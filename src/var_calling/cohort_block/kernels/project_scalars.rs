//! Phase A.1 layer 2 — per-(sample, allele) `AlleleSupportStats`
//! projection.
//!
//! Native column-native port of the row-shape `project_scalars` in
//! [`per_group_merger`](crate::var_calling::per_group_merger). Reads
//! the source pointers built by layer 1's
//! [`unify_alleles_columnar`](super::unify_alleles::unify_alleles_columnar)
//! and the raw `AlleleSupportStats` from the chunk's per-sample
//! columns, then folds them into a flat sample-major scalars matrix
//! `scalars[sample_idx * n_kept + allele_idx]` plus a per-sample
//! `other_scalars` pool of cap-dropped contributions.
//!
//! Four passes (same shape as the row-shape kernel):
//! 1. Sum per-position scalars for non-compound alleles.
//! 2. Pool the OTHER pool's sources into per-sample `other_scalars`.
//! 3. Project compound scalars (chain-anchored compounds get
//!    `count = |chain intersect|`, q-sum from the min-mean-q
//!    constituent, and bias counts scaled by `count / num_obs`).
//! 4. Subtract each compound's claim from its constituents to avoid
//!    double-counting at the group level.

use ahash::AHashMap;
use thiserror::Error;

use crate::pileup_record::AlleleSupportStats;
use crate::var_calling::cohort_block::columns::{MaterialisedChunk, SampleColumns};
use crate::var_calling::cohort_block::kernels::unify_alleles::UnifiedAllelesColumns;
use crate::var_calling::cohort_block::partition::WindowPartition;

/// Sample-major flat scalars matrix produced by
/// [`project_scalars_columnar`]. Same logical shape as
/// `MergedRecord.scalars` and `.other_scalars` so the downstream
/// likelihood kernel can consume either side without a branch.
///
/// **Layout.**
/// - [`Self::scalars`] — row-major sample-major. Cell for
///   `(sample_idx, allele_idx)` is
///   `scalars[sample_idx * n_kept + allele_idx]`. Length is
///   `n_samples * n_kept`.
/// - [`Self::other_scalars`] — one entry per sample; the pooled
///   contributions from alleles dropped by the max-alleles cap.
///   Length is `n_samples`.
// Mi1: `#[non_exhaustive]` — kernel-output columns.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectedScalarsColumns {
    pub n_samples: usize,
    pub n_kept: usize,
    pub scalars: Vec<AlleleSupportStats>,
    pub other_scalars: Vec<AlleleSupportStats>,
}

impl ProjectedScalarsColumns {
    pub fn empty() -> Self {
        Self {
            n_samples: 0,
            n_kept: 0,
            scalars: Vec::new(),
            other_scalars: Vec::new(),
        }
    }

    /// Reset every column to empty while preserving allocated
    /// capacity.
    pub fn clear(&mut self) {
        self.n_samples = 0;
        self.n_kept = 0;
        self.scalars.clear();
        self.other_scalars.clear();
    }

    /// Per-sample slice into [`Self::scalars`].
    pub fn scalars_row(&self, sample_idx: usize) -> &[AlleleSupportStats] {
        let start = sample_idx * self.n_kept;
        &self.scalars[start..start + self.n_kept]
    }
}

impl Default for ProjectedScalarsColumns {
    fn default() -> Self {
        Self::empty()
    }
}

/// Reusable scratch for [`project_scalars_columnar`]. Holds the
/// constituent-lookup map; cleared per call so the capacity carries
/// across groups.
#[derive(Debug, Default)]
pub struct ProjectScalarsScratch {
    /// `(sample_idx, record_idx_in_group, local_allele_idx) →
    /// allele_idx` for every non-compound entry's per-sample source.
    /// Built by [`build_source_index_columnar`] for the compound
    /// subtraction pass.
    pub(crate) source_index: AHashMap<(usize, usize, usize), usize>,
}

impl ProjectScalarsScratch {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn clear(&mut self) {
        self.source_index.clear();
    }
}

/// Errors surfaced by [`project_scalars_columnar`].
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum ProjectScalarsError {
    /// A chain-anchored compound's constituent record had
    /// `num_obs == 0`, which the compound's quality / subtraction
    /// math cannot handle. Signals an upstream invariant break — the
    /// walker should never emit a chain-id-touching allele with
    /// zero observations.
    #[error(
        "compound constituent with zero observations during {phase} at \
         group {chrom_id}:{group_start}..{group_end} \
         sample {sample_idx} allele {allele_idx} record {record_idx} \
         local_allele {local_allele_idx}"
    )]
    ZeroObservationConstituent {
        chrom_id: u32,
        group_start: u32,
        group_end: u32,
        phase: CompoundPhase,
        sample_idx: usize,
        allele_idx: usize,
        record_idx: usize,
        local_allele_idx: usize,
    },

    /// A compound allele had a non-zero `chain_anchor_count` but no
    /// constituent contributed a usable mean-q. Should be unreachable
    /// if `ZeroObservationConstituent` is not raised first.
    #[error(
        "chain-anchored compound at group {chrom_id}:{group_start}..{group_end} \
         sample {sample_idx} allele {allele_idx} has chain_count={inter} but \
         no constituent supplied a quality value"
    )]
    NoQualityForChainAnchoredCompound {
        chrom_id: u32,
        group_start: u32,
        group_end: u32,
        sample_idx: usize,
        allele_idx: usize,
        inter: u32,
    },
}

/// Phase tag for [`ProjectScalarsError::ZeroObservationConstituent`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompoundPhase {
    QualityGather,
    ConstituentSubtraction,
}

impl std::fmt::Display for CompoundPhase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CompoundPhase::QualityGather => f.write_str("quality-gather"),
            CompoundPhase::ConstituentSubtraction => f.write_str("constituent-subtraction"),
        }
    }
}

/// Phase A.1 layer 2 — full per-(sample, allele) scalar projection.
/// Composes the four sub-passes (per-position sum + OTHER pool +
/// compound projection + compound-constituent subtraction).
#[allow(clippy::too_many_arguments)]
pub fn project_scalars_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    unified: &UnifiedAllelesColumns,
    scratch: &mut ProjectScalarsScratch,
    out: &mut ProjectedScalarsColumns,
) -> Result<(), ProjectScalarsError> {
    let n_samples = unified.n_samples;
    let n_kept = unified.n_alleles();
    out.clear();
    out.n_samples = n_samples;
    out.n_kept = n_kept;
    out.scalars
        .resize(n_samples * n_kept, AlleleSupportStats::default());
    out.other_scalars
        .resize(n_samples, AlleleSupportStats::default());
    scratch.clear();

    sum_per_position_scalars_columnar(chunk, partition, group_idx, unified, &mut out.scalars);
    pool_dropped_other_scalars_columnar(
        chunk,
        partition,
        group_idx,
        unified,
        &mut out.other_scalars,
    );
    project_compound_scalars_columnar(chunk, partition, group_idx, unified, &mut out.scalars)?;
    build_source_index_columnar(unified, &mut scratch.source_index);
    subtract_compound_from_constituents_columnar(
        chunk,
        partition,
        group_idx,
        unified,
        &scratch.source_index,
        &mut out.scalars,
    )?;
    Ok(())
}

/// Pass 1: sum per-position-projection sources for every non-compound
/// allele into `scalars[sample × allele]`.
fn sum_per_position_scalars_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    unified: &UnifiedAllelesColumns,
    scalars: &mut [AlleleSupportStats],
) {
    let n_kept = unified.n_alleles();
    let n_samples = unified.n_samples;
    let group_position_lo = partition.position_range_for_group(group_idx).start;

    for allele_idx in 0..n_kept {
        if unified.is_compound[allele_idx] {
            continue;
        }
        for sample_idx in 0..n_samples {
            let (lo, hi) = unified.source_range(allele_idx, sample_idx);
            if lo == hi {
                continue;
            }
            let mut bucket = AlleleSupportStats::default();
            for k in lo..hi {
                let record_idx_in_group = unified.source_record_idx[k] as usize;
                let local_allele_idx = unified.source_local_allele_idx[k] as usize;
                if let Some(stats) = read_support_stats(
                    chunk,
                    partition,
                    group_position_lo,
                    record_idx_in_group,
                    sample_idx,
                    local_allele_idx,
                ) {
                    add_support(&mut bucket, &stats);
                }
            }
            scalars[sample_idx * n_kept + allele_idx] = bucket;
        }
    }
}

/// Pass 2: pool every OTHER-pool source's stats into per-sample
/// `other_scalars`. Cap-inactive groups land here as a fast no-op
/// because `other_offsets` is all-zero.
fn pool_dropped_other_scalars_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    unified: &UnifiedAllelesColumns,
    other_scalars: &mut [AlleleSupportStats],
) {
    let n_samples = unified.n_samples;
    let group_position_lo = partition.position_range_for_group(group_idx).start;
    #[allow(clippy::needless_range_loop)]
    for sample_idx in 0..n_samples {
        let (lo, hi) = unified.other_range(sample_idx);
        if lo == hi {
            continue;
        }
        for k in lo..hi {
            let record_idx_in_group = unified.other_record_idx[k] as usize;
            let local_allele_idx = unified.other_local_allele_idx[k] as usize;
            if let Some(stats) = read_support_stats(
                chunk,
                partition,
                group_position_lo,
                record_idx_in_group,
                sample_idx,
                local_allele_idx,
            ) {
                add_support(&mut other_scalars[sample_idx], &stats);
            }
        }
    }
}

/// Pass 3: for each chain-anchored compound × sample, fill the
/// `scalars[sample × allele]` cell from the constituent records'
/// stats. `num_obs` is the chain-id intersection count; `q_sum` uses
/// the homogeneous-quality approximation (min mean-q across
/// constituents × intersection). Bias counts (fwd, placed_left,
/// placed_start, mapq_sum, mapq_sum_sq) are scaled from the
/// pooled-constituent values by `count / bias_basis`.
fn project_compound_scalars_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    unified: &UnifiedAllelesColumns,
    scalars: &mut [AlleleSupportStats],
) -> Result<(), ProjectScalarsError> {
    let n_kept = unified.n_alleles();
    let n_samples = unified.n_samples;
    let group_position_lo = partition.position_range_for_group(group_idx).start;
    let group_start = partition.group_starts[group_idx];
    let group_end = partition.group_ends[group_idx];
    let chrom_id = chunk.chrom_id;

    for allele_idx in 0..n_kept {
        if !unified.is_compound[allele_idx] {
            continue;
        }
        for sample_idx in 0..n_samples {
            let inter = unified.chain_anchor_count(allele_idx, sample_idx);
            if inter == 0 {
                continue;
            }
            let (lo, hi) = unified.source_range(allele_idx, sample_idx);
            let mut min_mean_q: Option<f64> = None;
            let mut bias_fwd: u64 = 0;
            let mut bias_left: u64 = 0;
            let mut bias_start: u64 = 0;
            let mut bias_basis: u64 = 0;
            let mut bias_mapq_sum: u64 = 0;
            let mut bias_mapq_sum_sq: u128 = 0;
            for k in lo..hi {
                let record_idx_in_group = unified.source_record_idx[k] as usize;
                let local_allele_idx = unified.source_local_allele_idx[k] as usize;
                let Some(stats) = read_support_stats(
                    chunk,
                    partition,
                    group_position_lo,
                    record_idx_in_group,
                    sample_idx,
                    local_allele_idx,
                ) else {
                    continue;
                };
                if stats.num_obs == 0 {
                    return Err(ProjectScalarsError::ZeroObservationConstituent {
                        chrom_id,
                        group_start,
                        group_end,
                        phase: CompoundPhase::QualityGather,
                        sample_idx,
                        allele_idx,
                        record_idx: record_idx_in_group,
                        local_allele_idx,
                    });
                }
                let mean_q = stats.q_sum / f64::from(stats.num_obs);
                min_mean_q = Some(match min_mean_q {
                    Some(curr) => curr.min(mean_q),
                    None => mean_q,
                });
                bias_fwd += u64::from(stats.fwd);
                bias_left += u64::from(stats.placed_left);
                bias_start += u64::from(stats.placed_start);
                bias_basis += u64::from(stats.num_obs);
                bias_mapq_sum += u64::from(stats.mapq_sum);
                bias_mapq_sum_sq += stats.mapq_sum_sq as u128;
            }
            let count = inter;
            let mean_q =
                min_mean_q.ok_or(ProjectScalarsError::NoQualityForChainAnchoredCompound {
                    chrom_id,
                    group_start,
                    group_end,
                    sample_idx,
                    allele_idx,
                    inter,
                })?;
            let q_sum = mean_q * f64::from(count);
            let scale = if bias_basis == 0 {
                0.0
            } else {
                f64::from(count) / (bias_basis as f64)
            };
            scalars[sample_idx * n_kept + allele_idx] = AlleleSupportStats {
                num_obs: count,
                q_sum,
                fwd: ((bias_fwd as f64) * scale).round() as u32,
                placed_left: ((bias_left as f64) * scale).round() as u32,
                placed_start: ((bias_start as f64) * scale).round() as u32,
                mapq_sum: ((bias_mapq_sum as f64) * scale).round() as u32,
                mapq_sum_sq: ((bias_mapq_sum_sq as f64) * scale).round() as u64,
            };
        }
    }
    Ok(())
}

/// Build the `(sample, record_idx_in_group, local_allele_idx) →
/// allele_idx` map for every non-compound entry of `unified`. Used
/// by the compound-subtraction pass to find each constituent's kept
/// per-position allele in O(1).
fn build_source_index_columnar(
    unified: &UnifiedAllelesColumns,
    source_index: &mut AHashMap<(usize, usize, usize), usize>,
) {
    let n_kept = unified.n_alleles();
    let n_samples = unified.n_samples;
    source_index.clear();
    for allele_idx in 0..n_kept {
        if unified.is_compound[allele_idx] {
            continue;
        }
        for sample_idx in 0..n_samples {
            let (lo, hi) = unified.source_range(allele_idx, sample_idx);
            for k in lo..hi {
                let record_idx = unified.source_record_idx[k] as usize;
                let local_allele_idx = unified.source_local_allele_idx[k] as usize;
                source_index.insert((sample_idx, record_idx, local_allele_idx), allele_idx);
            }
        }
    }
}

/// Pass 4: in every chain-anchored sample, subtract the compound's
/// claim from each constituent's per-position scalars. Mirrors the
/// row-shape `subtract_compound_from_constituents` line by line.
fn subtract_compound_from_constituents_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    unified: &UnifiedAllelesColumns,
    source_index: &AHashMap<(usize, usize, usize), usize>,
    scalars: &mut [AlleleSupportStats],
) -> Result<(), ProjectScalarsError> {
    let n_kept = unified.n_alleles();
    let n_samples = unified.n_samples;
    let group_position_lo = partition.position_range_for_group(group_idx).start;
    let group_start = partition.group_starts[group_idx];
    let group_end = partition.group_ends[group_idx];
    let chrom_id = chunk.chrom_id;

    for allele_idx in 0..n_kept {
        if !unified.is_compound[allele_idx] {
            continue;
        }
        for sample_idx in 0..n_samples {
            let inter = unified.chain_anchor_count(allele_idx, sample_idx);
            if inter == 0 {
                continue;
            }
            let (lo, hi) = unified.source_range(allele_idx, sample_idx);
            for k in lo..hi {
                let record_idx_in_group = unified.source_record_idx[k] as usize;
                let local_allele_idx = unified.source_local_allele_idx[k] as usize;
                let Some(&constituent_idx) =
                    source_index.get(&(sample_idx, record_idx_in_group, local_allele_idx))
                else {
                    continue;
                };
                let Some(support) = read_support_stats(
                    chunk,
                    partition,
                    group_position_lo,
                    record_idx_in_group,
                    sample_idx,
                    local_allele_idx,
                ) else {
                    continue;
                };
                if support.num_obs == 0 {
                    return Err(ProjectScalarsError::ZeroObservationConstituent {
                        chrom_id,
                        group_start,
                        group_end,
                        phase: CompoundPhase::ConstituentSubtraction,
                        sample_idx,
                        allele_idx,
                        record_idx: record_idx_in_group,
                        local_allele_idx,
                    });
                }
                let dest_idx = sample_idx * n_kept + constituent_idx;
                let scale = f64::from(inter) / f64::from(support.num_obs);
                let clamped = scale.min(1.0);
                let mut to_subtract = AlleleSupportStats {
                    num_obs: inter.min(support.num_obs),
                    q_sum: support.q_sum * clamped,
                    fwd: ((f64::from(support.fwd)) * clamped).round() as u32,
                    placed_left: ((f64::from(support.placed_left)) * clamped).round() as u32,
                    placed_start: ((f64::from(support.placed_start)) * clamped).round() as u32,
                    mapq_sum: ((f64::from(support.mapq_sum)) * clamped).round() as u32,
                    mapq_sum_sq: ((support.mapq_sum_sq as f64) * clamped).round() as u64,
                };
                let dest = scalars[dest_idx];
                to_subtract.fwd = to_subtract.fwd.min(dest.fwd);
                to_subtract.placed_left = to_subtract.placed_left.min(dest.placed_left);
                to_subtract.placed_start = to_subtract.placed_start.min(dest.placed_start);
                to_subtract.mapq_sum = to_subtract.mapq_sum.min(dest.mapq_sum);
                to_subtract.mapq_sum_sq = to_subtract.mapq_sum_sq.min(dest.mapq_sum_sq);
                subtract_support(&mut scalars[dest_idx], &to_subtract);
            }
        }
    }
    Ok(())
}

/// Look up `(sample_idx, record_idx_in_group, local_allele_idx)`'s
/// `AlleleSupportStats` in the chunk's per-sample columns. Returns
/// `None` if the named sample has no record at the named position
/// in the group — the row-shape kernel skips silently in the same
/// case (the source pointer would be missing in
/// `per_sample_sources`).
pub(crate) fn read_support_stats(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_position_lo: usize,
    record_idx_in_group: usize,
    sample_idx: usize,
    local_allele_idx: usize,
) -> Option<AlleleSupportStats> {
    let row_idx = locate_sample_row_idx(
        partition,
        group_position_lo,
        record_idx_in_group,
        sample_idx,
    )?;
    let sample = &chunk.per_sample[sample_idx];
    Some(read_support_stats_from_columns(
        sample,
        row_idx,
        local_allele_idx,
    ))
}

/// Find `sample_idx`'s row index in `chunk.per_sample[sample_idx]` for
/// the position located at `group_position_lo + record_idx_in_group`.
/// Returns `None` if the sample has no record at that position.
pub(crate) fn locate_sample_row_idx(
    partition: &WindowPartition,
    group_position_lo: usize,
    record_idx_in_group: usize,
    sample_idx: usize,
) -> Option<usize> {
    let p = group_position_lo + record_idx_in_group;
    let sample_range = partition.sample_range_for_position(p);
    // M32: pin the `partition_window` invariant that
    // `samples_at_pos[sample_range]` is strictly sorted ascending by
    // sample_idx (each cohort sample contributes at most one record
    // per cohort position). The linear scan below relies on early-
    // return; without uniqueness a regression in `partition_window`'s
    // dedup could silently return the first row from a duplicated
    // sample without ever surfacing the error.
    debug_assert!(
        sample_range
            .clone()
            .map(|k| partition.samples_at_pos[k])
            .collect::<Vec<_>>()
            .windows(2)
            .all(|w| w[0] < w[1]),
        "locate_sample_row_idx: samples_at_pos[sample_range] must be strictly ascending",
    );
    // Linear scan over `samples_at_pos` at this position to find
    // `sample_idx` → row_idx. The slice is small (one entry per
    // sample-with-a-record-at-this-position) so the scan is fast.
    for k in sample_range {
        if partition.samples_at_pos[k] as usize == sample_idx {
            return Some(partition.rows_at_pos[k] as usize);
        }
    }
    None
}

/// Number of alleles emitted at `sample`'s `row_idx` row, derived
/// from the CSR `allele_offsets`. Used by the chain-broken-compound
/// path of layer 3 to iterate the row's per-allele stats.
pub(crate) fn n_local_alleles_at_row(sample: &SampleColumns, row_idx: usize) -> usize {
    let lo = sample.allele_offsets[row_idx] as usize;
    let hi = sample.allele_offsets[row_idx + 1] as usize;
    hi - lo
}

/// Extract the 7-component `AlleleSupportStats` for the allele at
/// `(row_idx, local_allele_idx)` from `sample`'s parallel column
/// vectors.
pub(crate) fn read_support_stats_from_columns(
    sample: &SampleColumns,
    row_idx: usize,
    local_allele_idx: usize,
) -> AlleleSupportStats {
    let allele_lo = sample.allele_offsets[row_idx] as usize;
    let k = allele_lo + local_allele_idx;
    AlleleSupportStats::new(
        sample.allele_num_obs[k],
        sample.allele_q_sum[k],
        sample.allele_fwd[k],
        sample.allele_placed_left[k],
        sample.allele_placed_start[k],
        sample.allele_mapq_sum[k],
        sample.allele_mapq_sum_sq[k],
    )
}

fn add_support(into: &mut AlleleSupportStats, src: &AlleleSupportStats) {
    let AlleleSupportStats {
        num_obs,
        q_sum,
        fwd,
        placed_left,
        placed_start,
        mapq_sum,
        mapq_sum_sq,
    } = *src;
    into.num_obs = into.num_obs.saturating_add(num_obs);
    into.q_sum += q_sum;
    into.fwd = into.fwd.saturating_add(fwd);
    into.placed_left = into.placed_left.saturating_add(placed_left);
    into.placed_start = into.placed_start.saturating_add(placed_start);
    into.mapq_sum = into.mapq_sum.saturating_add(mapq_sum);
    into.mapq_sum_sq = into.mapq_sum_sq.saturating_add(mapq_sum_sq);
}

fn subtract_support(into: &mut AlleleSupportStats, src: &AlleleSupportStats) {
    let AlleleSupportStats {
        num_obs,
        q_sum,
        fwd,
        placed_left,
        placed_start,
        mapq_sum,
        mapq_sum_sq,
    } = *src;
    into.num_obs = into.num_obs.saturating_sub(num_obs);
    // q_sum stays non-positive — see the row-shape kernel's note.
    into.q_sum = (into.q_sum - q_sum).min(0.0);
    into.fwd = into.fwd.saturating_sub(fwd);
    into.placed_left = into.placed_left.saturating_sub(placed_left);
    into.placed_start = into.placed_start.saturating_sub(placed_start);
    into.mapq_sum = into.mapq_sum.saturating_sub(mapq_sum);
    into.mapq_sum_sq = into.mapq_sum_sq.saturating_sub(mapq_sum_sq);
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::fasta::ChromRefFetcher;
    use crate::fasta::fetcher::ChromRefFetchError;
    use crate::pileup_record::{AlleleObservation, PileupRecord};
    use crate::var_calling::cohort_block::kernels::unify_alleles::{
        UnifiedAllelesColumns, UnifyAllelesScratch, unify_alleles_columnar,
    };
    use crate::var_calling::cohort_block::partition::{
        PartitionScratch, WindowPartition, partition_window,
    };
    use crate::var_calling::cohort_block::test_helpers::{loaded_chunk, record, ref_plus_alt};
    use crate::var_calling::per_group_merger::{
        PerGroupMerger, PerGroupMergerConfig, SharedRefFetcher,
    };
    use crate::var_calling::variant_grouping::{GrouperError, OverlappingVariantGroup};

    #[derive(Clone)]
    struct MockRef {
        seq: Vec<u8>,
        base_offset: u32,
    }

    impl crate::fasta::fetcher::sealed::Sealed for MockRef {}
    impl ChromRefFetcher for MockRef {
        fn length(&self) -> u32 {
            self.base_offset.saturating_sub(1) + self.seq.len() as u32
        }
        fn fetch(&self, start_1based: u32, length: u32) -> Result<Vec<u8>, ChromRefFetchError> {
            let start_idx = (start_1based - self.base_offset) as usize;
            Ok(self.seq[start_idx..start_idx + length as usize].to_vec())
        }
        fn iter_bases<'a>(
            &'a self,
        ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>
        {
            Ok(Box::new(self.seq.iter().copied().map(Ok)))
        }
    }

    fn shared_mock(seq: &[u8], base_offset: u32) -> SharedRefFetcher {
        Arc::new(MockRef {
            seq: seq.to_vec(),
            base_offset,
        })
    }

    fn group_fixture(
        per_sample_records: Vec<Vec<PileupRecord>>,
        ref_seq_full: &[u8],
        window: std::ops::Range<u32>,
        max_group_span: u32,
    ) -> (
        crate::var_calling::cohort_block::columns::MaterialisedChunk,
        WindowPartition,
        Vec<u8>,
    ) {
        let n = per_sample_records.len();
        let (chunk, _) = loaded_chunk(0, 1..200, per_sample_records);
        let mut scratch = PartitionScratch::with_n_samples(n);
        let mut partition = WindowPartition::empty();
        partition_window(
            &chunk,
            &window,
            &[],
            max_group_span,
            &mut scratch,
            &mut partition,
        )
        .unwrap();
        let group_start = partition.group_starts[0];
        let group_end = partition.group_ends[0];
        let ref_seq = ref_seq_full[(group_start - 1) as usize..group_end as usize].to_vec();
        (chunk, partition, ref_seq)
    }

    /// Pull `(MergedRecord.scalars, .other_scalars)` from the
    /// row-shape kernel — the byte-identity oracle for layer 2.
    #[allow(clippy::type_complexity)]
    fn scalars_via_existing_kernel(
        chunk: &MaterialisedChunk,
        partition: &WindowPartition,
        group_idx: usize,
        ref_fetcher: SharedRefFetcher,
        max_alleles: usize,
    ) -> (Vec<AlleleSupportStats>, Vec<AlleleSupportStats>) {
        let group = crate::var_calling::cohort_block::test_helpers::build_overlapping_variant_group(
            chunk,
            partition,
            group_idx,
            chunk.n_samples(),
            chunk.chrom_id,
        );
        let config = PerGroupMergerConfig::new(2, max_alleles, 64, 32).expect("merger config");
        let iter: Vec<Result<OverlappingVariantGroup, GrouperError>> = vec![Ok(group)];
        let mut merger = PerGroupMerger::with_config(iter.into_iter(), ref_fetcher, config);
        let record = merger
            .next()
            .expect("merger yields at least one item")
            .expect("merger succeeded");
        (record.scalars, record.other_scalars)
    }

    fn project_columnar(
        chunk: &MaterialisedChunk,
        partition: &WindowPartition,
        group_idx: usize,
        ref_seq: &[u8],
        max_alleles: usize,
    ) -> ProjectedScalarsColumns {
        let mut unify_scratch = UnifyAllelesScratch::new();
        let mut unified = UnifiedAllelesColumns::empty();
        unify_alleles_columnar(
            chunk,
            partition,
            group_idx,
            ref_seq,
            chunk.n_samples(),
            max_alleles,
            &mut unify_scratch,
            &mut unified,
        )
        .expect("unify succeeded");

        let mut scratch = ProjectScalarsScratch::new();
        let mut out = ProjectedScalarsColumns::empty();
        project_scalars_columnar(
            chunk,
            partition,
            group_idx,
            &unified,
            &mut scratch,
            &mut out,
        )
        .expect("project scalars succeeded");
        out
    }

    #[test]
    fn project_scalars_no_compounds_matches_row_shape() {
        // Two samples, SNP A→T at pos 50.
        let s0 = vec![record(50, ref_plus_alt(3, 4))];
        let s1 = vec![record(50, ref_plus_alt(5, 6))];
        let ref_full = vec![b'A'; 200];
        let (chunk, partition, ref_seq) = group_fixture(vec![s0, s1], &ref_full, 1..200, 50);

        let out = project_columnar(&chunk, &partition, 0, &ref_seq, 16);
        let (oracle_scalars, oracle_other) =
            scalars_via_existing_kernel(&chunk, &partition, 0, shared_mock(&ref_full, 1), 16);

        assert_eq!(out.scalars, oracle_scalars);
        assert_eq!(out.other_scalars, oracle_other);
    }

    #[test]
    fn project_scalars_multi_position_no_compound_matches_row_shape() {
        // Same MNP+SNP fixture used in layer 1's byte-identity test.
        let mnp_record = PileupRecord::new(
            0,
            50,
            vec![
                AlleleObservation::new(
                    b"ACA".to_vec(),
                    AlleleSupportStats::new(4, -1.0, 4, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    b"ACG".to_vec(),
                    AlleleSupportStats::new(5, -1.0, 5, 0, 0, 0, 0),
                    Vec::new(),
                ),
            ],
        );
        let snp_record = record(52, ref_plus_alt(2, 3));
        let s0 = vec![mnp_record, snp_record];
        let s1 = vec![record(50, ref_plus_alt(6, 0))];
        let mut ref_full = vec![b'A'; 200];
        ref_full[49] = b'A';
        ref_full[50] = b'C';
        ref_full[51] = b'A';
        let (chunk, partition, ref_seq) = group_fixture(vec![s0, s1], &ref_full, 1..200, 50);
        assert_eq!(ref_seq.as_slice(), b"ACA");

        let out = project_columnar(&chunk, &partition, 0, &ref_seq, 16);
        let (oracle_scalars, oracle_other) =
            scalars_via_existing_kernel(&chunk, &partition, 0, shared_mock(&ref_full, 1), 16);

        assert_eq!(out.scalars, oracle_scalars);
        assert_eq!(out.other_scalars, oracle_other);
    }

    #[test]
    fn project_scalars_with_compound_matches_row_shape() {
        // Same compound fixture as layer 1.2.
        let mnp_rec = PileupRecord::new(
            0,
            100,
            vec![
                AlleleObservation::new(
                    b"ACA".to_vec(),
                    AlleleSupportStats::new(4, -1.0, 4, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    b"ACC".to_vec(),
                    AlleleSupportStats::new(3, -1.0, 3, 0, 0, 0, 0),
                    vec![100],
                ),
            ],
        );
        let snp_100 = PileupRecord::new(
            0,
            100,
            vec![
                AlleleObservation::new(
                    b"A".to_vec(),
                    AlleleSupportStats::new(5, -1.0, 5, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    b"T".to_vec(),
                    AlleleSupportStats::new(4, -1.0, 4, 0, 0, 0, 0),
                    vec![1],
                ),
            ],
        );
        let snp_102 = PileupRecord::new(
            0,
            102,
            vec![
                AlleleObservation::new(
                    b"A".to_vec(),
                    AlleleSupportStats::new(5, -1.0, 5, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    b"G".to_vec(),
                    AlleleSupportStats::new(4, -1.0, 4, 0, 0, 0, 0),
                    vec![1],
                ),
            ],
        );
        let s0 = vec![snp_100, snp_102];
        let s1 = vec![mnp_rec];
        let mut ref_full = vec![b'A'; 200];
        ref_full[99] = b'A';
        ref_full[100] = b'C';
        ref_full[101] = b'A';
        let (chunk, partition, ref_seq) = group_fixture(vec![s0, s1], &ref_full, 1..200, 50);

        let out = project_columnar(&chunk, &partition, 0, &ref_seq, 16);
        let (oracle_scalars, oracle_other) =
            scalars_via_existing_kernel(&chunk, &partition, 0, shared_mock(&ref_full, 1), 16);

        assert_eq!(out.scalars, oracle_scalars);
        assert_eq!(out.other_scalars, oracle_other);
    }

    #[test]
    fn project_scalars_with_cap_drop_pools_other_matches_row_shape() {
        // Five-ALT fixture with max_alleles=2 → 3 ALTs go to OTHER.
        let alts = vec![
            AlleleObservation::new(
                b"A".to_vec(),
                AlleleSupportStats::new(10, -1.0, 10, 0, 0, 0, 0),
                Vec::new(),
            ),
            AlleleObservation::new(
                b"T".to_vec(),
                AlleleSupportStats::new(8, -1.0, 8, 0, 0, 0, 0),
                Vec::new(),
            ),
            AlleleObservation::new(
                b"C".to_vec(),
                AlleleSupportStats::new(3, -1.0, 3, 0, 0, 0, 0),
                Vec::new(),
            ),
            AlleleObservation::new(
                b"G".to_vec(),
                AlleleSupportStats::new(2, -1.0, 2, 0, 0, 0, 0),
                Vec::new(),
            ),
            AlleleObservation::new(
                b"N".to_vec(),
                AlleleSupportStats::new(1, -1.0, 1, 0, 0, 0, 0),
                Vec::new(),
            ),
        ];
        let s0 = vec![PileupRecord::new(0, 50, alts)];
        let ref_full = vec![b'A'; 200];
        let (chunk, partition, ref_seq) = group_fixture(vec![s0], &ref_full, 1..200, 50);

        let out = project_columnar(&chunk, &partition, 0, &ref_seq, 2);
        let (oracle_scalars, oracle_other) =
            scalars_via_existing_kernel(&chunk, &partition, 0, shared_mock(&ref_full, 1), 2);

        assert_eq!(out.scalars, oracle_scalars);
        assert_eq!(out.other_scalars, oracle_other);
    }
}
