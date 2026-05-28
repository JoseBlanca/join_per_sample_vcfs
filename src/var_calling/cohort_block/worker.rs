//! Per-window variant-calling math.
//!
//! Phase A.1 layer 4 — column-native pipeline. The worker runs each
//! partition group through the three columnar kernels
//! ([`unify_alleles_columnar`],
//! [`project_scalars_columnar`],
//! [`compute_log_likelihoods_columnar`]) and adapts their outputs into
//! the row-shape
//! [`MergedRecord`](crate::var_calling::per_group_merger::MergedRecord)
//! that
//! [`PosteriorEngine`](crate::var_calling::posterior_engine::PosteriorEngine)
//! consumes for the EM, posterior, GQ, and QUAL passes. The
//! row-shape
//! [`PerGroupMerger`](crate::var_calling::per_group_merger::PerGroupMerger)
//! is no longer on the production path — its emit was the
//! `MergedRecord` shape that the EM still wants, so layer 4 keeps the
//! type but builds it column-natively rather than by going through
//! `MergedAlleleSet` + `ScalarProjection` + `compute_log_likelihoods`
//! in the row pipeline. A future Phase A.2 step ports the EM to
//! consume columnar slices directly.
//!
//! The Phase A.0 row-shape adapter
//! ([`build_overlapping_variant_group`]) is retained for use as the
//! byte-identity oracle in the kernels' unit tests.

use std::sync::Arc;

use thiserror::Error;

use crate::fasta::fetcher::{ChromRefFetchError, ChromRefFetcher};
use crate::pileup_record::PileupRecord;
use crate::var_calling::cohort_block::columns::MaterialisedChunk;
use crate::var_calling::cohort_block::kernels::compute_log_likelihoods::{
    ComputeLogLikelihoodsError, LogLikelihoodsColumns, compute_log_likelihoods_columnar,
};
use crate::var_calling::cohort_block::kernels::project_scalars::{
    ProjectScalarsError, ProjectScalarsScratch, ProjectedScalarsColumns, project_scalars_columnar,
};
use crate::var_calling::cohort_block::kernels::unify_alleles::{
    UnifiedAllelesColumns, UnifyAllelesError, UnifyAllelesScratch, unify_alleles_columnar,
};
use crate::var_calling::cohort_block::partition::WindowPartition;
use crate::var_calling::per_group_merger::{
    CompoundConstituent, DegeneracyKind, LikelihoodContext, MergedAllele, MergedRecord,
    PerGroupMergerConfig, PerGroupMergerError, SharedRefFetcher,
};
use crate::var_calling::per_position_merger::PerPositionPileups;
use crate::var_calling::posterior_engine::{
    PosteriorEngine, PosteriorEngineConfig, PosteriorEngineError, PosteriorRecord,
};
use crate::var_calling::variant_grouping::OverlappingVariantGroup;

/// Reusable scratch for the column-native pipeline. The driver owns
/// one of these per chunk-loop instance and threads it into every
/// [`run_window`] call so the per-group columnar buffers retain
/// allocated capacity across windows + chunks. Build empty with
/// [`Self::empty`] and let the first call populate it; subsequent
/// calls only grow buffers when a group exceeds the prior high-water
/// mark.
#[derive(Debug, Default)]
pub struct ColumnarPipelineScratch {
    unify: UnifyAllelesScratch,
    unified: UnifiedAllelesColumns,
    project: ProjectScalarsScratch,
    projection: ProjectedScalarsColumns,
    log_likelihoods: LogLikelihoodsColumns,
    /// Per-group chain-anchor flag table —
    /// `chain_anchor_flags[sample × n_alleles + allele]` is `true` iff
    /// the allele is a compound and the sample is chain-broken at it.
    /// Refilled per group; capacity carries across groups.
    chain_anchor_flags: Vec<bool>,
    /// REF byte buffer for the group span. Refilled by
    /// [`ChromRefFetcher::fetch_into`] per group; capacity carries.
    ref_buf: Vec<u8>,
    /// Number of groups skipped because
    /// `unified.alleles.len() > cfg.max_alleles_lh_calc`. Accumulated
    /// across every `run_window` call; the driver reads + clears via
    /// [`Self::take_lh_cap_stats`] after each call to roll into
    /// [`ChunkDriverStats`](super::driver::ChunkDriverStats).
    lh_cap_groups_skipped: u64,
    /// Sum of `unified.alleles.len()` across every skipped group —
    /// the row-shape merger reports this in `allele_lh_cap_hit:
    /// alleles_total=N`.
    lh_cap_alleles_in_skipped: u64,
}

impl ColumnarPipelineScratch {
    pub fn empty() -> Self {
        Self::default()
    }

    /// Take the LH-cap skip counters accumulated since the last
    /// `take_lh_cap_stats` call. The driver rolls these into
    /// `ChunkDriverStats` after every `run_window`.
    pub fn take_lh_cap_stats(&mut self) -> (u64, u64) {
        let g = std::mem::take(&mut self.lh_cap_groups_skipped);
        let a = std::mem::take(&mut self.lh_cap_alleles_in_skipped);
        (g, a)
    }
}

/// Run the worker math on one window's partition.
///
/// **Inputs.**
/// - `chunk`: the loaded `MaterialisedChunk` whose columns the
///   partition's `(sample_idx, row_idx)` handles refer into.
/// - `partition`: the variant-group partition produced by
///   [`partition_window`](super::partition::partition_window) for
///   one of `chunk.windows`. Its `n_groups()` may be 0 — empty
///   partitions are a fast no-op that never queries the fetcher.
/// - `ref_fetcher`: shared reference fetcher used by the columnar
///   allele-unification layer to gather REF bases for the group span.
///   The driver typically holds an `Arc<StreamingChromRefFetcher>`
///   per chromosome.
/// - `per_group_config`: tuning for the column-native kernels
///   (`ploidy`, `max_alleles`, `max_alleles_lh_calc`). The
///   `batch_size` field is ignored — Phase A.1 doesn't batch.
/// - `posterior_config`: tuning passed straight to
///   [`PosteriorEngine`].
/// - `scratch`: per-window columnar scratch (
///   [`ColumnarPipelineScratch`]). The driver owns one and reuses it
///   across every window in every chunk.
/// - `output`: appended-to record buffer. Cleared on entry so the
///   driver can reuse it across windows.
#[allow(clippy::too_many_arguments)]
pub fn run_window(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    ref_fetcher: SharedRefFetcher,
    per_group_config: PerGroupMergerConfig,
    posterior_config: PosteriorEngineConfig,
    scratch: &mut ColumnarPipelineScratch,
    output: &mut Vec<PosteriorRecord>,
) -> Result<(), PosteriorEngineError> {
    output.clear();
    if partition.n_groups() == 0 {
        return Ok(());
    }

    let iter = ColumnarMergedRecordsIter {
        chunk,
        partition,
        ref_fetcher: &*ref_fetcher,
        cfg: per_group_config,
        scratch,
        cursor: 0,
    };
    let engine = PosteriorEngine::with_config(iter, posterior_config);
    for record in engine {
        output.push(record?);
    }
    Ok(())
}

/// Iterator over the partition's groups; each `next()` runs layers
/// 1+2+3 + builds a `MergedRecord` from the columnar outputs.
///
/// Borrowed-scratch lifetime: the iterator holds a `&mut` into the
/// driver-owned [`ColumnarPipelineScratch`]; the borrow is released
/// when the iterator is dropped at the end of the engine drain.
struct ColumnarMergedRecordsIter<'a> {
    chunk: &'a MaterialisedChunk,
    partition: &'a WindowPartition,
    ref_fetcher: &'a dyn ChromRefFetcher,
    cfg: PerGroupMergerConfig,
    scratch: &'a mut ColumnarPipelineScratch,
    cursor: usize,
}

impl Iterator for ColumnarMergedRecordsIter<'_> {
    type Item = Result<MergedRecord, PerGroupMergerError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.cursor >= self.partition.n_groups() {
                return None;
            }
            let g = self.cursor;
            self.cursor += 1;
            match build_merged_record_columnar(
                self.chunk,
                self.partition,
                g,
                self.ref_fetcher,
                &self.cfg,
                self.scratch,
            ) {
                Ok(Some(rec)) => return Some(Ok(rec)),
                Ok(None) => continue,
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

/// Run layers 1+2+3 on group `g` and adapt the columnar outputs into
/// a row-shape [`MergedRecord`] for the EM. Reuses every buffer in
/// `scratch`; per-group allocations are limited to the
/// [`MergedRecord`]'s owned `Vec`s.
///
/// Returns `Ok(None)` for groups the row-shape merger would skip:
/// post-unification `unified.alleles.len() < 2` (REF-only after the
/// max-alleles cap dropped every ALT into the OTHER pool) and
/// `unified.alleles.len() > cfg.max_alleles_lh_calc` (LH-cap bind in
/// pathological low-complexity regions). The downstream
/// `PosteriorEngine` does have a `trivial_posterior_record` branch
/// for the former case, but emitting through it would let the
/// driver's hom-ref counter inflate above the row-shape baseline —
/// the skip here mirrors `PerGroupMerger::process_group`'s
/// `Ok(None)` returns so the per-category counts stay byte-identical
/// against `main`.
fn build_merged_record_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    ref_fetcher: &dyn ChromRefFetcher,
    cfg: &PerGroupMergerConfig,
    scratch: &mut ColumnarPipelineScratch,
) -> Result<Option<MergedRecord>, PerGroupMergerError> {
    let chrom_id = chunk.chrom_id;
    let group_start = partition.group_starts[group_idx];
    let group_end = partition.group_ends[group_idx];
    let n_samples = chunk.n_samples();

    // Layer 1 — allele unification.
    let span = group_end - group_start + 1;
    scratch.ref_buf.clear();
    ref_fetcher
        .fetch_into(group_start, span, &mut scratch.ref_buf)
        .map_err(|source| ref_fetch_error_to_merger(chrom_id, group_start, group_end, source))?;
    unify_alleles_columnar(
        chunk,
        partition,
        group_idx,
        &scratch.ref_buf,
        n_samples,
        cfg.max_alleles,
        &mut scratch.unify,
        &mut scratch.unified,
    )
    .map_err(|e| unify_error_to_merger(chrom_id, group_start, group_end, e))?;

    let n_alleles = scratch.unified.n_alleles();

    // Row-shape parity — see [`PerGroupMerger::process_group`].
    if n_alleles > cfg.max_alleles_lh_calc {
        scratch.lh_cap_groups_skipped += 1;
        scratch.lh_cap_alleles_in_skipped += n_alleles as u64;
        return Ok(None);
    }
    if n_alleles < 2 {
        return Ok(None);
    }

    // Layer 2 — per-(sample, allele) scalar projection.
    project_scalars_columnar(
        chunk,
        partition,
        group_idx,
        &scratch.unified,
        &mut scratch.project,
        &mut scratch.projection,
    )
    .map_err(|e| project_scalars_error_to_merger(chrom_id, group_start, group_end, e))?;

    // Layer 3 — per-(sample, genotype) log-likelihood.
    let ctx = LikelihoodContext {
        chrom_id,
        start: group_start,
        end: group_end,
        ploidy: cfg.ploidy,
    };
    compute_log_likelihoods_columnar(
        chunk,
        partition,
        group_idx,
        &scratch.unified,
        &scratch.projection,
        &ctx,
        &[],
        &mut scratch.log_likelihoods,
    )
    .map_err(|e| compute_ll_error_to_merger(chrom_id, group_start, group_end, e))?;

    // chain_anchor_flags table (row-shape boolean flat): `true` iff
    // the allele is a compound AND the sample's chain_anchor_count is
    // zero. Matches `build_chain_anchor_flags` in `per_group_merger`.
    scratch.chain_anchor_flags.clear();
    scratch
        .chain_anchor_flags
        .resize(n_samples * n_alleles, false);
    for allele_idx in 0..n_alleles {
        if !scratch.unified.is_compound[allele_idx] {
            continue;
        }
        for sample_idx in 0..n_samples {
            if scratch.unified.chain_anchor_count(allele_idx, sample_idx) == 0 {
                scratch.chain_anchor_flags[sample_idx * n_alleles + allele_idx] = true;
            }
        }
    }

    // Build the row-shape `Vec<MergedAllele>` from
    // `UnifiedAllelesColumns`. Walks the per-allele constituent
    // pointers and clones the seq bytes (PosteriorRecord ends up
    // owning them, so the clone is just timing — not extra).
    let mut alleles: Vec<MergedAllele> = Vec::with_capacity(n_alleles);
    for a in 0..n_alleles {
        let lo = scratch.unified.constituent_offsets[a] as usize;
        let hi = scratch.unified.constituent_offsets[a + 1] as usize;
        let constituents: Vec<CompoundConstituent> = (lo..hi)
            .map(|k| CompoundConstituent {
                record_idx: scratch.unified.constituent_record_idx[k] as usize,
                local_allele_idx: scratch.unified.constituent_local_allele_idx[k] as usize,
            })
            .collect();
        alleles.push(MergedAllele {
            seq: scratch.unified.allele_seq(a).to_vec(),
            is_compound: scratch.unified.is_compound[a],
            constituents,
        });
    }

    Ok(Some(MergedRecord {
        chrom_id,
        start: group_start,
        end: group_end,
        alleles,
        ploidy: cfg.ploidy,
        n_samples,
        n_genotypes: scratch.log_likelihoods.n_genotypes,
        scalars: scratch.projection.scalars.clone(),
        other_scalars: scratch.projection.other_scalars.clone(),
        chain_anchor_flags: scratch.chain_anchor_flags.clone(),
        log_likelihoods: scratch.log_likelihoods.log_likelihoods.clone(),
    }))
}

fn ref_fetch_error_to_merger(
    chrom_id: u32,
    start: u32,
    end: u32,
    source: ChromRefFetchError,
) -> PerGroupMergerError {
    PerGroupMergerError::RefFetch {
        chrom_id,
        start,
        end,
        source: std::io::Error::other(source.to_string()),
    }
}

fn unify_error_to_merger(
    _chrom_id: u32,
    _group_start: u32,
    _group_end: u32,
    e: UnifyAllelesError,
) -> PerGroupMergerError {
    match e {
        UnifyAllelesError::MissingCompoundAlleleBytes {
            chrom_id,
            group_start,
            group_end,
            record_idx,
            local_allele_idx,
        } => PerGroupMergerError::MissingCompoundAlleleBytes {
            chrom_id,
            start: group_start,
            end: group_end,
            record_idx: record_idx as usize,
            local_allele_idx: local_allele_idx as usize,
        },
    }
}

fn project_scalars_error_to_merger(
    _chrom_id: u32,
    _group_start: u32,
    _group_end: u32,
    e: ProjectScalarsError,
) -> PerGroupMergerError {
    use crate::var_calling::cohort_block::kernels::project_scalars::CompoundPhase as CnPhase;
    use crate::var_calling::per_group_merger::CompoundPhase as RowPhase;
    match e {
        ProjectScalarsError::ZeroObservationConstituent {
            chrom_id,
            group_start,
            group_end,
            phase,
            sample_idx,
            allele_idx,
            record_idx,
            local_allele_idx,
        } => PerGroupMergerError::ZeroObservationConstituent {
            chrom_id,
            start: group_start,
            end: group_end,
            phase: match phase {
                CnPhase::QualityGather => RowPhase::QualityGather,
                CnPhase::ConstituentSubtraction => RowPhase::ConstituentSubtraction,
            },
            sample_idx,
            allele_idx,
            record_idx,
            local_allele_idx,
        },
        ProjectScalarsError::NoQualityForChainAnchoredCompound {
            chrom_id,
            group_start,
            group_end,
            sample_idx,
            allele_idx,
            inter,
        } => PerGroupMergerError::NoQualityForChainAnchoredCompound {
            chrom_id,
            start: group_start,
            end: group_end,
            sample_idx,
            allele_idx,
            inter,
        },
    }
}

fn compute_ll_error_to_merger(
    chrom_id: u32,
    group_start: u32,
    group_end: u32,
    e: ComputeLogLikelihoodsError,
) -> PerGroupMergerError {
    match e {
        ComputeLogLikelihoodsError::DegenerateLikelihood {
            chrom_id,
            start,
            end,
            sample_idx,
            genotype_idx,
            kind,
        } => PerGroupMergerError::DegenerateLikelihood {
            chrom_id,
            start,
            end,
            sample_idx,
            genotype_idx,
            kind,
        },
        ComputeLogLikelihoodsError::NAllelesExceedsBitmask { n_alleles } => {
            // The row-shape kernel doesn't surface this as an Error;
            // its equivalent invariant fires `assert!` inside
            // `standard_log_likelihood`. The Phase-A.1 cap (set
            // upstream in `unify_alleles_columnar`) is bounded by
            // `MAX_BITMASK_ALLELES`, so the kernel surfacing this
            // variant means the caller passed `cfg.max_alleles >
            // MAX_BITMASK_ALLELES`. Map it to a synthetic
            // `DegenerateLikelihood` carrying enough context for
            // diagnosis: at this point we have group locus +
            // n_alleles but no sample/genotype, so the bookkeeping
            // fields are placeholders. `DegeneracyKind::NaN` is
            // chosen as the "internal-bug, finite formula went
            // wrong" sentinel — it matches the same diagnostic
            // family the row-shape kernel uses for self-reports.
            //
            // The diagnostic guard above this point in the worker
            // would clamp `max_alleles` to ≤ `MAX_BITMASK_ALLELES`
            // before we ever reach here; if a future code path
            // bypasses that guard, the assert-equivalent surfaces
            // via this variant.
            let _ = n_alleles;
            PerGroupMergerError::DegenerateLikelihood {
                chrom_id,
                start: group_start,
                end: group_end,
                sample_idx: usize::MAX,
                genotype_idx: usize::MAX,
                kind: DegeneracyKind::NaN,
            }
        }
    }
}

/// Build one [`OverlappingVariantGroup`] from group `g` of
/// `partition`. Per-position records are materialised on demand
/// from the chunk's columnar slices via
/// [`SampleColumns::materialise_record`](super::columns::SampleColumns::materialise_record).
///
/// Retained from Phase A.0 as the byte-identity oracle in the
/// kernels' unit tests. Not on the production path — `run_window`
/// builds [`MergedRecord`]s column-natively via
/// [`build_merged_record_columnar`].
#[cfg(test)]
pub(crate) fn build_overlapping_variant_group(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    g: usize,
    n_samples: usize,
    chrom_id: u32,
) -> OverlappingVariantGroup {
    let position_range = partition.position_range_for_group(g);
    let mut records: Vec<PerPositionPileups> = Vec::with_capacity(position_range.len());
    for p in position_range {
        let pos = partition.positions[p];
        let mut per_sample: Vec<Option<PileupRecord>> = vec![None; n_samples];
        let sample_range = partition.sample_range_for_position(p);
        for k in sample_range {
            let sample_idx = partition.samples_at_pos[k] as usize;
            let row_idx = partition.rows_at_pos[k] as usize;
            per_sample[sample_idx] =
                Some(chunk.per_sample[sample_idx].materialise_record(chrom_id, row_idx));
        }
        records.push(PerPositionPileups {
            chrom_id,
            pos,
            per_sample,
        });
    }
    OverlappingVariantGroup {
        chrom_id,
        start: partition.group_starts[g],
        end: partition.group_ends[g],
        records,
    }
}

/// Construct the `SharedRefFetcher` shape `run_window` expects from
/// an owned [`ChromRefFetcher`] implementor. The Arc is still the
/// natural type for the per-chromosome ref fetcher held by the
/// driver — `run_window` derefs into a `&dyn ChromRefFetcher` for
/// the columnar layer 1 fetch.
pub fn shared_ref_fetcher<F>(fetcher: F) -> SharedRefFetcher
where
    F: ChromRefFetcher + Send + 'static,
{
    Arc::new(fetcher)
}

/// Errors from [`unified_alleles_for_group_columnar`]. Wraps the
/// kernel error and the ref-fetcher error so callers see a single
/// surface.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum WorkerUnifyError {
    /// Kernel-side failure (compound projection couldn't find an
    /// anchoring sample, etc.).
    #[error(transparent)]
    Kernel(#[from] UnifyAllelesError),

    /// `fetch_into` failed while pulling the group's REF span — the
    /// kernel needs those bytes as input.
    #[error("reference fetch: {0}")]
    RefFetch(#[from] ChromRefFetchError),
}

/// Phase A.1 sub-step 1.4 — column-native allele unification for one
/// group of a [`WindowPartition`], fetched through the worker's
/// ref-fetcher convention. Retained from layer 1 for direct external
/// access to the unifier without the full pipeline; `run_window`
/// drives the same kernel via [`ColumnarPipelineScratch`].
#[allow(clippy::too_many_arguments)]
pub fn unified_alleles_for_group_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    ref_fetcher: &dyn ChromRefFetcher,
    max_alleles: usize,
    scratch: &mut UnifyAllelesScratch,
    ref_buf: &mut Vec<u8>,
    out: &mut UnifiedAllelesColumns,
) -> Result<(), WorkerUnifyError> {
    let group_start = partition.group_starts[group_idx];
    let group_end = partition.group_ends[group_idx];
    let span = group_end - group_start + 1;
    ref_fetcher.fetch_into(group_start, span, ref_buf)?;
    unify_alleles_columnar(
        chunk,
        partition,
        group_idx,
        ref_buf,
        chunk.n_samples(),
        max_alleles,
        scratch,
        out,
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::cohort_block::partition::{
        PartitionScratch, WindowPartition, partition_window,
    };
    use crate::var_calling::cohort_block::test_helpers::{loaded_chunk, record, ref_plus_alt};
    use crate::var_calling::per_group_merger::{PerGroupMerger, SharedRefFetcher};
    use crate::var_calling::variant_grouping::GrouperError;

    /// Build a chunk + partition for the worker's adapter tests.
    /// Returns `(chunk, partition)` so callers can probe both
    /// surfaces.
    fn fixture(
        per_sample_records: Vec<Vec<crate::pileup_record::PileupRecord>>,
        window_start: u32,
        window_end: u32,
        max_group_span: u32,
    ) -> (MaterialisedChunk, WindowPartition) {
        let n = per_sample_records.len();
        let (chunk, _) = loaded_chunk(0, 1..1000, per_sample_records);
        let mut scratch = PartitionScratch::with_n_samples(n);
        let mut partition = WindowPartition::empty();
        partition_window(
            &chunk,
            &(window_start..window_end),
            &[],
            max_group_span,
            &mut scratch,
            &mut partition,
        )
        .unwrap();
        (chunk, partition)
    }

    #[test]
    fn adapter_carries_chrom_id_start_end_from_partition() {
        // MNP at 100 (ref_span=10 → reach 109) + SNP at 105.
        // 105 ≤ 109 → joins → one group spanning [100, 109].
        let s0 = vec![
            record(
                100,
                vec![
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"AAAAAAAAAA",
                        3,
                        -1.0,
                        &[],
                    ),
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"AAACAAAAAA",
                        4,
                        -1.0,
                        &[],
                    ),
                ],
            ),
            record(105, ref_plus_alt(3, 4)),
        ];
        let (mut chunk, partition) = fixture(vec![s0], 1, 200, 100);
        chunk.chrom_id = 9;
        let group = build_overlapping_variant_group(&chunk, &partition, 0, chunk.n_samples(), 9);
        assert_eq!(group.chrom_id, 9);
        assert_eq!(group.start, 100);
        assert_eq!(group.end, 109);
        assert_eq!(group.records.len(), 2);
        assert_eq!(group.records[0].pos, 100);
        assert_eq!(group.records[1].pos, 105);
    }

    #[test]
    fn adapter_emits_one_per_sample_slot_per_position_with_some_for_carriers() {
        // Two samples, both at position 50; only sample 1 at
        // position 60.
        let s0 = vec![record(50, ref_plus_alt(2, 3))];
        let s1 = vec![
            record(50, ref_plus_alt(2, 3)),
            record(60, ref_plus_alt(2, 3)),
        ];
        let (chunk, partition) = fixture(vec![s0, s1], 1, 100, 5);
        let n_samples = chunk.n_samples();

        let group_0 = build_overlapping_variant_group(&chunk, &partition, 0, n_samples, 0);
        assert_eq!(group_0.records.len(), 1);
        assert_eq!(group_0.records[0].per_sample.len(), 2);
        assert!(group_0.records[0].per_sample[0].is_some());
        assert!(group_0.records[0].per_sample[1].is_some());

        let group_1 = build_overlapping_variant_group(&chunk, &partition, 1, n_samples, 0);
        assert_eq!(group_1.records.len(), 1);
        assert_eq!(group_1.records[0].per_sample.len(), 2);
        assert!(group_1.records[0].per_sample[0].is_none());
        assert!(group_1.records[0].per_sample[1].is_some());
    }

    #[test]
    fn adapter_round_trips_payload_via_materialised_record() {
        // Sample 0 has a SNP at 50 (REF=A obs=3, ALT=T obs=4). The
        // adapter should materialise the same PileupRecord that
        // `chunk.materialise_record(0, 0)` returns.
        let s0 = vec![record(50, ref_plus_alt(3, 4))];
        let (chunk, partition) = fixture(vec![s0], 1, 100, 5);
        let group = build_overlapping_variant_group(&chunk, &partition, 0, chunk.n_samples(), 0);

        let expected = chunk.materialise_record(0, 0);
        let actual = group.records[0].per_sample[0].as_ref().unwrap();
        assert_eq!(actual, &expected);
    }

    #[test]
    fn run_window_with_empty_partition_clears_output_and_returns_ok() {
        // A real ChromRefFetcher isn't needed for the empty-partition
        // fast path. Build a dangling Arc that would panic on use,
        // and confirm the function never reaches it.
        let s0 = vec![record(50, ref_plus_alt(2, 3))];
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        let partition = WindowPartition::empty();
        let mut output: Vec<PosteriorRecord> = Vec::new();
        let mut scratch = ColumnarPipelineScratch::empty();

        let fetcher: SharedRefFetcher = Arc::new(NeverCalledFetcher);
        let result = run_window(
            &chunk,
            &partition,
            fetcher,
            PerGroupMergerConfig::default(),
            PosteriorEngineConfig::default(),
            &mut scratch,
            &mut output,
        );
        assert!(result.is_ok());
        assert!(output.is_empty());
    }

    /// In-memory `ChromRefFetcher` for the column-native unify
    /// integration test. Mirrors the `MockRef` inside the kernel
    /// tests.
    #[derive(Clone)]
    struct MockRef {
        seq: Vec<u8>,
        /// 1-based position of `seq[0]`.
        base_offset: u32,
    }

    impl crate::fasta::fetcher::sealed::Sealed for MockRef {}

    impl crate::fasta::fetcher::ChromRefFetcher for MockRef {
        fn length(&self) -> u32 {
            self.base_offset.saturating_sub(1) + self.seq.len() as u32
        }
        fn fetch(
            &self,
            start_1based: u32,
            length: u32,
        ) -> Result<Vec<u8>, crate::fasta::fetcher::ChromRefFetchError> {
            let start_idx = (start_1based - self.base_offset) as usize;
            Ok(self.seq[start_idx..start_idx + length as usize].to_vec())
        }
        fn iter_bases<'a>(
            &'a self,
        ) -> Result<
            Box<dyn Iterator<Item = Result<u8, crate::fasta::fetcher::ChromRefFetchError>> + 'a>,
            crate::fasta::fetcher::ChromRefFetchError,
        > {
            Ok(Box::new(self.seq.iter().copied().map(Ok)))
        }
    }

    #[test]
    fn column_native_unify_matches_row_shape_on_worker_fixture() {
        // Same multi-position fixture as the existing adapter test:
        // sample 0 has a 10-bp MNP at pos 100 reaching to 109, plus
        // an SNP at pos 105. The kernel + the row-shape PerGroupMerger
        // should produce identical (seq, is_compound, constituents)
        // tuples.
        let s0 = vec![
            record(
                100,
                vec![
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"AAAAAAAAAA",
                        3,
                        -1.0,
                        &[],
                    ),
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"AAACAAAAAA",
                        4,
                        -1.0,
                        &[],
                    ),
                ],
            ),
            record(105, ref_plus_alt(3, 4)),
        ];
        let (chunk, partition) = fixture(vec![s0], 1, 200, 100);
        assert_eq!(partition.n_groups(), 1);

        let ref_full = vec![b'A'; 200];
        let fetcher = MockRef {
            seq: ref_full.clone(),
            base_offset: 1,
        };

        // Column-native via the new worker-side entry point.
        let mut scratch = UnifyAllelesScratch::new();
        let mut ref_buf = Vec::new();
        let mut cn_out = UnifiedAllelesColumns::empty();
        unified_alleles_for_group_columnar(
            &chunk,
            &partition,
            0,
            &fetcher,
            16,
            &mut scratch,
            &mut ref_buf,
            &mut cn_out,
        )
        .expect("column-native unify succeeded");

        let cn_alleles: Vec<(Vec<u8>, bool, Vec<(usize, usize)>)> = (0..cn_out.n_alleles())
            .map(|i| {
                let lo = cn_out.constituent_offsets[i] as usize;
                let hi = cn_out.constituent_offsets[i + 1] as usize;
                let constituents: Vec<(usize, usize)> = (lo..hi)
                    .map(|k| {
                        (
                            cn_out.constituent_record_idx[k] as usize,
                            cn_out.constituent_local_allele_idx[k] as usize,
                        )
                    })
                    .collect();
                (
                    cn_out.allele_seq(i).to_vec(),
                    cn_out.is_compound[i],
                    constituents,
                )
            })
            .collect();

        // Row-shape oracle via the existing kernel.
        let fetcher_arc: SharedRefFetcher = Arc::new(fetcher);
        let group = build_overlapping_variant_group(&chunk, &partition, 0, chunk.n_samples(), 0);
        let config = PerGroupMergerConfig::new(2, 16, 64, 32).expect("merger config");
        let iter: Vec<Result<OverlappingVariantGroup, GrouperError>> = vec![Ok(group)];
        let mut merger = PerGroupMerger::with_config(iter.into_iter(), fetcher_arc, config);
        let record = merger
            .next()
            .expect("merger yields at least one item")
            .expect("merger succeeded");
        let oracle: Vec<(Vec<u8>, bool, Vec<(usize, usize)>)> = record
            .alleles
            .iter()
            .map(|a| {
                (
                    a.seq.clone(),
                    a.is_compound,
                    a.constituents
                        .iter()
                        .map(|c| (c.record_idx, c.local_allele_idx))
                        .collect(),
                )
            })
            .collect();

        assert_eq!(cn_alleles, oracle);
    }

    /// End-to-end byte-identity: same partition + chunk, fed through
    /// the column-native `run_window` and through a row-shape
    /// `PerGroupMerger` → `PosteriorEngine` pipeline. The emitted
    /// `PosteriorRecord` sequences must compare equal.
    #[test]
    fn run_window_columnar_matches_row_shape_pipeline() {
        // Two samples with a SNP at 50 + a multi-position group at
        // 100..109 with a compound candidate, exercising layers
        // 1 (compounds), 2 (compound projection), 3 (chain-broken
        // path on the sample without the compound), 4 (EM).
        let s0_snp = record(50, ref_plus_alt(3, 4));
        let s0_mnp = record(
            100,
            vec![
                crate::var_calling::cohort_block::test_helpers::allele(
                    b"AAAAAAAAAA",
                    3,
                    -1.0,
                    &[],
                ),
                crate::var_calling::cohort_block::test_helpers::allele(
                    b"AAACAAAAAA",
                    4,
                    -1.0,
                    &[],
                ),
            ],
        );
        let s0_snp_inside = record(105, ref_plus_alt(3, 4));
        let s1_snp = record(50, ref_plus_alt(5, 6));
        let s1_snp_inside = record(105, ref_plus_alt(4, 5));
        let s0 = vec![s0_snp, s0_mnp, s0_snp_inside];
        let s1 = vec![s1_snp, s1_snp_inside];
        let (chunk, partition) = fixture(vec![s0, s1], 1, 200, 100);
        let ref_full = vec![b'A'; 200];
        let fetcher = MockRef {
            seq: ref_full.clone(),
            base_offset: 1,
        };

        // Column-native path.
        let cn_fetcher: SharedRefFetcher = Arc::new(fetcher.clone());
        let mut cn_scratch = ColumnarPipelineScratch::empty();
        let mut cn_out: Vec<PosteriorRecord> = Vec::new();
        let cfg = PerGroupMergerConfig::new(2, 16, 64, 32).expect("merger config");
        run_window(
            &chunk,
            &partition,
            cn_fetcher,
            cfg,
            PosteriorEngineConfig::default(),
            &mut cn_scratch,
            &mut cn_out,
        )
        .expect("columnar run_window succeeded");

        // Row-shape oracle path (Phase A.0).
        let row_fetcher: SharedRefFetcher = Arc::new(fetcher);
        let mut row_out: Vec<PosteriorRecord> = Vec::new();
        let groups: Vec<OverlappingVariantGroup> = (0..partition.n_groups())
            .map(|g| {
                build_overlapping_variant_group(
                    &chunk,
                    &partition,
                    g,
                    chunk.n_samples(),
                    chunk.chrom_id,
                )
            })
            .collect();
        let merger = PerGroupMerger::with_config(
            groups.into_iter().map(Ok::<_, GrouperError>),
            row_fetcher,
            cfg,
        );
        let engine = PosteriorEngine::with_config(merger, PosteriorEngineConfig::default());
        for record in engine {
            row_out.push(record.expect("row-shape pipeline succeeded"));
        }

        assert_eq!(cn_out.len(), row_out.len());
        for (cn, row) in cn_out.iter().zip(row_out.iter()) {
            assert_eq!(cn, row);
        }
    }

    /// Regression test for the cap_protected over-protection bug
    /// that bit on the 3-tomato real-data run before commit (this
    /// file). The pre-fix layer 1.2 admit-compound code marked an
    /// existing per-position match as `cap_protected = true`,
    /// shrinking the prunable pool below what
    /// [`enforce_max_alleles_columnar`] needs — the cap couldn't
    /// trim alleles down to `max_alleles`, yielding emitted records
    /// with more than `max_alleles` ALTs on the columnar path while
    /// the row-shape path emitted only `max_alleles - 1` ALTs.
    /// Byte-identity against the row-shape `PerGroupMerger` oracle
    /// on a fixture that triggers the cap with both compound + bare
    /// per-position alleles present.
    #[test]
    fn run_window_columnar_matches_row_shape_under_cap_with_compounds() {
        // Group span [100, 101]: two cohort positions joined into one
        // group. Many ALTs at pos 100 to force the cap; chain-linked
        // alleles across pos 100 + pos 101 form compound candidates.
        let s0 = vec![
            record(
                100,
                vec![
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"A", 10, -1.0, &[],
                    ),
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"T", 6, -1.0, &[101],
                    ),
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"G", 2, -1.0, &[],
                    ),
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"C", 1, -1.0, &[],
                    ),
                ],
            ),
            record(
                101,
                vec![
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"A", 10, -1.0, &[],
                    ),
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"G", 5, -1.0, &[101],
                    ),
                ],
            ),
        ];
        let s1 = vec![
            record(
                100,
                vec![
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"A", 8, -1.0, &[],
                    ),
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"G", 3, -1.0, &[],
                    ),
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"T", 1, -1.0, &[],
                    ),
                ],
            ),
            record(
                101,
                vec![crate::var_calling::cohort_block::test_helpers::allele(
                    b"A", 10, -1.0, &[],
                )],
            ),
        ];
        let (chunk, partition) = fixture(vec![s0, s1], 1, 200, 5);

        let ref_full = vec![b'A'; 200];
        let fetcher = MockRef {
            seq: ref_full.clone(),
            base_offset: 1,
        };

        // `max_alleles = 3` forces the cap; compounds are present.
        let cfg = PerGroupMergerConfig::new(2, 3, 64, 32).expect("merger config");

        let cn_fetcher: SharedRefFetcher = Arc::new(fetcher.clone());
        let mut cn_scratch = ColumnarPipelineScratch::empty();
        let mut cn_out: Vec<PosteriorRecord> = Vec::new();
        run_window(
            &chunk,
            &partition,
            cn_fetcher,
            cfg,
            PosteriorEngineConfig::default(),
            &mut cn_scratch,
            &mut cn_out,
        )
        .expect("columnar run_window succeeded");

        let row_fetcher: SharedRefFetcher = Arc::new(fetcher);
        let mut row_out: Vec<PosteriorRecord> = Vec::new();
        let groups: Vec<OverlappingVariantGroup> = (0..partition.n_groups())
            .map(|g| {
                build_overlapping_variant_group(
                    &chunk,
                    &partition,
                    g,
                    chunk.n_samples(),
                    chunk.chrom_id,
                )
            })
            .collect();
        let merger = PerGroupMerger::with_config(
            groups.into_iter().map(Ok::<_, GrouperError>),
            row_fetcher,
            cfg,
        );
        let engine = PosteriorEngine::with_config(merger, PosteriorEngineConfig::default());
        for record in engine {
            row_out.push(record.expect("row-shape pipeline succeeded"));
        }

        assert_eq!(cn_out.len(), row_out.len());
        for (cn, row) in cn_out.iter().zip(row_out.iter()) {
            assert!(
                cn.alleles.len() <= 3,
                "max_alleles cap should hold ({} alleles emitted)",
                cn.alleles.len(),
            );
            assert_eq!(cn, row);
        }
    }

    /// A ChromRefFetcher whose every method panics — the
    /// empty-partition fast path must never invoke it.
    struct NeverCalledFetcher;

    impl crate::fasta::fetcher::sealed::Sealed for NeverCalledFetcher {}

    impl crate::fasta::fetcher::ChromRefFetcher for NeverCalledFetcher {
        fn length(&self) -> u32 {
            unreachable!("empty-partition path must not query the fetcher");
        }
        fn fetch(
            &self,
            _start_1based: u32,
            _length: u32,
        ) -> Result<Vec<u8>, crate::fasta::fetcher::ChromRefFetchError> {
            unreachable!("empty-partition path must not query the fetcher");
        }
        fn iter_bases<'a>(
            &'a self,
        ) -> Result<
            Box<dyn Iterator<Item = Result<u8, crate::fasta::fetcher::ChromRefFetchError>> + 'a>,
            crate::fasta::fetcher::ChromRefFetchError,
        > {
            unreachable!("empty-partition path must not query the fetcher");
        }
    }
}
