//! Phase A.1 layer 3 — per-(sample, genotype) log-likelihood
//! computation.
//!
//! Native column-native port of the row-shape `compute_log_likelihoods`
//! in [`per_group_merger`](crate::var_calling::per_group_merger). For
//! every (sample, genotype) pair the kernel picks one of two paths:
//!
//! - **Standard.** Closed-form Freebayes multinomial likelihood read
//!   straight from the per-sample slice of
//!   [`ProjectedScalarsColumns`].
//! - **Chain-broken-compound fallback.** When the genotype contains a
//!   compound allele the sample is chain-broken at (i.e. its
//!   `chain_anchor_count` for that (allele, sample) is zero), the
//!   likelihood degenerates to a per-constituent-position multinomial
//!   read from the chunk's per-sample columns directly.
//!
//! Output is a flat sample-major `log_likelihoods[sample × n_genotypes]`
//! `Vec<f64>` plus the resolved `n_genotypes` for the group.
//!
//! Byte-identity-tested against the row-shape kernel's
//! `MergedRecord.log_likelihoods` on the same fixture patterns used
//! for layer 1 / layer 2.

use thiserror::Error;

use crate::var_calling::cohort_block::columns::MaterialisedChunk;
use crate::var_calling::cohort_block::kernels::project_scalars::{
    ProjectedScalarsColumns, locate_sample_row_idx, n_local_alleles_at_row,
    read_support_stats_from_columns,
};
use crate::var_calling::cohort_block::kernels::unify_alleles::UnifiedAllelesColumns;
use crate::var_calling::cohort_block::partition::WindowPartition;
use crate::var_calling::per_group_merger::{
    DegeneracyKind, LikelihoodContext, MAX_BITMASK_ALLELES, genotype_order, ln_factorial, xlogy,
};

/// Flat sample-major log-likelihood table produced by
/// [`compute_log_likelihoods_columnar`]. Same logical shape as the
/// row-shape `(log_likelihoods, n_genotypes)` tuple returned by
/// `compute_log_likelihoods`.
///
/// **Layout.**
/// - [`Self::log_likelihoods`] — row-major sample-major. Cell for
///   `(sample_idx, genotype_idx)` is
///   `log_likelihoods[sample_idx * n_genotypes + genotype_idx]`.
///   Length is `n_samples * n_genotypes`.
/// - [`Self::n_samples`] — `unified.n_samples`.
/// - [`Self::n_genotypes`] — number of genotypes resolved from
///   `(ploidy, n_alleles)` against `genotype_tables` (or the
///   `genotype_order` fallback).
// Mi1: `#[non_exhaustive]` — kernel-output columns.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct LogLikelihoodsColumns {
    pub n_samples: usize,
    pub n_genotypes: usize,
    pub log_likelihoods: Vec<f64>,
}

impl LogLikelihoodsColumns {
    pub fn empty() -> Self {
        Self {
            n_samples: 0,
            n_genotypes: 0,
            log_likelihoods: Vec::new(),
        }
    }

    /// Reset to empty while preserving allocated capacity.
    pub fn clear(&mut self) {
        self.n_samples = 0;
        self.n_genotypes = 0;
        self.log_likelihoods.clear();
    }

    /// Per-sample slice into [`Self::log_likelihoods`].
    pub fn row(&self, sample_idx: usize) -> &[f64] {
        let start = sample_idx * self.n_genotypes;
        &self.log_likelihoods[start..start + self.n_genotypes]
    }
}

impl Default for LogLikelihoodsColumns {
    fn default() -> Self {
        Self::empty()
    }
}

/// Errors surfaced by [`compute_log_likelihoods_columnar`].
#[non_exhaustive]
#[derive(Error, Debug, PartialEq)]
pub enum ComputeLogLikelihoodsError {
    /// Closed-form multinomial produced NaN or +∞ for a particular
    /// `(sample, genotype)` pair. Carries the same locus + kind
    /// fields as the row-shape `PerGroupMergerError::DegenerateLikelihood`
    /// variant so the worker can surface the same diagnostic.
    #[error(
        "degenerate likelihood ({kind:?}) at \
         {chrom_id}:{start}..{end} sample {sample_idx} genotype {genotype_idx}"
    )]
    DegenerateLikelihood {
        chrom_id: u32,
        start: u32,
        end: u32,
        sample_idx: usize,
        genotype_idx: usize,
        kind: DegeneracyKind,
    },

    /// The unifier passed in more than `MAX_BITMASK_ALLELES` alleles —
    /// the `--max-alleles-lh-calc` cap should have prevented this.
    /// Signals an upstream invariant break.
    #[error(
        "n_alleles={n_alleles} exceeds MAX_BITMASK_ALLELES={MAX_BITMASK_ALLELES}; \
         max_alleles_lh_calc cap is misconfigured upstream"
    )]
    NAllelesExceedsBitmask { n_alleles: usize },
}

/// Phase A.1 layer 3 — full per-(sample, genotype) log-likelihood
/// computation. Mirrors `compute_log_likelihoods` in
/// [`per_group_merger`](crate::var_calling::per_group_merger) but
/// reads its inputs from the columnar layer 1 / layer 2 outputs and
/// the chunk's per-sample columns.
#[allow(clippy::too_many_arguments)]
pub fn compute_log_likelihoods_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    unified: &UnifiedAllelesColumns,
    projection: &ProjectedScalarsColumns,
    ctx: &LikelihoodContext,
    genotype_tables: &[Vec<Vec<u8>>],
    out: &mut LogLikelihoodsColumns,
) -> Result<(), ComputeLogLikelihoodsError> {
    let n_samples = unified.n_samples;
    let n_alleles = unified.n_alleles();
    let n_kept = projection.n_kept;
    let ploidy = ctx.ploidy;

    if n_alleles > MAX_BITMASK_ALLELES {
        return Err(ComputeLogLikelihoodsError::NAllelesExceedsBitmask { n_alleles });
    }

    // Resolve the genotype list — either from the precomputed cache or
    // via `genotype_order` fallback. Same pattern as the row-shape
    // kernel.
    let fallback;
    let genotypes: &[Vec<u8>] = if n_alleles < genotype_tables.len() {
        &genotype_tables[n_alleles]
    } else {
        fallback = genotype_order(ploidy, n_alleles);
        &fallback
    };
    let n_genotypes = genotypes.len();

    out.clear();
    out.n_samples = n_samples;
    out.n_genotypes = n_genotypes;
    out.log_likelihoods.resize(n_samples * n_genotypes, 0.0);

    let group_position_lo = partition.position_range_for_group(group_idx).start;

    for sample_idx in 0..n_samples {
        let scalars_row = projection.scalars_row(sample_idx);
        let other_scalars = &projection.other_scalars[sample_idx];
        let out_row_start = sample_idx * n_genotypes;

        for (g_idx, genotype) in genotypes.iter().enumerate() {
            // Chain-broken-compound check: pick the first compound
            // allele in `genotype` that this sample is chain-broken
            // at. Columnar equivalent of the row-shape
            // `unified.alleles[a].is_compound && flags_row[a]` test
            // (where `flags_row[a]` was true iff the allele is a
            // compound and the sample's chain_anchor_count is zero).
            let chain_broken_compound = genotype.iter().find_map(|&a| {
                let a = a as usize;
                if unified.is_compound[a] && unified.chain_anchor_count(a, sample_idx) == 0 {
                    Some(a)
                } else {
                    None
                }
            });

            let value = if let Some(compound_idx) = chain_broken_compound {
                chain_broken_log_likelihood_columnar(
                    chunk,
                    partition,
                    group_position_lo,
                    unified,
                    compound_idx,
                    sample_idx,
                    genotype,
                    ploidy,
                )
            } else {
                standard_log_likelihood_columnar(
                    scalars_row,
                    other_scalars,
                    genotype,
                    n_kept,
                    ploidy,
                )
            };

            if value.is_nan() {
                return Err(ComputeLogLikelihoodsError::DegenerateLikelihood {
                    chrom_id: ctx.chrom_id,
                    start: ctx.start,
                    end: ctx.end,
                    sample_idx,
                    genotype_idx: g_idx,
                    kind: DegeneracyKind::NaN,
                });
            }
            if value.is_infinite() && value > 0.0 {
                return Err(ComputeLogLikelihoodsError::DegenerateLikelihood {
                    chrom_id: ctx.chrom_id,
                    start: ctx.start,
                    end: ctx.end,
                    sample_idx,
                    genotype_idx: g_idx,
                    kind: DegeneracyKind::PositiveInfinity,
                });
            }

            out.log_likelihoods[out_row_start + g_idx] = value;
        }
    }

    Ok(())
}

/// Columnar equivalent of the row-shape `standard_log_likelihood`.
/// Reads from the per-sample slice of the layer 2 projection.
fn standard_log_likelihood_columnar(
    scalars: &[crate::pileup_record::AlleleSupportStats],
    other_scalars: &crate::pileup_record::AlleleSupportStats,
    genotype: &[u8],
    n_alleles: usize,
    ploidy: u8,
) -> f64 {
    // The kernel-level guard (`NAllelesExceedsBitmask`) already keeps
    // this invariant; the assert is a release-build backstop matching
    // the row-shape kernel.
    assert!(
        n_alleles <= MAX_BITMASK_ALLELES,
        "standard_log_likelihood_columnar invariant violated: n_alleles={n_alleles} > \
         MAX_BITMASK_ALLELES={MAX_BITMASK_ALLELES}",
    );

    let g_bits: u64 = genotype.iter().fold(0u64, |acc, &a| acc | (1u64 << a));

    let mut p_counts: [u32; MAX_BITMASK_ALLELES] = [0; MAX_BITMASK_ALLELES];
    for &a in genotype {
        p_counts[a as usize] += 1;
    }

    let mut log_l = 0.0;
    for (a, stats) in scalars.iter().enumerate().take(n_alleles) {
        if (g_bits >> a) & 1 == 0 {
            log_l += stats.q_sum;
        }
    }
    log_l += other_scalars.q_sum;

    let mut n_total: u64 = 0;
    let mut sum_ln_n_fact = 0.0;
    let mut sum_n_log_p = 0.0;
    for a in 0..n_alleles {
        if (g_bits >> a) & 1 == 0 {
            continue;
        }
        let n_i = scalars[a].num_obs as u64;
        n_total += n_i;
        sum_ln_n_fact += ln_factorial(n_i);
        let p_i = (p_counts[a] as f64) / (ploidy as f64);
        sum_n_log_p += xlogy(n_i as f64, p_i);
    }

    log_l += ln_factorial(n_total) - sum_ln_n_fact + sum_n_log_p;
    log_l
}

/// Columnar equivalent of the row-shape `chain_broken_log_likelihood`.
/// Walks the compound's constituents, reads each constituent
/// position's per-allele stats directly from the chunk's per-sample
/// columns, and accumulates the per-position multinomial.
#[allow(clippy::too_many_arguments)]
fn chain_broken_log_likelihood_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_position_lo: usize,
    unified: &UnifiedAllelesColumns,
    compound_idx: usize,
    sample_idx: usize,
    genotype: &[u8],
    ploidy: u8,
) -> f64 {
    let constituent_lo = unified.constituent_offsets[compound_idx] as usize;
    let constituent_hi = unified.constituent_offsets[compound_idx + 1] as usize;
    if constituent_lo == constituent_hi {
        return 0.0;
    }

    let compound_slot = compound_idx as u8;
    let n_slots = ploidy as usize;
    let mut log_l = 0.0;

    let sample = &chunk.per_sample[sample_idx];

    for k in constituent_lo..constituent_hi {
        let record_idx_in_group = unified.constituent_record_idx[k] as usize;
        let constituent_local_allele = unified.constituent_local_allele_idx[k] as usize;

        let Some(row_idx) = locate_sample_row_idx(
            partition,
            group_position_lo,
            record_idx_in_group,
            sample_idx,
        ) else {
            // No record from this sample at this constituent position
            // ⇒ zero scalars; the per-position multinomial reduces to
            // 0.0 in log space.
            continue;
        };
        let n_local_alleles = n_local_alleles_at_row(sample, row_idx);

        debug_assert!(
            n_local_alleles <= MAX_BITMASK_ALLELES,
            "chain_broken_log_likelihood_columnar assumes n_local_alleles ≤ \
             {MAX_BITMASK_ALLELES}; got {n_local_alleles}",
        );

        let mut per_pos_counts: [u32; MAX_BITMASK_ALLELES] = [0; MAX_BITMASK_ALLELES];
        for &slot in genotype {
            let local_idx = if slot == compound_slot {
                constituent_local_allele
            } else {
                0 // REF
            };
            per_pos_counts[local_idx] += 1;
        }

        // Error cost at this position: Σ q_sum over local alleles not
        // in per_pos_counts.
        for (a, &count) in per_pos_counts.iter().enumerate().take(n_local_alleles) {
            if count == 0 {
                let stats = read_support_stats_from_columns(sample, row_idx, a);
                log_l += stats.q_sum;
            }
        }

        // Multinomial at this position.
        let mut n_total: u64 = 0;
        let mut sum_ln_n_fact = 0.0;
        let mut sum_n_log_p = 0.0;
        for (a, &count) in per_pos_counts.iter().enumerate().take(n_local_alleles) {
            if count == 0 {
                continue;
            }
            let stats = read_support_stats_from_columns(sample, row_idx, a);
            let n_i = stats.num_obs as u64;
            n_total += n_i;
            sum_ln_n_fact += ln_factorial(n_i);
            let p_i = (count as f64) / (n_slots as f64);
            sum_n_log_p += xlogy(n_i as f64, p_i);
        }
        log_l += ln_factorial(n_total) - sum_ln_n_fact + sum_n_log_p;
    }

    log_l
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::fasta::ChromRefFetcher;
    use crate::fasta::fetcher::ChromRefFetchError;
    use crate::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};
    use crate::var_calling::cohort_block::columns::MaterialisedChunk;
    use crate::var_calling::cohort_block::kernels::project_scalars::{
        ProjectScalarsScratch, ProjectedScalarsColumns, project_scalars_columnar,
    };
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
    ) -> (MaterialisedChunk, WindowPartition, Vec<u8>) {
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

    /// Pull `(MergedRecord.log_likelihoods, .n_genotypes, ploidy)` from
    /// the row-shape kernel — the byte-identity oracle for layer 3.
    fn log_likelihoods_via_existing_kernel(
        chunk: &MaterialisedChunk,
        partition: &WindowPartition,
        group_idx: usize,
        ref_fetcher: SharedRefFetcher,
        max_alleles: usize,
    ) -> (Vec<f64>, usize, u32, u32, u32) {
        let group = crate::var_calling::cohort_block::worker::build_overlapping_variant_group(
            chunk,
            partition,
            group_idx,
            chunk.n_samples(),
            chunk.chrom_id,
        );
        let group_start = partition.group_starts[group_idx];
        let group_end = partition.group_ends[group_idx];
        let chrom_id = chunk.chrom_id;
        let config = PerGroupMergerConfig::new(2, max_alleles, 64, 32).expect("merger config");
        let iter: Vec<Result<OverlappingVariantGroup, GrouperError>> = vec![Ok(group)];
        let mut merger = PerGroupMerger::with_config(iter.into_iter(), ref_fetcher, config);
        let rec = merger
            .next()
            .expect("merger yields at least one item")
            .expect("merger succeeded");
        (
            rec.log_likelihoods,
            rec.n_genotypes,
            chrom_id,
            group_start,
            group_end,
        )
    }

    fn compute_columnar(
        chunk: &MaterialisedChunk,
        partition: &WindowPartition,
        group_idx: usize,
        ref_seq: &[u8],
        max_alleles: usize,
    ) -> LogLikelihoodsColumns {
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

        let mut proj_scratch = ProjectScalarsScratch::new();
        let mut projection = ProjectedScalarsColumns::empty();
        project_scalars_columnar(
            chunk,
            partition,
            group_idx,
            &unified,
            &mut proj_scratch,
            &mut projection,
        )
        .expect("project scalars succeeded");

        let group_start = partition.group_starts[group_idx];
        let group_end = partition.group_ends[group_idx];
        let ctx = LikelihoodContext {
            chrom_id: chunk.chrom_id,
            start: group_start,
            end: group_end,
            ploidy: 2,
        };
        let mut out = LogLikelihoodsColumns::empty();
        compute_log_likelihoods_columnar(
            chunk,
            partition,
            group_idx,
            &unified,
            &projection,
            &ctx,
            &[],
            &mut out,
        )
        .expect("compute log-likelihoods succeeded");
        out
    }

    #[test]
    fn log_likelihoods_no_compound_match_row_shape() {
        let s0 = vec![record(50, ref_plus_alt(3, 4))];
        let s1 = vec![record(50, ref_plus_alt(5, 6))];
        let ref_full = vec![b'A'; 200];
        let (chunk, partition, ref_seq) = group_fixture(vec![s0, s1], &ref_full, 1..200, 50);

        let out = compute_columnar(&chunk, &partition, 0, &ref_seq, 16);
        let (oracle_lls, oracle_n_genotypes, _, _, _) = log_likelihoods_via_existing_kernel(
            &chunk,
            &partition,
            0,
            shared_mock(&ref_full, 1),
            16,
        );

        assert_eq!(out.n_genotypes, oracle_n_genotypes);
        assert_eq!(out.log_likelihoods, oracle_lls);
    }

    #[test]
    fn log_likelihoods_multi_position_no_compound_match_row_shape() {
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

        let out = compute_columnar(&chunk, &partition, 0, &ref_seq, 16);
        let (oracle_lls, oracle_n_genotypes, _, _, _) = log_likelihoods_via_existing_kernel(
            &chunk,
            &partition,
            0,
            shared_mock(&ref_full, 1),
            16,
        );

        assert_eq!(out.n_genotypes, oracle_n_genotypes);
        assert_eq!(out.log_likelihoods, oracle_lls);
    }

    #[test]
    fn log_likelihoods_with_compound_match_row_shape() {
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

        let out = compute_columnar(&chunk, &partition, 0, &ref_seq, 16);
        let (oracle_lls, oracle_n_genotypes, _, _, _) = log_likelihoods_via_existing_kernel(
            &chunk,
            &partition,
            0,
            shared_mock(&ref_full, 1),
            16,
        );

        assert_eq!(out.n_genotypes, oracle_n_genotypes);
        assert_eq!(out.log_likelihoods, oracle_lls);
    }

    #[test]
    fn log_likelihoods_with_cap_drop_match_row_shape() {
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

        let out = compute_columnar(&chunk, &partition, 0, &ref_seq, 2);
        let (oracle_lls, oracle_n_genotypes, _, _, _) = log_likelihoods_via_existing_kernel(
            &chunk,
            &partition,
            0,
            shared_mock(&ref_full, 1),
            2,
        );

        assert_eq!(out.n_genotypes, oracle_n_genotypes);
        assert_eq!(out.log_likelihoods, oracle_lls);
    }
}
