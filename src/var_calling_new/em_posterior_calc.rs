//! Caller step 2 — per-group merge + per-record EM (appendix §D.2), plus the
//! `VariantCaller` worker that composes grouping (§D.1) with the math here.
//!
//! *(today: `var_calling::per_group_merger`, `posterior_engine`, `worker`)*
//!
//! `OverlappingPileupRecords` → `Vec<Variant>`, one group at a time:
//!
//! - per-group merge via the copied
//!   [`per_group_merger`](crate::var_calling_new::per_group_merger) — its
//!   row-shape `MergedRecord` already carries the **flat SoA**
//!   `log_likelihoods` buffer (`[sample * n_genotypes + g]`), so it feeds the
//!   SIMD EM directly through `MergedAllelesView` — **no SIMD loss**;
//! - per-record EM via the copied
//!   [`posterior_engine`](crate::var_calling_new::posterior_engine) — the
//!   `InterpUnivariateSimdMath` lane-of-4 `ln`/`exp` backend, kept verbatim.
//!
//! > **Design note (verified Phase 0).** The columnar `kernels/` chain
//! > (`unify_alleles_columnar` / `project_scalars_columnar` /
//! > `compute_log_likelihoods_columnar`) and its `MaterialisedChunk` /
//! > `WindowPartition` dependencies are **not** transplanted: they were a
//! > byte-identical reimplementation of the row `per_group_merger`, whose
//! > closed-form log-likelihood is scalar in *both* paths (the SIMD is only in
//! > the EM iteration, which fills its lanes from the flat SoA buffer either
//! > path produces). This is the design's §D.2 conditional — "row-shape
//! > `MergedRecord` is fine as long as it hands SoA likelihood buffers to the
//! > EM" — with its predicate confirmed against the code.
//!
//! Phase 3 builds the `VariantCaller` worker here.

use crate::pileup_record::AlleleSupportStats;
use crate::var_calling_new::per_group_merger::{
    MergeGroupOutcome, MergedRecord, PerGroupMergerConfig, PerGroupMergerError,
    build_genotype_tables, merge_group_with_ref,
};
use crate::var_calling_new::pileup_overlaps::overlapping_groups;
use crate::var_calling_new::posterior_engine::{
    PosteriorEngine, PosteriorEngineConfig, PosteriorEngineError,
};
use crate::var_calling_new::types::{CallStats, CalledChunk, PileupCohortChunk, Variant};
use crate::var_calling_new::variant_grouping::{GrouperConfig, GrouperError};

/// Errors the caller can surface (all pathological-input only — the EM math
/// itself is total on well-formed records).
#[derive(Debug, thiserror::Error)]
pub enum CallerError {
    /// Variant grouping failed (e.g. a group exceeding `max_variant_group_span`).
    #[error("grouping failed: {0}")]
    Group(#[from] GrouperError),
    /// Per-group merge failed.
    #[error("per-group merge failed: {0}")]
    Merge(#[from] PerGroupMergerError),
    /// The EM / posterior calculation failed.
    #[error("posterior calculation failed: {0}")]
    Em(#[from] PosteriorEngineError),
}

/// Section 2 — the cohort variant caller (appendix §D). One per worker thread;
/// turns a single [`PileupCohortChunk`] into a [`CalledChunk`] of [`Variant`]s.
///
/// Composes the two record-based steps back-to-back, group by group:
/// 1. **group** the chunk's records into overlapping variant groups
///    ([`overlapping_groups`], §D.1);
/// 2. **merge + call** each group — slice its REF from the chunk's `RefSpan`,
///    [`merge_group_with_ref`] into a `MergedRecord` (SoA likelihoods), apply
///    the pre-EM `min_alt_obs_per_sample` filter, then run the SIMD EM
///    ([`PosteriorEngine`], §D.2) → a [`Variant`].
///
/// The grouped records are created and consumed **here** — they never cross a
/// queue. Every chunk yields exactly one `CalledChunk` (empty `records` when no
/// group survives), keeping `chunk_order` gapless for the writer's reorder.
pub struct VariantCaller {
    grouper_cfg: GrouperConfig,
    merger_cfg: PerGroupMergerConfig,
    posterior_cfg: PosteriorEngineConfig,
    min_alt_obs: u32,
    /// Per-(ploidy, n_alleles) genotype enumeration cache, built once from
    /// `merger_cfg` (identical to what `PerGroupMerger` builds internally).
    genotype_tables: Vec<Vec<Vec<u8>>>,
}

impl VariantCaller {
    pub fn new(
        grouper_cfg: GrouperConfig,
        merger_cfg: PerGroupMergerConfig,
        posterior_cfg: PosteriorEngineConfig,
        min_alt_obs: u32,
    ) -> Self {
        let genotype_tables = build_genotype_tables(&merger_cfg);
        Self {
            grouper_cfg,
            merger_cfg,
            posterior_cfg,
            min_alt_obs,
            genotype_tables,
        }
    }

    /// Call one chunk → its [`CalledChunk`] (records may be empty).
    pub fn call_chunk(&self, chunk: PileupCohortChunk) -> Result<CalledChunk, CallerError> {
        let PileupCohortChunk {
            chunk_order,
            records,
            ref_span,
        } = chunk;
        let mut stats = CallStats::default();

        // Steps 1 + 2a: group, then merge each group (REF sliced from the
        // chunk's span). Surviving `MergedRecord`s feed the EM.
        let mut merged: Vec<MergedRecord> = Vec::new();
        for group in overlapping_groups(records, self.grouper_cfg) {
            let group = group?;
            let ref_slice = ref_span.slice(group.start, group.end);
            match merge_group_with_ref(group, ref_slice, &self.merger_cfg, &self.genotype_tables)? {
                MergeGroupOutcome::Merged(record) => {
                    // Pre-EM `min_alt_obs_per_sample` filter (matches main's
                    // columnar `columnar_passes_min_alt_obs`, applied on the
                    // projected, pre-prune scalars).
                    if passes_min_alt_obs(
                        &record.scalars,
                        record.n_samples,
                        record.alleles.len(),
                        self.min_alt_obs,
                    ) {
                        merged.push(record);
                    } else {
                        stats.records_dropped_low_alt_obs += 1;
                    }
                }
                MergeGroupOutcome::SkippedRefOnly => {
                    stats.groups_skipped_post_unify_ref_only += 1;
                }
                MergeGroupOutcome::SkippedLhCap { n_alleles } => {
                    stats.lh_cap_groups_skipped += 1;
                    stats.lh_cap_alleles_in_skipped += n_alleles as u64;
                }
            }
        }

        // Step 2b: the per-record SIMD EM (InterpUnivariateSimdMath by default).
        let engine = PosteriorEngine::with_config(
            merged
                .into_iter()
                .map(Ok::<MergedRecord, PerGroupMergerError>),
            self.posterior_cfg.clone(),
        );
        let mut called: Vec<Variant> = Vec::new();
        for variant in engine {
            called.push(variant?);
        }

        Ok(CalledChunk {
            chunk_order,
            records: called,
            stats,
        })
    }
}

/// Pre-EM `min_alt_obs_per_sample` predicate over the projected per-allele
/// scalars (`[sample * n_alleles + allele]`). Keep the record iff some ALT
/// column's max `num_obs` across samples reaches `min_obs`. Copied verbatim
/// from `worker::columnar_passes_min_alt_obs` — the byte-identity filter.
fn passes_min_alt_obs(
    scalars: &[AlleleSupportStats],
    n_samples: usize,
    n_alleles: usize,
    min_obs: u32,
) -> bool {
    if min_obs <= 1 || n_alleles <= 1 {
        return true;
    }
    (1..n_alleles).any(|allele_idx| {
        let mut max_obs = 0_u32;
        for sample_idx in 0..n_samples {
            let obs = scalars[sample_idx * n_alleles + allele_idx].num_obs;
            if obs > max_obs {
                max_obs = obs;
            }
        }
        max_obs >= min_obs
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::fetcher::{ChromRefFetchError, ChromRefFetcher};
    use crate::pileup_record::{AlleleObservation, PileupRecord};
    use crate::var_calling::test_helpers::allele;
    use crate::var_calling_new::per_group_merger::{PerGroupMerger, SharedRefFetcher};
    use crate::var_calling_new::pileup_overlaps::overlapping_groups;
    use crate::var_calling_new::types::{CohortPileupRecord, PileupCohortChunk, RefSpan};
    use std::sync::Arc;

    /// In-memory all-`A` reference, base 1 — matches every fixture's REF
    /// allele bytes so unification is trivial. Mirrors the worker tests' MockRef.
    struct MockRef {
        seq: Vec<u8>,
    }
    impl crate::fasta::fetcher::sealed::Sealed for MockRef {}
    impl ChromRefFetcher for MockRef {
        fn length(&self) -> u32 {
            self.seq.len() as u32
        }
        fn fetch(&self, start_1based: u32, length: u32) -> Result<Vec<u8>, ChromRefFetchError> {
            let lo = (start_1based - 1) as usize;
            Ok(self.seq[lo..lo + length as usize].to_vec())
        }
        fn iter_bases<'a>(
            &'a self,
        ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>
        {
            Ok(Box::new(self.seq.iter().copied().map(Ok)))
        }
    }

    /// One cohort position: per-sample `Some(alleles)` or `None` (no record).
    fn cohort_rec(pos: u32, per_sample: Vec<Option<Vec<AlleleObservation>>>) -> CohortPileupRecord {
        CohortPileupRecord {
            chrom_id: 0,
            pos,
            per_sample: per_sample
                .into_iter()
                .map(|opt| opt.map(|alleles| PileupRecord::new(0, pos, alleles)))
                .collect(),
        }
    }

    /// A SNP record: REF `A` + ALT `T`, both with the given obs.
    fn snp(ref_obs: u32, alt_obs: u32) -> Vec<AlleleObservation> {
        vec![
            allele(b"A", ref_obs, -1.0, &[]),
            allele(b"T", alt_obs, -2.0, &[]),
        ]
    }

    /// Two-sample fixture: SNPs at 10, 25 (well-supported) and 40 (singleton
    /// ALT — `min_alt_obs ≥ 2` drops it). Positions are far apart ⇒ one group
    /// each.
    fn fixture() -> Vec<CohortPileupRecord> {
        vec![
            cohort_rec(10, vec![Some(snp(5, 4)), Some(snp(6, 3))]),
            cohort_rec(25, vec![Some(snp(7, 5)), None]),
            cohort_rec(40, vec![Some(snp(9, 1)), Some(snp(8, 1))]),
        ]
    }

    fn chunk_of(records: Vec<CohortPileupRecord>) -> PileupCohortChunk {
        PileupCohortChunk {
            chunk_order: 7,
            records,
            // Over-covering all-A span [1, 200].
            ref_span: RefSpan {
                genomic_start: 1,
                bytes: vec![b'A'; 200],
            },
        }
    }

    fn caller(min_alt_obs: u32) -> VariantCaller {
        VariantCaller::new(
            GrouperConfig::default(),
            PerGroupMergerConfig::default(),
            PosteriorEngineConfig::default(),
            min_alt_obs,
        )
    }

    /// The canonical fetcher-based row pipeline (the independent reference):
    /// same grouper, but `PerGroupMerger` (fetch via MockRef) + `PosteriorEngine`.
    fn row_pipeline(records: Vec<CohortPileupRecord>) -> Vec<Variant> {
        let groups = overlapping_groups(records, GrouperConfig::default());
        let fetcher: SharedRefFetcher = Arc::new(MockRef {
            seq: vec![b'A'; 200],
        });
        let merger = PerGroupMerger::with_config(groups, fetcher, PerGroupMergerConfig::default());
        let engine = PosteriorEngine::with_config(merger, PosteriorEngineConfig::default());
        engine.map(|r| r.expect("row pipeline")).collect()
    }

    #[test]
    fn call_chunk_matches_fetcher_row_pipeline() {
        let records = fixture();
        // min_alt_obs inert (0) so the new caller doesn't filter — the row
        // pipeline (PerGroupMerger) doesn't either; outputs must be identical.
        let got = caller(0).call_chunk(chunk_of(records.clone())).unwrap();
        let want = row_pipeline(records);
        assert_eq!(got.chunk_order, 7);
        assert_eq!(
            got.records, want,
            "RefSpan-sliced caller == fetcher row pipeline"
        );
        assert!(!got.records.is_empty());
    }

    #[test]
    fn min_alt_obs_drops_singletons_and_counts_stat() {
        // Threshold 2: the pos-40 singleton-ALT group is dropped pre-EM.
        let got = caller(2).call_chunk(chunk_of(fixture())).unwrap();
        assert_eq!(got.stats.records_dropped_low_alt_obs, 1);
        let positions: Vec<u32> = got.records.iter().map(|r| r.locus.start).collect();
        assert!(positions.contains(&10) && positions.contains(&25));
        assert!(!positions.contains(&40), "singleton-ALT site filtered out");
        // Inert threshold keeps all three.
        let all = caller(0).call_chunk(chunk_of(fixture())).unwrap();
        assert_eq!(all.records.len(), 3);
        assert_eq!(all.stats.records_dropped_low_alt_obs, 0);
    }

    #[test]
    fn empty_chunk_yields_empty_called_chunk() {
        let got = caller(0).call_chunk(chunk_of(Vec::new())).unwrap();
        assert_eq!(got.chunk_order, 7);
        assert!(got.records.is_empty());
        assert_eq!(got.stats, CallStats::default());
    }

    #[test]
    fn ref_only_positions_skip_without_calling() {
        // Both samples pure-REF (no ALT obs) at a single position ⇒ no variant
        // group is even seeded, so nothing is called and no skip is recorded
        // (the grouper drops pure-REF seeds before the merger sees them).
        let recs = vec![cohort_rec(
            15,
            vec![
                Some(vec![allele(b"A", 5, -1.0, &[])]),
                Some(vec![allele(b"A", 6, -1.0, &[])]),
            ],
        )];
        let got = caller(0).call_chunk(chunk_of(recs)).unwrap();
        assert!(got.records.is_empty());
    }
}
