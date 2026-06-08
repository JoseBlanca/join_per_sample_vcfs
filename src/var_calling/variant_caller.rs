//! Caller section 2 — the `VariantCaller` worker (appendix §D). Composes
//! grouping (§D.1) with the per-group merge + per-record EM (§D.2): one
//! `VariantCaller` per worker thread turns a `RawCohortChunk` into a
//! `CalledChunk`.
//!
//! `OverlappingPileupRecords` → `Vec<Variant>`, one group at a time:
//!
//! - per-group merge via [`per_group_merger`](crate::var_calling::per_group_merger)
//!   — its row-shape `MergedRecord` already carries the **flat SoA**
//!   `log_likelihoods` buffer (`[sample * n_genotypes + g]`), so it feeds the
//!   SIMD EM directly through `MergedAllelesView` — **no SIMD loss**;
//! - per-record EM via [`posterior_engine`](crate::var_calling::posterior_engine)
//!   — the `InterpUnivariateSimdMath` lane-of-4 `ln`/`exp` backend.
//!
//! > **Design note.** There is no separate columnar likelihood path: the row
//! > `per_group_merger`'s closed-form log-likelihood is scalar, and the SIMD
//! > lives only in the EM iteration, which fills its lanes from the flat SoA
//! > buffer the row merger produces. So the record-based path runs the SIMD EM
//! > with no SIMD loss.

use crate::pileup_record::{AlleleSupportStats, PileupRecord};
use crate::psp::PspReadError;
use crate::var_calling::per_group_merger::{
    MergeGroupOutcome, MergedRecord, PerGroupMergerConfig, PerGroupMergerError,
    build_genotype_tables, merge_group_with_ref,
};
use crate::var_calling::per_position_merger::{PerPositionMerger, PerPositionMergerError};
use crate::var_calling::pileup_overlaps::overlapping_groups;
use crate::var_calling::posterior_engine::{
    PosteriorEngine, PosteriorEngineConfig, PosteriorEngineError,
};
use crate::var_calling::sample_reader::SamplePspChunk;
use crate::var_calling::types::{
    CallStats, CalledChunk, CohortPileupRecord, RawCohortChunk, RefSpan, Variant,
};
use crate::var_calling::variant_grouping::{GrouperConfig, GrouperError};

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
    /// The per-position merge (columns→records, moved onto the caller) failed.
    #[error("per-position merge failed: {0}")]
    PerPosition(#[from] PerPositionMergerError),
}

/// Wrap an owned [`PileupRecord`] as the `Result` item the
/// [`PerPositionMerger`] consumes. A named `fn` (not a closure) so every
/// per-sample iterator has the *same* type and they collect into a `Vec<I>`.
pub(crate) fn ok_record(r: PileupRecord) -> Result<PileupRecord, PspReadError> {
    Ok(r)
}

pub(crate) type KeptRecordIter = std::iter::Map<
    std::vec::IntoIter<PileupRecord>,
    fn(PileupRecord) -> Result<PileupRecord, PspReadError>,
>;

/// Reconstruct a chunk's [`CohortPileupRecord`]s from the per-sample compacted
/// columns — the columns→records conversion + per-position merge the caller
/// runs before grouping (the record-building lever runs here, off the single
/// producer thread).
///
/// Reuses the byte-identity-critical [`PerPositionMerger`] verbatim; the
/// per-sample [`records_all`](SamplePspChunk::records_all) feeds it exactly
/// the records the producer's compacted columns hold, in the same sample and
/// position order, so the merged output is identical to a producer-side merge.
/// `sample_names` is merger metadata only (diagnostics) — any cohort-length
/// list yields identical records.
pub(crate) fn merge_compacted_samples(
    per_sample: &[SamplePspChunk],
    sample_names: &[String],
) -> Result<Vec<CohortPileupRecord>, PerPositionMergerError> {
    let iters: Vec<KeptRecordIter> = per_sample
        .iter()
        .map(|c| {
            c.records_all()
                .into_iter()
                .map(ok_record as fn(PileupRecord) -> Result<PileupRecord, PspReadError>)
        })
        .collect();
    let merger = PerPositionMerger::new(iters, sample_names.to_vec(), Vec::new())?;
    let mut records = Vec::new();
    for item in merger {
        let pp = item?;
        records.push(CohortPileupRecord {
            chrom_id: pp.chrom_id,
            pos: pp.pos,
            per_sample: pp.per_sample,
        });
    }
    Ok(records)
}

/// Section 2 — the cohort variant caller (appendix §D). One per worker thread;
/// turns a single [`RawCohortChunk`] into a [`CalledChunk`] of [`Variant`]s.
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
    /// Cohort sample names, in sample order — passed to the per-position
    /// merger (diagnostics metadata) now that the columns→records merge runs
    /// on the caller. Immutable, so `&VariantCaller` stays `Sync`.
    sample_names: Vec<String>,
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
        sample_names: Vec<String>,
    ) -> Self {
        let genotype_tables = build_genotype_tables(&merger_cfg);
        Self {
            grouper_cfg,
            merger_cfg,
            posterior_cfg,
            min_alt_obs,
            sample_names,
            genotype_tables,
        }
    }

    /// Call one chunk → its [`CalledChunk`] (records may be empty).
    pub fn call_chunk(&self, chunk: RawCohortChunk) -> Result<CalledChunk, CallerError> {
        let RawCohortChunk {
            chunk_order,
            per_sample,
            ref_span,
        } = chunk;
        // Columns → records + per-position merge — the record-building work
        // moved off the single producer thread onto this parallel caller.
        let records = merge_compacted_samples(&per_sample, &self.sample_names)?;
        self.call_records(chunk_order, records, ref_span)
    }

    /// Group + per-group merge + EM over already-merged cohort `records`. Split
    /// from [`call_chunk`](Self::call_chunk) so this byte-identity-critical path
    /// can be tested directly against the fetcher row pipeline (the `call_chunk` front-end
    /// only adds the columns→records reconstruction).
    pub fn call_records(
        &self,
        chunk_order: u64,
        records: Vec<CohortPileupRecord>,
        ref_span: RefSpan,
    ) -> Result<CalledChunk, CallerError> {
        let mut stats = CallStats::default();

        // Steps 1 + 2a: group, then merge each group (REF sliced from the
        // chunk's span). Surviving `MergedRecord`s feed the EM.
        let mut merged: Vec<MergedRecord> = Vec::new();
        for group in overlapping_groups(records, self.grouper_cfg) {
            let group = group?;
            let ref_slice = ref_span.slice(group.start, group.end);
            match merge_group_with_ref(group, ref_slice, &self.merger_cfg, &self.genotype_tables)? {
                MergeGroupOutcome::Merged(record) => {
                    // Pre-EM `min_alt_obs_per_sample` filter, applied on the
                    // projected, pre-prune scalars (byte-identical to the
                    // pre-rewrite filter; verified out-of-tree).
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
/// column's max `num_obs` across samples reaches `min_obs`. Byte-identical to
/// the pre-rewrite filter (verified out-of-tree).
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
    use crate::var_calling::per_group_merger::{PerGroupMerger, SharedRefFetcher};
    use crate::var_calling::pileup_overlaps::overlapping_groups;
    use crate::var_calling::test_helpers::allele;
    use crate::var_calling::types::{CohortPileupRecord, RefSpan};
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

    /// Over-covering all-`A` span `[1, 200]` (matches every fixture's REF).
    fn test_ref_span() -> RefSpan {
        RefSpan {
            genomic_start: 1,
            bytes: vec![b'A'; 200],
        }
    }

    fn caller(min_alt_obs: u32) -> VariantCaller {
        VariantCaller::new(
            GrouperConfig::default(),
            PerGroupMergerConfig::default(),
            PosteriorEngineConfig::default(),
            min_alt_obs,
            vec!["S0".to_string(), "S1".to_string()],
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
        let got = caller(0)
            .call_records(7, records.clone(), test_ref_span())
            .unwrap();
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
        let got = caller(2)
            .call_records(7, fixture(), test_ref_span())
            .unwrap();
        assert_eq!(got.stats.records_dropped_low_alt_obs, 1);
        let positions: Vec<u32> = got.records.iter().map(|r| r.locus.start).collect();
        assert!(positions.contains(&10) && positions.contains(&25));
        assert!(!positions.contains(&40), "singleton-ALT site filtered out");
        // Inert threshold keeps all three.
        let all = caller(0)
            .call_records(7, fixture(), test_ref_span())
            .unwrap();
        assert_eq!(all.records.len(), 3);
        assert_eq!(all.stats.records_dropped_low_alt_obs, 0);
    }

    #[test]
    fn empty_chunk_yields_empty_called_chunk() {
        let got = caller(0)
            .call_records(7, Vec::new(), test_ref_span())
            .unwrap();
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
        let got = caller(0).call_records(7, recs, test_ref_span()).unwrap();
        assert!(got.records.is_empty());
    }
}
