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

use crate::paralog::SampleObservation;
use crate::pileup_record::{AlleleSupportStats, PileupRecord};
use crate::psp::PspReadError;
use crate::var_calling::paralog_filter::calibrate::ParalogScoringContext;
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
    CallStats, CalledChunk, CohortPileupRecord, LocusWindowCoverage, RawCohortChunk, RefSpan,
    Variant,
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
    /// The cohort-constant hidden-paralog scoring context when the filter is on;
    /// `None` when it is off. Read-only (so `&VariantCaller` stays `Sync`), it
    /// lets [`call_chunk`](Self::call_chunk) score each locus's LR inline —
    /// scoring runs on the parallel workers instead of a serial calibrate/write
    /// pass (arch `hidden_paralog_inline_scoring.md`).
    paralog: Option<ParalogScoringContext>,
}

impl VariantCaller {
    pub(crate) fn new(
        grouper_cfg: GrouperConfig,
        merger_cfg: PerGroupMergerConfig,
        posterior_cfg: PosteriorEngineConfig,
        min_alt_obs: u32,
        sample_names: Vec<String>,
        paralog: Option<ParalogScoringContext>,
    ) -> Self {
        let genotype_tables = build_genotype_tables(&merger_cfg);
        Self {
            grouper_cfg,
            merger_cfg,
            posterior_cfg,
            min_alt_obs,
            sample_names,
            genotype_tables,
            paralog,
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
        let mut called = self.call_records(chunk_order, records, ref_span)?;
        // Gather each called locus's per-sample window coverage from the
        // (still-owned) per-sample chunks, then score its hidden-paralog LR
        // inline — the worker already holds each locus's window coverage (just
        // gathered) + its AD (the called record), and the score depends on
        // nothing genome-wide, so this runs the expensive transcendental once, in
        // parallel, instead of a serial calibrate + write recompute (arch: inline
        // scoring). Both stay empty when the filter is off.
        called.window_coverage = gather_window_coverage(&per_sample, &called.records);
        called.paralog_lr = self.score_paralog_lrs(&called.records, &called.window_coverage);
        Ok(called)
    }

    /// Score each called locus's hidden-paralog LR from its record + gathered
    /// window coverage, one `f64` per record in `records` order (`f64::NAN` =
    /// unscored → kept). Returns an **empty** vec when the filter is off
    /// (`self.paralog` is `None`). Split from [`call_chunk`](Self::call_chunk) so
    /// the worker-scoring path is unit-testable without materialising a full
    /// [`RawCohortChunk`]; `obs_buf` is a scratch buffer reused across the chunk's
    /// loci (allocation-free per locus, memory `O(samples)`).
    fn score_paralog_lrs(
        &self,
        records: &[Variant],
        window_coverage: &[LocusWindowCoverage],
    ) -> Vec<f64> {
        let Some(paralog) = &self.paralog else {
            return Vec::new();
        };
        let mut obs_buf: Vec<Option<SampleObservation>> = Vec::new();
        records
            .iter()
            .zip(window_coverage)
            .map(|(record, window)| paralog.score(record, window, &mut obs_buf))
            .collect()
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
            // Both filled by `call_chunk` (which has the per-sample chunks to
            // gather from + the scoring context); left empty on this direct path
            // — see the field docs.
            window_coverage: Vec::new(),
            paralog_lr: Vec::new(),
            stats,
        })
    }
}

/// Gather each called locus's per-sample centred-window coverage from the
/// source chunk's per-sample windowed columns, one [`LocusWindowCoverage`] per
/// record in `records` order.
///
/// Keys each locus on its anchor `record.locus.start`; the join itself is a
/// per-sample forward merge (see [`gather_window_coverage_at_anchors`]).
fn gather_window_coverage(
    per_sample: &[SamplePspChunk],
    records: &[Variant],
) -> Vec<LocusWindowCoverage> {
    let anchors: Vec<u32> = records.iter().map(|r| r.locus.start).collect();
    gather_window_coverage_at_anchors(per_sample, &anchors)
}

/// Per-sample centred-window coverage for a list of locus `anchors`, in anchor
/// order, via a **forward per-sample merge-join**.
///
/// Both `anchors` (called-record positions) and each sample's `positions()`
/// (compacted rows) are coordinate-sorted ascending, so one monotone cursor per
/// sample matches every anchor in a single linear sweep — `O(anchors + positions)`
/// per sample, replacing the per-anchor `O(log positions)` binary search that
/// re-probed each sample from scratch for every locus. A hit copies that sample's
/// windowed GC + coverage; a miss leaves `NaN` (the sample has no covered record
/// there, so the score skips it — the "keep, don't flag" behaviour the retired
/// window join gave via `None`). Per-position by construction: the value is the
/// window centred on this exact locus, not a 500 bp tile mean.
///
/// Bit-identical to the old per-anchor binary search: `positions` are unique and
/// sorted, so the cursor lands on the same row index a binary search would find
/// (and a miss yields `NaN` on both). The `debug_assert`s pin the ascending-input
/// precondition the linear sweep relies on (both hold by construction — the
/// producer emits rows in coordinate order and the caller emits loci in genomic
/// order).
///
/// The windowed vectors are aligned 1:1 with `positions` by every chunk
/// constructor, so a hit's index is in bounds — but the reads go through `.get()`
/// so a hypothetical length skew degrades to `NaN` (sample absent → skipped),
/// never an out-of-bounds panic.
fn gather_window_coverage_at_anchors(
    per_sample: &[SamplePspChunk],
    anchors: &[u32],
) -> Vec<LocusWindowCoverage> {
    debug_assert!(
        anchors.windows(2).all(|w| w[0] <= w[1]),
        "anchors must be ascending for the forward merge-join"
    );
    let n_samples = per_sample.len();
    let mut out: Vec<LocusWindowCoverage> = anchors
        .iter()
        .map(|_| LocusWindowCoverage {
            gc: vec![f32::NAN; n_samples],
            coverage: vec![f32::NAN; n_samples],
        })
        .collect();
    for (s, chunk) in per_sample.iter().enumerate() {
        let positions = chunk.positions();
        debug_assert!(
            positions.windows(2).all(|w| w[0] <= w[1]),
            "sample positions must be ascending for the forward merge-join"
        );
        let gc = chunk.windowed_gc();
        let coverage = chunk.windowed_coverage();
        let mut cursor = 0usize;
        for (a, &anchor) in anchors.iter().enumerate() {
            while cursor < positions.len() && positions[cursor] < anchor {
                cursor += 1;
            }
            if cursor < positions.len() && positions[cursor] == anchor {
                out[a].gc[s] = gc.get(cursor).copied().unwrap_or(f32::NAN);
                out[a].coverage[s] = coverage.get(cursor).copied().unwrap_or(f32::NAN);
            }
        }
    }
    out
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
            None,
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

    /// The window-coverage gather aligns each sample's windowed values to a
    /// locus by exact position, leaves `NaN` where a sample has no record
    /// there, and indexes by cohort sample order.
    #[test]
    fn window_coverage_gather_aligns_by_position() {
        // Sample 0 covers 10, 25, 40; sample 1 covers 10, 40 (no record at 25).
        let per_sample = vec![
            SamplePspChunk::from_windowed_for_test(
                0,
                vec![10, 25, 40],
                vec![0.30, 0.31, 0.32],
                vec![100.0, 101.0, 102.0],
            ),
            SamplePspChunk::from_windowed_for_test(
                0,
                vec![10, 40],
                vec![0.40, 0.42],
                vec![200.0, 202.0],
            ),
        ];

        // One merge-join sweep over three ascending anchors (25, 40, 99).
        let got = gather_window_coverage_at_anchors(&per_sample, &[25, 40, 99]);

        // Locus at 25: sample 0 present, sample 1 absent (NaN).
        let at25 = &got[0];
        assert_eq!(at25.gc[0].to_bits(), 0.31f32.to_bits());
        assert_eq!(at25.coverage[0].to_bits(), 101.0f32.to_bits());
        assert!(at25.gc[1].is_nan(), "sample 1 has no record at 25");
        assert!(at25.coverage[1].is_nan());

        // Locus at 40: both present (the cursor advanced past 25 for sample 0
        // and past 10 for sample 1).
        let at40 = &got[1];
        assert_eq!(at40.gc[1].to_bits(), 0.42f32.to_bits());
        assert_eq!(at40.coverage[1].to_bits(), 202.0f32.to_bits());
        assert_eq!(at40.coverage[0].to_bits(), 102.0f32.to_bits());

        // A locus no sample covers → all NaN.
        let at99 = &got[2];
        assert!(at99.gc.iter().all(|v| v.is_nan()));
        assert!(at99.coverage.iter().all(|v| v.is_nan()));
    }

    /// The merge-join cursor handles an anchor that falls **before** a sample's
    /// first covered position and **between** two covered positions (a gap),
    /// leaving `NaN` at both without desyncing later matches — the cases a
    /// per-anchor binary search got for free but the linear cursor must get right.
    #[test]
    fn window_coverage_gather_cursor_handles_gaps_and_leading_miss() {
        // Sample 0 covers 20, 30, 40 (nothing at or before 10, a gap at 25).
        let per_sample = vec![SamplePspChunk::from_windowed_for_test(
            0,
            vec![20, 30, 40],
            vec![0.20, 0.30, 0.40],
            vec![200.0, 300.0, 400.0],
        )];
        // Anchors: 10 (before first pos → miss), 25 (in the 20..30 gap → miss),
        // 30 (hit), 40 (hit).
        let got = gather_window_coverage_at_anchors(&per_sample, &[10, 25, 30, 40]);
        assert!(
            got[0].coverage[0].is_nan(),
            "10 is before the first position"
        );
        assert!(got[1].coverage[0].is_nan(), "25 is in the 20..30 gap");
        assert_eq!(got[2].coverage[0].to_bits(), 300.0f32.to_bits());
        assert_eq!(got[3].coverage[0].to_bits(), 400.0f32.to_bits());
    }

    /// The gather keys on the locus **anchor** (`record.locus.start`), not the
    /// group span: a sample covered only at a grouped variant's interior — with
    /// no record at the anchor — is left `NaN` (absent → skipped, the safe
    /// keep-don't-flag direction). Pins the per-position anchor semantics so a
    /// future change to span-matching is a deliberate, test-visible choice.
    #[test]
    fn window_coverage_gather_keys_on_anchor_not_group_span() {
        // Sample 0 covers 10 and 11; sample 1 covers only 11 (no row at 10).
        let per_sample = vec![
            SamplePspChunk::from_windowed_for_test(
                0,
                vec![10, 11],
                vec![0.30, 0.31],
                vec![100.0, 101.0],
            ),
            SamplePspChunk::from_windowed_for_test(0, vec![11], vec![0.40], vec![200.0]),
        ];
        // A variant anchored at 10 (a group's leftmost position): sample 1 has no
        // covered record there, so it is absent even though it is covered at 11.
        let at10 = &gather_window_coverage_at_anchors(&per_sample, &[10])[0];
        assert_eq!(at10.gc[0].to_bits(), 0.30f32.to_bits());
        assert!(
            at10.gc[1].is_nan() && at10.coverage[1].is_nan(),
            "sample 1 is not covered at the anchor 10, so it is absent"
        );
    }

    /// A chunk from a legacy `.psp` without the windowed columns carries NaN-
    /// filled windowed vectors (still aligned 1:1 with positions). A hit yields a
    /// `NaN` value (sample skipped), never a panic or a search miss.
    #[test]
    fn window_coverage_gather_nan_for_legacy_windowless_chunk() {
        let chunk =
            SamplePspChunk::from_windowed_for_test(0, vec![10], vec![f32::NAN], vec![f32::NAN]);
        let at10 = &gather_window_coverage_at_anchors(&[chunk], &[10])[0];
        assert!(at10.gc[0].is_nan() && at10.coverage[0].is_nan());
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

    /// The worker-scoring path (`score_paralog_lrs`, the fill run by `call_chunk`)
    /// produces one LR per record in record order, stores a finite LR for a
    /// scorable biallelic-SNP locus and `NaN` for an unscored one (no usable
    /// window), and each stored LR bit-equals a standalone
    /// `ParalogScoringContext::score` — the load-bearing M7 mapping the
    /// byte-identity guarantee rests on. `None`-context returns an empty vec.
    #[test]
    fn score_paralog_lrs_fills_one_lr_per_record() {
        use crate::paralog::ParalogModelParams;
        use crate::var_calling::paralog_filter::calibrate::ParalogScoringContext;
        use crate::var_calling::paralog_filter::test_support::{
            prepass, snp_record, window_coverage,
        };

        let n = 4;
        let f = 0.02;
        let ctx = || ParalogScoringContext::new(prepass(n), f, &ParalogModelParams::default());
        let names: Vec<String> = (0..n).map(|s| format!("S{s}")).collect();
        let with_ctx = VariantCaller::new(
            GrouperConfig::default(),
            PerGroupMergerConfig::default(),
            PosteriorEngineConfig::default(),
            0,
            names.clone(),
            Some(ctx()),
        );

        // A scorable ~2× paralog-like locus (finite LR) and an unscored locus
        // (every sample's window coverage NaN → no usable observation).
        let scorable = snp_record(
            10,
            n,
            &|_| (20, 10),
            window_coverage(n, 0.5, &|_| Some(40.0)),
        );
        let unscored = snp_record(20, n, &|_| (20, 10), window_coverage(n, 0.5, &|_| None));
        let records = vec![scorable.record.clone(), unscored.record.clone()];
        let windows = vec![
            scorable.window_coverage.clone(),
            unscored.window_coverage.clone(),
        ];

        let lrs = with_ctx.score_paralog_lrs(&records, &windows);
        assert_eq!(lrs.len(), records.len(), "one LR per record");
        assert!(
            lrs[0].is_finite(),
            "scorable locus → finite LR, got {}",
            lrs[0]
        );
        assert!(lrs[1].is_nan(), "no usable window → unscored NaN");

        // Each stored LR bit-equals a standalone score from an identically-built
        // context (pins worker-score == the pure scorer, bit-for-bit).
        let standalone = ctx();
        let mut buf = Vec::new();
        for (i, (record, window)) in records.iter().zip(&windows).enumerate() {
            let want = standalone.score(record, window, &mut buf);
            assert_eq!(
                lrs[i].to_bits(),
                want.to_bits(),
                "record {i}: fill LR must bit-equal a standalone score",
            );
        }

        // Filter off (`None` context) → empty, so the sink pads to NaN (kept).
        let without = VariantCaller::new(
            GrouperConfig::default(),
            PerGroupMergerConfig::default(),
            PosteriorEngineConfig::default(),
            0,
            names,
            None,
        );
        assert!(without.score_paralog_lrs(&records, &windows).is_empty());
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
