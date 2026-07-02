//! Shared data types for the re-architected cohort `.psp` → VCF pipeline
//! (appendix §C of the re-architecture design).
//!
//! These are the structures that flow between the three runtime sections
//! (producer → caller → writer):
//!
//! - [`CohortPileupRecord`] — one variable cohort position with every
//!   sample's pileup data at it (the producer's emit unit, the grouper's
//!   input). Revived from `per_position_merger::PerPositionPileups`.
//! - [`RawCohortChunk`] — the producer→caller work-unit: a safe-gap-bounded
//!   per-sample columnar slice (`SamplePspChunk`s) + its [`RefSpan`], tagged
//!   `chunk_order`; the caller rebuilds the [`CohortPileupRecord`]s from it.
//! - [`RefSpan`] — the chunk's contiguous REF bytes, fetched once on the
//!   producer (monotonic-forward) and sliced per group by the caller.
//! - [`Variant`] — the final emitted record (today's `PosteriorRecord`).
//! - [`CalledChunk`] — the caller→writer payload (today's `BlockResult`),
//!   carrying [`CallStats`] (today's `WindowRunStats`).

use crate::pileup_record::PileupRecord;
use crate::var_calling::posterior_engine::PosteriorRecord;
use crate::var_calling::sample_reader::SamplePspChunk;

/// One **variable cohort position** with every sample's pileup data at it
/// (appendix §C, `[REVIVE — cf. PerPositionPileups]`).
///
/// The producer emits one of these per position the cohort fold deems
/// variable (AC / `min_alt_obs`, dust applied); the caller's grouper
/// consumes them. Same shape as
/// [`per_position_merger::PerPositionPileups`](crate::var_calling::per_position_merger::PerPositionPileups),
/// but only the variable positions are ever built into one — the ~96 % of
/// non-variant positions are dropped at the per-position merge and never
/// materialised cohort-wide (the memory invariant, appendix §B).
///
/// Each sample's [`PileupRecord`] carries its alleles + per-allele
/// obs/stats, which is enough both for grouping (reach is the REF allele's
/// `ref_span`, derivable from `alleles[0]`) and for the downstream
/// per-group merge + EM.
#[derive(Debug, Clone, PartialEq)]
pub struct CohortPileupRecord {
    pub chrom_id: u32,
    /// 1-based anchor position. Matches [`PileupRecord::pos`].
    pub pos: u32,
    /// Indexed by sample order. Samples without a record at this position
    /// carry `None`. `per_sample.len()` equals the cohort's sample count.
    pub per_sample: Vec<Option<PileupRecord>>,
}

/// The chunk's contiguous REF bytes (appendix §C, `[NEW]`).
///
/// Covers `[genomic_start, genomic_start + bytes.len())` in 1-based contig
/// coordinates, i.e. `bytes[0]` is the base at 1-based `genomic_start`.
/// Fetched **once per chunk on the producer** with a single
/// monotonic-forward read (the only thread that may, since
/// `StreamingChromRefFetcher` panics on non-monotonic access); the parallel
/// callers then [`slice`](Self::slice) per group by offset arithmetic — no
/// halo, because the safe-gap chunk cut makes every group's reach fall
/// inside the chunk's own span. Replaces today's per-group `Vec<Vec<u8>>`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RefSpan {
    /// 1-based contig coordinate of `bytes[0]`.
    pub genomic_start: u32,
    pub bytes: Vec<u8>,
}

impl RefSpan {
    /// The REF bases over the 1-based **inclusive** span
    /// `[group_start, group_end]` — exactly the per-group REF slice the
    /// caller's allele unification reads, with no padding (matching the old
    /// `prefetch_window_ref_bytes` group span `group_end - group_start + 1`).
    ///
    /// # Panics
    ///
    /// Panics if the requested span is not fully inside
    /// `[genomic_start, genomic_start + bytes.len())`.
    pub fn slice(&self, group_start: u32, group_end: u32) -> &[u8] {
        debug_assert!(
            group_start >= self.genomic_start,
            "RefSpan::slice: group_start {group_start} precedes span start {}",
            self.genomic_start,
        );
        let lo = (group_start - self.genomic_start) as usize;
        let hi = (group_end - self.genomic_start) as usize + 1;
        &self.bytes[lo..hi]
    }
}

/// The producer→caller **work-unit** (appendix §C, `[NEW]` — replaces the
/// columnar `MaterialisedChunk`).
///
/// Carries, per sample, a [`SamplePspChunk`] **compacted to this chunk's
/// variable rows only** (positions the cohort fold kept, dust applied) — the
/// columnar form, *not* yet rebuilt into [`PileupRecord`]s. The expensive
/// columns→records conversion + per-position merge (which used to run on the
/// single producer thread, ~½ its time) is deferred to the parallel caller,
/// which calls [`SamplePspChunk::records_all`] + the per-position merger to
/// reconstruct the [`CohortPileupRecord`]s before grouping. Chunk boundaries
/// are still the producer's safe-gap cuts, so each chunk groups in isolation.
///
/// `per_sample.len()` equals the cohort's sample count; a sample with no
/// variable row in this chunk carries an empty compacted chunk.
#[derive(Debug, Clone, PartialEq)]
pub struct RawCohortChunk {
    /// Monotonic, genomic-order production index stamped by the producer
    /// (the writer's reorder ticket).
    pub chunk_order: u64,
    /// Per sample, in cohort sample order: the chunk's variable rows in
    /// columnar form, ready for [`SamplePspChunk::records_all`].
    pub per_sample: Vec<SamplePspChunk>,
    pub ref_span: RefSpan,
}

/// The final emitted record.
///
/// A type alias onto the kernel output [`PosteriorRecord`] — keeping it the
/// *same* type the SIMD EM produces is how byte-identity is preserved (no
/// conversion at the caller→writer boundary). `Variant` is the name the
/// pipeline uses for it.
pub type Variant = PosteriorRecord;

/// Per-chunk run statistics carried alongside the calls (appendix §E; revival
/// of today's `worker::WindowRunStats`).
///
/// Pure run diagnostics — **not** part of the VCF byte-identity contract;
/// the writer rolls these into the run-level stats. The producer/caller fill
/// them; the four counters mirror `WindowRunStats` field-for-field so the
/// aggregated reporting stays unchanged.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct CallStats {
    /// Groups skipped because `unified.alleles.len() > max_alleles_lh_calc`.
    pub lh_cap_groups_skipped: u64,
    /// Sum of `unified.alleles.len()` across every skipped group.
    pub lh_cap_alleles_in_skipped: u64,
    /// Groups whose post-unification allele set was REF-only (every ALT
    /// absorbed into the OTHER pool by the max-alleles cap, or none survived
    /// the per-position projection dedup).
    pub groups_skipped_post_unify_ref_only: u64,
    /// Groups dropped at the merger→EM boundary because no ALT allele reached
    /// `--min-alt-obs-per-sample`.
    pub records_dropped_low_alt_obs: u64,
}

/// The caller→writer payload (appendix §E, `[RENAME BlockResult]`).
///
/// `records` may be **empty** — a chunk that produced no [`Variant`]s after
/// grouping/merge still ships exactly one `CalledChunk`, to keep
/// `chunk_order` gapless for the writer's reorder (an empty `Vec` does not
/// allocate, and the writer treats it uniformly).
#[derive(Debug, Clone, PartialEq)]
pub struct CalledChunk {
    /// Carried through untouched from the source [`RawCohortChunk`].
    pub chunk_order: u64,
    pub records: Vec<Variant>,
    /// Per-record centred-window coverage for the hidden-paralog score, in the
    /// **same order** as `records` (so `window_coverage[i]` describes
    /// `records[i]`). Gathered by [`VariantCaller::call_chunk`] from the source
    /// chunk's per-sample windowed columns; **empty** for a `CalledChunk`
    /// produced by [`VariantCaller::call_records`] directly (the test path,
    /// which has no per-sample chunks to gather from). The spill sink reads it
    /// alongside each record.
    pub window_coverage: Vec<LocusWindowCoverage>,
    pub stats: CallStats,
}

/// Per-sample centred-window coverage for one called locus, gathered from the
/// psp `windowed_gc`/`windowed_coverage` columns by matching each sample's
/// compacted-chunk row at the locus position.
///
/// Both vectors are indexed by cohort sample order and hold `NaN` where the
/// sample has no covered record at the locus — the score skips those samples,
/// matching the retired window join's `None`. Per-sample GC (not a single
/// shared reference GC) because each sample's window spans **its** covered
/// positions; storing it per sample is what lets var-calling stay window-free.
#[derive(Debug, Clone, Default)]
pub struct LocusWindowCoverage {
    /// Per-sample centred-window GC fraction; `NaN` = sample absent at locus.
    pub gc: Vec<f32>,
    /// Per-sample centred-window mean coverage; `NaN` = sample absent at locus.
    pub coverage: Vec<f32>,
}

/// Byte-identity equality: the windowed values carry `NaN` sentinels, so they
/// are compared on their bit pattern (a derived `PartialEq` would make
/// `NaN != NaN` and break `CalledChunk` equality in tests). No `Eq`/`Hash`.
impl PartialEq for LocusWindowCoverage {
    fn eq(&self, other: &Self) -> bool {
        let Self { gc, coverage } = self;
        gc.len() == other.gc.len()
            && gc
                .iter()
                .zip(&other.gc)
                .all(|(a, b)| a.to_bits() == b.to_bits())
            && coverage.len() == other.coverage.len()
            && coverage
                .iter()
                .zip(&other.coverage)
                .all(|(a, b)| a.to_bits() == b.to_bits())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn span() -> RefSpan {
        // Bases for 1-based coordinates 10..=14.
        RefSpan {
            genomic_start: 10,
            bytes: b"ACGTN".to_vec(),
        }
    }

    #[test]
    fn ref_span_slice_full_span() {
        assert_eq!(span().slice(10, 14), b"ACGTN");
    }

    #[test]
    fn ref_span_slice_single_base_is_inclusive() {
        // [12, 12] is one base — the 'G' at the 3rd position.
        assert_eq!(span().slice(12, 12), b"G");
    }

    #[test]
    fn ref_span_slice_interior_inclusive_range() {
        assert_eq!(span().slice(11, 13), b"CGT");
    }

    #[test]
    fn ref_span_slice_at_span_start() {
        assert_eq!(span().slice(10, 11), b"AC");
    }

    #[test]
    #[should_panic]
    fn ref_span_slice_past_end_panics() {
        // 15 is one past the last base (14).
        let _ = span().slice(14, 15);
    }

    #[test]
    #[should_panic]
    fn ref_span_slice_below_span_start_panics() {
        // group_start 5 precedes the span's genomic_start (10): the documented
        // precondition violation (a `debug_assert!`; in release the lower-bound
        // subtraction would underflow into an out-of-range index).
        let _ = span().slice(5, 12);
    }

    #[test]
    fn default_ref_span_is_zeroed() {
        let e = RefSpan::default();
        assert_eq!(e.genomic_start, 0);
        assert!(e.bytes.is_empty());
    }
}
