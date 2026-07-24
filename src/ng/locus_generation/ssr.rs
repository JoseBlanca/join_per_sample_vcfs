//! The STR locus generator — one microsatellite tract → one locus.
//!
//! The first [`LocusGenerator`](super::LocusGenerator): it consumes an `SsrSegment`,
//! fetches the reads over the tract, aligns each to read off its repeat, and tallies the
//! answers into a [`SampleLocusObservations`](super::SampleLocusObservations). A port of
//! production `src/ssr/pileup/`, adapted at two seams (the split of coordinates from bases,
//! and a widened admission gate). See `doc/devel/ng/spec/locus_generation_ssr.md` (design)
//! and `doc/devel/ng/arch/locus_generation_ssr.md` (types & interfaces).
//!
//! This module lands across the STR generator plan. So far: the config, counts, cap
//! constant, the working input ([`SsrLocus`] + its margin fetch), the reservoir cap, the
//! per-locus read fetch, and the read-region CIGAR mapping. The delimiter is the alignment
//! module's [`SsrFlatGapAligner`](crate::ng::alignment::ssr_best_path_flat_gap); the classify
//! pipeline, tally, and the `LocusGenerator` impl land next.

use super::LocusGenerationError;
use crate::bam::alignment_input::MappedRead;
use crate::ng::read::input::SampleReads;
use crate::ng::ref_seq::{ContigTable, RawRefSeq, RefSeq, RefSeqError};
use crate::ng::region_typing::segment_criteria::SsrSegment;
use crate::ng::types::{Bp, ContigId, GenomeRegion, Position};

/// An STR locus ready to align against: the segment plus the reference bases the aligner
/// aligns the reads to.
///
/// The ng counterpart of production's `Locus`, **split** so the coordinates come from region
/// typing (`SsrSegment`) and the bases are fetched here (spec §2). Fetching them makes the
/// port an *adaptation, not a lift* — the most likely place for the port to go subtly wrong,
/// which is why the flank lengths are always **measured** from the clamped span, never
/// assumed to be `flank_bp`.
#[derive(Debug, Clone, PartialEq)]
pub struct SsrLocus {
    /// The tract's coordinates, motif and purity — from region typing.
    pub segment: SsrSegment,
    /// The tract plus its query margin, canonical `{A,C,G,T,N}` bases, **clamped at contig
    /// ends** — so this may be shorter than `2 * flank_bp + tract`, and each flank must be
    /// measured, never assumed (spec §2).
    pub tract_with_margin_bases: Box<[u8]>,
    /// 1-based position of `tract_with_margin_bases[0]`.
    pub margin_start: Position,
}

impl SsrLocus {
    /// Fetch the tract ± `flank_bp` for `segment` into an [`SsrLocus`], clamped at the contig
    /// ends, using the reused `buffer` as scratch (spec §2).
    ///
    /// `contig` is the tract's contig (from the region handed to `begin_segment`); its length
    /// is read from `reference`'s contig table to clamp the right margin, and `fetch_into`
    /// validates the final window. The bases are canonical — the aligner compares against
    /// them, and production upper-cases the same way.
    pub fn fetch<R: RefSeq + ContigTable>(
        reference: &R,
        contig: ContigId,
        segment: SsrSegment,
        flank_bp: Bp,
        buffer: &mut Vec<u8>,
    ) -> Result<Self, RefSeqError> {
        let flank = flank_bp.get();
        let tract_start = segment.start(); // 1-based inclusive
        let tract_end = segment.end();
        // The left margin cannot run below base 1; the right cannot run past the contig.
        let margin_start = tract_start.saturating_sub(flank).max(1);
        let margin_end = match reference.contigs().entries.get(contig.get() as usize) {
            // Clamp the right margin to the contig — but only when the tract actually fits
            // inside it. A tract reaching past the contig end is a broken segment: leaving
            // the window unclamped makes `fetch_into` reject it as out-of-bounds, rather than
            // fetching a window that is missing the tract (which would underflow the derived
            // right flank). An unknown contig (`None`) is likewise left for `fetch_into`.
            Some(entry) if tract_end <= entry.length => (tract_end + flank).min(entry.length),
            _ => tract_end + flank,
        };
        let length = (margin_end + 1).saturating_sub(margin_start);

        reference.fetch_into(contig, margin_start, length, buffer)?;
        Ok(Self {
            segment,
            tract_with_margin_bases: buffer.as_slice().into(),
            margin_start: Position(margin_start),
        })
    }

    /// The left flank's length in bases — **measured** from the clamp (`tract_start −
    /// margin_start`), so a tract near the contig start reports a short flank, not `flank_bp`.
    pub fn left_flank_len(&self) -> usize {
        (self.segment.start() - self.margin_start.get()) as usize
    }

    /// The right flank's length in bases — measured as what remains after the left flank and
    /// the tract, so a tract near the contig end reports a short flank.
    pub fn right_flank_len(&self) -> usize {
        self.tract_with_margin_bases.len()
            - self.left_flank_len()
            - self.segment.tract_len() as usize
    }
}

/// The STR generator's knobs — owned and taken at construction
/// ([shared config discipline](super); spec §4). Each generator owns its own knobs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SsrGeneratorConfig {
    /// The flank fetched either side of the tract — the aligner's anchor and the read
    /// query margin. Must be `<= bundle_threshold`, the radius region typing guarantees
    /// repeat-free (checked by [`Self::check_flank_within`]); equal by default (spec §4).
    pub flank_bp: Bp,
    /// Reads kept per locus, reservoir-sampled. `None` = no cap (spec §4).
    pub max_reads_per_locus: Option<u32>,
}

/// ng's own per-locus read cap — **not** production's `MAX_READS_PER_LOCUS`, so the two can
/// diverge. Starts at 1000 (matching production) but is never-measured and soft: to be set
/// by experiment (spec §4).
pub const DEFAULT_SSR_MAX_READS_PER_LOCUS: u32 = 1000;

impl Default for SsrGeneratorConfig {
    /// The flank equals the bundle threshold's default (spec §4, "equal by default"), and
    /// the cap starts at [`DEFAULT_SSR_MAX_READS_PER_LOCUS`].
    fn default() -> Self {
        Self {
            flank_bp: Bp(crate::ng::region_typing::segment_criteria::DEFAULT_BUNDLE_THRESHOLD),
            max_reads_per_locus: Some(DEFAULT_SSR_MAX_READS_PER_LOCUS),
        }
    }
}

impl SsrGeneratorConfig {
    /// Check the cross-config invariant `flank_bp <= bundle_threshold`.
    ///
    /// A wider flank than the radius region typing guarantees repeat-free would let the read
    /// query hit a neighbouring repeat, leaving the aligner's anchor no longer clean (spec
    /// §4). It is a relation between two configs, so no newtype can hold it — the generator's
    /// constructor calls this. `bundle_threshold` is [`SsrSegmentCriteria::bundle_threshold`]
    /// ([`crate::ng::region_typing::segment_criteria`]).
    pub fn check_flank_within(&self, bundle_threshold: Bp) -> Result<(), SsrGeneratorConfigError> {
        if self.flank_bp.get() > bundle_threshold.get() {
            return Err(SsrGeneratorConfigError::FlankExceedsBundleThreshold {
                flank_bp: self.flank_bp.get(),
                bundle_threshold: bundle_threshold.get(),
            });
        }
        Ok(())
    }
}

/// A malformed STR generator configuration. `#[non_exhaustive]`; raised at construction so a
/// bad knob combination never reaches a locus.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum SsrGeneratorConfigError {
    /// The flank is wider than region typing's clean-flank guarantee — the read query could
    /// reach a neighbouring repeat and spoil the aligner's anchor (spec §4).
    #[error(
        "flank_bp ({flank_bp}) exceeds bundle_threshold ({bundle_threshold}): the read query \
         would reach past the repeat-free radius region typing guarantees"
    )]
    FlankExceedsBundleThreshold {
        flank_bp: u64,
        bundle_threshold: u64,
    },
}

/// Run-level STR counts, reported alongside the shared
/// [`LocusCounts`](super::LocusCounts). The locus records *that* reads yielded nothing
/// ([`reads_without_observation`](super::SampleLocusObservations::reads_without_observation));
/// **why** is this generator's to report, because the reasons are specific to how it reads a
/// tract and mean nothing to a pileup (spec §4).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SsrGeneratorCounts {
    /// Reads fetched over the tract-plus-margin query span.
    pub reads_fetched: u64,
    /// Reads the reservoir cap dropped.
    pub reads_discarded_by_cap: u64,
    /// Observations pinning the tract length (both borders seen).
    pub observations_complete: u64,
    /// Observations proving a lower bound (one border, ran off the other end).
    pub observations_partial: u64,
    /// Reads that reached the aligner and anchored no border (`read_preparation_ssr.md` §4).
    pub no_border_anchored: u64,
    /// Reads dropped for low base quality.
    pub low_quality: u64,
    /// Reads whose allele stayed flank-truncated even after widening.
    pub window_truncated: u64,
}

// ---------------------------------------------------------------------
// The per-locus read cap — a faithful port of production's reservoir sampler
// (src/ssr/pileup/fetch_reads.rs), keyed to ng's own seed and cap constant.
// ---------------------------------------------------------------------

/// A tiny deterministic PRNG (SplitMix64) — seeded per locus so the depth-cap subsample is
/// reproducible and thread-count-invariant, with no external RNG whose stream could shift.
/// Ported verbatim from production; the constants are load-bearing for byte parity.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e3779b97f4a7615);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
}

/// Deterministic per-locus reservoir seed from a contig name and a **0-based** tract start:
/// FNV-1a over the name bytes, folded with the start. Ported verbatim from production so the
/// kept set matches byte-for-byte (the parity oracle depends on it).
///
/// **The trap (spec §4):** the seed is over the contig **name** and the **0-based** start.
/// ng speaks `ContigId` and 1-based positions, so seeding from the id or the 1-based start
/// silently produces a *different* kept set — deterministic, so the parity test fails looking
/// like an aligner bug. Callers seed through [`seed_for_segment`], which does the conversion.
fn locus_seed(chrom: &str, start_0based: u32) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x100000001b3;
    let mut h = FNV_OFFSET;
    for &b in chrom.as_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h ^= start_0based as u64;
    h.wrapping_mul(FNV_PRIME)
}

/// The reservoir seed for an STR segment — the one place the seed trap is discharged: the
/// contig **name** and the tract's **0-based** start (`start() - 1`, since `SsrSegment` is
/// 1-based).
///
/// The `as u32` matches production's seed domain (its `Locus.start` is `u32`); per-contig
/// positions are far below `2^32` (the largest chromosome is ~250 Mb), so the cast never
/// truncates in practice and parity holds. `start() - 1` cannot underflow: `SsrSegment::new`
/// enforces `1 <= start`.
pub fn seed_for_segment(segment: &SsrSegment) -> u64 {
    locus_seed(segment.chrom(), (segment.start() - 1) as u32)
}

/// Reservoir sampler (Algorithm R) — an effectively-uniform sample of up to `capacity` items
/// from a stream of unknown length, in one pass with `O(capacity)` memory. Ported verbatim
/// from production: the eviction index `next_u64() % seen` carries a modulo bias bounded by
/// `seen / 2^64` (negligible at any real depth), accepted deliberately because the draw is
/// deterministic and thread-count-invariant — an unbiased reduction would change the kept set
/// and break that. The caller must `offer` admitted reads in a fixed total order —
/// `SampleReads`' merge order (spec §4).
pub struct Reservoir<T> {
    capacity: usize,
    held: Vec<T>,
    /// Admitted items offered so far (the `i` of Algorithm R).
    seen: u64,
    rng: SplitMix64,
}

impl<T> Reservoir<T> {
    /// A reservoir of `capacity` items seeded by `seed` (from [`seed_for_segment`]).
    pub fn new(capacity: usize, seed: u64) -> Self {
        Self {
            capacity,
            held: Vec::with_capacity(capacity),
            seen: 0,
            rng: SplitMix64::new(seed),
        }
    }

    /// Offer one admitted item. Keeps the first `capacity`; for the `i`-th item
    /// (`i > capacity`) keeps it with probability `capacity / i`, evicting one held item
    /// uniformly at random if kept.
    pub fn offer(&mut self, item: T) {
        self.seen += 1;
        if self.held.len() < self.capacity {
            self.held.push(item);
        } else {
            let j = (self.rng.next_u64() % self.seen) as usize;
            if j < self.capacity {
                self.held[j] = item;
            }
        }
    }

    /// The admitted depth — total items offered (the reservoir sees only admitted reads).
    pub fn seen(&self) -> u64 {
        self.seen
    }

    /// Consume the reservoir, yielding the sampled items (≤ `capacity`).
    pub fn into_held(self) -> Vec<T> {
        self.held
    }
}

// ---------------------------------------------------------------------
// The per-locus read fetch (D1): fetch over the tract+margin span, depth-cap.
// ---------------------------------------------------------------------

/// The reads kept for one locus after depth-capping, and how many were fetched — the caller
/// records `reads_discarded_by_cap = fetched - kept.len()`.
pub struct CappedReads {
    pub kept: Vec<MappedRead>,
    pub fetched: u64,
}

/// Fetch the reads over a locus's query span (tract + margin) and depth-cap them with the
/// reservoir, seeded per locus (spec §2, §4).
///
/// Admission is **relevance** — overlap with `query_span`, which `SampleReads` already applies
/// — **not** spanning: production's `reaches_locus` gate is deliberately not ported, because a
/// read that runs off mid-tract is exactly what a partial observation is made of (spec §2). The
/// cap sits between fetch and align; `None` keeps every read. `make_reference` is the per-query
/// reference factory `reads_in_region` needs for the mismatch-fraction filter.
pub fn fetch_capped_reads<R: RawRefSeq>(
    reads: &SampleReads,
    query_span: GenomeRegion,
    seed: u64,
    max_reads_per_locus: Option<u32>,
    make_reference: impl FnMut() -> R,
) -> Result<CappedReads, LocusGenerationError> {
    let stream = reads.reads_in_region(query_span, make_reference)?;
    let mut fetched = 0u64;
    let kept = match max_reads_per_locus {
        Some(cap) => {
            let mut reservoir = Reservoir::new(cap as usize, seed);
            for read in stream {
                reservoir.offer(read?);
                fetched += 1;
            }
            reservoir.into_held()
        }
        None => {
            let mut all = Vec::new();
            for read in stream {
                all.push(read?);
                fetched += 1;
            }
            all
        }
    };
    Ok(CappedReads { kept, fetched })
}

// ---------------------------------------------------------------------
// D2a — the read-region mapping: CIGAR → the read-coordinate slice covering the locus
// window, plus the region quality gate. Ported from production
// (`src/ssr/pileup/footprint.rs`, `alignment.rs`), adapted to take window coordinates
// directly and a plain scratch buffer. Consumed by the classify pipeline (D2b).
// ---------------------------------------------------------------------
mod read_region {
    #![allow(dead_code)] // consumed by D2b (the classify pipeline)

    use crate::pileup::walker::CigarOp;
    use std::ops::Range;

    /// The lower-quartile base-quality floor a delimited tract must clear (production's
    /// `MIN_REGION_Q1`).
    pub(super) const MIN_REGION_Q1: u8 = 15;

    /// A read's reference footprint from its CIGAR — the aligned span (0-based, half-open)
    /// and the soft-clip on each end (the optimistic long-allele reach).
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub(super) struct Footprint {
        pub(super) ref_start: u32,
        pub(super) ref_end: u32,
        pub(super) leading_clip: u32,
        pub(super) trailing_clip: u32,
    }

    /// First soft-clip length at an end (skipping a hard-clip), or `0` if the end is aligned.
    fn end_soft_clip(op: &CigarOp) -> Option<u32> {
        match op {
            CigarOp::HardClip(_) => None,
            CigarOp::SoftClip(n) => Some(*n),
            _ => Some(0),
        }
    }

    /// The read's reference footprint from its CIGAR and 1-based mapping position.
    pub(super) fn read_footprint(cigar: &[CigarOp], pos: u64) -> Footprint {
        let ref_start = pos.saturating_sub(1) as u32;
        let ref_span: u32 = cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n)
                | CigarOp::Deletion(n)
                | CigarOp::Skip(n)
                | CigarOp::SeqMatch(n)
                | CigarOp::SeqMismatch(n) => *n,
                _ => 0,
            })
            .sum();
        Footprint {
            ref_start,
            ref_end: ref_start + ref_span,
            leading_clip: cigar.iter().find_map(end_soft_clip).unwrap_or(0),
            trailing_clip: cigar.iter().rev().find_map(end_soft_clip).unwrap_or(0),
        }
    }

    /// Read coordinate of a 0-based reference position `target` inside the aligned span. A
    /// dual-cursor CIGAR walk; a `target` inside a deletion maps to the read position the
    /// deletion sits at.
    fn ref_to_read(cigar: &[CigarOp], ref_start: u32, leading_clip: u32, target: u32) -> usize {
        let mut ref_cur = ref_start;
        let mut read_cur = leading_clip as usize;
        for op in cigar {
            match op {
                CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                    if target < ref_cur + n {
                        return read_cur + (target - ref_cur) as usize;
                    }
                    ref_cur += n;
                    read_cur += *n as usize;
                }
                CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                    if target < ref_cur + n {
                        return read_cur;
                    }
                    ref_cur += n;
                }
                CigarOp::Insertion(n) => read_cur += *n as usize,
                CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
            }
        }
        read_cur
    }

    /// The read-coordinate span covering the locus window `[window_start, window_start +
    /// window_len)` (0-based reference), where the whole soft-clip is included on a side the
    /// window opens past the aligned span — that is where a long allele's extra tract lives,
    /// and the aligner realigns within the grabbed bases (spec §2).
    pub(super) fn extract_region(
        cigar: &[CigarOp],
        fp: Footprint,
        read_len: usize,
        window_start: u32,
        window_len: u32,
    ) -> Range<usize> {
        let w_end = window_start + window_len;
        let r_start = if window_start < fp.ref_start {
            0 // window opens left of the alignment → take the full leading clip
        } else {
            ref_to_read(cigar, fp.ref_start, fp.leading_clip, window_start)
        };
        let r_end = if w_end > fp.ref_end {
            read_len // window closes right of the alignment → take the full trailing clip
        } else {
            ref_to_read(cigar, fp.ref_start, fp.leading_clip, w_end)
        };
        // Clamp into `[0, read_len]` — defense-in-depth against a length-inconsistent CIGAR.
        let r_start = r_start.min(read_len);
        r_start..r_end.min(read_len).max(r_start)
    }

    /// Whether a delimited tract looks truncated by the extraction window rather than by the
    /// read genuinely ending — the long-allele recovery trigger. A side is window-bounded when
    /// [`extract_region`] did not reach the read edge there; the tract is suspicious when a
    /// window-bounded side carries fewer flank bytes than the locus declares. `tract` is
    /// region-relative.
    pub(super) fn flank_truncated(
        region: &Range<usize>,
        tract: &Range<usize>,
        read_len: usize,
        left_flank_len: usize,
        right_flank_len: usize,
    ) -> bool {
        let region_len = region.end - region.start;
        let left_flank_bytes = tract.start;
        let right_flank_bytes = region_len - tract.end;
        let left_window_bounded = region.start > 0;
        let right_window_bounded = region.end < read_len;
        (left_window_bounded && left_flank_bytes < left_flank_len)
            || (right_window_bounded && right_flank_bytes < right_flank_len)
    }

    /// Widen a read-coordinate region by one full reference flank each side, clamped to the
    /// read — in read, not reference, coordinates, sidestepping the unreliable CIGAR of a
    /// mis-aligned long-allele read.
    pub(super) fn widen_region(
        region: Range<usize>,
        read_len: usize,
        left_flank_len: usize,
        right_flank_len: usize,
    ) -> Range<usize> {
        let start = region.start.saturating_sub(left_flank_len);
        let end = (region.end + right_flank_len).min(read_len);
        start..end
    }

    /// Whether the tract's base qualities clear the gate: the nearest-rank lower quartile (the
    /// element at sorted index `⌊(n-1)/4⌋`) is at least `threshold`. `buffer` is reused
    /// scratch. An empty tract passes vacuously.
    pub(super) fn passes_quality_gate(quals: &[u8], threshold: u8, buffer: &mut Vec<u8>) -> bool {
        if quals.is_empty() {
            return true;
        }
        buffer.clear();
        buffer.extend_from_slice(quals);
        let k = (buffer.len() - 1) / 4;
        let (_, q1, _) = buffer.select_nth_unstable(k);
        *q1 >= threshold
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        /// A window fully inside an all-Match read maps to the same read offsets as reference
        /// offsets (shifted by the mapping position).
        #[test]
        fn extract_region_maps_a_window_inside_an_all_match_read() {
            // read at pos 11 (0-based 10), 30M; window [15, 25) (0-based).
            let cigar = vec![CigarOp::Match(30)];
            let fp = read_footprint(&cigar, 11);
            let region = extract_region(&cigar, fp, 30, 15, 10);
            // ref 15 → read 5 (15 - 10), ref 25 → read 15.
            assert_eq!(region, 5..15);
        }

        /// An internal deletion shifts the window's right edge left in read coordinates.
        #[test]
        fn extract_region_maps_across_an_internal_deletion() {
            // pos 11 (0-based 10): 10M 4D 16M — ref span 30, read len 26. Window [15, 28).
            let cigar = vec![CigarOp::Match(10), CigarOp::Deletion(4), CigarOp::Match(16)];
            let fp = read_footprint(&cigar, 11);
            assert_eq!(fp.ref_end, 40);
            let region = extract_region(&cigar, fp, 26, 15, 13);
            // ref 15 → read 5; ref 28 is past the deletion (ref 20..24 deleted), read = 10 + (28-24) = 14.
            assert_eq!(region, 5..14);
        }

        /// A window opening left of the alignment takes the full leading soft-clip.
        #[test]
        fn extract_region_takes_the_leading_softclip_when_the_window_opens_left() {
            // pos 21 (0-based 20): 5S 20M. Window [18, 30) opens left of ref_start 20.
            let cigar = vec![CigarOp::SoftClip(5), CigarOp::Match(20)];
            let fp = read_footprint(&cigar, 21);
            assert_eq!(fp.leading_clip, 5);
            let region = extract_region(&cigar, fp, 25, 18, 12);
            assert_eq!(region.start, 0, "the full leading clip is grabbed");
        }

        #[test]
        fn flank_truncated_flags_a_window_bounded_short_flank_only() {
            // region [3, 40) within a 50-base read (both sides window-bounded), tract [5, 32)
            // → left flank 5, right flank region_len(37) - 32 = 5.
            assert!(
                flank_truncated(&(3..40), &(5..32), 50, 6, 6),
                "a flank of 5 below the declared 6 is truncated"
            );
            assert!(
                !flank_truncated(&(3..40), &(5..32), 50, 5, 5),
                "flanks exactly the declared length are not truncated"
            );
            // A region reaching both read edges is never flagged (genuine allele ≥ read length).
            assert!(!flank_truncated(&(0..50), &(0..45), 50, 6, 6));
        }

        #[test]
        fn widen_region_extends_a_flank_each_side_clamped_to_the_read() {
            assert_eq!(widen_region(20..40, 100, 10, 10), 10..50);
            // Clamp at both ends.
            assert_eq!(widen_region(5..95, 100, 10, 10), 0..100);
        }

        #[test]
        fn the_quality_gate_keys_on_the_nearest_rank_lower_quartile() {
            let mut buffer = Vec::new();
            assert!(
                passes_quality_gate(&[], MIN_REGION_Q1, &mut buffer),
                "empty passes"
            );
            // len 1 → k = 0 (the lowest); 15 passes, 14 fails.
            assert!(passes_quality_gate(&[15], MIN_REGION_Q1, &mut buffer));
            assert!(!passes_quality_gate(&[14], MIN_REGION_Q1, &mut buffer));
            // len 8 → k = ⌊7/4⌋ = 1 (the 2nd-lowest): one base below the floor still passes.
            assert!(passes_quality_gate(
                &[10, 20, 20, 20, 20, 20, 20, 20],
                MIN_REGION_Q1,
                &mut buffer
            ));
            assert!(!passes_quality_gate(
                &[10, 10, 20, 20, 20, 20, 20, 20],
                MIN_REGION_Q1,
                &mut buffer
            ));
        }
    }
}

// ---------------------------------------------------------------------
// D2b — classify one read against a locus: delimit (the alignment module's flat-gap
// aligner), recover a long allele by widening, gate the tract quality, and map the result
// to an observation or a no-observation reason. The per-read pipeline ported from
// production's `classify_read`, over the new delimiter. Consumed by the tally + generator.
// ---------------------------------------------------------------------
mod classify {
    #![allow(dead_code)] // consumed by the tally + generator (D2c, D3)

    use super::SsrLocus;
    use super::read_region::{
        MIN_REGION_Q1, extract_region, flank_truncated, passes_quality_gate, read_footprint,
        widen_region,
    };
    use crate::bam::alignment_input::MappedRead;
    use crate::ng::alignment::ssr_best_path_flat_gap::{SsrFlatGapAligner, ViterbiScratch};
    use crate::ng::alignment::{
        BestPathAligner, Emission, ReadBases, RepeatContext, RepeatGeometry, RepeatSpan,
        StutterModel,
    };
    use crate::ng::locus_generation::ReadCoverage;
    use crate::ng::types::Bp;
    use std::ops::Range;

    /// Why a read yielded no usable observation — the tally increments the matching
    /// `SsrGeneratorCounts` reason (and `reads_without_observation`).
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub(super) enum NoObservationReason {
        /// Neither flank anchored (the read lies wholly inside the repeat), or a malformed read.
        NoBorderAnchored,
        /// The delimited tract failed the base-quality gate.
        LowQuality,
        /// The allele stayed flank-truncated even after widening.
        WindowTruncated,
    }

    /// What one read contributes to a locus: an observed tract (complete or partial), or a
    /// reason it observed nothing.
    #[derive(Debug, Clone, PartialEq)]
    pub(super) enum Classified {
        Observed {
            bases: Box<[u8]>,
            coverage: ReadCoverage,
            /// The read's base-quality error mass over the tract, `Σ ln(P_err)`, in log-error
            /// space (freebayes' `q_sum` convention). The **BQ** support moment the spec names
            /// (spec §3, "strand/BQ/MAPQ moments") — computed here because the tract base
            /// qualities are already sliced, a free by-product; the tally folds it into
            /// [`ObservedSequence::q_sum`](crate::ng::locus_generation::ObservedSequence). MAPQ
            /// is carried separately, off the [`MappedRead`], in the tally. **Soft** — filled,
            /// unconsumed today.
            q_sum: f64,
        },
        NoObservation(NoObservationReason),
    }

    /// Classify one read against `locus` using `aligner`. `stutter` feeds the context (the
    /// flat-gap aligner ignores it); `align_scratch` and `qual_buffer` are reused across reads.
    ///
    /// The pipeline mirrors production: extract the read's slice over the locus window, delimit
    /// it, and — if a *complete* tract looks truncated by the window (a mapper-collapsed long
    /// allele) — widen the slice and re-delimit, giving up as `WindowTruncated` if it stays
    /// truncated. A complete tract is quality-gated; **partials are kept as lower bounds**
    /// without the gate (spec §3, the new behaviour production discards).
    pub(super) fn classify_read<E: Emission>(
        read: &MappedRead,
        locus: &SsrLocus,
        aligner: &SsrFlatGapAligner<E>,
        stutter: &StutterModel,
        align_scratch: &mut ViterbiScratch,
        qual_buffer: &mut Vec<u8>,
    ) -> Classified {
        // A malformed record (qual length ≠ seq length) cannot be sliced safely.
        if read.qual.len() != read.seq.len() {
            return Classified::NoObservation(NoObservationReason::NoBorderAnchored);
        }
        let read_len = read.seq.len();
        let reference = &*locus.tract_with_margin_bases;
        let left = locus.left_flank_len();
        let right = locus.right_flank_len();
        let geometry = RepeatGeometry {
            left_flank_len: Bp(left as u64),
            right_flank_len: Bp(right as u64),
            motif: locus.segment.motif(),
        };
        let context = RepeatContext {
            geometry: &geometry,
            stutter,
        };

        let fp = read_footprint(&read.cigar, read.pos);
        let window_start = (locus.margin_start.get() - 1) as u32; // 0-based
        let window_len = reference.len() as u32;
        let region = extract_region(&read.cigar, fp, read_len, window_start, window_len);

        let span = delimit(aligner, read, &region, reference, context, align_scratch);
        match span {
            RepeatSpan::Between(tract)
                if flank_truncated(&region, &to_usize(&tract), read_len, left, right) =>
            {
                // A complete tract whose flanks look eaten by the window: widen and retry.
                let wide = widen_region(region, read_len, left, right);
                let wspan = delimit(aligner, read, &wide, reference, context, align_scratch);
                match wspan {
                    RepeatSpan::Between(wtract)
                        if !flank_truncated(&wide, &to_usize(&wtract), read_len, left, right) =>
                    {
                        complete_or_low_quality(read, &wide, &wtract, qual_buffer)
                    }
                    _ => Classified::NoObservation(NoObservationReason::WindowTruncated),
                }
            }
            RepeatSpan::Between(tract) => {
                complete_or_low_quality(read, &region, &tract, qual_buffer)
            }
            RepeatSpan::FromLeft(tract) => {
                partial(read, &region, &tract, ReadCoverage::PartialLeft)
            }
            RepeatSpan::FromRight(tract) => {
                partial(read, &region, &tract, ReadCoverage::PartialRight)
            }
            RepeatSpan::Unanchored => {
                Classified::NoObservation(NoObservationReason::NoBorderAnchored)
            }
        }
    }

    /// Align the read's `region` slice against the reference frame.
    fn delimit<E: Emission>(
        aligner: &SsrFlatGapAligner<E>,
        read: &MappedRead,
        region: &Range<usize>,
        reference: &[u8],
        context: RepeatContext<'_>,
        scratch: &mut ViterbiScratch,
    ) -> RepeatSpan {
        let bases = ReadBases::try_new(&read.seq[region.clone()], &read.qual[region.clone()])
            .expect("seq and qual are the same slice, hence equal length");
        aligner.align(bases, reference, context, scratch)
    }

    fn to_usize(range: &Range<u64>) -> Range<usize> {
        range.start as usize..range.end as usize
    }

    /// A complete tract, if it clears the base-quality gate; else `LowQuality`. `tract` is
    /// relative to the `region` slice.
    fn complete_or_low_quality(
        read: &MappedRead,
        region: &Range<usize>,
        tract: &Range<u64>,
        qual_buffer: &mut Vec<u8>,
    ) -> Classified {
        let region_seq = &read.seq[region.clone()];
        let region_qual = &read.qual[region.clone()];
        let tract = to_usize(tract);
        let tract_qual = &region_qual[tract.clone()];
        if passes_quality_gate(tract_qual, MIN_REGION_Q1, qual_buffer) {
            Classified::Observed {
                bases: region_seq[tract].into(),
                coverage: ReadCoverage::Complete,
                q_sum: ln_p_err_sum(tract_qual),
            }
        } else {
            Classified::NoObservation(NoObservationReason::LowQuality)
        }
    }

    /// `Σ ln(P_err)` over base qualities — the per-read base-quality error mass on the tract, in
    /// log-error space (freebayes' `q_sum` convention). Negative and monotone: higher quality →
    /// less error mass. The BQ support moment (spec §3); there is no BAQ on this path and MAPQ is
    /// carried separately, so this is the base-quality term alone.
    fn ln_p_err_sum(quals: &[u8]) -> f64 {
        const LN10_OVER_10: f64 = std::f64::consts::LN_10 / 10.0;
        quals.iter().map(|&q| -(q as f64) * LN10_OVER_10).sum()
    }

    /// A partial (lower-bound) observation: the tract bases the read showed, tagged with which
    /// border held and how far it reached (in read coordinates, clamped to `u16`).
    fn partial(
        read: &MappedRead,
        region: &Range<usize>,
        tract: &Range<u64>,
        coverage: fn(u16) -> ReadCoverage,
    ) -> Classified {
        let region_seq = &read.seq[region.clone()];
        let region_qual = &read.qual[region.clone()];
        let tract = to_usize(tract);
        // `reach` is the observed tract length in **read** coordinates. `ReadCoverage`'s reach
        // is nominally in locus positions, which diverge from read bases under stutter; the
        // approximation is bounded because `num_obs_along_locus` clamps the reach to the locus
        // length (a long partial correctly saturates to full coverage).
        let reach = (tract.end - tract.start).min(u16::MAX as usize) as u16;
        Classified::Observed {
            bases: region_seq[tract.clone()].into(),
            coverage: coverage(reach),
            // Partials are kept without the quality gate (spec §3), but still carry the BQ moment.
            q_sum: ln_p_err_sum(&region_qual[tract]),
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::ng::alignment::PerQualityEmission;
        use crate::ng::region_typing::segment_criteria::{Motif, SsrSegment};
        use crate::ng::types::Position;
        use crate::pileup::walker::CigarOp;

        // Reference frame: 6-base flanks around a CACACA tract → "GGGGGGCACACATTTTTT".
        const FRAME: &[u8] = b"GGGGGGCACACATTTTTT";

        fn locus() -> SsrLocus {
            SsrLocus {
                segment: SsrSegment::new("chr1".into(), 7, 12, Motif::new(b"CA").unwrap(), 1.0)
                    .unwrap(),
                tract_with_margin_bases: FRAME.into(),
                margin_start: Position(1),
            }
        }

        fn read(seq: &[u8], qual_value: u8) -> MappedRead {
            MappedRead {
                qname: b"r".to_vec(),
                flag: 0,
                ref_id: 0,
                pos: 1,
                mapq: 60,
                cigar: vec![CigarOp::Match(seq.len() as u32)],
                seq: seq.to_vec(),
                qual: vec![qual_value; seq.len()],
                mate_ref_id: None,
                mate_pos: None,
                adaptor_boundary: None,
                source_file_index: 0,
            }
        }

        fn classify(read: &MappedRead) -> Classified {
            let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
            let stutter = StutterModel::hipstr_shipped();
            let mut scratch = ViterbiScratch::new();
            let mut qual_buffer = Vec::new();
            classify_read(
                read,
                &locus(),
                &aligner,
                &stutter,
                &mut scratch,
                &mut qual_buffer,
            )
        }

        /// A read spanning the whole frame measures the tract exactly — a complete observation
        /// of the reference allele.
        #[test]
        fn a_spanning_read_yields_a_complete_observation() {
            match classify(&read(FRAME, 40)) {
                Classified::Observed {
                    bases,
                    coverage: ReadCoverage::Complete,
                    q_sum,
                } => {
                    assert_eq!(&*bases, b"CACACA");
                    // Six Q40 bases: Σ ln(P_err) = 6 · 40 · (−ln10/10), all negative.
                    assert!(q_sum < 0.0, "the BQ error mass is negative, got {q_sum}");
                }
                other => panic!("expected a complete observation, got {other:?}"),
            }
        }

        /// A read whose tract quality is below the gate is `LowQuality`, not an observation.
        #[test]
        fn a_low_quality_tract_is_rejected() {
            // Whole read at quality 5 (< MIN_REGION_Q1 = 15).
            let classified = classify(&read(FRAME, 5));
            assert_eq!(
                classified,
                Classified::NoObservation(NoObservationReason::LowQuality)
            );
        }

        /// A read that anchors the left flank but ends inside the tract is a partial (left)
        /// observation — a lower bound production would have discarded.
        #[test]
        fn a_read_running_off_the_right_is_a_left_partial() {
            // left flank + 4 tract bases, then the read ends (no right flank).
            let classified = classify(&read(b"GGGGGGCACA", 40));
            match classified {
                Classified::Observed {
                    bases,
                    coverage: ReadCoverage::PartialLeft(reach),
                    ..
                } => {
                    assert_eq!(&*bases, b"CACA");
                    assert_eq!(reach, 4);
                }
                other => panic!("expected a left partial, got {other:?}"),
            }
        }

        /// A read that anchors the right flank but begins inside the tract is a partial
        /// (right) — the mirror of the left case, which a `PartialLeft`/`PartialRight` swap
        /// would fail.
        #[test]
        fn a_read_running_off_the_left_is_a_right_partial() {
            // 4 tract bases + full right flank, mapped at the tract's 3rd base (0-based 8).
            let mut read = read(b"CACATTTTTT", 40);
            read.pos = 9;
            match classify(&read) {
                Classified::Observed {
                    bases,
                    coverage: ReadCoverage::PartialRight(reach),
                    ..
                } => {
                    assert_eq!(&*bases, b"CACA");
                    assert_eq!(reach, 4);
                }
                other => panic!("expected a right partial, got {other:?}"),
            }
        }

        /// A read wholly inside the tract anchors neither flank — no per-read fact, so no
        /// observation.
        #[test]
        fn a_read_wholly_inside_the_tract_is_unanchored() {
            let mut read = read(b"CACACA", 40); // exactly the tract, no flanks
            read.pos = 7;
            assert_eq!(
                classify(&read),
                Classified::NoObservation(NoObservationReason::NoBorderAnchored)
            );
        }

        /// A long-allele read the mapper laid down all-Match pushes its far flank out of the
        /// ref-sized window; the first delimit sees a truncated flank, and widening the read
        /// slice recovers the full long allele (spec §2, the long-allele window recovery).
        #[test]
        fn a_window_truncated_long_allele_is_recovered_by_widening() {
            // 5 CA units (10 bp tract) instead of the reference's 3, mapped 22M.
            match classify(&read(b"GGGGGGCACACACACATTTTTT", 40)) {
                Classified::Observed {
                    bases,
                    coverage: ReadCoverage::Complete,
                    ..
                } => assert_eq!(&*bases, b"CACACACACA"),
                other => panic!("expected the recovered long allele, got {other:?}"),
            }
        }
    }
}

// ---------------------------------------------------------------------
// D2c — the tally: fold each kept read's `Classified` into the locus's observed sequences,
// deduping by `(bases, read_coverage)`, accumulating the support moments, and counting the
// no-observation reasons. A port of production's `tally` (`src/ssr/pileup/locus_tally.rs`),
// extended with partial observations and the strand/BQ/MAPQ moments the shared type carries
// (spec §3, §6). Consumed by the generator (D3).
// ---------------------------------------------------------------------
mod tally {
    #![allow(dead_code)] // consumed by the generator (D3)

    use super::SsrGeneratorCounts;
    use super::classify::{Classified, NoObservationReason};
    use crate::bam::alignment_input::{FLAG_REVERSE_STRAND, MappedRead};
    use crate::ng::locus_generation::{ObservedSequence, ReadCoverage};
    use std::collections::HashMap;

    /// The per-locus tally the generator folds onto the `SampleLocusObservations`: the deduped
    /// observed sequences and how many reads reached the aligner but yielded nothing.
    /// (`reads_fetched` / `reads_discarded_by_cap` are the caller's to set — they come from the
    /// cap, not the outcomes.)
    pub(super) struct SsrTally {
        pub(super) observed_sequences: Vec<ObservedSequence>,
        pub(super) reads_without_observation: u32,
    }

    /// Accumulated support for one distinct `(bases, read_coverage)` bucket — the moments summed
    /// as reads fold in, materialised into an [`ObservedSequence`] at the end.
    #[derive(Default)]
    struct Support {
        num_obs: u32,
        num_fwd: u32,
        q_sum: f64,
        mapq_sum: u32,
        mapq_sum_sq: u64,
    }

    /// Fold each kept read's classification into the locus tally.
    ///
    /// Observations dedup by **`(bases, read_coverage)`** — a `Complete` and a partial of the
    /// same bases are different evidence and stay separate rows, two identical partials merge
    /// (spec §3). Each bucket accumulates the strand (`num_fwd` off the reverse-strand flag),
    /// BQ (`q_sum`, off the read's tract error mass) and MAPQ moments; `chain_ids` stays empty
    /// because the STR path does not phase. `observed_sequences` is sorted by `(bases,
    /// coverage)`, so — like production's `tally` — the bases, the counts and the integer moments
    /// are independent of the order reads were folded. `q_sum` is the one exception: it sums in
    /// `f64`, which is commutative but not associative, so a bucket's `q_sum` can differ in its
    /// low bits under a different fold order. That is immaterial while `q_sum` is soft and
    /// unconsumed, and the parity oracle checks only bytes and counts (spec §6).
    ///
    /// `counts` carries the **run-level** totals (complete/partial observations and the three
    /// no-observation reasons); the returned `reads_without_observation` is this locus's own
    /// total.
    pub(super) fn tally<'a>(
        reads_and_outcomes: impl IntoIterator<Item = (&'a MappedRead, Classified)>,
        counts: &mut SsrGeneratorCounts,
    ) -> SsrTally {
        let mut buckets: HashMap<(Box<[u8]>, ReadCoverage), Support> = HashMap::new();
        let mut reads_without_observation = 0u32;
        for (read, outcome) in reads_and_outcomes {
            match outcome {
                Classified::Observed {
                    bases,
                    coverage,
                    q_sum,
                } => {
                    match coverage {
                        ReadCoverage::Complete => counts.observations_complete += 1,
                        ReadCoverage::PartialLeft(_) | ReadCoverage::PartialRight(_) => {
                            counts.observations_partial += 1
                        }
                    }
                    let support = buckets.entry((bases, coverage)).or_default();
                    support.num_obs += 1;
                    if read.flag & FLAG_REVERSE_STRAND == 0 {
                        support.num_fwd += 1;
                    }
                    support.q_sum += q_sum;
                    let mapq = u32::from(read.mapq);
                    support.mapq_sum += mapq;
                    support.mapq_sum_sq += u64::from(mapq) * u64::from(mapq);
                }
                Classified::NoObservation(reason) => {
                    reads_without_observation += 1;
                    match reason {
                        NoObservationReason::NoBorderAnchored => counts.no_border_anchored += 1,
                        NoObservationReason::LowQuality => counts.low_quality += 1,
                        NoObservationReason::WindowTruncated => counts.window_truncated += 1,
                    }
                }
            }
        }

        let mut observed_sequences: Vec<ObservedSequence> = buckets
            .into_iter()
            .map(|((bases, read_coverage), support)| ObservedSequence {
                bases,
                read_coverage,
                num_obs: support.num_obs,
                num_fwd: support.num_fwd,
                q_sum: support.q_sum,
                mapq_sum: support.mapq_sum,
                mapq_sum_sq: support.mapq_sum_sq,
                // The STR path does not phase — there is no chain id to fold (spec §3).
                chain_ids: Vec::new(),
            })
            .collect();
        observed_sequences.sort_unstable_by(|a, b| {
            a.bases
                .cmp(&b.bases)
                .then_with(|| coverage_order(a.read_coverage).cmp(&coverage_order(b.read_coverage)))
        });

        SsrTally {
            observed_sequences,
            reads_without_observation,
        }
    }

    /// A total order over `ReadCoverage` (which is not `Ord`) for a deterministic tie-break when
    /// two observations share bases: complete first, then left partials, then right partials,
    /// each ordered by reach.
    fn coverage_order(coverage: ReadCoverage) -> (u8, u16) {
        match coverage {
            ReadCoverage::Complete => (0, 0),
            ReadCoverage::PartialLeft(n) => (1, n),
            ReadCoverage::PartialRight(n) => (2, n),
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::pileup::walker::CigarOp;

        /// A `MappedRead` with a given strand flag and MAPQ — the only fields the tally reads
        /// off the read; the sequence and qualities live in the `Classified` handed alongside.
        fn read(flag: u16, mapq: u8) -> MappedRead {
            MappedRead {
                qname: b"r".to_vec(),
                flag,
                ref_id: 0,
                pos: 1,
                mapq,
                cigar: vec![CigarOp::Match(6)],
                seq: b"CACACA".to_vec(),
                qual: vec![40; 6],
                mate_ref_id: None,
                mate_pos: None,
                adaptor_boundary: None,
                source_file_index: 0,
            }
        }

        fn observed(bases: &[u8], coverage: ReadCoverage, q_sum: f64) -> Classified {
            Classified::Observed {
                bases: bases.into(),
                coverage,
                q_sum,
            }
        }

        /// A `Complete` and a `PartialLeft` of the **same** bases are different evidence, so they
        /// stay as two separate rows (spec §3) — the property the `(bases, read_coverage)` dedup
        /// key rests on.
        #[test]
        fn a_complete_and_a_partial_of_the_same_bases_stay_separate() {
            let fwd = read(0, 60);
            let outcomes = vec![
                (&fwd, observed(b"CACACA", ReadCoverage::Complete, -1.0)),
                (
                    &fwd,
                    observed(b"CACACA", ReadCoverage::PartialLeft(6), -1.0),
                ),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, &mut counts);
            assert_eq!(result.observed_sequences.len(), 2);
            assert_eq!(counts.observations_complete, 1);
            assert_eq!(counts.observations_partial, 1);
        }

        /// Two identical partials (same bases, same coverage) are the identical constraint, so
        /// they merge into one row with `num_obs == 2` (spec §3).
        #[test]
        fn two_identical_partials_merge_into_one_count() {
            let fwd = read(0, 60);
            let outcomes = vec![
                (&fwd, observed(b"CACA", ReadCoverage::PartialLeft(4), -1.0)),
                (&fwd, observed(b"CACA", ReadCoverage::PartialLeft(4), -1.0)),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, &mut counts);
            assert_eq!(result.observed_sequences.len(), 1);
            assert_eq!(result.observed_sequences[0].num_obs, 2);
        }

        /// The strand, BQ and MAPQ moments accumulate across the reads folded into a bucket:
        /// `num_fwd` counts the forward-strand reads, `q_sum` sums the per-read BQ mass,
        /// `mapq_sum` / `mapq_sum_sq` sum MAPQ and MAPQ². `chain_ids` stays empty.
        #[test]
        fn the_support_moments_accumulate() {
            let fwd = read(0, 60);
            let rev = read(FLAG_REVERSE_STRAND, 30);
            let outcomes = vec![
                (&fwd, observed(b"CACACA", ReadCoverage::Complete, -2.0)),
                (&rev, observed(b"CACACA", ReadCoverage::Complete, -3.0)),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, &mut counts);
            let obs = &result.observed_sequences[0];
            assert_eq!(obs.num_obs, 2);
            assert_eq!(obs.num_fwd, 1, "one forward, one reverse");
            assert_eq!(obs.q_sum, -5.0, "the BQ masses sum");
            assert_eq!(obs.mapq_sum, 90, "60 + 30");
            assert_eq!(obs.mapq_sum_sq, 60 * 60 + 30 * 30);
            assert!(obs.chain_ids.is_empty(), "the STR path does not phase");
        }

        /// Each no-observation reason lands in its own run-level counter, and every such read is
        /// counted in the locus's `reads_without_observation` total.
        #[test]
        fn no_observation_reasons_are_counted_by_reason_and_in_total() {
            let r = read(0, 60);
            let outcomes = vec![
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::NoBorderAnchored),
                ),
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::LowQuality),
                ),
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::WindowTruncated),
                ),
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::LowQuality),
                ),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, &mut counts);
            assert!(result.observed_sequences.is_empty());
            assert_eq!(result.reads_without_observation, 4);
            assert_eq!(counts.no_border_anchored, 1);
            assert_eq!(counts.low_quality, 2);
            assert_eq!(counts.window_truncated, 1);
        }

        /// `observed_sequences` is sorted by bytes, then by coverage — so the record is identical
        /// regardless of the order reads folded in (production's order-independence, extended to
        /// the coverage tie-break).
        #[test]
        fn observed_is_sorted_by_bases_then_coverage() {
            let r = read(0, 60);
            let outcomes = vec![
                (&r, observed(b"GG", ReadCoverage::Complete, -1.0)),
                (&r, observed(b"AA", ReadCoverage::PartialRight(2), -1.0)),
                (&r, observed(b"AA", ReadCoverage::Complete, -1.0)),
                (&r, observed(b"AA", ReadCoverage::PartialLeft(2), -1.0)),
            ];
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, &mut counts);
            let order: Vec<(&[u8], ReadCoverage)> = result
                .observed_sequences
                .iter()
                .map(|o| (o.bases.as_ref(), o.read_coverage))
                .collect();
            assert_eq!(
                order,
                vec![
                    (b"AA".as_ref(), ReadCoverage::Complete),
                    (b"AA".as_ref(), ReadCoverage::PartialLeft(2)),
                    (b"AA".as_ref(), ReadCoverage::PartialRight(2)),
                    (b"GG".as_ref(), ReadCoverage::Complete),
                ]
            );
        }

        /// The integer moments of a **multi-read bucket** are order-independent: the same reads
        /// (differing strand and MAPQ) folded into one `(bases, coverage)` bucket in two orders
        /// give the same `num_obs` / `num_fwd` / `mapq_sum` / `mapq_sum_sq`. This is the case the
        /// singleton buckets of `tally_is_order_independent` never exercise, and it isolates
        /// `q_sum` as the sole order-sensitive field (its f64 sum is not asserted here).
        #[test]
        fn a_multi_read_bucket_accumulates_integer_moments_order_independently() {
            let fwd = read(0, 60);
            let rev = read(FLAG_REVERSE_STRAND, 30);
            let moments = |outcomes: Vec<(&MappedRead, Classified)>| {
                let mut counts = SsrGeneratorCounts::default();
                let obs = tally(outcomes, &mut counts).observed_sequences;
                assert_eq!(obs.len(), 1, "all reads share one bucket");
                let o = &obs[0];
                (o.num_obs, o.num_fwd, o.mapq_sum, o.mapq_sum_sq)
            };
            let forward_first = moments(vec![
                (&fwd, observed(b"CACACA", ReadCoverage::Complete, -2.0)),
                (&rev, observed(b"CACACA", ReadCoverage::Complete, -3.0)),
            ]);
            let reverse_first = moments(vec![
                (&rev, observed(b"CACACA", ReadCoverage::Complete, -3.0)),
                (&fwd, observed(b"CACACA", ReadCoverage::Complete, -2.0)),
            ]);
            assert_eq!(forward_first, reverse_first);
            assert_eq!(forward_first, (2, 1, 90, 60 * 60 + 30 * 30));
        }

        /// The tally is order-independent: the same multiset of outcomes folded in two different
        /// orders yields the same observed sequences and the same counts. (Its buckets are
        /// singletons, so no `q_sum` f64 reordering occurs — the multi-read integer case is
        /// `a_multi_read_bucket_accumulates_integer_moments_order_independently`.)
        #[test]
        fn tally_is_order_independent() {
            let r = read(0, 60);
            let run = |outcomes: Vec<(&MappedRead, Classified)>| {
                let mut counts = SsrGeneratorCounts::default();
                let result = tally(outcomes, &mut counts);
                (result.observed_sequences, counts)
            };
            let a = run(vec![
                (&r, observed(b"CACACA", ReadCoverage::Complete, -1.0)),
                (&r, observed(b"CACA", ReadCoverage::Complete, -1.0)),
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::LowQuality),
                ),
            ]);
            let b = run(vec![
                (
                    &r,
                    Classified::NoObservation(NoObservationReason::LowQuality),
                ),
                (&r, observed(b"CACA", ReadCoverage::Complete, -1.0)),
                (&r, observed(b"CACACA", ReadCoverage::Complete, -1.0)),
            ]);
            assert_eq!(a, b);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::ng::region_typing::segment_criteria::Motif;
    use std::collections::HashSet;

    /// A 100-base contig `chr1` with a known repeating pattern, so a fetched span can be
    /// checked byte-for-byte.
    fn reference_100() -> (InMemoryRefSeq, Vec<u8>) {
        let bases: Vec<u8> = b"ACGT".iter().cycle().take(100).copied().collect();
        (
            InMemoryRefSeq::from_named_contigs(vec![("chr1".to_string(), bases.clone())]),
            bases,
        )
    }

    fn tract(start: u64, end: u64) -> SsrSegment {
        SsrSegment::new("chr1".into(), start, end, Motif::new(b"AC").unwrap(), 1.0).unwrap()
    }

    /// A tract with room either side: both flanks are exactly `flank_bp`, the margin span is
    /// `tract + 2·flank_bp`, and the fetched bytes are exactly that slice of the reference.
    #[test]
    fn a_mid_contig_tract_fetches_equal_flanks_and_the_exact_bytes() {
        let (reference, bases) = reference_100();
        let mut buffer = Vec::new();
        // tract [40, 49] (10 bases), flank 10 → margin [30, 59].
        let locus =
            SsrLocus::fetch(&reference, ContigId(0), tract(40, 49), Bp(10), &mut buffer).unwrap();

        assert_eq!(locus.margin_start, Position(30));
        assert_eq!(locus.left_flank_len(), 10);
        assert_eq!(locus.right_flank_len(), 10);
        assert_eq!(locus.tract_with_margin_bases.len(), 30); // 10 + 10 + 10
        // margin [30, 59] 1-based == bases[29..59] 0-based.
        assert_eq!(&*locus.tract_with_margin_bases, &bases[29..59]);
    }

    /// A tract near the contig **start**: the left margin clamps to base 1, so the left flank
    /// is short while the right is full — the flanks are unequal, measured not assumed.
    #[test]
    fn a_tract_near_the_contig_start_has_a_short_left_flank() {
        let (reference, bases) = reference_100();
        let mut buffer = Vec::new();
        // tract [5, 14], flank 10 → margin_start clamps to 1, margin_end 24.
        let locus =
            SsrLocus::fetch(&reference, ContigId(0), tract(5, 14), Bp(10), &mut buffer).unwrap();

        assert_eq!(locus.margin_start, Position(1));
        assert_eq!(locus.left_flank_len(), 4, "clamped: 5 - 1");
        assert_eq!(locus.right_flank_len(), 10, "full");
        assert_ne!(locus.left_flank_len(), locus.right_flank_len());
        assert_eq!(&*locus.tract_with_margin_bases, &bases[0..24]);
    }

    /// A tract near the contig **end**: the right margin clamps to the contig length, so the
    /// right flank is short — again unequal.
    #[test]
    fn a_tract_near_the_contig_end_has_a_short_right_flank() {
        let (reference, bases) = reference_100();
        let mut buffer = Vec::new();
        // contig length 100; tract [92, 96], flank 10 → margin [82, 100].
        let locus =
            SsrLocus::fetch(&reference, ContigId(0), tract(92, 96), Bp(10), &mut buffer).unwrap();

        assert_eq!(locus.margin_start, Position(82));
        assert_eq!(locus.left_flank_len(), 10, "full");
        assert_eq!(locus.right_flank_len(), 4, "clamped: 100 - 96");
        assert_eq!(locus.tract_with_margin_bases.len(), 19); // 10 + 5 + 4
        assert_eq!(&*locus.tract_with_margin_bases, &bases[81..100]);
    }

    /// A tract reaching past the contig end is a broken segment: `fetch` leaves the window
    /// unclamped so `fetch_into` rejects it, rather than returning a locus whose derived
    /// right flank would underflow.
    #[test]
    fn a_tract_past_the_contig_end_is_rejected_not_silently_truncated() {
        let (reference, _bases) = reference_100();
        let mut buffer = Vec::new();
        // contig length 100; tract [98, 110] runs past the end.
        assert!(
            SsrLocus::fetch(&reference, ContigId(0), tract(98, 110), Bp(10), &mut buffer).is_err()
        );
    }

    #[test]
    fn default_flank_equals_the_bundle_threshold_default_and_caps_at_1000() {
        let config = SsrGeneratorConfig::default();
        assert_eq!(
            config.flank_bp,
            Bp(crate::ng::region_typing::segment_criteria::DEFAULT_BUNDLE_THRESHOLD),
            "flank equals the bundle threshold by default (spec §4)"
        );
        assert_eq!(
            config.max_reads_per_locus,
            Some(DEFAULT_SSR_MAX_READS_PER_LOCUS)
        );
    }

    #[test]
    fn the_flank_check_accepts_within_and_at_the_threshold() {
        // Equal is allowed (the default case), and strictly narrower.
        assert!(
            SsrGeneratorConfig {
                flank_bp: Bp(30),
                max_reads_per_locus: None,
            }
            .check_flank_within(Bp(30))
            .is_ok()
        );
        assert!(
            SsrGeneratorConfig {
                flank_bp: Bp(20),
                max_reads_per_locus: None,
            }
            .check_flank_within(Bp(30))
            .is_ok()
        );
    }

    #[test]
    fn the_flank_check_rejects_a_flank_wider_than_the_bundle_threshold() {
        let error = SsrGeneratorConfig {
            flank_bp: Bp(50),
            max_reads_per_locus: None,
        }
        .check_flank_within(Bp(30))
        .expect_err("50 > 30 must be refused");
        assert!(matches!(
            error,
            SsrGeneratorConfigError::FlankExceedsBundleThreshold {
                flank_bp: 50,
                bundle_threshold: 30,
            }
        ));
    }

    // --- the reservoir port, oracle'd against production's own tests ---------

    #[test]
    fn keeps_everything_when_offered_at_most_capacity() {
        let mut r = Reservoir::new(5, locus_seed("chr1", 100));
        for x in [10u32, 20, 30] {
            r.offer(x);
        }
        assert_eq!(r.seen(), 3);
        assert_eq!(r.into_held(), vec![10, 20, 30]); // first-K kept, in order
    }

    #[test]
    fn caps_at_capacity_and_counts_all_offers() {
        let mut r = Reservoir::new(10, locus_seed("chr1", 100));
        for x in 1..=100u32 {
            r.offer(x);
        }
        assert_eq!(r.seen(), 100);
        assert_eq!(r.into_held().len(), 10);
    }

    #[test]
    fn reservoir_is_deterministic_for_a_fixed_seed_and_order() {
        let run = || {
            let mut r = Reservoir::new(10, locus_seed("chr7", 4242));
            for x in 1..=100u32 {
                r.offer(x);
            }
            r.into_held()
        };
        assert_eq!(run(), run());
    }

    #[test]
    fn reservoir_keeps_a_deterministic_subset_far_past_capacity() {
        let run = || {
            let mut r = Reservoir::new(8, locus_seed("chrX", 7));
            for x in 1..=10_000u32 {
                r.offer(x);
            }
            (r.seen(), r.into_held())
        };
        let (seen, held) = run();
        assert_eq!(seen, 10_000);
        assert_eq!(held.len(), 8);
        assert!(
            held.iter().all(|x| (1..=10_000).contains(x)),
            "kept set is a subset of the stream"
        );
        assert_eq!(run().1, held, "kept set is identical across runs");
    }

    #[test]
    fn different_loci_sample_differently() {
        let sample = |chrom, start| {
            let mut r = Reservoir::new(10, locus_seed(chrom, start));
            for x in 1..=100u32 {
                r.offer(x);
            }
            r.into_held()
        };
        assert_ne!(sample("chr1", 100), sample("chr1", 101));
        assert_ne!(sample("chr1", 100), sample("chr2", 100));
    }

    #[test]
    fn every_item_can_be_selected_no_structural_exclusion() {
        let mut covered = HashSet::new();
        for seed in 0..500u64 {
            let mut r = Reservoir::new(10, seed);
            for x in 1..=100u32 {
                r.offer(x);
            }
            covered.extend(r.into_held());
        }
        assert_eq!(covered.len(), 100);
    }

    #[test]
    fn locus_seed_is_deterministic_and_distinguishes_loci() {
        assert_eq!(locus_seed("chr1", 100), locus_seed("chr1", 100));
        assert_ne!(locus_seed("chr1", 100), locus_seed("chr1", 101));
        assert_ne!(locus_seed("chr1", 100), locus_seed("chr2", 100));
    }

    /// The seed trap: [`seed_for_segment`] seeds from the contig **name** and the **0-based**
    /// start (`start - 1`), not the 1-based start — feeding the 1-based start would produce a
    /// different, deterministic kept set that fails parity looking like an aligner bug (spec §4).
    #[test]
    fn seed_for_segment_uses_the_contig_name_and_the_zero_based_start() {
        let segment = tract(101, 110); // 1-based start 101 → 0-based 100
        assert_eq!(seed_for_segment(&segment), locus_seed("chr1", 100));
        assert_ne!(
            seed_for_segment(&segment),
            locus_seed("chr1", 101),
            "the 1-based start is the trap the conversion avoids"
        );
    }

    /// **The parity oracle for the port**: ng's seed and reservoir must produce output
    /// identical to frozen production (`src/ssr/pileup/fetch_reads.rs`), byte for byte. The
    /// self-consistency tests above would survive a drifted constant; this one would not —
    /// it is what makes "byte-faithful port" a checked claim rather than an asserted one.
    /// (Calling production as a test-only oracle mirrors region typing's `build_loci`
    /// differential; ng does not depend on production at run time.)
    #[test]
    fn ng_seed_and_reservoir_match_frozen_production_byte_for_byte() {
        use crate::ssr::pileup::fetch_reads as production;

        for (chrom, start) in [
            ("chr1", 0u32),
            ("chr1", 100),
            ("chrX", 4242),
            ("scaffold_7", 999_999),
        ] {
            assert_eq!(
                locus_seed(chrom, start),
                production::locus_seed(chrom, start),
                "seed for ({chrom}, {start})"
            );
        }

        // The eviction branch, far past capacity — the kept set is where a drifted PRNG
        // constant would show.
        let kept = |seed| {
            let mut r = Reservoir::new(8, seed);
            for x in 1..=10_000u32 {
                r.offer(x);
            }
            r.into_held()
        };
        let prod_kept = |seed| {
            let mut r = production::Reservoir::new(8, seed);
            for x in 1..=10_000u32 {
                r.offer(x);
            }
            r.into_held()
        };
        let seed = locus_seed("chrX", 7);
        assert_eq!(
            kept(seed),
            prod_kept(seed),
            "the kept set must be identical"
        );
    }

    // --- D1: the per-locus read fetch + cap ---------------------------------

    fn span(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// A `RawRefSeq` over the fixture contigs (all `A`s), for `reads_in_region`'s
    /// mismatch-fraction filter — matches the fixture reference the reads are opened against.
    fn fixture_ref_bases() -> InMemoryRefSeq {
        use crate::ng::read::input::test_fixtures::matching_contigs;
        InMemoryRefSeq::from_contigs(
            matching_contigs()
                .iter()
                .map(|(_, len, _)| vec![b'A'; *len])
                .collect(),
        )
    }

    /// Open a `SampleReads` over one indexed BAM of `records`, against the fixture reference.
    fn sample_reads_with(
        records: &[noodles_sam::alignment::RecordBuf],
    ) -> (tempfile::TempDir, tempfile::TempDir, SampleReads) {
        use crate::ng::read::filtering::ReadFilterConfig;
        use crate::ng::read::input::test_fixtures::{
            fixture_reference, header, indexed_bam, matching_contigs,
        };
        let (reference_dir, reference) = fixture_reference(false);
        let (bam_dir, bam) = indexed_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some("NA12878"))],
            ),
            records,
        );
        let reads = SampleReads::open(&[bam], &reference, ReadFilterConfig::default(), false)
            .expect("the fixture sample opens");
        (reference_dir, bam_dir, reads)
    }

    /// The cap keeps `max_reads_per_locus` of the fetched reads and counts them all; `None`
    /// keeps everything (spec §4).
    #[test]
    fn fetch_caps_reads_and_counts_the_fetched() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let records: Vec<_> = (0..5)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 1 + i * 10, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);

        let capped = fetch_capped_reads(&reads, span(1, 100), 42, Some(3), fixture_ref_bases)
            .expect("fetch succeeds");
        assert_eq!(capped.fetched, 5);
        assert_eq!(capped.kept.len(), 3, "capped to 3 of the 5 fetched");

        let uncapped = fetch_capped_reads(&reads, span(1, 100), 42, None, fixture_ref_bases)
            .expect("fetch succeeds");
        assert_eq!(uncapped.fetched, 5);
        assert_eq!(uncapped.kept.len(), 5, "no cap keeps every read");
    }

    /// The kept subset is a deterministic function of the seed (the thread-invariance the seed
    /// buys, spec §4).
    #[test]
    fn the_capped_kept_set_is_deterministic_for_a_seed() {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        let records: Vec<_> = (0..20)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 1 + i * 4, 20))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let kept_positions = |seed| {
            let mut positions: Vec<u64> =
                fetch_capped_reads(&reads, span(1, 100), seed, Some(5), fixture_ref_bases)
                    .expect("fetch succeeds")
                    .kept
                    .iter()
                    .map(|read| read.pos)
                    .collect();
            positions.sort_unstable();
            positions
        };
        assert_eq!(
            kept_positions(7),
            kept_positions(7),
            "same seed → same kept set"
        );
    }
}
