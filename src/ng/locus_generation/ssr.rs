//! The STR locus generator — one microsatellite tract → one locus.
//!
//! The first [`LocusGenerator`](super::LocusGenerator): it consumes an `SsrSegment`,
//! fetches the reads over the tract, aligns each to read off its repeat, and tallies the
//! answers into a [`SampleLocusObservations`](super::SampleLocusObservations). A port of
//! production `src/ssr/pileup/`, adapted at two seams (the split of coordinates from bases,
//! and a widened admission gate). See `doc/devel/ng/spec/locus_generation_ssr.md` (design)
//! and `doc/devel/ng/arch/locus_generation_ssr.md` (types & interfaces).
//!
//! This module lands across the STR generator plan; so far the config, counts, cap constant,
//! and the working input ([`SsrLocus`] + its margin fetch). The align → tally rest of the
//! transform is gated on the ng STR aligner (`read_preparation_ssr.md`), unbuilt.

use crate::ng::ref_seq::{ContigTable, RefSeq, RefSeqError};
use crate::ng::region_typing::segment_criteria::SsrSegment;
use crate::ng::types::{Bp, ContigId, Position};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::ng::region_typing::segment_criteria::Motif;

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
}
