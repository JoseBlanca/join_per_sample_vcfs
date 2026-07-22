//! The STR locus generator — one microsatellite tract → one locus.
//!
//! The first [`LocusGenerator`](super::LocusGenerator): it consumes an `SsrSegment`,
//! fetches the reads over the tract, aligns each to read off its repeat, and tallies the
//! answers into a [`SampleLocusObservations`](super::SampleLocusObservations). A port of
//! production `src/ssr/pileup/`, adapted at two seams (the split of coordinates from bases,
//! and a widened admission gate). See `doc/devel/ng/spec/locus_generation_ssr.md` (design)
//! and `doc/devel/ng/arch/locus_generation_ssr.md` (types & interfaces).
//!
//! This module lands across the STR generator plan; so far its config, counts, and cap
//! constant. The transform (fetch → align → tally) is gated on the ng STR aligner
//! (`read_preparation_ssr.md`), unbuilt.

use crate::ng::types::Bp;

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
