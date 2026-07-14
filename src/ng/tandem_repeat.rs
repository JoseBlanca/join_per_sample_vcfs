//! A clean-room, use-agnostic short-period **tandem-repeat scanner** — a sequence
//! primitive shared by the STR catalog and (later) the ng snp/str caller. Design:
//! spec [`doc/devel/ng/spec/ssr_repeat_scanner.md`], arch
//! [`doc/devel/ng/arch/ssr_repeat_scanner.md`].
//!
//! It exposes **two interfaces** over one detection core (a lag-`p` self-comparison
//! plus a Ruzzo–Tompa maximal-scoring-segment pass):
//!
//! - a low-level **interval finder** — `find_tandem_repeats(seq, PeriodRange,
//!   &ScanParams) -> Vec<RepeatInterval>` — raw, possibly-overlapping intervals, for
//!   consumers that resolve overlaps themselves (the STR catalog); and
//! - a high-level **region seam** — `RegionScanner`, a streaming iterator that yields
//!   a gap-free repeat/satellite/unique tiling of the reference, for the caller's
//!   router.
//!
//! The scanner holds **no consumer policy** (no purity floor, homopolymer rule, or
//! period ceiling): a consumer passes the period range it wants, and that plus the two
//! scoring weights is the whole configuration surface.
//!
//! **Clean-room note:** derived only from Benson (1999)'s published idea (compare to a
//! period-shifted copy; grow runs of matches) and our own code — **not** from the
//! AGPL-v3 `TRF-mod` source, which was not read.
//!
//! Build status (incremental, per the impl plan): **Milestone A — types + scaffold.**
//! The interval finder (Milestone B) and the `RegionScanner` region seam (Milestone C)
//! land as later increments; this file currently carries only the type vocabulary.
//!
//! Visibility note: types are `pub` (matching the sibling ng modules, e.g. `types.rs`),
//! not the arch doc's illustrative `pub(crate)` — the ng convention exposes its module
//! vocabulary as reachable API, which also lets Milestone-A types be defined ahead of
//! their Milestone-B/C consumers without tripping `dead_code`.

use crate::fasta::ChromRefFetchError;

// ---------------------------------------------------------------------
// Inputs — the scope and scoring knobs
// ---------------------------------------------------------------------

/// The inclusive range of repeat **periods** (motif lengths, in bp) to scan for.
///
/// The caller chooses the range; the algorithm imposes no ceiling of its own. This is
/// the one **constrained** type in the module — a period of 0 is meaningless and an
/// empty range (`min > max`) is a caller bug — so its fields are private and every
/// value is built through the checked [`PeriodRange::new`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PeriodRange {
    min: u8,
    max: u8,
}

impl PeriodRange {
    /// Build a period range, rejecting `min == 0` and `min > max`. These are caller
    /// bugs (period 0 has no meaning; an empty range scans nothing), so they fail
    /// loudly rather than silently coercing.
    pub fn new(min: u8, max: u8) -> Result<Self, PeriodRangeError> {
        if min == 0 {
            return Err(PeriodRangeError::ZeroMin);
        }
        if min > max {
            return Err(PeriodRangeError::MinExceedsMax { min, max });
        }
        Ok(Self { min, max })
    }

    /// The smallest period scanned (inclusive).
    #[inline]
    pub fn min(self) -> u8 {
        self.min
    }

    /// The largest period scanned (inclusive).
    #[inline]
    pub fn max(self) -> u8 {
        self.max
    }
}

/// A [`PeriodRange`] could not be built — a caller bug, never a data condition.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PeriodRangeError {
    /// `min` was 0; period 0 is meaningless.
    #[error("period range min must be >= 1 (period 0 is meaningless)")]
    ZeroMin,
    /// `min` exceeded `max`; the range would scan nothing.
    #[error("period range min ({min}) exceeds max ({max})")]
    MinExceedsMax { min: u8, max: u8 },
}

/// The scanner's general default period range, 1..=6 — a named constant, not a magic
/// literal. Each consumer overrides as it sees fit (the catalog / ng snp-str caller
/// uses 2..=6); the value may change after empirical tuning.
pub const DEFAULT_PERIODS: (u8, u8) = (1, 6);

/// Default match reward: a lag-`p` match adds this to the running score.
pub const DEFAULT_MATCH_REWARD: i32 = 2;
/// Default mismatch penalty: a lag-`p` mismatch (or a non-ACGT base) subtracts this.
/// The ratio `mismatch_penalty / match_reward` is the mismatch tolerance — a ratio `r`
/// holds a tract to purity `r / (1 + r)`; 7/2 = 3.5 holds tracts to ≈ 0.78 (spec §3.3).
/// A starting point, tuned from experiments.
pub const DEFAULT_MISMATCH_PENALTY: i32 = 7;
/// Default minimum copy count for a segment to be emitted — a repeat needs at least
/// two copies to exist at all.
pub const DEFAULT_MIN_COPIES: u32 = 2;

/// The lag-`p` self-comparison scoring plus the minimum-copies emission floor.
///
/// A tandem repeat of period `p` shows up as a high-scoring run of the signal that adds
/// [`match_reward`](Self::match_reward) where `seq[j] == seq[j-p]` (both ACGT,
/// case-insensitive) and subtracts [`mismatch_penalty`](Self::mismatch_penalty)
/// otherwise. `Default` is the scanner's general starting weights (`2 / 7 / 2`), which a
/// consumer may override with weights matched to its own purity target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanParams {
    /// Score added on a lag-`p` match. Expected `> 0`.
    pub match_reward: i32,
    /// Score subtracted on a lag-`p` mismatch (or non-ACGT base). Expected `> 0`.
    pub mismatch_penalty: i32,
    /// A segment shorter than `min_copies · period` is not emitted (a permissive floor;
    /// a consumer applies its own, stricter, per-period floors afterward).
    pub min_copies: u32,
}

impl Default for ScanParams {
    fn default() -> Self {
        Self {
            match_reward: DEFAULT_MATCH_REWARD,
            mismatch_penalty: DEFAULT_MISMATCH_PENALTY,
            min_copies: DEFAULT_MIN_COPIES,
        }
    }
}

/// Above this merged-repeat length a region is [`Region::Satellite`], not a genotypeable
/// STR (satellite DNA; spec §3.6). Also bounds the streaming window halo. 1 kb is
/// comfortably above any real STR locus.
pub const DEFAULT_MAX_REPEAT_LEN: u32 = 1000;
/// Default streaming window core (bp): a memory/I/O-granularity knob only — the region
/// tiling is invariant to it.
pub const DEFAULT_WINDOW_BP: u32 = 100_000;
/// Default unique-gap bridging (0 = off): smoothing for the region tiling.
pub const DEFAULT_MERGE_GAP: u32 = 0;
/// Default minimum repeat-region length (0 = off): smoothing for the region tiling.
pub const DEFAULT_MIN_REPEAT_LEN: u32 = 0;

/// Region-seam shaping and streaming knobs — separate from [`ScanParams`] because they
/// shape the **partition and the walk**, not the detection.
///
/// `Default` = `{ 1000, 100_000, 0, 0 }`: the 1 kb satellite cap, a 100 kb streaming
/// window, and both smoothing knobs off (a pure coverage partition).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SegmentOptions {
    /// Merged repeat coverage longer than this is emitted as [`Region::Satellite`].
    pub max_repeat_len: u32,
    /// The streaming window core in bp; a memory knob, region-count-invariant.
    pub window_bp: u32,
    /// Smoothing: bridge unique gaps shorter than this into the flanking repeat.
    pub merge_gap: u32,
    /// Smoothing: reclassify a repeat region shorter than this as unique.
    pub min_repeat_len: u32,
}

impl Default for SegmentOptions {
    fn default() -> Self {
        Self {
            max_repeat_len: DEFAULT_MAX_REPEAT_LEN,
            window_bp: DEFAULT_WINDOW_BP,
            merge_gap: DEFAULT_MERGE_GAP,
            min_repeat_len: DEFAULT_MIN_REPEAT_LEN,
        }
    }
}

// ---------------------------------------------------------------------
// Output — the interval finder's natural shape
// ---------------------------------------------------------------------

/// One detected tandem-repeat interval.
///
/// Coordinates are **0-based half-open** (`[start, end)`); `period` is the lag it was
/// found at; `score` is the Ruzzo–Tompa segment total (an alignment-style
/// length-and-purity proxy). Nothing consumer-specific here — it happens to carry
/// exactly the four fields the STR post-filter reads, which is convenience, not coupling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RepeatInterval {
    /// Tract start, inclusive (0-based).
    pub start: u32,
    /// Tract end, exclusive (0-based).
    pub end: u32,
    /// Repeat period (motif length, in bp).
    pub period: u8,
    /// Ruzzo–Tompa segment score.
    pub score: i32,
}

// ---------------------------------------------------------------------
// Output — the region seam's tiling
// ---------------------------------------------------------------------

/// A half-open, 0-based coordinate pair (`[start, end)`) — the plain span the
/// [`Region::Satellite`] / [`Region::Unique`] variants carry and the union
/// [`RepeatRegion::span`] holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RegionSpan {
    /// Start, inclusive (0-based).
    pub start: u32,
    /// End, exclusive (0-based).
    pub end: u32,
}

/// A genotypeable tandem repeat: merged repeat coverage no longer than
/// [`SegmentOptions::max_repeat_len`]. It carries the merged span **and** the
/// constituent intervals (periods may differ), so the caller can read the repeat
/// structure without re-scanning.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatRegion {
    /// The union of `intervals`' spans.
    pub span: RegionSpan,
    /// The overlapping intervals merged into this region — non-empty, coordinate-ordered.
    pub intervals: Box<[RepeatInterval]>,
}

/// One tile of the reference produced by [`RegionScanner`]. The tiles a scan yields are
/// coordinate-ordered, pairwise non-overlapping, and cover the scanned span exactly;
/// consecutive tiles never share a kind (each is maximal).
///
/// The three kinds are the caller's three routes: `Repeat` → the STR path, `Unique` →
/// the SNP path, and `Satellite` → neither (mask/skip).
///
/// (`RegionScanner` lands in Milestone C; this enum is its yielded item, defined here
/// with the rest of the vocabulary.)
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Region {
    /// A genotypeable tandem repeat (coverage ≤ `max_repeat_len`).
    Repeat(RepeatRegion),
    /// Repeat coverage longer than `max_repeat_len` — satellite DNA. Repetitive (so not
    /// unique sequence) but too long to be a genotypeable STR, hence excluded from STR
    /// analysis.
    Satellite(RegionSpan),
    /// No tandem-repeat coverage — unique sequence.
    Unique(RegionSpan),
}

/// A region-seam iteration error. Only a reference-read failure surfaces here — the
/// detection itself is infallible; a mid-scan [`ChromRefFetchError`] (corrupt/truncated
/// reference) is fatal to the scan and fuses the iterator (spec §8). Construction-time
/// validation of the period range is [`PeriodRangeError`], separately.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum ScanError {
    /// Reading the reference failed mid-scan.
    #[error("reading the reference failed mid-scan")]
    Fetch {
        #[source]
        source: ChromRefFetchError,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- PeriodRange -------------------------------------------------------

    #[test]
    fn period_range_accepts_valid_ranges() {
        assert_eq!(
            PeriodRange::new(1, 6).map(|r| (r.min(), r.max())),
            Ok((1, 6))
        );
        // A single-period range is valid.
        assert_eq!(
            PeriodRange::new(6, 6).map(|r| (r.min(), r.max())),
            Ok((6, 6))
        );
        assert_eq!(
            PeriodRange::new(2, 6).map(|r| (r.min(), r.max())),
            Ok((2, 6))
        );
    }

    #[test]
    fn period_range_rejects_zero_min() {
        assert_eq!(PeriodRange::new(0, 6), Err(PeriodRangeError::ZeroMin));
        // Zero-min is checked before the min>max rule, so (0, 0) is ZeroMin, not empty.
        assert_eq!(PeriodRange::new(0, 0), Err(PeriodRangeError::ZeroMin));
    }

    #[test]
    fn period_range_rejects_min_above_max() {
        assert_eq!(
            PeriodRange::new(3, 2),
            Err(PeriodRangeError::MinExceedsMax { min: 3, max: 2 })
        );
    }

    // ---- Defaults ----------------------------------------------------------

    #[test]
    fn default_periods_is_one_to_six() {
        assert_eq!(DEFAULT_PERIODS, (1, 6));
        // The named default range is itself a valid PeriodRange.
        assert!(PeriodRange::new(DEFAULT_PERIODS.0, DEFAULT_PERIODS.1).is_ok());
    }

    #[test]
    fn scan_params_default_is_two_seven_two() {
        let p = ScanParams::default();
        assert_eq!(p.match_reward, 2);
        assert_eq!(p.mismatch_penalty, 7);
        assert_eq!(p.min_copies, 2);
    }

    #[test]
    fn segment_options_default_is_cap_window_no_smoothing() {
        let o = SegmentOptions::default();
        assert_eq!(o.max_repeat_len, 1000);
        assert_eq!(o.window_bp, 100_000);
        assert_eq!(o.merge_gap, 0);
        assert_eq!(o.min_repeat_len, 0);
    }
}
