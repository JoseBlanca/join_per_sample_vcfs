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
//! Build status (incremental, per the impl plan): **Milestones A–B done** — the type
//! vocabulary and the `find_tandem_repeats` interval finder (lag-`p` scoring + a
//! Ruzzo–Tompa maximal-scoring-segment pass). The `RegionScanner` region seam
//! (Milestone C) lands as a later increment.
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

// ---------------------------------------------------------------------
// The interval finder (Milestone B) — lag-p scoring + Ruzzo–Tompa
// ---------------------------------------------------------------------

/// Canonicalise an ACGT base to an uppercase byte, or `None` for any non-ACGT byte.
/// A tandem-repeat match requires **both** compared bases to canonicalise equal (so an
/// `N == N` never matches — the two `None`s are not equal by this function's use).
#[inline]
fn canonical_base(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(b'A'),
        b'C' | b'c' => Some(b'C'),
        b'G' | b'g' => Some(b'G'),
        b'T' | b't' => Some(b'T'),
        _ => None,
    }
}

/// One open segment in the Ruzzo–Tompa working stack: `l`/`r` are the cumulative score
/// strictly before the segment and through its end (so the segment's own score is
/// `r - l`); `start`/`end` are inclusive indices into the score sequence.
struct RtSeg {
    l: i64,
    r: i64,
    start: usize,
    end: usize,
}

/// Find every **maximal scoring subsequence** of a score sequence (Ruzzo & Tompa 1999),
/// calling `emit(start, end, score)` once per segment — `start`/`end` inclusive indices
/// into the input, `score = r - l > 0`. Segments are emitted in ascending `start` order.
///
/// A maximal scoring subsequence is a contiguous run of positive total that cannot be
/// extended or trimmed without lowering its score; substitutions and short interruptions
/// stay inside one segment as long as the surrounding score outweighs them (spec §3.4).
///
/// Runs online: when a new positive element finds **no** stacked segment with a smaller
/// `l` (Ruzzo–Tompa "rule 2"), every stacked segment is provably unabsorbable by any
/// future element — a later absorption chain hits the new element's `l` floor and stops —
/// so the stack is flushed and reset. On the sparse-positive signal this scanner produces
/// the working stack therefore stays small (≈ open segments, not O(n)); the emitted set is
/// identical to offline Ruzzo–Tompa (pinned by a brute-force property test).
fn maximal_scoring_subsequences(
    scores: impl Iterator<Item = i64>,
    mut emit: impl FnMut(usize, usize, i64),
) {
    let mut stack: Vec<RtSeg> = Vec::new();
    let mut total: i64 = 0;
    for (i, s) in scores.enumerate() {
        let l = total;
        total += s;
        let r = total;
        if s <= 0 {
            continue; // non-positive scores never start a segment; they advance the total
        }
        let mut cur = RtSeg {
            l,
            r,
            start: i,
            end: i,
        };
        loop {
            // The rightmost stacked segment whose `l` is smaller than `cur.l`.
            match stack.iter().rposition(|seg| seg.l < cur.l) {
                None => {
                    // Rule 2: no left-smaller segment → flush the (unabsorbable) stack.
                    for seg in stack.drain(..) {
                        emit(seg.start, seg.end, seg.r - seg.l);
                    }
                    break;
                }
                Some(j) => {
                    if stack[j].r >= cur.r {
                        break; // Rule 3: `cur` stays a separate segment to the right.
                    }
                    // Rule 4: `cur` absorbs `stack[j..]` — extend its left edge and re-search.
                    cur.l = stack[j].l;
                    cur.start = stack[j].start;
                    stack.truncate(j);
                }
            }
        }
        stack.push(cur);
    }
    for seg in stack.drain(..) {
        emit(seg.start, seg.end, seg.r - seg.l);
    }
}

/// Find every tandem-repeat interval in `seq` whose period lies in `periods` (arch §2.1).
///
/// For each period `p` it scores position `j` (`j >= p`) `+match_reward` when `seq[j]`
/// and `seq[j-p]` canonicalise to the same ACGT base and `-mismatch_penalty` otherwise
/// (a non-ACGT base never matches), takes every Ruzzo–Tompa maximal scoring segment, maps
/// it to the 0-based half-open tract it certifies, and emits it when its implied copy
/// count clears `params.min_copies`. Pure and total over arbitrary bytes; deterministic.
///
/// Returns raw, possibly-overlapping intervals — one region can match at several periods,
/// and the finder does **not** de-duplicate periods or resolve overlaps; that is a
/// consumer's job (spec §1, §3.5). Per-period results are start-sorted; across periods
/// they are concatenated (period-ascending).
pub fn find_tandem_repeats(
    seq: &[u8],
    periods: PeriodRange,
    params: &ScanParams,
) -> Vec<RepeatInterval> {
    let mut out = Vec::new();
    let n = seq.len();
    for period in periods.min()..=periods.max() {
        let p = period as usize;
        if p >= n {
            continue; // no position has a partner `p` bases back
        }
        // Score index `k` (0-based over `p..n`) corresponds to position `j = k + p`.
        let scores = (p..n).map(
            |j| match (canonical_base(seq[j]), canonical_base(seq[j - p])) {
                (Some(a), Some(b)) if a == b => i64::from(params.match_reward),
                _ => -i64::from(params.mismatch_penalty),
            },
        );
        maximal_scoring_subsequences(scores, |k0, k1, score| {
            // Segment [k0, k1] → tract [k0, k1 + p + 1): the earliest base involved is
            // `j0 - p = k0`, the latest is `j1 = k1 + p`, so the exclusive end is `k1+p+1`.
            let tract_start = k0 as u32;
            let tract_end = (k1 + p + 1) as u32;
            let copies = (tract_end - tract_start) / u32::from(period);
            if copies >= params.min_copies {
                out.push(RepeatInterval {
                    start: tract_start,
                    end: tract_end,
                    period,
                    score: i32::try_from(score).unwrap_or(i32::MAX),
                });
            }
        });
    }
    out
}

// ---------------------------------------------------------------------
// The region seam (Milestone C) — repeat / satellite / unique tiling
// ---------------------------------------------------------------------

/// Tile a resident sequence into an ordered, gap-free run of [`Region`]s (spec §3.6): the
/// coverage-merge core shared by [`RegionScanner::over_slice`] and (later) the windowed
/// streaming path — so it is also the oracle for window-count invariance.
///
/// Steps: run [`find_tandem_repeats`], sort by start, union overlapping/abutting intervals
/// into merged repeat spans (grouping their intervals), apply the [`SegmentOptions`]
/// smoothing (`merge_gap` bridges short unique gaps, `min_repeat_len` drops short repeat
/// blips back to unique), classify each surviving span `Repeat` vs `Satellite` by
/// `max_repeat_len`, and emit `Unique` for every gap. The output tiles `[0, seq.len())`
/// exactly, in order, with no two same-kind tiles adjacent.
fn tile(
    seq: &[u8],
    periods: PeriodRange,
    params: &ScanParams,
    opts: &SegmentOptions,
) -> Vec<Region> {
    let n = seq.len() as u32;
    let mut out = Vec::new();
    if n == 0 {
        return out;
    }

    // 1. Detect, then sort globally by start (find_tandem_repeats concatenates per-period).
    let mut intervals = find_tandem_repeats(seq, periods, params);
    intervals.sort_by_key(|iv| (iv.start, iv.end));

    // 2. Union overlapping/abutting intervals into merged repeat spans, keeping the
    //    constituent intervals of each.
    struct Merged {
        start: u32,
        end: u32,
        intervals: Vec<RepeatInterval>,
    }
    let mut spans: Vec<Merged> = Vec::new();
    for iv in intervals {
        match spans.last_mut() {
            Some(last) if iv.start <= last.end => {
                last.end = last.end.max(iv.end);
                last.intervals.push(iv);
            }
            _ => spans.push(Merged {
                start: iv.start,
                end: iv.end,
                intervals: vec![iv],
            }),
        }
    }

    // 3a. Smoothing — bridge unique gaps shorter than `merge_gap` (coalesce the spans).
    if opts.merge_gap > 0 {
        let mut bridged: Vec<Merged> = Vec::new();
        for span in spans {
            match bridged.last_mut() {
                Some(last) if span.start - last.end < opts.merge_gap => {
                    last.end = span.end;
                    last.intervals.extend(span.intervals);
                }
                _ => bridged.push(span),
            }
        }
        spans = bridged;
    }

    // 3b. Smoothing — drop repeat spans shorter than `min_repeat_len` (they become unique).
    if opts.min_repeat_len > 0 {
        spans.retain(|s| s.end - s.start >= opts.min_repeat_len);
    }

    // 4. Emit the tiling: Unique for each gap, Repeat/Satellite for each span.
    let mut pos = 0u32;
    for span in spans {
        if span.start > pos {
            out.push(Region::Unique(RegionSpan {
                start: pos,
                end: span.start,
            }));
        }
        let region_span = RegionSpan {
            start: span.start,
            end: span.end,
        };
        if span.end - span.start > opts.max_repeat_len {
            out.push(Region::Satellite(region_span));
        } else {
            out.push(Region::Repeat(RepeatRegion {
                span: region_span,
                intervals: span.intervals.into_boxed_slice(),
            }));
        }
        pos = span.end;
    }
    if pos < n {
        out.push(Region::Unique(RegionSpan { start: pos, end: n }));
    }
    out
}

/// Streams a repeat/satellite/unique tiling of a reference as an iterator of [`Region`]s —
/// the routing seam the snp/ssr caller consumes (spec §3.6). Yields `Result<Region,
/// ScanError>`; the resident [`RegionScanner::over_slice`] path is infallible (always
/// `Ok`), while the windowed streaming path (a later increment) surfaces a mid-scan
/// reference-read failure as a terminal `Err`.
pub struct RegionScanner {
    regions: std::vec::IntoIter<Region>,
}

impl RegionScanner {
    /// Tile a small, already-resident sequence (spec §3.6). Precomputes the whole tiling;
    /// the windowed, memory-bounded constructor over a `ChromRefFetcher` lands separately.
    pub fn over_slice(
        seq: &[u8],
        periods: PeriodRange,
        params: &ScanParams,
        opts: &SegmentOptions,
    ) -> Self {
        Self {
            regions: tile(seq, periods, params, opts).into_iter(),
        }
    }
}

impl Iterator for RegionScanner {
    type Item = Result<Region, ScanError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.regions.next().map(Ok)
    }
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

    // ---- Ruzzo–Tompa vs a brute-force oracle -------------------------------

    /// Collect `maximal_scoring_subsequences` into a start-sorted `(start, end, score)`
    /// vector for comparison.
    fn rt_segments(scores: &[i64]) -> Vec<(usize, usize, i64)> {
        let mut out = Vec::new();
        maximal_scoring_subsequences(scores.iter().copied(), |a, b, s| out.push((a, b, s)));
        out.sort_unstable();
        out
    }

    /// Brute-force all maximal scoring subsequences straight from the Ruzzo–Tompa
    /// definition (O(n⁴), for small test inputs only — the independent oracle). An
    /// inclusive range `[i, j]` qualifies iff: (1) its score is positive; (2) every
    /// *proper* sub-range scores strictly less (so it is internally optimal); and (3) it
    /// is not properly contained in another range satisfying (1)+(2) (so it is maximal).
    fn brute_segments(scores: &[i64]) -> Vec<(usize, usize, i64)> {
        let n = scores.len();
        let mut pre = vec![0i64; n + 1];
        for i in 0..n {
            pre[i + 1] = pre[i] + scores[i];
        }
        let score = |i: usize, j: usize| pre[j + 1] - pre[i]; // inclusive [i, j]

        // Candidates: positive score, every proper sub-range strictly smaller.
        let mut cands: Vec<(usize, usize, i64)> = Vec::new();
        for i in 0..n {
            for j in i..n {
                let sc = score(i, j);
                if sc <= 0 {
                    continue;
                }
                let internally_optimal =
                    (i..=j).all(|i2| (i2..=j).all(|j2| (i2, j2) == (i, j) || score(i2, j2) < sc));
                if internally_optimal {
                    cands.push((i, j, sc));
                }
            }
        }

        // Keep only those not properly contained in another candidate (maximality).
        let mut res: Vec<(usize, usize, i64)> = cands
            .iter()
            .copied()
            .filter(|&(i, j, _)| {
                !cands
                    .iter()
                    .any(|&(a, b, _)| a <= i && j <= b && (a, b) != (i, j))
            })
            .collect();
        res.sort_unstable();
        res
    }

    #[test]
    fn ruzzo_tompa_matches_brute_force_on_crafted_cases() {
        let cases: &[&[i64]] = &[
            &[],
            &[-7, -7, -7],
            &[2, 2, 2, 2],
            &[2, -7, 2],  // two singletons — the interruption is too costly to bridge
            &[2, -3, 2],  // one merged segment — the dip is bridged
            &[5, -10, 8], // dip below the left L → two segments
            &[5, -3, 8],  // dip stays above → one segment
            &[2, 2, -7, 2, 2],
            &[-7, 2, 2, -7, 2, -7, 2, 2, -7],
            &[1, 1, 1, -2, 1, 1, 1, -10, 5, 5],
        ];
        for case in cases {
            assert_eq!(rt_segments(case), brute_segments(case), "case {case:?}");
        }
    }

    #[test]
    fn ruzzo_tompa_matches_brute_force_on_pseudo_random() {
        // Deterministic LCG over ±-style scores; no rand dependency.
        let mut state: u64 = 0x9E3779B97F4A7C15;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        for _ in 0..400 {
            let len = (next() % 24) as usize;
            let scores: Vec<i64> = (0..len)
                .map(|_| if next() % 4 == 0 { 2 } else { -7 }) // 25% match-like
                .collect();
            assert_eq!(
                rt_segments(&scores),
                brute_segments(&scores),
                "scores {scores:?}"
            );
        }
    }

    // ---- find_tandem_repeats ----------------------------------------------

    fn p(min: u8, max: u8) -> PeriodRange {
        PeriodRange::new(min, max).unwrap()
    }

    /// A clean (CAG)*8 tract detected at period 3 as one exact interval.
    #[test]
    fn finds_a_clean_perfect_tract() {
        let seq = b"CAGCAGCAGCAGCAGCAGCAGCAG"; // 8 copies, 24 bp
        let got = find_tandem_repeats(seq, p(3, 3), &ScanParams::default());
        assert_eq!(
            got,
            vec![RepeatInterval {
                start: 0,
                end: 24,
                period: 3,
                score: 2 * (24 - 3), // match_reward × matches
            }]
        );
    }

    /// A single substitution inside a tract keeps it **one** interval (spec §3.4).
    #[test]
    fn a_substitution_keeps_one_interval() {
        // (CAG)*5, then one base flipped, then (CAG)*5.
        let mut seq = Vec::new();
        for _ in 0..5 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend_from_slice(b"CAT"); // the interrupting copy (G→T)
        for _ in 0..5 {
            seq.extend_from_slice(b"CAG");
        }
        let got = find_tandem_repeats(&seq, p(3, 3), &ScanParams::default());
        let period3: Vec<_> = got.iter().filter(|r| r.period == 3).collect();
        assert_eq!(
            period3.len(),
            1,
            "one span across the substitution: {got:?}"
        );
        assert_eq!(period3[0].start, 0);
        assert_eq!(period3[0].end, seq.len() as u32);
    }

    /// A single indel (a bounded mismatch burst) also keeps the tract **one** interval —
    /// the lag-p comparison is local, so downstream copies resume in phase (spec §3.4).
    #[test]
    fn a_single_indel_keeps_one_interval() {
        // (CAG)*6, one inserted base, (CAG)*6 — phase shifts but the burst is bounded.
        let mut seq = Vec::new();
        for _ in 0..6 {
            seq.extend_from_slice(b"CAG");
        }
        seq.push(b'T'); // a 1 bp insertion
        for _ in 0..6 {
            seq.extend_from_slice(b"CAG");
        }
        let got = find_tandem_repeats(&seq, p(3, 3), &ScanParams::default());
        let period3: Vec<_> = got.iter().filter(|r| r.period == 3).collect();
        assert_eq!(
            period3.len(),
            1,
            "the single indel is a bounded burst, not a split: {got:?}"
        );
        assert_eq!(period3[0].start, 0);
        assert_eq!(period3[0].end, seq.len() as u32);
    }

    /// A run of `N` yields nothing — `N == N` is not a match (spec §3.5).
    #[test]
    fn an_n_run_emits_nothing() {
        let seq = b"NNNNNNNNNNNNNNNNNNNN";
        assert!(find_tandem_repeats(seq, p(1, 6), &ScanParams::default()).is_empty());
    }

    /// A soft-masked (lowercase) tract is still found — comparison is case-insensitive.
    #[test]
    fn a_soft_masked_tract_is_found() {
        let seq = b"cagcagcagcagcagcagcagcag"; // lowercase (CAG)*8
        let got = find_tandem_repeats(seq, p(3, 3), &ScanParams::default());
        assert_eq!(got.len(), 1);
        assert_eq!((got[0].start, got[0].end, got[0].period), (0, 24, 3));
    }

    /// An (AT)*n tract matches at period 2 **and** its multiples 4 and 6 — all emitted;
    /// the finder does not de-duplicate periods (spec §3.5).
    #[test]
    fn multiples_of_the_fundamental_period_are_all_emitted() {
        let seq = b"ATATATATATATATATATAT"; // (AT)*10, 20 bp
        let got = find_tandem_repeats(seq, p(2, 6), &ScanParams::default());
        let periods: std::collections::BTreeSet<u8> = got.iter().map(|r| r.period).collect();
        assert!(
            periods.contains(&2),
            "fundamental period 2 present: {got:?}"
        );
        assert!(periods.contains(&4), "multiple period 4 present: {got:?}");
        assert!(periods.contains(&6), "multiple period 6 present: {got:?}");
        assert!(
            !periods.contains(&3) && !periods.contains(&5),
            "no match at non-divisor periods 3/5: {got:?}"
        );
    }

    /// The `min_copies` floor drops sub-threshold segments (e.g. a lone coincidental match).
    #[test]
    fn the_min_copies_floor_drops_short_segments() {
        // One coincidental period-2 match embedded in non-repetitive sequence.
        let seq = b"ACGTACAGTCTGAC";
        let got = find_tandem_repeats(seq, p(2, 2), &ScanParams::default());
        assert!(
            got.iter().all(|r| (r.end - r.start) / 2 >= 2),
            "every emitted interval has >= min_copies copies: {got:?}"
        );
    }

    // ---- RegionScanner::over_slice (the region seam, C1) -------------------

    /// A region's `(start, end, kind)` for the tiling-contract check — kind is 0 Repeat,
    /// 1 Satellite, 2 Unique.
    fn triple(r: &Region) -> (u32, u32, u8) {
        match r {
            Region::Repeat(rr) => (rr.span.start, rr.span.end, 0),
            Region::Satellite(s) => (s.start, s.end, 1),
            Region::Unique(s) => (s.start, s.end, 2),
        }
    }

    fn regions_of(seq: &[u8], periods: PeriodRange, opts: &SegmentOptions) -> Vec<Region> {
        RegionScanner::over_slice(seq, periods, &ScanParams::default(), opts)
            .map(|r| r.unwrap())
            .collect()
    }

    /// Every tiling is coordinate-ordered, gap-free over `[0, n)`, non-overlapping, and has
    /// no two same-kind tiles adjacent.
    fn assert_tiles(regions: &[Region], n: u32) {
        let mut pos = 0u32;
        let mut prev_kind: Option<u8> = None;
        for r in regions {
            let (s, e, kind) = triple(r);
            assert_eq!(
                s, pos,
                "tile starts at {s}, expected {pos} (gap or overlap)"
            );
            assert!(e > s, "empty tile [{s}, {e})");
            assert_ne!(Some(kind), prev_kind, "two same-kind tiles adjacent at {s}");
            prev_kind = Some(kind);
            pos = e;
        }
        assert_eq!(pos, n, "tiling must cover exactly [0, {n})");
    }

    #[test]
    fn empty_input_yields_no_regions() {
        assert!(regions_of(b"", p(1, 6), &SegmentOptions::default()).is_empty());
    }

    #[test]
    fn a_lone_repeat_tiles_as_unique_repeat_unique() {
        // Non-repetitive flanks around a clean (CAG)*8 tract.
        let mut seq = b"TTGGA".to_vec();
        for _ in 0..8 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend_from_slice(b"GGTTA");
        let regions = regions_of(&seq, p(3, 3), &SegmentOptions::default());
        assert_tiles(&regions, seq.len() as u32);
        let kinds: Vec<u8> = regions.iter().map(|r| triple(r).2).collect();
        assert_eq!(kinds, vec![2, 0, 2], "Unique, Repeat, Unique: {regions:?}");
        // The single Repeat covers the tract.
        if let Region::Repeat(rr) = &regions[1] {
            assert_eq!((rr.span.start, rr.span.end), (5, 5 + 24));
        } else {
            panic!("middle tile is the repeat");
        }
    }

    #[test]
    fn overlapping_periods_merge_into_one_repeat_carrying_all() {
        // (AT)*10 matches at periods 2, 4 and 6 — three overlapping intervals, one Repeat.
        let seq = b"ATATATATATATATATATAT"; // 20 bp
        let regions = regions_of(seq, p(2, 6), &SegmentOptions::default());
        assert_tiles(&regions, seq.len() as u32);
        let repeats: Vec<_> = regions
            .iter()
            .filter_map(|r| match r {
                Region::Repeat(rr) => Some(rr),
                _ => None,
            })
            .collect();
        assert_eq!(repeats.len(), 1, "one merged repeat region: {regions:?}");
        let periods: std::collections::BTreeSet<u8> =
            repeats[0].intervals.iter().map(|iv| iv.period).collect();
        assert!(
            periods.contains(&2) && periods.contains(&4) && periods.contains(&6),
            "the merged region carries all overlapping periods: {:?}",
            repeats[0].intervals
        );
    }

    #[test]
    fn a_repeat_over_the_cap_is_a_satellite() {
        // The same (CAG)*8 = 24 bp tract, but with a 10 bp satellite cap → Satellite.
        let mut seq = b"TTGGA".to_vec();
        for _ in 0..8 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend_from_slice(b"GGTTA");
        let opts = SegmentOptions {
            max_repeat_len: 10,
            ..SegmentOptions::default()
        };
        let regions = regions_of(&seq, p(3, 3), &opts);
        assert_tiles(&regions, seq.len() as u32);
        assert_eq!(
            regions.iter().map(|r| triple(r).2).collect::<Vec<_>>(),
            vec![2, 1, 2],
            "Unique, Satellite, Unique: {regions:?}"
        );
    }

    #[test]
    fn merge_gap_bridges_two_nearby_repeats() {
        // Two (CAG)*5 tracts separated by a 12 bp non-repetitive gap → two intervals.
        let mut seq = Vec::new();
        for _ in 0..5 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend_from_slice(b"TGACGTTGACGT"); // 12 bp gap the finder won't bridge
        for _ in 0..5 {
            seq.extend_from_slice(b"CAG");
        }
        // Off: two separate Repeats, a Unique gap between them.
        let off = regions_of(&seq, p(3, 3), &SegmentOptions::default());
        let repeats_off = off.iter().filter(|r| triple(r).2 == 0).count();
        assert_eq!(repeats_off, 2, "without merge_gap, two repeats: {off:?}");
        assert_tiles(&off, seq.len() as u32);
        // On: a merge_gap wider than the 12 bp gap coalesces them into one Repeat.
        let opts = SegmentOptions {
            merge_gap: 20,
            ..SegmentOptions::default()
        };
        let on = regions_of(&seq, p(3, 3), &opts);
        assert_eq!(
            on.iter().filter(|r| triple(r).2 == 0).count(),
            1,
            "with merge_gap, one bridged repeat: {on:?}"
        );
        assert_tiles(&on, seq.len() as u32);
    }

    #[test]
    fn min_repeat_len_reclassifies_a_short_blip_as_unique() {
        // An isolated (CAG)*3 = 9 bp tract, flanked.
        let mut seq = b"TTGGA".to_vec();
        for _ in 0..3 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend_from_slice(b"GGTTA");
        // Off: the 9 bp tract is a Repeat.
        let off = regions_of(&seq, p(3, 3), &SegmentOptions::default());
        assert!(
            off.iter().any(|r| triple(r).2 == 0),
            "a repeat is present: {off:?}"
        );
        // On: min_repeat_len = 20 drops it → the whole contig is one Unique tile.
        let opts = SegmentOptions {
            min_repeat_len: 20,
            ..SegmentOptions::default()
        };
        let on = regions_of(&seq, p(3, 3), &opts);
        assert_eq!(
            on.iter().map(|r| triple(r).2).collect::<Vec<_>>(),
            vec![2],
            "the short repeat is reclassified unique, gaps merged: {on:?}"
        );
        assert_tiles(&on, seq.len() as u32);
    }
}
