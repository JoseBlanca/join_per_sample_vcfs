//! A clean-room, use-agnostic short-period **tandem-repeat scanner** — a sequence
//! primitive shared by the STR catalog and (later) the ng snp/str caller. Design:
//! spec [`doc/devel/ng/spec/ssr_repeat_scanner.md`], arch
//! [`doc/devel/ng/arch/ssr_repeat_scanner.md`].
//!
//! It exposes **three interfaces** over one detection core (a lag-`p` self-comparison
//! plus a Ruzzo–Tompa maximal-scoring-segment pass):
//!
//! - a low-level **interval finder** — `find_tandem_repeats(seq, PeriodRange,
//!   &ScanParams) -> Vec<RepeatInterval>` — raw, possibly-overlapping intervals, for
//!   consumers that resolve overlaps themselves (the STR catalog); and
//! - a mid-level **windowed scan** — [`scan_windowed`], which yields one
//!   [`ScannedWindow`] at a time (the core + margin geometry, every detection in the
//!   fetched slice, and the two window rules — coverage clipped to cores, intervals
//!   attributed by start — as methods the consumer applies to *its own* interval set)
//!   and takes **no routing decision**; for a consumer that must inject its own policy
//!   between detection and any merge — the typed-region generator, which is why this
//!   layer is public (`typed_regions.md` §6.1); and
//! - a high-level **region seam** — `RegionScanner`, a streaming iterator that yields
//!   a gap-free repeat/satellite/unique tiling of the reference, built on the above.
//!
//! The scanner holds **no consumer policy** (no purity floor, homopolymer rule, or
//! period ceiling): a consumer passes the period range it wants, and that plus the two
//! scoring weights is the whole configuration surface.
//!
//! **Clean-room note:** derived only from Benson (1999)'s published idea (compare to a
//! period-shifted copy; grow runs of matches) and our own code — **not** from the
//! AGPL-v3 `TRF-mod` source, which was not read.
//!
//! Build status (incremental, per the impl plan): **the scanner plan's Milestones A–D are
//! done** — the type vocabulary, the `find_tandem_repeats` interval finder (lag-`p`
//! scoring + a Ruzzo–Tompa maximal-scoring-segment pass), the `RegionScanner` region seam
//! (coverage-merge tiling, resident via `over_slice` and windowed-streaming via `stream`),
//! and trf-mod parity. **Plus typed-regions B1**: [`scan_windowed`] promoted from a private
//! whole-contig-eager helper to a public per-window stream (`typed_regions.md` §6.1).
//!
//! Visibility note: types are `pub` (matching the sibling ng modules, e.g. `types.rs`),
//! not the arch doc's illustrative `pub(crate)` — the ng convention exposes its module
//! vocabulary as reachable API, which also lets Milestone-A types be defined ahead of
//! their Milestone-B/C consumers without tripping `dead_code`.

use crate::fasta::{ChromRefFetchError, ChromRefFetcher};

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
    ///
    /// `u32`, not `u64`: this is a **count**, not a coordinate or a length, so
    /// spec §4's widening does not reach it (the same reason ids stay `u32`, and
    /// matching `region_typing::admission::MinCopies`).
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
pub const DEFAULT_MAX_REPEAT_LEN: u64 = 1000;
/// Default streaming window core (bp): a memory/I/O-granularity knob only — the region
/// tiling is invariant to it.
pub const DEFAULT_WINDOW_BP: u64 = 100_000;
/// Default unique-gap bridging (0 = off): smoothing for the region tiling.
pub const DEFAULT_MERGE_GAP: u64 = 0;
/// Default minimum repeat-region length (0 = off): smoothing for the region tiling.
pub const DEFAULT_MIN_REPEAT_LEN: u64 = 0;

/// Region-seam shaping and streaming knobs — separate from [`ScanParams`] because they
/// shape the **partition and the walk**, not the detection.
///
/// `Default` = `{ 1000, 100_000, 0, 0 }`: the 1 kb satellite cap, a 100 kb streaming
/// window, and both smoothing knobs off (a pure coverage partition).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SegmentOptions {
    /// Merged repeat coverage longer than this is emitted as [`Region::Satellite`].
    pub max_repeat_len: u64,
    /// The streaming window core in bp; a memory knob, region-count-invariant.
    pub window_bp: u64,
    /// Smoothing: bridge unique gaps shorter than this into the flanking repeat.
    pub merge_gap: u64,
    /// Smoothing: reclassify a repeat region shorter than this as unique.
    pub min_repeat_len: u64,
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
///
/// **0-based, deliberately, where the rest of ng is 1-based** (spec §4 says 1-based
/// *everywhere*; this is the reasoned exception, owner-agreed 2026-07-16 at B2).
/// This is a **detector/slice type**: its coordinates index the byte slice
/// `find_tandem_repeats` was handed, and `seq[iv.start..iv.end]` is what every
/// producer and consumer of it actually does. Making it 1-based would put a `- 1`
/// at every one of those sites to buy consistency with a *genetic* coordinate it
/// never expresses. The 1-based boundary belongs where the genetics starts —
/// `Locus`, `GenomeRegion` — which is exactly where §4's payoff lands. `u64` is
/// the width regardless.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RepeatInterval {
    /// Tract start, inclusive (0-based).
    pub start: u64,
    /// Tract end, exclusive (0-based).
    pub end: u64,
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
    pub start: u64,
    /// End, exclusive (0-based).
    pub end: u64,
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
/// This is the textbook offline pass: a working stack of open segments, drained at the end.
/// Its correctness is pinned by a brute-force property test (`ruzzo_tompa_matches_brute_force`).
/// The stack is O(number of open segments) — bounded per window on the region seam's streaming
/// path; on a whole-contig `find_tandem_repeats` (the catalog's use) it can grow to O(n) on a
/// pathological input, an online-finalisation optimisation deferred behind that property test.
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
                None => break, // Rule 2: no left-smaller segment → `cur` is separate.
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
            let tract_start = k0 as u64;
            let tract_end = (k1 + p + 1) as u64;
            let copies = (tract_end - tract_start) / u64::from(period);
            if copies >= u64::from(params.min_copies) {
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

/// A merged repeat span (the union of overlapping/abutting intervals) plus the intervals
/// that formed it — the intermediate between raw intervals and emitted [`Region`]s.
/// Union a bag of half-open `[start, end)` spans into disjoint, start-sorted spans. The
/// repeat **coverage** — the positions some repeat covers — independent of which intervals
/// produced it (so it merges identically whether fed whole-contig intervals or windowed,
/// clipped fragments).
fn union_spans(mut spans: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    spans.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::new();
    for (s, e) in spans {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

/// Build the final region tiling of `[0, n)` from two inputs kept **deliberately separate**:
/// the repeat `coverage` (span-level, which classifies Repeat/Satellite and fixes the
/// boundaries) and the exact `str_intervals` (which decorate each `Repeat` with its
/// constituents). Splitting them is what makes the windowed path invariant: coverage
/// coalesces satellite fragments across windows, while STR intervals — always ≤
/// `max_repeat_len`, so captured whole by one window — are attached by start.
///
/// Steps: union the coverage into spans, apply the [`SegmentOptions`] smoothing (`merge_gap`
/// bridges short unique gaps, `min_repeat_len` drops short repeat blips), classify each
/// surviving span `Repeat` vs `Satellite` by `max_repeat_len` (attaching the `str_intervals`
/// whose start lies in a `Repeat` span), and fill every gap with `Unique`. The output tiles
/// `[0, n)` exactly, coordinate-ordered, with no two same-kind tiles adjacent.
fn build_regions(
    coverage: Vec<(u64, u64)>,
    mut str_intervals: Vec<RepeatInterval>,
    n: u64,
    opts: &SegmentOptions,
) -> Vec<Region> {
    let mut spans = union_spans(coverage);

    // Smoothing — bridge unique gaps shorter than `merge_gap`.
    if opts.merge_gap > 0 {
        let mut bridged: Vec<(u64, u64)> = Vec::new();
        for (s, e) in spans {
            match bridged.last_mut() {
                Some(last) if s - last.1 < opts.merge_gap => last.1 = e,
                _ => bridged.push((s, e)),
            }
        }
        spans = bridged;
    }
    // Smoothing — drop repeat spans shorter than `min_repeat_len` (they become unique).
    if opts.min_repeat_len > 0 {
        spans.retain(|&(s, e)| e - s >= opts.min_repeat_len);
    }

    str_intervals.sort_by_key(|iv| (iv.start, iv.end));
    let mut iv_idx = 0usize; // a single forward cursor into the start-sorted intervals

    let mut out = Vec::new();
    let mut pos = 0u64;
    for (start, end) in spans {
        if start > pos {
            out.push(Region::Unique(RegionSpan {
                start: pos,
                end: start,
            }));
        }
        let region_span = RegionSpan { start, end };
        if end - start > opts.max_repeat_len {
            out.push(Region::Satellite(region_span));
        } else {
            // Attach the intervals whose start lies in this span. Skipped intervals (in
            // gaps, satellites, or dropped spans) are stepped over by the shared cursor.
            while iv_idx < str_intervals.len() && str_intervals[iv_idx].start < start {
                iv_idx += 1;
            }
            let lo = iv_idx;
            while iv_idx < str_intervals.len() && str_intervals[iv_idx].start < end {
                iv_idx += 1;
            }
            out.push(Region::Repeat(RepeatRegion {
                span: region_span,
                intervals: str_intervals[lo..iv_idx].to_vec().into_boxed_slice(),
            }));
        }
        pos = end;
    }
    if pos < n {
        out.push(Region::Unique(RegionSpan { start: pos, end: n }));
    }
    out
}

/// Tile a resident sequence into an ordered, gap-free run of [`Region`]s (spec §3.6). The
/// **oracle** the windowed streaming path (`stream`) must reproduce (window-count invariance):
/// its coverage is every interval's span and its STR intervals are all of them.
fn tile(
    seq: &[u8],
    periods: PeriodRange,
    params: &ScanParams,
    opts: &SegmentOptions,
) -> Vec<Region> {
    let n = seq.len() as u64;
    if n == 0 {
        return Vec::new();
    }
    let intervals = find_tandem_repeats(seq, periods, params);
    let coverage = intervals.iter().map(|iv| (iv.start, iv.end)).collect();
    build_regions(coverage, intervals, n, opts)
}

/// One window's worth of scan — what [`scan_windowed`] yields per step.
///
/// **The data is one field; the window *rules* are the two methods.** A window
/// yields every detection in the slice it fetched, and the two derivations that
/// make a windowed scan add up to a whole-contig one —
/// [`coverage_in_core`](Self::coverage_in_core) and
/// [`starting_in_core`](Self::starting_in_core) — are methods **taking the
/// intervals as an argument** rather than pre-computed fields.
///
/// That shape is what lets the two consumers differ where they must and agree
/// where they must not. `RegionScanner` applies the rules to the raw detections;
/// the typed-region walk applies them to its **pre-filtered** ones (`admission::
/// prefilter`, which the scanner must not know about — this module holds no
/// consumer policy). Pre-computed fields would have forced the walk to either
/// accept raw-derived coverage — capping detector noise, which spec §2.4 forbids —
/// or to re-implement the clipping, which is exactly the invariant duplication
/// `typed_regions.md` §6.1 says to avoid.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScannedWindow {
    /// The window's core, `[start, end)`, 0-based — **not** the fetched slice
    /// ([`fetched`](Self::fetched)), which is the core plus a `max_repeat_len`
    /// margin each side. Cores tile the contig exactly: each starts where the last
    /// ended, so core-clipped coverage from adjacent windows abuts and a later union
    /// can rejoin a satellite across every window it spans.
    pub core: RegionSpan,
    /// The slice actually fetched and scanned, `[start, end)`, 0-based: the core
    /// plus a `max_repeat_len` margin each side, **clamped to the contig**.
    ///
    /// Public because the walk fetches these same bases for itself — admission needs
    /// the sequence, which a scan result does not carry — and must not compute the
    /// range independently: that arithmetic is this module's rule, and two copies of
    /// it are two copies that can drift.
    pub fetched: RegionSpan,
    /// **Every** detection in the fetched slice, in contig coordinates — including
    /// those starting in the margin, which is what makes them useful: a consumer's
    /// set-wise policy (redundancy elimination, the bundle flank test) needs an
    /// in-core interval's *neighbours* to reach the same verdict a whole-contig run
    /// would, and the `max_repeat_len` margin is exactly the radius that guarantees
    /// it.
    ///
    /// **Two caveats, both handled by the methods** rather than by the reader:
    ///
    /// - a detection is whole *unless* it runs off the fetched slice, where it is
    ///   **truncated**. Harmless by construction: to be truncated it must be longer
    ///   than `max_repeat_len` (i.e. satellite coverage, which
    ///   [`coverage_in_core`](Self::coverage_in_core) rejoins across windows anyway),
    ///   and it cannot then start in this core, so
    ///   [`starting_in_core`](Self::starting_in_core) never hands one out as exact.
    /// - concatenating these across windows **double-counts** the margins. Use
    ///   [`starting_in_core`](Self::starting_in_core) to attribute each to exactly
    ///   one window.
    ///
    /// **Order: period-major, then ascending start** — inherited from
    /// [`find_tandem_repeats`], which scans one period at a time. It is **not**
    /// start-sorted, and neither is the concatenation across windows. A consumer
    /// that needs coordinate order must sort; `admit` and `prefilter` both do, and
    /// `build_regions` re-sorts. Stated because the type is public and the ordering
    /// is not what a reader would assume.
    pub detections: Vec<RepeatInterval>,
}

impl ScannedWindow {
    /// `intervals`, **clipped to this window's core** — the coverage rule.
    ///
    /// Cores tile the contig, so unioning this across every window reconstructs the
    /// whole-contig coverage of the same intervals, and a satellite (necessarily
    /// longer than a window's reach) rejoins across all the windows it spans. This
    /// is why a truncated detection is harmless: whatever it loses lies outside
    /// some core, and that core's own window contributes it.
    ///
    /// `intervals` is a parameter rather than `self.detections` because the caller's
    /// policy decides *which* intervals count — see the type's docs.
    pub fn coverage_in_core(&self, intervals: &[RepeatInterval]) -> Vec<(u64, u64)> {
        intervals
            .iter()
            .filter_map(|iv| {
                let start = iv.start.max(self.core.start);
                let end = iv.end.min(self.core.end);
                (start < end).then_some((start, end))
            })
            .collect()
    }

    /// The intervals this window **owns**: kept whole, and only those whose *start*
    /// lands in the core — the attribution rule.
    ///
    /// Cores tile and do not overlap, so every interval is owned by exactly one
    /// window, whichever holds its start. An owned interval is never truncated: it
    /// starts at least `max_repeat_len` inside the fetched slice on the left and,
    /// were it long enough to reach the right edge, it would be satellite coverage
    /// rather than an exact repeat.
    pub fn starting_in_core<'a>(
        &self,
        intervals: &'a [RepeatInterval],
    ) -> impl Iterator<Item = &'a RepeatInterval> {
        let core = self.core;
        intervals
            .iter()
            .filter(move |iv| iv.start >= core.start && iv.start < core.end)
    }
}

/// One window's geometry: the core it owns, and the slice that must be read to scan it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WindowPlan {
    /// The core, `[start, end)` 0-based — the bases this window is responsible for.
    pub core: RegionSpan,
    /// The core plus a `max_repeat_len` margin each side, clamped to the contig.
    pub fetched: RegionSpan,
}

/// Walks a range of a contig as windows — **the geometry rule, and the only copy of it.**
///
/// # Why this is a separate type from [`scan_windowed`]
///
/// Because its two consumers cannot share an iterator. `scan_windowed` yields windows it
/// fetched through a *borrowed* byte source; [`crate::ng::region_typing::TypedRegionIterator`]
/// **owns** its reference (it moves onto a producer thread, `typed_regions.md` §7) and so
/// cannot hold an iterator borrowing it — that is a self-referential struct, which safe
/// Rust does not have. A cursor it can step by hand is the shape that fits both, and
/// `typed_regions.md` §6.1 is explicit that the thing not to do is duplicate this
/// arithmetic ("a subtle invariant that can then drift").
///
/// # Cores tile `cores`; margins reach outside it
///
/// `cores` is the range to walk and `contig_len` bounds the **margins**, which is not the
/// same number the moment a caller walks part of a contig (a BED region — spec §2.5): the
/// window at the range's edge must still read its margin from the sequence *beyond* the
/// range, or the repeat lying across that edge is measured as a fragment.
pub struct WindowCursor {
    cores: RegionSpan,
    contig_len: u64,
    window: u64,
    margin: u64,
    next_core: u64,
}

impl WindowCursor {
    /// Walk `cores` (0-based `[start, end)`) as windows, taking margins from a contig
    /// `contig_len` long. `window_bp = 0` coerces to 1 — a zero-length core would step
    /// forever, and a hang is the one failure a `should_panic` test cannot catch.
    pub fn new(cores: RegionSpan, contig_len: u64, opts: &SegmentOptions) -> Self {
        Self {
            cores,
            contig_len,
            window: opts.window_bp.max(1),
            margin: opts.max_repeat_len,
            next_core: cores.start,
        }
    }

    /// The next window, or `None` at the end of the range.
    pub fn next_window(&mut self) -> Option<WindowPlan> {
        if self.next_core >= self.cores.end {
            return None;
        }
        let core = RegionSpan {
            start: self.next_core,
            end: self
                .next_core
                .saturating_add(self.window)
                .min(self.cores.end),
        };
        let fetched = RegionSpan {
            start: core.start.saturating_sub(self.margin),
            end: core.end.saturating_add(self.margin).min(self.contig_len),
        };
        self.next_core = core.end;
        Some(WindowPlan { core, fetched })
    }

    /// Stop the walk: the next [`Self::next_window`] yields `None`. For a caller that has
    /// hit a fatal error and must not scan on — continuing would scan a hole.
    pub fn halt(&mut self) {
        self.next_core = self.cores.end;
    }
}

/// Scan one window's already-fetched bases into a [`ScannedWindow`] — **the lift, and the
/// only copy of it**: detections come out of [`find_tandem_repeats`] as offsets into
/// `bases`, and everything downstream speaks contig coordinates.
///
/// `bases` must be exactly `plan.fetched`, which is [`WindowCursor`]'s to say and not the
/// caller's to work out.
pub fn scan_window(
    plan: &WindowPlan,
    bases: &[u8],
    periods: PeriodRange,
    params: &ScanParams,
) -> ScannedWindow {
    debug_assert_eq!(
        bases.len() as u64,
        plan.fetched.end - plan.fetched.start,
        "bases must be exactly the planned slice"
    );
    ScannedWindow {
        core: plan.core,
        fetched: plan.fetched,
        detections: find_tandem_repeats(bases, periods, params)
            .into_iter()
            .map(|iv| RepeatInterval {
                start: iv.start + plan.fetched.start,
                end: iv.end + plan.fetched.start,
                period: iv.period,
                score: iv.score,
            })
            .collect(),
    }
}

/// Scan a contig in windows through a `ChromRefFetcher`, **yielding one
/// [`ScannedWindow`] at a time** — memory-bounded (peak ≈ `window_bp +
/// 2·max_repeat_len`, not the contig length; spec §3.6).
///
/// # Why this is public, and why it streams
///
/// It is the primitive the **typed-region generator** stands on
/// (`typed_regions.md` §6.1): that walk cannot use [`RegionScanner`], because
/// `RegionScanner` merges coverage and classifies satellites **before any
/// admission policy can run**, on raw permissive intervals, with nowhere to inject
/// the pre-filter between detection and merge. The layer underneath it fits
/// exactly — this one — because it does the genuinely hard part (core + margin,
/// coverage clipped to cores, intervals attributed by start) and takes **no
/// routing decision at all**.
///
/// Streaming rather than returning two contig-wide `Vec`s is what makes it usable
/// there: the walk holds one window plus a few open coordinates, never a contig
/// (`typed_regions.md` §2.6, §7). `RegionScanner::stream` collects it right back
/// up, which is fine — it was always whole-contig-eager.
///
/// # Why the byte source is a closure and not a `ChromRefFetcher`
///
/// Because its two consumers need **different bytes**, and neither can have the
/// other's:
///
/// - `RegionScanner` reads through production's `ChromRefFetcher`, which
///   **canonicalises** (`{A,C,G,T,N}`) as it fills its buffer.
/// - The **typed-region walk** needs **raw** bytes: the STR catalog reads the FASTA
///   verbatim and `Locus` compares *by value*, so canonical bytes would make every
///   IUPAC-carrying locus compare unequal to the catalog's and silently break the
///   `.cat` oracle (`typed_regions.md` §6). That is what
///   [`crate::ng::ref_seq::WindowedRefSeq`]'s raw path exists for — and it is a
///   `RawRefSeq`, not a `ChromRefFetcher`, so it cannot be handed to a fetcher-bound
///   scan at all.
///
/// Binding this to either source would force the other to reimplement the
/// windowing — which `typed_regions.md` §6.1 names as exactly the thing not to do
/// (*"duplicates a subtle invariant that can then drift"*) — or to hold **two**
/// readers, which §6's 14.6 GB lesson forbids. A closure lets both pass their own
/// bytes through **one** copy of the window rules.
///
/// Only *admission* needs raw, incidentally: detection is case-insensitive and
/// non-ACGT never matches, so `find_tandem_repeats` gives the same answer either
/// way. That is why nothing before the walk tripped on this.
///
/// The error type is the **caller's** (`E`): this function never inspects it, only
/// propagates it, so `RegionScanner` can yield a `ScanError` and the walk a
/// `TypedRegionError` without either bending to the other. `contig_len` likewise
/// comes from the caller — the source knows its own length, and for the walk that
/// number must come from the reference's contig table (`typed_regions.md` §2.6's
/// provenance obligation), which a fetcher-derived length could not guarantee.
///
/// A fetch failure is yielded once as `Some(Err(_))`, and the iterator then stops:
/// continuing would silently scan a hole, which is the failure this module's
/// callers are least able to see. The stop is a latch ([`WindowCursor::halt`]),
/// **and the returned iterator is `FusedIterator`** — [`ScanError`]'s own doc
/// already promises this, and `std::iter::FromFn` does not provide it, so the
/// promise is declared rather than merely kept. Mirrors the ng sibling
/// `ReadFilter`, which states the same contract the same way.
///
/// This function is now a thin loop: [`WindowCursor`] holds the geometry and
/// [`scan_window`] the lift, so the walk — which cannot use this iterator at all, see
/// [`WindowCursor`] — still runs the same two rules rather than a second copy of them.
/// What it adds is the **fetch loop and the failure latch**, and that is the whole of it.
///
/// Each window fetches its core plus a `max_repeat_len` margin on **each** side and runs
/// `find_tandem_repeats` on that slice; the two derivations that make the result
/// window-invariant are [`ScannedWindow`]'s methods. The margins are load-bearing:
/// Ruzzo–Tompa segmentation is context-dependent, so a shorter margin would let a window
/// re-segment a long repeat mid-tract and emit a spurious in-core fragment.
pub fn scan_windowed<F, E>(
    contig_len: u64,
    periods: PeriodRange,
    params: &ScanParams,
    opts: &SegmentOptions,
    mut fetch: F,
) -> impl std::iter::FusedIterator<Item = Result<ScannedWindow, E>>
where
    F: FnMut(u64, u64, &mut Vec<u8>) -> Result<(), E>,
{
    // Both config structs are `Copy` and are read out here, so tying the stream to
    // their borrows would constrain it for nothing.
    let params = *params;
    let mut cursor = WindowCursor::new(
        RegionSpan {
            start: 0,
            end: contig_len,
        },
        contig_len,
        opts,
    );
    // One buffer, reused across windows — the caller's `fetch` fills it. This is
    // why `fetch` writes into a `&mut Vec<u8>` rather than returning one: a window
    // is ~102 kb, and re-allocating it per step would be the memory profile the
    // streaming shape exists to avoid.
    let mut bases: Vec<u8> = Vec::new();

    let scan = std::iter::from_fn(move || {
        let plan = cursor.next_window()?;
        // `fetch` is 1-based.
        if let Err(e) = fetch(
            plan.fetched.start + 1,
            plan.fetched.end - plan.fetched.start,
            &mut bases,
        ) {
            // Stop after reporting: a fetch failure is fatal, and continuing would
            // silently scan a hole.
            cursor.halt();
            return Some(Err(e));
        }
        Some(Ok(scan_window(&plan, &bases, periods, &params)))
    });
    scan.fuse()
}

/// Streams a repeat/satellite/unique tiling of a reference as an iterator of [`Region`]s —
/// the routing seam the snp/ssr caller consumes (spec §3.6). Yields `Result<Region,
/// ScanError>`.
///
/// Two constructors: [`over_slice`](Self::over_slice) tiles a resident sequence, and
/// [`stream`](Self::stream) tiles a whole contig windowed through a `ChromRefFetcher`
/// without holding it in memory. `stream` is **window-count invariant for the STR tiling** —
/// it yields exactly the `Repeat`/`Unique` structure `over_slice` would, regardless of
/// `window_bp`. A `Satellite` (a repeat longer than `max_repeat_len`, hence longer than a
/// window's reach) may segment differently between the two, because its internal structure
/// depends on global context a bounded window cannot see; that is harmless, since satellites
/// are mask/skip regions, not genotyping targets (see [`scan_windowed`]).
pub struct RegionScanner {
    regions: std::vec::IntoIter<Region>,
}

impl RegionScanner {
    /// Tile a small, already-resident sequence (spec §3.6). Infallible (yields only `Ok`).
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

    /// Tile a whole contig windowed through `fetcher`, memory-bounded (spec §3.6). The
    /// contig length is taken from `fetcher.length()`. A reference-read failure during the
    /// windowed scan is surfaced here as `Err` (fail-fast), keeping iteration infallible.
    pub fn stream(
        fetcher: impl ChromRefFetcher,
        periods: PeriodRange,
        params: &ScanParams,
        opts: &SegmentOptions,
    ) -> Result<Self, ScanError> {
        let contig_len = u64::from(fetcher.length());
        // Collected straight back up: this consumer was always whole-contig-eager
        // (`build_regions` needs every window's coverage before it can union a
        // satellite). The streaming shape exists for the typed-region walk, which
        // cannot afford that — see `scan_windowed`.
        //
        // This closure is the **canonical** byte source: production's fetcher folds
        // to `{A,C,G,T,N}` as it fills, which is right for a routing tiling and
        // wrong for the walk (see `scan_windowed`'s docs on why the source is a
        // parameter at all).
        //
        // PANIC-FREE: every coordinate `scan_windowed` produces is bounded by
        // `contig_len`, which came from this fetcher's own `u32` length, so both
        // casts fit by construction. `expect` over `as`: were that reasoning ever
        // to lapse, a truncated fetch would silently scan the WRONG WINDOW — the
        // failure this module is least able to see.
        let mut coverage = Vec::new();
        let mut str_intervals = Vec::new();
        let windows = scan_windowed(contig_len, periods, params, opts, |start, len, dst| {
            let start = u32::try_from(start).expect("start <= contig_len, itself a u32");
            let len = u32::try_from(len).expect("len <= contig_len, itself a u32");
            let bytes = fetcher
                .fetch(start, len)
                .map_err(|source| ScanError::Fetch { source })?;
            dst.clear();
            dst.extend_from_slice(&bytes);
            Ok(())
        });
        for window in windows {
            let window = window?;
            // The two window rules, applied to the **raw** detections: this consumer
            // has no policy to insert between detection and merge, which is exactly
            // why it cannot serve the typed-region walk and `scan_windowed` is public
            // (`typed_regions.md` §6.1). The walk applies the same two rules to its
            // pre-filtered set.
            coverage.extend(window.coverage_in_core(&window.detections));
            str_intervals.extend(window.starting_in_core(&window.detections).copied());
        }
        let regions = build_regions(coverage, str_intervals, contig_len, opts);
        Ok(Self {
            regions: regions.into_iter(),
        })
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
        for _ in 0..2000 {
            let len = (next() % 160) as usize;
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
        assert_eq!(period3[0].end, seq.len() as u64);
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
        assert_eq!(period3[0].end, seq.len() as u64);
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
    fn triple(r: &Region) -> (u64, u64, u8) {
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
    fn assert_tiles(regions: &[Region], n: u64) {
        let mut pos = 0u64;
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
        assert_tiles(&regions, seq.len() as u64);
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
        assert_tiles(&regions, seq.len() as u64);
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
        assert_tiles(&regions, seq.len() as u64);
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
        assert_tiles(&off, seq.len() as u64);
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
        assert_tiles(&on, seq.len() as u64);
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
        assert_tiles(&on, seq.len() as u64);
    }

    // ---- scan_windowed (the streamed primitive, B1) -------------------------

    /// `scan_windowed` over a `ChromRefFetcher` — the canonical source, which is
    /// what `RegionScanner` uses and what these tests exercise. The typed-region
    /// walk passes a **raw** source instead; that the two share one copy of the
    /// window rules is the reason the source is a parameter (see `scan_windowed`).
    fn scan_over<'f>(
        fetcher: &'f impl ChromRefFetcher,
        periods: PeriodRange,
        params: &'f ScanParams,
        opts: &'f SegmentOptions,
    ) -> impl std::iter::FusedIterator<Item = Result<ScannedWindow, ScanError>> + 'f {
        scan_windowed(
            u64::from(fetcher.length()),
            periods,
            params,
            opts,
            |start, len, dst| {
                let bytes = fetcher
                    .fetch(start as u32, len as u32)
                    .map_err(|source| ScanError::Fetch { source })?;
                dst.clear();
                dst.extend_from_slice(&bytes);
                Ok(())
            },
        )
    }

    /// The streamed scan must say exactly what the whole-contig-eager form said —
    /// B1 reshapes the seam, it does not change the answer. `RegionScanner::stream`
    /// collecting it back up and still passing its own tests is one half of that;
    /// this is the other, stated directly against the primitive.
    #[test]
    fn scan_windowed_concatenates_to_the_whole_contig_scan() {
        let seq = mixed_contig();
        let (_dir, fetcher) = fetcher_over(&seq);
        let opts = SegmentOptions {
            window_bp: 40,
            ..SegmentOptions::default()
        };

        let windows: Vec<_> = scan_over(&fetcher, p(2, 6), &ScanParams::default(), &opts)
            .map(|w| w.expect("fetch"))
            .collect();
        assert!(windows.len() > 1, "the fixture must actually span windows");

        // The exact intervals, concatenated across windows, are the SAME SET a
        // single whole-contig scan finds — each attributed to exactly one core,
        // none lost, none doubled.
        //
        // **As a set, not a sequence, and the difference is real** (review found
        // this): `find_tandem_repeats` is period-major, so a whole-contig scan
        // groups by period while the stream groups by window. On this fixture at
        // `p(2, 6)` the stream yields periods [2,3,4,6,2,4,4] and the resident
        // scan [2,2,3,4,4,4,6] — same intervals, different order. An ordered
        // assertion here would be true only for a single period, which is exactly
        // how it was first written and why it passed.
        let mut streamed: Vec<RepeatInterval> = windows
            .iter()
            .flat_map(|w| w.starting_in_core(&w.detections).copied())
            .collect();
        let mut resident = find_tandem_repeats(&seq, p(2, 6), &ScanParams::default());
        assert!(!resident.is_empty(), "the fixture must find repeats");
        let key = |i: &RepeatInterval| (i.start, i.end, i.period);
        streamed.sort_by_key(key);
        resident.sort_by_key(key);
        assert_eq!(
            streamed, resident,
            "the streamed intervals must be exactly the resident scan's set"
        );
    }

    /// The cores tile `[0, contig_len)` exactly — contiguous, non-overlapping,
    /// complete. That is what lets a consumer trust that concatenating per-window
    /// coverage reconstructs the contig's, and it is the property the typed-region
    /// walk's own partition invariant will rest on.
    #[test]
    fn scan_windowed_cores_tile_the_contig_exactly() {
        let seq = mixed_contig();
        let (_dir, fetcher) = fetcher_over(&seq);
        for window_bp in [1u64, 7, 40, 1000, 100_000] {
            let opts = SegmentOptions {
                window_bp,
                ..SegmentOptions::default()
            };
            let cores: Vec<RegionSpan> =
                scan_over(&fetcher, p(3, 3), &ScanParams::default(), &opts)
                    .map(|w| w.expect("fetch").core)
                    .collect();

            assert_eq!(cores[0].start, 0, "window_bp = {window_bp}: starts at 0");
            assert_eq!(
                cores.last().unwrap().end,
                seq.len() as u64,
                "window_bp = {window_bp}: ends at the contig end"
            );
            for pair in cores.windows(2) {
                assert_eq!(
                    pair[0].end, pair[1].start,
                    "window_bp = {window_bp}: cores must abut, not gap or overlap"
                );
            }
            for c in &cores {
                assert!(c.start < c.end, "window_bp = {window_bp}: no empty core");
            }
        }
    }

    /// Coverage is clipped to the core and owned intervals are kept whole — the two
    /// derivations that make the scan window-invariant, and the reason they are
    /// **methods over a caller-chosen interval set** rather than fields.
    #[test]
    fn scan_windowed_clips_coverage_to_the_core_but_keeps_intervals_whole() {
        let seq = mixed_contig();
        let (_dir, fetcher) = fetcher_over(&seq);
        // A window small enough that the (CAG)*15 tract must straddle a core edge.
        let opts = SegmentOptions {
            window_bp: 8,
            ..SegmentOptions::default()
        };
        let windows: Vec<_> = scan_over(&fetcher, p(3, 3), &ScanParams::default(), &opts)
            .map(|w| w.expect("fetch"))
            .collect();

        for w in &windows {
            for (s, e) in w.coverage_in_core(&w.detections) {
                assert!(
                    s >= w.core.start && e <= w.core.end,
                    "coverage {s}..{e} escapes core {:?}",
                    w.core
                );
            }
            for iv in w.starting_in_core(&w.detections) {
                assert!(
                    iv.start >= w.core.start && iv.start < w.core.end,
                    "an interval is attributed to the core holding its START"
                );
            }
        }

        // Whole, not clipped: at least one owned interval must reach past its own
        // core, or this fixture is not exercising the property it claims to.
        assert!(
            windows.iter().any(|w| w
                .starting_in_core(&w.detections)
                .any(|iv| iv.end > w.core.end)),
            "the fixture must place a tract across a core edge"
        );
    }

    /// **The margin is visible, and it is what a set-wise consumer stands on.**
    ///
    /// `detections` carries the margin's repeats too, which is the whole reason the
    /// walk can pre-filter a window and reach the verdict a whole-contig run would:
    /// redundancy elimination and the bundle flank test are both *neighbour* rules,
    /// and the neighbours live in the margin. The old shape threw them away.
    ///
    /// Also pins the fetched span, which the walk re-fetches its bases from (it must
    /// not compute that range itself).
    #[test]
    fn a_window_carries_its_margins_detections_and_says_what_it_fetched() {
        let seq = mixed_contig();
        let (_dir, fetcher) = fetcher_over(&seq);
        let opts = SegmentOptions {
            window_bp: 8,
            // A margin the fixture can actually show: the 1 kb default would swallow
            // this whole contig, and every window would be the whole contig.
            max_repeat_len: 12,
            ..SegmentOptions::default()
        };
        let windows: Vec<_> = scan_over(&fetcher, p(3, 3), &ScanParams::default(), &opts)
            .map(|w| w.expect("fetch"))
            .collect();

        for w in &windows {
            assert_eq!(
                w.fetched.start,
                w.core.start.saturating_sub(12),
                "the fetched slice reaches a margin left of the core, clamped at 0"
            );
            assert_eq!(
                w.fetched.end,
                (w.core.end + 12).min(seq.len() as u64),
                "and a margin right of it, clamped at the contig end"
            );
            for iv in &w.detections {
                assert!(
                    iv.start >= w.fetched.start && iv.end <= w.fetched.end,
                    "a detection lies inside the slice it was found in"
                );
            }
        }

        // The property that matters: some window sees a repeat it does NOT own —
        // the neighbour context. Without this, `detections` would be `starting_in_core`
        // by another name and the walk's pre-filter would have no more context than
        // the old shape gave it.
        assert!(
            windows.iter().any(|w| {
                let owned = w.starting_in_core(&w.detections).count();
                w.detections.len() > owned
            }),
            "the margins must actually carry detections the core does not own"
        );
    }

    /// **A fetch failure is reported once and stops the scan.**
    ///
    /// New behaviour in B1, and untested until mutation said so: the eager form
    /// used `?` and returned `Err` immediately, so there was no "and then what?".
    /// An iterator has to answer it, and the wrong answers are both silent — keep
    /// yielding `Ok(empty)` and the consumer scans a **hole** it cannot see; keep
    /// going after the `Err` and it gets partial data past a failure.
    ///
    /// Driven by a `.fai` that overstates the contig: the fetcher believes there
    /// are 4 kb to read and the file holds ~100 bytes, so a window past the real
    /// EOF fails the way a truncated or corrupt reference would.
    #[test]
    fn scan_windowed_reports_a_fetch_failure_once_then_stops() {
        use std::io::Write;
        std::fs::create_dir_all("tmp").unwrap();
        let dir = tempfile::tempdir_in("tmp").unwrap();
        let fa = dir.path().join("ref.fa");
        let fai = dir.path().join("ref.fa.fai");
        let seq = mixed_contig();
        {
            let mut f = std::fs::File::create(&fa).unwrap();
            f.write_all(b">chr\n").unwrap();
            f.write_all(&seq).unwrap();
            f.write_all(b"\n").unwrap();
        }
        {
            // The lie: claim 4 kb where the file holds `seq.len()`.
            let mut f = std::fs::File::create(&fai).unwrap();
            writeln!(f, "chr\t4000\t5\t4000\t4001").unwrap();
        }
        let fetcher = crate::fasta::StreamingChromRefFetcher::for_contig(&fa, "chr").unwrap();
        let opts = SegmentOptions {
            window_bp: 64,
            max_repeat_len: 8,
            ..SegmentOptions::default()
        };

        let items: Vec<_> = scan_over(&fetcher, p(3, 3), &ScanParams::default(), &opts).collect();

        let errs = items.iter().filter(|i| i.is_err()).count();
        assert_eq!(
            errs, 1,
            "the failure is reported exactly once, not per window"
        );
        assert!(
            items.last().unwrap().is_err(),
            "the Err is the LAST item: the scan must stop, not carry on past a \
             window it could not read"
        );
        assert!(
            items[..items.len() - 1].iter().all(|i| i.is_ok()),
            "windows before the failure are still reported"
        );
    }

    /// **`window_bp = 0` must terminate.** The `.max(1)` is the only thing between
    /// this and `core_end == core`, which never advances the cursor: an iterator
    /// yielding empty windows **forever** — no panic, no memory growth, just a
    /// hang. Worse than a crash, and `window_bp` is a `pub` unvalidated field on a
    /// now-`pub` function.
    ///
    /// Untested until review: the other sweep starts at 1, and the pre-existing
    /// invariance helper does its own `.max(1)` over a sweep whose smallest value
    /// is 3 — so zero never reached the function from any direction. Coercion (not
    /// rejection) is the inherited behaviour and is kept: `window_bp` is a pure
    /// memory knob, so 0 → 1 changes the memory profile and not the answer, which
    /// the equality below states.
    #[test]
    fn scan_windowed_terminates_at_window_bp_zero() {
        let seq = mixed_contig();
        let (_dir, fetcher) = fetcher_over(&seq);
        let degenerate = SegmentOptions {
            window_bp: 0,
            ..SegmentOptions::default()
        };
        let windows: Vec<_> = scan_over(&fetcher, p(2, 6), &ScanParams::default(), &degenerate)
            .map(|w| w.expect("fetch"))
            .collect();

        // It terminates (reaching this line at all is the assertion), and coerces
        // to a 1 bp core: one window per base.
        assert_eq!(
            windows.len(),
            seq.len(),
            "window_bp = 0 coerces to 1, so there is one core per base"
        );
        assert!(windows.iter().all(|w| w.core.end == w.core.start + 1));

        // And it is still the same answer — a memory knob must not move the result.
        let mut streamed: Vec<RepeatInterval> = windows
            .iter()
            .flat_map(|w| w.starting_in_core(&w.detections).copied())
            .collect();
        let mut resident = find_tandem_repeats(&seq, p(2, 6), &ScanParams::default());
        let key = |i: &RepeatInterval| (i.start, i.end, i.period);
        streamed.sort_by_key(key);
        resident.sort_by_key(key);
        assert_eq!(
            streamed, resident,
            "window_bp = 0 must not change the answer"
        );
    }

    /// A repeat-free contig still tiles: every window yields, with nothing in it.
    #[test]
    fn scan_windowed_yields_empty_windows_for_repeat_free_sequence() {
        let seq = b"ACGTTGCAAGCTTGCAACGTTGCAAGCTTGCA".repeat(4);
        let (_dir, fetcher) = fetcher_over(&seq);
        let opts = SegmentOptions {
            window_bp: 16,
            ..SegmentOptions::default()
        };
        let windows: Vec<_> = scan_over(&fetcher, p(3, 3), &ScanParams::default(), &opts)
            .map(|w| w.expect("fetch"))
            .collect();
        assert!(windows.len() > 1);
        assert!(
            windows.iter().all(|w| w.detections.is_empty()),
            "no repeats to find, but every window must still be reported"
        );
    }

    // ---- RegionScanner::stream (windowed, C2) ------------------------------

    /// Build a single-contig `StreamingChromRefFetcher` over `seq` via a project-local temp
    /// FASTA + a hand-computed `.fai` (single-line record, so `offset` is the fixed 5-byte
    /// `>chr\n` header). The returned `TempDir` must be kept alive for the file to persist.
    fn fetcher_over(seq: &[u8]) -> (tempfile::TempDir, crate::fasta::StreamingChromRefFetcher) {
        use std::io::Write;
        std::fs::create_dir_all("tmp").unwrap();
        let dir = tempfile::tempdir_in("tmp").unwrap();
        let fa = dir.path().join("ref.fa");
        let fai = dir.path().join("ref.fa.fai");
        {
            let mut f = std::fs::File::create(&fa).unwrap();
            f.write_all(b">chr\n").unwrap();
            f.write_all(seq).unwrap();
            f.write_all(b"\n").unwrap();
        }
        {
            let mut f = std::fs::File::create(&fai).unwrap();
            // name  length  offset  line_bases  line_width
            writeln!(f, "chr\t{}\t5\t{}\t{}", seq.len(), seq.len(), seq.len() + 1).unwrap();
        }
        let fetcher = crate::fasta::StreamingChromRefFetcher::for_contig(&fa, "chr").unwrap();
        (dir, fetcher)
    }

    /// The mixed test contig: two flanks, a (CAG)*15 tract, a non-repetitive gap, an
    /// (AT)*20 tract, and a trailing flank — repeats, unique, and multi-period overlap.
    fn mixed_contig() -> Vec<u8> {
        let mut s = b"TGCATGCA".to_vec();
        for _ in 0..15 {
            s.extend_from_slice(b"CAG");
        }
        s.extend_from_slice(b"TTGACGTACGT"); // 11 bp non-repetitive gap
        for _ in 0..20 {
            s.extend_from_slice(b"AT");
        }
        s.extend_from_slice(b"GGCCTTAACC");
        s
    }

    fn stream_regions(seq: &[u8], periods: PeriodRange, opts: &SegmentOptions) -> Vec<Region> {
        let (_dir, fetcher) = fetcher_over(seq);
        RegionScanner::stream(fetcher, periods, &ScanParams::default(), opts)
            .unwrap()
            .map(|r| r.unwrap())
            .collect()
    }

    /// The streamed tiling equals the resident tiling for every `window_bp` — the
    /// boundary-halo / attribution correctness gate (spec §3.6; the analogue of the cohort
    /// DUST split-invariance test).
    fn assert_window_invariant(seq: &[u8], periods: PeriodRange, base: &SegmentOptions) {
        let oracle = regions_of(seq, periods, base);
        assert_tiles(&oracle, seq.len() as u64);
        for wbp in [seq.len() as u64 + 5, 100, 40, 16, 7, 3] {
            let opts = SegmentOptions {
                window_bp: wbp.max(1),
                ..*base
            };
            let streamed = stream_regions(seq, periods, &opts);
            assert_eq!(
                streamed, oracle,
                "window_bp={wbp}: streamed tiling must equal the resident oracle"
            );
        }
    }

    /// Two isolated STRs (each ≤ `max_repeat_len`) separated by an aperiodic gap wider than
    /// the window margin — so a small `max_repeat_len` genuinely truncates the per-window
    /// fetch, exercising real cross-window coverage assembly (not the whole-contig-per-window
    /// case the default cap degenerates to on a tiny contig).
    fn spaced_str_contig() -> Vec<u8> {
        let mut s = Vec::new();
        let mut st: u64 = 0x1234_5678_9abc_def0;
        let mut filler = |n: usize, out: &mut Vec<u8>| {
            for _ in 0..n {
                st = st
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                out.push(b"ACGT"[((st >> 40) % 4) as usize]);
            }
        };
        filler(15, &mut s); // left flank
        for _ in 0..8 {
            s.extend_from_slice(b"CAG"); // STR 1: 24 bp
        }
        filler(45, &mut s); // aperiodic gap, wider than the 30 bp margin below
        for _ in 0..10 {
            s.extend_from_slice(b"AT"); // STR 2: 20 bp
        }
        filler(15, &mut s); // right flank
        s
    }

    #[test]
    fn stream_is_window_invariant_default() {
        assert_window_invariant(&mixed_contig(), p(2, 6), &SegmentOptions::default());
    }

    #[test]
    fn stream_is_window_invariant_with_smoothing() {
        let opts = SegmentOptions {
            merge_gap: 15, // bridges the 11 bp gap between the two tracts
            min_repeat_len: 5,
            ..SegmentOptions::default()
        };
        assert_window_invariant(&mixed_contig(), p(2, 6), &opts);
    }

    #[test]
    fn stream_is_window_invariant_across_truncated_windows() {
        // A 30 bp cap (hence 30 bp margins) < the contig, so small windows really truncate the
        // fetch. Both STRs are ≤ 30 bp and the gap between them exceeds the margin, so each is
        // captured whole by its own window and segmented identically to the resident scan.
        let opts = SegmentOptions {
            max_repeat_len: 30,
            ..SegmentOptions::default()
        };
        assert_window_invariant(&spaced_str_contig(), p(2, 6), &opts);
    }

    #[test]
    fn stream_detects_a_satellite_within_a_window() {
        // A small cap makes the whole mixed contig one satellite. Windowed streaming
        // reproduces the resident tiling when the window contains the satellite; a satellite
        // longer than a window's reach may segment differently (context-dependent, and it is
        // mask/skip not an analysis target — see `scan_windowed`), so this asserts the
        // fits-in-window case only.
        let seq = mixed_contig();
        let opts = SegmentOptions {
            max_repeat_len: 20,
            ..SegmentOptions::default()
        };
        let oracle = regions_of(&seq, p(2, 6), &opts);
        assert!(
            oracle.iter().any(|r| matches!(r, Region::Satellite(_))),
            "the capped contig is a satellite: {oracle:?}"
        );
        for wbp in [seq.len() as u64 + 5, seq.len() as u64, 150] {
            let o = SegmentOptions {
                window_bp: wbp,
                ..opts
            };
            assert_eq!(stream_regions(&seq, p(2, 6), &o), oracle, "window_bp={wbp}");
        }
    }
}
