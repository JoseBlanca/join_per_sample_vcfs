//! Stage 3 — sdust low-complexity filter.
//!
//! Wraps a `PerPositionMerger` stream and silently drops items whose
//! reference position is low-complexity according to the symmetric
//! DUST algorithm (Morgulis et al. 2006). The filter sits between
//! the multi-way per-sample merge and the variant grouper:
//!
//! ```text
//! .psp_i ─┐
//! .psp_j ─┼─► PerPositionMerger ─► DustFilter ─► VariantGrouper ─► …
//! .psp_k ─┘                            ▲
//!                                      │
//!                          reference FASTA ──┘
//! ```
//!
//! See `doc/devel/implementation_plans/dust_filter.md` and
//! `doc/devel/specs/calling_pipeline_architecture.md` §"Stage 3 —
//! low-complexity filter".
//!
//! ## Algorithm
//!
//! The core is `sdust_mask`, a pure function that ports
//! `sdust_core` from `lh3/sdust` (vendored at `sdust/sdust.c`) line
//! by line. Given a byte slice of ACGT/N reference bases, it returns
//! the sorted, non-overlapping list of low-complexity intervals
//! (BED-style, 0-based half-open). The condition that triggers
//! masking is `10·s > T·L` (strictly greater), where `s` is the
//! triplet-pair score of a candidate subinterval and `L` is its
//! triplet count.
//!
//! Lowercase ACGT is treated identically to uppercase. N (and any
//! other non-ACGT byte) breaks the triplet stream — perfect
//! intervals never span an N.
//!
//! ## Streaming model
//!
//! [`DustFilter`] does not integrate sdust's per-base state machine
//! with the upstream `(chrom_id, pos)` iterator. Instead, on the
//! first emission of a new chromosome it fetches the full reference
//! sequence, runs `sdust_mask` once over it, keeps the resulting
//! masked-interval list together with a monotonic sweep pointer, and
//! answers pass/skip in `O(1)` amortised per upstream position.
//!
//! This trades the spec's stated `O(w)` memory bound (see
//! `doc/devel/specs/calling_pipeline_architecture.md` §"Stage 3")
//! for simplicity and a direct one-to-one port of `sdust_core`. See
//! `doc/devel/implementation_plans/dust_filter.md` §"Streaming
//! model" for the rationale.
//!
//! ## Coordinates
//!
//! [`PileupRecord::pos`](crate::per_sample_pileup::pileup::PileupRecord)
//! is 1-based; the internal sdust mask is 0-based half-open. The
//! conversion happens at the sweep-lookup boundary
//! ([`DustFilter::next`]) and nowhere else. `DustFilter` assumes
//! the upstream emits positions in strict monotonic order per
//! chromosome (as [`PerPositionMerger`](crate::var_calling::per_position_merger::PerPositionMerger)
//! guarantees); behaviour on contract violation is undefined and is
//! not the filter's responsibility to detect.

use std::collections::VecDeque;

use thiserror::Error;

use crate::per_sample_pileup::pileup::RefSeqFetcher;
use crate::per_sample_pileup::psp::header::ParsedChromosome;
use crate::var_calling::per_position_merger::{PerPositionMergerError, PerPositionPileups};

/// Default sdust window size (`W`). Matches the `lh3/sdust` and
/// minimap2 defaults — the maximum subinterval length the algorithm
/// considers.
pub const DEFAULT_DUST_WINDOW: u32 = 64;

/// Default sdust score-density threshold (`T`). Matches the
/// `lh3/sdust` and minimap2 defaults. A subinterval is flagged
/// low-complexity when `10·score > T·triplet_count`, i.e. the
/// score-density per triplet exceeds `T/10 = 2.0` at the default.
pub const DEFAULT_DUST_THRESHOLD: u32 = 20;

/// Triplet length, in bases. Fixed by the algorithm — every
/// reference to "triplet" in this module means a 3-base word.
/// Doubles as the minimum window the validator accepts — a window
/// smaller than the triplet length leaves no room for an interval.
pub const SD_WLEN: u32 = 3;

/// Largest [`DustFilterConfig::window`] the validator accepts.
/// Anchored at the value beyond which the `u32` running scores
/// `SdustState::rw` / `rv` could overflow: the maximum score for a
/// window of `W` triplets is `W·(W−1)/2`, which crosses
/// `u32::MAX` at `W ≈ 92682`. We round down to `u16::MAX` for a
/// clean number well below that boundary; any further use case
/// needing a larger window must first switch the score fields to
/// `u64`.
pub const MAX_DUST_WINDOW: u32 = u16::MAX as u32;

/// Number of distinct ACGT triplets (`4³`). Sizes the per-triplet
/// count arrays.
const SD_WTOT: usize = 1 << (SD_WLEN * 2);

/// Mask isolating the lower `SD_WLEN * 2 = 6` bits of a packed
/// triplet word.
const SD_WMSK: u32 = (SD_WTOT as u32) - 1;

/// ACGT-to-2bit table: `A/a→0, C/c→1, G/g→2, T/t→3`. Every other
/// byte (including `N/n`, ambiguity codes, and whitespace) maps to
/// `4`, which the algorithm treats as a triplet-stream break.
///
/// Ported verbatim from `sdust/sdust.c:22-39`.
const SEQ_NT4: [u8; 256] = build_seq_nt4_table();

const fn build_seq_nt4_table() -> [u8; 256] {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
}

/// Sorted, non-overlapping low-complexity intervals on a single
/// reference slice. Each `(start, finish)` is BED-style: bases at
/// indices `start..finish` (0-based, half-open) are masked.
/// Private — the public API surfaces masking decisions only through
/// [`DustFilter`]'s `Iterator` impl.
type SdustIntervals = Vec<(u32, u32)>;

// ============================================================
// Pure algorithmic core — `sdust_mask` and its helpers
// ============================================================

/// One candidate maximal subinterval inside the current window.
/// Coordinates are *base* coordinates (0-based half-open `finish`);
/// `r` is the triplet-pair score and `l` is the triplet count.
///
/// Ported from `sdust/sdust.c:12-15`.
#[derive(Debug, Clone, Copy)]
struct PerfInterval {
    start: u32,
    finish: u32,
    r: u32,
    l: u32,
}

/// Working state of the sdust scanner, lifted out of `sdust_core`
/// for readability. Field roles mirror the C variables one-to-one.
struct SdustState {
    /// Triplet deque, up to `W − SD_WLEN + 1 = W − 2` entries. Each
    /// entry is a packed 6-bit word with the high bits zeroed.
    window: VecDeque<u8>,
    /// Whole-window per-triplet counts.
    cw: [u32; SD_WTOT],
    /// Active-suffix per-triplet counts.
    cv: [u32; SD_WTOT],
    /// Whole-window score `Σ cw[t]·(cw[t]−1)/2`, kept incrementally.
    rw: u32,
    /// Active-suffix score, kept incrementally.
    rv: u32,
    /// Active-suffix length, in triplets. `L` in the C source.
    big_l: u32,
    /// Candidate maximal subintervals, sorted by descending `start`
    /// then ascending `finish` (the property `find_perfect` relies
    /// on for its dominance check).
    perf: Vec<PerfInterval>,
    /// Emitted masked intervals so far.
    masked_intervals: SdustIntervals,
    /// `W` parameter (window cap, in bases).
    window_cap: u32,
    /// `T` parameter (score-density threshold).
    threshold: u32,
}

/// `true` iff `10·score > threshold · length`, computed with `u64`
/// widening so the comparison stays correct even when `threshold` is
/// large. Mirrors the score-density gate from `sdust_core` (the
/// `cv[t] * 10 > T << 1`, `new_r * 10 > T * new_l`, and
/// `rw * 10 > L * T` comparisons at `sdust/sdust.c:78, 111, 148`).
#[inline]
fn dust_density_exceeds(score: u32, length: u32, threshold: u32) -> bool {
    (score as u64) * 10 > (threshold as u64) * (length as u64)
}

impl SdustState {
    fn new(window_cap: u32, threshold: u32) -> Self {
        // Capacity is the max triplet count for a window of `window_cap`
        // bases: `window_cap − SD_WLEN + 1`. Safe because the public
        // entry point validates `window_cap ≥ SD_WLEN`.
        let triplet_cap = window_cap.saturating_sub(SD_WLEN - 1) as usize;
        Self {
            window: VecDeque::with_capacity(triplet_cap.max(1)),
            cw: [0; SD_WTOT],
            cv: [0; SD_WTOT],
            rw: 0,
            rv: 0,
            big_l: 0,
            perf: Vec::new(),
            masked_intervals: Vec::new(),
            window_cap,
            threshold,
        }
    }

    /// Port of `shift_window` from `sdust/sdust.c:65-85`. Pushes
    /// triplet `t` into the window, evicts the leftmost when full,
    /// and trims the active suffix from the left until `10·cv[t]`
    /// no longer exceeds `2T`.
    fn shift_window(&mut self, t: u8) {
        // Eviction when the deque is at capacity (`W − SD_WLEN + 1`
        // triplets).
        let cap = (self.window_cap - SD_WLEN + 1) as usize;
        if self.window.len() >= cap {
            // PANIC-FREE: `window.len() >= cap` was just tested and
            // `cap >= 1` (because `window_cap >= SD_WLEN = 3` is
            // validated by `DustFilterConfig::new`), so the deque is
            // guaranteed non-empty here.
            let s = self.window.pop_front().expect("deque non-empty at cap") as usize;
            // INVARIANT: `cw[s]` was incremented when triplet `s`
            // entered the window via an earlier `shift_window` push;
            // it has not been decremented since because each triplet
            // is popped at most once on eviction. Same reasoning for
            // the symmetric `cv[s]` decrement below.
            self.cw[s] -= 1;
            self.rw -= self.cw[s];
            if self.big_l as usize > self.window.len() {
                self.big_l -= 1;
                self.cv[s] -= 1;
                self.rv -= self.cv[s];
            }
        }
        self.window.push_back(t);
        self.big_l += 1;
        self.rw += self.cw[t as usize];
        self.cw[t as usize] += 1;
        self.rv += self.cv[t as usize];
        self.cv[t as usize] += 1;
        // Suffix-trim invariant from `sdust.c:78`. With density
        // `length = 2`, the check is exactly the C `cv[t] * 10 > T<<1`.
        if dust_density_exceeds(self.cv[t as usize], 2, self.threshold) {
            loop {
                // INVARIANT: `big_l <= window.len()` (by the
                // suffix-trim contract — the suffix is a suffix of
                // the window), so `window.len() - big_l` does not
                // underflow.
                debug_assert!(self.big_l as usize <= self.window.len());
                let idx = self.window.len() - self.big_l as usize;
                let s = self.window[idx] as usize;
                // INVARIANT: `cv[s]` was last incremented when the
                // triplet at `idx` entered the suffix; the trim loop
                // shrinks the suffix from the left, so each `cv[s]`
                // decrement matches a prior increment.
                self.cv[s] -= 1;
                self.rv -= self.cv[s];
                self.big_l -= 1;
                if s == t as usize {
                    break;
                }
            }
        }
    }

    /// Port of `save_masked_regions` from `sdust/sdust.c:87-101`.
    /// Emits the smallest-start candidate interval (merging with the
    /// previous emitted interval if overlapping or adjacent), then
    /// truncates `perf` to drop every entry whose `start` has fallen
    /// out of the current window.
    fn save_masked_regions(&mut self, start_threshold: u32) {
        // `perf` is sorted by descending `start`, so the *smallest*
        // `start` is at the back — `Vec::last` is the candidate.
        let Some(&candidate) = self.perf.last() else {
            return;
        };
        if candidate.start >= start_threshold {
            return;
        }
        // Emit, merging into the previous emitted interval when
        // overlapping or adjacent.
        if let Some(last) = self.masked_intervals.last_mut()
            && candidate.start <= last.1
        {
            last.1 = last.1.max(candidate.finish);
        } else {
            self.masked_intervals
                .push((candidate.start, candidate.finish));
        }
        // Drop every entry whose `start` is past the window's left
        // edge. Because the list is sorted by descending start,
        // these are exactly the suffix of `perf`.
        while self.perf.last().is_some_and(|p| p.start < start_threshold) {
            self.perf.pop();
        }
    }

    /// Port of `find_perfect` from `sdust/sdust.c:103-127`. Scans
    /// the window backward from the suffix boundary, computing the
    /// triplet-pair score of each candidate ending at the window's
    /// right edge, and inserts new dominant intervals into `perf`.
    fn find_perfect(&mut self, start: u32) {
        let mut c = self.cv;
        let mut r = self.rv;
        let mut max_r: u32 = 0;
        let mut max_l: u32 = 0;
        let win_len = self.window.len() as u32;
        // INVARIANT: `shift_window` maintains `big_l <= win_len`
        // through the suffix-trim loop. A `debug_assert` here makes a
        // future refactor that breaks the invariant panic loudly
        // under tests instead of silently producing wrong masks.
        debug_assert!(
            self.big_l <= win_len,
            "shift_window must keep big_l <= window.len(); big_l={} win_len={}",
            self.big_l,
            win_len,
        );
        // Iterate over the prefix that precedes the active suffix,
        // from rightmost (window.len() - big_l - 1) down to 0. If
        // `big_l` covers the whole window the range is empty and
        // the loop is skipped — same effect as the C source's
        // `(long)kdq_size(w) - L - 1 >= 0` guard.
        let prefix_len = win_len.saturating_sub(self.big_l);
        for i in (0..prefix_len).rev() {
            let t = self.window[i as usize] as usize;
            r += c[t];
            c[t] += 1;
            let new_r = r;
            // INVARIANT: `i < prefix_len <= win_len`, so
            // `win_len - i - 1` does not underflow.
            let new_l = win_len - i - 1;
            if dust_density_exceeds(new_r, new_l, self.threshold) {
                // Find the insertion position `j`: the first entry
                // whose `start < i + start`. Along the way, track
                // the strongest dominator within the candidate's
                // span.
                let candidate_base_start = i + start;
                let mut j = 0usize;
                while j < self.perf.len() && self.perf[j].start >= candidate_base_start {
                    let p = self.perf[j];
                    // C: `p->r * max_l > max_r * p->l` — compare
                    // densities by cross-multiplication.
                    if max_r == 0 || (p.r as u64) * (max_l as u64) > (max_r as u64) * (p.l as u64) {
                        max_r = p.r;
                        max_l = p.l;
                    }
                    j += 1;
                }
                // Insert unless dominated by an existing better
                // candidate. C: `new_r * max_l >= max_r * new_l`.
                if max_r == 0 || (new_r as u64) * (max_l as u64) >= (max_r as u64) * (new_l as u64)
                {
                    max_r = new_r;
                    max_l = new_l;
                    self.perf.insert(
                        j,
                        PerfInterval {
                            start: candidate_base_start,
                            finish: win_len + (SD_WLEN - 1) + start,
                            r: new_r,
                            l: new_l,
                        },
                    );
                }
            }
        }
    }
}

/// Slice-based convenience wrapper around [`sdust_mask_streaming`].
/// Used by the in-module unit tests (which build the contig as a
/// `Vec<u8>` fixture); production drives the streaming form directly
/// via [`DustFilter::ensure_mask_for`]. Behaviour identical to
/// driving the streaming form with `seq.iter().copied().map(Ok)`.
#[cfg(test)]
fn sdust_mask(seq: &[u8], window: u32, threshold: u32) -> SdustIntervals {
    sdust_mask_streaming(
        seq.iter().copied().map(Ok),
        seq.len() as u32,
        window,
        threshold,
    )
    .expect("slice-backed iterator never yields Err")
}

/// Run sdust over a stream of bases and return the sorted,
/// non-overlapping list of low-complexity intervals as BED-style
/// 0-based half-open base coordinates.
///
/// Direct port of `sdust_core` from `sdust/sdust.c:129-159`, made
/// streaming so the caller never materialises the whole contig.
/// The masking condition is `10·score > T·triplet_count` (strictly
/// greater); `window` maps to `W` and `threshold` maps to `T`.
///
/// `bases` yields each contig base as `io::Result<u8>` (ACGT/N,
/// either case accepted — uppercased upstream in the production
/// fetcher). The first `Err` short-circuits the scan and is
/// returned to the caller. `length` is the declared contig length
/// in bases; the function consumes exactly that many bases plus the
/// trailing EOF-sentinel iteration that flushes any candidates
/// still in `perf`.
///
/// Private — every external entry point goes through
/// [`DustFilter::ensure_mask_for`], which passes a validated
/// `(window, threshold)` from [`DustFilterConfig`]. The
/// preconditions below are stated for internal callers only:
///
/// - `window >= SD_WLEN = 3` (otherwise no triplet fits in a
///   window). [`DustFilterConfig::new`] enforces this and is the
///   only path to a valid config; the [`debug_assert!`] below
///   converts an internal contract violation into a loud panic
///   under tests.
/// - In production `length` comes from `ParsedChromosome.length`
///   which is `u32`, so base coordinates fit naturally.
pub(crate) fn sdust_mask_streaming<I>(
    bases: I,
    length: u32,
    window: u32,
    threshold: u32,
) -> std::io::Result<SdustIntervals>
where
    I: IntoIterator<Item = std::io::Result<u8>>,
{
    debug_assert!(
        window >= SD_WLEN,
        "sdust_mask_streaming: window must be >= {SD_WLEN}; got {window}"
    );

    let mut state = SdustState::new(window, threshold);
    let mut l: u32 = 0;
    let mut t: u32 = 0;

    let len = length as usize;
    let mut iter = bases.into_iter();
    // Loop bound mirrors the C `i <= l_seq` — the trailing iteration
    // with `b = 4` flushes any candidates still in `perf`.
    for i in 0..=len {
        let b: u32 = if i < len {
            match iter.next() {
                Some(Ok(byte)) => SEQ_NT4[byte as usize] as u32,
                Some(Err(e)) => return Err(e),
                None => {
                    // Stream ended before `length` bases were yielded.
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::UnexpectedEof,
                        format!(
                            "sdust_mask_streaming: stream ended after {i} bases (declared {length})"
                        ),
                    ));
                }
            }
        } else {
            4
        };
        if b < 4 {
            l += 1;
            t = ((t << 2) | b) & SD_WMSK;
            if l >= SD_WLEN {
                // Base coord of the first base in the current ACGT run.
                // INVARIANT: `l` only grows when an ACGT base at
                // index `i` is consumed, so `l <= i + 1` always
                // holds when this branch runs and the subtraction
                // does not underflow.
                let run_origin = (i as u32 + 1) - l;
                // Window's left edge in base coords. Saturates at the
                // run origin until the run is longer than `window`.
                let start = if l > window {
                    (l - window) + run_origin
                } else {
                    run_origin
                };
                state.save_masked_regions(start);
                state.shift_window(t as u8);
                // Whole-window density gate — only scan for perfect
                // intervals when the window itself is dense enough.
                if dust_density_exceeds(state.rw, state.big_l, threshold) {
                    state.find_perfect(start);
                }
            }
        } else {
            // N or end-of-sequence. The C source uses `l - W + 1` as
            // the starting threshold (one past the in-loop value) and
            // advances it until `perf` is empty.
            let run_origin = (i as u32 + 1) - l;
            let mut start = if l + 1 > window {
                (l + 1 - window) + run_origin
            } else {
                run_origin
            };
            while !state.perf.is_empty() {
                state.save_masked_regions(start);
                start += 1;
            }
            // Reset the run state. The C source does **not** reset
            // `window`, `cw`, `cv`, `rw`, `rv`, `big_l`; stale state
            // is shifted out one base at a time as the next ACGT run
            // pushes new triplets. We match that behaviour exactly so
            // the committed golden vectors hold.
            l = 0;
            t = 0;
        }
    }

    Ok(state.masked_intervals)
}

// ============================================================
// Iterator-adaptor layer — `DustFilter`
// ============================================================

/// Validated DUST filter configuration. Fields are private so the
/// only way to obtain a value is through [`DustFilterConfig::new`]
/// or [`DustFilterConfig::default`], both of which guarantee the
/// window/threshold pair is in range.
///
/// "Make illegal states unrepresentable" — code downstream that
/// holds a `DustFilterConfig` does not have to re-validate it.
/// Mirrors the direction the Stage 6 review took for
/// `PosteriorEngineConfig` (Mi12).
#[derive(Debug, Clone, Copy)]
pub struct DustFilterConfig {
    /// sdust `W` — maximum subinterval length the algorithm
    /// considers. Defaults to [`DEFAULT_DUST_WINDOW`]; validated to
    /// `SD_WLEN..=MAX_DUST_WINDOW` (`3..=65535`) by [`Self::new`].
    window: u32,
    /// sdust `T` — score-density threshold. The masking condition
    /// is `10·s > T·L`; the density bar is `T/10` triplet-pairs per
    /// triplet. Defaults to [`DEFAULT_DUST_THRESHOLD`]. Unbounded; a
    /// very large value effectively disables masking.
    threshold: u32,
}

impl DustFilterConfig {
    /// Build a validated config. `window` must lie in
    /// `SD_WLEN..=MAX_DUST_WINDOW`; `threshold` is unbounded (a
    /// very large value is the effective "do not mask anything"
    /// setting).
    ///
    /// # Errors
    ///
    /// Returns [`DustFilterError::InvalidWindow`] if `window` is
    /// outside the supported range.
    pub fn new(window: u32, threshold: u32) -> Result<Self, DustFilterError> {
        if !(SD_WLEN..=MAX_DUST_WINDOW).contains(&window) {
            return Err(DustFilterError::InvalidWindow { window });
        }
        Ok(Self { window, threshold })
    }

    /// Validated sdust window size (`W`).
    pub fn window(&self) -> u32 {
        self.window
    }

    /// Validated sdust score-density threshold (`T`). The masking
    /// condition is `10·s > T·L`; see the module-level doc.
    pub fn threshold(&self) -> u32 {
        self.threshold
    }
}

impl Default for DustFilterConfig {
    /// `(DEFAULT_DUST_WINDOW, DEFAULT_DUST_THRESHOLD)` are
    /// known-valid, so this is infallible. The struct literal
    /// bypasses [`Self::new`]'s validation; if validation tightens
    /// in the future, defaults that pass today still construct.
    fn default() -> Self {
        Self {
            window: DEFAULT_DUST_WINDOW,
            threshold: DEFAULT_DUST_THRESHOLD,
        }
    }
}

/// Errors the filter can surface to its consumer.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum DustFilterError {
    /// Upstream merger surfaced an error. Boxed for the same
    /// `Result`-size reason as
    /// `PerPositionMergerError::PerSampleReader`. The source's
    /// `Display` is **not** interpolated here; let chain renderers
    /// (`anyhow`, `std::error::Error::source`) surface the cause.
    #[error("upstream merger failed")]
    Upstream(#[source] Box<PerPositionMergerError>),

    /// Reference fetch failed when loading a chromosome's bases.
    /// The cause is exposed via `source()`, not interpolated into
    /// this variant's `Display`.
    #[error("reference fetch failed for chrom {chrom_id}")]
    RefFetch {
        chrom_id: u32,
        #[source]
        source: std::io::Error,
    },

    /// Upstream yielded a `chrom_id` not present in the chromosomes
    /// table the filter was constructed with. Defensive — in
    /// production the table comes from the same merger that emits
    /// the records, so this fires only on pathological mocks.
    #[error("upstream emitted unknown chrom_id {chrom_id}")]
    UnknownChromId { chrom_id: u32 },

    /// Upstream emitted a `PileupRecord` with `pos == 0`, which
    /// breaks the documented 1-based-pos contract. Surfaced as an
    /// error (rather than left to silent `u32` underflow) so a
    /// contract violation is loud, not silent. The filter latches
    /// after this — subsequent `next()` calls return `None`.
    #[error("upstream emitted invalid 1-based pos 0 on chrom_id {chrom_id}")]
    InvalidPos { chrom_id: u32 },

    /// `DustFilterConfig::window` is outside the supported range.
    #[error(
        "invalid sdust window {window}: must be in {}..={MAX_DUST_WINDOW}",
        SD_WLEN
    )]
    InvalidWindow { window: u32 },
}

/// Cached mask state for a single chromosome. Bundling these three
/// fields together makes their "all consistent together or all
/// absent" invariant a type-level property: there is no way to have
/// a `chrom_id` set without also having the matching mask and a
/// reset sweep pointer.
struct LoadedChrom {
    chrom_id: u32,
    mask: SdustIntervals,
    sweep: usize,
}

/// Streaming low-complexity filter over a per-position pileup stream.
///
/// Wraps an upstream iterator yielding
/// [`PerPositionPileups`] / [`PerPositionMergerError`] and silently
/// drops items whose reference base is inside an sdust-masked
/// interval. Emits the surviving items in the same order as the
/// upstream.
///
/// On the first item for each chromosome, the filter fetches the
/// whole chromosome's reference via the supplied [`RefSeqFetcher`]
/// and runs `sdust_mask` over it once. Subsequent items on the
/// same chromosome are answered by a sweep over the cached
/// interval list (`O(1)` amortised per upstream position).
pub struct DustFilter<I, F>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    F: RefSeqFetcher,
{
    upstream: I,
    fetcher: F,
    chromosomes: Vec<ParsedChromosome>,
    config: DustFilterConfig,
    /// Currently-loaded chromosome's mask and sweep pointer.
    /// `None` before the first item, and reassigned atomically on
    /// every chromosome change.
    loaded: Option<LoadedChrom>,
    /// Latches `true` after the upstream signals end or any error,
    /// so subsequent `next()` calls return `None`. Matches the
    /// upstream merger's contract.
    is_finished: bool,
}

impl<I, F> DustFilter<I, F>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    F: RefSeqFetcher,
{
    /// Wrap `upstream`. `chromosomes` is the same table the upstream
    /// merger holds — pass `merger.chromosomes().to_vec()`. Infallible
    /// because `config` has already been validated by
    /// [`DustFilterConfig::new`].
    pub fn new(
        upstream: I,
        fetcher: F,
        chromosomes: Vec<ParsedChromosome>,
        config: DustFilterConfig,
    ) -> Self {
        Self {
            upstream,
            fetcher,
            chromosomes,
            config,
            loaded: None,
            is_finished: false,
        }
    }

    /// Returns the validated configuration in effect.
    pub fn config(&self) -> DustFilterConfig {
        self.config
    }

    /// Load (or reload) the mask for `chrom_id` if it isn't already
    /// the current one. Latches `is_finished` and returns the
    /// appropriate error on any failure. Streams the contig's
    /// bases through `sdust_mask_streaming` — no whole-contig
    /// buffer is materialised when the fetcher implements
    /// `iter_bases` natively (e.g. `StreamingChromRefFetcher`).
    fn ensure_mask_for(&mut self, chrom_id: u32) -> Result<(), DustFilterError> {
        if self.loaded.as_ref().is_some_and(|l| l.chrom_id == chrom_id) {
            return Ok(());
        }
        let entry = match self.chromosomes.get(chrom_id as usize) {
            Some(e) => e,
            None => {
                self.is_finished = true;
                return Err(DustFilterError::UnknownChromId { chrom_id });
            }
        };
        let bases = match self.fetcher.iter_bases(chrom_id, entry.length) {
            Ok(it) => it,
            Err(source) => {
                self.is_finished = true;
                return Err(DustFilterError::RefFetch { chrom_id, source });
            }
        };
        let mask = match sdust_mask_streaming(
            bases,
            entry.length,
            self.config.window,
            self.config.threshold,
        ) {
            Ok(m) => m,
            Err(source) => {
                self.is_finished = true;
                return Err(DustFilterError::RefFetch { chrom_id, source });
            }
        };
        // Atomic install: the three previously-co-dependent fields
        // (chrom_id, mask, sweep) move into the new `LoadedChrom` in
        // a single assignment, so they can never be out of sync.
        self.loaded = Some(LoadedChrom {
            chrom_id,
            mask,
            sweep: 0,
        });
        Ok(())
    }

    /// Returns `true` if the 0-based reference position `pos0` is
    /// inside an sdust-masked interval on the currently-loaded
    /// chromosome. Advances the cached sweep pointer past intervals
    /// whose `finish` is `≤ pos0` (they cannot match any future,
    /// larger position).
    ///
    /// Precondition: `ensure_mask_for` succeeded immediately before
    /// this call, so `self.loaded` is `Some(_)`. The
    /// [`debug_assert!`] makes a contract violation panic loudly
    /// under tests; in release, the `None` branch passes the record
    /// through unmasked rather than producing wrong output.
    fn is_masked(&mut self, pos0: u32) -> bool {
        let Some(loaded) = self.loaded.as_mut() else {
            debug_assert!(false, "is_masked called without a loaded chromosome");
            return false;
        };
        while loaded.sweep < loaded.mask.len() && loaded.mask[loaded.sweep].1 <= pos0 {
            loaded.sweep += 1;
        }
        match loaded.mask.get(loaded.sweep) {
            Some(&(start, finish)) => start <= pos0 && pos0 < finish,
            None => false,
        }
    }
}

impl<I, F> Iterator for DustFilter<I, F>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    F: RefSeqFetcher,
{
    type Item = Result<PerPositionPileups, DustFilterError>;

    /// Yields the next upstream item that is not low-complexity.
    ///
    /// Assumes the upstream iterator emits positions in strict
    /// monotonic order per chromosome, as
    /// [`PerPositionMerger`](crate::var_calling::per_position_merger::PerPositionMerger)
    /// guarantees. Out-of-order input would silently produce wrong
    /// pass/skip decisions; detecting that contract violation is
    /// the merger's responsibility, not this filter's.
    fn next(&mut self) -> Option<Self::Item> {
        if self.is_finished {
            return None;
        }
        loop {
            let pileups = match self.upstream.next() {
                None => {
                    self.is_finished = true;
                    return None;
                }
                Some(Err(err)) => {
                    self.is_finished = true;
                    return Some(Err(DustFilterError::Upstream(Box::new(err))));
                }
                Some(Ok(p)) => p,
            };
            if let Err(err) = self.ensure_mask_for(pileups.chrom_id) {
                return Some(Err(err));
            }
            // `PileupRecord::pos` is 1-based; the internal sdust
            // mask is 0-based half-open. A `pos == 0` is a contract
            // violation — surface it as a typed error rather than
            // letting the subtraction silently wrap in release.
            let pos0 = match pileups.pos.checked_sub(1) {
                Some(p) => p,
                None => {
                    self.is_finished = true;
                    return Some(Err(DustFilterError::InvalidPos {
                        chrom_id: pileups.chrom_id,
                    }));
                }
            };
            if self.is_masked(pos0) {
                // Drop and continue to the next upstream item.
                continue;
            }
            return Some(Ok(pileups));
        }
    }
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: build a random-ish high-complexity ACGT string of
    // length `n`. Uses a fixed seed so tests are deterministic and
    // do not depend on `rand`.
    fn high_complexity_bases(n: usize) -> Vec<u8> {
        bases_from_seed(0xdead_beef, n, b"ACGT")
    }

    /// Deterministic xorshift-based byte generator used by the
    /// invariant proptest below in addition to `high_complexity_bases`.
    /// Output depends only on `seed` and `n`, so failures reproduce
    /// across hosts.
    fn bases_from_seed(seed: u64, n: usize, alphabet: &[u8]) -> Vec<u8> {
        let mut state: u64 = seed;
        let len = alphabet.len();
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                alphabet[(state as usize) % len]
            })
            .collect()
    }

    // ---- sdust_mask: algorithmic core ----

    #[test]
    fn sdust_empty_input() {
        assert_eq!(sdust_mask(b"", 64, 20), Vec::<(u32, u32)>::new());
    }

    #[test]
    fn sdust_high_complexity_passes() {
        let seq = high_complexity_bases(256);
        let intervals = sdust_mask(&seq, 64, 20);
        // A random ACGT string at default threshold should produce
        // no masked intervals. If this flakes on a different seed,
        // tighten the seed rather than the threshold.
        assert!(
            intervals.is_empty(),
            "expected no intervals on high-complexity input, got {:?}",
            intervals
        );
    }

    #[test]
    fn sdust_homopolymer_is_masked() {
        // 128-bp pure poly-A flanked by nothing. The whole stretch
        // should land in one interval.
        let seq = vec![b'A'; 128];
        let intervals = sdust_mask(&seq, 64, 20);
        assert_eq!(intervals.len(), 1, "got {:?}", intervals);
        let (s, f) = intervals[0];
        // The masked region should at least cover the dense interior.
        // Exact bounds are pinned by the golden-vector test against
        // the cloned binary.
        assert_eq!(s, 0);
        assert_eq!(f, 128);
    }

    #[test]
    fn sdust_dinucleotide_repeat_is_masked() {
        let mut seq = Vec::with_capacity(128);
        for _ in 0..64 {
            seq.extend_from_slice(b"AT");
        }
        let intervals = sdust_mask(&seq, 64, 20);
        assert!(!intervals.is_empty(), "AT repeat should mask, got empty");
        // A perfect dinucleotide repeat across the whole input should
        // produce one interval spanning everything dense; exact
        // endpoints are pinned by the golden-vector test.
        assert_eq!(intervals[0].0, 0);
    }

    #[test]
    fn sdust_trinucleotide_repeat_is_masked() {
        let mut seq = Vec::with_capacity(126);
        for _ in 0..42 {
            seq.extend_from_slice(b"ATG");
        }
        let intervals = sdust_mask(&seq, 64, 20);
        assert!(!intervals.is_empty(), "ATG repeat should mask, got empty");
    }

    #[test]
    fn sdust_lowercase_treated_as_uppercase() {
        let upper = vec![b'A'; 128];
        let mut lower = upper.clone();
        for b in &mut lower {
            *b = b.to_ascii_lowercase();
        }
        assert_eq!(sdust_mask(&upper, 64, 20), sdust_mask(&lower, 64, 20));
    }

    #[test]
    fn sdust_n_breaks_the_run() {
        // Two homopolymer runs separated by an N. Each should be
        // independently masked; perfect intervals never span the N.
        let mut seq = vec![b'A'; 80];
        seq.push(b'N');
        seq.extend(std::iter::repeat_n(b'A', 80));
        let intervals = sdust_mask(&seq, 64, 20);
        // At minimum: two distinct masked spans, neither covering
        // the N at index 80.
        assert!(intervals.len() >= 2, "got {:?}", intervals);
        for (s, f) in &intervals {
            assert!(
                *f <= 80 || *s >= 81,
                "interval {:?} spans the N at 80",
                (s, f)
            );
        }
    }

    #[test]
    fn sdust_other_non_acgt_also_breaks_run() {
        // Any non-ACGT byte should behave like N. Test with a
        // whitespace byte — `SEQ_NT4` maps everything to 4 except
        // ACGT.
        let mut seq = vec![b'A'; 80];
        seq.push(b' ');
        seq.extend(std::iter::repeat_n(b'A', 80));
        let intervals = sdust_mask(&seq, 64, 20);
        for (s, f) in &intervals {
            assert!(*f <= 80 || *s >= 81);
        }
    }

    #[test]
    fn sdust_high_threshold_disables_masking() {
        let seq = vec![b'A'; 128];
        // T = 5000 → density bar of 500 triplet-pairs per triplet,
        // which a 64-base window cannot reach.
        assert_eq!(sdust_mask(&seq, 64, 5000), Vec::<(u32, u32)>::new());
    }

    #[test]
    fn sdust_intervals_are_sorted_and_disjoint() {
        // Two distant homopolymer runs in a single sequence should
        // produce intervals in increasing-start order with no overlap.
        let mut seq = high_complexity_bases(64);
        seq.extend(std::iter::repeat_n(b'A', 80));
        seq.extend(high_complexity_bases(64));
        seq.extend(std::iter::repeat_n(b'T', 80));
        seq.extend(high_complexity_bases(64));
        let intervals = sdust_mask(&seq, 64, 20);
        for win in intervals.windows(2) {
            assert!(win[0].1 <= win[1].0, "overlap: {:?}", win);
        }
    }

    #[test]
    #[should_panic(expected = "window must be")]
    fn sdust_mask_debug_asserts_on_tiny_window() {
        // The debug_assert only fires in debug; tests always run in
        // debug mode, so this catches any future caller bypassing
        // DustFilterConfig::new and passing window < SD_WLEN.
        let _ = sdust_mask(b"AAAA", 2, 20);
    }

    // ---- Golden-vector tests ----
    //
    // The `expected` intervals were generated once at implementation
    // time (2026-05-17) by running `lh3/sdust -w 64 -t 20` at
    // upstream commit 89c42cb41ba598e9cfa07c2ef99ae8c08f769b3e
    // against each snippet's bases and copying its BED-style
    // `(start, finish)` output verbatim. From that point on the
    // values committed below are the source of truth: the test
    // asserts that `sdust_mask` continues to produce them.
    //
    // The lh3/sdust binary itself is not a build- or test-time
    // dependency of this project. If you change a snippet's bases
    // (or add a new one), regenerating the matching `expected` is a
    // one-off, off-repo task — set up `lh3/sdust` independently,
    // run it on the snippet, paste the output. Code review is the
    // safety net for new-snippet additions.

    struct GoldenSnippet {
        name: &'static str,
        bases: &'static [u8],
        expected: &'static [(u32, u32)],
    }

    const GOLDEN_SNIPPETS: &[GoldenSnippet] = &[
        GoldenSnippet {
            name: "homopolymer_with_flanks",
            bases: b"ACGTCAGTACGATCAGTAGCATGCAGTAGCATCAGTACGAGCATCAGCAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGAGCTAGCATGCAGTGACAGTCACGCATAGCAGTCACGTAGCATCAGAC",
            expected: &[(50, 100)],
        },
        GoldenSnippet {
            name: "dinucleotide_repeat_with_flanks",
            bases: b"ACGTCAGTACGATCAGTAGCATGCAGTAGCATCAGTACGAGCATCAGCAGATATATATATATATATATATATATATATATATATATATATATATATATATATATATTGAGCTAGCATGCAGTGACAGTCACGCATAGCAGTCACGTAGCATCAGAC",
            expected: &[(50, 106)],
        },
        GoldenSnippet {
            name: "trinucleotide_repeat_with_flanks",
            bases: b"ACGTCAGTACGATCAGTAGCATGCAGTAGCATCAGTACGAGCATCAGCAGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGTGAGCTAGCATGCAGTGACAGTCACGCATAGCAGTCACGTAGCATCAGAC",
            expected: &[(49, 110)],
        },
        GoldenSnippet {
            name: "mixed_boundary",
            bases: b"ACGTCAGTACGATCAGTAGCATGCAGTAGCATCAGTACGAGCATCAGCAGAAAAAAAAAAAAAAAAAAAAATGAGCTAGCATGCAGTGACAGTCACGCATAGCAGTCACGTAGCATCAGAC",
            expected: &[(50, 71)],
        },
        GoldenSnippet {
            name: "n_inside_low_complexity",
            bases: b"ACGTCAGTACGATCAGTAGCATGCAGTAGCATCAGTACGAGCATCAGCAGAAAAAAAAAAAAAAAAAANNNNNNNNAAAAAAAAAAAAAAAAAATGAGCTAGCATGCAGTGACAGTCACGCATAGCAGTCACGTAGCATCAGAC",
            expected: &[(50, 68), (106, 140)],
        },
        GoldenSnippet {
            name: "low_complexity_long_tract",
            bases: b"ACGTCAGTACGATCAGTAGCATGCAGTAGCATCAGTACGAGCATCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTGAGCTAGCATGCAGTGACAGTCACGCATAGCAGTCACGTAGCATCAGAC",
            expected: &[(3, 182)],
        },
    ];

    #[test]
    fn sdust_matches_vendored_binary_on_golden_snippets() {
        for snippet in GOLDEN_SNIPPETS {
            let got = sdust_mask(snippet.bases, 64, 20);
            let want: Vec<(u32, u32)> = snippet.expected.to_vec();
            assert_eq!(
                got, want,
                "snippet {:?}: port output disagrees with lh3/sdust",
                snippet.name
            );
        }
    }

    // ---- DustFilterConfig validation ----

    #[test]
    fn config_new_accepts_defaults() {
        let cfg = DustFilterConfig::new(64, 20).unwrap();
        assert_eq!(cfg.window(), 64);
        assert_eq!(cfg.threshold(), 20);
    }

    #[test]
    fn config_new_accepts_boundaries() {
        assert!(DustFilterConfig::new(SD_WLEN, 20).is_ok());
        assert!(DustFilterConfig::new(MAX_DUST_WINDOW, 20).is_ok());
    }

    #[test]
    fn config_new_rejects_too_small_window() {
        match DustFilterConfig::new(SD_WLEN - 1, 20) {
            Err(DustFilterError::InvalidWindow { window }) if window == SD_WLEN - 1 => {}
            other => panic!("expected InvalidWindow below SD_WLEN, got {:?}", other),
        }
    }

    #[test]
    fn config_new_rejects_too_large_window() {
        let too_large = MAX_DUST_WINDOW + 1;
        match DustFilterConfig::new(too_large, 20) {
            Err(DustFilterError::InvalidWindow { window }) if window == too_large => {}
            other => panic!("expected InvalidWindow above MAX, got {:?}", other),
        }
    }

    #[test]
    fn config_new_accepts_threshold_zero_and_max() {
        // The threshold knob intentionally takes any u32 — `0`
        // masks aggressively, `u32::MAX` effectively disables.
        assert!(DustFilterConfig::new(64, 0).is_ok());
        assert!(DustFilterConfig::new(64, u32::MAX).is_ok());
    }

    #[test]
    fn config_default_passes_validation() {
        let cfg = DustFilterConfig::default();
        assert_eq!(cfg.window(), DEFAULT_DUST_WINDOW);
        assert_eq!(cfg.threshold(), DEFAULT_DUST_THRESHOLD);
    }

    // ---- DustFilter: iterator plumbing ----

    use crate::per_sample_pileup::pileup::{AlleleObservation, AlleleSupportStats, PileupRecord};
    use std::cell::RefCell;
    use std::io;

    /// Stub fetcher backed by per-chromosome byte vectors. Records
    /// the number of `fetch` calls per chromosome so tests can pin
    /// that the filter loads each chromosome at most once.
    struct StubFetcher {
        seqs: Vec<Vec<u8>>,
        fail_for: Option<u32>,
        call_counts: RefCell<Vec<u32>>,
    }

    impl StubFetcher {
        fn new(seqs: Vec<Vec<u8>>) -> Self {
            let n = seqs.len();
            Self {
                seqs,
                fail_for: None,
                call_counts: RefCell::new(vec![0; n]),
            }
        }

        fn with_failure_for(mut self, chrom_id: u32) -> Self {
            self.fail_for = Some(chrom_id);
            self
        }

        fn calls_for(&self, chrom_id: u32) -> u32 {
            self.call_counts.borrow()[chrom_id as usize]
        }
    }

    impl RefSeqFetcher for StubFetcher {
        fn fetch(
            &self,
            chrom_id: u32,
            _start_1based: u32,
            length: u32,
        ) -> Result<Vec<u8>, io::Error> {
            self.call_counts.borrow_mut()[chrom_id as usize] += 1;
            if self.fail_for == Some(chrom_id) {
                return Err(io::Error::other("stub fetch failure"));
            }
            let seq = &self.seqs[chrom_id as usize];
            assert_eq!(
                seq.len() as u32,
                length,
                "stub: caller asked for {} bases, have {}",
                length,
                seq.len()
            );
            Ok(seq.clone())
        }
    }

    fn make_chrom(name: &str, length: u32) -> ParsedChromosome {
        ParsedChromosome {
            name: name.to_string(),
            length,
            md5: "0".repeat(32),
        }
    }

    fn empty_pileup(chrom_id: u32, pos: u32) -> PerPositionPileups {
        let ref_allele =
            AlleleObservation::new(b"A".to_vec(), AlleleSupportStats::default(), Vec::new());
        PerPositionPileups {
            chrom_id,
            pos,
            per_sample: vec![Some(PileupRecord::new(chrom_id, pos, vec![ref_allele]))],
        }
    }

    fn positions_iter(
        items: Vec<(u32, u32)>,
    ) -> impl Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>> {
        items.into_iter().map(|(c, p)| Ok(empty_pileup(c, p)))
    }

    #[test]
    fn filter_passes_through_high_complexity_bases() {
        let seq = high_complexity_bases(256);
        let chromosomes = vec![make_chrom("chr1", seq.len() as u32)];
        let fetcher = StubFetcher::new(vec![seq]);
        let upstream = positions_iter((1..=10).map(|p| (0, p)).collect());
        let cfg = DustFilterConfig::default();
        let filter = DustFilter::new(upstream, fetcher, chromosomes, cfg);
        let out: Vec<u32> = filter.map(|r| r.unwrap().pos).collect();
        assert_eq!(out, (1..=10).collect::<Vec<_>>());
    }

    #[test]
    fn filter_drops_homopolymer_positions() {
        // Whole chromosome is 128-bp poly-A. Every upstream position
        // should be dropped.
        let seq = vec![b'A'; 128];
        let chromosomes = vec![make_chrom("chr1", seq.len() as u32)];
        let fetcher = StubFetcher::new(vec![seq]);
        let upstream = positions_iter((1..=20).map(|p| (0, p)).collect());
        let cfg = DustFilterConfig::default();
        let filter = DustFilter::new(upstream, fetcher, chromosomes, cfg);
        let out: Vec<_> = filter.collect();
        assert!(
            out.is_empty(),
            "expected all positions dropped, got {} survivors",
            out.len()
        );
    }

    #[test]
    fn filter_fetcher_call_count_pinned() {
        let chr1 = high_complexity_bases(128);
        let chr2 = high_complexity_bases(128);
        let chromosomes = vec![make_chrom("chr1", 128), make_chrom("chr2", 128)];
        // We need to inspect the fetcher after the filter exhausts.
        // Approach: wrap StubFetcher in &-borrow via the
        // `impl RefSeqFetcher for &T` blanket impl, so the original
        // remains accessible.
        let fetcher = StubFetcher::new(vec![chr1, chr2]);
        let positions: Vec<_> = (1..=50)
            .map(|p| (0, p))
            .chain((1..=50).map(|p| (1, p)))
            .collect();
        let cfg = DustFilterConfig::default();
        let filter = DustFilter::new(positions_iter(positions), &fetcher, chromosomes, cfg);
        let _drain: Vec<_> = filter.collect();
        assert_eq!(fetcher.calls_for(0), 1);
        assert_eq!(fetcher.calls_for(1), 1);
    }

    #[test]
    fn filter_resets_mask_across_chromosomes() {
        // chr0 is a homopolymer (drops everything); chr1 is high
        // complexity (passes everything). After switching to chr1,
        // a position with the same pos as a dropped chr0 position
        // must pass.
        let chr0 = vec![b'A'; 128];
        let chr1 = high_complexity_bases(128);
        let chromosomes = vec![make_chrom("chr1", 128), make_chrom("chr2", 128)];
        let fetcher = StubFetcher::new(vec![chr0, chr1]);
        let upstream = positions_iter(vec![(0, 5), (0, 6), (1, 5), (1, 6)]);
        let cfg = DustFilterConfig::default();
        let filter = DustFilter::new(upstream, fetcher, chromosomes, cfg);
        let out: Vec<_> = filter.map(|r| r.unwrap()).collect();
        assert_eq!(out.len(), 2);
        assert_eq!((out[0].chrom_id, out[0].pos), (1, 5));
        assert_eq!((out[1].chrom_id, out[1].pos), (1, 6));
    }

    #[test]
    fn filter_surfaces_upstream_error_and_latches() {
        let chromosomes = vec![make_chrom("chr1", 128)];
        let fetcher = StubFetcher::new(vec![high_complexity_bases(128)]);
        // First item is Ok, second is an OutOfOrder error. After
        // the error, subsequent next() must be None.
        let items: Vec<Result<PerPositionPileups, PerPositionMergerError>> = vec![
            Ok(empty_pileup(0, 5)),
            Err(PerPositionMergerError::OutOfOrder {
                sample_idx: 0,
                sample_name: "s0".into(),
                regressing_chrom_id: 0,
                regressing_pos: 3,
                last_emitted_chrom_id: 0,
                last_emitted_pos: 5,
            }),
            Ok(empty_pileup(0, 6)),
        ];
        let cfg = DustFilterConfig::default();
        let mut filter = DustFilter::new(items.into_iter(), fetcher, chromosomes, cfg);
        assert!(matches!(filter.next(), Some(Ok(_))));
        assert!(matches!(
            filter.next(),
            Some(Err(DustFilterError::Upstream(_)))
        ));
        assert!(filter.next().is_none(), "iterator must latch after error");
        assert!(filter.next().is_none(), "still latched on subsequent call");
    }

    #[test]
    fn filter_surfaces_ref_fetch_error_and_latches() {
        let chromosomes = vec![make_chrom("chr1", 128)];
        let fetcher = StubFetcher::new(vec![vec![b'A'; 128]]).with_failure_for(0);
        let upstream = positions_iter(vec![(0, 5), (0, 6)]);
        let cfg = DustFilterConfig::default();
        let mut filter = DustFilter::new(upstream, fetcher, chromosomes, cfg);
        match filter.next() {
            Some(Err(DustFilterError::RefFetch {
                chrom_id: 0,
                source: _,
            })) => {}
            other => panic!("expected RefFetch error, got {:?}", other),
        }
        assert!(filter.next().is_none());
    }

    #[test]
    fn filter_surfaces_unknown_chrom_id_and_latches() {
        // Upstream emits chrom_id=5, but the chromosomes table only
        // has chrom_id=0.
        let chromosomes = vec![make_chrom("chr1", 128)];
        let fetcher = StubFetcher::new(vec![high_complexity_bases(128)]);
        let upstream = positions_iter(vec![(5, 1)]);
        let cfg = DustFilterConfig::default();
        let mut filter = DustFilter::new(upstream, fetcher, chromosomes, cfg);
        match filter.next() {
            Some(Err(DustFilterError::UnknownChromId { chrom_id: 5 })) => {}
            other => panic!("expected UnknownChromId, got {:?}", other),
        }
        assert!(filter.next().is_none());
    }

    #[test]
    fn filter_empty_upstream_does_no_work() {
        let chromosomes = vec![make_chrom("chr1", 128)];
        let fetcher = StubFetcher::new(vec![high_complexity_bases(128)]);
        let upstream = positions_iter(Vec::<(u32, u32)>::new());
        let cfg = DustFilterConfig::default();
        let filter = DustFilter::new(upstream, &fetcher, chromosomes, cfg);
        let out: Vec<_> = filter.collect();
        assert!(out.is_empty());
        // The fetcher was never asked for anything.
        assert_eq!(fetcher.calls_for(0), 0);
    }

    // ---- Boundary and contract-violation tests added on review ----

    #[test]
    fn sdust_threshold_is_strictly_greater_at_density_boundary() {
        // Six identical bases form four triplets — the largest
        // candidate the algorithm considers has density exactly
        // `T/10 = 2.0`. With the strict `>` comparison this must NOT
        // mask; a future refactor that relaxed the comparison to `>=`
        // would mask here and fail this assertion. Seven identical
        // bases push past the boundary (density 2.5 over the
        // 5-triplet candidate) and must mask. Together these pin the
        // comparison direction in a way the existing
        // `sdust_high_threshold_disables_masking` smoke does not.
        assert_eq!(
            sdust_mask(b"AAAAAA", 64, 20),
            Vec::<(u32, u32)>::new(),
            "6 same bases sit at density == T/10; strict > must not mask",
        );
        assert_eq!(
            sdust_mask(b"AAAAAAA", 64, 20),
            vec![(0, 7)],
            "7 same bases cross the density boundary and must mask",
        );
    }

    #[test]
    fn sdust_returns_empty_on_seq_shorter_than_wlen() {
        // Sequences too short to fit one triplet must produce nothing
        // and must not panic. Covers the boundary where the inner
        // `if l >= SD_WLEN` branch never fires.
        assert!(sdust_mask(b"A", 64, 20).is_empty());
        assert!(sdust_mask(b"AC", 64, 20).is_empty());
    }

    #[test]
    fn sdust_handles_seq_of_exactly_wlen_bases() {
        // Exactly one triplet: the inner block runs once. Structural
        // check only — no specific output is pinned.
        for &(s, f) in &sdust_mask(b"AAA", 64, 20) {
            assert!(f <= 3, "interval beyond seq.len()={}: {}..{}", 3, s, f);
        }
    }

    #[test]
    fn sdust_returns_empty_on_all_n() {
        // 256 N's. The triplet stream never accumulates ACGT bases,
        // so no PerfInterval is ever inserted and the loop's N
        // branch fires every iteration. Must not panic.
        assert!(sdust_mask(&vec![b'N'; 256], 64, 20).is_empty());
    }

    #[test]
    fn sdust_handles_long_n_run_after_low_complexity() {
        // Low-complexity prefix then a long N gap exercises the
        // N-flush path with a non-empty `perf`, then iterates the N
        // branch for the rest of the sequence without producing more
        // intervals.
        let mut seq = vec![b'A'; 64];
        seq.extend(std::iter::repeat_n(b'N', 256));
        for &(_s, f) in &sdust_mask(&seq, 64, 20) {
            assert!(
                f <= 64,
                "interval should not extend into the N gap: finish={}",
                f
            );
        }
    }

    #[test]
    fn sdust_at_minimum_window_does_not_panic() {
        // `window = SD_WLEN = 3` exercises the boundary of
        // `triplet_cap = window − SD_WLEN + 1 = 1` at the
        // `SdustState::new` construction site. Structural assertions
        // only.
        let out = sdust_mask(&[b'A'; 32], 3, 20);
        for &(s, f) in &out {
            assert!(s < f);
            assert!(f <= 32);
        }
    }

    #[test]
    fn sdust_threshold_zero_masks_aggressively() {
        // T = 0 makes `10·s > 0` fire for any positive score. Even
        // random ACGT input should mask at least one interval.
        let out = sdust_mask(&high_complexity_bases(256), 64, 0);
        assert!(
            !out.is_empty(),
            "T=0 must mask at least one interval on any non-trivial input",
        );
    }

    #[test]
    fn sdust_threshold_max_does_not_overflow() {
        // u32::MAX exercises the u64 widening at every density gate.
        // Asserts the call returns rather than overflowing; output
        // shape is not pinned.
        let _ = sdust_mask(&[b'A'; 128], 64, u32::MAX);
    }

    #[test]
    fn sdust_invariants_hold_on_random_seeded_inputs() {
        // Stand-in for proptest (the project does not yet depend on
        // `proptest`): exercise `sdust_mask` over a fixed sweep of
        // seeded random inputs and assert the structural invariants
        // the iterator layer relies on. Catches any future refactor
        // of `shift_window` / `save_masked_regions` / `find_perfect`
        // that violates "sorted, non-overlapping, in-bounds, never
        // covers a non-ACGT byte" on inputs outside the golden set.
        let alphabets: &[&[u8]] = &[b"ACGT", b"ACGTN", b"acgtACGT", b"AN", b"ACGTNacgtn "];
        for seed in 0..32u64 {
            for &alphabet in alphabets {
                for len in [0usize, 1, 2, 3, 7, 31, 64, 128, 257, 511] {
                    for window in [3u32, 7, 32, 64, 128] {
                        for threshold in [1u32, 5, 20, 50, 200] {
                            let seq = bases_from_seed(seed.wrapping_mul(2654435761), len, alphabet);
                            let out = sdust_mask(&seq, window, threshold);
                            // Structural invariants the iterator layer
                            // relies on: each interval is non-empty
                            // and the list is sorted with no overlaps.
                            //
                            // Two things the assertions intentionally
                            // do *not* check:
                            //
                            // 1. `finish <= seq.len()`. The reported
                            //    `finish` is the algorithm's window
                            //    right edge, which can outrun
                            //    `seq.len()` after an N break (the
                            //    deque keeps stale triplets while the
                            //    run-origin resets). `lh3/sdust`
                            //    itself does this — the overshoot is a
                            //    port-fidelity property, not a bug.
                            // 2. Masked intervals never cover a non-
                            //    ACGT byte. They can — same reason as
                            //    (1). The spec's "perfect intervals
                            //    never span an N" rule is about
                            //    internal candidate `PerfInterval`s,
                            //    not the emitted BED intervals.
                            //
                            // Neither overshoot matters for the only
                            // consumer (`DustFilter::is_masked`): it
                            // compares positions against the interval
                            // endpoints, it never slices `seq` using
                            // them. Upstream never emits positions
                            // past chromosome length, so the
                            // out-of-bounds interval coordinates are
                            // unobservable in production.
                            for &(s, f) in &out {
                                assert!(s < f, "empty interval {}..{}", s, f);
                            }
                            for w in out.windows(2) {
                                assert!(
                                    w[0].1 <= w[1].0,
                                    "overlapping or out-of-order: {:?} {:?} (seed={} alphabet={:?} len={} w={} t={})",
                                    w[0],
                                    w[1],
                                    seed,
                                    alphabet,
                                    len,
                                    window,
                                    threshold,
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn filter_half_open_boundary_pos_at_start_masked_pos_at_finish_passes() {
        // Pin the half-open semantics of `is_masked` at both the
        // `start` and `finish` endpoints. Uses the
        // `homopolymer_with_flanks` golden snippet (mask = (50, 100)).
        let bases = GOLDEN_SNIPPETS[0].bases.to_vec();
        let len = bases.len() as u32;
        let chromosomes = vec![make_chrom("chr1", len)];
        let fetcher = StubFetcher::new(vec![bases]);
        // 1-based probes: pos = 50  → 0-based 49, just before mask → pass.
        //                  pos = 51  → 0-based 50, first masked base → drop.
        //                  pos = 100 → 0-based 99, last masked base  → drop.
        //                  pos = 101 → 0-based 100, just past mask    → pass.
        let upstream = positions_iter(vec![(0, 50), (0, 51), (0, 100), (0, 101)]);
        let cfg = DustFilterConfig::default();
        let filter = DustFilter::new(upstream, fetcher, chromosomes, cfg);
        let kept: Vec<u32> = filter.map(|r| r.unwrap().pos).collect();
        assert_eq!(
            kept,
            vec![50, 101],
            "half-open boundary: pos=51 (start) and pos=100 (finish-1) must drop; pos=50 and pos=101 must pass"
        );
    }

    #[test]
    fn filter_surfaces_invalid_pos_zero_and_latches() {
        // `pos == 0` violates the 1-based-pos contract. Must surface
        // `InvalidPos` (not silently underflow or panic) and latch.
        let chromosomes = vec![make_chrom("chr1", 128)];
        let fetcher = StubFetcher::new(vec![high_complexity_bases(128)]);
        let upstream = positions_iter(vec![(0, 0)]);
        let cfg = DustFilterConfig::default();
        let mut filter = DustFilter::new(upstream, fetcher, chromosomes, cfg);
        match filter.next() {
            Some(Err(DustFilterError::InvalidPos { chrom_id: 0 })) => {}
            other => panic!("expected InvalidPos {{ chrom_id: 0 }}, got {:?}", other),
        }
        assert!(
            filter.next().is_none(),
            "iterator must latch after InvalidPos"
        );
        assert!(filter.next().is_none(), "still latched on subsequent call");
    }

    #[test]
    fn filter_latches_after_upstream_exhaustion() {
        // Pinned behaviour: once the upstream returns None, the
        // filter latches and stays None even when polled again.
        let chromosomes = vec![make_chrom("chr1", 128)];
        let fetcher = StubFetcher::new(vec![high_complexity_bases(128)]);
        let upstream = positions_iter(vec![(0, 1)]);
        let cfg = DustFilterConfig::default();
        let mut filter = DustFilter::new(upstream, fetcher, chromosomes, cfg);
        assert!(matches!(filter.next(), Some(Ok(_))));
        assert!(filter.next().is_none(), "upstream exhausted");
        assert!(
            filter.next().is_none(),
            "polling after exhaustion must remain None"
        );
    }
}
