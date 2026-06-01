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
//! [`PileupRecord::pos`](crate::pileup_record::PileupRecord)
//! is 1-based; the internal sdust mask is 0-based half-open. The
//! conversion happens at the sweep-lookup boundary
//! ([`DustFilter::next`]) and nowhere else. `DustFilter` assumes
//! the upstream emits positions in strict monotonic order per
//! chromosome (as [`PerPositionMerger`](crate::var_calling::per_position_merger::PerPositionMerger)
//! guarantees); behaviour on contract violation is undefined and is
//! not the filter's responsibility to detect.

use std::collections::VecDeque;
use std::io;

use thiserror::Error;

use crate::fasta::{ChromRefFetchError, ChromRefFetcher};
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
/// via [`DustFilter::ensure_mask_loaded`]. Behaviour identical to
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
/// [`DustFilter::ensure_mask_loaded`], which passes a validated
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

/// Default minimum halo (bases) loaded on each side of an analysed span
/// before running DUST. Empirically large enough that the per-span mask
/// is byte-identical to a whole-contig scan at every position inside the
/// span for all >= 1 kb spans on whole-genome data (validated on the
/// tomato SL4.0 reference incl. the repeat-rich unplaced contig). The
/// `>= window` reset-barrier check in [`sdust_mask_for_span`] is the
/// correctness guarantee; this constant only sets where the search for
/// that barrier starts.
pub(crate) const MIN_DUST_HALO: u32 = 10_000;

/// True iff `buf[lo..hi]` contains a run of `>= w` consecutive `false`
/// (unmasked) entries — a DUST "reset barrier": the sliding window has
/// fully turned over with non-masking content, so the masked verdict to
/// its far side is independent of anything beyond the barrier.
fn has_unmasked_run(buf: &[bool], lo: usize, hi: usize, w: usize) -> bool {
    let mut run = 0usize;
    for &m in &buf[lo..hi.min(buf.len())] {
        if m {
            run = 0;
        } else {
            run += 1;
            if run >= w {
                return true;
            }
        }
    }
    false
}

/// Compute the DUST low-complexity mask for the analysed span
/// `[span_start, span_end)` (1-based, half-open) **without scanning the
/// whole contig** — the per-chunk replacement for a whole-chromosome
/// `sdust_mask_streaming` pre-pass.
///
/// DUST's masked verdict at a position depends on a left/right
/// neighbourhood that is bounded in practice but not by a fixed constant
/// (this sdust implementation extends masked runs across its window and
/// does not reset on sequence breaks). So we load the span plus a `halo`
/// on each side and *verify* the halo is large enough: each side's
/// padding must contain a reset barrier (a `>= window` run of unmasked
/// bases) or sit against the contig edge. If a side lacks a barrier (the
/// span is inside a low-complexity tract longer than the halo) the halo
/// is doubled and the span re-DUSTed. The starting halo is
/// `max(min_halo, window)`; on whole-genome data at >= 1 kb spans it
/// never needs to grow.
///
/// Returns `(mask, buffer_start)`: `mask` is the 1-based half-open
/// intervals **clipped to the span**, sorted and non-overlapping — the
/// shape [`partition_window`] expects; the verdict at every position
/// inside the span is byte-identical to a whole-contig
/// `sdust_mask_streaming`. `buffer_start` is the leftmost 1-based
/// position the computation actually read (after any halo expansion) —
/// the caller may `evict_before(buffer_start)` on its reference buffer
/// knowing nothing to the left of it was needed.
///
/// `fetch_bases(start_1based, len)` must return exactly `len` uppercased
/// bases of the contig starting at `start_1based`; `contig_len` is the
/// contig's total base count.
pub(crate) fn sdust_mask_for_span<F>(
    span_start: u32,
    span_end: u32,
    contig_len: u32,
    window: u32,
    threshold: u32,
    min_halo: u32,
    mut fetch_bases: F,
) -> std::io::Result<(Vec<std::ops::Range<u32>>, u32)>
where
    F: FnMut(u32, u32) -> std::io::Result<Vec<u8>>,
{
    debug_assert!(span_start >= 1, "span is 1-based");
    // The chunk's `safe_end` may be extended past the contig end (the
    // last-chunk logical extension); clamp to one-past-the-contig since
    // positions beyond the contig have no reference base (and no kept
    // records), so they carry no mask.
    let span_end = span_end.min(contig_len + 1);
    if span_start >= span_end {
        return Ok((Vec::new(), span_start));
    }
    let w = window as usize;
    let mut halo = min_halo.max(window);

    loop {
        // 1-based buffer bounds [a, b), clamped to [1, contig_len + 1).
        let a = span_start.saturating_sub(halo).max(1);
        let b = span_end.saturating_add(halo).min(contig_len + 1);
        let buf = fetch_bases(a, b - a)?;
        debug_assert_eq!(
            buf.len(),
            (b - a) as usize,
            "fetch_bases returned wrong length"
        );

        // sdust over the buffer; intervals are 0-based, buffer-relative.
        // The buffer is in memory so every base yields `Ok`; the Result
        // is propagated only for type-correctness and never errs here.
        let blen = buf.len() as u32;
        let intervals = sdust_mask_streaming(buf.iter().copied().map(Ok), blen, window, threshold)?;
        let mut buf_masked = vec![false; buf.len()];
        for (s, e) in &intervals {
            for p in *s..(*e).min(blen) {
                buf_masked[p as usize] = true;
            }
        }

        // Padding offsets within the buffer: left [0, span_start-a),
        // right [span_end-a, blen).
        let left_pad_end = (span_start - a) as usize;
        let right_pad_start = (span_end - a) as usize;
        let left_clean = a == 1 || has_unmasked_run(&buf_masked, 0, left_pad_end, w);
        let right_clean =
            b == contig_len + 1 || has_unmasked_run(&buf_masked, right_pad_start, buf.len(), w);

        if (left_clean && right_clean) || (a == 1 && b == contig_len + 1) {
            // Coalesce masked runs within [span_start, span_end), emit as
            // 1-based genomic half-open intervals (offset `off` in the
            // buffer is genomic position `a + off`).
            let lo = (span_start - a) as usize;
            let hi = (span_end - a) as usize;
            let mut out: Vec<std::ops::Range<u32>> = Vec::new();
            let mut run_start: Option<usize> = None;
            for (off, &masked) in buf_masked.iter().enumerate().take(hi).skip(lo) {
                if masked {
                    run_start.get_or_insert(off);
                } else if let Some(rs) = run_start.take() {
                    out.push((rs as u32 + a)..(off as u32 + a));
                }
            }
            if let Some(rs) = run_start {
                out.push((rs as u32 + a)..(hi as u32 + a));
            }
            return Ok((out, a));
        }
        halo = halo.saturating_mul(2);
    }
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

    /// Reference fetch failed when loading the contig's bases.
    /// The cause is exposed via `source()`, not interpolated into
    /// this variant's `Display`. After the
    /// `cohort_chrom_ref_fetcher_migration` plan landed, the source
    /// is `ChromRefFetchError` (typed) rather than the previous
    /// `io::Error` — `DustFilter` is bound to a single contig at
    /// construction, so the fetcher is a `ChromRefFetcher` impl and
    /// surfaces its typed error directly.
    #[error("reference fetch failed for chrom {chrom_id}")]
    RefFetch {
        chrom_id: u32,
        #[source]
        source: ChromRefFetchError,
    },

    /// Upstream emitted a `PileupRecord` on a different chrom than
    /// the one this `DustFilter` was constructed for. The filter is
    /// single-contig after the
    /// `cohort_chrom_ref_fetcher_migration` plan — multi-chrom
    /// records would mean the upstream wasn't chunked correctly.
    #[error("expected records on chrom_id {expected}, got chrom_id {got}")]
    ChromIdMismatch { expected: u32, got: u32 },

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

/// Streaming low-complexity filter over a per-position pileup stream.
///
/// Wraps an upstream iterator yielding
/// [`PerPositionPileups`] / [`PerPositionMergerError`] and silently
/// drops items whose reference base is inside an sdust-masked
/// interval. Emits the surviving items in the same order as the
/// upstream.
///
/// Bound to one contig at construction
/// (`bound_chrom_id`). On the first record, the filter streams the
/// contig's bases through `sdust_mask_streaming` to build the mask;
/// subsequent records on the same contig are answered by a sweep
/// over the cached interval list (`O(1)` amortised per upstream
/// position). Records on a different `chrom_id` are surfaced as
/// `DustFilterError::ChromIdMismatch` — upstream is expected to
/// have chunked records per chrom before reaching this filter (see
/// `process_one_chromosome` for the cohort path and
/// `PerChromRecordsIter` for the from-bam path).
pub struct DustFilter<I, F>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    F: ChromRefFetcher,
{
    upstream: I,
    fetcher: F,
    /// Chrom this filter is bound to. Every incoming record must
    /// match; mismatch is `DustFilterError::ChromIdMismatch`.
    bound_chrom_id: u32,
    config: DustFilterConfig,
    /// sdust mask for the bound contig. `None` before the first
    /// record arrives (mask is built lazily); `Some` after.
    mask: Option<SdustIntervals>,
    /// Sweep pointer into `mask` advancing monotonically with the
    /// upstream's per-position records. Reset to 0 when `mask` is
    /// installed.
    sweep: usize,
    /// Latches `true` after the upstream signals end or any error,
    /// so subsequent `next()` calls return `None`. Matches the
    /// upstream merger's contract.
    is_finished: bool,
}

impl<I, F> DustFilter<I, F>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    F: ChromRefFetcher,
{
    /// Wrap `upstream` with a filter bound to `bound_chrom_id`. The
    /// fetcher must be bound to the same contig (verified
    /// implicitly via the typed errors at fetch time). Infallible
    /// because `config` has already been validated by
    /// [`DustFilterConfig::new`].
    pub fn new(upstream: I, fetcher: F, bound_chrom_id: u32, config: DustFilterConfig) -> Self {
        Self {
            upstream,
            fetcher,
            bound_chrom_id,
            config,
            mask: None,
            sweep: 0,
            is_finished: false,
        }
    }

    /// Returns the validated configuration in effect.
    pub fn config(&self) -> DustFilterConfig {
        self.config
    }

    /// Build the mask if it isn't already built. Latches
    /// `is_finished` and returns the appropriate error on any
    /// failure. Streams the contig's bases through
    /// `sdust_mask_streaming` — no whole-contig buffer is
    /// materialised when the fetcher implements `iter_bases`
    /// natively (every production `ChromRefFetcher` impl does).
    fn ensure_mask_loaded(&mut self) -> Result<(), DustFilterError> {
        if self.mask.is_some() {
            return Ok(());
        }
        let length = self.fetcher.length();
        let bases = match self.fetcher.iter_bases() {
            Ok(it) => it,
            Err(source) => {
                self.is_finished = true;
                return Err(DustFilterError::RefFetch {
                    chrom_id: self.bound_chrom_id,
                    source,
                });
            }
        };
        // The inner iter yields `Result<u8, ChromRefFetchError>`,
        // but `sdust_mask_streaming` expects `Result<u8, io::Error>`.
        // Map errors at the iterator level so the streaming scan
        // surfaces ChromRefFetchError verbatim through `RefFetch`.
        //
        // We materialise the bases through a `.map(...)` because
        // sdust_mask_streaming's signature takes
        // `IntoIterator<Item = io::Result<u8>>`. Wrapping is cheap
        // (one .map adapter); the byte stream still flows lazily.
        let bases_io = bases.map(|res| res.map_err(|e| io::Error::other(format!("{e}"))));
        let mask =
            match sdust_mask_streaming(bases_io, length, self.config.window, self.config.threshold)
            {
                Ok(m) => m,
                Err(source) => {
                    self.is_finished = true;
                    return Err(DustFilterError::RefFetch {
                        chrom_id: self.bound_chrom_id,
                        source: ChromRefFetchError::Io {
                            chrom_name: format!("chrom_id {}", self.bound_chrom_id),
                            source,
                        },
                    });
                }
            };
        self.mask = Some(mask);
        self.sweep = 0;
        Ok(())
    }

    /// Returns `true` if the 0-based reference position `pos0` is
    /// inside an sdust-masked interval. Advances the cached sweep
    /// pointer past intervals whose `finish` is `≤ pos0` (they
    /// cannot match any future, larger position).
    ///
    /// Precondition: `ensure_mask_loaded` succeeded immediately
    /// before this call, so `self.mask` is `Some(_)`. The
    /// [`debug_assert!`] makes a contract violation panic loudly
    /// under tests; in release, the `None` branch passes the record
    /// through unmasked rather than producing wrong output.
    fn is_masked(&mut self, pos0: u32) -> bool {
        let Some(mask) = self.mask.as_ref() else {
            debug_assert!(false, "is_masked called without a loaded mask");
            return false;
        };
        while self.sweep < mask.len() && mask[self.sweep].1 <= pos0 {
            self.sweep += 1;
        }
        match mask.get(self.sweep) {
            Some(&(start, finish)) => start <= pos0 && pos0 < finish,
            None => false,
        }
    }
}

impl<I, F> Iterator for DustFilter<I, F>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    F: ChromRefFetcher,
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
            if pileups.chrom_id != self.bound_chrom_id {
                self.is_finished = true;
                return Some(Err(DustFilterError::ChromIdMismatch {
                    expected: self.bound_chrom_id,
                    got: pileups.chrom_id,
                }));
            }
            if let Err(err) = self.ensure_mask_loaded() {
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

    use crate::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};
    use std::cell::Cell;
    use std::io;

    /// Stub `ChromRefFetcher` backed by one contig's bytes. Counts
    /// `iter_bases` calls so tests can pin that the filter loads
    /// the mask at most once per construction.
    struct StubChromRefFetcher {
        seq: Vec<u8>,
        fail_iter: bool,
        iter_bases_calls: Cell<u32>,
    }

    impl StubChromRefFetcher {
        fn new(seq: Vec<u8>) -> Self {
            Self {
                seq,
                fail_iter: false,
                iter_bases_calls: Cell::new(0),
            }
        }

        fn with_iter_bases_failure(mut self) -> Self {
            self.fail_iter = true;
            self
        }

        fn iter_bases_calls(&self) -> u32 {
            self.iter_bases_calls.get()
        }
    }

    impl crate::fasta::fetcher::sealed::Sealed for StubChromRefFetcher {}
    impl ChromRefFetcher for StubChromRefFetcher {
        fn length(&self) -> u32 {
            self.seq.len() as u32
        }

        fn fetch(&self, start_1based: u32, length: u32) -> Result<Vec<u8>, ChromRefFetchError> {
            if start_1based == 0 {
                return Err(ChromRefFetchError::InvalidStart);
            }
            let s = (start_1based - 1) as usize;
            let e = s + length as usize;
            if e > self.seq.len() {
                return Err(ChromRefFetchError::OutOfBounds {
                    chrom_name: "stub".into(),
                    chrom_length: self.seq.len() as u32,
                    start: start_1based,
                    end: start_1based + length,
                });
            }
            Ok(self.seq[s..e].to_vec())
        }

        fn iter_bases<'a>(
            &'a self,
        ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>
        {
            self.iter_bases_calls.set(self.iter_bases_calls.get() + 1);
            if self.fail_iter {
                return Err(ChromRefFetchError::Io {
                    chrom_name: "stub".into(),
                    source: io::Error::other("stub iter_bases failure"),
                });
            }
            Ok(Box::new(self.seq.iter().copied().map(Ok)))
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
        let fetcher = StubChromRefFetcher::new(seq);
        let upstream = positions_iter((1..=10).map(|p| (0, p)).collect());
        let cfg = DustFilterConfig::default();
        let filter = DustFilter::new(upstream, fetcher, 0, cfg);
        let out: Vec<u32> = filter.map(|r| r.unwrap().pos).collect();
        assert_eq!(out, (1..=10).collect::<Vec<_>>());
    }

    #[test]
    fn filter_drops_homopolymer_positions() {
        // Whole contig is 128-bp poly-A. Every upstream position
        // should be dropped.
        let fetcher = StubChromRefFetcher::new(vec![b'A'; 128]);
        let upstream = positions_iter((1..=20).map(|p| (0, p)).collect());
        let cfg = DustFilterConfig::default();
        let filter = DustFilter::new(upstream, fetcher, 0, cfg);
        let out: Vec<_> = filter.collect();
        assert!(
            out.is_empty(),
            "expected all positions dropped, got {} survivors",
            out.len()
        );
    }

    #[test]
    fn filter_iter_bases_called_exactly_once() {
        // After the migration to single-contig DUST, the mask is
        // built lazily on the first record and reused for every
        // subsequent record. Pin that with a counter on the stub
        // fetcher: many upstream positions, one `iter_bases` call.
        let fetcher = StubChromRefFetcher::new(high_complexity_bases(128));
        let upstream = positions_iter((1..=50).map(|p| (0, p)).collect());
        let cfg = DustFilterConfig::default();
        let filter = DustFilter::new(upstream, &fetcher, 0, cfg);
        let _drain: Vec<_> = filter.collect();
        assert_eq!(fetcher.iter_bases_calls(), 1);
    }

    #[test]
    fn filter_surfaces_chrom_id_mismatch_and_latches() {
        // The filter is bound to chrom 0; upstream emits a record
        // on chrom 1. Must surface `ChromIdMismatch` and latch.
        let fetcher = StubChromRefFetcher::new(high_complexity_bases(128));
        let upstream = positions_iter(vec![(1, 5)]);
        let cfg = DustFilterConfig::default();
        let mut filter = DustFilter::new(upstream, fetcher, 0, cfg);
        match filter.next() {
            Some(Err(DustFilterError::ChromIdMismatch {
                expected: 0,
                got: 1,
            })) => {}
            other => panic!("expected ChromIdMismatch {{0,1}}, got {:?}", other),
        }
        assert!(
            filter.next().is_none(),
            "iterator must latch after mismatch"
        );
    }

    #[test]
    fn filter_surfaces_upstream_error_and_latches() {
        let fetcher = StubChromRefFetcher::new(high_complexity_bases(128));
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
        let mut filter = DustFilter::new(items.into_iter(), fetcher, 0, cfg);
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
        let fetcher = StubChromRefFetcher::new(vec![b'A'; 128]).with_iter_bases_failure();
        let upstream = positions_iter(vec![(0, 5), (0, 6)]);
        let cfg = DustFilterConfig::default();
        let mut filter = DustFilter::new(upstream, fetcher, 0, cfg);
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
    fn filter_empty_upstream_does_no_work() {
        let fetcher = StubChromRefFetcher::new(high_complexity_bases(128));
        let upstream = positions_iter(Vec::<(u32, u32)>::new());
        let cfg = DustFilterConfig::default();
        let filter = DustFilter::new(upstream, &fetcher, 0, cfg);
        let out: Vec<_> = filter.collect();
        assert!(out.is_empty());
        // The fetcher was never asked for anything.
        assert_eq!(fetcher.iter_bases_calls(), 0);
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
        let fetcher = StubChromRefFetcher::new(bases);
        // 1-based probes: pos = 50  → 0-based 49, just before mask → pass.
        //                  pos = 51  → 0-based 50, first masked base → drop.
        //                  pos = 100 → 0-based 99, last masked base  → drop.
        //                  pos = 101 → 0-based 100, just past mask    → pass.
        let upstream = positions_iter(vec![(0, 50), (0, 51), (0, 100), (0, 101)]);
        let cfg = DustFilterConfig::default();
        let filter = DustFilter::new(upstream, fetcher, 0, cfg);
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
        let fetcher = StubChromRefFetcher::new(high_complexity_bases(128));
        let upstream = positions_iter(vec![(0, 0)]);
        let cfg = DustFilterConfig::default();
        let mut filter = DustFilter::new(upstream, fetcher, 0, cfg);
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
        let fetcher = StubChromRefFetcher::new(high_complexity_bases(128));
        let upstream = positions_iter(vec![(0, 1)]);
        let cfg = DustFilterConfig::default();
        let mut filter = DustFilter::new(upstream, fetcher, 0, cfg);
        assert!(matches!(filter.next(), Some(Ok(_))));
        assert!(filter.next().is_none(), "upstream exhausted");
        assert!(
            filter.next().is_none(),
            "polling after exhaustion must remain None"
        );
    }

    /// `sdust_mask_for_span` over a span buried inside a low-complexity
    /// tract LONGER than the starting halo must (a) expand the halo until
    /// it finds an unmasked reset barrier in the high-complexity flanks,
    /// and (b) reproduce the whole-sequence mask's verdict bit-for-bit at
    /// every position in the span.
    #[test]
    fn sdust_mask_for_span_matches_whole_seq_with_expansion() {
        // high(200) + "AT"*150 (=300 bp low-complexity) + high(200).
        let mut seq = high_complexity_bases(200);
        for _ in 0..150 {
            seq.extend_from_slice(b"AT");
        }
        seq.extend(high_complexity_bases(200));
        let contig_len = seq.len() as u32; // 700
        let window = 64u32;
        let threshold = 20u32;

        // Ground truth: whole-sequence mask -> per-base bitvec (0-based).
        let gt = sdust_mask(&seq, window, threshold);
        let mut truth = vec![false; seq.len()];
        for (s, e) in &gt {
            for p in *s..*e {
                truth[p as usize] = true;
            }
        }
        // Sanity: the AT tract interior really is masked (else the test
        // would pass trivially).
        assert!(truth[300], "expected AT-tract interior to be masked");

        let fetch = |start_1based: u32, len: u32| -> std::io::Result<Vec<u8>> {
            let a = (start_1based - 1) as usize;
            Ok(seq[a..a + len as usize].to_vec())
        };

        // Span deep inside the AT tract (1-based [301, 401)); min_halo=32
        // is far smaller than the 300 bp tract, forcing 32->...->256
        // expansion to reach the unmasked flanks.
        let span_start = 301u32;
        let span_end = 401u32;
        let (mask, buf_start) = sdust_mask_for_span(
            span_start, span_end, contig_len, window, threshold, 32, fetch,
        )
        .unwrap();
        // Expansion reached back into the high-complexity left flank
        // (well before the span), so the watermark precedes the tract.
        assert!(
            buf_start < span_start,
            "expected halo expansion to read left of the span"
        );

        for pos in span_start..span_end {
            let masked_here = mask.iter().any(|r| r.start <= pos && pos < r.end);
            let truth_here = truth[(pos - 1) as usize];
            assert_eq!(
                masked_here, truth_here,
                "verdict mismatch at 1-based pos {pos}"
            );
        }

        // A high-complexity span (in the right flank) yields no mask and
        // is clean immediately (contig edge on the right, barrier left).
        let (hc, _) =
            sdust_mask_for_span(601, 701, contig_len, window, threshold, 32, fetch).unwrap();
        for pos in 601..701 {
            assert!(
                !hc.iter().any(|r| r.start <= pos && pos < r.end),
                "high-complexity pos {pos} should be unmasked"
            );
            assert!(
                !truth[(pos - 1) as usize],
                "ground truth disagrees at {pos}"
            );
        }
    }

    /// PROBE (2026-05-30 perf redesign): confirm that computing the
    /// DUST mask per-segment over `[seg.start - halo, seg.end + halo]`
    /// and clipping to the segment reproduces the whole-contig mask's
    /// masked/unmasked verdict at **every** position. This is the
    /// byte-identity gate for moving DUST from a whole-chromosome
    /// pre-pass to a per-chunk computation.
    ///
    /// Reads a real reference FASTA (env-gated, so it's a no-op in CI):
    ///   DUST_PROBE_FASTA=<path.fa> DUST_PROBE_CONTIG=<name> \
    ///   cargo test --release -- --ignored --nocapture dust_per_segment
    /// Optional: DUST_PROBE_W, DUST_PROBE_T, DUST_PROBE_SEG (segment bp).
    #[test]
    #[ignore = "reads a real reference FASTA; run manually (see doc comment)"]
    fn dust_per_segment_halo_matches_whole_contig() {
        use crate::fasta::StreamingChromRefFetcher;
        use crate::fasta::fetcher::ChromRefFetcher;

        let Ok(fasta) = std::env::var("DUST_PROBE_FASTA") else {
            eprintln!("DUST_PROBE_FASTA unset; skipping probe");
            return;
        };
        let contig = std::env::var("DUST_PROBE_CONTIG").unwrap_or_else(|_| "SL4.0ch01".into());
        let env_u32 = |k: &str, d: u32| {
            std::env::var(k)
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(d)
        };
        let window = env_u32("DUST_PROBE_W", DEFAULT_DUST_WINDOW);
        let threshold = env_u32("DUST_PROBE_T", DEFAULT_DUST_THRESHOLD);
        let seg_size = env_u32("DUST_PROBE_SEG", 100_000) as usize;

        let fetcher =
            StreamingChromRefFetcher::for_contig(std::path::Path::new(&fasta), &contig).unwrap();
        let length = fetcher.length();
        eprintln!("probe: contig={contig} length={length} W={window} T={threshold} seg={seg_size}");
        let seq = fetcher.fetch(1, length).unwrap();
        assert_eq!(seq.len(), length as usize);

        // Ground truth: whole-contig mask -> per-base bitvec.
        let gt = sdust_mask(&seq, window, threshold);
        let mut masked = vec![false; length as usize];
        for (s, e) in &gt {
            for p in *s..*e {
                masked[p as usize] = true;
            }
        }

        // FIXED-HALO sweep: when DUST_PROBE_FIXED_HALO is set, test a
        // sweep of *fixed* halos (no expansion) and report mismatches
        // for each. Answers "is a fixed N-kb halo byte-identical
        // everywhere?". Returns early (skips the barrier path).
        if std::env::var("DUST_PROBE_FIXED_HALO").is_ok() {
            let len = length as usize;
            for halo in [1000usize, 2000, 5000, 10000, 20000] {
                let mut mismatches = 0usize;
                let mut first: Vec<(usize, bool, bool)> = Vec::new();
                let mut s = 0usize;
                while s < len {
                    let e = (s + seg_size).min(len);
                    let a = s.saturating_sub(halo);
                    let b = (e + halo).min(len);
                    let seg_iv = sdust_mask(&seq[a..b], window, threshold);
                    let buf_len = (b - a) as u32;
                    let mut buf_masked = vec![false; b - a];
                    for (is, ie) in &seg_iv {
                        for p in *is..(*ie).min(buf_len) {
                            buf_masked[p as usize] = true;
                        }
                    }
                    for p in s..e {
                        if buf_masked[p - a] != masked[p] {
                            mismatches += 1;
                            if first.len() < 6 {
                                first.push((p, masked[p], buf_masked[p - a]));
                            }
                        }
                    }
                    s = e;
                }
                eprintln!(
                    "probe(fixed halo={halo:>6}): mismatches={mismatches} / {length}  first={first:?}"
                );
            }
            return;
        }

        // Barrier-termination criterion. A run of >= W consecutive
        // *unmasked* bases is a genuine sdust reset point: the window
        // (<= W bases) has fully turned over with non-masking content,
        // so positions to its right are independent of anything to its
        // left. For each segment, expand the halo until BOTH the left
        // padding [a, s) and the right padding (e, b] contain such a
        // barrier (or hit the contig edge, which is trivially clean),
        // then accept the segment's verdict. Self-terminating: the
        // whole contig is always clean.
        let len = length as usize;
        // True iff buf_masked[lo..hi] contains >= w consecutive false.
        let has_unmasked_run = |buf: &[bool], lo: usize, hi: usize, w: usize| -> bool {
            let mut run = 0usize;
            for &m in &buf[lo..hi] {
                if m {
                    run = 0;
                } else {
                    run += 1;
                    if run >= w {
                        return true;
                    }
                }
            }
            false
        };
        let w = window as usize;
        // COMBINED rule: a fixed minimum halo (kills the false-early-
        // stop failure of the pure barrier) PLUS the ≥W-unmasked-gap
        // barrier verification (catches tracts longer than the floor).
        let min_halo = env_u32("DUST_PROBE_MIN_HALO", 10_000) as usize;

        let mut mismatches = 0usize;
        let mut first: Vec<(usize, bool, bool)> = Vec::new();
        let mut max_halo = 0usize;
        let mut total_expansions = 0usize;
        let mut contig_fallbacks = 0usize; // segments that expanded to a contig edge
        let mut dust_bases: u128 = 0;
        let mut s = 0usize;
        while s < len {
            let e = (s + seg_size).min(len);
            let mut halo = w.max(min_halo); // fixed floor, then barrier-expand
            let (a, b, buf_masked) = loop {
                let a = s.saturating_sub(halo);
                let b = (e + halo).min(len);
                let seg_iv = sdust_mask(&seq[a..b], window, threshold);
                dust_bases += (b - a) as u128;
                let buf_len = (b - a) as u32;
                let mut buf_masked = vec![false; b - a];
                for (is, ie) in &seg_iv {
                    for p in *is..(*ie).min(buf_len) {
                        buf_masked[p as usize] = true;
                    }
                }
                let left_clean = a == 0 || has_unmasked_run(&buf_masked, 0, s - a, w);
                let right_clean = b == len || has_unmasked_run(&buf_masked, e - a, b - a, w);
                let full = a == 0 && b == len;
                if (left_clean && right_clean) || full {
                    if full && !(left_clean && right_clean) {
                        contig_fallbacks += 1;
                    }
                    break (a, b, buf_masked);
                }
                halo = halo.saturating_mul(2);
                total_expansions += 1;
            };
            let _ = b;
            max_halo = max_halo.max(halo);
            for p in s..e {
                if buf_masked[p - a] != masked[p] {
                    mismatches += 1;
                    if first.len() < 10 {
                        first.push((p, masked[p], buf_masked[p - a]));
                    }
                }
            }
            s = e;
        }
        let n_segs = len.div_ceil(seg_size);
        eprintln!(
            "probe(combined min_halo={min_halo}): mismatches={mismatches} / {length}  \
             max_halo={max_halo}  expansions={total_expansions} over {n_segs} segs  \
             contig_fallbacks={contig_fallbacks}  \
             dusted_bases={dust_bases} ({:.2}x contig)  first={first:?}",
            dust_bases as f64 / length as f64
        );
        assert_eq!(
            mismatches, 0,
            "combined per-segment DUST diverged from whole-contig mask"
        );
    }
}
