//! The coverage-by-GC accumulator.
//!
//! Folds the Stage-1 per-covered-position stream into the raw 2-D
//! histogram of tiled windows that becomes a [`CoverageByGcHistogram`]
//! (architecture `doc/devel/architecture/hidden_paralog_psp_integration.md`,
//! Premises 1b/2/3). The walker visits covered reference positions in
//! coordinate order on one thread, so this is a single-pass, no-allocation
//! (in steady state) fold: one running per-tile accumulator plus the count
//! matrix.
//!
//! Tiling is non-overlapping windows of `window_bp` anchored at genome
//! coordinates: a 1-based position `p` on chromosome `c` belongs to tile
//! `(c, (p - 1) / window_bp)`. A tile is finalised into one
//! `(GC fraction, covered-bases mean depth)` sample when the stream
//! crosses into a different tile (or chromosome) or at [`finish`].
//!
//! Both quantities are over the tile's **GC-defined covered positions**:
//! a position whose reference base is `N` has no GC and is excluded from
//! the covered count, the GC count, and the depth sum. A tile with no such
//! position yields no sample and is counted in `n_skipped_tiles`.
//!
//! [`finish`]: CoverageByGcAccumulator::finish

use std::collections::VecDeque;

use super::CoverageByGcHistogram;

/// Bin scheme for the coverage-by-GC histogram. These are the tuning
/// knobs (architecture Premise 2/3); the `pileup` CLI supplies them
/// (C1) and validates them at that boundary. Every field must be
/// positive (`depth_bin_width` finite and `> 0`); [`assert_valid`] asserts
/// this in both debug and release rather than silently using a scheme the
/// caller never chose.
///
/// [`assert_valid`]: CoverageBinScheme::assert_valid
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoverageBinScheme {
    /// Tile width in bp (the GC covariate window = analysis window).
    pub window_bp: u32,
    /// Number of GC bins, uniform over `[0, 1]`.
    pub gc_bins: u32,
    /// Width of one depth bin in covered-bases-mean-depth units.
    pub depth_bin_width: f64,
    /// Number of regular depth bins; one overflow bin follows.
    pub depth_bins: u32,
}

impl CoverageBinScheme {
    /// Cells per GC row: the regular depth bins plus one overflow bin.
    fn depth_cols(&self) -> usize {
        self.depth_bins as usize + 1
    }

    /// Total histogram cells = `gc_bins * (depth_bins + 1)`.
    fn n_cells(&self) -> usize {
        self.gc_bins as usize * self.depth_cols()
    }

    /// Panic (in both debug and release) if any field is non-positive
    /// (`depth_bin_width` not finite or `<= 0`). The scheme is producer config
    /// validated at the CLI boundary (C1); an invalid scheme reaching an
    /// accumulator is a programmer error, so it fails loudly rather than
    /// silently substituting a different scheme. Shared by both accumulators'
    /// constructors.
    pub(crate) fn assert_valid(&self) {
        assert!(self.window_bp >= 1, "window_bp must be >= 1");
        assert!(self.gc_bins >= 1, "gc_bins must be >= 1");
        assert!(self.depth_bins >= 1, "depth_bins must be >= 1");
        assert!(
            self.depth_bin_width.is_finite() && self.depth_bin_width > 0.0,
            "depth_bin_width must be finite and > 0, got {}",
            self.depth_bin_width,
        );
    }

    /// Map a window's `(GC fraction in [0, 1], mean depth >= 0)` to a
    /// row-major cell index. GC saturates to the last GC bin at 1.0; depth at
    /// or above `depth_bins * depth_bin_width` lands in the overflow column.
    /// `mean_depth` is always `>= 0` (a sum of non-negative counts over a
    /// positive divisor), so `0` maps to depth bin `0`; the `as usize` cast also
    /// saturates a hypothetical negative to `0`.
    fn cell_index(&self, gc_frac: f64, mean_depth: f64) -> usize {
        let gc_bins = self.gc_bins as usize;
        let gc_bin = ((gc_frac * gc_bins as f64) as usize).min(gc_bins - 1);
        let depth_bins = self.depth_bins as usize;
        let depth_bin = ((mean_depth / self.depth_bin_width) as usize).min(depth_bins);
        gc_bin * self.depth_cols() + depth_bin
    }
}

/// The running per-tile fold state. `covered` excludes `N` positions, so
/// `gc / covered` and `depth_sum / covered` are over GC-defined covered
/// positions only.
#[derive(Debug, Clone, Copy)]
struct OpenTile {
    /// `(chrom_id, tile_index)` this state belongs to.
    key: (u32, u32),
    /// GC-defined covered positions in the tile.
    covered: u64,
    /// Of `covered`, how many have a `G`/`C` reference base.
    gc: u64,
    /// Sum of per-position depth over the `covered` positions.
    depth_sum: u64,
}

/// Accumulates the coverage-by-GC histogram across one sample's
/// per-covered-position stream. Build with [`new`], feed each covered
/// position with [`observe`], and reduce with [`finish`].
///
/// [`new`]: Self::new
/// [`observe`]: Self::observe
/// [`finish`]: Self::finish
#[derive(Debug, Clone)]
pub struct CoverageByGcAccumulator {
    scheme: CoverageBinScheme,
    /// Row-major `[gc_bin][depth_bin]` counts, length `scheme.n_cells()`.
    counts: Vec<u32>,
    /// The tile currently being folded; `None` before the first position.
    open: Option<OpenTile>,
    /// `(chrom_id, pos)` of the previous [`observe`](Self::observe) call,
    /// used to debug-assert the walker's coordinate-order invariant.
    last_observed: Option<(u32, u32)>,
    /// Non-skipped tiles folded so far. In this (legacy, tiled) accumulator a
    /// tile spans many positions, so this is a tile count — it feeds the
    /// histogram's position-named `n_positions` field, which in the live
    /// sliding-window model is genuinely one-per-position.
    n_tiles: u64,
    n_skipped_tiles: u64,
    /// Running grand total of GC-defined (non-`N`) covered positions across
    /// every folded tile — the sample's callable-position count, kept here
    /// because the histogram cells lose each tile's covered count once the
    /// tile collapses to one `(GC, mean depth)` sample.
    callable_positions: u64,
}

impl CoverageByGcAccumulator {
    /// Construct an empty accumulator for the given bin scheme.
    ///
    /// # Panics
    ///
    /// Panics (in both debug and release) if any scheme field is
    /// non-positive (`depth_bin_width` not finite or `<= 0`). The scheme
    /// is producer config validated at the CLI boundary (C1); an invalid
    /// scheme here is a programmer error, so it fails loudly rather than
    /// silently substituting a different scheme.
    pub fn new(scheme: CoverageBinScheme) -> Self {
        scheme.assert_valid();
        Self {
            counts: vec![0; scheme.n_cells()],
            scheme,
            open: None,
            last_observed: None,
            n_tiles: 0,
            n_skipped_tiles: 0,
            callable_positions: 0,
        }
    }

    /// Fold one covered reference position. `ref_base` is the position's
    /// reference base (any case); `depth` is its total fragment depth
    /// (Σ allele observations).
    ///
    /// Positions must arrive in non-decreasing `(chrom_id, pos)` order —
    /// the walker's coordinate-order invariant — so a tile, once left, is
    /// never re-entered. This is debug-asserted; out-of-order input is a
    /// programmer error that would double-count tiles.
    pub fn observe(&mut self, chrom_id: u32, pos: u32, ref_base: u8, depth: u32) {
        debug_assert!(
            self.last_observed
                .is_none_or(|prev| (chrom_id, pos) >= prev),
            "coverage observe out of order: {:?} after {:?}",
            (chrom_id, pos),
            self.last_observed,
        );
        self.last_observed = Some((chrom_id, pos));

        // Tile index from the 1-based position; pos 0 is not a valid
        // 1-based coordinate but saturates to tile 0 rather than wrapping.
        let tile_index = pos.saturating_sub(1) / self.scheme.window_bp;
        let key = (chrom_id, tile_index);

        match self.open {
            Some(tile) if tile.key == key => {}
            _ => {
                // Crossed into a new tile (or first position): finalise the
                // one we were folding, then open this tile fresh.
                if let Some(tile) = self.open.take() {
                    self.finalise_tile(tile);
                }
                self.open = Some(OpenTile {
                    key,
                    covered: 0,
                    gc: 0,
                    depth_sum: 0,
                });
            }
        }

        // `N` reference bases have no defined GC, so they contribute to no
        // count — not even `covered`. Everything else is a GC-defined
        // covered position.
        if ref_base.eq_ignore_ascii_case(&b'N') {
            return;
        }
        // PANIC-FREE: the match arm above sets `self.open` to `Some` on
        // every path (it either matched an existing open tile or just
        // assigned a fresh one), and nothing between there and here clears
        // it, so the tile is present.
        let tile = self
            .open
            .as_mut()
            .expect("open tile was just set above for this position");
        tile.covered += 1;
        if matches!(ref_base.to_ascii_uppercase(), b'G' | b'C') {
            tile.gc += 1;
        }
        // `depth` is a per-position fragment count; `saturating_add` makes
        // the (practically unreachable) accumulation overflow explicit.
        tile.depth_sum = tile.depth_sum.saturating_add(u64::from(depth));
    }

    /// Finalise the last open tile and return the histogram.
    pub fn finish(mut self) -> CoverageByGcHistogram {
        if let Some(tile) = self.open.take() {
            self.finalise_tile(tile);
        }
        CoverageByGcHistogram {
            window_bp: self.scheme.window_bp,
            gc_bins: self.scheme.gc_bins,
            depth_bin_width: self.scheme.depth_bin_width,
            depth_bins: self.scheme.depth_bins,
            n_positions: self.n_tiles,
            n_skipped_tiles: self.n_skipped_tiles,
            callable_positions: self.callable_positions,
            counts: self.counts,
        }
    }

    /// Bin a completed tile's `(GC fraction, mean depth)` into the matrix,
    /// or count it as skipped when it has no GC-defined covered position.
    fn finalise_tile(&mut self, tile: OpenTile) {
        if tile.covered == 0 {
            self.n_skipped_tiles += 1;
            return;
        }
        // Sum the tile's covered positions into the callable grand total
        // before it collapses to a single histogram cell. Skipped tiles
        // (handled above) contribute nothing.
        self.callable_positions = self.callable_positions.saturating_add(tile.covered);
        let gc_frac = tile.gc as f64 / tile.covered as f64;
        let mean_depth = tile.depth_sum as f64 / tile.covered as f64;
        let cell = self.scheme.cell_index(gc_frac, mean_depth);
        self.counts[cell] += 1;
        self.n_tiles += 1;
    }
}

// ---------------------------------------------------------------------------
// Sliding-window (centred) coverage — the paralog per-position column source.
// ---------------------------------------------------------------------------

/// One covered position's centred-window coverage summary: the mean read depth
/// and GC fraction over the `window_bp`-wide window **centred** on the position,
/// both over the window's GC-defined (non-`N`) covered positions.
///
/// Emitted by [`SlidingWindowCoverageAccumulator`] as each position's window
/// becomes complete (i.e. once the stream has advanced `window_bp / 2` past it).
/// This is the value stored per position in the `.psp` and looked up, unbinned,
/// by the hidden-paralog score — the same value that is also folded into the
/// coverage-by-GC histogram used to fit the coverage model, so the yardstick
/// (σ₀) and the observation share one window definition
/// (`doc/devel/architecture/hidden_paralog_pileup_window_coverage.md`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowCoverage {
    /// Chromosome the position (and its window) lie on.
    pub chrom_id: u32,
    /// 1-based reference position the window is centred on. Always a
    /// GC-defined covered position (`N` positions emit no sample).
    pub pos: u32,
    /// Mean read depth over the window's covered positions.
    pub mean_depth: f32,
    /// GC fraction (`G`/`C` reference bases / covered positions) over the window.
    pub gc_fraction: f32,
}

/// A covered position retained in the sliding-window buffer.
#[derive(Debug, Clone, Copy)]
struct CoveredPosition {
    pos: u32,
    depth: u32,
    is_gc: bool,
}

/// Accumulates, for every GC-defined covered position, the mean depth and GC
/// fraction of a **centred** window of width `window_bp` (`± window_bp / 2`) over
/// the covered positions in that window — the paralog score's per-position
/// coverage input.
///
/// Fed the same one-thread, coordinate-ordered per-covered-position stream as
/// [`CoverageByGcAccumulator`], via [`observe`](Self::observe). A position `p` is
/// **finalised** — its centred window is complete — once the stream reaches
/// `p + window_bp / 2`; at that point its [`WindowCoverage`] is queued for the
/// caller to [`pop_ready`](Self::pop_ready) *and* folded into the coverage-by-GC
/// histogram (so a single pass serves both the per-position column and the model
/// fit). [`finish`](Self::finish) finalises the tail (windows truncated at the
/// last observed position / a contig end) and returns it plus the histogram in
/// one consuming call, so the tail can never be silently dropped.
///
/// Memory is `O(window_bp)`: the buffer holds only positions within the active
/// window span. Windows never span chromosomes — a contig change finalises the
/// previous contig's tail and resets the sliding state.
///
/// # Invariants (the two-pointer window state, held between `observe` calls)
///
/// The covered positions of the current contig are appended to `buf` in
/// coordinate order; `buf` holds the suffix still needed — from the left edge of
/// the next centre's window to the newest observed position. Two cursors index
/// into it, and the running sums cover a contiguous prefix:
///
/// - `centre_offset` — index in `buf` of the next centre to finalise.
/// - `summed_offset` — one past the last `buf` entry folded into
///   `sum_depth`/`sum_gc`/`count`; the sums equal `buf[0..summed_offset]`.
/// - Ordering: `centre_offset <= summed_offset <= buf.len()`.
/// - When the shrink-left step pops the front of `buf` (a position now left of
///   every remaining centre's window), it decrements **both** cursors in
///   lockstep so they keep pointing at the same logical entries — safe because a
///   front is popped only when it sits strictly left of the current centre (so
///   it was already summed and already before the centre), keeping both cursors
///   `>= 1` at the decrement.
#[derive(Debug, Clone)]
pub struct SlidingWindowCoverageAccumulator {
    scheme: CoverageBinScheme,
    /// Half-window in bp (`window_bp / 2`); the window centred on `p` spans the
    /// covered positions in `[p − half, p + half]`.
    half: u32,
    /// Row-major `[gc_bin][depth_bin]` counts, length `scheme.n_cells()`.
    counts: Vec<u32>,
    /// GC-defined covered positions finalised so far (one per emitted window).
    /// In the sliding model this is *both* the histogram's `n_positions` and its
    /// `callable_positions` — every finalised centre is exactly one covered
    /// position — so a single counter feeds both.
    n_covered_positions: u64,
    /// Contig currently being folded; `None` before the first position.
    chrom_id: Option<u32>,
    /// Covered positions retained (see the type's `# Invariants`).
    buf: VecDeque<CoveredPosition>,
    /// Index in `buf` of the next centre to finalise (see `# Invariants`).
    centre_offset: usize,
    /// One past the last `buf` entry folded into the running sums
    /// (see `# Invariants`).
    summed_offset: usize,
    /// Running sums over `buf[0..summed_offset]` — the current centre's window
    /// once the shrink-left step has removed entries left of it.
    sum_depth: u64,
    sum_gc: u64,
    count: u64,
    /// Finalised windows awaiting the caller.
    ready: VecDeque<WindowCoverage>,
    /// `(chrom_id, pos)` of the previous [`observe`](Self::observe), for the
    /// coordinate-order debug assert.
    last_observed: Option<(u32, u32)>,
}

impl SlidingWindowCoverageAccumulator {
    /// Construct an empty accumulator for the given bin scheme. Panics on a
    /// non-positive scheme field (via [`CoverageBinScheme::assert_valid`]),
    /// exactly as [`CoverageByGcAccumulator::new`].
    pub fn new(scheme: CoverageBinScheme) -> Self {
        scheme.assert_valid();
        Self {
            half: scheme.window_bp / 2,
            counts: vec![0; scheme.n_cells()],
            scheme,
            n_covered_positions: 0,
            chrom_id: None,
            buf: VecDeque::new(),
            centre_offset: 0,
            summed_offset: 0,
            sum_depth: 0,
            sum_gc: 0,
            count: 0,
            ready: VecDeque::new(),
            last_observed: None,
        }
    }

    /// Fold one covered reference position. `ref_base` is the position's
    /// reference base (any case); `depth` is its total fragment depth.
    ///
    /// Positions must arrive in non-decreasing `(chrom_id, pos)` order (the
    /// walker's coordinate-order invariant), debug-asserted. An `N` reference
    /// base is *not* a covered position: it does not become a centre and does
    /// not contribute to any window's sums (matching the tile accumulator's
    /// GC-defined rule), but it still advances the finalisation frontier.
    /// Finalised windows are queued for [`pop_ready`](Self::pop_ready).
    pub fn observe(&mut self, chrom_id: u32, pos: u32, ref_base: u8, depth: u32) {
        debug_assert!(
            self.last_observed
                .is_none_or(|prev| (chrom_id, pos) >= prev),
            "sliding coverage observe out of order: {:?} after {:?}",
            (chrom_id, pos),
            self.last_observed,
        );
        self.last_observed = Some((chrom_id, pos));

        // Contig change: no window spans chromosomes, so drain the previous
        // contig's remaining centres (truncated) and reset the sliding state.
        if self.chrom_id != Some(chrom_id) {
            self.finalise_all();
            self.reset_contig();
            self.chrom_id = Some(chrom_id);
        }

        // A non-`N` covered position joins the buffer as a centre + a
        // sum contributor; an `N` position only advances the frontier below.
        if !ref_base.eq_ignore_ascii_case(&b'N') {
            self.buf.push_back(CoveredPosition {
                pos,
                depth,
                is_gc: matches!(ref_base.to_ascii_uppercase(), b'G' | b'C'),
            });
        }

        // Finalise every centre whose right edge `centre + half` the frontier
        // `pos` has now reached — no later covered position can fall in its
        // window, so it is complete.
        while self.centre_offset < self.buf.len()
            && self.buf[self.centre_offset].pos.saturating_add(self.half) <= pos
        {
            self.finalise_centre();
        }
    }

    /// Pop one finalised window, oldest first, or `None` if none is ready.
    pub fn pop_ready(&mut self) -> Option<WindowCoverage> {
        self.ready.pop_front()
    }

    /// Finalise the stream: finalise every remaining centre (windows truncated
    /// at the last observed position / a contig end) and return, in one
    /// consuming call, the still-queued finalised windows (those not yet drained
    /// via [`pop_ready`](Self::pop_ready), plus the tail just finalised) together
    /// with the coverage-by-GC histogram. Consuming `self` and returning the tail
    /// means it can never be silently dropped — a caller cannot forget to drain.
    pub fn finish(mut self) -> (Vec<WindowCoverage>, CoverageByGcHistogram) {
        self.finalise_all();
        let tail: Vec<WindowCoverage> = self.ready.into_iter().collect();
        let histogram = CoverageByGcHistogram {
            window_bp: self.scheme.window_bp,
            gc_bins: self.scheme.gc_bins,
            depth_bin_width: self.scheme.depth_bin_width,
            depth_bins: self.scheme.depth_bins,
            // Every finalised centre is one covered position — there are no
            // "skipped tiles" in the sliding model, so `n_positions == callable ==
            // n_covered_positions` and `n_skipped_tiles == 0`.
            n_positions: self.n_covered_positions,
            n_skipped_tiles: 0,
            callable_positions: self.n_covered_positions,
            counts: self.counts,
        };
        (tail, histogram)
    }

    /// Finalise the centre at `centre_offset`: complete its window (extend the
    /// running sums right to `centre + half`, shrink them left past
    /// `centre − half`), emit the `(GC, mean depth)` sample, and fold it into
    /// the histogram. `count >= 1` always — the centre itself lies in its window.
    fn finalise_centre(&mut self) {
        let centre = self.buf[self.centre_offset].pos;
        // `saturating_add` clamps a centre within `half` of `u32::MAX` to
        // `u32::MAX` — unreachable at genomic scale (`u32::MAX` ≈ 4.29 Gbp is
        // larger than any contig), matching the `pos == 0` saturation elsewhere.
        let hi = centre.saturating_add(self.half);
        let lo = centre.saturating_sub(self.half);

        // Extend right: pull in every buffered position at or before the
        // window's right edge (all are already observed — the frontier is here).
        while self.summed_offset < self.buf.len() && self.buf[self.summed_offset].pos <= hi {
            let e = self.buf[self.summed_offset];
            self.sum_depth += u64::from(e.depth);
            self.sum_gc += u64::from(e.is_gc);
            self.count += 1;
            self.summed_offset += 1;
        }
        // Shrink left: drop every buffered position before the window's left
        // edge (also before every later centre's window → never needed again).
        // The `-=` cannot underflow: a popped front was summed by the extend
        // step above and sits left of the centre, so `count`/`summed_offset`/
        // `centre_offset` are all `>= 1` here (see the type's `# Invariants`).
        while let Some(front) = self.buf.front() {
            if front.pos < lo {
                self.sum_depth -= u64::from(front.depth);
                self.sum_gc -= u64::from(front.is_gc);
                self.count -= 1;
                self.buf.pop_front();
                self.summed_offset -= 1;
                self.centre_offset -= 1;
            } else {
                break;
            }
        }

        let gc_frac = self.sum_gc as f64 / self.count as f64;
        let mean_depth = self.sum_depth as f64 / self.count as f64;
        self.ready.push_back(WindowCoverage {
            // PANIC-FREE: `observe` sets `chrom_id = Some(_)` before buffering
            // any position, and `reset_contig` deliberately keeps it set, so a
            // centre is only ever finalised while a contig is active.
            chrom_id: self.chrom_id.expect("centre finalised within a contig"),
            pos: centre,
            // The stored `f32` holds a *mean* depth (bounded by per-position
            // depth, far under `f32`'s ~1.6e7 exact-integer limit), so the
            // narrowing is lossless in practice; the `.psp` column is `f32` too.
            mean_depth: mean_depth as f32,
            gc_fraction: gc_frac as f32,
        });
        self.counts[self.scheme.cell_index(gc_frac, mean_depth)] += 1;
        self.n_covered_positions += 1;
        self.centre_offset += 1;
    }

    /// Finalise all remaining centres (truncated windows).
    fn finalise_all(&mut self) {
        while self.centre_offset < self.buf.len() {
            self.finalise_centre();
        }
    }

    /// Reset the per-contig sliding state (called after draining a contig).
    /// `chrom_id` is deliberately *not* cleared — the caller sets it to the new
    /// contig immediately after.
    fn reset_contig(&mut self) {
        self.buf.clear();
        self.centre_offset = 0;
        self.summed_offset = 0;
        self.sum_depth = 0;
        self.sum_gc = 0;
        self.count = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scheme() -> CoverageBinScheme {
        CoverageBinScheme {
            window_bp: 10,
            gc_bins: 2,
            depth_bin_width: 1.0,
            depth_bins: 4,
        }
    }

    /// A single tile of uniform depth and 50% GC lands in exactly one
    /// cell: GC bin 1 (0.5 of 2 bins) and depth bin 6 -> overflow? No,
    /// depth 6 with width 1, depth_bins 4 -> overflow bin (index 4).
    #[test]
    fn single_tile_lands_in_expected_cell() {
        let mut acc = CoverageByGcAccumulator::new(scheme());
        // 4 positions in tile (0,0): A, C, G, T -> 50% GC; depth 2 each.
        for (p, b) in [(1u32, b'A'), (2, b'C'), (3, b'G'), (4, b'T')] {
            acc.observe(0, p, b, 2);
        }
        let h = acc.finish();
        assert_eq!(h.n_positions, 1);
        assert_eq!(h.n_skipped_tiles, 0);
        // GC = 0.5 -> gc_bin = floor(0.5 * 2) = 1. depth = 2.0, width 1 ->
        // depth_bin = 2. cell = 1 * (4 + 1) + 2 = 7.
        let total: u32 = h.counts.iter().sum();
        assert_eq!(total, 1);
        assert_eq!(h.counts[7], 1);
    }

    /// `callable_positions` is the grand total of GC-defined (non-`N`)
    /// covered positions across every folded tile — independent of how the
    /// tiles bin into the histogram. Here: two tiles, 3 non-`N` + 2 non-`N`
    /// covered positions (the interleaved `N`s do not count), so the total
    /// is 5.
    #[test]
    fn callable_positions_sum_covered_non_n() {
        let mut acc = CoverageByGcAccumulator::new(scheme());
        // Tile (0,0): positions 1..=10. 3 non-N (A, C, G) + 1 N.
        acc.observe(0, 1, b'A', 1);
        acc.observe(0, 2, b'N', 9);
        acc.observe(0, 3, b'C', 1);
        acc.observe(0, 4, b'G', 1);
        // Tile (0,1): positions 11..=20. 2 non-N (T, a) + 1 lowercase n.
        acc.observe(0, 11, b'T', 1);
        acc.observe(0, 12, b'n', 9);
        acc.observe(0, 13, b'a', 1);
        let h = acc.finish();
        assert_eq!(h.callable_positions, 5);
        assert_eq!(h.n_positions, 2);
    }

    /// The callable total sums covered positions across *every* tile,
    /// including the final still-open tile drained by [`finish`]. Three tiles
    /// with distinct counts (2 + 3 + 1) separate "sum all tiles" from
    /// "sum last tile only" or "drop the final open tile" — a bug in any of
    /// those would miscount here.
    #[test]
    fn callable_positions_counts_every_covered_position_across_many_tiles() {
        let mut acc = CoverageByGcAccumulator::new(scheme()); // window_bp 10
        // tile 0: 2 covered (pos 1,2); tile 1: 3 covered (11,12,13);
        // tile 2: 1 covered (21) — left OPEN, drained by finish().
        for (c, p) in [(0u32, 1u32), (0, 2), (0, 11), (0, 12), (0, 13), (0, 21)] {
            acc.observe(c, p, b'A', 1);
        }
        let h = acc.finish();
        assert_eq!(h.callable_positions, 6);
        assert_eq!(h.n_positions, 3);
    }

    /// A tile with only `N` covered positions contributes nothing to the
    /// callable total (and is counted as skipped), so an all-`N` sample has
    /// `callable_positions == 0`.
    #[test]
    fn all_n_tile_adds_zero_callable() {
        let mut acc = CoverageByGcAccumulator::new(scheme());
        acc.observe(0, 1, b'N', 10);
        acc.observe(0, 2, b'N', 10);
        let h = acc.finish();
        assert_eq!(h.callable_positions, 0);
        assert_eq!(h.n_skipped_tiles, 1);
    }

    /// `N` positions are excluded from covered/gc/depth, shifting the GC
    /// fraction and mean depth of the surviving positions.
    #[test]
    fn n_positions_are_excluded() {
        let mut acc = CoverageByGcAccumulator::new(scheme());
        // One G (depth 4) and two N (depth 100) -> covered = 1, gc = 1,
        // mean depth = 4, GC fraction = 1.0.
        acc.observe(0, 1, b'G', 4);
        acc.observe(0, 2, b'N', 100);
        acc.observe(0, 3, b'n', 100); // lowercase N too
        let h = acc.finish();
        assert_eq!(h.n_positions, 1);
        // GC = 1.0 -> gc_bin clamps to 1. depth 4.0 -> depth_bin 4 (overflow).
        // cell = 1 * 5 + 4 = 9.
        assert_eq!(h.counts[9], 1);
        assert_eq!(h.counts.iter().sum::<u32>(), 1);
    }

    /// A tile whose only covered positions are `N` yields no sample and is
    /// counted as skipped.
    #[test]
    fn all_n_tile_is_skipped() {
        let mut acc = CoverageByGcAccumulator::new(scheme());
        acc.observe(0, 1, b'N', 10);
        acc.observe(0, 2, b'N', 10);
        let h = acc.finish();
        assert_eq!(h.n_positions, 0);
        assert_eq!(h.n_skipped_tiles, 1);
        assert_eq!(h.counts.iter().sum::<u32>(), 0);
    }

    /// Crossing into a new tile (or chromosome) finalises the previous one.
    #[test]
    fn tile_and_chromosome_boundaries_finalise() {
        let mut acc = CoverageByGcAccumulator::new(scheme());
        // window_bp = 10: positions 1..=10 are tile 0; 11 is tile 1.
        acc.observe(0, 1, b'A', 1); // tile (0,0)
        acc.observe(0, 11, b'A', 1); // tile (0,1) -> finalises (0,0)
        acc.observe(1, 1, b'A', 1); // chrom 1 tile (1,0) -> finalises (0,1)
        let h = acc.finish(); // finalises (1,0)
        assert_eq!(h.n_positions, 3);
    }

    /// Depth at or above `depth_bins * depth_bin_width` lands in the
    /// overflow column; below the bottom edge lands in bin 0.
    #[test]
    fn depth_overflow_and_zero_bins() {
        let mut acc = CoverageByGcAccumulator::new(scheme());
        // tile 0: depth 0 -> depth_bin 0.
        acc.observe(0, 1, b'A', 0);
        // tile 1: depth 1000 -> overflow depth_bin 4.
        acc.observe(0, 11, b'A', 1000);
        let h = acc.finish();
        // GC = 0 for both -> gc_bin 0. cells: 0*5+0 = 0 and 0*5+4 = 4.
        assert_eq!(h.counts[0], 1);
        assert_eq!(h.counts[4], 1);
    }

    /// The finished histogram passes the B1 shape/validate invariants.
    #[test]
    fn finished_histogram_is_valid() {
        let mut acc = CoverageByGcAccumulator::new(scheme());
        acc.observe(0, 1, b'G', 3);
        acc.observe(0, 2, b'C', 3);
        let h = acc.finish();
        assert_eq!(
            h.counts.len(),
            (h.gc_bins as usize) * (h.depth_bins as usize + 1)
        );
        // Wrap into a SampleSummary and validate end-to-end.
        let summary = super::super::SampleSummary {
            version: super::super::SAMPLE_SUMMARY_VERSION,
            coverage_by_gc: h,
            heterozygosity: super::super::HetCounts {
                n_het_sites: 0,
                n_hom_alt_sites: 0,
                n_ambiguous_sites: 0,
                n_variant_sites: 0,
                min_depth: 1,
                error_rate: 0.02,
                lr_margin: std::f64::consts::LN_10,
                strand_bias_z: 3.0,
            },
        };
        summary
            .validate()
            .expect("accumulator output must validate");
    }

    /// An empty accumulator (no positions) finishes to an all-zero, valid
    /// histogram with no tiles.
    #[test]
    fn empty_accumulator_finishes_clean() {
        let h = CoverageByGcAccumulator::new(scheme()).finish();
        assert_eq!(h.n_positions, 0);
        assert_eq!(h.n_skipped_tiles, 0);
        assert_eq!(h.counts.iter().sum::<u32>(), 0);
        assert_eq!(h.counts.len(), 2 * 5);
    }

    /// Mean depth exactly on the top regular bin edge
    /// (`depth_bins * depth_bin_width`) lands in the overflow column, not
    /// the last regular bin — pins the `min(.., depth_bins)` boundary.
    #[test]
    fn depth_exactly_on_top_edge_is_overflow() {
        // window 10, depth_bins 4, width 1.0 -> top edge = 4.0.
        let mut acc = CoverageByGcAccumulator::new(scheme());
        acc.observe(0, 1, b'A', 4); // mean depth 4.0 (single position)
        let h = acc.finish();
        // GC 0 -> gc_bin 0; depth 4.0 -> overflow bin 4. cell = 4.
        assert_eq!(h.counts[4], 1);
        assert_eq!(h.counts.iter().sum::<u32>(), 1);
    }

    /// `pos == 0` (not a valid 1-based coordinate) saturates to tile 0
    /// rather than wrapping, and still bins.
    #[test]
    fn pos_zero_saturates_to_tile_zero() {
        let mut acc = CoverageByGcAccumulator::new(scheme());
        acc.observe(0, 0, b'A', 1);
        let h = acc.finish();
        assert_eq!(h.n_positions, 1);
        assert_eq!(h.counts.iter().sum::<u32>(), 1);
    }

    /// Out-of-order input violates the walker's coordinate-order invariant
    /// and panics in debug (where tests run), surfacing the programmer
    /// error rather than silently double-counting tiles.
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "out of order")]
    fn out_of_order_observe_panics_in_debug() {
        let mut acc = CoverageByGcAccumulator::new(scheme());
        acc.observe(0, 20, b'A', 1); // tile 1
        acc.observe(0, 5, b'A', 1); // tile 0 — backwards
    }

    /// A non-positive scheme is a programmer error: `new` panics in both
    /// debug and release rather than silently substituting a scheme.
    #[test]
    #[should_panic(expected = "gc_bins must be >= 1")]
    fn zero_gc_bins_scheme_panics() {
        let mut s = scheme();
        s.gc_bins = 0;
        let _ = CoverageByGcAccumulator::new(s);
    }

    // -- sliding-window accumulator ----------------------------------------

    /// A `window_bp = 4` scheme (`half = 2`): the window centred on `p` spans
    /// covered positions in `[p − 2, p + 2]`. Coarse GC/depth bins keep the
    /// per-value asserts about the emitted `WindowCoverage`, not the histogram.
    fn sliding_scheme() -> CoverageBinScheme {
        CoverageBinScheme {
            window_bp: 4,
            gc_bins: 2,
            depth_bin_width: 1.0,
            depth_bins: 40,
        }
    }

    /// Feed a `(pos, ref_base, depth)` stream on one contig and collect every
    /// emitted window (drained after each observe + a final `finish_stream`),
    /// plus the finished histogram.
    fn run_sliding(
        scheme: CoverageBinScheme,
        chrom_id: u32,
        stream: &[(u32, u8, u32)],
    ) -> (Vec<WindowCoverage>, CoverageByGcHistogram) {
        let mut acc = SlidingWindowCoverageAccumulator::new(scheme);
        let mut out = Vec::new();
        for &(pos, base, depth) in stream {
            acc.observe(chrom_id, pos, base, depth);
            while let Some(wc) = acc.pop_ready() {
                out.push(wc);
            }
        }
        let (tail, hist) = acc.finish();
        out.extend(tail);
        (out, hist)
    }

    /// A ramp `depth = pos` over positions 1..=10, all `G`, makes every centred
    /// window's mean hand-computable — and on a symmetric window over a linear
    /// ramp the interior means equal the centre, the edges are truncated.
    #[test]
    fn sliding_ramp_means_match_hand_computed() {
        let stream: Vec<(u32, u8, u32)> = (1..=10u32).map(|p| (p, b'G', p)).collect();
        let (windows, hist) = run_sliding(sliding_scheme(), 0, &stream);

        // One window per covered position, in position order.
        let positions: Vec<u32> = windows.iter().map(|w| w.pos).collect();
        assert_eq!(positions, (1..=10).collect::<Vec<_>>());

        // Hand-computed centred means over [p-2, p+2] ∩ [1,10]:
        let expected = [2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0];
        for (w, &e) in windows.iter().zip(expected.iter()) {
            assert!(
                (w.mean_depth - e).abs() < 1e-6,
                "pos {} mean {} != {e}",
                w.pos,
                w.mean_depth,
            );
            assert!((w.gc_fraction - 1.0).abs() < 1e-6, "all-G → gc 1.0");
        }
        // Every finalised centre is one covered position; no skipped tiles.
        assert_eq!(hist.n_positions, 10);
        assert_eq!(hist.callable_positions, 10);
        assert_eq!(hist.n_skipped_tiles, 0);
    }

    /// A uniform stream: every window mean equals the constant depth and GC,
    /// so all samples land in a single histogram cell (a flat sample).
    #[test]
    fn sliding_uniform_all_one_cell() {
        let stream: Vec<(u32, u8, u32)> = (1..=20u32).map(|p| (p, b'C', 12)).collect();
        let (windows, hist) = run_sliding(sliding_scheme(), 3, &stream);
        assert_eq!(windows.len(), 20);
        for w in &windows {
            assert!((w.mean_depth - 12.0).abs() < 1e-6);
            assert!((w.gc_fraction - 1.0).abs() < 1e-6);
            assert_eq!(w.chrom_id, 3);
        }
        let nonzero: Vec<u32> = hist.counts.iter().copied().filter(|&c| c > 0).collect();
        assert_eq!(nonzero, vec![20], "all 20 windows in one cell");
    }

    /// An `N` reference base is not a centre and contributes to no window's sum,
    /// but it still advances the finalisation frontier. Positions 1..=5 all
    /// depth 10; position 3 is `N`. The windows are over the four covered
    /// positions {1,2,4,5}; pos 3 emits nothing.
    #[test]
    fn sliding_n_position_excluded_but_advances_frontier() {
        let stream = [
            (1u32, b'G', 10u32),
            (2, b'G', 10),
            (3, b'N', 999), // huge depth, must not affect any mean
            (4, b'G', 10),
            (5, b'G', 10),
        ];
        let (windows, hist) = run_sliding(sliding_scheme(), 0, &stream);
        let positions: Vec<u32> = windows.iter().map(|w| w.pos).collect();
        assert_eq!(positions, vec![1, 2, 4, 5], "pos 3 (N) emits no window");
        for w in &windows {
            assert!(
                (w.mean_depth - 10.0).abs() < 1e-6,
                "the N position's depth 999 must not enter any mean (got {})",
                w.mean_depth,
            );
        }
        assert_eq!(hist.callable_positions, 4);
    }

    /// Windows never span chromosomes: a contig change finalises the previous
    /// contig's tail (truncated) and starts the next fresh. Two contigs, each
    /// positions 1..=3 at a distinct uniform depth.
    #[test]
    fn sliding_window_does_not_span_contigs() {
        let stream = [
            (1u32, b'G', 10u32),
            (2, b'G', 10),
            (3, b'G', 10),
            // contig 1 (chrom_id changes): a different depth
            (1, b'G', 20),
            (2, b'G', 20),
            (3, b'G', 20),
        ];
        let mut acc = SlidingWindowCoverageAccumulator::new(sliding_scheme());
        let mut out = Vec::new();
        for (i, &(pos, base, depth)) in stream.iter().enumerate() {
            let chrom = if i < 3 { 0 } else { 1 };
            acc.observe(chrom, pos, base, depth);
            while let Some(wc) = acc.pop_ready() {
                out.push(wc);
            }
        }
        let (tail, _) = acc.finish();
        out.extend(tail);
        // Six windows, three per contig; each carries only its contig's depth
        // (no cross-contig averaging).
        assert_eq!(out.len(), 6);
        for w in &out {
            let expected = if w.chrom_id == 0 { 10.0 } else { 20.0 };
            assert!(
                (w.mean_depth - expected).abs() < 1e-6,
                "chrom {} pos {} mean {} != {expected}",
                w.chrom_id,
                w.pos,
                w.mean_depth,
            );
        }
    }

    /// The accumulator is deterministic: the same stream yields identical
    /// windows and histogram.
    #[test]
    fn sliding_is_deterministic() {
        let stream: Vec<(u32, u8, u32)> = (1..=15u32).map(|p| (p, b'G', (p % 4) + 1)).collect();
        let a = run_sliding(sliding_scheme(), 0, &stream);
        let b = run_sliding(sliding_scheme(), 0, &stream);
        assert_eq!(a.0, b.0);
        assert_eq!(a.1, b.1);
    }

    /// An odd `window_bp` (`half` truncates): `window_bp = 5 → half = 2`, so the
    /// centred window is `[p-2, p+2]` (same as `window_bp = 4` here). Pins the
    /// documented `half = window_bp / 2` convention against silent drift.
    #[test]
    fn sliding_odd_window_uses_floor_half() {
        let scheme = CoverageBinScheme {
            window_bp: 5,
            ..sliding_scheme()
        };
        let stream: Vec<(u32, u8, u32)> = (1..=6u32).map(|p| (p, b'G', p)).collect();
        let (windows, _) = run_sliding(scheme, 0, &stream);
        // half = 2 → pos 3 window [1,5] = {1,2,3,4,5}, mean 3.0.
        let w3 = windows.iter().find(|w| w.pos == 3).unwrap();
        assert!((w3.mean_depth - 3.0).abs() < 1e-6, "half=2 window");
    }

    /// An empty stream yields no windows and an empty histogram.
    #[test]
    fn sliding_empty_stream_is_empty() {
        let (windows, hist) = run_sliding(sliding_scheme(), 0, &[]);
        assert!(windows.is_empty());
        assert_eq!(hist.callable_positions, 0);
        assert_eq!(hist.n_positions, 0);
    }

    /// `window_bp = 1` → `half = 0`: the window is `[p, p]` (the centre alone),
    /// the degenerate boundary where the two-pointer extents collapse. Each
    /// window's mean is exactly that position's own depth.
    #[test]
    fn sliding_window_bp_one_emits_self_only_window() {
        let scheme = CoverageBinScheme {
            window_bp: 1,
            ..sliding_scheme()
        };
        let stream: Vec<(u32, u8, u32)> = (1..=5u32).map(|p| (p, b'G', p * 10)).collect();
        let (windows, hist) = run_sliding(scheme, 0, &stream);
        assert_eq!(windows.len(), 5);
        for w in &windows {
            assert!(
                (w.mean_depth - (w.pos * 10) as f32).abs() < 1e-6,
                "pos {} window must contain only itself, got mean {}",
                w.pos,
                w.mean_depth,
            );
        }
        assert_eq!(hist.callable_positions, 5);
    }

    /// Consecutive covered positions farther apart than `window_bp`: every
    /// window is a singleton and the buffer fully drains between centres — the
    /// shrink-left / full-turnover path. A stale sum carried across the gap would
    /// contaminate a singleton mean with a far-away depth.
    #[test]
    fn sliding_large_gaps_yield_singleton_windows() {
        // window_bp 4 → half 2; positions ≫ 4 apart never share a window.
        let stream = [(10u32, b'G', 5u32), (100, b'G', 50), (1000, b'G', 500)];
        let (windows, _) = run_sliding(sliding_scheme(), 0, &stream);
        assert_eq!(windows.len(), 3);
        assert!((windows[0].mean_depth - 5.0).abs() < 1e-6);
        assert!((windows[1].mean_depth - 50.0).abs() < 1e-6);
        assert!((windows[2].mean_depth - 500.0).abs() < 1e-6);
    }

    /// A single covered position on the contig: `count == 1`, window truncated
    /// both sides — exercises the divisor-is-one path in isolation (the docs
    /// claim `count >= 1` always).
    #[test]
    fn sliding_single_covered_position_emits_one_self_window() {
        let (windows, hist) = run_sliding(sliding_scheme(), 7, &[(42u32, b'C', 9u32)]);
        assert_eq!(windows.len(), 1);
        assert_eq!(windows[0].pos, 42);
        assert!((windows[0].mean_depth - 9.0).abs() < 1e-6);
        assert!((windows[0].gc_fraction - 1.0).abs() < 1e-6);
        assert_eq!(windows[0].chrom_id, 7);
        assert_eq!(hist.callable_positions, 1);
    }

    /// `finish` returns the still-un-finalised tail (and the un-drained ready
    /// queue), so a caller that never calls `pop_ready` still receives every
    /// window — the tail cannot be silently dropped (the fix for the two-call
    /// `into_histogram` footgun).
    #[test]
    fn sliding_finish_returns_the_undrained_tail() {
        let mut acc = SlidingWindowCoverageAccumulator::new(sliding_scheme());
        for p in 1..=5u32 {
            acc.observe(0, p, b'G', 10);
            // Deliberately do NOT drain `pop_ready` during the stream.
        }
        let (windows, hist) = acc.finish();
        assert_eq!(windows.len(), 5, "every window returned, none dropped");
        assert_eq!(hist.callable_positions, 5);
        assert_eq!(hist.n_positions, 5);
    }
}
