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

use super::CoverageByGcHistogram;

/// Bin scheme for the coverage-by-GC histogram. These are the tuning
/// knobs (architecture Premise 2/3); the `pileup` CLI supplies them
/// (C1) and validates them at that boundary. Every field must be
/// positive (`depth_bin_width` finite and `> 0`); [`CoverageByGcAccumulator::new`]
/// asserts this in both debug and release rather than silently using a
/// scheme the caller never chose.
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
    n_tiles: u64,
    n_skipped_tiles: u64,
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
        assert!(scheme.window_bp >= 1, "window_bp must be >= 1");
        assert!(scheme.gc_bins >= 1, "gc_bins must be >= 1");
        assert!(scheme.depth_bins >= 1, "depth_bins must be >= 1");
        assert!(
            scheme.depth_bin_width.is_finite() && scheme.depth_bin_width > 0.0,
            "depth_bin_width must be finite and > 0, got {}",
            scheme.depth_bin_width,
        );
        Self {
            counts: vec![0; scheme.n_cells()],
            scheme,
            open: None,
            last_observed: None,
            n_tiles: 0,
            n_skipped_tiles: 0,
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
            n_tiles: self.n_tiles,
            n_skipped_tiles: self.n_skipped_tiles,
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
        let gc_frac = tile.gc as f64 / tile.covered as f64;
        let mean_depth = tile.depth_sum as f64 / tile.covered as f64;
        let cell = self.cell_index(gc_frac, mean_depth);
        self.counts[cell] += 1;
        self.n_tiles += 1;
    }

    /// Map a tile's `(GC fraction in [0, 1], mean depth >= 0)` to a
    /// row-major cell index. GC saturates to the last GC bin at 1.0; depth
    /// at or above `depth_bins * depth_bin_width` lands in the overflow
    /// column. `mean_depth` is always `>= 0` (it is a sum of non-negative
    /// counts over a positive divisor), so `0` maps to depth bin `0`; the
    /// `as usize` cast also saturates a hypothetical negative to `0`.
    fn cell_index(&self, gc_frac: f64, mean_depth: f64) -> usize {
        let gc_bins = self.scheme.gc_bins as usize;
        let gc_bin = ((gc_frac * gc_bins as f64) as usize).min(gc_bins - 1);

        let depth_bins = self.scheme.depth_bins as usize;
        let depth_bin = ((mean_depth / self.scheme.depth_bin_width) as usize).min(depth_bins);
        gc_bin * self.scheme.depth_cols() + depth_bin
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
        assert_eq!(h.n_tiles, 1);
        assert_eq!(h.n_skipped_tiles, 0);
        // GC = 0.5 -> gc_bin = floor(0.5 * 2) = 1. depth = 2.0, width 1 ->
        // depth_bin = 2. cell = 1 * (4 + 1) + 2 = 7.
        let total: u32 = h.counts.iter().sum();
        assert_eq!(total, 1);
        assert_eq!(h.counts[7], 1);
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
        assert_eq!(h.n_tiles, 1);
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
        assert_eq!(h.n_tiles, 0);
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
        assert_eq!(h.n_tiles, 3);
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
                n_variant_sites: 0,
                min_depth: 1,
                het_vaf_lo: 0.3,
                het_vaf_hi: 0.7,
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
        assert_eq!(h.n_tiles, 0);
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
        assert_eq!(h.n_tiles, 1);
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
}
