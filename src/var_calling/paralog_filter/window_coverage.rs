//! S3 — the per-sample analysis-window **mean-depth** accumulator.
//!
//! The paralog score needs, per locus per sample, that sample's 500 bp
//! analysis-window relative copy number =
//! `window_mean_depth / expected_single_copy_depth(window_gc)`. The two inputs
//! are sourced separately (settled 2026-07-01):
//!
//! - **`window_gc`** is a property of the *reference* over the window, shared by
//!   all samples — computed once from the reference in
//!   [`super::window_gc`], not here.
//! - **`window_mean_depth`** is **per-sample** and is what this module folds,
//!   straight from the caller pass's **light columns** (`depth = ref_obs +
//!   nonref_obs`, both light) — no sequence decode, no reference dependency.
//!
//! [`WindowMeanDepthAccumulator`] tiles one sample's coordinate-ordered
//! covered-position stream into `window_bp` windows and emits each window's mean
//! depth as the stream crosses a boundary.
//!
//! **`N`-handling (a deliberate approximation).** The `.psp` coverage histogram
//! the model is fit on excludes reference-`N` positions from a window's covered
//! count. This accumulator counts **every** covered position, because knowing
//! `N`-ness per position would require the reference base in the
//! throughput-critical fold (a per-position reference lookup or an
//! `O(interval)` `N` bitset — memory the design avoids). A sample's covered
//! positions are almost never reference-`N` (reads seldom place over
//! assembly-gap `N` runs), so the window mean depth is materially unchanged; the
//! `window_gc` term still excludes `N` (it is a serial reference walk where the
//! exclusion is free). The end-to-end effect is validated in T1.
//!
//! Memory-flat: one open window per sample, `O(samples)` running state, freed as
//! the coordinate-ordered stream advances — no per-window map, no genome-wide
//! structure.

/// A closed analysis window's per-sample mean depth.
///
/// Stored as `f32` to match the spill's `window_mean_depth` payload; the scorer
/// widens back to `f64` at use.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct WindowMeanDepth {
    /// `(chrom_id, tile_index)` the window covers, `tile_index =
    /// (pos − 1) / window_bp`.
    pub key: (u32, u32),
    /// Mean depth (`depth_sum / covered`) over the window's covered positions.
    pub mean_depth: f32,
}

/// The running per-window fold state for one sample.
#[derive(Debug, Clone, Copy)]
struct OpenWindow {
    key: (u32, u32),
    /// Covered positions folded into the window.
    covered: u64,
    /// Sum of per-position depth over `covered`.
    depth_sum: u64,
}

impl OpenWindow {
    /// Finalise to a [`WindowMeanDepth`], or `None` if the window has no covered
    /// position (nothing was folded — it cannot happen for an opened window,
    /// but the guard keeps the mean well-defined).
    fn finalise(self) -> Option<WindowMeanDepth> {
        if self.covered == 0 {
            return None;
        }
        Some(WindowMeanDepth {
            key: self.key,
            mean_depth: (self.depth_sum as f64 / self.covered as f64) as f32,
        })
    }
}

/// Tiles one sample's coordinate-ordered covered-position stream into
/// `window_bp` windows, emitting each window's mean depth as the stream leaves
/// it.
///
/// Feed every covered position with [`observe`] in non-decreasing
/// `(chrom_id, pos)` order; a boundary crossing returns the just-closed
/// window's mean depth. Call [`finish`] to drain the last open window.
///
/// [`observe`]: WindowMeanDepthAccumulator::observe
/// [`finish`]: WindowMeanDepthAccumulator::finish
#[derive(Debug, Clone)]
pub(crate) struct WindowMeanDepthAccumulator {
    /// Window width in bp (`>= 1`; a `0` is clamped to `1`).
    window_bp: u32,
    /// The window currently being folded; `None` before the first position.
    open: Option<OpenWindow>,
    /// `(chrom_id, pos)` of the previous [`observe`], to debug-assert the
    /// coordinate-order invariant the tiling relies on.
    last_observed: Option<(u32, u32)>,
}

impl WindowMeanDepthAccumulator {
    /// A fresh accumulator tiling into `window_bp`-wide windows.
    pub(crate) fn new(window_bp: u32) -> Self {
        Self {
            window_bp: window_bp.max(1),
            open: None,
            last_observed: None,
        }
    }

    /// The tile key for a 1-based position: `(chrom_id, (pos − 1) / window_bp)`.
    /// `pos == 0` saturates to tile `0` rather than wrapping — matching
    /// [`crate::sample_summary::coverage`]'s histogram tiling.
    fn tile_key(&self, chrom_id: u32, pos: u32) -> (u32, u32) {
        (chrom_id, pos.saturating_sub(1) / self.window_bp)
    }

    /// The window currently being folded (its key), if any.
    pub(crate) fn current_key(&self) -> Option<(u32, u32)> {
        self.open.as_ref().map(|w| w.key)
    }

    /// Fold one covered reference position. `depth` is the position's total
    /// fragment depth (`ref_obs + nonref_obs` = Σ allele observations, both
    /// light columns).
    ///
    /// Returns `Some(window)` when this position crosses into a new window (or
    /// chromosome), carrying the just-closed window's mean depth; `None` when it
    /// stays in the same window.
    ///
    /// Positions must arrive in non-decreasing `(chrom_id, pos)` order (the
    /// producer's per-sample coordinate order); this is debug-asserted.
    pub(crate) fn observe(
        &mut self,
        chrom_id: u32,
        pos: u32,
        depth: u32,
    ) -> Option<WindowMeanDepth> {
        debug_assert!(
            self.last_observed
                .is_none_or(|prev| (chrom_id, pos) >= prev),
            "window observe out of order: {:?} after {:?}",
            (chrom_id, pos),
            self.last_observed,
        );
        self.last_observed = Some((chrom_id, pos));

        let key = self.tile_key(chrom_id, pos);
        let mut closed = None;
        match self.open {
            Some(w) if w.key == key => {}
            _ => {
                if let Some(w) = self.open.take() {
                    closed = w.finalise();
                }
                self.open = Some(OpenWindow {
                    key,
                    covered: 0,
                    depth_sum: 0,
                });
            }
        }

        // PANIC-FREE: the match above sets `self.open` to `Some` on every path.
        let w = self
            .open
            .as_mut()
            .expect("open window was just set for this position");
        w.covered += 1;
        w.depth_sum = w.depth_sum.saturating_add(u64::from(depth));

        closed
    }

    /// Finalise the last open window, returning its mean depth. After this the
    /// accumulator holds no open window.
    pub(crate) fn finish(&mut self) -> Option<WindowMeanDepth> {
        self.open.take().and_then(OpenWindow::finalise)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The emitted mean depth matches an independent per-window computation over
    /// a hand-built stream (two windows).
    #[test]
    fn emits_expected_mean_depth() {
        let mut acc = WindowMeanDepthAccumulator::new(10);
        // Window (0,0): pos 1,3,5 with depths 4, 6, 2 → covered 3, sum 12,
        // mean 4.0.
        assert_eq!(acc.observe(0, 1, 4), None);
        assert_eq!(acc.observe(0, 3, 6), None);
        assert_eq!(acc.observe(0, 5, 2), None);
        // Crossing into window (0,1) closes (0,0).
        let closed = acc.observe(0, 11, 8).expect("window 0 closes");
        assert_eq!(closed.key, (0, 0));
        assert!((closed.mean_depth - 4.0).abs() < 1e-6);
        // Window (0,1): depth 8 so far. finish() drains it.
        let last = acc.finish().expect("window 1 drains");
        assert_eq!(last.key, (0, 1));
        assert_eq!(last.mean_depth, 8.0);
        assert!(acc.finish().is_none(), "nothing left after drain");
    }

    /// A chromosome change closes the open window even within the same tile
    /// index.
    #[test]
    fn chromosome_change_closes_window() {
        let mut acc = WindowMeanDepthAccumulator::new(10);
        acc.observe(0, 1, 3);
        let closed = acc.observe(1, 1, 7).expect("chrom 0 window closes");
        assert_eq!(closed.key, (0, 0));
        assert_eq!(closed.mean_depth, 3.0);
        assert_eq!(acc.current_key(), Some((1, 0)));
    }

    /// Several positions in one window average correctly; zero-depth positions
    /// count toward the denominator (a covered position with no reads is still
    /// covered).
    #[test]
    fn zero_depth_positions_count_as_covered() {
        let mut acc = WindowMeanDepthAccumulator::new(500);
        acc.observe(0, 1, 10);
        acc.observe(0, 2, 0);
        acc.observe(0, 3, 20);
        let last = acc.finish().expect("one window");
        // covered 3, sum 30, mean 10.0.
        assert_eq!(last.mean_depth, 10.0);
    }

    /// An empty accumulator drains to nothing.
    #[test]
    fn empty_accumulator_drains_to_none() {
        let mut acc = WindowMeanDepthAccumulator::new(500);
        assert!(acc.finish().is_none());
        assert_eq!(acc.current_key(), None);
    }

    /// `window_bp = 0` is clamped to 1 (every position its own window).
    #[test]
    fn zero_window_bp_is_clamped() {
        let mut acc = WindowMeanDepthAccumulator::new(0);
        acc.observe(0, 1, 2);
        let closed = acc.observe(0, 2, 4).expect("pos 1 window closes");
        assert_eq!(closed.key, (0, 0));
        assert_eq!(closed.mean_depth, 2.0);
    }

    /// The mean depth over a many-window stream matches an independent
    /// per-window Σdepth/count — the tiling boundaries line up with
    /// `(pos − 1) / window_bp`.
    #[test]
    fn matches_independent_per_window_mean() {
        use std::collections::BTreeMap;
        let window_bp = 500u32;
        let mut stream: Vec<(u32, u32, u32)> = Vec::new();
        for chrom in 0..2u32 {
            for pos in 1..=2600u32 {
                let depth = 1 + (pos % 37);
                stream.push((chrom, pos, depth));
            }
        }
        // Independent reference: (chrom, tile) → (sum, count).
        let mut reference: BTreeMap<(u32, u32), (u64, u64)> = BTreeMap::new();
        for &(c, p, d) in &stream {
            let key = (c, (p - 1) / window_bp);
            let e = reference.entry(key).or_insert((0, 0));
            e.0 += u64::from(d);
            e.1 += 1;
        }

        let mut acc = WindowMeanDepthAccumulator::new(window_bp);
        let mut emitted: BTreeMap<(u32, u32), f32> = BTreeMap::new();
        for &(c, p, d) in &stream {
            if let Some(w) = acc.observe(c, p, d) {
                emitted.insert(w.key, w.mean_depth);
            }
        }
        if let Some(w) = acc.finish() {
            emitted.insert(w.key, w.mean_depth);
        }

        assert_eq!(emitted.len(), reference.len());
        for (key, &(sum, count)) in &reference {
            let expected = (sum as f64 / count as f64) as f32;
            assert_eq!(emitted[key], expected, "window {key:?}");
        }
    }

    /// Out-of-order input violates the coordinate-order invariant and panics in
    /// debug.
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "out of order")]
    fn out_of_order_observe_panics_in_debug() {
        let mut acc = WindowMeanDepthAccumulator::new(10);
        acc.observe(0, 20, 1);
        acc.observe(0, 5, 1);
    }
}
