//! S6c — the window-spill builder (Approach A, write side).
//!
//! Turns the producer's per-sample covered-position stream into the **window
//! spill**: one [`WindowSpillRecord`] per variant window, in coordinate (tile)
//! order, carrying the window's shared reference GC and each sample's mean depth
//! over it. The calibrate / write passes then join it to each locus record by
//! tile key, so the paralog score's coverage inputs are supplied without
//! carrying per-sample window data in the (per-locus, all-samples) record spill.
//!
//! The producer feeds this once-through as it consumes each chunk's positions
//! `[chunk_start, cut)` (never in the re-folding fold), one sample at a time:
//!
//! - each sample's [`WindowMeanDepthAccumulator`] tiles its positions and emits
//!   a window's mean depth as the stream leaves it;
//! - a window is **closed for the whole cohort** once every sample's open window
//!   is past it (all samples were fed to the same `cut`), so windows with tile
//!   `< min(open tile)` are final — those that contain a variant are written and
//!   evicted.
//!
//! Memory-flat: per-sample `O(1)` accumulator state + a bounded map of
//! closed-but-not-yet-flushed windows near the frontier (freed as the stream
//! advances). Accumulators + GC are reset per contig; a window straddling a
//! covered-interval gap is folded continuously (the gap positions are simply
//! uncovered and never fed), matching the `.psp` histogram's continuous tiling.

use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;

use super::spill::{SpillError, WindowSpillRecord, WindowSpillWriter};
use super::window_coverage::WindowMeanDepthAccumulator;
use super::window_gc::ReferenceWindowGc;

/// Builds the window spill from the producer's per-sample position stream.
pub(crate) struct WindowSpillBuilder<W: Write> {
    n_samples: usize,
    window_bp: u32,
    /// Per cohort sample.
    accumulators: Vec<WindowMeanDepthAccumulator>,
    writer: WindowSpillWriter<W>,
    /// Contig currently being folded; `None` before the first interval.
    chrom_id: Option<u32>,
    /// Per-window reference GC for the current contig.
    gc: Option<ReferenceWindowGc>,
    /// Per-tile per-sample mean depth for windows closed but not yet flushed
    /// (keyed by tile within the current contig).
    pending: BTreeMap<u32, Vec<Option<f32>>>,
    /// Variant tiles (current contig) not yet flushed — only these are written
    /// (the join looks up variant windows).
    variant_tiles: BTreeSet<u32>,
}

impl<W: Write> WindowSpillBuilder<W> {
    pub(crate) fn new(n_samples: usize, window_bp: u32, writer: WindowSpillWriter<W>) -> Self {
        Self {
            n_samples,
            window_bp: window_bp.max(1),
            accumulators: (0..n_samples)
                .map(|_| WindowMeanDepthAccumulator::new(window_bp))
                .collect(),
            writer,
            chrom_id: None,
            gc: None,
            pending: BTreeMap::new(),
            variant_tiles: BTreeSet::new(),
        }
    }

    /// Begin an interval on `chrom_id`. On a **contig change**, the previous
    /// contig's remaining windows are drained + flushed and the accumulators
    /// reset; consecutive intervals of the *same* contig keep folding
    /// continuously (so a window straddling a covered gap is one window). `gc` is
    /// the current contig's per-window reference GC.
    pub(crate) fn begin_interval(
        &mut self,
        chrom_id: u32,
        gc: ReferenceWindowGc,
    ) -> Result<(), SpillError> {
        if self.chrom_id != Some(chrom_id) {
            self.finish_contig()?;
            for acc in &mut self.accumulators {
                *acc = WindowMeanDepthAccumulator::new(self.window_bp);
            }
            self.pending.clear();
            self.variant_tiles.clear();
            self.chrom_id = Some(chrom_id);
            self.gc = Some(gc);
        }
        Ok(())
    }

    /// Fold one covered position for one sample. `depth = ref_obs + nonref_obs`.
    /// Records a closed window's per-sample mean depth into `pending`.
    pub(crate) fn observe(&mut self, sample_idx: usize, pos: u32, depth: u32) {
        let chrom_id = self.chrom_id.expect("observe before begin_interval");
        if let Some(window) = self.accumulators[sample_idx].observe(chrom_id, pos, depth) {
            self.pending
                .entry(window.key.1)
                .or_insert_with(|| vec![None; self.n_samples])[sample_idx] =
                Some(window.mean_depth);
        }
    }

    /// Mark the tiles of this chunk's variant positions so their windows are
    /// written when they close.
    pub(crate) fn mark_variant_positions(&mut self, positions: &[u32]) {
        for &pos in positions {
            self.variant_tiles
                .insert(pos.saturating_sub(1) / self.window_bp);
        }
    }

    /// Flush every window now closed for the whole cohort — those with tile
    /// strictly below the tile of `fed_up_to`, the **exclusive** upper bound of
    /// the positions fed so far (the chunk's `cut`). Every sample was fed the
    /// same `[chunk_start, cut)` span and every later chunk feeds positions
    /// `>= cut`, so no future position can land in a tile below
    /// `tile_of(fed_up_to)` — those windows are final regardless of whether a
    /// given (lagging or gapped) sample happened to cover them.
    ///
    /// Using the **cohort-fed frontier** (not each sample's last-observed open
    /// tile) is what keeps `pending` bounded near the frontier: a sample with a
    /// long coverage gap would otherwise pin the flush point at its stale open
    /// tile and let `pending` grow to `O(contig tiles × samples)` (S6c review).
    /// Call after each chunk's positions are fed.
    ///
    /// A sample whose open window lies **below** the frontier is final (it will
    /// get no more positions there), but its mean is still sitting in the
    /// accumulator, not in `pending`. Force-close those first, else flushing the
    /// window would drop that sample's depth (write `None` for a sample that did
    /// cover the tile).
    pub(crate) fn flush_ready(&mut self, fed_up_to: u32) -> Result<(), SpillError> {
        let frontier_tile = fed_up_to.saturating_sub(1) / self.window_bp;
        let n_samples = self.n_samples;
        for (sample_idx, acc) in self.accumulators.iter_mut().enumerate() {
            if acc
                .current_key()
                .is_some_and(|(_, tile)| tile < frontier_tile)
                && let Some(window) = acc.finish()
            {
                self.pending
                    .entry(window.key.1)
                    .or_insert_with(|| vec![None; n_samples])[sample_idx] = Some(window.mean_depth);
            }
        }
        self.flush_below(frontier_tile)
    }

    /// Write + evict every pending variant window with tile `< limit`
    /// (`u32::MAX` at contig end flushes all).
    fn flush_below(&mut self, limit: u32) -> Result<(), SpillError> {
        // Collect the tiles to flush first (can't mutate `pending` while ranging).
        let tiles: Vec<u32> = self.pending.range(..limit).map(|(&t, _)| t).collect();
        for tile in tiles {
            let depths = self.pending.remove(&tile).expect("tile just listed");
            if self.variant_tiles.remove(&tile) {
                // A variant window's reference GC is normally defined; an all-`N`
                // window (impossible for a real variant) falls back to NaN, which
                // makes the locus score non-finite → unscored → kept.
                let gc = self
                    .gc
                    .as_ref()
                    .and_then(|g| g.gc_at_tile(tile))
                    .unwrap_or(f32::NAN);
                self.writer.append(&WindowSpillRecord {
                    chrom_id: self.chrom_id.expect("contig set when flushing"),
                    tile,
                    gc,
                    depths,
                })?;
            }
        }
        Ok(())
    }

    /// Drain the accumulators (closing each sample's last window) and flush all
    /// remaining variant windows for the current contig.
    fn finish_contig(&mut self) -> Result<(), SpillError> {
        if self.chrom_id.is_none() {
            return Ok(());
        }
        for sample_idx in 0..self.n_samples {
            if let Some(window) = self.accumulators[sample_idx].finish() {
                self.pending
                    .entry(window.key.1)
                    .or_insert_with(|| vec![None; self.n_samples])[sample_idx] =
                    Some(window.mean_depth);
            }
        }
        self.flush_below(u32::MAX)
    }

    /// Drain the last contig and flush the window-spill writer, returning it.
    pub(crate) fn finish(mut self) -> Result<W, SpillError> {
        self.finish_contig()?;
        self.writer.finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::paralog_filter::spill::WindowSpillReader;
    use std::io::Cursor;

    fn gc_all(window_bp: u32, tiles: u32, gc: f32) -> ReferenceWindowGc {
        // A reference where every position is 'G' (gc = 1.0) or built to a target
        // gc via a mix; simplest: all G except we want a specific fraction. For
        // the test we just need a defined per-tile gc, so build all-'G' (gc 1.0)
        // unless gc==0 (all 'A').
        let base = if gc >= 0.5 { b'G' } else { b'A' };
        let bases = vec![base; (tiles * window_bp) as usize];
        let _ = gc;
        ReferenceWindowGc::from_bases(&bases, window_bp)
    }

    /// Two samples, three windows: the builder emits one window-spill record per
    /// variant window, in tile order, with per-sample mean depths + the shared
    /// GC, and skips a non-variant window.
    #[test]
    fn builds_variant_windows_in_tile_order() {
        let window_bp = 10u32;
        let n = 2;
        let gc = gc_all(window_bp, 4, 1.0); // gc 1.0 for every tile
        let mut b = WindowSpillBuilder::new(
            n,
            window_bp,
            WindowSpillWriter::new(Cursor::new(Vec::new())),
        );
        b.begin_interval(0, gc).unwrap();

        // Feed both samples' positions across windows 0,1,2 (pos 1..30).
        // Window 0 (pos 1..10): s0 depth 10 each at 1,2,3; s1 depth 20 at 1,2.
        // Window 1 (pos 11..20): both covered. Window 2 (pos 21..30): both.
        for s in 0..n {
            let d = if s == 0 { 10 } else { 20 };
            for pos in [1u32, 2, 3] {
                b.observe(s, pos, d);
            }
            for pos in [11u32, 12] {
                b.observe(s, pos, d);
            }
            for pos in [21u32, 22] {
                b.observe(s, pos, d);
            }
        }
        // Variants in windows 0 and 2 only (tiles 0 and 2).
        b.mark_variant_positions(&[3, 22]);
        // Fed up to pos 22 (chunk cut 23 → frontier tile 2): windows 0 and 1
        // are below the frontier and flush; window 2 stays open until finish.
        b.flush_ready(23).unwrap();
        let bytes = b.finish().unwrap().into_inner();

        let mut reader = WindowSpillReader::new(Cursor::new(bytes));
        let recs: Vec<WindowSpillRecord> = std::iter::from_fn(|| reader.next_record())
            .map(|r| r.unwrap())
            .collect();
        // Only variant windows (tiles 0 and 2), in tile order.
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].tile, 0);
        assert_eq!(recs[1].tile, 2);
        // Per-sample mean depths: s0=10, s1=20; gc 1.0.
        assert_eq!(recs[0].depths, vec![Some(10.0), Some(20.0)]);
        assert_eq!(recs[0].gc, 1.0);
        assert_eq!(recs[1].depths, vec![Some(10.0), Some(20.0)]);
    }

    /// A sample absent from a window (no covered position there) gets `None`
    /// depth for that window.
    #[test]
    fn absent_sample_gets_none_depth() {
        let window_bp = 10u32;
        let n = 2;
        let gc = gc_all(window_bp, 2, 1.0);
        let mut b = WindowSpillBuilder::new(
            n,
            window_bp,
            WindowSpillWriter::new(Cursor::new(Vec::new())),
        );
        b.begin_interval(0, gc).unwrap();
        // Sample 0 covers window 0; sample 1 covers only window 1.
        b.observe(0, 1, 15);
        b.observe(0, 2, 15);
        b.observe(1, 11, 30);
        b.mark_variant_positions(&[1]); // variant in window 0
        // Fed up to pos 11 (cut 12 → frontier tile 1): sample 0's window-0 stays
        // open (it never crossed a boundary) but is below the frontier, so it is
        // force-closed and its depth captured — the S6c-review regression: a
        // lagging sample's contribution must not be lost when the cohort
        // frontier advances past its open window.
        b.flush_ready(12).unwrap();
        let bytes = b.finish().unwrap().into_inner();
        let mut reader = WindowSpillReader::new(Cursor::new(bytes));
        let recs: Vec<WindowSpillRecord> = std::iter::from_fn(|| reader.next_record())
            .map(|r| r.unwrap())
            .collect();
        // Window 0: sample 0 present (15), sample 1 absent (None) — written by
        // flush_ready (not finish), proving the force-close captured sample 0.
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].tile, 0);
        assert_eq!(recs[0].depths, vec![Some(15.0), None]);
    }

    /// A window straddling two intervals of the same contig is folded as ONE
    /// window (accumulators are not reset between same-contig intervals).
    #[test]
    fn same_contig_intervals_fold_continuously() {
        let window_bp = 100u32;
        let n = 1;
        let gc = gc_all(window_bp, 2, 1.0);
        let mut b = WindowSpillBuilder::new(
            n,
            window_bp,
            WindowSpillWriter::new(Cursor::new(Vec::new())),
        );
        // Interval 1 of chrom 0: positions 1, 2 (window 0).
        b.begin_interval(0, gc).unwrap();
        b.observe(0, 1, 10);
        b.observe(0, 2, 20);
        b.mark_variant_positions(&[1]);
        b.flush_ready(3).unwrap(); // cut 3 → frontier tile 0: window 0 still open
        // Interval 2 of chrom 0 (a second begin_interval, SAME chrom → no reset):
        b.begin_interval(0, gc_all(window_bp, 2, 1.0)).unwrap();
        b.observe(0, 50, 30); // still window 0
        b.flush_ready(51).unwrap(); // cut 51 → frontier tile 0: still open
        let bytes = b.finish().unwrap().into_inner();
        let mut reader = WindowSpillReader::new(Cursor::new(bytes));
        let recs: Vec<WindowSpillRecord> = std::iter::from_fn(|| reader.next_record())
            .map(|r| r.unwrap())
            .collect();
        // One window (tile 0), mean over all three positions: (10+20+30)/3 = 20.
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].tile, 0);
        assert_eq!(recs[0].depths, vec![Some(20.0)]);
    }
}
