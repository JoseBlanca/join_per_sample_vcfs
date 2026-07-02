//! S3 — the shared per-window reference GC.
//!
//! A window's GC fraction is a property of the **reference** over that window,
//! so it is the same for every sample at a locus (settled 2026-07-01). Rather
//! than decompress each sample's `.psp` REF column, we read the GC once from the
//! reference FASTA — the same memory-flat streaming path the rest of
//! var-calling already uses ([`crate::fasta::fetcher::StreamingChromRefFetcher`],
//! a ~1 MB sliding buffer, never the whole contig).
//!
//! [`ReferenceWindowGc`] holds one contig's per-window GC fractions, indexed by
//! window tile so a locus's GC is a lookup `gc_at(pos)`. GC is the fraction of
//! **non-`N`** reference bases in the window that are `G`/`C` — the same
//! `N`-exclusion the `.psp` coverage histogram used (here it is free, because
//! this is a serial reference walk).
//!
//! Memory: `O(contig / window_bp)` — a few MB for a whole mammalian chromosome,
//! independent of sample count (the same scale as the DUST mask), so it does not
//! reintroduce the full-contig footprint the streaming fetcher avoids.

/// One contig's per-window reference GC, keyed by window tile
/// (`tile = (pos − 1) / window_bp`). Build with [`from_bases`] /
/// [`from_base_iter`]; query with [`gc_at`].
///
/// [`from_bases`]: ReferenceWindowGc::from_bases
/// [`from_base_iter`]: ReferenceWindowGc::from_base_iter
/// [`gc_at`]: ReferenceWindowGc::gc_at
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ReferenceWindowGc {
    window_bp: u32,
    /// GC fraction per window tile, densely indexed from tile 0. `None` where
    /// the window had no non-`N` reference base (all-`N`, or off the contig).
    gc: Vec<Option<f32>>,
}

impl ReferenceWindowGc {
    /// The GC fraction of the window containing 1-based `pos`, or `None` if that
    /// window had no non-`N` reference base (or `pos` is past the contig). The
    /// production join keys by tile directly ([`gc_at_tile`](Self::gc_at_tile));
    /// this position-keyed convenience is exercised only by the unit tests.
    #[cfg(test)]
    pub(crate) fn gc_at(&self, pos: u32) -> Option<f32> {
        let tile = (pos.saturating_sub(1) / self.window_bp) as usize;
        self.gc.get(tile).copied().flatten()
    }

    /// The GC fraction of window `tile` directly (the builder keys by tile).
    /// `None` if the tile is past the contig or was all-`N`.
    pub(crate) fn gc_at_tile(&self, tile: u32) -> Option<f32> {
        self.gc.get(tile as usize).copied().flatten()
    }

    /// Number of window tiles covered. (Test-only; production reads GC by tile.)
    #[cfg(test)]
    pub(crate) fn n_windows(&self) -> usize {
        self.gc.len()
    }

    /// Build from a contiguous slice of reference bases for the contig, where
    /// `bases[i]` is the base at 1-based position `i + 1`. Any case; `N` excluded
    /// from the GC denominator. Production streams the contig via
    /// [`from_base_iter`](Self::from_base_iter); this slice form is a test
    /// convenience.
    #[cfg(test)]
    pub(crate) fn from_bases(bases: &[u8], window_bp: u32) -> Self {
        Self::from_base_iter(bases.iter().copied(), window_bp)
    }

    /// Build by folding a forward iterator of reference bases starting at
    /// 1-based position 1 (the whole contig). Streaming-friendly — nothing but
    /// the running per-window counts and the output vector is held, so the S6
    /// caller can drive it from
    /// [`crate::fasta::MultiChromRefFetcher::iter_bases`] without materialising
    /// the contig. Note that `iter_bases` yields `Result<u8, _>`, so the caller
    /// unwraps each base (surfacing a mid-contig fetch error as a run failure)
    /// before feeding the bare bytes here.
    pub(crate) fn from_base_iter(bases: impl Iterator<Item = u8>, window_bp: u32) -> Self {
        let window_bp = window_bp.max(1);
        let mut gc: Vec<Option<f32>> = Vec::new();
        // Running counts for the open tile.
        let mut gc_count: u64 = 0;
        let mut non_n: u64 = 0;
        let mut pos_in_tile: u32 = 0;

        let flush = |gc_out: &mut Vec<Option<f32>>, gc_count: u64, non_n: u64| {
            gc_out.push(if non_n == 0 {
                None
            } else {
                Some((gc_count as f64 / non_n as f64) as f32)
            });
        };

        for base in bases {
            if pos_in_tile == window_bp {
                flush(&mut gc, gc_count, non_n);
                gc_count = 0;
                non_n = 0;
                pos_in_tile = 0;
            }
            if !base.eq_ignore_ascii_case(&b'N') {
                non_n += 1;
                if matches!(base.to_ascii_uppercase(), b'G' | b'C') {
                    gc_count += 1;
                }
            }
            pos_in_tile += 1;
        }
        // Flush the final (possibly partial) tile if any base was seen.
        if pos_in_tile > 0 {
            flush(&mut gc, gc_count, non_n);
        }

        Self { window_bp, gc }
    }
}

/// Whether the reference base the FASTA walk reports at a locus matches the
/// REF allele the caller decoded there — the coordinate-consistency guard.
///
/// The `.psp` REF allele's first base *is* the reference base at that position,
/// so a mismatch means the window-key / coordinate arithmetic has drifted
/// (e.g. an off-by-one against the FASTA) and the GC lookups would be silently
/// wrong. Compared case-insensitively; an `N` on either side is treated as a
/// match (soft-masked / ambiguous reference is not a coordinate error).
pub(crate) fn reference_base_matches(fasta_base: u8, decoded_ref_base: u8) -> bool {
    fasta_base.eq_ignore_ascii_case(&b'N')
        || decoded_ref_base.eq_ignore_ascii_case(&b'N')
        || fasta_base.eq_ignore_ascii_case(&decoded_ref_base)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// GC is the fraction of non-`N` bases that are G/C, per window, and the
    /// lookup maps a position to its window.
    #[test]
    fn computes_per_window_gc_fraction() {
        // window_bp 4. Window 0 (pos 1..=4): A C G T → 2/4 = 0.5.
        //             Window 1 (pos 5..=8): G G C C → 4/4 = 1.0.
        //             Window 2 (pos 9..=10): A A → 0/2 = 0.0.
        let bases = b"ACGTGGCCAA";
        let gc = ReferenceWindowGc::from_bases(bases, 4);
        assert_eq!(gc.n_windows(), 3);
        assert_eq!(gc.gc_at(1), Some(0.5));
        assert_eq!(gc.gc_at(4), Some(0.5));
        assert_eq!(gc.gc_at(5), Some(1.0));
        assert_eq!(gc.gc_at(8), Some(1.0));
        assert_eq!(gc.gc_at(9), Some(0.0));
        assert_eq!(gc.gc_at(10), Some(0.0));
    }

    /// `N` bases are excluded from the GC denominator (any case).
    #[test]
    fn excludes_n_from_denominator() {
        // window_bp 4, window 0: G N n C → non-N = {G, C} → 2/2 = 1.0.
        let bases = b"GNnC";
        let gc = ReferenceWindowGc::from_bases(bases, 4);
        assert_eq!(gc.gc_at(1), Some(1.0));
    }

    /// An all-`N` window has no non-`N` base, so its GC is `None`.
    #[test]
    fn all_n_window_is_none() {
        let bases = b"NNNN";
        let gc = ReferenceWindowGc::from_bases(bases, 4);
        assert_eq!(gc.n_windows(), 1);
        assert_eq!(gc.gc_at(1), None);
    }

    /// A position past the contig yields `None` rather than panicking.
    #[test]
    fn out_of_range_position_is_none() {
        let gc = ReferenceWindowGc::from_bases(b"GCGC", 4);
        assert_eq!(gc.gc_at(100), None);
    }

    /// The streaming `from_base_iter` matches the slice builder over a long
    /// pseudo-sequence with interspersed `N`.
    #[test]
    fn base_iter_matches_slice_builder() {
        let bases: Vec<u8> = (0..5000u32)
            .map(|i| match i % 7 {
                0 => b'N',
                1 | 2 => b'G',
                3 => b'C',
                _ => b'A',
            })
            .collect();
        let from_slice = ReferenceWindowGc::from_bases(&bases, 500);
        let from_iter = ReferenceWindowGc::from_base_iter(bases.iter().copied(), 500);
        assert_eq!(from_slice, from_iter);
    }

    /// An empty reference yields no windows.
    #[test]
    fn empty_reference_has_no_windows() {
        let gc = ReferenceWindowGc::from_bases(b"", 4);
        assert_eq!(gc.n_windows(), 0);
        assert_eq!(gc.gc_at(1), None);
    }

    /// The coordinate-consistency guard matches on equal bases (any case) and on
    /// an `N` either side, and rejects a genuine mismatch.
    #[test]
    fn reference_base_consistency_guard() {
        assert!(reference_base_matches(b'A', b'A'));
        assert!(reference_base_matches(b'a', b'A')); // case-insensitive
        assert!(reference_base_matches(b'N', b'A')); // N reference is a match
        assert!(reference_base_matches(b'G', b'n')); // N allele is a match
        assert!(!reference_base_matches(b'A', b'C')); // genuine mismatch
    }
}
