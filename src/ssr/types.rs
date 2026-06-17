//! Shared SSR domain types — the spine every stage of the SSR caller passes.
//!
//! This module is built incrementally: it currently carries only the types
//! Stage 0 (the catalog builder) produces — [`Motif`] and [`Locus`]. The allele
//! representation (the on-/off-ladder hybrid used by Stages 1–2) is added when
//! those stages are built. The full design is in
//! `doc/devel/architecture/ssr_shared_types.md`.

use std::fmt;
use std::ops::Range;

/// SSR scope: a repeat unit (period) is between 1 and this many bases.
///
/// The bound is the microsatellite period ceiling (period ≤ 6) that scopes the
/// whole caller — see `doc/devel/specs/ssr_genotyping.md` (glossary: *SSR*,
/// *period*).
pub(crate) const MAX_MOTIF_LEN: usize = 6;

/// A tandem-repeat unit — the repeat's period, 1..=[`MAX_MOTIF_LEN`] bases.
///
/// The bytes are stored **verbatim**: the reference-strand, phase-faithful unit
/// exactly as it tiles the locus (e.g. `CAG`), *not* canonicalized. Rotating to
/// a canonical form (`CAG` → `AGC`) would break tiling, and reconstruction reads
/// phase-correct bytes off the reference anyway; the canonical *class* used for
/// stutter pooling is derived on demand downstream, never stored here.
///
/// Inline and `Copy`: a fixed 6-byte buffer plus a length, never heap-allocated.
/// Unused tail bytes are zero, so the derived `Eq`/`Hash` compare only the live
/// prefix (`0` is not a valid base, so it cannot collide with one).
///
/// INVARIANT (relied on by the derived `Eq`/`Hash`): every constructor MUST
/// zero-initialize the unused tail of `buf`. `new` is currently the only one;
/// any future constructor must uphold this or the derived impls will treat
/// equal motifs as distinct.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Motif {
    buf: [u8; MAX_MOTIF_LEN],
    len: u8,
}

/// A motif's bytes were not a valid SSR period.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub(crate) enum MotifError {
    /// Length is `0` or above [`MAX_MOTIF_LEN`] — outside the SSR period range.
    #[error("motif length {len} is outside the SSR range 1..={MAX_MOTIF_LEN}")]
    BadLength { len: usize },
}

impl Motif {
    /// Build a motif from its bytes, validating the SSR period range.
    ///
    /// The bytes are taken verbatim (no canonicalization); the caller is
    /// responsible for supplying the phase-faithful, reference-strand unit.
    pub fn new(bytes: &[u8]) -> Result<Self, MotifError> {
        let len = bytes.len();
        if len == 0 || len > MAX_MOTIF_LEN {
            return Err(MotifError::BadLength { len });
        }
        let mut buf = [0u8; MAX_MOTIF_LEN];
        buf[..len].copy_from_slice(bytes);
        Ok(Self {
            buf,
            len: len as u8,
        })
    }

    /// The motif bytes.
    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        &self.buf[..self.len as usize]
    }

    /// The period (motif length, in bases).
    #[inline]
    pub fn period(&self) -> usize {
        self.len as usize
    }
}

impl fmt::Debug for Motif {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Bases are ASCII; render as text for readable test output / logs,
        // falling back to bytes if a motif ever held non-UTF-8.
        match std::str::from_utf8(self.as_bytes()) {
            Ok(s) => write!(f, "Motif({s:?})"),
            Err(_) => write!(f, "Motif({:?})", self.as_bytes()),
        }
    }
}

/// A [`Locus`] could not be built because its inputs broke a documented invariant.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub(crate) enum LocusError {
    /// Tract coordinates are not ordered `ref_bytes_start <= start <= end`.
    #[error(
        "locus coordinates out of order: require ref_bytes_start ({ref_bytes_start}) \
         <= start ({start}) <= end ({end})"
    )]
    BadCoordinates {
        ref_bytes_start: u32,
        start: u32,
        end: u32,
    },
    /// The tract `end` runs past the embedded reference span.
    #[error(
        "tract end ({end}) exceeds embedded reference span \
         [{ref_bytes_start}, {ref_bytes_end})"
    )]
    TractBeyondRefBytes {
        end: u32,
        ref_bytes_start: u32,
        ref_bytes_end: u64,
    },
    /// `purity_fraction` is not a finite value in `[0.0, 1.0]`.
    #[error("purity fraction {purity_fraction} is not finite in [0.0, 1.0]")]
    BadPurity { purity_fraction: f32 },
}

/// One catalog locus — a single SSR, carrying its own local reference bases.
///
/// Embedding the reference (`ref_bytes`) is what lets Stages 1–2 reconstruct any
/// allele without ever opening the reference FASTA (the catalog is the only
/// reference-bearing artifact). `chrom` is the contig **name**: the catalog is a
/// name-keyed text artifact, matching the project's name-based contig model
/// (`fasta::ContigEntry`).
///
/// Coordinates are 0-based half-open. `ref_bytes` spans `[ref_bytes_start,
/// ref_bytes_start + ref_bytes.len())` — the tract `[start, end)` plus a flank
/// margin each side (clamped at contig ends), upper-cased. The fields are private
/// and the invariant `ref_bytes_start <= start <= end <= ref_bytes_start +
/// ref_bytes.len()` (plus `purity_fraction` finite ∈ `[0.0, 1.0]`) is enforced by
/// [`Locus::new`], so the accessors below are infallible by construction.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Locus {
    /// Contig name.
    chrom: Box<str>,
    /// Tract start (0-based).
    start: u32,
    /// Tract end (exclusive).
    end: u32,
    /// The repeat unit (verbatim, phase-faithful; `period = motif.period()`).
    motif: Motif,
    /// Fraction of the tract matching a perfect motif tiling, in `[0.0, 1.0]`
    /// (`1.0` = perfect). A *degree*, not a flag — see [`Self::is_perfect`].
    purity_fraction: f32,
    /// Embedded reference bases: the tract plus a flank margin each side,
    /// upper-cased, clamped at contig ends.
    ref_bytes: Box<[u8]>,
    /// Genomic coordinate of `ref_bytes[0]`.
    ref_bytes_start: u32,
}

impl Locus {
    /// Build a catalog locus, validating its coordinate and purity invariants.
    ///
    /// The inputs come from external tooling (TRF output + reference bytes), so
    /// the ordering and bounds are checked here rather than assumed: a bad
    /// catalog row becomes a typed [`LocusError`] instead of a release-mode
    /// underflow or out-of-bounds slice in the accessors below.
    pub fn new(
        chrom: Box<str>,
        start: u32,
        end: u32,
        motif: Motif,
        purity_fraction: f32,
        ref_bytes: Box<[u8]>,
        ref_bytes_start: u32,
    ) -> Result<Self, LocusError> {
        if !(ref_bytes_start <= start && start <= end) {
            return Err(LocusError::BadCoordinates {
                ref_bytes_start,
                start,
                end,
            });
        }
        let ref_bytes_end = ref_bytes_start as u64 + ref_bytes.len() as u64;
        if end as u64 > ref_bytes_end {
            return Err(LocusError::TractBeyondRefBytes {
                end,
                ref_bytes_start,
                ref_bytes_end,
            });
        }
        if !purity_fraction.is_finite() || !(0.0..=1.0).contains(&purity_fraction) {
            return Err(LocusError::BadPurity { purity_fraction });
        }
        Ok(Self {
            chrom,
            start,
            end,
            motif,
            purity_fraction,
            ref_bytes,
            ref_bytes_start,
        })
    }

    /// Contig name.
    #[inline]
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    /// Tract start (0-based, inclusive).
    #[inline]
    pub fn start(&self) -> u32 {
        self.start
    }

    /// Tract end (0-based, exclusive).
    #[inline]
    pub fn end(&self) -> u32 {
        self.end
    }

    /// The repeat unit (verbatim, phase-faithful).
    #[inline]
    pub fn motif(&self) -> Motif {
        self.motif
    }

    /// Fraction of the tract matching a perfect motif tiling, in `[0.0, 1.0]`.
    #[inline]
    pub fn purity_fraction(&self) -> f32 {
        self.purity_fraction
    }

    /// Embedded reference bases (tract plus flank margins, upper-cased).
    #[inline]
    pub fn ref_bytes(&self) -> &[u8] {
        &self.ref_bytes
    }

    /// Genomic coordinate of `ref_bytes()[0]`.
    #[inline]
    pub fn ref_bytes_start(&self) -> u32 {
        self.ref_bytes_start
    }

    /// The period (motif length, in bases).
    #[inline]
    pub fn period(&self) -> usize {
        self.motif.period()
    }

    /// Whether the tract is a perfect (uninterrupted) motif tiling.
    ///
    /// Derived from `purity_fraction`, never stored: the perfect/imperfect split
    /// is a threshold on the degree.
    #[inline]
    pub fn is_perfect(&self) -> bool {
        self.purity_fraction == 1.0
    }

    /// The tract's byte range within `ref_bytes`.
    ///
    /// Both bounds hold by construction (enforced in [`Locus::new`]); the
    /// `debug_assert` documents the lower-bound precondition the subtraction
    /// relies on.
    #[inline]
    fn tract_range(&self) -> Range<usize> {
        debug_assert!(
            self.ref_bytes_start <= self.start,
            "ref_bytes precedes tract start"
        );
        let offset = (self.start - self.ref_bytes_start) as usize;
        offset..offset + (self.end - self.start) as usize
    }

    /// The reference tract bytes — the REF allele's sequence.
    pub fn ref_tract(&self) -> &[u8] {
        &self.ref_bytes[self.tract_range()]
    }

    /// The left flank: embedded reference bases before the tract.
    pub fn left_flank(&self) -> &[u8] {
        &self.ref_bytes[..self.tract_range().start]
    }

    /// The right flank: embedded reference bases after the tract.
    pub fn right_flank(&self) -> &[u8] {
        &self.ref_bytes[self.tract_range().end..]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn motif_roundtrips_bytes_and_period() {
        let m = Motif::new(b"CAG").unwrap();
        assert_eq!(m.as_bytes(), b"CAG");
        assert_eq!(m.period(), 3);
    }

    #[test]
    fn motif_accepts_full_period_range() {
        for len in 1..=MAX_MOTIF_LEN {
            let bytes = vec![b'A'; len];
            assert_eq!(Motif::new(&bytes).unwrap().period(), len);
        }
    }

    #[test]
    fn motif_rejects_out_of_range_lengths() {
        assert_eq!(Motif::new(b""), Err(MotifError::BadLength { len: 0 }));
        assert_eq!(
            Motif::new(b"ACGTACG"),
            Err(MotifError::BadLength { len: 7 })
        );
    }

    #[test]
    fn motif_equality_ignores_unused_tail() {
        // Two motifs of different length must differ; same bytes must match
        // regardless of the zeroed buffer tail.
        assert_eq!(Motif::new(b"AT").unwrap(), Motif::new(b"AT").unwrap());
        assert_ne!(Motif::new(b"AT").unwrap(), Motif::new(b"ATA").unwrap());
        assert_ne!(Motif::new(b"AT").unwrap(), Motif::new(b"AC").unwrap());
    }

    #[test]
    fn motif_debug_renders_ascii_then_falls_back_on_non_utf8() {
        // ASCII bases render as text.
        assert_eq!(
            format!("{:?}", Motif::new(b"CAG").unwrap()),
            "Motif(\"CAG\")"
        );
        // A non-UTF-8 motif (constructible only via raw bytes) falls back to the
        // byte rendering rather than panicking.
        let non_utf8 = Motif::new(&[0xFF, 0xFE]).unwrap();
        assert_eq!(format!("{non_utf8:?}"), "Motif([255, 254])");
    }

    #[test]
    fn motif_is_hashable_and_dedups_in_a_set() {
        use std::collections::HashSet;
        let mut set = HashSet::new();
        set.insert(Motif::new(b"AT").unwrap());
        set.insert(Motif::new(b"AT").unwrap()); // duplicate collapses
        set.insert(Motif::new(b"ATA").unwrap()); // distinct length stays
        assert_eq!(set.len(), 2);
        assert!(set.contains(&Motif::new(b"AT").unwrap()));
        assert!(set.contains(&Motif::new(b"ATA").unwrap()));
    }

    /// A perfect dinucleotide locus: ref_bytes = 3bp left flank + (CA)x3 tract +
    /// 3bp right flank, tract at genomic [13, 19).
    fn sample_locus() -> Locus {
        Locus::new(
            "chr1".into(),
            13,
            19,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGCACACATTT").into(), // GGG | CACACA | TTT
            10,
        )
        .unwrap()
    }

    #[test]
    fn locus_ref_tract_and_flanks_slice_correctly() {
        let l = sample_locus();
        assert_eq!(l.ref_tract(), b"CACACA");
        assert_eq!(l.left_flank(), b"GGG");
        assert_eq!(l.right_flank(), b"TTT");
        assert_eq!(l.period(), 2);
    }

    #[test]
    fn locus_is_perfect_thresholds_on_purity() {
        assert!(sample_locus().is_perfect());
        let imperfect = Locus::new(
            "chr1".into(),
            13,
            19,
            Motif::new(b"CA").unwrap(),
            0.9,
            (*b"GGGCACACATTT").into(),
            10,
        )
        .unwrap();
        assert!(!imperfect.is_perfect());
    }

    #[test]
    fn locus_new_rejects_ref_bytes_start_after_start() {
        // ref_bytes_start (14) > start (13): in release this would underflow
        // `start - ref_bytes_start` in the accessors. Must be a typed error.
        let err = Locus::new(
            "chr1".into(),
            13,
            19,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGCACACATTT").into(),
            14,
        )
        .unwrap_err();
        assert_eq!(
            err,
            LocusError::BadCoordinates {
                ref_bytes_start: 14,
                start: 13,
                end: 19,
            }
        );
    }

    #[test]
    fn locus_new_rejects_start_after_end() {
        let err = Locus::new(
            "chr1".into(),
            20,
            19,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGCACACATTT").into(),
            10,
        )
        .unwrap_err();
        assert_eq!(
            err,
            LocusError::BadCoordinates {
                ref_bytes_start: 10,
                start: 20,
                end: 19,
            }
        );
    }

    #[test]
    fn locus_new_rejects_end_beyond_ref_bytes() {
        // ref_bytes spans [10, 22); end = 23 would slice past the buffer end.
        let err = Locus::new(
            "chr1".into(),
            13,
            23,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGCACACATTT").into(),
            10,
        )
        .unwrap_err();
        assert_eq!(
            err,
            LocusError::TractBeyondRefBytes {
                end: 23,
                ref_bytes_start: 10,
                ref_bytes_end: 22,
            }
        );
    }

    #[test]
    fn locus_new_rejects_non_finite_or_out_of_range_purity() {
        for bad in [f32::NAN, 1.5, -0.1, f32::INFINITY] {
            let err = Locus::new(
                "chr1".into(),
                13,
                19,
                Motif::new(b"CA").unwrap(),
                bad,
                (*b"GGGCACACATTT").into(),
                10,
            )
            .unwrap_err();
            // NaN != NaN, so match the variant rather than compare by value.
            assert!(
                matches!(err, LocusError::BadPurity { .. }),
                "expected BadPurity for purity {bad}, got {err:?}"
            );
        }
    }

    #[test]
    fn locus_new_accepts_boundary_purity() {
        for ok in [0.0, 1.0] {
            assert!(
                Locus::new(
                    "chr1".into(),
                    13,
                    19,
                    Motif::new(b"CA").unwrap(),
                    ok,
                    (*b"GGGCACACATTT").into(),
                    10,
                )
                .is_ok(),
                "purity {ok} should be accepted"
            );
        }
    }

    #[test]
    fn locus_ref_tract_returns_tract_when_left_flank_empty() {
        // ref_bytes_start == start (tract clamped at the contig start): the left
        // flank is empty and the tract is the whole leading slice.
        let l = Locus::new(
            "chr1".into(),
            10,
            16,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"CACACATTT").into(), // CACACA | TTT, no left flank
            10,
        )
        .unwrap();
        assert_eq!(l.left_flank(), b"");
        assert_eq!(l.ref_tract(), b"CACACA");
        assert_eq!(l.right_flank(), b"TTT");
    }

    #[test]
    fn locus_right_flank_empty_when_clamped_at_contig_end() {
        // end == ref_bytes_start + ref_bytes.len() (tract clamped at the contig
        // end): the right flank is empty and slicing must not panic.
        let l = Locus::new(
            "chr1".into(),
            13,
            19,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGCACACA").into(), // GGG | CACACA, no right flank
            10,
        )
        .unwrap();
        assert_eq!(l.left_flank(), b"GGG");
        assert_eq!(l.ref_tract(), b"CACACA");
        assert_eq!(l.right_flank(), b"");
    }

    #[test]
    fn locus_is_perfect_false_on_nan_purity() {
        // Documents the `== 1.0` behavior: a NaN purity is rejected at
        // construction, but if one ever reached `is_perfect` it returns false.
        assert!(
            Locus::new(
                "chr1".into(),
                13,
                19,
                Motif::new(b"CA").unwrap(),
                f32::NAN,
                (*b"GGGCACACATTT").into(),
                10,
            )
            .is_err()
        );
    }
}
