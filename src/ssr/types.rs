//! Shared SSR domain types — the spine every stage of the SSR caller passes.
//!
//! This module is built incrementally: it currently carries only the types
//! Stage 0 (the catalog builder) produces — [`Motif`] and [`Locus`]. The allele
//! representation (the on-/off-ladder hybrid used by Stages 1–2) is added when
//! those stages are built. The full design is in
//! `doc/devel/architecture/ssr_shared_types.md`.

use std::fmt;

/// SSR scope: a repeat unit (period) is between 1 and this many bases.
pub const MAX_MOTIF_LEN: usize = 6;

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
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Motif {
    buf: [u8; MAX_MOTIF_LEN],
    len: u8,
}

/// A motif's bytes were not a valid SSR period.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum MotifError {
    /// Length is `0` or above [`MAX_MOTIF_LEN`] — outside the SSR period range.
    #[error("motif length {0} is outside the SSR range 1..={MAX_MOTIF_LEN}")]
    BadLength(usize),
}

impl Motif {
    /// Build a motif from its bytes, validating the SSR period range.
    ///
    /// The bytes are taken verbatim (no canonicalization); the caller is
    /// responsible for supplying the phase-faithful, reference-strand unit.
    pub fn new(bytes: &[u8]) -> Result<Self, MotifError> {
        let len = bytes.len();
        if len == 0 || len > MAX_MOTIF_LEN {
            return Err(MotifError::BadLength(len));
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
/// margin each side (clamped at contig ends), upper-cased. Invariants, upheld at
/// construction: `ref_bytes_start <= start <= end <= ref_bytes_start +
/// ref_bytes.len()`.
#[derive(Clone, Debug, PartialEq)]
pub struct Locus {
    /// Contig name.
    pub chrom: Box<str>,
    /// Tract start (0-based).
    pub start: u32,
    /// Tract end (exclusive).
    pub end: u32,
    /// The repeat unit (verbatim, phase-faithful; `period = motif.period()`).
    pub motif: Motif,
    /// Fraction of the tract matching a perfect motif tiling, in `[0.0, 1.0]`
    /// (`1.0` = perfect). A *degree*, not a flag — see [`Self::is_perfect`].
    pub purity_fraction: f32,
    /// Embedded reference bases: the tract plus a flank margin each side,
    /// upper-cased, clamped at contig ends.
    pub ref_bytes: Box<[u8]>,
    /// Genomic coordinate of `ref_bytes[0]`.
    pub ref_bytes_start: u32,
}

impl Locus {
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

    /// Offset of the tract within `ref_bytes`.
    #[inline]
    fn tract_offset(&self) -> usize {
        debug_assert!(
            self.ref_bytes_start <= self.start,
            "ref_bytes precedes tract start"
        );
        (self.start - self.ref_bytes_start) as usize
    }

    /// The reference tract bytes — the REF allele's sequence.
    pub fn ref_tract(&self) -> &[u8] {
        let start = self.tract_offset();
        let end = start + (self.end - self.start) as usize;
        &self.ref_bytes[start..end]
    }

    /// The left flank: embedded reference bases before the tract.
    pub fn left_flank(&self) -> &[u8] {
        &self.ref_bytes[..self.tract_offset()]
    }

    /// The right flank: embedded reference bases after the tract.
    pub fn right_flank(&self) -> &[u8] {
        let end = self.tract_offset() + (self.end - self.start) as usize;
        &self.ref_bytes[end..]
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
        assert_eq!(Motif::new(b""), Err(MotifError::BadLength(0)));
        assert_eq!(Motif::new(b"ACGTACG"), Err(MotifError::BadLength(7)));
    }

    #[test]
    fn motif_equality_ignores_unused_tail() {
        // Two motifs of different length must differ; same bytes must match
        // regardless of the zeroed buffer tail.
        assert_eq!(Motif::new(b"AT").unwrap(), Motif::new(b"AT").unwrap());
        assert_ne!(Motif::new(b"AT").unwrap(), Motif::new(b"ATA").unwrap());
        assert_ne!(Motif::new(b"AT").unwrap(), Motif::new(b"AC").unwrap());
    }

    /// A perfect dinucleotide locus: ref_bytes = 3bp left flank + (CA)x3 tract +
    /// 3bp right flank, tract at genomic [13, 19).
    fn sample_locus() -> Locus {
        Locus {
            chrom: "chr1".into(),
            start: 13,
            end: 19,
            motif: Motif::new(b"CA").unwrap(),
            purity_fraction: 1.0,
            ref_bytes: (*b"GGGCACACATTT").into(), // GGG | CACACA | TTT
            ref_bytes_start: 10,
        }
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
        let mut l = sample_locus();
        assert!(l.is_perfect());
        l.purity_fraction = 0.9;
        assert!(!l.is_perfect());
    }
}
