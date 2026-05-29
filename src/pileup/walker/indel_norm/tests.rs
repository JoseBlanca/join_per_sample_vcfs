//! Unit tests for indel left-alignment (the pure CIGAR-rewrite routine).
//!
//! Each case fixes a reference window, a read sequence, and the CIGAR a
//! short-read aligner *could* have emitted (indel placed at an arbitrary
//! offset in a repeat), then asserts the normalized CIGAR puts the indel
//! at its leftmost equivalent position — the `bcftools norm` form.

use super::*;
use crate::pileup::walker::CigarOp::{self, *};

/// Left-align with `read_start = 0` (reference window begins at the read's
/// first aligned base — what the walker fetches per read).
fn la(cigar: &[CigarOp], ref_bases: &str, read: &str) -> LeftAlignResult {
    left_align_cigar(cigar, ref_bases.as_bytes(), read.as_bytes(), 0)
}

#[test]
fn passthrough_when_no_indel() {
    // A pure-match read is returned unchanged, no allocation surprises.
    let r = la(&[Match(5)], "ACGTA", "ACGTA");
    assert_eq!(r.cigar, vec![Match(5)]);
    assert_eq!(r.leading_deletion_bases_removed, 0);
}

#[test]
fn homopolymer_deletion_shifts_left() {
    // ref GAAAAT, read GAAAT (one A deleted). Aligner placed the deletion
    // at the rightmost A (4M1D1M); canonical form is the leftmost A.
    let r = la(&[Match(4), Deletion(1), Match(1)], "GAAAAT", "GAAAT");
    assert_eq!(r.cigar, vec![Match(1), Deletion(1), Match(4)]);
    assert_eq!(r.leading_deletion_bases_removed, 0);
}

#[test]
fn homopolymer_insertion_shifts_left() {
    // ref GAAAT, read GAAAAT (one A inserted). Aligner placed the
    // insertion at the right (4M1I1M); canonical form is the leftmost.
    let r = la(&[Match(4), Insertion(1), Match(1)], "GAAAT", "GAAAAT");
    assert_eq!(r.cigar, vec![Match(1), Insertion(1), Match(4)]);
    assert_eq!(r.leading_deletion_bases_removed, 0);
}

#[test]
fn dinucleotide_ssr_deletion_shifts_left() {
    // ref CATATAT = C(AT)(AT)(AT), read CATAT = C(AT)(AT) — one AT deleted.
    // Rightmost placement 5M2D; canonical is the leftmost AT unit.
    let r = la(&[Match(5), Deletion(2)], "CATATAT", "CATAT");
    assert_eq!(r.cigar, vec![Match(1), Deletion(2), Match(4)]);
    assert_eq!(r.leading_deletion_bases_removed, 0);
}

#[test]
fn already_canonical_unique_context_unchanged() {
    // Deletion of a non-repeated base cannot move.
    let r = la(&[Match(2), Deletion(1), Match(2)], "GACTT", "GATT");
    assert_eq!(r.cigar, vec![Match(2), Deletion(1), Match(2)]);
    assert_eq!(r.leading_deletion_bases_removed, 0);
}

#[test]
fn read_side_mismatch_blocks_overshift() {
    // ref GAAAAT but the read reads C where the second base aligns
    // (read GCAAT, a deletion of one A). Reference-only normalization
    // would roll to the first A; the read's C must stop it earlier, so
    // the indel is left-aligned only as far as the read still supports.
    let r = la(&[Match(4), Deletion(1), Match(1)], "GAAAAT", "GCAAT");
    // The deletion may not cross read position 1 (the C).
    assert_eq!(r.leading_deletion_bases_removed, 0);
    // Reconstruct the deletion's reference offset from the CIGAR.
    let mut ref_off = 0u32;
    let mut del_ref_pos = None;
    for op in &r.cigar {
        if let Deletion(_) = op {
            del_ref_pos = Some(ref_off);
            break;
        }
        ref_off += match op {
            Match(n) | Deletion(n) | SeqMatch(n) | SeqMismatch(n) | Skip(n) => *n,
            _ => 0,
        };
    }
    // 0-based ref offset of the deleted base must be >= 2 (cannot reach
    // the first A at offset 1, which the read contradicts with C).
    assert!(
        del_ref_pos.is_some_and(|p| p >= 2),
        "deletion over-shifted past the read mismatch: {:?}",
        r.cigar
    );
}

#[test]
fn deletion_shifted_to_read_start_becomes_leading() {
    // ref AAAAT, read AAAT — deletion rolls all the way to the read start
    // and is stripped, bumping the read's alignment_start by 1.
    let r = la(&[Match(3), Deletion(1), Match(1)], "AAAAT", "AAAT");
    assert_eq!(r.cigar, vec![Match(4)]);
    assert_eq!(r.leading_deletion_bases_removed, 1);
}

#[test]
fn insertion_seq_rotates_implicitly_via_read_offset() {
    // The routine moves only CIGAR boundaries; verify the inserted bases
    // a cursor would read from the shifted I op are the rotated/canonical
    // ones. ref GCACACT, read GCACACACT (one CA inserted in the CA SSR).
    // Canonical insertion is the leftmost CA unit.
    let r = la(&[Match(6), Insertion(2), Match(1)], "GCACACT", "GCACACACT");
    // Leftmost CA: G | CA(ins) | CACAC T  -> 1M 2I 6M
    assert_eq!(r.cigar, vec![Match(1), Insertion(2), Match(6)]);
}

#[test]
fn no_change_preserves_match_op_count() {
    // Indel already leftmost: unchanged.
    let r = la(&[Match(1), Deletion(1), Match(4)], "GAAAAT", "GAAAT");
    assert_eq!(r.cigar, vec![Match(1), Deletion(1), Match(4)]);
}

#[test]
fn soft_clips_are_preserved_and_not_crossed() {
    // Leading soft clip: read TT then GAAAT aligned. Deletion in the
    // A-run left-aligns within the aligned region but not into the clip.
    // ref GAAAAT, read (with 2S) = TT GAAAT.
    let r = left_align_cigar(
        &[SoftClip(2), Match(4), Deletion(1), Match(1)],
        "GAAAAT".as_bytes(),
        "TTGAAAT".as_bytes(),
        0,
    );
    assert_eq!(r.cigar, vec![SoftClip(2), Match(1), Deletion(1), Match(4)]);
    assert_eq!(r.leading_deletion_bases_removed, 0);
}
