//! Unit tests for indel left-alignment (the pure CIGAR-rewrite routine).
//!
//! Each case fixes a reference window, a read sequence, and the CIGAR a
//! short-read aligner *could* have emitted (indel placed at an arbitrary
//! offset in a repeat), then asserts the normalized CIGAR puts the indel
//! at its leftmost equivalent position — the `bcftools norm` form.

use super::*;
use crate::pileup::walker::CigarOp::{self, *};

/// Left-align with `read_start = 0` (reference window begins at the read's
/// first aligned base) and end-deletion stripping on (the canonical form).
fn la(cigar: &[CigarOp], ref_bases: &str, read: &str) -> LeftAlignResult {
    left_align_cigar(cigar, ref_bases.as_bytes(), read.as_bytes(), 0, true)
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
        true,
    );
    assert_eq!(r.cigar, vec![SoftClip(2), Match(1), Deletion(1), Match(4)]);
    assert_eq!(r.leading_deletion_bases_removed, 0);
}

#[test]
fn deletion_to_read_start_kept_when_not_removing_ends() {
    // Same input as `deletion_shifted_to_read_start_becomes_leading`, but
    // with end-deletion stripping off (the per-sample read-prep mode): the
    // leading deletion stays in the CIGAR and alignment_start is *not*
    // bumped, preserving the coordinate-sort invariant. The cursor rejects
    // this first-op deletion downstream.
    let r = left_align_cigar(
        &[Match(3), Deletion(1), Match(1)],
        "AAAAT".as_bytes(),
        "AAAT".as_bytes(),
        0,
        false,
    );
    assert_eq!(r.cigar, vec![Deletion(1), Match(4)]);
    assert_eq!(r.leading_deletion_bases_removed, 0);
}

#[test]
fn malformed_read_consumption_returns_input_unchanged() {
    // Untrusted-input guard: a CIGAR whose read-consuming ops (4+1+1 = 6)
    // disagree with seq.len() (5) is left untouched rather than producing a
    // wrong-length CIGAR. The downstream walker length check then rejects
    // the read, exactly as without normalization.
    let cigar = [Match(4), Insertion(1), Match(1)];
    let r = left_align_cigar(&cigar, "AAAAAA".as_bytes(), "AAAAA".as_bytes(), 0, false);
    assert_eq!(r.cigar, cigar.to_vec());
    assert_eq!(r.leading_deletion_bases_removed, 0);
}

#[test]
fn two_reads_with_same_deletion_converge() {
    // The core recall fix: the same biological deletion of one A from a
    // homopolymer, placed by the aligner at opposite ends (5M1D1M vs
    // 1M1D5M), must canonicalise to the *same* CIGAR so the cohort merge
    // buckets the two reads onto one allele.
    let right = la(&[Match(5), Deletion(1), Match(1)], "GAAAAAC", "GAAAAC");
    let left = la(&[Match(1), Deletion(1), Match(5)], "GAAAAAC", "GAAAAC");
    assert_eq!(right.cigar, left.cigar);
    assert_eq!(right.cigar, vec![Match(1), Deletion(1), Match(5)]);
}

#[test]
fn multi_base_ssr_insertion_shifts_to_leftmost_unit() {
    // ref GCACAC (G + CA + CA), read GCACACAC (one extra CA unit). The
    // inserted CA left-aligns to the first unit.
    let r = la(&[Match(2), Insertion(2), Match(4)], "GCACAC", "GCACACAC");
    assert_eq!(r.cigar, vec![Match(1), Insertion(2), Match(5)]);
}

// --- build_cigar: canonical-form assembly (review M7) ----------------
//
// The right-to-left emitter can produce an insertion before a deletion and
// adjacent same-ops; `build_cigar` is what makes the output canonical, and
// it is the logic that lets identical events bucket together. Tested
// directly because the end-to-end paths rarely produce a colliding run.

#[test]
fn build_cigar_orders_deletion_before_insertion() {
    // An I-before-D run must come out D-before-I (canonical order).
    let r = build_cigar(&[Match(2), Insertion(1), Deletion(1), Match(2)], false);
    assert_eq!(r.cigar, vec![Match(2), Deletion(1), Insertion(1), Match(2)],);
}

#[test]
fn build_cigar_merges_adjacent_identical_ops_and_drops_zeros() {
    // Adjacent same-ops merge; zero-length ops are dropped.
    assert_eq!(
        build_cigar(&[Match(1), Match(0), Match(4)], false).cigar,
        vec![Match(5)],
    );
    assert_eq!(
        build_cigar(&[Deletion(1), Deletion(2)], false).cigar,
        vec![Deletion(3)],
    );
}

#[test]
fn build_cigar_strips_leading_deletion_only_when_requested() {
    // remove_deletions_at_ends gates the leading-deletion strip + its
    // reported reference length.
    let kept = build_cigar(&[Deletion(2), Match(4)], false);
    assert_eq!(kept.cigar, vec![Deletion(2), Match(4)]);
    assert_eq!(kept.leading_deletion_bases_removed, 0);

    let stripped = build_cigar(&[Deletion(2), Match(4)], true);
    assert_eq!(stripped.cigar, vec![Match(4)]);
    assert_eq!(stripped.leading_deletion_bases_removed, 2);
}

// --- left_align_indels: the in-place wrapper the reader calls ---------
//
// These drive `left_align_indels(&mut cigar, seq, ref_seq)` — the
// `(&mut Vec<CigarOp>, seq, ref_seq)` drop-in for the reader's F3 site
// (`read_start = 0`, end-deletions kept). They exercise the wrapper's
// fast path and its debug mismatch-count invariant (which `left_align_cigar`
// callers bypass), and pin the `seq`-vs-`ref_seq` argument order: a swapped
// pair would trip the invariant on these repeat fixtures.

#[test]
fn wrapper_left_aligns_homopolymer_deletion_in_place() {
    let mut cigar = vec![Match(4), Deletion(1), Match(1)];
    left_align_indels(&mut cigar, b"GAAAT", b"GAAAAT");
    assert_eq!(cigar, vec![Match(1), Deletion(1), Match(4)]);
}

#[test]
fn wrapper_left_aligns_homopolymer_insertion_in_place() {
    let mut cigar = vec![Match(4), Insertion(1), Match(1)];
    left_align_indels(&mut cigar, b"GAAAAT", b"GAAAT");
    assert_eq!(cigar, vec![Match(1), Insertion(1), Match(4)]);
}

#[test]
fn wrapper_two_reads_converge_in_place() {
    // Same biological deletion, opposite placements → identical CIGAR.
    let mut right = vec![Match(5), Deletion(1), Match(1)];
    let mut left = vec![Match(1), Deletion(1), Match(5)];
    left_align_indels(&mut right, b"GAAAAC", b"GAAAAAC");
    left_align_indels(&mut left, b"GAAAAC", b"GAAAAAC");
    assert_eq!(right, left);
    assert_eq!(right, vec![Match(1), Deletion(1), Match(5)]);
}

#[test]
fn wrapper_keeps_leading_deletion_and_is_noop_without_indel() {
    // End-deletion not stripped (alignment_start stays fixed); and a
    // pure-match CIGAR is untouched (fast path).
    let mut rolled = vec![Match(3), Deletion(1), Match(1)];
    left_align_indels(&mut rolled, b"AAAT", b"AAAAT");
    assert_eq!(rolled, vec![Deletion(1), Match(4)]);

    let mut plain = vec![Match(5)];
    left_align_indels(&mut plain, b"ACGTA", b"ACGTA");
    assert_eq!(plain, vec![Match(5)]);
}
