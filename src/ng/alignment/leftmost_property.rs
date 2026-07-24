//! **The leftmost property checker — the normalizers' oracle.**
//!
//! This module's normalizers (algorithms 1a/1b/1c, Milestones B–C) are graded against a
//! *definition*, not against each other. "Leftmost" is one: for any output, take each indel,
//! shift it one base to the left, and ask whether the result still represents the same read
//! against the same reference without introducing a mismatch. If it does, the output was **not**
//! leftmost. This is the only comparison in the whole alignment module that has ground truth
//! without a rival implementation, which is why the plan builds it **first**, before any
//! normalizer exists, and lets it **outrank agreement** between the three (spec §6, arch §Test &
//! bench shape; plan step A2).
//!
//! # Why this is not a call into `normalize_alleles`
//!
//! Production's [`normalize_alleles`](crate::norm_seqs::normalize_alleles) *is* a left-aligner,
//! and algorithm 1a will port it. If this checker were built on it, it would grade 1a against
//! itself — an oracle has to be independent of the thing it grades. So the shift test here is
//! written straight from the definition: it looks at each indel and the aligned column
//! immediately to its left, and applies the one-base-left shift condition directly.
//!
//! # The shift condition, exactly
//!
//! An indel can slide one base left only by **crossing an aligned column that is a real match**
//! — left-alignment is a *lossless* re-placement that never moves a mismatch (production's port
//! asserts the read's mismatch count is unchanged, [indel_norm.rs](crate::pileup::walker)). So an
//! indel with the aligned column `(read[q-1] ↔ reference[r-1])` immediately to its left can shift
//! one base left iff that column matches **and** the base leaving the indel on the left equals the
//! base entering it from the right:
//!
//! - **deletion** of `reference[s..e)`: `reference[s-1] == read[q-1]` (crossed column matches) and
//!   `read[q-1] == reference[e-1]` (far edge preserved) — together, `reference[s-1] ==
//!   reference[e-1]`, crossed through a matching read base;
//! - **insertion** of `read[a..b)`: `reference[r-1] == read[a-1]` and `read[a-1] == read[b-1]`.
//!
//! This is exactly the pair of predicates production's `normalize_alleles` tests at each step —
//! `next_base_on_left_is_same` and `last_base_on_right_is_same` — arrived at from the definition
//! rather than borrowed from the code, so the two agreeing is a *result*, not an assumption.
//!
//! A subtle consequence worth stating: an indel sitting in a homopolymer whose left neighbour is
//! a **mismatch** is already leftmost, because the mismatch cannot be crossed — normalization is
//! placement within a fixed match/mismatch structure, not re-alignment. The tests pin that case.
//!
//! `#[cfg(test)]`: this is test support for the normalizers. It ships no production code.

use super::Alignment;
use crate::pileup::walker::CigarOp;

/// Which kind of indel can slide left. A closed set of two, so the sibling normalizer test
/// modules that match on a witness get an exhaustive `match` and the compiler catches a
/// mistyped or unhandled kind — where a `&'static str` `"insertion"`/`"deletion"` would not.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum IndelKind {
    Insertion,
    Deletion,
}

impl IndelKind {
    /// The word used to name the offender in a test-failure message.
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::Insertion => "insertion",
            Self::Deletion => "deletion",
        }
    }
}

/// An indel that can be shifted one base to the left — i.e. a witness that an alignment is
/// **not** leftmost. Carries enough to name the offender in a test failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ShiftableIndel {
    /// Index of the offending op in the alignment's `cigar`.
    pub op_index: usize,
    /// Whether the shiftable indel is an insertion or a deletion.
    pub kind: IndelKind,
    /// The equal base whose repetition permits the leftward slide (the crossed base).
    pub crossed_base: u8,
}

/// Find an indel in `alignment` that can shift one base left over a matching column — the
/// witness that the alignment is not leftmost — or `None` if every indel is already as far left
/// as a lossless slide allows.
///
/// `read` is the full read sequence; `reference` is the whole reference stretch the alignment
/// was computed against, with [`Alignment::reference_offset`] the base the placement starts at.
/// Returns the **first** shiftable indel in read order; one witness is enough to condemn an
/// output.
///
/// **Precondition — well-formed placement.** The cigar must consume exactly `read.len()` read
/// bases, and the placement (`reference_offset` plus the reference the cigar consumes) must fit
/// inside `reference`. This holds for anything a normalizer produces from the read it was given.
/// A malformed alignment (a normalizer under test emitting a cigar inconsistent with its read)
/// is caught by the debug assertions below rather than by a bare slice index-out-of-bounds — the
/// oracle's job is good diagnostics, so a malformation should fail *named*.
pub(crate) fn find_shiftable_indel(
    alignment: &Alignment,
    read: &[u8],
    reference: &[u8],
) -> Option<ShiftableIndel> {
    debug_assert_eq!(
        consumed_read(&alignment.cigar),
        read.len(),
        "malformed alignment: cigar consumes {} read bases, read has {}",
        consumed_read(&alignment.cigar),
        read.len(),
    );
    debug_assert!(
        alignment.reference_offset as usize + consumed_reference(&alignment.cigar)
            <= reference.len(),
        "malformed alignment: placement (offset {} + {} reference bases) runs past the {}-base \
         reference stretch",
        alignment.reference_offset,
        consumed_reference(&alignment.cigar),
        reference.len(),
    );

    // `reference_position` walks the reference stretch; `read_position` walks the read. Both
    // are cursors *into the slices*, so `reference_position` starts at the placement's offset.
    let mut reference_position = alignment.reference_offset as usize;
    let mut read_position = 0usize;
    // Whether the op just before this one was an aligned column with at least one base — the
    // only kind of neighbour an indel can slide across. Reset by every non-aligned op so an
    // indel abutting a clip, a skip, or another indel is never treated as slidable.
    let mut left_neighbour_is_aligned = false;

    for (op_index, &op) in alignment.cigar.iter().enumerate() {
        match op {
            CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                let n = n as usize;
                reference_position += n;
                read_position += n;
                left_neighbour_is_aligned = n > 0;
            }
            CigarOp::Insertion(n) => {
                let n = n as usize;
                if left_neighbour_is_aligned
                    && n >= 1
                    && read_position >= 1
                    && reference_position >= 1
                {
                    let crossed_reference = reference[reference_position - 1];
                    let crossed_read = read[read_position - 1];
                    let far_edge = read[read_position + n - 1];
                    if crossed_reference == crossed_read && crossed_read == far_edge {
                        return Some(ShiftableIndel {
                            op_index,
                            kind: IndelKind::Insertion,
                            crossed_base: crossed_read,
                        });
                    }
                }
                read_position += n;
                left_neighbour_is_aligned = false;
            }
            CigarOp::Deletion(n) => {
                let n = n as usize;
                if left_neighbour_is_aligned
                    && n >= 1
                    && read_position >= 1
                    && reference_position >= 1
                {
                    let crossed_reference = reference[reference_position - 1];
                    let crossed_read = read[read_position - 1];
                    let far_edge = reference[reference_position + n - 1];
                    if crossed_reference == crossed_read && crossed_read == far_edge {
                        return Some(ShiftableIndel {
                            op_index,
                            kind: IndelKind::Deletion,
                            crossed_base: crossed_read,
                        });
                    }
                }
                reference_position += n;
                left_neighbour_is_aligned = false;
            }
            CigarOp::SoftClip(n) => {
                read_position += n as usize;
                left_neighbour_is_aligned = false;
            }
            CigarOp::Skip(n) => {
                reference_position += n as usize;
                left_neighbour_is_aligned = false;
            }
            CigarOp::HardClip(_) | CigarOp::Padding(_) => {
                left_neighbour_is_aligned = false;
            }
        }
    }
    None
}

/// Whether `alignment` is left-aligned: no indel can shift one base left over a matching column.
pub(crate) fn is_left_aligned(alignment: &Alignment, read: &[u8], reference: &[u8]) -> bool {
    find_shiftable_indel(alignment, read, reference).is_none()
}

/// How many read bases a cigar consumes (M/=/X/I/S). For the well-formedness precondition only.
fn consumed_read(cigar: &[CigarOp]) -> usize {
    cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Insertion(n)
            | CigarOp::SoftClip(n) => *n as usize,
            CigarOp::Deletion(_)
            | CigarOp::Skip(_)
            | CigarOp::HardClip(_)
            | CigarOp::Padding(_) => 0,
        })
        .sum()
}

/// How many reference bases a cigar consumes (M/=/X/D/N). For the well-formedness precondition
/// only.
fn consumed_reference(cigar: &[CigarOp]) -> usize {
    cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Deletion(n)
            | CigarOp::Skip(n) => *n as usize,
            CigarOp::Insertion(_)
            | CigarOp::SoftClip(_)
            | CigarOp::HardClip(_)
            | CigarOp::Padding(_) => 0,
        })
        .sum()
}

/// Assert `alignment` is left-aligned, panicking with the offending indel named if it is not.
/// The oracle the normalizer tests call on every output.
#[track_caller]
pub(crate) fn assert_left_aligned(alignment: &Alignment, read: &[u8], reference: &[u8]) {
    if let Some(witness) = find_shiftable_indel(alignment, read, reference) {
        panic!(
            "alignment is not leftmost: the {} at cigar op {} can shift one base left over the \
             matching base {:?}\n  alignment: {:?}\n  read: {:?}\n  reference: {:?}",
            witness.kind.label(),
            witness.op_index,
            witness.crossed_base as char,
            alignment,
            std::str::from_utf8(read).unwrap_or("<non-utf8>"),
            std::str::from_utf8(reference).unwrap_or("<non-utf8>"),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::walker::CigarOp::{Deletion, Insertion, Match, Skip, SoftClip};

    fn alignment(reference_offset: u64, cigar: Vec<CigarOp>) -> Alignment {
        Alignment {
            reference_offset,
            cigar,
        }
    }

    // ---- No indel: leftmost by having nothing to shift -------------------------------

    #[test]
    fn a_pure_match_alignment_is_trivially_leftmost() {
        let a = alignment(0, vec![Match(5)]);
        assert!(is_left_aligned(&a, b"ACGTA", b"ACGTA"));
        assert_left_aligned(&a, b"ACGTA", b"ACGTA");
    }

    #[test]
    fn an_empty_cigar_is_leftmost() {
        // The smallest input: no ops, nothing to shift. Guards against a regression that
        // assumed at least one op or indexed `cigar[0]`.
        let a = alignment(0, vec![]);
        assert!(is_left_aligned(&a, b"", b""));
        assert_left_aligned(&a, b"", b"");
    }

    // ---- Homopolymer deletion --------------------------------------------------------

    #[test]
    fn a_rightmost_homopolymer_deletion_is_detected_as_shiftable() {
        // ref GAAAAT, read GAAAT: the deletion sits at the rightmost A (4M1D1M). It can slide
        // left through the A-run, so this is NOT leftmost.
        let a = alignment(0, vec![Match(4), Deletion(1), Match(1)]);
        assert!(!is_left_aligned(&a, b"GAAAT", b"GAAAAT"));
        let witness = find_shiftable_indel(&a, b"GAAAT", b"GAAAAT").expect("shiftable");
        assert_eq!(witness.kind, IndelKind::Deletion);
        assert_eq!(witness.op_index, 1);
        assert_eq!(witness.crossed_base, b'A');
    }

    #[test]
    fn a_leftmost_homopolymer_deletion_is_accepted() {
        // Same variant placed at the leftmost A (1M1D4M): the base left of the deletion is the
        // G anchor, not another A, so it cannot slide — leftmost.
        let a = alignment(0, vec![Match(1), Deletion(1), Match(4)]);
        assert!(is_left_aligned(&a, b"GAAAT", b"GAAAAT"));
        assert_left_aligned(&a, b"GAAAT", b"GAAAAT");
    }

    // ---- Homopolymer insertion -------------------------------------------------------

    #[test]
    fn a_rightmost_homopolymer_insertion_is_detected_as_shiftable() {
        // ref GAAAT, read GAAAAT: the inserted A sits at the right (4M1I1M) and can slide left.
        let a = alignment(0, vec![Match(4), Insertion(1), Match(1)]);
        assert!(!is_left_aligned(&a, b"GAAAAT", b"GAAAT"));
        let witness = find_shiftable_indel(&a, b"GAAAAT", b"GAAAT").expect("shiftable");
        assert_eq!(witness.kind, IndelKind::Insertion);
        assert_eq!(witness.op_index, 1);
        assert_eq!(witness.crossed_base, b'A');
    }

    #[test]
    fn a_leftmost_homopolymer_insertion_is_accepted() {
        // Same insertion at the leftmost A (1M1I4M): the base to its left is the G anchor.
        let a = alignment(0, vec![Match(1), Insertion(1), Match(4)]);
        assert!(is_left_aligned(&a, b"GAAAAT", b"GAAAT"));
    }

    #[test]
    fn a_trailing_insertion_at_the_read_end_is_evaluated_in_bounds() {
        // ref GAAA, read GAAAA (4M1I): the inserted A is the last read base, so the far-edge
        // index `read[read_position + n - 1]` lands on the final read byte. An off-by-one-high
        // far edge would read past the read and panic; this pins the boundary and that a
        // trailing shiftable insertion is still found.
        let a = alignment(0, vec![Match(4), Insertion(1)]);
        assert!(!is_left_aligned(&a, b"GAAAA", b"GAAA"));
    }

    // ---- Period-2 repeat -------------------------------------------------------------

    #[test]
    fn a_rightmost_period_two_deletion_is_detected_as_shiftable() {
        // ref CACACAG, read CACAG: one CA unit deleted, placed at the right (4M2D1M). The 2-base
        // deletion can slide one base left through the CA repeat.
        let a = alignment(0, vec![Match(4), Deletion(2), Match(1)]);
        assert!(!is_left_aligned(&a, b"CACAG", b"CACACAG"));
        let witness = find_shiftable_indel(&a, b"CACAG", b"CACACAG").expect("shiftable");
        assert_eq!(witness.kind, IndelKind::Deletion);
    }

    #[test]
    fn a_rightmost_period_two_insertion_is_detected_as_shiftable() {
        // ref CACAG, read CACACAG: one CA unit inserted at the right (4M2I1M). The far-edge base
        // is read[5]='A'; a mutation that read read[read_position]='C' (dropping the `+ n - 1`)
        // would wrongly report leftmost — only a multi-base insertion catches it, which the
        // single-base insertion tests cannot.
        let a = alignment(0, vec![Match(4), Insertion(2), Match(1)]);
        assert!(!is_left_aligned(&a, b"CACACAG", b"CACAG"));
        let witness = find_shiftable_indel(&a, b"CACACAG", b"CACAG").expect("shiftable");
        assert_eq!(witness.kind, IndelKind::Insertion);
        assert_eq!(witness.crossed_base, b'A');
    }

    #[test]
    fn a_leftmost_period_two_deletion_is_accepted() {
        // The same unit deleted at the very start (2D5M): nothing precedes it, so it cannot
        // slide — leftmost.
        let a = alignment(0, vec![Deletion(2), Match(5)]);
        assert!(is_left_aligned(&a, b"CACAG", b"CACACAG"));
        assert_left_aligned(&a, b"CACAG", b"CACACAG");
    }

    #[test]
    fn a_one_short_period_two_deletion_still_shifts() {
        // ref CACACAG, read CACAG, placed as 1M2D4M (delete AC at ref[1..3]). It is one base off
        // leftmost and must still be reported shiftable — the single-step test catches it.
        let a = alignment(0, vec![Match(1), Deletion(2), Match(4)]);
        assert!(!is_left_aligned(&a, b"CACAG", b"CACACAG"));
    }

    // ---- The match-required clause (the subtle one) ----------------------------------

    #[test]
    fn a_deletion_in_a_run_whose_left_neighbour_is_a_mismatch_is_leftmost() {
        // ref CAA, read CT, aligned 2M1D: read C↔ref C (match), read T↔ref A (MISMATCH), then
        // delete the last A. The A-run would *look* slidable, but the column left of the
        // deletion is a mismatch, and left-alignment never crosses a mismatch — so this is
        // already leftmost. This is the clause that separates "homopolymer" from "shiftable".
        let a = alignment(0, vec![Match(2), Deletion(1)]);
        assert!(is_left_aligned(&a, b"CT", b"CAA"));
        assert_left_aligned(&a, b"CT", b"CAA");
    }

    #[test]
    fn a_deletion_whose_far_edge_matches_but_left_column_mismatches_is_leftmost() {
        // ref TCA, read TA, aligned 2M1D: read[1]='A' ↔ ref[1]='C' is a MISMATCH, then delete
        // ref[2]='A'. Here the far edge (ref[2]='A') *equals* the crossed read base (read[1]=
        // 'A'), so the second clause of the shift condition passes — only the crossed-column-
        // match clause keeps this leftmost. Deleting that clause (crossed_reference ==
        // crossed_read) would report it shiftable. The mismatch-run test above cannot catch that,
        // because there both clauses are false; this one pins the match clause alone.
        let a = alignment(0, vec![Match(2), Deletion(1)]);
        assert!(is_left_aligned(&a, b"TA", b"TCA"));
        assert_left_aligned(&a, b"TA", b"TCA");
    }

    // ---- Boundary and cursor handling ------------------------------------------------

    #[test]
    fn an_indel_at_the_very_start_cannot_shift() {
        // A leading insertion has no aligned column to its left.
        let a = alignment(0, vec![Insertion(2), Match(3)]);
        assert!(is_left_aligned(&a, b"AACGT", b"CGT"));
    }

    #[test]
    fn an_indel_abutting_a_preceding_indel_is_not_treated_as_slidable() {
        // The `left_neighbour_is_aligned = false` reset in the indel arms is load-bearing: a
        // second indel whose left neighbour is another indel (not an aligned column) must never
        // be judged slidable. Here the first M is a MISMATCH, so the first indel is not
        // shiftable and read order does not stop there — leaving the second indel's gate as the
        // only thing keeping the alignment leftmost.
        //
        // ref ATTG, read TG, 1M2D1M spelled as two Deletion(1): M(read T ↔ ref A, mismatch),
        // delete ref[1]=T, delete ref[2]=T, M(read G ↔ ref G). If the deletion arm forgot the
        // reset, the second deletion would false-positive over the T-run and fail this
        // correct-and-leftmost alignment.
        let two_deletions = alignment(0, vec![Match(1), Deletion(1), Deletion(1), Match(1)]);
        assert!(is_left_aligned(&two_deletions, b"TG", b"ATTG"));

        // Symmetric guard for the insertion arm: M(mismatch), I, D, M. If the insertion arm
        // forgot the reset, the following deletion would false-positive.
        let insertion_then_deletion =
            alignment(0, vec![Match(1), Insertion(1), Deletion(1), Match(1)]);
        assert!(is_left_aligned(&insertion_then_deletion, b"ACG", b"CCG"));
    }

    #[test]
    fn a_shiftable_deletion_is_found_using_the_reference_offset_not_the_stretch_start() {
        // The placement starts at index 4 (the G) of CCCCGAAAAT. The rightmost-A deletion
        // (4M1D1M) is shiftable ONLY when read from reference[4..]; reading from reference[0..]
        // lands the crossed column on 'C' (a mismatch against read[3]='A') and wrongly reports
        // leftmost. This fixture's prefix base differs from the run, so dropping `reference_offset`
        // flips the verdict — unlike a homopolymer prefix, which would hide the bug.
        let a = alignment(4, vec![Match(4), Deletion(1), Match(1)]);
        assert!(!is_left_aligned(&a, b"GAAAT", b"CCCCGAAAAT"));
        let witness = find_shiftable_indel(&a, b"GAAAT", b"CCCCGAAAAT").expect("shiftable");
        assert_eq!(witness.crossed_base, b'A');
    }

    #[test]
    fn a_second_deletion_uses_the_first_deletions_reference_advance() {
        // ref GCAACT, read GAAT, 1M1D2M1D1M. The FIRST deletion (of ref[1]='C') is not shiftable
        // (G ≠ C). The SECOND deletion (of ref[4]='C') sits at reference_position 4 only because
        // the first deletion advanced the cursor; there its far edge is ref[4]='C' ≠ 'A', so it
        // is leftmost. If the deletion arm forgot `reference_position += n`, the second deletion
        // would be evaluated at reference_position 3, where its far edge is ref[3]='A' — a match
        // — and it would false-positive. So this pins the deletion reference-advance.
        let a = alignment(
            0,
            vec![Match(1), Deletion(1), Match(2), Deletion(1), Match(1)],
        );
        assert!(is_left_aligned(&a, b"GAAT", b"GCAACT"));
    }

    #[test]
    fn a_soft_clip_advances_the_read_cursor_so_the_crossed_read_base_is_correct() {
        // 2S4M1D1M; read = XX GCAA T, ref = GCAAAT. The crossed read base is read[5]='A' only
        // when the clip advanced the read cursor; without the advance it is read[3]='C' (a
        // mismatch against ref[3]='A') and the shiftable deletion is missed. A homopolymer
        // fixture cannot catch that, because the clipped read bytes would still be 'A'.
        let a = alignment(0, vec![SoftClip(2), Match(4), Deletion(1), Match(1)]);
        assert!(!is_left_aligned(&a, b"XXGCAAT", b"GCAAAT"));

        // The clip is not an aligned column, so an indel immediately after it cannot shift.
        let right_after_clip = alignment(0, vec![SoftClip(2), Deletion(1), Match(5)]);
        assert!(is_left_aligned(&right_after_clip, b"XXAAAAT", b"AAAAAT"));
    }

    #[test]
    fn a_skip_advances_the_reference_cursor() {
        // 2M3N4M1D1M with an intron (Skip) between the flanks. The post-skip deletion is a
        // shiftable homopolymer deletion, but only when the Skip advanced reference_position by
        // 3 — without it, the crossed base lands on 'G' (a mismatch) and the deletion is missed.
        // ref = TT CCC GAAA A T ; read = TT GAAA T.
        let a = alignment(0, vec![Match(2), Skip(3), Match(4), Deletion(1), Match(1)]);
        assert!(!is_left_aligned(&a, b"TTGAAAT", b"TTCCCGAAAAT"));
    }

    #[test]
    fn assert_left_aligned_panics_on_a_shiftable_alignment() {
        let a = alignment(0, vec![Match(4), Deletion(1), Match(1)]);
        let result = std::panic::catch_unwind(|| assert_left_aligned(&a, b"GAAAT", b"GAAAAT"));
        assert!(
            result.is_err(),
            "expected a panic on a non-leftmost alignment"
        );
    }
}
