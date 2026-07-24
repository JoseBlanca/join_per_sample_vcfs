//! Algorithm 1b — **repeated simple passes**, freebayes' shape, an [`AlignmentNormalizer`].
//!
//! Where algorithm 1a ([`StructuredLeftAligner`](super::left_align_structured)) does the whole job
//! in one structured right-to-left walk, 1b takes freebayes' approach: run a **simple** shift pass,
//! and re-run it until nothing moves, bounded by a cap (freebayes uses 20,
//! `LeftAlign.h:118`). Implemented **from the description, not by transliterating** freebayes'
//! `LeftAlign.cpp` (spec §8), and deliberately kept independent of 1a and of production's
//! `left_align_indels` — it is a rival, so it may not lean on what it competes with.
//!
//! # What a "simple pass" is here
//!
//! One pass slides **every shiftable indel one base to the left** and stops. An indel can slide one
//! base left only across an aligned column that is a real **match**, with the base leaving the indel
//! on the left equal to the base entering from the right — the same definition the property oracle
//! grades against (spec §6), implemented here independently because a rival must not borrow the
//! oracle. Sliding a gap through a run therefore takes one pass per base, which is why the passes
//! are *repeated*; a pass that moves nothing is the fixed point.
//!
//! **A recorded consequence of "one base per pass".** freebayes slides each indel *maximally* per
//! pass, so its cap bounds the number of cross-indel re-slides. This one-base pass instead makes the
//! cap bound the total shift **distance**: an indel more than `max_passes` bases from its leftmost
//! home reports exhaustion where 1a would place it in a single pass. That is a real, deliberate
//! difference — it is exactly the kind of thing the differ-at-all screen (Milestone D) exists to
//! surface — and it makes the exhaustion path reachable and testable rather than only theoretical.
//!
//! **A second recorded difference from 1a: 1b does not trim a complex indel.** 1b canonicalizes its
//! output the way 1a does — merging same-kind gaps, ordering a deletion before an abutting insertion,
//! dropping zero-length ops — so on a *pure* indel (the common case, and essentially all real reads)
//! it produces exactly 1a's spelling. It does **not**, however, trim a deletion/insertion overlap
//! whose bases cancel (`D2 I1` → `D1`), which 1a inherits from `normalize_alleles`. So for a complex
//! indel the two can differ in more than placement. Both remain leftmost (the property oracle grades
//! placement, not parsimony); the difference is left un-eliminated deliberately — trimming is a
//! separate normalization axis this "simple repeated shift" step does not own, and the screen
//! surfaces the difference rather than hiding it. **Recorded for the owner at Checkpoint C**: whether
//! 1b should also trim, for full parity with 1a, is a scope call, not a silent omission.
//!
//! # Exhaustion is reported, not swallowed
//!
//! freebayes computes a convergence flag and its own caller ignores it, so non-convergent results
//! ship silently (spec §6). Here the flag is the **return value** of [`RepeatedLeftAligner::left_align`]
//! ([`ConvergenceReport`]): a caller cannot drop it by accident. The [`AlignmentNormalizer`] trait's
//! `normalize` is a uniform interface that cannot carry a report, so it discards it *visibly* — but
//! the reporting method is the real API, and it is what the tests and the screen use. (Algorithm 1c
//! makes the opposite choice on the cap: it **panics** rather than returning a half-normalised
//! alignment.)

use super::{Alignment, AlignmentNormalizer};
use crate::pileup::walker::CigarOp;

/// The iteration cap: the most passes [`RepeatedLeftAligner`] runs before reporting exhaustion.
/// freebayes' `stablyLeftAlign` uses 20 (`LeftAlign.h:118`); adopted as the default.
pub const MAX_PASSES: u32 = 20;

/// Whether a repeated-pass left-alignment reached a fixed point within its cap — the report 1b
/// returns instead of swallowing (contrast freebayes, whose caller ignores the same fact).
///
/// `#[must_use]` is the compile-time half of "a caller cannot drop it by accident": a bare
/// `aligner.left_align(…);` statement warns unless the report is handled or discarded explicitly.
/// The trait's `normalize` is the one sanctioned discard (an explicit `let _ =`).
#[must_use = "an ExhaustedCap result is not known to be leftmost; handle it or discard it explicitly"]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvergenceReport {
    /// A pass moved nothing: the alignment is left-aligned. `passes` counts every pass run,
    /// including the final confirming no-move pass (so an already-leftmost input reports 1).
    Converged { passes: u32 },
    /// The cap was reached with a pass still moving the alignment. The result is **not** known to
    /// be leftmost; a caller must not treat it as canonical.
    ExhaustedCap,
}

impl ConvergenceReport {
    /// Whether a fixed point was reached (as opposed to hitting the cap).
    #[must_use]
    pub fn converged(self) -> bool {
        matches!(self, Self::Converged { .. })
    }
}

/// Algorithm 1b: left-align by re-running a simple one-base shift pass until nothing moves, capped.
///
/// Stateless apart from its cap. `max_passes` defaults to [`MAX_PASSES`]; a smaller value is used
/// in tests to drive the exhaustion path deliberately.
#[derive(Debug, Clone, Copy)]
pub struct RepeatedLeftAligner {
    max_passes: u32,
}

impl Default for RepeatedLeftAligner {
    fn default() -> Self {
        Self {
            max_passes: MAX_PASSES,
        }
    }
}

impl RepeatedLeftAligner {
    /// The default aligner, capped at [`MAX_PASSES`].
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// An aligner with a custom cap — the constructor tests use to force exhaustion.
    ///
    /// A `max_passes` of 0 is legal but degenerate: it runs no passes and reports
    /// [`ConvergenceReport::ExhaustedCap`] without touching the alignment (the honest "not known to
    /// be leftmost" answer for "no work permitted"). Real use wants at least 1.
    #[must_use]
    pub fn with_max_passes(max_passes: u32) -> Self {
        Self { max_passes }
    }

    /// The cap this aligner runs to.
    #[must_use]
    pub fn max_passes(self) -> u32 {
        self.max_passes
    }

    /// Left-align `alignment.cigar` in place with repeated simple passes, **returning** whether it
    /// converged within the cap.
    ///
    /// This is the reporting entry point — the fact freebayes drops is 1b's return value. `read`
    /// and `reference` are the full read and the whole reference stretch; `reference_offset` is
    /// where the placement starts (the offset is never moved — a deletion rolling to the read start
    /// stays a first-op `Deletion`, which is already leftmost).
    ///
    /// # Panics
    ///
    /// If `reference_offset > reference.len()` — the placement runs past the reference stretch. That
    /// is a fail-fast: a silently mis-sliced reference would be a wrong variant. A debug assertion
    /// names the precondition; the slice's own bounds check is what still panics in release.
    pub fn left_align(
        &self,
        alignment: &mut Alignment,
        read: &[u8],
        reference: &[u8],
    ) -> ConvergenceReport {
        let offset = alignment.reference_offset as usize;
        debug_assert!(
            offset <= reference.len(),
            "reference_offset {offset} exceeds reference length {}",
            reference.len(),
        );
        let reference_from_placement = &reference[offset..];
        for pass in 1..=self.max_passes {
            if !shift_pass(&mut alignment.cigar, read, reference_from_placement) {
                return ConvergenceReport::Converged { passes: pass };
            }
        }
        ConvergenceReport::ExhaustedCap
    }
}

impl AlignmentNormalizer for RepeatedLeftAligner {
    /// The uniform-interface entry point. It runs [`Self::left_align`] and **discards the
    /// convergence report** — the interface cannot carry it. The discard is deliberate and visible;
    /// a caller that needs to know whether the result is leftmost calls `left_align` directly. 1b
    /// never panics on exhaustion — that is 1c's contract, and the two must differ.
    fn normalize(&self, alignment: &mut Alignment, read: &[u8], reference: &[u8]) {
        let _ = self.left_align(alignment, read, reference);
    }
}

/// Run one simple pass: slide every shiftable indel one base to the left, then put the result in
/// canonical form. Returns whether the cigar changed (an unchanged pass is the fixed point).
///
/// Canonicalization runs on **every** pass, including the one where nothing is shiftable — an
/// already-leftmost input may still be spelled non-canonically (split same-kind gaps, an insertion
/// before an abutting deletion, a zero-length op), and 1b must consolidate it exactly as 1a does, or
/// two reads carrying the same event scatter across cohort-merge buckets. Movement is therefore
/// measured against the canonical result: a pass that only merges a non-canonical spelling counts as
/// motion (one canonicalising pass, then a confirming no-op).
fn shift_pass(cigar: &mut Vec<CigarOp>, read: &[u8], reference: &[u8]) -> bool {
    let shiftable = shiftable_indices(cigar, read, reference);
    // Rebuild: an alignment op immediately before a shiftable indel gives up its last base (that
    // base becomes an aligned column on the indel's right after the slide), spelled as a `Match(1)`
    // emitted right after the indel. With nothing shiftable this is a straight copy that
    // `canonicalize` still tidies.
    let mut shifted: Vec<CigarOp> = Vec::with_capacity(cigar.len() + shiftable.len() + 2);
    for (index, &op) in cigar.iter().enumerate() {
        if shiftable.binary_search(&(index + 1)).is_ok() {
            // The next op is a shiftable indel, so this (necessarily alignment) op loses one base.
            shifted.push(with_len(op, op_len(op) - 1));
        } else {
            shifted.push(op);
        }
        if shiftable.binary_search(&index).is_ok() {
            shifted.push(CigarOp::Match(1));
        }
    }
    let canonical = canonicalize(&shifted);
    let moved = canonical != *cigar;
    *cigar = canonical;
    moved
}

/// The cigar indices of every indel that can slide one base to the left: an insertion or deletion
/// whose immediately preceding op is an aligned column that **matches**, with the base leaving the
/// indel on the left equal to the base entering from the right. This is the same definition the
/// property oracle uses (spec §6), implemented independently because 1b is a rival to what the
/// oracle grades. Indices are returned in ascending order.
fn shiftable_indices(cigar: &[CigarOp], read: &[u8], reference: &[u8]) -> Vec<usize> {
    let mut shiftable = Vec::new();
    let mut reference_position = 0usize;
    let mut read_position = 0usize;
    let mut left_neighbour_is_aligned = false;

    for (index, &op) in cigar.iter().enumerate() {
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
                        shiftable.push(index);
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
                        shiftable.push(index);
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
    shiftable
}

/// Put a cigar into canonical form: drop zero-length ops, order a deletion before an abutting
/// insertion, and merge adjacent same-kind ops. Mirrors production's `build_cigar` (`indel_norm.rs`)
/// minus the end-deletion stripping 1b never does — duplicated because production's is private and
/// frozen, and a rival must not lean on it.
///
/// **It does *not* trim a deletion/insertion overlap.** Production's 1a reduces an abutting `D`/`I`
/// whose bases cancel (via `normalize_alleles`' parsimony trim), turning e.g. `D2 I1` into `D1`; 1b
/// does not, so for a *complex* indel the two can differ in more than placement (see the module doc's
/// note). Both spellings are still leftmost — the property oracle grades placement, not parsimony —
/// and the difference is itself something the differ-at-all screen legitimately surfaces.
fn canonicalize(ops: &[CigarOp]) -> Vec<CigarOp> {
    let mut nonzero: Vec<CigarOp> = ops.iter().copied().filter(|op| op_len(*op) > 0).collect();

    // Within each maximal run of consecutive indels, put the (single, merged) deletion before the
    // (single, merged) insertion.
    let mut i = 0;
    while i < nonzero.len() {
        if is_indel(nonzero[i]) {
            let mut j = i;
            while j < nonzero.len() && is_indel(nonzero[j]) {
                j += 1;
            }
            let mut deletion_len = 0u32;
            let mut insertion_len = 0u32;
            for op in &nonzero[i..j] {
                match op {
                    CigarOp::Deletion(n) => deletion_len += n,
                    CigarOp::Insertion(n) => insertion_len += n,
                    // UNREACHABLE: `j` advanced only across `is_indel` ops.
                    _ => unreachable!("indel run holds only insertions and deletions"),
                }
            }
            let mut replacement = Vec::with_capacity(2);
            if deletion_len > 0 {
                replacement.push(CigarOp::Deletion(deletion_len));
            }
            if insertion_len > 0 {
                replacement.push(CigarOp::Insertion(insertion_len));
            }
            let replacement_len = replacement.len();
            nonzero.splice(i..j, replacement);
            i += replacement_len;
        } else {
            i += 1;
        }
    }

    // Merge adjacent same-kind ops.
    let mut merged: Vec<CigarOp> = Vec::with_capacity(nonzero.len());
    for op in nonzero {
        if let Some(last) = merged.last_mut()
            && std::mem::discriminant(last) == std::mem::discriminant(&op)
        {
            *last = with_len(*last, op_len(*last) + op_len(op));
        } else {
            merged.push(op);
        }
    }
    merged
}

// The three helpers below (and `canonicalize` above) duplicate private twins in production's
// `pileup::walker::indel_norm` — `op_len`, `with_len`, `is_indel`, `build_cigar`. They are copied,
// not imported, because those are `fn`-private and production is frozen (arch §5): a rival must not
// lean on what it competes with. Debt: the four move together, and the port-back of this module is
// where they unify.

/// The length carried by any cigar op.
fn op_len(op: CigarOp) -> u32 {
    match op {
        CigarOp::Match(n)
        | CigarOp::Insertion(n)
        | CigarOp::Deletion(n)
        | CigarOp::Skip(n)
        | CigarOp::SoftClip(n)
        | CigarOp::HardClip(n)
        | CigarOp::Padding(n)
        | CigarOp::SeqMatch(n)
        | CigarOp::SeqMismatch(n) => n,
    }
}

/// The same op kind with a new length.
fn with_len(op: CigarOp, n: u32) -> CigarOp {
    match op {
        CigarOp::Match(_) => CigarOp::Match(n),
        CigarOp::Insertion(_) => CigarOp::Insertion(n),
        CigarOp::Deletion(_) => CigarOp::Deletion(n),
        CigarOp::Skip(_) => CigarOp::Skip(n),
        CigarOp::SoftClip(_) => CigarOp::SoftClip(n),
        CigarOp::HardClip(_) => CigarOp::HardClip(n),
        CigarOp::Padding(_) => CigarOp::Padding(n),
        CigarOp::SeqMatch(_) => CigarOp::SeqMatch(n),
        CigarOp::SeqMismatch(_) => CigarOp::SeqMismatch(n),
    }
}

/// Whether an op is an insertion or deletion.
///
/// Written as an exhaustive `match`, not `matches!`, deliberately: `canonicalize`'s `unreachable!`
/// relies on this covering exactly the indel ops, so a future `CigarOp` variant must force a
/// decision here at compile time rather than be silently classified as "not an indel".
fn is_indel(op: CigarOp) -> bool {
    match op {
        CigarOp::Insertion(_) | CigarOp::Deletion(_) => true,
        CigarOp::Match(_)
        | CigarOp::Skip(_)
        | CigarOp::SoftClip(_)
        | CigarOp::HardClip(_)
        | CigarOp::Padding(_)
        | CigarOp::SeqMatch(_)
        | CigarOp::SeqMismatch(_) => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::alignment::left_align_structured::StructuredLeftAligner;
    use crate::ng::alignment::leftmost_property::{assert_left_aligned, is_left_aligned};
    use crate::pileup::walker::CigarOp::{Deletion, Insertion, Match};

    fn align(reference_offset: u64, cigar: Vec<CigarOp>) -> Alignment {
        Alignment {
            reference_offset,
            cigar,
        }
    }

    /// Normalize (default cap) and return the resulting cigar.
    fn normalized_cigar(
        reference_offset: u64,
        cigar: Vec<CigarOp>,
        read: &[u8],
        reference: &[u8],
    ) -> Vec<CigarOp> {
        let mut alignment = align(reference_offset, cigar);
        RepeatedLeftAligner::new().normalize(&mut alignment, read, reference);
        alignment.cigar
    }

    #[test]
    fn a_pure_match_alignment_is_returned_unchanged_and_converges_in_one_pass() {
        let mut alignment = align(0, vec![Match(5)]);
        let report = RepeatedLeftAligner::new().left_align(&mut alignment, b"ACGTA", b"ACGTA");
        assert_eq!(alignment.cigar, vec![Match(5)]);
        assert_eq!(report, ConvergenceReport::Converged { passes: 1 });
    }

    #[test]
    fn a_homopolymer_deletion_is_left_aligned() {
        let cigar = normalized_cigar(
            0,
            vec![Match(4), Deletion(1), Match(1)],
            b"GAAAT",
            b"GAAAAT",
        );
        assert_eq!(cigar, vec![Match(1), Deletion(1), Match(4)]);
    }

    #[test]
    fn a_homopolymer_insertion_is_left_aligned() {
        let cigar = normalized_cigar(
            0,
            vec![Match(4), Insertion(1), Match(1)],
            b"GAAAAT",
            b"GAAAT",
        );
        assert_eq!(cigar, vec![Match(1), Insertion(1), Match(4)]);
    }

    #[test]
    fn the_pass_count_grows_with_the_shift_distance() {
        // The deletion sits three bases right of its leftmost home, so three one-base passes move
        // it and a fourth confirms — the defining property of "one base per pass".
        let mut alignment = align(0, vec![Match(4), Deletion(1), Match(1)]);
        let report = RepeatedLeftAligner::new().left_align(&mut alignment, b"GAAAT", b"GAAAAT");
        assert_eq!(alignment.cigar, vec![Match(1), Deletion(1), Match(4)]);
        assert_eq!(report, ConvergenceReport::Converged { passes: 4 });
    }

    #[test]
    fn a_deletion_rolling_to_the_read_start_keeps_a_nonzero_reference_offset() {
        // reference[3..] = "AAAAT", read "AAAT": the deletion rolls to a leading first-op Deletion;
        // 1b, like 1a, never moves the offset.
        let mut alignment = align(3, vec![Match(3), Deletion(1), Match(1)]);
        RepeatedLeftAligner::new().normalize(&mut alignment, b"AAAT", b"TTTAAAAT");
        assert_eq!(alignment.reference_offset, 3);
        assert_eq!(alignment.cigar, vec![Deletion(1), Match(4)]);
    }

    // ---- Exhaustion is reported, not swallowed (plan step C1) ------------------------

    #[test]
    fn a_cap_too_small_to_finish_reports_exhaustion() {
        // The deletion needs three one-base passes, but the cap is two: it must report
        // ExhaustedCap, not silently ship a half-shifted (non-leftmost) result.
        let mut alignment = align(0, vec![Match(4), Deletion(1), Match(1)]);
        let report =
            RepeatedLeftAligner::with_max_passes(2).left_align(&mut alignment, b"GAAAT", b"GAAAAT");
        assert_eq!(report, ConvergenceReport::ExhaustedCap);
        assert!(!report.converged());
        // The alignment did move (it is not left at its input), but it is NOT yet leftmost — the
        // whole point of reporting exhaustion.
        assert_ne!(alignment.cigar, vec![Match(4), Deletion(1), Match(1)]);
        assert!(!is_left_aligned(&alignment, b"GAAAT", b"GAAAAT"));
    }

    #[test]
    fn a_cap_exactly_large_enough_converges() {
        // Three passes move the deletion, a fourth confirms: a cap of four is exactly enough.
        let mut alignment = align(0, vec![Match(4), Deletion(1), Match(1)]);
        let report =
            RepeatedLeftAligner::with_max_passes(4).left_align(&mut alignment, b"GAAAT", b"GAAAAT");
        assert_eq!(report, ConvergenceReport::Converged { passes: 4 });
        assert!(is_left_aligned(&alignment, b"GAAAT", b"GAAAAT"));
    }

    #[test]
    fn a_cap_equal_to_the_shift_distance_reports_exhaustion_even_though_leftmost() {
        // The deliberate, surprising boundary: a deletion three bases from home needs three passes
        // that each move, but the cap of three leaves no room for the confirming no-move pass — so
        // the result IS leftmost yet the report is ExhaustedCap. This is why MAX_PASSES = 20 confirms
        // only indels needing <= 19 shifts; a refactor of the `1..=max_passes` bound would silently
        // move this boundary.
        let mut alignment = align(0, vec![Match(4), Deletion(1), Match(1)]);
        let report =
            RepeatedLeftAligner::with_max_passes(3).left_align(&mut alignment, b"GAAAT", b"GAAAAT");
        assert_eq!(report, ConvergenceReport::ExhaustedCap);
        assert!(is_left_aligned(&alignment, b"GAAAT", b"GAAAAT"));
    }

    #[test]
    fn a_zero_cap_does_nothing_and_reports_exhaustion() {
        // The degenerate cap: no passes run, so an already-leftmost input is left untouched and the
        // honest report is ExhaustedCap ("not known to be leftmost — no work permitted").
        let mut alignment = align(0, vec![Match(4), Deletion(1), Match(1)]);
        let report =
            RepeatedLeftAligner::with_max_passes(0).left_align(&mut alignment, b"GAAAT", b"GAAAAT");
        assert_eq!(report, ConvergenceReport::ExhaustedCap);
        assert_eq!(alignment.cigar, vec![Match(4), Deletion(1), Match(1)]);
    }

    // ---- The property oracle grades this rival independently (spec §6) ----------------

    #[test]
    fn every_converged_output_passes_the_leftmost_property() {
        struct Case {
            offset: u64,
            cigar: Vec<CigarOp>,
            read: &'static [u8],
            reference: &'static [u8],
        }
        let cases = [
            Case {
                offset: 0,
                cigar: vec![Match(4), Deletion(1), Match(1)],
                read: b"GAAAT",
                reference: b"GAAAAT",
            },
            Case {
                offset: 0,
                cigar: vec![Match(4), Insertion(1), Match(1)],
                read: b"GAAAAT",
                reference: b"GAAAT",
            },
            Case {
                offset: 0,
                cigar: vec![Match(5), Deletion(2)],
                read: b"CATAT",
                reference: b"CATATAT",
            },
            Case {
                offset: 0,
                cigar: vec![Match(6), Insertion(2), Match(1)],
                read: b"GCACACACT",
                reference: b"GCACACT",
            },
            Case {
                offset: 0,
                cigar: vec![Match(2), Deletion(1), Match(2)],
                read: b"GATT",
                reference: b"GACTT",
            },
            Case {
                offset: 2,
                cigar: vec![Match(4), Deletion(1), Match(1)],
                read: b"GAAAT",
                reference: b"TTGAAAAT",
            },
            Case {
                offset: 0,
                cigar: vec![Match(1), Deletion(1), Match(1), Deletion(1), Match(2)],
                read: b"AAAT",
                reference: b"AAAAAT",
            },
        ];
        for (index, case) in cases.iter().enumerate() {
            let mut alignment = align(case.offset, case.cigar.clone());
            let report =
                RepeatedLeftAligner::new().left_align(&mut alignment, case.read, case.reference);
            assert!(report.converged(), "case {index}: unexpectedly hit the cap");
            assert!(
                is_left_aligned(&alignment, case.read, case.reference),
                "case {index}: converged output is not leftmost: {:?}",
                alignment.cigar,
            );
            assert_left_aligned(&alignment, case.read, case.reference);
        }
    }

    // ---- Agreement with algorithm 1a, the other converging normalizer -----------------

    #[test]
    fn converged_output_agrees_with_the_structured_pass() {
        // When both converge, the two rivals must place the indel identically — this is what the
        // differ-at-all screen (Milestone D) measures at scale. Any disagreement here would be a
        // finding; there is none on these fixtures.
        struct Case {
            offset: u64,
            cigar: Vec<CigarOp>,
            read: &'static [u8],
            reference: &'static [u8],
        }
        let cases = [
            Case {
                offset: 0,
                cigar: vec![Match(4), Deletion(1), Match(1)],
                read: b"GAAAT",
                reference: b"GAAAAT",
            },
            Case {
                offset: 0,
                cigar: vec![Match(4), Insertion(1), Match(1)],
                read: b"GAAAAT",
                reference: b"GAAAT",
            },
            Case {
                offset: 0,
                cigar: vec![Match(5), Deletion(2)],
                read: b"CATAT",
                reference: b"CATATAT",
            },
            Case {
                offset: 0,
                cigar: vec![Match(1), Deletion(1), Match(1), Deletion(1), Match(2)],
                read: b"AAAT",
                reference: b"AAAAAT",
            },
            Case {
                offset: 3,
                cigar: vec![Match(3), Deletion(1), Match(1)],
                read: b"AAAT",
                reference: b"TTTAAAAT",
            },
            // Non-canonical but already-leftmost input: split same-kind insertions. Both must
            // consolidate to `[I3, M2]` — the case the pre-fix skip-canonicalize path got wrong.
            Case {
                offset: 0,
                cigar: vec![Insertion(1), Insertion(1), Insertion(1), Match(1), Match(1)],
                read: b"AAAAC",
                reference: b"AC",
            },
        ];
        for (index, case) in cases.iter().enumerate() {
            let repeated =
                normalized_cigar(case.offset, case.cigar.clone(), case.read, case.reference);

            let mut structured = align(case.offset, case.cigar.clone());
            StructuredLeftAligner.normalize(&mut structured, case.read, case.reference);

            assert_eq!(
                repeated, structured.cigar,
                "case {index}: 1b and 1a disagree",
            );
        }
    }

    #[test]
    fn an_already_leftmost_but_non_canonical_input_is_consolidated() {
        // Two-plus leading insertions, already leftmost (nothing to shift) but spelled split. The
        // fixed `shift_pass` canonicalizes on the no-shift path too, so 1b merges them — a
        // regression to the early-return would leave the input verbatim and scatter cohort buckets.
        let cigar = normalized_cigar(
            0,
            vec![Insertion(1), Insertion(1), Insertion(1), Match(1), Match(1)],
            b"AAAAC",
            b"AC",
        );
        assert_eq!(cigar, vec![Insertion(3), Match(2)]);
    }

    #[test]
    fn a_complex_indel_overlap_is_not_trimmed_and_may_differ_from_1a_in_spelling() {
        // The recorded difference from 1a: 1b canonicalizes but does NOT trim an overlapping
        // deletion/insertion. Here 1b converges to a leftmost `[D2, I1, M3]` while 1a trims it to
        // `[D1, M4]`. Both are leftmost (the oracle grades placement, not parsimony); the spelling
        // differs. This test pins the difference so it stays deliberate, not accidental.
        let repeated = normalized_cigar(
            0,
            vec![
                Match(1),
                Deletion(1),
                Deletion(1),
                Insertion(1),
                Match(1),
                Match(1),
            ],
            b"AAAC",
            b"AAAAC",
        );
        assert_eq!(repeated, vec![Deletion(2), Insertion(1), Match(3)]);

        let mut structured = align(
            0,
            vec![
                Match(1),
                Deletion(1),
                Deletion(1),
                Insertion(1),
                Match(1),
                Match(1),
            ],
        );
        StructuredLeftAligner.normalize(&mut structured, b"AAAC", b"AAAAC");
        assert_eq!(structured.cigar, vec![Deletion(1), Match(4)]);

        // Both are leftmost, so the divergence is spelling (parsimony), not placement.
        assert_ne!(repeated, structured.cigar);
        assert!(is_left_aligned(&align(0, repeated), b"AAAC", b"AAAAC"));
        assert!(is_left_aligned(&structured, b"AAAC", b"AAAAC"));
    }

    #[test]
    fn the_output_always_consumes_exactly_the_read() {
        // A reconstruction round-trip invariant the leftmost oracle does NOT check: `shift_pass`
        // must never drop or add a read base. A surgery bug that changed read consumption could
        // still leave nothing shiftable and pass the property test, so pin it directly.
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
        let inputs: [(Vec<CigarOp>, &[u8], &[u8]); 4] = [
            (vec![Match(4), Deletion(1), Match(1)], b"GAAAT", b"GAAAAT"),
            (vec![Match(4), Insertion(1), Match(1)], b"GAAAAT", b"GAAAT"),
            (
                vec![
                    Match(1),
                    Deletion(1),
                    Deletion(1),
                    Insertion(1),
                    Match(1),
                    Match(1),
                ],
                b"AAAC",
                b"AAAAC",
            ),
            (
                vec![Insertion(1), Insertion(1), Insertion(1), Match(1), Match(1)],
                b"AAAAC",
                b"AC",
            ),
        ];
        for (index, (cigar, read, reference)) in inputs.iter().enumerate() {
            let out = normalized_cigar(0, cigar.clone(), read, reference);
            assert_eq!(
                consumed_read(&out),
                read.len(),
                "case {index}: output consumes the wrong number of read bases: {out:?}",
            );
        }
    }

    // ---- Boundary and degenerate inputs ---------------------------------------------

    #[test]
    fn an_empty_cigar_converges_and_is_untouched() {
        let mut alignment = align(0, vec![]);
        let report = RepeatedLeftAligner::new().left_align(&mut alignment, b"", b"");
        assert_eq!(report, ConvergenceReport::Converged { passes: 1 });
        assert!(alignment.cigar.is_empty());
    }

    #[test]
    #[should_panic(expected = "reference_offset")]
    fn an_offset_past_the_reference_panics() {
        let mut alignment = align(10, vec![Match(4), Deletion(1), Match(1)]);
        RepeatedLeftAligner::new().normalize(&mut alignment, b"GAAAT", b"GAAAAT");
    }

    #[test]
    fn an_already_leftmost_deletion_is_unchanged_and_converges_immediately() {
        let mut alignment = align(0, vec![Match(1), Deletion(1), Match(4)]);
        let report = RepeatedLeftAligner::new().left_align(&mut alignment, b"GAAAT", b"GAAAAT");
        assert_eq!(alignment.cigar, vec![Match(1), Deletion(1), Match(4)]);
        assert_eq!(report, ConvergenceReport::Converged { passes: 1 });
    }
}
