//! ng alignment — lining a read up against a reference stretch.
//!
//! This is a **folder, not a file**, because it holds competing implementations that get
//! compared: the traits will live here in `mod.rs` and each algorithm sits in its own file
//! beside its siblings (`doc/devel/ng/arch/alignment.md` §Module home). File names say
//! which trait an algorithm implements and whether it is repeat-aware, because both are
//! load-bearing — the best-path and marginal families implement *different* traits, and
//! every repeat-aware algorithm has an affine counterpart it must not be confused with.
//!
//! **It is not a pipeline step, and it knows none.** Step 2 (read preparation) and step 7
//! (read likelihood) both call it; it knows about neither. It takes sequences and returns
//! alignments or probabilities, which is why it sits outside the step-per-module rule
//! rather than breaking it (`doc/devel/ng/spec/alignment.md` §1).
//!
//! Landed so far: [`Alignment`], the output shape of the general-purpose affine aligner.
//! No trait and no algorithm yet — those arrive with the later steps of the plan
//! `doc/devel/ng/impl_plan/alignment_best_path.md`.

use crate::pileup::walker::CigarOp;

/// A read placed against a reference stretch: where the placement starts, and the
/// operations that spell it out from there.
///
/// This is the output shape of the **affine** best-path aligner — the general-purpose
/// one. The repeat-aware aligners return where a read's repeat sits instead, which is a
/// different shape and a different type; that is why the aligner trait's output is an
/// associated type rather than this struct (arch §2.1, §3).
///
/// # Reusing production's `CigarOp`
///
/// The operations are production's [`CigarOp`], **reused rather than re-minted**. A
/// parallel operation type would buy nothing and cost a conversion at every boundary
/// where an alignment meets code that already speaks CIGAR — and unlike ng's [`Motif`],
/// which had to be ported because production's is `pub(crate)` and would have leaked out
/// of a `pub` signature, `CigarOp` is already fully public, so reuse here costs
/// production nothing (arch §5).
///
/// Worth naming, because the import path misleads: `CigarOp` lives inside
/// `pileup::walker`, so this module appears to depend on the pileup walker, which it does
/// not — it knows no callers at all. `CigarOp` is really crate-wide CIGAR vocabulary
/// (`bam/`, `ssr/pileup/`, `pop_var_caller/cli/`, the benches and ng all consume it), and
/// its home inside a pipeline stage is a **pre-existing misplacement**. Lifting it to a
/// top-level peer is a production edit, and production is frozen (owner, 2026-07-16), so
/// this is recorded as known debt rather than fixed here; porting this module back is the
/// natural moment to do it.
///
/// [`Motif`]: crate::ng::region_typing::segment_criteria::Motif
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    /// Where the placement starts, measured in bases from the start of **the reference
    /// stretch the aligner was given** — not a genome position. The aligner is handed a
    /// slice and knows nothing about where that slice came from, so it cannot return a
    /// coordinate; turning this into a genome position is the caller's job, and it needs
    /// the stretch's own start to do it.
    ///
    /// `u64`, not `usize`: ng speaks one width for coordinates and lengths, so nothing
    /// narrows and no off-by-width bug is possible ([`Bp`], and
    /// [`RepeatInterval`] — the same shape, an offset into the slice its
    /// producer was handed, `u64` for the same reason).
    ///
    /// [`Bp`]: crate::ng::types::Bp
    /// [`RepeatInterval`]: crate::ng::tandem_repeat::RepeatInterval
    pub reference_offset: u64,
    /// The operations from `reference_offset` onward, in read order.
    ///
    /// Named `cigar` rather than `ops` to match how every other site in the crate spells
    /// this concept (`cigar` / `cigar_ops`); the architecture doc's sketch says `ops`, but
    /// it states that its signatures are illustrative and the contract is what ships.
    ///
    /// **An affine aligner inhabits only part of [`CigarOp`].** The enum is production's
    /// *parser* vocabulary and must represent everything a BAM can contain, including
    /// `Skip`, `HardClip` and `Padding`, which describe an input record rather than a
    /// placement this module computed. Exactly which subset the affine aligner emits is
    /// **not settled here** — its algorithm lands in Milestone E and is currently gated,
    /// and arch §6 still has the aligner's output shape open. The obligation is recorded
    /// rather than guessed: the step that writes the producer states the subset, and only
    /// then may a consumer treat the remaining variants as unreachable.
    pub cigar: Vec<CigarOp>,
}

// A note on derives: `Hash` is not derived because [`CigarOp`] does not derive it, and
// adding it there is a production edit. If an `Alignment` ever needs to be a map key,
// that is the blocker to clear first.

#[cfg(test)]
mod tests {
    use super::*;

    /// A0 lands no algorithm, so this test's job is to be a **compile-time anchor**, and
    /// it is written to fail if the things it anchors change: that the type keeps both
    /// fields, derives structural equality over *both* of them, and still builds its
    /// operations out of production's [`CigarOp`].
    ///
    /// The whole-struct comparisons are the part that carries weight — a hand-written
    /// `PartialEq` that ignored either field would pass field-by-field assertions and
    /// fail these.
    ///
    /// TODO(Milestone E): the two contract properties stated on the fields have no test
    /// and cannot have one until a producer exists. The step that writes the affine
    /// aligner owes `align_returns_offset_relative_to_the_reference_stretch` (the
    /// stretch-versus-genome confusion is a wrong answer, not a panic) and
    /// `align_returns_operations_in_read_order`.
    #[test]
    fn alignment_distinguishes_offset_and_operation_order() {
        let alignment = Alignment {
            reference_offset: 7,
            cigar: vec![CigarOp::Match(10), CigarOp::Insertion(2), CigarOp::Match(5)],
        };

        assert_eq!(alignment.reference_offset, 7);
        assert_eq!(
            alignment.cigar,
            vec![CigarOp::Match(10), CigarOp::Insertion(2), CigarOp::Match(5)]
        );

        // Equality must see the offset: same operations, different placement.
        assert_ne!(
            alignment,
            Alignment {
                reference_offset: 8,
                cigar: alignment.cigar.clone(),
            }
        );

        // Equality must see operation order: same operations, permuted.
        assert_ne!(
            alignment,
            Alignment {
                reference_offset: 7,
                cigar: vec![CigarOp::Insertion(2), CigarOp::Match(10), CigarOp::Match(5)],
            }
        );
    }
}
