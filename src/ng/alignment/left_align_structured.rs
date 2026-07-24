//! Algorithm 1a — the **structured left-alignment pass**, an [`AlignmentNormalizer`].
//!
//! This is a **port, not a re-implementation**: it wraps production's
//! [`left_align_indels`](crate::pileup::walker::indel_norm) (itself a port of GATK's
//! `AlignmentUtils.leftAlignIndels`), and supplies only the trait wrapper. Nothing about
//! left-alignment is re-derived here — the plan's "reuse over rewrite" principle (arch §5, plan
//! step B1).
//!
//! # The two behaviours that make it "structured", and why they are the point
//!
//! A naive single left-to-right pass — "look at each indel, slide it one base left if the
//! neighbouring bases match" — is *assumed* to be what production does, and *assumed* to have one
//! failure: shifting one gap opens room for the next, and a single pass misses it. Production's pass
//! is not naive, and the two behaviours it has are exactly the ones that failure needs:
//!
//! - **it merges consecutive indels** — a deletion abutting an insertion (or a run of same-kind
//!   gaps) is collapsed into one canonical `D`-before-`I` event, so equal events spell the same way;
//! - **it propagates an indel across alignment blocks** rather than stopping at the first — the
//!   right-to-left walk carries a pending indel through a match block and keeps shifting into the
//!   next, so a gap slides as far left as the sequences allow, not just to the nearest boundary.
//!
//! Wrapping [`left_align_indels`] inherits both. Writing a fresh single pass here would reintroduce
//! precisely the bug those two behaviours exist to prevent, which is why this step is a wrapper and
//! its test asserts **byte-parity** against calling the production function directly.
//!
//! # It does not move `reference_offset`
//!
//! Production calls its left-aligner with end-deletion stripping **off**
//! (`remove_deletions_at_ends = false`): a deletion that rolls all the way to the read start stays a
//! first-op `Deletion`, and the alignment's start is left where it is. So this normalizer never
//! moves [`Alignment::reference_offset`]. The whole-`Alignment` signature exists because a normalizer
//! *may* need to (arch §3/§5); this one, faithfully to production, does not. A leading `Deletion` is
//! already leftmost — it has no aligned column to its left to slide across — so the property oracle
//! accepts it.
//!
//! # Isolated because its failure is silent
//!
//! A misplaced indel is a **wrong variant, not a crash**. So 1a lands as its own commit with the
//! property test green and byte-parity established, so that if indel placement ever moves, `git
//! bisect` lands here rather than in a bundle (plan step B1).
//!
//! # Algorithm 1c also lives here
//!
//! [`FixpointLeftAligner`] (1c) is a thin wrapper that applies 1a to a fixpoint, with an iteration
//! cap as a safety net that **fails loudly** (panics) rather than returning a half-normalised
//! alignment. It shares this file with 1a because it *is* 1a in a loop — no shifting of its own
//! (arch §Module home). It is the deliberate contrast with algorithm 1b, which *reports* cap
//! exhaustion instead of panicking on it.

use super::{Alignment, AlignmentNormalizer};
// A behavioural back-reference into the pipeline-stage module `pileup::walker`, deliberate under
// the plan's reuse-over-rewrite rule (arch §5): 1a *is* production's `left_align_indels`. It is the
// same kind of debt already recorded on `Alignment` for the `CigarOp` reuse — lifting `indel_norm`
// to a caller-agnostic peer is a production edit, and production is frozen; the port-back of this
// module is the moment to do it. `left_align_indels` is `pub(crate)`, so nothing leaks publicly.
use crate::pileup::walker::indel_norm::left_align_indels;

/// Algorithm 1a: left-align every indel in an alignment to its leftmost equivalent position with
/// **one structured pass**, by wrapping production's [`left_align_indels`].
///
/// Stateless — the pass reads only its arguments, so one value normalizes every read. Constructed
/// with no configuration; the sample group's error rate and the like are irrelevant to
/// normalization, which shifts gaps by matching bases alone (it is quality-blind).
#[derive(Debug, Clone, Copy, Default)]
pub struct StructuredLeftAligner;

impl AlignmentNormalizer for StructuredLeftAligner {
    /// Left-align `alignment.cigar` in place, leaving `alignment.reference_offset` untouched.
    ///
    /// `read` is the full read sequence; `reference` is the whole reference stretch, and
    /// `reference_offset` is where the placement starts within it — so the reference the production
    /// function needs (its `ref_seq[0]` at the read's first aligned base) is `reference` from the
    /// offset onward.
    ///
    /// **Fast path for the common case.** A read with no indel — the overwhelming majority — costs a
    /// single cigar scan and returns: `left_align_indels` early-returns before touching the
    /// reference, allocating nothing. So the caller does **not** need to pre-filter no-indel reads;
    /// this is already essentially free on them.
    ///
    /// **Two preconditions, handled differently.** The cigar's consistency with `read` (its
    /// read-consuming ops summing to `read.len()`) is upheld *gracefully* — production's own guard
    /// leaves a malformed cigar untouched rather than corrupting it. But `reference_offset <=
    /// reference.len()` is upheld by this call: a larger offset slices out of bounds and **panics**,
    /// which is the intended fail-fast (a silently mis-sliced reference would be a wrong variant).
    /// The debug assertion names the precondition so a dev/test build fails with the invariant it
    /// broke rather than an opaque slice-index message; the slice's own bounds check is what still
    /// panics in release. Anything a real placement satisfies.
    fn normalize(&self, alignment: &mut Alignment, read: &[u8], reference: &[u8]) {
        let offset = alignment.reference_offset as usize;
        debug_assert!(
            offset <= reference.len(),
            "reference_offset {offset} exceeds reference length {}",
            reference.len(),
        );
        let reference_from_placement = &reference[offset..];
        left_align_indels(&mut alignment.cigar, read, reference_from_placement);
    }
}

/// The most times [`FixpointLeftAligner`] applies its inner normalizer before giving up.
///
/// A correct [`StructuredLeftAligner`] reaches its fixpoint in **at most two** applications — one to
/// left-align, one to confirm nothing moved (`left_align_indels` is idempotent on already-leftmost
/// input). The cap sits well above that as a safety net: reaching it means the inner normalizer is
/// **not converging**, which is a defect, and 1c panics rather than return a half-normalised
/// alignment.
pub const FIXPOINT_MAX_ITERATIONS: u32 = 8;

// Compile-time guard: the cap must leave room for the two applications a correct 1a needs (one to
// left-align, one to confirm). Lowering it below that would panic on every real shifting read — a
// regression the runtime tests, which use explicit small caps, would not all catch.
const _: () = assert!(
    FIXPOINT_MAX_ITERATIONS >= 2,
    "FIXPOINT_MAX_ITERATIONS must allow the two applications a one-pass fixpoint needs"
);

/// Algorithm 1c: apply algorithm 1a repeatedly until the alignment stops changing, with the
/// iteration cap as a safety net that **fails loudly**.
///
/// It is a **thin wrapper over 1a**, not a third implementation of the shifting: [`drive_to_fixpoint`]
/// does no left-alignment of its own — it only re-applies [`StructuredLeftAligner`] and watches for
/// a fixpoint. Its whole reason to exist is what it does when the fixpoint is *not* reached within
/// the cap: it **panics**, rather than returning a half-normalised alignment the way freebayes'
/// capped loop silently does (spec §6). That is the deliberate contrast with algorithm 1b, which
/// *reports* exhaustion as a return value; 1c *fails loudly* on it. Because 1a is a proper one-pass
/// fixpoint, 1c converges in two iterations on every real input and never trips the cap — the panic
/// is a guard against a broken inner normalizer, not an expected outcome.
#[derive(Debug, Clone, Copy, Default)]
pub struct FixpointLeftAligner;

impl AlignmentNormalizer for FixpointLeftAligner {
    /// Left-align to a fixpoint by re-applying 1a. **Panics** if the alignment has not stopped
    /// changing after [`FIXPOINT_MAX_ITERATIONS`] applications — see the type docs.
    fn normalize(&self, alignment: &mut Alignment, read: &[u8], reference: &[u8]) {
        drive_to_fixpoint(
            &StructuredLeftAligner,
            alignment,
            read,
            reference,
            FIXPOINT_MAX_ITERATIONS,
        );
    }
}

/// Apply `inner` to `alignment` until it stops changing, up to `max_iterations` times; **panic** if
/// it has not converged by then.
///
/// This is the fixpoint machinery 1c is built from, factored out and generic so the cap-exhaustion
/// path can be driven in tests with a normalizer that never stabilises — the real inner (1a) never
/// reaches it. It performs **no** left-alignment itself: it clones the whole [`Alignment`] (so a
/// change in `reference_offset`, not just the cigar, counts as movement), re-applies `inner`, and
/// stops the instant an application changes nothing. Failing to converge is a **loud** panic, never
/// a quietly half-normalised return.
fn drive_to_fixpoint<N: AlignmentNormalizer>(
    inner: &N,
    alignment: &mut Alignment,
    read: &[u8],
    reference: &[u8],
    max_iterations: u32,
) {
    for _ in 0..max_iterations {
        let before = alignment.clone();
        inner.normalize(alignment, read, reference);
        if *alignment == before {
            return;
        }
    }
    panic!(
        "left-alignment did not reach a fixpoint within {max_iterations} iterations — the inner \
         normalizer is not converging; refusing to return a half-normalised alignment\n  \
         alignment: {:?}\n  read: {:?}\n  reference: {:?}",
        alignment,
        std::str::from_utf8(read).unwrap_or("<non-utf8>"),
        std::str::from_utf8(reference).unwrap_or("<non-utf8>"),
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::alignment::leftmost_property::{assert_left_aligned, is_left_aligned};
    use crate::pileup::walker::CigarOp::{self, Deletion, Insertion, Match};

    fn align(reference_offset: u64, cigar: Vec<CigarOp>) -> Alignment {
        Alignment {
            reference_offset,
            cigar,
        }
    }

    /// Normalize a copy and return its cigar — the common shape for the tests below.
    fn normalized_cigar(
        reference_offset: u64,
        cigar: Vec<CigarOp>,
        read: &[u8],
        reference: &[u8],
    ) -> Vec<CigarOp> {
        let mut alignment = align(reference_offset, cigar);
        StructuredLeftAligner.normalize(&mut alignment, read, reference);
        alignment.cigar
    }

    /// One table-driven fixture: an input placement (offset + cigar) and the two sequences it
    /// relates. Shared by the two table tests below.
    struct Case {
        offset: u64,
        cigar: Vec<CigarOp>,
        read: &'static [u8],
        reference: &'static [u8],
    }

    #[test]
    fn a_pure_match_alignment_is_returned_unchanged() {
        let cigar = normalized_cigar(0, vec![Match(5)], b"ACGTA", b"ACGTA");
        assert_eq!(cigar, vec![Match(5)]);
    }

    #[test]
    fn a_homopolymer_deletion_is_left_aligned() {
        // ref GAAAAT, read GAAAT: the rightmost-A deletion (4M1D1M) rolls to the leftmost A.
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
        // ref GAAAT, read GAAAAT: the rightmost inserted A (4M1I1M) rolls to the leftmost.
        let cigar = normalized_cigar(
            0,
            vec![Match(4), Insertion(1), Match(1)],
            b"GAAAAT",
            b"GAAAT",
        );
        assert_eq!(cigar, vec![Match(1), Insertion(1), Match(4)]);
    }

    #[test]
    fn two_placements_of_the_same_deletion_converge() {
        // The recall fix: the same one-A deletion placed at opposite ends of the homopolymer must
        // canonicalise to the SAME cigar, so the cohort merge buckets both reads onto one allele.
        let right = normalized_cigar(
            0,
            vec![Match(5), Deletion(1), Match(1)],
            b"GAAAAC",
            b"GAAAAAC",
        );
        let left = normalized_cigar(
            0,
            vec![Match(1), Deletion(1), Match(5)],
            b"GAAAAC",
            b"GAAAAAC",
        );
        assert_eq!(right, left);
        assert_eq!(right, vec![Match(1), Deletion(1), Match(5)]);
    }

    #[test]
    fn the_reference_offset_is_left_untouched() {
        // With a non-zero offset the placement starts partway into the stretch. 1a uses
        // remove_deletions_at_ends = false, so a rolled deletion stays a first-op Deletion and the
        // offset never moves.
        let mut alignment = align(2, vec![Match(4), Deletion(1), Match(1)]);
        StructuredLeftAligner.normalize(&mut alignment, b"GAAAT", b"TTGAAAAT");
        assert_eq!(alignment.reference_offset, 2);
        assert_eq!(alignment.cigar, vec![Match(1), Deletion(1), Match(4)]);
    }

    #[test]
    fn an_already_leftmost_deletion_comes_back_byte_identical() {
        // 1M1D4M on GAAAT/GAAAAT is already leftmost; it must return *unchanged*, not be
        // re-canonicalised to a different (also leftmost) spelling — idempotence is what keeps the
        // cohort merge's buckets stable. The property test alone would not catch a re-spelling.
        let cigar = normalized_cigar(
            0,
            vec![Match(1), Deletion(1), Match(4)],
            b"GAAAT",
            b"GAAAAT",
        );
        assert_eq!(cigar, vec![Match(1), Deletion(1), Match(4)]);
    }

    #[test]
    fn a_deletion_rolling_to_the_read_start_keeps_a_nonzero_reference_offset() {
        // reference[3..] = "AAAAT", read "AAAT": the deletion rolls all the way to a leading
        // first-op Deletion. With remove_deletions_at_ends = false, reference_offset must stay at 3
        // even though the deletion is now the first op — the crux of "1a never moves the offset".
        let mut alignment = align(3, vec![Match(3), Deletion(1), Match(1)]);
        StructuredLeftAligner.normalize(&mut alignment, b"AAAT", b"TTTAAAAT");
        assert_eq!(alignment.reference_offset, 3);
        assert_eq!(alignment.cigar, vec![Deletion(1), Match(4)]);
    }

    #[test]
    fn two_deletions_across_a_block_merge_and_left_align() {
        // 1M1D1M1D2M on AAAAAT/AAAT: two single-A deletions separated by a Match block. Left-
        // aligning carries the right deletion *across* the intervening block and merges both into
        // one leading Deletion(2) — the two load-bearing behaviours (cross-block propagation and
        // merge) exercised together, not just inherited by assumption.
        let mut alignment = align(
            0,
            vec![Match(1), Deletion(1), Match(1), Deletion(1), Match(2)],
        );
        StructuredLeftAligner.normalize(&mut alignment, b"AAAT", b"AAAAAT");
        assert_eq!(alignment.cigar, vec![Deletion(2), Match(4)]);
        assert!(is_left_aligned(&alignment, b"AAAT", b"AAAAAT"));
    }

    #[test]
    fn an_empty_cigar_is_left_untouched() {
        // The smallest legal input: empty read, empty cigar. Must no-op, not panic or mutate.
        let mut alignment = align(0, vec![]);
        StructuredLeftAligner.normalize(&mut alignment, b"", b"");
        assert_eq!(alignment.reference_offset, 0);
        assert!(alignment.cigar.is_empty());
    }

    #[test]
    fn an_offset_equal_to_the_reference_length_does_not_panic() {
        // Boundary: reference_offset == reference.len() yields an empty reference-from-placement
        // slice, which is legal (not out of bounds). With no indel the production function
        // early-returns; nothing should panic.
        let cigar = normalized_cigar(6, vec![], b"", b"GAAAAT");
        assert!(cigar.is_empty());
    }

    #[test]
    #[should_panic(expected = "reference_offset")]
    fn an_offset_past_the_reference_panics() {
        // Precondition reference_offset <= reference.len() violated. Fail-fast (panic) is the
        // intended failure mode — a silently clamped slice would be a wrong variant. The panic is
        // real in release too (the slice's own bounds check), not only via the debug assertion.
        let mut alignment = align(10, vec![Match(4), Deletion(1), Match(1)]);
        StructuredLeftAligner.normalize(&mut alignment, b"GAAAT", b"GAAAAT");
    }

    // ---- The property oracle grades the port (plan step B1; spec §6) -----------------
    //
    // The A2 property checker outranks agreement between the normalizers: wherever the port's
    // output is not leftmost, that is a finding. Run it on outputs spanning homopolymer,
    // period-2 SSR (deletion and insertion), and the already-leftmost / unique-context cases.

    #[test]
    fn every_normalized_output_passes_the_leftmost_property() {
        let cases = [
            // Homopolymer deletion, rightmost placement.
            Case {
                offset: 0,
                cigar: vec![Match(4), Deletion(1), Match(1)],
                read: b"GAAAT",
                reference: b"GAAAAT",
            },
            // Homopolymer insertion, rightmost placement.
            Case {
                offset: 0,
                cigar: vec![Match(4), Insertion(1), Match(1)],
                read: b"GAAAAT",
                reference: b"GAAAT",
            },
            // Period-2 SSR deletion (one AT unit), rightmost placement.
            Case {
                offset: 0,
                cigar: vec![Match(5), Deletion(2)],
                read: b"CATAT",
                reference: b"CATATAT",
            },
            // Period-2 SSR insertion (one CA unit) inside a longer repeat.
            Case {
                offset: 0,
                cigar: vec![Match(6), Insertion(2), Match(1)],
                read: b"GCACACACT",
                reference: b"GCACACT",
            },
            // Already leftmost, unique context — must stay put and still pass.
            Case {
                offset: 0,
                cigar: vec![Match(2), Deletion(1), Match(2)],
                read: b"GATT",
                reference: b"GACTT",
            },
            // Non-zero offset.
            Case {
                offset: 2,
                cigar: vec![Match(4), Deletion(1), Match(1)],
                read: b"GAAAT",
                reference: b"TTGAAAAT",
            },
        ];
        for (index, case) in cases.iter().enumerate() {
            let mut alignment = align(case.offset, case.cigar.clone());
            StructuredLeftAligner.normalize(&mut alignment, case.read, case.reference);
            assert!(
                is_left_aligned(&alignment, case.read, case.reference),
                "case {index}: port output is not leftmost: {:?}",
                alignment.cigar,
            );
            // Also exercise the panicking form, so a failure names the offender.
            assert_left_aligned(&alignment, case.read, case.reference);
        }
    }

    // ---- Byte-parity against the production function it ports (the port anchor) -------
    //
    // The wrapper must produce exactly what `left_align_indels` produces on the same read and the
    // same reference-from-offset. This is what proves the two load-bearing behaviours (merge +
    // cross-block propagation) are inherited rather than reimplemented, and it fails loudly if the
    // read/reference arguments are ever transposed.

    #[test]
    fn output_matches_left_align_indels_byte_for_byte() {
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
            // A colliding indel run: I-before-D that production reorders to D-before-I and merges.
            // This is the "merge consecutive indels" behaviour — inherited, not reimplemented.
            Case {
                offset: 0,
                cigar: vec![Match(2), Insertion(1), Deletion(1), Match(2)],
                read: b"GATT",
                reference: b"GACTT",
            },
            // Deletion rolling to the read start stays a leading Deletion (offset kept).
            Case {
                offset: 0,
                cigar: vec![Match(3), Deletion(1), Match(1)],
                read: b"AAAT",
                reference: b"AAAAT",
            },
            // Non-zero offset: parity must use reference-from-offset.
            Case {
                offset: 2,
                cigar: vec![Match(4), Deletion(1), Match(1)],
                read: b"GAAAT",
                reference: b"TTGAAAAT",
            },
            // Cross-block propagation + merge: two deletions separated by a Match block.
            Case {
                offset: 0,
                cigar: vec![Match(1), Deletion(1), Match(1), Deletion(1), Match(2)],
                read: b"AAAT",
                reference: b"AAAAAT",
            },
        ];
        for (index, case) in cases.iter().enumerate() {
            let normalized =
                normalized_cigar(case.offset, case.cigar.clone(), case.read, case.reference);

            let mut expected = case.cigar.clone();
            left_align_indels(
                &mut expected,
                case.read,
                &case.reference[case.offset as usize..],
            );

            assert_eq!(
                normalized, expected,
                "case {index}: wrapper diverged from left_align_indels",
            );
        }
    }

    // ---- Algorithm 1c: the fixpoint wrapper (plan step C2) ----------------------------

    /// A stand-in normalizer that **never** reaches a fixpoint — it changes the alignment on every
    /// application (by advancing the offset). Used to drive [`drive_to_fixpoint`]'s cap-exhaustion panic
    /// deliberately; the real inner (1a) converges and never reaches that path.
    struct NonConvergingNormalizer;
    impl AlignmentNormalizer for NonConvergingNormalizer {
        fn normalize(&self, alignment: &mut Alignment, _read: &[u8], _reference: &[u8]) {
            alignment.reference_offset += 1;
        }
    }

    #[test]
    fn the_fixpoint_wrapper_left_aligns_and_stays_leftmost() {
        let mut alignment = align(0, vec![Match(4), Deletion(1), Match(1)]);
        FixpointLeftAligner.normalize(&mut alignment, b"GAAAT", b"GAAAAT");
        assert_eq!(alignment.cigar, vec![Match(1), Deletion(1), Match(4)]);
        assert!(is_left_aligned(&alignment, b"GAAAT", b"GAAAAT"));
        assert_left_aligned(&alignment, b"GAAAT", b"GAAAAT");
    }

    #[test]
    fn an_already_leftmost_input_reaches_the_fixpoint_without_shifting() {
        // The one-iteration path: an input already leftmost is a fixpoint on the first application,
        // so 1c returns it unchanged and does not panic.
        let mut alignment = align(0, vec![Match(1), Deletion(1), Match(4)]);
        FixpointLeftAligner.normalize(&mut alignment, b"GAAAT", b"GAAAAT");
        assert_eq!(alignment.cigar, vec![Match(1), Deletion(1), Match(4)]);
    }

    #[test]
    fn the_fixpoint_wrapper_agrees_with_the_single_pass() {
        // 1c is 1a in a loop, and 1a is already a fixpoint, so the two produce identical output.
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
                offset: 2,
                cigar: vec![Match(4), Deletion(1), Match(1)],
                read: b"GAAAT",
                reference: b"TTGAAAAT",
            },
        ];
        for (index, case) in cases.iter().enumerate() {
            let mut fixpoint = align(case.offset, case.cigar.clone());
            FixpointLeftAligner.normalize(&mut fixpoint, case.read, case.reference);

            let mut structured = align(case.offset, case.cigar.clone());
            StructuredLeftAligner.normalize(&mut structured, case.read, case.reference);

            assert_eq!(fixpoint, structured, "case {index}: 1c and 1a disagree");
        }
    }

    #[test]
    fn a_correct_inner_converges_within_two_iterations() {
        // 1a needs exactly two iterations (one to left-align, one to confirm no change), so a cap of
        // two is enough and must NOT panic.
        let mut alignment = align(0, vec![Match(4), Deletion(1), Match(1)]);
        drive_to_fixpoint(
            &StructuredLeftAligner,
            &mut alignment,
            b"GAAAT",
            b"GAAAAT",
            2,
        );
        assert_eq!(alignment.cigar, vec![Match(1), Deletion(1), Match(4)]);
    }

    #[test]
    #[should_panic(expected = "did not reach a fixpoint")]
    fn a_deliberately_capped_run_surfaces_the_failure_rather_than_swallowing_it() {
        // The whole point of 1c: a run that cannot converge within its cap **panics**, it does not
        // return a half-normalised alignment. Driven with an inner that never stabilises.
        let mut alignment = align(0, vec![Match(3)]);
        drive_to_fixpoint(
            &NonConvergingNormalizer,
            &mut alignment,
            b"ACG",
            b"ACG",
            FIXPOINT_MAX_ITERATIONS,
        );
    }

    #[test]
    #[should_panic(expected = "did not reach a fixpoint")]
    fn a_cap_of_one_cannot_confirm_convergence_and_panics() {
        // A single iteration left-aligns but cannot confirm stability, so the safety net fires even
        // with the real, correct 1a — the cap guards the confirmation, not just non-convergence.
        let mut alignment = align(0, vec![Match(4), Deletion(1), Match(1)]);
        drive_to_fixpoint(
            &StructuredLeftAligner,
            &mut alignment,
            b"GAAAT",
            b"GAAAAT",
            1,
        );
    }
}
