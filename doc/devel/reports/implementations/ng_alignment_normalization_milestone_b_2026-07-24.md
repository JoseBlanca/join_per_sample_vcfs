# ng alignment — normalization (plan 3 of 3), Milestone B

**Date:** 2026-07-24
**Plan:** [alignment_normalization.md](../../ng/impl_plan/alignment_normalization.md) — step B1.
**Spec/arch:** [spec §6](../../ng/spec/alignment.md), [arch §5](../../ng/arch/alignment.md).
**Status:** implemented → reviewed → fixes-applied. At **Checkpoint B**.

---

## 1. Plan

Milestone B lands **algorithm 1a — the structured left-alignment pass** ([`StructuredLeftAligner`](../../../../src/ng/alignment/left_align_structured.rs)), the first `AlignmentNormalizer` implementation. It is a **port**: it wraps production's `left_align_indels` and supplies only the trait wrapper (arch §5's reuse-over-rewrite). Its own commit, because a misplaced indel is a wrong variant not a crash — the failure is silent, so it is kept bisectable (plan B1).

## 2. Assumptions / recorded deviations

- **1a does not move `reference_offset`.** Production calls its left-aligner with `remove_deletions_at_ends = false`, so a deletion that rolls to the read start stays a first-op `Deletion` and the alignment start is fixed. The whole-`Alignment` signature exists for a normalizer that *may* move the offset; this one, faithfully to production, does not. Recorded.
- **Behavioural back-reference into `pileup::walker`.** Calling `left_align_indels` is a real dependency from the caller-agnostic `ng::alignment` into a pipeline-stage module — deliberate under reuse-over-rewrite, recorded as debt alongside the existing `CigarOp` note (production is frozen; lifting `indel_norm` is the port-back moment).
- **`StructuredLeftAligner` name** kept over `...Normalizer` — "left-aligner" is standard genomics vocabulary and the whole normalizer sub-family shares it (`RepeatedLeftAligner` for 1b).

## 3. Changes made

- New [left_align_structured.rs](../../../../src/ng/alignment/left_align_structured.rs): `pub struct StructuredLeftAligner` (stateless unit struct) implementing `AlignmentNormalizer` by slicing `&reference[reference_offset..]` and calling `left_align_indels(&mut cigar, read, ref_from_placement)`. A `debug_assert` names the `reference_offset <= reference.len()` precondition. Registered `pub mod` in mod.rs.

## 4. Tests (13)

- Behaviour: pure-match no-op, homopolymer deletion/insertion left-aligned, two-placement convergence (the recall fix), offset-untouched, idempotence on already-leftmost, leading-deletion-at-nonzero-offset, cross-block propagation + merge, empty cigar, offset==len boundary, offset>len panic.
- **The two verifications the plan mandates:** the A2 **property oracle** on every normalized output (`every_normalized_output_passes_the_leftmost_property` — the oracle and production agree, no divergence), and **byte-parity** of the resulting operations against calling `left_align_indels` directly (`output_matches_left_align_indels_byte_for_byte`, recomputing `expected` from the raw inputs so a transposition/dropped-offset diverges).

## 5. Validation

- `cargo fmt --check` clean · `cargo clippy --all-targets --all-features -- -D warnings` clean · `cargo test --lib` → 2305 passed / 4 ignored (container).

## 6. Review and fixes

- Review [ng_normalizer_b1](../reviews/ng_normalizer_b1_2026-07-24.md): 0 Bl / 2 Maj / 6 Min / 2 Nit — all test-hardening or diagnostics; both anchor oracles confirmed genuine (not vacuous). Fixes [fixes_applied_ng_normalizer_b1](../reviews/fixes_applied_ng_normalizer_b1_2026-07-24.md): 8 applied (6 new tests + `debug_assert` + debt note), 1 Nit hoisted, 1 Nit won't-fix (`as usize`, ng convention).

## 7. Tradeoffs and follow-ups

- Milestone C (1b — freebayes' repeated passes, its own file; and 1c — the fixpoint wrapper over 1a, in *this* file) is next after the Checkpoint B pause.
