# ng alignment — normalization (plan 3 of 3), Milestone A

**Date:** 2026-07-24
**Plan:** [alignment_normalization.md](../../ng/impl_plan/alignment_normalization.md) — steps A1, A2.
**Spec/arch:** [spec §6/§7](../../ng/spec/alignment.md), [arch §3/§5/§Test & bench shape](../../ng/arch/alignment.md).
**Status:** implemented → reviewed → fixes-applied (per step). At **Checkpoint A**.

---

## 1. Plan

Milestone A lands the **interface and the oracle, before any implementation** (plan Principle:
"the oracle before the implementations"): the `AlignmentNormalizer` trait (A1) and the leftmost
property checker that grades the three normalizers against a definition (A2). No normalizer exists
yet — that is deliberate: the oracle is built and tested first so it grades every implementation
from its first commit and **outranks agreement** between the three (spec §6).

## 2. Assumptions / recorded deviations

- **`AlignmentNormalizer: Sized`** (A1) — the arch §3 sketch shows a plain trait, but the shipped
  sibling traits (`BestPathAligner`, `MarginalAligner`) both carry `: Sized` for the module's
  decided static-dispatch-only rule (spec §7, arch §4). Added it for consistency, with a parallel
  rationale doc comment. Recorded, not a design change.
- **A1 ships no `# Examples` doc test** — the two sibling aligner traits carry none either
  (prose + `#[cfg(test)]` anchors is the module pattern), and A1 has no implementor to write one
  against. Recorded exemption (review Mi2).
- **The oracle is derived from the definition, not from `normalize_alleles`** (A2) — building it on
  production's left-aligner would grade algorithm 1a against itself. The shift condition (crossed
  column must be a real match; far edge preserved) is written straight from the definition, and it
  independently reproduces the two predicates production's `normalize_alleles` tests — recorded as
  a *result*, not an assumption.

## 3. Changes made

- **A1** — `pub trait AlignmentNormalizer: Sized { fn normalize(&self, &mut Alignment, read: &[u8],
  reference: &[u8]); }` in [alignment/mod.rs](../../../../src/ng/alignment/mod.rs). Takes the whole
  `Alignment` (so a leading-deletion strip can move `reference_offset`); read/reference are bare
  `&[u8]` (quality-blind). No production impl; two `#[cfg(test)]` compile anchors
  (`IdentityNormalizer`, `LeadingDeletionStripper`) + a generic `normalize_via` driver.
- **A2** — new `#[cfg(test)]` module
  [alignment/leftmost_property.rs](../../../../src/ng/alignment/leftmost_property.rs): `IndelKind`,
  `ShiftableIndel`, `find_shiftable_indel`, `is_left_aligned`, `assert_left_aligned`, and
  well-formedness helpers. Registered in mod.rs.

## 4. Tests

- A1: `alignment_normalizer_leaves_an_already_leftmost_alignment_untouched`,
  `alignment_normalizer_can_move_the_reference_offset`, plus the two review-added guard tests
  (non-leading-deletion, empty-cigar), all driven through the generic `normalize_via`.
- A2: 20 unit tests — homopolymer + period-2 shiftable/leftmost (deletion and insertion), the
  match-required clause (two fixtures), and cursor/boundary handling (`reference_offset`, soft
  clip, `Skip`, adjacent indels, empty cigar, trailing insertion, chained reference-advance). Every
  discriminating fixture was hand-traced to confirm it fails under the specific mutation it targets.

## 5. Validation

- `cargo fmt --check` → clean · `cargo clippy --all-targets --all-features -- -D warnings` → clean ·
  `cargo test --lib` → 2292 passed / 0 failed / 4 ignored (container, `./scripts/dev.sh`).

## 6. Reviews and fixes

- A1 review [ng_normalizer_a1](../reviews/ng_normalizer_a1_2026-07-24.md): 0 Bl / 0 Maj / 2 Min / 2 Nit; fixes [fixes_applied_ng_normalizer_a1](../reviews/fixes_applied_ng_normalizer_a1_2026-07-24.md).
- A2 review [ng_normalizer_a2](../reviews/ng_normalizer_a2_2026-07-24.md): **1 Bl / 4 Maj** / 2 Min / 3 Nit — every Blocker/Major was a *test that could not fail* (the recurring ng lesson), all applied with hand-verified discriminating fixtures; fixes [fixes_applied_ng_normalizer_a2](../reviews/fixes_applied_ng_normalizer_a2_2026-07-24.md).

## 7. Tradeoffs and follow-ups

- The oracle is `#[cfg(test)]` and ships no production code; it is the grading harness for
  Milestones B–D.
- Milestone B (algorithm 1a — the port) is next after the Checkpoint A pause.
