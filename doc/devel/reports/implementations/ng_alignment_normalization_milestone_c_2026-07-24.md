# ng alignment — normalization (plan 3 of 3), Milestone C

**Date:** 2026-07-24
**Plan:** [alignment_normalization.md](../../ng/impl_plan/alignment_normalization.md) — steps C1, C2.
**Spec/arch:** [spec §6/§8](../../ng/spec/alignment.md), [arch §5/§Module home](../../ng/arch/alignment.md).
**Status:** implemented → reviewed → fixes-applied (per step). At **Checkpoint C**.

---

## 1. Plan

Milestone C lands the two remaining normalizers, both contrasting with 1a on how they treat the
iteration cap: **1b reports exhaustion; 1c fails loudly on it.**

- **C1 — algorithm 1b — `RepeatedLeftAligner`** ([left_align_repeated.rs](../../../../src/ng/alignment/left_align_repeated.rs)): freebayes' shape — repeated *simple* one-base shift passes until nothing moves, capped (`MAX_PASSES = 20`). A genuinely **independent** left-aligner (leans on neither 1a nor production), implemented from the description not by transliteration. **Exhaustion is reported** as a `#[must_use] ConvergenceReport`, the fact freebayes computes then drops.
- **C2 — algorithm 1c — `FixpointLeftAligner`** (in [left_align_structured.rs](../../../../src/ng/alignment/left_align_structured.rs), beside 1a per arch §Module home): 1a applied to a fixpoint via the generic `drive_to_fixpoint`, with `FIXPOINT_MAX_ITERATIONS = 8` as a safety net that **panics** rather than return a half-normalised alignment. A thin wrapper doing no shifting itself.

## 2. Assumptions / recorded deviations

- **C1: each pass slides one base** (freebayes slides maximally). This makes the cap bound total shift *distance* and the exhaustion path reachable/testable. Recorded in the module doc.
- **C1: 1b does not trim a complex D/I overlap** (1a does, via `normalize_alleles`). After canonicalization, this is the *only* residual 1a-vs-1b difference; both spellings remain leftmost (the oracle grades placement, not parsimony). Full trim reimplements `normalize_alleles`' trim+shift interleaving — substantial, risky, beyond C1's shift+cap+report scope, and rare on real reads. Corrected the doc, pinned the difference with a test, and **↓ flagged for owner decision (see §7).**
- **C2: `FIXPOINT_MAX_ITERATIONS = 8`** — a correct 1a is a one-pass idempotent fixpoint (converges in 2 iterations), so 8 is generous headroom; hitting the cap means the inner is broken. A compile-time assert pins the floor.

## 3. Changes made

- C1: new `left_align_repeated.rs` — `RepeatedLeftAligner`, `ConvergenceReport`, `MAX_PASSES`; independent `shift_pass`/`shiftable_indices`/`canonicalize` + small duplicated CIGAR helpers (production's are private + frozen).
- C2: `FixpointLeftAligner`, the generic `drive_to_fixpoint`, `FIXPOINT_MAX_ITERATIONS` added to `left_align_structured.rs`.

## 4. Tests

- C1 (17): homopolymer/period-2 left-alignment, pass-count grows with distance, exhaustion reported at a too-small cap (result non-leftmost), cap boundaries, `max_passes=0`, the A2 property oracle on converged output, agreement with 1a (incl. a non-canonical input), a complex-indel documented difference, a read-consumption round-trip invariant. (A reviewer's 50k-case fuzz confirmed no corruption + genuine independence.)
- C2 (5 for 1c): 1c left-aligns and stays leftmost (A2), agrees with 1a, converges in two iterations, an already-leftmost one-iteration path, and **two `#[should_panic]` tests** — a deliberately-capped non-converging run and a cap-of-one confirmation boundary — both confirmed genuine.

## 5. Validation

- `cargo fmt --check` clean · `cargo clippy --all-targets --all-features -- -D warnings` clean · `cargo test --lib` → 2328 passed / 4 ignored (container).

## 6. Reviews and fixes

- C1 review [ng_normalizer_c1](../reviews/ng_normalizer_c1_2026-07-24.md): 0 Bl / 3 Maj / 3 Min — `#[must_use]`, canonicalize-on-no-shift, and the documented-trim-difference; fixes [fixes_applied_ng_normalizer_c1](../reviews/fixes_applied_ng_normalizer_c1_2026-07-24.md).
- C2 review [ng_normalizer_c2](../reviews/ng_normalizer_c2_2026-07-24.md): 0 Bl / 0 Maj / 5 Min — fixpoint/panic logic and both `#[should_panic]` tests confirmed genuine; fixes [fixes_applied_ng_normalizer_c2](../reviews/fixes_applied_ng_normalizer_c2_2026-07-24.md).

## 7. Tradeoffs and follow-ups — **owner decision at Checkpoint C**

- **1b's complex-indel trim (C1 M3).** Should 1b also trim an overlapping D/I for full spelling parity with 1a, or is the documented difference acceptable? It matters for how the D1 differ-at-all screen (Milestone D) is read: an untrimmed complex indel counts as a 1a/1b disagreement — legitimately a "differently-spelled equivalent" the screen exists to surface, but attributable to 1b's omission rather than a placement ambiguity. Recommendation: leave as-is (rare on real reads; the screen surfaces it honestly), and revisit only if D1 shows complex-indel disagreements dominating. Not blocking; flagged for your call before D1.
