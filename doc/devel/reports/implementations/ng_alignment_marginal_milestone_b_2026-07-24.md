# ng alignment — marginal aligners, Milestone B (algorithm 5, the sequence marginal)

**Date:** 2026-07-24
**Plan:** [alignment_marginal.md](../../ng/impl_plan/alignment_marginal.md) (plan 2 of 3), steps B1–B3, at **Checkpoint B**.
**Design authority:** spec [alignment.md](../../ng/spec/alignment.md) §5.1, §7; arch [alignment.md](../../ng/arch/alignment.md) §3, §5; reuse target production's `align_subst` ([pair_hmm.rs](../../../../src/ssr/cohort/pair_hmm.rs)).
**Method:** plan-driven — one implement → review → apply-fixes → commit loop per step; **B2 and B3 are own-commit "silent-failure" steps, not bundled** (plan mandate). Per-category review outputs left in gitignored `tmp/review_2026-07-24_ng-marginal-b{1,2,3}/`.

---

## 1. What landed

Algorithm 5 — `SsrSequenceMarginal` in [ssr_marginal_sequence.rs](../../../../src/ng/alignment/ssr_marginal_sequence.rs), the first `MarginalAligner` implementation: the probability that one sequence produced another, summed over line-ups, gaps confined to the ends. A **faithful port of production's `align_subst` — this function alone**, not the `HipstrModel` stutter half (which is the genotyping's; spec §5.1 boundary note). Built in three own-commit steps:

- **B1** (`2a7f217`) — the equal-length case: one line-up, base-by-base under a flat rate ε; all-match `(1−ε)^len` fast path + the substitution product. Checked constructor `try_new` holding ε.
- **B2** (`fd551cd`) — the unequal-length banded forward + `linear_probability` dispatcher, gaps confined within `FLANK_SLOP`(=2) of either end. The **interior-gap restriction is load-bearing** (keeps base error and stutter separately identifiable).
- **B3** (this commit) — the logarithm boundary: `impl MarginalAligner`, one `.ln()` turning the linear result into `LogProb`, `-∞` for an unreachable line-up. The normalisation defect recorded in the impl doc.

## 2. The one design reconciliation (recorded, not a stop-and-ask)

**Algorithm 5 works in linear space and holds ε directly, not the module's log-space `FlatEmission` component** — a departure from the plan's precondition note ("reuses the Emission component"), driven by the design authority the plan itself cites. Required because (a) B2's forward *sums* probabilities, cheap in linear space as in the source (log-space would need `logsumexp` and round differently), and (b) B3's parity oracle is production's *linear* `align_subst`; reusing `FlatEmission` would force `exp()` round-trips that break byte-parity. The tie to the emission component is the shared validation (`DomainError::ErrorRate`, the same check `FlatEmission::try_new` makes) and the shared flat-rate concept. This is "adapt to the reused API's real shape" latitude, recorded in the module doc and the B1 commit.

## 3. Verification — the port is byte-identical to production

The load-bearing result: **ng's linear scorer is bit-for-bit identical to production's `align_subst`.** Pinned three ways in B3:
- `linear_probability_matches_production_align_subst_bit_for_bit` — all three arms (exact, equal-length substitution, insertion, deletion, a length-3 difference) × a sweep of ε including the endpoints 0.0 and 1.0.
- `linear_probability_matches_align_subst_on_random_sequences` — a **proptest** over random byte strings × ε in `[0,1]`, guarding the band/flank edges the example cases might all agree on.
- `the_marginal_is_the_logarithm_of_production_align_subst` — end to end: `marginal_probability(...).get() == align_subst(...).ln()`.

Plus the interior-gap and boundary guards: an interior indel scores ≪ an end gap; the band floor `|m−n|` admits a difference beyond `FLANK_SLOP` (drop it and larger stutter silently returns 0); the forward is symmetric (drives the deletion side production's own caller never reaches); one `.ln()` (not zero, not two); an unreachable line-up is `-∞` not 0.

## 4. Review outcomes (per step)

- **B1** (4 sub-agents / 8 checklists): 0 Blocker / 0 Major shipped; parity confirmed sound. Applied 2 tests (ε endpoints degenerate without flooring; fast-path == substitution product). Adapted (not applied): a `#[should_panic]` on the equal-length `debug_assert!` was declined — it passes in the test build while proving nothing about release (the project's recorded false-confidence trap); the real guarantee is B2's routing.
- **B2** (3 sub-agents / 7 checklists): no correctness bugs, port byte-faithful; both test-validity questions resolved positively. Applied 3 coverage tests (deletion-side symmetry; band-floor beyond the slop; empty-read boundary) + the `#[allow(dead_code)]` → `#[cfg_attr(not(test), expect(dead_code))]` **compiler backstop** (see §5).
- **B3** (3 sub-agents / 7 checklists): 0 Blocker / 0 Major; all four correctness questions verified. Applied: ε-sweep in the parity test, the end-to-end marginal-vs-`align_subst().ln()` test, the proptest, a dropped redundant assertion, and a doc-anchor fix.

Every finding across the three steps was Applied or Adapted-with-reason; none disputed as wrong, none deferred as design.

## 5. The dead-code backstop — a mechanical guarantee, discovered

B1–B2 built the "ported function" (`linear_probability` and its callees) before B3 wired the `MarginalAligner` impl, leaving it dead in the plain lib build. Rather than a manual "remove the `#[allow]` at B3" checklist item, B2 switched to `#[cfg_attr(not(test), expect(dead_code, reason=…))]` on **`linear_probability` alone** — the single entry point of the dead chain, since every callee counts as *used* through it. `expect` (not `allow`) fails the build via `unfulfilled_lint_expectations` the moment B3 makes the chain live; `cfg(not(test))` because the test build already uses it. **B3 confirmed it works**: wiring the trait forced removal of the attribute (the file now has zero dead-code suppressions), exactly the mechanical enforcement intended. (The reviewer's first suggestion, an unconditional `#[expect]`, was corrected — it is unfulfilled in the `--tests` build.)

## 6. Validation (container, via `./scripts/dev.sh`)

- `cargo fmt --check` → clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` → clean.
- `cargo test --lib` → **2247 passed, 0 failed**, 4 ignored.

The two pre-existing, unrelated red project-wide commands (`cargo test --all-targets` bench panic; `cargo doc --no-deps` intra-doc links) are the documented Standing items; untouched.

## 7. Follow-ups

- **Milestone C** — algorithm 6 (the whole-read forward, fixed synthetic qualities): the remaining `MarginalAligner` implementation, whose `Context` is the repeat geometry (the borrowed-context path A2's GAT exists for).
- The 5-vs-6 bake-off is **out of this plan** — it spans the genotyping (spec §10.3); Milestone C only proves each algorithm computes what it claims.

**Checkpoint B** — algorithm 5 matches the ported function on equal-length inputs, sums correctly over end gaps on unequal-length ones, refuses interior gaps, and returns logarithms. Hard pause for human review before Milestone C.
