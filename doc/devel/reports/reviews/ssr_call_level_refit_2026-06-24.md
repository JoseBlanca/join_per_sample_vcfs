# Code Review: ssr-call per-locus stutter-rate refit (Step I2)

**Date:** 2026-06-24
**Reviewer:** rust-code-review skill (orchestrator; inline synthesis)
**Scope:** commit `fd53330` on branch `ssr-cohort` — the per-locus stutter-level (rate) multiplier in the genotype EM
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** PR diff, commit `fd53330` (Step I2).
- **In-scope files:** [em.rs](../../../../src/ssr/cohort/em.rs) — `attribute_locus` (replacing `attribute_locus_slips`) + `LocusSlipFit`; `refit_level_multiplier` + `LEVEL_MULT_MAX`; `compute_data_ll`'s `level_mult` parameter; the combined θ+level refit loop; `EmCfg` `level_shrink`/`level_tol`; 3 new/updated tests.
- **Out of scope:** I1 (reviewed `46df902`); `refine_theta_locus` / `SlipProfile` (B2); the cohort-level `reduce_level` (E1).
- **Categories considered:** reliability, errors, naming, idiomatic, defaults, smells, refactor_safety, extras (determinism, hot path).

## 2. Verdict

**Approve-with-changes.** The per-locus rate multiplier is correctly formulated (`observed / expected` slips, shrunk toward 1), local, and byte-identity-preserving; it collapses to the group rate on thin/clean loci and the existing cohort `level_refit` test is unaffected. No Blocker/Major. Minors: the loop bound is still named `theta_max_rounds` though it now governs both refits, and the multiplier↔burn-in-level interaction deserves a documented note.

## 3. Execution status

- `cargo fmt --check` clean · `cargo clippy --lib -D warnings` clean · `cargo test --lib` = **1278 passed, 0 failed, 2 ignored** (+2 I2 tests; byte-identity + level_refit unchanged).
- Sub-agent fan-out not used (overload); reviewed inline. "Needs verification": 0.

## 4. Open questions and assumptions

1. **The multiplier is estimated from parent lengths, applied to candidate lengths.** `expected_slipped` sums the group level at each read's *parent* allele length; the multiplier then scales the group level at each *candidate* length in `compute_data_ll`. This is the intended "overall rate nudge, group keeps the length-dependence," and is correct, but it assumes the locus's rate deviation is length-independent (a scalar). Fine for v1; a per-length per-locus rate is out of scope.

## 5. Top 3 priorities
1. **Mi1** — rename `theta_max_rounds` → `refit_max_rounds` (it now bounds the θ *and* level refit).
2. **Mi2** — document that the per-locus multiplier and the burn-in group level form a hierarchy (group = prior mean estimated from genotypes; `m` = per-locus deviation), so they don't double-count.
3. Nits below.

## 6. Findings

### Minor

**Mi1: [em.rs](../../../../src/ssr/cohort/em.rs) — `theta_max_rounds` now bounds both refits.**
The loop refits `θ_locus` *and* the level multiplier and converges on both, but the bound is still named for θ only. Confidence: High (clarity).
*Fix:* rename `EmCfg::theta_max_rounds` → `refit_max_rounds` (the per-`θ`/`level` `*_shrink`/`*_tol` stay name-accurate). Touches the field, `dev_default`, the loop, and the `theta_max_rounds_zero_…` test.

**Mi2: [em.rs](../../../../src/ssr/cohort/em.rs) `attribute_locus` / `refit_level_multiplier` — document the multiplier ↔ group-level hierarchy.**
The burn-in's `reduce_level` estimates the group level from the called genotypes (not from `m`), while genotyping scales it per locus by `m`; at convergence `m ≈ 1` on loci that match the group, so the two do not double-count. This is correct but non-obvious. Confidence: Medium.
*Fix:* a short doc note on `refit_level_multiplier` stating that the group level is the prior mean (estimated upstream from genotypes) and `m` is the per-locus deviation around it.

### Nits
- `expected_slipped` is an `f64` sum; it is summed in fixed per-locus order and a locus is never split across threads, so byte-identity holds — already noted in the `attribute_locus` doc. No action.
- `LEVEL_MULT_MAX = 10.0` is a documented provisional constant (pinned in F2), consistent with the module's other dev defaults. No action.
- `attribute_locus` discards `_eps` from `sample_chemistry`; fine (it only needs the level), but a small `sample_level(...)` helper would avoid computing/ignoring `eps`. Optional.

## 7. Out of scope observations
None.

## 8. Missing tests to add now
None required — `refit_level_multiplier_collapses_to_one_…` covers the matched / no-data / stuttery cases, `attribute_locus_*` cover the sufficient statistics, and `level_multiplier_recovers_calls_under_an_underestimated_group_level` covers the end-to-end adaptation. (A strict "multiplier > 1 on this locus" assertion via an internal probe would be redundant with the unit test.)

## 9. What's good
- `LocusSlipFit` bundles the I1 and I2 sufficient statistics from a single attribution pass — no second read sweep, and the shape/level refits share one consistent hard-labelling.
- The multiplier form `(slipped + strength)/(expected + strength)` gives the anti-oscillation guarantee *and* a safe denominator in one expression: `expected ≈ 0` cannot divide-by-zero, and a thin locus collapses to `1`.
- Applying `m` as a multiplier on `group_level(candidate)` (not a replacement) preserves the group's length-dependence while letting the locus nudge the overall rate — exactly the design intent, and the `[0,1]` clamp on the resulting per-read level bounds any runaway.

## 10. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --lib --all-features -- -D warnings` · `cargo test --lib ssr::cohort::em`

### Author response convention
Address Mi1, Mi2, Nits by id; answer open question 1 first.
