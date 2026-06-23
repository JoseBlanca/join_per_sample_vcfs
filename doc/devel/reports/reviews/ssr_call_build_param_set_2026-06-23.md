# Code Review: ssr-call build_param_set (Step H1)

**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator; inline synthesis — see Execution status)
**Scope:** commit `4168331` on branch `ssr-cohort` — `build_param_set` (pre-pass → frozen `ParamSet`) + `SsrCallError::UnresolvedSamples`
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** PR diff, commit `4168331` (Step H1).
- **Reviewed against:** `4168331`, branch `ssr-cohort`.
- **In-scope files:** [driver.rs](../../../../src/ssr/cohort/driver.rs) — `build_param_set`, the `UnresolvedSamples` error variant, 3 new tests + test helpers.
- **Out of scope:** the rest of the driver (TSV dump path, unchanged); the consumed types (`EstimatedParams`/`GroupedParams`/`ParamSet`, reviewed in their own steps).
- **Categories considered:** reliability, errors, naming, idiomatic, defaults, refactor_safety, extras (diff-vs-intent). `unsafe_concurrency` n/a (pure function). `module_structure`/`tooling` n/a (single-file, no `Cargo.toml`).

## 2. Verdict

**Approve-with-changes.** `build_param_set` is a correct, pure assembly with the decision-E hard error and the Mi1 `G₀` backfill both implemented as designed. No Blocker/Major. One Minor (the backfill universe is the pre-pass-characterized periods, not the genotyped loci's periods — a documented boundary) and two test-completeness Minors/Nits.

## 3. Execution status

- `cargo fmt --check` — clean (exit 0).
- `cargo clippy --lib --all-features -- -D warnings` — clean (exit 0).
- `cargo test --lib` — **1263 passed, 0 failed, 2 ignored** (+3 H1 tests).
- Sub-agent fan-out **not used** — the parallel dispatch overloaded the server on the G1 review, and this diff is small (one pure function + one error variant), so the orchestrator reviewed inline. "Needs verification": 0.

## 4. Open questions and assumptions

1. **Backfill universe.** The Mi1 resolution backfills `G₀` for periods in `est.shape_by_period` (periods the pre-pass characterized), not for every period in the genotyped loci. A period whose loci yielded confident genotypes but **zero slips** (no shape entry) is therefore not backfilled and falls to the EM's own `0.5`. This is documented in the function and is benign (the EM default equals `G0FitCfg::dev_default().fallback_p`), but it means "a `p` for every period present" (arch §9) holds for *characterized* periods only. See M-Mi1. Affects whether H4 should pass the loci-derived period set.

## 5. Top 3 priorities

1. **Mi1** — decide whether the `shape_by_period` backfill universe is sufficient or H4 should supply the genotyped-loci period set (close the zero-slip-period gap).
2. **Mi2** — test that a `g0`-only period (present in `g0_by_period`, absent from `shape_by_period`) survives the assembly.
3. **Nits** — degenerate `n_samples == 0`; assert `stutter_shape_by_cell` in the happy-path test.

## 6. Findings

### Minor

**Mi1: [driver.rs](../../../../src/ssr/cohort/driver.rs#L97-L107) — `G₀` backfill universe is `shape_by_period`, not the genotyped loci's periods.**
A period with confident genotypes but no recorded slips has no `shape_by_period` entry, so it is neither fitted nor backfilled here and falls to the EM's hardcoded `0.5`. Benign today (`fallback_p` defaults to the same `0.5`), and documented in the function doc, but the "single source of truth" goal is only fully met for characterized periods. Confidence: High (behavioural boundary, verified against the code).
*Fix (choose at H4):* either accept the documented boundary as-is, or have H4 pass the set of periods actually present in the (subset/cohort) loci into `build_param_set` and backfill over *that* set. Record the choice in the H4 step. No code change required in H1 if the boundary is accepted.

**Mi2: missing test — a `g0`-only period survives assembly.**
All three tests put `g0_periods ⊆ shape_periods`, so the `est.g0_by_period.clone()` base of the backfill (which carries periods absent from `shape_by_period`) is never exercised. Confidence: High.
*Fix:* see section 8.

### Nits

- [driver.rs](../../../../src/ssr/cohort/driver.rs#L83-L95) — `n_samples == 0` yields an `Ok` `ParamSet` with an empty `group_of_sample`. Unreachable in practice (the merger rejects an empty cohort upstream with `NoInputs`), so not worth a guard, but a one-line test pinning the degenerate case would document the contract.
- [driver.rs](../../../../src/ssr/cohort/driver.rs#L227) — the happy-path test asserts most fields but not `stutter_shape_by_cell`; add one assertion so a future mismapping there is caught.
- A `grouped.group_of_sample` key `≥ n_samples` would be silently ignored (only `0..n_samples` is read). Contractually the pre-pass only emits `0..n_samples` indices, so this is correct by construction; no action.

## 7. Out of scope observations

None.

## 8. Missing tests to add now

Under `driver.rs`:

- `build_param_set_carries_a_g0_only_period` — `est` with `shape_by_period = {2}` and `g0_by_period = {2, 4}` (period 4 fitted but uncharacterized for shape); assert `pseudocount_decay_per_loci_group[&4]` equals the fitted value (proves the clone-first base, not just the shape-driven backfill).
- (Nit) `build_param_set_on_an_empty_cohort_is_ok_and_empty` — `n_samples = 0` → `Ok` with empty `group_of_sample`; documents the degenerate contract.

## 9. What's good

- The hard error returns **before** any `ParamSet` is built and collects *all* offending indices in one pass — the user sees every unresolved sample at once, not one-at-a-time. [driver.rs:83-95](../../../../src/ssr/cohort/driver.rs#L83-L95)
- The clone-then-`or_insert` backfill makes the fitted value authoritative and the fallback secondary in three lines, with the rationale documented inline. [driver.rs:99-107](../../../../src/ssr/cohort/driver.rs#L99-L107)
- Full struct-literal construction of `ParamSet` (no `..Default`) means a future field addition is a compile error here, not a silent default.

## 10. Commands to re-verify

- `cargo fmt --check`
- `cargo clippy --lib --all-features -- -D warnings`
- `cargo test --lib ssr::cohort::driver`

### Author response convention
Address Mi1, Mi2, and the Nits by identifier; answer open question 1 first.
