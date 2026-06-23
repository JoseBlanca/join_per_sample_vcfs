# Fix Application Report: ssr_call_build_param_set_2026-06-23.md

**Date:** 2026-06-23
**Source review:** `doc/devel/reports/reviews/ssr_call_build_param_set_2026-06-23.md`
**Source state reviewed against:** commit `4168331`, branch `ssr-cohort`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 2 (Mi1, Mi2) · Nits: 3

### Outcome totals
- Applied: 1 (Mi2) · Applied with adaptation: 0 · Already fixed: 0
- Deferred: 1 (Mi1 → H4 decision) · Disputed: 0 · Failed validation: 0
- Nits: 2 applied (n_samples==0 test, shape_by_cell assertion), 1 no-action (out-of-range key — correct by construction)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, **1265 passed, 0 failed, 2 ignored** (+2 from new tests)
- `cargo doc --no-deps` / `cargo audit` → not run (no doc-link or dependency change)
- Performance check → not applicable (pure assembly, not benched)

### Unresolved high-priority findings
- None. (Mi1 is a deliberate H4 design decision, recorded in the plan.)

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|
| Mi1 | Minor | `G₀` backfill universe = `shape_by_period`, not genotyped-loci periods | Defer | Deferred | None (plan note) | N/A | H4 decision |
| Mi2 | Minor | Missing `g0`-only-period test | Apply | Applied | `driver.rs` | Pass | No |
| Nit-a | Nit | `n_samples == 0` degenerate untested | Apply | Applied | `driver.rs` | Pass | No |
| Nit-b | Nit | Happy path doesn't assert `stutter_shape_by_cell` | Apply | Applied | `driver.rs` | Pass | No |
| Nit-c | Nit | Out-of-range `group_of_sample` key ignored | — | No action | None | N/A | No (correct by construction) |

## 3. Questions asked and answers

None — open question 1 (backfill universe) is a design decision deferred to H4, recorded in the plan; no user input was needed to land the H1 fixes.

## 4. Per-finding log

### Mi1 — `G₀` backfill universe
- **Final status:** Deferred. Benign today (the EM default `0.5` equals `G0FitCfg::dev_default().fallback_p`) and documented in the function. The clean close — passing the genotyped-loci period set into `build_param_set` — belongs in H4, which owns the loci stream. Added a carry-over note to the H4 plan step.
- **Files changed:** `doc/devel/implementation_plans/ssr_call_driver.md` (plan note).
- **Residual risk:** a period with confident genotypes but zero slips uses the EM's `0.5` rather than `G0FitCfg.fallback_p`; identical values today.

### Mi2 — `g0`-only-period test
- **Final status:** Applied. Added `build_param_set_carries_a_g0_only_period` (`shape = {2}`, `g0 = {2,4}` → period 4's fitted value survives), exercising the `est.g0_by_period.clone()` base of the backfill.
- **Files changed:** `src/ssr/cohort/driver.rs`. **Tests added:** the named test.
- **Validation:** `cargo test --lib ssr::cohort::driver` → 0, 10 passed.

### Nit-a — degenerate `n_samples == 0`
- **Final status:** Applied. Added `build_param_set_on_an_empty_cohort_is_ok_and_empty` pinning the `Ok`/empty contract.
- **Files changed:** `src/ssr/cohort/driver.rs`.

### Nit-b — assert `stutter_shape_by_cell`
- **Final status:** Applied. Added a `stutter_shape_by_cell[&(SampleGroupId(0), 2)]` assertion to the happy-path test.
- **Files changed:** `src/ssr/cohort/driver.rs`.

### Nit-c — out-of-range key
- **Final status:** No action. `grouped.group_of_sample` only ever holds `0..n_samples` indices (pre-pass contract), so reading exactly that range is correct; a guard would be dead code.

## 5. Deferred findings to carry forward
- Mi1 — decide at H4 whether to backfill `G₀` over the genotyped-loci period set instead of `shape_by_period`.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no Apply touched perf-sensitive (benched) code (pure assembly + tests).

## 10. Commands run
- `cargo fmt --check`
- `cargo clippy --lib --all-features -- -D warnings`
- `cargo test --lib`

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 1265 passed / 0 failed / 2 ignored

## 12. Notes
- Review and fixes both done inline (sub-agent fan-out overloaded earlier); the small pure-function scope made inline review reliable.
