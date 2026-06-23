# Fix Application Report: ssr_call_g0_fit_2026-06-23.md

**Date:** 2026-06-23
**Source review:** `doc/devel/reports/reviews/ssr_call_g0_fit_2026-06-23.md`
**Source state reviewed against:** commit `640a1e6`, branch `ssr-cohort`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 0
- Minors: 3 (Mi1, Mi2, Mi3)
- Nits: 4 (naming, redundant `.max(1)`, implicit NaN-absorption, optional helper extraction)

### Outcome totals
- Applied: 2 (Mi2, Mi3)
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 1 (Mi1 → Step H1)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0
- Nits: 2 applied (redundant `.max(1)` dropped via Mi2; implicit NaN-absorption made explicit), 2 won't-fix (naming, helper extraction)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean (after switching the divisibility check to `is_multiple_of`)
- `cargo test --lib` → 0, **1260 passed, 0 failed, 2 ignored** (+2 from the new G1 tests)
- `cargo doc --no-deps` → not run (doc-affecting change is comment-only; no intra-doc links added)
- `cargo audit` → not run (no dependency change)
- Performance check → not applicable (pre-pass accumulation is not covered by a `benches/` harness; the change is a debug-only assert + a divisor guard + tests)

### Unresolved high-priority findings
- None. (Mi1 deferred to H1 by design — it is the consumer's responsibility; recorded in the plan.)

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| Mi1 | Minor | Two unsynchronised G₀ fallbacks; periods w/o variable loci bypass `fallback_p` | Defer | Deferred | No | None (plan note added) | N/A | Step H1 (`build_param_set` fills every present period) |
| Mi2 | Minor | Integer multiplicity `ploidy / peaks.len()` assumes divisibility | Apply | Applied | No | `prepass.rs` | Pass | No |
| Mi3 | Minor | Missing het-multiplicity + asymmetric-spread tests | Apply | Applied | No | `prepass.rs` | Pass | No |
| Nit-a | Nit | `AlleleSpreadAccum` naming vs `*Profile`/`*Stats` | — | Won't fix | No | None | N/A | No |
| Nit-b | Nit | `peaks.len().max(1)` redundant | Apply | Applied | No | `prepass.rs` | Pass | No |
| Nit-c | Nit | Implicit NaN-absorption in `fit_g0_decay` | Apply | Applied | No | `prepass.rs` | Pass | No |
| Nit-d | Nit | Extract post-loop G₀ block into a helper | — | Won't fix | No | None | N/A | No |

## 3. Questions asked and answers

None — open question 1 (who owns "a `p` for every period present") is resolved by deferring to H1's `build_param_set`, the natural owner; no user decision was required.

## 4. Per-finding log

### Mi1 — two unsynchronised G₀ fallbacks
- **Severity:** Minor
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** No crash risk (`em.rs::period_decay` already falls back to `p=0.5`). The clean fix belongs in `build_param_set` (Step H1, not yet built), which knows exactly which periods the genotyping pass will need and can fill each with `G0FitCfg.fallback_p`. Implementing a period-tracking mechanism in G1 instead would duplicate that responsibility and pick a policy the consumer should own.
- **Implementation summary:** None in code. Added an explicit carry-over note to the H1 step of the plan ([ssr_call_driver.md plan, Step H1](../implementation_plans/ssr_call_driver.md)).
- **Files changed:** `doc/devel/implementation_plans/ssr_call_driver.md` (plan note only).
- **Follow-up:** H1 — fill `pseudocount_decay_per_loci_group` from `G0FitCfg.fallback_p` for every period present.
- **Residual risk:** Until H1, a period with only monomorphic loci silently uses the EM's hardcoded `0.5` rather than `G0FitCfg.fallback_p`. Both are `0.5` today, so behaviour is unchanged.

### Mi2 — integer multiplicity assumes divisibility
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** One-line, clearly-correct guard that makes the diploid assumption visible before polyploidy lands.
- **Implementation summary:** Added `debug_assert!((ploidy as usize).is_multiple_of(peaks.len()), …)` with an explanatory comment in `accumulate_locus`, and dropped the now-redundant `.max(1)` (a `Confident(Peaks)` always has ≥1 peak — this also resolves Nit-b).
- **Review suggestion used verbatim?:** No (the review suggested `% == 0`; clippy required `is_multiple_of`).
- **Adaptation:** Used `is_multiple_of` per `clippy::manual_is_multiple_of`.
- **Verification performed:** clippy clean; all prepass + lib tests pass.
- **Files changed:** `src/ssr/cohort/prepass.rs`.
- **Tests added or modified:** None (debug-only assert; the existing diploid tests exercise the path).
- **Validation:** `cargo clippy --lib … -D warnings` → 0, clean; `cargo test --lib` → 0, 1260 passed.

### Mi3 — missing tests
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The het (2-peak) multiplicity branch and the modal-mode distance reference were uncovered; both are core to the fit.
- **Implementation summary:** Added `g0_spread_counts_het_alleles_as_one_copy_each` (8 separated 6/10 hets → `total_copies == 16`, `distance_weighted == 32`) and `g0_spread_measures_distance_from_an_asymmetric_mode` (6 hom@10 + 3 hom@13 → `total_copies == 18`, `distance_weighted == 18`). Hand-computed values matched on first run.
- **Review suggestion used verbatim?:** No — adapted to the module's `sim`/`collect_loci` helpers and `het_spec`.
- **Files changed:** `src/ssr/cohort/prepass.rs`.
- **Tests added or modified:** the two named tests.
- **Validation:** `cargo test --lib ssr::cohort::prepass` → 0, 14 passed.

### Nit-b — redundant `.max(1)`
- **Final status:** Applied (folded into Mi2 — the `.max(1)` was removed when the divisibility assert made the ≥1-peak contract explicit).

### Nit-c — implicit NaN-absorption in `fit_g0_decay`
- **Final status:** Applied
- **Implementation summary:** Changed the thin-period guard to `accum.total_copies < cfg.min_copies.max(1)`, guaranteeing `total_copies ≥ 1` before the `K̄` division so a `0/0` NaN is impossible even if a caller sets `min_copies = 0`; added a comment. No behaviour change for the real call path (entries always carry ≥2 copies; default `min_copies = 30`).
- **Files changed:** `src/ssr/cohort/prepass.rs`.
- **Validation:** clippy clean; tests pass.

### Nit-a — `AlleleSpreadAccum` naming
- **Final status:** Won't fix. The review marked it defensible; it is a pair of sums (not a profile/histogram), and the doc comment justifies the shape. Renaming is churn without a correctness or clarity gain.

### Nit-d — extract post-loop G₀ block into a helper
- **Final status:** Won't fix. Optional readability refactor; the block is comment-delimited and `accumulate_locus` remains cohesive. Minimal-diff discipline.

## 5. Deferred findings to carry forward
- Mi1 — fill every present period's `G₀` decay from `G0FitCfg.fallback_p` in `build_param_set` (Step H1).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** No — the pre-pass accumulation path is not reachable from any `benches/` harness; the changes are a debug-only assert, a divisor guard, and tests.
- **Outcome:** Skipped — no Apply touched perf-sensitive (benched) code.

## 10. Commands run
- `cargo fmt --check`
- `cargo clippy --lib --all-features -- -D warnings`
- `cargo test --lib`
- `cargo test --lib ssr::cohort::prepass`

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 1260 passed / 0 failed / 2 ignored
- `cargo test --lib ssr::cohort::prepass` → 0, 14 passed

## 12. Notes
- Sub-agent fan-out for the review was unavailable (transient server overload); the review was synthesised inline. The fixes here are unaffected.
- Mi1's deferral is a genuine ownership decision (the consumer `build_param_set` is the right place), not a punt — recorded in the H1 plan step so it cannot be lost.
