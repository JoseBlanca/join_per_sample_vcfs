# Fix Application Report: ssr_call_theta_locus_2026-06-24.md

**Date:** 2026-06-24
**Source review:** `doc/devel/reports/reviews/ssr_call_theta_locus_2026-06-24.md`
**Source state reviewed against:** commit `46df902`, branch `ssr-cohort`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 2 (Mi1, Mi2) · Nits: 2

### Outcome totals
- Applied: 2 (Mi1 shared `add_slip`; the `theta_max_rounds=0` pin test) · Deferred: 1 (Mi2 warm-start optimization) · No action: 1 Nit (`compute_data_ll` arg count)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, **1276 passed, 0 failed, 2 ignored**
- `cargo doc --no-deps` / `cargo audit` → not run (no public-doc-link or dependency change)
- Performance check → not applicable (per-locus EM is not covered by a `benches/` harness; the cost note is Mi2)

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | `add_slip` duplicated in em + prepass | Apply | Applied | `param_estimation.rs`, `em.rs`, `prepass.rs` | Pass |
| Mi2 | Minor | Per-locus cost grows with the θ loop | Defer (note) | Deferred | None | N/A |
| Nit-a | Nit | `compute_data_ll` 10 args | — | No action | None | N/A |
| Nit-b | Nit | Pin `theta_max_rounds=0` == pre-I1 | Apply | Applied | `em.rs` | Pass |

## 3. Questions asked and answers
- **OQ1 (hard-label valley attribution):** kept as the documented hard-label approximation (deterministic via `min_by_key`, consistent with `reduce_level`); the soft split is the deferred refinement. No change.

## 4. Per-finding log

### Mi1 — shared `add_slip`
- **Final status:** Applied. Lifted `SlipProfile::add_slip(&mut self, delta: i32, count: u64)` into `param_estimation.rs` (next to `SlipProfile`); removed the two private copies and called the method from `prepass::accumulate_locus` (`*count as u64`) and `em::attribute_locus_slips`. Dropped the now-unused `MAX_SLIP` import from `em.rs`.
- **Files changed:** `src/ssr/cohort/param_estimation.rs`, `src/ssr/cohort/em.rs`, `src/ssr/cohort/prepass.rs`.
- **Validation:** prepass + em + full lib tests pass (1276); byte-identity test still green (the binning is unchanged).

### Mi2 — per-locus cost
- **Final status:** Deferred. The cost (data-ll recomputed + π restarted from `pi0` each θ round, bounded by `theta_max_rounds=3` + early `shapes_close` break) is already documented in the θ-loop comment and the module doc. Warm-starting π across θ rounds is a clean later optimization that would not change converged results; it folds naturally into the Milestone-J throughput work. No code change.

### Nit-a — arg count
- **Final status:** No action. `compute_data_ll`'s 10 args are an internal extraction under `#[allow]`; a context struct is premature.

### Nit-b — pin test
- **Final status:** Applied. Added `theta_max_rounds_zero_reproduces_the_seed_shape_result` (the refit disabled → still calls the clean truth), pinning that the refit is purely additive over the pre-I1 EM.

## 5. Deferred findings to carry forward
- Mi2 — warm-start π across θ rounds (a throughput optimization, candidate for Milestone J).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — the per-locus EM is not reachable from a `benches/` harness. The Mi2 cost is bounded and documented.

## 10. Commands run
- `cargo fmt --check`
- `cargo clippy --lib --all-features -- -D warnings`
- `cargo test --lib`

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 1276 passed / 0 failed / 2 ignored

## 12. Notes
- The shared `add_slip` is now the single home for `MAX_SLIP`-bounded signed-delta binning — used by both the pre-pass `θ`/per-sample accumulators and the per-locus `θ_locus` M-step.
