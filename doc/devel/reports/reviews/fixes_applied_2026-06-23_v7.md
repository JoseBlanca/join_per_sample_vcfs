# Fix Application Report: ssr-call genotyping+pre-pass — Milestone F review

**Date:** 2026-06-23
**Source review:** `doc/devel/reports/reviews/ssr_call_genotyping_milestone_f_2026-06-23.md`
**Source state reviewed against:** commit `0c2cd7b`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 1 · Nits: 0

### Outcome totals
- Applied: 1 (Mi1, documentation)

### Validation summary
- `cargo fmt --check` → pass
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-features` → 0, **1253 lib pass**, integration + doctests green
- Performance check → Skipped (no `benches/` harness covers Stage-2 cohort code).

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | `PrepassStats` struct order not canonical | Apply | Applied (doc) | `prepass.rs` | Pass |

## 3. Questions asked and answers
None.

## 4. Per-finding log

### Mi1 — determinism layer
- **Final status:** Applied (documentation) — added a `DETERMINISM LAYER` note to
  `PrepassStats::merge` clarifying that the *values* (and the derived `EstimatedParams`)
  are byte-identical across thread counts, while the struct's `by_length` `Vec` order is
  not canonical, so determinism tests compare `EstimatedParams` / `CohortCalls`, not the
  raw stats. No behavior change.

## 5–8. Deferred / Disputed / Failed / Blocked
None.

## 9. Performance check
Skipped — no `Apply` touched perf-sensitive code (documentation only).

## 10–11. Commands
- `cargo fmt` · `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean ·
  `cargo test --all-features` → 0, 1253 lib pass.

## 12. Notes
- **Milestone F shipped — all plan milestones A–F are now implemented, reviewed, and
  fixed.** The algorithmic SSR cohort caller is complete, composes end-to-end, and is
  parallel + byte-identical across thread counts. The remaining work is empirical
  (real-data calibration of the provisional constants) and the documented algorithmic
  refinements, plus wiring the genotyping path into the `ssr-call` driver to replace the
  Phase-1 TSV dump.
