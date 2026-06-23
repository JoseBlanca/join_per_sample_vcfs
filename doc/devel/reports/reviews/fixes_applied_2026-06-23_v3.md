# Fix Application Report: ssr-call genotyping+pre-pass — Milestone C review

**Date:** 2026-06-23
**Source review:** `doc/devel/reports/reviews/ssr_call_genotyping_milestone_c_2026-06-23.md`
**Source state reviewed against:** commits `19b1c61` (C1), `35555a5` (C2–C4)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 2 · Nits: 1 · Missing tests: 1

### Outcome totals
- Applied: 3 (Mi1, Mi2, MT-1) · Deferred: 1 (the Nit — perf, F1)

### Validation summary
- `cargo fmt --check` → pass
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-features` → 0, **1232 lib pass** (+1), integration + doctests green
- Performance check → Skipped — no `benches/` harness covers Stage-2 cohort code.

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | `run_locus_em` non-diploid panic | Apply | Applied (doc) | `em.rs` | Pass |
| Mi2 | Minor | first-of-length sequence selection | Apply | Applied (doc) | `em_init.rs`, `candidate_set.rs` | Pass |
| Nit | Nit | `read_likelihood` runs negligible-`S_θ` slips | Defer | Deferred | — | N/A |
| MT-1 | — | moderate-stutter EM test | Apply | Applied | `em.rs` | Pass |

## 3. Questions asked and answers
None.

## 4. Per-finding log

### Mi1 — non-diploid panic
- **Final status:** Applied (documentation). Kept the loud `assert_eq!` panic — a
  non-diploid call in a diploid-v1 build is a programming error, and loud-fail beats a
  silent no-call (the no-silent-fallback invariant, plan §5). Documented the contract on
  `run_locus_em`. No behavior change.
- **Files:** [em.rs](../../../../src/ssr/cohort/em.rs)

### Mi2 — first-of-length sequence selection
- **Final status:** Applied (documentation). Documented on both `candidate_of_length`
  (π⁰ tally) and `cohort_representative` (C1 nomination) that v1 tracks a single
  representative sequence per length — exact for pure tracts, an approximation for impure
  rungs (a follow-up tied to S2 impure alleles). No behavior change.
- **Files:** [em_init.rs](../../../../src/ssr/cohort/em_init.rs),
  [candidate_set.rs](../../../../src/ssr/cohort/candidate_set.rs)

### Nit — negligible-`S_θ` slips
- **Final status:** Deferred — a small-`S_θ` skip is a perf optimization for F1, not a
  correctness issue; recorded in the review.

### MT-1 — moderate-stutter EM test
- **Final status:** Applied — `em_calls_correct_genotypes_under_moderate_stutter` re-runs
  the checkpoint genotypes at stutter level `baseline 0.15` (vs 0.05) and asserts the homs
  and separated hets still call correctly. Guards against an EM that only works near zero
  stutter.
- **Files:** [em.rs](../../../../src/ssr/cohort/em.rs)

## 5–8. Deferred / Disputed / Failed / Blocked
- Deferred: the perf Nit (F1). Others: none.

## 9. Performance check
Skipped — no `Apply` touched perf-sensitive (`benches/`-reachable) code.

## 10–11. Commands
- `cargo fmt` · `cargo clippy --all-targets --all-features -- -D warnings` →
  0, clean · `cargo test --all-features` → 0, 1232 lib pass.

## 12. Notes
- Milestone C is now shipped; checkpoint 1 signed off. Next: Milestone D (the parameter
  pre-pass → checkpoint 2).
