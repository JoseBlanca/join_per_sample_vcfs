# Fix Application Report: ssr-call genotyping+pre-pass — Milestone A review

**Date:** 2026-06-23
**Source review:** `doc/devel/reports/reviews/ssr_call_genotyping_milestone_a_2026-06-23.md`
**Source state reviewed against:** commits `2229b16` (A1), `4e74b78` (A2)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 2 · Nits: 2

### Outcome totals
- Applied: 4 (Mi1, Mi2, both Nits) · Deferred: 0 · Disputed: 0 · Already fixed: 0

### Validation summary
- `cargo fmt --check` → pass
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-features` → 0, **1184 lib pass** (+3), integration + doctests green
- `cargo doc --no-deps` → not run (doc-comment-only change, no new intra-doc links beyond existing)
- `cargo audit` → not run (no dependency change)
- Performance check → Skipped — no `Apply` touched a `benches/`-covered hot path (Stage-2
  cohort code has no bench harness yet; the simulator is `#[cfg(test)]`).

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | `FixedPointAccum::add` swallows non-finite/extreme inputs | Apply | Applied | `param_estimation.rs` | Pass |
| Mi2 | Minor | Missing sim tests (separated-het, per-group shape) | Apply | Applied | `sim.rs` | Pass |
| Nit-1 | Nit | `below()` reads as a predicate | Apply | Applied | `sim.rs` | Pass |
| Nit-2 | Nit | `build_tract().len()` allocates for a length | Apply | Applied | `sim.rs` | Pass |

## 3. Questions asked and answers
None.

## 4. Per-finding log

### Mi1 — `FixedPointAccum::add` silently swallows non-finite/extreme inputs
- **Final status:** Applied
- **Implementation:** added `debug_assert!(x.is_finite(), …)` to `add`; documented the
  magnitude contract on the type and the method. Test-first: added
  `fixed_point_accum_rejects_non_finite_in_debug` (`#[should_panic(expected = "non-finite")]`)
  — confirmed it panics in debug (test passes).
- **Files:** [src/ssr/cohort/param_estimation.rs](../../../../src/ssr/cohort/param_estimation.rs)
- **Validation:** `cargo test --lib param_estimation` → 5 pass (incl. the new guard test).

### Mi2 — Missing simulator tests
- **Final status:** Applied
- **Implementation:** added `separated_het_deposits_support_at_both_allele_lengths`
  (both 4- and 9-unit tracts supported >50 reads) and
  `differing_group_shape_changes_the_slip_magnitude_distribution` (steeper `decay` ⇒ larger
  mean `|Δ|`), per the review's `## Missing tests`.
- **Files:** [src/ssr/cohort/sim.rs](../../../../src/ssr/cohort/sim.rs)
- **Validation:** `cargo test --lib ssr::cohort::sim` → 8 pass.

### Nit-1 — `below()` naming
- **Final status:** Applied — renamed `SplitMix64::below` → `index_below` (3 call sites,
  mechanical `sed`).
- **Files:** [src/ssr/cohort/sim.rs](../../../../src/ssr/cohort/sim.rs)

### Nit-2 — alloc-for-length in the record builder
- **Final status:** Applied — replaced
  `build_tract(&locus.motif, locus.ref_units).len()` with
  `locus.motif.period() * locus.ref_units as usize`.
- **Files:** [src/ssr/cohort/sim.rs](../../../../src/ssr/cohort/sim.rs)

## 5. Deferred findings to carry forward
None.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no `Apply` touched perf-sensitive code (`#[cfg(test)]` simulator + a debug-only
assert; no `benches/` harness covers Stage-2 cohort code).

## 10. Commands run
- `cargo fmt`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-features`

## 11. Command results
- `cargo clippy …` → 0, clean
- `cargo test --all-features` → 0, 1184 lib pass

## 12. Notes
- The pre-existing `benches/psp_writer_perf.rs:386` panic (out of scope, surfaced in the
  review) remains for a separate fix — not introduced by Milestone A.
