# Fix Application Report: ssr-call genotyping+pre-pass — Milestone D D3 review

**Date:** 2026-06-23
**Source review:** `doc/devel/reports/reviews/ssr_call_genotyping_milestone_d_d3_2026-06-23.md`
**Source state reviewed against:** commit `620d091`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 2 · Nits: 1 · Missing tests: 1

### Outcome totals
- Applied: 4 (Mi1, Mi2, Nit, MT-1)

### Validation summary
- `cargo fmt --check` → pass
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-features` → 0, **1242 lib pass** (+1), integration + doctests green
- Performance check → Skipped (no `benches/` harness covers Stage-2 cohort code).

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | hidden reference-depth constant | Apply | Applied | `sample_groups.rs` | Pass |
| Mi2 | Minor | single-linkage chaining | Apply | Applied (doc) | `sample_groups.rs` | Pass |
| Nit | Nit | precision-deflation direction unexplained | Apply | Applied (doc) | `sample_groups.rs` | Pass |
| MT-1 | — | clustering determinism test | Apply | Applied | `sample_groups.rs` | Pass |

## 3. Questions asked and answers
None.

## 4. Per-finding log

### Mi1 — hidden reference-depth constant
- **Final status:** Applied — added `precision_reference_depth` (default 100.0) to
  `ClusterCfg`; `scaled_distance` uses it instead of the inline `100.0`. The `1.0`
  saturation cap is kept with an explaining comment. No behavior change at the default.

### Mi2 — single-linkage chaining
- **Final status:** Applied (documentation) — documented the single-linkage chaining
  limitation on `cluster`, pointing at the deferred BIC split test (F). No behavior
  change (well-separated protocols are correct; the continuum case is the deferral).

### Nit — precision direction
- **Final status:** Applied — a "why" comment on the precision deflation (thinner sample
  ⇒ smaller distance ⇒ merges readily; saturates above the reference depth).

### MT-1 — clustering determinism test
- **Final status:** Applied — `group_samples_is_deterministic`: the M3 cohort grouped
  twice yields equal `GroupedParams` (guards the union-find + group-id labelling order).

## 5–8. Deferred / Disputed / Failed / Blocked
None.

## 9. Performance check
Skipped — no `Apply` touched perf-sensitive code.

## 10–11. Commands
- `cargo fmt` · `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean ·
  `cargo test --all-features` → 0, 1242 lib pass.

## 12. Notes
- Milestone D complete (D1+D2+D3) and shipped. Next: Milestone E (genotyping
  completeness — the `F`/level outer loop, allele-balance, full VCF).
