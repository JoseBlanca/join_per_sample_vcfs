# Fix Application Report: ssr-call genotyping+pre-pass — Milestone D (D1+D2) review

**Date:** 2026-06-23
**Source review:** `doc/devel/reports/reviews/ssr_call_genotyping_milestone_d_d1d2_2026-06-23.md`
**Source state reviewed against:** commit `4c3d756`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 2 · Nits: 2 · Missing tests: 2

### Outcome totals
- Applied: 6 (Mi1, Mi2, both Nits, MT-1, MT-2) · Deferred/Disputed: 0

### Validation summary
- `cargo fmt --check` → pass
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-features` → 0, **1238 lib pass** (+2), integration + doctests green
- Performance check → Skipped (no `benches/` harness covers Stage-2 cohort code).

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | het `ε` upward bias | Apply | Applied (doc) | `prepass.rs` | Pass |
| Mi2 | Minor | unclamped level line | Apply | Applied (doc) | `prepass.rs` | Pass |
| Nit-1 | Nit | inline `Rungs` path | Apply | Applied | `prepass.rs` | Pass |
| Nit-2 | Nit | `|Δ|>MAX_SLIP` cap unnoted | Apply | Applied (doc) | `prepass.rs` | Pass |
| MT-1 | — | het contribution test | Apply | Applied | `prepass.rs` | Pass |
| MT-2 | — | determinism test | Apply | Applied | `prepass.rs` | Pass |

## 3. Questions asked and answers
None.

## 4. Per-finding log

### Mi1 — het `ε` bias
- **Final status:** Applied (documentation) — a `BIAS NOTE` on `accumulate_locus`
  explaining the small upward `ε` bias off het faithful peaks and that the deferred
  soft-EM reduce removes it. No behavior change.

### Mi2 — unclamped level line
- **Final status:** Applied (documentation) — `fit_level` documents that the line is a
  seed clamped at evaluation (preserved raw for the E1 refit). No behavior change.

### Nit-1 / Nit-2
- **Final status:** Applied — imported `Rungs` (dropped the inline fully-qualified path);
  noted the `|Δ| > MAX_SLIP` drop is the deliberate spec cap.

### MT-1 — het contribution test
- **Final status:** Applied — `separated_hets_contribute_two_length_bins`: a 6/10 het
  cohort deposits faithful/slipped stats at *both* allele lengths and still recovers `ε`.

### MT-2 — determinism test
- **Final status:** Applied — `run_prepass_is_deterministic`: two runs over the same loci
  yield equal `EstimatedParams` (guards against an accidental order-dependence).

## 5–8. Deferred / Disputed / Failed / Blocked
None.

## 9. Performance check
Skipped — no `Apply` touched perf-sensitive code.

## 10–11. Commands
- `cargo fmt` · `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean ·
  `cargo test --all-features` → 0, 1238 lib pass.

## 12. Notes
- D1+D2 shipped; checkpoint 2 signed off. Next: **D3** (clustering + per-`(group,
  period)` shape + the `ε`-freeze check — the M3 milestone).
