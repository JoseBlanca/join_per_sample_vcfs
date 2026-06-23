# Fix Application Report: ssr-call genotyping+pre-pass — Milestone E review

**Date:** 2026-06-23
**Source review:** `doc/devel/reports/reviews/ssr_call_genotyping_milestone_e_2026-06-23.md`
**Source state reviewed against:** commits `ab12603` (E1), `7df023e` (E2)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 2 · Nits: 2 · Missing tests: 1

### Outcome totals
- Applied: 5 (Mi1, Mi2, both Nits, MT-1)

### Validation summary
- `cargo fmt --check` → pass
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-features` → 0, **1251 lib pass** (+1), integration + doctests green
- Performance check → Skipped (no `benches/` harness covers Stage-2 cohort code).

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | convergence ignores `Δlevel` | Apply | Applied | `inbreeding.rs` | Pass |
| Mi2 | Minor | duplicated `add_bin` | Apply | Applied | `prepass.rs`, `inbreeding.rs` | Pass |
| Nit-1 | Nit | masked group-index mismatch | Apply | Applied | `inbreeding.rs` | Pass |
| Nit-2 | Nit | keep `site_qual` deferral prominent | Apply | Already documented | `vcf_out.rs` | N/A |
| MT-1 | — | end-to-end pipeline test | Apply | Applied | `inbreeding.rs` | Pass |

## 3. Questions asked and answers
None.

## 4. Per-finding log

### Mi1 — convergence ignores `Δlevel`
- **Final status:** Applied — added `level_tol` to `OuterCfg`; `run_cohort_em` now also
  tracks the largest level-coefficient change and requires both `Δf < f_tol` and
  `Δlevel < level_tol` before breaking.

### Mi2 — duplicated `add_bin`
- **Final status:** Applied — `prepass::add_bin` made `pub(crate)`; `inbreeding.rs` imports
  and reuses it; the local copy is deleted (one definition).

### Nit-1 — masked group-index mismatch
- **Final status:** Applied — a `debug_assert` that the sample group is in range before the
  `.min(n_groups - 1)` clamp, so a `group_of_sample` / `level_seed` size mismatch surfaces
  in debug instead of hiding.

### Nit-2 — `site_qual` deferral
- **Final status:** Already documented — the proxy nature + the exact-AF kernel pointer are
  on the function and the module header. No change.

### MT-1 — end-to-end pipeline test
- **Final status:** Applied — `full_pipeline_calls_and_emits_a_variant_vcf_line`: a clean
  cohort runs `run_prepass_stats → estimate → group_samples → build ParamSet →
  run_cohort_em → apply_fp_control → site_qual → format_vcf_record`, asserting a `PASS`
  line with an ALT and correct `6/10` het calls. The first test that composes the whole
  caller.

## 5–8. Deferred / Disputed / Failed / Blocked
None.

## 9. Performance check
Skipped — no `Apply` touched perf-sensitive code.

## 10–11. Commands
- `cargo fmt` · `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean ·
  `cargo test --all-features` → 0, 1251 lib pass.

## 12. Notes
- Milestone E shipped. The whole caller now composes end-to-end (pre-pass → genotyping →
  FP control → VCF line), proven by the integration test. Next: Milestone F (parallelism +
  byte-identity, calibration) — the last milestone.
