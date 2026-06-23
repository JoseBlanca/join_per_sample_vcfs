# Code Review: ssr-call genotyping+pre-pass — Milestone E
**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator, focused inline pass)
**Scope:** Milestone E (commits `ab12603` E1, `7df023e` E2)
**Status:** Approve-with-changes

---

## 1. Scope
- **In-scope:** [inbreeding.rs](../../../../src/ssr/cohort/inbreeding.rs) (E1), the
  `em.rs` surgery (`run_locus_em_with` / `posterior_hom`),
  [vcf_out.rs](../../../../src/ssr/cohort/vcf_out.rs) FP-control + output (E2).
- **Categories:** reliability, errors, naming, idiomatic, refactor_safety, smells, extras
  (numerical correctness / determinism). Skipped unsafe_concurrency, tooling.

## 2. Verdict
**Approve-with-changes.** The `F` loop separates inbred from outbred samples and is
deterministic (`FixedPointAccum`); the allele-balance defence no-calls a depth-inflated
false het — the SNP-path blind spot — exactly as intended. Findings are a convergence
gap, a small duplication, a masked size-mismatch, and a valuable end-to-end test.

## 3. Execution status
- `cargo fmt --check` → pass · `cargo clippy --all-targets --all-features -- -D warnings`
  → pass · `cargo test --all-features` → **1250 lib pass**.
- Needs-verification findings: 0.

## 4. Top 3 priorities
1. **Mi1** — `run_cohort_em` convergence checks only `Δf`, not `Δlevel`; the loop can stop
   before the level settles (the plan's `converged(Δf, Δlevel)`).
2. **MT-1** — add an end-to-end integration test (prepass → cluster → `run_cohort_em` →
   FP control → `format_vcf_record`) — nothing yet exercises the whole caller composed.
3. **Mi2** — `inbreeding::add_bin` duplicates `prepass`'s; expose and reuse one.

## 5. Findings

### Minor

- `src/ssr/cohort/inbreeding.rs` (`run_cohort_em`) — **[Minor]** convergence ignores `Δlevel`
- **Confidence:** High
- **Problem:** The loop breaks on `max|Δf| < f_tol` only. The per-group level is refit
  each round but never gates termination, so the loop can stop with the level still
  moving (or, conversely, iterate `max_rounds` when only the level is unsettled).
- **Why it matters:** The plan specifies `converged(Δf, Δlevel)`; a level still in motion
  at exit slightly mis-weights the final calls.
- **Suggested fix:** track the previous `level_per_group` and also require the largest
  `baseline`/`slope` change to be below a `level_tol`.

- `src/ssr/cohort/inbreeding.rs` (`add_bin`) — **[Minor]** duplicated helper
- **Confidence:** High
- **Problem:** `add_bin` is a verbatim copy of the `prepass` one.
- **Suggested fix:** make `prepass::add_bin` `pub(crate)` and import it (one definition).

### Nits
- `reduce_level` clamps the group index with `.min(n_groups - 1)`, silently masking a
  `group_of_sample` / `level_seed` size mismatch. A `debug_assert` that the group is in
  range would surface a real wiring bug instead of hiding it.
- `site_qual` is a documented posterior proxy; keep the doc pointer to the exact-AF kernel
  (verify-fix #7b) prominent so it isn't mistaken for the final QUAL.

## 6. Out of scope observations
- `benches/psp_writer_perf.rs:386` — pre-existing bench panic, unchanged.

## 7. Missing tests to add now
- `full_pipeline_calls_and_emits_a_variant_vcf_line` — **input:** a clean simulated
  cohort. **Flow:** `run_prepass_stats` → `estimate` → `group_samples` → build a `ParamSet`
  from the grouped output → `run_cohort_em` → `apply_fp_control` per locus →
  `format_vcf_record`. **Bug it catches:** an interface mismatch between the pieces (the
  modules each test in isolation; nothing composes them). **Body:** assert the emitted
  line is `PASS`, carries the expected REF/ALT, and the high-depth genotypes are correct.

## 8. What's good
- The `F` estimator is the textbook excess-homozygosity form, reduced through
  `FixedPointAccum` for determinism, and the test demonstrates real separation
  ([inbreeding.rs](../../../../src/ssr/cohort/inbreeding.rs)).
- The EM surgery is genuinely minimal-ripple: `run_locus_em` stays a wrapper, so every
  C4 call site is untouched ([em.rs](../../../../src/ssr/cohort/em.rs)).
- The allele-balance test encodes the depth-inflation blind spot as an executable
  regression (a 95/5 "het" → no-call) ([vcf_out.rs](../../../../src/ssr/cohort/vcf_out.rs)).

## 9. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --all-targets --all-features -- -D warnings` ·
  `cargo test --all-features`

### Author response convention
Address each finding by ID (Mi1, Mi2, Nits, MT-1).
