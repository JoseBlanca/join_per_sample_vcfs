# Code Review + Fixes: SSR marginalized-prior wiring (Phase 3.4)

**Date:** 2026-07-07
**Reviewer:** rust-code-review skill (orchestrator, 3 category sub-agents)
**Scope:** `src/ssr/cohort/em.rs` — wire the marginalize+LOO DM prior into the SSR
per-locus EM behind an `EmCfg.marginalized_prior` toggle (default = plug-in). New:
`MarginalizedFit`, `run_pi_em_marginalized`, `final_calls_marginalized`,
`genotype_pass` dispatcher; the two `run_locus_em_with` call sites route through
`genotype_pass`.
**Status:** Approve-with-changes → fixes applied.

## 1. Execution status

- `cargo test -p pop_var_caller ssr::cohort` → **240 passed** (13 marginalized
  tests), 0 failed. clippy clean on the new code.
- **Correctness verdict (reliability agent, cross-checked against the plug-in
  `run_pi_em`/`final_calls` and the SNP `e_step_cohort_loo` + `run_em_loop`):** the
  algorithm is correct — flat-iteration-0/LOO-after, the convergence guard, state
  threading, and the driver delta all faithfully mirror the SNP engine.
- **Plug-in path byte-identical (refactor_safety agent):** `genotype_pass`'s
  else-branch passes identical arguments to the unchanged `run_pi_em`/`final_calls`;
  `EmCfg` has no `Default` impl and a single constructor, so the new field cannot
  silently default anywhere. All pre-existing tests run through this branch.

## 2. Categories dispatched

`reliability`, `refactor_safety`, `idiomatic`. Audit trail:
`tmp/review_2026-07-07_ssr-marginalized-wiring/`.

## 3. Findings and disposition

| ID | Sev | Cat | Finding | Disposition |
|---|---|---|---|---|
| M1 | Major | reliability | Multi-iteration LOO untested (only clean high-depth, where LOO is inert) | **Applied** — `..._leave_one_out_moves_an_ambiguous_sample` (hand-built cohort; flat-only vs full must differ) |
| M2 | Major | reliability | `F>0` IBD branch never exercised through the loop | **Applied** — `marginalized_inbreeding_is_threaded_and_raises_posterior_homozygosity` (F=0 vs 0.9, directional) |
| M3 | Major (low conf) | reliability | Empty cohort (`n_samples==0`) → NaN `π`, diverges from plug-in | **Applied** — `debug_assert!(n_samples > 0)` documenting the PASS-locus invariant (n≥1 is NaN-free) |
| Mi1 | Minor | reliability | No isolated tests of the two new functions | **Applied** — `..._single_sample_stays_finite_and_collapses_to_g0` |
| Mi2 | Minor | reliability | No marginalized-vs-plug-in agreement test | **Applied** — `marginalized_and_plugin_agree_on_clean_high_depth` |
| Mi3 | Minor | idiomatic | `final_calls_marginalized` duplicates `final_calls`'s ~18-line tail | **Deferred to 3.5** — the benchmark likely deletes one path (idiomatic agent's own recommendation); added the PANIC-FREE rationale to the marginalized copy so it isn't lost |
| N1 | Nit | reliability | Comments say "iteration 1" but the loop is 0-indexed | **Applied** — reworded to "first iteration (`iteration == 0`)" |
| N2/N3 | Nit | idiomatic | Per-iteration `Vec<Vec>` alloc; 3-tuple return | **Won't fix now** — per-locus scalar loop; matches the plug-in `run_pi_em` convention; revisit if 3.5 keeps this path and it profiles hot |
| — | obs | reliability | Changed `π` semantics (raw vs `G₀`-regularised) under the toggle — downstream QUAL/AF | **Noted for 3.5** — toggle defaults off; confirm downstream tolerance when the benchmark runs it on |

## 4. Validation

`cargo test -p pop_var_caller ssr::cohort` → 240 passed after fixes. The two
load-bearing new tests: `..._leave_one_out_moves_an_ambiguous_sample` fails if the
LOO threading breaks; `marginalized_inbreeding_...` fails if `F` is dropped — the
exact regressions M1/M2 warned about.

## 5. What's good

- The loop is a faithful port of the SNP cohort-LOO structure, verified line-by-line.
- The toggle keeps the plug-in path byte-identical (no `Default` impl → no silent
  field default), a clean, low-risk staging for the benchmark.
- SSR keeps its own mode-centred `G₀` seed and its slim loop; the SNP engine is
  untouched.
