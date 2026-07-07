# Code Review + Fixes: SSR marginalized DM prior (Phase 3.2)

**Date:** 2026-07-07
**Reviewer:** rust-code-review skill (orchestrator, 3 category sub-agents)
**Scope:** additive new code in `src/ssr/cohort/em.rs` — the marginalize+LOO+`F`
Dirichlet-multinomial genotype prior for SSR (`diploid_dm_inputs`,
`leave_one_out_alpha`, `marginalized_genotype_log_priors`) + unit tests. The SNP
engine and SSR's plug-in `genotype_prior` are untouched; the new prior is unwired
(behind a later `EmCfg` toggle).
**Status:** Approve-with-changes → fixes applied.

*(Consolidated review + fix-application report.)*

## 1. Execution status

- `cargo test -p pop_var_caller ssr::cohort::em` → **36 passed** (9 marginalized-prior
  tests), 0 failed. Build + clippy clean on the new code (only the pre-existing
  `vcf/writer.rs:384` lifetime warning, out of scope).
- **Math verdict (reliability agent, cross-checked against the SNP `e_step` and
  `genetics::dirichlet_multinomial_log_priors`): the Wright-`F` mixture is
  mathematically identical to the SNP engine.** No wrong-result finding — every
  finding was test-coverage or robustness.

## 2. Categories dispatched

`reliability`, `idiomatic`, `naming`. Audit trail:
`tmp/review_2026-07-07_ssr-marginalized-prior/`.

## 3. Findings and disposition

| ID | Sev | Categories | Finding | Disposition |
|---|---|---|---|---|
| M1 | Major | reliability | Mixture tested only at `f=0`/`f=1`; the `−log_sum_alpha` IBD normalization is unpinned (invisible at both boundaries) | **Applied** — added `..._at_f_one_equals_ibd_marginal_log_frequency` (value pin) + `..._mixes_dm_and_ibd_at_intermediate_f` (full-formula reconstruction) |
| Mi1 | Minor | reliability | `diploid_dm_inputs` shape only tested transitively | **Applied** — `diploid_dm_inputs_encodes_homozygote_as_double_count_and_het_as_two_singles` |
| Mi2 | Minor | reliability | Out-of-range `g.i`/`g.j` silently corrupts a neighbouring row | **Applied** — `debug_assert!(g.i < k && g.j < k)` in the builder |
| Mi3 | Minor | reliability | No `k=1` / empty-`genotypes` boundary tests | **Applied** — two boundary tests added |
| Mi4 | Minor | naming | `diploid_dm_shape` reads as "dimensions", not the `(counts, coeffs)` it returns | **Applied** — renamed to `diploid_dm_inputs` |
| Mi5 | Minor | idiomatic, reliability | `softmax` + tests inserted between two `use` lines, stranding `use HashMap` | **Applied** — moved the import up with the others |
| N1 | Nit | reliability | `log_sum_exp(&[a,b])` vs a dedicated `log_sum_exp_2` | **Won't fix** — behaviourally identical, reuses the SSR-local helper (DRY) |
| N2 | Nit | idiomatic | `k` derivable from `alpha.len()` | **Won't fix** — explicit `k` mirrors the shared `dirichlet_multinomial_log_priors(.., n_alleles, ..)` API |
| — | obs | reliability | `α > 0` only `debug_assert` in the DM primitive | **Already handled** — the `G₀` seed is floored at `G0_FLOOR = 1e-12 > 0` and LOO only adds `≥ 0`, so `α'_s > 0` by construction; matches the primitive's contract |

## 4. Validation

`cargo test -p pop_var_caller ssr::cohort::em` → 36 passed after fixes; fmt +
clippy clean on the changed file. The two M1 tests are the load-bearing ones: they
would fail if the `−log_sum_alpha` term or the DM/IBD balance regressed at
`0 < f < 1`, the exact term this port exists to carry over from the SNP engine.

## 5. What's good

- The math port is faithful — `log_p_effective[i] = ln(α_i) − ln(Σα)` and the
  homozygote/het mixture match the SNP `e_step` verbatim, using the already-shared
  `dirichlet_multinomial_log_priors` primitive.
- The seed stays SSR's own (mode-centred `G₀`) — the SNP REF-privileged seed is
  never imported (the §Q1 guardrail).
- SNP untouched — no risk to the SIMD perf gains (the §Q4 decision).
