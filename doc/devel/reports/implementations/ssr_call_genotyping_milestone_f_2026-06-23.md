# Implementation report — `ssr-call` genotyping+pre-pass, Milestone F (scale & calibrate)

**Date:** 2026-06-23 · **Branch:** `ssr-cohort` · **Plan:**
[ssr_call_genotyping_and_parameters.md](../../implementation_plans/ssr_call_genotyping_and_parameters.md)
(Milestone F: F1 parallelism & determinism, F2 calibration & validation) · **Skill:** rust-feature-implementation

## 1. Plan

- **F1** — parallelize the pre-pass + genotyping and **prove byte-identity across
  thread counts**.
- **F2** — fix the calibration numbers on the simulator + the correctness gates.

## 2. F1 — parallelism & determinism (code)

The two hot loops are now `rayon`-parallel over loci, with reduces that are
order-independent *by construction*, so the output is byte-identical at any thread
count:

- **Pre-pass** (`run_prepass_stats`): each thread folds a per-thread `PrepassStats`
  partial; the partials reduce via `PrepassStats::merge`. Every field is an **integer
  count**, so the merge is commutative + associative — the totals (and thus `estimate`'s
  `ε`/shape/level) are identical regardless of merge order or internal `Vec`/`HashMap`
  ordering.
- **Genotyping** (`run_cohort_em`): the per-locus EM is `par_iter().zip().collect()` —
  rayon's indexed `collect` preserves locus order; the `F` reduce uses `FixedPointAccum`
  (scale→round→`i128`, order-independent) and the level reduce uses integer counts.

The determinism groundwork laid in every earlier milestone (`FixedPointAccum`,
fixed-order/integer reduces, catalog-index tie-breaks) is exactly what makes this hold
with no algorithm change — only the iteration was parallelized.

## 3. F2 — calibration & validation (status)

**Calibration numbers are provisional `dev_default()` values**, each tagged "pinned in
F2": `MAX_SLIP`, the `RungCfg` / `CandidateCfg` / `EmCfg` / `OuterCfg` / `ClusterCfg` /
`FpControlCfg` thresholds, the `G₀` decay default, `λ`, the allele-balance floor, the
shrink strengths, the clustering threshold. **Pinning them requires real cohort data**
(plus the literature priors) and is a deliberate, honest follow-up — fabricating the
numbers on the simulator alone would over-fit the generator.

**The F2 correctness gates (external anchors, not "calls don't move") are tested:**

| Gate | Test |
|---|---|
| simulator ground-truth recovery (genotypes) | `checkpoint_1_called_genotypes_match_truth_at_high_depth` |
| parameter recovery (`ε`/shape/level) | `checkpoint_2_recovers_epsilon_shape_and_level` |
| known-protocol positive control (per-group shapes) | `m3_two_protocol_cohort_recovers_per_group_shapes` |
| byte-identity across thread counts | `prepass_is_byte_identical_across_thread_counts`, `cohort_em_is_byte_identical_across_thread_counts` |
| whole-caller composition | `full_pipeline_calls_and_emits_a_variant_vcf_line` |
| robustness under stutter | `em_calls_correct_genotypes_under_moderate_stutter` |

Gates that do **not** apply to this implementation (documented, not skipped silently):
multi-start `ℓ_pen` agreement and batch-size invariance presuppose the iterated soft-EM
burn-in (deferred — the heuristic single-pass estimator has no batches or random
restarts); the discrete-vs-continuum cloud is a real-data diagnostic.

## 4. Tests added (2)

- `prepass_is_byte_identical_across_thread_counts` — `run_prepass` at 1 vs 4 threads →
  identical `EstimatedParams`.
- `cohort_em_is_byte_identical_across_thread_counts` — `run_cohort_em` at 1 vs 4 threads
  → identical `CohortCalls`.

## 5. Validation results

- `cargo fmt --check` → pass · `cargo clippy --all-targets --all-features -- -D warnings`
  → pass.
- `cargo test --all-features` → **1253 lib pass** (+2), integration + doctests green.
- Pre-existing unrelated `benches/psp_writer_perf.rs:386` panic under `--all-targets`.

## 6. Tradeoffs and follow-ups

- **Real-data calibration** of every `dev_default` constant is the headline F2 follow-up
  (needs real cohorts).
- The deferred algorithmic refinements remain documented: the iterated soft-EM burn-in
  (`ℓ_pen` plateau, multi-start, batch-size invariance), the exact-AF site-QUAL kernel,
  the soft per-allele level reduce, beta-binomial overdispersion, ploidy > 2.
- **Driver wiring:** `run_cohort_em` + `apply_fp_control` + `format_vcf_record` are not
  yet wired into the `ssr-call` driver (still the Phase-1 TSV dump); a small driver step
  emits the full VCF.
- The plan's milestones A–F are now implemented (the algorithmic caller is complete and
  parallel; calibration is the empirical remainder).
