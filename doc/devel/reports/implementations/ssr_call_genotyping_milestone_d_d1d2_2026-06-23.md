# Implementation report — `ssr-call` genotyping+pre-pass, Milestone D (D1+D2 → checkpoint 2)

**Date:** 2026-06-23 · **Branch:** `ssr-cohort` · **Plan:**
[ssr_call_genotyping_and_parameters.md](../../implementation_plans/ssr_call_genotyping_and_parameters.md)
(Milestone D: D1 confident-genotype fitting, D2 burn-in/measure) · **Skill:** rust-feature-implementation

> This report covers **D1 + D2**, which reach **checkpoint 2** (the pre-pass recovers
> the simulator's known parameters). **D3** (sample-group clustering + per-`(group,
> period)` shape + the `ε`-freeze check — the M3 milestone) comes *after* the
> checkpoint-2 human gate.

## 1. Plan

The parameter pre-pass: estimate `ε`, the stutter shape, and the stutter level from
confident genotypes, so the genotyping half (C) no longer needs them supplied.

- **D1** — accumulate confident-genotype sufficient statistics (per-period
  `SlipProfile`, per-sample `SampleStutterStats`).
- **D2** — measure: `ε`, the per-period shape, the per-sample level line.

## 2. Assumptions / decisions

- **Gate = B1's heuristic resolution** (param-free). Because the gate doesn't depend
  on the params, the estimator is a **single pass** — no burn-in iteration needed to
  recover params on clean data. The model-based 1-vs-2-peak BIC gate (Q-P7), the soft
  full-cohort EM responsibility reduce (C2 amendment), the `ℓ_pen`-plateau burn-in
  with multi-start, and the `FixedPointAccum` float reduce are the refinements that
  matter once the gate co-evolves with the params or the data is messy — they arrive
  with the parallel machinery (F1). Documented in the module header.
- **Hard-label het inner-valley reads** to the nearer allele (a documented
  approximation; the soft responsibility split is the refinement).
- **Estimators:** `ε` = within-tract base-mismatch fraction off faithful reads; shape
  = direction split + geometric-decay MLE `(mean|Δ| − 1)/mean|Δ|` (the inverse of the
  simulator's geometric draw); level = weighted-least-squares line `baseline +
  slope·length` over the per-length faithful/slipped rates.
- **Statistics are integer counts** (`SlipProfile`/`SampleStutterStats`), so the
  single-pass reduce is order-independent by construction; the float `FixedPointAccum`
  is exercised by the deferred soft-EM/`ℓ_pen` layer.

## 3. Changes made

- [prepass.rs](../../../../src/ssr/cohort/prepass.rs) (new) — `PrepassStats`,
  `EstimatedParams`, `accumulate_locus` (D1), `estimate`/`estimate_shape`/`fit_level`
  (D2), `run_prepass`.
- [mod.rs](../../../../src/ssr/cohort/mod.rs) — wired `prepass`.

## 4. Tests added (4)

- **checkpoint 2** — `checkpoint_2_recovers_epsilon_shape_and_level`: a cohort of
  cohort-recurrent homozygotes spread over lengths 7..=12 under a known chemistry;
  recovered `ε` (±0.002), shape (contraction bias, direction split ±0.08, decay
  ±0.06), and mean level baseline (±0.03) all track the injected truth.
- `recovers_a_higher_epsilon` (ε = 0.02), `recovers_a_higher_stutter_level`
  (level = 0.25) — the estimators track changes, not just one operating point.
- `estimate_is_empty_without_confident_genotypes` — a thin cohort yields no stats.

## 5. Validation results

- `cargo fmt --check` → pass · `cargo clippy --all-targets --all-features -- -D warnings`
  → pass.
- `cargo test --all-features` → **1236 lib pass** (+4), integration + doctests green.
- ▶ **Checkpoint 2 holds:** the pre-pass recovers the simulator's `ε` / shape / level.
- Pre-existing unrelated `benches/psp_writer_perf.rs:386` panic under `--all-targets`.

## 6. Tradeoffs and follow-ups

- **D3** (clustering into sample groups, per-`(group, period)` shape, `ε`-freeze BIC
  check — M3) is next, after the checkpoint-2 gate.
- The deferred refinements (model-based gate, soft-EM reduce, `ℓ_pen`-plateau
  burn-in + multi-start, `FixedPointAccum` determinism) land with E/F.
- The pre-pass output is not yet wired into the genotyping pass (E1 consumes it via
  the outer `F`/level loop).
- **Checkpoint 2 is a human gate** — the loop halts here for sign-off before D3.
