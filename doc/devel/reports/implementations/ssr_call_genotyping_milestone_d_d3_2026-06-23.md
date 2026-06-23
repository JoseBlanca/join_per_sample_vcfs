# Implementation report — `ssr-call` genotyping+pre-pass, Milestone D D3 (clustering, M3)

**Date:** 2026-06-23 · **Branch:** `ssr-cohort` · **Plan:**
[ssr_call_genotyping_and_parameters.md](../../implementation_plans/ssr_call_genotyping_and_parameters.md)
(Milestone D step D3 — the M3 milestone) · **Skill:** rust-feature-implementation

## 1. Plan

Complete Milestone D: cluster the per-sample `(ε, level)` into data-driven sample
groups, fit the per-`(group, period)` stutter shape shrunk to the `θ_period` parent,
and run the `ε`-freeze check. **M3 milestone:** on a deliberately group-divergent
simulator, recover the per-group shapes; on a single protocol, collapse to one group.

## 2. Assumptions / decisions

- **Per-sample slip profiles added to D1** (`PrepassStats.slip_by_sample_period`) so
  D3 can re-aggregate the slip profile per sample group. Additive extension to the
  shipped D1+D2; the cohort `slip_by_period` parent is unchanged.
- **Clustering = deterministic union-find** over close neighbours in scaled
  `(ε, level)` space (no k-means / random init); the group **count** falls out of the
  distance threshold, ties resolve on the sample catalog index. Distances are deflated
  by the less-certain sample's precision (`√depth`).
- **Per-`(group, period)` shape** = the group's aggregated slip profile shrunk toward
  the `θ_period` parent via `refine_theta_locus` (reused from B2) — so a truly
  invariant group collapses to the shared shape and a divergent group keeps its own.
- **`ε`-freeze check** = a binomial BIC comparison of frozen-per-group `ε` vs.
  per-sample `ε`; the freeze stands when the richer model doesn't earn its parameters.
- **Deferred to F:** uncertainty-scaled distances beyond the basic depth weighting, a
  BIC-selected group *count* (here a calibrated threshold), the singleton-group
  shrinkage corner.

## 3. Changes made

- [prepass.rs](../../../../src/ssr/cohort/prepass.rs) — added
  `slip_by_sample_period` + `add_slip`; exposed `sample_eps` / `fit_level` /
  `merge_sample_stats`; split out `run_prepass_stats`.
- [sample_groups.rs](../../../../src/ssr/cohort/sample_groups.rs) (new) — `ClusterCfg`,
  `GroupedParams`, `cluster` (union-find), `group_samples`, `eps_freeze_check`.
- [mod.rs](../../../../src/ssr/cohort/mod.rs) — wired `sample_groups`.

## 4. Tests added (3)

- **M3** — `m3_two_protocol_cohort_recovers_per_group_shapes`: two protocols (distinct
  `ε`/level/decay) cluster into 2 groups; per-group `ε` separates in the injected
  direction; the per-`(group, period)` decays stay divergent (one shallow `< 0.2`, one
  heavy-tailed `> 0.3`).
- `single_protocol_cohort_collapses_to_one_group`.
- `frozen_epsilon_is_justified_for_a_single_protocol`.

## 5. Validation results

- `cargo fmt --check` → pass · `cargo clippy --all-targets --all-features -- -D warnings`
  → pass.
- `cargo test --all-features` → **1241 lib pass** (+3), integration + doctests green.
- ▶ **M3 holds:** per-group shapes recovered on a group-divergent cohort; single
  protocol collapses.

## 6. Tradeoffs and follow-ups

- Milestone D is now complete (D1+D2+D3). The pre-pass produces frozen `ε` (per group),
  the `θ_period` parent + per-`(group, period)` shapes, and the per-group level seeds.
- The grouped params are not yet wired into genotyping — **E1** consumes them (the
  outer `F` + per-group level loop), and **E2** adds allele-balance + full VCF.
- Calibration of the clustering threshold / shrink strength / `ε`-freeze penalty is F2.
