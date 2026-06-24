# SSR Stage-2 read-likelihood model bake-off — implementation report

**Date:** 2026-06-24
**Branch:** `ssr-cohort`
**Plan:** [ssr_stutter_scoring_model_bakeoff.md](../../implementation_plans/ssr_stutter_scoring_model_bakeoff.md)
**Background:** [ssr_delimiter_gap_penalty_2026-06-24.md](../research/ssr_delimiter_gap_penalty_2026-06-24.md)

## 1. Plan

Decide how Stage-2 (`ssr-call`) computes the read likelihood `Qᵣ(obs | candidate)` — the
probability an observed repeat tract arose from a candidate allele via stutter (length
change) + base error — by implementing three candidate models behind one swappable trait
and **baking them off on synthetic data with known truth**, then productionizing the
winner and removing the losers. The three models:

- **A — HipSTR**: explicit in-frame + out-of-frame geometric stutter (mirrors
  `StutterModel::log_stutter_pmf`).
- **B — fix-current**: the existing `Σ_Δ S_θ(Δ)·align_subst` decomposition.
- **C — two-penalty pair-HMM** (the user's idea): a single forward DP with a cheap
  whole-motif-unit gap and a stiff single-base gap applied uniformly along the alignment.

Executed in the plan's order of work: trait + Model B refactor → harness + metrics on G1
→ Model C → Model A → G2/G3 + full 3×3 sweep → decide + productionize.

## 2. Assumptions / declared choices

- **Bake-off scores on truth-derived parameters**, not the pre-pass estimator
  ([`param_set_from_truth`](../../../../src/ssr/cohort/bakeoff.rs)), so a model is judged on
  the *shape* of its `Qᵣ`, not on estimator recovery (plan §9). Each model's parameters are
  derived from the *same* `ReadScoringContext`, so the comparison is of models, not
  estimators.
- **A's out-of-frame mass** has no in-frame-context counterpart, so it is a fixed small
  fraction (`OUT_FRAME_REL = 0.05`) of the in-frame slip mass — a declared estimator for
  the bake-off. A real per-period out-of-frame estimator from the pre-pass is the
  productionization follow-up (see §6).
- **C's gap costs** are derived from `level` (split by the shape's direction bias) and `ε`;
  first-class frozen `g_unit`/`g_base` were never reached (C lost).
- **The bake-off relaxes the locus-admission periodicity gate** so the read-likelihood
  models are exercised on messy reads (see §5 finding), since that gate is upstream of
  `Qᵣ` and would otherwise filter the G2/G3 loci identically for every model.
- **Per-period parent shape seed**: when truth carries divergent per-group shapes,
  `param_set_from_truth` seeds each period from the first group's shape (single-group
  cohorts have no divergence).

## 3. Changes made

- **New `read_model/` module** ([mod.rs](../../../../src/ssr/cohort/read_model/mod.rs)):
  the `ReadLikelihoodModel` trait + per-call `ReadScoringContext` (motif, θ, level, ε);
  model-frozen parameters live on each model value.
- **Model B** ([classic.rs](../../../../src/ssr/cohort/read_model/classic.rs)):
  `ClassicStutterModel`, a pure adapter over the existing `read_likelihood`. Now
  **test-only**, retained as the bake-off's reference comparator.
- **Model A** ([hipstr.rs](../../../../src/ssr/cohort/read_model/hipstr.rs)):
  `HipstrModel`, mirroring `log_stutter_pmf` (in/out-of-frame split), reusing B's
  `reach_variants` + `align_subst` for the substitution-only align. **Now the production
  model.**
- **Model C** (removed): `two_penalty.rs` — a period-aware two-penalty forward DP with
  self-affinity normalization. Implemented, baked off, eliminated, and deleted.
- **EM threading** ([em.rs](../../../../src/ssr/cohort/em.rs)): `run_locus_em_with` /
  `compute_data_ll` are generic over `M: ReadLikelihoodModel`; `read_likelihood` became
  `model.q_r(...)` with identical inputs. Production call sites
  ([driver.rs](../../../../src/ssr/cohort/driver.rs),
  [inbreeding.rs](../../../../src/ssr/cohort/inbreeding.rs), the `run_locus_em` wrapper)
  use `HipstrModel`.
- **Generative axis** ([sim.rs](../../../../src/ssr/cohort/sim.rs)): `GenerativeNoise`
  (out-of-frame rate, impurity rate, extra substitution) + `simulate_with`; G1 =
  `none()`, G2 = `hipstr_like()`, G3 = `messy()`. Truth genotypes are unchanged by noise.
- **Bake-off harness** ([bakeoff.rs](../../../../src/ssr/cohort/bakeoff.rs)):
  `run_bakeoff_with`, `BakeoffMetrics` (genotype concordance, allele-length-error
  distribution, Expected Calibration Error, scoring wall), and `format_results_table`.

## 4. Results — the decision table

3×3 scoring × generative, truth-derived params, 16-sample / 3-locus sweep cohort
(`sweep_spec`, seed 2026), depth 60. Scoring ms is debug-build (the *ratio* is the model
cost):

```
model × generative    concordance  mean|Δ|  max|Δ|     ECE    scoring ms
B-classic / G1-clean       1.0000   0.0000       0  0.0000        1740
A-hipstr  / G1-clean       1.0000   0.0000       0  0.0000           9
C-2penalty/ G1-clean       0.6667   0.3229       2  0.2981         339
B-classic / G2-oof         1.0000   0.0000       0  0.0000        2367
A-hipstr  / G2-oof         1.0000   0.0000       0  0.0000          10
C-2penalty/ G2-oof         0.0208   5.8750      11  0.9769         682
B-classic / G3-messy       1.0000   0.0000       0  0.0166        8522
A-hipstr  / G3-messy       1.0000   0.0000       0  0.0001          34
C-2penalty/ G3-messy       0.0417   5.5000      11  0.9382        1948
```

Depth sweep on G3-messy: A and B both ≈0.98–1.0 by depth 20; C ≈0.05 at every depth.

**Decision (plan §7 priority: G3 accuracy → calibration → runtime → maintenance):**

- **C eliminated.** Folding stutter into the alignment lets out-of-frame reads pull the
  genotype. Self-affinity normalization (dividing by `forward(cand|cand)`) fixed an initial
  catastrophic long-allele bias — without it C called every genotype as the longest allele
  — but a residual +1 bias on clean data and a full collapse on G2/G3 remained.
- **A (HipSTR) chosen over B.** Both tie at perfect concordance across G1/G2/G3; A is far
  better calibrated (G3 ECE 0.0001 vs 0.017) and ~100–250× cheaper per `Qᵣ` (A scores one
  `bp_diff` term + a closed-form substitution; B does a 21-way Δ-sum with a banded DP per
  term — making B, the incumbent, the *slowest* model). Since `Qᵣ` is ~75% of pileup
  self-time, that runtime gap is material.
- **Productionized A; removed C; kept B** behind the trait as the bake-off baseline.
  Swapping production from B to A required **zero test re-baselining** — the two agree on
  every existing cohort/integration test (both accurate; calls + VCF byte-identical).

## 5. Tests added/updated

- **Trait/Model B**: `q_r == read_likelihood` bit-for-bit (the behaviour-preserving refactor).
- **Model A** (8): a direct `stutter_pmf` ↔ `log_stutter_pmf` term-by-term check; faithful
  dominated by `equal`; contraction-bias + geometric ordering; in-frame unit beats
  out-of-frame base; impure placements; out-of-frame resize; long allele; scratch reuse.
- **Model C** (while it existed): faithful normalizes to 1.0; long-candidate-bias removed
  by normalization; slip monotonicity; chained unit gaps; the in-frame-vs-out-of-frame
  property; impure / mononucleotide / long allele; scratch reuse. (Removed with C.)
- **Generative axis**: G2 produces non-unit lengths; G3 interruptions break the tiling.
- **Harness**: concordance/allele-error/calibration bookkeeping (`score_genotype`); ECE
  zero for a perfectly-calibrated set; metric determinism; the survivor sweep + depth
  robustness emitting the table.

## 6. Validation

- `cargo fmt --check`: clean.
- `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --lib`: **1325 passed, 0 failed, 2 ignored.**

## 7. Tradeoffs and follow-ups

- **A's out-of-frame parameter is a fixed `OUT_FRAME_REL` constant** (the only declared
  bake-off estimator). **Follow-up:** estimate out-of-frame mass per period from the
  pre-pass (out-of-frame reads binned separately, as HipSTR does), replacing the constant.
- **Separate Stage-1.5 finding (not this bake-off's scope):** the locus-admission
  periodicity gate ([`is_periodic`](../../../../src/ssr/cohort/candidate_set.rs)) is
  **presence-based** — it counts distinct observed lengths with no support floor, so a
  handful of out-of-frame reads injecting singleton odd lengths trips `NotPeriodic` and
  filters an otherwise-callable locus. The bake-off relaxed it to exercise the `Qᵣ` models;
  making it support-weighted is a worthwhile robustness fix.
- **Runtime numbers are debug-build.** The model-cost *ratio* (A ≪ C < B) is structural and
  trustworthy; absolute timings need a release bench (none exists for the SSR path yet).
- **Model B retained** behind the trait as the reference comparator; if future swaps stop
  needing it, B + `read_likelihood` + `s_theta` can be inlined/removed.
- **Real-data calibration** of absolute parameter values remains the next effort (plan §1
  scoped this bake-off to model *shape* on synthetic data).
