# General Dirichlet-multinomial SFS genotype prior — GIAB validation gate

**Date:** 2026-07-04
**Scope:** end-to-end GIAB per-sample validation of the **general
Dirichlet-multinomial** SFS genotype prior (the generalization of the
biallelic-only grid prior to all `(ploidy, n_alleles)` shapes), wired into the
posterior engine on branch `sfs-genotype-prior` (Steps G1–G4e). This is the gate
required **before** deleting the grid + HWE plug-in paths (Step G5).
**Verdict:** clean pass — matches/beats the grid at every depth, zero
precision/recall/FP cost → proceed to G5 deletion.

---

## What changed since the grid validation

The biallelic-only grid prior (`SfsGenotypePrior`, `e_step_sfs_biallelic`) is
replaced by the **Dirichlet-multinomial** closed form (`crate::genetics`) wired
into the standard `e_step` / `e_step_simd` for **every** shape:

- concentration `α_ref = 1`, `α_alt(a) = θ̂/(n_alleles−1)` from the measured
  cohort diversity θ̂ (the "C + choice-2" decision, arch §9.2);
- per-sample inbreeding `F` applied via the engine's existing Wright mixture
  (`fixation_index_overrides`, fed from `DiversityEstimate`);
- QUAL unchanged (separate Beta-pseudocount AF path).

The **key empirical question** this gate answers: the DM prior deliberately does
**not** reproduce the grid's inflated hom-ref weight (0.878) — it uses the
genetically-correct genome-wide variant fraction `3θ/2` (hom-ref ≈ 0.9985). The
concern was that the stronger hom-ref prior might drag weak-evidence
low-coverage samples toward hom-ref and regress the 5× concordance the grid won.
**It does not** (below).

## Setup

`POP_VAR_CALLER_BIN=target/release/pop_var_caller PRESET=high-recall
benchmarks/giab/src/run_ours_per_sample.sh {5,10,15,30}x` (HG002/3/4,
single-sample), DM prior live via the G4e pipeline flip. Host `cargo build
--release` (the script execs the host binary, not the container build).
Analysis: `tmp/gt_confusion_new.py` (GT concordance), `tmp/prec_recall.py`
(precision/recall/FP) — the same scripts as the grid validation.

## SNP genotype concordance (= 1 − GT-mismatch / TP)

| coverage | grid (prev) | **DM (new)** | freebayes | `1/1→0/1` overcall (new) |
|---|---|---|---|---|
| 5× | 94.6% | **94.9%** | 93.4% | 4 |
| 10× | 98.9% | **99.0%** | 98.3% | 6 |
| 15× | 99.7% | **99.8%** | 99.2% | 4 |
| 30× | 99.9% | **≈100%** (1 mismatch) | 99.4% | 2 |

The general DM prior **matches or slightly beats the grid at every depth** and
beats freebayes everywhere. The low-coverage het over-call stays fixed: the 5×
`1/1→0/1` transition is 4 (grid was 8). The residual 5× mismatch is dominated by
`0/1→1/1` (70) — true hets sequenced all-alt at low depth, the same
sampling-limited error freebayes has (75), not a prior artefact. **No
higher-depth regression.**

## Precision / recall / FP — unchanged

The genotype prior moves genotypes among *already-emitted* variants (het ↔
hom-alt both emit the ALT), so the allele-level metrics are identical to the
grid prior:

| cov | metric | grid (prev) | **DM (new)** |
|---|---|---|---|
| 5× | TP / FP / FN | 1455 / 13 / 606 | 1455 / **12** / 606 |
| 10× | TP / FP / FN | 1881 / 16 / 180 | 1881 / **14** / 180 |
| 15× | TP / FP / FN | 1978 / 8 / 83 | 1978 / **7** / 83 |
| 30× | TP / FP / FN | 2034 / 14 / 27 | 2034 / **12** / 27 |

TP and FN are identical at every depth; FP is within ±2 (pre-existing QUAL
non-determinism). SNP F1 beats freebayes at every depth (5×: 0.825 vs 0.738;
30×: 0.991 vs 0.962).

## The C + choice-2 decision is validated

The DM prior's hom-ref weight (0.9985) is the genetically-correct genome-wide
value, ~80× smaller variant fraction than the grid's un-normalised-sum artifact
(0.122). GIAB confirms this is **correct, not too conservative**: 5× concordance
did not regress (94.9% ≥ grid 94.6%). The reserved fallback (an independent
monomorphic-mass term to lower hom-ref without touching α_alt) is **not needed**.

## General shapes (multiallelic / polyploid)

GIAB is biallelic-diploid only, so it does not exercise the shapes the grid could
not cover. Those are validated by the engine unit tests
(`dirichlet_prior_matches_pure_value_{triallelic_diploid,tetraploid_biallelic}`,
which pin the engine posterior to `genetics::dirichlet_multinomial_log_priors`),
and a synthetic multiallelic/polyploid end-to-end pass is scheduled for G6.

## Conclusion

The general Dirichlet-multinomial SFS prior is a strict improvement on the grid:
matches/beats SNP genotyping at every depth, zero precision/recall/FP cost, and
covers every shape. **Go for G5** — delete `e_step_sfs_biallelic`, the grid
`SfsGenotypePrior` module, the `sfs_prior_tables` config path, and the HWE(p̂)
plug-in genotype prior. Owed in G6: a tomato (inbred, `F>0`) sanity pass +
synthetic multiallelic/polyploid end-to-end tests.

Reproduce: `tmp/gt_confusion_new.py`, `tmp/prec_recall.py` (host release build,
`PRESET=high-recall`).
