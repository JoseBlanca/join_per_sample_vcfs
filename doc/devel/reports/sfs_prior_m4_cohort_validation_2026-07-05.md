# SFS genotype prior — Milestone 4 (empirical-Bayes leave-one-out cohort prior) validation

**Date:** 2026-07-05
**Branch:** `sfs-genotype-prior` (tip `ae6aa03`)
**Scope:** validate the leave-one-out (LOO) cohort genotype prior on GIAB (single-sample,
truth-based) and tomato2 (59-sample real cohort, no truth). Decide go/no-go for merge.

**Outcome:** **M4 accepted** (owner sign-off 2026-07-05). Single-sample is byte-identical to
the shipped G5 DM prior; the cohort behaviour is a defensible recall change, not a demonstrated
regression. Three attempted refinements (F-aware first step, multi-start, read-based frequency)
were **tried and reverted** — the investigation showed they were solving a mis-diagnosed
problem.

---

## 1. GIAB (single-sample, truth-based) — PASS, byte-identical

The LOO prior collapses to the species prior when there is exactly one sample
(`E[cohort counts] − E[own] = 0`), so `dispatch_e_step` routes `n_samples == 1` through the
unchanged G5 path. GIAB is a per-sample benchmark, so it is **byte-identical to G5 by
construction**, and the run confirms it:

| depth | SNP GT concordance | precision | recall | FP |
|---|---|---|---|---|
| 5× | **94.9 %** (74/1455 mismatch) | 0.992 | 0.706 | 12 |
| 10× | 99.0 % | 0.993 | 0.913 | 14 |
| 15× | 99.8 % | 0.996 | 0.960 | 7 |
| 30× | 100 % | 0.994 | 0.987 | 12 |

Matches the G5 baseline exactly; beats freebayes on recall/F1 at every depth. **No
single-sample regression.**

## 2. tomato2 (59-sample cohort, no truth) — cohort het increase, investigated

Same 59 `.psp`, same flags, only the engine differs (G5 baseline binary built from the parent
commit `a219726` in a worktree). Under **standard filtering** (default `--min-qual`,
allele-balance filter, DUST, coverage-paralog filter all on), PASS records:

| metric | G5 | M4 | M4/G5 |
|---|---|---|---|
| sites emitted | 229,585 | 238,580 | +3.9 % |
| het calls, PASS (all GQ) | 216,863 | 682,526 | **3.15×** |
| het calls, PASS + GQ≥20 | 142,381 | 229,360 | **1.61×** |
| het calls, PASS + GQ≥30 | 105,128 | 183,365 | **1.74×** |

The raw ~3× tripling is dominated by a low-GQ tail; at a normal genotype threshold (GQ≥20) the
increase is **1.6×**. The investigation into *what* that increase is:

### 2.1 It is not frequency over-estimation

The initial hypothesis (the cohort EM inflates the allele frequency via genotype→prior
feedback) was **refuted by measurement**. The **called genotypes carry *less* alt than the
reads contain**: at DP≤5, `AC/AN = 0.075` vs read-fraction `0.094` (ratio 0.80×). The
genotyper is conservative relative to the reads — the alt is genuinely in the reads
(~9 %, far above the ~1 % error floor), and the genotyper calls less het than that would
justify. So this was never a genotype-prior bug.

### 2.2 The low-GQ, M4-only tail is correlated-mismapping artifact

The 11,081 sites M4 emits that G5 does not carry a clear **mismapping fingerprint**: mean
MQDiff **−10.2** (vs −3.4 shared), mean MQDiffT **−8.2** (vs −2.0), 56 % with MQDiff<−3, and
slightly elevated coverage (401 vs 369). These are reads mismapped from elsewhere (alt reads
map far worse than ref). Cohort pooling **amplifies** them because mismapping is *correlated*
across samples (same repeat structure → same mismap in every sample), so it accumulates like a
real shared variant — the "correlated errors" caveat, at cohort scale. **Note:** this is a
*MAPQ* class, **not** the coverage/het-excess class the hidden-paralog filter scores
(`PARALOG_POST` was only 0.037 on these sites) — the coverage-paralog filter does not catch
them. They are, however, **low-GQ** and are removed by the GQ≥20 filter.

### 2.3 The confident (GQ≥20) excess is not artifact — plausible recall

A post-hoc MAPQ-diff-t sweep does **not** close the confident-het gap:

| MQDiffT cut | G5 | M4 | M4/G5 |
|---|---|---|---|
| none | 142,381 | 229,360 | 1.61× |
| ≥ −3 | 110,823 | 170,006 | 1.53× |
| ≥ 0 (alt maps ≥ ref) | 56,446 | 88,278 | 1.56× |

Even keeping only sites with **no** mismapping signal (MQDiffT ≥ 0), M4 is still 1.56× G5. So
the confident excess is **not** the mismapping class — it is at sites **shared** with G5, where
M4 confidently recovers additional low-depth carriers G5 left as hom-ref. That is
strength-borrowing working as designed (the cohort frequency sharpens a low-depth carrier's
posterior to a confident het). It reads as a **recall improvement**; it cannot be *proven* real
without a cohort truth set or an independent caller (not available here).

## 3. Refinements tried and reverted

All three were built on the (later-refuted) premise that the cohort frequency was over-estimated:

- **F-aware neutral first EM step** — kept inbreeding F on iteration 1. Moved the het count by
  0.005 %. The inflation was steady-state, not an init artifact.
- **Multi-start + evidence selection** — ran the EM from a conservative and a permissive seed,
  kept the higher-evidence fixed point (data fit + SFS frequency penalty). Reduced hets ~1 %
  (661 k → 654 k). The SFS location penalty (~5 log-units) is dwarfed by the multi-sample data
  fit, so it flipped only ~5 k of ~14 k sites.
- **Read-based frequency** (frequency from observed allele read counts, not genotypes) — a
  measurement pre-check killed it before implementation: read-AF (0.094) is *higher* than the
  genotype frequency at low depth, so it would produce *more* hets, and raw read counts lose the
  base-quality weighting the likelihoods provide.

All three reverted to the committed `ae6aa03`; the branch ships the plain LOO prior.

## 4. Decision and residual risk

**Ship M4.** Rationale: clean where we have truth (single-sample byte-identical); the cohort
low-GQ artifacts are correlated-mismapping and filter out on GQ; the confident cohort excess
survives both the GQ and MAPQ-diff filters and reads as strength-borrowing recall.

**Residual risk (documented, not blocking):** the confident cohort het excess (~1.6×) is
*plausibly* real recall but *unproven* — a subtler artifact cannot be excluded without a cohort
truth set. Follow-ups if desired: an independent-caller (freebayes/GATK) cross-check on the same
cohort; a MAPQ-diff *filter* (only annotated today, not applied) for repeat-rich cohorts that
want the low-GQ mismapping tail gone at the source.

## 5. Commands

- GIAB: `POP_VAR_CALLER_BIN=target/release/pop_var_caller PRESET=high-recall bash benchmarks/giab/src/run_ours_per_sample.sh {5,10,15,30}x` → `uv run tmp/gt_confusion_new.py` + `uv run tmp/prec_recall.py`.
- tomato2: `POP_VAR_CALLER_BIN=… OUT=…/cohort_m4_std.vcf.gz` via `var-calling --reference … --regions … <psps>` (standard filters); G5 baseline from the `a219726` worktree binary.
- Analysis: `bcftools isec` (G5 vs M4), `bcftools query` on `AC/AN`, `AD`, `MQDiff`, `MQDiffT`,
  `PARALOG_POST` (host, regenerable; result VCFs under `benchmarks/tomato2/results/`).
