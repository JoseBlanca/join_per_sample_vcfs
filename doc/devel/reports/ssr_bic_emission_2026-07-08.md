# SSR emission by model selection (BIC) — in-caller implementation + benchmark

**Date:** 2026-07-08 · **Branch:** `ssr-bic-emission` (off `ssr-drop-fate`) ·
**Bench:** `benchmarks/ssr_tomato1` (51-sample *S. lycopersicum* selfer, median ~3 reads/plant)

## What this is

The SSR caller's emission decision was a heuristic chain — genotype each plant,
then *emit-iff-variable* plus an allele-balance no-call — that leaks false
positives (systematic per-chemistry stutter mis-read as an allele) and can't be
tuned cleanly. This replaces it, behind a toggle, with a **principled model
selection**: per locus, is the cohort's reads better explained by a *monomorphic*
model (one fixed allele + stutter) or a *polymorphic* one (an alt allele
segregates)? Emit iff the polymorphic model wins by **BIC**.

This is the in-caller realisation of the offline prototype
(`benchmarks/ssr_tomato1/scripts/bic_prototype.py`), now using the caller's real
per-chemistry-anchored stutter read model instead of a toy kernel. It is being
compared against a freebayes-style SFS-prior variant on a sibling branch, scored
on the identical silver-standard core.

## How it works (design)

- The EM already computes, per plant, the read log-likelihood of every diploid
  genotype (`data_ll`, from the Qᵣ stutter model — with the per-locus
  `refine_theta_locus` shape and level baked in). We reuse it; no new read model.
- `run_locus_em_with` now also returns `EmissionEvidence`
  ([em.rs](../../../src/ssr/cohort/em.rs)):
  - `ln_marginal` — Σ over plants of `log Σ_G P(G | π̂, F) P(reads | G)`, the
    maximised marginal likelihood under the polymorphic model (genotypes
    integrated out, frequency fit by the EM, inbreeding `F` in the prior).
  - `ln_monomorphic` — Σ of `P(reads | hom aa)` for the single fixed allele that
    best explains the cohort.
- In `genotype_locus` ([driver.rs](../../../src/ssr/cohort/driver.rs)), when
  `PVC_SSR_EMIT_MODEL=bic` (was `PVC_SSR_BIC_EMIT=1` before the emission models were
  unified behind one selector; see [freebayes spec](../specs/ssr_freebayes_marginal_emission.md)):
  emit iff `2·(ln_marginal − ln_monomorphic) > n_alt·ln(N)
  + 2·margin` (BIC, `N` = present plants), and only with a variant genotype so no
  monomorphic row is ever written. FP-control is **not** run — the model choice
  replaces it. `PVC_SSR_BIC_MARGIN` trades precision for recall. Default
  (`PVC_SSR_EMIT_MODEL=heuristic`) → **byte-identical** (verified by `diff`). 243
  `ssr::cohort` tests pass, plus a
  unit test that `EmissionEvidence` separates a monomorphic from a segregating locus.

## Benchmark (silver-standard confident core: 561 true100 / 8850 false100)

Recall = confident-real loci emitted; FP rate = confident-monomorphic loci emitted.

| emission model | recall | precision | FP rate |
|---|--:|--:|--:|
| baseline heuristic gates | 80.7% | 97.0% | 0.16% |
| BIC (plug-in EM), margin 0 | 83.4% | 94.2% | 0.33% |
| BIC (plug-in EM), margin 10 | 76.3% | 96.4% | 0.18% |
| **MARG prior + BIC, margin 0** | **87.3%** | 94.0% | 0.35% |
| MARG prior + BIC, margin 5 | 83.8% | 95.7% | 0.24% |
| **MARG prior + BIC, margin 10** | **81.1%** | 96.8% | **0.17%** |
| MARG prior + BIC, margin 20 | 75.0% | 98.6% | 0.07% |

## Read

- **The marginalized prior is essential to the combination.** Plain BIC caps at
  ~83% recall because it still emits the EM's MAP genotypes, so it cannot recover
  the loci the plug-in EM collapses to homozygous-reference. The marginalized
  leave-one-out prior re-genotypes those as variant; BIC then decides emission.
  MARG+BIC reaches 87% recall — 6+ points above the heuristic — where plain BIC
  and the heuristic both plateau.
- **At the heuristic's operating point it matches it** (MARG+BIC margin 10 =
  81.1% / 0.17% ≈ baseline 80.7% / 0.16%) — but as **one principled test with a
  clean knob**, not a stack of hand-tuned gates. And it is a *curve*: the same
  model reaches 87% recall (margin 0) or 98.6% precision (margin 20), operating
  points the fixed heuristic cannot express.
- **It does not dominate the heuristic on the raw trade-off** at this depth — the
  median-3-reads floor is real, and the residual true/false overlap is genuine
  low-depth ambiguity no model erases. The value is principle + tunability +
  the recovered-recall ceiling, not a free lunch.

## Verdict

The BIC emission model is a **competitive, principled replacement** for the
heuristic gate stack: it matches the baseline operating point, extends the
reachable recall/precision curve, and removes the ad-hoc corroboration knobs. Its
best form is **MARG prior + BIC**. Whether it or the freebayes-style SFS-prior
variant becomes the caller's emission model is the next comparison; both are
scored on this same core with `silver_standard.py`.

## Reproduce

```
# BIC emission on (reuses existing .ssr.psp); sweep the margin. The emission model is
# now selected by PVC_SSR_EMIT_MODEL = heuristic (default) | bic | freebayes; orthogonal
# to PVC_SSR_MARGINALIZED_PRIOR (which shapes the EM genotypes every model consumes).
PVC_SSR_MARGINALIZED_PRIOR=1 PVC_SSR_EMIT_MODEL=bic PVC_SSR_BIC_MARGIN=10 \
  pop_var_caller ssr-call <psp>/*.ssr.psp --catalog <cat> --output out.vcf
# score on the confident core
uv run --no-project benchmarks/ssr_tomato1/scripts/silver_standard.py --score out.vcf
```
