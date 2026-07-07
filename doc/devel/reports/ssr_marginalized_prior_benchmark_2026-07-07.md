# Benchmark: SSR marginalized DM prior vs plug-in (Phase 3.5)

**Date:** 2026-07-07 · **Branch:** `em-convergence-criterion` · **Verdict:
INCONCLUSIVE on tomato — HipSTR is not a fair reference here (see §Caveat). Keep
the marginalized prior opt-in; decision deferred to an outbred (human) benchmark.**

## What was tested

Phase 3 gave the SSR caller the SNP path's improved *marginalized + leave-one-out
Dirichlet-multinomial* genotype prior (seeded by SSR's own mode-centred `G₀`),
behind the `EmCfg.marginalized_prior` toggle (env `PVC_SSR_MARGINALIZED_PRIOR=1`).
Phase 3.5 benchmarks it against the current plug-in HWE(π) prior on **ssr_tomato1**
(51-sample *S. lycopersicum* cohort, ~15k SSR catalog), scored as genotype
concordance vs **HipSTR** on HipSTR's whole-unit length-polymorphic loci (no truth
set; HipSTR is the reference).

Both runs use the **same** catalog + `.ssr.psp` inputs; only `ssr-call` differs.

## Result

| metric (vs HipSTR, 1645 HipSTR length-polymorphic loci) | plug-in (default) | marginalized |
|---|---:|---:|
| total loci emitted | 74 573 | 73 993 |
| HipSTR-polymorphic loci **also emitted by ours** | **774 (47.1%)** | **539 (32.8%)** |
| HipSTR-polymorphic loci we call mono/no-call/filtered | 871 | **1106** |
| genotype concordance on shared sample cells | 96.5% (27899/28914) | 96.1% (11186/11639) |

**The marginalized prior drops ~30% of the comparable (HipSTR-polymorphic) loci
(774 → 539) while concordance is flat (96.5% → 96.1%).** It calls 235 more of these
loci length-monomorphic.

## Caveat — HipSTR is not a fair reference for a prior change on a selfer

Whether that extra conservatism is *correct* or *over-shoot* cannot be judged
against HipSTR here, because **HipSTR has no inbreeding (`F`) correction and tomato
is an extreme selfer (mean F_IS ≈ 0.82).** A biologically-correct caller on a selfer
should collapse most heterozygotes to homozygotes and call fewer loci polymorphic;
HipSTR, being `F`-blind, over-calls polymorphism/heterozygosity. So "recall vs
HipSTR" structurally penalises any `F`-aware prior, and the loci the marginalized
prior "misses" may be HipSTR false-polymorphisms. Both of our priors are `F`-aware
(the plug-in `genotype_prior` also takes `f`); the marginalized one is simply *more*
conservative, and against an `F`-blind reference that difference is unadjudicable.

**Conclusion: the tomato/HipSTR result is inconclusive, not a refutation.** The
decision needs an **outbred** cohort (F ≈ 0), where the `F` correction is inert and
the marginalization effect can be judged on its own — the human benchmark (below).

## Why the two priors diverge here (mechanism, still relevant)

- **The `G₀`-as-DM-concentration scale is the suspect** (the Phase-1 §Q3 risk we
  flagged to watch). `G₀ = p^|Δ|` was tuned as a *plug-in regularizer* (small
  counts nudging `π`); as a DM *concentration* the same small values make the
  frequency belief very diffuse, so marginalizing pulls hard toward the mode and
  suppresses minority alleles → fewer loci called polymorphic.
- **High inbreeding amplifies it.** The cohort's mean F_IS ≈ 0.82 (tomato is a
  selfer; `ssr-call` warns on it). The Wright-`F` branch strongly favours
  homozygotes; combined with the diffuse marginalized prior, borderline
  polymorphic loci collapse to monomorphic and drop out of the variable-only
  emission.
- This is the regime Q-G2 and Phase-1 predicted: the marginalize-vs-plug-in fix
  helped SNPs at *low coverage / small cohorts*; SSR loci here are multi-allelic
  and well-covered, so the fix doesn't help and the diffuse-prior side effect
  dominates.

Not a bug: the prior math was verified identical to the SNP engine and the loop
to the SNP `e_step_cohort_loo` (Phase 3.2/3.4 reviews); the unit tests (incl.
LOO-moves-a-sample and F-threading) pass. This is a genuine model property.

## Decision (deferred)

- **Verdict deferred, not refuted.** The default stays plug-in for now, but *only*
  because the outbred test that could justify flipping it hasn't run — not because
  tomato showed the marginalized prior is worse.
- The marginalized path stays **opt-in** behind `PVC_SSR_MARGINALIZED_PRIOR=1`.
  SSR production is unchanged (toggle off by default; plug-in byte-identical).
- **Next test — outbred human data (F ≈ 0).** Interim: the GIAB trio (thin — ~100
  regions, 3 samples), where the key diagnostic is whether the plug-in-vs-marginalized
  divergence *collapses* relative to tomato (which would confirm the tomato gap was
  the selfer/`F` effect). Definitive: a larger outbred human SSR cohort (in
  preparation) so the cohort leave-one-out prior is actually exercised.
- Secondary open item if it is ever revived at scale: the `G₀`→DM-concentration
  mapping (§Q3) may still want re-tuning.

## GIAB interim attempt (human trio) — unusable, too thin

Ran the SSR pipeline on the GIAB HG002/3/4 trio (whole-genome GRCh38 catalog =
515 352 loci; the M5-reheadered trio BAMs) to get an *outbred* read where the `F`
confound should vanish. It does **not** work as a benchmark:

- Only **20 loci are genotyped** (PASS). The trio BAMs are sliced to ~100
  SNP-benchmark regions, so 515 260 / 515 284 emitted records are `lowDepth`
  (no reads).
- **F_IS is mis-estimated at 0.99** — a pre-pass artifact of the 515k empty
  genome-wide loci plus low-depth het-dropout on the handful of covered ones. So
  this is *not* the clean `F ≈ 0` test; it carries the same high-`F` confound as
  tomato.
- On the 20 genotyped loci, plug-in and marginalized make **identical genotype
  calls** (16/20 identical; 4 differ only in GQ). A weak positive that the
  marginalized prior does not change human genotype calls here — but 20 loci with
  a broken `F` estimate cannot decide anything.

**The cohort SSR caller needs a real cohort of many covered polymorphic loci
(for the `F`/stutter/`G₀` pre-pass and for concordance). The GIAB 100-region trio
cannot provide that.** The decision genuinely awaits the dedicated outbred human
SSR benchmark (in preparation).

## Reproduce

```
PVC_SSR_MARGINALIZED_PRIOR=1  # (or unset for plug-in)
ssr-call --catalog ssr_tomato1.ssr.catalog --output out.vcf --threads 4 <*.ssr.psp>
benchmarks/lib/ssr_concordance.py --ours out.vcf --hipstr <hipstr.vcf.gz>
```
Inputs: `benchmarks/ssr_tomato1/results_ssr15k/` (psp + HipSTR); catalog at
`benchmarks/ssr_tomato1/results/ours/ssr_tomato1.ssr.catalog`.
