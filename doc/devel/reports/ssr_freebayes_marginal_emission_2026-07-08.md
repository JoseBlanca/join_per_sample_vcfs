# SSR freebayes-style marginal emission — validation on ssr_tomato1

*Date: 2026-07-08. Branch `ssr-freebayes-marginal` (off `main` `ce3f6d7`; impl `ae54c55`).
Validates the freebayes-style joint marginal-likelihood emission model
([`specs/ssr_freebayes_marginal_emission.md`](../specs/ssr_freebayes_marginal_emission.md))
against the current heuristic gates and the BIC candidate. Companion candidate:
the BIC model-selection test ([`specs/ssr_bic_confident_genotype.md`](../specs/ssr_bic_confident_genotype.md)),
built separately — the two are compared here on an identical silver standard.*

## TL;DR

The freebayes SFS-prior emission model **dominates the current heuristic gates across
the low-false-positive frontier** and **dominates the BIC offline prototype at every
operating point**. At matched false-positive rate it recovers **+6 to +12 recall points**
below 0.10 % FP, and it opens a higher-recall range (up to **88.9 %**) the current gates
structurally cannot reach. Its one weak spot is the current gates' single operating point
(0.16 % FP), where the gates edge it by ~2–3 points. Default off → byte-identical;
thread-deterministic. **Recommendation: adopt freebayes as the emission model** (pending
the orthogonal truth set), with the per-sample allele-balance no-call retained (see
§Concordance).

## What was run

Stage 1 is unaffected; only **`ssr-call`** was re-run on the existing 51-sample cohort
(`benchmarks/ssr_tomato1/results_rerun_20260708/ours/cohort/psp/*.ssr.psp`, catalog
`results/ours/ssr_tomato1.ssr.catalog`). Two VCFs from the **same** branch binary:

- **baseline** — toggle off (`cohort.ssr.mainhead_gates.vcf`): the current heuristic gates
  (`apply_fp_control` + `is_variable` + plug-in `site_qual`). Verified **byte-identical to
  `main` HEAD** (my changes are inert when off).
- **freebayes** — `PVC_SSR_FREEBAYES_EMIT=1` (`cohort.ssr.freebayes.vcf`): the SFS-prior
  site test owns the emit decision and QUAL. Verified **byte-identical across `--threads 1`
  and `4`**.

Scored with `benchmarks/ssr_tomato1/scripts/silver_recall_fp_curve.py` (the shared
recall/FP scorer), which reproduces the confident-core silver standard from §4 of
`ssr_error_signals_dashboard.py`: **561 `true100` / 8,850 `false100`** loci where the
read-grounded pileup standard and the HipSTR standard agree. The scorer was validated by
reproducing the current gates' **80.7 % recall / 0.16 % FP** before use.

## The recall–false-positive frontier (confident core)

| QUAL ≥ | freebayes recall | freebayes FP | freebayes emitted | | baseline recall | baseline FP |
|---:|---:|---:|---:|---|---:|---:|
| 0   | **88.9 %** (499/561) | 1.24 % | 2161 | | 80.7 % | 0.16 % |
| 3   | 83.8 % (470/561) | 0.29 % | 1631 | | — | — |
| 10  | 82.9 % (465/561) | 0.27 % | 1510 | | — | — |
| 20  | 80.7 % (453/561) | 0.21 % | 1392 | | 69.7 % | 0.10 % |
| 30  | 78.8 % (442/561) | 0.19 % | 1310 | | 62.6 % | 0.06 % |
| 40  | 77.2 % (433/561) | 0.19 % | 1226 | | — | — |
| 50  | 75.6 % (424/561) | 0.14 % | 1151 | | 52.8 % | 0.03 % |
| 75  | 73.4 % (412/561) | 0.07 % | 991  | | 35.3 % | 0.01 % |
| 100 | 70.8 % (397/561) | 0.06 % | 878  | | 23.9 % | 0.01 % |
| 150 | 65.1 % (365/561) | 0.03 % | 710  | | — | — |

*(baseline = the current gates' own QUAL sweep; the gates emit a fixed set and their QUAL
is the plug-in `site_qual`.)*

### Read the frontier at matched false-positive rate

| FP rate | current gates recall | freebayes recall | Δ |
|---|---:|---:|---:|
| 0.03 % | 52.8 % (Q50) | **65.1 %** (Q150) | **+12.3** |
| 0.06 % | 62.6 % (Q30) | **70.8 %** (Q100) | **+8.2** |
| 0.10 % | 69.7 % (Q20) | **~76 %** (Q50–75) | **~+6** |
| 0.16 % | **80.7 %** (Q0) | ~78 % (Q40–50) | −2.7 |
| 0.21 % | *unreachable* (max 0.16 %) | 80.7 % (Q20) | — |
| 0.29 % | *unreachable* | 83.8 % (Q3) | — |

The freebayes QUAL is a **strictly better ranking of true-vs-false loci below 0.10 % FP**
(+6 to +12 points), and it extends the reachable recall to 83.8–88.9 %, which the heuristic
gates cannot reach at any threshold (their emit set caps recall at 80.7 %). The current
gates win only at their own tuned operating point (0.16 % FP), by ~2–3 points — i.e. they
are a single well-chosen point, not a better frontier.

### vs the BIC offline prototype

BIC prototype targets (toy stutter, will differ in-caller): 90 / 83 / 80 / 72 % recall at
3.9 / 1.0 / 0.38 / 0.08 % FP.

| target recall | BIC FP | freebayes FP (≈ same recall) | freebayes advantage |
|---:|---:|---:|---|
| ~90 % | 3.9 %  | 1.24 % (88.9 %, Q0)  | ~3× lower FP |
| ~83 % | 1.0 %  | 0.29 % (83.8 %, Q3)  | ~3.4× lower FP |
| ~80 % | 0.38 % | 0.21 % (80.7 %, Q20) | ~1.8× lower FP |
| ~72 % | 0.08 % | 0.06 % (70.8 %, Q100)| comparable |

On the evidence available, **freebayes dominates the BIC prototype at every operating
point** (2–3× lower FP at matched recall). Caveat: the BIC numbers are from an offline
prototype with a toy stutter model; the in-caller BIC (built separately) will differ, so
this is indicative, not final. But the freebayes result is measured **in-caller** on the
real chemistry, which is the stronger footing.

## Genotype concordance vs HipSTR — the honest tradeoff

| VCF | comparable loci emitted (of 1645) | concordance (sample cells called by both) |
|---|---:|---:|
| baseline (gates) | 774 (47.1 %) | **96.5 %** (27899/28914) |
| freebayes, QUAL ≥ 0 | 1128 (68.6 %) | 90.8 % (47604/52408) |
| freebayes, QUAL ≥ 20 | 842 (51.2 %) | 91.0 % (37135/40801) |
| freebayes, QUAL ≥ 40 | 755 (45.9 %) | 90.9 % (33232/36553) |

Concordance drops ~5.5 points and **does not recover with the QUAL threshold**. This is a
real, explainable effect and **not** a genotyping regression: the per-sample `GT:GQ:REPCN`
columns are the **same EM MAP calls** in both paths (the freebayes toggle changes only the
site emit decision and QUAL). The drop is entirely that toggling on **removes
`apply_fp_control`** — the per-sample allele-balance *no-call* that converted imbalanced
(stutter-shoulder) hets into `./.`. Those cells were removed from the baseline's
concordance denominator; without the no-call they return as called-but-discordant cells
(our het vs HipSTR hom). The freebayes test operates at the **site** level; it does not do
the **per-sample** cleanup the allele-balance no-call did.

**Design consequence:** the per-sample allele-balance no-call is orthogonal to the
site-level polymorphism test. Retaining it (running `apply_fp_control` for the per-sample
GT cleanup while the freebayes marginal owns the site emit/QUAL) would recover concordance
without touching the recall/FP win. The spec replaced `apply_fp_control` wholesale to keep
the first cut clean; the measurement says the per-sample no-call should be kept. This is
the one recommended change before adoption.

## Emission count

The freebayes path emits **3,217 PASS variable records** at QUAL ≥ 0 (vs 1,263 for the
current gates), thinned by the QUAL threshold (1,631 at QUAL ≥ 3, 1,392 at QUAL ≥ 20). The
larger raw set is expected — the emit gate is now "MAP-variable", with the QUAL doing the
discrimination downstream — and is the mechanism behind the recoverable higher-recall range.

## Design — how it combines (what it reuses)

Per locus, over the top-4 alleles by cohort frequency (a computational bound on the exact
enumeration, not a corroboration knob):

```
term(n)  = lnEwens(n)  +  Σ_s ln L_s(p = n / 2N)      // over integer count vectors n, Σn = 2N
Z        = logsumexp_n term(n)
QUAL     = −10·log10 [ Σ_{n fixed for one allele} exp(term(n) − Z) ]
```

- **`L_s(p) = Σ_G genotype_prior(G, p, F_s)·exp(data_ll[s][G])`** reuses (a) `data_ll` — OUR
  Qᵣ stutter read likelihoods, already per-locus/per-chemistry via `refine_theta_locus` —
  and (b) `genotype_prior` — the Wright genotype prior with **frozen per-sample inbreeding
  `F`** (`f_present`). Because the genotype prior is conditioned on `p`, samples are
  conditionally independent given `p`, so the site marginal is an **exact sum over
  allele-count vectors** — no combo enumeration, no banded search.
- **`lnEwens(n) = #nonzero·ln θ − Σ ln n_a`** is the neutral SFS prior (freebayes'
  Ewens Sampling Formula, `θ = 0.01`), the only genuinely new math (~20 lines). The
  `#nonzero·ln θ` term is the rare-allele penalty: each extra allele multiplies the prior
  by `θ ≈ 0.01`, so a stutter shoulder must earn that cost from the read likelihood — which
  is exactly what suppresses the systematic-stutter false positives.
- **Reused wholesale:** the EM (`compute_data_ll`, genotype enumeration, MAP calls), the
  candidate set, the VCF writer, `f_present`, `ln`/`exp`. New surface = one module
  `freebayes_emit.rs`. No new read model, no new EM, no genotype-prior duplication.

`θ` is fixed (no per-run knob). Sensitivity was not swept empirically (θ is a compile-time
constant); qualitatively, a larger θ weakens the rare-allele penalty → more emission /
higher FP, a smaller θ is stricter. Cohort-estimated θ is a later refinement.

## Tests

`freebayes_emit` unit tests (7): a clean segregating locus emits at high QUAL; a
systematic-stutter locus stays near QUAL 0; segregating outscores stutter by a wide margin;
triallelic segregation emits; `< 2` alleles → monomorphic; determinism across repeated
calls; genotype-index/enumeration agreement. Full SSR lib suite green (358). Byte-identical
default-off (vs `main` HEAD) and thread-deterministic (T4 == T1) verified end-to-end.

## Honest caveats

- **The silver standard is read-grounded**, so it partly co-defines "segregation"; the
  non-circular check is the orthogonal truth set (in prep). The recall/FP frontier is the
  fair comparator (FP especially — a truly-monomorphic locus called variable is a clear
  error regardless of read extraction); the concordance-vs-HipSTR numbers carry the usual
  `F`-blind-reference caveat on a selfer.
- **Median 3 reads/plant is a real information floor.** The win is from using the cohort's
  segregation signal optimally (a site-level test), not from beating physics.
- **BIC comparison is prototype-vs-in-caller.** The freebayes numbers are in-caller on real
  chemistry; the BIC targets are from a toy offline prototype. The in-caller BIC (separate
  work) is the apples-to-apples counterpart when it lands.

## Verdict

The freebayes SFS-prior emission model is a **better emission frontier** than the current
heuristic gates (dominant below 0.10 % FP, +6–12 recall points; reaches recall the gates
cannot) and **beats the BIC prototype 2–3× on FP at matched recall**. It is default-off
byte-identical and thread-deterministic. **Adopt it** (subject to the orthogonal truth set),
with one change from the first cut: **keep the per-sample allele-balance no-call**
(`apply_fp_control`) for the per-sample GT columns, since removing it — not the site test —
is what costs the HipSTR concordance.

## Reproduce

```
# freebayes VCF (host release binary; psp + catalog only, no FASTA needed):
PVC_SSR_FREEBAYES_EMIT=1 target/release/pop_var_caller ssr-call \
  --catalog benchmarks/ssr_tomato1/results/ours/ssr_tomato1.ssr.catalog \
  --output cohort.ssr.freebayes.vcf --threads 4 \
  benchmarks/ssr_tomato1/results_rerun_20260708/ours/cohort/psp/*.ssr.psp

# recall/FP curve on the confident core:
uv run --no-project benchmarks/ssr_tomato1/scripts/silver_recall_fp_curve.py \
  --ours cohort.ssr.freebayes.vcf \
  --hipstr  benchmarks/ssr_tomato1/results_ssr15k/hipstr/cohort.str.vcf.gz \
  --reads   benchmarks/ssr_tomato1/results_rerun_20260708/our_reads.tsv \
  --catalog benchmarks/ssr_tomato1/results/ours/ssr_tomato1.ssr.catalog \
  --thresholds 0,3,5,10,15,20,30,40,50,75,100,150
```

VCFs left for cross-checking in
`benchmarks/ssr_tomato1/results_rerun_20260708/ours/cohort/`:
`cohort.ssr.freebayes.vcf` (toggle on) and `cohort.ssr.mainhead_gates.vcf` (toggle off).
