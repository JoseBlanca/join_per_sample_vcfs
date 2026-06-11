# Why our FP QUAL inflates with depth — diagnosis & fix options

**Date:** 2026-06-10 · **Branch:** `qual-analysis` · Follow-up to
[qual_distribution_by_depth_2026-06-10.md](qual_distribution_by_depth_2026-06-10.md).

> **Outcome (2026-06-11):** the `combo` fix (allele balance + strand/position
> bias) is now the **production variant QUAL** — `src/vcf/qual_refine.rs`,
> applied unconditionally at VCF-encode time (QUAL column only; genotypes,
> GQ and AF unchanged). The other prototypes (overdisp, fbformula,
> balance-only, bias-only) and the `PVC_QUAL_EXPERIMENT` env scaffolding
> were removed. The balance test compares against the *called-genotype*
> expected alt fraction (0.5 for a single het, cohort AF otherwise), so it
> is cohort-safe and skips hom-alt (`p_exp ≥ 0.9`) so a hom-alt's few
> error reads aren't mistaken for a deficit. Validated to reproduce the
> approved single-sample behaviour (SNP FP median → 0 at every depth, TP
> preserved: 76 / 555 / 5781 at 5×/30×/full).

## Observation

In the per-tool box plots (SNPs, y fixed to 0–100), **ours** shows the
median **false-positive** QUAL *climbing with coverage* — 5x→full the FP
QUAL median goes 1 → 3 → 150 and the box balloons past 100. freebayes
does **not**: its FP QUAL stays pinned near 0 at every depth while only
its TP QUAL rises. A caller whose *false* positives look *more* confident
as you sequence deeper is mis-calibrated — the opposite of what we want.

## Mechanism (confirmed in code + data)

QUAL is `−10·log10 P(K=0 | data)` (no non-ref allele), computed in
`compute_qual_via_exact_af` from per-sample **genotype likelihoods**
(`posterior_engine.rs`). Those likelihoods come from
`standard_log_likelihood` (`per_group_merger.rs:1948`):

```
logL(G) = Σ_{a ∉ G} q_sum[a]            # every read of an allele not in G is an "error"
        + multinomial(reads in G; p = copies/ploidy)
```

`q_sum[a]` is the sum over allele-`a` reads of `ln(P_err)`, and per read
(`open_record.rs:779`) `ln(P_err) = max(bq_ln_perr, mq_log_err)` — the
worse of base- and mapping-quality error.

So the het-vs-hom-ref log-ratio is dominated by `N_alt · (−ln ε)`:
**linear in the number of alt-supporting reads**, with a *tiny* fixed
ε. This is a pure **i.i.d. binomial/multinomial error model** — each alt
read is an independent rare error, so evidence against hom-ref piles up
linearly with depth.

### Empirical confirmation (ours, SNPs, TP=0003 vs FP=0001 vs truth)

| | depth | n | QUAL med | DP med | alt-VAF med | N_alt med |
|--|------:|--:|--------:|------:|-----------:|---------:|
| TP | 5    | 5156 | 86 | 5 | 1.00 | 4 |
| TP | 30   | 6730 | 571 | 28 | 0.57 | 17 |
| TP | full | 6687 | 5862 | 286 | 0.515 | 164 |
| FP | 5    | 435 | 1 | 6 | 0.43 | 2 |
| FP | 30   | 1351 | 3 | 23 | **0.235** | 5 |
| FP | full | 1403 | 150 | 227 | **0.196** | 45 |

The FP sites carry a **persistent ~20 % alt fraction** (VAF p25/p50/p75
at full = 0.17/0.20/0.25) that does *not* decay toward the error rate as
depth grows — so `N_alt` rises linearly with coverage (2→5→45) and QUAL
rises with it. A true diploid genotype is VAF 0.5 (het) or 1.0
(hom-alt); **0.20 is consistent with no real genotype**, yet the model
calls it het because the only question it asks is "is the alt fraction
bigger than the sequencing error rate ε?" — and at ε≈10⁻³ the answer is a
resounding, depth-amplified yes.

### What it is *not*

- **Not** a mapping-quality problem. `MQDiff` (MQ_ref − MQ_alt) median is
  **0.0** for both TP and FP; only ~2 % of FPs have |MQDiff|>3. The FP
  alt reads carry the same MAPQ≈70 and high base quality as TP reads, so
  `q_sum`'s existing `max(bq, mq)` error term gives ε≈10⁻³ either way.
  MAPQ-aware tweaks won't move these.
- **Not** caught by the existing gates: `--min-alt-obs-per-sample` is a
  raw count (a 20 %-VAF site easily clears it at depth) and
  `--min-mapq-diff-t` keys on a MAPQ difference that isn't there.

These are systematic artifacts (paralogs that map uniquely, recurrent
context errors, or reference-assembly issues) that produce a stable
minority allele of high-quality, well-mapped reads — exactly the case an
i.i.d. binomial cannot down-weight.

## Fix options (ranked)

1. **Overdispersed / explicit systematic-error term (principled).** Stop
   treating the per-site error rate as a fixed ε. Either (a) a
   **beta-binomial** site error rate (Dirichlet/Beta nuisance on the
   error allele fraction), or (b) an explicit **error-allele-fraction**
   parameter the hom-ref hypothesis is allowed to spend on a recurrent
   minority allele. Both make the alt-read evidence **saturate** instead
   of growing linearly, so a sustained ~20 % VAF is explained as artifact
   rather than forcing a het — capping the depth-driven QUAL inflation.
   This is the standard, calibration-correct remedy.

2. **Allele-balance-aware QUAL (surgical).** A germline diploid het is
   VAF≈0.5; penalise QUAL by a binomial allele-balance test of the best
   genotype's observed VAF against its expected value (GATK's old `AB`,
   freebayes' balance behaviour). Directly deflates the VAF~0.2 FPs.
   Caveat: would also penalise genuine low-VAF (mosaic/somatic) calls —
   acceptable for a germline diploid caller, but make it
   ploidy/mode-aware.

3. **Strand / read-position bias into QUAL.** We already collect
   per-allele `fwd` (strand) and `placed_left`/`placed_start` (position)
   bias counts in the merger, but they currently feed nothing in QUAL.
   Systematic artifacts are often strand- or position-biased; a bias
   penalty would target them. *Not yet measured here* — worth checking
   whether these FPs are biased (MAPQ wasn't, allele balance was).

4. **Depth-robust / capped QUAL (stopgap).** Cap per-read evidence or
   report a depth-normalised confidence. Stops the runaway but doesn't
   fix the underlying calibration; use only as a bridge.

**Recommendation:** prototype (1) — an error-allele-fraction / beta-binomial
term in the hom-ref likelihood — as the root-cause fix, with (2) as a
cheaper interim lever to validate that allele balance is indeed the
discriminator (the data says it is). Re-run the depth sweep after each to
confirm the FP-QUAL-vs-depth curve flattens while TP QUAL is preserved.

## Prototype results (2026-06-10)

All three options were prototyped as **reductions on the engine QUAL**
in `src/vcf/qual_experiment.rs` (gated by `PVC_QUAL_EXPERIMENT`; default
unset → byte-identical), then `var-calling` was re-run on each depth's
existing `.psp` per mode and re-scored against truth
(`benchmarks/lib/run_qual_experiments.sh`). Pileup/GATK/freebayes are
untouched, so the sweep is seconds per depth. Results overlay in the
"QUAL re-calibration experiments" section of `qual_depth_dashboard.py`.

**SNP median QUAL by mode across depth (5x → 301x):**

| | 5 | 10 | 20 | 30 | 50 | 100 | 301 |
|--|--:|--:|--:|--:|--:|--:|--:|
| **FP** engine | 0.9 | 0.9 | 1.5 | 3.5 | 12.5 | 44.9 | **150** |
| **FP** balance | 0.4 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | **0.0** |
| **FP** bias | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | **0.0** |
| **FP** overdisp | 0.3 | 0.2 | 0.5 | 1.8 | 9.6 | 41.8 | **134** |
| **TP** engine | 86 | 177 | 372 | 571 | 965 | 1937 | **5862** |
| **TP** balance | 85 | 177 | 372 | 570 | 964 | 1935 | **5862** |
| **TP** bias | 76 | 162 | 358 | 556 | 943 | 1905 | **5785** |
| **TP** overdisp | 15 | 40 | 89 | 132 | 209 | 360 | **813** |

**Verdict (SNPs):**

- **Balance (opt 2) — the winner.** Flattens FP QUAL to ~0 at *every*
  depth while leaving TP QUAL essentially unchanged (5862 = 5862 at full;
  ≤1 % anywhere). Eliminates the inflation with no TP cost.
- **Bias (opt 3) — also excellent**, and works even at 5–10x where the
  balance binomial has little power; tiny TP cost (~1–3 %). That it works
  at all confirms these FPs are *also* strand/position-biased — a
  complementary signal to allele balance.
- **Overdisp (opt 1) — fails as prototyped.** A unimodal Beta-Binomial
  barely dents the inflation (134 vs 150 at full) yet craters TP (813 vs
  5862, ~7×). It cannot tolerate a 20 % artifact without also tolerating a
  50 % het — the predicted tuning tension, confirmed across the param
  sweep (`conc` 12/60/150 all trade TP for little FP gain).

### Does freebayes' *formula* explain its flat FP QUAL? (No.)

We implemented freebayes' actual scoring (`PVC_QUAL_EXPERIMENT=fbformula`).
Reading the vendored source (`freebayes/src/DataLikelihood.cpp`), its
standard genotype likelihood is **identical to ours** — same error term
(`prodQout = Σ max(ln_bq, ln_mq)` over out-of-genotype reads = our
`q_sum`) and the same `multinomialSamplingProbLn` — except for one
ingredient: the **Read Dependence Factor** (`-D`, default 0.9), which
down-weights successive error reads as correlated:
`prodQout *= (1 + (countOut−1)·RDF)/countOut`. We applied exactly that.

| SNP FP median QUAL | 5 | 30 | 100 | 301 | FP count @301 |
|--|--:|--:|--:|--:|--:|
| engine | 0.9 | 3.5 | 44.9 | **150** | 1403 |
| fbformula (RDF 0.9) | 0.0 | 0.0 | 31.7 | **114** | 1403 |
| **freebayes (actual)** | 0.0 | 0.0 | 0.0 | **0.0** | **8650** |
| balance | 0.4 | 0.0 | 0.0 | **0.0** | 1403 |

freebayes' formula on *our* calls (`fbformula`) only nudges the inflation
(114 vs 150 at full; RDF=0.9 is too weak) — it does **not** reproduce
freebayes' flat-zero FP QUAL. The reason is the FP **population**, not
the formula: freebayes emits a *different and ~6× larger* FP set (8650 vs
our 1403 at full), dominated by low-count noise it reports at QUAL≈0, so
its median sits at 0. Our systematic ~20 %-VAF artifacts are a small
slice of freebayes' FPs and don't move its median. **Takeaway:**
freebayes' advantage on these sites is *what it calls*, not *how it
scores* — and our targeted `balance`/`bias` penalties control our own
artifact FPs better than freebayes' RDF does, at lower TP cost
(fbformula also shaves ~10 % off TP; balance shaves ~0 %).

**Indels** are a separate, harder story: balance/bias clear the
mid-depth FPs but at 100–301x indel FPs stay high under *every* mode
(balance FP 752/1003 vs engine 786/1021) — those FPs are not low-VAF or
biased the same way (likely alignment/repeat ambiguity at VAF≈0.5), so
allele balance can't see them. overdisp again craters TP.

**Recommendation:** ship **allele balance** as the SNP QUAL fix, with
**strand/position bias** as a complementary low-depth lever (the two can
combine). Promote it from the encode-time prototype into the QUAL
finaliser, make it ploidy/mode-aware (don't penalise legitimately
low-VAF calls in somatic/mosaic modes), and re-check **recall** (the
medians say TP is preserved, but confirm the TP lower tail isn't gated
out). Drop the unimodal-Beta-Binomial overdisp; a *mixture* error model
(a dedicated "artifact" component) would be the principled version if
option 1 is revisited.

## Reproduce

The fix now lives in `src/vcf/qual_refine.rs` (always-on; no flags). The
per-tool box plots and depth views are in `qual_depth_dashboard.py`.

The mode-comparison sweep that selected `combo` (a one-off
`run_qual_experiments.sh` that re-scored each candidate formula behind a
`PVC_QUAL_EXPERIMENT` env switch) was **retired once `combo` shipped** —
the numbers in this report are its record. The TP/FP field
characterisation (VAF, MQ, N_alt) was an ad-hoc `bcftools isec` of each
depth's ours VCF against the truth SNPs, querying
`%QUAL %INFO/DP %INFO/AF %INFO/MQDiff [%AD]`.
