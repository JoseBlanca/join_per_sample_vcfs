# The two rough callers — what production does, and what the literature offers instead

**Date:** 2026-07-23
**Scope:** ng step 4, the *parameter pre-pass* ([`ng_proposal.md`](../../ng/spec/ng_proposal.md) §4) —
i.e. the **rough** callers whose only job is to hand the real, cohort-Bayesian caller its
frozen parameters. Both paths are covered: the per-site SNP het caller
([`src/sample_summary/het.rs`](../../../../src/sample_summary/het.rs)) and the STR
confident-genotype gate ([`src/ssr/cohort/prepass.rs`](../../../../src/ssr/cohort/prepass.rs) +
[`rung_ladder.rs`](../../../../src/ssr/cohort/rung_ladder.rs)).
**Predecessor:** [`rough_snp_calling_heuristics_2026-07-07.md`](rough_snp_calling_heuristics_2026-07-07.md)
surveyed the *heuristic-caller* literature (VarScan, LoFreq, bcftools, freebayes gates) and its
recommendations 3.1–3.4 are now largely implemented. This report deliberately looks somewhere
else: the **population-genetics estimation** literature and the **parameter-calibration** code of
GATK/HipSTR, which answer a question the heuristic survey never asked.

---

## 0. The short version

The two rough callers share one design decision, and it is the decision the literature
disagrees with:

> **They call genotypes, threshold them for confidence, and then *count* the survivors.**

The SNP path counts confident hets; the STR path pools reads from confident genotypes. Both
therefore estimate a parameter from a **non-random subset** of the data — the easy sites — and
both inherit the threshold's biases. Every serious estimator in the literature does the same job
by **summing over the unknown genotype** instead of picking one: the genotype is a nuisance
variable to be marginalised, not an intermediate result to be computed and filtered.

That is not a more complicated algorithm. For the SNP path it is *the same three binomial
log-likelihoods we already compute*, added instead of compared; the parameters we want (the error
rate, the het rate) then fall out of a small maximisation. GATK's own STR calibration
(`CalibrateDragstrModel`) is literally this estimator, and it is ~40 lines of arithmetic.

Ranked recommendations for ng, value ÷ cost:

| # | Change | Path | Cost |
|---|---|---|---|
| **R1** | Estimate the het rate and error rate by **maximising the marginal likelihood** over the three genotype states — no threshold, no counting | SNP (and STR, same estimator) | small |
| **R2** | Get the inbreeding coefficient `F` from the **fraction of the genome in runs of homozygosity** (a 2-state HMM over windows), not from `1 − Hobs/Hexp` | SNP | small–medium |
| **R3** | Estimate the STR chemistry with **latent genotypes** — HipSTR's EM (soft) or DRAGstr's marginal ML (no genotypes at all) — instead of the confident-genotype gate | STR | medium |
| **R4** | If a confident-genotype gate is kept, replace the **BIC penalty with an explicit prior odds** on heterozygosity | STR | trivial |
| **R5** | Add the cohort's **excess-heterozygosity + allele-ratio-deviation** signal (HDplot) as a second, orthogonal paralog rejector | both | small |

And one experiment that decides all of it, described in §6: **downsample HG002 and check the
estimate does not move.** A rough caller whose het-rate estimate slides with depth is
mis-estimating `F` in exact proportion to how uneven our cohort's coverage is.

---

## 1. What the two rough callers do today

### 1.1 The SNP path — three binomials, two thresholds, one count

For every covered position that carries a non-reference allele,
[`HetAccumulator::observe_site`](../../../../src/sample_summary/het.rs) scores three models of
the alt count `k` out of depth `n` (`ε̂` is the site's effective error rate, taken from the reads'
own base/mapping qualities and clamped to `[0.02, 0.4]`):

```text
hom-ref :  k ~ Binomial(n, ε̂)        het :  k ~ Binomial(n, ½)        hom-alt :  k ~ Binomial(n, 1−ε̂)
```

then applies a confidence margin `M` twice — once to admit the site as variant at all, once to
split het from hom-alt — and increments one of three counters. A strand-bias veto demotes hets
whose alt allele is strand-confined. Downstream,
[`inbreeding.rs`](../../../../src/paralog/inbreeding.rs) forms
`obs_het = n_het_sites / callable_positions` and `F = clip(1 − obs_het/Hexp, 0, 0.99)`. `F` then
feeds two consumers: the hidden-paralog filter's Wright prior and the cohort caller's
inbreeding-aware genotype prior.

So the whole rough caller is graded on one number, `n_het_sites`, and its two known failure
directions are false hets (paralogs, mismapping) inflating it and threshold losses deflating it.

### 1.2 The STR path — a confident-genotype gate feeding the chemistry

[`resolve_confident_genotype`](../../../../src/ssr/cohort/rung_ladder.rs) runs a BIC
one-allele-vs-two-allele test per (sample, locus) over a cohort burn-in, using the production read
likelihood `Lr`. Genotypes that pass become **CG-seeds**: every one of their reads is hard-labelled
against its nearest parent allele, and those labelled reads are the entire substrate for the
frozen parameters — the per-base error rate `ε`, the per-period stutter shape, the per-sample
stutter level, the `G₀` spread, the sample-group clustering.
[`prepass.rs`](../../../../src/ssr/cohort/prepass.rs) is explicit that this is the hard-label
approximation and that the soft, responsibility-weighted version is the deferred refinement.

### 1.3 The property both share

Neither estimator is wrong about any individual site. The problem is what the *set* of surviving
sites looks like:

- **Threshold ascertainment.** A margin `M` admits a site only when the data are decisive. What is
  decisive depends on depth and quality, so the surviving sites are the deep, clean ones. The het
  rate measured on them is not the genome's het rate, and — critically — the gap between the two
  **depends on the sample's coverage**. In a cohort with 3× coverage spread, a spurious `F`
  gradient across samples is produced by nothing but the threshold.
- **The confident subset teaches the parameters.** The STR pre-pass's own rule — *only confident
  calls may teach the parameters* — is what keeps the bootstrap honest, but it also means `ε` and
  the stutter level are learned from short, pure, well-separated loci and then applied to long,
  interrupted, ambiguous ones. The same bias with a different name.
- **Counting throws away the uncertainty we computed.** An "ambiguous" site is not information-free;
  it is a site with, say, a 0.7 posterior of being het. The accumulator computes that number and
  then discards it in favour of a tally.

---

## 2. The estimator the literature actually uses (and GATK already ships)

The alternative is one line of restructuring: **do not choose a genotype — sum over it.** For a
site with depth `n` and alt count `k`, and with `π_het`, `π_hom_alt` the (unknown) genotype
frequencies:

```text
P(site | ε, π) = π_hom_ref · ε^k (1−ε)^(n−k)
               + π_het     · (½)^n
               + π_hom_alt · (1−ε)^k ε^(n−k)
```

Multiply over sites, maximise over `(ε, π_het)`. The maximiser **is** the het-rate estimate —
no threshold, no counting, no discarded ambiguity — and it comes with the error rate estimated
from the data rather than assumed from base qualities.

**This is not a research proposal; it is GATK's STR calibration verbatim.**
[`DragstrParametersEstimator.java`](../../../../../pop_var_caller/gatk/src/main/java/org/broadinstitute/hellbender/tools/dragstr/DragstrParametersEstimator.java)
(vendored) computes, for each locus case of `(depth n, indel-supporting reads k)`:

```java
log10ProbFunc = log10SumLog10(
    log10PHomRef + k*log10PError + (n-k)*log10PCorrect,
    log10PHet    + n*LOG10_ONE_HALF,
    n == k ? log10PHomVar : NEGATIVE_INFINITY);
```

— the same three states, summed — and grid-searches the error parameter (`GP`) and the variant
prior (`API`) to maximise the total. There is no rough genotyper anywhere in DRAGEN's calibration:
the genotype is integrated out, and the "how often is this a real variant" parameter is estimated
*alongside* the error rate rather than being supplied by a first-pass caller. Three details of
their implementation are worth copying wholesale:

1. **Stratify, then pool the thin strata.** Parameters are fit per (period × repeat-length) cell;
   cells with fewer than 50 loci are merged with their neighbours until they have enough data.
2. **Constrain the fit to be monotone.** Error rate must not fall as the repeat gets longer; a cell
   that violates it is merged back with its neighbour and re-fit. This is how they get stable
   estimates in sparse strata without a smoothing model.
3. **Cheap by construction.** Sites collapse to a `(n, k)` histogram before the fit, so the
   likelihood surface is evaluated over a few thousand cells, not a few hundred million sites.

The same estimator appears independently across the population-genetics literature, which is the
real endorsement:

- **Lynch (2008)** removed the assumption of a known error rate by estimating it *jointly* with
  heterozygosity from read counts; **Maruki & Lynch (2015, 2017)** built the ML genotype/allele
  frequency estimators on top.
- **Bryc, Patterson & Reich (2013)**, *A novel approach to estimating heterozygosity from
  low-coverage genome sequence* — estimates genome-wide heterozygosity "without an intermediate step
  that calls genotypes", jointly learning the **sequencing error distribution and the reference
  bias**. Their reference-bias term is worth stealing on its own (§4.2).
- **ANGSD / realSFS** (Korneliussen 2014) — the standard low-coverage toolkit; individual
  heterozygosity is read off the one-sample site-frequency spectrum computed from genotype
  likelihoods, explicitly "to avoid the bias inherent to hard-calling genotypes".
- **ngsF** (Vieira et al. 2013, *Genome Research*) — per-individual inbreeding coefficients by EM
  over genotype likelihoods, and their paper's point is precisely that ignoring the uncertainty
  biases both genotype calls and allele frequencies.
- **ATLAS** (Kousathanas et al. 2017) and **mlRho** (Haubold 2010) — same family, for ancient and
  low-coverage samples.

**What it costs us.** Almost nothing, and the plumbing already exists: the walker emits a record
for *every* covered position, and [`pileup_to_psp.rs`](../../../../src/pileup/per_sample/pileup_to_psp.rs)
already sees them — it just skips the pure-reference ones before reaching the accumulator. Those
`k = 0` cells are exactly what pins `ε` down. Accumulating a `(depth, alt-count)` histogram instead
of three counters is a smaller object than what we keep today.

---

## 3. `F` deserves a better estimator than `1 − Hobs/Hexp`

Even with a perfect het rate, `F = 1 − obs_het/Hexp` is fragile for our use, for three reasons
that compound: `Hexp` is one cohort-wide scalar applied to every sample regardless of its
ancestry; a uniform artifact floor of false hets biases every sample's `F` downward; and the
estimate has no way to distinguish "slightly outbred everywhere" from "mostly selfed with one
recent outcross".

The modern answer in the animal- and plant-breeding literature is **F_ROH**: the fraction of the
genome lying in runs of homozygosity. For a selfing crop this is not an approximation of `F`, it is
closer to the definition. Two implementations are directly relevant, and both are single-sample
and cheap:

- **ROHan** (Renaud et al. 2019, *Genetics*) — estimates a local heterozygosity rate in windows
  directly from the BAM, then runs a **2-state HMM** (in-ROH / out-of-ROH) genome-wide, reporting
  both the ROH segments and the global het rate *outside* them, with confidence intervals. It works
  down to 4–5× coverage on modern samples.
- **ngsF-HMM** (Vieira et al. 2016) — the same idea from genotype likelihoods, modelling inbreeding
  as an HMM along the genome rather than an independent per-site parameter.

**Why this specifically helps us.** The estimator's output becomes *the fraction of windows in the
low-het state* rather than *the absolute het rate*. A uniform floor of false hets from mismapping
and collapsed paralogs raises the fitted het rate of **both** HMM states and largely cancels out of
the fraction — whereas today it goes straight into `obs_het` and pushes `F` down. It also gives a
free diagnostic for exactly the pathology recorded in the tomato baseline (one strongly
het-inflated sample that no read-quality heuristic could reach): the HMM says whether that
sample's excess heterozygosity is **uniform** (artifact) or **segmental** (a real outcross or an
introgression), which the current single scalar cannot.

Cost: forward–backward over a two-state HMM on ~1 Mb windows of het counts. A few hundred lines,
no new data.

---

## 4. Two smaller SNP-path refinements the literature is clear about

### 4.1 Base qualities are the wrong place to get `ε`

Our `ε̂` is the geometric mean of the reads' own error probabilities. That is a real improvement
over a flat constant (report 2026-07-07 §3.1), but it inherits the instrument's calibration, and
the entire existence of BQSR is evidence that reported qualities are systematically off. Two
options, in increasing cost:

- **Free (recommended):** in the R1 estimator, let `ε` be a *fitted* parameter, and use the
  quality-derived `ε̂` only to **stratify** sites into a handful of quality bins, fitting one `ε`
  per bin. This is BQSR's core idea — empirical error rate by covariate — at 1 % of its machinery,
  and it is exactly DRAGstr's stratification pattern.
- **Not now:** full covariate recalibration. For a non-model organism there is no known-variants
  database, so GATK's BQSR needs bootstrapping; **Lacer** (Chung & Chen 2017) is the published
  alternative that recalibrates without a variant database. Both are large; neither is justified
  for a rough caller unless the stratified fit shows a big residual.

### 4.2 The het model's ½ is not ½

Reads carrying the alternative allele map slightly less often than reference-carrying ones, so a
true het sits near 0.47–0.49, not 0.50. Bryc et al. fit exactly this term. In our formulation it is
one extra free parameter (`k ~ Binomial(n, r)` with `r` fitted once per sample), and it matters
most where our errors are: a systematically low alt fraction makes marginal hets look like
hom-ref, which biases `obs_het` down in a depth-dependent way.

A related and probably larger effect: **our binomials assume independent errors**, but the false
hets we care about (collapsed paralogs, mismapped repeat copies) are systematic and correlated —
that is why 40 identical reads should not count as 40 independent observations (htslib's
`errmod.c` geometric down-weighting is the field's cheap fix, noted in the earlier report §3.5). In
the marginal-likelihood framework the principled version is a **beta-binomial** in place of the
binomial, adding one overdispersion parameter fitted from the data. Worth trying *after* R1, since
R1 makes it a one-parameter change rather than a new model.

### 4.3 A cohort signal we are not using: HDplot

**HDplot** (McKinney et al. 2017) diagnoses collapsed paralogs from two per-site cohort statistics:
the **proportion of samples called heterozygous** (paralogs show excess het — often near 1) and the
**deviation of the allele read ratio in those hets** from 1:1. It is the standard tool for exactly
our problem in non-model organisms, it needs no reference annotation, and it is orthogonal to our
existing coverage-based paralog filter. It is a *cohort* statistic, which fits ng's single-phase
design better than it fits production's per-sample Stage 1 — so ng is the right place to test
whether it catches what the coverage signal misses.

---

## 5. The STR pre-pass — three concrete alternatives to the confident-genotype gate

### 5.1 HipSTR's EM: the same gate, but soft

HipSTR does not gate. [`em_stutter_genotyper.cpp`](../../../../../pop_var_caller/HipSTR/src/em_stutter_genotyper.cpp)
(vendored) runs a plain EM in which the genotypes are latent: the E-step computes, for every read,
a posterior over (allele 1, allele 2, which haplotype it came from) given the current stutter
model; the M-step re-fits the stutter model and the allele frequencies from those
**responsibilities**. Every read contributes to the parameter estimate, weighted by how sure we are
about its origin — an ambiguous locus contributes a little, a clean one contributes a lot, and no
locus is excluded for being hard.

This is the refinement `prepass.rs` already names as deferred. The literature's verdict is that it
is not a refinement but the default: HipSTR's stutter models are locus-specific and learned this
way, from the same reads it later genotypes.

### 5.2 DRAGstr's route: no genotypes at all

§2's estimator applies unchanged to the STR path, with the stratification axis being
(period × tract length) — which is what our stutter level line (`baseline + slope·length`) is
already trying to capture, but fitted from labelled reads rather than from the marginal likelihood.
Two properties of the GATK implementation are worth adopting even if we keep our own likelihood:
**adaptive pooling** of strata with too little data, and the **monotonicity constraint** across
tract lengths, which together replace a parametric line with a shape-constrained fit that cannot
be dragged around by a sparse cell.

**The architectural payoff for ng.** R1 (SNP) and this (STR) are *one component* with a
stratification axis chosen by the marker type. ng's step-4 interface — `SampleSummarizer` +
`CohortEstimator` in [`ng_step_interfaces.md`](../../ng/arch/ng_step_interfaces.md) §4 — would then
have a single implementation serving both paths, which is the "SNP, indel and STR at one level"
property the proposal asks for, arriving for free at exactly the step where production is most
duplicated.

### 5.3 If the gate stays: BIC is the wrong penalty

BIC penalises the second allele by `k·ln n` with `n` = the read count. Two objections, both
standard:

- **Reads are not independent observations.** Stutter makes them correlated, so the effective
  sample size is well below `n` and BIC over-penalises heterozygotes at high depth — the deep,
  informative loci are the ones it is most conservative on.
- **There is a real prior available.** The question "is this sample het at this locus?" has a known
  answer rate: the heterozygosity expected at an STR of this period and length. Replacing
  `k·ln n` with the **log prior odds of het vs hom** is the same arithmetic with a principled,
  locus-specific constant, and STR mutation rates are well characterised as increasing with tract
  length (Sun et al. 2012; Willems et al. 2014; MUTEA, Gymrek et al. 2017). A gate that already
  knows long tracts are more polymorphic will stop treating a long-tract het as an extravagance.

### 5.4 The parameters we do not estimate at all

GangSTR's profile pass sweeps the **insert-size distribution, coverage, GC bins and read length**
before genotyping. We estimate none of these. They only pay off for alleles longer than a read
(spanning/FRR read classes), which our current STR path does not attempt — so this is a note for
scope, not a recommendation: adopting GangSTR-style long-allele recovery *requires* that pre-pass.

---

## 6. How to decide — the experiment that settles it

One measurement discriminates every proposal here, and it is cheap because the data exist
(`benchmarks/ssr_hg002/`, 300× HG002):

> **Downsample one sample to 5×, 10×, 20×, 40×, 80×, 300× and plot the estimated het rate (and `F`)
> against depth. A correct estimator gives a flat line.**

The current counting estimator cannot be flat — the margin `M` admits a depth-dependent subset —
and the size of the slope over 5–40× *is* the size of the spurious `F` gradient across a real
cohort. If the slope turns out to be small at our working depths, R1 is a nice-to-have and this
report's ranking collapses to R2 and R3. If it is large, R1 is the highest-value change available
in either rough caller, and it costs a histogram and a grid search.

Two supporting checks:

- **Synthetic (the only source of a known answer).** Simulate at a known het rate, known `ε`,
  known `F`, known stutter; require each estimator to recover them. This is the one place a
  parameter estimator can be graded directly rather than through downstream call quality — and
  ng's proposal already commits to building the generator (§2).
- **Tomato cohort (the diagnostic that costs nothing).** Correlate each sample's current `F` with
  its mean coverage. A strong correlation is the ascertainment bias made visible on real data, on a
  cohort we already have; a weak one bounds how much R1 can buy.

---

## 7. What not to build

- **A genotype-likelihood cohort EM in the rough pass.** That is the real caller; the rough pass
  must stay robust *without* the parameters it is estimating, or the bootstrap becomes circular.
  Everything above is single-pass or per-sample.
- **Full base-quality recalibration.** Stratified `ε` (§4.1) captures most of it for a fraction of
  the machinery; Lacer/BQSR is a separate tool-sized project.
- **Per-locus HMM refinements or phasing in the rough pass.** Both belong to the real caller.
- **A fourth "artifact" genotype state with a free alt fraction.** Tempting (it would absorb
  collapsed paralogs directly), but it is fitted on the same data as the het state and the two
  compete; the ROH HMM (R2) and HDplot (R5) get at the same artifacts with far less identifiability
  risk.

---

## Sources

**Vendored code read for this report** (paths under `/Users/jose/devel/pop_var_caller/`):
`gatk/src/main/java/org/broadinstitute/hellbender/tools/dragstr/{DragstrParametersEstimator,DragstrHyperParameters}.java`;
`HipSTR/src/em_stutter_genotyper.{h,cpp}`; `htslib/errmod.c`.

**Literature.**
Vieira et al. 2013, *Estimating inbreeding coefficients from NGS data* — [Genome Research 23:1852](https://genome.cshlp.org/content/23/11/1852.full) ([ngsF](https://github.com/fgvieira/ngsF));
Vieira et al. 2016, *Estimating IBD tracts from low coverage NGS data* — [Bioinformatics 32:2096](https://academic.oup.com/bioinformatics/article/32/14/2096/1743296);
Bryc, Patterson & Reich 2013, *A novel approach to estimating heterozygosity from low-coverage genome sequence* — [Genetics 195:553](https://academic.oup.com/genetics/article/195/2/553/5935399);
Kousathanas et al. 2017, *Inferring heterozygosity from ancient and low coverage genomes* — [Genetics 205:317](https://academic.oup.com/genetics/article/205/1/317/6095585);
Maruki & Lynch 2017, *Genotype calling from population-genomic sequencing data* — [G3 7:1393](https://academic.oup.com/g3journal/article/7/5/1393/6028285);
Renaud et al. 2019, *Joint estimates of heterozygosity and runs of homozygosity* — [Genetics 212:587](https://pubmed.ncbi.nlm.nih.gov/31088861/) ([ROHan](http://grenaud.github.io/ROHan/));
ANGSD heterozygosity / realSFS — [popgen.dk](https://www.popgen.dk/angsd/index.php/Heterozygosity);
McKinney et al. 2017, HDplot — [implementation](https://github.com/edgardomortiz/paralog-finder);
Chung & Chen 2017, *Lacer* — [bioRxiv 130732](https://www.biorxiv.org/content/10.1101/130732v2.full);
GATK [CalibrateDragstrModel](https://gatk.broadinstitute.org/hc/en-us/articles/27008029541787-CalibrateDragstrModel);
[HipSTR](https://hipstr-tool.github.io/HipSTR/); [GangSTR](https://github.com/gymreklab/GangSTR).
