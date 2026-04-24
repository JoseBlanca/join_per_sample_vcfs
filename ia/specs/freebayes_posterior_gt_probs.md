# freebayes: from BAM reads to genotype posterior probabilities

This document walks, step by step, through how freebayes turns aligned reads
in a BAM file into per-sample genotype posterior probabilities. It is written
for geneticists, not statisticians: every formula is accompanied by an
intuitive explanation, and every step is anchored to a specific location in
the freebayes C++ source code (under `freebayes/src/`).

freebayes is described in Garrison & Marth, *Haplotype-based variant
detection from short-read sequencing*, arXiv:1207.3907 (2012). The
implementation here is the reference implementation maintained at
<https://github.com/freebayes/freebayes>.

## Call flow

```
main()                                     freebayes.cpp:57
  AlleleParser::getNextAlleles()           AlleleParser.cpp    <- reads BAM, builds observations
  AlleleParser::buildHaplotypeAlleles()    AlleleParser.cpp    <- builds haplotype-level alleles
  getGenotypesByPloidy()                   Genotype.cpp        <- enumerates possible genotypes
  calculateSampleDataLikelihoods()         DataLikelihood.cpp  <- P(reads | genotype), per sample
    probObservedAllelesGivenGenotype()     DataLikelihood.cpp
  convergentGenotypeComboSearch()          Genotype.cpp        <- search over joint genotypes
    GenotypeCombo::calculatePosteriorProbability()  Genotype.cpp  <- likelihood × prior
      probabilityGivenAlleleFrequencyln()  Genotype.cpp        <- prior piece 1: P(G|AF)
      alleleFrequencyProbabilityln()       Ewens.cpp           <- prior piece 2: P(AF), Ewens
      hweProbGenotypeFrequencyln()         Genotype.cpp        <- prior piece 3 (optional): HWE
      binomialProbln / multinomialSampling Genotype.cpp        <- prior piece 4 (optional): obs bias
  marginalGenotypeLikelihoods()            Marginals.cpp       <- sum combos to per-sample posteriors
  ResultData::vcf()                        ResultData.cpp      <- emit GT, GQ, GL, QUAL
```

## Key source files

All under `freebayes/src/`:

| File | Role |
|---|---|
| `freebayes.cpp` | Main loop: one iteration per candidate site |
| `AlleleParser.cpp` | BAM parsing, per-read → per-allele observations, haplotype construction |
| `Allele.cpp`, `Sample.cpp` | Allele/observation data structures, per-sample read accumulators |
| `DataLikelihood.cpp` | `P(reads | genotype)` per sample |
| `Genotype.cpp` | Genotype enumeration, joint-combo search, prior and posterior assembly |
| `Ewens.cpp` | Ewens' sampling formula prior on allele frequencies |
| `Marginals.cpp` | Sum over combos to obtain per-sample marginal posteriors |
| `ResultData.cpp` | VCF record construction (GT, GQ, GL, QUAL, depth fields) |

---

## The mathematical model

freebayes uses Bayes' theorem, **not** an EM iteration. For each candidate site
it enumerates a joint assignment of genotypes across all samples (a "genotype
combo"), computes likelihood × prior for that combo, then marginalises back to
per-sample posteriors by summing over all combos. The marginalisation step
converts the joint distribution (one probability per full cohort assignment)
into the per-sample distributions that a VCF actually needs: to get the
probability that sample `S` has genotype `G`, you add up the posteriors of
every combo that assigned `G` to `S`, regardless of what the other samples
were assigned.

```
P(genotype combo | reads) ∝ P(reads | genotype combo) × P(genotype combo)
```

All arithmetic is done in natural log space (`long double` ln-probabilities) to
avoid underflow. The master equation is assembled in
[Genotype.cpp:1499-1605](../../freebayes/src/Genotype.cpp#L1499-L1605)
(`GenotypeCombo::calculatePosteriorProbability`):

```cpp
priorProb     = priorProbG_Af + priorProbAf
              + priorProbObservations + priorProbGenotypesGivenHWE;
posteriorProb = priorProb + probObsGivenGenotypes;
```

The rest of this document explains each of those terms.

---

## Three ways to solve the same problem: freebayes, GATK, MCMC

Before diving into the mechanics, it is worth placing freebayes' approach in
context. The core statistical goal is shared by every modern joint variant
caller:

> For each sample `s`, compute `P(Genotype_sample | all reads in the cohort)`.

The difficulty is a chicken-and-egg dependency:

- To compute a per-sample posterior you need a **prior**, and the correct
  prior at a locus depends on the **allele frequency** in the cohort.
- To know the allele frequency you need everyone's **genotypes**.
- But the genotypes are exactly what you are trying to compute.

Three standard strategies break this circle. freebayes uses one; GATK
GenotypeGVCFs uses another; a third exists in principle but is rarely used
in production pipelines.

### 1. Joint-combo marginalisation (freebayes)

Enumerate (or search over) joint assignments of genotypes to all samples
simultaneously. Under any fixed joint assignment, the allele frequency is
just whatever that assignment implies, so the prior and the likelihood are
both well-defined and the per-combo posterior can be computed directly.

To recover per-sample answers, sum the combo posteriors into each sample's
bucket — the marginalisation step described earlier. This is a numerical
approximation to the full Bayesian integration

```
P(G_s | reads) = Σ_{all other-sample assignments} P(whole combo | reads)
```

where the sum runs implicitly over all possible allele frequencies as well,
because each combo implies an allele frequency. freebayes approximates the
sum by the high-posterior combos found during its convergent search. This
is the closest of the three approaches to a full Bayesian answer, and it is
what the rest of this document describes in detail.

### 2. Empirical-Bayes EM (GATK GenotypeGVCFs)

Instead of summing over joint assignments, treat the allele frequency `p` as
a **latent parameter** and estimate a single best value `p̂` by
Expectation-Maximisation:

- **E step.** Given the current `p`, compute soft per-sample genotype
  posteriors for every sample (easy — for each sample this is a
  one-variable Bayes calculation).
- **M step.** Using those soft posteriors as weights, recount the expected
  allele copies across the cohort and update `p`.
- Iterate until `p` stops moving (usually 3–5 rounds).

After convergence, report the E-step per-sample posteriors evaluated at
`p̂`. This is *empirical Bayes* — the prior (driven by `p̂`) has been tuned
to the data, and the resulting per-sample posteriors are conditional on
that single tuned value.

GATK uses this approach because it can compute everything from the PL
fields already stored in each sample's gVCF, without re-reading the BAMs.
See [posterior_gt_probs.md](posterior_gt_probs.md) for the full GATK
walk-through.

### 3. Full Bayesian sampling (MCMC)

The textbook "gold standard" is to draw Monte Carlo samples from the joint
posterior `P(all genotypes, p | reads)` — typically by Gibbs sampling or
Hamiltonian Monte Carlo. Each draw is a complete `(genotypes, p)` tuple; by
tallying many draws you recover the per-sample marginal posteriors and
also a distribution over `p`, with correctly propagated uncertainty.

No mainstream production caller does this — it is much slower than the
other two approaches — but it is the reference point both freebayes and
GATK are ultimately approximating.

### How the three compare

| Approach | What is integrated over | How uncertainty in allele frequency is handled |
|---|---|---|
| Marginalisation (freebayes) | Joint genotype combos | Implicitly propagated — each combo implies its own AF, combos are weighted by posterior |
| EM (GATK) | Nothing; `p` is a point estimate | Discarded — per-sample posteriors conditioned on `p̂` as if it were known |
| MCMC | Everything, including `p` | Fully propagated by sampling |

**Practical consequence.** In large cohorts with well-determined allele
frequencies, all three agree closely because `p̂` is tight and the
integration collapses to plugging it in. In small cohorts, at rare
variants, or at sites near the detection threshold, the three can
disagree — EM tends to be slightly over-confident (it ignores uncertainty
in `p`), marginalisation is closer to the true Bayesian posterior, and
MCMC is the reference. This is why freebayes and GATK can produce
different per-sample calls on the same BAM data for borderline variants:
they are genuinely computing slightly different things, not just
differently-implemented versions of the same thing.

---

## Stage 1: From BAM to per-sample allele observations

Before any probability is computed, freebayes has to turn the reads in the BAM
file into a structured view of "at this reference position, sample S has N
reads supporting each of these alleles, at these base qualities."

### 1a. Streaming the BAM, one reference position at a time

Location: `AlleleParser::getNextAlleles()` in [AlleleParser.cpp](../../freebayes/src/AlleleParser.cpp),
driven from the main loop at [freebayes.cpp:117](../../freebayes/src/freebayes.cpp#L117).

freebayes walks the reference genome position by position. At each position it
maintains a **sliding window** of overlapping reads. For each read it:

1. Filters the read by mapping quality (`-m`, default 1), CIGAR validity,
   duplicate flags, and minimum alignment score.
2. Traverses the CIGAR string to decide what "allele observation" the read
   contributes at the current reference position — one of:
   - **reference match** (single base, matches reference)
   - **SNP** (single base, mismatches reference)
   - **MNP** (multi-nucleotide substitution)
   - **insertion** or **deletion** (indel, carried as a single
     variable-length allele)
   - **complex** (mixture of substitutions and indels)
3. Stores each observation with its base quality `Q` (Phred), its mapping
   quality `MQ`, its strand, and its position within the read.

Each `Allele` object thus records not only the base sequence but the quality
metadata that will drive the likelihood calculation later.

### 1b. Candidate allele discovery

At each site, freebayes looks at every allele that has been observed in at
least one sample, and keeps those that pass the minimum-count / minimum-fraction
filters ([freebayes.cpp:227](../../freebayes/src/freebayes.cpp#L227),
`sufficientAlternateObservations()`, controlled by `--min-alternate-count` and
`--min-alternate-fraction`, defaults 2 and 0.05). Everything else is treated as
noise.

**Genetics viewpoint:** freebayes is not evaluating every theoretical nucleotide
change; it only evaluates alleles the reads actually support. That's why
freebayes VCFs often have multiple ALTs at a single site — these are the
alternate alleles observed above threshold, not a pre-enumerated list.

### 1c. Haplotype-based allele definition

Location: `AlleleParser::buildHaplotypeAlleles()` at [freebayes.cpp:284](../../freebayes/src/freebayes.cpp#L284).

This is the step that makes freebayes "haplotype-based":

- If multiple candidate alleles overlap (e.g. a SNP 3 bp away from an
  insertion), freebayes extends the window to cover all of them and rewrites
  the alleles as **haplotype alleles**: each one is the full variant sequence
  across the extended window.
- A 3 bp insertion plus a nearby SNP becomes a single alt allele, not two
  alts.
- Every read is re-assigned to the haplotype allele it is most consistent
  with, with partial support tracked for reads that only cover part of the
  window.

**Genetics viewpoint:** this is why freebayes handles complex variants and
indels cleanly. It refuses to treat a 3 bp insertion + SNP as two independent
events. Instead, the genotype at this site is phrased in terms of
*haplotypes*: "sample S has one copy of the reference haplotype and one copy
of the ACGT-insertion-plus-SNP haplotype."

### 1d. Candidate genotypes per sample

Location: `getGenotypesByPloidy()` at [freebayes.cpp:331](../../freebayes/src/freebayes.cpp#L331).

Given the candidate allele set and each sample's ploidy (default 2, diploid;
controllable via `--ploidy` or `--cnv-map` per-sample), freebayes enumerates
*all possible unordered multisets of ploidy alleles*. For diploid with alleles
`{ref, alt1, alt2}` that gives six genotypes:

```
ref/ref, ref/alt1, ref/alt2, alt1/alt1, alt1/alt2, alt2/alt2
```

Each sample gets the same candidate-genotype set at this site. They will
differ in the likelihoods assigned to each.

---

## Stage 2: The likelihood — P(reads | genotype), one sample at a time

Location: `probObservedAllelesGivenGenotype()` in
[DataLikelihood.cpp:5-166](../../freebayes/src/DataLikelihood.cpp#L5-L166).

For one sample and one hypothesised genotype `G`, we want the probability of
the reads we observed. Reads are treated as independent, so this is a product
(sum of logs) over reads. For each read freebayes asks two questions:

### 2a. "Is this read consistent with genotype G?"

Every read supports some allele (or partially supports several). If that
allele is **in** the genotype `G`, the read is consistent: the read was
generated from one of the haplotypes in `G`. If the allele is **not** in `G`,
the read must be a sequencing error.

The per-base probability of error comes from the base quality. If `Q` is the
Phred score, then `error_prob = 10^(-Q/10)`. Mapping quality `MQ` is combined
in the same way; freebayes uses the **worst** of the two (max in log space,
[DataLikelihood.cpp:34](../../freebayes/src/DataLikelihood.cpp#L34)):

```
P(read is wrong) = max( 10^(-Q/10), 10^(-MQ/10) )
P(read is right) = 1 − P(read is wrong)
```

- Reads **inconsistent** with `G` contribute `log(1 − qual)` — the probability
  that they are in fact errors. Their contribution accumulates into
  `prodQout`.
- Reads **consistent** with `G` contribute an allele sampling probability
  (next subsection).

### 2b. "If the read is consistent, which haplotype did it come from?"

For a heterozygote `A/B`, half the reads should come from the A haplotype
and half from the B haplotype. freebayes encodes this with the **allele
sampling probability** `alleleSamplingProb(obs)` — the fraction of the
genotype's copies that carry the observed allele. For diploid:

```
hom ref     A/A    asampl(A) = 1    asampl(B) = 0
het         A/B    asampl(A) = 0.5  asampl(B) = 0.5
hom alt     B/B    asampl(A) = 0    asampl(B) = 1
```

This sampling probability is multiplied per read (summed in log space). A
read supporting `B` at an `A/A` site contributes `asampl = 0`, which in log
space is `-∞` — a zero-likelihood disaster — so the code instead routes such
reads through the "read is wrong" branch above, which yields a finite
penalty proportional to the base/mapping quality.

**Allele balance / observation counts (standard GL mode):**
[DataLikelihood.cpp:155](../../freebayes/src/DataLikelihood.cpp#L155) adds
one more term:

```
multinomialSamplingProbLn(alleleProbs, observationCounts)
```

This is the multinomial probability of seeing exactly these counts of each
allele given the genotype's sampling probabilities. Concretely: under a
clean heterozygote we expect a 50/50 count split, and a 10/0 split becomes
less likely the larger the depth. This term disciplines the likelihood
against genotypes that are theoretically possible but produce implausible
allele balances.

### 2c. Read-dependence correction (RDF)

[DataLikelihood.cpp:148](../../freebayes/src/DataLikelihood.cpp#L148):

```cpp
prodQout *= (1 + (countOut - 1) * RDF) / countOut;
```

If 20 reads all disagree with the genotype, they are probably not 20
independent errors — they more likely represent one PCR duplicate family or
a systematic mapping issue. freebayes therefore **down-weights the joint
evidence** from multiple contradicting reads. `RDF` defaults to 0.9, so the
10th contradicting read contributes much less evidence than the 1st.

**Genetics viewpoint:** this is freebayes' hedge against correlated errors.
Without it, a single duplicate cluster of 30 reads with a common sequencing
artefact would utterly dominate the likelihood.

### 2d. Contamination and reference bias

[DataLikelihood.cpp:116-130](../../freebayes/src/DataLikelihood.cpp#L116-L130).
If contamination estimates are provided (`--contamination-estimates`), the
sampling probability for the reference allele is rescaled by
`probRefGivenHet / 0.5`. Without a contamination file this is a no-op —
default behaviour is as above.

### 2e. The resulting per-sample data likelihood vector

After processing all reads, for each sample freebayes has, for *every*
candidate genotype, a log-likelihood value:

```
L_s(G) = log P(reads_s | G)
```

These are stored as `SampleDataLikelihood` objects sorted best-first
([DataLikelihood.cpp:270](../../freebayes/src/DataLikelihood.cpp#L270)).

**These values become the GL field** in the VCF (after conversion to log₁₀
and max-normalisation to 0; see Stage 6 below). At this stage the prior has
not yet been applied.

---

## Stage 3: Joint genotyping across samples

The per-sample likelihoods alone are not enough to call genotypes. Consider
a site where one sample has weak heterozygous evidence: we cannot tell if
it is a real het or sequencing noise. The answer depends on what the
*other* samples say. If many samples carry the alt allele the site is
plausibly a real common variant; if nobody else does, it is probably error.

freebayes handles this by computing the likelihood × prior for a **joint
assignment of genotypes across all samples**, called a `GenotypeCombo`
([Genotype.h:166](../../freebayes/src/Genotype.h#L166)).

### 3a. What a combo is

A genotype combo is a complete choice of one genotype per sample at this
site. For two samples with three candidate genotypes each, there are 3×3 = 9
combos; for 100 samples there are 3¹⁰⁰ — astronomically many, which is why
freebayes does not enumerate them all.

A combo carries:

- the chosen genotype for each sample
- the joint data likelihood `Σ_s L_s(G_s)`
- **cached allele frequencies across all samples** (`alleleCounters`)
- **cached observation totals, strand bias, position bias** for each allele
  pooled across all samples

Those cached counts are what make the population-wide prior computable.

### 3b. The joint data likelihood

[DataLikelihood.cpp:250-256](../../freebayes/src/DataLikelihood.cpp#L250-L256).
Samples are independent given their genotypes, so:

```
log P(reads | combo) = Σ_s L_s(G_s)
```

This is just addition (in log space) of each sample's likelihood for its
assigned genotype.

### 3c. The combo-level prior

The prior has up to four additive components, all in log space. The assembly
happens in `calculatePosteriorProbability()`
([Genotype.cpp:1499-1605](../../freebayes/src/Genotype.cpp#L1499-L1605)).

#### Component 1 — P(genotypes | allele frequency): the multinomial coefficient

[Genotype.cpp:1391-1405](../../freebayes/src/Genotype.cpp#L1391-L1405),
`probabilityGivenAlleleFrequencyln()`.

Given that the population contains some specific allele counts (derived from
the combo), how many ways can those allele copies be *distributed into
genotypes*? A heterozygote can be written `A/B` or `B/A` (two arrangements);
a homozygote only one way. The multinomial coefficient counts these
arrangements:

```
log P(G|AF) = log(permutations) − log multinomial_coefficient(n_total_alleles, allele_counts)
```

**Geneticist intuition:** at a given allele frequency, more heterozygotes
means more ways to arrange the alleles, so combos rich in hets get a slight
*a priori* boost over combos that concentrate the same allele counts into
homozygotes. This is the term that encodes the combinatorial aspect of
Hardy–Weinberg before any population-genetics law is invoked.

The `diffusionPriorScalar` parameter (default 1.0) divides this term if set,
dampening its effect when you have very many samples and a strong belief
that the site is under selection.

#### Component 2 — P(allele frequency): Ewens' sampling formula

[Ewens.cpp](../../freebayes/src/Ewens.cpp), called at
[Genotype.cpp:1582](../../freebayes/src/Genotype.cpp#L1582) when
`--no-population-priors` is *not* set.

Given a population mutation rate `θ` (theta, controlled by `-T`, default
0.001, optionally scaled by haplotype length at
[freebayes.cpp:309](../../freebayes/src/freebayes.cpp#L309)), Ewens' formula
gives the probability that a sample of `M` chromosomes has the observed
allele-frequency *partition*. The implementation
([Ewens.cpp:32-51](../../freebayes/src/Ewens.cpp#L32-L51)) computes:

```
             M!       θ^k₁       θ^k₂
P(partition) = ──── · ─────── · ─────── · ...
             θ · Π(θ+h)  1^k₁ · k₁!   2^k₂ · k₂!
```

where `k_f` is the number of alleles with frequency `f` in the sample, and
`M` is the total number of chromosomes.

**Geneticist intuition:** Ewens' formula is the neutral-theory expectation
for how common vs. rare alleles should be distributed in a finite sample of
chromosomes drawn from an infinite population. With `θ = 0.001` it strongly
prefers allele-frequency partitions dominated by one major allele with a
few rare ones — exactly the pattern real SNPs tend to show. A 50/50 split
is penalised much more than a 98/2 split. This is the mechanism by which
freebayes bakes "rare variants are a priori more plausible than common
ones" into its calls.

#### Component 3 — P(genotypes | HWE) (optional)

[Genotype.cpp:1417-1493](../../freebayes/src/Genotype.cpp#L1417-L1493),
enabled by default (`hwePriors = true`).

For each sample's genotype, freebayes computes the probability of observing
that genotype under Hardy–Weinberg equilibrium given the combo-wide allele
frequencies:

```
log P(G_s | HWE) = log multinomial(ploidy, allele_counts_in_G_s)
                 + log multinomial(total_genotypes, genotype_counts_across_samples)
                 − log multinomial(total_alleles, allele_counts_across_samples)
```

These per-sample contributions sum to the combo-level HWE prior. Combos
whose per-sample genotype counts (n_AA, n_AB, n_BB) match HWE expectations
get a boost; combos showing excess heterozygosity or homozygosity are
penalised.

#### Component 4 — Observation-level bias priors (optional)

[Genotype.cpp:1531-1569](../../freebayes/src/Genotype.cpp#L1531-L1569).

If a candidate alt allele is supported only by forward-strand reads, or only
by reads near one end of the fragment, it is much more likely to be a
systematic artefact than a real variant. freebayes checks:

- `binomialObsPriors` (default on): binomial `P(forward_strand_count | total_obs, 0.5)`,
  and similarly for read placement left/right and fragment start/end.
- `alleleBalancePriors` (default on): multinomial fit of observed allele
  counts against expected genotype sampling probabilities.

These are added to `priorProbObservations`. A combo that implies the alt
allele exists, but whose observations exhibit gross strand or positional
bias, is down-weighted.

### 3d. The combo posterior

Putting it together ([Genotype.cpp:1595-1596](../../freebayes/src/Genotype.cpp#L1595-L1596)):

```cpp
priorProb     = priorProbG_Af + priorProbAf + priorProbObservations + priorProbGenotypesGivenHWE;
posteriorProb = priorProb + probObsGivenGenotypes;   // likelihood × prior in log space
```

`posteriorProb` is the unnormalised log-posterior of this specific combo.

---

## Stage 4: Searching the combo space

There are too many combos to enumerate. freebayes therefore uses
`convergentGenotypeComboSearch()`
([Genotype.cpp:1101-1267](../../freebayes/src/Genotype.cpp#L1101-L1267)),
invoked from [freebayes.cpp:488](../../freebayes/src/freebayes.cpp#L488).

### 4a. Starting point: the GL-max combo

The initial combo picks each sample's *best* genotype by data likelihood
alone — that is, each sample's `rank 0` genotype. This is the
"genotype-likelihood maximum" combo and it is the combo that would be
returned by a naive per-sample caller with a flat prior.

### 4b. Banded local exploration

From the GL-max combo, freebayes explores nearby combos by demoting a
handful of samples to their 2nd, 3rd, ... best genotypes ("banded"
exploration, `bandedGenotypeCombinations()` at
[Genotype.cpp:965-1099](../../freebayes/src/Genotype.cpp#L965-L1099)):

- **bandwidth** — how many samples to deviate simultaneously
- **banddepth** — how many rank-positions down each sample may drop
- at polyallelic sites with many alleles, the default is **exhaustive local
  search** (bandwidth = 0, banddepth = 0) to ensure GQ values are properly
  normalised.

### 4c. Iterate until convergence

After each exploration step the posterior is recomputed for every generated
combo. The top-ranked combo from this iteration becomes the seed for the
next iteration. Iteration stops when the best combo is the same two rounds
in a row (or when a per-site iteration cap is reached,
`--genotyping-max-iterations`, default 1000).

At convergence freebayes has a list of high-posterior combos. Everything
downstream uses this list.

**Geneticist intuition:** the search is *not* EM. There is no latent
allele-frequency variable being re-estimated each round. The iteration is
simply because the combo space is too large to enumerate, and the search
gradually refines the set of high-posterior combos.

---

## Stage 5: Marginalising to per-sample posteriors

Location: `marginalGenotypeLikelihoods()` in
[Marginals.cpp:41-92](../../freebayes/src/Marginals.cpp#L41-L92).

The joint posterior tells us which *combo* is most probable. But a VCF
expresses genotypes sample by sample. For sample `s`, what is the
probability of genotype `G`? Sum over all combos that assigned `G` to
sample `s`:

```
P(G_s = G | reads_all)  =  Σ_{combos c with G_s = G} P(c | reads_all)
```

Implementation ([Marginals.cpp:62](../../freebayes/src/Marginals.cpp#L62)):

```cpp
rmgs[sdl.genotype] = log( safe_exp(rmgsItr->second) + safe_exp(gc->posteriorProb) );
```

This is `logsumexp` — numerically stable addition of probabilities that are
stored as logs.

The sum is then normalised across all genotypes for that sample
([Marginals.cpp:81-83](../../freebayes/src/Marginals.cpp#L81-L83)):

```cpp
long double normalizer = logsumexp_probs(rawprobs);
long double newmarginal = marginals[sdl->genotype] - normalizer;
```

so the marginals sum to 1 (i.e. log-marginals exponentiate to a valid
probability distribution). The result lives on each `SampleDataLikelihood`
as `sdl->marginal`.

**These per-sample marginal posteriors drive GT and GQ.**

---

## Stage 6: Emitting VCF fields

Location: `ResultData::vcf()` in [ResultData.cpp](../../freebayes/src/ResultData.cpp).

### GT (genotype)
[ResultData.cpp:511](../../freebayes/src/ResultData.cpp#L511). The genotype
assigned to the sample in the **best combo** (the one with the highest
posterior, after marginal re-scoring). Written as `0/1`, `0/2`, etc.,
relative to the ref/alt alleles actually emitted at this site.

### GQ (genotype quality)
[ResultData.cpp:513-519](../../freebayes/src/ResultData.cpp#L513-L519):

```cpp
double val = big2phred(1 - big_exp(sampleLikelihoods.front().marginal));
```

GQ is a Phred-scaled measure of `P(genotype_call is wrong)` using the
**marginal posterior** of the best genotype:

```
GQ = -10 · log₁₀( 1 − P_marginal(best_genotype) )
```

High marginal posterior → very small `1 − p` → large GQ. A GQ of 30 means a
1-in-1000 chance the call is wrong (under the model).

### GL (genotype log-likelihoods)
[ResultData.cpp:534-631](../../freebayes/src/ResultData.cpp#L534-L631).
Emits, in the VCF-standard ordering `F(j/k) = k(k+1)/2 + j`, the log₁₀ of
the **data likelihood** per genotype:

```cpp
genotypeLikelihoods[o->second] = ln2log10(g->prob);
...
genotypeLikelihoodsOutput[g->first] = convert(g->second - maxGL);
```

Note `g->prob` is the raw data likelihood from Stage 2, *not* the marginal
posterior. The values are max-normalised so the best genotype has GL = 0 and
the others are negative.

**Geneticist viewpoint:** GL in freebayes is the prior-free likelihood. If
you want to re-do downstream joint analysis with a different prior (e.g.
merging freebayes output into a larger cohort), the GL field is the honest,
reusable quantity. GT and GQ already bake in the site-level population prior
and should not be re-interpreted outside the called cohort.

### QUAL (site quality)
Around [freebayes.cpp:534-546](../../freebayes/src/freebayes.cpp#L534-L546):

```
pVar = 1 − Σ_{homozygous-reference combos} P(combo | reads)
QUAL = -10 · log₁₀(1 − pVar)       (approximately)
```

i.e. QUAL answers "what is the posterior probability that this site is
variant?" computed as 1 minus the total posterior weight of all combos in
which every sample is homozygous reference.

### Depth fields
`DP`, `AD`, `RO`, `AO`, `QR`, `QA` are set from the per-sample observation
tallies kept alongside the probability calculations
([ResultData.cpp:521-532](../../freebayes/src/ResultData.cpp#L521-L532)).
These are pure counts and qualities; they do not use any probabilities.

---

## End-to-end walkthrough: a 3-sample diploid SNP

Suppose we are at position 100 of chromosome 1 and three samples have the
following observations (base quality ~30, mapping quality ~60 for all
reads):

```
sample   A reads   T reads
S1       20          0       # looks like A/A
S2       10         10       # looks like A/T
S3        0         18       # looks like T/T
```

### What happens at this site

1. **BAM → observations (Stage 1).** Each read becomes an `Allele`
   observation. The candidate allele set is `{A, T}` because both pass
   `--min-alternate-count` and `--min-alternate-fraction`. No haplotype
   extension is needed — this is a simple SNP. Candidate genotypes per
   sample: `A/A`, `A/T`, `T/T`.
2. **Per-sample likelihoods (Stage 2).** For each sample and each genotype
   we compute `log P(reads_s | G)`:
   - S1 strongly favours `A/A` — all 20 reads consistent, product of their
     `1 − 10⁻³` ≈ near-zero log penalty. Likelihoods for `A/T` and `T/T`
     are much lower.
   - S2 favours `A/T` — the 50/50 balance matches the heterozygote
     sampling probability 0.5, and the multinomial allele-balance term
     [DataLikelihood.cpp:155](../../freebayes/src/DataLikelihood.cpp#L155)
     peaks exactly there.
   - S3 strongly favours `T/T` — symmetric to S1.
3. **Combo construction (Stage 3).** There are 3³ = 27 possible combos.
   The convergent search starts from the GL-max combo
   `(A/A, A/T, T/T)`, which also happens to be the posterior max.
   Nearby combos (e.g. `(A/A, A/A, T/T)`, which drops S2 to `A/A`) are
   also evaluated.
4. **Combo priors.**
   - For `(A/A, A/T, T/T)`: allele counts across samples are A:3, T:3.
     The Ewens prior is moderate (not extreme because the frequency is
     0.5, not rare). The HWE prior is fine — with AF = 0.5 we expect
     genotype frequencies 0.25/0.5/0.25, and our (1, 1, 1) matches that
     reasonably well. Multinomial-coefficient prior for the het S2 gets
     a small boost. Strand/placement priors depend on the reads but are
     likely neutral.
   - For `(A/A, A/A, T/T)`: allele counts A:4, T:2. Ewens prior is
     *better* (rarer alt allele), but the HWE term punishes this combo
     because we now have no heterozygotes despite f(T) = 1/3.
     Critically, the likelihood for S2 at `A/A` is much worse than at
     `A/T` — it has to explain away 10 T reads as errors, each costing
     ~log(10⁻³).
   - Net: the correct combo `(A/A, A/T, T/T)` wins.
5. **Marginalisation (Stage 5).** For S1, summing over all combos that
   assign `A/A` to S1 (most of the high-posterior combos) concentrates
   nearly all the posterior mass on `A/A`. Same for S3 on `T/T`. For S2,
   `A/T` gets nearly all the mass; a tiny fraction goes to `A/A` and
   `T/T` from low-posterior combos. GQ for S1 and S3 is very high; GQ
   for S2 is also high, but slightly lower.
6. **VCF output (Stage 6).**
   - QUAL is high: the posterior weight on "all hom-ref" combos is tiny.
   - GT: `0/0`, `0/1`, `1/1` for S1/S2/S3.
   - GQ: large for all three.
   - GL: for S1, roughly `(0, large_negative, very_large_negative)`
     because the A/A likelihood dominates.

### What if S2 had only 1 T read out of 20?

The per-sample likelihood for S2 now favours `A/A`, but weakly — one T
read is not impossible under `A/A` if we believe in base errors. With S1
and S3 as above, the combo `(A/A, A/A, T/T)` is now plausible. But the
Ewens prior still prefers combos where the alt allele is real and rare
(because S3 is T/T, freebayes already "knows" T exists at this site).
The outcome depends on the balance between the T-read penalty in S2's
likelihood and the Ewens prior on how common T is — this is exactly the
kind of tension the joint model is designed to resolve.

---

## freebayes vs GATK (for cross-reference)

| Aspect | freebayes | GATK GenotypeGVCFs |
|---|---|---|
| How priors are iterated | Joint search over genotype combos; no EM | EM over latent allele frequencies |
| Allele-frequency prior | Ewens' sampling formula, `θ` (default 0.001) | Dirichlet on allele frequencies, pseudocounts derived from heterozygosity |
| Per-sample likelihoods | Computed from BAM reads, at the moment of calling | Pre-computed at sample time in HaplotypeCaller, shipped as PL in gVCF |
| Unit of variation | Haplotype allele (complex variants preserved as single allele) | Site-level alleles, reassembled from HaplotypeCaller output |
| GL field meaning | Same: log₁₀ P(reads | G), no prior | Same |
| GQ field meaning | Same: Phred of 1 − P_marginal(best) | Same |
| Joint calling boundary | Operates per site, per cohort, using raw reads | Operates per site, per cohort, using pre-computed PLs |

Both approaches converge on the same *conceptual* answer — a per-sample
posterior over genotypes conditional on all reads in the cohort — but by
very different machinery. The key practical consequence for pipelines is
that freebayes re-reads the BAMs each run, while GATK reuses per-sample
gVCF PLs. freebayes therefore needs the BAMs at joint-calling time; GATK
does not.

---

## The data freebayes actually consumes

This section is a reference inventory: it lists, at the finest useful
granularity, every piece of per-sample / per-allele / per-read quantity
that freebayes reads out of the BAMs and flows into the posterior
calculation. If you wanted to reproduce freebayes' answer without
re-reading BAMs — for example, by caching a per-sample summary file like
the `.psf` proposed in
[calling_pipeline_architecture.md](calling_pipeline_architecture.md) — this
is the list of quantities that summary must either store directly or be
able to reconstruct.

### Per-read quantities (the most granular)

For each read that overlaps a candidate locus in a given sample, the
likelihood calculation at
[DataLikelihood.cpp:5-166](../../freebayes/src/DataLikelihood.cpp#L5-L166)
uses:

| Quantity | Source | How it is used |
|---|---|---|
| base quality `Q` at the variant position | BAM `QUAL` string, per position | `lnquality = ln(10^(-Q/10))` — probability this base is a sequencing error |
| mapping quality `MQ` | BAM `MAPQ` field | `lnmapQuality` — probability the read is mismapped; combined with `Q` via `max()` in log space |
| which allele the read supports | CIGAR + sequence vs. reference | decides whether the read is "consistent with `G`" or "an error under `G`" |
| partial-support haplotype assignment | read re-alignment to candidate haplotypes | when a read only covers part of the haplotype window, its evidence is distributed across the haplotypes it is consistent with |
| strand (forward / reverse) | BAM flag 0x10 | fed into per-combo strand-bias tallies (Component 4 prior) |
| position in read / fragment placement | read coordinate arithmetic | fed into per-combo placement-bias tallies |
| read group ID | BAM `RG` tag | looked up against the contamination table for sample-specific ref-bias correction |
| fragment identity (paired-end linkage) | BAM `QNAME` / flags | used by the RDF correction to down-weight correlated errors from the same fragment |

**Note.** freebayes does not need the read's full sequence outside the
haplotype window — only the allele assignment, the two qualities, and the
strand/position metadata. This is the source of the compressibility that
the PSP per-read record
([per_sample_pileup_format.md](per_sample_pileup_format.md)) exploits by
packing `(BQ, MQ, position_in_read, strand)` into 3 bytes.

### Per-sample, per-allele aggregates

Once per sample per candidate allele, the likelihood stage aggregates the
per-read data into:

| Quantity | Where computed | Where used |
|---|---|---|
| observation count per allele | `Sample` accumulator during BAM parsing | `multinomialSamplingProbLn(alleleProbs, observationCounts)` — the allele-balance term in the likelihood at [DataLikelihood.cpp:155](../../freebayes/src/DataLikelihood.cpp#L155) |
| sum over consistent reads of `log(asampl · scale)` | accumulated in `prodSample` | the per-haplotype sampling probability contribution |
| sum over inconsistent reads of `log(1 − qual)` | accumulated in `prodQout` | the "these reads must be errors" contribution |

These aggregates are per sample per genotype hypothesis, so they are
recomputed each time a new candidate genotype is evaluated.

### Per-cohort, per-allele aggregates (drive the priors)

The combo-level prior at
[Genotype.cpp:1499-1605](../../freebayes/src/Genotype.cpp#L1499-L1605)
depends on pools of observations across **all samples** in the current
combo. The `AlleleCounter` struct ([Genotype.h](../../freebayes/src/Genotype.h))
maintains, per candidate allele in the combo:

| Quantity | Drives which prior |
|---|---|
| `frequency` — total allele copies across all samples' genotypes | Component 1 (multinomial coefficient) and Component 2 (Ewens) and Component 3 (HWE) |
| `observations` — total reads supporting this allele across all samples | Component 4 (observation-bias and allele-balance) |
| `forwardStrand` / `reverseStrand` counts | `binomialProbln(forwardStrand, obs, 0.5)` — strand-bias prior |
| `placedLeft` / `placedRight` counts | `binomialProbln(placedLeft, obs, 0.5)` — read-placement bias |
| `placedStart` / `placedEnd` counts | `binomialProbln(placedStart, obs, 0.5)` — fragment-position bias |

**This is the critical subset for cohort-level posterior re-calculation.**
If you have per-sample, per-allele observation summaries with these
counts, you can reconstruct the combo-level `AlleleCounter` fields by
summing across samples — no BAM re-read required. This is exactly why
the `.psf` contract proposes storing strand, placement, and start counts
per allele per sample.

Derived combo-level quantities (computed from the above, not stored):

- `alleleFrequencyCounts` — the partition `{allele_count → number_of_alleles_with_that_count}`, the input to Ewens' formula at [Ewens.cpp:32-51](../../freebayes/src/Ewens.cpp#L32-L51).
- `allele_counts_in_each_sample's_genotype` — needed for the HWE
  multinomial coefficients at
  [Genotype.cpp:1417-1493](../../freebayes/src/Genotype.cpp#L1417-L1493).

### Per-sample metadata

| Quantity | Source | Where used |
|---|---|---|
| sample name | BAM `SM` read-group tag | VCF column, sample lookups |
| ploidy | `--ploidy` or `--cnv-map` | `getGenotypesByPloidy()` — determines the candidate genotype set |
| population label (optional) | `--populations` file | allows independent combo searches per population, priors computed within-population |

### Per-site metadata

| Quantity | Source | Where used |
|---|---|---|
| reference base(s) in the haplotype window | reference FASTA | likelihood (allele matching) and prior (the `<REF>` anchor) |
| haplotype window length | `AlleleParser::lastHaplotypeLength` | scales `θ` in the Ewens prior at [freebayes.cpp:309](../../freebayes/src/freebayes.cpp#L309): `theta = parameters.TH * lastHaplotypeLength` |
| candidate allele set and their haplotype sequences | discovered at site from pooled sample observations | the genotype space that everything else is defined over |

### What freebayes does *not* use

Worth explicitly listing, because it prunes the `.psf` contract:

- Read sequence outside the haplotype window — only the allele
  assignment + qualities matter.
- Pre-filter read counts (duplicates, low-MAPQ, secondary alignments are
  dropped before anything reaches the likelihood calculation).
- CIGAR string detail beyond what is needed to assign the read to an
  allele and compute its base quality at the variant position.
- Insert size, template orientation, PE pairing geometry beyond strand.

### Compactness implications — important subtlety about BQ and MQ

A naive reading of the tables above suggests that per-read base and
mapping qualities could be replaced by two scalars per allele:
`sum(BQ)` and `sum(MQ)`. **This is wrong.** freebayes combines BQ and
MQ *per read* via nonlinear functions before aggregating:

- Standard-GL path: `prodQout += max(ln_BQ_read, ln_MQ_read)` per read
  that disagrees with the genotype hypothesis
  ([DataLikelihood.cpp:34](../../freebayes/src/DataLikelihood.cpp#L34)).
- Non-standard-GL path: `qual_read = (1 − exp(ln_BQ_read)) · (1 − exp(ln_MQ_read))`
  per read, then `prodQout += log(1 − qual_read)`
  ([DataLikelihood.cpp:77,135](../../freebayes/src/DataLikelihood.cpp#L77)).

In both cases the per-read contribution is a nonlinear function of the
pair `(BQ, MQ)`. You cannot recover it from the per-allele marginal sums
of BQ and MQ separately. That nonlinearity is why the naive summary fails.

There are two clean options that do work:

**Option A — keep per-read records.** Store `(BQ, MQ, strand, position_in_read)`
per read per allele per sample (the PSP v1 approach at
[per_sample_pileup_format.md](per_sample_pileup_format.md), 3 bytes per
read). Any likelihood formula — freebayes' current one, a future revision,
or a different model entirely — remains computable.

**Option B — pre-compute a freebayes-specific scalar per allele.** Store
`S_allele = Σ over reads supporting this allele of max(ln_BQ_read, ln_MQ_read)`
(standard-GL path) or the corresponding sum of `log(1 − qual_read)` for
the non-standard path. Given that single scalar plus the observation
count, freebayes' `prodQout` contribution is reconstructable exactly. But
the summary is committed to freebayes' formula — switching to a different
likelihood model later would require per-read data you no longer have.

### Minimum sufficient per-sample data

Given the choice above, the minimum per-sample per-variant-locus
quantities to reproduce freebayes' posterior-driving aggregates *without
re-reading the BAM* are:

1. The candidate alleles at the locus (discovered locally from this
   sample's reads).
2. Per-allele bias counts: observation count, forward-strand count,
   placed-left count, placed-start count. These are simple sums — they
   aggregate linearly across reads, samples, and combos, so they safely
   summarise to a scalar per allele per sample.
3. Per-allele quality summary, chosen per Option A or Option B above:
   - **Option A:** per-read `(BQ, MQ, strand, position_in_read)` list.
   - **Option B:** pre-computed scalar `Σ max(ln_BQ, ln_MQ)` per allele
     (freebayes-specific).
4. Per-allele haplotype context (sequence + anchor offsets).
5. Sample ploidy.

For candidate genotype likelihoods (the thing that becomes GL in the
VCF), additionally one PL value per genotype against the sample's local
allele set — these are what freebayes computes and writes as GL today.

The bias counts (item 2) reconstruct Components 1–4 of the combo prior
exactly regardless of Option A vs B. Item 3 reconstructs the per-sample
likelihood (`prodQout` and `prodSample`). Items 4–5 let the merger
decide, at cohort-assembly time, whether a compound haplotype allele
introduced by other samples is consistent with this sample's
observations.

This is the key trade-off behind Open decision 1 in
[calling_pipeline_architecture.md](calling_pipeline_architecture.md):
flexibility for future modelling (Option A) vs. minimum file size
(Option B).

---

## Summary

freebayes goes from BAM to posterior probabilities in six conceptual
stages:

1. **Parse reads** into per-sample allele observations with base, strand,
   and quality metadata.
2. **Compute a per-sample likelihood** `P(reads | G)` for every candidate
   genotype `G`, using base and mapping qualities plus an allele-balance
   multinomial term, dampened by an RDF correction for correlated errors.
3. **Enumerate joint genotype combos** across all samples and compute the
   combo-level prior: P(genotypes | allele frequency) × Ewens P(AF) × HWE
   × observation-bias priors. Multiply by the joint likelihood.
4. **Search** the combo space (banded, convergent) rather than enumerate
   it, because it is too large.
5. **Marginalise** over combos to obtain a per-sample posterior
   distribution over genotypes.
6. **Emit VCF fields**: GT from the best combo, GQ from the marginal, GL
   as the prior-free per-sample log-likelihoods, QUAL from the total
   posterior weight on non-reference combos.

The whole pipeline stays in log space. The population-genetic content of
the calls lives almost entirely in the priors — especially Ewens' formula
(rare alleles are more plausible) and the HWE term (genotype counts should
match the implied allele frequencies). Without those priors the calls
collapse to naive per-sample maximum-likelihood, which cannot distinguish
real rare variants from stochastic sequencing error.
