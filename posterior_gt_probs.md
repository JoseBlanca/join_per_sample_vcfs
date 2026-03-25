# GATK GenotypeGVCFs: Genotype Posterior Probability Algorithm

This document describes how GATK calculates genotype and variant posterior
probabilities during the joint genotyping step (GenotypeGVCFs), where individual
per-sample gVCFs are merged into a final multi-sample VCF.

## Call Flow

```
GenotypeGVCFs.apply()
  -> GenotypeGVCFsEngine.callRegion()
    -> merge per-sample gVCFs
    -> regenotypeVC()
      -> GenotypingEngine.calculateGenotypes()
        -> AlleleFrequencyCalculator.calculate()   <- core math
```

## Key Source Files

All under `src/main/java/org/broadinstitute/hellbender/` in the
[GATK repo](https://github.com/broadinstitute/gatk):

| File | Role |
|---|---|
| `tools/walkers/GenotypeGVCFs.java` | Tool entry point |
| `tools/walkers/GenotypeGVCFsEngine.java` | Orchestrates merge + regenotyping |
| `tools/walkers/genotyper/GenotypingEngine.java` | Base class with `calculateGenotypes()` |
| `tools/walkers/genotyper/afcalc/AlleleFrequencyCalculator.java` | The core algorithm |
| `utils/Dirichlet.java` | Dirichlet distribution for allele frequency priors |

## The Mathematical Model

Uses a Bayesian model with a Dirichlet prior on allele frequencies and EM
(Expectation-Maximization) iteration.

### 1. Dirichlet Prior (pseudocounts)

With defaults (`snpHeterozygosity = 0.001`, `heterozygosityStandardDeviation = 0.01`,
`indelHeterozygosity = 1.25e-4`):

```
refPseudocount   = snpHet / hetStdDev^2  = 0.001 / 0.0001 = 10
snpPseudocount   = snpHet * refPseudo    = 0.001 * 10      = 0.01
indelPseudocount = indelHet * refPseudo  = 1.25e-4 * 10    = 0.00125
```

The prior is heavily biased toward the reference allele.

### 2. EM Iteration

**Initialize**: flat allele frequencies (`1/numAlleles` for all alleles).

**E-step**: for each sample, compute genotype posteriors given current allele
frequencies:

```
log10_posterior(g) = log10_multinomial_coeff(g)
                   + log10_likelihood(g)           <- from PL field in gVCF
                   + SUM(count_i * log10(freq_i))  <- HW prior given frequencies
```

- **Multinomial coefficient** = `ploidy! / (count_1! * count_2! * ...)`
  (e.g., diploid het = 2, hom = 1).
- **Likelihood** comes directly from the PL field stored in each sample's gVCF
  (computed earlier by HaplotypeCaller).

**M-step**: compute effective allele counts by summing posterior-weighted allele
counts across all samples:

```
effective_count[a] = SUM_samples SUM_genotypes ( posterior(g) * copies_of_a_in_g )
```

Update allele frequencies via Dirichlet conjugate update:

```
posterior_pseudocounts = prior_pseudocounts + effective_allele_counts
freq[a] = posterior_pseudocounts[a] / SUM(posterior_pseudocounts)
```

**Convergence**: iterate until max change in any allele count < 0.1. The first
iteration uses flat frequencies (not the prior) to avoid suppressing real
variants before they accumulate evidence.

### 3. Variant Quality (QUAL)

```
log10_P(no variant) = SUM_samples log10_posterior(hom_ref)
QUAL = -10 * log10_P(no variant)
```

### 4. Individual Genotype Assignment

Each sample's final genotype is the one with the highest posterior probability
(PL-based, `PREFER_PLS` method), after subsetting to the final allele set.

## EM Iteration: Intuitive Explanation

### What problem is the EM solving?

When you run HaplotypeCaller on a single sample, it produces genotype
likelihoods (the PL field) at each site. These likelihoods tell you how well
the sequencing reads match each possible genotype (0/0, 0/1, 1/1, etc.), but
they do not tell you the final genotype call. To make a confident call, you
also need to know how common each allele is in the population being sequenced.

For example, if a site has weak evidence for a het (0/1) in one sample, you
cannot tell whether it is a real variant or a sequencing artifact. But if 20
other samples in your cohort also show evidence for the same alt allele, even
weak evidence in your sample becomes meaningful. The EM algorithm is the
mechanism that pools evidence across all samples to estimate allele frequencies,
and then uses those frequencies to refine each sample's genotype call.

### The chicken-and-egg problem

There is a circular dependency:

- To know each sample's genotype, you need to know the allele frequencies in
  the cohort (because a common allele is more likely to be real).
- To know the allele frequencies, you need to know the genotypes of all
  samples.

The EM algorithm breaks this circle by alternating between two steps, each time
improving the estimates, until they stabilize.

### Step by step with a concrete example

Imagine a site where the reference allele is A and a possible alt allele is T.
You have 100 diploid samples (200 chromosomes total).

#### Starting point

You do not yet know the allele frequency of T, so you start with a naive guess:
assume A and T are equally common (50% each). This is intentionally flat to
avoid biasing the first round.

#### E-step ("Expectation"): assign provisional genotypes

Using the current allele frequency guess, compute for each sample the
probability of each possible genotype:

- **P(0/0 | data, freq)**: probability of being hom-ref, considering both the
  sequencing evidence (PL) and how common the ref allele currently appears.
- **P(0/1 | data, freq)**: probability of being het.
- **P(1/1 | data, freq)**: probability of being hom-alt.

These are not hard calls, but soft probabilities. A sample with ambiguous reads
might get P(0/0) = 0.6, P(0/1) = 0.35, P(1/1) = 0.05.

The key formula combines three things:

1. **Genotype likelihood from sequencing** (the PL field from the gVCF): how
   well do the reads match this genotype?
2. **Hardy-Weinberg expectation given current frequencies**: if freq(T) = 0.05,
   then P(0/1) = 2 * 0.95 * 0.05 = 0.095, P(1/1) = 0.05^2 = 0.0025. This is
   where the allele frequency acts as a prior: rare alleles make het and
   hom-alt genotypes less likely a priori.
3. **Normalize**: scale so the three probabilities sum to 1.

#### M-step ("Maximization"): update allele frequencies

Now count how many copies of each allele are present across the cohort, but
using the soft genotype probabilities instead of hard calls.

For each sample, the expected number of T alleles is:

```
expected_T = P(0/1) * 1 + P(1/1) * 2
```

A sample with P(0/0) = 0.6, P(0/1) = 0.35, P(1/1) = 0.05 contributes
0.35 * 1 + 0.05 * 2 = 0.45 copies of T.

Sum these expected counts across all 100 samples and divide by 200 (total
chromosomes) to get the updated allele frequency for T. If the sum is 9.0, then
freq(T) = 9.0 / 200 = 0.045.

At this point, the Dirichlet prior also contributes its pseudocounts (10 for
ref, 0.01 for each SNP alt), slightly pulling the estimate toward rare-variant
expectations. In practice, with 100+ samples, the data overwhelms the prior.

#### Repeat until stable

Go back to the E-step with the new freq(T) = 0.045, recompute all genotype
probabilities, then recount alleles (M-step), and so on. Typically converges in
a few iterations.

After the first round with flat frequencies, the algorithm might find freq(T) =
0.045. In the next E-step, samples with weak T evidence will now get lower
P(0/1) (because T is rarer than the initial 50% guess), while samples with
strong T evidence will keep high P(0/1). The M-step then gives a more accurate
freq(T). After 3-5 rounds, the estimates settle.

### What happens after convergence

Once allele frequencies and genotype posteriors are stable:

1. **QUAL score**: multiply together all samples' P(0/0) values. If every
   sample has high P(0/0), the product is close to 1 and QUAL is low (no
   variant). If even a few samples have low P(0/0), the product drops and QUAL
   is high. QUAL = -10 * log10(product of P(hom-ref) across all samples).

2. **Per-sample genotype (GT)**: each sample is assigned the genotype with the
   highest posterior probability.

3. **Genotype quality (GQ)**: reflects how confident the call is for that
   individual sample, derived from the difference between the best and
   second-best genotype posteriors.

### Why this matters in practice

- **Rare variants in large cohorts**: a variant seen weakly in 3 out of 1000
  samples can be called confidently because the EM pools evidence. Without joint
  calling, each of those 3 samples might be called hom-ref due to insufficient
  individual evidence.

- **Common variants are easy**: if freq(T) = 0.30, the HW prior strongly
  favors het and hom-alt genotypes, so even modest read evidence suffices.

- **Sequencing errors are suppressed**: random errors do not correlate across
  samples, so they do not accumulate in the allele frequency estimate. A true
  variant produces a consistent signal across multiple samples; an error does
  not.

## Summary

The PL values from individual gVCFs are combined implicitly through the EM
procedure: each sample's PLs contribute independently to effective allele
counts, and the shared allele frequency estimate ties samples together. This is
equivalent to assuming samples are drawn independently from a population with
shared allele frequencies (independent-samples model under Hardy-Weinberg).
