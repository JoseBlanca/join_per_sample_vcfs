# GATK vs. freebayes — approach comparison

A side-by-side comparison of the two reference approaches to multi-sample
variant calling this project draws from. Written to make the design
trade-offs legible when deciding what to keep from each.

For the full algorithmic details, see:

- [posterior_gt_probs.md](posterior_gt_probs.md) — GATK `GenotypeGVCFs`
  posterior math
- [freebayes_posterior_gt_probs.md](freebayes_posterior_gt_probs.md) —
  freebayes end-to-end BAM → posterior pipeline

This document compares them; it does not re-explain either in detail.

## The shared goal

Both tools compute, per sample `s` at each candidate site, the posterior
distribution over genotypes `G`:

```
P(G_s | reads from all samples in the cohort)
```

They agree on the end product (GT, GQ, QUAL, GL/PL in the VCF). They
differ radically in *how* they get there.

## Pipeline shape

| | GATK | freebayes |
|---|---|---|
| Per-sample stage | HaplotypeCaller (HC) → per-sample gVCF | none |
| Joint stage | GenotypeGVCFs (GGVCF) reads gVCFs → multi-sample VCF | reads all BAMs directly → multi-sample VCF |
| BAMs needed at joint calling time | no | **yes** |
| Per-sample artefact reusable across runs | yes (the gVCF) | no |
| Typical deployment | split: HC run once per sample, GGVCF run per cohort | single invocation across the cohort |

This shape difference is the root of everything else.

## Per-sample stage

### GATK HaplotypeCaller

1. **Active region detection.** Scans for regions with any sign of
   variation (reads showing mismatches, soft-clips, indels).
2. **Local haplotype assembly.** Builds a de Bruijn graph from reads in
   the region and traverses it to enumerate candidate haplotypes.
3. **PairHMM realignment.** For each read against each candidate
   haplotype, computes `P(read | haplotype)` with a pair hidden Markov
   model that models base-call errors, insertion opens/extensions,
   deletion opens/extensions, and mapping quality cap. This is the most
   compute-intensive part.
4. **Per-variant genotype likelihoods.** Marginalises the per-read
   per-haplotype likelihoods to per-read per-allele, combines those into
   per-sample per-genotype likelihoods against the local allele set.
5. **gVCF emission.** Writes one record per variant site plus
   reference-confidence blocks between variants. Every variant record
   includes:
   - alleles `[REF, alt1, ..., altN, <NON_REF>]`
   - PL values for every genotype against those alleles, **including
     genotypes involving `<NON_REF>`** — the likelihood assuming the
     sample has some unspecified alternate allele

### freebayes

There is no dedicated per-sample stage. freebayes walks the BAMs of all
samples simultaneously. At each candidate site it:

1. Pools per-read observations across samples to discover candidate
   alleles.
2. Extends the window into haplotype alleles if events overlap (Stage 1c
   in the freebayes doc).
3. Computes per-sample likelihoods on the fly using a simpler per-read
   error model: `max(ln_BQ, ln_MQ)` in log space, plus allele-balance
   multinomial.

### Key differences

| | GATK | freebayes |
|---|---|---|
| Read-haplotype likelihood model | PairHMM (HMM with indel states) | log-space max of base-qual and mapping-qual error probs |
| Indel modelling | explicit via HMM state transitions | implicit — indels are just another allele |
| Requires pre-alignment | yes | yes |
| Persistent per-sample output | yes (gVCF) | no |
| Compute per sample | ~30–60 min for WGS | n/a (fused into joint stage) |

## Allele representation

This is where the two most visibly diverge.

### GATK

HaplotypeCaller assembles haplotypes internally but **decomposes them
into site-level variants** when writing the gVCF. A haplotype carrying a
3 bp insertion plus a downstream SNP becomes two independent records in
the gVCF: one for the insertion, one for the SNP. The fact that they
were on the same haplotype is lost at gVCF boundary.

At joint calling, `GenotypeGVCFs` unifies alleles **per site** across
samples. It has no mechanism to reconstruct a compound haplotype that
spans originally-separate variant records.

### freebayes

Haplotype-level alleles are preserved end-to-end. A 3 bp insertion + SNP
on the same haplotype is a single allele throughout the pipeline. At
joint calling freebayes can reason about compound haplotype alleles
natively because they never got decomposed.

### Consequence

At merge-heavy sites — overlapping deletions, deletion+SNP combinations,
adjacent indels — **freebayes calls compound variants correctly;
`GenotypeGVCFs` does not**. This is the well-known limitation that
motivates merging rather than joining.

## Handling allele-set mismatches at joint calling

### GATK's `<NON_REF>` trick

At joint calling time, sample A's gVCF may contain an alt allele that
sample B's gVCF does not (because B's reads didn't support it). To assign
a likelihood to "sample B has that allele" without re-reading B's BAM,
`GenotypeGVCFs` uses B's PL value against the `<NON_REF>` symbolic allele
as a proxy: "the likelihood of B's reads assuming B has some unspecified
allele at this site."

This is good enough at simple SNP/indel sites. It is the main reason
GATK's per-sample → join workflow is feasible at all.

### freebayes

Not applicable. The problem doesn't exist because freebayes never
separates per-sample analysis from joint calling — all reads are in
memory together, so the allele set is unified before likelihoods are
computed.

## Prior structure

### GATK (`AlleleFrequencyCalculator`)

| Component | Form |
|---|---|
| Allele frequency prior | Dirichlet with pseudocounts (`ref=10, alt_snp=0.01, alt_indel=0.00125`) |
| Genotype prior given AF | multinomial (implicit via Dirichlet-multinomial conjugacy) |
| HWE | not explicit; flows from the Dirichlet-multinomial |
| Observation-bias priors (strand, placement) | **none in the posterior**; handled as post-hoc QC annotations (SOR, FS, etc.) |

### freebayes (`calculatePosteriorProbability`)

| Component | Form |
|---|---|
| Allele frequency prior | Ewens' sampling formula, `θ ≈ 0.001` scaled by haplotype length |
| Genotype prior given AF | explicit multinomial coefficient |
| HWE | optional, on by default |
| Observation-bias priors | **in the posterior**: binomial priors on strand balance, placement left/right, fragment-start position |

### Key differences

- Both encode "rare variants are *a priori* more plausible than common
  ones." GATK does it with a tiny alt Dirichlet pseudocount; freebayes
  does it with Ewens' neutral-theory partition distribution. For rare
  variants the two behave similarly.
- freebayes' observation-bias priors are a genuine structural addition:
  a variant whose alt-supporting reads are all on one strand gets
  penalised in the posterior itself, not just at post-hoc filtering.
  This matters most at borderline sites.

## Breaking the chicken-and-egg

The circular dependency — per-sample posteriors depend on allele
frequency, allele frequency depends on per-sample genotypes — is solved
differently.

### GATK: Expectation-Maximisation

Allele frequency `p` is treated as a latent parameter.

- **E-step.** For current `p`, compute soft per-sample genotype
  posteriors: `P(G_s | PL_s, p)`.
- **M-step.** From those soft posteriors, recount expected allele copies
  and update `p` via the Dirichlet conjugate update.
- **Iterate** until `p` stabilises (usually 3–5 rounds).
- **Output.** Plug the final `p̂` into the E-step formula to get the
  per-sample posteriors reported in the VCF.

This is **empirical Bayes**: a single best allele frequency is estimated
from the data and used as a plug-in prior.

### freebayes: joint combo search + marginalisation

Enumerate / search joint assignments of genotypes to all samples (each
such assignment is a "genotype combo"). For each combo:

- The implied allele frequency is whatever that combo's genotypes say.
- The prior and likelihood are both well-defined.
- Compute the combo's posterior = likelihood × prior.

To recover per-sample posteriors, sum combo posteriors into each
sample's buckets: `P(G_s = G | reads) = Σ_{combos assigning G to s} P(combo | reads)`.

This is an approximate **full Bayesian integration**: the sum over combos
implicitly integrates over allele frequency (each combo implies its own
`p`).

### Consequences

| | GATK EM | freebayes combo search |
|---|---|---|
| Allele frequency treatment | point estimate `p̂` | implicitly integrated over |
| Calibration at common variants | excellent | excellent |
| Calibration at rare variants with small N | over-confident | better |
| Calibration with large N | indistinguishable | indistinguishable |
| Iterations per site | 3–5 | 10–50 |
| Scaling with sample count | linear | sub-exponential |
| Works with thousands of samples | yes | no |

At your target of thousands of samples, `p̂` becomes extremely tight and
the two approaches converge to practically identical numbers. At small
cohorts (< ~100) freebayes is better-calibrated.

See [freebayes_posterior_gt_probs.md §Three ways to solve the same problem](freebayes_posterior_gt_probs.md#three-ways-to-solve-the-same-problem-freebayes-gatk-mcmc)
for a longer treatment including MCMC as a third reference point.

## Output fields

Both tools emit the same VCF fields with the same semantics. The
differences are in what each field refers to:

| Field | GATK | freebayes |
|---|---|---|
| GT | argmax of marginal posterior | argmax of marginal posterior |
| GQ | Phred of `1 − P_marginal(best)` | Phred of `1 − P_marginal(best)` |
| PL / GL | log-likelihoods against sample-local alleles + `<NON_REF>`, Phred-scaled | log₁₀-likelihoods against the cohort-merged allele set |
| QUAL | from total posterior weight on non-reference outcomes | same, from combo posteriors |

A subtle but important difference: **GATK's PLs are reusable at future
join time** (that's their whole purpose). **freebayes' GLs are not**,
because they are defined against whatever allele set happened at call
time. Re-calling a cohort with a different sample set changes the
allele set and invalidates the GLs.

## How each tool uses PL / GL

GATK and freebayes both emit genotype likelihoods in the VCF, but they
play fundamentally different roles in each pipeline. This is arguably
the single biggest architectural consequence of the split vs. joint
topology.

### In GATK: PL is the load-bearing intermediate

PLs are **computed once** at HaplotypeCaller time and **consumed
repeatedly** at GenotypeGVCFs time. The entire joint-calling step
operates on PLs alone; the BAMs are never opened again.

The lifecycle:

```
reads ──(PairHMM)──► P(read|haplotype)
                           │
                           ▼
                 per-sample genotype likelihood
                           │
                           ▼
                    normalise, Phred-scale
                           │
                           ▼
                 PL in gVCF (against [REF,alt1,...,altN,<NON_REF>])
                           │
                           │  ← gVCF is stored; HaplotypeCaller done
                           │
                           ▼
                 read by GenotypeGVCFs at join time
                           │
                           ▼
           EM combines PL + HW prior + allele-freq estimate
                           │
                           ▼
                 per-sample posterior → GT, GQ, QUAL
                           │
                           ▼
                 PL is also *copied verbatim* into the final VCF
```

Notable points:

- **PL is the sufficient statistic.** Once PLs are written, HC's
  PairHMM work, read alignments, and quality scores are all discarded as
  far as joint calling is concerned. Everything downstream is a function
  of PL.
- **PL survives the joint stage unchanged.** The final multi-sample VCF
  carries the same PL values that were in the gVCF. Only GT, GQ, and
  QUAL change, computed from PL + the EM-inferred allele frequency.
- **PL is defined against a sample-local allele set**, extended by
  `<NON_REF>`. At join time, when alleles from other samples get added
  to the unified allele set, `<NON_REF>` PLs stand in as likelihoods for
  the new alleles.
- **Reusable across runs.** Re-genotype the same cohort with different
  filters, or add samples and recall — no BAMs needed, no HC re-run.

### In freebayes: GL is a write-only output

freebayes computes per-sample likelihoods internally at joint-call time
from the BAMs, uses them inside the combo search, and writes them out
as GL. But **freebayes never reads a GL field in**. GLs only appear in
the output VCF; they are not part of any input format freebayes accepts.

The lifecycle:

```
BAMs of all samples  ──►  freebayes (single invocation)
                              │
            ┌─────────────────┼─────────────────┐
            │                 │                 │
            ▼                 ▼                 ▼
   pool observations    compute per-sample      combo search
   across samples       likelihoods from        (prior × likelihood)
   (haplotype-level     base/map quality                │
   allele set)          + allele balance                │
                                                        ▼
                                              per-combo posterior
                                                        │
                                                        ▼
                                              marginalise to per-sample
                                                        │
                                                        ▼
                                          GT, GQ, QUAL  +  GL (= per-sample
                                                            likelihoods,
                                                            log₁₀, max=0)
```

Notable points:

- **GL is derived, not durable.** It is a projection of the transient
  in-memory likelihoods computed from reads during the run.
- **GL is defined against the cohort-merged allele set.** This set
  depends on which samples were in the run. Re-run with a different
  cohort composition and the allele set changes, so the GL indexing
  changes too.
- **Not reusable.** Feeding freebayes' VCF back to freebayes is not a
  supported workflow. Any re-call requires the BAMs.
- **GL is prior-free.** It is just the data likelihood. The posterior
  (which drove GT and GQ) combined GL with the Ewens / HWE / bias priors
  internally, but those priors don't appear in the GL itself.

### Numerical conventions

Both represent the same mathematical object — the per-sample data
likelihood per genotype — but with different conventions:

| | GATK PL | freebayes GL |
|---|---|---|
| Scale | Phred = `−10 · log₁₀` | log₁₀ |
| Normalisation | best genotype → 0; others are positive Phred penalties | best genotype → 0; others are negative log₁₀ values |
| Conversion | `PL = −10 · GL` (after max-normalisation) | — |
| VCF field name | `PL` | `GL` |
| Ordering | VCF standard: `F(j/k) = k(k+1)/2 + j` | same |

They are the same quantity encoded differently.

### What the different roles mean for cohort recalling

- **GATK.** Add sample 10 001 to a cohort of 10 000: run HC on the new
  sample, then re-run GenotypeGVCFs on all 10 001 gVCFs. No existing
  BAMs are touched.
- **freebayes.** Add sample 10 001: re-run freebayes on all 10 001 BAMs.
  Every previous PL computation is redone.

This is why GATK dominates at scale even when freebayes is technically
better-calibrated at individual sites. PL's durability is the whole
product.

### Implication for this project

The architecture proposed in
[calling_pipeline_architecture.md](calling_pipeline_architecture.md)
takes a different route from both GATK and freebayes: **it stores
neither PL nor GL in the per-sample file, only the five per-allele
observation scalars**.

The rationale is that under freebayes-style likelihood computation, the
PL values are *exactly derivable* from the scalars — they are a
pre-computed projection of the same information, not an independent
quantity. GATK has to store PL because its gVCF is PL-centric: there
is no other stored per-sample evidence. Our `.psf` stores the
observation summary directly, which is strictly more information (it
can be used to compute a likelihood against *any* allele, not just the
sample-local set or the `<NON_REF>` catch-all), so adding PLs on top
would be redundant.

The scalars play the role `<NON_REF>` plays in gVCF — they give the
merger a usable likelihood for any new allele introduced by other
samples — and generalise it: the same scalars support compound
haplotype allele merging, which `<NON_REF>` alone cannot.

In short: **GATK's gVCF stores PL as the per-sample summary;
freebayes' VCF emits GL as an output only; our `.psf` stores the
observation scalars, from which any PL can be reconstructed on demand.**

## Scaling

| | GATK | freebayes |
|---|---|---|
| Per-sample compute | heavy (PairHMM) but trivially parallelisable across samples | n/a |
| Joint call compute | light (EM from PLs) | heavy (combo search over all samples) |
| Joint call memory | O(samples × sites_in_chunk) | O(all reads of all samples in window) |
| Scales to 1 000 samples | yes | no, typically |
| Scales to 10 000 samples | yes | no |
| Re-call after adding a sample | cheap (new HC on that sample + GGVCF only) | expensive (re-run on entire cohort's BAMs) |

This is the single largest reason production pipelines at scale use
GATK's topology even when they'd prefer freebayes' calling behaviour.

## What each does better

### GATK is better at

- **Scaling.** Thousands-of-samples joint calling is routine. Adding a
  sample costs one HaplotypeCaller run, not a full cohort recall.
- **Likelihood accuracy.** PairHMM is a more principled per-read model
  than freebayes' max-of-two-qualities heuristic, especially for indels.
- **Caching.** The gVCF is a durable per-sample artefact. You can
  re-genotype a cohort with different filters or samples without going
  back to BAMs.

### freebayes is better at

- **Compound / haplotype variants.** Because alleles are never
  decomposed, adjacent and overlapping events are called as single
  haplotype alleles. This is where GATK's join-not-merge limitation
  shows.
- **Observation-bias priors in the posterior.** Strand, position, and
  fragment-start biases contribute directly to the posterior rather than
  being post-hoc filters.
- **Calibration at rare variants in small cohorts.** The full Bayesian
  integration over allele frequency is better than EM's point estimate
  when the cohort is small.

## Implications for this project

This project's design
([calling_pipeline_architecture.md](calling_pipeline_architecture.md))
takes **GATK's pipeline shape** (per-sample first, joint second, BAM
read only once, persistent per-sample artefact) but replaces GATK's
per-site join with a **freebayes-style haplotype-aware merge**. In more
detail:

| Design aspect | Choice | Rationale |
|---|---|---|
| Pipeline topology | split per-sample / joint, like GATK | scales to thousands |
| Per-sample artefact | custom `.psf` with the five per-allele scalars freebayes needs | enables freebayes-equivalent likelihoods at joint time without BAM re-read |
| Allele representation in `.psf` | haplotype-level, with anchor positions (Open decision 2) | lets the merger reconstruct compound alleles across samples |
| Merge vs. join | **merge** — compound haplotype alleles assembled across samples | addresses GATK's main limitation |
| Per-sample likelihood model | freebayes-style (max(ln_BQ, ln_MQ) per read, pre-summed) | simpler than PairHMM; sufficient for our target data; flexible upgrade path |
| Joint chicken-and-egg breaker | EM (GATK-style) | linear scaling beats combo-search at 1000+ samples |
| Posterior prior | Dirichlet on allele frequency (GATK-style) | conjugate with EM; well-understood |
| Observation-bias priors | optional future addition, using the per-allele scalars already stored | keeps freebayes' posterior-level bias handling available if needed |
| `<NON_REF>` equivalent | not stored explicitly; the per-allele scalars serve the same purpose | the scalars let the merger compute a likelihood against any new allele, not just a catch-all |
| Per-sample PLs in the file | **not stored** — derivable from the scalars | avoids redundancy; scalars are strictly more informative |
| BAQ application | in-process during Stage 1, parallelised across reads | avoids the `samtools calmd` serial-processing bottleneck |

The result is intended to combine GATK's scalability with freebayes'
haplotype-aware merging — keeping EM because at our target cohort size
the approximation it makes is invisible, and keeping the freebayes-style
per-allele observation summary because it is the smallest sufficient
statistic that supports compound-allele merging without BAM re-read.
