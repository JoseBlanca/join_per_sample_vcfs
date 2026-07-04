# Genotype prior from the site-frequency spectrum

**Status:** draft, 2026-07-04. Model-intent spec — *what* genotype prior the
posterior engine should use, and *why*, to fix a systematic genotyping error at
low coverage in small cohorts. Grown from a GIAB benchmark investigation
([../reports/gt_concordance_vs_giab_2026-07-04.md](../reports/gt_concordance_vs_giab_2026-07-04.md))
and a long design discussion; the empirical prototype lives in
`tmp/prior_experiment.py`. The module/struct shape (the *how*) and the engine
integration are in the companion architecture doc
[../architecture/sfs_genotype_prior.md](../architecture/sfs_genotype_prior.md);
the build order is in
[../implementation_plans/sfs_genotype_prior.md](../implementation_plans/sfs_genotype_prior.md).

**Scope.** This is about the **SNP genotype prior** — how the caller weighs
hom-ref / het / hom-alt *before* reading a site's data. It does **not** address
the indel genotype errors found in the same investigation; those come from
upstream allele mis-assignment in repeats (a separate problem, see the report).

---

## 0. Vocabulary

The genetics you already know; the statistics terms are the ones worth pinning
down, so they get the fuller explanations. Skip whichever group you don't need.

### Genetics / model terms

- **Site-frequency spectrum (SFS).** The distribution of allele frequencies
  across variable sites in a population. Under neutral mutation–drift
  equilibrium most variant alleles are **rare** — the density of a derived
  allele at frequency `x` falls off as `θ/x`.
- **Nucleotide diversity `θ` (≈ `π`).** How much variation a species carries —
  the expected fraction of sites at which two randomly-drawn chromosomes differ.
  Ranges ~100-fold across life: ~0.001 in humans, ~0.02 in maize/Drosophila, up
  to ~0.05 in very diverse species. It is the *scale* of the SFS.
- **HWE — Hardy–Weinberg equilibrium.** The relation between allele frequency
  `p` and genotype frequencies. We use the **inbreeding-adjusted** (Wright) form:
  `P(het) = 2pq(1−F)`, `P(hom-alt) = p² + pqF`, `P(hom-ref) = q² + pqF`.
- **Inbreeding coefficient `F`.** A number in `[0, 1]` for how self-fertilising
  an individual is. `F = 0` outbred; `F → 1` a near-pure selfer. It packs alleles
  into homozygotes, so a selfer rarely carries a variant as a heterozygote. It is
  a *different axis* from `θ`: `θ` is how much variation exists, `F` is how that
  variation is packaged into genotypes.
- **Invariant (non-variant) site.** A genomic position that carries no variation
  in the population — the vast majority of the genome. Distinct from a **hom-ref
  genotype**, which is one individual carrying two reference alleles at a site
  that *may* be variable in others (§4).

### Statistics terms

- **Prior, likelihood, posterior.** For each candidate genotype the caller
  multiplies a **prior** (how probable that genotype is before seeing this site's
  reads) by a **likelihood** (how well that genotype explains the reads), then
  normalises to a **posterior** and takes the winner. The winner's probability
  sets the genotype; how far ahead it is sets the GQ.
- **Plug-in vs marginalize.** Two ways to handle the unknown allele frequency
  `p`. **Plug-in**: estimate a single `p̂`, then apply HWE as if `p̂` were a known
  fact. **Marginalize**: never commit to one `p`; average the genotype
  probability over the whole prior on `p` (weighted by the SFS). The difference
  is the heart of this spec (§3).
- **Empirical Bayes.** Estimating a prior's hyperparameter (here `θ`) from the
  data itself, rather than fixing it by hand.
- **Ewens sampling formula.** The probability of a sample's allelic
  configuration under neutral mutation–drift equilibrium — the SFS in its
  "partition of a sample" form, with the frequencies already integrated out. It
  is what freebayes uses (§3).

---

## 1. The problem, in one sentence

At low coverage in a single sample (or small cohort), the caller systematically
turns a true homozygous-alternate genotype into a heterozygote — `1/1 → 0/1` —
because its genotype prior favours the heterozygote far too strongly.

### The evidence (GIAB per-sample benchmark, HG002/3/4)

Genotype concordance at true-positive SNPs, ours vs freebayes:

| coverage | ours | freebayes | our dominant error |
|---:|---:|---:|---|
| 5× | 83.6% | ~95% | `1/1 → 0/1` (214 of 239), GQ ≈ 5 |
| 10× | 97.8% | 98% | `1/1 → 0/1` |
| ≥15× | ≥99.6% | ≥99% | — (resolved) |

The error is concentrated at 5–10×, is almost entirely `1/1 → 0/1`, and carries
**low GQ** (the caller is not confident — it is being pushed by the prior, not
the reads). It vanishes by 15× once there are enough alt reads to overrule the
prior. A worked example: a true hom-alt SNP sequenced as **0 ref, 2 alt reads**
is called `0/1` with GQ 7; freebayes calls `1/1`.

---

## 2. Root cause — a population prior applied to one sample

The posterior engine's prior is built in two stacked layers
([../reports/gt_concordance_vs_giab_2026-07-04.md](../reports/gt_concordance_vs_giab_2026-07-04.md)
has the full trace; `src/var_calling/posterior_engine.rs`):

1. **Allele-frequency layer (a Dirichlet "pseudocount" prior).** The per-record
   EM estimates the site's alt frequency `p̂` from the read counts **plus** a
   Dirichlet pseudocount of **10 reference / 0.01 alt** (`ref_pseudocount = 10.0`,
   `snp_alt_pseudocount = 0.01`). This encodes "reference is usually the common
   allele."
2. **Genotype layer (HWE).** Given `p̂`, HWE sets the genotype prior:
   `P(het) = 2·p̂·(1−p̂)`, `P(hom-alt) = p̂²`.

For a **single sample**, the frequency estimate is dominated by the pseudocount:
10 phantom reference observations against ~2 real reads pin `p̂ ≈ 0.08`
regardless of the true genotype. (Fingerprint: every het-called site emits the
*identical* `AF = 0.0834 = (1 + α_alt)/(2 + α_ref + α_alt)` — the pseudocount,
not the data, is setting it.) HWE then reads `p̂ ≈ 0.08` as "the alt is rare, so
carriers are heterozygotes," giving `P(het):P(hom-alt) ≈ 22:1`. The two reads
favour hom-alt only ~4:1. **The prior wins; the call flips to het.**

The failure is specific to the small-`n` regime. For a large cohort the
frequency estimate is real external information (many *other* samples inform
`p`), the pseudocount is negligible, and conditional HWE converges to the truth.

---

## 3. Why freebayes does not have this skew — and what the right answer is

Freebayes never estimates a per-site frequency and never applies `p²`/`2pq`. It
scores each genotype directly with the **Ewens** polymorphism prior
(`freebayes/src/Ewens.cpp`, default `θ = 0.01`). For one diploid sample the
Ewens prior ratio is:

| genotype | Ewens weight |
|---|---|
| `0/0` (monomorphic) | `1/(θ+1)` |
| `1/1` (monomorphic) | `1/(θ+1)` — **equal to `0/0`** |
| `0/1` (segregating) | `θ/(θ+1)` |

So freebayes *penalises* the heterozygote by `θ ≈ 0.01` (a segregating site
needs a mutation, which is rare) and treats the two homozygotes **symmetrically**
— there is no reference pseudocount and no `2pq` term. It has nothing that pulls
toward het, so at low depth it follows the reads to hom-alt. (Confirmed on
freebayes' own numbers: its likelihoods favour hom-alt ~4:1, same as ours, yet
it calls `1/1`.)

But freebayes over-corrects: its `θ:1` (≈ 1:100) leans *anti*-het, and at low
depth it over-calls homozygous (its dominant 5× SNP error is the mirror image,
`0/1 → 1/1`). The two callers have **opposite** low-depth biases.

**The theoretically-correct answer sits between them.** Your genetic intuition —
SFS from mutation–drift equilibrium, then HWE+F — is right, *provided you
marginalize over the frequency instead of plugging in a point estimate*.
Integrating the genotype over the SFS density `θ/x`:

- `P(het) = ∫ 2x(1−x)·(θ/x) dx ∝ θ`
- `P(hom-alt) = ∫ x²·(θ/x) dx ∝ θ/2`

gives **het : hom-alt = 2 : 1** — mildly pro-het, as the biology says, but a
factor of ~10 less than our plug-in produces. Against a 4:1 read likelihood, a
2:1 prior still yields a (low-GQ) hom-alt call — the correct answer.

| prior | het : hom-alt | verdict at 2 alt reads |
|---|---|---|
| ours (plug-in HWE, α_ref=10) | ≈ 22 : 1 | het (**wrong**) |
| SFS-marginalized HWE+F | ≈ 2 : 1 | hom-alt (correct, low GQ) |
| freebayes (Ewens θ) | ≈ 1 : 100 | hom-alt (correct, over-confident) |

**The fix is a plug-in→marginalize change, not a pseudocount tweak.** Our error
is conditioning the genotype on a fabricated `p̂`; marginalizing over the SFS
removes the fabrication.

---

## 4. The model — three levels

The prior is a hierarchy. Naming each level, and what is *estimated* versus
*assumed*, is the whole design:

```
θ  (species diversity)         ── top:   how much variation exists
│                                         estimate from data; prior as fallback
▼
p  (site allele frequency)     ── middle: a mixture prior set by θ —
│                                         point mass at p=0 (invariant site)
│                                         + SFS density θ/x (segregating)
▼
genotype                       ── bottom: HWE+F(p) per sample
```

### 4a. Genotype level — HWE+F, marginalized over `p`

Each sample's genotype prior is Wright HWE+F evaluated at `p`, but **averaged
over the posterior on `p`**, not conditioned on a point estimate:

```
prior(genotype_s) = ∫ HWE_F(genotype_s | p) · posterior(p | cohort reads) dp
```

The two regimes fall out of this one formula:

- **Single sample / small cohort** — the reads barely inform `p`, so
  `posterior(p)` ≈ the SFS prior, and the integral collapses to the fixed
  **2:1** marginal (§3). This is the case the prototype validated.
- **Large cohort** — `posterior(p)` is sharp (data-dominated), the integral
  ≈ `HWE_F(genotype | p̂)`, and we recover today's behaviour where it is already
  correct. No discontinuity between the regimes.

The **het:hom-alt ratio is diversity-independent** (2:1 at every `θ`), so this
level fixes the low-coverage error *regardless of the species*.

### 4b. Site-frequency level — the invariant mass is the false-positive knob

The prior on `p` is a mixture: a point mass at `p = 0` (the site is invariant)
plus the SFS density `θ/x` (the site segregates). Empirically (prototype sweep),
the invariant mass **does not touch** the het/hom-alt decision — that is set by
the segregating part's 2:1 — but it **controls how readily a weak-evidence site
is called a variant at all**. It is the precision/recall knob, cleanly separable
from the genotype decision.

Crucially, the invariant mass belongs at the **site** level, integrated **once**
using all samples — **not** re-added to every sample's genotype prior. Conflating
"hom-ref genotype" with "invariant site" is harmless for one sample (they are the
same event) but double-counts across a cohort (a site can be variant in others
while this sample is hom-ref). The engine already has a site-level home for this
quantity — the site QUAL is the Phred of `Π_s P(hom-ref)_s`, literally "is this a
variant at all."

### 4c. Diversity level — `θ` is estimated, with a species-range prior as fallback

What sizes the invariant mass is the species' diversity `θ`, and `θ` swings the
variant-vs-invariant odds ~100-fold across life:

| species | `θ` | P(invariant) | variant : invariant odds |
|---|---:|---:|---:|
| selfing crop line | 0.0005 | 0.9993 | 1 : 1332 |
| human | 0.001 | 0.9985 | 1 : 666 |
| maize / Drosophila | 0.02 | 0.970 | 1 : 32 |
| very diverse | 0.05 | 0.925 | 1 : 12 |

A single fixed `θ` (freebayes' 0.01) is ~10× too permissive for humans and
several-fold too strict for maize. The right answer is species-aware. And `θ` is
one of the most directly **estimable** quantities in population genetics — it
*is* the data — so we estimate it (§5) rather than guess, and fall back to a
weakly-informative prior over the species range only when the data are too thin.

---

## 5. Estimating `θ` before the calling phase

`θ` must be in hand *before* the posterior engine runs (it is the engine's
hyperprior), so it cannot come from the calling EM. The per-sample pileup already
records the sufficient statistics, in the `.psp` metadata
(`src/sample_summary/`): **`callable_positions`** and the rough single-sample
genotype counts **`n_het_sites` / `n_hom_alt_sites`**.

**Estimator — count allele copies, not genotypes (F-free):**

```
θ̂  =  (n_het + 2·n_hom_alt) / (2 · callable_positions)          per sample,
                                                                averaged over the cohort
```

This is F-independent by construction — you count alt-allele copies (het = 1,
hom-alt = 2) over chromosomes surveyed, so the individual's mating system drops
out (the reason `π`/Watterson estimators are F-robust while raw heterozygosity
`n_het/callable = θ(1−F)` is not). The same two counts also yield **F per
sample** via `F = (2−R)/(2+R)`, `R = n_het/n_hom_alt`.

**Accuracy.** The counts come from *rough* single-sample genotyping, so `θ̂`
inherits two biases: low coverage misses variants (biases `θ̂` down), sequencing
error invents hets (biases up). But `θ` only needs to be order-of-magnitude
right — it calibrates a boundary, it does not decide clear calls — so a rough
estimate in the right decade is enough. The coverage-by-GC histogram (also
already stored) is the covariate for an optional low-depth detection-bias
correction (a later tier).

**Architecture win:** the estimate is arithmetic on statistics the var-calling
phase already reads (the hidden-paralog filter opens the same summaries), so it
needs **no new pileup work and no new pass** — a one-shot pre-estimate handed to
the engine. Kept one-shot (not iterated into calling) on purpose, to keep the
coupling and the pipeline structure simple.

---

## 6. What is decided, and what is open

**Decided (validated on the GIAB data — prototype `tmp/prior_experiment.py`):**

- Replace the plug-in HWE prior with an **SFS-marginalized HWE+F** genotype
  prior. On GIAB it lifts 5× SNP concordance **83.6% → 94.5%** (freebayes range),
  cuts `1/1→0/1` from 214 → 5, with **no new hom-ref damage** and no cost at
  higher depth. Robust across error-rate and `θ` (the 2:1 ratio is
  `θ`-independent).
- The **2:1 genotype ratio and the `θ` invariant mass are decoupled knobs**: the
  genotype fix is locked in regardless of `θ`; `θ` only moves precision/recall.
- Estimate `θ` (F-free copy-counting) from the `.psp` summaries at var-calling
  start; species-range prior as the thin-data fallback; `F` on its own axis.

**Open (for the architecture and implementation to settle):**

- Exact vs approximate marginalization over `p` inside the EM hot loop
  (performance) — a `p`-grid integral (as the paralog H1 hypothesis already does)
  versus a closed-form small-`n` special case.
- Whether the SFS prior replaces or wraps the existing per-site EM
  frequency estimate.
- Default on/off and roll-out: this **changes calls by design** (not
  byte-identical), so it lands behind a flag with the full precision/recall panel
  as the gate before it becomes default.
- The low-depth detection-bias correction for `θ̂` (optional second tier), and
  whether a small per-sample alt-fraction histogram in pileup is worth adding for
  an F-robust, coverage-robust frequency-based `θ`.
- **Indels are out of scope** — a separate upstream allele-assignment fix.
