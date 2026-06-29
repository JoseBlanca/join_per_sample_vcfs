# Filtering SNP/indel calls caused by hidden (reference-collapsed) paralogs

**Status:** draft, 2026-06-28. Model-intent spec — *what* we compute to remove
false calls produced by duplications that are present in the sequenced genomes but
**collapsed into a single copy in the reference assembly**, and *why*. Grown one
agreed premise at a time during an exploratory study; the empirical work and all
figures live in [`benchmarks/tomato2/`](../../../benchmarks/tomato2/) (59-sample
*S. lycopersicum* WGS sliced to 160×200 kb random regions, ~6× per-sample depth).
The module/struct/signature shape (the *how*) and the production integration are
deferred to a later architecture doc; §8 sketches the intended architecture.

---

## 1. The problem, in one statement

> A genomic segment duplicated in the sample but **single-copy in the reference**
> collects reads from all its copies onto the one reference locus. Where the copies
> differ, this manufactures **artefactual heterozygotes** (and excess coverage). We
> must remove these calls **without removing introgressions** — divergent alleles
> that are genuinely single-copy.

The two look superficially identical at the variant level (a divergent ALT allele),
so the whole problem is finding a signal that separates them.

| | hidden paralog (REMOVE) | introgression (KEEP) |
|---|---|---|
| copy number at locus | **>1** (collapsed copies) | 1 (orthologous) |
| coverage | **excess** | normal |
| ALT allele | from a divergent collapsed copy | a real, divergent allele |
| genotype | pseudo-het (forced) | real (can be hom-alt) |

## 2. Core principle — coverage is the introgression-safe discriminator

**Only coverage separates the two classes**, because it is the one axis on which
they differ (extra copies → extra reads). Everything else (divergence, heterozygote
excess, allele balance) is *shared* by both classes and therefore cannot be a
primary signal. Empirically (tomato2, [`second_look.py`](../../../benchmarks/tomato2/src/second_look.py)):
among high-het loci, the normal-coverage ("introgression-like") class is **more**
divergent (mean MQDiff −12.5) than the excess-coverage paralog class (−7.3). So a
divergence-*only* filter would be worse than useless: at any MQDiff threshold a
**larger fraction of introgressions clears it than paralogs**, i.e. it would delete
introgressions at a *higher rate* than the very paralogs it is meant to remove. Coverage-based
measures are introgression-safe **by construction**: an introgression is depth-normal,
so any statistic that only fires on excess coverage leaves it untouched.

Two corollaries that shaped everything below:
- **Heterozygote excess must stay diagnostic / secondary.** Legitimate high-het
  loci exist (balancing selection, S-locus, recent introgression). Het alone is not
  a paralog signal; it is a tell-tale only *in combination with coverage excess*.
- **Reference-based mappability is the wrong tool here.** A collapse is *unique in
  the reference* (the duplicate isn't assembled), so k-mer mappability scores it as
  fully mappable — it is blind to exactly the loci we hunt. Mappability is also
  confounded with the CNV signal itself (normalising it away erases real CNVs;
  Janevski 2012, GENSENG). We do **not** use it. (Low-MAPQ reads are already
  filtered upstream, which would also make a mappability covariate inconsistent.)

## 3. The signals and what each can / cannot do

- **Coverage excess** — primary, introgression-safe. Must be measured **per sample,
  GC-corrected, in windows** (§4); per-base coverage at ~6× has no power to resolve
  a 2× carrier from noise.
- **Per-sample observed heterozygosity** — the autogamy/inbreeding level of each
  individual. Used inside the model (§5): in a selfer a real variant is almost
  always homozygous, so a het carrier is anomalous *unless that individual is itself
  heterozygous overall*. We track the **observed** per-sample het rate directly
  (not the derived inbreeding coefficient F = 1 − Hobs/Hexp, though they are
  equivalent); it must be **per sample** because a panel mixes near-homozygous
  inbreds with more-heterozygous wild relatives.
- **Divergence (MQDiff / MQDiffT)** — ALT-vs-REF mean mapping quality. Negative ⇒
  ALT reads map worse ⇒ divergent. Shared with introgression, so **only usable
  coverage-gated** (§5.4). **SNP-only**: indel-carrying reads get low MAPQ from the
  aligner regardless of paralogy, so MQDiff is unreliable for indels and must be
  disabled there.
- **Allele balance + hom-alt** — a collapse keeps the 2 orthologous REF alleles, so
  its VAF is capped below 1; a **confident hom-alt** sample (VAF≈1) is therefore
  near-impossible under the paralog hypothesis and is strong evidence for a real
  variant (§5, the hom-alt veto).

## 4. Measuring coverage (the part that needs real work)

Per-base depth at ~6× cannot tell a 2× carrier (≈12 reads) from an up-fluctuating
single-copy sample. So coverage is taken **per sample, in fixed windows, GC-corrected**.

1. **Windowed depth.** Per sample, mean read depth in fixed `W`-bp windows from the
   `.psp` (post-filter, deduplicated, mate-overlap-resolved fragment coverage — the
   *right* depth for copy number; better than raw BAM depth).
2. **Per-sample GC correction.** Fit each sample's own depth-vs-GC curve (robust /
   binned-median, so paralog windows don't inflate it) and normalise:
   `gc_rel = depth ÷ (sample scale × curve(GC))`, with `1.0 = single copy`.
   **Per-sample, not pooled** — GC bias is a PCR/library property; the curves fan
   out even within one project ([`build_gc_normalization.py`](../../../benchmarks/tomato2/src/build_gc_normalization.py):
   pooled leaves a systematic residual, per-sample flattens it).
3. **Window size.** Smaller windows dilute partial duplications less but have a
   noisier single-copy mode; the count-threshold is a stronger noise filter than
   window size, and **W = 500 bp** sits at a good operating point on tomato2 (single-
   copy mode robust-SD ≈ 0.23; 2× cleanly separable per window).

**Refinement owed:** GC bias acts at the *fragment* scale (Benjamini & Speed 2012),
so a ~insert-size GC window predicts coverage better than the analysis window —
likely tightens the single-copy mode further (§9).

## 5. The per-locus statistic — likelihood ratio of two hypotheses

Per locus, per sample `s` we use coverage `c_s` (windowed `gc_rel`) and the SNP
allele reads `(k_s, n_s)` (ALT, total = AD). Compare:

**H1 — real variant (single copy).** Per-sample genotype from inbreeding-adjusted
Hardy–Weinberg using that sample's observed het level; coverage `~ Normal(1, σ₀)`
**independent of genotype**; VAF ∈ {ε, ½, 1−ε}.

**H2 — hidden paralog (copy-number model).** A carrier has total copies `T` and
`m` mutant copies → coverage `T/2`, **VAF = m/T** (so high copy ⇒ low VAF — a single
mutated copy in a tandem array shows as a *low-VAF* alt). Configs capped at `T = 8`
(coverage 4×); per `T` keep only the single-PSV (`m=1`) and balanced (`m≈T/2`) modes
(the crowded intermediate low-VAF states are not separately identifiable). Carrier
dosage and the non-carrier prior also use the per-sample observed het (in a selfer a
carrier is usually homozygous for the dup).

**5.1 Hom-alt veto** — falls out for free: H2's VAF is capped below 1, so a confident
hom-alt sample has ~zero likelihood under H2 → pushes toward H1.

**5.2 Marginal (Laplace), not maximised.** Both hypotheses are **integrated** over
their parameters (allele freq for H1; carrier freq × configs for H2) with flat
priors, not maximised. Maximising rewards H2's extra flexibility and inflates the
apparent paralog fraction; marginalising applies the Occam penalty. (On tomato2 this
moved the estimated paralog fraction 12.3% → 11.1%; the bulk of the rate is genuine
flagging, not flexibility — see §7.)

**5.3 `LR = logP(data|H2) − logP(data|H1)`.** `>0` favours paralog, `<0` favours
real variant — but the decision cut is **not** 0 (§6).

**5.4 Coverage-gated mqdiff term.** Add `+ w·(−min(MQDiff, 0))` (w≈0.25, SNP-only).
A *diverging* paralog that survives into the callset has only **partial** coverage
excess (its most-divergent reads are MAPQ-filtered; if all were filtered, no variant
is called) plus a negative MQDiff. An introgression has negative MQDiff but **no**
excess. The term is **coverage-gated by construction**: with no excess, H2's coverage
penalty keeps the locus below the cut no matter how strong the term. Simulation
([`paralog_simulation.py`](../../../benchmarks/tomato2/src/paralog_simulation.py)):
single-carrier diverging-paralog recall **32% → 98%** at w≈0.2–0.3 while planted
introgressions stay **0% flagged** at every weight.

## 6. Calibration — empirical-Bayes posterior and the cut

Because the LR is a (marginal) log-likelihood ratio, the per-locus posterior is

> **P(paralog | data) = σ( LR + log(π/(1−π)) )**

with the paralog fraction **π estimated from the data** by EM over the genome-wide
LRs ([`build_paralog_eb.py`](../../../benchmarks/tomato2/src/build_paralog_eb.py)).
This answers "where is the cut": **not at LR = 0**. Paralogs are rare a priori, so
the 50%-posterior boundary sits at a **positive** LR (`log((1−π)/π)`), and at LR = 0
a locus only carries the base-rate posterior (~π). Set the operating point by **FDR
(q-value)**: high-confidence removal ≈ FDR 1% (posterior ≳ 0.9), aggressive ≈ FDR
10–20% (posterior ≳ 0.3–0.5). Default to high-confidence given the introgression-
preservation goal.

## 7. Key findings / design decisions (the *why*)

- **Coverage is the only introgression-safe discriminator** (§2); confirmed on
  tomato2 where divergence is *more* extreme on the introgression-like class.
- **Power requires windowed, per-sample, GC-corrected coverage** (§4); per-base at
  6× cannot resolve 2×.
- **Per-sample** GC curves *and* observed-het, not cohort-wide — library/PCR and
  inbreeding both vary per individual (tomato2 observed-het-derived F ranged
  0.00–0.97).
- **Segregating vs fixed paralogs.** A *segregating* dup is found by cross-sample
  coverage variation (carriers excess, others normal); a dup *fixed in the cohort*
  has no variation, so it is found by **absolute** coverage (vs the sample's own
  single-copy level) plus the heterozygote-excess / read-ratio signature (HWE has an
  absolute anchor the coverage-vs-cohort comparison lacks). The model uses absolute
  per-sample coverage, so it covers both.
- **Single-sample detection works only with confident coverage excess** — a lone het
  at a confidently 2× sample is flagged (within-sample coverage↔allele coupling); a
  lone het at *normal* coverage is kept (introgression-safe). At marginal (1.2–1.5×)
  excess a singleton is not confidently separable from noise.
- **The diverging-paralog blind spot** (coverage false-negatives whose divergent
  reads are MAPQ-filtered) is the one place coverage alone fails; the coverage-gated
  MQDiff term (§5.4) is the fix.
- **The cut is FDR-calibrated and positive** (§6), never LR = 0.
- **Absolute FDR is not yet anchored** — π≈11% is plausible (duplication-rich,
  autogamous genome) and the flagged set looks paralog-like, but only real ground
  truth (read-level simulation / known dups) can confirm the rate (§9).

## 8. Intended production architecture

This is **not** a VCF-only per-locus filter. The introgression-safe core (coverage)
needs the GC-corrected windowed coverage model, which is **not derivable from the
VCF** (the VCF's per-site `DP` is the weak, noisy signal we abandoned). The honest
shape is **pre-pass + per-locus scoring**:

```
PRE-PASS  (per sample, over .psp + reference FASTA)
    - per-sample observed heterozygosity (overall het rate)
    - depth-vs-GC curve + single-copy coverage scale
    - per-W-bp window depths  ->  gc_rel
        ↓
PER-LOCUS SCORING
    coverage gc_rel (pre-pass)  +  GT/AD/MQDiff (per-locus, from .psp/VCF)
        ->  H1-vs-H2 marginal LR  +  coverage-gated MQDiff  ->  posterior  ->  FDR cut
```

- **The `.psp` is already in the pipeline** (`pileup → .psp → var-calling`) and
  `var-calling` already reads every sample's `.psp`, so the pre-pass adds computation,
  not a new data dependency — the same pre-pass pattern as SSR param estimation and
  [contamination estimation](contamination_estimation.md).
- **Better, fold the coverage model into Stage 1 (pileup).** The pileup already
  walks every position per sample, so emitting per-window depth + the per-sample
  GC curve + observed het there is near-free and removes the separate pre-pass.
- **Per-locus, VCF-derivable** (cheap add-ons): observed het *aggregation* aside,
  the per-locus terms are `MQDiff`, allele balance, hom-alt, genotypes.

## 9. Open items / owed work

- **Ground truth** to anchor the absolute FDR: read-level simulation through the real
  pileup, or known tomato segmental dups. (The within-model
  [`paralog_simulation.py`](../../../benchmarks/tomato2/src/paralog_simulation.py)
  validated *mechanisms*, not the absolute rate.)
- **Fragment-scale GC** (Benjamini & Speed) instead of the analysis-window GC.
- **A per-locus production filter** (coverage gc_rel + observed het + MQDiff)
  calibrated against this LR model as the "teacher", with the introgression-safe vs
  aggressive modes selectable; indel handling drops MQDiff.
- **Indels**: covered by the coverage axis (allele-agnostic); MQDiff disabled.
- Optional: a negative-binomial coverage likelihood (coverage is ~1.8× over-dispersed
  vs Poisson) and the joint latent-factor (gCNV-style) formulation as the principled
  long-term model.

## 10. Empirical testbed (where the evidence lives)

All under [`benchmarks/tomato2/src/`](../../../benchmarks/tomato2/src/):
`build_psp.sh`, `run_varcalling.sh` (max-retention callset); `build_feature_table.py`
(per-site features); `build_window_coverage.sh` + `assemble_window_table.py` +
`build_gc_normalization.py` (windowed GC-corrected coverage); `build_window_locus_features.py`;
`build_paralog_lr.py` (the LR model); `build_paralog_eb.py` (empirical-Bayes);
`paralog_extended_test.py` / `paralog_simulation.py` (synthetic validation);
`dashboard.py` (marimo EDA over everything). The deep-research synthesis on coverage
covariates (mappability confound, fragment-scale GC, PCR vs PCR-free) is summarised
in the project notes.
