# Architecture — SFS genotype prior + diversity estimation

**Status:** design proposal, 2026-07-04. Turns the model in
[../specs/sfs_genotype_prior.md](../specs/sfs_genotype_prior.md) into a concrete
shape: which modules compute what, where each piece plugs into the existing
posterior engine and pileup summaries, and the CLI surface. Build order is in
[../implementation_plans/sfs_genotype_prior.md](../implementation_plans/sfs_genotype_prior.md).

This is a *proposal*, not a settled architecture — §7 lists the open decisions
that the first implementation steps will resolve. Nothing here is built yet.

---

## 0. Why (plain English, no formulas yet)

The caller currently decides a genotype by estimating how common the alt allele
is at the site, then asking Hardy–Weinberg "given that frequency, which genotype
is expected." For a single sample that frequency estimate is fiction — it is set
by a built-in "reference is common" pseudocount, not the reads — and the fiction
makes the caller expect heterozygotes, so at low depth it turns true
homozygous-alt sites into hets ([../specs/sfs_genotype_prior.md](../specs/sfs_genotype_prior.md) §2).

The fix is to stop committing to a single frequency. Instead the caller averages
the genotype over the whole population picture — the site-frequency spectrum,
scaled by the species' diversity — so no fabricated point estimate can push it.
Three quantities carry this: the species **diversity** (estimated once, before
calling), a **frequency prior** built from it, and the per-sample **genotype
prior** that reads out of it. This doc places each of the three in the codebase.

---

## 1. Glossary (software terms; the model terms are in the spec)

- **Posterior engine.** Stage 6, `src/var_calling/posterior_engine.rs` — the
  per-record EM that produces per-sample genotype posteriors, `best_genotype`,
  `gq_phred`, and the site `qual_phred`.
- **`.psp` sample summary.** `src/sample_summary/` — the per-sample metadata
  block the pileup phase writes into each `.psp` (coverage-by-GC histogram +
  rough heterozygosity counts). Read by var-calling at start-up.
- **Paralog pre-pass.** `src/var_calling/paralog_filter/prepass.rs` — the
  existing cohort-level pass that already opens every `.psp`'s sample summary
  before the main calling loop. Our `θ` estimate piggybacks on this read.
- **Folded-SFS marginalization.** The grid integral over allele frequency the
  hidden-paralog `H1` hypothesis already performs
  (`src/paralog/`, `score_locus_for_paralogy`) — the reusable primitive for §4.

---

## 2. The three components and where they live

```
 pileup phase (per sample, unchanged)
   └─ .psp SampleSummary { coverage_by_gc.callable_positions,
                           heterozygosity.{n_het_sites, n_hom_alt_sites} }
                                    │
 var-calling start-up               ▼
   ┌──────────────────────────────────────────────────────────┐
   │ (A) DiversityEstimate  ── src/var_calling/diversity.rs     │  NEW
   │     θ̂ = (n_het + 2·n_hom_alt)/(2·callable), avg over cohort │
   │     + per-sample F;  or --diversity override / prior       │
   └──────────────────────────────────────────────────────────┘
                                    │  θ, F
                                    ▼
   ┌──────────────────────────────────────────────────────────┐
   │ (B) SfsGenotypePrior  ── src/var_calling/sfs_prior.rs      │  NEW
   │     frequency mixture: point mass at p=0 (size ← θ)        │
   │       + SFS density θ/x   → marginalized HWE+F genotype     │
   │     reuses the paralog folded-SFS grid primitive           │
   └──────────────────────────────────────────────────────────┘
                                    │  prior tables
                                    ▼
   ┌──────────────────────────────────────────────────────────┐
   │ (C) posterior_engine.rs e_step  ── MODIFIED                │
   │     genotype prior comes from (B) instead of HWE(p̂)+Dirichlet│
   └──────────────────────────────────────────────────────────┘
```

---

## 3. Component A — diversity estimation (`DiversityEstimate`)

**New module** `src/var_calling/diversity.rs`. A small struct plus a pure
constructor from the per-sample summaries the paralog pre-pass already collects.

```rust
/// Cohort diversity and per-sample inbreeding, estimated from the rough
/// single-sample genotype counts in the `.psp` sample summaries. Handed to
/// the posterior engine as the SFS prior's hyperparameters, before calling.
pub struct DiversityEstimate {
    /// Nucleotide diversity θ (≈ π), F-free copy-count estimator averaged
    /// over the cohort. The scale of the site-frequency spectrum.
    pub theta: f64,
    /// Per-sample inbreeding coefficient F, from the het:hom-alt ratio.
    /// `None` when the CLI supplies a cohort-wide F instead.
    pub inbreeding: Option<Vec<f64>>,
    /// How θ was obtained (estimated / CLI-supplied / prior-fallback), for
    /// provenance in the VCF header.
    pub source: DiversitySource,
}

impl DiversityEstimate {
    /// θ̂ = mean_s (n_het_s + 2·n_hom_alt_s) / (2·callable_s).
    /// Falls back to `prior_theta` when the cohort is too thin to estimate
    /// (few callable positions / few variant sites).
    pub fn from_summaries(
        summaries: &[SampleSummary],
        prior_theta: f64,
        cli_override: Option<f64>,
    ) -> Self { /* … */ }
}
```

- **Input:** the `Vec<SampleSummary>` the paralog pre-pass already reads — no new
  I/O, no new pass. When the paralog filter is off, the summaries are still
  present in the `.psp` and read cheaply.
- **F-free by construction** (spec §5): counts alt-allele *copies*, so the
  mating system cancels. Per-sample F comes from the same two counts.
- **Thin-data fallback:** when there are too few callable/variant sites to pin
  `θ` (e.g. a tiny region, or a single very-low-coverage sample), fall back
  toward `prior_theta` (the species-range default, §6). This is where the
  "prior over diversities" lives — a weakly-informative regulariser, not a
  hard-coded constant.
- **Provenance:** the chosen `θ`, its source, and per-sample F go into the VCF
  header (`##diversity=…`), matching how the paralog filter records its `π`.

**Optional tier 2 (deferred):** a low-depth detection-bias correction using the
coverage-by-GC histogram — inflate `θ̂` for hets undetectable at the observed
depths. Kept out of the first cut (spec §5).

---

## 4. Component B — the SFS genotype prior (`SfsGenotypePrior`)

**New module** `src/var_calling/sfs_prior.rs`. Turns `(θ, F)` into the
per-genotype prior the E-step needs, by marginalizing HWE+F over the frequency
mixture. This is the mathematical core.

```rust
/// The genotype prior obtained by averaging Wright HWE+F over an SFS prior on
/// the allele frequency — the fix for the single-sample het over-call.
pub struct SfsGenotypePrior {
    /// Frequency grid + SFS weights θ/x (shared across sites of a ploidy).
    grid: FrequencyGrid,
    /// Invariant-site point mass at p=0, sized by θ. The FP/recall knob;
    /// does not affect the het:hom-alt ratio (spec §4b).
    invariant_mass: f64,
}

impl SfsGenotypePrior {
    /// Per-genotype log-prior for a sample, marginalized over p:
    ///   ∫ HWE_F(g | p) · posterior(p | reads) dp
    /// For the small-cohort regime posterior(p) ≈ the SFS prior and this is a
    /// fixed table; for a large cohort it sharpens toward HWE_F(g | p̂).
    fn log_prior(&self, ploidy: u8, n_alleles: usize, f: f64,
                 p_posterior: &FrequencyPosterior) -> Vec<f64>;
}
```

**Reuse.** The `p`-grid marginalization (weights `θ/x`, Wright HWE+F genotype
probabilities per grid point, log-sum-exp reduction) is exactly the machinery the
hidden-paralog `H1` hypothesis already runs in `src/paralog/`. Component B should
lift that primitive (or a shared version of it) rather than re-derive it.

**Two evaluation modes** (the §7 open decision on exact-vs-approximate):

1. **Small-cohort / marginal** — `posterior(p)` ≈ SFS prior, so the genotype
   prior is a *fixed table* per `(ploidy, n_alleles, F, θ)`. Cheap: computed once
   per record shape, splatted across samples (mirrors the engine's existing
   `fill_log_prior_per_g_homogeneous` fast path).
2. **Large-cohort / data-informed** — fold the cohort read evidence into
   `posterior(p)` before marginalizing, so the prior sharpens to `HWE_F(p̂)`. This
   is where it meets the existing EM `p̂`; §7 decides whether to reuse the EM
   estimate as `posterior(p)`'s centre or to integrate directly on the grid.

---

## 5. Component C — engine integration (the `e_step` change)

The change is localised to how the genotype prior enters the E-step; the
likelihood path, the QUAL path, and the record/VCF shape are untouched.

Today (`posterior_engine.rs`):

- `m_step_p_hat` (line ~2828) estimates `p̂` with the Dirichlet pseudocounts.
- `fill_log_indep_per_g` (~3296) and `fill_log_prior_per_g_homogeneous` (~3320)
  build the per-genotype HWE+F log-prior from `p̂`.
- `e_step` (~2531) adds that log-prior to the per-genotype likelihood.

Proposed:

- The per-genotype **log-prior source** switches from "HWE(p̂) built from the
  Dirichlet-regularised `p̂`" to **`SfsGenotypePrior::log_prior`** (Component B),
  parameterised by the cohort `θ` (Component A) and `F`.
- The **invariant mass** feeds the **site-level QUAL** path (`qual_phred`,
  `Π_s P(hom-ref)_s`) — its natural home (spec §4b) — not the per-sample genotype
  prior. This keeps hom-ref-genotype and invariant-site distinct.
- The Dirichlet `ref_pseudocount = 10` / `snp_alt_pseudocount` knobs are
  **retired** for the SNP genotype prior (their job — "reference is common /
  variants are rare" — moves to `θ` and the invariant mass, where it is
  species-calibrated instead of a fixed 10). The EM frequency estimate may stay
  for AF reporting and the large-cohort mode (§7).

**Selection.** A new prior mode (`--genotype-prior sfs|hwe-dirichlet`, default
`hwe-dirichlet` until validated) gates the switch, so the change is opt-in and
the old path stays bit-identical while the new one is measured.

---

## 6. CLI surface

- `--diversity <θ>` — override the estimated `θ`. Default: **estimated** from the
  summaries (Component A). Accepts a value for users who know their organism.
- `--diversity-prior <θ>` — the species-range fallback used when the cohort is
  too thin to estimate (default a low, human-ish `~1e-3`, documented as
  conservative; a diverse-organism user raises it).
- `--genotype-prior sfs|hwe-dirichlet` — selects the prior (§5). Default
  `hwe-dirichlet` (current behaviour) until the validation gate passes, then
  flips to `sfs`.
- `--inbreeding-coefficient <F>` — unchanged; when supplied it overrides the
  per-sample F estimate (keeps `θ` and `F` on separate axes, spec §4/§5).

All recorded in the VCF header for provenance.

---

## 7. Open design decisions (for the first implementation steps to settle)

1. **Exact vs approximate marginalization in the hot loop.** The `p`-grid
   integral per record has a cost; the small-cohort fixed-table mode (§4 mode 1)
   is cheap but is an approximation for mid-size cohorts. Decide the crossover and
   whether a closed form exists for biallelic-diploid (the dominant case).
2. **Reuse the EM `p̂` or integrate on the grid.** In the large-cohort mode, does
   `posterior(p)` re-use the existing EM estimate as its centre, or does the grid
   marginalization subsume the EM? Affects whether `m_step_p_hat` stays.
3. **Multiallelic and non-diploid shapes.** The prototype validated biallelic
   diploid. Confirm the folded-SFS grid extends cleanly to the engine's general
   `(ploidy, n_alleles)` (the paralog primitive is biallelic-oriented).
4. **Default roll-out.** This changes calls by design (not byte-identical), so it
   ships behind `--genotype-prior sfs` and becomes default only after the
   GIAB precision/recall/FP panel confirms no regression at higher depth and a
   real gain at 5–10× — and after a check on the tomato/plant (inbred, `F>0`)
   cohorts the caller actually targets.
5. **`θ̂` bias correction.** Whether the coverage-histogram detection-bias
   correction (Component A tier 2) is needed, decided by how far the rough `θ̂`
   lands from a cohort frequency-based estimate on real data.

---

## 9. Generalizing to all genotype shapes — the Dirichlet-multinomial prior

**Status:** design settled 2026-07-04 (resolves the §7.3 open item). Implemented
biallelic-diploid first (`SfsGenotypePrior` grid + `e_step_sfs_biallelic`,
GIAB-validated: 5× SNP concordance 83.6→94.6, zero precision/recall cost). This
section is the **next phase** — generalize to arbitrary `(ploidy, n_alleles)` and
**delete the HWE(p̂) + Dirichlet plug-in prior entirely**. Its own implementation
plan: [../implementation_plans/sfs_genotype_prior_generalization.md](../implementation_plans/sfs_genotype_prior_generalization.md).

### 9.1 Why a closed form, not a grid

The biallelic grid (integrate the Wright genotype prior over a `θ/p` frequency
grid) does not generalise: a multiallelic site's frequency is a point on the
`(k−1)`-simplex, and a grid over a simplex explodes with the number of alleles;
polyploidy needs the general Wright formula too. There is a closed form instead.

**Marginalising the multinomial genotype prior over a Dirichlet frequency prior
is the Dirichlet-multinomial.** For a genotype `g` with allele-count vector
`k(g)` (how many of each allele, summing to the ploidy `m`) at concentration
`α = (α_0, …, α_{k−1})` (`α_0` = REF):

```text
log P_random(g) = log_multinomial_coeff(g)
                + Σ_a [ lgamma(α_a + k_a(g)) − lgamma(α_a) ]
                − [ lgamma(Σα + m) − lgamma(Σα) ]      ← genotype-independent, cancels in normalisation
```

`log_multinomial_coeff(g)` and `k(g)` are exactly what the engine's
[`GenotypeShape`] already precomputes (`log_multinomial_coeffs`,
`genotype_allele_counts`) for every `(ploidy, n_alleles)`. So the general prior is
`lgamma` over data the engine already has — fast, exact, any shape.

### 9.2 The SFS enters as the concentration `α`

The site-frequency spectrum *is* a Dirichlet on allele frequency. The mapping,
**settled with the owner 2026-07-04** (the "C + choice-2" thread), is deliberately
simple:

- **`α_ref = 1`** — a fixed constant, not a calibration target. It is the value
  that makes the biallelic het:hom-alt ratio `2·α_ref/(α_alt+1)` approach the
  defensible **2:1** as `α_alt → 0`.
- **`α_alt(a) = θ̂ / (n_alleles − 1)`** — the *measured* cohort diversity θ̂
  ([`DiversityEstimate`]), split evenly across the ALT alleles (total alt
  concentration = θ̂, independent of allele count). **No calibration constant.**

**Why this is exactly the clean population-genetics answer (choice 2).** Marginal­
ising Hardy–Weinberg over the neutral SFS density `θ/x` gives the finite, exact
per-individual marginals `P(het) = ∫2x(1−x)(θ/x)dx = θ`, `P(hom-alt) =
∫x²(θ/x)dx = θ/2`, so the variant fraction is `3θ/2` and the monomorphic weight is
`1 − 3θ/2` — spec §3's derivation, no grid, no cutoff. The Dirichlet-multinomial
with `α_ref = 1, α_alt = θ̂` **reproduces these marginals to first order in θ**
(at θ = 1e-3: hom-ref 0.99850, het 0.000998 ≈ θ, hom-alt 0.000500 ≈ θ/2). So the
`α_ref = 1` reference weight *is* the monomorphic mass of §4b — but only because
we adopt the correct variant fraction. No separate invariant-mass term is needed
at the default.

**The 2:1 ratio stays θ-independent at every realistic diversity**, because
`α_alt = θ̂` is always small (even a very diverse organism at θ = 0.02 gives
`α_alt = 0.02` → ratio 1.96 ≈ 2:1). This is what fixes the low-coverage het
over-call regardless of species (spec §4a). The earlier worry that `α_alt ∝ θ`
would drag the ratio anti-het only applied to a *large* proportionality constant
(the discarded "match the grid's variant fraction" option); with the correct
`α_alt = θ̂` it does not arise.

**This intentionally does not reproduce the old grid's hom-ref weight (0.878).**
That figure was an artifact — the grid summed `θ/x` over 200 points without the
`1/200` spacing factor, inflating its variant fraction ~80× above the correct
`3θ/2`. The grid still validated on GIAB because the genotype prior's variant
fraction does not flow into QUAL/calling (that path uses the Beta pseudocounts,
§5) and, at an already-called-variant site with reads, the 2:1 packaging is what
decides the genotype. The new prior's stronger hom-ref weight (0.9985) is the
genetically-correct genome-wide value; **GIAB is the gate** for whether it is too
conservative at the very-low-coverage tail (a single alt read at a cohort-variant
site).

**Reserved knob (build only if GIAB regresses).** If the pure-marginal variant
fraction proves too conservative at the low-coverage tail, the fix is the "C"
separation held in reserve: an independent monomorphic-mass term on the
all-reference genotype, sized from θ̂, so the hom-ref weight can be lowered
*without* changing `α_alt` — keeping the 2:1 ratio intact. Not built at the
default because `α_ref = 1` already supplies the correct monomorphic mass.

**The Dirichlet pseudocounts become `α`.** `ref_pseudocount` / `snp_alt_pseudocount`
are already a Dirichlet on allele frequency — today the engine plugs in their MAP
estimate (the bug). The generalisation *marginalises* them, with `α_alt` sourced
from the measured θ̂ rather than a fixed `0.01`.

### 9.3 Inbreeding `F` — reuse the Wright mixture already in the engine

Wright's model is `P(g) = (1−F)·P_random(g) + F·P_IBD(g)`, where a fully-IBD
individual is homozygous for a single allele drawn at its marginal frequency.
This is **exactly** the mixture `fill_log_prior_per_g_homogeneous` already
computes (`log_sum_exp_2` of an independent term and an IBD homozygote term). The
generalisation keeps that structure verbatim and only swaps the two inputs:

- independent term: `log P_random(g)` (the Dirichlet-multinomial, §9.1) instead of
  the plug-in `Σ_a k_a(g)·log p̂_a`;
- IBD marginal frequency for allele `a`: `log(α_a / Σα)` instead of `log p̂_a`.

So `F` works cleanly at any ploidy — load-bearing for the inbred plant-cohort
target — with no new IBD machinery.

### 9.4 Engine integration and what gets deleted

- **`fill_log_indep_per_g`**: change the per-genotype independent term from
  `Σ k_a·log_p_effective[a]` to the Dirichlet-multinomial `lgamma` form. The
  Wright-`F` wrapper (`fill_log_prior_per_g_homogeneous` / the heterogeneous
  branch) is unchanged except that `log_p_effective[a]` now holds `log(α_a/Σα)`.
- **Applies to every shape** → the biallelic special case
  (`e_step_sfs_biallelic`, `SfsGenotypePrior` grid, `config.sfs_prior_tables`, the
  `EmContext`/routing added for it) is **removed**; the standard `e_step` /
  `e_step_simd` now carry the SFS prior for all records.
- **Delete the plug-in HWE genotype prior**: the genotype prior no longer reads
  the EM `p̂`. Keep `m_step_p_hat` **only** if AF (`INFO/AF`) reporting still wants
  the frequency estimate; otherwise drop it too. The `ref_pseudocount` /
  `snp_alt_pseudocount` CLI knobs are replaced by the θ-derived `α` (retire or
  repurpose them).
- **Open details to settle in implementation:** (a) **compound alleles** — the
  chain-linked `f̂_C` path substitutes a cohort frequency for `p̂[compound]`; give
  compound alleles their own `α_compound` in the Dirichlet-multinomial (or keep
  the `f̂_C` estimate feeding `α_compound`). (b) **`lgamma`** — Rust `std` has no
  `lgamma`; the crate has `ln_factorial` (integer only). Non-integer `α` needs a
  real `lgamma` (e.g. `libm`/`statrs`, or a vendored Lanczos approximation) — a
  small dependency decision. (c) the exact `α_ref` and `θ→α_alt` scale (§9.2),
  fixed by the GIAB re-validation.

### 9.5 Validation

The biallelic numbers shift slightly (the hom-ref weight is computed differently
from the grid), so the GIAB per-sample panel (SNP GT concordance +
precision/recall/FP) is re-run as the gate — the 2:1 ratio is preserved by
construction, but the 94.6 % needs re-confirming. New synthetic tests cover the
shapes GIAB does not exercise: **multiallelic diploid** and **polyploid**
(biallelic + multiallelic), each checking the Dirichlet-multinomial against a
hand-computed value and the `F`-mixture limit (`F=0` random, `F=1` homozygotes).

---

## 8. What this does not touch

- **Indels.** The indel genotype errors are upstream allele mis-assignment in
  repeats, not a prior problem — separate work (report §"Indels").
- **The likelihood, QUAL formula, contamination path, and VCF/record shapes** —
  unchanged. Only the genotype-prior *source* moves.
- **The pileup phase** — no new statistics required for the first cut; it already
  writes the counts Component A needs.

---

## 10. Milestone 4 — empirical-Bayes large-cohort sharpening (design, 2026-07-05)

**Status:** design **signed off by the owner 2026-07-05**; implementation in
progress. Implements Milestone 4 of
[../implementation_plans/sfs_genotype_prior_generalization.md](../implementation_plans/sfs_genotype_prior_generalization.md).
This section settles the four decisions the plan flags: the `α'` update mechanics,
the EM initialisation/convergence fix, the `m_step` interaction, and
biallelic-vs-general behaviour. **Sign-off decisions:** leave-one-out is the
mechanism (single-sample byte-identical to §9 by construction); the permissive
flat-start is the first cut (tomato2 gates spurious het inflation); the large-N
shared-prior compute shortcut is deferred to the perf review.

### 10.0 The problem, in plain English

The prior we ship today (§9) is *species-only*. It asks one question of every
sample in isolation: "given that this species carries variants at rate θ, and I
know nothing else about this site, what genotype do these reads imply?" For a
single sample that is exactly right — it is what stopped the low-coverage
het over-call. But it throws away a real source of evidence: **the other samples
at the same site.**

Picture five samples, each with shallow coverage — say two reads, one reference
and one alternate. Alone, each one is ambiguous: one alt read is as easily a
sequencing error as a true heterozygous site, and the species prior ("alt alleles
are rare") tips each one to homozygous-reference. But all five showing the *same*
alt allele is not a coincidence you expect from independent errors. Together they
are good evidence that the site really is polymorphic, and each individual is
probably a true heterozygote. The caller should let the cohort talk each sample
out of the conservative call. That is **strength borrowing**, and the species-only
prior cannot do it because it never looks at more than one sample at a time.

The fix keeps the exact same machine (§9) and only changes *what frequency the
prior is built from*: instead of the fixed species rate θ, use the frequency the
**cohort's own reads** imply at this site.

### 10.1 The model — a prior that the cohort updates

The site-frequency spectrum is a Dirichlet *prior* on the site's allele
frequency. The cohort's reads are data. Bayes' rule turns the prior into a
**posterior** frequency, and because the Dirichlet is conjugate to the multinomial
allele-count likelihood, that posterior is again a Dirichlet — with its
concentration bumped by the observed allele copies:

```text
species prior      α_species = alpha_from_diversity(k, θ̂)        =  [1, θ̂/(k−1), …]
cohort update      α'        = α_species + E[cohort allele copies]
genotype prior     log P(g)  = dirichlet_multinomial_log_priors(shape, α')   ← §9.1, α → α'
```

`E[cohort allele copies]` is the posterior-weighted allele-count vector the EM
already computes every iteration — `accumulate_expected_counts` inside
`m_step_p_hat` (`scratch.expected_counts`). So the genotype prior is *the same
Dirichlet-multinomial formula as §9*, evaluated at a concentration that has been
nudged by the cohort. No new prior implementation, no grid.

The interpolation is automatic and is the whole point:

- **Few samples** → `E[counts]` is small → `α' ≈ α_species` → the single-sample
  fix from §9 is untouched.
- **Large cohort agreeing** → `E[counts]` dominates the small `α_species` → the
  prior sharpens toward the cohort's actual frequency → strength borrowing.

`F` (inbreeding) rides on top exactly as in §9.3 — the Wright mixture is
unchanged; it now just reads its independent term and its IBD marginal
`log(α'_a/Σα')` off the *updated* `α'` instead of the fixed one.

### 10.2 Why the naive version fails, and the two things that fix it

Coupling the prior back to `E[counts]` re-introduces two classic pitfalls. Both
have clean, standard resolutions, and together they are the real content of this
milestone.

#### (a) A sample must not vote for its own prior — **leave-one-out**

If sample *s*'s genotype prior is built from an `E[counts]` that *includes sample
s's own reads*, the sample reinforces itself: one alt read raises `α'_alt`, the
raised prior makes the alt read look more like a real het, which raises `α'_alt`
further. A **single** low-coverage sample would talk *itself* into a
heterozygote — precisely the over-call §9 was built to kill. That would regress
GIAB, which is a non-negotiable.

The standard fix is **leave-one-out** (LOO): the prior handed to sample *s*
excludes *s*'s own expected counts.

```text
α'_s = α_species + ( E[total counts]  −  E[own counts of s] )
```

This is not a tuning knob; it is what makes the scheme correct, and it earns two
guarantees for free:

1. **Single-sample runs are byte-identical to §9.** With one sample,
   `E[total] − E[own] = 0`, so `α'_s = α_species` on every iteration. The
   genotype prior, the final posteriors, GQ, and the calls are exactly today's.
   **GIAB — a per-sample benchmark — passes by construction, not by re-measuring.**
   (The §9 parity tests, all single-sample, stay green for the same reason.)
2. **No self-reinforcement at any cohort size.** Each sample is judged against
   what the *rest* of the cohort says, never against an echo of itself.

The cost: the genotype prior is now **per-sample** (each `α'_s` differs), so the
homogeneous-fixation hoist that computes one shared prior row per iteration (H4 in
the perf review) no longer applies for `n_samples > 1`. Each sample pays its own
Dirichlet-multinomial `lgamma` evaluation per iteration. For `n_samples = 1` the
prior collapses to `α_species` and today's fast path is preserved. This is
expected per-iteration `lgamma` work; §10.5 flags it for the perf review.

#### (b) Starting at "everyone hom-ref" traps the cohort — **neutral first step**

The EM finds a *local* optimum. If it starts from the species prior (everyone
looks hom-ref → `E[counts]_alt ≈ 0` → `α'_alt` stays tiny → the prior keeps
everyone hom-ref), the cohort of weak hets is stuck at the wrong answer — the
counts never grow, so the prior never sharpens. Walking through the iterations by
hand confirms it: five samples at one-ref-one-alt started from `α_species` creep
*back* toward hom-ref, never toward het. This is the trap the plan warns about.

The escape is to let the **first** E-step run with a **flat genotype prior**
(likelihood only — no species prior, no cohort prior). The reads, unopposed, place
each sample at the genotype they actually favour; the resulting `E[counts]` is an
honest, cohort-wide read on the frequency; and from iteration 2 onward the
leave-one-out EB prior takes over. Concretely:

```text
iteration 1     genotype prior = flat            (posterior ∝ likelihood)
iteration ≥ 2   genotype prior = LOO EB prior    α'_s = α_species + (E[total] − E[own_s])
convergence     unchanged — still on p̂            (m_step / pseudocounts untouched)
```

The flat first step is neutral, not pro-het: a genuinely monomorphic site has
every sample's likelihood pinned at hom-ref, so `E[counts]_alt` comes out ≈ 0 on
iteration 1 and the site stays monomorphic. It only lets *ambiguous* sites float
up to where the reads point before the prior weighs in.

Crucially, **leave-one-out still protects a single sample even with the flat
start.** Iteration 1 may float a lone weak sample toward het, but with one sample
`E[total] − E[own] = 0`, so iteration 2 snaps it straight back to `α_species` →
hom-ref. The flat start changes the *cohort* trajectory without loosening the
single-sample guarantee.

### 10.3 What this does to small cohorts — the one behaviour to watch

The flat first step sets an implicit *threshold*: how many corroborating samples
it takes to overturn the species prior and call het. Hand-tracing the two-read
`(1 ref, 1 alt)` case: one sample → hom-ref (guaranteed), and a handful of
agreeing samples is enough to reach het. A strict marginal-likelihood accounting
(pay the SFS "variants are rare" penalty once, earn the per-sample likelihood gain
N times) puts the crossover at a small N that grows as the per-sample evidence
weakens. The flat-start + LOO scheme lands on the *permissive* side of that
crossover — it will call het for small agreeing cohorts.

This is a **model choice, not a bug**, and it is the one place the owner should
weigh in:

- **Permissive (this proposal).** Two-plus independent samples showing the same
  allele is real corroboration; calling het is defensible and is exactly the
  joint-calling benefit. Simplest; passes the cohort test; GIAB-safe by LOO.
- **Stricter.** Raise the crossover by making the first step *warm* instead of
  flat (a down-weighted species prior on iteration 1), or by a
  multi-start + evidence-selection EM (run from both the hom-ref and the flat
  seed, keep the higher cohort evidence). More faithful to the marginal-likelihood
  optimum; more code and compute.

**Recommendation:** ship the permissive flat-start + LOO scheme first. It is the
minimal correct change, byte-identical on GIAB, and passes the cohort unit test.
Use the **tomato cohort** (real, inbred, `F > 0`) as the accuracy gate: selfing
already damps heterozygosity, and the callset-coherence analysis will show whether
the small-cohort threshold produces *spurious het inflation*. If it does, escalate
to the warm-start or multi-start variant — the machinery (α' update, LOO) is
identical; only the iteration-1 seed changes.

### 10.4 `m_step`, pseudocounts, and AF — unchanged

`α'` is built from `α_species` (`= [1, θ̂/(k−1), …]`) **plus the EM's real
expected counts** — *not* from the `ref = 10` / SNP / indel Dirichlet
pseudocounts. Those pseudocounts keep their existing job on a separate axis:

- `m_step_p_hat` still forms `p̂` with the pseudocounts → drives `INFO/AF`.
- the `α_ref` / `α_alt` Beta pair still feeds the exact-AF QUAL marginalisation.

Neither path reads the genotype prior, so neither changes. The genotype prior and
the AF/QUAL estimate stay decoupled — the genotype prior marginalises the cohort
frequency *posterior*; AF/QUAL keep their pseudocount-regularised point estimate.
This preserves the §9 separation and means the QUAL/precision/recall profile
GIAB validated is untouched at `n = 1`.

### 10.5 Shapes, correctness, and the implementation surface

- **Any `(ploidy, n_alleles)`.** `E[counts]`, `α_species`, and the DM primitive
  are all defined per allele over the general genotype shape; the LOO subtraction
  is per allele. Biallelic diploid is just the `k = 2`, `ploidy = 2` case. The new
  synthetic multiallelic + polyploid end-to-end tests (plan G6) exercise the rest.
- **Positivity.** `α'_a` is `α_species_a` (strictly positive — `α_ref = 1`,
  `α_alt ≥ MIN_ALT_CONCENTRATION`) plus a non-negative LOO count, so it stays
  `> 0`; `lgamma` is always finite. (Floating-point subtraction of the own-count
  could dip a hair below the total; clamp `E[total] − E[own]` at 0 to be safe.)
- **Data flow.** The LOO prior for sample *s* in iteration *t* needs (i) the
  cohort total `E[counts]` from iteration *t−1* (already in
  `scratch.expected_counts` after the `m_step`) and (ii) *s*'s own expected counts
  from iteration *t−1*, recomputed on the fly from `scratch.posteriors[s]`
  **before** the E-step overwrites that row. No new per-sample storage.
- **Where the code changes.** `run_em_loop` gains an iteration counter feeding an
  "iteration 1 = flat prior" branch; the per-sample E-step
  (`fill_log_indep_per_g` + the Wright mixture) moves inside the sample loop for
  `n > 1` to consume the per-sample `α'_s` (via `log_p_effective`,
  `lgamma_alpha`, `alpha`, now refreshed per sample rather than once per record).
  `m_step_p_hat` / `accumulate_expected_counts` / the QUAL path are untouched.
- **Perf.** Per-sample, per-iteration `lgamma` for `n > 1`, replacing §9's
  once-per-record fill. Flagged for the perf review; the biallelic-diploid case is
  ~6 `lgamma` per sample per iteration. `n = 1` keeps the fast path.

### 10.6 Validation plan (gate)

1. Flip `#[ignore]` off
   `cohort_evidence_overcomes_rare_allele_prior_when_all_samples_agree` — the
   unit-level proof that five agreeing weak samples now call het (and AF > 0.15).
2. **GIAB** (`run_ours_per_sample.sh`, 5/10/15/30×): byte-identical to §9 by the
   LOO single-sample guarantee — re-run only to confirm the guarantee holds in
   practice (concordance ≈ 94.9 %, no precision/recall/FP move).
3. **tomato2** (59-sample real cohort, inbred): profile-coherence / callset-change
   analysis — cohort evidence sharpens low-depth-consistent sites toward het where
   warranted **without spurious het inflation**; π/profile stays coherent. Selfing
   `F > 0` damps the effect — a sanity gate, not a concordance number.
4. Add the synthetic multiallelic + polyploid end-to-end tests (plan G6).
