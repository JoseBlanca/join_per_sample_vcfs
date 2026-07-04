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
settled with the owner (tie `α_alt` to the estimated diversity):

- **`α_alt(a) ∝ θ`** — the estimated cohort diversity ([`DiversityEstimate`],
  already wired). A rare-alt SFS is a small `α_alt`; a diverse organism gets a
  larger one, so the variant-vs-invariant balance is species-aware by
  construction. For a multiallelic site split `θ` across the ALTs (e.g.
  `α_alt(a) = θ / (n_alleles − 1)`) — a design detail to pin during
  implementation.
- **`α_ref ≈ 1`** — the reference weight (the invariant-site mass of §4b, in
  Dirichlet form).

This reproduces the biallelic win by construction: at diploid-biallelic the
Dirichlet-multinomial het:hom-alt ratio is `2·α_ref / (α_alt + 1) ≈ 2:1` for
`α_ref ≈ 1, α_alt → 0` — the same defensible 2:1, θ-independent in the ratio,
θ-scaled in the variant mass. (The exact `α_ref` that best matches the validated
grid's hom-ref weight is a calibration target for the GIAB re-run.)

**The Dirichlet pseudocounts become `α`.** `ref_pseudocount` / `snp_alt_pseudocount`
are already a Dirichlet on allele frequency — today the engine plugs in their MAP
estimate (the bug). The generalisation *marginalises* them, with `α_alt` sourced
from θ rather than a fixed `0.01`.

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
