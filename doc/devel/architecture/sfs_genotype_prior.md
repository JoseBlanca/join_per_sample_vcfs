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

## 8. What this does not touch

- **Indels.** The indel genotype errors are upstream allele mis-assignment in
  repeats, not a prior problem — separate work (report §"Indels").
- **The likelihood, QUAL formula, contamination path, and VCF/record shapes** —
  unchanged. Only the genotype-prior *source* moves.
- **The pileup phase** — no new statistics required for the first cut; it already
  writes the counts Component A needs.
