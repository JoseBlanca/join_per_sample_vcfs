# SFS genotype prior — generalization to all shapes + HWE removal — implementation plan

**Status:** draft, 2026-07-04. Builds the design in
[../architecture/sfs_genotype_prior.md](../architecture/sfs_genotype_prior.md) §9.
Follows the completed biallelic-diploid phase (Steps A/B0/B/wiring, GIAB-validated:
5× SNP concordance 83.6→94.6, zero precision/recall cost — see
[../reports/sfs_prior_giab_validation_2026-07-04.md](../reports/sfs_prior_giab_validation_2026-07-04.md)).
This phase replaces the grid-based biallelic prior with the **Dirichlet-multinomial**
closed form so the SFS prior covers **every** `(ploidy, n_alleles)`, then
**deletes the HWE(p̂) + Dirichlet plug-in genotype prior entirely**.

Intended for a fresh session (the design one accumulated heavy context). Each step
is types-first, reviewed (implement → review → apply-fixes → commit), and paused
between milestones.

## Non-negotiables

- **The genotype prior no longer reads the EM `p̂`.** It is the Dirichlet-multinomial
  over the genotype shape, `α` from θ. `m_step_p_hat` survives only if `INFO/AF`
  needs it; the prior does not.
- **`F` reuses the existing Wright mixture** (`(1−F)·random + F·IBD`), not new IBD
  machinery — so it works at any ploidy for the inbred-cohort target.
- **`α_alt` is tied to the estimated θ** (owner decision), not a fixed pseudocount.
- **GIAB is the gate before deleting the old path** — re-confirm the SNP win and
  zero precision/recall cost; the 2:1 ratio is preserved by construction but the
  hom-ref weight changes.
- **Indels stay out of scope** (separate AD-driven fix).

---

## Milestone 1 — the general prior primitive (pure, no engine change)

### Step G1 — `lgamma`
- Decide the `lgamma` source (`libm`, `statrs`, or a vendored Lanczos). Add it,
  unit-test against known values. Keep `ln_factorial` for the integer paths.
- **Deliverable:** a tested `lgamma(x)` available to the prior.

### Step G2 — Dirichlet-multinomial genotype log-priors
- **Types first:** a pure function
  `dirichlet_multinomial_log_priors(shape: &GenotypeShape, alpha: &[f64]) -> Vec<f64>`
  (one log-prior per genotype), in `src/genetics.rs` beside the Wright primitives.
- **Build:** `log_multinomial_coeff(g) + Σ_a[lgamma(α_a+k_a) − lgamma(α_a)]` over
  the shape's `genotype_allele_counts`; the `Σα` term cancels (document that it is
  omitted).
- **Tests:** biallelic-diploid reproduces het:hom-alt `2α_ref/(α_alt+1)` and, at
  `α_ref≈1, α_alt→0`, the 2:1 ratio (tie to the Step-B numbers); a hand-computed
  triallelic-diploid case; a polyploid (tetraploid biallelic) case; normalisation.
- **Deliverable:** the general random-mating prior, matching the biallelic result.

### Step G3 — `α` from θ
- **Build:** `alpha_from_diversity(n_alleles, theta, ...) -> Vec<f64>` — `α_ref`
  (reference weight) + `α_alt(a)` from θ (single ALT: `α_alt = f(θ)`; multiallelic:
  split θ across ALTs). Calibrate `α_ref` / the θ scale so diploid-biallelic
  matches the validated grid at the human default.
- **Tests:** biallelic α reproduces the Step-B genotype prior within tolerance;
  θ scales the variant mass; multiallelic split sums correctly.
- **Deliverable:** the θ→α mapping; the biallelic prior reproduced from α.

**Checkpoint:** the general prior + α mapping exist and reproduce the biallelic
win, no engine change yet. Pause for review.

---

## Milestone 2 — wire into the engine for all shapes; delete the biallelic special case

### Step G4 — Wright-`F` mixture over the Dirichlet-multinomial
- **Build:** in `fill_log_indep_per_g`, replace the independent term
  `Σ k_a·log_p_effective[a]` with the Dirichlet-multinomial term (Step G2); set
  `log_p_effective[a] = log(α_a/Σα)` for the IBD term so
  `fill_log_prior_per_g_homogeneous` / the heterogeneous branch are unchanged.
- Thread `α` (from θ + n_alleles) into `EmContext`/`RecordScratch` in place of the
  per-record `p̂`-derived frequencies.
- **Handle compound alleles:** give them an `α_compound` (from the `f̂_C` path or a
  fixed prior — decide here, arch §9.4 detail a).
- **Tests:** the engine's per-genotype prior matches the pure Step-G2/G4 value for
  biallelic, multiallelic, and polyploid records; `F=0`/`F=1` limits.
- **Deliverable:** the SFS prior applies to **all** record shapes via the standard
  `e_step` / `e_step_simd`.

### Step G5 — remove the biallelic special case + the HWE plug-in
- **Delete:** `e_step_sfs_biallelic`, the grid `SfsGenotypePrior` +
  `src/var_calling/sfs_prior.rs` (or reduce to the α helper if reused),
  `config.sfs_prior_tables` + its `EmContext` field + routing, the driver's
  per-sample-table build. Delete the plug-in genotype-prior use of `p̂`; drop
  `m_step_p_hat` unless `INFO/AF` needs it (decide + document). Retire/repurpose the
  `ref_pseudocount`/`snp_alt_pseudocount` CLI knobs.
- **Tests:** the deleted-path tests are removed or repointed; the suite stays green.
- **Deliverable:** one genotype-prior path (Dirichlet-multinomial SFS), old prior gone.

**Checkpoint:** single prior path, HWE deleted. Pause + full-suite + clippy/fmt.

---

## Milestone 3 — validation gate

### Step G6 — GIAB re-validation + shape tests
- Rebuild the **host** binary (`cargo build --release`; the GIAB script execs the
  host build, not the container one — see the validation report). Re-run the GIAB
  per-sample panel at 5×/10×/15×/30×.
- **Gate:** SNP GT concordance still ≈ the biallelic-phase result (≥94 % at 5×,
  beats freebayes at every depth); precision/recall/FP unchanged; no regression at
  higher depth. A tomato (inbred, `F>0`) sanity pass confirms the F-mixture behaves.
- Add the synthetic **multiallelic** and **polyploid** end-to-end tests (shapes
  GIAB does not cover).
- **Deliverable:** a validation report; go/no-go on the deletion (already done in
  G5 — this confirms it).

---

## Notes on ordering and risk

- **Pure primitive before wiring** (M1 before M2): the Dirichlet-multinomial +
  α-from-θ are the risky maths, testable in isolation against the validated
  biallelic numbers before touching the engine core.
- **Wire before delete** (G4 before G5): confirm the general path matches the pure
  values on all shapes, then remove the old path.
- **`lgamma` dependency** (G1) is the one external decision; keep it small.
- **Compound alleles** and **`m_step_p̂`/AF** are the two integration details most
  likely to surprise — decide them explicitly in G4/G5, don't let them ride.
