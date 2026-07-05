# SFS genotype prior — implementation plan

**Status:** draft, 2026-07-04. Builds the design in
[../architecture/sfs_genotype_prior.md](../architecture/sfs_genotype_prior.md)
(model intent in [../specs/sfs_genotype_prior.md](../specs/sfs_genotype_prior.md)).
Types-first, pause between steps, each step reviewed and committed on its own.
The prototype (`tmp/prior_experiment.py`) has already validated the *model* on
the GIAB data; this plan lands it in the engine, behind a flag, and re-runs the
benchmark as the gate before it becomes default.

## Non-negotiables

- **Opt-in until validated.** The new prior ships behind
  `--genotype-prior sfs`; the default path stays **bit-identical** to today until
  the validation gate (Milestone 3) passes. No silent behaviour change.
- **`θ` estimated before calling, one-shot.** Component A reads the `.psp`
  summaries the paralog pre-pass already opens — **no new pileup work, no new
  pass** — and hands a fixed `θ` to the engine. It is not iterated back into
  calling.
- **`θ` and `F` on separate axes.** Diversity from allele-copy counts (F-free);
  inbreeding from the het:hom-alt ratio or `--inbreeding-coefficient`. Never
  conflate them.
- **The 2:1 genotype fix and the `θ` invariant mass are decoupled** and stay so
  in code: the genotype prior fixes het:hom-alt; the invariant mass feeds the
  site QUAL / FP calibration only.
- **Indels are out of scope** — a separate upstream fix.

---

## Milestone 0 — model validated (done)

The prototype `tmp/prior_experiment.py` established, on GIAB HG002/3/4 true-
positive SNPs, that the SFS-marginalized HWE+F prior lifts 5× GT concordance
83.6% → 94.5%, cuts `1/1→0/1` from 214 → 5, with no new hom-ref damage and no
higher-depth cost; robust across error-rate and `θ`; the 2:1 ratio is
`θ`-independent. This plan does not re-litigate the model — it productionises it.

---

## Milestone 1 — the two pure primitives (no engine wiring yet)

Types and pure functions first, unit-tested in isolation. Nothing calls them yet.

### Step A — `DiversityEstimate` from the `.psp` summaries

- **Types first:** `DiversityEstimate { theta, inbreeding, source }` +
  `DiversitySource` enum in a new `src/var_calling/diversity.rs`.
- **Build:** `from_summaries(&[SampleSummary], prior_theta, cli_override)`
  implementing `θ̂ = mean_s (n_het_s + 2·n_hom_alt_s)/(2·callable_s)`; per-sample
  `F = (2−R)/(2+R)`, `R = n_het/n_hom_alt`; thin-data fallback toward
  `prior_theta`.
- **Tests:** F-free property (θ̂ invariant when het↔hom-alt shift at fixed copy
  count under changing F); known-value cases (human-ish counts → θ ≈ 1e-3);
  thin-cohort → fallback; zero-callable / zero-variant corners; the F recovery
  from a selfing-vs-outbred count pair.
- **Deliverable:** a pure estimator, no I/O, no engine dependency. Full suite
  green.

### Step B — `SfsGenotypePrior` (marginalized HWE+F)

- **Types first:** `SfsGenotypePrior { grid, invariant_mass }` +
  `FrequencyGrid` in `src/var_calling/sfs_prior.rs`. Decide whether to lift the
  paralog folded-SFS grid primitive into a shared location or wrap it (arch §7.2).
- **Build:** `log_prior(ploidy, n_alleles, F, p_posterior)` — mode-1
  (fixed-table, `posterior(p)` ≈ SFS prior) first; the biallelic-diploid closed
  form if one exists, else the grid integral.
- **Tests:** reproduce the prototype's `[0.827, 0.115, 0.058]` and the **2:1**
  het:hom-alt at F=0; `θ`-independence of the ratio; F pushes mass to
  homozygotes; invariant-mass sweep moves only the hom-ref weight, not the
  het:hom-alt ratio (the decoupling non-negotiable, as a test).
- **Deliverable:** a pure prior table generator matching the validated prototype.
  Full suite green.

**Checkpoint:** both primitives exist and are unit-tested against the prototype
numbers, with no engine change yet. Pause for review.

---

## Milestone 2 — wire into the engine behind a flag

### Step C — read `θ`/`F` at var-calling start-up

- **Build:** in the var-calling driver, after the `.psp` summaries are loaded
  (the paralog pre-pass path), construct `DiversityEstimate`. Thread `θ` and the
  per-sample `F` into `PosteriorEngineConfig` (new fields, validating setters,
  matching the existing config style).
- **Tests:** the summaries are read once; a run with the paralog filter off still
  gets `θ` (summaries present); provenance recorded.
- **Deliverable:** `θ`/`F` available to the engine; no behaviour change yet
  (nothing consumes them).

### Step D — the `e_step` prior switch

- **Types first:** `--genotype-prior {sfs, hwe-dirichlet}` on the CLI →
  `PosteriorEngineConfig` mode field (default `hwe-dirichlet`).
- **Build:** in `e_step` / `fill_log_prior_per_g_*`, when mode is `sfs`, source
  the per-genotype log-prior from `SfsGenotypePrior` (Step B) instead of
  `HWE(p̂)`; route the invariant mass to the QUAL path, not the per-sample prior.
  Leave the `hwe-dirichlet` path untouched.
- **Tests:** `hwe-dirichlet` mode is **byte-identical** to `main` on the existing
  var-calling integration fixtures; `sfs` mode changes low-coverage single-sample
  hom-alt calls in the expected direction on a small synthetic fixture; QUAL still
  well-formed under `sfs`.
- **Deliverable:** the new prior runs end-to-end under `--genotype-prior sfs`; the
  default remains bit-identical. CLI + header provenance wired.

**Checkpoint:** the engine can call with either prior; default unchanged. Pause
for review + a full byte-identity check on the default path.

---

## Milestone 3 — the validation gate

### Step E — GIAB precision/recall/FP + GT-concordance panel

- **Build:** run `--genotype-prior sfs` end-to-end through the GIAB per-sample
  benchmark at 5×/10×/15×/30×/50×/300× (the pipeline in
  `benchmarks/giab/src/run_ours_per_sample.sh`), and re-score the freebayes
  comparison dashboard: the GT-concordance panel **and** precision/recall/FP.
- **Gate (all must hold before `sfs` can become default):**
  - 5–10× SNP GT concordance rises to freebayes range (~94%+ at 5×), matching the
    prototype;
  - **no regression** in GT concordance at ≥15×;
  - **no material precision/recall/FP regression** at any depth (the invariant
    mass / `θ` calibration is the knob to check here — sweep it if needed);
  - a sanity pass on an inbred (`F>0`) cohort (tomato) confirming `θ`/`F`
    estimation behaves and calls do not degrade.
- **Deliverable:** a validation report in `doc/devel/reports/` with the panels,
  the chosen `θ` handling, and a go/no-go on the default flip.

### Step F — flip the default (only if Step E passes)

- **Build:** default `--genotype-prior` to `sfs`; keep `hwe-dirichlet` available.
  Update the benchmark preset docs.
- **Deliverable:** SFS prior is the default SNP genotype prior; old path retained
  behind the flag.

**Checkpoint:** default flipped, or documented reasons it stayed opt-in. Pause.

---

## Deferred (post-gate, only if the data ask for it)

- **`θ̂` detection-bias correction** using the coverage-by-GC histogram
  (architecture Component A tier 2) — if the rough `θ̂` proves too coverage-biased
  on real cohorts.
- **F-robust frequency-based `θ`** — add a small per-sample alt-fraction
  histogram to the pileup summary and estimate `θ` (Watterson/π) in the cohort
  pre-pass, if the copy-count estimator is insufficient for diverse/inbred edges.
- **Large-cohort data-informed marginalization** (architecture §4 mode 2) — fold
  cohort read evidence into `posterior(p)` so the prior sharpens to `HWE(p̂)`;
  only needed if the fixed-table mode mis-serves mid-size cohorts.

---

## Notes on ordering and risk

- **Primitives before wiring** (Milestone 1 before 2): the maths is the risky
  part and it is testable in isolation against the prototype numbers — lock it
  down before touching the engine hot loop.
- **Opt-in before default** (Milestone 2 before 3): the default path stays
  bit-identical, so the change carries zero risk to existing runs until the gate
  passes.
- **The `θ` estimate is cheap and non-invasive** — reuses an existing read, adds
  no pass — so Component A is low-risk despite touching the driver.
- **The indel problem is independent** and can proceed in parallel without
  interaction (different failure mode, different stage).
