# Hidden-paralog score — single-individual reformulation — implementation plan

**Status:** draft, 2026-07-02, branch `tomato2-paralog-filter`. Builds the design
in
[hidden_paralog_single_sample_scoring.md](../architecture/hidden_paralog_single_sample_scoring.md).
Three code changes (drop `min_samples`; `F` from the cohort inbreeding
coefficient; remove `Hexp`) plus a profile-based re-validation. Types-first,
pause between steps, each step reviewed and committed on its own.

The statistic and its two-pass wiring already exist (Milestones Q, S, T). This
plan **subtracts** cohort machinery — it is mostly deletion — and then re-checks
the tomato2 profile, since removing `Hexp` and the gate will move the drop set.

## Non-negotiables

- **No new cohort-wide data quantity.** After this, the per-locus score reads only
  per-sample state (coverage model + the locus data) + a single cohort `F`
  constant; the sole global estimate is π (across loci).
- **`--no-paralog-filter` stays a literal no-op** over the single-pass path
  (unchanged).
- **Re-validation is a gate, not a footnote** (Step D). No truth set → judge on
  profile coherence (drop set keeps the coverage+het-excess signature, freebayes
  still emits it, π sane).

---

## Step A — remove the `min_samples` gate

The likelihood ratio self-gates (arch §3): an under-powered locus collapses to
LR≈0 and is kept, so the count gate is redundant.

- **Types first:** drop `min_samples` from `CalibrationConfig` and delete
  `DEFAULT_MIN_SAMPLES_FOR_SCORE` (`calibrate.rs`); drop the `min_samples`
  parameter from `score_spilled_locus` / `score_joined_locus` / `calibrate` /
  `run_write_pass` and their `pipeline.rs` call sites.
- **Build:** delete the `usable < min_samples` early-return in
  `score_spilled_locus` — a locus with any usable samples is scored; a locus with
  *zero* usable samples still returns `None` (nothing to score), unchanged.
- **Tests:** update the calibrate/write-pass unit tests that asserted the
  `min_samples` skip; add one that a low-usable-sample locus now scores to LR≈0
  and is *not* flagged at the default FDR (the self-gate). The
  `calibrate_separates_paralog_from_normal` fixtures still separate.
- **Deliverable:** the filter scores every biallelic-SNP locus; under-powered ones
  fold in near the null and are kept. Full suite green.

---

## Step B — `F` from the cohort inbreeding coefficient; remove `Hexp`

Replace the per-sample `Hexp`-derived `F` with the caller's single cohort
inbreeding coefficient, and delete the accumulator (arch §5, §8). This is the bulk
of the change and is mostly deletion.

- **Types first:**
  - `ParalogPrePass` / `ParalogSampleModel` (`prepass.rs`): drop `obs_het` and
    `callable_positions`; the pre-pass now holds *only* each sample's
    `SingleCopyCoverageModel` (+ its σ₀ accessor). Delete `HexpAccumulator`
    wholesale, plus `callable_reference` and `window_bp`-adjacent het plumbing
    that only fed `Hexp`.
  - Replace `inbreeding_by_sample(prepass, hexp) -> Vec<f64>` with
    `cohort_inbreeding(n_samples, f) -> Vec<f64>` (= `vec![f; n_samples]`), where
    `f` is the caller's `--inbreeding-coefficient`. Keep the `Vec<f64>` shape so
    `ParalogScorePrecompute::new(params, &inbreeding)` is untouched.
  - `calibrate` / `run_write_pass`: drop the `hexp` parameter; take the cohort
    `f` instead (or the pre-built slice).
- **Build:**
  - `pipeline.rs`: `SpillSink` stops accumulating `Hexp` (it only spills now);
    `SinkOutput::Spill` no longer carries a `HexpAccumulator`; the `SinkOutput`
    enum can collapse toward a plain "spill written" signal. Build the `F` slice
    from `cohort.inbreeding_coefficient` once, up front (it needs no pass).
  - Remove the `Hexp` `eprintln` diagnostics; log the cohort `F` used instead.
- **Tests:** update `prepass.rs` tests (no `obs_het`/`Hexp`/`callable_reference`),
  `calibrate.rs` / `write_pass.rs` (no `hexp` arg — pass a cohort `F`), the
  `pipeline.rs` `SpillSink` test (no `Hexp` accumulation assertion). Add a test
  that the cohort `F` reaches the scorer (a high-`F` cohort makes a balanced het
  slightly more paralog-like than an `F=0` cohort, all else equal).
- **Also:** update `examples/paralog_score_parity.rs` (the R1 harness) to use a
  single cohort `F` rather than its per-sample `Hexp` derivation, so it tracks the
  production model.
- **Deliverable:** no `Hexp`, no `obs_het`, no accumulator; `F` is one up-front
  constant; the two-pass main pass only spills (no accumulation). Full suite green,
  clippy clean.

**Checkpoint (after B):** the score is per-individual + per-locus + cohort-`F` +
across-loci π. Pause.

---

## Step C — re-validation on tomato2 (the gate)

Removing `Hexp` and the gate moves the drop set; confirm the paralog *profile*
survives (no truth set — arch §7).

- **Deliverable (extend the T1 report):**
  - re-run tomato2 (59 samples, summary-bearing `.psp`) filter-on vs
    `--no-paralog-filter`;
  - confirm the drop set still splits on **coverage excess + het excess** (T1's
    mean-DP / het-fraction table) and that **freebayes still emits** the dropped
    loci;
  - report π, the drop count, and a **new-vs-old drop-set diff** (Jaccard) to
    quantify the shift from the `F` change + the gone gate;
  - a **high-coverage single-sample** run (one deep tomato sample) to show the
    coverage-only path fires on clear paralogs and stays quiet elsewhere — the
    n=1 graceful-degradation evidence;
  - **re-pin `--paralog-fdr`'s default** if the profile calls for it.
- *Depends:* A, B. *Source:* arch §6, §7.

**Checkpoint (after C):** the reformulated score reproduces the paralog profile at
cohort scale and does something sensible at n=1. Pause for the owner.

---

## Deferred (not in this plan; flagged in the arch)

- **Calibrate the SFS prior to the cohort's observed spectrum** (arch §4) — an
  optional refinement, cohort-dependent, defaults to the theoretical prior at n=1;
  wire + measure only if the profile asks for it.
- **Soft-flag (INFO annotation) instead of hard-drop at very low n** (arch §6) —
  an operator policy; decide after the n=1 evidence in Step C.

## Risk register

- **Profile shift (Step C).** The new `F` is a weak knob (< 10 % historically),
  but combined with the gone gate the drop set will move; the gate is the
  coverage+het-excess profile, not the raw count.
- **π dilution from ungated loci** (arch §3). Conservative direction (stricter
  cut); watch it, do not re-add a hard gate.
- **The R1 harness** (`paralog_score_parity.rs`) must track the `F` change or its
  parity numbers drift from production — update it in Step B.
