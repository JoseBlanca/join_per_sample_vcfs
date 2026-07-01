# Hidden-paralog filter — the consumer model + var-calling wiring — implementation plan

**Status:** draft, 2026-07-01, branch `tomato2-paralog-filter`. Turns the
settled consumer architecture
([hidden_paralog_locus_statistic.md](../architecture/hidden_paralog_locus_statistic.md),
Premises 0–6) into build order. Design intent lives in the spec
[hidden_paralog_filter.md](../specs/hidden_paralog_filter.md) (§4–§7); the
validated reference maths live in
[`benchmarks/tomato2/`](../../../benchmarks/tomato2/)
(`build_paralog_lr.py`, `build_paralog_eb.py`, `build_gc_normalization.py`).

This is the **follow-on** to
[paralog_psp_summaries.md](paralog_psp_summaries.md) (the producer side —
coverage-by-GC histogram + het counts in the `.psp`, Milestones A1–D3
shipped). That plan stopped at "the data is produced, stored, and readable";
this one builds the model that consumes it and the var-calling integration
that emits a `FILTER` verdict.

## Domain intent

At each locus the filter asks which of two stories better explains all the
samples' coverage and allele counts: **H1** a real single-copy variant, or
**H2** a hidden (reference-collapsed) paralog. The verdict is a per-locus
**likelihood ratio**, turned into a probability via a data-estimated prior
(empirical Bayes) and cut at a false-discovery rate that preserves
introgressions. Coverage — measured per sample against the sample's own
single-copy level — is the introgression-safe primary signal; the allele
counts and the per-sample inbreeding coefficient are supporting.

## Scope / non-goals

In scope: one producer amendment (the callable-position total), the pure
statistics module `src/paralog/` (coverage-model fit, the H1-vs-H2 marginal
LR, the inbreeding scalar, the empirical-Bayes prior + FDR), a data-first
validation harness, and the var-calling wiring that scores every locus,
calibrates, and writes the `FILTER`/INFO into the VCF.

**Out of scope:** anchoring the *absolute* FDR to ground truth (read-level
simulation / known duplications — spec §9); the cheap teacher-trained
per-locus production filter (spec §9); indels beyond "they use the same
coverage-only score" (MQDiff already disabled there).

## Principles (how the order was chosen)

- **Types first, then implementation, within every step** (project rule);
  pause between milestones (incremental rule).
- **Pure core before wiring.** The three statistics pieces (coverage model,
  scorer, prior/FDR) are pure functions with no cohort-plumbing dependency
  (arch Premise 0), so they are built and unit-tested in isolation, then the
  var-calling wiring depends *inward* on a settled interface.
- **Validate on data, not on the prototype.** The Python prototype is a
  *reference draft* the production model has deliberately moved past (spec §3
  MQDiff dropped, §5.2 SFS prior, §5 het rate, §7 Normal kept). Validation is
  the Rust code on the tomato2 data — π in range, flagged set paralog-like,
  single-copy peak at `relative_copy_number = 1.0` — with a **loose
  Python↔Rust LR correlation** as a porting sanity-check, not a bit-exact gate
  (arch Premise 5).
- **Bounded RAM.** The prior/FDR are global (no streaming form), so the
  wiring spills per-locus records to an **ephemeral** temporary file and
  derives π + the `q(LR)` curve from a fixed-size streaming LR histogram —
  never a genome-wide in-memory vector (arch Premise 6). The temp file is
  scratch, deleted at run end; it is *not* a `.psp`-style durable artefact.

---

## Milestone P — producer amendment: the callable-position total

**P1. Store the per-sample callable-position total.** ☐
Observed heterozygosity is the het **rate** `n_het / callable_positions`
(spec §3; *not* `n_het/(n_het+n_hom_alt)`, which tracks reference divergence
and inverts). The het count is already stored; the denominator is not.
- Types first: add one `u64` grand-total to `CoverageByGcAccumulator` (it
  already counts `covered` per tile — accumulate the total) and a
  `callable_positions` field on the stored summary; extend the TOML payload
  and the reader accessor.
- **Deliverable:** the accumulator sums the total across tiles/regions; the
  field round-trips through the `.psp`; a test asserts it equals the summed
  covered (non-`N`) positions on a small fixture.
- *Depends:* the shipped producer (paralog_psp_summaries.md C2/D1).
  *Source:* arch Premise 3; [`src/sample_summary/coverage.rs`](../../../src/sample_summary/coverage.rs).

---

## Milestone Q — the pure statistics module `src/paralog/`

New top-level module, sibling to `src/sample_summary/`, depending on
`sample_summary` and nothing in `var_calling` (arch Premise 0). Layout:
`mod.rs`, `coverage_model.rs`, `locus_score.rs`, `prior.rs`.

**Q1. Module skeleton + `ParalogModelParams`.** ☐
Types first: the module tree and the fixed model constants/grids struct with
defaults — `pseudocount_vaf`, `max_relative_copy_number`, `carrier_copy_numbers`
`{3,4,6,8}`, `allele_freq_prior` (folded SFS `1/(p(1−p))` on `[1/2N, 1−1/2N]`),
`carrier_freq_grid` (40 pts, `0.004..0.6`, flat), the hom-alt-veto thresholds.
No divergence/MQDiff field (dropped — spec §3).
- **Deliverable:** the module compiles; `ParalogModelParams::default()` +
  a test pinning the default values. *Source:* arch Premise 2.

**Q2. `SingleCopyCoverageModel::fit`.** ☐
Types first: the model (`single_copy_scale`, `gc_bias_curve`,
`single_copy_depth_sd`, `window_bp`), `CoverageFitConfig`, `CoverageModelError`.
Then the fit from the stored `CoverageByGcHistogram`:
- `single_copy_scale` = the **global mode** of the depth distribution (not the
  first/lowest peak); **reject** (`CoverageModelError`) if `mode/median ∉
  [0.5, 1.5]` (the mode landed on a wrong copy-number peak — arch Premise 1,
  F-9 guard). Keep `depth_bin_width` a small fraction of the single-copy level.
- `gc_bias_curve` = per-GC-bin weighted median of `rel` over the single-copy
  band `0.4 < rel < 1.6`, gap-filled + smoothed.
- `single_copy_depth_sd` (σ₀) = `1.4826·MAD` of single-copy `rel` — **fit from
  the data**, not hardcoded (≈ 0.28 on tomato2).
- `relative_copy_number(gc, depth) = depth / (scale · curve(gc))`.
- **Deliverable:** unit tests (synthetic histogram → known scale/curve;
  wrong-peak histogram rejected; N-only/degenerate rejected). Data check in R1.
- *Depends:* Q1, P1. *Source:* arch Premise 1; spec §4.

**Q3. The scoring function `score_locus_for_paralogy`.** ☐
Types first: `SampleObservation` `(relative_copy_number, alt_reads,
total_reads, inbreeding_coefficient)`, `LocusObservations` (samples only — no
divergence field), `ParalogScore`. Then the **pure** LR (arch Premise 2, spec
§5), faithful to `build_paralog_lr.py`:
- **H1:** marginalise allele freq `p` under the folded-SFS prior; per-sample
  Wright HWE genotype prior with `F`; coverage `~ Normal(1, σ₀)` independent of
  genotype; binomial allele factor. Marginal by log-sum-exp over the grid.
- **H2:** per carrier config `(T, m)` (single-PSV `m=1` and balanced `m≈T/2`),
  coverage `~ Normal(T/2, σ₀·√(T/2))`, `vaf = m/T`; non-carrier single-copy;
  Wright carrier prior in `q`; marginal over `config × q`.
- **Result:** `paralog_log_likelihood_ratio = logL_H2 − logL_H1`, a *pure*
  likelihood ratio (nothing added). Report `confident_homalt_carriers` (the
  hom-alt veto falls out of H2's capped VAF).
- **Deliverable:** unit tests — a confident hom-alt sample forces negative LR
  (veto); excess-coverage carriers force positive; normal-coverage real variant
  stays negative regardless of allele pattern (introgression-safe); a fixed set
  of `(rel, alt, total, F, σ₀)` reproduces a hand-checked LR.
- *Depends:* Q1. *Source:* arch Premise 2; spec §5.

**Q4. The inbreeding scalar `F` from obs het + `Hexp`.** ☐
Types first: `obs_het(het_counts, callable_positions) = n_het/callable` (pure,
single-sample) and `inbreeding_coefficient(obs_het, hexp) = clip(1 −
obs_het/hexp, 0, 0.99)` (pure, given the cohort scalar `Hexp`). `Hexp` itself
is accumulated by the wiring (S1), not here.
- **Deliverable:** unit tests on the two pure functions; the data check
  (Spearman ≈ 0.86 vs the prototype's cohort F) lands in R1.
- *Depends:* Q1, P1. *Source:* arch Premise 3; spec §3.

**Q5. `ParalogPrior` — empirical-Bayes prior + FDR.** ☐
Types first: `ParalogPrior { prior_probability }`, `EmConfig`, and an **LR
histogram** accumulator (fixed bins) so the estimate never needs the per-locus
vector. Then (faithful to `build_paralog_eb.py`):
- `estimate`: EM on the histogram — `π ← Σ_bin w_bin·sigmoid(LR_bin + logit π)`,
  iterate to `|Δ| < tol` (concave → unique fixed point, start-independent).
- `paralog_posterior(LR) = sigmoid(LR + logit π)`; `half_posterior_ratio`.
- `q_of_lr`: the monotone tail-FDR-as-a-function-of-LR curve, computed from the
  histogram (so a locus's q-value is `q_of_lr(LR_i)` — no global sort).
- **Deliverable:** unit tests — EM recovers π on a synthetic two-component
  mixture; `q_of_lr` monotone decreasing in LR; histogram-derived π/q match a
  brute-force full-vector computation within bin tolerance.
- *Depends:* Q1. *Source:* arch Premise 4; spec §6.

---

## Milestone R — data-first validation harness

**R1. `examples/paralog_score_parity.rs`.** ☐
A runnable example (mirroring `examples/dump_sample_summary.rs`) that runs
Q2–Q5 on the tomato2 `.psp` summaries + cohort VCF and reports the
data-first checks:
- coverage model: the single-copy peak sits at `relative_copy_number = 1.0`;
  σ₀ ≈ 0.28; the mode/median guard passes on real samples;
- `F`: Spearman ≈ 0.86 vs the prototype's cohort F;
- scoring/prior: π ≈ 9% (folded-SFS), the FDR-flagged set carries the paralog
  profile (elevated het + coverage excess vs the kept set);
- **loose Python↔Rust LR correlation** ≥ (say) 0.98 on shared loci — a porting
  sanity-check, not an exact-parity assertion.
- **Deliverable:** the example + a short note in
  `doc/devel/reports/implementations/`; any correlation shortfall traced to a real
  bug, not tolerated as "close enough."
- *Depends:* Q2–Q5. *Source:* arch Premise 5; spec §7, §10.

---

## Milestone S — var-calling wiring (bounded-RAM, spill-to-disk)

Wire into `var_calling::pipeline::run_var_calling`. The data-flow architecture
is settled in
[hidden_paralog_varcalling_wiring.md](../architecture/hidden_paralog_varcalling_wiring.md)
(supersedes the Premise 6 sketch): the pass structure + `Hexp` timing, the
per-window coverage source (a per-sample body pass), the ephemeral spill format
(reuse the `.psp` block-writer), and the owner-settled decisions — **hard
removal** of flagged loci and **filter on by default** (`--paralog-fdr ≈ 0.01`,
`--no-paralog-filter` to disable).

**S1. Pre-pass: per-sample models + cohort `Hexp`.** ☐
Where the `PspReader`s are open: read each sample's summary → fit
`SingleCopyCoverageModel`; form `obs_het`. Accumulate the cohort
`Hexp = mean 2pq` **inside the existing per-locus genotype pass** (not a
separate pass); then `F_s = 1 − obs_het_s/Hexp`. Assemble a `ParalogPrePass`.
- **Deliverable:** the pre-pass value built for the tomato2 cohort; per-sample
  models + F available; a test that `Hexp` matches an independent sum of `2pq`.
- *Depends:* Q2, Q4. *Source:* arch Premise 6 step 1.

**S2. Pass 1 — score every locus, spill + histogram.** ☐
After each `PosteriorRecord` is built: measure the locus tile's `(gc,
mean_depth)` from the `.psp` body → `relative_copy_number` via the per-sample
model; build `LocusObservations`; call `score_locus_for_paralogy`. Write
`(position, LR, + the fields the final VCF needs)` to an **ephemeral temporary
file** (parquet or a `.psp`-style framing; columnar preferred), and fold the LR
into the fixed-size streaming histogram (Q5). (MQDiff still written to the VCF
INFO as today — not consumed by the score.)
- **Deliverable:** pass-1 run over tomato2 produces the temp file + the
  populated LR histogram; a test that the histogram's summary stats match the
  spilled LRs; the temp file lives in scratch and is registered for cleanup.
- *Depends:* Q3, S1. *Source:* arch Premise 6 steps 2–3.

**S3. Calibrate from the histogram.** ☐
From the streaming LR histogram: `ParalogPrior::estimate` → π; build the
`q_of_lr` curve and resolve the operating threshold for the target FDR
(default high-confidence ≈ 1%, introgression-safe; CLI-overridable).
- **Deliverable:** π + threshold computed with no re-read of the temp file; a
  test that they match a full-vector reference on a fixture.
- *Depends:* Q5, S2. *Source:* arch Premise 6 step 3; spec §6.

**S4. Pass 2 — write the final VCF.** ☐
Read the temp file **once**; for each locus stamp `paralog_posterior` /
`FILTER` (and any INFO) from `q_of_lr(LR_i)` vs the threshold; write the final
VCF once (no in-place FILTER patching). Declare the new FILTER/INFO in the
header. Delete the temp file (on success **and** on failure).
- **Deliverable:** end-to-end `var-calling` run on tomato2 emits a VCF with the
  paralog FILTER/INFO; header declarations present; pre-existing FILTER values
  preserved; the temp file is gone afterward.
- *Depends:* S3. *Source:* arch Premise 6 step 3.

**S5. CLI knob.** ☐
A `var-calling` flag for the paralog-filter FDR / mode (e.g.
`--paralog-fdr <q>` with an introgression-safe default; an off switch), threaded
through and recorded in the VCF header for provenance.
- **Deliverable:** flag parsed/validated/defaulted; help text; a test that a
  non-default value reaches the cut and is recorded.
- *Depends:* S4. *Source:* spec §6; `cli.rs`.

---

## Milestone T — end-to-end validation + cost

**T1. End-to-end behaviour on tomato2.** ☐
Confirm the FILTER'd callset matches the expected profile: flagged set is
paralog-like (het excess + coverage excess), planted/known introgression-like
loci are preserved, π ≈ 9%. Compare against the R1 numbers to confirm the
wiring reproduces the isolated statistics.
- **Deliverable:** a short report in `doc/devel/reports/implementations/`.
- *Depends:* S4. *Source:* spec §7.

**T2. Cost / RSS check.** ☐
Confirm the spill-to-disk + streaming histogram keeps RAM **flat in variant
count and sample count** (no genome-wide LR vector); record the wall/RSS delta
and the temp-file size on a representative run.
- **Deliverable:** before/after numbers in the T1 report; a fix only if RAM is
  not flat. *Depends:* S4. *Source:* arch Premise 6; [[low_memory_mode_branch]].

---

## Checkpoints

- After **P1**: the het rate is computable downstream (callable total stored,
  round-trips). Pause.
- After **Q5**: the three statistics pieces are built and unit-tested in
  isolation — the pure core is done. Pause.
- After **R1**: the isolated statistics reproduce the data-first checks on
  tomato2 (peak at 1.0, F Spearman ≈ 0.86, π ≈ 9%, LR correlates with the
  prototype). The maths is trustworthy before any wiring. Pause.
- After **S4**: a real `var-calling` run emits a FILTER'd VCF with bounded RAM
  and a cleaned-up temp file — first end-to-end artefact. Pause.
- After **T2**: cost confirmed flat; the filter is production-shaped.

## Open items carried from the architecture / spec

- **Absolute-FDR ground truth** (read-level simulation / known duplications) —
  the honest ceiling on any absolute-rate claim; not buildable here (spec §9).
- Temp-file format (parquet vs `.psp`-style) + LR-histogram bin resolution
  (arch Premise 6).
- Whether σ₀ is per-sample-fit or a single cohort value in production
  (default per-sample fit — Q2; arch Premise 1).
- `n_ambiguous`-weighted shrinkage of `obs_het` for low-coverage samples
  (arch Premise 3, later).
- Re-check the sample-independence assumption on a highly-related panel (F2 /
  family); down-weight to an effective sample size only if it bites (spec §7).
