# SSR interrupted-repeat recall — implementation plan

**Status:** draft, 2026-07-06, branch `ssr-interruptions`. The build order for making
`ssr-call` call the interrupted repeats it currently drops. Design is settled in spec
[ssr_interrupted_repeat_recall.md](../specs/ssr_interrupted_repeat_recall.md); the
code-layout companion is
[architecture/ssr_interrupted_repeat_recall.md](../architecture/ssr_interrupted_repeat_recall.md).
This turns that design into ordered steps and checkpoints. It sits *inside* the shipped
`ssr-call` stage (roadmap [ssr_call_roadmap.md](ssr_call_roadmap.md) Milestones A–F), not
before it.

---

## How the order was chosen (principles)

- **Types first, then implementation**, within every step (project rule).
- **Phase 1 before Phase 2, always.** Phase 1 (structural, sequence-keyed alleles) is the
  measurement instrument Phase 2 (statistical, per-allele stutter) fits on — there are no
  per-allele impure skirts to regress, and no purity covariate to identify, until impure
  alleles are being called (spec §6.4). Each phase ships and validates on its own.
- **Measure before building Phase 2.** Phase 2 is measurement-gated: its model *form*
  (purity separable from length? level or decay?) is decided from Phase 1's data, not
  assumed (spec §6.5 / §8 step 5). If purity moves nothing, Phase 2 stops — the honest
  outcome.
- **Conformance, not re-architecture.** S3 (`Qᵣ`) already scores sequences; the EM,
  `G₀`, enumeration, `is_variable`, ALT emission are already candidate-index keyed. The
  work is converting the four length-keyed holdouts (nomination, seed, slip attribution,
  allele balance) to key on the sequence — verify the "already conformant" pieces stay
  untouched.
- **Determinism throughout.** Every new sufficient statistic is an integer count or a
  fixed-order float reduce; every tie breaks on sequence bytes / candidate index. Prove
  byte-identity across thread counts is preserved (it is a hard requirement for SSR output).
- **Validate on the unit fixture *and* the benchmark.** The `ssr_tomato1` concordance +
  the SNP caller's end-to-end tests are the regression gate (pre-alpha, no byte-identity of
  SSR output required — but cross-thread determinism is).

## Preconditions (already in place)

- The `ssr-call` stage is built end-to-end (roadmap A–F): reader → candidates → `Qᵣ` → EM
  → VCF, with the pre-pass, the `F`/level outer loop, and FP control.
- **S3 is already sequence-aware** — `compute_data_ll` scores each observed sequence against
  each candidate *sequence*; two same-length candidates already get different likelihoods
  (spec §4, the linchpin).
- **The rung ladder already stores the per-length sequence set** (`seqs_by_length`,
  byte-sorted) — Phase 1 exposes it, it does not build it.
- Diagnostic tools exist on this branch: `examples/ssr_psp_seqdump.rs`,
  `examples/ssr_psp_dump.rs`, and the ours-vs-HipSTR dashboard `ssr_vs_hipstr_dashboard.py`.

---

## The steps

### Milestone P1 — Phase 1: sequence-keyed alleles (recover the drops)

**P1.1. Types & thresholds.**  ☐ types
Turn the per-length representative into a per-length *set* and add the admission knobs, no
logic yet: `rung_ladder.rs` gains the per-(length, sequence) **distinct-sample tally** (the
one new sufficient statistic — widen `SeqCount` or a side `BTreeMap`, Q-I4); `candidate_set.rs`
`CandidateCfg` gains `min_same_length_reads` / `min_same_length_samples` /
`min_same_length_fraction` (dev 8 / 3 / 0.10); `SampleCall` gains `allele_support` (the
deconvolved per-allele responsibility, arch §3). Settle the `cohort_alleles` and
`nearest_called_by_sequence` / `allele_responsibilities` signatures before touching callers.
*Depends:* — (first step). *Source:* spec §5.1/§5.2; arch §2.1–§2.5.

**P1.2. Set-valued nomination + §5.2 admission.**  ☐ arch ☐ plan
`candidate_set.rs`: replace `cohort_representative` with `cohort_alleles`; union every
same-length sequence clearing the three-way `&&` bar; reference allele seeded
unconditionally; `MAX_CANDIDATE_ALLELES` still caps. `build_rungs` fills the distinct-sample
tally. *Depends:* P1.1. *Source:* spec §5.1/§5.2; arch §2.1/§2.2.
*Tests:* a pure-vs-interrupted same-length fixture yields **two** candidates and admits; the
singleton noise variants fail on every count; a same-length ALT equal to REF is not
duplicated; the `4279322`-shaped case (26 reads / 17 samples / ~0.10) clears.

**P1.3. Sequence-aware seed.**  ☐ arch ☐ plan
`em_init.rs`: `candidate_of_length` → match the sample's representative sequence (fallback:
closest by `align_subst`). *Depends:* P1.1. *Source:* spec §5.3; arch §2.3.
*Tests:* a same-length homozygote seeds its own composition; a same-length **het** seeds a
hom for the majority (the documented limitation — recovery is P1.4's EM path, asserted in P1.5).

**P1.4. Sequence-aware attribution + allele balance → first interrupted call.**  ☐ arch ☐ plan
The core issue-2 fix. `attribution.rs`: add `nearest_called_by_sequence` (hard, integer slip
stats) and `allele_responsibilities` (soft), both over `align_subst`. `em.rs`:
`attribute_locus` uses the hard primitive (slip stats effect-neutral for same-length alleles,
but length-keying gone); `final_calls` computes the deconvolved `(n_A, n_B)` for the called
genotype and stores `allele_support`. `vcf_out.rs`: `allele_balance` triggers on distinct
candidate **indices** (drop the `genotype_units`-equal short-circuit) and consumes
`allele_support` — no length attribution at the VCF stage.
**Milestone: a same-length interruption het emits as a variable VCF row end-to-end, with a
correct allele-balance penalty.** *Depends:* P1.2, P1.3. *Source:* spec §5.3/§5.4 + Mark-2
§6 amendment; arch §2.4/§2.5/§3.
*Tests:* the balance term is *active* (not `1.0`) for a same-length het and penalizes an
80/20 split; determinism — byte-identical output across thread counts on a same-length locus.

**P1.5. Validate Phase 1.**  ☐ plan
The exit gate. *Depends:* P1.4. *Source:* spec §9.
- **Unit** — a cohort fixture with a same-length interruption polymorphism, **including a
  genuine same-length *het* carrier at low depth (~2–5 reads)**, must yield two candidates,
  a variable call, **and that carrier called *het*** (exercises the EM recovery the seed
  can't do, spec §5.3).
- **Benchmark** — re-run `ssr-call` on `ssr_tomato1` + the dashboard. **Exit criteria:** the
  ~23 % of the 765 drops carrying an interior interruption recover as emitted variable loci,
  **and** the 847 already-common PASS loci hold ≈96.7 % genotype concordance. Report
  per-sample **het** concordance (the undercall tail, spec §5.3), scored on the **sequence**
  genotype (compareSTR `metric-conc-seq`, statSTR `--use-length` off — spec §9/§5.4).
- **Regression watch** — no currently-`Pass` locus flips to `TooManyAlleles` (spec §10,
  issue-4 watch-item); the SNP caller's end-to-end tests still pass.
- **Threshold sweep** — sweep the §5.2 defaults (8 / 3 / 0.10) against recovered-loci vs new
  false positives; pick the knee.

**▶ Checkpoint 1 (after P1.4):** an interrupted locus comes out of real evidence as a variable
row, with an active sequence-aware allele balance. Proves the structural spine.
**▶ Checkpoint 2 (after P1.5):** the benchmark recovers the interruption drops without
regressing the common loci. Phase 1 is done and shippable.

### Milestone P2 — Phase 2: per-allele stutter (HipSTR-parity), measurement-gated

Only after P1.5 is clean.

**P2.0. Measure first (the gate).**  ☐ plan
With Phase 1 calling impure alleles, read out the observed per-allele slip rates (P1.4's
attribution already records them) and settle the model form before any type is written:
(i) **does purity lower the rate at *fixed length*** — on Phase 1's same-length pure-vs-impure
pairs, so purity is identified free of the length confound (spec §6.5a); (ii) **does purity
move the overall level or the slip-size *decay*** — compare the |Δ| distribution
pure-vs-impure (spec §6.5b); (iii) is there per-allele variation left after
length + period + motif + purity? **Gate:** keep a covariate only where it moves the data
(drop purity if it fails (i)); attach the purity term per (ii); add a per-allele deviation
term only if a residual survives (iii). **If nothing moves the rate, stop — Phase 2 is
unnecessary.** *Depends:* P1.5. *Source:* spec §6.5/§8 step 5.

**P2.1. Types.**  ☐ types
The per-allele purity measure (form from P2.0, Q-I2); the covariate default
`f(length, period, motif, purity)` and a per-allele rate field in `param_estimation.rs`,
shrinking toward `f`; the chosen attach point (level / decay / both, from P2.0, Q-I3).
*Depends:* P2.0. *Source:* spec §6.3; arch §5.

**P2.2. Estimate + apply.**  ☐ arch ☐ plan
Fit `f` on the pooled skirts and the per-allele shrinkage in the pre-pass; apply the
per-allele term at the read-model per-candidate site (level or decay per P2.0). Collapses to
Phase 1 identically where no impure alleles exist. *Depends:* P2.1. *Source:* spec §6.3/§6.4;
arch §5.

**P2.3. Validate Phase 2.**  ☐ plan
A fixture where the impure allele's known-lower slip rate must be recovered; re-run the
benchmark. **Exit criterion:** genotype-quality calibration on impure loci improves (or
matches HipSTR) with **no regression on pure loci and no loss of the Phase 1 recall gain**.
*Depends:* P2.2. *Source:* spec §9.

---

## Cross-cutting gates (hold the line)

- **Determinism** — byte-identity across thread counts on same-length loci, checked at P1.4
  and again after P2.2 (arch §4). Non-negotiable for SSR output.
- **Paralog FP audit** — until the same-length paralog guard exists (out of scope, tracked in
  `doc/devel/TODO.txt` + arch §6), P1.5 must include an audit flagging cohort-universal
  same-length hets as suspected paralog artifacts (spec §10). This gates whether the §5.2
  defaults are safe, not just tuned.
- **"Already conformant" review** — confirm `g0_pseudocounts`, `enumerate_diploid_genotypes`,
  `genotype_prior`, `run_pi_em`, `is_variable`, `site_qual`, and ALT emission are untouched
  (arch §2.6). A change there is a scope-creep smell.

---

## What this plan is *not*

- **Not** the Mark-2 BIC confident-genotype gate (roadmap D1's model-based form) — Phase 1 is
  its prerequisite, but landing Phase 1 does not remove the `ε` contamination the length-based
  gate leaves (spec §5.5); tracked separately (TODO).
- **Not** the paralog filter — the sharpest new FP mode, adapted from
  `src/var_calling/paralog_filter/`, is its own effort (TODO).
- **Not** the rare *length*-alt suppression behind the other portion of the 765 drops (spec
  §11) — a prior/threshold matter tracked with the short-period effort.
