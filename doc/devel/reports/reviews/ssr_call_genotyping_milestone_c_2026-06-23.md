# Code Review: ssr-call genotyping+pre-pass — Milestone C
**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator, focused inline pass)
**Scope:** Milestone C (commits `19b1c61` C1, `35555a5` C2–C4)
**Status:** Approve-with-changes

---

## 1. Scope
- **In-scope:** [candidate_set.rs](../../../../src/ssr/cohort/candidate_set.rs),
  [likelihood.rs](../../../../src/ssr/cohort/likelihood.rs),
  [allele_freq_prior.rs](../../../../src/ssr/cohort/allele_freq_prior.rs),
  [em_init.rs](../../../../src/ssr/cohort/em_init.rs), [em.rs](../../../../src/ssr/cohort/em.rs),
  [vcf_out.rs](../../../../src/ssr/cohort/vcf_out.rs), the `Rungs::modal_length` +
  `sample_clear_peaks` additions, and the `CohortLocus.ref_tract` extension.
- **Categories:** reliability, errors, naming, idiomatic, refactor_safety, smells, extras
  (numerical correctness / probability normalization). Skipped unsafe_concurrency,
  tooling.

## 2. Verdict
**Approve-with-changes.** The skeleton is correct: checkpoint 1 holds (genotypes match
truth), the genotype prior is provably normalized for any `F`, the likelihood sums to a
probability, and `G₀` keeps every `π > 0`. Findings are robustness/clarity Minors and a
moderate-stutter robustness test — no Blockers/Majors.

## 3. Execution status
- `cargo fmt --check` → pass · `cargo clippy --all-targets --all-features -- -D warnings`
  → pass · `cargo test --all-features` → **1231 lib pass**.
- Needs-verification findings: 0.

## 4. Top 3 priorities
1. **Mi1** — `run_locus_em` panics (`assert_eq!`) on non-diploid ploidy; confirm this is
   the intended loud-failure behavior and document it (resolved below in favor of the
   panic, per the no-silent-fallback invariant).
2. **Mi2** — `candidate_of_length` / `cohort_representative` pick the *first* sequence at
   a length; document that this is exact for pure tracts and an approximation when a rung
   holds multiple (impure) variants.
3. **MT-1** — add an EM robustness test under *moderate stutter* (clean calls under a
   higher level), not only the near-noiseless checkpoint cohort.

## 5. Findings

### Minor

- `src/ssr/cohort/em.rs` (`run_locus_em`) — **[Minor]** panics on non-diploid ploidy
- **Confidence:** High
- **Problem:** `assert_eq!(ploidy, 2, …)` aborts on any non-diploid call. For a library
  function a panic on an argument value is sharp, but a *silent* no-call would violate the
  project's no-silent-fallback invariant (cross-cutting §5).
- **Why it matters:** The behavior should be a deliberate, documented choice, not an
  accident of the diploid-only implementation.
- **Suggested fix:** keep the loud panic (a non-diploid call in a diploid-v1 build is a
  programming error, and loud-fail beats silent no-call), but document it as the chosen
  contract on the function and flag ploidy generalization as the follow-up.

- `src/ssr/cohort/em_init.rs` / `candidate_set.rs` — **[Minor]** first-of-length sequence selection
- **Confidence:** Medium
- **Problem:** `candidate_of_length` (π⁰ tally) and `cohort_representative` (C1 nomination)
  resolve a length to the *first* / most-supported single sequence. When a rung legitimately
  holds several same-length sequences (substitution/interruption variants — first-class per
  spec §5/§7), only one is tracked.
- **Why it matters:** Exact for the pure tracts the simulator emits and the common case;
  an approximation for impure loci that a later milestone should revisit.
- **Suggested fix:** document the single-representative-per-length simplification on both
  functions; multi-variant-per-rung handling is a follow-up (ties to S2 impure alleles).

### Nits
- `src/ssr/cohort/likelihood.rs` (`read_likelihood`) — every `Δ` with `S_θ > 0` runs
  `reach_variants` + `align_subst`, including far slips whose kernel weight is negligible.
  A small-`S_θ` skip would save work; deferred to F1 (perf), noted here.

## 6. Out of scope observations
- `benches/psp_writer_perf.rs:386` — pre-existing bench panic, unchanged.

## 7. Missing tests to add now
- `em_calls_correct_genotypes_under_moderate_stutter` — **input:** the checkpoint cohort
  re-run with a higher stutter level (e.g. `baseline 0.15`) and matching params.
  **Bug it catches:** an EM that only works at near-zero stutter (e.g. a mis-weighted
  outlier term or a too-greedy MAP) would pass the clean checkpoint but fail here.
  **Body:** assert the homozygotes and separated hets still call correctly.

## 8. What's good
- The genotype prior is normalized for *any* `F` (`Σ_G P(G) = F·1 + (1−F)·1 = 1`), and the
  EM exploits fixed `ε`/`θ`/level to precompute the per-genotype data-likelihood so only
  `π` iterates — a clean, fast structure ([em.rs](../../../../src/ssr/cohort/em.rs)).
- Checkpoint 1 drives the *real* merger end-to-end, not a shortcut, so it exercises the
  whole spine ([em.rs](../../../../src/ssr/cohort/em.rs) tests).
- `G₀`'s explicit floor is unit-tested to keep a length-200 candidate `> 0`
  ([allele_freq_prior.rs](../../../../src/ssr/cohort/allele_freq_prior.rs)).
- The Q-G2 decision is documented at the decision site with its rationale
  ([em.rs](../../../../src/ssr/cohort/em.rs) header).

## 9. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --all-targets --all-features -- -D warnings` ·
  `cargo test --all-features`

### Author response convention
Address each finding by ID (Mi1, Mi2, Nit, MT-1).
