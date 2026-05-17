# Implementation report — Stage 6 posterior engine (no-contamination v1)

**Date:** 2026-05-16
**Plan:** [posterior_engine.md](../../implementation_plans/posterior_engine.md)
**Skill:** rust-feature-implementation
**Status:** implemented (review pending)

---

## 1. Plan

Build the no-contamination posterior engine described in the plan: a
streaming iterator that consumes
`Result<MergedRecord, PerGroupMergerError>` from Stage 5 and emits one
`PosteriorRecord` per upstream record carrying the per-sample,
per-genotype posterior table, the converged per-allele frequency `p̂`,
the per-compound cohort frequency `f̂_C`, the per-sample best genotype
and GQ, and the site-level QUAL. Each record runs an independent EM
loop with closed-form M-steps on `p̂` (Dirichlet) and `f̂_C` (Beta),
under an HWE-with-`F` prior generalised to arbitrary ploidy via the
Wright–Fisher partition formula.

## 2. Assumptions (silent choices made during implementation)

The plan has a *"Decisions to confirm before implementation"* section
with eight open questions. The three that affect v1 scope or contract
were taken to the user; the rest were resolved by adopting the plan's
recommended defaults.

User-confirmed choices:

- **Polyploid support.** Engine reads `MergedRecord.ploidy` and
  evaluates the HWE-with-`F` prior at that ploidy via the
  Wright–Fisher partition mixture. Specialises cleanly to the spec's
  diploid biallelic formula at `ploidy = 2`.
- **Non-convergence behaviour.** Hard error
  `PosteriorEngineError::DidNotConverge`. Matches the Stage 4/5
  typed-error precedent; closed-form M-step non-convergence almost
  certainly signals an internal bug rather than a data condition.
- **`PosteriorRecord` shape.** Forwards `scalars`, `other_scalars`,
  and `chain_anchor_flags` from `MergedRecord` so the VCF writer
  needs no parallel join for DP/AD-style FORMAT fields.

Plan-default choices adopted silently (each consistent with the
plan's "Assumptions / silent choices" section):

- **Natural-log throughout**; Phred conversions cross to log10 only
  at the GQ/QUAL boundary.
- **Flat `p̂` and `f̂_C` initialisation** (`1 / n_alleles`).
- **`α_compound`'s Beta partner is `1 − α_compound`**, placing the
  per-compound estimator's prior mean exactly at `α_compound`
  (`= 0.001` default). The plan names `α_compound` but not
  `α_NotC`; this is the simplest interpretable Beta partner choice.
  One CLI knob covers both.
- **`PosteriorEngineConfig::approximate_posterior_calculation`** is
  wired (the plan asks for it) but the LUT implementation is
  out-of-scope for v1; setting it to `true` returns
  `PosteriorEngineError::ApproximateModeNotYetImplemented` so the
  flag is never silently ignored.
- **`PosteriorEngineConfig::contamination`** is wired as
  `Option<ContaminationConfig>`; `Some(_)` returns
  `PosteriorEngineError::ContaminationModeNotYetImplemented`.
  `ContaminationConfig` is an empty `#[non_exhaustive]` placeholder
  for the follow-up plan.
- **`fixation_index_default: f64` + `fixation_index_overrides:
  Option<Vec<f64>>`** instead of the plan's `Vec<f64>`. Engine
  internally builds the per-sample vector from the scalar on first
  record; per-sample F via a supplied file becomes a UI change
  filling `Some(_)` later — no engine refactor.
- **Compound's `f̂_C` substitutes for `p̂[compound]` in the prior for
  every sample, regardless of chain-anchor status.** The plan's
  wording on `ca_flags` is ambiguous; the spec's "VCF encoding of
  the flag" section pins it: both chain-evident and chain-broken
  samples use `f̂_C` as the compound's prior frequency. The
  chain-anchor flag affects only the likelihood, which Stage 5 has
  already baked into `log_likelihoods`.
- **Single-allele record** (`alleles.len() < 2`) returns a trivial
  `PosteriorRecord` with `QUAL = 0`, every sample hom-REF. Stage 5
  already filters these out; this is belt-and-braces.
- **GQ clamped at `MAX_GQ_PHRED = 99`** to keep an exact
  `posterior = 1` from yielding `+∞`.
- **Site QUAL is left as `f64::INFINITY`** when at least one sample
  has `P(hom-ref) = 0` exactly. The VCF writer can cap if it wants
  to; no in-engine cap.

## 3. Changes made

- **New file**: [src/var_calling/posterior_engine.rs](../../../src/var_calling/posterior_engine.rs)
  — `PosteriorEngine<I>`, `PosteriorEngineConfig` (with `DEFAULT_*`
  constants), `PosteriorRecord` (flat row-major posteriors +
  scalars + chain_anchor_flags), `EmDiagnostics`,
  `ContaminationConfig` (empty placeholder), `PosteriorEngineError`
  (six variants), `NonFiniteKind`, the EM loop (`run_em_for_record`),
  the HWE-with-`F` Wright–Fisher prior (`e_step`), the closed-form
  M-step on `p̂` and `f̂_C` (`m_step`), the per-sample GQ / site QUAL
  derivation (`summarise_posteriors`), allele-class routing for
  Dirichlet pseudocounts (`classify_alleles`), and a full
  `#[cfg(test)]` module.
- **Updated** [src/var_calling/mod.rs](../../../src/var_calling/mod.rs)
  to export `posterior_engine`.
- **New file**:
  [tests/posterior_engine_integration.rs](../../../tests/posterior_engine_integration.rs)
  — five integration tests exercising the engine over synthetic
  `MergedRecord` streams.

## 4. Tests added/updated

### Unit tests (28 in `src/var_calling/posterior_engine.rs`)

Numerical helpers:

- `log_sum_exp_2_handles_neg_infinity_arguments`
- `log_factorial_returns_zero_for_zero_and_one`
- `log_multinomial_coefficient_matches_closed_form`
- `homozygous_allele_returns_index_for_pure_genotypes`
- `classify_alleles_routes_pseudocounts_correctly`

EM correctness — strong-evidence calls:

- `single_sample_strong_ref_returns_homozygous_ref_and_low_qual`
- `single_sample_strong_alt_returns_homozygous_alt_and_high_qual`
- `single_sample_uniform_likelihood_leans_on_prior_and_pseudocount`
- `two_samples_opposite_evidence_pulls_p_hat_toward_one_half`

Cohort-pooling behaviour:

- `cohort_pooling_pulls_weak_het_sample_toward_ref_when_pseudocount_dominates`
- `em_converges_within_a_handful_of_iterations`
- `hard_iteration_cap_returns_did_not_converge_error`

Prior shape:

- `inbreeding_one_zeroes_heterozygote_posterior` — `F = 1` ⇒ no
  het mass.
- `snp_vs_indel_pseudocount_routing_pulls_indel_alt_down_more` —
  asserts the per-class pseudocount routing.
- `triploid_run_succeeds_with_polyploid_hwe_prior` — exercises the
  Wright–Fisher partition formula at `ploidy = 3`.

Compound handling:

- `compound_allele_emits_compound_frequency_in_posterior_record`
- `chain_broken_sample_uses_likelihood_unchanged_and_still_emits_compound_frequency`
- `all_samples_chain_broken_at_compound_remains_numerically_stable`

GQ / QUAL closed-form checks:

- `site_qual_matches_closed_form_on_weak_uniform_cohort`
- `gq_phred_matches_closed_form`

Config-mode error variants:

- `approximate_mode_returns_not_yet_implemented_error`
- `contamination_mode_returns_not_yet_implemented_error`
- `inbreeding_override_length_mismatch_surfaces_typed_error`

Iterator semantics:

- `error_latches_engine_after_first_failure`
- `upstream_exhaustion_is_clean_termination_not_error`

Property tests (proptest, 64 cases each):

- `posteriors_sum_to_one_per_sample`
- `p_hat_is_valid_simplex`
- `sample_permutation_preserves_p_hat_and_qual`
- `larger_ref_pseudocount_cannot_increase_p_alt` (monotonicity)

### Integration tests (5 in `tests/posterior_engine_integration.rs`)

- `streams_three_synthetic_records_in_genomic_order` — multi-record
  ordering.
- `cohort_prior_pulls_lone_weak_alt_back_to_ref` — pooling pulls a
  lone weak het to REF.
- `cohort_evidence_overcomes_rare_allele_prior_when_all_samples_agree`
  — pooled het evidence overcomes the rare-allele pseudocount.
- `engine_output_is_bit_identical_across_runs` — reproducibility.
- `config_override_propagates_through_engine` — config plumbing
  past the iterator boundary.

End-to-end PspReader→…→PosteriorEngine integration is **not** in
this commit; building a real `MergedRecord` stream needs a small
ref-fasta fixture and exercises Stage 1–5 plumbing rather than
Stage 6 correctness. Add when the cohort CLI lands and there is a
realistic pipeline to drive end-to-end.

## 5. Validation results

Commands run inside the project sandbox (cargo direct; the dev
container wrapper `./scripts/dev.sh` is unavailable in this
environment — podman is not installed).

- `cargo fmt --check` — passes after applying `cargo fmt`.
- `cargo build --all-targets` — passes.
- `cargo test --all-targets --all-features` — **600 lib tests pass,
  0 fail** (including the 28 new posterior-engine unit tests and 4
  proptests at 64 cases each); 5 new integration tests pass; all
  pre-existing integration suites still pass.
- `cargo clippy --all-targets --all-features` (without
  `-D warnings`) — **zero warnings on the new code**. Two
  pre-existing `needless_range_loop` warnings remain in
  [src/var_calling/per_group_merger.rs:1786 + 1796](../../../src/var_calling/per_group_merger.rs#L1786)
  (the chain-broken-compound fallback) and one pre-existing
  `manual_clamp` / `needless_range_loop` in
  `benches/var_calling_perf.rs`; both predate this implementation
  and are outside the scope of this work. Left for an unrelated
  clippy-cleanup pass.

## 6. Tradeoffs and follow-ups

### Deliberate non-goals (per the plan §"Out of scope")

- **Contamination machinery.** Both Algorithm 3 (in-engine
  scratch-file replay) and Algorithm 5 (post-hoc VCF-replay
  subcommand) are deferred. The engine surface accommodates them
  via `ContaminationConfig`.
- **VCF emission.** `PosteriorRecord` carries everything the
  writer needs (incl. forwarded scalars); the writer is a separate
  module.
- **Rayon parallelism across records.** Each record's EM is
  independent — trivial future rayon-iter. Defer until profiling
  motivates it.
- **CLI wiring.** `--inbreeding`, `--ref-pseudocount`,
  `--snp-alt-pseudocount`, `--indel-alt-pseudocount`,
  `--compound-alt-pseudocount`. Land with the cohort CLI subcommand.
- **Approximate-LUT inner loop.**
  `--approximate-posterior-calculation` flag is wired but its
  implementation is deferred; setting it `true` errors. See the
  plan §"Approximation via precomputed lookup tables" for the
  evaluation methodology.
- **Per-sample F file format.** `fixation_index_overrides` is
  reserved; the CLI today exposes only the scalar. Filling
  `Some(_)` from a supplied file is purely UI, no engine refactor.
- **End-to-end PspReader→PosteriorEngine integration test.** Lands
  with the cohort CLI subcommand.

### Open items surfaced by this implementation

- **Pre-existing clippy warnings in `per_group_merger.rs`** at
  lines 1786 + 1796 (and one in `benches/var_calling_perf.rs`).
  Out of scope for this commit; flag for a future clippy-cleanup
  pass.

### Tradeoffs taken

- **`fixation_index_default` + `Option<Vec<f64>>` overrides**
  instead of the plan's `Vec<f64>`. Cleaner per-record handling
  (engine builds the per-sample vector lazily from the scalar)
  while preserving the plan's "per-sample from day one" upgrade
  path.
- **Compound's `f̂_C` is computed for every record carrying a
  compound allele**, with the result reported in
  `PosteriorRecord.compound_frequencies` as `Some(_)`. The
  per-compound estimator runs even when no sample is chain-broken,
  matching the spec's "VCF encoding of the flag" semantics.
- **Final E-step after the last M-step** so the emitted posteriors
  reflect the converged `(p̂, f̂_C)` rather than the parameters
  that produced the last in-loop posteriors buffer. Small extra
  cost (one pass over samples × genotypes per record); avoids a
  subtle off-by-one between reported posteriors and reported
  `allele_frequencies`.
