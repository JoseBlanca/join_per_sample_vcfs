# Fixes applied — posterior_engine review 2026-05-16

**Date:** 2026-05-16
**Reviewer pass:** [posterior_engine_2026-05-16.md](posterior_engine_2026-05-16.md)
**Scope:** [src/var_calling/posterior_engine.rs](../../../src/var_calling/posterior_engine.rs), [tests/posterior_engine_integration.rs](../../../tests/posterior_engine_integration.rs)
**Status:** fixes-applied (review-followup)

---

## Open questions

The review surfaced four open questions; resolutions applied here:

1. **`EmDiagnostics.converged` contract.** Adopted "always converged on success, never converged on failure" by dropping the field — see Mi1. The `Iterator` impl returns `DidNotConverge` on the negative branch; a `PosteriorRecord` is only emitted when EM converged, so the redundant flag carried no information. The follow-up `EmOutcome` enum (soft-cap variant) lands when there is a real need for it.
2. **Engine-side vs CLI-side config validation.** **Deferred.** This needs a product call on whether the engine should reject ill-formed configs or whether the CLI is the validation boundary. The engine still validates the *internal-record* shape it depends on (M2 — `MalformedMergedRecord`) so that bad upstream records can no longer panic in release builds; the `f`/pseudocount-range validation (M3) waits for the CLI plan.
3. **`f̂_C` substitution policy for chain-evident samples.** Locked in by test M9 — the new `m_step_compound_frequency_matches_closed_form_with_no_observed_compound_counts` asserts the formula with both samples chain-evident, anchoring the spec's "always substitute" reading.
4. **`1 − α_compound` Beta-partner choice.** Documented in the `DEFAULT_COMPOUND_ALT_PSEUDOCOUNT` const doc with a source pointer to the implementation plan §"Step 4 — M-step on `f̂_C`". A single CLI knob remains correct; documentation is now sufficient to surface the implicit coupling.

## Per-finding disposition

Anchors refer to the unified review report's identifiers.

### Blocker

| ID | Disposition | Notes |
|---|---|---|
| **B1** | fixed | Added `e_step_returns_non_finite_posterior_when_likelihood_row_is_all_negative_infinity` and `classify_nonfinite_distinguishes_nan_pos_inf_neg_inf` unit tests. Both now lock in the only signal that catches malformed upstream data. |
| **B2** | fixed | Reworked `trivial_posterior_record` to fill `posteriors[s][HOM_REF_GENOTYPE_IDX] = 1.0` and zero elsewhere (was `vec![1.0; n_samples * n_genotypes]` — would have violated the row-sum invariant for any `n_genotypes > 1`). Added `trivial_posterior_record_returns_hom_ref_per_sample_when_only_ref_allele_present` and `trivial_posterior_record_returns_empty_allele_frequencies_when_alleles_empty`. |

### Major

| ID | Disposition | Notes |
|---|---|---|
| **M1** | fixed | `NonFiniteKind` gained `NegativeInfinity`; `classify_nonfinite` rewritten as the explicit three-way classifier. `NonFiniteKind` also gained a `Display` impl so the error message reads `kind: -infinity` instead of `kind: NegativeInfinity`. |
| **M2** | fixed | New `MalformedMergedRecord { locus, reason }` variant with a `MalformedMergedRecordReason` enum (GenotypeCountMismatch / LogLikelihoodLengthMismatch / ChainAnchorLengthMismatch). `validate_record_shape` runs unconditionally before EM. Three new tests pin each branch. |
| **M3** | deferred | Engine-side validation of `F`/pseudocount knobs needs Open Question 2 resolved. Internal-record validation (M2) covers the release-only panic risk; user-input validation can land with the cohort CLI subcommand. |
| **M4** | fixed | New `fixation_index_overrides_apply_per_sample_when_length_matches_record` test. Likelihoods strengthened so the F=0 sample calls het and the F=1 sample is forced to a homozygote — exercises both the happy path of the override resolver and the prior-shape difference between F values. |
| **M5** | fixed | `PosteriorEngine` holds `config: PosteriorEngineConfig` by value (was `Arc<PosteriorEngineConfig>` with no second owner). Debug impl no longer reborrows through `&**`. |
| **M6** | fixed | Introduced `RecordLocus` (chrom_id + start + end + `Display`) and `EmContext<'a>` (per-record dims + precomputed tables). `e_step`, `m_step_p_hat`, `m_step_f_hat_compound`, `run_em_loop` now take `(ctx, …)` and `trivial_posterior_record` takes a `TrivialRecordInputs` bundle. All four `#[allow(clippy::too_many_arguments)]` are removed. Locus also collapsed the `chrom_id/start/end` triple in every error variant. |
| **M7** | fixed | Added `posteriors_row_returns_correct_slice_for_last_sample_index`, `scalars_row_uses_n_alleles_stride_not_n_genotypes`, `chain_anchor_flags_row_returns_n_alleles_long_slice`. All three `*_row` accessors gained `# Panics` doc sections. |
| **M8** | fixed | Added `summarise_posteriors_returns_infinite_qual_when_any_sample_has_zero_hom_ref_posterior`. The combined assertion also covers the `MAX_GQ_PHRED` clamp branch. |
| **M9** | fixed | Added `m_step_compound_frequency_matches_closed_form_with_no_observed_compound_counts` — asserts `f̂_C = α / (α + (1 − α) + 4)` on the controlled 2-sample / 4-chromosome fixture. |
| **M10** | fixed | New `proptest_record_strategy()` walks `ploidy ∈ {2..=4}` × `n_alleles ∈ {2..=3}` × `n_samples ∈ {1..=4}`. The three sum-to-1 / simplex / permutation-invariance proptests now run at every shape in that matrix instead of diploid-biallelic only. |
| **M11** | fixed | Added `PosteriorEngineConfig::with_project_defaults()` whose rustdoc enumerates every default; `impl Default` delegates. Field docs on `PosteriorEngineConfig` now reference each `DEFAULT_*` constant by name. `PosteriorEngine::new` doc points at `with_project_defaults`. The `new` constructor is retained so test/integration call sites don't need to change. |

### Minor

| ID | Disposition | Notes |
|---|---|---|
| **Mi1** | fixed | Dropped `EmDiagnostics.converged` — the field was structurally always `true` in emitted records. The follow-up `EmOutcome` enum (soft-cap variant) lands when motivated. |
| **Mi2** | fixed (composed with M6) | `run_em_for_record` now orchestrates `validate_record_shape` → `resolve_fixation_indices` → setup → `run_em_loop` → final E-step → `summarise_posteriors`. The EM-loop body is its own function. |
| **Mi3** | fixed (composed with M6/Mi4) | `e_step`, `m_step_p_hat`, `summarise_posteriors` now use `chunks_exact` / `chunks_exact_mut` over the flat row-major buffers instead of manual stride arithmetic. |
| **Mi4** | fixed | `EmScratch { p_effective, log_p_effective, log_post_unnorm, expected_counts }` is allocated once in `run_em_for_record` and reused across iterations. `e_step` now fuses the `p_effective` / `log_p_effective` write into a single pass over `n_alleles`. |
| **Mi5** | deferred | The criterion bench is a substantive standalone follow-up that would benefit from cohort-CLI fixtures. Tracked as open work on the Stage 6 block in `PROJECT_STATUS.md`. |
| **Mi6** | deferred | The golden test follows the same logic as Mi5 — defer until the cohort CLI lands and there is a real pipeline to drive fixture generation. Tracked as open work. |
| **Mi7** | fixed | Each `DEFAULT_*` pseudocount constant + `MAX_GQ_PHRED` doc comment now carries a `Source:` paragraph pointing at the implementation plan / architecture spec section. |
| **Mi8** | fixed | `PosteriorEngine::config` gained a doc comment. `PosteriorEngine` struct doc gained `# Errors` paragraph and a runnable `# Examples` doctest. `PosteriorEngine::new` and `PosteriorEngine::with_config` carry doc comments. |
| **Mi9** | fixed | `resolve_fixation_indices` takes `Option<&[f64]>` (was `&Option<Vec<f64>>`); body rewritten with `let ... else` and uses `to_vec()` instead of `.clone()`. Caller updates from `&config.fixation_index_overrides` to `config.fixation_index_overrides.as_deref()`. |
| **Mi10** | fixed | `PosteriorEngine.done` renamed to `is_latched`. Updated the four use sites and the Debug destructure. |
| **Mi11** | fixed | `genotype_allele_counts` is now a flat `Vec<u32>` of length `n_genotypes * n_alleles`. Walked via `chunks_exact(n_alleles)` in `e_step` / `m_step_p_hat`. |
| **Mi12** | deferred | Moving the config-mode invariant checks into a `with_config -> Result<Self, _>` is a signature change that benefits from being landed with M3 (config validation). Cost today is one `bool` check + one `Option::is_some` per record — negligible. |
| **Mi13** | deferred | Test-fixture consolidation across crate / integration crate boundary needs a `pub(crate)` test-support module (or feature flag). Defer; two copies remain manageable. |

### Nits

Applied:

- `value` → `posterior` in `e_step`.
- `total` → `dirichlet_denominator` in `m_step_p_hat`.
- `ov` → `overrides` (collapsed alongside Mi9).
- `e_n_c` → `expected_compound_count` in `m_step_f_hat_compound`.
- `HOM_REF_GENOTYPE_IDX: usize = 0` const promoted to module scope; used in both `summarise_posteriors` and `trivial_posterior_record`.
- `PHRED_SCALE: f64 = -10.0` const promoted to module scope; used in both `summarise_posteriors`'s GQ and QUAL conversions.
- Tightened `hard_iteration_cap_returns_did_not_converge_error` to assert `locus.chrom_id`, `locus.start`, `locus.end`, `last_delta.is_finite()` in addition to `max_iterations`.
- `Display` impl on `NonFiniteKind` so error messages no longer use `{kind:?}`.
- `genotype_idx: Option<usize>` on `NonFinitePosterior` (was hard-coded to `0` for the `log_z` branch).

Not applied:

- `approx` helper → `float-cmp` crate dep: adding a dependency for a 3-line test helper isn't worth it.

### Convergent / disputed

- **`#[allow(clippy::field_reassign_with_default)]` workaround in `tests/posterior_engine_integration.rs`** — disputed. The `idiomatic` agent claimed `..Default::default()` works across crates on `#[non_exhaustive]`. Verified empirically against `cargo build` (E0639: "cannot create non-exhaustive struct using struct expression") — FRU is rejected across crates on `#[non_exhaustive]` structs. The workaround is load-bearing; comment retained, wording sharpened to "`#[non_exhaustive]` on `PosteriorEngineConfig` disallows struct-update syntax across crates".

## Validation

All commands run in the project sandbox (no podman in this environment).

- `cargo fmt --check` — passes (no diff).
- `cargo build --all-targets` — passes.
- `cargo test --all-targets --all-features` — **622 lib tests pass** (up from 600; +22 new tests for B1, B2, M2 × 3, M4, M7 × 3, M8, M9, M10 × 3, Mi8, plus the `engine_config_accessor`, `posterior_engine_debug`, `engine_emits_all_successful_records_before_latched_upstream_error`, `log_sum_exp_slice_returns_neg_infinity_when_every_entry_is_neg_infinity`, `log_sum_exp_slice_returns_neg_infinity_for_empty_slice`, `log_sum_exp_slice_matches_log_sum_exp_2_on_two_element_inputs`, `safe_ln_returns_neg_infinity_for_zero_and_negative_inputs`, `max_abs_diff_returns_largest_absolute_difference`, `pseudocount_for_returns_class_specific_value_from_config`, `classify_nonfinite_distinguishes_nan_pos_inf_neg_inf` helpers). **7 integration tests pass** (up from 5; added `upstream_error_propagates_as_posterior_engine_upstream_error`, `tetraploid_record_emits_posteriors_with_correct_n_genotypes`). All pre-existing suites green.
- `cargo clippy --all-targets --all-features` (without `-D warnings`) — **zero warnings on the in-scope files**. The pre-existing warnings in `src/var_calling/per_group_merger.rs:1786 + 1796` and `benches/var_calling_perf.rs` remain, out of scope per the original review.

## Open follow-ups for the next pass

- **M3 — config validation** (CLI-vs-engine boundary decision).
- **Mi1 — soft-cap `EmOutcome`** when there is a real product need to emit records that didn't fully converge.
- **Mi5 — `benches/posterior_engine_perf.rs`** with a regression-threshold harness, ideally with cohort-CLI fixtures.
- **Mi6 — golden `tests/golden/posterior_engine/*` fixtures** locking the emitted f64 values across a (single-sample, weak-het cohort, triploid, compound chain-evident, compound chain-broken) matrix.
- **Mi12 — `with_config -> Result<Self, _>`** lands with M3.
- **Mi13 — fixture consolidation** across the unit / integration crate boundary.
