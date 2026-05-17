# Code Review: contamination_estimation
**Date:** 2026-05-17
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** commit `0441600` — `var_calling: contamination_estimation side-pass v1 + Stage 6 mixture E-step`
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** single commit `0441600` on `main` of `/home/jose/devel/join_per_sample_vcfs`; the contamination-estimation side-pass module + Stage 6 wiring to consume its output.
- **Reviewed against:** commit `0441600` (HEAD of `main`).
- **In-scope files:**
  - [src/var_calling/contamination_estimation.rs](../../../src/var_calling/contamination_estimation.rs) — NEW, ~1200 lines incl. tests
  - [src/var_calling/posterior_engine.rs](../../../src/var_calling/posterior_engine.rs) — *only the changes in this commit*: `PosteriorEngineConfig.contamination` type change, removal of `ContaminationConfig` placeholder + `ContaminationModeNotYetImplemented` variant, addition of `ContaminationCohortSizeMismatch`, new `compute_mixture_log_likelihoods` + `map_to_contam_class` helpers, mixture branch at the top of `run_em_for_record`, new `try_single` test helper, 3 hardened proptests + 2 replacement tests
  - [src/var_calling/mod.rs](../../../src/var_calling/mod.rs) — one-line `pub mod contamination_estimation;`
  - [tests/contamination_estimation_integration.rs](../../../tests/contamination_estimation_integration.rs) — NEW, 4 integration tests
- **Deliberately out of scope:**
  - [src/var_calling/per_group_merger.rs](../../../src/var_calling/per_group_merger.rs) lines 1786 / 1796 — pre-existing `clippy::needless_range_loop` fix only, semantic-preserving.
  - [benches/var_calling_perf.rs](../../../benches/var_calling_perf.rs) lines ~123 / ~365 — same kind of pre-existing clippy fix.
  - The bulk of `posterior_engine.rs` not changed in this commit.
  - The spec / plan / impl-report Markdown files.
- **Categories dispatched** (8 of 10 — `unsafe_concurrency` skipped because the new code has no `unsafe`/`Arc`/`Mutex`/atomics/`async`/thread spawning; `tooling` skipped because no `Cargo.toml` changes):
  - `reliability` (always)
  - `errors` (always)
  - `naming` (always)
  - `defaults` (new public API + config defaults)
  - `idiomatic` (always)
  - `refactor_safety` (always)
  - `smells` (always)
  - `extras` (parser-adjacent: cohort-position filter; numerical correctness path; on cohort-CLI critical path)

## 2. Verdict

**Request-changes.**

Two Blocker findings (test-coverage gap on the central new numerical function; unvalidated public constructor on the CLI ingestion path) and sixteen Majors. The Blockers must land before this can be considered "shipped"; the Majors are not deferable to a separate follow-up without explicit justification. Several Majors are mechanical and trivially fixable (the test-only `..Default::default()` fixes; the placeholder-field shape change; the const-promotion + deduplication of `MIN_BASE_ERROR`/`MAX_BASE_ERROR`).

## 3. Execution status

Verification commands quoted in the dispatch context from the implementation report; **not re-run by the reviewer** (the implementation run already established a clean state and no in-scope file has changed since):

| Command | Result |
|---|---|
| `cargo fmt --check` | exit 0, clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0, clean |
| `cargo test --all-targets --all-features` | 645 lib + 4 contamination integration + 104 other integration tests pass; 0 failed |

Commands not run by the reviewer: `cargo doc --no-deps`, `cargo audit`. Run them at fix-application time if not part of CI already.

Findings labeled "Needs verification": 0. Two findings (`M3 try_single` rejection rate; `M11 InsufficientSites` recommendation cap branches) reference a verification step the author may want to run; both are filed at High confidence on the structural issue.

## 4. Open questions and assumptions

1. **Is `from_user_supplied` reachable from external CLI input today, or only from a future cohort subcommand?** The impl report says "CLI parser bindings deferred ... until then the side-pass is library-only and tested via the public API." If the constructor genuinely has no external caller today, B1 can land as Major with a "before any external caller wires it up" deferral. If the future cohort subcommand is already in flight, B1 is a true Blocker. **Affects:** B1, Mi5.
2. **Is the v1 intent that `compute_mixture_log_likelihoods` be the production mixture path, or a placeholder for a more sophisticated reconstruction later?** The impl report Assumption 4 sketches the per-allele uniform-error model as "the standard assumption shared with Stage 5's own-DNA likelihood." If a richer reconstruction is planned, the missing value-witness tests (B2) can lock the *current* behaviour as a regression baseline until the upgrade lands. **Affects:** B2, M1, M5, several Minor tests.
3. **Should `min_batch_size_for_contamination = 1` be a hard CLI error or a silent coercion to 2?** The current `.max(2)` in `finalise` silently coerces; M9 prefers explicit rejection. The spec doesn't pick. **Affects:** M9.
4. **What is the project's policy on `Default::default()` vs `with_project_defaults()` for `#[non_exhaustive]` config types?** The existing `posterior_engine.rs` ships both; this review's refactor-safety findings (M13) assume the project convention is "name every field" but the existing code base mixes both. Confirm before mechanical fix. **Affects:** M12, M13.
5. **Is the silent `γ = 0` fallback on degenerate mixture an actual unreachable branch (provable under positive Dirichlet pseudocounts) or a real recovery point?** The error finding (M2c) argues either it should be a `debug_assert!` (if unreachable) or a counted typed error (if reachable). Calibration data may answer. **Affects:** M2.

## 5. Top 3 priorities

1. **B2 — `compute_mixture_log_likelihoods` has no value-witness test.** The only new numerical code path in the consumer side is uncovered by any assertion that would fail if the math were wrong (sign error in `(1-c_s)` vs `c_s`, swapped `g_major_class` vs `class_idx`, transposed `genotype_allele_counts`, wrong `non_genotype_classes` denominator, etc.). Existing tests either short-circuit on `c_s = 0` or only assert posterior-row-sum=1 which holds regardless of formula correctness.
2. **B1 — `ContaminationEstimates::from_user_supplied` accepts arbitrary unvalidated input.** A `c_s = 5.0`, an unnormalised `q_b`, or a length-mismatched `sample_to_batch` is laundered into Stage 6's hot loop and surfaces as either silently wrong posteriors or an `index out of bounds` panic with no `locus`. The constructor is the documented CLI entry point for `--contamination-estimates`.
3. **M1 + M14 — `InformativeObservation { sample_idx: 0 }` placeholder pattern + `OnlineEmState` god-struct.** Both are architectural defects in the central new module surfaced by 5 of 8 reviewers convergent. The placeholder pattern is a load-bearing correctness invariant enforced only by physical proximity in `filter_and_update`; the god-struct's three field groups (running stats / frozen parameters / stability) are touched by disjoint method sets, so the layout actively obstructs further extension (adaptive tolerance, partition-parallel split). Both have mechanical fixes (drop the field; extract sub-structs) that simultaneously fix several Minor findings each.

## 6. Findings

### Blocker

**B1: [src/var_calling/contamination_estimation.rs:305](../../../src/var_calling/contamination_estimation.rs#L305) — `ContaminationEstimates::from_user_supplied` performs zero validation**
- **Categories:** errors. **Convergent evidence:** reliability (M-severity finding on the accessor-side panic propagation).
- **Confidence:** High
- **Problem:** The constructor accepts user-supplied `c_s_per_sample`, `q_b_per_batch`, and `sample_to_batch` and stores them verbatim. It does not check that (a) `c_s_per_sample.len() == sample_to_batch.len()`, (b) every `Some(c)` is finite and in `[0, 1]`, (c) every `q_b` row sums to ~1 or is the explicit floored-zero vector, (d) every `sample_to_batch[s] < q_b_per_batch.len()`, or (e) values are finite. The doc comment names this as the CLI's `--contamination-estimates` / `--contamination-source-distributions` entry point, so the input is genuinely untrusted.
- **Why it matters:** Stage 6 consumes these as frozen mixture parameters in `(1 − c_s) · p_own + c_s · p_contam`. A user supplying `c_s = 5.0` or an unnormalised `q_b` silently feeds degenerate parameters into every per-record EM. `effective_c_s` indexes `c_s_per_sample[sample_idx]` directly (line 323), and `q_b_for_sample` chains two raw indexes through `sample_to_batch` (line 329), so a length mismatch *panics* inside Stage 6's hot loop rather than producing a typed error at the boundary. This is the canonical "errors never pass silently" violation — a malformed user input is laundered into wrong results or a hot-path panic instead of a typed error at the construction site.
- **Suggested fix:** Make the constructor fallible and validate at the boundary; the same validation pattern `estimate_contamination` already uses:
  ```rust
  pub fn from_user_supplied(
      c_s_per_sample: Vec<Option<f64>>,
      q_b_per_batch: Vec<[f64; N_ALLELE_CLASSES]>,
      sample_to_batch: Vec<usize>,
  ) -> Result<Self, ContaminationEstimationError> {
      if c_s_per_sample.len() != sample_to_batch.len() {
          return Err(ContaminationEstimationError::BadInput(format!(
              "c_s_per_sample length {} does not match sample_to_batch length {}",
              c_s_per_sample.len(), sample_to_batch.len(),
          )));
      }
      for (s, c_opt) in c_s_per_sample.iter().enumerate() {
          if let Some(c) = c_opt && !(c.is_finite() && (0.0..=1.0).contains(c)) {
              return Err(ContaminationEstimationError::BadInput(format!(
                  "c_s_per_sample[{s}] = {c} is not a finite value in [0, 1]",
              )));
          }
      }
      for (b, q) in q_b_per_batch.iter().enumerate() {
          let s: f64 = q.iter().sum();
          let all_zero = q.iter().all(|v| *v == 0.0);
          if !all_zero && (s - 1.0).abs() > 1e-6 {
              return Err(ContaminationEstimationError::BadInput(format!(
                  "q_b_per_batch[{b}] sums to {s}; expected 1.0 (or all zeros for floored batches)",
              )));
          }
          for (a, v) in q.iter().enumerate() {
              if !v.is_finite() || !(0.0..=1.0).contains(v) {
                  return Err(ContaminationEstimationError::BadInput(format!(
                      "q_b_per_batch[{b}][{a}] = {v} is not a finite value in [0, 1]",
                  )));
              }
          }
      }
      if let Some(&max_b) = sample_to_batch.iter().max()
          && max_b >= q_b_per_batch.len()
      {
          return Err(ContaminationEstimationError::BadInput(format!(
              "sample_to_batch references batch_idx {max_b} >= q_b_per_batch length {}",
              q_b_per_batch.len(),
          )));
      }
      Ok(Self { c_s_per_sample, q_b_per_batch, sample_to_batch,
                source: ContaminationEstimateSource::UserSupplied })
  }
  ```
  Promoting `BadInput(String)` to a structured variant per M6 would tighten this further.

**B2: [src/var_calling/posterior_engine.rs:675](../../../src/var_calling/posterior_engine.rs#L675) — `compute_mixture_log_likelihoods` has no direct or value-checking test**
- **Categories:** reliability, extras (convergent).
- **Confidence:** High
- **Problem:** `compute_mixture_log_likelihoods` (lines 675–776) is the only place in the engine that consumes `c_s > 0` and produces the mixture log-likelihood table; the rest of the EM is unchanged. The two new in-module tests do not exercise its body in a meaningful way:
  - `contamination_estimates_zero_matches_no_contamination_path` (lines 1939–1958) uses `ContaminationEstimates::zero(...)`, every sample has `c_s_per_sample[s] = None`, `effective_c_s` returns `0.0`, and the function takes the `c_s <= 0.0` early-return branch (line 711) and never executes the inner `g_idx × a` loop.
  - `contamination_cohort_size_mismatch_surfaces_typed_error` (lines 1962–1976) errors out before the helper is called.
  - Integration test #2 `side_pass_output_feeds_stage_6_without_error` ([tests/contamination_estimation_integration.rs:166-215](../../../tests/contamination_estimation_integration.rs#L166-L215)) *does* run the body but only asserts `posteriors_row(s).sum() == 1.0 ± 1e-9` — a property that holds for almost any finite log-likelihood input because the EM normalises out. It cannot witness any value computed by the new formula.

  Result: the entire formula (`p_own`, `p_contam`, `mix`, `mix.ln()`, the genotype-allele-counts table build, the `mean_err` reconstruction from `q_sum / num_obs`, the compound→`None` pathway, the `mix <= 0.0` error path) is uncovered by an assertion that would fail if the math were wrong.
- **Why it matters:** This function is the entire contamination-aware likelihood path. A sign error, off-by-one in `ploidy_f`, swapped `(1 - c_s)` / `c_s`, transposed indexing into `scalars`, or wrong genotype enumeration would not be caught by any current test. The severity rubric's Blocker clause "absence of tests for any code path that, if broken, would produce wrong results without panicking" fits exactly. The impl report claims integration test #2 covers the mixture branch's "public API for the first time," but the assertion cannot distinguish "mixture path ran" from "mixture path was silently bypassed."
- **Suggested fix:** Add the value-witness tests itemised in §8 (notably `compute_mixture_log_likelihoods_matches_hand_calculation_at_c_s_5_percent`, `..._returns_non_finite_error_when_mix_is_zero`, `..._maps_compound_allele_to_zero_q_contam`, and a baseline-vs-contam contrast in the integration test). Concretely fix integration test #2 first:
  ```rust
  let baseline_record = merged_record(0, 100, alleles_seqs.clone(), 2, scalars_per_sample.clone(), lls_per_sample.clone());
  let baseline: Vec<_> = PosteriorEngine::new(std::iter::once(
      Ok::<MergedRecord, PerGroupMergerError>(baseline_record),
  )).collect();
  let pr_c = with_contam[0].as_ref().expect("engine should emit");
  let pr_n = baseline[0].as_ref().expect("baseline should emit");
  let mut any_diff = false;
  for s in 0..pr_c.n_samples {
      for (a, b) in pr_c.posteriors_row(s).iter().zip(pr_n.posteriors_row(s).iter()) {
          if (a - b).abs() > 1e-9 { any_diff = true; }
      }
  }
  assert!(any_diff, "mixture branch produced output identical to c_s=0 path");
  ```

### Major

**M1: [src/var_calling/contamination_estimation.rs:620](../../../src/var_calling/contamination_estimation.rs#L620) — `InformativeObservation { sample_idx: 0 }` placeholder filled in by caller (illegal intermediate state)**
- **Categories:** idiomatic, smells, reliability, naming, errors (convergent across 5).
- **Confidence:** High
- **Problem:** `step_1b_sample_filter` returns the struct with `sample_idx: 0` and a comment "filled in by the caller (we don't have sample_idx here)"; `filter_and_update` then mutates the field via `obs.sample_idx = sample_idx`. The struct briefly carries a value (`0`) that is silently wrong for every sample except sample 0. Any code reading the partial value during the gap — or any new code added between the two lines — silently mishandles every sample but the first. `InformativeObservation` is `Copy`, so the `mut` rebinding makes the placeholder lifetime even more confusing.
- **Why it matters:** "Make illegal states unrepresentable" — the sentinel-and-patch pattern is exactly the failure mode the rule guards. A future caller that forgets the patch silently double-credits sample 0's running statistics with every sample's reads; the bug is invisible at the type level, and no current test would catch it.
- **Suggested fix:** Drop the field and pass `sample_idx` to `apply_observation` directly:
  ```rust
  struct InformativeObservation {
      counts_by_class: [u32; N_ALLELE_CLASSES],
      q_sum_by_class: [f64; N_ALLELE_CLASSES],
      g_major_class: AlleleClass,
  }
  fn apply_observation(
      state: &mut OnlineEmState,
      sample_idx: usize,
      obs: &InformativeObservation,
      sample_to_batch: &[usize],
  ) { /* s = sample_idx */ }
  // in filter_and_update:
  let Some(obs) = step_1b_sample_filter(record, config) else { continue };
  apply_observation(state, sample_idx, &obs, sample_to_batch);
  ```

**M2: [src/var_calling/contamination_estimation.rs:746](../../../src/var_calling/contamination_estimation.rs#L746) — Silent fallbacks at three numerical guards swallow model-invariant signals**
- **Categories:** errors. **Convergent:** reliability (cross-category note); extras (Minor on subset).
- **Confidence:** Medium (each individual site; High that the *pattern* is wrong)
- **Problem:** Three silent-recovery sites in the EM hot path:
  - (a) Line 746 — `let mean_error = mean_log_error.exp().clamp(MIN_BASE_ERROR, MAX_BASE_ERROR);`. The comment says the floor exists because "Some BQ aggregates can underflow to near zero" and the ceiling is "a safety net" — both are recovery points; if a real input ever exceeds either bound, the side-pass silently mutates and continues. The matching site at [posterior_engine.rs:734](../../../src/var_calling/posterior_engine.rs#L734) (`per_read_log_err.exp().clamp(1e-12, 0.5)`) has the same critique *and* inlines the magic numbers (see M4).
  - (b) Line 766 — `let gamma = if mix > f64::MIN_POSITIVE { (c_s * p_contam) / mix } else { 0.0 };`. Comment correctly identifies that `p_own = p_contam = 0` is the only way to reach the branch under the uniform-error model, but after the `mean_error` clamp this branch is *structurally unreachable*. The silent γ = 0 produces a downward bias on `c_s` without any signal to the user.
  - (c) Line 800 — `*slot = (s_val / n_val as f64).clamp(0.0, 1.0);`. Each `γ_r ∈ [0, 1]` by construction, so the ratio is bounded in [0, 1] mathematically. If the clamp ever changes the value, an upstream invariant is violated; silently clipping to the boundary produces plausible-looking `c_s` while hiding the bug.
- **Why it matters:** Project memory note: "For silent-fallback / invariant-violation sites, prefer typed-error variants over `tracing::warn!`." The current three sites all silently rescue conditions that the math says shouldn't happen. Combined with B2 (no value-witness tests), the EM can drift into a regime the spec rules out and nobody notices.
- **Suggested fix:** For each site, decide whether the bound is a mathematical invariant or a real recovery point.
  - (a) If `MIN_BASE_ERROR`/`MAX_BASE_ERROR` are real safety nets, add a counter in `OnlineEmState` and surface it on `ContaminationEstimates`; if they're invariants, `debug_assert!` before the clamp.
  - (b) Promote to typed error variant (`ContaminationEstimationError::DegenerateMixture { sample_idx, batch_idx, class_idx }`) or at minimum a `debug_assert!(mix > 0.0, ...)`.
  - (c)
    ```rust
    if n_val > 0 {
        let ratio = s_val / n_val as f64;
        debug_assert!(
            (0.0..=1.0).contains(&ratio) && ratio.is_finite(),
            "c_s ratio {ratio} out of [0, 1] for sample (s_val={s_val}, n_val={n_val})"
        );
        *slot = ratio.clamp(0.0, 1.0); // release-mode safety net
    }
    ```

**M3: [src/var_calling/posterior_engine.rs:2473](../../../src/var_calling/posterior_engine.rs#L2473) — Hardened proptests silently skip *all* engine errors, not just `DidNotConverge`**
- **Categories:** errors, reliability, refactor_safety (convergent).
- **Confidence:** High
- **Problem:** The three hardened proptests (`posteriors_sum_to_one_at_any_ploidy_and_n_alleles`, `p_hat_is_a_valid_simplex_at_any_ploidy_and_n_alleles`, `sample_permutation_preserves_p_hat_and_qual_at_any_ploidy_and_n_alleles`) use the pattern `let Ok(pr) = try_single(record) else { return Ok(()); };`. The justifying comment ("EM can legitimately return `DidNotConverge` ... not a bug") names only `DidNotConverge`, but the `let-else` discards *every* `PosteriorEngineError` variant — `MalformedMergedRecord`, `NonFinitePosterior`, `ContaminationCohortSizeMismatch`, `Upstream`, etc. — all of which would be real bugs if surfaced by the proptest strategy. Two additional problems:
  - `let Ok(_) else { return Ok(()) }` trivially passes the case without recording it. `prop_assume!` would track the rejection rate and warn when it gets high; the current pattern hides degeneracy if `DidNotConverge` ever becomes the common outcome.
  - The strategy at [posterior_engine.rs:2456](../../../src/var_calling/posterior_engine.rs#L2456) draws likelihoods uniformly from `-50.0..0.0` for every genotype cell, which routinely produces conflicting per-genotype likelihoods near the edges of the range. Without a measured rejection-rate, the proptests can degenerate to no-ops.
- **Why it matters:** Errors rule: "No discarded `Result` ... without a comment naming the *specific* failure mode being discarded and why it is safe to ignore." The comment names `DidNotConverge` only. Reliability rule: a property test that can silently pass without testing the property is a regression hazard. The impl report admits this as in-scope hardening surfaced by the validation gate; finishing the hardening properly is in scope.
- **Suggested fix:**
  ```rust
  let result = try_single(record);
  prop_assume!(matches!(result, Ok(_) | Err(PosteriorEngineError::DidNotConverge { .. })));
  let pr = match result {
      Ok(pr) => pr,
      Err(_) => return Ok(()), // DidNotConverge — documented
  };
  ```
  And verify the rejection rate: `PROPTEST_CASES=512 ./scripts/dev.sh cargo test -- posteriors_sum_to_one_at_any_ploidy_and_n_alleles` and read the rejection-rate line. If it exceeds proptest's default threshold (typically 50 %), narrow the strategy's likelihood range away from `-50.0` or stratify by ploidy.

**M4: [src/var_calling/posterior_engine.rs:734](../../../src/var_calling/posterior_engine.rs#L734) — Mixture-clamp magic `1e-12` / `0.5` duplicates `MIN_BASE_ERROR` / `MAX_BASE_ERROR` with no shared name**
- **Categories:** defaults, naming, errors (convergent across 3).
- **Confidence:** High
- **Problem:** The Stage 6 mixture branch writes `mean_err[a] = per_read_log_err.exp().clamp(1e-12, 0.5);` inline. The side-pass at [contamination_estimation.rs:746](../../../src/var_calling/contamination_estimation.rs#L746) uses the same numbers behind the named `MIN_BASE_ERROR` / `MAX_BASE_ERROR` constants (lines 969 / 974) with a doc comment justifying each. The two sites must move together — the floor/ceiling are part of one model assumption — but here the literals are anonymous, so a future tuning of one site silently desyncs from the other.
- **Why it matters:** Defaults rule "No magic numbers": every behavioral default lives as a named `pub const` with a doc explaining choice and source. Naming rule "the literal would need to be kept in sync with anything else (another code site)" is hit. The same critique applies to the silent-recovery concern in M2(a).
- **Suggested fix:** Promote the side-pass constants to `pub` and reuse them at the engine site:
  ```rust
  // contamination_estimation.rs (promote)
  pub const MIN_BASE_ERROR: f64 = 1e-12;
  pub const MAX_BASE_ERROR: f64 = 0.5;

  // posterior_engine.rs:734
  use crate::var_calling::contamination_estimation::{MIN_BASE_ERROR, MAX_BASE_ERROR};
  mean_err[a] = per_read_log_err.exp().clamp(MIN_BASE_ERROR, MAX_BASE_ERROR);
  ```

**M5: [src/var_calling/contamination_estimation.rs:288](../../../src/var_calling/contamination_estimation.rs#L288) — `ContaminationEstimates::zero` enforces its invariant only in debug builds**
- **Categories:** errors, idiomatic, extras (convergent).
- **Confidence:** High
- **Problem:** The constructor uses `debug_assert_eq!(sample_to_batch.len(), n_samples);` and is silent in release. There is no check that every entry in `sample_to_batch` is `< n_batches`. A release-build caller with a malformed `sample_to_batch` silently builds an estimate whose `c_s_per_sample.len() == n_samples` but whose `sample_to_batch.len()` may differ — the new `q_b_for_sample` accessor and the Stage 6 `compute_mixture_log_likelihoods` rely on both being length-aligned. `from_user_supplied` (B1) has the same shape and the same missing check.
- **Why it matters:** This is the documented "no-contamination" constructor (the CLI takes this path when `--contamination-batches` is not supplied), so it lives on the default code path. Release-only divergence from debug behaviour is exactly the failure mode `debug_assert!` is not appropriate for; constructor validation runs once per cohort, not per-read, so the cost of a runtime check is negligible.
- **Suggested fix:** Replace with a runtime check, ideally a fallible constructor that returns `Result`. At minimum:
  ```rust
  pub fn zero(n_samples: usize, sample_to_batch: Vec<usize>, n_batches: usize) -> Self {
      assert_eq!(sample_to_batch.len(), n_samples,
          "sample_to_batch length {} does not match n_samples {}",
          sample_to_batch.len(), n_samples);
      if let Some(&max_b) = sample_to_batch.iter().max() {
          assert!(max_b < n_batches, "sample_to_batch contains batch_idx {max_b} >= n_batches {n_batches}");
      }
      Self { /* ... */ }
  }
  ```

**M6: [src/var_calling/contamination_estimation.rs:392](../../../src/var_calling/contamination_estimation.rs#L392) — `BadInput(String)` is a stringly-typed catch-all**
- **Categories:** errors.
- **Confidence:** High
- **Problem:** The variant funnels three structurally distinct invariant violations through one string. The three call sites (lines 421, 430, 436) format different fields into the message; callers cannot match on the specific failure mode without parsing English. The unit tests at lines 1394–1443 already only `matches!(... BadInput(_))` because the inner string is not structurally addressable.
- **Why it matters:** A downstream caller that wants different remediations for "sample_to_batch length wrong" vs "batch_idx out of range" vs "block_size = 0" cannot do so. The errors rule explicitly calls out `Other(String)` / mechanism-named variants.
- **Suggested fix:** Split into structured variants:
  ```rust
  #[error("sample_to_batch length {got} does not match n_samples {expected}")]
  SampleToBatchLengthMismatch { expected: usize, got: usize },
  #[error("sample_to_batch contains batch_idx {batch_idx} >= n_batches {n_batches}")]
  BatchIndexOutOfRange { batch_idx: usize, n_batches: usize },
  #[error("block_size must be > 0")]
  ZeroBlockSize,
  ```

**M7: [src/var_calling/contamination_estimation.rs:357](../../../src/var_calling/contamination_estimation.rs#L357) — `Upstream(#[from] PerPositionMergerError)` erases side-pass context at every `?`**
- **Categories:** errors.
- **Confidence:** High
- **Problem:** The single `?` site at line 446 (`let pileups = item?;`) silently drops the side-pass's own state (position number, per-block site count, partial-progress flag). When the side-pass fails mid-stream the user sees only the inner merger error and cannot tell whether the failure happened at block 1 or block 100, or whether partial estimates were already accumulated. The recommendation machinery does not fire for the `Upstream` path either.
- **Why it matters:** A transient reader hiccup at block 99/100 looks indistinguishable from a fatal header mismatch on the first record.
- **Suggested fix:** Replace the blanket `#[from]` with a structured variant and convert explicitly at the `?`:
  ```rust
  #[error("upstream merger failed after {sites_processed} informative sites")]
  Upstream {
      sites_processed: u32,
      #[source]
      source: PerPositionMergerError,
  },
  // at the ? site:
  let pileups = item.map_err(|source| ContaminationEstimationError::Upstream {
      sites_processed, source,
  })?;
  ```

**M8: [src/var_calling/contamination_estimation.rs:322,328](../../../src/var_calling/contamination_estimation.rs#L322) — `effective_c_s` and `q_b_for_sample` panic on out-of-range index with no typed error**
- **Categories:** errors, reliability (convergent).
- **Confidence:** High
- **Problem:** Both accessors index directly without bounds-checking and without a `// PANIC-FREE:` annotation. Stage 6 calls both inside the per-record loop ([posterior_engine.rs:710](../../../src/var_calling/posterior_engine.rs#L710) and L722); the engine's cohort-size guard at [posterior_engine.rs:862](../../../src/var_calling/posterior_engine.rs#L862) checks only `c_s_per_sample.len() != n_samples`, not that `sample_to_batch.len()` matches or that every `sample_to_batch[s] < q_b_per_batch.len()`. A `ContaminationEstimates` built via the unvalidated `from_user_supplied` (B1) slips past the guard and panics deep in the EM hot loop with `index out of bounds` and no `locus`.
- **Why it matters:** Direct slice indexing without a documented invariant is morally equivalent to `.unwrap()`. The project memory note "no logs — promote to typed errors" applies.
- **Suggested fix:** (a) Once B1 lands, document the invariant with `// PANIC-FREE:` comments referencing the validating constructors; or (b) tighten the Stage 6 guard at [posterior_engine.rs:862](../../../src/var_calling/posterior_engine.rs#L862) to also validate `sample_to_batch.len() == n_samples` and `sample_to_batch.iter().all(|&b| b < q_b_per_batch.len())`, returning typed `ContaminationCohortSampleToBatchMismatch` / `ContaminationBatchIndexOutOfRange` with `locus`.

**M9: [src/var_calling/contamination_estimation.rs:867](../../../src/var_calling/contamination_estimation.rs#L867) — `.max(2)` silently overrides user-set `min_batch_size_for_contamination`**
- **Categories:** defaults. **Convergent:** extras (cross-category).
- **Confidence:** High
- **Problem:** `n >= config.min_batch_size_for_contamination.max(2)` silently coerces a user-set floor of `1` to `2`. The struct-field doc at lines 199–204 notes the rule, but the constant doc at lines 80–82, `with_project_defaults` doc, and `ContaminationEstimates` provenance do not surface it; the effective floor is unrecoverable from a running instance.
- **Why it matters:** Defaults rule "Inspectable at runtime" — every effective default must be recoverable from a running instance. A user who reads `DEFAULT_MIN_BATCH_SIZE_FOR_CONTAMINATION = 5` and sets `min_batch_size_for_contamination = 1` to disable the floor silently gets a floor of 2; no error, no warning, no inspectable value.
- **Suggested fix:** Validate at the entry gate:
  ```rust
  if config.min_batch_size_for_contamination < 2 {
      return Err(ContaminationEstimationError::BadInput(
          "min_batch_size_for_contamination must be >= 2 (singletons are \
           always floored — c_s and q_b are unidentifiable in a batch of one)"
              .to_string(),
      ));
  }
  ```
  Then drop the `.max(2)`. Pair with the M6 structural-variant fix.

**M10: [src/var_calling/contamination_estimation.rs:89](../../../src/var_calling/contamination_estimation.rs#L89) — `DEFAULT_CS_INIT = 0.02` is behavior-affecting but not surfaced on the public config**
- **Categories:** defaults. **Naming sub-finding:** `CS` vs `c_s` spelling inconsistency (see naming notes).
- **Confidence:** High
- **Problem:** `DEFAULT_CS_INIT` seeds `c_s_frozen` and `c_s_last_snapshot` in `OnlineEmState::new` (lines 552, 554), determining the first block's γ_r computation and therefore the entire online-EM trajectory. It is not a field on `ContaminationEstimationConfig`, so a caller cannot override or inspect it, and the field table in the config's struct docs does not mention it. `OnlineEmState::new` already takes `_config: &ContaminationEstimationConfig` but does not read it.

  Naming sub-finding: every other reference uses `c_s` (`c_s_per_sample`, `c_s_frozen`, etc.), but the default constant uses the acronym `CS`. Rename to `DEFAULT_C_S_INIT` for consistency once the field promotion lands.
- **Why it matters:** A caller using `with_project_defaults()` has no way to see, from the type alone, that 0.02 (not 0.0) seeds the run. Silent priors on EM init can mask convergence bugs.
- **Suggested fix:**
  ```rust
  pub struct ContaminationEstimationConfig {
      // ...
      /// Initial per-sample c_s seeded into the first block's frozen
      /// snapshot. Defaults to [`DEFAULT_C_S_INIT`].
      pub c_s_init: f64,
  }
  // with_project_defaults: c_s_init: DEFAULT_C_S_INIT,
  // OnlineEmState::new: c_s_frozen: vec![config.c_s_init; n_samples],
  ```

**M11: [src/var_calling/contamination_estimation.rs:553](../../../src/var_calling/contamination_estimation.rs#L553) — `q_b` init uniform prior `1.0 / N_ALLELE_CLASSES` is a hidden default**
- **Categories:** defaults.
- **Confidence:** High
- **Problem:** `OnlineEmState::new` seeds `q_b_frozen` with `[1.0 / N_ALLELE_CLASSES as f64; N_ALLELE_CLASSES]` inline, with no named constant, no doc justifying the choice, no surface on the public config. Determines the first block's mixture denominator before the first refresh.
- **Why it matters:** Reader cannot tell whether the literal is a placeholder or a deliberate prior; not tunable; not inspectable.
- **Suggested fix:**
  ```rust
  /// Initial per-batch q_b seeded into the first block's frozen
  /// snapshot — a uniform-Dirichlet-mean prior so the first block's
  /// mixture denominator is well-defined for every allele class
  /// before any reads contribute.
  pub const DEFAULT_Q_B_INIT_PER_CLASS: f64 = 1.0 / N_ALLELE_CLASSES as f64;
  // OnlineEmState::new:
  q_b_frozen: vec![[DEFAULT_Q_B_INIT_PER_CLASS; N_ALLELE_CLASSES]; n_batches],
  ```

**M12: [src/var_calling/contamination_estimation.rs:214](../../../src/var_calling/contamination_estimation.rs#L214) — `Default for ContaminationEstimationConfig` is the exact anti-pattern the defaults rule forbids**
- **Categories:** defaults.
- **Confidence:** High
- **Problem:** `impl Default { fn default() -> Self { Self::with_project_defaults() } }` makes the project defaults reachable via `Default::default()`, hiding every behaviorally-significant value behind a zero-argument constructor. The integration test at [tests/contamination_estimation_integration.rs:125](../../../tests/contamination_estimation_integration.rs#L125) already uses `with_project_defaults()` — the `Default` impl is at odds with that convention.
- **Why it matters:** Once the impl is available, future call sites will write `ContaminationEstimationConfig::default()` or `..Default::default()`, defeating the explicit-defaults contract `with_project_defaults` was named to enforce.
- **Suggested fix:** Remove the `Default` impl entirely.
  ```rust
  // Delete the entire impl Default block at lines 214–218.
  ```
  Force every caller to choose `with_project_defaults()` (visible name) or a hand-built literal.

**M13: 14 sites — Test struct literals use `..Default::default()` / `..with_project_defaults()` tails**
- **Categories:** refactor_safety.
- **Confidence:** High
- **Problem:** Two new posterior-engine tests use `..Default::default()` ([posterior_engine.rs:1950](../../../src/var_calling/posterior_engine.rs#L1950), [posterior_engine.rs:1970](../../../src/var_calling/posterior_engine.rs#L1970)), and 11 new in-module side-pass tests use `..ContaminationEstimationConfig::with_project_defaults()` ([contamination_estimation.rs:1165, 1189, 1216, 1249, 1274, 1303, 1328, 1349, 1436, 1469, 1497, 1525](../../../src/var_calling/contamination_estimation.rs#L1165)). When `PosteriorEngineConfig` or `ContaminationEstimationConfig` grow a new field — the spec already anticipates per-class minimum-coverage knobs — every site silently picks up the new default rather than forcing the test author to decide. Some of these are sensitivity tests (`online_em_recovers_three_percent_contamination` asserts ±0.02 on `c_s ≈ 0.03`) where a new default-value field could quietly shift the asserted number.
- **Why it matters:** Project rule: "Make the compiler flag refactors — no `..Default::default()`, exhaustive destructures." The test suite must fail to compile until each test considers the new field.
- **Suggested fix:** Replace the tails with explicit per-field spellings, or share a `test_cfg(...)` helper that names every field:
  ```rust
  fn test_cfg(stopping_mode: StoppingMode, block_size: u32, min_batch: u32) -> ContaminationEstimationConfig {
      ContaminationEstimationConfig {
          stopping_mode, block_size,
          min_depth: DEFAULT_MIN_DEPTH,
          min_major_fraction: DEFAULT_MIN_MAJOR_FRACTION,
          min_cohort_minor_count: DEFAULT_MIN_COHORT_MINOR_COUNT,
          min_cohort_minor_fraction: DEFAULT_MIN_COHORT_MINOR_FRACTION,
          min_batch_size_for_contamination: min_batch,
          ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT,
          snp_alt_pseudocount: DEFAULT_SNP_ALT_PSEUDOCOUNT,
          indel_alt_pseudocount: DEFAULT_INDEL_ALT_PSEUDOCOUNT,
      }
  }
  ```
  The integration-test `fixed_n_config` at [tests/contamination_estimation_integration.rs:124-130](../../../tests/contamination_estimation_integration.rs#L124-L130) already names every override explicitly — it's the right shape. The two `let mut config = Default::default(); config.contamination = ...` sites at [tests/contamination_estimation_integration.rs:199](../../../tests/contamination_estimation_integration.rs#L199) and L277 leak the same silent-default-pickup risk and should be replaced with literals naming every field of `PosteriorEngineConfig`.

**M14: [src/var_calling/contamination_estimation.rs:520](../../../src/var_calling/contamination_estimation.rs#L520) — `OnlineEmState` god-struct (3 disjoint concerns, 10 fields)**
- **Categories:** smells.
- **Confidence:** High
- **Problem:** `OnlineEmState` bundles three disjoint concerns:
  - Running sufficient stats: `s_per_sample`, `n_per_sample`, `s_per_batch`, `n_per_batch` (mutated by `apply_observation`).
  - Frozen parameter snapshots: `c_s_frozen`, `q_b_frozen` (read by `apply_observation`, rewritten by `refresh_parameters_and_snapshot`).
  - Stability bookkeeping: `c_s_last_snapshot`, `last_block_deltas`, `last_block_max_delta`, `consecutive_within_tolerance` (only `refresh_parameters_and_snapshot` and the block-boundary branch of `estimate_contamination`; `apply_observation` never reads).
  `apply_observation` operates only on groups 1+2 — the stability fields are dead state from its perspective. Classic "different methods operate on disjoint subsets of fields" trigger plus "more than ~7 fields with no clear grouping."
- **Why it matters:** The current layout makes the EM math impossible to read in isolation from the stopping-criterion bookkeeping. Splitting moves each piece next to the method that operates on it, lets the stability tracker grow (adaptive tolerance, per-sample convergence) without touching the EM core, and gives `apply_observation` a tighter signature that makes its non-touching of stability a compile-time guarantee.
- **Suggested fix:**
  ```rust
  struct EmRunningStats {
      s_per_sample: Vec<f64>,
      n_per_sample: Vec<u64>,
      s_per_batch: Vec<[f64; N_ALLELE_CLASSES]>,
      n_per_batch: Vec<u64>,
  }
  struct EmFrozenParameters {
      c_s: Vec<f64>,
      q_b: Vec<[f64; N_ALLELE_CLASSES]>,
  }
  struct StabilityTracker {
      c_s_last_snapshot: Vec<f64>,
      last_block_deltas: Vec<f64>,
      last_block_max_delta: f64,
      consecutive_within_tolerance: u32,
  }
  struct OnlineEmState {
      stats: EmRunningStats,
      frozen: EmFrozenParameters,
      stability: StabilityTracker,
  }
  ```
  Follow-up: lift the per-block-boundary branch from `estimate_contamination` onto `impl OnlineEmState` (`fn on_block_boundary(&mut self, config) -> BoundaryOutcome`); the main loop becomes one decision per line.

**M15: [src/var_calling/contamination_estimation.rs:585](../../../src/var_calling/contamination_estimation.rs#L585) — `step_1a_cohort_filter` allocates per-position on the cohort-CLI hot path**
- **Categories:** extras. **Convergent:** smells, errors, reliability, idiomatic (cross-category notes).
- **Confidence:** High
- **Problem:** `step_1a_cohort_filter` runs once per cohort position and allocates a fresh `Vec<(Vec<u8>, u32)>` with capacity 4, then `obs.seq.clone()`s each novel allele's bytes. At cohort scale this is one outer `Vec` + N small `Vec<u8>`s per position, every position, for the entire side-pass.
- **Why it matters:** The spec calls this out as a critical-path scan ("CPU is dominated by the multi-way per-position scan over `.psp` files", spec §"Cost and parallelism"). Even with convergence-mode early-stop, every position before early-stop pays the cost. SNP byte-clones are small but go through the system allocator; indel clones can be longer.
- **Suggested fix:** Hoist the scratch vector outside the loop and use borrowed `&[u8]` keys (the `obs` references live as long as `pileups`):
  ```rust
  // In estimate_contamination(): allocate once, outside the loop.
  let mut cohort_alleles: Vec<(&[u8], u32)> = Vec::with_capacity(4);
  ```
  Pass `&mut cohort_alleles` into `step_1a_cohort_filter` and `cohort_alleles.clear()` at the top. If `find()` over a small vector ever becomes interesting, `ArrayVec<[(&[u8], u32); 8]>` is a drop-in replacement.

**M16: [src/var_calling/contamination_estimation.rs:721](../../../src/var_calling/contamination_estimation.rs#L721) — Integer saturation/overflow inconsistency across three sites in the EM hot path**
- **Categories:** reliability. **Convergent:** extras (Minor).
- **Confidence:** Medium
- **Problem:** Three sites use three different policies for the same kind of operation:
  - `step_1a_cohort_filter` line 590: `entry.1 = entry.1.saturating_add(obs.support.num_obs)` (cohort sum, saturating).
  - `step_1b_sample_filter` line 673: `counts_by_class[class as usize].saturating_add(obs.support.num_obs)` (per-class, saturating).
  - `apply_observation` line 731: `let total_reads: u32 = obs.counts_by_class.iter().sum();` (non-saturating; panics in debug, wraps in release on overflow).
  The depth sum at [step_1b_sample_filter:646-651](../../../src/var_calling/contamination_estimation.rs#L646-L651) uses `Iterator::sum`, also non-saturating. A saturated upstream `counts_by_class[a] = u32::MAX` would overflow before reaching `apply_observation`. None of the existing tests exercise high-coverage or large-cohort inputs (synthetic streams use `n_reads ∈ 20..80`, `n_samples ∈ 2..=4`).
- **Why it matters:** Realistic per-sample depths fit in `u32`, so this is unlikely to fire — but if `.psp` ever permits a corrupted aggregate with `num_obs` near `u32::MAX`, the silent saturation upstream and the panic-or-wrap downstream silently produce different bugs. Pick one policy.
- **Suggested fix:** Either (a) move all three sites to saturating arithmetic and accumulate in `u64` for the downstream sum:
  ```rust
  // apply_observation:
  let total_reads: u64 = obs.counts_by_class.iter().map(|&n| u64::from(n)).sum();
  ```
  Or (b) move all three to `checked_add(...).expect("psp count overflowed u32")` so the data condition surfaces as a panic at the source. Add a regression test `apply_observation_does_not_overflow_at_u32_max_reads` per the missing-tests list.

### Minor

- **Mi1: [src/var_calling/contamination_estimation.rs:741](../../../src/var_calling/contamination_estimation.rs#L741) — Dead `if total_reads > 0` after early return.** Convergent: errors, idiomatic, reliability, smells. The function returns at line 733 when `total_reads == 0`, so the `else { 0.0 }` branch at line 745 is unreachable. Drop the `if`: `let mean_log_error = total_q_sum / total_reads as f64;`.

- **Mi2: [src/var_calling/posterior_engine.rs:674](../../../src/var_calling/posterior_engine.rs#L674) — `#[allow(clippy::too_many_arguments)]` on `compute_mixture_log_likelihoods` lacks the 3-part justification.** Convergent: idiomatic, smells. Of the 9 parameters, `n_samples`, `n_genotypes`, and `n_alleles` are derivable from the other inputs (the function already calls `genotype_order` internally). Either derive sizes inside the function and drop the `allow`, or annotate the `allow` with why it stays (caller-side destructure) and what would trigger removal.

- **Mi3: [src/var_calling/contamination_estimation.rs:567](../../../src/var_calling/contamination_estimation.rs#L567) — `CohortFilterOutcome` is a single-bool wrapper.** Convergent: naming, idiomatic, smells. Return `bool` directly, or add a second field that justifies the struct (e.g. `cohort_minor_fraction` for the TSV report).

- **Mi4: [src/var_calling/contamination_estimation.rs:338](../../../src/var_calling/contamination_estimation.rs#L338) — `NextRunRecommendation` is half-structured + missing `#[non_exhaustive]` + missing field docs.** Convergent: errors, idiomatic, extras. The struct carries an actionable `suggested_num_sites` next to a free-form English `message`; tests already only assert that `suggested_num_sites.is_some() || !message.is_empty()`. Add `#[non_exhaustive]`, add `///` field docs, and consider a structured `RecommendationAction` enum to make the actions programmatically addressable.

- **Mi5: [src/var_calling/posterior_engine.rs:727](../../../src/var_calling/posterior_engine.rs#L727) — `compute_mixture_log_likelihoods` slice indexes lack precondition assertions.** The function indexes `scalars[sample_idx * n_alleles..(sample_idx+1) * n_alleles]` and `fallback_log_likelihoods[sample_idx * n_genotypes..(sample_idx+1) * n_genotypes]` without checking the slices are large enough. Internal callers validate `MergedRecord` shape upstream, but this helper is private and could be reused later from a non-validating path. Add `debug_assert_eq!(scalars.len(), n_samples * n_alleles)` at the top + a `// PANIC-FREE:` doc comment.

- **Mi6: [src/var_calling/contamination_estimation.rs:267,328](../../../src/var_calling/contamination_estimation.rs#L267) — `ContaminationEstimates` mixes public fields with effective-value accessors.** Field docs are `pub` but the docstring on `effective_c_s` says "callers should not branch on the `Option` themselves." The contract is unenforceable while the field is `pub`. Make the fields `pub(crate)` (preferred) or drop the accessors. The two existing call sites inside Stage 6 already go through the accessors.

- **Mi7: [src/var_calling/contamination_estimation.rs:906,944](../../../src/var_calling/contamination_estimation.rs#L906) — `recommend_for_non_convergence` / `recommend_for_insufficient_sites` should be `NextRunRecommendation` constructors.** Move into `impl NextRunRecommendation` as `for_non_convergence` / `for_insufficient_sites`. Discoverable via the type's doc page; same call-site shape.

- **Mi8: [src/var_calling/contamination_estimation.rs:1123](../../../src/var_calling/contamination_estimation.rs#L1123) — `synth_two_sample_contam_stream` builds three samples.** The integration-test copy at [tests/contamination_estimation_integration.rs:54](../../../tests/contamination_estimation_integration.rs#L54) is correctly named `synth_three_sample_stream`. Rename the in-module copy to match.

- **Mi9: Test fixture builders duplicated between in-module and integration tests** (`obs`, `pileup_record`/`rec`, `position`/`pos_pileups`, `synth_*`). Two divergences are already visible (the q_sum formula; the `s0/s1/s2` order). Consolidate into a `#[cfg(test)]` shared module or re-export under `pub use #[cfg(test)]`.

- **Mi10: [src/var_calling/contamination_estimation.rs:660](../../../src/var_calling/contamination_estimation.rs#L660) — `step_1b_sample_filter` divides by `depth` without enforcing `min_depth >= 1`.** A future config that exposes `min_depth = 0` without further validation could enter a state where a record with empty allele observations (alleles non-empty but every `num_obs == 0`) returns `Some(...)` with NaN scalars. Add `config.min_depth.max(1)` or validate at the entry point.

- **Mi11: [src/var_calling/contamination_estimation.rs:906](../../../src/var_calling/contamination_estimation.rs#L906) — `recommend_for_non_convergence` has no regression tests for the 10× cap, the zero-delta branch, or the infinite-delta branch.** Add `recommend_for_non_convergence_returns_capped_suggestion`, `..._returns_none_on_zero_delta`, `..._returns_finite_suggestion_on_infinite_delta` (specs in §8).

- **Mi12: [src/var_calling/posterior_engine.rs:655](../../../src/var_calling/posterior_engine.rs#L655) — `map_to_contam_class` has no direct test.** A future allele class added without remembering to update the bridge would silently route to `None` (treat-as-impossible), corrupting the mixture. Add `map_to_contam_class_round_trips_all_three_non_compound_classes` and `..._returns_none_for_compound`.

- **Mi13: [src/var_calling/posterior_engine.rs:698,727](../../../src/var_calling/posterior_engine.rs#L698) — `row` reused for two unrelated row-slices in `compute_mixture_log_likelihoods`.** Once for a genotype's per-allele copy counts, once for a sample's per-allele scalars. Rename:
  ```rust
  let genotype_copy_counts = &mut genotype_allele_counts[g_idx * n_alleles..(g_idx + 1) * n_alleles];
  let sample_allele_scalars = &scalars[sample_idx * n_alleles..(sample_idx + 1) * n_alleles];
  ```

- **Mi14: [src/var_calling/posterior_engine.rs:707,716](../../../src/var_calling/posterior_engine.rs#L707) — `out` is a generic noun for the returned log-likelihood matrix.** The function's docstring is explicit that the value is "a fresh `n_samples × n_genotypes` row-major `log_likelihoods` table"; the bare `out` discards that meaning at every use site. Rename to `mixture_log_likelihoods`.

- **Mi15: [src/var_calling/posterior_engine.rs:705](../../../src/var_calling/posterior_engine.rs#L705) — `non_genotype_classes` describes a structural role, not the divisor it holds.** Better:
  ```rust
  let other_allele_error_denom = (n_alleles as f64 - 1.0).max(1.0);
  ```

- **Mi16: [tests/contamination_estimation_integration.rs:124](../../../tests/contamination_estimation_integration.rs#L124) — `fixed_n_config` hides intent (also mutates `block_size` and `min_batch_size_for_contamination`).** Rename to `small_cohort_fixed_n_config` (carries every relaxation in the name) or split into single-purpose helpers.

- **Mi17: [src/var_calling/contamination_estimation.rs:851](../../../src/var_calling/contamination_estimation.rs#L851) — `finalise` uncovered against `min_batch_size_for_contamination = 1`.** A future refactor that drops the `.max(2)` would silently let singleton batches keep their unidentifiable estimates. Add `singleton_batch_with_floor_of_one_still_returns_none`. (Becomes moot after M9 lands at the validation gate.)

- **Mi18: [src/var_calling/contamination_estimation.rs:411](../../../src/var_calling/contamination_estimation.rs#L411) — `estimate_contamination`'s "did-not-tick if observed_count == 0" semantics are not tested.** A future refactor flipping to "count any Step-1a pass" would silently change the spec contract. Add `site_counter_does_not_tick_when_step_1b_finds_no_samples` (spec in §8).

- **Mi19: [src/var_calling/contamination_estimation.rs:64,69,73,77,82](../../../src/var_calling/contamination_estimation.rs#L64) — Five `DEFAULT_*` constants cite "VerifyBamID" or "Spec" without page/section/version.** The pseudocount constant at line 96 is the template ("Mirrors the posterior engine's REF pseudocount … per the spec's […] line"). Add a spec-section reference or a VerifyBamID release/flag-name to each.

- **Mi20: [src/var_calling/contamination_estimation.rs:794,811](../../../src/var_calling/contamination_estimation.rs#L794) — Redundant `.clone()` of `c_s_new` and `per_sample_deltas` in `refresh_parameters_and_snapshot`.** Both are returned by value or moved one line later. Restructure so the caller writes back into `state`, removing the per-call clone.

- **Mi21: [src/var_calling/contamination_estimation.rs:410](../../../src/var_calling/contamination_estimation.rs#L410) — `estimate_contamination` is a long flat state machine (~70 lines, ~9 decision points, 4 nesting levels).** Borderline on the function-length rule. Pairs with the M14 god-struct refactor — extracting the block-boundary branch onto `OnlineEmState::on_block_boundary` shrinks the main loop to a clear dispatch.

- **Mi22: [tests/contamination_estimation_integration.rs:199,277](../../../tests/contamination_estimation_integration.rs#L199) — `let mut config = PosteriorEngineConfig::default(); config.contamination = Some(...);` pattern leaks the same silent-default risk as M13.** Spell every field of `PosteriorEngineConfig` in one literal, or share a helper.

- **Mi23: [src/var_calling/contamination_estimation.rs:546](../../../src/var_calling/contamination_estimation.rs#L546) — `OnlineEmState::new` takes `_config: &ContaminationEstimationConfig` it never reads.** Convergent: idiomatic, smells, reliability, naming (cross-category), extras. Drop the parameter (resolves the nit at the call site) or wire `c_s_init`/`q_b_init` through it (resolves M10 + M11 in one move).

### Nits

Grouped — none enumerated:

- [contamination_estimation.rs:443](../../../src/var_calling/contamination_estimation.rs#L443): `let mut at_block_boundary;` declared outside the loop without initializer; can be a normal `let at_block_boundary = …` inside the loop.
- [contamination_estimation.rs:927-932,946-949](../../../src/var_calling/contamination_estimation.rs#L927): `match Option { Some => format!, None => String::new() }` reads cleaner as `.map_or(String::new(), |n| format!(...))`.
- [contamination_estimation.rs:1453-1454](../../../src/var_calling/contamination_estimation.rs#L1453): dead `.clamp(1, 3).min(n_samples)` in a proptest where the proptest input is already in `2usize..4`.
- [contamination_estimation.rs:920](../../../src/var_calling/contamination_estimation.rs#L920): redundant `.min(u32::MAX as f64)` since the `as` cast became saturating in Rust 1.45.
- [contamination_estimation.rs:678](../../../src/var_calling/contamination_estimation.rs#L678): the `sample_idx: 0, // filled in by the caller` placeholder comment becomes redundant once M1 is applied.
- [contamination_estimation.rs:130-138](../../../src/var_calling/contamination_estimation.rs#L130): `AlleleClass` variants lack per-variant `///` docs.
- [contamination_estimation.rs:338-341](../../../src/var_calling/contamination_estimation.rs#L338): `NextRunRecommendation`'s public fields lack `///` docs (rolls into Mi4).
- [contamination_estimation.rs:475-476](../../../src/var_calling/contamination_estimation.rs#L475): redundant assignments to `state.last_block_deltas` / `state.last_block_max_delta` already set by `refresh_parameters_and_snapshot`.
- [contamination_estimation.rs:283](../../../src/var_calling/contamination_estimation.rs#L283): `ContaminationEstimates::zero` doc would benefit from one example call given the three positional `usize`/`Vec<usize>`/`usize` arguments are easy to swap.
- [contamination_estimation.rs:794,811](../../../src/var_calling/contamination_estimation.rs#L794): brief comment ("init from the previous snapshot so unobserved samples keep their last value") to prevent a future reader from "optimising" the clones away — once Mi20 lands the comment can drop.
- [contamination_estimation.rs:585](../../../src/var_calling/contamination_estimation.rs#L585): positional-tuple `cohort_alleles[*].1` accesses could use a 2-field struct for readability — moot once M15's borrow refactor lands.
- [contamination_estimation.rs:1110](../../../src/var_calling/contamination_estimation.rs#L1110): `obs` binding in tests shadows the module-local `fn obs` fixture and the loop variable `obs` in `step_1b_sample_filter`; rename to `informative_obs`.
- [contamination_estimation.rs:1397,1413,1429](../../../src/var_calling/contamination_estimation.rs#L1397): the three `BadInput` tests assert only `matches!(... BadInput(_))` without checking message content. After M6's structural-variant fix the match becomes specific; until then, a `contains("sample_to_batch")` etc. would document the operator-facing string.
- [posterior_engine.rs:1602](../../../src/var_calling/posterior_engine.rs#L1602): the `try_single`-helper comment ("`DidNotConverge` … documented engine outcome on adversarial likelihood matrices, not a bug") should live on the rustdoc of `PosteriorEngineError::DidNotConverge`, not in a test helper.

## 7. Out of scope observations

- **[per_group_merger.rs:1786-1796](../../../src/var_calling/per_group_merger.rs#L1786-L1796)** — pre-existing `clippy::needless_range_loop` lines were resolved in this commit as a side-effect of clearing the `-D warnings` gate. Suggested follow-up: PROJECT_STATUS already noted these as a separate "clippy-cleanup pass" follow-up; that note is now resolved by this commit, no further action needed.
- **[benches/var_calling_perf.rs:~123,~365](../../../benches/var_calling_perf.rs)** — same as above; resolved in this commit.
- **[src/var_calling/posterior_engine.rs:2456](../../../src/var_calling/posterior_engine.rs#L2456)** — pre-existing `proptest_record_strategy` draws likelihoods uniformly from `-50.0..0.0` for every genotype cell. This strategy contributes to M3 by routinely producing `DidNotConverge` cases. Tightening the strategy (narrow the range, stratify by ploidy) is a separate follow-up not blocking this commit; it's the upstream cause of the `try_single` hardening having to exist at all.

## 8. Missing tests to add now

Grouped by function under test. The reliability sub-agent's challenge-tests pass feeds this section.

### `compute_mixture_log_likelihoods` ([posterior_engine.rs:675](../../../src/var_calling/posterior_engine.rs#L675))

- **`compute_mixture_log_likelihoods_matches_hand_calculation_at_c_s_5_percent`** — single sample, single batch, ploidy = 2, two alleles (`b"A"`, `b"C"`), `c_s = 0.05`, `q_b = [0.95, 0.05, 0.0]`, scalars `(num_obs=50, q_sum=50·ln(0.001))` for REF and default for ALT. Hand calc: `genotype RR` (k=2) → `p_own = 0.999`, `mix = 0.95·0.999 + 0.05·0.95 ≈ 0.99655`, `ll(RR) = 50·ln(0.99655) ≈ -0.1729`. Catches: swapped `(1-c_s)`/`c_s`, wrong q_b indexing, wrong `non_genotype_classes` denominator, mis-built `genotype_allele_counts`.
- **`compute_mixture_log_likelihoods_falls_back_to_precomputed_when_c_s_is_zero`** — sample with `c_s_per_sample[0] = None` (`effective_c_s = 0`). Provide recognisable fallback `vec![1.0, 2.0, 3.0]`; assert `out == fallback`. Catches: off-by-one in the fallback stride copy.
- **`compute_mixture_log_likelihoods_returns_non_finite_error_when_mix_is_zero`** — sample with `c_s = 1.0`, `q_b = [1.0, 0.0, 0.0]`, record where allele 1 (SNP-alt) has `num_obs > 0`. Then `p_contam = 0` for allele 1, `(1-c_s) = 0`, `mix = 0`, error must fire. Catches: silent NaN propagation if the guard is dropped.
- **`compute_mixture_log_likelihoods_maps_compound_allele_to_zero_q_contam`** — 3 alleles with `alleles[2].is_compound = true`. Two runs: one with `q_b = [1/3,1/3,1/3]` and one with `q_b = [0.0, 0.0, 1.0]`. The compound-allele likelihood term must be identical between runs (both yield `p_contam = 0`). Catches: M12-style regression mapping `CompoundAlt → SnpAlt`.
- **`compute_mixture_log_likelihoods_is_invariant_under_q_b_normalisation_at_zero_reads`** — sample with `n_obs[a] == 0` for every allele. Catches: a stray `+= mix.ln()` outside the `if n_a == 0 { continue }` branch.

### `map_to_contam_class` ([posterior_engine.rs:655](../../../src/var_calling/posterior_engine.rs#L655))

- **`map_to_contam_class_round_trips_all_three_non_compound_classes`** — assert each `AlleleClass::{Ref, SnpAlt, IndelAlt}` maps to the matching `ContamAlleleClass::*`. Catches: future allele class added without updating this fn (silent `None`).
- **`map_to_contam_class_returns_none_for_compound`** — assert `map_to_contam_class(AlleleClass::CompoundAlt).is_none()`.

### `apply_observation` ([contamination_estimation.rs:721](../../../src/var_calling/contamination_estimation.rs#L721))

- **`apply_observation_does_not_overflow_at_u32_max_reads`** — `counts_by_class: [u32::MAX/2 + 1, u32::MAX/2 + 1, 0]`, `q_sum_by_class: [-1e6, -1e6, 0.0]`. Build minimal state and config, call, assert no panic. Pair with M16.
- **`apply_observation_returns_early_on_zero_total_reads`** — `counts_by_class: [0; 3]`. Assert `state.s_per_sample[s]` unchanged.

### `refresh_parameters_and_snapshot` ([contamination_estimation.rs:787](../../../src/var_calling/contamination_estimation.rs#L787))

- **`refresh_parameters_keeps_init_c_s_for_samples_with_zero_reads`** — sample 1 has `n_per_sample[1] = 0` but sample 0 has accumulated stats. Assert `state.c_s_frozen[1] == DEFAULT_C_S_INIT`. Catches: regression that collapses unobserved samples to 0.0 on the first block.
- **`refresh_parameters_q_b_uses_uniform_when_batch_has_no_reads`** — `n_per_batch[0] = 0`. Assert `q_b_frozen[0] == [DEFAULT_Q_B_INIT_PER_CLASS; 3]` (assuming M11 lands).

### `step_1a_cohort_filter` ([contamination_estimation.rs:576](../../../src/var_calling/contamination_estimation.rs#L576))

- **`step_1a_returns_not_qualifies_on_all_none_per_sample`** — `pileups.per_sample = vec![None, None, None]`. Assert `!outcome.qualifies`.
- **`step_1a_returns_not_qualifies_on_all_zero_obs`** — cohort with two allele types each at `num_obs = 0`. Catches: a refactor dropping the `total_reads == 0` check at line 605 → `0/0 = NaN`.

### `step_1b_sample_filter` ([contamination_estimation.rs:640](../../../src/var_calling/contamination_estimation.rs#L640))

- **`step_1b_returns_none_on_empty_alleles`** — `let rec = pileup_record(100, vec![]); assert!(step_1b_sample_filter(&rec, &cfg).is_none());`
- **`step_1b_picks_documented_tie_winner_on_count_tie`** — two alleles with equal `num_obs >= min_depth/2` and major_fraction ≥ threshold. Assert `g_major_class` matches the `max_by_key`-returns-last semantics. Catches: silent change of tie-breaker in a refactor.

### `finalise` ([contamination_estimation.rs:851](../../../src/var_calling/contamination_estimation.rs#L851))

- **`singleton_batch_with_floor_of_one_still_returns_none`** — `min_batch_size_for_contamination: 1`, three singleton batches. Same assertions as the existing singleton test but with floor = 1. Catches: `.max(2)` removed → singletons keep estimates. Becomes moot after M9 lands.

### `estimate_contamination` ([contamination_estimation.rs:410](../../../src/var_calling/contamination_estimation.rs#L410))

- **`site_counter_does_not_tick_when_step_1b_finds_no_samples`** — cohort polymorphic (Step 1a passes) but every sample has `depth < min_depth` (Step 1b fails for all). Feed 500 such positions in FixedSites mode with `num_sites = 1`; expect `InsufficientSites { found: 0, .. }`.
- **`estimate_contamination_propagates_first_upstream_error`** — stream yielding `Err(PerPositionMergerError::...)` on the second item. Assert `matches!(result, Err(ContaminationEstimationError::Upstream { source: _, sites_processed: 1, .. }))` (post-M7).

### `recommend_for_non_convergence` / `recommend_for_insufficient_sites` ([contamination_estimation.rs:906,944](../../../src/var_calling/contamination_estimation.rs#L906))

- **`recommend_for_non_convergence_returns_capped_suggestion`** — `(sites_processed=100, observed_delta=1.0, tolerance=1e-3)`. Raw scaling `10⁶` → cap at `10× sites` = `1000`. Assert `r.suggested_num_sites == Some(1000)`.
- **`recommend_for_non_convergence_returns_none_on_zero_delta`** — `(100, 0.0, 1e-3)`. Branch at line 911 returns `None`.
- **`recommend_for_non_convergence_returns_finite_suggestion_on_infinite_delta`** — `(100, f64::INFINITY, 1e-3)`. `ratio = INF, scaled = INF, capped = 1000.0`, `is_finite && > 0` ✓. Assert `Some(1000)`.
- **`recommend_for_insufficient_sites_suggests_observed_count`** — `(found=42, requested=10_000)`. Assert `r.suggested_num_sites == Some(42)`.

### `AlleleClass::classify` ([contamination_estimation.rs:130](../../../src/var_calling/contamination_estimation.rs#L130))

- **`classify_returns_ref_on_byte_equal_inputs`** — `b"AC" == b"AC"` → `Ref`.
- **`classify_returns_snp_alt_on_same_length_different_bytes`** — `b"AC"` vs `b"AT"` → `SnpAlt`.
- **`classify_returns_indel_alt_on_different_lengths`** — both insertion (`b"A"` vs `b"ACT"`) and deletion (`b"ACT"` vs `b"A"`).
- **`classify_returns_ref_on_both_empty`** — `b"" == b""`. Pin the edge-case contract.

### `ContaminationEstimates::zero` / `from_user_supplied` / `effective_c_s` / `q_b_for_sample`

- **`from_user_supplied_returns_bad_input_on_length_mismatch`** — paired with B1.
- **`from_user_supplied_returns_bad_input_on_c_s_out_of_range`** — `Some(-0.1)` and `Some(2.0)` → `BadInput`.
- **`from_user_supplied_returns_bad_input_on_q_b_not_simplex`** — `[0.5, 0.3, 0.3]` → `BadInput`.
- **`from_user_supplied_returns_bad_input_on_batch_idx_out_of_range`** — `sample_to_batch = vec![0, 5]` with `q_b_per_batch.len() == 2`.
- **`zero_constructor_rejects_sample_to_batch_length_mismatch_in_release`** — paired with M5 (no `#[should_panic]` needed once the runtime check lands).
- **`q_b_for_sample_returns_the_batch_distribution_for_the_sample`** — `from_user_supplied` with `sample_to_batch = vec![1, 0, 1]` and distinct `q_b_per_batch`; assert `q_b_for_sample(0) == q_b_per_batch[1]`.

### Engine wiring ([posterior_engine.rs:862](../../../src/var_calling/posterior_engine.rs#L862))

- **`engine_returns_typed_error_when_contamination_sample_to_batch_is_malformed`** (depends on M8 fix) — `ContaminationEstimates` with `c_s_per_sample.len() == n_samples` (passes current guard) but `sample_to_batch.len() != n_samples`. Currently panics with `index out of bounds`; should be a typed `PosteriorEngineError::*`.
- **`engine_uses_mixture_for_c_s_positive_sample_and_fallback_for_floored_sample_in_same_record`** — 2-sample record, `c_s_per_sample = vec![Some(0.05), None]`. Verify sample 0 takes the mixture branch and sample 1 takes the fallback copy. Catches: a misplaced `continue` at line 720 causing sample 1 to inherit garbage from `out` initialised to 0.0.

## 9. What's good

- **Type-level mutual exclusivity for the stopping mode** ([contamination_estimation.rs:148-166](../../../src/var_calling/contamination_estimation.rs#L148-L166)). `StoppingMode::{Convergence, FixedSites}` as an enum makes "both at once" unrepresentable, demoting the CLI-parse-time error from a runtime check to a type-system invariant — exactly the "make illegal states unrepresentable" pattern the codebase aims for.
- **Provenance enum on `ContaminationEstimateSource`** ([contamination_estimation.rs:236-258](../../../src/var_calling/contamination_estimation.rs#L236-L258)). Carrying `{SidePass{mode, sites_processed}, UserSupplied, Mixed, Disabled}` on every output makes the construction path inspectable from downstream code without parsing log messages.
- **Spec-citation discipline on the new architecture-level decisions.** The mixture-likelihood doc comment on `compute_mixture_log_likelihoods` and the design choices in the module-level docs of [contamination_estimation.rs:1-36](../../../src/var_calling/contamination_estimation.rs#L1-L36) explicitly reference the spec and plan and name the "site-counting interpretation B" choice — exactly the level of cross-reference that lets a reader trust the code matches the design.
- **`StoppingMode::default()` returning `Convergence` at project defaults** ([contamination_estimation.rs:168-175](../../../src/var_calling/contamination_estimation.rs#L168-L175)). Right shape: the default mode is the one with the convergence guarantee, not the one with the silent "process N sites and return whatever".

## 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --all-targets --all-features`
- `PROPTEST_CASES=512 ./scripts/dev.sh cargo test -- posteriors_sum_to_one_at_any_ploidy_and_n_alleles --nocapture` (M3 verification step — read the rejection-rate line).
- New invocations the review introduced: the tests itemised in §8 above; run individually as
  `./scripts/dev.sh cargo test -- compute_mixture_log_likelihoods_matches_hand_calculation_at_c_s_5_percent --exact` etc.

### Author response convention

Address each finding by its identifier (e.g., "B2", "M5") with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer open questions from §4 first.
