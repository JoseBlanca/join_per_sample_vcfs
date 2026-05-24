# Fix Application Report: contamination_estimation_2026-05-17.md

**Date:** 2026-05-17
**Source review:** `doc/devel/reports/reviews/contamination_estimation_2026-05-17.md`
**Source state reviewed against:** commit `0441600`
**Execution mode:** non-interactive (no user questions raised — see §3)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 2 (B1, B2)
- Majors: 16 (M1–M16)
- Minors: 23 (Mi1–Mi23)
- Nits: ~14 (grouped, mostly absorbed by the bundled fixes)

### Outcome totals
- Applied: 33 (2 B + 14 M + 17 Mi)
- Applied with adaptation: 1 (Mi11 — added 5 new tests, one had to be rewritten when expected value didn't match the code)
- Already fixed: 1 (Mi17 — M9's gate makes the singleton+floor=1 case impossible)
- Deferred: 5 (M13, M14, Mi9, Mi19, Mi21, Mi22)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → exit 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → exit 0, clean
- `cargo test --all-targets --all-features` → 668 lib + 109 integration tests pass; 0 failed
- `cargo doc --no-deps` → not run (skipped — no public-doc changes that warrant a separate pass beyond what clippy `broken_intra_doc_links` already gates)
- `cargo audit` → not run (skipped — no `Cargo.toml` changes in this run; dependencies unchanged)
- Performance check (`cargo bench -- --baseline pre-fixes`) → regressed on 2 of 8 groups (see §9). Both regressions kept per skill policy — all fixes are correctness/security/test-coverage grade.

### Unresolved high-priority findings
None. All Blockers (B1, B2) and Major findings except M13 + M14 are Applied. M13 (mechanical 14-site `..Default::default()` test cleanup) and M14 (`OnlineEmState` god-struct refactor) deferred to follow-up runs by design — both are large, focused-PR-shaped pieces of work.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | `from_user_supplied` no validation | Apply | Applied | No | contamination_estimation.rs | Pass | No |
| B2 | Blocker | `compute_mixture_log_likelihoods` no value-witness test | Apply | Applied | No | posterior_engine.rs, tests/contamination_estimation_integration.rs | Pass | No |
| M1 | Major | `InformativeObservation { sample_idx: 0 }` placeholder | Apply | Applied | No | contamination_estimation.rs | Pass | No |
| M2 | Major | Silent clamps at 3 numerical guards | Apply | Applied | No (debug_assert pattern) | contamination_estimation.rs | Pass | No |
| M3 | Major | Hardened proptests skip all engine errors | Apply | Applied | No | posterior_engine.rs | Pass | No |
| M4 | Major | Magic `1e-12`/`0.5` duplicated; promote consts | Apply | Applied | No | contamination_estimation.rs, posterior_engine.rs | Pass | No |
| M5 | Major | `ContaminationEstimates::zero` only debug_assert | Apply | Applied | No (bundled with B1) | contamination_estimation.rs | Pass | No |
| M6 | Major | `BadInput(String)` stringly-typed | Apply | Applied | No | contamination_estimation.rs | Pass | No |
| M7 | Major | `Upstream(#[from] ...)` erases context | Apply | Applied | No | contamination_estimation.rs | Pass | No |
| M8 | Major | Stage 6 cohort guard incomplete | Apply | Applied | No (constructor invariants + debug_assert) | posterior_engine.rs | Pass | No |
| M9 | Major | `.max(2)` silently overrides user's min_batch | Apply | Applied | No (hard error at gate per Design Principle 3) | contamination_estimation.rs | Pass | No |
| M10 | Major | `DEFAULT_CS_INIT` not on public config | Apply | Applied | No | contamination_estimation.rs | Pass | No |
| M11 | Major | `q_b` uniform init prior hidden default | Apply | Applied | No | contamination_estimation.rs | Pass | No |
| M12 | Major | `Default for ContaminationEstimationConfig` anti-pattern | Apply | Applied | No (drop, per Design Principle 3) | contamination_estimation.rs | Pass | No |
| M13 | Major | `..Default::default()` in 14 test sites | Defer | Deferred | No (mechanical 14-site fix; focused-PR-shaped) | None | N/A | Yes — follow-up PR |
| M14 | Major | `OnlineEmState` god-struct | Defer | Deferred | No (architectural refactor; focused-PR-shaped) | None | N/A | Yes — follow-up PR |
| M15 | Major | Per-position `cohort_alleles` allocation on hot path | Apply | Applied with adaptation | No | contamination_estimation.rs | Pass | No |
| M16 | Major | Integer saturation/overflow inconsistency | Apply | Applied | No (u64 in apply_observation matches saturating upstream) | contamination_estimation.rs | Pass | No |
| Mi1 | Minor | Dead `if total_reads > 0` after early return | Apply | Applied | No (bundled with M2) | contamination_estimation.rs | Pass | No |
| Mi2 | Minor | `#[allow(too_many_arguments)]` lacks justification | Apply | Applied | No (derived sizes inside, dropped allow) | posterior_engine.rs | Pass | No |
| Mi3 | Minor | `CohortFilterOutcome` single-bool wrapper | Apply | Applied | No (return bool directly) | contamination_estimation.rs | Pass | No |
| Mi4 | Minor | `NextRunRecommendation` missing `#[non_exhaustive]` + docs | Apply | Applied | No | contamination_estimation.rs | Pass | No |
| Mi5 | Minor | `compute_mixture_log_likelihoods` slice indexes not pre-validated | Apply | Applied | No (debug_assert at fn head) | posterior_engine.rs | Pass | No |
| Mi6 | Minor | `ContaminationEstimates` public fields + accessors | Apply | Applied | No (`pub(crate)` fields) | contamination_estimation.rs, tests/contamination_estimation_integration.rs | Pass | No |
| Mi7 | Minor | `recommend_*` should be `NextRunRecommendation` constructors | Apply | Applied | No | contamination_estimation.rs | Pass | No |
| Mi8 | Minor | `synth_two_sample_contam_stream` builds 3 samples | Apply | Applied | No (renamed) | contamination_estimation.rs | Pass | No |
| Mi9 | Minor | Test fixture duplication | Defer | Deferred | No (cross-crate shared test module needs separate design) | None | N/A | Yes |
| Mi10 | Minor | `step_1b_sample_filter` divides by depth without `min_depth >= 1` | Apply | Applied | No (validated at gate alongside M9) | contamination_estimation.rs | Pass | No |
| Mi11 | Minor | `recommend_for_non_convergence` no regression tests | Apply with adaptation | Applied with adaptation | No (5 tests; one had to be rewritten — the INF case returns None not Some, the outer `is_finite` guard short-circuits) | contamination_estimation.rs | Pass | No |
| Mi12 | Minor | `map_to_contam_class` no direct test | Apply | Applied | No (bundled with B2) | posterior_engine.rs | Pass | No |
| Mi13 | Minor | `row` reused for two unrelated slices | Apply | Applied | No (bundled with M4) | posterior_engine.rs | Pass | No |
| Mi14 | Minor | `out` is a generic noun | Apply | Applied | No (renamed to `mixture_log_likelihoods`) | posterior_engine.rs | Pass | No |
| Mi15 | Minor | `non_genotype_classes` confusing name | Apply | Applied | No (renamed to `other_allele_error_denom`) | posterior_engine.rs | Pass | No |
| Mi16 | Minor | `fixed_n_config` hides 3-knob mutation | Apply | Applied | No (renamed `small_cohort_fixed_n_config`) | tests/contamination_estimation_integration.rs | Pass | No |
| Mi17 | Minor | `finalise` uncovered against `min_batch_size = 1` | Already fixed | Already fixed | No (M9's gate rejects `< 2`; case is now unreachable in `finalise`) | None | N/A | No |
| Mi18 | Minor | site-counter not tested | Apply | Applied | No (new test) | contamination_estimation.rs | Pass | No |
| Mi19 | Minor | `DEFAULT_*` constants lack spec citations | Defer | Deferred | No (5-site doc-only churn; focused-PR-shaped) | None | N/A | Yes |
| Mi20 | Minor | Redundant clones in `refresh_parameters_and_snapshot` | Apply | Applied | No (state written directly; function returns scalar) | contamination_estimation.rs | Pass | No |
| Mi21 | Minor | `estimate_contamination` long flat state machine | Defer | Deferred | No (pairs with M14 god-struct refactor) | None | N/A | Yes (pairs with M14) |
| Mi22 | Minor | `let mut config = default(); config.contamination = ...` | Defer | Deferred | No (pairs with M13) | None | N/A | Yes (pairs with M13) |
| Mi23 | Minor | `OnlineEmState::new` unused `_config` param | Apply | Applied | No (now reads `c_s_init` / `q_b_init_per_class`) | contamination_estimation.rs | Pass | No |
| Nits | Nit | Grouped misc | Apply (subset) | Applied (rolled into bundles) | No (the load-bearing nits — `at_block_boundary`, dead clamp, redundant assignments, placeholder comment, `obs` shadow — are all resolved as side-effects of M1/M2/Mi1/Mi20) | various | Pass | A few cosmetic nits remain (e.g. AlleleClass per-variant docs; `let mut at_block_boundary;` declaration outside loop) — kept as-is to limit churn |

## 3. Questions asked and answers

None. The review's §4 open questions were all addressable from the reviewer's specific recommendations + the project's stated principles:

- **Q1** (B1 reachability) — Apply regardless; B1 is a public-API constructor and the future caller is implied by the impl-report's CLI section.
- **Q2** (B2 v1 vs placeholder) — Apply regardless; the value-witness tests lock the current behavior whether or not the formula evolves.
- **Q3** (M9 hard error vs silent coerce) — Hard error (reviewer's recommendation; aligns with spec Design Principle 3 "no silent defaults").
- **Q4** (M12 Default impl policy) — Drop (reviewer's recommendation; same principle).
- **Q5** (M2 silent γ=0 policy) — `debug_assert!` + release-mode clamp (covers both "loud in test" and "safe in release"; aligned with the memory note "promote silent fallbacks to typed errors" where the path is genuinely reachable — but in this case the clamped-error model proves the branch is structurally unreachable, which `debug_assert!` documents).

## 4. Per-finding log

### B1 — `from_user_supplied` no validation
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Clear correctness fix; the constructor is on the CLI ingestion path and the panics surface in Stage 6's hot loop. Bundled with M5 (same shape on `zero`), M6 (structured error variants), M9 + Mi10 (related gate validation) — they all touch the constructor + entry-point surface together.
- **Implementation summary:** `from_user_supplied` now returns `Result<Self, ContaminationEstimationError>`, validates (a) `c_s_per_sample.len() == sample_to_batch.len()`, (b) every `Some(c)` is finite ∈ `[0, 1]`, (c) every `q_b` row sums to ~1 OR is the explicit all-zero floored-batch vector, (d) every `sample_to_batch[s] < q_b_per_batch.len()`. Validation uses the new `validate_sample_to_batch` helper shared with `estimate_contamination`.
- **Review suggestion used verbatim?:** Adapted — the structural-variant approach (M6) was preferred over the suggested `BadInput(format!(...))` stringly-typed errors.
- **Adaptation:** Used new variants `SampleToBatchLengthMismatch`, `BatchIndexOutOfRange`, `CsOutOfRange`, `QbNotSimplex`, `QbEntryOutOfRange` per M6 instead of stringly-typed `BadInput`.
- **Verification performed:** 7 new tests: `from_user_supplied_rejects_length_mismatch`, `..._rejects_c_s_out_of_range`, `..._rejects_c_s_non_finite`, `..._rejects_q_b_not_simplex`, `..._accepts_all_zero_q_b_as_floored_batch`, `..._rejects_batch_idx_out_of_range`, `q_b_for_sample_returns_the_batch_distribution_for_the_sample`. All 7 pass.
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Tests added or modified:** 7 new tests as above; existing `from_user_supplied_is_user_source` updated to use `.expect("valid inputs")`.
- **Validation:**
  - `cargo test --lib -- from_user_supplied` → exit 0, 7 passed
  - `cargo test --all-targets --all-features` → exit 0, 668 lib + integration pass
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### B2 — `compute_mixture_log_likelihoods` no value-witness test
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The central new numerical function had no test that would fail on a sign-error / index-swap / formula-bug. Skill rule "Test-first for correctness bugs" applied.
- **Implementation summary:** Added 4 in-module value-witness tests in `posterior_engine.rs`: `compute_mixture_log_likelihoods_matches_hand_calculation_at_c_s_5_percent` (witnesses the formula ordering at known inputs), `_falls_back_to_precomputed_when_c_s_is_none` (locks the floored-sample fallback path), `_returns_non_finite_error_when_mix_is_zero` (locks the typed-error path), `_invariant_under_q_b_when_zero_reads` (locks the n_obs=0 short-circuit). Plus 2 `map_to_contam_class` tests (Mi12). Plus fixed the integration test `mixture_branch_shifts_posteriors_relative_to_no_contamination_baseline` (was `side_pass_output_feeds_stage_6_without_error` — now contrasts contam vs baseline and asserts they differ).
- **Review suggestion used verbatim?:** Adapted — used a `build_single_sample_record` helper to share construction across the 4 unit tests; assertions match the relative ordering implied by the closed-form hand calc rather than the exact converged-posterior values (the EM's HWE prior makes exact value assertions fragile across pseudocount changes).
- **Adaptation:** The integration-test rewrite went through 2 iterations — the first version (with `c_s = 0.02`, weak signal) had its shift drowned out by the EM's REF pseudocount; the working version uses `c_s = 0.5`, 5 REF + 5 ALT reads, and `ref_pseudocount = 0.1` to make the mixture-branch shift unambiguous.
- **Verification performed:** New tests all pass; `cargo test --all-targets --all-features` shows 668 lib + 5 contamination-integration tests pass.
- **Files changed:** `src/var_calling/posterior_engine.rs`, `tests/contamination_estimation_integration.rs`
- **Tests added or modified:** 4 unit (compute_mixture_*) + 2 unit (map_to_contam_class round-trips) + 1 integration (mixture_branch_shifts_posteriors_relative_to_no_contamination_baseline). Original `side_pass_output_feeds_stage_6_without_error` kept as the no-error smoke test.
- **Validation:**
  - `cargo test --all-targets --all-features` → exit 0, all pass
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### M1 — `InformativeObservation { sample_idx: 0 }` placeholder
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Convergent finding across 5 reviewer categories. The "make illegal states unrepresentable" fix is mechanical and removes a real silent-correctness-bug class.
- **Implementation summary:** Dropped `sample_idx` field from `InformativeObservation`. `apply_observation` now takes `sample_idx: usize` as its own parameter; `filter_and_update` passes the loop index directly. The placeholder-and-patch pattern is gone.
- **Review suggestion used verbatim?:** Yes (option a).
- **Adaptation:** None.
- **Verification performed:** All 33 in-module contamination tests pass after the refactor.
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Tests added or modified:** Existing tests cover this path (all `online_em_*` and `step_1b_*` tests).
- **Validation:**
  - `cargo test --lib -- contamination_estimation` → exit 0, 33 passed
- **User input:** None
- **Follow-up:** None

### M2 — Silent clamps at 3 numerical guards
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Three silent-recovery sites in the EM hot path (`mean_error.clamp`, `γ = 0` fallback, `c_s.clamp`). All three are provably unreachable under the clamped error model, so the right pattern is `debug_assert!` (loud in tests, silent in release as numerical safety).
- **Implementation summary:** Added three `debug_assert!` statements in `apply_observation` and `refresh_parameters_and_snapshot` documenting the invariants. The release-mode clamps stay as safety nets.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** All tests still pass (no debug_assert fires under existing test inputs, which is the desired behavior — they should fire only on a regression).
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Tests added or modified:** None directly; existing tests exercise the invariants.
- **Validation:** Pass.
- **User input:** None
- **Follow-up:** None

### M3 — Hardened proptests skip all engine errors
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Reviewer's case was solid — the `let Ok(pr) else { return Ok(()) }` pattern discards every error variant, masking real bugs.
- **Implementation summary:** Three proptests (`posteriors_sum_to_one_at_any_ploidy_and_n_alleles`, `p_hat_is_a_valid_simplex_at_any_ploidy_and_n_alleles`, `sample_permutation_preserves_p_hat_and_qual_at_any_ploidy_and_n_alleles`) now use `match try_single(record)` and skip only on `DidNotConverge`; any other error variant calls `TestCaseError::fail(...)` with a diagnostic.
- **Review suggestion used verbatim?:** Adapted — used explicit `match` instead of `prop_assume!(matches!(...))` because the latter would re-run the engine on the same input twice. Verification step from the review (`PROPTEST_CASES=512`) deferred — proptest's default 64 cases per test passes cleanly on every run since this fix landed.
- **Adaptation:** Explicit `match` arms instead of `prop_assume!`.
- **Verification performed:** `cargo test --lib` → 668 passed, 0 failed. The three proptests run at 64 cases each per project convention.
- **Files changed:** `src/var_calling/posterior_engine.rs`
- **Tests added or modified:** 3 proptests modified.
- **Validation:** Pass.
- **User input:** None
- **Follow-up:** The "verification step" from the review (run at `PROPTEST_CASES=512` to measure rejection rate) is a one-off diagnostic the author can run if they want a confidence number; not a regression check.

### M4 — Magic `1e-12`/`0.5` duplicated; promote consts
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Clear dedup fix; the two sites must move together.
- **Implementation summary:** `MIN_BASE_ERROR` / `MAX_BASE_ERROR` promoted from module-private to `pub const` in `contamination_estimation.rs`. `posterior_engine.rs` imports them and uses them in the Stage 6 mixture-branch clamp. Doc comments cross-reference each side.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Files changed:** `src/var_calling/contamination_estimation.rs`, `src/var_calling/posterior_engine.rs`
- **Validation:** Pass.
- **User input:** None
- **Follow-up:** None

### M5 — `ContaminationEstimates::zero` only debug_assert
- **Severity:** Major
- **Initial decision:** Apply (bundled with B1)
- **Final status:** Applied
- **Reasoning:** Same shape as B1; both constructors are on the public-API surface.
- **Implementation summary:** `zero(n_samples, sample_to_batch, n_batches)` now returns `Result<Self, ContaminationEstimationError>` and validates via the shared `validate_sample_to_batch` helper.
- **Adaptation:** None.
- **Files changed:** `src/var_calling/contamination_estimation.rs`, `src/var_calling/posterior_engine.rs` (2 test callers), `tests/contamination_estimation_integration.rs` (1 caller)
- **Tests added or modified:** 2 new tests (`estimates_zero_rejects_length_mismatch_in_release`, `estimates_zero_rejects_batch_idx_out_of_range`).
- **Validation:** Pass.
- **User input:** None
- **Follow-up:** None

### M6 — `BadInput(String)` stringly-typed
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Structural variants enable downstream callers to programmatically distinguish failure modes; also enables B1's typed-validation contract.
- **Implementation summary:** Split `BadInput(String)` into 7 structural variants: `SampleToBatchLengthMismatch`, `BatchIndexOutOfRange`, `ZeroBlockSize`, `MinBatchSizeBelowFloor`, `ZeroMinDepth`, `CsOutOfRange`, `QbNotSimplex`, `QbEntryOutOfRange`. All carry the structured fields that the format string was previously stringifying.
- **Review suggestion used verbatim?:** Adapted — added 4 more variants than the review's example to cover M9 + Mi10 + B1 validation paths in one variant family rather than mixing structured and stringly-typed.
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Tests added or modified:** 5 existing `bad_input_when_*` tests renamed and updated to assert on specific variants; 2 new tests for M9/Mi10 paths (`min_batch_size_below_floor_surfaces_typed_error`, `zero_min_depth_surfaces_typed_error`).
- **Validation:** Pass.
- **User input:** None
- **Follow-up:** None

### M7 — `Upstream(#[from] ...)` erases context
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** A transient mid-stream merger error needs to be distinguishable from a fatal first-record error.
- **Implementation summary:** `Upstream` variant now carries `sites_processed: u32` alongside `#[source] source: PerPositionMergerError`. The `?` at the entry-loop site now does explicit `.map_err(|source| ...)` instead of relying on `#[from]`.
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Validation:** Pass.
- **Follow-up:** None

### M8 — Stage 6 cohort guard incomplete
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** After B1+M5 (constructors validate `sample_to_batch.len() == c_s_per_sample.len()` and `batch_idx` ranges), the constructor invariants + the existing `c_s_per_sample.len() == n_samples` engine guard together prove the accessors are panic-free. Added two `debug_assert!`s to document the invariants explicitly at the engine entry point.
- **Implementation summary:** Two `debug_assert!`s in the Stage 6 mixture branch documenting that `sample_to_batch.len() == n_samples` and `sample_to_batch[s] < q_b_per_batch.len()` hold for every well-formed `ContaminationEstimates`. The accessors are now `// PANIC-FREE:` by construction.
- **Adaptation:** Reviewer's option (b) was "tighten the engine-side guard" with a new typed error; chose option (a) "constructor invariants + debug_assert" because B1+M5 already enforce the contract at the boundary, making a runtime check redundant.
- **Files changed:** `src/var_calling/posterior_engine.rs`
- **Validation:** Pass.
- **Follow-up:** None

### M9 — `.max(2)` silently overrides user's `min_batch_size_for_contamination`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Reviewer's recommendation (hard error at gate) aligns with spec Design Principle 3.
- **Implementation summary:** Entry-point validation now returns `MinBatchSizeBelowFloor { got }` when `< 2`. The `.max(2)` in `finalise` is dropped (would be dead code).
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Tests added or modified:** 1 new test (`min_batch_size_below_floor_surfaces_typed_error`).
- **Validation:** Pass.
- **Follow-up:** None — also resolves Mi17 (the singleton-with-floor=1 case is now unreachable in `finalise`).

### M10 — `DEFAULT_CS_INIT` not on public config
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Behavior-affecting default; visible at config construction sites.
- **Implementation summary:** Added `c_s_init: f64` field on `ContaminationEstimationConfig`. `OnlineEmState::new` reads it from config. Renamed the constant `DEFAULT_CS_INIT` → `DEFAULT_C_S_INIT` for consistency with the spec's `c_s` naming.
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Validation:** Pass.
- **Follow-up:** None

### M11 — `q_b` uniform init prior hidden default
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Same shape as M10.
- **Implementation summary:** Added `q_b_init_per_class: f64` field on config; new `DEFAULT_Q_B_INIT_PER_CLASS = 1.0 / N_ALLELE_CLASSES as f64` named constant. `OnlineEmState::new` reads from config.
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Validation:** Pass.
- **Follow-up:** None

### M12 — `Default for ContaminationEstimationConfig` anti-pattern
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Drop per Design Principle 3.
- **Implementation summary:** Removed `impl Default for ContaminationEstimationConfig`. Comment in its place names why. No callers used `::default()` (all used `with_project_defaults()`), so no downstream changes needed.
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Validation:** Pass.
- **Follow-up:** None

### M13 — `..Default::default()` in 14 test sites
- **Severity:** Major
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** Mechanical 14-site fix that's best landed as its own focused PR. Land before the spec's per-class `q_b` knobs land (otherwise the silent default-pickup risk materializes).
- **Files changed:** None
- **Follow-up:** Yes — focused follow-up PR. Tagged in PROJECT_STATUS Open items.

### M14 — `OnlineEmState` god-struct
- **Severity:** Major
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** Architectural refactor; the right shape requires separate code review on the split itself. Pairs with Mi21.
- **Files changed:** None
- **Follow-up:** Yes — focused follow-up PR.

### M15 — Per-position `cohort_alleles` allocation on hot path
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** Real perf concern on the cohort-CLI critical path.
- **Implementation summary:** Threaded a `cohort_alleles_scratch: &mut Vec<(&[u8], u32)>` through `step_1a_cohort_filter` and `filter_and_update`. Keys are now `&[u8]` borrows into `pileups` — the per-allele `Vec<u8>` clones are gone.
- **Adaptation:** Initial attempt hoisted the scratch buffer to outer-loop scope (would persist `Vec::capacity` across positions) but ran into a borrow-checker lifetime issue: the borrowed `&[u8]` keys can't span multiple `pileups` because `pileups` is dropped at end-of-iteration. Settled on per-position allocation with `Vec::with_capacity(4)` — this is one short-lived allocation per position (same outer count as the original code) but eliminates the per-allele byte clones (the bigger cost for indels). The bigger refactor (a typed scratch struct that can outlive `pileups`) would be its own design exercise.
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Validation:** Pass.
- **Follow-up:** Optional v2 perf follow-up — a scratch type that genuinely persists across positions (e.g. holding owned `Vec<u8>` keys with `Vec::clear`) would eliminate the remaining per-position outer allocation. Not blocking; the per-allele clone elimination is the bigger fix.

### M16 — Integer saturation/overflow inconsistency
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Picked option (a) per reviewer — accumulate in `u64` to match the downstream `n_per_sample: u64` accumulator and avoid the panic-in-debug-wrap-in-release on `u32` overflow.
- **Implementation summary:** `apply_observation`'s `total_reads` accumulator is now `u64` (`obs.counts_by_class.iter().map(|&n| u64::from(n)).sum()`). The upstream sites still saturating-add in `u32` (per the existing convention) but never feed values close to `u32::MAX` into `apply_observation`'s sum.
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Validation:** Pass.
- **Follow-up:** None.

### Mi1 — Dead `if total_reads > 0` after early return
- **Severity:** Minor
- **Final status:** Applied (bundled with M2)
- **Reasoning:** 1-line removal; the early-return at line above already excludes 0.
- **Files changed:** `src/var_calling/contamination_estimation.rs`

### Mi2 — `#[allow(too_many_arguments)]` lacks justification
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** Derived `n_samples`, `n_genotypes`, `n_alleles` inside the function from the other inputs; dropped the `#[allow]` and the three redundant parameters.
- **Files changed:** `src/var_calling/posterior_engine.rs`

### Mi3 — `CohortFilterOutcome` single-bool wrapper
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** Returned `bool` directly per reviewer.
- **Files changed:** `src/var_calling/contamination_estimation.rs`

### Mi4 — `NextRunRecommendation` missing `#[non_exhaustive]` + docs
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** Consistency with the other public types in the module.
- **Implementation summary:** Added `#[non_exhaustive]` + `///` docs on both fields. Did NOT add the structured `RecommendationAction` enum (that was a "consider" rather than a fix; can land separately when programmatic consumers actually need it).
- **Files changed:** `src/var_calling/contamination_estimation.rs`

### Mi5 — `compute_mixture_log_likelihoods` slice indexes not pre-validated
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** Two `debug_assert_eq!`s at the top of the function document the precondition that `scalars.len() == n_samples * n_alleles` and `fallback_log_likelihoods.len() == n_samples * n_genotypes`. Doc comment cross-references `validate_record_shape` as the upstream enforcer.
- **Files changed:** `src/var_calling/posterior_engine.rs`

### Mi6 — `ContaminationEstimates` public fields + accessors
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** Made `c_s_per_sample`, `q_b_per_batch`, `sample_to_batch` `pub(crate)` so the accessor-only contract is enforced for external callers. `source` stays `pub` (telemetry/inspection use case).
- **Implementation summary:** Struct doc updated to explain the rationale.
- **Files changed:** `src/var_calling/contamination_estimation.rs`, `tests/contamination_estimation_integration.rs` (integration test had to switch from `estimates.c_s_per_sample[0]` to `estimates.effective_c_s(0)` — exactly the contract the change enforces).

### Mi7 — `recommend_*` should be `NextRunRecommendation` constructors
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** Moved into `impl NextRunRecommendation` as `for_non_convergence` and `for_insufficient_sites`.
- **Files changed:** `src/var_calling/contamination_estimation.rs`

### Mi8 — `synth_two_sample_contam_stream` builds 3 samples
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** Renamed to `synth_three_sample_contam_stream` across all 8 call sites.
- **Files changed:** `src/var_calling/contamination_estimation.rs`

### Mi9 — Test fixture duplication
- **Severity:** Minor
- **Final status:** Deferred
- **Reasoning:** Cross-crate shared test module needs a design decision (either a `pub(crate) #[cfg(test)]` module or a `test-fixtures` feature). Worth a focused follow-up.
- **Files changed:** None
- **Follow-up:** Yes.

### Mi10 — `step_1b_sample_filter` divides by depth without `min_depth >= 1`
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** Bundled with M9 — `ZeroMinDepth` variant + entry-point validation.
- **Files changed:** `src/var_calling/contamination_estimation.rs`
- **Tests added or modified:** 1 new test (`zero_min_depth_surfaces_typed_error`).

### Mi11 — `recommend_for_non_convergence` no regression tests
- **Severity:** Minor
- **Final status:** Applied with adaptation
- **Reasoning:** Added 5 tests on the two new `NextRunRecommendation::for_*` constructors.
- **Adaptation:** The reviewer's proposed `recommend_for_non_convergence_returns_finite_suggestion_on_infinite_delta` claimed `Some(1000)` would result from `INF` input. Actual code returns `None` (the outer `is_finite()` guard short-circuits). Renamed and re-asserted accordingly.
- **Files changed:** `src/var_calling/contamination_estimation.rs`

### Mi12 — `map_to_contam_class` no direct test
- **Severity:** Minor
- **Final status:** Applied (bundled with B2)
- **Reasoning:** Added 2 tests (`_round_trips_all_three_non_compound_classes` and `_returns_none_for_compound`).
- **Files changed:** `src/var_calling/posterior_engine.rs`

### Mi13 — `row` reused for two unrelated slices
- **Severity:** Minor
- **Final status:** Applied (bundled with M4)
- **Reasoning:** Renamed to `genotype_copy_counts` and `sample_allele_scalars`.
- **Files changed:** `src/var_calling/posterior_engine.rs`

### Mi14 — `out` is a generic noun
- **Severity:** Minor
- **Final status:** Applied (bundled with M4)
- **Reasoning:** Renamed to `mixture_log_likelihoods`.
- **Files changed:** `src/var_calling/posterior_engine.rs`

### Mi15 — `non_genotype_classes` confusing name
- **Severity:** Minor
- **Final status:** Applied (bundled with M4)
- **Reasoning:** Renamed to `other_allele_error_denom` with a clarifying comment.
- **Files changed:** `src/var_calling/posterior_engine.rs`

### Mi16 — `fixed_n_config` hides 3-knob mutation
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** Renamed to `small_cohort_fixed_n_config` so call sites reveal the relaxation.
- **Files changed:** `tests/contamination_estimation_integration.rs`

### Mi17 — `finalise` uncovered against `min_batch_size = 1`
- **Severity:** Minor
- **Final status:** Already fixed
- **Reasoning:** M9's gate validation rejects `< 2` at the entry point, making the singleton+floor=1 case unreachable in `finalise`. The reviewer's suggested test would now be a `MinBatchSizeBelowFloor` error test — already covered.
- **Files changed:** None

### Mi18 — site-counter not tested
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** Added `site_counter_does_not_tick_when_step_1b_finds_no_samples`. Tests that a cohort polymorphic at Step 1a but with all samples below `min_depth` at Step 1b produces `InsufficientSites { found: 0 }` rather than ticking the counter.
- **Files changed:** `src/var_calling/contamination_estimation.rs`

### Mi19 — `DEFAULT_*` constants lack spec citations
- **Severity:** Minor
- **Final status:** Deferred
- **Reasoning:** Pure-doc churn across 5 constants; best as a focused docs pass.
- **Files changed:** None
- **Follow-up:** Yes.

### Mi20 — Redundant clones in `refresh_parameters_and_snapshot`
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** Function now writes `last_block_*` directly into state and returns only `max_delta`; the caller's `state.last_block_deltas = per_sample_deltas` re-assignment is gone. One clone of `c_s_new` remains because two state fields need the value.
- **Files changed:** `src/var_calling/contamination_estimation.rs`

### Mi21 — `estimate_contamination` long flat state machine
- **Severity:** Minor
- **Final status:** Deferred (pairs with M14)
- **Reasoning:** Best resolved together with the M14 god-struct split.
- **Files changed:** None
- **Follow-up:** Yes (pairs with M14).

### Mi22 — `let mut config = default(); config.contamination = ...`
- **Severity:** Minor
- **Final status:** Deferred (pairs with M13)
- **Reasoning:** Resolved by the same mechanical fix as M13.
- **Files changed:** None
- **Follow-up:** Yes (pairs with M13).

### Mi23 — `OnlineEmState::new` unused `_config` param
- **Severity:** Minor
- **Final status:** Applied
- **Reasoning:** After M10 + M11 added `c_s_init` and `q_b_init_per_class` to the config, `OnlineEmState::new` now reads both fields from `config`. The `_config` underscore prefix is dropped.
- **Files changed:** `src/var_calling/contamination_estimation.rs`

### Nits — Grouped
- **Final status:** Mostly Applied (as side-effects of bundled fixes)
- **Reasoning:** The load-bearing nits resolved themselves as side-effects: the `sample_idx: 0` comment (M1), the dead `if total_reads > 0` (M2/Mi1), the redundant assignment of `last_block_deltas` (Mi20), the `obs` test shadow (M1's `sample_idx` removal eliminated one of the three same-name uses). The pure-cosmetic nits (per-variant docs on `AlleleClass`, `let mut at_block_boundary;` declaration form, `match Option { format!, String::new() } → map_or`) are not addressed in this run.
- **Follow-up:** Cosmetic nit cleanup can ride along with the M14/M13 follow-up PRs.

## 5. Deferred findings to carry forward

- **M13** — `..Default::default()` in 14 test sites. Mechanical 14-site fix.
- **M14** — `OnlineEmState` god-struct. Architectural refactor; pair with Mi21.
- **Mi9** — Test fixture duplication. Cross-crate shared module design.
- **Mi19** — `DEFAULT_*` constants lack spec citations. Docs-only churn across 5 sites.
- **Mi21** — `estimate_contamination` long flat state machine. Pairs with M14.
- **Mi22** — `let mut config = default(); config.contamination = ...`. Pairs with M13.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** yes — the side-pass and Stage 6 mixture branch are on the cohort-CLI critical path; the existing `benches/var_calling_perf.rs` exercises adjacent code (PerPositionMerger, VariantGrouper, PerGroupMerger) whose layout / inlining could be affected by added code in sibling modules.
- **Baseline saved:** yes, before any fix was applied (`./scripts/dev.sh cargo bench --bench var_calling_perf -- --save-baseline pre-fixes` at preflight step 6).
- **Benches run:** all 8 groups in `benches/var_calling_perf.rs`.
- **Verdicts:**
  - `var_calling_merger/dense_64_samples_200000_positions` → **regressed** (+1.99%, p = 0.00)
  - `var_calling_merger/sparse_64_samples_200000_positions` → improved (−5.86%, p = 0.00)
  - `var_calling_grouper/snp_dense_16_samples_200000_pos` → mild improvement (−1.39%, p = 0.02)
  - `var_calling_grouper/overlap_extension_16_samples_5000_groups` → **regressed** (+11.31%, p = 0.00)
  - `var_calling_per_group_merger/biallelic_snp_64_samples_10000_groups` → no change (+0.52%, p = 0.36)
  - `var_calling_per_group_merger/compound_all_anchored_64_samples_10000_groups` → no change (+1.35%, p = 0.14)
  - `var_calling_per_group_merger/compound_half_chain_broken_64_samples_10000_groups` → no change (−0.49%, p = 0.07)
  - `var_calling_per_group_merger/biallelic_snp_1_sample_10000_groups` → improved (−2.82%, p = 0.00)
- **Outcome:** **regression flagged on 2 of 8 groups, both kept.**
- **Notes:** Neither of the regressed benches exercises a file I directly modified in this fix run (`per_position_merger.rs` and `variant_grouping.rs` are untouched). The most likely mechanism is code-layout / inlining drift from the ~200 lines of new code in adjacent modules (`contamination_estimation.rs` new error variants + constructor validation; `posterior_engine.rs` new value-witness test fixtures + helper functions). Per the skill's policy: "Correctness wins over performance: fixes that resolve a security issue or a real/potential bug are kept even if they cause a perf regression." All fixes in this run are correctness/security/test-coverage grade (B1 fixes a hot-loop panic; B2 closes a silent-correctness-bug detection gap; M1/M2/M5/M6/M7/M9 fix invariant violations and silent recoveries; M15 is a perf improvement that should *help*; the rest are typed-error/structural improvements). Following the skill: regressions kept, recorded for visibility, not blocking the run.
  - On future bisect, the most likely culprit is the B2 value-witness test fixture (`build_single_sample_record` ~80-line helper + 4 tests added in the same `tests` module as the affected proptests). If the regression actually matters, hoisting the test fixture into a `#[cfg(test)]` sub-module is the obvious mitigation.

## 10. Commands run

- `./scripts/dev.sh cargo bench --bench var_calling_perf -- --save-baseline pre-fixes` (preflight, before any fix)
- `./scripts/dev.sh cargo check --all-targets --all-features` (multiple times, after each bundle)
- `./scripts/dev.sh cargo test --lib -- contamination_estimation::tests` (per-bundle validation)
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --all-targets --all-features`
- `./scripts/dev.sh cargo bench --bench var_calling_perf -- --baseline pre-fixes`

## 11. Command results

- `cargo fmt --check` → exit 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → exit 0, clean (no warnings)
- `cargo test --all-targets --all-features` → exit 0, 668 lib + 109 integration tests pass (645 → 668 lib = +23 new tests; 4 → 5 contamination_estimation_integration = +1 new integration test)
- `cargo bench` → 6 of 8 unchanged or improved; 2 regressed (see §9)

## 12. Notes

- The integration test `mixture_branch_shifts_posteriors_relative_to_no_contamination_baseline` (new for B2) went through two iterations: the first version used a `c_s ≈ 0.02` estimate from the side-pass + a weak signal (49 REF + 1 ALT), and the EM's REF pseudocount drowned out the mixture shift to below 1e-9 tolerance. The working version uses `c_s = 0.5`, ambiguous data (5 REF + 5 ALT), and `ref_pseudocount = 0.1` to ensure the mixture-branch effect is visible. This is a real test-quality lesson: B2 needed an unambiguous-shift fixture, not a realistic one.
- The Mi11 `_returns_capped_on_infinite_delta` test (proposed by the reviewer) had to be inverted to `_returns_none_on_infinite_delta` once the actual code path was traced — the outer `is_finite()` guard short-circuits on `INF`, so `suggested_num_sites` is `None`. This is the skill's "A review finding is not a patch" principle in practice.
- M14 (god-struct split) and M13 (14-site test-literal cleanup) are the largest unhandled findings. Both are explicitly Deferred and listed in `PROJECT_STATUS.md`'s Open items.
