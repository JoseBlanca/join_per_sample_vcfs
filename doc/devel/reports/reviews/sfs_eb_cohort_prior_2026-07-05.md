# Code Review: SFS empirical-Bayes leave-one-out cohort prior (Milestone 4)

**Date:** 2026-07-05
**Reviewer:** rust-code-review skill (orchestrator + 4 parallel category sub-agents)
**Scope:** working-tree diff on branch `sfs-genotype-prior` — the EM engine change adding
empirical-Bayes large-cohort sharpening (`src/var_calling/posterior_engine.rs`,
`tests/posterior_engine_integration.rs`).
**Status:** Approve-with-changes → **all findings applied** in the same session.

---

## 1. Scope

- **Reviewed:** the Milestone-4 diff — the leave-one-out (LOO) cohort genotype prior
  (arch [sfs_genotype_prior.md §10](../../architecture/sfs_genotype_prior.md)).
- **In-scope files:** [posterior_engine.rs](../../../../src/var_calling/posterior_engine.rs),
  [posterior_engine_integration.rs](../../../../tests/posterior_engine_integration.rs).
- **Out of scope:** the design note (docs); untouched engine code; the pre-existing
  `vcf/writer.rs:384` `needless_lifetimes` lint.
- **Categories dispatched (4 parallel sub-agents):** reliability; refactor_safety;
  idiomatic + smells; errors + naming. Audit trail in
  `tmp/review_2026-07-05_sfs-eb-cohort-prior/`.

## 2. Verdict

**Approve-with-changes.** No Blockers. No Majors on the production logic; the two Majors
were both missing-test coverage on silently-wrong-result paths. All findings applied.

## 3. Execution status

- `cargo build --lib` — OK.
- `cargo test --lib` — 1566 passed, 0 failed, 2 ignored (pre-existing).
- `cargo test --test posterior_engine_integration` — 11 passed, 0 failed.
- `cargo clippy --all-targets` — clean for the changed code (only the pre-existing
  `vcf/writer.rs:384` lint remains).
- `rustfmt --check` — clean on both changed files.

## 4. Findings and resolutions

### Major (both missing-coverage on silently-wrong paths)

- **M1 — LOO subtraction depends on an untested cross-path count invariant.**
  `expected_counts` (cohort total) uses `accumulate_expected_counts`'s hand-unrolled
  biallelic-diploid fast path, while `own_counts` uses the general `accumulate_row_counts`;
  the LOO counts `total − own` are only correct if the two agree per allele. **Applied:**
  new unit test `accumulate_expected_counts_fast_path_matches_general_path_for_biallelic_diploid`
  pins the two paths equal for a `(2,2)` multi-sample posterior matrix.

- **M2 — the multi-sample non-finite guards were untested** (the only non-finite test is
  single-sample and hits the species path). **Applied:** new test
  `cohort_flat_step_surfaces_non_finite_posterior_for_bad_non_first_sample` drives a cohort
  with an all-`-∞` row on sample index 1 and asserts `NonFinitePosterior { sample_idx: 1 }`
  — exercising `e_step_flat`'s own guard and its per-sample index.

### Minor

- **Mi1 — `.max(0.0)` LOO clamp silently masks a non-noise desync** (errors). **Applied:**
  added a `debug_assert!(raw_loo > -1e-6, …)` before the clamp so a materially-negative
  `total − own` (a future count-path desync) fails loudly in debug while release keeps the
  clamp.
- **Mi2 — `fill_log_indep_per_g_from` had an unchecked three-slice length contract**
  (refactor_safety). **Applied:** two `debug_assert_eq!` pinning
  `out.len() == shape.nonzero_pairs_offsets.len() == shape.log_multinomial_coeffs.len()`.
- **Mi3 — the log-normalise-and-write-back tail was duplicated in four E-step bodies**
  (smells). **Applied:** extracted `normalise_logits_into_posteriors`; `e_step`, the
  `e_step_simd` scalar tail, `e_step_flat`, and `e_step_cohort_loo` now share it. The
  single-sample byte-identity is preserved (same ops, same order; full suite green).
- **Mi4 — cohort could report `converged: true` on the flat first iteration**, so the LOO
  prior producing the emitted genotypes never enters the convergence test (reliability).
  **Applied:** `run_em_loop` now refuses to converge on iteration 1 when `n_samples > 1`, so
  a cohort always runs ≥1 steady-state (LOO) iteration; single-sample convergence is
  unchanged (byte-identical). New test
  `cohort_runs_at_least_one_leave_one_out_iteration_before_converging`.
- **Mi5 — per-sample permutation-invariance untested** (the existing proptest checks only
  aggregate p̂/QUAL) (reliability). **Applied:** deterministic test
  `cohort_loo_per_sample_calls_follow_sample_permutation` — sharp, tie-free per-sample
  evidence, reversed sample order, asserts per-sample calls are equivariant. (A strict
  proptest multiset assertion was *not* added: at an exact genotype-probability tie, the
  order-dependent `expected_counts` sum can flip a call at ULP scale, so exact discrete
  invariance is not guaranteed — only the tie-free case is.)
- **Convergent doc fix (Mi, flagged by refactor_safety + smells + naming):** the
  `fill_log_indep_per_g_from` docstring had a duplicated sentence and falsely claimed the
  cohort caller passes a "per-sample equivalent" buffer (it reuses the shared
  `scratch.log_indep_per_g`). **Applied:** rewritten to the truth.

### Nits applied

- `own_counts.fill(0.0)` replaces the manual clear loop.
- Dropped the redundant `[..n_alleles]` slice on `alpha_prime.iter().sum()`.
- `accumulate_row_counts` unit test extended with an assertion touching allele index 2 (the
  `AG` row), covering the full `n_alleles == 3` inner loop.

### Not applied (with rationale)

- **Single-sample golden byte-identity test** (reliability Mi, Medium confidence): a
  hardcoded-value golden is fragile (QUAL is depth-dependent) *and* would not actually catch
  an accidental `n=1 → cohort` re-routing, because LOO with one sample is *mathematically*
  `α_species`, so a numeric parity test passes either way. The real guards are (a) the
  **structural** dispatch (`dispatch_e_step` routes `n_samples == 1` to the exact species
  path — verified byte-identical by the refactor_safety sub-agent) and (b) the **GIAB
  per-sample benchmark** in Milestone-4 validation, which is the true byte-identity gate.
- **`loo` glossary anchor / identifier rename** (naming Mi): `e_step_cohort_loo` spells out
  "leave-one-out" in its own docstring's first line; the local was renamed `raw_loo`. Kept
  the function name (documented in place).

## 5. What's good (from the sub-agents)

- The single-sample path is **provably byte-identical to HEAD** — `dispatch_e_step` routes
  `n_samples == 1` to the exact `if M::HAS_LANE_4 { e_step_simd } else { e_step }` calls
  independent of phase, so the GIAB non-regression is structural, not measured.
- Allocation-free per-sample scratch reuse (`alpha_prime` / `lgamma_alpha_prime` /
  `log_p_effective_prime` / `own_counts` resized once, overwritten per sample, kept distinct
  from the species buffers) keeps the cohort E-step off the allocator.
- `EmStepPhase` + `dispatch_e_step` centralise E-step selection behind one exhaustive match.
- Both helper extractions (`fill_log_indep_per_g_from`, `accumulate_row_counts`) are
  behaviour-identical to the pre-refactor inline bodies, with the biallelic fast path
  untouched.

## 6. Deferred (perf review, owner-approved)

- The cohort path runs scalar per-sample (no SIMD/shared-prior sharing) because the LOO prior
  is per-sample. The large-N shared-prior shortcut (when removing one of N samples is below
  numerical noise) is deferred to the performance review per the design sign-off (arch
  §10.2a / §10.5).
