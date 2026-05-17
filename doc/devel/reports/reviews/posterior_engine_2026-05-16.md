# Code Review: posterior_engine
**Date:** 2026-05-16
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** Stage 6 posterior engine — new module + integration test
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** snippet (a freshly added module).
- **Reviewed against:** working tree on `main` (uncommitted: `src/var_calling/posterior_engine.rs`, `tests/posterior_engine_integration.rs`, `src/var_calling/mod.rs` one-line export).
- **In-scope files:**
  - `src/var_calling/posterior_engine.rs`
  - `tests/posterior_engine_integration.rs`
- **Deliberately out of scope:**
  - `src/var_calling/per_group_merger.rs` — Stage 5, unchanged here. Two pre-existing clippy `needless_range_loop` warnings at lines 1786 + 1796 were observed but are not posterior-engine concerns.
  - `src/var_calling/mod.rs` — only a one-line `pub mod posterior_engine;` addition.
  - `benches/var_calling_perf.rs` — pre-existing `manual_clamp` / `needless_range_loop` warnings unrelated to this work.
- **Categories dispatched** (one sub-agent per, in parallel):
  - `reliability` — always.
  - `errors` — always.
  - `naming` — always.
  - `defaults` — public config struct with `Default` impl + eleven behavioral knobs.
  - `idiomatic` — always.
  - `refactor_safety` — always.
  - `smells` — always.
  - `extras` — hot path per the plan; stable downstream contract (`PosteriorRecord`); no public-crate / parser-of-untrusted-input concerns.
  - **Skipped:** `unsafe_concurrency` (no `unsafe`, no `async`, no shared mutability beyond a single `Arc<Config>` — itself flagged as cosmetic by the `idiomatic` agent); `tooling` (no `Cargo.toml` changes).

### 2. Verdict

**Request-changes.** Two Blocker-class test gaps leave the engine's internal-bug guards entirely unexercised; the closed-form math and the `PosteriorRecord` row-sum invariant can silently regress without any test catching the drift. Several Major findings also surface a recurring theme — the engine accepts numerically invalid configuration and reports the downstream symptom as an "internal bug" rather than the actual root cause.

### 3. Execution status

Commands run (in the project sandbox; the `./scripts/dev.sh` container wrapper is unavailable in this environment because podman is not installed):

- `cargo fmt --check` — passes (no diff), exit 0.
- `cargo build --all-targets` — passes, exit 0.
- `cargo test --all-targets --all-features` — passes, exit 0. 600 lib tests, 5 new integration tests, all pre-existing suites green.
- `cargo clippy --all-targets --all-features` (without `-D warnings`) — zero warnings on the in-scope files. Pre-existing warnings noted in §1 are out of scope.

Commands not run:

- `cargo doc --no-deps` — not run; no behavior-impacting doc cross-references in this change.
- `cargo audit` — not run; no dependency edits in this change.

Findings labeled "Needs verification": 0. Two Major findings carry `Confidence: Medium` because they depend on whether downstream layers validate first (M3, M10) — see the per-finding `Verification step`.

### 4. Open questions and assumptions

1. **Is "EM converged on success, never converged on failure" the real `EmDiagnostics.converged` contract?** Affects M11. The current shape never emits `converged: false`, so the field carries no information for any downstream consumer. Two reasonable fixes — drop the field, or model the divergence path explicitly so the type rules out the impossible.
2. **Should `PosteriorEngineConfig` validate its knobs at construction, or is validation the CLI layer's job?** Affects M3 (and indirectly the missing test in B1). If the engine is allowed to receive an unvalidated config, every "non-finite posterior" today blames an internal bug for a user-input problem. If validation is engine-side, M3 is a real bug; if it's CLI-side, M3 downgrades to a doc note.
3. **Should `f̂_C` substitute for `p̂[compound]` in the prior for chain-evident samples as well, or only chain-broken ones?** Affects how M6's compound-frequency test should be written. The implementation report explicitly takes the spec's "VCF encoding of the flag" position (always substitute), and the implementation matches; locking this in a test that explicitly references the spec section anchors the contract.
4. **Is the `1 − α_compound` Beta partner choice for `f̂_C` a temporary calibration choice or the long-term shape?** Affects M9 / M14. The plan names `α_compound` but not its Beta partner; the implementation uses `1 − α_compound` to put the prior mean at `α_compound` exactly. The user-facing API today exposes only one knob, which is convenient but invisibly couples the two.

### 5. Top 3 priorities

1. **B1 — Add the missing `NonFinitePosterior` test.** Two-line `merged_record_simple` fixture; locks in the engine's only protection against malformed upstream data.
2. **B2 — Add direct tests for `trivial_posterior_record` and reconsider its row-sum policy.** The branch fills the posterior buffer with `1.0`, which violates the `PosteriorRecord` row-sum contract for `n_genotypes > 1`. Today no Stage-5 emission can trigger that path — but the branch exists precisely as a belt-and-braces for "if it ever does," and right now it would silently emit a malformed VCF.
3. **M1 — Fix `classify_nonfinite` to distinguish `-∞` from `+∞` (and tighten the `match` to make a third state a compile error).** Convergent finding across `errors`, `idiomatic`, and `smells`; the enum is `#[non_exhaustive]` so the fix is purely additive.

### 6. Findings

#### Blocker

##### B1: src/var_calling/posterior_engine.rs:776 — `NonFinitePosterior` error path is entirely untested

- **Categories:** reliability (primary), extras (Minor), errors (Nit)
- **Confidence:** High
- **Problem:** Both branches of the `NonFinitePosterior` guard inside `e_step` — `!log_z.is_finite()` at lines 776-786 and the per-genotype `!value.is_finite()` at lines 790-800 — are unreachable from every test in `src/var_calling/posterior_engine.rs` and `tests/posterior_engine_integration.rs`. The error variant is also the only diagnostic that catches an EM that produced a malformed posterior; if it regresses (e.g. someone replaces `!log_z.is_finite()` with `log_z.is_nan()`), the engine will start emitting `NaN`-contaminated posteriors silently. The helper `classify_nonfinite` is likewise untested.
- **Why it matters:** Per the severity rubric, "absence of tests for any code path that, if broken, would produce wrong results without panicking" is Blocker. Without this test the only signal that catches malformed upstream data can silently rot.
- **Suggested fix:**
  ```rust
  #[test]
  fn e_step_returns_non_finite_posterior_when_likelihood_row_is_all_negative_infinity() {
      let record = merged_record_simple(
          7, 42, vec![b"A", b"C"], 2,
          vec![vec![f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY]],
      );
      let out = engine_for(record);
      match &out[0] {
          Err(PosteriorEngineError::NonFinitePosterior {
              chrom_id, start, end, sample_idx, ..
          }) => {
              assert_eq!(*chrom_id, 7);
              assert_eq!(*start, 42);
              assert_eq!(*end, 42);
              assert_eq!(*sample_idx, 0);
          }
          other => panic!("expected NonFinitePosterior, got {other:?}"),
      }
  }

  #[test]
  fn classify_nonfinite_returns_nan_for_nan_and_distinguishes_pos_vs_neg_infinity() {
      // After M1 is applied; before then, lock in the current behavior so any change is loud.
      assert_eq!(classify_nonfinite(f64::NAN), NonFiniteKind::NaN);
      assert_eq!(classify_nonfinite(f64::INFINITY), NonFiniteKind::PositiveInfinity);
      // Either NegativeInfinity (after M1) or PositiveInfinity (before M1) — lock the variant.
  }
  ```

##### B2: src/var_calling/posterior_engine.rs:924-975 — `trivial_posterior_record` has zero tests and may emit posteriors that violate the row-sum contract

- **Categories:** reliability (primary), errors (Minor), defaults (cross-category), idiomatic (cross-category)
- **Confidence:** High
- **Problem:** `run_em_for_record` short-circuits to `trivial_posterior_record` when the upstream record has fewer than two alleles (line 489). The trivial path returns a hand-built record with `posteriors = vec![1.0; n_samples * n_genotypes]` (every entry `1.0`, regardless of `n_genotypes`). For `n_alleles == 1`, `genotype_order(_, 1).len() == 1` and the per-sample row trivially sums to 1; for any other `n_genotypes` the per-sample row would sum to `n_genotypes`, violating the `PosteriorRecord::posteriors` docstring at line 270 ("Each per-sample row sums to 1 within floating-point tolerance"). No test exercises any of this — neither the `n_alleles == 1` path nor an `alleles.is_empty()` defensive case.
- **Why it matters:** This branch is documented as defensive against Stage 5 letting a REF-only group slip through. Without tests, its output shape can silently regress; a row-sum violation will produce a malformed VCF rather than a panic. Per the rubric, missing tests for a path that produces wrong results without panicking is Blocker.
- **Suggested fix:** Add the tests, and fix the posteriors initialisation so the row-sum invariant holds for any `n_genotypes`. At minimum:
  ```rust
  // In trivial_posterior_record (line 940):
  let mut posteriors = vec![0.0_f64; n_samples * n_genotypes];
  for s in 0..n_samples {
      // Hom-REF is at canonical position 0 — see summarise_posteriors.
      posteriors[s * n_genotypes] = 1.0;
  }

  // New unit test:
  #[test]
  fn trivial_posterior_record_returns_hom_ref_per_sample_when_only_ref_allele_present() {
      let alleles = simple_alleles(&[b"A"]);
      let record = MergedRecord {
          chrom_id: 3, start: 7, end: 7,
          alleles, ploidy: 2, n_samples: 2, n_genotypes: 1,
          scalars: vec![AlleleSupportStats::default(); 2],
          other_scalars: vec![AlleleSupportStats::default(); 2],
          chain_anchor_flags: vec![false; 2],
          log_likelihoods: vec![0.0; 2],
      };
      let pr = single_ok(record);
      assert_eq!(pr.best_genotype, vec![0, 0]);
      assert_eq!(pr.qual_phred, 0.0);
      for s in 0..pr.n_samples {
          let sum: f64 = pr.posteriors_row(s).iter().sum();
          assert!((sum - 1.0).abs() < 1e-12, "row sum = {sum}");
      }
  }
  ```

#### Major

##### M1: src/var_calling/posterior_engine.rs:1049 — `classify_nonfinite` silently mis-labels `-∞` as `PositiveInfinity`

- **Categories:** errors, idiomatic, smells (Nit)
- **Confidence:** High
- **Problem:** `classify_nonfinite` distinguishes only `NaN` vs `PositiveInfinity`. `log_sum_exp_slice` legitimately returns `f64::NEG_INFINITY` when every input is `-∞`, and `e_step` documents that case explicitly at lines 770-775. When it fires the error variant reports `PositiveInfinity` for what is actually `-∞`; the `NonFiniteKind` enum has no `NegativeInfinity` variant at all.
- **Why it matters:** The whole point of `NonFinitePosterior` is to give operators an honest fingerprint of an internal bug. Reporting the wrong sign sends triage in the wrong direction. The enum is `#[non_exhaustive]`, so the fix is additive.
- **Suggested fix:** Add the missing variant and use `match` so adding a state is a compile error rather than silently absorbed:
  ```rust
  #[derive(Debug, Clone, Copy, PartialEq, Eq)]
  #[non_exhaustive]
  pub enum NonFiniteKind {
      NaN,
      PositiveInfinity,
      NegativeInfinity,
  }

  fn classify_nonfinite(value: f64) -> NonFiniteKind {
      if value.is_nan() {
          NonFiniteKind::NaN
      } else if value == f64::INFINITY {
          NonFiniteKind::PositiveInfinity
      } else {
          NonFiniteKind::NegativeInfinity
      }
  }
  ```

##### M2: src/var_calling/posterior_engine.rs:514-519 — `n_genotypes` mismatch is only checked in `debug_assert_eq!`

- **Categories:** errors, extras (cross-category)
- **Confidence:** Medium
- **Assumptions:** Assumes Stage 5 is the only producer of `MergedRecord` today. If future callers (tests, fixtures) hand-build `MergedRecord`, the silent miscompute in release builds becomes more likely.
- **Problem:** `n_genotypes` is taken from `MergedRecord` rather than recomputed from `genotype_order(ploidy, n_alleles)`. The agreement is checked only via `debug_assert_eq!`. In release builds, a mismatched `n_genotypes` will either panic on slice OOB inside `e_step` (`log_likelihoods[sample_idx * n_genotypes..]`) or, worse, read stride-incorrect rows. The whole engine is built on the assumption that this invariant holds.
- **Why it matters:** Debug-only checks of a load-bearing invariant are exactly the silent-fallback shape the project's "No logs — promote to typed errors" memory flags. The release-build symptom is either a release-only panic or stride-incorrect reads; both fail the engine's "errors latch as typed errors" contract.
- **Suggested fix:** Promote to a typed variant, check unconditionally. Also check `log_likelihoods.len() == n_samples * n_genotypes` and `chain_anchor_flags.len() == n_samples * n_alleles` (also dereferenced via raw slicing below).
  ```rust
  if genotypes.len() != n_genotypes {
      return Err(PosteriorEngineError::GenotypeCountMismatch {
          chrom_id, start, end,
          expected: genotypes.len(),
          got: n_genotypes,
      });
  }
  ```

##### M3: src/var_calling/posterior_engine.rs:89-149 — `PosteriorEngineConfig` accepts numerically invalid knobs silently

- **Categories:** errors, reliability (out-of-range F)
- **Confidence:** Medium
- **Assumptions:** Depends on whether CLI/upstream layers validate first. If they don't, the engine misdirects triage; if they do, this downgrades to a doc note.
- **Problem:** `fixation_index_default`, `fixation_index_overrides[*]`, every `*_pseudocount` field, `max_iterations`, `convergence_threshold`, `max_gq_phred` are bare numeric fields with no constructor-side validation. Out-of-range inputs propagate quietly:
  - `fixation_index_default ∈ {-0.5, 1.5}`: `safe_ln` silently maps to `-∞` (for `≤ 0.0`) or a valid log (for `> 1.0`), producing a mathematically-defined but meaningless prior.
  - `fixation_index_default = NaN`: `safe_ln(NaN) = NaN` (because `NaN <= 0.0` is false, so `safe_ln` falls through to `x.ln()`); contaminates the prior → posteriors go NaN → caught by `NonFinitePosterior`, but the error blames "internal bug" rather than "you gave me F=NaN".
  - `ref_pseudocount = 0.0` together with zero expected counts: `total = 0` in `m_step` → `(e + a) / total = NaN`.
- **Why it matters:** Per project memory `feedback_no_logs_use_errors.md` — invariant violations should be typed errors, not silently propagated to an unrelated error site. The current shape conflates configuration errors with engine bugs and points operators at the wrong layer.
- **Suggested fix:** Add a `PosteriorEngineConfig::validate(&self) -> Result<(), PosteriorEngineError>` (or a new `ConfigValidationError`) and call it at the top of `with_config` (or in a `try_new`). Minimum coverage: `F ∈ [0, 1]`, all pseudocounts strictly positive and finite, `max_iterations > 0`, `convergence_threshold > 0`, `max_gq_phred > 0`. Per-sample F values in `fixation_index_overrides` validated entry-wise.

##### M4: src/var_calling/posterior_engine.rs:111-117,669-692 — `fixation_index_overrides` happy path is untested

- **Categories:** reliability
- **Confidence:** High
- **Problem:** `resolve_fixation_indices` has two branches: `None` (covered transitively by every test) and `Some(_)`. The mismatch error path is covered (line 1473-1493); the matching-length happy path — the entire reason the field exists — is never exercised. A regression that silently ignored `overrides` and always returned `vec![default; n_samples]` would not fail any current test.
- **Why it matters:** When the file-based per-sample F loader lands, this path is load-bearing for genotype calls in inbred lines — exactly where the project distinguishes itself from GATK.
- **Suggested fix:**
  ```rust
  #[test]
  fn fixation_index_overrides_apply_per_sample_when_length_matches_record() {
      let record = merged_record_simple(
          1, 100, vec![b"A", b"C"], 2,
          vec![vec![-2.0, 0.0, -2.0], vec![-2.0, 0.0, -2.0]],
      );
      let config = PosteriorEngineConfig {
          fixation_index_overrides: Some(vec![0.0, 1.0]),
          ..Default::default()
      };
      let pr = engine_for_with_config(record, config)
          .into_iter().next().unwrap().expect("posterior");
      assert_eq!(pr.best_genotype[0], 1);   // F=0: het allowed
      assert_ne!(pr.best_genotype[1], 1);   // F=1: het forbidden
  }
  ```

##### M5: src/var_calling/posterior_engine.rs:338 — `Arc<PosteriorEngineConfig>` is cosmetic; no sharing happens

- **Categories:** idiomatic, smells (Minor)
- **Confidence:** High
- **Problem:** `PosteriorEngine.config` is wrapped in `Arc<_>`, but the only owner is the engine itself; no `.clone()` of the Arc happens, no second owner exists, and the config never crosses thread boundaries (v1 is sequential per the plan). `config()` returns `&PosteriorEngineConfig`, `process` calls `run_em_for_record(record, &self.config)`, Debug just reborrows.
- **Why it matters:** Heap allocation + atomic refcount bump for no semantic gain. The `Arc` also misleads readers — it suggests "this might be shared" and forces the reader to grep to confirm it isn't.
- **Suggested fix:** Hold the config by value.
  ```rust
  pub struct PosteriorEngine<I> where … {
      upstream: I,
      config: PosteriorEngineConfig,
      done: bool,
  }
  ```
  If a future rayon-across-records pass needs sharing, reintroduce the `Arc` at that point with a typed reason.

##### M6: src/var_calling/posterior_engine.rs:702 + :808 + :924 — Three unjustified `#[allow(clippy::too_many_arguments)]`

- **Categories:** smells (primary), extras (Nit), refactor_safety (cross-category), idiomatic (Nit)
- **Confidence:** High
- **Problem:** `e_step` (13 args), `m_step` (9 args), and `trivial_posterior_record` (11 args) each suppress `clippy::too_many_arguments` with a bare `#[allow(...)]` and no inline justification. The project's smells rule requires every `#[allow]` to carry a comment naming the lint, why the suppression is correct here, and the removal condition. The integration test at line 170-173 demonstrates the right shape.
- **Why it matters:** Unjustified `#[allow]`s rot. The recurring `(chrom_id, start, end)` triple + per-record dims + scratch buffers point at the real fix: extract a `RecordLocus` newtype and an `EmContext` struct.
- **Suggested fix:** Introduce two small types and drop the allows:
  ```rust
  #[derive(Debug, Clone, Copy)]
  struct RecordLocus { chrom_id: u32, start: u32, end: u32 }

  struct EmContext<'a> {
      locus: RecordLocus,
      n_samples: usize,
      n_genotypes: usize,
      n_alleles: usize,
      ploidy: u8,
      genotype_allele_counts: &'a [u32],   // see M13
      log_multinomial_coeffs: &'a [f64],
      compound_mask: &'a [bool],
      pseudocounts: &'a [f64],
      fixation_indices: &'a [f64],
      compound_pseudocount: f64,
  }
  ```
  `e_step(ctx, log_likelihoods, p_hat, f_hat_compound, &mut posteriors)` and `m_step(ctx, &posteriors)` then lose every allow. The error variants also collapse `(chrom_id, start, end)` to a single `locus: RecordLocus` field with a `Display` impl that renders `chrom {chrom_id} {start}-{end}` once.

##### M7: src/var_calling/posterior_engine.rs:295-315 — `PosteriorRecord::*_row` accessors lack boundary-index tests and stride-correctness tests

- **Categories:** reliability (primary), errors (Minor — undocumented panic surface)
- **Confidence:** High
- **Problem:** All three slice helpers are public API. Existing tests touch them only for `sample_idx = 0`. No test exercises `sample_idx = n_samples - 1`, and no test confirms `scalars_row`/`chain_anchor_flags_row` use the `n_alleles` stride rather than the `n_genotypes` stride. A regression that swapped strides between the three helpers would only fail asymmetrically depending on layout.
- **Why it matters:** These helpers are how every downstream consumer reads `PosteriorRecord`. A stride bug would silently mis-align posteriors to samples, mis-attributing genotype calls and GQs without a panic.
- **Suggested fix:**
  ```rust
  #[test]
  fn posteriors_row_returns_correct_slice_for_last_sample_index() {
      let record = merged_record_simple(
          1, 100, vec![b"A", b"C"], 2,
          vec![
              vec![0.0, -50.0, -50.0],   // s0 → RR
              vec![-50.0, 0.0, -50.0],   // s1 → RA
              vec![-50.0, -50.0, 0.0],   // s2 → AA
          ],
      );
      let pr = single_ok(record);
      let last = pr.posteriors_row(pr.n_samples - 1);
      assert_eq!(last.len(), 3);
      assert!(last[2] > 0.99);
  }

  #[test]
  fn scalars_row_uses_n_alleles_stride_not_n_genotypes() {
      let record = merged_record_simple(
          1, 100, vec![b"A", b"C", b"G"], 2,
          vec![vec![0.0; 6], vec![0.0; 6]],
      );
      let pr = single_ok(record);
      assert_eq!(pr.scalars_row(1).len(), 3); // n_alleles, not n_genotypes
  }
  ```
  Also add a `# Panics` note on each `*_row` doc comment.

##### M8: src/var_calling/posterior_engine.rs:864-922 — `summarise_posteriors` boundary cases untested

- **Categories:** reliability
- **Confidence:** High
- **Problem:** The defensive `p_best.min(1.0 - f64::EPSILON)` clamp (line 901) and the `p_hom_ref > 0.0` branch driving `qual_phred = f64::INFINITY` (lines 904-908, 917-921) have inline comments documenting intent but no tests pin them down. `gq_phred_matches_closed_form` only checks the formula on a strong-RR fixture where `p_best ≠ 1.0`.
- **Why it matters:** "QUAL = +∞ when at least one sample is certainly variant" is the contract the VCF writer will rely on to gate sites. A regression from `> 0.0` to `>= 0.0` silently makes QUAL a finite garbage value at the very sites where the engine is most confident.
- **Suggested fix:**
  ```rust
  #[test]
  fn summarise_posteriors_returns_infinite_qual_when_any_sample_has_zero_hom_ref_posterior() {
      let record = merged_record_simple(
          1, 100, vec![b"A", b"C"], 2,
          vec![vec![-1000.0, -1000.0, 0.0]],
      );
      let pr = single_ok(record);
      if pr.posteriors_row(0)[0] == 0.0 {
          assert!(pr.qual_phred.is_infinite() && pr.qual_phred > 0.0);
      }
      // GQ clamps regardless.
      assert!(pr.gq_phred[0].is_finite());
      assert!(pr.gq_phred[0] <= MAX_GQ_PHRED + 1e-9);
  }
  ```

##### M9: src/var_calling/posterior_engine.rs:809-862 — `f̂_C` closed-form M-step never numerically verified

- **Categories:** reliability
- **Confidence:** High
- **Problem:** Existing compound tests (`compound_allele_emits_compound_frequency_in_posterior_record`, `chain_broken_sample_uses_likelihood_unchanged_and_still_emits_compound_frequency`) only check that `f̂_C` is `Some(_)` and in `(0, 1)`. The actual formula `(α + E[n_C]) / (α + (1 - α) + ploidy · n_samples)` is never asserted against a known input. A regression that, say, swapped `(1 - α)` with `1.0` would still pass.
- **Why it matters:** `f̂_C` is one of two cohort-frequency outputs Stage 6 was built to produce. Its closed form is small enough to assert exactly under controlled inputs.
- **Suggested fix:**
  ```rust
  #[test]
  fn m_step_compound_frequency_matches_closed_form_with_no_observed_compound_counts() {
      let record = merged_record_with_compound(
          1, 100, vec![b"A", b"C"], &[1],
          vec![vec![false, false], vec![false, false]],
          2,
          vec![vec![0.0, -50.0, -50.0], vec![0.0, -50.0, -50.0]],
      );
      let pr = single_ok(record);
      let f_c = pr.compound_frequencies[1].expect("compound 1 frequency");
      let alpha = DEFAULT_COMPOUND_ALT_PSEUDOCOUNT;
      let expected = alpha / (alpha + (1.0 - alpha) + 2.0 * 2.0);
      assert!((f_c - expected).abs() < 1e-6, "f̂_C = {f_c}, expected {expected}");
  }
  ```

##### M10: src/var_calling/posterior_engine.rs:1605-1696 — Proptest coverage limited to diploid biallelic

- **Categories:** reliability
- **Confidence:** High
- **Problem:** `proptest_likelihood_row` generates a 3-element vector exactly. Every property test inherits that strategy. The user-confirmed "polyploid via Wright–Fisher" contract is never property-tested; only one unit test (`triploid_run_succeeds_with_polyploid_hwe_prior`) touches polyploid and it only checks row sum.
- **Why it matters:** Polyploid is one of the three user-confirmed v1 contracts; core invariants (simplex sum, permutation invariance, monotonicity in `α_ref`) should hold at every (ploidy, n_alleles) combination the engine claims to support.
- **Suggested fix:** Parameterise the strategy:
  ```rust
  fn proptest_record_strategy()
      -> impl Strategy<Value = (u8, Vec<&'static [u8]>, Vec<Vec<f64>>)>
  {
      (2u8..=4u8, 2usize..=3usize, 1usize..=4usize).prop_flat_map(
          |(ploidy, n_alleles, n_samples)| {
              let alleles: Vec<&[u8]> = b"ACGT".iter().take(n_alleles)
                  .map(|b| std::slice::from_ref(b)).collect();
              let n_genotypes = genotype_order(ploidy, n_alleles).len();
              let row = prop::collection::vec(-50.0_f64..0.0_f64, n_genotypes);
              let likelihoods = prop::collection::vec(row, n_samples);
              (Just(ploidy), Just(alleles), likelihoods)
          },
      )
  }

  proptest! {
      #[test]
      fn posteriors_sum_to_one_at_any_ploidy_and_n_alleles(
          (ploidy, alleles, likelihoods) in proptest_record_strategy()
      ) {
          let record = merged_record_simple(1, 100, alleles, ploidy, likelihoods);
          let pr = single_ok(record);
          for s in 0..pr.n_samples {
              let sum: f64 = pr.posteriors_row(s).iter().sum();
              prop_assert!((sum - 1.0).abs() < 1e-9);
          }
      }
  }
  ```

##### M11: src/var_calling/posterior_engine.rs:134 + :367 + :89-132 — `Default` impl hides every behavioral knob; `PosteriorEngine::new` lets call sites silently opt into the project priors

- **Categories:** defaults
- **Confidence:** High
- **Problem:** Three converging issues in the defaults story:
  1. `impl Default for PosteriorEngineConfig` returns a config with REF pseudocount 10, SNP pseudocount 0.01, indel pseudocount 0.00125, compound pseudocount 0.001, F=0, max-iter 50, GQ cap 99 — every one of those is behavior-shaping, not behaviorally inert. There is no `PosteriorEngineConfig::with_project_defaults()` whose rustdoc enumerates the defaults inline.
  2. `PosteriorEngine::new(upstream)` is the `Client::new(url)` anti-pattern: a one-arg constructor that silently selects every default. Both the test suite and the integration test take this path.
  3. Field docs on `PosteriorEngineConfig` describe what each field is but never reference the `DEFAULT_*` constant the `Default` impl uses, so the doc and the value can drift independently.
- **Why it matters:** Calibration of these defaults is the engine's defining choice. Reviewers seeing `PosteriorEngine::new(upstream)` at a call site cannot tell which prior was chosen without jumping into this file.
- **Suggested fix:** All three at once:
  ```rust
  impl PosteriorEngineConfig {
      /// Construct the engine config with the project defaults.
      ///
      /// # Defaults
      ///
      /// | Field | Constant | Value |
      /// |---|---|---|
      /// | `convergence_threshold`     | [`DEFAULT_CONVERGENCE_THRESHOLD`]     | 1e-4 |
      /// | `max_iterations`            | [`DEFAULT_MAX_ITERATIONS`]            | 50   |
      /// | `ref_pseudocount`           | [`DEFAULT_REF_PSEUDOCOUNT`]           | 10.0 |
      /// | `snp_alt_pseudocount`       | [`DEFAULT_SNP_ALT_PSEUDOCOUNT`]       | 0.01 |
      /// | `indel_alt_pseudocount`     | [`DEFAULT_INDEL_ALT_PSEUDOCOUNT`]     | 0.00125 |
      /// | `compound_alt_pseudocount`  | [`DEFAULT_COMPOUND_ALT_PSEUDOCOUNT`]  | 0.001 |
      /// | `fixation_index_default`    | [`DEFAULT_INBREEDING_COEFFICIENT`]    | 0.0 |
      /// | `max_gq_phred`              | [`MAX_GQ_PHRED`]                      | 99.0 |
      pub fn with_project_defaults() -> Self { Self::default() }
  }
  ```
  Append `Defaults to [\`DEFAULT_*\`].` to each field doc. Rename `PosteriorEngine::new` to `with_project_default_config` (or remove it and leave only `with_config`); update the test and integration-test call sites.

#### Minor

##### Mi1: src/var_calling/posterior_engine.rs:317 — `EmDiagnostics.converged` is structurally always `true` in emitted records

- **Categories:** idiomatic
- **Confidence:** High
- **Problem:** The EM loop returns `DidNotConverge` on the negative branch, so `PosteriorRecord` is never built with `converged: false`. The field therefore carries zero information for any downstream consumer and is co-dependent with `iterations` and `final_max_delta_p`.
- **Suggested fix:** Either drop the field (if "always converged on success" is the real contract), or replace `EmDiagnostics` with an enum that explicitly models a "soft cap" outcome. See Open question 1.
  ```rust
  // Either:
  pub struct EmDiagnostics { pub iterations: u32, pub final_max_delta_p: f64 }
  // Or:
  pub enum EmOutcome {
      Converged { iterations: u32, final_max_delta_p: f64 },
      HitIterationCap { iterations: u32, final_max_delta_p: f64 },
  }
  ```

##### Mi2: src/var_calling/posterior_engine.rs:466 — `run_em_for_record` runs ~200 lines across five phases

- **Categories:** smells
- **Confidence:** High
- **Problem:** Destructure + trivial short-circuit, setup (classes, pseudocounts, genotype tables, initial estimates), EM main loop, post-loop final E-step, summarise + record assembly — all inlined. The smells rule flags functions >~75 lines or with >~4 distinct decision points.
- **Suggested fix:** Compose with M6: extract `build_em_context(...) -> EmContext<'_>` and `run_em_loop(ctx, …) -> Result<EmDiagnostics, _>`, leaving the orchestrator at 30–40 lines.

##### Mi3: src/var_calling/posterior_engine.rs:723 + :820 + :886 — Index loops over `sample_idx`/`g_idx` could be `chunks_exact` iterator pipelines

- **Categories:** idiomatic
- **Confidence:** Medium
- **Problem:** Three hot loops use `for sample_idx in 0..n_samples { let row = &arr[s*n..(s+1)*n]; ... }`. The idiomatic shape for "walk a flat array in fixed-stride chunks" is `chunks_exact(stride)` paired with `enumerate()` where needed.
- **Suggested fix:** Convert `e_step`, `m_step`, and `summarise_posteriors` to `chunks_exact`/`chunks_exact_mut` over the flat buffers. Done alongside M6 the refactor is natural.

##### Mi4: src/var_calling/posterior_engine.rs:723-736 + :820-859 — Per-iteration `Vec` allocations inside the EM loop

- **Categories:** idiomatic, errors (cross-category), smells (cross-category), extras (cross-category — pending Mi5)
- **Confidence:** High
- **Problem:** Each EM iteration allocates `p_effective`, `log_p_effective`, `log_post_unnorm`, `new_p_hat`, `new_f_hat_compound`, and `expected_counts`. All are sized by record-scoped constants and could be lifted into scratch buffers reused across iterations. The 2026-05-16 Stage 5 perf review identified the same pattern as the dominant allocator-bound bottleneck on Stage 5.
- **Why it matters:** Engine runs once per merged record × millions of records on cohort-scale runs; allocator pressure is the known killer.
- **Suggested fix:** Lift buffers into `run_em_for_record` and pass `&mut`s through the EM context (M6). Fold the `p_effective → log_p_effective` two-pass into a single pass:
  ```rust
  for (a, (slot_p, slot_log)) in p_effective.iter_mut()
      .zip(log_p_effective.iter_mut()).enumerate()
  {
      let p = if compound_mask[a] { f_hat_compound[a] } else { p_hat[a] };
      *slot_p = p;
      *slot_log = safe_ln(p);
  }
  ```
  Land alongside Mi5's benchmark to defend the change.

##### Mi5: src/var_calling/posterior_engine.rs:466 — No `criterion` performance guard for the per-record EM hot path

- **Categories:** extras
- **Confidence:** High
- **Problem:** Scope explicitly identifies this engine as the cohort hot path. Six existing benches cover Stages 1–5; `posterior_engine_perf.rs` does not exist. There is no baseline against which an allocator regression (Mi4) or an algorithmic regression would be detected.
- **Suggested fix:** Add `benches/posterior_engine_perf.rs` with at least two groups: `posterior_engine/per_record_em` (single record, `n_samples ∈ {1, 10, 100}` × `ploidy ∈ {2, 4}`) and `posterior_engine/streaming` (10k synthetic records end-to-end). Build inputs outside the timed region; register in `Cargo.toml` next to `var_calling_perf`. Annotate each `bench_function` with a `// REGRESSION THRESHOLD: N%` comment.

##### Mi6: src/var_calling/posterior_engine.rs:246 — `PosteriorRecord` has no golden test

- **Categories:** extras
- **Confidence:** High
- **Problem:** Today's only stability check is `engine_output_is_bit_identical_across_runs` — that catches non-determinism, not numeric drift. A tweak to `DEFAULT_REF_PSEUDOCOUNT`, a reorder of the final E-step, or a change to `safe_ln` on the `x == 0.0` boundary would all pass and silently change downstream VCFs.
- **Suggested fix:** Commit `tests/golden/posterior_engine/*.json` (or `.ron` — anything lossless for `f64`) capturing the full `PosteriorRecord` minus `diagnostics` for a small representative matrix: (single-sample biallelic SNP), (cohort biallelic SNP with weak het), (triploid biallelic), (multi-allelic with one compound — chain-evident and chain-broken variants). A new test loads each fixture, runs the engine, and asserts deep equality with explicit bit-level `f64` comparison. Reblessing requires an explicit env var (`UPDATE_GOLDEN=1`).

##### Mi7: src/var_calling/posterior_engine.rs:68-87 — `DEFAULT_*` pseudocount and `MAX_GQ_PHRED` constants are not source-anchored

- **Categories:** defaults
- **Confidence:** High
- **Problem:** Each pseudocount doc comment is a one-line tagline (`DEFAULT_REF_PSEUDOCOUNT` says "GATK-style" without naming the GATK parameter; the `10 : 0.01 : 0.00125 : 0.001` ratio is presented without a citation). `MAX_GQ_PHRED = 99.0` claims "standard VCF convention" but doesn't cite a spec. Without sources, a future maintainer has no basis to decide whether bumping `DEFAULT_SNP_ALT_PSEUDOCOUNT` to `0.005` is a tuning improvement or a regression.
- **Suggested fix:** Point each constant at either the GATK source it was lifted from (with the exact GATK parameter name) or the calibration plan section that picked the value. Mark TBD values explicitly:
  ```rust
  /// Dirichlet pseudocount on a SNP/MNP alt allele.
  ///
  /// Source: matches GATK HaplotypeCaller's `--heterozygosity` default
  /// of 0.001 mapped onto the `α_ref : α_alt = 10 : 0.01` shape per the
  /// project plan §"Stage 6 — posterior engine, priors". Revisit against
  /// the calibration set in doc/devel/implementation_plans/posterior_engine.md.
  pub const DEFAULT_SNP_ALT_PSEUDOCOUNT: f64 = 0.01;
  ```

##### Mi8: src/var_calling/posterior_engine.rs:380 — `PosteriorEngine::config` is `pub` but undocumented; no `# Errors` section on the `Iterator` impl

- **Categories:** idiomatic, extras
- **Confidence:** High
- **Problem:** `pub fn config(&self) -> &PosteriorEngineConfig` carries no `///` comment. Its sibling constructors at lines 367 and 372 do. The file has zero `# Errors` / `# Panics` / `# Examples` markers across any pub item, and `PosteriorEngine` returns `Result<_, PosteriorEngineError>` via its `Iterator` impl with no enumeration of error cases for callers.
- **Suggested fix:** Add a one-line doc on `config()`. Add one runnable `# Examples` doctest at the `PosteriorEngine` struct doc that shows `PosteriorEngine::new(upstream).collect::<Result<Vec<_>, _>>()?`. Add a `# Errors` paragraph (also at the struct doc) enumerating the variants the engine can latch on. Enforcing `#![deny(missing_docs)]` repo-wide is out of scope for this review.

##### Mi9: src/var_calling/posterior_engine.rs:669 — `resolve_fixation_indices` takes `&Option<Vec<f64>>`

- **Categories:** idiomatic
- **Confidence:** High
- **Problem:** Classic less-general parameter shape that `clippy::ref_option` flags. The body matches on `None`/`Some(ov)`, checks the length, then returns `ov.clone()`. The function name promises an owned `Vec<f64>`, but the input type carries an unnecessary borrow-of-Option layer.
- **Suggested fix:**
  ```rust
  fn resolve_fixation_indices(
      …,
      overrides: Option<&[f64]>,
      default: f64,
  ) -> Result<Vec<f64>, PosteriorEngineError> {
      let Some(ov) = overrides else {
          return Ok(vec![default; n_samples]);
      };
      if ov.len() != n_samples {
          return Err(PosteriorEngineError::InbreedingOverrideLengthMismatch { … });
      }
      Ok(ov.to_vec())
  }
  ```
  Caller updates from `&config.fixation_index_overrides` to `config.fixation_index_overrides.as_deref()`.

##### Mi10: src/var_calling/posterior_engine.rs:339 — `PosteriorEngine.done: bool` is participle form

- **Categories:** naming
- **Confidence:** High
- **Problem:** Bare adjective/participle boolean names are forbidden by the project rule. `done` reads fine inside `if self.done { … }` because the surrounding `if` provides the predicate; the field itself does not.
- **Suggested fix:** Rename to `is_latched` — predicate form *and* aligned with the doc-comment vocabulary ("Errors latch …" at line 331). Update the four call sites in `next()` and the Debug destructure.

##### Mi11: src/var_calling/posterior_engine.rs:532 — `genotype_allele_counts: Vec<Vec<u32>>` should be flat row-major

- **Categories:** idiomatic
- **Confidence:** High
- **Problem:** Built with one inner `Vec<u32>` per genotype, then carried through `e_step` / `m_step` as `&[Vec<u32>]`. Inner Vecs are tiny (length = `n_alleles`, typically 2–4). The rest of the module commits to flat row-major precisely to avoid one-allocation-per-row.
- **Suggested fix:** Mirror the row-major layout. `let mut genotype_allele_counts = vec![0_u32; n_genotypes * n_alleles];`; helpers walk via `chunks_exact(n_alleles)`. Compose with M6 (the slice fits naturally into `EmContext`).

##### Mi12: src/var_calling/posterior_engine.rs:385-390 — Config-invariant checks repeated on every record

- **Categories:** smells
- **Confidence:** High
- **Problem:** `process` re-checks `approximate_posterior_calculation` and `contamination.is_some()` on every record. Both are config-level invariants that cannot change after engine construction.
- **Suggested fix:** Move the checks into `with_config` and have it return `Result<Self, PosteriorEngineError>`. `new` either delegates or stays infallible (the default config statically passes both checks). When the deferred features land, the checks just go away.

##### Mi13: src/var_calling/posterior_engine.rs:1070-1119 vs tests/posterior_engine_integration.rs:16-55 — Fixture builders duplicated

- **Categories:** smells
- **Confidence:** High
- **Problem:** `simple_alleles` and `merged_record_simple` are copy-pasted between the unit-test module and the integration test (the unit-test version asserts likelihood-row length, the integration version doesn't — already a small drift).
- **Suggested fix:** Promote to a small `pub(crate)` test-support module gated by `#[cfg(any(test, feature = "test-fixtures"))]`. Alternative: accept the duplication for two copies and add a one-line note in each pointing at the other.

#### Nits

- `src/var_calling/posterior_engine.rs:789` — local `value` is the per-genotype posterior; the rule lists `value` as a forbidden generic noun. Rename to `posterior`.
- `src/var_calling/posterior_engine.rs:834` — local `total` is the Dirichlet denominator; rename to `dirichlet_denominator` (the surrounding comment already uses that abstraction level).
- `src/var_calling/posterior_engine.rs:679` — parameter `ov` is an abbreviation of `overrides`; destructure into `Some(overrides_slice)` or rename when applying Mi9.
- `src/var_calling/posterior_engine.rs:852` — local `e_n_c` is three-letter; rename to `expected_compound_count`.
- `src/var_calling/posterior_engine.rs:884` — `let hom_ref_idx = 0_usize;` could be a module-level `const HOM_REF_GENOTYPE_IDX: usize = 0;` so the convention is auditable in one place.
- `src/var_calling/posterior_engine.rs:900,916` — the Phred scale factor `-10.0` repeats; consider a `PHRED_SCALE: f64 = -10.0` const.
- `src/var_calling/posterior_engine.rs:1167` — homegrown `approx` helper repeated across tests; consider pulling `float-cmp` or `approx`.
- `src/var_calling/posterior_engine.rs:1287-1301` — `hard_iteration_cap_returns_did_not_converge_error` only destructures `max_iterations`; tighten the match to also assert `chrom_id`, `start`, `end`, `last_delta` so a regression that zeroes those fields surfaces.
- `src/var_calling/posterior_engine.rs:184-185` — `NonFinitePosterior`'s `Display` uses `{kind:?}`. Implement `Display` on the enum and use `{kind}`.
- `src/var_calling/posterior_engine.rs:781-784` — `genotype_idx: 0` hard-coded in the `log_z` non-finite branch. There is no per-genotype context at that site; consider `genotype_idx: Option<usize>`.

### 7. Out of scope observations

- `src/var_calling/per_group_merger.rs:1786 + 1796` — pre-existing `needless_range_loop` clippy warnings flagged during this work's clippy run. Not introduced by Stage 6; suggest a follow-up clippy-cleanup pass on the Stage 5 module.
- `benches/var_calling_perf.rs` — pre-existing `manual_clamp` / `needless_range_loop` warnings (same situation as above).

### 8. Missing tests to add now

Grouped by function under test. Test names use `function_returns_expected_on_condition` form; the inline test bodies are inlined in B1, B2, M4, M7, M8, M9, M10 above and not repeated here.

#### `e_step` (src/var_calling/posterior_engine.rs:702-806)
- `e_step_returns_non_finite_posterior_when_likelihood_row_is_all_negative_infinity` — input: a sample row of `log_likelihoods` uniformly `f64::NEG_INFINITY`. Bug caught: regression that turns the `!log_z.is_finite()` guard into `log_z.is_nan()`. (Body in B1.)
- `e_step_handles_partial_negative_infinity_likelihoods_without_nan_contamination` — input: `[f64::NEG_INFINITY, 0.0, f64::NEG_INFINITY]`. Bug caught: regression that drops the `k == 0 ⇒ 0.0` short-circuit and computes `0.0 * NEG_INFINITY = NaN`.

#### `m_step` (src/var_calling/posterior_engine.rs:808-862)
- `m_step_compound_frequency_matches_closed_form_with_no_observed_compound_counts` — input: 2 samples, REF + 1 compound alt, strong-RR likelihoods. Bug caught: regression to the Beta denominator. (Body in M9.)
- `m_step_p_hat_matches_dirichlet_posterior_mean_on_controlled_counts` — input: 4 samples × strong RR; `E[n_ref] = 8`, `E[n_alt] = 0`. Bug caught: drops ploidy from the denominator, or normalises by `expected_counts.sum()` instead of `(expected + pseudocount).sum()`.

#### `summarise_posteriors` (src/var_calling/posterior_engine.rs:864-922)
- `summarise_posteriors_returns_infinite_qual_when_any_sample_has_zero_hom_ref_posterior` — see M8.
- `summarise_posteriors_clamps_gq_when_best_posterior_is_exactly_one` — see M8.
- `summarise_posteriors_best_genotype_breaks_ties_by_first_occurrence` — pins the `>` vs `>=` argmax fold.

#### `trivial_posterior_record` (src/var_calling/posterior_engine.rs:924-975)
- `trivial_posterior_record_returns_hom_ref_per_sample_when_only_ref_allele_present` — see B2.
- `trivial_posterior_record_returns_empty_allele_frequencies_when_alleles_empty` — see B2.

#### `resolve_fixation_indices` (src/var_calling/posterior_engine.rs:669-692)
- `fixation_index_overrides_apply_per_sample_when_length_matches_record` — see M4.

#### `PosteriorEngine::next` and error latching (src/var_calling/posterior_engine.rs:401-422)
- `engine_emits_all_successful_records_before_latched_upstream_error` — input: `[Ok, Ok, Ok, Err, Ok]`. Bug caught: a regression where `done = true` is set on the first record.
- `engine_emits_no_records_after_did_not_converge_error_on_third_record` — input: stream where record #3 hits the iteration cap. Bug caught: regression where `done = true` is set only for `Upstream` variants and not for engine-internal errors.

#### `PosteriorRecord::*_row` (src/var_calling/posterior_engine.rs:295-315)
- `posteriors_row_returns_correct_slice_for_last_sample_index` — see M7.
- `scalars_row_uses_n_alleles_stride_not_n_genotypes` — see M7.
- `chain_anchor_flags_row_returns_n_alleles_long_slice` — see M7.

#### `PosteriorEngine::config` (src/var_calling/posterior_engine.rs:380)
- `engine_config_accessor_returns_the_value_used_at_construction` — input: `with_config` with a custom `max_iterations`; assert `eng.config().max_iterations == 7`.

#### `Debug for PosteriorEngine` (src/var_calling/posterior_engine.rs:342-360)
- `posterior_engine_debug_includes_config_and_done_fields` — assert the rendered string mentions `PosteriorEngine`, `config`, `done`.

#### `log_sum_exp_slice` (src/var_calling/posterior_engine.rs:996-1003)
- `log_sum_exp_slice_returns_neg_infinity_when_every_entry_is_neg_infinity`.
- `log_sum_exp_slice_returns_neg_infinity_for_empty_slice`.
- `log_sum_exp_slice_matches_log_sum_exp_2_on_two_element_inputs`.

#### `log_sum_exp_2` (src/var_calling/posterior_engine.rs:985-994)
- `log_sum_exp_2_is_symmetric_under_argument_swap` — proptest over finite `f64` pairs.

#### Property tests (src/var_calling/posterior_engine.rs:1605-1696)
- `posteriors_sum_to_one_at_any_ploidy_and_n_alleles` — see M10.
- `p_hat_is_a_valid_simplex_at_any_ploidy_and_n_alleles` — same strategy as above.
- `sample_permutation_preserves_p_hat_and_qual_at_any_ploidy_and_n_alleles`.

#### Misc small helpers
- `safe_ln_returns_neg_infinity_for_zero_and_negative_inputs`.
- `max_abs_diff_returns_largest_absolute_difference`.
- `genotype_counts_tallies_each_allele_index_with_multiplicity`.
- `pseudocount_for_returns_class_specific_value_from_config`.
- `classify_nonfinite_distinguishes_nan_pos_inf_neg_inf` (after M1).

#### Integration tests (tests/posterior_engine_integration.rs)
- `upstream_error_propagates_as_posterior_engine_upstream_error` — `[Ok, Err]` stream; assert second item is `Err(PosteriorEngineError::Upstream(_))`.
- `tetraploid_record_emits_posteriors_with_n_genotypes` — single tetraploid record; assert row sum = 1 and `n_genotypes == genotype_order(4, 3).len()`.

### 9. What's good

- **Module-level documentation pulls its weight.** The module doc at `src/var_calling/posterior_engine.rs:1-44` opens with the algorithmic contract (Wright–Fisher partition for the HWE-with-`F` prior at arbitrary ploidy, closed-form M-steps on `p̂` and `f̂_C`, error-latching), then names the v1 / deferred-mode split. A new reader gets the whole picture in one screen.
- **Flat row-major storage is consistent.** `PosteriorRecord.posteriors`, `scalars`, and `chain_anchor_flags` mirror Stage 5's `MergedRecord` layout (per the 2026-05-16 perf review) with matching `*_row` accessors. Layout discipline carries forward without bespoke duplication.
- **Hand-written `Debug` exhaustively destructures `Self`.** `src/var_calling/posterior_engine.rs:342-360` — `let Self { upstream: _, config, done } = self;` makes a new field a compile error and surfaces the omission deliberately. Worth keeping as a pattern.
- **Polyploid is unit-tested as a sanity case.** `triploid_run_succeeds_with_polyploid_hwe_prior` at line 1524 — explicit (ploidy=3, n_alleles=2) exercise of the Wright–Fisher branch. (Property-test extension covered by M10.)
- **Typed-error deferrals don't silently no-op.** `ApproximateModeNotYetImplemented` and `ContaminationModeNotYetImplemented` surface as distinct error variants rather than ignored config flags — matching the project memory `feedback_no_logs_use_errors.md`.

### 10. Commands to re-verify

Re-run to confirm the baseline:
- `cargo fmt --check`
- `cargo build --all-targets`
- `cargo test --all-targets --all-features`
- `cargo clippy --all-targets --all-features`

New commands the review introduces (run after applying fixes):
- The two new `#[test]` blocks from B1 and B2 must pass under `cargo test --lib --tests`.
- `cargo bench --bench posterior_engine_perf -- --baseline` once Mi5 lands (no baseline today).

### Author response convention

Address each finding by its identifier (e.g., "B2", "M5") with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer open questions from §4 first.
