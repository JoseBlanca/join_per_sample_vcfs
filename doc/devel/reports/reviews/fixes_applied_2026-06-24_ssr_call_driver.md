# Fix Application Report: ssr_call_driver_2026-06-24.md

**Date:** 2026-06-24
**Source review:** `doc/devel/reports/reviews/ssr_call_driver_2026-06-24.md`
**Source state reviewed against:** branch `ssr-cohort`, HEAD `ce91077` (post ia→ai/doc move)
**Execution mode:** interactive
**Overall status:** Completed (Blocker + all Majors Applied; 8 Minors Deferred as follow-ups)

---

## 1. Executive summary

### Review totals
- Blockers: 1 (B1)
- Majors: 6 (M1–M6)
- Minors: 16 (Mi1–Mi16)
- Nits: grouped

### Outcome totals
- Applied: 15 (B1; M1–M6; Mi6, Mi7, Mi8, Mi9, Mi10, Mi11, Mi12, Mi15)
- Applied with adaptation: 0 (B1/M2/M5/M6 adapted the review's *suggestion detail* but met intent — noted per finding)
- Already fixed: 0
- Deferred: 8 (Mi1, Mi2, Mi3, Mi4, Mi5, Mi13, Mi14, Mi16)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary (final full gate)
- `cargo fmt --check` → 0
- `cargo clippy --all-targets --all-features -- -D warnings` → 0
- `cargo test --lib --all-features` → 0, `1287 passed; 0 failed; 2 ignored` (+7 vs the 1280 baseline)
- `cargo test --all-targets --all-features` → 101 — the **only** failure is the pre-existing `benches/psp_writer_perf.rs:386` panic (baseline; not this work). All lib + integration tests pass.
- `cargo doc --no-deps` → 0
- `cargo audit` → not run (cargo-audit not installed; diff adds no dependencies).
- Performance check: **skipped** — no bench under `benches/` references `ssr`/`cohort::`, so no Apply touches a bench-covered hot path.

### Unresolved high-priority findings
- None. Every Blocker and Major is Applied. The 8 Deferred are all Minors (refactors / bench-gated allocation levers / policy choices) — see §5.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| B1 | Blocker | present-order VCF columns | Apply (test-first) | Applied | driver.rs, vcf_out.rs, inbreeding.rs | Pass |
| M1 | Major | `"?"` contig fallback | Apply | Applied | driver.rs | Pass |
| M2 | Major | `sample_chemistry` silent defaults | Apply | Applied | em.rs | Pass |
| M3 | Major | duplicate `G₀` fallback const | Apply | Applied | param_estimation.rs, em.rs | Pass |
| M4 | Major | `partial_cmp().unwrap()` NaN-argmax | Apply | Applied | em.rs, vcf_out.rs | Pass |
| M5 | Major | unguarded contig/sample names | Apply (test-first) | Applied | driver.rs | Pass |
| M6 | Major | untested filtered-locus emit | Apply (test-only) | Applied | driver.rs (test) | Pass |
| Mi1 | Minor | duplicated attribution helper | Defer | — | — | — |
| Mi2 | Minor | `LocusModel` bundle (11-arg) | Defer | — | — | — |
| Mi3 | Minor | per-round alloc churn | Defer | — | — | — |
| Mi4 | Minor | per-locus `f_present` alloc | Defer | — | — | — |
| Mi5 | Minor | tuple bins primitive obsession | Defer | — | — | — |
| Mi6 | Minor | bare `#[allow(too_many_arguments)]` | Apply | Applied | driver.rs, em.rs | Pass |
| Mi7 | Minor | `..Default::default()` in test literals | Apply | Applied | em.rs (test) | Pass |
| Mi8 | Minor | `level_mult`→`level_multiplier` | Apply | Applied | em.rs | Pass |
| Mi9 | Minor | `FrozenParams.params`→`chemistry` | Apply | Applied | driver.rs | Pass |
| Mi10 | Minor | stale "sweep is serial" comment | Apply | Applied | driver.rs (test comment) | Pass |
| Mi11 | Minor | `EmCfg.inbreeding_f` test-only | Apply | Applied | em.rs (doc) | Pass |
| Mi12 | Minor | unsurfaced threads/queue_depth coercion | Apply | Applied | driver.rs (doc) | Pass |
| Mi13 | Minor | `from_utf8_lossy` alleles | Defer | — | — | — |
| Mi14 | Minor | `QUAL=.` for variable-but-zero locus | Defer | — | — | — |
| Mi15 | Minor | refit non-convergence untested | Apply (test-only) | Applied | em.rs (test) | Pass |
| Mi16 | Minor | once-per-run `level_per_group.clone()` | Defer | — | — | — |

## 3. Questions asked and answers

None — every Applied finding had one clearly-correct implementation path; the Deferred set is held for design/policy decisions (Mi13/Mi14) or a missing SSR bench (Mi3/Mi4), not blocked on a question.

## 4. Per-finding log

### B1 — present-order VCF columns
- **Severity:** Blocker
- **Initial decision:** Apply (test-first)
- **Final status:** Applied
- **Reasoning:** Verified against current code: `format_vcf_record` iterated `call.calls` (present-only) producing `present_count` columns; the cursor yields `Ok(None)` for an absent sample and the merger sparse-omits only all-absent loci, so partial-coverage loci reach the formatter. Confirmed reachable by a test that demonstrably failed (15 cols vs 21).
- **Implementation summary:** `format_vcf_record` now takes `n_samples` and builds `sample_fields` dense at width `n_samples`, placing each present call at `locus.present[k]` and leaving absent samples as the `./.:.:.` placeholder. `n_samples` threaded from `run` → `write_genotyped_chunk` → `genotype_locus` → `format_vcf_record`. `genotype_locus` gained a justified `#[allow(clippy::too_many_arguments)]` (8 args; pre-empts Mi6 for this fn).
- **Review suggestion used verbatim?:** No (adapted — review sketched a `Vec<Option<&SampleCall>>`; used a `Vec<String>` placeholder fill keyed by `present`, fewer allocations).
- **Adaptation:** placeholder-string fill rather than an intermediate `Option` vector.
- **Verification performed:** new test `run_emits_dense_sample_columns_for_a_partial_coverage_locus` failed pre-fix (15 vs 21 cols), passes post-fix; full `ssr::cohort` suite green.
- **Files changed:** `src/ssr/cohort/vcf_out.rs`, `src/ssr/cohort/driver.rs`, `src/ssr/cohort/inbreeding.rs` (test call-site arg).
- **Tests added or modified:** `run_emits_dense_sample_columns_for_a_partial_coverage_locus` (new); `format_vcf_record` call-sites in vcf_out/inbreeding tests updated for the new arg.
- **Validation:**
  - `cargo test --lib --all-features ssr::cohort` → 0, `143 passed; 0 failed`
  - `cargo fmt --check` → 0
  - `cargo clippy --lib --all-features -- -D warnings` → 0
- **Follow-up:** None.
- **Residual risk:** None — `present[k] < n_samples` holds by construction (present indices are cohort sample ids `< n_samples`).

### M1 — `"?"` contig fallback
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The silent `"?"` fallback masked a merger-invariant break with a structurally-valid-but-wrong VCF row. The merger hard-errors `UnknownCatalogChrom` upstream so the id is always in range — making this a "cannot happen" path; the crate convention for broken internal invariants is a loud panic (cf. the ploidy `assert_eq!`). Chose panic over threading a `Result` to keep the change minimal and `genotype_locus`'s `Option` shape (the emit/drop signal) intact.
- **Implementation summary:** replaced `.unwrap_or("?")` with `.unwrap_or_else(|| panic!(...))` naming the bad `chrom_id` and the table size, with a `// PANIC-FREE:` comment citing the merger guarantee.
- **Review suggestion used verbatim?:** Yes (the review's minimal panic form).
- **Adaptation:** None.
- **Verification performed:** `ssr::cohort` suite green; no test constructs an out-of-range `chrom_id`, so the panic path is not exercised (it is provably unreachable given the merger contract).
- **Files changed:** `src/ssr/cohort/driver.rs`
- **Tests added or modified:** None (unreachable-by-contract path; a test would need a deliberately corrupt `CohortLocus` the merger cannot produce).
- **Validation:**
  - `cargo test --lib --all-features ssr::cohort` → 0, `143 passed`
  - `cargo fmt --check` → 0
  - `cargo clippy --lib --all-features -- -D warnings` → 0
- **Follow-up:** None.
- **Residual risk:** None.

### M2 — `sample_chemistry` silent defaults
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The three `unwrap_or` fallbacks (group-0 / ε=0.01 / level baseline 0.05) silently genotype a sample on fabricated chemistry if the frozen-`ParamSet` density invariant breaks — the exact silent-cohort-default that decision E + `UnresolvedSamples` exist to prevent. Made the lookups total so a broken invariant fails loud, matching the crate's no-silent-default style. Verified all current callers pass dense vectors (`clean_params`/`build_param_set`), so the change is behaviour-preserving on the supported path.
- **Implementation summary:** replaced the three `.get(...).unwrap_or(...)` with `.get(...).expect(...)` carrying decision-E messages; added a `// PANIC-FREE:` doc paragraph. Removed the now-unused inline `SampleGroupId(0)` fallback path.
- **Review suggestion used verbatim?:** No (review offered direct `[]` indexing or `.expect`; chose `.expect` for clearer diagnostics).
- **Adaptation:** `.expect` with domain messages rather than bare `[]` index.
- **Verification performed:** added `should_panic` test pinning the loud-failure; full `ssr::cohort` suite green.
- **Files changed:** `src/ssr/cohort/em.rs`
- **Tests added or modified:** `sample_chemistry_panics_on_a_sample_missing_its_group` (new, `#[should_panic(expected = "frozen sample group")]`).
- **Validation:**
  - `cargo test --lib --all-features ssr::cohort::em` → 0, `17 passed`
  - `cargo fmt --check` → 0
  - `cargo clippy --lib --all-features -- -D warnings` → 0
- **Follow-up:** None.
- **Residual risk:** None — the panic fires only on a broken decision-E invariant, which `build_param_set` rejects upstream.

### M3 — duplicate `G₀` fallback const
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `em::period_decay`'s `FALLBACK { p: 0.5 }` was an independent copy of `G0FitCfg::fallback_p`'s `0.5`; the driver's backfill doc claims `fallback_p` is "the single fallback source of truth," which the EM literal silently contradicts for any unobserved period. Promoted to one named const.
- **Implementation summary:** added `pub(crate) const DEFAULT_G0_FALLBACK_P: f64 = 0.5` in `param_estimation.rs`; `G0FitCfg::dev_default` and `em::period_decay`'s `FALLBACK` both reference it. Imported the const in `em.rs`.
- **Review suggestion used verbatim?:** Yes (the review's named-const form).
- **Adaptation:** None.
- **Verification performed:** value unchanged (still `0.5`), so behaviour-identical; `ssr::cohort` suite green.
- **Files changed:** `src/ssr/cohort/param_estimation.rs`, `src/ssr/cohort/em.rs`
- **Tests added or modified:** None (pure single-source-of-truth refactor; existing `period_decay`/`G0FitCfg` tests cover the value).
- **Validation:**
  - `cargo test --lib --all-features ssr::cohort` → 0, `144 passed`
  - `cargo fmt --check` → 0
  - `cargo clippy --lib --all-features -- -D warnings` → 0
- **Follow-up:** None.
- **Residual risk:** None.

### M4 — `partial_cmp().unwrap()` NaN-argmax on the parallel emit path
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Both `final_calls` (em.rs) and `site_qual` (vcf_out.rs) selected an argmax via `partial_cmp(...).unwrap()`, which panics a rayon worker on a NaN with no context. `total_cmp` is a total order (NaN-safe) and removes the `unwrap`; for the established finite invariants it gives the identical argmax, so behaviour is preserved (the genotype/QUAL-asserting tests still pass).
- **Implementation summary:** `final_calls`: `partial_cmp(...).unwrap()` → `total_cmp(...)` and the trailing `.unwrap()` → `.expect("a PASS locus enumerates ≥1 genotype")`, with a `// PANIC-FREE:` comment. `site_qual`: `partial_cmp(...).unwrap()` → `total_cmp(...)` (empty handled by `.unwrap_or(ref_idx)`), with a comment.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** all genotype/QUAL-asserting `ssr::cohort` tests pass unchanged; grep confirmed these were the only two `partial_cmp().unwrap()` sites in scope.
- **Files changed:** `src/ssr/cohort/em.rs`, `src/ssr/cohort/vcf_out.rs`
- **Tests added or modified:** None (existing tests pin the argmax; the change removes a panic on an unreachable NaN).
- **Validation:**
  - `cargo test --lib --all-features ssr::cohort` → 0, `144 passed`
  - `cargo fmt --check` → 0
  - `cargo clippy --lib --all-features -- -D warnings` → 0
- **Follow-up:** None.
- **Residual risk:** None — `total_cmp`/`partial_cmp` agree for finite/−∞ inputs; ties keep the first element under `max_by` in both.

### M5 — unguarded contig/sample names
- **Severity:** Major
- **Initial decision:** Apply (test-first)
- **Final status:** Applied
- **Reasoning:** Contig names (untrusted `.ssr.psp` headers) and sample names (input basenames) were written verbatim into the header; the writer doc claims the caller validates VCF-cleanliness but only uniqueness was checked. A tab/newline/structured char would silently corrupt the VCF. Added a typed-error validation pass before the burn-in (fails fast, before any output).
- **Implementation summary:** new `SsrCallError::InvalidVcfName { kind, name, reason }`; `check_vcf_safe_names` (called in `run` right after `check_unique_sample_names`) + `vcf_name_violation` helper rejecting `\t`/`\n`/`\r` in any name and `,`/`<`/`>` additionally in contig names; added the `ParsedChromosome` import.
- **Review suggestion used verbatim?:** No (review suggested validate-or-escape; implemented reject-with-typed-error, the loud-failure style the crate prefers).
- **Adaptation:** typed error rather than escaping; validation lives in `run` (the writer doc already delegates this to the caller).
- **Verification performed:** test-first via direct unit tests on `check_vcf_safe_names` (a run()-level test is entangled with catalog/contig-name reconciliation, so the function is tested directly per the review's "assert the specific error variant").
- **Files changed:** `src/ssr/cohort/driver.rs`
- **Tests added or modified:** `check_vcf_safe_names_rejects_a_tab_in_a_sample_name`, `check_vcf_safe_names_rejects_structured_chars_and_newlines_in_a_contig`, `check_vcf_safe_names_accepts_clean_names` (all new).
- **Validation:**
  - `cargo test --lib --all-features ssr::cohort` → 0, `147 passed`
  - `cargo fmt --check` → 0
  - `cargo clippy --lib --all-features -- -D warnings` → 0
- **Follow-up:** Mi13 (allele-byte UTF-8 validation) is the remaining untrusted-output gap, deferred.
- **Residual risk:** Low — the char set rejected is the VCF-structural minimum; a stricter VCF 4.4 contig-ID grammar could be layered later if needed.

### M6 — untested filtered-locus emit branch
- **Severity:** Major
- **Initial decision:** Apply (test-only)
- **Final status:** Applied
- **Reasoning:** The `candidates.admit != Pass` emit branch (filtered locus kept with its reason) was uncovered end-to-end; a widened/inverted drop guard would silently truncate the VCF with no failing test. Added a `run`-level regression test.
- **Implementation summary:** new test driving a cohort whose locus B is covered only by sample 0 at cohort depth 3 (< the `min_cohort_depth = 10` floor) → `Admission::LowDepth`, asserting it is emitted with `FILTER == lowDepth`, dense `9 + 12` columns all `./.:.:.`, alongside the PASS variant at locus A. Production code unchanged — the branch was already correct; this closes the coverage gap (and incidentally pins B1's dense-row fix for filtered loci).
- **Review suggestion used verbatim?:** No (adapted — used a single-sample partial-coverage locus B to hit the `LowDepth` floor reliably given `reads_for(&[])` depth 3 × the 10-floor).
- **Adaptation:** single-coverage `LowDepth` fixture rather than an all-thin locus (which would exceed the cohort-depth floor).
- **Verification performed:** new test passes; full `ssr::cohort` suite green.
- **Files changed:** `src/ssr/cohort/driver.rs` (test only).
- **Tests added or modified:** `run_emits_a_filtered_locus_with_its_reason` (new).
- **Validation:**
  - `cargo test --lib --all-features ssr::cohort` → 0, `148 passed`
  - `cargo fmt --check` → 0
  - `cargo clippy --lib --all-features -- -D warnings` → 0
- **Follow-up:** None.
- **Residual risk:** None.

### Mi8 — `level_mult` → `level_multiplier` & Mi6 — bare `#[allow(too_many_arguments)]`
- **Severity:** Minor (both)
- **Initial decision:** Apply (both)
- **Final status:** Applied (both)
- **Reasoning:** Bundled — both are small mechanical em.rs/driver.rs cleanups touching the same `compute_data_ll`/`run_locus_em_with` region (`mult` is a non-standard abbreviation while the spelled form is used everywhere else; the three remaining bare `#[allow]` lacked the project-mandated justification comment — `genotype_locus`'s allow already got one in B1).
- **Implementation summary:** Mi8 — renamed the `level_mult` local/param/`new_level_mult` to `level_multiplier`/`new_level_multiplier` across `run_locus_em_with`, `compute_data_ll`, and the doc comments (left `refit_level_multiplier`/`level_multiplier_recovers…` untouched). Mi6 — added a one-line justification comment above the `#[allow(clippy::too_many_arguments)]` on `compute_data_ll`, `run_locus_em_with` (em.rs) and `write_genotyped_chunk` (driver.rs).
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** rename is symbol-local (compiler-checked); `ssr::cohort` suite + clippy green.
- **Files changed:** `src/ssr/cohort/em.rs`, `src/ssr/cohort/driver.rs`
- **Tests added or modified:** None.
- **Validation:**
  - `cargo test --lib --all-features ssr::cohort` → 0, `148 passed`
  - `cargo fmt --check` → 0
  - `cargo clippy --lib --all-features -- -D warnings` → 0
- **Follow-up:** None.
- **Residual risk:** None.

### Mi9 — `FrozenParams.params` → `chemistry`
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `frozen.params` (a bare generic noun beside the descriptive `f_per_sample`/`level_per_group`) read as a tautology; the field's doc already calls it "Frozen chemistry."
- **Implementation summary:** renamed the `FrozenParams.params` field to `chemistry`; updated the constructor in `run` and the two `&frozen.params` reads in `genotype_locus`.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** symbol-local rename (compiler-checked); suite + clippy green.
- **Files changed:** `src/ssr/cohort/driver.rs`
- **Tests added or modified:** None.
- **Validation:**
  - `cargo test --lib --all-features ssr::cohort` → 0, `148 passed`
  - `cargo fmt --check` → 0
  - `cargo clippy --lib --all-features -- -D warnings` → 0
- **Follow-up:** None.
- **Residual risk:** None.

### Mi7 — `..Default::default()` in `LocusSlipFit` test literals
- **Severity:** Minor — **Final status:** Applied
- **Reasoning/implementation:** the two `refit_level_multiplier_collapses_to_one…` literals used `..Default::default()`, which would silently test a future `LocusSlipFit` field at its default; spelled `profile: SlipProfile::default()` explicitly (zero behaviour change, makes a future field a compile error).
- **Files changed:** `src/ssr/cohort/em.rs` (test). **Tests:** existing test edited.
- **Validation:** `cargo test --lib ssr::cohort` → `149 passed`; fmt/clippy 0.

### Mi10 — stale "sweep is serial" comment
- **Severity:** Minor — **Final status:** Applied
- **Reasoning/implementation:** the `run_is_byte_identical_across_thread_counts` comment said "the sweep is serial"; post-J it is chunk-parallel. Reworded to describe the order-preserving `par_iter` collect over a pure per-locus map.
- **Files changed:** `src/ssr/cohort/driver.rs` (test comment). **Tests:** none.
- **Validation:** as above.

### Mi11 — `EmCfg.inbreeding_f` is test-wrapper-only
- **Severity:** Minor — **Final status:** Applied
- **Reasoning/implementation:** documented on the field that it seeds only the `run_locus_em` convenience wrapper; the production `run_locus_em_with` path takes explicit per-sample `F` and never reads it (a latent misuse trap, now flagged in the doc).
- **Files changed:** `src/ssr/cohort/em.rs` (doc). **Tests:** none.
- **Validation:** as above.

### Mi12 — unsurfaced `threads`/`queue_depth` coercions
- **Severity:** Minor — **Final status:** Applied
- **Reasoning/implementation:** documented on the `SsrCallConfig` fields that `threads == 0` is coerced to a single-threaded pool and `queue_depth == 0` is the unset sentinel replaced by `DEFAULT_SWEEP_CHUNK` (output is chunk-size-invariant). Doc-only; no runtime logging added (kept minimal — the values are visible at the call site and in the docs).
- **Files changed:** `src/ssr/cohort/driver.rs` (doc). **Tests:** none.
- **Validation:** as above.

### Mi15 — refit non-convergence untested
- **Severity:** Minor — **Final status:** Applied
- **Reasoning/implementation:** added `capped_refit_returns_calls_consistent_with_the_final_round` — `refit_max_rounds: 1` recovers the truth at high depth and agrees with the converged (`refit_max_rounds: 3`) genotypes, pinning that a capped/non-converged run returns the last-round-recomputed calls, never stale pre-loop calls.
- **Files changed:** `src/ssr/cohort/em.rs` (test). **Tests:** `capped_refit_returns_calls_consistent_with_the_final_round` (new).
- **Validation:** `cargo test --lib ssr::cohort` → `149 passed`; fmt/clippy/doc 0.

*Commit grouping note:* Mi7/Mi11/Mi15 all live in `em.rs` and Mi10/Mi12 in `driver.rs`; since interactive `git add -p` is unavailable to split a single file across commits, the `em.rs` trivials share one commit and the `driver.rs` trivials another. Each finding is logged separately above.

## 5. Deferred findings to carry forward

- **Mi1** — extract the shared read-to-nearest-allele attribution helper (3–4 call sites across `em.rs`/`prepass.rs`/`vcf_out.rs`). Cross-file refactor; tie-breaks agree today, so it is a maintainability/drift cost, not a bug. Best done alongside the planned soft-split change.
- **Mi2** — bundle `compute_data_ll`'s 11-arg signature into a `LocusModel`. Signature refactor; the justification comments (Mi6) cover the lint in the interim.
- **Mi3** — hoist `compute_data_ll`'s per-round `Vec<Vec<f64>>`/`obs_qr` allocations into reused scratch. Allocation lever the review says to bench-gate; **no SSR bench exists**, so deferred until one does (and the wall cost is measured).
- **Mi4** — per-locus `f_present` allocation in the parallel sweep → `map_init` scratch. Same bench-gated rationale as Mi3; small (× locus count).
- **Mi5** — replace the `(u16,u64,u64)`/`(u16,u64)` tuple bins with named structs. Multi-site refactor (`add_bin`/`merge_sample_stats`/`fit_level`/`add_allele_copies`); cosmetic, low value relative to churn.
- **Mi13** — validate allele bytes are ACGTN/ASCII (vs `from_utf8_lossy`). A different untrusted-boundary than M5 (allele tracts, not names); needs a design choice on where to validate. Low real-world likelihood.
- **Mi14** — `QUAL=.` vs `0.0` for a variable-but-zero-confidence emitted locus. The review flags this as an "either/or" **policy choice** (`.` is valid VCF per spec §4.5); deferred for a product decision rather than guessed.
- **Mi16** — restructure the once-per-run `grouped.level_per_group.clone()` to a move. One-time burn-in cost; negligible, restructure not worth the churn now.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** No — no bench under `benches/` references `ssr`/`cohort::` (`grep -rl "ssr\|cohort::" benches/` empty), so no Apply changed a bench-covered hot path.
- **Baseline saved:** N/A.
- **Outcome:** Skipped — no Apply touched perf-sensitive code.

## 10. Commands run

- `cargo test --lib --all-features ssr::cohort` (per finding)
- `cargo fmt --check` / `cargo fmt`
- `cargo clippy --lib --all-features -- -D warnings` (per finding)
- `cargo clippy --all-targets --all-features -- -D warnings` (final)
- `cargo test --lib --all-features` (final)
- `cargo test --all-targets --all-features` (final)
- `cargo doc --no-deps` (final)

## 11. Command results

- `cargo fmt --check` → 0
- `cargo clippy --all-targets --all-features -- -D warnings` → 0
- `cargo test --lib --all-features` → 0, `1287 passed; 0 failed; 2 ignored`
- `cargo test --all-targets --all-features` → 101 (pre-existing `psp_writer_perf` bench panic only)
- `cargo doc --no-deps` → 0

## 12. Notes

- All Applied fixes are behaviour-preserving on the supported path except the intended ones: B1 changes the data-row width (the bug fix); M1/M2 convert silent fallbacks to loud panics (broken-invariant paths only); M5 adds a new typed error for malformed names; M4 swaps `partial_cmp`→`total_cmp` (identical argmax for finite inputs).
- Two commits group same-file trivials (`em.rs`: Mi7/Mi11/Mi15; `driver.rs`: Mi10/Mi12) because interactive `git add -p` is unavailable to split one file across commits; each finding is logged individually above.
- The deferred Minors are recorded as `Open:` items in `PROJECT_STATUS.md`'s Driver-wiring block.
