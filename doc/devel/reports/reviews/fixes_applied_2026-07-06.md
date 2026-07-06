# Fix Application Report: ssr_interrupted_repeat_recall_2026-07-06.md

**Date:** 2026-07-06
**Source review:** `doc/devel/reports/reviews/ssr_interrupted_repeat_recall_2026-07-06.md`
**Source state reviewed against:** branch `ssr-interruptions`, `fb27e76`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 1
- Minors: 11
- Nits: (grouped, ~9)

### Outcome totals
- Applied: 8 (M1, Mi1–Mi7)
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 4 (Mi8, Mi9, Mi10, Mi11) + Nits
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

New tests added (findings + review §8 missing-tests): **9** (`ssr::cohort` 206 → 215).

### Validation summary
- `cargo fmt --check` → **0** on in-scope files (crate-wide shows only pre-existing `src/paralog/` drift + the untracked `examples/ssr_psp_seqdump.rs`).
- `cargo clippy --all-targets --all-features -- -D warnings` → **101**, but the *only* error is the pre-existing `src/vcf/writer.rs:384` (out of scope, SNP path). In-scope `ssr/cohort` + `examples/ssr_slip_dump.rs` are clippy-clean.
- `cargo test --lib` → **101**: `1601 passed; 1 failed; 3 ignored`. The single failure is the pre-existing out-of-scope `var_calling::posterior_engine::tests::larger_ref_pseudocount_cannot_increase_p_alt` — **unchanged** (no new failures; +9 tests over the 1592 baseline).
- `cargo doc --no-deps` → fails only on pre-existing links in `read_model/mod.rs` + `var_calling/diversity.rs` (out of scope); no in-scope doc error.
- `cargo audit` → not run (cargo-audit not installed).
- Performance check → **Skipped** — no `benches/` harness reaches the changed SSR genotyping code (the SSR bench was deleted). Instead, an end-to-end benchmark re-run confirmed no recall regression (below).

### Real-data safety check (M1 changes a sufficient statistic)
Re-ran `ssr-call` on the ssr_tomato1 cohort (post-fix binary): **1633 PASS / 74574 total — identical** to the pre-fix VCF. 70 of 74574 data rows differ, all per-sample **GQ ± 1** tweaks (no genotype flips, no FILTER/emission changes) — the θ_locus attribution correction taking effect on equidistant-read loci, exactly what M1 targets. No recall loss.

### Unresolved high-priority findings
- None. (Deferred items are all Minor/Nit — see §5.)

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | `nearest_called_by_sequence` unconditional align_subst (perf + equidistant slip-sign flip) | Apply | Applied | No | attribution.rs | Pass (test-first) | No |
| Mi1 | Minor | `interruption_count` over-scores first-unit interruption | Apply | Applied | No | param_estimation.rs | Pass (test-first) | No |
| Mi2 | Minor | `allele_balance` silently disables FP defence on empty support | Apply | Applied | No | vcf_out.rs | Pass | No |
| Mi3 | Minor | `candidate_level` recomputed per (obs, candidate) | Apply | Applied | No | em.rs | Pass (byte-identical) | No |
| Mi4 | Minor | Missing `// PANIC-FREE:` marker on swapped `.expect` | Apply | Applied | No | em.rs | Pass (doc) | No |
| Mi5 | Minor | `em_init::sample_eps` panic path untested | Apply | Applied | No | em_init.rs | Pass | No |
| Mi6 | Minor | `MIN_PURITY_FIT_READS` doc/value mismatch | Apply | Applied | No | prepass.rs | Pass (doc) | No |
| Mi7 | Minor | `PURITY_FACTOR_FLOOR` magnitude unsourced | Apply | Applied | No | prepass.rs | Pass (doc) | No |
| Mi8 | Minor | §5.2 defaults duplicated in prose + `dev_default()` | Defer | Deferred | No | None | N/A | Yes |
| Mi9 | Minor | Purity-fit gates non-sweepable consts | Defer | Deferred | No | None | N/A | Yes |
| Mi10 | Minor | P2.0 diagnostic harness placement | Defer | Deferred | No | None | N/A | Yes |
| Mi11 | Minor | Duplicated test fixtures across 6 modules | Defer | Deferred | No | None | N/A | Yes |
| Nits | Nit | naming / casts / tie-break dup / etc. | Defer | Deferred | No | None | N/A | Yes |

## 3. Questions asked and answers

None — no user question was necessary. The three review open questions were resolved by the chosen implementations:
1. **(M1 — is the equidistant different-length composition tie-break intended?)** Resolved by the conservative, behavior-preserving fix: composition breaks *only same-length* ties; different-length equidistant ties keep `nearest_parent`'s lowest-index rule. This makes the code match the documented "effect-neutral for the θ_locus stats" claim rather than change it silently.
3. **(Mi1 — does the first-unit-phase precondition hold?)** Resolved by removing the dependency: the measure is now phase-robust (per-phase majority), so no unenforced precondition remains.
Open question 2 (purity gates sweepable?) is carried forward with **Mi9**.

## 4. Per-finding log

### M1 — `nearest_called_by_sequence` scores composition on the length-nearest allele unconditionally
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** Both facets (hot-path DP on single-nearest slip reads; equidistant slip-sign flip) share one root cause and one fix. The conservative fix preserves the old behavior for different-length ties (matching `nearest_parent`) while keeping the intended same-length composition tie-break, so it resolves the correctness concern without inventing new policy.
- **Implementation summary:** rewrote `nearest_called_by_sequence` to (1) take the single length-nearest allele directly (no `align_subst`), and (2) use the composition (`align_subst`) tie-break **only** when ≥2 length-nearest alleles share a length; different-length equidistant ties fall to the lowest index. Updated the doc comment.
- **Review suggestion used verbatim?:** No (adapted — used a `SmallVec` of tied indices + an `all_same_length` guard rather than the review's sketch).
- **Verification performed:** test-first — added `equidistant_different_length_tie_falls_to_lower_index_not_composition`, confirmed it FAILED pre-fix (`left: Some((1, -1))`, expected `Some((0, 1))`), passed post-fix. All prior attribution tests (same-length composition, length-separated, slip-Δ) still pass.
- **Files changed:** src/ssr/cohort/attribution.rs
- **Tests added:** `equidistant_different_length_tie_falls_to_lower_index_not_composition`
- **Validation:** `cargo test --lib ssr::cohort::attribution` → 0, 10 passed. `cargo test --lib ssr::cohort` → 0, 207 passed. Byte-identity e2e (`same_length_interruption_locus_...`) → 0, passed. Benchmark: no emission change (§1 safety check).
- **Follow-up:** None. **Residual risk:** the correctness half was Medium-confidence (bounded to equidistant reads on a nuisance param); the fix restores the documented behavior, so the remaining behavior change on real data is the intended correction (70 loci, GQ ±1).

### Mi1 — `interruption_count` over-scores a first-unit interruption
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** the first-unit-anchored `seq[i % period]` measure miscounts a corrupted first unit as ~`units−1`; a phase-robust majority measure counts each interruption once and needs no unenforced precondition.
- **Implementation summary:** rewrote `interruption_count` to score each phase position against its **majority base** across the tract (not the first unit). Behavior-identical for well-phased tracts; correct for first-unit-corrupted ones. Kept the `period < 2` / `len < period` short-circuits.
- **Review suggestion used verbatim?:** No (chose the reviewer's option 2 — a phase-robust measure — over a bare `debug_assert`, since it *fixes* the defect rather than documenting it).
- **Verification performed:** test-first — `interruption_count_scores_a_first_unit_interruption_once` (`"TACACACACACA"` → 1, was 5); existing pure / one-interruption / 4279322 / period-1 assertions unchanged; added a sub-period boundary test.
- **Files changed:** src/ssr/cohort/param_estimation.rs
- **Tests added:** `interruption_count_scores_a_first_unit_interruption_once`, `interruption_count_is_zero_for_a_tract_shorter_than_one_period`
- **Validation:** `cargo test --lib ssr::cohort::param_estimation::tests::interruption` → 0, 3 passed. **Follow-up:** None.

### Mi2 — `allele_balance` silently disables the FP defence on empty support
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Implementation summary:** added a `debug_assert!` in `apply_fp_control` that a called het carries ≥2 `allele_support` entries (fails loud in debug/test if a future path forgets `fill_allele_support`); the release path keeps the defensive `1.0`.
- **Review suggestion used verbatim?:** Yes (debug_assert + the defensive-branch test).
- **Files changed:** src/ssr/cohort/vcf_out.rs · **Tests added:** `allele_balance_returns_one_when_support_absent`
- **Validation:** `cargo test --lib ssr::cohort::vcf_out` → 0, 14 passed (existing fp_control tests fill support, so the assert does not fire). **Follow-up:** None.

### Mi3 — `candidate_level` recomputed per (obs, candidate)
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Reasoning:** the level is `obs`-independent; hoisting to a per-candidate `Vec` once per sample removes `k × n_obs` redundant `powi`/`clamp` on the genotyping hot path.
- **Implementation summary:** in `compute_data_ll` and `fill_allele_support`, precompute `cand_levels` / `called_levels` before the read loop and index them inside. Byte-identical (same values).
- **Review suggestion used verbatim?:** N/A (mechanical hoist).
- **Verification performed:** the T1-vs-T4 byte-identity e2e test passes, and the benchmark's 70 differing rows are attributable to M1's attribution change alone (Mi3 changes no values) — confirming byte-identity.
- **Files changed:** src/ssr/cohort/em.rs · **Tests added:** none (covered by existing determinism + fill-support tests). **Follow-up:** None.

### Mi4 — missing `// PANIC-FREE:` marker
- **Severity:** Minor · **Final status:** Applied · **Implementation summary:** added the marker comment documenting the non-empty-`called` invariant at the `attribute_locus` `.expect`. Doc-only. **Files changed:** src/ssr/cohort/em.rs · **Validation:** compiles; covered by the cohort suite. **Follow-up:** None.

### Mi5 — `em_init::sample_eps` panic path untested
- **Severity:** Minor · **Final status:** Applied · **Implementation summary:** added `sample_eps_panics_on_a_sample_missing_its_group` (`#[should_panic(expected = "frozen sample group")]`), mirroring `em::sample_chemistry`'s regression test. **Files changed:** src/ssr/cohort/em_init.rs · **Validation:** in the 215-test cohort run. **Follow-up:** None.

### Mi6 — `MIN_PURITY_FIT_READS` doc/value mismatch
- **Severity:** Minor · **Final status:** Applied · **Implementation summary:** corrected the doc to say cell-total (not "each side"), and to explain the 50-vs-100 choice + its F2-provisional status. Doc-only. **Files changed:** src/ssr/cohort/prepass.rs · **Follow-up:** None.

### Mi7 — `PURITY_FACTOR_FLOOR` magnitude unsourced
- **Severity:** Minor · **Final status:** Applied · **Implementation summary:** documented it as a conservative (not measured) guard, revisited in F2. Doc-only. **Files changed:** src/ssr/cohort/prepass.rs · **Follow-up:** None.

### Mi8 — §5.2 defaults duplicated in prose + `dev_default()`
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** the "value in the field prose *and* in `dev_default()`" pattern is a **module-wide convention** (`RungCfg`, `G0FitCfg`, `CandidateCfg` all do it); changing it is a convention decision touching sibling configs, not a one-field edit — out of minimal-diff scope for this run. **Follow-up:** carry forward.

### Mi9 — purity-fit gates non-sweepable consts
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** whether to lift `MIN_PURITY_FIT_READS`/`PURITY_FACTOR_FLOOR` into a `PurityFitCfg` depends on open question 2 (are they meant to be swept?), which is itself tied to the deferred "refit purity from the broad genotypes" follow-up (`doc/devel/TODO.txt`). Decide when a purity sweep is scheduled. **Follow-up:** carry forward.

### Mi10 — P2.0 diagnostic harness placement
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** promoting the ~185-line `driver.rs::mod tests` harness to an example (or refactoring `measure_allele_slips` to reuse `run`'s passes) is a maintainability restructure, not a correctness fix; the harness works and the P2.0 gate is settled. **Follow-up:** carry forward (add a removal condition or promote to `examples/`).

### Mi11 — duplicated test fixtures across 6 modules
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** consolidating `ca_seq`/`pure6`/… into a shared `#[cfg(test)]` module is a cross-module refactor beyond minimal-diff scope. **Follow-up:** carry forward.

### Nits — naming / casts / tie-break duplication / etc.
- **Severity:** Nit · **Final status:** Deferred · **Reasoning:** cosmetic; out of scope for a conservative apply run. **Follow-up:** optional cleanup.

### Review §8 missing tests
- **Added (with the findings above or standalone):** `interruption_count_is_zero_for_a_tract_shorter_than_one_period`, `interruption_count_scores_a_first_unit_interruption_once` (Mi1), `equidistant_different_length_tie...` (M1), `candidate_for_sequence_returns_none_when_no_candidate_sits_at_the_length`, `fill_allele_support_skips_a_read_matching_no_called_allele`, `fit_purity_level_floors_an_extreme_contrast`, `fit_purity_level_fits_multiple_interruption_levels_through_origin`, `allele_balance_returns_one_when_support_absent` (Mi2), `sample_eps_panics_on_a_sample_missing_its_group` (Mi5).
- **Deferred:** `run_is_byte_identical_across_threads_with_a_purity_contrast` — requires constructing a confident-genotype cohort strong enough to activate `purity_level != none()` (non-trivial); the purity path's determinism is already traced (unsafe_concurrency: No findings) and `fit_purity_level`'s sorted-order reduce is unit-tested. Carry forward.

## 5. Deferred findings to carry forward
- Mi8 — §5.2 default duplication (module-wide convention decision).
- Mi9 — purity-fit gates → `PurityFitCfg` (gated on the purity-sweep decision / Open Q2).
- Mi10 — P2.0 harness → example / removal condition.
- Mi11 — shared test-fixtures module.
- §8 — `run_is_byte_identical_across_threads_with_a_purity_contrast`.
- Nits — naming (`purity_slip`, `PurityLevel`), cast consistency, tie-break duplication, `purity_slip` tuple→struct.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — no `benches/` harness reaches the changed SSR genotyping code (`grep` over `benches/` found no reference to `ssr::cohort` / the changed files; the SSR Mark-2 bench was deleted).
- **Baseline saved:** no (not applicable).
- **Outcome:** Skipped — no Apply touched perf-sensitive code covered by `benches/`. Instead, an end-to-end `ssr-call` re-run on ssr_tomato1 confirmed identical emission (1633 PASS / 74574 total) and only 70 GQ-±1 row differences from M1's correctness fix. M1 *removes* a per-slip-read pair-HMM DP, so its perf effect is a net improvement (unmeasured for lack of a bench).

## 10. Commands run
- `./scripts/dev.sh cargo test --lib ssr::cohort::attribution` (M1 test-first)
- `./scripts/dev.sh cargo test --lib ssr::cohort` (215 passed)
- `./scripts/dev.sh cargo test --lib` (1601 passed, 1 pre-existing failure)
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features`
- `./scripts/dev.sh cargo build --release --bin pop_var_caller` + `ssr-call` on ssr_tomato1

## 11. Command results
- `cargo test --lib ssr::cohort` → 0, `215 passed; 0 failed; 1 ignored`.
- `cargo test --lib` → 101, `1601 passed; 1 failed` (pre-existing out-of-scope `var_calling::posterior_engine`).
- `cargo fmt --check` → in-scope clean.
- `cargo clippy --all-targets --all-features -- -D warnings` → 101, only pre-existing `vcf/writer.rs:384`.
- `ssr-call` on ssr_tomato1 → identical emission; 70 GQ-±1 row diffs.

## 12. Notes
- All 8 Applied findings are validated; no new test failures were introduced (the one failing lib test and the clippy/doc reds are all pre-existing and out of scope, in `src/var_calling/` + `src/vcf/` + `read_model/`).
- M1 changes θ_locus attribution on equidistant reads; the benchmark confirms this only refines GQ (±1) on 70 loci with no recall loss.
