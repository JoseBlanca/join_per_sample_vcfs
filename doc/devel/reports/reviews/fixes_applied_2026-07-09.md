# Fix Application Report: ssr_emission_null_toggles_2026-07-09.md

**Date:** 2026-07-09
**Source review:** `doc/devel/reports/reviews/ssr_emission_null_toggles_2026-07-09.md`
**Source state reviewed against:** branch `ssr-bic-emission` (worktree), uncommitted working tree
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

Two user decisions (recorded in §3) reshaped the fix set: **delete** the two empirically-dead
toggles (idea A `bic_null_refit`, genotyper-B `null_in_genotyper`), and **correct**
`null_hom_frac` `0.8 → 0.75` to match the silver-standard scorer. The deletion superseded six
findings by removal. The surviving feature is emission-only `null_from_homs` + its two validated
knobs + `keep_fp_control`.

### Review totals
- Blockers: 0 · Majors: 12 · Minors: 8 · Nits: 5

### Outcome totals
- Applied: 10 · Applied with adaptation: 2 · Already fixed: 0 · Deferred: 4 · Disputed: 0
- Failed validation: 0 · Blocked by context mismatch: 0 · Superseded (by deletion): 7 · Awaiting user answer: 0

### Validation summary (host — the dev container mounts the main checkout, not this worktree)
- `cargo fmt --check` (in-scope files) → **clean**. (Pre-existing out-of-scope fmt diffs left untouched; a stray whole-crate `cargo fmt` was reverted on the 4 out-of-scope files.)
- `cargo clippy --lib -- -D warnings` → **pass** (0). `cargo clippy --lib --tests` → my 3 files clean; 4 remaining warnings are **pre-existing** `freebayes_emit.rs` test code (out of scope; baseline clippy was `--lib` only).
- `cargo test --lib ssr::cohort` → **258 passed; 0 failed; 2 ignored** (was 250 → +8 new tests).
- `cargo doc --no-deps` → not re-run; the one new-code doc issue (the `[[wikilink]]`) was deleted with `null_in_genotyper`. Pre-existing link errors remain out of scope.
- `cargo audit` → not run.
- Performance check → **skipped** — no `Apply` touches the default hot path (deletions are off-by-default branches; default VCF re-verified byte-identical), so no bench baseline was warranted.

### Extra end-to-end evidence
- **Default byte-identity:** `ssr-call` (all toggles off) diff-clean vs the pre-fix baseline `cohort.ssr.vcf`.
- **M6 thread-determinism:** `PVC_SSR_NULL_FROM_HOMS=1` VCF byte-identical across `--threads 1` and `8`.
- **M7/Mi5 fail-loud:** `PVC_SSR_NULL_HOM_FRAC=nan|2.5` and `PVC_SSR_NULL_HOM_MINDEPTH=0` each abort with a typed `InvalidNullKnob` error naming the var, value, and reason.

### Unresolved high-priority findings
- None. All Majors are Applied, Applied-with-adaptation, or Superseded-by-deletion.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation |
|---|---|---|---|---|---|---|---|
| M1 | Major | `attribute_clean_homs` untested | Apply | Applied | No | em.rs | 6 tests pass |
| M2 | Major | `null_in_genotyper` untested | — | Superseded | Yes (delete) | em.rs | path deleted |
| M3 | Major | `null_from_homs` emission untested | Apply | Applied | No | em.rs | 1 test pass |
| M4 | Major | `bic_null_refit` monotonicity untested | — | Superseded | Yes (delete) | em.rs | path deleted |
| M5 | Major | `keep_fp_control` emit-invariance untested | Apply | Applied | No | driver.rs | 1 test pass (BIC+fb) |
| M6 | Major | thread-determinism untested | Apply | Applied w/ adaptation | No | em.rs | e2e T1==T8 + pinned float |
| M7 | Major | `NULL_HOM_FRAC` accepts NaN/out-of-range | Apply | Applied | No | driver.rs | fail-loud verified |
| M8 | Major | silent-fallback vs fail-loud invariant | Apply | Applied w/ adaptation | No | driver.rs | numeric fail-loud; bools deferred |
| M9 | Major | `0.8` magic / diverges from scorer `0.75` | Ask | Applied | Yes (→0.75) | em.rs | consts + provenance |
| M10 | Major | co-dependent bools → enum | — | Superseded | Yes (delete) | em.rs | one bool remains |
| M11 | Major | two dead toggles kept | Ask | Applied | Yes (delete) | em.rs, driver.rs | deleted both |
| M12 | Major | `..dev_default()` silent-default trap | Apply | Applied | No | driver.rs | comment at both sites |
| Mi1 | Minor | dup best-monomorphic loop | — | Superseded | Yes (delete) | em.rs | in deleted A block |
| Mi2 | Minor | `attribute_clean_homs` dup `attribute_locus` | Defer | Deferred | No | — | follow-up note |
| Mi3 | Minor | long function / inline probe | — | Superseded | Yes (delete) | em.rs | A block gone |
| Mi4 | Minor | toggles not logged / no precedence test | Defer | Deferred | No | — | see §5 |
| Mi5 | Minor | `MINDEPTH=0` accepted | Apply | Applied | No | driver.rs | rejected (≥1) |
| Mi6 | Minor | codename locals | Apply | Applied | No | em.rs | renamed |
| Mi7 | Minor | `m_b`/`m0` abbreviations | Apply | Applied | No | em.rs | folded into Mi6 |
| Mi8 | Minor | `null_from_homs` vs `null_in_genotyper` axis | — | Superseded | Yes (delete) | em.rs | one flag remains |
| N-fmt | Nit | fmt em.rs:1209 | Apply | Applied | No | em.rs | fmt clean |
| N-wiki | Nit | `[[wikilink]]` doc | — | Superseded | Yes (delete) | em.rs | doc deleted |
| N-waste | Nit | recompute under off-label combo | Defer | Deferred | No | — | benign; see §5 |
| N-doc | Nit | doc-comment volume | Apply | Applied | No | em.rs | shrank via deletion |
| N-name | Nit | `null_hom_frac` naming | Defer | Deferred | No | — | optional |

## 3. Questions asked and answers

1. **M11 / M10** — Delete idea A + genotyper-B, or keep them behind toggles as documented dead-ends?
   - **Answer:** Delete both.
2. **M9** — Is `null_hom_frac = 0.8` deliberate (vs the scorer's `0.75`), or an unintended divergence?
   - **Answer:** Correct to 0.75.

## 4. Per-finding log

### M11 — Delete the two dead-end toggles (bundled with M2, M4, M10, Mi1, Mi3, Mi8, N-wiki)
- **Final status:** Applied. **Review suggestion used:** Yes (deletion steps followed).
- **Implementation:** Removed the `bic_null_refit` field/doc/default, its ~46-line emission block, and its env parse; removed the `null_in_genotyper` field/doc/default, its genotyping-seed branch, the `refit_rounds` override, and its env parse. The genotyping θ/level init reverts to `let mut theta = seed.theta0; let mut level_multiplier = 1.0;` and the loop to `for _ in 0..cfg.refit_max_rounds`. The emission-B guard simplified from `if cfg.null_from_homs && !cfg.null_in_genotyper` to `if cfg.null_from_homs`.
- **Superseded by this deletion:** M2 & M4 (tests for deleted paths), M10 (only one null strategy remains → no enum needed), Mi1 (dup loop was in the A block), Mi3 (inline probe block gone), Mi8 (`null_in_genotyper` gone → `null_from_homs` stands alone), N-wiki (`[[wikilink]]` was in the deleted `null_in_genotyper` doc).
- **Validation:** default VCF byte-identical; `cargo build` clean (no unused imports/vars); 258 tests pass.
- **Files changed:** em.rs, driver.rs. **Residual risk:** None — off-by-default branches; default path unchanged.

### M1 — `attribute_clean_homs` tests
- **Final status:** Applied. Added 6 tests: `attribute_clean_homs_bins_off_modal_reads_as_slips`, `_excludes_a_plant_below_hom_frac`, `_excludes_a_plant_below_min_depth`, `_empty_pool_round_trips_to_the_prior`, `_tie_break_is_deterministic`, `_counts_far_slip_in_slipped_but_not_profile`. Adapted the reliability agent's specs to the real `CohortLocus::new`/`SampleEvidence`/`clean_params` helpers and signature `attribute_clean_homs(&locus, period, &params, &params.level_seed, min_depth, hom_frac)`. The binning test additionally pins `expected_slipped` (5.6) — the float sum whose fixed order underpins thread-determinism (partial M6).
- **Validation:** all 6 pass. **Files:** em.rs.

### M3 — `null_from_homs` leaves genotypes unchanged
- **Final status:** Applied. Added `null_from_homs_leaves_genotypes_unchanged`: runs `run_locus_em_with` on the `checkpoint_spec` sim cohort (clean homs present so the null actually shifts) with the toggle off vs on, asserts `calls` and `pi` are equal. **Validation:** pass. **Files:** em.rs.

### M5 — `keep_fp_control` cannot gate emission
- **Final status:** Applied. Added `keep_fp_control_cleans_genotypes_without_changing_model_emit_decision`: a BIC/freebayes-polymorphic locus whose lone lopsided (190/10) het `apply_fp_control` no-calls still emits (`drop == false`, `EmittedVariable`) while that het's GT is cleaned to a no-call. Covers both `EmitModel::Bic` and `EmitModel::Freebayes`. **Validation:** pass. **Files:** driver.rs.

### M6 — thread-determinism under the toggle
- **Final status:** Applied with adaptation. Rather than a bespoke multi-thread unit harness (disproportionate for a Medium finding), covered by (a) an e2e check: `PVC_SSR_NULL_FROM_HOMS=1` VCF byte-identical across `--threads 1` and `8`, and (b) pinning the exact `expected_slipped` float in the M1 binning test so a sum-reordering regression fails. **Files:** em.rs (assertion). **Residual risk:** low — determinism is structural (serial per-locus sums in fixed sample→observation order).

### M7 / Mi5 — validate the numeric knobs (fail loud)
- **Final status:** Applied. Added `SsrCallError::InvalidNullKnob { var, value, reason }` and a `parse_validated_env` helper; `PVC_SSR_NULL_HOM_FRAC` must be finite in `[0.0, 1.0]`, `PVC_SSR_NULL_HOM_MINDEPTH` must be `>= 1`. An unset var still returns the default (the only silent path), matching the `PVC_SSR_EMIT_MODEL` fail-loud precedent. **Validation:** `nan`, `2.5`, and `0` each abort with the typed error. **Files:** driver.rs.

### M8 — silent-fallback consistency
- **Final status:** Applied with adaptation. The *numeric* knobs now fail loud (M7/Mi5). The four **booleans** keep the house `is_ok_and(|v| v == "1")` pattern shared by every existing `PVC_SSR_*` boolean (`MARGINALIZED_PRIOR`, `COHORT_FP_CONTROL`, …); changing only these two would be inconsistent the other way, so that half is deferred as a repo-wide convention question (§5). **Files:** driver.rs.

### M9 — magic numbers / scorer divergence
- **Final status:** Applied. Hoisted `DEFAULT_NULL_HOM_MIN_DEPTH = 4` and `DEFAULT_NULL_HOM_FRAC = 0.75` (was 0.8) as named consts with a provenance comment pointing at `silver_standard.py`'s `MIN_DEPTH`/`HOM_FRAC`; field docs now reference the consts instead of restating literals. **Note:** correcting 0.8→0.75 slightly changes the `null_from_homs` benchmark numbers reported earlier in the emission-model comparison; the comparison report should be re-scored if those exact figures are cited going forward. **Files:** em.rs.

### M12 — `..dev_default()` silent-default trap
- **Final status:** Applied (minimal). Added a comment at both production struct-update sites (`em_cfg`, `fp_cfg`) recording that a *new* field defaults silently there (no compile error at the `..dev_default()` tail) and must be reviewed when added. The heavier exhaustive-destructure was judged out of proportion for a benign-today mechanism. **Files:** driver.rs.

### Mi6 / Mi7 — local naming
- **Final status:** Applied. In the surviving emission-B block: `data_ll_b → data_ll_clean_hom_null`, `fit_b → clean_hom_fit`, `theta_b → clean_hom_theta`, `m_b → clean_hom_level_multiplier`. (`fit0`/`m0` were in the deleted A block.) **Files:** em.rs.

### N-fmt / N-doc — formatting & doc volume
- **Final status:** Applied. `cargo fmt` on the 3 in-scope files (a stray whole-crate run was reverted on 4 out-of-scope files). Doc-comment volume shrank naturally with the dead-toggle deletion.

## 5. Deferred findings to carry forward
- **M8 (booleans half)** — the four `PVC_SSR_*` booleans keep the silent `is_ok_and(v=="1")` house pattern; making them fail-loud is a repo-wide convention decision (all `PVC_SSR_*` bools), not this diff's to set unilaterally.
- **Mi2** — factor the shared `(delta, count, parent_level) → fit` slip accumulation out of `attribute_clean_homs`/`attribute_locus` if a third caller appears; a two-caller note today.
- **Mi4** — a resolved-toggle startup dump + env-override test; consistent with the untested/unlogged existing knobs, and env-var override tests are process-global/flaky.
- **N-waste** — under `PVC_SSR_EMIT_MODEL=heuristic` + `PVC_SSR_NULL_FROM_HOMS=1` the clean-hom `data_ll` recompute runs though heuristic never reads the evidence; benign (an off-label combo), gate on `emit_model` if it ever matters.
- **N-name** — `null_hom_frac` → `null_hom_modal_frac` (optional).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** No. No `Apply` touches the default hot path — the deletions remove off-by-default branches, the renames are semantics-preserving, and the validation/consts live in startup config. The default VCF is re-verified byte-identical, so a bench baseline was not warranted.
- **Outcome:** Skipped — no Apply touched perf-sensitive default-path code.

## 10. Commands run
- `cargo build --release`
- `cargo fmt --check` / `cargo fmt -- <in-scope files>`
- `cargo clippy --lib -- -D warnings` ; `cargo clippy --lib --tests`
- `cargo test --lib ssr::cohort`
- `ssr-call` default byte-identity, T1-vs-T8 determinism, and fail-loud validation (host, on the ssr_tomato1 psp).

## 11. Command results
- `cargo build --release` → 0, clean.
- `cargo fmt --check` (in-scope) → 0, clean.
- `cargo clippy --lib -- -D warnings` → 0, clean.
- `cargo clippy --lib --tests` → in-scope clean; 4 pre-existing `freebayes_emit.rs` test warnings (out of scope).
- `cargo test --lib ssr::cohort` → 0, `258 passed; 0 failed; 2 ignored`.
- default byte-identity → identical; T1==T8 → identical; bad knobs → abort with `InvalidNullKnob`.

## 12. Notes
- The surviving experimental surface is now just `null_from_homs` (+ `null_hom_min_depth`/`null_hom_frac`) and `keep_fp_control` — both keepers per the emission-model investigation. Idea A and genotyper-B are gone from the compiled path (their findings in the emission-model comparison memory/report stand as the record).
- Correcting `null_hom_frac` to 0.75 (M9) slightly shifts the `null_from_homs` benchmark numbers vs the earlier comparison report; re-score if those exact figures are quoted downstream.
