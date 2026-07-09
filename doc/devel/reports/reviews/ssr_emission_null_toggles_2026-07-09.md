# Code Review: ssr-emission-null-toggles
**Date:** 2026-07-09
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** Uncommitted working-tree diff on branch `ssr-bic-emission` — four env-gated experimental toggles on the SSR Stage-2 emission path (`ssr-call`)
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** PR-like uncommitted diff (3 files, +246/−5).
- **Reviewed against:** branch `ssr-bic-emission` (worktree `/Users/jose/devel/pop_var_caller-ssr-bic-emission`), working tree as-provided.
- **In-scope files:**
  - [em.rs](../../../../src/ssr/cohort/em.rs) — `EmCfg` fields + defaults; `null_in_genotyper` genotyper path; idea-B emission `data_ll_b` block; idea-A `bic_null_refit` block; new `attribute_clean_homs`.
  - [driver.rs](../../../../src/ssr/cohort/driver.rs) — env parsing in `run`; `decide_emission` BIC/Freebayes `keep_fp_control` calls.
  - [vcf_out.rs](../../../../src/ssr/cohort/vcf_out.rs) — `FpControlCfg::keep_fp_control` field + default.
- **Deliberately out of scope:** all other files; pre-existing fmt/doc issues in `examples/ssr_psp_seqdump.rs`, `src/paralog/locus_score.rs`, `src/var_calling/paralog_filter/write_pass.rs`.
- **Categories dispatched:** reliability, errors, naming, defaults, idiomatic, refactor_safety, module_structure, smells, extras (9). `unsafe_concurrency` skipped (new code adds no unsafe/Arc/Mutex/channels/async; thread-determinism folded into reliability). `tooling` handled by the orchestrator from the verification output.

The four toggles: **(1)** `PVC_SSR_KEEP_FP_CONTROL` (keeper) — run per-sample `apply_fp_control` no-call on emitted GTs under BIC/freebayes without gating the site decision; **(2)** `PVC_SSR_BIC_NULL_REFIT` (dead end per author) — BIC monomorphic model fits its own stutter rate; **(3)** `PVC_SSR_NULL_FROM_HOMS` (keeper) + knobs `_MINDEPTH`/`_FRAC` — per-locus stutter null from clean-hom plants (`attribute_clean_homs`), emission-only; **(4)** `PVC_SSR_NULL_GENOTYPER` (dead end per author) — same null in the genotyping EM.

## 2. Verdict

**Approve-with-changes.** No Blockers. The logic is correct and its four stated invariants were verified by code-tracing (below). The changes are Major-worthy only on hardening/coverage/design-hygiene grounds: zero tests for five new runtime-reachable paths, one real input-validation foot-gun, and two known-dead toggles plus a co-dependent-bool cluster that should be simplified before merge.

## 3. Execution status

Commands run on the **host** (the dev container mounts the main checkout, not this worktree):
- `cargo fmt --check` → **fail**; in-scope hit at `em.rs:1209` (`attribute_clean_homs` `parent_level` over-wrap). Remaining diffs are pre-existing out-of-scope files.
- `cargo clippy --lib -- -D warnings` → **pass** (clean, new code included).
- `cargo test --lib ssr::cohort` → **pass**, `250 passed; 0 failed; 2 ignored` — *same count as before the diff → zero new tests*.
- `cargo doc --no-deps --lib` → 7 unresolved-link errors, all pre-existing/out-of-scope; the new `[[wikilink]]` doc comment (em.rs:124) renders as literal text.
- `cargo audit` → not run.

Findings labeled "Needs verification": 0 (all Major/Minor findings are High/Medium confidence against read code).

**Four intent claims verified by tracing (extras + naming + refactor_safety converge):** `bic_null_refit` mutates only `evidence.ln_monomorphic`, read only by the BIC arm → no freebayes leak; `keep_fp_control` reads `variable` before `apply_fp_control` → emit decision invariant, GT-only; `null_from_homs` feeds only the emission evidence → GTs untouched; all-off byte-identity holds structurally (`data_ll_emit == &data_ll`, `(seed.theta0, 1.0)`, `refit_rounds == cfg.refit_max_rounds`).

## 4. Open questions and assumptions

1. **Disposition of the two dead-end toggles** (affects M11, M10, Mi1, Mi3, Mi6, Mi7, Mi8, and test findings for the deleted paths): the author found `bic_null_refit` (idea A) and `null_in_genotyper` (genotyper-B) empirically worthless. Delete them, or keep them behind toggles as documented dead-ends? Deletion resolves roughly half the findings by removal.
2. **`null_hom_frac = 0.8` vs the silver scorer's `HOM_FRAC = 0.75`** (affects M9): deliberate tightening, or an unintended transcription divergence to be corrected to 0.75?

## 5. Top 3 priorities

1. **M1–M6 — zero test coverage for five new runtime-reachable statistical paths** (reliability). They fail *silently* (wrong numbers, no panic) when broken. Ready-to-use tests supplied in §8.
2. **M11 + M10 — delete the two dead-end toggles, collapse the co-dependent-bool cluster** (smells). Resolves the largest chunk of design-hygiene findings and shrinks the diff.
3. **M7 — `PVC_SSR_NULL_HOM_FRAC` accepts NaN/negative, silently pooling every heterozygote into the stutter null** (errors). A real, silent scientific-correctness foot-gun on a research knob.

## 6. Findings

### Major

**M1: [em.rs:1177](../../../../src/ssr/cohort/em.rs#L1177) — `attribute_clean_homs` has zero tests.** *(reliability, High)* A new non-trivial pure function (per-donor length histogram, modal pick, `min_depth`/`hom_frac` gates, `MAX_SLIP`-truncated slip binning, float `expected_slipped` sum) with no coverage of any input class. Its sibling `attribute_locus` has three dedicated tests. Wrong modal pick / gate / binning → silently wrong stutter null feeding every downstream model. **Fix:** the six challenge tests in §8.

**M2: [em.rs:548](../../../../src/ssr/cohort/em.rs#L548) — `null_in_genotyper` genotyper path untested.** *(reliability, High)* The `refit_rounds = 0` + clean-hom-seed branch has no test that it changes GTs, or that an empty clean-hom pool reproduces the frozen-seed result. A regression making it a silent no-op passes all 250 tests. **Fix:** `null_in_genotyper_with_empty_clean_hom_pool_matches_frozen_seed` + an end-to-end shoulder-het resolution test. *(Moot if OQ1 = delete.)*

**M3: [em.rs:654](../../../../src/ssr/cohort/em.rs#L654) — `null_from_homs` emission path untested.** *(reliability, High)* No test covers the `data_ll_b.as_deref().unwrap_or(&data_ll)` selection (the load-bearing line for "emission-only, GTs untouched"), the `&& !null_in_genotyper` skip, or that GTs stay put while the evidence moves. **Fix:** `null_from_homs_leaves_genotypes_unchanged_but_shifts_emission_evidence`.

**M4: [em.rs:687](../../../../src/ssr/cohort/em.rs#L687) — `bic_null_refit` monotonicity invariant untested.** *(reliability, High)* The `ln_monomorphic.max(ln_mono_refit)` doc-claims the null "only ever gets stronger, never manufactures a false positive"; no regression test pins it. **Fix:** `bic_null_refit_never_lowers_ln_monomorphic`. *(Moot if OQ1 = delete.)*

**M5: [driver.rs:721](../../../../src/ssr/cohort/driver.rs#L721) — `keep_fp_control` emit-invariance untested.** *(reliability, High)* The `variable`-read-before-`apply_fp_control` ordering is the crux of the emit-vs-GT separation the diff exists to provide; asserted nowhere. A future edit moving the read below the no-call would silently reintroduce FP-control gating under BIC/freebayes. **Fix:** `keep_fp_control_cleans_genotypes_without_changing_bic_emit_decision` (BIC + freebayes).

**M6: [em.rs:548](../../../../src/ssr/cohort/em.rs#L548) — no thread-determinism test for the new float sums under the toggles.** *(reliability, Medium)* `attribute_clean_homs`'s `expected_slipped` sum and the `ln_mono_refit` reduction run inside a rayon `par_iter`; thread-invariance rests on fixed per-locus sum order (it is correct today) but no test asserts byte-identity across `--threads` with a toggle on. **Fix:** a determinism test under `null_from_homs=1`.

**M7: [driver.rs:255](../../../../src/ssr/cohort/driver.rs#L255) — `PVC_SSR_NULL_HOM_FRAC` accepts out-of-range/NaN/inf.** *(errors, High; convergent with extras)* `.parse::<f64>().ok()` lets `2.5`, `-1`, `nan`, `inf` through with no `[0,1]` check. A negative or NaN `hom_frac` makes the `modal_count < hom_frac*total` gate always-false → **every plant, including pure hets, is pooled as a clean-hom stutter donor**, silently corrupting the null (and, under `null_in_genotyper`, the genotypes). **Fix:** validate the domain and fail loud (mirror `InvalidEmitModel`) or `clamp(0.0,1.0)`.

**M8: [driver.rs:246,250,261,293](../../../../src/ssr/cohort/driver.rs#L246) — silent-fallback inconsistent with the fail-loud invariant.** *(errors, High)* The four booleans use `.is_ok_and(|v| v == "1")`, so `=true`/`=yes`/`=TRUE`/`=0` all silently mean off; `null_hom_min_depth` silently falls back on bad input — while `PVC_SSR_EMIT_MODEL` right above fails loud with the comment "the no-silent-fallback invariant." A sweep launched with `PVC_SSR_NULL_FROM_HOMS=true` silently runs the default path, mislabeled "feature on." **Fix:** treat set-but-unparseable as a hard error and/or `tracing::warn!` the ignored value; one error contract for the whole toggle group.

**M9: [em.rs:113,115,145,146](../../../../src/ssr/cohort/em.rs#L145) — `4`/`0.8` magic numbers; `0.8` silently diverges from the scorer's `0.75`.** *(defaults, High doc-gap / Medium intent)* Nothing records that these are the clean-hom donor thresholds borrowed from `silver_standard.py` (`MIN_DEPTH=4`, `HOM_FRAC=0.75`), and the code uses `0.8 ≠ 0.75` with no note of intent. Doc restates the literal in prose, so doc and value can drift. **Fix:** hoist to named consts with provenance comments; document the deliberate deviation, or correct to 0.75 (OQ2).

**M10: [em.rs:104,112,125](../../../../src/ssr/cohort/em.rs#L104) — three co-dependent bools are one enum.** *(smells High + idiomatic Minor, convergent)* `bic_null_refit`/`null_from_homs`/`null_in_genotyper` encode mutually-exclusive null strategies; exclusivity is hand-enforced by scattered `&& !other` guards + prose precedence. The 2³ truth table has meaningless silently-collapsed rows. **Fix:** one `NullStrategy` enum decoded once in the driver, `null_hom_*` kept as fields. *(Becomes unnecessary if OQ1 = delete — only `FromHoms` survives.)*

**M11: [em.rs:654-732](../../../../src/ssr/cohort/em.rs#L654), [driver.rs:246-261](../../../../src/ssr/cohort/driver.rs#L246) — two known-dead toggles kept with no cleanup plan.** *(smells, High)* `bic_null_refit` (~46-line block) and `null_in_genotyper` (genotyper branch + refit override) are empirically dead ends yet fully wired, with no `TODO`/removal condition (unlike the sibling `MARGINALIZED_PRIOR` knob). **Fix (recommended):** delete both; only `null_from_homs` + its two knobs + `keep_fp_control` survive, and M10/Mi1/Mi3 largely vanish. Depends on OQ1.

**M12: [driver.rs:262,277](../../../../src/ssr/cohort/driver.rs#L262) — `..dev_default()` re-introduces the silent-default trap.** *(refactor_safety, High mechanism / benign now)* Adding the six fields forced a compile error only at the `dev_default()` body; ~11 struct-update sites (8 EmCfg test literals, `bakeoff.rs:231`, 2 FpControlCfg test literals, and *both production sites*, which still end in `..dev_default()`) silently absorbed the new fields. Safe only because the defaults are no-ops — not because the compiler enforced it; the PROJECT_STATUS "flags every construction site" claim is only partially true. **Fix:** a comment at the two production sites recording that new fields default silently there (minimal), or correct the PROJECT_STATUS wording; full exhaustive-destructure is heavier and optional.

### Minor

- **Mi1: [em.rs:722](../../../../src/ssr/cohort/em.rs#L722) vs [em.rs:776](../../../../src/ssr/cohort/em.rs#L776) — duplicated best-monomorphic loop.** *(idiomatic + smells, convergent)* Verbatim copy of the `emission_evidence` reduction. Extract `best_monomorphic_ll(data_ll, genotypes, k)` as an iterator `fold(NEG_INFINITY, f64::max)`. *(Vanishes if OQ1 = delete.)*
- **Mi2: [em.rs:1177](../../../../src/ssr/cohort/em.rs#L1177) — `attribute_clean_homs` duplicates `attribute_locus`'s slip core.** *(smells, Medium)* Follow-up note; factor the `(delta,count,parent_level)→fit` accumulation if a third caller appears.
- **Mi3: [em.rs:681-732](../../../../src/ssr/cohort/em.rs#L681) — ~50-line inline probe block lengthens `run_locus_em_with`.** *(smells)* Hoist to a helper if retained; deleting idea A removes it. *(Depends on OQ1.)*
- **Mi4: [driver.rs:256-293](../../../../src/ssr/cohort/driver.rs#L256) — resolved toggles neither logged nor precedence-tested.** *(defaults)* Add a one-line dump per resolved toggle + a default-vs-override test.
- **Mi5: [driver.rs:251](../../../../src/ssr/cohort/driver.rs#L251) — `PVC_SSR_NULL_HOM_MINDEPTH=0` disables the depth gate.** *(errors)* `total < 0` never fires. Reject `0` alongside the M7 validation.
- **Mi6: [em.rs:654,655,663,664,703,711](../../../../src/ssr/cohort/em.rs#L654) — codename locals (`data_ll_b`, `fit_b`, `theta_b`, `m_b`, `fit0`, `m0`).** *(naming)* Rename to clean-hom / mono-refit domain names. *(A-block ones vanish if OQ1 = delete.)*
- **Mi7: [em.rs:664,711](../../../../src/ssr/cohort/em.rs#L664) — `m_b`/`m0` abbreviate the crate-standard `level_multiplier`.** *(naming)* Fold into Mi6.
- **Mi8: [em.rs:86,99](../../../../src/ssr/cohort/em.rs#L86) — `null_from_homs` vs `null_in_genotyper` named on inconsistent axes** (source vs consumer). *(naming, Medium)* Vary only a scope suffix. *(Moot if OQ1 = delete `null_in_genotyper`.)*

### Nits

- `cargo fmt` fails at `em.rs:1209` (in-scope, `attribute_clean_homs`) — mechanical `cargo fmt` (convergent across ~all categories).
- `em.rs:124` doc comment uses `[[wikilink]]` syntax that rustdoc renders as literal text.
- `bic_null_refit`/`null_from_homs` recompute runs regardless of `emit_model` → wasted `compute_data_ll` under `Heuristic`/`Freebayes` (benign; extras + module_structure).
- Four toggle doc-comments run 6–9 `///` lines each; shrinks with the dead-toggle deletion.
- `null_hom_frac` could be `null_hom_modal_frac` (naming, optional).

## 7. Out of scope observations

- Pre-existing `cargo fmt` diffs in `examples/ssr_psp_seqdump.rs` and `src/paralog/locus_score.rs`; pre-existing `cargo doc` unresolved-link errors (7) — separate cleanup.
- Pre-existing unvalidated numeric knobs (`sfs_theta`, `bic_margin`, `freebayes_min_qual`, `min_corroboration`) share the silent-`.parse().ok()` pattern; M7/M8 are about the *new* bounded-domain knobs, but a house-wide validation pass would be a reasonable follow-up.

## 8. Missing tests to add now

From the reliability agent (adapt to the crate's actual `CohortLocus::new`/`SampleEvidence`/`clean_params` test helpers; verify each fails on the mutation it targets):

- `attribute_clean_homs_bins_off_modal_reads_as_slips` — clean-hom donor with Δ=−1/+1/−2 reads; catches wrong modal pick or slip sign/magnitude binning.
- `attribute_clean_homs_excludes_a_plant_below_hom_frac` — 50/50 het; catches admitting hets as clean homs (the circularity the function avoids).
- `attribute_clean_homs_excludes_a_plant_below_min_depth` — `total < min_depth`; catches the depth-gate off-by-one (add a `total == min_depth` admits-boundary case).
- `attribute_clean_homs_returns_empty_fit_when_no_clean_homs` — all donors filtered; catches the empty-pool → frozen-prior degradation contract.
- `attribute_clean_homs_tie_break_is_deterministic` — two lengths tied at modal (use `hom_frac=0.5`); pins the `max_by_key`-returns-last choice guarding thread-invariance.
- `attribute_clean_homs_counts_far_slip_in_slipped_but_not_profile` — `|Δ| > MAX_SLIP`; locks the `slipped`-vs-`profile` divergence.
- `keep_fp_control_cleans_genotypes_without_changing_bic_emit_decision` — BIC-polymorphic locus whose sample `apply_fp_control` would no-call; assert `drop == false` + `EmittedVariable` while GT is no-called. Repeat for `Freebayes`.
- `null_from_homs_leaves_genotypes_unchanged_but_shifts_emission_evidence`.
- *(If OQ1 = keep)* `null_in_genotyper_with_empty_clean_hom_pool_matches_frozen_seed`; `bic_null_refit_never_lowers_ln_monomorphic`.

## 9. What's good

- **Off-path byte-identity is structural, not incidental:** `data_ll_emit = data_ll_b.as_deref().unwrap_or(&data_ll)` and the `(seed.theta0, 1.0)` / `refit_rounds` fallbacks reduce to the exact pre-diff bindings when off — verified by three independent agents.
- **The emit-vs-GT separation is correctly ordered:** `variable` is read before `apply_fp_control`, so `keep_fp_control` cleans genotypes without gating emission (the whole point of the toggle).
- **Thread-determinism preserved:** every new float sum is serial within a single-locus `par_iter` item, in fixed sample→observation order over integer slip counts.
- **Module boundaries respected:** no new `use` edges; `attribute_clean_homs` correctly sits beside `attribute_locus`; env parsing stays centralized in `driver::run`; the `EmCfg`/`FpControlCfg` split tracks the consuming module.
- **`bic_null_refit` is accurately scoped** — it mutates only `ln_monomorphic`, read only by the BIC arm; no freebayes leak.

## 10. Commands to re-verify

- `cargo fmt --check` (expect the `em.rs:1209` in-scope hit until fixed).
- `cargo clippy --lib -- -D warnings` (clean).
- `cargo test --lib ssr::cohort` (250 pass; will rise once §8 tests land).
- New: the §8 tests once added.

### Author response convention
Answer OQ1/OQ2 first, then address findings by ID (`fixed in <commit>` / `disputed because …` / `deferred to …` / `won't fix because …`).

Audit trail: `tmp/review_2026-07-09_ssr-emission-null-toggles/` (per-category files).
