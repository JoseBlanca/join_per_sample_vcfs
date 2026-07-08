# Code Review: rough_snp_het_caller
**Date:** 2026-07-07
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** Quality-aware ε̂ + strand-bias veto in the Stage-1 rough SNP het caller (commit `87cc19b`, branch `rough-snp-heuristics`)
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the PR diff of commit `87cc19b` (`git diff HEAD~1 HEAD`).
- **Reviewed against:** commit `87cc19b` on branch `rough-snp-heuristics`.
- **In-scope files:**
  - [src/sample_summary/het.rs](../../../../src/sample_summary/het.rs) — primary: per-site ε̂, strand veto, `AlleleGroupStats`/`SiteCounts` reshape, `HetClassifyParams`, constructor asserts, tests, two rejected-approach notes.
  - [src/sample_summary/mod.rs](../../../../src/sample_summary/mod.rs) — `DEFAULT_HET_STRAND_BIAS_Z`.
  - [src/pileup/per_sample/pileup_to_psp.rs](../../../../src/pileup/per_sample/pileup_to_psp.rs) — the feed (builds `AlleleGroupStats` from the record).
  - [src/pop_var_caller/cli.rs](../../../../src/pop_var_caller/cli.rs) — constructs `HetClassifyParams`; `PVC_HET_STRAND_Z` env override.
  - [examples/paralog_score_parity.rs](../../../../examples/paralog_score_parity.rs) — mirrors the feed.
  - [examples/het_baseline.rs](../../../../examples/het_baseline.rs) — new measurement example.
- **Deliberately out of scope:** the research report `doc/devel/reports/research/rough_snp_calling_heuristics_2026-07-07.md` (not code); pre-existing gate failures in untouched files (see §7).
- **Categories dispatched:** reliability, errors, naming, defaults, idiomatic, refactor_safety, module_structure, smells, tooling, extras. `unsafe_concurrency` skipped (no `unsafe`/`Arc`/`Mutex`/atomics/channels/`async`/threads in the diff).

## 2. Verdict

**Approve-with-changes.** The caller logic is correct and the change does what its commit message claims (verified by the `extras` diff-matches-intent pass: ε̂ only rises from the floor; the strand veto touches only het calls; the MAPQ/min-alt items are genuinely notes-only). No Blockers. The changes to make before merge are hardening (env-override validation, one unenforced public-API invariant), reproducibility (record the veto threshold in the `.psp`), gate hygiene (`cargo fmt`; unblock clippy so the new code is actually linted), and a set of regression tests for the new load-bearing guards.

## 3. Execution status

Run in the container per the dev workflow. Verbatim tails in `tmp/review_2026-07-07_rough-het/verification.txt`.

- `cargo fmt --check` → **exit 1.** In-scope violations introduced by this diff: `examples/het_baseline.rs:63`, `src/pileup/per_sample/pileup_to_psp.rs:25`. (Other hunks are pre-existing drift in untouched files.)
- `cargo clippy --all-targets --all-features -- -D warnings` → **exit 101**, aborting on a **pre-existing** `needless_lifetimes` at `src/vcf/writer.rs:384` (untouched). Consequence: the in-scope files were **never reached by clippy** (see M4).
- `cargo test --all-targets --all-features` → **1621 passed, 1 failed.** The one failure, `var_calling::posterior_engine::tests::larger_ref_pseudocount_cannot_increase_p_alt`, is **pre-existing and out of scope** (tracked in project memory). **All `sample_summary::het` tests pass.**
- `cargo doc --no-deps` → **exit 101**, on **pre-existing** broken intra-doc links in `src/var_calling/diversity.rs` (untouched).
- `cargo audit` → not run: `cargo-audit` not installed (deferred, no dependency changes in this diff).

Findings labeled "Needs verification": 0 (all locations were read; M5's exploit path is reachability-caveated, not unverified).

## 4. Open questions and assumptions

1. **Is `PVC_HET_STRAND_Z` intended to survive as a shipped knob, or is it session-only measurement scaffolding?** If shipped, it needs validation + provenance (M1, M2); if it is only for the sweep and should not ship, remove it and the whole question dissolves. Affects **M1, M2**.
2. **Should the reshaped `AlleleGroupStats`/`SiteCounts` stay `pub` with public fields, or move behind a validating constructor?** The `fwd ≤ obs` invariant (M5) and the non-exhaustive ALT construction (Mi1) both hinge on this. A `SiteCounts::from_record` constructor resolves M5, Mi1, and the naming concern Mi4 at once. Affects **M5, Mi1, Mi4**.

## 5. Top 3 priorities

1. **M1 — validate the `PVC_HET_STRAND_Z` env override at the CLI boundary.** Today a typo silently falls back to the default (a wrong sweep with no signal), and `0`/`-1`/`nan` sail through to the constructor `assert!` and panic inside a function that otherwise returns `Result`.
2. **M2 — record `strand_bias_z` in the `.psp` `SampleSummary`.** Three of the four het knobs (`min_depth`, `error_rate`, `lr_margin`) are already persisted in `HetCounts`; the veto threshold is not, so a `.psp` built under an override is metadata-indistinguishable from a default one — its het counts and F are unreproducible.
3. **M3 + M4 — make the gates green for this code.** `cargo fmt` the two touched files, and unblock clippy (the new ε̂/strand code has currently received **zero** clippy coverage because the gate aborts on a pre-existing error before reaching it).

## 6. Findings

### Major

**M1: [src/pop_var_caller/cli.rs:461](../../../../src/pop_var_caller/cli.rs#L461) — `PVC_HET_STRAND_Z` env override has no validation (silent fallback + panic-on-bad-value)**
**Categories:** errors, defaults (convergent). **Confidence:** High.
`strand_bias_z: std::env::var("PVC_HET_STRAND_Z").ok().and_then(|s| s.parse::<f64>().ok()).unwrap_or(DEFAULT_HET_STRAND_BIAS_Z)` has two silent failure modes:
1. A malformed/typo'd/non-UTF-8 value is swallowed by `.ok()` and silently replaced by the default — a threshold sweep is then quietly wrong with no diagnostic (violates the "recovery points emit a log" rule).
2. A parseable-but-invalid value (`0`, `-1`, or `nan` → `Ok(NaN)`) flows unchecked into `HetAccumulator::new`, whose `assert!` ([het.rs:217](../../../../src/sample_summary/het.rs#L217)) then panics — turning a runtime config typo into a panic inside `run_pileup`, which otherwise returns `Result<_, PileupCliError>`.
**Fix:** parse and validate at the CLI boundary, distinguishing "unset" from "invalid", and return a `PileupCliError` on a bad value rather than deferring to the constructor assert. (The constructor `assert!`s remain the right guard for `error_rate`/`lr_margin`, which are only ever fed compile-time constants — the problem is specifically the new *runtime* path reaching the same assert.)

**M2: [src/pop_var_caller/cli.rs:461](../../../../src/pop_var_caller/cli.rs#L461) / [src/sample_summary/het.rs:300](../../../../src/sample_summary/het.rs#L300) — the strand-veto threshold is not recorded in the `.psp` `SampleSummary`**
**Categories:** defaults, extras (convergent). **Confidence:** High.
`HetCounts` (via `finish()`) persists `min_depth`, `error_rate`, and `lr_margin`, but **not** `strand_bias_z`. Because the veto directly changes `n_het_sites` (and thus `obs_het` → F), two `.psp` files produced with different `PVC_HET_STRAND_Z` values are byte-for-byte indistinguishable in their metadata — the parameter that produced the het counts is unrecoverable. This is an inconsistent omission relative to the three knobs already persisted, and it undermines the pipeline's byte-identity/reproducibility discipline.
**Fix:** add `strand_bias_z` to `HetCounts` and thread it through `finish()`, beside the three existing knobs. (At minimum, emit a startup log when the override is set — but persisting it matches the existing pattern and is the real fix.)

**M3: [examples/het_baseline.rs:63](../../../../examples/het_baseline.rs#L63), [src/pileup/per_sample/pileup_to_psp.rs:25](../../../../src/pileup/per_sample/pileup_to_psp.rs#L25) — diff committed non-`cargo fmt`-clean**
**Categories:** tooling (+ noted by idiomatic, errors, smells, reliability, extras). **Confidence:** High (verified: `cargo fmt --check` exit 1).
Two in-scope hunks fail `rustfmt`: the `println!` header args in the example (must be one per line) and the multi-line `use` block in the feed (rustfmt collapses it to one line). `cargo fmt --check` is a merge-blocking gate; this diff turns it red on its own account.
**Fix:** `cargo fmt` scoped to the two touched files. Do **not** sweep the pre-existing drift in untouched files (`examples/ssr_psp_seqdump.rs`, `src/paralog/locus_score.rs`, `src/var_calling/paralog_filter/write_pass.rs`, `src/var_calling/vcf_writer.rs`) into this PR.

**M4: gate health — the in-scope code was never clippy-linted**
**Categories:** tooling. **Confidence:** High.
`cargo clippy --all-targets --all-features -- -D warnings` aborts at `src/vcf/writer.rs:384` (a pre-existing `needless_lifetimes`, untouched by this diff) before it ever compiles the changed files, so the new ε̂/strand-veto code and the new example have received **zero** clippy coverage. The root cause is out of scope, but the *consequence* — unlinted new code merging — is in scope.
**Fix:** apply the one-line lifetime elision at `src/vcf/writer.rs:384` so the gate reaches the diff (smallest unblock), or run clippy scoped to the changed paths, and confirm the new code is clean before merge.

**M5: [src/sample_summary/het.rs:332](../../../../src/sample_summary/het.rs#L332) — `strand_biased` trusts an unenforced `fwd ≤ obs` invariant on a public struct**
**Categories:** reliability, errors (convergent). **Confidence:** Medium (production-unreachable; reachable via the public API).
`let (rr, ar) = (nr - rf, na - af);` subtracts on `u64` under the comment "`fwd ≤ obs` by construction." The shipping pileup feed does maintain that invariant ([open_record.rs](../../../../src/pileup/walker/open_record.rs) accumulates `fwd` under an underflow guard). But `AlleleGroupStats` has independent `pub obs` / `pub fwd` fields, derives `Default`, and validates nothing, so a directly-constructed `AlleleGroupStats { obs: 15, fwd: 20, .. }` (examples, tests, any downstream) underflows: a debug panic, or in release a wrap to ~`u64::MAX` that passes the margin guard and yields a garbage `z` → **silently-wrong het/ambiguous classification.**
**Fix:** `debug_assert!(rf <= nr && af <= na, "fwd must not exceed obs")` before the subtraction, plus `saturating_sub` for the release path, plus the negative test in §8. Alternatively, a validating `SiteCounts::from_record` constructor with private fields (see Mi1) removes the class of bug.

**M6: [src/sample_summary/het.rs:206](../../../../src/sample_summary/het.rs#L206), [het.rs:217](../../../../src/sample_summary/het.rs#L217) — the new/changed constructor guards have no regression tests**
**Categories:** reliability. **Confidence:** High.
Two load-bearing constructor invariants are unverified: (a) the new `strand_bias_z` `assert!` (`NaN`/`≤ 0`) — if regressed, a `NaN` threshold makes `z.abs() > NaN` always `false`, silently disabling the veto and over-counting hets with no panic; (b) the `error_rate` upper bound was tightened from `< 1.0` to `< MAX_SITE_ERROR_RATE` (0.4) — a real behavior change whose only test still sets the *lower* bound (`0.0`). The tightened bound is the sole thing preventing `clamp(min=0.5, max=0.4)` from panicking inside the hot `observe_site`.
**Fix:** the `#[should_panic]` tests in §8 (`new_panics_on_strand_bias_z_{nan,zero}`, `error_rate_at_ceiling_panics`).

### Minor

**Mi1: [src/pileup/per_sample/pileup_to_psp.rs:97](../../../../src/pileup/per_sample/pileup_to_psp.rs#L97), [examples/paralog_score_parity.rs:331](../../../../examples/paralog_score_parity.rs#L331) — duplicated feed logic + non-exhaustive ALT-group construction**
**Categories:** smells, refactor_safety, idiomatic (convergent). **Confidence:** High.
The "build `SiteCounts` from a `PileupRecord`" block (REF group from `alleles[0]`, ALT pooled over `alleles[1..]`, then `observe_site`) is near-verbatim at both sites — and the parity example exists precisely to reproduce the production feed, so drift there fails silently. Separately, within each block the REF group is an *exhaustive* struct literal (compiler-checked) but the ALT group is `AlleleGroupStats::default()` + a field-wise `+=` loop: if a field is later added to `AlleleGroupStats`, the compiler flags the REF literal but stays silent on the ALT loop → the new field silently fails to aggregate on the ALT side, biasing the REF-vs-ALT veto.
**Fix:** one `SiteCounts::from_record(&PileupRecord) -> Option<SiteCounts>` (or `AlleleGroupStats::from_support` returning an exhaustive literal) used by both call sites — dedups the logic *and* makes the ALT aggregation exhaustively field-checked in a single move.

**Mi2: [src/sample_summary/het.rs:198](../../../../src/sample_summary/het.rs#L198) — stale `HetAccumulator::new` panic doc**
**Categories:** smells. **Confidence:** High.
The `/// # Panics` section still documents `error_rate` as "not finite in `(0, 1)`" (the code now asserts `< MAX_SITE_ERROR_RATE` = 0.4) and omits the new `strand_bias_z` panic entirely. The panic *message* was updated; the doc comment was not.
**Fix:** update the Panics section to the current `(0, MAX_SITE_ERROR_RATE)` bound and add the `strand_bias_z` clause.

**Mi3: [src/sample_summary/het.rs:561](../../../../src/sample_summary/het.rs#L561) — `strand_biased` is a `-> bool` predicate not named as one**
**Categories:** naming. **Confidence:** Medium.
The crate uses the `is_`/`has_` predicate convention; a `-> bool` named `strand_biased` reads as a state, not a question.
**Fix:** rename to `is_strand_biased`.

**Mi4: [src/sample_summary/het.rs:372](../../../../src/sample_summary/het.rs#L372) — `AlleleGroupStats.obs` renames the canonical `num_obs`**
**Categories:** naming. **Confidence:** Medium.
The record's `AlleleSupportStats.num_obs` ([pileup_record.rs:46](../../../../src/pileup_record.rs#L46)) is the canonical name for the identical fragment-count concept, and the feed even writes `obs: u64::from(s.num_obs)`. Same concept, two names.
**Fix:** rename `obs` → `num_obs` (folds naturally into the `from_record` constructor of Mi1).

**Mi5: [src/sample_summary/mod.rs:83](../../../../src/sample_summary/mod.rs#L83) — the `DEFAULT_HET_*` default family is split incoherently**
**Categories:** module_structure, naming (convergent). **Confidence:** Medium.
`DEFAULT_HET_MIN_DEPTH`/`_ERROR_RATE`/`_LR_MARGIN` are literals in the `sample_summary` façade, but the fourth default is defined as `het::DEFAULT_STRAND_BIAS_Z` and re-exported from the façade under a *second* name, `DEFAULT_HET_STRAND_BIAS_Z`. One value, two public spellings, no rule for where a default's source of truth lives.
**Fix:** pick one — move the `3.0` literal up beside its three siblings, or re-export the others through `het::` too. (Either is fine; the point is consistency.)

**Mi6: [src/sample_summary/het.rs:233](../../../../src/sample_summary/het.rs#L233) — `observe_site` NaN-freedom silently depends on `min_depth ≥ 1`**
**Categories:** errors. **Confidence:** Medium.
The "no NaN" guarantee relies on `n > 0` at the `log_error_sum / n` division; `n > 0` holds only because `min_depth ≥ 1` (default 4) makes the `n_total < min_depth` guard reject `n = 0`. The constructor asserts the other three params but not `min_depth`. Currently unreachable (not CLI-exposed), but the invariant is unenforced.
**Fix:** add `assert!(params.min_depth >= 1, ...)` to `HetAccumulator::new`, or a comment pinning the dependency.

**Mi7: [src/sample_summary/het.rs:370](../../../../src/sample_summary/het.rs#L370) — `SiteCounts` grew to 48 bytes but still derives `Copy`**
**Categories:** idiomatic. **Confidence:** Low.
`AlleleGroupStats` is 24 bytes, `SiteCounts` now 48 — above the ~16-byte `Copy`-discipline guideline. `Copy` is not load-bearing (`observe_site` takes by value/move, `strand_biased` borrows).
**Fix:** consider dropping `Copy` (keep `Clone`); low priority.

### Nits

- **[src/sample_summary/het.rs:344](../../../../src/sample_summary/het.rs#L344) — unreachable `variance <= 0.0` branch with a misleading comment** (reliability, extras). Given the preceding every-margin-≥-2 guard, `fwd_total, rev_total ≥ 2` and `fwd_total + rev_total = nr + na`, so `p_pooled ∈ (0,1)` strictly and `variance > 0` always. The branch is harmless defense-in-depth but its comment ("Both alleles fully on one strand") describes a state the guard already excludes. Either annotate it as unreachable-given-the-guard or drop it.
- **`f64::INFINITY`-as-disabled sentinel** for the veto threshold (idiomatic, smells) — correct and documented, but `Option<f64>` would make the disabled state unrepresentable-when-wrong. Deferred (ripple cost across params + construction sites).
- The `.parse().ok()` in M1 conflates "unset" and "garbage" — folded into M1's fix.

## 7. Out of scope observations

Pre-existing issues in untouched files, surfaced by the verification run. None introduced by this diff; each is a separate follow-up.

- `src/vcf/writer.rs:384` — `needless_lifetimes` fails `cargo clippy -D warnings` (a one-line elision). Blocks the clippy gate crate-wide; see M4 for the in-scope consequence.
- `src/var_calling/diversity.rs:6-7` — broken intra-doc links (`doc/devel/...` paths written as `[`...`]`) fail `cargo doc --no-deps`.
- `var_calling::posterior_engine::tests::larger_ref_pseudocount_cannot_increase_p_alt` — pre-existing test failure (tracked in project memory: EM convergence pseudocount leak).
- `cargo fmt` drift in `examples/ssr_psp_seqdump.rs`, `src/paralog/locus_score.rs`, `src/var_calling/paralog_filter/write_pass.rs`, `src/var_calling/vcf_writer.rs`.
- `cargo-audit` not installed (deferred item Mi12; no dependency changes here).

## 8. Missing tests to add now

Grouped by function. Bodies are in `tmp/review_2026-07-07_rough-het/reliability.md`.

**`HetAccumulator::new`** (guards M6):
- `new_panics_on_strand_bias_z_nan`, `new_panics_on_strand_bias_z_zero` — `strand_bias_z ∈ {NaN, 0.0}`; catches regression of the veto-disabling guard.
- `error_rate_at_ceiling_panics` — `error_rate = MAX_SITE_ERROR_RATE` (0.4) and `0.5`; catches the tightened bound regressing into a `clamp(min>max)` panic in `observe_site`.
- `new_panics_on_zero_min_depth` (Mi6) — `min_depth = 0`; pins the NaN-freedom precondition.

**`strand_biased` / `observe_site`** (guards M5 + the "hets only" contract):
- `strand_biased_rejects_fwd_exceeding_obs` — directly-constructed `fwd > obs`; catches the release-mode wrap.
- `strand_confined_hom_alt_is_not_vetoed` — confident hom-alt with strand-confined ALT stays hom-alt; catches the veto being hoisted out of the het branch.
- `strand_skew_on_non_het_does_not_touch_counts` — a lone-alt (hom-ref) strand-skewed site increments nothing; confirms the veto is not reached for non-hets.
- `infinite_log_error_sum_clamps_to_floor_no_nan` — `log_error_sum = -inf`; pins ε̂ clamp NaN-safety.

**`SiteCounts::log_error_sum`**:
- `log_error_sum_adds_ref_and_alt_groups` — direct call with distinct REF/ALT sums; catches a dropped/double-counted component.

## 9. What's good

- **Numeric robustness of ε̂ is real, not assumed** — the clamp handles `log_error_sum` of `0`, slightly-positive floating error, and `-inf` without producing `NaN`; `n > 0` is genuinely guaranteed by the min-depth early return ([het.rs:233-262](../../../../src/sample_summary/het.rs#L233)).
- **The veto is correctly scoped to het calls only** — `strand_biased` runs inside the `log_lr > lr_margin` branch, so hom-alt/ambiguous are provably untouched, matching the commit's "hets only" claim ([het.rs:287](../../../../src/sample_summary/het.rs#L287)).
- **Determinism is preserved** — the ALT-group aggregation is a fixed allele-index loop with no HashMap iteration or parallel float reduction, so the `.psp` het counts stay reproducible.
- **The two new tuning constants are exemplary** — `MAX_SITE_ERROR_RATE` and `DEFAULT_STRAND_BIAS_Z` carry their value choice, units, the `z ≈ 3 ⇒ p ≈ 0.003` derivation, research-report provenance, and the `INFINITY`-disables convention, at both the const and the field ([het.rs:317-326](../../../../src/sample_summary/het.rs#L317)).
- **The MAPQ removal is clean and the rejected-veto note is a legitimate *why*-comment** — `mapq_sum`/`mapq_sum_sq` remain live on `AlleleSupportStats` for the downstream paralog filter/merger/VCF; nothing dangles ([het.rs, note above `AlleleGroupStats`](../../../../src/sample_summary/het.rs)).

## 10. Commands to re-verify

Reviewer ran (in the container): `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --all-targets --all-features`, `cargo doc --no-deps`. Re-run after fixes.

New invocations the review introduces:
- After M3: `cargo fmt --check` should pass for the two in-scope files.
- After M4: `cargo clippy --all-targets --all-features -- -D warnings` should compile through to the in-scope files (fix the pre-existing `vcf/writer.rs:384` elision first) and report clean on them.
- After M6/§8: `cargo test --lib sample_summary::het` should show the new `#[should_panic]` and veto-contract tests passing.

### Author response convention

Address each finding by its identifier (M1–M6, Mi1–Mi7) with `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`, answering the §4 open questions first.
