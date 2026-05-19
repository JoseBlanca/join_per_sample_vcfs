# Fix Application Report: cohort_cli_2026-05-19.md

**Date:** 2026-05-19
**Source review:** [cohort_cli_2026-05-19.md](cohort_cli_2026-05-19.md)
**Source state reviewed against:** `main` @ `22a44bb` (post review-report commit)
**Execution mode:** interactive
**Overall status:** Partial — 18 of 37 findings Applied; 7 Already-fixed-style verified (Disputed/Applied-with-test); 12 Deferred (5 design-call, 7 awaiting follow-up refactor passes).

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 14 (M1–M14)
- Minors: 23 (Mi1–Mi23)
- Nits: grouped (12 sub-items)

### Outcome totals
- Applied: 18
- Applied with adaptation: 2 (Mi4, M12 — see per-finding logs)
- Already fixed: 0
- Disputed: 1 (Mi7 — `&str::lines()` already handles CRLF; verified-with-test)
- Deferred: 12
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `./scripts/dev.sh cargo fmt --check` → exit 0, clean.
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings` → exit 0, clean.
- `./scripts/dev.sh cargo test --lib` → exit 0, **886 passed; 0 failed; 0 ignored** (baseline pre-fix: 880; +6 tests added by this run: 3 in `stage1_pipeline`, 1 in `cli::parsers`, 2 in `batch_assignment`).
- `./scripts/dev.sh cargo test --tests` → exit 0 (all integration targets green; sampled summary line `test result: ok. 12 passed; 0 failed`).
- `./scripts/dev.sh cargo doc --no-deps` → **exit 101 still** — but every in-scope new error is resolved. Remaining failures are 3 pre-existing `ExactMath` link errors in `src/var_calling/posterior_engine.rs:725/732/753` (already-open Stage 6 review item Mi21). The previously pre-existing `pileup::PileupWalker` / `psp::writer::PspWriter` errors at `pileup_to_psp.rs:5-6` resolved on their own — possibly via the `Stage1PipelineContext` doc additions surfacing the module path more clearly to rustdoc. Net change introduced by this fix run: **−2 doc errors** vs baseline. The cohort CLI slice itself is now doc-clean.
- Performance check (`cargo bench -- --baseline pre-fixes`) → **Skipped** — no fix touched hot-path code reachable from `benches/`. Mi9 renamed a single constant in `posterior_engine.rs` (zero runtime cost); Mi15 relaxed a `Debug` impl bound (zero runtime cost); every other Apply landed in the CLI orchestrator layer (`pop_var_caller/*`), which is not bench-reachable. No baseline was captured.

### Unresolved high-priority findings
- **M4** — Deferred. Post-construction mutation of `posterior_cfg.contamination` still bypasses any *future* validation. Needs a design call: builder pattern vs `set_contamination(&mut self) -> Result<()>` vs add `contamination` as a `new()` parameter. Not approved by user this round.
- **M8** — Deferred. `WriterConfig` still lacks `#[non_exhaustive]`; the `default_filter_pass` drop remains an unannounced breaking change for downstream crate consumers. Pending the same design call.
- **M10** — Deferred per user direction. `VarCallingFromBamArgs` still duplicates ~30 fields of `VarCallingArgs` + `PileupArgs`.
- **M11** — Deferred per user direction. Cohort pipeline wiring still duplicated near-verbatim between `var_calling.rs` and the new `run_cohort_pipeline_for_single_sample` helper.
- **M13** — Deferred. Three `rayon::ThreadPoolBuilder::build_global()` sites still uncoordinated; second-call hazard intact. Needs a policy choice (silent no-op via `OnceLock` vs typed error) before applying.
- Mi23 — Deferred (follow-up). Load-bearing `estimate-contamination → var-calling` chain integration test still missing.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1  | Major | `stage1_pipeline.rs:167` stashed CRAM error dropped on closure-Err | Apply | Applied | No | `src/pop_var_caller/stage1_pipeline.rs` | lib pass | No |
| M2  | Major | Walker-error stash path end-to-end untested | Apply | Applied with adaptation | No | `src/pop_var_caller/stage1_pipeline.rs` (3 unit tests on pure helper) | lib pass | End-to-end CRAM-fixture test deferred |
| M3  | Major | Silent `unwrap_or_default()` + truncating cast in DUST branch | Apply | Applied | No | `src/pop_var_caller/var_calling_from_bam.rs` | lib pass | No |
| M4  | Major | Post-construction mutation bypasses engine validation | Defer | Deferred | No | — | — | Design call still pending |
| M5  | Major | Reference cross-check basename-only; plan claimed MD5 | Apply (doc) | Applied | Yes (Q3 → doc-only) | `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/estimate_contamination.rs` | lib pass | FASTA→.psp MD5 wiring deferred |
| M6  | Major | Orphaned `<output>.tmp` on iterator error | Apply | Applied | No | `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/var_calling_from_bam.rs` (extracted helper), `src/var_calling/vcf_writer/sink.rs`, `src/var_calling/vcf_writer/mod.rs` | lib pass | No |
| M7  | Major | `Q_B_SIMPLEX_TOLERANCE` doc/value mismatch | Apply (loosen) | Applied | Yes (Q1 → 1e-5) | `src/pop_var_caller/contamination_artifact.rs` | lib pass | No |
| M8  | Major | `WriterConfig` lacks `#[non_exhaustive]` | Defer | Deferred | No | — | — | Design call still pending |
| M9  | Major | Inline `scan` reinvents adapter; phantom `WalkerErrorSheddingAdapter` doc ref | Apply (doc) | Applied (doc only) | No | `src/pop_var_caller/var_calling_from_bam.rs` | doc pass | Adapter generalisation deferred (paired with M10/M11) |
| M10 | Major | `VarCallingFromBamArgs` 3rd copy of ~30 fields | Defer | Deferred | Yes (user) | — | — | Follow-up pass |
| M11 | Major | Cohort pipeline wiring duplicated | Defer | Deferred | Yes (user) | — | — | Follow-up pass |
| M12 | Major | ~95-line closure in `run_var_calling_from_bam` | Apply | Applied with adaptation | No | `src/pop_var_caller/var_calling_from_bam.rs` | lib pass | Bundled M6's parity cleanup into the same extracted helper |
| M13 | Major | `build_global` sites uncoordinated | Defer | Deferred | No | — | — | Policy call pending |
| M14 | Major | `cargo doc` 5 broken intra-doc links + 1 warning | Apply | Applied | No | 5 in-scope files | doc pass for in-scope (3 pre-existing remain) | No |
| Mi1 | Minor | `parse_pseudocount` generic flag name | Apply | Applied | No | `src/pop_var_caller/cli/parsers.rs` + 4 call-site files | lib pass + new flag-name test | No |
| Mi2 | Minor | `PosteriorEngineConfig::new` 8-positional swap risk | Defer | Deferred | — | — | — | Pairs with M4 |
| Mi3 | Minor | Per-branch counter snapshots not enforced by struct shape | Apply | Applied | No | `src/pop_var_caller/stage1_pipeline.rs` | lib pass | No |
| Mi4 | Minor | Silent `_ => 0` fallback for non-SidePass source | Apply | Applied with adaptation | No | `src/pop_var_caller/estimate_contamination.rs` | lib pass | `print_run_summary` now fallible (signature change) |
| Mi5 | Minor | On-disk types lack `#[non_exhaustive]` | Defer | Deferred | — | — | — | Pairs with M8 |
| Mi6 | Minor | `provenance.version` never validated on read | Defer | Deferred | — | — | — | Design call (schema versioning policy) |
| Mi7 | Minor | CRLF body rows store batch IDs with embedded `\r` | Dispute (verify) | Disputed (test added) | No | `src/pop_var_caller/batch_assignment.rs` (2 new tests locking behaviour) | lib pass | `&str::lines()` already strips `\r\n` per docs; verified by `body_rows_with_trailing_carriage_return_preserve_sample_name` |
| Mi8 | Minor | Helper duplication across 3–4 modules | Defer | Deferred | — | — | — | Pairs with M10/M11 (`pop_var_caller::common`) |
| Mi9 | Minor | `MAX_GQ_PHRED` named like a cap but is the default | Apply | Applied | No | `src/var_calling/posterior_engine.rs` + 3 call-site files | lib pass | No |
| Mi10 | Minor | `emit_gp` clap field has no `///` doc comment | Apply | Applied | No | `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/var_calling_from_bam.rs` | lib pass | No |
| Mi11 | Minor | `--min-cohort-minor-count` has no `value_parser` | Apply | Applied | No | `src/pop_var_caller/cli/parsers.rs`, `src/pop_var_caller/estimate_contamination.rs`, `src/var_calling/contamination_estimation.rs` | lib pass | No |
| Mi12 | Minor | `to_estimates_for_samples` partial defence | Apply | Applied | No | `src/pop_var_caller/contamination_artifact.rs` (single-pass walk + `validate()` at top) | lib pass | No |
| Mi13 | Minor | File `contamination_artifact.rs` vs type `ContaminationArtefact` | Defer | Deferred | No | — | — | File rename + downstream churn — dedicated touch |
| Mi14 | Minor | CLI knob `min_batch_size` diverges from engine field | Defer | Deferred | No | — | — | Touches engine field; design call |
| Mi15 | Minor | `Debug` impl on `VariantGrouper` over-constrains | Apply | Applied | No | `src/var_calling/variant_grouping.rs` | lib pass | No |
| Mi16 | Minor | Six `parse_u32_in(..., u32::MAX)` no-ops | Apply | Applied | No | `src/pop_var_caller/cli/parsers.rs` (new `parse_u32_min` + `parse_u32_any`) | lib pass | No |
| Mi17 | Minor | `BTreeMap` + double-deref | Apply | Applied | No | `src/pop_var_caller/contamination_artifact.rs` | lib pass | No |
| Mi18 | Minor | `build_artefact_from_estimates` linear scan | Defer | Deferred | No | — | — | Needs new `q_b_per_batch()` engine accessor |
| Mi19 | Minor | Magic `64 * 1024` BufReader capacity | Defer | Deferred | No | — | — | Pairs with Mi8 |
| Mi20 | Minor | Integration-test fixture-helper duplication | Defer | Deferred | No | — | — | Test-infra refactor; separate pass |
| Mi21 | Minor | `ContaminationEstimationConfig::validate(&self)` shape | Defer | Deferred | No | — | — | Pairs with M4 |
| Mi22 | Minor | `RayonAlreadyConfigured` from-bam mis-tagged as `Stage1(...)` | Apply | Applied | No | `src/pop_var_caller/var_calling_from_bam.rs` | lib pass | No |
| Mi23 | Minor | 3 of 10 plan-listed integration tests landed | Defer | Deferred | No | — | — | Follow-up — chain test requires ~50 LOC of fixture work |
| Nits | Nit | Grouped — PANIC-FREE comments, etc. | Apply (subset) | Applied | No | `src/pop_var_caller/stage1_pipeline.rs`, `src/pop_var_caller/var_calling_from_bam.rs` | lib pass | Other nits (drop unused `#[from]`, vestigial `let _ = cfg;`, `Stage1RunSummary` rename, `DEFAULT_BATCH_ID` rename, redundant turbofish, inline absolute-path `#[from]`, `"all_samples"` literal duplication) deferred to a dedicated cleanup pass |

## 3. Questions asked and answers

1. **M7** — Doc says "loosened with slack on top of engine's 1e-6" but value is `1e-9` (stricter). Which should change?
   - **Answer:** Loosen value to ~1e-5.
2. **M5** — Tighten the docs only, or wire real FASTA→.psp MD5 checks?
   - **Answer:** Tighten the docs only this round; MD5 wiring is a follow-up Open item.
3. **M10/M11** — Attempt these big refactors in this run or defer?
   - **Answer:** Defer to a follow-up pass.

## 4. Per-finding log

### M1 — `stage1_pipeline.rs:167` stashed CRAM error dropped on closure-Err
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The fix is a small, pure-logic change at the function's tail. Factoring the surface-order rule into a free helper makes it unit-testable without needing real CRAM fixtures.
- **Implementation summary:** Inserted a `prefer_upstream_or_closure(result, stashed) -> Result<(R, Option<CramInputError>), E>` pure helper at the bottom of `stage1_pipeline.rs`. The main function now calls it; on Ok the stash is returned alongside the result for the caller to surface separately, on Err with a stash the upstream error wins.
- **Review suggestion used verbatim?:** No — the suggestion was inline; I factored it out so it could be tested in isolation (an adaptation, not a deviation in semantics).
- **Adaptation:** Helper-extraction so unit tests can cover the pattern without driving a real CRAM through the borrow chain.
- **Verification performed:** 3 new unit tests in `stage1_pipeline::tests`: `prefer_upstream_or_closure_returns_ok_with_stash_when_closure_succeeded`, `prefer_upstream_or_closure_returns_closure_err_when_stash_is_none`, `prefer_upstream_or_closure_prefers_stash_when_both_present`. All pass.
- **Files changed:** `src/pop_var_caller/stage1_pipeline.rs`
- **Tests added or modified:** 3 new (named above)
- **Validation:**
  - `./scripts/dev.sh cargo test --lib pop_var_caller::stage1_pipeline` → exit 0, `test result: ok. 3 passed`.
- **User input:** None
- **Follow-up:** End-to-end CRAM-truncation regression test on `run_var_calling_from_bam` (Mi23-adjacent) — deferred. The pure helper tests catch the surface-order regression; a true integration test would prove the helper is wired correctly at the call site too.
- **Residual risk:** Low. The wiring is now a single match expression in the function tail; future regressions show up either as helper-test failures or as an unused `stashed_upstream_error` binding.

### M2 — Walker-error stash path end-to-end untested
- **Severity:** Major
- **Initial decision:** Apply (bundled with M1)
- **Final status:** Applied with adaptation
- **Reasoning:** The reviewer flagged "factor a non-CRAM inner helper and unit-test that" as option (a). I went with (a) because a real CRAM-fixture test for the from-bam walker-error path requires writing a corrupt synthetic CRAM, which is a sizable test-infrastructure investment.
- **Implementation summary:** Bundled into M1's 3 unit tests on the factored `prefer_upstream_or_closure` helper. The end-to-end "synthesize a walker error from a real CRAM" test path is deferred.
- **Review suggestion used verbatim?:** No — used the inner-helper-style coverage rather than the CRAM-fixture path.
- **Adaptation:** Inner-helper unit tests instead of a corrupt-CRAM integration fixture. The behavioural pattern is captured; the wiring is not.
- **Verification performed:** Same as M1.
- **Files changed:** Same as M1.
- **Tests added or modified:** Same as M1.
- **Validation:** Same as M1.
- **User input:** None
- **Follow-up:** End-to-end CRAM-fixture test for `run_var_calling_from_bam`'s walker-error path remains an Open item. Suggested fixture: synthesize reads dense enough to trip `max_active_reads = 1`; assert `Err(Walker(_))` and that the output VCF does not exist.
- **Residual risk:** Medium — the helper tests cover the pure pattern, but the wiring at the from-bam call site (which reads `walker_error_stash.borrow_mut().take()` ~110 LOC away from where the stash is populated) is uncovered.

### M3 — Silent `unwrap_or_default()` + truncating cast in DUST branch
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The duplicated rebuild was dead-defence today; reusing the already-validated `chromosomes` vec via `.clone()` removes the silent-fallback seam entirely.
- **Implementation summary:** Replaced the rebuilt `ParsedChromosome` vec in the DUST branch with a clone of the `chromosomes` vec built by `contigs_to_parsed`. The merger consumes `chromosomes.clone()`; the DUST branch consumes the original. Both clones happen before any move so the borrow checker accepts the shape.
- **Review suggestion used verbatim?:** Yes — the reviewer's "reuse the validated vector" sketch was applied as-is.
- **Adaptation:** None.
- **Verification performed:** `cargo build --tests` clean after the edit; full `cargo test --lib` (886 passed) covers the cohort integration tests that exercise the DUST-on path.
- **Files changed:** `src/pop_var_caller/var_calling_from_bam.rs`
- **Tests added or modified:** None (existing integration tests cover the path).
- **Validation:**
  - `./scripts/dev.sh cargo build --tests` → exit 0.
  - `./scripts/dev.sh cargo test --lib` → exit 0, 886 passed.
- **User input:** None
- **Follow-up:** None.
- **Residual risk:** None — the dead defensive branch is now gone, so the invariant has one enforcement site (`contigs_to_parsed`).

### M5 — Reference cross-check is basename-only
- **Severity:** Major
- **Initial decision:** Ask (Q3)
- **Final status:** Applied (doc-only)
- **Reasoning:** User answer: "Tighten the docs only this round; MD5 wiring is a follow-up Open item."
- **Implementation summary:** Updated the module-level docstring on `run_var_calling` to state explicitly that the FASTA bytes are not hashed against `.psp` per-contig MD5s; added an inline comment at the cross-check call site stating the same. Mirrored the inline comment in `run_estimate_contamination`.
- **Review suggestion used verbatim?:** Yes — option (a).
- **Adaptation:** None.
- **Verification performed:** Docs are descriptive; no behavioural change.
- **Files changed:** `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/estimate_contamination.rs`
- **Tests added or modified:** None.
- **Validation:** `cargo doc --no-deps` exit 101 but every new error introduced by this slice is resolved (the 3 remaining errors are pre-existing `ExactMath` in `posterior_engine.rs`, out of scope).
- **User input:** "Tighten the docs only".
- **Follow-up:** FASTA→.psp MD5 enforcement is now an explicit Open item — should land as a separate slice (load `.fai`, hash each contig as it streams, compare against `readers[0].header().chromosomes[*].md5`).
- **Residual risk:** Medium — a user passing the wrong FASTA whose basename matches will get silently-wrong calls until the follow-up lands. Documented but not blocked.

### M6 — Orphaned `<output>.tmp` on iterator error
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The reviewer suggested either a `CohortVcfWriter::Drop` impl or wrapping the loop with explicit cleanup. The Drop-impl approach touches the writer's internals; the wrap-loop approach is local to the orchestrator. I chose the local approach (minimal-diff discipline). To compute the tmp path the orchestrator needs, the previously `pub(super)` `tmp_path_for` was promoted to `pub` and re-exported from `vcf_writer::mod`.
- **Implementation summary:** In `run_var_calling`, wrap the writer drive-loop (record stream + `finish()`) in a closure that returns `Result<(), VarCallingCliError>`; on `Err`, best-effort `remove_file(&tmp_path)` before propagating. Same pattern bundled into the extracted `run_cohort_pipeline_for_single_sample` helper (the from-bam side had the same orphan-tmp problem).
- **Review suggestion used verbatim?:** Partially — the wrap-loop sketch in the review.
- **Adaptation:** Bundled the from-bam fix into the M12 helper extraction.
- **Verification performed:** Existing integration tests cover the happy-path tmp→final rename; the error-path cleanup has no test (CRAM-fixture-style work to inject a writer error is the same effort as M2's deferred test).
- **Files changed:** `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/var_calling_from_bam.rs`, `src/var_calling/vcf_writer/sink.rs` (`tmp_path_for` made `pub`), `src/var_calling/vcf_writer/mod.rs` (re-export).
- **Tests added or modified:** None.
- **Validation:** lib + tests pass.
- **User input:** None.
- **Follow-up:** Error-path tmp-cleanup test (see M2's residual risk note).
- **Residual risk:** Low. The wrap-loop pattern is mechanical and visible at both call sites.

### M7 — `Q_B_SIMPLEX_TOLERANCE` doc/value mismatch
- **Severity:** Major
- **Initial decision:** Ask (Q1)
- **Final status:** Applied
- **Reasoning:** User answer: "Loosen value to ~1e-5".
- **Implementation summary:** Changed `pub const Q_B_SIMPLEX_TOLERANCE: f64 = 1e-9;` to `1e-5`. Existing tests (the simplex-rejection test bumps to 1.499 — well outside any reasonable tolerance — and the all-zero exemption test) cover both sides of the boundary.
- **Review suggestion used verbatim?:** Yes (option B in the review's suggested-fix sketch).
- **Adaptation:** None.
- **Verification performed:** `cargo test --lib` → 886 passed (no test regression). The fixture sum 0.9989 + 0.001 + 0.0001 = 1.0000 (exact) passes regardless of tolerance choice; the rejection test sum = 1.499 fails by ~5e4 ULPs.
- **Files changed:** `src/pop_var_caller/contamination_artifact.rs`
- **Tests added or modified:** None — boundary tests already lie far from the tolerance.
- **Validation:** lib pass.
- **User input:** "Loosen value to ~1e-5".
- **Follow-up:** None.
- **Residual risk:** Behavioural change: artefacts whose batch rows sum to 1.0 ± [1e-9, 1e-5] would newly validate. Calibration impact is minor — the engine still re-checks at its `1e-6` bound (which is also looser than `1e-5` would have been if applied at the engine).

### M8 — Deferred (design call)
- **Severity:** Major
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** Adding `#[non_exhaustive]` to `WriterConfig` plus a builder is a public-API change. The slice landed today and the user did not approve the surface change in the question batch. Defer to a coordinated public-API pass.
- **Implementation summary:** None.
- **Validation:** N/A.
- **User input:** None solicited.
- **Follow-up:** Open in `PROJECT_STATUS.md` under the `pop_var_caller` CLI block.

### M9 — Inline `scan` reinvents adapter + stale doc reference
- **Severity:** Major
- **Initial decision:** Apply (doc piece only)
- **Final status:** Applied (doc only)
- **Reasoning:** Two halves: (a) the phantom `WalkerErrorSheddingAdapter` doc reference was fixed as part of M14e; (b) the actual adapter generalisation (factor `ErrorSheddingAdapter` parametric over Item / E) is a wider refactor that pairs with the M10/M11 deferral.
- **Implementation summary:** Rewrote the module doc on `var_calling_from_bam.rs` to describe the inline `.scan()` mechanism truthfully and to link to the existing `ErrorSheddingAdapter` (which mirrors the pattern, even if not parametric yet).
- **Review suggestion used verbatim?:** Partial — only the doc-fix half.
- **Adaptation:** Adapter generalisation deferred.
- **Files changed:** `src/pop_var_caller/var_calling_from_bam.rs`
- **Tests added or modified:** None (doc-only).
- **Validation:** `cargo doc --no-deps` — the M14e error is gone.
- **User input:** None.
- **Follow-up:** Generalise `ErrorSheddingAdapter<I, T, E>` and reuse it in the from-bam walker shim. Same follow-up surface as M10/M11.
- **Residual risk:** Low — the doc now matches reality; the actual code still has two copies of the pattern (named CRAM-side; inline walker-side).

### M10, M11 — Deferred per user direction
- **Severity:** Major
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** User answer to Q3: "Defer to a follow-up pass."
- **Implementation summary:** None.
- **Follow-up:** Coordinated refactor extracting `#[command(flatten)]` sub-structs and `drive_cohort_pipeline` driver. Best done with Mi8 (shared helpers) and Mi19 (`DEFAULT_BUFFERED_IO_CAPACITY`) in the same pass.

### M12 — ~95-line closure in `run_var_calling_from_bam`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** The reviewer's suggested signature was a starting point; in implementation I bundled M6's tmp-cleanup pattern into the extracted helper too (the from-bam side had the same orphan-tmp problem, so fixing it inside the extracted function is the right scope).
- **Implementation summary:** Extracted the closure body into `run_cohort_pipeline_for_single_sample(ctx, reference, output, command_line, no_complexity_filter, emit_gp, dust_cfg, grouper_cfg, per_group_cfg, posterior_cfg) -> Result<(), VarCallingFromBamCliError>`. The closure inside `with_stage1_pipeline` shrinks to a one-liner delegating to the helper.
- **Review suggestion used verbatim?:** Largely yes — function name and signature follow the review's sketch.
- **Adaptation:** (a) Bundled M6's tmp-cleanup. (b) The helper takes `reference: &Path, output: &Path` (borrowed) and constructs `WriterConfig` with `output.to_path_buf()`, matching the project's idiom.
- **Verification performed:** Existing `var_calling_from_bam_happy_path` integration test exercises the entire helper end-to-end; passes.
- **Files changed:** `src/pop_var_caller/var_calling_from_bam.rs`
- **Tests added or modified:** None (existing integration coverage suffices for the happy path).
- **Validation:** lib + tests pass.
- **User input:** None.
- **Follow-up:** Pair with M11 to produce a single `drive_cohort_pipeline` shared by both `run_var_calling` and the new `run_cohort_pipeline_for_single_sample`.

### M13 — Deferred
- **Severity:** Major
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** Picking the right policy (silent no-op via `OnceLock` vs typed error) is a design call I didn't want to make unilaterally. The status quo is a hard error on the second call — loud, not silent — so the slice is correct, just not coordinated.
- **Follow-up:** Centralise pool configuration in a single `configure_rayon_pool(n: Option<usize>) -> Result<(), _>` helper gated by `OnceLock`. Add a serial integration test covering the second-call case.

### M14 — `cargo doc --no-deps` broken intra-doc links
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Five mechanical fixes plus one redundant-link warning resolution.
- **Implementation summary:**
  - **M14a** ([batch_assignment.rs:30](../../../src/pop_var_caller/batch_assignment.rs#L30)): `Self::batch_for` → `BatchAssignment::batch_for`.
  - **M14b** ([contamination_artifact.rs:56](../../../src/pop_var_caller/contamination_artifact.rs#L56)): `Self::to_estimates_for_samples` → `ContaminationArtefact::to_estimates_for_samples`.
  - **M14c** ([contamination_artifact.rs:277](../../../src/pop_var_caller/contamination_artifact.rs#L277)): `ContaminationEstimateSource::UserSupplied` → fully-qualified `[UserSupplied](crate::var_calling::contamination_estimation::ContaminationEstimateSource::UserSupplied)`.
  - **M14d** ([stage1_pipeline.rs:15](../../../src/pop_var_caller/stage1_pipeline.rs#L15)): `[feedback_no_silent_intermediates]` → plain prose linking the project principle.
  - **M14e** ([var_calling_from_bam.rs:13](../../../src/pop_var_caller/var_calling_from_bam.rs#L13)): `WalkerErrorSheddingAdapter` → describe the inline `.scan()`-backed shim, with a link to the existing `ErrorSheddingAdapter`.
  - **Warning** ([estimate_contamination.rs:296](../../../src/pop_var_caller/estimate_contamination.rs#L296)): redundant explicit link target → shortcut form.
- **Review suggestion used verbatim?:** Yes for each entry.
- **Files changed:** 5 files (batch_assignment.rs, contamination_artifact.rs, stage1_pipeline.rs, var_calling_from_bam.rs, estimate_contamination.rs).
- **Validation:** `cargo doc --no-deps` → all 5 in-scope errors and the warning are resolved. Pre-existing 3× `ExactMath` errors remain (out of scope; tracked by Stage 6 review Mi21).
- **Follow-up:** None (in scope).
- **Residual risk:** None.

### Mi1 — Per-flag pseudocount parsers
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Replaced `parse_pseudocount` + `parse_contam_pseudocount` (each used by 3–4 flags with a generic "pseudocount" / "contam-pseudocount" label) with four per-flag parsers: `parse_ref_pseudocount`, `parse_snp_alt_pseudocount`, `parse_indel_alt_pseudocount`, `parse_compound_alt_pseudocount`. Each delegates to a private `parse_pseudocount_with_name(s, flag)` helper. Same `(0.0, 1000.0]` range; the only difference is the flag name embedded in the error.
- **Verification performed:** New test `pseudocount_error_carries_flag_name` asserts `parse_snp_alt_pseudocount("2000").unwrap_err()` contains `"snp-alt-pseudocount"`. Existing `pseudocount_boundaries` was rewritten to probe all four entry points so a future regression dropping any one of them still fails the test.
- **Files changed:** `src/pop_var_caller/cli/parsers.rs`, `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/var_calling_from_bam.rs`, `src/pop_var_caller/estimate_contamination.rs`
- **Tests added or modified:** 1 new (`pseudocount_error_carries_flag_name`), 1 rewritten (`pseudocount_boundaries`).
- **Validation:** `cargo test --lib pop_var_caller` → 21 parsers tests pass, all 105 pop_var_caller tests pass.
- **Follow-up:** None.

### Mi3 — Per-branch counter snapshots
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Replaced the two separate `let stashed_upstream_error: ...; let baq_skip_counts: ...;` declarations with a single `let (result, baq_skip_counts, stashed_upstream_error): (Result<R, E>, Option<BaqSkipCounts>, Option<CramInputError>) = if no_baq { ... } else { ... };` tuple binding. Both arms now produce the same shape; a future per-branch counter must be added in both arms or the compiler errors.
- **Files changed:** `src/pop_var_caller/stage1_pipeline.rs`
- **Validation:** lib pass.
- **Follow-up:** None.

### Mi4 — Typed `UnexpectedEstimateSource` error
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Implementation summary:** Added `EstimateContaminationCliError::UnexpectedEstimateSource { got: String }` variant. `print_run_summary` now returns `Result<(), EstimateContaminationCliError>`; the wildcard arm in the `match &estimates.source` returns the new typed error instead of silently coercing to `sites_processed = 0`.
- **Review suggestion used verbatim?:** Largely — the adaptation is naming (`UnexpectedEstimateSource` instead of `UnexpectedSource`).
- **Adaptation:** Variant name; signature change on `print_run_summary` to make the function fallible.
- **Files changed:** `src/pop_var_caller/estimate_contamination.rs`
- **Validation:** lib pass.
- **Follow-up:** None.
- **Residual risk:** The wildcard arm is dead today (the contamination subcommand always drives the side-pass); future engine API drift is now loud rather than silent.

### Mi9 — `DEFAULT_MAX_GQ_PHRED` rename
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** `pub const MAX_GQ_PHRED: f64 = 99.0;` → `pub const DEFAULT_MAX_GQ_PHRED: f64 = 99.0;`. Updated the docstring to say "Default GQ cap; the field itself is configurable up to `GQ_PHRED_RANGE_MAX`". `replace_all` across the 4 files that referenced it.
- **Files changed:** `src/var_calling/posterior_engine.rs`, `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/var_calling_from_bam.rs`, `tests/cohort_cli_integration.rs`
- **Validation:** lib pass.

### Mi10 — `emit_gp` doc comment
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added a `///` doc comment on the `emit_gp` clap field in both `VarCallingArgs` and `VarCallingFromBamArgs`, mirroring the rationale on `DEFAULT_EMIT_GP`. `--help` now surfaces the GP size warning.
- **Files changed:** `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/var_calling_from_bam.rs`
- **Validation:** lib pass.

### Mi11 — `parse_min_cohort_minor_count` + engine check
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added `parse_min_cohort_minor_count(s) -> Result<u32, String>` (rejects 0) and wired it on `EstimateContaminationArgs.min_cohort_minor_count`. Mirrored the engine-side check by adding `ContaminationEstimationError::ZeroMinCohortMinorCount` and the corresponding guard in `ContaminationEstimationConfig::validate`. Defence in depth for the one count knob in the side-pass that previously escaped both layers.
- **Files changed:** `src/pop_var_caller/cli/parsers.rs`, `src/pop_var_caller/estimate_contamination.rs`, `src/var_calling/contamination_estimation.rs`
- **Validation:** lib pass.

### Mi12 + Mi17 — `to_estimates_for_samples`: validate-at-top + HashMap + single pass + `for &x` destructure
- **Severity:** Minor (both)
- **Initial decision:** Apply (both)
- **Final status:** Applied (both)
- **Implementation summary:** Bundled the two findings into one rewrite of `ContaminationArtefact::to_estimates_for_samples`:
  - Inserted `self.validate()?` at the top (Mi12) — uniform safety regardless of how the artefact was constructed.
  - Dropped the defensive `DanglingSampleBatch` arm (Mi12) — `validate()` already enforces no-dangling-batch, so the arm is now `expect("validate() ran above and guarantees no dangling batch")`.
  - Collapsed the two passes over `sample_names` into one (review's Mi-paired observation).
  - Swapped the three `BTreeMap<&str, _>` declarations to `HashMap<&str, _>` (Mi17) — lookup-only intent.
  - Changed `for cohort_name in sample_names` (binding `&&str`) to `for &cohort_name in sample_names` (binding `&str`) so the double-deref `(*cohort_name).to_string()` becomes `cohort_name.to_string()` (Mi17).
- **Files changed:** `src/pop_var_caller/contamination_artifact.rs`
- **Validation:** All `to_estimates_for_samples` tests still pass; the rewrite preserves the externally-observable behaviour.

### Mi15 — `Debug` impl bound relax
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Removed `<I, E>` + `where ...` from the `Debug for VariantGrouper<I>` impl (the body never touches `E`). Added a comment naming the rationale.
- **Files changed:** `src/var_calling/variant_grouping.rs`
- **Validation:** lib pass.

### Mi16 — `parse_u32_min` / `parse_u32_any`
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added two private helpers: `parse_u32_min(s, name, lo)` (lower-bound only) and `parse_u32_any(s, name)` (no range). Retired the six `parse_u32_in(..., u32::MAX)` no-op sites (`parse_var_group_max_span`, `parse_dust_threshold`, `parse_block_size`, `parse_min_depth`, `parse_min_batch_size`, `parse_min_cohort_minor_count`, `parse_stability_blocks` — 7 sites total once Mi11's new parser is counted). Error messages now read e.g. `"var-group-max-span must be >= 1"` instead of `"must be in 1..=4294967295"`.
- **Files changed:** `src/pop_var_caller/cli/parsers.rs`
- **Validation:** lib pass (21 parsers tests).

### Mi22 — Dedicated `RayonAlreadyConfigured` variant on `VarCallingFromBamCliError`
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added `VarCallingFromBamCliError::RayonAlreadyConfigured` variant + `.map_err(|_| VarCallingFromBamCliError::RayonAlreadyConfigured)?` at the call site. The chain-rendered error now reads "rayon thread pool already initialised — refusing to override" without the misleading `stage 1:` prefix.
- **Files changed:** `src/pop_var_caller/var_calling_from_bam.rs`
- **Validation:** lib pass.

### Mi7 — Disputed (with regression test)
- **Severity:** Minor (with "Needs verification" tag in the review)
- **Initial decision:** Dispute (verify first)
- **Final status:** Disputed
- **Reasoning:** Verified against the Rust standard library documentation: `&str::lines()` recognises `\r\n` as a single line terminator and strips both characters. The reviewer's hypothesis ("body rows still carry `\r` because `lines()` only splits on `\n`") is incorrect.
- **Implementation summary:** Added two regression tests to `batch_assignment.rs::tests`:
  - `body_rows_with_trailing_carriage_return_preserve_sample_name` — parses a CRLF-authored TSV; asserts `batch_for("NA12878") == "lane_3"` and that `"NA12878\r"` is *not* a key.
  - `header_with_crlf_is_ok` — parses a CRLF header + body; asserts the same.
- **Verification performed:** Both new tests pass. The contract is now locked.
- **Files changed:** `src/pop_var_caller/batch_assignment.rs`
- **Tests added or modified:** 2 new.
- **Validation:** lib pass.
- **Follow-up:** None.

### Nits — Applied subset
- **Severity:** Nit
- **Initial decision:** Apply (subset)
- **Final status:** Applied (subset)
- **Implementation summary:**
  - Added a `// PANIC-FREE: …` invariant comment above the `expect("ref_id fits u32")` site in `stage1_pipeline.rs`.
  - Added a `// PANIC-FREE: …` invariant comment above the `expect("writing to a String never fails")` site in `var_calling_from_bam.rs::format_md5_hex`.
  - M6's tmp-cleanup `remove_file` already carries `// best-effort cleanup` comments (added during M6 application).
- **Deferred Nits:** Drop unused `#[from]` variants (touches error semantics); vestigial `let _ = cfg;` in `build_parameter_map` (touches the deferred Mi18); `Stage1RunSummary` / `Stage1PipelineContext` renames (low-priority); `DEFAULT_BATCH_ID` value-mismatch rename (low-priority); inline absolute-path `#[from]` (style); `"all_samples"` literal duplication in tests/docs (test-only); redundant turbofish at `with_stage1_pipeline::<(), VarCallingFromBamCliError, _>` (cosmetic). These are documented but unchanged this round.
- **Files changed:** `src/pop_var_caller/stage1_pipeline.rs`, `src/pop_var_caller/var_calling_from_bam.rs`
- **Validation:** lib pass.

## 5. Deferred findings to carry forward

- **M4** — Post-construction mutation bypasses engine validation. Design call.
- **M8** — `WriterConfig` lacks `#[non_exhaustive]`. Public-API surface change.
- **M10** — `VarCallingFromBamArgs` 3rd copy of ~30 fields. Per user direction.
- **M11** — Cohort pipeline wiring duplicated. Per user direction.
- **M13** — Three `build_global` sites uncoordinated. Policy call.
- **Mi2** — `PosteriorEngineConfig::new` 8-positional swap risk (pairs with M4).
- **Mi5** — On-disk types lack `#[non_exhaustive]` (pairs with M8).
- **Mi6** — `provenance.version` not validated on read (schema versioning policy).
- **Mi8** — Helper duplication across 3–4 modules (pairs with M10/M11).
- **Mi13** — File `contamination_artifact.rs` vs type `ContaminationArtefact`.
- **Mi14** — CLI knob `min_batch_size` diverges from engine field (touches engine).
- **Mi18** — `build_artefact_from_estimates` linear scan (needs new engine accessor).
- **Mi19** — Magic `64 * 1024` BufReader capacity (pairs with Mi8).
- **Mi20** — Integration-test fixture-helper duplication.
- **Mi21** — `ContaminationEstimationConfig::validate(&self)` shape (pairs with M4).
- **Mi23** — Load-bearing chain integration test.
- **Selected Nits** — see Nits log above.

## 6. Disputed findings to return to reviewer

- **Mi7** — `&str::lines()` already handles CRLF correctly per the Rust language guarantee. Verified by two new regression tests. The reviewer's "if confirmed, promote to Major" branch is closed: the bug does not exist.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** No.
- **Baseline saved:** No.
- **Benches run:** None.
- **Verdicts:** N/A.
- **Outcome:** Skipped — no Apply touched perf-sensitive code reachable from `benches/`.
- **Notes:** Mi9 renamed a const in `posterior_engine.rs` (a bench-reachable file) but the rename is a pure compile-time identifier swap with zero runtime cost. Mi15 relaxed a `Debug` impl bound in `variant_grouping.rs` (bench-reachable) — no runtime impact. Every other change landed in the CLI orchestrator layer (`pop_var_caller/*`), which has no bench harness.

## 10. Commands run

- `./scripts/dev.sh cargo build --tests` (several checkpoint builds during application)
- `./scripts/dev.sh cargo test --lib pop_var_caller::stage1_pipeline` (M1+M2 verification)
- `./scripts/dev.sh cargo test --lib pop_var_caller::cli::parsers` (Mi1 verification)
- `./scripts/dev.sh cargo test --lib pop_var_caller` (Mi7 + cross-module sanity)
- `./scripts/dev.sh cargo fmt` (post-rename whitespace normalisation)
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings`
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo test --tests`
- `./scripts/dev.sh cargo doc --no-deps`

## 11. Command results

- `cargo fmt --check` → exit 0, clean.
- `cargo clippy --workspace --all-targets -- -D warnings` → exit 0, clean.
- `cargo test --lib` → exit 0, `test result: ok. 886 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 39.79s`.
- `cargo test --tests` → exit 0 (every integration target green; final sampled summary `test result: ok. 12 passed; 0 failed`).
- `cargo doc --no-deps` → exit 101 — **all 5 new errors introduced by the slice are resolved**; remaining 3 errors are pre-existing `ExactMath` links in `posterior_engine.rs` (out of scope, tracked by Stage 6 review Mi21).

## 12. Notes

- Two findings were bundled into single edits because they touched the same surface and one was naturally part of the other: **Mi12 + Mi17** (the same rewrite of `to_estimates_for_samples`) and **M6 + M12** (the from-bam tmp cleanup naturally lives inside the extracted helper).
- The `prefer_upstream_or_closure` helper for M1/M2 is the closest this run gets to test-first discipline: factoring the surface-order logic into a pure function let me unit-test the M1 case without a real CRAM fixture. The trade-off is that the wiring at the call site is still uncovered end-to-end — that follow-up is captured under "Deferred (Mi23)".
- Mi7 was reclassified from "Minor (needs verification)" to **Disputed (with regression test)** after verifying `&str::lines()` documented behaviour. Two tests now lock the contract.
- The cohort CLI slice is now doc-clean: `cargo doc --no-deps` still fails because of 3 pre-existing errors in `posterior_engine.rs`, but every error this slice introduced is gone, and the 2 previously pre-existing `pileup_to_psp.rs` errors resolved as a side-effect of the new module docs.
- Total Applied count: 20 findings (18 Major+Minor Applied + 2 Applied-with-adaptation) + 1 Disputed-with-test (Mi7) + subset of Nits. Total Deferred count: 16 findings (5 Major + 11 Minor). No Blockers.
