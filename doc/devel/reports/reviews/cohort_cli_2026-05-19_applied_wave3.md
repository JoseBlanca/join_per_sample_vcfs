# Fix Application Report: cohort_cli_2026-05-19.md — Wave 3 (Shared-infrastructure refactor)

**Date:** 2026-05-19
**Source review:** `doc/devel/reports/reviews/cohort_cli_2026-05-19.md`
**Source plan:** `doc/devel/implementation_plans/pop_var_caller_cohort_cli_followup.md` (Wave 3)
**Source state reviewed against:** branch `main` at commit `f44c086` (the Wave 2 docs commit)
**Execution mode:** interactive
**Overall status:** Completed

> Wave 3 of five. The largest mechanical wave (~1500 LOC).
> Covers Mi8, Mi19, M9 follow-up, M10, M11, Mi18 from the
> 2026-05-19 review.

---

## 1. Executive summary

### Wave 3 findings totals (subset of the source review)

- Blockers: 0
- Majors: 2 — **M10**, **M11**
- Minors: 3 — **Mi8**, **Mi18**, **Mi19**
- Doc follow-ups: 1 — **M9 follow-up**
- Nits: 0

### Outcome totals

- Applied: 6 of 6 (M10, M11, Mi8, Mi18, Mi19, M9 follow-up)
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0

### Validation summary

- `./scripts/dev.sh cargo fmt --check` → **exit 0**, clean (after a single `cargo fmt` apply pass that rustfmt'd the wave's new code).
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings` → **exit 0**, clean.
- `./scripts/dev.sh cargo test --all-targets` → **exit 0**, **887** lib tests pass (was 889 end of Wave 2; net **−2** from helper-test deduplication into `pop_var_caller::common::tests`) plus every integration / bench / dhat / example target green.
- `./scripts/dev.sh cargo doc --no-deps` → **exit 101**, but **every error is identical to Wave 1's end-state** (2 in `pileup_to_psp.rs`, 3 `ExactMath` in `posterior_engine.rs`). **Zero errors introduced by Wave 3** after fixing the three `[CohortVcfWriter]` / `[DustFilter]` / `[CramInputError]` doc-link errors that briefly surfaced when their imports moved.
- `cargo audit` → not re-run (no `Cargo.toml` / `Cargo.lock` change).
- Performance check → **not triggered**: the wave is structural code-movement only. The cohort pipeline's runtime data path is byte-for-byte identical to Wave 2 — `drive_cohort_pipeline` runs the same iterator chain, just defined in one place instead of two.

### Unresolved high-priority findings

None for Wave 3. The remaining Deferred set is down to: M13 (rayon
policy, Wave 4) + Mi20 / Mi23 / M1+M2-followup / M5-followup
(Wave 5). 

## 2. Findings table

| ID  | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|-----|----------|-------|------------------|--------------|------------|---------------|------------|-----------|
| Mi8 | Minor    | Helper duplication (basename / format_md5_hex / current_command_line / rfc3339_now / civil_from_days) across 3-4 files | Apply | Applied | No | New `src/pop_var_caller/common.rs`; deleted duplicates and migrated imports across `cli.rs`, `var_calling.rs`, `var_calling_from_bam.rs`, `estimate_contamination.rs` | fmt+clippy+test pass | No |
| Mi19| Minor    | `64 * 1024` BufReader/BufWriter capacity duplicated at 6+ sites | Apply | Applied | No | `src/pop_var_caller/common.rs` (`DEFAULT_BUFFERED_IO_CAPACITY`); migrated 5 in-scope sites (cli.rs, var_calling.rs, estimate_contamination.rs, psp_to_pileup.rs ×2, vcf_writer/sink.rs) | fmt+clippy+test pass | No |
| M9 (followup) | Major (follow-up) | Generalise `ErrorSheddingAdapter` so the from-bam walker shim reuses it | Apply | Applied | No | `src/pop_var_caller/cli/error_bridge.rs` (made generic over `<I, T, E>`); `src/pop_var_caller/var_calling_from_bam.rs` (replaced open-coded `.scan()` with the parametric adapter) | adapter tests pass + 1 new parametric-types regression test | No |
| Mi18| Minor    | `build_artefact_from_estimates` linear-scans for a representative sample per batch | Apply | Applied | No | `src/var_calling/contamination_estimation.rs` (new `q_b_per_batch()` accessor); `src/pop_var_caller/estimate_contamination.rs` (rewrote to index by batch directly) | artefact tests pass | No |
| M10 | Major    | `VarCallingFromBamArgs` duplicates ~30 fields of `VarCallingArgs` + `PileupArgs` | Apply | Applied | No | New `src/pop_var_caller/cli/shared_args.rs` (`Stage1Args` + `CohortPipelineArgs`); all 3 args structs flattened via `#[command(flatten)]`; helpers `cram/baq/walker_config_from_args` unified into `cli.rs` taking `&Stage1Args`; var_calling.rs's body uses `args.cohort.foo`; var_calling_from_bam.rs's body uses `args.stage1.foo` + `args.cohort.foo`; 4 test fixtures migrated | fmt+clippy+test pass | No |
| M11 | Major    | Cohort pipeline wiring duplicated between `var_calling.rs` and `run_cohort_pipeline_for_single_sample` | Apply | Applied | No | New `src/pop_var_caller/cohort_driver.rs` (`CohortPipelineParams` + `drive_cohort_pipeline`); both `run_var_calling` and `run_cohort_pipeline_for_single_sample` delegate | fmt+clippy+test pass | No |

## 3. Questions asked and answers

None for Wave 3 — every finding's shape was locked in the
plan-review on 2026-05-19. Wave 3 is direct execution.

## 4. Per-finding log

### Mi8 + Mi19 — `pop_var_caller::common` module
- **Severity:** Minor + Minor
- **Initial decision:** Apply (both)
- **Final status:** Applied (both)
- **Reasoning:** Five helpers (`basename`, `format_md5_hex`,
  `current_command_line`, `rfc3339_now`, `civil_from_days`)
  duplicated 2-4 places each across `cli.rs`, `var_calling.rs`,
  `var_calling_from_bam.rs`, `estimate_contamination.rs`. Plus the
  64 KiB `BufReader` / `BufWriter` capacity literal at 6+ sites.
  One extracted module covers both.
- **Implementation summary:**
  - New `src/pop_var_caller/common.rs` with `pub(crate)` versions
    of all five helpers and `pub(crate) const DEFAULT_BUFFERED_IO_CAPACITY: usize = 64 * 1024;`
    plus consolidated `#[cfg(test)] mod tests` covering every
    helper.
  - Deleted duplicates from the four files; their test blocks
    lost the helper-specific tests (now centralised in
    `common::tests`).
  - Migrated 5 in-scope `BufReader::with_capacity(64 * 1024, ...)` /
    `BufWriter::with_capacity(64 * 1024, ...)` sites: `cli.rs`,
    `var_calling.rs`, `estimate_contamination.rs`,
    `psp_to_pileup.rs` (×2), `vcf_writer/sink.rs`.
- **Files changed:** 1 new (`common.rs`), 4 trimmed
  (cli.rs / var_calling.rs / var_calling_from_bam.rs /
  estimate_contamination.rs), 2 migrated to the constant
  (`psp_to_pileup.rs`, `vcf_writer/sink.rs`), 1 wired
  (`mod.rs`).
- **Tests added or modified:** 5 helper tests now live exclusively
  in `common::tests` (was 8 spread across 4 files). Plus one new
  `basename_falls_back_to_full_path_on_trailing_slash` smoke test
  that exercises the previously-unguarded fallback branch.
- **Validation:** see Wave-3 summary.

### M9 follow-up — Generalise `ErrorSheddingAdapter`
- **Severity:** Major (follow-up)
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Adapter was hard-coded to
  `Iterator<Item = Result<PreparedRead, CramInputError>>`. The
  from-bam walker shim needed the same shape with different
  types (`PileupRecord` / `WalkerError`) and was open-coding it
  with an inline `.scan()` + `Rc<RefCell<Option<WalkerError>>>`
  stash.
- **Implementation summary:**
  - `ErrorSheddingAdapter<I>` → `ErrorSheddingAdapter<I, T, E>`
    where `I: Iterator<Item = Result<T, E>>`; the handle is now
    `ErrorHandle<E>` (manual `Clone` impl — the auto-derive
    would have added an unwarranted `E: Clone` bound that
    `Rc::clone` doesn't need).
  - Stage 1 CRAM-input call sites (`stage1_pipeline.rs:130, 159`)
    unchanged; Rust infers `T = PreparedRead, E = CramInputError`
    from the input iterator.
  - From-bam walker shim in `var_calling_from_bam.rs` replaced
    its open-coded `.scan()` with
    `ErrorSheddingAdapter::<_, PileupRecord, WalkerError>::new(ctx.walker)`,
    then `.map(Ok)` to wrap the shed-stream items into the
    `Result<_, PspReadError>` shape `PerPositionMerger` wants.
  - Dropped `use std::cell::RefCell` and `use std::rc::Rc` from
    var_calling_from_bam.rs — no longer needed.
- **Tests:** added `parametric_over_item_and_error_types` test
  with `T = i32, E = String` to lock the generalisation in
  (decoupled from any specific domain types).

### Mi18 — `q_b_per_batch` accessor on `ContaminationEstimates`
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The artefact builder needed batch-side access to
  `q_b` rows but the engine only exposed `q_b_for_sample(idx)` —
  sample-side. The builder worked around this by linear-scanning
  `sample_to_batch.iter().position(|b| *b == batch_idx)` to find
  a representative sample, then reading its `q_b`. Plus a
  defensive `[0.0; 3]` fall-through for "no representative
  sample" — impossible per `build_dense_batches`' contract but
  defensible without an engine-side accessor.
- **Implementation summary:**
  - New
    `pub fn q_b_per_batch(&self) -> &[[f64; N_ALLELE_CLASSES]]`
    on `ContaminationEstimates`. Trivial getter that returns the
    `pub(crate)` `q_b_per_batch` field's slice. Doc references
    Mi18.
  - Rewrote `build_artefact_from_estimates` to index `q_b_rows
    [batch_idx]` directly via the new accessor — the linear-scan
    workaround and the `[0.0; 3]` fall-through are both gone.
    The parallel `batch_id_for_idx` and `q_b_rows` lengths
    are guaranteed equal by `build_dense_batches`, so direct
    indexing is total.

### M10 — Shared `Stage1Args` + `CohortPipelineArgs`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Three args structs (`PileupArgs`,
  `VarCallingArgs`, `VarCallingFromBamArgs`) duplicated ~30 fields
  verbatim — three places for the same flag-set. clap-derive's
  `#[command(flatten)]` makes this transparent at the user
  level.
- **Implementation summary:**
  - New `src/pop_var_caller/cli/shared_args.rs`:
    - `Stage1Args` (16 fields: CRAM-input filters, BAQ HMM,
      pileup walker) — flattened into `PileupArgs` and
      `VarCallingFromBamArgs`.
    - `CohortPipelineArgs` (13 fields: DUST, variant grouping,
      per-group merger, posterior engine, ploidy, VCF writer
      `emit_gp`) — flattened into `VarCallingArgs` and
      `VarCallingFromBamArgs`.
  - `PileupArgs` reduced from 22 fields → 5 (reference,
    output, threads, crams, `stage1`).
  - `VarCallingArgs` reduced from 20 fields → 7 (reference,
    output, threads, psp_files, contamination_estimates,
    no_complexity_filter, `cohort`).
  - `VarCallingFromBamArgs` reduced from 34 fields → 7
    (reference, output, threads, crams, no_complexity_filter,
    `stage1`, `cohort`).
  - The `cram_config_from_args`, `baq_config_from_args`,
    `walker_config_from_args` helpers in `cli.rs` switched
    signature from `&PileupArgs` → `&Stage1Args` and gained
    `pub(super)` visibility. The duplicate copies in
    var_calling_from_bam.rs (~40 LOC) are gone — the from-bam
    orchestrator calls the cli.rs helpers with `&args.stage1`.
  - `effective_parameters` and `build_writer_header` also moved
    to `&Stage1Args`.
  - Body migrations: `args.foo` → `args.stage1.foo` or
    `args.cohort.foo` in both orchestrators. Local rebinds
    (`let stage1 = &args.stage1; let cohort = &args.cohort;`)
    at the top of each orchestrator keep the body diff small.
  - 4 integration / unit test fixtures migrated to the nested
    struct shape (cohort_cli_integration.rs ×3,
    pileup_cli_integration.rs ×1, cli.rs `default_args` builder).
- **Adaptation:** Helper unification (cram / baq / walker
  config_from_args) went beyond the plan's "look one field
  deeper" prescription — the unified `&Stage1Args` signature
  dropped a second axis of duplication that fell out naturally.
- **Files changed:** 1 new (`shared_args.rs`); 5 modified
  (cli.rs, var_calling.rs, var_calling_from_bam.rs, mod.rs,
  cohort_cli_integration.rs, pileup_cli_integration.rs).

### M11 — `drive_cohort_pipeline` shared driver
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Same iterator chain (optional DUST → grouper →
  per-group merger → posterior engine → VCF writer + tmp-rename
  cleanup) duplicated between `run_var_calling` and
  `run_cohort_pipeline_for_single_sample`.
- **Implementation summary:**
  - New `src/pop_var_caller/cohort_driver.rs` with:
    - `pub(crate) struct CohortPipelineParams` — bundle of
      per-stage configs + shared resources (`SharedRefFetcher`,
      chromosomes table, `no_complexity_filter` flag).
    - `pub(crate) fn drive_cohort_pipeline<M, E>(...) -> Result<u64, E>` —
      generic over the merger type (so .psp-based and
      walker-based mergers share the same code) and the typed
      error enum (`VarCallingCliError` /
      `VarCallingFromBamCliError`); both enums already have the
      five required `From<...>` impls.
  - Lifetime quirk: the boxed `Iterator + '_` ties the DUST
    branch's borrow of the local `fetcher` to the function
    body's scope rather than the merger's lifetime. Documented
    inline.
  - Both `run_var_calling` (var_calling.rs) and
    `run_cohort_pipeline_for_single_sample` (var_calling_from_bam.rs)
    now construct `CohortPipelineParams` and delegate to the
    driver. The from-bam wrapper retains the walker-error
    surfacing step (which is from-bam-specific).
- **Validation:** see Wave-3 summary. The drive_cohort_pipeline
  function exercises end-to-end via the existing happy-path
  integration tests (`var_calling_happy_path_three_samples`,
  `var_calling_from_bam_happy_path`) — both pass post-refactor.

## 5. Deferred findings to carry forward

None for Wave 3. Wave 4 (rayon policy, M13) and Wave 5
(test infrastructure + missing coverage) carry the remaining 5
findings.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** No.
- **Baseline saved:** N/A.
- **Outcome:** Skipped — the wave is structural code-movement
  only. `drive_cohort_pipeline` runs the same iterator chain that
  was inlined in two places before. The cohort pipeline's runtime
  data path is byte-for-byte identical. clap's
  `#[command(flatten)]` is a compile-time expansion; the parsed
  arg shape is unchanged at runtime.

## 10. Commands run

- `./scripts/dev.sh cargo check --tests` (×7 — once after each task)
- `./scripts/dev.sh cargo fmt`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings`
- `./scripts/dev.sh cargo test --all-targets`
- `./scripts/dev.sh cargo doc --no-deps`

## 11. Command results

- `cargo check --tests` (post each task) → exit 0
- `cargo fmt --check` → exit 0 after one rustfmt-apply pass
- `cargo clippy --workspace --all-targets -- -D warnings` → exit 0
- `cargo test --all-targets` → exit 0, **887 lib tests pass** + every integration / bench / dhat / example target green
- `cargo doc --no-deps` → exit 101 with the same 5 pre-existing errors as Wave 1 / Wave 2 (2 `pileup_to_psp.rs` + 3 `ExactMath`); zero introduced by Wave 3 after fixing the 3 transient `[CohortVcfWriter]` / `[DustFilter]` / `[CramInputError]` doc-link errors that surfaced when their imports moved

## 12. Notes

- **Commit shape.** Working tree splits cleanly into five
  per-task commits + one docs commit, matching Waves 1 / 2's
  pattern:
  1. Task 3.1 — `pop_var_caller::common` extraction (Mi8 + Mi19).
  2. Task 3.2 — `ErrorSheddingAdapter<I, T, E>` generalisation
     (M9 follow-up).
  3. Task 3.3 — `q_b_per_batch()` engine accessor + artefact
     builder rewrite (Mi18).
  4. Task 3.4 — `Stage1Args` + `CohortPipelineArgs` via clap
     flatten + helper unification (M10).
  5. Task 3.5 — `drive_cohort_pipeline` shared driver (M11).
  6. Docs: this report + PROJECT_STATUS update.
- **No user-visible CLI break.** clap-derive's
  `#[command(flatten)]` is transparent at the help-text and
  parsing layers; `--help` output and CLI grammar are unchanged.
- **Cargo.toml unchanged.** No new dev-dependencies; no Cargo
  surface change.
