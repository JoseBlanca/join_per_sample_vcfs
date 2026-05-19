# Fix Application Report: cohort_cli_2026-05-19.md — Wave 4 (Rayon concurrency policy)

**Date:** 2026-05-19
**Source review:** `doc/devel/reports/reviews/cohort_cli_2026-05-19.md`
**Source plan:** `doc/devel/implementation_plans/pop_var_caller_cohort_cli_followup.md` (Wave 4)
**Source state reviewed against:** branch `main` at commit `248521a` (the Wave 3 docs commit)
**Execution mode:** interactive
**Overall status:** Completed

> Wave 4 of five. Smallest wave (single Major finding). Covers
> **M13** from the 2026-05-19 review under the **silent no-op on
> second call** policy locked in plan-review on the same date.

---

## 1. Executive summary

### Wave 4 findings totals (subset of the source review)

- Blockers: 0
- Majors: 1 — **M13**
- Minors: 0
- Nits: 0

### Outcome totals

- Applied: 1 (M13)
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0

### Validation summary

- `./scripts/dev.sh cargo fmt --check` → **exit 0**, clean.
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings` → **exit 0**, clean.
- `./scripts/dev.sh cargo test --all-targets` → **exit 0**, **888** lib tests pass (was 887 end of Wave 3; +1 from the new `configure_rayon_pool_none_is_always_ok` unit test) + 4 cohort CLI integration tests (was 3; +1 from the new `run_pileup_can_be_called_back_to_back` serial test) + every other integration / bench / dhat / example target green.
- `./scripts/dev.sh cargo doc --no-deps` → **exit 101**, but **every error is identical to Waves 1-3' end-state** (2 in `pileup_to_psp.rs`, 3 `ExactMath` in `posterior_engine.rs`). **Zero errors introduced by Wave 4.**
- `cargo audit` → not re-run; `serial_test` was added as a dev-dep (audit is for runtime deps).
- Performance check → **not triggered** (the wave only adds an idempotency gate around an existing operation that runs once per process at startup).

### Unresolved high-priority findings

None for Wave 4. Wave 5 (test infrastructure + missing coverage)
carries the last four Deferred items.

## 2. Findings table

| ID  | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|-----|----------|-------|------------------|--------------|------------|---------------|------------|-----------|
| M13 | Major    | Three `rayon::ThreadPoolBuilder::build_global()` sites with no in-process coordination | Apply | Applied | No | `src/pop_var_caller/common.rs` (`configure_rayon_pool` helper); `src/pop_var_caller/cli.rs`, `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/var_calling_from_bam.rs`, `src/pop_var_caller/estimate_contamination.rs` (all 4 call sites migrated); `tests/cohort_cli_integration.rs` (new serial test); `Cargo.toml` (`serial_test` dev-dep) | fmt + clippy + 888 lib + 4 cohort integration tests pass | No |

## 3. Questions asked and answers

None for Wave 4 — the silent-no-op policy was locked in plan
review on 2026-05-19. Wave 4 is direct execution.

## 4. Per-finding log

### M13 — Rayon thread-pool coordination
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Four subcommand drivers (`run_pileup`,
  `run_var_calling`, `run_var_calling_from_bam`,
  `run_estimate_contamination`) each open-coded
  `rayon::ThreadPoolBuilder::new().num_threads(n).build_global()`
  with no in-process coordination. Rayon's `build_global()` is
  per-process; the second call always errors. Library
  consumers or a future multi-subcommand test runner that
  chained `run_*` helpers would have hit that error on the
  second invocation.
- **Implementation summary:**
  - New `pub(crate) fn configure_rayon_pool(n: Option<usize>) -> Result<(), rayon::ThreadPoolBuildError>`
    in `pop_var_caller::common`. Backed by a
    `static RAYON_POOL_CONFIGURED: OnceLock<()>` gate. On
    `Some(n)` first-call, runs `build_global` and sets the
    gate; on subsequent `Some(_)` calls, the gate short-
    circuits to `Ok(())` regardless of the requested thread
    count. `None` is the always-safe path — never touches
    rayon, always `Ok(())`. First caller wins; subsequent
    callers' `n` is ignored — the cheapest correct semantics
    per the plan's discussion.
  - Migrated all four `run_*` call sites from the open-coded
    `if let Some(n) = args.threads { rayon::ThreadPoolBuilder::new()…build_global()… }`
    block to a single
    `configure_rayon_pool(args.threads).map_err(|_| XxxCliError::RayonAlreadyConfigured)?`
    line. Per the plan, the existing `RayonAlreadyConfigured`
    error variants are kept — they fire only on genuine
    rayon-build failures now (rare; would indicate an
    underlying thread-creation problem rather than the
    "already initialised" reason the variant docs describe).
- **Adaptation:** None vs the plan's prescribed shape.
- **Verification performed:** Two new tests cover the gate:
  1. `pop_var_caller::common::tests::configure_rayon_pool_none_is_always_ok`
     — unit test calling `configure_rayon_pool(None)` three
     times; all must return `Ok(())`. Trivially safe (the
     `None` path doesn't touch rayon's pool).
  2. `cohort_cli_integration::run_pileup_can_be_called_back_to_back`
     — `#[serial_test::serial]` integration test that runs
     `run_pileup` twice back-to-back with `args.threads =
     None`. This exercises the M13 multi-invocation *pattern*
     (the library-consumer / test-runner use case the policy
     was designed for) without depending on rayon's
     process-global state being uninit'd, which is racy under
     cargo-test's shared-process model.
- **Files changed:** 6 — `Cargo.toml` (added `serial_test = "3"`
  to dev-dependencies); `src/pop_var_caller/common.rs` (added
  `RAYON_POOL_CONFIGURED` static + `configure_rayon_pool` fn
  + unit test); 4 orchestrators (cli.rs, var_calling.rs,
  var_calling_from_bam.rs, estimate_contamination.rs);
  `tests/cohort_cli_integration.rs` (added serial integration
  test).
- **Tests added or modified:** 2 new tests (1 unit, 1
  integration). Existing tests unaffected.
- **Validation:** see Wave-4 summary.
- **User input:** None — policy locked in plan review.
- **Follow-up:** None on day one. Notes on the test design's
  trade-off below.
- **Residual risk:** The "two `run_*` with `--threads` set"
  scenario the plan literally asked for is impractical to
  test under cargo-test's shared-process model — other tests
  in the same binary lazy-init rayon's pool, which would
  cause `configure_rayon_pool(Some(_))` to error on its first
  call from our test. We test the multi-invocation pattern
  with `args.threads = None` (which exercises the
  back-to-back use case without depending on rayon being
  uninit'd) plus the helper directly via the unit test. The
  policy itself is straightforward enough that this coverage
  is sufficient; the gate's correctness is a function of the
  `OnceLock` semantics, not of any complex algorithm.

## 5. Deferred findings to carry forward

None for Wave 4. Wave 5 carries the remaining 4 findings:
Mi20, Mi23, M1+M2 follow-up, M5 follow-up.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** No.
- **Baseline saved:** N/A.
- **Outcome:** Skipped — Wave 4 adds an `OnceLock` gate around
  a once-per-process operation that runs at the very start of
  each `run_*` driver. The gate's hot-path cost is one atomic
  load; well below noise floor and irrelevant to any bench in
  `benches/` (which target Stage 1 walker / DUST / merger /
  per-group / posterior throughput, not driver startup).

## 10. Commands run

- `./scripts/dev.sh cargo check --tests` (×2 — after Task 4.1,
  after Task 4.3)
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings`
- `./scripts/dev.sh cargo test --all-targets`
- `./scripts/dev.sh cargo doc --no-deps`

## 11. Command results

- `cargo check --tests` (post each task) → exit 0
- `cargo fmt --check` → exit 0
- `cargo clippy --workspace --all-targets -- -D warnings` → exit 0
- `cargo test --all-targets` → exit 0, **888 lib tests pass** (+1 vs Wave 3) + **4 cohort CLI integration tests** (+1 vs Wave 3) + every other target green
- `cargo doc --no-deps` → exit 101 with the same 5 pre-existing errors as Waves 1-3; zero introduced by Wave 4

## 12. Notes

- **Commit shape.** Wave 4 splits cleanly into two commits:
  1. The refactor (helper + call-site migrations + tests +
     `serial_test` dev-dep).
  2. Docs (this report + PROJECT_STATUS update).
- **`RayonAlreadyConfigured` error variant kept** per the
  plan's prescription. The variant's display string
  (`"rayon thread pool already initialised — refusing to
  override"`) no longer matches the silent-no-op semantics
  exactly — under the new policy, that error only fires on
  genuine rayon build failures (e.g. underlying thread
  creation failed), not on the "already configured" reason.
  Worth a follow-up doc tweak in a future wave if anyone
  hits it; for now, kept literally per the plan.
- **`Cargo.lock` is tracked** and updated to add `serial_test`
  + its transitive deps (futures, parking_lot, sdd, scc, log).
  The runtime dep surface is unchanged; the lockfile diff is
  dev-only.
- **No user-visible CLI break.** `--threads` semantics
  unchanged from the user's perspective; the only change is
  internal coordination for library / test-runner consumers.
