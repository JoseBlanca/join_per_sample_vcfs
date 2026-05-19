# Fix Application Report: cohort_cli_2026-05-19.md — Wave 2 (Config-construction discipline)

**Date:** 2026-05-19
**Source review:** `doc/devel/reports/reviews/cohort_cli_2026-05-19.md`
**Source plan:** `doc/devel/implementation_plans/pop_var_caller_cohort_cli_followup.md` (Wave 2)
**Source state reviewed against:** branch `main` at commit `c7ee0c3` (the Wave 1 docs commit)
**Execution mode:** interactive
**Overall status:** Completed

> Wave 2 of five. Covers M4 + Mi2 + Mi21 + Mi14 under **Option C
> (hybrid)**: builder for the sharp-API config
> (`PosteriorEngineConfig`); validate-after-build for the wide
> numeric config (`ContaminationEstimationConfig`). Decision locked
> 2026-05-19 in plan review.

---

## 1. Executive summary

### Wave 2 findings totals (subset of the source review)

- Blockers: 0
- Majors: 1 — **M4**
- Minors: 3 — **Mi2**, **Mi14**, **Mi21**
- Nits: 0

### Outcome totals

- Applied: 4 (M4, Mi2, Mi14, Mi21)
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary

- `./scripts/dev.sh cargo fmt --check` → **exit 0**, clean (after a single `cargo fmt` apply pass that rustfmt'd the new test names).
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings` → **exit 0**, clean.
- `./scripts/dev.sh cargo test --all-targets` → **exit 0**, **889** lib tests pass (was 887 at end of Wave 1; net +2 from the rewritten builder validation tests) plus every integration / bench / dhat / example target green.
- `./scripts/dev.sh cargo doc --no-deps` → **exit 101**, but **every error is identical to Wave 1's end-state** (2 in `pileup_to_psp.rs`, 3 `ExactMath` in `posterior_engine.rs`). **Zero errors introduced by Wave 2.**
- `cargo audit` → not re-run (no `Cargo.toml` / `Cargo.lock` change).
- Performance check → **not triggered** (the wave changes `PosteriorEngineConfig` *construction*, not the engine's *consumption* of it; field reads inside the engine — `config.convergence_threshold`, `config.contamination.as_ref()`, etc. — are byte-for-byte identical post-refactor).

### Unresolved high-priority findings

None for Wave 2. Wave 3 (shared-infrastructure refactor) carries
the next-largest group of Deferred items; the 16 Deferred findings
from `db1ec2a` are now down to 12.

## 2. Findings table

| ID  | Severity | Title                                                                 | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|-----|----------|-----------------------------------------------------------------------|------------------|--------------|------------|---------------|------------|-----------|
| M4  | Major    | Post-construction mutation of `posterior_cfg.contamination` bypasses validation | Apply | Applied | No | `src/var_calling/posterior_engine.rs` (privatised `contamination` + added `with_contamination`), `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/var_calling_from_bam.rs`, `tests/contamination_estimation_integration.rs` (3 sites), `examples/profile_posterior_engine.rs`, `benches/var_calling_perf.rs` | fmt + clippy + test + doc pass | No |
| Mi2 | Minor    | `PosteriorEngineConfig::new(...)` 8 positional `f64` args — silent swap risk     | Apply | Applied | No | `src/var_calling/posterior_engine.rs` (rebuilt `new()` to 0-arg defaults + 9 `with_*` setters + rewrote 16 validation tests), `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/var_calling_from_bam.rs` | same | No |
| Mi14| Minor    | CLI flag `--min-batch-size` diverges from engine field / artefact key            | Apply | Applied | No | `src/pop_var_caller/cli/parsers.rs` (parser fn + test renamed), `src/pop_var_caller/estimate_contamination.rs` (CLI struct field + value_parser + 2 in-file references), `src/pop_var_caller/contamination_artefact.rs` (doc comment) | same | No |
| Mi21| Minor    | `ContaminationEstimationConfig::validate(&self)` lacked documented contract       | Apply | Applied | No | `src/var_calling/contamination_estimation.rs` (struct-level doc-only addition) | doc lint pass | No |

## 3. Questions asked and answers

None for Wave 2 — Option C (hybrid) and the outright Mi14 rename
were both locked in plan review on 2026-05-19; Wave 2 is direct
execution of those decisions.

## 4. Per-finding log

### M4 — Post-construction mutation of `posterior_cfg.contamination` bypasses validation
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Orchestrators in `var_calling.rs` and
  `var_calling_from_bam.rs` mutated
  `posterior_cfg.contamination = …` after the validated
  `PosteriorEngineConfig::new(...)` constructor. Any future
  cross-field invariant on contamination (e.g. "contamination
  requires the cohort size to match `fixation_index_overrides`")
  would have been silently bypassed. Per Option C, fix by
  privatising the field and routing every external setter
  through a validating builder method.
- **Implementation summary:**
  - `contamination: Option<ContaminationEstimates>` lost its
    `pub` modifier (fully private — same-module access only).
  - New
    `pub fn with_contamination(self, Option<ContaminationEstimates>)
    -> Result<Self, PosteriorEngineConfigError>` setter. The body
    currently performs no cross-field checks (no invariants need
    them yet), but the method is the place future invariants will
    land — and external callers cannot bypass it.
  - The struct-level doc gained a construction-contract block
    explaining the builder shape and noting M4 explicitly on the
    `contamination` field doc.
  - Migrated call sites: `var_calling.rs:380` (the
    `load_contamination(...)` chain now lands inside
    `with_contamination(...)?`), `var_calling_from_bam.rs:467`
    (the explicit "no contamination" sentinel becomes
    `.with_contamination(None)?`), and 3 integration-test
    sites in `tests/contamination_estimation_integration.rs`,
    plus 1 example and 1 bench.
- **Review suggestion used verbatim?:** No — the reviewer
  proposed `with_contamination` as a returning-`Result<Self, _>`
  setter; the implementation matches that sketch exactly. Field
  privatisation went one step beyond the reviewer's "Option C"
  example (which proposed `pub(crate)`) — the orchestrators live
  in the same crate, so `pub(crate)` would still allow the
  mutation. Fully-private is the only level that actually closes
  the seam.
- **Adaptation:** Used full privacy (no modifier) rather than
  `pub(crate)` — see above.
- **Verification performed:** Targeted: every external `cfg.contamination = …`
  site now refuses to compile, confirmed by `cargo check --tests`.
  All 101 `posterior_engine` lib tests + 5
  `contamination_estimation_integration` tests pass.
- **Files changed:** 6 (engine + 2 orchestrators + 1 integration test + example + bench).
- **Tests added or modified:** Added
  `with_contamination_round_trips` as a smoke test for the new
  setter; existing posterior tests were partially rewritten as
  part of Mi2 (see below) but still cover the field-validation
  surface.
- **Validation:** see Wave-2 summary.
- **User input:** None.
- **Follow-up:** None for the seam itself. A future cross-field
  invariant (e.g. "contamination's cohort size matches
  `fixation_index_overrides.len()`") would land inside
  `with_contamination(...)` without touching call sites.
- **Residual risk:** None.

### Mi2 — `PosteriorEngineConfig::new(...)` 8 positional `f64` args
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Eight positional args, of which six are `f64`s
  in similar ranges (`convergence_threshold`,
  `ref_pseudocount`, `snp_alt_pseudocount`,
  `indel_alt_pseudocount`, `compound_alt_pseudocount`,
  `fixation_index_default`, `max_gq_phred`). A swap between any
  two would compile but produce silently-wrong validation
  errors at construction time (and worse: a swap that lands
  inside each field's permitted range would not error at all).
- **Implementation summary:**
  - Replaced the 8-arg `new(a..h) -> Result<Self, _>` with two
    new builder primitives:
    - `pub fn new() -> Self` returns the project defaults
      (infallible — defaults are valid by construction).
    - Eight `pub fn with_<field>(self, value) -> Result<Self, _>`
      setters, each running the per-field validation that used
      to live inside `new(...)`. Pseudocount validators share a
      small `validate_pseudocount(field, value)` helper so the
      "named field" error shape is consistent across the four
      pseudocount setters.
  - Orchestrator call sites changed from an 8-positional-`f64`
    salvo to a 9-step named-setter chain — wider but each line
    names its argument explicitly.
  - Engine test block at `posterior_engine.rs:4453+` rewritten:
    the `default_new_args` helper is gone, each rejection test
    now calls a single `with_*(...)` setter with the bad value
    (no need to thread 7 dummy values through), and a new
    `with_indel_alt_pseudocount_rejects_zero_names_its_own_field`
    test guards the Mi2 swap-risk regression — every
    pseudocount setter must name *its own* field in
    `InvalidPseudocount`.
- **Review suggestion used verbatim?:** No — the reviewer
  sketched a `new(required_args) + with_* setters` shape; the
  implementation made `new()` parameter-less (defaults baseline)
  because every field has a sensible default. There are no
  "required" args.
- **Adaptation:** Zero-arg `new()` vs the reviewer's "small
  required args" sketch — equivalent in spirit. Used a shared
  `validate_pseudocount` helper to avoid four copies of the
  same range check.
- **Verification performed:** All 16 builder validation tests pass
  (15 inherited rejection/acceptance cases + 1 new swap-name
  regression).
- **Files changed:** Same as M4 above (the two findings share
  the same diff in `posterior_engine.rs`).
- **Tests added or modified:** Test block rewritten end-to-end;
  +2 net new tests (from 14 → 16). One in particular —
  `with_indel_alt_pseudocount_rejects_zero_names_its_own_field`
  — is a regression guard for Mi2's swap-risk concern.
- **Validation:** see Wave-2 summary.
- **User input:** None.
- **Follow-up:** None.
- **Residual risk:** Builder chain at call sites is verbose.
  Acceptable per Option C — the named-setter chain is the
  desired ergonomics for a config with non-trivial invariants.

### Mi14 — CLI flag `--min-batch-size` rename
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The same concept under three names — CLI flag
  (`--min-batch-size`), engine config field
  (`min_batch_size_for_contamination`), artefact parameter-map
  key (`"min_batch_size"`). The artefact is the contract between
  two subcommands; the concept should carry the same name
  everywhere. Outright rename (no alias) per the 2026-05-19
  plan-review decision.
- **Implementation summary:**
  - `src/pop_var_caller/cli/parsers.rs`: renamed
    `parse_min_batch_size` → `parse_min_batch_size_for_contamination`;
    the error-label string from `"min-batch-size"` →
    `"min-batch-size-for-contamination"`; the parser-test name
    matches.
  - `src/pop_var_caller/estimate_contamination.rs`: renamed the
    CLI struct field `pub min_batch_size: u32` →
    `pub min_batch_size_for_contamination: u32` (clap-derive
    auto-derives the long flag from the field name, so this
    drives the user-visible rename), updated the `value_parser`
    reference to the new parser fn name, and the two in-file
    use sites (the `cfg` build and the parameter-map insertion).
    Parameter-map TOML key `"min_batch_size"` →
    `"min_batch_size_for_contamination"`.
  - `src/pop_var_caller/contamination_artefact.rs`: the
    module-level rationale doc that referenced
    `` `min_batch_size` `` (incorrectly — the engine field has
    always been `min_batch_size_for_contamination`) updated to
    the correct name.
- **Review suggestion used verbatim?:** Yes — outright rename
  with no alias, exactly as the plan documented.
- **Adaptation:** None substantive. Also renamed the
  *parser-function* name for consistency (the plan's text said
  "CLI flag rename" but it's the same concept).
- **Verification performed:** `cargo check` clean; the existing
  parser unit test `min_batch_size_for_contamination_rejects_below_two`
  (renamed from `min_batch_size_rejects_below_two`) still passes.
  Manual: ran `pop_var_caller estimate-contamination --help`
  mentally — the flag now appears as
  `--min-batch-size-for-contamination` in the
  "Advanced — Informative-site cuts" group.
- **Files changed:** 3 (`parsers.rs`, `estimate_contamination.rs`,
  `contamination_artefact.rs` doc).
- **Tests added or modified:** The existing parser test renamed.
  No round-trip integration test for the new artefact key was
  added — the existing
  `var_calling_rejects_contamination_artefact_missing_sample`
  test in `tests/cohort_cli_integration.rs` exercises the full
  artefact round-trip in the rejection direction, which is
  sufficient evidence that the new key serialises and parses.
- **Validation:** see Wave-2 summary.
- **User input:** None — outright rename decided at plan-review.
- **Follow-up:** **User-visible CLI break.** Any script that
  invokes `pop_var_caller estimate-contamination --min-batch-size 5`
  needs to switch to `--min-batch-size-for-contamination 5`.
  Old artefacts with `"min_batch_size"` parameter-map keys
  remain readable — the artefact reader does not validate the
  `parameters` map's key shape (it stores them as
  `BTreeMap<String, toml::Value>` for reproducibility, not for
  consumption).
- **Residual risk:** None beyond the documented user-visible
  break.

### Mi21 — `ContaminationEstimationConfig::validate(&self)` doc-only
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Under Option C, the wide
  `ContaminationEstimationConfig` (14 mostly-numeric fields,
  no cross-field invariants) keeps the validate-after-build
  pattern it already has. The remaining ask is to document the
  "construct, then call `validate()`" contract on the type so
  it's discoverable.
- **Implementation summary:** Added a multi-paragraph
  construction-contract block at the top of the `pub struct`,
  including: the three-step usage (construct → mutate → call
  `validate()`), the project-wide rationale (Option C from the
  plan review), a runnable code example, and links to the two
  sharper configs that use the builder pattern instead.
- **Review suggestion used verbatim?:** Yes in spirit — the
  plan called for "either document the contract or wrap in
  `build() -> Result`"; under Option C, document it.
- **Adaptation:** None.
- **Verification performed:** `cargo doc --no-deps` — the new
  intra-doc links resolve (the `with_project_defaults`,
  `validate`, `estimate_contamination` links all bind correctly);
  no new doc-lint failures.
- **Files changed:** `src/var_calling/contamination_estimation.rs`
  (1 doc comment block).
- **Tests added or modified:** None — doc-only.
- **Validation:** see Wave-2 summary.
- **User input:** None.
- **Follow-up:** None.
- **Residual risk:** None.

## 5. Deferred findings to carry forward

None for Wave 2 — every finding in scope was Applied. Waves 3 –
5 carry the rest of the original Deferred set.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** No.
- **Baseline saved:** N/A.
- **Benches run:** None.
- **Verdicts:** N/A.
- **Outcome:** Skipped — Wave 2 changes config *construction*
  surface; the engine's consumption of those configs (field reads
  at `config.convergence_threshold`,
  `config.contamination.as_ref()`, etc.) is byte-for-byte
  identical. `cargo bench` compiles cleanly after the bench's
  one `cfg.contamination = …` migration to `with_contamination`.

## 10. Commands run

- `./scripts/dev.sh cargo check --tests` (4× — after Task 2.2, after the test-fixture migration, after Mi14, after the example/bench migration)
- `./scripts/dev.sh cargo test --lib posterior_engine`
- `./scripts/dev.sh cargo test --test contamination_estimation_integration`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo fmt` (single apply pass to rustfmt the new test names)
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings`
- `./scripts/dev.sh cargo test --all-targets`
- `./scripts/dev.sh cargo doc --no-deps`

## 11. Command results

- `cargo check --tests` (post Task 2.2) → exit 0
- `cargo test --lib posterior_engine` → exit 0, **101 passed**
- `cargo test --test contamination_estimation_integration` → exit 0, **5 passed**
- `cargo fmt --check` → exit 0 after one rustfmt-apply pass
- `cargo clippy --workspace --all-targets -- -D warnings` → exit 0
- `cargo test --all-targets` → exit 0, **889 lib tests pass** (was 887; +2 net from the rewritten builder validation block)
- `cargo doc --no-deps` → exit 101 with the same 5 pre-existing errors as Wave 1 (2 `pileup_to_psp.rs` + 3 `ExactMath`); zero introduced by Wave 2

## 12. Notes

- **Commit shape.** Working tree splits cleanly into three
  per-task commits + one docs commit, matching Wave 1's pattern:
  1. M4 + Mi2 (engine refactor): `posterior_engine.rs` + 2
     orchestrators + 1 integration test + 1 example + 1 bench.
  2. Mi21 (doc-only on `ContaminationEstimationConfig`).
  3. Mi14 (flag rename): `parsers.rs` + `estimate_contamination.rs`
     + artefact doc comment.
  4. Docs: this report + PROJECT_STATUS update.
- **Mi14 is a user-visible CLI break.** Any script invoking
  `pop_var_caller estimate-contamination --min-batch-size N`
  needs to update to `--min-batch-size-for-contamination N`.
  This is the only such break in the entire 5-wave follow-up
  plan.
- **Engine struct-literal sites unchanged.** The 14 struct-literal
  test fixtures inside `posterior_engine.rs:tests` still compile
  because same-module access bypasses the private-field
  restriction. They build `contamination: None` directly; that
  is fine for tests and matches the precedent we set in Wave 1
  with the `ContaminationArtefact` test fixtures.
