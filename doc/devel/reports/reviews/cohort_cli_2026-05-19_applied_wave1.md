# Fix Application Report: cohort_cli_2026-05-19.md — Wave 1 (Public-API hygiene)

**Date:** 2026-05-19
**Source review:** `doc/devel/reports/reviews/cohort_cli_2026-05-19.md`
**Source plan:** `doc/devel/implementation_plans/pop_var_caller_cohort_cli_followup.md` (Wave 1)
**Source state reviewed against:** branch `main` at commit `db1ec2a` (the prior fix-application pass)
**Execution mode:** interactive
**Overall status:** Completed

> Per the follow-up plan's §Validation note, the deferred fixes ship
> as five per-wave reports (`_wave1` … `_wave5`). This report covers
> Wave 1 only: M8, Mi5, Mi6, Mi13.

---

## 1. Executive summary

### Wave 1 findings totals (subset of the source review)

- Blockers: 0
- Majors: 1 — **M8**
- Minors: 3 — **Mi5**, **Mi6**, **Mi13**
- Nits: 0

### Outcome totals

- Applied: 4 (M8, Mi5, Mi6, Mi13)
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary

- `./scripts/dev.sh cargo fmt --check` → **exit 0**, clean.
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings` → **exit 0**, clean.
- `./scripts/dev.sh cargo test --all-targets` → **exit 0**, **887** lib tests pass (was 886 before Wave 1; the new
  `read_rejects_artefact_with_unknown_version` brings the count up by one) plus every integration / bench / dhat / example target green.
- `./scripts/dev.sh cargo doc --no-deps` → **exit 101**, but **every error is pre-existing** (confirmed by a stash-based bisection of the unmodified `db1ec2a` tree):
  - 2 × `pileup_to_psp.rs:5-6` — `[pileup::PileupWalker]` / `[psp::writer::PspWriter]` (pre-existing; previous fix-application run claimed they had auto-resolved, but the stash bisection in this run shows they are present on the un-edited tree).
  - 3 × `posterior_engine.rs:725 / 732 / 753` — `[ExactMath]` (pre-existing; tracked by Stage 6 review Mi21).
  - **No errors introduced by Wave 1.** The two `[Self::validate]` / `[Self::read]` errors briefly introduced by an initial draft of the `SUPPORTED_VERSIONS` doc were caught and fixed before the wave concluded (the const sits at module scope, so `Self::` resolves to the const, not the impl block).
- `cargo audit` → not run this wave; no Cargo.toml or Cargo.lock change so the dependency surface is unchanged from `db1ec2a`.
- Performance check → **not triggered** (no Apply touches code reachable from `benches/`; the wave is `#[non_exhaustive]` attributes, a builder setter, a one-time `validate()` gate, and a file rename).

### Unresolved high-priority findings

None for Wave 1. The 12 still-Deferred findings carried forward
from `db1ec2a` (5 Majors + 7 paired Minors and follow-ups for
Waves 2 – 5) remain Deferred under their own per-wave reports.

## 2. Findings table

| ID  | Severity | Title                                                                 | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|-----|----------|-----------------------------------------------------------------------|------------------|--------------|------------|---------------|------------|-----------|
| M8  | Major    | `WriterConfig` lacks `#[non_exhaustive]`; add builder for `emit_gp`   | Apply            | Applied      | No         | `src/var_calling/vcf_writer/mod.rs`, `src/var_calling/vcf_writer/header.rs`, `src/var_calling/vcf_writer/record_encode.rs`, `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/var_calling_from_bam.rs`, `tests/cohort_vcf_writer_integration.rs` | fmt + clippy + test + doctest pass | No |
| Mi5 | Minor    | `#[non_exhaustive]` missing on `ContaminationArtefact` + 4 sub-structs | Apply            | Applied      | No         | `src/pop_var_caller/contamination_artefact.rs`, `tests/cohort_cli_integration.rs` | fmt + clippy + test pass | No |
| Mi6 | Minor    | `provenance.version` recorded but never validated on read              | Apply            | Applied      | No         | `src/pop_var_caller/contamination_artefact.rs` (added `SUPPORTED_VERSIONS`, `UnsupportedVersion` variant, gate in `validate()`, new test `read_rejects_artefact_with_unknown_version`) | new test pass + 19 pre-existing artefact tests still pass | No |
| Mi13| Minor    | File `contamination_artifact.rs` (American) holds British type         | Apply            | Applied      | No         | `git mv` of `src/pop_var_caller/contamination_artifact.rs` → `src/pop_var_caller/contamination_artefact.rs`; `use` sweep in `src/pop_var_caller/mod.rs`, `src/pop_var_caller/var_calling.rs`, `src/pop_var_caller/estimate_contamination.rs` | full build + test pass | No |

## 3. Questions asked and answers

None for Wave 1. The four design questions (Q1 Wave-2 shape,
Q2 Mi14 flag rename, Q3 Wave-4 rayon policy, Q4 M5 follow-up
scope) resolved during plan review (2026-05-19) and govern later
waves. Wave 1 has no open policy questions.

## 4. Per-finding log

### M8 — `WriterConfig` lacks `#[non_exhaustive]`; add builder for `emit_gp`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The cohort CLI slice dropped `default_filter_pass`
  from `WriterConfig` without a `#[non_exhaustive]` annotation —
  an unannounced breaking change for any downstream crate using
  struct-literal construction. Fix is additive (`#[non_exhaustive]`
  + `with_emit_gp(...)` setter) and forces in-tree callers off
  struct-literal form so the next field drop / addition is a
  compile error rather than a silent break.
- **Implementation summary:**
  - Added `#[non_exhaustive]` to `pub struct WriterConfig` and a
    `pub fn with_emit_gp(self, bool) -> Self` builder setter.
  - Updated the type's module-level doc to call out the new
    construction pattern, the rationale for `#[non_exhaustive]`,
    and a runnable doctest:
    ```rust
    let cfg = WriterConfig::new(PathBuf::from("out.vcf"))
                  .with_emit_gp(true);
    ```
  - Migrated **every** in-tree struct-literal `WriterConfig {…}`
    site to `WriterConfig::new(out)[.with_emit_gp(b)]`. Sites:
    `src/pop_var_caller/var_calling.rs:397`,
    `src/pop_var_caller/var_calling_from_bam.rs:585`,
    `src/var_calling/vcf_writer/record_encode.rs:553`,
    `src/var_calling/vcf_writer/header.rs:296`,
    and four `tests/cohort_vcf_writer_integration.rs` fixtures
    (lines 161, 228, 261, 300).
- **Review suggestion used verbatim?:** No — the review sketched
  a `#[non_exhaustive]` + builder shape; this run adopted that
  shape and added the doctest the plan called for.
- **Adaptation:** None substantive.
- **Verification performed:** `cargo check --tests` clean after
  the migration; full `cargo test --all-targets` clean.
- **Files changed:** 6 files (1 in `vcf_writer` core + 2 fixtures
  in `vcf_writer` test modules + 2 orchestrators + 1 integration
  test file).
- **Tests added or modified:** Doctest added on
  `WriterConfig::new` exercising the builder pattern; 4 existing
  struct-literal `WriterConfig` test fixtures migrated.
- **Validation:**
  - `./scripts/dev.sh cargo check --tests` → 0, clean
  - `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings` → 0, clean
  - `./scripts/dev.sh cargo test --all-targets` → 0, 887 lib tests + every integration target pass
- **User input:** None.
- **Follow-up:** None. The builder pattern is the precedent the
  Wave 2 Option C builder for `PosteriorEngineConfig` will follow.
- **Residual risk:** None — `#[non_exhaustive]` on a struct only
  blocks construction by struct literal from *outside* the
  defining crate. Inside the crate, struct literals and field
  mutation continue to compile, so test fixtures and internal
  builders are unaffected; only external consumers are forced
  through the public constructor.

### Mi5 — `#[non_exhaustive]` missing on `ContaminationArtefact` + sub-structs
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The five on-disk artefact structs
  (`ContaminationArtefact`, `Provenance`, `ProvenanceInputs`,
  `BatchEntry`, `SampleEntry`) are all `pub` with all-`pub`
  fields. Without `#[non_exhaustive]`, adding a single new field
  (e.g. `provenance.git_sha`, `parameters.dust_window`) would
  silently break every downstream struct-literal constructor.
  Pairs with M8; same fix shape.
- **Implementation summary:**
  - Added `#[non_exhaustive]` to all five `pub struct`s in
    `src/pop_var_caller/contamination_artefact.rs`.
  - Updated `ContaminationArtefact`'s doc to note that
    inside-crate construction continues to work; only foreign
    crates are forced through public constructors / TOML
    deserialisation.
  - Migrated the one external struct-literal site in
    `tests/cohort_cli_integration.rs:354` to a TOML-literal +
    `fs::write` shape. This is also the **more realistic** test
    path — `run_var_calling` always loads its artefact from
    disk via `ContaminationArtefact::read`, so a TOML literal
    mirrors what a user editing the file by hand produces.
    Dropped 5 now-unused imports (`BatchEntry`,
    `ContaminationArtefact`, `Provenance`, `ProvenanceInputs`,
    `SampleEntry`) plus the `std::collections::BTreeMap` import.
- **Review suggestion used verbatim?:** No — the review
  proposed `#[non_exhaustive]` annotations only; this run also
  performed the necessary migration of the one external
  struct-literal site that the annotation would have broken.
  The integration-test refactor to TOML literal is adaptation:
  it preserves test semantics while reducing test-side coupling
  to the on-disk struct surface.
- **Adaptation:** The integration test was refactored to write
  a TOML literal rather than build the struct in code. This
  preserves the test's behaviour (the artefact arrives with
  no `NA00001` sample → `run_var_calling` must surface
  `ContaminationSampleMissing`) and removes the test-side
  dependency on the now-`#[non_exhaustive]` struct surface.
- **Verification performed:** Re-ran the full integration test
  suite after the migration; the chained
  `var_calling_rejects_contamination_artefact_missing_sample`
  case still passes.
- **Files changed:** `src/pop_var_caller/contamination_artefact.rs`,
  `tests/cohort_cli_integration.rs`.
- **Tests added or modified:** No new tests; one existing test
  migrated to TOML-literal input shape.
- **Validation:** see Wave-1 summary.
- **User input:** None.
- **Follow-up:** None — the same migration pattern (test fixture
  → TOML literal) is the natural path for any future test that
  needs to hand-craft a synthetic artefact.
- **Residual risk:** None.

### Mi6 — `provenance.version` recorded but never validated on read
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The artefact stores `provenance.version`
  (`"0.1.0"` today) but `validate()` never inspects it. A future
  `v0.2.0` format with a renamed or reshaped field would
  deserialize against the v0.1.0 reader with no error. Strict
  allow-list (`SUPPORTED_VERSIONS = &["0.1.0"]`) is two lines
  + one test and catches exactly the failure mode the field
  exists to catch.
- **Implementation summary:**
  - Added
    `pub const SUPPORTED_VERSIONS: &[&str] = &["0.1.0"];`
    next to `Q_B_SIMPLEX_TOLERANCE`. Public so callers can
    introspect supported versions.
  - Added the gate at the top of
    `ContaminationArtefact::validate`, before any per-batch or
    per-sample check — fast-fail the file before any other
    invariant work.
  - Added the new error variant:
    ```rust
    #[error(
        "contamination artefact: provenance.version `{got}` is not in the \
         supported set {supported:?}"
    )]
    UnsupportedVersion {
        got: String,
        supported: &'static [&'static str],
    },
    ```
  - Added the test `read_rejects_artefact_with_unknown_version`
    next to the existing file-level read tests. The test writes
    a structurally valid TOML with `version = "99.0.0"` and
    asserts `Err(UnsupportedVersion { got: "99.0.0", supported: SUPPORTED_VERSIONS })`.
- **Review suggestion used verbatim?:** No — the review sketched
  the variant shape; the implemented variant matches that
  sketch exactly.
- **Adaptation:** None.
- **Verification performed:** Targeted unit test
  `read_rejects_artefact_with_unknown_version` exercises the
  full file → `validate()` → error variant path and asserts on
  both fields of the variant. Existing tests
  (`round_trips_through_disk`, `validate_accepts_fixture`,
  …) confirm `"0.1.0"` continues to pass.
- **Files changed:** `src/pop_var_caller/contamination_artefact.rs`.
- **Tests added or modified:** New test
  `pop_var_caller::contamination_artefact::tests::read_rejects_artefact_with_unknown_version`.
- **Validation:**
  - `./scripts/dev.sh cargo test --lib contamination_artefact` →
    20 tests pass (was 19; new test is included).
  - Full `cargo test --all-targets` → clean.
- **User input:** None.
- **Follow-up:** None on day one. The gate is dead code today
  (only v0.1.0 has ever existed); its value lands the day a
  future producer writes anything else. The strict allow-list
  choice (rather than semver-aware leniency) was confirmed by
  the plan-review discussion as the right shape for a v0.x
  format.
- **Residual risk:** None.

### Mi13 — File `contamination_artifact.rs` (American) holds British type
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The file path used American "artifact" while
  the type, error, and every consumer used British "artefact".
  Mechanical mismatch; renaming the file (lower-blast-radius
  direction) makes the module-path/type-name spelling
  uniformly British.
- **Implementation summary:**
  - `git mv src/pop_var_caller/contamination_artifact.rs src/pop_var_caller/contamination_artefact.rs`.
  - Updated the `pub mod` declaration in
    `src/pop_var_caller/mod.rs:8` and the `pub use` re-export at
    `src/pop_var_caller/mod.rs:17`.
  - Updated `use crate::pop_var_caller::contamination_artifact::{…}`
    → `…::contamination_artefact::{…}` in
    `src/pop_var_caller/var_calling.rs` and
    `src/pop_var_caller/estimate_contamination.rs`.
  - Updated one doc-link reference at
    `src/pop_var_caller/estimate_contamination.rs:10`.
  - The Wave-1 Task 1.2 refactor of
    `tests/cohort_cli_integration.rs` (TOML-literal migration)
    incidentally removed the only test-side import of the
    artefact module, so no test-side `use` sweep was needed
    for the rename.
- **Review suggestion used verbatim?:** Yes — the suggested
  fix was a mechanical file rename + import sweep.
- **Adaptation:** None.
- **Verification performed:** `cargo check --tests` clean
  post-rename; full `cargo test --all-targets` clean. `git mv`
  preserved file history (rename is tracked).
- **Files changed:** Rename of one file + 3 source files with
  `use` line edits (mod.rs, var_calling.rs, estimate_contamination.rs).
- **Tests added or modified:** None.
- **Validation:** see Wave-1 summary.
- **User input:** None.
- **Follow-up:** None. The American spelling of "artifact" no
  longer appears in `src/` or `tests/`. Historical documents
  under `doc/devel/reports/` and the `pop_var_caller_cohort_cli.md`
  / `pop_var_caller_cohort_cli_followup.md` plan docs retain
  the old path for accuracy of the historical record.
- **Residual risk:** None.

## 5. Deferred findings to carry forward

None for Wave 1 — every finding in scope was Applied.

Waves 2 – 5 carry the rest of the deferred set; see the follow-up
plan at
`doc/devel/implementation_plans/pop_var_caller_cohort_cli_followup.md`
for the per-wave breakdown.

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
- **Outcome:** Skipped — no Apply touched perf-sensitive code
  (Wave 1 is `#[non_exhaustive]` attributes on already-existing
  public types, a builder `with_emit_gp` setter, a `validate()`
  gate that runs once at read time, and a file rename).

## 10. Commands run

- `./scripts/dev.sh cargo check --tests` (twice — once after the
  WriterConfig migration, once after the artefact-struct migration)
- `./scripts/dev.sh cargo test --lib contamination_artefact`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings`
- `./scripts/dev.sh cargo test --all-targets`
- `./scripts/dev.sh cargo doc --no-deps`
- `git -C . mv src/pop_var_caller/contamination_artifact.rs src/pop_var_caller/contamination_artefact.rs`

## 11. Command results

- `cargo check --tests` (post Task 1.1) → exit 0, `Finished dev profile`.
- `cargo check --tests` (post Task 1.2) → exit 0, `Finished dev profile`.
- `cargo test --lib contamination_artefact` (post Task 1.3) → exit 0, **20 passed**, 0 failed.
- `cargo fmt --check` → exit 0, no diff.
- `cargo clippy --workspace --all-targets -- -D warnings` → exit 0, no warnings.
- `cargo test --all-targets` → exit 0, **887** lib tests pass + integration / bench / dhat / example targets pass.
- `cargo doc --no-deps` → exit 101 with **only pre-existing** errors (2 in `pileup_to_psp.rs`, 3 `ExactMath` in `posterior_engine.rs`); zero errors introduced by Wave 1. Confirmed by stash-based bisection.
- `git mv` of the artefact file → completed; rename tracked.

## 12. Notes

- **Commit shape.** All four Wave-1 tasks landed in the
  working tree before commits. The user's plan-review decision
  was "per-task commits behind one plan per wave" (~3 commits
  per wave). On request the diff splits cleanly into:
  1. Tasks 1.1 + integration-test struct-literal migrations
     (`WriterConfig` hygiene).
  2. Task 1.2 + Task 1.3 + integration-test TOML-literal
     migration (artefact struct hygiene + schema-version gate).
  3. Task 1.4 (file rename + import sweep).
  Or a single Wave-1 commit, at the user's discretion — the
  protocol's accounting is via this report + validation, not
  via commit count.
- **Cargo doc pre-existing errors.** The previous fix-application
  report claimed the two `pileup_to_psp.rs` intra-doc-link
  errors had auto-resolved after the cohort CLI slice. The stash
  bisection in this run shows they are still present on the
  un-edited `db1ec2a` tree, so either the previous report's
  inference was incorrect or the lint surface changed in a way
  this run did not chase down. Either way, the errors are
  pre-existing relative to Wave 1 and out of scope; they should
  fold into the same Stage 1 doc-cleanup follow-up that Mi21
  tracks for `ExactMath`.
- **No `Cargo.toml` / `Cargo.lock` change.** `cargo audit` was
  not re-run for that reason.
