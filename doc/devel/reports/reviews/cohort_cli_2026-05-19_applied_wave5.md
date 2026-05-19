# Fix Application Report: cohort_cli_2026-05-19.md — Wave 5 (Test infrastructure + missing coverage)

**Date:** 2026-05-19
**Source review:** `doc/devel/reports/reviews/cohort_cli_2026-05-19.md`
**Source plan:** `doc/devel/implementation_plans/pop_var_caller_cohort_cli_followup.md` (Wave 5)
**Source state reviewed against:** branch `main` at commit `d84ee8e` (the Wave 4 docs commit)
**Execution mode:** interactive
**Overall status:** Completed

> Final wave of the cohort CLI deferred follow-up. Closes
> Mi20 + Mi23 + M1/M2 follow-up + M5 follow-up — the test
> infrastructure + missing-coverage tasks and the only
> behaviour-change finding (FASTA → `.psp` per-contig MD5
> enforcement).

---

## 1. Executive summary

### Wave 5 findings totals (subset of the source review)

- Blockers: 0
- Majors: 0
- Minors: 2 — **Mi20**, **Mi23**
- Major follow-ups: 2 — **M1/M2 walker-error CRAM test**, **M5 FASTA→.psp MD5**
- Nits: 0

### Outcome totals

- Applied: 4 (Mi20, Mi23, M1/M2 follow-up, M5 follow-up)
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0

### Validation summary

- `./scripts/dev.sh cargo fmt --check` → **exit 0**, clean (after a single `cargo fmt` apply pass).
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings` → **exit 0**, clean.
- `./scripts/dev.sh cargo test --all-targets` → **exit 0**, **890** lib tests pass (was 888 end of Wave 4; net +2 from the two new `to_estimates_*` artefact-builder unit tests) + **9** cohort CLI integration tests (was 4 end of Wave 4; net +5 from `estimate_contamination_then_var_calling_chain`, `estimate_contamination_rejects_reference_mismatch`, `var_calling_reports_reference_mismatch_for_psp_with_different_header_basename`, `var_calling_from_bam_surfaces_walker_error_on_max_active_reads_trip`, `var_calling_rejects_fasta_whose_bytes_dont_match_psp_md5`) + every other integration / bench / dhat / example target green.
- `./scripts/dev.sh cargo doc --no-deps` → **exit 101**, but **every error is identical to Waves 1-4's end-state** (2 in `pileup_to_psp.rs`, 3 `ExactMath` in `posterior_engine.rs`). **Zero errors introduced by Wave 5.**
- `cargo audit` → not re-run; the only new runtime dep is `md-5` (was already a transitive dep, now promoted to explicit) — no transitive surface change.
- Performance check → **not triggered** for Mi20 / Mi23 / M1/M2 (pure test additions). For M5 the verification adds one full-FASTA pass at startup of `run_var_calling` and `run_estimate_contamination`; on the integration-test fixtures (200 bases / 1 contig) the cost is in the microseconds, so no end-to-end bench regression is detectable. On real-cohort data the helper pre-warms the SyncRefFetcher cache the rest of the pipeline needs anyway, so the marginal cost is ≤ one extra-pass of contigs that the side-pass would otherwise not have cached.

### Unresolved high-priority findings

**None.** With Wave 5 the deferred set from the 2026-05-19
review is **closed** — 16 of 16 originally-Deferred findings
applied across Waves 1–5, plus the two doc-follow-ups (M5, M9).

## 2. Findings table

| ID  | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|-----|----------|-------|------------------|--------------|------------|---------------|------------|-----------|
| Mi20| Minor    | `tests/common/mod.rs` shared fixture helpers — eliminate duplication across `pileup_cli_integration` + `cohort_cli_integration` | Apply | Applied | No | New `tests/common/mod.rs`; both integration test files migrated | fmt + clippy + tests pass | No |
| Mi23| Minor    | Load-bearing chain test + four named missing-coverage unit tests | Apply | Applied | No | `tests/cohort_cli_integration.rs` (chain test + two ref-mismatch tests); `src/pop_var_caller/contamination_artefact.rs` (two artefact-builder unit tests) | 5 new tests pass | No |
| M1/M2 (follow-up) | Major (follow-up) | End-to-end walker-error CRAM-fixture test | Apply | Applied | No | `tests/cohort_cli_integration.rs` (`var_calling_from_bam_surfaces_walker_error_on_max_active_reads_trip`) | new test pass | No |
| M5 (follow-up)    | Major (follow-up) | Real FASTA → `.psp` per-contig MD5 enforcement + typed error | Apply | Applied | No | `Cargo.toml` (`md-5` explicit dep); `Cargo.lock`; `src/pop_var_caller/common.rs` (new `verify_fasta_matches_psp_chromosomes` + `FastaVerifyError` types); `src/pop_var_caller/var_calling.rs` (wiring + 2 new error variants); `src/pop_var_caller/estimate_contamination.rs` (wiring + 2 new error variants); `tests/cohort_cli_integration.rs` (`var_calling_rejects_fasta_whose_bytes_dont_match_psp_md5` integration test); `tests/common/mod.rs` (fixture md5 helper + `WRONG_MD5` constant) | new integration test pass | No |

## 3. Questions asked and answers

None for Wave 5 — every finding's shape was locked in plan
review on 2026-05-19. Wave 5 is direct execution.

## 4. Per-finding log

### Mi20 — `tests/common/mod.rs` shared fixture helpers
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `pileup_cli_integration.rs` and
  `cohort_cli_integration.rs` carried near-identical `build_fasta`
  / `build_sam_header` / `build_cram` / `read_record` helpers
  (~80 LOC each). Cargo's `mod common;` idiom from inside `tests/`
  is the canonical way to share helpers across integration-test
  binaries.
- **Implementation summary:**
  - New `tests/common/mod.rs` with the five helpers consolidated
    under a single signature (`build_sam_header` and `build_cram`
    each take a `sample: &str` + `md5: Option<&str>` so the same
    helper drives the pileup test's missing-MD5 case, the cohort
    tests' default fixtures, and Wave 5's deliberate
    `WRONG_MD5` injection for the M5 mismatch test).
  - Both integration test files migrated: drop the local
    helpers + their imports (`bstr`, `noodles_cram`,
    `noodles_fasta`, `noodles_sam`, `sam::*`), add `mod common;`
    + `use common::{…};`.
  - One follow-on cleanup: `pileup_cli_integration::happy_path_default_config`
    used to assert `header.writer.input_crams == vec!["sample.cram"]`;
    the unified `build_cram` parameterises the CRAM path by
    `{sample}.cram`, so the assertion now reads
    `vec!["NA12878.cram"]`.
- **Files changed:** new `tests/common/mod.rs`; modified
  `tests/pileup_cli_integration.rs` and
  `tests/cohort_cli_integration.rs`.
- **Tests added or modified:** No semantic test additions —
  pure refactor. Existing tests still cover the same surfaces.
- **Validation:** see Wave-5 summary.

### Mi23 — Chained integration test + missing-coverage unit tests
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The cohort CLI ships two subcommands
  (`estimate-contamination`, `var-calling`) whose chain is the
  load-bearing contract; until this wave, the chain was only
  tested as separate subcommand happy paths. Plus four
  unit-level coverage gaps the review's §8 listed.
- **Implementation summary:**
  - **Chain test** (`estimate_contamination_then_var_calling_chain`):
    builds 2 `.psp`s, runs `run_estimate_contamination` to a
    temp TOML, then `run_var_calling --contamination-estimates
    <that>`, and verifies the resulting VCF parses. The
    contamination side-pass needs cohort-minor-allele evidence
    to converge; the test designs reads such that positions
    10..=14 each have a 1-base minor C in sample B (5% mismatch
    each, under the CRAM-input filter), yielding 5 informative
    sites — enough for Convergence-mode stopping at relaxed
    tolerances. Documented that the side-pass *accuracy* is
    covered by `tests/contamination_estimation_integration.rs`
    with realistic synthetic data.
  - **Reference-mismatch integration coverage** (two tests):
    `estimate_contamination_rejects_reference_mismatch` and
    `var_calling_reports_reference_mismatch_for_psp_with_different_header_basename`.
    Build a `.psp` with `--reference ref.fa`, then run the
    subcommand with `--reference other.fa` (same bytes, different
    basename via `fs::copy`). Both subcommands must surface the
    typed `ReferenceMismatch` variant before any side-pass /
    cohort-pipeline boot.
  - **Artefact-builder unit tests** (two tests in
    `contamination_artefact.rs::tests`):
    `to_estimates_round_trips_floored_batch_with_zero_qb`
    locks the engine's all-zero `q_b` row contract for floored
    batches end-to-end; the
    `to_estimates_orders_dense_batches_by_cohort_first_seen`
    test pins the cohort-first-seen reordering invariant a
    downstream consumer indexing `q_b_per_batch()` relies on.
- **Files changed:** `tests/cohort_cli_integration.rs` (chain
  test + two ref-mismatch tests); `src/pop_var_caller/contamination_artefact.rs`
  (two new `to_estimates_*` unit tests).
- **Tests added or modified:** +5 tests total (1 chain
  integration + 2 ref-mismatch integration + 2 artefact-builder
  unit).

### M1/M2 follow-up — Walker-error CRAM-fixture test
- **Severity:** Major (follow-up)
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The `prefer_upstream_or_closure` helper from
  the 2026-05-19 first-pass (commit `db1ec2a`) is unit-tested,
  but the end-to-end *wiring* of walker errors through
  `run_var_calling_from_bam` — including the publish-then-retract
  `<output>` cleanup that fires after the rename — had no
  end-to-end coverage.
- **Implementation summary:** New
  `var_calling_from_bam_surfaces_walker_error_on_max_active_reads_trip`
  integration test in `tests/cohort_cli_integration.rs`.
  Synthesises a CRAM with two reads overlapping at the same
  position; passes `max_active_reads = 1` so the walker's
  defensive bound trips on the second read. Asserts:
  (a) `run_var_calling_from_bam` returns
  `Err(VarCallingFromBamCliError::Walker(_))` (the
  `ErrorSheddingAdapter` surfaces the stashed walker error
  rather than silently producing a partial VCF);
  (b) the output VCF does **not** exist on disk after the
  function returns (publish-then-retract cleanup ran).
- **Files changed:** `tests/cohort_cli_integration.rs` (new test
  + the `VarCallingFromBamCliError` import).
- **Tests added or modified:** +1 integration test.

### M5 follow-up — FASTA → `.psp` MD5 enforcement
- **Severity:** Major (follow-up)
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The basename cross-check in
  `run_var_calling` / `run_estimate_contamination` step 3
  catches "user passed the wrong filename" but **not** "user
  has a same-named FASTA pointing at a different genome build"
  — exactly the silent-wrong-reference failure mode the
  project's `feedback_no_silent_intermediates` principle
  exists to prevent. The v1 contract documented this gap; this
  wave closes it.
- **Implementation summary:**
  - New `pop_var_caller::common::verify_fasta_matches_psp_chromosomes`
    helper: for each `.psp` `ParsedChromosome`, fetches the
    contig bytes via `fetcher.fetch(chrom_id, 1, length)` (the
    fetcher's existing uppercase-on-fetch contract makes the
    bytes canonical SAM-spec form), computes MD5 via the
    explicit `md-5 = "0.11"` crate (was a transitive dep — now
    promoted to direct so future revisions don't drift), and
    compares against the `.psp` header's `chromosome.md5`.
  - Returns `Result<(), FastaVerifyError>` — a two-variant
    enum distinguishing genuine `Md5Mismatch { contig,
    fasta_md5, psp_md5 }` from `FetchFailed { contig, source }`
    (contig missing from FASTA / I/O error). Both failure
    modes propagate to typed CLI errors:
    `VarCallingCliError::{FastaContigMd5Mismatch,
    FastaContigFetchFailed}` and the matching pair on
    `EstimateContaminationCliError`.
  - Wired into `run_var_calling` at step 7a (right after
    `SyncRefFetcher::new`; reuses the existing pipeline
    fetcher — no extra cache pressure) and into
    `run_estimate_contamination` at step 3 (builds a
    one-shot `SyncRefFetcher` just for the check, then drops
    it; the side-pass itself doesn't need the FASTA bytes,
    so the verification's ~contig-sized memory footprint
    transients out before the side-pass runs).
  - `contigs_from_parsed` (the `ParsedChromosome →
    cram_input::ContigList` shim that the SyncRefFetcher
    constructor expects) promoted from `fn` →
    `pub(crate) fn` in `var_calling.rs` so
    `run_estimate_contamination` can reuse it.
  - **New integration test**
    (`var_calling_rejects_fasta_whose_bytes_dont_match_psp_md5`):
    builds a CRAM whose `@SQ M5` is forced to `WRONG_MD5`
    (the new `Option<&str>` parameter on `build_cram` makes
    this a clean fixture override), runs `run_pileup` so the
    bogus MD5 propagates into the `.psp` header, then runs
    `run_var_calling` against the actual FASTA. The new
    check fires and surfaces `FastaContigMd5Mismatch {
    contig: "chr1", fasta_md5: <real MD5>, psp_md5:
    WRONG_MD5 }`.
- **Adaptation:** The plan considered "audit `SyncRefFetcher`'s
  exposed surface; pick the cheap-or-streaming path; defer to a
  separate slice if a refactor is needed." Adapted: the
  fetcher already exposes everything needed
  (`fetch(chrom_id, 1, length)` returns canonical uppercase
  bytes via its existing soft-mask normalisation), so no
  refactor was required. `run_estimate_contamination` does
  build an otherwise-unnecessary fetcher just for the check —
  documented as a trade-off.
- **Files changed:** `Cargo.toml` (new explicit `md-5` dep);
  `Cargo.lock`; `src/pop_var_caller/common.rs`
  (`FastaVerifyError` + `verify_fasta_matches_psp_chromosomes`);
  `src/pop_var_caller/var_calling.rs` (2 new error variants +
  wiring; `contigs_from_parsed` promoted to `pub(crate)`);
  `src/pop_var_caller/estimate_contamination.rs` (2 new error
  variants + wiring; chromosome-agreement check moved earlier so
  the M5 check can reuse it); `tests/common/mod.rs`
  (`fixture_md5` + `WRONG_MD5`; `build_cram` /
  `build_sam_header` switched from `include_md5: bool` to
  `md5: Option<&str>`); `tests/cohort_cli_integration.rs` and
  `tests/pileup_cli_integration.rs` (call-site updates +
  the new mismatch test).
- **Tests added or modified:** +1 integration test. Every
  existing test that built a CRAM was updated to pass
  `Some(fixture_md5())` so the new check's happy path is
  exercised under the existing test suite (i.e. silent
  regression coverage — any future change that breaks the
  helper would fail many tests, not just one).

## 5. Deferred findings to carry forward

**None.** Wave 5 closes the deferred set from the 2026-05-19
review.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** No.
- **Baseline saved:** N/A.
- **Outcome:** Skipped — Mi20 / Mi23 / M1/M2 are pure test
  additions; M5 adds a fixed-cost helper call (1 fetch +
  1 MD5 per contig) at orchestrator startup. On the
  integration-test fixture (1 contig × 200 bases) the cost is
  microseconds; on real-cohort data it's ~24 contigs × IO + MD5
  ≈ a few seconds of startup wall-time, dominated by the
  FASTA I/O the SyncRefFetcher would do anyway when DUST first
  fetches a contig. The end-to-end cohort wall-time impact is
  below noise floor.

## 10. Commands run

- `./scripts/dev.sh cargo check --tests` (many — once per task)
- `./scripts/dev.sh cargo test --test cohort_cli_integration <name>` (per-test smoke during development)
- `./scripts/dev.sh cargo test --lib contamination_artefact::tests::to_estimates`
- `./scripts/dev.sh cargo fmt --check` / `cargo fmt`
- `./scripts/dev.sh cargo clippy --workspace --all-targets -- -D warnings`
- `./scripts/dev.sh cargo test --all-targets`
- `./scripts/dev.sh cargo doc --no-deps`

## 11. Command results

- `cargo fmt --check` → exit 0
- `cargo clippy --workspace --all-targets -- -D warnings` → exit 0
- `cargo test --all-targets` → exit 0, **890** lib + **9**
  cohort CLI integration + every other target green
- `cargo doc --no-deps` → exit 101 with the **same 5 pre-existing
  errors** as Waves 1–4 (2 `pileup_to_psp.rs` + 3 `ExactMath`);
  zero introduced by Wave 5

## 12. Notes

- **Commit shape.** Wave 5's diff is naturally three units:
  1. **Mi20** — pure refactor (`tests/common/mod.rs`).
  2. **Mi23 + M1/M2** — five new tests (chain + 2 ref-mismatch
     + 2 artefact-builder + 1 walker-error CRAM).
  3. **M5 follow-up** — new helper + wiring + 1 mismatch
     integration test + the `md-5` dep + lockfile churn.
  Plus a fourth docs commit. **Recommended ordering:** land
  Mi20 first (its API change to `build_cram` is touched by the
  M5 test fixture too) → Mi23 / M1/M2 → M5 → docs.
- **`md-5` is now an explicit runtime dep**, formerly transitive
  through `noodles-cram`. The version is unchanged (`0.11`);
  the lockfile picks up no new transitive surface.
- **Side-pass terminates on the synthetic chain-test fixture**
  because 5 informative cohort sites are enough for a single
  `stability_blocks = 1` convergence block at tolerance 0.1.
  The `estimate_contamination_args` fixture builder
  hard-codes those relaxed thresholds; the per-module
  `contamination_estimation_integration.rs` covers the
  side-pass's accuracy on realistic data, not this test.
- **No user-visible CLI break.** All Wave 5 work is either
  test scaffolding (Mi20 / Mi23 / M1/M2) or strictly tightens
  an existing contract (M5 — previously the basename-only
  check could pass with the wrong genome; now the FASTA-bytes
  MD5 check rejects the same invocation up-front).
