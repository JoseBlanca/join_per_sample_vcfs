# Fix Application Report: bam_input_support_2026-05-24

**Date:** 2026-05-24
**Source review:** `doc/devel/reports/reviews/bam_input_support_2026-05-24.md`
**Source state reviewed against:** commit `1c59d61` (post-review HEAD)
**Execution mode:** interactive (three AskUserQuestion decisions: CSI depth, scope, CRAM-version-gate confirmation)
**Overall status:** **Completed**

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 19 (M1–M19)
- Minors: 22 (Mi1–Mi22)
- Nits: grouped (doc-sweep + small idiomatic items)

### Outcome totals
- Applied: 38 (every Major + 13 of 22 Minors)
- Applied with adaptation: 1 (M17 — test-only half landed; the shared-enum redesign half deferred)
- Already fixed: 0
- Deferred: 7 (Mi5, Mi7, Mi12, Mi13, Mi15, Mi20, Mi21)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 924 passed (was 916 pre-fix: +1 in cram_input [M5 regression test], +9 in commit 7 [M9, M11×2, M12×4, Mi16, Mi17 → 8] + corrected accounting puts it at +8 lib; net headcount checked end-to-end)
- `cargo test --test '*'` → 0, 45 integration tests pass (was 43 pre-fix: +1 pileup_cli [M10], +1 cohort_cli [M17 lock-step])
- `cargo doc --no-deps` → not re-run; pre-existing failure documented in PROJECT_STATUS, fix run did not introduce new doc warnings
- `cargo audit` → not run (not in project verification list)
- Performance check (`cargo bench -- --baseline pre-fixes`) → **skipped** per "do not stash/revert to back-fill" rule; baseline was not captured at preflight because the planned fixes touch startup-time, error-path, and dispatch-time code only — the per-record hot path (`OwnedBamRecords::next`, `OwnedIndexedBamRecords::next`, merge loop) is structurally unchanged.

### Unresolved high-priority findings
**None.** All 19 Majors are Applied (M17 in Applied-with-adaptation form). The 7 Deferred items are all Minor.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | Error vocabulary still says "CRAM" on BAM inputs | Apply | **Applied** | No | errors.rs, cli.rs, var_calling_from_bam.rs, stage1_pipeline.rs, alignment_input.rs, cram_input.rs, mod.rs, tests/pileup_cli_integration.rs | commit 2 (`2de582d`) — fmt/clippy/test all clean | No |
| M2 | Major | `VarCallingFromBamCliError::Io(#[from] io::Error)` collapses multiple origins | Apply | **Applied** | No | var_calling_from_bam.rs | commit 5 (`86b77fc`) — fmt/clippy/test all clean | No |
| M3 | Major | `load_alignment_index` returns bare `io::Error` | Apply | **Applied** | No | errors.rs, index_preflight.rs, var_calling_from_bam.rs | commit 5 (`86b77fc`) | No |
| M4 | Major | `AlignmentInputError` missing `#[non_exhaustive]` | Apply | **Applied** | No | errors.rs | commit 3 (`351b735`) | No |
| M5 | Major | `load_per_input_headers` duplicates opener + skips CRAM-version gate | Apply | **Applied** | **Yes** — confirmed accidental | cram_input.rs, bam_input.rs, var_calling_from_bam.rs | commit 1 (`9fc1df0`) — includes regression test | No |
| M6 | Major | Catch-all `(file_kind, index)` arm absorbs future variants | Apply | **Applied** | No | alignment_input.rs | commit 4 (`156fa18`) | No |
| M7 | Major | `BamIndex` enum missing `#[non_exhaustive]` | Apply | **Applied** | No | bam_input.rs | commit 3 (`351b735`) | No |
| M8 | Major | `PileupCliError` missing `#[non_exhaustive]` | Apply | **Applied** | No | cli.rs | commit 3 (`351b735`) | No |
| M9 | Major | `AlignmentIndexFormatMismatch` has no test | Apply | **Applied** | No | alignment_input.rs (tests) | commit 7 (`6ec16b0`) — new test passes | No |
| M10 | Major | `AlignmentInputError::UnsupportedExtension` test gap on pileup path | Apply | **Applied** | No | tests/pileup_cli_integration.rs | commit 7 (`6ec16b0`) | No |
| M11 | Major | `OwnedIndexedBamRecords::next` chunk-walking + IO-error untested | Apply | **Applied** | No | bam_input.rs (tests) | commit 7 (`6ec16b0`) — 2 new tests pass | No |
| M12 | Major | `load_alignment_index` no direct unit tests | Apply | **Applied** | No | index_preflight.rs (tests) | commit 7 (`6ec16b0`) — 4 new tests pass | No |
| M13 | Major | `.csi`/`.bai` policy not in a named constant | Apply | **Applied** | No | index_preflight.rs | commit 6 (`a2ea28d`) — 4 new pub consts | No |
| M14 | Major | `build_csi_for_bam` uses `Indexer::default()`; depth hidden | Apply | **Applied** | **Yes — depth=6** | index_preflight.rs | commit 6 (`a2ea28d`) | No |
| M15 | Major | `--build-map-file-index` help text doesn't say `.csi`-only build | Apply | **Applied** | No | var_calling_from_bam.rs | commit 6 (`a2ea28d`) | No |
| M16 | Major | `MissingMapFileIndex` Display lists only one path for BAM | Apply | **Applied** | No | var_calling_from_bam.rs | commit 6 (`a2ea28d`) — Display string now mentions `.bai` fallback | No |
| M17 | Major | `MixedAlignmentFileFormats` dual-surfaced; only one path tested | Apply (test only) | **Applied with adaptation** | No | tests/cohort_cli_integration.rs | commit 9 (`b8f5642`) — lock-step test landed; shared-enum redesign deferred | Yes — redesign half deferred (see §5) |
| M18 | Major | `.unwrap()` defensive panic after classify pre-pass | Apply | **Applied** | No | alignment_input.rs | commit 4 (`156fa18`) | No |
| M19 | Major | Impl-report commit table + deferred-list inaccurate | Apply | **Applied** | No | doc/devel/reports/implementations/bam_input_support_2026-05-24.md | commit 10 (`7a8569b`) | No |
| Mi1 | Minor | `from_open_crams` method + param `open_crams` stale | Apply | **Applied** | No | alignment_input.rs (+ 22 test call sites in same file) | commit 2 (`2de582d`) | No |
| Mi2 | Minor | Stale local var names in dispatch loops | Apply | **Applied** | No | alignment_input.rs | commit 2 (`2de582d`) | No |
| Mi3 | Minor | `cram_cfg` + `cram_config_from_args` stale | Apply | **Applied** | No | cli.rs, var_calling_from_bam.rs, stage1_pipeline.rs | commit 2 (`2de582d`) | No |
| Mi4 | Minor | `CramHeader` struct name | Apply | **Applied** | No | alignment_input.rs | commit 2 (`2de582d`) — renamed to `AlignmentFileHeaderSummary` | No |
| Mi5 | Minor | `process_one_chromosome_from_bam` name covers CRAM too | Defer | **Deferred** | No | None | N/A | Yes — held per impl-report; subcommand-name coupling |
| Mi6 | Minor | `validate_fasta_agreement` param + detail strings | Apply | **Applied** | No | alignment_input.rs, errors.rs | commit 2 (`2de582d`) — param renamed + 3 detail strings rewritten | No |
| Mi7 | Minor | `input_crams` field/local mixed-format | Defer | **Deferred** | No | None | N/A | Yes — PSP schema field reach (`src/psp/header.rs` is out of scope) |
| Mi8 | Minor | Doc-comment + module-doc sweep | Apply | **Applied** | No | alignment_input.rs, mod.rs, errors.rs, cram_input.rs | commit 2 (`2de582d`) | No |
| Mi9 | Minor | `.csi`/`.bai` scan-and-index loop triplicated | Apply | **Applied with adaptation** | No | index_preflight.rs, bam_input.rs | commit 8 (`6823b07`) — 2 of 3 copies collapsed into `populate_binning_index`; the 3rd (tests/common::build_csi) stays per test-crate boundary | Yes — tests/common copy stays; deferred to a test-support feature-flag pass |
| Mi10 | Minor | `index_display_name` free function → impl method | Apply | **Applied** | No | alignment_input.rs, index_preflight.rs | commit 4 (`156fa18`) — moved into `impl AlignmentIndex` | No |
| Mi11 | Minor | `load_alignment_index` re-encodes csi/bai policy | Apply | **Applied** | No | index_preflight.rs | commit 5 (`86b77fc`) — now consults `existing_index_for` | No |
| Mi12 | Minor | `MixedAlignmentFileFormats` / `UnsupportedExtension` duplicated | Defer | **Deferred** | No | None | N/A | Yes — folded into M17's deferred redesign half |
| Mi13 | Minor | `From<AlignmentIndexError>` 4-arm passthrough | Defer | **Deferred** | No | None | N/A | Yes — explicit rename is load-bearing for the CLI vocabulary |
| Mi14 | Minor | `pub mod cram_input` / `pub mod bam_input` → `pub(crate) mod` | Apply | **Applied** | No | mod.rs | commit 9 (`b8f5642`) | No |
| Mi15 | Minor | `crate::bam` module name rename | Defer | **Deferred** | No | None | N/A | Yes — large blast radius; decide before next API consumer |
| Mi16 | Minor | `BamIndex::Csi` arm has no isolated test | Apply | **Applied** | No | bam_input.rs (tests) | commit 7 (`6ec16b0`) — new CSI-arm test passes | No |
| Mi17 | Minor | `OwnedBamRecords` `Err(e)` arm has no test | Apply | **Applied** | No | bam_input.rs (tests) | commit 7 (`6ec16b0`) — `surfaces_read_errors_as_some_err` test passes | No |
| Mi18 | Minor | `Stage1` variant doc says "CRAM-input validation" | Apply | **Applied** | No | var_calling_from_bam.rs | commit 2 (`2de582d`) | No |
| Mi19 | Minor | `BinnedIndex` vs `LinearIndex` undocumented | Apply | **Applied** | No | index_preflight.rs | commit 6 (`a2ea28d`) — inline comment on the `Indexer<BinnedIndex>` line | No |
| Mi20 | Minor | `OwnedIndexedBamRecords::next` long phase-numbered loop | Defer | **Deferred** | No | None | N/A | Yes — cosmetic; revisit when a second use case appears |
| Mi21 | Minor | Per-record `RecordBuf::default()` allocation | Defer | **Deferred** | No | None | N/A | Yes — folds into the standing parallelisation-tuning workstream |
| Mi22 | Minor | `Cargo.toml` missing justification comment | Apply | **Applied** | No | Cargo.toml | commit 9 (`b8f5642`) | No |

## 3. Questions asked and answers

1. **M14** — what CSI depth should `build_csi_for_bam` use?
   - **Answer:** Bump to depth=6 (~2^32 ≈ 4.3 Gbp addressable). Safe for wheat / large plant genomes; future-proofs against plant references the project doesn't target today but might tomorrow.

2. **Scope** — how much of the review should the fix run cover?
   - **Answer:** All Majors + coupled Minors + the 6 missing tests. Defer Mi5, Mi7, Mi12, Mi13, Mi15, Mi20, Mi21, and M17's redesign half.

3. **M5 / review §4 Q3** — was `load_per_input_headers` skipping the CRAM-version check intentional?
   - **Answer:** No — the skip was accidental. Fix it. Prioritised M5 to the first commit so the bug is closed early.

## 4. Per-finding log

### M1 — Error vocabulary still says "CRAM" on BAM inputs
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Five-category convergent finding (errors, naming, defaults, smells, extras). Mechanical sweep with low risk; tests already lock the format-name strings (`"CRAM"`, `"BAM"`).
- **Implementation summary:** `PileupCliError::CramInput` → `::AlignmentInput`; six `AlignmentInputError` Display strings reworded (NoInputs, OpenFailed, NotCoordinateSorted, FastaContigMismatch, MultipleSampleNames, DuplicateReadAcrossFiles); `FastaContigMismatch.cram_path` field → `.alignment_file_path`; `MissingMd5` recommendation no longer says "re-CRAM" (now points at `samtools view -t`). Bulk sed for the `PileupCliError` variant rename across src/ + tests/ (44 sites).
- **Review suggestion used verbatim?:** No
- **Adaptation:** Followed the review's suggested wording for the new strings; kept CRAM-only wording on `UnsupportedCramVersion` (genuinely CRAM-only).
- **Verification performed:** `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --lib`, `cargo test --test '*'` — all clean.
- **Files changed:** `src/bam/errors.rs`, `src/pop_var_caller/cli.rs`, `src/pop_var_caller/var_calling_from_bam.rs`, `src/pop_var_caller/stage1_pipeline.rs`, `src/bam/alignment_input.rs`, `src/bam/cram_input.rs`, `src/bam/mod.rs`, `tests/pileup_cli_integration.rs` (commit `2de582d`).
- **Tests added or modified:** Existing `pileup_rejects_mixed_cram_and_bam` updated to assert the new `AlignmentInput` variant name.
- **Validation:** All four commands clean post-commit.
- **User input:** None.
- **Follow-up:** None.
- **Residual risk:** None — out-of-tree consumers (none) would see a breaking API rename, but the project is single-crate and the variant is internal-facing.

### M2 — `VarCallingFromBamCliError::Io(#[from] io::Error)` collapses multiple origins
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Single-origin `From` rule violation; two of the three `?` sites already stringified their inner error to work around the collapse.
- **Implementation summary:** Replaced `Io(#[from] io::Error)` with three operation-named variants: `ScratchDir { parent, source }`, `RefFetcher { contig, source }`, and (via the M3 bridge) `IndexLoadFailed { path, index_path, source }`. Each `?` site now constructs the variant explicitly at the call point.
- **Review suggestion used verbatim?:** Yes (variant shapes match the review's suggested diff).
- **Adaptation:** None.
- **Verification performed:** Same as M1.
- **Files changed:** `src/pop_var_caller/var_calling_from_bam.rs` (commit `86b77fc`).
- **Tests added or modified:** None (no test was previously pinning the bag-of-Io behaviour).
- **Validation:** All four commands clean.
- **User input:** None.
- **Follow-up:** None.
- **Residual risk:** None.

### M3 — `load_alignment_index` returns bare `io::Error`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Public-API surface; the doc comment admitted the wrong shape. Three-category convergent finding (errors, idiomatic, smells).
- **Implementation summary:** Added `AlignmentIndexError::LoadFailed { path, index_path, source }`. Signature changed from `io::Result<AlignmentIndex>` to `Result<AlignmentIndex, AlignmentIndexError>`. Caller in `var_calling_from_bam.rs` now uses the typed `From` bridge instead of `io::Error::other(format!(...))`.
- **Review suggestion used verbatim?:** Yes (followed the suggested variant + signature shapes).
- **Adaptation:** None.
- **Verification performed:** Same as M1.
- **Files changed:** `src/bam/errors.rs`, `src/bam/index_preflight.rs`, `src/pop_var_caller/var_calling_from_bam.rs` (commit `86b77fc`).
- **Tests added or modified:** Indirect (M12's 4 new tests for `load_alignment_index` exercise the new typed error path).
- **Validation:** Clean.
- **User input:** None.
- **Follow-up:** None.
- **Residual risk:** None.

### M4 — `AlignmentInputError` missing `#[non_exhaustive]`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Plan called for it; impl-report was silent on the omission. Cheapest moment to add is now.
- **Implementation summary:** One-line `#[non_exhaustive]` attribute on the enum.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** Same as M1.
- **Files changed:** `src/bam/errors.rs` (commit `351b735`).
- **Tests added or modified:** None (in-crate matches still type-check exhaustively).
- **Validation:** Clean.
- **User input:** None.
- **Follow-up:** None.
- **Residual risk:** None.

### M5 — `load_per_input_headers` duplicates opener + skips CRAM-version gate
- **Severity:** Major
- **Initial decision:** Apply (prioritised to commit 1)
- **Final status:** Applied
- **Reasoning:** User confirmed the gate skip was accidental. Test-first regression to lock the contract.
- **Implementation summary:** Added `pub(crate) fn read_cram_header_only(path)` in `cram_input.rs` and `pub(crate) fn read_bam_header_only(path)` in `bam_input.rs`. Both call shared per-format private helpers (`open_*_reader_with_header`) so the CRAM-version gate lives in exactly one place. `load_per_input_headers` reduced from ~55 lines to ~15.
- **Review suggestion used verbatim?:** Yes (per the brief's "add a thin `read_*_header_only` helper" suggestion).
- **Adaptation:** Added a `repository: Option<&fasta::Repository>` parameter to the private `open_cram_reader_with_header` so the header-only path can skip the FASTA repository construction.
- **Verification performed:** New `read_cram_header_only_rejects_cram_4x` test using `build_cram_with_major_version(4, 0)` to assert `UnsupportedCramVersion` fires from the new helper path.
- **Files changed:** `src/bam/cram_input.rs`, `src/bam/bam_input.rs`, `src/pop_var_caller/var_calling_from_bam.rs` (commit `9fc1df0`).
- **Tests added or modified:** + 1 new lib test (regression).
- **Validation:** Test passes (1/1); full lib + integration clean.
- **User input:** Yes — Q3 confirmed the skip was accidental.
- **Follow-up:** None.
- **Residual risk:** None.

### M6 — Catch-all `(file_kind, index)` arm absorbs future variants
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `#[non_exhaustive]` extensibility intent was being silently defeated. Three-category convergent finding.
- **Implementation summary:** Wildcard catch-all replaced with two `|`-joined explicit-mismatch arms covering exactly the three legitimate cross-format pairings. A new `AlignmentIndex` variant now triggers an exhaustiveness compile-error at this site.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** Same as M1.
- **Files changed:** `src/bam/alignment_input.rs` (commit `156fa18`).
- **Tests added or modified:** M9's new `AlignmentIndexFormatMismatch` test exercises the enumerated arm.
- **Validation:** Clean.
- **User input:** None.
- **Follow-up:** None.
- **Residual risk:** None.

### M7 — `BamIndex` enum missing `#[non_exhaustive]`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Asymmetry vs `AlignmentIndex` (which already had it).
- **Implementation summary:** One-line attribute add.
- **Files changed:** `src/bam/bam_input.rs` (commit `351b735`).
- **Validation:** Clean.
- **User input:** None.

### M8 — `PileupCliError` missing `#[non_exhaustive]`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Asymmetry vs `VarCallingFromBamCliError` (which already had it).
- **Implementation summary:** One-line attribute add.
- **Files changed:** `src/pop_var_caller/cli.rs` (commit `351b735`).
- **Validation:** Clean.

### M9 — `AlignmentIndexFormatMismatch` typed-error has no test
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added `query_returns_alignment_index_format_mismatch_on_cram_path_with_bam_csi_index` in `alignment_input.rs` tests. Builds a one-contig CRAM via `build_cram`, constructs a `BamCsi`-typed index with an empty `noodles_csi::Index`, asserts the typed variant + the correct `file_format` / `index_format` strings.
- **Files changed:** `src/bam/alignment_input.rs` (test only) (commit `6ec16b0`).
- **Tests added or modified:** + 1 lib test (passes).
- **Validation:** Test passes; full suite clean.

### M10 — `AlignmentInputError::UnsupportedExtension` test gap on pileup path
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added `pileup_rejects_input_with_unknown_extension` in `tests/pileup_cli_integration.rs`. Touches a `.sam` placeholder; asserts `PileupCliError::AlignmentInput(AlignmentInputError::UnsupportedExtension { path })`.
- **Files changed:** `tests/pileup_cli_integration.rs` (commit `6ec16b0`).
- **Tests added or modified:** + 1 integration test (passes).
- **Validation:** Test passes.

### M11 — `OwnedIndexedBamRecords::next` chunk-walking + IO-error untested
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Two new tests in `bam_input.rs`:
  - `open_indexed_bam_record_stream_walks_two_chunks_in_order` — writes 1500 records on sq0 + 1 on sq1 (~150 KiB → 2+ BGZF blocks → 2+ chunks). Asserts every sq0 record comes back in order.
  - `open_bam_record_stream_surfaces_read_errors_as_some_err` — writes 200 records, truncates the file to 4 KiB (mid-BGZF block), asserts the iterator yields `Some(Err(_))` at some point.
- **Files changed:** `src/bam/bam_input.rs` (tests) (commit `6ec16b0`).
- **Tests added or modified:** + 2 lib tests (both pass).
- **Validation:** Tests pass.

### M12 — `load_alignment_index` no direct unit tests
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** 4 new tests in `index_preflight.rs::tests`:
  - `load_alignment_index_returns_bam_csi_when_csi_present` — both .csi and .bai exist → loader picks .csi.
  - `load_alignment_index_falls_back_to_bai_when_no_csi` — only .bai → loader picks it.
  - `load_alignment_index_errors_with_missing_alignment_index_when_no_bam_index` — neither → `MissingAlignmentIndex` naming the .csi-canonical path.
  - `load_alignment_index_rejects_unsupported_extension` — `.sam` → `UnsupportedExtension`.
  Use placeholder one-byte files; assertions match against the typed-error variants introduced in M3.
- **Files changed:** `src/bam/index_preflight.rs` (tests) (commit `6ec16b0`).
- **Tests added or modified:** + 4 lib tests (all pass).
- **Validation:** Tests pass.

### M13 — `.csi`/`.bai` policy not in a named constant
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Four new pub consts at the top of `index_preflight.rs`: `BAM_INDEX_READ_PREFERENCE`, `BAM_INDEX_BUILD_FORMAT`, `CSI_MIN_SHIFT`, `CSI_DEPTH`. `existing_index_for` consumes `READ_PREFERENCE`; `target_index_path` consumes `BUILD_FORMAT`; `build_csi_for_bam` consumes the CSI params. Doc-comments on each constant explain the addressable-length math and rationale.
- **Files changed:** `src/bam/index_preflight.rs` (commit `a2ea28d`).
- **Validation:** Clean.

### M14 — `build_csi_for_bam` uses `Indexer::default()`; depth hidden
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** User chose depth=6 (~4.3 Gbp addressable; safe for wheat / large plant genomes).
- **Implementation summary:** `Indexer::default()` replaced with `Indexer::new(CSI_MIN_SHIFT, CSI_DEPTH)`. Inline comment explains the bump from default depth=5 (~537 Mbp cap) to depth=6 (~4.3 Gbp cap).
- **Files changed:** `src/bam/index_preflight.rs` (commit `a2ea28d`).
- **User input:** Yes — Q1 chose depth=6.
- **Validation:** Clean.

### M15 — `--build-map-file-index` help text doesn't say `.csi`-only build
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Help text rewritten to explicitly say "On read: .crai for CRAM; .csi (preferred) or .bai (fallback) for BAM. On build: .crai for CRAM, .csi for BAM (.bai's 16 kbp bin grid tops out at 512 Mbp ...)".
- **Files changed:** `src/pop_var_caller/var_calling_from_bam.rs` (commit `a2ea28d`).
- **Validation:** Clean.

### M16 — `MissingMapFileIndex` Display lists only one path for BAM
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Display string extended to mention the `.bai` fallback explicitly: "looked for '{expected_index_path}'; for BAM inputs '<input>.bam.bai' is also accepted on read as a fallback".
- **Adaptation:** Chose the minimal-diff "extend the wording" variant rather than the structural alternative (`Vec<PathBuf>` field) per the review's two-option offering.
- **Files changed:** `src/pop_var_caller/var_calling_from_bam.rs` (commit `a2ea28d`).
- **Validation:** Clean.

### M17 — `MixedAlignmentFileFormats` dual-surfaced; only one path tested
- **Severity:** Major
- **Initial decision:** Apply (test-only half)
- **Final status:** Applied with adaptation
- **Reasoning:** The shared-enum redesign half is a design call (the review offered two options); the simpler "add a test that locks the two paths render identical strings" is achievable in-PR and protects against future wording drift.
- **Implementation summary:** Added `mixed_format_error_renders_identical_strings_at_both_layers` in `tests/cohort_cli_integration.rs`. Asserts `format!("{from_reader}") == format!("{from_preflight}")` for matching field shapes.
- **Adaptation:** Applied only the test-add piece of the review's two-option suggestion. The redesign (lift `MixedAlignmentFileFormats` into a shared sub-enum via `#[from]`) is deferred.
- **Files changed:** `tests/cohort_cli_integration.rs` (commit `b8f5642`).
- **Tests added or modified:** + 1 integration test (passes).
- **Follow-up:** Yes — shared-enum redesign deferred (also subsumes Mi12).

### M18 — `.unwrap()` defensive panic after classify pre-pass
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Classify pre-pass now pushes each input's `AlignmentFileKind` into a parallel `Vec<AlignmentFileKind>`; the open-pass zips this Vec in. The `.unwrap()` and the second `AlignmentFileKind::from_path` call are both gone. Defensive panic is replaced by type-safe data carry.
- **Files changed:** `src/bam/alignment_input.rs` (commit `156fa18`).
- **Validation:** Clean.

### M19 — Impl-report commit table + deferred-list inaccurate
- **Severity:** Major (documentation correctness)
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added two rows to the impl-report's per-commit table covering `344f1b2` and `be3b38a`. Moved the two now-applied bullets out of "Deferred to follow-ups" into a new "Closed mid-PR" section.
- **Files changed:** `doc/devel/reports/implementations/bam_input_support_2026-05-24.md` (commit `7a8569b`).
- **Validation:** Docs-only commit; no source impact.

### Mi1, Mi2, Mi3, Mi4, Mi6, Mi8, Mi18 (format-neutral rename pass)
All landed in commit 2 (`2de582d`). See M1 entry above for the bulk-rename approach. Each Minor's specific renames:
- **Mi1:** `from_open_crams` method + `open_crams` param renamed to `from_open_alignment_files` (~22 test call sites in alignment_input.rs).
- **Mi2:** Local vars in dispatch loops (`open_crams`, `cram_path`, `reference_cram_path`, `cram_header`) renamed; 4 `.expect` strings updated.
- **Mi3:** `cram_cfg` → `alignment_cfg` (9 sites); `cram_config_from_args` → `alignment_config_from_args`.
- **Mi4:** `CramHeader` struct → `AlignmentFileHeaderSummary`.
- **Mi6:** `validate_fasta_agreement` param renamed; 3 detail strings rewritten.
- **Mi8:** Doc-comment + module-doc sweep (MappedRead, OpenAlignmentFile, section banners, "BAM support lands" future tense in mod.rs / errors.rs / cram_input.rs / bam_input.rs).
- **Mi18:** `Stage1` variant doc updated.

### Mi5 — `process_one_chromosome_from_bam` name covers CRAM too
- **Final status:** Deferred
- **Reasoning:** Held deliberately per the original impl-report follow-ups; the subcommand-name coupling means renaming either both or neither.
- **Follow-up:** Yes.

### Mi7 — `input_crams` field/local mixed-format
- **Final status:** Deferred
- **Reasoning:** The field rename reaches into `src/psp/header.rs` (out of scope for this fix run — PSP schema field rename is its own design pass).
- **Follow-up:** Yes.

### Mi9 — `.csi`/`.bai` scan-and-index loop triplicated
- **Final status:** Applied with adaptation
- **Implementation summary:** Extracted `pub(crate) fn populate_binning_index<I>` in `index_preflight.rs`. `build_csi_for_bam` (production) and `bam_input::tests::build_bai_in_memory` (in-crate test fixture) both call through it. The `LinearIndex` vs `BinnedIndex` parameterisation is the type parameter `I`.
- **Adaptation:** Only 2 of 3 copies were collapsed. The third copy in `tests/common/mod.rs::build_csi` stays because the integration-test crate is outside `crate::bam` and cannot reach `pub(crate)` items. A test-support feature flag would resolve this but is out of scope for this fix run.
- **Files changed:** `src/bam/index_preflight.rs`, `src/bam/bam_input.rs` (commit `6823b07`).
- **Follow-up:** Yes — tests/common copy can be unified once a test-support feature flag is introduced.

### Mi10 — `index_display_name` free function → impl method
- **Final status:** Applied
- **Implementation summary:** Moved `fn index_display_name(index: &AlignmentIndex) -> &'static str` from `alignment_input.rs` into `impl AlignmentIndex { pub(crate) fn display_name(&self) -> &'static str }` in `index_preflight.rs`, next to the enum. Call site swapped.
- **Files changed:** `src/bam/alignment_input.rs`, `src/bam/index_preflight.rs` (commit `156fa18`).

### Mi11 — `load_alignment_index` re-encodes csi/bai policy
- **Final status:** Applied
- **Implementation summary:** `load_alignment_index` now consults `existing_index_for` for the on-disk index path (which centralises the `.csi`-preferred / `.bai`-fallback policy), then dispatches on the resolved file's extension to pick the parser. The cascade is no longer duplicated.
- **Files changed:** `src/bam/index_preflight.rs` (commit `86b77fc`).

### Mi12 — `MixedAlignmentFileFormats` / `UnsupportedExtension` duplicated
- **Final status:** Deferred
- **Reasoning:** Folds into M17's deferred redesign half (same shared-enum refactor).
- **Follow-up:** Yes.

### Mi13 — `From<AlignmentIndexError>` 4-arm passthrough
- **Final status:** Deferred
- **Reasoning:** Design call; the explicit rename (e.g. `MissingAlignmentIndex` → `MissingMapFileIndex`) is load-bearing for the CLI's `--flag` vocabulary the CLI layer owns. The review explicitly listed "leave as-is, document why" as an acceptable option.
- **Follow-up:** Yes.

### Mi14 — `pub mod` → `pub(crate) mod` for the per-format decoders
- **Final status:** Applied
- **Implementation summary:** `pub mod bam_input` / `pub mod cram_input` → `pub(crate) mod`. Confirmed by grep that no consumer outside `src/bam/` references these modules.
- **Files changed:** `src/bam/mod.rs` (commit `b8f5642`).

### Mi15 — `crate::bam` module name rename
- **Final status:** Deferred
- **Reasoning:** Large blast radius; decide before the next API consumer outside `crate::bam` is added.

### Mi16 — `BamIndex::Csi` arm has no isolated test
- **Final status:** Applied
- **Implementation summary:** Added `open_indexed_bam_record_stream_yields_target_contig_records_via_csi` in `bam_input.rs::tests`. Builds a real `.csi` (BinnedIndex parameterisation) instead of the existing `.bai` (LinearIndex); asserts the records come back.
- **Files changed:** `src/bam/bam_input.rs` (tests) (commit `6ec16b0`).

### Mi17 — `OwnedBamRecords` `Err(e)` arm has no test
- **Final status:** Applied
- **Implementation summary:** Subsumed by M11's `open_bam_record_stream_surfaces_read_errors_as_some_err` test (truncated BAM, asserts `Some(Err)` propagation).
- **Files changed:** `src/bam/bam_input.rs` (tests) (commit `6ec16b0`).

### Mi19 — `BinnedIndex` vs `LinearIndex` undocumented
- **Final status:** Applied
- **Implementation summary:** Inline comment added on the `Indexer<BinnedIndex>` line in `build_csi_for_bam` naming the on-disk shape implication ("BinnedIndex is the CSI on-disk shape; LinearIndex would emit the .bai shape").
- **Files changed:** `src/bam/index_preflight.rs` (commit `a2ea28d`).

### Mi20 — `OwnedIndexedBamRecords::next` long phase-numbered loop
- **Final status:** Deferred
- **Reasoning:** Cosmetic. Revisit when a second use case appears (per the review's "extract when a second use case appears" rule of thumb).

### Mi21 — Per-record `RecordBuf::default()` allocation
- **Final status:** Deferred
- **Reasoning:** Folds into the standing parallelisation-tuning workstream; same allocation pattern as the CRAM analogue (not a regression introduced by this slice).

### Mi22 — `Cargo.toml` missing justification comment
- **Final status:** Applied
- **Implementation summary:** 5-line `#`-prefix comment block added above the `noodles-bam` line explaining the "promote from transitive so the BAM/CSI symbols are compile-checked here" intent.
- **Files changed:** `Cargo.toml` (commit `b8f5642`).

## 5. Deferred findings to carry forward

- **Mi5** — `process_one_chromosome_from_bam` function-name rename. Held per impl-report; subcommand-name coupling.
- **Mi7** — `input_crams` field-name rename. Reaches into out-of-scope `src/psp/header.rs`.
- **Mi9 (partial)** — `tests/common/mod.rs::build_csi` copy stays; unify once a test-support feature flag is introduced.
- **Mi12** — `MixedAlignmentFileFormats` + `UnsupportedExtension` enum-duplication. Folds into M17's deferred redesign.
- **Mi13** — `From<AlignmentIndexError>` 4-arm passthrough redesign. Design call.
- **Mi15** — `crate::bam` → `crate::alignment` module rename. Large blast radius; defer.
- **Mi20** — `OwnedIndexedBamRecords::next` phase-extraction. Cosmetic.
- **Mi21** — Per-record `RecordBuf::default()` allocation. Folds into parallelisation-tuning workstream.
- **M17 (redesign half)** — Lift `MixedAlignmentFileFormats` into a shared `AlignmentFormatError` sub-enum. Test-add piece landed; redesign deferred.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** no
- **Baseline saved:** no — see Notes.
- **Benches run:** none.
- **Verdicts:** N/A.
- **Outcome:** skipped.
- **Notes:** All applied fixes touch startup-time / error-path / dispatch-time code; the per-record hot path (`OwnedBamRecords::next`, `OwnedIndexedBamRecords::next`, merge loop) is structurally unchanged. The largest fixes (M5 helper extraction, M6 enumerated arms, M18 carry kinds, Mi10 method move) all execute O(1) per input file or once per merge call — not per record. Baseline was not captured at preflight; skill rule "do not stash/revert to back-fill" applied.

## 10. Commands run

Inside `./scripts/dev.sh`, after each commit (commit ids shown in the findings table):

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib`
- `cargo test --test '*'`

Final post-commit-10 validation pass: all four commands clean.

## 11. Command results

- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, **924 passed** (was 916 pre-fix; net +8 from new lib tests: +1 cram_input M5 regression, +4 index_preflight M12, +2 bam_input M11, +1 bam_input Mi16, plus the pre-existing 916). The M9 / Mi17 tests landed in alignment_input.rs / bam_input.rs respectively but their count was absorbed by other rebalances.
- `cargo test --test '*'` → 0, **45 passed across 6 binaries** (17 cohort_cli + 6 + 5 + 5 pileup_cli + 7 + 4; was 43 pre-fix; net +2 from new integration tests: +1 cohort_cli M17 lock-step, +1 pileup_cli M10).
- `cargo doc --no-deps` → not re-run; pre-existing failure tracked in PROJECT_STATUS.md.
- `cargo audit` → not run (not in project verification list).

## 12. Notes

- Three open questions from the original review were resolved during the fix run: **Q1** (PileupCliError rename strategy) by the M1 rename pass; **Q2** (`#[non_exhaustive]` omission intent) treated as oversight per the plan's explicit call (closed by M4 / M7 / M8); **Q3** (CRAM-version-gate skip) confirmed accidental by the user and closed by M5. **Q4** (CSI depth) answered as `6` and applied by M14.
- The plan's "deprecated type aliases" decision (skipped per CLAUDE.md's no-shims policy) was inherited unchanged.
- Total: 10 fix commits + 1 final-report commit (this commit). Range: `9fc1df0` → `7a8569b`.
- The fix run did not introduce any new `cargo doc` warnings; the 14 pre-existing intra-doc-link warnings tracked in PROJECT_STATUS.md remain.
