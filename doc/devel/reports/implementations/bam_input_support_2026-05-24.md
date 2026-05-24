# BAM input support — implementation report (2026-05-24)

Adds `.bam` ingestion alongside the existing `.cram` pipeline on
both `pop_var_caller pileup` and `pop_var_caller var-calling-from-bam`.
Stage 1's per-file decoder seam grows a BAM-side sibling; the merge,
filter cascade, header validation, and per-chromosome parallel
driver stay unchanged.

Plan: [bam_input_support.md](../../implementation_plans/bam_input_support.md).

## Plan

The plan from the linked document was executed as 7 commits (1
docs + 6 code) on `main`:

| Commit | Hash | Purpose |
|---|---|---|
| docs: plan | `b87ec89` | The plan itself + the `git mv` from cram_input.rs → alignment_input.rs |
| 1 | `18a9b9e` | Rename `Cram*` types → `Alignment*` family across the source tree |
| 2 | `266e79a` | Make `noodles-bam` + `noodles-csi` explicit Cargo deps |
| 3a | `4ad1e04` | Extract CRAM decoder bits into a sibling `cram_input.rs` module |
| 3b | `630ac7c` | Add `bam_input.rs` (`OwnedBamRecords`, `OwnedIndexedBamRecords`, open helpers + 4 unit tests) |
| 4 | `a4d1f6d` | Per-extension dispatch in `alignment_input.rs`; BAM branches + mixed-format error in `index_preflight.rs` |
| 5 | `bee6bc1` | `tests/common/mod.rs` BAM fixture helpers (`build_bam`, `build_csi`, `build_bai`) |
| 6 | `d0af049` | CLI positional rename `crams` → `alignment_files`; BAM integration-test siblings on both subcommands |
| 7 (this commit) | `344f1b2` | Internal helper-parameter renames `crams: &[PathBuf]` → `alignment_files: &[PathBuf]` (4 sites); `AlignmentInputError::PerInputHandleCountMismatch.crams` → `.inputs` field rename + Display update |
| 8 (docs) | `be3b38a` | `PROJECT_STATUS.md` Open-list update: close the two items above |

The plan's "commit 1" was split into two for tighter review: a
pure type rename (`18a9b9e`) and a structural split (`4ad1e04`).
The deprecated type aliases the plan had reserved as a transitional
shim were skipped per the project's `CLAUDE.md` policy on
backwards-compatibility shims (in-tree only; all in-tree callers
migrated in the same commit).

## Assumptions

These were silent choices the plan did not pin and that ended up
mattering during implementation:

- **`noodles-bam` version pinning.** Pinned to `0.89.0` and
  `noodles-csi` to `0.56.0` — the versions `noodles-cram 0.93.0`
  was already pulling transitively. No version bump to the
  noodles-* suite was needed.
- **No EOF-latch on `OwnedBamRecords`.** Verified via the
  noodles-bam source: `read_record_buf` calls `read_block_size`
  which uses `read_exact_or_eof` — cleanly idempotent at EOF
  (returns 0 on every post-EOF call, never an error). Documented
  this on the struct's doc comment with a pointer to the CRAM
  analogue's latch field, so a future reader doesn't "fix" a
  non-bug. Confirmed by a dedicated regression test
  `open_bam_record_stream_returns_none_at_eof_idempotently`.
- **`MalformedBam` shape.** Folded into the existing
  `AlignmentInputError::MalformedRecord` variant rather than a
  separate `MalformedBam` variant — the noodles-bam errors are
  `io::Error` shaped and the existing variant already carries
  `path` + optional `qname` + source.
- **`open_*_record_stream` symmetry.** Both `cram_input` and
  `bam_input` export the same factory shape: `(sam::Header,
  AlignmentRecordsIter)` for the linear opener and
  `AlignmentRecordsIter` for the indexed opener. The CRAM
  opener takes an extra `fasta::Repository` (CRAM decode needs
  it for sequence reconstruction); the BAM opener does not. The
  caller threads the repository through anyway because the
  downstream F1 mismatch-fraction filter on the merge side
  still needs it.
- **`AlignmentFileKind`** was promoted from private to
  `pub(crate)` so the merged-reader dispatch and the pre-flight
  classifier use the same enum. The `display_name` method went
  with it.
- **`index_display_name`** is a new private helper in
  `alignment_input.rs` that maps `&AlignmentIndex` → `"CRAI"` /
  `"CSI"` / `"BAI"` for the format-mismatch error message
  (`AlignmentIndexFormatMismatch`).
- **`load_per_input_headers`** in
  [`src/pop_var_caller/var_calling_from_bam.rs`](../../../src/pop_var_caller/var_calling_from_bam.rs)
  was found mid-implementation to be hard-coded to
  `noodles_cram::io::reader::Builder`. The BAM integration tests
  caught it on first run (every BAM-side test failed with
  "invalid CRAM header"). Fix is in the same commit (6) as the
  CLI rename: it now classifies each input via
  `AlignmentFileKind::from_path` and reads the header via the
  matching noodles-* reader.

## Changes made

### Modules / files

| File | Change |
|---|---|
| [src/bam/cram_input.rs](../../../src/bam/cram_input.rs) | NEW — `OwnedCramRecords`, `OwnedIndexedCramRecords`, `open_cram_record_stream`, `open_indexed_cram_record_stream`. Moved out of `alignment_input.rs` in commit 3a. |
| [src/bam/bam_input.rs](../../../src/bam/bam_input.rs) | NEW — `OwnedBamRecords`, `OwnedIndexedBamRecords`, `BamIndex` enum (`Bai` / `Csi`), `open_bam_record_stream`, `open_indexed_bam_record_stream`. The indexed iterator walks `BinningIndex::query` chunks manually so it owns the reader (`'static + Send` — same constraint as the CRAM analogue). |
| [src/bam/alignment_input.rs](../../../src/bam/alignment_input.rs) | EXISTING (renamed from `cram_input.rs` in the plan-commit) — gained per-extension dispatch in `AlignmentMergedReader::new` and `::query`, the new `AlignmentRecordsIter` shared type alias, `AlignmentPeekable` rename of the old `CramPeekable`, and the `index_display_name` helper. |
| [src/bam/index_preflight.rs](../../../src/bam/index_preflight.rs) | EXTENDED — `AlignmentFileKind::Bam` variant; `AlignmentIndex::BamCsi` + `::BamBai` variants; classification pass rejects mixed CRAM+BAM; `existing_index_for` BAM picks `.csi`-first / `.bai`-fallback; `target_index_path` BAM → `.csi`; `build_index` BAM → noodles-csi `Indexer<BinnedIndex>` pass + `noodles_csi::fs::write`. 7 new unit tests. |
| [src/bam/errors.rs](../../../src/bam/errors.rs) | EXTENDED — `AlignmentIndexError::MixedAlignmentFileFormats`, `AlignmentInputError::UnsupportedExtension` + `::MixedAlignmentFileFormats` + `::AlignmentIndexFormatMismatch`. |
| [src/bam/mod.rs](../../../src/bam/mod.rs) | Re-exports `cram_input` and `bam_input` as sibling modules under the `bam` root. |
| [src/pop_var_caller/cli.rs](../../../src/pop_var_caller/cli.rs) | `PileupArgs.crams` → `.alignment_files` + help text updated. Test helpers reflect the rename. |
| [src/pop_var_caller/var_calling_from_bam.rs](../../../src/pop_var_caller/var_calling_from_bam.rs) | `VarCallingFromBamArgs.crams` → `.alignment_files`; `VarCallingFromBamCliError::MixedAlignmentFormats` variant; `From<AlignmentIndexError>` bridge for `MixedAlignmentFileFormats`; `load_per_input_headers` dispatches per-extension. |
| [tests/common/mod.rs](../../../tests/common/mod.rs) | `build_bam(dir, sample, md5, records)`, `build_csi(bam_path)`, `build_bai(bam_path)` — parallel to the existing `build_cram` helper, reusing `build_sam_header` and `read_record`. |
| [tests/pileup_cli_integration.rs](../../../tests/pileup_cli_integration.rs) | + `happy_path_default_config_bam`, + `pileup_rejects_mixed_cram_and_bam`. |
| [tests/cohort_cli_integration.rs](../../../tests/cohort_cli_integration.rs) | + `var_calling_from_bam_happy_path_bam`, + `from_bam_errors_when_bam_index_missing_and_flag_unset`, + `from_bam_builds_missing_bam_index_when_flag_set`, + `var_calling_from_bam_accepts_bai_when_no_csi`, + `from_bam_rejects_mixed_cram_and_bam_inputs`. |
| [Cargo.toml](../../../Cargo.toml) | + `noodles-bam = "0.89.0"`, + `noodles-csi = "0.56.0"` (both already pulled transitively by `noodles-cram = 0.93.0`). |

### Behaviour changes

- `pop_var_caller pileup [ALIGNMENT_FILES]...` and
  `pop_var_caller var-calling-from-bam [ALIGNMENT_FILES]...` now
  accept `.bam` inputs as well as `.cram`. Same per-read filter
  cascade applies (BAM doesn't change the filter contract).
- A single invocation that mixes `.cram` and `.bam` inputs hard-
  errors before any reader opens, naming both offending paths and
  format strings. Surfaces as
  `AlignmentIndexError::MixedAlignmentFileFormats` from the
  pre-flight (var-calling-from-bam) or as
  `AlignmentInputError::MixedAlignmentFileFormats` from the
  reader (pileup, which does not run pre-flight).
- BAM index pre-flight is `.csi`-preferred / `.bai`-fallback on
  read, `.csi`-only on build (`--build-map-file-index`). `.bai`'s
  512 Mbp contig-length cap motivates the build-side choice.
- Per-chromosome parallelism on `var-calling-from-bam` works
  identically against `.bam` + `.csi` (or `.bam` + `.bai`) inputs
  as it does against `.cram` + `.crai`.

### Renames (public)

- `CramMergedReader` → `AlignmentMergedReader`.
- `CramMergedReaderConfig` → `AlignmentMergedReaderConfig`.
- `CramInputError` → `AlignmentInputError`.
- `OpenCram` → `OpenAlignmentFile` (`pub(crate)`).
- `bam/cram_input.rs` (4878 lines) → `bam/alignment_input.rs`
  (refactored) + new sibling `bam/cram_input.rs` (CRAM-only
  decoder bits) + new sibling `bam/bam_input.rs` (BAM-only
  decoder bits).
- `PileupArgs.crams` → `.alignment_files`.
- `VarCallingFromBamArgs.crams` → `.alignment_files`.

No deprecated type aliases retained. Internal helper params
named `crams: &[PathBuf]` in
[`stage1_pipeline.rs`](../../../src/pop_var_caller/stage1_pipeline.rs)
and
[`var_calling_from_bam.rs`](../../../src/pop_var_caller/var_calling_from_bam.rs)
were not renamed — they are local function-parameter names with
no API impact, and renaming the public `args.crams` field flows
through cleanly without touching them.

## Tests added / updated

| File | Count (delta) | What it validates |
|---|---|---|
| [src/bam/bam_input.rs](../../../src/bam/bam_input.rs) | +4 | `open_bam_record_stream` happy path (3 records in order); empty-BAM EOF-idempotency; indexed query yields only target-contig records; indexed query on a contig with no records returns nothing across repeated polls. |
| [src/bam/index_preflight.rs](../../../src/bam/index_preflight.rs) | +7 | `.csi` accepted; `.bai` accepted as fallback; `.csi` preferred when both present; build-canonical target index path is `.csi`-shaped; missing BAM index error names the `.csi` path; mixed CRAM+BAM rejected (both orderings). |
| [tests/pileup_cli_integration.rs](../../../tests/pileup_cli_integration.rs) | +2 | BAM happy-path mirror of the existing CRAM test (with `header.writer.input_crams == ["NA12878.bam"]`); mixed CRAM+BAM rejected via `PileupCliError::CramInput(AlignmentInputError::Mixed...)`. |
| [tests/cohort_cli_integration.rs](../../../tests/cohort_cli_integration.rs) | +5 | BAM happy-path mirror of `var_calling_from_bam_happy_path`; missing-BAM-index error names `.csi`; `--build-map-file-index` writes `.csi`-only; `.bai`-fallback works end-to-end through the per-chromosome driver; mixed CRAM+BAM rejected via `VarCallingFromBamCliError::MixedAlignmentFormats`. |

Net: +18 new tests over the pre-implementation baseline.

## Validation results

All commands run inside the project's dev container
(`./scripts/dev.sh`):

| Command | Result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib` | **915 passed** (904 pre-implementation baseline + 4 `bam_input` + 7 `index_preflight` bam-side tests) |
| `cargo test --test '*'` (all 6 integration binaries) | **43 passed** (36 pre-implementation baseline + 2 pileup + 5 cohort BAM-side tests) |

Each commit was validated end-to-end before the next one began;
no commit on the branch is broken on its own.

## Tradeoffs and follow-ups

### Decisions taken (vs the plan's open list)

- **Mixed CRAM+BAM rejection** lives in both the pre-flight
  (`AlignmentIndexError::MixedAlignmentFileFormats`) AND the
  merged reader's `new` (`AlignmentInputError::MixedAlignmentFileFormats`)
  because the `pileup` subcommand does not run pre-flight.
  Same display-name strings (`"CRAM"`, `"BAM"`) on both sides so
  the error messages match.
- **CSI vs BAI build choice**: only `.csi` is built by
  `--build-map-file-index`. `.bai` is read-side-only fallback. The
  build-side test `from_bam_builds_missing_bam_index_when_flag_set`
  asserts no `.bai` ever appears on disk.
- **No alignment-file extension sniff**: rely on `.cram` / `.bam`
  file extension; an unknown extension surfaces
  `UnsupportedExtension`. Same policy as the pre-existing
  CRAM-only world had.

### Closed mid-PR (commit `344f1b2`)

Two items the original plan had reserved as deferred follow-ups
landed inside this same PR scope after additional discussion
during the impl run:

- **Internal helper-parameter renames** from `crams: &[PathBuf]`
  to `alignment_files: &[PathBuf]` in
  [`stage1_pipeline.rs`](../../../src/pop_var_caller/stage1_pipeline.rs)
  and in three private helpers inside
  [`var_calling_from_bam.rs`](../../../src/pop_var_caller/var_calling_from_bam.rs).
  Internal-only; no API impact.
- **`PerInputHandleCountMismatch.crams: usize` field** in
  [`AlignmentInputError`](../../../src/bam/errors.rs) → renamed
  to `.inputs: usize` + Display string update. Public-error-type
  breaking change for any out-of-tree consumer (none in tree).

### Deferred to follow-ups

- **Wall-time validation on real multi-chrom BAMs.** No real-data
  BAM run analogous to the cohort H1's 3.85× at T=13 figure on
  tomato CRAMs. The
  `examples/profile_from_bam_e2e.rs` /
  `benches/from_bam_e2e_perf.rs` infrastructure (commit 5 of the
  earlier per-chromosome plan) was scoped out then and remains
  out of scope now. The standing parallelisation-tuning pass
  picks this up.
- **Lift the no-mixing restriction** if a real workload appears
  with a per-sample need to merge CRAM and BAM. The pre-flight
  gate is the only place the restriction lives; the merge core
  is already format-agnostic.
