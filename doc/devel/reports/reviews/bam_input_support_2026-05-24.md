# Code Review: bam_input_support
**Date:** 2026-05-24
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** BAM-file input support slice, commits b87ec89 → be3b38a (excluding the unrelated `4d0da97` project_history doc commit)
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** A 10-commit slice on `main` adding BAM-file ingestion alongside CRAM in `pop_var_caller`. PR-shaped multi-commit slice (no remote PR yet).
- **Reviewed against:** local `main` at commit `be3b38a`.
- **In-scope files:**
  - [src/bam/alignment_input.rs](../../../src/bam/alignment_input.rs) (focus on lines 244-272, 480-520, 565-905 — the AlignmentRecordsIter alias, dispatch loops, format-mismatch helper; the bulk of the file is pre-existing CRAM-only logic, renamed only)
  - [src/bam/bam_input.rs](../../../src/bam/bam_input.rs) (new — 458 lines)
  - [src/bam/cram_input.rs](../../../src/bam/cram_input.rs) (new sibling — extracted decoder + 2 open helpers)
  - [src/bam/index_preflight.rs](../../../src/bam/index_preflight.rs) (extended)
  - [src/bam/errors.rs](../../../src/bam/errors.rs) (new variants)
  - [src/bam/mod.rs](../../../src/bam/mod.rs) (re-exports)
  - [src/pop_var_caller/var_calling_from_bam.rs](../../../src/pop_var_caller/var_calling_from_bam.rs) (dispatch + CLI rename + error bridge)
  - [src/pop_var_caller/cli.rs](../../../src/pop_var_caller/cli.rs) (CLI rename only)
  - [src/pop_var_caller/stage1_pipeline.rs](../../../src/pop_var_caller/stage1_pipeline.rs) (internal rename only)
  - [tests/common/mod.rs](../../../tests/common/mod.rs) (BAM fixture helpers)
  - [tests/pileup_cli_integration.rs](../../../tests/pileup_cli_integration.rs) (BAM-side tests)
  - [tests/cohort_cli_integration.rs](../../../tests/cohort_cli_integration.rs) (BAM-side tests)
  - [Cargo.toml](../../../Cargo.toml) (deps)
- **Deliberately out of scope:**
  - Mechanical sed-renames in `benches/baq_perf.rs`, `examples/dhat_baq.rs`, `src/pileup/per_sample/*`, `src/pileup/walker/*`, `src/pop_var_caller/cli/error_bridge.rs`, `src/pop_var_caller/cli/shared_args.rs` — type-name change only, no logic.
  - Pre-existing CRAM-only logic in `alignment_input.rs` (merge loop, filter cascade, header validators, BAQ adapters, walker tests).
  - Commit `4d0da97` (project_history.md by user) — unrelated.
  - `doc/devel/*` files except the plan + impl report.
- **Categories dispatched (11, all applied):** reliability, errors, naming, defaults, idiomatic, refactor_safety, module_structure, unsafe_concurrency, smells, tooling, extras.

## 2. Verdict

**Approve-with-changes.** The slice's architecture is sound: clean per-format module split, `'static + Send` owned iterators, no `unsafe`, no shared mutable state. The Major findings cluster around four themes: (a) the rename pass left CRAM vocabulary in error displays + variant name `PileupCliError::CramInput` + several locals (mislabels every BAM error); (b) missing `#[non_exhaustive]` on three public/internal enums that should grow; (c) `load_per_input_headers` reimplements opener work without CRAM-version gating; (d) defaults around BAM index format (.csi/.bai policy + CSI depth) are body-buried rather than named constants. None of these is a correctness blocker for the BAM happy path, but the rename-debt items mislead users on every BAM error and the missing `#[non_exhaustive]` lock the public API into a breaking-change tax on future error variants.

## 3. Execution status

- `cargo fmt --check` — exit 0, clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0, clean.
- `cargo test --all-targets --all-features` — all pass: 915 lib tests; 17 + 6 + 5 + 4 + 7 + 4 = 43 integration; 3 doctests; 6 zero-test bench-compile targets.
- `cargo doc --no-deps --lib` — **fails** on 14 pre-existing `broken_intra_doc_links = "deny"` warnings; none introduced by this slice (`git blame` on `src/bam/alignment_input.rs:135` shows the only in-scope-file warning dates to commit `4e187f0a`, 2026-05-12). PROJECT_STATUS.md already notes this as a pre-existing issue.
- `cargo audit` — not run; not in project verification list.
- Findings labeled "Needs verification": 0.

## 4. Open questions and assumptions

1. **`PileupCliError::CramInput` rename strategy.** The variant wraps `AlignmentInputError` and is matched by name at `var_calling_from_bam.rs:495, 504, 510, 516, 526, 532, 680, 777` and `cli.rs:218`. Should the rename be applied this PR (mechanical sed; ~9 sites) or deferred to a coordinated naming pass? Affects **M1, M2, M3**.

2. **`#[non_exhaustive]` on `AlignmentInputError`.** The plan called for it (`bam_input_support.md:209-213`); the impl-report is silent on why it landed without. Was the omission intentional, or an oversight? Affects **M4**.

3. **`load_per_input_headers` CRAM-version-gating skip.** The driver-side helper reads `read_file_definition()` but does not check `version.major() == 3`, while `open_cram_record_stream` does. Was this intentional (deferred version check happens later) or accidental? Affects **M5**.

4. **`.csi` build-side depth.** `Indexer::default()` picks depth=5 → ~537 Mbp addressable, marginally above `.bai`'s 512 Mbp. The README rationale for choosing `.csi` over `.bai` cites "no contig-length limit", but at the default depth that's misleading. Is depth=5 the right pick for the project's target species, or should it be bumped (e.g. depth=6 → ~2^32 bp)? Affects **M14**.

## 5. Top 3 priorities

1. **M5** ([var_calling_from_bam.rs:486-545](../../../src/pop_var_caller/var_calling_from_bam.rs#L486-L545)) — `load_per_input_headers` reimplements per-format header opening AND skips the CRAM-version gating that `open_cram_record_stream` enforces. A CRAM 4.x file passes header-load and is only rejected later by `AlignmentMergedReader::query`. Fix: have the loader call into `open_cram_record_stream` / `open_bam_record_stream` (discarding the iterator), or factor out `read_*_header_only` helpers.

2. **M4** ([errors.rs:12](../../../src/bam/errors.rs#L12)) — `AlignmentInputError` is `pub` but not `#[non_exhaustive]`; sibling `AlignmentIndexError` is. This slice added three variants. Cheapest moment to add the attribute is now; retrofitting later is itself a breaking change for any consumer that exhaustively matched.

3. **M1** ([cli.rs:121-124](../../../src/pop_var_caller/cli.rs#L121-L124) + [errors.rs:14-99](../../../src/bam/errors.rs#L14-L99)) — Every BAM error currently prints as `stage 1: CRAM input: failed to open CRAM '/path/sample.bam': …`. Three layers of CRAM-mislabel on a file that is not CRAM. Fix is mechanical: rename `PileupCliError::CramInput` → `AlignmentInput`, drop `"CRAM"` from the six `AlignmentInputError` variant display strings that are now format-agnostic, rename `FastaContigMismatch.cram_path` field. ~30-line sed + ~10 call-site updates.

## 6. Findings

### Major

- **M1: [src/bam/errors.rs:14-99](../../../src/bam/errors.rs#L14-L99) + [src/pop_var_caller/cli.rs:121-124](../../../src/pop_var_caller/cli.rs#L121-L124) — Error vocabulary still says "CRAM" on BAM inputs**
- **Categories:** errors, naming, defaults, smells, extras (5-way convergent)
- **Confidence:** High
- **Problem:** `PileupCliError::CramInput(#[from] AlignmentInputError)` (cli.rs:123-124) with `#[error("CRAM input: {0}")]` now wraps BAM errors after the dispatch additions in `var_calling_from_bam.rs:495,504,510,516,526,532,680,777`. Six `AlignmentInputError` variants embed "CRAM" verbatim in their Display strings even though `bam_input.rs:189,195,222,231` raises them for BAM:
  - `NoInputs` ("at least one CRAM input is required", errors.rs:14)
  - `OpenFailed` ("failed to open CRAM '{path}'", errors.rs:17)
  - `NotCoordinateSorted` ("CRAM '{path}' is not coordinate-sorted", errors.rs:32)
  - `FastaContigMismatch` (errors.rs:42; field also named `cram_path`)
  - `MultipleSampleNames` ("multiple sample names across CRAMs", errors.rs:56)
  - `DuplicateReadAcrossFiles` ("duplicate read across CRAMs", errors.rs:91)
  Plus `MissingMd5` recommends "re-CRAM with a tool that emits M5" — meaningless for BAM (cli.rs:128-131 and var_calling_from_bam.rs:165-169).
- **Why it matters:** Every BAM error prints with three layers of mislabel. User sees `stage 1: CRAM input: failed to open CRAM '/path/sample.bam': ...`. Likely cause of confused bug reports. The fix is mechanical and consistent across the codebase.
- **Suggested fix:**
  ```rust
  // cli.rs
  #[error("alignment input: {0}")]
  AlignmentInput(#[from] AlignmentInputError),

  // errors.rs
  #[error("at least one alignment-file input is required")]
  NoInputs,
  #[error("failed to open alignment file '{path}': {source}")]
  OpenFailed { /* … */ },
  #[error("alignment file '{path}' is not coordinate-sorted (found SO:'{sort_order}')")]
  NotCoordinateSorted { /* … */ },
  #[error("FASTA '{fasta_path}' disagrees with alignment file '{alignment_file_path}': {detail}")]
  FastaContigMismatch { fasta_path: PathBuf, alignment_file_path: PathBuf, detail: String },
  #[error("multiple sample names across alignment files: …")]
  MultipleSampleNames { /* … */ },
  #[error("duplicate read across alignment files: …")]
  DuplicateReadAcrossFiles { /* … */ },
  ```
  Then sweep call sites (`var_calling_from_bam.rs:495,504,510,516,526,532,680,777`; `cli.rs:218`; `alignment_input.rs:432,449,460,666,984` for the `cram_path` field-rename). Keep CRAM-only wording on the genuinely-CRAM variants (`UnsupportedCramVersion`, `MalformedMd5`).

- **M2: [src/pop_var_caller/var_calling_from_bam.rs:128-129, 376, 554-559, 795-798](../../../src/pop_var_caller/var_calling_from_bam.rs#L128-L129) — `VarCallingFromBamCliError::Io(#[from] io::Error)` collapses multiple origins**
- **Categories:** errors
- **Confidence:** High
- **Problem:** `Io(#[from] io::Error)` is constructed from at least four unrelated sites: (a) `TempDir::new_in`, (b) the index-loader bridge at lines 554-559 (which already string-flattens the underlying error to fit), (c) the ref-fetcher constructor at 795-798 (same flattening), and (d) any `?` that produces an `io::Error`. The "Implement `From` only when one origin exists" rule is broken; two of the four sites already work around it by `io::Error::other(format!(…))`, which is exactly the leak-of-context the rule prevents.
- **Why it matters:** A user seeing `io: failed to load alignment index for '…': …` cannot programmatically distinguish "tempdir create" from "load index" from "ref fetcher build". Add a fifth `io::Error` site and the typed channel is gone for the rest too.
- **Suggested fix:** Replace the catch-all with operation-named variants:
  ```rust
  #[error("failed to create scratch directory under '{parent}': {source}")]
  ScratchDir { parent: PathBuf, #[source] source: io::Error },
  #[error("failed to load alignment index for '{path}': {source}")]
  LoadIndex { path: PathBuf, #[source] source: io::Error },
  #[error("failed to build reference fetcher for '{contig}': {source}")]
  RefFetcher { contig: String, #[source] source: io::Error },
  ```
  Map at each `?` site explicitly; drop the bare `Io` variant. This lets `load_per_input_indexes` carry the underlying error natively instead of flattening to a string.

- **M3: [src/bam/index_preflight.rs:72-112](../../../src/bam/index_preflight.rs#L72-L112) — `load_alignment_index` returns bare `io::Error` despite typed `AlignmentIndexError` variants existing for the same conditions**
- **Categories:** errors, idiomatic, smells (3-way convergent)
- **Confidence:** High
- **Problem:** The `pub fn` signature is `std::io::Result<AlignmentIndex>`. The "no BAM index" arm at lines 94-102 manufactures `std::io::Error::new(NotFound, format!(…))`; the "unsupported extension" arm at 104-110 does the same with `InvalidInput`. Sibling `preflight_alignment_indexes` in the same file uses typed `AlignmentIndexError::{MissingAlignmentIndex, UnsupportedExtension}` for exactly these conditions. The doc-comment at lines 65-71 explicitly acknowledges the wrong shape ("prefer running `preflight_alignment_indexes` first so the failure surfaces at the typed `AlignmentIndexError` layer rather than as a bare `io::Error` here"). The caller (`var_calling_from_bam.rs:554-558`) then wraps the bare `io::Error` in another `io::Error::other(format!(...))`, double-stringifying the cause.
- **Why it matters:** Part of the public API. Forces callers to format-and-stringify rather than match on the failure mode. The same information is already structured next door.
- **Suggested fix:** Return `Result<AlignmentIndex, AlignmentIndexError>`. Add an `AlignmentIndexError::LoadFailed { path, source: io::Error }` variant for the underlying `noodles_cram::crai::fs::read` / `noodles_csi::fs::read` / `noodles_bam::bai::fs::read` failures. Update the bridge in `var_calling_from_bam.rs:213-260` for the new variant.

- **M4: [src/bam/errors.rs:12](../../../src/bam/errors.rs#L12) — `AlignmentInputError` missing `#[non_exhaustive]`**
- **Categories:** errors, refactor_safety, extras (3-way convergent)
- **Confidence:** High
- **Problem:** `AlignmentInputError` is public and gained three new variants in this slice (`UnsupportedExtension`, `MixedAlignmentFileFormats`, `AlignmentIndexFormatMismatch`); sibling `AlignmentIndexError` at line 204 carries `#[non_exhaustive]`. The asymmetry is not principled — both surfaces are equally public. The plan (`bam_input_support.md:209-213`) explicitly called for `#[non_exhaustive]`; the impl-report's "Assumptions" section is silent on the omission.
- **Why it matters:** Every variant added to a non-`#[non_exhaustive]` `pub enum` is a semver-breaking change. Retrofitting `#[non_exhaustive]` later is itself a break for any consumer who exhaustively-matched.
- **Suggested fix:**
  ```rust
  #[derive(Error, Debug)]
  #[non_exhaustive]
  pub enum AlignmentInputError { /* … */ }
  ```

- **M5: [src/pop_var_caller/var_calling_from_bam.rs:486-542](../../../src/pop_var_caller/var_calling_from_bam.rs#L486-L542) — `load_per_input_headers` duplicates open-helper work + skips CRAM-version gating**
- **Categories:** module_structure, smells, idiomatic (3-way convergent)
- **Confidence:** High
- **Problem:** The driver classifies via `AlignmentFileKind::from_path` and then inlines the CRAM-specific `build_from_path` / `read_file_definition` / `read_file_header` sequence and the BAM-specific `build_from_path` / `read_header` sequence. The exact same per-format step sequence lives in `cram_input.rs::open_cram_record_stream` (lines 264-298) and `bam_input.rs::open_bam_record_stream` (lines 184-203). The header is even already returned from those helpers as the first tuple element. The map_err pattern `PileupCliError::CramInput(AlignmentInputError::OpenFailed { path: path.clone(), source })` repeats five times in 50 lines. Worst: `open_cram_record_stream` checks `version.major() == 3` (cram_input.rs:285-291) but the inline duplicate at `load_per_input_headers:509` skips the check — a CRAM 4.x file passes header-load and is only rejected later by `AlignmentMergedReader::query`.
- **Why it matters:** Two consequences. (1) Any future change to per-format open sequence has to land in two places, and they have already drifted (the version check). (2) Five identical map_err closures obscure the actual control flow.
- **Suggested fix:** Add `pub(crate) fn read_cram_header_only(path: &Path) -> Result<sam::Header, AlignmentInputError>` in `cram_input.rs` and `read_bam_header_only` in `bam_input.rs`, each calling the existing `open_*_record_stream` and discarding the iterator. The driver becomes a 10-line dispatch loop and the version-gating drift disappears.

- **M6: [src/bam/alignment_input.rs:897](../../../src/bam/alignment_input.rs#L897) — Catch-all `(file_kind, index) =>` arm in `query` dispatch silently absorbs future `AlignmentIndex` variants**
- **Categories:** refactor_safety, idiomatic, smells (3-way convergent)
- **Confidence:** High
- **Problem:** `AlignmentIndex` is `#[non_exhaustive]` with three variants today; `AlignmentFileKind` is a 2-variant internal enum. The match enumerates the three legal pairings then routes everything else through a wildcard `(file_kind, index) => Err(AlignmentIndexFormatMismatch …)`. Today the catch-all hits only the three genuine mismatches. Tomorrow when `AlignmentIndex::BamTabix` lands, the wildcard silently absorbs `(Bam, BamTabix)` as "format mismatch" — defeating the `#[non_exhaustive]` extensibility intent the slice's own plan documents.
- **Why it matters:** The whole point of `#[non_exhaustive]` on `AlignmentIndex` plus a parallel-growth `AlignmentFileKind` is that the compiler flags every dispatch site on extension. This one site cancels that guarantee.
- **Suggested fix:** Enumerate the genuine-mismatch arms explicitly:
  ```rust
  let owned_records = match (file_kind, alignment_index) {
      (AlignmentFileKind::Cram, AlignmentIndex::Crai(idx)) => { /* … */ }
      (AlignmentFileKind::Bam,  AlignmentIndex::BamCsi(idx)) => { /* … */ }
      (AlignmentFileKind::Bam,  AlignmentIndex::BamBai(idx)) => { /* … */ }
      (AlignmentFileKind::Cram, AlignmentIndex::BamCsi(_) | AlignmentIndex::BamBai(_))
      | (AlignmentFileKind::Bam, AlignmentIndex::Crai(_)) => {
          return Err(AlignmentInputError::AlignmentIndexFormatMismatch { /* … */ });
      }
      // REVIEW ON UPGRADE: AlignmentIndex is #[non_exhaustive];
      // adding a variant requires a new arm here.
  };
  ```
  A new variant then triggers an exhaustiveness check that points at this site.

- **M7: [src/bam/bam_input.rs:48](../../../src/bam/bam_input.rs#L48) — `BamIndex` enum missing `#[non_exhaustive]`**
- **Categories:** refactor_safety
- **Confidence:** High
- **Problem:** `BamIndex` is `pub(super)` and sits immediately next to the `#[non_exhaustive]` `AlignmentIndex` in `index_preflight.rs`. A maintainer reading `BamIndex` today vs `AlignmentIndex` immediately above will reasonably conclude the BAM-side set is closed. That is exactly the asymmetry that produces silent assumption-breakage in three months.
- **Suggested fix:**
  ```rust
  #[derive(Clone)]
  #[non_exhaustive]
  pub(super) enum BamIndex { /* … */ }
  ```
  Policy-only (no behaviour change for in-crate matches), but consistent with the surrounding type's openness signal.

- **M8: [src/pop_var_caller/cli.rs:121](../../../src/pop_var_caller/cli.rs#L121) — `PileupCliError` missing `#[non_exhaustive]`**
- **Categories:** refactor_safety
- **Confidence:** High
- **Problem:** `PileupCliError` is `pub` and the documented return type of `run_pileup`, and is wrapped by `VarCallingFromBamCliError::Stage1` (`#[from]` at var_calling_from_bam.rs:126). The sibling `VarCallingFromBamCliError` correctly carries `#[non_exhaustive]`; `PileupCliError` is the odd one out.
- **Suggested fix:** Add `#[non_exhaustive]`. Verify in-crate matches (cli.rs:218 and others) destructure specific variants — unaffected.

- **M9: [src/bam/alignment_input.rs:897](../../../src/bam/alignment_input.rs#L897), [errors.rs:188-201](../../../src/bam/errors.rs#L188-L201) — `AlignmentIndexFormatMismatch` typed error has no test**
- **Categories:** reliability, errors (2-way convergent)
- **Confidence:** High
- **Problem:** New variant in this slice; only safety net against driver-side index/file-format mismatches. `grep -rn AlignmentIndexFormatMismatch tests/ src/` returns only the definition and the construction site — zero coverage.
- **Why it matters:** The whole point of promoting a defensive panic into a typed error (per `feedback_no_logs_use_errors`) is exercised never. A future refactor could swap the display-name strings or drop one match arm without test failure.
- **Suggested fix:** See **Missing tests §8.1**.

- **M10: [src/bam/alignment_input.rs:600-618](../../../src/bam/alignment_input.rs#L600-L618), [:867](../../../src/bam/alignment_input.rs#L867) — `AlignmentInputError::UnsupportedExtension` test coverage gap on the pileup path**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** `pileup` does not run through pre-flight (the comment at alignment_input.rs:592-598 explicitly says so), so the `AlignmentMergedReader::new`-side `UnsupportedExtension` construction is the only gate. `grep -rn 'AlignmentInputError::UnsupportedExtension' tests/` returns zero. The parallel `AlignmentIndexError::UnsupportedExtension` is covered (`preflight_errors_on_unsupported_extension`) but exercises a different code path.
- **Suggested fix:** See **Missing tests §8.2**.

- **M11: [src/bam/bam_input.rs:120-172](../../../src/bam/bam_input.rs#L120-L172) — `OwnedIndexedBamRecords::next` chunk-walking has no test for multi-chunk traversal or IO-error propagation**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** Three branches uncovered:
  1. "Advance to next chunk" path (line 145, `if vp >= chunk_end { current_chunk_end = None; continue }`).
  2. "EOF mid-chunk" path (lines 152-158, `Ok(0)` with chunks still pending).
  3. "Read returned `Err`" path (line 169).
  The 4 existing tests use small enough record sets that one chunk likely covers everything. None forces a chunk boundary or surfaces a read error.
- **Suggested fix:** See **Missing tests §8.3 and §8.4**.

- **M12: [src/bam/index_preflight.rs:72-112](../../../src/bam/index_preflight.rs#L72-L112) — `pub fn load_alignment_index` has no direct unit tests**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** Part of the public API. Three primary paths (CRAM, BAM-csi-preferred, BAM-bai-fallback) and three error paths. Only coverage is end-to-end via the integration tests (each setting up one index format). A regression that reversed the csi/bai preference would fail no test, only the integration tests would catch it diffusely.
- **Suggested fix:** See **Missing tests §8.5**.

- **M13: [src/bam/index_preflight.rs:243](../../../src/bam/index_preflight.rs#L243) — `.csi`-preferred / `.bai`-fallback policy not in a named constant; encoded in two function bodies**
- **Categories:** defaults
- **Confidence:** High
- **Problem:** Preference order is encoded twice — in `existing_index_for` (lines 244-258) and again in `load_alignment_index` (lines 78-103) — and only mentioned in doc-comment prose ("`.csi` preferred over `.bai`. If both exist, `.csi` wins" at line 79). No `pub const BAM_INDEX_READ_PREFERENCE`, no `# Defaults` section on `load_alignment_index`'s docs. Two future maintainers can drift the two helpers apart.
- **Suggested fix:**
  ```rust
  /// Preference order when both `.csi` and `.bai` exist next to a
  /// BAM input. `.csi` wins because it has no 512 Mbp contig-length
  /// cap (some plant references exceed this); `.bai` is accepted on
  /// read for compatibility with existing data sets. Build always
  /// emits `.csi`.
  pub const BAM_INDEX_READ_PREFERENCE: &[&str] = &["csi", "bai"];
  pub const BAM_INDEX_BUILD_FORMAT: &str = "csi";
  ```
  Both helpers read the order from one place; `AlignmentIndex` enum docs reference the constant.

- **M14: [src/bam/index_preflight.rs:309-342](../../../src/bam/index_preflight.rs#L309-L342) — `build_csi_for_bam` uses `Indexer::default()`, hiding min-shift=14 / depth=5 — addressable contig length ~537 Mbp, barely above `.bai`'s 512 Mbp cap**
- **Categories:** defaults
- **Confidence:** High
- **Problem:** Line 318 reads `let mut indexer: Indexer<BinnedIndex> = Indexer::default();`. CSI v1 has two behaviorally significant parameters (`min_shift` defaulting to 14 → 16 kbp bins, `depth` defaulting to 5 → addressable ~2^29 ≈ 537 Mbp). The function's docstring at lines 304-308 describes what but not the addressable-length implications. The rationale at `target_index_path:263` ("`.csi` has no contig-length limit while `.bai`'s … tops out at 512 Mbp") is misleading in light of these defaults — at depth=5 `.csi` tops out around 537 Mbp.
- **Why it matters:** The first time someone indexes a contig past 537 Mbp this silently fails or produces a malformed index. The defaults rubric says "`Default` is discouraged — fields that drive runtime behavior must not hide behind a no-argument constructor."
- **Suggested fix:** Add named constants and explicit `Indexer` construction (if noodles' API exposes a 2-arg constructor; otherwise builder-style). See finding's full diff in the per-category file. Bump depth if the project's plant references warrant it.

- **M15: [src/pop_var_caller/var_calling_from_bam.rs:95](../../../src/pop_var_caller/var_calling_from_bam.rs#L95) — `--build-map-file-index` help text says "build it in place" without revealing `.csi`-only emission**
- **Categories:** defaults
- **Confidence:** High
- **Problem:** Help text reads "(`.crai` for CRAM; `.bai` or `.csi` for BAM)" suggesting either BAM index might be built. The build path actually only emits `.csi` (per `target_index_path:265-269` and `build_csi_for_bam`). User who expects `.bai` to be produced will find an unexpected `.csi` appear.
- **Suggested fix:** Rewrite the help block to spell out read-vs-build:
  ```rust
  /// If any input lacks its index, build it in place. On read we
  /// accept `.crai` for CRAM and `.csi` or `.bai` for BAM (.csi
  /// preferred); on build we always emit `.crai` for CRAM and
  /// `.csi` for BAM (.bai's 16 kbp bin grid tops out at 512 Mbp
  /// and is not future-proof). …
  ```

- **M16: [src/pop_var_caller/var_calling_from_bam.rs:185](../../../src/pop_var_caller/var_calling_from_bam.rs#L185) + [src/bam/errors.rs:207-214](../../../src/bam/errors.rs#L207-L214) — `MissingMapFileIndex` / `MissingAlignmentIndex` Display lists only one path but BAM loader consults two**
- **Categories:** defaults
- **Confidence:** High
- **Problem:** For BAM the loader (`existing_index_for`) consults `.csi` then `.bai`. The `expected_index_path` field is filled from `target_index_path` which returns only the `.csi`-canonical path. An operator who reads the error sees "looked for 'sample.bam.csi'" and may chase a `.csi` build pipeline when they already have a `.bai` that would have been accepted. The lower-level `load_alignment_index:94-102` DOES list both — inconsistent.
- **Suggested fix:** Either change `MissingAlignmentIndex` to carry `looked_at: Vec<PathBuf>`, or extend the Display string to mention the `.bai` fallback explicitly. The latter is the minimal-diff version.

- **M17: [src/bam/errors.rs:172-177, 240-245](../../../src/bam/errors.rs#L172-L177) — `MixedAlignmentFileFormats` exists in both `AlignmentInputError` and `AlignmentIndexError` with identical wording; only one path tested**
- **Categories:** errors, smells (2-way convergent)
- **Confidence:** Medium
- **Problem:** Pre-flight path raises `AlignmentIndexError::MixedAlignmentFileFormats` and is bridged; merged-reader path raises `AlignmentInputError::MixedAlignmentFileFormats` and flows through `PileupCliError::CramInput(...)`. The integration test `from_bam_rejects_mixed_cram_and_bam_inputs` exercises only the pre-flight path; `tests/pileup_cli_integration.rs:264` exercises the merged-reader path for the `pileup` caller. The display strings drift if either path's wording is updated and the other isn't (which the M1 fix will trigger when applied unevenly).
- **Suggested fix:** Lift the variant into a small shared `AlignmentFormatError` that both enums embed via `#[from]`, or have one re-export. Add a test that asserts both paths render identical strings (`format!("{e}")`) to lock against drift.

- **M18: [src/bam/alignment_input.rs:636-641](../../../src/bam/alignment_input.rs#L636-L641) — `.unwrap()` re-classifies a path that was just classified; defensive panic with no test**
- **Categories:** reliability, idiomatic, smells, errors (4-way convergent)
- **Confidence:** High
- **Problem:** Two-pass loop: classify-pass at lines 599-618 classifies every input, then drops the result; open-pass at line 636 re-runs `AlignmentFileKind::from_path(cram_path).unwrap()` with a comment-guarded `unwrap`. The invariant lives in a comment 30 lines away from the panic site. No `// PANIC-FREE:` marker per project convention.
- **Why it matters:** Future refactor that moves the classify pre-pass behind a conditional, or adds a new extension, breaks the invariant silently. Defensive panics in input-validation paths violate `feedback_no_silent_bugs`.
- **Suggested fix:** Harvest the kinds in the first pass and zip them in the second:
  ```rust
  let mut input_kinds: Vec<AlignmentFileKind> = Vec::with_capacity(alignment_files.len());
  // … in classify pass, push `kind` instead of dropping ...
  for (cram_path, kind) in alignment_files.iter().zip(input_kinds.iter().copied()) {
      let (header, records) = match kind {
          AlignmentFileKind::Cram => open_cram_record_stream(cram_path, repository.clone())?,
          AlignmentFileKind::Bam  => open_bam_record_stream(cram_path)?,
      };
      // …
  }
  ```
  Eliminates the `.unwrap()` and the duplicated `from_path` call.

- **M19: [doc/devel/reports/implementations/bam_input_support_2026-05-24.md:17-26, 199-210](../../../doc/devel/reports/implementations/bam_input_support_2026-05-24.md#L17-L26) — Impl-report commit table + "Deferred follow-ups" claims items that landed in 344f1b2 as still deferred**
- **Categories:** extras (diff-matches-stated-intent)
- **Confidence:** High
- **Problem:** The report's commit table at lines 17-26 lists 8 commits ending at `d0af049`. Two more in scope are missing: `344f1b2` (the rename + field-rename) and `be3b38a` (the status update). The "Deferred to follow-ups" list at lines 199-210 still describes the `PerInputHandleCountMismatch.crams → .inputs` field rename and the internal-helper-param renames as **deferred**, when both landed in `344f1b2` within the same review scope.
- **Why it matters:** A consumer reading the impl report ends up confused on both fronts. The "Deferred" claim could mislead someone reviewing the public API to think the field rename has not yet happened.
- **Suggested fix:** Add the two missing commit rows; move the two now-applied bullets out of "Deferred to follow-ups" into a new "Closed mid-PR" section (or into the existing "Decisions taken" sub-bullet).

### Minor

- **Mi1: [src/bam/alignment_input.rs:717, 697, 911](../../../src/bam/alignment_input.rs#L717) — `from_open_crams` method + parameter `open_crams` are stale**
  - **Categories:** naming, refactor_safety, idiomatic (3-way convergent)
  - **Confidence:** Medium (deliberately left in this slice per the impl-report follow-ups; recording here per the brief's request to weigh against that deferral).
  - Method takes `Vec<OpenAlignmentFile>` but is named `from_open_crams`. Five positional args at the call sites at 697, 911, plus ~20 test-code call sites. Fix: rename to `from_open_alignment_files`, parameter `open_alignment_files`. Mechanical.

- **Mi2: [src/bam/alignment_input.rs:622, 627, 642, 649, 654, 663, 669, 671, 678-681, 685, 687, 694, 860, 866, 868, 905-908](../../../src/bam/alignment_input.rs#L622) — Stale local var names (`open_crams`, `cram_path`, `reference_cram_path`, `cram_header`) in format-agnostic dispatch loops**
  - **Categories:** naming, idiomatic, refactor_safety, module_structure, smells (5-way convergent)
  - **Confidence:** High (also deliberately left this PR per the impl-report follow-ups). Same staleness as Mi1, at narrower scope. Rename to `open_alignment_files`, `input_path`, `reference_input_path`, `extracted_header`. Update the 4 `.expect(...)` strings to match. Body is fully self-contained.

- **Mi3: [src/pop_var_caller/var_calling_from_bam.rs:651, 425, 294, 678; cli.rs:160, 182, 236](../../../src/pop_var_caller/var_calling_from_bam.rs#L651) — Local var `cram_cfg` + factory `cram_config_from_args` carry an `AlignmentMergedReaderConfig` for both formats**
  - **Categories:** naming (cross-cat from idiomatic, module_structure)
  - **Confidence:** High. Rename `cram_cfg` → `alignment_cfg`, `cram_config_from_args` → `alignment_config_from_args`. Mechanical, ~9 sites.

- **Mi4: [src/bam/alignment_input.rs:272, 281, 307](../../../src/bam/alignment_input.rs#L272) — `CramHeader` struct still says "Cram" in a format-agnostic header summary**
  - **Categories:** naming
  - **Confidence:** High. Private to the module. `extract_header(path, sam_header) -> CramHeader` is called from the dispatch loop for both formats. Rename to `AlignmentFileHeaderSummary` or `PerInputHeader`.

- **Mi5: [src/pop_var_caller/var_calling_from_bam.rs:642](../../../src/pop_var_caller/var_calling_from_bam.rs#L642) — `process_one_chromosome_from_bam` function name covers CRAM too**
  - **Categories:** naming
  - **Confidence:** High (deliberately not renamed because the subcommand is `var-calling-from-bam`). Recording per the brief. The body has no BAM-specific code; the function's job is "per-chromosome alignment-file → VCF worker". Either rename to `process_one_chromosome_from_alignment_files` (keep subcommand) or rename the subcommand (more invasive). Pick at follow-up time.

- **Mi6: [src/bam/alignment_input.rs:411-415, 432-434, 448-451, 459-462](../../../src/bam/alignment_input.rs#L411-L415) — `validate_fasta_agreement` parameter `cram_path` + three "in CRAM" detail strings**
  - **Categories:** naming
  - **Confidence:** High. Helper called from `new()` at 691-695 where the dispatch loop has already accepted both formats. The path passed in may be a `.bam`. Rename parameter and update the three detail strings ("CRAM has" / "in CRAM"). Pairs with M1's `FastaContigMismatch.cram_path` field rename.

- **Mi7: [src/pop_var_caller/cli.rs:336-355](../../../src/pop_var_caller/cli.rs#L336-L355) — `input_crams` local + `WriterProvenance.input_crams` field carry mixed-format paths**
  - **Categories:** naming
  - **Confidence:** Medium. The field rename reaches into `src/psp/header.rs` (out of scope for this slice; the test at `tests/pileup_cli_integration.rs:226-228` already acknowledges the asymmetry in a comment). Rename the local + the parameter in `build_writer_header`; defer the public-PSP-schema field rename to a coordinated naming pass.

- **Mi8: [src/bam/alignment_input.rs:60-62, 68-72, 202, 235-244, 271, 408-410, 488, 499, 540-543, 762-783; mod.rs:6, 8, 14-18; bam_input.rs:14-15; cram_input.rs:14-15; index_preflight.rs:1-2; errors.rs:1-3](../../../src/bam/alignment_input.rs#L60-L72) — Doc-comments stale: "CRAM only" / "BAM support lands as" future-tense**
  - **Categories:** naming, idiomatic, module_structure, errors, extras (5-way convergent)
  - **Confidence:** High. Multiple doc-comments in alignment_input.rs and module-level docs in mod.rs / cram_input.rs / bam_input.rs / index_preflight.rs / errors.rs still say "CRAM-only" or "BAM support lands as a sibling submodule" when BAM has already landed. Mechanical doc-sweep pass: change "currently CRAM-only" → "CRAM + BAM" / re-tense "lands as" → "lives at"/"is provided by".

- **Mi9: [src/bam/index_preflight.rs:309-342, src/bam/bam_input.rs:325-358, tests/common/mod.rs:196-235](../../../src/bam/index_preflight.rs#L309-L342) — `.csi`/`.bai` scan-and-index loop triplicated**
  - **Categories:** smells, idiomatic, reliability (3-way convergent)
  - **Confidence:** High. Three near-identical bodies for "open BAM, walk records, build CSI/BAI Indexer". `tests/common/mod.rs:191-195` doc already admits the duplication. Extract a generic `pub(crate) fn build_binning_index<I: ReferenceSequenceIndex>(&Path) -> io::Result<Index<I>>` in `index_preflight.rs`; the test fixtures + library helper all call it.

- **Mi10: [src/bam/alignment_input.rs:479-485](../../../src/bam/alignment_input.rs#L479-L485) — `index_display_name` free function with a single call site; should be an inherent method on the enum**
  - **Categories:** smells
  - **Confidence:** High. Defined in alignment_input.rs, used exactly once at line 901. The symmetric helper `AlignmentFileKind::display_name` is right next to the enum it describes (index_preflight.rs:234-239). Move `index_display_name` body to `impl AlignmentIndex` in `index_preflight.rs`; swap the call site to `index.display_name()`; delete the free fn.

- **Mi11: [src/bam/index_preflight.rs:79](../../../src/bam/index_preflight.rs#L79) — `.csi`-then-`.bai` policy re-encoded in `load_alignment_index`; `existing_index_for` already centralises it**
  - **Categories:** smells
  - **Confidence:** High. The two policies match today by convention. Make `load_alignment_index` consult `existing_index_for` for the path, then dispatch on the resolved file's extension. Removes the second source of truth.

- **Mi12: [src/bam/errors.rs:155-162, 168-177, 225-231, 236-245](../../../src/bam/errors.rs#L155-L162) — `MixedAlignmentFileFormats` and `UnsupportedExtension` duplicated verbatim across `AlignmentInputError` and `AlignmentIndexError`**
  - **Categories:** smells, errors (2-way convergent)
  - **Confidence:** High. Same fields, same `&'static str` types, word-for-word `#[error]` messages. Each downstream consumer pays two `From`/`match` arms (see the `From<AlignmentIndexError>` bridge at var_calling_from_bam.rs:229-259). Extract `AlignmentFormatError { UnsupportedExtension, MixedFormats }` and embed via `#[from]` in both enums; one source of wording.

- **Mi13: [src/pop_var_caller/var_calling_from_bam.rs:229-259](../../../src/pop_var_caller/var_calling_from_bam.rs#L229-L259) — `From<AlignmentIndexError>` 4-arm passthrough**
  - **Categories:** smells
  - **Confidence:** Medium. Each arm is a 1:1 rename except `MissingAlignmentIndex` → `MissingMapFileIndex`, which adds remediation text in the receiving variant. Either replace four CLI variants with `#[error(transparent)] AlignmentIndex(#[from] AlignmentIndexError)` + Display wrapper for the remediation case, or add a one-line comment that the explicit rename is load-bearing for the `--flag` vocabulary the CLI owns.

- **Mi14: [src/bam/mod.rs:33-37](../../../src/bam/mod.rs#L33-L37) — `pub mod cram_input` / `pub mod bam_input` should be `pub(crate) mod`**
  - **Categories:** module_structure
  - **Confidence:** High. Confirmed via `grep -rn "crate::bam::{cram,bam}_input"` outside `src/bam/` — zero hits. The items inside expose their cross-format-visible symbols with `pub(super)` already, so `pub mod` over-claims the public API.

- **Mi15: [src/bam/mod.rs:1-31](../../../src/bam/mod.rs#L1-L31) — Module name `crate::bam` hosts format-agnostic alignment input; rename to `crate::alignment` is a follow-up**
  - **Categories:** module_structure
  - **Confidence:** Medium. The mod.rs doc owns the choice explicitly. Path strings still imply BAM at use-sites (`crate::bam::alignment_input::AlignmentMergedReaderConfig`). Defer the rename; before the next API consumer outside `crate::bam` is added, decide.

- **Mi16: [src/bam/bam_input.rs:402-457](../../../src/bam/bam_input.rs#L402-L457) — Indexed BAM unit tests build `.bai` only; no isolated test for the `BamIndex::Csi` arm**
  - **Categories:** reliability
  - **Confidence:** High. CSI is the production-preferred branch yet has no isolated unit-level coverage (only end-to-end). Duplicate the two existing indexed-BAM tests with a CSI-built index or parametrise.

- **Mi17: [src/bam/bam_input.rs:73-89](../../../src/bam/bam_input.rs#L73-L89) — `OwnedBamRecords` `Err(e)` arm of `next` has no test**
  - **Categories:** reliability
  - **Confidence:** Medium. The "no-EOF-latch-needed" design is load-bearing on noodles' idempotency contract; a truncated-BAM test would lock in the contract. Modest belt-and-braces.

- **Mi18: [src/pop_var_caller/var_calling_from_bam.rs:125-126](../../../src/pop_var_caller/var_calling_from_bam.rs#L125-L126) — `Stage1` variant doc still says "CRAM-input validation"**
  - **Categories:** errors
  - **Confidence:** High. Update to "alignment-input validation, reader open, header parse".

- **Mi19: [src/bam/index_preflight.rs:312, src/bam/bam_input.rs:331](../../../src/bam/index_preflight.rs#L312) — `BinnedIndex` vs `LinearIndex` choice undocumented**
  - **Categories:** defaults
  - **Confidence:** High. Production CSI build uses `Indexer<BinnedIndex>`; test helper `build_bai_in_memory` uses `Indexer<LinearIndex>`. Neither names the on-disk format. Add a one-line comment on each `Indexer<…>` line.

- **Mi20: [src/bam/bam_input.rs:120-173](../../../src/bam/bam_input.rs#L120-L173) — `OwnedIndexedBamRecords::next` is a 50-line numbered-phase loop; extraction would help future tuning**
  - **Categories:** idiomatic
  - **Confidence:** Medium. The state-machine design is correct for the `'static + Send` constraint. The three phases (advance + seek / bound check / read + filter) are candidates for small helpers (`advance_to_next_chunk`, `read_one_record_in_current_chunk`). Defer until a second use case appears.

- **Mi21: [src/bam/bam_input.rs:77, 150](../../../src/bam/bam_input.rs#L77) — Per-record `RecordBuf::default()` allocation in the indexed BAM hot path**
  - **Categories:** extras (hot-path)
  - **Confidence:** Medium. Construct a stable scratch field and clone into the yielded value. Same allocation pattern as the CRAM analogue, so not a regression — fold into the standing parallelisation-tuning workstream.

- **Mi22: [Cargo.toml:58, 62](../../../Cargo.toml#L58) — Missing justification comment on newly explicit `noodles-bam` / `noodles-csi` lines**
  - **Categories:** tooling
  - **Confidence:** Medium. Adjacent deps carry rationale comments (`bytemuck`, `mimalloc`, `wide`). The "promoted from transitive so BAM/CSI is compile-checked" intent should live next to the lines, not only in the impl plan.

### Nits

Grouped per skill convention. Apply via mechanical sweeps where applicable.

- **Doc-comment / string-literal sweep**: `MappedRead` doc still says "A read decoded out of a CRAM" (alignment_input.rs:69-72). Section banner "Per-CRAM lazy record stream" at 235-238. `mismatch_bq_floor` doc says "raw value from the CRAM" at 202. `ReadFingerprint` doc references "Duplicate-read detection across CRAMs" at 499. `query` doc 762-783 says "each input CRAM fresh". Bulk replace "CRAM" → "alignment file" wherever format-agnostic; keep CRAM where genuinely CRAM-specific.
- **Module-doc stale future tense**: `mod.rs:8` "BAM record-stream support lands as a sibling submodule", `mod.rs:14-18` "currently CRAM-only, inlined here", `errors.rs:1-3` "CRAM today; BAM support tracked in …", `cram_input.rs:17` "BAM support lands as a sibling `bam_input` module", `alignment_input.rs:257-265` "When BAM support lands, a `crate::bam::bam_input` sibling will host the BAM-side analogues" — re-tense to present.
- **`index_preflight.rs:1-2`**: lede still says "detect `.crai` next to each input CRAM" — mention `.csi`/`.bai`.
- **`tests/pileup_cli_integration.rs:78` vs `:203`**: 70-line near-clones (CRAM happy-path vs BAM happy-path). Acceptable as long as a third clone for a future format doesn't land.
- **`Cargo.toml`** exact-pin: consistent with project policy; no action.
- Other small idiomatic items (`first_seen: Option<(usize, AlignmentFileKind)>` could be a named struct; three `crai_path_for` / `csi_path_for` / `bai_path_for` wrappers could fold into an enum-dispatch helper) — acceptable as-is.

## 7. Out of scope observations

Pre-existing issues surfaced during this review but not in the BAM-input scope. File as follow-up.

- **`cargo doc --no-deps --lib` fails on 14 intra-doc-link warnings.** The deny gate `broken_intra_doc_links = "deny"` at `Cargo.toml:21` is firing across the codebase; the one in-scope-file warning (`src/bam/alignment_input.rs:135` `BaqSkipReason`) is **pre-existing** per `git blame` (commit `4e187f0a`, 2026-05-12). PROJECT_STATUS.md already tracks this. The full sweep is a candidate follow-up for a doc-link cleanup pass.
- **`src/bam/cram_input.rs:142-249` (`OwnedIndexedCramRecords::next`)** has the same per-branch test gap as the BAM analogue (no multi-chunk, no IO-error propagation, no empty-index-cursor unit test). Pre-existing logic moved verbatim. Symmetry follow-up.
- **`src/bam/index_preflight.rs:187` `eprintln!`** from a library function. Pre-existing CRAM-side behaviour (the build progress message). Logging policy belongs in the caller; pre-existing.
- **`src/bam/index_preflight.rs:166-168, 175-182` first-missing-only reporting.** With N missing indexes the user re-runs N times. Pre-existing; out-of-scope UX item.
- **`--build-map-file-index` symlink write risk.** The build writes a `.csi` (or `.crai`) next to the user-provided input via `noodles_csi::fs::write`. A symlink pointing outside the project tree would be honoured. Pre-existing CRAM-side pattern (`noodles_cram::crai::fs::write` at index_preflight.rs:294-302); not introduced by this slice.
- **`src/bam/alignment_input.rs:653-687` `.expect("…is set when … is Some")` patterns** (pre-existing; the invariants come from the first-iteration pattern). API-design refactor candidate (use type-level state instead of `Option`).
- **`WriterProvenance.input_crams` field name** in `src/psp/header.rs` — out-of-scope (PSP schema rename has reach beyond this slice).
- **CRAM-version-3 gate** at `cram_input.rs:285-291` is hard-coded; no `pub const SUPPORTED_CRAM_MAJOR: u8 = 3`. Pre-existing; marginal.

## 8. Missing tests to add now

Per the reliability sub-agent's challenge pass.

### 8.1 `query_returns_alignment_index_format_mismatch_on_cram_path_with_bam_index`
- **Target:** [src/bam/alignment_input.rs:897-903](../../../src/bam/alignment_input.rs#L897-L903) — `AlignmentIndexFormatMismatch` arm.
- **Input class:** `.cram` path paired with `AlignmentIndex::BamCsi` (or `.bam` with `Crai`).
- **Bug caught:** Regression in the match arms of `query` that swallows the mismatch; typo in display-name strings.
- **Spec:** Build a 1-contig CRAM via `tests/common::build_cram`, construct an empty `noodles_csi::Index` (the match fires before any seek), assert:
  ```rust
  let err = AlignmentMergedReader::query(
      &[cram_path], &fasta, contigs, "sample".into(),
      &[Arc::new(header)], &[AlignmentIndex::BamCsi(Arc::new(empty_csi))],
      "chr1", AlignmentMergedReaderConfig::default(),
  ).expect_err("mismatch must error");
  assert!(matches!(err,
      AlignmentInputError::AlignmentIndexFormatMismatch {
          file_format: "CRAM", index_format: "CSI", ..
      }));
  ```

### 8.2 `pileup_rejects_input_with_unknown_extension`
- **Target:** [src/bam/alignment_input.rs:600-618](../../../src/bam/alignment_input.rs#L600-L618) — `AlignmentInputError::UnsupportedExtension` on the pileup path.
- **Input class:** input vec containing one `.sam` (or extensionless) path.
- **Bug caught:** Future refactor that drops the classify pre-pass and lets the second-pass `.unwrap()` panic; classifier that starts treating unintended extensions as Cram/Bam.
- **Spec:**
  ```rust
  #[test]
  fn pileup_rejects_input_with_unknown_extension() {
      let dir = TempDir::new().unwrap();
      let p = dir.path().join("sample.sam");
      std::fs::write(&p, b"x").unwrap();
      let fasta = build_fasta(dir.path());
      let args = default_args(fasta, dir.path().join("out.psp"), vec![p.clone()]);
      let err = run_pileup(&args).expect_err("must reject sam");
      assert!(matches!(err,
          PileupCliError::CramInput(AlignmentInputError::UnsupportedExtension { path }) if path == p));
  }
  ```

### 8.3 `indexed_bam_iterator_walks_two_chunks_in_order`
- **Target:** [src/bam/bam_input.rs:130-172](../../../src/bam/bam_input.rs#L130-L172) — chunk-walking advance branch.
- **Input class:** BAM whose CSI assigns at least two chunks to one contig (~32 KiB of records to cross a BGZF block boundary, or hand-craft chunks via `Indexer::add_record`).
- **Bug caught:** Regression in the "fall through to next chunk" branch — e.g., setting `current_chunk_end = None` but forgetting to `continue`.
- **Spec:** Construct a synthetic BAM with enough records that BGZF emits at least two blocks; build the CSI; drive `open_indexed_bam_record_stream`; assert yielded records' positions form the full ordered set on the target contig.

### 8.4 `indexed_bam_iterator_surfaces_read_errors_as_some_err`
- **Target:** [src/bam/bam_input.rs:169](../../../src/bam/bam_input.rs#L169) — `Err(e)` arm.
- **Input class:** truncated BAM whose first indexed chunk decodes partially then fails.
- **Bug caught:** Future refactor that converted the `Err(e) => return Some(Err(e))` arm into silent loop continuation.
- **Spec:** Write a complete BAM, build its CSI, truncate the body bytes (preserve header + first BGZF block, cut the rest); open via `open_indexed_bam_record_stream`; assert at least one `Some(Err(_))` is yielded.

### 8.5 `load_alignment_index` triplet
- **Targets:** [src/bam/index_preflight.rs:72-112](../../../src/bam/index_preflight.rs#L72-L112). Three tests:
  - `load_alignment_index_returns_bam_csi_when_csi_present`
  - `load_alignment_index_falls_back_to_bai_when_no_csi`
  - `load_alignment_index_errors_with_not_found_when_no_bam_index`
- **Input classes:** (a) BAM + both .csi/.bai; (b) BAM + only .bai; (c) BAM + neither.
- **Bugs caught:** Reversed csi/bai preference order; silent fallback to a default variant; error path that names the wrong path.
- **Spec:** Use the existing `touch` pattern (`preflight_accepts_existing_csi`). For (c), assert `io::ErrorKind::NotFound` and that the message contains both `csi_path` and `bai_path` strings.

### 8.6 Bonus: `linear_bam_iterator_surfaces_read_errors_as_some_err`
- **Target:** [src/bam/bam_input.rs:86](../../../src/bam/bam_input.rs#L86) — `Err(e)` arm of `OwnedBamRecords::next`.
- **Input class:** BAM truncated mid-record after the header.
- **Bug caught:** Regression that swallowed `read_record_buf` error and returned `None` (silent termination), letting the merge believe the file ended cleanly.
- **Spec:** Mirror of §8.4 against `open_bam_record_stream`.

## 9. What's good

- **`unsafe_concurrency` came back clean.** No `unsafe`, no `static mut`, no interior mutability across thread boundaries. `Send`/`Sync` correctness is established structurally via two compile-time bounds (the `Box<dyn Iterator + Send>` trait object at [alignment_input.rs:253-254](../../../src/bam/alignment_input.rs#L253-L254) + the `rayon::par_iter` `Sync` capture at [var_calling_from_bam.rs:413-435](../../../src/pop_var_caller/var_calling_from_bam.rs#L413-L435)). The chunk-walking state is held by-value with no aliasing.
- **`OwnedBamRecords` EOF-no-latch design** ([bam_input.rs:80-95](../../../src/bam/bam_input.rs#L80-L95)) is well-documented with explicit rationale (`read_record_buf` calls `read_exact_or_eof`, idempotent at EOF) and a regression test (`open_bam_record_stream_returns_none_at_eof_idempotently`) that locks the contract.
- **Mixed-format rejection is enforced at both layers** ([alignment_input.rs:599-618](../../../src/bam/alignment_input.rs#L599-L618) for pileup, [index_preflight.rs:146-200](../../../src/bam/index_preflight.rs#L146-L200) for the indexed-query path). Pre-flight's classify-pass runs entirely before any disk write, so a `--build-map-file-index` run with one mixed-format input does not partially build the others. Integration test `from_bam_rejects_mixed_cram_and_bam_inputs` covers the pre-flight path end-to-end.
- **`AlignmentIndexError` is correctly `#[non_exhaustive]`** ([errors.rs:204](../../../src/bam/errors.rs#L204)) and the `From<AlignmentIndexError>` bridge in `var_calling_from_bam.rs:229-258` exhaustively handles all 4 variants — adding a 5th will force a compile error in the bridge.
- **The mid-implementation catch of `load_per_input_headers` being CRAM-only** (called out in the impl report's Assumptions) demonstrates the integration-test-first discipline catches latent dispatch bugs that unit tests wouldn't.

## 10. Commands to re-verify

Inside `./scripts/dev.sh`:

- `cargo fmt --check` (re-run after the M1 doc/error-string sweep)
- `cargo clippy --all-targets --all-features -- -D warnings` (re-run after each Major fix)
- `cargo test --lib` (re-run after M9-M12 + Mi16-Mi17 land — expect +6 to +8 new tests)
- `cargo test --test '*'` (re-run after §8.2 test lands in pileup_cli_integration)
- New invocations the review introduces:
  - `cargo test --lib bam::alignment_input::tests::query_returns_alignment_index_format_mismatch_on_cram_path_with_bam_index` (M9 / §8.1)
  - `cargo test --test pileup_cli_integration pileup_rejects_input_with_unknown_extension` (M10 / §8.2)
  - `cargo test --lib bam::bam_input::tests::indexed_bam_iterator_walks_two_chunks_in_order` (M11 / §8.3)
  - `cargo test --lib bam::index_preflight::tests::load_alignment_index_prefers_csi_over_bai` (M12 / §8.5)

### Author response convention

Address each finding by its identifier (e.g., "M5") with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the four open questions in §4 first.
