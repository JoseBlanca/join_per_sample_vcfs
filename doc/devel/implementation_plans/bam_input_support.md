# BAM input support

Implementation plan for accepting BAM files as a Stage-1 mapped-read
input alongside the existing CRAM pipeline. Today every subcommand
that reads alignments (`pileup`, `var-calling-from-bam`) is
CRAM-only; this plan delivers `.bam` ingestion behind the same
interfaces, with the same per-read filter cascade and the same
per-chromosome parallel driver.

Companion specs: [per_sample_caller.md](../specs/per_sample_caller.md)
§"Inputs" and §"Read filters"; the existing CRAM slice plan
[per_sample_caller_cram_input.md](per_sample_caller_cram_input.md)
sets the policy this plan is preserving.

The module tree at [src/bam/](../../src/bam/) is already named for
BAM's eventual arrival — the doc comments in
[src/bam/mod.rs](../../src/bam/mod.rs),
[src/bam/cram_input.rs](../../src/bam/cram_input.rs), and
[src/bam/index_preflight.rs](../../src/bam/index_preflight.rs)
explicitly anticipate this work and call out the growth points.
This plan executes against those growth points.

## Scope

In:

- A new BAM record stream that decodes a coordinate-sorted `.bam`
  file into `sam::alignment::RecordBuf`s and plugs into the existing
  `OpenAlignmentFile { path, records: Box<dyn Iterator<Item=
  io::Result<sam::alignment::RecordBuf>> + Send> }` seam (renamed
  from `OpenCram`).
- Header validation for BAM mirroring the CRAM path: `@HD
  SO:coordinate`, `@SQ` cross-file consistency, FASTA contig
  agreement, `@RG SM` single-sample harvest. The validators already
  take a `sam::Header` and are format-agnostic — they get reused
  unchanged.
- `.csi` / `.bai` index detection and opt-in build via
  `--build-map-file-index`, slotted into `index_preflight`.
- Per-chromosome random-access query via `noodles_bam` index seek,
  plugged into the existing per-chromosome parallel driver
  (`var_calling_from_bam`).
- One refactor pass: rename `cram_input.rs` →
  `alignment_input.rs`, lift `AlignmentMergedReader` (formerly
  `CramMergedReader`) into the shared module, and split CRAM-only
  decoder bits into a `cram_input.rs` that only contains
  `OwnedCramRecords` / `OwnedIndexedCramRecords` + open helpers.
- CLI rename: the positional argument exposed today as `crams:
  Vec<PathBuf>` (help text `[CRAMS]...`) becomes `alignment_files:
  Vec<PathBuf>` (help text `[ALIGNMENT_FILES]...`) on both
  `pileup` and `var-calling-from-bam`. Hard-error when a single
  invocation mixes `.cram` and `.bam` inputs.
- Integration-test coverage: every existing CRAM-path integration
  test gets a `_bam` sibling that exercises the same scenario
  against a synthetic `.bam` fixture.

Out:

- Mixed CRAM+BAM inputs in one invocation. The pre-flight will
  fail with `MixedAlignmentFileFormats` and name both offending
  paths. (See "Design decisions" §1.)
- BAM-format CLI subcommand renames (`var-calling-from-bam` stays
  named as-is — the historical name already covered the
  yet-to-come BAM case).
- SAM (plain-text `.sam`) input. Stays unsupported; the
  pre-flight already rejects unknown extensions.
- BAM reading via stdin or a URL. File paths only, as for CRAM.
- BAM-specific perf-tuning. The CRAM path's per-chrom parallel
  ceiling (`var_calling_from_bam_per_chromosome` impl report —
  3.85× at T=13, gated by per-chrom workload imbalance) is the
  baseline; BAM should land at the same or better, but no
  dedicated bench commit is included. The standing
  parallelisation-tuning pass picks this up.
- The deferred `examples/profile_from_bam_e2e.rs` /
  `benches/from_bam_e2e_perf.rs` infrastructure (commit 5 of the
  earlier per-chromosome plan — scoped out then; still scoped
  out now). Once that lands it covers both formats.

## Design decisions

### 1. No mixing of CRAM and BAM inputs per invocation

Hard-error if a single command line names both formats. Rationale:
the merge produces one coordinate-sorted record stream that
downstream stages consume as one sample; a single sample carrying
some lanes as CRAM and others as BAM is a corner case the
project does not need to optimise for, and the all-or-nothing
guarantee simplifies both the pre-flight code and the user-facing
error reporting (one format keyword in the error message, no
ambiguity about which file's reference panel was used). If a real
need surfaces, lift the restriction in a follow-up — the merge
machinery is already format-agnostic at the `OpenAlignmentFile`
boundary.

The pre-flight raises a new typed error,
`AlignmentIndexError::MixedAlignmentFileFormats { first_path,
first_kind, other_path, other_kind }`. The CLI bridges this into a
user-facing `VarCallingFromBamCliError::MixedAlignmentFormats` /
`PileupCliError::MixedAlignmentFormats` variant whose `Display`
names both paths and the policy ("CRAM and BAM cannot be mixed in
one invocation").

### 2. `.csi` preferred, `.bai` fallback on read; `.csi` only on build

On the read path (`load_alignment_index` + the per-input handle
loader in `var_calling_from_bam`), pre-flight looks for
`<input>.bam.csi` first, then `<input>.bam.bai`, and loads
whichever exists. On the build path (`--build-map-file-index`),
`build_index` always emits `<input>.bam.csi`. Rationale: `.bai`
cannot index contigs longer than 512 Mb (the project does not
target those today, but several plant references — the project's
target species — are close to the limit; `.csi` removes the cap
permanently); meanwhile, existing pipelines often ship `.bai`
indexes already, so we accept those without forcing a rebuild.

The `AlignmentIndex` enum grows two variants:

```rust
#[non_exhaustive]
pub enum AlignmentIndex {
    Crai(Arc<noodles_cram::crai::Index>),
    BamCsi(Arc<noodles_csi::Index>),
    BamBai(Arc<noodles_bam::bai::Index>),
}
```

(Both BAM variants behave equivalently as opaque payloads handed
to a BAM `Reader::query`; the variant distinguishes only the
on-disk format.)

### 3. Renaming + shared-core extraction

[src/bam/cram_input.rs](../../src/bam/cram_input.rs) splits into:

- `src/bam/alignment_input.rs` (NEW) — `AlignmentMergedReader`
  (formerly `CramMergedReader`), `OpenAlignmentFile` (formerly
  `OpenCram`), `MappedRead`, `FilterCounts`,
  `AlignmentMergedReaderConfig` (formerly
  `CramMergedReaderConfig`), the merge loop, the filter cascade,
  the cross-file duplicate detector, all flag constants. Header
  validation helpers (`extract_header`, `validate_fasta_agreement`)
  live here too — they take a `sam::Header` and a `&Path` and care
  nothing about the source format.
- `src/bam/cram_input.rs` (REDUCED) — `OwnedCramRecords`,
  `OwnedIndexedCramRecords`, `open_cram_record_stream`,
  `open_indexed_cram_record_stream`, and the CRAM-specific
  `Repository`-bearing decoder seam. No `Reader` impls move out.
- `src/bam/bam_input.rs` (NEW) — `OwnedBamRecords`,
  `OwnedIndexedBamRecords`, `open_bam_record_stream`,
  `open_indexed_bam_record_stream`. Symmetric with
  `cram_input.rs`. No `Reader` impl, no merge logic.

The merge-side public API only mentions
`AlignmentMergedReader::new` / `::query` and takes `&[PathBuf]`; it
dispatches per-path on extension to either the CRAM or BAM open
helper, then plugs the resulting `Box<dyn Iterator<Item=
io::Result<sam::alignment::RecordBuf>>>` into the unchanged merge.
The per-input format is recorded once at open time and used to
choose which `extract_header` opener to call — for CRAM we still
need the `cram::Reader::read_file_header` path; for BAM we use
`noodles_bam::io::Reader::read_header`. Both return a
`sam::Header`, so `extract_header` (the validator) is unchanged.

### 4. CLI rename `crams` → `alignment_files`

Affected structs:

- `PileupArgs.crams` →
  [PileupArgs.alignment_files](../../src/pop_var_caller/cli.rs).
- `VarCallingFromBamArgs.crams` →
  [VarCallingFromBamArgs.alignment_files](../../src/pop_var_caller/var_calling_from_bam.rs).

Help text changes from `[CRAMS]...  Coordinate-sorted CRAM file(s)
for one sample` to `[ALIGNMENT_FILES]...  Coordinate-sorted CRAM
or BAM file(s) for one sample (all inputs in one invocation must
share the same format).`

The subcommand `var-calling-from-bam` keeps its name — the
historical naming was already anticipating BAM, and renaming the
subcommand would break callers' invocations without delivering
any new function. The error type
`VarCallingFromBamCliError` keeps its name for the same reason.

Field references inside the driver and Stage-1 pipeline
(`args.crams`, `&args.crams`, etc.) all rename. Integration tests
build their `Args` structs via field-name updates.

## Module layout

```
src/bam/
    mod.rs                  (REWRITE: re-exports change shape; doc
                             comment loses its "BAM is planned"
                             paragraph and gains a "format dispatch
                             happens in alignment_input.rs" pointer)
    alignment_input.rs      (NEW: AlignmentMergedReader, MappedRead,
                             FilterCounts, AlignmentMergedReaderConfig,
                             OpenAlignmentFile, format-agnostic
                             header validators, merge loop)
    cram_input.rs           (REDUCED: only OwnedCramRecords +
                             OwnedIndexedCramRecords + open_cram_*
                             helpers + CRAM-specific extract_header
                             header-shim wrapper)
    bam_input.rs            (NEW: OwnedBamRecords, OwnedIndexedBamRecords,
                             open_bam_record_stream,
                             open_indexed_bam_record_stream)
    index_preflight.rs      (EXTEND: AlignmentFileKind::Bam variant;
                             load_alignment_index handles .csi/.bai;
                             build_index writes .csi)
    errors.rs               (EXTEND: rename CramInputError →
                             AlignmentInputError with #[non_exhaustive]
                             and a deprecated type alias; new
                             AlignmentIndexError::MixedAlignmentFileFormats
                             variant)
```

[src/lib.rs](../../src/lib.rs) — no public-module change. The
`bam::cram_input` and `bam::bam_input` modules are private
implementation; only `bam::alignment_input::AlignmentMergedReader`
and `bam::index_preflight` need to be `pub use`'d from `bam/mod.rs`
for downstream callers.

## Dependencies (Cargo.toml)

Add two crates compatible with the existing noodles suite
(`noodles-cram = 0.93.0`, `noodles-sam = 0.85.0`,
`noodles-fasta = 0.61.0`):

- `noodles-bam` — BAM decoder, header reader, `Reader::query`.
- `noodles-csi` — `.csi` index reader/writer (likely already pulled
  transitively by `noodles-cram`; confirm and pin explicitly if so).

`noodles-bam`'s `bai` submodule provides `.bai` reading; we use it
through the same crate so no separate `noodles-bai` dep is needed.

Pin the exact versions at implementation time and record them in
this plan and in the impl report. Use the same noodles family as
`noodles-cram = 0.93`; consult `cargo search noodles-bam` to find
the matching release.

## Public API changes

### Renamed types (with deprecated aliases for one release)

In [src/bam/alignment_input.rs](../../src/bam/alignment_input.rs):

```rust
pub struct AlignmentMergedReader { /* fields, all format-agnostic */ }
pub struct AlignmentMergedReaderConfig { /* fields unchanged */ }

#[deprecated(note = "renamed to AlignmentMergedReader; BAM support
                     made the CRAM-specific name misleading")]
pub type CramMergedReader = AlignmentMergedReader;

#[deprecated(note = "renamed to AlignmentMergedReaderConfig")]
pub type CramMergedReaderConfig = AlignmentMergedReaderConfig;
```

`AlignmentMergedReader::new` and `::query` keep the same parameter
shape (the input paths slice is still `&[PathBuf]`; the helper
caller decides the format per-path inside).

The deprecated aliases get removed in the impl-report follow-up
once every in-tree caller has migrated (they will have, in this
PR — the aliases are belt-and-braces for any out-of-tree user;
project guideline §"Backwards-compatibility shims" says we may
delete them outright. **Decision: include the aliases in the
commit that does the rename, delete them in the same PR's
clean-up commit once tests are green.**)

`OpenCram` → `OpenAlignmentFile`. No deprecated alias — it is
`pub(crate)`.

### `bam::index_preflight` extensions

```rust
#[non_exhaustive]
pub enum AlignmentIndex {
    Crai(Arc<noodles_cram::crai::Index>),
    BamCsi(Arc<noodles_csi::Index>),
    BamBai(Arc<noodles_bam::bai::Index>),
}

// load_alignment_index dispatches on path extension; for .bam it
// tries .bam.csi first, then .bam.bai, returning the first that
// loads cleanly. If neither exists, the caller (preflight) is
// responsible for surfacing AlignmentIndexError::MissingAlignmentIndex.

// preflight_alignment_indexes: behaviour unchanged for CRAM. For
// BAM, the "missing" check accepts either .csi or .bai. The
// "build" path always writes .csi.
```

### Error-type renames

`CramInputError` →  `AlignmentInputError`. The variant names
(`OpenFailed`, `MissingFastaIndex`, `UnsupportedCramVersion`,
`NotCoordinateSorted`, `ContigListMismatch`, `FastaContigMismatch`,
`MissingSampleTag`, `MultipleSampleNames`,
`MultipleSampleNamesInFile`, `OutOfOrderRead`,
`DuplicateReadAcrossFiles`, `MalformedMd5`, `MalformedRecord`,
`Io`, `ContigNotInList`, `PerInputHandleCountMismatch`) stay as
they are — none of them name the format in a way that BAM
contradicts, except `UnsupportedCramVersion`, which stays
CRAM-specific (BAM is a single-version format and a BAM-side
analogue would be redundant).

Add one new variant:

```rust
#[error(
    "input file '{path}' is not a valid BAM file: {detail}"
)]
MalformedBam { path: PathBuf, detail: String },
```

(or fold this into `MalformedRecord` if the noodles-bam error
shape lets us — confirm at implementation time.)

A new deprecated type alias `CramInputError = AlignmentInputError`
is exported for one release per the same policy as the reader
rename.

`AlignmentIndexError::MixedAlignmentFileFormats { first_path,
first_kind, other_path, other_kind }` is added (variant name
self-explanatory). The two `_kind: AlignmentFileKind` fields carry
the variant; `AlignmentFileKind` is promoted from private to
`pub(crate)` for this purpose (still not re-exported).

## Implementation slices (commit plan)

Each numbered slice is one commit. All slices land on a feature
branch; PR opens after slice 6 is green.

### Commit 1 — rename + extract shared core (no new function)

Mechanical refactor; **no behaviour change**. Splits
`bam/cram_input.rs` into `bam/alignment_input.rs` +
`bam/cram_input.rs` (reduced). Renames `CramMergedReader` →
`AlignmentMergedReader`, `CramMergedReaderConfig` →
`AlignmentMergedReaderConfig`, `CramInputError` →
`AlignmentInputError`, `OpenCram` → `OpenAlignmentFile`. Adds
`#[deprecated]` type aliases for the public renames.

Updates every caller in-tree to use the new names (the deprecated
aliases exist so an out-of-tree consumer would still compile, but
in-tree we eat the rename).

Acceptance: `cargo fmt --check` clean; `cargo clippy --lib
--tests --all-features -D warnings` clean; `cargo test
--all-targets --all-features` — all 898 lib + every integration
test still pass with zero changes.

Net diff: ~+200 / -200 (most of the existing
`cram_input.rs` body moves verbatim to `alignment_input.rs`; the
move is rename-detected by git so the diff stays small).

### Commit 2 — `noodles-bam` + `noodles-csi` Cargo deps, no code change

Just `Cargo.toml` + `Cargo.lock`. Pin the exact versions that pair
with the rest of the noodles family. Sanity check: `cargo build`
inside the container compiles the new crates; `cargo test --lib`
unchanged.

This commit isolates the dependency churn from the code change so
a `git bisect` can attribute any "deps broke compilation" issue
without confounding by the new modules.

### Commit 3 — `bam_input.rs`: BAM record streams

Adds [src/bam/bam_input.rs](../../src/bam/bam_input.rs) with:

- `pub(crate) struct OwnedBamRecords` — owns
  `noodles_bam::io::Reader<bgzf::Reader<File>>` (exact wrapper type
  per noodles-bam API) + a `sam::Header`. Implements `Iterator<Item =
  io::Result<sam::alignment::RecordBuf>>`. Mirrors the EOF-latch
  pattern from `OwnedCramRecords` if noodles-bam's `read_record`
  has the same idempotency hazard; otherwise document the
  difference inline.
- `pub(crate) struct OwnedIndexedBamRecords` — owns a BAM `Reader`,
  an `Arc<AlignmentIndex::Bam*>`, and a
  `target_reference_sequence_id: usize`. Implements `Iterator`,
  driving `Reader::query(&header, &Region)` for the target contig
  and filtering decoded records to that contig.
- `pub(crate) fn open_bam_record_stream(path: &Path) ->
  Result<(sam::Header, Box<dyn Iterator<Item =
  io::Result<sam::alignment::RecordBuf>> + Send>),
  AlignmentInputError>` — opens the BAM, reads the header, builds
  an `OwnedBamRecords`, returns both.
- `pub(crate) fn open_indexed_bam_record_stream(path: &Path,
  header: Arc<sam::Header>, index: &AlignmentIndex,
  target_ref_seq_id: usize) -> Result<Box<dyn Iterator<Item =
  io::Result<sam::alignment::RecordBuf>> + Send>,
  AlignmentInputError>` — opens the BAM, seeks via the
  pre-loaded index, builds an `OwnedIndexedBamRecords`.

Adds unit tests in the module:

- A synthetic-BAM fixture (write 3 reads on 1 contig with the
  test helper from commit 5), open the stream, assert the
  emitted `RecordBuf` count + the qnames in order.
- Empty-BAM fixture: header only, zero records — assert the
  stream returns `None` immediately and does not re-enter
  noodles past EOF (the analogue of the CRAM EOF-latch test).
- Indexed variant: write a BAM with reads on contigs A and B,
  build a `.csi`, open the indexed stream for contig B only,
  assert only contig-B reads come back.

Acceptance: `cargo test --lib bam::bam_input` passes
(3+ new tests); no other test touches the new module.

### Commit 4 — dispatch in `alignment_input.rs` + index_preflight extension

Changes:

- `index_preflight.rs`: `AlignmentFileKind::Bam` variant added;
  `existing_index_for(&path, AlignmentFileKind::Bam)` looks for
  `.csi` then `.bai`; `target_index_path(&path,
  AlignmentFileKind::Bam) = <path>.csi`; `build_index(&path,
  AlignmentFileKind::Bam)` writes a `.csi`.
- `index_preflight.rs`: pre-flight gains a "no mixing" pass — if
  the first input's kind differs from any subsequent input's kind,
  return `AlignmentIndexError::MixedAlignmentFileFormats`. Done in
  the first classification pass, before the missing-index pass, so
  the error fires before any build is attempted.
- `index_preflight.rs`: `AlignmentIndex` enum gets `BamCsi` and
  `BamBai` variants; `load_alignment_index(&Path)` dispatches on
  extension and (for `.bam`) tries `.csi` then `.bai`.
- `alignment_input.rs`: the per-path open helper inside
  `AlignmentMergedReader::new` and `::query` dispatches on
  extension. For `.cram`, call the existing
  `open_cram_record_stream` (renamed in commit 1 if needed);
  for `.bam`, call `open_bam_record_stream`. The merge,
  header-validation, and filter cascade run unchanged.
- `alignment_input.rs`: `AlignmentMergedReader::query` accepts the
  `headers` slice as `&[Arc<sam::Header>]` and the `indexes` slice
  as `&[AlignmentIndex]` unchanged; per-input dispatch picks
  `open_indexed_bam_record_stream` vs
  `open_indexed_cram_record_stream` based on the same per-path
  extension.

Tests added:

- `index_preflight` unit tests for the BAM cases (existing CRAM
  tests stay green): preflight accepts `.bam.csi`; preflight
  accepts `.bam.bai`; preflight prefers `.csi` when both
  present; preflight builds `.csi` when missing-with-flag;
  preflight errors on missing-without-flag (names `.csi` first
  in the "looked for" message); preflight errors on mixed
  `.cram` + `.bam` inputs and names both offending paths.
- `alignment_input` unit tests: open a BAM via
  `AlignmentMergedReader::new`, iterate to exhaustion, assert
  the harvested sample name + contig list match what the test
  fixture wrote.

Acceptance: lib tests + every existing integration test green;
new tests pass.

### Commit 5 — test-fixture helpers for BAM in `tests/common/mod.rs`

Adds parallel helpers:

- `pub fn build_bam(dir: &Path, fasta_path: &Path, sample: &str,
  md5: Option<&str>, records: &[RecordBuf]) -> PathBuf` — writes
  a coordinate-sorted BAM using
  `noodles_bam::io::writer::Builder`. Reuses
  `build_sam_header(sample, md5)` (format-agnostic).
- `pub fn build_csi(bam_path: &Path) -> PathBuf` — builds a
  `.csi` next to the BAM via `noodles_bam::bai::fs::index` or the
  equivalent `noodles_csi` builder.
- `pub fn build_bai(bam_path: &Path) -> PathBuf` — same for
  `.bai`, used in one fallback-coverage test.

These helpers live next to the existing `build_cram` and reuse the
same `build_fasta` + `read_record` + `fixture_md5` underneath.

No production code touched in this commit. New helpers are
exercised in commit 6.

### Commit 6 — CLI rename + integration-test sweep

Renames:

- `PileupArgs.crams` → `PileupArgs.alignment_files` (field +
  every call site).
- `VarCallingFromBamArgs.crams` →
  `VarCallingFromBamArgs.alignment_files` (field + every call
  site).
- Help text in
  [src/pop_var_caller/cli.rs](../../src/pop_var_caller/cli.rs)
  updated to `Coordinate-sorted CRAM or BAM file(s) for one
  sample (all inputs in one invocation must share the same
  format).`
- Driver internals (`load_per_input_headers`,
  `load_per_input_indexes`) keep their names but rename their
  local args to match.

Integration tests:

- Every test in
  [tests/cohort_cli_integration.rs](../../tests/cohort_cli_integration.rs)
  that exercises `run_var_calling_from_bam` gets a `_bam` sibling
  built from `build_bam` + `build_csi` instead of `build_cram`.
  Cover at least: happy path, walker-error surfacing,
  missing-index-without-flag (now expects the message to name
  `.csi` first), missing-index-with-flag-builds (asserts a
  `.csi` appeared on disk).
- Every test in
  [tests/pileup_cli_integration.rs](../../tests/pileup_cli_integration.rs)
  that exercises `run_pileup` gets a `_bam` sibling.
- One new mixed-format integration test in
  `tests/cohort_cli_integration.rs`: pass one `.cram` and one
  `.bam` to `var-calling-from-bam`, assert
  `VarCallingFromBamCliError::MixedAlignmentFormats` (or the
  equivalent CLI bridge variant) with both paths in the
  message.
- One new `.bai`-fallback integration test: build a BAM with a
  `.bai` (no `.csi`), assert pre-flight accepts it, assert the
  pipeline runs to completion.

Acceptance: `cargo fmt --check` clean; `cargo clippy --lib
--tests --all-features -D warnings` clean; `cargo test
--all-targets --all-features` — all prior tests + the new BAM
siblings pass. Expect ~+12-15 new integration tests and
~+5-8 new lib unit tests across commits 3–6.

Wall-time validation on real multi-chrom BAMs is deferred — see
"Out of scope". The CRAM path's 3.85×-at-T=13 ceiling is the
working assumption; if BAM lands materially worse the standing
parallelisation-tuning pass will catch it.

## Tests (summary)

| Where | Count (delta) | What it covers |
|---|---|---|
| `src/bam/bam_input.rs` (mod tests) | +3-4 | Synthetic BAM open + iterate; empty BAM EOF-latch; indexed query yields only target-contig reads. |
| `src/bam/index_preflight.rs` (mod tests) | +6 | `.csi` accepted; `.bai` accepted; `.csi` preferred when both present; build emits `.csi`; missing-without-flag names `.csi` first; mixed CRAM+BAM rejected. |
| `src/bam/alignment_input.rs` (mod tests) | +1-2 | BAM-input round-trip via `AlignmentMergedReader::new` (sample name + contigs harvested correctly). |
| `tests/common/mod.rs` | n/a (helpers) | `build_bam`, `build_csi`, `build_bai`. |
| `tests/pileup_cli_integration.rs` | +N (mirror of existing CRAM tests) | Every CRAM `pileup` integration test gets a `_bam` sibling. |
| `tests/cohort_cli_integration.rs` | +N + 2 | Every CRAM `var-calling-from-bam` integration test gets a `_bam` sibling; plus 1 mixed-format-rejected test; plus 1 `.bai`-fallback test. |

`cargo doc --no-deps --lib` should stay at its existing
pre-existing-error count (no new warnings introduced by this
slice).

## Open questions

- **`noodles-bam` version pinning.** The exact compatible release
  pairs with `noodles-cram = 0.93.0` need to be confirmed at
  implementation time (cargo search + a build pass). Record the
  pinned versions in the impl report.
- **EOF-latch on `OwnedBamRecords`.** Check whether
  `noodles_bam::io::Reader::read_record` is idempotent at EOF (the
  CRAM analogue famously is not; see the
  [latch fix on OwnedCramRecords](../reports/implementations/var_calling_from_bam_per_chromosome_2026-05-24.md)).
  If BAM has the same hazard, mirror the latch with the same
  comment block; if not, document the asymmetry inline so a future
  reader does not "fix" a non-bug.
- **`AlignmentInputError::MalformedBam` shape.** Resolve at
  implementation time whether noodles-bam's malformed-record
  errors map cleanly into the existing `MalformedRecord` variant
  (which already carries a `path` + optional `qname` + source
  `io::Error`) or whether a separate `MalformedBam` variant gives
  better diagnostics. Default: fold into `MalformedRecord` unless
  a concrete diagnostic gain shows up.
- **Subcommand naming.** `var-calling-from-bam` keeps its name.
  Confirmed in the design pass — the historical name is already
  format-neutral in intent.

## Follow-ups (not in this PR)

- Wall-time validation on real multi-chrom BAMs (analogous to the
  cohort H1's 3.85×-at-T=13 figure on tomato CRAMs). Picked up by
  the deferred `examples/profile_from_bam_e2e.rs` /
  `benches/from_bam_e2e_perf.rs` infrastructure (commit 5 of the
  earlier per-chromosome plan, scoped out then).
- Delete the deprecated `CramMergedReader` /
  `CramMergedReaderConfig` / `CramInputError` type aliases once
  the change has shipped one release. (See "Public API changes"
  §Renamed types — current decision is to ship aliases in this
  PR and delete in the same PR's clean-up commit; revisit if any
  out-of-tree caller is in scope.)
- BAM-side perf tuning if measurements show BAM trails CRAM
  meaningfully. The merge / filter / writer pipeline is shared,
  so any divergence is in the noodles decoder cost — out-of-scope
  here, in-scope for the standing parallelisation-tuning pass.
- Lift the no-mixing restriction if a real workload surfaces that
  needs to merge CRAM and BAM inputs into a single sample. The
  pre-flight gate is the only place the restriction lives; the
  merge core is already format-agnostic.
