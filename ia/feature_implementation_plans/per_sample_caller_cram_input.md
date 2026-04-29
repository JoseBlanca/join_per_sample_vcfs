# Per-sample caller — CRAM input slice

Implementation plan for the first slice of Stage 1 of the multi-sample
calling pipeline: turning N coordinate-sorted CRAMs + a reference FASTA
into a single coordinate-sorted stream of decoded reads. Specified in
[per_sample_caller.md](../specs/per_sample_caller.md) §"Inputs",
§"Multi-CRAM ingestion", and §"Errors".

This slice covers header validation, decoding, the peek-and-scan
merge, and the per-read filter cascade from spec §"Read filters" —
all the cheap rejection-based filters that decide which reads enter
the pipeline. It does *not* include BAQ, the pileup walker, phase
chains, or `.psf` writing — those are separate slices.

## Scope

What this slice produces:

- A new module tree under [src/per_sample_caller/](../../src/per_sample_caller/).
- A `CramMergedReader` type that owns N CRAM decoders, validates their
  headers against each other and against a reference FASTA, applies
  the per-read filter cascade from spec §"Read filters", and yields
  surviving reads in coordinate order.
- A `CramMergedReaderConfig` struct with `Default` providing the spec's
  defaults (`min_mapq = Some(20)`, all flag-based drops on,
  `min_read_length = Some(30)`). Both thresholds are `Option`s so
  `None` cleanly expresses "no minimum".
- Enough of `MappedRead` to carry every field downstream stages need
  (CIGAR, SEQ, BQ, MAPQ, flags, mate position, source-file index).
- Hard errors for every failure mode listed in
  [per_sample_caller.md §Errors](../specs/per_sample_caller.md).

Per-read filters absorbed into this slice (Option X — see
[spec §Read filters](../specs/per_sample_caller.md)):

- Unmapped (`flag & 0x4`)
- Secondary (`flag & 0x100`)
- Supplementary (`flag & 0x800`)
- QC fail (`flag & 0x200`)
- Duplicate (`flag & 0x400`)
- MAPQ < `config.min_mapq`
- Read length 0 / no SEQ (via `config.min_read_length`)

What it does **not** include:

- The "first or last CIGAR op is `I`/`D`" rule. Spec says this only
  suppresses the *indel* contribution, not the whole read; it lives
  in the pileup walker, not as a filter.
- BAQ.
- Pileup walking, allele extraction, scalar accumulation.
- Phase-chain slot allocation.
- `.psf` writing.
- Parallel decoder threads (deferred per spec §Parallelism — start
  single-threaded, thread later behind the same iterator API).
- `--region` / `.crai`-driven seeking (no-`--region` streaming path
  only, per spec §CLI).

There is **no separate `read_filter.rs` module**. The filter cascade
is a property of which reads enter the pipeline, and its natural home
is the pipeline entry. The spec's `read_filter.rs` line in the module
layout is superseded by this plan.

## Dependencies (Cargo.toml)

Add:

- `noodles-cram` — CRAM decoder.
- `noodles-sam` — header types, alignment record types.
- `noodles-fasta` — FASTA + `.fai` reader.
- `noodles-core` — `Position`, region types (transitive but worth
  pinning).
- `anyhow` — used at the orchestrator boundary in
  `per_sample_caller/mod.rs`.

`thiserror` is already in the manifest.

Pin to the latest stable noodles release at implementation time and
record the version in this plan after the fact.

## Module layout

```
src/per_sample_caller/
    mod.rs             — re-exports + future CLI orchestration
    cram_input.rs      — this slice: header validation, peek-and-scan
                         merge, per-read filter cascade
    errors.rs          — CramInputError (thiserror)
```

[src/lib.rs](../../src/lib.rs) gets one new line: `pub mod per_sample_caller;`.

The `read_filter.rs` line in spec §"Module layout" is superseded by
this plan — see Scope above.

## Public API

### Default constants

All tunable defaults are named `pub const` items at module scope —
the Rust analogue of Python module-level `DEF_*` constants, but
compile-time and zero-cost. No magic numbers in `Default` impls or
in code; every threshold has a self-documenting name a reader can
follow.

```rust
/// Reads with MAPQ strictly below this are dropped. Matches bcftools'
/// default and the spec recommendation in §"Read filters".
pub const DEFAULT_MIN_MAPQ: u8 = 20;

/// Decoded SEQ length below this is dropped. Reads shorter than this
/// rarely contribute reliable alignments at the project's coverage
/// targets.
pub const DEFAULT_MIN_READ_LENGTH: u32 = 30;
```

If a future filter needs a tunable, add it here and reference it from
`Default`.

### `CramMergedReaderConfig`

`min_mapq` and `min_read_length` are `Option`s: `None` is the explicit
"no minimum" state. Using `0` as a disable sentinel would make
`Some(0)` and `None` redundant representations of the same behaviour;
`Option` makes the absent state structurally distinct from any
specific threshold.

```rust
pub struct CramMergedReaderConfig {
    /// `None` = no minimum; `Some(n)` = drop reads with MAPQ < n.
    pub min_mapq: Option<u8>,
    /// `None` = no minimum; `Some(n)` = drop reads with decoded SEQ
    /// length < n.
    pub min_read_length: Option<u32>,
    pub drop_unmapped: bool,        // flag & 0x4
    pub drop_secondary: bool,       // flag & 0x100
    pub drop_supplementary: bool,   // flag & 0x800
    pub drop_qc_fail: bool,         // flag & 0x200
    pub drop_duplicate: bool,       // flag & 0x400
}

impl Default for CramMergedReaderConfig {
    fn default() -> Self {
        Self {
            min_mapq: Some(DEFAULT_MIN_MAPQ),
            min_read_length: Some(DEFAULT_MIN_READ_LENGTH),
            drop_unmapped: true,
            drop_secondary: true,
            drop_supplementary: true,
            drop_qc_fail: true,
            drop_duplicate: true,
        }
    }
}
```

Call-site ergonomics via struct update syntax:

```rust
// All defaults
CramMergedReader::new(&crams, fa, CramMergedReaderConfig::default())?;

// Override one field — raise the threshold
CramMergedReader::new(&crams, fa, CramMergedReaderConfig {
    min_mapq: Some(30),
    ..Default::default()
})?;

// Disable the MAPQ filter entirely
CramMergedReader::new(&crams, fa, CramMergedReaderConfig {
    min_mapq: None,
    ..Default::default()
})?;
```

### `CramMergedReader`

Two constructors. The public `new` is the convenience that opens
files; the crate-internal `from_open_crams` is the testable core
that takes pre-built per-CRAM record streams. The split pushes I/O
to the program edge: production callers use `new`, tests inject
canned record streams through `from_open_crams` and skip the
synthetic-CRAM dance for everything except the I/O-glue tests.

```rust
pub struct CramMergedReader {
    // private: peekers, contig list, sample name, dup-detection window,
    //          config, per-file order tracking, filter counters.
}

impl CramMergedReader {
    /// Public, convenient. Opens every CRAM, validates headers against
    /// each other and the FASTA, sets up the noodles fasta repository,
    /// and delegates the merge/filter logic to `from_open_crams`.
    pub fn new(
        crams: &[PathBuf],
        fasta: &Path,
        config: CramMergedReaderConfig,
    ) -> Result<Self, CramInputError>;

    /// Crate-internal core. Takes pre-built per-CRAM record streams
    /// plus the canonical contig list and sample name (the caller is
    /// responsible for cross-file header validation — `new` does this
    /// before calling here). All I/O for opening lives outside this
    /// function, which is why it is the single point of entry tests
    /// drive against.
    pub(crate) fn from_open_crams(
        open_crams: Vec<OpenCram>,
        contigs: ContigList,
        sample_name: String,
        config: CramMergedReaderConfig,
    ) -> Result<Self, CramInputError>;

    pub fn sample_name(&self) -> &str;
    pub fn contigs(&self) -> &ContigList;

    /// Counts of reads dropped by each filter. Logged to stderr by the
    /// orchestrator at end-of-run per spec §"Errors" / Warnings.
    pub fn filter_counts(&self) -> &FilterCounts;
}

/// One CRAM's worth of already-decoded records, plus the path string
/// used for error messages. The `records` iterator is consumed lazily
/// by the merge — it is *not* drained eagerly. `Box<dyn Iterator>`
/// keeps the signature concrete so `Vec<OpenCram>` is a plain owned
/// type with no generic parameters or lifetimes.
pub(crate) struct OpenCram {
    /// Used in error messages naming the source of an offending read.
    pub path_for_errors: PathBuf,
    /// Lazy stream of decoded records. In production this wraps a
    /// `noodles_cram::io::Reader::records()`; in tests it wraps a
    /// `Vec<RecordBuf>::into_iter().map(Ok)`.
    pub records: Box<dyn Iterator<Item = io::Result<noodles_sam::alignment::RecordBuf>> + Send>,
}

pub struct FilterCounts {
    pub unmapped: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub qc_fail: u64,
    pub duplicate: u64,
    pub low_mapq: u64,
    pub too_short: u64,
}

impl Iterator for CramMergedReader {
    type Item = Result<MappedRead, CramInputError>;
    fn next(&mut self) -> Option<Self::Item> { ... }
}
```

`OpenCram` exists at `pub(crate)` visibility — it carries
`noodles::sam::alignment::RecordBuf`, and exposing that publicly
would re-couple our public API to noodles. Tests live in the same
crate, so `pub(crate)` is enough.

Header-level errors surface from `new()`. Per-record errors (malformed
record, out-of-order within a file, duplicate read across files) surface
during iteration.

Filter precedence at each pulled record. Per-check cost is roughly
uniform across the flag-bit and MAPQ filters (one mask+compare or
one byte compare, all on already-decoded record metadata), so the
optimization is on *hit rate*: filters that drop more reads run
first, so subsequent filters skip the most work. Approximate hit
rates in a typical short-read CRAM are noted in parentheses; they
shift with sample and aligner but the relative ordering is stable.

1. Duplicate (~10–30%).
2. Low MAPQ — `min_mapq.is_some_and(|m| read.mapq < m)` (~5–15%).
3. Supplementary (~1–5%).
4. Secondary (~1–5%).
5. Unmapped (~0–5%, often already filtered upstream).
6. QC fail (<1%).
7. Decoded SEQ length < `min_read_length` — structurally last
   because it requires SEQ decoding.

Each flag-bit step is gated by its `drop_*` config bool. A dropped
read increments its bucket in `FilterCounts` and is not yielded; the
iterator pulls again.

### `MappedRead`

Owned, decoupled from noodles. All fields public.

```rust
pub struct MappedRead {
    pub qname: Vec<u8>,
    pub flag: u16,
    pub ref_id: usize,         // index into ContigList
    pub pos: u32,              // 1-based leftmost mapped position
    pub mapq: u8,
    pub cigar: Vec<CigarOp>,
    pub seq: Vec<u8>,          // ACGTN, uppercase
    pub qual: Vec<u8>,         // raw BQ per base (BAQ runs later)
    pub mate_ref_id: Option<usize>,
    pub mate_pos: Option<u32>,
    pub source_file_index: usize,
}

pub enum CigarOp {
    Match(u32), Insertion(u32), Deletion(u32), Skip(u32),
    SoftClip(u32), HardClip(u32), Padding(u32),
    SeqMatch(u32), SeqMismatch(u32),
}
```

### `ContigList`

```rust
pub struct ContigList {
    pub entries: Vec<ContigEntry>,
}

pub struct ContigEntry {
    pub name: String,
    pub length: u32,
    pub md5: Option<[u8; 16]>,   // from @SQ M5 if present
}
```

### `CramInputError`

`thiserror` enum. Variants (one per failure mode, each carrying
enough context to point at the offending file/read):

- `OpenFailed { path, source }` — generic file-open / parse error
  from noodles.
- `MissingFastaIndex { fasta_path }`.
- `UnsupportedCramVersion { path, major, minor }`.
- `NotCoordinateSorted { path, sort_order }`.
- `ContigListMismatch { reference_path, other_path, detail }` —
  detail names the field that disagreed (`name`, `length`, `md5`).
- `FastaContigMismatch { fasta_path, cram_path, detail }`.
- `MissingSampleTag { path, read_group_id }`.
- `MultipleSampleNames { path_a, sm_a, path_b, sm_b }`.
- `OutOfOrderRead { path, qname, prev_pos, this_pos }`.
- `DuplicateReadAcrossFiles { qname, path_a, path_b, ref_id, pos }`.
- `MalformedRecord { path, qname, source }` — generic decode error.
- `Io { path, source: io::Error }`.

The Stage 1 orchestrator (later, in `mod.rs`) uses `anyhow::Result`
and adds `.with_context(|| ...)` around calls into this module.

## Algorithms

### Pre-flight validation (in `new()`)

1. Open the FASTA and load its `.fai`. If `.fai` is missing →
   `MissingFastaIndex`.
2. For each CRAM in input order:
   - Open via `noodles_cram::io::Reader`.
   - Read the file header (CRAM container header) and the embedded SAM
     header.
   - Check CRAM major version == 3. If 4+, return
     `UnsupportedCramVersion`.
   - Check `@HD SO == coordinate`. If absent or different →
     `NotCoordinateSorted`.
   - Extract the `@SQ` list into a local `ContigList`.
   - Extract every `@RG SM`. Empty/missing SM → `MissingSampleTag`.
3. Validate cross-file invariants:
   - All per-file `ContigList`s equal (same length, same name+length+md5
     in the same order). First mismatch → `ContigListMismatch` naming
     the first file as `reference_path` and the offender.
   - All distinct `SM` strings collapse to a single value. Otherwise →
     `MultipleSampleNames` naming the two files that disagree.
4. Validate FASTA agreement: contigs in `.fai` must match the canonical
   `ContigList` in name+length, in the same order. We trust the
   `.fai` per the user decision; we do **not** recompute MD5 from the
   FASTA bytes. If a CRAM `@SQ` carries `M5` and the `.fai` does not,
   we keep the CRAM's MD5 in `ContigEntry.md5` for downstream stages
   (BAQ may want it).
5. Build the canonical `ContigList` and the canonical `sample_name`.
   Store on `CramMergedReader`.

### Decoding via noodles

CRAM is a complex container format (slice-based, reference-driven,
multiple codec versions). We do not decode it ourselves; this slice
sits on top of `noodles_cram`. The integration is small and lives
entirely inside `cram_input.rs` — none of it leaks into the public
API.

**One reference repository, shared by all decoders.** A
`noodles_fasta::Repository` is built once from the indexed FASTA and
passed by reference to every CRAM reader. CRAM blocks decoded
against the same reference share the cache; opening N CRAMs does not
load the FASTA N times.

```rust
let fasta_repository = noodles_fasta::repository::Builder::default()
    .add_adapter(noodles_fasta::indexed_reader::Builder::default()
        .build_from_path(fasta_path)?)
    .build();
```

**One reader per input CRAM.** Each input path is opened via the
noodles builder, configured to use the shared repository:

```rust
let cram_reader = noodles_cram::io::reader::Builder::default()
    .set_reference_sequence_repository(fasta_repository.clone())
    .build_from_path(cram_path)?;
let _header = cram_reader.read_header()?;   // also used for validation
let records_iter = cram_reader.records();   // Iterator<Item = io::Result<RecordBuf>>
```

`RecordBuf` is the *owned* form of a SAM/CRAM record — it does not
borrow from the decoder's internal buffers. That matters for the
merge: peekers can hold heads across multiple `peek()` calls and the
records we hand downstream do not pin any decoder lifetime. Using
the borrowed `Record` form would force lifetime parameters all the
way through the iterator API.

The records iterator is what gets wrapped in `BufferedPeekable` for
the merge loop (see Streaming merge below).

**Per-record conversion** from `RecordBuf` into our owned
`MappedRead` is a one-shot copy in a private helper
`record_buf_to_mapped_read(buf, source_file_index) -> MappedRead`:

| `MappedRead` field   | Source on `RecordBuf`                                   |
|----------------------|---------------------------------------------------------|
| `qname`              | `buf.name()` → bytes, cloned                            |
| `flag`               | `buf.flags().bits()`                                    |
| `ref_id`             | `buf.reference_sequence_id()` (already filtered for unmapped) |
| `pos`                | `buf.alignment_start().unwrap().get() as u32`           |
| `mapq`               | `buf.mapping_quality().map(u8::from).unwrap_or(0)`      |
| `cigar`              | walk `buf.cigar().as_ref()`, convert each op (see below)|
| `seq`                | `buf.sequence().as_ref().to_vec()` (uppercased to ACGTN)|
| `qual`               | `buf.quality_scores().as_ref().to_vec()`                |
| `mate_ref_id`        | `buf.mate_reference_sequence_id()`                      |
| `mate_pos`           | `buf.mate_alignment_start().map(\|p\| p.get() as u32)`    |
| `source_file_index`  | passed in by the merge loop                             |

After conversion, the `MappedRead` is fully owned and disconnected
from the decoder; downstream stages can move it freely.

The exact `RecordBuf` accessor names follow the pinned noodles
release; this table is the contract, the method names are
implementation detail and may change with a noodles upgrade.

### Streaming merge

State on `CramMergedReader`:

- `peekers: Vec<BufferedPeekable<noodles record iterator, ...>>`
  with `buffer_size = 1`.
- `paths: Vec<PathBuf>` (parallel to `peekers`, for error messages).
- `prev_per_file: Vec<Option<(ref_id, pos)>>` for per-file order
  validation.
- `current_start_window: Vec<MappedReadKey>` — `(qname, flag, ref_id,
  pos)` of every read accepted at the current `(ref_id, pos)`. Cleared
  whenever the next pulled read advances past it.
- `config: CramMergedReaderConfig`.
- `filter_counts: FilterCounts`.

Per `next()`:

1. Refill peekers' heads. For each peeker, peek; apply the pre-decode
   filter cascade against the head record in hit-rate order
   (duplicate → low MAPQ → supplementary → secondary → unmapped →
   qc_fail). All checks read only fields already available in the
   alignment record's fixed-position header, so no SEQ/CIGAR decode
   work is paid for a dropped read. If the head fails any filter,
   increment the matching `FilterCounts` bucket, consume the head,
   and re-peek. If the head is `Err`, surface it (mapped to
   `CramInputError`).
2. Choose the smallest surviving head by `(ref_id, pos, file_index)`.
   `file_index` is the deterministic tiebreaker.
3. Validate per-file order: if the chosen peeker's head regresses
   relative to its `prev_per_file` entry, return `OutOfOrderRead`.
4. Maintain the duplicate-detection window: if the chosen `(ref_id,
   pos)` advanced past the window's anchor, clear the window. If the
   chosen read's `(qname, flag, ref_id, pos)` matches anything still
   in the window, return `DuplicateReadAcrossFiles`.
5. Consume the chosen peeker's head. Convert the noodles record into a
   `MappedRead` (clone bytes, decode CIGAR into our enum).
6. Apply the post-decode filter: if
   `config.min_read_length.is_some_and(|min| merged.seq.len() < min as usize)`,
   increment `filter_counts.too_short`, skip, pull next. (Done after
   decode because SEQ length isn't reliably available pre-decode in
   CRAM.)
7. Push the read's key into the window. Update `prev_per_file`.
8. Return `Ok(MappedRead)`.

When all peekers are exhausted, return `None`.

### CIGAR conversion

One-shot conversion from noodles' `Cigar` to `Vec<CigarOp>`. No
allocations beyond the `Vec` itself; CIGARs are short (typically
< 20 ops for short reads).

### Sequence and qualities

`SEQ` decoded by noodles into a byte buffer of ACGTN. Quality scores
likewise — both as raw bytes (Phred 0–93). We clone both into the
`MappedRead`. If `SEQ` is `*` (no sequence), `seq` and `qual` are
empty `Vec`s; the post-decode `min_read_length` filter rejects them.

## TDD order

Tests split into three groups:

- **Pure-type tests.** Touch neither constructor, just exercise types
  in isolation.
- **Group A — logic tests via `from_open_crams`.** Inject `Vec<RecordBuf>`
  as the per-CRAM record stream. No filesystem, no noodles writer, no
  FASTA on disk. Covers every non-I/O concern: merge order, tiebreaker,
  out-of-order detection, dedup window, the entire filter cascade.
  This is where most of the test mass lives.
- **Group B — I/O-glue tests via `new`.** Real synthetic CRAMs and a
  real FASTA on disk. Each test exercises one rule of header
  validation, plus a single end-to-end smoke. Smaller set; covers
  the noodles wiring and `new`'s file-opening / header-parsing logic.

Each numbered step within a group is one implementation unit: write
the failing test first, then make it pass.

### Pure-type tests

#### P1 — `ContigList` + equality
- Derived `PartialEq` with a small twist: when one side's `md5` is
  `None`, it acts as a wildcard against the other's `Some`.
- Test: two equal lists compare equal; mismatched name / length /
  both have `Some(md5)` but differ → not equal; one `None` MD5 vs
  one `Some` MD5 → equal.

### Group A — via `from_open_crams`

These tests build `Vec<OpenCram>` directly. A test helper
`open_cram_from_records(path_for_errors, records: Vec<RecordBuf>)`
wraps a `Vec` in the `Box<dyn Iterator<...>>`. No file is touched.

#### A1 — Single-stream pass-through
- One `OpenCram` with three records at positions 100, 200, 300.
  Iterate; expect three `MappedRead`s in order, fields preserved
  (CIGAR / SEQ / flag / MAPQ / mate fields).

#### A2 — Multi-stream merge order
- Two `OpenCram`s: A with records at 100, 300, 500; B with records at
  150, 200, 400. Iterate; expect 100, 150, 200, 300, 400, 500 with
  `source_file_index` alternating correctly.

#### A3 — Tiebreaker on equal coordinates
- A at (chr1, 100, qname=R1); B at (chr1, 100, qname=R2). Expect A
  first (lower file_index), B second.

#### A4 — Out-of-order within a single stream
- A's records: position 200, then position 100. Iterate; first pull
  yields pos=200; second pull returns `OutOfOrderRead` naming the
  source path.

#### A5 — Duplicate read across streams
- A and B both contain `(qname=R1, flag=0x0, ref_id=0, pos=100)`.
  Iterate; first pull yields the read; second pull returns
  `DuplicateReadAcrossFiles` naming both source paths.

#### A6 — Duplicate-window clears on advance
- A: `(R1, pos=100)`, `(R1, pos=200)`. The second R1 is at a
  different position so it's not a duplicate. Iterate; both yield
  `Ok` (window cleared between positions).

#### A7 — `min_mapq` filter
- A: three records with MAPQ values straddling `DEFAULT_MIN_MAPQ`
  (one well below, one just below, one above). Construct with
  `CramMergedReaderConfig::default()`. Iterate; expect only the
  above-threshold record yielded. `filter_counts().low_mapq == 2`.
- A separate sanity test asserts `DEFAULT_MIN_MAPQ == 20` so a
  future edit to the constant is caught loudly.

#### A8 — `min_mapq = None` disables the MAPQ filter
- Same input as A7, `min_mapq: None`. Iterate; expect all three
  records. `filter_counts().low_mapq == 0`.

#### A9 — Flag-bit drops, one filter at a time
- One stream containing records tagged respectively unmapped
  (0x4), secondary (0x100), supplementary (0x800), qc_fail (0x200),
  duplicate (0x400), plus one clean record.
- For each filter `F`: construct with all defaults except every
  other drop bool set to `false` so only `F` is active. Iterate;
  expect exactly the record tagged with `F` to be missing, and the
  matching `FilterCounts` field == 1.
- Then with full defaults: expect only the clean record yielded;
  every drop bucket == 1.

#### A10 — All flag drops disabled passes everything
- Same input as A9, all `drop_*` bools `false`, `min_mapq: None`,
  `min_read_length: None`. Iterate; expect all six records.

#### A11 — `min_read_length` drops short and empty-SEQ records
- A: long record of length `DEFAULT_MIN_READ_LENGTH + 20`, short
  record of length `DEFAULT_MIN_READ_LENGTH - 5`, record with
  `SEQ = *`. Defaults. Iterate; expect only the long record.
  `filter_counts().too_short == 2`.
- Sanity test asserts `DEFAULT_MIN_READ_LENGTH == 30`.
- A second test with `min_read_length: None` yields all three.

#### A12 — Filter precedence is hit-rate-ordered
- A record flag-marked duplicate *and* unmapped *and* MAPQ 0. With
  all drops on and `min_mapq = Some(20)`, only `duplicate` should
  increment (the highest-hit-rate filter runs first; subsequent
  ones are short-circuited).
- A record unmapped *and* MAPQ 0, but not duplicate: `low_mapq`
  increments (next in order), not `unmapped`.
- A record unmapped only: `unmapped` increments. Together these pin
  the cascade order at duplicate → MAPQ → unmapped without
  enumerating every link.

#### A13 — Empty stream
- One `OpenCram` with an empty record vector. Iterate; expect
  immediate `None`.

#### A14 — Mixed empty + non-empty streams
- A empty, B with two records. Iterate; expect B's two records in
  order.

### Group B — via `new`

These tests build real synthetic CRAMs on disk via
`test_fixtures::build_cram`, plus a real FASTA + `.fai`.

#### B1 — Header parsing on a known-good CRAM
- One CRAM with one contig, one `@RG SM:s1`, two records. Open via
  `new`; assert `sample_name() == "s1"` and `contigs()` matches the
  expected single-contig list.

#### B2 — CRAM 4.x rejection
- Synthesise a header byte string with major version 4 (full file
  not needed — the version check happens before record decode).
- Open via `new`; expect `UnsupportedCramVersion { major: 4, .. }`.

#### B3 — Sort-order rejection
- CRAM with `@HD SO:queryname`. Open via `new`; expect
  `NotCoordinateSorted`.

#### B4 — Sample-tag handling
- One CRAM, one `@RG SM:foo` → `sample_name() == "foo"`.
- One CRAM, one `@RG` with no `SM` → `MissingSampleTag`.
- Two CRAMs with `SM:foo` and `SM:bar` → `MultipleSampleNames`.

#### B5 — Contig-list mismatch across CRAMs
- Two CRAMs differing in `@SQ` order → `ContigListMismatch{detail: "name"}`.
- Differing in length → `ContigListMismatch{detail: "length"}`.
- Conflicting MD5s where both have one → `ContigListMismatch{detail: "md5"}`.

#### B6 — FASTA agreement
- CRAM with contig `chr1:1000`; FASTA `.fai` listing `chr1\t1000` → OK.
- FASTA missing `.fai` → `MissingFastaIndex`.
- FASTA listing a different length → `FastaContigMismatch`.

#### B7 — End-to-end smoke
- One CRAM with three records (positions 100, 200, 300), valid
  reference and `.fai`. Open via `new`, iterate, assert three
  `MappedRead`s with the expected fields. This is the only test
  whose purpose is to confirm the noodles wiring and the
  `new → from_open_crams` plumbing produce the same answer the
  Group A tests get from injected streams.

## Test-fixture helpers

Two helper modules under `src/per_sample_caller/`, both gated
`#[cfg(test)]`. The split mirrors the Group A / Group B test split:
most tests use `record_specs` and never touch the disk; only Group B
needs `cram_files`.

### `record_specs.rs` — for Group A (no I/O)

- `record_spec(name, flag, ref_id, pos, mapq, cigar, seq, qual,
  mate_ref_id, mate_pos) -> RecordBuf` — concise builder for one
  noodles `RecordBuf`. Defaults for fields not relevant to the test.
- `open_cram_from_records(path_for_errors: &str, records:
  Vec<RecordBuf>) -> OpenCram` — wraps a `Vec` in the
  `Box<dyn Iterator<...>> + Send` shape `OpenCram::records` expects.
- `default_contigs() -> ContigList` — a single-contig
  (`chr1`, length 100_000) list for tests that don't care about the
  contig list.

These primitives let a Group A test fit in 5–10 lines:

```rust
let cram_a = open_cram_from_records("a.cram", vec![
    record_spec("R1", 0, 0, 100, 60, "50M", "A".repeat(50), ...),
    record_spec("R2", 0, 0, 300, 60, "50M", ...),
]);
let cram_b = open_cram_from_records("b.cram", vec![
    record_spec("R3", 0, 0, 150, 60, "50M", ...),
]);
let reader = CramMergedReader::from_open_crams(
    vec![cram_a, cram_b],
    default_contigs(),
    "sample".into(),
    CramMergedReaderConfig::default(),
)?;
```

### `cram_files.rs` — for Group B (real I/O)

- `build_fasta(contigs: &[(name, sequence)]) -> (TempDir, PathBuf)` —
  writes a FASTA + `.fai` to a tempdir, returns the path. Keeps
  the `TempDir` alive for the test.
- `build_cram(fasta_path, header_overrides, records: &[RecordBuf]) ->
  (TempDir, PathBuf)` — uses `noodles_cram::io::Writer` to build a
  coordinate-sorted CRAM with the given header tweaks (sort order,
  `@RG SM`, `@SQ` overrides) and records. Returns the CRAM path.

Both helpers stay in `src/per_sample_caller/` so later slices (BAQ,
pileup walker) can reuse them. `record_spec` in particular is
expected to be the workhorse for Group A in every later slice.

## Validation

Run after each step and at the end:

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo test per_sample_caller::cram_input` (targeted)

## Tradeoffs and follow-ups

- **Single-threaded decoders.** Spec calls for one decoder thread per
  CRAM. Deferred: the `Iterator` API stays stable when we later move
  decoding off-thread. Add behind the scenes when profiling shows
  decode-bound runs.
- **Per-record allocations.** Each `MappedRead` clones `qname`, `seq`,
  `qual`, and the CIGAR vector. At 30× coverage this is billions of
  allocations across a human genome. If profiling later shows this is
  hot, we can switch to an arena or an `Arc<[u8]>` shared with the
  decoder buffer. Decoupling from noodles was the user-driven
  decision; we accept the cost.
- **MD5 trust.** We trust `.fai` length and CRAM `@SQ M5` and do not
  recompute MD5 from the FASTA. If a malformed CRAM lies about `M5`,
  BAQ will still get correct sequence from the FASTA — the lie just
  silently goes unflagged. Accepted risk.
- **All per-read filters live here (Option X).** The spec's
  module-layout sketch lists a separate `read_filter.rs`. We absorb
  the entire cascade into the reader: the filter decides which reads
  enter the pipeline, and the natural home for that decision is the
  pipeline entry. It also lets every flag/MAPQ-based drop short-circuit
  before any SEQ/CIGAR decode work. The two filters that need decoded
  data — read length (`min_read_length`) and the "first/last CIGAR op
  is `I`/`D`" indel-only suppression — split: read length stays here
  as a post-decode pass; the CIGAR-shape rule moves to the pileup
  walker, where it semantically belongs (it suppresses an *allele*,
  not a read).
- **Two-constructor split (`new` + `from_open_crams`).** Pushes I/O
  to the program edge per ports-and-adapters reasoning. Public
  callers pay no cost — they call `new(paths, ...)` as before. The
  `pub(crate) from_open_crams` exists so the merge / order / dedup /
  filter logic can be tested with injected `Vec<RecordBuf>` streams,
  without round-tripping through the noodles writer and a tempdir.
  `OpenCram` carries `noodles::sam::alignment::RecordBuf` and is
  `pub(crate)` to avoid leaking noodles into the public API. Going
  further (e.g. abstracting filesystem I/O behind a trait) was
  considered and rejected — single-implementor abstractions add
  generic noise without proportional benefit in Rust.
- **`--region` / `.crai`-driven seeking.** Out of scope. When added,
  `CramMergedReader::with_region(...)` will be a sibling constructor;
  the iterator API stays the same.
