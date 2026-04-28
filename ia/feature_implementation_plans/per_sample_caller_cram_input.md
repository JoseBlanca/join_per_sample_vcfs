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
- A `CramReaderConfig` struct with `Default` providing the spec's
  defaults (`min_mapq = 20`, all flag-based drops on, `min_read_length
  = 1`).
- Enough of `MergedRead` to carry every field downstream stages need
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

### `CramReaderConfig`

```rust
pub struct CramReaderConfig {
    /// Reads with MAPQ strictly below this are dropped. 0 disables.
    pub min_mapq: u8,
    /// Reads whose decoded SEQ length is strictly below this are dropped.
    /// 0 disables.
    pub min_read_length: u32,
    pub drop_unmapped: bool,        // flag & 0x4
    pub drop_secondary: bool,       // flag & 0x100
    pub drop_supplementary: bool,   // flag & 0x800
    pub drop_qc_fail: bool,         // flag & 0x200
    pub drop_duplicate: bool,       // flag & 0x400
}

impl Default for CramReaderConfig {
    fn default() -> Self {
        Self {
            min_mapq: DEFAULT_MIN_MAPQ,
            min_read_length: DEFAULT_MIN_READ_LENGTH,
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
CramMergedReader::new(&crams, fa, CramReaderConfig::default())?;

// Override one field
CramMergedReader::new(&crams, fa, CramReaderConfig {
    min_mapq: 30,
    ..Default::default()
})?;
```

### `CramMergedReader`

```rust
pub struct CramMergedReader {
    // private: peekers, contig list, sample name, dup-detection window,
    //          config, per-file order tracking, filter counters.
}

impl CramMergedReader {
    /// Open every CRAM, validate headers against each other and the
    /// FASTA, apply the per-read filter cascade per `config`, and
    /// prepare a merged coordinate-sorted stream of surviving reads.
    pub fn new(
        crams: &[PathBuf],
        fasta: &Path,
        config: CramReaderConfig,
    ) -> Result<Self, CramInputError>;

    pub fn sample_name(&self) -> &str;
    pub fn contigs(&self) -> &ContigList;

    /// Counts of reads dropped by each filter. Logged to stderr by the
    /// orchestrator at end-of-run per spec §"Errors" / Warnings.
    pub fn filter_counts(&self) -> &FilterCounts;
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
    type Item = Result<MergedRead, CramInputError>;
    fn next(&mut self) -> Option<Self::Item> { ... }
}
```

Header-level errors surface from `new()`. Per-record errors (malformed
record, out-of-order within a file, duplicate read across files) surface
during iteration.

Filter precedence at each pulled record (cheapest checks first to
short-circuit decode work):

1. Flag-bit drops in this order: unmapped, secondary, supplementary,
   qc_fail, duplicate. Each enabled by its config bool.
2. MAPQ < `min_mapq`.
3. Decoded SEQ length < `min_read_length`.

A dropped read increments its bucket in `FilterCounts` and is not
yielded; the iterator pulls again.

### `MergedRead`

Owned, decoupled from noodles. All fields public.

```rust
pub struct MergedRead {
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

### Streaming merge

State on `CramMergedReader`:

- `peekers: Vec<BufferedPeekable<noodles record iterator, ...>>`
  with `buffer_size = 1`.
- `paths: Vec<PathBuf>` (parallel to `peekers`, for error messages).
- `prev_per_file: Vec<Option<(ref_id, pos)>>` for per-file order
  validation.
- `current_start_window: Vec<MergedReadKey>` — `(qname, flag, ref_id,
  pos)` of every read accepted at the current `(ref_id, pos)`. Cleared
  whenever the next pulled read advances past it.
- `config: CramReaderConfig`.
- `filter_counts: FilterCounts`.

Per `next()`:

1. Refill peekers' heads. For each peeker, peek; apply the filter
   cascade against the head record (flag bits, then MAPQ — both are
   in the alignment record's fixed-position header, so no SEQ/CIGAR
   decode work is paid for a dropped read). If the head fails any
   filter, increment the matching `FilterCounts` bucket, consume the
   head, and re-peek. If the head is `Err`, surface it (mapped to
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
   `MergedRead` (clone bytes, decode CIGAR into our enum).
6. Apply the post-decode filter: if `merged.seq.len() <
   config.min_read_length`, increment `filter_counts.too_short`,
   skip, pull next. (Done after decode because SEQ length isn't
   reliably available pre-decode in CRAM.)
7. Push the read's key into the window. Update `prev_per_file`.
8. Return `Ok(MergedRead)`.

When all peekers are exhausted, return `None`.

### CIGAR conversion

One-shot conversion from noodles' `Cigar` to `Vec<CigarOp>`. No
allocations beyond the `Vec` itself; CIGARs are short (typically
< 20 ops for short reads).

### Sequence and qualities

`SEQ` decoded by noodles into a byte buffer of ACGTN. Quality scores
likewise — both as raw bytes (Phred 0–93). We clone both into the
`MergedRead`. If `SEQ` is `*` (no sequence), `seq` and `qual` are
empty `Vec`s; the post-decode `min_read_length` filter rejects them.

## TDD order

Each numbered step is one implementation unit: write the failing test
first, then make it pass.

### Step 1 — `ContigList` + equality
- New type with derived `PartialEq` (modulo MD5 absence: `None == Some`
  on either side is treated as a match per spec, the full match check
  is only when both sides have `Some`).
- Test: two equal lists compare equal; mismatched name / length / both
  have `Some(md5)` but differ → not equal.

### Step 2 — Header parsing on a known-good CRAM
- Build a tiny CRAM in-memory using noodles (one contig, two reads),
  read its header through our `read_header_metadata` helper.
- Test: returns expected contig list, sort order, sample name.

### Step 3 — CRAM 4.x rejection
- Synthesise a fake header byte string with major version 4 (we don't
  need a fully valid file; the version check happens before record
  decode).
- Test: `UnsupportedCramVersion` returned with the right major.

### Step 4 — Sort-order rejection
- Build a CRAM with `@HD SO:queryname`.
- Test: `NotCoordinateSorted`.

### Step 5 — Sample-tag handling
- Build a CRAM with one `@RG` and `SM:foo`, expect `Ok("foo")`.
- Build a CRAM with one `@RG` missing SM, expect `MissingSampleTag`.
- Build two CRAMs with `SM:foo` and `SM:bar`, run `new()`, expect
  `MultipleSampleNames`.

### Step 6 — Contig-list mismatch
- Two CRAMs with different `@SQ` orders → `ContigListMismatch{detail:
  "name"}`.
- Two CRAMs differing on length → `ContigListMismatch{detail: "length"}`.
- Two CRAMs with conflicting MD5s where both have one →
  `ContigListMismatch{detail: "md5"}`.

### Step 7 — FASTA agreement
- CRAM with contig `chr1:1000`; FASTA `.fai` listing `chr1\t1000`. OK.
- FASTA missing `.fai` → `MissingFastaIndex`.
- FASTA listing a different length → `FastaContigMismatch`.

### Step 8 — Single-CRAM streaming
- One CRAM with three reads at positions 100, 200, 300. Iterate;
  expect three `MergedRead`s in order, with the right CIGAR / SEQ /
  flag / mate fields.

### Step 9 — Multi-CRAM merge order
- Two CRAMs: file A with reads at 100, 300, 500; file B with reads
  at 150, 200, 400.
- Iterate; expect 100, 150, 200, 300, 400, 500 with `source_file_index`
  alternating correctly.

### Step 10 — Tiebreaker on equal coordinates
- File A read at (chr1, 100) qname=R1; file B read at (chr1, 100)
  qname=R2. Expect file A first (lower file_index), file B second.

### Step 11 — Out-of-order within a single file
- Build a CRAM that violates coordinate sort (read at 200 before read
  at 100; we have to construct this manually since noodles' writer
  may sort for us). Iterate; expect `OutOfOrderRead` on the second
  pull.

### Step 12 — Duplicate read across files
- File A and file B both contain the same `(qname=R1, flag=0x0,
  ref_id=0, pos=100)`. Iterate; first pull yields the read, second
  pull returns `DuplicateReadAcrossFiles` naming both files.

### Step 13 — Duplicate-window clears on advance
- File A: `(R1, pos=100)`, `(R1, pos=200)`. The second R1 is at a
  different position so it's not a duplicate. Iterate; both should
  yield `Ok` (the window cleared between positions).

### Step 14 — `min_mapq` filter
- File A: three reads with MAPQ values straddling
  `DEFAULT_MIN_MAPQ` (one well below, one just below, one above).
  Construct with `CramReaderConfig::default()`. Iterate; expect only
  the above-threshold read yielded. `filter_counts().low_mapq == 2`.
- A separate sanity test asserts `DEFAULT_MIN_MAPQ == 20` so a future
  edit to the constant is caught loudly.

### Step 15 — `min_mapq = 0` disables the MAPQ filter
- Same input as Step 14, `min_mapq: 0`. Iterate; expect all three
  reads. `filter_counts().low_mapq == 0`.

### Step 16 — Flag-bit drops, one filter at a time
- One file containing reads tagged respectively unmapped (0x4),
  secondary (0x100), supplementary (0x800), qc_fail (0x200),
  duplicate (0x400), plus one clean read.
- For each filter `F`: construct with all defaults except every other
  drop bool set to `false` so only `F` is active. Iterate; expect
  exactly the read tagged with `F` to be missing, and the matching
  `FilterCounts` field to be 1.
- Then with full defaults: expect only the clean read yielded; counts
  reflect each drop bucket = 1.

### Step 17 — All flag drops disabled passes everything
- Same input as Step 16, all `drop_*` bools `false`, `min_mapq: 0`,
  `min_read_length: 0`. Iterate; expect all six reads.

### Step 18 — `min_read_length` drops short and empty-SEQ reads
- File A: a long read of length `DEFAULT_MIN_READ_LENGTH + 20`, a
  short read of length `DEFAULT_MIN_READ_LENGTH - 5`, and a read
  with `SEQ = *` (no sequence). Construct with
  `CramReaderConfig::default()`. Iterate; expect only the long read.
  `filter_counts().too_short == 2`.
- A sanity test asserts `DEFAULT_MIN_READ_LENGTH == 30`.
- A second test with `min_read_length: 0` yields all three reads.

### Step 19 — Filter precedence is cheap-first
- A single read is unmapped *and* has MAPQ 0 *and* is also flag-marked
  duplicate. With all drops on, only `unmapped` should increment
  (confirms short-circuit ordering, not just final outcome).

### Step 20 — Empty CRAM
- A CRAM with a header but no records. Iterate; expect immediate `None`.

### Step 21 — Mixed empty + non-empty
- File A empty, file B with two reads. Iterate; expect file B's two
  reads in order.

## Test-fixture helpers

Tests need synthetic CRAMs. Build a small helper module
`src/per_sample_caller/test_fixtures.rs` (gated `#[cfg(test)]`):

- `build_fasta(contigs: &[(name, sequence)]) -> (TempDir, PathBuf)` —
  writes a FASTA + `.fai` to a tempdir, returns the path. Keeps the
  TempDir alive for the test.
- `build_cram(fasta_path, header_overrides, reads) -> (TempDir,
  PathBuf)` — uses `noodles_cram::io::Writer` to build a coordinate-
  sorted CRAM with the given header tweaks (sort order, RG/SM, SQ
  list overrides) and reads. Returns the path.
- `read_spec(name, flag, ref_id, pos, mapq, cigar, seq, qual,
  mate_ref_id, mate_pos)` — concise constructor for record specs.

These helpers are also useful for the later slices (filter, BAQ,
pileup walker), so they go in a shared location now.

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
- **Per-record allocations.** Each `MergedRead` clones `qname`, `seq`,
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
- **`--region` / `.crai`-driven seeking.** Out of scope. When added,
  `CramMergedReader::with_region(...)` will be a sibling constructor;
  the iterator API stays the same.
