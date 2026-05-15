# Per-sample pileup — `.psp` writer and reader

Implementation plan for Stage 1's serialiser and the inverse reader
that Stages 3–6 will consume. Specified in
[per_sample_pileup_format.md](../specs/per_sample_pileup_format.md);
the high-level contract lives in
[calling_pipeline_architecture.md §"Stage 2"](../specs/calling_pipeline_architecture.md).

The walker (already implemented at
[`src/per_sample_caller/pileup/`](../../src/per_sample_caller/pileup/))
pushes `PileupRecord`s through a `SyncSender<PileupRecord>` channel
(see [pileup/walker.rs:27](../../src/per_sample_caller/pileup/walker.rs#L27)).
The writer is the consumer of that channel. The reader is the inverse:
file path in, `PileupRecord` stream out. **Both halves ship in the same
slice** — they share the wire format and validate each other through
round-trip tests; testing one without the other would force
hand-crafted byte fixtures for every case and isn't worth the cost.

## Scope

What this slice produces:

- A new module tree under [`src/per_sample_caller/psp/`](../../src/per_sample_caller/psp/).
- A `PspWriter` that takes a `PileupRecord` stream and writes a v1.0
  `.psp` to a `Write` sink. Single-pass, streaming, no back-patching.
- A `PspReader` that opens a `.psp` from a `Read + Seek` source and
  yields `PileupRecord`s. Sequential streaming path; a separate
  `PspReader::seek_to_region(chrom_id, start, end)` for tail-index
  random access (`Stage 3` doesn't need this, but `psp dump` will, and
  it's cheap to wire up alongside the index parser).
- A central column-tag registry — the **single source of truth**
  shared by writer, reader, and (later) the `psp_spec_dump` helper
  that emits the spec's column table from the same data.
- Typed errors at module boundary (`PspWriteError`, `PspReadError`);
  `anyhow` context at the orchestrator edge.
- Hard errors for every check listed in spec
  §"Header-binary consistency: required reader checks" and
  §"Block invariants".

What it does **not** include:

- The CLI surface that pumps records from the walker into the writer.
  This plan documents the writer's `Write`-sink API; wiring it to the
  pipeline's actual `--output` flag is a separate, smaller slice
  alongside the eventual `per_sample_caller::cli` module.
- `psp head` / `psp dump` utilities. The reader exposes everything
  they need; these tools are a few hundred lines of presentation on
  top, deferred until the reader has shipped.
- BCF / VCF export. Out of scope per
  [architecture doc §"Why custom binary rather than BCF"](../specs/calling_pipeline_architecture.md).
- Multi-threaded column decompression on the reader side. The reader
  decodes columns sequentially. Per-column rayon parallelism is
  trivially derivable later from the existing code shape (columns are
  independent zstd frames) but is not pursued until profiling shows
  it matters.
- Sidecar indexes, schema evolution beyond v1.0, capability tables.
  All deferred per the spec.

## Cross-references

- [per_sample_pileup_format.md](../specs/per_sample_pileup_format.md)
  — byte-layout spec. Authoritative; this plan defers every wire-level
  decision to it.
- [calling_pipeline_architecture.md §"Stage 2"](../specs/calling_pipeline_architecture.md)
  — high-level contract: per-position records, five scalars, phase
  chain slot lifecycle, columnar zstd, schema-evolution policy.
- [pileup/mod.rs](../../src/per_sample_caller/pileup/mod.rs) — the
  walker's emit types (`PileupRecord`, `AlleleObservation`,
  `AlleleSupportStats`). The writer consumes these; the reader
  reconstructs them.
- [design_principles.md](../specs/design_principles.md) — especially
  §3 "Errors must not pass silently" (load-bearing for every
  consistency check the reader performs) and §4 "Push side effects to
  the edges" (motivates the `Write`-sink / `Read + Seek`-source
  constructors split).

## Dependencies (Cargo.toml)

Add:

- `toml` — TOML v1.0.0 parser and writer. Same crate Cargo itself
  uses; mature, no surprises, tiny.
- `zstd` — Rust bindings to the zstd C library. Already exposed via a
  high-level `encode`/`decode` API; we use it at level 9.
- `xxhash-rust` (feature `xxh3`) — XXH3-64 for the trailer's
  `index_checksum` (spec §"File trailer"). Same hash family zstd uses
  internally; pulling it in directly gives us a sync, no-alloc
  `xxh3_64` for ~5 KB of compiled code.

Already in the manifest: `anyhow`, `thiserror`, `serde` (transitively
via `toml`).

Pin to the latest stable releases at implementation time; record the
versions in the implementation report after the fact.

## Module layout

```
src/per_sample_caller/psp/
    mod.rs              — re-exports, top-level documentation
    errors.rs           — PspWriteError, PspReadError (thiserror)
    registry.rs         — column-tag registry: the single source of
                          truth shared by writer, reader, and (later)
                          psp_spec_dump. v1.0 columns as `const`.
    varint.rs           — LEB128 encode/decode for u64 + zig-zag svarint.
                          Tiny; pulled out so tests can hammer the
                          edge cases (max-width, partial buffer).
    header.rs           — TOML header build, validate, parse.
                          Per-field rule enforcement.
    block.rs            — Block accumulator (writer side), Block decoder
                          (reader side). Column buffer types.
    index.rs            — Block-index entry struct, encode/decode,
                          XXH3-64 checksum.
    trailer.rs          — Trailer struct (32 B), encode/decode.
    writer.rs           — PspWriter public API.
    reader.rs           — PspReader public API.
    tests.rs            — Integration tests; round-trip, fixtures.
```

The split is by concern, not by file size. Each file should be small
enough to read in one pass (cf. existing pileup module, where
`active_read_set.rs`, `cigar_cursor.rs`, etc. each own one idea).

[src/per_sample_caller/mod.rs](../../src/per_sample_caller/mod.rs)
gets one new line: `pub mod psp;` plus a brief sentence in its
crate-level doc pointing at the new module.

The spec's earlier `psp_writer.rs` name is superseded; the module is
`psp`, and the writer entry point inside it is `writer.rs`.

## The column registry is the single source of truth

Per spec §"Required columns in v1.0":

> This table **is** the v1.0 column-tag registry. It is also what a
> v1.0 writer emits verbatim into the file header's `[[column]]` array
> (§"File header — binary schema"). The spec and the in-file schema
> are two views of the same data; readers verify they agree.

In Rust, that means **the registry lives in one place** and every
other piece of code keys off it:

```rust
// src/per_sample_caller/psp/registry.rs

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cardinality { PerRecord, PerAllele }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Shape { Scalar, List, Bytes }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElementType { U8, U16, U32, U64, I32, I64, F32, F64, Varint, Svarint, Bool }

#[derive(Debug, Clone, Copy)]
pub struct ColumnDef {
    pub tag: u16,                       // Spec writes these in hex
                                        // (0x01 etc.); kept as u16
                                        // because the v1.0 registry
                                        // fits and a wider tag would
                                        // mean a major bump.
    pub name: &'static str,             // e.g. "allele-q-sum-log"
    pub cardinality: Cardinality,
    pub shape: Shape,
    pub element_type: Option<ElementType>,  // None iff shape == Bytes
    pub length_column: Option<&'static str>,// Some iff shape == Bytes
    pub required: bool,
    pub description: &'static str,      // Verbatim from the spec
}

/// The v1.0 registry. Order matches spec §"Required columns in v1.0"
/// (tag-ascending; writer emits in this order).
pub const V1_0_COLUMNS: &[ColumnDef] = &[
    ColumnDef { tag: 0x01, name: "delta-pos", /* ... */ },
    ColumnDef { tag: 0x02, name: "n-alleles", /* ... */ },
    // ... twelve entries total
];
```

Three things fall out of this single declaration:

1. **The writer emits the in-file TOML schema by walking
   `V1_0_COLUMNS`** and pretty-printing each `ColumnDef`. No
   parallel source of truth; the writer's output is the registry
   serialised.
2. **The reader validates the in-file schema against `V1_0_COLUMNS`**
   field by field (spec consistency check #6). Any disagreement is a
   hard error — the file claims something different from what the
   spec at its declared `format-version` says.
3. **A future `psp_spec_dump` binary** can emit the spec's table by
   walking the same constant. The spec doc and the writer can't
   drift because the doc is generated.

Internal helper:

```rust
pub fn lookup_by_tag(tag: u16) -> Option<&'static ColumnDef>;
pub fn lookup_by_name(name: &str) -> Option<&'static ColumnDef>;
```

Linear scan is fine — twelve entries.

## Errors

Typed enums at the module boundary, `thiserror`-derived, named after
the specific failure they represent. No `Other(String)` catch-alls.

```rust
// src/per_sample_caller/psp/errors.rs

#[derive(thiserror::Error, Debug)]
pub enum PspWriteError {
    #[error("I/O error writing {context}: {source}")]
    Io { context: &'static str, #[source] source: std::io::Error },

    #[error("invalid TOML field {key:?}: {reason}")]
    InvalidHeaderField { key: String, reason: String },

    #[error("zstd compression failed for column {column}: {source}")]
    Compression { column: &'static str, #[source] source: std::io::Error },

    #[error("non-finite float in column {column} at record {record} allele {allele}: {value}")]
    NonFiniteFloat { column: &'static str, record: usize, allele: usize, value: f64 },

    #[error("allele sequence length {len} for record {record} allele {allele} exceeds cap {cap}")]
    AlleleSeqLenOutOfBounds { record: usize, allele: usize, len: u64, cap: u64 },

    // ... one variant per validation failure the writer can flag.
}

#[derive(thiserror::Error, Debug)]
pub enum PspReadError {
    #[error("I/O error reading {context}: {source}")]
    Io { context: &'static str, #[source] source: std::io::Error },

    #[error("not a .psp file: header magic mismatch")]
    BadMagic,

    #[error("toml_body_length {got} outside allowed range [1, {max}]")]
    BadHeaderLength { got: u64, max: u64 },

    #[error("TOML header parse failed: {source}")]
    HeaderToml { #[source] source: toml::de::Error },

    #[error("invalid TOML field {key:?}: {reason}")]
    InvalidHeaderField { key: String, reason: String },

    #[error("header length prefix and sentinel disagree at offset {offset}")]
    SentinelMismatch { offset: u64 },

    #[error("in-file schema disagrees with v1.0 registry on column {name}: {field} = {got} (expected {expected})")]
    SchemaMismatch { name: String, field: &'static str, got: String, expected: String },

    #[error("required column {name} (tag {tag:#x}) is not recognised by this reader")]
    UnknownRequiredColumn { name: String, tag: u16 },

    #[error("block {block} has out-of-order or duplicate column tags in manifest")]
    BlockManifestOrder { block: u64 },

    #[error("block {block} column {column}: uncompressed_len {got} disagrees with schema-predicted {expected}")]
    UncompressedLenMismatch { block: u64, column: String, got: u64, expected: u64 },

    #[error("block {block}: non-finite float in column {column} at record {record} allele {allele}")]
    NonFiniteFloat { block: u64, column: String, record: usize, allele: usize },

    #[error("block index checksum mismatch: stored {stored:#010x}, computed {computed:#010x}")]
    IndexChecksum { stored: u32, computed: u32 },

    #[error("block {block}: phase-chain active set inconsistency at record {record}: {detail}")]
    PhaseChainConsistency { block: u64, record: usize, detail: String },

    // ... etc.
}
```

The error names are deliberately spec-aligned: when a check from
§"Header-binary consistency" fails, its error variant says so.

## Writer public API

```rust
// src/per_sample_caller/psp/writer.rs

pub struct PspWriter<W: Write> {
    // private: sink, header (cached so we can re-emit / validate),
    //          current block accumulator, block index in progress,
    //          n_blocks written, zstd encoder pool, position
    //          state for delta-pos enforcement, phase-chain snapshot
    //          tracker for active_chain_slots_at_block_start.
}

/// Required up-front metadata the writer needs to assemble the
/// header. Everything Stage 1 already knows by the time the walker
/// starts emitting records — no field here depends on per-record
/// content.
pub struct PspWriterHeader {
    pub sample: String,
    pub reference: String,
    pub chromosomes: Vec<ChromosomeEntry>,    // matches input @SQ
    pub writer: WriterProvenance,
}

pub struct ChromosomeEntry {
    pub name: String,
    pub length: u32,
    pub md5: [u8; 16],                        // 16 raw bytes; stringified
                                              // to 32 hex chars on emit.
}

pub struct WriterProvenance {
    pub tool: String,                         // e.g. "join_per_sample_vcfs"
    pub version: String,                      // e.g. env!("CARGO_PKG_VERSION")
    pub subcommand: String,                   // e.g. "per-sample"
    pub input_crams: Vec<String>,             // basenames only
    pub input_fasta: String,                  // basename only
    pub parameters: BTreeMap<String, ParameterValue>,
}

pub enum ParameterValue {
    Integer(i64),
    Float(f64),
    Boolean(bool),
    String(String),
}

impl<W: Write> PspWriter<W> {
    /// Opens the writer and emits the file header immediately. The
    /// header is fully determined by `header` plus the v1.0 column
    /// registry — no record-dependent state goes into the header.
    /// Validation against the per-field rules (spec §"Per-field
    /// character set") happens here; any failure is a write-time
    /// error before any block is written.
    pub fn new(sink: W, header: PspWriterHeader) -> Result<Self, PspWriteError>;

    /// Append one record. Buffered into the current block; the
    /// writer flushes a block when the projected uncompressed
    /// payload reaches the configured target (default 16 MiB) or
    /// when the next record's chrom_id differs from the current
    /// block's.
    ///
    /// Returns the number of bytes written to `sink` as a side effect
    /// of any flush triggered by this call (0 if no flush).
    pub fn write_record(&mut self, record: &PileupRecord) -> Result<u64, PspWriteError>;

    /// Finalise the file: flush the current block if any records are
    /// buffered, write the block index, write the trailer. After
    /// `finish` returns, the sink is at a complete `.psp` boundary;
    /// any subsequent write to `sink` would corrupt the file.
    ///
    /// Consumes the writer to make double-finish unrepresentable.
    pub fn finish(self) -> Result<W, PspWriteError>;
}
```

### Constants

```rust
/// Target uncompressed bytes per block. Spec §"Block sizing". Not on
/// the CLI; format-internal.
pub const DEFAULT_TARGET_BLOCK_BYTES: u64 = 16 * 1024 * 1024;

/// zstd compression level. Spec §"Compression". Not on the CLI;
/// format-internal.
pub const ZSTD_COMPRESSION_LEVEL: i32 = 9;

/// Hard upper bound for the TOML header body length. Spec §"File
/// header / Layout" / Constraints.
pub const MAX_HEADER_BODY_BYTES: u64 = 1_048_547;  // 1 MiB - 29 B framing

/// Per-allele sequence length cap. Spec §"Required columns in v1.0"
/// row 0x03 (allele_seq_len).
pub const MAX_ALLELE_SEQ_LEN: u64 = 10_000;
```

### Writer-side invariants and validation

The writer enforces every invariant the reader will later check. The
goal is "if `PspWriter::finish()` returned `Ok`, the file is valid";
the reader is for files written elsewhere (other tool versions, on
other machines) — but every check the reader runs, the writer also
runs at emission time, so a buggy writer fails its own
`write_record` call rather than silently producing a corrupt file
that only a stricter reader catches months later.

Enforced at `write_record` time:

- `record.alleles` non-empty (REF always present); enforces
  `n_total_alleles ≥ n_records`.
- Every `allele.seq.len()` in `[1, MAX_ALLELE_SEQ_LEN]`.
- Every byte of every `allele.seq` is in `{A, C, G, T, N}`.
- `record.alleles[0].support.q_sum` and every alt's `q_sum` is
  finite (`f64::is_finite()`).
- `record.chrom_id` is within the header's chromosome table.
- `record.pos ≥ 1` and `record.pos ≤ chromosome.length`.
- `record.pos > previous_record.pos` when same chrom; strictly
  increasing.
- Strict monotonic non-decreasing in `chrom_id` (chromosome
  boundaries can only move forward).

Enforced at block-flush time:

- `n_records ≥ 1` (the writer never flushes an empty block).
- `active_chain_slots_at_block_start` snapshot is correctly computed
  from the running active set the writer maintains.
- Column manifest entries are emitted in tag-ascending order.

Enforced at `finish` time:

- The last block (if any) has been flushed.
- The block index XXH3-64 is computed over the emitted index bytes.
- The trailer's `index_offset + index_byte_length` equals the
  trailer's own offset.

## Reader public API

```rust
// src/per_sample_caller/psp/reader.rs

pub struct PspReader<R: Read + Seek> {
    // private: source, parsed header, parsed block index (kept in
    //          memory; small), block iterator state, sequential-mode
    //          phase-chain active set running cursor.
}

/// Parsed and validated file header. Returned by the reader so
/// callers don't have to re-parse it themselves.
#[derive(Debug, Clone)]
pub struct ParsedHeader {
    pub format_version: (u16, u16),
    pub sample: String,
    pub reference: String,
    pub created: chrono::DateTime<chrono::Utc>,
    pub chromosomes: Vec<ParsedChromosome>,
    pub writer: ParsedWriter,
    pub columns: Vec<ParsedColumn>,                 // in-file binary schema
}

impl<R: Read + Seek> PspReader<R> {
    /// Open: tail-reads the trailer, decodes and checksum-verifies
    /// the block index, then reads and validates the header. All
    /// header-binary consistency checks that don't require block
    /// payloads (spec §"Header-binary consistency" items 1–6) run
    /// here. Block-level checks run as blocks are read.
    pub fn new(source: R) -> Result<Self, PspReadError>;

    pub fn header(&self) -> &ParsedHeader;

    /// Sequential streaming: yields records in genomic order across
    /// all blocks. Internally walks the block index. Maintains the
    /// inter-block phase-chain continuity check (consistency check
    /// #11) by carrying the running active set across block
    /// boundaries.
    pub fn records(&mut self) -> RecordsIter<'_, R>;

    /// Random-access region query. Uses the block index to find
    /// blocks overlapping `[start, end]` on `chrom_id` and yields
    /// records inside. Cannot perform the inter-block continuity
    /// check (no preceding-block active set in memory) — relies on
    /// each block's `active_chain_slots_at_block_start` snapshot,
    /// which is exactly what that snapshot exists for (spec
    /// §"Phase-chain state across blocks").
    pub fn region_records(
        &mut self,
        chrom_id: u32,
        start: u32,
        end: u32,
    ) -> RecordsIter<'_, R>;
}

pub struct RecordsIter<'r, R: Read + Seek> { /* private */ }

impl<'r, R: Read + Seek> Iterator for RecordsIter<'r, R> {
    type Item = Result<PileupRecord, PspReadError>;
    fn next(&mut self) -> Option<Self::Item>;
}
```

### Reader-side validation order

Spec §"Header-binary consistency" enumerates the checks; the reader
runs them in this order, eager-fail at every step:

1. (in `new`) magic, length prefix range, header parse, per-field
   rules, schema-vs-registry agreement, sentinel cross-check.
2. (in `new`) trailer magic, trailer arithmetic, block-index
   XXH3-64 checksum.
3. (in `new`) block-index sanity: `n_blocks` matches index entries;
   index entries are coordinate-monotonic.
4. (per block, lazily) block header well-formedness, block
   invariants, per-block manifest ordering, `uncompressed_len`
   a-priori prediction (spec consistency check #7), decompressed
   size match, structural correctness during decode, NaN check,
   phase-chain active-set consistency, inter-block continuity (for
   sequential reads only).

The phase-chain `active_chain_slots_at_block_start` snapshot is
trusted on random-access path and cross-checked on sequential path
(spec check #11). The implementation:

- For sequential `records()`: maintain a `BTreeSet<SlotId>`
  carried across blocks. At each block boundary, compare the set
  to the new block's snapshot; mismatch → hard error.
- For random-access `region_records()`: at each block, initialise
  the running set from the block's snapshot directly. No
  cross-check possible.

## Block layout — write path

The writer maintains one `BlockAccumulator` at a time. Per record,
the accumulator:

1. Validates the record (see §"Writer-side invariants").
2. Computes `delta_pos = record.pos - last_pos` (0 for first
   record of block).
3. Pushes one entry to each per-record column buffer: `delta_pos`,
   `n_alleles`, `new_chain_slots`, `expired_chain_slots`.
4. Pushes one entry per allele to each per-allele column buffer:
   `allele_seq_len`, `allele_seq` (bytes), the five scalar columns,
   `allele_chain_slots`.
5. Updates the running active-slot set (`+= new_chains`,
   `-= expired_chains`).
6. Accumulates the projected uncompressed-payload total.

When the projected total crosses `DEFAULT_TARGET_BLOCK_BYTES`, or
when the next record's `chrom_id` differs from the current block's,
the accumulator flushes:

1. Writes the uncompressed block header (varint stream).
2. For each column, in tag-ascending order:
   - Serialises the buffered values into the column's wire form
     (delegates to per-shape helpers in `block.rs`).
   - zstd-compresses the bytes at level 9 into a `Vec<u8>`.
   - Appends `(column_tag, compressed_len, uncompressed_len)` to
     an in-memory manifest.
3. Emits the manifest entries (varint stream).
4. Emits each compressed payload back-to-back.
5. Records the block's `(chrom_id, first_pos, last_pos, n_records,
   block_offset)` into the in-memory block-index buffer.
6. Resets the accumulator; updates the snapshot of currently-active
   slots for the next block's `active_chain_slots_at_block_start`.

The accumulator is bounded in memory by the target block size plus
some slack for the column buffers themselves; on a 16 MiB target
that's tens of MiB peak, acceptable.

### Per-column wire encoding (writer side, mirrors spec §"Encoding rules per cardinality/shape")

```rust
// src/per_sample_caller/psp/block.rs

fn encode_per_record_scalar<T: WireScalar>(values: &[T], out: &mut Vec<u8>);
fn encode_per_allele_scalar<T: WireScalar>(values: &[T], out: &mut Vec<u8>);
fn encode_per_record_list<T: WireScalar>(lists: &[&[T]], out: &mut Vec<u8>);
fn encode_per_allele_list<T: WireScalar>(lists: &[&[T]], out: &mut Vec<u8>);
fn encode_per_allele_bytes(values: &[&[u8]], lengths_out: &mut Vec<u8>, bytes_out: &mut Vec<u8>);
```

`WireScalar` is a small private trait implemented for `u8`/`u16`/
`u32`/`u64`/`i32`/`i64`/`f32`/`f64`/`bool` (fixed-width LE) plus
varint/svarint wrappers. Floats serialise via `to_le_bytes`; the
finite-value check happens before this call.

### Block index and trailer

```rust
// src/per_sample_caller/psp/index.rs

pub struct BlockIndexEntry {
    pub chrom_id: u32,
    pub first_pos: u32,
    pub last_pos: u32,
    pub n_records: u32,
    pub block_offset: u64,
}

pub fn encode_index(entries: &[BlockIndexEntry]) -> Vec<u8>;
pub fn decode_index(bytes: &[u8], expected_n_blocks: u64) -> Result<Vec<BlockIndexEntry>, PspReadError>;
pub fn xxh3_64_trunc32(bytes: &[u8]) -> u32;        // truncate to low 32 bits
```

```rust
// src/per_sample_caller/psp/trailer.rs

pub const TRAILER_BYTES: usize = 32;
pub const TRAILER_MAGIC: [u8; 4] = *b"PSPE";

pub struct Trailer {
    pub index_offset: u64,
    pub index_byte_length: u64,
    pub n_blocks: u64,
    pub index_checksum: u32,
}

pub fn encode_trailer(t: &Trailer) -> [u8; TRAILER_BYTES];
pub fn decode_trailer(bytes: &[u8; TRAILER_BYTES]) -> Result<Trailer, PspReadError>;
```

### Header build (writer side)

```rust
// src/per_sample_caller/psp/header.rs

pub fn build_header_toml(
    header: &PspWriterHeader,
    created: chrono::DateTime<chrono::Utc>,
) -> Result<Vec<u8>, PspWriteError>;
```

Implementation: use `toml::Value` / `toml_edit::DocumentMut` to
construct the document programmatically (key order matches the spec's
example), then serialise. Validate every field against the per-field
rules **before** serialising — we don't trust the user-supplied
`PspWriterHeader` strings to be ASCII-printable, sample-name-clean,
etc.

Writer order (matches the spec example):

1. Top-level keys (`format-version`, `sample`, `reference`, `created`).
2. `[[chromosome]]` array, in `@SQ` order.
3. `[writer]` table with provenance.
4. `[writer.parameters]` sub-table.
5. `[[column]]` array, walking `registry::V1_0_COLUMNS`.

The serialised TOML body, the length prefix, the magic, the sentinel
are then emitted in one `sink.write_all` if possible (the whole header
is small).

### Header parse (reader side)

```rust
pub fn parse_header(toml_body: &str) -> Result<ParsedHeader, PspReadError>;
pub fn validate_header(parsed: &ParsedHeader) -> Result<(), PspReadError>;
```

Two-step:

1. `toml::from_str` into a `Deserialize`-derived `RawHeader` struct.
   Parse failures map to `PspReadError::HeaderToml`.
2. Walk every field, run its per-field rule (spec §"Per-field
   character set"). Map failures to `PspReadError::InvalidHeaderField`
   with the key and a human-readable reason.
3. Cross-check `[[column]]` entries against `V1_0_COLUMNS` field by
   field (consistency check #6).

## Block layout — read path

Mirror image. Per block:

1. `read_exact` the next bytes through the manifest (uncompressed
   varint stream).
2. Validate block invariants:
   - `n_records ≥ 1`.
   - `n_total_alleles ≥ n_records`.
   - `chrom_id` in range.
   - `first_pos` in `[1, chromosome.length]`.
   - `active_chain_slots_at_block_start`: ascending, distinct,
     empty if first block on the chromosome.
   - Manifest is tag-ascending with no duplicates.
3. For each declared column, predict `uncompressed_len` from the
   schema where possible (consistency check #7), reject on
   mismatch.
4. `read_exact` each `compressed_len`-byte payload; zstd-decode it;
   confirm output size matches `uncompressed_len`.
5. Decode column values into in-memory buffers (mirroring the
   writer's encoders).
6. Walk the per-record and per-allele buffers in tandem,
   reconstructing `PileupRecord`s. Phase-chain active-set
   consistency and float-finite checks fire here.

## Tests

Three layers, each addressing a specific failure mode without
depending on the others.

### Group V — wire primitives (no I/O, no zstd)

Unit tests in `varint.rs`, `index.rs`, `trailer.rs`.

#### V1 — varint round-trips
- Encode/decode `0`, `1`, `127`, `128`, `16383`, `16384`,
  `u32::MAX as u64`, `u64::MAX`. Each must round-trip.
- Truncated buffer (one byte short of a multi-byte value) returns
  an `IncompleteVarint` error.

#### V2 — svarint round-trips
- `0`, `±1`, `±63`, `±64`, `i32::MIN as i64`, `i64::MIN`.

#### V3 — trailer round-trip
- Build a trailer, encode, decode, assert equality.
- Bad magic byte → `BadMagic`.
- `index_offset + index_byte_length` overflow → `BadTrailer`.

#### V4 — index round-trip with checksum
- Build N entries (1, 100, 10_000), encode, compute XXH3-64-trunc,
  decode, assert all entries equal and `n_blocks` matches.
- Flip one byte in the encoded index → `IndexChecksum`.

### Group H — header build/parse/validate (no zstd, no blocks)

Unit tests in `header.rs`.

#### H1 — round-trip a minimal valid header
- Build with one chromosome, one input CRAM, two parameters.
- Build → bytes; parse → struct; equality.

#### H2 — every per-field rule, both sides
For each rule in spec §"Per-field character set":
- `sample` with a non-ASCII byte → writer error
  `InvalidHeaderField{key: "sample"}`; reader error on a hand-
  crafted file that smuggles the bad value through.
- `chromosome.name` with `<` (illegal per `@SQ SN` regex) →
  writer error; reader error.
- `chromosome.md5` with 31 hex chars → both errors.
- `chromosome.length = 0` → both errors.
- `writer.input-crams[*]` containing `/` → both errors.
- `column.element-type = "u128"` (not in the closed enum) → both
  errors.
- `column.shape = "bytes"` without `length-column` → both errors.
- `column.shape = "scalar"` *with* `length-column` → both errors.

#### H3 — schema-vs-registry disagreement
- Hand-craft a header whose `allele-q-sum-log` says
  `element-type = "f32"`. Reader rejects with `SchemaMismatch`
  naming the column.
- Hand-craft a header whose `new-chain-slots` says
  `cardinality = "per-allele"`. Same.
- Hand-craft a header whose `[[column]]` array omits `allele-seq`
  (required). Reader rejects with `MissingRequiredColumn`.

#### H4 — unknown top-level keys are tolerated
- Add `future-feature = "hi"` at top level. Parser accepts; reader
  proceeds.

#### H5 — TOML body length prefix and sentinel
- Build a file with a deliberately wrong `toml_body_length` (one
  byte too short). Reader gets a `SentinelMismatch`.
- Build a file whose body legally contains a multi-line string
  with embedded `\n` (impossible under the v1.0 per-field rules,
  but achievable with a misbehaving writer). Length prefix
  navigates around it correctly; sentinel cross-check passes.
- A `toml_body_length` of 0 or >1_048_547 → `BadHeaderLength`.

### Group B — block round-trip (zstd in play)

Integration tests in `tests.rs`.

#### B1 — single record, single allele, single block
- Build a `PileupRecord` with chrom=0, pos=100, one REF allele,
  five scalars set to known values, no chains.
- Write → bytes; read → records; assert one record matches.

#### B2 — multi-allele record
- Two records on chr1: pos=100 with REF=A, ALT=C (SNP);
  pos=200 with REF=ACGT, ALT1=A (deletion), ALT2=ACGTAG (insertion).
- Round-trip; assert allele sequences and scalars survive.

#### B3 — multi-block file
- 100k records, target block size set to 64 KiB so several blocks
  emit. Round-trip; assert record count and last record match.

#### B4 — multi-chromosome
- 10 records on chr1, then 10 on chr2. Two blocks (chromosome
  boundary forces flush). Round-trip; assert per-chrom record
  counts.

#### B5 — phase chains across block boundary
- Construct a stream where slot 7 is opened in block N's last
  record and consumed in block N+1's first record. Writer's
  `active_chain_slots_at_block_start` snapshot must list `7` on
  block N+1; reader's sequential continuity check must pass.
- Mutate the test fixture to corrupt the snapshot (drop `7`).
  Reader's continuity check fires with
  `PhaseChainConsistency`.

#### B6 — empty file (zero blocks)
- Write a file with no `write_record` calls. Round-trip; reader
  yields zero records; trailer arithmetic still validates;
  index has zero entries.

#### B7 — random-access region query
- 1000 records across chr1:1–100k, chr2:1–100k. Open with
  `region_records(chrom_id=0, start=50_000, end=60_000)`; assert
  only matching records returned.
- Same with a region wholly outside any block (chrom_id=0,
  start=200_000, end=300_000): zero records, no error.

### Group F — failure-mode catches

Negative tests; each constructs a deliberately-corrupt file and
asserts the matching `PspReadError` variant.

#### F1 — header magic flipped
File starts with `PSQ\n` instead of `PSP\n`. Reader returns
`BadMagic` immediately.

#### F2 — truncated mid-block
Truncate a known-good file mid-payload. Reader returns
`Io { context: "block payload" }`.

#### F3 — torn zstd frame
Flip one byte inside a compressed column payload. Reader returns
`Compression { source: <zstd decode error> }`.

#### F4 — corrupted block index
Flip one byte inside the index region. Reader returns
`IndexChecksum`.

#### F5 — NaN injected into `allele-q-sum-log`
Hand-craft a block whose `allele-q-sum-log` payload contains a
NaN bit pattern. Reader returns `NonFiniteFloat`.

#### F6 — out-of-order column tags
Hand-craft a block manifest with `[0x02, 0x01]`. Reader returns
`BlockManifestOrder`.

#### F7 — `uncompressed_len` lying about fixed-width-scalar size
Hand-craft a manifest where `allele-obs-count`'s
`uncompressed_len = n_total_alleles × 5` instead of `× 4`. Reader
returns `UncompressedLenMismatch` *before* decompression.

#### F8 — first block of chromosome carries a non-empty snapshot
Set `active_chain_slots_at_block_start = [3]` on block 0. Reader
returns `PhaseChainConsistency` naming the violation.

## Test fixtures

Two helper modules under `src/per_sample_caller/psp/`, both gated
`#[cfg(test)]`.

### `test_records.rs` — synthetic `PileupRecord` builders

- `record(chrom_id, pos, alleles: &[(seq, support, chain_slots)],
   new_chains, expired_chains) -> PileupRecord` — concise builder.
  Sets sensible defaults for fields not relevant to the test.
- `support_stats(num_obs, q_sum, fwd, placed_left, placed_start) ->
   AlleleSupportStats`.
- `default_header(n_chroms: usize) -> PspWriterHeader` —
  single-sample, generic contig names, all-default parameters.

Most tests fit in 5–10 lines:

```rust
let records = vec![
    record(0, 100, &[(b"A", support_stats(4, -1.0, 2, 1, 0), vec![1])],
           vec![1], vec![]),
    record(0, 200, &[(b"C", support_stats(3, -0.8, 1, 0, 0), vec![1])],
           vec![], vec![1]),
];

let mut writer = PspWriter::new(Cursor::new(Vec::new()), default_header(1))?;
for r in &records {
    writer.write_record(r)?;
}
let bytes = writer.finish()?.into_inner();

let mut reader = PspReader::new(Cursor::new(bytes))?;
let read_back: Vec<_> = reader.records().collect::<Result<_, _>>()?;
assert_eq!(records, read_back);
```

### `corrupt_fixtures.rs` — byte-level mutation helpers

- `corrupt_byte(bytes: &mut [u8], offset: usize, new_value: u8)` —
  documented assertion + mutation for Group F tests.
- `replace_pattern(bytes: &mut Vec<u8>, find: &[u8], replace: &[u8])`
  — used by F1.
- `find_block_index_start(bytes: &[u8]) -> usize` — reads the
  trailer's `index_offset` to find a block-index byte to flip in
  F4.
- `find_allele_q_sum_log_payload(bytes: &[u8], block: usize) ->
   std::ops::Range<usize>` — used by F5 to inject a NaN into a
  known location.

These live close to the tests so anyone writing a new negative test
has a small library of mutations to draw from.

## Validation

Run after each step and at the end:

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo test per_sample_caller::psp` (targeted)

A representative round-trip test (Group B) is the closest thing this
slice has to a smoke test: writer → reader → equality over a few
thousand records exercises every column type, every per-record /
per-allele cardinality, every phase-chain marker, the index path,
and the trailer.

## Implementation order

Suggested slicing to keep each commit's diff focused:

1. **`registry.rs` + `errors.rs` + `varint.rs`.** No dependencies on
   anything beyond `std`. Tests: Group V (V1, V2).
2. **`trailer.rs` + `index.rs`.** Adds `xxhash-rust`. Tests: V3, V4.
3. **`header.rs`.** Adds `toml`. Tests: Group H entirely. Group H
   is large; landing it here gives the writer and reader's hardest
   shared concern (TOML schema agreement) full coverage early.
4. **`block.rs`.** Adds `zstd`. Implements the per-shape encoder /
   decoder primitives. Unit tests for each (`encode_per_record_scalar`
   round-trips through `decode_per_record_scalar`, etc.).
5. **`writer.rs`.** Wires registry + header + block + index + trailer
   into a streaming writer. No new dependencies.
6. **`reader.rs`.** Symmetric. Adds the sequential and random-access
   iterator entry points. Tests: Groups B and F.

Each step ends green (`cargo test` passes) before the next begins.

## Tradeoffs and follow-ups

- **Block index in memory, not streaming.** The writer accumulates
  index entries in `Vec<BlockIndexEntry>` until `finish()`. On a 5×
  WGS that's ~3 k entries (spec §"Size estimate") — under 100 KB
  peak. If a future use case needs the writer to handle pathological
  block counts (millions), spill-to-disk is easy to add without
  changing the on-disk format.
- **No parallel column decompression.** Columns inside a block are
  independent zstd frames, so `rayon::scope` over them is trivial.
  Deferred to a follow-up; the sequential decode is simpler to get
  right first, and benchmarks haven't shown it to be the bottleneck.
- **No mmap.** Reader uses ordinary `pread` via `Read + Seek`.
  `mmap` is plausibly faster for random-access region queries but
  brings file-descriptor / safety / signal-handling subtleties that
  the project doesn't need before Stage 3 actually exercises the
  random-access path on large files. Revisit if/when profiling
  shows it.
- **`PspWriter::write_record` runs the walker's consumer half on the
  caller's thread.** The walker today emits over a `SyncSender`;
  pulling from the receiver and feeding `write_record` is a thin
  wrapper around `for record in receiver { writer.write_record(&record)?; }`.
  Moving the writer to its own thread (so walker → channel → writer
  is true producer/consumer parallelism) is one line of `thread::spawn`
  in the orchestration layer, and lives there rather than in this
  module.
- **`chrono` vs `time` for `created`.** TOML's offset date-time maps
  cleanly to `chrono::DateTime<Utc>`. The `time` crate is the modern
  alternative; sticking with `chrono` because the rest of the
  ecosystem this project depends on uses it.
- **`zstd` crate vs `zstd-safe`.** The high-level `zstd` crate
  wraps `zstd-safe`. We use the high-level one; the encode/decode
  functions are sufficient and a streaming API isn't needed
  per-column (each column fits in memory at 16 MiB target).
- **Column emission strategy.** The writer assembles each column's
  buffer in full before zstd-compressing. An alternative would be to
  stream column bytes into the zstd encoder as records arrive. For
  the per-record scalar columns this is straightforward; for
  per-allele lists and bytes columns it's more involved. The
  buffered approach uses one extra block's worth of memory (~16 MiB)
  but keeps the per-shape encoders straightforward; revisit if peak
  memory becomes the bottleneck.
- **`PileupRecord` ↔ wire mapping is direct.** The walker's
  in-memory record type matches the spec's logical record one-to-one,
  so the writer is "validate then serialise" and the reader is
  "parse then validate". No intermediate IR.
