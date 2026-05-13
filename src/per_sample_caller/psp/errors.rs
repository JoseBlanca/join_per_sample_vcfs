//! Typed error enums for the `.psp` writer and reader.
//!
//! Two top-level enums — `PspWriteError` for the writer side,
//! `PspReadError` for the reader — match the producer / consumer
//! split in the spec's §"Header-binary consistency: required reader
//! checks" and §"Block invariants". Sub-enums (e.g. `VarintError`)
//! describe failures intrinsic to a primitive and get wrapped into
//! the top-level enums with context at the call site.
//!
//! Each variant is one concrete failure mode. No `Other(String)`
//! catch-alls — callers that need to react to a specific failure
//! match on the variant; tests assert on the exact variant rather
//! than stringly-comparing messages.
//!
//! New variants land as later modules need them; this file currently
//! covers what the registry / varint / wire-primitive slice requires.

use thiserror::Error;

/// Failures intrinsic to LEB128 / zig-zag-LEB128 decoding. Always
/// wrapped by [`PspReadError::Varint`] (the read side is the only
/// side that decodes — the writer constructs varints by construction
/// and cannot fail mid-encode).
#[derive(Error, Debug, Clone, Copy, PartialEq, Eq)]
pub enum VarintError {
    /// Buffer ran out mid-varint: a continuation bit was set on the
    /// last byte available, but no more bytes were present.
    #[error("varint truncated: continuation bit set on the final available byte")]
    Truncated,

    /// More than 10 LEB128 continuation bytes were consumed. The
    /// spec caps the varint width at 10 bytes (covers `u64`); a
    /// longer encoding cannot represent a valid `u64` and indicates
    /// either corruption or a writer bug.
    #[error("varint overflow: continuation bytes exceeded the 10-byte cap")]
    Overflow,
}

/// Failures the `WireScalar` decoders can produce. Generalises
/// [`VarintError`] with fixed-width-only failure modes
/// (buffer-too-short, illegal bool byte).
#[derive(Error, Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScalarDecodeError {
    #[error("buffer truncated mid-decode")]
    Truncated,

    #[error("varint exceeded 10-byte cap")]
    VarintOverflow,

    /// `bool` columns use a single byte; only `0` (false) and `1`
    /// (true) are legal.
    #[error("invalid bool byte {0:#04x}; expected 0 or 1")]
    InvalidBool(u8),
}

impl From<VarintError> for ScalarDecodeError {
    fn from(e: VarintError) -> Self {
        match e {
            VarintError::Truncated => Self::Truncated,
            VarintError::Overflow => Self::VarintOverflow,
        }
    }
}

/// Errors the `.psp` reader can emit. Variants land as the
/// corresponding slice ships; the registry / varint slice
/// contributes the [`Self::Varint`] case.
#[derive(Error, Debug)]
pub enum PspReadError {
    #[error("varint decode failed{}: {source}", match context {
        Some(c) => format!(" while reading {c}"),
        None => String::new(),
    })]
    Varint {
        /// What was being decoded when the failure occurred — e.g.
        /// `"block header n_records"`. Optional so callers can pass
        /// the bare primitive error up where context is added higher.
        context: Option<&'static str>,
        #[source]
        source: VarintError,
    },

    #[error("I/O error reading {context}: {source}")]
    Io {
        context: &'static str,
        #[source]
        source: std::io::Error,
    },

    /// Trailer's magic bytes did not equal `PSPE`. Either the file is
    /// not a `.psp`, has been truncated past the trailer's start, or
    /// is corrupted at the tail.
    #[error("trailer magic mismatch: got {got:02x?}, expected {expected:02x?}")]
    BadTrailerMagic { got: [u8; 4], expected: [u8; 4] },

    /// The XXH3-64 hash (low 32 bits) computed over the index bytes
    /// does not match the value stored in the trailer. Spec
    /// §"Header-binary consistency" rule 11 / spec §"File trailer".
    /// Hard error — the index region is the one file region not
    /// covered by a zstd frame checksum, and silently trusting a
    /// corrupt index would redirect every subsequent block read to
    /// wrong bytes.
    #[error("block index checksum mismatch: stored {stored:#010x}, computed {computed:#010x}")]
    IndexChecksum { stored: u32, computed: u32 },

    /// Failed to decode a varint inside a block-index entry.
    #[error("block index entry {entry} field {field}: {source}")]
    IndexEntryDecode {
        entry: usize,
        field: &'static str,
        #[source]
        source: VarintError,
    },

    /// A block-index varint decoded successfully but exceeds the
    /// representable range of the field's owning type (`u32` for the
    /// genomic-coordinate-ish fields).
    #[error("block index entry {entry} field {field}: value {value} exceeds the field's u32 range")]
    IndexFieldOverflow {
        entry: usize,
        field: &'static str,
        value: u64,
    },

    /// The block-index buffer ran out before a fixed-width field
    /// (currently only `block_offset: u64`) could be read.
    #[error("block index ran out of bytes while reading entry {entry} field {field}")]
    IndexTruncated { entry: usize, field: &'static str },

    /// The block-index buffer carried extra bytes past the last
    /// entry. Either the trailer's `index_byte_length` is wrong, the
    /// writer over-counted, or there is corruption between the index
    /// and the trailer.
    #[error(
        "block index has {trailing_bytes} trailing bytes after the last entry; index_byte_length and the writer disagree"
    )]
    IndexTrailingBytes { trailing_bytes: usize },

    /// The decoded entry count disagrees with the trailer's
    /// `n_blocks` field. Either the index is short of entries (and a
    /// varint truncation should have been caught earlier) or the
    /// trailer's count is wrong.
    #[error("block index has {got} entries, trailer claims {expected}")]
    IndexEntryCountMismatch { got: usize, expected: u64 },

    /// The first four bytes of the file are not `PSP\n` — either the
    /// file is not a `.psp` or the head magic is corrupt.
    #[error("file head magic mismatch: got {got:02x?}, expected {expected:02x?}")]
    BadHeadMagic { got: [u8; 4], expected: [u8; 4] },

    /// `toml_body_length` (the 8-byte u64 length prefix between the
    /// magic and the TOML body) is outside the legal `[1, 1048547]`
    /// range. Validated before any buffer is allocated so a malicious
    /// or corrupt file cannot drive a large allocation off this
    /// field.
    #[error("toml_body_length {got} outside allowed range [{min}, {max}]")]
    BadHeaderLength { got: u64, min: u64, max: u64 },

    /// The TOML body parsed successfully as TOML but its content
    /// violates a per-field rule (§"Per-field character set" or a
    /// schema invariant). `key` names the offending TOML key; `reason`
    /// is a human-readable explanation of which rule failed.
    #[error("invalid header field {key:?}: {reason}")]
    InvalidHeaderField { key: String, reason: String },

    /// The TOML body failed to parse as TOML at all. Wraps the
    /// underlying `toml::de::Error` so the line/column information
    /// survives.
    #[error("TOML header parse failed: {source}")]
    HeaderToml {
        #[source]
        source: toml::de::Error,
    },

    /// The 17 bytes immediately after the TOML body do not equal the
    /// `---END-HEADER---\n` sentinel. Either the writer is buggy
    /// (`toml_body_length` and the actual TOML body do not match) or
    /// the file is corrupt.
    #[error(
        "header sentinel mismatch at byte offset {offset}: length prefix and sentinel disagree"
    )]
    SentinelMismatch { offset: u64 },

    /// A column declared in the file's `[[column]]` array has
    /// `required = true` but the reader's column-tag registry does
    /// not recognise its `name` (or `tag`). The file's contract is
    /// "this required column must be honoured"; an unknown reader
    /// must refuse rather than proceed without it.
    #[error("required column {name:?} (tag {tag:#x}) is not recognised by this reader")]
    UnknownRequiredColumn { name: String, tag: u16 },

    /// A column the v1.0 registry requires is absent from the file's
    /// `[[column]]` array.
    #[error("required column {name:?} (tag {tag:#x}) is missing from the file's [[column]] array")]
    MissingRequiredColumn { name: String, tag: u16 },

    /// A column whose `name`/`tag` the reader knows disagrees with
    /// the registry on a structural field (cardinality, shape,
    /// element-type, length-column, or required). The file claims
    /// something different from what the spec at its declared
    /// `format-version` says.
    #[error(
        "in-file schema disagrees with v1.0 registry on column {name:?}: {field} = {got:?} (expected {expected:?})"
    )]
    SchemaMismatch {
        name: String,
        field: &'static str,
        got: String,
        expected: String,
    },

    /// `format-version` is syntactically `MAJOR.MINOR` but the
    /// declared major version exceeds what this reader supports.
    /// Hard error — a higher-major reader's content may rest on
    /// layout assumptions this reader cannot make.
    #[error(
        "file was written by format-version {file_major}.{file_minor}; this reader supports up to {reader_major}.{reader_minor}"
    )]
    UnsupportedFormatVersion {
        file_major: u16,
        file_minor: u16,
        reader_major: u16,
        reader_minor: u16,
    },

    /// A column payload ran out before the expected entry count was
    /// decoded. `decoded` entries had been produced when the cursor
    /// hit the end of the buffer; `expected` is the cardinality-
    /// derived count.
    #[error("column {column:?} payload truncated: decoded {decoded} entries, expected {expected}")]
    ColumnTruncated {
        column: String,
        decoded: usize,
        expected: usize,
    },

    /// A column payload carried bytes past the last expected entry —
    /// the writer wrote more than the cardinality demands.
    #[error("column {column:?} payload has {trailing} trailing bytes after the last entry")]
    ColumnTrailingBytes { column: String, trailing: usize },

    /// zstd decompression failed for a column payload. Almost always
    /// means corruption inside the compressed bytes (a single-byte
    /// flip catches here via zstd's built-in frame checksum).
    #[error("zstd decompression failed for {context}: {source}")]
    Zstd {
        context: String,
        #[source]
        source: std::io::Error,
    },

    /// A column's decompressed payload size did not match the
    /// per-block manifest's `uncompressed_len`. Catches torn frames
    /// whose internal checksum happens to pass — should not happen
    /// in practice but the check is cheap.
    #[error("column {column:?} decompressed to {got} bytes, manifest claimed {expected}")]
    UncompressedLenMismatch {
        column: String,
        got: u64,
        expected: u64,
    },

    /// A column's manifest-claimed `uncompressed_len` disagrees with
    /// what the schema predicts a-priori for fixed-width or
    /// length-column-bounded shapes. Caught *before* decompression
    /// (spec §"Header-binary consistency" check #7).
    #[error("column {column:?}: uncompressed_len {got} disagrees with schema-predicted {expected}")]
    UncompressedLenSchemaMismatch {
        column: String,
        got: u64,
        expected: u64,
    },

    /// A block header field could not be decoded; usually a varint
    /// truncation inside the header byte stream.
    #[error("block header field {field}: {source}")]
    BlockHeaderField {
        field: &'static str,
        #[source]
        source: VarintError,
    },

    /// A single element inside a column payload failed to decode.
    /// `entry` is the per-cardinality index (record index for
    /// per-record columns, allele index for per-allele).
    #[error("column {column:?} entry {entry}: {source}")]
    ColumnElementDecode {
        column: String,
        entry: usize,
        #[source]
        source: ScalarDecodeError,
    },

    /// A block-header structural invariant was violated. The
    /// `reason` text names the invariant; the writer caught one of
    /// the per-block invariants from spec §"Block sizing / Block
    /// invariants" before emitting, or the reader caught a
    /// corrupted block on the way in.
    #[error("block header invariant violation: {reason}")]
    BlockHeaderInvariant { reason: String },
}

/// Errors the `.psp` writer can emit. Variants land as the
/// corresponding slice ships.
#[derive(Error, Debug)]
pub enum PspWriteError {
    #[error("I/O error writing {context}: {source}")]
    Io {
        context: &'static str,
        #[source]
        source: std::io::Error,
    },

    /// A field the writer was about to emit violates one of the spec's
    /// per-field rules (§"Per-field character set"). Writer-side
    /// pre-emit validation catches the failure before any byte is
    /// written.
    #[error("invalid header field {key:?}: {reason}")]
    InvalidHeaderField { key: String, reason: String },

    /// The TOML serializer failed (typically because a field cannot be
    /// represented in TOML at all — should not happen for our schema,
    /// but the failure surface exists). Carries the underlying
    /// message as a string because `toml::ser::Error`'s exact type is
    /// not part of the `toml` crate's public surface.
    #[error("TOML header serialisation failed: {message}")]
    HeaderToml { message: String },

    /// The serialised TOML body exceeded the 1 MiB-minus-framing cap
    /// (or somehow came out empty). Spec §"File header / Layout"
    /// pins the bounds.
    #[error("header TOML body length {got} outside allowed range [{min}, {max}]")]
    HeaderBodyTooLarge { got: u64, min: u64, max: u64 },

    /// A record handed to `write_record` violates a spec invariant.
    /// `reason` names which rule failed; `record_index` is the
    /// zero-based position in the writer's input stream so the
    /// error message can point at the offending input.
    #[error("invalid record at index {record_index}: {reason}")]
    InvalidRecord { record_index: u64, reason: String },

    /// A record's `chrom_id` is outside the writer's chromosome
    /// table (declared via [`WriterHeader::chromosomes`] in `new`).
    #[error(
        "record at index {record_index} has chrom_id {chrom_id}; header declared only {n_chroms} chromosomes"
    )]
    UnknownChromId {
        record_index: u64,
        chrom_id: u32,
        n_chroms: u32,
    },

    /// A record's `pos` is `0` (1-based) or exceeds the contig's
    /// declared length.
    #[error(
        "record at index {record_index}: pos {pos} out of [1, {chrom_length}] for chrom_id {chrom_id}"
    )]
    PosOutOfRange {
        record_index: u64,
        chrom_id: u32,
        pos: u32,
        chrom_length: u32,
    },

    /// Two consecutive records on the same chromosome have
    /// non-strictly-increasing positions, or `chrom_id` decreased
    /// from the previous record.
    #[error(
        "record at index {record_index} regresses: ({this_chrom}, {this_pos}) follows ({prev_chrom}, {prev_pos})"
    )]
    OutOfOrderRecord {
        record_index: u64,
        prev_chrom: u32,
        prev_pos: u32,
        this_chrom: u32,
        this_pos: u32,
    },

    /// A phase-chain marker is inconsistent with the running
    /// active-slot set: `new_chains` references a slot already
    /// active, or `expired_chains` references a slot not active.
    #[error("record at index {record_index}: phase-chain marker inconsistency: {reason}")]
    PhaseChainMarkerInconsistency { record_index: u64, reason: String },

    /// `finish` has been called but the writer-side encoded block
    /// header rejected our own data — bug in the writer.
    #[error("block emission failed at index {block_index}: {source}")]
    BlockEmission {
        block_index: u64,
        #[source]
        source: PspReadError,
    },
}
