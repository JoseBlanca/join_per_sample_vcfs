//! Typed error enums for the `.psp` writer and reader.
//!
//! Two top-level enums — `PspWriteError` for the writer side,
//! `PspReadError` for the reader — match the producer / consumer
//! split in the spec's §"Header-binary consistency: required reader
//! checks" and §"Block invariants". Sub-enums describe failures
//! intrinsic to a primitive (varint decoding, fixed-width scalar
//! decoding) or to a cross-cut category (block-header invariants,
//! record-validation rules, phase-chain marker rules), and get
//! wrapped into the top-level enums with context at the call site.
//!
//! Each variant is one concrete failure mode. No `Other(String)`
//! catch-alls — callers that need to react to a specific failure
//! match on the variant; tests assert on the exact variant rather
//! than stringly-comparing messages.
//!
//! Variants are organised by spec section: header framing
//! ([`PspReadError::BadHeadMagic`], [`PspReadError::SentinelMismatch`], …),
//! block-header invariants ([`PspReadError::BlockHeaderField`],
//! [`PspReadError::BlockHeaderInvariant`]), column payloads
//! ([`PspReadError::ColumnTruncated`], …), and trailer/index
//! integrity ([`PspReadError::BadTrailerMagic`],
//! [`PspReadError::IndexChecksum`], …). Both top-level enums carry
//! `#[non_exhaustive]` so future variants can land as additive
//! changes.

use thiserror::Error;

/// Failures intrinsic to LEB128 / zig-zag-LEB128 decoding. Always
/// wrapped by a more specific [`PspReadError`] variant
/// (`IndexEntryDecode`, `BlockHeaderField`, `ColumnElementDecode`)
/// that carries the file-region context.
#[non_exhaustive]
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
#[non_exhaustive]
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

/// Reasons a block-header value can violate the structural
/// invariants the spec pins. Shared between the writer (which
/// checks before emission) and the reader (which checks on decode);
/// each side wraps the same kind in its own top-level error variant
/// (see [`PspReadError::BlockHeaderInvariant`] and
/// [`PspWriteError::BlockEmission`]). M8.
#[non_exhaustive]
#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum BlockHeaderInvariantKind {
    #[error("n_records must be >= 1 (empty blocks are forbidden)")]
    EmptyBlock,

    /// Retained for wire/format compatibility but **no longer raised**
    /// by the generic block-header validator (architecture §10):
    /// `n_total_alleles >= n_records` is a SNP semantic, not a
    /// container one — the SSR schema legitimately stores loci with
    /// zero per-record entries. SNP integrity is now enforced where it
    /// belongs: the writer rejects zero-allele records
    /// (`validate_record`) and the reader cross-checks the per-allele
    /// column counts against `n_total_alleles` (see
    /// `block::validate_block_header_invariants`). This variant is kept
    /// (not deleted) because `BlockHeaderInvariantKind` is part of the
    /// stable error surface; do not re-wire it into the generic
    /// validator without re-introducing the SNP-only coupling.
    #[error(
        "n_total_alleles {n_total_alleles} < n_records {n_records} \
         (every record has at least one allele)"
    )]
    AllelesLessThanRecords {
        n_records: u32,
        n_total_alleles: u32,
    },

    #[error("manifest tags not strictly ascending: {prev:#x} then {next:#x}")]
    ManifestTagsNotAscending { prev: u16, next: u16 },

    #[error("{field}: value {value} exceeds u32 range")]
    FieldExceedsU32 { field: &'static str, value: u64 },

    #[error("{field}: value {value} exceeds u16 range")]
    FieldExceedsU16 { field: &'static str, value: u64 },

    /// The decoded `n-alleles` column sums to a value different
    /// from the block header's `n_total_alleles`. Spec
    /// §"Per-block manifest agreement": any over- or under-run is
    /// a hard error. Without this check, an over-run drives an
    /// index-out-of-bounds panic on per-allele indexing, and an
    /// under-run silently truncates the trailing alleles.
    #[error(
        "n-alleles column sums to {sum_n_alleles}, block header \
         declares n_total_alleles = {n_total_alleles}"
    )]
    NAllelesSumMismatch {
        n_total_alleles: u32,
        sum_n_alleles: u64,
    },

    /// A record's decoded `n-alleles` entry is `0`. The writer forbids
    /// zero-allele records (the REF allele is always present —
    /// [`InvalidRecordKind::ZeroAlleles`]), so a `0` here marks a
    /// corrupt block. Without this check the per-allele CSR offsets
    /// alias the next record's range: a trailing zero-allele record
    /// indexes one past the offsets array (panic), an interior one
    /// silently mis-attributes the following record's allele span.
    #[error("record {record_index} has zero alleles (every record has at least one allele)")]
    ZeroAlleleRecord { record_index: u32 },
}

/// Rules a record handed to `write_record` can violate. Carried as
/// the `kind` of [`PspWriteError::InvalidRecord`]. Mi10.
#[non_exhaustive]
#[derive(Error, Debug, Clone, PartialEq)]
pub enum InvalidRecordKind {
    #[error("record has zero alleles (REF always present)")]
    ZeroAlleles,

    #[error("allele {allele_index} sequence length {length} outside [1, {max}]")]
    AlleleSeqLen {
        allele_index: usize,
        length: usize,
        max: u64,
    },

    #[error(
        "allele {allele_index} byte {byte_offset} = {byte:#04x} \
         (only A/C/G/T/N allowed)"
    )]
    InvalidAlleleByte {
        allele_index: usize,
        byte_offset: usize,
        byte: u8,
    },

    #[error("allele {allele_index} q_sum is non-finite ({q_sum})")]
    NonFiniteQSum { allele_index: usize, q_sum: f64 },

    #[error("allele {allele_index} chain_ids not strictly ascending")]
    AlleleChainIdsNotAscending { allele_index: usize },
}

/// Wrapper around `toml::ser::Error` so the writer's error chain
/// keeps the `#[source]` link without leaking the foreign type's
/// shape into our public API. Mi12.
#[derive(Debug)]
pub struct TomlSerError(toml::ser::Error);

impl TomlSerError {
    pub(crate) fn new(e: toml::ser::Error) -> Self {
        Self(e)
    }
}

impl std::fmt::Display for TomlSerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl std::error::Error for TomlSerError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.0.source()
    }
}

/// Errors the `.psp` reader can emit. New variants are added per
/// spec section as each slice lands.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum PspReadError {
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

    /// The header's `kind` schema-family tag (architecture §10.3) is
    /// not one this reader knows. The `kind` selects which column
    /// registry the `[[column]]` array is cross-checked against
    /// (`"snp"` today; `"ssr"` later); an unknown value means the file
    /// was written by a newer/foreign producer and cannot be decoded.
    #[error("unknown .psp kind {kind:?} (this reader knows: {known})")]
    UnknownKind { kind: String, known: &'static str },

    /// The typed reader was instantiated for a schema (`expected`) that
    /// does not match the file's header `kind` (`found`) — e.g.
    /// `records_of::<SnpKind>()` on an `.ssr.psp` file. Decoding columns
    /// under the wrong schema would yield garbage records, so the
    /// iterator refuses (architecture §10.3). Distinct from
    /// [`Self::UnknownKind`], which is raised at header parse when the
    /// `kind` itself is unrecognised; this is raised on the typed read
    /// path when the caller-chosen schema disagrees with a *known* kind.
    #[error("requested schema {expected:?} does not match the file's kind {found:?}")]
    KindMismatch {
        expected: &'static str,
        found: String,
    },

    /// A per-block structural invariant specific to a schema's columnar
    /// layout was violated on decode — a corrupt or foreign block, not
    /// an I/O fault. `context` names the violated invariant (e.g. the
    /// SSR parallel-CSR columns disagree on row boundaries, or the
    /// per-record grouping count over-runs the decoded entry total).
    /// Kept separate from [`Self::Io`] so callers can distinguish a
    /// torn read (retryable) from deterministic corruption (fatal).
    #[error("block structural invariant violated: {context}")]
    BlockStructureInvalid { context: &'static str },

    /// The SSR per-record grouping column (`n-spanning`) sums to a
    /// different total than the decoded per-profile entry count
    /// (`n_total_alleles`). `got < expected` would silently drop
    /// trailing profiles; `got > expected` would over-run the CSR
    /// offsets. Either way the block is malformed/foreign.
    #[error(
        "ssr profile count mismatch: n-spanning sums to {got} but the block declares {expected} \
         per-profile entries"
    )]
    SsrProfileCountMismatch { expected: u32, got: u64 },

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

    /// A column the v1.0 registry requires is absent from a
    /// per-block manifest. Symmetric with
    /// [`Self::MissingRequiredColumn`] but checked at block-decode
    /// time (the TOML `[[column]]` array is validated at header
    /// parse; the per-block manifest is an independent layer and
    /// must be cross-checked separately — B1).
    #[error(
        "required column {name:?} (tag {tag:#x}) is missing from a \
         per-block manifest"
    )]
    MissingRequiredColumnInManifest { name: String, tag: u16 },

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

    /// The block-header grow-on-incomplete read loop hit its size
    /// cap (currently 64 KiB) without `decode_block_header`
    /// succeeding. Either the header is genuinely larger than the
    /// cap (writer bug or unsupported future format) or the file's
    /// bytes never form a valid header. Distinct from
    /// [`Self::BlockHeaderField`] with `VarintError::Overflow`,
    /// which signals a single-varint overflow at the decoder layer.
    /// M3.
    #[error(
        "block header exceeds {cap}-byte read cap (read {consumed} \
         bytes without a successful decode)"
    )]
    BlockHeaderExceedsCap { cap: usize, consumed: usize },

    /// EOF was hit mid-decode while reading more bytes into the
    /// block-header buffer. Distinct from
    /// [`Self::BlockHeaderField`]'s field-specific truncation —
    /// here we ran out of source bytes before any decode attempt
    /// could discover where in the header it failed. `offset` is
    /// the start of this block on the source; `consumed` is the
    /// number of bytes successfully pulled before EOF. M4.
    #[error(
        "block header truncated at offset {offset} after reading \
         {consumed} bytes"
    )]
    BlockHeaderTruncated { offset: u64, consumed: usize },

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

    /// A block-header structural invariant was violated; the reader
    /// caught a corrupted block on the way in. The cause carries the
    /// specific invariant; see [`BlockHeaderInvariantKind`].
    #[error("block header invariant violation")]
    BlockHeaderInvariant {
        #[source]
        kind: BlockHeaderInvariantKind,
    },

    /// A `block_offset` in the decoded block index is out of order
    /// or out of range. `PspReader` rejects this before any block
    /// is read so a tampered index cannot redirect block reads to
    /// arbitrary file offsets.
    #[error(
        "block index entry {block}: block_offset {offset} \
         outside expected range [{min}, {max})"
    )]
    BlockIndexOffsetInvalid {
        block: usize,
        offset: u64,
        min: u64,
        max: u64,
    },

    /// A block-index entry's `chrom_id` is outside the file's
    /// declared chromosome table. The header parse must succeed
    /// before this check runs.
    #[error(
        "block index entry {block}: chrom_id {chrom_id} \
         outside the file's chromosome table (n_chroms = {n_chroms})"
    )]
    BlockIndexChromOutOfRange {
        block: usize,
        chrom_id: u32,
        n_chroms: u32,
    },

    /// A block-index entry's `first_pos` or `last_pos` lies outside
    /// `[1, chromosomes[chrom_id].length]`, or `first_pos >
    /// last_pos`.
    #[error(
        "block index entry {block}: pos {pos} outside [1, {chromosome_length}] \
         for chrom_id {chrom_id}"
    )]
    BlockIndexPosOutOfRange {
        block: usize,
        chrom_id: u32,
        pos: u32,
        chromosome_length: u32,
    },

    /// An `f32` or `f64` column carried a non-finite value (NaN or
    /// ±∞). The writer enforces finite-only on its side; this
    /// fires if a hand-crafted or torn file slips a non-finite
    /// past the zstd frame check. Spec check #9.
    #[error("column {column:?} entry {entry}: non-finite float {value}")]
    NonFiniteFloat {
        column: String,
        entry: usize,
        value: f64,
    },
}

/// Errors the `.psp` writer can emit. Variants land as the
/// corresponding slice ships.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum PspWriteError {
    /// I/O error during a write. Optional `block_index` and
    /// `column_tag` fields identify the file region that was being
    /// written when the failure occurred — populated for per-block
    /// and per-column writes, `None` for file-level writes (header,
    /// trailer, final flush). Mi9.
    #[error(
        "I/O error writing {context}{}: {source}",
        format_block_column(*block_index, *column_tag),
    )]
    Io {
        context: &'static str,
        block_index: Option<u64>,
        column_tag: Option<u16>,
        #[source]
        source: std::io::Error,
    },

    /// A field the writer was about to emit violates one of the spec's
    /// per-field rules (§"Per-field character set"). Writer-side
    /// pre-emit validation catches the failure before any byte is
    /// written.
    #[error("invalid header field {key:?}: {reason}")]
    InvalidHeaderField { key: String, reason: String },

    /// The TOML serializer failed (typically because a field cannot
    /// be represented in TOML at all — should not happen for our
    /// schema, but the failure surface exists). The wrapper newtype
    /// preserves the `Error::source()` chain.
    #[error("TOML header serialisation failed")]
    HeaderToml {
        #[source]
        source: TomlSerError,
    },

    /// The serialised TOML body exceeded the 1 MiB-minus-framing cap
    /// (or somehow came out empty). Spec §"File header / Layout"
    /// pins the bounds.
    #[error("header TOML body length {got} outside allowed range [{min}, {max}]")]
    HeaderBodyTooLarge { got: u64, min: u64, max: u64 },

    /// A record handed to `write_record` violates a spec invariant.
    /// `record_index` is the zero-based position in the writer's
    /// input stream so the error message can point at the offending
    /// input; `kind` names which rule failed.
    #[error("invalid record at index {record_index}")]
    InvalidRecord {
        record_index: u64,
        #[source]
        kind: InvalidRecordKind,
    },

    /// A record's `chrom_id` is outside the writer's chromosome
    /// table (declared via
    /// [`super::header::WriterHeader::chromosomes`] in `new`).
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

    /// An interval record's `end` is out of range. The locus interval
    /// is half-open `[start, end)`, so `end` is legal in
    /// `(start, chrom_length + 1]` — it may sit one past the last
    /// 1-based position. This variant carries both bounds so the
    /// message states the true accepted range (unlike
    /// [`Self::PosOutOfRange`], whose `[1, length]` text is wrong for
    /// an exclusive `end`), and distinguishes the `end <= start`
    /// (empty/negative span) failure from `end` past the contig.
    #[error(
        "record at index {record_index}: locus end {end} out of ({start}, {chrom_length}+1] \
         for chrom_id {chrom_id}"
    )]
    LocusEndOutOfRange {
        record_index: u64,
        chrom_id: u32,
        start: u32,
        end: u32,
        chrom_length: u32,
    },

    /// An SSR spanning-read profile carries a `NaN` log-probability.
    /// `-inf` is permitted (a legitimate `log(0)` profile weight), but
    /// `NaN` is never a valid probability and would poison the
    /// downstream EM sum if it round-tripped silently (the SSR decode
    /// path runs no finite sweep — `amb-logliks` is
    /// `finite_constraint: false` because `-inf` is legal).
    #[error(
        "record at index {record_index}: profile {profile_index} carries a NaN log-probability"
    )]
    NonFiniteLoglik {
        record_index: u64,
        profile_index: usize,
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

    /// `finish` has been called but the writer-side block-header
    /// invariant check rejected our own data — bug in the writer.
    /// `kind` carries the specific invariant violated; M8.
    #[error("block emission failed at index {block_index}")]
    BlockEmission {
        block_index: u64,
        #[source]
        kind: BlockHeaderInvariantKind,
    },

    /// A column the writer just encoded has a byte length that
    /// disagrees with the schema-predicted length. M5: this
    /// previously fired only via `debug_assert!`; it is now a real
    /// runtime check so a writer bug surfaces at the producer, not
    /// downstream at read time.
    #[error("column {column:?} self-check failed: encoded {got} bytes, schema predicts {expected}")]
    ColumnSizeSelfCheck {
        column: &'static str,
        got: usize,
        expected: usize,
    },
}

/// Helper for [`PspWriteError::Io`]'s `Display`: renders the
/// optional `(block_index, column_tag)` context inline.
fn format_block_column(block_index: Option<u64>, column_tag: Option<u16>) -> String {
    match (block_index, column_tag) {
        (Some(b), Some(t)) => format!(" (block {b}, column {t:#x})"),
        (Some(b), None) => format!(" (block {b})"),
        (None, Some(t)) => format!(" (column {t:#x})"),
        (None, None) => String::new(),
    }
}
