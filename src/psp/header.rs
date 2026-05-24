//! `.psp` file header — build, parse, validate.
//!
//! The header sits at the very start of every `.psp` file. Layout
//! (spec §"File header / Layout"):
//!
//! ```text
//! +-----------------------+
//! | PSP\n                 |  4 bytes — head magic
//! +-----------------------+
//! | toml_body_length u64  |  8 bytes — little-endian, authoritative
//! +-----------------------+    for the body boundary
//! | <TOML body>           |  exactly toml_body_length bytes
//! +-----------------------+
//! | ---END-HEADER---\n    |  17 bytes — sentinel, cross-check
//! +-----------------------+
//! ```
//!
//! The TOML body carries:
//!
//! - `format-version` (string `"MAJOR.MINOR"`),
//! - `sample`, `reference`, `created`,
//! - `[[chromosome]]` array (≥ 1 entry),
//! - `[writer]` table with provenance and `[writer.parameters]`,
//! - `[[column]]` array — the file's authoritative binary schema.
//!
//! Two layers of types live here. **Wire types** (`Wire*`,
//! `pub(super)` and not exposed) are serde-derived and round-trip
//! through TOML directly. **Public types** (`WriterHeader` on the
//! build side, `ParsedHeader` on the read side) carry strong types
//! (`(u16, u16)` for format version, [`registry::Cardinality`] for
//! column cardinality, etc.) and are what callers see.
//!
//! All per-field rules from spec §"Encoding conventions / Per-field
//! character set" live here as small validators; the validation
//! pass runs them once on the build side and once on the parse
//! side. Spec §"Errors must not pass silently" requires that any
//! rule violation halt — and they do.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use toml::value::Datetime;

use super::errors::{PspReadError, PspWriteError, TomlSerError};
use super::registry::{self, Cardinality, ColumnDef, ElementType, Shape, V1_0_COLUMNS};

// ---------------------------------------------------------------------
// Framing constants
// ---------------------------------------------------------------------

/// The 4-byte head magic — printable ASCII so the trailing newline
/// opens the TOML body on the next line for `head file.psp` output.
pub const HEAD_MAGIC: [u8; 4] = *b"PSP\n";

/// The 17-byte sentinel line that closes the header. Verified by the
/// reader as a structural cross-check against the length prefix.
pub const HEAD_SENTINEL: &[u8; 17] = b"---END-HEADER---\n";

/// Byte length of [`HEAD_SENTINEL`].
pub const HEAD_SENTINEL_LEN: usize = 17;

/// Total fixed framing overhead: magic + length prefix + sentinel.
pub const HEADER_FRAMING_BYTES: usize = 4 + 8 + HEAD_SENTINEL_LEN;

/// Minimum allowed `toml_body_length`. At least one byte of TOML —
/// the empty body cannot be valid because `format-version` is
/// required.
pub const MIN_HEADER_BODY_BYTES: u64 = 1;

/// Maximum allowed `toml_body_length`. The hard 1 MiB cap minus the
/// framing overhead. Validated before any buffer is allocated so a
/// malicious or corrupt file cannot drive a large allocation off
/// this field.
pub const MAX_HEADER_BODY_BYTES: u64 = (1024 * 1024) - HEADER_FRAMING_BYTES as u64;

/// Format version this reader writes and consumes — the v1.0
/// baseline.
pub const READER_FORMAT_VERSION: (u16, u16) = (1, 0);

// ---------------------------------------------------------------------
// Public types — caller-facing
// ---------------------------------------------------------------------

/// Everything the writer needs to assemble a header. The `[[column]]`
/// array is **not** carried here — the writer emits the v1.0 column
/// registry ([`V1_0_COLUMNS`]) verbatim.
#[derive(Debug, Clone)]
pub struct WriterHeader {
    /// `(major, minor)`. The writer accepts exactly
    /// [`READER_FORMAT_VERSION`]; any other value is rejected by
    /// [`build_header_bytes`] with a `PspWriteError::InvalidHeaderField`.
    pub format_version: (u16, u16),
    pub sample: String,
    pub reference: String,
    pub created: Datetime,
    /// Must have at least one entry.
    pub chromosomes: Vec<ChromosomeEntry>,
    pub writer: WriterProvenance,
}

/// One row of the `[[chromosome]]` array. Matches `@SQ` order across
/// the input CRAMs; the zero-based index here is the `chrom_id`
/// referenced by the binary body.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChromosomeEntry {
    pub name: String,
    pub length: u32,
    /// 32 lowercase hex characters — MD5 of the uppercase reference
    /// sequence for this contig (SAM `@SQ M5` semantics).
    pub md5: String,
}

/// Provenance: what produced the file, in a form that supports
/// reproducibility on any host without leaking the producer's
/// directory structure or username.
#[derive(Debug, Clone, PartialEq)]
pub struct WriterProvenance {
    pub tool: String,
    pub version: String,
    pub subcommand: String,
    /// **Basenames only** — no directory component. Spec §"File
    /// header — provenance and path stripping".
    pub input_crams: Vec<String>,
    /// **Basename only**.
    pub input_fasta: String,
    /// One entry per CLI-exposed knob, in any order (`BTreeMap`
    /// gives a deterministic order on the wire).
    pub parameters: BTreeMap<String, ParameterValue>,
}

/// A parameter value: TOML-typed.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum ParameterValue {
    /// TOML integer. `i64` to match TOML's signed integer width.
    Integer(i64),
    /// TOML float. Must be finite (validator rejects NaN / ±∞).
    Float(f64),
    Boolean(bool),
    String(String),
}

/// Output of [`parse_header_bytes`] — the in-file header in strong
/// types after validation.
#[derive(Debug, Clone)]
pub struct ParsedHeader {
    pub format_version: (u16, u16),
    pub sample: String,
    pub reference: String,
    pub created: Datetime,
    pub chromosomes: Vec<ParsedChromosome>,
    pub writer: ParsedWriter,
    /// The file's binary schema, parsed into strong types. The
    /// number and order matches the in-file `[[column]]` array.
    pub columns: Vec<ParsedColumn>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedChromosome {
    pub name: String,
    pub length: u32,
    pub md5: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ParsedWriter {
    pub tool: String,
    pub version: String,
    pub subcommand: String,
    pub input_crams: Vec<String>,
    pub input_fasta: String,
    pub parameters: BTreeMap<String, ParameterValue>,
}

/// A `[[column]]` entry as it appears in the file, after validation
/// against the enum closed-vocabularies. The reader cross-checks
/// every recognised column's fields against [`V1_0_COLUMNS`]
/// separately (spec consistency check #6).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedColumn {
    pub tag: u16,
    pub name: String,
    pub cardinality: Cardinality,
    /// Per-payload-shape data; mirrors [`ColumnDef::payload`] but
    /// uses `String` for the length-column reference because parse
    /// inputs are owned. M17.
    pub payload: ParsedColumnPayload,
    pub required: bool,
    pub description: String,
}

/// Owned counterpart to [`super::registry::ColumnPayload`] —
/// parse-side variant where `length_column` is a `String`. M17.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParsedColumnPayload {
    Scalar { element_type: ElementType },
    List { element_type: ElementType },
    Bytes { length_column: String },
}

impl ParsedColumnPayload {
    pub fn shape(&self) -> Shape {
        match self {
            Self::Scalar { .. } => Shape::Scalar,
            Self::List { .. } => Shape::List,
            Self::Bytes { .. } => Shape::Bytes,
        }
    }

    pub fn element_type(&self) -> Option<ElementType> {
        match self {
            Self::Scalar { element_type } | Self::List { element_type } => Some(*element_type),
            Self::Bytes { .. } => None,
        }
    }

    pub fn length_column(&self) -> Option<&str> {
        match self {
            Self::Bytes { length_column } => Some(length_column.as_str()),
            Self::Scalar { .. } | Self::List { .. } => None,
        }
    }
}

// ---------------------------------------------------------------------
// Wire types — private; what serde rounds-trips through TOML
// ---------------------------------------------------------------------

/// Top-level wire-form of the TOML header. Field order here is what
/// the serializer emits, which is also what the spec example shows:
/// scalars first (format-version, sample, reference, created), then
/// the `[[chromosome]]` array, then `[writer]`, then `[[column]]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct WireHeader {
    #[serde(rename = "format-version")]
    format_version: String,
    sample: String,
    reference: String,
    created: Datetime,
    #[serde(rename = "chromosome", default)]
    chromosomes: Vec<WireChromosome>,
    writer: WireWriter,
    #[serde(rename = "column", default)]
    columns: Vec<WireColumn>,
    /// Unknown top-level keys are captured here so they survive the
    /// round-trip and the reader can apply the "unknown keys are
    /// skipped" rule without losing them silently.
    #[serde(flatten)]
    extras: BTreeMap<String, toml::Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct WireChromosome {
    name: String,
    length: u32,
    md5: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct WireWriter {
    tool: String,
    version: String,
    subcommand: String,
    #[serde(rename = "input-crams")]
    input_crams: Vec<String>,
    #[serde(rename = "input-fasta")]
    input_fasta: String,
    parameters: BTreeMap<String, ParameterValue>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct WireColumn {
    tag: u16,
    name: String,
    cardinality: String,
    shape: String,
    #[serde(
        rename = "element-type",
        default,
        skip_serializing_if = "Option::is_none"
    )]
    element_type: Option<String>,
    #[serde(
        rename = "length-column",
        default,
        skip_serializing_if = "Option::is_none"
    )]
    length_column: Option<String>,
    required: bool,
    description: String,
}

// ---------------------------------------------------------------------
// Build path: WriterHeader → bytes
// ---------------------------------------------------------------------

/// Serialise a [`WriterHeader`] into the TOML body bytes only — no
/// magic, no length prefix, no sentinel. The caller usually wants
/// [`build_header_bytes`]; this lower-level form exists so the
/// writer slice can frame the body without rebuilding it.
pub fn build_header_toml(header: &WriterHeader) -> Result<Vec<u8>, PspWriteError> {
    validate_writer_header(header)?;
    let wire = wire_from_writer_header(header);
    let s = toml::to_string(&wire).map_err(|e| PspWriteError::HeaderToml {
        source: TomlSerError::new(e),
    })?;
    Ok(s.into_bytes())
}

/// Serialise a [`WriterHeader`] into the full, framed header
/// (magic + length prefix + TOML body + sentinel). The first byte
/// of the returned `Vec` is the file's first byte.
pub fn build_header_bytes(header: &WriterHeader) -> Result<Vec<u8>, PspWriteError> {
    let body = build_header_toml(header)?;
    let body_len = u64::try_from(body.len()).expect("header body fits in u64");
    if !(MIN_HEADER_BODY_BYTES..=MAX_HEADER_BODY_BYTES).contains(&body_len) {
        return Err(PspWriteError::HeaderBodyTooLarge {
            got: body_len,
            min: MIN_HEADER_BODY_BYTES,
            max: MAX_HEADER_BODY_BYTES,
        });
    }
    let mut out = Vec::with_capacity(HEADER_FRAMING_BYTES + body.len());
    out.extend_from_slice(&HEAD_MAGIC);
    out.extend_from_slice(&body_len.to_le_bytes());
    out.extend_from_slice(&body);
    out.extend_from_slice(HEAD_SENTINEL);
    Ok(out)
}

fn wire_from_writer_header(header: &WriterHeader) -> WireHeader {
    let format_version = format!("{}.{}", header.format_version.0, header.format_version.1);
    let chromosomes = header
        .chromosomes
        .iter()
        .map(|c| WireChromosome {
            name: c.name.clone(),
            length: c.length,
            md5: c.md5.clone(),
        })
        .collect();
    let writer = WireWriter {
        tool: header.writer.tool.clone(),
        version: header.writer.version.clone(),
        subcommand: header.writer.subcommand.clone(),
        input_crams: header.writer.input_crams.clone(),
        input_fasta: header.writer.input_fasta.clone(),
        parameters: header.writer.parameters.clone(),
    };
    let columns = V1_0_COLUMNS.iter().map(wire_column_from_def).collect();
    WireHeader {
        format_version,
        sample: header.sample.clone(),
        reference: header.reference.clone(),
        created: header.created,
        chromosomes,
        writer,
        columns,
        extras: BTreeMap::new(),
    }
}

fn wire_column_from_def(c: &ColumnDef) -> WireColumn {
    WireColumn {
        tag: c.tag,
        name: c.name.to_string(),
        cardinality: c.cardinality.as_str().to_string(),
        shape: c.shape().as_str().to_string(),
        element_type: c.element_type().map(|e| e.as_str().to_string()),
        length_column: c.length_column().map(|s| s.to_string()),
        required: c.required,
        description: c.description.to_string(),
    }
}

// ---------------------------------------------------------------------
// Parse path: bytes → ParsedHeader
// ---------------------------------------------------------------------

/// Parse a fully-framed header starting at the beginning of `bytes`.
/// Returns the parsed header and the total number of bytes consumed
/// (magic + length prefix + body + sentinel), so the caller can
/// position itself at the first byte of block 0.
pub fn parse_header_bytes(bytes: &[u8]) -> Result<(ParsedHeader, usize), PspReadError> {
    // 1. Magic.
    if bytes.len() < 4 {
        return Err(PspReadError::BadHeadMagic {
            got: [0; 4],
            expected: HEAD_MAGIC,
        });
    }
    let mut got_magic = [0u8; 4];
    got_magic.copy_from_slice(&bytes[0..4]);
    if got_magic != HEAD_MAGIC {
        return Err(PspReadError::BadHeadMagic {
            got: got_magic,
            expected: HEAD_MAGIC,
        });
    }

    // 2. Length prefix.
    if bytes.len() < 12 {
        return Err(PspReadError::Io {
            context: "header length prefix",
            source: std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "less than 12 bytes available",
            ),
        });
    }
    let body_len = u64::from_le_bytes(bytes[4..12].try_into().unwrap());
    if !(MIN_HEADER_BODY_BYTES..=MAX_HEADER_BODY_BYTES).contains(&body_len) {
        return Err(PspReadError::BadHeaderLength {
            got: body_len,
            min: MIN_HEADER_BODY_BYTES,
            max: MAX_HEADER_BODY_BYTES,
        });
    }
    let body_len_usize = body_len as usize;

    // 3. TOML body bytes available?
    let body_start = 12;
    let body_end = body_start + body_len_usize;
    let sentinel_end = body_end + HEAD_SENTINEL_LEN;
    if bytes.len() < sentinel_end {
        return Err(PspReadError::Io {
            context: "header TOML body / sentinel",
            source: std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "header truncated"),
        });
    }

    // 4. Parse TOML.
    let body_bytes = &bytes[body_start..body_end];
    let body_str =
        std::str::from_utf8(body_bytes).map_err(|e| PspReadError::InvalidHeaderField {
            key: "<toml body>".to_string(),
            reason: format!(
                "TOML body is not valid UTF-8 at byte offset {}",
                e.valid_up_to()
            ),
        })?;
    let wire: WireHeader =
        toml::from_str(body_str).map_err(|source| PspReadError::HeaderToml { source })?;

    // 5. Validate + convert wire → parsed.
    let parsed = parsed_from_wire(wire)?;

    // 6. Sentinel cross-check.
    let sentinel_bytes = &bytes[body_end..sentinel_end];
    if sentinel_bytes != HEAD_SENTINEL.as_slice() {
        return Err(PspReadError::SentinelMismatch {
            offset: body_end as u64,
        });
    }

    Ok((parsed, sentinel_end))
}

/// Lower-level variant of [`parse_header_bytes`] that takes the TOML
/// body alone — already extracted between length prefix and
/// sentinel by the caller — and returns the parsed header. Useful
/// for unit tests that hand-craft a body without re-implementing the
/// framing.
pub fn parse_header_toml(body: &[u8]) -> Result<ParsedHeader, PspReadError> {
    let body_str = std::str::from_utf8(body).map_err(|e| PspReadError::InvalidHeaderField {
        key: "<toml body>".to_string(),
        reason: format!(
            "TOML body is not valid UTF-8 at byte offset {}",
            e.valid_up_to()
        ),
    })?;
    let wire: WireHeader =
        toml::from_str(body_str).map_err(|source| PspReadError::HeaderToml { source })?;
    parsed_from_wire(wire)
}

fn parsed_from_wire(wire: WireHeader) -> Result<ParsedHeader, PspReadError> {
    // Per-field rules first — every value validated before structural
    // checks fire, so the error message names the lowest-level rule
    // that failed.
    let format_version = parse_format_version(&wire.format_version)?;
    validate_printable_ascii("sample", &wire.sample)?;
    validate_printable_ascii("reference", &wire.reference)?;

    // M6: reject any file whose major version exceeds the reader's.
    // Higher-major files may rest on layout assumptions a v1.x
    // reader cannot make — spec §"Versioning policy".
    if format_version.0 > READER_FORMAT_VERSION.0 {
        return Err(PspReadError::UnsupportedFormatVersion {
            file_major: format_version.0,
            file_minor: format_version.1,
            reader_major: READER_FORMAT_VERSION.0,
            reader_minor: READER_FORMAT_VERSION.1,
        });
    }

    if wire.chromosomes.is_empty() {
        return Err(PspReadError::InvalidHeaderField {
            key: "chromosome".to_string(),
            reason: "[[chromosome]] array must have at least one entry".to_string(),
        });
    }
    let chromosomes: Vec<ParsedChromosome> = wire
        .chromosomes
        .into_iter()
        .enumerate()
        .map(|(i, c)| parse_chromosome(i, c))
        .collect::<Result<_, _>>()?;

    let writer = parse_writer(wire.writer)?;

    let columns: Vec<ParsedColumn> = wire
        .columns
        .into_iter()
        .enumerate()
        .map(|(i, c)| parse_column(i, c))
        .collect::<Result<_, _>>()?;
    // Tag uniqueness inside the file.
    {
        let mut seen = std::collections::HashSet::new();
        for c in &columns {
            if !seen.insert(c.tag) {
                return Err(PspReadError::InvalidHeaderField {
                    key: format!("column[{}]", c.name),
                    reason: format!("duplicate column tag {:#x}", c.tag),
                });
            }
        }
    }

    // Schema-vs-registry cross-check + required-column coverage.
    cross_check_against_registry(&columns)?;

    Ok(ParsedHeader {
        format_version,
        sample: wire.sample,
        reference: wire.reference,
        created: wire.created,
        chromosomes,
        writer,
        columns,
    })
}

fn parse_chromosome(index: usize, c: WireChromosome) -> Result<ParsedChromosome, PspReadError> {
    let key_prefix = format!("chromosome[{index}]");
    validate_chrom_name(&format!("{key_prefix}.name"), &c.name)?;
    if c.length == 0 || c.length > i32::MAX as u32 {
        return Err(PspReadError::InvalidHeaderField {
            key: format!("{key_prefix}.length"),
            reason: format!("{} is outside [1, 2^31 - 1]; SAM @SQ LN bounds", c.length),
        });
    }
    validate_md5_hex(&format!("{key_prefix}.md5"), &c.md5)?;
    Ok(ParsedChromosome {
        name: c.name,
        length: c.length,
        md5: c.md5,
    })
}

fn parse_writer(w: WireWriter) -> Result<ParsedWriter, PspReadError> {
    validate_printable_ascii_no_path_sep("writer.tool", &w.tool)?;
    validate_printable_ascii_no_path_sep("writer.version", &w.version)?;
    validate_printable_ascii_no_path_sep("writer.subcommand", &w.subcommand)?;
    for (i, c) in w.input_crams.iter().enumerate() {
        validate_basename(&format!("writer.input-crams[{i}]"), c)?;
    }
    validate_basename("writer.input-fasta", &w.input_fasta)?;
    for (k, v) in w.parameters.iter() {
        validate_bare_key(&format!("writer.parameters.{k}"), k)?;
        if let ParameterValue::Float(f) = v
            && !f.is_finite()
        {
            return Err(PspReadError::InvalidHeaderField {
                key: format!("writer.parameters.{k}"),
                reason: format!("non-finite float {f}"),
            });
        }
    }
    Ok(ParsedWriter {
        tool: w.tool,
        version: w.version,
        subcommand: w.subcommand,
        input_crams: w.input_crams,
        input_fasta: w.input_fasta,
        parameters: w.parameters,
    })
}

fn parse_column(index: usize, c: WireColumn) -> Result<ParsedColumn, PspReadError> {
    let key_prefix = format!("column[{index}]");
    validate_column_name(&format!("{key_prefix}.name"), &c.name)?;
    let cardinality = Cardinality::parse_token(&c.cardinality).ok_or_else(|| {
        PspReadError::InvalidHeaderField {
            key: format!("{key_prefix}.cardinality"),
            reason: format!("{:?} is not a valid cardinality token", c.cardinality),
        }
    })?;
    let shape = Shape::parse_token(&c.shape).ok_or_else(|| PspReadError::InvalidHeaderField {
        key: format!("{key_prefix}.shape"),
        reason: format!("{:?} is not a valid shape token", c.shape),
    })?;
    let payload =
        cross_check_shape_consistency(&key_prefix, shape, c.element_type, c.length_column)?;
    Ok(ParsedColumn {
        tag: c.tag,
        name: c.name,
        cardinality,
        payload,
        required: c.required,
        description: c.description,
    })
}

/// Cross-validate the `(shape, element_type, length_column)` triple
/// — only three of the twelve combinations are legal. Extracted
/// from `parse_column` per Mi29 to keep the cross-rule logic in one
/// place. Constructs the typed [`ParsedColumnPayload`] on success.
fn cross_check_shape_consistency(
    key_prefix: &str,
    shape: Shape,
    element_type: Option<String>,
    length_column: Option<String>,
) -> Result<ParsedColumnPayload, PspReadError> {
    let bad = |field: &str, reason: String| PspReadError::InvalidHeaderField {
        key: format!("{key_prefix}.{field}"),
        reason,
    };
    match (shape, element_type.as_deref(), length_column.as_deref()) {
        (Shape::Bytes, None, Some(lc)) => Ok(ParsedColumnPayload::Bytes {
            length_column: lc.to_string(),
        }),
        (Shape::Bytes, Some(_), _) => Err(bad(
            "element-type",
            "bytes-shape columns must not carry element-type".to_string(),
        )),
        (Shape::Bytes, None, None) => Err(bad(
            "length-column",
            "bytes-shape columns require length-column".to_string(),
        )),
        (Shape::Scalar | Shape::List, None, _) => Err(bad(
            "element-type",
            format!("non-bytes shape {:?} requires element-type", shape.as_str()),
        )),
        (Shape::Scalar | Shape::List, Some(_), Some(_)) => Err(bad(
            "length-column",
            format!(
                "non-bytes shape {:?} must not carry length-column",
                shape.as_str()
            ),
        )),
        (Shape::Scalar, Some(et_str), None) => {
            let element_type = ElementType::parse_token(et_str).ok_or_else(|| {
                bad(
                    "element-type",
                    format!("{et_str:?} is not a valid element-type token"),
                )
            })?;
            Ok(ParsedColumnPayload::Scalar { element_type })
        }
        (Shape::List, Some(et_str), None) => {
            let element_type = ElementType::parse_token(et_str).ok_or_else(|| {
                bad(
                    "element-type",
                    format!("{et_str:?} is not a valid element-type token"),
                )
            })?;
            Ok(ParsedColumnPayload::List { element_type })
        }
    }
}

/// Cross-check every column in the file's `[[column]]` array against
/// [`V1_0_COLUMNS`]: required columns the reader knows must match
/// the registry on cardinality / shape / element-type / length-column
/// / required; required columns in the registry must appear in the
/// file; required columns the reader does not recognise abort.
fn cross_check_against_registry(columns: &[ParsedColumn]) -> Result<(), PspReadError> {
    for c in columns {
        if let Some(def) = registry::lookup_by_name(&c.name) {
            check_match(c, def)?;
        } else if c.required {
            return Err(PspReadError::UnknownRequiredColumn {
                name: c.name.clone(),
                tag: c.tag,
            });
        }
        // Unknown optional columns: tolerated and skipped, per the
        // unknown-key / minor-bump-additions-are-optional rules
        // (spec §"Versioning policy").
    }

    // Every required column in the registry must be present in the
    // file.
    for def in V1_0_COLUMNS {
        if !def.required {
            continue;
        }
        if !columns.iter().any(|c| c.name == def.name) {
            return Err(PspReadError::MissingRequiredColumn {
                name: def.name.to_string(),
                tag: def.tag,
            });
        }
    }

    Ok(())
}

fn check_match(file: &ParsedColumn, def: &ColumnDef) -> Result<(), PspReadError> {
    if file.tag != def.tag {
        return Err(PspReadError::SchemaMismatch {
            name: file.name.clone(),
            field: "tag",
            got: format!("{:#x}", file.tag),
            expected: format!("{:#x}", def.tag),
        });
    }
    if file.cardinality != def.cardinality {
        return Err(PspReadError::SchemaMismatch {
            name: file.name.clone(),
            field: "cardinality",
            got: file.cardinality.as_str().to_string(),
            expected: def.cardinality.as_str().to_string(),
        });
    }
    if file.payload.shape() != def.shape() {
        return Err(PspReadError::SchemaMismatch {
            name: file.name.clone(),
            field: "shape",
            got: file.payload.shape().as_str().to_string(),
            expected: def.shape().as_str().to_string(),
        });
    }
    if file.payload.element_type() != def.element_type() {
        return Err(PspReadError::SchemaMismatch {
            name: file.name.clone(),
            field: "element-type",
            got: file
                .payload
                .element_type()
                .map(|e| e.as_str().to_string())
                .unwrap_or_else(|| "<none>".to_string()),
            expected: def
                .element_type()
                .map(|e| e.as_str().to_string())
                .unwrap_or_else(|| "<none>".to_string()),
        });
    }
    if file.payload.length_column() != def.length_column() {
        return Err(PspReadError::SchemaMismatch {
            name: file.name.clone(),
            field: "length-column",
            got: file
                .payload
                .length_column()
                .map(|s| s.to_string())
                .unwrap_or_else(|| "<none>".to_string()),
            expected: def
                .length_column()
                .map(|s| s.to_string())
                .unwrap_or_else(|| "<none>".to_string()),
        });
    }
    if file.required != def.required {
        return Err(PspReadError::SchemaMismatch {
            name: file.name.clone(),
            field: "required",
            got: file.required.to_string(),
            expected: def.required.to_string(),
        });
    }
    Ok(())
}

// ---------------------------------------------------------------------
// Header-field validators (M16: shared generic helper)
// ---------------------------------------------------------------------
//
// Header-field validation rules are identical on the read and write
// sides — only the outer error type differs. M16 collapses the
// previously-duplicated `validate_*` wrappers, the raw `check_*`
// layer, and the ad-hoc inline `validate_writer_header` block into
// a single generic-over-`E` helper `wrap`, where `E` carries the
// `invalid_field(key, reason)` constructor via the
// [`HeaderFieldError`] trait. Adding a new field rule now means
// adding one `check_*` raw function plus one `wrap(...)` call on
// each side — never two parallel branches that can drift.

/// Outer error types that carry an "invalid header field" variant
/// with `{key, reason}` shape. Both [`PspReadError`] and
/// [`PspWriteError`] implement this. M16.
trait HeaderFieldError {
    fn invalid_field(key: String, reason: String) -> Self;
}

impl HeaderFieldError for PspReadError {
    fn invalid_field(key: String, reason: String) -> Self {
        PspReadError::InvalidHeaderField { key, reason }
    }
}

impl HeaderFieldError for PspWriteError {
    fn invalid_field(key: String, reason: String) -> Self {
        PspWriteError::InvalidHeaderField { key, reason }
    }
}

/// Map a `Result<(), String>` from a raw `check_*` function into a
/// typed outer error, tagging it with the TOML-key prefix.
fn wrap<E: HeaderFieldError>(key: &str, r: Result<(), String>) -> Result<(), E> {
    r.map_err(|reason| E::invalid_field(key.to_string(), reason))
}

/// Run every per-field rule on a [`WriterHeader`]. Used by the
/// writer to refuse malformed callers before any byte is emitted.
fn validate_writer_header(header: &WriterHeader) -> Result<(), PspWriteError> {
    validate_header_fields::<PspWriteError>(header)
}

/// Generic per-field validation pass shared between the writer
/// (`validate_writer_header`) and the reader (after the wire types
/// have parsed back into a struct-shaped form). All rules go through
/// [`wrap`] so the error type is uniform.
fn validate_header_fields<E: HeaderFieldError>(header: &WriterHeader) -> Result<(), E> {
    // M6: the writer only emits files at the single
    // `READER_FORMAT_VERSION` it understands.
    if header.format_version != READER_FORMAT_VERSION {
        return Err(E::invalid_field(
            "format-version".to_string(),
            format!(
                "writer only emits {}.{}; got {}.{}",
                READER_FORMAT_VERSION.0,
                READER_FORMAT_VERSION.1,
                header.format_version.0,
                header.format_version.1,
            ),
        ));
    }
    wrap::<E>("sample", check_printable_ascii(&header.sample))?;
    wrap::<E>("reference", check_printable_ascii(&header.reference))?;
    if header.chromosomes.is_empty() {
        return Err(E::invalid_field(
            "chromosome".to_string(),
            "[[chromosome]] array must have at least one entry".to_string(),
        ));
    }
    for (i, c) in header.chromosomes.iter().enumerate() {
        wrap::<E>(&format!("chromosome[{i}].name"), check_chrom_name(&c.name))?;
        if c.length == 0 || c.length > i32::MAX as u32 {
            return Err(E::invalid_field(
                format!("chromosome[{i}].length"),
                format!("{} is outside [1, 2^31 - 1]", c.length),
            ));
        }
        wrap::<E>(&format!("chromosome[{i}].md5"), check_md5_hex(&c.md5))?;
    }
    wrap::<E>(
        "writer.tool",
        check_printable_ascii_no_path_sep(&header.writer.tool),
    )?;
    wrap::<E>(
        "writer.version",
        check_printable_ascii_no_path_sep(&header.writer.version),
    )?;
    wrap::<E>(
        "writer.subcommand",
        check_printable_ascii_no_path_sep(&header.writer.subcommand),
    )?;
    for (i, c) in header.writer.input_crams.iter().enumerate() {
        wrap::<E>(&format!("writer.input-crams[{i}]"), check_basename(c))?;
    }
    wrap::<E>(
        "writer.input-fasta",
        check_basename(&header.writer.input_fasta),
    )?;
    for (k, v) in header.writer.parameters.iter() {
        wrap::<E>(&format!("writer.parameters.{k}"), check_bare_key(k))?;
        if let ParameterValue::Float(f) = v
            && !f.is_finite()
        {
            return Err(E::invalid_field(
                format!("writer.parameters.{k}"),
                format!("non-finite float {f}"),
            ));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------
// Per-field validators (read-side wrappers around `wrap`)
// ---------------------------------------------------------------------

/// Parse and validate a `format-version` string in one pass —
/// returns `(major, minor)` on success or a typed error on any of
/// the rejection conditions (missing `.`, empty halves, non-digit
/// bytes, half exceeds `u16::MAX`). M9 folded the prior
/// separate `validate_format_version_str` step into this function
/// so the call site no longer has a doc-only "caller has validated"
/// invariant backing an `.unwrap()`.
fn parse_format_version(s: &str) -> Result<(u16, u16), PspReadError> {
    let bad = |reason: String| PspReadError::InvalidHeaderField {
        key: "format-version".to_string(),
        reason,
    };
    let Some((a, b)) = s.split_once('.') else {
        return Err(bad(format!("{s:?} does not match [0-9]+\\.[0-9]+")));
    };
    if a.is_empty() || b.is_empty() {
        return Err(bad(format!("{s:?} does not match [0-9]+\\.[0-9]+")));
    }
    for half in [a, b] {
        if !half.bytes().all(|b| b.is_ascii_digit()) {
            return Err(bad(format!("{s:?} contains a non-digit")));
        }
    }
    let major = a
        .parse::<u16>()
        .map_err(|_| bad(format!("major component {a:?} exceeds u16 range")))?;
    let minor = b
        .parse::<u16>()
        .map_err(|_| bad(format!("minor component {b:?} exceeds u16 range")))?;
    Ok((major, minor))
}

// Read-side validators are one-line wrappers over `wrap` + the raw
// `check_*` functions, exactly mirroring the build side. M16
// collapsed the prior per-rule wrappers into these calls.

fn validate_printable_ascii(key: &str, s: &str) -> Result<(), PspReadError> {
    wrap(key, check_printable_ascii(s))
}

fn validate_printable_ascii_no_path_sep(key: &str, s: &str) -> Result<(), PspReadError> {
    wrap(key, check_printable_ascii_no_path_sep(s))
}

fn validate_basename(key: &str, s: &str) -> Result<(), PspReadError> {
    wrap(key, check_basename(s))
}

fn validate_md5_hex(key: &str, s: &str) -> Result<(), PspReadError> {
    wrap(key, check_md5_hex(s))
}

fn validate_chrom_name(key: &str, s: &str) -> Result<(), PspReadError> {
    wrap(key, check_chrom_name(s))
}

fn validate_bare_key(key: &str, s: &str) -> Result<(), PspReadError> {
    wrap(key, check_bare_key(s))
}

fn validate_column_name(key: &str, s: &str) -> Result<(), PspReadError> {
    wrap(key, check_column_name(s))
}

// ---------------------------------------------------------------------
// Raw character-set checks. Each returns `Ok(())` or `Err(reason)`,
// where `reason` is a single sentence explaining which rule failed.
// ---------------------------------------------------------------------

/// ASCII printable: bytes in `[0x20, 0x7E]`. Empty strings are
/// rejected — none of our fields are legitimately empty.
fn check_printable_ascii(s: &str) -> Result<(), String> {
    if s.is_empty() {
        return Err("empty string".to_string());
    }
    for (i, b) in s.bytes().enumerate() {
        if !is_printable_ascii(b) {
            return Err(format!(
                "byte {b:#04x} at position {i} is not ASCII printable"
            ));
        }
    }
    Ok(())
}

/// ASCII printable AND no `/` or `\` (forbidden in tool / version /
/// subcommand fields — those don't carry paths).
fn check_printable_ascii_no_path_sep(s: &str) -> Result<(), String> {
    check_printable_ascii(s)?;
    for (i, b) in s.bytes().enumerate() {
        if b == b'/' || b == b'\\' {
            return Err(format!(
                "path separator {:?} at position {i} not allowed",
                b as char
            ));
        }
    }
    Ok(())
}

/// ASCII printable filename: no `/` or `\` (basename only).
fn check_basename(s: &str) -> Result<(), String> {
    check_printable_ascii_no_path_sep(s)
}

/// 32 lowercase hex characters.
fn check_md5_hex(s: &str) -> Result<(), String> {
    if s.len() != 32 {
        return Err(format!("md5 is {} chars long, expected 32", s.len()));
    }
    for (i, b) in s.bytes().enumerate() {
        if !matches!(b, b'0'..=b'9' | b'a'..=b'f') {
            return Err(format!("non-hex character {:?} at position {i}", b as char));
        }
    }
    Ok(())
}

/// SAM `@SQ SN` regex:
/// `[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*`.
fn check_chrom_name(s: &str) -> Result<(), String> {
    let bytes = s.as_bytes();
    if bytes.is_empty() {
        return Err("empty name".to_string());
    }
    if !is_chrom_first_byte(bytes[0]) {
        return Err(format!(
            "first character {:?} not allowed in @SQ SN",
            bytes[0] as char
        ));
    }
    for (i, &b) in bytes.iter().enumerate().skip(1) {
        if !is_chrom_tail_byte(b) {
            return Err(format!(
                "character {:?} at position {i} not allowed in @SQ SN",
                b as char
            ));
        }
    }
    Ok(())
}

/// TOML bare-key character set: `[A-Za-z0-9_-]+`.
fn check_bare_key(s: &str) -> Result<(), String> {
    if s.is_empty() {
        return Err("empty key".to_string());
    }
    for (i, b) in s.bytes().enumerate() {
        if !(b.is_ascii_alphanumeric() || b == b'_' || b == b'-') {
            return Err(format!(
                "character {:?} at position {i} not allowed in TOML bare key",
                b as char
            ));
        }
    }
    Ok(())
}

/// Column name: lowercase kebab-case `[a-z][a-z0-9-]*`.
fn check_column_name(s: &str) -> Result<(), String> {
    let bytes = s.as_bytes();
    if bytes.is_empty() {
        return Err("empty name".to_string());
    }
    if !bytes[0].is_ascii_lowercase() {
        return Err(format!(
            "first character {:?} must be a lowercase letter",
            bytes[0] as char
        ));
    }
    for (i, &b) in bytes.iter().enumerate().skip(1) {
        if !(b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-') {
            return Err(format!(
                "character {:?} at position {i} not allowed in kebab-case name",
                b as char
            ));
        }
    }
    Ok(())
}

#[inline]
fn is_printable_ascii(b: u8) -> bool {
    (0x20..=0x7E).contains(&b)
}

#[inline]
fn is_chrom_first_byte(b: u8) -> bool {
    b.is_ascii_alphanumeric()
        || matches!(
            b,
            b'!' | b'#'
                | b'$'
                | b'%'
                | b'&'
                | b'+'
                | b'.'
                | b'/'
                | b':'
                | b';'
                | b'?'
                | b'@'
                | b'^'
                | b'_'
                | b'|'
                | b'~'
                | b'-'
        )
}

#[inline]
fn is_chrom_tail_byte(b: u8) -> bool {
    is_chrom_first_byte(b) || b == b'*' || b == b'='
}

#[cfg(test)]
mod tests {
    use super::*;

    // Mi20: the realistic-header fixture moved to the shared
    // `super::super::test_fixtures` module — see that file for the
    // canonical builder.
    use super::super::test_fixtures::realistic_writer_header as minimal_writer_header;

    /// (Plan group H — H1.) Build with one chromosome and a couple
    /// of parameters, frame the bytes, parse them back, assert the
    /// round-trip preserves every field.
    #[test]
    fn round_trip_minimal_header() {
        let header = minimal_writer_header();
        let bytes = build_header_bytes(&header).expect("build should succeed");
        // Magic + length + body + sentinel must be at least 30 bytes.
        assert!(bytes.len() > HEADER_FRAMING_BYTES);
        // First 4 bytes are the magic.
        assert_eq!(&bytes[0..4], &HEAD_MAGIC);
        // Last 17 bytes are the sentinel.
        assert_eq!(&bytes[bytes.len() - HEAD_SENTINEL_LEN..], HEAD_SENTINEL);

        let (parsed, consumed) = parse_header_bytes(&bytes).expect("parse should succeed");
        assert_eq!(consumed, bytes.len());

        assert_eq!(parsed.format_version, (1, 0));
        assert_eq!(parsed.sample, "NA12878");
        assert_eq!(parsed.reference, "GRCh38.fa");
        assert_eq!(parsed.chromosomes.len(), 1);
        assert_eq!(parsed.chromosomes[0].name, "chr1");
        assert_eq!(parsed.chromosomes[0].length, 248956422);
        assert_eq!(parsed.writer.tool, "join_per_sample_vcfs");
        assert_eq!(parsed.writer.parameters.len(), 2);
        // All 12 v1.0 columns appear in the parsed header.
        assert_eq!(parsed.columns.len(), V1_0_COLUMNS.len());
        assert_eq!(parsed.columns[0].name, "delta-pos");
        assert_eq!(parsed.columns.last().unwrap().name, "allele-chain-ids");
    }

    /// (M3.) Pin the literal kebab-case TOML wire keys produced by
    /// the serde renames on `WireHeader`, `WireWriter`, `WireColumn`.
    /// Dropping a `#[serde(rename = ...)]` annotation would cause
    /// every test in this file to keep passing (because the reader
    /// uses the same wire definition) but break file compatibility
    /// with separately-built consumers. This test fails loudly in
    /// that case.
    #[test]
    fn wire_keys_are_pinned() {
        let header = minimal_writer_header();
        let body = String::from_utf8(build_header_toml(&header).unwrap()).unwrap();
        for key in [
            "format-version =",
            "[[chromosome]]",
            "[writer]",
            "input-crams =",
            "input-fasta =",
            "[[column]]",
            "element-type =",
            "length-column =",
        ] {
            assert!(body.contains(key), "wire key {key:?} absent from TOML body");
        }
        // Converse: the snake_case forms must not leak — catches a
        // dropped rename annotation that would otherwise be invisible
        // because both serializer and deserializer would agree on
        // the wrong name.
        for forbidden in [
            "format_version",
            "input_crams",
            "input_fasta",
            "element_type",
            "length_column",
        ] {
            assert!(
                !body.contains(forbidden),
                "snake_case wire key {forbidden:?} leaked into TOML body"
            );
        }
    }

    /// (Plan group H — H1.) Output is roundtrip-stable: parse a
    /// header, write it back, parse again, get identical strong-typed
    /// values. The byte representation may or may not be
    /// byte-identical (TOML serializer cosmetics) but the data must
    /// be.
    #[test]
    fn round_trip_is_data_stable() {
        let header = minimal_writer_header();
        let bytes_1 = build_header_bytes(&header).unwrap();
        let (parsed_1, _) = parse_header_bytes(&bytes_1).unwrap();
        // Reconstruct a WriterHeader from parsed_1 and re-emit.
        let header_2 = WriterHeader {
            format_version: parsed_1.format_version,
            sample: parsed_1.sample.clone(),
            reference: parsed_1.reference.clone(),
            created: parsed_1.created,
            chromosomes: parsed_1
                .chromosomes
                .iter()
                .map(|c| ChromosomeEntry {
                    name: c.name.clone(),
                    length: c.length,
                    md5: c.md5.clone(),
                })
                .collect(),
            writer: WriterProvenance {
                tool: parsed_1.writer.tool.clone(),
                version: parsed_1.writer.version.clone(),
                subcommand: parsed_1.writer.subcommand.clone(),
                input_crams: parsed_1.writer.input_crams.clone(),
                input_fasta: parsed_1.writer.input_fasta.clone(),
                parameters: parsed_1.writer.parameters.clone(),
            },
        };
        let bytes_2 = build_header_bytes(&header_2).unwrap();
        let (parsed_2, _) = parse_header_bytes(&bytes_2).unwrap();
        assert_eq!(parsed_1.format_version, parsed_2.format_version);
        assert_eq!(parsed_1.sample, parsed_2.sample);
        assert_eq!(parsed_1.chromosomes, parsed_2.chromosomes);
        assert_eq!(parsed_1.writer.parameters, parsed_2.writer.parameters);
        assert_eq!(parsed_1.columns.len(), parsed_2.columns.len());
    }

    /// `head file.psp` should start with `PSP` on its own first line.
    /// Test that the head magic places `PSP\n` at byte 0.
    #[test]
    fn head_command_shows_psp_first() {
        let header = minimal_writer_header();
        let bytes = build_header_bytes(&header).unwrap();
        assert_eq!(&bytes[0..4], b"PSP\n");
    }

    // ----- H2: per-field rule rejections, both sides --------------

    /// Helper that builds-then-parses a "broken" wire-form header.
    /// Used when the bug we want to test is structurally fine for the
    /// builder to emit (so we can't trigger it via WriterHeader) but
    /// must still fail at parse time. Returns the parse error.
    fn parse_from_wire(wire: &WireHeader) -> PspReadError {
        let body = toml::to_string(wire).expect("wire serialises");
        parse_header_toml(body.as_bytes()).expect_err("should fail")
    }

    fn valid_wire() -> WireHeader {
        WireHeader {
            format_version: "1.0".to_string(),
            sample: "NA12878".to_string(),
            reference: "GRCh38.fa".to_string(),
            created: "2026-05-13T10:00:00Z".parse().unwrap(),
            chromosomes: vec![WireChromosome {
                name: "chr1".to_string(),
                length: 100,
                md5: "6aef897c3d6ff0c78aff06ac189178dd".to_string(),
            }],
            writer: WireWriter {
                tool: "join_per_sample_vcfs".to_string(),
                version: "0.3.0".to_string(),
                subcommand: "per-sample".to_string(),
                input_crams: vec!["a.cram".to_string()],
                input_fasta: "GRCh38.fa".to_string(),
                parameters: BTreeMap::new(),
            },
            columns: V1_0_COLUMNS.iter().map(wire_column_from_def).collect(),
            extras: BTreeMap::new(),
        }
    }

    #[test]
    fn writer_rejects_non_ascii_sample() {
        let mut header = minimal_writer_header();
        header.sample = "NA12878\u{e9}".to_string();
        let err = build_header_bytes(&header).expect_err("non-ASCII byte must fail");
        match err {
            PspWriteError::InvalidHeaderField { key, .. } => assert_eq!(key, "sample"),
            other => panic!("expected InvalidHeaderField{{sample}}, got {other:?}"),
        }
    }

    // ----- M12: writer-side per-field validator coverage ---------
    //
    // Each test mirrors an `h2_reader_rejects_*` sibling above.
    // Writer-side validation is the only line of defence against a
    // buggy producer emitting a malformed file; the build-side
    // validator path was previously only exercised by the
    // `non_ascii_sample` test above.

    /// `writer.reference` must satisfy printable-ASCII.
    #[test]
    fn writer_rejects_non_ascii_reference() {
        let mut header = minimal_writer_header();
        header.reference = "ref\u{e9}.fa".to_string();
        let err = build_header_bytes(&header).expect_err("must fail");
        match err {
            PspWriteError::InvalidHeaderField { key, .. } => assert_eq!(key, "reference"),
            other => panic!("expected InvalidHeaderField{{reference}}, got {other:?}"),
        }
    }

    #[test]
    fn writer_rejects_empty_chromosomes() {
        let mut header = minimal_writer_header();
        header.chromosomes.clear();
        let err = build_header_bytes(&header).expect_err("must fail");
        match err {
            PspWriteError::InvalidHeaderField { key, .. } => assert_eq!(key, "chromosome"),
            other => panic!("expected InvalidHeaderField{{chromosome}}, got {other:?}"),
        }
    }

    #[test]
    fn writer_rejects_invalid_chrom_name() {
        let mut header = minimal_writer_header();
        header.chromosomes[0].name = "ch<".to_string();
        let err = build_header_bytes(&header).expect_err("must fail");
        match err {
            PspWriteError::InvalidHeaderField { key, .. } => {
                assert_eq!(key, "chromosome[0].name");
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    #[test]
    fn writer_rejects_zero_chrom_length() {
        let mut header = minimal_writer_header();
        header.chromosomes[0].length = 0;
        let err = build_header_bytes(&header).expect_err("must fail");
        match err {
            PspWriteError::InvalidHeaderField { key, .. } => {
                assert_eq!(key, "chromosome[0].length");
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    #[test]
    fn writer_rejects_short_md5() {
        let mut header = minimal_writer_header();
        header.chromosomes[0].md5 = "deadbeef".to_string();
        let err = build_header_bytes(&header).expect_err("must fail");
        match err {
            PspWriteError::InvalidHeaderField { key, .. } => {
                assert_eq!(key, "chromosome[0].md5");
            }
            other => panic!("expected InvalidHeaderField{{md5}}, got {other:?}"),
        }
    }

    #[test]
    fn writer_rejects_path_separator_in_input_crams() {
        let mut header = minimal_writer_header();
        header.writer.input_crams[0] = "subdir/a.cram".to_string();
        let err = build_header_bytes(&header).expect_err("must fail");
        match err {
            PspWriteError::InvalidHeaderField { key, .. } => {
                assert_eq!(key, "writer.input-crams[0]");
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    #[test]
    fn writer_rejects_path_separator_in_input_fasta() {
        let mut header = minimal_writer_header();
        header.writer.input_fasta = "subdir/ref.fa".to_string();
        let err = build_header_bytes(&header).expect_err("must fail");
        match err {
            PspWriteError::InvalidHeaderField { key, .. } => {
                assert_eq!(key, "writer.input-fasta");
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    #[test]
    fn writer_rejects_path_separator_in_tool() {
        let mut header = minimal_writer_header();
        header.writer.tool = "bin/tool".to_string();
        let err = build_header_bytes(&header).expect_err("must fail");
        match err {
            PspWriteError::InvalidHeaderField { key, .. } => {
                assert_eq!(key, "writer.tool");
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    #[test]
    fn writer_rejects_invalid_parameter_bare_key() {
        let mut header = minimal_writer_header();
        // Space is not in the TOML bare-key set [A-Za-z0-9_-].
        header
            .writer
            .parameters
            .insert("bad key".to_string(), ParameterValue::Integer(0));
        let err = build_header_bytes(&header).expect_err("must fail");
        match err {
            PspWriteError::InvalidHeaderField { key, .. } => {
                assert!(key.starts_with("writer.parameters."));
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    #[test]
    fn writer_rejects_non_finite_parameter_float() {
        let mut header = minimal_writer_header();
        header
            .writer
            .parameters
            .insert("frac".to_string(), ParameterValue::Float(f64::NAN));
        let err = build_header_bytes(&header).expect_err("must fail");
        match err {
            PspWriteError::InvalidHeaderField { key, reason } => {
                assert_eq!(key, "writer.parameters.frac");
                assert!(reason.contains("non-finite"));
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_non_ascii_sample() {
        let mut wire = valid_wire();
        wire.sample = "NA12878\u{e9}".to_string();
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::InvalidHeaderField { key, .. } => assert_eq!(key, "sample"),
            other => panic!("expected InvalidHeaderField{{sample}}, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_invalid_chrom_name() {
        let mut wire = valid_wire();
        // `<` is not in @SQ SN's allowed set.
        wire.chromosomes[0].name = "ch<".to_string();
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::InvalidHeaderField { key, reason } => {
                assert_eq!(key, "chromosome[0].name");
                assert!(reason.contains("not allowed"));
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_short_md5() {
        let mut wire = valid_wire();
        wire.chromosomes[0].md5 = "deadbeef".to_string();
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::InvalidHeaderField { key, .. } => {
                assert_eq!(key, "chromosome[0].md5");
            }
            other => panic!("expected InvalidHeaderField{{md5}}, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_zero_chrom_length() {
        let mut wire = valid_wire();
        wire.chromosomes[0].length = 0;
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::InvalidHeaderField { key, .. } => {
                assert_eq!(key, "chromosome[0].length");
            }
            other => panic!("expected InvalidHeaderField{{length}}, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_path_separator_in_input_crams() {
        let mut wire = valid_wire();
        wire.writer.input_crams[0] = "subdir/a.cram".to_string();
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::InvalidHeaderField { key, .. } => {
                assert_eq!(key, "writer.input-crams[0]");
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_invalid_column_element_type() {
        let mut wire = valid_wire();
        // Find allele-obs-count and lie about its element-type.
        let c = wire
            .columns
            .iter_mut()
            .find(|c| c.name == "allele-obs-count")
            .unwrap();
        c.element_type = Some("u128".to_string());
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::InvalidHeaderField { key, .. } => {
                assert!(key.contains("element-type"));
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_bytes_shape_without_length_column() {
        let mut wire = valid_wire();
        let c = wire
            .columns
            .iter_mut()
            .find(|c| c.name == "allele-seq")
            .unwrap();
        c.length_column = None;
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::InvalidHeaderField { key, .. } => {
                assert!(key.contains("length-column"));
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_scalar_shape_with_length_column() {
        let mut wire = valid_wire();
        let c = wire
            .columns
            .iter_mut()
            .find(|c| c.name == "allele-obs-count")
            .unwrap();
        c.length_column = Some("allele-seq-len".to_string());
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::InvalidHeaderField { key, .. } => {
                assert!(key.contains("length-column"));
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    // ----- H3: schema-vs-registry disagreement --------------------

    #[test]
    fn reader_rejects_q_sum_log_as_f32() {
        let mut wire = valid_wire();
        let c = wire
            .columns
            .iter_mut()
            .find(|c| c.name == "allele-q-sum-log")
            .unwrap();
        c.element_type = Some("f32".to_string());
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::SchemaMismatch {
                name,
                field,
                got,
                expected,
            } => {
                assert_eq!(name, "allele-q-sum-log");
                assert_eq!(field, "element-type");
                assert_eq!(got, "f32");
                assert_eq!(expected, "f64");
            }
            other => panic!("expected SchemaMismatch, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_n_alleles_cardinality_flipped_to_per_allele() {
        let mut wire = valid_wire();
        let c = wire
            .columns
            .iter_mut()
            .find(|c| c.name == "n-alleles")
            .unwrap();
        // n-alleles is per-record; force it to per-allele to
        // trigger a SchemaMismatch.
        c.cardinality = "per-allele".to_string();
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::SchemaMismatch { name, field, .. } => {
                assert_eq!(name, "n-alleles");
                assert_eq!(field, "cardinality");
            }
            other => panic!("expected SchemaMismatch, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_missing_required_column() {
        let mut wire = valid_wire();
        // Drop allele-seq.
        wire.columns.retain(|c| c.name != "allele-seq");
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::MissingRequiredColumn { name, .. } => {
                assert_eq!(name, "allele-seq");
            }
            other => panic!("expected MissingRequiredColumn, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_unknown_required_column() {
        let mut wire = valid_wire();
        wire.columns.push(WireColumn {
            tag: 0x77,
            name: "future-required".to_string(),
            cardinality: "per-record".to_string(),
            shape: "scalar".to_string(),
            element_type: Some("u32".to_string()),
            length_column: None,
            required: true,
            description: "imagined future column".to_string(),
        });
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::UnknownRequiredColumn { name, tag } => {
                assert_eq!(name, "future-required");
                assert_eq!(tag, 0x77);
            }
            other => panic!("expected UnknownRequiredColumn, got {other:?}"),
        }
    }

    #[test]
    fn reader_tolerates_unknown_optional_column() {
        let mut wire = valid_wire();
        wire.columns.push(WireColumn {
            tag: 0x33,
            name: "future-optional".to_string(),
            cardinality: "per-record".to_string(),
            shape: "scalar".to_string(),
            element_type: Some("u32".to_string()),
            length_column: None,
            required: false,
            description: "imagined future optional column".to_string(),
        });
        let body = toml::to_string(&wire).unwrap();
        let parsed = parse_header_toml(body.as_bytes()).expect("optional unknown should pass");
        // The optional unknown column is parsed and carried through.
        assert!(parsed.columns.iter().any(|c| c.name == "future-optional"));
    }

    // ----- H4: unknown top-level keys ------------------------------

    #[test]
    fn reader_tolerates_unknown_top_level_keys() {
        let mut wire = valid_wire();
        wire.extras.insert(
            "future-feature".to_string(),
            toml::Value::String("hi".into()),
        );
        let body = toml::to_string(&wire).unwrap();
        let parsed = parse_header_toml(body.as_bytes()).expect("unknown top-level should pass");
        assert_eq!(parsed.sample, "NA12878");
    }

    // ----- H5: length prefix and sentinel framing -----------------

    #[test]
    fn reader_rejects_zero_body_length() {
        // Construct: magic + length=0 + sentinel — no body at all.
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&HEAD_MAGIC);
        bytes.extend_from_slice(&0u64.to_le_bytes());
        bytes.extend_from_slice(HEAD_SENTINEL);
        let err = parse_header_bytes(&bytes).expect_err("zero-length body must fail");
        match err {
            PspReadError::BadHeaderLength { got, .. } => assert_eq!(got, 0),
            other => panic!("expected BadHeaderLength, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_huge_body_length() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&HEAD_MAGIC);
        bytes.extend_from_slice(&u64::MAX.to_le_bytes());
        // Don't actually allocate the body — the length check must
        // fire before any allocation.
        let err = parse_header_bytes(&bytes).expect_err("huge length must fail");
        match err {
            PspReadError::BadHeaderLength { got, .. } => assert_eq!(got, u64::MAX),
            other => panic!("expected BadHeaderLength, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_mismatched_sentinel() {
        let header = minimal_writer_header();
        let mut bytes = build_header_bytes(&header).unwrap();
        // Corrupt one byte of the trailing sentinel.
        let sentinel_start = bytes.len() - HEAD_SENTINEL_LEN;
        bytes[sentinel_start] = b'X';
        let err = parse_header_bytes(&bytes).expect_err("bad sentinel must fail");
        match err {
            PspReadError::SentinelMismatch { offset } => {
                assert_eq!(offset, sentinel_start as u64);
            }
            other => panic!("expected SentinelMismatch, got {other:?}"),
        }
    }

    #[test]
    fn reader_rejects_bad_magic() {
        let header = minimal_writer_header();
        let mut bytes = build_header_bytes(&header).unwrap();
        bytes[0] = b'X';
        let err = parse_header_bytes(&bytes).expect_err("bad magic must fail");
        match err {
            PspReadError::BadHeadMagic { .. } => {}
            other => panic!("expected BadHeadMagic, got {other:?}"),
        }
    }

    #[test]
    fn writer_emits_length_prefix_matching_body() {
        let header = minimal_writer_header();
        let bytes = build_header_bytes(&header).unwrap();
        let claimed_len = u64::from_le_bytes(bytes[4..12].try_into().unwrap());
        let actual_body_len = bytes.len() - 4 - 8 - HEAD_SENTINEL_LEN;
        assert_eq!(claimed_len as usize, actual_body_len);
    }

    // ----- M6: format-version refusal on both sides --------------

    /// Writer rejects any `format_version` other than the single
    /// version it actually produces (`READER_FORMAT_VERSION`).
    #[test]
    fn writer_rejects_non_v1_0_format_version() {
        let mut header = minimal_writer_header();
        header.format_version = (2, 0);
        let err =
            build_header_bytes(&header).expect_err("writer must reject non-v1.0 format version");
        match err {
            PspWriteError::InvalidHeaderField { key, .. } => {
                assert_eq!(key, "format-version");
            }
            other => panic!("expected InvalidHeaderField, got {other:?}"),
        }
    }

    /// Reader rejects any file declaring a higher major version
    /// than it supports — guards against a future v2.0 writer's
    /// output silently flowing through a v1.0 reader.
    #[test]
    fn reader_rejects_higher_major_format_version() {
        let mut wire = valid_wire();
        wire.format_version = "2.0".to_string();
        let err = parse_from_wire(&wire);
        match err {
            PspReadError::UnsupportedFormatVersion {
                file_major,
                reader_major,
                ..
            } => {
                assert_eq!(file_major, 2);
                assert_eq!(reader_major, 1);
            }
            other => panic!("expected UnsupportedFormatVersion, got {other:?}"),
        }
    }

    /// Reader still accepts higher *minor* versions of the same
    /// major (forward-compat per spec §"Versioning policy").
    #[test]
    fn reader_accepts_higher_minor_format_version() {
        let mut wire = valid_wire();
        wire.format_version = "1.99".to_string();
        let body = toml::to_string(&wire).unwrap();
        let parsed = parse_header_toml(body.as_bytes())
            .expect("higher minor must parse — minor bumps are forward-compatible");
        assert_eq!(parsed.format_version, (1, 99));
    }

    // ----- Per-field check unit tests (table-driven) -------------

    #[test]
    fn printable_ascii_checks() {
        assert!(check_printable_ascii("NA12878").is_ok());
        assert!(check_printable_ascii("name with spaces").is_ok());
        assert!(check_printable_ascii("").is_err());
        assert!(check_printable_ascii("\u{e9}").is_err());
        assert!(check_printable_ascii("name\twith\ttab").is_err()); // 0x09 is below 0x20
    }

    #[test]
    fn md5_hex_checks() {
        assert!(check_md5_hex(&"a".repeat(32)).is_ok());
        assert!(check_md5_hex(&"0".repeat(32)).is_ok());
        assert!(check_md5_hex(&"f".repeat(32)).is_ok());
        assert!(check_md5_hex(&"A".repeat(32)).is_err()); // uppercase not allowed
        assert!(check_md5_hex(&"a".repeat(31)).is_err()); // wrong length
        assert!(check_md5_hex(&"g".repeat(32)).is_err()); // non-hex
    }

    #[test]
    fn chrom_name_checks() {
        assert!(check_chrom_name("chr1").is_ok());
        assert!(check_chrom_name("GL000192.1").is_ok());
        assert!(check_chrom_name("chrUn_KI270448v1").is_ok());
        assert!(check_chrom_name("").is_err());
        assert!(check_chrom_name("*chr1").is_err()); // can't start with `*`
        assert!(check_chrom_name("=chr1").is_err()); // can't start with `=`
        assert!(check_chrom_name("ch<").is_err()); // `<` is forbidden
        assert!(check_chrom_name("name with space").is_err());
    }

    #[test]
    fn column_name_checks() {
        assert!(check_column_name("delta-pos").is_ok());
        assert!(check_column_name("a1").is_ok());
        assert!(check_column_name("").is_err());
        assert!(check_column_name("1abc").is_err()); // can't start with digit
        assert!(check_column_name("Delta-pos").is_err()); // uppercase not allowed
        assert!(check_column_name("delta_pos").is_err()); // underscore not allowed
    }

    #[test]
    fn format_version_checks() {
        assert_eq!(parse_format_version("1.0").unwrap(), (1, 0));
        assert_eq!(parse_format_version("1.10").unwrap(), (1, 10));
        assert_eq!(parse_format_version("100.200").unwrap(), (100, 200));
        assert!(parse_format_version("1").is_err());
        // "1.0.0" splits on the first '.' into ("1", "0.0"); the
        // second half contains a non-digit '.', so the digit check
        // rejects it.
        assert!(parse_format_version("1.0.0").is_err());
        assert!(parse_format_version("v1.0").is_err());
        assert!(parse_format_version(".").is_err());
        assert!(parse_format_version("").is_err());
        // Each half must fit in u16 (M9 regression).
        assert!(parse_format_version("65536.0").is_err());
    }
}
