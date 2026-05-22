// Mi5: this module is `pub(crate)`; `lookup_by_tag` is reader-side
// surface awaiting the not-yet-built `PspReader`.
#![allow(dead_code)]
//! The v1.0 column-tag registry — the single source of truth shared
//! by writer, reader, and (later) the `psp_spec_dump` helper.
//!
//! The byte-layout spec's §"Required columns in v1.0" table and the
//! in-file TOML `[[column]]` array carry the same information; this
//! constant is what they both serialise / deserialise against. Any
//! disagreement between a file's declared schema and this table is a
//! hard error at file-open (spec §"Header-binary consistency: required
//! reader checks", item 6).

/// How a column's entry count is determined inside a block.
///
/// `PerRecord` columns have one entry per per-position record
/// (count = `n_records`); `PerAllele` columns have one entry per
/// allele across all records (count = `n_total_alleles`).
///
/// **Not `#[non_exhaustive]` by design.** Every match on this enum
/// must be exhaustive so adding a new cardinality in v1.x or v2.0
/// is a compile error at every dispatch site (writer encoders,
/// `as_str` / `parse_token`, etc.).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cardinality {
    PerRecord,
    PerAllele,
}

impl Cardinality {
    /// String form matching the TOML enum surface (spec §"Encoding
    /// conventions / Per-field character set" → `column.cardinality`).
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::PerRecord => "per-record",
            Self::PerAllele => "per-allele",
        }
    }

    pub fn parse_token(s: &str) -> Option<Self> {
        match s {
            "per-record" => Some(Self::PerRecord),
            "per-allele" => Some(Self::PerAllele),
            _ => None,
        }
    }
}

/// Per-entry wire shape inside a column's payload.
///
/// - `Scalar`: one `element-type` value per entry.
/// - `List`: per entry, a `varint` count followed by that many
///   `element-type` values.
/// - `Bytes`: a flat byte stream chunked by a paired length column
///   (whose name lives in [`ColumnDef::length_column`]).
///
/// **Not `#[non_exhaustive]` by design** — same reason as
/// [`Cardinality`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Shape {
    Scalar,
    List,
    Bytes,
}

impl Shape {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Scalar => "scalar",
            Self::List => "list",
            Self::Bytes => "bytes",
        }
    }

    pub fn parse_token(s: &str) -> Option<Self> {
        match s {
            "scalar" => Some(Self::Scalar),
            "list" => Some(Self::List),
            "bytes" => Some(Self::Bytes),
            _ => None,
        }
    }
}

/// Per-element type inside a scalar or list column. Absent on
/// `Bytes`-shaped columns.
///
/// **Not `#[non_exhaustive]` by design.** Adding a new variant must
/// be a compile error at every dispatch site
/// ([`Self::fixed_byte_width`], [`Self::as_str`], the writer's
/// `predict_uncompressed_len`, etc.).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElementType {
    U8,
    U16,
    U32,
    U64,
    I32,
    I64,
    F32,
    F64,
    Varint,
    Svarint,
    Bool,
}

impl ElementType {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::U8 => "u8",
            Self::U16 => "u16",
            Self::U32 => "u32",
            Self::U64 => "u64",
            Self::I32 => "i32",
            Self::I64 => "i64",
            Self::F32 => "f32",
            Self::F64 => "f64",
            Self::Varint => "varint",
            Self::Svarint => "svarint",
            Self::Bool => "bool",
        }
    }

    pub fn parse_token(s: &str) -> Option<Self> {
        match s {
            "u8" => Some(Self::U8),
            "u16" => Some(Self::U16),
            "u32" => Some(Self::U32),
            "u64" => Some(Self::U64),
            "i32" => Some(Self::I32),
            "i64" => Some(Self::I64),
            "f32" => Some(Self::F32),
            "f64" => Some(Self::F64),
            "varint" => Some(Self::Varint),
            "svarint" => Some(Self::Svarint),
            "bool" => Some(Self::Bool),
            _ => None,
        }
    }

    /// Byte width of one element on the wire, if it is fixed.
    ///
    /// `None` for `Varint` and `Svarint` (variable width); `Some(n)`
    /// for every fixed-width integer, float, and `Bool` (which is one
    /// byte). Used by the reader to a-priori predict the
    /// `uncompressed_len` of a scalar column from the block-level
    /// counts (spec §"Header-binary consistency", check #7).
    pub const fn fixed_byte_width(self) -> Option<usize> {
        match self {
            Self::U8 | Self::Bool => Some(1),
            Self::U16 => Some(2),
            Self::U32 | Self::I32 | Self::F32 => Some(4),
            Self::U64 | Self::I64 | Self::F64 => Some(8),
            Self::Varint | Self::Svarint => None,
        }
    }
}

/// Hard cap on the bytes of a single allele's sequence (the
/// `allele_seq_len` row of the registry). Spec §"Required columns
/// in v1.0" → 0x03 says: real biological alleles are far shorter,
/// and a value over this cap is a producer bug.
///
/// **Enforcement.** The writer rejects oversized alleles in
/// [`PspWriter::write_record`](crate::per_sample_pileup::psp::writer::PspWriter::write_record).
/// On the read side this cap is honoured by
/// [`decode_bytes_split`](super::block::decode_bytes_split) when
/// the caller passes `Some(MAX_ALLELE_SEQ_LEN)` as its
/// `max_entry_len` argument — the eventual `PspReader` is required
/// to do so for the `allele-seq` column so the registry's
/// single-source-of-truth cap is enforced symmetrically.
pub const MAX_ALLELE_SEQ_LEN: u64 = 10_000;

/// Exhaustive enumeration of every v1.0 column the registry knows
/// about, by static identity. Used by the writer to dispatch
/// per-column encoders without the compile-time hole a `match` on
/// `u16` tag literals leaves (M4). Adding a v1.x column means
/// adding a variant here, which forces every dispatch match to be
/// updated as a compile error.
///
/// **Not `#[non_exhaustive]` by design** — see the comment on
/// [`Cardinality`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ColumnKey {
    DeltaPos,
    NAlleles,
    AlleleSeqLen,
    AlleleSeq,
    AlleleObsCount,
    AlleleQSumLog,
    AlleleFwdCount,
    AllelePlacedLeftCount,
    AllelePlacedStartCount,
    AlleleMapqSum,
    AlleleMapqSumSq,
    AlleleChainIds,
}

/// Per-payload-shape data, pinning the valid combinations of
/// `(shape, element_type, length_column)` that previously lived as
/// three flat `ColumnDef` fields (M17). Only the three combinations
/// reachable in v1.0 are representable.
///
/// **Not `#[non_exhaustive]` by design** — same rationale as the
/// neighbouring enums.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnPayload {
    /// One `element_type` value per entry.
    Scalar { element_type: ElementType },
    /// Per entry: a `varint` count followed by that many
    /// `element_type` values.
    List { element_type: ElementType },
    /// A flat byte stream chunked by the named length column.
    Bytes { length_column: &'static str },
}

impl ColumnPayload {
    /// Recover the legacy `Shape` token for the wire / TOML schema.
    pub const fn shape(self) -> Shape {
        match self {
            Self::Scalar { .. } => Shape::Scalar,
            Self::List { .. } => Shape::List,
            Self::Bytes { .. } => Shape::Bytes,
        }
    }

    /// `Some(element_type)` for `Scalar`/`List`; `None` for `Bytes`.
    pub const fn element_type(self) -> Option<ElementType> {
        match self {
            Self::Scalar { element_type } | Self::List { element_type } => Some(element_type),
            Self::Bytes { .. } => None,
        }
    }

    /// `Some(name)` for `Bytes`; `None` otherwise.
    pub const fn length_column(self) -> Option<&'static str> {
        match self {
            Self::Bytes { length_column } => Some(length_column),
            Self::Scalar { .. } | Self::List { .. } => None,
        }
    }
}

/// One row of the column-tag registry. Identical information appears
/// in the spec's §"Required columns in v1.0" table and in every
/// v1.0 file's TOML `[[column]]` array; this struct is the canonical
/// in-code representation.
#[derive(Debug, Clone, Copy)]
pub struct ColumnDef {
    /// Numeric tag the per-block manifest references. Spec writes
    /// these in hex (`0x01`, `0x11`, ...); the value here is the
    /// same integer.
    pub tag: u16,

    /// Static identity used for compile-time-exhaustive dispatch
    /// (writer encoders, decoders) without going through the numeric
    /// `tag`. Adding a registry entry forces every dispatch match to
    /// be updated.
    pub key: ColumnKey,

    /// Lowercase kebab-case identifier, e.g. `"allele-q-sum-log"`.
    /// Stable across versions.
    pub name: &'static str,

    pub cardinality: Cardinality,

    /// Per-payload-shape data — pins the valid `(shape, element_type,
    /// length_column)` combinations at the type level (M17).
    pub payload: ColumnPayload,

    pub required: bool,

    /// When `true`, the reader sweeps every decoded float entry of
    /// this column and rejects NaN / ±∞ (spec check #9). Only
    /// meaningful for `ColumnPayload::Scalar { element_type: F32 |
    /// F64 }` (and future float-list payloads); silently `false`
    /// otherwise.
    ///
    /// Lifted from a hard-coded reader-side `match` so the set of
    /// finite-constrained columns is *registry-defined* — a future
    /// F-typed column flips this flag instead of remembering to
    /// teach the reader's `decode_block_payload` about it. M11.
    pub finite_constraint: bool,

    /// Canonical one-line description, reproduced verbatim by the
    /// writer into the in-file `[[column]]` array.
    pub description: &'static str,
}

impl ColumnDef {
    /// Convenience accessor for the legacy `shape` field.
    pub const fn shape(&self) -> Shape {
        self.payload.shape()
    }

    /// Convenience accessor for the legacy `element_type` field.
    pub const fn element_type(&self) -> Option<ElementType> {
        self.payload.element_type()
    }

    /// Convenience accessor for the legacy `length_column` field.
    pub const fn length_column(&self) -> Option<&'static str> {
        self.payload.length_column()
    }
}

/// The v1.0 column registry. Order matches the spec's §"Required
/// columns in v1.0" table (tag-ascending); writers emit columns in
/// this order, readers walk the in-file schema and cross-check
/// every column against the matching entry here by `name`.
pub const V1_0_COLUMNS: &[ColumnDef] = &[
    ColumnDef {
        tag: 0x01,
        key: ColumnKey::DeltaPos,
        name: "delta-pos",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::Varint,
        },
        required: true,
        finite_constraint: false,
        description: "Distance to the previous record's position. \
            First record in a block has value 0 and uses the block \
            header's first-pos.",
    },
    ColumnDef {
        tag: 0x02,
        key: ColumnKey::NAlleles,
        name: "n-alleles",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::Varint,
        },
        required: true,
        finite_constraint: false,
        description: "Number of alleles in this record (>= 1). Sum \
            across records equals n_total_alleles.",
    },
    ColumnDef {
        tag: 0x03,
        key: ColumnKey::AlleleSeqLen,
        name: "allele-seq-len",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::Varint,
        },
        required: true,
        finite_constraint: false,
        description: "Byte length of each allele's sequence. Drives \
            chunking of allele-seq.",
    },
    ColumnDef {
        tag: 0x04,
        key: ColumnKey::AlleleSeq,
        name: "allele-seq",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Bytes {
            length_column: "allele-seq-len",
        },
        required: true,
        finite_constraint: false,
        description: "Allele sequence bytes (uppercase ASCII over \
            {A,C,G,T,N}). REF is the first allele in each record.",
    },
    ColumnDef {
        tag: 0x10,
        key: ColumnKey::AlleleObsCount,
        name: "allele-obs-count",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Observation count: reads supporting this allele.",
    },
    ColumnDef {
        tag: 0x11,
        key: ColumnKey::AlleleQSumLog,
        name: "allele-q-sum-log",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::F64,
        },
        required: true,
        finite_constraint: true,
        description: "Sum over supporting reads of max(ln_BQ_BAQ, \
            ln_MQ). The freebayes per-read combination of base and \
            mapping quality.",
    },
    ColumnDef {
        tag: 0x12,
        key: ColumnKey::AlleleFwdCount,
        name: "allele-fwd-count",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Forward-strand count: reads on the forward \
            strand supporting this allele.",
    },
    ColumnDef {
        tag: 0x13,
        key: ColumnKey::AllelePlacedLeftCount,
        name: "allele-placed-left-count",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Reads whose 5' end is to the left of the \
            record's position (freebayes' placedLeft).",
    },
    ColumnDef {
        tag: 0x14,
        key: ColumnKey::AllelePlacedStartCount,
        name: "allele-placed-start-count",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Reads whose 5' end is the record's position \
            (freebayes' placedStart).",
    },
    ColumnDef {
        tag: 0x15,
        key: ColumnKey::AlleleMapqSum,
        name: "allele-mapq-sum",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Σ mapping quality over reads supporting this \
            allele. Combined with allele-obs-count this gives the \
            mean MAPQ; combined with allele-mapq-sum-sq it gives the \
            sample variance for the Welch's t multi-mapper filter.",
    },
    ColumnDef {
        tag: 0x16,
        key: ColumnKey::AlleleMapqSumSq,
        name: "allele-mapq-sum-sq",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U64,
        },
        required: true,
        finite_constraint: false,
        description: "Σ mapq² over reads supporting this allele. With \
            allele-mapq-sum and allele-obs-count, computes sample \
            variance: var = (sum_sq − sum²/n)/(n − 1). u64 to allow \
            extreme-depth amplicons without overflow.",
    },
    ColumnDef {
        tag: 0x22,
        key: ColumnKey::AlleleChainIds,
        name: "allele-chain-ids",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::List {
            element_type: ElementType::U64,
        },
        required: true,
        finite_constraint: false,
        description: "Phase-chain identifiers contributing to each \
            allele observation. Ascending fixed-width little-endian \
            u64s (relying on zstd to compress the leading zero bytes \
            on early ids). Identifiers are unique within the .psp \
            file and never recycled, so two observations sharing an \
            identifier came from the same read or read-pair in this \
            sample.",
    },
];

/// Linear-scan lookup by tag. Twelve entries — no hash table needed.
pub fn lookup_by_tag(tag: u16) -> Option<&'static ColumnDef> {
    V1_0_COLUMNS.iter().find(|c| c.tag == tag)
}

/// Linear-scan lookup by canonical name.
pub fn lookup_by_name(name: &str) -> Option<&'static ColumnDef> {
    V1_0_COLUMNS.iter().find(|c| c.name == name)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// (M17.) The `ColumnPayload` enum *is* the shape invariant —
    /// only `Scalar`/`List` carry an element-type and only `Bytes`
    /// carries a length-column. This test is now mostly redundant
    /// with the type system, but kept as a smoke check.
    #[test]
    fn shape_invariants_per_entry() {
        for c in V1_0_COLUMNS {
            match c.payload {
                ColumnPayload::Bytes { length_column } => {
                    assert!(!length_column.is_empty(), "{}: empty length-column", c.name);
                }
                ColumnPayload::Scalar { .. } | ColumnPayload::List { .. } => {}
            }
        }
    }

    /// Every `length_column` reference resolves to a column that
    /// shares the parent's `cardinality`, has `shape = Scalar`, and
    /// has `element_type = Varint`. (Spec consistency check #3.)
    #[test]
    fn length_column_references_are_well_typed() {
        for c in V1_0_COLUMNS {
            if let Some(len_name) = c.length_column() {
                let target = lookup_by_name(len_name).unwrap_or_else(|| {
                    panic!("{}: length-column '{}' not in registry", c.name, len_name)
                });
                assert_eq!(
                    target.cardinality, c.cardinality,
                    "{}: cardinality mismatch with length-column",
                    c.name
                );
                assert_eq!(
                    target.shape(),
                    Shape::Scalar,
                    "{}: length-column must be Scalar",
                    c.name
                );
                assert_eq!(
                    target.element_type(),
                    Some(ElementType::Varint),
                    "{}: length-column must be Varint",
                    c.name
                );
            }
        }
    }

    /// Tags are unique and strictly ascending in array order (the
    /// spec/writer commit to emitting columns tag-ascending, and
    /// having the registry pre-sorted lets the writer iterate it
    /// directly).
    #[test]
    fn tags_are_unique_and_ascending() {
        let mut prev: Option<u16> = None;
        for c in V1_0_COLUMNS {
            if let Some(p) = prev {
                assert!(
                    c.tag > p,
                    "tags not strictly ascending: {p:#x} then {:#x}",
                    c.tag
                );
            }
            prev = Some(c.tag);
        }
    }

    /// Tags fall in the ranges the spec reserves for v1.0: per-record
    /// structural in `[0x01, 0x0F]`, per-allele scalars in
    /// `[0x10, 0x1F]`, phase-chain in `[0x20, 0x2F]`. `0x00` is
    /// forbidden. (Spec §"Reserved column-tag ranges".)
    #[test]
    fn tags_fall_in_v1_ranges() {
        for c in V1_0_COLUMNS {
            assert_ne!(c.tag, 0x00, "{}: tag 0x00 is reserved", c.name);
            let in_v1_range = (0x01..=0x2F).contains(&c.tag);
            assert!(
                in_v1_range,
                "{}: tag {:#x} outside v1.0 ranges",
                c.name, c.tag
            );
        }
    }

    #[test]
    fn lookup_by_tag_round_trip() {
        for c in V1_0_COLUMNS {
            let looked_up = lookup_by_tag(c.tag).expect("tag should resolve");
            assert_eq!(looked_up.name, c.name);
        }
        assert!(lookup_by_tag(0x00).is_none());
        assert!(lookup_by_tag(0xFFFF).is_none());
    }

    #[test]
    fn lookup_by_name_round_trip() {
        for c in V1_0_COLUMNS {
            let looked_up = lookup_by_name(c.name).expect("name should resolve");
            assert_eq!(looked_up.tag, c.tag);
        }
        assert!(lookup_by_name("nonexistent-column").is_none());
    }

    #[test]
    fn enum_str_round_trip() {
        for c in [Cardinality::PerRecord, Cardinality::PerAllele] {
            assert_eq!(Cardinality::parse_token(c.as_str()), Some(c));
        }
        for s in [Shape::Scalar, Shape::List, Shape::Bytes] {
            assert_eq!(Shape::parse_token(s.as_str()), Some(s));
        }
        for e in [
            ElementType::U8,
            ElementType::U16,
            ElementType::U32,
            ElementType::U64,
            ElementType::I32,
            ElementType::I64,
            ElementType::F32,
            ElementType::F64,
            ElementType::Varint,
            ElementType::Svarint,
            ElementType::Bool,
        ] {
            assert_eq!(ElementType::parse_token(e.as_str()), Some(e));
        }
        assert_eq!(Cardinality::parse_token("bogus"), None);
        assert_eq!(Shape::parse_token("bogus"), None);
        assert_eq!(ElementType::parse_token("bogus"), None);
    }

    #[test]
    fn fixed_byte_widths_match_spec() {
        assert_eq!(ElementType::U8.fixed_byte_width(), Some(1));
        assert_eq!(ElementType::Bool.fixed_byte_width(), Some(1));
        assert_eq!(ElementType::U16.fixed_byte_width(), Some(2));
        assert_eq!(ElementType::U32.fixed_byte_width(), Some(4));
        assert_eq!(ElementType::I32.fixed_byte_width(), Some(4));
        assert_eq!(ElementType::F32.fixed_byte_width(), Some(4));
        assert_eq!(ElementType::U64.fixed_byte_width(), Some(8));
        assert_eq!(ElementType::I64.fixed_byte_width(), Some(8));
        assert_eq!(ElementType::F64.fixed_byte_width(), Some(8));
        assert_eq!(ElementType::Varint.fixed_byte_width(), None);
        assert_eq!(ElementType::Svarint.fixed_byte_width(), None);
    }

    /// Specific spot-check: tag, name, type triples for the load-bearing
    /// per-allele scalars match what the spec's example header carries
    /// and what the implementation produces.
    #[test]
    fn key_columns_have_expected_types() {
        let q = lookup_by_name("allele-q-sum-log").unwrap();
        assert_eq!(q.tag, 0x11);
        assert_eq!(q.element_type(), Some(ElementType::F64));

        let chain_ids = lookup_by_name("allele-chain-ids").unwrap();
        assert_eq!(chain_ids.tag, 0x22);
        assert_eq!(chain_ids.cardinality, Cardinality::PerAllele);
        assert_eq!(chain_ids.shape(), Shape::List);
        assert_eq!(chain_ids.element_type(), Some(ElementType::U64));

        let seq = lookup_by_name("allele-seq").unwrap();
        assert_eq!(seq.shape(), Shape::Bytes);
        assert_eq!(seq.length_column(), Some("allele-seq-len"));
    }

    /// (M4.) `ColumnKey` is exhaustive over every registry entry.
    /// Adding a `ColumnDef` row without giving it a `key` is a
    /// compile error; reusing an existing key for a new row would
    /// fail this uniqueness test.
    #[test]
    fn column_keys_are_unique_and_cover_v1_0() {
        let mut seen = std::collections::HashSet::new();
        for c in V1_0_COLUMNS {
            assert!(
                seen.insert(c.key),
                "duplicate ColumnKey {:?} in V1_0_COLUMNS",
                c.key
            );
        }
        assert_eq!(seen.len(), V1_0_COLUMNS.len());
    }
}
