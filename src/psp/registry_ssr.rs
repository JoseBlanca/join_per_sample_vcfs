// New SSR `.psp` schema (architecture §10.4). Net-new code with no
// production consumers yet (the Stage-1 SSR driver is deferred); the
// round-trip test exercises it. `#![allow(dead_code)]` mirrors the rest
// of the pre-consumer SSR module.
#![allow(dead_code)]
//! The `ssr` `.psp` schema — the second [`PspKind`] (after `snp`).
//!
//! Maps Stage-1's per-locus evidence ([`crate::ssr::pileup`]'s all-CSR
//! `SsrLocusRecord`) onto the generic container. The fit is the SNP
//! 2-level shape (architecture §10.2): a **locus is a record**, a
//! **spanning-read profile is an entry** (`n-spanning` is the per-record
//! grouping count, exactly as `n-alleles` groups alleles in SNP), and
//! each profile carries two parallel per-profile CSR ragged lists —
//! `amb-lengths` (u16 repeat counts) and `amb-logliks` (f32 renormalized
//! log-probabilities). No `hist_*` (all-CSR, §11 revision); off-ladder
//! columns are deferred until off-ladder candidate generation lands.
//!
//! Every column lands on a codec [`super::block`] already has — varint /
//! fixed-width scalar / CSR ragged-list — so this is "a table + a record
//! mapping," no new wire code (§10.4).

use super::block::{
    BlockHeader, decode_list_column_csr, decode_scalar_column_pod, decode_varint_column,
    encode_list_column_csr, encode_scalar_column, encode_varint_column,
};
use super::errors::{PspReadError, PspWriteError};
use super::kind::{BlockAccumulator, BlockDecoder, PspKind};
use super::reader::{read_and_inflate_column, read_compressed_blob};
use super::registry::{Cardinality, ColumnDef, ColumnPayload, ElementType};

/// The SSR schema-family `kind` tag (header §10.3). Selected by
/// [`super::registry::columns_for_kind`].
pub const SSR_KIND: &str = "ssr";

/// `#[allow]`-free initial hints for the per-locus / per-profile / CSR
/// buffers in [`SsrBlock`]. SSR blocks hold far fewer records than SNP
/// (loci are sparse), so modest fixed hints suffice — they grow if a
/// dense catalog needs more.
const INITIAL_LOCI_HINT: usize = 256;
const INITIAL_PROFILES_HINT: usize = 4096;
const INITIAL_PAIRS_HINT: usize = 8192;

/// Static identity for the SSR columns — compile-time-exhaustive
/// dispatch in the encode/decode, mirroring SNP's `ColumnKey`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SsrColumnKey {
    DeltaStart,
    Span,
    NSpanning,
    Depth,
    NFlanking,
    NFrr,
    NFiltered,
    MappedReads,
    AmbLengths,
    AmbLogliks,
}

impl SsrColumnKey {
    /// The SSR tag for this key. Exhaustive `match Self`, so adding a
    /// variant forces an arm — the source of the key↔tag pairing (the
    /// schema-agnostic [`ColumnDef`] carries no key).
    const fn tag(self) -> u16 {
        match self {
            Self::DeltaStart => 0x01,
            Self::Span => 0x02,
            Self::NSpanning => 0x03,
            Self::Depth => 0x04,
            Self::NFlanking => 0x05,
            Self::NFrr => 0x06,
            Self::NFiltered => 0x07,
            Self::MappedReads => 0x08,
            Self::AmbLengths => 0x10,
            Self::AmbLogliks => 0x11,
        }
    }

    fn from_tag(tag: u16) -> Option<Self> {
        [
            Self::DeltaStart,
            Self::Span,
            Self::NSpanning,
            Self::Depth,
            Self::NFlanking,
            Self::NFrr,
            Self::NFiltered,
            Self::MappedReads,
            Self::AmbLengths,
            Self::AmbLogliks,
        ]
        .into_iter()
        .find(|k| k.tag() == tag)
    }
}

/// The SSR column registry. Per-locus structural/QC columns in
/// `[0x01, 0x0F]`; per-profile CSR lists in `[0x10, 0x1F]` (mirrors the
/// SNP per-record / per-allele tag-range split).
pub const SSR_COLUMNS: &[ColumnDef] = &[
    ColumnDef {
        tag: 0x01,
        name: "delta-start",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::Varint,
        },
        required: true,
        finite_constraint: false,
        description: "Distance to the previous locus's start. First locus \
            in a block has value 0 and uses the block header's first-pos.",
    },
    ColumnDef {
        tag: 0x02,
        name: "span",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::Varint,
        },
        required: true,
        finite_constraint: false,
        description: "Locus length: end - start. The locus interval is \
            [start, start + span).",
    },
    ColumnDef {
        tag: 0x03,
        name: "n-spanning",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::Varint,
        },
        required: true,
        finite_constraint: false,
        description: "Number of spanning-read profiles at this locus. \
            Groups the per-profile entries (amb-* CSR) into loci, as \
            n-alleles groups alleles in the snp schema. Sum across loci \
            equals the per-profile total.",
    },
    ColumnDef {
        tag: 0x04,
        name: "depth",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Total reads considered at the locus (pre-triage).",
    },
    ColumnDef {
        tag: 0x05,
        name: "n-flanking",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Reads anchored on one side of the repeat but not \
            spanning it.",
    },
    ColumnDef {
        tag: 0x06,
        name: "n-frr",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Fully-repeat reads (in-repeat, no unique anchor).",
    },
    ColumnDef {
        tag: 0x07,
        name: "n-filtered",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Reads dropped by the admission gate (low-MAPQ / \
            duplicate / clipped).",
    },
    ColumnDef {
        tag: 0x08,
        name: "mapped-reads",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Mapped reads overlapping the locus window — the \
            normalized-depth denominator.",
    },
    ColumnDef {
        tag: 0x10,
        name: "amb-lengths",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::List {
            element_type: ElementType::U16,
        },
        required: true,
        finite_constraint: false,
        description: "Per spanning-read profile: the on-ladder repeat unit \
            counts the read is consistent with, ascending. Parallel to \
            amb-logliks (same per-profile lengths).",
    },
    ColumnDef {
        tag: 0x11,
        name: "amb-logliks",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::List {
            element_type: ElementType::F32,
        },
        required: true,
        finite_constraint: false,
        description: "Per spanning-read profile: the renormalized log- \
            probabilities (linear-space sum 1) aligned with amb-lengths.",
    },
];

/// One sample's evidence at one locus, **container form** — the
/// chrom_id-keyed mirror of the chrom-name-keyed in-memory
/// [`crate::ssr::pileup`] record (the Stage-1 driver adapts between
/// them, resolving the chromosome name against the header table, exactly
/// as Stage-1 SNP records become `chrom_id`-keyed `PileupRecord`s). The
/// per-spanning-read profiles are the all-CSR storage; `n_spanning` is
/// `spanning.len()` (the grouping count), not stored separately.
#[derive(Debug, Clone, PartialEq)]
pub struct SsrLocusRecord {
    pub chrom_id: u32,
    /// 1-based locus start.
    pub start: u32,
    /// Locus end (exclusive); the interval is `[start, end)`.
    pub end: u32,
    pub depth: u32,
    pub n_flanking: u32,
    pub n_frr: u32,
    pub n_filtered: u32,
    pub mapped_reads: u32,
    /// One renormalized `Qᵣ` profile per spanning read: `(on-ladder
    /// length, log-probability)` pairs, ascending by length.
    pub spanning: Vec<Vec<(u16, f32)>>,
}

// ---------------------------------------------------------------------
// SsrKind
// ---------------------------------------------------------------------

/// The SSR `.psp` schema.
pub struct SsrKind;

impl PspKind for SsrKind {
    type Record = SsrLocusRecord;
    type Block = SsrBlock;
    type Decoder = SsrDecoder;
    const KIND: &'static str = SSR_KIND;

    fn columns() -> &'static [ColumnDef] {
        SSR_COLUMNS
    }

    fn encode_column(
        def: &ColumnDef,
        block: &SsrBlock,
        out: &mut Vec<u8>,
    ) -> Result<(), PspWriteError> {
        // `ColumnDef` is schema-agnostic; recover the SSR dispatch key
        // from the tag (always `Some` — `def` is an `SSR_COLUMNS` row).
        let key = SsrColumnKey::from_tag(def.tag)
            .expect("encode_column called with a non-SSR column tag");
        match key {
            SsrColumnKey::DeltaStart => encode_varint_column(&block.delta_start, out),
            SsrColumnKey::Span => encode_varint_column(&block.span, out),
            SsrColumnKey::NSpanning => encode_varint_column(&block.n_spanning, out),
            SsrColumnKey::Depth => encode_scalar_column(&block.depth, out),
            SsrColumnKey::NFlanking => encode_scalar_column(&block.n_flanking, out),
            SsrColumnKey::NFrr => encode_scalar_column(&block.n_frr, out),
            SsrColumnKey::NFiltered => encode_scalar_column(&block.n_filtered, out),
            SsrColumnKey::MappedReads => encode_scalar_column(&block.mapped_reads, out),
            // amb-lengths / amb-logliks share the per-profile offsets.
            SsrColumnKey::AmbLengths => {
                encode_list_column_csr(&block.amb_lengths_data, &block.amb_offsets, out)
            }
            SsrColumnKey::AmbLogliks => {
                encode_list_column_csr(&block.amb_logliks_data, &block.amb_offsets, out)
            }
        }
        Ok(())
    }

    fn record_interval(record: &SsrLocusRecord) -> (u32, u32, u32) {
        // SSR loci are real intervals; `end` is already exclusive.
        (record.chrom_id, record.start, record.end)
    }
}

// ---------------------------------------------------------------------
// SsrBlock — the write-side accumulator
// ---------------------------------------------------------------------

/// SSR per-block accumulator. Per-locus column buffers + the two
/// per-profile CSR slabs (sharing one offsets vec) + structural
/// metadata. Mirror of [`super::writer::SnpBlock`].
pub struct SsrBlock {
    chrom_id: u32,
    first_pos: u32,
    last_pos: u32,
    prev_start: u32,
    // Per-locus columns.
    delta_start: Vec<u64>,
    span: Vec<u64>,
    n_spanning: Vec<u64>,
    depth: Vec<u32>,
    n_flanking: Vec<u32>,
    n_frr: Vec<u32>,
    n_filtered: Vec<u32>,
    mapped_reads: Vec<u32>,
    // Per-profile CSR. `amb_offsets` (len = n_profiles + 1) indexes both
    // data slabs identically.
    amb_offsets: Vec<u32>,
    amb_lengths_data: Vec<u16>,
    amb_logliks_data: Vec<f32>,
    projected_bytes: usize,
}

impl BlockAccumulator for SsrBlock {
    type Record = SsrLocusRecord;

    fn new_block(chrom_id: u32, first_pos: u32) -> Self {
        let mut amb_offsets = Vec::with_capacity(INITIAL_PROFILES_HINT + 1);
        amb_offsets.push(0);
        Self {
            chrom_id,
            first_pos,
            last_pos: first_pos,
            prev_start: first_pos,
            delta_start: Vec::with_capacity(INITIAL_LOCI_HINT),
            span: Vec::with_capacity(INITIAL_LOCI_HINT),
            n_spanning: Vec::with_capacity(INITIAL_LOCI_HINT),
            depth: Vec::with_capacity(INITIAL_LOCI_HINT),
            n_flanking: Vec::with_capacity(INITIAL_LOCI_HINT),
            n_frr: Vec::with_capacity(INITIAL_LOCI_HINT),
            n_filtered: Vec::with_capacity(INITIAL_LOCI_HINT),
            mapped_reads: Vec::with_capacity(INITIAL_LOCI_HINT),
            amb_offsets,
            amb_lengths_data: Vec::with_capacity(INITIAL_PAIRS_HINT),
            amb_logliks_data: Vec::with_capacity(INITIAL_PAIRS_HINT),
            projected_bytes: 0,
        }
    }

    fn append(&mut self, record: &SsrLocusRecord) {
        let is_first = self.delta_start.is_empty();
        let delta = if is_first {
            0u64
        } else {
            (record.start - self.prev_start) as u64
        };
        self.delta_start.push(delta);
        self.span.push((record.end - record.start) as u64);
        self.n_spanning.push(record.spanning.len() as u64);
        self.depth.push(record.depth);
        self.n_flanking.push(record.n_flanking);
        self.n_frr.push(record.n_frr);
        self.n_filtered.push(record.n_filtered);
        self.mapped_reads.push(record.mapped_reads);

        for profile in &record.spanning {
            for &(length, loglik) in profile {
                self.amb_lengths_data.push(length);
                self.amb_logliks_data.push(loglik);
            }
            self.amb_offsets.push(self.amb_lengths_data.len() as u32);
        }

        self.prev_start = record.start;
        self.last_pos = self.last_pos.max(record.end);

        // Rough projection: per-locus scalars + per-pair (u16 + f32) +
        // per-profile varint count. Doesn't need to be precise.
        let pairs: usize = record.spanning.iter().map(|p| p.len()).sum();
        self.projected_bytes += 2 + 1 + 4 * 5 // deltas/span/n-spanning varints + five u32s
            + record.spanning.len() // per-profile CSR count varints
            + pairs * (2 + 4); // u16 length + f32 loglik per pair
    }

    fn chrom_id(&self) -> u32 {
        self.chrom_id
    }
    fn first_pos(&self) -> u32 {
        self.first_pos
    }
    fn last_pos(&self) -> u32 {
        self.last_pos
    }
    fn n_records(&self) -> u32 {
        self.delta_start.len() as u32
    }
    fn n_entries(&self) -> u32 {
        // Per-profile total = number of CSR rows = offsets.len() - 1.
        (self.amb_offsets.len() - 1) as u32
    }
    fn projected_bytes(&self) -> usize {
        self.projected_bytes
    }
}

// ---------------------------------------------------------------------
// SsrDecoder — the read-side decoder
// ---------------------------------------------------------------------

/// SSR per-block decoder. Owns the decoded per-locus columns + the two
/// per-profile CSR slabs (reused across blocks) + the per-block record
/// cursor. Mirror of [`super::reader::SnpDecoder`].
pub struct SsrDecoder {
    chrom_id: u32,
    first_pos: u32,
    n_records: u32,
    // Per-locus columns.
    delta_start: Vec<u64>,
    span: Vec<u64>,
    n_spanning: Vec<u64>,
    depth: Vec<u32>,
    n_flanking: Vec<u32>,
    n_frr: Vec<u32>,
    n_filtered: Vec<u32>,
    mapped_reads: Vec<u32>,
    // Per-profile CSR (reused across blocks).
    lengths_data: Vec<u16>,
    lengths_offsets: Vec<u32>,
    logliks_data: Vec<f32>,
    logliks_offsets: Vec<u32>,
    // Cursor.
    next_record: u32,
    next_profile: u32,
    last_start: u32,
    loaded: bool,
}

impl BlockDecoder for SsrDecoder {
    type Record = SsrLocusRecord;

    fn new_decoder() -> Self {
        Self {
            chrom_id: 0,
            first_pos: 0,
            n_records: 0,
            delta_start: Vec::new(),
            span: Vec::new(),
            n_spanning: Vec::new(),
            depth: Vec::new(),
            n_flanking: Vec::new(),
            n_frr: Vec::new(),
            n_filtered: Vec::new(),
            mapped_reads: Vec::new(),
            lengths_data: Vec::new(),
            lengths_offsets: Vec::new(),
            logliks_data: Vec::new(),
            logliks_offsets: Vec::new(),
            next_record: 0,
            next_profile: 0,
            last_start: 0,
            loaded: false,
        }
    }

    fn decode_block<R: std::io::Read>(
        &mut self,
        source: &mut R,
        header: &BlockHeader,
        budget: u64,
        decompressor: &mut zstd::bulk::Decompressor<'static>,
        compressed_scratch: &mut Vec<u8>,
        decompressed_scratch: &mut Vec<u8>,
    ) -> Result<(), PspReadError> {
        // B1: every required SSR column must appear in this block's manifest.
        for def in SSR_COLUMNS {
            if def.required && !header.manifest.iter().any(|e| e.tag == def.tag) {
                return Err(PspReadError::MissingRequiredColumnInManifest {
                    name: def.name.to_string(),
                    tag: def.tag,
                });
            }
        }

        let n_records = header.n_records as usize;
        let n_profiles = header.n_total_alleles as usize;
        let mut remaining = budget;

        for entry in &header.manifest {
            let def = SSR_COLUMNS.iter().find(|d| d.tag == entry.tag);
            let Some(def) = def else {
                // Unknown (future/optional) column — read past it so the
                // source cursor stays aligned, then skip.
                let _ = read_compressed_blob(source, entry, remaining)?;
                remaining = remaining.saturating_sub(entry.compressed_len as u64);
                continue;
            };
            read_and_inflate_column(
                source,
                entry,
                remaining,
                def.name,
                decompressor,
                compressed_scratch,
                decompressed_scratch,
            )?;
            remaining = remaining.saturating_sub(entry.compressed_len as u64);
            let bytes: &[u8] = decompressed_scratch;
            let key = SsrColumnKey::from_tag(def.tag)
                .expect("decode_block resolved a non-SSR column tag");
            match key {
                SsrColumnKey::DeltaStart => {
                    self.delta_start = decode_varint_column(bytes, n_records, def.name)?;
                }
                SsrColumnKey::Span => {
                    self.span = decode_varint_column(bytes, n_records, def.name)?;
                }
                SsrColumnKey::NSpanning => {
                    self.n_spanning = decode_varint_column(bytes, n_records, def.name)?;
                }
                SsrColumnKey::Depth => {
                    self.depth = decode_scalar_column_pod::<u32>(bytes, n_records, def.name)?;
                }
                SsrColumnKey::NFlanking => {
                    self.n_flanking = decode_scalar_column_pod::<u32>(bytes, n_records, def.name)?;
                }
                SsrColumnKey::NFrr => {
                    self.n_frr = decode_scalar_column_pod::<u32>(bytes, n_records, def.name)?;
                }
                SsrColumnKey::NFiltered => {
                    self.n_filtered = decode_scalar_column_pod::<u32>(bytes, n_records, def.name)?;
                }
                SsrColumnKey::MappedReads => {
                    self.mapped_reads =
                        decode_scalar_column_pod::<u32>(bytes, n_records, def.name)?;
                }
                SsrColumnKey::AmbLengths => {
                    decode_list_column_csr::<u16>(
                        bytes,
                        n_profiles,
                        def.name,
                        &mut self.lengths_data,
                        &mut self.lengths_offsets,
                    )?;
                }
                SsrColumnKey::AmbLogliks => {
                    decode_list_column_csr::<f32>(
                        bytes,
                        n_profiles,
                        def.name,
                        &mut self.logliks_data,
                        &mut self.logliks_offsets,
                    )?;
                }
            }
        }

        // The two per-profile CSR columns must agree on row boundaries
        // (they were written from the same offsets); a disagreement is a
        // malformed/foreign file.
        if self.lengths_offsets != self.logliks_offsets {
            return Err(PspReadError::Io {
                context: "amb-lengths / amb-logliks CSR offsets disagree",
                source: std::io::Error::other("ssr profile column mismatch"),
            });
        }

        self.chrom_id = header.chrom_id;
        self.first_pos = header.first_pos;
        self.n_records = header.n_records;
        self.next_record = 0;
        self.next_profile = 0;
        self.last_start = header.first_pos;
        self.loaded = true;
        Ok(())
    }

    fn next_record(&mut self) -> Option<Result<SsrLocusRecord, PspReadError>> {
        if !self.loaded || self.next_record >= self.n_records {
            return None;
        }
        let i = self.next_record as usize;
        // delta-start[0] == 0 by writer invariant ⇒ start = first_pos.
        let delta = self.delta_start[i];
        let delta_u32 = match u32::try_from(delta) {
            Ok(d) => d,
            Err(_) => {
                return Some(Err(PspReadError::ColumnElementDecode {
                    column: "delta-start".to_string(),
                    entry: i,
                    source: super::errors::ScalarDecodeError::VarintOverflow,
                }));
            }
        };
        let start = if i == 0 {
            self.first_pos
        } else {
            self.last_start.saturating_add(delta_u32)
        };
        let end = start.saturating_add(self.span[i] as u32);

        let n_here = self.n_spanning[i] as usize;
        let profile_start = self.next_profile as usize;
        let profile_end = profile_start + n_here;
        // Bounds proof for the per-profile loop (offsets carry
        // n_profiles + 1 entries; B1 + the column-count checks make this
        // hold on well-formed input).
        if profile_end >= self.lengths_offsets.len() {
            return Some(Err(PspReadError::Io {
                context: "ssr profile index past CSR offsets",
                source: std::io::Error::other("n-spanning exceeds decoded profiles"),
            }));
        }

        let mut spanning = Vec::with_capacity(n_here);
        for j in profile_start..profile_end {
            let lo = self.lengths_offsets[j] as usize;
            let hi = self.lengths_offsets[j + 1] as usize;
            let profile: Vec<(u16, f32)> = self.lengths_data[lo..hi]
                .iter()
                .copied()
                .zip(self.logliks_data[lo..hi].iter().copied())
                .collect();
            spanning.push(profile);
        }

        self.last_start = start;
        self.next_record += 1;
        self.next_profile += n_here as u32;
        Some(Ok(SsrLocusRecord {
            chrom_id: self.chrom_id,
            start,
            end,
            depth: self.depth[i],
            n_flanking: self.n_flanking[i],
            n_frr: self.n_frr[i],
            n_filtered: self.n_filtered[i],
            mapped_reads: self.mapped_reads[i],
            spanning,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::psp::PspReader;
    use crate::psp::test_fixtures::writer_header;
    use crate::psp::writer::PspWriter;
    use std::io::Cursor;

    fn rec(chrom_id: u32, start: u32, end: u32, spanning: Vec<Vec<(u16, f32)>>) -> SsrLocusRecord {
        SsrLocusRecord {
            chrom_id,
            start,
            end,
            depth: 30,
            n_flanking: 3,
            n_frr: 1,
            n_filtered: 2,
            mapped_reads: 33,
            spanning,
        }
    }

    /// SSR columns are well-formed: unique ascending tags, `key_ssr`
    /// covers every entry, ranges sane.
    #[test]
    fn ssr_columns_are_well_formed() {
        let mut prev = 0u16;
        let mut seen = std::collections::HashSet::new();
        for c in SSR_COLUMNS {
            assert!(
                c.tag > prev,
                "tags strictly ascending: {prev:#x} {:#x}",
                c.tag
            );
            prev = c.tag;
            let key = SsrColumnKey::from_tag(c.tag)
                .unwrap_or_else(|| panic!("{} (tag {:#x}) has no key", c.name, c.tag));
            assert_eq!(key.tag(), c.tag, "{}: from_tag/tag round-trip", c.name);
            assert!(seen.insert(key), "duplicate key {key:?}");
            assert!(
                (0x01..=0x2F).contains(&c.tag),
                "{}: tag out of range",
                c.name
            );
        }
        assert_eq!(seen.len(), SSR_COLUMNS.len());
    }

    /// Full round-trip: write a multi-block `.ssr.psp` via the generic
    /// writer, read it back through the typed reader, assert per-record
    /// equality. Exercises SsrBlock encode + SsrDecoder decode end-to-end
    /// over the shared container (the first §10.6-step-5 hook).
    #[test]
    fn ssr_round_trips_through_the_container() {
        // Window 100 + generous byte cap so the genomic grid drives the
        // block cuts deterministically; multiple loci per/over windows
        // plus a chromosome boundary → multiple blocks.
        let records = vec![
            rec(0, 50, 56, vec![vec![(12, 0.0)]]),
            rec(
                0,
                60,
                70,
                vec![vec![(11, -0.7f32), (12, -0.69f32)], vec![(12, 0.0)]],
            ),
            rec(0, 150, 162, vec![]), // a locus with no spanning reads
            rec(0, 270, 276, vec![vec![(8, -0.2f32), (9, -1.6f32)]]),
            rec(
                1,
                5,
                11,
                vec![
                    vec![(20, 0.0)],
                    vec![(20, 0.0)],
                    vec![(19, -3.0f32), (20, -0.05f32)],
                ],
            ),
        ];

        let sink = Cursor::new(Vec::<u8>::new());
        let mut writer = PspWriter::<_, SsrKind>::new_ssr_with_block_layout(
            sink,
            writer_header(2),
            16 * 1024 * 1024,
            100,
        )
        .unwrap();
        for r in &records {
            writer.write_locus(r).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();

        // Header declares the ssr kind.
        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.header().kind, "ssr");
        // Multiple blocks were cut (windows 0, 1, 2 on chr0 + chr1).
        assert!(reader.block_index().len() >= 4);

        let read_back: Vec<SsrLocusRecord> = reader
            .records_of::<SsrKind>()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(read_back, records);
    }
}
