// SSR `.ssr.psp` schema — Mark 2 (the empirical-candidate model). The second
// [`PspKind`] (after `snp`). Net-new wire shape with one consumer (the Stage-1
// `ssr` driver). `#![allow(dead_code)]` mirrors the rest of the pre-consumer SSR
// surface.
#![allow(dead_code)]
//! The `ssr` `.psp` schema — **observed repeat-region sequences + counts**.
//!
//! Maps Stage-1's per-locus evidence ([`crate::ssr::pileup`]'s `SsrLocusObs`) onto
//! the generic container. The fit is the SNP 2-level shape (architecture §10.2): a
//! **locus is a record**, an **observed sequence is an entry** (`n-obs` is the
//! per-record grouping count, exactly as `n-alleles` groups alleles in SNP), and
//! each observation carries a `count` (varint), a `seq-len` (varint), and its
//! sequence `bytes` — the last two paired exactly as the SNP `allele-seq` /
//! `allele-seq-len` columns, so this rides existing wire codecs (varint / fixed
//! scalar / bytes-by-length-column), no new block code.
//!
//! Mark-2 has no on/off-ladder and no per-read likelihood profiles: a repeat
//! allele is just a byte sequence, the signal is the observation count, and
//! quality was gated in Stage 1 (model `ssr_ladder_model.md`).
//!
//! Glossary:
//! - **obs** — a distinct observed repeat-region sequence at a locus, with the
//!   number of reads that showed it (`obs-count`).

use super::block::{
    BlockHeader, decode_bytes_split, decode_scalar_column_pod, decode_varint_column,
    encode_bytes_concat, encode_scalar_column, encode_varint_column,
};
use super::errors::{PspReadError, PspWriteError};
use super::kind::{BlockAccumulator, BlockDecoder, PspKind};
use super::reader::{read_and_inflate_column, read_compressed_blob};
use super::registry::{Cardinality, ColumnDef, ColumnPayload, ElementType, MAX_ALLELE_SEQ_LEN};

/// The SSR schema-family `kind` tag (header §10.3).
pub const SSR_KIND: &str = "ssr";

/// Capacity hints for the per-locus / per-observation / byte buffers. SSR blocks
/// hold far fewer records than SNP (loci are sparse); modest fixed hints suffice.
const INITIAL_LOCI_HINT: usize = 256;
const INITIAL_OBS_HINT: usize = 1024; // ~256 loci × ~4 distinct sequences
const INITIAL_OBS_BYTES_HINT: usize = 32 * 1024; // ~1024 obs × ~32 bp tracts

super::kind::column_key! {
    /// Static identity for the SSR columns — compile-time-exhaustive dispatch in
    /// the encode/decode, mirroring SNP's `ColumnKey`. Per-locus columns in
    /// `[0x01, 0x0F]`; per-observation columns in `[0x10, 0x1F]`.
    pub enum SsrColumnKey {
        DeltaStart = 0x01,
        Span = 0x02,
        NObs = 0x03,
        Depth = 0x04,
        NFiltered = 0x05,
        MappedReads = 0x06,
        NLowQuality = 0x07,
        NBorderOffEnd = 0x08,
        ObsCount = 0x10,
        ObsSeqLen = 0x11,
        ObsSeq = 0x12,
    }
}

/// The SSR column registry.
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
        name: "n-obs",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::Varint,
        },
        required: true,
        finite_constraint: false,
        description: "Number of distinct observed repeat-region sequences at \
            this locus. Groups the per-observation entries into loci, as \
            n-alleles groups alleles in the snp schema.",
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
        description: "Usable primary reads considered at the locus.",
    },
    ColumnDef {
        tag: 0x05,
        name: "n-filtered",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Reads dropped by the admission gate (low-MAPQ / \
            duplicate / qc-fail / short).",
    },
    ColumnDef {
        tag: 0x06,
        name: "mapped-reads",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Dup-free primary coverage denominator.",
    },
    ColumnDef {
        tag: 0x07,
        name: "n-low-quality",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Reads dropped by the first-quartile repeat-region \
            quality gate.",
    },
    ColumnDef {
        tag: 0x08,
        name: "n-border-off-end",
        cardinality: Cardinality::PerRecord,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::U32,
        },
        required: true,
        finite_constraint: false,
        description: "Reads where a flank ran off the read end (allele \
            >= read length).",
    },
    ColumnDef {
        tag: 0x10,
        name: "obs-count",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::Varint,
        },
        required: true,
        finite_constraint: false,
        description: "Per observed sequence: the number of reads that showed \
            it. Parallel to obs-seq (same per-observation order).",
    },
    ColumnDef {
        tag: 0x11,
        name: "obs-seq-len",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Scalar {
            element_type: ElementType::Varint,
        },
        required: true,
        finite_constraint: false,
        description: "Byte length of each observed sequence. Drives chunking \
            of obs-seq.",
    },
    ColumnDef {
        tag: 0x12,
        name: "obs-seq",
        cardinality: Cardinality::PerAllele,
        payload: ColumnPayload::Bytes {
            length_column: "obs-seq-len",
        },
        required: true,
        finite_constraint: false,
        description: "Observed repeat-region sequence bytes, ascending by \
            bytes within each locus.",
    },
];

/// One sample's evidence at one locus, **container form** — the chrom_id-keyed,
/// 1-based mirror of the chrom-name-keyed in-memory [`crate::ssr::pileup`] record.
/// `n_obs` is `observed.len()` (the grouping count), not stored separately.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SsrLocusRecord {
    pub chrom_id: u32,
    /// 1-based locus start.
    pub start: u32,
    /// Locus end (exclusive); the interval is `[start, end)`.
    pub end: u32,
    pub depth: u32,
    pub n_filtered: u32,
    pub mapped_reads: u32,
    pub n_low_quality: u32,
    pub n_border_off_end: u32,
    /// Distinct observed repeat-region sequences → observation count, ascending
    /// by bytes.
    pub observed: Vec<(Box<[u8]>, u32)>,
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
        // UNREACHABLE: `def` is an `SSR_COLUMNS` row, so its tag is an `SsrColumnKey`.
        let key = SsrColumnKey::from_tag(def.tag)
            .expect("encode_column called with a non-SSR column tag");
        match key {
            SsrColumnKey::DeltaStart => encode_varint_column(&block.delta_start, out),
            SsrColumnKey::Span => encode_varint_column(&block.span, out),
            SsrColumnKey::NObs => encode_varint_column(&block.n_obs, out),
            SsrColumnKey::Depth => encode_scalar_column(&block.depth, out),
            SsrColumnKey::NFiltered => encode_scalar_column(&block.n_filtered, out),
            SsrColumnKey::MappedReads => encode_scalar_column(&block.mapped_reads, out),
            SsrColumnKey::NLowQuality => encode_scalar_column(&block.n_low_quality, out),
            SsrColumnKey::NBorderOffEnd => encode_scalar_column(&block.n_border_off_end, out),
            SsrColumnKey::ObsCount => encode_varint_column(&block.obs_count, out),
            SsrColumnKey::ObsSeqLen => encode_varint_column(&block.obs_seq_len, out),
            // The per-entry chunking lives in obs-seq-len; the payload is the
            // concatenated bytes as-is.
            SsrColumnKey::ObsSeq => encode_bytes_concat(&[&block.obs_seq_bytes[..]], out),
        }
        Ok(())
    }

    fn record_interval(record: &SsrLocusRecord) -> (u32, u32, u32) {
        (record.chrom_id, record.start, record.end)
    }
}

// ---------------------------------------------------------------------
// SsrBlock — the write-side accumulator
// ---------------------------------------------------------------------

/// SSR per-block accumulator: per-locus column buffers + the per-observation
/// count/seq-len vecs + the concatenated sequence bytes. Mirror of
/// [`super::writer::SnpBlock`].
pub struct SsrBlock {
    chrom_id: u32,
    first_pos: u32,
    last_pos: u32,
    prev_start: u32,
    // Per-locus columns.
    delta_start: Vec<u64>,
    span: Vec<u64>,
    n_obs: Vec<u64>,
    depth: Vec<u32>,
    n_filtered: Vec<u32>,
    mapped_reads: Vec<u32>,
    n_low_quality: Vec<u32>,
    n_border_off_end: Vec<u32>,
    // Per-observation columns.
    obs_count: Vec<u64>,
    obs_seq_len: Vec<u64>,
    obs_seq_bytes: Vec<u8>,
    projected_bytes: usize,
}

impl BlockAccumulator for SsrBlock {
    type Record = SsrLocusRecord;

    fn new_block(chrom_id: u32, first_pos: u32) -> Self {
        Self {
            chrom_id,
            first_pos,
            last_pos: first_pos,
            prev_start: first_pos,
            delta_start: Vec::with_capacity(INITIAL_LOCI_HINT),
            span: Vec::with_capacity(INITIAL_LOCI_HINT),
            n_obs: Vec::with_capacity(INITIAL_LOCI_HINT),
            depth: Vec::with_capacity(INITIAL_LOCI_HINT),
            n_filtered: Vec::with_capacity(INITIAL_LOCI_HINT),
            mapped_reads: Vec::with_capacity(INITIAL_LOCI_HINT),
            n_low_quality: Vec::with_capacity(INITIAL_LOCI_HINT),
            n_border_off_end: Vec::with_capacity(INITIAL_LOCI_HINT),
            obs_count: Vec::with_capacity(INITIAL_OBS_HINT),
            obs_seq_len: Vec::with_capacity(INITIAL_OBS_HINT),
            obs_seq_bytes: Vec::with_capacity(INITIAL_OBS_BYTES_HINT),
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
        self.n_obs.push(record.observed.len() as u64);
        self.depth.push(record.depth);
        self.n_filtered.push(record.n_filtered);
        self.mapped_reads.push(record.mapped_reads);
        self.n_low_quality.push(record.n_low_quality);
        self.n_border_off_end.push(record.n_border_off_end);

        for (seq, count) in &record.observed {
            self.obs_count.push(*count as u64);
            self.obs_seq_len.push(seq.len() as u64);
            self.obs_seq_bytes.extend_from_slice(seq);
        }

        self.prev_start = record.start;
        // M1: store the inclusive last-touched position (the block index treats
        // `last_pos` as an inclusive max). `record.end` is exclusive; `end > start
        // >= 1` is guaranteed by `validate_locus`, so `end - 1` never underflows.
        self.last_pos = self.last_pos.max(record.end - 1);

        let total_seq_bytes: usize = record.observed.iter().map(|(s, _)| s.len()).sum();
        self.projected_bytes += 2 + 1 + 4 * 5 // deltas/span/n-obs varints + five u32s
            + record.observed.len() * (1 + 1) // per-obs count + seq-len varints
            + total_seq_bytes; // per-obs sequence bytes
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
        self.obs_count.len() as u32
    }
    fn projected_bytes(&self) -> usize {
        self.projected_bytes
    }
}

// ---------------------------------------------------------------------
// SsrDecoder — the read-side decoder
// ---------------------------------------------------------------------

/// SSR per-block decoder. Owns the decoded per-locus columns + the per-observation
/// vecs + the per-block record cursor. Mirror of [`super::reader::SnpDecoder`].
pub struct SsrDecoder {
    chrom_id: u32,
    first_pos: u32,
    n_records: u32,
    // Per-locus columns.
    delta_start: Vec<u64>,
    span: Vec<u64>,
    n_obs: Vec<u64>,
    depth: Vec<u32>,
    n_filtered: Vec<u32>,
    mapped_reads: Vec<u32>,
    n_low_quality: Vec<u32>,
    n_border_off_end: Vec<u32>,
    // Per-observation columns.
    obs_count: Vec<u64>,
    obs_seq_len: Vec<u64>,
    obs_seq: Vec<Vec<u8>>,
    // Cursor.
    next_record: u32,
    next_obs: u32,
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
            n_obs: Vec::new(),
            depth: Vec::new(),
            n_filtered: Vec::new(),
            mapped_reads: Vec::new(),
            n_low_quality: Vec::new(),
            n_border_off_end: Vec::new(),
            obs_count: Vec::new(),
            obs_seq_len: Vec::new(),
            obs_seq: Vec::new(),
            next_record: 0,
            next_obs: 0,
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
        let n_obs_total = header.n_total_alleles as usize;
        let mut remaining = budget;

        for entry in &header.manifest {
            let def = SSR_COLUMNS.iter().find(|d| d.tag == entry.tag);
            let Some(def) = def else {
                // Unknown (future/optional) column — read past it so the source
                // cursor stays aligned, then skip.
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
            // UNREACHABLE: `def` resolved via `SSR_COLUMNS.iter().find(...)`.
            let key = SsrColumnKey::from_tag(def.tag)
                .expect("decode_block resolved a non-SSR column tag");
            match key {
                SsrColumnKey::DeltaStart => {
                    self.delta_start = decode_varint_column(bytes, n_records, def.name)?;
                }
                SsrColumnKey::Span => {
                    self.span = decode_varint_column(bytes, n_records, def.name)?;
                }
                SsrColumnKey::NObs => {
                    self.n_obs = decode_varint_column(bytes, n_records, def.name)?;
                }
                SsrColumnKey::Depth => {
                    self.depth = decode_scalar_column_pod::<u32>(bytes, n_records, def.name)?;
                }
                SsrColumnKey::NFiltered => {
                    self.n_filtered = decode_scalar_column_pod::<u32>(bytes, n_records, def.name)?;
                }
                SsrColumnKey::MappedReads => {
                    self.mapped_reads =
                        decode_scalar_column_pod::<u32>(bytes, n_records, def.name)?;
                }
                SsrColumnKey::NLowQuality => {
                    self.n_low_quality =
                        decode_scalar_column_pod::<u32>(bytes, n_records, def.name)?;
                }
                SsrColumnKey::NBorderOffEnd => {
                    self.n_border_off_end =
                        decode_scalar_column_pod::<u32>(bytes, n_records, def.name)?;
                }
                SsrColumnKey::ObsCount => {
                    self.obs_count = decode_varint_column(bytes, n_obs_total, def.name)?;
                }
                SsrColumnKey::ObsSeqLen => {
                    self.obs_seq_len = decode_varint_column(bytes, n_obs_total, def.name)?;
                }
                SsrColumnKey::ObsSeq => {
                    // obs-seq-len (tag 0x11) precedes obs-seq (0x12) in the
                    // tag-ordered manifest, so it is already decoded here.
                    self.obs_seq = decode_bytes_split(
                        bytes,
                        &self.obs_seq_len,
                        def.name,
                        Some(MAX_ALLELE_SEQ_LEN),
                    )?;
                }
            }
        }

        // M2: the per-record grouping column (`n-obs`) must sum to exactly the
        // decoded per-observation entry count (`n_total_alleles`).
        let declared: u64 = self.n_obs.iter().sum();
        if declared != header.n_total_alleles as u64 {
            return Err(PspReadError::SsrProfileCountMismatch {
                expected: header.n_total_alleles,
                got: declared,
            });
        }

        self.chrom_id = header.chrom_id;
        self.first_pos = header.first_pos;
        self.n_records = header.n_records;
        self.next_record = 0;
        self.next_obs = 0;
        self.last_start = header.first_pos;
        self.loaded = true;
        Ok(())
    }

    fn next_record(&mut self) -> Option<Result<SsrLocusRecord, PspReadError>> {
        if !self.loaded || self.next_record >= self.n_records {
            return None;
        }
        let i = self.next_record as usize;
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
        let span = match u32::try_from(self.span[i]) {
            Ok(s) => s,
            Err(_) => {
                return Some(Err(PspReadError::ColumnElementDecode {
                    column: "span".to_string(),
                    entry: i,
                    source: super::errors::ScalarDecodeError::VarintOverflow,
                }));
            }
        };
        let end = start.saturating_add(span);

        let n_here = self.n_obs[i] as usize;
        let obs_start = self.next_obs as usize;
        let obs_end = obs_start + n_here;
        // Bounds proof: M2's decode-time `sum(n_obs) == n_total_alleles` check
        // makes this unreachable on a block that passed `decode_block`; kept as a
        // defensive guard.
        if obs_end > self.obs_seq.len() {
            return Some(Err(PspReadError::BlockStructureInvalid {
                context: "ssr n-obs over-runs the decoded observations",
            }));
        }

        let mut observed = Vec::with_capacity(n_here);
        for j in obs_start..obs_end {
            let count = match u32::try_from(self.obs_count[j]) {
                Ok(c) => c,
                Err(_) => {
                    return Some(Err(PspReadError::ColumnElementDecode {
                        column: "obs-count".to_string(),
                        entry: j,
                        source: super::errors::ScalarDecodeError::VarintOverflow,
                    }));
                }
            };
            observed.push((self.obs_seq[j].clone().into_boxed_slice(), count));
        }

        self.last_start = start;
        self.next_record += 1;
        self.next_obs += n_here as u32;
        Some(Ok(SsrLocusRecord {
            chrom_id: self.chrom_id,
            start,
            end,
            depth: self.depth[i],
            n_filtered: self.n_filtered[i],
            mapped_reads: self.mapped_reads[i],
            n_low_quality: self.n_low_quality[i],
            n_border_off_end: self.n_border_off_end[i],
            observed,
        }))
    }

    fn unload(&mut self) {
        // Free the exhausted block's columns and mark the decoder unloaded.
        self.loaded = false;
        self.delta_start = Vec::new();
        self.span = Vec::new();
        self.n_obs = Vec::new();
        self.depth = Vec::new();
        self.n_filtered = Vec::new();
        self.mapped_reads = Vec::new();
        self.n_low_quality = Vec::new();
        self.n_border_off_end = Vec::new();
        self.obs_count = Vec::new();
        self.obs_seq_len = Vec::new();
        self.obs_seq = Vec::new();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::psp::PspReader;
    use crate::psp::test_fixtures::writer_header;
    use crate::psp::writer::PspWriter;
    use std::io::Cursor;

    fn obs(pairs: &[(&[u8], u32)]) -> Vec<(Box<[u8]>, u32)> {
        pairs
            .iter()
            .map(|(s, c)| (s.to_vec().into_boxed_slice(), *c))
            .collect()
    }

    fn rec(chrom_id: u32, start: u32, end: u32, observed: Vec<(Box<[u8]>, u32)>) -> SsrLocusRecord {
        SsrLocusRecord {
            chrom_id,
            start,
            end,
            depth: 30,
            n_filtered: 2,
            mapped_reads: 33,
            n_low_quality: 1,
            n_border_off_end: 3,
            observed,
        }
    }

    /// SSR columns are well-formed: unique ascending tags, total key coverage.
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

    /// Full round-trip: a multi-block `.ssr.psp` written via the generic writer,
    /// read back through the typed reader, per-record equal. Covers a locus with
    /// zero observations, one with several, and a chromosome boundary.
    #[test]
    fn ssr_round_trips_through_the_container() {
        let records = vec![
            rec(0, 50, 56, obs(&[(b"CACACA", 12)])),
            rec(
                0,
                60,
                70,
                obs(&[(b"CACACA", 8), (b"CACAACA", 1), (b"CACACACA", 3)]),
            ),
            rec(0, 150, 162, obs(&[])), // a locus with no observations
            rec(0, 270, 276, obs(&[(b"GTGTGT", 20)])),
            rec(1, 5, 11, obs(&[(b"AT", 4), (b"ATAT", 2)])),
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

        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.header().kind, "ssr");
        assert!(reader.block_index().len() >= 4);

        let read_back: Vec<SsrLocusRecord> = reader
            .records_of::<SsrKind>()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(read_back, records);
    }

    /// M1: a locus whose half-open interval ends one past the contig's last base
    /// (`end == length + 1`) must round-trip.
    #[test]
    fn ssr_locus_at_contig_end_round_trips() {
        let len = 1_000_000u32;
        let records = vec![rec(0, len - 5, len + 1, obs(&[(b"CACACA", 10)]))];
        let sink = Cursor::new(Vec::<u8>::new());
        let mut writer = PspWriter::<_, SsrKind>::new_ssr(sink, writer_header(2)).unwrap();
        for r in &records {
            writer.write_locus(r).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();
        let mut reader =
            PspReader::new(Cursor::new(bytes)).expect("end-of-contig locus must reopen");
        let read_back: Vec<SsrLocusRecord> = reader
            .records_of::<SsrKind>()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(read_back, records);
    }

    /// M4: reading an `.ssr.psp` through the default SNP path must yield a single
    /// `KindMismatch` then poison — not silently decode SSR columns as SNP.
    #[test]
    fn records_of_wrong_kind_yields_kind_mismatch_then_poisons() {
        let records = vec![rec(0, 50, 56, obs(&[(b"CACACA", 12)]))];
        let sink = Cursor::new(Vec::<u8>::new());
        let mut writer = PspWriter::<_, SsrKind>::new_ssr(sink, writer_header(2)).unwrap();
        for r in &records {
            writer.write_locus(r).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();
        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let mut it = reader.records(); // defaults to SnpKind; file kind is "ssr"
        match it.next() {
            Some(Err(PspReadError::KindMismatch { expected, found })) => {
                assert_eq!(expected, "snp");
                assert_eq!(found, "ssr");
            }
            other => panic!("expected KindMismatch, got {other:?}"),
        }
        assert!(
            it.next().is_none(),
            "iterator must poison after the mismatch"
        );
    }
}
