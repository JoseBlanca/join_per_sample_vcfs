//! Per-sample reading (appendix §A) — `SamplePspReader` + `SamplePspChunk`.
//!
//! *(today: `var_calling::column_span_reader`, backed by `psp::PspReader`)*
//!
//! All per-sample responsibility lives here:
//!
//! - [`SamplePspReader`] — one per sample, created from a
//!   [`psp::PspReader`](crate::psp::PspReader) + a chromosome region.
//!   **No dust** (the producer is the sole dust consumer, §2.3). Hands out
//!   one [`SamplePspChunk`] per psp **segment** (the natural read unit), via
//!   [`next_chunk`](SamplePspReader::next_chunk); the producer takes
//!   `min(peek_next_span)` across readers to advance the cohort in lockstep.
//!   One reader is owned per sample; the producer decodes the per-sample
//!   readers in parallel via disjoint `&mut` (`par_iter_mut`), so `R: Send`
//!   is required (appendix §A).
//! - [`SamplePspChunk`] — one sample's columns for **one psp segment**.
//!   `new` decodes the *light* fold columns eagerly
//!   ([`positions`](SamplePspChunk::positions) /
//!   [`nonref_obs`](SamplePspChunk::nonref_obs) /
//!   [`ref_spans`](SamplePspChunk::ref_spans), cached); the typed getters
//!   ([`take_seq`](SamplePspChunk::take_seq) /
//!   [`take_chain_ids`](SamplePspChunk::take_chain_ids) /
//!   [`take_fixed`](SamplePspChunk::take_fixed)) move out their heavy
//!   column(s) for the `keep` rows only.
//!
//! Phase 1 uses a **simple decode**: `SamplePspReader` decodes the whole
//! segment up front (via the shared [`BlockColumnReader`]) and the chunk
//! owns a copy of every column; the take-getters slice the `keep` rows out
//! of that copy. The column-selective *skip-the-rest* lever (decode only the
//! light columns to decide variable, then only the heavy columns of the
//! variable rows, freeing the compressed bytes) is Phase 5 — the take-getter
//! signatures (`&mut self`, move-out, once-only) are already shaped for it.

use std::io::{Read, Seek};

use crate::pileup_record::{AlleleObservation, AlleleSupportStats, ChainId, PileupRecord};
use crate::psp::ScalarDecodeError;
use crate::psp::{BlockColumnReader, BlockColumns, BlockIndexEntry, PspReadError, PspReader};

// ---------------------------------------------------------------------------
// Heavy per-allele column carriers (copied from the old columnar reader).
//
// These are the byte-identity-sensitive destination shapes the appendix §A
// take-getters return. `var_calling::columns` was deliberately not
// transplanted (the columnar EM block is gone), and referencing it would
// break at the P7 swap when the old package is deleted — so the minimal
// subset the reader needs is copied here.
// ---------------------------------------------------------------------------

/// Per-(record, allele) fixed-width scalar columns — the 7 components of
/// [`AlleleSupportStats`], one cell per allele.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PerAlleleFixed {
    pub num_obs: Vec<u32>,
    pub q_sum: Vec<f64>,
    pub fwd: Vec<u32>,
    pub placed_left: Vec<u32>,
    pub placed_start: Vec<u32>,
    pub mapq_sum: Vec<u32>,
    pub mapq_sum_sq: Vec<u64>,
}

impl PerAlleleFixed {
    /// Total per-allele cells across all records.
    pub fn len(&self) -> usize {
        self.num_obs.len()
    }

    /// True when no per-allele cells are stored.
    pub fn is_empty(&self) -> bool {
        self.num_obs.is_empty()
    }

    /// Append `src`'s cells in `range` to `self`.
    pub fn extend_from_range(&mut self, src: &Self, range: std::ops::Range<usize>) {
        self.num_obs.extend_from_slice(&src.num_obs[range.clone()]);
        self.q_sum.extend_from_slice(&src.q_sum[range.clone()]);
        self.fwd.extend_from_slice(&src.fwd[range.clone()]);
        self.placed_left
            .extend_from_slice(&src.placed_left[range.clone()]);
        self.placed_start
            .extend_from_slice(&src.placed_start[range.clone()]);
        self.mapq_sum
            .extend_from_slice(&src.mapq_sum[range.clone()]);
        self.mapq_sum_sq.extend_from_slice(&src.mapq_sum_sq[range]);
    }

    /// Build an [`AlleleSupportStats`] from the cell at `allele_idx`.
    pub fn support_at(&self, allele_idx: usize) -> AlleleSupportStats {
        AlleleSupportStats::new(
            self.num_obs[allele_idx],
            self.q_sum[allele_idx],
            self.fwd[allele_idx],
            self.placed_left[allele_idx],
            self.placed_start[allele_idx],
            self.mapq_sum[allele_idx],
            self.mapq_sum_sq[allele_idx],
        )
    }
}

/// Per-allele variable-length sequence-bytes column, nested-CSR.
/// `offsets[k]..offsets[k + 1]` is allele `k`'s slice of [`Self::bytes`];
/// the leading `0` keeps `offsets.len() == n_alleles + 1`.
#[derive(Debug, Clone, PartialEq)]
pub struct PerAlleleSeq {
    pub offsets: Vec<u32>,
    pub bytes: Vec<u8>,
}

impl PerAlleleSeq {
    /// Empty columns, CSR sentinel seeded.
    pub fn empty() -> Self {
        Self {
            offsets: vec![0],
            bytes: Vec::new(),
        }
    }

    /// Append one allele's sequence bytes and push the CSR offset.
    pub fn push(&mut self, seq: &[u8]) {
        self.bytes.extend_from_slice(seq);
        self.offsets.push(u32_from_usize(self.bytes.len()));
    }

    /// Append `src`'s alleles in `allele_range` to `self`, copying bytes and
    /// rewriting offsets `self`-relative.
    pub fn extend_from_range(&mut self, src: &Self, allele_range: std::ops::Range<usize>) {
        let seq_lo = src.offsets[allele_range.start] as usize;
        let seq_hi = src.offsets[allele_range.end] as usize;
        let bytes_base = self.bytes.len();
        self.bytes.extend_from_slice(&src.bytes[seq_lo..seq_hi]);
        for k in (allele_range.start + 1)..=allele_range.end {
            let inner_offset = src.offsets[k] as usize - seq_lo;
            self.offsets.push(u32_from_usize(bytes_base + inner_offset));
        }
    }

    /// Number of alleles stored.
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// True when no alleles are stored.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Allele `allele_idx`'s sequence bytes.
    pub fn slice_at(&self, allele_idx: usize) -> &[u8] {
        let lo = self.offsets[allele_idx] as usize;
        let hi = self.offsets[allele_idx + 1] as usize;
        &self.bytes[lo..hi]
    }
}

/// Per-allele variable-length chain-id column, nested-CSR (mirrors
/// [`PerAlleleSeq`] for [`ChainId`]s).
#[derive(Debug, Clone, PartialEq)]
pub struct PerAlleleChainIds {
    pub offsets: Vec<u32>,
    pub ids: Vec<ChainId>,
}

impl PerAlleleChainIds {
    /// Empty columns, CSR sentinel seeded.
    pub fn empty() -> Self {
        Self {
            offsets: vec![0],
            ids: Vec::new(),
        }
    }

    /// Append one allele's chain ids and push the CSR offset.
    pub fn push(&mut self, chain_ids: &[ChainId]) {
        self.ids.extend_from_slice(chain_ids);
        self.offsets.push(u32_from_usize(self.ids.len()));
    }

    /// Append `src`'s alleles in `allele_range` to `self`.
    pub fn extend_from_range(&mut self, src: &Self, allele_range: std::ops::Range<usize>) {
        let id_lo = src.offsets[allele_range.start] as usize;
        let id_hi = src.offsets[allele_range.end] as usize;
        let ids_base = self.ids.len();
        self.ids.extend_from_slice(&src.ids[id_lo..id_hi]);
        for k in (allele_range.start + 1)..=allele_range.end {
            let inner_offset = src.offsets[k] as usize - id_lo;
            self.offsets.push(u32_from_usize(ids_base + inner_offset));
        }
    }

    /// Number of alleles stored.
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// True when no alleles are stored.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Allele `allele_idx`'s chain ids.
    pub fn slice_at(&self, allele_idx: usize) -> &[ChainId] {
        let lo = self.offsets[allele_idx] as usize;
        let hi = self.offsets[allele_idx + 1] as usize;
        &self.ids[lo..hi]
    }
}

fn u32_from_usize(v: usize) -> u32 {
    debug_assert!(v <= u32::MAX as usize, "column offset {v} exceeds u32");
    v as u32
}

// ---------------------------------------------------------------------------
// SamplePspChunk
// ---------------------------------------------------------------------------

/// One sample's columns for **one psp segment**, clamped to the reader's
/// region (appendix §A `[NEW]`).
///
/// `new` decodes and caches the light fold columns
/// ([`positions`](Self::positions) / [`nonref_obs`](Self::nonref_obs) /
/// [`ref_spans`](Self::ref_spans), one value per record); the heavy
/// per-allele columns are owned and handed out by the move-out take-getters
/// for the `keep` rows. Record indices used by the getters and by `keep` are
/// `0..len()` over this chunk's (region-clamped) records.
#[derive(Debug, Clone, PartialEq)]
pub struct SamplePspChunk {
    chrom_id: u32,
    /// 1-based positions, one per record. Light, cached at `new`.
    positions: Vec<u32>,
    /// Per-record sum of non-REF allele obs (alleles `[1..]`). Light.
    nonref_obs: Vec<u32>,
    /// Per-record REF-allele (`alleles[0]`) reference span in bases. Light.
    ref_spans: Vec<u32>,
    /// CSR record→allele offsets (`len() == n_records + 1`); maps a record
    /// index to its allele range in the heavy columns.
    allele_offsets: Vec<u32>,
    /// Heavy columns, owned; moved out (and emptied) by the take-getters.
    fixed: PerAlleleFixed,
    seq: PerAlleleSeq,
    chain_ids: PerAlleleChainIds,
}

impl SamplePspChunk {
    /// Decode one block's region-clamped records into an owned chunk.
    ///
    /// Records with `pos < floor` or `pos > region_end` are dropped (only
    /// the region-edge blocks carry any). Returns `None` when no record of
    /// the block falls in `[floor, region_end]`.
    fn from_block(
        cols: &BlockColumns<'_>,
        chrom_id: u32,
        floor: u32,
        region_end: u32,
    ) -> Result<Option<Self>, PspReadError> {
        let n = cols.n_records as usize;

        // Absolute positions for the whole block (delta-decoded).
        let mut abs = Vec::with_capacity(n);
        let mut p = cols.first_pos;
        for i in 0..n {
            if i > 0 {
                p = advance_pos(p, cols.delta_pos[i], i)?;
            }
            abs.push(p);
        }

        // Region clamp: ascending positions ⇒ a contiguous record range.
        let r_lo = abs.partition_point(|&pos| pos < floor);
        let r_hi = abs.partition_point(|&pos| pos <= region_end);
        if r_lo >= r_hi {
            return Ok(None);
        }

        // Allele range of the clamped record range [r_lo, r_hi): `a_lo` is
        // the cumulative allele count before record `r_lo`.
        let a_lo: usize = cols.n_alleles[..r_lo].iter().map(|&n| n as usize).sum();
        let a_hi = a_lo
            + cols.n_alleles[r_lo..r_hi]
                .iter()
                .map(|&n| n as usize)
                .sum::<usize>();

        // CSR record→allele offsets, chunk-relative.
        let mut allele_offsets = Vec::with_capacity(r_hi - r_lo + 1);
        allele_offsets.push(0u32);
        let mut cum = 0u32;
        for &na in &cols.n_alleles[r_lo..r_hi] {
            cum += na as u32;
            allele_offsets.push(cum);
        }

        // Heavy columns for the clamped allele range.
        let mut fixed = PerAlleleFixed::default();
        fixed
            .num_obs
            .extend_from_slice(&cols.allele_obs_count[a_lo..a_hi]);
        fixed
            .q_sum
            .extend_from_slice(&cols.allele_q_sum_log[a_lo..a_hi]);
        fixed
            .fwd
            .extend_from_slice(&cols.allele_fwd_count[a_lo..a_hi]);
        fixed
            .placed_left
            .extend_from_slice(&cols.allele_placed_left_count[a_lo..a_hi]);
        fixed
            .placed_start
            .extend_from_slice(&cols.allele_placed_start_count[a_lo..a_hi]);
        fixed
            .mapq_sum
            .extend_from_slice(&cols.allele_mapq_sum[a_lo..a_hi]);
        fixed
            .mapq_sum_sq
            .extend_from_slice(&cols.allele_mapq_sum_sq[a_lo..a_hi]);

        let mut seq = PerAlleleSeq::empty();
        let mut chain_ids = PerAlleleChainIds::empty();
        for k in a_lo..a_hi {
            let s = cols.allele_seq_offsets[k] as usize;
            let e = cols.allele_seq_offsets[k + 1] as usize;
            seq.push(&cols.allele_seq_data[s..e]);
            let cs = cols.allele_chain_ids_offsets[k] as usize;
            let ce = cols.allele_chain_ids_offsets[k + 1] as usize;
            chain_ids.push(&cols.allele_chain_ids_data[cs..ce]);
        }

        // Light per-record columns over the clamped range.
        let positions = abs[r_lo..r_hi].to_vec();
        let mut nonref_obs = Vec::with_capacity(r_hi - r_lo);
        let mut ref_spans = Vec::with_capacity(r_hi - r_lo);
        let mut base = a_lo;
        for r in r_lo..r_hi {
            let lo = base;
            let hi = lo + cols.n_alleles[r] as usize;
            base = hi;
            // Non-REF obs: sum over alleles [1..]; REF is allele 0.
            let nro = cols.allele_obs_count[(lo + 1).min(hi)..hi]
                .iter()
                .fold(0u32, |acc, &c| acc.saturating_add(c));
            nonref_obs.push(nro);
            // REF span: byte length of allele 0's sequence.
            let ref_span = cols.allele_seq_offsets[lo + 1] - cols.allele_seq_offsets[lo];
            ref_spans.push(ref_span);
        }

        Ok(Some(Self {
            chrom_id,
            positions,
            nonref_obs,
            ref_spans,
            allele_offsets,
            fixed,
            seq,
            chain_ids,
        }))
    }

    /// Chromosome id of every record in this chunk.
    pub fn chrom_id(&self) -> u32 {
        self.chrom_id
    }

    /// Number of (region-clamped) records.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// True when the chunk has no records.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// 1-based positions, one per record (light, cached).
    pub fn positions(&self) -> &[u32] {
        &self.positions
    }

    /// Per-record sum of non-REF allele obs (light, cached). The
    /// per-sample input to the cohort `min_alt_obs` fold.
    pub fn nonref_obs(&self) -> &[u32] {
        &self.nonref_obs
    }

    /// Per-record REF-allele reference span in bases (light, cached). The
    /// per-sample reach input to the cohort fold / grouping.
    pub fn ref_spans(&self) -> &[u32] {
        &self.ref_spans
    }

    /// Number of alleles at record `record_idx`.
    pub fn n_alleles_at(&self, record_idx: usize) -> usize {
        let lo = self.allele_offsets[record_idx] as usize;
        let hi = self.allele_offsets[record_idx + 1] as usize;
        hi - lo
    }

    /// Reconstruct owned [`PileupRecord`]s for the `keep` records — the
    /// columnar→record boundary at the producer (§2.2 step 3). **Non-consuming**
    /// (copies via `slice_at`/`support_at`), so a segment straddling a chunk
    /// cut survives in the buffer to serve the next chunk; the move-out
    /// [`take_seq`](Self::take_seq) / `take_*` getters are the Phase-5
    /// column-selective seam instead. `keep.len()` must equal [`len`](Self::len).
    pub fn records_for(&self, keep: &[bool]) -> Vec<PileupRecord> {
        debug_assert_eq!(keep.len(), self.len(), "keep mask must cover every record");
        let mut out = Vec::new();
        for (r, &k) in keep.iter().enumerate() {
            if !k {
                continue;
            }
            let lo = self.allele_offsets[r] as usize;
            let hi = self.allele_offsets[r + 1] as usize;
            let mut alleles = Vec::with_capacity(hi - lo);
            for a in lo..hi {
                alleles.push(AlleleObservation::new(
                    self.seq.slice_at(a).to_vec(),
                    self.fixed.support_at(a),
                    self.chain_ids.slice_at(a).to_vec(),
                ));
            }
            out.push(PileupRecord::new(self.chrom_id, self.positions[r], alleles));
        }
        out
    }

    /// Move out the fixed-scalar columns for the `keep` records (appendix
    /// §A typed getter). `keep.len()` must equal [`len`](Self::len);
    /// the column is emptied, so a second call yields nothing.
    pub fn take_fixed(&mut self, keep: &[bool]) -> PerAlleleFixed {
        debug_assert_eq!(keep.len(), self.len(), "keep mask must cover every record");
        let src = std::mem::take(&mut self.fixed);
        let mut out = PerAlleleFixed::default();
        for (r, &k) in keep.iter().enumerate() {
            if k {
                let lo = self.allele_offsets[r] as usize;
                let hi = self.allele_offsets[r + 1] as usize;
                out.extend_from_range(&src, lo..hi);
            }
        }
        out
    }

    /// Move out the allele-sequence column for the `keep` records.
    ///
    /// Call-once (per the appendix §A "fetched at most once" contract): the
    /// column is moved out, so a second call on the same chunk is misuse.
    pub fn take_seq(&mut self, keep: &[bool]) -> PerAlleleSeq {
        debug_assert_eq!(keep.len(), self.len(), "keep mask must cover every record");
        let src = std::mem::replace(&mut self.seq, PerAlleleSeq::empty());
        let mut out = PerAlleleSeq::empty();
        for (r, &k) in keep.iter().enumerate() {
            if k {
                let lo = self.allele_offsets[r] as usize;
                let hi = self.allele_offsets[r + 1] as usize;
                out.extend_from_range(&src, lo..hi);
            }
        }
        out
    }

    /// Move out the allele-chain-id column for the `keep` records.
    pub fn take_chain_ids(&mut self, keep: &[bool]) -> PerAlleleChainIds {
        debug_assert_eq!(keep.len(), self.len(), "keep mask must cover every record");
        let src = std::mem::replace(&mut self.chain_ids, PerAlleleChainIds::empty());
        let mut out = PerAlleleChainIds::empty();
        for (r, &k) in keep.iter().enumerate() {
            if k {
                let lo = self.allele_offsets[r] as usize;
                let hi = self.allele_offsets[r + 1] as usize;
                out.extend_from_range(&src, lo..hi);
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// SamplePspReader
// ---------------------------------------------------------------------------

/// Per-sample segment reader bound to one chromosome region (appendix §A
/// `[RENAME ColumnSpanReader]`). Serves the records of `[region_start,
/// region_end]` (1-based inclusive) one **psp segment at a time**, in
/// ascending position order.
pub struct SamplePspReader<R: Read + Seek> {
    blocks: BlockColumnReader<R>,
    chrom_id: u32,
    /// Inclusive upper region bound; records past it are not served.
    region_end: u32,
    /// Inclusive lower region bound; earlier records (region-edge block
    /// only) are skipped.
    floor: u32,
}

impl<R: Read + Seek> SamplePspReader<R> {
    /// Open a reader over `[region_start, region_end]` of `chrom_id`,
    /// taking ownership of `reader`.
    pub fn new(reader: PspReader<R>, chrom_id: u32, region_start: u32, region_end: u32) -> Self {
        let mut blocks = reader.into_column_blocks();
        blocks.seek_to(chrom_id, region_start, region_end);
        Self {
            blocks,
            chrom_id,
            region_end,
            floor: region_start,
        }
    }

    /// Re-point at a new region, reusing the decode buffers (the cohort
    /// producer resets one reader per sample at every covered-interval /
    /// chromosome boundary so the zstd context + column slabs persist).
    pub fn reset(&mut self, chrom_id: u32, region_start: u32, region_end: u32) {
        self.blocks.seek_to(chrom_id, region_start, region_end);
        self.chrom_id = chrom_id;
        self.region_end = region_end;
        self.floor = region_start;
    }

    /// The sample's on-disk block index — the producer unions these across
    /// samples to find the covered intervals per chromosome.
    pub fn block_index(&self) -> &[BlockIndexEntry] {
        self.blocks.index()
    }

    /// Inclusive end of the next segment this reader serves **without a new
    /// decode** (the next block's `last_pos`), clamped to `region_end`.
    /// `None` once the region is exhausted. The producer takes `min` of
    /// this across samples to advance the cohort in lockstep.
    pub fn peek_next_span(&self) -> Option<u32> {
        let entry = self.blocks.peek_block()?;
        if entry.chrom_id != self.chrom_id || entry.first_pos > self.region_end {
            return None;
        }
        Some(entry.last_pos.min(self.region_end))
    }

    /// Decode and return the next in-region segment as a [`SamplePspChunk`],
    /// advancing the cursor past it. `None` once the region is exhausted.
    /// Blocks whose records all fall below `floor` are skipped transparently.
    pub fn next_chunk(&mut self) -> Result<Option<SamplePspChunk>, PspReadError> {
        loop {
            match self.blocks.peek_block() {
                Some(entry)
                    if entry.chrom_id == self.chrom_id && entry.first_pos <= self.region_end => {}
                _ => return Ok(None),
            }
            let loaded = self.blocks.load_current()?;
            debug_assert!(loaded, "peek_block was Some, so the cursor is in range");
            let chunk = {
                let cols = self.blocks.columns().expect("block just loaded");
                SamplePspChunk::from_block(&cols, self.chrom_id, self.floor, self.region_end)?
            };
            self.blocks.advance();
            match chunk {
                Some(c) => return Ok(Some(c)),
                // Whole block was below the region floor — keep scanning.
                None => continue,
            }
        }
    }
}

/// Advance a delta-encoded position, surfacing a `u32`-overflowing delta as
/// the same typed error the row reader raises (rather than a silent wrap).
fn advance_pos(prev: u32, delta: u64, entry: usize) -> Result<u32, PspReadError> {
    let delta_u32 = u32::try_from(delta).map_err(|_| PspReadError::ColumnElementDecode {
        column: "delta-pos".to_string(),
        entry,
        source: ScalarDecodeError::VarintOverflow,
    })?;
    Ok(prev.saturating_add(delta_u32))
}

#[cfg(test)]
mod tests {
    //! The light accessors and the typed getters must equal a full row
    //! decode (`PspReader::records`) of the same segment — the P1 gate.

    use std::io::Cursor;

    use super::*;
    use crate::pileup_record::PileupRecord;
    use crate::psp::test_fixtures::writer_header;
    use crate::psp::writer::PspWriter;
    use crate::var_calling::test_helpers::{allele, record};

    /// Multi-allele records at ascending positions, dense enough that a
    /// small block target yields several segments.
    fn fixture_records() -> Vec<PileupRecord> {
        (0..200u32)
            .map(|i| {
                let pos = 10 + i * 7;
                match i % 3 {
                    0 => record(pos, vec![allele(b"A", 12, -1.5, &[])]),
                    1 => record(
                        pos,
                        vec![
                            allele(b"A", 9, -1.0, &[]),
                            allele(b"T", 4, -2.0, &[(i as u64) + 1]),
                        ],
                    ),
                    _ => record(
                        pos,
                        vec![
                            allele(b"AC", 7, -0.5, &[]),
                            allele(b"AG", 3, -3.0, &[(i as u64) + 1, (i as u64) + 2]),
                            allele(b"A", 5, -1.2, &[]),
                        ],
                    ),
                }
            })
            .collect()
    }

    fn psp_bytes(records: &[PileupRecord], block_target: usize) -> Vec<u8> {
        let mut writer = PspWriter::new_with_block_target(
            Cursor::new(Vec::new()),
            writer_header(1),
            block_target,
        )
        .expect("writer opens");
        for r in records {
            writer.write_record(r).expect("write_record");
        }
        writer.finish().expect("finish").into_inner()
    }

    /// Reference: full row decode of `[start, end]`.
    fn rows(bytes: &[u8], chrom: u32, start: u32, end: u32) -> Vec<PileupRecord> {
        let mut reader = PspReader::new(Cursor::new(bytes.to_vec())).expect("reader opens");
        reader
            .region_records(chrom, start, end)
            .map(|r| r.expect("record decodes"))
            .collect()
    }

    /// Drain every chunk of a reader over `[start, end]`, concatenating the
    /// light columns and the `keep`-all heavy columns into flat buffers
    /// matched against the row reference.
    fn drain_all(
        bytes: &[u8],
        chrom: u32,
        start: u32,
        end: u32,
    ) -> (
        Vec<u32>,
        Vec<u32>,
        Vec<u32>,
        PerAlleleFixed,
        PerAlleleSeq,
        PerAlleleChainIds,
    ) {
        let reader = PspReader::new(Cursor::new(bytes.to_vec())).expect("reader opens");
        let mut sr = SamplePspReader::new(reader, chrom, start, end);

        let mut positions = Vec::new();
        let mut nonref_obs = Vec::new();
        let mut ref_spans = Vec::new();
        let mut fixed = PerAlleleFixed::default();
        let mut seq = PerAlleleSeq::empty();
        let mut chain = PerAlleleChainIds::empty();

        while let Some(mut chunk) = sr.next_chunk().expect("next_chunk") {
            positions.extend_from_slice(chunk.positions());
            nonref_obs.extend_from_slice(chunk.nonref_obs());
            ref_spans.extend_from_slice(chunk.ref_spans());
            let keep = vec![true; chunk.len()];
            let cf = chunk.take_fixed(&keep);
            let cs = chunk.take_seq(&keep);
            let cc = chunk.take_chain_ids(&keep);
            let n = cf.len();
            fixed.extend_from_range(&cf, 0..n);
            seq.extend_from_range(&cs, 0..cs.len());
            chain.extend_from_range(&cc, 0..cc.len());
        }
        (positions, nonref_obs, ref_spans, fixed, seq, chain)
    }

    /// Build the expected flat columns straight from the row records.
    fn expected_from_rows(
        recs: &[PileupRecord],
    ) -> (
        Vec<u32>,
        Vec<u32>,
        Vec<u32>,
        PerAlleleFixed,
        PerAlleleSeq,
        PerAlleleChainIds,
    ) {
        let mut positions = Vec::new();
        let mut nonref_obs = Vec::new();
        let mut ref_spans = Vec::new();
        let mut fixed = PerAlleleFixed::default();
        let mut seq = PerAlleleSeq::empty();
        let mut chain = PerAlleleChainIds::empty();
        for r in recs {
            positions.push(r.pos);
            let nro = r.alleles[1..]
                .iter()
                .fold(0u32, |acc, a| acc.saturating_add(a.support.num_obs));
            nonref_obs.push(nro);
            ref_spans.push(r.alleles[0].seq.len() as u32);
            for a in &r.alleles {
                fixed.num_obs.push(a.support.num_obs);
                fixed.q_sum.push(a.support.q_sum);
                fixed.fwd.push(a.support.fwd);
                fixed.placed_left.push(a.support.placed_left);
                fixed.placed_start.push(a.support.placed_start);
                fixed.mapq_sum.push(a.support.mapq_sum);
                fixed.mapq_sum_sq.push(a.support.mapq_sum_sq);
                seq.push(&a.seq);
                chain.push(&a.chain_ids);
            }
        }
        (positions, nonref_obs, ref_spans, fixed, seq, chain)
    }

    #[test]
    fn whole_region_matches_row_decode() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        // Sanity: multi-segment.
        let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
        assert!(
            reader.block_index().len() >= 3,
            "fixture should be multi-block"
        );

        let got = drain_all(&bytes, 0, 1, 1_000_000);
        let want = expected_from_rows(&rows(&bytes, 0, 1, 1_000_000));
        assert_eq!(got.0, want.0, "positions");
        assert_eq!(got.1, want.1, "nonref_obs");
        assert_eq!(got.2, want.2, "ref_spans");
        assert_eq!(got.3, want.3, "fixed");
        assert_eq!(got.4, want.4, "seq");
        assert_eq!(got.5, want.5, "chain_ids");
    }

    #[test]
    fn clamped_region_matches_row_decode() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        let (start, end) = (137, 941);
        let got = drain_all(&bytes, 0, start, end);
        let want = expected_from_rows(&rows(&bytes, 0, start, end));
        assert_eq!(got.0, want.0, "positions");
        assert_eq!(got.3, want.3, "fixed");
        assert_eq!(got.4, want.4, "seq");
        assert_eq!(got.5, want.5, "chain_ids");
        // The clamp really dropped records outside the window.
        assert!(got.0.len() < recs.len());
        assert_eq!(*got.0.first().unwrap(), *want.0.first().unwrap());
        assert!(*got.0.first().unwrap() >= start && *got.0.last().unwrap() <= end);
    }

    #[test]
    fn take_with_keep_subset_selects_those_records() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 4096); // single-ish block
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let mut sr = SamplePspReader::new(reader, 0, 1, 1_000_000);

        let mut chunk = sr.next_chunk().expect("next_chunk").expect("a chunk");
        // Keep only even-indexed records of this chunk.
        let keep: Vec<bool> = (0..chunk.len()).map(|i| i % 2 == 0).collect();
        let kept_positions: Vec<u32> = chunk
            .positions()
            .iter()
            .enumerate()
            .filter_map(|(i, &p)| (i % 2 == 0).then_some(p))
            .collect();
        let expected_alleles: usize = (0..chunk.len())
            .filter(|i| i % 2 == 0)
            .map(|i| chunk.n_alleles_at(i))
            .sum();

        let fixed = chunk.take_fixed(&keep);
        let seq = chunk.take_seq(&keep);
        let chain = chunk.take_chain_ids(&keep);
        assert_eq!(fixed.len(), expected_alleles, "fixed allele count");
        assert_eq!(seq.len(), expected_alleles, "seq allele count");
        assert_eq!(chain.len(), expected_alleles, "chain allele count");
        // Cross-check against a fresh decode of just those positions.
        let rows_all = rows(&psp_bytes(&recs, 4096), 0, 1, 1_000_000);
        let want: Vec<u32> = kept_positions.clone();
        let got: Vec<u32> = rows_all
            .iter()
            .filter(|r| want.contains(&r.pos))
            .map(|r| r.pos)
            .collect();
        assert_eq!(got, want);
    }

    #[test]
    fn peek_reports_ascending_segment_ends() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let mut sr = SamplePspReader::new(reader, 0, 1, 1_000_000);
        let mut peeks = Vec::new();
        while let Some(end) = sr.peek_next_span() {
            peeks.push(end);
            sr.next_chunk().unwrap();
        }
        assert!(peeks.len() >= 3, "multi-segment");
        assert!(
            peeks.windows(2).all(|w| w[0] < w[1]),
            "segment ends ascending"
        );
    }

    #[test]
    fn empty_region_yields_no_chunks() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        // Positions are 10, 17, 24, ...; (11, 16) holds none.
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let mut sr = SamplePspReader::new(reader, 0, 11, 16);
        assert!(sr.next_chunk().unwrap().is_none());
    }
}
