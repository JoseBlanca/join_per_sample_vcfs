//! Streaming reader over a `.psp` byte source.
//!
//! Composes the decoder primitives in [`super::header`],
//! [`super::trailer`], [`super::index`], and [`super::block`] into
//! an end-to-end open + iterate pipeline. The orchestration
//! responsibilities — not in the primitives — are:
//!
//! - **Open sequence:** tail-first read of the trailer, decode of
//!   the block index (with XXH3-64 checksum check), then head-first
//!   read + parse of the header. Each step gates the next on data
//!   already known to be well-formed.
//! - **Cross-region cross-checks:** every `block_offset` falls in
//!   `[header_end, trailer.index_offset)` and is strictly
//!   increasing; every `chrom_id` is `< header.chromosomes.len()`;
//!   every `first_pos` / `last_pos` falls in
//!   `[1, chromosomes[chrom_id].length]`.
//! - **Sequential iteration:** per-block decode + per-record
//!   materialisation with running phase-chain bookkeeping and the
//!   block-start snapshot continuity check (spec check #11).
//! - **Region iteration:** binary-search the block index for the
//!   first overlapping block; snapshot-init the active-slot set at
//!   each block (no cross-block continuity check on random-access
//!   reads); clamp records to the requested window.
//!
//! See `ia/feature_implementation_plans/psp_reader.md` for the
//! design rationale and the test plan (Groups B, F, and R).

use std::collections::BTreeSet;
use std::io::{Read, Seek, SeekFrom};

use super::block::{
    BlockHeader, ColumnManifestEntry, decode_block_header, decode_bytes_split, decode_list_column,
    decode_scalar_column, decode_varint_column, zstd_decompress,
};
use super::errors::{
    BlockHeaderInvariantKind, PhaseChainConsistencyKind, PspReadError, VarintError,
};
use super::header::{
    HEAD_SENTINEL_LEN, HEADER_FRAMING_BYTES, MAX_HEADER_BODY_BYTES, MIN_HEADER_BODY_BYTES,
    ParsedHeader, parse_header_bytes,
};
use super::index::{BlockIndexEntry, checksum_index, decode_index};
use super::registry::{ColumnKey, MAX_ALLELE_SEQ_LEN, lookup_by_tag};
use super::trailer::{TRAILER_BYTES, Trailer, decode_trailer};
use crate::per_sample_caller::pileup::{
    AlleleObservation, AlleleSupportStats, PileupRecord, SlotId,
};

/// Cap on how many bytes the reader will pull off the source in
/// pursuit of a single block header before giving up. The header is
/// a varint stream with no outer length prefix, so the decoder is a
/// grow-on-incomplete loop; this cap stops a malformed file from
/// inducing an unbounded read. 64 KiB is far more than any sane
/// block header (the writer's own headers run in the hundreds of
/// bytes even with 12 columns).
const BLOCK_HEADER_READ_CAP: usize = 64 * 1024;

/// Initial chunk size for the block-header grow-on-incomplete read
/// loop. Doubles up to [`BLOCK_HEADER_READ_CAP`] before giving up.
const BLOCK_HEADER_INITIAL_CHUNK: usize = 4 * 1024;

/// Reader for a `.psp` byte source. Constructed by [`Self::new`],
/// which performs the open sequence (trailer → block index →
/// header) and the cross-region integrity checks. The owned source
/// is positioned at end-of-header on success — ready for an
/// immediate sequential [`Self::records`] call.
pub struct PspReader<R: Read + Seek> {
    source: R,
    header: ParsedHeader,
    trailer: Trailer,
    /// Decoded block index. Empty for a zero-block file.
    index: Vec<BlockIndexEntry>,
    /// Absolute file offset of the first byte past the framed
    /// header. Used by `new` as the lower bound for `block_offset`
    /// validation; Slice 5's region init may also use it.
    #[allow(dead_code)]
    header_end_offset: u64,
}

impl<R: Read + Seek> std::fmt::Debug for PspReader<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PspReader")
            .field("sample", &self.header.sample)
            .field("n_blocks", &self.index.len())
            .field("trailer", &self.trailer)
            .finish()
    }
}

impl<R: Read + Seek> PspReader<R> {
    /// Open a `.psp` source. Performs the full validation pipeline
    /// described in spec §"Header-binary consistency: required
    /// reader checks", steps 1–6 (header), the XXH3-64 index
    /// checksum, and the block-index cross-region checks
    /// (block offsets in range + strictly increasing; chrom_id and
    /// position within the declared chromosome table).
    pub fn new(mut source: R) -> Result<Self, PspReadError> {
        // 1. Size the file. A `.psp` cannot be shorter than the
        //    framed header plus the 32-byte trailer.
        let file_len = source.seek(SeekFrom::End(0)).map_err(io_err("file size"))?;
        let min_len = (HEADER_FRAMING_BYTES + TRAILER_BYTES) as u64;
        if file_len < min_len {
            return Err(PspReadError::Io {
                context: "trailer",
                source: std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    format!("file is {file_len} bytes; minimum for a valid .psp is {min_len}"),
                ),
            });
        }

        // 2. Read + decode the trailer.
        source
            .seek(SeekFrom::End(-(TRAILER_BYTES as i64)))
            .map_err(io_err("trailer"))?;
        let mut trailer_bytes = [0u8; TRAILER_BYTES];
        source
            .read_exact(&mut trailer_bytes)
            .map_err(io_err("trailer"))?;
        let trailer = decode_trailer(&trailer_bytes)?;

        // Cross-check trailer arithmetic against the file size. The
        // index plus the 32-byte trailer must fit in the file
        // exactly — anything else means the trailer's pointer math
        // is wrong (corruption) or a write was truncated.
        let trailer_start = file_len - TRAILER_BYTES as u64;
        let computed_end = trailer.index_offset.checked_add(trailer.index_byte_length);
        if computed_end != Some(trailer_start) {
            return Err(PspReadError::IndexTrailingBytes {
                trailing_bytes: trailer_start
                    .saturating_sub(trailer.index_offset + trailer.index_byte_length)
                    as usize,
            });
        }

        // 3. Read + decode the block index, then verify the XXH3-64
        //    checksum the trailer stamped over it.
        source
            .seek(SeekFrom::Start(trailer.index_offset))
            .map_err(io_err("block index"))?;
        let mut index_bytes = vec![0u8; trailer.index_byte_length as usize];
        source
            .read_exact(&mut index_bytes)
            .map_err(io_err("block index"))?;
        let computed_checksum = checksum_index(&index_bytes);
        if computed_checksum != trailer.index_checksum {
            return Err(PspReadError::IndexChecksum {
                stored: trailer.index_checksum,
                computed: computed_checksum,
            });
        }
        let index = decode_index(&index_bytes, trailer.n_blocks)?;

        // 4. Sanity-check `block_offset`s: every entry's offset must
        //    be strictly less than `trailer.index_offset` (blocks
        //    end before the index starts), and the sequence must be
        //    strictly increasing. The lower-bound check (≥
        //    end-of-header) happens after the header parse in
        //    step 6.
        for (i, entry) in index.iter().enumerate() {
            if entry.block_offset >= trailer.index_offset {
                return Err(PspReadError::BlockIndexOffsetInvalid {
                    block: i,
                    offset: entry.block_offset,
                    min: 0,
                    max: trailer.index_offset,
                });
            }
            if i > 0 {
                let prev = &index[i - 1];
                if entry.block_offset <= prev.block_offset {
                    return Err(PspReadError::BlockIndexOffsetInvalid {
                        block: i,
                        offset: entry.block_offset,
                        min: prev.block_offset + 1,
                        max: trailer.index_offset,
                    });
                }
            }
        }

        // 5. Read + parse the framed header.
        source.seek(SeekFrom::Start(0)).map_err(io_err("header"))?;
        let mut head_prefix = [0u8; 12];
        source
            .read_exact(&mut head_prefix)
            .map_err(io_err("header magic + length prefix"))?;
        // Range-check the body length before allocating, so a
        // tampered length cannot drive an unbounded read.
        let body_len = u64::from_le_bytes(head_prefix[4..12].try_into().unwrap());
        if !(MIN_HEADER_BODY_BYTES..=MAX_HEADER_BODY_BYTES).contains(&body_len) {
            return Err(PspReadError::BadHeaderLength {
                got: body_len,
                min: MIN_HEADER_BODY_BYTES,
                max: MAX_HEADER_BODY_BYTES,
            });
        }
        let body_and_sentinel_len = body_len as usize + HEAD_SENTINEL_LEN;
        let mut header_bytes = Vec::with_capacity(12 + body_and_sentinel_len);
        header_bytes.extend_from_slice(&head_prefix);
        header_bytes.resize(12 + body_and_sentinel_len, 0);
        source
            .read_exact(&mut header_bytes[12..])
            .map_err(io_err("header body + sentinel"))?;
        let (header, header_end_usize) = parse_header_bytes(&header_bytes)?;
        let header_end_offset = header_end_usize as u64;

        // 6. Cross-bind the index against the parsed header.
        if let Some(first) = index.first()
            && first.block_offset < header_end_offset
        {
            return Err(PspReadError::BlockIndexOffsetInvalid {
                block: 0,
                offset: first.block_offset,
                min: header_end_offset,
                max: trailer.index_offset,
            });
        }
        let n_chroms = header.chromosomes.len() as u32;
        for (i, entry) in index.iter().enumerate() {
            if entry.chrom_id >= n_chroms {
                return Err(PspReadError::BlockIndexChromOutOfRange {
                    block: i,
                    chrom_id: entry.chrom_id,
                    n_chroms,
                });
            }
            let chrom = &header.chromosomes[entry.chrom_id as usize];
            for pos in [entry.first_pos, entry.last_pos] {
                if pos == 0 || pos > chrom.length {
                    return Err(PspReadError::BlockIndexPosOutOfRange {
                        block: i,
                        chrom_id: entry.chrom_id,
                        pos,
                        chromosome_length: chrom.length,
                    });
                }
            }
            if entry.first_pos > entry.last_pos {
                return Err(PspReadError::BlockIndexPosOutOfRange {
                    block: i,
                    chrom_id: entry.chrom_id,
                    pos: entry.first_pos,
                    chromosome_length: chrom.length,
                });
            }
        }

        // 7. Position the cursor at the start of block 0 so a
        //    subsequent sequential `records()` call needs no extra
        //    seek. For a zero-block file this is end-of-header
        //    (== trailer.index_offset).
        source
            .seek(SeekFrom::Start(header_end_offset))
            .map_err(io_err("post-header seek"))?;

        Ok(Self {
            source,
            header,
            trailer,
            index,
            header_end_offset,
        })
    }

    /// Parsed file header. Available immediately after [`Self::new`]
    /// returns — does not require any block reads.
    pub fn header(&self) -> &ParsedHeader {
        &self.header
    }

    /// Decoded block index, in genomic order. Empty for a
    /// zero-block file. Exposed for tests, the future `psp dump`
    /// utility, and the cohort-merge driver that wants per-sample
    /// block layout up front (e.g. to schedule reads).
    pub fn block_index(&self) -> &[BlockIndexEntry] {
        &self.index
    }

    /// Trailer in strong-typed form. Hidden from rustdoc — the
    /// trailer fields are derived from `block_index()` plus the
    /// file size, and downstream consumers should not need them.
    /// Kept for tests and the future `psp dump`.
    #[doc(hidden)]
    pub fn trailer(&self) -> Trailer {
        self.trailer
    }

    /// Sequential iterator over every record in genomic order.
    ///
    /// Slice 2 lands the single-block decode loop; Slice 3 the
    /// multi-block continuation; Slice 4 the phase-chain
    /// bookkeeping. Region iteration ships in Slice 5.
    pub fn records(&mut self) -> RecordsIter<'_, R> {
        RecordsIter::new(self, RangeClamp::None)
    }

    /// Iterator over records inside `[start, end]` (inclusive
    /// both ends) on `chrom_id`.
    ///
    /// Position the iterator at the first block whose range
    /// overlaps the window, init the active-slot set from that
    /// block's snapshot (random-access mode skips spec check #11,
    /// per the plan), and let `RecordsIter` clamp per-record:
    /// records before `start` are skipped; the first record past
    /// `end` ends iteration. Empty window or no-overlap returns
    /// `None` from the first `next()` call.
    pub fn region_records(&mut self, chrom_id: u32, start: u32, end: u32) -> RecordsIter<'_, R> {
        let first = first_block_overlapping(&self.index, chrom_id, start, end);
        let clamp = RangeClamp::Window {
            chrom_id,
            start,
            end,
        };
        let mut it = RecordsIter::new(self, clamp);
        match first {
            Some(idx) => it.cur_block_idx = idx,
            None => it.cur_block_idx = it.reader.index.len(),
        }
        it
    }
}

/// Binary-search the block index for the first entry whose
/// `chrom_id` matches and whose `[first_pos, last_pos]` overlaps
/// `[start, end]`. Returns `None` if no entry overlaps.
fn first_block_overlapping(
    index: &[BlockIndexEntry],
    chrom_id: u32,
    start: u32,
    end: u32,
) -> Option<usize> {
    // Block index is sorted by `(chrom_id, first_pos)`. The
    // first block whose range overlaps `[start, end]` is the
    // first one with `chrom_id == chrom_id && last_pos >= start
    // && first_pos <= end`.
    //
    // `partition_point` lands at the first entry whose
    // `(chrom_id, first_pos)` strictly exceeds `(chrom_id,
    // start)`. The overlap candidate is at `pivot - 1` (its
    // `first_pos ≤ start`, possibly with `last_pos ≥ start`) or
    // at `pivot` (`first_pos > start`, may still be `≤ end`).
    // Check both, then fall forward.
    let pivot = index.partition_point(|e| (e.chrom_id, e.first_pos) <= (chrom_id, start));
    let candidates = pivot.saturating_sub(1)..(pivot + 1).min(index.len());
    for i in candidates {
        let e = &index[i];
        if e.chrom_id == chrom_id && e.last_pos >= start && e.first_pos <= end {
            return Some(i);
        }
    }
    // No candidate at the pivot pair — scan forward in case the
    // window starts before any block on this chromosome.
    for (i, e) in index.iter().enumerate().skip(pivot) {
        if e.chrom_id != chrom_id || e.first_pos > end {
            return None;
        }
        if e.last_pos >= start {
            return Some(i);
        }
    }
    None
}

/// Materialised columns for one decoded block. Constructed by
/// [`RecordsIter::load_next_block`]; consumed record-by-record by
/// [`RecordsIter::next`].
#[derive(Debug)]
struct DecodedBlock {
    chrom_id: u32,
    first_pos: u32,
    n_records: u32,
    #[allow(dead_code)] // sanity-checked against the manifest in Slice 3+
    n_total_alleles: u32,
    /// Active phase-chain slots at block start; consumed by
    /// Slice 4's continuity check.
    #[allow(dead_code)]
    active_chain_slots_snapshot: Vec<SlotId>,
    // Per-record columns
    delta_pos: Vec<u64>,
    n_alleles: Vec<u64>,
    new_chain_slots: Vec<Vec<SlotId>>,
    expired_chain_slots: Vec<Vec<SlotId>>,
    // Per-allele columns
    allele_seqs: Vec<Vec<u8>>,
    allele_obs_count: Vec<u32>,
    allele_q_sum_log: Vec<f64>,
    allele_fwd_count: Vec<u32>,
    allele_placed_left_count: Vec<u32>,
    allele_placed_start_count: Vec<u32>,
    allele_chain_slots: Vec<Vec<SlotId>>,
}

/// Region-iteration window. `None` for sequential iteration
/// (the whole file, with full cross-block continuity); `Window`
/// for region iteration (snapshot-init the active set at each
/// block, no continuity check, clamp records to `[start, end]`
/// on `chrom_id`).
#[derive(Debug, Clone, Copy)]
enum RangeClamp {
    None,
    Window {
        chrom_id: u32,
        /// 1-based inclusive lower bound.
        start: u32,
        /// 1-based inclusive upper bound.
        end: u32,
    },
}

/// Iterator over `PspReader`'s records. Owns a mutable borrow of
/// the source, so only one iterator lives at a time per
/// `PspReader`. **Not `Send`** by design — see
/// `ia/feature_implementation_plans/psp_reader.md` §"Tradeoffs"
/// for the rationale.
pub struct RecordsIter<'r, R: Read + Seek> {
    reader: &'r mut PspReader<R>,
    /// Index of the next block to load. Equals
    /// `reader.index.len()` once iteration is exhausted.
    cur_block_idx: usize,
    /// Decoded payload of the current block. `None` before the
    /// first block has been pulled from disk; dropped when the
    /// block is exhausted to free its column buffers.
    cur_block: Option<DecodedBlock>,
    /// Index of the next record (within the current block) to
    /// materialise. Resets to 0 on each new block.
    next_record_in_block: u32,
    /// Cumulative allele offset into the per-allele columns,
    /// scoped to the current block. Resets to 0 on each new
    /// block.
    next_allele_in_block: u32,
    /// Running 1-based reference position of the last-materialised
    /// record. Reset to `cur_block.first_pos` on the first record
    /// of every block.
    last_pos: u32,
    /// Running phase-chain active-slot set, maintained across
    /// records and across blocks for sequential iteration.
    /// Sequential mode cross-checks every block-start snapshot
    /// against this set (spec check #11); region mode resets it
    /// from the snapshot at every block (Slice 5).
    active_chain_slots: BTreeSet<SlotId>,
    /// `RangeClamp::None` for sequential iteration; `Window` for
    /// region iteration. Tells `load_next_block` whether to run
    /// the cross-block continuity check (sequential) or to
    /// snapshot-init from the block header (region).
    clamp: RangeClamp,
    /// Sticky: once `next()` has yielded `Some(Err(_))`, future
    /// calls return `None`. Required for `Iterator` correctness
    /// when the source state after an error is unknown.
    poisoned: bool,
}

impl<'r, R: Read + Seek> RecordsIter<'r, R> {
    fn new(reader: &'r mut PspReader<R>, clamp: RangeClamp) -> Self {
        Self {
            reader,
            cur_block_idx: 0,
            cur_block: None,
            next_record_in_block: 0,
            next_allele_in_block: 0,
            last_pos: 0,
            active_chain_slots: BTreeSet::new(),
            clamp,
            poisoned: false,
        }
    }

    /// Seek to the next un-decoded block and decode it into
    /// `self.cur_block`. Returns `Ok(false)` when no more blocks
    /// remain.
    ///
    /// In sequential mode (the only mode in Slices 2–4) runs spec
    /// check #11: the new block's `active_chain_slots_at_block_start`
    /// snapshot must equal the active set carried in from the
    /// previous block. On a chromosome change the active set must
    /// have been emptied by the writer's chain-closure invariant,
    /// so both sides drop to empty.
    fn load_next_block(&mut self) -> Result<bool, PspReadError> {
        if self.cur_block_idx >= self.reader.index.len() {
            return Ok(false);
        }
        let entry = self.reader.index[self.cur_block_idx];
        self.reader
            .source
            .seek(SeekFrom::Start(entry.block_offset))
            .map_err(io_err("block start seek"))?;
        let (block_header, _consumed) = read_block_header(&mut self.reader.source)?;

        match self.clamp {
            RangeClamp::None => {
                // Spec check #11: snapshot must match the
                // carried-in set (sequential mode only).
                let snapshot: BTreeSet<SlotId> =
                    block_header.active_chain_slots.iter().copied().collect();
                if snapshot != self.active_chain_slots {
                    return Err(PspReadError::BlockHeaderInvariant {
                        kind: BlockHeaderInvariantKind::SnapshotMismatch {
                            block: self.cur_block_idx,
                            snapshot_len: snapshot.len(),
                            expected_len: self.active_chain_slots.len(),
                        },
                    });
                }
            }
            RangeClamp::Window { .. } => {
                // Random-access mode trusts the snapshot — the
                // spec defers check #11 to sequential reads. The
                // running set is reinitialised at every block.
                self.active_chain_slots.clear();
                self.active_chain_slots
                    .extend(block_header.active_chain_slots.iter().copied());
            }
        }

        let decoded = decode_block_payload(&mut self.reader.source, &block_header)?;
        self.cur_block = Some(decoded);
        self.next_record_in_block = 0;
        self.next_allele_in_block = 0;
        // `delta_pos[0] = 0` by writer invariant, so the first
        // record's `pos = first_pos + delta_pos[0] = first_pos`.
        // Initialise `last_pos` to make that fall out of the same
        // addition path used for the rest of the block.
        self.last_pos = block_header.first_pos;
        Ok(true)
    }

    /// Construct the next `PileupRecord` from the currently-loaded
    /// block. Caller must ensure a block is loaded and has records
    /// left. Also runs spec check #10: phase-chain marker
    /// consistency against the running active-slot set.
    fn materialise_next_record(&mut self) -> Result<PileupRecord, PspReadError> {
        let block = self
            .cur_block
            .as_ref()
            .expect("materialise_next_record requires a loaded block");
        let i = self.next_record_in_block as usize;

        // Spec check #10, part 1: every expired slot must currently
        // be in the active set; no new slot may already be active.
        for &slot in &block.expired_chain_slots[i] {
            if !self.active_chain_slots.contains(&slot) {
                return Err(PspReadError::PhaseChainConsistency {
                    block: self.cur_block_idx,
                    record_in_block: self.next_record_in_block,
                    kind: PhaseChainConsistencyKind::ExpiredNotActive { slot },
                });
            }
        }
        for &slot in &block.new_chain_slots[i] {
            if self.active_chain_slots.contains(&slot) {
                return Err(PspReadError::PhaseChainConsistency {
                    block: self.cur_block_idx,
                    record_in_block: self.next_record_in_block,
                    kind: PhaseChainConsistencyKind::NewAlreadyActive { slot },
                });
            }
        }
        // Apply markers (expired then new) — matches the writer's
        // semantics in `apply_record_to_block`.
        for &slot in &block.expired_chain_slots[i] {
            self.active_chain_slots.remove(&slot);
        }
        for &slot in &block.new_chain_slots[i] {
            self.active_chain_slots.insert(slot);
        }

        let delta = block.delta_pos[i];
        let pos = if i == 0 {
            block.first_pos
        } else {
            self.last_pos.saturating_add(delta as u32)
        };

        let n_alleles_here = block.n_alleles[i] as usize;
        let allele_start = self.next_allele_in_block as usize;
        let allele_end = allele_start + n_alleles_here;

        let mut alleles = Vec::with_capacity(n_alleles_here);
        for j in allele_start..allele_end {
            // Spec check #10, part 2: every slot referenced by an
            // allele must be in the post-application active set.
            for &slot in &block.allele_chain_slots[j] {
                if !self.active_chain_slots.contains(&slot) {
                    return Err(PspReadError::PhaseChainConsistency {
                        block: self.cur_block_idx,
                        record_in_block: self.next_record_in_block,
                        kind: PhaseChainConsistencyKind::AlleleReferencesUnknownSlot {
                            allele_index: j - allele_start,
                            slot,
                        },
                    });
                }
            }
            alleles.push(AlleleObservation {
                seq: block.allele_seqs[j].clone(),
                support: AlleleSupportStats {
                    num_obs: block.allele_obs_count[j],
                    q_sum: block.allele_q_sum_log[j],
                    fwd: block.allele_fwd_count[j],
                    placed_left: block.allele_placed_left_count[j],
                    placed_start: block.allele_placed_start_count[j],
                },
                chain_slots: block.allele_chain_slots[j].clone(),
            });
        }

        let record = PileupRecord {
            chrom_id: block.chrom_id,
            pos,
            new_chains: block.new_chain_slots[i].clone(),
            expired_chains: block.expired_chain_slots[i].clone(),
            alleles,
        };

        self.last_pos = pos;
        self.next_record_in_block += 1;
        self.next_allele_in_block += n_alleles_here as u32;
        Ok(record)
    }
}

impl<R: Read + Seek> Iterator for RecordsIter<'_, R> {
    type Item = Result<PileupRecord, PspReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.poisoned {
            return None;
        }
        loop {
            if let Some(block) = &self.cur_block
                && self.next_record_in_block < block.n_records
            {
                let record = match self.materialise_next_record() {
                    Ok(r) => r,
                    Err(e) => {
                        self.poisoned = true;
                        return Some(Err(e));
                    }
                };
                // Region clamp: skip records before `start`,
                // terminate the moment we step past `end`.
                if let RangeClamp::Window {
                    chrom_id,
                    start,
                    end,
                } = self.clamp
                {
                    if record.chrom_id != chrom_id || record.pos > end {
                        return None;
                    }
                    if record.pos < start {
                        continue; // pre-window record, drop and ask for the next
                    }
                }
                return Some(Ok(record));
            }
            // Current block exhausted (or never loaded). Drop its
            // buffers, advance, and load the next.
            if self.cur_block.is_some() {
                self.cur_block = None;
                self.cur_block_idx += 1;
            }
            // For region iteration, peek at the next block's
            // range before paying the decode cost: if its
            // `first_pos > end` (or `chrom_id` no longer matches),
            // we're done.
            if let RangeClamp::Window {
                chrom_id,
                start: _,
                end,
            } = self.clamp
                && self.cur_block_idx < self.reader.index.len()
            {
                let entry = &self.reader.index[self.cur_block_idx];
                if entry.chrom_id != chrom_id || entry.first_pos > end {
                    return None;
                }
            }
            match self.load_next_block() {
                Ok(true) => continue,
                Ok(false) => return None,
                Err(e) => {
                    self.poisoned = true;
                    return Some(Err(e));
                }
            }
        }
    }
}

// ---------------------------------------------------------------------
// Block-header read: grow-on-incomplete loop
// ---------------------------------------------------------------------

/// Read a block header from the source by growing a buffer until
/// `decode_block_header` succeeds. Returns the decoded header and
/// the number of bytes it consumed; positions the source cursor at
/// the first byte after the header.
fn read_block_header<R: Read + Seek>(source: &mut R) -> Result<(BlockHeader, usize), PspReadError> {
    let header_start = source
        .stream_position()
        .map_err(io_err("block header start position"))?;
    let mut buf: Vec<u8> = Vec::with_capacity(BLOCK_HEADER_INITIAL_CHUNK);
    let mut chunk = BLOCK_HEADER_INITIAL_CHUNK;
    loop {
        let prev_len = buf.len();
        buf.resize(prev_len + chunk, 0);
        let read = source
            .read(&mut buf[prev_len..])
            .map_err(io_err("block header"))?;
        if read == 0 {
            return Err(PspReadError::BlockHeaderField {
                field: "header truncated mid-decode",
                source: VarintError::Truncated,
            });
        }
        buf.truncate(prev_len + read);

        match decode_block_header(&buf) {
            Ok((header, consumed)) => {
                // We over-read into `buf`; rewind the source so
                // the caller's cursor lands at the start of the
                // column payloads.
                source
                    .seek(SeekFrom::Start(header_start + consumed as u64))
                    .map_err(io_err("post-block-header seek"))?;
                return Ok((header, consumed));
            }
            Err(PspReadError::BlockHeaderField {
                source: VarintError::Truncated,
                ..
            }) => {
                // Not enough bytes — feed more and retry.
            }
            Err(other) => return Err(other),
        }

        if buf.len() >= BLOCK_HEADER_READ_CAP {
            return Err(PspReadError::BlockHeaderField {
                field: "header exceeds read cap",
                source: VarintError::Overflow,
            });
        }
        chunk = (chunk * 2).min(BLOCK_HEADER_READ_CAP - buf.len());
    }
}

// ---------------------------------------------------------------------
// Block-payload decode: walk the manifest, decode each column
// ---------------------------------------------------------------------

/// Decode every column in a block. The source cursor must be
/// positioned at the first byte of the column payloads (i.e.
/// immediately after the block header).
fn decode_block_payload<R: Read>(
    source: &mut R,
    header: &BlockHeader,
) -> Result<DecodedBlock, PspReadError> {
    let n_records = header.n_records as usize;
    let n_total_alleles = header.n_total_alleles as usize;

    // Per-record slots.
    let mut delta_pos: Option<Vec<u64>> = None;
    let mut n_alleles: Option<Vec<u64>> = None;
    let mut new_chain_slots: Option<Vec<Vec<SlotId>>> = None;
    let mut expired_chain_slots: Option<Vec<Vec<SlotId>>> = None;

    // Per-allele slots.
    let mut allele_seq_len: Option<Vec<u64>> = None;
    let mut allele_seqs: Option<Vec<Vec<u8>>> = None;
    let mut allele_obs_count: Option<Vec<u32>> = None;
    let mut allele_q_sum_log: Option<Vec<f64>> = None;
    let mut allele_fwd_count: Option<Vec<u32>> = None;
    let mut allele_placed_left_count: Option<Vec<u32>> = None;
    let mut allele_placed_start_count: Option<Vec<u32>> = None;
    let mut allele_chain_slots: Option<Vec<Vec<SlotId>>> = None;

    for entry in &header.manifest {
        let column = decode_one_column(
            source,
            entry,
            n_records,
            n_total_alleles,
            allele_seq_len.as_deref(),
        )?;
        match column {
            DecodedColumn::DeltaPos(v) => delta_pos = Some(v),
            DecodedColumn::NAlleles(v) => n_alleles = Some(v),
            DecodedColumn::AlleleSeqLen(v) => allele_seq_len = Some(v),
            DecodedColumn::AlleleSeq(v) => allele_seqs = Some(v),
            DecodedColumn::AlleleObsCount(v) => allele_obs_count = Some(v),
            DecodedColumn::AlleleQSumLog(v) => {
                // Spec check #9: every f64 entry must be finite.
                // The writer enforces this on the produce side
                // ([`InvalidRecordKind::NonFiniteQSum`]); the
                // reader re-checks against torn-frame or
                // hand-crafted attacks.
                for (i, &q) in v.iter().enumerate() {
                    if !q.is_finite() {
                        return Err(PspReadError::NonFiniteFloat {
                            column: "allele-q-sum-log".to_string(),
                            entry: i,
                            value: q,
                        });
                    }
                }
                allele_q_sum_log = Some(v);
            }
            DecodedColumn::AlleleFwdCount(v) => allele_fwd_count = Some(v),
            DecodedColumn::AllelePlacedLeftCount(v) => allele_placed_left_count = Some(v),
            DecodedColumn::AllelePlacedStartCount(v) => allele_placed_start_count = Some(v),
            DecodedColumn::NewChainSlots(v) => new_chain_slots = Some(v),
            DecodedColumn::ExpiredChainSlots(v) => expired_chain_slots = Some(v),
            DecodedColumn::AlleleChainSlots(v) => allele_chain_slots = Some(v),
            DecodedColumn::Unknown => {} // optional minor-bump column: skipped
        }
    }

    // The header parse's `cross_check_against_registry` guarantees
    // every v1.0 required column is present in the manifest, so
    // each `expect` below is structurally unreachable on a
    // schema-valid file.
    Ok(DecodedBlock {
        chrom_id: header.chrom_id,
        first_pos: header.first_pos,
        n_records: header.n_records,
        n_total_alleles: header.n_total_alleles,
        active_chain_slots_snapshot: header.active_chain_slots.clone(),
        delta_pos: delta_pos.expect("delta-pos column required by v1.0 schema"),
        n_alleles: n_alleles.expect("n-alleles column required by v1.0 schema"),
        new_chain_slots: new_chain_slots.expect("new-chain-slots column required by v1.0 schema"),
        expired_chain_slots: expired_chain_slots
            .expect("expired-chain-slots column required by v1.0 schema"),
        allele_seqs: allele_seqs.expect("allele-seq column required by v1.0 schema"),
        allele_obs_count: allele_obs_count
            .expect("allele-obs-count column required by v1.0 schema"),
        allele_q_sum_log: allele_q_sum_log
            .expect("allele-q-sum-log column required by v1.0 schema"),
        allele_fwd_count: allele_fwd_count
            .expect("allele-fwd-count column required by v1.0 schema"),
        allele_placed_left_count: allele_placed_left_count
            .expect("allele-placed-left-count column required by v1.0 schema"),
        allele_placed_start_count: allele_placed_start_count
            .expect("allele-placed-start-count column required by v1.0 schema"),
        allele_chain_slots: allele_chain_slots
            .expect("allele-chain-slots column required by v1.0 schema"),
    })
}

/// One decoded v1.0 column — keyed by [`ColumnKey`] so the
/// per-column dispatch is compile-time exhaustive against the
/// registry.
enum DecodedColumn {
    DeltaPos(Vec<u64>),
    NAlleles(Vec<u64>),
    AlleleSeqLen(Vec<u64>),
    AlleleSeq(Vec<Vec<u8>>),
    AlleleObsCount(Vec<u32>),
    AlleleQSumLog(Vec<f64>),
    AlleleFwdCount(Vec<u32>),
    AllelePlacedLeftCount(Vec<u32>),
    AllelePlacedStartCount(Vec<u32>),
    NewChainSlots(Vec<Vec<SlotId>>),
    ExpiredChainSlots(Vec<Vec<SlotId>>),
    AlleleChainSlots(Vec<Vec<SlotId>>),
    /// Optional column the registry does not know about — the
    /// "minor-bump-additions-are-optional" rule. Slice 2 consumes
    /// the bytes without decoding so the source cursor stays in
    /// sync; the column contents are dropped.
    Unknown,
}

/// Read + decompress + decode one column. `allele_seq_len` is
/// `Some(...)` once the `allele-seq-len` column has been decoded
/// (manifest order is tag-ascending and 0x03 < 0x04); used to
/// chunk the `allele-seq` bytes column.
fn decode_one_column<R: Read>(
    source: &mut R,
    entry: &ColumnManifestEntry,
    n_records: usize,
    n_total_alleles: usize,
    allele_seq_len: Option<&[u64]>,
) -> Result<DecodedColumn, PspReadError> {
    // Unknown tag = optional future column; read past it so the
    // source cursor stays aligned.
    let def = match lookup_by_tag(entry.tag) {
        Some(d) => d,
        None => {
            let mut sink = vec![0u8; entry.compressed_len as usize];
            source
                .read_exact(&mut sink)
                .map_err(io_err("unknown optional column payload"))?;
            return Ok(DecodedColumn::Unknown);
        }
    };
    let column_name = def.name;

    // Spec check #7 case 1: predict + verify before decompression
    // for fixed-width scalar columns.
    if let Some(predicted) = predict_fixed_uncompressed_len(def, n_records, n_total_alleles)
        && predicted as u32 != entry.uncompressed_len
    {
        return Err(PspReadError::UncompressedLenSchemaMismatch {
            column: column_name.to_string(),
            got: entry.uncompressed_len as u64,
            expected: predicted as u64,
        });
    }

    // Spec check #7 case 2 for `allele-seq`: predict against the
    // sum of `allele-seq-len` already decoded in this block.
    if def.key == ColumnKey::AlleleSeq
        && let Some(lens) = allele_seq_len
    {
        let expected: u64 = lens.iter().sum();
        if expected != entry.uncompressed_len as u64 {
            return Err(PspReadError::UncompressedLenSchemaMismatch {
                column: column_name.to_string(),
                got: entry.uncompressed_len as u64,
                expected,
            });
        }
    }

    // Read + decompress + verify decompressed length.
    let mut compressed = vec![0u8; entry.compressed_len as usize];
    source
        .read_exact(&mut compressed)
        .map_err(io_err("block column payload"))?;
    let bytes = zstd_decompress(&compressed).map_err(|source| PspReadError::Zstd {
        context: format!("column {column_name:?}"),
        source,
    })?;
    if bytes.len() as u64 != entry.uncompressed_len as u64 {
        return Err(PspReadError::UncompressedLenMismatch {
            column: column_name.to_string(),
            got: bytes.len() as u64,
            expected: entry.uncompressed_len as u64,
        });
    }

    // Per-column decode, exhaustive on `ColumnKey` so a registry
    // addition forces a new arm here as a compile error.
    let decoded = match def.key {
        ColumnKey::DeltaPos => {
            DecodedColumn::DeltaPos(decode_varint_column(&bytes, n_records, column_name)?)
        }
        ColumnKey::NAlleles => {
            DecodedColumn::NAlleles(decode_varint_column(&bytes, n_records, column_name)?)
        }
        ColumnKey::AlleleSeqLen => {
            let lens = decode_varint_column(&bytes, n_total_alleles, column_name)?;
            // Enforce per-entry cap symmetrically with the writer.
            for (i, &len) in lens.iter().enumerate() {
                if len > MAX_ALLELE_SEQ_LEN {
                    return Err(PspReadError::ColumnElementDecode {
                        column: column_name.to_string(),
                        entry: i,
                        source: super::errors::ScalarDecodeError::Truncated,
                    });
                }
            }
            DecodedColumn::AlleleSeqLen(lens)
        }
        ColumnKey::AlleleSeq => {
            let lens = allele_seq_len.expect(
                "allele-seq-len decoded before allele-seq per manifest tag-ascending order",
            );
            DecodedColumn::AlleleSeq(decode_bytes_split(
                &bytes,
                lens,
                column_name,
                Some(MAX_ALLELE_SEQ_LEN),
            )?)
        }
        ColumnKey::AlleleObsCount => DecodedColumn::AlleleObsCount(decode_scalar_column::<u32>(
            &bytes,
            n_total_alleles,
            column_name,
        )?),
        ColumnKey::AlleleQSumLog => DecodedColumn::AlleleQSumLog(decode_scalar_column::<f64>(
            &bytes,
            n_total_alleles,
            column_name,
        )?),
        ColumnKey::AlleleFwdCount => DecodedColumn::AlleleFwdCount(decode_scalar_column::<u32>(
            &bytes,
            n_total_alleles,
            column_name,
        )?),
        ColumnKey::AllelePlacedLeftCount => DecodedColumn::AllelePlacedLeftCount(
            decode_scalar_column::<u32>(&bytes, n_total_alleles, column_name)?,
        ),
        ColumnKey::AllelePlacedStartCount => DecodedColumn::AllelePlacedStartCount(
            decode_scalar_column::<u32>(&bytes, n_total_alleles, column_name)?,
        ),
        ColumnKey::NewChainSlots => DecodedColumn::NewChainSlots(decode_list_column::<SlotId>(
            &bytes,
            n_records,
            column_name,
        )?),
        ColumnKey::ExpiredChainSlots => {
            DecodedColumn::ExpiredChainSlots(decode_list_column::<SlotId>(
                &bytes,
                n_records,
                column_name,
            )?)
        }
        ColumnKey::AlleleChainSlots => DecodedColumn::AlleleChainSlots(
            decode_list_column::<SlotId>(&bytes, n_total_alleles, column_name)?,
        ),
    };
    Ok(decoded)
}

/// Predict the column's uncompressed byte length from the schema
/// alone (cardinality + fixed-width element type). `None` for
/// columns whose uncompressed length is genuinely variable
/// (varint scalars, list columns, bytes columns). Used by Spec
/// check #7 case 1.
fn predict_fixed_uncompressed_len(
    def: &super::registry::ColumnDef,
    n_records: usize,
    n_total_alleles: usize,
) -> Option<usize> {
    use super::registry::{Cardinality, ColumnPayload};
    match def.payload {
        ColumnPayload::Scalar { element_type } => {
            let width = element_type.fixed_byte_width()?;
            let count = match def.cardinality {
                Cardinality::PerRecord => n_records,
                Cardinality::PerAllele => n_total_alleles,
            };
            Some(count * width)
        }
        ColumnPayload::List { .. } | ColumnPayload::Bytes { .. } => None,
    }
}

/// Wrap an `io::Error` into `PspReadError::Io` with a static
/// context label, so the open-sequence code is not littered with
/// repeated closure literals.
fn io_err(context: &'static str) -> impl Fn(std::io::Error) -> PspReadError {
    move |source| PspReadError::Io { context, source }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::per_sample_caller::pileup::{AlleleObservation, AlleleSupportStats, PileupRecord};
    use crate::per_sample_caller::psp::header::{
        ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance,
    };
    use crate::per_sample_caller::psp::test_fixtures::writer_header;
    use crate::per_sample_caller::psp::writer::PspWriter;
    use std::collections::BTreeMap;
    use std::io::Cursor;

    fn finish_empty_writer(header: WriterHeader) -> Vec<u8> {
        let writer = PspWriter::new(Cursor::new(Vec::new()), header).unwrap();
        writer.finish().unwrap().into_inner()
    }

    fn one_allele_record(chrom_id: u32, pos: u32, seq: &[u8]) -> PileupRecord {
        PileupRecord {
            chrom_id,
            pos,
            new_chains: Vec::new(),
            expired_chains: Vec::new(),
            alleles: vec![AlleleObservation {
                seq: seq.to_vec(),
                support: AlleleSupportStats {
                    num_obs: 1,
                    q_sum: -1.0,
                    fwd: 1,
                    placed_left: 0,
                    placed_start: 1,
                },
                chain_slots: Vec::new(),
            }],
        }
    }

    /// (R1) Zero-record file round-trips through `PspReader::new`
    /// successfully; `header()` matches; `block_index().is_empty()`;
    /// `records()` and `region_records()` yield `None` immediately.
    #[test]
    fn r1_empty_file_round_trip() {
        let bytes = finish_empty_writer(writer_header(1));
        let mut reader = PspReader::new(Cursor::new(bytes)).expect("empty file opens");
        assert_eq!(reader.header().sample, "sample");
        assert!(reader.block_index().is_empty());
        assert_eq!(reader.trailer().n_blocks, 0);
        assert!(reader.records().next().is_none());
        assert!(reader.region_records(0, 1, 100_000).next().is_none());
    }

    /// (R2) The header is accessible before any iteration: opening
    /// plus `header()` is a complete operation that requires no
    /// block reads.
    #[test]
    fn r2_header_accessible_before_iteration() {
        let bytes = finish_empty_writer(writer_header(2));
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.header().chromosomes.len(), 2);
        assert_eq!(reader.header().sample, "sample");
        // No iteration; drop.
    }

    /// (R3) Header round-trip equality: every non-default writer
    /// header field survives the write → read pipeline byte-stable.
    #[test]
    fn r3_header_round_trip_equality() {
        let mut params = BTreeMap::new();
        params.insert("min-mapq".to_string(), ParameterValue::Integer(42));
        params.insert("min-bq".to_string(), ParameterValue::Integer(13));
        params.insert("drop-duplicate".to_string(), ParameterValue::Boolean(true));
        params.insert(
            "tag".to_string(),
            ParameterValue::String("calibration".to_string()),
        );
        let written = WriterHeader {
            format_version: (1, 0),
            sample: "S".repeat(64),
            reference: "GRCh38.fa".to_string(),
            created: "2026-05-13T11:22:33Z".parse().unwrap(),
            chromosomes: vec![
                ChromosomeEntry {
                    name: "chr1".to_string(),
                    length: 248_956_422,
                    md5: "6aef897c3d6ff0c78aff06ac189178dd".to_string(),
                },
                ChromosomeEntry {
                    name: "chr2".to_string(),
                    length: 242_193_529,
                    md5: "f98db672eb0993dcfdabafe2a882905c".to_string(),
                },
            ],
            writer: WriterProvenance {
                tool: "join_per_sample_vcfs".to_string(),
                version: "0.3.0".to_string(),
                subcommand: "per-sample".to_string(),
                input_crams: vec![
                    "a.cram".to_string(),
                    "b.cram".to_string(),
                    "c.cram".to_string(),
                ],
                input_fasta: "GRCh38.fa".to_string(),
                parameters: params,
            },
        };
        let bytes = finish_empty_writer(written.clone());
        let reader = PspReader::new(Cursor::new(bytes)).expect("round-trip opens");
        let parsed = reader.header();
        assert_eq!(parsed.format_version, written.format_version);
        assert_eq!(parsed.sample, written.sample);
        assert_eq!(parsed.reference, written.reference);
        assert_eq!(parsed.chromosomes.len(), written.chromosomes.len());
        for (got, want) in parsed.chromosomes.iter().zip(written.chromosomes.iter()) {
            assert_eq!(got.name, want.name);
            assert_eq!(got.length, want.length);
            assert_eq!(got.md5, want.md5);
        }
        assert_eq!(parsed.writer.tool, written.writer.tool);
        assert_eq!(parsed.writer.version, written.writer.version);
        assert_eq!(parsed.writer.subcommand, written.writer.subcommand);
        assert_eq!(parsed.writer.input_crams, written.writer.input_crams);
        assert_eq!(parsed.writer.input_fasta, written.writer.input_fasta);
        assert_eq!(parsed.writer.parameters, written.writer.parameters);
    }

    /// (R4) `block_index()` matches the writer's emission for a
    /// multi-block fixture. Uses `PspWriter::new_with_block_target`
    /// to force several flushes with a small record count.
    #[test]
    fn r4_block_index_matches_writer_emission() {
        let header = writer_header(1);
        let mut writer =
            PspWriter::new_with_block_target(Cursor::new(Vec::new()), header, 4 * 1024).unwrap();
        for i in 1u32..=600 {
            writer.write_record(&one_allele_record(0, i, b"A")).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let entries = reader.block_index();
        assert!(
            entries.len() >= 2,
            "fixture should produce ≥2 blocks; got {}",
            entries.len()
        );
        for entry in entries {
            assert_eq!(entry.chrom_id, 0);
            assert!(entry.first_pos >= 1);
            assert!(entry.last_pos >= entry.first_pos);
            assert!(entry.last_pos <= 1_000_000);
            assert!(entry.n_records >= 1);
        }
        for w in entries.windows(2) {
            assert!(w[0].block_offset < w[1].block_offset);
            assert!(w[0].last_pos < w[1].first_pos);
        }
    }

    fn allele(
        seq: &[u8],
        num_obs: u32,
        q_sum: f64,
        fwd: u32,
        placed_left: u32,
        placed_start: u32,
        chain_slots: &[SlotId],
    ) -> AlleleObservation {
        AlleleObservation {
            seq: seq.to_vec(),
            support: AlleleSupportStats {
                num_obs,
                q_sum,
                fwd,
                placed_left,
                placed_start,
            },
            chain_slots: chain_slots.to_vec(),
        }
    }

    /// (B1) Single record, single allele, single block — every
    /// scalar survives the round-trip exactly.
    #[test]
    fn b1_single_record_single_allele_round_trip() {
        let want = PileupRecord {
            chrom_id: 0,
            pos: 100,
            new_chains: Vec::new(),
            expired_chains: Vec::new(),
            alleles: vec![allele(b"A", 17, -42.5, 9, 3, 7, &[])],
        };
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer.write_record(&want).unwrap();
        let bytes = writer.finish().unwrap().into_inner();

        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let got: Vec<PileupRecord> = reader
            .records()
            .collect::<Result<_, _>>()
            .expect("records iterate cleanly");
        assert_eq!(got.len(), 1);
        let got = &got[0];
        assert_eq!(got.chrom_id, want.chrom_id);
        assert_eq!(got.pos, want.pos);
        assert_eq!(got.new_chains, want.new_chains);
        assert_eq!(got.expired_chains, want.expired_chains);
        assert_eq!(got.alleles.len(), want.alleles.len());
        let g = &got.alleles[0];
        let w = &want.alleles[0];
        assert_eq!(g.seq, w.seq);
        assert_eq!(g.support.num_obs, w.support.num_obs);
        assert_eq!(g.support.q_sum, w.support.q_sum);
        assert_eq!(g.support.fwd, w.support.fwd);
        assert_eq!(g.support.placed_left, w.support.placed_left);
        assert_eq!(g.support.placed_start, w.support.placed_start);
        assert_eq!(g.chain_slots, w.chain_slots);
    }

    /// (B2) Multi-allele records: a SNP at pos=100 and a
    /// deletion + insertion at pos=200. Allele sequences and
    /// scalars survive.
    #[test]
    fn b2_multi_allele_round_trip() {
        let snp = PileupRecord {
            chrom_id: 0,
            pos: 100,
            new_chains: Vec::new(),
            expired_chains: Vec::new(),
            alleles: vec![
                allele(b"A", 30, -50.0, 16, 0, 30, &[]),
                allele(b"C", 8, -8.0, 4, 0, 8, &[]),
            ],
        };
        let multi = PileupRecord {
            chrom_id: 0,
            pos: 200,
            new_chains: Vec::new(),
            expired_chains: Vec::new(),
            alleles: vec![
                allele(b"ACGT", 25, -40.0, 13, 0, 25, &[]),
                allele(b"A", 5, -5.0, 2, 0, 5, &[]), // deletion of CGT
                allele(b"ACGTAG", 3, -3.0, 1, 0, 3, &[]), // insertion of AG
            ],
        };
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer.write_record(&snp).unwrap();
        writer.write_record(&multi).unwrap();
        let bytes = writer.finish().unwrap().into_inner();

        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let got: Vec<PileupRecord> = reader.records().collect::<Result<_, _>>().unwrap();
        assert_eq!(got.len(), 2);
        // SNP
        assert_eq!(got[0].pos, 100);
        assert_eq!(got[0].alleles.len(), 2);
        assert_eq!(got[0].alleles[0].seq, b"A");
        assert_eq!(got[0].alleles[1].seq, b"C");
        // Multi
        assert_eq!(got[1].pos, 200);
        assert_eq!(got[1].alleles.len(), 3);
        assert_eq!(got[1].alleles[0].seq, b"ACGT");
        assert_eq!(got[1].alleles[1].seq, b"A");
        assert_eq!(got[1].alleles[2].seq, b"ACGTAG");
        assert_eq!(got[1].alleles[0].support.num_obs, 25);
        assert_eq!(got[1].alleles[2].support.num_obs, 3);
    }

    /// (B3) Multi-block file: enough records to force several
    /// flushes at a small block target. Record count and the last
    /// record's coordinates round-trip exactly.
    #[test]
    fn b3_multi_block_round_trip() {
        let header = writer_header(1);
        // 4 KiB target → ~130 records per block at ~30 bytes each.
        let mut writer =
            PspWriter::new_with_block_target(Cursor::new(Vec::new()), header, 4 * 1024).unwrap();
        let n = 1000u32;
        for i in 1..=n {
            writer.write_record(&one_allele_record(0, i, b"A")).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();

        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        assert!(
            reader.block_index().len() >= 3,
            "fixture should produce ≥3 blocks; got {}",
            reader.block_index().len()
        );

        let got: Vec<PileupRecord> = reader.records().collect::<Result<_, _>>().unwrap();
        assert_eq!(got.len() as u32, n);
        // Positions are strictly increasing 1..=n.
        for (i, r) in got.iter().enumerate() {
            assert_eq!(r.pos as usize, i + 1);
            assert_eq!(r.chrom_id, 0);
            assert_eq!(r.alleles.len(), 1);
            assert_eq!(r.alleles[0].seq, b"A");
        }
    }

    /// (B4) Multi-chromosome: 10 records on chr1 then 10 on chr2.
    /// The chrom_id-change flush boundary produces two blocks;
    /// per-chrom record counts round-trip.
    #[test]
    fn b4_multi_chromosome_round_trip() {
        let header = writer_header(2);
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), header).unwrap();
        for i in 1u32..=10 {
            writer.write_record(&one_allele_record(0, i, b"A")).unwrap();
        }
        for i in 1u32..=10 {
            writer.write_record(&one_allele_record(1, i, b"C")).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();

        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        assert_eq!(
            reader.block_index().len(),
            2,
            "chrom change forces a flush; expected 2 blocks"
        );
        let got: Vec<PileupRecord> = reader.records().collect::<Result<_, _>>().unwrap();
        assert_eq!(got.len(), 20);
        let n_chr1 = got.iter().filter(|r| r.chrom_id == 0).count();
        let n_chr2 = got.iter().filter(|r| r.chrom_id == 1).count();
        assert_eq!(n_chr1, 10);
        assert_eq!(n_chr2, 10);
        // Chrom 1 records carry seq "A", chrom 2 records "C".
        for r in got.iter().filter(|r| r.chrom_id == 0) {
            assert_eq!(r.alleles[0].seq, b"A");
        }
        for r in got.iter().filter(|r| r.chrom_id == 1) {
            assert_eq!(r.alleles[0].seq, b"C");
        }
    }

    /// Build a record that opens slot `slot_id` at this position.
    fn record_open_slot(chrom_id: u32, pos: u32, slot_id: SlotId) -> PileupRecord {
        PileupRecord {
            chrom_id,
            pos,
            new_chains: vec![slot_id],
            expired_chains: Vec::new(),
            alleles: vec![allele(b"A", 1, -1.0, 1, 0, 1, &[slot_id])],
        }
    }

    /// Record that consumes slot `slot_id` (references it but does
    /// not close it).
    fn record_use_slot(chrom_id: u32, pos: u32, slot_id: SlotId) -> PileupRecord {
        PileupRecord {
            chrom_id,
            pos,
            new_chains: Vec::new(),
            expired_chains: Vec::new(),
            alleles: vec![allele(b"A", 1, -1.0, 1, 0, 1, &[slot_id])],
        }
    }

    /// Record that closes slot `slot_id`.
    fn record_close_slot(chrom_id: u32, pos: u32, slot_id: SlotId) -> PileupRecord {
        PileupRecord {
            chrom_id,
            pos,
            new_chains: Vec::new(),
            expired_chains: vec![slot_id],
            alleles: vec![allele(b"A", 1, -1.0, 1, 0, 1, &[])],
        }
    }

    /// (B5 / R9) Slot 7 opens in block N's last record, gets
    /// referenced (consumed) in block N+1's first record, and
    /// closes later. The reader's cross-block snapshot check
    /// passes, and slot 7 is present in the per-record
    /// `chain_slots` on both sides of the boundary.
    #[test]
    fn b5_r9_phase_chain_across_block_boundary() {
        let header = writer_header(1);
        // Tiny block target so a few records fill a block.
        let mut writer =
            PspWriter::new_with_block_target(Cursor::new(Vec::new()), header, 1024).unwrap();

        // Block 0: 100 plain records, then one that opens slot 7.
        for i in 1u32..=100 {
            writer.write_record(&one_allele_record(0, i, b"A")).unwrap();
        }
        writer.write_record(&record_open_slot(0, 101, 7)).unwrap();
        // Block 1: a record that uses slot 7, then one that closes it.
        writer.write_record(&record_use_slot(0, 102, 7)).unwrap();
        writer.write_record(&record_close_slot(0, 103, 7)).unwrap();
        let bytes = writer.finish().unwrap().into_inner();

        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        assert!(
            reader.block_index().len() >= 2,
            "fixture should produce ≥2 blocks; got {}",
            reader.block_index().len()
        );
        let got: Vec<PileupRecord> = reader
            .records()
            .collect::<Result<_, _>>()
            .expect("phase-chain continuity holds across block boundaries");

        // Slot 7 appears on the open record, the use record, and
        // the close record's `expired_chains`.
        let opener = got.iter().find(|r| r.pos == 101).unwrap();
        let user = got.iter().find(|r| r.pos == 102).unwrap();
        let closer = got.iter().find(|r| r.pos == 103).unwrap();
        assert!(opener.new_chains.contains(&7));
        assert!(opener.alleles[0].chain_slots.contains(&7));
        assert!(user.alleles[0].chain_slots.contains(&7));
        assert!(closer.expired_chains.contains(&7));
    }

    /// (R10) Tight 4-record fixture engineered so block 1's
    /// `active_chain_slots_at_block_start` snapshot is `[7]`. The
    /// test then mutates that snapshot byte (slot 7 → slot 9),
    /// breaking the carried-in active set continuity, and asserts
    /// the reader rejects the file at the snapshot check.
    #[test]
    fn r10_snapshot_mismatch_at_block_boundary() {
        let header = writer_header(1);
        // Target = 64 B: each ~30-byte record fills a block fast,
        // so the second flush boundary lands between the
        // open-slot record (in block 0) and the use-slot record
        // (in block 1). Block 1's snapshot is therefore `[7]`.
        let mut writer =
            PspWriter::new_with_block_target(Cursor::new(Vec::new()), header, 64).unwrap();
        writer.write_record(&one_allele_record(0, 1, b"A")).unwrap();
        writer.write_record(&record_open_slot(0, 2, 7)).unwrap();
        writer.write_record(&record_use_slot(0, 3, 7)).unwrap();
        writer.write_record(&record_close_slot(0, 4, 7)).unwrap();
        let mut bytes = writer.finish().unwrap().into_inner();

        // Block 1's header at this fixture: 5 single-byte varints
        // (chrom_id=0, first_pos=3, n_records=2,
        // n_total_alleles=2, slot_count=1), then the snapshot
        // slot LE u16 starting at `block1_offset + 5`. Slot 7
        // (0x07 0x00) → slot 9 (0x09 0x00) breaks the snapshot.
        let (block1_offset, n_blocks) = {
            let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
            (
                reader.block_index()[1].block_offset as usize,
                reader.block_index().len(),
            )
        };
        assert!(
            n_blocks >= 2,
            "fixture must produce ≥2 blocks; got {n_blocks}"
        );
        let snapshot_slot_offset = block1_offset + 5;
        assert_eq!(
            bytes[snapshot_slot_offset], 0x07,
            "fixture invariant: block 1's snapshot slot LE u16 starts with 0x07"
        );
        assert_eq!(bytes[snapshot_slot_offset + 1], 0x00);
        bytes[snapshot_slot_offset] = 0x09;

        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let err = reader
            .records()
            .find_map(|r| r.err())
            .expect("snapshot mismatch must surface as an error");
        match err {
            PspReadError::BlockHeaderInvariant {
                kind: BlockHeaderInvariantKind::SnapshotMismatch { .. },
            } => {}
            PspReadError::PhaseChainConsistency { .. } => {
                // Acceptable: the mismatch may surface first as a
                // per-record consistency violation if the snapshot
                // happens to be the *right* size but wrong slot id.
            }
            other => panic!("expected SnapshotMismatch or PhaseChainConsistency, got {other:?}"),
        }
    }

    /// Build a 3-block fixture: 3 blocks on chrom 0 each holding
    /// 100 records, covering positions 1–100, 101–200, 201–300.
    fn three_block_chrom0_fixture(target: usize) -> Vec<u8> {
        let mut writer =
            PspWriter::new_with_block_target(Cursor::new(Vec::new()), writer_header(1), target)
                .unwrap();
        for i in 1u32..=300 {
            writer.write_record(&one_allele_record(0, i, b"A")).unwrap();
        }
        writer.finish().unwrap().into_inner()
    }

    /// (B7) Random-access region query across multiple blocks.
    /// 1000 records on two chromosomes; `region_records(0,
    /// 50_000, 60_000)` (out-of-range on a fixture whose max pos
    /// is 1000) yields zero records. Inside-range case:
    /// `region_records(0, 100, 200)` returns positions 100..=200.
    #[test]
    fn b7_random_access_region_query() {
        let mut writer =
            PspWriter::new_with_block_target(Cursor::new(Vec::new()), writer_header(2), 8 * 1024)
                .unwrap();
        for i in 1u32..=500 {
            writer.write_record(&one_allele_record(0, i, b"A")).unwrap();
        }
        for i in 1u32..=500 {
            writer.write_record(&one_allele_record(1, i, b"C")).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();
        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        // In-range: positions 100..=200 on chr0.
        let in_range: Vec<PileupRecord> = reader
            .region_records(0, 100, 200)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(in_range.len(), 101); // inclusive both ends
        for r in &in_range {
            assert_eq!(r.chrom_id, 0);
            assert!((100..=200).contains(&r.pos));
        }
        assert_eq!(in_range.first().unwrap().pos, 100);
        assert_eq!(in_range.last().unwrap().pos, 200);
        // Out-of-range: positions 50_000..=60_000 (no coverage).
        let out_of_range: Vec<PileupRecord> = reader
            .region_records(0, 50_000, 60_000)
            .collect::<Result<_, _>>()
            .unwrap();
        assert!(out_of_range.is_empty());
    }

    /// (R5) Region query that spans multiple blocks: positions
    /// 50..=250 across 3 blocks each covering 100 positions.
    #[test]
    fn r5_region_across_multiple_blocks() {
        // Choose target so each block holds ~100 records.
        // ~30 bytes per record, 100 records = ~3000 bytes.
        let bytes = three_block_chrom0_fixture(3 * 1024);
        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        assert_eq!(
            reader.block_index().len(),
            3,
            "fixture should produce exactly 3 blocks; got {}",
            reader.block_index().len()
        );
        let got: Vec<PileupRecord> = reader
            .region_records(0, 50, 250)
            .collect::<Result<_, _>>()
            .unwrap();
        // Records inclusive on both ends — 50..=250 = 201 records.
        assert_eq!(got.len(), 201);
        for (k, r) in got.iter().enumerate() {
            assert_eq!(r.pos as usize, 50 + k);
            assert_eq!(r.chrom_id, 0);
        }
    }

    /// (R6) Region wholly outside any block's coverage returns no
    /// records and no error.
    #[test]
    fn r6_region_wholly_outside_coverage() {
        let bytes = three_block_chrom0_fixture(3 * 1024);
        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let got: Vec<PileupRecord> = reader
            .region_records(0, 1_000_000, 2_000_000)
            .collect::<Result<_, _>>()
            .unwrap();
        assert!(got.is_empty());
    }

    /// (R7) Dropping a mid-stream iterator and taking a fresh one
    /// restarts from the top — iterators do not share progress
    /// with the source.
    #[test]
    fn r7_iterator_dropped_mid_stream_restarts() {
        let header = writer_header(1);
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), header).unwrap();
        for i in 1u32..=100 {
            writer.write_record(&one_allele_record(0, i, b"A")).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();
        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        // Consume the first 20 records, then drop.
        {
            let mut it = reader.records();
            for _ in 0..20 {
                let _ = it.next().unwrap().unwrap();
            }
        }
        // A fresh iterator yields all 100.
        let all: Vec<_> = reader.records().collect::<Result<_, _>>().unwrap();
        assert_eq!(all.len(), 100);
        assert_eq!(all[0].pos, 1);
        assert_eq!(all[99].pos, 100);
    }

    /// (R8) Sequential and region iterators produce the same
    /// records inside any window.
    #[test]
    fn r8_sequential_and_region_agree_inside_window() {
        let bytes = three_block_chrom0_fixture(3 * 1024);
        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();

        let sequential: Vec<PileupRecord> = reader.records().collect::<Result<_, _>>().unwrap();
        let in_window: Vec<&PileupRecord> = sequential
            .iter()
            .filter(|r| r.chrom_id == 0 && (75..=175).contains(&r.pos))
            .collect();

        let via_region: Vec<PileupRecord> = reader
            .region_records(0, 75, 175)
            .collect::<Result<_, _>>()
            .unwrap();

        assert_eq!(in_window.len(), via_region.len());
        for (s, r) in in_window.iter().zip(via_region.iter()) {
            assert_eq!(s.chrom_id, r.chrom_id);
            assert_eq!(s.pos, r.pos);
            assert_eq!(s.alleles.len(), r.alleles.len());
        }
    }

    // ---------- Slice 6: F1–F8 corruption / negative tests ------

    fn small_multi_block_fixture() -> Vec<u8> {
        let mut writer =
            PspWriter::new_with_block_target(Cursor::new(Vec::new()), writer_header(1), 512)
                .unwrap();
        for i in 1u32..=50 {
            writer.write_record(&one_allele_record(0, i, b"A")).unwrap();
        }
        writer.finish().unwrap().into_inner()
    }

    /// (F1) Head magic flipped — `PspReader::new` returns
    /// `BadHeadMagic` immediately.
    #[test]
    fn f1_head_magic_flipped() {
        let mut bytes = finish_empty_writer(writer_header(1));
        bytes[1] = b'Q'; // PSP\n → PQP\n
        let err = PspReader::new(Cursor::new(bytes)).expect_err("bad magic must fail");
        assert!(matches!(err, PspReadError::BadHeadMagic { .. }));
    }

    /// (F2) Truncate the file mid-block — iteration surfaces an
    /// `Io` (or `Zstd`) error, or open-time validation already
    /// fails because the trailer is missing/garbled.
    #[test]
    fn f2_truncated_mid_block() {
        let bytes = small_multi_block_fixture();
        let truncate_at = {
            let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
            reader.block_index()[0].block_offset as usize + 50
        };
        // Lop everything past `truncate_at`, then patch the
        // trailer back in so the file at least *looks* opened.
        // The trailer still points the index past EOF so open
        // will fail at index validation — that's the truncation
        // surfacing.
        let mut truncated = bytes[..truncate_at].to_vec();
        truncated.extend_from_slice(&bytes[bytes.len() - TRAILER_BYTES..]);
        if let Ok(mut reader) = PspReader::new(Cursor::new(truncated)) {
            let err = reader
                .records()
                .find_map(|r| r.err())
                .expect("truncation must surface during iteration");
            assert!(matches!(
                err,
                PspReadError::Io { .. } | PspReadError::Zstd { .. }
            ));
        }
        // else: open-time variant accepted
    }

    /// (F3) Flip one byte inside a compressed column payload —
    /// the reader surfaces `Zstd` (or one of the downstream
    /// length / decode errors when the corruption survives the
    /// frame check).
    #[test]
    fn f3_torn_zstd_frame() {
        let mut bytes = small_multi_block_fixture();
        let block_offset = {
            let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
            reader.block_index()[0].block_offset as usize
        };
        let flip_at = block_offset + 200;
        assert!(flip_at < bytes.len() - TRAILER_BYTES);
        bytes[flip_at] ^= 0xFF;
        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let err = reader
            .records()
            .find_map(|r| r.err())
            .expect("torn zstd must surface as an error");
        assert!(matches!(
            err,
            PspReadError::Zstd { .. }
                | PspReadError::UncompressedLenMismatch { .. }
                | PspReadError::ColumnTruncated { .. }
                | PspReadError::ColumnTrailingBytes { .. }
                | PspReadError::ColumnElementDecode { .. }
        ));
    }

    /// (F4) Flip one byte inside the block index region —
    /// `PspReader::new` returns `IndexChecksum`.
    #[test]
    fn f4_corrupted_block_index() {
        let mut bytes = small_multi_block_fixture();
        let index_offset = {
            let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
            reader.trailer().index_offset as usize
        };
        let flip_at = index_offset + 4;
        assert!(flip_at < bytes.len() - TRAILER_BYTES);
        bytes[flip_at] ^= 0xFF;
        let err = PspReader::new(Cursor::new(bytes)).expect_err("corrupted index must fail");
        assert!(matches!(err, PspReadError::IndexChecksum { .. }));
    }

    /// (F6) Block manifest with tags out of ascending order is
    /// rejected by `encode_block_header` (the writer side; the
    /// reader path inherits the same validator via
    /// `decode_block_header`).
    #[test]
    fn f6_out_of_order_column_tags() {
        use super::super::block::{BlockHeader, ColumnManifestEntry, encode_block_header};
        let bad = BlockHeader {
            chrom_id: 0,
            first_pos: 1,
            n_records: 1,
            n_total_alleles: 1,
            active_chain_slots: vec![],
            manifest: vec![
                ColumnManifestEntry {
                    tag: 0x02,
                    compressed_len: 0,
                    uncompressed_len: 0,
                },
                ColumnManifestEntry {
                    tag: 0x01,
                    compressed_len: 0,
                    uncompressed_len: 0,
                },
            ],
        };
        let mut buf = Vec::new();
        let err = encode_block_header(&bad, &mut buf).unwrap_err();
        assert!(matches!(
            err,
            BlockHeaderInvariantKind::ManifestTagsNotAscending { .. }
        ));
    }

    /// (F7) Manifest `uncompressed_len` lies about a fixed-width
    /// scalar column's size — the reader fires
    /// `UncompressedLenSchemaMismatch` *before* decompression.
    #[test]
    fn f7_uncompressed_len_schema_mismatch() {
        use super::super::block::{BlockHeader, ColumnManifestEntry, encode_block_header};
        let bad = BlockHeader {
            chrom_id: 0,
            first_pos: 1,
            n_records: 1,
            n_total_alleles: 1,
            active_chain_slots: vec![],
            manifest: vec![ColumnManifestEntry {
                tag: 0x10, // allele-obs-count, u32 per allele
                compressed_len: 0,
                uncompressed_len: 5, // lies; truth is 4 for 1 allele × 4 bytes
            }],
        };
        let mut buf = Vec::new();
        encode_block_header(&bad, &mut buf).unwrap();
        let mut src = Cursor::new(buf);
        let (header, _) = read_block_header(&mut src).unwrap();
        let err = decode_block_payload(&mut src, &header).unwrap_err();
        match err {
            PspReadError::UncompressedLenSchemaMismatch {
                column,
                got,
                expected,
            } => {
                assert_eq!(column, "allele-obs-count");
                assert_eq!(got, 5);
                assert_eq!(expected, 4);
            }
            other => panic!("expected UncompressedLenSchemaMismatch, got {other:?}"),
        }
    }

    /// (F8) Forge block 0's snapshot to declare slot_count = 1
    /// (originally 0). Sequential iteration with an empty
    /// carried-in active set fires `SnapshotMismatch` (or, if
    /// the forged byte stream mangles the rest of the header,
    /// some other block-header variant).
    #[test]
    fn f8_first_block_non_empty_snapshot() {
        let mut bytes = small_multi_block_fixture();
        let block0_offset = {
            let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
            reader.block_index()[0].block_offset as usize
        };
        // For this fixture every leading field encodes to one
        // byte, so slot_count sits at `block0_offset + 4`.
        let slot_count_offset = block0_offset + 4;
        assert_eq!(
            bytes[slot_count_offset], 0,
            "fixture invariant: block 0 starts with empty snapshot"
        );
        bytes[slot_count_offset] = 1;
        let mut reader = match PspReader::new(Cursor::new(bytes)) {
            Ok(r) => r,
            Err(_) => return, // forging mangled the file at open
        };
        let err = reader
            .records()
            .find_map(|r| r.err())
            .expect("forged non-empty first-block snapshot must error");
        assert!(matches!(
            err,
            PspReadError::BlockHeaderInvariant { .. } | PspReadError::BlockHeaderField { .. }
        ));
    }
}
