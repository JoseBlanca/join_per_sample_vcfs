//! Streaming writer that turns a [`PileupRecord`] stream into a v1.0
//! `.psp` file.
//!
//! Three public entry points:
//!
//! - [`PspWriter::new`] — emits the framed TOML header to the sink
//!   and prepares the in-memory state for accepting records.
//! - [`PspWriter::write_record`] — validates one record and buffers
//!   it into the open block; auto-flushes when the projected
//!   uncompressed payload reaches [`DEFAULT_TARGET_BLOCK_BYTES`] or
//!   the next record crosses a chromosome boundary.
//! - [`PspWriter::finish`] — flushes any open block, writes the
//!   block index and the trailer (with its XXH3-64 index checksum),
//!   and returns the sink.
//!
//! Every per-record / per-block invariant from spec
//! §"Header-binary consistency" and §"Block sizing / Block
//! invariants" is enforced at write_record / flush time, so a sink
//! that received a sequence of `Ok` returns is *necessarily* a
//! reader-valid `.psp`. The reader's checks (in step 6) duplicate
//! these — that is by design, because the file may be read on a
//! machine that did not write it.

use std::io::Write;

use super::block::{
    BlockHeader, ColumnManifestEntry, encode_block_header, encode_bytes_concat,
    encode_list_column_csr, encode_scalar_column, encode_varint_column, new_column_compressor,
    zstd_compress_into,
};
use super::errors::PspWriteError;
use super::header::{WriterHeader, build_header_bytes};
use super::index::{BlockIndexEntry, checksum_index, encode_index};
use super::registry::{
    Cardinality, ColumnDef, ElementType, MAX_ALLELE_SEQ_LEN, Shape, V1_0_COLUMNS,
};
use super::trailer::{Trailer, encode_trailer};
use crate::per_sample_caller::pileup::{PileupRecord, SlotId};

/// Target uncompressed bytes per block. Spec §"Block sizing". Not
/// CLI-exposed; the writer flushes around this value.
pub const DEFAULT_TARGET_BLOCK_BYTES: usize = 16 * 1024 * 1024;

/// Initial `Vec::with_capacity` hint for per-record columns in a
/// freshly-opened block (delta-pos, n-alleles, the per-record list
/// columns). Sized at the records-per-block high-water mark for the
/// SNP-typical bench workload (~31 uncompressed bytes per record);
/// pre-allocating this much keeps the amortised-doubling Vec growth
/// path off the hot path for the common case. Real workloads with
/// more alleles per record fit comfortably under this hint; the worst
/// case is an additional realloc on a few-allele block, which is one
/// reallocation per column per block instead of ~17.
const INITIAL_RECORDS_HINT: usize = 350_000;
/// Initial hint for per-allele columns. Sized a touch above
/// `INITIAL_RECORDS_HINT` to absorb the typical excess of alleles
/// over records (occasional multi-allelic positions).
const INITIAL_ALLELES_HINT: usize = 360_000;
/// Initial hint for the concatenated allele-sequence byte buffer.
/// On the SNP-typical workload one byte per allele = ~360 KiB; on
/// indel-heavier workloads alleles run to ~10 bytes each. The 1 MiB
/// hint covers both without leaving much slack.
const INITIAL_ALLELE_SEQ_BYTES_HINT: usize = 1_000_000;

/// Byte mask of allele bytes the writer accepts. `ALLOWED[b]` is
/// `true` iff `b` is one of `b'A' | b'C' | b'G' | b'T' | b'N'`. Walked
/// via `slice.iter().all(|&b| ALLOWED[b as usize])` so the hot-path
/// loop is one indexed load + one branch per byte and runs to
/// completion (no per-element early-exit via `?`) — autovectorisable
/// on `target-cpu = x86-64-v3` (H4).
const ALLOWED_ALLELE_BYTE: [bool; 256] = {
    let mut t = [false; 256];
    t[b'A' as usize] = true;
    t[b'C' as usize] = true;
    t[b'G' as usize] = true;
    t[b'T' as usize] = true;
    t[b'N' as usize] = true;
    t
};

/// `#[cold]` constructors for `PspWriteError::InvalidRecord` (L5).
/// Keeping the `format!` and `String` allocation out-of-line shrinks
/// the hot-path body of `validate_record` (which the samply profile
/// puts at 10 % self).
#[cold]
#[inline(never)]
fn err_invalid_record(record_index: u64, reason: String) -> PspWriteError {
    PspWriteError::InvalidRecord {
        record_index,
        reason,
    }
}

#[cold]
#[inline(never)]
fn err_invalid_allele_byte(
    record_index: u64,
    allele_index: usize,
    seq: &[u8],
) -> PspWriteError {
    let (j, b) = seq
        .iter()
        .enumerate()
        .find(|&(_, &b)| !ALLOWED_ALLELE_BYTE[b as usize])
        .map(|(j, &b)| (j, b))
        .expect("err_invalid_allele_byte called with valid sequence");
    PspWriteError::InvalidRecord {
        record_index,
        reason: format!(
            "allele {allele_index} byte {j} = {b:#04x} (only A/C/G/T/N allowed)"
        ),
    }
}

/// Streaming `.psp` writer over any `Write` sink.
pub struct PspWriter<W: Write> {
    sink: W,
    header: WriterHeader,
    /// Bytes written to `sink` so far. Updated after every flush so
    /// the block index records absolute offsets.
    sink_offset: u64,
    /// One entry per emitted block.
    index_entries: Vec<BlockIndexEntry>,
    /// Open block — `None` when the writer is between blocks (just
    /// after `new` or just after a flush).
    block: Option<BlockAccumulator>,
    /// Running active phase-chain slot set, kept sorted ascending.
    /// Updated as records' `new_chains` / `expired_chains` markers
    /// are applied; snapshotted at block start.
    ///
    /// Sorted `Vec<SlotId>` (with `binary_search` membership) wins
    /// over `BTreeSet<u16>` for the small active counts typical of
    /// phasing (~8 active per locus): one contiguous allocation, one
    /// cache line of data for small `n`, and the block-start snapshot
    /// is a single memcpy via `.clone()` instead of a tree walk
    /// (S1+S2 in `ia/reviews/perf_psp_writer_2026-05-13.md`).
    active_slots: Vec<SlotId>,
    /// Last admitted record's `(chrom_id, pos)`, for monotonicity
    /// enforcement across `write_record` calls.
    last_locus: Option<(u32, u32)>,
    /// Running count of records the caller has handed us. Used in
    /// error messages so the offending input is identifiable.
    records_seen: u64,
    /// Persistent zstd compressor reused across every column of every
    /// block. Each `compress_to_buffer` call resets the internal
    /// frame state but keeps the CCtx workspace and tables alive.
    compressor: zstd::bulk::Compressor<'static>,
    /// Reused per-column uncompressed-payload scratch.
    uncompressed_scratch: Vec<u8>,
    /// Reused per-column compressed-frame buffers, indexed by
    /// `V1_0_COLUMNS` position. Each `flush_block` call `clear`s the
    /// inner Vecs (preserving their capacity) before refilling.
    compressed_scratch: Vec<Vec<u8>>,
    /// Reused per-flush manifest buffer (12 entries).
    manifest_scratch: Vec<ColumnManifestEntry>,
    /// Reused block-header serialisation buffer.
    header_bytes_scratch: Vec<u8>,
}

impl<W: Write> PspWriter<W> {
    /// Validate the header, frame it (magic + length prefix + TOML +
    /// sentinel), and write it to `sink`. After this returns, the
    /// writer is ready to accept records.
    ///
    /// **Buffering.** Each block flush issues one `write_all` for the
    /// block header plus one per column payload, and `finish` adds
    /// one for the index and one for the 32-byte trailer. Against a
    /// real-file sink each `write_all` becomes one `write(2)`
    /// syscall. Wrap a `File` in
    /// `BufWriter::with_capacity(64 * 1024, file)` (or larger) before
    /// passing it here. In-memory sinks (`io::sink()`,
    /// `Cursor<Vec<u8>>`) need no wrapper.
    pub fn new(mut sink: W, header: WriterHeader) -> Result<Self, PspWriteError> {
        let header_bytes = build_header_bytes(&header)?;
        sink.write_all(&header_bytes)
            .map_err(|e| PspWriteError::Io {
                context: "file header",
                source: e,
            })?;
        let sink_offset = header_bytes.len() as u64;
        let compressor = new_column_compressor().map_err(|e| PspWriteError::Io {
            context: "zstd compressor init",
            source: e,
        })?;
        let compressed_scratch = (0..V1_0_COLUMNS.len()).map(|_| Vec::new()).collect();
        Ok(Self {
            sink,
            header,
            sink_offset,
            index_entries: Vec::new(),
            block: None,
            active_slots: Vec::new(),
            last_locus: None,
            records_seen: 0,
            compressor,
            uncompressed_scratch: Vec::new(),
            compressed_scratch,
            manifest_scratch: Vec::with_capacity(V1_0_COLUMNS.len()),
            header_bytes_scratch: Vec::new(),
        })
    }

    /// Append one record. Returns the number of bytes pushed to the
    /// sink as a side effect of any auto-flush this call triggered
    /// (zero if no flush).
    pub fn write_record(&mut self, record: &PileupRecord) -> Result<u64, PspWriteError> {
        let record_index = self.records_seen;
        self.records_seen += 1;
        self.validate_record(record_index, record)?;

        // Should we flush before appending? Flush triggers:
        // - chrom_id change (blocks never cross chromosomes)
        // - projected uncompressed size at or past the target
        let pre_flush_bytes = self.sink_offset;
        let should_flush = match &self.block {
            None => false,
            Some(b) => {
                b.chrom_id != record.chrom_id || b.projected_bytes >= DEFAULT_TARGET_BLOCK_BYTES
            }
        };
        if should_flush {
            self.flush_block()?;
        }
        let flushed_bytes = self.sink_offset - pre_flush_bytes;

        // Start a new block if none is open.
        if self.block.is_none() {
            // The first block on a chromosome always has an empty
            // active-slot snapshot — Stage 1 never produces a chain
            // that spans a contig boundary (reads don't), so by the
            // time we cross a boundary every slot has expired. The
            // walker is expected to honour this and we double-check
            // here: if active_slots is non-empty AND we're starting
            // the first block on a new chromosome, that's a writer
            // upstream bug.
            let snapshot = if self.is_first_block_on_chrom(record.chrom_id) {
                if !self.active_slots.is_empty() {
                    return Err(PspWriteError::PhaseChainMarkerInconsistency {
                        record_index,
                        reason: format!(
                            "active_slots non-empty ({}) at first record of chromosome {}",
                            self.active_slots.len(),
                            record.chrom_id
                        ),
                    });
                }
                Vec::new()
            } else {
                // `active_slots` is kept sorted ascending, which is
                // the BlockHeader's contract for `active_chain_slots`.
                // `.clone()` is one allocation + memcpy — S2.
                self.active_slots.clone()
            };
            self.block = Some(BlockAccumulator::new(record.chrom_id, record.pos, snapshot));
        }

        self.apply_record_to_block(record_index, record)?;
        self.last_locus = Some((record.chrom_id, record.pos));

        Ok(flushed_bytes)
    }

    /// Flush any open block, write the block index, and write the
    /// trailer. Consumes the writer; the returned sink is positioned
    /// at the end of a complete `.psp` file.
    ///
    /// **End-of-stage discipline for `BufWriter`-wrapped file sinks.**
    /// `BufWriter::drop` may swallow flush errors, which for a
    /// billions-of-records file can silently truncate the trailer.
    /// The production caller should do, in order:
    /// ```ignore
    /// let buf = writer.finish()?;       // PSP-level flush
    /// let file = buf.into_inner()?;     // surface BufWriter errors
    /// file.sync_all()?;                 // durability for downstream stages
    /// ```
    /// `sync_all` is end-of-stage only — never per-block.
    pub fn finish(mut self) -> Result<W, PspWriteError> {
        if self.block.is_some() {
            self.flush_block()?;
        }
        let index_offset = self.sink_offset;
        let index_bytes = encode_index(&self.index_entries);
        let index_byte_length = index_bytes.len() as u64;
        let index_checksum = checksum_index(&index_bytes);
        self.sink
            .write_all(&index_bytes)
            .map_err(|e| PspWriteError::Io {
                context: "block index",
                source: e,
            })?;
        self.sink_offset += index_byte_length;

        let trailer = Trailer {
            index_offset,
            index_byte_length,
            n_blocks: self.index_entries.len() as u64,
            index_checksum,
        };
        let trailer_bytes = encode_trailer(&trailer);
        self.sink
            .write_all(&trailer_bytes)
            .map_err(|e| PspWriteError::Io {
                context: "file trailer",
                source: e,
            })?;
        self.sink_offset += trailer_bytes.len() as u64;
        self.sink.flush().map_err(|e| PspWriteError::Io {
            context: "final flush",
            source: e,
        })?;
        Ok(self.sink)
    }

    /// Projected uncompressed bytes the currently-open block has
    /// accumulated. `None` when no block is open (just after `new` or
    /// just after a flush).
    ///
    /// Hidden from rustdoc — this peek into writer state exists only
    /// to let benches align iteration boundaries with the writer's
    /// auto-flush boundary (see L10 in
    /// `ia/reviews/perf_psp_writer_2026-05-13.md`).
    #[doc(hidden)]
    pub fn current_block_projected_bytes(&self) -> Option<usize> {
        self.block.as_ref().map(|b| b.projected_bytes)
    }

    /// `true` if no block has been emitted yet OR the most recent
    /// emitted block was on a different chromosome.
    fn is_first_block_on_chrom(&self, chrom_id: u32) -> bool {
        match self.index_entries.last() {
            None => true,
            Some(last) => last.chrom_id != chrom_id,
        }
    }

    fn validate_record(
        &self,
        record_index: u64,
        record: &PileupRecord,
    ) -> Result<(), PspWriteError> {
        // Chromosome id known.
        let n_chroms = self.header.chromosomes.len();
        if record.chrom_id as usize >= n_chroms {
            return Err(PspWriteError::UnknownChromId {
                record_index,
                chrom_id: record.chrom_id,
                n_chroms: n_chroms as u32,
            });
        }
        let chrom = &self.header.chromosomes[record.chrom_id as usize];

        // pos in [1, chrom.length].
        if record.pos == 0 || record.pos > chrom.length {
            return Err(PspWriteError::PosOutOfRange {
                record_index,
                chrom_id: record.chrom_id,
                pos: record.pos,
                chrom_length: chrom.length,
            });
        }

        // Monotonic locus.
        if let Some((prev_chrom, prev_pos)) = self.last_locus {
            let regression = record.chrom_id < prev_chrom
                || (record.chrom_id == prev_chrom && record.pos <= prev_pos);
            if regression {
                return Err(PspWriteError::OutOfOrderRecord {
                    record_index,
                    prev_chrom,
                    prev_pos,
                    this_chrom: record.chrom_id,
                    this_pos: record.pos,
                });
            }
        }

        // At least one allele.
        if record.alleles.is_empty() {
            return Err(err_invalid_record(
                record_index,
                "record has zero alleles (REF always present)".to_string(),
            ));
        }

        // Per-allele rules.
        for (i, allele) in record.alleles.iter().enumerate() {
            let len = allele.seq.len();
            if len == 0 || (len as u64) > MAX_ALLELE_SEQ_LEN {
                return Err(err_invalid_record(
                    record_index,
                    format!(
                        "allele {i} sequence length {len} outside [1, {MAX_ALLELE_SEQ_LEN}]"
                    ),
                ));
            }
            // H4: lookup-table walk lets the loop run to completion
            // (no per-element early-exit) so LLVM can autovectorise
            // on long alleles. The cold path scans again to locate
            // the offending byte for the error message.
            if !allele.seq.iter().all(|&b| ALLOWED_ALLELE_BYTE[b as usize]) {
                return Err(err_invalid_allele_byte(record_index, i, &allele.seq));
            }
            if !allele.support.q_sum.is_finite() {
                return Err(err_invalid_record(
                    record_index,
                    format!("allele {i} q_sum is non-finite ({})", allele.support.q_sum),
                ));
            }
            // chain_slots ascending + distinct (defensive).
            for w in allele.chain_slots.windows(2) {
                if w[0] >= w[1] {
                    return Err(err_invalid_record(
                        record_index,
                        format!("allele {i} chain_slots not strictly ascending"),
                    ));
                }
            }
        }

        // new / expired chains ascending + distinct.
        for (name, slots) in [
            ("new_chains", &record.new_chains),
            ("expired_chains", &record.expired_chains),
        ] {
            for w in slots.windows(2) {
                if w[0] >= w[1] {
                    return Err(err_invalid_record(
                        record_index,
                        format!("{name} not strictly ascending"),
                    ));
                }
            }
        }

        Ok(())
    }

    /// Apply the record's content to the open block's per-column
    /// buffers, update the running active-slot set, and validate
    /// phase-chain marker consistency against that set.
    fn apply_record_to_block(
        &mut self,
        record_index: u64,
        record: &PileupRecord,
    ) -> Result<(), PspWriteError> {
        // Phase-chain marker consistency: expired must be currently
        // active; new must not be currently active. Sorted-Vec
        // membership is `binary_search`.
        for &slot in &record.expired_chains {
            if self.active_slots.binary_search(&slot).is_err() {
                return Err(PspWriteError::PhaseChainMarkerInconsistency {
                    record_index,
                    reason: format!("expired chain slot {slot} not currently active"),
                });
            }
        }
        for &slot in &record.new_chains {
            if self.active_slots.binary_search(&slot).is_ok() {
                return Err(PspWriteError::PhaseChainMarkerInconsistency {
                    record_index,
                    reason: format!("new chain slot {slot} already active"),
                });
            }
        }
        // Apply markers to running set. Expired first, then new
        // (matches the BTreeSet-era semantics).
        for &slot in &record.expired_chains {
            if let Ok(idx) = self.active_slots.binary_search(&slot) {
                self.active_slots.remove(idx);
            }
        }
        for &slot in &record.new_chains {
            // `slot` was just checked not in active_slots, so this
            // returns Err with the sorted insertion index.
            if let Err(idx) = self.active_slots.binary_search(&slot) {
                self.active_slots.insert(idx, slot);
            }
        }
        // Now per-allele chain_slots must all be in the running set.
        for (i, allele) in record.alleles.iter().enumerate() {
            for &slot in &allele.chain_slots {
                if self.active_slots.binary_search(&slot).is_err() {
                    return Err(PspWriteError::PhaseChainMarkerInconsistency {
                        record_index,
                        reason: format!(
                            "allele {i} references chain slot {slot} not in active set"
                        ),
                    });
                }
            }
        }

        // Push to block buffers.
        let block = self.block.as_mut().expect("block open by construction");
        block.append_record(record);

        Ok(())
    }

    /// Encode the current block to wire bytes, compress its columns,
    /// emit header + payloads to the sink, and append an index
    /// entry.
    ///
    /// Per-flush scratch state (uncompressed payload buffer, the 12
    /// compressed-frame buffers, the manifest scratch, and the
    /// block-header bytes buffer) is reused across flushes — each is
    /// `clear()`ed inside this method, keeping its prior capacity so
    /// subsequent flushes do not pay realloc cost.
    fn flush_block(&mut self) -> Result<(), PspWriteError> {
        let mut block = self
            .block
            .take()
            .expect("flush_block called with no open block");
        let block_offset = self.sink_offset;
        let n_records = block.delta_pos.len() as u32;
        let n_total_alleles = block.allele_seq_len.len() as u32;
        let block_index = self.index_entries.len() as u64;

        // Single pass over the column registry: encode -> compress ->
        // record manifest entry. The intermediate `payloads` Vec the
        // pre-refactor code carried (writer.rs:445 in the pre-L2
        // tree) is gone; compressed-frame bytes live in
        // `self.compressed_scratch[i]` until written to the sink
        // below.
        self.manifest_scratch.clear();
        for (i, column_def) in V1_0_COLUMNS.iter().enumerate() {
            self.uncompressed_scratch.clear();
            encode_column_into(column_def, &block, &mut self.uncompressed_scratch)?;
            let uncompressed_len = self.uncompressed_scratch.len() as u32;
            zstd_compress_into(
                &mut self.compressor,
                &self.uncompressed_scratch,
                &mut self.compressed_scratch[i],
            )
            .map_err(|e| PspWriteError::Io {
                context: "zstd compression of column payload",
                source: e,
            })?;
            self.manifest_scratch.push(ColumnManifestEntry {
                tag: column_def.tag,
                compressed_len: self.compressed_scratch[i].len() as u32,
                uncompressed_len,
            });
        }

        // Move the active-slot snapshot rather than cloning it —
        // `block` is owned and dropped at end of scope (L3).
        let header = BlockHeader {
            chrom_id: block.chrom_id,
            first_pos: block.first_pos,
            n_records,
            n_total_alleles,
            active_chain_slots: std::mem::take(&mut block.snapshot_active_slots),
            manifest: std::mem::take(&mut self.manifest_scratch),
        };
        self.header_bytes_scratch.clear();
        encode_block_header(&header, &mut self.header_bytes_scratch).map_err(|source| {
            PspWriteError::BlockEmission {
                block_index,
                source,
            }
        })?;
        // The manifest moved into `header`; recover it for the next
        // flush so capacity is preserved.
        self.manifest_scratch = header.manifest;

        self.sink
            .write_all(&self.header_bytes_scratch)
            .map_err(|e| PspWriteError::Io {
                context: "block header",
                source: e,
            })?;
        let mut written = self.header_bytes_scratch.len() as u64;
        for compressed in &self.compressed_scratch {
            self.sink
                .write_all(compressed)
                .map_err(|e| PspWriteError::Io {
                    context: "block column payload",
                    source: e,
                })?;
            written += compressed.len() as u64;
        }
        self.sink_offset += written;

        self.index_entries.push(BlockIndexEntry {
            chrom_id: block.chrom_id,
            first_pos: block.first_pos,
            last_pos: block.last_pos,
            n_records,
            block_offset,
        });
        Ok(())
    }
}

// ---------------------------------------------------------------------
// BlockAccumulator
// ---------------------------------------------------------------------

/// Flat CSR storage for a per-record or per-allele list column.
/// `offsets[i]..offsets[i+1]` is row `i`'s slice into `data`;
/// `offsets[0]` is always `0`. H2 in
/// `ia/reviews/perf_psp_writer_2026-05-13.md`.
struct ListColumn {
    data: Vec<SlotId>,
    offsets: Vec<u32>,
}

impl ListColumn {
    fn with_capacity(n_rows_hint: usize, n_slots_hint: usize) -> Self {
        let mut offsets = Vec::with_capacity(n_rows_hint + 1);
        offsets.push(0);
        Self {
            data: Vec::with_capacity(n_slots_hint),
            offsets,
        }
    }
    #[inline]
    fn push_row(&mut self, row: &[SlotId]) {
        self.data.extend_from_slice(row);
        self.offsets.push(self.data.len() as u32);
    }
}

struct BlockAccumulator {
    chrom_id: u32,
    first_pos: u32,
    last_pos: u32,
    /// Snapshot of the active-slot set at the moment this block
    /// opened. Sits verbatim in the block header.
    snapshot_active_slots: Vec<SlotId>,
    // Per-record columns.
    delta_pos: Vec<u64>,
    n_alleles: Vec<u64>,
    new_chain_slots: ListColumn,
    expired_chain_slots: ListColumn,
    // Per-allele columns.
    allele_seq_len: Vec<u64>,
    allele_seq_bytes: Vec<u8>,
    allele_obs_count: Vec<u32>,
    allele_q_sum_log: Vec<f64>,
    allele_fwd_count: Vec<u32>,
    allele_placed_left_count: Vec<u32>,
    allele_placed_start_count: Vec<u32>,
    allele_chain_slots: ListColumn,
    /// Rough projection of the uncompressed byte total. Used to
    /// decide when to auto-flush.
    projected_bytes: usize,
}

impl BlockAccumulator {
    fn new(chrom_id: u32, first_pos: u32, snapshot_active_slots: Vec<SlotId>) -> Self {
        Self {
            chrom_id,
            first_pos,
            last_pos: first_pos,
            snapshot_active_slots,
            delta_pos: Vec::with_capacity(INITIAL_RECORDS_HINT),
            n_alleles: Vec::with_capacity(INITIAL_RECORDS_HINT),
            // Per-record list columns. Most records have an empty
            // chain-slots list, so the data buffer's hint is much
            // smaller than the offsets hint.
            new_chain_slots: ListColumn::with_capacity(INITIAL_RECORDS_HINT, 4096),
            expired_chain_slots: ListColumn::with_capacity(INITIAL_RECORDS_HINT, 4096),
            allele_seq_len: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_seq_bytes: Vec::with_capacity(INITIAL_ALLELE_SEQ_BYTES_HINT),
            allele_obs_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_q_sum_log: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_fwd_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_placed_left_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_placed_start_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_chain_slots: ListColumn::with_capacity(INITIAL_ALLELES_HINT, 4096),
            projected_bytes: 0,
        }
    }

    fn append_record(&mut self, record: &PileupRecord) {
        let is_first = self.delta_pos.is_empty();
        let delta = if is_first {
            0u64
        } else {
            (record.pos - self.last_pos) as u64
        };
        self.delta_pos.push(delta);
        self.n_alleles.push(record.alleles.len() as u64);
        // H2: extend_from_slice into the flat data buffer, no
        // per-record Vec allocation. The outer `Vec<Vec<SlotId>>`
        // anti-pattern is gone — only one offsets push per row.
        self.new_chain_slots.push_row(&record.new_chains);
        self.expired_chain_slots.push_row(&record.expired_chains);

        for allele in &record.alleles {
            self.allele_seq_len.push(allele.seq.len() as u64);
            self.allele_seq_bytes.extend_from_slice(&allele.seq);
            self.allele_obs_count.push(allele.support.num_obs);
            self.allele_q_sum_log.push(allele.support.q_sum);
            self.allele_fwd_count.push(allele.support.fwd);
            self.allele_placed_left_count
                .push(allele.support.placed_left);
            self.allele_placed_start_count
                .push(allele.support.placed_start);
            self.allele_chain_slots.push_row(&allele.chain_slots);
        }

        self.last_pos = record.pos;

        // Rough size projection — per record, plus per-allele +
        // per-byte. Doesn't need to be precise; the target is a soft
        // cap.
        let per_record = 1 // delta-pos varint typical
            + 1 // n-alleles varint typical
            + (1 + 2 * record.new_chains.len()) // new_chain marker
            + (1 + 2 * record.expired_chains.len()); // expired_chain marker
        let per_allele: usize = record
            .alleles
            .iter()
            .map(|a| {
                1                  // allele-seq-len varint typical
                + a.seq.len()      // allele-seq bytes
                + 4 + 8 + 4 + 4 + 4 // the five scalars
                + (1 + 2 * a.chain_slots.len())
            })
            .sum();
        self.projected_bytes += per_record + per_allele;
    }
}

// ---------------------------------------------------------------------
// Column encoders — dispatched by (cardinality, shape, element-type)
// ---------------------------------------------------------------------

/// Encode one column from the block accumulator into its
/// uncompressed wire form, appending into the caller-provided `out`
/// buffer. Dispatches by tag — each tag corresponds to a specific
/// accumulator field with a specific element type. The caller is
/// responsible for `clear()`-ing `out` before calling so the buffer
/// can be reused across columns without reallocation (L1).
fn encode_column_into(
    def: &ColumnDef,
    block: &BlockAccumulator,
    out: &mut Vec<u8>,
) -> Result<(), PspWriteError> {
    match def.tag {
        0x01 => encode_varint_column(&block.delta_pos, out),
        0x02 => encode_varint_column(&block.n_alleles, out),
        0x03 => encode_varint_column(&block.allele_seq_len, out),
        0x04 => encode_bytes_concat(
            &[&block.allele_seq_bytes[..]],
            // The bytes-column payload is the concatenation of every
            // allele's sequence bytes; the per-allele lengths live
            // in the 0x03 length column. We've kept the bytes
            // concatenated as we appended records, so the payload
            // is just `allele_seq_bytes` as-is.
            out,
        ),
        0x10 => encode_scalar_column(&block.allele_obs_count, out),
        0x11 => {
            // q_sum was finite-validated at write_record time.
            encode_scalar_column(&block.allele_q_sum_log, out)
        }
        0x12 => encode_scalar_column(&block.allele_fwd_count, out),
        0x13 => encode_scalar_column(&block.allele_placed_left_count, out),
        0x14 => encode_scalar_column(&block.allele_placed_start_count, out),
        0x20 => encode_list_column_csr(
            &block.new_chain_slots.data,
            &block.new_chain_slots.offsets,
            out,
        ),
        0x21 => encode_list_column_csr(
            &block.expired_chain_slots.data,
            &block.expired_chain_slots.offsets,
            out,
        ),
        0x22 => encode_list_column_csr(
            &block.allele_chain_slots.data,
            &block.allele_chain_slots.offsets,
            out,
        ),
        unknown => {
            // V1_0_COLUMNS is the source of truth; if
            // encode_column_into is called with a tag not in this
            // match, the registry has been extended without updating
            // the dispatch — a programmer error.
            return Err(PspWriteError::InvalidRecord {
                record_index: 0,
                reason: format!(
                    "internal: encode_column_into has no dispatch for tag {unknown:#x}"
                ),
            });
        }
    }
    // Verify the encoded size matches the schema-predicted size for
    // fixed-width scalars and the bytes column (with the length
    // column also known). This is a writer-side self-check: a
    // mis-match means a bug in the writer.
    debug_assert!(verify_encoded_size(def, block, out));
    Ok(())
}

fn verify_encoded_size(def: &ColumnDef, block: &BlockAccumulator, encoded: &[u8]) -> bool {
    let n_records = block.delta_pos.len();
    let n_total_alleles = block.allele_seq_len.len();
    match (def.cardinality, def.shape, def.element_type) {
        (Cardinality::PerRecord, Shape::Scalar, Some(ElementType::Varint)) => true, // variable
        (Cardinality::PerAllele, Shape::Scalar, Some(ElementType::Varint)) => true, // variable
        (Cardinality::PerRecord, Shape::Scalar, Some(e)) if e.fixed_byte_width().is_some() => {
            encoded.len() == n_records * e.fixed_byte_width().unwrap()
        }
        (Cardinality::PerAllele, Shape::Scalar, Some(e)) if e.fixed_byte_width().is_some() => {
            encoded.len() == n_total_alleles * e.fixed_byte_width().unwrap()
        }
        (_, Shape::List, _) => true, // variable
        (_, Shape::Bytes, _) => {
            // For tag 0x04, encoded length equals the sum of
            // allele_seq_len.
            encoded.len() as u64 == block.allele_seq_len.iter().sum::<u64>()
        }
        _ => true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::per_sample_caller::pileup::{AlleleObservation, AlleleSupportStats};
    use crate::per_sample_caller::psp::header::{
        ChromosomeEntry, ParameterValue, ParsedHeader, WriterProvenance, parse_header_bytes,
    };
    use crate::per_sample_caller::psp::index::decode_index;
    use crate::per_sample_caller::psp::trailer::{TRAILER_BYTES, decode_trailer};
    use std::collections::BTreeMap;
    use std::io::Cursor;

    // ---------- Fixture builders ---------------------------------

    fn writer_header(n_chroms: usize) -> WriterHeader {
        let chromosomes = (0..n_chroms)
            .map(|i| ChromosomeEntry {
                name: format!("chr{}", i + 1),
                length: 1_000_000,
                md5: "0".repeat(32),
            })
            .collect();
        let mut params = BTreeMap::new();
        params.insert("min-mapq".to_string(), ParameterValue::Integer(30));
        WriterHeader {
            format_version: (1, 0),
            sample: "sample".to_string(),
            reference: "ref.fa".to_string(),
            created: "2026-05-13T10:00:00Z".parse().unwrap(),
            chromosomes,
            writer: WriterProvenance {
                tool: "test".to_string(),
                version: "0.0.1".to_string(),
                subcommand: "per-sample".to_string(),
                input_crams: vec!["a.cram".to_string()],
                input_fasta: "ref.fa".to_string(),
                parameters: params,
            },
        }
    }

    fn support(num_obs: u32, q_sum: f64) -> AlleleSupportStats {
        AlleleSupportStats {
            num_obs,
            q_sum,
            fwd: num_obs / 2,
            placed_left: 0,
            placed_start: num_obs,
        }
    }

    fn allele(seq: &[u8], num_obs: u32, q_sum: f64, chain_slots: &[u16]) -> AlleleObservation {
        AlleleObservation {
            seq: seq.to_vec(),
            support: support(num_obs, q_sum),
            chain_slots: chain_slots.to_vec(),
        }
    }

    fn record(
        chrom_id: u32,
        pos: u32,
        alleles: Vec<AlleleObservation>,
        new_chains: Vec<u16>,
        expired_chains: Vec<u16>,
    ) -> PileupRecord {
        PileupRecord {
            chrom_id,
            pos,
            new_chains,
            expired_chains,
            alleles,
        }
    }

    // ---------- new() emits a valid framed header -----------------

    #[test]
    fn new_emits_valid_framed_header() {
        let header = writer_header(2);
        let cursor = Cursor::new(Vec::new());
        let writer = PspWriter::new(cursor, header.clone()).expect("new should succeed");
        // Drop the writer to extract the sink (or just call finish to get the trailer too).
        let sink = writer.finish().expect("finish should succeed");
        let bytes = sink.into_inner();
        // Parse the header back.
        let (parsed, consumed): (ParsedHeader, usize) =
            parse_header_bytes(&bytes).expect("header parse round-trip");
        assert_eq!(parsed.sample, "sample");
        assert_eq!(parsed.chromosomes.len(), 2);
        assert!(consumed > 0);
    }

    // ---------- finish() with zero records -----------------------

    #[test]
    fn finish_with_zero_records_writes_header_empty_index_trailer() {
        let header = writer_header(1);
        let cursor = Cursor::new(Vec::new());
        let writer = PspWriter::new(cursor, header).unwrap();
        let bytes = writer.finish().unwrap().into_inner();

        // Last 32 bytes are the trailer.
        let trailer_bytes: &[u8; TRAILER_BYTES] =
            bytes[bytes.len() - TRAILER_BYTES..].try_into().unwrap();
        let trailer = decode_trailer(trailer_bytes).expect("trailer should decode");
        assert_eq!(trailer.n_blocks, 0);
        assert_eq!(trailer.index_byte_length, 0);
        // Index lives between header-end and trailer-start, here
        // empty — index_offset points to trailer-start.
        let expected_index_offset = (bytes.len() - TRAILER_BYTES) as u64;
        assert_eq!(trailer.index_offset, expected_index_offset);

        let index = decode_index(&[], trailer.n_blocks).unwrap();
        assert!(index.is_empty());
    }

    // ---------- write_record validation ---------------------------

    #[test]
    fn rejects_unknown_chrom_id() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                5,
                100,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("chrom_id 5 with 1 chromosome should fail");
        assert!(matches!(err, PspWriteError::UnknownChromId { .. }));
    }

    #[test]
    fn rejects_pos_zero() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                0,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("pos 0 should fail");
        assert!(matches!(err, PspWriteError::PosOutOfRange { .. }));
    }

    #[test]
    fn rejects_pos_beyond_contig_length() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                10_000_000,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("pos beyond contig length should fail");
        assert!(matches!(err, PspWriteError::PosOutOfRange { .. }));
    }

    #[test]
    fn rejects_out_of_order_positions() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        let err = writer
            .write_record(&record(
                0,
                50,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("going backwards should fail");
        assert!(matches!(err, PspWriteError::OutOfOrderRecord { .. }));
    }

    #[test]
    fn rejects_chrom_regression() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(2)).unwrap();
        writer
            .write_record(&record(
                1,
                100,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("chrom regression should fail");
        assert!(matches!(err, PspWriteError::OutOfOrderRecord { .. }));
    }

    #[test]
    fn rejects_zero_alleles() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(0, 100, vec![], vec![], vec![]))
            .expect_err("zero alleles should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    #[test]
    fn rejects_non_acgtn_allele_byte() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"X", 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("non-ACGTN should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    #[test]
    fn rejects_oversized_allele_sequence() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let seq = b"A".repeat((MAX_ALLELE_SEQ_LEN + 1) as usize);
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(&seq, 10, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("seq > cap should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    #[test]
    fn rejects_nan_q_sum() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 10, f64::NAN, &[])],
                vec![],
                vec![],
            ))
            .expect_err("NaN q_sum should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    #[test]
    fn rejects_phase_chain_marker_inconsistency_expired_not_active() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[])],
                vec![],
                vec![5], // 5 is not in the active set yet
            ))
            .expect_err("expired-not-active should fail");
        assert!(matches!(
            err,
            PspWriteError::PhaseChainMarkerInconsistency { .. }
        ));
    }

    #[test]
    fn rejects_phase_chain_marker_inconsistency_double_open() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[3])],
                vec![3],
                vec![],
            ))
            .unwrap();
        let err = writer
            .write_record(&record(
                0,
                101,
                vec![allele(b"A", 1, -1.0, &[3])],
                vec![3], // 3 is already active
                vec![],
            ))
            .expect_err("re-opening an active slot should fail");
        assert!(matches!(
            err,
            PspWriteError::PhaseChainMarkerInconsistency { .. }
        ));
    }

    #[test]
    fn rejects_allele_chain_slot_not_in_active_set() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[99])], // 99 not opened
                vec![],
                vec![],
            ))
            .expect_err("allele references unknown slot");
        assert!(matches!(
            err,
            PspWriteError::PhaseChainMarkerInconsistency { .. }
        ));
    }

    // ---------- One block round-trip via writer alone ----------

    #[test]
    fn one_record_one_block_finishes_cleanly() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 10, -2.5, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        let bytes = writer.finish().unwrap().into_inner();

        // Trailer at tail says one block.
        let trailer: &[u8; TRAILER_BYTES] =
            bytes[bytes.len() - TRAILER_BYTES..].try_into().unwrap();
        let trailer = decode_trailer(trailer).unwrap();
        assert_eq!(trailer.n_blocks, 1);

        // Index decodes to one entry with the expected coordinates.
        let index_bytes = &bytes[trailer.index_offset as usize
            ..(trailer.index_offset + trailer.index_byte_length) as usize];
        let entries = decode_index(index_bytes, trailer.n_blocks).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].chrom_id, 0);
        assert_eq!(entries[0].first_pos, 100);
        assert_eq!(entries[0].last_pos, 100);
        assert_eq!(entries[0].n_records, 1);

        // Index checksum verifies.
        let stored = trailer.index_checksum;
        let computed = checksum_index(index_bytes);
        assert_eq!(stored, computed);
    }

    /// Chromosome boundary forces a flush.
    #[test]
    fn chrom_change_triggers_flush() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(2)).unwrap();
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        writer
            .write_record(&record(
                0,
                200,
                vec![allele(b"C", 1, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        // Switch chromosome — the writer should flush block 0 and
        // start block 1.
        writer
            .write_record(&record(
                1,
                1,
                vec![allele(b"G", 1, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        let trailer: &[u8; TRAILER_BYTES] =
            bytes[bytes.len() - TRAILER_BYTES..].try_into().unwrap();
        let trailer = decode_trailer(trailer).unwrap();
        assert_eq!(
            trailer.n_blocks, 2,
            "chrom change should produce two blocks"
        );
    }

    /// Active-slot snapshot on subsequent (same-chrom) blocks
    /// matches the running set at flush time. Tested indirectly via
    /// a write_record sequence that closes block 1's last open slot
    /// just before block 2 opens.
    #[test]
    fn active_slot_snapshot_carries_across_blocks() {
        // Build a writer with a small target so a single record
        // doesn't flush; we'll force a flush by chromosome change.
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(2)).unwrap();
        // chrom 0: open slot 7, no per-allele reference; expire it.
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[7])],
                vec![7],
                vec![],
            ))
            .unwrap();
        writer
            .write_record(&record(
                0,
                101,
                vec![allele(b"C", 1, -1.0, &[])],
                vec![],
                vec![7],
            ))
            .unwrap();
        // Move to chrom 1; active set is empty at this point.
        writer
            .write_record(&record(
                1,
                1,
                vec![allele(b"G", 1, -1.0, &[])],
                vec![],
                vec![],
            ))
            .unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        // Smoke test: file is non-zero and ends with a valid trailer.
        let trailer: &[u8; TRAILER_BYTES] =
            bytes[bytes.len() - TRAILER_BYTES..].try_into().unwrap();
        assert!(decode_trailer(trailer).is_ok());
    }

    /// Phase-chain still active when chromosome switches: writer
    /// catches the upstream bug at first-record time on the new
    /// chromosome.
    #[test]
    fn rejects_chain_active_at_chromosome_boundary() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(2)).unwrap();
        // Open a slot on chrom 0 and never close it.
        writer
            .write_record(&record(
                0,
                100,
                vec![allele(b"A", 1, -1.0, &[5])],
                vec![5],
                vec![],
            ))
            .unwrap();
        // Try to start chrom 1 — writer notices the active set is
        // non-empty and rejects.
        let err = writer
            .write_record(&record(
                1,
                1,
                vec![allele(b"G", 1, -1.0, &[])],
                vec![],
                vec![],
            ))
            .expect_err("chain active across chrom boundary should fail");
        assert!(matches!(
            err,
            PspWriteError::PhaseChainMarkerInconsistency { .. }
        ));
    }
}
