//! Streaming writer that turns a [`PileupRecord`] stream into a v1.0
//! `.psp` file.
//!
//! Three public entry points:
//!
//! - [`PspWriter::new`] — emits the framed TOML header to the sink
//!   and prepares the in-memory state for accepting records.
//! - [`PspWriter::write_record`] — validates one record and buffers
//!   it into the open block; auto-flushes when the projected
//!   uncompressed payload reaches [`TARGET_BLOCK_BYTES`] or
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
use super::errors::{InvalidRecordKind, PspWriteError};
use super::header::{WriterHeader, build_header_bytes};
use super::index::{BlockIndexEntry, checksum_index, encode_index};
use super::registry::{
    ColumnDef, ColumnKey, ColumnPayload, ElementType, MAX_ALLELE_SEQ_LEN, V1_0_COLUMNS,
};
use super::trailer::{Trailer, encode_trailer};
use crate::per_sample_pileup::pileup::{ChainId, PileupRecord};

/// Target uncompressed bytes per block. Spec §"Block sizing" pins
/// this; it is not a tunable. The writer auto-flushes when an open
/// block's projected payload reaches this value. Tests can override
/// via [`PspWriter::new_with_block_target`] (Mi1, M14).
pub const TARGET_BLOCK_BYTES: usize = 16 * 1024 * 1024;

/// Initial `Vec::with_capacity` hint for per-record columns in a
/// freshly-opened block (delta-pos, n-alleles, the per-record list
/// columns). Sized at the records-per-block high-water mark for the
/// SNP-typical bench workload (~31 uncompressed bytes per record),
/// where the 16 MiB block target fills at ~530k records; the hint
/// pre-allocates for the bulk so column `Vec`s do not realloc on
/// the hot path.
///
/// **Indel-heavy / multi-allelic workloads** carry more bytes per
/// record, so blocks flush at fewer records (5k–50k typical). In
/// that regime this hint *over-allocates* by ~5–10×; the excess
/// capacity is held for the block's lifetime and freed at flush.
/// No reallocation occurs in either direction — the trade is
/// "spend some peak RAM, save the hot-path doubling cost on SNP
/// runs."
const INITIAL_RECORDS_HINT: usize = 350_000;
/// Initial hint for per-allele columns. Sized a touch above
/// `INITIAL_RECORDS_HINT` to absorb the typical excess of alleles
/// over records (occasional multi-allelic positions). Same
/// SNP-vs-indel trade as `INITIAL_RECORDS_HINT`.
const INITIAL_ALLELES_HINT: usize = 360_000;
/// Initial hint for the concatenated allele-sequence byte buffer.
/// On the SNP-typical workload one byte per allele = ~360 KiB; on
/// indel-heavier workloads alleles run to ~10 bytes each. The 1 MiB
/// hint covers both without leaving much slack.
const INITIAL_ALLELE_SEQ_BYTES_HINT: usize = 1_000_000;
/// Initial chain-id-data capacity for the list-shaped chain-id
/// column (Mi4). Most records carry zero chain ids; the
/// SNP-typical workload's high-water mark sits well under this hint.
const INITIAL_CHAIN_IDS_HINT: usize = 4096;

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

/// `#[cold]` constructor for `PspWriteError::InvalidRecord` (L5).
/// Keeping the construction out-of-line shrinks the hot-path body
/// of `validate_record` (which the samply profile puts at 10 % self).
#[cold]
#[inline(never)]
fn err_invalid_record(record_index: u64, kind: InvalidRecordKind) -> PspWriteError {
    PspWriteError::InvalidRecord { record_index, kind }
}

#[cold]
#[inline(never)]
fn err_invalid_allele_byte(record_index: u64, allele_index: usize, seq: &[u8]) -> PspWriteError {
    // PANIC-FREE: validate_record only calls this helper after the
    // `seq.iter().all(|&b| ALLOWED_ALLELE_BYTE[b as usize])` check
    // at writer.rs:399 has returned `false` — so at least one
    // invalid byte exists in `seq`. find() returning None here
    // would mean the two walks disagree, which is a codebase bug.
    let (byte_offset, byte) = seq
        .iter()
        .enumerate()
        .find(|&(_, &b)| !ALLOWED_ALLELE_BYTE[b as usize])
        .map(|(j, &b)| (j, b))
        .expect("invariant: caller observed an invalid byte but find() returned None");
    PspWriteError::InvalidRecord {
        record_index,
        kind: InvalidRecordKind::InvalidAlleleByte {
            allele_index,
            byte_offset,
            byte,
        },
    }
}

/// Running per-record state shared across `write_record` calls.
/// M15 extracts this from [`PspWriter`] so the record-ingest path
/// owns its own bookkeeping without touching scratch buffers or
/// the sink.
struct IngestState {
    /// Last admitted record's `(chrom_id, pos)`, for monotonicity
    /// enforcement across `write_record` calls.
    last_locus: Option<(u32, u32)>,
    /// Running count of records the caller has handed us. Used in
    /// error messages so the offending input is identifiable.
    records_seen: u64,
}

impl IngestState {
    fn new() -> Self {
        Self {
            last_locus: None,
            records_seen: 0,
        }
    }
}

/// Reused per-flush scratch state. M15 hoists the zstd compressor
/// and the four buffer-sized scratches out of [`PspWriter`] so the
/// god-struct can shrink and the flush phases can be written as
/// free functions taking `&mut WriterScratch`.
struct WriterScratch {
    /// Persistent zstd compressor reused across every column of
    /// every block. Each `compress_to_buffer` call resets the
    /// internal frame state but keeps the CCtx workspace and tables
    /// alive.
    compressor: zstd::bulk::Compressor<'static>,
    /// Reused per-column uncompressed-payload scratch.
    uncompressed: Vec<u8>,
    /// Reused per-column compressed-frame buffers, indexed by
    /// `V1_0_COLUMNS` position. Each flush `clear()`s the inner
    /// Vecs (preserving their capacity) before refilling.
    compressed: Vec<Vec<u8>>,
    /// Reused per-flush manifest buffer (12 entries).
    manifest: Vec<ColumnManifestEntry>,
    /// Reused block-header serialisation buffer.
    header_bytes: Vec<u8>,
}

impl WriterScratch {
    fn new() -> Result<Self, PspWriteError> {
        let compressor = new_column_compressor().map_err(|e| PspWriteError::Io {
            context: "zstd compressor init",
            block_index: None,
            column_tag: None,
            source: e,
        })?;
        Ok(Self {
            compressor,
            uncompressed: Vec::new(),
            compressed: (0..V1_0_COLUMNS.len()).map(|_| Vec::new()).collect(),
            manifest: Vec::with_capacity(V1_0_COLUMNS.len()),
            header_bytes: Vec::new(),
        })
    }
}

/// Streaming `.psp` writer over any `Write` sink. M15: down from
/// 14 fields to 7 by extracting [`IngestState`] and
/// [`WriterScratch`].
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
    /// Running per-record ingest state.
    ingest: IngestState,
    /// Reused per-flush scratch buffers + zstd compressor.
    scratch: WriterScratch,
    /// Auto-flush trigger: the writer flushes when an open block's
    /// projected uncompressed bytes reach this value. Defaults to
    /// [`TARGET_BLOCK_BYTES`]; tests override via
    /// [`PspWriter::new_with_block_target`] (M14, Mi26).
    target_block_bytes: usize,
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
    pub fn new(sink: W, header: WriterHeader) -> Result<Self, PspWriteError> {
        Self::new_with_block_target(sink, header, TARGET_BLOCK_BYTES)
    }

    /// Like [`Self::new`] but overrides the auto-flush projected-byte
    /// threshold. Exposed (hidden from rustdoc) so tests can force
    /// size-driven flushes with realistic record counts. Not for
    /// production use — the spec pins the block target at
    /// [`TARGET_BLOCK_BYTES`].
    #[doc(hidden)]
    pub fn new_with_block_target(
        mut sink: W,
        header: WriterHeader,
        target_block_bytes: usize,
    ) -> Result<Self, PspWriteError> {
        let header_bytes = build_header_bytes(&header)?;
        sink.write_all(&header_bytes)
            .map_err(|e| PspWriteError::Io {
                context: "file header",
                block_index: None,
                column_tag: None,
                source: e,
            })?;
        let sink_offset = header_bytes.len() as u64;
        let scratch = WriterScratch::new()?;
        Ok(Self {
            sink,
            header,
            sink_offset,
            index_entries: Vec::new(),
            block: None,
            ingest: IngestState::new(),
            scratch,
            target_block_bytes,
        })
    }

    /// Append one record. Returns the number of bytes pushed to the
    /// sink as a side effect of any auto-flush this call triggered
    /// (zero if no flush).
    pub fn write_record(&mut self, record: &PileupRecord) -> Result<u64, PspWriteError> {
        let record_index = self.ingest.records_seen;
        self.ingest.records_seen += 1;
        self.validate_record(record_index, record)?;

        // Should we flush before appending? Flush triggers:
        // - chrom_id change (blocks never cross chromosomes)
        // - projected uncompressed size at or past the target
        let pre_flush_bytes = self.sink_offset;
        let should_flush = match &self.block {
            None => false,
            Some(b) => {
                b.chrom_id != record.chrom_id || b.projected_bytes >= self.target_block_bytes
            }
        };
        if should_flush {
            self.flush_block()?;
        }
        let flushed_bytes = self.sink_offset - pre_flush_bytes;

        // Start a new block if none is open. With unique-per-file
        // chain ids there is no active-slot snapshot to copy; the
        // block header only carries structural metadata
        // (chrom_id, first_pos, n_records, n_total_alleles, manifest).
        if self.block.is_none() {
            self.block = Some(BlockAccumulator::new(record.chrom_id, record.pos));
        }

        self.apply_record_to_block(record_index, record)?;
        self.ingest.last_locus = Some((record.chrom_id, record.pos));

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
                block_index: None,
                column_tag: None,
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
                block_index: None,
                column_tag: None,
                source: e,
            })?;
        self.sink_offset += trailer_bytes.len() as u64;
        self.sink.flush().map_err(|e| PspWriteError::Io {
            context: "final flush",
            block_index: None,
            column_tag: None,
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
        if let Some((prev_chrom, prev_pos)) = self.ingest.last_locus {
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
                InvalidRecordKind::ZeroAlleles,
            ));
        }

        // Per-allele rules.
        for (i, allele) in record.alleles.iter().enumerate() {
            let len = allele.seq.len();
            if len == 0 || (len as u64) > MAX_ALLELE_SEQ_LEN {
                return Err(err_invalid_record(
                    record_index,
                    InvalidRecordKind::AlleleSeqLen {
                        allele_index: i,
                        length: len,
                        max: MAX_ALLELE_SEQ_LEN,
                    },
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
                    InvalidRecordKind::NonFiniteQSum {
                        allele_index: i,
                        q_sum: allele.support.q_sum,
                    },
                ));
            }
            // chain_ids ascending + distinct (defensive).
            for w in allele.chain_ids.windows(2) {
                if w[0] >= w[1] {
                    return Err(err_invalid_record(
                        record_index,
                        InvalidRecordKind::AlleleChainIdsNotAscending { allele_index: i },
                    ));
                }
            }
        }

        Ok(())
    }

    /// Apply the record's content to the open block's per-column
    /// buffers.
    fn apply_record_to_block(
        &mut self,
        record_index: u64,
        record: &PileupRecord,
    ) -> Result<(), PspWriteError> {
        // With unique-per-file chain ids there's no active-set
        // bookkeeping or lifecycle-marker validation. Each
        // `allele.chain_ids` is just a list of `u64`s that must
        // be strictly ascending (per-record well-formedness — pinned
        // in `validate_record`); identifier collisions are
        // structurally impossible because the chain-id allocator
        // monotonically mints distinct values per file.
        //
        // PANIC-FREE: write_record opens `self.block` on every path
        // that reaches this call. The Option shape is an artefact
        // of the flush/open cycle, not a real "may be absent"
        // condition here.
        let _ = record_index;
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
        // PANIC-FREE: flush_block is called only from `write_record`
        // (gated by `should_flush`, which requires `self.block` to
        // already be `Some`) and from `finish` (gated by
        // `if self.block.is_some()`). Both call sites ensure a block
        // is open.
        let block = self
            .block
            .take()
            .expect("flush_block called with no open block");
        let block_offset = self.sink_offset;
        let n_records = block.delta_pos.len() as u32;
        let n_total_alleles = block.allele_seq_len.len() as u32;
        let block_index = self.index_entries.len() as u64;

        // Mi30: split flush into three phases. Each is a free
        // function taking `&mut WriterScratch` (and, where needed,
        // the open block, the sink, and `block_index` for error
        // context).
        encode_and_compress_columns(&block, &mut self.scratch, block_index)?;
        assemble_block_header(
            &block,
            &mut self.scratch,
            block_index,
            n_records,
            n_total_alleles,
        )?;
        let written = emit_block_to_sink(&mut self.sink, &self.scratch, block_index)?;
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
// Flush phases (Mi30)
// ---------------------------------------------------------------------

/// Phase 1: walk the column registry, encoding each column's
/// uncompressed bytes, compressing them into `scratch.compressed[i]`,
/// and recording a manifest entry. `block_index` is used only for
/// error context.
fn encode_and_compress_columns(
    block: &BlockAccumulator,
    scratch: &mut WriterScratch,
    block_index: u64,
) -> Result<(), PspWriteError> {
    scratch.manifest.clear();
    for (i, column_def) in V1_0_COLUMNS.iter().enumerate() {
        scratch.uncompressed.clear();
        encode_column_into(column_def, block, &mut scratch.uncompressed)?;
        let uncompressed_len = scratch.uncompressed.len() as u32;
        zstd_compress_into(
            &mut scratch.compressor,
            &scratch.uncompressed,
            &mut scratch.compressed[i],
        )
        .map_err(|e| PspWriteError::Io {
            context: "zstd compression of column payload",
            block_index: Some(block_index),
            column_tag: Some(column_def.tag),
            source: e,
        })?;
        scratch.manifest.push(ColumnManifestEntry {
            tag: column_def.tag,
            compressed_len: scratch.compressed[i].len() as u32,
            uncompressed_len,
        });
    }
    Ok(())
}

/// Phase 2: build the [`BlockHeader`] struct (moving the manifest
/// scratch into it) and serialise it into `scratch.header_bytes`.
/// The manifest scratch is restored at the end so capacity carries
/// over to the next flush.
fn assemble_block_header(
    block: &BlockAccumulator,
    scratch: &mut WriterScratch,
    block_index: u64,
    n_records: u32,
    n_total_alleles: u32,
) -> Result<(), PspWriteError> {
    let header = BlockHeader {
        chrom_id: block.chrom_id,
        first_pos: block.first_pos,
        n_records,
        n_total_alleles,
        manifest: std::mem::take(&mut scratch.manifest),
    };
    scratch.header_bytes.clear();
    encode_block_header(&header, &mut scratch.header_bytes)
        .map_err(|kind| PspWriteError::BlockEmission { block_index, kind })?;
    // Recover the manifest for the next flush so capacity is
    // preserved.
    scratch.manifest = header.manifest;
    Ok(())
}

/// Phase 3: write the encoded block header followed by every
/// compressed column payload to the sink. Returns the total number
/// of bytes written so the caller can update `sink_offset`.
fn emit_block_to_sink<W: Write>(
    sink: &mut W,
    scratch: &WriterScratch,
    block_index: u64,
) -> Result<u64, PspWriteError> {
    sink.write_all(&scratch.header_bytes)
        .map_err(|e| PspWriteError::Io {
            context: "block header",
            block_index: Some(block_index),
            column_tag: None,
            source: e,
        })?;
    let mut written = scratch.header_bytes.len() as u64;
    for (i, compressed) in scratch.compressed.iter().enumerate() {
        sink.write_all(compressed).map_err(|e| PspWriteError::Io {
            context: "block column payload",
            block_index: Some(block_index),
            column_tag: Some(V1_0_COLUMNS[i].tag),
            source: e,
        })?;
        written += compressed.len() as u64;
    }
    Ok(written)
}

// ---------------------------------------------------------------------
// BlockAccumulator
// ---------------------------------------------------------------------

/// Flat CSR storage for a per-record or per-allele list column.
/// `offsets[i]..offsets[i+1]` is row `i`'s slice into `data`;
/// `offsets[0]` is always `0`. H2 in
/// `ia/reviews/perf_psp_writer_2026-05-13.md`.
struct ListColumn {
    data: Vec<ChainId>,
    offsets: Vec<u32>,
}

impl ListColumn {
    fn with_capacity(n_rows_hint: usize, n_chain_ids_hint: usize) -> Self {
        let mut offsets = Vec::with_capacity(n_rows_hint + 1);
        offsets.push(0);
        Self {
            data: Vec::with_capacity(n_chain_ids_hint),
            offsets,
        }
    }
    #[inline]
    fn push_row(&mut self, row: &[ChainId]) {
        self.data.extend_from_slice(row);
        self.offsets.push(self.data.len() as u32);
    }
}

struct BlockAccumulator {
    chrom_id: u32,
    first_pos: u32,
    last_pos: u32,
    // Per-record columns.
    delta_pos: Vec<u64>,
    n_alleles: Vec<u64>,
    // Per-allele columns.
    allele_seq_len: Vec<u64>,
    allele_seq_bytes: Vec<u8>,
    allele_obs_count: Vec<u32>,
    allele_q_sum_log: Vec<f64>,
    allele_fwd_count: Vec<u32>,
    allele_placed_left_count: Vec<u32>,
    allele_placed_start_count: Vec<u32>,
    allele_mapq_sum: Vec<u32>,
    allele_mapq_sum_sq: Vec<u64>,
    allele_chain_ids: ListColumn,
    /// Rough projection of the uncompressed byte total. Used to
    /// decide when to auto-flush.
    projected_bytes: usize,
}

impl BlockAccumulator {
    fn new(chrom_id: u32, first_pos: u32) -> Self {
        Self {
            chrom_id,
            first_pos,
            last_pos: first_pos,
            delta_pos: Vec::with_capacity(INITIAL_RECORDS_HINT),
            n_alleles: Vec::with_capacity(INITIAL_RECORDS_HINT),
            allele_seq_len: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_seq_bytes: Vec::with_capacity(INITIAL_ALLELE_SEQ_BYTES_HINT),
            allele_obs_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_q_sum_log: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_fwd_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_placed_left_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_placed_start_count: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_mapq_sum: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_mapq_sum_sq: Vec::with_capacity(INITIAL_ALLELES_HINT),
            allele_chain_ids: ListColumn::with_capacity(
                INITIAL_ALLELES_HINT,
                INITIAL_CHAIN_IDS_HINT,
            ),
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
            self.allele_mapq_sum.push(allele.support.mapq_sum);
            self.allele_mapq_sum_sq.push(allele.support.mapq_sum_sq);
            self.allele_chain_ids.push_row(&allele.chain_ids);
        }

        self.last_pos = record.pos;

        // Rough size projection — per record, plus per-allele +
        // per-byte. Doesn't need to be precise; the target is a soft
        // cap. Chain ids are u64 little-endian (8 bytes/id);
        // zstd compresses the high-order zero bytes effectively.
        let per_record = 1 // delta-pos varint typical
            + 1; // n-alleles varint typical
        let per_allele: usize = record
            .alleles
            .iter()
            .map(|a| {
                1                       // allele-seq-len varint typical
                + a.seq.len()           // allele-seq bytes
                + 4 + 8 + 4 + 4 + 4     // num_obs / q_sum / fwd / placed_left / placed_start
                + 4 + 8                 // mapq_sum (u32) + mapq_sum_sq (u64)
                + (1 + 8 * a.chain_ids.len()) // varint count + u64 ids
            })
            .sum();
        self.projected_bytes += per_record + per_allele;
    }
}

// ---------------------------------------------------------------------
// Column encoders — exhaustive dispatch on `ColumnKey` (M4)
// ---------------------------------------------------------------------

/// Encode one column from the block accumulator into its
/// uncompressed wire form, appending into the caller-provided `out`
/// buffer. Dispatches exhaustively on `def.key`, so adding a
/// `ColumnKey` variant forces an arm here as a compile error — no
/// runtime fall-through. The caller is responsible for `clear()`-ing
/// `out` before calling so the buffer can be reused across columns
/// without reallocation (L1).
fn encode_column_into(
    def: &ColumnDef,
    block: &BlockAccumulator,
    out: &mut Vec<u8>,
) -> Result<(), PspWriteError> {
    match def.key {
        ColumnKey::DeltaPos => encode_varint_column(&block.delta_pos, out),
        ColumnKey::NAlleles => encode_varint_column(&block.n_alleles, out),
        ColumnKey::AlleleSeqLen => encode_varint_column(&block.allele_seq_len, out),
        ColumnKey::AlleleSeq => encode_bytes_concat(
            // The bytes-column payload is the concatenation of every
            // allele's sequence bytes; the per-allele lengths live
            // in the `allele-seq-len` column. We kept the bytes
            // concatenated as we appended records, so the payload
            // is just `allele_seq_bytes` as-is.
            &[&block.allele_seq_bytes[..]],
            out,
        ),
        ColumnKey::AlleleObsCount => encode_scalar_column(&block.allele_obs_count, out),
        ColumnKey::AlleleQSumLog => {
            // q_sum was finite-validated at write_record time.
            encode_scalar_column(&block.allele_q_sum_log, out)
        }
        ColumnKey::AlleleFwdCount => encode_scalar_column(&block.allele_fwd_count, out),
        ColumnKey::AllelePlacedLeftCount => {
            encode_scalar_column(&block.allele_placed_left_count, out)
        }
        ColumnKey::AllelePlacedStartCount => {
            encode_scalar_column(&block.allele_placed_start_count, out)
        }
        ColumnKey::AlleleMapqSum => encode_scalar_column(&block.allele_mapq_sum, out),
        ColumnKey::AlleleMapqSumSq => encode_scalar_column(&block.allele_mapq_sum_sq, out),
        ColumnKey::AlleleChainIds => encode_list_column_csr(
            &block.allele_chain_ids.data,
            &block.allele_chain_ids.offsets,
            out,
        ),
    }
    // M5: promote the writer-side column-size self-check from a
    // debug-only assertion to a real check. Catches a writer bug
    // (encoded length disagrees with the schema-predicted length)
    // at the producer, not weeks later when a consumer rejects the
    // file. The cost is one length comparison per column per block
    // flush — well below the noise floor of the block's
    // compression cost.
    let predicted = predict_uncompressed_len(def, block);
    if let Some(expected) = predicted
        && expected != out.len()
    {
        return Err(PspWriteError::ColumnSizeSelfCheck {
            column: def.name,
            got: out.len(),
            expected,
        });
    }
    Ok(())
}

/// Predict the encoded byte length for a column whose shape gives
/// us an a-priori bound. Returns `None` for columns whose encoded
/// length is genuinely variable (varint scalars, list columns).
///
/// The exhaustive match on `ColumnPayload` (and on the inner
/// `ElementType` for fixed-width scalars) forces a compile-time
/// update when a new variant lands — replacing the prior
/// `_ => true` catch-all that silently approved any future
/// combination.
fn predict_uncompressed_len(def: &ColumnDef, block: &BlockAccumulator) -> Option<usize> {
    use super::registry::Cardinality;
    let n_records = block.delta_pos.len();
    let n_total_alleles = block.allele_seq_len.len();
    match def.payload {
        ColumnPayload::Bytes { .. } => {
            // `MAX_ALLELE_SEQ_LEN = 10_000`; with `u32::MAX` alleles
            // per block the sum still fits comfortably in `usize` on
            // 64-bit targets, but be defensive.
            let total: u64 = block.allele_seq_len.iter().sum();
            usize::try_from(total).ok()
        }
        ColumnPayload::List { .. } => None,
        ColumnPayload::Scalar { element_type: et } => {
            let count = match def.cardinality {
                Cardinality::PerRecord => n_records,
                Cardinality::PerAllele => n_total_alleles,
            };
            match et {
                ElementType::Varint | ElementType::Svarint => None,
                ElementType::U8
                | ElementType::U16
                | ElementType::U32
                | ElementType::U64
                | ElementType::I32
                | ElementType::I64
                | ElementType::F32
                | ElementType::F64
                | ElementType::Bool => Some(
                    count
                        * et.fixed_byte_width()
                            .expect("fixed-width path: et was checked above"),
                ),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::per_sample_pileup::pileup::{AlleleObservation, AlleleSupportStats};
    use crate::per_sample_pileup::psp::header::{ParsedHeader, parse_header_bytes};
    use crate::per_sample_pileup::psp::index::decode_index;
    use crate::per_sample_pileup::psp::trailer::{TRAILER_BYTES, decode_trailer};
    use std::io::Cursor;

    // ---------- Fixture builders ---------------------------------
    //
    // Mi20: the `writer_header` fixture lives in the shared
    // `super::super::test_fixtures` module so a new mandatory
    // `WriterHeader` field updates one place.
    use super::super::test_fixtures::writer_header;

    fn support(num_obs: u32, q_sum: f64) -> AlleleSupportStats {
        AlleleSupportStats {
            num_obs,
            q_sum,
            fwd: num_obs / 2,
            placed_left: 0,
            placed_start: num_obs,
        
            mapq_sum: 0,
            mapq_sum_sq: 0,
        }
    }

    fn allele(seq: &[u8], num_obs: u32, q_sum: f64, chain_ids: &[u64]) -> AlleleObservation {
        AlleleObservation {
            seq: seq.to_vec(),
            support: support(num_obs, q_sum),
            chain_ids: chain_ids.to_vec(),
        }
    }

    fn record(chrom_id: u32, pos: u32, alleles: Vec<AlleleObservation>) -> PileupRecord {
        PileupRecord {
            chrom_id,
            pos,
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
            .write_record(&record(5, 100, vec![allele(b"A", 10, -1.0, &[])]))
            .expect_err("chrom_id 5 with 1 chromosome should fail");
        assert!(matches!(err, PspWriteError::UnknownChromId { .. }));
    }

    #[test]
    fn rejects_pos_zero() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(0, 0, vec![allele(b"A", 10, -1.0, &[])]))
            .expect_err("pos 0 should fail");
        assert!(matches!(err, PspWriteError::PosOutOfRange { .. }));
    }

    #[test]
    fn rejects_pos_beyond_contig_length() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(0, 10_000_000, vec![allele(b"A", 10, -1.0, &[])]))
            .expect_err("pos beyond contig length should fail");
        assert!(matches!(err, PspWriteError::PosOutOfRange { .. }));
    }

    #[test]
    fn rejects_out_of_order_positions() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer
            .write_record(&record(0, 100, vec![allele(b"A", 10, -1.0, &[])]))
            .unwrap();
        let err = writer
            .write_record(&record(0, 50, vec![allele(b"A", 10, -1.0, &[])]))
            .expect_err("going backwards should fail");
        assert!(matches!(err, PspWriteError::OutOfOrderRecord { .. }));
    }

    #[test]
    fn rejects_chrom_regression() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(2)).unwrap();
        writer
            .write_record(&record(1, 100, vec![allele(b"A", 10, -1.0, &[])]))
            .unwrap();
        let err = writer
            .write_record(&record(0, 100, vec![allele(b"A", 10, -1.0, &[])]))
            .expect_err("chrom regression should fail");
        assert!(matches!(err, PspWriteError::OutOfOrderRecord { .. }));
    }

    #[test]
    fn rejects_zero_alleles() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(0, 100, vec![]))
            .expect_err("zero alleles should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    #[test]
    fn rejects_non_acgtn_allele_byte() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(0, 100, vec![allele(b"X", 10, -1.0, &[])]))
            .expect_err("non-ACGTN should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    #[test]
    fn rejects_oversized_allele_sequence() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let seq = b"A".repeat((MAX_ALLELE_SEQ_LEN + 1) as usize);
        let err = writer
            .write_record(&record(0, 100, vec![allele(&seq, 10, -1.0, &[])]))
            .expect_err("seq > cap should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    #[test]
    fn rejects_nan_q_sum() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        let err = writer
            .write_record(&record(0, 100, vec![allele(b"A", 10, f64::NAN, &[])]))
            .expect_err("NaN q_sum should fail");
        assert!(matches!(err, PspWriteError::InvalidRecord { .. }));
    }

    // Phase-chain marker inconsistency tests removed: the writer no
    // longer tracks an active-chain set, so the marker-vs-active
    // validation it used to enforce is gone. Chain ids are unique
    // per `.psp` file (u64), so the same family of violations is
    // structurally impossible. The remaining per-record well-
    // formedness check (chain_ids strictly ascending) is pinned
    // by the iterator-style assertions inside the round-trip tests.

    // ---------- One block round-trip via writer alone ----------

    #[test]
    fn one_record_one_block_finishes_cleanly() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer
            .write_record(&record(0, 100, vec![allele(b"A", 10, -2.5, &[])]))
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
            .write_record(&record(0, 100, vec![allele(b"A", 1, -1.0, &[])]))
            .unwrap();
        writer
            .write_record(&record(0, 200, vec![allele(b"C", 1, -1.0, &[])]))
            .unwrap();
        // Switch chromosome — the writer should flush block 0 and
        // start block 1.
        writer
            .write_record(&record(1, 1, vec![allele(b"G", 1, -1.0, &[])]))
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

    /// (M14.) The size-driven auto-flush trigger fires when an open
    /// block's projected bytes cross `target_block_bytes`. We use
    /// the `#[cfg(test)]`-exposed knob to override the 16 MiB
    /// default so a modestly-sized fixture exercises the branch.
    #[test]
    fn auto_flushes_on_projected_size_boundary() {
        // Target ~8 KiB worth of projected bytes — each simple
        // record contributes ~30 bytes, so ~270 records per block.
        let target = 8 * 1024;
        let mut writer =
            PspWriter::new_with_block_target(Cursor::new(Vec::new()), writer_header(1), target)
                .unwrap();
        for i in 1u32..=2000 {
            writer
                .write_record(&record(0, i, vec![allele(b"A", 1, -1.0, &[])]))
                .unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();
        let trailer: &[u8; TRAILER_BYTES] =
            bytes[bytes.len() - TRAILER_BYTES..].try_into().unwrap();
        let trailer = decode_trailer(trailer).unwrap();
        assert!(
            trailer.n_blocks >= 2,
            "size-based auto-flush should produce multiple blocks, got {}",
            trailer.n_blocks
        );
    }

    /// (Mi26.) `write_record` returns the number of bytes pushed to
    /// the sink as a side effect of any auto-flush it triggered.
    /// Defending the contract against a future refactor of the
    /// `pre_flush_bytes` accounting.
    #[test]
    fn write_record_returns_flushed_byte_count() {
        // Tiny block target so the second batch of records forces a
        // flush.
        let target = 4 * 1024;
        let mut writer =
            PspWriter::new_with_block_target(Cursor::new(Vec::new()), writer_header(1), target)
                .unwrap();
        // Records up to ~the target are absorbed with no flush.
        let mut saw_flush = false;
        for i in 1u32..=2000 {
            let pushed = writer
                .write_record(&record(0, i, vec![allele(b"A", 1, -1.0, &[])]))
                .unwrap();
            if pushed > 0 {
                saw_flush = true;
                break;
            }
        }
        assert!(
            saw_flush,
            "expected at least one auto-flush before record 2000"
        );
        // The very next record after a flush is on a fresh block — it
        // must report 0 bytes pushed.
        let post_flush = writer
            .write_record(&record(0, 10_000, vec![allele(b"A", 1, -1.0, &[])]))
            .unwrap();
        assert_eq!(post_flush, 0, "post-flush record must not push bytes");
    }

    /// (Mi27.) Multi-block index entries match the input order and
    /// coordinates — defends against a bug that reverses entries,
    /// drops one, or computes a wrong `last_pos` / `block_offset`.
    #[test]
    fn index_entries_match_input_order_and_coordinates() {
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(2)).unwrap();
        writer
            .write_record(&record(0, 100, vec![allele(b"A", 1, -1.0, &[])]))
            .unwrap();
        writer
            .write_record(&record(0, 200, vec![allele(b"C", 1, -1.0, &[])]))
            .unwrap();
        writer
            .write_record(&record(1, 1, vec![allele(b"G", 1, -1.0, &[])]))
            .unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        let trailer: &[u8; TRAILER_BYTES] =
            bytes[bytes.len() - TRAILER_BYTES..].try_into().unwrap();
        let trailer = decode_trailer(trailer).unwrap();
        let index_bytes = &bytes[trailer.index_offset as usize
            ..(trailer.index_offset + trailer.index_byte_length) as usize];
        let entries = decode_index(index_bytes, trailer.n_blocks).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].chrom_id, 0);
        assert_eq!(entries[0].first_pos, 100);
        assert_eq!(entries[0].last_pos, 200);
        assert_eq!(entries[0].n_records, 2);
        assert_eq!(entries[1].chrom_id, 1);
        assert_eq!(entries[1].first_pos, 1);
        assert_eq!(entries[1].last_pos, 1);
        assert_eq!(entries[1].n_records, 1);
        assert!(
            entries[1].block_offset > entries[0].block_offset,
            "second block's offset must come after first"
        );
    }

    // The `active_slot_snapshot_carries_across_blocks` and
    // `rejects_chain_active_at_chromosome_boundary` tests are gone:
    // block headers no longer carry an active-slot snapshot, and
    // there's no active-set bookkeeping for the writer to validate
    // against at a chromosome boundary. Chain ids are unique per
    // file, so chains "still active across blocks" is just a normal
    // sequence of ids referenced by alleles in both blocks — nothing
    // special to check.
}
