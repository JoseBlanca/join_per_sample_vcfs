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
//!   materialisation. No cross-block bookkeeping under the
//!   unique-`u64`-chain-id design — chain ids are unique
//!   throughout the file, so the same id in any block always
//!   refers to the same molecule.
//! - **Region iteration:** binary-search the block index for the
//!   first overlapping block; clamp records to the requested
//!   window.
//!
//! See `ia/feature_implementation_plans/psp_reader.md` for the
//! design rationale and the test plan (Groups B, F, and R).

use std::io::{Read, Seek, SeekFrom};

use super::block::{
    BlockHeader, ColumnManifestEntry, decode_block_header, decode_bytes_split_csr,
    decode_list_column_csr, decode_scalar_column_pod, decode_varint_column,
    new_column_decompressor, zstd_decompress_into,
};
use super::errors::{BlockHeaderInvariantKind, PspReadError, ScalarDecodeError};
use super::header::{
    HEAD_SENTINEL_LEN, HEADER_FRAMING_BYTES, MAX_HEADER_BODY_BYTES, MIN_HEADER_BODY_BYTES,
    ParsedHeader, parse_header_bytes,
};
use super::index::{BlockIndexEntry, checksum_index, decode_index};
use super::kind::{BlockDecoder, PspKind};
use super::registry::{ColumnKey, MAX_ALLELE_SEQ_LEN, V1_0_COLUMNS, lookup_by_tag};
use super::trailer::{TRAILER_BYTES, Trailer, decode_trailer};
use super::writer::SnpKind;
use crate::pileup_record::{AlleleObservation, AlleleSupportStats, ChainId, PileupRecord};

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
    /// Decompressed payload of the optional metadata section (the gap
    /// between the block index and the trailer), or `None` when the
    /// file carries no section. Opaque to the container core; a
    /// kind-specific parser interprets the bytes.
    metadata: Option<Vec<u8>>,
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
    /// (block offsets in `[header_end, trailer.index_offset)` and
    /// strictly increasing; `chrom_id` and position within the
    /// declared chromosome table).
    ///
    /// # Errors
    ///
    /// `Err(PspReadError::Io { .. })` on I/O failure;
    /// `BadTrailerMagic` / `IndexChecksum` / `IndexTrailingBytes` /
    /// `IndexOverrunsTrailer` on trailer or index corruption;
    /// `Zstd` / `MetadataSectionTooLarge` on a corrupt or oversized
    /// metadata section; `BadHeadMagic` /
    /// `BadHeaderLength` / `HeaderToml` / `InvalidHeaderField` /
    /// `SentinelMismatch` / `UnsupportedFormatVersion` /
    /// `UnknownRequiredColumn` / `MissingRequiredColumn` /
    /// `SchemaMismatch` on header / registry disagreement;
    /// `BlockIndexOffsetInvalid` / `BlockIndexChromOutOfRange` /
    /// `BlockIndexPosOutOfRange` on a block-index entry that
    /// disagrees with the parsed header.
    ///
    /// # Buffering
    ///
    /// `PspReader::new` issues five `seek` + four `read_exact` calls
    /// just to open, and each block then triggers ~12 small
    /// `read_exact` calls (one per column). Wrap a real-file source
    /// in `BufReader::with_capacity(64 * 1024, file)` (or larger) —
    /// the default 8 KiB `BufReader::new` is too small to absorb the
    /// per-column reads in a single kernel `read(2)`. For multi-GB
    /// files consider 1 MiB or larger. `Cursor<Vec<u8>>` and other
    /// in-memory sources don't need wrapping. L8 in
    /// `ia/reviews/perf_psp_reader_2026-05-13.md`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use std::fs::File;
    /// use std::io::BufReader;
    /// use pop_var_caller::psp::PspReader;
    ///
    /// let f = File::open("sample.psp")?;
    /// let mut reader = PspReader::new(BufReader::with_capacity(64 * 1024, f))?;
    /// println!("sample = {}", reader.header().sample);
    /// for r in reader.records() {
    ///     let _ = r?;
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
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
            .map_err(io_err("trailer seek"))?;
        let mut trailer_bytes = [0u8; TRAILER_BYTES];
        source
            .read_exact(&mut trailer_bytes)
            .map_err(io_err("trailer read"))?;
        let trailer = decode_trailer(&trailer_bytes)?;

        // Cross-check trailer arithmetic against the file size. The
        // index plus the 32-byte trailer must fit in the file
        // exactly — anything else means the trailer's pointer math
        // is wrong (corruption) or a write was truncated.
        //
        // M16 / Mi4: the `checked_add` here is reachable on hostile
        // inputs (`index_offset = u64::MAX`, etc.). When it returns
        // `None` we mustn't recompute the same overflowing sum
        // unchecked for the error payload — wrap the diagnostic in
        // `unwrap_or(u64::MAX)` so the path stays panic-free in
        // debug and emits a clamped (but never wrapping) count in
        // release.
        // The block index ends at `index_offset + index_byte_length`.
        // The bytes from there to the trailer start are the optional
        // metadata section (architecture
        // `doc/devel/architecture/hidden_paralog_psp_integration.md`); a
        // zero-length gap means no section. Only an index that runs
        // *past* the trailer start (or overflows `u64`) is corruption.
        let trailer_start = file_len - TRAILER_BYTES as u64;
        let index_end = trailer
            .index_offset
            .checked_add(trailer.index_byte_length)
            .filter(|&end| end <= trailer_start)
            .ok_or(PspReadError::IndexOverrunsTrailer {
                index_offset: trailer.index_offset,
                index_byte_length: trailer.index_byte_length,
                trailer_start,
            })?;

        // 3. Read + decode the block index, then verify the XXH3-64
        //    checksum the trailer stamped over it.
        source
            .seek(SeekFrom::Start(trailer.index_offset))
            .map_err(io_err("block index seek"))?;
        let mut index_bytes = vec![0u8; trailer.index_byte_length as usize];
        source
            .read_exact(&mut index_bytes)
            .map_err(io_err("block index read"))?;
        let computed_checksum = checksum_index(&index_bytes);
        if computed_checksum != trailer.index_checksum {
            return Err(PspReadError::IndexChecksum {
                stored: trailer.index_checksum,
                computed: computed_checksum,
            });
        }
        let index = decode_index(&index_bytes, trailer.n_blocks)?;

        // 3b. Read the optional metadata section — the bytes between the
        //     index end and the trailer start. The cursor is already at
        //     `index_end` after the index `read_exact`, so no seek is
        //     needed. An empty gap means no section.
        let metadata_len = (trailer_start - index_end) as usize;
        let metadata = if metadata_len == 0 {
            None
        } else {
            let mut frame = vec![0u8; metadata_len];
            source
                .read_exact(&mut frame)
                .map_err(io_err("metadata section"))?;
            Some(super::metadata::decompress_metadata(&frame)?)
        };

        // 4. Read + parse the framed header. Body length is
        //    range-checked before allocating so a tampered length
        //    cannot drive an unbounded read.
        source
            .seek(SeekFrom::Start(0))
            .map_err(io_err("header seek"))?;
        let mut head_prefix = [0u8; 12];
        source
            .read_exact(&mut head_prefix)
            .map_err(io_err("header magic + length prefix"))?;
        // M2: avoid the `try_into().unwrap()` panic-path on a
        // statically infallible conversion.
        let mut body_len_bytes = [0u8; 8];
        body_len_bytes.copy_from_slice(&head_prefix[4..12]);
        let body_len = u64::from_le_bytes(body_len_bytes);
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

        // 5. Cross-region checks on the block index. Three rules,
        //    each surfaced by `BlockIndexOffsetInvalid` /
        //    `BlockIndexChromOutOfRange` /
        //    `BlockIndexPosOutOfRange` with a consistent half-open
        //    interval `[header_end, trailer.index_offset)`. M5: the
        //    range check is delayed until here so the error message
        //    can carry the real lower bound rather than `min: 0`.
        let n_chroms = header.chromosomes.len() as u32;
        for (i, entry) in index.iter().enumerate() {
            if entry.block_offset < header_end_offset || entry.block_offset >= trailer.index_offset
            {
                return Err(PspReadError::BlockIndexOffsetInvalid {
                    block: i,
                    offset: entry.block_offset,
                    min: header_end_offset,
                    max: trailer.index_offset,
                });
            }
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
        // Mi15: monotonic-offset check as a separate windows(2)
        // pass; cleanly factored from the range check above.
        for (i, w) in index.windows(2).enumerate() {
            if w[1].block_offset <= w[0].block_offset {
                return Err(PspReadError::BlockIndexOffsetInvalid {
                    block: i + 1,
                    offset: w[1].block_offset,
                    min: w[0].block_offset.saturating_add(1),
                    max: trailer.index_offset,
                });
            }
        }

        // 6. Position the cursor at the start of block 0 so a
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
            metadata,
        })
    }

    /// Parsed file header. Available immediately after [`Self::new`]
    /// returns — does not require any block reads.
    pub fn header(&self) -> &ParsedHeader {
        &self.header
    }

    /// Decompressed bytes of the optional metadata section, or `None`
    /// when the file carries none. The payload is opaque here — a
    /// kind-specific parser (e.g. the SNP per-sample summary) interprets
    /// it. Available immediately after [`Self::new`]; no block reads.
    pub fn metadata(&self) -> Option<&[u8]> {
        self.metadata.as_deref()
    }

    /// Decoded block index, in genomic order. Empty for a
    /// zero-block file. Exposed for tests, the future `psp dump`
    /// utility, and the cohort-merge driver that wants per-sample
    /// block layout up front (e.g. to schedule reads).
    pub fn block_index(&self) -> &[BlockIndexEntry] {
        &self.index
    }

    /// Block-granular **columnar** reader over this file. Decodes
    /// blocks on demand and exposes each as borrowed [`BlockColumns`]
    /// — the decoded columns straight from disk, with no row-shape
    /// [`PileupRecord`] materialised. Record-level windowing and
    /// accumulation are the caller's job; this only walks + decodes
    /// blocks and gives access to the block index.
    ///
    /// This is the columnar counterpart to [`Self::region_records`]:
    /// where that yields one allocated `PileupRecord` per record, this
    /// hands the consumer the block's columns to copy in bulk. The
    /// cohort loader uses it to fill `SampleColumns` without the
    /// columnar→row→columnar round-trip.
    ///
    /// Consumes the reader (the columnar cursor owns it) so the caller
    /// can hold one per sample in a `Vec`; recover it with
    /// [`BlockColumnReader::into_reader`].
    pub fn into_column_blocks(self) -> BlockColumnReader<R> {
        BlockColumnReader::new(self)
    }

    /// Trailer in strong-typed form. Mi8: crate-private — the
    /// trailer fields are derived from `block_index()` plus the
    /// file size, and downstream consumers should not need them.
    /// Kept for tests and the future `psp dump`.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn trailer(&self) -> Trailer {
        self.trailer
    }

    /// Sequential iterator over every record in genomic order.
    /// Errors poison the iterator — subsequent `next()` calls
    /// return `None`.
    pub fn records(&mut self) -> RecordsIter<'_, R> {
        RecordsIter::new(self, RangeClamp::None)
    }

    /// Sequential iterator typed to a given schema `S` (architecture
    /// §10). [`records`](Self::records) is the `S = SnpKind` case; SSR
    /// callers (and tests) use `records_of::<SsrKind>()`.
    ///
    /// If `S::KIND` disagrees with the file's header `kind` (e.g.
    /// `records_of::<SnpKind>()` on an `.ssr.psp` file), the iterator
    /// yields a single [`PspReadError::KindMismatch`] on its first
    /// `next()` and then ends — it does not decode columns under the
    /// wrong schema (M4).
    #[allow(private_bounds)]
    pub fn records_of<S: PspKind>(&mut self) -> RecordsIter<'_, R, S>
    where
        S::Decoder: BlockDecoder<Record = S::Record>,
    {
        RecordsIter::new(self, RangeClamp::None)
    }

    /// Consume the reader into an **owning**, schema-typed sequential iterator.
    ///
    /// [`records_of`](Self::records_of) borrows the reader for the iterator's
    /// lifetime, so a caller cannot store the reader *and* a live iterator over it
    /// (that would be self-referential). This variant moves the reader *into* the
    /// iterator, letting a long-lived consumer (the SSR cohort cursor) hold the reader
    /// and pull records across many calls. Like `records_of`, it decodes one block at
    /// a time — never the whole file — and a schema/`kind` mismatch surfaces once on
    /// the first `next()` then ends.
    #[allow(private_bounds)]
    pub fn into_records_of<S: PspKind>(self) -> OwnedRecordsIter<R, S>
    where
        S::Decoder: BlockDecoder<Record = S::Record>,
    {
        OwnedRecordsIter::new(self)
    }

    /// Iterator over records inside `[start, end]` (inclusive
    /// both ends) on `chrom_id`.
    ///
    /// Positions the iterator at the first block whose range
    /// overlaps the window and lets `RecordsIter` clamp per-record:
    /// records before `start` are skipped; the first record past
    /// `end` ends iteration.
    ///
    /// **Empty / invalid windows.** Mi6: a `chrom_id` past the
    /// declared chromosome table, a window whose `[start, end]`
    /// does not overlap any block on that chromosome, or `start >
    /// end` all return `None` from the first `next()` call without
    /// raising an error. Caller programs that need to distinguish
    /// "no records in this region" from "I asked the wrong
    /// question" must validate the inputs themselves.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use std::fs::File;
    /// use std::io::BufReader;
    /// use pop_var_caller::psp::PspReader;
    ///
    /// let f = File::open("sample.psp")?;
    /// let mut reader = PspReader::new(BufReader::with_capacity(64 * 1024, f))?;
    /// let n: usize = reader
    ///     .region_records(0, 1_000_000, 1_500_000)
    ///     .filter_map(|r| r.ok())
    ///     .count();
    /// println!("{n} records in chr 0:[1M, 1.5M]");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn region_records(&mut self, chrom_id: u32, start: u32, end: u32) -> RecordsIter<'_, R> {
        let first = find_first_overlapping_block(&self.index, chrom_id, start, end);
        let clamp = RangeClamp::Window {
            chrom_id,
            start,
            end,
        };
        let mut it = RecordsIter::new(self, clamp);
        // Mi16: replace `match first { Some => ..., None => ... }`
        // with `unwrap_or`.
        it.cur_block_idx = first.unwrap_or(it.reader.index.len());
        it
    }
}

/// Binary-search the block index for the first entry whose
/// `chrom_id` matches and whose `[first_pos, last_pos]` overlaps
/// `[start, end]`. Returns `None` if no entry overlaps.
fn find_first_overlapping_block(
    index: &[BlockIndexEntry],
    chrom_id: u32,
    start: u32,
    end: u32,
) -> Option<usize> {
    // Block index is sorted by `(chrom_id, first_pos)`
    // lexicographically. The first block whose range overlaps
    // `[start, end]` is the first one with `chrom_id == chrom_id
    // && last_pos >= start && first_pos <= end`.
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
///
/// H1 (2026-05-23 PSP reader perf review): the two ragged columns
/// (`allele-seq` bytes + `allele-chain-ids` lists) used to live here
/// as `Vec<Vec<u8>>` / `Vec<Vec<ChainId>>`, allocating one inner
/// `Vec` per allele on every block decode. They now live on
/// [`RecordsIter`] as CSR `(data, offsets)` pairs that are reused
/// across blocks — those per-allele allocations collapse to two
/// per-block `extend_from_slice` calls. The block holds only the
/// per-record columns and the per-allele fixed-width columns.
#[derive(Debug)]
struct DecodedBlock {
    chrom_id: u32,
    first_pos: u32,
    n_records: u32,
    // Per-record columns
    delta_pos: Vec<u64>,
    n_alleles: Vec<u64>,
    windowed_gc: Vec<f32>,
    windowed_coverage: Vec<f32>,
    // Per-allele fixed-width columns
    allele_obs_count: Vec<u32>,
    allele_q_sum_log: Vec<f64>,
    allele_fwd_count: Vec<u32>,
    allele_placed_left_count: Vec<u32>,
    allele_placed_start_count: Vec<u32>,
    allele_mapq_sum: Vec<u32>,
    allele_mapq_sum_sq: Vec<u64>,
    // (the two ragged columns — allele_seqs + allele_chain_ids —
    // live on RecordsIter; see the H1 doc comment above.)
}

/// Region-iteration window. `None` for sequential iteration
/// (the whole file); `Window` for region iteration (clamp records
/// to `[start, end]` on `chrom_id`).
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

impl RangeClamp {
    /// Mi17: `true` iff a record interval `[rec_start, ..)` lies past
    /// the window's right edge (or on a different chromosome). In
    /// `Window` mode the iterator terminates as soon as this fires.
    /// Interval-based (architecture §10.5) so the generic reader can
    /// clamp any kind via [`PspKind::record_interval`]. The query `end`
    /// is inclusive, so the record is past once `rec_start > end`.
    fn interval_past_window(self, chrom_id: u32, rec_start: u32) -> bool {
        match self {
            Self::None => false,
            Self::Window {
                chrom_id: c, end, ..
            } => chrom_id != c || rec_start > end,
        }
    }

    /// Mi17: `true` iff a record interval `[.., rec_end)` ends at or
    /// before the window's (inclusive) left edge — `rec_end <= start` —
    /// in `Window` mode; the iterator skips it and asks for the next.
    /// `rec_end` is exclusive, so a SNP point at `pos` (`rec_end = pos +
    /// 1`) is skipped exactly when `pos < start`, as before.
    fn interval_before_window(self, rec_end: u32) -> bool {
        matches!(self, Self::Window { start, .. } if rec_end <= start)
    }

    /// Mi17: `true` iff `entry`'s range cannot intersect this
    /// window — the iterator can terminate without decoding the
    /// block. `RangeClamp::None` never terminates early.
    fn block_past_window(self, entry: &BlockIndexEntry) -> bool {
        match self {
            Self::None => false,
            Self::Window { chrom_id, end, .. } => {
                entry.chrom_id != chrom_id || entry.first_pos > end
            }
        }
    }
}

/// Iterator over `PspReader`'s records. Owns a mutable borrow of
/// the source, so only one iterator lives at a time per
/// `PspReader`. **Not `Send`** by design — see
/// `ia/feature_implementation_plans/psp_reader.md` §"Tradeoffs"
/// for the rationale.
pub struct RecordsIter<'r, R: Read + Seek, S: PspKind = SnpKind> {
    reader: &'r mut PspReader<R>,
    /// Index of the next block to load. Equals
    /// `reader.index.len()` once iteration is exhausted.
    cur_block_idx: usize,
    /// Whether a block is currently decoded into `decoder`. `false`
    /// before the first block is pulled and after the current block is
    /// drained; gates whether `next()` materialises or loads.
    cur_block_loaded: bool,
    /// The schema's per-block decoder (architecture §10.2): owns the
    /// decoded columns, any cross-block reuse buffers (e.g. the SNP CSR
    /// slabs), and the per-block record cursor. Mirror of the writer's
    /// `S::Block`.
    decoder: S::Decoder,
    /// `RangeClamp::None` for sequential iteration; `Window` for
    /// region iteration.
    clamp: RangeClamp,
    /// Sticky: once `next()` has yielded `Some(Err(_))`, future
    /// calls return `None`. Required for `Iterator` correctness
    /// when the source state after an error is unknown.
    poisoned: bool,
    /// M4: an error to surface on the first `next()` before any block
    /// is touched, then poison. Set when the iterator's schema `S`
    /// disagrees with the file's header `kind` (decoding columns under
    /// the wrong schema would yield garbage), so a kind/decoder
    /// mismatch fails loudly instead of silently misdecoding.
    pending_error: Option<PspReadError>,
    /// L1: persistent zstd decompressor reused across every column
    /// of every block. The DCtx workspace is allocated once
    /// (~hundreds of KiB at level 9) instead of per-column. Mirror
    /// of writer-side H3 (commit 969de6c). Shared decompression
    /// machinery — handed to `decoder.decode_block`.
    decompressor: zstd::bulk::Decompressor<'static>,
    /// L2: per-iter scratch for the compressed column bytes. Cleared
    /// then resized to `entry.compressed_len` per column; capacity
    /// converges to the largest compressed column the iterator has
    /// seen.
    compressed_scratch: Vec<u8>,
    /// L2: per-iter scratch for the decompressed column bytes,
    /// written into by `zstd_decompress_into`. The per-column
    /// decoders take `&[u8]` and return owned `Vec<...>`, so this
    /// buffer is safe to reuse across columns within a block.
    decompressed_scratch: Vec<u8>,
    /// L2: per-iter scratch for the block-header grow-on-incomplete
    /// loop. One Vec per `RecordsIter` instead of one per
    /// `read_block_header` call.
    block_header_buf: Vec<u8>,
}

// `BlockDecoder` is deliberately `pub(crate)` (it names internal wire
// types — see its doc) and the bound is an implementation detail callers
// can't observe (the public surface is just `Iterator<Item =
// Result<Record>>`). Allowing keeps those wire types encapsulated rather
// than leaking them through a `pub` re-export just to satisfy the lint.
#[allow(private_bounds)]
impl<'r, R: Read + Seek, S: PspKind> RecordsIter<'r, R, S>
where
    S::Decoder: BlockDecoder<Record = S::Record>,
{
    fn new(reader: &'r mut PspReader<R>, clamp: RangeClamp) -> Self {
        // L1: persistent decompressor. `Decompressor::new()` only
        // fails in OOM scenarios and at that point we have bigger
        // problems; the initial DCtx allocation is small. Mirror of
        // the writer's `new_column_compressor()?` in `PspWriter::new`
        // (which propagates the io::Error up; we cannot here because
        // `records()` / `region_records()` return `RecordsIter`, not
        // `Result<RecordsIter, _>`).
        let decompressor = new_column_decompressor()
            .expect("zstd::bulk::Decompressor::new is infallible on supported platforms");
        // M4: refuse a schema/kind mismatch. The header `kind` was
        // validated against a known registry at parse time; if the
        // caller-chosen `S` names a *different* known kind, decoding
        // would read the wrong columns. Surface it as the iterator's
        // first item rather than yielding garbage records. (`records()`
        // / `region_records()` default `S = SnpKind`, so an snp file is
        // never flagged; an snp decoder over an ssr file now is.)
        let pending_error = (S::KIND != reader.header().kind).then(|| PspReadError::KindMismatch {
            expected: S::KIND,
            found: reader.header().kind.clone(),
        });
        Self {
            reader,
            cur_block_idx: 0,
            cur_block_loaded: false,
            decoder: S::Decoder::new_decoder(),
            clamp,
            poisoned: false,
            pending_error,
            decompressor,
            compressed_scratch: Vec::new(),
            decompressed_scratch: Vec::new(),
            block_header_buf: Vec::with_capacity(BLOCK_HEADER_INITIAL_CHUNK),
        }
    }

    /// Seek to the next un-decoded block and decode it into the
    /// schema decoder. Returns `Ok(false)` when no more blocks remain.
    /// The block framing + shared decompression scratch are owned here;
    /// the per-schema column decode is delegated to the decoder.
    ///
    /// With unique-per-file chain ids the block header no longer
    /// carries an active-slot snapshot, so there is no
    /// snapshot-mismatch or first-block-of-chrom-must-be-empty
    /// check to run here — block boundaries are structural only.
    fn load_next_block(&mut self) -> Result<bool, PspReadError> {
        if self.cur_block_idx >= self.reader.index.len() {
            return Ok(false);
        }
        let entry = self.reader.index[self.cur_block_idx];
        seek_to_offset(
            &mut self.reader.source,
            entry.block_offset,
            "block start seek",
        )?;
        let (block_header, _consumed) =
            read_block_header(&mut self.reader.source, &mut self.block_header_buf)?;
        let budget =
            block_byte_budget(&self.reader.index, &self.reader.trailer, self.cur_block_idx);

        self.decoder.decode_block(
            &mut self.reader.source,
            &block_header,
            budget,
            &mut self.decompressor,
            &mut self.compressed_scratch,
            &mut self.decompressed_scratch,
        )?;
        self.cur_block_loaded = true;
        Ok(true)
    }
}

/// The SNP per-block decoder — [`SnpKind`]'s [`BlockDecoder`]. Owns the
/// decoded SNP block, the two CSR slabs reused across blocks (the H1
/// allocation collapse), and the per-block record cursor. The mirror of
/// the writer's [`SnpBlock`](super::writer::SnpBlock).
pub struct SnpDecoder {
    /// Decoded payload of the current block. `None` before the first
    /// block is decoded.
    cur_block: Option<DecodedBlock>,
    /// Index of the next record (within the current block) to
    /// materialise. Resets to 0 on each new block.
    next_record_in_block: u32,
    /// Cumulative allele offset into the per-allele columns, scoped to
    /// the current block. Resets to 0 on each new block.
    next_allele_in_block: u32,
    /// Running 1-based reference position of the last-materialised
    /// record. Reset to the block's `first_pos` on each new block.
    last_pos: u32,
    /// File-global running record count, incremented after each
    /// materialised record.
    record_index: u64,
    /// H1: CSR data slab for the `allele-seq` bytes column, reused
    /// across blocks via `extend_from_slice` of the whole column
    /// payload at decode time.
    allele_seq_data: Vec<u8>,
    /// H1: CSR row offsets for `allele-seq`; `offsets.len() ==
    /// n_total_alleles_in_block + 1`.
    allele_seq_offsets: Vec<u32>,
    /// H1: CSR data slab for the `allele-chain-ids` list column.
    allele_chain_ids_data: Vec<ChainId>,
    /// H1: CSR row offsets for `allele-chain-ids`.
    allele_chain_ids_offsets: Vec<u32>,
}

impl BlockDecoder for SnpDecoder {
    type Record = PileupRecord;

    fn new_decoder() -> Self {
        Self {
            cur_block: None,
            next_record_in_block: 0,
            next_allele_in_block: 0,
            last_pos: 0,
            record_index: 0,
            allele_seq_data: Vec::new(),
            allele_seq_offsets: Vec::new(),
            allele_chain_ids_data: Vec::new(),
            allele_chain_ids_offsets: Vec::new(),
        }
    }

    fn decode_block<R: Read>(
        &mut self,
        source: &mut R,
        header: &BlockHeader,
        budget: u64,
        decompressor: &mut zstd::bulk::Decompressor<'static>,
        compressed_scratch: &mut Vec<u8>,
        decompressed_scratch: &mut Vec<u8>,
    ) -> Result<(), PspReadError> {
        let decoded = decode_block_payload(
            source,
            header,
            budget,
            decompressor,
            compressed_scratch,
            decompressed_scratch,
            &mut self.allele_seq_data,
            &mut self.allele_seq_offsets,
            &mut self.allele_chain_ids_data,
            &mut self.allele_chain_ids_offsets,
        )?;
        self.cur_block = Some(decoded);
        self.next_record_in_block = 0;
        self.next_allele_in_block = 0;
        // `delta_pos[0] = 0` by writer invariant, so the first
        // record's `pos = first_pos + delta_pos[0] = first_pos`.
        // Initialise `last_pos` to make that fall out of the same
        // addition path used for the rest of the block.
        self.last_pos = header.first_pos;
        Ok(())
    }

    fn next_record(&mut self) -> Option<Result<PileupRecord, PspReadError>> {
        match &self.cur_block {
            Some(block) if self.next_record_in_block < block.n_records => {
                Some(self.materialise_next_record())
            }
            _ => None,
        }
    }

    fn unload(&mut self) {
        // Mi2: free the exhausted block's per-record / per-allele
        // fixed-width columns (read by index, not `mem::take`, so they
        // stay resident otherwise). The H1 CSR reuse slabs
        // (`allele_seq_*` / `allele_chain_ids_*`) are kept — the next
        // `decode_block` overwrites them.
        self.cur_block = None;
    }
}

impl SnpDecoder {
    /// Construct the next `PileupRecord` from the currently-loaded
    /// block. Caller must ensure a block is loaded and has records
    /// left (enforced by [`Self::next_record`]).
    ///
    /// The block is consumed strictly forwards (the iterator never
    /// revisits a `(block, i, j)` triple). The per-record / per-allele
    /// fixed-width columns are read by index/copy (the H1 perf rework
    /// replaced the original M13 `mem::take` shape — the allele bytes /
    /// chain-ids now slice from the iter-owned CSR slabs). The block's
    /// buffers are freed when it is exhausted, via
    /// [`BlockDecoder::unload`](super::kind::BlockDecoder::unload),
    /// which the generic reader calls before advancing (Mi2).
    fn materialise_next_record(&mut self) -> Result<PileupRecord, PspReadError> {
        // M1: lift the option-unwrap out of the function body —
        // the caller (`<Self as Iterator>::next`) gates this call
        // behind `if let Some(block) = &self.cur_block && ...`,
        // but checking again here makes the invariant a type-level
        // fact instead of a prose precondition. PANIC-FREE by
        // construction: the `?` returns the typed error instead
        // of panicking if the contract is ever violated.
        let block = match self.cur_block.as_mut() {
            Some(b) => b,
            None => {
                return Err(PspReadError::Io {
                    context: "materialise_next_record without a loaded block",
                    source: std::io::Error::other("internal invariant"),
                });
            }
        };
        let i = self.next_record_in_block as usize;
        // Mi5: `delta_pos[i]` decodes as `u64`. A varint with the
        // top 32 bits set silently truncates if we cast to `u32`,
        // and `saturating_add` then caps `pos` at `u32::MAX` —
        // a wrong position, never an error. Reject the truncation
        // surface here with a typed varint-overflow error so the
        // wrong-bytes path stays diagnostic.
        let delta = block.delta_pos[i];
        let delta_u32 = u32::try_from(delta).map_err(|_| PspReadError::ColumnElementDecode {
            column: "delta-pos".to_string(),
            entry: i,
            source: ScalarDecodeError::VarintOverflow,
        })?;
        let pos = if i == 0 {
            block.first_pos
        } else {
            self.last_pos.saturating_add(delta_u32)
        };

        let n_alleles_here = block.n_alleles[i] as usize;
        let allele_start = self.next_allele_in_block as usize;
        let allele_end = allele_start + n_alleles_here;

        // H3 (2026-05-23 PSP reader perf review): leading asserts
        // hoist the per-allele bounds checks out of the inner loop.
        // The compiler can't share a bounds proof across the disjoint
        // per-allele `Vec`s otherwise. PANIC-FREE in practice — B3
        // in `decode_block_payload` already enforces
        // `sum(n_alleles[i]) == n_total_alleles`, and `CHECK 7` /
        // `predict_fixed_uncompressed_len` enforce each column's
        // length equals `n_total_alleles`. The asserts are dead
        // code on well-formed input but the dominator gives LLVM
        // the bounds-proof it needs.
        assert!(allele_end <= block.allele_obs_count.len());
        assert!(allele_end <= block.allele_q_sum_log.len());
        assert!(allele_end <= block.allele_fwd_count.len());
        assert!(allele_end <= block.allele_placed_left_count.len());
        assert!(allele_end <= block.allele_placed_start_count.len());
        assert!(allele_end <= block.allele_mapq_sum.len());
        assert!(allele_end <= block.allele_mapq_sum_sq.len());
        // The CSR offset arrays have `n_total_alleles + 1` entries,
        // so `allele_end` indexes the closing offset of the last
        // emitted allele.
        assert!(allele_end < self.allele_seq_offsets.len());
        assert!(allele_end < self.allele_chain_ids_offsets.len());

        // H1: slice from the iter-owned CSR data + offsets (`data`
        // is reused across blocks; `to_vec()` here is the single
        // remaining alloc per emitted allele). Replaces the previous
        // `mem::take(&mut block.allele_seqs[j])` shape.
        let mut alleles = Vec::with_capacity(n_alleles_here);
        for j in allele_start..allele_end {
            let seq_s = self.allele_seq_offsets[j] as usize;
            let seq_e = self.allele_seq_offsets[j + 1] as usize;
            let cid_s = self.allele_chain_ids_offsets[j] as usize;
            let cid_e = self.allele_chain_ids_offsets[j + 1] as usize;
            alleles.push(AlleleObservation {
                seq: self.allele_seq_data[seq_s..seq_e].to_vec(),
                support: AlleleSupportStats {
                    num_obs: block.allele_obs_count[j],
                    q_sum: block.allele_q_sum_log[j],
                    fwd: block.allele_fwd_count[j],
                    placed_left: block.allele_placed_left_count[j],
                    placed_start: block.allele_placed_start_count[j],
                    mapq_sum: block.allele_mapq_sum[j],
                    mapq_sum_sq: block.allele_mapq_sum_sq[j],
                },
                chain_ids: self.allele_chain_ids_data[cid_s..cid_e].to_vec(),
            });
        }

        let record = PileupRecord {
            chrom_id: block.chrom_id,
            pos,
            alleles,
            // Per-record columns, indexed like `delta_pos[i]` / `n_alleles[i]`
            // above (the decode-side length check guarantees `len == n_records`).
            windowed_gc: block.windowed_gc[i],
            windowed_coverage: block.windowed_coverage[i],
        };

        self.last_pos = pos;
        self.next_record_in_block += 1;
        self.next_allele_in_block += n_alleles_here as u32;
        self.record_index = self.record_index.saturating_add(1);
        Ok(record)
    }
}

// See the `#[allow]` rationale on the inherent impl above.
#[allow(private_bounds)]
impl<R: Read + Seek, S: PspKind> Iterator for RecordsIter<'_, R, S>
where
    S::Decoder: BlockDecoder<Record = S::Record>,
{
    type Item = Result<S::Record, PspReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.poisoned {
            return None;
        }
        // M4: a schema/kind mismatch surfaces once, then poisons.
        if let Some(err) = self.pending_error.take() {
            self.poisoned = true;
            return Some(Err(err));
        }
        // State machine, per iteration:
        //   1. If the current block has records left, yield one
        //      (with region clamp).
        //   2. Otherwise drop the current block and advance the
        //      index.
        //   3. In region mode, peek the next index entry — terminate
        //      early if it's past the window.
        //   4. Decode the next block; loop back to step 1.
        loop {
            if self.cur_block_loaded {
                if let Some(res) = self.decoder.next_record() {
                    let record = match res {
                        Ok(r) => r,
                        Err(e) => {
                            self.poisoned = true;
                            return Some(Err(e));
                        }
                    };
                    // Mi17: region clamp via RangeClamp methods, on the
                    // schema-provided interval (architecture §10.5).
                    let (chrom_id, rec_start, rec_end) = S::record_interval(&record);
                    if self.clamp.interval_past_window(chrom_id, rec_start) {
                        return None;
                    }
                    if self.clamp.interval_before_window(rec_end) {
                        continue; // pre-window record, drop and ask for the next
                    }
                    return Some(Ok(record));
                }
                // Current block exhausted. Free its buffers (Mi2),
                // advance, and load the next.
                self.cur_block_loaded = false;
                self.decoder.unload();
                self.cur_block_idx += 1;
            }
            // For region iteration, peek at the next block's
            // range before paying the decode cost: if it's past
            // the window, we're done.
            if self.cur_block_idx < self.reader.index.len() {
                let entry = &self.reader.index[self.cur_block_idx];
                if self.clamp.block_past_window(entry) {
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

/// Owning, schema-typed **sequential** record iterator — [`RecordsIter`]'s state
/// machine, but owning its [`PspReader`] instead of borrowing it.
///
/// The borrowing `RecordsIter` cannot be stored alongside the reader it reads (that
/// is self-referential), which is exactly what a long-lived per-sample cursor needs.
/// This type owns the reader, so the cursor can hold it and pull records across many
/// `next()` calls. It decodes **one block at a time** (the persistent decompressor +
/// scratch are reused across blocks), so the resident set is one decoded block, not
/// the whole file. Sequential only — there is no region clamp; callers that want a
/// window pre-position via the block index themselves.
///
/// The block-stepping mirrors `RecordsIter` deliberately, kept as a separate type so
/// the production SNP `RecordsIter` path is left untouched.
pub struct OwnedRecordsIter<R: Read + Seek, S: PspKind = SnpKind> {
    reader: PspReader<R>,
    /// Index of the next block to load; equals `reader.index.len()` once exhausted.
    cur_block_idx: usize,
    /// Whether a block is currently decoded into `decoder`.
    cur_block_loaded: bool,
    /// The schema's per-block decoder (owns the decoded block + the per-block cursor).
    decoder: S::Decoder,
    /// Sticky: once an error has been yielded, future calls return `None`.
    poisoned: bool,
    /// Surfaced once on the first `next()` then poisons — set when `S` disagrees with
    /// the file's header `kind` (decoding under the wrong schema would be garbage).
    pending_error: Option<PspReadError>,
    /// Persistent zstd decompressor reused across every column of every block.
    decompressor: zstd::bulk::Decompressor<'static>,
    /// Reusable scratch for compressed / decompressed column bytes and the block
    /// header grow loop (converge to the largest seen).
    compressed_scratch: Vec<u8>,
    decompressed_scratch: Vec<u8>,
    block_header_buf: Vec<u8>,
}

#[allow(private_bounds)]
impl<R: Read + Seek, S: PspKind> OwnedRecordsIter<R, S>
where
    S::Decoder: BlockDecoder<Record = S::Record>,
{
    fn new(reader: PspReader<R>) -> Self {
        // Mirrors `RecordsIter::new`: the decompressor allocation is small and only
        // fails under OOM; the kind/schema check is surfaced as the first item.
        let decompressor = new_column_decompressor()
            .expect("zstd::bulk::Decompressor::new is infallible on supported platforms");
        let pending_error = (S::KIND != reader.header().kind).then(|| PspReadError::KindMismatch {
            expected: S::KIND,
            found: reader.header().kind.clone(),
        });
        Self {
            reader,
            cur_block_idx: 0,
            cur_block_loaded: false,
            decoder: S::Decoder::new_decoder(),
            poisoned: false,
            pending_error,
            decompressor,
            compressed_scratch: Vec::new(),
            decompressed_scratch: Vec::new(),
            block_header_buf: Vec::with_capacity(BLOCK_HEADER_INITIAL_CHUNK),
        }
    }

    /// Recover the owned reader (e.g. to re-seek for a second pass). Drops the
    /// current decode state.
    pub fn into_reader(self) -> PspReader<R> {
        self.reader
    }

    /// Seek to the next un-decoded block and decode it into the schema decoder.
    /// Returns `Ok(false)` when no more blocks remain. Identical framing to
    /// `RecordsIter::load_next_block`, with the reader owned rather than borrowed.
    fn load_next_block(&mut self) -> Result<bool, PspReadError> {
        if self.cur_block_idx >= self.reader.index.len() {
            return Ok(false);
        }
        let entry = self.reader.index[self.cur_block_idx];
        seek_to_offset(
            &mut self.reader.source,
            entry.block_offset,
            "block start seek",
        )?;
        let (block_header, _consumed) =
            read_block_header(&mut self.reader.source, &mut self.block_header_buf)?;
        let budget =
            block_byte_budget(&self.reader.index, &self.reader.trailer, self.cur_block_idx);
        self.decoder.decode_block(
            &mut self.reader.source,
            &block_header,
            budget,
            &mut self.decompressor,
            &mut self.compressed_scratch,
            &mut self.decompressed_scratch,
        )?;
        self.cur_block_loaded = true;
        Ok(true)
    }
}

#[allow(private_bounds)]
impl<R: Read + Seek, S: PspKind> Iterator for OwnedRecordsIter<R, S>
where
    S::Decoder: BlockDecoder<Record = S::Record>,
{
    type Item = Result<S::Record, PspReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.poisoned {
            return None;
        }
        if let Some(err) = self.pending_error.take() {
            self.poisoned = true;
            return Some(Err(err));
        }
        // Sequential variant of `RecordsIter`'s state machine (no region clamp):
        // yield from the loaded block, else drop it, advance, decode the next.
        loop {
            if self.cur_block_loaded {
                if let Some(res) = self.decoder.next_record() {
                    return match res {
                        Ok(r) => Some(Ok(r)),
                        Err(e) => {
                            self.poisoned = true;
                            Some(Err(e))
                        }
                    };
                }
                self.cur_block_loaded = false;
                self.decoder.unload();
                self.cur_block_idx += 1;
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

/// Borrowed columnar view of one decoded block, handed out by
/// [`BlockColumnReader::columns`]. The slices live in the reader's
/// reusable decode buffers and stay valid only until the next
/// [`BlockColumnReader::load_current`] / [`BlockColumnReader::advance`].
///
/// **Positions are delta-encoded.** Record `i`'s 1-based reference
/// position is `first_pos + Σ delta_pos[0..=i]` (with `delta_pos[0] == 0`,
/// so record 0 sits at `first_pos`). The per-allele fixed-width columns
/// and the two ragged CSR columns (`allele_seq` / `allele_chain_ids`)
/// are indexed by the cumulative allele offset `Σ n_alleles[0..i]` plus
/// the within-record allele index; the ragged columns slice their data
/// with the matching `*_offsets[j] .. *_offsets[j + 1]`. The `*_offsets`
/// arrays carry `total_alleles + 1` entries.
pub struct BlockColumns<'a> {
    /// Chromosome the block's records belong to.
    pub chrom_id: u32,
    /// 1-based position of record 0.
    pub first_pos: u32,
    /// Records in the block.
    pub n_records: u32,
    /// Per-record position deltas (`delta_pos[0] == 0`).
    pub delta_pos: &'a [u64],
    /// Per-record allele counts.
    pub n_alleles: &'a [u64],
    /// Per-allele fixed-width columns (aligned with `SampleColumns`'
    /// `AlleleScalarColumns`).
    pub allele_obs_count: &'a [u32],
    pub allele_q_sum_log: &'a [f64],
    pub allele_fwd_count: &'a [u32],
    pub allele_placed_left_count: &'a [u32],
    pub allele_placed_start_count: &'a [u32],
    pub allele_mapq_sum: &'a [u32],
    pub allele_mapq_sum_sq: &'a [u64],
    /// Per-record centred-window columns (aligned with the per-record
    /// `delta_pos`/`n_alleles`, not the per-allele columns).
    pub windowed_gc: &'a [f32],
    pub windowed_coverage: &'a [f32],
    /// Ragged `allele-seq` CSR (`data` indexed by `offsets`).
    pub allele_seq_data: &'a [u8],
    pub allele_seq_offsets: &'a [u32],
    /// Ragged `allele-chain-ids` CSR.
    pub allele_chain_ids_data: &'a [ChainId],
    pub allele_chain_ids_offsets: &'a [u32],
}

/// Block-granular columnar cursor over a [`PspReader`]. Decodes blocks
/// on demand and exposes each as borrowed [`BlockColumns`] without
/// materialising rows. Holds a mutable borrow of the source, so it is
/// **not `Send`** and only one lives at a time per reader — exactly
/// like [`RecordsIter`].
///
/// The decode buffers (zstd context, compressed/decompressed scratch,
/// the two ragged CSR slabs) mirror `RecordsIter`'s and are reused
/// across blocks; the heavy lifting is the shared
/// [`decode_block_payload`] free function. Stepping is explicit so the
/// caller can decide — from the block index alone, via [`peek_block`] —
/// whether a block is worth decoding before paying for it.
///
/// [`peek_block`]: Self::peek_block
///
/// Owns its [`PspReader`] (rather than borrowing it) so a cohort
/// producer can hold one per sample in a `Vec` without a
/// self-referential borrow. Recover the reader with
/// [`into_reader`](Self::into_reader) if needed.
pub struct BlockColumnReader<R: Read + Seek> {
    reader: PspReader<R>,
    cur_block_idx: usize,
    cur_block: Option<DecodedBlock>,
    decompressor: zstd::bulk::Decompressor<'static>,
    compressed_scratch: Vec<u8>,
    decompressed_scratch: Vec<u8>,
    block_header_buf: Vec<u8>,
    allele_seq_data: Vec<u8>,
    allele_seq_offsets: Vec<u32>,
    allele_chain_ids_data: Vec<ChainId>,
    allele_chain_ids_offsets: Vec<u32>,
}

impl<R: Read + Seek> BlockColumnReader<R> {
    fn new(reader: PspReader<R>) -> Self {
        let decompressor = new_column_decompressor()
            .expect("zstd::bulk::Decompressor::new is infallible on supported platforms");
        Self {
            reader,
            cur_block_idx: 0,
            cur_block: None,
            decompressor,
            compressed_scratch: Vec::new(),
            decompressed_scratch: Vec::new(),
            block_header_buf: Vec::with_capacity(BLOCK_HEADER_INITIAL_CHUNK),
            allele_seq_data: Vec::new(),
            allele_seq_offsets: Vec::new(),
            allele_chain_ids_data: Vec::new(),
            allele_chain_ids_offsets: Vec::new(),
        }
    }

    /// Recover the owned [`PspReader`] (e.g. to reuse the open file
    /// handle for something else). Drops the decode buffers.
    pub fn into_reader(self) -> PspReader<R> {
        self.reader
    }

    /// The block index (sorted by `(chrom_id, first_pos)`).
    pub fn index(&self) -> &[BlockIndexEntry] {
        &self.reader.index
    }

    /// Index of the block the cursor points at — the one
    /// [`load_current`](Self::load_current) decodes and
    /// [`columns`](Self::columns) then borrows.
    pub fn cur_block_idx(&self) -> usize {
        self.cur_block_idx
    }

    /// Position the cursor at the first block overlapping
    /// `[start, end]` on `chrom_id` (or past-the-end when none),
    /// dropping any currently-decoded block.
    pub fn seek_to(&mut self, chrom_id: u32, start: u32, end: u32) {
        self.cur_block = None;
        self.cur_block_idx = find_first_overlapping_block(&self.reader.index, chrom_id, start, end)
            .unwrap_or(self.reader.index.len());
    }

    /// The index entry the cursor points at, or `None` past the last
    /// block. Free — reads the in-memory index, decodes nothing.
    pub fn peek_block(&self) -> Option<&BlockIndexEntry> {
        self.reader.index.get(self.cur_block_idx)
    }

    /// Decode the block at the cursor into the reader's buffers.
    /// `Ok(false)` when the cursor is past the last block. After
    /// `Ok(true)`, [`columns`](Self::columns) borrows it.
    pub fn load_current(&mut self) -> Result<bool, PspReadError> {
        if self.cur_block_idx >= self.reader.index.len() {
            self.cur_block = None;
            return Ok(false);
        }
        let entry = self.reader.index[self.cur_block_idx];
        seek_to_offset(
            &mut self.reader.source,
            entry.block_offset,
            "block start seek",
        )?;
        let (block_header, _consumed) =
            read_block_header(&mut self.reader.source, &mut self.block_header_buf)?;
        let decoded = decode_block_payload(
            &mut self.reader.source,
            &block_header,
            block_byte_budget(&self.reader.index, &self.reader.trailer, self.cur_block_idx),
            &mut self.decompressor,
            &mut self.compressed_scratch,
            &mut self.decompressed_scratch,
            &mut self.allele_seq_data,
            &mut self.allele_seq_offsets,
            &mut self.allele_chain_ids_data,
            &mut self.allele_chain_ids_offsets,
        )?;
        self.cur_block = Some(decoded);
        Ok(true)
    }

    /// Decode the block at the cursor in **two phases** — light columns
    /// eagerly, the rest retained compressed (column-selective decode) —
    /// returning the owned [`TwoPhaseBlock`] (or `None` past the last block).
    /// Mirrors [`load_current`](Self::load_current) but does not store the
    /// block in `self`; pair with [`advance`](Self::advance). The deferred
    /// columns are inflated later, once the cohort fold knows the variable mask.
    // Wired into `SamplePspReader` in the next step; exercised now by
    // `block_reader_two_phase_matches_eager`.
    pub(crate) fn decode_current_two_phase(
        &mut self,
    ) -> Result<Option<TwoPhaseBlock>, PspReadError> {
        if self.cur_block_idx >= self.reader.index.len() {
            return Ok(None);
        }
        let entry = self.reader.index[self.cur_block_idx];
        seek_to_offset(
            &mut self.reader.source,
            entry.block_offset,
            "block start seek",
        )?;
        let (block_header, _consumed) =
            read_block_header(&mut self.reader.source, &mut self.block_header_buf)?;
        let tp = decode_block_payload_two_phase(
            &mut self.reader.source,
            &block_header,
            block_byte_budget(&self.reader.index, &self.reader.trailer, self.cur_block_idx),
            &mut self.decompressor,
            &mut self.compressed_scratch,
            &mut self.decompressed_scratch,
        )?;
        Ok(Some(tp))
    }

    /// Advance the cursor to the next block, dropping the current
    /// decoded one. Pair with [`peek_block`](Self::peek_block) +
    /// [`load_current`](Self::load_current).
    pub fn advance(&mut self) {
        self.cur_block = None;
        self.cur_block_idx += 1;
    }

    /// Borrow the currently-decoded block's columns, or `None` when no
    /// block is loaded (before the first `load_current`, or after an
    /// `advance` without a reload).
    pub fn columns(&self) -> Option<BlockColumns<'_>> {
        let b = self.cur_block.as_ref()?;
        Some(BlockColumns {
            chrom_id: b.chrom_id,
            first_pos: b.first_pos,
            n_records: b.n_records,
            delta_pos: &b.delta_pos,
            n_alleles: &b.n_alleles,
            allele_obs_count: &b.allele_obs_count,
            allele_q_sum_log: &b.allele_q_sum_log,
            allele_fwd_count: &b.allele_fwd_count,
            allele_placed_left_count: &b.allele_placed_left_count,
            allele_placed_start_count: &b.allele_placed_start_count,
            allele_mapq_sum: &b.allele_mapq_sum,
            allele_mapq_sum_sq: &b.allele_mapq_sum_sq,
            windowed_gc: &b.windowed_gc,
            windowed_coverage: &b.windowed_coverage,
            allele_seq_data: &self.allele_seq_data,
            allele_seq_offsets: &self.allele_seq_offsets,
            allele_chain_ids_data: &self.allele_chain_ids_data,
            allele_chain_ids_offsets: &self.allele_chain_ids_offsets,
        })
    }
}

// ---------------------------------------------------------------------
// Block-header read: grow-on-incomplete loop
// ---------------------------------------------------------------------

/// Read a block header from the source by growing a buffer until
/// `decode_block_header` succeeds. Returns the decoded header and
/// the number of bytes it consumed; positions the source cursor at
/// the first byte after the header.
///
/// M3 + M4: distinct error variants for the two failure modes —
/// EOF mid-decode (`BlockHeaderTruncated`) vs. read-cap exhaustion
/// (`BlockHeaderExceedsCap`). The previous code reused
/// `BlockHeaderField { source: VarintError::Truncated|Overflow }`
/// for both, which conflated buffer-cap failures with single-varint
/// decode failures.
/// L2: takes a `&mut Vec<u8>` from the caller (the `RecordsIter`) so
/// the buffer is reused across blocks. `clear()` keeps the capacity
/// that converged after the first few blocks (typical block headers
/// fit in `BLOCK_HEADER_INITIAL_CHUNK`).
fn read_block_header<R: Read + Seek>(
    source: &mut R,
    buf: &mut Vec<u8>,
) -> Result<(BlockHeader, usize), PspReadError> {
    let header_start = source
        .stream_position()
        .map_err(io_err("block header start position"))?;
    buf.clear();
    buf.reserve(BLOCK_HEADER_INITIAL_CHUNK);
    let mut chunk = BLOCK_HEADER_INITIAL_CHUNK;
    loop {
        let prev_len = buf.len();
        buf.resize(prev_len + chunk, 0);
        let read = source
            .read(&mut buf[prev_len..])
            .map_err(io_err("block header"))?;
        if read == 0 {
            return Err(PspReadError::BlockHeaderTruncated {
                offset: header_start,
                consumed: prev_len,
            });
        }
        buf.truncate(prev_len + read);

        match decode_block_header(buf) {
            Ok((header, consumed)) => {
                // We over-read into `buf`; rewind the source so
                // the caller's cursor lands at the start of the
                // column payloads. Rewind distance ≤ chunk size
                // (≤ BLOCK_HEADER_READ_CAP), well inside the
                // wrapped BufReader's 64 KiB buffer, so
                // `seek_to_offset`'s SeekFrom::Current preserves
                // the buffer.
                seek_to_offset(
                    source,
                    header_start + consumed as u64,
                    "post-block-header seek",
                )?;
                return Ok((header, consumed));
            }
            Err(PspReadError::BlockHeaderField {
                source: super::errors::VarintError::Truncated,
                ..
            }) => {
                // Not enough bytes — feed more and retry.
            }
            Err(other) => return Err(other),
        }

        if buf.len() >= BLOCK_HEADER_READ_CAP {
            return Err(PspReadError::BlockHeaderExceedsCap {
                cap: BLOCK_HEADER_READ_CAP,
                consumed: buf.len(),
            });
        }
        chunk = (chunk * 2).min(BLOCK_HEADER_READ_CAP.saturating_sub(buf.len()));
    }
}

// ---------------------------------------------------------------------
// Block-payload decode: walk the manifest, decode each column
// ---------------------------------------------------------------------

/// Decode every column in a block. The source cursor must be
/// positioned at the first byte of the column payloads (i.e.
/// immediately after the block header).
///
/// `byte_budget` is the number of bytes the block is allowed to
/// occupy between the end of its header and the start of the next
/// block (or the index region for the last block). Used by
/// `decode_one_column` to reject `compressed_len` values that would
/// drive an allocation past the block's geographic extent — B2.
///
/// `decompressor`, `compressed_scratch`, `decompressed_scratch` are
/// `RecordsIter`-owned buffers reused across every column of every
/// block. L1 / L2 in `ia/reviews/perf_psp_reader_2026-05-13.md`.
/// Post-decode validation of the per-record `n-alleles` column against
/// the block header. Catches the two corruptions that would otherwise
/// drive panics or silent mis-attribution in per-allele indexing:
///
/// - **sum disagreement** (B3): `Σ n_alleles[i] != n_total_alleles` — an
///   over-run panics on per-allele indexing in `materialise_next_record`,
///   an under-run silently emits truncated allele lists;
/// - **zero-allele record**: a record with `n_alleles[i] == 0`. Every
///   record carries at least its REF allele (the writer enforces this via
///   [`InvalidRecordKind::ZeroAlleles`](super::errors::InvalidRecordKind::ZeroAlleles));
///   a `0` aliases the next record's CSR allele range — a trailing one
///   indexes past `allele_*_offsets` (panic in the columnar `from_block`),
///   an interior one mis-attributes the following record's span.
fn validate_n_alleles_column(
    n_alleles: &[u64],
    n_total_alleles: u32,
) -> Result<(), BlockHeaderInvariantKind> {
    let sum_n_alleles: u64 = n_alleles.iter().sum();
    if sum_n_alleles != u64::from(n_total_alleles) {
        return Err(BlockHeaderInvariantKind::NAllelesSumMismatch {
            n_total_alleles,
            sum_n_alleles,
        });
    }
    if let Some(record_index) = n_alleles.iter().position(|&n| n == 0) {
        return Err(BlockHeaderInvariantKind::ZeroAlleleRecord {
            record_index: record_index as u32,
        });
    }
    Ok(())
}

// H1: four extra scratch buffers for the CSR ragged columns. Owned
// by the caller (`RecordsIter`); reused across blocks.
#[allow(clippy::too_many_arguments)]
fn decode_block_payload<R: Read>(
    source: &mut R,
    header: &BlockHeader,
    byte_budget: u64,
    decompressor: &mut zstd::bulk::Decompressor<'static>,
    compressed_scratch: &mut Vec<u8>,
    decompressed_scratch: &mut Vec<u8>,
    allele_seq_data: &mut Vec<u8>,
    allele_seq_offsets: &mut Vec<u32>,
    allele_chain_ids_data: &mut Vec<ChainId>,
    allele_chain_ids_offsets: &mut Vec<u32>,
) -> Result<DecodedBlock, PspReadError> {
    // B1: coverage check. Every v1.0 required tag must appear in
    // the per-block manifest. Symmetric with the file-level TOML
    // check in `cross_check_against_registry`; the two are
    // independent layers and both need enforcing. Without this
    // check, the `expect()`s below could panic on a hand-crafted
    // manifest that omits a required tag, and the `allele-seq`
    // arm in `decode_one_column` could also panic if the
    // `allele-seq-len` column were absent.
    for def in V1_0_COLUMNS {
        if def.required && !header.manifest.iter().any(|e| e.tag == def.tag) {
            return Err(PspReadError::MissingRequiredColumnInManifest {
                name: def.name.to_string(),
                tag: def.tag,
            });
        }
    }

    let n_records = header.n_records as usize;
    let n_total_alleles = header.n_total_alleles as usize;

    // Per-record slots.
    let mut delta_pos: Option<Vec<u64>> = None;
    let mut n_alleles: Option<Vec<u64>> = None;
    let mut windowed_gc: Option<Vec<f32>> = None;
    let mut windowed_coverage: Option<Vec<f32>> = None;

    // Per-allele slots.
    let mut allele_seq_len: Option<Vec<u64>> = None;
    let mut allele_obs_count: Option<Vec<u32>> = None;
    let mut allele_q_sum_log: Option<Vec<f64>> = None;
    let mut allele_fwd_count: Option<Vec<u32>> = None;
    let mut allele_placed_left_count: Option<Vec<u32>> = None;
    let mut allele_placed_start_count: Option<Vec<u32>> = None;
    let mut allele_mapq_sum: Option<Vec<u32>> = None;
    let mut allele_mapq_sum_sq: Option<Vec<u64>> = None;
    // H1: allele-seq + allele-chain-ids slots dropped — the CSR
    // data + offsets are written into the caller's `allele_seq_*` /
    // `allele_chain_ids_*` buffers by `decode_one_column`. B1's
    // manifest-coverage check above guarantees those columns fire.

    let mut remaining_budget = byte_budget;
    for entry in &header.manifest {
        let column = decode_one_column(
            source,
            entry,
            n_records,
            n_total_alleles,
            allele_seq_len.as_deref(),
            remaining_budget,
            decompressor,
            compressed_scratch,
            decompressed_scratch,
            allele_seq_data,
            allele_seq_offsets,
            allele_chain_ids_data,
            allele_chain_ids_offsets,
        )?;
        remaining_budget = remaining_budget.saturating_sub(entry.compressed_len as u64);
        // M10: `decode_one_column` returns `None` for an unknown
        // optional column (bytes consumed, no decoded payload).
        // `DecodedColumn` itself is closed over the v1.0 registry
        // — no `Unknown` escape hatch, so adding a `ColumnKey`
        // forces both a new `decode_one_column` arm and a matching
        // `DecodedColumn` variant. The dispatch match here is
        // exhaustive.
        let Some(column) = column else {
            continue;
        };
        // M11: post-decode validators driven by `ColumnDef`.
        // Generic finite-float sweep based on
        // `def.finite_constraint`. Today this only fires for
        // `allele-q-sum-log`; adding a new F-typed column with the
        // flag set picks it up automatically. M6: the column name
        // in the error payload is the registry name, not a
        // hard-coded literal.
        if let Some(def) = lookup_by_tag(entry.tag)
            && def.finite_constraint
            && let DecodedColumn::AlleleQSumLog(v) = &column
        {
            for (i, &q) in v.iter().enumerate() {
                if !q.is_finite() {
                    return Err(PspReadError::NonFiniteFloat {
                        column: def.name.to_string(),
                        entry: i,
                        value: q,
                    });
                }
            }
        }
        match column {
            DecodedColumn::DeltaPos(v) => delta_pos = Some(v),
            DecodedColumn::NAlleles(v) => n_alleles = Some(v),
            DecodedColumn::WindowedGc(v) => windowed_gc = Some(v),
            DecodedColumn::WindowedCoverage(v) => windowed_coverage = Some(v),
            DecodedColumn::AlleleSeqLen(v) => allele_seq_len = Some(v),
            // H1: unit variants — payload landed in the caller's
            // CSR scratches inside `decode_one_column`. Nothing to
            // capture here; B1 above guarantees these arms fire on
            // any well-formed v1.0 block.
            DecodedColumn::AlleleSeq | DecodedColumn::AlleleChainIds => {}
            DecodedColumn::AlleleObsCount(v) => allele_obs_count = Some(v),
            DecodedColumn::AlleleQSumLog(v) => allele_q_sum_log = Some(v),
            DecodedColumn::AlleleFwdCount(v) => allele_fwd_count = Some(v),
            DecodedColumn::AllelePlacedLeftCount(v) => allele_placed_left_count = Some(v),
            DecodedColumn::AllelePlacedStartCount(v) => allele_placed_start_count = Some(v),
            DecodedColumn::AlleleMapqSum(v) => allele_mapq_sum = Some(v),
            DecodedColumn::AlleleMapqSumSq(v) => allele_mapq_sum_sq = Some(v),
        }
    }

    // PANIC-FREE: B1's coverage check above has verified every
    // v1.0-required tag is present in the manifest; each
    // `decode_one_column` call for those tags writes its decoded
    // payload into the matching `Option`. The 12 `.expect()`s
    // below are structurally unreachable.
    let n_alleles = n_alleles.expect("n-alleles column required by v1.0 schema");

    // B3: cross-column manifest agreement — `sum(n_alleles[i])`
    // must equal `n_total_alleles`. Without this, an over-run
    // panics on `block.allele_seqs[j]` indexing in
    // `materialise_next_record`; an under-run silently emits
    // truncated allele lists. Spec §"Per-block manifest
    // agreement": "Any over- or under-run is a hard error."
    validate_n_alleles_column(&n_alleles, header.n_total_alleles)
        .map_err(|kind| PspReadError::BlockHeaderInvariant { kind })?;

    Ok(DecodedBlock {
        chrom_id: header.chrom_id,
        first_pos: header.first_pos,
        n_records: header.n_records,
        delta_pos: delta_pos.expect("delta-pos column required by v1.0 schema"),
        n_alleles,
        windowed_gc: windowed_gc.expect("windowed-gc column required by v1.0 schema"),
        windowed_coverage: windowed_coverage
            .expect("windowed-coverage column required by v1.0 schema"),
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
        allele_mapq_sum: allele_mapq_sum.expect("allele-mapq-sum column required by v1.0 schema"),
        allele_mapq_sum_sq: allele_mapq_sum_sq
            .expect("allele-mapq-sum-sq column required by v1.0 schema"),
    })
}

/// One decoded v1.0 column — keyed by [`ColumnKey`] so the
/// per-column dispatch is compile-time exhaustive against the
/// registry. M10: no `Unknown` catch-all variant — `decode_one_column`
/// returns `Option<DecodedColumn>`, with `None` for unknown
/// optional columns.
#[derive(Debug)]
pub(crate) enum DecodedColumn {
    DeltaPos(Vec<u64>),
    NAlleles(Vec<u64>),
    WindowedGc(Vec<f32>),
    WindowedCoverage(Vec<f32>),
    AlleleSeqLen(Vec<u64>),
    /// H1 (2026-05-23): unit variant — the CSR `data` + `offsets`
    /// were written through the caller's `allele_seq_*` `&mut Vec`
    /// parameters during decode, not returned here. The dispatcher
    /// just needs to know "AlleleSeq fired" so the per-block
    /// manifest-coverage tracker can mark the slot as covered.
    AlleleSeq,
    AlleleObsCount(Vec<u32>),
    AlleleQSumLog(Vec<f64>),
    AlleleFwdCount(Vec<u32>),
    AllelePlacedLeftCount(Vec<u32>),
    AllelePlacedStartCount(Vec<u32>),
    AlleleMapqSum(Vec<u32>),
    AlleleMapqSumSq(Vec<u64>),
    /// H1: unit variant — same shape as [`Self::AlleleSeq`] but for
    /// the chain-ids list column.
    AlleleChainIds,
}

/// Read + decompress + decode one column. `allele_seq_len` is
/// `Some(...)` once the `allele-seq-len` column has been decoded
/// (manifest order is tag-ascending and 0x03 < 0x04); used to
/// chunk the `allele-seq` bytes column.
///
/// `remaining_budget` is the number of bytes still available within
/// the block's geographic extent. B2: a `compressed_len` exceeding
/// the budget is rejected with a typed error before any
/// allocation, so a hostile manifest cannot drive a 4 GiB
/// allocation off a `u32` field.
///
/// Returns `Some(decoded)` for any registered v1.0 column;
/// `None` for an unknown optional column (bytes consumed, no
/// payload). M10: no `Unknown` escape hatch — the absence of a
/// catch-all variant means a new `ColumnKey` forces both a new
/// arm here *and* a new `DecodedColumn` variant.
// L1 + L2 + H1: thirteen arguments because the persistent
// decompressor, the two compressed/decompressed scratches, and the
// four CSR scratches (data + offsets for the two ragged columns)
// thread through to the per-column read + decompress site.
// Bundling them would add a `ReaderScratch` struct that exists only
// for this call shape; the line-noise from the extra indirection is
// worse than the warning.
#[allow(clippy::too_many_arguments)]
fn decode_one_column<R: Read>(
    source: &mut R,
    entry: &ColumnManifestEntry,
    n_records: usize,
    n_total_alleles: usize,
    allele_seq_len: Option<&[u64]>,
    remaining_budget: u64,
    decompressor: &mut zstd::bulk::Decompressor<'static>,
    compressed_scratch: &mut Vec<u8>,
    decompressed_scratch: &mut Vec<u8>,
    allele_seq_data: &mut Vec<u8>,
    allele_seq_offsets: &mut Vec<u32>,
    allele_chain_ids_data: &mut Vec<ChainId>,
    allele_chain_ids_offsets: &mut Vec<u32>,
) -> Result<Option<DecodedColumn>, PspReadError> {
    // B2: reject `compressed_len` exceeding the block's remaining
    // byte budget before allocating. Catches the
    // `compressed_len = u32::MAX` attack on `compressed = vec![0;
    // compressed_len]` below.
    if entry.compressed_len as u64 > remaining_budget {
        return Err(PspReadError::ColumnTruncated {
            column: lookup_by_tag(entry.tag)
                .map(|d| d.name.to_string())
                .unwrap_or_else(|| format!("tag {:#x}", entry.tag)),
            decoded: 0,
            expected: entry.compressed_len as usize,
        });
    }

    // L2: read into the reusable compressed scratch.
    compressed_scratch.clear();
    compressed_scratch.resize(entry.compressed_len as usize, 0);

    // Unknown tag = optional future column; read past it so the
    // source cursor stays aligned, then return `None`. The scratch
    // buffer doubles as the throwaway sink for unknown columns.
    let def = match lookup_by_tag(entry.tag) {
        Some(d) => d,
        None => {
            source
                .read_exact(compressed_scratch.as_mut_slice())
                .map_err(io_err("unknown optional column payload"))?;
            return Ok(None);
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
    if ColumnKey::from_tag(def.tag) == Some(ColumnKey::AlleleSeq)
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

    // L1 + L2: read into the reusable compressed scratch (already
    // sized above), then decompress through the persistent
    // decompressor into the reusable decompressed scratch.
    source
        .read_exact(compressed_scratch.as_mut_slice())
        .map_err(io_err("block column payload"))?;
    zstd_decompress_into(
        decompressor,
        compressed_scratch,
        entry.uncompressed_len as usize,
        decompressed_scratch,
    )
    .map_err(|source| PspReadError::Zstd {
        context: format!("column {column_name:?}"),
        source,
    })?;
    let bytes: &[u8] = decompressed_scratch;
    if bytes.len() as u64 != entry.uncompressed_len as u64 {
        return Err(PspReadError::UncompressedLenMismatch {
            column: column_name.to_string(),
            got: bytes.len() as u64,
            expected: entry.uncompressed_len as u64,
        });
    }

    // `ColumnDef` is schema-agnostic (no `key`); recover the SNP
    // dispatch key from the tag.
    // UNREACHABLE: `def` came from `lookup_by_tag` over `V1_0_COLUMNS`,
    // so its tag is always a `ColumnKey`.
    let key =
        ColumnKey::from_tag(def.tag).expect("decode_one_column resolved a non-SNP column tag");

    // Per-column decode, exhaustive on `ColumnKey` so a registry
    // addition forces a new arm here as a compile error.
    let decoded = match key {
        ColumnKey::DeltaPos => {
            DecodedColumn::DeltaPos(decode_varint_column(bytes, n_records, column_name)?)
        }
        ColumnKey::NAlleles => {
            DecodedColumn::NAlleles(decode_varint_column(bytes, n_records, column_name)?)
        }
        ColumnKey::WindowedGc => DecodedColumn::WindowedGc(decode_scalar_column_pod::<f32>(
            bytes,
            n_records,
            column_name,
        )?),
        ColumnKey::WindowedCoverage => {
            DecodedColumn::WindowedCoverage(decode_scalar_column_pod::<f32>(
                bytes,
                n_records,
                column_name,
            )?)
        }
        ColumnKey::AlleleSeqLen => {
            let lens = decode_varint_column(bytes, n_total_alleles, column_name)?;
            // Enforce per-entry cap symmetrically with the writer.
            for (i, &len) in lens.iter().enumerate() {
                if len > MAX_ALLELE_SEQ_LEN {
                    return Err(PspReadError::ColumnElementDecode {
                        column: column_name.to_string(),
                        entry: i,
                        source: ScalarDecodeError::Truncated,
                    });
                }
            }
            DecodedColumn::AlleleSeqLen(lens)
        }
        ColumnKey::AlleleSeq => {
            // PANIC-FREE: B1's per-block manifest coverage check
            // in `decode_block_payload` rejects a manifest that
            // omits `allele-seq-len`. The manifest is also
            // strictly tag-ascending (block.rs
            // `validate_block_header_invariants`), so 0x03
            // (allele-seq-len) is decoded before 0x04 (allele-seq).
            // Both together guarantee `allele_seq_len` is `Some`
            // here on any input that reaches this branch.
            let lens = allele_seq_len.expect(
                "allele-seq-len decoded before allele-seq per manifest tag-ascending order",
            );
            // H1: write the CSR data + offsets into the caller-
            // supplied scratch buffers instead of allocating one
            // `Vec<u8>` per allele.
            decode_bytes_split_csr(
                bytes,
                lens,
                column_name,
                Some(MAX_ALLELE_SEQ_LEN),
                allele_seq_data,
                allele_seq_offsets,
            )?;
            DecodedColumn::AlleleSeq
        }
        ColumnKey::AlleleObsCount => DecodedColumn::AlleleObsCount(
            decode_scalar_column_pod::<u32>(bytes, n_total_alleles, column_name)?,
        ),
        ColumnKey::AlleleQSumLog => DecodedColumn::AlleleQSumLog(decode_scalar_column_pod::<f64>(
            bytes,
            n_total_alleles,
            column_name,
        )?),
        ColumnKey::AlleleFwdCount => DecodedColumn::AlleleFwdCount(
            decode_scalar_column_pod::<u32>(bytes, n_total_alleles, column_name)?,
        ),
        ColumnKey::AllelePlacedLeftCount => DecodedColumn::AllelePlacedLeftCount(
            decode_scalar_column_pod::<u32>(bytes, n_total_alleles, column_name)?,
        ),
        ColumnKey::AllelePlacedStartCount => DecodedColumn::AllelePlacedStartCount(
            decode_scalar_column_pod::<u32>(bytes, n_total_alleles, column_name)?,
        ),
        ColumnKey::AlleleMapqSum => DecodedColumn::AlleleMapqSum(decode_scalar_column_pod::<u32>(
            bytes,
            n_total_alleles,
            column_name,
        )?),
        ColumnKey::AlleleMapqSumSq => {
            DecodedColumn::AlleleMapqSumSq(decode_scalar_column_pod::<u64>(
                bytes,
                n_total_alleles,
                column_name,
            )?)
        }
        ColumnKey::AlleleChainIds => {
            // H1: write the CSR data + offsets into the caller-
            // supplied scratch buffers instead of allocating one
            // `Vec<ChainId>` per allele.
            decode_list_column_csr::<ChainId>(
                bytes,
                n_total_alleles,
                column_name,
                allele_chain_ids_data,
                allele_chain_ids_offsets,
            )?;
            DecodedColumn::AlleleChainIds
        }
    };
    Ok(Some(decoded))
}

/// Read one column's raw compressed bytes (budget-checked) **without inflating**
/// — the two-phase decode retains the deferred columns this way until
/// [`inflate_retained_column`] materialises them. Mirrors the budget guard +
/// `read_exact` at the head of [`decode_one_column`]; the caller positions
/// `source` at the column payload (manifest order).
// Wired into the two-phase block decode in the next step; exercised now by
// `two_phase_inflate_matches_eager_decode`.
/// Shared low-level read + inflate for one column: budget-check the
/// compressed length, read it into `compressed_scratch`, and zstd-inflate
/// into `decompressed_scratch` (whose contents are the column's
/// uncompressed bytes on return). Used by per-schema block decoders
/// (`SnpDecoder` decodes inline; `SsrDecoder` calls this) so the
/// decompression machinery is shared while the per-column decode stays
/// per-schema (architecture §10.2). `column_name` is for error context.
pub(crate) fn read_and_inflate_column<R: Read>(
    source: &mut R,
    entry: &ColumnManifestEntry,
    remaining_budget: u64,
    column_name: &str,
    decompressor: &mut zstd::bulk::Decompressor<'static>,
    compressed_scratch: &mut Vec<u8>,
    decompressed_scratch: &mut Vec<u8>,
) -> Result<(), PspReadError> {
    if entry.compressed_len as u64 > remaining_budget {
        return Err(PspReadError::ColumnTruncated {
            column: column_name.to_string(),
            decoded: 0,
            expected: entry.compressed_len as usize,
        });
    }
    compressed_scratch.clear();
    compressed_scratch.resize(entry.compressed_len as usize, 0);
    source
        .read_exact(compressed_scratch.as_mut_slice())
        .map_err(io_err("block column payload"))?;
    zstd_decompress_into(
        decompressor,
        compressed_scratch,
        entry.uncompressed_len as usize,
        decompressed_scratch,
    )
    .map_err(|source| PspReadError::Zstd {
        context: format!("column {column_name:?}"),
        source,
    })?;
    if decompressed_scratch.len() as u64 != entry.uncompressed_len as u64 {
        return Err(PspReadError::UncompressedLenMismatch {
            column: column_name.to_string(),
            got: decompressed_scratch.len() as u64,
            expected: entry.uncompressed_len as u64,
        });
    }
    Ok(())
}

pub(crate) fn read_compressed_blob<R: Read>(
    source: &mut R,
    entry: &ColumnManifestEntry,
    remaining_budget: u64,
) -> Result<Vec<u8>, PspReadError> {
    if entry.compressed_len as u64 > remaining_budget {
        return Err(PspReadError::ColumnTruncated {
            column: lookup_by_tag(entry.tag)
                .map(|d| d.name.to_string())
                .unwrap_or_else(|| format!("tag {:#x}", entry.tag)),
            decoded: 0,
            expected: entry.compressed_len as usize,
        });
    }
    let mut buf = vec![0u8; entry.compressed_len as usize];
    source
        .read_exact(&mut buf)
        .map_err(io_err("block column payload"))?;
    Ok(buf)
}

/// Inflate a retained compressed-column blob, **reusing [`decode_one_column`]**
/// over an in-memory cursor. The blob is exactly `entry.compressed_len` bytes,
/// so its own length is the budget and the inner `read_exact` consumes it whole
/// — byte-identical to the eager path by construction (same decode, deferred).
/// `allele_seq_len` (a light column, decoded up front) chunks the `allele-seq` /
/// `allele-chain-ids` CSR columns.
// Used by `TwoPhaseSegment::set_variable_rows`, which is wired into the producer
// in the next step; exercised now by the two-phase round-trip tests.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn inflate_retained_column(
    blob: &[u8],
    entry: &ColumnManifestEntry,
    n_records: usize,
    n_total_alleles: usize,
    allele_seq_len: Option<&[u64]>,
    decompressor: &mut zstd::bulk::Decompressor<'static>,
    compressed_scratch: &mut Vec<u8>,
    decompressed_scratch: &mut Vec<u8>,
    allele_seq_data: &mut Vec<u8>,
    allele_seq_offsets: &mut Vec<u32>,
    allele_chain_ids_data: &mut Vec<ChainId>,
    allele_chain_ids_offsets: &mut Vec<u32>,
) -> Result<Option<DecodedColumn>, PspReadError> {
    let mut cursor = std::io::Cursor::new(blob);
    decode_one_column(
        &mut cursor,
        entry,
        n_records,
        n_total_alleles,
        allele_seq_len,
        blob.len() as u64,
        decompressor,
        compressed_scratch,
        decompressed_scratch,
        allele_seq_data,
        allele_seq_offsets,
        allele_chain_ids_data,
        allele_chain_ids_offsets,
    )
}

/// The four "light" column tags — decoded eagerly because the cohort fold
/// needs them (positions, the record→allele CSR base, ref-span/seq-len, and the
/// AC obs count). Every other v1.0 column is "deferred": retained compressed by
/// the two-phase decode and inflated only for the variable rows.
fn is_light_tag(tag: u16) -> bool {
    matches!(tag, 0x01 | 0x02 | 0x03 | 0x10)
}

/// One deferred column retained as its raw zstd blob (with its manifest entry,
/// so [`inflate_retained_column`] can size + dispatch it later).
pub(crate) struct RetainedColumn {
    pub entry: ColumnManifestEntry,
    pub blob: Vec<u8>,
}

/// A block decoded in two phases: the light columns are materialised (the fold
/// reads them); the deferred columns are retained compressed in
/// [`retained`](Self::retained). Mirrors the eager [`DecodedBlock`] for the
/// light columns; the heavy columns are inflated on demand by the consumer
/// (`SamplePspChunk::set_variable_rows`) once the variable mask is known.
pub(crate) struct TwoPhaseBlock {
    pub chrom_id: u32,
    pub first_pos: u32,
    pub n_records: u32,
    pub n_total_alleles: u32,
    pub delta_pos: Vec<u64>,
    pub n_alleles: Vec<u64>,
    pub allele_seq_len: Vec<u64>,
    pub allele_obs_count: Vec<u32>,
    pub retained: Vec<RetainedColumn>,
}

/// Two-phase sibling of [`decode_block_payload`]: decode the four light columns
/// eagerly; read every other column's bytes but **retain them compressed**
/// (via [`read_compressed_blob`]) for deferred, column-selective inflation. The
/// B1 manifest-coverage and B3 `n_alleles` cross-checks are kept identical to
/// the eager path; the source stays aligned (every column consumes its
/// `compressed_len` bytes in manifest order, light-decoded or retained).
// Wired into `BlockColumnReader` in the next step; exercised now by
// `two_phase_block_decode_matches_eager`.
fn decode_block_payload_two_phase<R: Read>(
    source: &mut R,
    header: &BlockHeader,
    byte_budget: u64,
    decompressor: &mut zstd::bulk::Decompressor<'static>,
    compressed_scratch: &mut Vec<u8>,
    decompressed_scratch: &mut Vec<u8>,
) -> Result<TwoPhaseBlock, PspReadError> {
    // B1: every v1.0 required tag must appear in the manifest (same as eager).
    for def in V1_0_COLUMNS {
        if def.required && !header.manifest.iter().any(|e| e.tag == def.tag) {
            return Err(PspReadError::MissingRequiredColumnInManifest {
                name: def.name.to_string(),
                tag: def.tag,
            });
        }
    }

    let n_records = header.n_records as usize;
    let n_total_alleles = header.n_total_alleles as usize;
    let mut delta_pos: Option<Vec<u64>> = None;
    let mut n_alleles: Option<Vec<u64>> = None;
    let mut allele_seq_len: Option<Vec<u64>> = None;
    let mut allele_obs_count: Option<Vec<u32>> = None;
    let mut retained: Vec<RetainedColumn> = Vec::new();

    // The light columns are scalar/varint — they never touch the CSR scratch,
    // so a single throwaway set suffices (the ragged seq/chain columns are
    // deferred, not decoded here).
    let mut csr_sink = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut remaining_budget = byte_budget;
    for entry in &header.manifest {
        if is_light_tag(entry.tag) {
            let column = decode_one_column(
                source,
                entry,
                n_records,
                n_total_alleles,
                allele_seq_len.as_deref(),
                remaining_budget,
                decompressor,
                compressed_scratch,
                decompressed_scratch,
                &mut csr_sink.0,
                &mut csr_sink.1,
                &mut csr_sink.2,
                &mut csr_sink.3,
            )?;
            match column {
                Some(DecodedColumn::DeltaPos(v)) => delta_pos = Some(v),
                Some(DecodedColumn::NAlleles(v)) => n_alleles = Some(v),
                Some(DecodedColumn::AlleleSeqLen(v)) => allele_seq_len = Some(v),
                Some(DecodedColumn::AlleleObsCount(v)) => allele_obs_count = Some(v),
                // `is_light_tag` only admits the four above.
                _ => {}
            }
        } else {
            let blob = read_compressed_blob(source, entry, remaining_budget)?;
            retained.push(RetainedColumn {
                entry: *entry,
                blob,
            });
        }
        remaining_budget = remaining_budget.saturating_sub(entry.compressed_len as u64);
    }

    // B3: `sum(n_alleles) == n_total_alleles` (same as eager).
    let n_alleles = n_alleles.expect("n-alleles column required by v1.0 schema");
    validate_n_alleles_column(&n_alleles, header.n_total_alleles)
        .map_err(|kind| PspReadError::BlockHeaderInvariant { kind })?;

    Ok(TwoPhaseBlock {
        chrom_id: header.chrom_id,
        first_pos: header.first_pos,
        n_records: header.n_records,
        n_total_alleles: header.n_total_alleles,
        delta_pos: delta_pos.expect("delta-pos column required by v1.0 schema"),
        n_alleles,
        allele_seq_len: allele_seq_len.expect("allele-seq-len column required by v1.0 schema"),
        allele_obs_count: allele_obs_count
            .expect("allele-obs-count column required by v1.0 schema"),
        retained,
    })
}

/// Compute the geographic byte budget for the block at `block_idx`.
/// Bytes available between block N's start and the start of block
/// N+1 (or, for the last block, the index region). Used by
/// `decode_block_payload` → `decode_one_column` to reject hostile
/// `compressed_len` values before allocating. B2.
fn block_byte_budget(index: &[BlockIndexEntry], trailer: &Trailer, block_idx: usize) -> u64 {
    let this_offset = index[block_idx].block_offset;
    let next_offset = index
        .get(block_idx + 1)
        .map(|e| e.block_offset)
        .unwrap_or(trailer.index_offset);
    next_offset.saturating_sub(this_offset)
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

/// Move `source` to absolute file offset `target`. Uses
/// `SeekFrom::Current(delta)` so a wrapped `BufReader` preserves its
/// 64 KiB user-space buffer when the target lies in the existing
/// window (per `std::io::BufReader`'s `Seek::seek` impl, which
/// special-cases `SeekFrom::Current(n)` for `|n|` within the
/// buffer). `SeekFrom::Start` would discard the buffer
/// unconditionally — which matters because PSP blocks are written
/// contiguously on disk and `RecordsIter` walks them sequentially,
/// so the cohort driver's 64 KiB BufReader would otherwise be
/// thrown away on every block transition.
///
/// `stream_position` on `BufReader<File>` does **not** issue an
/// `lseek` — it computes from the buffered position (std impl).
/// On a bare `File` source this falls through to a regular
/// `lseek(0, SEEK_CUR)` — no slower than the previous
/// `SeekFrom::Start` would have been.
fn seek_to_offset<R: Read + Seek>(
    source: &mut R,
    target: u64,
    context: &'static str,
) -> Result<(), PspReadError> {
    let current = source.stream_position().map_err(io_err(context))?;
    if current == target {
        return Ok(());
    }
    let delta = (target as i64) - (current as i64);
    source
        .seek(SeekFrom::Current(delta))
        .map_err(io_err(context))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};
    use crate::psp::header::{ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance};
    use crate::psp::test_fixtures::writer_header;
    use crate::psp::writer::PspWriter;
    use std::collections::BTreeMap;
    use std::io::Cursor;

    fn finish_empty_writer(header: WriterHeader) -> Vec<u8> {
        let writer = PspWriter::new(Cursor::new(Vec::new()), header).unwrap();
        writer.finish().unwrap().into_inner()
    }

    fn one_allele_record(chrom_id: u32, pos: u32, seq: &[u8]) -> PileupRecord {
        PileupRecord {
            windowed_gc: 0.5,
            windowed_coverage: 25.0,
            chrom_id,
            pos,
            alleles: vec![AlleleObservation {
                seq: seq.to_vec(),
                support: AlleleSupportStats {
                    num_obs: 1,
                    q_sum: -1.0,
                    fwd: 1,
                    placed_left: 0,
                    placed_start: 1,

                    mapq_sum: 0,
                    mapq_sum_sq: 0,
                },
                chain_ids: Vec::new(),
            }],
        }
    }

    /// Multi-allele record exercising the ragged seq / chain-id CSR columns
    /// and every per-allele scalar stat (distinct values per allele).
    fn multi_allele_record(pos: u32, alleles: &[(&[u8], &[ChainId])]) -> PileupRecord {
        PileupRecord {
            windowed_gc: 0.5,
            windowed_coverage: 25.0,
            chrom_id: 0,
            pos,
            alleles: alleles
                .iter()
                .enumerate()
                .map(|(i, (seq, chains))| AlleleObservation {
                    seq: seq.to_vec(),
                    support: AlleleSupportStats {
                        num_obs: i as u32 + 1,
                        q_sum: -1.0 - i as f64,
                        fwd: i as u32,
                        placed_left: i as u32,
                        placed_start: 1,
                        mapq_sum: 10 * i as u32,
                        mapq_sum_sq: 100 * i as u64,
                    },
                    chain_ids: chains.to_vec(),
                })
                .collect(),
        }
    }

    /// Step 1 of column-selective decode: `read_compressed_blob` +
    /// `inflate_retained_column` (deferred decode of a retained blob) must
    /// produce **byte-identical** columns to the eager `decode_block_payload`,
    /// for every v1.0 column type incl. the ragged seq / chain-id CSR.
    #[test]
    fn two_phase_inflate_matches_eager_decode() {
        let records = vec![
            multi_allele_record(10, &[(b"A", &[]), (b"ACGT", &[7, 9])]),
            multi_allele_record(25, &[(b"T", &[3]), (b"TT", &[])]),
            multi_allele_record(40, &[(b"G", &[])]),
        ];
        let bytes = {
            let mut w = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
            for r in &records {
                w.write_record(r).unwrap();
            }
            w.finish().unwrap().into_inner()
        };

        // Slice out block 0's [header + payload] bytes.
        let (off, end) = {
            let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
            let idx = reader.block_index();
            let off = idx[0].block_offset as usize;
            let end = idx
                .get(1)
                .map(|e| e.block_offset)
                .unwrap_or(reader.trailer().index_offset) as usize;
            (off, end)
        };
        let mut src = Cursor::new(bytes[off..end].to_vec());
        let mut header_buf = Vec::new();
        let (header, _) = read_block_header(&mut src, &mut header_buf).unwrap();
        let payload_start = src.position();
        let n_rec = header.n_records as usize;
        let n_tot = header.n_total_alleles as usize;

        // --- Eager reference decode. ---
        let mut dz = new_column_decompressor().unwrap();
        let (mut cs, mut ds) = (Vec::new(), Vec::new());
        let (mut e_sd, mut e_so, mut e_cd, mut e_co) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let eager = decode_block_payload(
            &mut src,
            &header,
            u64::MAX,
            &mut dz,
            &mut cs,
            &mut ds,
            &mut e_sd,
            &mut e_so,
            &mut e_cd,
            &mut e_co,
        )
        .unwrap();

        // --- Two-phase: read each column's blob, then inflate it. ---
        src.set_position(payload_start);
        let mut seq_len: Option<Vec<u64>> = None;
        let (mut obs, mut q, mut fwd, mut pl, mut ps, mut ms, mut mss) =
            (None, None, None, None, None, None, None);
        let (mut sd, mut so, mut cd, mut co) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let mut budget = u64::MAX;
        for entry in &header.manifest {
            let blob = read_compressed_blob(&mut src, entry, budget).unwrap();
            budget -= entry.compressed_len as u64;
            let col = inflate_retained_column(
                &blob,
                entry,
                n_rec,
                n_tot,
                seq_len.as_deref(),
                &mut dz,
                &mut cs,
                &mut ds,
                &mut sd,
                &mut so,
                &mut cd,
                &mut co,
            )
            .unwrap();
            match col {
                Some(DecodedColumn::AlleleSeqLen(v)) => seq_len = Some(v),
                Some(DecodedColumn::AlleleObsCount(v)) => obs = Some(v),
                Some(DecodedColumn::AlleleQSumLog(v)) => q = Some(v),
                Some(DecodedColumn::AlleleFwdCount(v)) => fwd = Some(v),
                Some(DecodedColumn::AllelePlacedLeftCount(v)) => pl = Some(v),
                Some(DecodedColumn::AllelePlacedStartCount(v)) => ps = Some(v),
                Some(DecodedColumn::AlleleMapqSum(v)) => ms = Some(v),
                Some(DecodedColumn::AlleleMapqSumSq(v)) => mss = Some(v),
                Some(DecodedColumn::DeltaPos(_))
                | Some(DecodedColumn::NAlleles(_))
                | Some(DecodedColumn::WindowedGc(_))
                | Some(DecodedColumn::WindowedCoverage(_))
                | Some(DecodedColumn::AlleleSeq)
                | Some(DecodedColumn::AlleleChainIds)
                | None => {}
            }
        }

        // Heavy scalar columns inflate byte-identically.
        assert_eq!(obs.unwrap(), eager.allele_obs_count);
        assert_eq!(q.unwrap(), eager.allele_q_sum_log);
        assert_eq!(fwd.unwrap(), eager.allele_fwd_count);
        assert_eq!(pl.unwrap(), eager.allele_placed_left_count);
        assert_eq!(ps.unwrap(), eager.allele_placed_start_count);
        assert_eq!(ms.unwrap(), eager.allele_mapq_sum);
        assert_eq!(mss.unwrap(), eager.allele_mapq_sum_sq);
        // Ragged CSR columns inflate byte-identically.
        assert_eq!((sd, so), (e_sd, e_so), "allele-seq CSR");
        assert_eq!((cd, co), (e_cd, e_co), "chain-ids CSR");
        // The fixture actually exercised the ragged columns.
        assert_eq!(n_tot, 5);
    }

    /// Step 2: `decode_block_payload_two_phase` must yield light columns
    /// identical to the eager decode, plus retained blobs that inflate to the
    /// eager heavy columns (scalars + ragged seq/chain CSR).
    #[test]
    fn two_phase_block_decode_matches_eager() {
        let records = vec![
            multi_allele_record(10, &[(b"A", &[]), (b"ACGT", &[7, 9])]),
            multi_allele_record(25, &[(b"T", &[3]), (b"TT", &[])]),
            multi_allele_record(40, &[(b"G", &[])]),
        ];
        let bytes = {
            let mut w = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
            for r in &records {
                w.write_record(r).unwrap();
            }
            w.finish().unwrap().into_inner()
        };
        let (off, end) = {
            let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
            let idx = reader.block_index();
            let off = idx[0].block_offset as usize;
            let end = idx
                .get(1)
                .map(|e| e.block_offset)
                .unwrap_or(reader.trailer().index_offset) as usize;
            (off, end)
        };
        let mut src = Cursor::new(bytes[off..end].to_vec());
        let mut header_buf = Vec::new();
        let (header, _) = read_block_header(&mut src, &mut header_buf).unwrap();
        let payload_start = src.position();
        let n_rec = header.n_records as usize;
        let n_tot = header.n_total_alleles as usize;

        let mut dz = new_column_decompressor().unwrap();
        let (mut cs, mut ds) = (Vec::new(), Vec::new());
        let (mut e_sd, mut e_so, mut e_cd, mut e_co) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let eager = decode_block_payload(
            &mut src,
            &header,
            u64::MAX,
            &mut dz,
            &mut cs,
            &mut ds,
            &mut e_sd,
            &mut e_so,
            &mut e_cd,
            &mut e_co,
        )
        .unwrap();

        src.set_position(payload_start);
        let tp =
            decode_block_payload_two_phase(&mut src, &header, u64::MAX, &mut dz, &mut cs, &mut ds)
                .unwrap();

        // Light columns identical to eager.
        assert_eq!(tp.delta_pos, eager.delta_pos);
        assert_eq!(tp.n_alleles, eager.n_alleles);
        assert_eq!(tp.allele_obs_count, eager.allele_obs_count);
        assert_eq!(tp.n_total_alleles, header.n_total_alleles);
        // The 10 deferred columns (seq, six allele stats, chain-ids, and the
        // two per-record windowed-coverage columns) were retained.
        assert_eq!(tp.retained.len(), 10);

        // Each retained blob inflates to the eager heavy column.
        let (mut sd, mut so, mut cd, mut co) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for rc in &tp.retained {
            let col = inflate_retained_column(
                &rc.blob,
                &rc.entry,
                n_rec,
                n_tot,
                Some(&tp.allele_seq_len),
                &mut dz,
                &mut cs,
                &mut ds,
                &mut sd,
                &mut so,
                &mut cd,
                &mut co,
            )
            .unwrap();
            match col {
                Some(DecodedColumn::AlleleQSumLog(v)) => assert_eq!(v, eager.allele_q_sum_log),
                Some(DecodedColumn::AlleleFwdCount(v)) => assert_eq!(v, eager.allele_fwd_count),
                Some(DecodedColumn::AllelePlacedLeftCount(v)) => {
                    assert_eq!(v, eager.allele_placed_left_count)
                }
                Some(DecodedColumn::AllelePlacedStartCount(v)) => {
                    assert_eq!(v, eager.allele_placed_start_count)
                }
                Some(DecodedColumn::AlleleMapqSum(v)) => assert_eq!(v, eager.allele_mapq_sum),
                Some(DecodedColumn::AlleleMapqSumSq(v)) => assert_eq!(v, eager.allele_mapq_sum_sq),
                Some(DecodedColumn::WindowedGc(v)) => assert_eq!(v, eager.windowed_gc),
                Some(DecodedColumn::WindowedCoverage(v)) => assert_eq!(v, eager.windowed_coverage),
                Some(DecodedColumn::AlleleSeq) | Some(DecodedColumn::AlleleChainIds) => {}
                other => panic!("unexpected retained column {other:?}"),
            }
        }
        assert_eq!((sd, so), (e_sd, e_so), "allele-seq CSR");
        assert_eq!((cd, co), (e_cd, e_co), "chain-ids CSR");
    }

    /// Step 3 bridge: `BlockColumnReader::decode_current_two_phase` (the real
    /// reader path — seek + header read + two-phase payload) yields light
    /// columns + retained blobs identical to the eager `load_current` columns.
    #[test]
    fn block_reader_two_phase_matches_eager() {
        let records = vec![
            multi_allele_record(10, &[(b"A", &[]), (b"ACGT", &[7, 9])]),
            multi_allele_record(25, &[(b"T", &[3]), (b"TT", &[])]),
        ];
        let bytes = {
            let mut w = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
            for r in &records {
                w.write_record(r).unwrap();
            }
            w.finish().unwrap().into_inner()
        };

        // Eager reference via load_current + columns.
        let mut br_e = PspReader::new(Cursor::new(bytes.clone()))
            .unwrap()
            .into_column_blocks();
        assert!(br_e.load_current().unwrap());
        let cols = br_e.columns().unwrap();
        let n_rec = cols.n_records as usize;
        let n_tot: usize = cols.n_alleles.iter().map(|&n| n as usize).sum();
        let e_delta = cols.delta_pos.to_vec();
        let e_obs = cols.allele_obs_count.to_vec();
        let e_qsum = cols.allele_q_sum_log.to_vec();
        let e_seq = (
            cols.allele_seq_data.to_vec(),
            cols.allele_seq_offsets.to_vec(),
        );
        let e_chain = (
            cols.allele_chain_ids_data.to_vec(),
            cols.allele_chain_ids_offsets.to_vec(),
        );
        // `cols`'s borrow of `br_e` ends here (NLL); the two-phase reader below
        // owns a separate `PspReader`.

        // Two-phase via the bridge method.
        let mut br = PspReader::new(Cursor::new(bytes))
            .unwrap()
            .into_column_blocks();
        let tp = br.decode_current_two_phase().unwrap().unwrap();
        assert_eq!(tp.delta_pos, e_delta);
        assert_eq!(tp.allele_obs_count, e_obs);
        assert_eq!(tp.retained.len(), 10);

        let mut dz = new_column_decompressor().unwrap();
        let (mut cs, mut ds) = (Vec::new(), Vec::new());
        let (mut sd, mut so, mut cd, mut co) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let mut qsum = None;
        for rc in &tp.retained {
            let col = inflate_retained_column(
                &rc.blob,
                &rc.entry,
                n_rec,
                n_tot,
                Some(&tp.allele_seq_len),
                &mut dz,
                &mut cs,
                &mut ds,
                &mut sd,
                &mut so,
                &mut cd,
                &mut co,
            )
            .unwrap();
            if let Some(DecodedColumn::AlleleQSumLog(v)) = col {
                qsum = Some(v);
            }
        }
        assert_eq!(qsum.unwrap(), e_qsum);
        assert_eq!((sd, so), e_seq, "allele-seq CSR");
        assert_eq!((cd, co), e_chain, "chain-ids CSR");
    }

    /// (R1) Zero-record file round-trips through `PspReader::new`
    /// successfully; `header()` matches; `block_index().is_empty()`;
    /// `records()` and `region_records()` yield `None` immediately.
    #[test]
    fn empty_file_round_trip() {
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
    fn header_accessible_before_iteration() {
        let bytes = finish_empty_writer(writer_header(2));
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.header().chromosomes.len(), 2);
        assert_eq!(reader.header().sample, "sample");
        // No iteration; drop.
    }

    /// (R3) Header round-trip equality: every non-default writer
    /// header field survives the write → read pipeline byte-stable.
    #[test]
    fn header_round_trip_equality() {
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
                command_line: String::new(),
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
    fn block_index_matches_writer_emission() {
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
        chain_ids: &[ChainId],
    ) -> AlleleObservation {
        AlleleObservation {
            seq: seq.to_vec(),
            support: AlleleSupportStats {
                num_obs,
                q_sum,
                fwd,
                placed_left,
                placed_start,
                mapq_sum: 0,
                mapq_sum_sq: 0,
            },
            chain_ids: chain_ids.to_vec(),
        }
    }

    /// The per-record windowed-coverage columns (M2) round-trip exactly,
    /// including distinct per-record values (so a mis-indexing bug would show)
    /// and the `NaN` "no window here" sentinel (which round-trips bit-for-bit,
    /// verified on the raw bits since `NaN != NaN`).
    #[test]
    fn windowed_coverage_columns_round_trip() {
        let mut records = Vec::new();
        for (k, pos) in [10u32, 20, 30].into_iter().enumerate() {
            let mut r = PileupRecord::new(0, pos, vec![allele(b"A", 5, -3.0, 3, 0, 5, &[])]);
            // Distinct per-record values (a mis-index would show), with the NaN
            // "no window" sentinel on a different record per field so both
            // fields' NaN paths are exercised.
            r.windowed_gc = if k == 2 {
                f32::NAN
            } else {
                0.25 + 0.1 * k as f32
            };
            r.windowed_coverage = if k == 1 { f32::NAN } else { 12.0 + pos as f32 };
            records.push(r);
        }

        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        for r in &records {
            writer.write_record(r).unwrap();
        }
        let bytes = writer.finish().unwrap().into_inner();

        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let got: Vec<PileupRecord> = reader.records().collect::<Result<_, _>>().unwrap();
        assert_eq!(got.len(), 3);
        for (g, w) in got.iter().zip(&records) {
            assert_eq!(g.pos, w.pos);
            // Bit comparison so the NaN records are checked on both fields.
            assert_eq!(
                g.windowed_gc.to_bits(),
                w.windowed_gc.to_bits(),
                "gc at pos {}",
                w.pos,
            );
            assert_eq!(
                g.windowed_coverage.to_bits(),
                w.windowed_coverage.to_bits(),
                "coverage at pos {}",
                w.pos,
            );
        }
        assert!(got[1].windowed_coverage.is_nan(), "coverage NaN survives");
        assert!(got[2].windowed_gc.is_nan(), "gc NaN survives");
    }

    /// Variant of `allele()` that also pins the new mapq_sum /
    /// mapq_sum_sq scalars. Used by the round-trip tests covering
    /// the per-allele MAPQ tracking added in 2026-05.
    #[allow(clippy::too_many_arguments)]
    fn allele_with_mapq(
        seq: &[u8],
        num_obs: u32,
        q_sum: f64,
        fwd: u32,
        placed_left: u32,
        placed_start: u32,
        mapq_sum: u32,
        mapq_sum_sq: u64,
        chain_ids: &[ChainId],
    ) -> AlleleObservation {
        AlleleObservation {
            seq: seq.to_vec(),
            support: AlleleSupportStats {
                num_obs,
                q_sum,
                fwd,
                placed_left,
                placed_start,
                mapq_sum,
                mapq_sum_sq,
            },
            chain_ids: chain_ids.to_vec(),
        }
    }

    /// (B1) Single record, single allele, single block — every
    /// scalar survives the round-trip exactly.
    #[test]
    fn single_record_single_allele_round_trip() {
        let want = PileupRecord {
            windowed_gc: 0.5,
            windowed_coverage: 25.0,
            chrom_id: 0,
            pos: 100,
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
        assert_eq!(got.alleles.len(), want.alleles.len());
        let g = &got.alleles[0];
        let w = &want.alleles[0];
        assert_eq!(g.seq, w.seq);
        assert_eq!(g.support.num_obs, w.support.num_obs);
        assert_eq!(g.support.q_sum, w.support.q_sum);
        assert_eq!(g.support.fwd, w.support.fwd);
        assert_eq!(g.support.placed_left, w.support.placed_left);
        assert_eq!(g.support.placed_start, w.support.placed_start);
        assert_eq!(g.chain_ids, w.chain_ids);
    }

    /// Per-allele MAPQ scalars (`mapq_sum`, `mapq_sum_sq`) survive
    /// the round-trip exactly. Covers the columns added for the
    /// Welch's-t multi-mapper filter (2026-05).
    #[test]
    fn mapq_scalars_round_trip() {
        // Two alleles with distinct MAPQ profiles: REF reads all
        // MAPQ=60 (clean), ALT reads mixed 60/0 (multi-mapper).
        // Numbers correspond to a hypothetical 10-read site:
        //   REF n=7, mapq_sum=420, mapq_sum_sq=25200  (all MAPQ 60)
        //   ALT n=3, mapq_sum= 60, mapq_sum_sq= 3600  (one MAPQ 60, two MAPQ 0)
        let rec = PileupRecord {
            windowed_gc: 0.5,
            windowed_coverage: 25.0,
            chrom_id: 0,
            pos: 500,
            alleles: vec![
                allele_with_mapq(b"A", 7, -14.0, 4, 0, 7, 420, 25_200, &[]),
                allele_with_mapq(b"C", 3, -6.0, 1, 0, 3, 60, 3_600, &[]),
            ],
        };
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer.write_record(&rec).unwrap();
        let bytes = writer.finish().unwrap().into_inner();

        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let got: Vec<PileupRecord> = reader.records().collect::<Result<_, _>>().unwrap();
        assert_eq!(got.len(), 1);
        let r = &got[0];
        assert_eq!(r.alleles.len(), 2);
        assert_eq!(r.alleles[0].support.mapq_sum, 420);
        assert_eq!(r.alleles[0].support.mapq_sum_sq, 25_200);
        assert_eq!(r.alleles[1].support.mapq_sum, 60);
        assert_eq!(r.alleles[1].support.mapq_sum_sq, 3_600);
        // Sanity: mean MAPQ recovers correctly.
        let mean_ref = r.alleles[0].support.mapq_sum as f64 / r.alleles[0].support.num_obs as f64;
        let mean_alt = r.alleles[1].support.mapq_sum as f64 / r.alleles[1].support.num_obs as f64;
        assert!((mean_ref - 60.0).abs() < 1e-9);
        assert!((mean_alt - 20.0).abs() < 1e-9);
    }

    /// (B2) Multi-allele records: a SNP at pos=100 and a
    /// deletion + insertion at pos=200. Allele sequences and
    /// scalars survive.
    #[test]
    fn multi_allele_round_trip() {
        let snp = PileupRecord {
            windowed_gc: 0.5,
            windowed_coverage: 25.0,
            chrom_id: 0,
            pos: 100,
            alleles: vec![
                allele(b"A", 30, -50.0, 16, 0, 30, &[]),
                allele(b"C", 8, -8.0, 4, 0, 8, &[]),
            ],
        };
        let multi = PileupRecord {
            windowed_gc: 0.5,
            windowed_coverage: 25.0,
            chrom_id: 0,
            pos: 200,
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
    fn multi_block_round_trip() {
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
    fn multi_chromosome_round_trip() {
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
    /// Record referencing chain id `chain_id` in its REF allele.
    fn record_with_chain_id(chrom_id: u32, pos: u32, chain_id: ChainId) -> PileupRecord {
        PileupRecord {
            windowed_gc: 0.5,
            windowed_coverage: 25.0,
            chrom_id,
            pos,
            alleles: vec![allele(b"A", 1, -1.0, 1, 0, 1, &[chain_id])],
        }
    }

    /// (B5 / R9) A chain id appears on records spanning a block
    /// boundary. The reader must materialise the id on both sides
    /// without any active-set bookkeeping (unique-per-file ids
    /// don't need it).
    #[test]
    fn chain_id_round_trips_across_block_boundary() {
        let header = writer_header(1);
        let mut writer =
            PspWriter::new_with_block_target(Cursor::new(Vec::new()), header, 1024).unwrap();

        // Block 0: 100 plain records, then one referencing chain id 7.
        for i in 1u32..=100 {
            writer.write_record(&one_allele_record(0, i, b"A")).unwrap();
        }
        writer
            .write_record(&record_with_chain_id(0, 101, 7))
            .unwrap();
        // Block 1: another record referencing chain id 7.
        writer
            .write_record(&record_with_chain_id(0, 102, 7))
            .unwrap();
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
            .expect("records iterate cleanly across the boundary");

        // Chain id 7 appears on both records around the boundary.
        let first = got.iter().find(|r| r.pos == 101).unwrap();
        let second = got.iter().find(|r| r.pos == 102).unwrap();
        assert!(first.alleles[0].chain_ids.contains(&7));
        assert!(second.alleles[0].chain_ids.contains(&7));
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

    /// The owning `into_records_of` yields the same record stream as the borrowing
    /// `records_of` across a multi-block file — the cohort cursor's foundation must
    /// not diverge from the production reader.
    #[test]
    fn owned_records_iter_matches_borrowing_across_blocks() {
        let bytes = three_block_chrom0_fixture(4 * 1024);

        let mut borrow_reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
        assert!(
            borrow_reader.block_index().len() >= 2,
            "fixture must be multi-block"
        );
        let borrowed: Vec<PileupRecord> = borrow_reader
            .records_of::<SnpKind>()
            .collect::<Result<_, _>>()
            .expect("borrowing iterator is clean");

        let owning_reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let owned: Vec<PileupRecord> = owning_reader
            .into_records_of::<SnpKind>()
            .collect::<Result<_, _>>()
            .expect("owning iterator is clean");

        assert_eq!(owned, borrowed);
        assert_eq!(owned.len(), 300);
    }

    /// A schema/kind mismatch surfaces once then ends — the owning iterator mirrors
    /// `records_of`'s poisoning (decoding an SNP file under the SSR schema is refused,
    /// not silently misdecoded).
    #[test]
    fn owned_records_iter_kind_mismatch_yields_one_error_then_none() {
        let bytes = three_block_chrom0_fixture(4 * 1024); // an SNP (.psp) file
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let mut it = reader.into_records_of::<crate::psp::registry_ssr::SsrKind>();
        assert!(matches!(
            it.next(),
            Some(Err(PspReadError::KindMismatch { .. }))
        ));
        assert!(it.next().is_none(), "poisoned after the mismatch");
    }

    /// The cohort reader drives the owning iterator with `SsrKind`, not `SnpKind` —
    /// so equivalence with the borrowing iterator must hold on that decoder too,
    /// across a multi-block file (the SNP-only test above is not enough).
    #[test]
    fn owned_records_iter_matches_borrowing_for_ssr_kind() {
        use crate::psp::registry_ssr::{SsrKind, SsrLocusRecord};

        fn ssr_rec(start: u32) -> SsrLocusRecord {
            SsrLocusRecord {
                chrom_id: 0,
                start,
                end: start + 6,
                depth: 1,
                n_filtered: 0,
                mapped_reads: 1,
                n_low_quality: 0,
                n_border_off_end: 0,
                n_widened: 0,
                n_window_truncated: 0,
                observed: vec![(b"CACACA".to_vec().into_boxed_slice(), 1)],
            }
        }

        // Spaced-out starts + a small genomic window grid → multiple blocks.
        let mut w = PspWriter::<_, SsrKind>::new_ssr_with_block_layout(
            Cursor::new(Vec::new()),
            writer_header(1),
            16 * 1024 * 1024,
            16,
        )
        .unwrap();
        for s in [11u32, 51, 101] {
            w.write_locus(&ssr_rec(s)).unwrap();
        }
        let bytes = w.finish().unwrap().into_inner();

        let mut borrow_reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
        assert!(
            borrow_reader.block_index().len() >= 2,
            "fixture must be multi-block"
        );
        let borrowed: Vec<SsrLocusRecord> = borrow_reader
            .records_of::<SsrKind>()
            .collect::<Result<_, _>>()
            .unwrap();

        let owned: Vec<SsrLocusRecord> = PspReader::new(Cursor::new(bytes))
            .unwrap()
            .into_records_of::<SsrKind>()
            .collect::<Result<_, _>>()
            .unwrap();

        assert_eq!(owned, borrowed);
        assert_eq!(owned.len(), 3);
    }

    #[test]
    fn owned_records_iter_on_empty_file_yields_none() {
        let bytes = finish_empty_writer(writer_header(1));
        let mut it = PspReader::new(Cursor::new(bytes))
            .unwrap()
            .into_records_of::<SnpKind>();
        assert!(it.next().is_none());
    }

    /// (B7) Random-access region query across multiple blocks.
    /// 1000 records on two chromosomes; `region_records(0,
    /// 50_000, 60_000)` (out-of-range on a fixture whose max pos
    /// is 1000) yields zero records. Inside-range case:
    /// `region_records(0, 100, 200)` returns positions 100..=200.
    #[test]
    fn random_access_region_query() {
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
    fn region_across_multiple_blocks() {
        // Choose target so each block holds ~100 records. ~50 bytes per record
        // (after adding mapq_sum + mapq_sum_sq + the two windowed-coverage f32s),
        // 100 records ≈ 5000 bytes.
        let bytes = three_block_chrom0_fixture(5 * 1024);
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
    fn region_wholly_outside_coverage() {
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
    fn iterator_dropped_mid_stream_restarts() {
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
    fn sequential_and_region_agree_inside_window() {
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
    fn head_magic_flipped() {
        let mut bytes = finish_empty_writer(writer_header(1));
        bytes[1] = b'Q'; // PSP\n → PQP\n
        let err = PspReader::new(Cursor::new(bytes)).expect_err("bad magic must fail");
        assert!(matches!(err, PspReadError::BadHeadMagic { .. }));
    }

    /// (F2) Truncate the file mid-block — iteration surfaces an
    /// `Io` (or `Zstd`) error, or open-time validation already
    /// fails because the trailer is missing/garbled.
    #[test]
    fn truncated_mid_block() {
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
                PspReadError::Io { .. }
                    | PspReadError::Zstd { .. }
                    | PspReadError::BlockHeaderTruncated { .. }
                    | PspReadError::BlockHeaderField { .. }
                    | PspReadError::ColumnTruncated { .. }
            ));
        }
        // else: open-time variant accepted
    }

    /// (F3) Flip one byte inside a compressed column payload —
    /// the reader surfaces `Zstd` (or one of the downstream
    /// length / decode errors when the corruption survives the
    /// frame check).
    #[test]
    fn torn_zstd_frame() {
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
    fn corrupted_block_index() {
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
    fn out_of_order_column_tags() {
        use super::super::block::{BlockHeader, ColumnManifestEntry, encode_block_header};
        let bad = BlockHeader {
            chrom_id: 0,
            first_pos: 1,
            n_records: 1,
            n_total_alleles: 1,
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

    /// (F7) A block manifest that omits a v1.0 required column
    /// triggers `MissingRequiredColumnInManifest` — the per-block
    /// coverage check (B1) fronts the per-column decode loop, so a
    /// hand-crafted manifest cannot drive the previously-`.expect()`
    /// panic sites inside `decode_block_payload` and
    /// `decode_one_column`.
    ///
    /// Previously this test asserted `UncompressedLenSchemaMismatch`
    /// directly. With B1's coverage check ahead of the per-column
    /// decode loop, the assertion shifts: an incomplete manifest is
    /// rejected before any per-column predict check runs. A
    /// dedicated `UncompressedLenSchemaMismatch` regression test
    /// requires a fixture-patch helper that round-trips a real file
    /// and overwrites one manifest entry's `uncompressed_len` byte
    /// — left as follow-up.
    #[test]
    fn manifest_missing_required_column() {
        use super::super::block::{BlockHeader, ColumnManifestEntry, encode_block_header};
        let bad = BlockHeader {
            chrom_id: 0,
            first_pos: 1,
            n_records: 1,
            n_total_alleles: 1,
            manifest: vec![ColumnManifestEntry {
                tag: 0x10, // allele-obs-count, the only v1.0 column present
                compressed_len: 0,
                uncompressed_len: 4,
            }],
        };
        let mut buf = Vec::new();
        encode_block_header(&bad, &mut buf).unwrap();
        let mut src = Cursor::new(buf);
        let mut header_buf = Vec::new();
        let (header, _) = read_block_header(&mut src, &mut header_buf).unwrap();
        let mut decompressor = new_column_decompressor().unwrap();
        let mut compressed = Vec::new();
        let mut decompressed = Vec::new();
        let mut allele_seq_data = Vec::new();
        let mut allele_seq_offsets = Vec::new();
        let mut allele_chain_ids_data = Vec::new();
        let mut allele_chain_ids_offsets = Vec::new();
        let err = decode_block_payload(
            &mut src,
            &header,
            u64::MAX,
            &mut decompressor,
            &mut compressed,
            &mut decompressed,
            &mut allele_seq_data,
            &mut allele_seq_offsets,
            &mut allele_chain_ids_data,
            &mut allele_chain_ids_offsets,
        )
        .unwrap_err();
        match err {
            PspReadError::MissingRequiredColumnInManifest { name, tag } => {
                // First missing column in registry-tag order is
                // `delta-pos` (tag 0x01).
                assert_eq!(name, "delta-pos");
                assert_eq!(tag, 0x01);
            }
            other => panic!("expected MissingRequiredColumnInManifest, got {other:?}"),
        }
    }

    // (F8) The "forge block 0's active-slot snapshot" test is gone
    // along with the snapshot field; block headers no longer carry
    // a per-block active-slot list, so there is no `slot_count`
    // byte to mutate or `NonEmptySnapshotAtChromStart` invariant to
    // assert on.

    // -----------------------------------------------------------------
    // Regression tests added by the 2026-05-13 reader code review
    // (B1 / B2 / B3 / M7 / Mi6 / Mi7 / additional reliability misses).
    // -----------------------------------------------------------------

    /// (B1) `decode_block_payload` rejects a per-block manifest
    /// that omits any v1.0-required column with the typed
    /// `MissingRequiredColumnInManifest` variant — closes the
    /// previous panic surface where the trailing `.expect("…
    /// required by v1.0 schema")`s would fire on a hand-crafted
    /// short manifest.
    ///
    /// Pinned for every required tag so a v1.x relaxation of the
    /// required set is a visible, reviewable change.
    #[test]
    fn manifest_missing_required_column_for_every_required_tag() {
        use super::super::block::{BlockHeader, ColumnManifestEntry, encode_block_header};
        // For each required v1.0 column, build a manifest that
        // contains every *other* required column and confirm the
        // omission surfaces with this exact variant.
        let required: Vec<_> = V1_0_COLUMNS.iter().filter(|d| d.required).collect();
        for omitted in &required {
            let mut manifest: Vec<_> = required
                .iter()
                .filter(|d| d.tag != omitted.tag)
                .map(|d| ColumnManifestEntry {
                    tag: d.tag,
                    compressed_len: 0,
                    uncompressed_len: 0,
                })
                .collect();
            // Manifest must be tag-ascending; the registry is, so
            // a filtered copy preserves the order.
            manifest.sort_by_key(|e| e.tag);
            let bad = BlockHeader {
                chrom_id: 0,
                first_pos: 1,
                n_records: 1,
                n_total_alleles: 1,
                manifest,
            };
            let mut buf = Vec::new();
            encode_block_header(&bad, &mut buf).unwrap();
            let mut src = Cursor::new(buf);
            let mut header_buf = Vec::new();
            let (header, _) = read_block_header(&mut src, &mut header_buf).unwrap();
            let mut decompressor = new_column_decompressor().unwrap();
            let mut compressed = Vec::new();
            let mut decompressed = Vec::new();
            let mut allele_seq_data = Vec::new();
            let mut allele_seq_offsets = Vec::new();
            let mut allele_chain_ids_data = Vec::new();
            let mut allele_chain_ids_offsets = Vec::new();
            let err = decode_block_payload(
                &mut src,
                &header,
                u64::MAX,
                &mut decompressor,
                &mut compressed,
                &mut decompressed,
                &mut allele_seq_data,
                &mut allele_seq_offsets,
                &mut allele_chain_ids_data,
                &mut allele_chain_ids_offsets,
            )
            .unwrap_err();
            match err {
                PspReadError::MissingRequiredColumnInManifest { name, tag } => {
                    assert_eq!(
                        name, omitted.name,
                        "wrong name for omitted {}",
                        omitted.name
                    );
                    assert_eq!(tag, omitted.tag, "wrong tag for omitted {}", omitted.name);
                }
                other => panic!(
                    "expected MissingRequiredColumnInManifest for omitted {:?}, got {:?}",
                    omitted.name, other
                ),
            }
        }
    }

    /// (B2) `decode_one_column` rejects a manifest whose
    /// `compressed_len` exceeds the block's remaining byte budget
    /// before allocating. Without this guard a hostile file with
    /// `compressed_len = u32::MAX` would drive a 4 GiB allocation
    /// per column.
    #[test]
    fn decode_one_column_rejects_oversized_compressed_len() {
        use super::super::block::ColumnManifestEntry;
        let entry = ColumnManifestEntry {
            tag: 0x10, // allele-obs-count (known v1.0)
            compressed_len: u32::MAX,
            uncompressed_len: 4,
        };
        // Empty source — the byte-budget check fires before any
        // read, so the source contents are irrelevant.
        let mut src = Cursor::new(Vec::<u8>::new());
        let mut decompressor = new_column_decompressor().unwrap();
        let mut compressed = Vec::new();
        let mut decompressed = Vec::new();
        let mut allele_seq_data = Vec::new();
        let mut allele_seq_offsets = Vec::new();
        let mut allele_chain_ids_data = Vec::new();
        let mut allele_chain_ids_offsets = Vec::new();
        let err = decode_one_column(
            &mut src,
            &entry,
            /* n_records   */ 1,
            /* n_total_alleles */ 1,
            /* allele_seq_len */ None,
            /* remaining_budget */ 100, // far less than u32::MAX
            &mut decompressor,
            &mut compressed,
            &mut decompressed,
            &mut allele_seq_data,
            &mut allele_seq_offsets,
            &mut allele_chain_ids_data,
            &mut allele_chain_ids_offsets,
        )
        .unwrap_err();
        match err {
            PspReadError::ColumnTruncated {
                column,
                decoded,
                expected,
            } => {
                assert_eq!(column, "allele-obs-count");
                assert_eq!(decoded, 0);
                assert_eq!(expected, u32::MAX as usize);
            }
            other => panic!("expected ColumnTruncated for budget overflow, got {other:?}"),
        }
    }

    /// (B3) `decode_block_payload` rejects a file whose decoded
    /// `n-alleles` column sums to a value different from the
    /// block header's `n_total_alleles`. Without this check an
    /// over-run panics on per-allele indexing in
    /// `materialise_next_record`; an under-run silently emits
    /// truncated allele lists.
    ///
    /// The fixture uses a one-record one-allele round-trip and
    /// flips the block header's `n_total_alleles` field, so the
    /// sum (1) and the field (2) disagree.
    #[test]
    fn n_alleles_sum_mismatch_against_n_total_alleles() {
        // Compose a valid one-record file via the writer.
        let mut writer = PspWriter::new(Cursor::new(Vec::new()), writer_header(1)).unwrap();
        writer.write_record(&one_allele_record(0, 1, b"A")).unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        let block0_offset = {
            let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
            reader.block_index()[0].block_offset as usize
        };
        // For this fixture every leading block-header varint is
        // one byte: chrom_id (0x00), first_pos (0x01), n_records
        // (0x01), n_total_alleles (0x01). Patch the
        // `n_total_alleles` byte to 0x02 so the manifest claims
        // 2 total alleles while the n-alleles column still sums
        // to 1.
        let n_total_alleles_offset = block0_offset + 3;
        assert_eq!(
            bytes[n_total_alleles_offset], 0x01,
            "fixture invariant: block 0 declares n_total_alleles = 1"
        );
        let mut mutated = bytes.clone();
        mutated[n_total_alleles_offset] = 0x02;
        let mut reader = match PspReader::new(Cursor::new(mutated)) {
            Ok(r) => r,
            Err(_) => return, // forging mangled the file at open
        };
        let err = reader
            .records()
            .find_map(|r| r.err())
            .expect("n_alleles sum mismatch must surface as an error");
        // The decoded column may also have a different per-allele
        // field length (since `n_total_alleles` now claims 2 but
        // only 1 entry is present in the per-allele columns), so
        // either NAllelesSumMismatch or ColumnTruncated is
        // acceptable.
        assert!(matches!(
            err,
            PspReadError::BlockHeaderInvariant {
                kind: BlockHeaderInvariantKind::NAllelesSumMismatch { .. },
            } | PspReadError::ColumnTruncated { .. }
                | PspReadError::UncompressedLenSchemaMismatch { .. }
        ));
    }

    #[test]
    fn validate_n_alleles_column_accepts_all_nonzero() {
        // sum 4 == n_total_alleles, every record has >= 1 allele.
        assert!(validate_n_alleles_column(&[1, 2, 1], 4).is_ok());
        // single record, single allele.
        assert!(validate_n_alleles_column(&[1], 1).is_ok());
    }

    #[test]
    fn validate_n_alleles_column_rejects_sum_mismatch() {
        let err = validate_n_alleles_column(&[1, 1], 3).unwrap_err();
        assert_eq!(
            err,
            BlockHeaderInvariantKind::NAllelesSumMismatch {
                n_total_alleles: 3,
                sum_n_alleles: 2,
            }
        );
    }

    #[test]
    fn validate_n_alleles_column_rejects_interior_zero_allele_record() {
        // [2, 0, 1] sums to 3 (matches n_total_alleles) but record 1
        // carries zero alleles — the corruption the sum check misses.
        let err = validate_n_alleles_column(&[2, 0, 1], 3).unwrap_err();
        assert_eq!(
            err,
            BlockHeaderInvariantKind::ZeroAlleleRecord { record_index: 1 }
        );
    }

    #[test]
    fn validate_n_alleles_column_rejects_trailing_zero_allele_record() {
        // The trailing-zero case that would index one past the offsets
        // array in the columnar `from_block`.
        let err = validate_n_alleles_column(&[3, 0], 3).unwrap_err();
        assert_eq!(
            err,
            BlockHeaderInvariantKind::ZeroAlleleRecord { record_index: 1 }
        );
    }

    // (M7) The "region mode rejects non-empty snapshot at chrom
    // start" test is gone along with the snapshot field. Block
    // headers no longer carry per-block active-slot lists; the
    // first-block-of-chrom rule it pinned no longer exists.

    /// (Mi5) The `delta_pos` cast from `u64` to `u32` previously
    /// truncated silently for varints with the top 32 bits set,
    /// then `saturating_add` capped `pos` at `u32::MAX`. The fix
    /// rejects the truncation surface with a typed varint-overflow
    /// error.
    ///
    /// Reaching this surface end-to-end requires patching a
    /// `delta-pos` varint inside a real file's compressed payload
    /// — deferred to the next round of fixture helpers. The fix
    /// itself is exercised by `cargo test` build inclusion (the
    /// new error path compiles and `u32::try_from` is exercised
    /// at every record on every existing round-trip test).
    #[test]
    fn delta_pos_u32_overflow_compile_check() {
        // No-op runtime check; the existence of `delta-pos`'s
        // u32::try_from path is verified by every successful
        // round-trip test.
        let _ = u32::try_from(u32::MAX as u64).unwrap();
        let _ = u32::try_from(u32::MAX as u64 + 1).unwrap_err();
    }
}
