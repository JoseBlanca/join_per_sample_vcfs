//! A shared, thread-safe indexed-segment read source.
//!
//! Given a coordinate-sorted + indexed BAM/CRAM and a genomic segment
//! `[start, end]` (1-based inclusive) on a contig, [`AlignmentFile`]
//! returns the reads whose reference footprint overlaps that segment.
//! It exploits the coordinate-sort: index-seek to the segment, then
//! linear-scan in coordinate order until past the segment end.
//!
//! The primitive is built for the SSR Stage-1 fetcher, which queries
//! ~10⁶ tiny loci. The cost it removes versus re-running the existing
//! per-`query` scanners is the per-call file re-open + header re-parse
//! (and, for CRAM, container re-decode of the file header): each
//! [`AlignmentFile`] holds a **pool of idle open readers**, so the
//! Nth segment call on a thread reuses the reader the (N-1)th left
//! behind rather than opening the file again. The file is opened
//! roughly once per concurrent caller (≈ once per thread at the SSR
//! driver's "one segment per thread" rate), then reused.
//!
//! # Thread-safety
//!
//! [`AlignmentFile`] is `Sync`: the index and header map are shared
//! `Arc`s, the filter config is immutable `Copy`, and the reader pool
//! is a `Mutex<Vec<_>>` locked only for the pop/push at the iterator's
//! creation and drop — never during iteration. Many threads call
//! [`AlignmentFile::get_reads_from_segment`] concurrently for
//! different segments through a shared `&AlignmentFile`; each call
//! borrows its own reader, so there is no shared cursor and no lock
//! held across work.
//!
//! # What it is not
//!
//! - **Not a multi-file merge.** One file → its reads for a segment.
//!   A sample spread across several files is each queried separately
//!   and combined by the driver (SSR's reservoir is order-independent,
//!   so it just concatenates).
//! - **Not a clipper.** A read overlapping two queried segments is
//!   yielded *whole* by both calls; the consumer decides what to do
//!   with the part outside its segment.
//!
//! The per-record decode reuses the existing per-format helpers
//! ([`crate::bam::alignment_input`]'s `classify_pre_decode` /
//! `record_buf_to_mapped_read`, [`crate::bam::bam_input`]'s
//! `query_interval`, and `ContigInterval::overlaps_record`); only the
//! reader ownership + pooling is new here. The existing
//! `OwnedIndexed{Bam,Cram}Records` scanners are left in place — they
//! still back `AlignmentMergedReader::query` until the SNP `--regions`
//! path is retrofitted onto this primitive (a separate, measured
//! change).

// The primitive's only live consumers are its own tests until the SSR
// Stage-1 fetcher (increment #4) and the SNP `--regions` retrofit
// (increment #5) are wired onto it — both deferred to their own
// measured changes per the implementation plan. Without that wiring the
// public surface (`get_reads_from_segment`, the CRAM path, the
// dispatch enum) is reachable only from `#[cfg(test)]`, which reads as
// dead code to the non-test lib build. Allow it module-wide rather than
// scatter per-item attributes; remove when the consumer lands.
#![allow(dead_code)]

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, SeekFrom};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::vec;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_csi::binning_index::BinningIndex;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::alignment_input::{
    AlignmentMergedReaderConfig, ContigInterval, MappedRead, PreDecodeOutcome, classify_pre_decode,
    record_buf_to_mapped_read,
};
use super::bam_input::{BamIndex, open_bam_reader_with_header, query_interval};
use super::cram_input::open_cram_reader_with_header;
use super::errors::AlignmentInputError;
use super::index_preflight::{AlignmentFileKind, AlignmentIndex};

// ---------------------------------------------------------------------
// Header ref map
// ---------------------------------------------------------------------

/// The parsed SAM header plus a contig-name → reference-id lookup.
///
/// The header is the one the decoder needs to materialise records
/// (`read_record_buf` / `try_from_alignment_record` both take it); the
/// name map turns a caller's `chrom` string into the integer
/// `reference_sequence_id` the index query and the record's own
/// `reference_sequence_id()` use. That id equals the `@SQ` order
/// index, which equals the index into the canonical `ContigList` the
/// rest of the pipeline keys on, so [`MappedRead::ref_id`] needs no
/// remapping.
///
/// Built once in [`AlignmentFile::from_input`] and shared by `Arc`
/// across every pooled reader.
pub(crate) struct HeaderRefMap {
    header: Arc<sam::Header>,
    name_to_id: HashMap<String, usize>,
}

impl HeaderRefMap {
    fn from_header(header: Arc<sam::Header>) -> Self {
        let name_to_id = header
            .reference_sequences()
            .keys()
            .enumerate()
            .map(|(id, name)| (String::from_utf8_lossy(name.as_ref()).into_owned(), id))
            .collect();
        Self { header, name_to_id }
    }

    fn ref_id(&self, chrom: &str) -> Option<usize> {
        self.name_to_id.get(chrom).copied()
    }

    fn header(&self) -> &sam::Header {
        &self.header
    }
}

// ---------------------------------------------------------------------
// AlignmentFile
// ---------------------------------------------------------------------

/// A coordinate-sorted, indexed alignment file (BAM or CRAM) that
/// serves per-segment read queries. Immutable + `Sync`; share it by
/// `&` across threads. See the module docs for the pooling and
/// thread-safety design.
pub(crate) enum AlignmentFile {
    Bam(BamFile),
    Cram(CramFile),
}

impl AlignmentFile {
    /// Build an [`AlignmentFile`] from one input's already-validated
    /// handles: its `path`, the shared parsed `header`, the loaded
    /// `index`, the per-read filter `cfg`, and the
    /// `source_file_index` stamped into every [`MappedRead`] this file
    /// yields (the caller's position for this file in its input list).
    ///
    /// `repository` is the FASTA reference; **required for CRAM**
    /// (slice decoding consults it) and ignored for BAM. The caller
    /// builds it once via
    /// [`crate::bam::alignment_input::build_fasta_repository`] and
    /// shares the clone (an `Arc` bump).
    ///
    /// # Errors
    ///
    /// - [`AlignmentInputError::UnsupportedExtension`] — `path`'s
    ///   extension is neither `.cram` nor `.bam`.
    /// - [`AlignmentInputError::AlignmentIndexFormatMismatch`] — the
    ///   `index` variant's format disagrees with the file extension.
    /// - [`AlignmentInputError::MissingCramReference`] — a CRAM input
    ///   was given no `repository`.
    pub(crate) fn from_input(
        path: PathBuf,
        header: Arc<sam::Header>,
        index: AlignmentIndex,
        repository: Option<fasta::Repository>,
        cfg: AlignmentMergedReaderConfig,
        source_file_index: usize,
    ) -> Result<Self, AlignmentInputError> {
        let kind = AlignmentFileKind::from_path(&path)
            .ok_or_else(|| AlignmentInputError::UnsupportedExtension { path: path.clone() })?;
        let ref_map = Arc::new(HeaderRefMap::from_header(header));

        match (kind, index) {
            (AlignmentFileKind::Bam, AlignmentIndex::BamCsi(idx)) => Ok(Self::Bam(BamFile {
                path,
                index: BamIndex::Csi(idx),
                ref_map,
                cfg,
                source_file_index,
                readers_pool: Mutex::new(Vec::new()),
            })),
            (AlignmentFileKind::Bam, AlignmentIndex::BamBai(idx)) => Ok(Self::Bam(BamFile {
                path,
                index: BamIndex::Bai(idx),
                ref_map,
                cfg,
                source_file_index,
                readers_pool: Mutex::new(Vec::new()),
            })),
            (AlignmentFileKind::Cram, AlignmentIndex::Crai(idx)) => {
                let repository = repository.ok_or_else(|| {
                    AlignmentInputError::MissingCramReference { path: path.clone() }
                })?;
                Ok(Self::Cram(CramFile {
                    path,
                    index: idx,
                    repository,
                    ref_map,
                    cfg,
                    source_file_index,
                    readers_pool: Mutex::new(Vec::new()),
                }))
            }
            // Format / index disagreement — the driver's index loader
            // and the file extension say different things.
            (kind, index) => Err(AlignmentInputError::AlignmentIndexFormatMismatch {
                path,
                file_format: kind.display_name(),
                index_format: index.display_name(),
            }),
        }
    }

    /// Reads whose reference footprint overlaps `[start, end]`
    /// (1-based inclusive) on `chrom`, in coordinate order. Borrows a
    /// pooled reader, index-seeks it to the segment, and returns an
    /// iterator that streams until past the segment end; the reader
    /// returns to the pool when the iterator drops.
    ///
    /// Callable concurrently from many threads for different segments
    /// through a shared `&AlignmentFile`.
    ///
    /// # Errors
    ///
    /// - [`AlignmentInputError::InvalidSegment`] — not a valid 1-based
    ///   inclusive range.
    /// - [`AlignmentInputError::ContigNotInList`] — `chrom` is absent
    ///   from the header.
    /// - [`AlignmentInputError::OpenFailed`] / [`AlignmentInputError::Io`]
    ///   — opening a fresh pooled reader, or the BAM index query,
    ///   failed.
    pub(crate) fn get_reads_from_segment(
        &self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<MappedReadsInSegment<'_>, AlignmentInputError> {
        match self {
            Self::Bam(file) => Ok(MappedReadsInSegment::Bam(
                file.get_reads_from_segment(chrom, start, end)?,
            )),
            Self::Cram(file) => Ok(MappedReadsInSegment::Cram(
                file.get_reads_from_segment(chrom, start, end)?,
            )),
        }
    }
}

/// Validate a caller's `(chrom, start, end)` into a [`ContigInterval`]
/// plus the resolved reference id. Shared by both formats.
fn resolve_segment(
    ref_map: &HeaderRefMap,
    chrom: &str,
    start: u32,
    end: u32,
) -> Result<(usize, ContigInterval), AlignmentInputError> {
    if start == 0 || start > end {
        return Err(AlignmentInputError::InvalidSegment {
            chrom: chrom.to_string(),
            start,
            end,
        });
    }
    let target = ref_map
        .ref_id(chrom)
        .ok_or_else(|| AlignmentInputError::ContigNotInList {
            contig: chrom.to_string(),
            known_contigs: ref_map.name_to_id.len(),
        })?;
    Ok((target, ContigInterval { start, end }))
}

// ---------------------------------------------------------------------
// BAM
// ---------------------------------------------------------------------

/// A pooled, re-seekable BAM reader plus the metadata its segment
/// queries need. See the module docs.
pub(crate) struct BamFile {
    path: PathBuf,
    index: BamIndex,
    ref_map: Arc<HeaderRefMap>,
    cfg: AlignmentMergedReaderConfig,
    source_file_index: usize,
    readers_pool: Mutex<Vec<BamReaderHandle>>,
}

/// One idle open BAM reader, positioned past the header and ready to
/// be seeked to a chunk's start.
struct BamReaderHandle {
    reader: bam::io::Reader<bgzf::io::Reader<File>>,
}

impl BamFile {
    fn get_reads_from_segment(
        &self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<BamSegmentReads<'_>, AlignmentInputError> {
        let (target, segment) = resolve_segment(&self.ref_map, chrom, start, end)?;

        // Index → candidate chunks. The query is narrowed to the
        // segment; the per-record overlap filter below drops what the
        // bin-granular chunk set over-returns at the edges.
        let interval = query_interval(Some(segment));
        let chunks = match &self.index {
            BamIndex::Bai(idx) => idx.query(target, interval),
            BamIndex::Csi(idx) => idx.query(target, interval),
        }
        .map_err(|source| AlignmentInputError::Io {
            path: self.path.clone(),
            source,
        })?;

        let handle = self.borrow_handle()?;
        Ok(BamSegmentReads {
            file: self,
            handle: Some(handle),
            chunks: chunks.into_iter(),
            current_chunk_end: None,
            target_reference_sequence_id: target,
            segment,
            done: false,
        })
    }

    /// Pop an idle reader from the pool, or open a fresh one if the
    /// pool is empty (the only file-open; happens ~once per concurrent
    /// caller).
    fn borrow_handle(&self) -> Result<BamReaderHandle, AlignmentInputError> {
        if let Some(handle) = self
            .readers_pool
            .lock()
            .expect("reader pool not poisoned")
            .pop()
        {
            return Ok(handle);
        }
        let (reader, _header) = open_bam_reader_with_header(&self.path)?;
        Ok(BamReaderHandle { reader })
    }

    fn return_handle(&self, handle: BamReaderHandle) {
        if let Ok(mut pool) = self.readers_pool.lock() {
            pool.push(handle);
        }
    }

    #[cfg(test)]
    fn pool_len(&self) -> usize {
        self.readers_pool
            .lock()
            .expect("reader pool not poisoned")
            .len()
    }
}

/// Per-call BAM iterator. `(target, segment)` are fixed at creation;
/// the only mutable state is the chunk cursor and the reader handle.
/// On drop it returns the borrowed reader to its [`BamFile`]'s pool.
pub(crate) struct BamSegmentReads<'a> {
    file: &'a BamFile,
    /// `Some` while borrowed; taken in `Drop` to return to the pool.
    handle: Option<BamReaderHandle>,
    chunks: vec::IntoIter<Chunk>,
    /// Virtual position where the current chunk ends, or `None` to
    /// advance to the next chunk on the next poll.
    current_chunk_end: Option<bgzf::VirtualPosition>,
    target_reference_sequence_id: usize,
    segment: ContigInterval,
    /// Latched once we pass the segment end (sort order guarantees
    /// nothing later overlaps) or exhaust the chunks.
    done: bool,
}

impl BamSegmentReads<'_> {
    fn io_error(&self, source: io::Error) -> AlignmentInputError {
        AlignmentInputError::Io {
            path: self.file.path.clone(),
            source,
        }
    }
}

impl Iterator for BamSegmentReads<'_> {
    type Item = Result<MappedRead, AlignmentInputError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            // 1. Land inside a chunk (advance + seek if needed).
            let chunk_end = match self.current_chunk_end {
                Some(end) => end,
                None => {
                    let Some(chunk) = self.chunks.next() else {
                        self.done = true;
                        return None;
                    };
                    let reader = &mut self.handle.as_mut().expect("handle held").reader;
                    if let Err(e) = reader.get_mut().seek(chunk.start()) {
                        return Some(Err(self.io_error(e)));
                    }
                    self.current_chunk_end = Some(chunk.end());
                    chunk.end()
                }
            };

            // 2. If the previous read pushed us past the chunk end,
            //    fall through to the next chunk.
            {
                let reader = &self.handle.as_ref().expect("handle held").reader;
                if reader.get_ref().virtual_position() >= chunk_end {
                    self.current_chunk_end = None;
                    continue;
                }
            }

            // 3. Read one record at the current virtual position.
            let mut record = sam::alignment::RecordBuf::default();
            let read = {
                let reader = &mut self.handle.as_mut().expect("handle held").reader;
                reader.read_record_buf(self.file.ref_map.header(), &mut record)
            };
            match read {
                Ok(0) => {
                    // EOF mid-chunk — advance to the next chunk.
                    self.current_chunk_end = None;
                    continue;
                }
                Ok(_) => {
                    // Different contig (a chunk straddling contigs) — skip.
                    if record.reference_sequence_id() != Some(self.target_reference_sequence_id) {
                        continue;
                    }
                    // Past the segment end → nothing later overlaps.
                    if let Some(astart) = record.alignment_start()
                        && usize::from(astart) as u64 > u64::from(self.segment.end)
                    {
                        self.done = true;
                        return None;
                    }
                    // Footprint does not overlap (ends before start, or
                    // the chunk over-returned) — skip.
                    if !self.segment.overlaps_record(&record) {
                        continue;
                    }
                    // Pre-decode flag/mapq filter.
                    if let PreDecodeOutcome::Drop(_) = classify_pre_decode(&self.file.cfg, &record)
                    {
                        continue;
                    }
                    // Min read length (cheap, pre-conversion).
                    if let Some(min) = self.file.cfg.min_read_length
                        && (record.sequence().as_ref().len() as u32) < min
                    {
                        continue;
                    }
                    match record_buf_to_mapped_read(&record, self.file.source_file_index) {
                        Ok(mapped) => return Some(Ok(mapped)),
                        Err(source) => {
                            self.done = true;
                            return Some(Err(AlignmentInputError::MalformedRecord {
                                path: self.file.path.clone(),
                                qname: record
                                    .name()
                                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned()),
                                source,
                            }));
                        }
                    }
                }
                Err(e) => {
                    self.done = true;
                    return Some(Err(self.io_error(e)));
                }
            }
        }
    }
}

impl Drop for BamSegmentReads<'_> {
    fn drop(&mut self) {
        if let Some(handle) = self.handle.take() {
            self.file.return_handle(handle);
        }
    }
}

// ---------------------------------------------------------------------
// CRAM
// ---------------------------------------------------------------------

/// A pooled, re-seekable CRAM reader plus the metadata its segment
/// queries need. See the module docs.
pub(crate) struct CramFile {
    path: PathBuf,
    index: Arc<cram::crai::Index>,
    repository: fasta::Repository,
    ref_map: Arc<HeaderRefMap>,
    cfg: AlignmentMergedReaderConfig,
    source_file_index: usize,
    readers_pool: Mutex<Vec<CramReaderHandle>>,
}

/// One idle open CRAM reader, positioned past the file definition +
/// header and ready to be seeked to a container offset.
struct CramReaderHandle {
    reader: cram::io::Reader<File>,
}

impl CramFile {
    fn get_reads_from_segment(
        &self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<CramSegmentReads<'_>, AlignmentInputError> {
        let (target, segment) = resolve_segment(&self.ref_map, chrom, start, end)?;
        let handle = self.borrow_handle()?;
        Ok(CramSegmentReads {
            file: self,
            handle: Some(handle),
            next_index_record: 0,
            pending: Vec::new().into_iter(),
            target_reference_sequence_id: target,
            segment,
            done: false,
        })
    }

    fn borrow_handle(&self) -> Result<CramReaderHandle, AlignmentInputError> {
        if let Some(handle) = self
            .readers_pool
            .lock()
            .expect("reader pool not poisoned")
            .pop()
        {
            return Ok(handle);
        }
        let (reader, _header) = open_cram_reader_with_header(&self.path, Some(&self.repository))?;
        Ok(CramReaderHandle { reader })
    }

    fn return_handle(&self, handle: CramReaderHandle) {
        if let Ok(mut pool) = self.readers_pool.lock() {
            pool.push(handle);
        }
    }

    #[cfg(test)]
    fn pool_len(&self) -> usize {
        self.readers_pool
            .lock()
            .expect("reader pool not poisoned")
            .len()
    }
}

/// Per-call CRAM iterator. Walks the `.crai` for containers on the
/// target contig, seeks to each, decodes its records, and buffers the
/// segment-overlapping survivors in `pending`. On drop it returns the
/// borrowed reader to its [`CramFile`]'s pool.
pub(crate) struct CramSegmentReads<'a> {
    file: &'a CramFile,
    handle: Option<CramReaderHandle>,
    /// Cursor into the `.crai` (held by `Arc` on the file, walked by
    /// integer to keep the struct free of self-referential borrows).
    next_index_record: usize,
    /// Converted, segment-overlapping reads from the last decoded
    /// container, drained before the next container is read.
    pending: vec::IntoIter<MappedRead>,
    target_reference_sequence_id: usize,
    segment: ContigInterval,
    done: bool,
}

impl CramSegmentReads<'_> {
    fn io_error(&self, source: io::Error) -> AlignmentInputError {
        AlignmentInputError::Io {
            path: self.file.path.clone(),
            source,
        }
    }

    /// Decode the next target-contig container into `pending`. Returns
    /// `Ok(true)` once the index is exhausted.
    fn refill(&mut self) -> Result<bool, AlignmentInputError> {
        loop {
            // Clone the small index record so the `Arc<Index>` borrow
            // is released before we mutate the cursor / buffers.
            let record = match self.file.index.as_slice().get(self.next_index_record) {
                Some(r) => r.clone(),
                None => return Ok(true),
            };
            self.next_index_record += 1;

            if record.reference_sequence_id() != Some(self.target_reference_sequence_id) {
                continue;
            }

            // Container-level skip: a container whose footprint is
            // disjoint from the segment holds nothing we want, so we
            // can avoid the seek + decode entirely.
            if let Some(container_start) = record.alignment_start() {
                let container_start = usize::from(container_start) as u64;
                let span = record.alignment_span() as u64;
                if span > 0 {
                    let container_end = container_start + span - 1;
                    if container_end < u64::from(self.segment.start)
                        || container_start > u64::from(self.segment.end)
                    {
                        continue;
                    }
                }
            }

            // Seek + decode. The reader borrow is confined to this
            // block; container ownership outlives it.
            let mut container = cram::io::reader::Container::default();
            let decoded = {
                let reader = &mut self.handle.as_mut().expect("handle held").reader;
                reader.seek(SeekFrom::Start(record.offset())).map_err(|e| {
                    AlignmentInputError::Io {
                        path: self.file.path.clone(),
                        source: e,
                    }
                })?;
                reader
                    .read_container(&mut container)
                    .map_err(|e| AlignmentInputError::Io {
                        path: self.file.path.clone(),
                        source: e,
                    })?
            };
            if decoded == 0 {
                return Ok(true);
            }

            let header = self.file.ref_map.header();
            let compression_header = container
                .compression_header()
                .map_err(|e| self.io_error(e))?;
            let mut converted: Vec<MappedRead> = Vec::new();
            for slice_result in container.slices() {
                let slice = slice_result.map_err(|e| self.io_error(e))?;
                let (core_data_src, external_data_srcs) =
                    slice.decode_blocks().map_err(|e| self.io_error(e))?;
                let cram_records = slice
                    .records(
                        self.file.repository.clone(),
                        header,
                        &compression_header,
                        &core_data_src,
                        &external_data_srcs,
                    )
                    .map_err(|e| self.io_error(e))?;
                for cram_record in &cram_records {
                    let record =
                        sam::alignment::RecordBuf::try_from_alignment_record(header, cram_record)
                            .map_err(|e| self.io_error(e))?;
                    if record.reference_sequence_id() != Some(self.target_reference_sequence_id) {
                        continue;
                    }
                    if !self.segment.overlaps_record(&record) {
                        continue;
                    }
                    if let PreDecodeOutcome::Drop(_) = classify_pre_decode(&self.file.cfg, &record)
                    {
                        continue;
                    }
                    if let Some(min) = self.file.cfg.min_read_length
                        && (record.sequence().as_ref().len() as u32) < min
                    {
                        continue;
                    }
                    let mapped = record_buf_to_mapped_read(&record, self.file.source_file_index)
                        .map_err(|source| AlignmentInputError::MalformedRecord {
                            path: self.file.path.clone(),
                            qname: record
                                .name()
                                .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned()),
                            source,
                        })?;
                    converted.push(mapped);
                }
            }
            self.pending = converted.into_iter();
            return Ok(false);
        }
    }
}

impl Iterator for CramSegmentReads<'_> {
    type Item = Result<MappedRead, AlignmentInputError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(mapped) = self.pending.next() {
                return Some(Ok(mapped));
            }
            if self.done {
                return None;
            }
            match self.refill() {
                Ok(true) => {
                    self.done = true;
                    return None;
                }
                Ok(false) => continue,
                Err(e) => {
                    self.done = true;
                    return Some(Err(e));
                }
            }
        }
    }
}

impl Drop for CramSegmentReads<'_> {
    fn drop(&mut self) {
        if let Some(handle) = self.handle.take() {
            self.file.return_handle(handle);
        }
    }
}

// ---------------------------------------------------------------------
// The per-call iterator (format-dispatching wrapper)
// ---------------------------------------------------------------------

/// The reads overlapping a segment, streamed in coordinate order.
/// Yields `Result<MappedRead, _>`; a per-record decode failure surfaces
/// as `Some(Err(_))` and fuses the iterator. Borrows the
/// [`AlignmentFile`] for the call's duration and returns its reader to
/// the pool on drop.
pub(crate) enum MappedReadsInSegment<'a> {
    Bam(BamSegmentReads<'a>),
    Cram(CramSegmentReads<'a>),
}

impl Iterator for MappedReadsInSegment<'_> {
    type Item = Result<MappedRead, AlignmentInputError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Bam(reads) => reads.next(),
            Self::Cram(reads) => reads.next(),
        }
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::num::NonZero;
    use std::path::Path;

    use noodles_core::Position;
    use noodles_csi::binning_index::Indexer;
    use noodles_csi::binning_index::index::reference_sequence::index::BinnedIndex;
    use noodles_sam::alignment::Record as _;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::MappingQuality;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::ReferenceSequence;
    use rayon::prelude::*;
    use tempfile::TempDir;

    use super::*;
    use crate::bam::alignment_input::build_fasta_repository;
    use crate::pileup::per_sample::cram_files::{
        ContigSpec, HeaderOverrides, build_cram, build_fasta,
    };

    const CONTIG_LEN: usize = 200;

    /// Filters off — these overlap tests probe the seek + overlap
    /// logic, not the flag/mapq cascade (which has its own test).
    fn permissive_config() -> AlignmentMergedReaderConfig {
        AlignmentMergedReaderConfig {
            min_mapq: None,
            min_read_length: None,
            drop_qc_fail: false,
            drop_duplicate: false,
            max_read_mismatch_fraction: None,
            mismatch_bq_floor: 0,
        }
    }

    fn assert_sync<T: Sync>() {}

    #[test]
    fn alignment_file_is_sync() {
        // The whole point of the primitive: shareable by `&` across
        // threads. A regression that made a field non-`Sync` (e.g. a
        // `RefCell` in the pool) would fail to compile here.
        assert_sync::<AlignmentFile>();
    }

    // --- BAM fixtures -------------------------------------------------

    fn header_two_contigs() -> sam::Header {
        sam::Header::builder()
            .set_header(Default::default())
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZero::new(CONTIG_LEN).unwrap()),
            )
            .add_reference_sequence(
                "chr2",
                Map::<ReferenceSequence>::new(NonZero::new(CONTIG_LEN).unwrap()),
            )
            .build()
    }

    /// A mapped primary read of `len` `M`-ops (footprint
    /// `[start, start + len - 1]`) with sequence + qualities so the
    /// conversion to `MappedRead` has bytes to copy.
    fn aln_record(qname: &str, ref_id: usize, start: usize, len: usize, mapq: u8) -> RecordBuf {
        RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_reference_sequence_id(ref_id)
            .set_flags(Flags::default())
            .set_mapping_quality(MappingQuality::new(mapq).expect("mapq in range"))
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar([Op::new(Kind::Match, len)].into_iter().collect())
            .set_sequence(Sequence::from(vec![b'A'; len]))
            .set_quality_scores(QualityScores::from(vec![30u8; len]))
            .build()
    }

    fn write_bam(records: &[RecordBuf], header: &sam::Header) -> (TempDir, PathBuf) {
        let dir = TempDir::new().expect("tempdir");
        let bam_path = dir.path().join("sample.bam");
        let file = File::create(&bam_path).expect("create bam");
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(header).expect("write header");
        for record in records {
            writer
                .write_alignment_record(header, record)
                .expect("write record");
        }
        writer.try_finish().expect("finish writer");
        (dir, bam_path)
    }

    fn build_csi_in_memory(bam_path: &Path) -> noodles_csi::Index {
        let mut reader = bam::io::reader::Builder
            .build_from_path(bam_path)
            .expect("open bam");
        let parsed_header = reader.read_header().expect("read header");
        let mut indexer: Indexer<BinnedIndex> = Indexer::default();
        let mut chunk_start = reader.get_ref().virtual_position();
        let mut record = bam::Record::default();
        while reader.read_record(&mut record).expect("read") != 0 {
            let chunk_end = reader.get_ref().virtual_position();
            let alignment_context = match (
                record.reference_sequence_id().transpose().expect("ref"),
                record.alignment_start().transpose().expect("start"),
                record.alignment_end().transpose().expect("end"),
            ) {
                (Some(id), Some(start), Some(end)) => {
                    Some((id, start, end, !record.flags().is_unmapped()))
                }
                _ => None,
            };
            indexer
                .add_record(alignment_context, Chunk::new(chunk_start, chunk_end))
                .expect("add");
            chunk_start = chunk_end;
        }
        indexer.build(parsed_header.reference_sequences().len())
    }

    /// Build a `BAM`-backed [`AlignmentFile`] over `records`. Keeps the
    /// `TempDir` alive via the returned handle.
    fn bam_alignment_file(
        records: &[RecordBuf],
        cfg: AlignmentMergedReaderConfig,
    ) -> (TempDir, AlignmentFile) {
        let header = header_two_contigs();
        let (dir, bam_path) = write_bam(records, &header);
        let csi = build_csi_in_memory(&bam_path);
        let file = AlignmentFile::from_input(
            bam_path,
            Arc::new(header),
            AlignmentIndex::BamCsi(Arc::new(csi)),
            None,
            cfg,
            /* source_file_index = */ 7,
        )
        .expect("build alignment file");
        (dir, file)
    }

    fn pool_len(file: &AlignmentFile) -> usize {
        match file {
            AlignmentFile::Bam(f) => f.pool_len(),
            AlignmentFile::Cram(f) => f.pool_len(),
        }
    }

    /// Drain a segment query into the sorted list of `(pos, qname)` it
    /// yielded, asserting no error.
    fn drain(file: &AlignmentFile, chrom: &str, start: u32, end: u32) -> Vec<(u64, String)> {
        file.get_reads_from_segment(chrom, start, end)
            .expect("open segment")
            .map(|r| {
                let read = r.expect("decode read");
                (read.pos, String::from_utf8_lossy(&read.qname).into_owned())
            })
            .collect()
    }

    // --- BAM tests ----------------------------------------------------

    #[test]
    fn bam_returns_exactly_overlapping_reads_with_inclusive_edges() {
        // chr1 footprints (Match(4) → [start, start+3]):
        // a[1,4] b[6,9] c[8,11] d[20,23] e[40,43]; plus a chr2 read.
        let records = [
            aln_record("a", 0, 1, 4, 60),
            aln_record("b", 0, 6, 4, 60),
            aln_record("c", 0, 8, 4, 60),
            aln_record("d", 0, 20, 4, 60),
            aln_record("e", 0, 40, 4, 60),
            aln_record("z", 1, 5, 4, 60),
        ];
        let (_dir, file) = bam_alignment_file(&records, permissive_config());

        // [8,25]: b spans the start edge (ends at 9 >= 8), c & d inside;
        // a (ends 4) and e (starts 40) excluded.
        let got = drain(&file, "chr1", 8, 25);
        assert_eq!(
            got,
            vec![(6, "b".into()), (8, "c".into()), (20, "d".into())]
        );

        // Source-file index is stamped onto every read.
        let first = file
            .get_reads_from_segment("chr1", 8, 25)
            .unwrap()
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(first.source_file_index, 7);
    }

    #[test]
    fn bam_segment_boundaries_are_one_based_inclusive() {
        let records = [aln_record("a", 0, 1, 4, 60)]; // footprint [1,4]
        let (_dir, file) = bam_alignment_file(&records, permissive_config());

        // End exactly on the read's last base → included.
        assert_eq!(drain(&file, "chr1", 4, 4), vec![(1, "a".into())]);
        // One past the read's last base → excluded.
        assert!(drain(&file, "chr1", 5, 5).is_empty());
        // Start exactly on the read's first base → included.
        assert_eq!(drain(&file, "chr1", 1, 1), vec![(1, "a".into())]);
    }

    #[test]
    fn bam_read_spanning_two_segments_is_yielded_by_both_whole() {
        // One long read [10,49] (Match(40)).
        let records = [aln_record("long", 0, 10, 40, 60)];
        let (_dir, file) = bam_alignment_file(&records, permissive_config());

        let left = drain(&file, "chr1", 12, 15);
        let right = drain(&file, "chr1", 30, 35);
        assert_eq!(left, vec![(10, "long".into())]);
        assert_eq!(right, vec![(10, "long".into())]);
    }

    #[test]
    fn bam_empty_segment_yields_nothing_and_returns_handle() {
        let records = [aln_record("a", 0, 1, 4, 60), aln_record("b", 0, 100, 4, 60)];
        let (_dir, file) = bam_alignment_file(&records, permissive_config());

        // A gap between the two reads.
        assert!(drain(&file, "chr1", 50, 60).is_empty());
        // The borrowed reader came back to the pool.
        assert_eq!(pool_len(&file), 1);
    }

    #[test]
    fn bam_pool_opens_once_and_returns_to_resting_size() {
        let records = [aln_record("a", 0, 1, 4, 60), aln_record("b", 0, 20, 4, 60)];
        let (_dir, file) = bam_alignment_file(&records, permissive_config());

        // Fresh file: empty pool.
        assert_eq!(pool_len(&file), 0);

        for _ in 0..5 {
            {
                let reads = file.get_reads_from_segment("chr1", 1, 200).unwrap();
                // During the borrow the pool is empty (the reader is out).
                assert_eq!(pool_len(&file), 0);
                let count = reads.count();
                assert_eq!(count, 2);
            }
            // After drop, exactly one reader is resting — never grew.
            assert_eq!(pool_len(&file), 1);
        }
    }

    #[test]
    fn bam_parallel_segments_match_sequential() {
        // Reads every 5 bp across chr1; query overlapping windows.
        let records: Vec<RecordBuf> = (0..30)
            .map(|i| aln_record(&format!("r{i}"), 0, 1 + i * 5, 4, 60))
            .collect();
        let (_dir, file) = bam_alignment_file(&records, permissive_config());

        let segments: Vec<(u32, u32)> = (0..30).map(|i| (1 + i * 5, 10 + i * 5)).collect();

        let sequential: Vec<Vec<(u64, String)>> = segments
            .iter()
            .map(|&(s, e)| drain(&file, "chr1", s, e))
            .collect();

        // Shared `&AlignmentFile` across rayon workers — the proof the
        // `Sync` design holds and there is no shared cursor.
        let parallel: Vec<Vec<(u64, String)>> = segments
            .par_iter()
            .map(|&(s, e)| drain(&file, "chr1", s, e))
            .collect();

        assert_eq!(parallel, sequential);
        // Every read overlaps at least one window, so each query is
        // non-trivial somewhere.
        assert!(sequential.iter().any(|reads| !reads.is_empty()));
    }

    #[test]
    fn bam_target_contig_filter_excludes_other_contigs() {
        let records = [aln_record("a", 0, 10, 4, 60), aln_record("z", 1, 10, 4, 60)];
        let (_dir, file) = bam_alignment_file(&records, permissive_config());

        assert_eq!(drain(&file, "chr1", 1, 200), vec![(10, "a".into())]);
        assert_eq!(drain(&file, "chr2", 1, 200), vec![(10, "z".into())]);
    }

    #[test]
    fn bam_low_mapq_read_is_filtered_under_default_config() {
        let records = [
            aln_record("low", 0, 10, 4, 5),
            aln_record("high", 0, 12, 4, 60),
        ];
        // Isolate the MAPQ gate (min_mapq = 20) from the other filters.
        let cfg = AlignmentMergedReaderConfig {
            min_mapq: Some(20),
            ..permissive_config()
        };
        let (_dir, file) = bam_alignment_file(&records, cfg);

        assert_eq!(drain(&file, "chr1", 1, 200), vec![(12, "high".into())]);
    }

    #[test]
    fn bam_invalid_segment_is_rejected() {
        let records = [aln_record("a", 0, 10, 4, 60)];
        let (_dir, file) = bam_alignment_file(&records, permissive_config());

        // start > end.
        assert!(matches!(
            file.get_reads_from_segment("chr1", 20, 10),
            Err(AlignmentInputError::InvalidSegment { .. })
        ));
        // start == 0 (not 1-based).
        assert!(matches!(
            file.get_reads_from_segment("chr1", 0, 10),
            Err(AlignmentInputError::InvalidSegment { .. })
        ));
    }

    #[test]
    fn bam_unknown_contig_is_rejected() {
        let records = [aln_record("a", 0, 10, 4, 60)];
        let (_dir, file) = bam_alignment_file(&records, permissive_config());

        assert!(matches!(
            file.get_reads_from_segment("nope", 1, 10),
            Err(AlignmentInputError::ContigNotInList { .. })
        ));
    }

    // --- CRAM fixtures + tests ----------------------------------------

    /// Build a CRAM-backed [`AlignmentFile`] over `records` on a single
    /// contig `chr1`. Keeps both tempdirs (FASTA + CRAM) alive.
    fn cram_alignment_file(
        records: &[RecordBuf],
        cfg: AlignmentMergedReaderConfig,
    ) -> (TempDir, TempDir, AlignmentFile) {
        let contigs = [ContigSpec {
            name: "chr1".into(),
            length: CONTIG_LEN as u64,
        }];
        let (fasta_dir, fasta_path) = build_fasta(&contigs).expect("build fasta");
        let (cram_dir, cram_path) =
            build_cram(&fasta_path, &contigs, &HeaderOverrides::default(), records)
                .expect("build cram");

        let repository = build_fasta_repository(&fasta_path).expect("repository");
        let (_reader, header) =
            open_cram_reader_with_header(&cram_path, Some(&repository)).expect("read header");
        let index = cram::fs::index(&cram_path).expect("build crai");

        let file = AlignmentFile::from_input(
            cram_path,
            Arc::new(header),
            AlignmentIndex::Crai(Arc::new(index)),
            Some(repository),
            cfg,
            /* source_file_index = */ 3,
        )
        .expect("build alignment file");
        (fasta_dir, cram_dir, file)
    }

    #[test]
    fn cram_returns_exactly_overlapping_reads() {
        // chr1 footprints: a[1,4] b[6,9] c[8,11] d[20,23] e[40,43].
        let records = [
            aln_record("a", 0, 1, 4, 60),
            aln_record("b", 0, 6, 4, 60),
            aln_record("c", 0, 8, 4, 60),
            aln_record("d", 0, 20, 4, 60),
            aln_record("e", 0, 40, 4, 60),
        ];
        let (_fa, _cram, file) = cram_alignment_file(&records, permissive_config());

        let got = drain(&file, "chr1", 8, 25);
        assert_eq!(
            got,
            vec![(6, "b".into()), (8, "c".into()), (20, "d".into())]
        );

        let first = file
            .get_reads_from_segment("chr1", 8, 25)
            .unwrap()
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(first.source_file_index, 3);
    }

    #[test]
    fn cram_read_spanning_two_segments_is_yielded_by_both() {
        let records = [aln_record("long", 0, 10, 40, 60)];
        let (_fa, _cram, file) = cram_alignment_file(&records, permissive_config());

        assert_eq!(drain(&file, "chr1", 12, 15), vec![(10, "long".into())]);
        assert_eq!(drain(&file, "chr1", 30, 35), vec![(10, "long".into())]);
    }

    #[test]
    fn cram_empty_segment_yields_nothing_and_returns_handle() {
        let records = [aln_record("a", 0, 1, 4, 60), aln_record("b", 0, 100, 4, 60)];
        let (_fa, _cram, file) = cram_alignment_file(&records, permissive_config());

        assert!(drain(&file, "chr1", 50, 60).is_empty());
        assert_eq!(pool_len(&file), 1);
    }

    #[test]
    fn cram_pool_opens_once_and_returns_to_resting_size() {
        let records = [aln_record("a", 0, 1, 4, 60), aln_record("b", 0, 20, 4, 60)];
        let (_fa, _cram, file) = cram_alignment_file(&records, permissive_config());

        assert_eq!(pool_len(&file), 0);
        for _ in 0..5 {
            {
                let reads = file.get_reads_from_segment("chr1", 1, 200).unwrap();
                assert_eq!(pool_len(&file), 0);
                assert_eq!(reads.count(), 2);
            }
            assert_eq!(pool_len(&file), 1);
        }
    }

    #[test]
    fn cram_parallel_segments_match_sequential() {
        let records: Vec<RecordBuf> = (0..20)
            .map(|i| aln_record(&format!("r{i}"), 0, 1 + i * 5, 4, 60))
            .collect();
        let (_fa, _cram, file) = cram_alignment_file(&records, permissive_config());

        let segments: Vec<(u32, u32)> = (0..20).map(|i| (1 + i * 5, 10 + i * 5)).collect();
        let sequential: Vec<Vec<(u64, String)>> = segments
            .iter()
            .map(|&(s, e)| drain(&file, "chr1", s, e))
            .collect();
        let parallel: Vec<Vec<(u64, String)>> = segments
            .par_iter()
            .map(|&(s, e)| drain(&file, "chr1", s, e))
            .collect();

        assert_eq!(parallel, sequential);
    }

    #[test]
    fn cram_without_repository_is_rejected() {
        // A CRAM input with no reference repository is a programmer
        // error surfaced as a typed error, not a panic.
        let contigs = [ContigSpec {
            name: "chr1".into(),
            length: CONTIG_LEN as u64,
        }];
        let (_fa, fasta_path) = build_fasta(&contigs).expect("build fasta");
        let (_cram, cram_path) = build_cram(
            &fasta_path,
            &contigs,
            &HeaderOverrides::default(),
            &[aln_record("a", 0, 1, 4, 60)],
        )
        .expect("build cram");
        let repository = build_fasta_repository(&fasta_path).expect("repository");
        let (_reader, header) =
            open_cram_reader_with_header(&cram_path, Some(&repository)).expect("read header");
        let index = cram::fs::index(&cram_path).expect("build crai");

        assert!(matches!(
            AlignmentFile::from_input(
                cram_path,
                Arc::new(header),
                AlignmentIndex::Crai(Arc::new(index)),
                None,
                permissive_config(),
                0,
            ),
            Err(AlignmentInputError::MissingCramReference { .. })
        ));
    }

    #[test]
    fn from_input_rejects_format_index_mismatch() {
        // A BAM file paired with a CRAI index — driver-bug mismatch.
        let records = [aln_record("a", 0, 1, 4, 60)];
        let header = header_two_contigs();
        let (_dir, bam_path) = write_bam(&records, &header);
        let empty_crai = cram::crai::Index::default();

        assert!(matches!(
            AlignmentFile::from_input(
                bam_path,
                Arc::new(header),
                AlignmentIndex::Crai(Arc::new(empty_crai)),
                None,
                permissive_config(),
                0,
            ),
            Err(AlignmentInputError::AlignmentIndexFormatMismatch { .. })
        ));
    }
}
