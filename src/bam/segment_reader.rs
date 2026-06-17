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
// public surface (`get_reads_from_segment`, the CRAM path, the dispatch
// enum) is reachable only from `#[cfg(test)]`, which reads as dead code
// to the non-test lib build. Suppress that — but only in the non-test
// build, so the test build still flags any genuinely-dead helper added
// later. Remove the attribute entirely when the consumer lands.
#![cfg_attr(not(test), allow(dead_code))]

use std::collections::HashMap;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, SeekFrom};
use std::path::{Path, PathBuf};
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
    AlignmentMergedReaderConfig, ContigInterval, DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH,
    FilterCounts, MappedRead, PreDecodeOutcome, classify_pre_decode, record_buf_to_mapped_read,
};
use super::bam_input::{BamIndex, open_bam_reader_with_header, query_interval};
use super::cram_input::open_cram_reader_with_header;
use super::errors::AlignmentInputError;
use super::index_preflight::{AlignmentFileKind, AlignmentIndex};

// ---------------------------------------------------------------------
// Read filter (the cheap, pre-decode subset)
// ---------------------------------------------------------------------

/// The read filters the segment reader applies — deliberately only the
/// ones that are **cheap to evaluate before fully using a record**: SAM
/// flags (duplicate / QC-fail / secondary / supplementary / unmapped),
/// MAPQ, and decoded read length.
///
/// This is narrower than [`AlignmentMergedReaderConfig`] on purpose. The
/// expensive, reference-dependent filters that type also carries —
/// mismatch fraction (`max_read_mismatch_fraction` / `mismatch_bq_floor`),
/// bad-CIGAR, indel left-alignment — are intentionally **not**
/// representable here. A segment read is yielded whole and the consumer
/// applies those itself: the SSR fetcher realigns every spanning read
/// (so it distrusts the mapper's CIGAR that those filters lean on), and
/// the SNP `--regions` path already applies them in its parallel
/// read-processing stage. Making the reference-dependent fields
/// unrepresentable means a caller cannot set one and have it silently
/// ignored.
///
/// `min_mapq` / `min_read_length` are `Option`s: `None` is the explicit
/// "no minimum" state, distinct from any specific threshold.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SegmentReadFilter {
    /// `None` = no minimum; `Some(n)` = drop reads with MAPQ < n.
    pub min_mapq: Option<u8>,
    /// `None` = no minimum; `Some(n)` = drop reads with decoded SEQ
    /// length < n.
    pub min_read_length: Option<u32>,
    /// Drop reads with the QC-fail flag (`0x200`) set.
    pub drop_qc_fail: bool,
    /// Drop reads with the duplicate flag (`0x400`) set.
    pub drop_duplicate: bool,
}

impl Default for SegmentReadFilter {
    fn default() -> Self {
        Self {
            min_mapq: Some(DEFAULT_MIN_MAPQ),
            min_read_length: Some(DEFAULT_MIN_READ_LENGTH),
            drop_qc_fail: true,
            drop_duplicate: true,
        }
    }
}

impl From<&AlignmentMergedReaderConfig> for SegmentReadFilter {
    /// Project the merged reader's full config onto the cheap subset the
    /// segment reader applies. The reference-dependent fields
    /// (`max_read_mismatch_fraction` / `mismatch_bq_floor`) are dropped —
    /// the segment reader never applies them; the SNP `--regions` path's
    /// downstream read-processing stage does (and the SSR fetcher realigns).
    fn from(config: &AlignmentMergedReaderConfig) -> Self {
        Self {
            min_mapq: config.min_mapq,
            min_read_length: config.min_read_length,
            drop_qc_fail: config.drop_qc_fail,
            drop_duplicate: config.drop_duplicate,
        }
    }
}

impl SegmentReadFilter {
    /// The flag/MAPQ view consumed by
    /// [`classify_pre_decode`](super::alignment_input::classify_pre_decode),
    /// reusing the merged reader's single source of truth for the
    /// flag-drop policy. The reference-dependent fields are forced
    /// inert — the segment reader never applies them (see the type
    /// docs). A new field on [`AlignmentMergedReaderConfig`] breaks this
    /// literal, forcing a deliberate decision about whether it belongs
    /// in the cheap-filter set.
    fn pre_decode_config(&self) -> AlignmentMergedReaderConfig {
        AlignmentMergedReaderConfig {
            min_mapq: self.min_mapq,
            min_read_length: self.min_read_length,
            drop_qc_fail: self.drop_qc_fail,
            drop_duplicate: self.drop_duplicate,
            max_read_mismatch_fraction: None,
            mismatch_bq_floor: 0,
        }
    }
}

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
    /// `index`, the cheap-read `filter`, and the `source_file_index`
    /// stamped into every [`MappedRead`] this file yields (the caller's
    /// position for this file in its input list).
    ///
    /// `repository` is the FASTA reference; **required for CRAM**
    /// (slice decoding consults it) and ignored for BAM. The caller
    /// builds it once via
    /// [`crate::bam::alignment_input::build_fasta_repository`] and
    /// shares the clone (an `Arc` bump).
    ///
    /// # Filtering
    ///
    /// Only the cheap, pre-decode filters in [`SegmentReadFilter`] are
    /// applied (flags, MAPQ, read length). Reads are otherwise yielded
    /// whole; the reference-dependent filters (mismatch fraction,
    /// bad-CIGAR, indel left-alignment) are the consumer's job — see the
    /// [`SegmentReadFilter`] docs for why.
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
        filter: SegmentReadFilter,
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
                filter,
                source_file_index,
                readers_pool: Mutex::new(Vec::new()),
            })),
            (AlignmentFileKind::Bam, AlignmentIndex::BamBai(idx)) => Ok(Self::Bam(BamFile {
                path,
                index: BamIndex::Bai(idx),
                ref_map,
                filter,
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
                    filter,
                    source_file_index,
                    readers_pool: Mutex::new(Vec::new()),
                }))
            }
            // Genuine format / index disagreements — the driver's index
            // loader and the file extension say different things. These
            // are enumerated explicitly (rather than caught by a `_`
            // wildcard) so that adding a new `AlignmentIndex` variant —
            // the enum is `#[non_exhaustive]` — forces a new arm here
            // instead of being silently absorbed as a "format mismatch".
            // Mirrors `AlignmentMergedReader::query`'s convention.
            (
                AlignmentFileKind::Cram,
                index @ (AlignmentIndex::BamCsi(_) | AlignmentIndex::BamBai(_)),
            )
            | (AlignmentFileKind::Bam, index @ AlignmentIndex::Crai(_)) => {
                Err(AlignmentInputError::AlignmentIndexFormatMismatch {
                    path,
                    file_format: kind.display_name(),
                    index_format: index.display_name(),
                })
            }
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
    /// Reads failing the cheap [`SegmentReadFilter`] (flags / MAPQ /
    /// length) are dropped; everything else is yielded whole, with the
    /// reference-dependent filters left to the consumer.
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

    /// The input file's path — used by the multi-file merge for error
    /// messages (which file an out-of-order or duplicate read came from).
    pub(crate) fn path(&self) -> &Path {
        match self {
            Self::Bam(file) => &file.path,
            Self::Cram(file) => &file.path,
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

/// Apply the per-record cascade shared by both formats once a record has
/// been decoded: target-contig check, segment overlap, the cheap
/// pre-decode flag/MAPQ filter, and the min-read-length filter — then
/// convert to [`MappedRead`]. Returns `Ok(Some(read))` to keep,
/// `Ok(None)` to drop, and `Err` for a malformed record.
///
/// This is the single source of truth for filter membership and ordering
/// across the BAM and CRAM iterators, so a change to either applies to
/// both. The BAM-only sorted early-stop (`alignment_start > segment.end`)
/// stays in the BAM loop — it ends the scan rather than dropping one
/// record, and the CRAM path has no per-record sort guarantee to exploit.
///
/// `counts` accumulates the **filter** drops (flags / MAPQ / length) so a
/// consumer can report them (the SSR fetcher's `n_filtered` QC column).
/// Out-of-segment drops — wrong contig, or a footprint the bin-granular
/// index over-returns at the edges — are *not* counted: they are not reads
/// "at the locus", just chunk slop.
fn classify_segment_record(
    record: &sam::alignment::RecordBuf,
    target_reference_sequence_id: usize,
    segment: &ContigInterval,
    filter: &SegmentReadFilter,
    counts: &mut FilterCounts,
    source_file_index: usize,
    path: &Path,
) -> Result<Option<MappedRead>, AlignmentInputError> {
    // Different contig (a chunk/container straddling contigs) — drop,
    // uncounted (not a read at this segment).
    if record.reference_sequence_id() != Some(target_reference_sequence_id) {
        return Ok(None);
    }
    // Footprint does not overlap the segment — drop, uncounted (index
    // over-return at the chunk edges).
    if !segment.overlaps_record(record) {
        return Ok(None);
    }
    // Cheap pre-decode flag/MAPQ filter — a counted filter drop.
    if let PreDecodeOutcome::Drop(bucket) = classify_pre_decode(&filter.pre_decode_config(), record)
    {
        counts.record_drop(bucket);
        return Ok(None);
    }
    // Min read length — a counted filter drop.
    if let Some(min) = filter.min_read_length
        && (record.sequence().as_ref().len() as u32) < min
    {
        counts.too_short += 1;
        return Ok(None);
    }
    let mapped = record_buf_to_mapped_read(record, source_file_index).map_err(|source| {
        AlignmentInputError::MalformedRecord {
            path: path.to_path_buf(),
            qname: record
                .name()
                .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned()),
            source,
        }
    })?;
    Ok(Some(mapped))
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
    filter: SegmentReadFilter,
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
            filter_counts: FilterCounts::default(),
        })
    }

    /// Pop an idle reader from the pool, or open a fresh one if the
    /// pool is empty (the only file-open; happens ~once per concurrent
    /// caller).
    fn borrow_handle(&self) -> Result<BamReaderHandle, AlignmentInputError> {
        if let Some(handle) = self.lock_pool().pop() {
            return Ok(handle);
        }
        let (reader, _header) = open_bam_reader_with_header(&self.path)?;
        Ok(BamReaderHandle { reader })
    }

    fn return_handle(&self, handle: BamReaderHandle) {
        self.lock_pool().push(handle);
    }

    /// Lock the pool, recovering the guard if a previous panic poisoned
    /// the mutex. The lock only ever guards `Vec` pop/push/len — none of
    /// which can unwind — so poisoning cannot happen through this type's
    /// own code; recovering (rather than `expect`-panicking on borrow and
    /// silently dropping the reader on return) keeps a poison introduced
    /// elsewhere from cascading across every later query.
    fn lock_pool(&self) -> std::sync::MutexGuard<'_, Vec<BamReaderHandle>> {
        self.readers_pool
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    #[cfg(test)]
    fn pool_len(&self) -> usize {
        self.lock_pool().len()
    }
}

/// Per-call BAM iterator. `(target, segment)` are fixed at creation;
/// the only mutable state is the chunk cursor and the reader handle.
/// On drop it returns the borrowed reader to its [`BamFile`]'s pool.
pub(crate) struct BamSegmentReads<'a> {
    file: &'a BamFile,
    /// The borrowed reader. `Some` for the whole iterator lifetime;
    /// taken only in `Drop` to return it to the pool. `next` therefore
    /// always finds it present — the `.expect("handle held")` accesses
    /// below cannot fire, because `next` never runs during `Drop`.
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
    /// Running tally of reads dropped by the cheap filter (flags / MAPQ /
    /// length). Read via [`Self::filter_counts`] after draining.
    filter_counts: FilterCounts,
}

impl BamSegmentReads<'_> {
    fn io_error(&self, source: io::Error) -> AlignmentInputError {
        AlignmentInputError::Io {
            path: self.file.path.clone(),
            source,
        }
    }

    /// The filter-drop tally accumulated so far. See
    /// [`MappedReadsInSegment::filter_counts`].
    fn filter_counts(&self) -> &FilterCounts {
        &self.filter_counts
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
                    // BAM-only sorted early-stop: once a target-contig
                    // read starts past the segment end, nothing later in
                    // this linear scan can overlap, so end the iterator.
                    if record.reference_sequence_id() == Some(self.target_reference_sequence_id)
                        && let Some(astart) = record.alignment_start()
                        && usize::from(astart) as u64 > u64::from(self.segment.end)
                    {
                        self.done = true;
                        return None;
                    }
                    // Shared per-record cascade (contig / overlap / filter
                    // / convert).
                    match classify_segment_record(
                        &record,
                        self.target_reference_sequence_id,
                        &self.segment,
                        &self.file.filter,
                        &mut self.filter_counts,
                        self.file.source_file_index,
                        &self.file.path,
                    ) {
                        Ok(Some(mapped)) => return Some(Ok(mapped)),
                        Ok(None) => continue,
                        Err(e) => {
                            self.done = true;
                            return Some(Err(e));
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
    filter: SegmentReadFilter,
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
            filter_counts: FilterCounts::default(),
        })
    }

    fn borrow_handle(&self) -> Result<CramReaderHandle, AlignmentInputError> {
        if let Some(handle) = self.lock_pool().pop() {
            return Ok(handle);
        }
        let (reader, _header) = open_cram_reader_with_header(&self.path, Some(&self.repository))?;
        Ok(CramReaderHandle { reader })
    }

    fn return_handle(&self, handle: CramReaderHandle) {
        self.lock_pool().push(handle);
    }

    /// Lock the pool, recovering the guard on poison. See
    /// [`BamFile::lock_pool`] for the rationale (the lock guards only
    /// `Vec` operations, so a poison can only originate elsewhere, and
    /// recovering keeps it from cascading across later queries).
    fn lock_pool(&self) -> std::sync::MutexGuard<'_, Vec<CramReaderHandle>> {
        self.readers_pool
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    #[cfg(test)]
    fn pool_len(&self) -> usize {
        self.lock_pool().len()
    }
}

/// Per-call CRAM iterator. Walks the `.crai` for containers on the
/// target contig, seeks to each, decodes its records, and buffers the
/// segment-overlapping survivors in `pending`. On drop it returns the
/// borrowed reader to its [`CramFile`]'s pool.
pub(crate) struct CramSegmentReads<'a> {
    file: &'a CramFile,
    /// The borrowed reader. `Some` for the whole iterator lifetime;
    /// taken only in `Drop`. `refill` always finds it present — the
    /// `.expect("handle held")` access cannot fire, because `refill`
    /// never runs during `Drop`.
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
    /// Running tally of reads dropped by the cheap filter (flags / MAPQ /
    /// length). Read via [`Self::filter_counts`] after draining.
    filter_counts: FilterCounts,
}

/// Seek to `offset`, read one CRAM container, and decode all of its slices into
/// owned records. Returns `Ok(None)` at end-of-stream (`read_container` reads 0
/// — the EOF marker) and `Ok(Some(records))` for a decoded container (`records`
/// may be empty). Shared by the per-call [`CramSegmentReads::refill`] and the
/// decode-once [`CachingCramReader`] so the seek/decode path — and its EOF stop
/// semantics — cannot drift between the two readers.
fn decode_cram_container(
    reader: &mut cram::io::Reader<File>,
    repository: &fasta::Repository,
    header: &sam::Header,
    path: &Path,
    offset: u64,
) -> Result<Option<Vec<sam::alignment::RecordBuf>>, AlignmentInputError> {
    let to_io = |source| AlignmentInputError::Io {
        path: path.to_path_buf(),
        source,
    };
    reader.seek(SeekFrom::Start(offset)).map_err(&to_io)?;
    let mut container = cram::io::reader::Container::default();
    if reader.read_container(&mut container).map_err(&to_io)? == 0 {
        return Ok(None); // EOF — the index walk is exhausted.
    }
    let compression_header = container.compression_header().map_err(&to_io)?;
    let mut out = Vec::new();
    for slice_result in container.slices() {
        let slice = slice_result.map_err(&to_io)?;
        let (core_data_src, external_data_srcs) = slice.decode_blocks().map_err(&to_io)?;
        let cram_records = slice
            .records(
                repository.clone(),
                header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )
            .map_err(&to_io)?;
        for cram_record in &cram_records {
            out.push(
                sam::alignment::RecordBuf::try_from_alignment_record(header, cram_record)
                    .map_err(&to_io)?,
            );
        }
    }
    Ok(Some(out))
}

impl CramSegmentReads<'_> {
    /// The filter-drop tally accumulated so far. See
    /// [`MappedReadsInSegment::filter_counts`].
    fn filter_counts(&self) -> &FilterCounts {
        &self.filter_counts
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

            // Container-level narrowing using the coordinate-ordered
            // `.crai`: skip a container entirely before the segment, and
            // — crucially — *stop the whole walk* once a container starts
            // past the segment end, since every later container on this
            // contig starts no earlier and so cannot overlap. Without the
            // stop, a tiny locus near the start of a large contig would
            // scan that contig's entire index tail on every call.
            if let Some(container_start) = record.alignment_start() {
                let container_start = usize::from(container_start) as u64;
                if container_start > u64::from(self.segment.end) {
                    return Ok(true);
                }
                let span = record.alignment_span() as u64;
                if span > 0 {
                    let container_end = container_start + span - 1;
                    if container_end < u64::from(self.segment.start) {
                        continue;
                    }
                }
            }

            // Seek + decode this container via the shared helper, so the decode
            // path and its EOF stop semantics stay identical to the decode-once
            // `CachingCramReader`. The reader borrow is confined to this block.
            let records = {
                let reader = &mut self.handle.as_mut().expect("handle held").reader;
                decode_cram_container(
                    reader,
                    &self.file.repository,
                    self.file.ref_map.header(),
                    &self.file.path,
                    record.offset(),
                )?
            };
            let records = match records {
                None => return Ok(true), // EOF — index exhausted.
                Some(records) => records,
            };

            let mut converted: Vec<MappedRead> = Vec::new();
            for record_buf in &records {
                // Shared per-record cascade (contig / overlap / filter
                // / convert) — identical to the BAM path.
                if let Some(mapped) = classify_segment_record(
                    record_buf,
                    self.target_reference_sequence_id,
                    &self.segment,
                    &self.file.filter,
                    &mut self.filter_counts,
                    self.file.source_file_index,
                    &self.file.path,
                )? {
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
// CRAM — per-thread caching reader (decode-once)
// ---------------------------------------------------------------------

/// Default number of decoded CRAM containers a [`CachingCramReader`] keeps.
/// A locus window overlaps ≤ 2 containers (1 normally, 2 at a boundary); a
/// third covers forward jitter. See `doc/devel/architecture/ssr_pileup_read_buffer.md`.
///
/// Intentionally a fixed constant, not a CLI/runtime knob: it is a pure
/// speed/memory tradeoff and **does not affect the `.ssr.psp` output** (decode
/// count varies with the cap, decoded records do not), so it is neither exposed
/// nor recorded in the artefact header.
pub(crate) const DEFAULT_MAX_CACHED_CONTAINERS: usize = 3;

/// A tiny FIFO cache of decoded CRAM containers, keyed by container file
/// offset. Each entry is the container's **unfiltered** decoded records (all
/// of them — filtering happens per window in [`CachingCramReader`]), shared by
/// `Arc<[_]>` so a served window can outlive an eviction. `get`/`insert` scan
/// linearly: intentional at this `capacity` (a handful of entries), where a
/// linear scan over a `VecDeque` beats a map's hashing + indirection.
struct ContainerCache {
    capacity: usize,
    entries: VecDeque<(u64, Arc<[sam::alignment::RecordBuf]>)>,
}

impl ContainerCache {
    fn new(capacity: usize) -> Self {
        Self {
            capacity: capacity.max(1),
            entries: VecDeque::new(),
        }
    }

    fn get(&self, offset: u64) -> Option<Arc<[sam::alignment::RecordBuf]>> {
        self.entries
            .iter()
            .find(|(o, _)| *o == offset)
            .map(|(_, v)| Arc::clone(v))
    }

    fn insert(&mut self, offset: u64, value: Arc<[sam::alignment::RecordBuf]>) {
        if self.entries.iter().any(|(o, _)| *o == offset) {
            return;
        }
        if self.entries.len() >= self.capacity {
            self.entries.pop_front();
        }
        self.entries.push_back((offset, value));
    }
}

/// A per-thread CRAM reader that decodes each container **once** and serves
/// many nearby locus windows from a small FIFO cache. A drop-in for repeated
/// [`CramFile::get_reads_from_segment`] over nearby windows (the SSR Stage-1
/// fetch pattern): it returns the same `(reads, FilterCounts)` a per-call query
/// would, reusing [`classify_segment_record`] so categorization is identical.
///
/// **Not shared between threads** — it owns its file handle and cache, so it
/// needs no locking. Construct one per worker (see
/// [`CramFile::caching_reader`]); it shares the index / header via `Arc`, plus a
/// clone-shared `repository` and the `Copy` `filter`, with the originating
/// [`CramFile`]. Design + rationale:
/// `doc/devel/architecture/ssr_pileup_read_buffer.md`.
pub(crate) struct CachingCramReader {
    /// Opened lazily on the first decode so construction is infallible (the
    /// driver builds one per worker without I/O). `None` until first use.
    reader: Option<cram::io::Reader<File>>,
    path: PathBuf,
    index: Arc<cram::crai::Index>,
    repository: fasta::Repository,
    ref_map: Arc<HeaderRefMap>,
    filter: SegmentReadFilter,
    source_file_index: usize,
    cache: ContainerCache,
}

impl CachingCramReader {
    /// Reads overlapping `[start, end]` on `chrom`, plus the **complete** filter
    /// tally over that window — byte-identical to
    /// [`CramFile::get_reads_from_segment`] drained to exhaustion. Decodes each
    /// overlapping container once (cache miss) or reuses it (cache hit).
    pub(crate) fn fetch_mapped_reads(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<(Vec<MappedRead>, FilterCounts), AlignmentInputError> {
        let (target, segment) = resolve_segment(&self.ref_map, chrom, start, end)?;
        let mut reads = Vec::new();
        let mut counts = FilterCounts::default();

        // Walk the `.crai` exactly as `CramSegmentReads::refill` does: skip
        // off-contig / before-segment containers, and stop the whole walk once a
        // container starts past the segment end (coordinate-ordered `.crai`).
        for idx in 0..self.index.as_slice().len() {
            let record = self.index.as_slice()[idx].clone();
            if record.reference_sequence_id() != Some(target) {
                continue;
            }
            if let Some(container_start) = record.alignment_start() {
                let container_start = usize::from(container_start) as u64;
                if container_start > u64::from(segment.end) {
                    break;
                }
                let span = record.alignment_span() as u64;
                if span > 0 {
                    let container_end = container_start + span - 1;
                    if container_end < u64::from(segment.start) {
                        continue;
                    }
                }
            }

            let decoded = self.get_or_decode(record.offset())?;
            if decoded.is_empty() {
                // Empty decode == end-of-stream (`read_container` read 0) — stop
                // the walk, matching `CramSegmentReads::refill`'s `Ok(true)`. (A
                // decoded-but-record-less container is not a real CRAM shape, so
                // this never short-stops a live container.)
                break;
            }
            for record_buf in decoded.iter() {
                if let Some(mapped) = classify_segment_record(
                    record_buf,
                    target,
                    &segment,
                    &self.filter,
                    &mut counts,
                    self.source_file_index,
                    &self.path,
                )? {
                    reads.push(mapped);
                }
            }
        }

        Ok((reads, counts))
    }

    /// The decoded records of the container at `offset`, from cache or by
    /// decoding (and caching) it. An empty slice means end-of-stream (the shared
    /// helper read 0 — see [`decode_cram_container`]).
    fn get_or_decode(
        &mut self,
        offset: u64,
    ) -> Result<Arc<[sam::alignment::RecordBuf]>, AlignmentInputError> {
        if let Some(hit) = self.cache.get(offset) {
            return Ok(hit);
        }
        let decoded: Arc<[sam::alignment::RecordBuf]> =
            self.decode_container(offset)?.unwrap_or_default().into();
        self.cache.insert(offset, Arc::clone(&decoded));
        Ok(decoded)
    }

    /// Open the underlying reader if not yet open (lazy first-use), so
    /// construction stays infallible and I/O-free.
    fn ensure_open(&mut self) -> Result<(), AlignmentInputError> {
        if self.reader.is_none() {
            let (reader, _header) =
                open_cram_reader_with_header(&self.path, Some(&self.repository))?;
            self.reader = Some(reader);
        }
        Ok(())
    }

    /// Decode the container at `offset` into owned `RecordBuf`s — unfiltered (the
    /// caller classifies per window). `Ok(None)` at end-of-stream. Delegates to
    /// the shared [`decode_cram_container`] so the per-call and decode-once paths
    /// cannot drift.
    fn decode_container(
        &mut self,
        offset: u64,
    ) -> Result<Option<Vec<sam::alignment::RecordBuf>>, AlignmentInputError> {
        self.ensure_open()?;
        // PANIC-FREE: `ensure_open` populated `self.reader` immediately above.
        let reader = self.reader.as_mut().expect("opened by ensure_open");
        decode_cram_container(
            reader,
            &self.repository,
            self.ref_map.header(),
            &self.path,
            offset,
        )
    }
}

impl CramFile {
    /// Build a per-thread [`CachingCramReader`] over this file: its own open
    /// handle + a FIFO container cache, sharing the index / header / repository
    /// / filter via `Arc`. One per worker (not pooled).
    pub(crate) fn caching_reader(&self, max_cached_containers: usize) -> CachingCramReader {
        CachingCramReader {
            reader: None,
            path: self.path.clone(),
            index: Arc::clone(&self.index),
            repository: self.repository.clone(),
            ref_map: Arc::clone(&self.ref_map),
            filter: self.filter,
            source_file_index: self.source_file_index,
            cache: ContainerCache::new(max_cached_containers),
        }
    }
}

/// A per-worker read source: a caching CRAM reader (decode-once), or the shared
/// pooled path for BAM (lighter decode, no slice concept). Built by
/// [`AlignmentFile::worker_reader`] — one per worker via the driver's
/// `map_init` — and unifies fetching into `(reads, FilterCounts)`.
pub(crate) enum WorkerReader<'a> {
    Cram(CachingCramReader),
    /// BAM keeps the existing pooled per-call path; the borrow ties the reader
    /// to the shared [`AlignmentFile`] (which outlives the worker pool).
    Bam(&'a AlignmentFile),
}

impl WorkerReader<'_> {
    /// Reads overlapping `[start, end]` on `chrom` + the complete filter tally —
    /// the same `(reads, FilterCounts)` for both formats. For CRAM this hits the
    /// per-worker container cache; for BAM it drains the pooled iterator.
    pub(crate) fn fetch_mapped_reads(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<(Vec<MappedRead>, FilterCounts), AlignmentInputError> {
        match self {
            Self::Cram(reader) => reader.fetch_mapped_reads(chrom, start, end),
            Self::Bam(file) => {
                let mut iter = file.get_reads_from_segment(chrom, start, end)?;
                let mut reads = Vec::new();
                for read in iter.by_ref() {
                    reads.push(read?);
                }
                let counts = *iter.filter_counts();
                Ok((reads, counts))
            }
        }
    }
}

impl AlignmentFile {
    /// Build a per-worker [`WorkerReader`]: a caching CRAM reader (own handle +
    /// cache) for CRAM, or a borrow of `self` for BAM. Infallible (CRAM opens
    /// lazily on first decode). One per worker, via the driver's `map_init`.
    pub(crate) fn worker_reader(&self) -> WorkerReader<'_> {
        match self {
            Self::Cram(file) => {
                WorkerReader::Cram(file.caching_reader(DEFAULT_MAX_CACHED_CONTAINERS))
            }
            Self::Bam(_) => WorkerReader::Bam(self),
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

impl MappedReadsInSegment<'_> {
    /// The reads dropped by the cheap [`SegmentReadFilter`] (flags / MAPQ
    /// / length) during the scan so far, bucketed by reason. Reads dropped
    /// as out-of-segment (wrong contig, or index over-return at the chunk
    /// edges) are **not** counted — only genuine filter rejections are.
    ///
    /// The tally grows as the iterator is polled, so read it **after
    /// draining** for the final per-segment totals. The SSR fetcher uses
    /// it for the `n_filtered` QC column (the reader owns the filter, so it
    /// owns the count).
    pub(crate) fn filter_counts(&self) -> &FilterCounts {
        match self {
            Self::Bam(reads) => reads.filter_counts(),
            Self::Cram(reads) => reads.filter_counts(),
        }
    }
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
    fn permissive_config() -> SegmentReadFilter {
        SegmentReadFilter {
            min_mapq: None,
            min_read_length: None,
            drop_qc_fail: false,
            drop_duplicate: false,
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

    fn assert_send<T: Send>() {}

    #[test]
    fn worker_reader_is_send() {
        // `map_init` hands each `WorkerReader` to a rayon worker, so it must be
        // `Send`. A future non-`Send` field would otherwise surface as a deep
        // generic error at the driver call site; pin it to this named test.
        assert_send::<WorkerReader<'static>>();
    }

    #[test]
    fn container_cache_evicts_oldest_over_capacity() {
        let mut cache = ContainerCache::new(3);
        // Entry for offset N holds N records, so length identifies the entry.
        let entry = |n: u64| -> Arc<[RecordBuf]> {
            (0..n).map(|_| RecordBuf::default()).collect()
        };
        for off in 1..=4u64 {
            cache.insert(off, entry(off));
        }
        assert!(cache.get(1).is_none(), "oldest (offset 1) must be evicted");
        assert_eq!(cache.get(2).map(|v| v.len()), Some(2));
        assert_eq!(cache.get(3).map(|v| v.len()), Some(3));
        assert_eq!(cache.get(4).map(|v| v.len()), Some(4));
    }

    #[test]
    fn container_cache_insert_is_idempotent_on_repeat_offset() {
        let mut cache = ContainerCache::new(3);
        let first: Arc<[RecordBuf]> = Vec::new().into(); // 0 records
        let second: Arc<[RecordBuf]> = vec![RecordBuf::default()].into(); // 1 record
        cache.insert(7, Arc::clone(&first));
        cache.insert(7, Arc::clone(&second)); // ignored — offset already cached
        let got = cache.get(7).expect("offset 7 present");
        assert_eq!(got.len(), 0, "the first insert wins; the repeat is dropped");
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
        filter: SegmentReadFilter,
    ) -> (TempDir, AlignmentFile) {
        let header = header_two_contigs();
        let (dir, bam_path) = write_bam(records, &header);
        let csi = build_csi_in_memory(&bam_path);
        let file = AlignmentFile::from_input(
            bam_path,
            Arc::new(header),
            AlignmentIndex::BamCsi(Arc::new(csi)),
            None,
            filter,
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
    fn bam_filter_drops_are_counted_by_bucket_and_out_of_segment_reads_are_not() {
        use noodles_sam::alignment::record::Flags;
        let good_len = (DEFAULT_MIN_READ_LENGTH as usize) + 4; // 34, above the floor

        // A duplicate-flagged read (otherwise admissible) — drops as
        // `duplicate` under the default filter.
        let dup = RecordBuf::builder()
            .set_name(b"dup")
            .set_reference_sequence_id(0)
            .set_flags(Flags::DUPLICATE)
            .set_mapping_quality(MappingQuality::new(60).expect("mapq"))
            .set_alignment_start(Position::try_from(10).unwrap())
            .set_cigar([Op::new(Kind::Match, good_len)].into_iter().collect())
            .set_sequence(Sequence::from(vec![b'A'; good_len]))
            .set_quality_scores(QualityScores::from(vec![30u8; good_len]))
            .build();

        // All footprints start at 10 (overlap the [1,50] query) except
        // `faraway`, which sits past the segment end so the sorted
        // early-stop should drop it *before* the filter ever sees it.
        let records = [
            aln_record("good", 0, 10, good_len, 60),    // kept
            aln_record("lowmq", 0, 10, good_len, 5),    // low MAPQ (5 < 20)
            aln_record("short", 0, 10, 10, 60),         // SEQ len 10 < 30
            dup,                                        // duplicate flag
            aln_record("faraway", 0, 150, good_len, 5), // low MAPQ but out of segment
        ];
        // Default filter (not the permissive one) so the cascade bites.
        let (_dir, file) = bam_alignment_file(&records, SegmentReadFilter::default());

        let mut reads = file.get_reads_from_segment("chr1", 1, 50).expect("segment");
        let kept: Vec<String> = reads
            .by_ref()
            .map(|r| String::from_utf8_lossy(&r.expect("decode").qname).into_owned())
            .collect();
        let counts = reads.filter_counts();

        assert_eq!(kept, vec!["good".to_string()]);
        assert_eq!(counts.low_mapq, 1);
        assert_eq!(counts.too_short, 1);
        assert_eq!(counts.duplicate, 1);
        // `faraway` is out of segment (sorted early-stop) → never counted,
        // even though it would have failed the MAPQ filter.
        assert_eq!(counts.qc_fail, 0);
        assert_eq!(counts.secondary, 0);
        assert_eq!(counts.supplementary, 0);
        assert_eq!(counts.unmapped, 0);
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
    fn bam_pool_survives_a_poisoned_lock() {
        let records = [aln_record("a", 0, 10, 4, 60), aln_record("b", 0, 20, 4, 60)];
        let (_dir, file) = bam_alignment_file(&records, permissive_config());

        // Poison the pool mutex by panicking while holding its guard.
        let AlignmentFile::Bam(bam) = &file else {
            unreachable!("built a BAM file");
        };
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _guard = bam.readers_pool.lock().unwrap();
            panic!("poison the pool");
        }));
        assert!(bam.readers_pool.is_poisoned(), "pool should be poisoned");

        // A poisoned pool must not take down queries: borrow + return
        // recover the guard rather than panicking (the old `expect`) or
        // silently dropping the reader (the old `if let Ok`).
        assert_eq!(drain(&file, "chr1", 1, 200).len(), 2);
        assert_eq!(pool_len(&file), 1);
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
        let filter = SegmentReadFilter {
            min_mapq: Some(20),
            ..permissive_config()
        };
        let (_dir, file) = bam_alignment_file(&records, filter);

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
        filter: SegmentReadFilter,
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
            filter,
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

    /// The byte-identity gate for [`CachingCramReader`]: for a spread of
    /// repeated and overlapping windows, the caching reader must return exactly
    /// the same `(reads, FilterCounts)` as draining the per-call
    /// `get_reads_from_segment` path. Includes a low-MAPQ read so the filter
    /// tally is non-trivial, and repeats `8..25` so a cache hit is exercised.
    #[test]
    fn caching_cram_reader_matches_per_call_path() {
        let records = [
            aln_record("a", 0, 1, 4, 60),
            aln_record("b", 0, 6, 4, 60),
            aln_record("lowmq", 0, 8, 4, 5),
            aln_record("c", 0, 8, 4, 60),
            aln_record("d", 0, 20, 4, 60),
            aln_record("e", 0, 40, 4, 60),
        ];
        let filter = SegmentReadFilter {
            min_mapq: Some(20),
            min_read_length: None,
            drop_qc_fail: true,
            drop_duplicate: true,
        };
        let (_fa, _cram, file) = cram_alignment_file(&records, filter);
        let cram = match &file {
            AlignmentFile::Cram(c) => c,
            _ => unreachable!("cram_alignment_file builds a CRAM file"),
        };
        let mut caching = cram.caching_reader(DEFAULT_MAX_CACHED_CONTAINERS);

        for (s, e) in [
            (1u32, 50),
            (8, 25),
            (6, 9),
            (8, 11),
            (20, 23),
            (100, 110),
            (40, 43),
            (8, 25),
        ] {
            let mut iter = file.get_reads_from_segment("chr1", s, e).unwrap();
            let mut want_reads = Vec::new();
            for r in iter.by_ref() {
                want_reads.push(r.unwrap());
            }
            let want_counts = *iter.filter_counts();
            drop(iter);

            let (got_reads, got_counts) = caching.fetch_mapped_reads("chr1", s, e).unwrap();
            assert_eq!(got_reads, want_reads, "reads differ for window {s}..{e}");
            assert_eq!(got_counts, want_counts, "counts differ for window {s}..{e}");
        }
    }

    /// Exercises the container-level early-stop: when the container the
    /// `.crai` points at starts past the segment end, the walk stops and
    /// yields nothing — without dropping a container that *does* overlap.
    /// (The fixture writer packs these few reads into one container, so
    /// this checks the stop fires correctly on the observable single-
    /// container case; the multi-container stop is the same predicate.)
    #[test]
    fn cram_early_stops_when_container_starts_past_segment() {
        // All reads start at 100+ → the only container starts at 100.
        let records = [
            aln_record("a", 0, 100, 4, 60),
            aln_record("b", 0, 110, 4, 60),
        ];
        let (_fa, _cram, file) = cram_alignment_file(&records, permissive_config());

        // Segment entirely before the container start → early-stop, empty.
        assert!(drain(&file, "chr1", 1, 50).is_empty());
        assert_eq!(pool_len(&file), 1);

        // Segment reaching the container is still served (stop not too eager).
        assert_eq!(
            drain(&file, "chr1", 100, 130),
            vec![(100, "a".into()), (110, "b".into())]
        );
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

    #[test]
    fn from_input_rejects_cram_with_bam_index() {
        // The other mismatch arm: a CRAM file paired with a BAM CSI index.
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

        assert!(matches!(
            AlignmentFile::from_input(
                cram_path,
                Arc::new(header),
                AlignmentIndex::BamCsi(Arc::new(noodles_csi::Index::default())),
                Some(repository),
                permissive_config(),
                0,
            ),
            Err(AlignmentInputError::AlignmentIndexFormatMismatch { .. })
        ));
    }
}
