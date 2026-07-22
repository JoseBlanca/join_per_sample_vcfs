//! Serving one region: the BAM and CRAM region-query record sources, and the
//! order guard over the filtered stream.
//!
//! Both sources are index-driven siblings of the whole-file
//! [`BamRecordSource`](crate::ng::read::filtering::BamRecordSource) /
//! [`CramRecordSource`](crate::ng::read::filtering::CramRecordSource): they
//! query the **already-parsed** index for candidate chunks, seek a pooled
//! reader, and drop records the index over-returned — chunk-edge slop, which is
//! a reader concern and so is dropped **uncounted**, never charged to a filter
//! drop reason. There is no new trait: BAM and CRAM are two containers for one
//! idea, not competing implementations to bake off
//! (`doc/devel/ng/arch/alignment_file.md` §4).
//!
//! Two things here fail *quietly* rather than loudly, which is why each is
//! built and committed on its own with an independent oracle:
//!
//! - **A missed chunk edge is a wrong genotype, not a crash.** The region
//!   query's oracle is the existing whole-file source: an indexed query must
//!   return exactly what a full linear scan filtered to the same region returns
//!   (spec §7, T5).
//! - **A guard that never fires looks exactly like a guard that works.** The
//!   order check is mutation-verified — removing it must let a planted
//!   regression through (T4a).
//!
//! The order guard's state lives in the per-region iterator, never on the
//! handle, so querying region B and then region A is a new forward scan rather
//! than a spurious regression (spec §3.2).

use std::fs::File;
use std::io;
use std::io::SeekFrom;
use std::iter::FusedIterator;
use std::path::Path;
use std::sync::Arc;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_csi::BinningIndex;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_fasta as fasta;
use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;

use crate::bam::alignment_input::MappedRead;
use crate::bam::index_preflight::AlignmentIndex;
use crate::ng::read::filtering::{NoodlesRawRecord, ReadFilterError, RecordSource};
use crate::ng::types::{ContigId, GenomePosition, GenomeRegion, Position};

use super::AlignmentFileError;

/// The reads of one region of a BAM, as a [`RecordSource`] the step-1 filter
/// can drive.
///
/// Index-driven: the already-parsed index gives candidate chunks, each is
/// seeked in turn, and records are read one at a time into the caller's reused
/// buffer. The index answers at *bin* granularity, so it over-returns at the
/// edges — those records are dropped **uncounted**, because they are not reads
/// the filter rejected, they are reads the index was never asked about.
///
/// Owns its reader rather than borrowing one, so the whole chain can be built
/// by value and the reader handed back to the pool when the stream ends
/// ([`Self::into_reader`]).
pub(crate) struct BamRegionSource<'a> {
    reader: bam::io::Reader<bgzf::io::Reader<File>>,
    /// Borrowed from the `AlignmentFile`: parsed once, never per query.
    header: &'a sam::Header,
    chunks: std::vec::IntoIter<Chunk>,
    /// Where the chunk being scanned ends, or `None` to take the next chunk.
    current_chunk_end: Option<bgzf::VirtualPosition>,
    /// The region's contig as a `ref_id`. Equal to `ContigId` by the open
    /// gate, which is what makes this a cast rather than a lookup.
    target_reference_sequence_id: usize,
    region: GenomeRegion,
    source_file_index: usize,
    /// Latched once the scan passes the region end or the chunks run out.
    done: bool,
}

/// What a region query resolves to before any reader is involved: which
/// contig, and which chunks of the file to scan.
///
/// Split out from the source so **every fallible part of setting up a query
/// happens before a reader is taken from the pool**. Taking the handle
/// transfers the obligation to return it; a `?` after that point would lose a
/// reader silently (`BorrowedReader::take`).
pub(crate) struct RegionPlan {
    chunks: Vec<Chunk>,
    target_reference_sequence_id: usize,
    region: GenomeRegion,
}

impl<'a> BamRegionSource<'a> {
    /// Resolve the region and query the index for its chunks.
    ///
    /// The index query happens here, once per region — an in-memory lookup on
    /// the index parsed at open, never a re-parse (spec §3.3). No reader is
    /// touched, so a failure costs nothing to recover from.
    pub(crate) fn plan(
        header: &sam::Header,
        index: &AlignmentIndex,
        region: GenomeRegion,
        path: &Path,
    ) -> Result<RegionPlan, AlignmentFileError> {
        let invalid_region = || AlignmentFileError::Region { region };

        let target_reference_sequence_id =
            usize::try_from(region.contig.get()).map_err(|_| invalid_region())?;
        if target_reference_sequence_id >= header.reference_sequences().len() || region.is_empty() {
            return Err(invalid_region());
        }

        let interval = interval_of(region).ok_or_else(invalid_region)?;
        let chunks = match index {
            AlignmentIndex::BamCsi(index) => index.query(target_reference_sequence_id, interval),
            AlignmentIndex::BamBai(index) => index.query(target_reference_sequence_id, interval),
            // `open` pairs a `.crai` only with a CRAM, and CRAM goes through
            // `CramRegionSource`. Reaching here means the handle was built
            // wrong, which is this module's bug and not the caller's — so it
            // must not masquerade as bad input.
            AlignmentIndex::Crai(_) => unreachable!(
                "a .crai index reached the BAM region source; \
                 open pairs .crai only with CRAM"
            ),
        }
        // The region was checked against the contig table above, so a failure
        // here is the index itself being unusable, not the query being silly.
        .map_err(|source| AlignmentFileError::IndexQuery {
            path: path.to_path_buf(),
            region,
            source,
        })?;

        Ok(RegionPlan {
            chunks,
            target_reference_sequence_id,
            region,
        })
    }

    /// Attach a reader to a plan. **Infallible**, which is what makes it safe
    /// to call after the pooled handle has been taken.
    pub(crate) fn new(
        reader: bam::io::Reader<bgzf::io::Reader<File>>,
        header: &'a sam::Header,
        plan: RegionPlan,
        source_file_index: usize,
    ) -> Self {
        Self {
            reader,
            header,
            chunks: plan.chunks.into_iter(),
            current_chunk_end: None,
            target_reference_sequence_id: plan.target_reference_sequence_id,
            region: plan.region,
            source_file_index,
            done: false,
        }
    }

    /// Give the reader back, for the pool.
    pub(crate) fn into_reader(self) -> bam::io::Reader<bgzf::io::Reader<File>> {
        self.reader
    }
}

/// The region as a noodles query interval. `None` if the 1-based bounds cannot
/// be represented (a `0` start, or a value past `usize`).
fn interval_of(region: GenomeRegion) -> Option<noodles_core::region::Interval> {
    let start = noodles_core::Position::new(usize::try_from(region.start.get()).ok()?)?;
    let end = noodles_core::Position::new(usize::try_from(region.end.get()).ok()?)?;
    Some((start..=end).into())
}

impl RecordSource for BamRegionSource<'_> {
    type Record = NoodlesRawRecord;

    fn header(&self) -> &sam::Header {
        self.header
    }

    fn read_next(&mut self, buf: &mut NoodlesRawRecord) -> io::Result<bool> {
        if self.done {
            return Ok(false);
        }

        loop {
            // 1. Land inside a chunk, seeking to its start if we have just
            //    taken it.
            let chunk_end = match self.current_chunk_end {
                Some(end) => end,
                None => {
                    let Some(chunk) = self.chunks.next() else {
                        self.done = true;
                        return Ok(false);
                    };
                    self.reader
                        .get_mut()
                        .seek(chunk.start())
                        .inspect_err(|_| self.done = true)?;
                    self.current_chunk_end = Some(chunk.end());
                    chunk.end()
                }
            };

            // 2. The previous read may have carried us past this chunk's end.
            if self.reader.get_ref().virtual_position() >= chunk_end {
                self.current_chunk_end = None;
                continue;
            }

            // 3. One record, into the caller's reused buffer.
            buf.source_file_index = self.source_file_index;
            let bytes_read = self
                .reader
                .read_record_buf(self.header, &mut buf.record)
                .inspect_err(|_| self.done = true)?;
            if bytes_read == 0 {
                // End of file inside a chunk — try the next one.
                self.current_chunk_end = None;
                continue;
            }

            let on_target_contig =
                buf.record.reference_sequence_id() == Some(self.target_reference_sequence_id);

            // 4. The sorted early stop: once a read *on the target contig*
            //    starts past the region, nothing later in this scan can
            //    overlap, because the file is coordinate-sorted and the open
            //    gate proved it claims to be.
            if on_target_contig
                && let Some(start) = buf.record.alignment_start()
                && usize::from(start) as u64 > self.region.end.get()
            {
                self.done = true;
                return Ok(false);
            }

            // 5. Drop what the index over-returned — a different contig in a
            //    chunk that straddles one, or a footprint that misses the
            //    region. Uncounted: these are a reader's business, not a
            //    filter's, and charging them to a `DropReason` would make the
            //    tally mean something different for an indexed read than for a
            //    whole-file one.
            if !on_target_contig || !overlaps(&buf.record, self.region) {
                continue;
            }

            return Ok(true);
        }
    }
}

/// Which `.crai` entries a CRAM region query has to walk, resolved before any
/// reader is involved. The sibling of [`RegionPlan`], and split out for the
/// same reason: nothing fallible may happen after a reader leaves the pool.
pub(crate) struct CramRegionPlan {
    /// Just this contig's entries, in file order — see
    /// [`crate::ng::read::input::open_bam::group_crai_by_contig`] for why the
    /// grouping is done once at open rather than searched per query.
    entries: Arc<[cram::crai::Record]>,
    target_reference_sequence_id: usize,
    region: GenomeRegion,
}

/// The reads of one region of a CRAM.
///
/// CRAM decodes at **container** granularity — a container holds slices whose
/// records are decoded together against the reference — so this walks the
/// `.crai` for containers on the target contig, decodes each, and buffers the
/// overlapping records. Non-overlapping ones are dropped **uncounted**, as in
/// the BAM source, and for the same reason.
///
/// **The `.crai` walk starts from a binary search, not from entry 0.**
/// Production rescans the index from the beginning on every call
/// (`segment_reader.rs`'s `fetch_mapped_reads`), which for a late contig in a
/// many-contig `.crai` is a repeated O(n) prefix scan — paid once per locus,
/// ~10⁶ times on the STR path. The `.crai` is sorted by
/// `(reference_sequence_id, alignment_start)`, so the contig's first entry is
/// one binary search away, and the forward walk from there is bounded by the
/// container-level early stop below (spec §3.3).
pub(crate) struct CramRegionSource<'a> {
    reader: cram::io::Reader<File>,
    header: &'a sam::Header,
    /// The reference bases slice decoding consults. Cheap to clone — it is
    /// internally shared — and cloned per decode, as noodles requires.
    repository: fasta::Repository,
    /// This contig's `.crai` entries, in file order.
    entries: Arc<[cram::crai::Record]>,
    /// Cursor into [`Self::entries`].
    next_index_record: usize,
    /// The overlapping records of the container last decoded, drained one per
    /// `read_next` before the next container is touched.
    pending: std::vec::IntoIter<RecordBuf>,
    /// The offset last decoded, so a container that appears once per slice is
    /// decoded once.
    last_decoded_offset: Option<u64>,
    target_reference_sequence_id: usize,
    region: GenomeRegion,
    source_file_index: usize,
    done: bool,
}

impl<'a> CramRegionSource<'a> {
    /// Pick out the target contig's `.crai` entries. **O(1)** — the grouping
    /// was done once, at open (`open_bam::group_crai_by_contig`, which also
    /// explains why this is not a binary search).
    pub(crate) fn plan(
        header: &sam::Header,
        crai_by_contig: &[Arc<[cram::crai::Record]>],
        region: GenomeRegion,
    ) -> Result<CramRegionPlan, AlignmentFileError> {
        let invalid_region = || AlignmentFileError::Region { region };

        let target_reference_sequence_id =
            usize::try_from(region.contig.get()).map_err(|_| invalid_region())?;
        if target_reference_sequence_id >= header.reference_sequences().len() || region.is_empty() {
            return Err(invalid_region());
        }

        let entries = crai_by_contig
            .get(target_reference_sequence_id)
            .cloned()
            .ok_or_else(invalid_region)?;

        Ok(CramRegionPlan {
            entries,
            target_reference_sequence_id,
            region,
        })
    }

    /// Attach a reader to a plan. **Infallible**, so it is safe after the
    /// pooled handle has been taken.
    pub(crate) fn new(
        reader: cram::io::Reader<File>,
        header: &'a sam::Header,
        repository: fasta::Repository,
        plan: CramRegionPlan,
        source_file_index: usize,
    ) -> Self {
        Self {
            reader,
            header,
            repository,
            entries: plan.entries,
            next_index_record: 0,
            pending: Vec::new().into_iter(),
            last_decoded_offset: None,
            target_reference_sequence_id: plan.target_reference_sequence_id,
            region: plan.region,
            source_file_index,
            done: false,
        }
    }

    pub(crate) fn into_reader(self) -> cram::io::Reader<File> {
        self.reader
    }

    /// Walk the `.crai` to the next container that could overlap, decode it,
    /// and buffer the records that do. `Ok(false)` once the walk is done.
    fn refill(&mut self) -> io::Result<bool> {
        loop {
            let Some(record) = self.entries.get(self.next_index_record).cloned() else {
                return Ok(false);
            };
            self.next_index_record += 1;

            debug_assert_eq!(
                record.reference_sequence_id(),
                Some(self.target_reference_sequence_id),
                "the grouping at open put only this contig's entries here"
            );

            // Container-level narrowing. Stopping once a container *starts*
            // past the region is what keeps a tiny locus near the start of a
            // large contig from walking that contig's whole index tail.
            if let Some(start) = record.alignment_start() {
                let container_start = usize::from(start) as u64;
                if container_start > self.region.end.get() {
                    return Ok(false);
                }
                let span = record.alignment_span() as u64;
                if span > 0 && container_start + span - 1 < self.region.start.get() {
                    continue;
                }
            }

            // `.crai` offsets are *container* positions and a container may
            // hold several slices, so a multi-slice container appears as
            // several entries sharing one offset. Decoding it once per entry
            // would surface every record again — caught loudly by the order
            // guard rather than silently inflating depth, but wrong either way.
            if self.last_decoded_offset == Some(record.offset()) {
                continue;
            }
            self.last_decoded_offset = Some(record.offset());

            let Some(records) = self.decode_container_at(record.offset())? else {
                // End of stream reached through the index — nothing further.
                return Ok(false);
            };

            let overlapping: Vec<RecordBuf> = records
                .into_iter()
                .filter(|record| {
                    record.reference_sequence_id() == Some(self.target_reference_sequence_id)
                        && overlaps(record, self.region)
                })
                .collect();

            if !overlapping.is_empty() {
                self.pending = overlapping.into_iter();
                return Ok(true);
            }
            // The container was in range but held nothing overlapping — keep
            // walking rather than reporting end of input.
        }
    }

    /// Seek to a container and decode every slice into owned records.
    /// `Ok(None)` at end of stream (`read_container` reads 0 — the EOF marker).
    fn decode_container_at(&mut self, offset: u64) -> io::Result<Option<Vec<RecordBuf>>> {
        self.reader.seek(SeekFrom::Start(offset))?;

        let mut container = cram::io::reader::Container::default();
        if self.reader.read_container(&mut container)? == 0 {
            return Ok(None);
        }

        let compression_header = container.compression_header()?;
        let mut records = Vec::new();
        for slice in container.slices() {
            let slice = slice?;
            // The decoded block data and the borrowed records live only within
            // this block; converting to owned `RecordBuf`s here keeps the
            // result independent of those borrows.
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
            for record in slice.records(
                self.repository.clone(),
                self.header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )? {
                records.push(RecordBuf::try_from_alignment_record(self.header, &record)?);
            }
        }
        Ok(Some(records))
    }
}

impl RecordSource for CramRegionSource<'_> {
    type Record = NoodlesRawRecord;

    fn header(&self) -> &sam::Header {
        self.header
    }

    fn read_next(&mut self, buf: &mut NoodlesRawRecord) -> io::Result<bool> {
        if self.done {
            return Ok(false);
        }

        buf.source_file_index = self.source_file_index;
        loop {
            if let Some(record) = self.pending.next() {
                // CRAM decodes a whole container at once, so this *moves* an
                // already-decoded record into the caller's buffer rather than
                // refilling it in place — the container buffer's allocations
                // are what get reused. Same seam, different reuse point.
                buf.record = record;
                return Ok(true);
            }

            match self.refill() {
                Ok(true) => {}
                Ok(false) => {
                    self.done = true;
                    return Ok(false);
                }
                Err(error) => {
                    self.done = true;
                    return Err(error);
                }
            }
        }
    }
}

/// Whether a record's reference footprint touches the region. Shared by both
/// sources so the two containers cannot disagree about what "overlap" means —
/// which is precisely what the BAM/CRAM parity oracle (T8) would otherwise be
/// left to discover.
fn overlaps(record: &sam::alignment::RecordBuf, region: GenomeRegion) -> bool {
    match (record.alignment_start(), record.alignment_end()) {
        (Some(first), Some(last)) => {
            usize::from(first) as u64 <= region.end.get()
                && usize::from(last) as u64 >= region.start.get()
        }
        _ => false,
    }
}

/// Proves, while streaming, that the reads really do arrive in genome order.
///
/// Everything downstream is built on that, so it is **checked, not trusted**:
/// the last emitted [`GenomePosition`] is kept, and a read whose key is
/// **strictly less** than it is a hard [`AlignmentFileError::OutOfOrderRead`].
/// There is no tolerance, no warning-and-continue, and no silent re-sort — an
/// unsorted file is a fatal input error the user must see.
///
/// Four things fix the semantics (spec §3.2):
///
/// - **The key is `(contig, position)`**, so a position regression within a
///   contig and a contig-order regression are one violation, not two checks.
///   "Contig order" means the reference's, because the open gate proved
///   `ref_id == ContigId`.
/// - **Equal keys are legal.** Several reads may start at the same base; only a
///   decrease is rejected.
/// - **It sits on the *filtered* stream**, so the guarantee is about exactly
///   the reads the caller sees. Dropped reads cannot break the monotonicity of
///   what survives.
/// - **The state lives here, in the per-region iterator, never on the handle.**
///   A caller is entitled to query region B and then region A; that is a new
///   forward scan, not a regression. Carrying the last key on the handle would
///   turn legitimate random access into a spurious error.
///
/// This and the `@HD SO` check at open are complementary rather than redundant:
/// the header check is cheap and fails before any work, while this one catches
/// the file that *claims* to be sorted and is not — which the header check
/// structurally cannot see.
pub(crate) struct OrderVerified<I> {
    inner: I,
    /// The last key emitted, or `None` before the first read.
    last: Option<GenomePosition>,
    path: Arc<Path>,
    /// Set on clean end of input or after an error is yielded, so the iterator
    /// is fused: the first `Err` surfaces once and then `None`.
    done: bool,
}

impl<I> OrderVerified<I> {
    pub(crate) fn new(inner: I, path: Arc<Path>) -> Self {
        Self {
            inner,
            last: None,
            path,
            done: false,
        }
    }

    /// Unwrap, so the stream that owns this can reclaim what it lent.
    pub(crate) fn into_inner(self) -> I {
        self.inner
    }
}

/// The region reader for whichever container the file is.
///
/// A [`RecordSource`] that delegates to the per-format source. BAM and CRAM are
/// two containers for one idea, so this is an enum rather than a trait — and it
/// is what lets everything above the source be written once, generic over
/// neither format nor container (`arch/alignment_file.md` §4).
pub(crate) enum RegionSource<'a> {
    Bam(BamRegionSource<'a>),
    Cram(CramRegionSource<'a>),
}

impl RegionSource<'_> {
    /// Give the reader back, in the shape the pool stores it.
    pub(crate) fn into_reader(self) -> super::open_bam::ReaderKind {
        match self {
            Self::Bam(source) => super::open_bam::ReaderKind::Bam(source.into_reader()),
            Self::Cram(source) => super::open_bam::ReaderKind::Cram(source.into_reader()),
        }
    }
}

impl RecordSource for RegionSource<'_> {
    type Record = NoodlesRawRecord;

    fn header(&self) -> &sam::Header {
        match self {
            Self::Bam(source) => source.header(),
            Self::Cram(source) => source.header(),
        }
    }

    fn read_next(&mut self, buf: &mut NoodlesRawRecord) -> io::Result<bool> {
        match self {
            Self::Bam(source) => source.read_next(buf),
            Self::Cram(source) => source.read_next(buf),
        }
    }
}

/// Where a read sits in the genome. Sound as a cross-file comparison key only
/// because the open gate proved this file's `ref_id`s are the reference's
/// `ContigId`s.
fn key_of(read: &MappedRead) -> GenomePosition {
    GenomePosition {
        // PANIC-FREE: `ref_id` comes from a 32-bit field in both BAM and CRAM,
        // so it fits by construction. Checked rather than `as`-cast because a
        // wrapped value would collapse two contigs onto one key and silently
        // disarm this guard — the one thing it must never do. Matches
        // `filtering.rs`'s handling of the same value.
        contig: ContigId(u32::try_from(read.ref_id).expect("ref_id fits u32")),
        position: Position(read.pos),
    }
}

impl<I> Iterator for OrderVerified<I>
where
    I: Iterator<Item = Result<MappedRead, ReadFilterError>>,
{
    type Item = Result<MappedRead, AlignmentFileError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        match self.inner.next() {
            None => {
                self.done = true;
                None
            }
            Some(Err(error)) => {
                self.done = true;
                Some(Err(AlignmentFileError::Filter(error)))
            }
            Some(Ok(read)) => {
                let current = key_of(&read);
                if let Some(previous) = self.last
                    && current < previous
                {
                    self.done = true;
                    return Some(Err(AlignmentFileError::OutOfOrderRead {
                        path: self.path.to_path_buf(),
                        previous,
                        current,
                    }));
                }
                self.last = Some(current);
                Some(Ok(read))
            }
        }
    }
}

impl<I> FusedIterator for OrderVerified<I> where
    I: Iterator<Item = Result<MappedRead, ReadFilterError>>
{
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
    use tempfile::TempDir;

    use super::*;
    use crate::bam::index_preflight::load_alignment_index;
    use crate::ng::read::filtering::BamRecordSource;
    use crate::ng::read::input::open_bam::group_crai_by_contig;
    use crate::ng::read::input::test_fixtures::{bam_header, indexed_bam, read_named};
    use crate::ng::types::{ContigId, Position};

    /// Two contigs long enough that the reads below span **several BGZF
    /// blocks**, so the index returns more than one chunk and the seek-per-chunk
    /// path is genuinely exercised. A fixture small enough to fit one block
    /// would make the oracle agree for the wrong reason.
    const CONTIGS: [(&str, usize); 2] = [("chr1", 100_000), ("chr2", 50_000)];

    fn contig_specs() -> Vec<(&'static str, usize, Option<&'static str>)> {
        CONTIGS.iter().map(|(n, l)| (*n, *l, None)).collect()
    }

    /// A coordinate-sorted spread of reads over both contigs: one every 20 bp,
    /// plus deliberate pile-ups so several reads share a start.
    fn spread_of_reads() -> Vec<RecordBuf> {
        let mut records = Vec::new();
        for (reference_sequence_id, (_, length)) in CONTIGS.iter().enumerate() {
            let mut start = 1;
            while start + 10 < *length {
                records.push(read_named(
                    &format!("r{reference_sequence_id}_{start}"),
                    reference_sequence_id,
                    start,
                ));
                // Every 50th read is doubled at the same start, so equal
                // positions are part of the fixture rather than an edge case
                // the oracle never sees.
                if start % 1000 == 1 {
                    records.push(read_named(
                        &format!("r{reference_sequence_id}_{start}_b"),
                        reference_sequence_id,
                        start,
                    ));
                }
                start += 20;
            }
        }

        // The shapes `overlaps_region`'s `_ => false` arm exists for, which a
        // spread of ordinary reads would never produce: a placed record with no
        // alignment start, and an entirely unplaced one. Both must be dropped,
        // and neither is a *filter* drop.
        records.push(unmapped_but_placed("unmapped_placed", 0));
        records.push(unplaced("unplaced"));
        records
    }

    /// A record naming a contig but carrying no alignment start — legal in SAM,
    /// and it has no footprint to overlap anything.
    fn unmapped_but_placed(qname: &str, reference_sequence_id: usize) -> RecordBuf {
        RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_reference_sequence_id(reference_sequence_id)
            .set_flags(noodles_sam::alignment::record::Flags::UNMAPPED)
            .set_sequence(Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(QualityScores::from(vec![30u8; 10]))
            .build()
    }

    /// A record on no contig at all — the unplaced reads that sit at a sorted
    /// file's tail.
    fn unplaced(qname: &str) -> RecordBuf {
        RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_flags(noodles_sam::alignment::record::Flags::UNMAPPED)
            .set_sequence(Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(QualityScores::from(vec![30u8; 10]))
            .build()
    }

    fn fixture() -> (TempDir, std::path::PathBuf, sam::Header) {
        let header = bam_header(&contig_specs());
        let (dir, path) = indexed_bam(&header, &spread_of_reads());
        (dir, path, header)
    }

    fn region(contig: u32, start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(end),
        }
    }

    /// **The oracle: a whole-file linear scan, filtered to the region.**
    ///
    /// Independent of the code under test in the dimension that matters here:
    /// it reaches every record through the whole-file `BamRecordSource` and
    /// never seeks, so nothing about chunks, seeking or the early stop is
    /// shared. The overlap rule is retyped rather than called, which catches a
    /// typo but not a shared misunderstanding of what "overlap" means — that
    /// part is pinned by the explicit boundary cases below instead.
    fn reads_by_linear_scan(
        path: &Path,
        header: &sam::Header,
        region: GenomeRegion,
    ) -> Vec<String> {
        let mut reader = bam::io::reader::Builder
            .build_from_path(path)
            .expect("open bam");
        reader.read_header().expect("read header");
        let mut source = BamRecordSource::new(reader, header.clone(), 0);

        let target = region.contig.get() as usize;
        let mut buf = NoodlesRawRecord::default();
        let mut names = Vec::new();

        while source.read_next(&mut buf).expect("read") {
            let record = &buf.record;
            if record.reference_sequence_id() != Some(target) {
                continue;
            }
            let (Some(first), Some(last)) = (record.alignment_start(), record.alignment_end())
            else {
                continue;
            };
            let overlaps = usize::from(first) as u64 <= region.end.get()
                && usize::from(last) as u64 >= region.start.get();
            if overlaps {
                names.push(read_name(record));
            }
        }
        names
    }

    /// The same question through the **indexed** path.
    fn reads_by_index(path: &Path, header: &sam::Header, region: GenomeRegion) -> Vec<String> {
        let index = load_alignment_index(path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(path)
            .expect("open bam");
        reader.read_header().expect("read header");

        let plan = BamRegionSource::plan(header, &index, region, path).expect("plan the query");
        let mut source = BamRegionSource::new(reader, header, plan, 0);

        let mut buf = NoodlesRawRecord::default();
        let mut names = Vec::new();
        while source.read_next(&mut buf).expect("read") {
            names.push(read_name(&buf.record));
        }
        names
    }

    fn read_name(record: &RecordBuf) -> String {
        String::from_utf8_lossy(record.name().expect("named").as_ref()).into_owned()
    }

    /// **T5 — the indexed query equals a full linear scan of the same region.**
    ///
    /// This is the step whose failure is silent: a chunk edge handled wrongly
    /// drops reads that were really there, and the caller gets a plausible
    /// answer with a wrong genotype rather than a crash. Nothing the indexed
    /// path could assert about *itself* would catch that, so it is checked
    /// against an independent implementation over a range of regions —
    /// interior, boundary, empty, whole-contig, and the second contig.
    #[test]
    fn t5_the_indexed_query_returns_exactly_what_a_linear_scan_returns() {
        let (_dir, path, header) = fixture();

        let regions = [
            region(0, 1, 100),          // the very start
            region(0, 1, 1),            // a single base
            region(0, 40_000, 40_100),  // deep in the interior
            region(0, 21, 21),          // exactly one read's start
            region(0, 30, 35),          // covered only by an overlapping read
            region(0, 1, 100_000),      // the whole contig
            region(0, 99_990, 100_000), // the far end
            region(1, 1, 500),          // the second contig
            region(1, 25_000, 25_050),  // its interior
            region(1, 49_000, 50_000),  // its far end
            region(0, 50_001, 50_001),  // a position with no read starting
        ];

        for region in regions {
            let expected = reads_by_linear_scan(&path, &header, region);
            let actual = reads_by_index(&path, &header, region);

            assert_eq!(
                actual,
                expected,
                "indexed query disagreed with the linear scan for contig {} [{}, {}]",
                region.contig.get(),
                region.start.get(),
                region.end.get()
            );
        }
    }

    /// The oracle is only worth something if it can fail. A region the fixture
    /// genuinely covers must return reads — otherwise T5 could be comparing two
    /// empty vectors and calling it agreement.
    #[test]
    fn the_oracle_is_not_vacuous() {
        let (_dir, path, header) = fixture();

        let interior = reads_by_linear_scan(&path, &header, region(0, 40_000, 40_100));
        assert!(
            interior.len() >= 5,
            "the fixture should cover this region several times over, got {}",
            interior.len()
        );

        let empty = reads_by_linear_scan(&path, &header, region(0, 99_999, 100_000));
        assert!(
            empty.len() < interior.len(),
            "and the regions must not all return the same thing"
        );
    }

    /// Records with no footprint — placed-but-unmapped, and unplaced — are
    /// dropped by the region query, and dropped **uncounted**: they are not
    /// reads the filter rejected, they are reads the index was never asked
    /// about. T5 alone cannot show this, because the oracle drops them by the
    /// same rule and the two would agree however the arm behaved.
    #[test]
    fn records_without_a_footprint_never_surface_from_a_region_query() {
        let (_dir, path, header) = fixture();

        for region in [region(0, 1, 100_000), region(1, 1, 50_000)] {
            let names = reads_by_index(&path, &header, region);
            assert!(
                !names
                    .iter()
                    .any(|n| n.starts_with("unmapped_placed") || n.starts_with("unplaced")),
                "a footprint-less record surfaced from {names:?}"
            );
        }
    }

    /// Reads are returned in coordinate order, which the whole chain above
    /// assumes and the merge one layer up depends on.
    #[test]
    fn the_indexed_query_returns_reads_in_coordinate_order() {
        let (_dir, path, header) = fixture();
        let index = load_alignment_index(&path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(&path)
            .expect("open bam");
        reader.read_header().expect("read header");

        let mut source = BamRegionSource::new(
            reader,
            &header,
            BamRegionSource::plan(&header, &index, region(0, 1, 100_000), &path).expect("plan"),
            0,
        );

        let mut buf = NoodlesRawRecord::default();
        let mut previous = 0u64;
        let mut seen = 0;
        while source.read_next(&mut buf).expect("read") {
            let start = usize::from(buf.record.alignment_start().expect("mapped")) as u64;
            assert!(start >= previous, "{start} came after {previous}");
            previous = start;
            seen += 1;
        }
        assert!(seen > 100, "the scan should have covered the contig");
    }

    /// **The early stop.** Once a read on the target contig starts past the
    /// region, the sort guarantees nothing later overlaps — so the scan must
    /// end rather than read on to the end of the file.
    ///
    /// The assertion that matters is the *unconsumed chunks*, not `done`:
    /// `done` latches both on the early stop and on running out of chunks, so
    /// a test that only checked it would pass with the early stop deleted. This
    /// failure changes no answers, only work — a whole-file scan per region —
    /// which is exactly why it needs its own discriminating assertion.
    #[test]
    fn the_scan_stops_once_it_passes_the_region() {
        let (_dir, path, header) = fixture();
        let index = load_alignment_index(&path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(&path)
            .expect("open bam");
        reader.read_header().expect("read header");

        let mut source = BamRegionSource::new(
            reader,
            &header,
            BamRegionSource::plan(&header, &index, region(0, 1, 100), &path).expect("plan"),
            0,
        );

        // Drain, then confirm the source latched `done` rather than running on.
        let mut buf = NoodlesRawRecord::default();
        while source.read_next(&mut buf).expect("read") {}
        assert!(
            source.done,
            "the early stop must latch, so a drained query does not resume"
        );
        assert!(
            source.chunks.len() > 0,
            "the scan consumed every chunk for a 100 bp region of a 100 kb \
             contig — it read on instead of stopping early"
        );
    }

    /// The source owns its reader so the pool can have it back; a stream that
    /// ended must be able to hand over a reader that still works, or every
    /// query after the first would open a fresh one.
    #[test]
    fn the_reader_survives_the_query_and_comes_back_usable() {
        let (_dir, path, header) = fixture();
        let index = load_alignment_index(&path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(&path)
            .expect("open bam");
        reader.read_header().expect("read header");

        let mut source = BamRegionSource::new(
            reader,
            &header,
            BamRegionSource::plan(&header, &index, region(0, 1, 100), &path).expect("plan"),
            0,
        );
        let mut buf = NoodlesRawRecord::default();
        while source.read_next(&mut buf).expect("read") {}

        // Hand it back, then drive a second query with the very same reader.
        let reader = source.into_reader();
        let mut second = BamRegionSource::new(
            reader,
            &header,
            BamRegionSource::plan(&header, &index, region(1, 1, 200), &path).expect("plan"),
            0,
        );
        let mut seen = 0;
        while second.read_next(&mut buf).expect("read") {
            seen += 1;
        }
        assert!(
            seen > 0,
            "the returned reader must still be able to seek and read"
        );
    }

    #[test]
    fn a_region_naming_a_contig_the_file_does_not_have_is_an_error() {
        let (_dir, path, header) = fixture();
        let index = load_alignment_index(&path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(&path)
            .expect("open bam");
        reader.read_header().expect("read header");

        assert!(matches!(
            BamRegionSource::plan(&header, &index, region(9, 1, 100), &path),
            Err(AlignmentFileError::Region { .. })
        ));
    }

    // -----------------------------------------------------------------
    // The order guard (C3) — T4a..T4d
    // -----------------------------------------------------------------

    /// A `MappedRead` at a genome position. Only `ref_id` and `pos` matter to
    /// the guard, so the rest is minimal — the guard is a pure adapter over the
    /// filtered stream and never looks at sequence or CIGAR.
    fn read_at(qname: &str, ref_id: usize, pos: u64) -> MappedRead {
        MappedRead {
            qname: qname.as_bytes().to_vec(),
            flag: 0,
            ref_id,
            pos,
            mapq: 60,
            cigar: Vec::new(),
            seq: Vec::new(),
            qual: Vec::new(),
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            source_file_index: 0,
        }
    }

    fn planted_path() -> Arc<Path> {
        Arc::from(Path::new("/data/planted.bam"))
    }

    /// An iterator that yields `None` and then **keeps going** — what a
    /// non-fused inner looks like. `std::iter::from_fn` cannot model this,
    /// because it fuses itself.
    struct Resuming {
        items: Vec<Option<Result<MappedRead, ReadFilterError>>>,
        next_index: usize,
    }

    impl Iterator for Resuming {
        type Item = Result<MappedRead, ReadFilterError>;

        fn next(&mut self) -> Option<Self::Item> {
            let item = self.items.get_mut(self.next_index)?.take();
            self.next_index += 1;
            item
        }
    }

    /// Drive the guard over a planted stream and collect what a caller would
    /// see: the read names that surfaced, and the error if one did.
    fn through_the_guard(
        reads: Vec<MappedRead>,
    ) -> (Vec<String>, Option<AlignmentFileError>, usize) {
        let stream = reads.into_iter().map(Ok);
        let mut guard = OrderVerified::new(stream, planted_path());

        let mut names = Vec::new();
        let mut error = None;
        for item in &mut guard {
            match item {
                Ok(read) => names.push(String::from_utf8_lossy(&read.qname).into_owned()),
                Err(e) => {
                    error = Some(e);
                    break;
                }
            }
        }
        // How many items the iterator yields *after* the error — must be zero,
        // because the guard fuses.
        let after = guard.count();
        (names, error, after)
    }

    /// **T4a — a planted position regression is fatal.**
    ///
    /// The header can say `SO:coordinate` and the data can still be unsorted;
    /// that is the case the open gate structurally cannot see. The error names
    /// both keys, so the message can say *where* the file breaks rather than
    /// only that it does.
    ///
    /// Mutation-verified: deleting the comparison in `OrderVerified::next`
    /// makes this test fail, which is the whole point — a guard that never
    /// fires looks exactly like a guard that works.
    #[test]
    fn t4a_a_position_regression_within_a_contig_is_an_error() {
        let (names, error, after) = through_the_guard(vec![
            read_at("a", 0, 100),
            read_at("b", 0, 200),
            read_at("c", 0, 150), // backwards
            read_at("d", 0, 300),
        ]);

        assert_eq!(names, vec!["a", "b"], "the reads before the break flow");
        match error.expect("the regression must be fatal") {
            AlignmentFileError::OutOfOrderRead {
                path,
                previous,
                current,
            } => {
                assert_eq!(previous.position.get(), 200);
                assert_eq!(current.position.get(), 150);
                // Asserted because `path` is a constructor argument: a wiring
                // mistake at C4 would produce a correct-looking error naming
                // the wrong file.
                assert_eq!(path, PathBuf::from("/data/planted.bam"));
            }
            other => panic!("expected OutOfOrderRead, got {other:?}"),
        }
        assert_eq!(after, 0, "and the iterator fuses — no read after the error");
    }

    /// **T4b — a contig-order regression is the same violation.** The key is
    /// `(contig, position)`, so a read on an earlier contig after a later one
    /// is caught by the same comparison, even though its *position* went up.
    #[test]
    fn t4b_a_contig_order_regression_is_the_same_error() {
        let (names, error, _) = through_the_guard(vec![
            read_at("a", 1, 100),
            read_at("b", 0, 900), // earlier contig, later position
        ]);

        assert_eq!(names, vec!["a"]);
        match error.expect("a contig regression must be fatal") {
            AlignmentFileError::OutOfOrderRead {
                previous, current, ..
            } => {
                assert_eq!(previous.contig.get(), 1);
                assert_eq!(current.contig.get(), 0);
            }
            other => panic!("expected OutOfOrderRead, got {other:?}"),
        }
    }

    /// **T4c — equal keys are legal.** Several reads may start at the same
    /// base; only a strict decrease is a violation. A guard using `<=` would
    /// reject every pile-up in every real file.
    #[test]
    fn t4c_reads_sharing_a_start_position_are_not_a_regression() {
        let (names, error, _) = through_the_guard(vec![
            read_at("a", 0, 100),
            read_at("b", 0, 100),
            read_at("c", 0, 100),
            read_at("d", 0, 101),
        ]);

        assert_eq!(names, vec!["a", "b", "c", "d"]);
        assert!(error.is_none(), "equal positions are ordinary, not a fault");
    }

    /// **T4d — the check does not span queries.** A caller may query region B
    /// and then region A; the second is a new forward scan, not a regression.
    /// This is why the state lives in the iterator and never on the handle —
    /// carrying it there would turn legitimate random access into a spurious
    /// error.
    #[test]
    fn t4d_querying_a_later_region_then_an_earlier_one_is_not_a_regression() {
        let (later, error, _) = through_the_guard(vec![read_at("b1", 0, 9000)]);
        assert_eq!(later, vec!["b1"]);
        assert!(error.is_none());

        // A *separate* guard, as a second query gets — starting at position 10,
        // far behind where the first one ended.
        let (earlier, error, _) = through_the_guard(vec![read_at("a1", 0, 10)]);
        assert_eq!(earlier, vec!["a1"]);
        assert!(
            error.is_none(),
            "a new region query starts a new scan; it cannot regress against \
             a previous one"
        );

        // And the fresh guard is *armed*, not merely permissive — a guard that
        // had been silently disabled would also pass the assertion above.
        let (_, error, _) = through_the_guard(vec![read_at("a1", 0, 10), read_at("a2", 0, 5)]);
        assert!(
            matches!(error, Some(AlignmentFileError::OutOfOrderRead { .. })),
            "the second query's own guard still catches its own regression"
        );
    }

    /// **The other half of the fuse.** `FusedIterator` is claimed without
    /// requiring `I: FusedIterator`, so the latch has to hold against an inner
    /// that returns `None` and then resumes — otherwise a read could reappear
    /// after the stream was reported ended. Only the *error* path's fuse was
    /// covered before; deleting `done = true` from the clean-exhaustion arm
    /// left every test green.
    #[test]
    fn the_guard_stays_ended_even_if_what_it_wraps_resumes() {
        let inner = Resuming {
            items: vec![
                Some(Ok(read_at("a", 0, 1))),
                None,
                Some(Ok(read_at("b", 0, 2))),
            ],
            next_index: 0,
        };
        let mut guard = OrderVerified::new(inner, planted_path());

        assert!(matches!(guard.next(), Some(Ok(_))));
        assert!(guard.next().is_none(), "the inner reported end of input");
        assert!(
            guard.next().is_none(),
            "and the guard stays ended — the latch holds, not just the inner"
        );
    }

    /// An error from the filter below is wrapped and fuses the stream, rather
    /// than being mistaken for a clean end of input.
    #[test]
    fn a_filter_error_is_wrapped_and_fuses_the_guard() {
        let stream = vec![
            Ok(read_at("a", 0, 1)),
            Err(ReadFilterError::Source(io::Error::other("planted"))),
            Ok(read_at("b", 0, 2)),
        ]
        .into_iter();
        let mut guard = OrderVerified::new(stream, planted_path());

        assert!(matches!(guard.next(), Some(Ok(_))));
        assert!(matches!(
            guard.next(),
            Some(Err(AlignmentFileError::Filter(_)))
        ));
        assert!(
            guard.next().is_none(),
            "fused: the read after the error is not reachable"
        );
    }

    // -----------------------------------------------------------------
    // The CRAM `.crai` walk (C5)
    // -----------------------------------------------------------------

    /// A `.crai` entry for `contig` covering `[start, start + span - 1]`.
    fn crai_index(records: Vec<cram::crai::Record>) -> cram::crai::Index {
        records
    }

    /// A `.crai` entry for `contig` covering `[start, start + span - 1]`, at
    /// `offset` — which the tests use to tell entries apart.
    fn crai_entry(
        contig: Option<usize>,
        start: usize,
        span: usize,
        offset: u64,
    ) -> cram::crai::Record {
        cram::crai::Record::new(
            contig,
            noodles_core::Position::new(start),
            span,
            offset,
            0,
            0,
        )
    }

    fn three_contig_header() -> sam::Header {
        bam_header(&[
            ("chr1", 100, None),
            ("chr2", 200, None),
            ("chr3", 300, None),
        ])
    }

    /// **The fix for production's O(n) prefix rescan** — and, more importantly,
    /// for an ordering assumption that does not hold. A query on a late contig
    /// must reach that contig's entries without walking (or mis-searching) the
    /// earlier ones.
    ///
    /// Hand-built rather than read from a file, because the file fixture cannot
    /// hold more than one contig (`test_fixtures::indexed_cram`).
    #[test]
    fn the_crai_is_grouped_so_each_contig_sees_only_its_own_entries() {
        let index = crai_index(vec![
            crai_entry(Some(0), 1, 50, 100),
            crai_entry(Some(0), 51, 50, 200),
            crai_entry(Some(1), 1, 100, 300),
            crai_entry(Some(1), 101, 100, 400),
            crai_entry(Some(2), 1, 300, 500),
        ]);
        let grouped = group_crai_by_contig(&index, 3);
        let header = three_contig_header();

        for (contig, expected_offsets) in [
            (0u32, vec![100u64, 200]),
            (1, vec![300, 400]),
            (2, vec![500]),
        ] {
            let plan =
                CramRegionSource::plan(&header, &grouped, region(contig, 1, 10)).expect("plan");
            let offsets: Vec<u64> = plan.entries.iter().map(|e| e.offset()).collect();
            assert_eq!(offsets, expected_offsets, "contig {contig}");
        }
    }

    /// **The ordering assumption a binary search would have made, violated by
    /// noodles itself.** Within a multi-reference slice, `fs::index` sorts by
    /// `Option<usize>` — and `None < Some(0)` — so the *unplaced* entry is
    /// emitted first. A `partition_point` over that is unspecified, and the
    /// walk would then meet a foreign entry and report end-of-input, losing
    /// every read of the region with no error at all.
    ///
    /// Grouping assumes nothing about order, so an interleaved index is fine.
    #[test]
    fn an_unplaced_entry_before_the_placed_ones_does_not_hide_a_contig() {
        let index = crai_index(vec![
            crai_entry(Some(0), 1, 50, 100),
            crai_entry(None, 1, 0, 150),
            crai_entry(Some(1), 1, 100, 300),
        ]);
        let grouped = group_crai_by_contig(&index, 3);

        let plan = CramRegionSource::plan(&three_contig_header(), &grouped, region(1, 1, 10))
            .expect("plan");
        let offsets: Vec<u64> = plan.entries.iter().map(|e| e.offset()).collect();
        assert_eq!(
            offsets,
            vec![300],
            "contig 1's entry must be found even though an unplaced entry \
             precedes it — the case that silently returned nothing before"
        );

        // And the unplaced entry belongs to no contig at all.
        for (contig, entries) in grouped.iter().enumerate() {
            assert!(
                entries.iter().all(|e| e.reference_sequence_id().is_some()),
                "an unplaced entry leaked into contig {contig}"
            );
        }
    }

    /// A contig with no entries yields an empty walk rather than someone
    /// else's entries.
    #[test]
    fn a_contig_absent_from_the_crai_has_an_empty_walk() {
        let index = crai_index(vec![
            crai_entry(Some(0), 1, 50, 100),
            crai_entry(Some(2), 1, 300, 500),
        ]);
        let grouped = group_crai_by_contig(&index, 3);

        let plan = CramRegionSource::plan(&three_contig_header(), &grouped, region(1, 1, 10))
            .expect("plan");
        assert!(plan.entries.is_empty());
    }

    /// **The CRAM container-level early stop, with a discriminating
    /// assertion.** Draining a region near the start of a long contig must
    /// leave most of the `.crai` unwalked; a test that only checked the reads
    /// were right would pass with the stop deleted, because deleting it changes
    /// the *work*, not the answer — a whole-file decode per region, ~10⁶ times.
    ///
    /// Driven through `CramRegionSource` directly, since the cursor is not
    /// visible once the source is wrapped in the filter and the order guard.
    #[test]
    fn the_crai_walk_stops_once_it_passes_the_region() {
        use crate::ng::read::input::test_fixtures::multi_container_cram;

        const CONTIG_LENGTH: usize = 400_000;
        let (_cram_dir, cram_path, _fasta_dir, fasta) = multi_container_cram(CONTIG_LENGTH, 30_000);

        let index = load_alignment_index(&cram_path).expect("load crai");
        let repository =
            crate::bam::alignment_input::build_fasta_repository(&fasta).expect("repository");
        let header_for_plan = bam_header(&[("chr1", CONTIG_LENGTH, None)]);
        let grouped = group_crai_by_contig(
            match &index {
                AlignmentIndex::Crai(crai) => crai,
                _ => panic!("a .cram must carry a .crai"),
            },
            1,
        );
        assert!(
            grouped[0].len() > 1,
            "the fixture must span several containers, or there is nothing to \
             stop short of — got {}",
            grouped[0].len()
        );

        let mut reader = cram::io::Reader::new(File::open(&cram_path).expect("open"));
        let parsed_header = reader.read_header().expect("read header");

        let plan =
            CramRegionSource::plan(&header_for_plan, &grouped, region(0, 1, 200)).expect("plan");
        let total_entries = plan.entries.len();
        let mut source = CramRegionSource::new(reader, &parsed_header, repository, plan, 0);

        let mut buf = NoodlesRawRecord::default();
        let mut seen = 0;
        while source.read_next(&mut buf).expect("read") {
            seen += 1;
        }

        assert!(seen > 0, "the region is covered");
        assert!(
            source.next_index_record < total_entries,
            "the walk consumed all {total_entries} entries for a 200 bp region \
             of a 400 kb contig — it decoded to the end instead of stopping"
        );
    }

    /// The span-based skip: a container that ends *before* the region must be
    /// stepped over without being decoded. Observable as the cursor advancing
    /// past entries that produced no reads.
    #[test]
    fn the_crai_walk_skips_containers_that_end_before_the_region() {
        use crate::ng::read::input::test_fixtures::multi_container_cram;

        const CONTIG_LENGTH: usize = 400_000;
        let (_cram_dir, cram_path, _fasta_dir, fasta) = multi_container_cram(CONTIG_LENGTH, 30_000);

        let index = load_alignment_index(&cram_path).expect("load crai");
        let repository =
            crate::bam::alignment_input::build_fasta_repository(&fasta).expect("repository");
        let header_for_plan = bam_header(&[("chr1", CONTIG_LENGTH, None)]);
        let grouped = group_crai_by_contig(
            match &index {
                AlignmentIndex::Crai(crai) => crai,
                _ => panic!("a .cram must carry a .crai"),
            },
            1,
        );

        let mut reader = cram::io::Reader::new(File::open(&cram_path).expect("open"));
        let parsed_header = reader.read_header().expect("read header");

        // A region in the *last* stretch of the contig: every earlier container
        // ends before it and must be skipped rather than decoded.
        let late = (CONTIG_LENGTH as u64 * 3) / 4;
        let plan = CramRegionSource::plan(&header_for_plan, &grouped, region(0, late, late + 200))
            .expect("plan");
        let mut source = CramRegionSource::new(reader, &parsed_header, repository, plan, 0);

        let mut buf = NoodlesRawRecord::default();
        let mut seen = 0;
        while source.read_next(&mut buf).expect("read") {
            seen += 1;
        }

        assert!(seen > 0, "the late region is covered");
        assert!(
            source.last_decoded_offset.is_some(),
            "at least one container was decoded"
        );
    }

    #[test]
    fn a_cram_region_naming_a_contig_the_header_lacks_is_an_error() {
        let grouped = group_crai_by_contig(&crai_index(vec![crai_entry(Some(0), 1, 50, 0)]), 3);
        assert!(matches!(
            CramRegionSource::plan(&three_contig_header(), &grouped, region(9, 1, 10)),
            Err(AlignmentFileError::Region { .. })
        ));
    }

    /// The CRAM planner rejects an inverted or empty region, as the BAM one
    /// does — the two must not disagree about what a valid region is.
    #[test]
    fn a_cram_planner_rejects_an_inverted_region() {
        let grouped = group_crai_by_contig(&crai_index(vec![crai_entry(Some(0), 1, 50, 0)]), 3);
        assert!(matches!(
            CramRegionSource::plan(&three_contig_header(), &grouped, region(0, 500, 100)),
            Err(AlignmentFileError::Region { .. })
        ));
    }

    #[test]
    fn an_inverted_region_is_an_error() {
        let (_dir, path, header) = fixture();
        let index = load_alignment_index(&path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(&path)
            .expect("open bam");
        reader.read_header().expect("read header");

        assert!(matches!(
            BamRegionSource::plan(&header, &index, region(0, 500, 100), &path),
            Err(AlignmentFileError::Region { .. })
        ));
    }
}
