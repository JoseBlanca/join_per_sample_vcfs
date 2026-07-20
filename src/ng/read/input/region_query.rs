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

// This module's items are built by `reads_in_region`, which lands in step C4;
// until then only its own tests construct them. `expect` rather than `allow`,
// so once C4 makes everything live the expectation goes unfulfilled and the
// build fails naming this line.
//
// Module-level rather than per-item because rustc reports the associated items
// as one group, so per-item expectations go unfulfilled individually. The cost
// is that a genuinely-dead *new* item would hide here for the one commit until
// C4; the enforcement that matters — this attribute cannot outlive C4 — holds.
#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "reads_in_region (C4) builds and drives this module"
    )
)]

use std::fs::File;
use std::io;
use std::iter::FusedIterator;
use std::path::Path;
use std::sync::Arc;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_csi::BinningIndex;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_sam as sam;

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

impl<'a> BamRegionSource<'a> {
    /// Query the index for the region's chunks and prepare to scan them.
    ///
    /// The index query happens here, once per region — an in-memory lookup on
    /// the index parsed at open, never a re-parse (spec §3.3).
    pub(crate) fn new(
        reader: bam::io::Reader<bgzf::io::Reader<File>>,
        header: &'a sam::Header,
        index: &AlignmentIndex,
        region: GenomeRegion,
        path: &Path,
        source_file_index: usize,
    ) -> Result<Self, AlignmentFileError> {
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

        Ok(Self {
            reader,
            header,
            chunks: chunks.into_iter(),
            current_chunk_end: None,
            target_reference_sequence_id,
            region,
            source_file_index,
            done: false,
        })
    }

    /// Give the reader back, for the pool.
    pub(crate) fn into_reader(self) -> bam::io::Reader<bgzf::io::Reader<File>> {
        self.reader
    }

    /// Whether this record's reference footprint touches the region.
    ///
    /// A record with no mapped position never overlaps. The comparison is on
    /// `[alignment_start, alignment_end]` — the *footprint*, not the start —
    /// because a read beginning before the region can still cover it.
    fn overlaps_region(&self, record: &sam::alignment::RecordBuf) -> bool {
        match (record.alignment_start(), record.alignment_end()) {
            (Some(first), Some(last)) => {
                usize::from(first) as u64 <= self.region.end.get()
                    && usize::from(last) as u64 >= self.region.start.get()
            }
            _ => false,
        }
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
            if !on_target_contig || !self.overlaps_region(&buf.record) {
                continue;
            }

            return Ok(true);
        }
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

        let mut source = BamRegionSource::new(reader, header, &index, region, path, 0)
            .expect("build the source");

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

        let mut source =
            BamRegionSource::new(reader, &header, &index, region(0, 1, 100_000), &path, 0)
                .expect("source");

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

        let mut source = BamRegionSource::new(reader, &header, &index, region(0, 1, 100), &path, 0)
            .expect("source");

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

        let mut source = BamRegionSource::new(reader, &header, &index, region(0, 1, 100), &path, 0)
            .expect("source");
        let mut buf = NoodlesRawRecord::default();
        while source.read_next(&mut buf).expect("read") {}

        // Hand it back, then drive a second query with the very same reader.
        let reader = source.into_reader();
        let mut second = BamRegionSource::new(reader, &header, &index, region(1, 1, 200), &path, 0)
            .expect("source");
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
            BamRegionSource::new(reader, &header, &index, region(9, 1, 100), &path, 0),
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

    #[test]
    fn an_inverted_region_is_an_error() {
        let (_dir, path, header) = fixture();
        let index = load_alignment_index(&path).expect("load index");
        let mut reader = bam::io::reader::Builder
            .build_from_path(&path)
            .expect("open bam");
        reader.read_header().expect("read header");

        assert!(matches!(
            BamRegionSource::new(reader, &header, &index, region(0, 500, 100), &path, 0),
            Err(AlignmentFileError::Region { .. })
        ));
    }
}
