//! BAM record-stream decoders. Sibling of [`crate::bam::cram_input`]
//! — see that module's docs for the design; this file is the BAM
//! analogue with the same factory shape so the format-agnostic
//! merge in [`crate::bam::alignment_input`] can consume either.
//!
//! - [`OwnedBamRecords`] — linear-order owned iterator over a single
//!   BAM's records. Unlike CRAM, BAM's `read_record_buf` is
//!   cleanly idempotent at EOF (returns 0 repeatedly), so no
//!   EOF-latch field is needed; the merge's `BufferedPeekable`
//!   post-EOF poll is harmless here.
//! - [`OwnedIndexedBamRecords`] — `.csi`- or `.bai`-driven owned
//!   iterator over a single BAM's records on one target contig.
//!   Walks chunks from `BinningIndex::query` and filters records
//!   on `reference_sequence_id` (a BAM file's chunks for one
//!   contig may technically contain records on other contigs at
//!   chunk boundaries).
//! - [`open_bam_record_stream`] / [`open_indexed_bam_record_stream`]
//!   — the open helpers `AlignmentMergedReader::{new, query}`
//!   call once per input BAM. Mirror of `cram_input`'s two
//!   helpers, returning the same `(sam::Header, AlignmentRecordsIter)`
//!   / `AlignmentRecordsIter` shapes.

use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Arc;
use std::vec;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::binning_index::BinningIndex;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_sam as sam;

use super::alignment_input::AlignmentRecordsIter;
use super::errors::AlignmentInputError;

// ---------------------------------------------------------------------
// BAM-side index payload
// ---------------------------------------------------------------------

/// BAM index payload handed to the indexed open helper. Either index
/// format works at noodles' [`BinningIndex`] layer; the caller (the
/// `index_preflight` layer, downstream) hands us whichever it
/// loaded.
#[derive(Clone)]
pub(super) enum BamIndex {
    /// `.bai` — the legacy BAM index. 16 kbp bins; cannot address
    /// contigs above 512 Mbp.
    Bai(Arc<bam::bai::Index>),
    /// `.csi` — coordinate-sorted index. Preferred on read and
    /// always produced on build, because it has no contig-length
    /// limit.
    Csi(Arc<noodles_csi::Index>),
}

// ---------------------------------------------------------------------
// Owned BAM record iterator (linear order)
// ---------------------------------------------------------------------

/// Owns a BAM `Reader<bgzf::Reader<File>>` plus its `sam::Header`
/// and yields decoded `RecordBuf`s. Equivalent of
/// [`crate::bam::cram_input::OwnedCramRecords`] for the BAM
/// decoder; the simpler shape (no compression-header / slice /
/// repository plumbing) is because BAM stores records flat and
/// inline-sequence, with no per-container batching.
struct OwnedBamRecords {
    reader: bam::io::Reader<bgzf::io::Reader<File>>,
    header: sam::Header,
}

impl Iterator for OwnedBamRecords {
    type Item = io::Result<sam::alignment::RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record_buf = sam::alignment::RecordBuf::default();
        match self.reader.read_record_buf(&self.header, &mut record_buf) {
            // Clean EOF. `read_record_buf` is idempotent at EOF (it
            // calls `read_block_size` -> `read_exact_or_eof`, which
            // gracefully signals 0 bytes both at the first EOF and
            // on every subsequent call), so no separate latch field
            // is needed — unlike the CRAM analogue.
            Ok(0) => None,
            Ok(_) => Some(Ok(record_buf)),
            Err(e) => Some(Err(e)),
        }
    }
}

// ---------------------------------------------------------------------
// Indexed (per-contig) owned BAM record iterator
// ---------------------------------------------------------------------

/// Owned iterator that walks a BAM via the chunks reported by its
/// `.csi` or `.bai` index, yielding only records that align to
/// `target_reference_sequence_id`. Mirror of
/// [`crate::bam::cram_input::OwnedIndexedCramRecords`] driven by
/// `BinningIndex::query` instead of by the CRAM index slice; the
/// chunk-walking loop reproduces what noodles' own
/// [`bam::io::reader::Query`] does internally, but with the reader
/// owned (not borrowed) so the iterator is `'static + Send`.
struct OwnedIndexedBamRecords {
    reader: bam::io::Reader<bgzf::io::Reader<File>>,
    /// Shared parsed header. We use the caller's already-loaded
    /// `Arc<sam::Header>` rather than re-parsing per worker — the
    /// driver has cross-validated every input's header is
    /// identical.
    header: Arc<sam::Header>,
    chunks: vec::IntoIter<Chunk>,
    /// Virtual-position where the current chunk ends. `None` means
    /// "advance to the next chunk on the next poll". Set to `Some`
    /// after a successful seek, cleared back to `None` once the
    /// chunk is exhausted (either because we read past its end or
    /// because we hit EOF mid-chunk on the last chunk).
    current_chunk_end: Option<bgzf::VirtualPosition>,
    target_reference_sequence_id: usize,
}

impl Iterator for OwnedIndexedBamRecords {
    type Item = io::Result<sam::alignment::RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // 1. Land inside a chunk (advance + seek if we weren't
            //    inside one already).
            let chunk_end = match self.current_chunk_end {
                Some(e) => e,
                None => {
                    let chunk = self.chunks.next()?;
                    if let Err(e) = self.reader.get_mut().seek(chunk.start()) {
                        return Some(Err(e));
                    }
                    self.current_chunk_end = Some(chunk.end());
                    chunk.end()
                }
            };

            // 2. If the previous record read pushed us past the
            //    current chunk's end, fall through to the next
            //    chunk on the next loop iteration. Mirrors the
            //    pre-read bound check in noodles'
            //    `csi::io::Query<R>::read` impl.
            if self.reader.get_ref().virtual_position() >= chunk_end {
                self.current_chunk_end = None;
                continue;
            }

            // 3. Read one record at the current virtual position.
            let mut record_buf = sam::alignment::RecordBuf::default();
            match self.reader.read_record_buf(&self.header, &mut record_buf) {
                Ok(0) => {
                    // EOF mid-chunk (the writer truncated, or the
                    // chunk end is exactly at file EOF). Move on
                    // to the next chunk if there is one.
                    self.current_chunk_end = None;
                    continue;
                }
                Ok(_) => {
                    if record_buf.reference_sequence_id() == Some(self.target_reference_sequence_id)
                    {
                        return Some(Ok(record_buf));
                    }
                    // Record on a different contig — happens at
                    // chunk boundaries when a chunk straddles
                    // contigs. Skip; the loop re-checks the
                    // virtual position and reads on.
                }
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

// ---------------------------------------------------------------------
// Open helpers
// ---------------------------------------------------------------------

/// Open `path` as a BAM, read the SAM header, and return the
/// header + an owned record-stream iterator ready to feed into
/// [`super::alignment_input::OpenAlignmentFile`]. BAM has no
/// separate file-definition / version step (the magic-byte check
/// is folded into `read_header`).
pub(super) fn open_bam_record_stream(
    path: &Path,
) -> Result<(sam::Header, AlignmentRecordsIter), AlignmentInputError> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .map_err(|source| AlignmentInputError::OpenFailed {
            path: path.to_path_buf(),
            source,
        })?;
    let header = reader
        .read_header()
        .map_err(|source| AlignmentInputError::OpenFailed {
            path: path.to_path_buf(),
            source,
        })?;
    let records: AlignmentRecordsIter = Box::new(OwnedBamRecords {
        reader,
        header: header.clone(),
    });
    Ok((header, records))
}

/// Open `path` as a BAM, advance past its header (the caller
/// already holds the validated header via `header`), ask the index
/// for the chunks covering the entire target contig, and return an
/// owned per-contig indexed record-stream iterator.
///
/// `target_reference_sequence_id` must be valid against the
/// canonical contig list the caller resolved against (the driver
/// does this once at startup).
pub(super) fn open_indexed_bam_record_stream(
    path: &Path,
    header: Arc<sam::Header>,
    index: BamIndex,
    target_reference_sequence_id: usize,
) -> Result<AlignmentRecordsIter, AlignmentInputError> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .map_err(|source| AlignmentInputError::OpenFailed {
            path: path.to_path_buf(),
            source,
        })?;
    // Advance past the header so any chunk-seek lands at a record
    // boundary. The parsed header is discarded — `header` already
    // carries the validated one.
    reader
        .read_header()
        .map_err(|source| AlignmentInputError::OpenFailed {
            path: path.to_path_buf(),
            source,
        })?;

    // Ask the index for the chunks covering the whole target contig
    // (interval = unbounded). `BinningIndex::query` is the polymorphic
    // entry point both `.bai` and `.csi` implement.
    let chunks = match &index {
        BamIndex::Bai(idx) => idx.query(target_reference_sequence_id, Interval::from(..)),
        BamIndex::Csi(idx) => idx.query(target_reference_sequence_id, Interval::from(..)),
    }
    .map_err(|source| AlignmentInputError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    let records: AlignmentRecordsIter = Box::new(OwnedIndexedBamRecords {
        reader,
        header,
        chunks: chunks.into_iter(),
        current_chunk_end: None,
        target_reference_sequence_id,
    });
    Ok(records)
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use noodles_core::Position;
    use noodles_csi::binning_index::Indexer;
    use noodles_sam::alignment::Record as _;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::ReferenceSequence;
    use tempfile::TempDir;

    use super::*;

    fn header_two_contigs() -> sam::Header {
        sam::Header::builder()
            .set_header(Default::default())
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(64).unwrap() }),
            )
            .add_reference_sequence(
                "sq1",
                Map::<ReferenceSequence>::new(const { NonZero::new(64).unwrap() }),
            )
            .build()
    }

    fn record_on(reference_sequence_id: usize, start: usize) -> RecordBuf {
        RecordBuf::builder()
            .set_reference_sequence_id(reference_sequence_id)
            .set_flags(Flags::default())
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
            .build()
    }

    /// Write a synthetic coordinate-sorted BAM with the given
    /// records to a temp dir; return the dir handle (so the path
    /// stays alive) and the file path.
    fn write_bam(records: &[RecordBuf], header: &sam::Header) -> (TempDir, std::path::PathBuf) {
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

    /// Build a `.bai` for the just-written BAM by scanning it. Mirror
    /// of the helper noodles-bam's own query tests use; we keep our
    /// own copy so the bam_input tests don't pull in noodles' test
    /// support feature.
    fn build_bai_in_memory(bam_path: &std::path::Path) -> bam::bai::Index {
        let mut reader = bam::io::reader::Builder
            .build_from_path(bam_path)
            .expect("open bam to index");
        let header = reader.read_header().expect("read header for index");

        let mut indexer: Indexer<
            noodles_csi::binning_index::index::reference_sequence::index::LinearIndex,
        > = Indexer::default();
        let mut chunk_start = reader.get_ref().virtual_position();
        let mut record = bam::Record::default();

        while reader.read_record(&mut record).expect("read record") != 0 {
            let chunk_end = reader.get_ref().virtual_position();
            let alignment_context = match (
                record.reference_sequence_id().transpose().expect("ref id"),
                record.alignment_start().transpose().expect("start"),
                record.alignment_end().transpose().expect("end"),
            ) {
                (Some(id), Some(start), Some(end)) => {
                    let is_mapped = !record.flags().is_unmapped();
                    Some((id, start, end, is_mapped))
                }
                _ => None,
            };
            let chunk = Chunk::new(chunk_start, chunk_end);
            indexer
                .add_record(alignment_context, chunk)
                .expect("add record");
            chunk_start = chunk_end;
        }

        indexer.build(header.reference_sequences().len())
    }

    #[test]
    fn open_bam_record_stream_streams_three_records_in_order() {
        let header = header_two_contigs();
        let records = [record_on(0, 5), record_on(0, 12), record_on(1, 1)];
        let (_dir, bam_path) = write_bam(&records, &header);

        let (returned_header, mut stream) =
            open_bam_record_stream(&bam_path).expect("open bam stream");

        assert_eq!(returned_header.reference_sequences().len(), 2);

        let actual: Vec<RecordBuf> = std::iter::from_fn(|| stream.next())
            .map(|r| r.expect("decode record"))
            .collect();
        assert_eq!(actual.len(), 3, "stream yielded {} records", actual.len());
        assert_eq!(actual[0], records[0]);
        assert_eq!(actual[1], records[1]);
        assert_eq!(actual[2], records[2]);
    }

    #[test]
    fn open_bam_record_stream_returns_none_at_eof_idempotently() {
        // Empty body (header-only BAM): the iterator must return
        // None on every poll past the first, mirroring the no-latch
        // contract documented on OwnedBamRecords.
        let header = header_two_contigs();
        let (_dir, bam_path) = write_bam(&[], &header);

        let (_header, mut stream) = open_bam_record_stream(&bam_path).expect("open empty bam");

        assert!(stream.next().is_none(), "empty BAM should yield no records");
        // Poll a few more times — must stay None, must not surface
        // an error (this is exactly the bug the CRAM EOF-latch
        // guards against in OwnedCramRecords).
        for _ in 0..3 {
            assert!(
                stream.next().is_none(),
                "post-EOF poll must not surface a spurious record or error"
            );
        }
    }

    #[test]
    fn open_indexed_bam_record_stream_yields_only_target_contig() {
        let header = header_two_contigs();
        // Three records: sq0 (irrelevant), sq1 at 1, sq1 at 8.
        // Target contig is sq1 (id=1).
        let records = [record_on(0, 1), record_on(1, 1), record_on(1, 8)];
        let (_dir, bam_path) = write_bam(&records, &header);

        let bai_index = build_bai_in_memory(&bam_path);
        let header_arc = Arc::new(header.clone());

        let mut stream = open_indexed_bam_record_stream(
            &bam_path,
            Arc::clone(&header_arc),
            BamIndex::Bai(Arc::new(bai_index)),
            /* target_reference_sequence_id = */ 1,
        )
        .expect("open indexed bam stream");

        let actual: Vec<RecordBuf> = std::iter::from_fn(|| stream.next())
            .map(|r| r.expect("decode record"))
            .collect();
        assert_eq!(actual.len(), 2, "expected the two sq1 records");
        assert!(
            actual.iter().all(|r| r.reference_sequence_id() == Some(1)),
            "every yielded record must be on the queried contig"
        );
    }

    #[test]
    fn open_indexed_bam_record_stream_yields_nothing_for_empty_contig() {
        let header = header_two_contigs();
        // All records on sq0; query sq1 -> no records.
        let records = [record_on(0, 1), record_on(0, 8)];
        let (_dir, bam_path) = write_bam(&records, &header);

        let bai_index = build_bai_in_memory(&bam_path);
        let header_arc = Arc::new(header.clone());

        let mut stream = open_indexed_bam_record_stream(
            &bam_path,
            Arc::clone(&header_arc),
            BamIndex::Bai(Arc::new(bai_index)),
            /* target_reference_sequence_id = */ 1,
        )
        .expect("open indexed bam stream");

        assert!(
            stream.next().is_none(),
            "no records on the queried contig should mean no yields"
        );
        // Idempotent at exhaustion.
        for _ in 0..3 {
            assert!(stream.next().is_none());
        }
    }
}
