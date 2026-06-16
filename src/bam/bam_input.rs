//! BAM record-stream decoders. Sibling of [`crate::bam::cram_input`]
//! — see that module's docs for the design; this file is the BAM
//! analogue with the same factory shape so the format-agnostic
//! merge in [`crate::bam::alignment_input`] can consume either.
//!
//! - [`OwnedBamRecords`] — linear-order owned iterator over a single
//!   BAM's records. Unlike CRAM, BAM's `read_record_buf` is
//!   cleanly idempotent at EOF (returns 0 repeatedly), so no
//!   EOF-latch field is needed.
//! - [`open_bam_record_stream`] — the open helper `AlignmentMergedReader::new`
//!   calls once per input BAM, returning `(sam::Header, AlignmentRecordsIter)`.
//! - [`open_bam_reader_with_header`] — opens + reads the header, used by
//!   `load_pileup_inputs` and the pooled [`crate::bam::segment_reader`].
//! - [`query_interval`] / [`BamIndex`] — the index-query shapes the pooled
//!   segment reader reuses.

use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Arc;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_core::region::Interval;
use noodles_sam as sam;

use super::alignment_input::{AlignmentRecordsIter, ContigInterval};
use super::errors::AlignmentInputError;

// ---------------------------------------------------------------------
// BAM-side index payload
// ---------------------------------------------------------------------

/// BAM index payload handed to the indexed open helper. Either index
/// format works at noodles' [`BinningIndex`] layer; the caller (the
/// `index_preflight` layer, downstream) hands us whichever it
/// loaded.
#[derive(Clone)]
#[non_exhaustive]
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
    let (reader, header) = open_bam_reader_with_header(path)?;
    let records: AlignmentRecordsIter = Box::new(OwnedBamRecords {
        reader,
        header: header.clone(),
    });
    Ok((header, records))
}

/// Shared BAM open + header-read sequence used by
/// `open_bam_record_stream` so the magic-byte check and the error
/// wrapping stay in one place. Also used by the pooled
/// [`crate::bam::segment_reader`] to open a fresh seekable reader
/// (discarding the re-parsed header — the caller already holds the
/// validated one), so the reader is positioned past the header and
/// any chunk seek lands at a record boundary.
pub(super) fn open_bam_reader_with_header(
    path: &Path,
) -> Result<(bam::io::Reader<bgzf::io::Reader<File>>, sam::Header), AlignmentInputError> {
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
    Ok((reader, header))
}

/// The `noodles` query interval for a `region`: the 1-based inclusive
/// range when `Some`, or the whole contig (unbounded) when `None`.
pub(super) fn query_interval(region: Option<ContigInterval>) -> Interval {
    match region {
        Some(iv) => {
            // `iv.start`/`iv.end` are >= 1 by construction, so the
            // `Position::new` calls cannot fail.
            let start = Position::new(iv.start as usize).expect("region start >= 1");
            let end = Position::new(iv.end as usize).expect("region end >= 1");
            Interval::from(start..=end)
        }
        None => Interval::from(..),
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use noodles_core::Position;
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

    /// Write a synthetic coordinate-sorted BAM with the given records to a
    /// temp dir; return the dir handle (so the path stays alive) and the path.
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
        // Empty body (header-only BAM): the iterator must return None on
        // every poll past the first, mirroring the no-latch contract.
        let header = header_two_contigs();
        let (_dir, bam_path) = write_bam(&[], &header);

        let (_header, mut stream) = open_bam_record_stream(&bam_path).expect("open empty bam");

        assert!(stream.next().is_none(), "empty BAM should yield no records");
        for _ in 0..3 {
            assert!(
                stream.next().is_none(),
                "post-EOF poll must not surface a spurious record or error"
            );
        }
    }

    /// The linear BAM iterator must surface a read error as `Some(Err(_))`
    /// rather than silently terminating with `None`. Trigger it by writing a
    /// small BAM then truncating it mid-stream so the BGZF reader fails.
    #[test]
    fn open_bam_record_stream_surfaces_read_errors_as_some_err() {
        use std::fs::OpenOptions;

        let header = header_two_contigs();
        let records: Vec<RecordBuf> = (0..200).map(|i| record_on(0, 1 + (i % 60))).collect();
        let (_dir, bam_path) = write_bam(&records, &header);

        let file = OpenOptions::new()
            .write(true)
            .open(&bam_path)
            .expect("open for truncate");
        file.set_len(4096).expect("truncate to 4 KiB");
        drop(file);

        let (_header, mut stream) = open_bam_record_stream(&bam_path).expect("open truncated bam");

        let mut saw_err = false;
        let mut polls = 0;
        for item in stream.by_ref() {
            polls += 1;
            if item.is_err() {
                saw_err = true;
                break;
            }
            if polls > 5000 {
                break;
            }
        }
        assert!(
            saw_err,
            "truncated BAM should surface a Some(Err(_)) (saw {polls} successful polls then None)"
        );
    }
}
