//! The production seam between the pileup walker and the `.psp`
//! writer.
//!
//! [`drive_pileup_to_psp`] pulls every record from a
//! [`pileup::PileupWalker`] and feeds it through a
//! [`psp::writer::PspWriter`], finalising the writer when the
//! walker exhausts. Errors from either side surface through the
//! combined [`PileupToPspError`] enum.
//!
//! This is the in-tool Stage 1 emission half — the wiring from
//! per-position records to the on-disk artefact downstream stages
//! consume. The CRAM input → BAQ → walker plumbing is upstream of
//! this and lives in sibling slices.

use std::io::Write;

use thiserror::Error;

use super::pileup::{PileupWalker, PreparedRead, RunSummary, WalkerError};
use crate::fasta::MultiChromRefFetcher;
use crate::psp::PspWriteError;
use crate::psp::writer::PspWriter;

/// Failure modes for [`drive_pileup_to_psp`]. Either the walker
/// surfaced an error mid-stream, or the writer rejected a record /
/// failed to finalise. The typed wrapping keeps each side's
/// context intact so operator triage can still tell *where* the
/// failure originated.
#[derive(Error, Debug)]
pub enum PileupToPspError {
    #[error("pileup walker failed: {0}")]
    Walker(#[from] WalkerError),
    #[error("psp writer failed: {0}")]
    Psp(#[from] PspWriteError),
}

/// Drive a [`PileupWalker`] end-to-end into a [`PspWriter`].
///
/// Returns the writer's underlying sink (post-`finish`, ready for a
/// caller-side `BufWriter::into_inner` + `File::sync_all` step where
/// applicable — see the doc-comment on [`PspWriter::finish`]) and the
/// walker's cumulative [`RunSummary`].
///
/// On the first error from either side, iteration stops immediately:
/// a walker error short-circuits before `write_record`, and a writer
/// error short-circuits before pulling the next record. In the
/// writer-error case the walker's open state is dropped without
/// flushing — there is no partial `.psp` artefact left in a
/// consistent shape; callers that need recovery should re-run from
/// scratch.
pub fn drive_pileup_to_psp<I, F, W>(
    mut walker: PileupWalker<I, F>,
    mut writer: PspWriter<W>,
) -> Result<(W, RunSummary), PileupToPspError>
where
    I: Iterator<Item = PreparedRead>,
    F: MultiChromRefFetcher,
    W: Write,
{
    for item in walker.by_ref() {
        let record = item?;
        writer.write_record(&record)?;
    }
    let summary = walker.summary();
    let sink = writer.finish()?;
    Ok((sink, summary))
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::super::pileup::tests::{MockFasta, snp_read};
    use super::super::pileup::{WalkerConfig, run};
    use super::*;
    use crate::psp::PspReader;
    use crate::psp::test_fixtures::writer_header;
    use crate::psp::writer::PspWriter;

    /// Drive a tiny pileup → psp run through the seam, then read the
    /// resulting bytes back through `PspReader` and confirm
    /// per-record parity. This is the contract the seam exists to
    /// preserve: whatever the walker emits, the writer encodes
    /// losslessly, and the reader yields the same records.
    ///
    /// Note: under the unique-chain-id design this test passes for
    /// real now. The earlier shape of this test pinned the walker /
    /// writer contract mismatch that was the proximate motivation
    /// for switching to unique ids.
    #[test]
    fn drive_pileup_to_psp_roundtrips_records_through_in_memory_sink() {
        let fa = MockFasta::new("ACGTA");
        let reads = vec![
            snp_read("r1", 1, b"ACGTA", &[30; 5]),
            snp_read("r2", 1, b"ACGTA", &[30; 5]),
        ];

        let walker = run(reads, &fa, &WalkerConfig::default());
        let writer = PspWriter::new(Cursor::new(Vec::<u8>::new()), writer_header(1))
            .expect("writer construction");

        let (sink, summary) = drive_pileup_to_psp(walker, writer).expect("seam should run cleanly");
        assert_eq!(summary.records_emitted, 5, "expect five emitted records");

        // Re-open the in-memory `.psp` artefact and pull every record
        // back. Parity is per-record equality against an independent
        // walker run over the same inputs.
        let bytes = sink.into_inner();
        let mut reader = PspReader::new(Cursor::new(bytes)).expect("reader open");
        let read_back: Vec<_> = reader
            .records()
            .map(|r| r.expect("reader yielded error"))
            .collect();

        let fa2 = MockFasta::new("ACGTA");
        let expected: Vec<_> = run(
            vec![
                snp_read("r1", 1, b"ACGTA", &[30; 5]),
                snp_read("r2", 1, b"ACGTA", &[30; 5]),
            ],
            &fa2,
            &WalkerConfig::default(),
        )
        .map(|r| r.expect("walker yielded error"))
        .collect();

        assert_eq!(
            read_back, expected,
            "records read back must equal records the walker emitted",
        );
    }

    /// Walker error surfaces as [`PileupToPspError::Walker`], not as
    /// a writer error and not as a panic. The fixture forces the
    /// walker to fail by giving it a reference shorter than the
    /// read (which trips `WalkerError::Fasta` inside the open-record
    /// fetch path).
    #[test]
    fn walker_error_surfaces_as_walker_variant() {
        let fa = MockFasta::new("AC"); // too short to satisfy a 4-base read
        let reads = vec![snp_read("r", 1, b"ACGT", &[30; 4])];
        let walker = run(reads, &fa, &WalkerConfig::default());
        let writer = PspWriter::new(Cursor::new(Vec::<u8>::new()), writer_header(1))
            .expect("writer construction");
        let err = drive_pileup_to_psp(walker, writer).expect_err("must surface walker error");
        match err {
            PileupToPspError::Walker(WalkerError::Fasta { .. }) => {}
            other => panic!("expected Walker(Fasta), got {other:?}"),
        }
    }

    /// Writer error surfaces as [`PileupToPspError::Psp`]. The
    /// fixture forces the writer to fail by giving the walker a
    /// `chrom_id` the writer header does not declare — chrom 0 is
    /// the only chromosome in `writer_header(1)`, while the read
    /// is on chrom 7. The walker accepts (the mock FASTA has eight
    /// chromosomes); the writer's `validate_record` rejects on
    /// first call.
    #[test]
    fn writer_error_surfaces_as_psp_variant() {
        let fa = MockFasta::with_chromosomes(&["A", "A", "A", "A", "A", "A", "A", "ACGTA"]);
        let mut r = snp_read("r", 1, b"ACGTA", &[30; 5]);
        r.chrom_id = 7;
        let walker = run(vec![r], &fa, &WalkerConfig::default());
        let writer = PspWriter::new(Cursor::new(Vec::<u8>::new()), writer_header(1))
            .expect("writer construction");
        let err = drive_pileup_to_psp(walker, writer).expect_err("must surface writer error");
        assert!(
            matches!(err, PileupToPspError::Psp(_)),
            "expected Psp variant, got {err:?}",
        );
    }
}
