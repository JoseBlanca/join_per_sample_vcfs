//! The production seam between the pileup walker and the `.psp`
//! writer.
//!
//! [`drive_pileup_to_psp`] pulls every record from a
//! [`PileupWalker`] and feeds it
//! through a [`PspWriter`], finalising the writer when the
//! walker exhausts. Errors from either side surface through the
//! combined [`PileupToPspError`] enum.
//!
//! This is the in-tool Stage 1 emission half — the wiring from
//! per-position records to the on-disk artefact downstream stages
//! consume. The CRAM input → BAQ → walker plumbing is upstream of
//! this and lives in sibling slices.

use std::io::Write;

use thiserror::Error;

use crate::fasta::MultiChromRefFetcher;
use crate::pileup::walker::{PileupWalker, PreparedRead, RunSummary, WalkerError};
use crate::pileup_record::PileupRecord;
use crate::psp::PspWriteError;
use crate::psp::writer::PspWriter;
use crate::sample_summary::coverage::{CoverageBinScheme, CoverageByGcAccumulator};
use crate::sample_summary::het::{HetAccumulator, HetClassifyParams, SiteCounts};
use crate::sample_summary::{SAMPLE_SUMMARY_VERSION, SampleSummary};

/// Bundles the two per-sample summary accumulators (coverage-by-GC and
/// observed heterozygosity) and owns the **record → counts shaping**:
/// the seam already sees every written record in coordinate order, so it
/// turns each one into the per-position depth/GC observations and the
/// per-variant-site `(ref, alt)` counts the accumulators consume. The
/// accumulators themselves do only the math
/// ([`crate::sample_summary`]). The bundle persists across regions (one
/// per sample run) and is reduced once, by the caller, before
/// `PspWriter::finish`.
pub struct SampleSummaryAccumulators {
    coverage: CoverageByGcAccumulator,
    het: HetAccumulator,
}

impl SampleSummaryAccumulators {
    /// Build from the coverage bin scheme and het classifier parameters
    /// (CLI-supplied; see `pileup`'s `--gc-window-bp` and the
    /// `sample_summary` defaults).
    pub fn new(coverage: CoverageBinScheme, het: HetClassifyParams) -> Self {
        Self {
            coverage: CoverageByGcAccumulator::new(coverage),
            het: HetAccumulator::new(het),
        }
    }

    /// Fold one finalised record into both summaries. Call once per
    /// written record, in coordinate order.
    pub fn observe_record(&mut self, record: &PileupRecord) {
        // Every record carries ≥ 1 allele with REF first (a `PileupRecord`
        // invariant); guard defensively so a malformed empty-allele record
        // is skipped rather than index-panicking.
        let Some(ref_allele) = record.alleles.first() else {
            return;
        };
        // Total fragment depth at the site = Σ observations over every
        // allele (REF + ALTs).
        let depth: u32 = record.alleles.iter().map(|a| a.support.num_obs).sum();

        // Coverage: ONE observation at the record's **anchor** position.
        // The walker emits a separate per-position record for every other
        // covered reference base, so attributing the record's whole REF
        // span here would double-count a multi-base (deletion / MNP)
        // record's interior — which the following per-position records
        // already cover — and feed the coordinate-ordered accumulator a
        // backwards position. The REF base is the first byte of the REF
        // allele's sequence (`N` for the impossible empty seq → GC-excluded).
        let ref_base = ref_allele.seq.first().copied().unwrap_or(b'N');
        self.coverage
            .observe(record.chrom_id, record.pos, ref_base, depth);

        // Het: one candidate site per record that carries a non-reference
        // allele. The binomial variant gate then decides whether it is a
        // real variant. Pure-REF columns (one allele) are skipped.
        if record.alleles.len() > 1 {
            let ref_obs = u64::from(ref_allele.support.num_obs);
            let alt_obs: u64 = record.alleles[1..]
                .iter()
                .map(|a| u64::from(a.support.num_obs))
                .sum();
            self.het.observe_site(SiteCounts { ref_obs, alt_obs });
        }
    }

    /// Reduce to the per-sample [`SampleSummary`] for the `.psp` metadata
    /// section.
    pub fn finish(self) -> SampleSummary {
        SampleSummary {
            version: SAMPLE_SUMMARY_VERSION,
            coverage_by_gc: self.coverage.finish(),
            heterozygosity: self.het.finish(),
        }
    }
}

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

/// Drive one region's [`PileupWalker`] into a **shared** [`PspWriter`],
/// writing only records whose position lies in the 1-based inclusive
/// range `[start, end]`. Returns the walker's [`RunSummary`].
///
/// Unlike [`drive_pileup_to_psp`], the writer is borrowed and **not**
/// finalised here: the region-driven pileup calls this once per region
/// against one shared writer, then finalises after the last region.
///
/// The walker runs over reads that *overlap* the region — including
/// reads starting before `start`, so BAQ keeps its flanking context —
/// so it can emit columns just outside the range; the clamp drops
/// those, leaving exactly the region's columns. Because analysis
/// regions are disjoint, every emitted column is written by exactly one
/// region, so concatenating the per-region writes yields a coordinate-
/// sorted, duplicate-free record stream.
pub fn drive_region_into_writer<I, F, W>(
    mut walker: PileupWalker<I, F>,
    writer: &mut PspWriter<W>,
    start: u32,
    end: u32,
    summary: &mut SampleSummaryAccumulators,
) -> Result<RunSummary, PileupToPspError>
where
    I: Iterator<Item = PreparedRead>,
    F: MultiChromRefFetcher,
    W: Write,
{
    for item in walker.by_ref() {
        let record = item?;
        // `start`/`end` are the region's 1-based inclusive bounds. The
        // summary accumulators see exactly the records that are written,
        // so the stored statistics match the `.psp` body.
        if (start..=end).contains(&record.pos) {
            summary.observe_record(&record);
            writer.write_record(&record)?;
        }
    }
    Ok(walker.summary())
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use crate::pileup::walker::tests::{MockFasta, snp_read};
    use crate::pileup::walker::{WalkerConfig, run};
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

    /// A multi-base (e.g. deletion-anchor) record at `pos` followed by a
    /// narrow record at `pos+1` — the shape that exposed the span-iteration
    /// coverage bug — must fold in coordinate order without double-counting.
    /// Anchor-only attribution means each record contributes one coverage
    /// observation at its own anchor; the walker's separate per-position
    /// records cover the wide record's interior.
    #[test]
    fn observe_record_handles_multibase_record_in_coordinate_order() {
        use crate::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};

        let rec = |pos: u32, ref_seq: &[u8], obs: u32| PileupRecord {
            chrom_id: 0,
            pos,
            alleles: vec![AlleleObservation {
                seq: ref_seq.to_vec(),
                support: AlleleSupportStats {
                    num_obs: obs,
                    ..Default::default()
                },
                chain_ids: Vec::new(),
            }],
        };

        let mut acc = test_summary_accumulators();
        // Wide REF span at pos 2 (4 bp), then a narrow record at pos 3.
        // The pre-fix span loop would `observe(2,3,4,5)` then `observe(3)`,
        // tripping the coordinate-order assert; anchor-only does not.
        acc.observe_record(&rec(2, b"ACGT", 6));
        acc.observe_record(&rec(3, b"C", 6));
        let h = acc.finish().coverage_by_gc;
        // Exactly two anchor observations (pos 2, pos 3); both at depth 6,
        // so no double-counting of the wide record's interior.
        assert_eq!(h.n_tiles, 1, "positions 2 and 3 fall in one window-3 tile");
    }

    /// Test-default summary accumulators (tiny GC window so the few test
    /// positions span more than one tile).
    fn test_summary_accumulators() -> SampleSummaryAccumulators {
        use crate::sample_summary::{
            DEFAULT_DEPTH_BIN_WIDTH, DEFAULT_DEPTH_BINS, DEFAULT_GC_BINS, DEFAULT_HET_ERROR_RATE,
            DEFAULT_HET_LR_MARGIN, DEFAULT_HET_MIN_DEPTH,
        };
        SampleSummaryAccumulators::new(
            CoverageBinScheme {
                window_bp: 3,
                gc_bins: DEFAULT_GC_BINS,
                depth_bin_width: DEFAULT_DEPTH_BIN_WIDTH,
                depth_bins: DEFAULT_DEPTH_BINS,
            },
            HetClassifyParams {
                min_depth: DEFAULT_HET_MIN_DEPTH,
                error_rate: DEFAULT_HET_ERROR_RATE,
                lr_margin: DEFAULT_HET_LR_MARGIN,
            },
        )
    }

    /// Two regions driven into one shared writer, each clamped to its
    /// range, produce a single `.psp` carrying exactly the in-range
    /// positions in coordinate order. The shared summary accumulator sees
    /// exactly the written records and its reduction attaches as the
    /// metadata section. Mirrors how the region-driven pileup reuses one
    /// writer + one accumulator across regions.
    #[test]
    fn drive_region_into_writer_clamps_and_shares_one_writer() {
        use crate::sample_summary::SampleSummary;

        let mut writer = PspWriter::new(Cursor::new(Vec::<u8>::new()), writer_header(1))
            .expect("writer construction");
        let mut summary = test_summary_accumulators();

        // Region A: positions [2, 4]. Reads cover the whole 1..=5 ref.
        let fa_a = MockFasta::new("ACGTA");
        let walker_a = run(
            vec![snp_read("r1", 1, b"ACGTA", &[30; 5])],
            &fa_a,
            &WalkerConfig::default(),
        );
        drive_region_into_writer(walker_a, &mut writer, 2, 4, &mut summary)
            .expect("region A drives");

        // Region B: a single position [5, 5] (disjoint from A).
        let fa_b = MockFasta::new("ACGTA");
        let walker_b = run(
            vec![snp_read("r2", 1, b"ACGTA", &[30; 5])],
            &fa_b,
            &WalkerConfig::default(),
        );
        drive_region_into_writer(walker_b, &mut writer, 5, 5, &mut summary)
            .expect("region B drives");

        // Reduce the shared accumulator and attach it as the metadata
        // section before finishing — the production finalisation order.
        let summary_doc = summary.finish();
        // Four written REF columns (2,3,4,5), all covered once -> tiles
        // carry depth; no ALT alleles -> no variant sites.
        assert!(summary_doc.coverage_by_gc.n_tiles >= 1);
        assert_eq!(summary_doc.heterozygosity.n_variant_sites, 0);
        writer
            .attach_metadata(summary_doc.to_toml_bytes().expect("serialise summary"))
            .expect("attach");

        let bytes = writer.finish().expect("finish").into_inner();
        let mut reader = PspReader::new(Cursor::new(bytes)).expect("reader open");
        // The metadata section round-trips back to the same summary.
        let parsed = SampleSummary::from_toml_bytes(reader.metadata().expect("section present"))
            .expect("parse summary");
        assert_eq!(parsed, summary_doc);
        let positions: Vec<u32> = reader
            .records()
            .map(|r| r.expect("record decodes").pos)
            .collect();
        // A contributes 2,3,4; B contributes 5 — clamped and in order,
        // no out-of-range (1) position and no duplicates.
        assert_eq!(positions, vec![2, 3, 4, 5]);
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
