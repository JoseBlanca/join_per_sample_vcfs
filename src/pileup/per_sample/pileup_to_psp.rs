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

use std::collections::VecDeque;
use std::io::Write;

use thiserror::Error;

use crate::fasta::MultiChromRefFetcher;
use crate::pileup::walker::{PileupWalker, PreparedRead, RunSummary, WalkerError};
use crate::pileup_record::PileupRecord;
use crate::psp::PspWriteError;
use crate::psp::writer::PspWriter;
use crate::sample_summary::coverage::{
    CoverageBinScheme, SlidingWindowCoverageAccumulator, WindowCoverage,
};
use crate::sample_summary::het::{
    AlleleGroupStats, HetAccumulator, HetClassifyParams, SiteCounts,
};
use crate::sample_summary::{SAMPLE_SUMMARY_VERSION, SampleSummary};

/// The REF base of a record — its REF allele's first byte, or `N` for the
/// impossible empty-seq case (GC/coverage-excluded, never a window centre).
fn record_ref_base(record: &PileupRecord) -> u8 {
    record
        .alleles
        .first()
        .and_then(|a| a.seq.first())
        .copied()
        .unwrap_or(b'N')
}

/// Bundles the per-sample summary accumulators (the sliding-window
/// coverage-by-GC fold and observed heterozygosity) **and** the write-side
/// look-ahead buffer that populates each record's per-position windowed
/// coverage before it is written.
///
/// The seam sees every written record in coordinate order. It feeds each one to
/// the accumulators (turning it into a per-position depth/GC observation and,
/// for variant sites, `(ref, alt)` counts) and stages it for writing. Because
/// the centred window around a position closes only once the stream has read
/// `window_bp / 2` past it, a record cannot be written the moment it arrives —
/// it is held in `pending` until its [`WindowCoverage`] is finalised, then
/// written with its `windowed_gc` / `windowed_coverage` filled in. The bundle
/// (and its buffer) persists across regions — a record near a region's end has
/// its window close in the next region — and is drained + reduced once, by the
/// caller, before `PspWriter::finish` (via [`finish`](Self::finish)).
pub struct SampleSummaryAccumulators {
    coverage: SlidingWindowCoverageAccumulator,
    het: HetAccumulator,
    /// Records staged but not yet written (their windows still open), in
    /// coordinate order. Bounded by `~window_bp` near the stream frontier.
    pending: VecDeque<PileupRecord>,
    /// Finalised window values awaiting correlation with `pending`, in position
    /// order — one per non-`N` covered position.
    ready: VecDeque<WindowCoverage>,
}

impl SampleSummaryAccumulators {
    /// Build from the coverage bin scheme and het classifier parameters
    /// (CLI-supplied; see `pileup`'s `--gc-window-bp` and the
    /// `sample_summary` defaults).
    pub fn new(coverage: CoverageBinScheme, het: HetClassifyParams) -> Self {
        Self {
            coverage: SlidingWindowCoverageAccumulator::new(coverage),
            het: HetAccumulator::new(het),
            pending: VecDeque::new(),
            ready: VecDeque::new(),
        }
    }

    /// Fold one record into the summaries and stage it for writing. Call once
    /// per record, in coordinate order. Records reach `writer` in the same
    /// order, lagged by up to the window half-width (a record is written once
    /// its centred window closes).
    pub fn stage_record<W: Write>(
        &mut self,
        record: PileupRecord,
        writer: &mut PspWriter<W>,
    ) -> Result<(), PspWriteError> {
        // Het: one candidate site per record that carries a non-reference
        // allele. Pure-REF columns (one allele) are skipped.
        if let Some(ref_allele) = record.alleles.first()
            && record.alleles.len() > 1
        {
            // REF allele group; ALT groups pooled over every non-REF allele.
            let s = &ref_allele.support;
            let reference = AlleleGroupStats {
                obs: u64::from(s.num_obs),
                fwd: u64::from(s.fwd),
                log_error_sum: s.q_sum,
            };
            let mut alt = AlleleGroupStats::default();
            for a in &record.alleles[1..] {
                let s = &a.support;
                alt.obs += u64::from(s.num_obs);
                alt.fwd += u64::from(s.fwd);
                alt.log_error_sum += s.q_sum;
            }
            self.het.observe_site(SiteCounts { reference, alt });
        }

        // Coverage: ONE observation at the record's **anchor** position. The
        // walker emits a separate per-position record for every other covered
        // reference base, so attributing a multi-base (deletion / MNP) record's
        // whole REF span here would double-count its interior. Total fragment
        // depth = Σ observations over every allele (REF + ALTs).
        let depth: u32 = record.alleles.iter().map(|a| a.support.num_obs).sum();
        self.coverage
            .observe(record.chrom_id, record.pos, record_ref_base(&record), depth);
        while let Some(wc) = self.coverage.pop_ready() {
            self.ready.push_back(wc);
        }

        self.pending.push_back(record);
        flush_ready(&mut self.pending, &mut self.ready, writer)
    }

    /// Finalise: flush the remaining staged records (their windows truncated at
    /// the stream end) and reduce to the per-sample [`SampleSummary`] for the
    /// `.psp` metadata section. Call once, after the last `stage_record`, before
    /// `PspWriter::attach_metadata` / `finish`.
    pub fn finish<W: Write>(
        self,
        writer: &mut PspWriter<W>,
    ) -> Result<SampleSummary, PspWriteError> {
        let Self {
            coverage,
            het,
            mut pending,
            mut ready,
        } = self;
        let (tail, coverage_by_gc) = coverage.finish();
        ready.extend(tail);
        flush_ready(&mut pending, &mut ready, writer)?;
        debug_assert!(
            pending.is_empty(),
            "every staged record must be written at finish ({} left)",
            pending.len(),
        );
        Ok(SampleSummary {
            version: SAMPLE_SUMMARY_VERSION,
            coverage_by_gc,
            heterozygosity: het.finish(),
        })
    }
}

/// Write every `pending` record whose per-position window is resolved, oldest
/// first, until the front record's window is still open. An `N`-ref record has
/// no window (it never becomes a centre) and is written immediately with its
/// `NaN` placeholder; a non-`N` record is written once its [`WindowCoverage`] is
/// `ready`, with `windowed_gc` / `windowed_coverage` filled in.
///
/// The correlation is exact: `ready` holds one window per non-`N` position in
/// coordinate order, and non-`N` records consume them in the same order (`N`
/// records consume none), so when the front record is non-`N` the front of
/// `ready` — if present — is its window (debug-asserted on position).
fn flush_ready<W: Write>(
    pending: &mut VecDeque<PileupRecord>,
    ready: &mut VecDeque<WindowCoverage>,
    writer: &mut PspWriter<W>,
) -> Result<(), PspWriteError> {
    while let Some(front) = pending.front() {
        if record_ref_base(front).eq_ignore_ascii_case(&b'N') {
            // N position: never a window centre → write with the NaN placeholder.
            let record = pending.pop_front().expect("peeked non-empty above");
            writer.write_record(&record)?;
        } else if ready.front().is_some() {
            let wc = ready.pop_front().expect("peeked Some above");
            let mut record = pending.pop_front().expect("peeked non-empty above");
            // Hard invariant (release-active, not debug-only): the front non-N
            // record pairs with the front window, on the **same contig and
            // position**. A mismatch means the pending/ready correlation drifted
            // — abort rather than silently write a wrong per-position coverage
            // into the `.psp`. `chrom_id` is load-bearing: positions restart at
            // 1 per contig, so a `pos`-only check would miss a cross-contig
            // mis-pair. The cost is two `u32` compares per record.
            assert_eq!(
                (wc.chrom_id, wc.pos),
                (record.chrom_id, record.pos),
                "window/record correlation drift",
            );
            record.windowed_gc = wc.gc_fraction;
            record.windowed_coverage = wc.mean_depth;
            writer.write_record(&record)?;
        } else {
            // The front record's centred window has not closed yet.
            break;
        }
    }
    Ok(())
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
///
/// The clamp filters on `record.pos` alone: the caller supplies **one contig
/// per region** (each region drives a walker built for a single contig), so a
/// record's `chrom_id` is implied by the region and a position filter is
/// sufficient. Records from more than one contig must never reach a single
/// call.
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
            // Feed the summaries and stage the record; it is written (with its
            // windowed coverage filled in) once its centred window closes —
            // possibly during a later region or at `SampleSummaryAccumulators::finish`.
            summary.stage_record(record, writer)?;
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

        // `PileupRecord`'s `PartialEq` is byte-identity (windowed floats compared
        // on bits), so the walker's `NaN` windowed placeholders round-trip and
        // compare equal.
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
    fn multibase_record_anchor_only_no_double_count() {
        use crate::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};

        let rec = |pos: u32, ref_seq: &[u8], obs: u32| PileupRecord {
            chrom_id: 0,
            pos,
            windowed_gc: f32::NAN,
            windowed_coverage: f32::NAN,
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
        let mut writer = PspWriter::new(Cursor::new(Vec::<u8>::new()), writer_header(1)).unwrap();
        // Wide REF span at pos 2 (4 bp), then a narrow record at pos 3. Coverage
        // is attributed to each record's anchor only (pos 2, pos 3) — the wide
        // record's interior is covered by the walker's own per-position records,
        // not re-counted here — so the coordinate-ordered accumulator never sees
        // a backwards position and each of the two covered positions is one
        // window centre.
        acc.stage_record(rec(2, b"ACGT", 6), &mut writer).unwrap();
        acc.stage_record(rec(3, b"C", 6), &mut writer).unwrap();
        let h = acc.finish(&mut writer).unwrap().coverage_by_gc;
        // Two covered positions (2 and 3), each finalised into one window
        // sample — no double-count of the wide record's interior.
        assert_eq!(
            h.n_positions, 2,
            "two covered positions → two window samples"
        );
        assert_eq!(h.callable_positions, 2);
    }

    /// Test-default summary accumulators (tiny GC window so the few test
    /// positions span more than one tile).
    fn test_summary_accumulators() -> SampleSummaryAccumulators {
        use crate::sample_summary::{
            DEFAULT_DEPTH_BIN_WIDTH, DEFAULT_DEPTH_BINS, DEFAULT_GC_BINS, DEFAULT_HET_ERROR_RATE,
            DEFAULT_HET_LR_MARGIN, DEFAULT_HET_MIN_DEPTH, DEFAULT_HET_STRAND_BIAS_Z,
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
                strand_bias_z: DEFAULT_HET_STRAND_BIAS_Z,
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

        // Reduce the shared accumulator (flushing the last staged records) and
        // attach it as the metadata section before finishing — the production
        // finalisation order.
        let summary_doc = summary
            .finish(&mut writer)
            .expect("finish flushes + reduces");
        // Four covered positions (2,3,4,5), each one window sample; no ALT
        // alleles -> no variant sites.
        assert_eq!(summary_doc.coverage_by_gc.n_positions, 4);
        assert_eq!(summary_doc.coverage_by_gc.callable_positions, 4);
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

    /// End-to-end seam: staged records are written in coordinate order, lagged
    /// by the window half-width, with each record's `windowed_coverage` filled
    /// in from the centred sliding window. An `N` reference position emits no
    /// window (it never becomes a centre and contributes to no neighbour's
    /// window) and is written with the `NaN` sentinel. `window_bp = 3` → the
    /// window centred on `p` spans covered positions in `[p-1, p+1]`.
    #[test]
    fn seam_fills_windowed_coverage_and_skips_n_positions() {
        use crate::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};

        let rec = |pos: u32, base: u8, depth: u32| PileupRecord {
            chrom_id: 0,
            pos,
            windowed_gc: f32::NAN,
            windowed_coverage: f32::NAN,
            alleles: vec![AlleleObservation {
                seq: vec![base],
                support: AlleleSupportStats {
                    num_obs: depth,
                    ..Default::default()
                },
                chain_ids: Vec::new(),
            }],
        };

        let mut acc = test_summary_accumulators(); // window_bp = 3 → half = 1
        let mut writer = PspWriter::new(Cursor::new(Vec::<u8>::new()), writer_header(1)).unwrap();
        // Positions 1,2,4,5 covered; position 3 is an `N` reference base (huge
        // depth, which must NOT leak into any window mean).
        for (pos, base, depth) in [
            (1u32, b'A', 10u32),
            (2, b'C', 20),
            (3, b'N', 999),
            (4, b'G', 40),
            (5, b'T', 50),
        ] {
            acc.stage_record(rec(pos, base, depth), &mut writer)
                .unwrap();
        }
        let summary = acc.finish(&mut writer).unwrap();
        // Four GC-defined covered positions (the N is excluded).
        assert_eq!(summary.coverage_by_gc.callable_positions, 4);

        let bytes = writer.finish().unwrap().into_inner();
        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let got: Vec<_> = reader.records().map(|r| r.unwrap()).collect();

        // Written in coordinate order, all five positions present.
        let positions: Vec<u32> = got.iter().map(|r| r.pos).collect();
        assert_eq!(positions, vec![1, 2, 3, 4, 5]);

        // Centred means over `[p-1, p+1]` ∩ covered (pos 3 excluded):
        //   pos1 {1,2}=15, pos2 {1,2}=15 (3 is N), pos4 {4,5}=45, pos5 {4,5}=45.
        let cov = |p: u32| got.iter().find(|r| r.pos == p).unwrap().windowed_coverage;
        assert!((cov(1) - 15.0).abs() < 1e-6, "pos1 {}", cov(1));
        assert!((cov(2) - 15.0).abs() < 1e-6, "pos2 {}", cov(2));
        assert!(cov(3).is_nan(), "N position has no window");
        assert!((cov(4) - 45.0).abs() < 1e-6, "pos4 {}", cov(4));
        assert!((cov(5) - 45.0).abs() < 1e-6, "pos5 {}", cov(5));
    }

    /// A single-`REF`-allele record at `pos` on `chrom` with reference base
    /// `base` and total depth `depth` (walker placeholders for the windowed
    /// fields).
    fn ref_record(chrom: u32, pos: u32, base: u8, depth: u32) -> PileupRecord {
        use crate::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};
        PileupRecord {
            chrom_id: chrom,
            pos,
            windowed_gc: f32::NAN,
            windowed_coverage: f32::NAN,
            alleles: vec![AlleleObservation {
                seq: vec![base],
                support: AlleleSupportStats {
                    num_obs: depth,
                    ..Default::default()
                },
                chain_ids: Vec::new(),
            }],
        }
    }

    /// The look-ahead buffer resolves each record's window on its **own** contig:
    /// when a new contig's first record arrives it finalises the previous
    /// contig's still-open windows, and `flush_ready` pairs each pending record
    /// with a window from the same contig — never averaging or mis-pairing across
    /// the boundary. (The last contig-0 record has an open window when contig 1
    /// starts — the boundary case the `chrom_id`-aware pairing guard protects.)
    #[test]
    fn seam_pairs_windows_within_contig_across_contig_change() {
        let mut acc = test_summary_accumulators(); // window_bp = 3 → half = 1
        let mut writer = PspWriter::new(Cursor::new(Vec::<u8>::new()), writer_header(2)).unwrap();
        for (c, p, d) in [
            (0u32, 1u32, 10u32),
            (0, 2, 10),
            (0, 3, 10),
            (1, 1, 100),
            (1, 2, 100),
            (1, 3, 100),
        ] {
            acc.stage_record(ref_record(c, p, b'G', d), &mut writer)
                .unwrap();
        }
        acc.finish(&mut writer).unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let got: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
        assert_eq!(got.len(), 6);
        // No cross-contig depth bleed: every contig-0 window mean is 10, every
        // contig-1 window mean is 100 (a mis-pair would pull the other contig's
        // depth in).
        for r in &got {
            let expected = if r.chrom_id == 0 { 10.0 } else { 100.0 };
            assert!(
                (r.windowed_coverage - expected).abs() < 1e-6,
                "chrom {} pos {} coverage {} bled across the contig boundary",
                r.chrom_id,
                r.pos,
                r.windowed_coverage,
            );
        }
    }

    /// A long all-`N` region does not accumulate in `pending`: `N` records are
    /// written immediately (never blocking on a window) and advance the
    /// finalisation frontier, so the covered flanks resolve and the buffer
    /// drains fully — `finish`'s `debug_assert!(pending.is_empty())` would fire
    /// otherwise. All records are written in coordinate order.
    #[test]
    fn seam_all_n_region_drains_and_does_not_grow_pending() {
        let mut acc = test_summary_accumulators(); // window_bp = 3 → half = 1
        let mut writer = PspWriter::new(Cursor::new(Vec::<u8>::new()), writer_header(1)).unwrap();
        // Covered flank, then 200 consecutive `N` positions (a run far longer
        // than window_bp), then a covered flank. The N records must write
        // immediately (never blocking on a window) and advance the frontier, so
        // the buffer never grows past ~window_bp and drains fully at finish.
        acc.stage_record(ref_record(0, 1, b'G', 10), &mut writer)
            .unwrap();
        for pos in 2..=201u32 {
            acc.stage_record(ref_record(0, pos, b'N', 999), &mut writer)
                .unwrap();
        }
        acc.stage_record(ref_record(0, 202, b'G', 100), &mut writer)
            .unwrap();
        // `finish`'s `debug_assert!(pending.is_empty())` fires if the run left
        // records stuck in the buffer.
        let summary = acc.finish(&mut writer).unwrap();
        // Only the two covered flanks are GC-defined covered positions (the N
        // run contributes nothing).
        assert_eq!(summary.coverage_by_gc.callable_positions, 2);

        let bytes = writer.finish().unwrap().into_inner();
        let mut reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let got: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
        assert_eq!(got.len(), 202);
        assert!(
            got.windows(2).all(|w| w[0].pos < w[1].pos),
            "coordinate order",
        );
        // Flanks carry finite window means; the N interior carries NaN.
        assert!((got[0].windowed_coverage - 10.0).abs() < 1e-6);
        assert!(got[100].windowed_coverage.is_nan(), "N interior → NaN");
        assert!((got[201].windowed_coverage - 100.0).abs() < 1e-6);
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
