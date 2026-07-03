//! S5 — the write pass: spill → read stored LR → apply cut → VCF.
//!
//! The **single** spill read. Each locus's LR was scored **once in the main
//! pass** — in the caller worker (`ParalogScoringContext::score`), stored on its
//! spill record (`ParalogSpillRecord::paralog_lr`) and folded into the histogram
//! the calibration was derived from (arch: inline scoring).
//! This pass streams the spill, looks up each **stored** LR's FDR q-value on the
//! calibration, and either **drops** the locus (a paralog, counted in
//! [`WriterStats::records_dropped_paralog`]) or routes the surviving record back
//! through the unchanged [`VcfWriter`] — so the existing downstream filters
//! (hom-ref / min-QUAL / allele-balance) and the byte-identity-critical VCF
//! formatting run exactly as in the single-pass path. The run's parameters go
//! in the VCF header ([`paralog_provenance`]) since a dropped call leaves no
//! per-record trace (arch §7).
//!
//! Applying the cut to the *stored* LR (rather than re-scoring) makes the cut
//! match the histogram it came from **by construction**, and holds only one
//! record at a time — RAM stays flat in the variant count.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;

use crate::fasta::{ChromRefFetchError, ChromRefFetcher, StreamingChromRefFetcher};
use crate::var_calling::posterior_engine::PosteriorRecord;
use crate::var_calling::vcf_writer::{FormatOutcome, VcfWriter, WriterError, WriterStats};

use super::calibrate::ParalogCalibration;
use super::spill::{ParalogSpillReader, SpillError};

/// Failure modes of the write pass.
#[derive(Debug, thiserror::Error)]
pub enum WritePassError {
    /// Reading the spill back failed.
    #[error("spill read failed")]
    Spill(#[from] SpillError),
    /// Writing the surviving VCF failed.
    #[error("VCF write failed")]
    Writer(#[from] WriterError),
    /// Fetching the reference base for the coordinate-consistency guard failed.
    #[error("reference fetch during the paralog write pass")]
    RefFetch(#[from] ChromRefFetchError),
    /// The FASTA base at a surviving locus disagreed with the decoded REF allele
    /// — the window-key / coordinate arithmetic has drifted (an off-by-one
    /// against the FASTA), so the GC lookups the score used were silently wrong.
    /// Fail loud rather than emit a mis-scored callset.
    #[error(
        "reference/decoded REF mismatch at chrom_id {chrom_id} pos {pos}: \
         FASTA '{fasta}' vs decoded REF '{decoded}' — window/coordinate drift"
    )]
    ReferenceMismatch {
        chrom_id: u32,
        pos: u32,
        fasta: char,
        decoded: char,
    },
}

/// Coordinate-consistency guard (owner's ask): the FASTA base at a surviving
/// locus must match the decoded REF allele's first base — the `.psp` REF
/// allele's first base *is* the reference base there, so a mismatch means the
/// window-key / coordinate arithmetic has drifted and the GC lookups that fed
/// the paralog score were silently wrong. Fails loud.
///
/// Rebuilds the per-contig fetcher on a chromosome change; survivors arrive in
/// strictly-increasing `(chrom_id, pos)` order, satisfying the streaming
/// fetcher's monotonic-forward contract.
fn check_reference_consistency(
    chrom_id: u32,
    locus_start: u32,
    decoded: Option<u8>,
    fetcher: &mut Option<(u32, StreamingChromRefFetcher)>,
    reference: &Path,
    chrom_names: &[String],
) -> Result<(), WritePassError> {
    // An alleleless record never reaches the writer, but stay total.
    let Some(decoded) = decoded else {
        return Ok(());
    };
    if fetcher.as_ref().map(|(c, _)| *c) != Some(chrom_id) {
        let name = &chrom_names[chrom_id as usize];
        *fetcher = Some((
            chrom_id,
            StreamingChromRefFetcher::for_contig(reference, name)?,
        ));
    }
    // UNREACHABLE: `fetcher` is `Some` on both branches above.
    let (_, f) = fetcher.as_ref().expect("fetcher just set for this contig");
    let fasta_base = f.fetch(locus_start, 1)?.first().copied().unwrap_or(b'N');
    if !reference_base_matches(fasta_base, decoded) {
        return Err(WritePassError::ReferenceMismatch {
            chrom_id,
            pos: locus_start,
            fasta: fasta_base as char,
            decoded: decoded as char,
        });
    }
    Ok(())
}

/// Whether the reference base the FASTA walk reports at a locus matches the
/// REF allele the caller decoded there — the coordinate-consistency guard used
/// by [`check_reference_consistency`].
///
/// The `.psp` REF allele's first base *is* the reference base at that position,
/// so a mismatch means the coordinate arithmetic has drifted (e.g. an off-by-one
/// against the FASTA). Compared case-insensitively; an `N` on either side is
/// treated as a match (soft-masked / ambiguous reference is not a coordinate
/// error).
fn reference_base_matches(fasta_base: u8, decoded_ref_base: u8) -> bool {
    fasta_base.eq_ignore_ascii_case(&b'N')
        || decoded_ref_base.eq_ignore_ascii_case(&b'N')
        || fasta_base.eq_ignore_ascii_case(&decoded_ref_base)
}

/// The `##paralogFilter=...` header provenance for a calibration: the target
/// FDR, the estimated `π`, the resolved LR cut, and whether the EM converged.
/// A reader can tell the filter ran and with what settings (the dropped records
/// carry no trace).
pub(crate) fn paralog_provenance(calibration: &ParalogCalibration) -> String {
    let cut = match calibration.lr_threshold {
        Some(lr) => format!("{lr:.4}"),
        None => "none".to_string(),
    };
    format!(
        "target_fdr={:.4};pi={:.6};lr_cut={};em_converged={}",
        calibration.target_fdr,
        calibration.prior.prior_probability,
        cut,
        calibration.prior.converged,
    )
}

/// Depth (per worker) of both write-pass hand-off queues. Bounds resident
/// records to `O(workers)` so the spill's flat-RAM property survives; `2` keeps
/// every worker fed across one hand-off.
const WRITE_PASS_QUEUE_DEPTH_PER_WORKER: usize = 2;

/// Stream the spill, drop the paralog-flagged loci, and write the survivors —
/// **parallel-format pipeline**. Returns the writer's stats with
/// [`WriterStats::records_dropped_paralog`] set.
///
/// The serial predecessor decoded, filtered, formatted, and bgzf-wrote each
/// survivor on one thread; profiling put ~70% of that in the *VCF text format*
/// (`encode` + noodles serialisation), which is pure per record. So this splits
/// into three stages inside a [`std::thread::scope`], hand-offs bounded for
/// back-pressure (RAM stays flat in the variant count):
///
/// 1. **reader** (one thread): decode each spill frame, stamp the paralog
///    posterior, apply the stored-LR cut, and forward survivors — tagged with a
///    contiguous `order` in spill (= genomic) order — to the workers;
/// 2. **`W` format workers**: run the downstream filters (hom-ref → min-QUAL →
///    allele-balance, order load-bearing) and serialise each survivor to its VCF
///    line ([`ParallelRecordFormatter`](crate::var_calling::vcf_writer::ParallelRecordFormatter)),
///    off the writer thread;
/// 3. **sink** (this thread): reorder by `order`, run the coordinate-consistency
///    guard, and commit the pre-formatted bytes to the bgzf sink in genomic order
///    — byte-identical VCF (same `encode`, same header, same order) and identical
///    run-summary counters to the serial path.
///
/// Each locus's LR was scored once in the main pass (arch: inline scoring); this
/// pass only applies the stored cut, no re-score.
pub(crate) fn run_write_pass<R: Read + Send>(
    spill: &mut ParalogSpillReader<R>,
    calibration: &ParalogCalibration,
    keep_dup_artifacts: bool,
    reference: &Path,
    chrom_names: &[String],
    mut writer: VcfWriter,
    worker_count: usize,
) -> Result<WriterStats, WritePassError> {
    let n_workers = worker_count.max(1);
    let cap = (WRITE_PASS_QUEUE_DEPTH_PER_WORKER * n_workers).max(1);
    let (work_tx, work_rx) = crossbeam_channel::bounded::<(u64, PosteriorRecord)>(cap);
    let (done_tx, done_rx) = crossbeam_channel::bounded::<(u64, FormatOutcome)>(cap);

    let records_dropped_paralog = std::thread::scope(
        |scope| -> Result<u64, WritePassError> {
            // Reader: decode + paralog cut, contiguous `order` to survivors.
            let reader = scope.spawn(move || -> Result<u64, WritePassError> {
                let mut order = 0u64;
                let mut dropped = 0u64;
                while let Some(record) = spill.next_record() {
                    let mut record = record?;
                    // Apply the cut to the stored LR (folded into the histogram this
                    // `calibration` came from — so the cut matches the value it was
                    // built from). `NaN` = unscored → never flagged → kept.
                    let flagged = calibration.flags(record.paralog_lr);
                    // Stamp the paralog posterior on every scored locus (drop or
                    // keep) so emitted records carry `PARALOG_POST`; `NaN` → `None`.
                    record.record.paralog_posterior = calibration.posterior(record.paralog_lr);
                    if flagged && !keep_dup_artifacts {
                        dropped += 1;
                        continue;
                    }
                    if work_tx.send((order, record.record)).is_err() {
                        break; // workers gone
                    }
                    order += 1;
                }
                Ok(dropped)
            });

            // Format workers: filter + serialise, off the writer thread.
            let mut worker_handles = Vec::with_capacity(n_workers);
            for _ in 0..n_workers {
                let work_rx = work_rx.clone();
                let done_tx = done_tx.clone();
                let mut formatter = writer.make_formatter();
                worker_handles.push(scope.spawn(move || -> Result<(), WritePassError> {
                    for (order, record) in work_rx {
                        let outcome = formatter.process(&record).map_err(WriterError::from)?;
                        if done_tx.send((order, outcome)).is_err() {
                            break; // sink gone
                        }
                    }
                    Ok(())
                }));
            }
            // Drop the main-thread channel handles so the channels disconnect once
            // the reader / workers finish (else the sink loop never ends).
            drop(work_rx);
            drop(done_tx);

            // Sink (this thread): reorder by `order`, guard, commit in genomic order.
            let mut ref_fetcher: Option<(u32, StreamingChromRefFetcher)> = None;
            let mut pending: BTreeMap<u64, FormatOutcome> = BTreeMap::new();
            let mut next = 0u64;
            let sink_result = (|| -> Result<(), WritePassError> {
                for (order, outcome) in done_rx.iter() {
                    pending.insert(order, outcome);
                    while let Some(outcome) = pending.remove(&next) {
                        emit_outcome(
                            outcome,
                            &mut writer,
                            &mut ref_fetcher,
                            reference,
                            chrom_names,
                        )?;
                        next += 1;
                    }
                }
                Ok(())
            })();

            // Join everyone; first error wins (sink → reader → workers). The scope
            // would auto-join, but we join explicitly to surface thread errors.
            let reader_result = reader.join().expect("write-pass reader thread panicked");
            let mut worker_err = None;
            for h in worker_handles {
                if let Err(e) = h.join().expect("write-pass format worker panicked") {
                    worker_err.get_or_insert(e);
                }
            }
            sink_result?;
            let dropped = reader_result?;
            if let Some(e) = worker_err {
                return Err(e);
            }
            Ok(dropped)
        },
    )?;

    let mut stats = writer.finish()?;
    stats.records_dropped_paralog = records_dropped_paralog;
    Ok(stats)
}

/// Commit one worker verdict on the sink thread: tally a drop, or (for a
/// survivor) run the coordinate-consistency guard and write the pre-formatted
/// bytes to the sink in genomic order.
fn emit_outcome(
    outcome: FormatOutcome,
    writer: &mut VcfWriter,
    ref_fetcher: &mut Option<(u32, StreamingChromRefFetcher)>,
    reference: &Path,
    chrom_names: &[String],
) -> Result<(), WritePassError> {
    match outcome {
        FormatOutcome::Drop(bucket) => writer.tally_drop(bucket),
        FormatOutcome::Emit {
            bytes,
            chrom_id,
            pos_1based,
            locus_start,
            decoded_ref,
            unconverged,
        } => {
            check_reference_consistency(
                chrom_id,
                locus_start,
                decoded_ref,
                ref_fetcher,
                reference,
                chrom_names,
            )?;
            if unconverged {
                writer.tally_unconverged();
            }
            writer.write_preformatted_line(&bytes, chrom_id, pos_1based)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::paralog::{ParalogLrHistogram, ParalogModelParams};
    use crate::psp::header::ParsedChromosome;
    use crate::var_calling::paralog_filter::calibrate::{
        CalibrationConfig, ParalogScoringContext, calibrate_from_histogram,
    };
    use crate::var_calling::paralog_filter::spill::{ParalogSpillReader, ParalogSpillRecord};
    use crate::var_calling::paralog_filter::test_support::{
        normal_locus, paralog_locus, prepass, write_spill,
    };
    use crate::var_calling::vcf_writer::{AlleleBalanceFilter, DownstreamFilters};
    use crate::vcf::{CohortMetadata, WriterConfig};
    use std::io::Cursor;
    use tempfile::tempdir;

    const COHORT_F: f64 = 0.02;

    /// Fill each record's `paralog_lr` as the worker does — score it from its own
    /// window coverage — and build the calibration from those LRs (the histogram
    /// the sink folds inline), exactly the production path. Mutates `records` in
    /// place so `write_spill(records)` then carries the stored LR the write pass
    /// reads. Returns the calibration for the target `fdr`.
    fn score_and_calibrate(
        records: &mut [ParalogSpillRecord],
        n: usize,
        fdr: f64,
    ) -> ParalogCalibration {
        let ctx = ParalogScoringContext::new(prepass(n), COHORT_F, &ParalogModelParams::default());
        let mut buf = Vec::new();
        let mut histogram = ParalogLrHistogram::with_defaults();
        for record in records.iter_mut() {
            let lr = ctx.score(&record.record, &record.window_coverage, &mut buf);
            record.paralog_lr = lr;
            histogram.push(lr);
        }
        calibrate_from_histogram(&histogram, fdr, &CalibrationConfig::default())
    }

    /// The single contig name the fixtures live on (matches [`metadata`]).
    fn chrom_names() -> Vec<String> {
        vec!["chr1".to_string()]
    }

    /// Write a minimal FASTA + `.fai` for contig `chr1`: `len` copies of `base`.
    /// The test loci all use REF `A`, so an all-`A` reference matches (and an
    /// all-`C` one triggers the coordinate-consistency guard).
    fn reference_fasta(len: usize, base: u8) -> (tempfile::TempDir, std::path::PathBuf) {
        use std::io::Write as _;
        let dir = tempdir().unwrap();
        let fasta = dir.path().join("ref.fa");
        let fai = dir.path().join("ref.fa.fai");
        let header = ">chr1\n";
        let mut fa = std::fs::File::create(&fasta).unwrap();
        fa.write_all(header.as_bytes()).unwrap();
        fa.write_all(&vec![base; len]).unwrap();
        fa.write_all(b"\n").unwrap();
        // name\tlength\toffset\tline_bases\tline_width
        std::fs::write(
            &fai,
            format!(
                "chr1\t{len}\t{off}\t{len}\t{width}\n",
                off = header.len(),
                width = len + 1
            ),
        )
        .unwrap();
        (dir, fasta)
    }

    fn metadata(n_samples: usize, provenance: String) -> CohortMetadata {
        CohortMetadata {
            sample_names: (0..n_samples).map(|s| format!("S{s}")).collect(),
            contigs: vec![ParsedChromosome {
                name: "chr1".into(),
                length: 1_000_000,
                md5: "0123456789abcdef0123456789abcdef".into(),
            }],
            tool_string: "pop_var_caller paralog-write-test".into(),
            command_line: "pop_var_caller cohort".into(),
            paralog_provenance: provenance,
        }
    }

    fn no_ab_filters() -> DownstreamFilters {
        DownstreamFilters {
            min_qual_phred: 30.0,
            allele_balance: AlleleBalanceFilter::Off,
        }
    }

    /// The write pass drops exactly the paralog-flagged loci and writes the rest;
    /// the surviving VCF has the kept positions, and the header records the
    /// provenance.
    #[test]
    fn drops_flagged_writes_survivors_and_records_provenance() {
        let n = 20;
        let mut records = Vec::new();
        for i in 0..80u32 {
            records.push(normal_locus(1000 + i, n));
        }
        for i in 0..20u32 {
            records.push(paralog_locus(5000 + i, n));
        }
        let cal = score_and_calibrate(&mut records, n, 0.05);
        let provenance = paralog_provenance(&cal);
        assert!(provenance.contains("target_fdr=0.05"));

        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let writer = VcfWriter::new(
            metadata(n, provenance.clone()),
            WriterConfig::new(out.clone()),
            no_ab_filters(),
        )
        .unwrap();

        let (_ref_dir, reference) = reference_fasta(6000, b'A');
        let stats = run_write_pass(
            &mut ParalogSpillReader::new(Cursor::new(write_spill(&records))),
            &cal,
            false, // keep_dup_artifacts: default drop mode
            &reference,
            &chrom_names(),
            writer,
            2, // format workers (exercise the parallel path)
        )
        .expect("write pass");

        // Every paralog-like locus is flagged; every normal locus survives.
        assert_eq!(stats.records_dropped_paralog, 20);
        assert_eq!(stats.records_written, 80);

        let text = std::fs::read_to_string(&out).unwrap();
        // Header provenance present.
        assert!(
            text.contains("##paralogFilter=") && text.contains("target_fdr=0.05"),
            "header must record the paralog provenance"
        );
        // The 20 paralog positions (5000..5020) are absent; a normal one is present.
        let data_positions: Vec<&str> = text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').nth(1).unwrap())
            .collect();
        assert_eq!(data_positions.len(), 80);
        assert!(data_positions.contains(&"1000"));
        assert!(!data_positions.iter().any(|p| p.starts_with("500")));
    }

    /// `--do-not-drop-dup-artifacts` keeps every flagged locus (nothing dropped)
    /// and stamps the `PARALOG_POST` posterior on the scored records, so the
    /// artifact class is present and annotated instead of removed.
    #[test]
    fn keep_dup_artifacts_retains_flagged_and_annotates_posterior() {
        let n = 20;
        let mut records = Vec::new();
        for i in 0..80u32 {
            records.push(normal_locus(1000 + i, n));
        }
        for i in 0..20u32 {
            records.push(paralog_locus(5000 + i, n));
        }
        let cal = score_and_calibrate(&mut records, n, 0.05);

        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        // Min-qual 0 so ONLY the paralog keep/drop decision governs — the
        // skewed-AB paralog loci would otherwise trip the refined-QUAL min-qual
        // filter downstream (the research callset pairs this flag with
        // `--min-qual 0` for exactly that reason).
        let filters = DownstreamFilters {
            min_qual_phred: 0.0,
            allele_balance: AlleleBalanceFilter::Off,
        };
        let writer = VcfWriter::new(
            metadata(n, paralog_provenance(&cal)),
            WriterConfig::new(out.clone()),
            filters,
        )
        .unwrap();

        let (_ref_dir, reference) = reference_fasta(6000, b'A');
        let stats = run_write_pass(
            &mut ParalogSpillReader::new(Cursor::new(write_spill(&records))),
            &cal,
            true, // keep_dup_artifacts: annotate but do not drop
            &reference,
            &chrom_names(),
            writer,
            2, // format workers (exercise the parallel path)
        )
        .expect("write pass");

        // Nothing dropped; all 100 loci written (80 normal + 20 flagged kept).
        assert_eq!(stats.records_dropped_paralog, 0);
        assert_eq!(stats.records_written, 100);

        let text = std::fs::read_to_string(&out).unwrap();
        // The PARALOG_POST INFO field is declared and emitted on scored loci.
        assert!(
            text.contains("##INFO=<ID=PARALOG_POST"),
            "header must declare PARALOG_POST"
        );
        let data: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data.len(), 100);
        // The flagged paralog positions (5000..) are RETAINED here.
        assert!(data.iter().any(|l| l.split('\t').nth(1) == Some("5000")));
        // Every scored biallelic-SNP record carries the posterior.
        assert!(
            data.iter().all(|l| l.contains("PARALOG_POST=")),
            "each scored locus must carry PARALOG_POST"
        );
    }

    /// A reference that disagrees with the decoded REF allele fails the write
    /// pass loud (coordinate/window drift), rather than emitting a mis-scored
    /// callset. Here every survivor's REF is `A` but the FASTA is all `C`.
    #[test]
    fn reference_mismatch_fails_loud() {
        let n = 10;
        let mut records: Vec<ParalogSpillRecord> =
            (0..5u32).map(|i| normal_locus(1000 + i, n)).collect();
        let cal = score_and_calibrate(&mut records, n, -1.0); // flags nothing → all survive

        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let writer = VcfWriter::new(
            metadata(n, paralog_provenance(&cal)),
            WriterConfig::new(out),
            no_ab_filters(),
        )
        .unwrap();
        let (_ref_dir, reference) = reference_fasta(2000, b'C'); // REF is A → mismatch
        let err = run_write_pass(
            &mut ParalogSpillReader::new(Cursor::new(write_spill(&records))),
            &cal,
            false, // keep_dup_artifacts: default drop mode
            &reference,
            &chrom_names(),
            writer,
            2, // format workers (exercise the parallel path)
        )
        .expect_err("mismatch must fail loud");
        match err {
            WritePassError::ReferenceMismatch { fasta, decoded, .. } => {
                assert_eq!(fasta, 'C');
                assert_eq!(decoded, 'A');
            }
            other => panic!("expected ReferenceMismatch, got {other:?}"),
        }
    }

    /// With the filter's cut unachievable (nothing flagged), the write pass
    /// writes every record — same as the pre-filter callset.
    #[test]
    fn nothing_flagged_writes_everything() {
        let n = 10;
        let mut records: Vec<ParalogSpillRecord> =
            (0..30u32).map(|i| normal_locus(1000 + i, n)).collect();
        // A negative FDR target is unachievable → flags nothing.
        let cal = score_and_calibrate(&mut records, n, -1.0);

        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let writer = VcfWriter::new(
            metadata(n, paralog_provenance(&cal)),
            WriterConfig::new(out),
            no_ab_filters(),
        )
        .unwrap();

        let (_ref_dir, reference) = reference_fasta(2000, b'A');
        let stats = run_write_pass(
            &mut ParalogSpillReader::new(Cursor::new(write_spill(&records))),
            &cal,
            false, // keep_dup_artifacts: default drop mode
            &reference,
            &chrom_names(),
            writer,
            2, // format workers (exercise the parallel path)
        )
        .unwrap();
        assert_eq!(stats.records_dropped_paralog, 0);
        assert_eq!(stats.records_written, 30);
    }

    /// An empty spill writes nothing and finishes cleanly (no stranded reorder
    /// chunk) — the zero-record boundary of the gapless-`chunk_order` invariant.
    #[test]
    fn writes_nothing_on_empty_spill() {
        let n = 4;
        let cal = score_and_calibrate(&mut [normal_locus(1, n)], n, 0.05);
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let writer = VcfWriter::new(
            metadata(n, paralog_provenance(&cal)),
            WriterConfig::new(out),
            no_ab_filters(),
        )
        .unwrap();
        let (_ref_dir, reference) = reference_fasta(100, b'A');
        let stats = run_write_pass(
            &mut ParalogSpillReader::new(Cursor::new(Vec::<u8>::new())),
            &cal,
            false, // keep_dup_artifacts: default drop mode
            &reference,
            &chrom_names(),
            writer,
            2, // format workers (exercise the parallel path)
        )
        .unwrap();
        assert_eq!(stats.records_written, 0);
        assert_eq!(stats.records_dropped_paralog, 0);
    }

    /// The write pass is deterministic: the same spill + calibration produce the
    /// identical VCF (the cut is a pure curve lookup on the stored LR).
    #[test]
    fn write_pass_is_deterministic() {
        let n = 15;
        let mut records = Vec::new();
        for i in 0..40u32 {
            records.push(normal_locus(1000 + i, n));
        }
        for i in 0..10u32 {
            records.push(paralog_locus(5000 + i, n));
        }
        let cal = score_and_calibrate(&mut records, n, 0.1);

        let (_ref_dir, reference) = reference_fasta(6000, b'A');
        let run = || {
            let dir = tempdir().unwrap();
            let out = dir.path().join("out.vcf");
            let writer = VcfWriter::new(
                metadata(n, paralog_provenance(&cal)),
                WriterConfig::new(out.clone()),
                no_ab_filters(),
            )
            .unwrap();
            run_write_pass(
                &mut ParalogSpillReader::new(Cursor::new(write_spill(&records))),
                &cal,
                false, // keep_dup_artifacts: default drop mode
                &reference,
                &chrom_names(),
                writer,
                2, // format workers (exercise the parallel path)
            )
            .unwrap();
            std::fs::read_to_string(&out).unwrap()
        };
        assert_eq!(run(), run(), "same inputs → identical VCF");
    }

    /// The coordinate-consistency guard matches on equal bases (any case) and on
    /// an `N` either side, and rejects a genuine mismatch.
    #[test]
    fn reference_base_consistency_guard() {
        assert!(reference_base_matches(b'A', b'A'));
        assert!(reference_base_matches(b'a', b'A')); // case-insensitive
        assert!(reference_base_matches(b'N', b'A')); // N reference is a match
        assert!(reference_base_matches(b'G', b'n')); // N allele is a match
        assert!(!reference_base_matches(b'A', b'C')); // genuine mismatch
    }
}
