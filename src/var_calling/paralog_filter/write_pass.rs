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

use std::io::Read;
use std::path::Path;

use crate::fasta::{ChromRefFetchError, ChromRefFetcher, StreamingChromRefFetcher};
use crate::var_calling::posterior_engine::PosteriorRecord;
use crate::var_calling::types::{CallStats, CalledChunk};
use crate::var_calling::vcf_writer::{VcfWriter, WriterError, WriterStats};

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
    record: &PosteriorRecord,
    fetcher: &mut Option<(u32, StreamingChromRefFetcher)>,
    reference: &Path,
    chrom_names: &[String],
) -> Result<(), WritePassError> {
    // An alleleless record never reaches the writer, but stay total.
    let Some(&decoded) = record.alleles.first().and_then(|a| a.seq.first()) else {
        return Ok(());
    };
    let chrom_id = record.locus.chrom_id;
    if fetcher.as_ref().map(|(c, _)| *c) != Some(chrom_id) {
        let name = &chrom_names[chrom_id as usize];
        *fetcher = Some((
            chrom_id,
            StreamingChromRefFetcher::for_contig(reference, name)?,
        ));
    }
    // UNREACHABLE: `fetcher` is `Some` on both branches above.
    let (_, f) = fetcher.as_ref().expect("fetcher just set for this contig");
    let fasta_base = f
        .fetch(record.locus.start, 1)?
        .first()
        .copied()
        .unwrap_or(b'N');
    if !reference_base_matches(fasta_base, decoded) {
        return Err(WritePassError::ReferenceMismatch {
            chrom_id,
            pos: record.locus.start,
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

/// Stream the spill, drop the paralog-flagged loci, and write the survivors
/// through `writer`. Returns the writer's stats with
/// [`WriterStats::records_dropped_paralog`] set.
///
/// Each locus's LR was **scored once in the main pass** (the worker) and stored
/// on the spill ([`ParalogSpillRecord::paralog_lr`]), and folded into the
/// histogram `calibration` was derived from — so this pass just reads the stored
/// LR and applies the cut, with no re-score (arch: inline scoring). Survivors are
/// fed to the writer in spill order (which is genomic order — the main pass
/// spills post-reorder), one record per `CalledChunk` with a gapless
/// `chunk_order`; the writer's reorder passes them straight through.
pub(crate) fn run_write_pass<R: Read>(
    spill: &mut ParalogSpillReader<R>,
    calibration: &ParalogCalibration,
    reference: &Path,
    chrom_names: &[String],
    mut writer: VcfWriter,
) -> Result<WriterStats, WritePassError> {
    let mut records_dropped_paralog = 0u64;
    let mut chunk_order = 0u64;
    // Per-contig reference fetcher for the coordinate-consistency guard, built
    // lazily on the first survivor of each contig.
    let mut ref_fetcher: Option<(u32, StreamingChromRefFetcher)> = None;

    while let Some(record) = spill.next_record() {
        let record = record?;
        // Apply the cut to the locus's stored LR (scored once in the worker,
        // folded into the histogram `calibration` came from — so the cut applies
        // to the exact value the calibration was built from). An unscored locus
        // carries `NaN` and is never flagged → kept.
        let flagged = calibration.flags(record.paralog_lr);

        if flagged {
            records_dropped_paralog += 1;
            continue;
        }

        // Survivor: verify the coordinate is consistent with the reference
        // before emitting (fails loud on window/coordinate drift).
        check_reference_consistency(&record.record, &mut ref_fetcher, reference, chrom_names)?;

        // Survivor: route through the unchanged writer (existing filters +
        // byte-identity-critical formatting). One record per chunk keeps the
        // reorder cursor gapless from 0.
        writer.handle(CalledChunk {
            chunk_order,
            // The VCF writer ignores the window coverage + LR (they rode the spill
            // only for the score, already applied above); empty vecs suffice here.
            window_coverage: Vec::new(),
            paralog_lr: Vec::new(),
            records: vec![record.record],
            stats: CallStats::default(),
        })?;
        chunk_order += 1;
    }

    let mut stats = writer.finish()?;
    stats.records_dropped_paralog = records_dropped_paralog;
    Ok(stats)
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
            &reference,
            &chrom_names(),
            writer,
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
            &reference,
            &chrom_names(),
            writer,
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
            &reference,
            &chrom_names(),
            writer,
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
            &reference,
            &chrom_names(),
            writer,
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
                &reference,
                &chrom_names(),
                writer,
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
