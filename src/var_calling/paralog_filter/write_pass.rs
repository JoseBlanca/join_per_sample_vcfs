//! S5 — the write pass: spill → recompute LR → apply cut → VCF.
//!
//! The second spill read. It streams the spill a second time, **recomputes**
//! each locus's LR from the same spilled inputs via the shared
//! [`score_spilled_locus`] (bit-identical to the calibrate pass, since the
//! scorer is pure — arch §3), looks up its FDR q-value on the calibration, and
//! either **drops** the locus (a paralog, counted in
//! [`WriterStats::records_dropped_paralog`]) or routes the surviving record back
//! through the unchanged [`VcfWriter`] — so the existing downstream filters
//! (hom-ref / min-QUAL / allele-balance) and the byte-identity-critical VCF
//! formatting run exactly as in the single-pass path. The run's parameters go
//! in the VCF header ([`paralog_provenance`]) since a dropped call leaves no
//! per-record trace (arch §7).
//!
//! Recompute-not-cache is the non-negotiable invariant: holding the calibrate
//! pass's LRs (or records) to skip this second read would make RAM grow with the
//! variant count — the exact thing the spill exists to prevent. This pass holds
//! one record + an `O(samples)` scratch buffer at a time.

use std::io::Read;
use std::path::Path;

use crate::fasta::{ChromRefFetchError, ChromRefFetcher, StreamingChromRefFetcher};
use crate::paralog::{ParalogModelParams, ParalogScorePrecompute};
use crate::var_calling::posterior_engine::PosteriorRecord;
use crate::var_calling::types::{CallStats, CalledChunk};
use crate::var_calling::vcf_writer::{VcfWriter, WriterError, WriterStats};

use super::calibrate::{ParalogCalibration, cohort_inbreeding, score_joined_locus};
use super::prepass::ParalogPrePass;
use super::spill::{ParalogSpillReader, SpillError, WindowJoin};
use super::window_gc::reference_base_matches;

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
        *fetcher = Some((chrom_id, StreamingChromRefFetcher::for_contig(reference, name)?));
    }
    // UNREACHABLE: `fetcher` is `Some` on both branches above.
    let (_, f) = fetcher.as_ref().expect("fetcher just set for this contig");
    let fasta_base = f.fetch(record.locus.start, 1)?.first().copied().unwrap_or(b'N');
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
/// `prepass`, `inbreeding_coefficient`, `window_bp`, and `params` must be the
/// **same** as the calibrate pass used (both `inbreeding` and
/// `single_copy_depth_sd` are derived from `prepass`), and `windows` must be a
/// fresh joiner over the **same** window spill, so the recomputed LRs match the
/// histogram the cut came from. Survivors are fed to the writer in spill order
/// (which is genomic order — the main pass spills post-reorder), one record per
/// `CalledChunk` with a gapless `chunk_order`; the writer's reorder passes them
/// straight through.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_write_pass<R: Read, WR: Read>(
    spill: &mut ParalogSpillReader<R>,
    windows: &mut WindowJoin<WR>,
    window_bp: u32,
    prepass: &ParalogPrePass,
    inbreeding_coefficient: f64,
    calibration: &ParalogCalibration,
    params: &ParalogModelParams,
    reference: &Path,
    chrom_names: &[String],
    mut writer: VcfWriter,
) -> Result<WriterStats, WritePassError> {
    let single_copy_depth_sd = prepass.single_copy_depth_sd();
    let inbreeding = cohort_inbreeding(single_copy_depth_sd.len(), inbreeding_coefficient);
    // Built once from the same params + cohort inbreeding the calibrate pass
    // used, so the recomputed LRs are bit-identical to the histogram's.
    let precompute = ParalogScorePrecompute::new(params, &inbreeding);
    let mut obs_buf = Vec::new();
    let mut records_dropped_paralog = 0u64;
    let mut chunk_order = 0u64;
    // Per-contig reference fetcher for the coordinate-consistency guard, built
    // lazily on the first survivor of each contig.
    let mut ref_fetcher: Option<(u32, StreamingChromRefFetcher)> = None;

    while let Some(record) = spill.next_record() {
        let record = record?;
        // Recompute the LR from the same joined inputs (bit-identical to
        // calibrate) and drop the locus if the calibration flags it. An unscored
        // locus (not a biallelic SNP, no usable samples, or no joined window) is
        // never flagged → kept.
        let flagged = score_joined_locus(
            &record.record,
            windows,
            window_bp,
            prepass,
            &inbreeding,
            &single_copy_depth_sd,
            &precompute,
            &mut obs_buf,
        )?
        .is_some_and(|lr| calibration.flags(lr));

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
            // The VCF writer ignores window coverage (it rode the spill only for
            // the score, already applied above); an empty vec suffices here.
            window_coverage: Vec::new(),
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
    use crate::psp::header::ParsedChromosome;
    use crate::var_calling::paralog_filter::calibrate::{CalibrationConfig, calibrate};
    use crate::var_calling::paralog_filter::spill::{
        ParalogSpillReader, ParalogSpillRecord, WindowSpillReader, WindowSpillRecord,
    };
    use crate::var_calling::paralog_filter::test_support::{
        TEST_WINDOW_BP, normal_locus, normal_window, paralog_locus, paralog_window, prepass,
        write_spill, write_window_spill,
    };
    use crate::var_calling::vcf_writer::{AlleleBalanceFilter, DownstreamFilters};
    use crate::vcf::{CohortMetadata, WriterConfig};
    use std::io::Cursor;
    use tempfile::tempdir;

    const COHORT_F: f64 = 0.02;

    /// A fresh [`WindowJoin`] over the given window fixtures.
    fn join(windows: &[WindowSpillRecord]) -> WindowJoin<Cursor<Vec<u8>>> {
        WindowJoin::new(WindowSpillReader::new(Cursor::new(write_window_spill(windows)))).unwrap()
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
            format!("chr1\t{len}\t{off}\t{len}\t{width}\n", off = header.len(), width = len + 1),
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

    fn calibration_for(
        records: &[ParalogSpillRecord],
        windows: &[WindowSpillRecord],
        n: usize,
        fdr: f64,
    ) -> ParalogCalibration {
        calibrate(
            &prepass(n),
            COHORT_F,
            &mut ParalogSpillReader::new(Cursor::new(write_spill(records))),
            &mut join(windows),
            TEST_WINDOW_BP,
            &ParalogModelParams::default(),
            fdr,
            &CalibrationConfig::default(),
        )
        .expect("calibrate")
    }

    /// The write pass drops exactly the paralog-flagged loci and writes the rest;
    /// the surviving VCF has the kept positions, and the header records the
    /// provenance.
    #[test]
    fn drops_flagged_writes_survivors_and_records_provenance() {
        let n = 20;
        let mut records = Vec::new();
        let mut windows = Vec::new();
        for i in 0..80u32 {
            records.push(normal_locus(1000 + i, n));
            windows.push(normal_window(1000 + i, n));
        }
        for i in 0..20u32 {
            records.push(paralog_locus(5000 + i, n));
            windows.push(paralog_window(5000 + i, n));
        }
        let cal = calibration_for(&records, &windows, n, 0.05);
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
            &mut join(&windows),
            TEST_WINDOW_BP,
            &prepass(n),
            COHORT_F,
            &cal,
            &ParalogModelParams::default(),
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
        let records: Vec<ParalogSpillRecord> =
            (0..5u32).map(|i| normal_locus(1000 + i, n)).collect();
        let windows: Vec<WindowSpillRecord> =
            (0..5u32).map(|i| normal_window(1000 + i, n)).collect();
        let cal = calibration_for(&records, &windows, n, -1.0); // flags nothing → all survive

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
            &mut join(&windows),
            TEST_WINDOW_BP,
            &prepass(n),
            COHORT_F,
            &cal,
            &ParalogModelParams::default(),
            &reference,
            &chrom_names(),
            writer,
        )
        .expect_err("mismatch must fail loud");
        match err {
            WritePassError::ReferenceMismatch {
                fasta, decoded, ..
            } => {
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
        let records: Vec<ParalogSpillRecord> =
            (0..30u32).map(|i| normal_locus(1000 + i, n)).collect();
        let windows: Vec<WindowSpillRecord> =
            (0..30u32).map(|i| normal_window(1000 + i, n)).collect();
        // A negative FDR target is unachievable → flags nothing.
        let cal = calibration_for(&records, &windows, n, -1.0);

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
            &mut join(&windows),
            TEST_WINDOW_BP,
            &prepass(n),
            COHORT_F,
            &cal,
            &ParalogModelParams::default(),
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
        let cal = calibration_for(&[normal_locus(1, n)], &[normal_window(1, n)], n, 0.05);
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
            &mut join(&[]),
            TEST_WINDOW_BP,
            &prepass(n),
            COHORT_F,
            &cal,
            &ParalogModelParams::default(),
            &reference,
            &chrom_names(),
            writer,
        )
        .unwrap();
        assert_eq!(stats.records_written, 0);
        assert_eq!(stats.records_dropped_paralog, 0);
    }

    /// The write pass is deterministic: the same spill + calibration produce the
    /// identical VCF (the LR recompute is pure).
    #[test]
    fn write_pass_is_deterministic() {
        let n = 15;
        let mut records = Vec::new();
        let mut windows = Vec::new();
        for i in 0..40u32 {
            records.push(normal_locus(1000 + i, n));
            windows.push(normal_window(1000 + i, n));
        }
        for i in 0..10u32 {
            records.push(paralog_locus(5000 + i, n));
            windows.push(paralog_window(5000 + i, n));
        }
        let cal = calibration_for(&records, &windows, n, 0.1);

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
                &mut join(&windows),
                TEST_WINDOW_BP,
                &prepass(n),
                COHORT_F,
                &cal,
                &ParalogModelParams::default(),
                &reference,
                &chrom_names(),
                writer,
            )
            .unwrap();
            std::fs::read_to_string(&out).unwrap()
        };
        assert_eq!(run(), run(), "same inputs → identical VCF");
    }
}
