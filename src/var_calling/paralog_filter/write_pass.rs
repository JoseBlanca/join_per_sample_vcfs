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

use crate::paralog::ParalogModelParams;
use crate::var_calling::types::{CallStats, CalledChunk};
use crate::var_calling::vcf_writer::{VcfWriter, WriterError, WriterStats};

use super::calibrate::{ParalogCalibration, inbreeding_by_sample, score_spilled_locus};
use super::prepass::ParalogPrePass;
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
/// `prepass`, `hexp`, `params`, and `min_samples` must be the **same** as the
/// calibrate pass used (both `inbreeding` and `single_copy_depth_sd` are derived
/// from `prepass`), so the recomputed LRs match the histogram the cut came from.
/// Survivors are fed to the writer in spill order (which is genomic order
/// — the main pass spills post-reorder), one record per `CalledChunk` with a
/// gapless `chunk_order`; the writer's reorder passes them straight through.
pub(crate) fn run_write_pass<R: Read>(
    spill: &mut ParalogSpillReader<R>,
    prepass: &ParalogPrePass,
    hexp: f64,
    calibration: &ParalogCalibration,
    params: &ParalogModelParams,
    min_samples: usize,
    mut writer: VcfWriter,
) -> Result<WriterStats, WritePassError> {
    let inbreeding = inbreeding_by_sample(prepass, hexp);
    let single_copy_depth_sd = prepass.single_copy_depth_sd();
    let mut obs_buf = Vec::new();
    let mut records_dropped_paralog = 0u64;
    let mut chunk_order = 0u64;

    while let Some(record) = spill.next_record() {
        let record = record?;
        // Recompute the LR from the spilled inputs (bit-identical to calibrate)
        // and drop the locus if the calibration flags it. An unscored locus
        // (not a biallelic SNP, too few samples) is never flagged → kept.
        let flagged = score_spilled_locus(
            &record,
            prepass,
            &inbreeding,
            &single_copy_depth_sd,
            params,
            min_samples,
            &mut obs_buf,
        )
        .is_some_and(|lr| calibration.flags(lr));

        if flagged {
            records_dropped_paralog += 1;
            continue;
        }

        // Survivor: route through the unchanged writer (existing filters +
        // byte-identity-critical formatting). One record per chunk keeps the
        // reorder cursor gapless from 0.
        writer.handle(CalledChunk {
            chunk_order,
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
    use crate::var_calling::paralog_filter::spill::ParalogSpillRecord;
    use crate::var_calling::paralog_filter::test_support::{
        normal_locus, paralog_locus, prepass, write_spill,
    };
    use crate::var_calling::vcf_writer::{AlleleBalanceFilter, DownstreamFilters};
    use crate::vcf::{CohortMetadata, WriterConfig};
    use std::io::Cursor;
    use tempfile::tempdir;

    const HEXP: f64 = 0.02;
    const MIN_SAMPLES: usize = 5;

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

    fn calibration_for(records: &[ParalogSpillRecord], n: usize, fdr: f64) -> ParalogCalibration {
        let bytes = write_spill(records);
        calibrate(
            &prepass(n),
            HEXP,
            &mut ParalogSpillReader::new(Cursor::new(bytes)),
            &ParalogModelParams::default(),
            fdr,
            &CalibrationConfig {
                min_samples: MIN_SAMPLES,
                ..Default::default()
            },
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
        for i in 0..80u32 {
            records.push(normal_locus(1000 + i, n));
        }
        for i in 0..20u32 {
            records.push(paralog_locus(5000 + i, n));
        }
        let cal = calibration_for(&records, n, 0.05);
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

        let bytes = write_spill(&records);
        let stats = run_write_pass(
            &mut ParalogSpillReader::new(Cursor::new(bytes)),
            &prepass(n),
            HEXP,
            &cal,
            &ParalogModelParams::default(),
            MIN_SAMPLES,
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

    /// With the filter's cut unachievable (nothing flagged), the write pass
    /// writes every record — same as the pre-filter callset.
    #[test]
    fn nothing_flagged_writes_everything() {
        let n = 10;
        let records: Vec<ParalogSpillRecord> =
            (0..30u32).map(|i| normal_locus(1000 + i, n)).collect();
        // A negative FDR target is unachievable → flags nothing.
        let cal = calibration_for(&records, n, -1.0);

        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let writer = VcfWriter::new(
            metadata(n, paralog_provenance(&cal)),
            WriterConfig::new(out),
            no_ab_filters(),
        )
        .unwrap();

        let bytes = write_spill(&records);
        let stats = run_write_pass(
            &mut ParalogSpillReader::new(Cursor::new(bytes)),
            &prepass(n),
            HEXP,
            &cal,
            &ParalogModelParams::default(),
            MIN_SAMPLES,
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
        let cal = calibration_for(&[normal_locus(1, n)], n, 0.05);
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let writer = VcfWriter::new(
            metadata(n, paralog_provenance(&cal)),
            WriterConfig::new(out),
            no_ab_filters(),
        )
        .unwrap();
        let stats = run_write_pass(
            &mut ParalogSpillReader::new(Cursor::new(Vec::<u8>::new())),
            &prepass(n),
            HEXP,
            &cal,
            &ParalogModelParams::default(),
            MIN_SAMPLES,
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
        for i in 0..40u32 {
            records.push(normal_locus(1000 + i, n));
        }
        for i in 0..10u32 {
            records.push(paralog_locus(5000 + i, n));
        }
        let cal = calibration_for(&records, n, 0.1);

        let run = || {
            let dir = tempdir().unwrap();
            let out = dir.path().join("out.vcf");
            let writer = VcfWriter::new(
                metadata(n, paralog_provenance(&cal)),
                WriterConfig::new(out.clone()),
                no_ab_filters(),
            )
            .unwrap();
            let bytes = write_spill(&records);
            run_write_pass(
                &mut ParalogSpillReader::new(Cursor::new(bytes)),
                &prepass(n),
                HEXP,
                &cal,
                &ParalogModelParams::default(),
                MIN_SAMPLES,
                writer,
            )
            .unwrap();
            std::fs::read_to_string(&out).unwrap()
        };
        assert_eq!(run(), run(), "same inputs → identical VCF");
    }
}
