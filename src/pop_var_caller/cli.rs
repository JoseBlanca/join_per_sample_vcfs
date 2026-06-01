//! `pop_var_caller pileup` — the Stage 1 CLI orchestrator.
//!
//! Glues `AlignmentMergedReader` → `BaqStream` (or BAQ-bypass passthrough)
//! → pileup walker → `.psp` writer. Argument parsing is delegated to
//! clap; everything else is in [`run_pileup`].
//!
//! Plan: `ia/feature_implementation_plans/pop_var_caller_pileup_cli.md`.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};

use clap::{Args, Parser, Subcommand};
use thiserror::Error;

use crate::bam::alignment_input::{AlignmentMergedReaderConfig, FilterCounts};
use crate::bam::errors::AlignmentInputError;
use crate::baq::BaqConfig;
use crate::fasta::ContigList;
use crate::pileup::per_sample::baq_stream::BaqSkipCounts;
use crate::pileup::per_sample::pileup_to_psp::{PileupToPspError, drive_pileup_to_psp};
use crate::pileup::walker::{RunSummary, WalkerConfig};
use crate::pop_var_caller::common::{
    DEFAULT_BUFFERED_IO_CAPACITY, basename, configure_rayon_pool, format_md5_hex, rfc3339_now,
};
use crate::psp::header::{ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance};
use crate::psp::writer::{PspWriter, TARGET_BLOCK_BYTES};

pub mod error_bridge;
pub mod parsers;
pub mod shared_args;

// ---------------------------------------------------------------------
// Clap surface
// ---------------------------------------------------------------------

/// Top-level CLI for the `pop_var_caller` binary.
#[derive(Debug, Parser)]
#[command(name = "pop_var_caller", version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub cmd: PopVarCallerCommand,
}

#[derive(Debug, Subcommand)]
pub enum PopVarCallerCommand {
    /// Stage 1: run BAQ + pileup over one sample's CRAMs, emit a .psp.
    Pileup(PileupArgs),

    /// Stream a .psp as samtools-mpileup-style text plus a trailing
    /// column with per-allele aggregates PSP carries.
    PspToPileup(super::psp_to_pileup::PspToPileupArgs),

    /// Estimate per-sample contamination from cohort `.psp` files and
    /// emit a TOML artefact consumed by `var-calling`.
    EstimateContamination(super::estimate_contamination::EstimateContaminationArgs),

    /// Call SNPs across a cohort: `.psp` files → multi-sample VCF.
    VarCalling(super::var_calling::VarCallingArgs),
}

/// Arguments accepted by the `pileup` subcommand. The struct is the
/// authoritative knob list; the orchestrator translates it into module
/// configs. Defaults are pulled from each module's `DEFAULT_*` consts
/// so they stay in sync.
///
/// Stage 1 knobs (CRAM-input filters, BAQ HMM, pileup walker) live in
/// the flattened [`Stage1Args`](shared_args::Stage1Args) sub-struct so
/// the `var-calling-from-bam` subcommand reuses the same surface
/// — M10 from the 2026-05-19 cohort CLI review.
#[derive(Debug, Args, Clone)]
pub struct PileupArgs {
    // ===== Common flags (visible in `-h`) =====================
    /// Reference FASTA. A sibling `.fai` is required.
    #[arg(long)]
    pub reference: PathBuf,

    /// Output .psp path. Must be a regular file; stdout is not
    /// supported (the .psp format is random-access on read).
    #[arg(long)]
    pub output: PathBuf,

    /// Worker threads for the BAQ stage and any other rayon work.
    /// If omitted, rayon's default (all logical cores) is used.
    #[arg(long)]
    pub threads: Option<usize>,

    /// One or more coordinate-sorted CRAM or BAM file(s) for one
    /// sample. All inputs in one invocation must share the same
    /// format (CRAM-only or BAM-only); mixed CRAM + BAM is rejected
    /// with a typed error.
    #[arg(required = true)]
    pub alignment_files: Vec<PathBuf>,

    /// PSP writer auto-flush threshold (uncompressed bytes per block).
    /// Range 16 KiB..=16 MiB. Default 1 MiB.
    ///
    /// Trade-off — peak memory vs on-disk size. The cohort
    /// var-calling driver opens one `PspReader` per sample per
    /// chromosome worker, each holding one decoded block live, so
    /// peak heap scales as `n_threads × N × per_block`. Smaller
    /// blocks → lower cohort-step memory, slightly worse zstd
    /// compression context (so larger .psp on disk).
    ///
    /// Sweep on tomato1 (N=18, T=4): 16 MiB → 2501 MB peak / 183 MB
    /// on disk; 1 MiB → 261 MB / 206 MB (+13% disk, default);
    /// 256 KiB → 108 MB / 246 MB (+34% disk); 64 KiB → 59 MB /
    /// 321 MB (+75% disk, near the knee). Wall time is flat across
    /// the range — no CPU penalty.
    ///
    /// Dial down (256 KiB, 64 KiB) for large-cohort joint
    /// genotyping that hits a memory ceiling. Dial up (4 MiB,
    /// 16 MiB) for single-sample archival where on-disk size
    /// matters more than cohort-step heap.
    #[arg(
        long,
        default_value_t = TARGET_BLOCK_BYTES,
        value_parser = parsers::parse_block_target_bytes,
        help_heading = "Advanced — PSP writer",
    )]
    pub block_target_bytes: usize,

    // ===== Stage 1 (shared with var-calling-from-bam) =========
    #[command(flatten)]
    pub stage1: shared_args::Stage1Args,
}

pub(crate) fn parse_mismatch_fraction(s: &str) -> Result<f32, String> {
    let v: f32 = s.parse().map_err(|e| format!("not a number: {e}"))?;
    if !v.is_finite() {
        return Err(format!("must be finite, got `{s}`"));
    }
    if !(0.0..=1.0).contains(&v) {
        return Err(format!("must be in [0.0, 1.0], got `{s}`"));
    }
    Ok(v)
}

// ---------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum PileupCliError {
    #[error("alignment input: {0}")]
    AlignmentInput(#[from] AlignmentInputError),
    #[error("pipeline: {0}")]
    Pipeline(#[from] PileupToPspError),
    #[error(
        "contig '{contig}' has no @SQ M5 checksum in the input alignment file(s); \
         re-create the file with a tool that emits @SQ M5 (e.g. samtools view -t)"
    )]
    MissingMd5 { contig: String },
    #[error("rayon thread pool already initialised — refusing to override")]
    RayonAlreadyConfigured,
    #[error("io: {0}")]
    Io(#[from] io::Error),
    #[error("contig '{name}' length {length} exceeds u32::MAX")]
    ContigLengthOverflow { name: String, length: u64 },
    #[error("internal: failed to format current timestamp as RFC3339")]
    TimestampFormat,
}

// ---------------------------------------------------------------------
// Top-level driver
// ---------------------------------------------------------------------

/// Construct the four stages from `args` and drive the pipeline.
/// Records on success are written to `<output>.tmp` and atomically
/// renamed to `<output>` on success. Stderr receives a one-shot
/// run-summary block after the pipeline exhausts.
///
/// The Stage 1 borrow chain (`reader → BAQ → adapter → walker`) is
/// extracted into
/// [`with_stage1_pipeline`](super::stage1_pipeline::with_stage1_pipeline)
/// so the same code path drives both this subcommand and (Task 8)
/// `var-calling-from-bam`.
pub fn run_pileup(args: &PileupArgs) -> Result<(), PileupCliError> {
    let stage1 = &args.stage1;

    // 1. Translate CLI → module configs.
    let alignment_cfg = alignment_config_from_args(stage1);
    let baq_cfg = baq_config_from_args(stage1);
    let walker_cfg = walker_config_from_args(stage1);

    // 2. Size rayon's global pool when --threads is given. The
    //    helper is idempotent (silent no-op on second call) so
    //    library consumers and test runners can invoke multiple
    //    `run_*` helpers back-to-back; the first caller wins. M13
    //    from the 2026-05-19 cohort CLI review.
    configure_rayon_pool(args.threads).map_err(|_| PileupCliError::RayonAlreadyConfigured)?;

    // 3. Tmp output path. The closure opens the file lazily so the
    //    .tmp file only appears once the Stage 1 chain is fully set
    //    up — failures before this point leave no garbage on disk.
    let tmp_path = tmp_path_for(&args.output);

    // 4. Drive the Stage 1 chain via the shared helper. The closure
    //    runs while the chain is alive; counter snapshots happen
    //    afterwards inside the helper.
    let outputs = super::stage1_pipeline::with_stage1_pipeline::<_, PileupCliError, _>(
        &args.alignment_files,
        &args.reference,
        alignment_cfg,
        baq_cfg,
        walker_cfg,
        stage1.baq_chunk_size,
        stage1.no_baq,
        |ctx| {
            // Build the writer header from the reader's metadata (the
            // M5 hard-error fires here).
            let header = build_writer_header(
                ctx.sample_name,
                &args.reference,
                ctx.contigs,
                &args.alignment_files,
                stage1,
                args.block_target_bytes,
            )?;
            // Open the .tmp output, build the writer, drive the walker.
            let file = File::create(&tmp_path).map_err(PileupCliError::Io)?;
            let buf = BufWriter::with_capacity(DEFAULT_BUFFERED_IO_CAPACITY, file);
            let psp_writer = PspWriter::new_with_block_target(buf, header, args.block_target_bytes)
                .map_err(PileupToPspError::from)?;
            let (sink, run_summary) = drive_pileup_to_psp(ctx.walker, psp_writer)?;
            Ok((sink, run_summary))
        },
    )?;

    let (buf_sink, summary) = outputs.result;
    let filter_counts_for_summary = outputs.run_summary.filter_counts;
    let baq_skip_counts_for_summary = outputs.run_summary.baq_skip_counts;

    finalise_output(buf_sink, &tmp_path, &args.output)?;

    // 5. Surface a deferred upstream error if one was stashed by the
    //    bridge. This is the only place upstream errors are
    //    reported — the walker would have exited cleanly from its
    //    perspective.
    if let Some(e) = outputs.stashed_upstream_error {
        let _ = fs::remove_file(&args.output); // best-effort cleanup
        return Err(PileupCliError::AlignmentInput(e));
    }

    // 6. Format the stderr run-summary block.
    print_run_summary(
        &outputs.sample_name,
        &filter_counts_for_summary,
        baq_skip_counts_for_summary.as_ref(),
        &summary,
    );

    Ok(())
}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

pub(crate) fn alignment_config_from_args(
    args: &shared_args::Stage1Args,
) -> AlignmentMergedReaderConfig {
    AlignmentMergedReaderConfig {
        min_mapq: if args.min_mapq == 0 {
            None
        } else {
            Some(args.min_mapq)
        },
        min_read_length: if args.min_read_length == 0 {
            None
        } else {
            Some(args.min_read_length)
        },
        drop_qc_fail: !args.keep_qc_fail,
        drop_duplicate: !args.keep_duplicates,
        max_read_mismatch_fraction: if args.max_read_mismatch_fraction == 0.0 {
            None
        } else {
            Some(args.max_read_mismatch_fraction)
        },
        mismatch_bq_floor: args.mismatch_bq_floor,
    }
}

pub(crate) fn baq_config_from_args(args: &shared_args::Stage1Args) -> BaqConfig {
    BaqConfig {
        gap_open_prob: args.baq_gap_open_prob,
        gap_extend_prob: args.baq_gap_extend_prob,
        band_half_width: args.baq_band_half_width,
    }
}

pub(crate) fn walker_config_from_args(args: &shared_args::Stage1Args) -> WalkerConfig {
    WalkerConfig {
        max_snp_column_depth: args.max_snp_column_depth,
        max_indel_column_depth: args.max_indel_column_depth,
        max_record_span: args.max_record_span,
        mate_lookup_window: args.mate_lookup_window,
        max_active_reads: args.max_active_reads,
    }
}

/// Append `.tmp` to the file name. Works for any path: `foo.psp` →
/// `foo.psp.tmp`, `out` → `out.tmp`.
fn tmp_path_for(output: &Path) -> PathBuf {
    let mut s = output.as_os_str().to_os_string();
    s.push(".tmp");
    PathBuf::from(s)
}

/// Recover the underlying File from the BufWriter the walker→writer
/// seam returned, fsync it, then atomically rename the tmp file
/// into place.
fn finalise_output(
    buf: BufWriter<File>,
    tmp_path: &Path,
    final_path: &Path,
) -> Result<(), PileupCliError> {
    let file = buf.into_inner().map_err(|e| {
        // `IntoInnerError` carries both the writer and the IO error;
        // we just want the io::Error for the user-facing message.
        PileupCliError::Io(e.into_error())
    })?;
    file.sync_all().map_err(PileupCliError::Io)?;
    fs::rename(tmp_path, final_path).map_err(PileupCliError::Io)?;
    Ok(())
}

/// Build a [`WriterHeader`] from the merged reader's metadata + CLI
/// inputs. Hard-errors when a contig has no `@SQ M5`.
fn build_writer_header(
    sample: &str,
    fasta_path: &Path,
    contigs: &ContigList,
    cram_paths: &[PathBuf],
    args: &shared_args::Stage1Args,
    block_target_bytes: usize,
) -> Result<WriterHeader, PileupCliError> {
    let mut chromosomes = Vec::with_capacity(contigs.entries.len());
    for entry in &contigs.entries {
        let length =
            u32::try_from(entry.length).map_err(|_| PileupCliError::ContigLengthOverflow {
                name: entry.name.clone(),
                length: entry.length,
            })?;
        let md5 = entry
            .md5
            .map(format_md5_hex)
            .ok_or_else(|| PileupCliError::MissingMd5 {
                contig: entry.name.clone(),
            })?;
        chromosomes.push(ChromosomeEntry {
            name: entry.name.clone(),
            length,
            md5,
        });
    }

    let created = rfc3339_now()
        .parse()
        .map_err(|_| PileupCliError::TimestampFormat)?;

    let input_crams: Vec<String> = cram_paths.iter().map(|p| basename(p)).collect();
    let input_fasta = basename(fasta_path);
    let parameters = effective_parameters(args, block_target_bytes);

    Ok(WriterHeader {
        format_version: (1, 0),
        sample: sample.to_string(),
        reference: input_fasta.clone(),
        created,
        chromosomes,
        writer: WriterProvenance {
            tool: "pop_var_caller".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            subcommand: "pileup".to_string(),
            input_crams,
            input_fasta,
            parameters,
        },
    })
}

/// Effective parameter map, written into `WriterProvenance.parameters`
/// so the .psp self-describes the run. Every CLI knob — including the
/// ones left at their defaults — is recorded.
fn effective_parameters(
    args: &shared_args::Stage1Args,
    block_target_bytes: usize,
) -> BTreeMap<String, ParameterValue> {
    let mut p = BTreeMap::new();
    // CRAM input filters
    p.insert(
        "min_mapq".into(),
        ParameterValue::Integer(args.min_mapq as i64),
    );
    p.insert(
        "min_read_length".into(),
        ParameterValue::Integer(args.min_read_length as i64),
    );
    p.insert(
        "drop_qc_fail".into(),
        ParameterValue::Boolean(!args.keep_qc_fail),
    );
    p.insert(
        "drop_duplicate".into(),
        ParameterValue::Boolean(!args.keep_duplicates),
    );
    p.insert(
        "max_read_mismatch_fraction".into(),
        ParameterValue::Float(args.max_read_mismatch_fraction as f64),
    );
    p.insert(
        "mismatch_bq_floor".into(),
        ParameterValue::Integer(args.mismatch_bq_floor as i64),
    );
    // BAQ
    p.insert("baq_enabled".into(), ParameterValue::Boolean(!args.no_baq));
    p.insert(
        "baq_gap_open_prob".into(),
        ParameterValue::Float(args.baq_gap_open_prob as f64),
    );
    p.insert(
        "baq_gap_extend_prob".into(),
        ParameterValue::Float(args.baq_gap_extend_prob as f64),
    );
    p.insert(
        "baq_band_half_width".into(),
        ParameterValue::Integer(args.baq_band_half_width as i64),
    );
    p.insert(
        "baq_chunk_size".into(),
        ParameterValue::Integer(args.baq_chunk_size as i64),
    );
    // Walker
    p.insert(
        "max_snp_column_depth".into(),
        ParameterValue::Integer(args.max_snp_column_depth as i64),
    );
    p.insert(
        "max_indel_column_depth".into(),
        ParameterValue::Integer(args.max_indel_column_depth as i64),
    );
    p.insert(
        "max_record_span".into(),
        ParameterValue::Integer(args.max_record_span as i64),
    );
    p.insert(
        "mate_lookup_window".into(),
        ParameterValue::Integer(args.mate_lookup_window as i64),
    );
    p.insert(
        "max_active_reads".into(),
        ParameterValue::Integer(args.max_active_reads as i64),
    );
    // Threads (effective count after rayon was sized).
    let threads = rayon::current_num_threads();
    p.insert("threads".into(), ParameterValue::Integer(threads as i64));
    // PSP writer.
    p.insert(
        "block_target_bytes".into(),
        ParameterValue::Integer(block_target_bytes as i64),
    );
    p
}

fn print_run_summary(
    sample: &str,
    filter_counts: &FilterCounts,
    baq_skip: Option<&BaqSkipCounts>,
    summary: &RunSummary,
) {
    eprintln!("sample: {sample}");
    eprintln!(
        "filter_counts: unmapped={} secondary={} supplementary={} qc_fail={} duplicate={} \
         low_mapq={} too_short={} high_mismatch_fraction={} bad_cigar={} baq_rejected={}",
        filter_counts.unmapped,
        filter_counts.secondary,
        filter_counts.supplementary,
        filter_counts.qc_fail,
        filter_counts.duplicate,
        filter_counts.low_mapq,
        filter_counts.too_short,
        filter_counts.high_mismatch_fraction,
        filter_counts.bad_cigar,
        filter_counts.baq_rejected,
    );
    match baq_skip {
        Some(b) => eprintln!("baq: total_skipped={}", b.total),
        None => eprintln!("baq: disabled"),
    }
    eprintln!(
        "walker: reads_admitted={} records_emitted={} record_widen_events={} \
         mate_overlap_positions={} active_reads_high_water={} mate_lookup_evictions={} \
         column_depth_truncations={}",
        summary.reads_admitted,
        summary.records_emitted,
        summary.record_widen_events,
        summary.mate_overlap_positions,
        summary.active_reads_high_water,
        summary.mate_lookup_evictions,
        summary.column_depth_truncations,
    );
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::{
        DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH,
        DEFAULT_MISMATCH_BQ_FLOOR,
    };
    use crate::baq::{
        SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH, SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
        SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
    };
    use crate::pileup::per_sample::baq_stream::DEFAULT_BAQ_CHUNK_SIZE;
    use crate::pileup::walker::{
        DEFAULT_MATE_LOOKUP_WINDOW, DEFAULT_MAX_ACTIVE_READS, DEFAULT_MAX_INDEL_COLUMN_DEPTH,
        DEFAULT_MAX_RECORD_SPAN, DEFAULT_MAX_SNP_COLUMN_DEPTH,
    };

    fn default_args(
        reference: PathBuf,
        output: PathBuf,
        alignment_files: Vec<PathBuf>,
    ) -> PileupArgs {
        PileupArgs {
            reference,
            output,
            threads: None,
            alignment_files,
            block_target_bytes: TARGET_BLOCK_BYTES,
            stage1: shared_args::Stage1Args {
                min_mapq: DEFAULT_MIN_MAPQ,
                no_baq: false,
                min_read_length: DEFAULT_MIN_READ_LENGTH,
                keep_qc_fail: false,
                keep_duplicates: false,
                max_read_mismatch_fraction: DEFAULT_MAX_READ_MISMATCH_FRACTION,
                mismatch_bq_floor: DEFAULT_MISMATCH_BQ_FLOOR,
                baq_gap_open_prob: SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
                baq_gap_extend_prob: SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
                baq_band_half_width: SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
                baq_chunk_size: DEFAULT_BAQ_CHUNK_SIZE,
                max_snp_column_depth: DEFAULT_MAX_SNP_COLUMN_DEPTH,
                max_indel_column_depth: DEFAULT_MAX_INDEL_COLUMN_DEPTH,
                max_record_span: DEFAULT_MAX_RECORD_SPAN,
                mate_lookup_window: DEFAULT_MATE_LOOKUP_WINDOW,
                max_active_reads: DEFAULT_MAX_ACTIVE_READS,
            },
        }
    }

    #[test]
    fn parse_mismatch_fraction_accepts_in_range() {
        assert_eq!(parse_mismatch_fraction("0.0").unwrap(), 0.0);
        assert_eq!(parse_mismatch_fraction("0.1").unwrap(), 0.1);
        assert_eq!(parse_mismatch_fraction("1.0").unwrap(), 1.0);
    }

    #[test]
    fn parse_mismatch_fraction_rejects_out_of_range() {
        assert!(parse_mismatch_fraction("-0.001").is_err());
        assert!(parse_mismatch_fraction("1.0001").is_err());
        assert!(parse_mismatch_fraction("nan").is_err());
        assert!(parse_mismatch_fraction("inf").is_err());
        assert!(parse_mismatch_fraction("hello").is_err());
    }

    #[test]
    fn tmp_path_appends_suffix() {
        assert_eq!(
            tmp_path_for(Path::new("/tmp/x.psp")),
            PathBuf::from("/tmp/x.psp.tmp")
        );
        assert_eq!(tmp_path_for(Path::new("out")), PathBuf::from("out.tmp"));
    }

    #[test]
    fn effective_parameters_records_every_knob() {
        let args = default_args(PathBuf::from("r.fa"), PathBuf::from("o.psp"), vec![]);
        let p = effective_parameters(&args.stage1, args.block_target_bytes);
        // Sanity: every knob shows up.
        for key in &[
            "min_mapq",
            "min_read_length",
            "drop_qc_fail",
            "drop_duplicate",
            "max_read_mismatch_fraction",
            "mismatch_bq_floor",
            "baq_enabled",
            "baq_gap_open_prob",
            "baq_gap_extend_prob",
            "baq_band_half_width",
            "baq_chunk_size",
            "max_snp_column_depth",
            "max_indel_column_depth",
            "max_record_span",
            "mate_lookup_window",
            "max_active_reads",
            "threads",
            "block_target_bytes",
        ] {
            assert!(p.contains_key(*key), "missing parameter: {key}");
        }
        // Default booleans land where we expect.
        assert_eq!(p.get("drop_qc_fail"), Some(&ParameterValue::Boolean(true)));
        assert_eq!(p.get("baq_enabled"), Some(&ParameterValue::Boolean(true)));
        // Block-target default tracks `TARGET_BLOCK_BYTES`.
        assert_eq!(
            p.get("block_target_bytes"),
            Some(&ParameterValue::Integer(TARGET_BLOCK_BYTES as i64))
        );
    }

    #[test]
    fn build_writer_header_errors_on_missing_md5() {
        use crate::fasta::{ContigEntry, ContigList};
        let contigs = ContigList {
            entries: vec![ContigEntry {
                name: "chr1".into(),
                length: 100,
                md5: None,
            }],
        };
        let args = default_args(PathBuf::from("r.fa"), PathBuf::from("o.psp"), vec![]);
        let err = build_writer_header(
            "s",
            Path::new("r.fa"),
            &contigs,
            &[],
            &args.stage1,
            args.block_target_bytes,
        )
        .expect_err("must error");
        assert!(matches!(err, PileupCliError::MissingMd5 { ref contig } if contig == "chr1"));
    }

    #[test]
    fn build_writer_header_strips_paths_to_basenames() {
        use crate::fasta::{ContigEntry, ContigList};
        let md5 = [0u8; 16];
        let contigs = ContigList {
            entries: vec![ContigEntry {
                name: "chr1".into(),
                length: 100,
                md5: Some(md5),
            }],
        };
        let args = default_args(
            PathBuf::from("/data/ref/grch38.fa"),
            PathBuf::from("o.psp"),
            vec![],
        );
        let alignment_files = [
            PathBuf::from("/scratch/run1/sample.cram"),
            PathBuf::from("./run2/sample2.cram"),
        ];
        let header = build_writer_header(
            "NA12878",
            Path::new("/data/ref/grch38.fa"),
            &contigs,
            &alignment_files,
            &args.stage1,
            args.block_target_bytes,
        )
        .expect("build_writer_header");
        assert_eq!(header.writer.input_fasta, "grch38.fa");
        assert_eq!(
            header.writer.input_crams,
            vec!["sample.cram", "sample2.cram"]
        );
        // MD5 round-trips as 32 hex chars.
        assert_eq!(header.chromosomes[0].md5.len(), 32);
        // Subcommand identifier matches the binary surface.
        assert_eq!(header.writer.subcommand, "pileup");
        assert_eq!(header.writer.tool, "pop_var_caller");
    }
}
