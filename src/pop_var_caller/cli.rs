//! `pop_var_caller pileup` — the Stage 1 CLI orchestrator.
//!
//! Glues `CramMergedReader` → `BaqStream` (or BAQ-bypass passthrough)
//! → pileup walker → `.psp` writer. Argument parsing is delegated to
//! clap; everything else is in [`run_pileup`].
//!
//! Plan: `ia/feature_implementation_plans/pop_var_caller_pileup_cli.md`.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{Args, Parser, Subcommand};
use thiserror::Error;

use crate::per_sample_caller::baq::{
    BaqConfig, BaqSkipCounts, BaqStream, DEFAULT_BAQ_CHUNK_SIZE, SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
    SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB, SAMTOOLS_ILLUMINA_GAP_OPEN_PROB, prepare_passthrough,
};
use crate::per_sample_caller::cram_input::{
    ContigList, CramMergedReader, CramMergedReaderConfig, DEFAULT_MAX_READ_MISMATCH_FRACTION,
    DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH, DEFAULT_MISMATCH_BQ_FLOOR, FilterCounts,
};
use crate::per_sample_caller::errors::CramInputError;
use crate::per_sample_caller::pileup;
use crate::per_sample_caller::pileup::{
    DEFAULT_MATE_LOOKUP_WINDOW, DEFAULT_MAX_ACTIVE_SLOTS, DEFAULT_MAX_INDEL_COLUMN_DEPTH,
    DEFAULT_MAX_RECORD_SPAN, DEFAULT_MAX_SNP_COLUMN_DEPTH, RunSummary, WalkerConfig,
};
use crate::per_sample_caller::pileup_to_psp::{PileupToPspError, drive_pileup_to_psp};
use crate::per_sample_caller::psp::header::{
    ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance,
};
use crate::per_sample_caller::psp::writer::PspWriter;
use crate::per_sample_caller::ref_fetcher::{ChromBoundaryRefFetcher, SyncRefFetcher};

pub mod error_bridge;
use error_bridge::ErrorSheddingAdapter;

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
}

/// Arguments accepted by the `pileup` subcommand. The struct is the
/// authoritative knob list; the orchestrator translates it into module
/// configs. Defaults are pulled from each module's `DEFAULT_*` consts
/// so they stay in sync.
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

    /// One or more coordinate-sorted CRAMs for the sample.
    #[arg(required = true)]
    pub crams: Vec<PathBuf>,

    /// Drop reads with MAPQ < N. `0` admits everything.
    #[arg(long, default_value_t = DEFAULT_MIN_MAPQ)]
    pub min_mapq: u8,

    /// Skip the BAQ HMM. `bq_baq` becomes a copy of the raw CRAM
    /// QUAL.
    #[arg(long)]
    pub no_baq: bool,

    // ===== Advanced flags (hidden from `-h`, shown in `--help`) =====
    /// Drop reads with decoded SEQ length < N. `0` admits any length.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MIN_READ_LENGTH,
        help_heading = "Advanced — CRAM input filters",
    )]
    pub min_read_length: u32,

    /// Keep reads with the QC-fail flag (0x200) set.
    #[arg(
        long,
        hide_short_help = true,
        help_heading = "Advanced — CRAM input filters"
    )]
    pub keep_qc_fail: bool,

    /// Keep reads with the duplicate flag (0x400) set.
    #[arg(
        long,
        hide_short_help = true,
        help_heading = "Advanced — CRAM input filters"
    )]
    pub keep_duplicates: bool,

    /// Drop reads whose M-op mismatch fraction exceeds X. Must be in
    /// [0.0, 1.0]; pass 0.0 to disable the filter entirely. Values
    /// outside the range are a hard error.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_READ_MISMATCH_FRACTION,
        value_parser = parse_mismatch_fraction,
        help_heading = "Advanced — CRAM input filters",
    )]
    pub max_read_mismatch_fraction: f32,

    /// BQ floor below which a mismatch does not count toward
    /// --max-read-mismatch-fraction. `0` makes every mismatch count.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MISMATCH_BQ_FLOOR,
        help_heading = "Advanced — CRAM input filters",
    )]
    pub mismatch_bq_floor: u8,

    /// BAQ gap-open probability. samtools/htslib default: 1e-3.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
        help_heading = "Advanced — BAQ HMM",
    )]
    pub baq_gap_open_prob: f32,

    /// BAQ gap-extension probability. samtools/htslib default: 0.1.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
        help_heading = "Advanced — BAQ HMM",
    )]
    pub baq_gap_extend_prob: f32,

    /// BAQ band half-width. samtools/htslib default: 7. The engine
    /// auto-widens per-read on long indels.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
        help_heading = "Advanced — BAQ HMM",
    )]
    pub baq_band_half_width: i32,

    /// BAQ batch size — reads processed in parallel per rayon chunk.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_BAQ_CHUNK_SIZE,
        help_heading = "Advanced — BAQ HMM",
    )]
    pub baq_chunk_size: usize,

    /// Max contributors folded at a pure-SNP/REF column.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_SNP_COLUMN_DEPTH,
        help_heading = "Advanced — Pileup walker",
    )]
    pub max_snp_column_depth: u32,

    /// Max contributors folded at a column carrying any indel
    /// observation.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_INDEL_COLUMN_DEPTH,
        help_heading = "Advanced — Pileup walker",
    )]
    pub max_indel_column_depth: u32,

    /// Hard cap on per-record reference span.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_RECORD_SPAN,
        help_heading = "Advanced — Pileup walker",
    )]
    pub max_record_span: u32,

    /// How far past a first mate the walker keeps a pending-mates
    /// entry before evicting and treating the first mate as solo.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MATE_LOOKUP_WINDOW,
        help_heading = "Advanced — Pileup walker",
    )]
    pub mate_lookup_window: u32,

    /// Hard cap on concurrently-active reads (defensive bound;
    /// exceeding it is a hard error).
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_ACTIVE_SLOTS,
        help_heading = "Advanced — Pileup walker",
    )]
    pub max_active_reads: u32,
}

fn parse_mismatch_fraction(s: &str) -> Result<f32, String> {
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
pub enum PileupCliError {
    #[error("CRAM input: {0}")]
    CramInput(#[from] CramInputError),
    #[error("pipeline: {0}")]
    Pipeline(#[from] PileupToPspError),
    #[error(
        "contig '{contig}' has no @SQ M5 checksum in the input CRAM(s); \
         re-CRAM with a tool that emits M5 (e.g. samtools)"
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
pub fn run_pileup(args: &PileupArgs) -> Result<(), PileupCliError> {
    // 1. Translate CLI → module configs.
    let cram_cfg = cram_config_from_args(args);
    let baq_cfg = baq_config_from_args(args);
    let walker_cfg = walker_config_from_args(args);

    // 2. Size rayon's global pool when --threads is given. If the
    //    pool is already configured (e.g. a prior call in the same
    //    process), we report it; the binary only calls this once.
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|_| PileupCliError::RayonAlreadyConfigured)?;
    }

    // 3. Open the merged reader, capture metadata for the writer header.
    let mut reader = CramMergedReader::new(&args.crams, &args.reference, cram_cfg)?;
    let sample_name = reader.sample_name().to_string();
    let contigs = reader.contigs().clone();

    // 4. Build writer header (this is where missing @SQ M5 hard-errors).
    let header = build_writer_header(&sample_name, &args.reference, &contigs, &args.crams, args)?;

    // 5. Reference fetchers — two of them. BAQ runs in parallel and
    //    needs a Sync fetcher; the walker is single-threaded and uses
    //    the chrom-boundary-evicting one. See [`SyncRefFetcher`]
    //    docs for the memory trade-off.
    let baq_fetcher =
        SyncRefFetcher::new(&args.reference, contigs.clone()).map_err(PileupCliError::Io)?;
    let walker_fetcher = ChromBoundaryRefFetcher::new(&args.reference, contigs.clone())
        .map_err(PileupCliError::Io)?;

    // 6. Open the output sink at `<output>.tmp`. Atomic rename on
    //    success — partial files never appear at `output`.
    let tmp_path = tmp_path_for(&args.output);
    let file = File::create(&tmp_path).map_err(PileupCliError::Io)?;
    let buf = BufWriter::with_capacity(64 * 1024, file);
    let psp_writer = PspWriter::new(buf, header).map_err(PileupToPspError::from)?;

    // 7. Assemble the pipeline. `reader` and `baq_stream` (when on)
    //    are kept by-ref through the chain so they survive until we
    //    pull `FilterCounts` / `BaqSkipCounts` for the run-summary
    //    block at the end. Borrow order is:
    //
    //        &mut reader → BaqStream<&mut reader> → &mut baq_stream
    //                  → ErrorSheddingAdapter      → &mut adapter
    //                  → PileupWalker              → drive_pileup_to_psp
    //
    //    Each adapter is dropped in reverse before the one below it
    //    is touched again, so the borrow checker accepts the chain.
    //    Both branches end with `buf_sink` ready for `finalise_output`,
    //    plus the three counter snapshots and any stashed upstream error.
    let summary: RunSummary;
    let stashed_error: Option<CramInputError>;
    let baq_skip_counts_for_summary: Option<BaqSkipCounts>;
    let buf_sink: BufWriter<File>;

    if args.no_baq {
        let prepared = reader.by_ref().map(|r| {
            r.map(|read| {
                let chrom_id = u32::try_from(read.ref_id).expect("ref_id fits u32");
                prepare_passthrough(read, chrom_id)
            })
        });
        let mut adapter = ErrorSheddingAdapter::new(prepared);
        let handle = adapter.error_handle();
        let walker = pileup::run(adapter.by_ref(), &walker_fetcher, &walker_cfg);
        let (sink, run_summary) = drive_pileup_to_psp(walker, psp_writer)?;
        summary = run_summary;
        stashed_error = handle.take();
        baq_skip_counts_for_summary = None;
        buf_sink = sink;
        // Drop adapter (and the map closure inside) so the &mut
        // borrow of `reader` is released — `filter_counts()` next.
        drop(adapter);
    } else {
        let mut baq_stream =
            BaqStream::new(reader.by_ref(), baq_cfg, &baq_fetcher, args.baq_chunk_size);
        let mut adapter = ErrorSheddingAdapter::new(baq_stream.by_ref());
        let handle = adapter.error_handle();
        let walker = pileup::run(adapter.by_ref(), &walker_fetcher, &walker_cfg);
        let (sink, run_summary) = drive_pileup_to_psp(walker, psp_writer)?;
        summary = run_summary;
        stashed_error = handle.take();
        buf_sink = sink;
        // Adapter borrows baq_stream; drop it first.
        drop(adapter);
        // Now baq_stream is accessible — snapshot skip counts before dropping.
        baq_skip_counts_for_summary = Some(*baq_stream.skip_counts());
        drop(baq_stream); // releases &mut reader for `filter_counts()` below
    }

    // Pull FilterCounts before the writer-tmp is renamed — the
    // counts live in the reader and we want them in the summary even
    // if the rename later fails.
    let filter_counts_for_summary: FilterCounts = *reader.filter_counts();

    finalise_output(buf_sink, &tmp_path, &args.output)?;

    // 8. Surface a deferred upstream error if one was stashed by the
    //    bridge. This is the only place upstream errors are
    //    reported — the walker would have exited cleanly from its
    //    perspective.
    if let Some(e) = stashed_error {
        let _ = fs::remove_file(&args.output); // best-effort cleanup
        return Err(PileupCliError::CramInput(e));
    }

    // 9. Format the stderr run-summary block.
    print_run_summary(
        &sample_name,
        &filter_counts_for_summary,
        baq_skip_counts_for_summary.as_ref(),
        &summary,
    );

    Ok(())
}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

fn cram_config_from_args(args: &PileupArgs) -> CramMergedReaderConfig {
    CramMergedReaderConfig {
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

fn baq_config_from_args(args: &PileupArgs) -> BaqConfig {
    BaqConfig {
        gap_open_prob: args.baq_gap_open_prob,
        gap_extend_prob: args.baq_gap_extend_prob,
        band_half_width: args.baq_band_half_width,
    }
}

fn walker_config_from_args(args: &PileupArgs) -> WalkerConfig {
    WalkerConfig {
        max_snp_column_depth: args.max_snp_column_depth,
        max_indel_column_depth: args.max_indel_column_depth,
        max_record_span: args.max_record_span,
        mate_lookup_window: args.mate_lookup_window,
        max_active_slots: args.max_active_reads,
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
    args: &PileupArgs,
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
    let parameters = effective_parameters(args);

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

fn basename(p: &Path) -> String {
    p.file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| p.to_string_lossy().into_owned())
}

/// Format a 16-byte MD5 as 32 lowercase hex characters, matching
/// SAM `@SQ M5` and the `.psp` `chromosome.md5` field.
fn format_md5_hex(bytes: [u8; 16]) -> String {
    let mut out = String::with_capacity(32);
    for b in bytes {
        use std::fmt::Write as _;
        write!(&mut out, "{b:02x}").expect("writing to a String never fails");
    }
    out
}

/// Effective parameter map, written into `WriterProvenance.parameters`
/// so the .psp self-describes the run. Every CLI knob — including the
/// ones left at their defaults — is recorded.
fn effective_parameters(args: &PileupArgs) -> BTreeMap<String, ParameterValue> {
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
         mate_overlap_positions={} slot_high_water={} mate_lookup_evictions={} \
         column_depth_truncations={}",
        summary.reads_admitted,
        summary.records_emitted,
        summary.record_widen_events,
        summary.mate_overlap_positions,
        summary.slot_high_water,
        summary.mate_lookup_evictions,
        summary.column_depth_truncations,
    );
}

// ---------------------------------------------------------------------
// RFC3339 timestamp (no date-crate dep)
// ---------------------------------------------------------------------

/// Format the current UTC time as a TOML-compatible RFC3339 string
/// (`YYYY-MM-DDTHH:MM:SSZ`). Done by hand so we don't pull in a date
/// crate just for this one call. Uses Howard Hinnant's civil-date
/// algorithm to turn days-since-epoch into (year, month, day).
fn rfc3339_now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = (secs / 86_400) as i64;
    let sod = secs % 86_400;
    let (y, m, d) = civil_from_days(days);
    let h = sod / 3600;
    let min = (sod % 3600) / 60;
    let s = sod % 60;
    format!("{y:04}-{m:02}-{d:02}T{h:02}:{min:02}:{s:02}Z")
}

/// Howard Hinnant's `civil_from_days` — translate days-since-epoch
/// (1970-01-01 = 0) into a proleptic Gregorian (year, month, day).
/// Verified by external reference; see
/// http://howardhinnant.github.io/date_algorithms.html
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64; // [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365; // [0, 399]
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32; // [1, 31]
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32; // [1, 12]
    let y = if m <= 2 { y + 1 } else { y };
    (y, m, d)
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn default_args(reference: PathBuf, output: PathBuf, crams: Vec<PathBuf>) -> PileupArgs {
        PileupArgs {
            reference,
            output,
            threads: None,
            crams,
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
            max_active_reads: DEFAULT_MAX_ACTIVE_SLOTS,
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
    fn format_md5_hex_is_32_lowercase() {
        let bytes = [
            0x6a, 0xef, 0x89, 0x7c, 0x3d, 0x6f, 0xf0, 0xc7, 0x8a, 0xff, 0x06, 0xac, 0x18, 0x91,
            0x78, 0xdd,
        ];
        assert_eq!(format_md5_hex(bytes), "6aef897c3d6ff0c78aff06ac189178dd");
    }

    #[test]
    fn basename_strips_directory() {
        assert_eq!(basename(Path::new("/tmp/a/b/sample.cram")), "sample.cram");
        assert_eq!(basename(Path::new("ref.fa")), "ref.fa");
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
    fn civil_from_days_matches_known_dates() {
        // Unix epoch
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        // 1970-12-31 (day 364, 1970 is not a leap year)
        assert_eq!(civil_from_days(364), (1970, 12, 31));
        // 2000-01-01 (day 10957 since epoch, common Y2K test)
        assert_eq!(civil_from_days(10957), (2000, 1, 1));
        // 2024-02-29 — leap day
        assert_eq!(civil_from_days(19782), (2024, 2, 29));
    }

    #[test]
    fn rfc3339_now_parses_as_toml_datetime() {
        let s = rfc3339_now();
        let _: toml::value::Datetime = s.parse().expect("must parse as toml Datetime");
    }

    #[test]
    fn effective_parameters_records_every_knob() {
        let args = default_args(PathBuf::from("r.fa"), PathBuf::from("o.psp"), vec![]);
        let p = effective_parameters(&args);
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
        ] {
            assert!(p.contains_key(*key), "missing parameter: {key}");
        }
        // Default booleans land where we expect.
        assert_eq!(p.get("drop_qc_fail"), Some(&ParameterValue::Boolean(true)));
        assert_eq!(p.get("baq_enabled"), Some(&ParameterValue::Boolean(true)));
    }

    #[test]
    fn build_writer_header_errors_on_missing_md5() {
        use crate::per_sample_caller::cram_input::{ContigEntry, ContigList};
        let contigs = ContigList {
            entries: vec![ContigEntry {
                name: "chr1".into(),
                length: 100,
                md5: None,
            }],
        };
        let args = default_args(PathBuf::from("r.fa"), PathBuf::from("o.psp"), vec![]);
        let err = build_writer_header("s", Path::new("r.fa"), &contigs, &[], &args)
            .expect_err("must error");
        assert!(matches!(err, PileupCliError::MissingMd5 { ref contig } if contig == "chr1"));
    }

    #[test]
    fn build_writer_header_strips_paths_to_basenames() {
        use crate::per_sample_caller::cram_input::{ContigEntry, ContigList};
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
        let crams = [
            PathBuf::from("/scratch/run1/sample.cram"),
            PathBuf::from("./run2/sample2.cram"),
        ];
        let header = build_writer_header(
            "NA12878",
            Path::new("/data/ref/grch38.fa"),
            &contigs,
            &crams,
            &args,
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
