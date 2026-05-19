//! `pop_var_caller var-calling-from-bam` — single-sample BAM → VCF.
//!
//! Convenience subcommand for one-off calls: pileup + var-calling fused
//! in one process, no intermediate `.psp` file. Single sample only —
//! multi-sample work always goes through the `.psp` route (which is
//! both faster and the only path that supports `--contamination-estimates`).
//!
//! Architecture: the Stage 1 borrow chain is built inside
//! [`with_stage1_pipeline`]; the closure receives the walker, wraps it
//! as a k=1 input to [`PerPositionMerger`], and chains the rest of the
//! cohort pipeline exactly as `var-calling` does. The walker's
//! `WalkerError` is stashed by an inline `.scan()`-based error-shedding
//! adapter (backed by `Rc<RefCell<Option<WalkerError>>>` — same shape
//! as [`ErrorSheddingAdapter`](crate::pop_var_caller::cli::error_bridge::ErrorSheddingAdapter)
//! for CRAM-input errors) so the merger sees a clean
//! `Result<_, PspReadError>` stream; any stashed walker error is
//! surfaced after the pipeline drains.
//!
//! Plan: `doc/devel/implementation_plans/pop_var_caller_cohort_cli.md`
//! §"Subcommand: var-calling-from-bam".

use std::cell::RefCell;
use std::ffi::OsString;
use std::io;
use std::path::PathBuf;
use std::rc::Rc;
use std::sync::Arc;

use clap::Args;
use thiserror::Error;

use crate::per_sample_pileup::baq::{
    BaqConfig, DEFAULT_BAQ_CHUNK_SIZE, SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
    SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB, SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
};
use crate::per_sample_pileup::cram_input::{
    ContigList, CramMergedReaderConfig, DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ,
    DEFAULT_MIN_READ_LENGTH, DEFAULT_MISMATCH_BQ_FLOOR,
};
use crate::per_sample_pileup::pileup::{
    DEFAULT_MATE_LOOKUP_WINDOW, DEFAULT_MAX_ACTIVE_READS, DEFAULT_MAX_INDEL_COLUMN_DEPTH,
    DEFAULT_MAX_RECORD_SPAN, DEFAULT_MAX_SNP_COLUMN_DEPTH, PileupRecord, WalkerConfig, WalkerError,
};
use crate::per_sample_pileup::psp::PspReadError;
use crate::per_sample_pileup::psp::header::ParsedChromosome;
use crate::per_sample_pileup::ref_fetcher::SyncRefFetcher;
use crate::pop_var_caller::cli::PileupCliError;
use crate::pop_var_caller::cli::parsers;
use crate::pop_var_caller::stage1_pipeline::{Stage1PipelineContext, with_stage1_pipeline};
use crate::var_calling::dust_filter::{
    DEFAULT_DUST_THRESHOLD, DEFAULT_DUST_WINDOW, DustFilter, DustFilterConfig, DustFilterError,
};
use crate::var_calling::per_group_merger::{
    DEFAULT_BATCH_SIZE, DEFAULT_MAX_ALLELES_PER_RECORD, DEFAULT_PLOIDY, PerGroupMerger,
    PerGroupMergerConfig, PerGroupMergerError, SharedRefFetcher,
};
use crate::var_calling::per_position_merger::{
    PerPositionMerger, PerPositionMergerError, PerPositionPileups,
};
use crate::var_calling::posterior_engine::{
    DEFAULT_COMPOUND_ALT_PSEUDOCOUNT, DEFAULT_CONVERGENCE_THRESHOLD,
    DEFAULT_INBREEDING_COEFFICIENT, DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MAX_GQ_PHRED,
    DEFAULT_MAX_ITERATIONS, DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT, PosteriorEngine,
    PosteriorEngineConfig, PosteriorEngineConfigError, PosteriorEngineError,
};
use crate::var_calling::variant_grouping::{
    DEFAULT_MAX_VARIANT_GROUP_SPAN, GrouperConfig, GrouperConfigError, GrouperError, VariantGrouper,
};
use crate::var_calling::vcf_writer::{
    CohortMetadata, CohortVcfWriter, DEFAULT_EMIT_GP, VcfWriteError, WriterConfig, tmp_path_for,
};

// ---------------------------------------------------------------------
// CLI surface
// ---------------------------------------------------------------------

/// Arguments accepted by `pop_var_caller var-calling-from-bam`. Single
/// sample's CRAM(s) → VCF, no `.psp` on disk. No contamination knob
/// (contamination correction requires the `.psp` route — see plan).
///
/// The fields duplicate `PileupArgs` and `VarCallingArgs` (minus
/// contamination); future cleanup may share them via
/// `#[command(flatten)]` sub-structs.
#[derive(Debug, Args, Clone)]
pub struct VarCallingFromBamArgs {
    // ===== Stage 1 — Common ===================================
    /// Reference FASTA. A sibling `.fai` is required.
    #[arg(long)]
    pub reference: PathBuf,

    /// Output VCF path. `.vcf.gz` / `.vcf.bgz` → bgzf; anything else
    /// → plain text. Atomic `<output>.tmp` → `fs::rename`.
    #[arg(long)]
    pub output: PathBuf,

    /// Worker threads for the BAQ stage and any other rayon work.
    #[arg(long)]
    pub threads: Option<usize>,

    /// One or more coordinate-sorted CRAMs for the sample.
    #[arg(required = true)]
    pub crams: Vec<PathBuf>,

    /// Drop reads with MAPQ < N. `0` admits everything.
    #[arg(long, default_value_t = DEFAULT_MIN_MAPQ)]
    pub min_mapq: u8,

    /// Skip the BAQ HMM. `bq_baq` becomes a copy of the raw CRAM QUAL.
    #[arg(long)]
    pub no_baq: bool,

    /// Skip the low-complexity (sdust) filter entirely.
    #[arg(long)]
    pub no_complexity_filter: bool,

    /// Cohort-wide ploidy.
    #[arg(long, default_value_t = DEFAULT_PLOIDY, value_parser = parsers::parse_ploidy)]
    pub ploidy: u8,

    // ===== Advanced — CRAM input filters ======================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MIN_READ_LENGTH,
        help_heading = "Advanced — CRAM input filters",
    )]
    pub min_read_length: u32,

    #[arg(
        long,
        hide_short_help = true,
        help_heading = "Advanced — CRAM input filters"
    )]
    pub keep_qc_fail: bool,

    #[arg(
        long,
        hide_short_help = true,
        help_heading = "Advanced — CRAM input filters"
    )]
    pub keep_duplicates: bool,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_READ_MISMATCH_FRACTION,
        value_parser = crate::pop_var_caller::cli::parse_mismatch_fraction,
        help_heading = "Advanced — CRAM input filters",
    )]
    pub max_read_mismatch_fraction: f32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MISMATCH_BQ_FLOOR,
        help_heading = "Advanced — CRAM input filters",
    )]
    pub mismatch_bq_floor: u8,

    // ===== Advanced — BAQ HMM =================================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
        help_heading = "Advanced — BAQ HMM",
    )]
    pub baq_gap_open_prob: f32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
        help_heading = "Advanced — BAQ HMM",
    )]
    pub baq_gap_extend_prob: f32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
        help_heading = "Advanced — BAQ HMM",
    )]
    pub baq_band_half_width: i32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_BAQ_CHUNK_SIZE,
        help_heading = "Advanced — BAQ HMM",
    )]
    pub baq_chunk_size: usize,

    // ===== Advanced — Pileup walker ===========================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_SNP_COLUMN_DEPTH,
        help_heading = "Advanced — Pileup walker",
    )]
    pub max_snp_column_depth: u32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_INDEL_COLUMN_DEPTH,
        help_heading = "Advanced — Pileup walker",
    )]
    pub max_indel_column_depth: u32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_RECORD_SPAN,
        help_heading = "Advanced — Pileup walker",
    )]
    pub max_record_span: u32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MATE_LOOKUP_WINDOW,
        help_heading = "Advanced — Pileup walker",
    )]
    pub mate_lookup_window: u32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_ACTIVE_READS,
        help_heading = "Advanced — Pileup walker",
    )]
    pub max_active_reads: u32,

    // ===== Advanced — DUST filter =============================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_DUST_WINDOW,
        value_parser = parsers::parse_dust_window,
        help_heading = "Advanced — DUST filter",
    )]
    pub complexity_window: u32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_DUST_THRESHOLD,
        value_parser = parsers::parse_dust_threshold,
        help_heading = "Advanced — DUST filter",
    )]
    pub complexity_threshold: u32,

    // ===== Advanced — Variant grouping ========================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_VARIANT_GROUP_SPAN,
        value_parser = parsers::parse_var_group_max_span,
        help_heading = "Advanced — Variant grouping",
    )]
    pub var_group_max_span: u32,

    // ===== Advanced — Per-group merger ========================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_ALLELES_PER_RECORD,
        value_parser = parsers::parse_max_alleles,
        help_heading = "Advanced — Per-group merger",
    )]
    pub max_alleles_per_var: usize,

    // ===== Advanced — Posterior engine ========================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_INBREEDING_COEFFICIENT,
        value_parser = parsers::parse_inbreeding_coefficient,
        help_heading = "Advanced — Posterior engine",
    )]
    pub inbreeding_coefficient: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_CONVERGENCE_THRESHOLD,
        value_parser = parsers::parse_em_convergence_threshold,
        help_heading = "Advanced — Posterior engine",
    )]
    pub em_convergence_threshold: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_ITERATIONS,
        value_parser = parsers::parse_em_max_iterations,
        help_heading = "Advanced — Posterior engine",
    )]
    pub em_max_iterations: u32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_REF_PSEUDOCOUNT,
        value_parser = parsers::parse_ref_pseudocount,
        help_heading = "Advanced — Posterior engine",
    )]
    pub ref_pseudocount: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_SNP_ALT_PSEUDOCOUNT,
        value_parser = parsers::parse_snp_alt_pseudocount,
        help_heading = "Advanced — Posterior engine",
    )]
    pub snp_alt_pseudocount: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_INDEL_ALT_PSEUDOCOUNT,
        value_parser = parsers::parse_indel_alt_pseudocount,
        help_heading = "Advanced — Posterior engine",
    )]
    pub indel_alt_pseudocount: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_COMPOUND_ALT_PSEUDOCOUNT,
        value_parser = parsers::parse_compound_alt_pseudocount,
        help_heading = "Advanced — Posterior engine",
    )]
    pub compound_alt_pseudocount: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_GQ_PHRED,
        value_parser = parsers::parse_max_gq_phred,
        help_heading = "Advanced — Posterior engine",
    )]
    pub max_gq_phred: f64,

    // ===== Advanced — VCF writer ==============================
    /// Emit `GP` (genotype posteriors) `FORMAT` per sample. Off by
    /// default — `GP` is `Number=G`, so the per-sample cell grows as
    /// `(ploidy + n_alleles − 1) choose ploidy` (21 floats at
    /// ploidy=2, n_alleles=6). Opt in when posteriors are wanted on
    /// disk.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_EMIT_GP,
        help_heading = "Advanced — VCF writer",
    )]
    pub emit_gp: bool,
}

// ---------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------

/// Errors surfaced by [`run_var_calling_from_bam`].
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum VarCallingFromBamCliError {
    /// Stage 1 setup / IO failure surfaced from `with_stage1_pipeline`.
    #[error("stage 1: {0}")]
    Stage1(#[from] PileupCliError),

    #[error("io: {0}")]
    Io(#[from] io::Error),

    #[error("walker: {0}")]
    Walker(#[from] WalkerError),

    #[error("psp shape adapter: {0}")]
    Psp(#[from] PspReadError),

    #[error("merger: {0}")]
    Merger(#[from] PerPositionMergerError),

    #[error("dust filter: {0}")]
    Dust(#[from] DustFilterError),

    #[error("grouper: {0}")]
    Grouper(#[from] GrouperError),

    #[error("per-group merger: {0}")]
    PerGroup(#[from] PerGroupMergerError),

    #[error("posterior engine: {0}")]
    Posterior(#[from] PosteriorEngineError),

    #[error("vcf writer: {0}")]
    Vcf(#[from] VcfWriteError),

    #[error("grouper config: {0}")]
    GrouperConfig(#[from] GrouperConfigError),

    #[error("per-group merger config: {0}")]
    PerGroupConfig(#[from] crate::var_calling::per_group_merger::PerGroupMergerConfigError),

    #[error("posterior engine config: {0}")]
    PosteriorConfig(#[from] PosteriorEngineConfigError),

    /// One contig in the CRAM header has no `@SQ M5` checksum.
    #[error(
        "contig '{contig}' has no @SQ M5 checksum in the input CRAM(s); \
         re-CRAM with a tool that emits M5 (e.g. samtools)"
    )]
    MissingMd5 { contig: String },

    /// Contig length exceeds `u32::MAX` (would overflow the cohort
    /// pipeline's `ParsedChromosome`).
    #[error("contig '{name}' length {length} exceeds u32::MAX")]
    ContigLengthOverflow { name: String, length: u64 },

    /// `rayon::ThreadPoolBuilder::build_global()` already ran in this
    /// process. The binary calls it at most once.
    #[error("rayon thread pool already initialised — refusing to override")]
    RayonAlreadyConfigured,
}

// ---------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------

/// Run pileup + var-calling fused, no `.psp` intermediate. Single
/// sample (one BAM-set in, one VCF out).
pub fn run_var_calling_from_bam(
    args: &VarCallingFromBamArgs,
) -> Result<(), VarCallingFromBamCliError> {
    // 1. Rayon pool — same once-per-process policy as the other
    //    subcommands.
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|_| VarCallingFromBamCliError::RayonAlreadyConfigured)?;
    }

    // 2. Build the Stage 1 configs from args.
    let cram_cfg = cram_config_from_args(args);
    let baq_cfg = baq_config_from_args(args);
    let walker_cfg = walker_config_from_args(args);

    // 3. Build per-stage configs (cohort side).
    let dust_cfg = DustFilterConfig::new(args.complexity_window, args.complexity_threshold)
        .map_err(VarCallingFromBamCliError::Dust)?;
    let grouper_cfg = GrouperConfig::new(args.var_group_max_span)?;
    let per_group_cfg =
        PerGroupMergerConfig::new(args.ploidy, args.max_alleles_per_var, DEFAULT_BATCH_SIZE)?;
    let mut posterior_cfg = PosteriorEngineConfig::new(
        args.em_convergence_threshold,
        args.em_max_iterations,
        args.ref_pseudocount,
        args.snp_alt_pseudocount,
        args.indel_alt_pseudocount,
        args.compound_alt_pseudocount,
        args.inbreeding_coefficient,
        args.max_gq_phred,
    )?;
    // No contamination — explicit "none" sentinel.
    posterior_cfg.contamination = None;

    // 4. Drive the Stage 1 chain via the shared helper. The named
    //    helper `run_cohort_pipeline_for_single_sample` takes the
    //    walker over and runs the rest of the cohort pipeline; its
    //    signature names every captured value so the data flow is
    //    legible (previously: ~95-line closure capturing six values).
    let no_complexity_filter = args.no_complexity_filter;
    let emit_gp = args.emit_gp;
    let reference = args.reference.clone();
    let output = args.output.clone();
    let command_line = current_command_line();

    let outputs = with_stage1_pipeline::<(), VarCallingFromBamCliError, _>(
        &args.crams,
        &args.reference,
        cram_cfg,
        baq_cfg,
        walker_cfg,
        args.baq_chunk_size,
        args.no_baq,
        |ctx| {
            run_cohort_pipeline_for_single_sample(
                ctx,
                &reference,
                &output,
                command_line,
                no_complexity_filter,
                emit_gp,
                dust_cfg,
                grouper_cfg,
                per_group_cfg,
                posterior_cfg,
            )
        },
    )?;

    // 5. Stage 1 summary — print after the closure returns, mirrors
    //    `run_pileup`'s pattern. Filter / BAQ counts come from the
    //    helper's snapshot.
    let filter_counts = outputs.run_summary.filter_counts;
    let baq_skip = outputs.run_summary.baq_skip_counts;
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

    // 6. Surface a deferred CRAM-input error if one was stashed.
    if let Some(e) = outputs.stashed_upstream_error {
        let _ = std::fs::remove_file(&args.output);
        return Err(VarCallingFromBamCliError::Stage1(
            PileupCliError::CramInput(e),
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

/// Drive Stages 3–6 over a single sample's Stage 1 walker.
///
/// Pulled out of the inline closure that used to live inside
/// [`run_var_calling_from_bam`] so the data flow is legible (every
/// captured value becomes a named parameter) and so an eventual
/// shared `drive_cohort_pipeline` (deferred review item M11) has a
/// clear extraction surface. The closure's outer captures map 1:1
/// onto the parameters here.
///
/// Walker errors are stashed via an inline `.scan()` adapter — the
/// merger sees a clean `Result<_, PspReadError>` stream and never
/// observes the walker's `Err` arm; the stash is read after the
/// pipeline drains and converted into [`VarCallingFromBamCliError::Walker`].
///
/// On any driver-loop failure (write_record / posterior emission /
/// finish), the writer's `<output>.tmp` is best-effort removed —
/// matching `run_pileup`'s discipline (M6 parity).
#[allow(clippy::too_many_arguments)]
fn run_cohort_pipeline_for_single_sample(
    ctx: Stage1PipelineContext<'_>,
    reference: &std::path::Path,
    output: &std::path::Path,
    command_line: String,
    no_complexity_filter: bool,
    emit_gp: bool,
    dust_cfg: DustFilterConfig,
    grouper_cfg: GrouperConfig,
    per_group_cfg: PerGroupMergerConfig,
    posterior_cfg: PosteriorEngineConfig,
) -> Result<(), VarCallingFromBamCliError> {
    // Convert ContigList → Vec<ParsedChromosome>. Single validated
    // source of truth for chromosomes; the merger consumes a clone
    // and the DUST branch consumes the original.
    let chromosomes = contigs_to_parsed(ctx.contigs)?;
    let sample_name_owned = ctx.sample_name.to_string();

    // Cohort metadata + writer config.
    let metadata = CohortMetadata {
        sample_names: vec![sample_name_owned.clone()],
        contigs: chromosomes.clone(),
        tool_string: format!("pop_var_caller {}", env!("CARGO_PKG_VERSION")),
        command_line,
    };
    let writer_cfg = WriterConfig::new(output.to_path_buf()).with_emit_gp(emit_gp);

    // Reference fetcher shared between DUST + PerGroupMerger.
    let fetcher_concrete = SyncRefFetcher::new(reference, ctx.contigs.clone())
        .map_err(VarCallingFromBamCliError::Io)?;
    let fetcher: SharedRefFetcher = Arc::new(fetcher_concrete);

    // Wrap the walker as a k=1 input for PerPositionMerger. The
    // shim stashes walker errors and reports end-of-stream,
    // mirroring `ErrorSheddingAdapter`'s pattern for CRAM input
    // errors.
    //
    // Invariant: the stash is touched only on the merger's main
    // thread. `PerGroupMerger::refill` collects upstream items into
    // a `Vec` before fanning rayon workers across it, so the walker
    // chain itself never crosses thread boundaries. If upstream of
    // `PerGroupMerger` ever becomes multi-threaded, this needs to
    // become `Arc<Mutex<...>>` (and the constraint should move to
    // the bridge module).
    let walker_error_stash: Rc<RefCell<Option<WalkerError>>> = Rc::new(RefCell::new(None));
    let stash_clone = walker_error_stash.clone();
    let walker_iter = ctx.walker.scan(stash_clone, |stash, r| match r {
        Ok(rec) => Some(Ok::<PileupRecord, PspReadError>(rec)),
        Err(e) => {
            *stash.borrow_mut() = Some(e);
            None
        }
    });
    let walker_iter: Box<dyn Iterator<Item = Result<PileupRecord, PspReadError>>> =
        Box::new(walker_iter);

    // Build merger / DUST / grouper / per-group / posterior. The
    // merger consumes a clone of `chromosomes`; the DUST branch
    // consumes the original (both branches converge on a single
    // boxed iterator yielding `Result<_, GrouperError>`).
    let merger = PerPositionMerger::new(
        vec![walker_iter],
        vec![sample_name_owned.clone()],
        chromosomes.clone(),
    )?;
    let upstream_for_grouper: Box<dyn Iterator<Item = Result<PerPositionPileups, GrouperError>>> =
        if no_complexity_filter {
            Box::new(merger.map(|r| r.map_err(GrouperError::from)))
        } else {
            let dust = DustFilter::new(merger, &*fetcher, chromosomes, dust_cfg);
            Box::new(dust.map(|r| r.map_err(GrouperError::from)))
        };
    let grouper = VariantGrouper::with_config(upstream_for_grouper, grouper_cfg);
    let per_group = PerGroupMerger::with_config(grouper, fetcher.clone(), per_group_cfg);
    let posterior = PosteriorEngine::with_config(per_group, posterior_cfg);

    // Open writer, stream records, finalise. On any failure inside
    // the loop, best-effort remove `<output>.tmp` — `finish()`
    // consumes `self`, so a `?` short-circuit leaves the tmp on
    // disk otherwise. Matches `run_var_calling`'s M6 cleanup
    // discipline.
    let tmp_path = tmp_path_for(output);
    let mut writer = CohortVcfWriter::new(metadata, writer_cfg)?;
    let mut records_written: u64 = 0;
    let drive_result: Result<(), VarCallingFromBamCliError> = (|| {
        for item in posterior {
            let record = item?;
            writer.write_record(&record)?;
            records_written += 1;
        }
        writer.finish()?;
        Ok(())
    })();
    if let Err(e) = drive_result {
        let _ = std::fs::remove_file(&tmp_path); // best-effort cleanup
        return Err(e);
    }

    // Surface walker errors stashed by the scan-adapter. The
    // pipeline ran to completion, so `<output>` has been renamed
    // from `<output>.tmp` — retract that publish (publish-then-
    // retract semantics: we'd rather no file than a truncated one).
    if let Some(e) = walker_error_stash.borrow_mut().take() {
        // best-effort cleanup; the rename already happened
        let _ = std::fs::remove_file(output);
        return Err(VarCallingFromBamCliError::Walker(e));
    }

    eprintln!(
        "var-calling-from-bam: sample={} records_emitted={}",
        sample_name_owned, records_written,
    );
    Ok(())
}

fn cram_config_from_args(args: &VarCallingFromBamArgs) -> CramMergedReaderConfig {
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

fn baq_config_from_args(args: &VarCallingFromBamArgs) -> BaqConfig {
    BaqConfig {
        gap_open_prob: args.baq_gap_open_prob,
        gap_extend_prob: args.baq_gap_extend_prob,
        band_half_width: args.baq_band_half_width,
    }
}

fn walker_config_from_args(args: &VarCallingFromBamArgs) -> WalkerConfig {
    WalkerConfig {
        max_snp_column_depth: args.max_snp_column_depth,
        max_indel_column_depth: args.max_indel_column_depth,
        max_record_span: args.max_record_span,
        mate_lookup_window: args.mate_lookup_window,
        max_active_reads: args.max_active_reads,
    }
}

/// Convert the CRAM-side [`ContigList`] into the .psp-shaped
/// [`ParsedChromosome`] vector the cohort pipeline consumes. Hard-errors
/// on missing M5 (the cohort merger requires it) or contig lengths
/// exceeding `u32::MAX`.
fn contigs_to_parsed(
    contigs: &ContigList,
) -> Result<Vec<ParsedChromosome>, VarCallingFromBamCliError> {
    contigs
        .entries
        .iter()
        .map(|e| {
            let length = u32::try_from(e.length).map_err(|_| {
                VarCallingFromBamCliError::ContigLengthOverflow {
                    name: e.name.clone(),
                    length: e.length,
                }
            })?;
            let md5 =
                e.md5
                    .map(format_md5_hex)
                    .ok_or_else(|| VarCallingFromBamCliError::MissingMd5 {
                        contig: e.name.clone(),
                    })?;
            Ok(ParsedChromosome {
                name: e.name.clone(),
                length,
                md5,
            })
        })
        .collect()
}

/// Format a 16-byte MD5 as 32 lowercase hex characters. Mirrors the
/// helper in `cli.rs`; duplicated rather than re-exported to keep the
/// subcommand self-contained.
fn format_md5_hex(bytes: [u8; 16]) -> String {
    use std::fmt::Write as _;
    let mut out = String::with_capacity(32);
    for b in bytes {
        // PANIC-FREE: `std::fmt::Write` on a `String` is infallible
        // — the only failure mode is OOM, which the runtime turns
        // into an abort, not an `Err`.
        write!(&mut out, "{b:02x}").expect("writing to a String never fails");
    }
    out
}

fn current_command_line() -> String {
    std::env::args_os()
        .map(|a: OsString| a.to_string_lossy().into_owned())
        .collect::<Vec<_>>()
        .join(" ")
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::per_sample_pileup::cram_input::ContigEntry;

    #[test]
    fn contigs_to_parsed_strips_md5_hex() {
        let md5 = [0x6a; 16];
        let contigs = ContigList {
            entries: vec![ContigEntry {
                name: "chr1".into(),
                length: 1000,
                md5: Some(md5),
            }],
        };
        let parsed = contigs_to_parsed(&contigs).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].name, "chr1");
        assert_eq!(parsed[0].length, 1000);
        assert_eq!(parsed[0].md5.len(), 32);
    }

    #[test]
    fn contigs_to_parsed_errors_on_missing_md5() {
        let contigs = ContigList {
            entries: vec![ContigEntry {
                name: "chr1".into(),
                length: 1000,
                md5: None,
            }],
        };
        let err = contigs_to_parsed(&contigs).unwrap_err();
        assert!(matches!(
            err,
            VarCallingFromBamCliError::MissingMd5 { ref contig } if contig == "chr1"
        ));
    }

    #[test]
    fn contigs_to_parsed_errors_on_overflow() {
        let contigs = ContigList {
            entries: vec![ContigEntry {
                name: "chr1".into(),
                length: (u32::MAX as u64) + 1,
                md5: Some([0u8; 16]),
            }],
        };
        let err = contigs_to_parsed(&contigs).unwrap_err();
        assert!(matches!(
            err,
            VarCallingFromBamCliError::ContigLengthOverflow { .. }
        ));
    }
}
