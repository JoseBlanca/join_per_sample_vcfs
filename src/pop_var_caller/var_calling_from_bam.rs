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
//! `WalkerError` is stashed by an
//! [`ErrorSheddingAdapter`](crate::pop_var_caller::cli::error_bridge::ErrorSheddingAdapter)
//! parameterised over `(PileupRecord, WalkerError)` — the same
//! adapter used at the Stage 1 CRAM-input seam, just bound to
//! different `T` and `E` — so the merger sees a clean
//! `Result<_, PspReadError>` stream; any stashed walker error is
//! surfaced after the pipeline drains.
//!
//! Plan: `doc/devel/implementation_plans/pop_var_caller_cohort_cli.md`
//! §"Subcommand: var-calling-from-bam".

use std::io;
use std::path::PathBuf;
use std::sync::Arc;

use clap::Args;
use tempfile::TempDir;
use thiserror::Error;

use crate::fasta::{ContigList, StreamingChromRefFetcher};
use crate::pileup::walker::WalkerError;
use crate::pileup_record::PileupRecord;
use crate::pop_var_caller::cli::PileupCliError;
use crate::pop_var_caller::cli::error_bridge::ErrorSheddingAdapter;
use crate::pop_var_caller::cohort_driver::{
    CohortDriveStats, CohortPipelineParams, drive_cohort_pipeline,
};
use crate::pop_var_caller::common::{configure_rayon_pool, current_command_line, format_md5_hex};
use crate::pop_var_caller::stage1_pipeline::{Stage1PipelineContext, with_stage1_pipeline};
use crate::psp::PspReadError;
use crate::psp::header::ParsedChromosome;
use crate::var_calling::dust_filter::{DustFilterConfig, DustFilterError};
use crate::var_calling::per_group_merger::{
    DEFAULT_BATCH_SIZE, PerGroupMergerConfig, PerGroupMergerError, SharedRefFetcher,
};
use crate::var_calling::per_position_merger::{PerPositionMerger, PerPositionMergerError};
use crate::var_calling::posterior_engine::{
    PosteriorEngineConfig, PosteriorEngineConfigError, PosteriorEngineError,
};
use crate::var_calling::variant_grouping::{GrouperConfig, GrouperConfigError, GrouperError};
use crate::vcf::concat::concat_fragments;
use crate::vcf::{CohortMetadata, VcfWriteError, WriterConfig, path_is_bgzf};

// ---------------------------------------------------------------------
// CLI surface
// ---------------------------------------------------------------------

/// Arguments accepted by `pop_var_caller var-calling-from-bam`. Single
/// sample's CRAM(s) → VCF, no `.psp` on disk. No contamination knob
/// (contamination correction requires the `.psp` route — see plan).
///
/// Flattens both [`Stage1Args`](crate::pop_var_caller::cli::shared_args::Stage1Args)
/// (CRAM-input filters, BAQ HMM, pileup walker — shared with the
/// `pileup` subcommand) and
/// [`CohortPipelineArgs`](crate::pop_var_caller::cli::shared_args::CohortPipelineArgs)
/// (DUST, grouper, per-group merger, posterior engine, ploidy, VCF
/// writer — shared with `var-calling`). M10 from the 2026-05-19
/// cohort CLI review unified the three previously-duplicated args
/// surfaces.
#[derive(Debug, Args, Clone)]
pub struct VarCallingFromBamArgs {
    // ===== Top-level (visible in `-h`) ========================
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

    /// Skip the low-complexity (sdust) filter entirely.
    #[arg(long)]
    pub no_complexity_filter: bool,

    /// If any input CRAM/BAM lacks its alignment index (`.crai` for
    /// CRAM; `.bai` or `.csi` for BAM), build it in place next to
    /// the source file before running. Without this flag, missing
    /// indexes are a hard error.
    ///
    /// Default off — per-chromosome parallelism requires an
    /// alignment index for each input, but the caller must opt in
    /// before we write files next to user inputs.
    #[arg(long)]
    pub build_map_file_index: bool,

    // ===== Stage 1 (shared with pileup) =======================
    #[command(flatten)]
    pub stage1: crate::pop_var_caller::cli::shared_args::Stage1Args,

    // ===== Cohort pipeline (shared with var-calling) ==========
    #[command(flatten)]
    pub cohort: crate::pop_var_caller::cli::shared_args::CohortPipelineArgs,
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

    /// An input alignment file has no index and `--build-map-file-index`
    /// was not passed. The message names both the flag and the
    /// `samtools index` recipe so the user can fix the run without
    /// digging.
    #[error(
        "input alignment file '{path}' has no index (looked for '{expected_index_path}')\n\
         per-chromosome parallelism requires an index for each input.\n\
         either:\n\
           - re-run with --build-map-file-index to have pop_var_caller build it, or\n\
           - run `samtools index {path}` yourself before invoking this command"
    )]
    MissingMapFileIndex {
        path: PathBuf,
        expected_index_path: PathBuf,
    },

    /// `--build-map-file-index` was set but the build failed.
    /// Commonly: the directory holding the source alignment file is
    /// read-only.
    #[error("failed to build alignment index for '{path}': {source}")]
    IndexBuildFailed {
        path: PathBuf,
        #[source]
        source: io::Error,
    },

    /// Input's extension is neither `.cram` nor `.bam`.
    #[error(
        "input alignment file '{path}' has an unsupported extension \
         (expected .cram or .bam)"
    )]
    UnsupportedAlignmentExtension { path: PathBuf },
}

impl From<crate::bam::errors::AlignmentIndexError> for VarCallingFromBamCliError {
    fn from(e: crate::bam::errors::AlignmentIndexError) -> Self {
        use crate::bam::errors::AlignmentIndexError;
        match e {
            AlignmentIndexError::MissingAlignmentIndex {
                path,
                expected_index_path,
            } => Self::MissingMapFileIndex {
                path,
                expected_index_path,
            },
            AlignmentIndexError::BuildFailed { path, source } => {
                Self::IndexBuildFailed { path, source }
            }
            AlignmentIndexError::UnsupportedExtension { path } => {
                Self::UnsupportedAlignmentExtension { path }
            }
        }
    }
}

// ---------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------

/// Run pileup + var-calling fused, no `.psp` intermediate. Single
/// sample (one BAM-set in, one VCF out).
pub fn run_var_calling_from_bam(
    args: &VarCallingFromBamArgs,
) -> Result<(), VarCallingFromBamCliError> {
    let stage1 = &args.stage1;
    let cohort = &args.cohort;

    // 0. Alignment-index pre-flight. Per-chromosome parallelism (the
    //    structural reason this driver exists) requires an index next
    //    to every CRAM/BAM input so each rayon worker can issue a
    //    contig-scoped `Reader::query(...)`. Runs before any other
    //    setup so a missing-index error fires before tempdirs or
    //    rayon pools are touched. Indexes are only created when the
    //    user opts in via `--build-map-file-index`.
    crate::bam::index_preflight::preflight_alignment_indexes(
        &args.crams,
        args.build_map_file_index,
    )?;

    // 1. Rayon pool — idempotent under the silent-no-op policy
    //    (M13, locked 2026-05-19); first caller wins, subsequent
    //    callers' --threads is silently ignored.
    configure_rayon_pool(args.threads)
        .map_err(|_| VarCallingFromBamCliError::RayonAlreadyConfigured)?;

    // 2. Build the Stage 1 configs from args. The helpers live in
    //    `cli.rs` and take `&Stage1Args` directly — same code path
    //    drives this subcommand and `pileup`, no duplication.
    let cram_cfg = crate::pop_var_caller::cli::cram_config_from_args(stage1);
    let baq_cfg = crate::pop_var_caller::cli::baq_config_from_args(stage1);
    let walker_cfg = crate::pop_var_caller::cli::walker_config_from_args(stage1);

    // 3. Build per-stage configs (cohort side).
    let dust_cfg = DustFilterConfig::new(cohort.complexity_window, cohort.complexity_threshold)
        .map_err(VarCallingFromBamCliError::Dust)?;
    let grouper_cfg = GrouperConfig::new(cohort.var_group_max_span)?;
    let per_group_cfg = PerGroupMergerConfig::new(
        cohort.ploidy,
        cohort.max_alleles_per_var,
        DEFAULT_BATCH_SIZE,
    )?;
    // Build the posterior-engine config via the validating
    // builder chain. The `from-bam` subcommand has no
    // `--contamination-estimates` surface, so `contamination`
    // stays at the project default (`None`); the explicit setter
    // call is kept as a "no contamination here" anchor for
    // future maintainers reading this flow alongside
    // `var_calling.rs`.
    let posterior_cfg = PosteriorEngineConfig::new()
        .with_convergence_threshold(cohort.em_convergence_threshold)?
        .with_max_iterations(cohort.em_max_iterations)?
        .with_ref_pseudocount(cohort.ref_pseudocount)?
        .with_snp_alt_pseudocount(cohort.snp_alt_pseudocount)?
        .with_indel_alt_pseudocount(cohort.indel_alt_pseudocount)?
        .with_compound_alt_pseudocount(cohort.compound_alt_pseudocount)?
        .with_fixation_index_default(cohort.inbreeding_coefficient)?
        .with_max_gq_phred(cohort.max_gq_phred)?
        .with_contamination(None)?;

    // 4. Drive the Stage 1 chain via the shared helper. The named
    //    helper `run_cohort_pipeline_for_single_sample` takes the
    //    walker over and runs the rest of the cohort pipeline; its
    //    signature names every captured value so the data flow is
    //    legible (previously: ~95-line closure capturing six values).
    let no_complexity_filter = args.no_complexity_filter;
    let emit_gp = cohort.emit_gp;
    let min_qual_phred = cohort.min_qual_phred;
    let min_alt_obs_per_sample = cohort.min_alt_obs_per_sample;
    let no_mapq_diff_filter = cohort.no_mapq_diff_filter;
    let min_mapq_diff_t = cohort.min_mapq_diff_t;
    let reference = args.reference.clone();
    let output = args.output.clone();
    let command_line = current_command_line();

    let outputs = with_stage1_pipeline::<(), VarCallingFromBamCliError, _>(
        &args.crams,
        &args.reference,
        cram_cfg,
        baq_cfg,
        walker_cfg,
        stage1.baq_chunk_size,
        stage1.no_baq,
        |ctx| {
            run_cohort_pipeline_for_single_sample(
                ctx,
                &reference,
                &output,
                command_line,
                no_complexity_filter,
                emit_gp,
                min_qual_phred,
                min_alt_obs_per_sample,
                no_mapq_diff_filter,
                min_mapq_diff_t,
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
/// matching `run_pileup`'s discipline (M6 parity). The shared
/// [`drive_cohort_pipeline`] handles the merger-onwards wiring and
/// the cleanup loop; this wrapper sets up the walker-shim and
/// metadata, then surfaces any stashed walker error after the
/// driver returns.
#[allow(clippy::too_many_arguments)]
fn run_cohort_pipeline_for_single_sample(
    ctx: Stage1PipelineContext<'_>,
    reference: &std::path::Path,
    output: &std::path::Path,
    command_line: String,
    no_complexity_filter: bool,
    emit_gp: bool,
    min_qual_phred: f64,
    min_alt_obs_per_sample: u32,
    no_mapq_diff_filter: bool,
    min_mapq_diff_t: f32,
    dust_cfg: DustFilterConfig,
    grouper_cfg: GrouperConfig,
    per_group_cfg: PerGroupMergerConfig,
    posterior_cfg: PosteriorEngineConfig,
) -> Result<(), VarCallingFromBamCliError> {
    // Convert ContigList → Vec<ParsedChromosome>. Single validated
    // source of truth for chromosomes; the per-chrom loop clones it
    // for each fragment.
    let chromosomes = contigs_to_parsed(ctx.contigs)?;
    let sample_name_owned = ctx.sample_name.to_string();

    // Cohort metadata + writer config template.
    let metadata = CohortMetadata {
        sample_names: vec![sample_name_owned.clone()],
        contigs: chromosomes.clone(),
        tool_string: format!("pop_var_caller {}", env!("CARGO_PKG_VERSION")),
        command_line,
    };
    let writer_cfg_template = WriterConfig::new(output.to_path_buf()).with_emit_gp(emit_gp);

    // Per-chromosome fragment paths under a TempDir placed next to the
    // final output — matches the cohort `var-calling` path so the
    // concat's atomic rename stays on one filesystem. The TempDir
    // RAII-drops on Ok or on panic; the fragments live just long
    // enough for `concat_fragments` to assemble them into `output`.
    let output_parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| std::path::Path::new("."));
    let frags_dir = TempDir::new_in(output_parent).map_err(VarCallingFromBamCliError::Io)?;
    let frag_ext = if path_is_bgzf(output) {
        "vcf.gz"
    } else {
        "vcf"
    };
    let fragment_paths: Vec<PathBuf> = chromosomes
        .iter()
        .enumerate()
        .map(|(cid, c)| {
            frags_dir
                .path()
                .join(format!("chr_{cid:03}_{name}.{frag_ext}", name = c.name))
        })
        .collect();

    // Wrap the walker as a clean `Result<_, PspReadError>` stream.
    // Walker errors are stashed via `ErrorSheddingAdapter` (the
    // merger sees a clean stream); the stash is surfaced after the
    // per-chrom loop drains.
    //
    // Invariant: the stash is touched only on the merger's main
    // thread. `PerGroupMerger::refill` collects upstream items into
    // a `Vec` before fanning rayon workers across it, so the walker
    // chain itself never crosses thread boundaries.
    let walker_adapter: ErrorSheddingAdapter<_, PileupRecord, WalkerError> =
        ErrorSheddingAdapter::new(ctx.walker);
    let walker_error_handle = walker_adapter.error_handle();
    let walker_iter: Box<dyn Iterator<Item = Result<PileupRecord, PspReadError>>> =
        Box::new(walker_adapter.map(Ok));

    // Drive the cohort pipeline per chrom. The walker emits records
    // in coord order across all chroms; `PerChromRecordsIter` chunks
    // that single stream by `chrom_id`. Each chrom builds its own
    // `StreamingChromRefFetcher::for_contig` (1 MB sliding buffer)
    // that drops at the end of the iteration — no contig is ever
    // wholly resident.
    let pipeline_params = CohortPipelineParams {
        no_complexity_filter,
        dust_cfg,
        grouper_cfg,
        per_group_cfg,
        posterior_cfg,
        min_qual_phred,
        min_alt_obs_per_sample,
        no_mapq_diff_filter,
        min_mapq_diff_t,
    };

    let mut per_chrom = PerChromRecordsIter::new(walker_iter);
    let mut next_walker_chrom = per_chrom.open_next_chrom();
    let mut total_stats = CohortDriveStats::default();

    for (cid, chrom_entry) in chromosomes.iter().enumerate() {
        let chrom_id = cid as u32;
        let walker_is_on_this_chrom = next_walker_chrom == Some(chrom_id);

        // Per-chrom reference fetcher: opens its own indexed FASTA
        // reader, serves bases through a 1 MB sliding buffer, drops
        // at end of scope.
        let streaming = StreamingChromRefFetcher::for_contig(reference, &chrom_entry.name)
            .map_err(|e| {
                VarCallingFromBamCliError::Io(io::Error::other(format!(
                    "ref fetcher construction failed: {e}"
                )))
            })?;
        // See cohort_driver::process_one_chromosome for the
        // `Arc<!Sync>` rationale; switching to `Rc` is a tracked
        // follow-up to the H1 perf fix.
        #[allow(clippy::arc_with_non_send_sync)]
        let fetcher: SharedRefFetcher = Arc::new(streaming);

        // Build the per-chrom record iter. When the walker has no
        // records on this chrom, feed an empty iter — the cohort
        // pipeline still writes a header-only fragment, matching the
        // multi-sample `var-calling` flow.
        let chrom_records: Box<dyn Iterator<Item = Result<PileupRecord, PspReadError>> + '_> =
            if walker_is_on_this_chrom {
                Box::new(per_chrom.consume_current_chrom())
            } else {
                Box::new(std::iter::empty())
            };
        let merger = PerPositionMerger::new(
            vec![chrom_records],
            vec![sample_name_owned.clone()],
            chromosomes.clone(),
        )?;

        let fragment_path = fragment_paths[cid].clone();
        let writer_cfg = WriterConfig {
            output: fragment_path.clone(),
            ..writer_cfg_template.clone()
        };

        let stats = drive_cohort_pipeline::<_, VarCallingFromBamCliError>(
            chrom_id,
            merger,
            pipeline_params.clone(),
            fetcher,
            &fragment_path,
            metadata.clone(),
            writer_cfg,
        )?;
        total_stats.records_written += stats.records_written;
        total_stats.records_unconverged += stats.records_unconverged;
        total_stats.records_dropped_hom_ref += stats.records_dropped_hom_ref;
        total_stats.records_dropped_low_qual += stats.records_dropped_low_qual;
        total_stats.records_dropped_low_alt_obs += stats.records_dropped_low_alt_obs;
        total_stats.records_dropped_low_mapq_diff_t += stats.records_dropped_low_mapq_diff_t;

        if walker_is_on_this_chrom {
            next_walker_chrom = per_chrom.open_next_chrom();
        }
    }

    // Concat fragments in contig-table order into the final output.
    // `concat_fragments` owns the atomic rename + parent directory
    // fsync; on Ok the file is durable at `output`.
    concat_fragments(output, &fragment_paths)?;
    // frags_dir RAII-drops the scratch TempDir + fragments here.

    // Surface walker errors stashed by the ErrorSheddingAdapter.
    // The pipeline ran to completion, so `<output>` has been
    // published — retract that publish (publish-then-retract
    // semantics: we'd rather no file than a truncated one).
    if let Some(e) = walker_error_handle.take() {
        let _ = std::fs::remove_file(output);
        return Err(VarCallingFromBamCliError::Walker(e));
    }

    let emnoconv_note = if total_stats.records_unconverged > 0 {
        format!(
            " records_emnoconv={} (FILTER=EMNoConv; EM iteration cap)",
            total_stats.records_unconverged,
        )
    } else {
        String::new()
    };
    eprintln!(
        "var-calling-from-bam: sample={} records_emitted={}{}",
        sample_name_owned, total_stats.records_written, emnoconv_note,
    );
    Ok(())
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

// ---------------------------------------------------------------------
// Per-chromosome worker (commit 3 of the
// var-calling-from-bam-per-chromosome plan; commit 4 wires it into
// `run_var_calling_from_bam` via `rayon::par_iter`).
// ---------------------------------------------------------------------

/// Drive Stages 1 + 3-6 for **one chromosome** of one sample.
///
/// Mirror of the cohort-side
/// [`crate::pop_var_caller::cohort_driver::process_one_chromosome`]
/// (PSP inputs) but the front of the chain is Stage 1 (CRAM →
/// BAQ → walker) instead of a `PspReader`. Each call owns its
/// reader set, BAQ stream, walker, reference fetcher, and writer
/// outright — no shared mutable state with sibling workers.
///
/// Returns `(chrom_id, CohortDriveStats)` on success. Errors are
/// surfaced through the same typed enum the serial driver uses;
/// the per-chrom orchestrator (commit 4) collects worker results
/// and returns the first error after every worker has joined.
///
/// Stage 1 borrow chain: the returned reader from
/// [`CramMergedReader::query`] is owned on this function's stack;
/// the BAQ stream (or no-BAQ passthrough) borrows from it; the
/// error-shedding adapter borrows from BAQ; the pileup walker
/// borrows from the adapter; the walker is wrapped in a second
/// error-shedding adapter so its `WalkerError` is stashed before
/// the merger sees the stream; the resulting
/// `Result<PileupRecord, PspReadError>` stream is fed into a
/// single-sample [`PerPositionMerger`] which is finally handed to
/// [`drive_cohort_pipeline`].
///
/// `headers` and `indexes` are per-input handles loaded once at
/// driver startup (commit 4) and shared across all workers via
/// `Arc`. `all_chromosomes` is the canonical contig table —
/// `chrom_id` indexes into it. `pipeline_params` is the Stages 3-6
/// config bundle shared across workers (cohort-side
/// `CohortPipelineParams` is `Clone`).
///
/// **No-Stage-1-summary**: per-chrom workers cannot meaningfully
/// snapshot the global `FilterCounts` / `BaqSkipCounts` because
/// the orchestrator runs N of them in parallel. The driver
/// (commit 4) will either drop the per-Stage-1 counter summary
/// for the parallel path, or aggregate counts across workers via
/// a separate accumulator.
// `dead_code` is intentional in commit 3 — this helper lands
// standalone so its diff is reviewable in isolation; commit 4 of
// the var-calling-from-bam-per-chromosome plan wires it into
// `run_var_calling_from_bam` via `rayon::par_iter`, at which point
// the attribute is removed.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)] // bundling would obscure the borrow chain
fn process_one_chromosome_from_bam(
    chrom_id: u32,
    crams: &[PathBuf],
    headers: &[Arc<noodles_sam::Header>],
    indexes: &[crate::bam::index_preflight::AlignmentIndex],
    canonical_contigs: ContigList,
    sample_name: String,
    all_chromosomes: Vec<ParsedChromosome>,
    reference: &std::path::Path,
    cram_cfg: crate::bam::cram_input::CramMergedReaderConfig,
    baq_cfg: crate::baq::BaqConfig,
    walker_cfg: crate::pileup::walker::WalkerConfig,
    baq_chunk_size: usize,
    no_baq: bool,
    fragment_path: PathBuf,
    metadata: crate::vcf::CohortMetadata,
    writer_cfg_template: crate::vcf::WriterConfig,
    pipeline_params: CohortPipelineParams,
) -> Result<(u32, CohortDriveStats), VarCallingFromBamCliError> {
    use crate::pileup::per_sample::baq_engine::prepare_passthrough;
    use crate::pileup::per_sample::baq_stream::BaqStream;
    use crate::pileup::walker::{self, PreparedRead};

    let chrom_entry = all_chromosomes.get(chrom_id as usize).expect(
        "chrom_id past contig table — caller (run_var_calling_from_bam) must validate the bound",
    );

    // 1. Per-worker CRAM record source, bounded to this chrom.
    let mut reader = crate::bam::cram_input::CramMergedReader::query(
        crams,
        reference,
        canonical_contigs.clone(),
        sample_name.clone(),
        headers,
        indexes,
        &chrom_entry.name,
        cram_cfg,
    )
    .map_err(PileupCliError::CramInput)?;

    // 2. Walker fetcher — multi-chrom by trait but only ever sees
    //    this one chrom's records, so the swap-on-chrom-transition
    //    path is exercised at most once.
    let walker_fetcher = crate::fasta::MultiChromStreamingRefFetcher::new(
        reference.to_path_buf(),
        canonical_contigs,
    );

    // 3. Build the BAQ → adapter → walker chain. Mirror of
    //    `with_stage1_pipeline`'s body, restricted to the
    //    per-chrom record stream.
    let (per_chrom_records, cram_error_handle): (
        Box<dyn Iterator<Item = Result<crate::pileup_record::PileupRecord, PspReadError>>>,
        _,
    ) = {
        // (a) Build the chain on this function's stack.
        // (b) Consume it by collecting into a `Vec` — the merger
        //     consumes the records by value anyway. This sidesteps
        //     the self-referential borrow that `with_stage1_pipeline`
        //     resolves via a callback, at the cost of materialising
        //     one chrom's worth of `PileupRecord`s in memory before
        //     handing them to the merger.
        //
        // Per-chrom record counts on real workloads are in the
        // 1-100 K range (sites that survive the filter cascade);
        // each record carries a small `Vec<AlleleObservation>`. The
        // peak memory cost is bounded and rebuilds at the next
        // worker — no cross-worker accumulation.
        let baq_skip = if no_baq {
            // No-BAQ branch: passthrough off the merged reader.
            let mut adapter = ErrorSheddingAdapter::new(reader.by_ref().map(|r| {
                r.map(|read| {
                    let chrom = u32::try_from(read.ref_id).expect("ref_id fits u32");
                    prepare_passthrough(read, chrom)
                })
            }));
            let cram_error_handle = adapter.error_handle();
            let input: Box<dyn Iterator<Item = PreparedRead> + '_> = Box::new(adapter.by_ref());
            let walker = walker::run(input, &walker_fetcher, &walker_cfg);

            let mut walker_adapter: ErrorSheddingAdapter<_, _, WalkerError> =
                ErrorSheddingAdapter::new(walker);
            let walker_error_handle = walker_adapter.error_handle();
            let records: Vec<crate::pileup_record::PileupRecord> =
                walker_adapter.by_ref().collect();
            drop(walker_adapter);
            drop(adapter);

            // If the walker errored mid-stream, surface that
            // immediately — the merger should never see a partial
            // chrom.
            if let Some(e) = walker_error_handle.take() {
                return Err(VarCallingFromBamCliError::Walker(e));
            }

            (records, cram_error_handle)
        } else {
            // BAQ-on: reader → BaqStream → adapter → walker.
            let mut baq_stream = BaqStream::new(
                reader.by_ref(),
                baq_cfg,
                reference.to_path_buf(),
                all_chromosomes_to_contig_list(&all_chromosomes),
                baq_chunk_size,
            );
            let mut adapter = ErrorSheddingAdapter::new(baq_stream.by_ref());
            let cram_error_handle = adapter.error_handle();
            let input: Box<dyn Iterator<Item = PreparedRead> + '_> = Box::new(adapter.by_ref());
            let walker = walker::run(input, &walker_fetcher, &walker_cfg);

            let mut walker_adapter: ErrorSheddingAdapter<_, _, WalkerError> =
                ErrorSheddingAdapter::new(walker);
            let walker_error_handle = walker_adapter.error_handle();
            let records: Vec<crate::pileup_record::PileupRecord> =
                walker_adapter.by_ref().collect();
            drop(walker_adapter);
            drop(adapter);
            drop(baq_stream);

            if let Some(e) = walker_error_handle.take() {
                return Err(VarCallingFromBamCliError::Walker(e));
            }

            (records, cram_error_handle)
        };
        let cram_handle = baq_skip.1;
        let iter: Box<dyn Iterator<Item = Result<_, PspReadError>>> =
            Box::new(baq_skip.0.into_iter().map(Ok));
        (iter, cram_handle)
    };

    // 4. Surface a stashed CRAM-input error from inside the chain
    //    before going downstream.
    if let Some(e) = cram_error_handle.take() {
        return Err(VarCallingFromBamCliError::Stage1(
            PileupCliError::CramInput(e),
        ));
    }

    // 5. Single-sample merger over this chrom's records.
    let merger = crate::var_calling::per_position_merger::PerPositionMerger::new(
        vec![per_chrom_records],
        vec![sample_name.clone()],
        all_chromosomes.clone(),
    )?;

    // 6. Per-worker cohort-side fetcher (DUST + PerGroupMerger).
    //    Mirror of `cohort_driver::process_one_chromosome`'s setup;
    //    constructor returns the project's `StreamingChromRefFetcher`
    //    bound to this chrom.
    let streaming =
        crate::fasta::StreamingChromRefFetcher::for_contig(reference, &chrom_entry.name).map_err(
            |e| {
                VarCallingFromBamCliError::Io(io::Error::other(format!(
                    "ref fetcher construction failed: {e}"
                )))
            },
        )?;
    #[allow(clippy::arc_with_non_send_sync)]
    let fetcher: SharedRefFetcher = Arc::new(streaming);

    // 7. Writer config inherits `emit_gp` and any future flags from
    //    the template; only the output path is per-fragment.
    let writer_cfg = crate::vcf::WriterConfig {
        output: fragment_path.clone(),
        ..writer_cfg_template
    };

    let stats = crate::pop_var_caller::cohort_driver::drive_cohort_pipeline::<
        _,
        VarCallingFromBamCliError,
    >(
        chrom_id,
        merger,
        pipeline_params,
        fetcher,
        &fragment_path,
        metadata,
        writer_cfg,
    )?;

    Ok((chrom_id, stats))
}

/// `BaqStream::new` wants a `ContigList`, not the cohort-side
/// `Vec<ParsedChromosome>`. This helper round-trips through the
/// project's [`ContigList`] shape without losing names or lengths.
///
/// Dead-code-allowed for the same commit-3-only reason as
/// [`process_one_chromosome_from_bam`].
#[allow(dead_code)]
fn all_chromosomes_to_contig_list(chromosomes: &[ParsedChromosome]) -> ContigList {
    ContigList {
        entries: chromosomes
            .iter()
            .map(|c| crate::fasta::ContigEntry {
                name: c.name.clone(),
                length: c.length as u64,
                md5: None,
            })
            .collect(),
    }
}

// ---------------------------------------------------------------------
// PerChromRecordsIter — chunk a coord-sorted record stream by chrom.
// ---------------------------------------------------------------------
//
// The Stage 1 walker emits `PileupRecord`s in coordinate order
// across the whole CRAM input (all chroms in sequence). To drive
// `drive_cohort_pipeline` per chrom — which the
// `cohort_chrom_ref_fetcher_migration` plan requires so the
// pipeline can use the single-contig `ChromRefFetcher` API — we
// need to split that one stream into per-chrom sub-streams.
//
// This helper does it without buffering: it peeks one record
// ahead from the upstream, and `consume_current_chrom()` yields a
// borrowed iterator that stops the moment the next record's
// `chrom_id` differs. The caller then loops, calling
// `open_next_chrom()` to advance to the next chrom and again
// `consume_current_chrom()` for that chrom's records, until
// `open_next_chrom()` returns `None`.
//
// Errors from the upstream propagate through the consumed iter
// — the upstream's `Err` is yielded as-is so the caller sees the
// same `Result<PileupRecord, PspReadError>` shape `PerPositionMerger`
// already expects.

/// Per-chrom view over a coord-sorted record stream. See the
/// module comment above for the lifecycle.
pub struct PerChromRecordsIter<I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    upstream: std::iter::Peekable<I>,
    /// Chrom currently being consumed by `consume_current_chrom`.
    /// `None` before the first `open_next_chrom()` call and after
    /// the upstream is exhausted.
    current_chrom: Option<u32>,
}

impl<I> PerChromRecordsIter<I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    pub fn new(upstream: I) -> Self {
        Self {
            upstream: upstream.peekable(),
            current_chrom: None,
        }
    }

    /// Peek the next upstream record; if there is one, set the
    /// "current chrom" to its `chrom_id` and return `Some(chrom_id)`.
    /// Returns `None` when the upstream is exhausted (clean EOF).
    ///
    /// If the peeked item is an `Err`, the error is **left in the
    /// peek slot** so the next `consume_current_chrom().next()`
    /// returns it to the caller. This keeps error surfacing
    /// uniform with how the underlying iterator behaves —
    /// `open_next_chrom` itself never reports errors directly.
    pub fn open_next_chrom(&mut self) -> Option<u32> {
        match self.upstream.peek() {
            Some(Ok(rec)) => {
                self.current_chrom = Some(rec.chrom_id);
                self.current_chrom
            }
            Some(Err(_)) => {
                // Leave the error in place for `consume_current_chrom`
                // to yield. We don't know the chrom_id of the next
                // valid record, so pretend we're "on" the most
                // recently seen chrom — the consume iter will yield
                // the error and then None, the caller loops back to
                // `open_next_chrom` which sees the iterator drained
                // and returns None.
                self.current_chrom
            }
            None => {
                self.current_chrom = None;
                None
            }
        }
    }

    /// Yield records while their `chrom_id` matches `current_chrom`.
    /// Returns `None` at the chrom transition, at EOF, or after an
    /// upstream error has been surfaced. The caller's loop should
    /// drain this iterator, then call `open_next_chrom()` to advance.
    pub fn consume_current_chrom(&mut self) -> ConsumeCurrentChromIter<'_, I> {
        ConsumeCurrentChromIter { parent: self }
    }
}

/// Borrowed iterator returned by
/// [`PerChromRecordsIter::consume_current_chrom`]. Yields records
/// for the current chrom only.
pub struct ConsumeCurrentChromIter<'a, I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    parent: &'a mut PerChromRecordsIter<I>,
}

impl<'a, I> Iterator for ConsumeCurrentChromIter<'a, I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    type Item = Result<PileupRecord, PspReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        let current = self.parent.current_chrom?;
        match self.parent.upstream.peek() {
            Some(Ok(rec)) if rec.chrom_id == current => {
                // Consume the peeked record and yield it.
                self.parent.upstream.next()
            }
            Some(Ok(_)) => {
                // Different chrom — stop here so the caller loops.
                None
            }
            Some(Err(_)) => {
                // Surface the error to the caller; consuming it here
                // also clears the peek slot so the next
                // `open_next_chrom` sees the upstream cleanly.
                self.parent.upstream.next()
            }
            None => None,
        }
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::ContigEntry;

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

    // ---------------------------------------------------------------
    // PerChromRecordsIter tests
    // ---------------------------------------------------------------

    use crate::pileup_record::PileupRecord;

    fn record_on(chrom_id: u32, pos: u32) -> PileupRecord {
        // Use the public constructor and patch chrom_id/pos.
        // PileupRecord::new takes (chrom_id, pos, alleles), so this
        // is just a thin convenience.
        PileupRecord::new(chrom_id, pos, Vec::new())
    }

    #[test]
    fn per_chrom_iter_yields_chunked_records() {
        // Three chroms, sizes 3 / 2 / 1 records. Two outer loops
        // confirm `open_next_chrom` walks chroms in order and
        // `consume_current_chrom` yields exactly the records on
        // each.
        let inputs: Vec<Result<PileupRecord, PspReadError>> = vec![
            Ok(record_on(0, 1)),
            Ok(record_on(0, 2)),
            Ok(record_on(0, 3)),
            Ok(record_on(1, 1)),
            Ok(record_on(1, 2)),
            Ok(record_on(2, 1)),
        ];
        let mut iter = PerChromRecordsIter::new(inputs.into_iter());

        let mut seen: Vec<(u32, Vec<u32>)> = Vec::new();
        while let Some(chrom_id) = iter.open_next_chrom() {
            let positions: Vec<u32> = iter
                .consume_current_chrom()
                .map(|r| r.unwrap().pos)
                .collect();
            seen.push((chrom_id, positions));
        }
        assert_eq!(
            seen,
            vec![(0, vec![1, 2, 3]), (1, vec![1, 2]), (2, vec![1])],
        );
    }

    #[test]
    fn per_chrom_iter_handles_single_chrom() {
        let inputs: Vec<Result<PileupRecord, PspReadError>> =
            vec![Ok(record_on(7, 100)), Ok(record_on(7, 200))];
        let mut iter = PerChromRecordsIter::new(inputs.into_iter());
        assert_eq!(iter.open_next_chrom(), Some(7));
        let positions: Vec<u32> = iter
            .consume_current_chrom()
            .map(|r| r.unwrap().pos)
            .collect();
        assert_eq!(positions, vec![100, 200]);
        assert_eq!(iter.open_next_chrom(), None);
    }

    #[test]
    fn per_chrom_iter_handles_empty_input() {
        let inputs: Vec<Result<PileupRecord, PspReadError>> = vec![];
        let mut iter = PerChromRecordsIter::new(inputs.into_iter());
        assert_eq!(iter.open_next_chrom(), None);
    }

    #[test]
    fn per_chrom_iter_propagates_upstream_error_mid_stream() {
        // An error mid-stream surfaces through the consumed iter.
        // The caller looks at the Result; here we just confirm the
        // error variant escapes the chunker.
        let inputs: Vec<Result<PileupRecord, PspReadError>> = vec![
            Ok(record_on(0, 1)),
            Err(PspReadError::Io {
                context: "test",
                source: std::io::Error::other("test"),
            }),
        ];
        let mut iter = PerChromRecordsIter::new(inputs.into_iter());
        assert_eq!(iter.open_next_chrom(), Some(0));
        let mut consumed = iter.consume_current_chrom();
        // First the OK record.
        assert_eq!(consumed.next().unwrap().unwrap().pos, 1);
        // Then the error.
        let err = consumed.next().unwrap();
        assert!(err.is_err());
        // Then EOF on the consume iter.
        assert!(consumed.next().is_none());
        // Re-borrow happens implicitly when `consumed` goes out of
        // scope; no need for an explicit `drop` call (ConsumeCurrentChromIter
        // is not Drop-impl'd, just a borrow).
        let _ = consumed;
        // And `open_next_chrom` returns None after the error path.
        assert_eq!(iter.open_next_chrom(), None);
    }

    #[test]
    fn per_chrom_iter_propagates_upstream_error_at_chrom_transition() {
        // Error happens right at a chrom transition. The first
        // chrom's records come through, then the error is yielded
        // in the next consume cycle.
        let inputs: Vec<Result<PileupRecord, PspReadError>> = vec![
            Ok(record_on(0, 1)),
            Err(PspReadError::Io {
                context: "test",
                source: std::io::Error::other("test"),
            }),
            Ok(record_on(1, 1)),
        ];
        let mut iter = PerChromRecordsIter::new(inputs.into_iter());
        assert_eq!(iter.open_next_chrom(), Some(0));
        let chrom0: Vec<_> = iter.consume_current_chrom().collect();
        // chr0 had one OK record; the error halts the chunk
        // boundary at the chrom-0 group.
        assert_eq!(chrom0.len(), 2); // OK + Err in the same group
        assert!(chrom0[0].is_ok());
        assert!(chrom0[1].is_err());
        // After draining, the next open advances to chrom 1.
        assert_eq!(iter.open_next_chrom(), Some(1));
        let chrom1_positions: Vec<u32> = iter
            .consume_current_chrom()
            .map(|r| r.unwrap().pos)
            .collect();
        assert_eq!(chrom1_positions, vec![1]);
    }
}
