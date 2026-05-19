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
use thiserror::Error;

use crate::per_sample_pileup::cram_input::ContigList;
use crate::per_sample_pileup::pileup::{PileupRecord, WalkerError};
use crate::per_sample_pileup::psp::PspReadError;
use crate::per_sample_pileup::psp::header::ParsedChromosome;
use crate::per_sample_pileup::ref_fetcher::SyncRefFetcher;
use crate::pop_var_caller::cli::PileupCliError;
use crate::pop_var_caller::cli::error_bridge::ErrorSheddingAdapter;
use crate::pop_var_caller::cohort_driver::{CohortPipelineParams, drive_cohort_pipeline};
use crate::pop_var_caller::common::{current_command_line, format_md5_hex};
use crate::pop_var_caller::stage1_pipeline::{Stage1PipelineContext, with_stage1_pipeline};
use crate::var_calling::dust_filter::{DustFilterConfig, DustFilterError};
use crate::var_calling::per_group_merger::{
    DEFAULT_BATCH_SIZE, PerGroupMergerConfig, PerGroupMergerError, SharedRefFetcher,
};
use crate::var_calling::per_position_merger::{PerPositionMerger, PerPositionMergerError};
use crate::var_calling::posterior_engine::{
    PosteriorEngineConfig, PosteriorEngineConfigError, PosteriorEngineError,
};
use crate::var_calling::variant_grouping::{GrouperConfig, GrouperConfigError, GrouperError};
use crate::var_calling::vcf_writer::{CohortMetadata, VcfWriteError, WriterConfig};

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

    // 1. Rayon pool — same once-per-process policy as the other
    //    subcommands.
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|_| VarCallingFromBamCliError::RayonAlreadyConfigured)?;
    }

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
    dust_cfg: DustFilterConfig,
    grouper_cfg: GrouperConfig,
    per_group_cfg: PerGroupMergerConfig,
    posterior_cfg: PosteriorEngineConfig,
) -> Result<(), VarCallingFromBamCliError> {
    // Convert ContigList → Vec<ParsedChromosome>. Single validated
    // source of truth for chromosomes; the merger consumes a clone
    // and the cohort driver takes the original through
    // CohortPipelineParams.
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
    // adapter stashes walker errors and reports end-of-stream; the
    // `.map(Ok)` lifts the bare `PileupRecord`s back into the
    // `Result<_, PspReadError>` shape that `PerPositionMerger`
    // wants (no walker error will ever flow through this Result —
    // it is stashed in the handle instead).
    //
    // Invariant: the stash is touched only on the merger's main
    // thread. `PerGroupMerger::refill` collects upstream items into
    // a `Vec` before fanning rayon workers across it, so the walker
    // chain itself never crosses thread boundaries. If upstream of
    // `PerGroupMerger` ever becomes multi-threaded, the stash
    // becomes `Arc<Mutex<...>>` (and the constraint should move to
    // `error_bridge`).
    let walker_adapter: ErrorSheddingAdapter<_, PileupRecord, WalkerError> =
        ErrorSheddingAdapter::new(ctx.walker);
    let walker_error_handle = walker_adapter.error_handle();
    let walker_iter: Box<dyn Iterator<Item = Result<PileupRecord, PspReadError>>> =
        Box::new(walker_adapter.map(Ok));

    // Build the merger from the walker-shim, then hand the rest of
    // the pipeline to the shared cohort driver.
    let merger = PerPositionMerger::new(
        vec![walker_iter],
        vec![sample_name_owned.clone()],
        chromosomes.clone(),
    )?;
    let pipeline_params = CohortPipelineParams {
        no_complexity_filter,
        dust_cfg,
        grouper_cfg,
        per_group_cfg,
        posterior_cfg,
        fetcher,
        chromosomes,
    };
    let records_written = drive_cohort_pipeline::<_, VarCallingFromBamCliError>(
        merger,
        pipeline_params,
        output,
        metadata,
        writer_cfg,
    )?;

    // Surface walker errors stashed by the ErrorSheddingAdapter.
    // The pipeline ran to completion, so `<output>` has been
    // renamed from `<output>.tmp` — retract that publish
    // (publish-then-retract semantics: we'd rather no file than a
    // truncated one).
    if let Some(e) = walker_error_handle.take() {
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
