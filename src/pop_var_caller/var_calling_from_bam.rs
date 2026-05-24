//! `pop_var_caller var-calling-from-bam` — single-sample BAM → VCF.
//!
//! Single-sample, single-VCF convenience subcommand: pileup + cohort
//! var-calling fused in one process, no intermediate `.psp` file.
//! Multi-sample work always goes through the `.psp` route (which is
//! both faster and the only path that supports
//! `--contamination-estimates`).
//!
//! Architecture: the driver pre-flights every input's alignment
//! index (and optionally builds one on demand via
//! `--build-map-file-index`), harvests the canonical contig list +
//! sample name once via [`crate::bam::alignment_input::AlignmentMergedReader::new`],
//! loads each input's `Arc<sam::Header>` and
//! `Arc<crai::Index>`, and dispatches one [`process_one_chromosome_from_bam`]
//! worker per chromosome via `rayon::par_iter`. Each worker owns
//! its own CRAM reader (via [`crate::bam::alignment_input::AlignmentMergedReader::query`]),
//! BAQ chunk pool, walker, reference fetchers, and per-fragment VCF
//! writer — no shared mutable state across workers. Fragments
//! concat in contig-table order via [`crate::vcf::concat::concat_fragments`],
//! which owns the atomic rename + parent-dir fsync.
//!
//! Plan: `doc/devel/implementation_plans/var_calling_from_bam_per_chromosome.md`.

use std::io;
use std::path::PathBuf;
use std::sync::Arc;

use clap::Args;
use tempfile::TempDir;
use thiserror::Error;

use crate::fasta::ContigList;
use crate::pileup::walker::WalkerError;
use crate::pop_var_caller::cli::PileupCliError;
use crate::pop_var_caller::cli::error_bridge::ErrorSheddingAdapter;
use crate::pop_var_caller::cohort_driver::{CohortDriveStats, CohortPipelineParams};
use crate::pop_var_caller::common::{configure_rayon_pool, current_command_line, format_md5_hex};
use crate::psp::PspReadError;
use crate::psp::header::ParsedChromosome;
use crate::var_calling::dust_filter::{DustFilterConfig, DustFilterError};
use crate::var_calling::per_group_merger::{
    DEFAULT_BATCH_SIZE, PerGroupMergerConfig, PerGroupMergerError, SharedRefFetcher,
};
use crate::var_calling::per_position_merger::PerPositionMergerError;
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

    /// One or more coordinate-sorted CRAM or BAM file(s) for one
    /// sample. All inputs in one invocation must share the same
    /// format (CRAM-only or BAM-only); mixed CRAM + BAM is rejected
    /// with a typed error.
    #[arg(required = true)]
    pub alignment_files: Vec<PathBuf>,

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
    /// Stage 1 setup / IO failure: alignment-input validation,
    /// reader open, header parse.
    #[error("stage 1: {0}")]
    Stage1(#[from] PileupCliError),

    /// Failed to create the per-run tempdir for the per-chrom
    /// VCF fragments. Commonly the output's parent directory is
    /// read-only or missing.
    #[error("failed to create scratch directory under '{parent}': {source}")]
    ScratchDir {
        parent: PathBuf,
        #[source]
        source: io::Error,
    },

    /// Failed to build the per-chrom reference fetcher. Commonly
    /// the FASTA's `.fai` is missing or unreadable, or the contig
    /// name is not present in the indexed FASTA.
    #[error("failed to build reference fetcher for contig '{contig}': {source}")]
    RefFetcher {
        contig: String,
        #[source]
        source: io::Error,
    },
    // NOTE: the previous `Io(#[from] io::Error)` catch-all variant
    // was removed by M2. New `io::Error` origins should be wrapped
    // in operation-named variants (like `ScratchDir` / `RefFetcher`
    // above) at the `?` site, not collapsed back into one bag.
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

    /// One contig in the alignment-file header has no `@SQ M5` checksum.
    #[error(
        "contig '{contig}' has no @SQ M5 checksum in the input alignment file(s); \
         re-create the file with a tool that emits @SQ M5 (e.g. samtools view -t)"
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

    /// Loading a previously-built alignment index from disk
    /// failed (corrupted index, permission denied, …). Surfaced
    /// for both CRAM (`.crai`) and BAM (`.csi`/`.bai`) inputs.
    #[error("failed to load alignment index '{index_path}' for '{path}': {source}")]
    IndexLoadFailed {
        path: PathBuf,
        index_path: PathBuf,
        #[source]
        source: io::Error,
    },

    /// Input's extension is neither `.cram` nor `.bam`.
    #[error(
        "input alignment file '{path}' has an unsupported extension \
         (expected .cram or .bam)"
    )]
    UnsupportedAlignmentExtension { path: PathBuf },

    /// Inputs to a single invocation mixed `.cram` and `.bam` files.
    /// One format per invocation only; the merge downstream consumes
    /// one sample's reads as one ordered stream.
    #[error(
        "mixed alignment-file formats are not supported in one \
         invocation: '{first_path}' is {first_format}, '{other_path}' is {other_format}"
    )]
    MixedAlignmentFormats {
        first_path: PathBuf,
        first_format: &'static str,
        other_path: PathBuf,
        other_format: &'static str,
    },
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
            AlignmentIndexError::LoadFailed {
                path,
                index_path,
                source,
            } => Self::IndexLoadFailed {
                path,
                index_path,
                source,
            },
            AlignmentIndexError::UnsupportedExtension { path } => {
                Self::UnsupportedAlignmentExtension { path }
            }
            AlignmentIndexError::MixedAlignmentFileFormats {
                first_path,
                first_format,
                other_path,
                other_format,
            } => Self::MixedAlignmentFormats {
                first_path,
                first_format,
                other_path,
                other_format,
            },
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
        &args.alignment_files,
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
    let alignment_cfg = crate::pop_var_caller::cli::alignment_config_from_args(stage1);
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

    // 4. Harvest the canonical contig list + sample name via the
    //    streaming `AlignmentMergedReader::new`. Validates per-file headers
    //    (sort order, `@RG SM`, `@SQ M5`), reconciles across CRAMs,
    //    and cross-checks against the FASTA's `.fai`. The reader is
    //    dropped immediately afterwards — per-chrom workers open
    //    their own readers via `AlignmentMergedReader::query`. The double
    //    open per CRAM (once here for validation, once per worker
    //    for query) is intentional: it keeps the validation surface
    //    centralised and is small relative to the per-contig decode.
    let canonical_contigs;
    let sample_name;
    {
        let reader = crate::bam::alignment_input::AlignmentMergedReader::new(
            &args.alignment_files,
            &args.reference,
            alignment_cfg,
        )
        .map_err(PileupCliError::AlignmentInput)?;
        canonical_contigs = reader.contigs().clone();
        sample_name = reader.sample_name().to_string();
    }
    let chromosomes = contigs_to_parsed(&canonical_contigs)?;

    // 5. Load per-input SAM headers + alignment indexes once. The
    //    headers parsed here, plus the `Arc<crai::Index>` payloads,
    //    are cloned-by-Arc across rayon workers below — one atomic
    //    per worker, no per-worker disk parse.
    let headers = load_per_input_headers(&args.alignment_files)?;
    let indexes = load_per_input_indexes(&args.alignment_files)?;

    // 6. Cohort metadata + writer config template. Workers clone
    //    these per fragment; the template's `output` field is
    //    overridden per worker so each writes to its own fragment
    //    path.
    let metadata = CohortMetadata {
        sample_names: vec![sample_name.clone()],
        contigs: chromosomes.clone(),
        tool_string: format!("pop_var_caller {}", env!("CARGO_PKG_VERSION")),
        command_line: current_command_line(),
    };
    let writer_cfg_template = WriterConfig::new(args.output.clone()).with_emit_gp(cohort.emit_gp);

    // 7. Per-chrom fragment paths under a TempDir placed next to the
    //    final output. Mirror of the cohort `var-calling` driver — the
    //    final atomic rename stays on one filesystem and the
    //    `TempDir` RAII-drops the fragments on Ok and on panic.
    let output_parent = args
        .output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| std::path::Path::new("."));
    let frags_dir =
        TempDir::new_in(output_parent).map_err(|source| VarCallingFromBamCliError::ScratchDir {
            parent: output_parent.to_path_buf(),
            source,
        })?;
    let frag_ext = if path_is_bgzf(&args.output) {
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

    // 8. Pipeline params (cloned per worker; `CohortPipelineParams`
    //    is `Clone` — cheap configs only).
    let pipeline_params = CohortPipelineParams {
        no_complexity_filter: args.no_complexity_filter,
        dust_cfg,
        grouper_cfg,
        per_group_cfg,
        posterior_cfg,
        min_qual_phred: cohort.min_qual_phred,
        min_alt_obs_per_sample: cohort.min_alt_obs_per_sample,
        no_mapq_diff_filter: cohort.no_mapq_diff_filter,
        min_mapq_diff_t: cohort.min_mapq_diff_t,
    };

    // 9. Parallel per-chrom drive. Each worker owns its own readers,
    //    BAQ chunk pool, walker, fetchers, and writer — no shared
    //    mutable state across workers. The rayon pool is the
    //    process-global one configured in step 1.
    use rayon::prelude::*;
    let per_chrom_results: Vec<Result<(u32, CohortDriveStats), VarCallingFromBamCliError>> =
        chromosomes
            .par_iter()
            .enumerate()
            .map(|(cid, _contig)| {
                process_one_chromosome_from_bam(
                    cid as u32,
                    &args.alignment_files,
                    &headers,
                    &indexes,
                    canonical_contigs.clone(),
                    sample_name.clone(),
                    chromosomes.clone(),
                    &args.reference,
                    alignment_cfg,
                    baq_cfg,
                    walker_cfg,
                    stage1.baq_chunk_size,
                    stage1.no_baq,
                    fragment_paths[cid].clone(),
                    metadata.clone(),
                    writer_cfg_template.clone(),
                    pipeline_params.clone(),
                )
            })
            .collect();

    // 10. Fail-fast surfacing of the first error after all workers
    //     joined. Mirror of the cohort path's discipline.
    let mut total_stats = CohortDriveStats::default();
    for r in per_chrom_results {
        let (_cid, stats) = r?;
        total_stats.records_written += stats.records_written;
        total_stats.records_unconverged += stats.records_unconverged;
        total_stats.records_dropped_hom_ref += stats.records_dropped_hom_ref;
        total_stats.records_dropped_low_qual += stats.records_dropped_low_qual;
        total_stats.records_dropped_low_alt_obs += stats.records_dropped_low_alt_obs;
        total_stats.records_dropped_low_mapq_diff_t += stats.records_dropped_low_mapq_diff_t;
    }

    // 11. Concat fragments in contig-table order. `concat_fragments`
    //     owns the atomic rename + parent directory fsync; on Ok the
    //     file is durable at `args.output`. `frags_dir` RAII-drops
    //     the scratch tempdir + fragments here.
    concat_fragments(&args.output, &fragment_paths)?;

    // 12. Run summary. Per-Stage-1 `FilterCounts` / `BaqSkipCounts`
    //     are no longer surfaced because there is no shared
    //     accumulator across the per-chrom workers; reintroducing
    //     them would require an explicit cross-worker reducer
    //     (tracked as a follow-up). The records-emitted summary is
    //     authoritative for what landed in the VCF.
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
        sample_name, total_stats.records_written, emnoconv_note,
    );
    Ok(())
}

/// Open each input once, read its SAM header, and return the
/// parsed headers wrapped in `Arc` for cross-worker sharing.
/// Mirrors what
/// [`crate::bam::alignment_input::AlignmentMergedReader::query`] does
/// per-worker, lifted to driver-startup time so workers only pay an
/// `Arc::clone`. Dispatches on the input file's extension via the
/// shared per-format helpers `read_cram_header_only` /
/// `read_bam_header_only` so the CRAM-version gate (and any future
/// format-specific magic-byte check) stays in one place.
fn load_per_input_headers(
    alignment_files: &[PathBuf],
) -> Result<Vec<Arc<noodles_sam::Header>>, VarCallingFromBamCliError> {
    use crate::bam::bam_input::read_bam_header_only;
    use crate::bam::cram_input::read_cram_header_only;
    use crate::bam::errors::AlignmentInputError;
    use crate::bam::index_preflight::AlignmentFileKind;

    let mut headers = Vec::with_capacity(alignment_files.len());
    for path in alignment_files {
        let kind = AlignmentFileKind::from_path(path).ok_or_else(|| {
            PileupCliError::AlignmentInput(AlignmentInputError::UnsupportedExtension {
                path: path.clone(),
            })
        })?;
        let header = match kind {
            AlignmentFileKind::Cram => {
                read_cram_header_only(path).map_err(PileupCliError::AlignmentInput)?
            }
            AlignmentFileKind::Bam => {
                read_bam_header_only(path).map_err(PileupCliError::AlignmentInput)?
            }
        };
        headers.push(Arc::new(header));
    }
    Ok(headers)
}

/// Load each input's alignment index from the canonical sibling
/// path. Pre-flight (step 0 of `run_var_calling_from_bam`) has
/// already confirmed each index exists or built one when
/// `--build-map-file-index` was set, so a failure here is genuinely
/// I/O. Failures route through the existing
/// `From<AlignmentIndexError>` bridge for `VarCallingFromBamCliError`
/// (M3 + M2: the loader is typed; the previous `io::Error::other`
/// string-flattening is gone).
fn load_per_input_indexes(
    alignment_files: &[PathBuf],
) -> Result<Vec<crate::bam::index_preflight::AlignmentIndex>, VarCallingFromBamCliError> {
    let mut indexes = Vec::with_capacity(alignment_files.len());
    for path in alignment_files {
        let index = crate::bam::index_preflight::load_alignment_index(path)?;
        indexes.push(index);
    }
    Ok(indexes)
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
/// [`AlignmentMergedReader::query`] is owned on this function's stack;
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
#[allow(clippy::too_many_arguments)] // bundling would obscure the borrow chain
fn process_one_chromosome_from_bam(
    chrom_id: u32,
    alignment_files: &[PathBuf],
    headers: &[Arc<noodles_sam::Header>],
    indexes: &[crate::bam::index_preflight::AlignmentIndex],
    canonical_contigs: ContigList,
    sample_name: String,
    all_chromosomes: Vec<ParsedChromosome>,
    reference: &std::path::Path,
    alignment_cfg: crate::bam::alignment_input::AlignmentMergedReaderConfig,
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
    let mut reader = crate::bam::alignment_input::AlignmentMergedReader::query(
        alignment_files,
        reference,
        canonical_contigs.clone(),
        sample_name.clone(),
        headers,
        indexes,
        &chrom_entry.name,
        alignment_cfg,
    )
    .map_err(PileupCliError::AlignmentInput)?;

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
            PileupCliError::AlignmentInput(e),
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
            |source| VarCallingFromBamCliError::RefFetcher {
                contig: chrom_entry.name.clone(),
                source: io::Error::other(format!("{source}")),
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
}
