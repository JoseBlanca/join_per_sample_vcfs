//! `pop_var_caller var-calling` — cohort `.psp` → multi-sample VCF.
//!
//! Wires Stages 3–6 of the pipeline:
//!
//!   PspReader(s) → PerPositionMerger → DustFilter → VariantGrouper
//!   → PerGroupMerger → PosteriorEngine → CohortVcfWriter
//!
//! The DUST filter is bypassed when `--no-complexity-filter` is set;
//! both branches converge on the same downstream chain via a
//! `Box<dyn Iterator>` adapter that lifts the upstream error into
//! [`GrouperError`] — the wiring lives in
//! [`crate::var_calling::driver::drive_cohort_chunked`].
//!
//! Contamination plumbing:
//!
//! - With `--contamination-estimates <FILE>`: the artefact is loaded
//!   and reconciled against the cohort sample names (extras are
//!   tolerated, absences are not), then handed to the posterior engine
//!   via [`PosteriorEngineConfig::contamination`].
//! - Without: the engine runs in "no contamination" mode and the
//!   `contamination` field stays `None`.
//!
//! Plan: `doc/devel/implementation_plans/pop_var_caller_cohort_cli.md`
//! §"Subcommand: var-calling".

use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};

use clap::Args;
use thiserror::Error;

use crate::pop_var_caller::common::{
    DEFAULT_BUFFERED_IO_CAPACITY, FastaVerifyError, basename, configure_rayon_pool,
    current_command_line, verify_fasta_matches_psp_chromosomes,
};
use crate::pop_var_caller::contamination_artefact::{
    ContaminationArtefact, ContaminationArtefactError,
};
use crate::psp::{PspReadError, PspReader};
use crate::var_calling::contamination_estimation::ContaminationEstimates;
use crate::var_calling::dust_filter::{DustFilterConfig, DustFilterError};
use crate::var_calling::per_group_merger::{
    DEFAULT_BATCH_SIZE, PerGroupMergerConfig, PerGroupMergerError,
};
use crate::var_calling::per_position_merger::{PerPositionMergerError, check_chromosome_agreement};
use crate::var_calling::posterior_engine::{
    PosteriorEngineConfig, PosteriorEngineConfigError, PosteriorEngineError,
};
use crate::var_calling::variant_grouping::{GrouperConfig, GrouperConfigError, GrouperError};
use crate::var_calling::{
    ChunkDriverError, ChunkDriverParams, ChunkDriverStats, ChunkSizingParams,
    DownstreamFilterParams, drive_cohort_chunked,
};
use crate::vcf::{CohortMetadata, VcfWriteError, WriterConfig};

/// Default desired number of variable variants per materialised block
/// (the `BlockIterator`'s accumulation target) when
/// `--target-variants-per-chunk` is left at `0`. Sets block size: the
/// iterator pulls .psp chunks until it has this many variable variants,
/// then cuts at a clean group boundary. Output is independent of it;
/// it trades block size for memory (resident ≈ target × n_samples).
/// Tunable; validate against the high-N synthetic scaling benchmark.
const DEFAULT_DESIRED_VARIANTS_PER_BLOCK: u32 = 1024;

// ---------------------------------------------------------------------
// CLI surface
// ---------------------------------------------------------------------

/// Arguments accepted by `pop_var_caller var-calling`.
///
/// Cohort-pipeline knobs (DUST, grouper, per-group merger,
/// posterior engine, ploidy, VCF writer) live in the flattened
/// [`CohortPipelineArgs`](crate::pop_var_caller::cli::shared_args::CohortPipelineArgs)
/// sub-struct so the `var-calling-from-bam` subcommand reuses the
/// same surface — M10 from the 2026-05-19 cohort CLI review.
#[derive(Debug, Args, Clone)]
pub struct VarCallingArgs {
    // ===== Common (visible in `-h`) ===========================
    /// Reference FASTA used to produce the `.psp` files. Basename
    /// must match the `reference` field every `.psp` header carries.
    #[arg(long)]
    pub reference: PathBuf,

    /// Output VCF path. Extension picks the sink kind:
    /// `.vcf.gz` / `.vcf.bgz` → bgzf, anything else → plain text.
    /// Written via atomic `<output>.tmp` → `fs::rename`.
    #[arg(long)]
    pub output: PathBuf,

    /// BED file restricting variant calling to the listed regions; the
    /// rest of the genome is ignored. Coordinates are 0-based half-open
    /// (standard BED). Without this flag every position in the `.psp`
    /// files is called. The `.psp` block index is used to seek to the
    /// regions, so a sparse BED is cheap.
    #[arg(long)]
    pub regions: Option<PathBuf>,

    /// Worker threads for the rayon pool. If omitted, rayon picks
    /// the default (all logical cores).
    #[arg(long)]
    pub threads: Option<usize>,

    /// Contamination-estimates artefact produced by
    /// `pop_var_caller estimate-contamination`. When omitted, every
    /// sample's `c_s` is treated as `0` (no contamination correction).
    #[arg(long)]
    pub contamination_estimates: Option<PathBuf>,

    /// Skip the low-complexity (sdust) filter entirely.
    #[arg(long)]
    pub no_complexity_filter: bool,

    /// Soft lower bound on the post-filter variant count per chunk.
    /// The loader grows each chunk's BP span adaptively until the
    /// kept-position count crosses this target — decoupling worker
    /// workload from PSP block size and per-region variant density.
    /// Pass `0` (the default) to keep the legacy BP-only loop (one
    /// pull of `--chunk-genomic-span` BP per chunk).
    #[arg(long, default_value_t = 0)]
    pub target_variants_per_chunk: u32,

    /// Low-memory mode: summarise each covered interval and re-read it to
    /// materialise only the kept positions, instead of folding every
    /// sample's window records into memory at once. Holds ≤ pool-width
    /// decoded blocks plus a compact per-position summary, trading a
    /// second pass over the `.psp` files for a much smaller peak RSS.
    /// Output is byte-identical to the default path.
    #[arg(long, default_value_t = false)]
    pub low_memory: bool,

    /// One or more cohort `.psp` files.
    #[arg(required = true)]
    pub psp_files: Vec<PathBuf>,

    // ===== Cohort pipeline (shared with var-calling-from-bam) ==
    #[command(flatten)]
    pub cohort: crate::pop_var_caller::cli::shared_args::CohortPipelineArgs,
}

// ---------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------

/// Errors surfaced by [`run_var_calling`]. Wraps every upstream module
/// error so the CLI's `format_error_chain` can render the cause.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum VarCallingCliError {
    #[error("io: {0}")]
    Io(#[from] io::Error),

    #[error("psp reader: {0}")]
    PspReader(#[from] PspReadError),

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

    /// Surfaced by [`drive_cohort_chunked`] — wraps every per-stage
    /// error from the chunk-loop driver (loader, pre-pass, partition,
    /// worker, vcf writer, etc.) so they all flow through one
    /// CLI-side variant.
    #[error("chunk driver: {0}")]
    ChunkDriver(#[from] ChunkDriverError),

    #[error("contamination artefact: {0}")]
    ContamArtefact(#[from] ContaminationArtefactError),

    #[error("regions BED: {0}")]
    Bed(#[from] crate::regions::BedError),

    #[error("grouper config: {0}")]
    GrouperConfig(#[from] GrouperConfigError),

    #[error("per-group merger config: {0}")]
    PerGroupConfig(#[from] crate::var_calling::per_group_merger::PerGroupMergerConfigError),

    #[error("posterior engine config: {0}")]
    PosteriorConfig(#[from] PosteriorEngineConfigError),

    // `DustFilterConfig::new` also returns `DustFilterError`, which
    // funnels through the `Dust` variant above — no separate
    // `DustConfig` slot needed.
    /// A `.psp` file's `reference` field doesn't match the basename
    /// of `--reference`. Surfaces before any record is read.
    #[error(
        "psp {psp}: reference mismatch — header has `{psp_ref}`, CLI passed `{supplied_ref}` \
         (basename comparison)"
    )]
    ReferenceMismatch {
        psp: PathBuf,
        psp_ref: String,
        supplied_ref: String,
    },

    /// `rayon::ThreadPoolBuilder::build_global()` already ran in this
    /// process. The binary calls it at most once.
    #[error("rayon thread pool already initialised — refusing to override")]
    RayonAlreadyConfigured,

    /// FASTA contig bytes don't match the `.psp` header's
    /// per-contig MD5. The basename cross-check (variant
    /// `ReferenceMismatch`) is necessary but not sufficient — this
    /// variant catches "right filename, wrong bytes" (e.g. the user
    /// has a same-named FASTA pointing at a different genome
    /// build). M5 follow-up from the 2026-05-19 cohort CLI review.
    #[error(
        "fasta `{contig}` bytes don't match the .psp header MD5 — \
         fasta has `{fasta_md5}`, .psp has `{psp_md5}`"
    )]
    FastaContigMd5Mismatch {
        contig: String,
        fasta_md5: String,
        psp_md5: String,
    },

    /// Failed to read a contig from the FASTA referenced by
    /// `--reference` — the contig declared by the `.psp` is missing
    /// (or unreadable) from the FASTA on disk. Distinguished from
    /// `FastaContigMd5Mismatch` because the user-facing remediation
    /// differs: the FASTA either doesn't contain the contig at all
    /// or the `.fai` index is stale.
    #[error("fasta `{contig}`: failed to read for MD5 verification — {source}")]
    FastaContigFetchFailed {
        contig: String,
        #[source]
        source: io::Error,
    },
}

// ---------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------

/// Run the cohort var-calling pipeline and emit a multi-sample VCF.
///
/// Pipeline (per the cohort CLI plan §"Subcommand: var-calling" +
/// the 2026-05-20 per-chromosome-parallel implementation plan):
///
/// 1. Size the global rayon pool when `--threads` is set.
/// 2. Open every `.psp` with [`PspReader`] for header validation
///    (the validated readers are then dropped — each per-chrom
///    worker re-opens its own set).
/// 3. Cross-check `header.reference` against the basename of
///    `--reference`. **Basename-only — the FASTA bytes are
///    additionally MD5-verified at step 7a.**
/// 4. Cross-check chromosomes across readers
///    ([`check_chromosome_agreement`]). Compares per-contig MD5s
///    *between readers* only.
/// 5. Load `--contamination-estimates` if supplied and reconcile
///    sample names; absence → `c_s = 0` for every sample.
/// 6. Build + validate every per-stage config from the CLI args.
/// 7. Verify the FASTA's per-contig bytes against the `.psp`
///    per-contig MD5s by streaming each contig through `Md5::update`
///    in parallel — peak memory is one 64 KiB window per worker, no
///    contig-sized allocations. No process-wide reference fetcher is
///    built — the chunk driver constructs a per-chromosome
///    `StreamingChromRefFetcher` as it advances.
/// 8. Cohort metadata + writer-config template.
/// 9. Drive the chunk loop ([`drive_cohort_chunked`]): a sequential
///    per-chromosome outer loop with a within-chromosome chunk loop
///    that materialises one chunk × N samples at a time, partitions it
///    into windows, and runs the per-group merger + posterior EM on
///    each window. Records stream straight into a single
///    [`CohortVcfWriter`](crate::vcf::CohortVcfWriter) in genomic order
///    — no per-chromosome fragments, no concat.
/// 10. Print a one-shot stderr run-summary block (includes
///     `effective_threads` + the `requested … capped by N chromosomes`
///     parenthetical when the soft cap bit).
pub fn run_var_calling(args: &VarCallingArgs) -> Result<(), VarCallingCliError> {
    // TEMP (P6 measurement, re-architecture): route to the new pipeline when
    // `PVC_NEW_PIPELINE` is set, so the same CLI invocation A/B-tests both paths
    // for byte-identity + RSS/wall. Removed at the P7 swap (which repoints the
    // CLI to the new entry directly). Off by default → no behaviour change.
    if std::env::var_os("PVC_NEW_PIPELINE").is_some() {
        crate::var_calling_new::pipeline::run_var_calling(args)
            .unwrap_or_else(|e| panic!("new pipeline failed: {e}"));
        return Ok(());
    }

    let cohort = &args.cohort;

    // 1. Rayon pool — idempotent under the silent-no-op policy
    //    (M13, locked 2026-05-19); first caller wins, subsequent
    //    callers' --threads is silently ignored.
    configure_rayon_pool(args.threads).map_err(|_| VarCallingCliError::RayonAlreadyConfigured)?;

    // 2. Open every .psp once for header validation. Per-chrom
    //    workers re-open their own readers below; keeping the
    //    parent's readers alive would multiply concurrent open
    //    file descriptors by ~N_chromosomes for no read-side win.
    let mut readers: Vec<PspReader<BufReader<File>>> = Vec::with_capacity(args.psp_files.len());
    for path in &args.psp_files {
        let file = File::open(path).map_err(VarCallingCliError::Io)?;
        let buf = BufReader::with_capacity(DEFAULT_BUFFERED_IO_CAPACITY, file);
        let reader = PspReader::new(buf)?;
        readers.push(reader);
    }

    // 3. Cross-check references. Basename comparison only — the
    //    FASTA bytes get an MD5 cross-check at step 7a.
    let supplied_ref = basename(&args.reference);
    for (path, reader) in args.psp_files.iter().zip(readers.iter()) {
        if reader.header().reference != supplied_ref {
            return Err(VarCallingCliError::ReferenceMismatch {
                psp: path.clone(),
                psp_ref: reader.header().reference.clone(),
                supplied_ref: supplied_ref.clone(),
            });
        }
    }

    // 4. Sample names + chromosome agreement.
    let sample_names: Vec<String> = readers.iter().map(|r| r.header().sample.clone()).collect();
    let chromosomes = check_chromosome_agreement(&readers)?;
    // Validation done; drop the parent readers so workers can own
    // their own file handles without inflating the concurrent fd
    // count. Plan §"RLIMIT_NOFILE bump" tracks the v1 ceiling at
    // N_samples × N_chromosomes.
    drop(readers);

    // 5. Build + validate every per-stage config.
    let dust_cfg = DustFilterConfig::new(cohort.complexity_window, cohort.complexity_threshold)?;
    let grouper_cfg = GrouperConfig::new(cohort.var_group_max_span)?;
    let per_group_cfg = PerGroupMergerConfig::new(
        cohort.ploidy,
        cohort.max_alleles_per_var,
        cohort.max_alleles_lh_calc,
        DEFAULT_BATCH_SIZE,
    )?;
    // 5/6. Posterior config + contamination wired through the named-
    //      setter builder chain so out-of-range values surface as
    //      typed errors at construction time.
    let posterior_cfg = PosteriorEngineConfig::new()
        .with_convergence_threshold(cohort.em_convergence_threshold)?
        .with_max_iterations(cohort.em_max_iterations)?
        .with_ref_pseudocount(cohort.ref_pseudocount)?
        .with_snp_alt_pseudocount(cohort.snp_alt_pseudocount)?
        .with_indel_alt_pseudocount(cohort.indel_alt_pseudocount)?
        .with_compound_alt_pseudocount(cohort.compound_alt_pseudocount)?
        .with_fixation_index_default(cohort.inbreeding_coefficient)?
        .with_max_gq_phred(cohort.max_gq_phred)?
        .with_contamination(load_contamination(
            args.contamination_estimates.as_deref(),
            &sample_names,
        )?)?;

    // 7a. Cross-check FASTA bytes against the .psp per-contig MD5s
    //     *before* building the runtime fetcher. Catches "right
    //     basename, wrong genome build". Streams each contig through
    //     `Md5::update` in 64 KiB windows, per-chrom in parallel —
    //     peak resident is one window per worker, zero contig-sized
    //     allocations. No fetcher cache pre-warm (Phase A of
    //     reference_fasta_streaming).
    match verify_fasta_matches_psp_chromosomes(&args.reference, &chromosomes) {
        Ok(()) => {}
        Err(FastaVerifyError::Md5Mismatch {
            contig,
            fasta_md5,
            psp_md5,
        }) => {
            return Err(VarCallingCliError::FastaContigMd5Mismatch {
                contig,
                fasta_md5,
                psp_md5,
            });
        }
        Err(FastaVerifyError::FetchFailed { contig, source }) => {
            return Err(VarCallingCliError::FastaContigFetchFailed { contig, source });
        }
    }

    // 7b. No global fetcher is built here. Each per-chrom worker
    //     constructs its own `StreamingChromRefFetcher` inside
    //     `process_one_chromosome` — no Arc shared across threads,
    //     no Mutex contention. The contig's sliding buffer lives
    //     exactly for the worker's lifetime; peak resident during
    //     the `par_iter` below is `min(threads, n_chroms) × 1 MB`,
    //     independent of contig size.

    // 8. Cohort metadata for the writer header + the writer-config
    //    template inherited by every per-chrom fragment writer.
    let metadata = CohortMetadata {
        sample_names: sample_names.clone(),
        contigs: chromosomes.clone(),
        tool_string: format!("pop_var_caller {}", env!("CARGO_PKG_VERSION")),
        command_line: current_command_line(),
    };
    let writer_cfg_template = WriterConfig::new(args.output.clone()).with_emit_gp(cohort.emit_gp);

    // 9. Drive the chunk loop. The within-chromosome chunk-parallel
    //    rewrite (Phase A — see
    //    `doc/devel/implementation_plans/cohort_within_chromosome_parallel.md`)
    //    runs the entire cohort end-to-end through one
    //    `CohortVcfWriter`; no per-chrom fragments, no concat. The
    //    new driver owns its own per-sample PSP readers + per-chrom
    //    `StreamingChromRefFetcher` + the persistent chunk-loop
    //    scratch.
    let driver_params = ChunkDriverParams {
        no_complexity_filter: args.no_complexity_filter,
        dust_cfg,
        grouper_cfg,
        per_group_cfg,
        posterior_cfg,
        sizing: {
            // Desired variable-variants per block — the BlockIterator's
            // accumulation target. `0` (default) uses the const default;
            // an explicit `--target-variants-per-chunk` overrides. Output
            // is independent of this (the iterator always cuts at clean
            // group boundaries); it trades block size for memory
            // (resident ≈ target × n_samples — the one-chunk-at-a-time
            // bound the rewrite trades genomic span for).
            let target_variants = if args.target_variants_per_chunk == 0 {
                DEFAULT_DESIRED_VARIANTS_PER_BLOCK
            } else {
                args.target_variants_per_chunk
            };
            ChunkSizingParams {
                target_variants_per_chunk: std::num::NonZeroU32::new(target_variants),
            }
        },
        downstream: DownstreamFilterParams {
            min_alt_obs_per_sample: args.cohort.min_alt_obs_per_sample,
            min_qual_phred: args.cohort.min_qual_phred,
            no_mapq_diff_filter: args.cohort.no_mapq_diff_filter,
            min_mapq_diff_t: args.cohort.min_mapq_diff_t,
        },
        low_memory: args.low_memory,
    };
    // Region set: the --regions BED resolved against the cohort
    // chromosomes, or None (whole genome — the driver runs the
    // identical, byte-for-byte path when this is None).
    let region_set = match &args.regions {
        Some(bed_path) => {
            let contig_bounds: Vec<crate::regions::ContigBounds> = chromosomes
                .iter()
                .map(|c| crate::regions::ContigBounds {
                    name: &c.name,
                    length: c.length,
                })
                .collect();
            Some(crate::regions::RegionSet::from_bed_path(
                bed_path,
                &contig_bounds,
            )?)
        }
        None => None,
    };

    let chunk_stats = drive_cohort_chunked(
        &args.psp_files,
        sample_names.clone(),
        chromosomes.clone(),
        &args.reference,
        &args.output,
        metadata,
        writer_cfg_template,
        driver_params,
        region_set.as_ref(),
    )?;

    // 11. Stderr summary.
    print_run_summary(&sample_names, chunk_stats, args.threads, chromosomes.len());
    Ok(())
}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

/// Load and reconcile a contamination-estimates artefact against the
/// cohort sample order. Returns `None` if `path` is `None`.
fn load_contamination(
    path: Option<&Path>,
    sample_names: &[String],
) -> Result<Option<ContaminationEstimates>, VarCallingCliError> {
    let path = match path {
        Some(p) => p,
        None => return Ok(None),
    };
    let artefact = ContaminationArtefact::read(path)?;
    let refs: Vec<&str> = sample_names.iter().map(String::as_str).collect();
    Ok(Some(artefact.to_estimates_for_samples(&refs)?))
}

/// Print the one-shot stderr run summary.
///
/// `effective_threads` is the per-chrom parallelism actually
/// realised — `min(rayon_pool_size, n_chromosomes)`. The
/// parenthetical `(requested R; capped by N chromosomes)` is
/// appended only when the cap actually bit (i.e. the user passed
/// `--threads R` with `R > n_chromosomes`). Matches the convention
/// in samtools / bcftools where over-provisioning the thread pool
/// is a warning-grade fact, not an error.
///
/// `records_emnoconv` is a tally of records whose posterior EM hit
/// the iteration cap and were emitted with `FILTER=EMNoConv`. Only
/// printed when at least one such record was produced; absence keeps
/// the happy-path summary tight.
fn print_run_summary(
    sample_names: &[String],
    stats: ChunkDriverStats,
    requested_threads: Option<usize>,
    n_chromosomes: usize,
) {
    let pool_size = rayon::current_num_threads();
    let effective_threads = pool_size.min(n_chromosomes.max(1));
    let cap_note = match requested_threads {
        Some(req) if req > n_chromosomes && n_chromosomes > 0 => {
            format!(" (requested {req}; capped by {n_chromosomes} chromosomes)")
        }
        _ => String::new(),
    };
    let emnoconv_note = if stats.records_unconverged > 0 {
        format!(
            " records_emnoconv={} (FILTER=EMNoConv; EM iteration cap)",
            stats.records_unconverged,
        )
    } else {
        String::new()
    };
    let dropped_hom_ref_note = if stats.records_dropped_hom_ref > 0 {
        format!(
            " records_dropped_hom_ref={} (no sample carries an ALT)",
            stats.records_dropped_hom_ref,
        )
    } else {
        String::new()
    };
    let dropped_low_qual_note = if stats.records_dropped_low_qual > 0 {
        format!(
            " records_dropped_low_qual={} (QUAL < --min-qual)",
            stats.records_dropped_low_qual,
        )
    } else {
        String::new()
    };
    let dropped_low_alt_obs_note = if stats.records_dropped_low_alt_obs > 0 {
        format!(
            " records_dropped_low_alt_obs={} (no ALT max(num_obs) >= --min-alt-obs-per-sample)",
            stats.records_dropped_low_alt_obs,
        )
    } else {
        String::new()
    };
    let dropped_low_mapq_diff_t_note = if stats.records_dropped_low_mapq_diff_t > 0 {
        format!(
            " records_dropped_low_mapq_diff_t={} (any ALT Welch's t < --min-mapq-diff-t)",
            stats.records_dropped_low_mapq_diff_t,
        )
    } else {
        String::new()
    };
    let lh_cap_note = if stats.lh_cap_groups_skipped > 0 {
        format!(
            " allele_lh_cap_hit: groups_skipped={} alleles_total={} (--max-alleles-lh-calc bound)",
            stats.lh_cap_groups_skipped, stats.lh_cap_alleles_in_skipped,
        )
    } else {
        String::new()
    };
    let chunks_loaded_note = if stats.chunks_loaded > 0 {
        let avg = stats.chunk_variants_total as f64 / stats.chunks_loaded as f64;
        format!(
            " chunks_loaded={} avg_variants_per_chunk={:.1}",
            stats.chunks_loaded, avg,
        )
    } else {
        String::new()
    };
    eprintln!(
        "var-calling: n_samples={} records_emitted={} effective_threads={}{}{}{}{}{}{}{}{}",
        sample_names.len(),
        stats.records_written,
        effective_threads,
        cap_note,
        emnoconv_note,
        dropped_hom_ref_note,
        dropped_low_qual_note,
        dropped_low_alt_obs_note,
        dropped_low_mapq_diff_t_note,
        lh_cap_note,
        chunks_loaded_note,
    );
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_contamination_returns_none_without_path() {
        let names = vec!["NA12878".to_string()];
        let result = load_contamination(None, &names).unwrap();
        assert!(result.is_none());
    }
}
