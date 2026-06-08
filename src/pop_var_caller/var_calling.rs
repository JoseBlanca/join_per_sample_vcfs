//! `pop_var_caller var-calling` — cohort `.psp` → multi-sample VCF.
//!
//! Wires Stages 3–6 of the pipeline:
//!
//!   PspReader(s) → PerPositionMerger → DustFilter → VariantGrouper
//!   → PerGroupMerger → PosteriorEngine → CohortVcfWriter
//!
//! This CLI layer validates inputs (chromosome agreement, FASTA-vs-`.psp`
//! MD5), loads any `--contamination-estimates`, and hands off to the
//! re-architected record-streaming pipeline that runs those stages — the
//! wiring lives in [`crate::var_calling::pipeline::run_var_calling`].
//!
//! Contamination plumbing:
//!
//! - With `--contamination-estimates <FILE>`: the artefact is loaded
//!   and reconciled against the cohort sample names (extras are
//!   tolerated, absences are not), then handed to the posterior engine.
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
    verify_fasta_matches_psp_chromosomes,
};
use crate::pop_var_caller::contamination_artefact::{
    ContaminationArtefact, ContaminationArtefactError,
};
use crate::psp::{PspReadError, PspReader};
use crate::var_calling::contamination_estimation::ContaminationEstimates;
use crate::var_calling::per_position_merger::{PerPositionMergerError, check_chromosome_agreement};
use crate::var_calling::pipeline::PipelineError;

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

    /// Worker threads. If omitted, defaults to all logical cores
    /// (`available_parallelism`). Also sizes the pipeline's caller-thread
    /// pool and the two bounded hand-off queues (depth ≈ 2 × threads), so
    /// it is the primary peak-memory knob.
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

    /// Soft target for the number of variable (cohort-variant) positions
    /// per chunk — the producer accumulates records until it reaches this
    /// many, then cuts at the next safe gap. Primarily a wall trade (finer
    /// chunks load-balance the producer→caller pipeline better); peak RSS is
    /// flat across chunk size at sane targets, so it is not a memory knob.
    /// The emitted VCF is independent of it (cuts always fall on clean group
    /// boundaries). `0` (the default) selects the built-in default of 128;
    /// any positive value overrides. The resolved value is printed in the
    /// startup log.
    #[arg(long, default_value_t = 0)]
    pub target_variants_per_chunk: u32,

    /// Low-memory mode: trade wall time for a smaller peak RSS by turning off
    /// the speed-for-memory optimisations the default path uses. Currently this
    /// disables the producer's straddler decode cache — a segment that spans a
    /// chunk boundary is re-decompressed for each chunk instead of being held
    /// decompressed in memory (whose cost scales with the cohort size). Expect a
    /// few percent slower for a peak-RSS reduction that grows with sample count.
    /// A bucket for further low-memory toggles. Output is byte-identical to the
    /// default path either way.
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

    /// Surfaced by the re-architected pipeline — wraps every per-stage error
    /// (config build, .psp decode, REF/dust, grouping/merge/EM, VCF write) so
    /// they flow through one CLI-side variant.
    #[error("pipeline: {0}")]
    Pipeline(#[from] PipelineError),

    #[error("contamination artefact: {0}")]
    ContamArtefact(#[from] ContaminationArtefactError),

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
/// 9. Drive the record-streaming pipeline
///    ([`crate::var_calling::pipeline::run_var_calling`]): a producer
///    that materialises one chunk × N samples at a time, partitions it
///    into windows, and runs the per-group merger + posterior EM on
///    each window. Records stream straight into a single
///    [`CohortVcfWriter`](crate::vcf::CohortVcfWriter) in genomic order
///    — no per-chromosome fragments, no concat.
/// 10. Print a one-shot stderr run-summary block (includes
///     `effective_threads` + the `requested … capped by N chromosomes`
///     parenthetical when the soft cap bit).
pub fn run_var_calling(args: &VarCallingArgs) -> Result<(), VarCallingCliError> {
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

    // 5. Load `--contamination-estimates` (reconciled to the cohort order;
    //    `None` → no contamination). The per-stage EM / grouper / merger
    //    configs are built inside the pipeline from the args.
    let contamination = load_contamination(args.contamination_estimates.as_deref(), &sample_names)?;

    // 6. Cross-check FASTA bytes against the .psp per-contig MD5s
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

    // 7. Run the re-architected pipeline (producer → W callers → writer). It
    //    opens its own per-sample readers, builds every per-stage config from
    //    the args, applies `--regions`, and writes the VCF; it returns the
    //    run-level counters for the summary.
    let stats = crate::var_calling::pipeline::run_var_calling(args, contamination)?;

    // 8. Stderr summary.
    print_run_summary(&sample_names, stats, args.threads, chromosomes.len());
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
    stats: crate::var_calling::vcf_writer::WriterStats,
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
    eprintln!(
        "var-calling: n_samples={} records_emitted={} effective_threads={}{}{}{}{}{}{}{}",
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
