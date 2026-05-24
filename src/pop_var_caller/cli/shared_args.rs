//! Shared clap-derive sub-structs for the `pop_var_caller`
//! subcommands. Closes **M10** from the 2026-05-19 cohort CLI
//! review — the three top-level args structs ([`PileupArgs`],
//! [`VarCallingArgs`], [`VarCallingFromBamArgs`]) used to duplicate
//! ~30 fields verbatim, three places for the same flag-set.
//!
//! Two shared sub-structs cover the duplication:
//!
//! - [`Stage1Args`] — every CRAM-input / BAQ / walker knob; flattened
//!   into `PileupArgs` and `VarCallingFromBamArgs` (the two
//!   subcommands that drive Stage 1).
//! - [`CohortPipelineArgs`] — every cohort-pipeline knob (DUST,
//!   variant grouping, per-group merger, posterior engine, VCF
//!   writer, ploidy); flattened into `VarCallingArgs` and
//!   `VarCallingFromBamArgs` (the two subcommands that drive Stages
//!   3 – 6).
//!
//! At the user-visible CLI surface, clap-derive's `#[command(flatten)]`
//! is transparent: the flags appear in `--help` exactly as before,
//! and parsing accepts the same CLI grammar. Only the access shape
//! changes internally: `args.foo` becomes `args.stage1.foo` or
//! `args.cohort.foo`.
//!
//! [`PileupArgs`]: super::super::cli::PileupArgs
//! [`VarCallingArgs`]: super::super::var_calling::VarCallingArgs
//! [`VarCallingFromBamArgs`]: super::super::var_calling_from_bam::VarCallingFromBamArgs

use clap::Args;

use crate::bam::cram_input::{
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
use crate::pop_var_caller::cli::parse_mismatch_fraction;
use crate::pop_var_caller::cli::parsers;
use crate::pop_var_caller::cohort_driver::{
    DEFAULT_MIN_ALT_OBS_PER_SAMPLE, DEFAULT_MIN_MAPQ_DIFF_T, DEFAULT_MIN_QUAL_PHRED,
};

/// Tiny shim so `clap`'s `default_value_t` (which needs a `Display`
/// expression evaluated at attribute time) can read
/// [`DEFAULT_MIN_MAPQ_DIFF_T`] without dragging the path through
/// every CLI helper that reads this constant.
const fn pop_var_caller_default_min_mapq_diff_t() -> f32 {
    DEFAULT_MIN_MAPQ_DIFF_T
}
use crate::var_calling::dust_filter::{DEFAULT_DUST_THRESHOLD, DEFAULT_DUST_WINDOW};
use crate::var_calling::per_group_merger::{DEFAULT_MAX_ALLELES_PER_RECORD, DEFAULT_PLOIDY};
use crate::var_calling::posterior_engine::{
    DEFAULT_COMPOUND_ALT_PSEUDOCOUNT, DEFAULT_CONVERGENCE_THRESHOLD,
    DEFAULT_INBREEDING_COEFFICIENT, DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MAX_GQ_PHRED,
    DEFAULT_MAX_ITERATIONS, DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT,
};
use crate::var_calling::variant_grouping::DEFAULT_MAX_VARIANT_GROUP_SPAN;
use crate::vcf::DEFAULT_EMIT_GP;

/// Stage 1 knobs (CRAM-input filters, BAQ HMM, pileup walker).
/// Flattened into [`PileupArgs`](super::super::cli::PileupArgs) and
/// [`VarCallingFromBamArgs`](super::super::var_calling_from_bam::VarCallingFromBamArgs).
#[derive(Debug, Args, Clone)]
pub struct Stage1Args {
    // ===== Common flags (visible in `-h`) =====================
    /// Drop reads with MAPQ < N. `0` admits everything.
    #[arg(long, default_value_t = DEFAULT_MIN_MAPQ)]
    pub min_mapq: u8,

    /// Skip the BAQ HMM. `bq_baq` becomes a copy of the raw CRAM
    /// QUAL.
    #[arg(long)]
    pub no_baq: bool,

    // ===== Advanced — CRAM input filters ======================
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

    // ===== Advanced — BAQ HMM =================================
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

    // ===== Advanced — Pileup walker ===========================
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
        default_value_t = DEFAULT_MAX_ACTIVE_READS,
        help_heading = "Advanced — Pileup walker",
    )]
    pub max_active_reads: u32,
}

/// Cohort-pipeline knobs (DUST filter, variant grouping, per-group
/// merger, posterior engine, ploidy, VCF writer). Flattened into
/// [`VarCallingArgs`](super::super::var_calling::VarCallingArgs) and
/// [`VarCallingFromBamArgs`](super::super::var_calling_from_bam::VarCallingFromBamArgs).
#[derive(Debug, Args, Clone)]
pub struct CohortPipelineArgs {
    // ===== Common (visible in `-h`) ===========================
    /// Cohort-wide ploidy.
    #[arg(long, default_value_t = DEFAULT_PLOIDY, value_parser = parsers::parse_ploidy)]
    pub ploidy: u8,

    // ===== Advanced — DUST filter ==============================
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

    // ===== Advanced — Variant grouping =========================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_VARIANT_GROUP_SPAN,
        value_parser = parsers::parse_var_group_max_span,
        help_heading = "Advanced — Variant grouping",
    )]
    pub var_group_max_span: u32,

    // ===== Advanced — Per-group merger =========================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_ALLELES_PER_RECORD,
        value_parser = parsers::parse_max_alleles,
        help_heading = "Advanced — Per-group merger",
    )]
    pub max_alleles_per_var: usize,

    // ===== Advanced — Posterior engine =========================
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

    // ===== Advanced — VCF writer ===============================
    /// Drop records with site-level `QUAL` strictly below this
    /// (phred). `0` disables the filter. Default matches GATK
    /// HaplotypeCaller's emission gate; caveat: our `QUAL` is a
    /// product over per-sample hom-ref posteriors, which inflates
    /// with cohort size, so the same threshold means different
    /// things at different N.
    #[arg(
        long = "min-qual",
        hide_short_help = true,
        default_value_t = DEFAULT_MIN_QUAL_PHRED,
        value_parser = parsers::parse_min_qual_phred,
        help_heading = "Advanced — VCF writer",
    )]
    pub min_qual_phred: f64,

    /// Drop records before the EM where no ALT allele has
    /// `max(num_obs across samples) >= N`. `0` and `1` disable the
    /// filter (every candidate has ≥ 1 supporting read by
    /// construction). Default `2` empirically cuts false positives
    /// against GATK by ~70% at < 1 pp recall cost — single-alt-read
    /// candidates dominate the noise tail (89.9% on the multichrom
    /// validation cohort). GATK's analogue is `--min-pruning = 2` on
    /// the assembly graph.
    #[arg(
        long = "min-alt-obs-per-sample",
        hide_short_help = true,
        default_value_t = DEFAULT_MIN_ALT_OBS_PER_SAMPLE,
        value_parser = parsers::parse_min_alt_obs_per_sample,
        help_heading = "Advanced — Per-group merger",
    )]
    pub min_alt_obs_per_sample: u32,

    /// Skip the Welch's-t MAPQ-difference drop entirely. The
    /// `INFO/MQRef`, `MQAlt`, `MQDiff`, and `MQDiffT` annotations
    /// still emit; setting this keeps records that would otherwise
    /// be dropped for having alt-supporting reads at systematically
    /// lower MAPQ than ref-supporting reads (the multi-mapper
    /// fingerprint).
    #[arg(
        long = "no-mapq-diff-filter",
        hide_short_help = true,
        default_value_t = false,
        help_heading = "Advanced — MAPQ filter"
    )]
    pub no_mapq_diff_filter: bool,

    /// Threshold for the Welch's-t MAPQ-difference drop. Records
    /// with at least one ALT whose cohort-pooled Welch's t (ALT vs
    /// REF MAPQ) is below this value AND both sides have ≥ 3
    /// supporting reads are dropped before reaching the writer.
    /// Pass `-inf` to disable; pass `--no-mapq-diff-filter` for the
    /// equivalent boolean. Default `-3.0` empirically catches 25 %
    /// of GATK-extreme multi-mapper sites at a 2.0 % false-flag
    /// rate on GATK-clean sites; see the per_allele_mapq_tracking
    /// implementation plan for the cross-validation.
    #[arg(
        long = "min-mapq-diff-t",
        hide_short_help = true,
        default_value_t = pop_var_caller_default_min_mapq_diff_t(),
        value_parser = parsers::parse_min_mapq_diff_t,
        help_heading = "Advanced — MAPQ filter",
    )]
    pub min_mapq_diff_t: f32,

    /// Emit `GP` (genotype posteriors) `FORMAT` per sample. Off by
    /// default — `GP` is `Number=G`, so the per-sample cell grows as
    /// `(ploidy + n_alleles − 1) choose ploidy` (21 floats at
    /// ploidy=2, n_alleles=6). Opt in when posteriors are wanted on
    /// disk; a 1000-sample × 100 K-variant cohort with `emit_gp` adds
    /// ~8 GiB to the VCF.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_EMIT_GP,
        help_heading = "Advanced — VCF writer",
    )]
    pub emit_gp: bool,
}
