//! Shared clap-derive sub-structs for the `pop_var_caller`
//! subcommands. Closes **M10** from the 2026-05-19 cohort CLI
//! review — the top-level args structs ([`PileupArgs`],
//! [`VarCallingArgs`]) used to duplicate ~30 fields verbatim, the
//! same flag-set spelled out in each place. (A third struct,
//! `VarCallingFromBamArgs`, also shared them before the direct
//! BAM→VCF path was removed.)
//!
//! Two shared sub-structs cover the duplication:
//!
//! - [`Stage1Args`] — every CRAM-input / BAQ / walker knob; flattened
//!   into `PileupArgs` (the subcommand that drives Stage 1).
//! - [`CohortPipelineArgs`] — every cohort-pipeline knob (DUST,
//!   variant grouping, per-group merger, posterior engine, VCF
//!   writer, ploidy); flattened into `VarCallingArgs` (the subcommand
//!   that drives Stages 3 – 6).
//!
//! At the user-visible CLI surface, clap-derive's `#[command(flatten)]`
//! is transparent: the flags appear in `--help` exactly as before,
//! and parsing accepts the same CLI grammar. Only the access shape
//! changes internally: `args.foo` becomes `args.stage1.foo` or
//! `args.cohort.foo`.
//!
//! [`PileupArgs`]: super::super::cli::PileupArgs
//! [`VarCallingArgs`]: super::super::var_calling::VarCallingArgs

use clap::Args;

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
use crate::pop_var_caller::cli::parse_mismatch_fraction;
use crate::pop_var_caller::cli::parsers;
use crate::var_calling::dust_filter::{DEFAULT_DUST_THRESHOLD, DEFAULT_DUST_WINDOW};
use crate::var_calling::per_group_merger::{
    DEFAULT_MAX_ALLELES_LH_CALC, DEFAULT_MAX_ALLELES_PER_RECORD, DEFAULT_PLOIDY,
};
use crate::var_calling::posterior_engine::{
    DEFAULT_COMPOUND_ALT_PSEUDOCOUNT, DEFAULT_CONVERGENCE_THRESHOLD,
    DEFAULT_INBREEDING_COEFFICIENT, DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MAX_GQ_PHRED,
    DEFAULT_MAX_ITERATIONS, DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT,
};
use crate::var_calling::variant_grouping::DEFAULT_MAX_VARIANT_GROUP_SPAN;
use crate::var_calling::{DEFAULT_MIN_ALT_OBS_PER_SAMPLE, DEFAULT_MIN_QUAL_PHRED};
use crate::vcf::DEFAULT_EMIT_GP;

/// Stage 1 knobs (CRAM-input filters, BAQ HMM, pileup walker).
/// Flattened into [`PileupArgs`](super::super::cli::PileupArgs).
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
/// [`VarCallingArgs`](super::super::var_calling::VarCallingArgs).
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

    /// Absolute ceiling on the unified allele set per variant group
    /// (REF + compounds + prunable), enforced *after*
    /// `--max-alleles-per-var`. Backstops the per-genotype likelihood
    /// routine's stack-allocated scratch + u64 genotype-membership
    /// bitmask. Normally inert; only binds in pathological
    /// low-complexity regions where chain-anchored compound alleles
    /// pile up. When it bites, the lowest-cohort_count alleles are
    /// dropped (REF excluded, compounds included) and the run summary
    /// reports `allele_lh_cap_hit: groups=N alleles_dropped=M`.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_ALLELES_LH_CALC,
        value_parser = parsers::parse_max_alleles_lh_calc,
        help_heading = "Advanced — Per-group merger",
    )]
    pub max_alleles_lh_calc: usize,

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

    /// Skip the allele-balance drop entirely (max sensitivity). The
    /// per-call allele balance is otherwise tested for biallelic het
    /// calls: a site is dropped only when every variant-carrying
    /// sample is a het whose observed ALT read fraction is
    /// inconsistent with the expected ~0.5 (the persistent-low-VAF
    /// artefact fingerprint). On by default.
    #[arg(
        long = "no-allele-balance-filter",
        hide_short_help = true,
        default_value_t = false,
        help_heading = "Advanced — Allele-balance filter"
    )]
    pub no_allele_balance_filter: bool,

    /// Drop threshold on the allele-balance log-likelihood ratio
    /// (Beta-Binomial fit of the genotype's expected balance vs a
    /// free fit). More negative = more inconsistent with a real
    /// genotype. The score is depth-aware, so a single threshold is
    /// self-silencing at low coverage. Default `-5.0` removes ~96 %
    /// of false positives at 300× for ~0.8 % true-call loss on the
    /// GIAB per_sample benchmark, and ~nothing below ~30×.
    #[arg(
        long = "min-allele-balance-log-lr",
        hide_short_help = true,
        default_value_t = crate::var_calling::allele_balance::DEFAULT_AB_MIN_LOG_LR,
        help_heading = "Advanced — Allele-balance filter",
    )]
    pub min_allele_balance_log_lr: f64,

    /// Beta-Binomial concentration (overdispersion) for the
    /// allele-balance test. Larger = tighter around the expected
    /// balance. This is the depth-honesty knob; the default was fit
    /// by MLE on high-coverage true hets.
    #[arg(
        long = "allele-balance-concentration",
        hide_short_help = true,
        default_value_t = crate::var_calling::allele_balance::DEFAULT_AB_CONCENTRATION,
        help_heading = "Advanced — Allele-balance filter",
    )]
    pub allele_balance_concentration: f64,

    /// Target false-discovery rate for the hidden-paralog filter, which
    /// drops candidate SNPs that a coverage + allele-balance model
    /// explains better as a reference-collapsed hidden paralog than as a
    /// real single-copy variant. On by default at a low (introgression-
    /// safe) rate: a real introgressed variant has normal coverage, so it
    /// never approaches the cut. `0.0` disables the filter (equivalent to
    /// `--no-paralog-filter`). See also `--no-paralog-filter`.
    #[arg(
        long = "paralog-fdr",
        hide_short_help = true,
        default_value_t = crate::var_calling::paralog_filter::DEFAULT_PARALOG_FDR,
        value_parser = parsers::parse_paralog_fdr,
        help_heading = "Advanced — Hidden-paralog filter",
    )]
    pub paralog_fdr: f64,

    /// Skip the hidden-paralog filter entirely, restoring the single-pass,
    /// direct-to-VCF behaviour (no spill, no coverage calibration). On by
    /// default the filter runs; use this (or `--paralog-fdr 0`) to turn it
    /// off.
    #[arg(
        long = "no-paralog-filter",
        hide_short_help = true,
        default_value_t = false,
        help_heading = "Advanced — Hidden-paralog filter"
    )]
    pub no_paralog_filter: bool,

    /// Run the hidden-paralog filter but do NOT drop the flagged loci: keep
    /// every scored locus and annotate it with the `PARALOG_POST` INFO field
    /// (the collapsed-paralog posterior). Use this to inspect the artifact
    /// class instead of removing it. No effect when the filter is off
    /// (`--no-paralog-filter` / `--paralog-fdr 0`).
    #[arg(
        long = "do-not-drop-dup-artifacts",
        hide_short_help = true,
        default_value_t = false,
        help_heading = "Advanced — Hidden-paralog filter"
    )]
    pub do_not_drop_dup_artifacts: bool,

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
