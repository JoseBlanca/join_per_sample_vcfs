//! Stage 6 — posterior engine.
//!
//! Consumes `Result<MergedRecord, PerGroupMergerError>` items from
//! Stage 5 and emits one [`PosteriorRecord`] per upstream record: a
//! per-sample, per-genotype posterior table plus the site-level
//! quantities (`p̂`, `f̂_C`, QUAL) and per-sample summaries (best
//! genotype, GQ) that the VCF writer needs. See
//! `doc/devel/implementation_plans/posterior_engine.md` and the
//! algorithmic contract in
//! `doc/devel/specs/calling_pipeline_architecture.md` §"Stage 6 —
//! posterior engine".
//!
//! Each upstream record is processed by an independent per-record EM
//! loop:
//!
//! 1. Flat `p̂[a] = 1 / n_alleles` and `f̂_C = 1 / n_alleles` initial
//!    estimates.
//! 2. E-step combines `MergedRecord.log_likelihoods` with the
//!    HWE-with-`F` prior (Wright–Fisher partition formula at
//!    arbitrary ploidy) into per-sample, per-genotype posteriors.
//! 3. M-step on `p̂` from posterior-weighted allele counts plus
//!    Dirichlet pseudocounts, and on `f̂_C` from posterior-weighted
//!    compound counts plus a Beta(α_compound, 1 − α_compound) prior.
//! 4. Convergence check on `max |Δp̂|`; hard cap at
//!    `config.max_iterations` returns
//!    [`PosteriorEngineError::DidNotConverge`].
//!
//! Compound alleles get the cohort-derived `f̂_C` substituted for
//! `p̂[compound]` everywhere the prior is evaluated, regardless of
//! whether the sample is chain-evident or chain-broken. The
//! chain-broken flag in [`MergedRecord::chain_anchor_flag`] only
//! affects the likelihood path, which Stage 5 has already baked into
//! `log_likelihoods`.
//!
//! Contamination correction is opt-in via
//! [`PosteriorEngineConfig::contamination`]: `Some(_)` activates the
//! mixture-likelihood E-step using the frozen `c_s` / `q_b` produced
//! by the side-pass at
//! [`crate::var_calling::contamination_estimation`]; `None` runs the
//! `c_s = 0` path with no mixture overhead. The
//! `--approximate-posterior-calculation` flag exists in the config so
//! the CLI surface can land independently, but the precomputed-LUT
//! implementation is deferred; setting the flag to `true` returns
//! [`PosteriorEngineError::ApproximateModeNotYetImplemented`].

pub mod backends;
mod interp;
mod shape;

use std::fmt;

use thiserror::Error;

use self::backends::{InterpUnivariateMath, MathBackend};
use self::shape::{GenotypeShape, shape_for};
use crate::per_sample_pileup::pileup::AlleleSupportStats;
use crate::var_calling::contamination_estimation::{
    AlleleClass as ContamAlleleClass, ContaminationEstimates, MAX_BASE_ERROR, MIN_BASE_ERROR,
};
use crate::var_calling::per_group_merger::{
    MergedAllele, MergedRecord, PerGroupMergerError, genotype_order,
};

/// Default EM convergence threshold on `max_a |p̂_new[a] − p̂_old[a]|`.
/// Inherited from the implementation plan; revisit against real data.
///
/// Source: implementation-plan `posterior_engine.md` §"Step 5 —
/// convergence check" picks `1e-4` as the per-record threshold; GATK's
/// `0.1` on raw allele counts is roughly equivalent at depth ≈ 20.
pub const DEFAULT_CONVERGENCE_THRESHOLD: f64 = 1e-4;

/// Default hard cap on per-record EM iterations. Exceeding it returns
/// [`PosteriorEngineError::DidNotConverge`]. Empirically the
/// closed-form M-steps converge in 3–5 rounds per the GATK reference;
/// 50 is the belt-and-braces ceiling.
///
/// Source: implementation-plan `posterior_engine.md` §"Step 5 —
/// convergence check".
pub const DEFAULT_MAX_ITERATIONS: u32 = 50;

/// Dirichlet pseudocount on the REF allele. GATK-style. Equivalent to
/// a pre-observation belief that REF is overwhelmingly the common
/// allele at a typical site.
///
/// Source: implementation-plan `posterior_engine.md` §"Configurable
/// parameters" and the architecture spec's prior shape. Mirrors the
/// `α_ref = 10`, `α_alt = 0.01` GATK-shaped Dirichlet recorded at
/// `calling_pipeline_architecture.md` §"From likelihood to posterior".
/// Revisit against the cohort calibration set.
pub const DEFAULT_REF_PSEUDOCOUNT: f64 = 10.0;

/// Dirichlet pseudocount on a non-compound SNP/MNP alt allele.
///
/// Source: matches the GATK-style `α_alt = 0.01` for SNP/MNP alts
/// recorded at `calling_pipeline_architecture.md` §"From likelihood
/// to posterior". Revisit against the cohort calibration set.
pub const DEFAULT_SNP_ALT_PSEUDOCOUNT: f64 = 0.01;

/// Dirichlet pseudocount on a non-compound indel alt allele.
///
/// Source: matches the GATK-style `α_alt = 0.00125` for indels
/// recorded at `calling_pipeline_architecture.md` §"From likelihood
/// to posterior". The `0.01 : 0.00125 = 8 : 1` SNP-to-indel ratio
/// mirrors GATK's per-class default. Revisit against the cohort
/// calibration set.
pub const DEFAULT_INDEL_ALT_PSEUDOCOUNT: f64 = 0.00125;

/// Beta-`α` pseudocount on a chain-anchored compound allele for the
/// per-compound `f̂_C` estimator. The Beta partner is implicitly
/// `1 − α_compound`, placing the prior mean at `α_compound` exactly.
///
/// Source: implementation-plan `posterior_engine.md` §"Step 4 —
/// M-step on `f̂_C`". The Beta-partner choice is documented in the
/// 2026-05-16 implementation report. Revisit against the cohort
/// calibration set.
pub const DEFAULT_COMPOUND_ALT_PSEUDOCOUNT: f64 = 0.001;

/// Default per-sample inbreeding coefficient. `0.0` recovers strict
/// HWE; `1.0` allows only homozygous genotypes. The project default
/// is outcrossing-friendly even though the plant-breeding cohort
/// target tends to want a non-zero override per the spec.
pub const DEFAULT_INBREEDING_COEFFICIENT: f64 = 0.0;

/// Phred cap on per-sample GQ. Matches the GATK and bcftools
/// convention (`GQ` capped at 99) so downstream tooling sees a
/// familiar range; the cap also prevents `+∞` GQ when EM yields
/// `P(best) = 1` exactly. The CLI exposes this as `--max-gq-phred`
/// once the cohort subcommand lands.
pub const MAX_GQ_PHRED: f64 = 99.0;

/// Phred conversion factor — `Phred = -10 · log10(p)`. The literal
/// `-10.0` appears in every Phred-conversion site; the named const
/// keeps the meaning in one place.
const PHRED_SCALE: f64 = -10.0;

/// Index of the hom-REF genotype in canonical [`genotype_order`].
/// REF/REF/.../REF is the all-zero allele tuple and lands at position
/// 0 for every (ploidy, n_alleles).
const HOM_REF_GENOTYPE_IDX: usize = 0;

/// Tunable knobs for the posterior engine.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct PosteriorEngineConfig {
    /// Convergence threshold on `max_a |p̂_new[a] − p̂_old[a]|`.
    /// Iteration stops as soon as the latest M-step on `p̂` moves all
    /// alleles by less than this. Defaults to
    /// [`DEFAULT_CONVERGENCE_THRESHOLD`].
    pub convergence_threshold: f64,
    /// Hard cap on per-record EM iterations. Defaults to
    /// [`DEFAULT_MAX_ITERATIONS`].
    pub max_iterations: u32,
    /// Dirichlet pseudocount for the REF allele. Defaults to
    /// [`DEFAULT_REF_PSEUDOCOUNT`].
    pub ref_pseudocount: f64,
    /// Dirichlet pseudocount for a non-compound SNP/MNP alt allele.
    /// Defaults to [`DEFAULT_SNP_ALT_PSEUDOCOUNT`].
    pub snp_alt_pseudocount: f64,
    /// Dirichlet pseudocount for a non-compound indel alt allele.
    /// Defaults to [`DEFAULT_INDEL_ALT_PSEUDOCOUNT`].
    pub indel_alt_pseudocount: f64,
    /// Beta-`α` for the per-compound `f̂_C` estimator. The Beta
    /// partner is implicitly `1 − α_compound` so the prior mean is
    /// `α_compound`. Defaults to [`DEFAULT_COMPOUND_ALT_PSEUDOCOUNT`].
    pub compound_alt_pseudocount: f64,
    /// Default per-sample inbreeding coefficient `F`. Used for every
    /// sample not covered by `fixation_index_overrides`. Defaults to
    /// [`DEFAULT_INBREEDING_COEFFICIENT`].
    pub fixation_index_default: f64,
    /// Per-sample inbreeding overrides. If `Some`, length equals the
    /// cohort size and entries replace `fixation_index_default`. Empty
    /// in v1; the CLI exposes a scalar today, and per-sample F via a
    /// file is a follow-up that flips this from `None` to `Some(_)`
    /// without an engine-side refactor.
    pub fixation_index_overrides: Option<Vec<f64>>,
    /// GQ Phred cap. Defaults to [`MAX_GQ_PHRED`].
    pub max_gq_phred: f64,
    /// Opt-in flag for the precomputed-LUT inner loop (multinomial
    /// coefficients, HWE-with-`F` quantised priors). **The LUT
    /// implementation is not yet in tree** — setting this `true`
    /// returns [`PosteriorEngineError::ApproximateModeNotYetImplemented`].
    /// The field exists so the CLI surface can be wired before the
    /// implementation lands; see the plan §"Approximation via
    /// precomputed lookup tables".
    pub approximate_posterior_calculation: bool,
    /// Frozen contamination estimates produced by the side-pass
    /// ([`estimate_contamination`]) or supplied by the user. `None`
    /// disables contamination correction entirely and the per-record
    /// EM uses [`MergedRecord::log_likelihoods`] verbatim. `Some(_)`
    /// activates the mixture-likelihood E-step using the carried
    /// `c_s` and `q_b` values.
    ///
    /// [`estimate_contamination`]: crate::var_calling::contamination_estimation::estimate_contamination
    pub contamination: Option<ContaminationEstimates>,
}

impl Default for PosteriorEngineConfig {
    fn default() -> Self {
        Self::with_project_defaults()
    }
}

impl PosteriorEngineConfig {
    /// Construct the engine config with the project defaults.
    ///
    /// # Defaults
    ///
    /// | Field                              | Constant                              | Value    |
    /// |---|---|---|
    /// | `convergence_threshold`            | [`DEFAULT_CONVERGENCE_THRESHOLD`]     | 1e-4     |
    /// | `max_iterations`                   | [`DEFAULT_MAX_ITERATIONS`]            | 50       |
    /// | `ref_pseudocount`                  | [`DEFAULT_REF_PSEUDOCOUNT`]           | 10.0     |
    /// | `snp_alt_pseudocount`              | [`DEFAULT_SNP_ALT_PSEUDOCOUNT`]       | 0.01     |
    /// | `indel_alt_pseudocount`            | [`DEFAULT_INDEL_ALT_PSEUDOCOUNT`]     | 0.00125  |
    /// | `compound_alt_pseudocount`         | [`DEFAULT_COMPOUND_ALT_PSEUDOCOUNT`]  | 0.001    |
    /// | `fixation_index_default`           | [`DEFAULT_INBREEDING_COEFFICIENT`]    | 0.0      |
    /// | `fixation_index_overrides`         | —                                     | `None`   |
    /// | `max_gq_phred`                     | [`MAX_GQ_PHRED`]                      | 99.0     |
    /// | `approximate_posterior_calculation`| —                                     | `false`  |
    /// | `contamination`                    | —                                     | `None`   |
    pub fn with_project_defaults() -> Self {
        Self {
            convergence_threshold: DEFAULT_CONVERGENCE_THRESHOLD,
            max_iterations: DEFAULT_MAX_ITERATIONS,
            ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT,
            snp_alt_pseudocount: DEFAULT_SNP_ALT_PSEUDOCOUNT,
            indel_alt_pseudocount: DEFAULT_INDEL_ALT_PSEUDOCOUNT,
            compound_alt_pseudocount: DEFAULT_COMPOUND_ALT_PSEUDOCOUNT,
            fixation_index_default: DEFAULT_INBREEDING_COEFFICIENT,
            fixation_index_overrides: None,
            max_gq_phred: MAX_GQ_PHRED,
            approximate_posterior_calculation: false,
            contamination: None,
        }
    }
}

/// Why a posterior calculation hit an unrecoverable value. Surfaced
/// through [`PosteriorEngineError::NonFinitePosterior`] because the
/// closed-form EM is finite for every valid `MergedRecord`; surfacing
/// signals an internal bug, not a data condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum NonFiniteKind {
    NaN,
    PositiveInfinity,
    NegativeInfinity,
}

impl fmt::Display for NonFiniteKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            NonFiniteKind::NaN => "NaN",
            NonFiniteKind::PositiveInfinity => "+infinity",
            NonFiniteKind::NegativeInfinity => "-infinity",
        };
        f.write_str(s)
    }
}

/// Reason a `MergedRecord` failed an internal-shape invariant. The
/// posterior engine trusts Stage 5 to produce well-shaped records;
/// this enum names what was wrong when the trust is misplaced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MalformedMergedRecordReason {
    /// `n_genotypes` disagrees with `genotype_order(ploidy, n_alleles).len()`.
    GenotypeCountMismatch { expected: usize, got: usize },
    /// `log_likelihoods.len()` is not `n_samples * n_genotypes`.
    LogLikelihoodLengthMismatch { expected: usize, got: usize },
    /// `chain_anchor_flags.len()` is not `n_samples * n_alleles`.
    ChainAnchorLengthMismatch { expected: usize, got: usize },
}

impl fmt::Display for MalformedMergedRecordReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MalformedMergedRecordReason::GenotypeCountMismatch { expected, got } => write!(
                f,
                "n_genotypes mismatch: expected {expected} per `genotype_order`, got {got}"
            ),
            MalformedMergedRecordReason::LogLikelihoodLengthMismatch { expected, got } => write!(
                f,
                "log_likelihoods length mismatch: expected {expected} (n_samples × n_genotypes), got {got}"
            ),
            MalformedMergedRecordReason::ChainAnchorLengthMismatch { expected, got } => write!(
                f,
                "chain_anchor_flags length mismatch: expected {expected} (n_samples × n_alleles), got {got}"
            ),
        }
    }
}

/// Genomic-position triple — chrom_id + 1-based inclusive `[start, end]`
/// span — used to annotate engine errors with where they happened.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RecordLocus {
    pub chrom_id: u32,
    pub start: u32,
    pub end: u32,
}

impl fmt::Display for RecordLocus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "chrom {} {}-{}", self.chrom_id, self.start, self.end)
    }
}

/// Errors surfaced by the posterior engine.
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum PosteriorEngineError {
    #[error("failed to consume upstream merged record")]
    Upstream(#[from] PerGroupMergerError),

    /// EM produced a non-finite per-sample posterior. The closed-form
    /// math is finite for every valid `MergedRecord`; surfacing this
    /// signals an internal bug, not a data condition.
    ///
    /// `genotype_idx` is `None` when the non-finite value is the
    /// per-sample normaliser `log_z` (no single genotype to blame);
    /// `Some(g)` when the value is the per-genotype posterior at
    /// index `g`.
    #[error(
        "non-finite posterior at {locus} for sample_idx {sample_idx} \
         (genotype_idx {genotype_idx:?}): {kind}"
    )]
    NonFinitePosterior {
        locus: RecordLocus,
        sample_idx: usize,
        genotype_idx: Option<usize>,
        kind: NonFiniteKind,
    },

    /// EM did not converge within `max_iterations`. A hard cap exists
    /// so that pathological records cannot stall the engine.
    #[error(
        "EM did not converge within {max_iterations} iterations at {locus} \
         (last delta: {last_delta})"
    )]
    DidNotConverge {
        locus: RecordLocus,
        max_iterations: u32,
        last_delta: f64,
    },

    /// `fixation_index_overrides` was supplied but its length does
    /// not match the upstream record's sample count.
    #[error(
        "fixation_index_overrides length {override_len} does not match \
         record sample count {record_samples} at {locus}"
    )]
    InbreedingOverrideLengthMismatch {
        locus: RecordLocus,
        override_len: usize,
        record_samples: usize,
    },

    /// Upstream record violates an internal shape invariant
    /// (`n_genotypes` / `log_likelihoods` / `chain_anchor_flags`
    /// dimensions). Stage 5 is supposed to maintain these; surfacing
    /// signals an upstream bug rather than a data condition.
    #[error("malformed merged record at {locus}: {reason}")]
    MalformedMergedRecord {
        locus: RecordLocus,
        reason: MalformedMergedRecordReason,
    },

    /// Approximate-posterior mode is configured but the LUT
    /// implementation is not yet in tree.
    #[error("approximate-posterior calculation is not yet implemented")]
    ApproximateModeNotYetImplemented,

    /// `PosteriorEngineConfig::contamination.sample_to_batch.len()`
    /// did not match the upstream record's `n_samples`. The frozen
    /// contamination estimates were produced for a different cohort
    /// than the one driving Stage 5.
    #[error(
        "contamination estimate cohort size {estimates_n_samples} does not match \
         record sample count {record_samples} at {locus}"
    )]
    ContaminationCohortSizeMismatch {
        locus: RecordLocus,
        estimates_n_samples: usize,
        record_samples: usize,
    },
}

/// One emitted item: the per-sample, per-genotype posterior for one
/// merged record plus the site-level frequencies, QUAL, and the
/// forwarded scalar table the VCF writer needs.
///
/// **Storage layout.** `posteriors`, `scalars`, and
/// `chain_anchor_flags` are row-major flat `Vec`s (sample-major), the
/// same layout as [`MergedRecord`]. Use [`Self::posteriors_row`],
/// [`Self::scalars_row`], and [`Self::chain_anchor_flags_row`] to
/// read them. The flat layout collapses what would otherwise be one
/// allocation per sample per row into one allocation per row, in
/// keeping with the 2026-05-16 perf review of Stage 5.
#[derive(Debug, Clone, PartialEq)]
pub struct PosteriorRecord {
    pub locus: RecordLocus,
    /// Merged allele set, forwarded from the upstream `MergedRecord`.
    /// `alleles[0]` is always REF.
    pub alleles: Vec<MergedAllele>,
    /// Cohort-wide ploidy, forwarded from the upstream `MergedRecord`.
    pub ploidy: u8,
    /// Number of samples this record carries posteriors for.
    pub n_samples: usize,
    /// Number of genotypes per sample in `posteriors`. Equals
    /// `genotype_order(ploidy, alleles.len()).len()`.
    pub n_genotypes: usize,
    /// Per-allele frequency estimate `p̂` after EM convergence.
    /// Length equals `alleles.len()`; sums to 1.
    pub allele_frequencies: Vec<f64>,
    /// Per-compound cohort frequency estimate `f̂_C`. `Some(f)` for
    /// compound alleles, `None` otherwise. Length equals
    /// `alleles.len()`.
    pub compound_frequencies: Vec<Option<f64>>,
    /// Flat row-major posterior table —
    /// `posteriors[sample_idx * n_genotypes + genotype_idx]`. Genotype
    /// order is canonical per [`genotype_order`]. Each per-sample row
    /// sums to 1 within floating-point tolerance.
    pub posteriors: Vec<f64>,
    /// Per-sample best genotype (argmax of `posteriors[sample_idx]`).
    /// Indexed by genotype-enumeration position.
    pub best_genotype: Vec<usize>,
    /// Per-sample genotype quality, Phred of `1 − P(best)`. Clamped
    /// at [`PosteriorEngineConfig::max_gq_phred`].
    pub gq_phred: Vec<f64>,
    /// Site-level QUAL: Phred of `Π_s P(hom-ref)_s`. Passed through
    /// as `f64::INFINITY` when at least one sample is certainly
    /// variant; the VCF writer caps if needed.
    pub qual_phred: f64,
    /// Forwarded from `MergedRecord.scalars`, same row-major layout.
    /// The VCF writer uses these for DP/AD-style FORMAT fields
    /// without a parallel join.
    pub scalars: Vec<AlleleSupportStats>,
    /// Forwarded from `MergedRecord.other_scalars`.
    pub other_scalars: Vec<AlleleSupportStats>,
    /// Forwarded from `MergedRecord.chain_anchor_flags`, same
    /// row-major layout.
    pub chain_anchor_flags: Vec<bool>,
    /// EM bookkeeping — iteration count, final delta.
    pub diagnostics: EmDiagnostics,
}

impl PosteriorRecord {
    /// Per-sample slice into the flat row-major `posteriors` table.
    ///
    /// # Panics
    ///
    /// Panics if `sample_idx >= self.n_samples` (slice out of bounds).
    pub fn posteriors_row(&self, sample_idx: usize) -> &[f64] {
        let start = sample_idx * self.n_genotypes;
        &self.posteriors[start..start + self.n_genotypes]
    }

    /// Per-sample slice into the flat row-major `scalars` table.
    ///
    /// # Panics
    ///
    /// Panics if `sample_idx >= self.n_samples` (slice out of bounds).
    pub fn scalars_row(&self, sample_idx: usize) -> &[AlleleSupportStats] {
        let n_alleles = self.alleles.len();
        let start = sample_idx * n_alleles;
        &self.scalars[start..start + n_alleles]
    }

    /// Per-sample slice into the flat row-major `chain_anchor_flags` table.
    ///
    /// # Panics
    ///
    /// Panics if `sample_idx >= self.n_samples` (slice out of bounds).
    pub fn chain_anchor_flags_row(&self, sample_idx: usize) -> &[bool] {
        let n_alleles = self.alleles.len();
        let start = sample_idx * n_alleles;
        &self.chain_anchor_flags[start..start + n_alleles]
    }
}

/// EM bookkeeping for a single record. A `PosteriorRecord` is only
/// emitted on convergence, so an emitted record's
/// `iterations <= config.max_iterations` and
/// `final_max_delta_p < config.convergence_threshold` always hold.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EmDiagnostics {
    /// Iterations actually run (≤ `config.max_iterations`).
    pub iterations: u32,
    /// `max_a |p̂_new[a] − p̂_old[a]|` at the last iteration.
    pub final_max_delta_p: f64,
}

/// Streaming posterior engine. Wraps an upstream `MergedRecord`
/// iterator; emits one [`PosteriorRecord`] per upstream record.
///
/// Errors latch: once `next()` surfaces any error variant, subsequent
/// calls return `None`. Matches the Stage 4 / 5 merger precedent.
///
/// # Errors
///
/// The `Iterator` impl can yield any [`PosteriorEngineError`] variant:
/// `Upstream` from Stage 5, `NonFinitePosterior` / `DidNotConverge` /
/// `MalformedMergedRecord` from per-record EM,
/// `InbreedingOverrideLengthMismatch` from
/// [`PosteriorEngineConfig::fixation_index_overrides`], and
/// `ApproximateModeNotYetImplemented` for the deferred LUT mode,
/// and `ContaminationCohortSizeMismatch` when
/// [`PosteriorEngineConfig::contamination`]'s cohort size disagrees
/// with the upstream record's `n_samples`.
///
/// # Examples
///
/// ```
/// use merge_per_sample_vcfs::var_calling::per_group_merger::{
///     MergedRecord, PerGroupMergerError,
/// };
/// use merge_per_sample_vcfs::var_calling::posterior_engine::{
///     PosteriorEngine, PosteriorEngineError, PosteriorRecord,
/// };
///
/// fn run_stage6<I>(upstream: I) -> Result<Vec<PosteriorRecord>, PosteriorEngineError>
/// where
///     I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
/// {
///     PosteriorEngine::new(upstream).collect()
/// }
/// ```
pub struct PosteriorEngine<I, M = InterpUnivariateMath>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
    M: MathBackend,
{
    upstream: I,
    config: PosteriorEngineConfig,
    math: M,
    is_latched: bool,
}

impl<I, M> fmt::Debug for PosteriorEngine<I, M>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
    M: MathBackend,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Exhaustive destructure: a new field on `PosteriorEngine`
        // fails to compile here so the omission from Debug output is
        // explicit rather than silent.
        let Self {
            upstream: _,
            config,
            math: _,
            is_latched,
        } = self;
        f.debug_struct("PosteriorEngine")
            .field("config", config)
            .field("is_latched", is_latched)
            .finish()
    }
}

impl<I> PosteriorEngine<I, InterpUnivariateMath>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
{
    /// Construct an engine with the [project defaults] and the default
    /// math backend ([`InterpUnivariateMath`], ~10 % faster than
    /// [`ExactMath`] on the contam-on bench with ~50× margin under the
    /// approximate-parity budget — see the implementation plan's
    /// step-7 decision).
    ///
    /// For bit-identical reproducibility against the unoptimised
    /// engine, construct with
    /// [`with_math_backend`](Self::with_math_backend) and pass
    /// [`ExactMath`] explicitly.
    ///
    /// [project defaults]: PosteriorEngineConfig::with_project_defaults
    pub fn new(upstream: I) -> Self {
        Self::with_config(upstream, PosteriorEngineConfig::with_project_defaults())
    }

    /// Construct an engine with explicit tuning and the default math
    /// backend ([`InterpUnivariateMath`]).
    pub fn with_config(upstream: I, config: PosteriorEngineConfig) -> Self {
        Self::with_math_backend(upstream, config, InterpUnivariateMath)
    }
}

impl<I, M> PosteriorEngine<I, M>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
    M: MathBackend,
{
    /// Construct an engine with explicit tuning and an explicit math
    /// backend. Use this to opt into an approximate backend (e.g.
    /// `InterpUnivariateMath`) or to opt back into [`ExactMath`] if a
    /// future release flips the default.
    pub fn with_math_backend(upstream: I, config: PosteriorEngineConfig, math: M) -> Self {
        Self {
            upstream,
            config,
            math,
            is_latched: false,
        }
    }

    /// Read-only view of the tuning this engine was constructed with.
    /// Stable across `next()` calls.
    pub fn config(&self) -> &PosteriorEngineConfig {
        &self.config
    }

    fn process(&mut self, record: MergedRecord) -> Result<PosteriorRecord, PosteriorEngineError> {
        if self.config.approximate_posterior_calculation {
            return Err(PosteriorEngineError::ApproximateModeNotYetImplemented);
        }
        run_em_for_record(record, &self.config, &self.math)
    }
}

impl<I, M> Iterator for PosteriorEngine<I, M>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
    M: MathBackend,
{
    type Item = Result<PosteriorRecord, PosteriorEngineError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_latched {
            return None;
        }
        match self.upstream.next() {
            None => {
                self.is_latched = true;
                None
            }
            Some(Err(e)) => {
                self.is_latched = true;
                Some(Err(e.into()))
            }
            Some(Ok(record)) => match self.process(record) {
                Ok(pr) => Some(Ok(pr)),
                Err(e) => {
                    self.is_latched = true;
                    Some(Err(e))
                }
            },
        }
    }
}

// ---------------------------------------------------------------------
// Per-record EM
// ---------------------------------------------------------------------

/// Allele class for Dirichlet-pseudocount routing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AlleleClass {
    Ref,
    SnpAlt,
    IndelAlt,
    CompoundAlt,
}

fn classify_alleles(alleles: &[MergedAllele]) -> Vec<AlleleClass> {
    let ref_len = alleles[0].seq.len();
    alleles
        .iter()
        .enumerate()
        .map(|(idx, allele)| {
            if idx == 0 {
                AlleleClass::Ref
            } else if allele.is_compound {
                AlleleClass::CompoundAlt
            } else if allele.seq.len() == ref_len {
                AlleleClass::SnpAlt
            } else {
                AlleleClass::IndelAlt
            }
        })
        .collect()
}

/// Map Stage 6's internal [`AlleleClass`] onto the contamination
/// side-pass's 3-variant [`ContamAlleleClass`]. Compound alleles do
/// not appear in the side-pass's `q_b` (the side-pass operates on raw
/// per-position observations), so v1 maps them to `None` — the
/// mixture-likelihood term for compound classes is then `q_b = 0`,
/// which treats contamination at compound alleles as impossible. This
/// is conservative; revisit if real data shows the assumption matters.
fn map_to_contam_class(class: AlleleClass) -> Option<ContamAlleleClass> {
    match class {
        AlleleClass::Ref => Some(ContamAlleleClass::Ref),
        AlleleClass::SnpAlt => Some(ContamAlleleClass::SnpAlt),
        AlleleClass::IndelAlt => Some(ContamAlleleClass::IndelAlt),
        AlleleClass::CompoundAlt => None,
    }
}

/// Mixture-likelihood reconstruction from Stage 5's per-(sample,
/// allele) scalars and the frozen contamination parameters. Returns a
/// fresh `n_samples × n_genotypes` row-major `log_likelihoods` table
/// suitable for the existing EM loop.
///
/// Per-genotype likelihood: `Σ_a n_{s,a} · log[(1 − c_s) · P(read with
/// allele a | G) + c_s · q_{b(s)}[class(a)]]`, where `P(read | G)` for
/// a genotype with allele copy counts `k_a` is `(k_a / ploidy) · (1 −
/// ε_{s,a}) + ((ploidy − k_a) / ploidy) · (ε_{s,a} / (n_alleles − 1))`.
/// Per-allele `ε_{s,a}` is derived from the `q_sum` aggregate.
///
/// # Panics
///
/// Debug-only invariants:
/// - `scalars.len() == n_samples * n_alleles`
/// - `fallback_log_likelihoods.len() == n_samples * n_genotypes`
///
/// In release builds these are caller responsibilities; the only
/// in-tree caller is [`run_em_for_record`], which validates the
/// merged record's shape via [`validate_record_shape`] before calling
/// here.
fn compute_mixture_log_likelihoods<M: MathBackend>(
    locus: RecordLocus,
    alleles: &[MergedAllele],
    ploidy: u8,
    scalars: &[AlleleSupportStats],
    estimates: &ContaminationEstimates,
    fallback_log_likelihoods: &[f64],
    math: &M,
) -> Result<Vec<f64>, PosteriorEngineError> {
    let n_alleles = alleles.len();
    let shape = shape_for(ploidy, n_alleles);
    let n_genotypes = shape.n_genotypes;
    let n_samples = scalars.len() / n_alleles;

    debug_assert_eq!(
        scalars.len(),
        n_samples * n_alleles,
        "scalars length does not match n_samples * n_alleles"
    );
    debug_assert_eq!(
        fallback_log_likelihoods.len(),
        n_samples * n_genotypes,
        "fallback_log_likelihoods length does not match n_samples * n_genotypes"
    );

    let allele_classes = classify_alleles(alleles);
    let contam_class_per_allele: Vec<Option<ContamAlleleClass>> = allele_classes
        .iter()
        .copied()
        .map(map_to_contam_class)
        .collect();
    let genotype_allele_counts = &shape.genotype_allele_counts;

    let ploidy_f = f64::from(ploidy);
    // Denominator for spreading per-base error mass across alleles
    // other than the one the genotype carries; floored at 1 so the
    // n_alleles == 1 degenerate case doesn't divide by zero.
    let other_allele_error_denom = (n_alleles as f64 - 1.0).max(1.0);

    let mut mixture_log_likelihoods = vec![0.0_f64; n_samples * n_genotypes];

    for sample_idx in 0..n_samples {
        let c_s = estimates.effective_c_s(sample_idx);
        if c_s <= 0.0 {
            // Sample is below-floor / singleton-batch: c_s = 0, the
            // mixture collapses to the own-DNA term. Fall back to
            // Stage 5's precomputed value to avoid re-deriving the
            // same number with floating-point drift.
            let dst = &mut mixture_log_likelihoods
                [sample_idx * n_genotypes..(sample_idx + 1) * n_genotypes];
            let src =
                &fallback_log_likelihoods[sample_idx * n_genotypes..(sample_idx + 1) * n_genotypes];
            dst.copy_from_slice(src);
            continue;
        }
        let q_b = estimates.q_b_for_sample(sample_idx);

        // Per-allele observation count and mean per-read error rate
        // for this sample at this site, recovered from the row of
        // AlleleSupportStats Stage 5 emitted.
        let sample_allele_scalars = &scalars[sample_idx * n_alleles..(sample_idx + 1) * n_alleles];
        let mut n_obs = vec![0_u32; n_alleles];
        let mut mean_err = vec![0.0_f64; n_alleles];
        for (a, support) in sample_allele_scalars.iter().enumerate() {
            n_obs[a] = support.num_obs;
            if support.num_obs > 0 {
                let per_read_log_err = support.q_sum / f64::from(support.num_obs);
                mean_err[a] = math.exp(per_read_log_err).clamp(MIN_BASE_ERROR, MAX_BASE_ERROR);
            }
        }

        for (g_idx, gt_counts) in genotype_allele_counts.chunks_exact(n_alleles).enumerate() {
            let mut ll: f64 = 0.0;
            for a in 0..n_alleles {
                let n_a = n_obs[a];
                if n_a == 0 {
                    continue;
                }
                let k_a = f64::from(gt_counts[a]);
                let eps = mean_err[a];
                // P(read carrying allele a | genotype G).
                let p_own = (k_a / ploidy_f) * (1.0 - eps)
                    + ((ploidy_f - k_a) / ploidy_f) * (eps / other_allele_error_denom);
                let p_contam = match contam_class_per_allele[a] {
                    Some(class) => q_b[class as usize],
                    None => 0.0, // compound allele class — see map_to_contam_class
                };
                let mix = (1.0 - c_s) * p_own + c_s * p_contam;
                if mix <= 0.0 {
                    // Defensive: this can only happen if both p_own
                    // and p_contam are 0, i.e. the read carries an
                    // allele the genotype could not produce and the
                    // contaminant cannot supply. Surface as a
                    // non-finite-posterior because the resulting log
                    // is −∞.
                    return Err(PosteriorEngineError::NonFinitePosterior {
                        locus,
                        sample_idx,
                        genotype_idx: Some(g_idx),
                        kind: NonFiniteKind::NegativeInfinity,
                    });
                }
                ll += f64::from(n_a) * math.ln(mix);
            }
            mixture_log_likelihoods[sample_idx * n_genotypes + g_idx] = ll;
        }
    }

    Ok(mixture_log_likelihoods)
}

fn pseudocount_for(class: AlleleClass, config: &PosteriorEngineConfig) -> f64 {
    match class {
        AlleleClass::Ref => config.ref_pseudocount,
        AlleleClass::SnpAlt => config.snp_alt_pseudocount,
        AlleleClass::IndelAlt => config.indel_alt_pseudocount,
        AlleleClass::CompoundAlt => config.compound_alt_pseudocount,
    }
}

/// Per-record EM inputs that don't change across iterations.
/// Bundled into one struct so `e_step` / `m_step` don't take 10+
/// parameters each.
struct EmContext<'a> {
    locus: RecordLocus,
    n_samples: usize,
    n_genotypes: usize,
    n_alleles: usize,
    ploidy: u8,
    /// Precomputed genotype-shape artefacts:
    /// `genotype_allele_counts`, `log_multinomial_coeffs`,
    /// `nonzero_pairs` / `nonzero_pairs_offsets`, and
    /// `homozygous_allele_for`. Shared across records of the same
    /// `(ploidy, n_alleles)` via the thread-local cache.
    shape: &'a GenotypeShape,
    compound_mask: &'a [bool],
    pseudocounts: &'a [f64],
    /// `safe_ln(f_s)` per sample, where `f_s` is the per-sample
    /// fixation index. Hoisted out of `e_step`'s per-iteration loop
    /// because `fixation_indices` is record-static.
    log_f_per_sample: &'a [f64],
    /// `safe_ln(1.0 - f_s)` per sample. Same record-static rationale as
    /// `log_f_per_sample`.
    log_one_minus_f_per_sample: &'a [f64],
    compound_pseudocount: f64,
}

/// Per-record EM scratch buffers, allocated once per record and
/// reused across iterations.
struct EmScratch {
    /// Effective per-allele frequency seen by the prior (compounds
    /// substituted with `f̂_C`). Length `n_alleles`.
    p_effective: Vec<f64>,
    /// Natural log of `p_effective`. Length `n_alleles`.
    log_p_effective: Vec<f64>,
    /// Per-genotype unnormalised log-posterior for the current
    /// sample. Length `n_genotypes`.
    log_post_unnorm: Vec<f64>,
    /// Posterior-weighted allele counts E[n_a]. Length `n_alleles`.
    expected_counts: Vec<f64>,
}

impl EmScratch {
    fn new(n_alleles: usize, n_genotypes: usize) -> Self {
        Self {
            p_effective: vec![0.0; n_alleles],
            log_p_effective: vec![0.0; n_alleles],
            log_post_unnorm: vec![0.0; n_genotypes],
            expected_counts: vec![0.0; n_alleles],
        }
    }
}

fn run_em_for_record<M: MathBackend>(
    record: MergedRecord,
    config: &PosteriorEngineConfig,
    math: &M,
) -> Result<PosteriorRecord, PosteriorEngineError> {
    let MergedRecord {
        chrom_id,
        start,
        end,
        alleles,
        ploidy,
        n_samples,
        n_genotypes,
        scalars,
        other_scalars,
        chain_anchor_flags,
        log_likelihoods,
    } = record;

    let locus = RecordLocus {
        chrom_id,
        start,
        end,
    };
    let n_alleles = alleles.len();

    // Mixture-likelihood branch: when contamination is configured,
    // replace the Stage-5-supplied (c_s = 0) log-likelihood table with
    // the mixture-recomputed version before the EM loop runs.
    // Everything downstream is unchanged.
    let log_likelihoods = if let Some(estimates) = config.contamination.as_ref() {
        if estimates.c_s_per_sample.len() != n_samples {
            return Err(PosteriorEngineError::ContaminationCohortSizeMismatch {
                locus,
                estimates_n_samples: estimates.c_s_per_sample.len(),
                record_samples: n_samples,
            });
        }
        // M8: the side-pass's fallible constructors
        // (`ContaminationEstimates::zero` / `::from_user_supplied`)
        // and the in-engine `finalise` are the only construction
        // sites; each enforces `sample_to_batch.len() ==
        // c_s_per_sample.len()` and `batch_idx < q_b_per_batch.len()`.
        // Together with the cohort-size check above, this makes
        // `effective_c_s` and `q_b_for_sample` indexing PANIC-FREE for
        // every well-formed `ContaminationEstimates`. Debug-only
        // sanity check.
        debug_assert_eq!(
            estimates.sample_to_batch.len(),
            n_samples,
            "ContaminationEstimates internal invariant violated"
        );
        debug_assert!(
            estimates
                .sample_to_batch
                .iter()
                .all(|&b| b < estimates.q_b_per_batch.len()),
            "ContaminationEstimates internal invariant violated"
        );
        if n_alleles >= 2 {
            compute_mixture_log_likelihoods(
                locus,
                &alleles,
                ploidy,
                &scalars,
                estimates,
                &log_likelihoods,
                math,
            )?
        } else {
            log_likelihoods
        }
    } else {
        log_likelihoods
    };

    // Belt-and-braces: Stage 5 drops REF-only groups (alleles.len() <
    // 2). If one slips through, emit a trivial posterior record per
    // the plan §"Edge cases" — every sample is hom-REF with QUAL 0.
    if n_alleles < 2 {
        return Ok(trivial_posterior_record(
            TrivialRecordInputs {
                locus,
                alleles,
                ploidy,
                n_samples,
                scalars,
                other_scalars,
                chain_anchor_flags,
            },
            config.max_gq_phred,
        ));
    }

    validate_record_shape(
        locus,
        ploidy,
        n_alleles,
        n_samples,
        n_genotypes,
        &log_likelihoods,
        &chain_anchor_flags,
    )?;

    let fixation_indices = resolve_fixation_indices(
        locus,
        n_samples,
        config.fixation_index_overrides.as_deref(),
        config.fixation_index_default,
    )?;

    // `safe_ln(f_s)` and `safe_ln(1 - f_s)` are functions of
    // `fixation_indices` only, so they're record-static. Precompute
    // once here instead of inside `e_step`'s per-iteration sample loop.
    let log_f_per_sample: Vec<f64> = fixation_indices
        .iter()
        .map(|&f_s| safe_ln(math, f_s))
        .collect();
    let log_one_minus_f_per_sample: Vec<f64> = fixation_indices
        .iter()
        .map(|&f_s| safe_ln(math, 1.0 - f_s))
        .collect();

    let allele_classes = classify_alleles(&alleles);
    let pseudocounts: Vec<f64> = allele_classes
        .iter()
        .map(|c| pseudocount_for(*c, config))
        .collect();
    let compound_mask: Vec<bool> = alleles.iter().map(|a| a.is_compound).collect();

    // Genotype-shape artefacts (allele-count table, multinomial
    // coefficients, non-zero-pair and homozygous-allele lookups) come
    // from a thread-local `(ploidy, n_alleles)`-keyed cache. Repeated
    // shapes (the common case — every biallelic-diploid record shares
    // `(2, 2)`) hit the cache instead of rebuilding.
    let shape = shape_for(ploidy, n_alleles);

    let ctx = EmContext {
        locus,
        n_samples,
        n_genotypes,
        n_alleles,
        ploidy,
        shape: &shape,
        compound_mask: &compound_mask,
        pseudocounts: &pseudocounts,
        log_f_per_sample: &log_f_per_sample,
        log_one_minus_f_per_sample: &log_one_minus_f_per_sample,
        compound_pseudocount: config.compound_alt_pseudocount,
    };

    let mut scratch = EmScratch::new(n_alleles, n_genotypes);
    let mut p_hat = vec![1.0 / n_alleles as f64; n_alleles];
    let mut f_hat_compound = vec![1.0 / n_alleles as f64; n_alleles];
    let mut posteriors = vec![0.0_f64; n_samples * n_genotypes];

    let diagnostics = run_em_loop(
        &ctx,
        config,
        math,
        &log_likelihoods,
        &mut p_hat,
        &mut f_hat_compound,
        &mut posteriors,
        &mut scratch,
    )?;

    // Final E-step under the converged (p̂, f̂_C) so the emitted
    // posteriors reflect the *post-final-M-step* parameters rather
    // than the parameters that produced the last `posteriors` buffer
    // contents.
    e_step(
        &ctx,
        math,
        &log_likelihoods,
        &p_hat,
        &f_hat_compound,
        &mut scratch,
        &mut posteriors,
    )?;

    let (best_genotype, gq_phred, qual_phred) =
        summarise_posteriors(n_samples, n_genotypes, &posteriors, config.max_gq_phred);

    let compound_frequencies: Vec<Option<f64>> = compound_mask
        .iter()
        .enumerate()
        .map(|(a, &is_compound)| {
            if is_compound {
                Some(f_hat_compound[a])
            } else {
                None
            }
        })
        .collect();

    Ok(PosteriorRecord {
        locus,
        alleles,
        ploidy,
        n_samples,
        n_genotypes,
        allele_frequencies: p_hat,
        compound_frequencies,
        posteriors,
        best_genotype,
        gq_phred,
        qual_phred,
        scalars,
        other_scalars,
        chain_anchor_flags,
        diagnostics,
    })
}

fn validate_record_shape(
    locus: RecordLocus,
    ploidy: u8,
    n_alleles: usize,
    n_samples: usize,
    n_genotypes: usize,
    log_likelihoods: &[f64],
    chain_anchor_flags: &[bool],
) -> Result<(), PosteriorEngineError> {
    let expected_n_genotypes = genotype_order(ploidy, n_alleles).len();
    if expected_n_genotypes != n_genotypes {
        return Err(PosteriorEngineError::MalformedMergedRecord {
            locus,
            reason: MalformedMergedRecordReason::GenotypeCountMismatch {
                expected: expected_n_genotypes,
                got: n_genotypes,
            },
        });
    }
    let expected_ll = n_samples * n_genotypes;
    if log_likelihoods.len() != expected_ll {
        return Err(PosteriorEngineError::MalformedMergedRecord {
            locus,
            reason: MalformedMergedRecordReason::LogLikelihoodLengthMismatch {
                expected: expected_ll,
                got: log_likelihoods.len(),
            },
        });
    }
    let expected_flags = n_samples * n_alleles;
    if chain_anchor_flags.len() != expected_flags {
        return Err(PosteriorEngineError::MalformedMergedRecord {
            locus,
            reason: MalformedMergedRecordReason::ChainAnchorLengthMismatch {
                expected: expected_flags,
                got: chain_anchor_flags.len(),
            },
        });
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_em_loop<M: MathBackend>(
    ctx: &EmContext<'_>,
    config: &PosteriorEngineConfig,
    math: &M,
    log_likelihoods: &[f64],
    p_hat: &mut Vec<f64>,
    f_hat_compound: &mut Vec<f64>,
    posteriors: &mut [f64],
    scratch: &mut EmScratch,
) -> Result<EmDiagnostics, PosteriorEngineError> {
    let mut iterations: u32 = 0;
    let mut last_delta = f64::INFINITY;
    let max_iterations = config.max_iterations;

    while iterations < max_iterations {
        iterations += 1;

        e_step(
            ctx,
            math,
            log_likelihoods,
            p_hat,
            f_hat_compound,
            scratch,
            posteriors,
        )?;

        let new_p_hat = m_step_p_hat(ctx, posteriors, scratch);
        let new_f_hat_compound = m_step_f_hat_compound(ctx, scratch);

        last_delta = max_abs_diff(p_hat, &new_p_hat);
        *p_hat = new_p_hat;
        *f_hat_compound = new_f_hat_compound;

        if last_delta < config.convergence_threshold {
            return Ok(EmDiagnostics {
                iterations,
                final_max_delta_p: last_delta,
            });
        }
    }

    Err(PosteriorEngineError::DidNotConverge {
        locus: ctx.locus,
        max_iterations,
        last_delta,
    })
}

fn resolve_fixation_indices(
    locus: RecordLocus,
    n_samples: usize,
    overrides: Option<&[f64]>,
    default: f64,
) -> Result<Vec<f64>, PosteriorEngineError> {
    let Some(overrides) = overrides else {
        return Ok(vec![default; n_samples]);
    };
    if overrides.len() != n_samples {
        return Err(PosteriorEngineError::InbreedingOverrideLengthMismatch {
            locus,
            override_len: overrides.len(),
            record_samples: n_samples,
        });
    }
    Ok(overrides.to_vec())
}

fn e_step<M: MathBackend>(
    ctx: &EmContext<'_>,
    math: &M,
    log_likelihoods: &[f64],
    p_hat: &[f64],
    f_hat_compound: &[f64],
    scratch: &mut EmScratch,
    posteriors: &mut [f64],
) -> Result<(), PosteriorEngineError> {
    // Effective per-allele frequency seen by the prior: `f̂_C`
    // substitutes for `p̂[compound]` in every compound slot, regardless
    // of any per-sample chain-broken flag (the chain-broken flag only
    // affects the likelihood, which Stage 5 has already baked into
    // `log_likelihoods`).
    for (a, (slot_p, slot_log)) in scratch
        .p_effective
        .iter_mut()
        .zip(scratch.log_p_effective.iter_mut())
        .enumerate()
    {
        let p = if ctx.compound_mask[a] {
            f_hat_compound[a]
        } else {
            p_hat[a]
        };
        *slot_p = p;
        *slot_log = safe_ln(math, p);
    }

    let ll_chunks = log_likelihoods.chunks_exact(ctx.n_genotypes);
    let post_chunks = posteriors.chunks_exact_mut(ctx.n_genotypes);

    for (sample_idx, (ll_row, post_row)) in ll_chunks.zip(post_chunks).enumerate() {
        // `log_f` and `log_one_minus_f` are precomputed in EmContext
        // because `fixation_indices` is record-static; this loop only
        // needs to read them indexed by `sample_idx`.
        let log_one_minus_f = ctx.log_one_minus_f_per_sample[sample_idx];
        let log_f = ctx.log_f_per_sample[sample_idx];

        for (g_idx, &ll) in ll_row.iter().enumerate() {
            // Inner sum walks only the non-zero (allele, count) pairs
            // for this genotype — the `if k == 0` skip-branch and the
            // per-cell `homozygous_allele` linear scan that the
            // un-shape-cached path used are both gone. For diploid each
            // pairs slice is 1 (homozygous) or 2 (heterozygous) entries.
            let (start, len) = ctx.shape.nonzero_pairs_offsets[g_idx];
            let pairs = &ctx.shape.nonzero_pairs[start as usize..start as usize + len as usize];

            let log_indep = ctx.shape.log_multinomial_coeffs[g_idx]
                + pairs
                    .iter()
                    .map(|&(a, k)| f64::from(k) * scratch.log_p_effective[a as usize])
                    .sum::<f64>();

            let log_prior = match ctx.shape.homozygous_allele_for[g_idx] {
                Some(homo_allele) => log_sum_exp_2(
                    math,
                    log_one_minus_f + log_indep,
                    log_f + scratch.log_p_effective[homo_allele as usize],
                ),
                // Heterozygous genotype: IBD component contributes 0.
                None => log_one_minus_f + log_indep,
            };
            scratch.log_post_unnorm[g_idx] = ll + log_prior;
        }

        let log_z = log_sum_exp_slice(math, &scratch.log_post_unnorm);
        // log_z is `-∞` iff every entry of `log_post_unnorm` is `-∞`,
        // which can only happen if the prior assigns zero mass to
        // every genotype the likelihood permits. Pseudocounts keep p̂
        // and f̂_C strictly positive, so the only realistic trigger is
        // a likelihood matrix that is all `-∞` — itself a malformed
        // upstream record. Surface it loudly.
        if !log_z.is_finite() {
            return Err(PosteriorEngineError::NonFinitePosterior {
                locus: ctx.locus,
                sample_idx,
                genotype_idx: None,
                kind: classify_nonfinite(log_z),
            });
        }

        for (g_idx, post_slot) in post_row.iter_mut().enumerate() {
            let posterior = math.exp(scratch.log_post_unnorm[g_idx] - log_z);
            if !posterior.is_finite() {
                return Err(PosteriorEngineError::NonFinitePosterior {
                    locus: ctx.locus,
                    sample_idx,
                    genotype_idx: Some(g_idx),
                    kind: classify_nonfinite(posterior),
                });
            }
            *post_slot = posterior;
        }
    }

    Ok(())
}

fn m_step_p_hat(ctx: &EmContext<'_>, posteriors: &[f64], scratch: &mut EmScratch) -> Vec<f64> {
    for slot in scratch.expected_counts.iter_mut() {
        *slot = 0.0;
    }

    for post_row in posteriors.chunks_exact(ctx.n_genotypes) {
        for (g_idx, gt_counts) in ctx
            .shape
            .genotype_allele_counts
            .chunks_exact(ctx.n_alleles)
            .enumerate()
        {
            let weight = post_row[g_idx];
            for (a, &k) in gt_counts.iter().enumerate() {
                if k != 0 {
                    scratch.expected_counts[a] += weight * k as f64;
                }
            }
        }
    }

    // M-step on p̂ — Dirichlet posterior mean over the simplex.
    let dirichlet_denominator: f64 = scratch
        .expected_counts
        .iter()
        .zip(ctx.pseudocounts.iter())
        .map(|(e, a)| e + a)
        .sum();
    scratch
        .expected_counts
        .iter()
        .zip(ctx.pseudocounts.iter())
        .map(|(e, a)| (e + a) / dirichlet_denominator)
        .collect()
}

fn m_step_f_hat_compound(ctx: &EmContext<'_>, scratch: &EmScratch) -> Vec<f64> {
    // M-step on f̂_C — Beta(α_compound, 1 − α_compound) posterior
    // mean. Each compound is independent of the others.
    let chromosomes_total = ctx.ploidy as f64 * ctx.n_samples as f64;
    let alpha_other = 1.0 - ctx.compound_pseudocount;
    (0..ctx.n_alleles)
        .map(|a| {
            if ctx.compound_mask[a] {
                let expected_compound_count = scratch.expected_counts[a];
                (ctx.compound_pseudocount + expected_compound_count)
                    / (ctx.compound_pseudocount + alpha_other + chromosomes_total)
            } else {
                0.0
            }
        })
        .collect()
}

fn summarise_posteriors(
    n_samples: usize,
    n_genotypes: usize,
    posteriors: &[f64],
    max_gq: f64,
) -> (Vec<usize>, Vec<f64>, f64) {
    let mut best_genotype = vec![0_usize; n_samples];
    let mut gq_phred = vec![0.0_f64; n_samples];
    let mut log10_p_hom_ref_sum = 0.0_f64;

    for (sample_idx, post_row) in posteriors.chunks_exact(n_genotypes).enumerate() {
        let (best_idx, best_p) = post_row.iter().enumerate().fold(
            (0_usize, f64::NEG_INFINITY),
            |(best_i, best_v), (i, &v)| {
                if v > best_v { (i, v) } else { (best_i, best_v) }
            },
        );
        best_genotype[sample_idx] = best_idx;

        // Clamp p_best below 1 by one ULP so the Phred calculation is
        // finite when EM produced an exact 1.0 (e.g. a singleton
        // genotype after softmax of `[finite, -∞, -∞]`).
        let p_best_clamped = best_p.min(1.0 - f64::EPSILON);
        let gq = PHRED_SCALE * (1.0 - p_best_clamped).log10();
        gq_phred[sample_idx] = gq.clamp(0.0, max_gq);

        let p_hom_ref = post_row[HOM_REF_GENOTYPE_IDX];
        if p_hom_ref > 0.0 {
            log10_p_hom_ref_sum += p_hom_ref.log10();
        } else {
            log10_p_hom_ref_sum = f64::NEG_INFINITY;
        }
    }

    // QUAL = -10 · log10(Π_s P(hom-ref)_s). If any sample has
    // P(hom-ref) = 0 exactly, the product is 0 and QUAL = +∞ — a
    // legitimate "site is certainly variant" signal; the VCF writer
    // can cap it if it wants to.
    let qual_phred = if log10_p_hom_ref_sum.is_finite() {
        PHRED_SCALE * log10_p_hom_ref_sum
    } else {
        f64::INFINITY
    };

    (best_genotype, gq_phred, qual_phred)
}

/// Inputs to [`trivial_posterior_record`]. Grouped into a struct
/// because the trivial-path fixture forwards seven distinct
/// `MergedRecord` fields verbatim and `too_many_arguments` would
/// otherwise fire.
struct TrivialRecordInputs {
    locus: RecordLocus,
    alleles: Vec<MergedAllele>,
    ploidy: u8,
    n_samples: usize,
    scalars: Vec<AlleleSupportStats>,
    other_scalars: Vec<AlleleSupportStats>,
    chain_anchor_flags: Vec<bool>,
}

fn trivial_posterior_record(inputs: TrivialRecordInputs, max_gq: f64) -> PosteriorRecord {
    let TrivialRecordInputs {
        locus,
        alleles,
        ploidy,
        n_samples,
        scalars,
        other_scalars,
        chain_anchor_flags,
    } = inputs;
    let n_alleles = alleles.len();
    let n_genotypes = genotype_order(ploidy, n_alleles).len();

    // Per-sample posterior row: all mass on hom-REF (genotype index 0).
    // For an empty allele set `n_genotypes == 0` and the buffer is
    // empty; otherwise each per-sample row sums to 1.
    let mut posteriors = vec![0.0_f64; n_samples * n_genotypes];
    if n_genotypes > 0 {
        for row in posteriors.chunks_exact_mut(n_genotypes) {
            row[HOM_REF_GENOTYPE_IDX] = 1.0;
        }
    }
    let best_genotype = vec![HOM_REF_GENOTYPE_IDX; n_samples];
    let gq_phred = vec![max_gq; n_samples];

    let allele_frequencies = if alleles.is_empty() {
        Vec::new()
    } else {
        let mut p = vec![0.0; n_alleles];
        p[0] = 1.0;
        p
    };
    let compound_frequencies = alleles
        .iter()
        .map(|a| if a.is_compound { Some(0.0) } else { None })
        .collect();

    PosteriorRecord {
        locus,
        alleles,
        ploidy,
        n_samples,
        n_genotypes,
        allele_frequencies,
        compound_frequencies,
        posteriors,
        best_genotype,
        gq_phred,
        qual_phred: 0.0,
        scalars,
        other_scalars,
        chain_anchor_flags,
        diagnostics: EmDiagnostics {
            iterations: 0,
            final_max_delta_p: 0.0,
        },
    }
}

// ---------------------------------------------------------------------
// Numeric helpers
// ---------------------------------------------------------------------

#[inline]
fn safe_ln<M: MathBackend>(math: &M, x: f64) -> f64 {
    if x <= 0.0 { f64::NEG_INFINITY } else { math.ln(x) }
}

#[inline]
fn log_sum_exp_2<M: MathBackend>(math: &M, a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let m = a.max(b);
    m + math.ln(math.exp(a - m) + math.exp(b - m))
}

#[inline]
fn log_sum_exp_slice<M: MathBackend>(math: &M, values: &[f64]) -> f64 {
    let m = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let sum: f64 = values.iter().map(|v| math.exp(*v - m)).sum();
    m + math.ln(sum)
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

fn classify_nonfinite(value: f64) -> NonFiniteKind {
    if value.is_nan() {
        NonFiniteKind::NaN
    } else if value == f64::INFINITY {
        NonFiniteKind::PositiveInfinity
    } else {
        NonFiniteKind::NegativeInfinity
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::per_group_merger::MergedAllele;

    // -----------------------------------------------------------------
    // Test fixtures
    // -----------------------------------------------------------------

    fn simple_alleles(seqs: &[&[u8]]) -> Vec<MergedAllele> {
        seqs.iter()
            .map(|s| MergedAllele {
                seq: s.to_vec(),
                is_compound: false,
                constituents: Vec::new(),
            })
            .collect()
    }

    /// Build a `MergedRecord` from per-sample, per-genotype log-likelihood
    /// vectors. `alleles[0]` is treated as REF; all alleles non-compound.
    fn merged_record_simple(
        chrom_id: u32,
        pos: u32,
        alleles: Vec<&[u8]>,
        ploidy: u8,
        likelihoods: Vec<Vec<f64>>,
    ) -> MergedRecord {
        let alleles = simple_alleles(&alleles);
        let n_samples = likelihoods.len();
        let n_alleles = alleles.len();
        let n_genotypes = genotype_order(ploidy, n_alleles).len();
        for (s, row) in likelihoods.iter().enumerate() {
            assert_eq!(
                row.len(),
                n_genotypes,
                "fixture sample {s} has {} genotypes; expected {n_genotypes}",
                row.len(),
            );
        }
        let ref_len = alleles[0].seq.len() as u32;
        let scalars = vec![AlleleSupportStats::default(); n_samples * n_alleles];
        let other_scalars = vec![AlleleSupportStats::default(); n_samples];
        let chain_anchor_flags = vec![false; n_samples * n_alleles];
        let log_likelihoods: Vec<f64> = likelihoods.into_iter().flatten().collect();
        MergedRecord {
            chrom_id,
            start: pos,
            end: pos + ref_len - 1,
            alleles,
            ploidy,
            n_samples,
            n_genotypes,
            scalars,
            other_scalars,
            chain_anchor_flags,
            log_likelihoods,
        }
    }

    /// Same as [`merged_record_simple`] but lets the caller mark
    /// specific alleles as compound and supply per-sample
    /// chain-anchor flags.
    fn merged_record_with_compound(
        chrom_id: u32,
        pos: u32,
        alleles: Vec<&[u8]>,
        compound_indices: &[usize],
        chain_anchor_flags_per_sample: Vec<Vec<bool>>,
        ploidy: u8,
        likelihoods: Vec<Vec<f64>>,
    ) -> MergedRecord {
        let mut record = merged_record_simple(chrom_id, pos, alleles, ploidy, likelihoods);
        for &i in compound_indices {
            record.alleles[i].is_compound = true;
        }
        let n_alleles = record.alleles.len();
        assert_eq!(chain_anchor_flags_per_sample.len(), record.n_samples);
        for (s, row) in chain_anchor_flags_per_sample.iter().enumerate() {
            assert_eq!(row.len(), n_alleles);
            for (a, &flag) in row.iter().enumerate() {
                record.chain_anchor_flags[s * n_alleles + a] = flag;
            }
        }
        record
    }

    // The semantic / numeric tests below assert tight tolerances
    // (often 1e-12) against hand-computed expected values. They pin
    // against `ExactMath` explicitly so the production default
    // (`InterpUnivariateMath`) can move without breaking them; the
    // backend-equivalence guarantee is enforced separately by
    // `tests_math_backend_accuracy` (see the accuracy harness).
    fn engine_for(record: MergedRecord) -> Vec<Result<PosteriorRecord, PosteriorEngineError>> {
        engine_for_with_config(record, PosteriorEngineConfig::with_project_defaults())
    }

    fn engine_for_with_config(
        record: MergedRecord,
        config: PosteriorEngineConfig,
    ) -> Vec<Result<PosteriorRecord, PosteriorEngineError>> {
        let upstream = std::iter::once(Ok::<_, PerGroupMergerError>(record));
        PosteriorEngine::with_math_backend(upstream, config, super::backends::ExactMath).collect()
    }

    fn single_ok(record: MergedRecord) -> PosteriorRecord {
        let mut out = engine_for(record);
        assert_eq!(out.len(), 1);
        out.remove(0).expect("posterior record")
    }

    /// Proptest-friendly variant of [`single_ok`]: returns the engine's
    /// `Result` so the caller can `prop_assume!` over pathological
    /// inputs (e.g. `DidNotConverge`, which is a documented engine
    /// error variant for adversarial likelihood matrices, not a bug).
    fn try_single(record: MergedRecord) -> Result<PosteriorRecord, PosteriorEngineError> {
        let mut out = engine_for(record);
        assert_eq!(out.len(), 1);
        out.remove(0)
    }

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // -----------------------------------------------------------------
    // Unit tests — no-contamination mode (per the plan §"Unit tests")
    // -----------------------------------------------------------------

    #[test]
    fn single_sample_strong_ref_returns_homozygous_ref_and_low_qual() {
        let record =
            merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![0.0, -50.0, -50.0]]);
        let pr = single_ok(record);
        assert_eq!(pr.best_genotype, vec![0]);
        assert!(pr.posteriors_row(0)[0] > 0.999);
        assert!(pr.qual_phred < 1.0, "qual_phred = {}", pr.qual_phred);
        assert!(pr.allele_frequencies[0] > 0.9);
        assert!(pr.allele_frequencies[1] < 0.1);
    }

    #[test]
    fn single_sample_strong_alt_returns_homozygous_alt_and_high_qual() {
        let record =
            merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![-50.0, -50.0, 0.0]]);
        let pr = single_ok(record);
        assert_eq!(pr.best_genotype, vec![2]);
        assert!(pr.posteriors_row(0)[2] > 0.999);
        assert!(pr.qual_phred > 20.0, "qual_phred = {}", pr.qual_phred);
    }

    #[test]
    fn single_sample_uniform_likelihood_leans_on_prior_and_pseudocount() {
        let record =
            merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![-1.0, -1.0, -1.0]]);
        let pr = single_ok(record);
        let row = pr.posteriors_row(0);
        let sum: f64 = row.iter().sum();
        assert!(approx(sum, 1.0, 1e-9), "row sum = {sum}");
        assert!(pr.allele_frequencies[0] > pr.allele_frequencies[1]);
        assert_eq!(pr.best_genotype, vec![0]);
    }

    #[test]
    fn two_samples_opposite_evidence_pulls_p_hat_toward_one_half() {
        let record = merged_record_simple(
            1,
            100,
            vec![b"A", b"C"],
            2,
            vec![
                vec![0.0, -50.0, -50.0], // S0: strong RR
                vec![-50.0, -50.0, 0.0], // S1: strong AA
            ],
        );
        let pr = single_ok(record);
        assert_eq!(pr.best_genotype, vec![0, 2]);
        assert!(pr.allele_frequencies[0] > 0.5);
        assert!(pr.allele_frequencies[1] > 0.1);
    }

    #[test]
    fn cohort_pooling_pulls_weak_het_sample_toward_ref_when_pseudocount_dominates() {
        let mut likelihoods: Vec<Vec<f64>> = (0..19).map(|_| vec![0.0, -50.0, -50.0]).collect();
        likelihoods.push(vec![-2.0, 0.0, -2.0]);
        let record = merged_record_simple(1, 100, vec![b"A", b"C"], 2, likelihoods);
        let pr = single_ok(record);
        assert!(
            pr.allele_frequencies[1] < 0.05,
            "p̂[alt] = {}",
            pr.allele_frequencies[1]
        );
        assert_eq!(pr.best_genotype[19], 0);
    }

    #[test]
    fn em_converges_within_a_handful_of_iterations() {
        let record = merged_record_simple(
            1,
            100,
            vec![b"A", b"C"],
            2,
            vec![
                vec![0.0, -50.0, -50.0],
                vec![-50.0, -50.0, 0.0],
                vec![-50.0, 0.0, -50.0],
                vec![0.0, -50.0, -50.0],
            ],
        );
        let pr = single_ok(record);
        assert!(
            pr.diagnostics.iterations <= 10,
            "EM took {} iterations",
            pr.diagnostics.iterations
        );
    }

    #[test]
    fn hard_iteration_cap_returns_did_not_converge_error() {
        let record = merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![0.0, -1.0, -2.0]]);
        let config = PosteriorEngineConfig {
            max_iterations: 1,
            convergence_threshold: 1e-12,
            ..Default::default()
        };
        let out = engine_for_with_config(record, config);
        match &out[0] {
            Err(PosteriorEngineError::DidNotConverge {
                locus,
                max_iterations,
                last_delta,
            }) => {
                assert_eq!(locus.chrom_id, 1);
                assert_eq!(locus.start, 100);
                assert_eq!(locus.end, 100);
                assert_eq!(*max_iterations, 1);
                assert!(
                    last_delta.is_finite() && *last_delta >= 0.0,
                    "last_delta = {last_delta}"
                );
            }
            other => panic!("expected DidNotConverge, got {other:?}"),
        }
    }

    #[test]
    fn inbreeding_one_zeroes_heterozygote_posterior() {
        let record = merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![-2.0, 0.0, -2.0]]);
        let config = PosteriorEngineConfig {
            fixation_index_default: 1.0,
            ..Default::default()
        };
        let out = engine_for_with_config(record, config);
        let pr = out.into_iter().next().unwrap().expect("posterior");
        assert!(pr.posteriors_row(0)[1].abs() < 1e-12);
        assert_ne!(pr.best_genotype[0], 1);
    }

    #[test]
    fn fixation_index_overrides_apply_per_sample_when_length_matches_record() {
        // Both samples have strong het likelihood. Sample 0's F=0
        // prior allows the het call; sample 1's F=1 zeroes het mass
        // entirely so the call must be a homozygote regardless of
        // likelihood strength.
        let record = merged_record_simple(
            1,
            100,
            vec![b"A", b"C"],
            2,
            vec![vec![-50.0, 0.0, -50.0], vec![-50.0, 0.0, -50.0]],
        );
        let config = PosteriorEngineConfig {
            fixation_index_overrides: Some(vec![0.0, 1.0]),
            ..Default::default()
        };
        let pr = engine_for_with_config(record, config)
            .into_iter()
            .next()
            .unwrap()
            .expect("posterior");
        // Sample 0 (F=0): heterozygote prior is finite, strong het
        // likelihood wins.
        assert_eq!(pr.best_genotype[0], 1);
        // Sample 1 (F=1): heterozygote prior is 0, forced to a homozygote
        // even though the likelihood prefers het.
        assert_ne!(pr.best_genotype[1], 1);
    }

    #[test]
    fn snp_vs_indel_pseudocount_routing_pulls_indel_alt_down_more() {
        let record = merged_record_simple(1, 100, vec![b"A", b"C", b"AT"], 2, vec![vec![0.0; 6]]);
        let pr = single_ok(record);
        assert!(
            pr.allele_frequencies[1] > pr.allele_frequencies[2],
            "p̂[snp_alt] = {} not greater than p̂[indel_alt] = {}",
            pr.allele_frequencies[1],
            pr.allele_frequencies[2]
        );
    }

    #[test]
    fn compound_allele_emits_compound_frequency_in_posterior_record() {
        let record = merged_record_with_compound(
            1,
            100,
            vec![b"A", b"C"],
            &[1],
            vec![vec![false, false], vec![false, false]],
            2,
            vec![vec![-50.0, -50.0, 0.0], vec![-50.0, -50.0, 0.0]],
        );
        let pr = single_ok(record);
        assert!(pr.compound_frequencies[0].is_none());
        let f_c = pr.compound_frequencies[1].expect("compound 1 frequency");
        assert!(f_c > 0.0 && f_c < 1.0, "f̂_C = {f_c}");
        assert_eq!(pr.best_genotype, vec![2, 2]);
    }

    #[test]
    fn m_step_compound_frequency_matches_closed_form_with_no_observed_compound_counts() {
        // 2 samples × ploidy 2 = 4 chromosomes; both samples strongly RR ⇒ E[n_C] ≈ 0.
        // f̂_C = α / (α + (1 − α) + 4) = α / 5.
        let record = merged_record_with_compound(
            1,
            100,
            vec![b"A", b"C"],
            &[1],
            vec![vec![false, false], vec![false, false]],
            2,
            vec![vec![0.0, -50.0, -50.0], vec![0.0, -50.0, -50.0]],
        );
        let pr = single_ok(record);
        let f_c = pr.compound_frequencies[1].expect("compound 1 frequency");
        let alpha = DEFAULT_COMPOUND_ALT_PSEUDOCOUNT;
        let expected = alpha / (alpha + (1.0 - alpha) + 2.0 * 2.0);
        assert!(
            approx(f_c, expected, 1e-6),
            "f̂_C = {f_c}, expected {expected}"
        );
    }

    #[test]
    fn chain_broken_sample_uses_likelihood_unchanged_and_still_emits_compound_frequency() {
        let record = merged_record_with_compound(
            1,
            100,
            vec![b"A", b"C"],
            &[1],
            vec![vec![false, true], vec![false, false]],
            2,
            vec![vec![-50.0, -50.0, 0.0], vec![-50.0, -50.0, 0.0]],
        );
        let pr = single_ok(record);
        assert!(pr.chain_anchor_flags_row(0)[1]);
        assert!(!pr.chain_anchor_flags_row(1)[1]);
        assert!(pr.compound_frequencies[1].is_some());
    }

    #[test]
    fn all_samples_chain_broken_at_compound_remains_numerically_stable() {
        let record = merged_record_with_compound(
            1,
            100,
            vec![b"A", b"C"],
            &[1],
            vec![vec![false, true], vec![false, true]],
            2,
            vec![vec![-10.0, -10.0, 0.0], vec![-10.0, -10.0, 0.0]],
        );
        let pr = single_ok(record);
        for s in 0..pr.n_samples {
            for v in pr.posteriors_row(s) {
                assert!(v.is_finite() && *v >= 0.0 && *v <= 1.0);
            }
        }
        assert!(pr.qual_phred.is_finite());
        assert!(pr.compound_frequencies[1].is_some());
    }

    #[test]
    fn site_qual_matches_closed_form_on_weak_uniform_cohort() {
        let likelihoods: Vec<Vec<f64>> = (0..5).map(|_| vec![0.0, -3.0, -3.0]).collect();
        let record = merged_record_simple(1, 100, vec![b"A", b"C"], 2, likelihoods);
        let pr = single_ok(record);
        let expected_qual: f64 = -10.0
            * (0..pr.n_samples)
                .map(|s| pr.posteriors_row(s)[0].log10())
                .sum::<f64>();
        assert!(
            approx(pr.qual_phred, expected_qual, 1e-9),
            "qual_phred={} expected={}",
            pr.qual_phred,
            expected_qual
        );
    }

    #[test]
    fn gq_phred_matches_closed_form() {
        let record =
            merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![0.0, -50.0, -50.0]]);
        let pr = single_ok(record);
        let p_best = pr.posteriors_row(0)[pr.best_genotype[0]];
        let p_best_clamped = p_best.min(1.0 - f64::EPSILON);
        let expected_gq = ((-10.0_f64) * (1.0 - p_best_clamped).log10()).clamp(0.0, MAX_GQ_PHRED);
        assert!(
            approx(pr.gq_phred[0], expected_gq, 1e-9),
            "gq_phred={} expected={}",
            pr.gq_phred[0],
            expected_gq
        );
    }

    #[test]
    fn summarise_posteriors_returns_infinite_qual_when_any_sample_has_zero_hom_ref_posterior() {
        // log-likelihood ratio so extreme that softmax of the alts
        // underflows P(RR) to exactly 0.
        let record = merged_record_simple(
            1,
            100,
            vec![b"A", b"C"],
            2,
            vec![vec![-1000.0, -1000.0, 0.0]],
        );
        let pr = single_ok(record);
        if pr.posteriors_row(0)[0] == 0.0 {
            assert!(
                pr.qual_phred.is_infinite() && pr.qual_phred > 0.0,
                "qual_phred = {}",
                pr.qual_phred
            );
        }
        // GQ clamps regardless of how extreme the posterior is.
        assert!(pr.gq_phred[0].is_finite());
        assert!(pr.gq_phred[0] <= MAX_GQ_PHRED + 1e-9);
    }

    #[test]
    fn approximate_mode_returns_not_yet_implemented_error() {
        let record =
            merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![0.0, -50.0, -50.0]]);
        let config = PosteriorEngineConfig {
            approximate_posterior_calculation: true,
            ..Default::default()
        };
        let out = engine_for_with_config(record, config);
        match &out[0] {
            Err(PosteriorEngineError::ApproximateModeNotYetImplemented) => {}
            other => panic!("expected ApproximateModeNotYetImplemented, got {other:?}"),
        }
    }

    #[test]
    fn contamination_estimates_zero_matches_no_contamination_path() {
        // Passing `ContaminationEstimates::zero(...)` (all c_s = None,
        // every batch zeroed) should produce a posterior record
        // identical to the `contamination = None` path because the
        // mixture-likelihood helper falls back to Stage 5's
        // precomputed log_likelihoods for every below-floor sample.
        let record =
            merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![0.0, -50.0, -50.0]]);
        let baseline = engine_for(record.clone());

        let config = PosteriorEngineConfig {
            contamination: Some(ContaminationEstimates::zero(1, vec![0], 1).expect("valid inputs")),
            ..Default::default()
        };
        let mirrored = engine_for_with_config(record, config);

        match (&baseline[0], &mirrored[0]) {
            (Ok(b), Ok(m)) => assert_eq!(b, m),
            other => panic!("expected matching Ok records, got {other:?}"),
        }
    }

    #[test]
    fn contamination_cohort_size_mismatch_surfaces_typed_error() {
        // ContaminationEstimates built for a 2-sample cohort, applied
        // to a 1-sample record → typed mismatch error.
        let record =
            merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![0.0, -50.0, -50.0]]);
        let config = PosteriorEngineConfig {
            contamination: Some(
                ContaminationEstimates::zero(2, vec![0, 0], 1).expect("valid inputs"),
            ),
            ..Default::default()
        };
        let out = engine_for_with_config(record, config);
        match &out[0] {
            Err(PosteriorEngineError::ContaminationCohortSizeMismatch {
                estimates_n_samples: 2,
                record_samples: 1,
                ..
            }) => {}
            other => panic!("expected ContaminationCohortSizeMismatch, got {other:?}"),
        }
    }

    // ---------- B2 value-witness tests for compute_mixture_log_likelihoods ----------

    /// Build a 1-sample, 2-allele (REF=A, SNP-alt=C) MergedRecord with
    /// the given per-allele scalars and a stub fallback log-likelihood
    /// table. The fallback is shape-correct but only used by the
    /// floored-sample branch (`c_s = 0`), so the value is irrelevant
    /// for c_s > 0 tests.
    fn build_single_sample_record(
        num_obs_ref: u32,
        num_obs_alt: u32,
        mean_err: f64,
        n_genotypes: usize,
    ) -> MergedRecord {
        let alleles = vec![
            MergedAllele {
                seq: b"A".to_vec(),
                is_compound: false,
                constituents: vec![],
            },
            MergedAllele {
                seq: b"C".to_vec(),
                is_compound: false,
                constituents: vec![],
            },
        ];
        let q_sum_ref = if num_obs_ref > 0 {
            mean_err.ln() * f64::from(num_obs_ref)
        } else {
            0.0
        };
        let q_sum_alt = if num_obs_alt > 0 {
            mean_err.ln() * f64::from(num_obs_alt)
        } else {
            0.0
        };
        let scalars = vec![
            AlleleSupportStats::new(num_obs_ref, q_sum_ref, num_obs_ref / 2, 0, 0),
            AlleleSupportStats::new(num_obs_alt, q_sum_alt, num_obs_alt / 2, 0, 0),
        ];
        MergedRecord {
            chrom_id: 1,
            start: 100,
            end: 100,
            alleles,
            ploidy: 2,
            n_samples: 1,
            n_genotypes,
            scalars,
            other_scalars: vec![AlleleSupportStats::default(); 1],
            chain_anchor_flags: vec![false; 2],
            log_likelihoods: vec![0.0; n_genotypes], // fallback stub
        }
    }

    #[test]
    fn compute_mixture_log_likelihoods_matches_hand_calculation_at_c_s_5_percent() {
        // Hand calc for: 1 sample, ploidy=2, alleles=[A, C], n_obs=[50, 0],
        // ε = 0.001 per read, c_s = 0.05, q_b = [0.95, 0.05, 0.0].
        // n_alleles=2 ⇒ non_genotype error denom = 1.0.
        //
        // For each genotype G with copy count k_REF in {0, 1, 2}:
        //   p_own(REF | G) = (k/2)·(1-ε) + ((2-k)/2)·(ε/1)
        //                  = (k/2)·0.999 + ((2-k)/2)·0.001
        //   q_b[Ref] = 0.95
        //   mix = 0.95·p_own + 0.05·0.95
        //   ll  = 50·ln(mix)         [n_obs only on REF; ALT contributes 0]
        //
        //   G=RR (k=2): p_own=0.999    mix = 0.95·0.999 + 0.05·0.95 ≈ 0.9985
        //   G=RA (k=1): p_own=0.5005   mix = 0.95·0.5005 + 0.05·0.95 ≈ 0.522975
        //   G=AA (k=0): p_own=0.001    mix = 0.95·0.001 + 0.05·0.95 ≈ 0.04845
        let estimates = ContaminationEstimates::from_user_supplied(
            vec![Some(0.05)],
            vec![[0.95, 0.05, 0.0]],
            vec![0],
        )
        .expect("valid inputs");
        let n_genotypes = genotype_order(2, 2).len();
        assert_eq!(n_genotypes, 3); // RR, RA, AA

        let record = build_single_sample_record(50, 0, 0.001, n_genotypes);
        let config = PosteriorEngineConfig {
            contamination: Some(estimates),
            ..Default::default()
        };
        let out = engine_for_with_config(record, config);
        let pr = out.into_iter().next().unwrap().expect("ok");

        // The EM converges these likelihoods into posteriors. We
        // can't easily hand-calc the EM-converged posteriors but we
        // CAN check the relative ordering implied by the
        // log-likelihoods (RR > RA > AA, monotonic in the witness
        // mix calc above).
        let post = pr.posteriors_row(0);
        assert!(
            post[0] > post[1] && post[1] > post[2],
            "mixture-likelihood ordering wrong: RR={}, RA={}, AA={} \
             (expected strictly decreasing)",
            post[0],
            post[1],
            post[2],
        );
        // RR posterior should dominate (the mix value at G=RR is
        // ~21× the mix value at G=RA, multiplied by 50 reads).
        assert!(post[0] > 0.95, "P(RR) should dominate; got {}", post[0]);
    }

    #[test]
    fn compute_mixture_log_likelihoods_falls_back_to_precomputed_when_c_s_is_none() {
        // When `effective_c_s` returns 0 (sample's `c_s_per_sample[s]`
        // is None — singleton/below-floor), the mixture function
        // copies Stage 5's `log_likelihoods` row through verbatim
        // rather than recomputing. We can witness this by passing a
        // recognisable fallback that produces a specific posterior,
        // then asserting the engine's output matches a baseline run
        // with `contamination = None` on the same record.
        let alleles = vec![
            MergedAllele {
                seq: b"A".to_vec(),
                is_compound: false,
                constituents: vec![],
            },
            MergedAllele {
                seq: b"C".to_vec(),
                is_compound: false,
                constituents: vec![],
            },
        ];
        let scalars = vec![
            AlleleSupportStats::new(49, 0.001_f64.ln() * 49.0, 24, 0, 0),
            AlleleSupportStats::new(1, 0.001_f64.ln(), 0, 0, 0),
        ];
        let n_genotypes = genotype_order(2, 2).len();
        // Use a recognisable hom-REF-leaning likelihood vector.
        let lls = vec![0.0_f64, -30.0, -50.0];
        let make_record = || MergedRecord {
            chrom_id: 1,
            start: 100,
            end: 100,
            alleles: alleles.clone(),
            ploidy: 2,
            n_samples: 1,
            n_genotypes,
            scalars: scalars.clone(),
            other_scalars: vec![AlleleSupportStats::default(); 1],
            chain_anchor_flags: vec![false; 2],
            log_likelihoods: lls.clone(),
        };

        // contamination = Some(estimates with c_s = None) → fallback path.
        let estimates = ContaminationEstimates::zero(1, vec![0], 1).expect("valid inputs");
        let config_contam = PosteriorEngineConfig {
            contamination: Some(estimates),
            ..Default::default()
        };
        let pr_contam = engine_for_with_config(make_record(), config_contam);
        let pr_baseline = engine_for(make_record());

        match (&pr_contam[0], &pr_baseline[0]) {
            (Ok(c), Ok(b)) => assert_eq!(
                c, b,
                "c_s=None should produce identical posteriors to contamination=None"
            ),
            other => panic!("expected matching Ok records, got {other:?}"),
        }
    }

    #[test]
    fn compute_mixture_log_likelihoods_returns_non_finite_error_when_mix_is_zero() {
        // c_s = 1.0 (all reads are contaminant), q_b = [1.0, 0.0, 0.0]
        // (contaminant always carries REF). A read carrying the SNP
        // alt (allele 1) has p_contam = 0 and p_own term is zeroed by
        // (1 - c_s) = 0 → mix = 0 → log = -∞.
        let estimates = ContaminationEstimates::from_user_supplied(
            vec![Some(1.0)],
            vec![[1.0, 0.0, 0.0]],
            vec![0],
        )
        .expect("valid inputs");
        let n_genotypes = genotype_order(2, 2).len();
        // 1 ALT read forces the mixture term on allele 1 to fire.
        let record = build_single_sample_record(0, 1, 0.001, n_genotypes);
        let config = PosteriorEngineConfig {
            contamination: Some(estimates),
            ..Default::default()
        };
        let out = engine_for_with_config(record, config);
        match &out[0] {
            Err(PosteriorEngineError::NonFinitePosterior {
                kind: NonFiniteKind::NegativeInfinity,
                sample_idx: 0,
                ..
            }) => {}
            other => panic!("expected NonFinitePosterior(NegativeInfinity), got {other:?}"),
        }
    }

    #[test]
    fn compute_mixture_log_likelihoods_invariant_under_q_b_when_zero_reads() {
        // A sample with `n_obs = [0, 0]` for every allele should
        // produce posteriors that are *independent* of `q_b` — only
        // the EM prior shapes the result. Verify by running with two
        // distinct `q_b` rows and asserting matching posteriors.
        let n_genotypes = genotype_order(2, 2).len();
        let make_record = || build_single_sample_record(0, 0, 0.001, n_genotypes);

        let est_a = ContaminationEstimates::from_user_supplied(
            vec![Some(0.05)],
            vec![[0.95, 0.05, 0.0]],
            vec![0],
        )
        .expect("valid inputs");
        let est_b = ContaminationEstimates::from_user_supplied(
            vec![Some(0.05)],
            vec![[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]],
            vec![0],
        )
        .expect("valid inputs");

        let pr_a = engine_for_with_config(
            make_record(),
            PosteriorEngineConfig {
                contamination: Some(est_a),
                ..Default::default()
            },
        );
        let pr_b = engine_for_with_config(
            make_record(),
            PosteriorEngineConfig {
                contamination: Some(est_b),
                ..Default::default()
            },
        );

        match (&pr_a[0], &pr_b[0]) {
            (Ok(a), Ok(b)) => assert_eq!(
                a.posteriors, b.posteriors,
                "posteriors should be independent of q_b when n_obs = [0; K]"
            ),
            other => panic!("expected matching Ok records, got {other:?}"),
        }
    }

    // ---------- B2 / Mi12: map_to_contam_class direct coverage ----------

    #[test]
    fn map_to_contam_class_round_trips_all_three_non_compound_classes() {
        assert_eq!(
            map_to_contam_class(AlleleClass::Ref),
            Some(ContamAlleleClass::Ref)
        );
        assert_eq!(
            map_to_contam_class(AlleleClass::SnpAlt),
            Some(ContamAlleleClass::SnpAlt)
        );
        assert_eq!(
            map_to_contam_class(AlleleClass::IndelAlt),
            Some(ContamAlleleClass::IndelAlt)
        );
    }

    #[test]
    fn map_to_contam_class_returns_none_for_compound() {
        assert!(map_to_contam_class(AlleleClass::CompoundAlt).is_none());
    }

    #[test]
    fn inbreeding_override_length_mismatch_surfaces_typed_error() {
        let record =
            merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![0.0, -50.0, -50.0]]);
        // 2-entry override vs 1-sample record ⇒ mismatch.
        let config = PosteriorEngineConfig {
            fixation_index_overrides: Some(vec![0.5, 0.5]),
            ..Default::default()
        };
        let out = engine_for_with_config(record, config);
        match &out[0] {
            Err(PosteriorEngineError::InbreedingOverrideLengthMismatch {
                override_len,
                record_samples,
                ..
            }) => {
                assert_eq!(*override_len, 2);
                assert_eq!(*record_samples, 1);
            }
            other => panic!("expected InbreedingOverrideLengthMismatch, got {other:?}"),
        }
    }

    #[test]
    fn error_latches_engine_after_first_failure() {
        let record1 =
            merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![0.0, -50.0, -50.0]]);
        let err = PerGroupMergerError::RefFetch {
            chrom_id: 1,
            start: 100,
            end: 100,
            source: std::io::Error::other("synthetic"),
        };
        let record2 =
            merged_record_simple(1, 200, vec![b"A", b"C"], 2, vec![vec![0.0, -50.0, -50.0]]);
        let upstream = vec![Ok(record1), Err(err), Ok(record2)].into_iter();
        let out: Vec<_> = PosteriorEngine::new(upstream).collect();
        assert!(out[0].is_ok());
        assert!(out[1].is_err());
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn engine_emits_all_successful_records_before_latched_upstream_error() {
        let r =
            |pos| merged_record_simple(1, pos, vec![b"A", b"C"], 2, vec![vec![0.0, -50.0, -50.0]]);
        let err = PerGroupMergerError::RefFetch {
            chrom_id: 1,
            start: 400,
            end: 400,
            source: std::io::Error::other("synthetic"),
        };
        let upstream = vec![Ok(r(100)), Ok(r(200)), Ok(r(300)), Err(err), Ok(r(500))].into_iter();
        let out: Vec<_> = PosteriorEngine::new(upstream).collect();
        assert_eq!(out.len(), 4);
        assert!(out[0].is_ok() && out[1].is_ok() && out[2].is_ok());
        assert!(matches!(out[3], Err(PosteriorEngineError::Upstream(_))));
        let starts: Vec<u32> = out
            .iter()
            .take(3)
            .map(|r| r.as_ref().unwrap().locus.start)
            .collect();
        assert_eq!(starts, vec![100, 200, 300]);
    }

    #[test]
    fn upstream_exhaustion_is_clean_termination_not_error() {
        let upstream: std::iter::Empty<Result<MergedRecord, PerGroupMergerError>> =
            std::iter::empty();
        let out: Vec<_> = PosteriorEngine::new(upstream).collect();
        assert!(out.is_empty());
    }

    #[test]
    fn triploid_run_succeeds_with_polyploid_hwe_prior() {
        // Triploid biallelic: genotype_order = [000, 001, 011, 111].
        let record = merged_record_simple(
            1,
            100,
            vec![b"A", b"C"],
            3,
            vec![vec![0.0, -50.0, -50.0, -50.0]],
        );
        let pr = single_ok(record);
        assert_eq!(pr.n_genotypes, 4);
        assert_eq!(pr.best_genotype, vec![0]);
        let s: f64 = pr.posteriors_row(0).iter().sum();
        assert!(approx(s, 1.0, 1e-9));
    }

    // -----------------------------------------------------------------
    // Internal-bug guards (NonFinitePosterior + MalformedMergedRecord)
    // -----------------------------------------------------------------

    #[test]
    fn e_step_returns_non_finite_posterior_when_likelihood_row_is_all_negative_infinity() {
        let record = merged_record_simple(
            7,
            42,
            vec![b"A", b"C"],
            2,
            vec![vec![
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
            ]],
        );
        let out = engine_for(record);
        match &out[0] {
            Err(PosteriorEngineError::NonFinitePosterior {
                locus,
                sample_idx,
                genotype_idx,
                kind,
            }) => {
                assert_eq!(locus.chrom_id, 7);
                assert_eq!(locus.start, 42);
                assert_eq!(locus.end, 42);
                assert_eq!(*sample_idx, 0);
                assert!(genotype_idx.is_none());
                assert_eq!(*kind, NonFiniteKind::NegativeInfinity);
            }
            other => panic!("expected NonFinitePosterior, got {other:?}"),
        }
    }

    #[test]
    fn malformed_merged_record_detected_for_wrong_n_genotypes() {
        // 2 alleles + ploidy 2 ⇒ genotype_order yields 3 genotypes;
        // here we claim 4. The engine must surface a typed error
        // rather than indexing past the log_likelihoods buffer.
        let record = MergedRecord {
            chrom_id: 1,
            start: 100,
            end: 100,
            alleles: simple_alleles(&[b"A", b"C"]),
            ploidy: 2,
            n_samples: 1,
            n_genotypes: 4,
            scalars: vec![AlleleSupportStats::default(); 2],
            other_scalars: vec![AlleleSupportStats::default(); 1],
            chain_anchor_flags: vec![false; 2],
            log_likelihoods: vec![0.0; 4],
        };
        let out = engine_for(record);
        match &out[0] {
            Err(PosteriorEngineError::MalformedMergedRecord {
                reason: MalformedMergedRecordReason::GenotypeCountMismatch { expected, got },
                ..
            }) => {
                assert_eq!(*expected, 3);
                assert_eq!(*got, 4);
            }
            other => panic!("expected MalformedMergedRecord::GenotypeCountMismatch, got {other:?}"),
        }
    }

    #[test]
    fn malformed_merged_record_detected_for_wrong_log_likelihoods_length() {
        let record = MergedRecord {
            chrom_id: 1,
            start: 100,
            end: 100,
            alleles: simple_alleles(&[b"A", b"C"]),
            ploidy: 2,
            n_samples: 1,
            n_genotypes: 3,
            scalars: vec![AlleleSupportStats::default(); 2],
            other_scalars: vec![AlleleSupportStats::default(); 1],
            chain_anchor_flags: vec![false; 2],
            log_likelihoods: vec![0.0; 5], // expected 3
        };
        let out = engine_for(record);
        match &out[0] {
            Err(PosteriorEngineError::MalformedMergedRecord {
                reason: MalformedMergedRecordReason::LogLikelihoodLengthMismatch { expected, got },
                ..
            }) => {
                assert_eq!(*expected, 3);
                assert_eq!(*got, 5);
            }
            other => {
                panic!("expected MalformedMergedRecord::LogLikelihoodLengthMismatch, got {other:?}")
            }
        }
    }

    #[test]
    fn malformed_merged_record_detected_for_wrong_chain_anchor_flags_length() {
        let record = MergedRecord {
            chrom_id: 1,
            start: 100,
            end: 100,
            alleles: simple_alleles(&[b"A", b"C"]),
            ploidy: 2,
            n_samples: 1,
            n_genotypes: 3,
            scalars: vec![AlleleSupportStats::default(); 2],
            other_scalars: vec![AlleleSupportStats::default(); 1],
            chain_anchor_flags: vec![false; 5], // expected 2
            log_likelihoods: vec![0.0; 3],
        };
        let out = engine_for(record);
        match &out[0] {
            Err(PosteriorEngineError::MalformedMergedRecord {
                reason: MalformedMergedRecordReason::ChainAnchorLengthMismatch { expected, got },
                ..
            }) => {
                assert_eq!(*expected, 2);
                assert_eq!(*got, 5);
            }
            other => {
                panic!("expected MalformedMergedRecord::ChainAnchorLengthMismatch, got {other:?}")
            }
        }
    }

    // -----------------------------------------------------------------
    // Trivial-record path (alleles.len() < 2)
    // -----------------------------------------------------------------

    #[test]
    fn trivial_posterior_record_returns_hom_ref_per_sample_when_only_ref_allele_present() {
        let record = MergedRecord {
            chrom_id: 3,
            start: 7,
            end: 7,
            alleles: simple_alleles(&[b"A"]),
            ploidy: 2,
            n_samples: 2,
            n_genotypes: 1,
            scalars: vec![AlleleSupportStats::default(); 2],
            other_scalars: vec![AlleleSupportStats::default(); 2],
            chain_anchor_flags: vec![false; 2],
            log_likelihoods: vec![0.0; 2],
        };
        let pr = single_ok(record);
        assert_eq!(pr.best_genotype, vec![0, 0]);
        assert_eq!(pr.qual_phred, 0.0);
        assert_eq!(pr.diagnostics.iterations, 0);
        for s in 0..pr.n_samples {
            let sum: f64 = pr.posteriors_row(s).iter().sum();
            assert!(approx(sum, 1.0, 1e-12), "row sum = {sum}");
        }
        let p_sum: f64 = pr.allele_frequencies.iter().sum();
        assert!(approx(p_sum, 1.0, 1e-12));
    }

    #[test]
    fn trivial_posterior_record_returns_empty_allele_frequencies_when_alleles_empty() {
        let record = MergedRecord {
            chrom_id: 0,
            start: 1,
            end: 1,
            alleles: Vec::new(),
            ploidy: 2,
            n_samples: 1,
            n_genotypes: 0,
            scalars: Vec::new(),
            other_scalars: vec![AlleleSupportStats::default(); 1],
            chain_anchor_flags: Vec::new(),
            log_likelihoods: Vec::new(),
        };
        let pr = single_ok(record);
        assert!(pr.allele_frequencies.is_empty());
        assert!(pr.compound_frequencies.is_empty());
        assert!(pr.posteriors.is_empty());
    }

    // -----------------------------------------------------------------
    // PosteriorRecord *_row stride correctness
    // -----------------------------------------------------------------

    #[test]
    fn posteriors_row_returns_correct_slice_for_last_sample_index() {
        let record = merged_record_simple(
            1,
            100,
            vec![b"A", b"C"],
            2,
            vec![
                vec![0.0, -50.0, -50.0],
                vec![-50.0, 0.0, -50.0],
                vec![-50.0, -50.0, 0.0],
            ],
        );
        let pr = single_ok(record);
        let last = pr.posteriors_row(pr.n_samples - 1);
        assert_eq!(last.len(), 3);
        assert!(last[2] > 0.99, "last row = {last:?}");
    }

    #[test]
    fn scalars_row_uses_n_alleles_stride_not_n_genotypes() {
        let record = merged_record_simple(
            1,
            100,
            vec![b"A", b"C", b"G"],
            2,
            vec![vec![0.0; 6], vec![0.0; 6]],
        );
        let pr = single_ok(record);
        assert_eq!(pr.alleles.len(), 3);
        let row = pr.scalars_row(1);
        assert_eq!(row.len(), 3, "scalars_row stride should equal n_alleles");
    }

    #[test]
    fn chain_anchor_flags_row_returns_n_alleles_long_slice() {
        let record = merged_record_simple(1, 100, vec![b"A", b"C", b"G"], 2, vec![vec![0.0; 6]]);
        let pr = single_ok(record);
        assert_eq!(pr.chain_anchor_flags_row(0).len(), pr.alleles.len());
    }

    // -----------------------------------------------------------------
    // Engine accessors and Debug
    // -----------------------------------------------------------------

    #[test]
    fn engine_config_accessor_returns_the_value_used_at_construction() {
        let upstream = std::iter::empty::<Result<MergedRecord, PerGroupMergerError>>();
        let config = PosteriorEngineConfig {
            max_iterations: 7,
            ..Default::default()
        };
        let eng = PosteriorEngine::with_config(upstream, config);
        assert_eq!(eng.config().max_iterations, 7);
    }

    #[test]
    fn posterior_engine_debug_includes_config_and_is_latched_fields() {
        let upstream = std::iter::empty::<Result<MergedRecord, PerGroupMergerError>>();
        let eng = PosteriorEngine::new(upstream);
        let dbg = format!("{eng:?}");
        assert!(dbg.contains("PosteriorEngine"));
        assert!(dbg.contains("config"));
        assert!(dbg.contains("is_latched"));
    }

    // -----------------------------------------------------------------
    // Numeric-helper unit tests
    // -----------------------------------------------------------------

    #[test]
    fn log_sum_exp_2_handles_neg_infinity_arguments() {
        let m = ExactMath;
        assert_eq!(log_sum_exp_2(&m, f64::NEG_INFINITY, 3.0), 3.0);
        assert_eq!(log_sum_exp_2(&m, 3.0, f64::NEG_INFINITY), 3.0);
        assert!(approx(log_sum_exp_2(&m, 0.0, 0.0), 2.0_f64.ln(), 1e-12));
    }

    #[test]
    fn log_sum_exp_slice_returns_neg_infinity_when_every_entry_is_neg_infinity() {
        let v = vec![f64::NEG_INFINITY; 5];
        assert_eq!(log_sum_exp_slice(&ExactMath, &v), f64::NEG_INFINITY);
    }

    #[test]
    fn log_sum_exp_slice_returns_neg_infinity_for_empty_slice() {
        assert_eq!(log_sum_exp_slice(&ExactMath, &[]), f64::NEG_INFINITY);
    }

    #[test]
    fn log_sum_exp_slice_matches_log_sum_exp_2_on_two_element_inputs() {
        let m = ExactMath;
        for &(a, b) in &[(0.0_f64, 1.0), (-3.0, 7.0), (-5.5, 5.5)] {
            assert!(approx(
                log_sum_exp_slice(&m, &[a, b]),
                log_sum_exp_2(&m, a, b),
                1e-12
            ));
        }
    }

    // `log_factorial`, `log_multinomial_coefficient`, and
    // `homozygous_allele` were moved alongside the genotype-shape cache
    // to `posterior_engine::shape`; their tests live in that module.

    #[test]
    fn classify_alleles_routes_pseudocounts_correctly() {
        let alleles = simple_alleles(&[b"A", b"C", b"AT"]);
        let classes = classify_alleles(&alleles);
        assert_eq!(
            classes,
            vec![AlleleClass::Ref, AlleleClass::SnpAlt, AlleleClass::IndelAlt]
        );
        let mut compound_alleles = simple_alleles(&[b"A", b"C"]);
        compound_alleles[1].is_compound = true;
        let classes = classify_alleles(&compound_alleles);
        assert_eq!(classes, vec![AlleleClass::Ref, AlleleClass::CompoundAlt]);
    }

    #[test]
    fn pseudocount_for_returns_class_specific_value_from_config() {
        let config = PosteriorEngineConfig {
            ref_pseudocount: 1.0,
            snp_alt_pseudocount: 2.0,
            indel_alt_pseudocount: 3.0,
            compound_alt_pseudocount: 4.0,
            ..Default::default()
        };
        assert_eq!(pseudocount_for(AlleleClass::Ref, &config), 1.0);
        assert_eq!(pseudocount_for(AlleleClass::SnpAlt, &config), 2.0);
        assert_eq!(pseudocount_for(AlleleClass::IndelAlt, &config), 3.0);
        assert_eq!(pseudocount_for(AlleleClass::CompoundAlt, &config), 4.0);
    }

    #[test]
    fn safe_ln_returns_neg_infinity_for_zero_and_negative_inputs() {
        let m = ExactMath;
        assert_eq!(safe_ln(&m, 0.0), f64::NEG_INFINITY);
        assert_eq!(safe_ln(&m, -0.0), f64::NEG_INFINITY);
        assert_eq!(safe_ln(&m, -1.0), f64::NEG_INFINITY);
        assert!(approx(safe_ln(&m, std::f64::consts::E), 1.0, 1e-12));
        assert!(safe_ln(&m, f64::NAN).is_nan());
    }

    #[test]
    fn max_abs_diff_returns_largest_absolute_difference() {
        assert_eq!(max_abs_diff(&[0.0, 0.0], &[0.0, 0.0]), 0.0);
        assert_eq!(max_abs_diff(&[1.0, 2.0], &[1.0, 2.5]), 0.5);
        assert_eq!(max_abs_diff(&[1.0], &[-1.0]), 2.0);
        assert_eq!(max_abs_diff(&[-1.0], &[1.0]), 2.0);
    }

    #[test]
    fn classify_nonfinite_distinguishes_nan_pos_inf_neg_inf() {
        assert_eq!(classify_nonfinite(f64::NAN), NonFiniteKind::NaN);
        assert_eq!(
            classify_nonfinite(f64::INFINITY),
            NonFiniteKind::PositiveInfinity
        );
        assert_eq!(
            classify_nonfinite(f64::NEG_INFINITY),
            NonFiniteKind::NegativeInfinity
        );
    }

    // -----------------------------------------------------------------
    // Math-backend approximate-parity harness
    //
    // For every fixture in the battery we run both `ExactMath` (the
    // bit-identical baseline) and `InterpUnivariateMath` (the approx
    // backend), then compare outputs field-by-field. Per-cell errors
    // are accumulated into distributions so the test asserts on max
    // *and* p99 (catches a backend that meets the max but has a fat
    // tail). `argmax_mismatch_above_margin` counts samples whose best
    // genotype differs *and* whose top posterior margin is above the
    // "borderline" threshold — i.e. cases where a near-tie is not the
    // excuse. Tolerances follow §"Approximate parity" of the
    // implementation plan.
    // -----------------------------------------------------------------

    use super::backends::{ExactMath, InterpUnivariateMath, MathBackend};

    const POSTERIOR_MAX_TOL: f64 = 1e-4;
    const POSTERIOR_P99_TOL: f64 = 5e-5;
    const PHRED_MAX_TOL: f64 = 1.0;
    const PHRED_P99_TOL: f64 = 0.5;
    const ALLELE_FREQ_MAX_TOL: f64 = 1e-4;
    const BORDERLINE_MARGIN: f64 = 0.05;

    /// Per-metric error summary. Fields are read via the `Debug`
    /// derive when the harness prints the report on assertion failure.
    #[derive(Debug)]
    #[allow(dead_code)]
    struct ErrorReport {
        label: &'static str,
        n: usize,
        max: f64,
        p99: f64,
        p95: f64,
        p50: f64,
        mean: f64,
        argmax_fixture: &'static str,
    }

    impl ErrorReport {
        fn from_samples(label: &'static str, mut samples: Vec<(f64, &'static str)>) -> Self {
            samples.sort_by(|a, b| {
                a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)
            });
            let n = samples.len();
            let pick = |q: f64| -> (f64, &'static str) {
                if n == 0 {
                    (0.0, "")
                } else {
                    let idx = ((n as f64) * q).floor() as usize;
                    samples[idx.min(n - 1)]
                }
            };
            let (max, argmax_fixture) = samples.last().copied().unwrap_or((0.0, ""));
            let (p99, _) = pick(0.99);
            let (p95, _) = pick(0.95);
            let (p50, _) = pick(0.50);
            let mean = if n == 0 {
                0.0
            } else {
                samples.iter().map(|s| s.0).sum::<f64>() / n as f64
            };
            ErrorReport {
                label,
                n,
                max,
                p99,
                p95,
                p50,
                mean,
                argmax_fixture,
            }
        }
    }

    #[derive(Default)]
    struct AccuracyAccumulator {
        posterior_errors: Vec<(f64, &'static str)>,
        phred_errors: Vec<(f64, &'static str)>,
        allele_freq_errors: Vec<(f64, &'static str)>,
        argmax_mismatch_above_margin: usize,
        argmax_mismatch_borderline_skipped: usize,
        fixtures_total: usize,
    }

    impl AccuracyAccumulator {
        fn ingest(
            &mut self,
            fixture: &'static str,
            exact: &PosteriorRecord,
            approx: &PosteriorRecord,
        ) {
            self.fixtures_total += 1;
            assert_eq!(exact.n_samples, approx.n_samples);
            assert_eq!(exact.n_genotypes, approx.n_genotypes);
            assert_eq!(exact.allele_frequencies.len(), approx.allele_frequencies.len());

            for (e, a) in exact
                .posteriors
                .iter()
                .zip(approx.posteriors.iter())
            {
                self.posterior_errors.push(((e - a).abs(), fixture));
            }
            for (e, a) in exact
                .allele_frequencies
                .iter()
                .zip(approx.allele_frequencies.iter())
            {
                self.allele_freq_errors.push(((e - a).abs(), fixture));
            }
            for (e, a) in exact.gq_phred.iter().zip(approx.gq_phred.iter()) {
                if e.is_finite() && a.is_finite() {
                    self.phred_errors.push(((e - a).abs(), fixture));
                }
            }
            if exact.qual_phred.is_finite() && approx.qual_phred.is_finite() {
                self.phred_errors
                    .push(((exact.qual_phred - approx.qual_phred).abs(), fixture));
            }

            // `best_genotype` is checked per sample only when the
            // exact top posterior dominates by ≥ BORDERLINE_MARGIN.
            // Closer races are tagged "borderline" and excluded — the
            // approximate posterior may flip a near-tied argmax without
            // implying a real divergence.
            for sample_idx in 0..exact.n_samples {
                let row = exact.posteriors_row(sample_idx);
                let mut sorted = row.to_vec();
                sorted.sort_by(|a, b| {
                    b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal)
                });
                let margin = if sorted.len() >= 2 {
                    sorted[0] - sorted[1]
                } else {
                    1.0
                };
                if margin < BORDERLINE_MARGIN {
                    self.argmax_mismatch_borderline_skipped += 1;
                    continue;
                }
                if exact.best_genotype[sample_idx] != approx.best_genotype[sample_idx] {
                    self.argmax_mismatch_above_margin += 1;
                }
            }
        }

        fn summarise(self) -> AccuracyReport {
            AccuracyReport {
                posterior_error: ErrorReport::from_samples("posterior", self.posterior_errors),
                phred_error: ErrorReport::from_samples("phred", self.phred_errors),
                allele_freq_error: ErrorReport::from_samples(
                    "allele_freq",
                    self.allele_freq_errors,
                ),
                argmax_mismatch_above_margin: self.argmax_mismatch_above_margin,
                argmax_mismatch_borderline_skipped: self.argmax_mismatch_borderline_skipped,
                fixtures_total: self.fixtures_total,
            }
        }
    }

    #[derive(Debug)]
    #[allow(dead_code)]
    struct AccuracyReport {
        posterior_error: ErrorReport,
        phred_error: ErrorReport,
        allele_freq_error: ErrorReport,
        argmax_mismatch_above_margin: usize,
        argmax_mismatch_borderline_skipped: usize,
        fixtures_total: usize,
    }

    /// Run a single fixture under both backends and ingest the per-cell
    /// errors into `accum`.
    fn drive_pair<M: MathBackend + Copy>(
        accum: &mut AccuracyAccumulator,
        fixture: &'static str,
        record: MergedRecord,
        approx_backend: M,
        config: Option<PosteriorEngineConfig>,
    ) {
        let cfg = config.unwrap_or_else(PosteriorEngineConfig::with_project_defaults);
        // The "exact" baseline is pinned to `ExactMath` regardless of
        // the engine's default backend — the harness's contract is "is
        // `approx_backend` close to the bit-identical reference?".
        let exact_engine = PosteriorEngine::with_math_backend(
            std::iter::once(Ok::<_, PerGroupMergerError>(record.clone())),
            cfg.clone(),
            ExactMath,
        );
        let approx_engine = PosteriorEngine::with_math_backend(
            std::iter::once(Ok::<_, PerGroupMergerError>(record)),
            cfg,
            approx_backend,
        );
        let exact: Vec<_> = exact_engine
            .collect::<Result<Vec<_>, _>>()
            .expect("ExactMath fixture errored");
        let approx: Vec<_> = approx_engine
            .collect::<Result<Vec<_>, _>>()
            .expect("approximate backend fixture errored");
        assert_eq!(exact.len(), 1);
        assert_eq!(approx.len(), 1);
        accum.ingest(fixture, &exact[0], &approx[0]);
    }

    fn cohort_likelihoods_biallelic_diploid(n_samples: usize) -> Vec<Vec<f64>> {
        // Half ALT-leaning, half REF-leaning. Likelihoods chosen so
        // every sample's argmax margin sits well above the borderline.
        (0..n_samples)
            .map(|s: usize| {
                if s.is_multiple_of(2) {
                    vec![-30.0, -10.0, -0.1]
                } else {
                    vec![-0.1, -10.0, -30.0]
                }
            })
            .collect()
    }

    fn build_accuracy_battery<M: MathBackend + Copy>(approx: M) -> AccuracyReport {
        let mut accum = AccuracyAccumulator::default();

        // Biallelic, single sample, strong-REF.
        drive_pair(
            &mut accum,
            "biallelic_single_sample_strong_ref",
            merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![-0.1, -10.0, -30.0]]),
            approx,
            None,
        );

        // Biallelic, single sample, strong-ALT.
        drive_pair(
            &mut accum,
            "biallelic_single_sample_strong_alt",
            merged_record_simple(1, 200, vec![b"A", b"C"], 2, vec![vec![-30.0, -10.0, -0.1]]),
            approx,
            None,
        );

        // Biallelic, 64 samples, mixed evidence.
        drive_pair(
            &mut accum,
            "biallelic_64_samples_mixed",
            merged_record_simple(
                1,
                300,
                vec![b"A", b"C"],
                2,
                cohort_likelihoods_biallelic_diploid(64),
            ),
            approx,
            None,
        );

        // Triallelic SNP, 8 samples.
        let trial_likelihoods: Vec<Vec<f64>> = (0..8)
            .map(|s| {
                // 6 genotypes for n_alleles=3, ploidy=2: AA, AB, BB, AC, BC, CC
                let mut row = vec![-20.0; 6];
                let idx = s % 6;
                row[idx] = -0.1;
                row
            })
            .collect();
        drive_pair(
            &mut accum,
            "triallelic_8_samples",
            merged_record_simple(
                1,
                400,
                vec![b"A", b"C", b"G"],
                2,
                trial_likelihoods,
            ),
            approx,
            None,
        );

        // Compound allele, chain-anchored, 8 samples.
        let compound_ll: Vec<Vec<f64>> = (0..8)
            .map(|s: usize| {
                if s.is_multiple_of(2) {
                    vec![-0.1, -10.0, -30.0]
                } else {
                    vec![-30.0, -5.0, -0.1]
                }
            })
            .collect();
        drive_pair(
            &mut accum,
            "compound_chain_anchored",
            merged_record_with_compound(
                1,
                500,
                vec![b"A", b"CT"],
                &[1],
                vec![vec![false, true]; 8],
                2,
                compound_ll,
            ),
            approx,
            None,
        );

        // Compound allele with a chain-broken sample.
        let compound_broken_ll: Vec<Vec<f64>> = (0..6)
            .map(|s: usize| {
                if s == 0 {
                    vec![-0.5, -2.0, -5.0]
                } else if s.is_multiple_of(2) {
                    vec![-0.1, -10.0, -30.0]
                } else {
                    vec![-30.0, -5.0, -0.1]
                }
            })
            .collect();
        let mut chain_flags = vec![vec![false, true]; 6];
        chain_flags[0] = vec![false, false];
        drive_pair(
            &mut accum,
            "compound_chain_broken_sample",
            merged_record_with_compound(
                1,
                600,
                vec![b"A", b"CT"],
                &[1],
                chain_flags,
                2,
                compound_broken_ll,
            ),
            approx,
            None,
        );

        // Non-zero fixation index default.
        let mut config = PosteriorEngineConfig::with_project_defaults();
        config.fixation_index_default = 0.3;
        drive_pair(
            &mut accum,
            "fixation_index_default_0_3",
            merged_record_simple(
                1,
                700,
                vec![b"A", b"C"],
                2,
                cohort_likelihoods_biallelic_diploid(32),
            ),
            approx,
            Some(config),
        );

        // Per-sample fixation index overrides.
        let mut config = PosteriorEngineConfig::with_project_defaults();
        config.fixation_index_overrides =
            Some((0..32).map(|s: usize| if s.is_multiple_of(2) { 0.0 } else { 0.4 }).collect());
        drive_pair(
            &mut accum,
            "fixation_index_per_sample_overrides",
            merged_record_simple(
                1,
                800,
                vec![b"A", b"C"],
                2,
                cohort_likelihoods_biallelic_diploid(32),
            ),
            approx,
            Some(config),
        );

        // Triploid (polyploid path).
        let triploid_ll: Vec<Vec<f64>> = (0..16)
            .map(|s| {
                // 4 genotypes for ploidy=3, n_alleles=2: 0/0/0, 0/0/1, 0/1/1, 1/1/1
                let mut row = vec![-20.0_f64; 4];
                row[s % 4] = -0.1;
                row
            })
            .collect();
        drive_pair(
            &mut accum,
            "triploid_biallelic",
            merged_record_simple(1, 900, vec![b"A", b"C"], 3, triploid_ll),
            approx,
            None,
        );

        // Contamination on, all samples c_s > 0.
        let mut contam_cfg = PosteriorEngineConfig::with_project_defaults();
        contam_cfg.contamination = Some(
            crate::var_calling::contamination_estimation::ContaminationEstimates::from_user_supplied(
                vec![Some(0.03_f64); 32],
                vec![[0.6, 0.3, 0.1]],
                vec![0; 32],
            )
            .expect("valid contamination fixture"),
        );
        drive_pair(
            &mut accum,
            "contam_on_all_samples",
            merged_record_simple(
                1,
                1000,
                vec![b"A", b"C"],
                2,
                cohort_likelihoods_biallelic_diploid(32),
            ),
            approx,
            Some(contam_cfg),
        );

        // Contamination on, mixed c_s = None and c_s > 0.
        let mut mixed_cfg = PosteriorEngineConfig::with_project_defaults();
        let c_s_mixed: Vec<Option<f64>> = (0..32)
            .map(|s: usize| if s.is_multiple_of(3) { None } else { Some(0.05_f64) })
            .collect();
        mixed_cfg.contamination = Some(
            crate::var_calling::contamination_estimation::ContaminationEstimates::from_user_supplied(
                c_s_mixed,
                vec![[0.6, 0.3, 0.1]],
                vec![0; 32],
            )
            .expect("valid mixed contamination"),
        );
        drive_pair(
            &mut accum,
            "contam_on_mixed_c_s",
            merged_record_simple(
                1,
                1100,
                vec![b"A", b"C"],
                2,
                cohort_likelihoods_biallelic_diploid(32),
            ),
            approx,
            Some(mixed_cfg),
        );

        // Contamination on with triallelic site.
        let mut tri_contam_cfg = PosteriorEngineConfig::with_project_defaults();
        tri_contam_cfg.contamination = Some(
            crate::var_calling::contamination_estimation::ContaminationEstimates::from_user_supplied(
                vec![Some(0.03_f64); 8],
                vec![[0.6, 0.3, 0.1]],
                vec![0; 8],
            )
            .expect("valid triallelic contam fixture"),
        );
        let tri_contam_ll: Vec<Vec<f64>> = (0..8)
            .map(|s| {
                let mut row = vec![-20.0_f64; 6];
                row[s % 6] = -0.1;
                row
            })
            .collect();
        drive_pair(
            &mut accum,
            "contam_on_triallelic",
            merged_record_simple(1, 1200, vec![b"A", b"C", b"G"], 2, tri_contam_ll),
            approx,
            Some(tri_contam_cfg),
        );

        accum.summarise()
    }

    #[test]
    fn interp_univariate_clears_approximate_parity_budget() {
        let report = build_accuracy_battery(InterpUnivariateMath);
        eprintln!("InterpUnivariateMath accuracy report:\n{report:#?}");

        assert!(
            report.posterior_error.max <= POSTERIOR_MAX_TOL,
            "posterior_error.max = {} > {}",
            report.posterior_error.max,
            POSTERIOR_MAX_TOL,
        );
        assert!(
            report.posterior_error.p99 <= POSTERIOR_P99_TOL,
            "posterior_error.p99 = {} > {}",
            report.posterior_error.p99,
            POSTERIOR_P99_TOL,
        );
        assert!(
            report.phred_error.max <= PHRED_MAX_TOL,
            "phred_error.max = {} > {}",
            report.phred_error.max,
            PHRED_MAX_TOL,
        );
        assert!(
            report.phred_error.p99 <= PHRED_P99_TOL,
            "phred_error.p99 = {} > {}",
            report.phred_error.p99,
            PHRED_P99_TOL,
        );
        assert!(
            report.allele_freq_error.max <= ALLELE_FREQ_MAX_TOL,
            "allele_freq_error.max = {} > {}",
            report.allele_freq_error.max,
            ALLELE_FREQ_MAX_TOL,
        );
        assert_eq!(
            report.argmax_mismatch_above_margin, 0,
            "{} samples picked a different best genotype with margin >= {}",
            report.argmax_mismatch_above_margin, BORDERLINE_MARGIN,
        );
    }

    /// Sanity check: `ExactMath` against itself produces all-zero
    /// errors. Guards the harness machinery — if this fails, the
    /// comparison logic is broken, not the backends.
    #[test]
    fn exact_against_exact_yields_zero_errors() {
        let report = build_accuracy_battery(ExactMath);
        assert_eq!(report.posterior_error.max, 0.0);
        assert_eq!(report.phred_error.max, 0.0);
        assert_eq!(report.allele_freq_error.max, 0.0);
        assert_eq!(report.argmax_mismatch_above_margin, 0);
    }

    // -----------------------------------------------------------------
    // Property tests (proptest)
    // -----------------------------------------------------------------

    use proptest::prelude::*;

    /// (ploidy, n_alleles, n_samples) → a `(ploidy, alleles, likelihoods)`
    /// tuple suitable for `merged_record_simple`. Walks the
    /// (ploidy ∈ {2..=4}, n_alleles ∈ {2..=3}, n_samples ∈ {1..=4})
    /// matrix so each property holds at every supported polyploid /
    /// multiallelic shape, not just diploid biallelic.
    fn proptest_record_strategy() -> impl Strategy<Value = (u8, Vec<&'static [u8]>, Vec<Vec<f64>>)>
    {
        (2u8..=4u8, 2usize..=3usize, 1usize..=4usize).prop_flat_map(
            |(ploidy, n_alleles, n_samples)| {
                static ACGT: &[&[u8]] = &[b"A", b"C", b"G", b"T"];
                let alleles: Vec<&'static [u8]> = ACGT.iter().take(n_alleles).copied().collect();
                let n_genotypes = genotype_order(ploidy, n_alleles).len();
                let row = prop::collection::vec(-50.0_f64..0.0_f64, n_genotypes);
                let likelihoods = prop::collection::vec(row, n_samples);
                (Just(ploidy), Just(alleles), likelihoods)
            },
        )
    }

    proptest! {
        #![proptest_config(ProptestConfig { cases: 64, ..ProptestConfig::default() })]

        #[test]
        fn posteriors_sum_to_one_at_any_ploidy_and_n_alleles(
            (ploidy, alleles, likelihoods) in proptest_record_strategy()
        ) {
            let record = merged_record_simple(1, 100, alleles, ploidy, likelihoods);
            // M3: skip only the specific `DidNotConverge` outcome (a
            // documented engine error on adversarial likelihood
            // matrices); any other error variant is a real bug and
            // should fail the proptest loudly.
            let pr = match try_single(record) {
                Ok(pr) => pr,
                Err(PosteriorEngineError::DidNotConverge { .. }) => return Ok(()),
                Err(other) => return Err(TestCaseError::fail(format!(
                    "unexpected engine error: {other}"
                ))),
            };
            for s in 0..pr.n_samples {
                let sum: f64 = pr.posteriors_row(s).iter().sum();
                prop_assert!((sum - 1.0).abs() < 1e-9, "ploidy={ploidy} sum={sum}");
            }
        }

        #[test]
        fn p_hat_is_a_valid_simplex_at_any_ploidy_and_n_alleles(
            (ploidy, alleles, likelihoods) in proptest_record_strategy()
        ) {
            let record = merged_record_simple(1, 100, alleles, ploidy, likelihoods);
            let pr = match try_single(record) {
                Ok(pr) => pr,
                Err(PosteriorEngineError::DidNotConverge { .. }) => return Ok(()),
                Err(other) => return Err(TestCaseError::fail(format!(
                    "unexpected engine error: {other}"
                ))),
            };
            for p in &pr.allele_frequencies {
                prop_assert!(*p >= 0.0 && *p <= 1.0);
            }
            let sum: f64 = pr.allele_frequencies.iter().sum();
            prop_assert!((sum - 1.0).abs() < 1e-9);
        }

        #[test]
        fn sample_permutation_preserves_p_hat_and_qual_at_any_ploidy_and_n_alleles(
            (ploidy, alleles, likelihoods) in proptest_record_strategy(),
            seed in any::<u64>(),
        ) {
            // Determined permutation of the sample order.
            let mut perm: Vec<usize> = (0..likelihoods.len()).collect();
            let mut s = seed | 1;
            for i in (1..perm.len()).rev() {
                s ^= s << 13; s ^= s >> 7; s ^= s << 17;
                let j = (s as usize) % (i + 1);
                perm.swap(i, j);
            }
            let permuted: Vec<Vec<f64>> = perm.iter().map(|&i| likelihoods[i].clone()).collect();
            let record_a = merged_record_simple(1, 100, alleles.clone(), ploidy, likelihoods);
            let record_b = merged_record_simple(1, 100, alleles, ploidy, permuted);
            // Both records must converge for the invariance assertion
            // to be meaningful; skip only on DidNotConverge.
            let pr_a = match try_single(record_a) {
                Ok(pr) => pr,
                Err(PosteriorEngineError::DidNotConverge { .. }) => return Ok(()),
                Err(other) => return Err(TestCaseError::fail(format!(
                    "unexpected engine error on record_a: {other}"
                ))),
            };
            let pr_b = match try_single(record_b) {
                Ok(pr) => pr,
                Err(PosteriorEngineError::DidNotConverge { .. }) => return Ok(()),
                Err(other) => return Err(TestCaseError::fail(format!(
                    "unexpected engine error on record_b: {other}"
                ))),
            };
            for (a, b) in pr_a.allele_frequencies.iter().zip(pr_b.allele_frequencies.iter()) {
                prop_assert!((a - b).abs() < 1e-6);
            }
            prop_assert!((pr_a.qual_phred - pr_b.qual_phred).abs() < 1e-6);
        }

        #[test]
        fn larger_ref_pseudocount_cannot_increase_p_alt(
            likelihoods in prop::collection::vec(
                prop::collection::vec(-50.0_f64..0.0_f64, 3..=3),
                1..6,
            )
        ) {
            let record_a = merged_record_simple(
                1, 100, vec![b"A", b"C"], 2, likelihoods.clone(),
            );
            let pr_a = single_ok(record_a);
            let config = PosteriorEngineConfig {
                ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT * 10.0,
                ..Default::default()
            };
            let record_b = merged_record_simple(1, 100, vec![b"A", b"C"], 2, likelihoods);
            let out_b = engine_for_with_config(record_b, config);
            let pr_b = out_b.into_iter().next().unwrap().expect("posterior");
            prop_assert!(pr_b.allele_frequencies[1] <= pr_a.allele_frequencies[1] + 1e-9);
        }
    }
}
