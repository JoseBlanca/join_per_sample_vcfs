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
//! 4. Convergence check on `max |Δp̂|`. Records that satisfy the
//!    threshold emit with `EmDiagnostics::converged = true`; records
//!    that hit `config.max_iterations` without satisfying it emit
//!    with `EmDiagnostics::converged = false` and are routed by the
//!    VCF writer to `FILTER=EMNoConv`. The
//!    [`PosteriorEngineError::DidNotConverge`] variant is reserved
//!    for future engine variants but is no longer raised by the
//!    exact-math EM loop — emit-with-flag preserves the record so a
//!    single hard site doesn't kill the whole cohort run.
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

use self::backends::{InterpUnivariateSimdMath, MathBackend};
use self::shape::{GenotypeShape, shape_for};
use crate::pileup_record::AlleleSupportStats;
use crate::var_calling::contamination_estimation::{
    AlleleClass as ContamAlleleClass, ContaminationEstimates, MAX_BASE_ERROR, MIN_BASE_ERROR,
};
use crate::var_calling::per_group_merger::{
    MAX_PLOIDY, MergedAllele, MergedRecord, PerGroupMergerError, genotype_order,
};

/// Default EM convergence threshold on `max_a |p̂_new[a] − p̂_old[a]|`.
///
/// Source: implementation-plan `posterior_engine.md` §"Step 5 —
/// convergence check" originally picked `1e-4`; GATK's `0.1` on raw
/// allele counts is roughly equivalent at depth ≈ 20.
///
/// Relaxed from `1e-4` to `1e-3` (2026-05-20) after tomato (SL4.0)
/// SRR7279725_small.psp surfaced hard records where the EM plateau'd
/// in the `5e-4 … 1e-3` band and hit
/// [`DEFAULT_MAX_ITERATIONS`] without crossing the original threshold
/// (e.g. SL4.0ch01:2025118, `last_delta = 7.33e-4` at iter 50).
/// `1e-3` is in the typical range for population-genetics EM and the
/// well-conditioned records that previously converged below `1e-4`
/// still do so well before the iteration cap — the change buys
/// robustness on hard records without giving up precision on easy
/// ones.
pub const DEFAULT_CONVERGENCE_THRESHOLD: f64 = 1e-3;

/// Default hard cap on per-record EM iterations. Records that exceed
/// it emit with `EmDiagnostics::converged = false` (routed by the VCF
/// writer to `FILTER=EMNoConv`) rather than hard-failing the cohort
/// run. Empirically the closed-form M-steps converge in 3–5 rounds
/// per the GATK reference; 50 is the belt-and-braces ceiling.
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

/// Default Phred cap on per-sample GQ. Matches the GATK and bcftools
/// convention (`GQ` capped at 99) so downstream tooling sees a
/// familiar range; the cap also prevents `+∞` GQ when EM yields
/// `P(best) = 1` exactly. The CLI exposes this as `--max-gq-phred`;
/// the actual validator upper bound is [`GQ_PHRED_RANGE_MAX`] (200.0).
pub const DEFAULT_MAX_GQ_PHRED: f64 = 99.0;

/// Validation upper bound on `convergence_threshold` (cohort CLI plan).
/// Loose-but-not-degenerate: a `1.0` threshold would always exit after
/// one EM iteration (max possible per-allele change is `1.0`).
pub const CONVERGENCE_THRESHOLD_RANGE_MAX: f64 = 0.1;

/// Validation upper bound on `max_iterations` (cohort CLI plan). EM
/// converges in 3–5 iterations on the GATK reference; `500` is 10× the
/// default of `50`, generous headroom for debugging without entering
/// runaway territory.
pub const MAX_ITERATIONS_RANGE_MAX: u32 = 500;

/// Validation upper bound on Dirichlet pseudocounts. Beyond `1000` the
/// prior swamps any realistic cohort-level evidence (forces the prior
/// unconditionally).
pub const PSEUDOCOUNT_RANGE_MAX: f64 = 1000.0;

/// Validation lower bound on `max_gq_phred` (exclusive). Phred values
/// below `10` are absurdly low confidence caps.
pub const GQ_PHRED_RANGE_MIN_EXCLUSIVE: f64 = 10.0;

/// Validation upper bound on `max_gq_phred`. Phred `200` corresponds to
/// `p = 1e-20` — already absurdly confident; anything higher has no
/// physical meaning.
pub const GQ_PHRED_RANGE_MAX: f64 = 200.0;

/// Phred conversion factor — `Phred = -10 · log10(p)`. The literal
/// `-10.0` appears in every Phred-conversion site; the named const
/// keeps the meaning in one place.
const PHRED_SCALE: f64 = -10.0;

/// Index of the hom-REF genotype in canonical [`genotype_order`].
/// REF/REF/.../REF is the all-zero allele tuple and lands at position
/// 0 for every (ploidy, n_alleles).
const HOM_REF_GENOTYPE_IDX: usize = 0;

/// Tunable knobs for the posterior engine.
///
/// **Construction.** Build via [`PosteriorEngineConfig::new`] (or
/// [`Default::default()`], or
/// [`with_project_defaults`](Self::with_project_defaults) — all
/// three yield the project defaults) then chain `with_*` setters
/// to override individual knobs:
///
/// ```ignore
/// let cfg = PosteriorEngineConfig::new()
///     .with_max_iterations(80)?
///     .with_max_gq_phred(60.0)?
///     .with_contamination(Some(estimates))?;
/// ```
///
/// Every `with_*` setter returns `Result<Self, PosteriorEngineConfigError>`
/// and validates the field's invariants on the way in. The
/// `contamination` field is **private** — it can only be set via
/// [`with_contamination`](Self::with_contamination), so any future
/// cross-field invariant added to that setter cannot be bypassed
/// by direct field assignment from outside the module.
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
    /// GQ Phred cap. Defaults to [`DEFAULT_MAX_GQ_PHRED`].
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
    /// **Private.** Set via
    /// [`with_contamination`](Self::with_contamination); the setter
    /// runs whatever cross-field validation lands in future revisions
    /// of this engine. Direct mutation from outside the module would
    /// bypass that gate (this is the M4 fix from the 2026-05-19
    /// cohort CLI review).
    ///
    /// [`estimate_contamination`]: crate::var_calling::contamination_estimation::estimate_contamination
    contamination: Option<ContaminationEstimates>,

    /// Precomputed per-sample biallelic-diploid genotype log-priors from the
    /// SFS-marginalized prior (`[ln P(0/0), ln P(0/1), ln P(1/1)]` per sample,
    /// cohort-length). `Some` → the engine uses this SFS-marginal prior in place
    /// of the HWE(p̂) + Dirichlet genotype prior for **biallelic-diploid**
    /// records (multiallelic / non-diploid records keep the HWE path — the SFS
    /// prior's general-shape support is deferred). `None` → the HWE(p̂) prior for
    /// every shape (the pre-SFS behaviour). Built once by the driver from the
    /// cohort diversity θ and per-sample inbreeding F; see
    /// [`crate::var_calling::sfs_prior`].
    ///
    /// **Private.** Set via
    /// [`with_sfs_prior_tables`](Self::with_sfs_prior_tables).
    sfs_prior_tables: Option<Vec<[f64; 3]>>,

    /// Cohort nucleotide diversity `θ̂` for the general SFS genotype prior. `Some`
    /// → every record's genotype prior is the **Dirichlet-multinomial** with
    /// concentration `α` from θ̂ ([`crate::genetics::alpha_from_diversity`]),
    /// replacing the HWE(p̂) + Dirichlet plug-in for **all** `(ploidy, n_alleles)`
    /// shapes. `None` → the plug-in HWE(p̂) prior (the pre-SFS behaviour). This is
    /// the general-shape successor to `sfs_prior_tables` (which covered only
    /// biallelic-diploid). If both are set (transitional only — G5 removes
    /// `sfs_prior_tables`), the legacy table path still takes precedence for
    /// biallelic-diploid records via the E-step dispatch, and the general
    /// Dirichlet-multinomial covers every other shape. The driver sets exactly
    /// one.
    ///
    /// **Private.** Set via
    /// [`with_nucleotide_diversity`](Self::with_nucleotide_diversity).
    nucleotide_diversity: Option<f64>,
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
    /// | `convergence_threshold`            | [`DEFAULT_CONVERGENCE_THRESHOLD`]     | 1e-3     |
    /// | `max_iterations`                   | [`DEFAULT_MAX_ITERATIONS`]            | 50       |
    /// | `ref_pseudocount`                  | [`DEFAULT_REF_PSEUDOCOUNT`]           | 10.0     |
    /// | `snp_alt_pseudocount`              | [`DEFAULT_SNP_ALT_PSEUDOCOUNT`]       | 0.01     |
    /// | `indel_alt_pseudocount`            | [`DEFAULT_INDEL_ALT_PSEUDOCOUNT`]     | 0.00125  |
    /// | `compound_alt_pseudocount`         | [`DEFAULT_COMPOUND_ALT_PSEUDOCOUNT`]  | 0.001    |
    /// | `fixation_index_default`           | [`DEFAULT_INBREEDING_COEFFICIENT`]    | 0.0      |
    /// | `fixation_index_overrides`         | —                                     | `None`   |
    /// | `max_gq_phred`                     | [`DEFAULT_MAX_GQ_PHRED`]                      | 99.0     |
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
            max_gq_phred: DEFAULT_MAX_GQ_PHRED,
            approximate_posterior_calculation: false,
            contamination: None,
            sfs_prior_tables: None,
            nucleotide_diversity: None,
        }
    }

    /// Construct with the project defaults; equivalent to
    /// [`with_project_defaults`](Self::with_project_defaults) and
    /// [`Default::default()`]. Defaults are documented on
    /// `with_project_defaults` and are guaranteed valid by
    /// construction, so [`Self::new`] is infallible.
    ///
    /// Chain `with_*` setters to override individual knobs; each
    /// setter validates and returns `Result<Self, _>`:
    ///
    /// ```ignore
    /// let cfg = PosteriorEngineConfig::new()
    ///     .with_max_iterations(80)?
    ///     .with_max_gq_phred(60.0)?
    ///     .with_contamination(Some(estimates))?;
    /// ```
    pub fn new() -> Self {
        Self::with_project_defaults()
    }

    /// Validating setter for `convergence_threshold` — finite, in
    /// `(0.0, CONVERGENCE_THRESHOLD_RANGE_MAX]` (i.e. `(0.0, 0.1]`).
    pub fn with_convergence_threshold(
        mut self,
        convergence_threshold: f64,
    ) -> Result<Self, PosteriorEngineConfigError> {
        if !(convergence_threshold.is_finite()
            && 0.0 < convergence_threshold
            && convergence_threshold <= CONVERGENCE_THRESHOLD_RANGE_MAX)
        {
            return Err(PosteriorEngineConfigError::InvalidConvergenceThreshold {
                got: convergence_threshold,
            });
        }
        self.convergence_threshold = convergence_threshold;
        Ok(self)
    }

    /// Validating setter for `max_iterations` — in
    /// `1..=MAX_ITERATIONS_RANGE_MAX` (i.e. `1..=500`).
    pub fn with_max_iterations(
        mut self,
        max_iterations: u32,
    ) -> Result<Self, PosteriorEngineConfigError> {
        if !(1..=MAX_ITERATIONS_RANGE_MAX).contains(&max_iterations) {
            return Err(PosteriorEngineConfigError::InvalidMaxIterations {
                got: max_iterations,
            });
        }
        self.max_iterations = max_iterations;
        Ok(self)
    }

    /// Validating setter for `ref_pseudocount` — finite, in
    /// `(0.0, PSEUDOCOUNT_RANGE_MAX]`.
    pub fn with_ref_pseudocount(
        mut self,
        ref_pseudocount: f64,
    ) -> Result<Self, PosteriorEngineConfigError> {
        validate_pseudocount("ref_pseudocount", ref_pseudocount)?;
        self.ref_pseudocount = ref_pseudocount;
        Ok(self)
    }

    /// Validating setter for `snp_alt_pseudocount` — finite, in
    /// `(0.0, PSEUDOCOUNT_RANGE_MAX]`.
    pub fn with_snp_alt_pseudocount(
        mut self,
        snp_alt_pseudocount: f64,
    ) -> Result<Self, PosteriorEngineConfigError> {
        validate_pseudocount("snp_alt_pseudocount", snp_alt_pseudocount)?;
        self.snp_alt_pseudocount = snp_alt_pseudocount;
        Ok(self)
    }

    /// Validating setter for `indel_alt_pseudocount` — finite, in
    /// `(0.0, PSEUDOCOUNT_RANGE_MAX]`.
    pub fn with_indel_alt_pseudocount(
        mut self,
        indel_alt_pseudocount: f64,
    ) -> Result<Self, PosteriorEngineConfigError> {
        validate_pseudocount("indel_alt_pseudocount", indel_alt_pseudocount)?;
        self.indel_alt_pseudocount = indel_alt_pseudocount;
        Ok(self)
    }

    /// Validating setter for `compound_alt_pseudocount` — finite, in
    /// `(0.0, PSEUDOCOUNT_RANGE_MAX]`.
    pub fn with_compound_alt_pseudocount(
        mut self,
        compound_alt_pseudocount: f64,
    ) -> Result<Self, PosteriorEngineConfigError> {
        validate_pseudocount("compound_alt_pseudocount", compound_alt_pseudocount)?;
        self.compound_alt_pseudocount = compound_alt_pseudocount;
        Ok(self)
    }

    /// Validating setter for `fixation_index_default` — finite, in
    /// `[0.0, 1.0]`.
    pub fn with_fixation_index_default(
        mut self,
        fixation_index_default: f64,
    ) -> Result<Self, PosteriorEngineConfigError> {
        if !(fixation_index_default.is_finite() && (0.0..=1.0).contains(&fixation_index_default)) {
            return Err(PosteriorEngineConfigError::InvalidFixationIndex {
                got: fixation_index_default,
            });
        }
        self.fixation_index_default = fixation_index_default;
        Ok(self)
    }

    /// Validating setter for `max_gq_phred` — finite, in
    /// `(GQ_PHRED_RANGE_MIN_EXCLUSIVE, GQ_PHRED_RANGE_MAX]`
    /// (i.e. `(10.0, 200.0]`).
    pub fn with_max_gq_phred(
        mut self,
        max_gq_phred: f64,
    ) -> Result<Self, PosteriorEngineConfigError> {
        if !(max_gq_phred.is_finite()
            && GQ_PHRED_RANGE_MIN_EXCLUSIVE < max_gq_phred
            && max_gq_phred <= GQ_PHRED_RANGE_MAX)
        {
            return Err(PosteriorEngineConfigError::InvalidMaxGqPhred { got: max_gq_phred });
        }
        self.max_gq_phred = max_gq_phred;
        Ok(self)
    }

    /// Validating setter for the private `contamination` field. The
    /// only public path to set the contamination state; closes
    /// **M4** from the 2026-05-19 cohort CLI review by forcing any
    /// future cross-field invariant through this setter rather than
    /// allowing it to be bypassed via direct field assignment.
    /// Currently no cross-field invariants are checked; the setter
    /// is the place to add them.
    pub fn with_contamination(
        mut self,
        contamination: Option<ContaminationEstimates>,
    ) -> Result<Self, PosteriorEngineConfigError> {
        self.contamination = contamination;
        Ok(self)
    }

    /// Validating setter for the private `sfs_prior_tables` field. Each entry
    /// is one sample's biallelic-diploid genotype log-priors
    /// `[ln P(0/0), ln P(0/1), ln P(1/1)]`; `−∞` is allowed (a floored genotype
    /// at full inbreeding) but `NaN`/`+∞` are rejected as a wiring error. The
    /// cohort-length invariant (`len == n_samples`) is checked per record in the
    /// engine, where the sample count is known.
    pub fn with_sfs_prior_tables(
        mut self,
        sfs_prior_tables: Option<Vec<[f64; 3]>>,
    ) -> Result<Self, PosteriorEngineConfigError> {
        if let Some(tables) = &sfs_prior_tables {
            for row in tables {
                for &v in row {
                    if v.is_nan() || v == f64::INFINITY {
                        return Err(PosteriorEngineConfigError::InvalidSfsPriorTable { got: v });
                    }
                }
            }
        }
        self.sfs_prior_tables = sfs_prior_tables;
        Ok(self)
    }

    /// Validating setter for the private `nucleotide_diversity` field (the
    /// general SFS genotype prior's `θ̂`). `Some(θ̂)` must be finite and `>= 0`
    /// (it comes from [`crate::var_calling::diversity::DiversityEstimate`], which
    /// guarantees both). `None` leaves the plug-in HWE(p̂) prior in place.
    pub fn with_nucleotide_diversity(
        mut self,
        nucleotide_diversity: Option<f64>,
    ) -> Result<Self, PosteriorEngineConfigError> {
        if let Some(theta) = nucleotide_diversity
            && !(theta.is_finite() && theta >= 0.0)
        {
            return Err(PosteriorEngineConfigError::InvalidNucleotideDiversity { got: theta });
        }
        self.nucleotide_diversity = nucleotide_diversity;
        Ok(self)
    }
}

/// Shared validator for the four `*_pseudocount` setters. Keeps the
/// "field name + value" error shape consistent across them.
fn validate_pseudocount(field: &'static str, value: f64) -> Result<(), PosteriorEngineConfigError> {
    if !(value.is_finite() && 0.0 < value && value <= PSEUDOCOUNT_RANGE_MAX) {
        return Err(PosteriorEngineConfigError::InvalidPseudocount { field, got: value });
    }
    Ok(())
}

/// Config-construction errors for [`PosteriorEngineConfig::new`].
/// Distinct from [`PosteriorEngineError`] (runtime errors during
/// iteration) — a `PosteriorEngineConfigError` is surfaced at
/// construction time before any upstream records have been pulled.
#[non_exhaustive]
#[derive(Error, Debug, PartialEq)]
pub enum PosteriorEngineConfigError {
    /// `convergence_threshold` was non-finite or outside `(0.0, 0.1]`.
    /// A threshold ≥ 1.0 would always exit after one EM iteration; the
    /// 0.1 upper bound is a strict "10 % per-allele change is loose but
    /// not absurd" line.
    #[error(
        "convergence_threshold must be finite and in (0.0, {CONVERGENCE_THRESHOLD_RANGE_MAX}], got {got}"
    )]
    InvalidConvergenceThreshold { got: f64 },

    /// `max_iterations` was outside `1..=500`. The default is 50; 500
    /// is 10× generous headroom.
    #[error("max_iterations must be in 1..={MAX_ITERATIONS_RANGE_MAX}, got {got}")]
    InvalidMaxIterations { got: u32 },

    /// One of the Dirichlet pseudocounts was non-finite or outside
    /// `(0.0, 1000.0]`. `field` names which one.
    #[error("{field} must be finite and in (0.0, {PSEUDOCOUNT_RANGE_MAX}], got {got}")]
    InvalidPseudocount { field: &'static str, got: f64 },

    /// `fixation_index_default` was non-finite or outside `[0.0, 1.0]`.
    /// `0.0` is outcrossing, `1.0` is full inbreeding.
    #[error("fixation_index_default must be finite and in [0.0, 1.0], got {got}")]
    InvalidFixationIndex { got: f64 },

    /// `max_gq_phred` was non-finite or outside `(10.0, 200.0]`. Phred
    /// 200 already corresponds to `p = 1e-20`; values below 10 are
    /// absurdly low confidence caps.
    #[error(
        "max_gq_phred must be finite and in ({GQ_PHRED_RANGE_MIN_EXCLUSIVE}, {GQ_PHRED_RANGE_MAX}], got {got}"
    )]
    InvalidMaxGqPhred { got: f64 },

    /// An `sfs_prior_tables` entry was `NaN` or `+∞`. The SFS genotype prior's
    /// log-priors are finite or `−∞` (a floored genotype); `NaN`/`+∞` signal a
    /// construction bug upstream.
    #[error("sfs_prior_tables contained a non-finite log-prior (NaN or +inf), got {got}")]
    InvalidSfsPriorTable { got: f64 },

    /// `nucleotide_diversity` (`θ̂`) was non-finite or negative. It comes from
    /// [`crate::var_calling::diversity::DiversityEstimate`], which guarantees a
    /// finite, non-negative value; a bad value signals a wiring error.
    #[error("nucleotide_diversity must be finite and >= 0, got {got}")]
    InvalidNucleotideDiversity { got: f64 },
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

    /// **Historical — no longer raised by the exact-math EM loop.**
    /// Records that hit `max_iterations` without satisfying
    /// `convergence_threshold` are now emitted with
    /// `EmDiagnostics::converged = false` and routed by the VCF
    /// writer to `FILTER=EMNoConv`. The variant is retained because
    /// the engine's `#[non_exhaustive]` `Result` shape is part of
    /// its public API and future engine variants (e.g. the deferred
    /// approximate-LUT mode) may still need a hard non-convergence
    /// surface; removing it would be a SemVer-breaking change for
    /// callers exhaustively matching on the enum.
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

    /// B4: the column-native log-likelihood kernel was handed more
    /// alleles than its bitmask-based per-genotype bitset can index.
    /// The Phase A.1 unifier caps `n_alleles` upstream at
    /// `MAX_BITMASK_ALLELES`, so surfacing this variant means the
    /// caller passed `cfg.max_alleles > MAX_BITMASK_ALLELES` — an
    /// upstream-invariant break, **not** a math-degeneracy
    /// condition. The locus and `n_alleles` are the only useful
    /// diagnostic fields (no single sample / genotype to blame), so
    /// the variant carries them directly rather than synthesising
    /// `usize::MAX` placeholders on a sibling variant.
    #[error(
        "column-native log-likelihood kernel was handed n_alleles={n_alleles} \
         at {locus}, exceeding the bitmask cap — `max_alleles_lh_calc` \
         (or `max_alleles`) misconfigured"
    )]
    NAllelesExceedsBitmask {
        locus: RecordLocus,
        n_alleles: usize,
    },

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

    /// `PosteriorEngineConfig::sfs_prior_tables` length did not match the
    /// upstream record's `n_samples`. The precomputed per-sample genotype prior
    /// tables were built for a different cohort than the one driving Stage 5.
    #[error(
        "sfs_prior_tables cohort size {table_samples} does not match \
         record sample count {record_samples} at {locus}"
    )]
    SfsPriorTableCohortSizeMismatch {
        locus: RecordLocus,
        table_samples: usize,
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
    /// Hidden-paralog posterior probability `P(paralog | data)`, stamped by the
    /// paralog write pass when the locus was scored (biallelic SNP with a finite
    /// LR). `None` for loci the paralog filter did not score (indels,
    /// multiallelic, or the whole filter disabled) → the `PARALOG_POST` INFO
    /// field is then omitted. Not part of the EM output; set downstream.
    pub paralog_posterior: Option<f64>,
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

    /// Drop ALT alleles that no sample's argmax genotype references and
    /// re-index every allele-keyed field, so the emitted VCF lists only
    /// alleles the cohort actually carries. Returns the number of ALT
    /// alleles removed.
    ///
    /// Why this exists: the upstream allele set is the read-supported
    /// *candidate* set for the whole variant group (the cohort union,
    /// capped at `max_alleles`), and the engine genotypes over all of
    /// it. The per-sample argmax can leave candidates with zero genotype
    /// support; emitting them yields VCF ALT alleles with `AC=0` (and,
    /// after `bcftools norm -m`, phantom `0/0` records). GATK's
    /// `GenotypeGVCFs` prunes the same way. REF (index 0) is always kept.
    ///
    /// Only multiallelic records (`alleles.len() >= 3`) can carry an
    /// unsupported ALT *alongside* a supported one. A biallelic record's
    /// single ALT is either supported, or the whole record is hom-ref
    /// everywhere and dropped upstream by [`Self::is_variant_call`] — so
    /// those short-circuit with no work or allocation.
    ///
    /// Re-indexing: `alleles`, `allele_frequencies` (renormalised to sum
    /// 1), `compound_frequencies`, and the per-sample `scalars` /
    /// `chain_anchor_flags` allele columns are compacted to the kept
    /// alleles; `posteriors` are restricted to the surviving genotypes
    /// and renormalised per sample; `best_genotype` and `gq_phred` are
    /// recomputed from the restricted rows (GQ capped at `max_gq_phred`,
    /// matching the engine's `summarise_posteriors`); `n_genotypes` is
    /// updated. A dropped allele's per-sample read support leaves
    /// `scalars`, so it no longer counts toward DP/AD — the same
    /// convention the `max_alleles` cap path uses for its OTHER pool.
    pub fn prune_unsupported_alleles(&mut self, max_gq_phred: f64) -> usize {
        let n_alleles = self.alleles.len();
        if n_alleles < 3 {
            return 0;
        }
        let ploidy = self.ploidy;
        let old_table = genotype_order(ploidy, n_alleles);

        // 1. Mark alleles referenced by some sample's argmax genotype.
        let mut supported = vec![false; n_alleles];
        supported[0] = true; // REF is always kept
        for &gt_idx in &self.best_genotype {
            for &allele in &old_table[gt_idx] {
                supported[allele as usize] = true;
            }
        }
        if supported.iter().all(|&keep| keep) {
            return 0;
        }
        // If no ALT survives, the site is hom-ref in every sample. Leave
        // it intact: the whole-site `is_variant_call` filter drops it,
        // and pruning to a REF-only record here would needlessly destroy
        // the EM's allele-frequency estimates for a discarded record.
        if !supported[1..].iter().any(|&keep| keep) {
            return 0;
        }

        // 2. kept[new] = old allele index (REF first); compaction map.
        let kept: Vec<usize> = (0..n_alleles).filter(|&i| supported[i]).collect();
        let new_n_alleles = kept.len();
        let removed = n_alleles - new_n_alleles;

        // 3. New genotype enumeration + an old-tuple -> old-index lookup.
        //    Each new genotype's allele indices map back to old via
        //    `kept` (monotone, so the tuple stays non-decreasing) and we
        //    look its old position up to copy the right posterior.
        let new_table = genotype_order(ploidy, new_n_alleles);
        let new_n_genotypes = new_table.len();
        let old_index_of: std::collections::HashMap<&[u8], usize> = old_table
            .iter()
            .enumerate()
            .map(|(i, gt)| (gt.as_slice(), i))
            .collect();
        let mut new_to_old_gt = vec![0usize; new_n_genotypes];
        let mut tuple_buf: Vec<u8> = Vec::with_capacity(ploidy as usize);
        for (j, gt_new) in new_table.iter().enumerate() {
            tuple_buf.clear();
            tuple_buf.extend(gt_new.iter().map(|&a_new| kept[a_new as usize] as u8));
            new_to_old_gt[j] = old_index_of[tuple_buf.as_slice()];
        }

        // 4. Per-allele vectors (compact to kept; renormalise AF).
        let mut alleles = Vec::with_capacity(new_n_alleles);
        let mut allele_frequencies = Vec::with_capacity(new_n_alleles);
        let mut compound_frequencies = Vec::with_capacity(new_n_alleles);
        for &old_i in &kept {
            alleles.push(self.alleles[old_i].clone());
            allele_frequencies.push(self.allele_frequencies[old_i]);
            compound_frequencies.push(self.compound_frequencies[old_i]);
        }
        let af_sum: f64 = allele_frequencies.iter().sum();
        if af_sum > 0.0 {
            for v in &mut allele_frequencies {
                *v /= af_sum;
            }
        }

        // 5. Per-sample per-allele tables (drop the removed columns).
        let mut scalars = Vec::with_capacity(self.n_samples * new_n_alleles);
        let mut chain_anchor_flags = Vec::with_capacity(self.n_samples * new_n_alleles);
        for s in 0..self.n_samples {
            let base = s * n_alleles;
            for &old_i in &kept {
                scalars.push(self.scalars[base + old_i]);
                chain_anchor_flags.push(self.chain_anchor_flags[base + old_i]);
            }
        }

        // 6. Posteriors restricted to surviving genotypes + renormalised;
        //    best_genotype + GQ recomputed from the restricted rows.
        let mut posteriors = vec![0.0_f64; self.n_samples * new_n_genotypes];
        let mut best_genotype = vec![0usize; self.n_samples];
        let mut gq_phred = vec![0.0_f64; self.n_samples];
        for s in 0..self.n_samples {
            let old_base = s * self.n_genotypes;
            let new_base = s * new_n_genotypes;
            let mut row_sum = 0.0;
            for j in 0..new_n_genotypes {
                let p = self.posteriors[old_base + new_to_old_gt[j]];
                posteriors[new_base + j] = p;
                row_sum += p;
            }
            let inv = if row_sum > 0.0 { 1.0 / row_sum } else { 1.0 };
            let mut best_j = 0usize;
            let mut best_p = f64::NEG_INFINITY;
            for j in 0..new_n_genotypes {
                let p = posteriors[new_base + j] * inv;
                posteriors[new_base + j] = p;
                if p > best_p {
                    best_p = p;
                    best_j = j;
                }
            }
            best_genotype[s] = best_j;
            let p_best_clamped = best_p.min(1.0 - f64::EPSILON);
            gq_phred[s] = (PHRED_SCALE * (1.0 - p_best_clamped).log10()).clamp(0.0, max_gq_phred);
        }

        // 7. Commit. `n_alleles` is derived from `alleles.len()`, so
        //    swapping the vector updates it; `n_genotypes` is explicit.
        self.alleles = alleles;
        self.allele_frequencies = allele_frequencies;
        self.compound_frequencies = compound_frequencies;
        self.scalars = scalars;
        self.chain_anchor_flags = chain_anchor_flags;
        self.posteriors = posteriors;
        self.best_genotype = best_genotype;
        self.gq_phred = gq_phred;
        self.n_genotypes = new_n_genotypes;

        removed
    }

    /// True iff at least one sample's argmax genotype carries a
    /// non-reference allele.
    ///
    /// The canonical [`genotype_order`] table places the all-REF
    /// genotype at index 0 for every
    /// `(ploidy, n_alleles)` (e.g. `0/0` for diploid biallelic,
    /// `AAAA` for tetraploid biallelic), so the check is just
    /// `best_genotype.iter().any(|&i| i != 0)`.
    ///
    /// Used by the cohort driver to drop records that the per-sample
    /// genotype caller resolved to hom-ref everywhere — the EM
    /// allele-frequency estimate `p̂` can still be non-zero at such
    /// sites (it's a likelihood-weighted estimate from raw read
    /// counts, not from argmax genotypes), and emitting them would
    /// produce VCF rows with `AC=0` across all samples.
    pub fn is_variant_call(&self) -> bool {
        self.best_genotype.iter().any(|&i| i != 0)
    }
}

// VcfWritable: connect PosteriorRecord to the VCF writer's input
// contract. Every method here is a one-line accessor for an existing
// field or a delegation to the per-sample-row helper, so monomorphised
// callers inline back to a direct field access.
impl crate::vcf::VcfWritable for PosteriorRecord {
    fn chrom_id(&self) -> u32 {
        self.locus.chrom_id
    }
    fn pos_1based(&self) -> u32 {
        self.locus.start
    }
    fn ploidy(&self) -> u8 {
        self.ploidy
    }
    fn n_samples(&self) -> usize {
        self.n_samples
    }
    fn n_genotypes(&self) -> usize {
        self.n_genotypes
    }
    fn n_alleles(&self) -> usize {
        self.alleles.len()
    }
    fn allele_seq(&self, allele_idx: usize) -> &[u8] {
        &self.alleles[allele_idx].seq
    }
    fn qual_phred(&self) -> f64 {
        self.qual_phred
    }
    fn converged(&self) -> bool {
        self.diagnostics.converged
    }
    fn allele_frequencies(&self) -> &[f64] {
        &self.allele_frequencies
    }
    fn compound_frequencies(&self) -> &[Option<f64>] {
        &self.compound_frequencies
    }
    fn best_genotype(&self) -> &[usize] {
        &self.best_genotype
    }
    fn gq_phred(&self) -> &[f64] {
        &self.gq_phred
    }
    fn posteriors_row(&self, sample_idx: usize) -> &[f64] {
        self.posteriors_row(sample_idx)
    }
    fn scalars_row(&self, sample_idx: usize) -> &[AlleleSupportStats] {
        self.scalars_row(sample_idx)
    }
    fn chain_anchor_flags_row(&self, sample_idx: usize) -> &[bool] {
        self.chain_anchor_flags_row(sample_idx)
    }
    fn posteriors_len(&self) -> usize {
        self.posteriors.len()
    }
    fn scalars_len(&self) -> usize {
        self.scalars.len()
    }
    fn chain_anchor_flags_len(&self) -> usize {
        self.chain_anchor_flags.len()
    }
    fn paralog_posterior(&self) -> Option<f64> {
        self.paralog_posterior
    }
}

/// EM bookkeeping for a single record. A `PosteriorRecord` is emitted
/// regardless of EM convergence — `converged` carries that bit so the
/// VCF writer can route non-converging records into a flagged
/// `FILTER` value (`EMNoConv`) rather than hard-erroring the whole
/// cohort run on a single problematic site. Hard-erroring on
/// non-convergence was the v0 behaviour; the emit-with-flag policy
/// matches GATK / bcftools convention and preserves
/// downstream-filterable information that would otherwise be lost.
///
/// On `converged == true`: `final_max_delta_p < config.convergence_threshold`.
/// On `converged == false`: the EM loop hit `config.max_iterations`
/// without satisfying the threshold; the emitted record's
/// `allele_frequencies` / `posteriors` are computed from the
/// un-converged p̂/f̂ values via the same final E-step the
/// converged path uses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EmDiagnostics {
    /// Iterations actually run (≤ `config.max_iterations`).
    pub iterations: u32,
    /// `max_a |p̂_new[a] − p̂_old[a]|` at the last iteration.
    pub final_max_delta_p: f64,
    /// `true` iff the EM loop exited via the convergence threshold
    /// rather than the iteration cap.
    pub converged: bool,
}

// ---------------------------------------------------------------------
// Phase A.2 step 1 — column-native EM input boundary
//
// `run_em_for_record(record: MergedRecord)` keeps its existing
// signature for the row-shape consumers (the in-module test surface +
// `var-calling-from-bam`'s streaming pipeline). It is a thin shim
// around the borrowed-slice [`run_em_columnar`] body, which reads its
// per-allele queries through the [`AllelesView`] trait. Two
// `AllelesView` impls live in this crate:
//
// - [`MergedAllelesView`] here — wraps `&[MergedAllele]`. Used by the
//   row-shape shim.
// - `ColumnarAllelesView` in
//   [`crate::var_calling::worker`] — wraps
//   `&UnifiedAllelesColumns`. Used by the column-native worker once
//   Phase A.2 step 2 lands.
//
// The hot EM loops never touch the view: they read
// `scratch.compound_mask: &[bool]` which is populated once per record
// via `inputs.alleles.is_compound(a)`. The view's `&dyn` virtual
// dispatch only costs at the per-record per-allele classify step.

/// Minimal accessor surface for the per-allele queries the EM body
/// makes on the allele set. Designed for cheap implementation over
/// both row-shape `&[MergedAllele]` and column-native
/// `&UnifiedAllelesColumns`.
///
/// The EM reads:
/// - [`Self::len`] — for `n_alleles` bound checks.
/// - [`Self::ref_len`] — for `classify_allele`'s SNP-vs-indel test.
/// - [`Self::seq_len`] — per-allele, also for `classify_allele`.
/// - [`Self::is_compound`] — for `compound_mask`,
///   `mixture_contam_class_per_allele`, and `compound_frequencies`.
///
/// The EM never asks for per-allele bytes or compound constituents —
/// those live on `PosteriorRecord` for the VCF writer's `ALT` /
/// compound-emit path, downstream of the EM.
pub(crate) trait AllelesView {
    fn len(&self) -> usize;
    fn ref_len(&self) -> usize;
    fn seq_len(&self, allele_idx: usize) -> usize;
    fn is_compound(&self, allele_idx: usize) -> bool;
}

/// Row-shape [`AllelesView`] over a `&[MergedAllele]` slice. Used by
/// the [`run_em_for_record`] row-shape shim.
pub(crate) struct MergedAllelesView<'a> {
    alleles: &'a [MergedAllele],
}

impl<'a> MergedAllelesView<'a> {
    pub(crate) fn new(alleles: &'a [MergedAllele]) -> Self {
        Self { alleles }
    }
}

impl AllelesView for MergedAllelesView<'_> {
    fn len(&self) -> usize {
        self.alleles.len()
    }
    fn ref_len(&self) -> usize {
        self.alleles.first().map_or(0, |a| a.seq.len())
    }
    fn seq_len(&self, allele_idx: usize) -> usize {
        self.alleles[allele_idx].seq.len()
    }
    fn is_compound(&self, allele_idx: usize) -> bool {
        self.alleles[allele_idx].is_compound
    }
}

/// Borrowed-slice view of one group's EM inputs. Same logical content
/// as a [`MergedRecord`] minus the per-allele bytes (carried through
/// `alleles` as an [`AllelesView`]) and minus the
/// `other_scalars` / `chain_anchor_flags` row-shape buffers that the
/// EM body never reads (those are forwarded straight to
/// [`PosteriorRecord`] by the caller).
///
/// The EM body itself reads:
/// - `locus`, `ploidy`, `n_samples`, `n_genotypes` — scalars.
/// - `alleles` via the [`AllelesView`] trait — per-allele classify.
/// - `scalars` — only the mixture pre-pass reads this (contamination on).
/// - `log_likelihoods` — the E-step's core input.
/// - `chain_anchor_flags_for_validation` — only [`validate_record_shape`]
///   reads this, as a length check. The bool table itself never feeds
///   any math; it is forwarded to `PosteriorRecord` by the caller.
pub(crate) struct EmInputs<'a> {
    pub locus: RecordLocus,
    pub ploidy: u8,
    pub n_samples: usize,
    pub n_genotypes: usize,
    pub alleles: &'a dyn AllelesView,
    pub scalars: &'a [AlleleSupportStats],
    pub log_likelihoods: &'a [f64],
    pub chain_anchor_flags_for_validation: &'a [bool],
}

/// Owned outputs produced by [`run_em_columnar`]. The caller joins
/// these with the per-record allele bytes + row-shape
/// `scalars` / `other_scalars` / `chain_anchor_flags` to form the
/// emitted [`PosteriorRecord`].
pub(crate) struct EmOutputs {
    pub n_genotypes: usize,
    pub allele_frequencies: Vec<f64>,
    pub compound_frequencies: Vec<Option<f64>>,
    pub posteriors: Vec<f64>,
    pub best_genotype: Vec<usize>,
    pub gq_phred: Vec<f64>,
    pub qual_phred: f64,
    pub diagnostics: EmDiagnostics,
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
/// use pop_var_caller::var_calling::per_group_merger::{
///     MergedRecord, PerGroupMergerError,
/// };
/// use pop_var_caller::var_calling::posterior_engine::{
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
pub struct PosteriorEngine<I, M = InterpUnivariateSimdMath>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
    M: MathBackend,
{
    upstream: I,
    config: PosteriorEngineConfig,
    math: M,
    is_latched: bool,
    /// Reusable per-record scratch buffers. Allocated empty here and
    /// resized once per record in `run_em_for_record`; steady-state
    /// records of the same shape allocate nothing. See
    /// [`RecordScratch`] for the layout and the post-H4 profile
    /// motivation.
    scratch: RecordScratch,
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
            scratch: _,
        } = self;
        f.debug_struct("PosteriorEngine")
            .field("config", config)
            .field("is_latched", is_latched)
            .finish()
    }
}

impl<I> PosteriorEngine<I, InterpUnivariateSimdMath>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
{
    /// Construct an engine with the [project defaults] and the default
    /// math backend ([`InterpUnivariateSimdMath`], ~25 % faster than
    /// [`ExactMath`](backends::ExactMath) on the contam-on bench with ~50× margin under the
    /// approximate-parity budget — see the 2026-05-18 SIMD analysis
    /// report).
    ///
    /// For bit-identical reproducibility against the unoptimised
    /// engine, construct with
    /// [`with_math_backend`](Self::with_math_backend) and pass
    /// [`ExactMath`](backends::ExactMath) explicitly.
    ///
    /// [project defaults]: PosteriorEngineConfig::with_project_defaults
    pub fn new(upstream: I) -> Self {
        Self::with_config(upstream, PosteriorEngineConfig::with_project_defaults())
    }

    /// Construct an engine with explicit tuning and the default math
    /// backend ([`InterpUnivariateSimdMath`]).
    pub fn with_config(upstream: I, config: PosteriorEngineConfig) -> Self {
        Self::with_math_backend(upstream, config, InterpUnivariateSimdMath)
    }
}

impl<I, M> PosteriorEngine<I, M>
where
    I: Iterator<Item = Result<MergedRecord, PerGroupMergerError>>,
    M: MathBackend,
{
    /// Construct an engine with explicit tuning and an explicit math
    /// backend. Use this to opt into a scalar-only approximate backend
    /// (`InterpUnivariateMath`) or to opt back into [`ExactMath`](backends::ExactMath) if
    /// bit-identical reproducibility against the unoptimised engine
    /// matters.
    pub fn with_math_backend(upstream: I, config: PosteriorEngineConfig, math: M) -> Self {
        Self {
            upstream,
            config,
            math,
            is_latched: false,
            scratch: RecordScratch::empty(),
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
        run_em_for_record(record, &self.config, &self.math, &mut self.scratch)
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

/// Single-allele classifier. Shared between [`run_em_for_record`]
/// (which writes the result into `scratch.allele_classes` per allele
/// without a Vec allocation) and the test-only [`classify_alleles`]
/// helper.
#[cfg(test)]
fn classify_allele(allele: &MergedAllele, ref_len: usize, is_first: bool) -> AlleleClass {
    if is_first {
        AlleleClass::Ref
    } else if allele.is_compound {
        AlleleClass::CompoundAlt
    } else if allele.seq.len() == ref_len {
        AlleleClass::SnpAlt
    } else {
        AlleleClass::IndelAlt
    }
}

#[cfg(test)]
fn classify_alleles(alleles: &[MergedAllele]) -> Vec<AlleleClass> {
    let ref_len = alleles[0].seq.len();
    alleles
        .iter()
        .enumerate()
        .map(|(idx, allele)| classify_allele(allele, ref_len, idx == 0))
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
/// Fill `scratch.mixture_log_likelihoods` with the
/// mixture-recomputed per-(sample, genotype) log-likelihood under
/// the configured contamination estimates. Scalar fallback path;
/// `compute_mixture_log_likelihoods_simd` is the AVX2 default.
///
/// Caller invariants: `scratch.mixture_contam_class_per_allele` has
/// been filled, and `scratch.mixture_log_likelihoods` /
/// `scratch.mixture_n_obs` / `scratch.mixture_mean_err` /
/// `scratch.mixture_c_s_all` are resized to the current record's
/// shape.
#[allow(clippy::too_many_arguments)]
fn compute_mixture_log_likelihoods<M: MathBackend>(
    locus: RecordLocus,
    shape: &GenotypeShape,
    n_alleles: usize,
    n_genotypes: usize,
    n_samples: usize,
    ploidy: u8,
    scalars: &[AlleleSupportStats],
    estimates: &ContaminationEstimates,
    fallback_log_likelihoods: &[f64],
    math: &M,
    scratch: &mut RecordScratch,
) -> Result<(), PosteriorEngineError> {
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

    let genotype_allele_counts = &shape.genotype_allele_counts;

    let ploidy_f = f64::from(ploidy);
    // Denominator for spreading per-base error mass across alleles
    // other than the one the genotype carries; floored at 1 so the
    // n_alleles == 1 degenerate case doesn't divide by zero.
    let other_allele_error_denom = (n_alleles as f64 - 1.0).max(1.0);

    for sample_idx in 0..n_samples {
        let c_s = estimates.effective_c_s(sample_idx);
        if c_s <= 0.0 {
            // Sample is below-floor / singleton-batch: c_s = 0, the
            // mixture collapses to the own-DNA term. Fall back to
            // Stage 5's precomputed value to avoid re-deriving the
            // same number with floating-point drift.
            let dst = &mut scratch.mixture_log_likelihoods
                [sample_idx * n_genotypes..(sample_idx + 1) * n_genotypes];
            let src =
                &fallback_log_likelihoods[sample_idx * n_genotypes..(sample_idx + 1) * n_genotypes];
            dst.copy_from_slice(src);
            continue;
        }
        let q_b = estimates.q_b_for_sample(sample_idx);

        // Per-allele observation count and mean per-read error rate
        // for this sample at this site, recovered from the row of
        // AlleleSupportStats Stage 5 emitted. Both scratch slots are
        // fully rewritten below, so no per-iter clear needed.
        let sample_allele_scalars = &scalars[sample_idx * n_alleles..(sample_idx + 1) * n_alleles];
        for (a, support) in sample_allele_scalars.iter().enumerate() {
            scratch.mixture_n_obs[a] = support.num_obs;
            scratch.mixture_mean_err[a] = if support.num_obs > 0 {
                let per_read_log_err = support.q_sum / f64::from(support.num_obs);
                math.exp(per_read_log_err)
                    .clamp(MIN_BASE_ERROR, MAX_BASE_ERROR)
            } else {
                0.0
            };
        }

        for (g_idx, gt_counts) in genotype_allele_counts.chunks_exact(n_alleles).enumerate() {
            let mut ll: f64 = 0.0;
            for (a, &k) in gt_counts.iter().enumerate() {
                let n_a = scratch.mixture_n_obs[a];
                if n_a == 0 {
                    continue;
                }
                let k_a = f64::from(k);
                let eps = scratch.mixture_mean_err[a];
                // P(read carrying allele a | genotype G).
                let p_own = (k_a / ploidy_f) * (1.0 - eps)
                    + ((ploidy_f - k_a) / ploidy_f) * (eps / other_allele_error_denom);
                let p_contam = match scratch.mixture_contam_class_per_allele[a] {
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
            scratch.mixture_log_likelihoods[sample_idx * n_genotypes + g_idx] = ll;
        }
    }

    Ok(())
}

/// Lane-of-4 SIMD variant of [`compute_mixture_log_likelihoods`].
/// Process samples in batches of 4 using `f64x4` math; the per-sample
/// outer loop, the per-(g, a) inner work, and the dominant
/// `math.ln(mix)` call all go SIMD via [`MathBackend::ln_x4`].
///
/// The samply profile from 2026-05-18 showed ~12 % of contam-on
/// self-time in the scalar mixture pre-pass (3.8 % on `math.ln(mix)`
/// at line 824 of the scalar version, plus ~8 % in `interp::*_approx`
/// driven from that call). The SIMD body replaces it.
///
/// Dispatched from [`run_em_for_record`] when `M::HAS_LANE_4` is
/// true. Scalar backends fall through to
/// [`compute_mixture_log_likelihoods`].
#[allow(clippy::too_many_arguments)]
fn compute_mixture_log_likelihoods_simd<M: MathBackend>(
    locus: RecordLocus,
    shape: &GenotypeShape,
    n_alleles: usize,
    n_genotypes: usize,
    n_samples: usize,
    ploidy: u8,
    scalars: &[AlleleSupportStats],
    estimates: &ContaminationEstimates,
    fallback_log_likelihoods: &[f64],
    math: &M,
    scratch: &mut RecordScratch,
) -> Result<(), PosteriorEngineError> {
    use wide::{CmpGt, CmpLe, f64x4};

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

    let genotype_allele_counts = &shape.genotype_allele_counts;

    let ploidy_f = f64::from(ploidy);
    let other_allele_error_denom = (n_alleles as f64 - 1.0).max(1.0);

    // Per-sample pre-pass: collect c_s and mean_err for every
    // (sample, allele) cell so the SIMD batch loop can lane-load
    // them with simple gathers. n_obs reads directly from `scalars`
    // (its `num_obs: u32` field is already where we'd put it).
    for sample_idx in 0..n_samples {
        scratch.mixture_c_s_all[sample_idx] = estimates.effective_c_s(sample_idx);
        let sample_allele_scalars = &scalars[sample_idx * n_alleles..(sample_idx + 1) * n_alleles];
        for (a, support) in sample_allele_scalars.iter().enumerate() {
            scratch.mixture_mean_err_all[sample_idx * n_alleles + a] = if support.num_obs > 0 {
                let per_read_log_err = support.q_sum / f64::from(support.num_obs);
                math.exp(per_read_log_err)
                    .clamp(MIN_BASE_ERROR, MAX_BASE_ERROR)
            } else {
                0.0
            };
        }
    }

    // Initialize output with the fallback (Stage-5 precomputed)
    // values. Samples with `c_s <= 0` (inactive in the mixture model)
    // keep their fallback row through the SIMD loop; active lanes
    // overwrite their per-genotype cell on the way out.
    scratch
        .mixture_log_likelihoods
        .copy_from_slice(fallback_log_likelihoods);

    let n_full_batches = n_samples / 4;
    let tail_start = n_full_batches * 4;

    for batch_idx in 0..n_full_batches {
        let s0 = batch_idx * 4;

        let c_s_v = f64x4::from([
            scratch.mixture_c_s_all[s0],
            scratch.mixture_c_s_all[s0 + 1],
            scratch.mixture_c_s_all[s0 + 2],
            scratch.mixture_c_s_all[s0 + 3],
        ]);
        let one_minus_c_s_v = f64x4::splat(1.0) - c_s_v;
        let active_mask = c_s_v.simd_gt(f64x4::splat(0.0));

        if !active_mask.any() {
            // All 4 lanes inactive — fallback rows already in place,
            // nothing more to do for this batch.
            continue;
        }

        for g_idx in 0..n_genotypes {
            let gt_offset = g_idx * n_alleles;
            let mut ll_v = f64x4::splat(0.0);

            for a in 0..n_alleles {
                // `n_a` per lane (different samples at fixed allele).
                let n_a_v = f64x4::from([
                    f64::from(scalars[s0 * n_alleles + a].num_obs),
                    f64::from(scalars[(s0 + 1) * n_alleles + a].num_obs),
                    f64::from(scalars[(s0 + 2) * n_alleles + a].num_obs),
                    f64::from(scalars[(s0 + 3) * n_alleles + a].num_obs),
                ]);
                let n_a_mask = n_a_v.simd_gt(f64x4::splat(0.0));

                // `eps` per lane.
                let eps_v = f64x4::from([
                    scratch.mixture_mean_err_all[s0 * n_alleles + a],
                    scratch.mixture_mean_err_all[(s0 + 1) * n_alleles + a],
                    scratch.mixture_mean_err_all[(s0 + 2) * n_alleles + a],
                    scratch.mixture_mean_err_all[(s0 + 3) * n_alleles + a],
                ]);

                // p_own depends on (eps, k_a, ploidy_f). The
                // coefficients are lane-uniform (depend on g, a, not
                // sample), so we splat them and combine with the
                // per-lane `eps`.
                let k_a = f64::from(genotype_allele_counts[gt_offset + a]);
                let coeff_own = k_a / ploidy_f;
                let coeff_err = (ploidy_f - k_a) / (ploidy_f * other_allele_error_denom);
                let p_own_v = f64x4::splat(coeff_own) * (f64x4::splat(1.0) - eps_v)
                    + f64x4::splat(coeff_err) * eps_v;

                // p_contam per lane: gather the per-sample-batch
                // q_b value for the current allele's contamination
                // class. Compound-allele class → `None` → 0.0 across
                // all lanes.
                let p_contam_v = match scratch.mixture_contam_class_per_allele[a] {
                    Some(class) => {
                        let ci = class as usize;
                        f64x4::from([
                            estimates.q_b_for_sample(s0)[ci],
                            estimates.q_b_for_sample(s0 + 1)[ci],
                            estimates.q_b_for_sample(s0 + 2)[ci],
                            estimates.q_b_for_sample(s0 + 3)[ci],
                        ])
                    }
                    None => f64x4::splat(0.0),
                };

                let mix_v = one_minus_c_s_v * p_own_v + c_s_v * p_contam_v;

                // Defensive: every active lane with `n_a > 0` must
                // produce `mix > 0`. If any such lane has `mix <= 0`,
                // surface the first offending lane as a
                // non-finite-posterior error (matches the scalar
                // body's behaviour exactly).
                let contributing_mask = active_mask & n_a_mask;
                let bad_mask = mix_v.simd_le(f64x4::splat(0.0)) & contributing_mask;
                if bad_mask.any() {
                    let mix_arr = mix_v.to_array();
                    for lane in 0..4 {
                        if scratch.mixture_c_s_all[s0 + lane] > 0.0
                            && scalars[(s0 + lane) * n_alleles + a].num_obs > 0
                            && mix_arr[lane] <= 0.0
                        {
                            return Err(PosteriorEngineError::NonFinitePosterior {
                                locus,
                                sample_idx: s0 + lane,
                                genotype_idx: Some(g_idx),
                                kind: NonFiniteKind::NegativeInfinity,
                            });
                        }
                    }
                }

                // Replace `mix` with 1.0 in lanes that won't
                // contribute, so `ln` doesn't see 0 / negative inputs
                // (which would produce `-∞` / `NaN` and poison the
                // sum via `0 * -∞ = NaN`).
                let mix_safe = contributing_mask.blend(mix_v, f64x4::splat(1.0));
                let log_mix_v = math.ln_x4(mix_safe);
                let contrib_v = n_a_v * log_mix_v;
                let contrib_masked = contributing_mask.blend(contrib_v, f64x4::splat(0.0));
                ll_v += contrib_masked;
            }

            // Scatter back: each active lane writes its cell;
            // inactive lanes keep the fallback already copied into
            // `scratch.mixture_log_likelihoods` at function entry.
            let ll_arr = ll_v.to_array();
            for (lane, &ll_lane) in ll_arr.iter().enumerate() {
                if scratch.mixture_c_s_all[s0 + lane] > 0.0 {
                    scratch.mixture_log_likelihoods[(s0 + lane) * n_genotypes + g_idx] = ll_lane;
                }
            }
        }
    }

    // Tail: process leftover samples (n_samples % 4) with scalar
    // math. The pre-pass already populated `scratch.mixture_mean_err_all`,
    // so the tail just reads from it.
    for sample_idx in tail_start..n_samples {
        let c_s = scratch.mixture_c_s_all[sample_idx];
        if c_s <= 0.0 {
            continue; // fallback row already in place
        }
        let q_b = estimates.q_b_for_sample(sample_idx);
        let sample_allele_scalars = &scalars[sample_idx * n_alleles..(sample_idx + 1) * n_alleles];

        for (g_idx, gt_counts) in genotype_allele_counts.chunks_exact(n_alleles).enumerate() {
            let mut ll: f64 = 0.0;
            for a in 0..n_alleles {
                let n_a = sample_allele_scalars[a].num_obs;
                if n_a == 0 {
                    continue;
                }
                let k_a = f64::from(gt_counts[a]);
                let eps = scratch.mixture_mean_err_all[sample_idx * n_alleles + a];
                let p_own = (k_a / ploidy_f) * (1.0 - eps)
                    + ((ploidy_f - k_a) / ploidy_f) * (eps / other_allele_error_denom);
                let p_contam = match scratch.mixture_contam_class_per_allele[a] {
                    Some(class) => q_b[class as usize],
                    None => 0.0,
                };
                let mix = (1.0 - c_s) * p_own + c_s * p_contam;
                if mix <= 0.0 {
                    return Err(PosteriorEngineError::NonFinitePosterior {
                        locus,
                        sample_idx,
                        genotype_idx: Some(g_idx),
                        kind: NonFiniteKind::NegativeInfinity,
                    });
                }
                ll += f64::from(n_a) * math.ln(mix);
            }
            scratch.mixture_log_likelihoods[sample_idx * n_genotypes + g_idx] = ll;
        }
    }

    Ok(())
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
/// Per-record EM context: scalar parameters + a shared reference to
/// the cached genotype shape.
///
/// **No slice fields.** The arrays that used to live here as borrowed
/// slices (`compound_mask`, `pseudocounts`, `log_f_per_sample`,
/// `log_one_minus_f_per_sample`) now live in [`RecordScratch`]; EM
/// reads them through the same `&mut RecordScratch` it already
/// receives, sidestepping the borrow-checker conflict between holding
/// `&EmContext<'_>` (which would borrow scratch) and `&mut
/// RecordScratch` simultaneously. `EmContext` is `Copy` so passing it
/// by value to each EM helper is free.
#[derive(Clone, Copy)]
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
    /// `(ploidy, n_alleles)` via the thread-local cache. The borrow
    /// is held for the duration of `run_em_for_record` against the
    /// `Arc<GenotypeShape>` returned by `shape_for(...)`.
    shape: &'a GenotypeShape,
    /// Every sample shares the same `f_s` value (e.g. the default
    /// config sets `fixation_index_default = 0.0` and leaves
    /// `fixation_index_overrides = None`). When `true`,
    /// `scratch.log_f_per_sample[0]` and
    /// `scratch.log_one_minus_f_per_sample[0]` are the cohort-wide
    /// values, and the entire per-genotype log-prior (including the
    /// homozygous-IBD `log_sum_exp_2` term) is sample-invariant
    /// within an EM iteration — `e_step` / `e_step_simd` precompute
    /// it once into `RecordScratch::log_prior_per_g` and splat from
    /// there. See H2 + H4 in
    /// `doc/devel/reports/reviews/perf_posterior_engine_2026-05-18.md`.
    homogeneous_fixation: bool,
    compound_pseudocount: f64,
    /// Per-sample biallelic-diploid genotype log-priors from the SFS-marginal
    /// prior, borrowed from `config.sfs_prior_tables`. `Some` **only** when the
    /// prior is configured *and* this record is biallelic-diploid
    /// (`ploidy == 2, n_alleles == 2`); the E-step then reads the fixed
    /// per-sample table instead of building the HWE(p̂) prior. `None` routes to
    /// the HWE(p̂) path (multiallelic/non-diploid records, or the prior off).
    sfs_prior_tables: Option<&'a [[f64; 3]]>,
    /// `true` when `config.nucleotide_diversity` is set: the genotype prior is
    /// the general Dirichlet-multinomial (concentration `α` from θ̂, precomputed
    /// per record into `scratch.alpha` / `scratch.lgamma_alpha` /
    /// `scratch.log_alpha_over_sum`) for **every** shape. The E-step then fills
    /// `log_p_effective` from `log(α_a/Σα)` (the IBD marginal) and
    /// `fill_log_indep_per_g` uses the `lgamma` form instead of the plug-in
    /// `Σ k_a·log p̂_a`. `false` keeps the HWE(p̂) plug-in prior.
    use_dirichlet_prior: bool,
}

/// Per-engine scratch buffers reused across every record.
///
/// `PosteriorEngine` owns one instance; [`run_em_for_record`] resizes
/// each `Vec` to the current record's shape on entry, then the EM
/// loop overwrites cells in place. Steady-state allocation count
/// drops to ~zero for the intermediate buffers — `Vec::resize` keeps
/// capacity when the new shape is ≤ the high-water mark and only
/// reallocates when growing past it. The May-18 post-H4 sampling
/// profile showed glibc allocator self-time at ~16 % of cycles
/// before this lift; see
/// `doc/devel/reports/implementations/posterior_engine_post_h4_profile_2026-05-18.md`.
///
/// Output buffers (`p_hat`, `f_hat_compound`, `posteriors`,
/// `best_genotype`, `gq_phred`, `compound_frequencies`) are filled
/// here and `.clone()`-d out into [`PosteriorRecord`] on the way
/// back — trades one malloc-then-zero-init per record for one
/// malloc-then-memcpy at the same size. Memcpy beats zero-init on
/// the small Vecs and ties on the large `posteriors` buffer; on the
/// small Vecs (`p_hat`, `f_hat_compound`, `compound_frequencies`)
/// where malloc overhead dominates over the byte cost, lifting them
/// into the scratch removes the alloc entirely from the EM
/// iteration loop (`m_step_*_p_hat` no longer collects a fresh `Vec`
/// every EM iter).
pub(crate) struct RecordScratch {
    // ===== EM working set (previously EmScratch). =====
    /// Effective per-allele frequency seen by the prior (compounds
    /// substituted with `f̂_C`). Length `n_alleles`.
    p_effective: Vec<f64>,
    /// Natural log of `p_effective`. Length `n_alleles`.
    log_p_effective: Vec<f64>,
    /// Dirichlet concentration `α` for the general SFS prior (`α_ref = 1`,
    /// `α_alt = θ̂/(n_alleles−1)`); filled once per record from
    /// `config.nucleotide_diversity` when `ctx.use_dirichlet_prior`. Length
    /// `n_alleles`. Empty/unused under the plug-in HWE prior.
    alpha: Vec<f64>,
    /// `lgamma(α_a)` per allele — the baseline the Dirichlet-multinomial
    /// independent term subtracts. Length `n_alleles`. Filled beside `alpha`.
    lgamma_alpha: Vec<f64>,
    /// `log(α_a / Σα)` per allele — the IBD marginal allele frequency the
    /// Dirichlet-multinomial prior feeds to the Wright-`F` mixture in place of
    /// `log p̂_a`. The E-step copies it into `log_p_effective` under
    /// `ctx.use_dirichlet_prior`. Length `n_alleles`.
    log_alpha_over_sum: Vec<f64>,
    /// Per-genotype `log_indep` = `log_multinomial_coeffs[g] + Σ k *
    /// log_p_effective[a]` over non-zero pairs. Sample-invariant
    /// within one EM iteration; `e_step` / `e_step_simd` rebuild it
    /// from the current `log_p_effective` at the top of each call
    /// (H2 in the perf review). Length `n_genotypes`.
    log_indep_per_g: Vec<f64>,
    /// Per-genotype log-prior under the homogeneous-fixation hoist.
    /// Filled at the top of each `e_step` call when
    /// `ctx.homogeneous_fixation` (H4). Length `n_genotypes`.
    log_prior_per_g: Vec<f64>,
    /// Per-genotype unnormalised log-posterior for the current
    /// sample. Length `n_genotypes`.
    log_post_unnorm: Vec<f64>,
    /// SIMD-lane buffer: per-genotype `f64x4` of unnormalised
    /// log-posteriors across a 4-sample batch. Length `n_genotypes`.
    log_post_unnorm_lane: Vec<wide::f64x4>,
    /// Posterior-weighted allele counts E[n_a]. Length `n_alleles`.
    expected_counts: Vec<f64>,

    // ===== Record-static intermediates (no copy-out). =====
    /// Per-sample fixation indices `f_s`. Length `n_samples`.
    fixation_indices: Vec<f64>,
    /// `safe_ln(f_s)` per sample; hoisted out of `e_step`'s
    /// per-iteration sample loop. Length `n_samples`.
    log_f_per_sample: Vec<f64>,
    /// `safe_ln(1.0 - f_s)` per sample. Length `n_samples`.
    log_one_minus_f_per_sample: Vec<f64>,
    /// Per-allele [`AlleleClass`]. Length `n_alleles`.
    allele_classes: Vec<AlleleClass>,
    /// Per-allele Dirichlet pseudocount routed from
    /// [`PosteriorEngineConfig`] by [`AlleleClass`]. Length
    /// `n_alleles`.
    pseudocounts: Vec<f64>,
    /// Per-allele compound-allele flag. Length `n_alleles`.
    compound_mask: Vec<bool>,
    /// Per-allele contamination class — `None` for compound alleles
    /// (the mixture's `p_contam` term collapses to 0). Length
    /// `n_alleles`. Only populated when contamination is configured.
    mixture_contam_class_per_allele: Vec<Option<ContamAlleleClass>>,

    // ===== EM state (current + next for the convergence test). =====
    /// Current allele-frequency estimate. Length `n_alleles`.
    p_hat: Vec<f64>,
    /// Next-iteration allele-frequency estimate. After the M-step
    /// fills it, `max_abs_diff(p_hat, p_hat_next)` drives convergence
    /// and `mem::swap` rotates the roles. Length `n_alleles`.
    p_hat_next: Vec<f64>,
    /// Current compound-frequency estimate. Length `n_alleles`.
    f_hat_compound: Vec<f64>,
    /// Per-`(sample, genotype)` posteriors, sample-major. Length
    /// `n_samples * n_genotypes`.
    posteriors: Vec<f64>,

    // ===== Output buffers (cloned into PosteriorRecord at the end). =====
    /// Per-sample argmax genotype index. Length `n_samples`.
    best_genotype: Vec<usize>,
    /// Per-sample GQ in phred units. Length `n_samples`.
    gq_phred: Vec<f64>,
    /// Per-allele `Option<f64>` of the compound-frequency estimate
    /// (`None` for non-compound alleles). Length `n_alleles`.
    compound_frequencies: Vec<Option<f64>>,

    // ===== Site-QUAL exact-AF buffers (used after EM converges). =====
    /// Per-sample log-likelihoods collapsed by total non-ref allele
    /// count. Layout: row-major `[sample_idx * (ploidy + 1) +
    /// non_ref_count]`. Filled by [`compute_qual_via_exact_af`] before
    /// the convolution. Length `n_samples * (ploidy + 1)`.
    per_sample_alt_log_lik: Vec<f64>,
    /// Rolling per-cohort-AC log-probability buffer for the exact-AF
    /// convolution. After processing the first `s` samples,
    /// `log_p_ac_curr[k]` holds `log P(seen data | AC = k)` for
    /// `k = 0..=s·ploidy`. After the prior step it holds the
    /// unnormalised log-posterior over `K = 0..=ploidy·n_samples`.
    /// Length `ploidy * n_samples + 1`.
    log_p_ac_curr: Vec<f64>,
    /// Swap target for [`Self::log_p_ac_curr`] during the convolution.
    /// Same length.
    log_p_ac_next: Vec<f64>,
    /// Left-padded (`ploidy` leading `0.0` slots) **linear-domain**
    /// rolling AC buffer for the convolution ([`convolve_ac_linear`]).
    /// The hot O(n_samples²·ploidy) fold runs here as a transcendental-
    /// free multiply-add FIR (the log-domain fold spent ~88 % of the
    /// consume path in scalar `exp`/`ln`); per-sample renormalisation
    /// keeps the values in range and the running log-scale recovers the
    /// log-domain result afterwards. Length `ploidy + ploidy*n_samples + 1`.
    log_p_ac_curr_lin: Vec<f64>,
    /// Swap target for [`Self::log_p_ac_curr_lin`]. Same length.
    log_p_ac_next_lin: Vec<f64>,

    // ===== Mixture pre-pass buffers (only used when contam-on). =====
    /// Mixture-recomputed log-likelihoods. Length
    /// `n_samples * n_genotypes`. Filled by the mixture pre-pass when
    /// contamination is configured; the EM loop reads from here
    /// instead of the upstream `MergedRecord.log_likelihoods`.
    mixture_log_likelihoods: Vec<f64>,
    /// Per-allele observation count for the current sample in the
    /// scalar mixture pre-pass. Length `n_alleles`.
    mixture_n_obs: Vec<u32>,
    /// Per-allele mean per-read error rate for the current sample
    /// in the scalar mixture pre-pass. Length `n_alleles`.
    mixture_mean_err: Vec<f64>,
    /// Per-sample contamination fraction for the SIMD mixture
    /// pre-pass. Length `n_samples`.
    mixture_c_s_all: Vec<f64>,
    /// Per-`(sample, allele)` mean per-read error rate for the SIMD
    /// mixture pre-pass. Length `n_samples * n_alleles`.
    mixture_mean_err_all: Vec<f64>,
}

impl RecordScratch {
    /// All Vecs empty; capacity grows on first `resize_to` call.
    pub(crate) fn empty() -> Self {
        Self {
            p_effective: Vec::new(),
            log_p_effective: Vec::new(),
            alpha: Vec::new(),
            lgamma_alpha: Vec::new(),
            log_alpha_over_sum: Vec::new(),
            log_indep_per_g: Vec::new(),
            log_prior_per_g: Vec::new(),
            log_post_unnorm: Vec::new(),
            log_post_unnorm_lane: Vec::new(),
            expected_counts: Vec::new(),
            fixation_indices: Vec::new(),
            log_f_per_sample: Vec::new(),
            log_one_minus_f_per_sample: Vec::new(),
            allele_classes: Vec::new(),
            pseudocounts: Vec::new(),
            compound_mask: Vec::new(),
            mixture_contam_class_per_allele: Vec::new(),
            p_hat: Vec::new(),
            p_hat_next: Vec::new(),
            f_hat_compound: Vec::new(),
            posteriors: Vec::new(),
            best_genotype: Vec::new(),
            gq_phred: Vec::new(),
            compound_frequencies: Vec::new(),
            per_sample_alt_log_lik: Vec::new(),
            log_p_ac_curr: Vec::new(),
            log_p_ac_next: Vec::new(),
            log_p_ac_curr_lin: Vec::new(),
            log_p_ac_next_lin: Vec::new(),
            mixture_log_likelihoods: Vec::new(),
            mixture_n_obs: Vec::new(),
            mixture_mean_err: Vec::new(),
            mixture_c_s_all: Vec::new(),
            mixture_mean_err_all: Vec::new(),
        }
    }
}

impl Default for RecordScratch {
    fn default() -> Self {
        Self::empty()
    }
}

impl RecordScratch {
    /// Size every buffer to the current record's shape. `Vec::resize`
    /// preserves capacity above the high-water mark, so steady-state
    /// records of the same shape allocate nothing.
    fn resize_to(&mut self, n_samples: usize, n_alleles: usize, n_genotypes: usize, ploidy: u8) {
        self.p_effective.resize(n_alleles, 0.0);
        self.log_p_effective.resize(n_alleles, 0.0);
        self.alpha.resize(n_alleles, 0.0);
        self.lgamma_alpha.resize(n_alleles, 0.0);
        self.log_alpha_over_sum.resize(n_alleles, 0.0);
        self.log_indep_per_g.resize(n_genotypes, 0.0);
        self.log_prior_per_g.resize(n_genotypes, 0.0);
        self.log_post_unnorm.resize(n_genotypes, 0.0);
        self.log_post_unnorm_lane
            .resize(n_genotypes, wide::f64x4::ZERO);
        self.expected_counts.resize(n_alleles, 0.0);

        self.fixation_indices.resize(n_samples, 0.0);
        self.log_f_per_sample.resize(n_samples, 0.0);
        self.log_one_minus_f_per_sample.resize(n_samples, 0.0);
        self.allele_classes.resize(n_alleles, AlleleClass::Ref);
        self.pseudocounts.resize(n_alleles, 0.0);
        self.compound_mask.resize(n_alleles, false);
        self.mixture_contam_class_per_allele.resize(n_alleles, None);

        self.p_hat.resize(n_alleles, 0.0);
        self.p_hat_next.resize(n_alleles, 0.0);
        self.f_hat_compound.resize(n_alleles, 0.0);
        self.posteriors.resize(n_samples * n_genotypes, 0.0);

        self.best_genotype.resize(n_samples, 0);
        self.gq_phred.resize(n_samples, 0.0);
        self.compound_frequencies.resize(n_alleles, None);

        // Exact-AF QUAL buffers. `(ploidy + 1)` slots per sample for
        // the collapsed-by-alt-count likelihoods; `ploidy * n_samples
        // + 1` slots for the rolling AC log-probability vectors.
        let p = ploidy as usize;
        self.per_sample_alt_log_lik
            .resize(n_samples * (p + 1), f64::NEG_INFINITY);
        let ac_len = n_samples * p + 1;
        self.log_p_ac_curr.resize(ac_len, f64::NEG_INFINITY);
        self.log_p_ac_next.resize(ac_len, f64::NEG_INFINITY);
        // Linear-domain convolution buffers: `p` leading `0.0` pad slots
        // (so the FIR's `curr[k-c]` read is branch-free at the left edge)
        // + the AC range.
        let lin_len = p + ac_len;
        self.log_p_ac_curr_lin.resize(lin_len, 0.0);
        self.log_p_ac_next_lin.resize(lin_len, 0.0);

        // Mixture buffers are only used when contam-on; resize
        // unconditionally so the buffer is the right shape if a
        // subsequent record toggles contam state. (`Vec::resize` is
        // O(0) when length == capacity == requested.)
        self.mixture_log_likelihoods
            .resize(n_samples * n_genotypes, 0.0);
        self.mixture_n_obs.resize(n_alleles, 0);
        self.mixture_mean_err.resize(n_alleles, 0.0);
        self.mixture_c_s_all.resize(n_samples, 0.0);
        self.mixture_mean_err_all.resize(n_samples * n_alleles, 0.0);
    }
}

fn run_em_for_record<M: MathBackend>(
    record: MergedRecord,
    config: &PosteriorEngineConfig,
    math: &M,
    scratch: &mut RecordScratch,
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

    let view = MergedAllelesView::new(&alleles);
    let inputs = EmInputs {
        locus,
        ploidy,
        n_samples,
        n_genotypes,
        alleles: &view,
        scalars: &scalars,
        log_likelihoods: &log_likelihoods,
        chain_anchor_flags_for_validation: &chain_anchor_flags,
    };
    let em_outputs = run_em_columnar(inputs, config, math, scratch)?;

    let mut record = PosteriorRecord {
        locus,
        alleles,
        ploidy,
        n_samples,
        n_genotypes: em_outputs.n_genotypes,
        allele_frequencies: em_outputs.allele_frequencies,
        compound_frequencies: em_outputs.compound_frequencies,
        posteriors: em_outputs.posteriors,
        best_genotype: em_outputs.best_genotype,
        gq_phred: em_outputs.gq_phred,
        qual_phred: em_outputs.qual_phred,
        scalars,
        other_scalars,
        chain_anchor_flags,
        diagnostics: em_outputs.diagnostics,
        // Set later by the paralog write pass (if the filter runs and scores it).
        paralog_posterior: None,
    };
    // Drop candidate ALT alleles no sample's argmax genotype uses, so
    // the VCF never carries AC=0 alleles (see the method docs). No-op
    // for biallelic records.
    record.prune_unsupported_alleles(config.max_gq_phred);
    Ok(record)
}

/// Column-native EM entry point — Phase A.2 step 1. Same body as
/// [`run_em_for_record`] but reads its inputs by borrow through
/// [`EmInputs`] + [`AllelesView`]. Used by the row-shape shim
/// (above) and, after Phase A.2 step 2, directly by the cohort
/// worker.
///
/// **Returns:** [`EmOutputs`] — the EM-specific fields. The caller
/// joins these with the per-record allele bytes + row-shape
/// `scalars` / `other_scalars` / `chain_anchor_flags` to form the
/// emitted [`PosteriorRecord`].
pub(crate) fn run_em_columnar<M: MathBackend>(
    inputs: EmInputs<'_>,
    config: &PosteriorEngineConfig,
    math: &M,
    scratch: &mut RecordScratch,
) -> Result<EmOutputs, PosteriorEngineError> {
    let EmInputs {
        locus,
        ploidy,
        n_samples,
        n_genotypes,
        alleles,
        scalars,
        log_likelihoods,
        chain_anchor_flags_for_validation,
    } = inputs;
    let n_alleles = alleles.len();

    // Belt-and-braces: Stage 5 drops REF-only groups (n_alleles < 2).
    // If one slips through, emit a trivial posterior per the plan
    // §"Edge cases" — every sample is hom-REF with QUAL 0. (Checked
    // before the mixture pre-pass and `scratch.resize_to` since the
    // trivial path doesn't touch scratch.)
    if n_alleles < 2 {
        return Ok(trivial_em_outputs(
            alleles,
            ploidy,
            n_samples,
            config.max_gq_phred,
        ));
    }

    validate_record_shape(
        locus,
        ploidy,
        n_alleles,
        n_samples,
        n_genotypes,
        log_likelihoods,
        chain_anchor_flags_for_validation,
    )?;

    // Resize the per-engine scratch to the current record's shape
    // before any further field accesses. Steady-state records of the
    // same shape allocate nothing here.
    scratch.resize_to(n_samples, n_alleles, n_genotypes, ploidy);

    // Fill the record-static scratch fields (no mallocs — overwrites
    // pre-sized buffers).
    resolve_fixation_indices(
        locus,
        n_samples,
        config.fixation_index_overrides.as_deref(),
        config.fixation_index_default,
        &mut scratch.fixation_indices,
    )?;

    // `safe_ln(f_s)` and `safe_ln(1 - f_s)` are functions of
    // `fixation_indices` only, so they're record-static. Precompute
    // once here instead of inside `e_step`'s per-iteration sample loop.
    for s in 0..n_samples {
        let f_s = scratch.fixation_indices[s];
        scratch.log_f_per_sample[s] = safe_ln(math, f_s);
        scratch.log_one_minus_f_per_sample[s] = safe_ln(math, 1.0 - f_s);
    }

    // Detect the homogeneous-fixation case (every sample shares the
    // same `f_s`). Default config — `fixation_index_default = 0.0`,
    // `fixation_index_overrides = None` — hits this path; a CLI
    // `--inbreeding` global override also hits it. The
    // heterogeneous case (per-sample overrides differing) stays on
    // the per-sample path. Comparing on `f_s` rather than on the
    // logged value sidesteps `NaN != NaN` (and `-INF == -INF` is
    // `true` for any IEEE-754 implementation, but the typed input
    // here is always finite).
    let homogeneous_fixation = scratch
        .fixation_indices
        .first()
        .is_some_and(|&first| scratch.fixation_indices.iter().all(|&f| f == first));

    // Per-allele classification, pseudocounts, compound mask, and
    // contamination class — all record-static. Fill into scratch so
    // the EM loop and the mixture pre-pass read them by slice rather
    // than rebuilding per call.
    let ref_len = alleles.ref_len();
    for a in 0..n_alleles {
        let class = classify_allele_via_view(alleles, ref_len, a);
        scratch.allele_classes[a] = class;
        scratch.pseudocounts[a] = pseudocount_for(class, config);
        scratch.compound_mask[a] = alleles.is_compound(a);
        scratch.mixture_contam_class_per_allele[a] = map_to_contam_class(class);
    }

    // General SFS genotype prior (Dirichlet-multinomial): when a cohort
    // diversity θ̂ is configured, derive this record's concentration
    // α = [1, θ̂/(n_alleles−1), …] (the in-place equivalent of
    // `crate::genetics::alpha_from_diversity`, avoiding a per-record alloc) plus
    // the two tables the E-step reads: `lgamma(α_a)` (the Dirichlet-multinomial
    // independent-term baseline) and `log(α_a/Σα)` (the IBD marginal allele
    // frequency the Wright-`F` mixture consumes in place of `log p̂_a`). Compound
    // alleles are treated as ordinary rare ALTs here — they take the same
    // θ̂/(n_alleles−1) share; only the likelihood (already baked upstream) and AF
    // reporting distinguish them.
    let use_dirichlet_prior = if let Some(theta) = config.nucleotide_diversity {
        scratch.alpha[0] = crate::genetics::ALPHA_REF;
        let n_alt = n_alleles - 1;
        if n_alt > 0 {
            let per_alt = (theta / n_alt as f64).max(crate::genetics::MIN_ALT_CONCENTRATION);
            for a in 1..n_alleles {
                scratch.alpha[a] = per_alt;
            }
        }
        let sum_alpha: f64 = scratch.alpha[..n_alleles].iter().sum();
        let log_sum_alpha = math.ln(sum_alpha);
        for a in 0..n_alleles {
            scratch.lgamma_alpha[a] = crate::genetics::lgamma(scratch.alpha[a]);
            scratch.log_alpha_over_sum[a] = math.ln(scratch.alpha[a]) - log_sum_alpha;
        }
        true
    } else {
        false
    };

    // Genotype-shape artefacts (allele-count table, multinomial
    // coefficients, non-zero-pair and homozygous-allele lookups) come
    // from a thread-local `(ploidy, n_alleles)`-keyed cache. Repeated
    // shapes (the common case — every biallelic-diploid record shares
    // `(2, 2)`) hit the cache instead of rebuilding.
    let shape = shape_for(ploidy, n_alleles);

    // Mixture-likelihood branch: when contamination is configured,
    // fill `scratch.mixture_log_likelihoods` with the c_s≠0 version,
    // then `mem::take` it into a local owner so the EM body can read
    // a slice that does NOT alias `scratch`. We restore the buffer
    // before returning so the next record reuses the capacity.
    let mut local_mixture_ll: Vec<f64> = Vec::new();
    let did_mixture = if let Some(estimates) = config.contamination.as_ref() {
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
        // `M::HAS_LANE_4` is a compile-time const; the unused branch
        // is dead-code-eliminated after monomorphisation.
        if M::HAS_LANE_4 {
            compute_mixture_log_likelihoods_simd(
                locus,
                &shape,
                n_alleles,
                n_genotypes,
                n_samples,
                ploidy,
                scalars,
                estimates,
                log_likelihoods,
                math,
                scratch,
            )?;
        } else {
            compute_mixture_log_likelihoods(
                locus,
                &shape,
                n_alleles,
                n_genotypes,
                n_samples,
                ploidy,
                scalars,
                estimates,
                log_likelihoods,
                math,
                scratch,
            )?;
        }
        // Move the mixture buffer out of scratch so the EM body can
        // read it via a slice that does NOT alias `scratch`.
        local_mixture_ll = std::mem::take(&mut scratch.mixture_log_likelihoods);
        true
    } else {
        false
    };
    let ll_for_em: &[f64] = if did_mixture {
        &local_mixture_ll
    } else {
        log_likelihoods
    };

    // SFS-marginal genotype prior: used only for biallelic-diploid records
    // (`ploidy == 2, n_alleles == 2`). For those, the per-sample genotype prior
    // is the fixed table `config.sfs_prior_tables[sample]`; the HWE(p̂) path is
    // bypassed. Other shapes keep the HWE path even when the prior is
    // configured. The cohort-length invariant is checked here where `n_samples`
    // is known.
    let sfs_prior_tables: Option<&[[f64; 3]]> = match config.sfs_prior_tables.as_deref() {
        Some(tables) if ploidy == 2 && n_alleles == 2 => {
            if tables.len() != n_samples {
                return Err(PosteriorEngineError::SfsPriorTableCohortSizeMismatch {
                    locus,
                    table_samples: tables.len(),
                    record_samples: n_samples,
                });
            }
            Some(tables)
        }
        _ => None,
    };

    let ctx = EmContext {
        locus,
        n_samples,
        n_genotypes,
        n_alleles,
        ploidy,
        shape: &shape,
        homogeneous_fixation,
        compound_pseudocount: config.compound_alt_pseudocount,
        sfs_prior_tables,
        use_dirichlet_prior,
    };

    // EM initial point: uniform allele frequencies and uniform
    // compound frequencies. Overwrites scratch in place — no malloc.
    let p_init = 1.0 / n_alleles as f64;
    for slot in scratch.p_hat.iter_mut() {
        *slot = p_init;
    }
    for slot in scratch.f_hat_compound.iter_mut() {
        *slot = p_init;
    }

    let diagnostics = run_em_loop(ctx, config, math, ll_for_em, scratch)?;

    // Final E-step under the converged (p̂, f̂_C) so the emitted
    // posteriors reflect the *post-final-M-step* parameters rather
    // than the parameters that produced the last `posteriors` buffer
    // contents.
    //
    // `M::HAS_LANE_4` is a compile-time const so the unused branch is
    // dead-code-eliminated after monomorphisation. The SFS-marginal prior takes
    // its own biallelic-diploid E-step (fixed per-sample table, no p̂).
    if ctx.sfs_prior_tables.is_some() {
        e_step_sfs_biallelic(ctx, math, ll_for_em, scratch)?;
    } else if M::HAS_LANE_4 {
        e_step_simd(ctx, math, ll_for_em, scratch)?;
    } else {
        e_step(ctx, math, ll_for_em, scratch)?;
    }

    summarise_posteriors(n_samples, n_genotypes, scratch, config.max_gq_phred);
    // Beta(α_alt, α_ref) prior shape for the AC marginalisation —
    // read off the per-allele pseudocounts the EM already routed
    // through `pseudocount_for` into `scratch.pseudocounts`. Sum over
    // ALTs (a ≥ 1) for the "any non-ref" collapse.
    let alpha_ref = scratch.pseudocounts[0];
    let alpha_alt: f64 = scratch.pseudocounts[1..n_alleles].iter().sum();
    let qual_phred = compute_qual_via_exact_af(
        math,
        n_samples,
        n_genotypes,
        n_alleles,
        ploidy,
        ll_for_em,
        &shape,
        alpha_ref,
        alpha_alt,
        scratch,
    );

    // Restore the mixture buffer (if we took it) so the next record's
    // `resize_to` + mixture pre-pass reuse the high-water-mark
    // capacity instead of reallocating from empty.
    if did_mixture {
        scratch.mixture_log_likelihoods = local_mixture_ll;
    }

    // Fill `scratch.compound_frequencies` from the converged f̂_C and
    // the compound mask, then clone-out below. Reuses the scratch
    // slot instead of `.collect()`ing a fresh Vec.
    for a in 0..n_alleles {
        scratch.compound_frequencies[a] = if scratch.compound_mask[a] {
            Some(scratch.f_hat_compound[a])
        } else {
            None
        };
    }

    Ok(EmOutputs {
        n_genotypes,
        // Clone the converged outputs out of scratch. The clone is a
        // single malloc-then-memcpy per Vec — cheaper than the
        // pre-lift malloc-then-zero-init since memcpy avoids the
        // zero-init pass for an already-initialised buffer.
        allele_frequencies: scratch.p_hat.clone(),
        compound_frequencies: scratch.compound_frequencies.clone(),
        posteriors: scratch.posteriors.clone(),
        best_genotype: scratch.best_genotype.clone(),
        gq_phred: scratch.gq_phred.clone(),
        qual_phred,
        diagnostics,
    })
}

/// Allele-classification helper that reads through the
/// [`AllelesView`] trait — the column-native equivalent of
/// [`classify_allele`]. Identical decision tree.
fn classify_allele_via_view(
    alleles: &dyn AllelesView,
    ref_len: usize,
    allele_idx: usize,
) -> AlleleClass {
    if allele_idx == 0 {
        AlleleClass::Ref
    } else if alleles.is_compound(allele_idx) {
        AlleleClass::CompoundAlt
    } else if alleles.seq_len(allele_idx) == ref_len {
        AlleleClass::SnpAlt
    } else {
        AlleleClass::IndelAlt
    }
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

fn run_em_loop<M: MathBackend>(
    ctx: EmContext<'_>,
    config: &PosteriorEngineConfig,
    math: &M,
    log_likelihoods: &[f64],
    scratch: &mut RecordScratch,
) -> Result<EmDiagnostics, PosteriorEngineError> {
    let mut iterations: u32 = 0;
    let mut last_delta = f64::INFINITY;
    let max_iterations = config.max_iterations;

    while iterations < max_iterations {
        iterations += 1;

        if ctx.sfs_prior_tables.is_some() {
            e_step_sfs_biallelic(ctx, math, log_likelihoods, scratch)?;
        } else if M::HAS_LANE_4 {
            e_step_simd(ctx, math, log_likelihoods, scratch)?;
        } else {
            e_step(ctx, math, log_likelihoods, scratch)?;
        }

        // M-step writes `p_hat_next` and overwrites `f_hat_compound`
        // in place. `p_hat` is double-buffered against `p_hat_next`
        // so the convergence test can read both before we advance.
        m_step_p_hat(ctx, scratch);
        m_step_f_hat_compound(ctx, scratch);

        last_delta = max_abs_diff(&scratch.p_hat, &scratch.p_hat_next);
        std::mem::swap(&mut scratch.p_hat, &mut scratch.p_hat_next);

        if last_delta < config.convergence_threshold {
            return Ok(EmDiagnostics {
                iterations,
                final_max_delta_p: last_delta,
                converged: true,
            });
        }
    }

    // Iteration cap reached. Emit the record with `converged: false`
    // so the caller (and ultimately the VCF writer) can flag it via
    // `FILTER=EMNoConv` rather than hard-erroring the whole cohort
    // run on a single problematic site. The un-converged p̂/f̂
    // values flow through the same post-loop final E-step the
    // converged path uses.
    Ok(EmDiagnostics {
        iterations: max_iterations,
        final_max_delta_p: last_delta,
        converged: false,
    })
}

/// Fill `out` with the per-sample fixation indices `f_s`, either
/// from the per-sample overrides (when configured) or by broadcasting
/// the cohort-wide default. `out` is assumed pre-sized to `n_samples`
/// by `RecordScratch::resize_to`.
fn resolve_fixation_indices(
    locus: RecordLocus,
    n_samples: usize,
    overrides: Option<&[f64]>,
    default: f64,
    out: &mut [f64],
) -> Result<(), PosteriorEngineError> {
    debug_assert_eq!(out.len(), n_samples);
    match overrides {
        Some(overrides) => {
            if overrides.len() != n_samples {
                return Err(PosteriorEngineError::InbreedingOverrideLengthMismatch {
                    locus,
                    override_len: overrides.len(),
                    record_samples: n_samples,
                });
            }
            out.copy_from_slice(overrides);
        }
        None => {
            for slot in out.iter_mut() {
                *slot = default;
            }
        }
    }
    Ok(())
}

fn e_step<M: MathBackend>(
    ctx: EmContext<'_>,
    math: &M,
    log_likelihoods: &[f64],
    scratch: &mut RecordScratch,
) -> Result<(), PosteriorEngineError> {
    // Effective per-allele frequency seen by the prior: `f̂_C`
    // substitutes for `p̂[compound]` in every compound slot, regardless
    // of any per-sample chain-broken flag (the chain-broken flag only
    // affects the likelihood, which Stage 5 has already baked into
    // `log_likelihoods`).
    for a in 0..ctx.n_alleles {
        let p = if scratch.compound_mask[a] {
            scratch.f_hat_compound[a]
        } else {
            scratch.p_hat[a]
        };
        scratch.p_effective[a] = p;
        // Under the Dirichlet-multinomial prior the frequency the mixture sees is
        // the IBD marginal `log(α_a/Σα)` (precomputed, p̂-independent), not the
        // plug-in `log p̂_a`.
        scratch.log_p_effective[a] = if ctx.use_dirichlet_prior {
            scratch.log_alpha_over_sum[a]
        } else {
            safe_ln(math, p)
        };
    }

    // H2: `log_indep` depends only on `log_p_effective` and the
    // genotype shape — sample-invariant within an EM iteration.
    // Precompute once and read inside the sample loop instead of
    // rebuilding per (sample, g_idx).
    fill_log_indep_per_g(ctx, scratch);

    // H4: when every sample shares the same `f_s`, the whole
    // per-genotype log-prior is also sample-invariant; fold it into
    // `scratch.log_prior_per_g` once and let the sample loop splat.
    // Heterogeneous-fixation records take the per-sample path below.
    if ctx.homogeneous_fixation {
        fill_log_prior_per_g_homogeneous(ctx, math, scratch);
    }

    let n_genotypes = ctx.n_genotypes;
    for sample_idx in 0..ctx.n_samples {
        let ll_row = &log_likelihoods[sample_idx * n_genotypes..(sample_idx + 1) * n_genotypes];

        if ctx.homogeneous_fixation {
            // log_prior_per_g[g] already carries the IBD term; just
            // add the sample's likelihood per genotype.
            for (g_idx, &ll) in ll_row.iter().enumerate() {
                scratch.log_post_unnorm[g_idx] = ll + scratch.log_prior_per_g[g_idx];
            }
        } else {
            // Heterogeneous: `log_f` / `log_one_minus_f` change per
            // sample, so the IBD branch has to run per-cell.
            let log_one_minus_f = scratch.log_one_minus_f_per_sample[sample_idx];
            let log_f = scratch.log_f_per_sample[sample_idx];

            for (g_idx, &ll) in ll_row.iter().enumerate() {
                let log_indep = scratch.log_indep_per_g[g_idx];
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

        let post_base = sample_idx * n_genotypes;
        for g_idx in 0..n_genotypes {
            let posterior = math.exp(scratch.log_post_unnorm[g_idx] - log_z);
            if !posterior.is_finite() {
                return Err(PosteriorEngineError::NonFinitePosterior {
                    locus: ctx.locus,
                    sample_idx,
                    genotype_idx: Some(g_idx),
                    kind: classify_nonfinite(posterior),
                });
            }
            scratch.posteriors[post_base + g_idx] = posterior;
        }
    }

    Ok(())
}

/// E-step for the biallelic-diploid **SFS-marginal** genotype prior.
///
/// The per-sample genotype log-prior is a fixed table
/// `[ln P(0/0), ln P(0/1), ln P(1/1)]` (from `ctx.sfs_prior_tables`) that does
/// **not** depend on the EM's allele-frequency estimate `p̂` — so the posterior
/// is just `likelihood × prior`, normalised, with no HWE(p̂) construction. Only
/// reached when `ctx.sfs_prior_tables` is `Some`, which the caller sets only for
/// `ploidy == 2, n_alleles == 2` (three genotypes in canonical order
/// `0/0, 0/1, 1/1`, matching the table's index order).
fn e_step_sfs_biallelic<M: MathBackend>(
    ctx: EmContext<'_>,
    math: &M,
    log_likelihoods: &[f64],
    scratch: &mut RecordScratch,
) -> Result<(), PosteriorEngineError> {
    let tables = ctx
        .sfs_prior_tables
        .expect("e_step_sfs_biallelic called without sfs_prior_tables");
    let n_genotypes = ctx.n_genotypes;
    debug_assert_eq!(
        n_genotypes, 3,
        "SFS-marginal E-step is biallelic-diploid only"
    );
    debug_assert_eq!(tables.len(), ctx.n_samples);

    for sample_idx in 0..ctx.n_samples {
        let ll_row = &log_likelihoods[sample_idx * n_genotypes..(sample_idx + 1) * n_genotypes];
        let prior = &tables[sample_idx];
        for g_idx in 0..n_genotypes {
            scratch.log_post_unnorm[g_idx] = ll_row[g_idx] + prior[g_idx];
        }

        let log_z = log_sum_exp_slice(math, &scratch.log_post_unnorm[..n_genotypes]);
        if !log_z.is_finite() {
            return Err(PosteriorEngineError::NonFinitePosterior {
                locus: ctx.locus,
                sample_idx,
                genotype_idx: None,
                kind: classify_nonfinite(log_z),
            });
        }

        let post_base = sample_idx * n_genotypes;
        for g_idx in 0..n_genotypes {
            let posterior = math.exp(scratch.log_post_unnorm[g_idx] - log_z);
            if !posterior.is_finite() {
                return Err(PosteriorEngineError::NonFinitePosterior {
                    locus: ctx.locus,
                    sample_idx,
                    genotype_idx: Some(g_idx),
                    kind: classify_nonfinite(posterior),
                });
            }
            scratch.posteriors[post_base + g_idx] = posterior;
        }
    }

    Ok(())
}

/// SIMD-aware variant of [`e_step`]. Processes samples in batches of
/// 4 lanes; falls back to the scalar `e_step` body for the tail when
/// `n_samples` is not a multiple of 4.
///
/// The dispatch in [`run_em_loop`] / [`run_em_for_record`] selects
/// this when the math backend's `HAS_LANE_4` const is `true`. For
/// scalar backends the const is `false` and this function is
/// dead-code-eliminated by monomorphisation.
fn e_step_simd<M: MathBackend>(
    ctx: EmContext<'_>,
    math: &M,
    log_likelihoods: &[f64],
    scratch: &mut RecordScratch,
) -> Result<(), PosteriorEngineError> {
    // `log_p_effective` build is per-allele scalar work — identical
    // to the scalar e_step.
    for a in 0..ctx.n_alleles {
        let p = if scratch.compound_mask[a] {
            scratch.f_hat_compound[a]
        } else {
            scratch.p_hat[a]
        };
        scratch.p_effective[a] = p;
        // See `e_step`: the Dirichlet-multinomial prior feeds `log(α_a/Σα)`, not
        // `log p̂_a`, to the mixture.
        scratch.log_p_effective[a] = if ctx.use_dirichlet_prior {
            scratch.log_alpha_over_sum[a]
        } else {
            safe_ln(math, p)
        };
    }

    // H2: cache `log_indep[g]` once per iteration; the batch loop
    // splats from it instead of recomputing per (batch, g_idx).
    fill_log_indep_per_g(ctx, scratch);

    // H4: when fixation is homogeneous (default config — `f_s = 0`
    // for every sample), the whole per-genotype log-prior is
    // sample-invariant. Fold it once into `scratch.log_prior_per_g`
    // and let the batch loop splat from there instead of doing a
    // 3-transcendental `log_sum_exp_2_x4` per (batch, homozygous g)
    // whose `log_f_v = -∞` lane makes the work redundant anyway.
    if ctx.homogeneous_fixation {
        fill_log_prior_per_g_homogeneous(ctx, math, scratch);
    }

    let n_genotypes = ctx.n_genotypes;
    let n_samples = ctx.n_samples;
    let n_full_batches = n_samples / 4;
    let tail_start = n_full_batches * 4;

    // Batch loop: 4 samples per iteration.
    for batch_idx in 0..n_full_batches {
        let s0 = batch_idx * 4;

        // `log_f_v` / `log_one_minus_f_v` only matter on the
        // heterogeneous-fixation path; the homogeneous case reads
        // its log-prior straight from `scratch.log_prior_per_g`.
        let (log_f_v, log_one_minus_f_v) = if ctx.homogeneous_fixation {
            (wide::f64x4::ZERO, wide::f64x4::ZERO)
        } else {
            (
                wide::f64x4::from([
                    scratch.log_f_per_sample[s0],
                    scratch.log_f_per_sample[s0 + 1],
                    scratch.log_f_per_sample[s0 + 2],
                    scratch.log_f_per_sample[s0 + 3],
                ]),
                wide::f64x4::from([
                    scratch.log_one_minus_f_per_sample[s0],
                    scratch.log_one_minus_f_per_sample[s0 + 1],
                    scratch.log_one_minus_f_per_sample[s0 + 2],
                    scratch.log_one_minus_f_per_sample[s0 + 3],
                ]),
            )
        };

        // Compute `log_post_unnorm[g]` as a `f64x4` across the 4
        // samples for every genotype.
        for g_idx in 0..n_genotypes {
            let log_prior_v = if ctx.homogeneous_fixation {
                wide::f64x4::splat(scratch.log_prior_per_g[g_idx])
            } else {
                let log_indep_v = wide::f64x4::splat(scratch.log_indep_per_g[g_idx]);
                match ctx.shape.homozygous_allele_for[g_idx] {
                    Some(homo_allele) => {
                        let term1 = log_one_minus_f_v + log_indep_v;
                        let term2 = log_f_v
                            + wide::f64x4::splat(scratch.log_p_effective[homo_allele as usize]);
                        log_sum_exp_2_x4(math, term1, term2)
                    }
                    None => log_one_minus_f_v + log_indep_v,
                }
            };

            // Lane-load the per-sample likelihood for this genotype.
            // `log_likelihoods` is `[s * n_genotypes + g]` row-major.
            let ll_v = wide::f64x4::from([
                log_likelihoods[s0 * n_genotypes + g_idx],
                log_likelihoods[(s0 + 1) * n_genotypes + g_idx],
                log_likelihoods[(s0 + 2) * n_genotypes + g_idx],
                log_likelihoods[(s0 + 3) * n_genotypes + g_idx],
            ]);

            scratch.log_post_unnorm_lane[g_idx] = ll_v + log_prior_v;
        }

        // Per-lane `log_z` via SIMD log_sum_exp_slice.
        let log_z_v = log_sum_exp_slice_x4(math, &scratch.log_post_unnorm_lane);
        let log_z_arr = log_z_v.to_array();

        // Surface non-finite `log_z` per lane.
        for (lane, &z) in log_z_arr.iter().enumerate() {
            if !z.is_finite() {
                return Err(PosteriorEngineError::NonFinitePosterior {
                    locus: ctx.locus,
                    sample_idx: s0 + lane,
                    genotype_idx: None,
                    kind: classify_nonfinite(z),
                });
            }
        }

        // Normalise into `scratch.posteriors`. Layout stays
        // sample-major, so we scatter each lane back to its
        // `(sample, genotype)` slot.
        for g_idx in 0..n_genotypes {
            let post_v = math.exp_x4(scratch.log_post_unnorm_lane[g_idx] - log_z_v);
            let post_arr = post_v.to_array();
            for (lane, &p) in post_arr.iter().enumerate() {
                if !p.is_finite() {
                    return Err(PosteriorEngineError::NonFinitePosterior {
                        locus: ctx.locus,
                        sample_idx: s0 + lane,
                        genotype_idx: Some(g_idx),
                        kind: classify_nonfinite(p),
                    });
                }
                scratch.posteriors[(s0 + lane) * n_genotypes + g_idx] = p;
            }
        }
    }

    // Tail: process any samples beyond the last full batch via the
    // scalar inner loop. Reads from the same `log_indep_per_g` /
    // `log_prior_per_g` scratch the batch loop populated above, so
    // no per-tail-sample recomputation.
    for sample_idx in tail_start..n_samples {
        let ll_row = &log_likelihoods[sample_idx * n_genotypes..(sample_idx + 1) * n_genotypes];

        if ctx.homogeneous_fixation {
            for (g_idx, &ll) in ll_row.iter().enumerate() {
                scratch.log_post_unnorm[g_idx] = ll + scratch.log_prior_per_g[g_idx];
            }
        } else {
            let log_one_minus_f = scratch.log_one_minus_f_per_sample[sample_idx];
            let log_f = scratch.log_f_per_sample[sample_idx];

            for (g_idx, &ll) in ll_row.iter().enumerate() {
                let log_indep = scratch.log_indep_per_g[g_idx];
                let log_prior = match ctx.shape.homozygous_allele_for[g_idx] {
                    Some(homo_allele) => log_sum_exp_2(
                        math,
                        log_one_minus_f + log_indep,
                        log_f + scratch.log_p_effective[homo_allele as usize],
                    ),
                    None => log_one_minus_f + log_indep,
                };
                scratch.log_post_unnorm[g_idx] = ll + log_prior;
            }
        }

        let log_z = log_sum_exp_slice(math, &scratch.log_post_unnorm);
        if !log_z.is_finite() {
            return Err(PosteriorEngineError::NonFinitePosterior {
                locus: ctx.locus,
                sample_idx,
                genotype_idx: None,
                kind: classify_nonfinite(log_z),
            });
        }
        let post_base = sample_idx * n_genotypes;
        for g_idx in 0..n_genotypes {
            let posterior = math.exp(scratch.log_post_unnorm[g_idx] - log_z);
            if !posterior.is_finite() {
                return Err(PosteriorEngineError::NonFinitePosterior {
                    locus: ctx.locus,
                    sample_idx,
                    genotype_idx: Some(g_idx),
                    kind: classify_nonfinite(posterior),
                });
            }
            scratch.posteriors[post_base + g_idx] = posterior;
        }
    }

    Ok(())
}

/// M-step on p̂. Writes the new estimate into `scratch.p_hat_next`;
/// the caller (`run_em_loop`) then computes `max_abs_diff(p_hat,
/// p_hat_next)` and `mem::swap`s the two buffers to advance the
/// iteration.
fn m_step_p_hat(ctx: EmContext<'_>, scratch: &mut RecordScratch) {
    accumulate_expected_counts(ctx, &scratch.posteriors, &mut scratch.expected_counts);

    // M-step on p̂ — Dirichlet posterior mean over the simplex.
    let dirichlet_denominator: f64 = scratch
        .expected_counts
        .iter()
        .zip(scratch.pseudocounts.iter())
        .map(|(e, a)| e + a)
        .sum();
    for a in 0..ctx.n_alleles {
        scratch.p_hat_next[a] =
            (scratch.expected_counts[a] + scratch.pseudocounts[a]) / dirichlet_denominator;
    }
}

/// Sum `posterior[s, g] · gt_counts[g, a]` over (s, g) into
/// `expected_counts[a]`, clearing the output first. Split out from
/// [`m_step_p_hat`] so the biallelic-diploid fast path can specialise
/// the inner triple loop without cluttering the surrounding Dirichlet
/// arithmetic.
fn accumulate_expected_counts(ctx: EmContext<'_>, posteriors: &[f64], expected_counts: &mut [f64]) {
    // Biallelic-diploid is the dominant production case: every
    // standard SNP and small indel hits `(n_alleles, n_genotypes) ==
    // (2, 3)`. `genotype_allele_counts` is then statically
    // `[2, 0, 1, 1, 0, 2]` (RR, RA, AA), so the per-sample inner
    // loop reduces to two FMAs with no data-dependent branch.
    if ctx.n_alleles == 2 && ctx.n_genotypes == 3 {
        let mut e_ref = 0.0_f64;
        let mut e_alt = 0.0_f64;
        for post_row in posteriors.chunks_exact(3) {
            // `chunks_exact(3)` already proves `post_row.len() == 3`,
            // so the indexed loads compile down to bounds-check-free
            // movsd instructions.
            let p_rr = post_row[0];
            let p_ra = post_row[1];
            let p_aa = post_row[2];
            // RR=[2,0] → contributes 2·P_RR to REF, 0 to ALT
            // RA=[1,1] → contributes   P_RA to REF,   P_RA to ALT
            // AA=[0,2] → contributes 0 to REF,      2·P_AA to ALT
            e_ref += 2.0 * p_rr + p_ra;
            e_alt += p_ra + 2.0 * p_aa;
        }
        expected_counts[0] = e_ref;
        expected_counts[1] = e_alt;
        return;
    }

    // General path: any `(ploidy, n_alleles)` shape.
    for slot in expected_counts.iter_mut() {
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
                    expected_counts[a] += weight * k as f64;
                }
            }
        }
    }
}

/// M-step on f̂_C. Writes the new estimate into
/// `scratch.f_hat_compound`. Unlike `p_hat`, `f_hat_compound` is
/// not part of the convergence test, so it can be overwritten in
/// place without double-buffering.
fn m_step_f_hat_compound(ctx: EmContext<'_>, scratch: &mut RecordScratch) {
    // M-step on f̂_C — Beta(α_compound, 1 − α_compound) posterior
    // mean. Each compound is independent of the others.
    let chromosomes_total = ctx.ploidy as f64 * ctx.n_samples as f64;
    let alpha_other = 1.0 - ctx.compound_pseudocount;
    let denominator = ctx.compound_pseudocount + alpha_other + chromosomes_total;
    for a in 0..ctx.n_alleles {
        scratch.f_hat_compound[a] = if scratch.compound_mask[a] {
            (ctx.compound_pseudocount + scratch.expected_counts[a]) / denominator
        } else {
            0.0
        };
    }
}

/// Compute per-sample best genotype, GQ phred, and site-level QUAL
/// from the EM-converged posteriors.
///
/// Writes the per-sample outputs (`best_genotype`, `gq_phred`) into
/// `scratch` so the caller can `.clone()` them into the output
/// record without a fresh malloc-then-zero-init per record. Returns
/// only the scalar `qual_phred`.
/// Walk each sample's posterior row and fill the per-sample output
/// buffers (`best_genotype` argmax and `gq_phred`). The site-level
/// QUAL is no longer derived here — see [`compute_qual_via_exact_af`].
fn summarise_posteriors(
    n_samples: usize,
    n_genotypes: usize,
    scratch: &mut RecordScratch,
    max_gq: f64,
) {
    for sample_idx in 0..n_samples {
        let post_row =
            &scratch.posteriors[sample_idx * n_genotypes..(sample_idx + 1) * n_genotypes];
        let (best_idx, best_p) = post_row.iter().enumerate().fold(
            (0_usize, f64::NEG_INFINITY),
            |(best_i, best_v), (i, &v)| {
                if v > best_v { (i, v) } else { (best_i, best_v) }
            },
        );
        scratch.best_genotype[sample_idx] = best_idx;

        // Clamp p_best below 1 by one ULP so the Phred calculation is
        // finite when EM produced an exact 1.0 (e.g. a singleton
        // genotype after softmax of `[finite, -∞, -∞]`).
        let p_best_clamped = best_p.min(1.0 - f64::EPSILON);
        let gq = PHRED_SCALE * (1.0 - p_best_clamped).log10();
        scratch.gq_phred[sample_idx] = gq.clamp(0.0, max_gq);
    }
}

/// Lanczos approximation to `ln Γ(x)` for `x > 0`. Used by
/// [`compute_qual_via_exact_af`] to evaluate the Beta-Binomial AC
/// prior — argument range is `≈ 0.001..few-thousand` (Dirichlet
/// pseudocounts on one end, cohort factorials on the other), which
/// the standard `g = 7, n = 9` form covers at ~15-digit accuracy.
/// For `x < 0.5` the reflection formula `Γ(z)·Γ(1 − z) = π / sin(πz)`
/// is applied; in our use the argument is always `> 0` (positive
/// pseudocounts and factorial offsets), but the reflection branch is
/// kept for robustness.
fn ln_gamma(x: f64) -> f64 {
    debug_assert!(x > 0.0, "ln_gamma called with non-positive x={x}");
    const G: f64 = 7.0;
    const COEFFS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_2,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().abs().ln()
            - ln_gamma(1.0 - x);
    }
    let z = x - 1.0;
    let mut a = COEFFS[0];
    for (i, &c) in COEFFS.iter().enumerate().skip(1) {
        a += c / (z + i as f64);
    }
    let t = z + G + 0.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + a.ln() + (z + 0.5) * t.ln() - t
}

/// Site-level QUAL via "any non-ref carrier" exact AC enumeration.
///
/// Computes `−10·log10 P(K = 0 | data)` where `K = Σ_s c_s` is the
/// total non-ref allele count across the cohort, `c_s ∈ {0..ploidy}`
/// is sample `s`'s non-ref allele count, and the prior on `K` is the
/// Beta-Binomial induced by a Dirichlet on (ref, any-non-ref) allele
/// frequency with `α_ref = config.ref_pseudocount` and
/// `α_alt = Σ_{a ≥ 1} pseudocount_for(allele_classes[a])`. The
/// "any non-ref" collapse means multi-allelic sites use a single QUAL
/// (one allele frequency over the union of ALTs); this matches the
/// VCF intent of `−10·log10 P(no variant)`.
///
/// In contrast to the previous `Π_s P(GT_s = hom-ref)` formula, this
/// is a properly normalised marginal posterior — adding hom-ref
/// samples adds essentially zero log-evidence against `K > 0` under
/// the sparse-het prior, so QUAL is bounded by what the few
/// non-hom-ref samples actually justify rather than scaling with the
/// cohort size. See the QUAL-semantics project memory for the
/// alignment with GATK (`P(AC = 0 | data)`) and FreeBayes
/// (`1 - P(K = 1 | reads)`).
///
/// Algorithm (log-domain throughout, `log_sum_exp_2` for stability):
///   1. **Collapse** per-sample per-genotype log-likelihoods into
///      per-sample-per-`c` (= non-ref-allele-count) cells:
///      `l_s[c] = logsumexp_{g : non_ref(g) = c} log_likelihoods[s, g]`.
///   2. **Convolution** sample-by-sample over the cohort AC axis:
///      `log_p_ac[k] ← logsumexp_{c} (log_p_ac_prev[k − c] + l_s[c])`.
///   3. **Apply log-prior** `log P(K = k)` from
///      Beta-Binomial(`ploidy · n_samples`, `α_alt`, `α_ref`).
///   4. **Normalise** with logsumexp over `K = 0..=ploidy·n_samples`
///      and read off `P(K = 0 | data)`.
#[allow(clippy::too_many_arguments)] // shape, math backend, scratch, and the four numeric record shape parameters are all distinct concerns — bundling adds indirection without clarity gain.
fn compute_qual_via_exact_af<M: MathBackend>(
    math: &M,
    n_samples: usize,
    n_genotypes: usize,
    n_alleles: usize,
    ploidy: u8,
    log_likelihoods: &[f64],
    shape: &GenotypeShape,
    alpha_ref: f64,
    alpha_alt: f64,
    scratch: &mut RecordScratch,
) -> f64 {
    let p = ploidy as usize;
    let max_k = n_samples * p;
    let max_k_p1 = max_k + 1;

    // Step 1: per-sample collapsed log-likelihoods.
    for slot in scratch.per_sample_alt_log_lik.iter_mut() {
        *slot = f64::NEG_INFINITY;
    }
    for s in 0..n_samples {
        let ll_row = &log_likelihoods[s * n_genotypes..(s + 1) * n_genotypes];
        let alt_row_base = s * (p + 1);
        for (g, &ll) in ll_row.iter().enumerate() {
            let counts = &shape.genotype_allele_counts[g * n_alleles..(g + 1) * n_alleles];
            // c = ploidy − count of REF (allele 0) = total non-ref copies.
            let c = (p as u32 - counts[0]) as usize;
            let slot = &mut scratch.per_sample_alt_log_lik[alt_row_base + c];
            *slot = log_sum_exp_2(math, *slot, ll);
        }
    }

    // Step 2: convolution. The log-domain fold spent ~88 % of the
    // consume path's self-time in scalar `exp`/`ln` and is
    // O(n_samples²·ploidy), so it runs in the linear domain as a
    // transcendental-free multiply-add FIR with per-sample
    // renormalisation. The result is written back into
    // `log_p_ac_curr[..max_k_p1]` in the log domain, so Steps 3–4 are
    // unchanged.
    convolve_ac_linear(
        n_samples,
        p,
        max_k_p1,
        &scratch.per_sample_alt_log_lik,
        &mut scratch.log_p_ac_curr_lin,
        &mut scratch.log_p_ac_next_lin,
        &mut scratch.log_p_ac_curr,
    );

    // Step 3: apply the Beta-Binomial log-prior on `K`.
    //
    // The (ref, any-non-ref) collapse on the allele-frequency simplex
    // turns the per-allele Dirichlet pseudocounts into a Beta prior
    // on `f_non_ref` with shape (`α_alt`, `α_ref`), which in turn
    // induces Beta-Binomial(max_k, α_alt, α_ref) on `K`.
    let log_beta_norm = ln_gamma(alpha_alt) + ln_gamma(alpha_ref) - ln_gamma(alpha_alt + alpha_ref);
    let max_k_f = max_k as f64;
    let log_max_k_fact = ln_gamma(max_k_f + 1.0);
    let log_denom = ln_gamma(alpha_alt + alpha_ref + max_k_f);

    // Overwrite `log_p_ac_next` with the unnormalised log-posterior
    // (likelihood + prior) so we can sweep it for the partition once.
    for k in 0..max_k_p1 {
        let k_f = k as f64;
        let log_binom = log_max_k_fact - ln_gamma(k_f + 1.0) - ln_gamma(max_k_f - k_f + 1.0);
        let log_beta_kk = ln_gamma(alpha_alt + k_f) + ln_gamma(alpha_ref + max_k_f - k_f);
        let log_prior_k = log_binom + log_beta_kk - log_denom - log_beta_norm;
        scratch.log_p_ac_next[k] = scratch.log_p_ac_curr[k] + log_prior_k;
    }

    // Step 4: normalise via logsumexp; derive QUAL from P(K = 0).
    let log_z = log_sum_exp_slice(math, &scratch.log_p_ac_next[..max_k_p1]);
    if log_z == f64::NEG_INFINITY {
        // Pathological — every (K, data) combination has zero mass.
        // Safer to surface as QUAL = 0 than to claim infinite
        // confidence; this branch should not be reachable from real
        // inputs.
        return 0.0;
    }
    let log_p_k0_post = scratch.log_p_ac_next[0] - log_z;
    if log_p_k0_post == f64::NEG_INFINITY {
        return f64::INFINITY;
    }
    // QUAL = −10 · log10(P(K = 0 | data)).
    PHRED_SCALE * (log_p_k0_post / std::f64::consts::LN_10)
}

/// Step 2 of [`compute_qual_via_exact_af`] — fold each sample's
/// `(ploidy + 1)`-tap collapsed-likelihood kernel into the rolling
/// allele-count distribution — run in the **linear domain**.
///
/// The mathematically-equivalent log-domain fold
/// (`next[k] = log Σ_c exp(curr[k-c] + l_s[c])`) spends its whole
/// O(n_samples²·ploidy) inner loop in `exp`/`ln` (~88 % of the consume
/// path's self-time at N=200). Here the inner loop is instead the plain
/// linear convolution `next[k] = Σ_c curr[k-c] · Lₛ[c]` — a fixed-tap
/// multiply-add FIR with **no transcendentals**, which autovectorises
/// to FMA. The only transcendentals left are O(n_samples): one `exp`
/// per kernel tap and one `ln` per sample for the renormalisation.
///
/// **Numerical scheme.** A scaled forward algorithm. After folding `s`
/// samples the invariant is `P(AC = k) = cur[p + k] · exp(log_scale)`:
/// each sample's kernel is divided by its max (folded into `log_scale`)
/// and the post-fold distribution is renormalised by its max (also
/// folded in), keeping the live values in `(0, 1]`. The absolute
/// log-probabilities `ln(cur[p+k]) + log_scale` are recovered at the
/// end. Tail entries that underflow to `0.0` map to `-∞`; they are
/// negligible in the partition `Z` and are never read individually.
///
/// **K = 0 is tracked exactly** in the log domain (`log P(AC = 0) =
/// Σ_s l_s[0]`, a trivial recurrence). It is the QUAL numerator and can
/// decay far below the linear range for a confident large cohort, so
/// reading it back from the renormalised linear buffer would lose it to
/// underflow (and report QUAL = +∞). The exact value overrides
/// `out_log[0]` and stays consistent with the bulk because `log_scale`
/// is the shared true scale.
///
/// `cur`/`nxt` are the left-padded (`ploidy` leading `0.0` slots) linear
/// buffers; `out_log` receives the log-domain result in `[..max_k_p1]`.
/// Transcendentals use the exact `std` `exp`/`ln` (not the approximate
/// math backend) — they are off the hot path now, so exact costs
/// nothing and keeps QUAL backend-independent.
fn convolve_ac_linear<'a>(
    n_samples: usize,
    p: usize,
    max_k_p1: usize,
    per_sample_alt_log_lik: &[f64],
    mut cur: &'a mut [f64],
    mut nxt: &'a mut [f64],
    out_log: &mut [f64],
) {
    // Empty-prefix distribution: mass 1 on `K = 0` (logical index 0 ->
    // padded index `p`).
    cur.fill(0.0);
    cur[p] = 1.0;
    let mut log_scale = 0.0f64;
    let mut log_pk0 = 0.0f64; // exact log P(AC = 0)

    for s in 0..n_samples {
        let l_s = &per_sample_alt_log_lik[s * (p + 1)..(s + 1) * (p + 1)];
        let kmax_after = (s + 1) * p;

        // Exact K=0 recurrence — independent of the linear bulk.
        log_pk0 += l_s[0];

        // Linear kernel, normalised by its max so the taps land in
        // `(0, 1]`; the max folds back into the running log-scale.
        let kmax_s = l_s.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mut lker = [0.0f64; MAX_PLOIDY as usize + 1];
        for (c, &ll) in l_s.iter().enumerate().take(p + 1) {
            lker[c] = if ll > f64::NEG_INFINITY {
                (ll - kmax_s).exp()
            } else {
                0.0
            };
        }

        // FIR: next[k] = Σ_c cur[k-c] · lker[c], tap-major so each pass
        // is a contiguous `dst[k] += src[k] · lc` axpy that vectorises.
        let live = kmax_after + 1;
        nxt[..p + live].fill(0.0);
        for (c, &lc) in lker.iter().enumerate().take(p + 1) {
            if lc == 0.0 {
                continue;
            }
            let src = &cur[(p - c)..(p - c) + live];
            let dst = &mut nxt[p..p + live];
            for (d, &sv) in dst.iter_mut().zip(src.iter()) {
                *d += sv * lc;
            }
        }

        // Renormalise the live region by its max; fold both the kernel
        // max and the renorm factor into the running log-scale.
        log_scale += kmax_s;
        let m = nxt[p..p + live].iter().copied().fold(0.0f64, f64::max);
        if m > 0.0 {
            let inv = 1.0 / m;
            for v in nxt[p..p + live].iter_mut() {
                *v *= inv;
            }
            log_scale += m.ln();
        }

        std::mem::swap(&mut cur, &mut nxt);
    }

    // Recover absolute log P(AC = k); override K=0 with the exact value.
    for (k, slot) in out_log[..max_k_p1].iter_mut().enumerate() {
        let v = cur[p + k];
        *slot = if v > 0.0 {
            v.ln() + log_scale
        } else {
            f64::NEG_INFINITY
        };
    }
    out_log[0] = log_pk0;
}

/// EM-output fields for the trivial path — `n_alleles < 2`. Reads
/// the allele set through the [`AllelesView`] trait so both the
/// row-shape shim and the column-native worker share one body. The
/// caller joins this with the per-record allele bytes + row-shape
/// `scalars` / `other_scalars` / `chain_anchor_flags` to form the
/// trivial [`PosteriorRecord`].
fn trivial_em_outputs(
    alleles: &dyn AllelesView,
    ploidy: u8,
    n_samples: usize,
    max_gq: f64,
) -> EmOutputs {
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

    let allele_frequencies = if n_alleles == 0 {
        Vec::new()
    } else {
        let mut p = vec![0.0; n_alleles];
        p[0] = 1.0;
        p
    };
    let compound_frequencies = (0..n_alleles)
        .map(|a| {
            if alleles.is_compound(a) {
                Some(0.0)
            } else {
                None
            }
        })
        .collect();

    EmOutputs {
        n_genotypes,
        allele_frequencies,
        compound_frequencies,
        posteriors,
        best_genotype,
        gq_phred,
        qual_phred: 0.0,
        diagnostics: EmDiagnostics {
            iterations: 0,
            final_max_delta_p: 0.0,
            // Trivial record — the closed-form output IS the
            // converged fixed point, no iteration needed.
            converged: true,
        },
    }
}

// ---------------------------------------------------------------------
// E-step prelude helpers
// ---------------------------------------------------------------------

/// Fill `scratch.log_indep_per_g` — the random-mating (no-inbreeding) genotype
/// log-prior term that the Wright-`F` mixture then bends by each sample's `F`.
///
/// Two forms, selected by `ctx.use_dirichlet_prior`:
///
/// - **Dirichlet-multinomial** (SFS prior, `ctx.use_dirichlet_prior == true`):
///   `log_multinomial_coeffs[g] + Σ_a [lgamma(α_a+k_a) − lgamma(α_a)]` over the
///   genotype's non-zero `(allele, count)` pairs — the exact marginal of the
///   multinomial genotype prior over the Dirichlet frequency prior with
///   concentration `α` (`crate::genetics::dirichlet_multinomial_log_priors`,
///   computed here in place from `scratch.alpha` / `scratch.lgamma_alpha` to
///   avoid a per-record allocation). p̂-independent — constant across EM
///   iterations.
/// - **Plug-in HWE** (`ctx.use_dirichlet_prior == false`, the legacy path):
///   `log_multinomial_coeffs[g] + Σ k_a · log_p̂_a` — plugs the EM frequency
///   estimate straight into the genotype prior.
///
/// Sample-invariant within an EM iteration; both `e_step` and `e_step_simd` call
/// this once per iteration so the sample / batch loop reads from the cache
/// instead of rebuilding (H2 in `perf_posterior_engine_2026-05-18.md`).
#[inline]
fn fill_log_indep_per_g(ctx: EmContext<'_>, scratch: &mut RecordScratch) {
    for g_idx in 0..ctx.n_genotypes {
        let (start, len) = ctx.shape.nonzero_pairs_offsets[g_idx];
        let pairs = &ctx.shape.nonzero_pairs[start as usize..start as usize + len as usize];
        let sum: f64 = if ctx.use_dirichlet_prior {
            // Σ_a [lgamma(α_a + k_a) − lgamma(α_a)]; a zero count contributes 0
            // (excluded from `nonzero_pairs`).
            pairs
                .iter()
                .map(|&(a, k)| {
                    crate::genetics::lgamma(scratch.alpha[a as usize] + f64::from(k))
                        - scratch.lgamma_alpha[a as usize]
                })
                .sum()
        } else {
            pairs
                .iter()
                .map(|&(a, k)| f64::from(k) * scratch.log_p_effective[a as usize])
                .sum()
        };
        scratch.log_indep_per_g[g_idx] = ctx.shape.log_multinomial_coeffs[g_idx] + sum;
    }
}

/// Fill `scratch.log_prior_per_g` from the current `log_indep_per_g`
/// under the homogeneous-fixation invariant (every sample has the
/// same `f_s`). Caller must have just run [`fill_log_indep_per_g`].
///
/// Reads the cohort-wide `log_f` / `log_one_minus_f` from
/// `scratch.log_f_per_sample[0]` /
/// `scratch.log_one_minus_f_per_sample[0]`. The scalar
/// `log_sum_exp_2` short-circuits on `-∞` arguments (default config
/// path: `f_s = 0.0` → `log_f = -∞`), so the homozygous branch
/// usually folds to `log_one_minus_f + log_indep` without any
/// transcendental call. H4 in the perf review.
#[inline]
fn fill_log_prior_per_g_homogeneous<M: MathBackend>(
    ctx: EmContext<'_>,
    math: &M,
    scratch: &mut RecordScratch,
) {
    debug_assert!(
        !scratch.log_f_per_sample.is_empty(),
        "homogeneous_fixation set but scratch.log_f_per_sample is empty"
    );
    let log_f = scratch.log_f_per_sample[0];
    let log_one_minus_f = scratch.log_one_minus_f_per_sample[0];
    for g_idx in 0..ctx.n_genotypes {
        let log_indep = scratch.log_indep_per_g[g_idx];
        scratch.log_prior_per_g[g_idx] = match ctx.shape.homozygous_allele_for[g_idx] {
            Some(homo_allele) => log_sum_exp_2(
                math,
                log_one_minus_f + log_indep,
                log_f + scratch.log_p_effective[homo_allele as usize],
            ),
            None => log_one_minus_f + log_indep,
        };
    }
}

// ---------------------------------------------------------------------
// Numeric helpers
// ---------------------------------------------------------------------

#[inline]
fn safe_ln<M: MathBackend>(math: &M, x: f64) -> f64 {
    if x <= 0.0 {
        f64::NEG_INFINITY
    } else {
        math.ln(x)
    }
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

/// Lane-of-4 version of [`log_sum_exp_2`]. The both-`-INF` lane case
/// (which the scalar guards via early returns) is left to produce
/// `NaN` here — in the EM hot path, log-priors are sums of finite
/// pseudocount-stabilised log-frequencies, so this case never arises
/// on realistic inputs.
#[inline]
fn log_sum_exp_2_x4<M: MathBackend>(math: &M, a: wide::f64x4, b: wide::f64x4) -> wide::f64x4 {
    let m = a.fast_max(b);
    let ea = math.exp_x4(a - m);
    let eb = math.exp_x4(b - m);
    let ln_sum = math.ln_x4(ea + eb);
    m + ln_sum
}

/// Lane-of-4 version of [`log_sum_exp_slice`]. Each lane independently
/// computes `log_sum_exp` over the genotype dimension (the slice).
///
/// If any lane's max is `-∞`, that lane's result is forced to `-∞`
/// (matching the scalar's all-`-∞` short circuit) rather than the
/// `NaN` the raw computation would produce.
#[inline]
fn log_sum_exp_slice_x4<M: MathBackend>(math: &M, values: &[wide::f64x4]) -> wide::f64x4 {
    use wide::CmpEq;
    let mut m = wide::f64x4::splat(f64::NEG_INFINITY);
    for &v in values {
        m = m.fast_max(v);
    }
    let mut sum = wide::f64x4::splat(0.0);
    for &v in values {
        sum += math.exp_x4(v - m);
    }
    let ln_sum = math.ln_x4(sum);
    let candidate = m + ln_sum;
    // Per-lane fallback: if a lane's `m == -∞`, force result to `-∞`.
    let neg_inf = wide::f64x4::splat(f64::NEG_INFINITY);
    let all_neg_mask = m.simd_eq(neg_inf);
    all_neg_mask.blend(neg_inf, candidate)
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

    /// Diploid triallelic record (A=REF, T, C) with per-sample argmax
    /// genotypes + posterior rows, for the allele-prune tests.
    /// `genotype_order(2, 3)` is `AA, AB, BB, AC, BC, CC` (indices 0..5).
    fn triallelic_record(best: Vec<usize>, posteriors: Vec<f64>) -> PosteriorRecord {
        let n_samples = best.len();
        let n_genotypes = genotype_order(2, 3).len();
        assert_eq!(posteriors.len(), n_samples * n_genotypes);
        let scalars: Vec<AlleleSupportStats> = (0..n_samples)
            .flat_map(|_| {
                [
                    AlleleSupportStats {
                        num_obs: 20,
                        ..Default::default()
                    },
                    AlleleSupportStats {
                        num_obs: 5,
                        ..Default::default()
                    },
                    AlleleSupportStats {
                        num_obs: 7,
                        ..Default::default()
                    },
                ]
            })
            .collect();
        PosteriorRecord {
            locus: RecordLocus {
                chrom_id: 0,
                start: 100,
                end: 101,
            },
            alleles: simple_alleles(&[b"A", b"T", b"C"]),
            ploidy: 2,
            n_samples,
            n_genotypes,
            allele_frequencies: vec![0.6, 0.25, 0.15],
            compound_frequencies: vec![None, None, None],
            posteriors,
            best_genotype: best,
            gq_phred: vec![0.0; n_samples],
            qual_phred: 50.0,
            scalars,
            other_scalars: vec![],
            chain_anchor_flags: vec![false; n_samples * 3],
            diagnostics: EmDiagnostics {
                iterations: 1,
                final_max_delta_p: 0.0,
                converged: true,
            },
            paralog_posterior: None,
        }
    }

    #[test]
    fn prune_drops_alt_no_genotype_supports_and_reindexes() {
        // One sample, best genotype = idx 3 = [0,2] = ref/C, so ALT T
        // (allele 1) has zero genotype support and must be dropped.
        // Posteriors also carry mass on idx0 (0/0) and idx5 (2/2) so the
        // renormalised survivor row is checkable.
        let post = vec![0.05, 0.05, 0.05, 0.70, 0.10, 0.05];
        let mut rec = triallelic_record(vec![3], post);
        let removed = rec.prune_unsupported_alleles(99.0);

        assert_eq!(removed, 1, "exactly the unsupported ALT T is dropped");
        assert_eq!(rec.alleles.len(), 2);
        assert_eq!(rec.alleles[0].seq, b"A");
        assert_eq!(rec.alleles[1].seq, b"C", "T dropped, C compacted to ALT 1");
        assert_eq!(rec.n_genotypes, 3, "2 alleles diploid -> 3 genotypes");
        // old [0,2] (ref/C) re-indexes to new [0,1] = genotype index 1.
        assert_eq!(rec.best_genotype, vec![1]);
        assert!(rec.is_variant_call());

        // Surviving genotypes 00/0C/CC = old posteriors 0.05/0.70/0.05,
        // renormalised over their 0.80 sum.
        let row = rec.posteriors_row(0);
        assert_eq!(row.len(), 3);
        assert!((row[0] - 0.05 / 0.80).abs() < 1e-9);
        assert!((row[1] - 0.70 / 0.80).abs() < 1e-9);
        assert!((row[2] - 0.05 / 0.80).abs() < 1e-9);

        // AF: keep A=0.6, C=0.15, renormalise over 0.75 -> 0.8, 0.2.
        assert!((rec.allele_frequencies[0] - 0.8).abs() < 1e-9);
        assert!((rec.allele_frequencies[1] - 0.2).abs() < 1e-9);

        // scalars compacted to the A and C columns (20, 7); T's 5 dropped.
        assert_eq!(rec.scalars.len(), 2);
        assert_eq!(rec.scalars[0].num_obs, 20);
        assert_eq!(rec.scalars[1].num_obs, 7);

        // GQ recomputed from the renormalised best p = 0.70/0.80 = 0.875.
        let expected_gq = (PHRED_SCALE * (1.0 - 0.875_f64).log10()).clamp(0.0, 99.0);
        assert!((rec.gq_phred[0] - expected_gq).abs() < 1e-6);
    }

    #[test]
    fn prune_is_noop_when_every_alt_is_supported() {
        // Two samples: s0 = 0/1 (carries T), s1 = 0/2 (carries C), so
        // every allele is referenced -> nothing to prune.
        let post = vec![
            0.1, 0.7, 0.05, 0.05, 0.05, 0.05, // s0 best = idx1 (0/1)
            0.05, 0.05, 0.05, 0.7, 0.1, 0.05, // s1 best = idx3 (0/2)
        ];
        let mut rec = triallelic_record(vec![1, 3], post.clone());
        let removed = rec.prune_unsupported_alleles(99.0);

        assert_eq!(removed, 0);
        assert_eq!(rec.alleles.len(), 3);
        assert_eq!(rec.n_genotypes, 6);
        assert_eq!(rec.best_genotype, vec![1, 3]);
        assert_eq!(rec.posteriors, post, "posteriors untouched on no-op");
    }

    #[test]
    fn prune_short_circuits_for_biallelic() {
        // 2 alleles -> the < 3 guard returns immediately, untouched.
        let mut rec = PosteriorRecord {
            locus: RecordLocus {
                chrom_id: 0,
                start: 1,
                end: 2,
            },
            alleles: simple_alleles(&[b"A", b"T"]),
            ploidy: 2,
            n_samples: 1,
            n_genotypes: 3,
            allele_frequencies: vec![0.7, 0.3],
            compound_frequencies: vec![None, None],
            posteriors: vec![0.1, 0.8, 0.1],
            best_genotype: vec![1],
            gq_phred: vec![30.0],
            qual_phred: 50.0,
            scalars: vec![
                AlleleSupportStats {
                    num_obs: 10,
                    ..Default::default()
                },
                AlleleSupportStats {
                    num_obs: 10,
                    ..Default::default()
                },
            ],
            other_scalars: vec![],
            chain_anchor_flags: vec![false, false],
            diagnostics: EmDiagnostics {
                iterations: 1,
                final_max_delta_p: 0.0,
                converged: true,
            },
            paralog_posterior: None,
        };
        assert_eq!(rec.prune_unsupported_alleles(99.0), 0);
        assert_eq!(rec.alleles.len(), 2);
        assert_eq!(rec.best_genotype, vec![1]);
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
    // (`InterpUnivariateSimdMath`) can move without breaking them; the
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

    /// End-to-end wiring check for the SFS-marginal genotype prior: a single
    /// sample at a biallelic-diploid site whose likelihood favours hom-alt ~4:1
    /// (the "0 ref, 2 alt reads" shape). With the default HWE(p̂) + Dirichlet
    /// prior the n=1 reference pseudocount pulls the call to het (the bug); with
    /// the SFS prior's 2:1 het:hom-alt table the reads win and it calls hom-alt.
    #[test]
    fn sfs_prior_flips_single_sample_low_coverage_homalt_call() {
        use crate::var_calling::sfs_prior::{
            DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS, SfsGenotypePrior,
        };
        let eps = 0.01_f64;
        // genotype log-likelihoods [0/0, 0/1, 1/1] for "both reads alt".
        let likelihoods = vec![vec![
            (eps * eps).ln(),
            0.25_f64.ln(),
            ((1.0 - eps) * (1.0 - eps)).ln(),
        ]];
        let record = || merged_record_simple(1, 100, vec![b"A", b"G"], 2, likelihoods.clone());

        // Default HWE(p̂)+Dirichlet prior → het (the low-coverage over-call).
        assert_eq!(single_ok(record()).best_genotype, vec![1]);

        // SFS-marginal prior (het:hom-alt 2:1 at F=0) → hom-alt (correct).
        let table =
            SfsGenotypePrior::new(1e-3, DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS)
                .biallelic_diploid_log_priors(0.0);
        let cfg = PosteriorEngineConfig::with_project_defaults()
            .with_sfs_prior_tables(Some(vec![table]))
            .unwrap();
        let out = engine_for_with_config(record(), cfg);
        assert_eq!(out[0].as_ref().unwrap().best_genotype, vec![2]);
    }

    // ---- General Dirichlet-multinomial SFS prior (G4) ----

    /// Softmax of a log-prior row (test helper).
    fn softmax_row(logs: &[f64]) -> Vec<f64> {
        let m = logs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let exps: Vec<f64> = logs.iter().map(|&l| (l - m).exp()).collect();
        let z: f64 = exps.iter().sum();
        exps.iter().map(|&e| e / z).collect()
    }

    /// With **flat** likelihoods the per-sample posterior IS the normalised
    /// genotype prior, so at `F = 0` the engine's Dirichlet-multinomial path must
    /// reproduce the pure `genetics::dirichlet_multinomial_log_priors` value for
    /// the record's shape — the G4 "engine matches the pure primitive" check.
    fn assert_dm_prior_matches_pure(ploidy: u8, alleles: Vec<&[u8]>, theta: f64) {
        let n_alleles = alleles.len();
        let n_genotypes = genotype_order(ploidy, n_alleles).len();
        let flat = vec![vec![0.0_f64; n_genotypes]]; // one sample, flat likelihood
        let record = merged_record_simple(1, 100, alleles, ploidy, flat);
        let cfg = PosteriorEngineConfig::with_project_defaults()
            .with_nucleotide_diversity(Some(theta))
            .unwrap();
        let mut out = engine_for_with_config(record, cfg);
        let rec = out.remove(0).expect("posterior record");

        let alpha = crate::genetics::alpha_from_diversity(n_alleles, theta);
        let shape = shape_for(ploidy, n_alleles);
        let dm = crate::genetics::dirichlet_multinomial_log_priors(
            &shape.genotype_allele_counts,
            &shape.log_multinomial_coeffs,
            n_alleles,
            &alpha,
        );
        let expected = softmax_row(&dm);
        assert_eq!(rec.posteriors.len(), n_genotypes);
        for g in 0..n_genotypes {
            assert!(
                (rec.posteriors[g] - expected[g]).abs() < 1e-12,
                "genotype {g}: engine {} vs pure {}",
                rec.posteriors[g],
                expected[g]
            );
        }
    }

    #[test]
    fn dirichlet_prior_matches_pure_value_biallelic_diploid() {
        assert_dm_prior_matches_pure(2, vec![b"A", b"G"], 1e-3);
    }

    #[test]
    fn dirichlet_prior_matches_pure_value_triallelic_diploid() {
        assert_dm_prior_matches_pure(2, vec![b"A", b"G", b"C"], 5e-3);
    }

    #[test]
    fn dirichlet_prior_matches_pure_value_tetraploid_biallelic() {
        assert_dm_prior_matches_pure(4, vec![b"A", b"G"], 2e-2);
    }

    /// Full inbreeding (`F = 1`): the Wright mixture floors every heterozygote and
    /// leaves the homozygotes carrying the IBD marginal `α_a/Σα`, so for a
    /// biallelic site AA:BB = α_ref:α_alt = 1:θ. Confirms the DM prior threads
    /// through the existing Wright-`F` mixture unchanged at the `F = 1` limit.
    #[test]
    fn dirichlet_prior_full_inbreeding_concentrates_on_homozygotes() {
        let theta = 1e-3;
        let flat = vec![vec![0.0_f64; 3]];
        let record = merged_record_simple(1, 100, vec![b"A", b"G"], 2, flat);
        let cfg = PosteriorEngineConfig::with_project_defaults()
            .with_nucleotide_diversity(Some(theta))
            .unwrap()
            .with_fixation_index_default(1.0)
            .unwrap();
        let mut out = engine_for_with_config(record, cfg);
        let rec = out.remove(0).expect("posterior record");
        // Heterozygote floored.
        assert!(
            rec.posteriors[1] < 1e-9,
            "het not floored: {}",
            rec.posteriors[1]
        );
        // AA:BB = α_ref:α_alt = 1:θ.
        let ratio = rec.posteriors[0] / rec.posteriors[2];
        assert!(
            (ratio / (1.0 / theta) - 1.0).abs() < 1e-6,
            "AA:BB ratio {ratio} vs {}",
            1.0 / theta
        );
    }

    /// The Dirichlet-multinomial prior reproduces the biallelic low-coverage
    /// fix end-to-end: the same "0 ref, 2 alt reads" record the grid prior flips
    /// to hom-alt is also called hom-alt under the general θ-driven path (the DM
    /// het:hom-alt ratio ≈ 2:1 overrides the weak likelihood).
    #[test]
    fn dirichlet_prior_flips_single_sample_low_coverage_homalt_call() {
        let eps = 0.01_f64;
        let likelihoods = vec![vec![
            (eps * eps).ln(),
            0.25_f64.ln(),
            ((1.0 - eps) * (1.0 - eps)).ln(),
        ]];
        let record = merged_record_simple(1, 100, vec![b"A", b"G"], 2, likelihoods);
        // Plug-in HWE prior (no θ) → het (the low-coverage over-call).
        assert_eq!(single_ok(record.clone()).best_genotype, vec![1]);
        // General DM SFS prior (θ set) → hom-alt (correct).
        let cfg = PosteriorEngineConfig::with_project_defaults()
            .with_nucleotide_diversity(Some(1e-3))
            .unwrap();
        let mut out = engine_for_with_config(record, cfg);
        assert_eq!(out.remove(0).unwrap().best_genotype, vec![2]);
    }

    /// The production **SIMD** backend (`InterpUnivariateSimdMath`, `HAS_LANE_4`)
    /// runs its own DM fill loop and 4-lane batch consumption, distinct from the
    /// scalar `e_step` the other DM tests exercise. With 5 samples the batch loop
    /// (samples 0..4) AND the scalar tail (sample 4) both run; the SIMD posteriors
    /// must match the exact scalar path (pinned to the pure DM value elsewhere)
    /// within the interpolation tolerance. Guards a silent typo in the SIMD DM
    /// branch — the path production actually takes.
    #[test]
    fn dirichlet_prior_simd_backend_matches_scalar_across_batch_and_tail() {
        let theta = 5e-3;
        let flat = vec![vec![0.0_f64; 3]; 5]; // 5 samples → batch of 4 + 1 tail
        let record = merged_record_simple(1, 100, vec![b"A", b"G"], 2, flat);
        let cfg = || {
            PosteriorEngineConfig::with_project_defaults()
                .with_nucleotide_diversity(Some(theta))
                .unwrap()
        };
        // Exact scalar reference.
        let mut scalar = engine_for_with_config(record.clone(), cfg());
        let scalar = scalar.remove(0).expect("scalar posterior");
        // Production SIMD backend.
        let simd_upstream = std::iter::once(Ok::<_, PerGroupMergerError>(record));
        let mut simd: Vec<_> = PosteriorEngine::with_math_backend(
            simd_upstream,
            cfg(),
            super::backends::InterpUnivariateSimdMath,
        )
        .collect();
        let simd = simd.remove(0).expect("simd posterior");
        assert_eq!(scalar.posteriors.len(), simd.posteriors.len());
        for (a, b) in scalar.posteriors.iter().zip(&simd.posteriors) {
            assert!(
                (a - b).abs() < 1e-4,
                "SIMD DM diverged from scalar: {a} vs {b}"
            );
        }
        assert_eq!(scalar.best_genotype, simd.best_genotype);
    }

    /// `with_nucleotide_diversity` accepts finite non-negative θ̂ (including 0)
    /// and rejects NaN / +∞ / negative — the wiring-error guard.
    #[test]
    fn with_nucleotide_diversity_validates_theta() {
        let base = PosteriorEngineConfig::with_project_defaults;
        assert!(base().with_nucleotide_diversity(Some(0.0)).is_ok());
        assert!(base().with_nucleotide_diversity(Some(1e-3)).is_ok());
        assert!(base().with_nucleotide_diversity(None).is_ok());
        assert_eq!(
            base().with_nucleotide_diversity(Some(-1e-3)).unwrap_err(),
            PosteriorEngineConfigError::InvalidNucleotideDiversity { got: -1e-3 }
        );
        assert!(matches!(
            base()
                .with_nucleotide_diversity(Some(f64::NAN))
                .unwrap_err(),
            PosteriorEngineConfigError::InvalidNucleotideDiversity { .. }
        ));
        assert!(matches!(
            base()
                .with_nucleotide_diversity(Some(f64::INFINITY))
                .unwrap_err(),
            PosteriorEngineConfigError::InvalidNucleotideDiversity { .. }
        ));
    }

    /// Under the DM prior, per-sample inbreeding F still routes through the
    /// heterogeneous Wright mixture: a selfing sample (F=0.9) carries a smaller
    /// heterozygote prior than an outbred one (F=0) at the same site. Exercises
    /// the per-sample (non-homogeneous) branch with the α-based IBD term.
    #[test]
    fn dirichlet_prior_heterogeneous_inbreeding_lowers_het_for_the_selfer() {
        let flat = vec![vec![0.0_f64; 3], vec![0.0_f64; 3]]; // 2 samples
        let record = merged_record_simple(1, 100, vec![b"A", b"G"], 2, flat);
        let mut cfg = PosteriorEngineConfig::with_project_defaults()
            .with_nucleotide_diversity(Some(1e-2))
            .unwrap();
        cfg.fixation_index_overrides = Some(vec![0.0, 0.9]);
        let mut out = engine_for_with_config(record, cfg);
        let rec = out.remove(0).expect("posterior record");
        // posteriors are sample-major: [s0: 0/0,0/1,1/1 | s1: 0/0,0/1,1/1].
        let het_outbred = rec.posteriors[1];
        let het_selfer = rec.posteriors[3 + 1];
        assert!(
            het_selfer < het_outbred,
            "selfer het {het_selfer} not < outbred het {het_outbred}"
        );
    }

    /// θ̂ = 0 (a fully-invariant cohort, or `--diversity 0`): α_alt floors to
    /// `MIN_ALT_CONCENTRATION` so the prior is effectively certain hom-ref rather
    /// than `NaN` — graceful degradation matching the biallelic grid path.
    #[test]
    fn dirichlet_prior_zero_theta_is_certain_homref() {
        let flat = vec![vec![0.0_f64; 3]];
        let record = merged_record_simple(1, 100, vec![b"A", b"G"], 2, flat);
        let cfg = PosteriorEngineConfig::with_project_defaults()
            .with_nucleotide_diversity(Some(0.0))
            .unwrap();
        let mut out = engine_for_with_config(record, cfg);
        let rec = out.remove(0).expect("posterior record");
        assert!(
            rec.posteriors[0] > 0.999_999_999,
            "homref {} not ~1",
            rec.posteriors[0]
        );
        assert_eq!(rec.best_genotype, vec![0]);
    }

    /// The SFS prior is bypassed for non-biallelic-diploid records even when
    /// configured: a triallelic record still genotypes via the HWE path (the
    /// table only covers 3 genotypes, so a cohort-size/shape mismatch must not
    /// occur — the engine routes multiallelic to HWE regardless).
    #[test]
    fn sfs_prior_ignored_for_multiallelic_record() {
        use crate::var_calling::sfs_prior::{
            DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS, SfsGenotypePrior,
        };
        // Triallelic diploid, one sample: 6 genotypes.
        let likelihoods = vec![vec![-10.0, -1.0, -10.0, -10.0, -10.0, -10.0]];
        let record = merged_record_simple(1, 100, vec![b"A", b"G", b"T"], 2, likelihoods);
        let table =
            SfsGenotypePrior::new(1e-3, DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS)
                .biallelic_diploid_log_priors(0.0);
        let cfg = PosteriorEngineConfig::with_project_defaults()
            .with_sfs_prior_tables(Some(vec![table]))
            .unwrap();
        // Must not panic or error on the shape mismatch — routed to the HWE path.
        let out = engine_for_with_config(record, cfg);
        assert!(out[0].is_ok());
    }

    /// Proptest-friendly variant of [`single_ok`]: returns the engine's
    /// `Result` so the caller can drop adversarial inputs whose EM
    /// hit the iteration cap (`pr.diagnostics.converged == false` —
    /// invariants like "permutation invariance" are only meaningful
    /// when both runs reached the same fixed point). The engine no
    /// longer surfaces non-convergence as `Err(DidNotConverge)`; it
    /// emits the record with `converged: false` instead.
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
    fn is_variant_call_false_when_every_sample_is_hom_ref() {
        // Two samples, both with strong-RR evidence — driver should
        // skip this record before it reaches the VCF writer.
        let record = merged_record_simple(
            1,
            100,
            vec![b"A", b"C"],
            2,
            vec![vec![0.0, -50.0, -50.0], vec![0.0, -50.0, -50.0]],
        );
        let pr = single_ok(record);
        assert_eq!(pr.best_genotype, vec![0, 0]);
        assert!(!pr.is_variant_call());
    }

    #[test]
    fn is_variant_call_true_when_any_sample_is_non_ref() {
        // One sample RR, one sample AA — exactly the case that today
        // emits AC>0 but historically motivated the filter on the
        // other side: a single non-ref carrier is enough.
        let record = merged_record_simple(
            1,
            100,
            vec![b"A", b"C"],
            2,
            vec![vec![0.0, -50.0, -50.0], vec![-50.0, -50.0, 0.0]],
        );
        let pr = single_ok(record);
        assert!(pr.is_variant_call());
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
    fn hard_iteration_cap_emits_record_with_converged_false() {
        let record = merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![0.0, -1.0, -2.0]]);
        let config = PosteriorEngineConfig {
            max_iterations: 1,
            convergence_threshold: 1e-12,
            ..Default::default()
        };
        let out = engine_for_with_config(record, config);
        let pr = match out.into_iter().next() {
            Some(Ok(pr)) => pr,
            other => panic!("expected Ok(PosteriorRecord), got {other:?}"),
        };
        assert_eq!(pr.locus.chrom_id, 1);
        assert_eq!(pr.locus.start, 100);
        assert_eq!(pr.locus.end, 100);
        assert!(
            !pr.diagnostics.converged,
            "iteration-cap exit must set converged=false"
        );
        assert_eq!(pr.diagnostics.iterations, 1);
        assert!(
            pr.diagnostics.final_max_delta_p.is_finite() && pr.diagnostics.final_max_delta_p >= 0.0,
            "final_max_delta_p = {}",
            pr.diagnostics.final_max_delta_p,
        );
        // Even though EM didn't converge, the emitted record must
        // still satisfy the basic invariants — the final E-step
        // normalises posteriors regardless of convergence.
        for s in 0..pr.n_samples {
            let sum: f64 = pr.posteriors_row(s).iter().sum();
            assert!((sum - 1.0).abs() < 1e-9, "row sum = {sum}");
        }
        let p_sum: f64 = pr.allele_frequencies.iter().sum();
        assert!((p_sum - 1.0).abs() < 1e-9, "p̂ sum = {p_sum}");
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
    fn site_qual_matches_closed_form_for_single_sample() {
        // Single diploid biallelic sample with arbitrary log-likelihoods
        // — the AC convolution collapses to the per-sample row, so the
        // expected QUAL is the Beta-Binomial(2, α_alt, α_ref) posterior
        // on `K = 0` computed in closed form. Catches sign / base / prior
        // mistakes in `compute_qual_via_exact_af` against the math.
        let l = vec![-1.5_f64, -0.2, -2.0];
        let record = merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![l.clone()]);
        let pr = single_ok(record);

        let alpha_alt = DEFAULT_SNP_ALT_PSEUDOCOUNT;
        let alpha_ref = DEFAULT_REF_PSEUDOCOUNT;
        // log B(α_alt, α_ref).
        let log_beta_norm =
            ln_gamma(alpha_alt) + ln_gamma(alpha_ref) - ln_gamma(alpha_alt + alpha_ref);
        // P(K = k | 2, α_alt, α_ref) for k = 0, 1, 2.
        let log_prior = |k: u32| -> f64 {
            let kf = k as f64;
            let log_binom = ln_gamma(3.0) - ln_gamma(kf + 1.0) - ln_gamma(2.0 - kf + 1.0);
            let log_beta_kk = ln_gamma(alpha_alt + kf) + ln_gamma(alpha_ref + 2.0 - kf);
            log_binom + log_beta_kk - ln_gamma(alpha_alt + alpha_ref + 2.0) - log_beta_norm
        };
        let log_post: [f64; 3] = [
            l[0] + log_prior(0),
            l[1] + log_prior(1),
            l[2] + log_prior(2),
        ];
        let m = log_post.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let log_z = m + (log_post.iter().map(|x| (x - m).exp()).sum::<f64>()).ln();
        let expected_qual = -10.0 * (log_post[0] - log_z) / std::f64::consts::LN_10;

        assert!(
            approx(pr.qual_phred, expected_qual, 1e-6),
            "qual_phred={} expected={}",
            pr.qual_phred,
            expected_qual,
        );
    }

    #[test]
    fn site_qual_does_not_inflate_when_no_sample_has_strong_alt_evidence() {
        // Pathology the new formula fixes: a cohort where every
        // sample has weak ambiguous evidence (mild hint toward het)
        // but no sample is a confident carrier. Under the old
        // `Π_s P(GT_s = hom-ref)` formula every such sample
        // contributed a few extra phred and QUAL grew linearly with
        // N — that's how the synthetic-10-duplicate cohort hit
        // QUAL ≈ 30 with no real signal. Under the cohort-AC
        // posterior with a sparse-het prior, weak per-sample
        // evidence translates to a weak posterior shift at the
        // cohort-AC level, bounded by the prior; QUAL must stay
        // small regardless of N.
        let weak = vec![0.0_f64, -1.0, -3.0]; // mild ref preference per sample

        let one = merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![weak.clone()]);
        let qual_one = single_ok(one).qual_phred;

        let many_likelihoods: Vec<Vec<f64>> = (0..20).map(|_| weak.clone()).collect();
        let many = merged_record_simple(1, 100, vec![b"A", b"C"], 2, many_likelihoods);
        let qual_many = single_ok(many).qual_phred;

        // The defaults' sparse-het prior keeps QUAL ≲ 5 phred for any
        // cohort under this evidence shape. Under the old formula
        // QUAL would have grown by ~10 phred between N=1 and N=20.
        assert!(
            qual_one < 5.0,
            "qual_one={qual_one} should be < 5 (weak evidence, no clear carrier)",
        );
        assert!(
            qual_many < 5.0,
            "qual_many={qual_many} should be < 5 (cohort growth must not inflate QUAL when no sample is a clear carrier)",
        );
    }

    #[test]
    fn gq_phred_matches_closed_form() {
        let record =
            merged_record_simple(1, 100, vec![b"A", b"C"], 2, vec![vec![0.0, -50.0, -50.0]]);
        let pr = single_ok(record);
        let p_best = pr.posteriors_row(0)[pr.best_genotype[0]];
        let p_best_clamped = p_best.min(1.0 - f64::EPSILON);
        let expected_gq =
            ((-10.0_f64) * (1.0 - p_best_clamped).log10()).clamp(0.0, DEFAULT_MAX_GQ_PHRED);
        assert!(
            approx(pr.gq_phred[0], expected_gq, 1e-9),
            "gq_phred={} expected={}",
            pr.gq_phred[0],
            expected_gq
        );
    }

    #[test]
    fn site_qual_is_large_but_finite_under_extreme_alt_evidence() {
        // Under the new `−10·log10 P(K = 0 | data)` formula QUAL is
        // bounded by the convolution arithmetic — it grows large but
        // does not reach `f64::INFINITY` unless the entire numerator
        // at `K = 0` underflows to `−∞` exactly. The VCF writer caps
        // at `QUAL_MAX = 9999.0` anyway. This test pins the
        // "large-finite, not ±∞ or NaN" contract for an extreme
        // single-sample alt-favoured input.
        let record = merged_record_simple(
            1,
            100,
            vec![b"A", b"C"],
            2,
            vec![vec![-1000.0, -1000.0, 0.0]],
        );
        let pr = single_ok(record);
        assert!(
            pr.qual_phred.is_finite() && pr.qual_phred > 100.0,
            "qual_phred = {} (expected large-finite)",
            pr.qual_phred,
        );
        assert!(pr.gq_phred[0].is_finite());
        assert!(pr.gq_phred[0] <= DEFAULT_MAX_GQ_PHRED + 1e-9);
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
            AlleleSupportStats::new(num_obs_ref, q_sum_ref, num_obs_ref / 2, 0, 0, 0, 0),
            AlleleSupportStats::new(num_obs_alt, q_sum_alt, num_obs_alt / 2, 0, 0, 0, 0),
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
            AlleleSupportStats::new(49, 0.001_f64.ln() * 49.0, 24, 0, 0, 0, 0),
            AlleleSupportStats::new(1, 0.001_f64.ln(), 0, 0, 0, 0, 0),
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

    use super::backends::{ExactMath, InterpUnivariateMath, InterpUnivariateSimdMath, MathBackend};

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
            samples.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
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
            assert_eq!(
                exact.allele_frequencies.len(),
                approx.allele_frequencies.len()
            );

            for (e, a) in exact.posteriors.iter().zip(approx.posteriors.iter()) {
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
                sorted.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
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
            merged_record_simple(1, 400, vec![b"A", b"C", b"G"], 2, trial_likelihoods),
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
        config.fixation_index_overrides = Some(
            (0..32)
                .map(|s: usize| if s.is_multiple_of(2) { 0.0 } else { 0.4 })
                .collect(),
        );
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
            .map(|s: usize| {
                if s.is_multiple_of(3) {
                    None
                } else {
                    Some(0.05_f64)
                }
            })
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

    #[test]
    fn interp_univariate_simd_clears_approximate_parity_budget() {
        let report = build_accuracy_battery(InterpUnivariateSimdMath);
        eprintln!("InterpUnivariateSimdMath accuracy report:\n{report:#?}");

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
            // Records that hit the iteration cap on adversarial
            // likelihood matrices are emitted with `converged: false`
            // (the engine no longer surfaces non-convergence as an
            // `Err`). Skip those — the per-row sum invariant holds
            // for them too, but separating the converged-only
            // population keeps this proptest a clean check on the
            // converged code path.
            let pr = match try_single(record) {
                Ok(pr) if pr.diagnostics.converged => pr,
                Ok(_) => return Ok(()),
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
                Ok(pr) if pr.diagnostics.converged => pr,
                Ok(_) => return Ok(()),
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
            // to be meaningful — un-converged runs stop at the
            // iteration cap with parameter values that depend on
            // path-history. Skip when either run lands at
            // `converged: false`.
            let pr_a = match try_single(record_a) {
                Ok(pr) if pr.diagnostics.converged => pr,
                Ok(_) => return Ok(()),
                Err(other) => return Err(TestCaseError::fail(format!(
                    "unexpected engine error on record_a: {other}"
                ))),
            };
            let pr_b = match try_single(record_b) {
                Ok(pr) if pr.diagnostics.converged => pr,
                Ok(_) => return Ok(()),
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

    // ---------- PosteriorEngineConfig builder validation tests -----

    #[test]
    fn config_new_returns_project_defaults() {
        let cfg = PosteriorEngineConfig::new();
        assert_eq!(cfg.convergence_threshold, DEFAULT_CONVERGENCE_THRESHOLD);
        assert_eq!(cfg.max_iterations, DEFAULT_MAX_ITERATIONS);
        assert_eq!(cfg.ref_pseudocount, DEFAULT_REF_PSEUDOCOUNT);
        assert_eq!(cfg.snp_alt_pseudocount, DEFAULT_SNP_ALT_PSEUDOCOUNT);
        assert_eq!(cfg.indel_alt_pseudocount, DEFAULT_INDEL_ALT_PSEUDOCOUNT);
        assert_eq!(
            cfg.compound_alt_pseudocount,
            DEFAULT_COMPOUND_ALT_PSEUDOCOUNT
        );
        assert_eq!(cfg.fixation_index_default, DEFAULT_INBREEDING_COEFFICIENT);
        assert!(cfg.fixation_index_overrides.is_none());
        assert_eq!(cfg.max_gq_phred, DEFAULT_MAX_GQ_PHRED);
        assert!(!cfg.approximate_posterior_calculation);
        assert!(cfg.contamination.is_none());
    }

    #[test]
    fn with_convergence_threshold_rejects_nan() {
        let err = PosteriorEngineConfig::new()
            .with_convergence_threshold(f64::NAN)
            .unwrap_err();
        assert!(matches!(
            err,
            PosteriorEngineConfigError::InvalidConvergenceThreshold { .. }
        ));
    }

    #[test]
    fn with_convergence_threshold_rejects_zero() {
        let err = PosteriorEngineConfig::new()
            .with_convergence_threshold(0.0)
            .unwrap_err();
        assert!(matches!(
            err,
            PosteriorEngineConfigError::InvalidConvergenceThreshold { .. }
        ));
    }

    #[test]
    fn with_convergence_threshold_rejects_just_above_max() {
        let err = PosteriorEngineConfig::new()
            .with_convergence_threshold(CONVERGENCE_THRESHOLD_RANGE_MAX + 1e-9)
            .unwrap_err();
        assert!(matches!(
            err,
            PosteriorEngineConfigError::InvalidConvergenceThreshold { .. }
        ));
    }

    #[test]
    fn with_convergence_threshold_accepts_max() {
        PosteriorEngineConfig::new()
            .with_convergence_threshold(CONVERGENCE_THRESHOLD_RANGE_MAX)
            .unwrap();
    }

    #[test]
    fn with_max_iterations_rejects_zero() {
        let err = PosteriorEngineConfig::new()
            .with_max_iterations(0)
            .unwrap_err();
        assert_eq!(
            err,
            PosteriorEngineConfigError::InvalidMaxIterations { got: 0 }
        );
    }

    #[test]
    fn with_max_iterations_rejects_above_max() {
        let err = PosteriorEngineConfig::new()
            .with_max_iterations(MAX_ITERATIONS_RANGE_MAX + 1)
            .unwrap_err();
        assert!(matches!(
            err,
            PosteriorEngineConfigError::InvalidMaxIterations { .. }
        ));
    }

    #[test]
    fn with_ref_pseudocount_rejects_zero() {
        let err = PosteriorEngineConfig::new()
            .with_ref_pseudocount(0.0)
            .unwrap_err();
        assert!(matches!(
            err,
            PosteriorEngineConfigError::InvalidPseudocount {
                field: "ref_pseudocount",
                ..
            }
        ));
    }

    #[test]
    fn with_ref_pseudocount_rejects_above_max() {
        let err = PosteriorEngineConfig::new()
            .with_ref_pseudocount(PSEUDOCOUNT_RANGE_MAX + 1.0)
            .unwrap_err();
        assert!(matches!(
            err,
            PosteriorEngineConfigError::InvalidPseudocount {
                field: "ref_pseudocount",
                ..
            }
        ));
    }

    #[test]
    fn with_indel_alt_pseudocount_rejects_zero_names_its_own_field() {
        // Mi2 swap-risk regression: every pseudocount setter must
        // surface the field name in the error, so an
        // out-of-range value in the chained call is unambiguous
        // about which knob the caller meant.
        let err = PosteriorEngineConfig::new()
            .with_indel_alt_pseudocount(0.0)
            .unwrap_err();
        assert!(matches!(
            err,
            PosteriorEngineConfigError::InvalidPseudocount {
                field: "indel_alt_pseudocount",
                ..
            }
        ));
    }

    #[test]
    fn with_fixation_index_default_rejects_above_one() {
        let err = PosteriorEngineConfig::new()
            .with_fixation_index_default(1.5)
            .unwrap_err();
        assert!(matches!(
            err,
            PosteriorEngineConfigError::InvalidFixationIndex { .. }
        ));
    }

    #[test]
    fn with_fixation_index_default_accepts_boundaries() {
        PosteriorEngineConfig::new()
            .with_fixation_index_default(0.0)
            .unwrap();
        PosteriorEngineConfig::new()
            .with_fixation_index_default(1.0)
            .unwrap();
    }

    #[test]
    fn with_max_gq_phred_rejects_at_min_exclusive() {
        let err = PosteriorEngineConfig::new()
            .with_max_gq_phred(GQ_PHRED_RANGE_MIN_EXCLUSIVE)
            .unwrap_err();
        assert!(matches!(
            err,
            PosteriorEngineConfigError::InvalidMaxGqPhred { .. }
        ));
    }

    #[test]
    fn with_max_gq_phred_rejects_above_max() {
        let err = PosteriorEngineConfig::new()
            .with_max_gq_phred(GQ_PHRED_RANGE_MAX + 0.1)
            .unwrap_err();
        assert!(matches!(
            err,
            PosteriorEngineConfigError::InvalidMaxGqPhred { .. }
        ));
    }

    #[test]
    fn with_max_gq_phred_accepts_max() {
        PosteriorEngineConfig::new()
            .with_max_gq_phred(GQ_PHRED_RANGE_MAX)
            .unwrap();
    }

    #[test]
    fn with_contamination_round_trips() {
        // M4: the only public path to set `contamination` is the
        // validating setter, which today is the identity but
        // reserves the slot for future cross-field checks.
        let cfg = PosteriorEngineConfig::new()
            .with_contamination(None)
            .unwrap();
        assert!(cfg.contamination.is_none());
    }
}
