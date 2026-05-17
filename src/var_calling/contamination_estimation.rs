//! Contamination-estimation side-pass.
//!
//! Estimates per-sample contamination fractions `c_s` and per-batch
//! contamination-source allele-class distributions `q_b` from raw
//! `.psp` data using block-wise online EM on a read-level mixture
//! model. The output is consumed by the Stage 6 posterior engine as
//! frozen inputs.
//!
//! Spec: `doc/devel/specs/contamination_estimation.md`. The plan that
//! drives this implementation is at
//! `doc/devel/implementation_plans/contamination_estimation.md`.
//!
//! Architecture: a side-pass, not a stage. It consumes
//! [`PerPositionPileups`] from a [`PerPositionMerger`] (or any
//! iterator with that shape) but does **not** go through Stage 3
//! (DUST filter), Stage 4 (grouping), or Stage 5 (per-group merger).
//! Bypassing them is the point: contamination is a per-read mixture
//! parameter and the raw `.psp` aggregates carry everything the
//! estimator needs.
//!
//! The estimator stops by exactly one of two stopping criteria,
//! selected at config-construction time via the [`StoppingMode`]
//! enum:
//!
//! - [`StoppingMode::Convergence`] — run until the per-block snapshot
//!   of `c_s` is stable across `stability_blocks` consecutive blocks
//!   within `tolerance`. Default.
//! - [`StoppingMode::FixedSites`] — process exactly `num_sites`
//!   informative positions and return whatever estimate the running
//!   sufficient statistics yield.
//!
//! The two modes are mutually exclusive at the type level (the enum
//! makes both-at-once unrepresentable). The future CLI parser is
//! responsible for mapping the two `--contamination-stability-tolerance`
//! / `--contamination-num-sites` flags into the enum, and for erroring
//! at parse time if both are passed.

use std::fmt;

use thiserror::Error;

use crate::per_sample_pileup::pileup::PileupRecord;
use crate::var_calling::per_position_merger::{PerPositionMergerError, PerPositionPileups};

// ---------------------------------------------------------------------
// Default constants
// ---------------------------------------------------------------------

/// Default convergence-mode tolerance on `max_s |c_s_new − c_s_prev|`
/// between consecutive block snapshots. Spec value.
pub const DEFAULT_STABILITY_TOLERANCE: f64 = 1e-3;

/// Default number of consecutive within-tolerance snapshots required
/// to declare convergence. Spec value.
pub const DEFAULT_STABILITY_BLOCKS: u32 = 3;

/// Default informative-site count per block heartbeat. The same knob
/// controls how often the online-EM parameters refresh and (in
/// convergence mode) how often the stability check fires. Spec value.
pub const DEFAULT_BLOCK_SIZE: u32 = 1000;

/// Default minimum read depth for a (sample, site) pair to be considered
/// for the per-sample hom-major call. Adopted from VerifyBamID. Spec value.
pub const DEFAULT_MIN_DEPTH: u32 = 10;

/// Default minimum observed major-allele fraction for a (sample, site)
/// pair to qualify as confidently homozygous-major. Adopted from
/// VerifyBamID. Spec value.
pub const DEFAULT_MIN_MAJOR_FRACTION: f64 = 0.95;

/// Default cohort-summed minimum minor-allele read count for a site
/// to be informative (Step 1a). Spec value.
pub const DEFAULT_MIN_COHORT_MINOR_COUNT: u32 = 2;

/// Default cohort-summed minimum minor-allele fraction for a site to
/// be informative (Step 1a). Spec value.
pub const DEFAULT_MIN_COHORT_MINOR_FRACTION: f64 = 0.005;

/// Default minimum batch size for contamination correction. Batches
/// below this floor (and all singleton batches) get `c_s = 0` regardless
/// of the EM output, per the spec.
pub const DEFAULT_MIN_BATCH_SIZE_FOR_CONTAMINATION: u32 = 5;

/// Initial per-sample `c_s` value used before the first block boundary
/// refresh. Set to the typical sequencing-batch contamination prior
/// per the spec rather than 0 so the first block's `γ_r` values are
/// computed against a plausible (rather than degenerate) parameter
/// snapshot.
pub const DEFAULT_CS_INIT: f64 = 0.02;

/// Default Dirichlet pseudocount on REF reads for the `q_b` update.
/// Mirrors the posterior engine's REF pseudocount; the side-pass
/// shares the value rather than re-introducing it as a separate knob
/// per the spec's "Dirichlet pseudocounts `α_a` from the per-position
/// pseudocounts" line.
pub const DEFAULT_REF_PSEUDOCOUNT: f64 = 10.0;

/// Default Dirichlet pseudocount on SNP-alt reads for the `q_b`
/// update. Mirrors the posterior engine value.
pub const DEFAULT_SNP_ALT_PSEUDOCOUNT: f64 = 0.01;

/// Default Dirichlet pseudocount on indel-alt reads for the `q_b`
/// update. Mirrors the posterior engine value.
pub const DEFAULT_INDEL_ALT_PSEUDOCOUNT: f64 = 0.00125;

// ---------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------

/// Allele class used to index per-batch `q_b`. The side-pass operates
/// on per-position observations only, so compounds (a Stage 5 concept)
/// never appear here — three classes is the full alphabet.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum AlleleClass {
    Ref = 0,
    SnpAlt = 1,
    IndelAlt = 2,
}

impl AlleleClass {
    /// Classify a per-position allele observation by comparing its
    /// bytes to the per-record REF observation's bytes. Equal bytes
    /// → REF; same length but different bytes → SNP/MNP; different
    /// length → indel (insertions and deletions both fold into
    /// `IndelAlt`).
    ///
    /// The walker invariant says `record.alleles[0]` is unconditionally
    /// REF, so the caller passes `record.alleles[0].seq` as `ref_seq`.
    pub fn classify(observation_seq: &[u8], ref_seq: &[u8]) -> Self {
        if observation_seq == ref_seq {
            AlleleClass::Ref
        } else if observation_seq.len() == ref_seq.len() {
            AlleleClass::SnpAlt
        } else {
            AlleleClass::IndelAlt
        }
    }
}

/// Number of slots in a `q_b` vector. Equals the number of variants
/// of [`AlleleClass`].
pub const N_ALLELE_CLASSES: usize = 3;

/// Stopping criterion for the side-pass. The two modes are mutually
/// exclusive by construction — the enum makes "both at once"
/// unrepresentable.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum StoppingMode {
    /// Convergence mode. Process informative sites until the per-block
    /// snapshot of `c_s` is stable across `stability_blocks`
    /// consecutive blocks within `tolerance`. If the upstream stream
    /// is exhausted first, the side-pass returns
    /// [`ContaminationEstimationError::DidNotConverge`].
    Convergence {
        tolerance: f64,
        stability_blocks: u32,
    },
    /// Fixed-N mode. Process exactly `num_sites` informative positions
    /// and return whatever estimate the running sufficient statistics
    /// yield. No convergence check. If the upstream stream is
    /// exhausted before `num_sites` are found, the side-pass returns
    /// [`ContaminationEstimationError::InsufficientSites`].
    FixedSites { num_sites: u32 },
}

impl Default for StoppingMode {
    fn default() -> Self {
        StoppingMode::Convergence {
            tolerance: DEFAULT_STABILITY_TOLERANCE,
            stability_blocks: DEFAULT_STABILITY_BLOCKS,
        }
    }
}

/// Tunable knobs for [`estimate_contamination`].
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct ContaminationEstimationConfig {
    /// Stopping criterion. Defaults to
    /// [`StoppingMode::Convergence`] at the project defaults.
    pub stopping_mode: StoppingMode,
    /// Sites per heartbeat — controls both the online-EM parameter
    /// refresh and (in convergence mode) the stability snapshot.
    pub block_size: u32,
    /// Per-sample minimum read depth at a site to consider the
    /// (sample, site) pair informative.
    pub min_depth: u32,
    /// Per-sample minimum observed major-allele fraction to call
    /// hom-major at a site.
    pub min_major_fraction: f64,
    /// Cohort-summed minimum minor-allele read count for a site to be
    /// informative (Step 1a).
    pub min_cohort_minor_count: u32,
    /// Cohort-summed minimum minor-allele fraction for a site to be
    /// informative (Step 1a).
    pub min_cohort_minor_fraction: f64,
    /// Batches strictly smaller than this floor get `c_s = 0` for
    /// every sample regardless of the EM. Singleton batches (size
    /// exactly 1) are always treated as below-floor regardless of
    /// this value because their `c_s` and `q_b` are unidentifiable
    /// from each other.
    pub min_batch_size_for_contamination: u32,
    /// Dirichlet pseudocount on REF reads in the `q_b` update. Mirrors
    /// the posterior engine's REF pseudocount by default.
    pub ref_pseudocount: f64,
    /// Dirichlet pseudocount on SNP-alt reads in the `q_b` update.
    pub snp_alt_pseudocount: f64,
    /// Dirichlet pseudocount on indel-alt reads in the `q_b` update.
    pub indel_alt_pseudocount: f64,
}

impl Default for ContaminationEstimationConfig {
    fn default() -> Self {
        Self::with_project_defaults()
    }
}

impl ContaminationEstimationConfig {
    /// Construct the side-pass config with the project defaults
    /// listed against each field on the struct docs.
    pub fn with_project_defaults() -> Self {
        Self {
            stopping_mode: StoppingMode::default(),
            block_size: DEFAULT_BLOCK_SIZE,
            min_depth: DEFAULT_MIN_DEPTH,
            min_major_fraction: DEFAULT_MIN_MAJOR_FRACTION,
            min_cohort_minor_count: DEFAULT_MIN_COHORT_MINOR_COUNT,
            min_cohort_minor_fraction: DEFAULT_MIN_COHORT_MINOR_FRACTION,
            min_batch_size_for_contamination: DEFAULT_MIN_BATCH_SIZE_FOR_CONTAMINATION,
            ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT,
            snp_alt_pseudocount: DEFAULT_SNP_ALT_PSEUDOCOUNT,
            indel_alt_pseudocount: DEFAULT_INDEL_ALT_PSEUDOCOUNT,
        }
    }
}

/// How a [`ContaminationEstimates`] was produced. Carried for
/// provenance and for the optional TSV artefact.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ContaminationEstimateSource {
    /// Produced by [`estimate_contamination`] with the given mode and
    /// site count.
    SidePass {
        mode: StoppingMode,
        sites_processed: u32,
    },
    /// Loaded from user-supplied files (`--contamination-estimates`
    /// and/or `--contamination-source-distributions`).
    UserSupplied,
    /// A combination of side-pass output and user overrides.
    Mixed,
    /// All zeros — the engine was built without
    /// `--contamination-batches`, so contamination correction is
    /// disabled.
    Disabled,
}

/// Output of the side-pass. Consumed by Stage 6 as a frozen input via
/// [`PosteriorEngineConfig::contamination`].
///
/// [`PosteriorEngineConfig::contamination`]: crate::var_calling::posterior_engine::PosteriorEngineConfig::contamination
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ContaminationEstimates {
    /// `c_s` per cohort-sample-idx. `None` for samples whose batch is
    /// singleton or below-floor — Stage 6 reads `None` as `c_s = 0`.
    pub c_s_per_sample: Vec<Option<f64>>,
    /// `q_b` per batch-idx. Each entry is a probability vector over
    /// [`AlleleClass`] indices. Singleton / below-floor batches carry
    /// a zero vector (unused — the samples in those batches have
    /// `c_s = None`).
    pub q_b_per_batch: Vec<[f64; N_ALLELE_CLASSES]>,
    /// Cohort-sample-idx → batch-idx mapping, as supplied at
    /// construction.
    pub sample_to_batch: Vec<usize>,
    /// Provenance.
    pub source: ContaminationEstimateSource,
}

impl ContaminationEstimates {
    /// All-zero estimates: every sample has `c_s = None` (equivalent
    /// to `c_s = 0`), every batch has a zero `q_b` vector. Used by
    /// Stage 6 callers that want to take the contamination-aware code
    /// path with a no-op estimate.
    pub fn zero(n_samples: usize, sample_to_batch: Vec<usize>, n_batches: usize) -> Self {
        debug_assert_eq!(sample_to_batch.len(), n_samples);
        Self {
            c_s_per_sample: vec![None; n_samples],
            q_b_per_batch: vec![[0.0; N_ALLELE_CLASSES]; n_batches],
            sample_to_batch,
            source: ContaminationEstimateSource::Disabled,
        }
    }

    /// Build from user-supplied `c_s` values (one per sample) and
    /// optional user-supplied `q_b` vectors. The provenance is
    /// [`ContaminationEstimateSource::UserSupplied`].
    ///
    /// Used by the CLI path
    /// `--contamination-estimates`/`--contamination-source-distributions`
    /// to skip the side-pass entirely.
    pub fn from_user_supplied(
        c_s_per_sample: Vec<Option<f64>>,
        q_b_per_batch: Vec<[f64; N_ALLELE_CLASSES]>,
        sample_to_batch: Vec<usize>,
    ) -> Self {
        Self {
            c_s_per_sample,
            q_b_per_batch,
            sample_to_batch,
            source: ContaminationEstimateSource::UserSupplied,
        }
    }

    /// Effective `c_s` for a sample — `0.0` whenever
    /// `c_s_per_sample[sample_idx]` is `None` (singleton / below-floor
    /// batch). Stage 6 uses this accessor; callers should not branch
    /// on the `Option` themselves.
    pub fn effective_c_s(&self, sample_idx: usize) -> f64 {
        self.c_s_per_sample[sample_idx].unwrap_or(0.0)
    }

    /// `q_b` for the batch the given sample belongs to. Returns a
    /// reference into [`Self::q_b_per_batch`].
    pub fn q_b_for_sample(&self, sample_idx: usize) -> &[f64; N_ALLELE_CLASSES] {
        &self.q_b_per_batch[self.sample_to_batch[sample_idx]]
    }
}

/// Recommendation surfaced alongside the two non-convergence error
/// variants. Carries a free-form message for the user plus an
/// optional concrete `--contamination-num-sites` value the user can
/// pass to recover.
#[derive(Debug, Clone, PartialEq)]
pub struct NextRunRecommendation {
    pub message: String,
    pub suggested_num_sites: Option<u32>,
}

impl fmt::Display for NextRunRecommendation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

/// Errors the side-pass can surface. Latches on the first error: the
/// caller should treat an `Err` return as terminal and re-invoke the
/// side-pass after addressing the cause.
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum ContaminationEstimationError {
    /// Upstream merger / reader error.
    #[error("upstream merger: {0}")]
    Upstream(#[from] PerPositionMergerError),

    /// Convergence mode, upstream exhausted before convergence fired.
    #[error(
        "contamination estimate did not converge after {sites_processed} \
         informative sites; max per-sample delta in last two blocks = {max_delta:e} \
         (tolerance = {tolerance:e}). {recommendation}"
    )]
    DidNotConverge {
        sites_processed: u32,
        max_delta: f64,
        tolerance: f64,
        /// Per-sample deltas at the last block snapshot. Length is the
        /// cohort size. NaN for samples that never received any reads
        /// (singleton / below-floor / no informative sites in the
        /// stream so far).
        per_sample_deltas: Vec<f64>,
        recommendation: NextRunRecommendation,
    },

    /// Fixed-N mode, upstream exhausted before `num_sites` informative
    /// sites were found.
    #[error(
        "contamination estimate found only {found} informative sites; \
         {requested} were requested. {recommendation}"
    )]
    InsufficientSites {
        requested: u32,
        found: u32,
        recommendation: NextRunRecommendation,
    },

    /// Caller-side input violates an invariant — e.g. `sample_to_batch`
    /// length doesn't match the cohort size, or a `batch_idx` is out
    /// of range.
    #[error("bad input: {0}")]
    BadInput(String),
}

// ---------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------

/// Drive the side-pass to completion against the supplied upstream
/// per-position stream. The function consumes the upstream iterator
/// up to its stopping criterion and returns a single
/// [`ContaminationEstimates`] hand-off value.
///
/// `upstream` is generic so callers can compose with a
/// [`PerPositionMerger`] in production or with hand-built fixtures in
/// tests. The merger is the production source.
///
/// [`PerPositionMerger`]: crate::var_calling::per_position_merger::PerPositionMerger
pub fn estimate_contamination<I>(
    upstream: I,
    n_samples: usize,
    sample_to_batch: Vec<usize>,
    n_batches: usize,
    config: ContaminationEstimationConfig,
) -> Result<ContaminationEstimates, ContaminationEstimationError>
where
    I: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
{
    if sample_to_batch.len() != n_samples {
        return Err(ContaminationEstimationError::BadInput(format!(
            "sample_to_batch length {} does not match n_samples {}",
            sample_to_batch.len(),
            n_samples
        )));
    }
    if let Some(&max_batch) = sample_to_batch.iter().max()
        && max_batch >= n_batches
    {
        return Err(ContaminationEstimationError::BadInput(format!(
            "sample_to_batch contains batch_idx {} >= n_batches {}",
            max_batch, n_batches
        )));
    }
    if config.block_size == 0 {
        return Err(ContaminationEstimationError::BadInput(
            "block_size must be > 0".to_string(),
        ));
    }

    let mut state = OnlineEmState::new(n_samples, n_batches, &config);
    let mut sites_processed: u32 = 0;
    let mut at_block_boundary;

    for item in upstream {
        let pileups = item?;
        let observed_count = filter_and_update(&pileups, &sample_to_batch, &mut state, &config);
        if observed_count == 0 {
            // Site passed Step 1a but no (sample, site) pair passed
            // Step 1b: it did not advance the estimate, so it does not
            // tick the informative-site counter. Per the plan's
            // "site-counting interpretation (B)" decision.
            continue;
        }
        sites_processed += 1;
        at_block_boundary = sites_processed.is_multiple_of(config.block_size);

        if at_block_boundary {
            let (max_delta, per_sample_deltas) =
                refresh_parameters_and_snapshot(&mut state, &config);

            match config.stopping_mode {
                StoppingMode::Convergence {
                    tolerance,
                    stability_blocks,
                } => {
                    if max_delta < tolerance {
                        state.consecutive_within_tolerance += 1;
                        if state.consecutive_within_tolerance >= stability_blocks {
                            return Ok(finalise(state, sample_to_batch, config, sites_processed));
                        }
                    } else {
                        state.consecutive_within_tolerance = 0;
                    }
                    state.last_block_deltas = per_sample_deltas;
                    state.last_block_max_delta = max_delta;
                }
                StoppingMode::FixedSites { num_sites } => {
                    if sites_processed >= num_sites {
                        return Ok(finalise(state, sample_to_batch, config, sites_processed));
                    }
                }
            }
        }
    }

    // Upstream exhausted without hitting the stopping criterion.
    match config.stopping_mode {
        StoppingMode::Convergence {
            tolerance,
            stability_blocks: _,
        } => Err(ContaminationEstimationError::DidNotConverge {
            sites_processed,
            max_delta: state.last_block_max_delta,
            tolerance,
            per_sample_deltas: state.last_block_deltas,
            recommendation: recommend_for_non_convergence(
                sites_processed,
                state.last_block_max_delta,
                tolerance,
            ),
        }),
        StoppingMode::FixedSites { num_sites } => {
            Err(ContaminationEstimationError::InsufficientSites {
                requested: num_sites,
                found: sites_processed,
                recommendation: recommend_for_insufficient_sites(sites_processed, num_sites),
            })
        }
    }
}

// ---------------------------------------------------------------------
// Online-EM state
// ---------------------------------------------------------------------

/// Running sufficient statistics + frozen parameter snapshots. The
/// only mutable state the side-pass keeps; size is
/// `O(n_samples + n_batches × N_ALLELE_CLASSES)`.
struct OnlineEmState {
    /// `S_s = Σ over reads r for sample s of γ_r`.
    s_per_sample: Vec<f64>,
    /// `N_s = total reads observed for sample s (at qualifying sites)`.
    n_per_sample: Vec<u64>,
    /// `S_b[a] = Σ over reads r of allele class a in batch b of γ_r`.
    s_per_batch: Vec<[f64; N_ALLELE_CLASSES]>,
    /// `N_b = total reads observed for batch b`.
    n_per_batch: Vec<u64>,
    /// `c_s_frozen` for the current block's γ_r computations.
    c_s_frozen: Vec<f64>,
    /// `q_b_frozen` for the current block's γ_r computations.
    q_b_frozen: Vec<[f64; N_ALLELE_CLASSES]>,
    /// Snapshot of `c_s` at the previous block boundary; used to
    /// compute per-sample deltas at the next boundary.
    c_s_last_snapshot: Vec<f64>,
    /// Per-sample deltas at the last block boundary. Reported in the
    /// `DidNotConverge` error.
    last_block_deltas: Vec<f64>,
    /// Largest per-sample delta at the last block boundary.
    last_block_max_delta: f64,
    /// Count of consecutive within-tolerance blocks (convergence mode).
    consecutive_within_tolerance: u32,
}

impl OnlineEmState {
    fn new(n_samples: usize, n_batches: usize, _config: &ContaminationEstimationConfig) -> Self {
        Self {
            s_per_sample: vec![0.0; n_samples],
            n_per_sample: vec![0; n_samples],
            s_per_batch: vec![[0.0; N_ALLELE_CLASSES]; n_batches],
            n_per_batch: vec![0; n_batches],
            c_s_frozen: vec![DEFAULT_CS_INIT; n_samples],
            q_b_frozen: vec![[1.0 / N_ALLELE_CLASSES as f64; N_ALLELE_CLASSES]; n_batches],
            c_s_last_snapshot: vec![DEFAULT_CS_INIT; n_samples],
            last_block_deltas: vec![0.0; n_samples],
            last_block_max_delta: f64::INFINITY,
            consecutive_within_tolerance: 0,
        }
    }
}

// ---------------------------------------------------------------------
// Step 1 — informative-site filter
// ---------------------------------------------------------------------

/// Outcome of Step 1a (cohort-wide filter) at a single position.
#[derive(Debug, Clone, Copy)]
struct CohortFilterOutcome {
    /// `true` iff the site passes Step 1a — i.e. is polymorphic in the
    /// cohort at the configured thresholds.
    qualifies: bool,
}

/// Sum allele counts across all samples present at a position, then
/// apply Step 1a's cohort-wide thresholds.
fn step_1a_cohort_filter(
    pileups: &PerPositionPileups,
    config: &ContaminationEstimationConfig,
) -> CohortFilterOutcome {
    // The position is shared across all samples; the allele alphabet
    // is per-sample. Aggregating across samples means summing
    // num_obs across entries with matching `seq`. The expected
    // alphabet size is small (typically 2–4), so a linear scan over a
    // small Vec is the right shape.
    let mut cohort_alleles: Vec<(Vec<u8>, u32)> = Vec::with_capacity(4);
    for slot in &pileups.per_sample {
        let Some(rec) = slot else { continue };
        for obs in &rec.alleles {
            if let Some(entry) = cohort_alleles.iter_mut().find(|(seq, _)| seq == &obs.seq) {
                entry.1 = entry.1.saturating_add(obs.support.num_obs);
            } else {
                cohort_alleles.push((obs.seq.clone(), obs.support.num_obs));
            }
        }
    }

    if cohort_alleles.len() < 2 {
        // Only one allele type observed → monomorphic in the cohort.
        return CohortFilterOutcome { qualifies: false };
    }

    // Sort descending by count to identify major and second-most-common.
    cohort_alleles.sort_by_key(|entry| std::cmp::Reverse(entry.1));
    let total_reads: u64 = cohort_alleles.iter().map(|(_, n)| u64::from(*n)).sum();
    if total_reads == 0 {
        return CohortFilterOutcome { qualifies: false };
    }
    let second_count = cohort_alleles[1].1;
    let cohort_minor_fraction = u64::from(second_count) as f64 / total_reads as f64;

    let qualifies = second_count >= config.min_cohort_minor_count
        && cohort_minor_fraction >= config.min_cohort_minor_fraction;
    CohortFilterOutcome { qualifies }
}

/// Per-sample observation that passes Step 1b. Each retained
/// `(sample, site)` pair produces one of these; the online-EM update
/// folds it into the running sufficient statistics immediately.
#[derive(Debug, Clone, Copy)]
struct InformativeObservation {
    sample_idx: usize,
    /// Per-allele-class read counts at this (sample, site). Reads of
    /// allele classes the sample's observation didn't carry contribute
    /// 0; only the classes actually seen here are non-zero.
    counts_by_class: [u32; N_ALLELE_CLASSES],
    /// Per-allele-class summed log-error scalars at this (sample,
    /// site). Used to derive the per-class mean base-error rate.
    q_sum_by_class: [f64; N_ALLELE_CLASSES],
    /// Which allele class is the major (g_major). Used to compute
    /// `P(read | own genotype)` in the mixture likelihood.
    g_major_class: AlleleClass,
}

/// Apply Step 1b to a single sample's record at a qualifying position
/// and return an `InformativeObservation` if the sample passes.
///
/// Allele classification uses `record.alleles[0].seq` as the REF
/// sequence — that slot is unconditionally REF by the walker
/// invariant documented on [`PileupRecord`].
fn step_1b_sample_filter(
    record: &PileupRecord,
    config: &ContaminationEstimationConfig,
) -> Option<InformativeObservation> {
    if record.alleles.is_empty() {
        return None;
    }

    let depth: u32 = record.alleles.iter().map(|a| a.support.num_obs).sum();
    if depth < config.min_depth {
        return None;
    }

    let (major_idx, major_count) = record
        .alleles
        .iter()
        .enumerate()
        .map(|(i, a)| (i, a.support.num_obs))
        .max_by_key(|&(_, n)| n)?;

    let major_fraction = u64::from(major_count) as f64 / u64::from(depth) as f64;
    if major_fraction < config.min_major_fraction {
        return None;
    }

    let ref_seq = &record.alleles[0].seq;
    let g_major_class = AlleleClass::classify(&record.alleles[major_idx].seq, ref_seq);

    let mut counts_by_class = [0_u32; N_ALLELE_CLASSES];
    let mut q_sum_by_class = [0.0_f64; N_ALLELE_CLASSES];
    for obs in &record.alleles {
        let class = AlleleClass::classify(&obs.seq, ref_seq);
        counts_by_class[class as usize] =
            counts_by_class[class as usize].saturating_add(obs.support.num_obs);
        q_sum_by_class[class as usize] += obs.support.q_sum;
    }

    Some(InformativeObservation {
        sample_idx: 0, // filled in by the caller (we don't have sample_idx here)
        counts_by_class,
        q_sum_by_class,
        g_major_class,
    })
}

/// Walk the per-sample slots at one qualifying position, accumulate
/// each passing observation into the online-EM state, and return how
/// many samples contributed (used to decide whether the site ticks
/// the informative-site counter).
fn filter_and_update(
    pileups: &PerPositionPileups,
    sample_to_batch: &[usize],
    state: &mut OnlineEmState,
    config: &ContaminationEstimationConfig,
) -> u32 {
    let cohort = step_1a_cohort_filter(pileups, config);
    if !cohort.qualifies {
        return 0;
    }

    let mut contributed: u32 = 0;
    for (sample_idx, slot) in pileups.per_sample.iter().enumerate() {
        let Some(record) = slot else { continue };
        let Some(mut obs) = step_1b_sample_filter(record, config) else {
            continue;
        };
        obs.sample_idx = sample_idx;
        apply_observation(state, &obs, sample_to_batch);
        contributed += 1;
    }
    contributed
}

// ---------------------------------------------------------------------
// Step 3 + Step 4 — γ computation and sufficient-statistic update
// ---------------------------------------------------------------------

/// Fold one informative observation into the running sufficient
/// statistics. Implements the per-observation E-step (γ at frozen
/// parameters) and the M-step contribution from this observation in
/// one pass — there is no per-block read buffer.
fn apply_observation(
    state: &mut OnlineEmState,
    obs: &InformativeObservation,
    sample_to_batch: &[usize],
) {
    let s = obs.sample_idx;
    let b = sample_to_batch[s];
    let c_s = state.c_s_frozen[s];
    let q_b = state.q_b_frozen[b];

    let total_reads: u32 = obs.counts_by_class.iter().sum();
    if total_reads == 0 {
        return;
    }

    // Depth-weighted mean per-base error for this sample at this site.
    // `q_sum` is in log-error space already (ln of per-read error
    // probability), so the per-read mean is exp(mean log error) =
    // exp(sum / count).
    let total_q_sum: f64 = obs.q_sum_by_class.iter().sum();
    let mean_log_error = if total_reads > 0 {
        total_q_sum / total_reads as f64
    } else {
        0.0
    };
    let mean_error = mean_log_error.exp().clamp(MIN_BASE_ERROR, MAX_BASE_ERROR);

    for (class_idx, &n_reads_class) in obs.counts_by_class.iter().enumerate() {
        if n_reads_class == 0 {
            continue;
        }
        // P(read carrying allele of this class | own genotype is
        // hom-major in g_major_class). Spec's "uniform-error" model:
        // the major class gets (1 − ε); the other classes share ε
        // uniformly.
        let p_own = if class_idx == obs.g_major_class as usize {
            1.0 - mean_error
        } else {
            mean_error / (N_ALLELE_CLASSES as f64 - 1.0)
        };
        let p_contam = q_b[class_idx];

        // Mixture denominator, floored to avoid numerical 0/0 in the
        // degenerate case where both p_own and p_contam are tiny.
        let mix = (1.0 - c_s) * p_own + c_s * p_contam;
        let gamma = if mix > f64::MIN_POSITIVE {
            (c_s * p_contam) / mix
        } else {
            // Defensive: at p_own=p_contam=0 gamma is undefined; treat
            // the contribution as zero so it does not perturb the
            // running statistics.
            0.0
        };

        let weight = gamma * f64::from(n_reads_class);
        state.s_per_sample[s] += weight;
        state.n_per_sample[s] += u64::from(n_reads_class);
        state.s_per_batch[b][class_idx] += weight;
        state.n_per_batch[b] += u64::from(n_reads_class);
    }
}

/// Recompute the parameter estimates from the running sufficient
/// statistics; replace the frozen snapshots for the next block; update
/// the last-block delta bookkeeping; return the largest per-sample
/// delta and the per-sample delta vector.
fn refresh_parameters_and_snapshot(
    state: &mut OnlineEmState,
    config: &ContaminationEstimationConfig,
) -> (f64, Vec<f64>) {
    // M-step on c_s: simple ratio per sample, with the init value
    // retained for samples that never received any reads (so no
    // division by zero and no spurious jump from init to a noisy ratio).
    let mut c_s_new = state.c_s_frozen.clone();
    for (slot, (&s_val, &n_val)) in c_s_new
        .iter_mut()
        .zip(state.s_per_sample.iter().zip(state.n_per_sample.iter()))
    {
        if n_val > 0 {
            *slot = (s_val / n_val as f64).clamp(0.0, 1.0);
        }
    }

    // M-step on q_b: Dirichlet-smoothed per-batch fractions.
    let alpha = [
        config.ref_pseudocount,
        config.snp_alt_pseudocount,
        config.indel_alt_pseudocount,
    ];
    let alpha_sum: f64 = alpha.iter().sum();
    let mut q_b_new = state.q_b_frozen.clone();
    for (q_b_slot, s_b_row) in q_b_new.iter_mut().zip(state.s_per_batch.iter()) {
        let s_sum: f64 = s_b_row.iter().sum();
        let denom = alpha_sum + s_sum;
        if denom > 0.0 {
            for (a, slot) in q_b_slot.iter_mut().enumerate() {
                *slot = (alpha[a] + s_b_row[a]) / denom;
            }
        }
        // else: keep the init uniform — no reads seen for this batch yet.
    }

    // Compute per-sample deltas against the previous snapshot.
    let mut per_sample_deltas = vec![0.0_f64; c_s_new.len()];
    let mut max_delta = 0.0_f64;
    for s in 0..c_s_new.len() {
        let delta = (c_s_new[s] - state.c_s_last_snapshot[s]).abs();
        per_sample_deltas[s] = delta;
        if delta > max_delta {
            max_delta = delta;
        }
    }

    // Promote new params to frozen snapshots and record the snapshot.
    state.c_s_frozen = c_s_new.clone();
    state.q_b_frozen = q_b_new;
    state.c_s_last_snapshot = c_s_new;
    state.last_block_max_delta = max_delta;
    state.last_block_deltas = per_sample_deltas.clone();

    (max_delta, per_sample_deltas)
}

// ---------------------------------------------------------------------
// Step 5 — finalisation
// ---------------------------------------------------------------------

/// Apply the singleton/below-floor floor and return the
/// `ContaminationEstimates` hand-off. Called once when the stopping
/// criterion fires.
fn finalise(
    state: OnlineEmState,
    sample_to_batch: Vec<usize>,
    config: ContaminationEstimationConfig,
    sites_processed: u32,
) -> ContaminationEstimates {
    let n_samples = state.c_s_frozen.len();
    let n_batches = state.q_b_frozen.len();

    // Identify batches at or above the floor (and not singletons).
    let mut batch_sizes = vec![0_u32; n_batches];
    for &b in &sample_to_batch {
        batch_sizes[b] = batch_sizes[b].saturating_add(1);
    }
    let batch_qualifies: Vec<bool> = batch_sizes
        .iter()
        .map(|&n| n >= config.min_batch_size_for_contamination.max(2))
        .collect();

    let c_s_per_sample: Vec<Option<f64>> = (0..n_samples)
        .map(|s| {
            let b = sample_to_batch[s];
            if batch_qualifies[b] {
                Some(state.c_s_frozen[s])
            } else {
                None
            }
        })
        .collect();

    let q_b_per_batch: Vec<[f64; N_ALLELE_CLASSES]> = (0..n_batches)
        .map(|b| {
            if batch_qualifies[b] {
                state.q_b_frozen[b]
            } else {
                [0.0; N_ALLELE_CLASSES]
            }
        })
        .collect();

    ContaminationEstimates {
        c_s_per_sample,
        q_b_per_batch,
        sample_to_batch,
        source: ContaminationEstimateSource::SidePass {
            mode: config.stopping_mode,
            sites_processed,
        },
    }
}

// ---------------------------------------------------------------------
// Recommendation builders
// ---------------------------------------------------------------------

fn recommend_for_non_convergence(
    sites_processed: u32,
    observed_delta: f64,
    tolerance: f64,
) -> NextRunRecommendation {
    let suggested = if observed_delta.is_finite() && observed_delta > 0.0 && tolerance > 0.0 {
        // Under the expected MLE scaling delta ∝ N^(−1/2), the N at
        // which the delta drops to `tolerance` is
        //   N_recommended ≈ N · (observed_delta / tolerance)^2.
        // Cap at 10× to keep the recommendation actionable.
        let ratio = observed_delta / tolerance;
        let scaled = f64::from(sites_processed) * ratio * ratio;
        let capped = scaled.min(f64::from(sites_processed) * 10.0);
        if capped.is_finite() && capped > 0.0 {
            Some(capped.round().min(u32::MAX as f64) as u32)
        } else {
            None
        }
    } else {
        None
    };
    let suggested_clause = match suggested {
        Some(n) => format!(
            " — try `--contamination-num-sites {}` or relax `--contamination-stability-tolerance`",
            n
        ),
        None => String::new(),
    };
    NextRunRecommendation {
        message: format!(
            "Recommended next run: relax `--contamination-stability-tolerance`, \
             lower `--contamination-min-major-fraction`, or switch to fixed-N mode{}",
            suggested_clause
        ),
        suggested_num_sites: suggested,
    }
}

fn recommend_for_insufficient_sites(found: u32, requested: u32) -> NextRunRecommendation {
    let suggested = if found > 0 { Some(found) } else { None };
    let suggested_clause = match suggested {
        Some(n) => format!(" — try `--contamination-num-sites {}`", n),
        None => " — try a smaller `--contamination-num-sites`".to_string(),
    };
    NextRunRecommendation {
        message: format!(
            "Recommended next run: lower `--contamination-num-sites` from \
             {} to ~{} (the observed count), switch to convergence mode \
             with `--contamination-stability-tolerance`, or relax filtering \
             thresholds (e.g. `--contamination-min-major-fraction`){}",
            requested, found, suggested_clause
        ),
        suggested_num_sites: suggested,
    }
}

// ---------------------------------------------------------------------
// Numerical guards
// ---------------------------------------------------------------------

/// Floor on per-base error rate. Some BQ aggregates can underflow to
/// near zero for very high-quality reads; clamping prevents the
/// mixture denominator from collapsing.
const MIN_BASE_ERROR: f64 = 1e-12;

/// Ceiling on per-base error rate. The model assumes hom-major reads,
/// so a sample with mean error near 1.0 would not have passed
/// Step 1b's `min_major_fraction` filter; this is a safety net.
const MAX_BASE_ERROR: f64 = 0.5;

// =====================================================================
// Tests
// =====================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::per_sample_pileup::pileup::{AlleleObservation, AlleleSupportStats, PileupRecord};
    use proptest::prelude::*;

    // ----- Fixture builders ----------------------------------------

    /// Build a synthetic allele observation. `seq` is the raw allele
    /// bytes; the REF/ALT distinction is positional inside a
    /// `PileupRecord` (slot 0 is always REF), so this helper does not
    /// take a `is_reference` flag.
    fn obs(seq: &[u8], num_obs: u32, mean_bq_err: f64) -> AlleleObservation {
        let q_sum = if num_obs == 0 {
            0.0
        } else {
            mean_bq_err.ln() * f64::from(num_obs)
        };
        AlleleObservation {
            seq: seq.to_vec(),
            support: AlleleSupportStats {
                num_obs,
                q_sum,
                fwd: num_obs / 2,
                placed_left: num_obs / 2,
                placed_start: 0,
            },
            chain_ids: Vec::new(),
        }
    }

    fn pileup_record(pos: u32, alleles: Vec<AlleleObservation>) -> PileupRecord {
        PileupRecord::new(0, pos, alleles)
    }

    fn position(pos: u32, per_sample: Vec<Option<PileupRecord>>) -> PerPositionPileups {
        PerPositionPileups {
            chrom_id: 0,
            pos,
            per_sample,
        }
    }

    // ----- Step 1a: cohort-wide filter -----------------------------

    #[test]
    fn step_1a_rejects_monomorphic_cohort() {
        // All samples report only REF.
        let pileups = position(
            100,
            vec![
                Some(pileup_record(100, vec![obs(b"A", 30, 0.001)])),
                Some(pileup_record(100, vec![obs(b"A", 25, 0.001)])),
                Some(pileup_record(100, vec![obs(b"A", 40, 0.001)])),
            ],
        );
        let cfg = ContaminationEstimationConfig::with_project_defaults();
        let outcome = step_1a_cohort_filter(&pileups, &cfg);
        assert!(!outcome.qualifies);
    }

    #[test]
    fn step_1a_accepts_polymorphic_cohort() {
        // Two samples REF-heavy, one sample alt-heavy → cohort
        // polymorphic with cohort-minor fraction well above default.
        let pileups = position(
            200,
            vec![
                Some(pileup_record(
                    200,
                    vec![obs(b"A", 30, 0.001), obs(b"C", 1, 0.001)],
                )),
                Some(pileup_record(
                    200,
                    vec![obs(b"A", 25, 0.001), obs(b"C", 1, 0.001)],
                )),
                Some(pileup_record(
                    200,
                    vec![obs(b"A", 5, 0.001), obs(b"C", 30, 0.001)],
                )),
            ],
        );
        let cfg = ContaminationEstimationConfig::with_project_defaults();
        let outcome = step_1a_cohort_filter(&pileups, &cfg);
        assert!(outcome.qualifies);
    }

    #[test]
    fn step_1a_rejects_below_min_cohort_minor_count() {
        // Two samples REF-only, one with a single odd alt read →
        // cohort minor count = 1, below the default threshold of 2.
        let pileups = position(
            300,
            vec![
                Some(pileup_record(300, vec![obs(b"A", 30, 0.001)])),
                Some(pileup_record(300, vec![obs(b"A", 25, 0.001)])),
                Some(pileup_record(
                    300,
                    vec![obs(b"A", 40, 0.001), obs(b"C", 1, 0.001)],
                )),
            ],
        );
        let cfg = ContaminationEstimationConfig::with_project_defaults();
        let outcome = step_1a_cohort_filter(&pileups, &cfg);
        assert!(!outcome.qualifies);
    }

    // ----- Step 1b: per-sample filter ------------------------------

    #[test]
    fn step_1b_rejects_low_depth() {
        let rec = pileup_record(100, vec![obs(b"A", 5, 0.001), obs(b"C", 0, 0.001)]);
        let cfg = ContaminationEstimationConfig::with_project_defaults();
        // depth = 5, default min_depth = 10.
        assert!(step_1b_sample_filter(&rec, &cfg).is_none());
    }

    #[test]
    fn step_1b_rejects_heterozygous_signal() {
        // 60 / 40 split — major fraction = 0.6, well below 0.95.
        let rec = pileup_record(100, vec![obs(b"A", 60, 0.001), obs(b"C", 40, 0.001)]);
        let cfg = ContaminationEstimationConfig::with_project_defaults();
        assert!(step_1b_sample_filter(&rec, &cfg).is_none());
    }

    #[test]
    fn step_1b_accepts_confident_hom_major() {
        // 99 REF, 1 alt → major fraction = 0.99, depth = 100.
        let rec = pileup_record(100, vec![obs(b"A", 99, 0.001), obs(b"C", 1, 0.001)]);
        let cfg = ContaminationEstimationConfig::with_project_defaults();
        let obs = step_1b_sample_filter(&rec, &cfg).expect("should qualify");
        assert_eq!(obs.g_major_class, AlleleClass::Ref);
        assert_eq!(obs.counts_by_class[AlleleClass::Ref as usize], 99);
        assert_eq!(obs.counts_by_class[AlleleClass::SnpAlt as usize], 1);
    }

    // ----- Online EM ground-truth tests ----------------------------

    /// Build a stream of N synthetic informative sites where sample 0
    /// is hom-REF with `n_reads` reads per site, of which a `c_true`
    /// fraction are switched to the SNP alt allele (simulating
    /// contamination from another sample). All other samples are
    /// hom-REF with no contamination.
    fn synth_two_sample_contam_stream(
        n_sites: u32,
        n_reads: u32,
        c_true: f64,
    ) -> Vec<Result<PerPositionPileups, PerPositionMergerError>> {
        let mut out = Vec::with_capacity(n_sites as usize);
        // Use deterministic rounding to keep the test reproducible.
        let n_contam = (f64::from(n_reads) * c_true).round() as u32;
        let n_own = n_reads - n_contam;
        let mean_err = 0.001;
        for i in 0..n_sites {
            // Sample 0: contaminated.
            let s0 = pileup_record(
                i + 1,
                vec![obs(b"A", n_own, mean_err), obs(b"C", n_contam, mean_err)],
            );
            // Sample 1: clean hom-REF.
            let s1 = pileup_record(i + 1, vec![obs(b"A", n_reads, mean_err)]);
            // Sample 2: clean alt-carrier (provides cohort polymorphism so
            // Step 1a accepts).
            let s2 = pileup_record(
                i + 1,
                vec![obs(b"A", 0, mean_err), obs(b"C", n_reads, mean_err)],
            );
            out.push(Ok(position(i + 1, vec![Some(s0), Some(s1), Some(s2)])));
        }
        out
    }

    #[test]
    fn online_em_recovers_zero_contamination() {
        // Three samples in one batch, no contamination injected.
        let stream = synth_two_sample_contam_stream(2000, 50, 0.0);
        let estimates = estimate_contamination(
            stream.into_iter(),
            3,
            vec![0, 0, 0],
            1,
            ContaminationEstimationConfig {
                stopping_mode: StoppingMode::FixedSites { num_sites: 1500 },
                block_size: 500,
                min_batch_size_for_contamination: 2,
                ..ContaminationEstimationConfig::with_project_defaults()
            },
        )
        .expect("should succeed");
        // Sample 0 (the only hom-REF carrier with informative reads) should
        // estimate c_s near 0. With base-error rate 0.001 and 50 reads/site
        // there is some residual signal so we use a loose tolerance.
        let c_s_0 = estimates.c_s_per_sample[0].expect("not floored");
        assert!(c_s_0 < 0.05, "c_s for clean sample = {}", c_s_0);
    }

    #[test]
    fn online_em_recovers_three_percent_contamination() {
        // Sample 0 contaminated at 3 %; should recover within tolerance.
        let stream = synth_two_sample_contam_stream(4000, 50, 0.03);
        let estimates = estimate_contamination(
            stream.into_iter(),
            3,
            vec![0, 0, 0],
            1,
            ContaminationEstimationConfig {
                stopping_mode: StoppingMode::FixedSites { num_sites: 3000 },
                block_size: 500,
                min_batch_size_for_contamination: 2,
                ..ContaminationEstimationConfig::with_project_defaults()
            },
        )
        .expect("should succeed");
        let c_s_0 = estimates.c_s_per_sample[0].expect("not floored");
        // The online estimator has some init-block bias; ±0.02 is
        // generous and matches the spec's expected accuracy window
        // for a 3 000-site subsample.
        assert!(
            (c_s_0 - 0.03).abs() < 0.02,
            "c_s for contaminated sample = {}; expected ~0.03",
            c_s_0
        );
    }

    #[test]
    fn fixed_n_mode_processes_exact_count() {
        let stream = synth_two_sample_contam_stream(5000, 50, 0.02);
        let estimates = estimate_contamination(
            stream.into_iter(),
            3,
            vec![0, 0, 0],
            1,
            ContaminationEstimationConfig {
                stopping_mode: StoppingMode::FixedSites { num_sites: 2000 },
                block_size: 500,
                min_batch_size_for_contamination: 2,
                ..ContaminationEstimationConfig::with_project_defaults()
            },
        )
        .expect("should succeed");
        match estimates.source {
            ContaminationEstimateSource::SidePass {
                mode: StoppingMode::FixedSites { num_sites },
                sites_processed,
            } => {
                assert_eq!(num_sites, 2000);
                assert_eq!(sites_processed, 2000);
            }
            other => panic!("unexpected source: {other:?}"),
        }
    }

    #[test]
    fn convergence_mode_returns_did_not_converge_on_short_stream() {
        // Very short stream — convergence cannot fire in 3 blocks of
        // 500 sites = 1500 sites, and we provide only 200 sites.
        let stream = synth_two_sample_contam_stream(200, 50, 0.02);
        let result = estimate_contamination(
            stream.into_iter(),
            3,
            vec![0, 0, 0],
            1,
            ContaminationEstimationConfig {
                stopping_mode: StoppingMode::Convergence {
                    tolerance: 1e-4,
                    stability_blocks: 3,
                },
                block_size: 500,
                min_batch_size_for_contamination: 2,
                ..ContaminationEstimationConfig::with_project_defaults()
            },
        );
        match result {
            Err(ContaminationEstimationError::DidNotConverge {
                sites_processed, ..
            }) => {
                assert!(sites_processed <= 200);
            }
            other => panic!("expected DidNotConverge, got {other:?}"),
        }
    }

    #[test]
    fn fixed_n_returns_insufficient_sites_on_short_stream() {
        let stream = synth_two_sample_contam_stream(100, 50, 0.02);
        let result = estimate_contamination(
            stream.into_iter(),
            3,
            vec![0, 0, 0],
            1,
            ContaminationEstimationConfig {
                stopping_mode: StoppingMode::FixedSites { num_sites: 5000 },
                block_size: 500,
                min_batch_size_for_contamination: 2,
                ..ContaminationEstimationConfig::with_project_defaults()
            },
        );
        match result {
            Err(ContaminationEstimationError::InsufficientSites {
                requested, found, ..
            }) => {
                assert_eq!(requested, 5000);
                assert!(found <= 100);
            }
            other => panic!("expected InsufficientSites, got {other:?}"),
        }
    }

    // ----- Step 5: floors and overrides ----------------------------

    #[test]
    fn singleton_batch_gets_c_s_none() {
        // Two samples, each in its own batch → both singleton.
        let stream = synth_two_sample_contam_stream(2000, 50, 0.03);
        let estimates = estimate_contamination(
            stream.into_iter(),
            3,
            vec![0, 1, 2], // three singleton batches
            3,
            ContaminationEstimationConfig {
                stopping_mode: StoppingMode::FixedSites { num_sites: 1500 },
                block_size: 500,
                min_batch_size_for_contamination: 5,
                ..ContaminationEstimationConfig::with_project_defaults()
            },
        )
        .expect("should succeed");
        for c in &estimates.c_s_per_sample {
            assert!(c.is_none(), "singleton batches should be floored");
        }
        for q in &estimates.q_b_per_batch {
            assert_eq!(*q, [0.0; N_ALLELE_CLASSES]);
        }
    }

    #[test]
    fn below_floor_batch_gets_c_s_none() {
        // All three samples in one batch with floor = 5 → below floor.
        let stream = synth_two_sample_contam_stream(2000, 50, 0.03);
        let estimates = estimate_contamination(
            stream.into_iter(),
            3,
            vec![0, 0, 0],
            1,
            ContaminationEstimationConfig {
                stopping_mode: StoppingMode::FixedSites { num_sites: 1500 },
                block_size: 500,
                min_batch_size_for_contamination: 5,
                ..ContaminationEstimationConfig::with_project_defaults()
            },
        )
        .expect("should succeed");
        for c in &estimates.c_s_per_sample {
            assert!(c.is_none(), "below-floor batches should be floored");
        }
    }

    #[test]
    fn cs_in_unit_interval_and_qb_is_simplex() {
        let stream = synth_two_sample_contam_stream(2000, 50, 0.04);
        let estimates = estimate_contamination(
            stream.into_iter(),
            3,
            vec![0, 0, 0],
            1,
            ContaminationEstimationConfig {
                stopping_mode: StoppingMode::FixedSites { num_sites: 1500 },
                block_size: 500,
                min_batch_size_for_contamination: 2,
                ..ContaminationEstimationConfig::with_project_defaults()
            },
        )
        .expect("should succeed");
        for c in estimates.c_s_per_sample.iter().flatten() {
            assert!((0.0..=1.0).contains(c), "c_s out of range: {c}");
        }
        for q in &estimates.q_b_per_batch {
            let sum: f64 = q.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-9 || sum == 0.0,
                "q_b not a simplex: sum = {sum}"
            );
            for v in q {
                assert!((0.0..=1.0).contains(v), "q_b entry out of range: {v}");
            }
        }
    }

    // ----- Provenance ----------------------------------------------

    #[test]
    fn estimates_zero_constructor_is_disabled_source() {
        let est = ContaminationEstimates::zero(3, vec![0, 0, 0], 1);
        assert_eq!(est.source, ContaminationEstimateSource::Disabled);
        assert!(est.c_s_per_sample.iter().all(|c| c.is_none()));
        for s in 0..3 {
            assert_eq!(est.effective_c_s(s), 0.0);
        }
    }

    #[test]
    fn from_user_supplied_is_user_source() {
        let est = ContaminationEstimates::from_user_supplied(
            vec![Some(0.03), None, Some(0.01)],
            vec![[0.9, 0.08, 0.02]],
            vec![0, 0, 0],
        );
        assert_eq!(est.source, ContaminationEstimateSource::UserSupplied);
        assert_eq!(est.effective_c_s(0), 0.03);
        assert_eq!(est.effective_c_s(1), 0.0);
    }

    // ----- Bad-input contract --------------------------------------

    #[test]
    fn bad_input_when_sample_to_batch_length_wrong() {
        let stream: Vec<Result<PerPositionPileups, PerPositionMergerError>> = vec![];
        let result = estimate_contamination(
            stream.into_iter(),
            3,
            vec![0, 0], // wrong length
            1,
            ContaminationEstimationConfig::with_project_defaults(),
        );
        assert!(matches!(
            result,
            Err(ContaminationEstimationError::BadInput(_))
        ));
    }

    #[test]
    fn bad_input_when_batch_idx_out_of_range() {
        let stream: Vec<Result<PerPositionPileups, PerPositionMergerError>> = vec![];
        let result = estimate_contamination(
            stream.into_iter(),
            3,
            vec![0, 1, 5], // batch_idx 5 with n_batches = 2
            2,
            ContaminationEstimationConfig::with_project_defaults(),
        );
        assert!(matches!(
            result,
            Err(ContaminationEstimationError::BadInput(_))
        ));
    }

    #[test]
    fn bad_input_when_block_size_zero() {
        let stream: Vec<Result<PerPositionPileups, PerPositionMergerError>> = vec![];
        let result = estimate_contamination(
            stream.into_iter(),
            1,
            vec![0],
            1,
            ContaminationEstimationConfig {
                block_size: 0,
                ..ContaminationEstimationConfig::with_project_defaults()
            },
        );
        assert!(matches!(
            result,
            Err(ContaminationEstimationError::BadInput(_))
        ));
    }

    // ----- Proptests -----------------------------------------------

    proptest! {
        #[test]
        fn proptest_c_s_in_unit_interval(
            n_sites in 600u32..1500,
            n_reads in 20u32..80,
            c_true in 0.0_f64..0.10,
            batch_size in 2usize..4,
        ) {
            let stream = synth_two_sample_contam_stream(n_sites, n_reads, c_true);
            let n_samples = 3;
            let effective_batch_size = batch_size.clamp(1, 3);
            let s_to_b: Vec<usize> = (0..n_samples).map(|i| i % effective_batch_size).collect();
            let n_batches = effective_batch_size.min(n_samples);
            let estimates = estimate_contamination(
                stream.into_iter(),
                n_samples,
                s_to_b,
                n_batches,
                ContaminationEstimationConfig {
                    stopping_mode: StoppingMode::FixedSites { num_sites: 500 },
                    block_size: 250,
                    min_batch_size_for_contamination: 2,
                    ..ContaminationEstimationConfig::with_project_defaults()
                },
            );
            // The estimator may produce DidNotConverge or InsufficientSites on
            // pathological samples; both are valid outcomes for the proptest.
            if let Ok(est) = estimates {
                for c in est.c_s_per_sample.iter().flatten() {
                    prop_assert!((0.0..=1.0).contains(c));
                }
            }
        }

        #[test]
        fn proptest_q_b_is_simplex(
            n_sites in 600u32..1500,
            n_reads in 20u32..80,
            c_true in 0.0_f64..0.10,
        ) {
            let stream = synth_two_sample_contam_stream(n_sites, n_reads, c_true);
            let estimates = estimate_contamination(
                stream.into_iter(),
                3,
                vec![0, 0, 0],
                1,
                ContaminationEstimationConfig {
                    stopping_mode: StoppingMode::FixedSites { num_sites: 500 },
                    block_size: 250,
                    min_batch_size_for_contamination: 2,
                    ..ContaminationEstimationConfig::with_project_defaults()
                },
            );
            if let Ok(est) = estimates {
                for q in &est.q_b_per_batch {
                    let sum: f64 = q.iter().sum();
                    prop_assert!((sum - 1.0).abs() < 1e-9 || sum == 0.0);
                    for v in q {
                        prop_assert!((0.0..=1.0).contains(v));
                    }
                }
            }
        }

        #[test]
        fn proptest_block_size_independence_within_tolerance(
            c_true in 0.01_f64..0.06,
        ) {
            // For the same input, different block sizes should
            // produce c_s estimates that agree within a generous
            // tolerance once both converge.
            let n_sites = 3000;
            let n_reads = 50;
            let s_to_b = vec![0_usize, 0, 0];
            let common_cfg = |bs: u32| ContaminationEstimationConfig {
                stopping_mode: StoppingMode::FixedSites { num_sites: 2000 },
                block_size: bs,
                min_batch_size_for_contamination: 2,
                ..ContaminationEstimationConfig::with_project_defaults()
            };
            let stream_a = synth_two_sample_contam_stream(n_sites, n_reads, c_true);
            let stream_b = synth_two_sample_contam_stream(n_sites, n_reads, c_true);
            let est_a = estimate_contamination(stream_a.into_iter(), 3, s_to_b.clone(), 1, common_cfg(500))
                .expect("a converges");
            let est_b = estimate_contamination(stream_b.into_iter(), 3, s_to_b, 1, common_cfg(1000))
                .expect("b converges");
            for s in 0..3 {
                let a = est_a.c_s_per_sample[s].unwrap_or(0.0);
                let b = est_b.c_s_per_sample[s].unwrap_or(0.0);
                prop_assert!(
                    (a - b).abs() < 0.02,
                    "block-size sensitivity at sample {}: {} vs {}",
                    s, a, b
                );
            }
        }
    }
}
