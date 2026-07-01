//! The hidden-paralog filter's consumer statistics.
//!
//! This module holds the *pure* statistics that decide, at each locus,
//! which of two stories better explains all the samples' coverage and
//! allele counts: **H1** a real single-copy variant, or **H2** a hidden
//! (reference-collapsed) paralog. It consumes the per-sample summaries the
//! Stage-1 pileup already stores in the `.psp`
//! ([`crate::sample_summary`]) and depends on nothing in `var_calling`
//! (architecture
//! `doc/devel/architecture/hidden_paralog_locus_statistic.md`, Premise 0),
//! so every function here is unit-testable and parity-checkable in
//! isolation before the var-calling wiring depends inward on it.
//!
//! Build order (implementation plan
//! `doc/devel/implementation_plans/paralog_filter_model.md`, Milestone Q):
//!
//! - **Q1 (this file):** the fixed model constants/grids —
//!   [`ParalogModelParams`].
//! - Q2 `coverage_model` — fit the per-sample single-copy depth scale +
//!   GC-bias curve → a window's `relative_copy_number`.
//! - Q3 `locus_score` — the pure H1-vs-H2 marginal likelihood ratio.
//! - Q4 the inbreeding scalar `F` from observed het + `Hexp`.
//! - Q5 `prior` — the empirical-Bayes prior + FDR that turn the raw
//!   likelihood ratio into a per-locus probability and a cut.
//!
//! The design intent lives in the spec
//! `doc/devel/specs/hidden_paralog_filter.md` (§4–§7); the validated maths
//! reference is the tomato2 prototype
//! (`benchmarks/tomato2/src/build_paralog_lr.py`,
//! `build_paralog_eb.py`, `build_gc_normalization.py`) — a *reference
//! draft*, not a bit-exact oracle.

pub mod coverage_model;
pub mod inbreeding;
pub mod locus_score;

pub use coverage_model::{
    CoverageFitConfig, CoverageModelError, ModeMedianRatioBounds, SingleCopyCoverageModel,
};
pub use inbreeding::{inbreeding_coefficient, obs_het};
pub use locus_score::{
    LocusObservations, ParalogScore, SampleObservation, score_locus_for_paralogy,
};

/// A per-read allele error floor (`ε`): the VAF assigned to a genotype's
/// "should be absent" allele (hom-ref / non-carrier), and, as `1 − ε`, the
/// VAF of a "should be fully present" allele (hom-alt). It keeps the
/// binomial allele factor away from the `ln 0` singularity. Prototype
/// `EPS`.
pub const DEFAULT_PSEUDOCOUNT_VAF: f64 = 0.01;

/// Winsor cap on `relative_copy_number` for the Normal coverage tail: any
/// window at or above this many copies reads as "looks like the top carrier
/// configuration" (`T/2 = 4` for `T = 8`). Above it the model no longer
/// distinguishes 4× from 8×+, which is fine for *detection* (all such
/// windows are paralog-like). If the carrier set or cap is retuned, keep the
/// cap above the top carrier mean + a few σ (arch Premise 2).
pub const DEFAULT_MAX_RELATIVE_COPY_NUMBER: f64 = 4.0;

/// Hidden-paralog total copy numbers `T` enumerated under H2 (coverage
/// `1.5×, 2×, 3×, 4×`; `T ≤ 8` so the top mean equals the winsor cap).
/// Prototype `T_CARRIER`.
pub const DEFAULT_CARRIER_COPY_NUMBERS: &[u32] = &[3, 4, 6, 8];

/// Grid resolution for marginalising the allele frequency `p` under the
/// folded site-frequency-spectrum prior. Dense because the `1/(p(1−p))`
/// prior concentrates weight at low `p`. Prototype `NGRID`.
pub const DEFAULT_ALLELE_FREQ_PRIOR_POINTS: usize = 200;

/// Lower endpoint of the default carrier-frequency `q` grid (rare).
/// Prototype `QEXT` lower bound.
pub const DEFAULT_CARRIER_FREQ_LO: f64 = 0.004;

/// Upper endpoint of the default carrier-frequency `q` grid (common; cannot
/// reach fixation — a fixed duplication is caught by absolute coverage, not
/// this term). Prototype `QEXT` upper bound.
pub const DEFAULT_CARRIER_FREQ_HI: f64 = 0.6;

/// Number of points in the default carrier-frequency `q` grid. Prototype
/// `QEXT` length.
pub const DEFAULT_CARRIER_FREQ_POINTS: usize = 40;

/// VAF threshold for a "confident hom-alt" sample in the hom-alt veto (a
/// sample this saturated cannot be a `VAF = m/T < 1` paralog carrier).
/// Prototype `af > 0.9`.
pub const DEFAULT_HOMALT_VAF_THRESHOLD: f64 = 0.9;

/// Total-read-depth threshold for the hom-alt veto — below it a high VAF is
/// not "confident". Prototype `nn >= 5`.
pub const DEFAULT_HOMALT_MIN_DEPTH: u32 = 5;

/// The folded site-frequency-spectrum prior on the allele frequency `p`,
/// under which H1 marginalises `p`. The shape is fixed — proportional to
/// `1/(p(1−p))` over `p ∈ [1/2N, 1−1/2N]`, where `N` is the number of
/// samples (real variants are mostly rare, and the `1/2N` floor — a
/// singleton, the finest frequency the panel can resolve — replaces an
/// arbitrary low-frequency cutoff; spec §5.2). The only free parameter is
/// the grid resolution; the `N`-dependent bounds and weights are derived at
/// scoring time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SfsPriorSpec {
    /// Number of `p` grid points on `[1/2N, 1−1/2N]` (uniform in `p`).
    pub n_points: usize,
}

impl Default for SfsPriorSpec {
    fn default() -> Self {
        Self {
            n_points: DEFAULT_ALLELE_FREQ_PRIOR_POINTS,
        }
    }
}

/// A uniform (flat) linear grid `[lo, hi]` with `n_points` samples,
/// inclusive of both ends. Used for the carrier-frequency `q` under H2.
///
/// A well-formed grid has `lo < hi` and `n_points >= 2` (so it spans the
/// interval rather than collapsing to a point / empty set); [`new`]
/// enforces that at the construction boundary. The fields stay public
/// because the scorer (Q3) reads `lo`/`hi`/`n_points` directly to lay out
/// the linspace — but callers building a grid from untrusted input should
/// go through [`new`] rather than a raw literal.
///
/// (Unlike [`SfsPriorSpec`], which has one canonical default and so gets a
/// [`Default`] impl, `GridSpec` is a generic grid with no universal
/// default; its one named default is [`GridSpec::DEFAULT_CARRIER_FREQ`], an
/// associated `const` usable in `const` context.)
///
/// [`new`]: GridSpec::new
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridSpec {
    /// Lower endpoint (included).
    pub lo: f64,
    /// Upper endpoint (included).
    pub hi: f64,
    /// Number of grid points (`>= 2` to span the interval).
    pub n_points: usize,
}

impl GridSpec {
    /// The default carrier-frequency grid: 40 flat points on `[0.004, 0.6]`
    /// (rare to common). The `0.6` ceiling cannot reach fixation, but a
    /// fully-fixed duplication is caught by absolute coverage, not the
    /// carrier-frequency term — raising it barely moves π (spec §7).
    /// Prototype `QEXT = linspace(0.004, 0.6, 40)`.
    pub const DEFAULT_CARRIER_FREQ: Self = Self {
        lo: DEFAULT_CARRIER_FREQ_LO,
        hi: DEFAULT_CARRIER_FREQ_HI,
        n_points: DEFAULT_CARRIER_FREQ_POINTS,
    };

    /// Build a grid, returning `None` for a degenerate spec (`lo >= hi`, a
    /// non-finite endpoint, or `n_points < 2`) that would give an empty or
    /// point-collapsed integration. The trusted [`DEFAULT_CARRIER_FREQ`]
    /// bypasses this (a compile-time literal); use `new` for anything
    /// derived from input.
    ///
    /// [`DEFAULT_CARRIER_FREQ`]: GridSpec::DEFAULT_CARRIER_FREQ
    pub fn new(lo: f64, hi: f64, n_points: usize) -> Option<Self> {
        (lo.is_finite() && hi.is_finite() && lo < hi && n_points >= 2).then_some(Self {
            lo,
            hi,
            n_points,
        })
    }
}

/// The fixed model constants and integration grids for the hidden-paralog
/// likelihood ratio. Every field defaults to the value the validated
/// tomato2 prototype used; [`ParalogModelParams::default`] reproduces that
/// configuration. The scoring function (Q3) takes this by reference so it
/// depends only on its arguments (arch Premise 2).
#[derive(Debug, Clone, PartialEq)]
pub struct ParalogModelParams {
    /// VAF error floor `ε` for hom-ref / non-carrier alleles (and `1 − ε`
    /// for hom-alt). See [`DEFAULT_PSEUDOCOUNT_VAF`].
    pub pseudocount_vaf: f64,
    /// Winsor cap on `relative_copy_number` for the Normal coverage tail.
    /// See [`DEFAULT_MAX_RELATIVE_COPY_NUMBER`].
    pub max_relative_copy_number: f64,
    /// Hidden-paralog total copy numbers `T` enumerated under H2. See
    /// [`DEFAULT_CARRIER_COPY_NUMBERS`].
    pub carrier_copy_numbers: Vec<u32>,
    /// The folded-SFS prior on the allele frequency `p` marginalised under
    /// H1.
    pub allele_freq_prior: SfsPriorSpec,
    /// The flat carrier-frequency `q` grid marginalised under H2.
    pub carrier_freq_grid: GridSpec,
    /// VAF threshold for the hom-alt veto. See
    /// [`DEFAULT_HOMALT_VAF_THRESHOLD`].
    pub homalt_vaf_threshold: f64,
    /// Depth threshold for the hom-alt veto. See
    /// [`DEFAULT_HOMALT_MIN_DEPTH`].
    pub homalt_min_depth: u32,
}

impl Default for ParalogModelParams {
    fn default() -> Self {
        Self {
            pseudocount_vaf: DEFAULT_PSEUDOCOUNT_VAF,
            max_relative_copy_number: DEFAULT_MAX_RELATIVE_COPY_NUMBER,
            carrier_copy_numbers: DEFAULT_CARRIER_COPY_NUMBERS.to_vec(),
            allele_freq_prior: SfsPriorSpec::default(),
            carrier_freq_grid: GridSpec::DEFAULT_CARRIER_FREQ,
            homalt_vaf_threshold: DEFAULT_HOMALT_VAF_THRESHOLD,
            homalt_min_depth: DEFAULT_HOMALT_MIN_DEPTH,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pin the default model configuration to the validated prototype
    /// constants (`build_paralog_lr.py`): `EPS`, `T_CARRIER`, the `q` grid
    /// `linspace(0.004, 0.6, 40)`, `NGRID = 200`, the `4×` cap, and the
    /// hom-alt veto `(af > 0.9, depth >= 5)`. A drift in any of these
    /// changes the score, so it is asserted rather than left implicit.
    #[test]
    fn default_params_match_prototype_constants() {
        let p = ParalogModelParams::default();
        assert_eq!(p.pseudocount_vaf, 0.01);
        assert_eq!(p.max_relative_copy_number, 4.0);
        assert_eq!(p.carrier_copy_numbers, vec![3, 4, 6, 8]);
        assert_eq!(p.allele_freq_prior.n_points, 200);
        assert_eq!(
            p.carrier_freq_grid,
            GridSpec {
                lo: 0.004,
                hi: 0.6,
                n_points: 40,
            }
        );
        assert_eq!(p.homalt_vaf_threshold, 0.9);
        assert_eq!(p.homalt_min_depth, 5);
    }

    /// The two grid specs default independently to their documented shapes.
    #[test]
    fn grid_specs_have_expected_defaults() {
        assert_eq!(SfsPriorSpec::default().n_points, 200);
        let q = GridSpec::DEFAULT_CARRIER_FREQ;
        assert_eq!((q.lo, q.hi, q.n_points), (0.004, 0.6, 40));
        assert!(q.lo < q.hi && q.n_points >= 2);
    }

    /// The carrier copy numbers and the winsor cap are not independent: the
    /// top configuration's coverage mean `max(T)/2` must equal
    /// `max_relative_copy_number`, or the top carrier is pinned at (or past)
    /// the clip boundary (arch Premise 2). A retune that raised `T` without
    /// the cap — or vice versa — would break the coverage model Q2/Q3 rely
    /// on; this pins the relationship the per-value tests miss.
    #[test]
    fn default_carrier_copy_numbers_ascending_and_top_mean_equals_cap() {
        let p = ParalogModelParams::default();
        assert!(
            p.carrier_copy_numbers.windows(2).all(|w| w[0] < w[1]),
            "carrier copy numbers must be strictly ascending"
        );
        let top_mean = f64::from(*p.carrier_copy_numbers.iter().max().unwrap()) / 2.0;
        assert_eq!(
            top_mean, p.max_relative_copy_number,
            "top carrier mean (max T / 2) must equal the winsor cap"
        );
    }

    /// `GridSpec::new` accepts a well-formed grid and rejects every
    /// degenerate spec (reversed / collapsed bounds, non-finite endpoint,
    /// fewer than two points).
    #[test]
    fn grid_spec_new_validates_bounds() {
        assert!(GridSpec::new(0.004, 0.6, 40).is_some()); // well-formed
        assert!(GridSpec::new(0.004, 0.6, 2).is_some()); // minimal valid
        assert!(GridSpec::new(0.6, 0.004, 40).is_none()); // lo > hi
        assert!(GridSpec::new(0.004, 0.004, 40).is_none()); // lo == hi
        assert!(GridSpec::new(0.004, 0.6, 1).is_none()); // < 2 points
        assert!(GridSpec::new(f64::NAN, 0.6, 40).is_none()); // non-finite
    }
}
