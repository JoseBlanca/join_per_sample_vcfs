//! The SFS-marginalized Hardy–Weinberg genotype prior.
//!
//! The fix for the single-sample low-coverage het over-call (spec
//! [`doc/devel/specs/sfs_genotype_prior.md`], arch
//! [`doc/devel/architecture/sfs_genotype_prior.md`] §4). Instead of estimating a
//! per-site allele frequency and plugging it into Hardy–Weinberg — which, for a
//! single sample, pins the frequency on a reference-heavy pseudocount and skews
//! ~22:1 toward heterozygotes — this **averages** the genotype prior over the
//! whole site-frequency spectrum, so no fabricated point estimate can push it.
//!
//! # What it computes
//!
//! For one diploid sample at a biallelic site, the prior probability of each
//! genotype is the Wright inbreeding-adjusted Hardy–Weinberg probability
//! `P(g | p, F)` averaged over a prior on the alt frequency `p`:
//!
//! ```text
//! prior(g) = invariant_site_mass·[g = hom-ref]  +  Σ_i  (θ / p_i) · P(g | p_i, F)
//! ```
//!
//! The frequency prior is a mixture: a point mass at `p = 0` (the site is
//! **invariant** — the whole cohort is reference) plus the neutral
//! mutation–drift spectrum `θ / p` over the segregating grid. Two properties
//! matter and both are proven as tests:
//!
//! - **The het:hom-alt ratio is `2(1−F)/(1+F)`, independent of `θ`.** On the
//!   symmetric grid the segregating sums give exactly this Wright ratio, so the
//!   genotype call is fixed regardless of the species' diversity — at `F = 0` it
//!   is the defensible 2:1 (mildly pro-het), a factor of ~10 less than the old
//!   plug-in prior.
//! - **The invariant mass only moves the hom-ref weight.** It is the
//!   false-positive/recall knob (how readily a weak-evidence site is called a
//!   variant at all), cleanly separated from the het-vs-hom-alt decision.
//!
//! # Scope (first cut)
//!
//! Biallelic diploid, "mode 1" of the architecture: the frequency posterior is
//! taken to be the SFS prior itself (the single-sample / small-cohort regime),
//! so the genotype prior is a fixed table per `(θ, invariant_site_mass, F)`. The
//! large-cohort mode (fold cohort read evidence into the frequency posterior so
//! the prior sharpens toward `HWE(p̂)`) and general `(ploidy, n_alleles)` shapes
//! are deferred (arch §7).

use crate::genetics::{sfs_grid_point, wright_genotype_log_priors};

/// Default number of frequency grid points for the marginalization — the same
/// 200 points the paralog filter uses for its allele-frequency prior. This
/// resolves the `θ/p` spectrum finely enough that the genotype prior is stable
/// to well under the precision the downstream call needs.
pub const DEFAULT_SFS_GRID_POINTS: usize = 200;

/// Default invariant-site point mass — **provisional**, pending the calibration
/// sweep in the validation step. A reference weight of `1.0` against the
/// `θ`-scaled segregating sums makes the prior odds of a site being variant at
/// all scale with `θ` (diverse species → more readily variant), which is the
/// intended shape; the magnitude is the false-positive/recall knob and does not
/// affect the het:hom-alt ratio.
pub const DEFAULT_INVARIANT_SITE_MASS: f64 = 1.0;

/// Genotype index into the biallelic-diploid log-prior triple returned by
/// [`SfsGenotypePrior::biallelic_diploid_log_priors`].
const HOM_REF: usize = 0;
const HET: usize = 1;
const HOM_ALT: usize = 2;

/// A folded site-frequency-spectrum frequency grid: the evaluation points `p_i`
/// in `(0, 1)` and their neutral mutation–drift prior weights `θ / p_i`.
///
/// The points are uniform on `[inset, 1 − inset]` with `inset = 1/(2n)`, so the
/// grid is symmetric about `0.5` and never lands on the degenerate frequencies
/// `0` or `1`. Internal to [`SfsGenotypePrior`].
#[derive(Debug, Clone)]
struct FrequencyGrid {
    /// Frequency evaluation points `p_i ∈ (0, 1)`.
    points: Vec<f64>,
    /// Neutral SFS prior weight `θ / p_i` at each point.
    weights: Vec<f64>,
}

impl FrequencyGrid {
    /// Build the folded neutral-SFS grid of `n_points` points with diversity
    /// `nucleotide_diversity` (`θ`). `n_points == 0` yields an empty grid (the
    /// genotype prior then reduces to the invariant mass alone — hom-ref).
    fn folded_sfs(n_points: usize, nucleotide_diversity: f64) -> Self {
        if n_points == 0 {
            return Self {
                points: Vec::new(),
                weights: Vec::new(),
            };
        }
        let inset = 1.0 / (2.0 * n_points as f64);
        let points: Vec<f64> = (0..n_points)
            .map(|i| sfs_grid_point(i, n_points, inset))
            .collect();
        let weights: Vec<f64> = points.iter().map(|&p| nucleotide_diversity / p).collect();
        Self { points, weights }
    }
}

/// The SFS-marginalized genotype prior. Holds the cohort-constant pieces (the
/// diversity-scaled frequency grid and the invariant-site mass); the per-sample
/// inbreeding `F` is supplied per call, since it varies by sample.
#[derive(Debug, Clone)]
pub struct SfsGenotypePrior {
    grid: FrequencyGrid,
    invariant_site_mass: f64,
}

impl SfsGenotypePrior {
    /// Build the prior for cohort diversity `nucleotide_diversity` (`θ`),
    /// invariant-site mass `invariant_site_mass`, over `n_points` frequency grid
    /// points.
    ///
    /// # Preconditions
    ///
    /// `nucleotide_diversity` and `invariant_site_mass` must be finite and
    /// `>= 0`. They originate from [`DiversityEstimate`] (θ is non-negative by
    /// construction) and a compile-time constant, so the debug assertions here
    /// catch a wiring mistake in development without making this pure prior
    /// fallible; a negative θ would otherwise drive genotype probabilities
    /// negative and yield `NaN` log-priors.
    ///
    /// [`DiversityEstimate`]: crate::var_calling::diversity::DiversityEstimate
    pub fn new(nucleotide_diversity: f64, invariant_site_mass: f64, n_points: usize) -> Self {
        debug_assert!(
            nucleotide_diversity.is_finite() && nucleotide_diversity >= 0.0,
            "nucleotide_diversity must be finite and non-negative, got {nucleotide_diversity}"
        );
        debug_assert!(
            invariant_site_mass.is_finite() && invariant_site_mass >= 0.0,
            "invariant_site_mass must be finite and non-negative, got {invariant_site_mass}"
        );
        Self {
            grid: FrequencyGrid::folded_sfs(n_points, nucleotide_diversity),
            invariant_site_mass,
        }
    }

    /// The biallelic-diploid genotype **log**-priors `[ln P(0/0), ln P(0/1),
    /// ln P(1/1)]` at inbreeding coefficient `inbreeding_coefficient` (`F`),
    /// normalised to sum to one.
    ///
    /// Marginalises the Wright HWE probabilities over the frequency grid and
    /// adds the invariant mass to the hom-ref term. A grid with no segregating
    /// mass (empty grid, or `θ = 0`) leaves only the invariant mass, i.e. the
    /// prior collapses to certain hom-ref.
    ///
    /// # Preconditions
    ///
    /// `inbreeding_coefficient` must be finite and in `[0, 1]` (checked with a
    /// debug assertion; it comes from `DiversityEstimate`, which clamps it).
    pub fn biallelic_diploid_log_priors(&self, inbreeding_coefficient: f64) -> [f64; 3] {
        debug_assert!(
            inbreeding_coefficient.is_finite() && (0.0..=1.0).contains(&inbreeding_coefficient),
            "inbreeding_coefficient must be finite and in [0, 1], got {inbreeding_coefficient}"
        );

        // Linear accumulation of Σ_i w_i · P(g | p_i, F); with valid inputs all
        // terms are non-negative and there is no underflow concern at these
        // scales.
        let mut acc = [0.0_f64; 3];
        for (&p, &w) in self.grid.points.iter().zip(&self.grid.weights) {
            let (log_hom_ref, log_het, log_hom_alt) =
                wright_genotype_log_priors(p, inbreeding_coefficient);
            acc[HOM_REF] += w * log_hom_ref.exp();
            acc[HET] += w * log_het.exp();
            acc[HOM_ALT] += w * log_hom_alt.exp();
        }
        acc[HOM_REF] += self.invariant_site_mass;

        let total = acc[HOM_REF] + acc[HET] + acc[HOM_ALT];
        // With valid inputs `total >= invariant_site_mass >= 0`. `!(total > 0.0)`
        // (rather than `total <= 0.0`) also catches a non-finite `total`, so the
        // degenerate "no segregating mass" case — and any slipped-through
        // invalid input in a release build — collapses to certain hom-ref
        // instead of propagating `NaN`.
        if !(total > 0.0) {
            return [0.0, f64::NEG_INFINITY, f64::NEG_INFINITY];
        }
        [
            (acc[HOM_REF] / total).ln(),
            (acc[HET] / total).ln(),
            (acc[HOM_ALT] / total).ln(),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn probs(prior: &SfsGenotypePrior, f: f64) -> [f64; 3] {
        let logs = prior.biallelic_diploid_log_priors(f);
        [logs[0].exp(), logs[1].exp(), logs[2].exp()]
    }

    /// The log-priors are a proper distribution: they exponentiate to a vector
    /// summing to one.
    #[test]
    fn log_priors_are_normalised() {
        let prior =
            SfsGenotypePrior::new(1e-3, DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS);
        for &f in &[0.0, 0.3, 0.9] {
            let s: f64 = probs(&prior, f).iter().sum();
            assert!((s - 1.0).abs() < 1e-12, "F={f} sum={s}");
        }
    }

    /// At no inbreeding the het:hom-alt ratio is exactly 2:1 — the defensible
    /// SFS-marginal answer (a symmetric grid makes Σ(1−p) = Σp, so the Wright
    /// segregating sums give 2:1 to floating-point precision).
    #[test]
    fn het_homalt_ratio_is_two_to_one_at_no_inbreeding() {
        let prior =
            SfsGenotypePrior::new(1e-3, DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS);
        let p = probs(&prior, 0.0);
        assert!(
            (p[HET] / p[HOM_ALT] - 2.0).abs() < 1e-9,
            "ratio={}",
            p[HET] / p[HOM_ALT]
        );
    }

    /// With inbreeding the ratio follows Wright: het:hom-alt = 2(1−F)/(1+F).
    #[test]
    fn het_homalt_ratio_follows_wright_with_inbreeding() {
        let prior =
            SfsGenotypePrior::new(1e-3, DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS);
        for &f in &[0.0, 0.3, 0.6, 0.9] {
            let p = probs(&prior, f);
            let expected = 2.0 * (1.0 - f) / (1.0 + f);
            assert!(
                (p[HET] / p[HOM_ALT] - expected).abs() < 1e-9,
                "F={f} got {} want {expected}",
                p[HET] / p[HOM_ALT]
            );
        }
    }

    /// The het:hom-alt ratio does not depend on the diversity θ — θ only sets
    /// the variant-vs-invariant balance, not the genotype split.
    #[test]
    fn ratio_is_diversity_independent() {
        let low = SfsGenotypePrior::new(1e-4, DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS);
        let high =
            SfsGenotypePrior::new(2e-2, DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS);
        let rl = probs(&low, 0.0);
        let rh = probs(&high, 0.0);
        assert!((rl[HET] / rl[HOM_ALT] - rh[HET] / rh[HOM_ALT]).abs() < 1e-9);
    }

    /// A larger invariant mass raises the hom-ref probability but leaves the
    /// het:hom-alt ratio untouched — the two are decoupled knobs.
    #[test]
    fn invariant_mass_raises_homref_not_the_ratio() {
        let small = SfsGenotypePrior::new(1e-3, 1.0, DEFAULT_SFS_GRID_POINTS);
        let large = SfsGenotypePrior::new(1e-3, 100.0, DEFAULT_SFS_GRID_POINTS);
        let ps = probs(&small, 0.0);
        let pl = probs(&large, 0.0);
        assert!(
            pl[HOM_REF] > ps[HOM_REF],
            "{} !> {}",
            pl[HOM_REF],
            ps[HOM_REF]
        );
        // ratio unchanged
        assert!((ps[HET] / ps[HOM_ALT] - pl[HET] / pl[HOM_ALT]).abs() < 1e-9);
    }

    /// Inbreeding moves mass off the heterozygote and onto the homozygotes.
    #[test]
    fn inbreeding_shifts_mass_to_homozygotes() {
        let prior =
            SfsGenotypePrior::new(1e-3, DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS);
        let outbred = probs(&prior, 0.0);
        let selfer = probs(&prior, 0.9);
        assert!(
            selfer[HET] < outbred[HET],
            "het did not fall under inbreeding"
        );
        assert!(selfer[HOM_ALT] > outbred[HOM_ALT], "hom-alt did not rise");
    }

    /// Reproduces the validated prototype at the human-diversity default:
    /// θ = 1e-3, mass = 1, F = 0 gives hom-ref ≈ 0.83, het ≈ 0.11, and the 2:1
    /// het:hom-alt split. Bounds are loose (the exact triple is grid-dependent);
    /// the ratio is the tight, load-bearing assertion above.
    #[test]
    fn matches_prototype_at_human_default() {
        let prior =
            SfsGenotypePrior::new(1e-3, DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS);
        let p = probs(&prior, 0.0);
        assert!((0.78..0.88).contains(&p[HOM_REF]), "hom-ref={}", p[HOM_REF]);
        assert!((0.08..0.15).contains(&p[HET]), "het={}", p[HET]);
    }

    /// Zero diversity → no segregating mass → certain hom-ref.
    #[test]
    fn zero_diversity_is_certain_homref() {
        let prior =
            SfsGenotypePrior::new(0.0, DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS);
        let logs = prior.biallelic_diploid_log_priors(0.0);
        assert_eq!(logs[HOM_REF], 0.0);
        assert_eq!(logs[HET], f64::NEG_INFINITY);
        assert_eq!(logs[HOM_ALT], f64::NEG_INFINITY);
    }

    /// Zero invariant mass with a real grid is still a valid distribution — the
    /// prior is then purely segregating (no monomorphic-site mass), so hom-ref
    /// is smaller but the 2:1 het:hom-alt ratio still holds and all logs finite.
    #[test]
    fn zero_invariant_mass_still_valid() {
        let prior = SfsGenotypePrior::new(1e-3, 0.0, DEFAULT_SFS_GRID_POINTS);
        let logs = prior.biallelic_diploid_log_priors(0.0);
        assert!(logs.iter().all(|l| l.is_finite()));
        let p = probs(&prior, 0.0);
        assert!((p.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!((p[HET] / p[HOM_ALT] - 2.0).abs() < 1e-9);
    }

    /// No segregating mass AND no invariant mass (empty grid, zero mass) hits the
    /// `total <= 0` guard and collapses to certain hom-ref rather than `NaN`.
    #[test]
    fn total_zero_collapses_to_homref() {
        let prior = SfsGenotypePrior::new(1e-3, 0.0, 0);
        let logs = prior.biallelic_diploid_log_priors(0.0);
        assert_eq!(logs[HOM_REF], 0.0);
        assert_eq!(logs[HET], f64::NEG_INFINITY);
        assert_eq!(logs[HOM_ALT], f64::NEG_INFINITY);
    }

    /// Full inbreeding (`F = 1`) floors the heterozygote: its probability is the
    /// probability floor, so its log-prior is finite and very negative, and the
    /// homozygotes carry ~all the mass.
    #[test]
    fn full_inbreeding_floors_het() {
        let prior =
            SfsGenotypePrior::new(1e-3, DEFAULT_INVARIANT_SITE_MASS, DEFAULT_SFS_GRID_POINTS);
        let logs = prior.biallelic_diploid_log_priors(1.0);
        assert!(logs.iter().all(|l| l.is_finite()), "logs={logs:?}");
        assert!(
            logs[HET] < -100.0,
            "het log-prior not floored: {}",
            logs[HET]
        );
        let p = probs(&prior, 1.0);
        assert!(p[HOM_REF] + p[HOM_ALT] > 0.999);
    }

    /// The folded-SFS grid is symmetric about 0.5 (Σp = Σ(1−p), the fact the
    /// exact 2:1 ratio rests on) and carries the neutral `θ/p` weights.
    #[test]
    fn folded_sfs_grid_is_symmetric_with_theta_over_p_weights() {
        let theta = 1e-3;
        let grid = FrequencyGrid::folded_sfs(200, theta);
        let sum_p: f64 = grid.points.iter().sum();
        let sum_1mp: f64 = grid.points.iter().map(|&p| 1.0 - p).sum();
        assert!((sum_p - sum_1mp).abs() < 1e-9, "grid not symmetric");
        for (&p, &w) in grid.points.iter().zip(&grid.weights) {
            assert!((w - theta / p).abs() < 1e-18);
            assert!(p > 0.0 && p < 1.0);
        }
    }
}
