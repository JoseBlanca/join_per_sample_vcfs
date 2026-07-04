//! Cohort nucleotide diversity `θ` and per-sample inbreeding `F`, estimated
//! from the rough single-sample genotype counts already stored in each `.psp`
//! sample summary — computed *before* the calling phase so it can seed the
//! SFS genotype prior's hyperparameters.
//!
//! See the spec [`doc/devel/specs/sfs_genotype_prior.md`] §5 and the
//! architecture doc [`doc/devel/architecture/sfs_genotype_prior.md`] §3. This
//! module is the pure estimator only: it reads the per-sample summaries and
//! returns the numbers; wiring it into the driver and the posterior engine is a
//! later step.
//!
//! # Why this shape
//!
//! Nucleotide diversity `θ` (≈ `π`) is the scale of the site-frequency
//! spectrum — how much variation the species carries. The SFS genotype prior
//! needs it as a hyperparameter, and needs it *before* genotyping starts, so it
//! cannot come out of the calling EM. Fortunately the per-sample pileup already
//! records the sufficient statistics in the `.psp` metadata: the number of
//! confidently-called heterozygous and homozygous-alt sites, and the count of
//! callable (confidently-covered) positions. This estimator is arithmetic on
//! those — no new pileup work, no extra pass over the data.
//!
//! ## The estimator counts allele *copies*, not genotypes
//!
//! The naive per-individual heterozygosity `n_het / callable` estimates
//! `θ·(1 − F)` — it is confounded by inbreeding, which is fatal for a selfing
//! crop line (almost no heterozygotes, yet the species may be diverse). Counting
//! alt-allele *copies* instead removes the confounding: a heterozygote carries
//! one alt copy, a homozygous-alt carries two, over `2 × callable` chromosomes
//! surveyed, so the individual's mating system cancels:
//!
//! ```text
//! θ̂ = Σ_s (n_het_s + 2·n_hom_alt_s)  /  (2 · Σ_s callable_s)
//! ```
//!
//! This is the standard reason `π`/Watterson estimators are inbreeding-robust
//! while raw heterozygosity is not. The same two counts also yield the
//! per-sample inbreeding coefficient `F = (2 − R)/(2 + R)` with `R = n_het /
//! n_hom_alt` (a het:hom-alt ratio of 2:1 is outbred `F = 0`; all-homozygous is
//! `F → 1`).
//!
//! ## Deliberate choices where the spec left a gap
//!
//! - **Ambiguous sites are excluded.** The summary also carries
//!   `n_ambiguous_sites` (variant sites where het vs hom-alt could not be
//!   resolved, common at low coverage). The copy count uses only the confident
//!   het/hom-alt calls, so ambiguous sites contribute no copies. This
//!   under-counts at low coverage and is the concrete mechanism of the
//!   "`θ̂` biases down at low coverage" limitation the spec already records; a
//!   coverage-based correction is a deferred tier.
//! - **Pooled, not a per-sample mean.** The cohort estimate is total copies over
//!   total chromosomes (callable-weighted), which down-weights thin, noisy
//!   samples — the correct cohort diversity estimator. The spec's "mean over the
//!   cohort" is realised as this pooled form.
//!
//! `θ` only needs to be right to an order of magnitude — it calibrates the
//! false-positive/recall boundary of the prior, it does not decide clear calls —
//! so a rough estimate in the right decade is enough.

use crate::paralog::inbreeding::MAX_INBREEDING_COEFFICIENT;
use crate::sample_summary::{CoverageByGcHistogram, HetCounts, SampleSummary};

/// Chromosomes per diploid position — the `2` in the copy-count denominator.
/// Named so the diploid-ploidy assumption is explicit and greppable if the
/// caller ever supports other ploidies.
const CHROMOSOMES_PER_DIPLOID: u64 = 2;

/// Alt-allele copies carried by a homozygous-alt genotype (both chromosomes).
/// A distinct constant from [`CHROMOSOMES_PER_DIPLOID`] despite the equal value:
/// this is "copies per hom-alt genotype", that is "chromosomes per position".
const ALT_COPIES_PER_HOM_ALT: u64 = 2;

/// Default species-range fallback diversity `θ`, used when the cohort is too
/// thin to estimate (see [`MIN_COPIES_FOR_ESTIMATE`]). A conservative,
/// low-diversity value (~human nucleotide diversity); a user working with a more
/// diverse organism should raise it. This is the weakly-informative prior the
/// estimator falls back toward, not a hard default for the estimate itself.
pub const DEFAULT_DIVERSITY_PRIOR: f64 = 1e-3;

/// Minimum total confident alt-allele copies across the whole cohort for the
/// diversity estimate to be trusted. Below this the sampling variance is too
/// high — the relative standard error of a copy-count rate is roughly
/// `1/√copies`, so ~30 copies is ~18 %, adequate for an order-of-magnitude `θ`
/// but the floor below which we defer to the prior instead of reporting noise.
pub const MIN_COPIES_FOR_ESTIMATE: u64 = 30;

/// How the [`DiversityEstimate::nucleotide_diversity`] value was obtained.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiversitySource {
    /// Estimated from the cohort's `.psp` sample summaries.
    Estimated,
    /// Supplied on the command line, overriding any estimate.
    CliOverride,
    /// The cohort was too thin to estimate (fewer than
    /// [`MIN_COPIES_FOR_ESTIMATE`] confident alt copies, or no callable
    /// positions); fell back to the species-range prior value.
    PriorFallback,
}

/// Cohort nucleotide diversity plus per-sample inbreeding, ready to seed the
/// SFS genotype prior. Built by [`DiversityEstimate::from_summaries`] from the
/// per-sample `.psp` summaries before the calling phase.
#[derive(Debug, Clone, PartialEq)]
pub struct DiversityEstimate {
    /// Nucleotide diversity `θ` (≈ `π`), the scale of the site-frequency
    /// spectrum. Finite and `>= 0`: the estimated path is provably non-negative,
    /// and the `CliOverride` / `PriorFallback` paths carry it through only under
    /// the caller preconditions documented on
    /// [`from_summaries`](Self::from_summaries).
    pub nucleotide_diversity: f64,
    /// Per-sample inbreeding coefficient `F ∈ [0, MAX_INBREEDING_COEFFICIENT]`,
    /// one entry per input summary in the same order, from the confident
    /// het:hom-alt ratio. Empty when `summaries` is empty.
    pub inbreeding_coefficients: Vec<f64>,
    /// Provenance of `nucleotide_diversity`, recorded in the VCF header
    /// downstream.
    pub source: DiversitySource,
}

impl DiversityEstimate {
    /// Estimate cohort `θ` and per-sample `F` from the `.psp` sample summaries.
    ///
    /// - `cli_override`: `Some(θ)` forces that `θ` (source [`CliOverride`]);
    ///   per-sample `F` is still estimated from the counts.
    /// - `prior_theta`: the species-range fallback used when the cohort carries
    ///   fewer than [`MIN_COPIES_FOR_ESTIMATE`] confident alt copies or has no
    ///   callable positions (source [`PriorFallback`]).
    ///
    /// An empty `summaries` slice yields the fallback `θ` with an empty
    /// per-sample `F` vector.
    ///
    /// # Preconditions
    ///
    /// `prior_theta` and any `cli_override` must be finite and `>= 0` — they are
    /// carried straight into [`nucleotide_diversity`](Self::nucleotide_diversity),
    /// so a `NaN`/negative value would silently poison the downstream prior.
    /// They originate from a compile-time species-range constant and a
    /// CLI-validated flag; the debug assertions here catch a wiring mistake in
    /// development without turning this pure estimator fallible.
    ///
    /// [`CliOverride`]: DiversitySource::CliOverride
    /// [`PriorFallback`]: DiversitySource::PriorFallback
    pub fn from_summaries(
        summaries: &[SampleSummary],
        prior_theta: f64,
        cli_override: Option<f64>,
    ) -> Self {
        debug_assert!(
            prior_theta.is_finite() && prior_theta >= 0.0,
            "prior_theta must be finite and non-negative, got {prior_theta}"
        );

        let inbreeding_coefficients: Vec<f64> = summaries
            .iter()
            .map(|s| estimate_inbreeding(&s.heterozygosity))
            .collect();

        if let Some(theta) = cli_override {
            debug_assert!(
                theta.is_finite() && theta >= 0.0,
                "cli_override theta must be finite and non-negative, got {theta}"
            );
            return Self {
                nucleotide_diversity: theta,
                inbreeding_coefficients,
                source: DiversitySource::CliOverride,
            };
        }

        // Pooled copy-count estimate: total alt-allele copies over total
        // chromosomes surveyed. Saturating adds match the summary accumulators'
        // convention and cannot realistically overflow (genome-scale counts).
        let mut total_copies: u64 = 0;
        let mut total_callable: u64 = 0;
        for summary in summaries {
            total_copies = total_copies.saturating_add(alt_allele_copies(&summary.heterozygosity));
            total_callable =
                total_callable.saturating_add(callable_position_count(&summary.coverage_by_gc));
        }

        if total_copies < MIN_COPIES_FOR_ESTIMATE || total_callable == 0 {
            return Self {
                nucleotide_diversity: prior_theta,
                inbreeding_coefficients,
                source: DiversitySource::PriorFallback,
            };
        }

        // Two chromosomes per diploid callable position.
        let nucleotide_diversity =
            total_copies as f64 / (CHROMOSOMES_PER_DIPLOID as f64 * total_callable as f64);
        Self {
            nucleotide_diversity,
            inbreeding_coefficients,
            source: DiversitySource::Estimated,
        }
    }
}

/// Alt-allele copies a sample contributes to the diversity estimate: one per
/// confident heterozygote, two per confident homozygous-alt (ambiguous sites
/// carry none — see the module docs). Exhaustively destructures `HetCounts` so
/// that adding a new copy-bearing count field there is a compile error here,
/// where the decision to include it must be made, rather than a silent bias.
fn alt_allele_copies(het: &HetCounts) -> u64 {
    let HetCounts {
        n_het_sites,
        n_hom_alt_sites,
        n_ambiguous_sites: _,
        n_variant_sites: _,
        min_depth: _,
        error_rate: _,
        lr_margin: _,
    } = het;
    n_het_sites.saturating_add(n_hom_alt_sites.saturating_mul(ALT_COPIES_PER_HOM_ALT))
}

/// Callable (confidently-covered) positions — the chromosome-count denominator.
/// Exhaustively destructures `CoverageByGcHistogram` so a renamed or split
/// callable-position field forces a decision here instead of silently changing
/// the denominator.
fn callable_position_count(coverage: &CoverageByGcHistogram) -> u64 {
    let CoverageByGcHistogram {
        callable_positions,
        window_bp: _,
        gc_bins: _,
        depth_bin_width: _,
        depth_bins: _,
        n_positions: _,
        n_skipped_tiles: _,
        counts: _,
    } = coverage;
    *callable_positions
}

/// Per-sample inbreeding coefficient from the confident het:hom-alt counts:
/// `F = (2 − R)/(2 + R)` with `R = n_het / n_hom_alt`, rearranged to
/// `F = (2·n_hom_alt − n_het) / (2·n_hom_alt + n_het)` so the `n_hom_alt = 0`
/// (`R → ∞`) case needs no special guard. Clamped to
/// `[0, MAX_INBREEDING_COEFFICIENT]`; a sample with no confident variant sites
/// is inbreeding-unidentifiable and reported as outbred (`0.0`). Exhaustively
/// destructures `HetCounts` for the same field-safety reason as
/// [`alt_allele_copies`].
fn estimate_inbreeding(het: &HetCounts) -> f64 {
    let HetCounts {
        n_het_sites,
        n_hom_alt_sites,
        n_ambiguous_sites: _,
        n_variant_sites: _,
        min_depth: _,
        error_rate: _,
        lr_margin: _,
    } = het;
    let n_het = *n_het_sites as f64;
    let n_hom_alt = *n_hom_alt_sites as f64;
    let denominator = 2.0 * n_hom_alt + n_het;
    if denominator == 0.0 {
        return 0.0;
    }
    // Reuse the paralog filter's inbreeding ceiling so the two F clamps in the
    // codebase (paralog scorer, this estimator) share one source of truth and
    // cannot drift — same `[0, 0.99]` inbreeding range, same quantity.
    ((2.0 * n_hom_alt - n_het) / denominator).clamp(0.0, MAX_INBREEDING_COEFFICIENT)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sample_summary::{CoverageByGcHistogram, SAMPLE_SUMMARY_VERSION};

    /// Build a minimal `SampleSummary` carrying just the fields the diversity
    /// estimator reads: the confident het/hom-alt counts and the callable
    /// total. The histogram matrix and unrelated het-model knobs are set to
    /// inert values.
    fn summary(n_het: u64, n_hom_alt: u64, callable: u64) -> SampleSummary {
        SampleSummary {
            version: SAMPLE_SUMMARY_VERSION,
            coverage_by_gc: CoverageByGcHistogram {
                window_bp: 500,
                gc_bins: 1,
                depth_bin_width: 1.0,
                depth_bins: 1,
                n_positions: callable,
                n_skipped_tiles: 0,
                callable_positions: callable,
                counts: vec![0, 0],
            },
            heterozygosity: HetCounts {
                n_het_sites: n_het,
                n_hom_alt_sites: n_hom_alt,
                n_ambiguous_sites: 0,
                n_variant_sites: n_het + n_hom_alt,
                min_depth: 1,
                error_rate: 0.01,
                lr_margin: std::f64::consts::LN_10,
            },
        }
    }

    /// The estimate counts alt-allele copies, so it is invariant to how those
    /// copies are packaged into genotypes — shifting all copies from
    /// heterozygous into homozygous-alt (i.e. changing `F`) at the same copy
    /// total and callable count leaves `θ̂` unchanged.
    #[test]
    fn theta_is_inbreeding_free() {
        // 200 alt copies either way: all-het (200 het) vs all-hom (100 hom-alt).
        let outbred = DiversityEstimate::from_summaries(&[summary(200, 0, 1_000_000)], 1e-3, None);
        let selfer = DiversityEstimate::from_summaries(&[summary(0, 100, 1_000_000)], 1e-3, None);
        assert_eq!(outbred.source, DiversitySource::Estimated);
        assert_eq!(selfer.source, DiversitySource::Estimated);
        assert!((outbred.nucleotide_diversity - selfer.nucleotide_diversity).abs() < 1e-12);
        assert!((outbred.nucleotide_diversity - 200.0 / 2_000_000.0).abs() < 1e-12);
        // ...but F does differ: outbred (only hets) vs a pure selfer (only homs).
        assert_eq!(outbred.inbreeding_coefficients[0], 0.0);
        assert_eq!(
            selfer.inbreeding_coefficients[0],
            MAX_INBREEDING_COEFFICIENT
        );
    }

    /// A human-scale count (θ ≈ 1e-3): 1000 het + 500 hom-alt = 2000 copies over
    /// 1e6 callable → 2000 / (2·1e6) = 1e-3.
    #[test]
    fn theta_matches_human_scale() {
        let est = DiversityEstimate::from_summaries(&[summary(1000, 500, 1_000_000)], 1e-4, None);
        assert_eq!(est.source, DiversitySource::Estimated);
        assert!(
            (est.nucleotide_diversity - 1e-3).abs() < 1e-9,
            "theta = {}",
            est.nucleotide_diversity
        );
    }

    /// Pooled (callable-weighted), not an unweighted mean of per-sample rates:
    /// a high-copy/high-callable sample and a thin sample combine by totals.
    #[test]
    fn theta_is_pooled_over_the_cohort() {
        // Sample 1: 2000 copies / 1e6 callable (rate 1e-3).
        // Sample 2:   40 copies / 1e5 callable (rate 2e-4).
        // Pooled: 2040 / (2·1.1e6) ≈ 9.27e-4, NOT the mean of the two rates (6e-4).
        let est = DiversityEstimate::from_summaries(
            &[summary(2000, 0, 1_000_000), summary(40, 0, 100_000)],
            1e-4,
            None,
        );
        let pooled = 2040.0 / (2.0 * 1_100_000.0);
        assert!((est.nucleotide_diversity - pooled).abs() < 1e-12);
        let mean_of_rates = (1e-3 + 2e-4) / 2.0;
        assert!((est.nucleotide_diversity - mean_of_rates).abs() > 1e-5);
    }

    /// Exactly `MIN_COPIES_FOR_ESTIMATE` copies estimates; one below falls back —
    /// pins the `<` boundary against an off-by-one flip.
    #[test]
    fn estimates_at_exactly_min_copies() {
        // 30 copies exactly (0 het + 15 hom-alt = 30), callable > 0.
        let at = DiversityEstimate::from_summaries(&[summary(0, 15, 1_000_000)], 1e-3, None);
        assert_eq!(at.source, DiversitySource::Estimated);
        // 29 copies (29 het) is one below the floor → fallback.
        let below = DiversityEstimate::from_summaries(&[summary(29, 0, 1_000_000)], 1e-3, None);
        assert_eq!(below.source, DiversitySource::PriorFallback);
    }

    /// Too few confident copies → defer to the prior rather than report noise.
    #[test]
    fn thin_cohort_falls_back_to_prior() {
        // 10 het + 5 hom-alt = 20 copies < MIN_COPIES_FOR_ESTIMATE (30).
        let est = DiversityEstimate::from_summaries(&[summary(10, 5, 1_000_000)], 5e-3, None);
        assert_eq!(est.source, DiversitySource::PriorFallback);
        assert_eq!(est.nucleotide_diversity, 5e-3);
        // F is still estimated from the (thin) counts.
        assert_eq!(est.inbreeding_coefficients.len(), 1);
    }

    /// No callable positions → fallback (avoids a divide-by-zero rate).
    #[test]
    fn zero_callable_falls_back() {
        let est = DiversityEstimate::from_summaries(&[summary(500, 500, 0)], 2e-3, None);
        assert_eq!(est.source, DiversitySource::PriorFallback);
        assert_eq!(est.nucleotide_diversity, 2e-3);
    }

    /// Empty cohort → fallback θ, empty per-sample F.
    #[test]
    fn empty_cohort_falls_back() {
        let est = DiversityEstimate::from_summaries(&[], 1e-3, None);
        assert_eq!(est.source, DiversitySource::PriorFallback);
        assert_eq!(est.nucleotide_diversity, 1e-3);
        assert!(est.inbreeding_coefficients.is_empty());
    }

    /// A CLI override wins over the estimate, but per-sample F is still computed.
    #[test]
    fn cli_override_wins_but_f_still_estimated() {
        let est =
            DiversityEstimate::from_summaries(&[summary(1000, 500, 1_000_000)], 1e-4, Some(0.02));
        assert_eq!(est.source, DiversitySource::CliOverride);
        assert_eq!(est.nucleotide_diversity, 0.02);
        assert_eq!(est.inbreeding_coefficients.len(), 1);
    }

    /// `nucleotide_diversity` is finite and non-negative on every source, given
    /// the documented preconditions (release build — debug builds assert first).
    #[cfg(not(debug_assertions))]
    #[test]
    fn theta_is_finite_and_nonnegative_for_all_sources() {
        for est in [
            DiversityEstimate::from_summaries(&[summary(1000, 500, 1_000_000)], 1e-3, None),
            DiversityEstimate::from_summaries(&[summary(1000, 500, 1_000_000)], 1e-3, Some(0.02)),
            DiversityEstimate::from_summaries(&[summary(1, 0, 1_000_000)], 1e-3, None),
        ] {
            assert!(est.nucleotide_diversity.is_finite() && est.nucleotide_diversity >= 0.0);
        }
    }

    /// Per-sample F stays aligned with its summary's position — a swap in the
    /// map/collect would show up here (existing F tests are single-sample only).
    #[test]
    fn preserves_inbreeding_order_across_samples() {
        // Sample 0: pure selfer (F ceiling). Sample 1: outbred 2:1 (F = 0).
        let est = DiversityEstimate::from_summaries(
            &[summary(0, 100, 1_000_000), summary(200, 100, 1_000_000)],
            1e-3,
            None,
        );
        assert_eq!(est.inbreeding_coefficients.len(), 2);
        assert_eq!(est.inbreeding_coefficients[0], MAX_INBREEDING_COEFFICIENT);
        assert_eq!(est.inbreeding_coefficients[1], 0.0);
    }

    /// F recovery from the het:hom-alt ratio across the mating-system range.
    #[test]
    fn inbreeding_from_het_hom_ratio() {
        // Outbred 2:1 het:hom → F = 0.
        assert_eq!(
            estimate_inbreeding(&summary(200, 100, 0).heterozygosity),
            0.0
        );
        // Pure selfer, only hom-alt → F clamped to the ceiling.
        assert_eq!(
            estimate_inbreeding(&summary(0, 100, 0).heterozygosity),
            MAX_INBREEDING_COEFFICIENT
        );
        // Excess heterozygosity (more het than 2:1) clamps at 0, never negative.
        assert_eq!(
            estimate_inbreeding(&summary(400, 100, 0).heterozygosity),
            0.0
        );
        // Intermediate: 1:1 het:hom → R = 1 → F = (2−1)/(2+1) = 1/3.
        let f = estimate_inbreeding(&summary(100, 100, 0).heterozygosity);
        assert!((f - 1.0 / 3.0).abs() < 1e-12, "F = {f}");
        // No confident variant sites → outbred by default.
        assert_eq!(estimate_inbreeding(&summary(0, 0, 0).heterozygosity), 0.0);
    }
}
