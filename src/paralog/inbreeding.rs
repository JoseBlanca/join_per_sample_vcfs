//! The per-sample inbreeding coefficient `F` (Q4).
//!
//! The Wright genotype prior in the H1/H2 score
//! ([`super::score_locus_for_paralogy`]) needs each sample's inbreeding
//! coefficient `F` (selfer ↔ outbred). It is formed from two ingredients:
//!
//! - **`obs_het`** — the sample's observed heterozygosity as a *rate*,
//!   `n_het_sites / callable_positions`, a single-sample observable from the
//!   Stage-1 pileup. It must be the het *rate*, **not**
//!   `n_het/(n_het+n_hom_alt)` (het over *variant* sites), which is dominated
//!   by the hom-alt count — i.e. by divergence from the reference — and
//!   *inverts* (Spearman −0.28 vs the true cohort `F` on tomato2; the het rate
//!   gives +0.86). Spec §3, arch Premise 3.
//! - **`Hexp`** — one cohort-wide expected heterozygosity scalar (`mean 2pq`),
//!   accumulated by the var-calling wiring (S1) inside its existing per-locus
//!   genotype pass, not here.
//!
//! Then `F = clip(1 − obs_het / Hexp, 0, 0.99)`. Both functions here are pure;
//! `Hexp` is supplied by the caller.

use crate::sample_summary::HetCounts;

/// The ceiling on `F`: a fully-selfing sample never reaches exactly 1 (spec
/// §3; matches the prototype's `np.clip(..., 0.0, 0.99)`).
pub const MAX_INBREEDING_COEFFICIENT: f64 = 0.99;

/// The sample's observed heterozygosity **rate**:
/// `n_het_sites / callable_positions`. Returns `0.0` when there are no
/// callable positions. A zero-callable sample thus maps to `F = 0.99`
/// (maximally inbred via [`inbreeding_coefficient`]), but such a sample has
/// no usable per-locus coverage/allele data, so its `F` is never consumed by
/// the scorer — the value is inert. (`n_het_sites` and `callable_positions`
/// come from different summary fields, so the rate can exceed `1.0` on
/// malformed input; that is the caller's data-quality problem, not clamped
/// here.)
pub fn obs_het(het_counts: &HetCounts, callable_positions: u64) -> f64 {
    if callable_positions == 0 {
        return 0.0;
    }
    het_counts.n_het_sites as f64 / callable_positions as f64
}

/// The inbreeding coefficient `F = clip(1 − obs_het / Hexp, 0, 0.99)` from the
/// sample's observed het rate and the cohort's expected heterozygosity `hexp`.
///
/// A non-finite `obs_het`, or a non-positive / non-finite `hexp` (no expected
/// variation, or an unformed accumulator) leaves `F` undefined; it returns
/// `0.0` (treat the sample as outbred) rather than propagating a `NaN` into
/// the Wright prior. This diverges from the prototype only at the
/// `obs_het == 0, hexp → 0` corner: the prototype floors the denominator at
/// `1e-9`, giving `F ≈ 0.99` there, whereas a cohort with no expected
/// variation is treated here as `F` undefined → `0`.
pub fn inbreeding_coefficient(obs_het: f64, hexp: f64) -> f64 {
    if !(obs_het.is_finite() && hexp.is_finite() && hexp > 0.0) {
        return 0.0;
    }
    (1.0 - obs_het / hexp).clamp(0.0, MAX_INBREEDING_COEFFICIENT)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn het_counts(n_het: u64) -> HetCounts {
        HetCounts {
            n_het_sites: n_het,
            n_hom_alt_sites: 0,
            n_ambiguous_sites: 0,
            n_variant_sites: n_het,
            min_depth: 4,
            error_rate: 0.02,
            lr_margin: std::f64::consts::LN_10,
        }
    }

    /// `obs_het` is the het count over *callable positions*, not over variant
    /// sites — 359 het / 1e6 callable ≈ 359 het/Mb.
    #[test]
    fn obs_het_is_rate_over_callable_positions() {
        let h = obs_het(&het_counts(359), 1_000_000);
        assert!((h - 359e-6).abs() < 1e-12);
    }

    /// No callable positions → a `0.0` rate rather than a divide-by-zero.
    #[test]
    fn obs_het_zero_callable_is_zero() {
        assert_eq!(obs_het(&het_counts(0), 0), 0.0);
        assert_eq!(obs_het(&het_counts(10), 0), 0.0);
    }

    /// `F = 1 − Hobs/Hexp`: an outbred sample (`Hobs ≈ Hexp`) gets `F ≈ 0`, a
    /// selfer (`Hobs ≪ Hexp`) gets `F` near the ceiling.
    #[test]
    fn inbreeding_coefficient_maps_het_deficit_to_f() {
        // Hobs == Hexp → F = 0 (outbred, e.g. an F2).
        assert!((inbreeding_coefficient(0.01, 0.01) - 0.0).abs() < 1e-12);
        // Hobs = ½ Hexp → F = 0.5.
        assert!((inbreeding_coefficient(0.005, 0.01) - 0.5).abs() < 1e-12);
        // Near-zero het against real expected variation → clipped at 0.99.
        assert_eq!(inbreeding_coefficient(0.0, 0.01), 0.99);
    }

    /// A het *excess* (`Hobs > Hexp`, e.g. a collapsed-paralog pseudo-het
    /// signal) clips `F` up to `0`, never negative.
    #[test]
    fn inbreeding_coefficient_clips_negative_to_zero() {
        assert_eq!(inbreeding_coefficient(0.02, 0.01), 0.0);
    }

    /// A degenerate `Hexp` (zero / negative / non-finite / `+inf`) yields
    /// `F = 0` rather than an inf/NaN.
    #[test]
    fn inbreeding_coefficient_degenerate_hexp_is_zero() {
        assert_eq!(inbreeding_coefficient(0.01, 0.0), 0.0);
        assert_eq!(inbreeding_coefficient(0.01, -1.0), 0.0);
        assert_eq!(inbreeding_coefficient(0.01, f64::NAN), 0.0);
        assert_eq!(inbreeding_coefficient(0.01, f64::INFINITY), 0.0);
    }

    /// A non-finite `obs_het` (e.g. a `NaN` leaking from an upstream fit) is
    /// guarded to `F = 0` rather than riding through `.clamp` as `NaN` and
    /// poisoning the Wright prior.
    #[test]
    fn inbreeding_coefficient_nan_obs_het_is_zero() {
        assert_eq!(inbreeding_coefficient(f64::NAN, 0.01), 0.0);
        assert_eq!(inbreeding_coefficient(f64::INFINITY, 0.01), 0.0);
    }

    /// The `obs_het == 0, hexp == 0` corner is `F = 0` here (the deliberate
    /// divergence from the prototype's `1e-9`-floored `≈ 0.99`).
    #[test]
    fn inbreeding_coefficient_zero_het_zero_hexp_is_zero() {
        assert_eq!(inbreeding_coefficient(0.0, 0.0), 0.0);
    }

    /// `obs_het` can exceed `1.0` on malformed input (het count from a
    /// different summary than the callable total); it is not clamped here.
    #[test]
    fn obs_het_can_exceed_one_on_malformed_counts() {
        assert!(obs_het(&het_counts(10), 5) > 1.0);
    }
}
