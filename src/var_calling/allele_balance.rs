//! Beta-binomial allele-balance score for the per-genotype allele-ratio filter.
//!
//! A genotype implies an expected alternate-allele fraction: a diploid het
//! should show ALT in ≈50% of reads, a hom-alt in ≈100%. A persistent
//! low-fraction signal (~20% VAF) fits *no* real genotype and is the
//! fingerprint of a paralog/systematic artefact the genotype likelihood cannot
//! see (its het emission is `0.5` per read, independent of the observed split).
//!
//! We score each variant-carrying sample by how well its observed (ref, alt)
//! read counts fit the called genotype's expected balance, as a Beta-Binomial
//! log-likelihood ratio against a free-fraction fit:
//!
//! ```text
//!   AB_LR = logBB(alt | n, mean = p0(genotype), s) − logBB(alt | n, mean = alt/n, s)
//! ```
//!
//! `AB_LR ≈ 0`  → the split is consistent with the genotype (keep).
//! `AB_LR ≪ 0`  → the split fits no real genotype balance (artefact).
//!
//! The overdispersion `s` is the depth-honesty knob. It is what makes the score
//! conservative at low coverage (the Beta-Binomial is wide, so even a skewed
//! split is not far from `p0` in likelihood → `AB_LR` stays near 0 → no false
//! drop) and only decisive at high coverage. The defaults below were validated
//! on the GIAB per_sample benchmark across 5×–300× (see
//! `benchmarks/giab/src/allele_balance_dashboard.py`): at `s = 150` a threshold
//! of `AB_LR < -5` removes ~96% of false positives at 300× for ~0.8% true-call
//! loss, and removes essentially nothing below ~30×.

/// Beta-Binomial overdispersion (concentration). Larger → tighter around the
/// genotype's expected balance. MLE on GIAB 300× true hets ≈ 150.
pub const DEFAULT_AB_CONCENTRATION: f64 = 150.0;

/// Default drop threshold on the allele-balance log-LR. Validated knee on the
/// GIAB coverage sweep (see module docs).
pub const DEFAULT_AB_MIN_LOG_LR: f64 = -5.0;

/// Expected VAF for a fully-alt genotype (hom-alt). Left below 1.0 so a couple
/// of reference-base sequencing errors among otherwise-all-alt reads do not
/// drive the Beta-Binomial likelihood to −∞.
pub const HOM_ALT_EXPECTED_VAF: f64 = 0.98;

/// Clamp keeping the free-fraction fit (and any p0) strictly inside (0, 1).
const VAF_CLAMP: f64 = 1e-6;

/// Expected alternate-allele fraction for a genotype carrying `alt_copies` of
/// `ploidy` total. `0` copies → hom-ref (no expectation to test; returns ~0);
/// a full house → [`HOM_ALT_EXPECTED_VAF`]; otherwise `alt_copies / ploidy`.
#[must_use]
pub fn expected_vaf(alt_copies: u8, ploidy: u8) -> f64 {
    debug_assert!(ploidy > 0, "ploidy must be positive");
    if alt_copies == 0 {
        return VAF_CLAMP;
    }
    if alt_copies >= ploidy {
        return HOM_ALT_EXPECTED_VAF;
    }
    f64::from(alt_copies) / f64::from(ploidy)
}

/// Allele-balance log-likelihood ratio for one sample's `(ref_obs, alt_obs)`
/// against the expected fraction `p0` of its called genotype, with
/// concentration `s`. Returns `0.0` when there is no read evidence
/// (`ref_obs + alt_obs == 0`) — nothing to test.
///
/// The result is `≤ 0` by construction (the free fit can only match the data
/// better than a fixed `p0`); the more negative, the worse the observed split
/// fits the genotype.
#[must_use]
pub fn allele_balance_log_lr(ref_obs: u32, alt_obs: u32, p0: f64, s: f64) -> f64 {
    let n = u64::from(ref_obs) + u64::from(alt_obs);
    if n == 0 {
        return 0.0;
    }
    let n_f = n as f64;
    let k = f64::from(alt_obs);
    let p0 = p0.clamp(VAF_CLAMP, 1.0 - VAF_CLAMP);
    let p_free = (k / n_f).clamp(VAF_CLAMP, 1.0 - VAF_CLAMP);
    beta_binom_logpmf(k, n_f, p0, s) - beta_binom_logpmf(k, n_f, p_free, s)
}

/// `log P(k alt | n total, mean, concentration s)` for a Beta-Binomial whose
/// underlying Beta has shape `(mean·s, (1−mean)·s)`. The `C(n, k)` term cancels
/// in the log-LR but is kept so this is a true log-pmf usable on its own.
fn beta_binom_logpmf(k: f64, n: f64, mean: f64, s: f64) -> f64 {
    let a = mean * s;
    let b = (1.0 - mean) * s;
    let log_choose = ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0);
    log_choose + ln_beta(k + a, n - k + b) - ln_beta(a, b)
}

#[inline]
fn ln_beta(a: f64, b: f64) -> f64 {
    ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
}

/// Lanczos `ln Γ(x)` for `x > 0`. Matches the implementation in
/// `posterior_engine` (kept local so this module is self-contained; a shared
/// math util could fold the copies together later).
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

#[cfg(test)]
mod tests {
    use super::*;

    const S: f64 = DEFAULT_AB_CONCENTRATION;

    #[test]
    fn expected_vaf_maps_genotype_to_balance() {
        assert_eq!(expected_vaf(0, 2), VAF_CLAMP); // hom-ref
        assert!((expected_vaf(1, 2) - 0.5).abs() < 1e-12); // het
        assert_eq!(expected_vaf(2, 2), HOM_ALT_EXPECTED_VAF); // hom-alt
        assert!((expected_vaf(1, 4) - 0.25).abs() < 1e-12); // tetraploid simplex
        assert_eq!(expected_vaf(4, 4), HOM_ALT_EXPECTED_VAF);
    }

    #[test]
    fn balanced_het_scores_near_zero() {
        // 130 alt / 130 ref at high depth: a textbook het.
        let lr = allele_balance_log_lr(130, 130, 0.5, S);
        assert!(lr > -1.0, "balanced deep het should score ~0, got {lr}");
    }

    #[test]
    fn skewed_het_at_high_depth_is_strongly_negative() {
        // ~20% VAF het at 300×: the artefact signature.
        let lr = allele_balance_log_lr(200, 50, 0.5, S);
        assert!(lr < -10.0, "deep 0.2-VAF het should be flagged, got {lr}");
    }

    #[test]
    fn same_skew_at_low_depth_is_not_flagged() {
        // ~20% VAF (1 alt / 5) at 5×: too little evidence to reject a het.
        let lr = allele_balance_log_lr(4, 1, 0.5, S);
        assert!(lr > -5.0, "shallow 0.2-VAF should NOT be flagged, got {lr}");
    }

    #[test]
    fn depth_monotonicity_for_fixed_skew() {
        // The same 0.2 fraction is more decisively rejected with more reads.
        let shallow = allele_balance_log_lr(8, 2, 0.5, S);
        let deep = allele_balance_log_lr(200, 50, 0.5, S);
        assert!(
            deep < shallow,
            "deeper same-skew must score lower: {deep} !< {shallow}"
        );
    }

    #[test]
    fn hom_alt_call_at_full_vaf_clears_drop_threshold() {
        // A clean 1/1 (0 ref / 250 alt) is mildly penalised vs the 0.98
        // expectation (the model anticipates ~2% ref errors), but stays well
        // above the drop threshold. In practice the filter only runs on het
        // calls — hom-alt artefacts sit at VAF≈1 and allele balance can't see
        // them — so this is a guard on the general function, not the gate.
        let lr = allele_balance_log_lr(0, 250, HOM_ALT_EXPECTED_VAF, S);
        assert!(
            lr > DEFAULT_AB_MIN_LOG_LR,
            "clean hom-alt must clear the drop threshold, got {lr}"
        );
    }

    #[test]
    fn no_reads_scores_zero() {
        assert_eq!(allele_balance_log_lr(0, 0, 0.5, S), 0.0);
    }

    #[test]
    fn log_lr_is_non_positive() {
        // Free fit can only match the data at least as well as a fixed p0.
        for (r, a) in [(157u32, 53u32), (100, 100), (10, 1), (0, 80)] {
            let lr = allele_balance_log_lr(r, a, 0.5, S);
            assert!(lr <= 1e-9, "log-LR must be ≤ 0, got {lr} for {r}/{a}");
        }
    }

    #[test]
    fn matches_scipy_reference_values() {
        // Parity with the offline analysis (scipy betaln/gammaln), s=150, p0=0.5.
        // Locks the numerics so a refactor of the gamma kernel can't drift.
        for (r, a, want) in [
            (200u32, 50u32, -19.285_174),
            (130, 130, 0.0),
            (4, 1, -0.934_461),
            (8, 2, -1.814_137),
        ] {
            let got = allele_balance_log_lr(r, a, 0.5, S);
            assert!(
                (got - want).abs() < 1e-5,
                "ref {r}/{a}: got {got}, want {want}"
            );
        }
    }
}
