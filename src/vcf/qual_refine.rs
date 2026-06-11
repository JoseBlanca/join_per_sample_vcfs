//! Variant-QUAL refinement — deflate site QUAL when the alternate-allele
//! support looks like a systematic artifact rather than a real variant.
//!
//! ## Why
//!
//! Our genotype likelihood scores each alt-supporting read as an
//! independent rare error, so QUAL grows ~linearly with the number of
//! alt reads. At a systematic-artifact site (mismapping, paralog,
//! recurrent context error) the alt reads recur at a steady fraction, so
//! QUAL inflates with depth even though the site is false. A QUAL that
//! rises with coverage at a false positive is mis-calibrated.
//!
//! The fix (validated on the GIAB HG002 depth sweep, 5–301x — see
//! `doc/devel/reports/qual_fp_depth_inflation_2026-06-10.md`) is to judge
//! the **shape** of the alt evidence, which exposes artifacts *more*
//! clearly with depth, instead of just its amount. Two complementary
//! shape tests, each a Phred reduction on the engine QUAL:
//!
//! * **allele balance** — a real call's alt-read fraction matches what
//!   the called genotypes imply (≈0.5 for a single het, the cohort
//!   allele frequency otherwise). We penalise an alt-read *deficit*
//!   relative to that expectation by the binomial improbability of the
//!   observed split. This sharpens with depth (20 % is unmistakably not
//!   50 % once you have many reads) and is ~0 for a well-balanced call.
//! * **strand / read-position bias** — a real variant's reads are a fair
//!   sample of all reads at the site (both strands, varied positions);
//!   an artifact's reads often pile on one strand or position. We
//!   penalise by the improbability of the alt reads' strand / placed-left
//!   fraction against the REF reads'.
//!
//! The penalties are summed (treated as independent evidence the site is
//! an artifact) and subtracted from QUAL. A clean, evenly-sampled call
//! pays ~0 on both; an artifact pays on whichever test exposes it.
//!
//! ## Scope
//!
//! This refines the **QUAL column only** — genotypes, GQ and allele
//! frequencies are untouched. It operates on the primary ALT (the
//! highest-observation non-REF allele) under a biallelic approximation,
//! which is exact for the common biallelic case and a close approximation
//! for multiallelic sites. The expected alt fraction is read from the
//! *called genotypes* (not the EM frequency, which would adapt to the
//! artifact), so it is cohort-correct: a common variant whose genotypes
//! explain its allele fraction is not penalised.

use super::writable::VcfWritable;

const PHRED: f64 = -10.0;

/// Refine the engine QUAL for one record. Returns `baseline_qual` unchanged
/// when there is no ALT support to reason about.
pub(super) fn refine_qual<R: VcfWritable>(
    record: &R,
    table: &[Vec<u8>],
    baseline_qual: f64,
) -> f64 {
    let n_alleles = record.n_alleles();
    if n_alleles < 2 {
        return baseline_qual;
    }
    let ploidy = f64::from(record.ploidy());
    if ploidy <= 0.0 {
        return baseline_qual;
    }
    let n_samples = record.n_samples();

    // Primary ALT = the non-REF allele with the most observations.
    let mut primary = 0usize;
    let mut best_obs = 0u64;
    for a in 1..n_alleles {
        let mut obs = 0u64;
        for s in 0..n_samples {
            obs += u64::from(record.scalars_row(s)[a].num_obs);
        }
        if obs > best_obs {
            best_obs = obs;
            primary = a;
        }
    }
    if primary == 0 || best_obs == 0 {
        return baseline_qual;
    }

    // Pool counts; accumulate the genotype-expected alt-read total.
    let (mut n_ref, mut fwd_ref, mut pl_ref) = (0.0_f64, 0.0, 0.0);
    let (mut n_alt, mut fwd_alt, mut pl_alt) = (0.0_f64, 0.0, 0.0);
    let mut n_total = 0.0_f64;
    let mut expected_alt_reads = 0.0_f64;
    let best_genotype = record.best_genotype();
    for s in 0..n_samples {
        let row = record.scalars_row(s);
        let mut depth = 0.0_f64;
        for stat in row {
            depth += f64::from(stat.num_obs);
        }
        n_total += depth;

        let r = &row[0];
        n_ref += f64::from(r.num_obs);
        fwd_ref += f64::from(r.fwd);
        pl_ref += f64::from(r.placed_left);

        let al = &row[primary];
        n_alt += f64::from(al.num_obs);
        fwd_alt += f64::from(al.fwd);
        pl_alt += f64::from(al.placed_left);

        // Copies of the primary ALT in this sample's called genotype.
        if let Some(gt) = best_genotype.get(s).and_then(|&g| table.get(g)) {
            let copies = gt.iter().filter(|&&x| usize::from(x) == primary).count();
            expected_alt_reads += (copies as f64 / ploidy) * depth;
        }
    }
    if n_alt < 1.0 || n_total < 1.0 {
        return baseline_qual;
    }

    // Allele-balance penalty: expected alt-read fraction from the called
    // genotypes (0.5 for a single het, cohort AF otherwise). Penalise only
    // an alt *deficit* — our artifacts present fewer alt reads than a real
    // call of this frequency would. Skip when the genotypes expect
    // near-all-alt (p_exp ≥ 0.9): a hom-alt's few missing alt reads are
    // sequencing error, not an artifact, and a binomial against p≈1 would
    // mistake that error for a deficit.
    let p_exp = (expected_alt_reads / n_total).clamp(1e-6, 0.999);
    let balance = if p_exp < 0.9 && n_alt < p_exp * n_total {
        tail_phred(n_alt, n_total, p_exp)
    } else {
        0.0
    };

    // Strand / read-position bias penalty: alt reads should match the REF
    // reads' forward-strand and placed-left fractions.
    let ref_fwd = if n_ref > 0.0 {
        (fwd_ref / n_ref).clamp(0.01, 0.99)
    } else {
        0.5
    };
    let ref_pl = if n_ref > 0.0 {
        (pl_ref / n_ref).clamp(0.01, 0.99)
    } else {
        0.5
    };
    let bias = tail_phred(fwd_alt, n_alt, ref_fwd).max(tail_phred(pl_alt, n_alt, ref_pl));

    (baseline_qual - balance - bias).max(0.0)
}

// ---- self-contained math (no external deps) ----------------------------

/// Lanczos ln Γ(x), x > 0. Good to ~1e-13 relative — ample for QUAL.
fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    // Published Lanczos g=7 coefficients, kept verbatim; the last digit
    // of a couple of them is below f64 resolution (clippy flags it) but
    // dropping it would not change the stored value.
    #[allow(clippy::excessive_precision)]
    const C: [f64; 9] = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        std::f64::consts::PI.ln() - (std::f64::consts::PI * x).sin().ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = C[0];
        let t = x + G + 0.5;
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

fn ln_binom_pmf(k: f64, n: f64, p: f64) -> f64 {
    if p <= 0.0 {
        return if k == 0.0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if p >= 1.0 {
        return if k == n { 0.0 } else { f64::NEG_INFINITY };
    }
    let ln_coeff = ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0);
    ln_coeff + k * p.ln() + (n - k) * (1.0 - p).ln()
}

/// Above this `n`, [`binom_two_sided_p`] switches from the exact
/// O(n) discrete sum to the O(1) normal approximation. Chosen so every
/// calibrated and realistic single-sample depth stays on the exact
/// path (the QUAL depth-sweep calibration ran at ≤301x; the cap leaves
/// wide margin) — preserving byte-identical QUAL there — while bounding
/// the per-record cost for large cohorts (where `n` is the cohort-wide
/// depth) and for adversarial / overflow inputs. The normal
/// approximation is most accurate precisely in the large-`n` regime
/// where it is used.
const EXACT_TAIL_MAX_N: u64 = 2_000;

/// Two-sided binomial-tail p-value for `k ~ Binom(n, p)`: total
/// probability of outcomes no more likely than the observed `k`.
///
/// Exact (discrete sum over all `n+1` outcomes) for `n ≤
/// EXACT_TAIL_MAX_N`; a continuity-corrected normal approximation above
/// that, so the cost is bounded for large cohorts and cannot blow up on
/// pathological depths (e.g. a corrupt `.psp` with `num_obs` near
/// `u32::MAX`, which would otherwise loop billions of times).
fn binom_two_sided_p(k: f64, n: f64, p: f64) -> f64 {
    if n < 1.0 {
        return 1.0;
    }
    let n_i = n.round() as u64;
    let k = k.round();
    if n_i > EXACT_TAIL_MAX_N {
        return binom_two_sided_p_normal(k, n, p);
    }
    let obs = ln_binom_pmf(k, n, p);
    let tol = 1e-7;
    let mut acc = 0.0_f64;
    for i in 0..=n_i {
        let lp = ln_binom_pmf(i as f64, n, p);
        if lp <= obs + tol {
            acc += lp.exp();
        }
    }
    acc.clamp(1e-300, 1.0)
}

/// Continuity-corrected normal approximation to [`binom_two_sided_p`],
/// used when `n` exceeds [`EXACT_TAIL_MAX_N`]. By the CLT, `Binom(n, p)
/// ≈ Normal(np, np(1-p))`; "outcomes no more likely than the observed
/// `k`" become the two tails at least as far from the mean as `k`, so
/// the two-sided tail mass is `2·Φ(−z) = erfc(z/√2)` with `z` the
/// continuity-corrected standardized deviation.
fn binom_two_sided_p_normal(k: f64, n: f64, p: f64) -> f64 {
    let mean = n * p;
    let var = n * p * (1.0 - p);
    if var <= 0.0 {
        // p∈{0,1}: the distribution is degenerate at np. (Callers
        // clamp p away from the bounds, so this is only a guard.)
        return if (k - mean).abs() < 0.5 { 1.0 } else { 1e-300 };
    }
    let sigma = var.sqrt();
    // Half-count continuity correction; never let the deviation go
    // negative (k at the mean → z = 0 → p = 1).
    let dev = ((k - mean).abs() - 0.5).max(0.0);
    let z = dev / sigma;
    erfc(z / std::f64::consts::SQRT_2).clamp(1e-300, 1.0)
}

/// Complementary error function, Numerical Recipes `erfcc` rational
/// approximation — fractional error < 1.2e-7 over the whole real line,
/// far tighter than QUAL's ~0.01-phred resolution needs.
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t
        * (-z * z - 1.265_512_23
            + t * (1.000_023_68
                + t * (0.374_091_96
                    + t * (0.096_784_18
                        + t * (-0.186_288_06
                            + t * (0.278_868_07
                                + t * (-1.135_203_98
                                    + t * (1.488_515_87
                                        + t * (-0.822_152_23 + t * 0.170_872_77)))))))))
            .exp();
    if x >= 0.0 { ans } else { 2.0 - ans }
}

/// Phred of a two-sided binomial tail — the balance / bias penalty.
fn tail_phred(k: f64, n: f64, p: f64) -> f64 {
    (PHRED * binom_two_sided_p(k, n, p).log10()).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ln_gamma_matches_factorial() {
        assert!((ln_gamma(6.0) - 120.0_f64.ln()).abs() < 1e-9);
        assert!(ln_gamma(1.0).abs() < 1e-9);
    }

    #[test]
    fn balance_penalty_grows_with_depth_at_fixed_low_vaf() {
        // VAF = 0.2 held constant vs an expected 0.5 — the penalty must
        // increase with depth (the deeper, the clearer it is not a het).
        let lo = tail_phred(2.0, 10.0, 0.5);
        let hi = tail_phred(40.0, 200.0, 0.5);
        assert!(hi > lo, "penalty should grow with depth: {lo} -> {hi}");
    }

    #[test]
    fn balanced_call_pays_no_penalty() {
        // 50/100 alt against an expected 0.5 is the modal outcome → tail
        // ~1 → ~0 penalty.
        assert!(tail_phred(50.0, 100.0, 0.5) < 1.0);
    }

    #[test]
    fn erfc_matches_known_values() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-6);
        // erfc(1) ≈ 0.157299, erfc(2) ≈ 0.004678 (tables).
        assert!((erfc(1.0) - 0.157_299).abs() < 1e-5);
        assert!((erfc(2.0) - 0.004_678).abs() < 1e-5);
        // Symmetry: erfc(-x) = 2 - erfc(x).
        assert!((erfc(-1.0) - (2.0 - erfc(1.0))).abs() < 1e-12);
    }

    /// Exact two-sided tail via the discrete sum — the reference the
    /// approximation is graded against (independent of the `n` cap).
    fn exact_two_sided(k: f64, n: f64, p: f64) -> f64 {
        let obs = ln_binom_pmf(k, n, p);
        let mut acc = 0.0_f64;
        for i in 0..=(n as u64) {
            let lp = ln_binom_pmf(i as f64, n, p);
            if lp <= obs + 1e-7 {
                acc += lp.exp();
            }
        }
        acc.clamp(1e-300, 1.0)
    }

    #[test]
    fn normal_approx_agrees_with_exact_in_meaningful_range() {
        // Where the penalty can actually move a QUAL (moderate
        // deviations, z ≲ 3 → phred ~10–25), the approximation must
        // track the exact discrete sum to within ~2 phred so the
        // switchover at the n cap is seamless. The deep tail (z ≫ 6),
        // where the normal approx diverges from the binomial, is checked
        // separately — there both saturate QUAL to 0, so the gap is moot.
        let n = EXACT_TAIL_MAX_N as f64;
        for &(k, p) in &[(950.0, 0.5), (935.0, 0.5), (560.0, 0.3), (550.0, 0.3)] {
            let exact_phred = -10.0 * exact_two_sided(k, n, p).log10();
            let approx_phred = -10.0 * binom_two_sided_p_normal(k, n, p).log10();
            let d = (exact_phred - approx_phred).abs();
            assert!(
                d < 2.0,
                "phred mismatch at k={k}, p={p}: exact={exact_phred} approx={approx_phred} (Δ={d})"
            );
        }
    }

    #[test]
    fn normal_approx_saturates_in_the_deep_tail() {
        // A gross deficit (k far below np) is "definitely an artifact"
        // under both methods: each yields a penalty larger than any
        // baseline QUAL, so QUAL clamps to 0 regardless of the exact
        // value. We only require both to be very large, not equal.
        let n = EXACT_TAIL_MAX_N as f64;
        let exact_phred = -10.0 * exact_two_sided(300.0, n, 0.3).log10();
        let approx_phred = -10.0 * binom_two_sided_p_normal(300.0, n, 0.3).log10();
        assert!(exact_phred > 100.0 && approx_phred > 100.0);
    }

    #[test]
    fn large_n_is_bounded_and_sane() {
        // The whole point of the fix: a pathological depth returns a
        // finite p-value in O(1) instead of looping billions of times.
        let p = binom_two_sided_p(u32::MAX as f64 / 5.0, u32::MAX as f64, 0.5);
        assert!((1e-300..=1.0).contains(&p), "p out of range: {p}");
        // A 20% VAF against an expected 0.5 at enormous depth is
        // astronomically improbable → p pinned at the clamp floor.
        assert_eq!(p, 1e-300);
    }
}
