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

/// Above this `n`, [`binom_two_sided_p`] switches from the exact discrete sum
/// to the closed-form incomplete-beta tail ([`binom_two_sided_p_beta`]). The
/// beta method is itself exact at every `n` and `p`, so this is a
/// *byte-identity* boundary, not an accuracy one: the discrete sum is kept
/// below it only to preserve bit-for-bit-identical QUAL at every
/// calibrated/realistic single-sample depth (the QUAL depth sweep ran at ≤301x;
/// the cap leaves wide margin). Above it the per-record cost is bounded to
/// O(log n) — for large cohorts (where `n` is the cohort-wide depth) and for
/// adversarial / overflow inputs alike.
const EXACT_TAIL_MAX_N: u64 = 2_000;

/// Two-sided binomial-tail p-value for `k ~ Binom(n, p)`: total probability of
/// outcomes no more likely than the observed `k` (the Sterne two-sided test).
///
/// Exact discrete sum over all `n + 1` outcomes for `n ≤ EXACT_TAIL_MAX_N`; the
/// exact incomplete-beta tail above that. Both are exact (the beta method
/// reproduces the discrete sum to floating point); the split exists only to
/// keep QUAL byte-identical at validated depths while bounding the cost to
/// O(log n) on large cohorts and pathological inputs (e.g. a corrupt `.psp`
/// with `num_obs` near `u32::MAX`, which the old `0..=n` sum looped billions of
/// times over).
fn binom_two_sided_p(k: f64, n: f64, p: f64) -> f64 {
    if n < 1.0 {
        return 1.0;
    }
    let n_i = n.round() as u64;
    let k = k.round();
    if n_i > EXACT_TAIL_MAX_N {
        return binom_two_sided_p_beta(k, n, p);
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

/// Lentz continued fraction for the incomplete beta (Numerical Recipes
/// `betacf`). Converges to ~1e-15 in a few dozen iterations regardless of the
/// magnitude of `a`/`b`, which is what makes the tail O(1) in read depth.
fn betacf(a: f64, b: f64, x: f64) -> f64 {
    const MAXIT: usize = 300;
    const EPS: f64 = 1e-15;
    const FPMIN: f64 = 1e-300;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < FPMIN {
        d = FPMIN;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=MAXIT {
        let m_f = m as f64;
        let m2 = 2.0 * m_f;
        // Even step.
        let aa = m_f * (b - m_f) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;
        // Odd step.
        let aa = -(a + m_f) * (qab + m_f) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS {
            break;
        }
    }
    h
}

/// Regularized incomplete beta `I_x(a, b)` (Numerical Recipes `betai`),
/// self-contained via [`ln_gamma`] + [`betacf`]. Good to ~1e-12 relative —
/// ample for a QUAL tail. This is the exact binomial CDF in closed form.
fn reg_incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let ln_front = ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln();
    let front = ln_front.exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * betacf(a, b, x) / a
    } else {
        1.0 - front * betacf(b, a, 1.0 - x) / b
    }
}

/// `P(X <= k)` for `X ~ Binom(n, p)`, exact via the incomplete beta:
/// `P(X <= k) = I_{1-p}(n - k, k + 1)`.
fn binom_cdf_le(k: u64, n: u64, p: f64) -> f64 {
    if k >= n {
        return 1.0;
    }
    reg_incomplete_beta((n - k) as f64, (k + 1) as f64, 1.0 - p)
}

/// `P(X >= k)` for `X ~ Binom(n, p)`.
fn binom_sf_ge(k: u64, n: u64, p: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    1.0 - binom_cdf_le(k - 1, n, p)
}

/// Exact two-sided Sterne tail via the incomplete-beta CDF/SF plus a binary
/// search for the opposite-flank cutoff, used when `n` exceeds
/// [`EXACT_TAIL_MAX_N`]. The binomial is unimodal, so the "no more likely than
/// `k`" set is the near tail through `k` plus a far tail beyond the
/// opposite-side outcome of equal pmf; each flank is monotone, so the far-tail
/// cutoff is one binary search over `ln_binom_pmf`. Cost: O(log n) pmf
/// evaluations + two beta evaluations — independent of read depth.
fn binom_two_sided_p_beta(k: f64, n: f64, p: f64) -> f64 {
    if n < 1.0 {
        return 1.0;
    }
    let n_u = n.round() as u64;
    let k_u = k.round().clamp(0.0, n_u as f64) as u64;
    let n_f = n_u as f64;
    let lpk = ln_binom_pmf(k_u as f64, n_f, p);
    let tol = 1e-7;
    // `floor((n+1)p)` is always an argmax, so `pmf(mode)` is the maximum.
    let mode = (((n_f + 1.0) * p).floor()).clamp(0.0, n_f) as u64;

    // Observed at the peak: every outcome is "no more likely", so p = 1.0.
    if k_u == mode {
        return 1.0;
    }

    let two_sided = if k_u < mode {
        let near = binom_cdf_le(k_u, n_u, p); // {i <= k}
        // Far tail {i >= j}: smallest j on the right flank `[mode, n]` whose pmf
        // is no greater than pmf(k). Empty if even pmf(n) exceeds pmf(k).
        let far = if ln_binom_pmf(n_f, n_f, p) > lpk + tol {
            0.0
        } else {
            let (mut lo, mut hi) = (mode, n_u);
            while lo < hi {
                let mid = lo + (hi - lo) / 2;
                if ln_binom_pmf(mid as f64, n_f, p) <= lpk + tol {
                    hi = mid;
                } else {
                    lo = mid + 1;
                }
            }
            binom_sf_ge(lo, n_u, p)
        };
        near + far
    } else {
        let near = binom_sf_ge(k_u, n_u, p); // {i >= k}
        // Far tail {i <= j}: largest j on the left flank `[0, mode]` whose pmf
        // is no greater than pmf(k).
        let far = if ln_binom_pmf(0.0, n_f, p) > lpk + tol {
            0.0
        } else {
            let (mut lo, mut hi) = (0u64, mode);
            while lo < hi {
                let mid = lo + (hi - lo).div_ceil(2);
                if ln_binom_pmf(mid as f64, n_f, p) <= lpk + tol {
                    lo = mid;
                } else {
                    hi = mid - 1;
                }
            }
            binom_cdf_le(lo, n_u, p)
        };
        near + far
    };
    two_sided.clamp(1e-300, 1.0)
}

/// Phred of a two-sided binomial tail — the balance / bias penalty.
fn tail_phred(k: f64, n: f64, p: f64) -> f64 {
    (PHRED * binom_two_sided_p(k, n, p).log10()).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

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

    // ----- verification of the incomplete-beta tail (above the cap) -----
    //
    // The exact discrete sum is the reference: it was validated against
    // scipy.stats.binomtest to 0.0000 phred (tmp/qual_binom_calib.py), so the
    // incomplete-beta method is correct iff it reproduces the sum wherever the
    // result can change QUAL.

    /// Exact two-sided Sterne tail by discrete enumeration — the reference.
    fn enumerated_two_sided_p(k: f64, n: f64, p: f64) -> f64 {
        if n < 1.0 {
            return 1.0;
        }
        let n_i = n.round() as u64;
        let k = k.round();
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

    fn phred_of(p: f64) -> f64 {
        (PHRED * p.log10()).max(0.0)
    }

    /// Core gate: the incomplete-beta tail must reproduce the enumeration across
    /// the whole operating grid — every `p` incl. the clamp extremes, depths up
    /// to where enumeration is affordable, `k` swept across both tails. Within
    /// the decision band (penalty ≤ 100 phred, where QUAL can change) the two
    /// must be numerically equal; deeper in the tail both saturate, so we only
    /// require both to stay past the band.
    #[test]
    fn beta_tail_matches_enumeration_across_grid() {
        let ps = [1e-6, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.99, 0.999];
        let ns = [1u64, 2, 5, 20, 50, 200, 1000, 2001, 4000];
        let mut worst_band = 0.0_f64;
        for &n in &ns {
            for &p in &ps {
                let nf = n as f64;
                let mode = (((nf + 1.0) * p).floor()).clamp(0.0, nf);
                let sd = (nf * p * (1.0 - p)).sqrt().max(1.0);
                let mut ks: Vec<u64> = vec![0, n, mode as u64];
                let mut z = -8.0;
                while z <= 8.0 {
                    let k = (mode + z * sd).round();
                    if (0.0..=nf).contains(&k) {
                        ks.push(k as u64);
                    }
                    z += 0.25;
                }
                ks.sort_unstable();
                ks.dedup();
                for &k in &ks {
                    let kf = k as f64;
                    let got = phred_of(binom_two_sided_p_beta(kf, nf, p));
                    let want = phred_of(enumerated_two_sided_p(kf, nf, p));
                    if want <= 100.0 {
                        let d = (got - want).abs();
                        worst_band = worst_band.max(d);
                        assert!(
                            d < 0.1,
                            "decision-band mismatch at n={n} p={p} k={k}: \
                             beta={got:.4} enum={want:.4} (Δ={d:.4})"
                        );
                    } else {
                        assert!(
                            got > 100.0,
                            "deep-tail disagreement at n={n} p={p} k={k}: \
                             beta={got:.4} enum={want:.4}"
                        );
                    }
                }
            }
        }
        assert!(
            worst_band < 1e-3,
            "band agreement weaker than expected: {worst_band}"
        );
    }

    /// Dispatch + byte-identity: the production entry point uses the exact
    /// discrete sum at/below the cap (bit-for-bit) and the beta method above it.
    #[test]
    fn production_dispatch_matches_path_by_cap() {
        let cases = [(20.0, 0.3), (50.0, 0.5), (180.0, 0.1), (995.0, 0.5)];
        for &n in &[100.0, 1000.0, EXACT_TAIL_MAX_N as f64] {
            for &(k, p) in &cases {
                if k <= n {
                    assert_eq!(
                        binom_two_sided_p(k, n, p),
                        enumerated_two_sided_p(k, n, p),
                        "≤cap path must be the exact sum at n={n} k={k} p={p}"
                    );
                }
            }
        }
        let n = (EXACT_TAIL_MAX_N + 1) as f64;
        for &(k, p) in &[(900.0, 0.5), (300.0, 0.3), (2001.0, 0.9)] {
            assert_eq!(
                binom_two_sided_p(k, n, p),
                binom_two_sided_p_beta(k.round(), n, p),
                ">cap path must be the beta method at n={n} k={k} p={p}"
            );
        }
    }

    #[test]
    fn reg_incomplete_beta_boundaries_and_symmetry() {
        assert_eq!(reg_incomplete_beta(2.0, 3.0, 0.0), 0.0);
        assert_eq!(reg_incomplete_beta(2.0, 3.0, 1.0), 1.0);
        // I_x(a,b) = 1 - I_{1-x}(b,a).
        for (a, b, x) in [(2.0, 3.0, 0.4), (7.5, 2.0, 0.6), (50.0, 20.0, 0.7)] {
            let lhs = reg_incomplete_beta(a, b, x);
            let rhs = 1.0 - reg_incomplete_beta(b, a, 1.0 - x);
            assert!((lhs - rhs).abs() < 1e-12, "symmetry broke at ({a},{b},{x})");
        }
    }

    #[test]
    fn reg_incomplete_beta_matches_scipy_golden() {
        // Pinned from scipy.special.betainc.
        let cases = [
            (2.0, 3.0, 0.4, 0.524_800_000_000_000),
            (5.0, 5.0, 0.5, 0.500_000_000_000_000),
            (10.0, 1.0, 0.9, 0.348_678_440_100_000),
            (0.5, 0.5, 0.3, 0.369_010_119_565_545),
            (50.0, 20.0, 0.7, 0.382_509_248_381_233),
            (3.0, 7.0, 0.2, 0.261_802_496_000_000),
        ];
        for (a, b, x, want) in cases {
            let got = reg_incomplete_beta(a, b, x);
            assert!(
                (got - want).abs() < 1e-10,
                "I_{x}({a},{b})={got}, want {want}"
            );
        }
    }

    #[test]
    fn binom_cdf_matches_direct_summation() {
        for (n, p) in [(10u64, 0.3), (25, 0.5), (40, 0.85)] {
            for k in 0..=n {
                let direct: f64 = (0..=k)
                    .map(|i| ln_binom_pmf(i as f64, n as f64, p).exp())
                    .sum();
                let beta = binom_cdf_le(k, n, p);
                assert!(
                    (direct - beta).abs() < 1e-9,
                    "CDF mismatch n={n} p={p} k={k}: direct={direct} beta={beta}"
                );
            }
        }
    }

    /// At depths the enumeration can't reach (millions), the beta path must
    /// still obey the structural invariants of a two-sided tail.
    #[test]
    fn beta_invariants_at_huge_depth() {
        for &p in &[1e-6, 0.01, 0.5, 0.99, 0.999] {
            let n = 2_000_000.0_f64;
            let mode = ((n + 1.0) * p).floor();
            assert!(
                binom_two_sided_p_beta(mode, n, p) > 0.999_999,
                "p at mode should be ~1 (p={p})"
            );
            let mut prev = binom_two_sided_p_beta(mode, n, p);
            let sd = (n * p * (1.0 - p)).sqrt().max(1.0);
            for step in 1..=8 {
                let k = (mode + step as f64 * sd).min(n);
                let cur = binom_two_sided_p_beta(k, n, p);
                assert!(cur.is_finite() && (0.0..=1.0).contains(&cur));
                assert!(cur <= prev + 1e-12, "not monotone right of mode (p={p})");
                prev = cur;
            }
        }
    }

    #[test]
    fn large_n_is_bounded_and_sane() {
        // The whole point of the fix: a pathological depth returns a
        // finite p-value in O(log n) instead of looping billions of times.
        let p = binom_two_sided_p(u32::MAX as f64 / 5.0, u32::MAX as f64, 0.5);
        assert!((1e-300..=1.0).contains(&p), "p out of range: {p}");
        // A 20% VAF against an expected 0.5 at enormous depth is
        // astronomically improbable → p pinned at the clamp floor.
        assert_eq!(p, 1e-300);
    }

    proptest! {
        /// Fuzz the grid gaps: a random (n, p, k) where enumeration is still
        /// affordable must match the beta method in the decision band.
        #[test]
        fn proptest_beta_matches_enumeration(
            n in 1u64..3000,
            p in 1e-6_f64..0.999_999,
            frac in 0.0_f64..1.0,
        ) {
            let k = (frac * n as f64).round();
            let got = phred_of(binom_two_sided_p_beta(k, n as f64, p));
            let want = phred_of(enumerated_two_sided_p(k, n as f64, p));
            if want <= 100.0 {
                prop_assert!(
                    (got - want).abs() < 0.1,
                    "n={n} p={p} k={k}: beta={got} enum={want}"
                );
            } else {
                prop_assert!(got > 100.0, "n={n} p={p} k={k}: beta={got} enum={want}");
            }
        }
    }
}
