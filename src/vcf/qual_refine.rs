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
        for a in 0..n_alleles {
            depth += f64::from(row[a].num_obs);
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
    let ref_fwd = if n_ref > 0.0 { (fwd_ref / n_ref).clamp(0.01, 0.99) } else { 0.5 };
    let ref_pl = if n_ref > 0.0 { (pl_ref / n_ref).clamp(0.01, 0.99) } else { 0.5 };
    let bias = tail_phred(fwd_alt, n_alt, ref_fwd).max(tail_phred(pl_alt, n_alt, ref_pl));

    (baseline_qual - balance - bias).max(0.0)
}

// ---- self-contained math (no external deps) ----------------------------

/// Lanczos ln Γ(x), x > 0. Good to ~1e-13 relative — ample for QUAL.
fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
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
        std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().ln()
            - ln_gamma(1.0 - x)
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

/// Two-sided binomial-tail p-value for `k ~ Binom(n, p)`: total
/// probability of outcomes no more likely than the observed `k`.
fn binom_two_sided_p(k: f64, n: f64, p: f64) -> f64 {
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
}
