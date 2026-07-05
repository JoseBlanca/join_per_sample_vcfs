//! Shared population-genetics primitives.
//!
//! Small, pure building blocks used by more than one part of the caller: the
//! **Wright inbreeding-adjusted Hardy–Weinberg genotype priors**, the
//! **site-frequency-spectrum frequency grid**, and the **Dirichlet-multinomial**
//! genotype prior with its `α`-from-`θ` mapping. Both the hidden-paralog filter
//! (`crate::paralog`) and the SFS genotype prior in the posterior engine
//! (`crate::var_calling::posterior_engine`) build on these — keeping one copy
//! here means they cannot drift.
//!
//! These functions moved out of `crate::paralog::locus_score` verbatim (same
//! floating-point operations, same order, same constants) so the paralog
//! filter's output is unchanged bit-for-bit.

/// A tiny floor clamping an analytically-positive probability away from
/// `ln(0) = −∞`. The clamp is equivalent to `f64` precision for any realistic
/// probability and never bites on the default grids.
pub const PROBABILITY_FLOOR: f64 = 1e-300;

/// Natural log of the gamma function, `ln Γ(x)`, for `x > 0`.
///
/// The Dirichlet-multinomial genotype prior needs `ln Γ` at **non-integer**
/// arguments (the concentration `α` is derived from the continuous diversity
/// `θ`), which the integer-only factorial helpers cannot supply. This is a thin
/// wrapper over `libm::lgamma` (the rust-lang port of musl's libm), kept here so
/// the rest of the caller depends on one project-owned name and the accuracy
/// tests live beside it.
///
/// Only ever called with `x > 0` in this crate (Dirichlet concentrations and
/// `α + k` are strictly positive). `libm::lgamma` is defined for other inputs,
/// but this wrapper documents the contract the callers actually rely on.
///
/// The `x > 0` precondition is checked with a debug assertion (zero-cost in
/// release): `libm::lgamma` does not panic on a non-positive argument, so a
/// wiring regression that let `α` reach `0` or go negative would otherwise
/// silently produce a corrupt log-prior. The bare `libm` call in release keeps
/// the hot path allocation- and branch-free.
pub fn lgamma(x: f64) -> f64 {
    debug_assert!(x > 0.0, "lgamma requires x > 0, got {x}");
    libm::lgamma(x)
}

/// The `i`-th of `n` points on the inclusive linear grid `[lo, hi]`. `n <= 1`
/// yields the midpoint (a degenerate grid still returns a usable point).
pub fn linear_grid_point(i: usize, n: usize, lo: f64, hi: f64) -> f64 {
    if n <= 1 {
        return 0.5 * (lo + hi);
    }
    lo + (hi - lo) * (i as f64) / ((n - 1) as f64)
}

/// The `i`-th of `n` folded-SFS grid points on `[grid_inset, 1−grid_inset]`,
/// uniform in `p`. A single-point grid collapses to the midpoint `0.5`.
/// `grid_inset` is the amount the endpoints are pulled in from `0` and `1`
/// (typically `1/(2n)`) so no grid point sits at a degenerate frequency.
pub fn sfs_grid_point(i: usize, n: usize, grid_inset: f64) -> f64 {
    linear_grid_point(i, n, grid_inset, 1.0 - grid_inset)
}

/// The Wright inbreeding-adjusted Hardy–Weinberg genotype log-priors
/// `(hom-ref, het, hom-alt)` at ALT frequency `p` and inbreeding coefficient
/// `f`: `P(het) = 2pq(1−f)`, homozygotes `q²+f·pq` / `p²+f·pq` (`q = 1−p`).
/// Each probability is floored at [`PROBABILITY_FLOOR`] before the `ln` so a
/// zero-probability genotype yields a finite, very negative log-prior rather
/// than `−∞`.
pub fn wright_genotype_log_priors(p: f64, f: f64) -> (f64, f64, f64) {
    let q = 1.0 - p;
    let het = 2.0 * p * q * (1.0 - f);
    let hom_ref = q * q + f * p * q;
    let hom_alt = p * p + f * p * q;
    (
        hom_ref.max(PROBABILITY_FLOOR).ln(),
        het.max(PROBABILITY_FLOOR).ln(),
        hom_alt.max(PROBABILITY_FLOOR).ln(),
    )
}

/// The random-mating (no-inbreeding) genotype log-priors for an arbitrary
/// `(ploidy, n_alleles)` shape, obtained by marginalising the multinomial
/// genotype prior over a Dirichlet frequency prior with concentration `alpha`
/// — the **Dirichlet-multinomial**. This is the general-shape replacement for
/// the biallelic-only [`wright_genotype_log_priors`] × frequency-grid integral,
/// and the site-frequency spectrum enters as `alpha` (small `α_alt` ⇒ a
/// rare-variant SFS; see the SFS-prior architecture doc §9).
///
/// For a genotype `g` with allele-count vector `k(g)` (how many copies of each
/// allele, summing to the ploidy `m`):
///
/// ```text
/// log P_random(g) = log_multinomial_coeff(g)
///                 + Σ_a [ lgamma(α_a + k_a) − lgamma(α_a) ]
/// ```
///
/// # Inputs (the flat arrays the engine's `GenotypeShape` already precomputes)
///
/// - `genotype_allele_counts`: row-major `n_genotypes × n_alleles` copy counts,
///   `[g * n_alleles + a]` = copies of allele `a` in genotype `g`.
/// - `log_multinomial_coeffs`: `ln C(ploidy; k(g))` per genotype (length
///   `n_genotypes`). Passed in rather than recomputed so this primitive cannot
///   drift from the engine's own coefficients.
/// - `n_alleles`: the number of alleles (row stride).
/// - `alpha`: the Dirichlet concentration, one strictly-positive entry per
///   allele (`α_0` = reference). Its scale carries the diversity `θ`.
///
/// Taking the flat arrays instead of the engine's `GenotypeShape` keeps this
/// shared primitive free of a back-reference into the posterior-engine module.
///
/// # Return value — unnormalised by a genotype-independent constant
///
/// The genotype-independent term `lgamma(Σα + m) − lgamma(Σα)` of the exact
/// Dirichlet-multinomial log-pmf is **omitted**: it is the same for every
/// genotype, so it cancels once the caller mixes in the Wright-`F` IBD term and
/// normalises. The returned values are therefore log-priors up to that additive
/// constant; `softmax` over them recovers the true random-mating genotype
/// distribution. One `f64` per genotype, in the input's genotype order.
///
/// # Preconditions
///
/// The **structural** invariants — `n_alleles > 0`, `alpha.len() == n_alleles`,
/// and `genotype_allele_counts.len() == log_multinomial_coeffs.len() *
/// n_alleles` — are enforced with hard assertions even in release: a shorter
/// `log_multinomial_coeffs` would otherwise let the `chunks_exact`/`zip` pipe
/// **silently truncate** the result to the wrong length, corrupting every
/// downstream genotype index without a panic. The **value** precondition (every
/// `α_a` finite and `> 0`) is a `debug_assert` only, matching [`lgamma`]: a bad
/// `α` degrades a log-prior but cannot mis-shape the output.
pub fn dirichlet_multinomial_log_priors(
    genotype_allele_counts: &[u32],
    log_multinomial_coeffs: &[f64],
    n_alleles: usize,
    alpha: &[f64],
) -> Vec<f64> {
    assert!(n_alleles > 0, "n_alleles must be > 0");
    assert_eq!(
        alpha.len(),
        n_alleles,
        "alpha must have one entry per allele"
    );
    assert_eq!(
        genotype_allele_counts.len(),
        log_multinomial_coeffs.len() * n_alleles,
        "genotype_allele_counts must be n_genotypes × n_alleles"
    );
    debug_assert!(
        alpha.iter().all(|&a| a > 0.0 && a.is_finite()),
        "every Dirichlet concentration must be finite and > 0, got {alpha:?}"
    );

    // lgamma(α_a) once per allele — the baseline each genotype's term subtracts.
    let lgamma_alpha: Vec<f64> = alpha.iter().map(|&a| lgamma(a)).collect();

    genotype_allele_counts
        .chunks_exact(n_alleles)
        .zip(log_multinomial_coeffs)
        .map(|(counts, &log_coeff)| {
            // Σ_a [ lgamma(α_a + k_a) − lgamma(α_a) ]. A zero count contributes
            // exactly zero (lgamma(α_a + 0) − lgamma(α_a) = 0), so skip it.
            counts.iter().zip(alpha).zip(&lgamma_alpha).fold(
                log_coeff,
                |acc, ((&k, &a), &lgamma_a)| {
                    if k == 0 {
                        acc
                    } else {
                        acc + lgamma(a + f64::from(k)) - lgamma_a
                    }
                },
            )
        })
        .collect()
}

/// The reference-allele Dirichlet concentration `α_ref`. Fixed at `1`, the value
/// that makes the biallelic-diploid het:hom-alt ratio `2·α_ref/(α_alt+1)`
/// approach the defensible **2:1** as `α_alt → 0`. It doubles as the
/// monomorphic-site weight (arch §9.2): with the small `α_alt = θ̂` from
/// [`alpha_from_diversity`] the Dirichlet-multinomial's hom-ref probability comes
/// out at the genetically-correct `1 − 3θ/2`, so no separate invariant mass is
/// needed at the default.
pub const ALPHA_REF: f64 = 1.0;

/// A tiny positive floor for each ALT concentration, so the
/// Dirichlet-multinomial's `lgamma(α_alt)` stays finite when the estimated
/// diversity is exactly zero (a fully invariant cohort, or `--diversity 0`). It
/// sits far below any real diversity (human `θ ≈ 1e-3`), so it never perturbs a
/// genuine estimate; at `θ = 0` it yields an effectively-certain hom-ref prior,
/// matching the biallelic grid path's `θ = 0` behaviour.
pub const MIN_ALT_CONCENTRATION: f64 = 1e-12;

/// The Dirichlet concentration `α = (α_ref, α_alt(1), …, α_alt(k−1))` for the SFS
/// genotype prior at estimated diversity `theta` (`θ̂`).
///
/// - `α_ref = ALPHA_REF = 1`.
/// - The total ALT concentration is `θ̂`, split evenly across the `n_alleles − 1`
///   ALT alleles: `α_alt(a) = θ̂ / (n_alleles − 1)`, floored at
///   [`MIN_ALT_CONCENTRATION`] so it stays strictly positive. Splitting keeps a
///   site's total polymorphism `θ̂` independent of how many ALT alleles it
///   carries.
///
/// Fed to [`dirichlet_multinomial_log_priors`], this yields the clean
/// population-genetics marginals for a biallelic-diploid site (`F = 0`): `P(het)
/// ≈ θ`, `P(hom-alt) ≈ θ/2`, monomorphic weight `≈ 1 − 3θ/2`, and a het:hom-alt
/// ratio that stays `≈ 2:1` at every realistic diversity (because `θ̂` — hence
/// `α_alt` — is always small). Per-sample inbreeding `F` is applied on top by the
/// engine's Wright mixture, not here. See the SFS-prior architecture doc §9.2 for
/// why this is the settled mapping (no calibration constant).
///
/// # Preconditions
///
/// `n_alleles >= 1` is a hard assertion (a zero-allele shape is impossible — every
/// site has a reference allele — and would flow a wrong-length `α` into the
/// Dirichlet-multinomial). A monomorphic shape (`n_alleles == 1`) has no ALT to
/// carry diversity and returns `[ALPHA_REF]`. `theta` finite and `>= 0` is a
/// `debug_assert` (a bad θ degrades the prior but cannot mis-shape `α`).
pub fn alpha_from_diversity(n_alleles: usize, theta: f64) -> Vec<f64> {
    assert!(n_alleles >= 1, "n_alleles must be >= 1");
    debug_assert!(
        theta.is_finite() && theta >= 0.0,
        "theta must be finite and non-negative, got {theta}"
    );

    let n_alt = n_alleles - 1;
    if n_alt == 0 {
        return vec![ALPHA_REF];
    }
    let per_alt = (theta / n_alt as f64).max(MIN_ALT_CONCENTRATION);
    let mut alpha = Vec::with_capacity(n_alleles);
    alpha.push(ALPHA_REF);
    alpha.resize(n_alleles, per_alt);
    alpha
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `ln` of the rising factorial `Π_{j=0}^{k−1} (a + j)` — the closed form of
    /// `lgamma(a + k) − lgamma(a)` for an integer `k`, computed independently
    /// (no `lgamma`) so it is a genuine cross-check of
    /// [`dirichlet_multinomial_log_priors`].
    fn pochhammer_ln(a: f64, k: u32) -> f64 {
        (0..k).map(|j| (a + f64::from(j)).ln()).sum()
    }

    /// The independent oracle for one genotype's Dirichlet-multinomial log-prior:
    /// `log_coeff + Σ_a pochhammer_ln(α_a, k_a)`.
    fn dm_log_prior_oracle(counts: &[u32], log_coeff: f64, alpha: &[f64]) -> f64 {
        log_coeff
            + counts
                .iter()
                .zip(alpha)
                .map(|(&k, &a)| pochhammer_ln(a, k))
                .sum::<f64>()
    }

    fn softmax(logs: &[f64]) -> Vec<f64> {
        let max = logs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let exps: Vec<f64> = logs.iter().map(|&l| (l - max).exp()).collect();
        let total: f64 = exps.iter().sum();
        exps.iter().map(|&e| e / total).collect()
    }

    // ---- Dirichlet-multinomial genotype log-priors (G2) ----

    /// Biallelic-diploid: each genotype's log-prior matches the independent
    /// Pochhammer oracle, and the het:hom-alt ratio is exactly `2α_ref/(α_alt+1)`
    /// (the closed form the biallelic grid prior must reproduce).
    #[test]
    fn dirichlet_multinomial_biallelic_diploid_matches_oracle_and_ratio() {
        // AA, AB, BB with ln C(2;k) = 0, ln2, 0.
        let counts = [2u32, 0, 1, 1, 0, 2];
        let coeffs = [0.0, 2.0_f64.ln(), 0.0];
        let alpha = [1.0_f64, 0.02];
        let got = dirichlet_multinomial_log_priors(&counts, &coeffs, 2, &alpha);

        for (g, chunk) in counts.chunks_exact(2).enumerate() {
            let want = dm_log_prior_oracle(chunk, coeffs[g], &alpha);
            assert!(
                (got[g] - want).abs() < 1e-12,
                "genotype {g}: {} vs {want}",
                got[g]
            );
        }
        // het:hom-alt = 2·α_ref/(α_alt+1).
        let ratio = (got[1] - got[2]).exp();
        let want = 2.0 * alpha[0] / (alpha[1] + 1.0);
        assert!((ratio - want).abs() < 1e-12, "ratio {ratio} vs {want}");
    }

    /// At `α_ref = 1, α_alt → 0` the het:hom-alt ratio is the defensible 2:1 —
    /// the SFS-marginal answer the whole feature rests on.
    #[test]
    fn dirichlet_multinomial_biallelic_ratio_is_two_to_one_at_rare_alt() {
        let counts = [2u32, 0, 1, 1, 0, 2];
        let coeffs = [0.0, 2.0_f64.ln(), 0.0];
        let alpha = [1.0_f64, 1e-9];
        let got = dirichlet_multinomial_log_priors(&counts, &coeffs, 2, &alpha);
        let ratio = (got[1] - got[2]).exp();
        assert!((ratio - 2.0).abs() < 1e-6, "ratio {ratio} not ~2:1");
    }

    /// Biallelic-diploid: the softmax of the returned (constant-omitted)
    /// log-priors equals the true Dirichlet-multinomial pmf, computed
    /// independently — confirms the omitted `lgamma(Σα+m)−lgamma(Σα)` term is
    /// genuinely the normaliser.
    #[test]
    fn dirichlet_multinomial_softmax_is_the_true_pmf() {
        let counts = [2u32, 0, 1, 1, 0, 2];
        let coeffs = [0.0, 2.0_f64.ln(), 0.0];
        let alpha = [1.3_f64, 0.7];
        let got = dirichlet_multinomial_log_priors(&counts, &coeffs, 2, &alpha);
        let probs = softmax(&got);

        // True DM pmf: C(m;k)·Π_a Poch(α_a,k_a) / Poch(Σα, m).
        let sum_alpha: f64 = alpha.iter().sum();
        let denom = pochhammer_ln(sum_alpha, 2);
        for (g, chunk) in counts.chunks_exact(2).enumerate() {
            let ln_pmf = dm_log_prior_oracle(chunk, coeffs[g], &alpha) - denom;
            assert!((probs[g] - ln_pmf.exp()).abs() < 1e-12, "pmf {g}");
        }
        assert!((probs.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    /// Triallelic-diploid (a shape the biallelic grid could not cover): every
    /// genotype matches the oracle and the softmax is a proper distribution.
    #[test]
    fn dirichlet_multinomial_triallelic_diploid_matches_oracle() {
        // Six genotypes: AA, AB, BB, AC, BC, CC (homozygotes 0, hets ln2).
        let counts = [
            2u32, 0, 0, // AA
            1, 1, 0, // AB
            0, 2, 0, // BB
            1, 0, 1, // AC
            0, 1, 1, // BC
            0, 0, 2, // CC
        ];
        let ln2 = 2.0_f64.ln();
        let coeffs = [0.0, ln2, 0.0, ln2, ln2, 0.0];
        let alpha = [1.0_f64, 0.02, 0.02];
        let got = dirichlet_multinomial_log_priors(&counts, &coeffs, 3, &alpha);
        for (g, chunk) in counts.chunks_exact(3).enumerate() {
            let want = dm_log_prior_oracle(chunk, coeffs[g], &alpha);
            assert!((got[g] - want).abs() < 1e-12, "genotype {g}");
        }
        assert!((softmax(&got).iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    /// Tetraploid biallelic (polyploid — the general Wright formula): five
    /// genotypes, each matching the oracle, softmax normalised.
    #[test]
    fn dirichlet_multinomial_tetraploid_biallelic_matches_oracle() {
        // k = 4,3,2,1,0 alt copies; ln C(4;k) = 0, ln4, ln6, ln4, 0.
        let counts = [4u32, 0, 3, 1, 2, 2, 1, 3, 0, 4];
        let coeffs = [0.0, 4.0_f64.ln(), 6.0_f64.ln(), 4.0_f64.ln(), 0.0];
        let alpha = [1.0_f64, 0.02];
        let got = dirichlet_multinomial_log_priors(&counts, &coeffs, 2, &alpha);
        for (g, chunk) in counts.chunks_exact(2).enumerate() {
            let want = dm_log_prior_oracle(chunk, coeffs[g], &alpha);
            assert!((got[g] - want).abs() < 1e-12, "genotype {g}");
        }
        assert!((softmax(&got).iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    // ---- α from diversity θ (G3) ----

    /// The biallelic-diploid genotype probabilities (F = 0) that
    /// [`alpha_from_diversity`] + [`dirichlet_multinomial_log_priors`] produce at
    /// diversity `theta`: `[P(hom-ref), P(het), P(hom-alt)]`.
    fn biallelic_probs_from_theta(theta: f64) -> [f64; 3] {
        let alpha = alpha_from_diversity(2, theta);
        // AA, AB, BB with ln C(2;k) = 0, ln2, 0.
        let counts = [2u32, 0, 1, 1, 0, 2];
        let coeffs = [0.0, 2.0_f64.ln(), 0.0];
        let logs = dirichlet_multinomial_log_priors(&counts, &coeffs, 2, &alpha);
        let p = softmax(&logs);
        [p[0], p[1], p[2]]
    }

    /// `alpha_from_diversity` sets `α_ref = 1` and splits `θ` across the ALTs:
    /// biallelic → `[1, θ]`, triallelic → `[1, θ/2, θ/2]`.
    #[test]
    fn alpha_from_diversity_sets_ref_one_and_splits_theta() {
        assert_eq!(alpha_from_diversity(2, 1e-3), vec![1.0, 1e-3]);
        let tri = alpha_from_diversity(3, 2e-3);
        assert_eq!(tri[0], 1.0);
        assert!((tri[1] - 1e-3).abs() < 1e-18 && (tri[2] - 1e-3).abs() < 1e-18);
        // Total ALT concentration is θ regardless of allele count.
        assert!((tri[1] + tri[2] - 2e-3).abs() < 1e-15);
    }

    /// At the human default θ = 1e-3 the biallelic prior gives the clean SFS
    /// marginals: het ≈ θ, hom-alt ≈ θ/2, monomorphic ≈ 1 − 3θ/2. This is the
    /// settled "choice-2" behaviour — NOT the old grid's inflated hom-ref 0.878.
    #[test]
    fn alpha_from_diversity_reproduces_clean_sfs_marginals() {
        let theta = 1e-3;
        let [homref, het, homalt] = biallelic_probs_from_theta(theta);
        assert!((het - theta).abs() < 5e-6, "het {het} vs θ {theta}");
        assert!(
            (homalt - theta / 2.0).abs() < 5e-6,
            "homalt {homalt} vs θ/2"
        );
        assert!(
            (homref - (1.0 - 1.5 * theta)).abs() < 5e-6,
            "homref {homref} vs 1−3θ/2"
        );
    }

    /// The het:hom-alt ratio stays ≈ 2:1 across the whole realistic diversity
    /// range — the θ-independence spec §4a requires — because `α_alt = θ` is
    /// always small (even a very diverse organism at θ = 0.02 stays near 2:1).
    #[test]
    fn alpha_from_diversity_ratio_is_two_to_one_across_theta() {
        for &theta in &[1e-4, 1e-3, 2e-3, 1e-2, 2e-2] {
            let [_, het, homalt] = biallelic_probs_from_theta(theta);
            let ratio = het / homalt;
            assert!(
                (ratio - 2.0).abs() < 0.05,
                "θ={theta}: ratio {ratio} not ≈ 2:1"
            );
        }
    }

    /// θ scales the variant mass monotonically: a more diverse cohort puts more
    /// prior weight on carrying a variant (less on hom-ref).
    #[test]
    fn alpha_from_diversity_variant_mass_grows_with_theta() {
        let variant = |theta: f64| {
            let [_, het, homalt] = biallelic_probs_from_theta(theta);
            het + homalt
        };
        assert!(variant(1e-4) < variant(1e-3));
        assert!(variant(1e-3) < variant(1e-2));
        // And tracks 3θ/2 at the human default.
        assert!((variant(1e-3) - 1.5e-3).abs() < 1e-5);
    }

    /// Zero diversity floors `α_alt` at [`MIN_ALT_CONCENTRATION`] rather than
    /// producing `α_alt = 0` (which would break the Dirichlet-multinomial's
    /// `lgamma`), yielding an effectively-certain hom-ref prior.
    #[test]
    fn alpha_from_diversity_zero_theta_is_floored_and_homref() {
        let alpha = alpha_from_diversity(2, 0.0);
        assert_eq!(alpha[0], 1.0);
        assert_eq!(alpha[1], MIN_ALT_CONCENTRATION);
        let [homref, _, _] = biallelic_probs_from_theta(0.0);
        assert!(homref > 0.999_999_999, "homref {homref} not ≈ 1");
    }

    /// A monomorphic shape (`n_alleles = 1`) has no ALT to carry diversity and
    /// returns just `[α_ref]`.
    #[test]
    fn alpha_from_diversity_monomorphic_shape_is_ref_only() {
        assert_eq!(alpha_from_diversity(1, 1e-3), vec![1.0]);
    }

    /// A legitimately tiny — but non-zero — diversity passes through unfloored:
    /// `MIN_ALT_CONCENTRATION = 1e-12` sits far below any real cohort's θ, so a
    /// low-diversity cohort at θ = 1e-8 keeps `α_alt = 1e-8`. Guards the doc
    /// claim that the floor "never perturbs a genuine estimate".
    #[test]
    fn alpha_from_diversity_tiny_real_theta_is_not_floored() {
        let alpha = alpha_from_diversity(2, 1e-8);
        assert_eq!(alpha[1], 1e-8, "a real θ=1e-8 was clamped by the floor");
        assert!(alpha[1] > MIN_ALT_CONCENTRATION);
    }

    /// A nonsensical θ > 1 (reachable via `--diversity`) is passed through as
    /// `[1, θ]` rather than rejected — documents the map has no upper clamp (the
    /// value precondition is only θ ≥ 0).
    #[test]
    fn alpha_from_diversity_theta_above_one_passes_through() {
        assert_eq!(alpha_from_diversity(2, 3.0), vec![1.0, 3.0]);
    }

    /// A zero-allele shape is impossible and trips the hard assertion rather than
    /// silently returning a wrong-length `α` (which would then mis-shape the
    /// Dirichlet-multinomial).
    #[test]
    #[should_panic(expected = "n_alleles must be >= 1")]
    fn alpha_from_diversity_panics_on_zero_alleles() {
        alpha_from_diversity(0, 1e-3);
    }

    /// Haploid (`ploidy = 1`) is structurally distinct — each genotype carries a
    /// single allele, `k_a ∈ {0, 1}`, `log_coeff = 0` — so its per-genotype
    /// log-prior is just `ln α_a`, matching the oracle.
    #[test]
    fn dirichlet_multinomial_haploid_matches_oracle() {
        // Two haploid "genotypes": A, B.
        let counts = [1u32, 0, 0, 1];
        let coeffs = [0.0, 0.0];
        let alpha = [1.0_f64, 0.3];
        let got = dirichlet_multinomial_log_priors(&counts, &coeffs, 2, &alpha);
        assert!((got[0] - alpha[0].ln()).abs() < 1e-12);
        assert!((got[1] - alpha[1].ln()).abs() < 1e-12);
    }

    /// An empty genotype set yields an empty result (no genotypes to score),
    /// not a panic.
    #[test]
    fn dirichlet_multinomial_empty_input_is_empty() {
        let got = dirichlet_multinomial_log_priors(&[], &[], 2, &[1.0, 0.02]);
        assert!(got.is_empty());
    }

    /// A wrong-length `log_multinomial_coeffs` (fewer coefficients than
    /// genotypes) trips the hard structural assertion rather than silently
    /// truncating the returned vector — the release-mode failure mode the
    /// `chunks_exact`/`zip` pipe would otherwise hide.
    #[test]
    #[should_panic(expected = "n_genotypes × n_alleles")]
    fn dirichlet_multinomial_panics_on_mismatched_coeff_length() {
        // Three genotypes' worth of counts, only two coefficients.
        let counts = [2u32, 0, 1, 1, 0, 2];
        dirichlet_multinomial_log_priors(&counts, &[0.0, 0.0], 2, &[1.0, 0.02]);
    }

    /// A degenerate (`n <= 1`) grid returns the interval midpoint rather than
    /// dividing by `n − 1 = 0`.
    #[test]
    fn linear_grid_point_degenerate_grid_returns_midpoint() {
        assert_eq!(linear_grid_point(0, 1, 0.2, 0.8), 0.5);
        assert_eq!(linear_grid_point(0, 0, 0.2, 0.8), 0.5);
        assert_eq!(sfs_grid_point(0, 1, 0.01), 0.5);
    }

    /// `lgamma` reproduces the log-factorial identity `ln Γ(n+1) = ln n!` at
    /// integer arguments and the classic `ln Γ(1/2) = ln √π` at a non-integer
    /// one — the case the integer factorial helpers cannot cover.
    #[test]
    fn lgamma_matches_known_values() {
        // ln Γ(1) = ln 0! = 0, ln Γ(2) = ln 1! = 0.
        assert!(lgamma(1.0).abs() < 1e-12);
        assert!(lgamma(2.0).abs() < 1e-12);
        // ln Γ(n+1) = ln n! for a few n.
        for (n, fact) in [(3u32, 6.0), (5, 120.0), (6, 720.0)] {
            let got = lgamma(n as f64 + 1.0);
            assert!(
                (got - (fact as f64).ln()).abs() < 1e-10,
                "lgamma({}) = {got}, want ln {fact}",
                n + 1
            );
        }
        // Half-integer absolute anchors — closed forms with √π, so a shared
        // systematic error (which the relative recurrence test cannot see)
        // would show up here.
        let half_ln_pi = std::f64::consts::PI.ln() / 2.0;
        // ln Γ(1/2) = ln √π.
        assert!((lgamma(0.5) - half_ln_pi).abs() < 1e-12, "Γ(1/2)");
        // ln Γ(3/2) = ½ln π − ln 2.
        assert!(
            (lgamma(1.5) - (half_ln_pi - 2.0_f64.ln())).abs() < 1e-12,
            "Γ(3/2) = {}",
            lgamma(1.5)
        );
        // ln Γ(5/2) = ln(3/4) + ½ln π.
        assert!(
            (lgamma(2.5) - ((3.0_f64 / 4.0).ln() + half_ln_pi)).abs() < 1e-12,
            "Γ(5/2) = {}",
            lgamma(2.5)
        );
    }

    /// The `x > 0` contract is enforced in debug builds: a non-positive argument
    /// trips the debug assertion rather than silently returning a corrupt value.
    #[test]
    #[should_panic(expected = "lgamma requires x > 0")]
    fn lgamma_panics_on_non_positive_in_debug() {
        lgamma(0.0);
    }

    /// `lgamma` satisfies the recurrence `ln Γ(x+1) = ln x + ln Γ(x)` at the
    /// small non-integer arguments the prior actually evaluates (`α ≈ θ`, small).
    #[test]
    fn lgamma_satisfies_recurrence_at_small_args() {
        for &x in &[1e-3, 0.01, 0.3, 1.7] {
            let lhs = lgamma(x + 1.0);
            let rhs = x.ln() + lgamma(x);
            assert!((lhs - rhs).abs() < 1e-10, "x={x}: {lhs} vs {rhs}");
        }
    }

    /// The Wright genotype priors are a proper distribution: `hom-ref + het +
    /// hom-alt = 1` for every `(p, F)` — a coefficient typo would break this.
    #[test]
    fn wright_genotype_priors_sum_to_one() {
        for &p in &[0.01, 0.2, 0.5, 0.9] {
            for &f in &[0.0, 0.3, 0.99] {
                let (a, b, c) = wright_genotype_log_priors(p, f);
                let sum = a.exp() + b.exp() + c.exp();
                assert!((sum - 1.0).abs() < 1e-12, "p={p} F={f} sum={sum}");
            }
        }
    }
}
