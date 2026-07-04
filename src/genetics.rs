//! Shared population-genetics primitives.
//!
//! Small, pure building blocks used by more than one part of the caller: the
//! **Wright inbreeding-adjusted Hardy–Weinberg genotype priors** and the
//! **site-frequency-spectrum frequency grid** they are evaluated on. Both the
//! hidden-paralog filter (`crate::paralog`) and the SFS genotype prior
//! (`crate::var_calling::sfs_prior`) need the same maths — keeping one copy here
//! means the two cannot drift.
//!
//! These functions moved out of `crate::paralog::locus_score` verbatim (same
//! floating-point operations, same order, same constants) so the paralog
//! filter's output is unchanged bit-for-bit.

/// A tiny floor clamping an analytically-positive probability away from
/// `ln(0) = −∞`. The clamp is equivalent to `f64` precision for any realistic
/// probability and never bites on the default grids.
pub const PROBABILITY_FLOOR: f64 = 1e-300;

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

#[cfg(test)]
mod tests {
    use super::*;

    /// A degenerate (`n <= 1`) grid returns the interval midpoint rather than
    /// dividing by `n − 1 = 0`.
    #[test]
    fn linear_grid_point_degenerate_grid_returns_midpoint() {
        assert_eq!(linear_grid_point(0, 1, 0.2, 0.8), 0.5);
        assert_eq!(linear_grid_point(0, 0, 0.2, 0.8), 0.5);
        assert_eq!(sfs_grid_point(0, 1, 0.01), 0.5);
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
