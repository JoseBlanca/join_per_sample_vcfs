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
