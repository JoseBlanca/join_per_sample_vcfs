//! Pluggable math backends for [`PosteriorEngine`].
//!
//! The trait surface (`MathBackend::ln` / `::exp`) is the only API the
//! engine uses for natural-log and exponential operations in its hot
//! path. Concrete backends live in this module so consumers can swap
//! one implementation for another via
//! [`PosteriorEngine::with_math_backend`].
//!
//! Two backends ship today:
//!
//! - [`ExactMath`] — calls `f64::ln` / `f64::exp` directly. Bit-identical
//!   to the unoptimised engine and the default for
//!   [`PosteriorEngine::new`] / [`PosteriorEngine::with_config`].
//! - Future: `InterpUnivariateMath` — interpolated lookup tables over a
//!   uniform grid in the IEEE 754 mantissa. Approximate but ~4× faster
//!   per call. Lands in a later step of the interpolated-LUTs plan.
//!
//! The trait is `Sync` because backends are constructed once and shared
//! across the EM loop's reads; no mutable state is required for either
//! backend (`ExactMath` is zero-sized, `InterpUnivariateMath` will hold
//! `'static` tables).
//!
//! [`PosteriorEngine`]: super::PosteriorEngine
//! [`PosteriorEngine::new`]: super::PosteriorEngine::new
//! [`PosteriorEngine::with_config`]: super::PosteriorEngine::with_config
//! [`PosteriorEngine::with_math_backend`]: super::PosteriorEngine::with_math_backend

/// Math operations the posterior engine routes through a backend.
///
/// Hot calls (`mix.ln()` in the mixture pre-pass, the per-cell
/// normalisation `exp()` in `e_step`, `log_sum_exp_2` /
/// `log_sum_exp_slice` internals) all dispatch through these
/// methods. End-of-record `log10` calls and the integer
/// `log_factorial` helper stay native — they're outside the inner EM
/// loop.
///
/// Implementations should be cheap to call by reference: the engine
/// inlines `&self` through monomorphisation, so zero-sized
/// [`ExactMath`] has no overhead.
///
/// `Sync` is sufficient for today's single-threaded engine and keeps
/// the door open for rayon-over-records parallelism (each thread
/// reads the same `'static` tables, no `Send`-only state).
pub trait MathBackend: Sync {
    /// Natural logarithm. Inputs `<= 0.0`, `NaN`, and `+∞` follow the
    /// IEEE 754 contract of `f64::ln` (NaN, NaN, +∞ respectively);
    /// callers that need a defined behaviour at `0.0` (the typical
    /// "log of a probability" case) wrap the call in `safe_ln`.
    fn ln(&self, x: f64) -> f64;

    /// Exponential. Inputs propagate IEEE 754 semantics directly.
    fn exp(&self, x: f64) -> f64;
}

/// Bit-identical baseline. Delegates `ln` / `exp` to `f64::ln` /
/// `f64::exp`. Use this when reproducibility against the unoptimised
/// engine matters (e.g. comparing against a prior cohort run).
///
/// Zero-sized — the field is purely a type-level tag at runtime.
#[derive(Debug, Default, Clone, Copy)]
pub struct ExactMath;

impl MathBackend for ExactMath {
    #[inline]
    fn ln(&self, x: f64) -> f64 {
        x.ln()
    }

    #[inline]
    fn exp(&self, x: f64) -> f64 {
        x.exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_math_ln_matches_native() {
        let m = ExactMath;
        for x in [0.5_f64, 1.0, std::f64::consts::E, 100.0] {
            assert_eq!(m.ln(x), x.ln());
        }
    }

    #[test]
    fn exact_math_exp_matches_native() {
        let m = ExactMath;
        for y in [-10.0_f64, 0.0, 1.0, 5.0] {
            assert_eq!(m.exp(y), y.exp());
        }
    }

    #[test]
    fn exact_math_ln_zero_is_neg_infinity() {
        assert_eq!(ExactMath.ln(0.0), f64::NEG_INFINITY);
    }

    #[test]
    fn exact_math_exp_neg_infinity_is_zero() {
        assert_eq!(ExactMath.exp(f64::NEG_INFINITY), 0.0);
    }
}
