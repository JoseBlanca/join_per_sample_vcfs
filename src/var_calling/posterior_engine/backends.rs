//! Pluggable math backends for [`PosteriorEngine`].
//!
//! The trait surface (`MathBackend::ln` / `::exp`) is the only API the
//! engine uses for natural-log and exponential operations in its hot
//! path. Concrete backends live in this module so consumers can swap
//! one implementation for another via
//! [`PosteriorEngine::with_math_backend`].
//!
//! Three backends ship today:
//!
//! - [`ExactMath`] — calls `f64::ln` / `f64::exp` directly. Bit-identical
//!   to the unoptimised engine. Opt-in via
//!   [`PosteriorEngine::with_math_backend`] when reproducibility against
//!   the unoptimised engine matters.
//! - [`InterpUnivariateMath`] — IEEE-decomposition + 1D linear-interp
//!   lookup tables. Approximate (~`1e-6` relative error per call) and
//!   ~10 % faster than [`ExactMath`] in pure scalar code. Opt-in; useful
//!   when the SIMD backend isn't suitable (e.g. testing the scalar
//!   approximation in isolation).
//! - [`InterpUnivariateSimdMath`] — lane-of-4 `ln` / `exp` via `wide`'s
//!   built-in polynomial approximation. The **default** for
//!   [`PosteriorEngine::new`] / [`PosteriorEngine::with_config`]: same
//!   parity-budget margin as [`InterpUnivariateMath`] and ~25 % faster
//!   than [`ExactMath`] on the contam-on bench. See
//!   `doc/devel/reports/implementations/posterior_engine_simd_analysis_2026-05-18.md`.
//!
//! The trait is `Sync` because backends are constructed once and shared
//! across the EM loop's reads; no mutable state is required for any
//! backend (all three are zero-sized — the interp backends read
//! `'static`-lifetime tables defined in the sibling `interp` module).
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

    /// Lane-of-4 natural log. Default falls back to four scalar `ln`
    /// calls; SIMD backends override with a lane-parallel
    /// implementation. The engine routes its hottest log/exp calls
    /// through `ln_x4` / `exp_x4` only on the
    /// [`Self::HAS_LANE_4`]-true code path.
    ///
    /// Takes and returns [`wide::f64x4`] directly so the SIMD call
    /// sites in `e_step_simd` / `log_sum_exp_*_x4` don't pay an array
    /// round-trip per call.
    #[inline]
    fn ln_x4(&self, x: wide::f64x4) -> wide::f64x4 {
        let a = x.to_array();
        wide::f64x4::from([self.ln(a[0]), self.ln(a[1]), self.ln(a[2]), self.ln(a[3])])
    }

    /// Lane-of-4 exponential. Default falls back to four scalar
    /// calls; SIMD backends override.
    #[inline]
    fn exp_x4(&self, x: wide::f64x4) -> wide::f64x4 {
        let a = x.to_array();
        wide::f64x4::from([
            self.exp(a[0]),
            self.exp(a[1]),
            self.exp(a[2]),
            self.exp(a[3]),
        ])
    }

    /// Whether this backend's lane methods are SIMD-accelerated and
    /// the engine should pick its lane-batched `e_step` body.
    /// Defaults to `false`; only [`InterpUnivariateSimdMath`]
    /// currently overrides.
    ///
    /// Resolved at monomorphisation time, so the engine's
    /// `if M::HAS_LANE_4 { … } else { … }` dispatch compiles to a
    /// dead-code-elimination of the unused branch.
    const HAS_LANE_4: bool = false;
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

/// Interpolated `ln` / `exp` via static lookup tables.
///
/// Per-call relative error on typical EM inputs is ~`1e-6` (set by
/// the 1024-bin sub-table resolution in [`super::interp`]); see the
/// module docs there for the IEEE-decomposition trick that linearises
/// the inputs.
///
/// **API stability.** The internals (table layout, resolution, the
/// precise sub-domain partition) may change without a major version
/// bump as the accuracy harness drives the resolution tuning. For
/// stable output across versions, construct the engine via
/// [`super::PosteriorEngine::new`] / [`super::PosteriorEngine::with_config`]
/// and don't name the backend explicitly.
#[derive(Debug, Default, Clone, Copy)]
pub struct InterpUnivariateMath;

impl MathBackend for InterpUnivariateMath {
    #[inline]
    fn ln(&self, x: f64) -> f64 {
        super::interp::ln_approx(x)
    }

    #[inline]
    fn exp(&self, x: f64) -> f64 {
        super::interp::exp_approx(x)
    }
}

/// SIMD-accelerated interpolated math. Same approximation contract as
/// [`InterpUnivariateMath`] for the scalar `ln` / `exp` methods, plus
/// genuine lane-of-4 implementations for `ln_x4` / `exp_x4` that the
/// engine uses on its hot path when this backend is selected.
///
/// On x86_64 this lowers to AVX2 (via `wide`'s `f64x4`); on aarch64
/// it lowers to two NEON `float64x2_t` pairs. Both paths share the
/// same source.
///
/// **API stability.** Same caveat as [`InterpUnivariateMath`] —
/// internals (lane width, table layout) may change without a major
/// version bump.
#[derive(Debug, Default, Clone, Copy)]
pub struct InterpUnivariateSimdMath;

impl MathBackend for InterpUnivariateSimdMath {
    #[inline]
    fn ln(&self, x: f64) -> f64 {
        super::interp::ln_approx(x)
    }

    #[inline]
    fn exp(&self, x: f64) -> f64 {
        super::interp::exp_approx(x)
    }

    // Lane-of-4 routes through `wide`'s built-in polynomial `f64x4::ln`
    // / `f64x4::exp` — fully SIMD via safe_arch AVX2 intrinsics +
    // bytemuck bit-cast. See
    // doc/devel/reports/implementations/posterior_engine_simd_analysis_2026-05-18.md.
    #[inline]
    fn ln_x4(&self, x: wide::f64x4) -> wide::f64x4 {
        x.ln()
    }

    #[inline]
    fn exp_x4(&self, x: wide::f64x4) -> wide::f64x4 {
        x.exp()
    }

    const HAS_LANE_4: bool = true;
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
