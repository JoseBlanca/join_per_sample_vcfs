//! Univariate interpolated approximations of `f64::ln` and `f64::exp`.
//!
//! Both functions exploit the IEEE 754 mantissa/exponent decomposition
//! so the actual tabulated function only varies over a small bounded
//! sub-domain:
//!
//! - `ln(x) = e * LN_2 + ln(m)` for `m ∈ [1, 2)`. We tabulate `ln(m)`.
//! - `exp(y) = 2^(y · LOG2_E)`. Let `t = y · LOG2_E`, `n = floor(t)`,
//!   `f = t - n`. Then `exp(y) = 2^n · 2^f`. We tabulate `2^f` for
//!   `f ∈ [0, 1)`. The `2^n` factor is constructed by writing the
//!   biased exponent directly into the IEEE 754 bit pattern.
//!
//! Inside each sub-domain the function's second derivative is mild
//! (`|f''(m)| ∈ [0.25, 1]` for `ln(m)`, varies by 2× for `2^f`), so
//! uniform spacing gives near-uniform output error for linear
//! interpolation. The plan §"A / Resolution is a tunable" calls for a
//! starting resolution of 1024 entries per sub-table — the accuracy
//! harness in step 5 drives the final choice.
//!
//! Out-of-domain inputs (`0.0` / negatives / `NaN` / `±∞` /
//! subnormals) fall back to the native `f64::ln` or `f64::exp` so the
//! IEEE semantics of those edge cases are preserved exactly.

use std::sync::LazyLock;

/// Number of *bins* in each sub-table (one less than the entry count
/// because we keep a trailing sentinel so linear interp doesn't need
/// a bounds check on `idx + 1`).
const TABLE_BINS: usize = 1024;

/// Number of mantissa bits we route to the sub-table index. With 1024
/// bins this is 10 bits, leaving 42 mantissa bits below for the
/// linear-interp fraction.
const TABLE_INDEX_BITS: u32 = 10;

const _: () = assert!(1 << TABLE_INDEX_BITS == TABLE_BINS);

/// `f64` significand width: 52 mantissa bits.
const MANTISSA_BITS: u32 = 52;
const MANTISSA_MASK: u64 = (1u64 << MANTISSA_BITS) - 1;
const FRAC_BITS: u32 = MANTISSA_BITS - TABLE_INDEX_BITS;
const FRAC_MASK: u64 = (1u64 << FRAC_BITS) - 1;
const FRAC_SCALE: f64 = 1.0_f64 / ((1u64 << FRAC_BITS) as f64);

/// `ln(m)` sampled at `m_i = 1 + i / TABLE_BINS` for `i ∈ [0, TABLE_BINS]`.
/// The `+1` trailing entry lets the linear-interp lookup read
/// `[idx + 1]` without a bounds check.
static LN_MANTISSA_TABLE: LazyLock<[f64; TABLE_BINS + 1]> = LazyLock::new(|| {
    let mut t = [0.0_f64; TABLE_BINS + 1];
    for (i, slot) in t.iter_mut().enumerate() {
        let m = 1.0 + (i as f64) / (TABLE_BINS as f64);
        *slot = m.ln();
    }
    t
});

/// `2^f` sampled at `f_i = i / TABLE_BINS` for `i ∈ [0, TABLE_BINS]`.
static EXP2_FRAC_TABLE: LazyLock<[f64; TABLE_BINS + 1]> = LazyLock::new(|| {
    let mut t = [0.0_f64; TABLE_BINS + 1];
    for (i, slot) in t.iter_mut().enumerate() {
        let f = (i as f64) / (TABLE_BINS as f64);
        *slot = f.exp2();
    }
    t
});

/// Approximate `f64::ln(x)` via the IEEE-decomposition + table lookup
/// described in the module docs.
#[inline]
pub(super) fn ln_approx(x: f64) -> f64 {
    // Fall back to native for non-positive, non-finite, and subnormal
    // inputs. The hot path's inputs are probabilities in `(0, 1]` and
    // intermediate likelihoods that stay in the normal range; the
    // fallback covers numerical edge cases (e.g. `safe_ln(0.0)`).
    if !(x.is_finite() && x > 0.0) {
        return x.ln();
    }
    let bits = x.to_bits();
    let biased_exp = ((bits >> MANTISSA_BITS) & 0x7FF) as i64;
    if biased_exp == 0 {
        // Subnormal — uncommon in the EM hot path; defer to native.
        return x.ln();
    }
    let e_unbiased = biased_exp - 1023;
    let mantissa_bits = bits & MANTISSA_MASK;
    let idx = (mantissa_bits >> FRAC_BITS) as usize;
    let frac = ((mantissa_bits & FRAC_MASK) as f64) * FRAC_SCALE;

    let table = &*LN_MANTISSA_TABLE;
    // Linear interp between adjacent samples. The trailing sentinel at
    // `table[TABLE_BINS]` makes `idx + 1` always in-bounds for
    // `idx ∈ [0, TABLE_BINS - 1]`.
    let ln_m = table[idx] + frac * (table[idx + 1] - table[idx]);

    (e_unbiased as f64) * std::f64::consts::LN_2 + ln_m
}

/// Approximate `f64::exp(y)` via `2^(y · LOG2_E)` + table lookup.
#[inline]
pub(super) fn exp_approx(y: f64) -> f64 {
    if !y.is_finite() {
        return y.exp();
    }
    // Compose the per-call constant `LOG2_E * TABLE_BINS` so we can do
    // one fmul instead of two, and bake the bin-index extraction into
    // a single floor + bitmask (TABLE_BINS is a power of two, so the
    // `% TABLE_BINS` and `/ TABLE_BINS` decompose to mask + shift —
    // no fdiv).
    const SCALE: f64 = std::f64::consts::LOG2_E * (TABLE_BINS as f64);
    let scaled = y * SCALE;
    // Range check on the un-scaled exponent. `t = scaled / TABLE_BINS`
    // must stay in `[-1022, 1023]` so the exponent injection below
    // doesn't escape the normal range.
    if !(-1022.0 * (TABLE_BINS as f64)..=1023.0 * (TABLE_BINS as f64)).contains(&scaled) {
        return y.exp();
    }
    let total = scaled.floor();
    let frac = scaled - total;
    let total_i = total as i32;
    let n_i = total_i >> TABLE_INDEX_BITS;
    let idx = (total_i as usize) & (TABLE_BINS - 1);

    let table = &*EXP2_FRAC_TABLE;
    let pow2_f = table[idx] + frac * (table[idx + 1] - table[idx]);

    // 2^n via direct IEEE 754 bit injection: a biased exponent of
    // `n + 1023` in bits 52..62 with all other bits zero encodes
    // `2^n` exactly.
    let pow2_n = f64::from_bits(((n_i + 1023) as u64) << MANTISSA_BITS);

    pow2_f * pow2_n
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Relative error tolerance per call. With 1024 bins on a smooth
    /// `ln(m)` over `[1, 2)`, linear interp gives ~`1e-6` worst-case
    /// relative error (the second derivative is bounded by 1 and the
    /// bin width is ~`1e-3`); we test against a slightly looser bound.
    const REL_TOL: f64 = 5e-6;

    fn assert_relative(a: f64, b: f64, tol: f64, label: &str) {
        let scale = b.abs().max(1.0);
        let err = (a - b).abs() / scale;
        assert!(
            err <= tol,
            "{label}: |{a} - {b}| / {scale} = {err} > {tol}"
        );
    }

    // ---- ln_approx ----

    #[test]
    fn ln_approx_matches_native_across_typical_inputs() {
        for &x in &[1e-12_f64, 0.001, 0.1, 0.5, 1.0, 2.0, 5.0, 100.0, 1e6, 1e12] {
            assert_relative(ln_approx(x), x.ln(), REL_TOL, &format!("ln_approx({x})"));
        }
    }

    #[test]
    fn ln_approx_zero_returns_neg_infinity() {
        assert_eq!(ln_approx(0.0), f64::NEG_INFINITY);
        assert_eq!(ln_approx(-0.0), f64::NEG_INFINITY);
    }

    #[test]
    fn ln_approx_negative_input_returns_nan() {
        assert!(ln_approx(-1.0).is_nan());
    }

    #[test]
    fn ln_approx_nan_propagates() {
        assert!(ln_approx(f64::NAN).is_nan());
    }

    #[test]
    fn ln_approx_positive_infinity_returns_positive_infinity() {
        assert_eq!(ln_approx(f64::INFINITY), f64::INFINITY);
    }

    #[test]
    fn ln_approx_one_returns_zero() {
        // `ln(1.0) == 0.0` exactly — checks the `e_unbiased = 0`,
        // `mantissa_bits = 0`, `idx = 0` corner.
        assert_eq!(ln_approx(1.0), 0.0);
    }

    // ---- exp_approx ----

    #[test]
    fn exp_approx_matches_native_across_typical_inputs() {
        for &y in &[-20.0_f64, -1.0, -0.1, 0.0, 0.5, 1.0, 5.0, 20.0] {
            assert_relative(exp_approx(y), y.exp(), REL_TOL, &format!("exp_approx({y})"));
        }
    }

    #[test]
    fn exp_approx_zero_returns_one() {
        assert_eq!(exp_approx(0.0), 1.0);
    }

    #[test]
    fn exp_approx_neg_infinity_returns_zero() {
        assert_eq!(exp_approx(f64::NEG_INFINITY), 0.0);
    }

    #[test]
    fn exp_approx_positive_infinity_returns_positive_infinity() {
        assert_eq!(exp_approx(f64::INFINITY), f64::INFINITY);
    }

    #[test]
    fn exp_approx_nan_propagates() {
        assert!(exp_approx(f64::NAN).is_nan());
    }

    // ---- Round-trip sanity ----

    #[test]
    fn ln_then_exp_is_close_to_identity() {
        for &x in &[0.001_f64, 0.1, 0.5, 1.0, 2.0, 10.0, 1000.0] {
            // Two interp calls compose: total relative error ~2x per-call.
            assert_relative(
                exp_approx(ln_approx(x)),
                x,
                2.0 * REL_TOL,
                &format!("exp(ln({x}))"),
            );
        }
    }
}
