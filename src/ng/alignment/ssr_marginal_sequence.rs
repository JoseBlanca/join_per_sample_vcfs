//! Algorithm 5 — the sequence-versus-sequence marginal.
//!
//! The probability that one sequence produced another, summed over the ways it could have,
//! with gaps allowed only near the ends (spec §5.1). It is a faithful port of production's
//! `align_subst` ([src/ssr/cohort/pair_hmm.rs](../../../ssr/cohort/pair_hmm.rs)) — **this
//! function alone**, not the `HipstrModel` that calls it: the stutter half of production's
//! read-likelihood model is the *genotyping's*, and porting it here would put a genetics
//! model inside a module that knows only about alignment (spec §5.1's boundary note).
//!
//! ## Two cases, and the common one is degenerate
//!
//! - **equal lengths** — there is exactly one way to line the two sequences up, so the sum
//!   has a single term and the result is a straight base-by-base comparison under the error
//!   rate. This is what lands here in B1.
//! - **different lengths** — the length difference is absorbed by gaps confined to a couple
//!   of bases of either end; the sum then runs over the few ways that can happen, a genuine
//!   banded forward pass. That lands in B2, and **interior gaps are forbidden** — an indel
//!   in the middle of the repeat is what the *stutter* model explains, so letting this
//!   algorithm explain it too would make base error and slippage indistinguishable.
//!
//! ## Why this works in linear space, and holds a flat rate rather than the `Emission` component
//!
//! It runs in **ordinary probabilities**, exactly as `align_subst` does, and converts to a
//! logarithm once at the trait boundary (B3) — which spec §7 explicitly permits ("run its
//! matrix in ordinary probabilities and take a single logarithm at the end, safe at repeat
//! sizes"). Two reasons make linear the right space here rather than the module's log-space
//! [`Emission`](super::Emission) component:
//!
//! - The B2 forward **sums** probabilities across paths. In linear space that is a
//!   multiply-and-add, as in the source; in log space each add is a `logsumexp`, dearer and
//!   a different rounding — which would break parity against the very function being ported.
//! - The oracle for this port is production's linear `align_subst` (B3 checks parity against
//!   it). Matching its arithmetic means matching its representation.
//!
//! So the flat per-base error rate ε is held **directly**, not as a
//! [`FlatEmission`](super::FlatEmission) (whose scores are logarithms). The tie to the
//! emission component is the shared *validation* — the checked constructor rejects the same
//! out-of-range rates [`FlatEmission::try_new`](super::FlatEmission::try_new) does, via the
//! same [`DomainError::ErrorRate`] — and the shared flat-rate *concept*, not the
//! representation. (Recorded reconciliation: the plan's preconditions note said algorithm 5
//! "reuses the Emission component"; the design authority it cites — spec §5.1, §7, and the
//! arch's "port this function … returns a linear probability, convert at the boundary" —
//! requires the linear form, which the log-space component cannot provide without `exp()`
//! round-trips that defeat parity.)

use crate::ng::types::DomainError;

/// The sequence-versus-sequence marginal aligner (algorithm 5).
///
/// A stateless value carrying only its flat per-base error rate ε, so it can be held by
/// value and shared across every read at a locus — the same shape production's
/// read-likelihood model has (arch §4). Everything that varies per read is a call argument;
/// there is nothing locus-specific, which is why its [`MarginalAligner::Context`] is `()`
/// (landed in B3).
///
/// [`MarginalAligner::Context`]: super::MarginalAligner::Context
#[derive(Debug, Clone, Copy)]
pub struct SsrSequenceMarginal {
    /// The flat per-base error rate ε, in `[0, 1]` — enforced by [`Self::try_new`].
    eps: f64,
}

impl SsrSequenceMarginal {
    /// Build the aligner from a flat per-base error rate ε.
    ///
    /// # Errors
    ///
    /// [`DomainError::ErrorRate`] if `error_rate` is not a finite probability in `[0, 1]` —
    /// the same check [`FlatEmission::try_new`](super::FlatEmission::try_new) makes, since ε
    /// is the sample group's fixed configuration and an out-of-range value is a setup error,
    /// not something to coerce on the hot path.
    pub fn try_new(error_rate: f64) -> Result<Self, DomainError> {
        if !error_rate.is_finite() || !(0.0..=1.0).contains(&error_rate) {
            return Err(DomainError::ErrorRate(error_rate));
        }
        Ok(Self { eps: error_rate })
    }

    /// The equal-length line-up's probability — **linear**, to be turned into a
    /// [`LogProb`](crate::ng::types::LogProb) only at the trait boundary (B3).
    ///
    /// When the two sequences are the same length there is exactly one line-up, so the sum
    /// has a single term: each read base is scored against the reference base in the same
    /// column, a match worth `1 − ε` and a mismatch `ε / 3` (the error split evenly across
    /// the three other bases). The all-match case is `(1 − ε)^len`, taken as a fast path
    /// exactly as the source does.
    ///
    /// # Panics (debug only)
    ///
    /// Debug-asserts the two slices are the same length — the unequal-length forward is
    /// B2's, and the routing between the two cases is added in B2 as well. A release build
    /// with a length mismatch here would score only the shared prefix, which is why the
    /// assertion names the precondition; before B3 wires a [`MarginalAligner`] impl, nothing
    /// external can reach this.
    ///
    /// **The guarantee is structural, not this assertion.** The `debug_assert_eq!` compiles
    /// out of the release build the project runs, so what actually keeps this arm off an
    /// unequal-length input is B2's length routing — which B2 owes a test for. A
    /// `#[should_panic]` test on the assertion is deliberately *not* added: it would pass in
    /// the test build (where `debug_assertions` is on) while proving nothing about release,
    /// the false-confidence trap this project has recorded more than once.
    ///
    /// [`MarginalAligner`]: super::MarginalAligner
    // `dead_code` in the plain lib build only: this is exercised by the B1 tests and wired
    // to the `MarginalAligner` impl in B3, which removes this allow. B2 extends it with the
    // unequal-length forward in between.
    #[allow(dead_code)]
    fn equal_length_probability(&self, read: &[u8], reference: &[u8]) -> f64 {
        debug_assert_eq!(
            read.len(),
            reference.len(),
            "equal_length_probability is the same-length case; the unequal-length forward is B2"
        );
        if read == reference {
            return (1.0 - self.eps).powi(read.len() as i32);
        }
        substitution_product(read, reference, self.eps)
    }
}

/// The closed-form substitution score for two equal-length sequences: a match scores
/// `1 − ε`, a mismatch `ε / 3`, multiplied across the columns. A free function taking ε,
/// mirroring the source's `substitution_product` so the port reads against it.
// `dead_code` in the plain lib build only until B3 wires the trait — see
// `equal_length_probability`.
#[allow(dead_code)]
fn substitution_product(read: &[u8], reference: &[u8], eps: f64) -> f64 {
    let mismatch = eps / 3.0;
    read.iter()
        .zip(reference)
        .map(|(a, b)| if a == b { 1.0 - eps } else { mismatch })
        .product()
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 0.01;

    fn aligner(eps: f64) -> SsrSequenceMarginal {
        SsrSequenceMarginal::try_new(eps).expect("eps in [0, 1]")
    }

    /// An exact match scores `(1 − ε)^len` — the single-term sum with every column a match,
    /// taken through the fast path.
    #[test]
    fn an_exact_match_scores_one_minus_eps_to_the_length() {
        let seq = b"CACACACA";
        let got = aligner(EPS).equal_length_probability(seq, seq);
        assert!((got - (1.0 - EPS).powi(8)).abs() < 1e-15, "got {got}");
    }

    /// One substitution costs **exactly** a factor `(ε/3)/(1−ε)` against the exact match —
    /// one column flips from a `1 − ε` match to an `ε/3` mismatch, and nothing else moves.
    #[test]
    fn one_substitution_costs_a_factor_of_eps_over_three_over_one_minus_eps() {
        let reference = b"CACACACA";
        let read = b"CACAGACA"; // one substitution (C→G) at an interior column
        let a = aligner(EPS);
        let got = a.equal_length_probability(read, reference);

        // Absolute form: seven matches, one mismatch.
        let expected = (1.0 - EPS).powi(7) * (EPS / 3.0);
        assert!(
            (got - expected).abs() < 1e-15,
            "got {got}, expected {expected}"
        );

        // Relative form: the ratio to the exact match is exactly the per-column factor.
        let exact = a.equal_length_probability(reference, reference);
        assert!((got / exact - (EPS / 3.0) / (1.0 - EPS)).abs() < 1e-12);
    }

    /// The empty sequence lines up with itself the one trivial way, probability 1
    /// (`(1 − ε)^0`) — a degenerate but legal input.
    #[test]
    fn two_empty_sequences_score_one() {
        assert_eq!(aligner(EPS).equal_length_probability(b"", b""), 1.0);
    }

    /// **At the ε endpoints the arithmetic degenerates, and the linear port does not floor**
    /// — unlike the log-space `FlatEmission`, whose match score is floored off zero. Both
    /// endpoints are admitted by `try_new`, so both are pinned:
    ///
    /// - ε = 0: a match is certain (`1 − 0 = 1`), so an exact line-up scores 1 and a single
    ///   mismatch (factor `0/3 = 0`) zeroes the whole product.
    /// - ε = 1: a match is impossible (`1 − 1 = 0`), so an exact line-up scores `0^len = 0`
    ///   through the fast path — and a mismatch (`1/3`) then *outscores* a match, exactly as
    ///   the source's un-floored linear arithmetic does.
    #[test]
    fn the_epsilon_endpoints_degenerate_without_flooring() {
        let zero = aligner(0.0);
        assert_eq!(zero.equal_length_probability(b"ACGT", b"ACGT"), 1.0);
        assert_eq!(zero.equal_length_probability(b"ACGT", b"ACGA"), 0.0);

        let one = aligner(1.0);
        assert_eq!(one.equal_length_probability(b"ACGT", b"ACGT"), 0.0);
        let all_mismatch = one.equal_length_probability(b"AAAA", b"CCCC");
        assert!(
            (all_mismatch - (1.0f64 / 3.0).powi(4)).abs() < 1e-18,
            "got {all_mismatch}"
        );
        // A mismatch outscores a match at ε = 1: (1/3)^4 > 0 = 0^4.
        assert!(all_mismatch > one.equal_length_probability(b"AAAA", b"AAAA"));
    }

    /// The `read == reference` fast path is an *optimization* of the general substitution
    /// product, not a different answer: on an exact match the two must agree (to within the
    /// last bit — `.powi` and repeated multiplication may round differently). Production
    /// never compares them on the same input; this pins that the fast path is exact.
    #[test]
    fn the_fast_path_agrees_with_the_substitution_product_on_an_exact_match() {
        let seq = b"CACAGTCA";
        let fast = aligner(EPS).equal_length_probability(seq, seq); // read == reference
        let product = substitution_product(seq, seq, EPS); // the general product, same input
        assert!(
            (fast - product).abs() < 1e-15,
            "fast {fast}, product {product}"
        );
    }

    /// Every column a mismatch is `(ε/3)^len` — the opposite extreme from the exact match,
    /// pinning that the substitution product runs over all columns, not just the differing
    /// prefix.
    #[test]
    fn an_all_mismatch_line_up_scores_eps_over_three_to_the_length() {
        // AAAA vs CCCC: four mismatches, no shared base.
        let got = aligner(EPS).equal_length_probability(b"AAAA", b"CCCC");
        assert!((got - (EPS / 3.0).powi(4)).abs() < 1e-18, "got {got}");
    }

    /// The checked constructor rejects the same out-of-range rates `FlatEmission::try_new`
    /// does — a non-finite or out-of-`[0, 1]` ε is a setup error, not something to coerce.
    #[test]
    fn try_new_rejects_rates_outside_the_unit_interval() {
        assert!(SsrSequenceMarginal::try_new(0.0).is_ok());
        assert!(SsrSequenceMarginal::try_new(1.0).is_ok());
        assert!(matches!(
            SsrSequenceMarginal::try_new(-0.01),
            Err(DomainError::ErrorRate(_))
        ));
        assert!(matches!(
            SsrSequenceMarginal::try_new(1.01),
            Err(DomainError::ErrorRate(_))
        ));
        assert!(SsrSequenceMarginal::try_new(f64::NAN).is_err());
        assert!(SsrSequenceMarginal::try_new(f64::INFINITY).is_err());
    }
}
