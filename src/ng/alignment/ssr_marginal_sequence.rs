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

use super::MarginalAligner;
use crate::ng::types::{DomainError, LogProb};

/// Bases of slop tolerated at each end of a sequence: a gap is admitted only within this
/// many positions of either end, the interior being diagonal-only. **Ported verbatim from
/// production's `FLANK_SLOP`** ([src/ssr/cohort/pair_hmm.rs](../../../ssr/cohort/pair_hmm.rs)),
/// and the restriction it enforces is **load-bearing, not an optimisation**: an indel in the
/// middle of the repeat is what the *stutter* model explains, so letting this algorithm
/// explain it too would make base error and slippage indistinguishable (spec §5.1). Confining
/// gaps to the ends is what keeps the base-error rate and the stutter rate separately
/// estimable.
const FLANK_SLOP: usize = 2;

/// The reused forward-matrix buffer for the unequal-length case (B2), so scoring a locus's
/// reads allocates once rather than per read — the caller owns it and hands it back each
/// call, exactly as production's `HmmScratch` does. Empty until first used and grown on
/// demand; it decides nothing about a result, so its `Default` is safe as the trait's
/// `Scratch` bound requires.
///
/// The equal-length case needs no matrix and ignores this.
#[derive(Debug, Default)]
pub struct SequenceMarginalScratch {
    /// The flattened `(m+1) × (n+1)` forward matrix, `m`/`n` the two sequence lengths.
    /// Compute is band-limited but the buffer is the full matrix — memory is not banded,
    /// which is fine at repeat sizes (production makes the same tradeoff).
    cells: Vec<f64>,
}

impl SequenceMarginalScratch {
    /// A fresh, empty buffer.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

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

    /// The full ported function's probability — **linear**, to become a
    /// [`LogProb`](crate::ng::types::LogProb) only at the trait boundary (B3). This is the
    /// dispatch `align_subst` itself performs
    /// ([src/ssr/cohort/pair_hmm.rs](../../../ssr/cohort/pair_hmm.rs)):
    ///
    /// - **equal lengths** → [`Self::equal_length_probability`], the single-term sum (which
    ///   itself splits the exact-match fast path from the substitution product);
    /// - **different lengths** → [`banded_forward`], summing the few ways the length
    ///   difference is absorbed by **end** gaps.
    ///
    /// **This routing is the real guarantee** that `equal_length_probability` only ever sees
    /// equal-length input — the thing its `debug_assert!` merely documents, since that
    /// assertion compiles out of release. Which is why the length test lives here, not on the
    /// assertion.
    fn linear_probability(
        &self,
        read: &[u8],
        reference: &[u8],
        scratch: &mut SequenceMarginalScratch,
    ) -> f64 {
        if read.len() == reference.len() {
            self.equal_length_probability(read, reference)
        } else {
            banded_forward(read, reference, self.eps, scratch)
        }
    }
}

/// The closed-form substitution score for two equal-length sequences: a match scores
/// `1 − ε`, a mismatch `ε / 3`, multiplied across the columns. A free function taking ε,
/// mirroring the source's `substitution_product` so the port reads against it.
fn substitution_product(read: &[u8], reference: &[u8], eps: f64) -> f64 {
    let mismatch = eps / 3.0;
    read.iter()
        .zip(reference)
        .map(|(a, b)| if a == b { 1.0 - eps } else { mismatch })
        .product()
}

/// Whether position `pos` (0-based) lies in either flank of a sequence of length `len` —
/// within [`FLANK_SLOP`] of either end. Ported verbatim from production's `in_flank`.
fn in_flank(pos: usize, len: usize) -> bool {
    pos < FLANK_SLOP || pos + FLANK_SLOP >= len
}

/// The banded forward sum for **unequal-length** sequences, gaps confined to the flanks.
///
/// A faithful port of production's `banded_forward`
/// ([src/ssr/cohort/pair_hmm.rs](../../../ssr/cohort/pair_hmm.rs)). `read` has length `m`,
/// `reference` length `n`; the length difference is absorbed by gaps, but **only within
/// [`FLANK_SLOP`] of either end** — a gap touching an interior position is not admitted, so
/// an interior indel is *not competed* against an end gap and scores far below one. That
/// restriction is the whole point (spec §5.1): it is what keeps a sequencing error and a
/// stutter slip separately identifiable, and porting it — rather than writing a general
/// forward that would run fine and silently destroy that separation — is why this is its own
/// step.
///
/// The forward runs in **linear** probabilities: each cell sums its reachable predecessors
/// times the transition (a diagonal emission `1 − ε` / `ε/3`, or an end gap `ε`), and the
/// answer is the bottom-right cell. Cells outside the band `|i − j| ≤ |m − n| + FLANK_SLOP`
/// stay zero. The band floor is arithmetic — a global alignment must consume both sequences,
/// so the path is forced `|m − n|` off the diagonal — plus the flank slop as headroom.
fn banded_forward(
    read: &[u8],
    reference: &[u8],
    eps: f64,
    scratch: &mut SequenceMarginalScratch,
) -> f64 {
    let m = read.len();
    let n = reference.len();
    let band = m.abs_diff(n) + FLANK_SLOP;
    let gap = eps; // a flank-slop base is error-like under the flat emission, as in the source.
    let mismatch = eps / 3.0;
    let emit = |a: u8, b: u8| if a == b { 1.0 - eps } else { mismatch };

    let width = n + 1;
    scratch.cells.clear();
    scratch.cells.resize((m + 1) * width, 0.0);
    let at = |i: usize, j: usize| i * width + j;
    scratch.cells[at(0, 0)] = 1.0;

    for i in 0..=m {
        for j in 0..=n {
            if i == 0 && j == 0 {
                continue;
            }
            if i.abs_diff(j) > band {
                continue;
            }
            let mut acc = 0.0;
            if i > 0 && j > 0 {
                acc += scratch.cells[at(i - 1, j - 1)] * emit(read[i - 1], reference[j - 1]);
            }
            // Insertion: consume read[i-1] against a gap in the reference — flank only.
            if i > 0 && in_flank(i - 1, m) {
                acc += scratch.cells[at(i - 1, j)] * gap;
            }
            // Deletion: consume reference[j-1] against a gap in the read — flank only.
            if j > 0 && in_flank(j - 1, n) {
                acc += scratch.cells[at(i, j - 1)] * gap;
            }
            scratch.cells[at(i, j)] = acc;
        }
    }
    scratch.cells[at(m, n)]
}

impl MarginalAligner for SsrSequenceMarginal {
    type Scratch = SequenceMarginalScratch;
    /// `()` — the sequence-versus-sequence marginal needs nothing locus-specific: it compares
    /// two bare sequences under the fixed error rate, with no geometry (contrast algorithm 6,
    /// whose context is the repeat geometry).
    type Context<'a> = ();

    /// **The logarithm boundary — the one place the linear result becomes a [`LogProb`].**
    ///
    /// `linear_probability` computes the answer in ordinary probabilities, faithfully
    /// to production's `align_subst`; this returns its **natural logarithm**, and that single
    /// `.ln()` is the whole of the conversion. It is isolated in its own step because the
    /// failure mode is silent: a missing `.ln()` would return a linear probability typed as a
    /// `LogProb`, and a doubled one `ln(ln(x))` — both plausible numbers, neither an error.
    ///
    /// **An unreachable line-up returns `LogProb(f64::NEG_INFINITY)`, not zero.** The linear
    /// probability of an impossible line-up is `0.0`, and `0f64.ln()` is `-∞` — exactly the
    /// distinction the [`LogProb`] contract exists to preserve: `-∞` is a value a caller can
    /// see and act on, where a linear `0` is indistinguishable from underflow (arch §1, spec
    /// §7). No special-casing is needed; `.ln()` already maps `0` to `-∞`.
    ///
    /// # A known normalisation defect, reproduced knowingly
    ///
    /// Over **equal-length** inputs this is a proper distribution (the substitution products
    /// sum to one). The unequal-length case adds end-gap paths *on top* of that, so summed
    /// over all differing-length outputs the total mass **slightly exceeds one** (spec §5.1).
    /// It is inherited from `align_subst` deliberately, not fixed here: it is harmless where
    /// this is used — the genotyping normalises per read — and the equal-length case dominates
    /// overwhelmingly. Recorded so it is reproduced knowingly rather than discovered as a
    /// surprise.
    fn marginal_probability(
        &self,
        read: &[u8],
        reference: &[u8],
        _context: (),
        scratch: &mut Self::Scratch,
    ) -> LogProb {
        LogProb(self.linear_probability(read, reference, scratch).ln())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

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

    // ---- B2: the unequal-length banded forward -----------------------------------------
    //
    // Driven **deliberately** through `linear_probability` (the dispatcher), because the
    // caller algorithm 5 was written for — production's read-likelihood model — always makes
    // the two sequences equal-length first, so it NEVER reaches this branch (spec §5.1). A
    // test built by imitating that caller would leave the forward pass untouched.

    /// **The assertion that pins the interior-gap restriction — the point of the whole step.**
    /// An indel in the *middle* of the tract must score far below an end gap of the same
    /// size. It is not that the interior indel loses a close contest: an interior gap is not
    /// admitted at all, so it is never competed against an end gap — it can only be reached
    /// as a run of substitutions, which costs orders of magnitude more. A general forward
    /// that allowed interior gaps would score the two comparably and silently destroy the
    /// base-error-vs-stutter separability this restriction exists to keep.
    #[test]
    fn an_interior_indel_scores_far_below_an_end_gap() {
        let a = aligner(EPS);
        let mut scratch = SequenceMarginalScratch::new();
        let reference = b"CACACACACACA"; // 12 bp
        let end_indel = b"CACACACACACAC"; // one extra base at the trailing flank
        let mut interior_indel = reference.to_vec();
        interior_indel.insert(6, b'G'); // one extra base in the middle of the tract
        assert_eq!(end_indel.len(), reference.len() + 1);
        assert_eq!(interior_indel.len(), reference.len() + 1);

        let end = a.linear_probability(end_indel, reference, &mut scratch);
        let interior = a.linear_probability(&interior_indel, reference, &mut scratch);
        assert!(
            interior < end * 1e-3,
            "interior indel ({interior}) must score far below an end gap ({end})"
        );
    }

    /// A single end (flank) indel *is* scored, and tracks the clean end-gap path: the
    /// diagonal over all of the reference, then one gap at the flank, `(1 − ε)^len · ε`.
    #[test]
    fn a_single_end_indel_is_scored() {
        let a = aligner(EPS);
        let mut scratch = SequenceMarginalScratch::new();
        let reference = b"CACACACA";
        let read = b"CACACACAC"; // one extra base at the trailing flank
        let got = a.linear_probability(read, reference, &mut scratch);
        let clean_end_gap_path = (1.0 - EPS).powi(8) * EPS;
        assert!(got > 0.0);
        assert!(
            (got / clean_end_gap_path - 1.0).abs() < 0.05,
            "end-indel score {got} should track the clean end-gap path {clean_end_gap_path}"
        );
    }

    /// **The dispatcher routes by length — the structural guarantee behind B1's debug-only
    /// assertion.** An equal-length pair goes to the single-term sum (identical to calling
    /// `equal_length_probability` directly); an unequal-length pair goes to the forward and
    /// comes back a valid probability — crucially *without* reaching (and tripping) the
    /// equal-length arm's precondition. This is what actually keeps that arm off unequal
    /// input, which is why the length guarantee is tested here rather than on the assertion.
    #[test]
    fn linear_probability_routes_equal_and_unequal_lengths_to_the_right_arm() {
        let a = aligner(EPS);
        let mut scratch = SequenceMarginalScratch::new();

        // Equal length → exactly the single-term sum.
        let via_dispatch = a.linear_probability(b"CACACACA", b"CACAGACA", &mut scratch);
        assert_eq!(
            via_dispatch,
            a.equal_length_probability(b"CACACACA", b"CACAGACA")
        );

        // Unequal length → the forward, a valid probability, no panic.
        let unequal = a.linear_probability(b"CACACACAC", b"CACACACA", &mut scratch);
        assert!(unequal > 0.0 && unequal <= 1.0, "unequal score {unequal}");
    }

    /// The unequal-length forward returns a probability in `[0, 1]` — the boundary-slop
    /// normalisation defect lets the total mass over *differing-length* outputs slightly
    /// exceed one, but any single score is still a probability.
    #[test]
    fn a_longer_read_returns_a_probability_in_the_unit_interval() {
        let a = aligner(EPS);
        let mut scratch = SequenceMarginalScratch::new();
        let p = a.linear_probability(b"CACACAC", b"CACACA", &mut scratch);
        assert!((0.0..=1.0).contains(&p), "p = {p}");
    }

    /// A reused scratch buffer gives a bit-identical result to a fresh one, even after being
    /// primed by an unrelated, larger alignment — the buffer decides nothing, which is what
    /// the trait's `Scratch: Default` contract requires and what keeps the cohort's
    /// byte-identity guarantee intact.
    #[test]
    fn a_reused_scratch_gives_a_bit_identical_result() {
        let a = aligner(EPS);
        let mut primed = SequenceMarginalScratch::new();
        let mut fresh = SequenceMarginalScratch::new();
        let reference = b"CACACACA";
        let read = b"CACACACAC";
        // Prime with an unrelated, larger alignment, then reuse.
        let _ = a.linear_probability(b"CACACACACACA", b"CACA", &mut primed);
        let reused = a.linear_probability(read, reference, &mut primed);
        let clean = a.linear_probability(read, reference, &mut fresh);
        assert_eq!(reused.to_bits(), clean.to_bits());
    }

    /// `in_flank` marks exactly the [`FLANK_SLOP`] positions at each end and nothing in
    /// between — the boundary that admits an end gap but forbids an interior one.
    #[test]
    fn in_flank_marks_only_the_ends() {
        // len 8, FLANK_SLOP 2: positions 0,1 and 6,7 are flank; 2..=5 are interior.
        let flags: Vec<bool> = (0..8).map(|pos| in_flank(pos, 8)).collect();
        assert_eq!(
            flags,
            vec![true, true, false, false, false, false, true, true]
        );
    }

    /// **The forward is symmetric in its two sequences, which drives the deletion side.** The
    /// insertion transition (`in_flank(i-1, m)`) and the deletion transition
    /// (`in_flank(j-1, n)`) are mirror images: scoring the longer sequence as the read
    /// exercises the first, scoring it as the reference exercises the second, and the two must
    /// agree. Without this, a regression that mishandled the *deletion* branch — read shorter
    /// than reference — would slip past every other test, which only ever makes the read the
    /// longer or equal sequence. Approximate, not bit-exact: the two fills add the insertion
    /// and deletion terms in opposite order, so they round differently by a ULP.
    #[test]
    fn the_forward_is_symmetric_so_the_deletion_side_is_scored() {
        let a = aligner(EPS);
        let mut scratch = SequenceMarginalScratch::new();
        let longer = b"GGACGTACGT"; // 10 bp: the 8-bp sequence with two leading-flank bases
        let shorter = b"ACGTACGT"; // 8 bp
        let read_longer = a.linear_probability(longer, shorter, &mut scratch); // insertion side
        let read_shorter = a.linear_probability(shorter, longer, &mut scratch); // deletion side
        assert!(read_longer > 0.0 && read_shorter > 0.0);
        assert!(
            (read_longer - read_shorter).abs() < 1e-15,
            "symmetric: read-longer {read_longer} vs read-shorter {read_shorter}"
        );
    }

    /// **The band floor is `|m − n|`, not `FLANK_SLOP` — and that is load-bearing.** A length
    /// difference of 3 (beyond `FLANK_SLOP` = 2) is reachable only because the band widens by
    /// the length difference: two end gaps at the leading flank and one at the trailing flank
    /// give a valid line-up, and the final cell sits `|m − n| = 3` off the diagonal. Drop the
    /// `|m − n|` term from the band — a plausible "simplification" — and the answer cell falls
    /// outside a `FLANK_SLOP`-only band and the score collapses to `0.0`, silently losing every
    /// stutter event larger than the slop. All the other unequal-length tests use a difference
    /// of 1, so only this one would catch that.
    #[test]
    fn the_band_floor_admits_a_length_difference_beyond_the_flank_slop() {
        let a = aligner(EPS);
        let mut scratch = SequenceMarginalScratch::new();
        let reference = b"ACGTACGT"; // 8 bp
        let read = b"GGACGTACGTA"; // 11 bp: "GG" prepended, "A" appended — three end gaps
        assert_eq!(read.len(), reference.len() + 3);
        let got = a.linear_probability(read, reference, &mut scratch);
        assert!(
            got > 0.0 && got <= 1.0,
            "a length-3 difference must score above zero (band floor); got {got}"
        );
    }

    /// An empty read against a reference longer than `2 · FLANK_SLOP` scores exactly `0`: it
    /// could only be produced by deleting the whole reference, but interior deletions are
    /// forbidden, so no line-up reaches the end — the interior-gap restriction at its extreme.
    #[test]
    fn an_empty_read_cannot_delete_a_long_references_interior() {
        let a = aligner(EPS);
        let mut scratch = SequenceMarginalScratch::new();
        let p = a.linear_probability(b"", b"ACGTACGT", &mut scratch); // 0 vs 8 bp
        assert_eq!(p, 0.0);
    }

    // ---- B3: the logarithm boundary + parity against the ported function ----------------

    /// **Parity against the ported function — the byte-level oracle for this port.** ng's
    /// linear scorer must reproduce production's `align_subst`
    /// ([src/ssr/cohort/pair_hmm.rs](../../../ssr/cohort/pair_hmm.rs)) **bit for bit**, across
    /// all three of its arms and both gap directions. This is the fidelity the logarithm
    /// boundary rests on: match the linear values and the logged ones match too. Production is
    /// a read-only oracle here, exactly as `delimit_parity` uses `delimit_read`.
    #[test]
    fn linear_probability_matches_production_align_subst_bit_for_bit() {
        use crate::ssr::cohort::pair_hmm::{HmmScratch, align_subst};
        let cases: &[(&[u8], &[u8])] = &[
            (b"CACACACA", b"CACACACA"),    // exact match — the fast path
            (b"CACAGACA", b"CACACACA"),    // equal-length substitution
            (b"CACACACAC", b"CACACACA"),   // read longer — insertion side
            (b"CACACACA", b"CACACACAC"),   // read shorter — deletion side
            (b"GGACGTACGTA", b"ACGTACGT"), // a length-3 difference through the band floor
        ];
        // Sweep ε, endpoints included, so the degenerate 0.0 / 1.0 cases are checked against
        // production and not only against B1's hand-computed values.
        for &eps in &[0.0, 0.001, EPS, 0.2, 1.0] {
            let a = aligner(eps);
            for &(read, reference) in cases {
                let mut ours_scratch = SequenceMarginalScratch::new();
                let mut prod_scratch = HmmScratch::new();
                let ours = a.linear_probability(read, reference, &mut ours_scratch);
                let prod = align_subst(read, reference, eps, &mut prod_scratch);
                assert_eq!(
                    ours.to_bits(),
                    prod.to_bits(),
                    "linear parity failed at eps {eps} for {read:?} vs {reference:?}: \
                     ours {ours}, prod {prod}"
                );
            }
        }
    }

    /// **The whole boundary, end to end, as one contract:** the marginal equals production's
    /// `align_subst` *logarithm*. The two bit-exact tests above (linear parity, and the
    /// marginal is one `.ln()` of the linear) already imply this transitively, but stating it
    /// as a single assertion documents the contract the module actually offers a caller.
    #[test]
    fn the_marginal_is_the_logarithm_of_production_align_subst() {
        use crate::ssr::cohort::pair_hmm::{HmmScratch, align_subst};
        let a = aligner(EPS);
        for &(read, reference) in &[
            (b"CACACACA".as_slice(), b"CACACACA".as_slice()), // exact
            (b"CACAGACA".as_slice(), b"CACACACA".as_slice()), // substitution
            (b"CACACACAC".as_slice(), b"CACACACA".as_slice()), // insertion
            (b"CACACACA".as_slice(), b"CACACACAC".as_slice()), // deletion
        ] {
            let mut ours = SequenceMarginalScratch::new();
            let mut prod = HmmScratch::new();
            let marginal = a.marginal_probability(read, reference, (), &mut ours);
            let prod_log = align_subst(read, reference, EPS, &mut prod).ln();
            assert_eq!(marginal.get().to_bits(), prod_log.to_bits());
        }
    }

    /// **The boundary is exactly one logarithm — not zero, not two.** The marginal is the
    /// natural log of the linear probability. A missing `.ln()` would return a linear
    /// probability typed as a `LogProb`; a doubled one would return `ln(ln(x))` — both
    /// plausible wrong numbers, neither an error, which is why this pins the conversion
    /// directly against the linear value.
    #[test]
    fn marginal_probability_is_the_single_logarithm_of_the_linear_probability() {
        let a = aligner(EPS);
        let cases: &[(&[u8], &[u8])] = &[
            (b"CACACACA", b"CACACACA"),
            (b"CACAGACA", b"CACACACA"),
            (b"CACACACAC", b"CACACACA"),
        ];
        for &(read, reference) in cases {
            let mut s1 = SequenceMarginalScratch::new();
            let mut s2 = SequenceMarginalScratch::new();
            let marginal = a.marginal_probability(read, reference, (), &mut s1);
            let linear = a.linear_probability(read, reference, &mut s2);
            assert_eq!(marginal.get().to_bits(), linear.ln().to_bits());
        }
    }

    /// **An unreachable line-up is `-∞`, not `0` — the distinction the `LogProb` contract
    /// exists for.** An empty read against a reference too long to reach by end gaps has linear
    /// probability exactly `0`; the marginal is `ln(0) = -∞`, a value a caller can see and act
    /// on, rather than a linear `0` indistinguishable from underflow.
    #[test]
    fn an_unreachable_line_up_is_negative_infinity_not_zero() {
        let a = aligner(EPS);
        let mut scratch = SequenceMarginalScratch::new();
        let marginal = a.marginal_probability(b"", b"ACGTACGT", (), &mut scratch);
        // `== NEG_INFINITY` already excludes 0, NaN, and +∞ and pins the sign.
        assert_eq!(marginal.get(), f64::NEG_INFINITY);
    }

    proptest! {
        /// **The port matches production over the whole sequence domain, not just five hand-
        /// picked pairs.** For any two short byte strings and any ε in `[0, 1]`, ng's linear
        /// scorer must be bit-for-bit identical to production's `align_subst`. This is what
        /// guards the band and flank edges the example cases might all happen to agree on — a
        /// divergence here would be a real port defect, surfaced as a shrunk counterexample.
        #[test]
        fn linear_probability_matches_align_subst_on_random_sequences(
            read in prop::collection::vec(any::<u8>(), 0..12),
            reference in prop::collection::vec(any::<u8>(), 0..12),
            eps in 0.0f64..=1.0,
        ) {
            use crate::ssr::cohort::pair_hmm::{HmmScratch, align_subst};
            let a = SsrSequenceMarginal::try_new(eps).expect("eps in [0, 1]");
            let mut ours_scratch = SequenceMarginalScratch::new();
            let mut prod_scratch = HmmScratch::new();
            let ours = a.linear_probability(&read, &reference, &mut ours_scratch);
            let prod = align_subst(&read, &reference, eps, &mut prod_scratch);
            prop_assert_eq!(
                ours.to_bits(),
                prod.to_bits(),
                "eps {}, read {:?}, reference {:?}: ours {}, prod {}",
                eps, read, reference, ours, prod
            );
        }
    }

    /// **C2 — the sequence marginal scores its own generating allele highest.** A read whose
    /// *measured repeat* is a *k*-unit tract (with two substitution errors) scores higher
    /// against the *k*-unit candidate than against a neighbour (k±1, k±2): a length mismatch is
    /// absorbed only as end gaps and pays for it, and a mismatch of more than the flank slop is
    /// unreachable altogether. The ordering is **stable across the flat error rate**. This is
    /// what "computes what it claims" means for algorithm 5 *in isolation* — composing the
    /// stutter model that would properly price the length change is the genotyping's, not this
    /// module's (spec §5.1, §10.3).
    #[test]
    fn the_sequence_marginal_scores_its_generating_allele_highest() {
        let unit = b"CAG";
        for truth in [6usize, 10] {
            // The read's measured repeat: a truth-unit tract with two *real* substitution
            // errors — tract[1] is `A` (→ `T`) and tract[5] is `G` (→ `T`).
            let mut read = unit.repeat(truth);
            read[1] = b'T';
            read[5] = b'T';

            for eps in [0.001, 0.05] {
                let a = SsrSequenceMarginal::try_new(eps).expect("eps in [0, 1]");
                let mut scratch = SequenceMarginalScratch::new();
                let mut truth_score = f64::NEG_INFINITY;
                let mut best_neighbour = f64::NEG_INFINITY;
                for units in (truth - 2)..=(truth + 2) {
                    let candidate = unit.repeat(units);
                    let s = a
                        .marginal_probability(&read, &candidate, (), &mut scratch)
                        .get();
                    if units == truth {
                        truth_score = s;
                    } else {
                        best_neighbour = best_neighbour.max(s);
                    }
                }
                // The equal-length candidate needs no gaps; a neighbour pays end-gap penalties
                // (k±1) or is unreachable beyond the flank slop (k±2 → −∞). Either way the
                // truth wins by a clear margin, not a rounding-width one.
                assert!(
                    truth_score > best_neighbour + 1.0,
                    "at eps {eps}, truth={truth}: {truth_score} must beat the best neighbour \
                     {best_neighbour} by a clear margin"
                );
            }
        }
    }
}
