//! Algorithm 6 — the whole-read forward.
//!
//! Scores the **whole read** against a **full flank-repeat-flank** reference in one forward
//! pass, summed over every line-up, returned as a [`LogProb`]. This is the alternative to
//! algorithm 5 ([`ssr_marginal_sequence`](super::ssr_marginal_sequence)): where algorithm 5
//! scores an already-*measured* repeat against a candidate, this scores the read verbatim and
//! **marginalizes the repeat length itself**, so it inherits no delimitation and carries
//! forward no mistake a ruler might have made (spec §9). It has **no production counterpart** —
//! it is new code (arch §5) — built here to settle whether the whole-read design wins before
//! anyone pays for the per-read storage it would otherwise require. The head-to-head against
//! algorithm 5 is **not** in this module; it spans the genotyping that composes it (spec
//! §10.3). This file only proves the algorithm computes what it claims.
//!
//! ## A two-regime forward, reusing the delimiter's model
//!
//! It is the **sum-reduction analog of the tract delimiter** (algorithm 3): the same two gap
//! regimes — a **stiff** gap in the flanks and a **soft** gap inside the tract — but summed
//! over all line-ups (a forward), not maximised with a traceback (Viterbi). The two regimes
//! are what make it repeat-aware and what let a read score highest against **its own** allele:
//! the flanks anchor (stiff gaps hold the junctions), while the tract is free to change length
//! (soft gaps), so a read whose repeat is *k* units aligns cheaply to a *k*-unit reference and
//! dearly to any other. **This is the opposite of algorithm 5's rule:** algorithm 5 *forbids*
//! interior gaps (they are the stutter model's to explain); algorithm 6 *permits* soft tract
//! gaps precisely so it can marginalize the length without a prior measurement.
//!
//! Reuse, not rewrite (arch §5):
//!
//! - **Transitions** — production's [`TransitionCosts`], the exact two-regime gap model
//!   algorithm 3 uses (and algorithm 4 shares), read through its accessors. So algorithm 6
//!   inherits algorithm 3's *provisional* gap-open calibration verbatim — including the known
//!   match→match inconsistency, which computes `ln_match_to_match` from the flank open and
//!   never recomputes it under the tract regime (spec §4.2). Reproduced knowingly by reuse,
//!   not re-derived.
//! - **Emission** — the module's log-space [`FlatEmission`]. Algorithm 6 scores with **one
//!   flat error rate and a fixed synthetic quality** (spec §9): the flat emission ignores
//!   per-base quality, so every base is scored as if it carried one synthetic quality. That is
//!   the whole point — no real per-read quality is read or stored, which keeps the per-locus
//!   observation table a tally rather than per-read records. (Unlike algorithm 5, which needed
//!   a *linear* representation for byte-parity against `align_subst`, algorithm 6 has no parity
//!   oracle and runs in **log space**, so it reuses the log-space component directly — the
//!   "reuse the Emission component" the plan's preconditions asked for.)
//!
//! ## Log space, and unbanded
//!
//! The forward runs in **log probabilities** with a stable log-sum-exp, so the result is a
//! [`LogProb`] with no separate conversion step: a whole read is longer than a repeat, and log
//! space keeps it from underflowing however long the read (spec §7). It fills the **full**
//! `(m+1)×(n+1)` matrix for each of the three states — **unbanded**, matching how algorithm 3
//! first shipped: correctness before the banding optimisation, which is a later, separately
//! proved change (spec §5, §10.3).

use super::emission::{Emission, FlatEmission};
use super::ssr_best_path_flat_gap::TransitionCosts;
use super::{MarginalAligner, RepeatContext, RepeatGeometry};
use crate::ng::types::{BaseQual, LogProb};

/// The fixed synthetic quality every base is scored at. Its numeric value is **immaterial**:
/// [`FlatEmission`] ignores per-base quality entirely, and passing a single fixed quality is
/// how "one flat error rate, no per-read qualities" is expressed at the emission interface
/// (spec §9). Named rather than inlined so the deliberate synthesised-quality choice is
/// visible, not hidden behind a bare literal.
const SYNTHETIC_QUALITY: BaseQual = BaseQual(40);

/// The reused three-state forward buffer (match / insertion / deletion), each a flattened
/// `(m+1)×(n+1)` log-probability matrix — so scoring a locus's reads allocates once rather
/// than per read. Empty until first used and grown on demand; it decides nothing about a
/// result, as the trait's `Scratch` bound requires. **Unbanded**: the full matrix is filled
/// (memory is not banded), which is fine at repeat-plus-flank sizes.
#[derive(Debug, Default)]
pub struct WholeReadMarginalScratch {
    match_state: Vec<f64>,
    insertion_state: Vec<f64>,
    deletion_state: Vec<f64>,
}

impl WholeReadMarginalScratch {
    /// A fresh, empty buffer.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

/// The whole-read forward aligner (algorithm 6).
///
/// A stateless value carrying its two configuration pieces — the flat emission and the
/// two-regime transition costs — both fixed for the run, so it can be held by value and shared
/// across every read at a locus (arch §4). Everything that varies per locus (the repeat
/// geometry) arrives per call through the [`RepeatContext`].
#[derive(Debug, Clone, Copy)]
pub struct SsrWholeReadMarginal {
    emission: FlatEmission,
    transitions: TransitionCosts,
}

impl SsrWholeReadMarginal {
    /// Build the aligner from a flat emission model. The emission carries the sample group's
    /// error rate (already validated by [`FlatEmission::try_new`]); the two-regime transition
    /// costs are the fixed, reused delimiter calibration.
    #[must_use]
    pub fn new(emission: FlatEmission) -> Self {
        Self {
            emission,
            transitions: TransitionCosts::new(),
        }
    }

    /// The forward log-probability of `read` against the full flank-repeat-flank `reference`,
    /// summed over every line-up. **Already a logarithm** — the trait returns it directly.
    ///
    /// # Panics (debug only)
    ///
    /// Debug-asserts the geometry fits the reference (`left + right <= reference.len()`); the
    /// caller upholds this ([`RepeatGeometry::fits_reference`]), and a release build with a
    /// geometry that overran would misclassify the tract, not crash.
    fn forward_log(
        &self,
        read: &[u8],
        reference: &[u8],
        geometry: &RepeatGeometry,
        scratch: &mut WholeReadMarginalScratch,
    ) -> f64 {
        let m = read.len();
        let n = reference.len();
        // Saturating, matching the sibling delimiter's handling of the same `Bp(u64)` geometry:
        // on a 32-bit target a bare `as usize` would wrap, and a wrapped flank length would
        // misclassify the tract rather than fail. `u64 → usize` is lossless on 64-bit.
        let left_flank = usize::try_from(geometry.left_flank_len.get()).unwrap_or(usize::MAX);
        let right_flank = usize::try_from(geometry.right_flank_len.get()).unwrap_or(usize::MAX);
        debug_assert!(
            left_flank + right_flank <= n,
            "geometry (left {left_flank} + right {right_flank}) must fit the reference ({n})"
        );
        // Tract = reference columns [left_flank, n - right_flank). A gap touching a column in
        // this range is priced at the soft tract rate; every other column keeps the stiff flank
        // rate. `tract_end_column` is exclusive.
        let tract_end_column = n.saturating_sub(right_flank);

        // Hoist every per-cell constant out of the loop: the emission scores (constant for a
        // flat model), and the four transition logs.
        let scores = self.emission.scores_for(SYNTHETIC_QUALITY);
        let emit_insert = self.emission.insert_ln();
        let match_to_match = self.transitions.ln_match_to_match();
        let gap_extend = self.transitions.ln_gap_extend();
        let gap_close = self.transitions.ln_gap_close();
        let gap_open_flank = self.transitions.ln_gap_open();
        let gap_open_tract = self.transitions.ln_gap_open_tract();

        let width = n + 1;
        let size = (m + 1) * width;
        for state in [
            &mut scratch.match_state,
            &mut scratch.insertion_state,
            &mut scratch.deletion_state,
        ] {
            state.clear();
            state.resize(size, f64::NEG_INFINITY);
        }
        let at = |i: usize, j: usize| i * width + j;
        // The start: in the match state, nothing consumed, log-probability 0 (= ln 1).
        scratch.match_state[at(0, 0)] = 0.0;

        for i in 0..=m {
            for j in 0..=n {
                if i == 0 && j == 0 {
                    continue;
                }
                // A gap at column j is priced by the reference base it touches — the one most
                // recently consumed, `reference[j - 1]`; before any reference (j = 0) it is a
                // leading-flank gap.
                let gap_open = if j >= 1 && is_tract_column(j - 1, left_flank, tract_end_column) {
                    gap_open_tract
                } else {
                    gap_open_flank
                };

                // Match: consume read[i-1] against reference[j-1], from any state on the
                // diagonal predecessor.
                if i >= 1 && j >= 1 {
                    let emit = scores.pick(read[i - 1], reference[j - 1]);
                    scratch.match_state[at(i, j)] = emit
                        + ln_sum3(
                            match_to_match + scratch.match_state[at(i - 1, j - 1)],
                            gap_close + scratch.insertion_state[at(i - 1, j - 1)],
                            gap_close + scratch.deletion_state[at(i - 1, j - 1)],
                        );
                }

                // Insertion: consume read[i-1] against a gap in the reference (column stays j).
                if i >= 1 {
                    scratch.insertion_state[at(i, j)] = emit_insert
                        + ln_sum2(
                            gap_open + scratch.match_state[at(i - 1, j)],
                            gap_extend + scratch.insertion_state[at(i - 1, j)],
                        );
                }

                // Deletion: consume reference[j-1] against a gap in the read (row stays i). A
                // deletion emits nothing — there is no read base to score.
                if j >= 1 {
                    scratch.deletion_state[at(i, j)] = ln_sum2(
                        gap_open + scratch.match_state[at(i, j - 1)],
                        gap_extend + scratch.deletion_state[at(i, j - 1)],
                    );
                }
            }
        }

        // The whole read and the whole reference consumed, ending in any state.
        ln_sum3(
            scratch.match_state[at(m, n)],
            scratch.insertion_state[at(m, n)],
            scratch.deletion_state[at(m, n)],
        )
    }
}

/// Whether reference column `col` (0-based) lies inside the tract — at or after the last
/// left-flank base and before the first right-flank base, i.e. `left_flank <= col <
/// tract_end_column` (with `tract_end_column = n − right_flank`, exclusive). The two junctions
/// are load-bearing: column `left_flank − 1` is the last flank column (stiff gaps) and
/// `left_flank` the first tract column (soft gaps); one column later, `tract_end_column − 1` is
/// the last tract column and `tract_end_column` the first right-flank column. A free function so
/// the boundary can be tested directly, as the sibling delimiter tests its own tract window.
fn is_tract_column(col: usize, left_flank: usize, tract_end_column: usize) -> bool {
    col >= left_flank && col < tract_end_column
}

/// Add two log-probabilities: `ln(exp(a) + exp(b))`, numerically stable, with `-∞` (log 0)
/// absorbed rather than turned into `NaN`.
fn ln_sum2(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    hi + (lo - hi).exp().ln_1p()
}

/// Add three log-probabilities.
fn ln_sum3(a: f64, b: f64, c: f64) -> f64 {
    ln_sum2(ln_sum2(a, b), c)
}

impl MarginalAligner for SsrWholeReadMarginal {
    type Scratch = WholeReadMarginalScratch;
    /// [`RepeatContext`] — algorithm 6 needs the repeat's **geometry** (where the flanks end)
    /// to price gaps by region. It ignores the context's stutter model, exactly as the
    /// flat-gap delimiter does; length change is what this algorithm measures, not something a
    /// stutter model pre-judges.
    type Context<'a> = RepeatContext<'a>;

    /// The forward log-probability of the whole read against the full flank-repeat-flank
    /// reference. Already a logarithm (the forward runs in log space), so no boundary
    /// conversion; an impossible line-up sums to `-∞`, the value the [`LogProb`] contract
    /// preserves.
    fn marginal_probability(
        &self,
        read: &[u8],
        reference: &[u8],
        context: RepeatContext<'_>,
        scratch: &mut Self::Scratch,
    ) -> LogProb {
        LogProb(self.forward_log(read, reference, context.geometry, scratch))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::alignment::StutterModel;
    use crate::ng::types::{Bp, Motif};

    const EPS: f64 = 0.01;

    fn aligner() -> SsrWholeReadMarginal {
        SsrWholeReadMarginal::new(FlatEmission::try_new(EPS).expect("eps in [0, 1]"))
    }

    /// Build a flank-repeat-flank reference and its geometry from three pieces.
    fn locus(left_flank: &[u8], tract: &[u8], right_flank: &[u8]) -> (Vec<u8>, RepeatGeometry) {
        let mut reference = Vec::new();
        reference.extend_from_slice(left_flank);
        reference.extend_from_slice(tract);
        reference.extend_from_slice(right_flank);
        let geometry = RepeatGeometry {
            left_flank_len: Bp(left_flank.len() as u64),
            right_flank_len: Bp(right_flank.len() as u64),
            motif: Motif::new(b"CAG").expect("CAG is a valid period-3 motif"),
        };
        (reference, geometry)
    }

    /// Score a read against a locus, through the public trait entry point.
    fn score(read: &[u8], reference: &[u8], geometry: &RepeatGeometry) -> f64 {
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry,
            stutter: &stutter,
        };
        let mut scratch = WholeReadMarginalScratch::new();
        aligner()
            .marginal_probability(read, reference, context, &mut scratch)
            .get()
    }

    const LEFT: &[u8] = b"ACGTACGT";
    const RIGHT: &[u8] = b"TGCATGCA";
    const UNIT: &[u8] = b"CAG";

    fn tract(units: usize) -> Vec<u8> {
        UNIT.repeat(units)
    }

    /// A geometry with the given flank lengths (motif immaterial here — the forward reads only
    /// the flank boundaries).
    fn geometry(left: u64, right: u64) -> RepeatGeometry {
        RepeatGeometry {
            left_flank_len: Bp(left),
            right_flank_len: Bp(right),
            motif: Motif::new(b"CAG").expect("CAG is a valid period-3 motif"),
        }
    }

    /// **Absolute anchor #1 — the recurrence combines the right emission and transition.** A
    /// one-base read against a one-base flank reference has exactly one line-up, so the forward
    /// is the *hand-computed* `match_ln + ln_match_to_match`: a match emits the flat match score
    /// and pays exactly the match→match transition, nothing else. Every other test is
    /// relational; this pins an actual number, so a misplaced emission or a swapped transition
    /// that preserves ordering cannot hide.
    #[test]
    fn a_single_cell_forward_is_the_match_emission_times_the_match_transition() {
        let got = score(b"A", b"A", &geometry(1, 0));
        let costs = TransitionCosts::new();
        let scores = FlatEmission::try_new(EPS)
            .expect("eps in [0, 1]")
            .scores_for(SYNTHETIC_QUALITY);
        let expected = scores.match_ln + costs.ln_match_to_match();
        assert!(
            (got - expected).abs() < 1e-12,
            "got {got}, expected {expected}"
        );
    }

    /// **Absolute anchor #2 — it is a SUM, not a max.** A one-base read against a *two*-base
    /// flank reference lines up two ways: match the first base and delete the second, or delete
    /// the first and match the second. A forward SUMS them —
    /// `match_ln + gap_open + ln_sum2(gap_close, ln_match_to_match)` — whereas a Viterbi (or an
    /// accidental `f64::max` in the forward's place) returns only the larger single path,
    /// `match_ln + gap_open + ln_match_to_match`, about 0.49 nats lower. This is the one
    /// property that separates a marginal from a best path, and **nothing else in the suite
    /// pins it** — a `max`-based regression passes every relational test but fails here.
    #[test]
    fn a_two_path_forward_sums_the_paths_rather_than_maxing_them() {
        let got = score(b"A", b"AA", &geometry(2, 0));

        let costs = TransitionCosts::new();
        let match_ln = FlatEmission::try_new(EPS)
            .expect("eps in [0, 1]")
            .scores_for(SYNTHETIC_QUALITY)
            .match_ln;
        let open = costs.ln_gap_open(); // both columns are flank
        let expected = match_ln + open + ln_sum2(costs.ln_gap_close(), costs.ln_match_to_match());
        assert!(
            (got - expected).abs() < 1e-12,
            "got {got}, expected {expected}"
        );

        // Strictly above the best single path — the defining forward-vs-Viterbi gap (~0.49).
        let best_single_path = match_ln + open + costs.ln_match_to_match();
        assert!(
            got > best_single_path + 0.4,
            "a forward must sum the two paths (got {got}, best single {best_single_path})"
        );
    }

    /// **The tract window's two junctions, pinned directly.** With `left_flank = 2` and
    /// `tract_end_column = 4`, exactly columns 2 and 3 are tract; column 1 (last left flank) and
    /// column 4 (first right flank) are not. An off-by-one at either junction — the mistake that
    /// survived the sibling delimiter's whole suite once — fails here.
    #[test]
    fn the_tract_window_starts_at_the_left_flank_and_ends_before_the_right_flank() {
        let (left_flank, tract_end_column) = (2, 4);
        assert!(!is_tract_column(1, left_flank, tract_end_column)); // last left-flank column
        assert!(is_tract_column(2, left_flank, tract_end_column)); // first tract column
        assert!(is_tract_column(3, left_flank, tract_end_column)); // last tract column
        assert!(!is_tract_column(4, left_flank, tract_end_column)); // first right-flank column
    }

    /// Degenerate empty inputs have defined, finite answers rather than a panic or `NaN`. Two
    /// empty sequences line up the one trivial way (log-probability 0); a read against an empty
    /// reference is entirely inserted; and `n = 0` does not trip the indexing.
    #[test]
    fn degenerate_empty_sequences_have_defined_finite_answers() {
        assert_eq!(score(b"", b"", &geometry(0, 0)), 0.0);
        assert!(score(b"AC", b"", &geometry(0, 0)).is_finite());
    }

    /// An exact match — the read *is* the reference — scores higher than any altered read, and
    /// finitely (a real log-probability, not `-∞`).
    #[test]
    fn an_exact_match_scores_highest_and_finite() {
        let (reference, geometry) = locus(LEFT, &tract(4), RIGHT);
        let exact = score(&reference, &reference, &geometry);
        assert!(exact.is_finite(), "exact match must be finite, got {exact}");

        // One substitution in the tract scores strictly lower.
        let mut mutated = reference.clone();
        mutated[LEFT.len() + 1] = b'T'; // flip a tract base
        let substituted = score(&mutated, &reference, &geometry);
        assert!(
            substituted < exact,
            "a substitution ({substituted}) must score below the exact match ({exact})"
        );
    }

    /// **The tract is free to change length; the read stays reachable.** A read whose repeat is
    /// one unit longer than the reference still scores finitely — soft tract gaps absorb the
    /// difference. This is the defining difference from algorithm 5, which *forbids* interior
    /// gaps and would score such a read far below any end-gap path.
    #[test]
    fn an_expanded_tract_is_reachable() {
        let (reference, geometry) = locus(LEFT, &tract(4), RIGHT);
        let mut read = Vec::new();
        read.extend_from_slice(LEFT);
        read.extend_from_slice(&tract(5)); // one extra unit in the tract
        read.extend_from_slice(RIGHT);
        let expanded = score(&read, &reference, &geometry);
        assert!(
            expanded.is_finite(),
            "an expanded tract must be reachable, got {expanded}"
        );
        // And a contracted tract likewise.
        let mut short = Vec::new();
        short.extend_from_slice(LEFT);
        short.extend_from_slice(&tract(3));
        short.extend_from_slice(RIGHT);
        assert!(score(&short, &reference, &geometry).is_finite());
    }

    /// **The two regimes: the same-size length change is cheaper in the tract than in a
    /// flank.** Two reads, both three bases longer than the reference — one with the extra
    /// bases inside the tract (soft gap), one with them inside the left flank (stiff gap). The
    /// tract insertion must score higher, because the flank gap is ~350× stiffer per open. This
    /// is what "repeat-aware" means, and it is the property that makes a read score highest
    /// against its own allele.
    #[test]
    fn a_tract_length_change_scores_higher_than_a_flank_length_change() {
        let (reference, geometry) = locus(LEFT, &tract(4), RIGHT);

        // Extra unit inside the tract.
        let mut in_tract = Vec::new();
        in_tract.extend_from_slice(LEFT);
        in_tract.extend_from_slice(&tract(5));
        in_tract.extend_from_slice(RIGHT);

        // Three extra bases inside the left flank — same total length, different location.
        let mut in_flank = Vec::new();
        in_flank.extend_from_slice(&LEFT[..4]);
        in_flank.extend_from_slice(b"GGG"); // three inserted flank bases
        in_flank.extend_from_slice(&LEFT[4..]);
        in_flank.extend_from_slice(&tract(4));
        in_flank.extend_from_slice(RIGHT);
        assert_eq!(in_tract.len(), in_flank.len());

        let tract_score = score(&in_tract, &reference, &geometry);
        let flank_score = score(&in_flank, &reference, &geometry);
        assert!(
            tract_score > flank_score,
            "a tract length change ({tract_score}) must score above a flank one ({flank_score})"
        );
    }

    /// An empty read against a non-empty reference is **reachable** (unlike algorithm 5): every
    /// reference base is deleted, and deletions are permitted in both regimes. A defined, finite
    /// answer, and lower than a real match.
    #[test]
    fn an_empty_read_is_reachable_by_deleting_the_whole_reference() {
        let (reference, geometry) = locus(LEFT, &tract(4), RIGHT);
        let empty = score(b"", &reference, &geometry);
        assert!(
            empty.is_finite(),
            "an all-deletion path must be reachable, got {empty}"
        );
        assert!(empty < score(&reference, &reference, &geometry));
    }

    /// A reused scratch gives a bit-identical result to a fresh one, even after being primed by
    /// an unrelated, larger alignment — the buffer decides nothing, which the trait's
    /// `Scratch: Default` contract requires and the cohort's byte-identity rests on.
    #[test]
    fn a_reused_scratch_gives_a_bit_identical_result() {
        let (reference, geometry) = locus(LEFT, &tract(4), RIGHT);
        let read = &reference;
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry: &geometry,
            stutter: &stutter,
        };
        let a = aligner();

        let mut primed = WholeReadMarginalScratch::new();
        // Prime with an unrelated, larger alignment.
        let (big_ref, big_geo) = locus(LEFT, &tract(9), RIGHT);
        let big_stutter = StutterModel::hipstr_shipped();
        let big_context = RepeatContext {
            geometry: &big_geo,
            stutter: &big_stutter,
        };
        let _ = a.marginal_probability(&big_ref, &big_ref, big_context, &mut primed);
        let reused = a.marginal_probability(read, &reference, context, &mut primed);

        let mut fresh = WholeReadMarginalScratch::new();
        let clean = a.marginal_probability(read, &reference, context, &mut fresh);
        assert_eq!(reused.get().to_bits(), clean.get().to_bits());
    }

    /// A flankless locus (the whole reference is tract) is handled — every column is priced at
    /// the tract rate, and an exact match is still finite.
    #[test]
    fn a_flankless_locus_prices_every_column_as_tract() {
        let (reference, geometry) = locus(b"", &tract(4), b"");
        assert_eq!(geometry.left_flank_len, Bp(0));
        assert!(score(&reference, &reference, &geometry).is_finite());
    }

    /// `ln_sum2` is a stable log-add that absorbs `-∞` (log 0) on either side and never
    /// produces `NaN`, which the forward relies on to seed and terminate its matrices.
    #[test]
    fn ln_sum_absorbs_negative_infinity() {
        let ni = f64::NEG_INFINITY;
        assert_eq!(ln_sum2(ni, -1.0), -1.0);
        assert_eq!(ln_sum2(-1.0, ni), -1.0);
        assert_eq!(ln_sum2(ni, ni), ni);
        // ln(exp(0) + exp(0)) = ln 2.
        assert!((ln_sum2(0.0, 0.0) - 2.0f64.ln()).abs() < 1e-15);
        // Three-way, with one impossible term ignored.
        assert!((ln_sum3(0.0, 0.0, ni) - 2.0f64.ln()).abs() < 1e-15);
    }

    /// **C2 — the whole-read forward scores its own generating allele highest.** A read
    /// generated from a *k*-unit allele (flanks + *k* copies, with two substitution errors in
    /// the tract) must score higher against the *k*-unit reference than against any neighbour
    /// (k±1, k±2), and that ordering must be **stable across the flat error rate**. This is
    /// what "computes what it claims" means for algorithm 6 in isolation. It is deliberately
    /// **not** a comparison against algorithm 5 — that head-to-head scores a measured repeat
    /// against a whole read and spans the genotyping, not this module (spec §10.3).
    #[test]
    fn the_whole_read_forward_scores_its_generating_allele_highest() {
        for truth in [6usize, 10] {
            // A read from the truth allele, with two *real* substitution errors in the tract:
            // `CAG` repeats, so tract[1] is `A` (→ `T`) and tract[5] is `G` (→ `T`).
            let (mut read, _) = locus(LEFT, &tract(truth), RIGHT);
            read[LEFT.len() + 1] = b'T';
            read[LEFT.len() + 5] = b'T';

            for eps in [0.001, 0.05] {
                let a =
                    SsrWholeReadMarginal::new(FlatEmission::try_new(eps).expect("eps in [0, 1]"));
                let stutter = StutterModel::hipstr_shipped();
                let mut scratch = WholeReadMarginalScratch::new();
                let mut truth_score = f64::NEG_INFINITY;
                let mut best_neighbour = f64::NEG_INFINITY;
                for units in (truth - 2)..=(truth + 2) {
                    let (reference, geometry) = locus(LEFT, &tract(units), RIGHT);
                    let context = RepeatContext {
                        geometry: &geometry,
                        stutter: &stutter,
                    };
                    let s = a
                        .marginal_probability(&read, &reference, context, &mut scratch)
                        .get();
                    if units == truth {
                        truth_score = s;
                    } else {
                        best_neighbour = best_neighbour.max(s);
                    }
                }
                // A strict margin, not merely argmax: one unit is a soft tract gap, several
                // nats, so a bare "highest" that a rounding-width tie could satisfy is not
                // enough to claim the algorithm discriminates alleles.
                assert!(
                    truth_score > best_neighbour + 1.0,
                    "at eps {eps}, truth={truth}: {truth_score} must beat the best neighbour \
                     {best_neighbour} by a clear margin"
                );
            }
        }
    }
}
