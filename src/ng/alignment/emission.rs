//! Emission — scoring one read base against one reference base.
//!
//! This is a **component, not a variant**: per-base quality and a flat error rate are two
//! configurations of the same algorithm, which is what makes comparing them a swap rather
//! than a fork (spec §3). Every aligner in this module takes its emission model as a type
//! parameter, never behind `dyn` — this is called once per matrix cell, so a virtual call
//! per cell is not affordable (arch §4). The `Sized` supertrait on [`Emission`] makes that
//! decision a compile error rather than a convention.
//!
//! Two implementations land here: [`PerQualityEmission`], which reads the read's own
//! quality score, and [`FlatEmission`], which ignores it in favour of one rate. The first
//! is the default, on the reasoning that the read already carries a confidence for every
//! base and throwing it away to align the read would discard information we paid for
//! (spec §4.1); the second is the quality-blind end of the comparison.
//!
//! ## The scores must stay bit-equal to production
//!
//! [`PerQualityEmission`]'s table is ported from `EMISSION_LN`
//! ([src/ssr/pileup/alignment.rs](../../../ssr/pileup/alignment.rs)), and the repeat-aware
//! aligner built on it has to reproduce production's measured repeats **byte for byte** —
//! that parity is this module's only hard oracle (spec §10.3). So these are not merely
//! "close enough" numbers: reformulating the arithmetic in a way that moves the last bit
//! moves every downstream score. `per_quality_table_is_bit_exact` exists to make such a
//! change fail loudly instead of silently.

use crate::ng::types::BaseQual;
use std::sync::LazyLock;

/// What a base is worth at one quality: the score if it agrees with the reference, and
/// the score if it does not. Both are natural logarithms.
///
/// This exists so quality can be resolved **once per matrix row** rather than once per
/// cell — see [`Emission::scores_for`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaseScores {
    /// Score for a read base equal to the reference base, `ln(1 − ε)`.
    pub match_ln: f64,
    /// Score for a read base differing from the reference base, `ln(ε / 3)`.
    pub mismatch_ln: f64,
}

impl BaseScores {
    /// Pick the match or mismatch score by comparing the two bases.
    ///
    /// **Comparison is raw byte equality** — see [`Emission::scores_for`] for the
    /// precondition that puts on the caller.
    #[inline]
    #[must_use]
    pub fn pick(&self, read_base: u8, reference_base: u8) -> f64 {
        if read_base == reference_base {
            self.match_ln
        } else {
            self.mismatch_ln
        }
    }
}

/// Scores one read base against one reference base, in log space.
///
/// # Contract: pure and total
///
/// An implementation is a pure function of its arguments and its own constructor state —
/// no hidden mutation, so a score never depends on call order or thread count, which the
/// cohort's byte-identity guarantee rests on.
///
/// **Total** is the sharp end: *no* argument may produce a non-finite score. In
/// particular a quality-zero base must not score `-inf`. Quality zero means an error
/// probability of 1, so the naive `ln(1 − ε)` is `ln(0)`, and a single `-inf` annihilates
/// **every path** through that cell — the whole read becomes unalignable because of one
/// worthless base. Production floors it precisely to stop that (spec §4.2), and both
/// implementations here inherit the floor. `emission_is_finite_at_every_quality` pins it.
///
/// A `+inf` is worse still, and is the reason totality is stated unconditionally rather
/// than only for well-formed input: a positive log-probability makes the event it scores
/// *infinitely preferred*, so one would not produce a slightly wrong alignment but a
/// uniformly wrong one, with nothing in the output to show for it.
///
/// # Precondition: the caller canonicalizes the bases
///
/// Bases are compared by **raw byte equality**, exactly as production does. Two
/// consequences the caller owns, because nothing here can check them without paying for
/// it on every cell:
///
/// - **Case matters.** A soft-masked reference spells its bases in lower case, and `b'a'`
///   is not `b'A'` — a soft-masked stretch would score as a mismatch at every base.
///   Soft-masking marks repeats, which is precisely where the repeat-aware aligner works,
///   so this is not a remote hazard. ng offers both shapes deliberately:
///   [`RefSeq::fetch`] canonicalizes (upper-cases ACGT, folds everything else to `N`)
///   while [`RawRefSeq`] returns bases verbatim, soft-mask intact. **Feed this the
///   canonical form.**
/// - **`N` is a base like any other here.** An `N` read base against an `N` reference
///   base scores a full-confidence *match*. Neither this component nor production's
///   treats ambiguity specially; if that should change it is a scoring-model decision for
///   the aligner steps, not a silent fix here.
///
/// Both behaviours are pinned by `emission_compares_bases_by_raw_byte_equality` so they
/// are recorded rather than discovered.
///
/// [`RefSeq::fetch`]: crate::ng::ref_seq::RefSeq::fetch
/// [`RawRefSeq`]: crate::ng::ref_seq::RawRefSeq
pub trait Emission: Sized {
    /// Resolve the two scores for one quality, so a caller can hoist the lookup out of
    /// its inner loop.
    ///
    /// **This is the primary method, and the shape is deliberate.** Arch §6 left open
    /// "whether quality arrives per call or as a pre-resolved row", to be settled once
    /// two implementations existed; they now do, and the answer is per row. A quality
    /// belongs to a *read base*, so it is constant along a whole row of the alignment
    /// matrix while the reference base varies along it — production's own loop resolves
    /// it once per row and leaves a compare-and-select in the inner loop
    /// ([src/ssr/pileup/alignment.rs](../../../ssr/pileup/alignment.rs)). Billing a table
    /// lookup per *cell* would put work in the hottest loop in the module that the
    /// structure of the problem does not require.
    ///
    /// No measurement backs the choice — nothing is built yet to measure. It rests on
    /// matching the loop shape of the code this module must reproduce byte for byte.
    fn scores_for(&self, quality: BaseQual) -> BaseScores;

    /// Score an **inserted** read base, which has no reference base to score against.
    ///
    /// The two implementations have **no common source** for this value, so it is worth
    /// stating where each comes from. Production's per-quality path scores an inserted
    /// base against a uniform base composition, `ln(1/4)`. Its flat path has no emission
    /// for an inserted base at all — the value it uses there is a *transition* cost
    /// (`gap = eps`, [src/ssr/cohort/pair_hmm.rs](../../../ssr/cohort/pair_hmm.rs)), not
    /// an emission. So [`FlatEmission`]'s value is a **decision, not a port** (arch §2.3);
    /// see its documentation for the reasoning.
    fn insert_ln(&self) -> f64;

    /// Score `read_base` against `reference_base` at one quality — the convenience form
    /// of [`Self::scores_for`], for callers not in an inner loop.
    #[inline]
    fn emit_ln(&self, read_base: u8, reference_base: u8, quality: BaseQual) -> f64 {
        self.scores_for(quality).pick(read_base, reference_base)
    }
}

/// The smallest positive value a floored log may take the logarithm of. Using
/// [`f64::MIN_POSITIVE`] rather than an arbitrary epsilon keeps the floored probability as
/// close to the true zero as the type allows, so a floored base is *merely* wildly
/// improbable (about `-708` in log space) instead of impossible.
const PROBABILITY_FLOOR: f64 = f64::MIN_POSITIVE;

/// An inserted read base has no reference base to be compared with, so it is scored
/// against a uniform base composition: any of the four bases, `1/4`. Production's
/// per-quality path uses exactly this (`INS_EMIT_LN`).
///
/// A `const` rather than a lazily-computed static: this is read on the hot path, and the
/// module's whole argument against `dyn` is that per-cell overhead is not affordable — an
/// atomic initialisation check would be the same kind of cost. `uniform_base_ln_is_ln_of_a_quarter`
/// pins the literal against `0.25f64.ln()`.
const UNIFORM_BASE_LN: f64 = -1.386_294_361_119_890_6;

/// Per-Phred-quality emission scores, indexed by the raw quality byte — all 256 of them,
/// so no quality a BAM can hold can index out of range.
///
/// The Dindel base-quality model: with `ε = 10^(−Q/10)`, a match scores `ln(1 − ε)` and a
/// mismatch `ln(ε / 3)`, the error split evenly across the three other bases. Ported from
/// production's `EMISSION_LN`
/// ([src/ssr/pileup/alignment.rs](../../../ssr/pileup/alignment.rs)) **including its
/// quality-zero floor**, which is the whole reason this is a table and not a formula at
/// the call site.
///
/// **One deliberate difference from the source, which changes no value.** Production
/// floors the *match* term only; this floors the mismatch term too, so the trait's
/// totality contract holds by construction rather than by arithmetic accident. The
/// mismatch floor can never bind: the smallest `ε / 3` any `u8` quality produces is about
/// `1.05e-26` (at Q255), roughly 282 orders of magnitude above `f64::MIN_POSITIVE`. The
/// two tables are therefore bit-identical, which
/// `mismatch_floor_never_binds_over_the_quality_domain` asserts rather than assumes.
static PER_QUALITY_LN: LazyLock<[BaseScores; 256]> = LazyLock::new(|| {
    let mut table = [BaseScores {
        match_ln: 0.0,
        mismatch_ln: 0.0,
    }; 256];
    for (quality, entry) in table.iter_mut().enumerate() {
        let error_rate = 10f64.powf(-(quality as f64) / 10.0);
        *entry = BaseScores {
            match_ln: (1.0 - error_rate).max(PROBABILITY_FLOOR).ln(),
            mismatch_ln: (error_rate / 3.0).max(PROBABILITY_FLOOR).ln(),
        };
    }
    table
});

/// Scores each base with the read's own quality — the default, and production's model.
///
/// Stateless: the scores live in a shared table built once, so constructing one of these
/// is free and holding one costs nothing. That is what lets an aligner take it by value
/// as a type parameter.
///
/// `Default` is written out rather than derived, deliberately: a derived `Default` is the
/// one construction path that would keep compiling if this gained a field, silently
/// zero-filling it, whereas [`Self::new`] would fail to compile and say so.
#[derive(Debug, Clone, Copy)]
pub struct PerQualityEmission;

impl PerQualityEmission {
    /// Build the per-quality model. Free — the scores are a shared table.
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}

impl Default for PerQualityEmission {
    fn default() -> Self {
        Self::new()
    }
}

impl Emission for PerQualityEmission {
    #[inline]
    fn scores_for(&self, quality: BaseQual) -> BaseScores {
        PER_QUALITY_LN[quality.get() as usize]
    }

    #[inline]
    fn insert_ln(&self) -> f64 {
        UNIFORM_BASE_LN
    }
}

/// Scores every base at one fixed error rate, ignoring the read's qualities — the
/// quality-blind end of the emission comparison (spec §3).
///
/// The rate is the sample group's per-base substitution rate, so it is the experiment's
/// configuration rather than a per-read fact; that is why it is constructor state while
/// the repeat geometry travels per call (arch §4).
///
/// Both scores are computed **once, at construction**, because a flat model's scores do
/// not depend on anything the matrix varies — which also keeps the per-cell work to one
/// comparison and one field read. There is deliberately **no `Default`**: a default error
/// rate would be exactly the kind of behaviourally-significant hidden value that has to be
/// visible at the call site.
#[derive(Debug, Clone, Copy)]
pub struct FlatEmission {
    scores: BaseScores,
}

impl FlatEmission {
    /// Build a flat model from a per-base error rate, `ε`, matching production's flat
    /// path (`align_subst`,
    /// [src/ssr/cohort/pair_hmm.rs](../../../ssr/cohort/pair_hmm.rs)): a match scores
    /// `1 − ε` and a mismatch `ε / 3`, converted to log space here.
    ///
    /// # Preconditions, and why they are not a `Result`
    ///
    /// `error_rate` must be a finite probability in `[0, 1]`. Nothing in this module
    /// returns a `Result` — a fallible emission would push error handling onto the
    /// caller's hottest path for a value that is the run's fixed configuration, not
    /// per-read data (arch §3). The precondition is therefore a documented one plus a
    /// debug assertion, and it is worth being blunt about what that means: **a debug
    /// assertion compiles out of the release build this project actually runs.**
    ///
    /// **This rate is not user-typeable today** — the only callers are tests. Arch §3
    /// says that if such a condition turns out to be reachable from untrusted input it
    /// becomes a *checked constructor on the context type*, and this type is that context
    /// type. So the moment a CLI flag or config field feeds this, `new` should become a
    /// `try_new` returning a `DomainError` — the shape [`MismatchFraction`] already uses.
    ///
    /// # Panics
    ///
    /// Panics in debug builds if `error_rate` is not a finite value in `[0, 1]`. In
    /// release builds the assertion is compiled out and the value is clamped instead
    /// (below), so the call never panics there.
    ///
    /// # Why a release clamp as well
    ///
    /// The debug assertion is the real check; the clamp is what keeps the **totality
    /// contract** true once it is gone. Without it the floor is a lower bound only, so
    /// `ε = ∞` yields `ln(∞) = +∞` for every mismatch — a positive log-probability that
    /// makes mismatching infinitely attractive, and a negative rate yields a match score
    /// above zero. Those are not slightly-wrong answers that surface downstream; they
    /// silently invert the model. A non-finite rate is treated as `1`, the
    /// no-information end of the scale, so garbage in cannot produce confident output.
    ///
    /// The `PROBABILITY_FLOOR` at each endpoint is separate and binds only at `ε = 0`
    /// (which would make every mismatch impossible) and `ε = 1` (every match impossible).
    /// Within the documented range it changes nothing.
    ///
    /// [`MismatchFraction`]: crate::ng::types::MismatchFraction
    #[must_use]
    pub fn new(error_rate: f64) -> Self {
        debug_assert!(
            error_rate.is_finite() && (0.0..=1.0).contains(&error_rate),
            "flat emission error rate {error_rate} is not a finite probability in [0, 1]"
        );
        Self::from_unchecked_rate(error_rate)
    }

    /// Everything [`Self::new`] does **except** the debug assertion — that is, exactly
    /// what a release build runs.
    ///
    /// It exists to be testable. The totality contract has to hold once the assertion is
    /// compiled out, but a test cannot reach that path through `new`, because the
    /// assertion fires in the test profile. Testing the clamp through this function is
    /// the only way to assert the release behaviour without conditionally compiling the
    /// test out of the build anyone actually runs.
    fn from_unchecked_rate(error_rate: f64) -> Self {
        let error_rate = if error_rate.is_nan() {
            1.0
        } else {
            error_rate.clamp(0.0, 1.0)
        };
        Self {
            scores: BaseScores {
                match_ln: (1.0 - error_rate).max(PROBABILITY_FLOOR).ln(),
                mismatch_ln: (error_rate / 3.0).max(PROBABILITY_FLOOR).ln(),
            },
        }
    }
}

impl Emission for FlatEmission {
    #[inline]
    fn scores_for(&self, _quality: BaseQual) -> BaseScores {
        self.scores
    }

    /// `ln(1/4)` — the same value [`PerQualityEmission`] uses, and **a decision rather
    /// than a port** (arch §2.3).
    ///
    /// The reasoning is that an inserted base has no reference base in *either* model, so
    /// there is nothing for `ε` to describe: `ε` is the chance of misreading a base whose
    /// true identity the reference supplies, and an insertion has no such base. Scoring
    /// it against a uniform composition is the same modelling assumption in both, so the
    /// two agree here and differ only where they are meant to — matched bases.
    ///
    /// The alternative would be to follow production's flat path, whose unequal-length
    /// route charges `gap = eps` for a base absorbed in the flank slop. That is a
    /// **transition** cost, not an emission, and putting it in this slot would price the
    /// same event twice once the aligner's own gap model is applied on top.
    #[inline]
    fn insert_ln(&self) -> f64 {
        UNIFORM_BASE_LN
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Phred quality to error probability, spelled out independently of the table under
    /// test so the test cannot pass by sharing the implementation's mistake.
    fn error_probability(quality: u8) -> f64 {
        10f64.powf(-f64::from(quality) / 10.0)
    }

    /// **Bit-exactness over the whole quality domain.** The repeat-aware aligner built on
    /// this table has to reproduce production's measured repeats byte for byte, so a
    /// reformulation that is merely *close* — `ln_1p`, `powi`, splitting the division —
    /// would move every downstream score while a tolerance test waved it through.
    /// `assert_eq!` on `f64` is deliberate here for that reason. The expectation is
    /// re-derived from the published model rather than hardcoded, so it stays portable.
    #[test]
    fn per_quality_table_is_bit_exact() {
        let emission = PerQualityEmission::new();
        for quality in 0..=u8::MAX {
            let error_rate = error_probability(quality);
            let expected = BaseScores {
                match_ln: (1.0 - error_rate).max(f64::MIN_POSITIVE).ln(),
                mismatch_ln: (error_rate / 3.0).max(f64::MIN_POSITIVE).ln(),
            };
            assert_eq!(
                emission.scores_for(BaseQual(quality)),
                expected,
                "table diverged from the Dindel model at Q{quality}"
            );
        }
    }

    /// The published Dindel model at qualities where the arithmetic is checkable by hand:
    /// Q10 → ε = 0.1, Q20 → 0.01, Q30 → 0.001. Tolerance is appropriate here — these are
    /// decimal literals a human can verify, not the bit-level contract, which
    /// `per_quality_table_is_bit_exact` carries.
    #[test]
    fn per_quality_emission_reproduces_the_dindel_model() {
        let emission = PerQualityEmission::new();
        for (quality, error_rate) in [(10u8, 0.1f64), (20, 0.01), (30, 0.001)] {
            assert!(
                (emission.emit_ln(b'A', b'A', BaseQual(quality)) - (1.0 - error_rate).ln()).abs()
                    < 1e-12,
                "match score at Q{quality}"
            );
            assert!(
                (emission.emit_ln(b'A', b'C', BaseQual(quality)) - (error_rate / 3.0).ln()).abs()
                    < 1e-12,
                "mismatch score at Q{quality}"
            );
        }
    }

    /// **The floor, which is the reason this is a table.** At Q0 the error probability is
    /// 1, so an unfloored `ln(1 − ε)` is `ln(0) = -inf`, and one such cell annihilates
    /// every path through it — the read stops being alignable because of a single
    /// worthless base. This test fails if the floor is dropped.
    #[test]
    fn per_quality_emission_floors_the_quality_zero_match_instead_of_annihilating_it() {
        let emission = PerQualityEmission::new();
        let match_ln = emission.emit_ln(b'A', b'A', BaseQual(0));

        assert!(
            match_ln.is_finite(),
            "a Q0 match must be finite, was {match_ln}"
        );
        // Checked against the value independently, not against `PROBABILITY_FLOOR` —
        // otherwise editing the constant would move both sides of the assertion. The
        // floor is documented as the smallest positive normal f64, ~2.2250738585072014e-308.
        assert_eq!(PROBABILITY_FLOOR, 2.225_073_858_507_201_4e-308);
        assert_eq!(match_ln, 2.225_073_858_507_201_4e-308f64.ln());
        assert!((match_ln - -708.396_418_532_264_1).abs() < 1e-9);
        // A Q0 mismatch (ε/3 = 1/3) stays an ordinary, unfloored score.
        assert!((emission.emit_ln(b'A', b'C', BaseQual(0)) - (1.0f64 / 3.0).ln()).abs() < 1e-12);
    }

    /// The port floors the mismatch term where production does not. That is safe only
    /// because the floor cannot bind there — assert it rather than assume it, since the
    /// claim is what makes the two tables bit-identical.
    #[test]
    fn mismatch_floor_never_binds_over_the_quality_domain() {
        for quality in 0..=u8::MAX {
            let smallest = error_probability(quality) / 3.0;
            assert!(
                smallest > PROBABILITY_FLOOR,
                "the mismatch floor would bind at Q{quality} ({smallest:e})"
            );
        }
        // The tightest case, at the top of the domain, with a wide margin.
        assert!(error_probability(u8::MAX) / 3.0 > 1e-27);
    }

    /// Totality, over the whole domain rather than at sampled points: every one of the
    /// 256 quality bytes a BAM can hold must give a finite score for both outcomes, in
    /// both implementations. This is the contract the trait states.
    #[test]
    fn emission_is_finite_at_every_quality() {
        let per_quality = PerQualityEmission::new();
        let flat = FlatEmission::new(0.01);
        for quality in 0..=u8::MAX {
            for score in [
                per_quality.emit_ln(b'A', b'A', BaseQual(quality)),
                per_quality.emit_ln(b'A', b'C', BaseQual(quality)),
                flat.emit_ln(b'A', b'A', BaseQual(quality)),
                flat.emit_ln(b'A', b'C', BaseQual(quality)),
            ] {
                assert!(score.is_finite(), "non-finite score at Q{quality}");
                assert!(score <= 0.0, "a log-probability above zero at Q{quality}");
            }
        }
        assert!(per_quality.insert_ln().is_finite());
        assert!(flat.insert_ln().is_finite());
    }

    /// **Totality must survive a violated precondition**, because the debug assertion
    /// that guards it is compiled out of the release build this project runs. Each of
    /// these rates would, without the clamp, produce a score that inverts the model
    /// rather than merely skewing it: `+inf` mismatches, or a match score above zero.
    ///
    /// Goes through `from_unchecked_rate` rather than `new` because `new`'s debug
    /// assertion fires in the test profile — the path under test is precisely the one
    /// that runs when that assertion is gone.
    #[test]
    fn flat_emission_stays_total_for_rates_outside_the_contract() {
        for error_rate in [
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
            -0.5,
            1.5,
            f64::MAX,
        ] {
            let flat = FlatEmission::from_unchecked_rate(error_rate);
            for score in [
                flat.emit_ln(b'A', b'A', BaseQual(30)),
                flat.emit_ln(b'A', b'C', BaseQual(30)),
            ] {
                assert!(score.is_finite(), "non-finite score for rate {error_rate}");
                assert!(
                    score <= 0.0,
                    "log-probability above zero ({score}) for rate {error_rate}"
                );
            }
        }
    }

    /// A match must outscore a mismatch at every usable quality. If that ever inverts,
    /// the aligner prefers disagreeing with the reference and every result is nonsense.
    ///
    /// **Q0 and Q1 are genuine exceptions, and the boundary is asserted rather than
    /// assumed.** A match beats a mismatch exactly when `1 − ε > ε / 3`, i.e. when
    /// `ε < 0.75`, which is `Q > 1.249`. So the model says a mismatch is the *likelier*
    /// reading at Q0 and Q1 — correctly: at ε ≈ 0.79 the base is very nearly noise, and
    /// there are three ways to disagree against one way to agree. This is a property of
    /// the Dindel model being ported, not a defect, but it is exactly the kind of edge a
    /// caller would assume away, so both sides of the boundary are pinned here.
    #[test]
    fn a_match_outscores_a_mismatch_above_the_quality_one_crossover() {
        let emission = PerQualityEmission::new();

        for quality in 2..=u8::MAX {
            assert!(
                emission.emit_ln(b'A', b'A', BaseQual(quality))
                    > emission.emit_ln(b'A', b'C', BaseQual(quality)),
                "match did not outscore mismatch at Q{quality}"
            );
        }

        for quality in [0u8, 1] {
            assert!(
                emission.emit_ln(b'A', b'A', BaseQual(quality))
                    < emission.emit_ln(b'A', b'C', BaseQual(quality)),
                "expected the mismatch to win below the crossover, at Q{quality}"
            );
        }
    }

    /// The same ordering for the flat model. Its `ε` is a **constructor parameter**, so
    /// unlike the quality table the crossover at `ε = 0.75` is reachable by
    /// configuration — and without this test, swapping the two structurally identical
    /// lines that build `match_ln` and `mismatch_ln` would pass the whole suite.
    #[test]
    fn flat_emission_orders_match_above_mismatch_below_the_crossover() {
        for error_rate in [0.0, 0.001, 0.01, 0.1, 0.5, 0.74] {
            let flat = FlatEmission::new(error_rate);
            assert!(
                flat.emit_ln(b'A', b'A', BaseQual(30)) > flat.emit_ln(b'A', b'C', BaseQual(30)),
                "match did not outscore mismatch at ε = {error_rate}"
            );
        }
        // Above the crossover the order reverses, for the same reason it does at Q0/Q1.
        for error_rate in [0.76, 0.9, 1.0] {
            let flat = FlatEmission::new(error_rate);
            assert!(
                flat.emit_ln(b'A', b'A', BaseQual(30)) < flat.emit_ln(b'A', b'C', BaseQual(30)),
                "expected the mismatch to win at ε = {error_rate}"
            );
        }
    }

    #[test]
    fn flat_emission_scores_every_quality_alike() {
        let flat = FlatEmission::new(0.02);
        let at_low = flat.emit_ln(b'A', b'A', BaseQual(2));
        let at_high = flat.emit_ln(b'A', b'A', BaseQual(60));
        assert_eq!(at_low, at_high);
        assert!((at_low - 0.98f64.ln()).abs() < 1e-12);
        assert!((flat.emit_ln(b'A', b'G', BaseQual(2)) - (0.02f64 / 3.0).ln()).abs() < 1e-12);
    }

    /// The flat model's degenerate endpoints are floored rather than annihilating, the
    /// same discipline the quality table applies at Q0.
    #[test]
    fn flat_emission_floors_both_degenerate_error_rates() {
        let perfect = FlatEmission::new(0.0);
        assert_eq!(perfect.emit_ln(b'A', b'A', BaseQual(30)), 0.0); // ln(1)
        assert_eq!(
            perfect.emit_ln(b'A', b'C', BaseQual(30)),
            PROBABILITY_FLOOR.ln()
        );

        let hopeless = FlatEmission::new(1.0);
        assert_eq!(
            hopeless.emit_ln(b'A', b'A', BaseQual(30)),
            PROBABILITY_FLOOR.ln()
        );
        assert!(hopeless.emit_ln(b'A', b'C', BaseQual(30)).is_finite());
    }

    /// Both implementations score an inserted base against a uniform base composition.
    /// They agree deliberately — see `FlatEmission::insert_ln` for why this is a decision
    /// rather than two ports that happened to coincide.
    #[test]
    fn both_implementations_score_an_inserted_base_against_a_uniform_composition() {
        assert_eq!(PerQualityEmission::new().insert_ln(), UNIFORM_BASE_LN);
        assert_eq!(FlatEmission::new(0.01).insert_ln(), UNIFORM_BASE_LN);
    }

    /// The literal is written out so it costs nothing on the hot path; this is what keeps
    /// it honest.
    #[test]
    fn uniform_base_ln_is_ln_of_a_quarter() {
        assert_eq!(UNIFORM_BASE_LN, 0.25f64.ln());
    }

    /// **The caller's precondition, pinned so it is recorded rather than discovered.**
    /// Bases are compared by raw byte equality, so a soft-masked (lower-case) reference
    /// mismatches everywhere, and `N` against `N` is a full-confidence match. Both follow
    /// production; neither is fixed here, because doing so would be a scoring-model change
    /// smuggled into a component.
    #[test]
    fn emission_compares_bases_by_raw_byte_equality() {
        let emission = PerQualityEmission::new();
        let quality = BaseQual(30);
        let matched = emission.emit_ln(b'A', b'A', quality);
        let mismatched = emission.emit_ln(b'A', b'C', quality);

        // Soft-masked reference: same base, different case, scored as a mismatch.
        assert_eq!(emission.emit_ln(b'A', b'a', quality), mismatched);
        assert_eq!(emission.emit_ln(b'a', b'a', quality), matched);

        // `N` is not special: it matches itself at full confidence.
        assert_eq!(emission.emit_ln(b'N', b'N', quality), matched);
        assert_eq!(emission.emit_ln(b'N', b'A', quality), mismatched);
    }

    /// Emission depends on whether the bases agree, never on which bases they are — the
    /// model has no per-base composition beyond that.
    #[test]
    fn emission_depends_only_on_whether_the_bases_agree() {
        let emission = PerQualityEmission::new();
        let quality = BaseQual(25);
        for (read_base, reference_base) in [(b'A', b'A'), (b'C', b'C'), (b'G', b'G'), (b'T', b'T')]
        {
            assert_eq!(
                emission.emit_ln(read_base, reference_base, quality),
                emission.emit_ln(b'A', b'A', quality)
            );
        }
        for (read_base, reference_base) in [(b'A', b'C'), (b'C', b'A'), (b'G', b'T'), (b'T', b'N')]
        {
            assert_eq!(
                emission.emit_ln(read_base, reference_base, quality),
                emission.emit_ln(b'A', b'C', quality)
            );
        }
    }

    /// The row-resolved and per-call forms must agree — `emit_ln` is a provided method
    /// over `scores_for`, and a caller that hoists must get the same scores as one that
    /// does not.
    #[test]
    fn row_resolved_scores_agree_with_the_per_call_form() {
        let per_quality = PerQualityEmission::new();
        let flat = FlatEmission::new(0.01);
        for quality in 0..=u8::MAX {
            let quality = BaseQual(quality);
            let row = per_quality.scores_for(quality);
            assert_eq!(
                row.pick(b'A', b'A'),
                per_quality.emit_ln(b'A', b'A', quality)
            );
            assert_eq!(
                row.pick(b'A', b'C'),
                per_quality.emit_ln(b'A', b'C', quality)
            );

            let flat_row = flat.scores_for(quality);
            assert_eq!(flat_row.pick(b'A', b'A'), flat.emit_ln(b'A', b'A', quality));
            assert_eq!(flat_row.pick(b'A', b'C'), flat.emit_ln(b'A', b'C', quality));
        }
    }
}
