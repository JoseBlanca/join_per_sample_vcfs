//! The stutter model — how likely each length change is, for one repeat.
//!
//! **Not an alignment algorithm**, and it lives here anyway because *two* consumers share
//! it: the two-penalty best-path aligner in this module (spec §4.2), and the genotyping
//! likelihood outside it. One description, so the two cannot drift apart (spec §5.2). It
//! is HipSTR's model.
//!
//! Everything here is a **linear probability**, not a logarithm — unlike HipSTR's own
//! fields, which are logs (`log_equal_`), and matching production's `stutter_pmf`, which
//! also returns a linear value. The module carries no `_ln` names for that reason; a value
//! that ever does become a logarithm must say so in its name (the crate-wide convention).
//!
//! ## The two regimes are the model's defining structure
//!
//! A read's length change is either a whole number of repeat units or it is not, and the
//! two mean different things:
//!
//! - **in frame** — a whole number of units. This is slippage, the common event, and its
//!   size is measured in **units**.
//! - **out of frame** — not a whole number of units. This is a sequencing indel or an
//!   interruption, not slippage; it is rarer, it gets its **own** parameters, and its size
//!   is measured in **base pairs**.
//!
//! Each regime splits again by direction, because stutter is asymmetric: **losing units is
//! more common than gaining them**.
//!
//! **Which regime applies is decided by arithmetic alone** — is the change a multiple of
//! the period? — and **never by what was actually inserted**. So an insertion that happens
//! to be period-sized is treated as slippage whether or not its bases are the repeat unit.
//! The composition is caught downstream instead, as a base-error mismatch against the
//! re-tiled candidate. The mis-routing is worst at period 1, where it catches roughly three
//! of every four single-base insertions — which is why the algorithm-3-versus-4 comparison
//! (Milestone D) must score an indel *of the repeat's own base* separately from an indel of
//! a different base, or the two effects cancel in the average (spec §4.2).

use std::num::NonZeroU8;

/// Largest slip this model scores at all; anything past it is **zero**, so an implausibly
/// large change is not explained away as stutter — such a read falls to the genotyping's
/// outlier handling instead (spec §5.2).
///
/// Copied from production's `MAX_SLIP` rather than imported: ng is a from-scratch caller
/// that does not depend on production (owner, 2026-07-16). **The value must stay equal to
/// production's** while the two models are meant to agree, and
/// `the_copied_cutoff_still_equals_productions` is what enforces that rather than trusting
/// this sentence. Production's is a compile-time array bound described there as a
/// provisional choice awaiting recalibration on the simulator.
///
/// # The one constant, two scales — inherited deliberately
///
/// Production applies this single number to the **unit** count in the in-frame branch and
/// to the re-indexed **base-pair** count in the out-of-frame branch. Those are different
/// scales, and spec §5.2 says to decide about that rather than inherit it silently.
///
/// **The decision here is to reproduce production's behaviour**, and to write down what it
/// costs. Out-of-frame changes are cut off about **`period − 1` times sooner** in real
/// terms: the re-indexed size is `Δ − Δ/period`, so a cutoff of 10 re-indexed steps admits
/// roughly `10 · period/(period − 1)` base pairs, against `10 · period` in frame. At period
/// 4 that is about 13 bp out of frame against 40 bp in frame; **at period 2 the effect
/// vanishes** (about 20 against 20), and it grows with the period. Changing this is a
/// **behaviour change to a model the genotyping shares**, and this plan's rule is transcribe
/// first, change separately with its own evidence. Recorded as a follow-up, not fixed here.
pub const MAX_SLIP: u32 = 10;

/// Lower bound on a geometric success probability, and the floor under
/// [`StutterModel::equal`]. Production uses one constant for both, so this does too.
const GEOM_MIN: f64 = 0.01;
/// Upper bound on a geometric success probability. See [`GEOM_MIN`].
const GEOM_MAX: f64 = 0.99;

/// The six stutter rates a [`StutterModel`] is built from, **named** rather than positional.
///
/// The names are load-bearing, not decoration. As six positional `f64` arguments, an
/// `in_up`/`in_down` or `in_geom`/`out_geom` transposition is invisible at every call site
/// *and* in a test suite whose fixtures happen to give the two members of a pair the same
/// value — which is exactly what the first version of this module's tests did. HipSTR keeps
/// the two geometrics genuinely independent (spec §5.2), so the pairing matters.
///
/// `equal` is **not** here: it is derived from these four masses, and taking it as an input
/// would let it disagree with the values that define it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StutterRates {
    /// Probability of an in-frame expansion, of any size.
    pub in_up: f64,
    /// Probability of an in-frame contraction, of any size. Usually exceeds `in_up` —
    /// stutter is contraction-biased.
    pub in_down: f64,
    /// How fast in-frame slip size decays, per **unit**.
    ///
    /// A geometric **success** probability, *not* a decay rate — see [`StutterModel::new`]
    /// for the trap.
    pub in_geom: f64,
    /// Probability of an out-of-frame expansion, of any size.
    pub out_up: f64,
    /// Probability of an out-of-frame contraction, of any size.
    pub out_down: f64,
    /// How fast out-of-frame size decays, per **base pair**. A geometric **success**
    /// probability. HipSTR keeps this independent of [`Self::in_geom`]; production ties the
    /// two to a single value, which is an undeclared placeholder rather than a fitted
    /// result (spec §5.2).
    pub out_geom: f64,
}

/// How likely each length change is, for one repeat.
///
/// # How it is built: per locus, from the *reference* allele
///
/// This model belongs to a **locus**, and is built from that locus's stutter shape and its
/// **reference** allele length — never from a candidate allele's length. The distinction is
/// not pedantic: stutter rises with allele length, so a per-*candidate* slip level is a real
/// and useful thing, but it belongs to the **genotyping**, which scores reads against
/// candidates. This module *measures*; letting the measurement's own model vary with the
/// candidate being tested would bias the ruler toward the answer under test (arch §5, §2.4).
///
/// Production derives its parameters per call from a per-locus stutter shape plus a
/// per-read stutter level. **ng has no stutter-shape type yet** — fitting one is the
/// genotyping's job (spec §5.2) — so this step lands the parameters and the distribution,
/// and the adapting constructor is owed once that type exists.
///
/// # Contract: `equal` is floored, so **do not test that the five sum to one**
///
/// The unchanged mass is *defined* as whatever the four direction masses leave,
/// `equal = 1 − in_up − in_down − out_up − out_down` — **but it is floored**. When the floor
/// binds, the five values sum to slightly more than one. That is deliberate: the floor is
/// what stops a hostile parameter combination producing a negative probability. It means
/// "the five sum to one" must **not** be written as a test, and
/// `the_five_masses_do_not_sum_to_one_when_the_floor_binds` exists to make that explicit
/// rather than merely absent.
///
/// (HipSTR instead *asserts* the masses sum below one at construction and never clamps.
/// Both disciplines are defensible; they differ, and only the floor matches the code being
/// ported — arch §2.4.)
///
/// # Why the fields are private
///
/// The clamps **are** the contract. Arch §2.4 sketches these as public fields, but public
/// fields would let a caller assemble a model whose geometrics sit at 0 or 1, or whose
/// `equal` is negative, and every guarantee above would be a comment rather than a fact.
/// Construction goes through [`StutterModel::new`], which applies them; the seven values are
/// readable through accessors. (Arch states its signatures are illustrative and the contract
/// is the deliverable.)
///
/// Not `Copy`: at seven `f64`s this is 56 bytes, well past the size where an implicit copy
/// is free, and the per-call repeat context holds it **by reference** by design — so `Copy`
/// would buy nothing and hide the cost of the cases it did serve.
#[derive(Debug, Clone, PartialEq)]
pub struct StutterModel {
    equal: f64,
    in_up: f64,
    in_down: f64,
    in_geom: f64,
    out_up: f64,
    out_down: f64,
    out_geom: f64,
}

impl StutterModel {
    /// Build a model from the six rates. `equal` is derived, as one minus the four direction
    /// masses, then floored.
    ///
    /// # The trap this type exists to make hard
    ///
    /// **`in_geom` and `out_geom` are geometric *success* probabilities, not decay rates.**
    /// A geometric with success probability `g` puts mass `g` on the first step and
    /// multiplies by `(1 − g)` for each step after, so a **larger** `g` concentrates mass on
    /// *single*-unit slips — HipSTR ships 0.95, meaning nineteen slips in twenty are exactly
    /// one unit. If a parameter arrives expressed as the probability of *continuing* to the
    /// next step (mean size `1/(1 − decay)`), it is the **complement**: `geom = 1 − decay`.
    /// Getting this backwards inverts the size distribution — large slips become common —
    /// and **nothing crashes** (spec §5.2, trap 1).
    /// `a_single_unit_slip_outweighs_a_larger_one` fails if it is inverted.
    ///
    /// # Panics
    ///
    /// Panics in debug builds if any rate is not a finite value in `[0, 1]`. In release the
    /// assertion is compiled out and the values are sanitized instead (below), so the call
    /// never panics there.
    ///
    /// # What sanitizing does, and why it is not the same thing as the clamps
    ///
    /// Two different jobs, easy to conflate. The **clamps on the geometrics and the floor
    /// under `equal` are ported model behaviour** — production does exactly this — and they
    /// bind only at degenerate values. The **mass sanitizing is a release-mode backstop** for
    /// a violated precondition, and it exists because without it this type would break its
    /// own stated contract: an unchecked mass of `2.0` makes `probability` return `1.9`, a
    /// negative mass returns a negative probability, and `NaN` poisons every score while
    /// *reporting a healthy floor* — because `f64::max` absorbs `NaN` and `f64::clamp` passes
    /// it straight through, so the constructor would look fine and only the results would be
    /// wrong. A non-finite rate becomes `0` (no stutter), the no-information end of the scale.
    #[must_use]
    pub fn new(rates: StutterRates) -> Self {
        debug_assert!(
            [
                rates.in_up,
                rates.in_down,
                rates.in_geom,
                rates.out_up,
                rates.out_down,
                rates.out_geom
            ]
            .iter()
            .all(|rate| rate.is_finite() && (0.0..=1.0).contains(rate)),
            "stutter rates must be finite probabilities in [0, 1]: {rates:?}"
        );
        Self::sanitized(rates)
    }

    /// Everything [`Self::new`] does **except** the debug assertion — that is, exactly what a
    /// release build runs.
    ///
    /// Separated so the release path is reachable from a test: the assertion fires in the
    /// test profile, so `new` cannot be used to check what happens once it is gone.
    fn sanitized(rates: StutterRates) -> Self {
        let in_up = sanitize_mass(rates.in_up);
        let in_down = sanitize_mass(rates.in_down);
        let out_up = sanitize_mass(rates.out_up);
        let out_down = sanitize_mass(rates.out_down);

        Self {
            equal: (1.0 - in_up - in_down - out_up - out_down).clamp(GEOM_MIN, 1.0),
            in_up,
            in_down,
            in_geom: sanitize_geom(rates.in_geom),
            out_up,
            out_down,
            out_geom: sanitize_geom(rates.out_geom),
        }
    }

    /// HipSTR's **shipped** default parameters, as a matched set.
    ///
    /// HipSTR has **two** parameter sets, and mixing them yields a pairing that exists
    /// nowhere — spec §5.2 records that an earlier draft of the spec did exactly that. These
    /// constructors exist so the two rows cannot be crossed by hand: this one is in-frame
    /// geom 0.95 with 0.05/0.05, out-of-frame geom 0.95 with 0.01/0.01.
    ///
    /// Note the shipped row makes expansion and contraction **equal**. That symmetry is a
    /// starting point, not a claim — HipSTR's *fitted* values are contraction-biased.
    #[must_use]
    pub fn hipstr_shipped() -> Self {
        Self::new(StutterRates {
            in_up: 0.05,
            in_down: 0.05,
            in_geom: 0.95,
            out_up: 0.01,
            out_down: 0.01,
            out_geom: 0.95,
        })
    }

    /// HipSTR's EM **starting point** — an initialisation its fitting immediately moves away
    /// from, and **not** a default. Provided as a matched set for the same reason as
    /// [`Self::hipstr_shipped`]: in-frame geom 0.9 with 0.1/0.1, out-of-frame geom 0.8 with
    /// 0.01/0.01.
    ///
    /// Note this row keeps the two geometrics **different** (0.9 against 0.8), which is the
    /// shape HipSTR actually uses — production ties them to one value, an undeclared
    /// placeholder (spec §5.2).
    #[must_use]
    pub fn hipstr_em_start() -> Self {
        Self::new(StutterRates {
            in_up: 0.1,
            in_down: 0.1,
            in_geom: 0.9,
            out_up: 0.01,
            out_down: 0.01,
            out_geom: 0.8,
        })
    }

    /// `P(length change)` for a change of `bp_diff` bases on a repeat of period `period`.
    ///
    /// Both regimes are a direction mass times a geometric over size, following spec §5.2
    /// term by term. Zero beyond [`MAX_SLIP`], so an implausibly large slip is not explained
    /// away.
    ///
    /// # Why out-of-frame sizes are re-indexed
    ///
    /// The out-of-frame geometric is indexed by `e = Δ − Δ/period` (truncated division), not
    /// by `Δ`. The reason is **not** double-counting — the two regimes are disjoint by
    /// construction, since a change is out-of-frame precisely when it is *not* a multiple of
    /// the period, so no length can reach both. What the re-indexing does is **compress the
    /// ranks**: it maps the out-of-frame values onto consecutive integers so the geometric's
    /// support has no gaps. At period 3 the out-of-frame values 1, 2, 4, 5, 7 map to
    /// e = 1, 2, 3, 4, 5. Without it the geometric would be evaluated at indices that skip
    /// the multiples, distorting the distribution.
    ///
    /// # Why `period` is a [`NonZeroU8`]
    ///
    /// The in-frame test divides by the period, so a zero would divide by zero — and a
    /// period-zero repeat is not a repeat in the first place. A `debug_assert!` would be the
    /// module's usual answer (arch §3), but here the illegal state is cheap to make
    /// **unrepresentable**, which is better than making it testable: there is no release
    /// path to get wrong, and no guard that compiles out. [`Motif`] already guarantees a
    /// period in `1..=6`, so the conversion at a real call site cannot fail.
    ///
    /// [`Motif`]: crate::ng::types::Motif
    /// [`NonZeroU8`]: std::num::NonZeroU8
    #[must_use]
    pub fn probability(&self, bp_diff: i64, period: NonZeroU8) -> f64 {
        let period = i64::from(period.get());

        if bp_diff % period == 0 {
            // In frame: slippage, sized in whole units.
            let units = bp_diff / period;
            if units == 0 {
                self.equal
            } else {
                regime(self.in_up, self.in_down, self.in_geom, units)
            }
        } else {
            // Out of frame: a sequencing indel or an interruption, sized in base pairs and
            // rank-compressed (see above). The size is never zero here — `bp_diff` is not a
            // multiple of the period, so truncated division cannot consume all of it.
            regime(
                self.out_up,
                self.out_down,
                self.out_geom,
                bp_diff - bp_diff / period,
            )
        }
    }

    /// Probability the read shows the allele's length unchanged. **Floored** — see the
    /// type's contract.
    #[must_use]
    pub fn equal(&self) -> f64 {
        self.equal
    }
    /// Probability of an in-frame expansion, of any size.
    #[must_use]
    pub fn in_up(&self) -> f64 {
        self.in_up
    }
    /// Probability of an in-frame contraction, of any size. Usually exceeds
    /// [`Self::in_up`] — stutter is contraction-biased.
    #[must_use]
    pub fn in_down(&self) -> f64 {
        self.in_down
    }
    /// How fast in-frame slip size decays, per **unit**. A geometric *success* probability.
    #[must_use]
    pub fn in_geom(&self) -> f64 {
        self.in_geom
    }
    /// Probability of an out-of-frame expansion, of any size.
    #[must_use]
    pub fn out_up(&self) -> f64 {
        self.out_up
    }
    /// Probability of an out-of-frame contraction, of any size.
    #[must_use]
    pub fn out_down(&self) -> f64 {
        self.out_down
    }
    /// How fast out-of-frame size decays, per **base pair**. A geometric *success*
    /// probability.
    #[must_use]
    pub fn out_geom(&self) -> f64 {
        self.out_geom
    }
}

/// One regime's whole answer: pick the direction by the sign of `steps`, drop anything past
/// the cutoff, and evaluate the geometric.
///
/// Written once and shared by both regimes because they are the same shape over different
/// parameters — which also makes the "size is at least one" precondition provable in one
/// place. `steps` is a **unit** count in frame and a re-indexed **base-pair** count out of
/// frame; [`MAX_SLIP`] is applied to whichever it is, which is the two-scale inheritance
/// recorded on that constant.
#[inline]
fn regime(up: f64, down: f64, geom: f64, steps: i64) -> f64 {
    debug_assert!(
        steps != 0,
        "a zero-size change is the `equal` case, not a regime"
    );
    let size = steps.unsigned_abs();
    if size > u64::from(MAX_SLIP) {
        return 0.0;
    }
    let mass = if steps < 0 { down } else { up };
    // `size >= 1` here, so the exponent cannot underflow. `unsigned_abs` also means
    // `i64::MIN` is safe, which a `-steps - 1` form would not be.
    mass * geom * (1.0 - geom).powi((size - 1) as i32)
}

/// A direction mass, made safe for release: `NaN` becomes zero (no stutter — the
/// no-information end), and anything outside `[0, 1]` is clamped into it.
#[inline]
fn sanitize_mass(mass: f64) -> f64 {
    if mass.is_nan() {
        0.0
    } else {
        mass.clamp(0.0, 1.0)
    }
}

/// A geometric success probability, held strictly inside `(0, 1)`. `NaN` becomes
/// [`GEOM_MIN`]; `f64::clamp` would otherwise pass it straight through.
#[inline]
fn sanitize_geom(geom: f64) -> f64 {
    if geom.is_nan() {
        GEOM_MIN
    } else {
        geom.clamp(GEOM_MIN, GEOM_MAX)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A repeat period, which the type system requires to be non-zero.
    fn period(bases: u8) -> NonZeroU8 {
        NonZeroU8::new(bases).expect("a test period is never zero")
    }

    /// A fixture whose six rates are **all different**, so a transposition of any pair
    /// changes an answer.
    ///
    /// This exists because the first version of these tests could not see one. Every fixture
    /// then used 0.05/0.05 for the in-frame masses, 0.01/0.01 for the out-of-frame ones, and
    /// 0.95 for *both* geometrics — so swapping `in_up` with `in_down`, or `in_geom` with
    /// `out_geom`, passed all twelve tests byte-identically. Contraction-biased in both
    /// regimes, as HipSTR's fitted values are.
    fn all_distinct() -> StutterModel {
        StutterModel::new(StutterRates {
            in_up: 0.03,
            in_down: 0.07,
            in_geom: 0.95,
            out_up: 0.004,
            out_down: 0.012,
            out_geom: 0.8,
        })
    }

    /// The published in-frame formula, term by term: `mass · geom · (1 − geom)^(n − 1)`,
    /// in **both** directions with different masses.
    #[test]
    fn the_in_frame_branch_reproduces_the_published_formula() {
        let model = all_distinct();
        let period = period(3);
        for units in 1..=5i64 {
            let decay = 0.05f64.powi((units - 1) as i32);
            let bp = units * i64::from(period.get());
            assert!((model.probability(bp, period) - 0.03 * 0.95 * decay).abs() < 1e-15);
            assert!((model.probability(-bp, period) - 0.07 * 0.95 * decay).abs() < 1e-15);
        }
    }

    /// The published out-of-frame formula, on the **re-indexed** size `e = Δ − Δ/period`,
    /// in **both** directions.
    ///
    /// The direction half is the point: with `out_up == out_down` — as every fixture in the
    /// first version of this file had — swapping the two out-of-frame masses passed the
    /// entire suite, because no assertion anywhere used a negative out-of-frame change.
    #[test]
    fn the_out_of_frame_branch_reproduces_the_published_formula_in_both_directions() {
        let model = all_distinct();
        let period = period(3);
        for bp_diff in [1i64, 2, 4, 5, 7] {
            let effective = bp_diff - bp_diff / i64::from(period.get());
            let decay = 0.2f64.powi((effective - 1) as i32);
            assert!((model.probability(bp_diff, period) - 0.004 * 0.8 * decay).abs() < 1e-15);
            assert!((model.probability(-bp_diff, period) - 0.012 * 0.8 * decay).abs() < 1e-15);
        }
    }

    /// **Direction asymmetry must be expressible in both regimes.** Stutter is
    /// contraction-biased, and a direction-symmetric model would be a step backwards from
    /// production's scoring, which already carries this (spec §4.2).
    #[test]
    fn a_contraction_outscores_an_expansion_of_the_same_size_in_both_regimes() {
        let model = all_distinct();
        let period = period(4);
        let unit = i64::from(period.get());

        // In frame.
        assert!(model.probability(-unit, period) > model.probability(unit, period));
        // Out of frame — the case nothing used to cover.
        assert!(model.probability(-1, period) > model.probability(1, period));
    }

    /// **Rank compression, the reason for the re-indexing.** At period 3 the out-of-frame
    /// changes 1, 2, 4, 5, 7 must map onto consecutive geometric steps 1, 2, 3, 4, 5 — so
    /// the support has no gaps. Indexing by Δ itself would skip the multiples and distort
    /// the distribution. Negatives mirror exactly, because Rust truncates toward zero.
    #[test]
    fn out_of_frame_sizes_compress_onto_consecutive_ranks() {
        let period = 3i64;
        let ranks: Vec<i64> = [1i64, 2, 4, 5, 7]
            .iter()
            .map(|&bp_diff| bp_diff - bp_diff / period)
            .collect();
        assert_eq!(ranks, vec![1, 2, 3, 4, 5]);

        let negative_ranks: Vec<i64> = [-1i64, -2, -4, -5, -7]
            .iter()
            .map(|&bp_diff| bp_diff - bp_diff / period)
            .collect();
        assert_eq!(negative_ranks, vec![-1, -2, -3, -4, -5]);
    }

    /// **The inverted-geometric trap** (spec §5.2, trap 1). A geometric *success*
    /// probability of 0.95 means nineteen slips in twenty are exactly one unit, so a
    /// one-unit slip must outweigh a two-unit slip by exactly `1/(1 − geom) = 20`. Read as a
    /// *decay* instead — the complement, 0.05 — the distribution inverts and large slips
    /// become common, and nothing crashes.
    ///
    /// **The quantitative assertion is the one doing the work.** Monotonicity alone would
    /// *not* catch the inversion: at geom = 0.05 the sequence still decreases, just far more
    /// slowly. The ratio is what pins it.
    #[test]
    fn a_single_unit_slip_outweighs_a_larger_one() {
        let model = StutterModel::hipstr_shipped();
        let period = period(4);
        let unit = i64::from(period.get());
        let one_unit = model.probability(unit, period);
        let two_units = model.probability(2 * unit, period);
        let three_units = model.probability(3 * unit, period);

        assert!(one_unit > two_units && two_units > three_units);
        assert!((one_unit / two_units - 20.0).abs() < 1e-9);
        assert!((two_units / three_units - 20.0).abs() < 1e-9);
    }

    /// **The `equal` floor, and why "the five sum to one" must not be a test** (arch §2.4).
    /// With four direction masses summing past 1, the derived unchanged mass would go
    /// negative; the floor stops that, and the consequence is that the five values then sum
    /// to *more* than one. Asserted here so the absence of a sums-to-one test reads as a
    /// decision rather than an oversight.
    #[test]
    fn the_five_masses_do_not_sum_to_one_when_the_floor_binds() {
        let hostile = StutterModel::new(StutterRates {
            in_up: 0.5,
            in_down: 0.5,
            in_geom: 0.95,
            out_up: 0.1,
            out_down: 0.1,
            out_geom: 0.95,
        });
        assert!(
            hostile.equal() > 0.0,
            "the floor must keep `equal` positive"
        );
        assert_eq!(hostile.equal(), 0.01);

        let total = hostile.equal()
            + hostile.in_up()
            + hostile.in_down()
            + hostile.out_up()
            + hostile.out_down();
        assert!(total > 1.0, "expected the floor to push the total past one");
    }

    /// With well-behaved masses the floor does not bind, and `equal` is exactly the
    /// remainder.
    #[test]
    fn equal_is_the_remainder_when_the_floor_does_not_bind() {
        let model = StutterModel::hipstr_shipped();
        assert!((model.equal() - (1.0 - 0.05 - 0.05 - 0.01 - 0.01)).abs() < 1e-15);
        assert_eq!(model.probability(0, period(4)), model.equal());
    }

    /// The geometrics are held strictly inside `(0, 1)`, so neither a certainty nor an
    /// impossibility is expressible (spec §5.2, trap 2). Asserted against the **contractual
    /// values** 0.01 and 0.99, not against the constants the implementation uses — otherwise
    /// editing a constant would move both sides of the assertion.
    #[test]
    fn the_geometrics_are_held_strictly_inside_zero_and_one() {
        let extreme = StutterModel::new(StutterRates {
            in_up: 0.05,
            in_down: 0.05,
            in_geom: 0.0,
            out_up: 0.01,
            out_down: 0.01,
            out_geom: 1.0,
        });
        assert_eq!(extreme.in_geom(), 0.01);
        assert_eq!(extreme.out_geom(), 0.99);
    }

    /// **Every ill-formed rate must still yield a probability**, because the debug assertion
    /// that rejects them is compiled out of the release build this project runs. Each of
    /// these produced a non-probability before the sanitizing: a mass of 2.0 gave 1.9, a
    /// negative mass gave a negative score, and `NaN` poisoned every score *while the
    /// constructor reported a healthy floor* — `f64::max` absorbs `NaN` and `f64::clamp`
    /// passes it through, so only the results were wrong.
    ///
    /// Goes through `sanitized`, which is `new` minus the debug assertion — the release path
    /// is unreachable through `new` in this profile.
    #[test]
    fn ill_formed_rates_still_yield_probabilities() {
        let well_formed = StutterRates {
            in_up: 0.05,
            in_down: 0.05,
            in_geom: 0.95,
            out_up: 0.01,
            out_down: 0.01,
            out_geom: 0.95,
        };
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.5, 2.0, 1e300] {
            for slot in 0..6 {
                let mut rates = well_formed;
                match slot {
                    0 => rates.in_up = bad,
                    1 => rates.in_down = bad,
                    2 => rates.in_geom = bad,
                    3 => rates.out_up = bad,
                    4 => rates.out_down = bad,
                    _ => rates.out_geom = bad,
                }
                let model = StutterModel::sanitized(rates);
                assert!(model.equal().is_finite(), "rate {bad} in slot {slot}");
                for period in (1..=6u8).map(period) {
                    for bp_diff in [-9i64, -4, -1, 0, 1, 4, 9] {
                        let probability = model.probability(bp_diff, period);
                        assert!(
                            probability.is_finite() && (0.0..=1.0).contains(&probability),
                            "rate {bad} in slot {slot} gave {probability} at Δ={bp_diff}"
                        );
                    }
                }
            }
        }
    }

    /// **Past the cutoff the answer is zero, not a small number** — an implausibly large
    /// change must not be explained away as stutter. Pinned on **both** sides of the
    /// boundary, in both regimes and both directions, so a `>` / `>=` slip cannot survive.
    #[test]
    fn slips_past_the_cutoff_score_zero() {
        let model = all_distinct();
        let period = period(4);
        let unit = i64::from(period.get());
        let last = i64::from(MAX_SLIP) * unit;

        // In frame: exactly at the cutoff is still scored; one unit past is zero.
        assert!(model.probability(last, period) > 0.0);
        assert!(model.probability(-last, period) > 0.0);
        assert_eq!(model.probability(last + unit, period), 0.0);
        assert_eq!(model.probability(-last - unit, period), 0.0);

        // Out of frame: the boundary is on the *re-indexed* size. Δ = 13 → e = 10 (scored),
        // Δ = 14 → e = 11 (zero).
        assert_eq!(13 - 13 / unit, 10);
        assert_eq!(14 - 14 / unit, 11);
        assert!(model.probability(13, period) > 0.0);
        assert!(model.probability(-13, period) > 0.0);
        assert_eq!(model.probability(14, period), 0.0);
        assert_eq!(model.probability(-14, period), 0.0);
    }

    /// **The one constant applied at two scales**, reproduced from production and pinned so
    /// the asymmetry is visible rather than surprising: in frame the cutoff counts *units*,
    /// out of frame it counts *re-indexed base pairs*, so out-of-frame changes are cut off
    /// about `period − 1` times sooner in real terms.
    #[test]
    fn the_cutoff_counts_units_in_frame_and_base_pairs_out_of_frame() {
        let model = all_distinct();

        // Period 4: 40 bp accepted in frame, but 14 bp already rejected out of frame.
        assert!(model.probability(40, period(4)) > 0.0);
        assert_eq!(model.probability(14, period(4)), 0.0);

        // Period 2: the effect very nearly vanishes — both scales admit about 20 bp, which
        // is why the rationale says `period − 1` and not `period`.
        assert!(model.probability(20, period(2)) > 0.0);
        assert_eq!(19 - 19 / 2, 10);
        assert!(model.probability(19, period(2)) > 0.0);
        assert_eq!(model.probability(21, period(2)), 0.0);
    }

    /// Period 1 is the case the comparison most needs, not the one to skip: **every**
    /// change is in frame there, so the out-of-frame branch is unreachable and the
    /// arithmetic regime split collapses. What does *not* collapse is direction, size decay
    /// and placement multiplicity (spec §4.2).
    #[test]
    fn every_change_is_in_frame_at_period_one() {
        let model = all_distinct();
        for bp_diff in [-3i64, -1, 1, 2, 5] {
            let mass = if bp_diff < 0 { 0.07 } else { 0.03 };
            let steps = bp_diff.unsigned_abs();
            let expected = mass * 0.95 * 0.05f64.powi((steps - 1) as i32);
            assert!((model.probability(bp_diff, period(1)) - expected).abs() < 1e-15);
        }
        // Direction asymmetry survives at period 1 — the one-penalty model cannot express
        // it at all, which is what Milestone D's comparison is about.
        assert!(model.probability(-1, period(1)) > model.probability(1, period(1)));
    }

    /// Every probability the model returns is a real probability — finite and in `[0, 1]` —
    /// across both regimes, both directions, several periods, and past the cutoff.
    #[test]
    fn every_probability_is_finite_and_within_zero_and_one() {
        for model in [
            StutterModel::hipstr_shipped(),
            StutterModel::hipstr_em_start(),
            all_distinct(),
        ] {
            for period in (1..=6u8).map(period) {
                for bp_diff in -60i64..=60 {
                    let probability = model.probability(bp_diff, period);
                    assert!(
                        probability.is_finite(),
                        "non-finite at Δ={bp_diff}, period {period:?}"
                    );
                    assert!(
                        (0.0..=1.0).contains(&probability),
                        "probability {probability} out of range at Δ={bp_diff}, period {period:?}"
                    );
                }
            }
        }
    }

    /// The two HipSTR parameter rows are matched sets, and mixing them yields a pairing that
    /// exists nowhere (spec §5.2 records an earlier draft of the spec doing exactly that).
    /// The named constructors are what stop that happening by hand; this pins their contents.
    #[test]
    fn the_two_hipstr_parameter_sets_are_kept_as_matched_rows() {
        let shipped = StutterModel::hipstr_shipped();
        assert_eq!(shipped.in_geom(), 0.95);
        assert_eq!(shipped.out_geom(), 0.95);
        assert_eq!(shipped.in_up(), 0.05);
        assert_eq!(shipped.in_down(), 0.05);
        assert_eq!(shipped.out_up(), 0.01);

        let em_start = StutterModel::hipstr_em_start();
        assert_eq!(em_start.in_geom(), 0.9);
        assert_eq!(em_start.out_geom(), 0.8);
        assert_eq!(em_start.in_up(), 0.1);
        assert_eq!(em_start.out_up(), 0.01);

        // The rows differ in their geometrics — which is the whole reason to keep them
        // apart, and also why HipSTR treats the two geometrics as independent.
        assert_ne!(em_start.in_geom(), em_start.out_geom());
    }

    /// **The copied cutoff must not drift from production's.** The doc on [`MAX_SLIP`]
    /// asserts they stay equal while the two models are meant to agree; this makes that true
    /// rather than aspirational. A **test-only** reference, so shipping ng code still depends
    /// on nothing in production.
    #[test]
    fn the_copied_cutoff_still_equals_productions() {
        assert_eq!(
            MAX_SLIP as usize,
            crate::ssr::cohort::param_estimation::MAX_SLIP,
            "ng's MAX_SLIP has drifted from production's"
        );
    }
}
