//! Stage-2 parameter & accumulator data model — the frozen [`ParamSet`] the
//! genotyping EM consumes, the sufficient-statistic accumulators the pre-pass
//! fills, and the order-independent [`FixedPointAccum`] both halves reduce floats
//! through.
//!
//! These are the A1 *nouns* (arch `ssr_call_parameters.md` §0/§3, the
//! implementation plan §3.1–§3.2). The estimation logic that fills them lives in
//! the pre-pass (`prepass.rs`, `sample_groups.rs`) and the stutter kernel
//! (`stutter.rs`); the only behaviour defined here is [`FixedPointAccum`]'s
//! accumulation, because its associativity is the determinism guarantee the rest
//! of the stage is built on.

use std::collections::HashMap;

/// Largest stutter slip `Δ` (in **repeat units**, either direction) the kernel and
/// the [`SlipProfile`] accumulator model.
///
/// This is a **compile-time array bound** (see [`SlipProfile`]), so it must hold a
/// concrete value now; `10` is a provisional choice recalibrated on the simulator
/// in F2.
pub(crate) const MAX_SLIP: usize = 10;

/// A data-driven cluster of samples sharing chemistry + provenance (arch §0).
///
/// Unobserved (no reliable user/SRA labels) — inferred in the pre-pass from each
/// sample's `ε` + stutter level. Used as a dense index into the per-group parameter
/// vectors of [`ParamSet`] and as half of the per-`(group, period)` shape key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct SampleGroupId(pub(crate) u16);

/// Per-base within-tract **substitution** rate `ε` (arch §0; C1: `align` admits no
/// sub-motif indel inside the tract, so `ε` moves composition, not length).
///
/// Frozen per sample group by the pre-pass — it lives inside `align`, making that a
/// pure, iteration-invariant function.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PerBaseError(pub(crate) f64);

/// Stutter **shape** — *where* a slip lands (arch §0).
///
/// The same geometric form is used at every granularity: the per-period parent
/// `θ_period`, the per-`(group, period)` cell, and the per-locus refinement
/// `θ_locus`. `down_rate` usually exceeds `up_rate` (contraction bias); `decay` is
/// the geometric falloff `ρ` per extra unit of `|Δ|`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct StutterShape {
    /// Mass placed on expansions (`+Δ`).
    pub(crate) up_rate: f64,
    /// Mass placed on contractions (`−Δ`); usually `> up_rate`.
    pub(crate) down_rate: f64,
    /// Geometric decay `ρ` applied per additional unit of `|Δ|`.
    pub(crate) decay: f64,
}

/// Stutter **level** — *how often* a slip happens, per sample group, as a line in
/// repeat length: `level(length) = baseline + slope · length` (arch §0).
///
/// Pre-pass-seeded then refined in the genotyping outer loop (it is an `S_θ`
/// re-weight outside `align`, so refining it rebuilds nothing).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct StutterLevel {
    /// Intercept — the stutter rate at the short-length end of the modeled range.
    pub(crate) baseline: f64,
    /// Per-unit increase in stutter rate (the length response).
    pub(crate) slope: f64,
}

/// `G₀` geometric pseudocount decay on the allele-frequency prior `π`, per loci
/// group (= period in v1; arch §0/§5).
///
/// `pseudocount(Δ) = max(p^|Δ|, FLOOR)` where `Δ` is the candidate's unit offset
/// from the per-locus modal allele.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct G0PseudocountDecay {
    /// The geometric base `p` of the per-unit decay.
    pub(crate) p: f64,
}

/// Controls for fitting the per-period `G₀` decay `p` from the cohort's confident
/// germline allele spread (arch `ssr_call_driver.md` §9). Pinned in F2.
///
/// The fit reflects "*given* a locus varies, how spread are its alleles," so it runs on
/// variable loci only; thin periods fall back to [`fallback_p`](Self::fallback_p), and
/// the fit is clamped no steeper than [`min_p`](Self::min_p) so a relatedness- or
/// structure-compressed cohort cannot over-tighten the prior.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct G0FitCfg {
    /// Minimum allele-copies for a period before the data-driven fit is trusted; below
    /// it the period takes the coded [`fallback_p`](Self::fallback_p).
    pub(crate) min_copies: u64,
    /// The coded fallback decay used for a thin period.
    pub(crate) fallback_p: f64,
    /// The steepest (smallest) `p` the fit may return — the over-tightening guard.
    pub(crate) min_p: f64,
}

/// Coded `G₀` decay `p` for a period with no data-driven fit — the single source of truth
/// for the fallback, shared by [`G0FitCfg::dev_default`]'s `fallback_p` and the genotyping
/// EM's own last-resort fallback (`em::period_decay`) so recalibrating one moves both
/// (review M3).
pub(crate) const DEFAULT_G0_FALLBACK_P: f64 = 0.5;

impl G0FitCfg {
    pub(crate) fn dev_default() -> Self {
        Self {
            min_copies: 30,
            fallback_p: DEFAULT_G0_FALLBACK_P,
            min_p: 0.1,
        }
    }
}

/// The **frozen** output of the pre-pass that the genotyping EM consumes (the
/// plan §1 interface).
///
/// `level_seed` is the SEED only — refined in the outer loop (E1), not frozen.
/// Everything else is fixed for the whole genotyping pass.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ParamSet {
    /// `ε`, frozen per sample group; indexed by [`SampleGroupId`].
    pub(crate) error_per_sample_group: Vec<PerBaseError>,
    /// `θ_period`: the cohort-per-period shape parent, keyed by period.
    pub(crate) stutter_shape_parent: HashMap<u8, StutterShape>,
    /// `θ_(group,period)`: the per-cell shape, shrunk toward the parent.
    pub(crate) stutter_shape_by_cell: HashMap<(SampleGroupId, u8), StutterShape>,
    /// Per-sample-group stutter-level seed (refined in E1, not frozen).
    pub(crate) level_seed: Vec<StutterLevel>,
    /// `G₀` decay, per loci group (= period in v1).
    pub(crate) pseudocount_decay_per_loci_group: HashMap<u8, G0PseudocountDecay>,
    /// Sample index → its sample group.
    pub(crate) group_of_sample: Vec<SampleGroupId>,
    /// `F⁰` seed for the Phase-3 prior-side `F_i` loop (estimated there, not here).
    pub(crate) f0_seed: f64,
}

/// Per `(group, period)` — the conditional profile of a slip (arch §3).
///
/// Integer counts, so the reduce is commutative and order-independent (the
/// determinism trick). `down[k]` / `up[k]` count slips of `−(k+1)` / `+(k+1)`
/// units.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct SlipProfile {
    /// Contraction counts: `down[k]` is the count of `Δ = −(k+1)`.
    pub(crate) down: [u64; MAX_SLIP],
    /// Expansion counts: `up[k]` is the count of `Δ = +(k+1)`.
    pub(crate) up: [u64; MAX_SLIP],
}

impl SlipProfile {
    /// Add a slip of signed size `delta` (count `count`) to the profile, respecting the
    /// `MAX_SLIP` cap (a `|Δ|` beyond the cap is dropped). Shared by the pre-pass `θ`
    /// accumulators and the per-locus `θ_locus` M-step.
    pub(crate) fn add_slip(&mut self, delta: i32, count: u64) {
        let magnitude = delta.unsigned_abs() as usize;
        if (1..=MAX_SLIP).contains(&magnitude) {
            if delta > 0 {
                self.up[magnitude - 1] += count;
            } else {
                self.down[magnitude - 1] += count;
            }
        }
    }
}

/// One repeat-length bin of a sample's stutter tally: at `length` repeat units, how many
/// reads sat faithfully on the peak vs slipped off it. The level line
/// `baseline + slope·length` is fit from these. A named struct so `bin.faithful += …`
/// reads itself rather than `bin.1 += …` (review Mi5); the `Vec` + linear scan is kept
/// deliberately — the distinct-length cardinality per sample is tiny.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct LengthBin {
    pub(crate) length: u16,
    pub(crate) faithful: u64,
    pub(crate) slipped: u64,
}

/// Per sample — feeds (a) the per-group level line, (b) the clustering, and
/// (c) — once groups are fixed — the per-`(group, period)` shape key (arch §3).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct SampleStutterStats {
    /// Per-repeat-length faithful/slipped bins → the level line
    /// `baseline + slope · length`. A separated het adds a bin at each allele
    /// length.
    pub(crate) by_length: Vec<LengthBin>,
    /// Within-tract base matches off the faithful peaks → `ε`.
    pub(crate) base_match: u64,
    /// Within-tract base mismatches off the faithful peaks → `ε` (substitutions).
    pub(crate) base_mismatch: u64,
    /// Dup-free primary depth — a clustering **precision** weight (`1/√depth`)
    /// only, never a feature (arch §4; dup-free by construction).
    pub(crate) read_depth: u64,
}

/// Per loci group (= period) — the count-weighted spread of confident **germline**
/// alleles around the per-locus modal allele, over **variable loci only**, used to fit
/// the `G₀` decay `p` (arch `ssr_call_driver.md` §9).
///
/// Two integer running sums, not a histogram: they are the exact sufficient statistic
/// for the closed-form mean distance `K̄ = distance_weighted / total_copies`, and being
/// integer the reduce is commutative + associative (byte-identical across threads).
/// Storing sums rather than a binned histogram also sidesteps any fixed distance bound
/// (a long tract can put an allele far from the mode).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct AlleleSpreadAccum {
    /// `Σ C[k]` — total allele-copies counted across the period's variable loci.
    pub(crate) total_copies: u64,
    /// `Σ k·C[k]` — distance-weighted copies, `k = |length − mode|` in repeat units.
    pub(crate) distance_weighted: u64,
}

/// Fixed-point scale: floats are rounded to multiples of `2⁻⁴⁰` before summing, so
/// the reduce is exact integer addition — associative, commutative, and identical
/// across thread counts.
const SCALE: f64 = (1u64 << 40) as f64;

/// Order-independent float reduce (arch §3, verify-fix #1): scale → round → sum
/// into an `i128`.
///
/// Used for the decision floats whose value must not depend on the order loci
/// complete in — `ℓ_pen` (the burn-in plateau stop), the BIC log-likelihoods, and
/// the `F` / level responsibility reduces. Integer addition is associative and
/// commutative, so partials reduced in any order, or merged in any grouping, yield
/// the identical bit pattern from [`value`](Self::value).
///
/// Inputs to [`add`](Self::add) are expected **finite and bounded** well within
/// `i128::MAX / 2⁴⁰ ≈ 1.5e26` (per-locus log-normalizers are `reads × O(few
/// nats)`). A non-finite or out-of-range input would otherwise be absorbed
/// silently (`NaN → 0`, overflow → saturation) and defeat the whole point — that a
/// non-monotone `ℓ_pen` is a *sharp* diagnostic — so `add` debug-asserts finiteness.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct FixedPointAccum {
    acc: i128,
}

impl FixedPointAccum {
    /// A zeroed accumulator.
    pub(crate) fn new() -> Self {
        Self { acc: 0 }
    }

    /// Add one float contribution (scaled and rounded to the fixed-point grid).
    ///
    /// Inputs must be finite (see the type's magnitude contract); a non-finite
    /// value is a caller bug and is caught in debug builds rather than silently
    /// absorbed as `0`.
    pub(crate) fn add(&mut self, x: f64) {
        debug_assert!(
            x.is_finite(),
            "FixedPointAccum::add got a non-finite value: {x}"
        );
        self.acc += (x * SCALE).round() as i128;
    }

    /// Fold another accumulator's running total in (associative + commutative).
    pub(crate) fn merge(&mut self, other: &Self) {
        self.acc += other.acc;
    }

    /// The accumulated value as an `f64` (a pure function of the integer total, so
    /// bit-identical for equal totals).
    pub(crate) fn value(&self) -> f64 {
        self.acc as f64 / SCALE
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn param_set_round_trips_its_fields() {
        let mut parent = HashMap::new();
        parent.insert(
            2u8,
            StutterShape {
                up_rate: 0.01,
                down_rate: 0.05,
                decay: 0.3,
            },
        );
        let mut by_cell = HashMap::new();
        by_cell.insert(
            (SampleGroupId(0), 2u8),
            StutterShape {
                up_rate: 0.02,
                down_rate: 0.06,
                decay: 0.3,
            },
        );
        let mut g0 = HashMap::new();
        g0.insert(2u8, G0PseudocountDecay { p: 0.25 });

        let params = ParamSet {
            error_per_sample_group: vec![PerBaseError(0.001), PerBaseError(0.002)],
            stutter_shape_parent: parent,
            stutter_shape_by_cell: by_cell,
            level_seed: vec![StutterLevel {
                baseline: 0.01,
                slope: 0.002,
            }],
            pseudocount_decay_per_loci_group: g0,
            group_of_sample: vec![SampleGroupId(0), SampleGroupId(1), SampleGroupId(0)],
            f0_seed: 0.05,
        };

        assert_eq!(params.error_per_sample_group[1], PerBaseError(0.002));
        assert_eq!(params.stutter_shape_parent[&2].down_rate, 0.05);
        assert_eq!(
            params.stutter_shape_by_cell[&(SampleGroupId(0), 2)].up_rate,
            0.02
        );
        assert_eq!(params.group_of_sample[2], SampleGroupId(0));
        assert_eq!(params.pseudocount_decay_per_loci_group[&2].p, 0.25);
        assert_eq!(params.f0_seed, 0.05);
    }

    #[test]
    fn slip_profile_default_is_zeroed() {
        let p = SlipProfile::default();
        assert!(p.down.iter().all(|&c| c == 0));
        assert!(p.up.iter().all(|&c| c == 0));
        assert_eq!(p.down.len(), MAX_SLIP);
    }

    #[test]
    fn fixed_point_accum_value_is_independent_of_add_order() {
        let xs = [0.1, -0.4, 12.75, 0.000_3, -7.2, 3.333_333, 1e6, -1e6 + 0.5];
        let mut forward = FixedPointAccum::new();
        for &x in &xs {
            forward.add(x);
        }
        let mut backward = FixedPointAccum::new();
        for &x in xs.iter().rev() {
            backward.add(x);
        }
        // Bit-identical, not merely approximately equal — this is the determinism
        // guarantee the whole stage relies on.
        assert_eq!(forward.value().to_bits(), backward.value().to_bits());
    }

    #[test]
    #[should_panic(expected = "non-finite")]
    fn fixed_point_accum_rejects_non_finite_in_debug() {
        // A NaN would otherwise be absorbed as 0 and hide a broken upstream reduce.
        let mut acc = FixedPointAccum::new();
        acc.add(f64::NAN);
    }

    #[test]
    fn fixed_point_accum_merge_matches_sequential_and_is_grouping_independent() {
        let xs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];

        let mut sequential = FixedPointAccum::new();
        for &x in &xs {
            sequential.add(x);
        }

        // Split into two partials and merge — must equal the sequential total.
        let mut left = FixedPointAccum::new();
        for &x in &xs[..2] {
            left.add(x);
        }
        let mut right = FixedPointAccum::new();
        for &x in &xs[2..] {
            right.add(x);
        }
        left.merge(&right);
        assert_eq!(left.value().to_bits(), sequential.value().to_bits());

        // A different grouping of the same partials reduces identically.
        let mut a = FixedPointAccum::new();
        for &x in &xs[..4] {
            a.add(x);
        }
        let mut b = FixedPointAccum::new();
        for &x in &xs[4..] {
            b.add(x);
        }
        b.merge(&a);
        assert_eq!(b.value().to_bits(), sequential.value().to_bits());
    }
}
