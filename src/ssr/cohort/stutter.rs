//! The stutter kernel and slip reachability (arch `ssr_call_genotyping.md` §3/§4,
//! spec §5.2/§6/§7, implementation plan §3.4, B2).
//!
//! Three pieces: [`s_theta`] is the kernel `S_θ(Δ) = level × shape(Δ)` — the
//! *scoring* counterpart of the simulator's *generative* slip model (the plan §B
//! contract: this and `sim.rs::slip_length` are one model from two ends);
//! [`reach_variants`] enumerates the placement-distinct realizations of
//! `candidate ⊕ Δ` (one for a pure tract, ≈ #interruptions + 1 for an impure one);
//! and [`refine_theta_locus`] is the per-locus shape M-step, shrunk toward the
//! `(group, period)` prior (which is itself shrunk toward `θ_period`).

use crate::ssr::cohort::param_estimation::{MAX_SLIP, SlipProfile, StutterShape};
use crate::ssr::types::Motif;

/// One placement-distinct realization of `candidate ⊕ Δ` (arch §4, verify-fix #3).
///
/// A pure tract has a single placement of a `Δ`-unit slip; an impure tract (runs
/// separated by interruptions) has one variant per run the slip could land in. The
/// likelihood sums `align` over these with a uniform position prior.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PlacementVariant {
    /// The realized tract sequence for this placement of the slip.
    pub(crate) seq: Box<[u8]>,
}

/// The stutter kernel `S_θ(Δ)` — the probability a read of a parent allele shows a
/// slip of `Δ` whole motif units, at the supplied per-read `level` (≈ `P(Δ ≠ 0)`).
///
/// Mirrors the simulator's forward draw (`sim.rs`): `Δ = 0` is faithful with mass
/// `1 − level`; a slip splits by direction (`up_rate : down_rate`) and its magnitude
/// `m = |Δ|` is geometric in `decay`, truncated at [`MAX_SLIP`] with the tail
/// absorbed into the last bin — so `Σ_{Δ=-MAX..MAX} S_θ(Δ) = 1` exactly.
pub(crate) fn s_theta(delta: i32, shape: &StutterShape, level: f64) -> f64 {
    if delta == 0 {
        return 1.0 - level;
    }
    let magnitude = delta.unsigned_abs() as usize;
    if magnitude > MAX_SLIP {
        return 0.0;
    }
    let direction_mass = shape.up_rate + shape.down_rate;
    let direction_fraction = if direction_mass > 0.0 {
        if delta > 0 {
            shape.up_rate / direction_mass
        } else {
            shape.down_rate / direction_mass
        }
    } else {
        0.5
    };
    level * direction_fraction * geometric_magnitude(magnitude, shape.decay)
}

/// `P(magnitude = m)` for `m ∈ 1..=MAX_SLIP`: geometric with continuation
/// probability `decay`, the final bin absorbing the truncated tail (matching the
/// simulator's `while m < MAX_SLIP && chance(decay)` loop, so the bins sum to 1).
fn geometric_magnitude(magnitude: usize, decay: f64) -> f64 {
    debug_assert!((1..=MAX_SLIP).contains(&magnitude));
    if magnitude < MAX_SLIP {
        (1.0 - decay) * decay.powi(magnitude as i32 - 1)
    } else {
        decay.powi(MAX_SLIP as i32 - 1)
    }
}

/// A candidate tract decomposed into repeat runs and the interruptions between them.
enum Segment {
    /// A maximal run of `units` whole motif copies.
    Run(usize),
    /// A stretch of bases that is not part of a motif copy (an interruption).
    Interruption(Vec<u8>),
}

/// Split `cand` into alternating motif runs and interruption segments (greedy
/// left-to-right motif matching). A pure tract yields a single [`Segment::Run`].
fn segment(cand: &[u8], motif: &[u8]) -> Vec<Segment> {
    let period = motif.len();
    let mut segments = Vec::new();
    let mut run_units = 0usize;
    let mut interruption: Vec<u8> = Vec::new();
    let mut i = 0;
    while i < cand.len() {
        if i + period <= cand.len() && &cand[i..i + period] == motif {
            if !interruption.is_empty() {
                segments.push(Segment::Interruption(std::mem::take(&mut interruption)));
            }
            run_units += 1;
            i += period;
        } else {
            if run_units > 0 {
                segments.push(Segment::Run(run_units));
                run_units = 0;
            }
            interruption.push(cand[i]);
            i += 1;
        }
    }
    if run_units > 0 {
        segments.push(Segment::Run(run_units));
    }
    if !interruption.is_empty() {
        segments.push(Segment::Interruption(interruption));
    }
    segments
}

/// Render `segments` to bytes, overriding the unit count of the run at
/// `target_run_index` to `new_units`.
fn render(
    segments: &[Segment],
    motif: &[u8],
    target_run_index: usize,
    new_units: usize,
) -> Box<[u8]> {
    let mut out = Vec::new();
    let mut run_seen = 0usize;
    for seg in segments {
        match seg {
            Segment::Run(units) => {
                let u = if run_seen == target_run_index {
                    new_units
                } else {
                    *units
                };
                for _ in 0..u {
                    out.extend_from_slice(motif);
                }
                run_seen += 1;
            }
            Segment::Interruption(bytes) => out.extend_from_slice(bytes),
        }
    }
    out.into_boxed_slice()
}

/// Enumerate the placement-distinct variants of `cand ⊕ delta` into `out` (cleared
/// first). A pure candidate yields one variant; an impure one yields up to one per
/// run, deduplicated. A run that cannot absorb the slip (would go below zero units)
/// contributes no variant.
pub(crate) fn reach_variants(
    cand: &[u8],
    motif: &Motif,
    delta: i32,
    out: &mut Vec<PlacementVariant>,
) {
    out.clear();
    if delta == 0 {
        out.push(PlacementVariant {
            seq: cand.to_vec().into_boxed_slice(),
        });
        return;
    }
    let motif_bytes = motif.as_bytes();
    let segments = segment(cand, motif_bytes);
    let run_count = segments
        .iter()
        .filter(|s| matches!(s, Segment::Run(_)))
        .count();
    for target in 0..run_count {
        // The unit count of the `target`-th run.
        let units = segments
            .iter()
            .filter_map(|s| match s {
                Segment::Run(u) => Some(*u),
                Segment::Interruption(_) => None,
            })
            .nth(target)
            .unwrap();
        let new_units = units as i32 + delta;
        if new_units < 0 {
            continue;
        }
        let seq = render(&segments, motif_bytes, target, new_units as usize);
        let variant = PlacementVariant { seq };
        if !out.contains(&variant) {
            out.push(variant);
        }
    }
}

/// The per-locus shape M-step: estimate `θ_locus` from a slip profile and shrink it
/// toward `prior` (the sample's `(group, period)` shape) with pseudo-count weight
/// `strength`. With no observed slips it collapses to `prior`.
pub(crate) fn refine_theta_locus(
    profile: &SlipProfile,
    prior: &StutterShape,
    strength: f64,
) -> StutterShape {
    let up_total: u64 = profile.up.iter().sum();
    let down_total: u64 = profile.down.iter().sum();
    let total = up_total + down_total;
    if total == 0 {
        return *prior;
    }

    // Direction fractions, shrunk toward the prior's direction split.
    let prior_mass = prior.up_rate + prior.down_rate;
    let prior_up_fraction = if prior_mass > 0.0 {
        prior.up_rate / prior_mass
    } else {
        0.5
    };
    let n = total as f64;
    let up_fraction = (up_total as f64 + strength * prior_up_fraction) / (n + strength);
    let down_fraction = 1.0 - up_fraction;

    // Geometric decay MLE from the magnitude distribution: for P(m) ∝ decay^(m-1),
    // the mean magnitude is 1/(1−decay), so decay = (mean − 1) / mean.
    let mut weighted_magnitude = 0.0;
    for k in 0..MAX_SLIP {
        let magnitude = (k + 1) as f64;
        weighted_magnitude += magnitude * (profile.up[k] + profile.down[k]) as f64;
    }
    let mean_magnitude = weighted_magnitude / n;
    let decay_estimate = ((mean_magnitude - 1.0) / mean_magnitude).clamp(0.0, 0.999);
    let decay = (n * decay_estimate + strength * prior.decay) / (n + strength);

    StutterShape {
        up_rate: up_fraction,
        down_rate: down_fraction,
        decay,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn motif(bytes: &[u8]) -> Motif {
        Motif::new(bytes).unwrap()
    }

    fn shape() -> StutterShape {
        StutterShape {
            up_rate: 1.0,
            down_rate: 2.0, // contraction bias
            decay: 0.2,
        }
    }

    #[test]
    fn placement_variant_holds_its_sequence() {
        let v = PlacementVariant {
            seq: Box::from(b"ATATATAT".as_slice()),
        };
        assert_eq!(&*v.seq, b"ATATATAT");
    }

    #[test]
    fn s_theta_normalizes_over_all_slips() {
        let level = 0.3;
        let sum: f64 = (-(MAX_SLIP as i32)..=(MAX_SLIP as i32))
            .map(|d| s_theta(d, &shape(), level))
            .sum();
        assert!((sum - 1.0).abs() < 1e-12, "kernel must sum to 1, got {sum}");
    }

    #[test]
    fn s_theta_faithful_mass_is_one_minus_level() {
        assert!((s_theta(0, &shape(), 0.3) - 0.7).abs() < 1e-12);
    }

    #[test]
    fn s_theta_respects_contraction_bias_and_decay() {
        let s = shape();
        // down (contraction) carries more mass than up at the same magnitude.
        assert!(s_theta(-1, &s, 0.3) > s_theta(1, &s, 0.3));
        // larger magnitude decays.
        assert!(s_theta(-1, &s, 0.3) > s_theta(-2, &s, 0.3));
        // beyond MAX_SLIP is unreachable.
        assert_eq!(s_theta(MAX_SLIP as i32 + 1, &s, 0.3), 0.0);
    }

    #[test]
    fn reach_variants_pure_allele_is_a_single_tiling() {
        let mut out = Vec::new();
        // (CA)x5 expanded by +2 → (CA)x7, one variant.
        reach_variants(b"CACACACACA", &motif(b"CA"), 2, &mut out);
        assert_eq!(out.len(), 1);
        assert_eq!(&*out[0].seq, b"CACACACACACACA"); // CA x 7
    }

    #[test]
    fn reach_variants_delta_zero_returns_the_candidate() {
        let mut out = Vec::new();
        reach_variants(b"CACACA", &motif(b"CA"), 0, &mut out);
        assert_eq!(out.len(), 1);
        assert_eq!(&*out[0].seq, b"CACACA");
    }

    #[test]
    fn reach_variants_one_interruption_gives_two_placements() {
        // (CA)x3 TT (CA)x2 — one interruption ⇒ two runs ⇒ two +1 placements.
        let mut out = Vec::new();
        reach_variants(b"CACACATTCACA", &motif(b"CA"), 1, &mut out);
        assert_eq!(out.len(), 2);
        // +1 in the first run, or +1 in the second run.
        assert!(out.iter().any(|v| &*v.seq == b"CACACACATTCACA"));
        assert!(out.iter().any(|v| &*v.seq == b"CACACATTCACACA"));
    }

    #[test]
    fn reach_variants_pure_contraction_shortens_the_tiling() {
        // (CA)x5 contracted by −2 → (CA)x3, one variant.
        let mut out = Vec::new();
        reach_variants(b"CACACACACA", &motif(b"CA"), -2, &mut out);
        assert_eq!(out.len(), 1);
        assert_eq!(&*out[0].seq, b"CACACA"); // CA x 3
    }

    #[test]
    fn reach_variants_skips_runs_that_cannot_contract() {
        // A 2-unit run cannot lose 3 units; with a single run that means no variant.
        let mut out = Vec::new();
        reach_variants(b"CACA", &motif(b"CA"), -3, &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn refine_theta_locus_collapses_to_prior_without_data() {
        let prior = shape();
        let refined = refine_theta_locus(&SlipProfile::default(), &prior, 5.0);
        assert_eq!(refined, prior);
    }

    #[test]
    fn refine_theta_locus_moves_toward_observed_direction_and_decay() {
        // Profile dominated by single-unit contractions: down direction, low decay.
        let mut profile = SlipProfile::default();
        profile.down[0] = 90; // Δ = −1
        profile.down[1] = 10; // Δ = −2
        profile.up[0] = 5; // Δ = +1
        // A weak prior so the data dominates.
        let prior = StutterShape {
            up_rate: 0.5,
            down_rate: 0.5,
            decay: 0.5,
        };
        let refined = refine_theta_locus(&profile, &prior, 1.0);
        assert!(
            refined.down_rate > refined.up_rate,
            "observed contraction bias should pull down_rate up"
        );
        // mean magnitude ≈ (90·1 + 10·2 + 5·1)/105 ≈ 1.10 → decay ≈ 0.09, well below
        // the 0.5 prior.
        assert!(
            refined.decay < 0.3,
            "decay should track the short-slip data"
        );
    }
}
