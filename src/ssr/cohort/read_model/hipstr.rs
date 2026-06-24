//! Model A — the HipSTR explicit in-frame + out-of-frame stutter model (bake-off plan
//! §2 Model A; mirrors `HipSTR/src/stutter_model.cpp::log_stutter_pmf`).
//!
//! **This is the production read-likelihood model** — the bake-off
//! (`ssr_stutter_scoring_model_bakeoff_2026-06-24.md`) picked it: tied top concordance
//! with Model B across clean/out-of-frame/messy data, but far better calibrated and
//! ~100× cheaper per `Qᵣ` (one `bp_diff` term + a closed-form substitution, vs B's
//! 21-way Δ-sum-with-DP).
//!
//! `Qᵣ = P_stutter(Δ_bp | period) · Σ_v (1/|v|) · subst(obs | resize(cand, Δ_bp))`,
//! where the read's length change `Δ_bp = |obs| − |cand|` is the *single* length the
//! substitution-only align must explain (a stutter resizes the candidate to the read's
//! length; the only residual is base composition). The stutter PMF splits two regimes
//! with **distinct** geometrics (HipSTR's defining structure):
//!
//! - **in-frame** (`Δ_bp` a multiple of `period`): direction mass `in_up`/`in_down` ×
//!   geometric in *units* (`in_geom`);
//! - **out-of-frame** (`Δ_bp` not a unit multiple, rarer): direction mass
//!   `out_up`/`out_down` × geometric in *bp* (`out_geom`), indexed by the HipSTR
//!   effective-difference `eff = Δ_bp − Δ_bp/period`.
//!
//! `log_equal = 1 − in_up − in_down − out_up − out_down` is the faithful mass. The
//! in-frame placements reuse B's [`reach_variants`]; out-of-frame resizing appends or
//! trims partial-motif bases at the tract end (a simpler placement than HipSTR's
//! `StutterAlignerClass`, sufficient at our tract sizes).
//!
//! **Step-4 parameterization.** Like Model C, the six HipSTR parameters are derived
//! from the same [`ReadScoringContext`] every model reads, so A is scored on identical
//! truth/pre-pass parameters (plan §9): the in-frame mass + geometric from `level` and
//! the shape's direction split + `decay`; `ε` drives `subst`. The **out-of-frame** mass
//! has no counterpart in the in-frame context, so it is a fixed small fraction
//! ([`OUT_FRAME_REL`]) of the in-frame slip mass — a declared estimator for the
//! bake-off. A real out-of-frame estimator (out-of-frame reads binned separately in the
//! pre-pass) is the productionization follow-up (Step 5).

use crate::ssr::cohort::pair_hmm::{HmmScratch, align_subst};
use crate::ssr::cohort::param_estimation::MAX_SLIP;
use crate::ssr::cohort::read_model::{ReadLikelihoodModel, ReadScoringContext};
use crate::ssr::cohort::stutter::{PlacementVariant, reach_variants};

/// Out-of-frame slip mass as a fraction of the in-frame slip mass (the Step-4 declared
/// estimator). Out-of-frame stutter is rare relative to whole-unit slips, so it carries
/// a small share; pinned to a real per-period estimate in Step 5.
const OUT_FRAME_REL: f64 = 0.05;
/// Bounds keeping each geometric success probability strictly inside `(0, 1)`.
const GEOM_MIN: f64 = 0.01;
const GEOM_MAX: f64 = 0.99;

/// Model A: HipSTR's explicit in-frame + out-of-frame geometric stutter. Stateless —
/// the six HipSTR parameters are derived per call from the [`ReadScoringContext`] (see
/// the module docs), so it shares B's parameter interface.
#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct HipstrModel;

/// The six HipSTR stutter parameters for one read-vs-candidate evaluation, plus the
/// faithful mass — the linear-probability form of `StutterModel`'s log fields.
#[derive(Debug, Clone, Copy)]
struct HipstrParams {
    /// Faithful mass `P(Δ = 0) = 1 − in_up − in_down − out_up − out_down`.
    equal: f64,
    /// In-frame expansion / contraction direction mass.
    in_up: f64,
    in_down: f64,
    /// In-frame geometric **success** probability (per unit; `1 − decay`).
    in_geom: f64,
    /// Out-of-frame expansion / contraction direction mass.
    out_up: f64,
    out_down: f64,
    /// Out-of-frame geometric **success** probability (per bp).
    out_geom: f64,
}

/// Reusable buffers: the in-frame placement-variant set, the align DP, and the
/// out-of-frame resized-candidate scratch — so scoring a locus's reads allocates once.
#[derive(Debug, Default)]
pub(crate) struct HipstrScratch {
    variants: Vec<PlacementVariant>,
    hmm: HmmScratch,
    resized: Vec<u8>,
}

impl ReadLikelihoodModel for HipstrModel {
    type Scratch = HipstrScratch;

    fn q_r(
        &self,
        obs: &[u8],
        cand: &[u8],
        ctx: &ReadScoringContext,
        scratch: &mut Self::Scratch,
    ) -> f64 {
        let params = hipstr_params(ctx);
        let period = ctx.motif.period() as i32;
        let bp_diff = obs.len() as i32 - cand.len() as i32;

        let p_stutter = stutter_pmf(bp_diff, period, &params);
        if p_stutter <= 0.0 {
            return 0.0;
        }

        let HipstrScratch {
            variants,
            hmm,
            resized,
        } = scratch;

        // Resize the candidate to the read's length, then score the residual composition
        // with the substitution-only align (equal length ⇒ no interior gaps).
        if bp_diff % period == 0 {
            let units = bp_diff / period; // in-frame: whole-unit placement(s)
            reach_variants(cand, ctx.motif, units, variants);
            if variants.is_empty() {
                return 0.0; // this Δ cannot be reached from `cand` (over-contraction)
            }
            let pr_v = 1.0 / variants.len() as f64;
            let align: f64 = variants
                .iter()
                .map(|v| pr_v * align_subst(obs, &v.seq, ctx.eps, hmm))
                .sum();
            p_stutter * align
        } else {
            resize_out_of_frame(cand, ctx.motif.as_bytes(), bp_diff, resized);
            p_stutter * align_subst(obs, resized, ctx.eps, hmm)
        }
    }
}

/// Derive the six HipSTR parameters from the per-call context (the Step-4 estimator):
/// in-frame mass + geometric from `level`, the shape's direction split, and `decay`;
/// out-of-frame mass a small fraction ([`OUT_FRAME_REL`]) of it, sharing the geometric.
fn hipstr_params(ctx: &ReadScoringContext) -> HipstrParams {
    let direction_mass = ctx.shape.up_rate + ctx.shape.down_rate;
    let up_fraction = if direction_mass > 0.0 {
        ctx.shape.up_rate / direction_mass
    } else {
        0.5
    };
    let down_fraction = 1.0 - up_fraction;

    let in_total = ctx.level;
    let out_total = ctx.level * OUT_FRAME_REL;
    // B's `decay` is the geometric *continuation* probability (mean magnitude
    // 1/(1−decay)); HipSTR's geom is the *success* probability (mean 1/geom), so the
    // matching conversion is geom = 1 − decay.
    let geom = (1.0 - ctx.shape.decay).clamp(GEOM_MIN, GEOM_MAX);

    let in_up = in_total * up_fraction;
    let in_down = in_total * down_fraction;
    let out_up = out_total * up_fraction;
    let out_down = out_total * down_fraction;
    let equal = (1.0 - in_up - in_down - out_up - out_down).max(GEOM_MIN);

    HipstrParams {
        equal,
        in_up,
        in_down,
        in_geom: geom,
        out_up,
        out_down,
        out_geom: geom,
    }
}

/// The stutter PMF `P(Δ_bp)` — the linear form of HipSTR's `log_stutter_pmf`. In-frame
/// when `Δ_bp` is a unit multiple (geometric in units), out-of-frame otherwise
/// (geometric in HipSTR's effective bp difference `eff = Δ_bp − Δ_bp/period`). Slips
/// beyond [`MAX_SLIP`] are dropped to `0` (the read falls to the EM's outlier floor).
fn stutter_pmf(bp_diff: i32, period: i32, params: &HipstrParams) -> f64 {
    if bp_diff % period != 0 {
        // Out of frame. `eff` reuses HipSTR's truncated-division re-indexing so the
        // out-of-frame geometric does not double-count the in-frame unit multiples.
        let eff = bp_diff - bp_diff / period;
        if eff.unsigned_abs() as usize > MAX_SLIP {
            return 0.0;
        }
        if eff < 0 {
            params.out_down * params.out_geom * (1.0 - params.out_geom).powi(-eff - 1)
        } else {
            params.out_up * params.out_geom * (1.0 - params.out_geom).powi(eff - 1)
        }
    } else {
        let units = bp_diff / period;
        if units == 0 {
            params.equal
        } else if units.unsigned_abs() as usize > MAX_SLIP {
            0.0
        } else if units < 0 {
            params.in_down * params.in_geom * (1.0 - params.in_geom).powi(-units - 1)
        } else {
            params.in_up * params.in_geom * (1.0 - params.in_geom).powi(units - 1)
        }
    }
}

/// Resize `cand` by `bp_diff` out-of-frame bases into `out` (cleared first): append
/// motif bases continuing the tract's phase for an expansion, or trim trailing bases for
/// a contraction. A single placement (the tract end), enough for the rare out-of-frame
/// term at our tract sizes; HipSTR's `StutterAlignerClass` places it more elaborately.
fn resize_out_of_frame(cand: &[u8], motif: &[u8], bp_diff: i32, out: &mut Vec<u8>) {
    out.clear();
    let period = motif.len();
    if bp_diff >= 0 {
        out.extend_from_slice(cand);
        // Continue the motif tiling from the tract's current phase.
        for t in 0..bp_diff as usize {
            out.push(motif[(cand.len() + t) % period]);
        }
    } else {
        let keep = cand.len().saturating_sub((-bp_diff) as usize);
        out.extend_from_slice(&cand[..keep]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::param_estimation::StutterShape;
    use crate::ssr::types::Motif;

    const EPS: f64 = 0.005;
    const LEVEL: f64 = 0.1;

    fn shape() -> StutterShape {
        StutterShape {
            up_rate: 1.0,
            down_rate: 2.0, // contraction bias
            decay: 0.2,
        }
    }

    fn ctx<'a>(motif: &'a Motif, shape: &'a StutterShape) -> ReadScoringContext<'a> {
        ReadScoringContext {
            motif,
            shape,
            level: LEVEL,
            eps: EPS,
        }
    }

    fn ca(units: usize) -> Vec<u8> {
        std::iter::repeat_n(*b"CA", units).flatten().collect()
    }

    fn score(obs: &[u8], cand: &[u8], motif: &Motif) -> f64 {
        let shape = shape();
        let mut scratch = HipstrScratch::default();
        HipstrModel.q_r(obs, cand, &ctx(motif, &shape), &mut scratch)
    }

    /// The PMF must reproduce HipSTR's `log_stutter_pmf` term-by-term (linear form).
    #[test]
    fn stutter_pmf_mirrors_the_hipstr_formula() {
        let p = hipstr_params(&ctx(&Motif::new(b"CA").unwrap(), &shape()));
        let period = 2;
        // Δ = 0 → equal.
        assert_eq!(stutter_pmf(0, period, &p), p.equal);
        // +1 unit (in-frame) → in_up · in_geom.
        assert!((stutter_pmf(2, period, &p) - p.in_up * p.in_geom).abs() < 1e-15);
        // −1 unit → in_down · in_geom.
        assert!((stutter_pmf(-2, period, &p) - p.in_down * p.in_geom).abs() < 1e-15);
        // +2 units → in_up · in_geom · (1−in_geom).
        let plus_two = p.in_up * p.in_geom * (1.0 - p.in_geom);
        assert!((stutter_pmf(4, period, &p) - plus_two).abs() < 1e-15);
        // +1 bp (out of frame) → out_up · out_geom (eff = 1).
        assert!((stutter_pmf(1, period, &p) - p.out_up * p.out_geom).abs() < 1e-15);
        // −1 bp (out of frame) → out_down · out_geom (eff = −1).
        assert!((stutter_pmf(-1, period, &p) - p.out_down * p.out_geom).abs() < 1e-15);
    }

    #[test]
    fn faithful_read_is_dominated_by_the_equal_term() {
        let motif = Motif::new(b"CA").unwrap();
        let cand = ca(8);
        let q = score(&cand, &cand, &motif);
        let p = hipstr_params(&ctx(&motif, &shape()));
        let faithful = p.equal * (1.0 - EPS).powi(16);
        assert!(q >= faithful, "q {q} ≥ equal·subst {faithful}");
        assert!(
            q < faithful * 1.05,
            "q {q} should be dominated by the equal term"
        );
    }

    #[test]
    fn contraction_bias_and_geometric_decay_order_the_in_frame_slips() {
        let motif = Motif::new(b"CA").unwrap();
        let cand = ca(8);
        let minus_one = score(&ca(7), &cand, &motif);
        let plus_one = score(&ca(9), &cand, &motif);
        let plus_two = score(&ca(10), &cand, &motif);
        // contraction bias: −1 beats +1.
        assert!(minus_one > plus_one, "−1 {minus_one} > +1 {plus_one}");
        // geometric decay: +1 beats +2.
        assert!(plus_one > plus_two, "+1 {plus_one} > +2 {plus_two}");
    }

    #[test]
    fn an_in_frame_unit_beats_an_out_of_frame_single_base() {
        // HipSTR's defining split: a whole-unit (in-frame) change is far likelier than a
        // single out-of-frame base, so the +1-unit read outscores a +1-bp out-of-frame
        // read (OUT_FRAME_REL ≪ 1).
        let motif = Motif::new(b"CA").unwrap();
        let cand = ca(6);
        let in_frame = score(&ca(7), &cand, &motif); // +2 bp, in frame
        let mut oof = cand.clone();
        oof.push(b'G'); // +1 bp, out of frame
        let out_of_frame = score(&oof, &cand, &motif);
        assert!(
            in_frame > out_of_frame,
            "in-frame unit {in_frame} should beat out-of-frame base {out_of_frame}"
        );
    }

    #[test]
    fn impure_candidate_sums_over_placements() {
        let motif = Motif::new(b"CA").unwrap();
        let cand = b"CACACATTCACA"; // (CA)3 TT (CA)2
        let q = score(cand, cand, &motif);
        assert!(q > 0.0 && q <= 1.0 + 1e-9);
        // A +1-unit slip is reachable in either run (in-frame placement sum).
        assert!(score(b"CACACACATTCACA", cand, &motif) > 0.0);
    }

    #[test]
    fn out_of_frame_resize_continues_then_trims_the_tract() {
        let motif = b"CA";
        let mut out = Vec::new();
        // +1 bp on CA×3 (len 6, phase 0 next) appends motif[6 % 2 = 0] = 'C'.
        resize_out_of_frame(b"CACACA", motif, 1, &mut out);
        assert_eq!(&out, b"CACACAC");
        // −1 bp trims one trailing base.
        resize_out_of_frame(b"CACACA", motif, -1, &mut out);
        assert_eq!(&out, b"CACAC");
    }

    #[test]
    fn long_allele_faithful_and_slipped_both_score() {
        let motif = Motif::new(b"CA").unwrap();
        let cand = ca(20);
        let faithful = score(&cand, &cand, &motif);
        let plus_two = score(&ca(22), &cand, &motif);
        assert!(faithful > 0.0 && plus_two > 0.0);
        assert!(faithful > plus_two);
    }

    #[test]
    fn scratch_reuse_gives_identical_results() {
        let motif = Motif::new(b"CA").unwrap();
        let shape = shape();
        let context = ctx(&motif, &shape);
        let mut reused = HipstrScratch::default();
        let mut fresh = HipstrScratch::default();
        // Prime `reused` with unrelated in-frame and out-of-frame alignments first.
        let _ = HipstrModel.q_r(&ca(20), &ca(4), &context, &mut reused);
        let mut oof = ca(8);
        oof.push(b'G');
        let _ = HipstrModel.q_r(&oof, &ca(8), &context, &mut reused);
        let a = HipstrModel.q_r(&ca(9), &ca(8), &context, &mut reused);
        let b = HipstrModel.q_r(&ca(9), &ca(8), &context, &mut fresh);
        assert_eq!(a.to_bits(), b.to_bits());
    }
}
