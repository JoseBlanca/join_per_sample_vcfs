//! Model B — the current production read likelihood, behind the
//! [`ReadLikelihoodModel`] trait (bake-off plan §2 Model B).
//!
//! `Qᵣ = Σ_Δ S_θ(Δ) · Σ_v (1/|v|) · align_subst(obs | cand ⊕ Δ)`: in-frame
//! (whole-unit) stutter `S_θ` marginalized *outside* a substitution-only tract align,
//! with the out-of-frame residual absorbed as flank slop and any larger mismatch
//! handled by the EM's `λ·(1/D)` outlier floor. This is a thin adapter over
//! [`read_likelihood`] so the existing likelihood unit tests keep pinning the same
//! function directly; the bake-off harness reaches the identical computation through
//! the trait.

use crate::ssr::cohort::likelihood::{LikelihoodScratch, read_likelihood};
use crate::ssr::cohort::read_model::{ReadLikelihoodModel, ReadScoringContext};

/// Model B: the current `Σ_Δ S_θ(Δ)·align_subst` likelihood. Stateless — every
/// parameter arrives through the [`ReadScoringContext`] (the EM computes the M-step
/// shape / level / ε), so the model value carries no frozen constants of its own.
#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct ClassicStutterModel;

impl ReadLikelihoodModel for ClassicStutterModel {
    type Scratch = LikelihoodScratch;

    fn q_r(
        &self,
        obs: &[u8],
        cand: &[u8],
        ctx: &ReadScoringContext,
        scratch: &mut Self::Scratch,
    ) -> f64 {
        read_likelihood(obs, cand, ctx.motif, ctx.shape, ctx.level, ctx.eps, scratch)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::param_estimation::StutterShape;
    use crate::ssr::types::Motif;

    fn shape() -> StutterShape {
        StutterShape {
            up_rate: 1.0,
            down_rate: 2.0,
            decay: 0.2,
        }
    }

    fn ca(units: usize) -> Vec<u8> {
        std::iter::repeat_n(*b"CA", units).flatten().collect()
    }

    /// The trait path must reproduce the free function bit-for-bit — Model B is a pure
    /// adapter, so the behaviour-preserving refactor is only correct if `q_r` and
    /// `read_likelihood` agree exactly (not just approximately) across obs/cand pairs.
    #[test]
    fn q_r_matches_read_likelihood_bit_for_bit() {
        let motif = Motif::new(b"CA").unwrap();
        let shape = shape();
        let model = ClassicStutterModel;
        let mut trait_scratch = LikelihoodScratch::new();
        let mut direct_scratch = LikelihoodScratch::new();

        for (obs_units, cand_units) in [(8, 8), (7, 8), (8, 3), (10, 8), (4, 9)] {
            let obs = ca(obs_units);
            let cand = ca(cand_units);
            let ctx = ReadScoringContext {
                motif: &motif,
                shape: &shape,
                level: 0.1,
                eps: 0.01,
            };
            let via_trait = model.q_r(&obs, &cand, &ctx, &mut trait_scratch);
            let via_fn =
                read_likelihood(&obs, &cand, &motif, &shape, 0.1, 0.01, &mut direct_scratch);
            assert_eq!(
                via_trait.to_bits(),
                via_fn.to_bits(),
                "q_r must equal read_likelihood for obs {obs_units}u / cand {cand_units}u"
            );
        }
    }

    /// An impure candidate exercises the `Σ_v` placement branch through the trait.
    #[test]
    fn q_r_matches_on_an_impure_candidate() {
        let motif = Motif::new(b"CA").unwrap();
        let shape = shape();
        let model = ClassicStutterModel;
        let mut trait_scratch = LikelihoodScratch::new();
        let mut direct_scratch = LikelihoodScratch::new();
        let cand = b"CACACATTCACA"; // (CA)3 TT (CA)2

        let ctx = ReadScoringContext {
            motif: &motif,
            shape: &shape,
            level: 0.1,
            eps: 0.01,
        };
        let via_trait = model.q_r(cand, cand, &ctx, &mut trait_scratch);
        let via_fn = read_likelihood(cand, cand, &motif, &shape, 0.1, 0.01, &mut direct_scratch);
        assert_eq!(via_trait.to_bits(), via_fn.to_bits());
    }
}
