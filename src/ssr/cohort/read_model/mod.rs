//! The swappable read-likelihood interface `Qᵣ(obs | cand)` (Stage-2 scoring
//! bake-off, plan `ssr_stutter_scoring_model_bakeoff.md` §3).
//!
//! `Qᵣ` is the probability an observed repeat tract arose from a candidate allele via
//! stutter (a whole-tract length change) plus base error (a composition change). It
//! drives the per-locus EM, so its shape sets genotyping accuracy and calibration.
//! Three credible models are implemented behind this one trait and baked off on
//! synthetic data with known truth; the production model is then chosen on measured
//! accuracy / calibration / robustness / speed.
//!
//! The bake-off (report `ssr_stutter_scoring_model_bakeoff_2026-06-24.md`) chose
//! **Model A (HipSTR)** as the production model: tied top genotype concordance with
//! Model B across the clean/out-of-frame/messy generative cases, but far better
//! calibrated and ~100× cheaper per `Qᵣ`. Model C (two-penalty pair-HMM) was eliminated
//! — even self-affinity-normalized it collapsed on out-of-frame reads — and its code is
//! removed. Model B is **kept behind the trait** as the bake-off's reference comparator.
//!
//! - [`HipstrModel`] (Model A): the production model — explicit in-frame + out-of-frame
//!   geometric stutter.
//! - [`ClassicStutterModel`] (Model B): the prior `Σ_Δ S_θ(Δ)·align_subst` code, now the
//!   test-only baseline the harness still scores A against.

mod classic;
mod hipstr;

pub(crate) use hipstr::HipstrModel;
// Model B is retained only as the bake-off's reference comparator (the harness scores A
// against it); production uses A, so B's handle is test-only.
#[cfg(test)]
pub(crate) use classic::ClassicStutterModel;

use crate::ssr::cohort::param_estimation::StutterShape;
use crate::ssr::types::Motif;

/// The per-call scoring inputs the genotyping EM computes for each
/// `(sample, locus, candidate)`: the locus motif, the per-locus stutter shape
/// `θ_locus`, the per-read stutter `level`, and the per-sample-group base error `ε`.
///
/// These are the quantities that genuinely vary call-to-call (the M-step refines
/// `shape`, `level` is `baseline + slope·len` re-weighted per candidate length, `ε`
/// is per sample group). The *model-frozen* parameters — HipSTR's out-of-frame
/// geometric, the two-penalty gaps — instead live on the model value itself, so every
/// model reads the same dynamic context while owning its own constants.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ReadScoringContext<'a> {
    /// The locus repeat motif (its `period()` is the unit size).
    pub(crate) motif: &'a Motif,
    /// `θ_locus` — *where* a slip lands (direction split + geometric decay).
    pub(crate) shape: &'a StutterShape,
    /// The per-read stutter `level` (≈ `P(Δ ≠ 0)`) for this candidate length.
    pub(crate) level: f64,
    /// The per-sample-group within-tract base error `ε`.
    pub(crate) eps: f64,
}

/// A swappable read-likelihood model: `Qᵣ(obs | cand)` under one treatment of stutter
/// (length change) and base error (composition change).
///
/// Implementations must be **pure** per call — a function of `obs`, `cand`, and the
/// [`ReadScoringContext`], with no hidden state beyond the reusable
/// [`Scratch`](Self::Scratch) — so the cohort stays byte-identical across thread
/// counts (the determinism contract).
pub(crate) trait ReadLikelihoodModel {
    /// Reusable per-read scratch (the placement-variant set, the DP buffers), so
    /// scoring a locus's reads allocates once. `Qᵣ` is ~75% of pileup self-time, so
    /// each model keeps its own scratch shape rather than a shared one.
    type Scratch: Default;

    /// `Qᵣ(obs | cand)` — the probability `obs` arose from `cand` under this model.
    fn q_r(
        &self,
        obs: &[u8],
        cand: &[u8],
        ctx: &ReadScoringContext,
        scratch: &mut Self::Scratch,
    ) -> f64;
}
