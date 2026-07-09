//! The per-locus genotype EM (arch `ssr_call_genotyping.md` §5, spec §4.2/§5.4,
//! implementation plan §3.6, C4).
//!
//! **Q-G2 — decided here: a slim SSR-specific EM, not a `posterior_engine` graft.**
//! The SSR genotype model is small and clean to express directly: genotypes are
//! size-`ploidy` multisets of candidate alleles, the prior is HWE with an IBD-`F`
//! autozygous branch, and the M-step is the per-candidate `G₀`-regularized
//! frequency update. Because `ε`, the stutter shape `θ`, and the level are *fixed*
//! for this milestone (supplied params), each genotype's **data** log-likelihood is
//! a one-time precompute and only `π` iterates — so the EM is a few cheap passes.
//! Bending the SNP engine (chain anchors, class pseudocounts, contamination) would
//! cost more than it saves. `ε` is frozen by the burn-in; the **stutter shape and rate
//! adapt per locus** (I1 + I2): an inner M-step attributes the locus's reads to their
//! called alleles and refits both the shape `θ_locus` (shrunk toward the frozen
//! `θ_period` seed) and a per-locus rate multiplier on the frozen group level (shrunk
//! toward 1) — so a thin locus stays at the frozen priors (no oscillation) — re-genotyping
//! until both settle. The group level keeps the length-dependence; the locus nudges the
//! overall rate. The soft per-read responsibility split is the deferred refinement.
//!
//! v1 is **diploid** (the simulator's ploidy); higher ploidy is a documented
//! follow-up (the genotype enumeration generalizes, the dosage priors need care).

use crate::ssr::cohort::allele_freq_prior::g0_pseudocounts;
use crate::ssr::cohort::attribution::{allele_responsibilities, nearest_called_by_sequence};
use crate::ssr::cohort::candidate_set::{Admission, CandidateSet};
use crate::ssr::cohort::em_init::LocusSeed;
use crate::ssr::cohort::likelihood::read_given_genotype;
use crate::ssr::cohort::pair_hmm::HmmScratch;
use crate::ssr::cohort::param_estimation::{
    DEFAULT_G0_FALLBACK_P, G0PseudocountDecay, ParamSet, PurityLevel, SlipProfile, StutterLevel,
    StutterShape, interruption_count,
};
use crate::ssr::cohort::read_model::{HipstrModel, ReadLikelihoodModel, ReadScoringContext};
use crate::ssr::cohort::rung_ladder::Rungs;
use crate::ssr::cohort::stutter::refine_theta_locus;
use crate::ssr::cohort::types::CohortLocus;
use smallvec::{SmallVec, smallvec};

/// EM controls (pinned in F2).
/// The per-locus **emission model**, selected by `PVC_SSR_EMIT_MODEL` (default
/// `Heuristic` → byte-identical). Orthogonal to [`EmCfg::marginalized_prior`], which
/// shapes the EM genotypes / `π` that *every* model consumes — any emission model
/// combines with either genotype prior (e.g. MARG + BIC, MARG + Freebayes).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EmitModel {
    /// Emit-iff-variable + allele-balance FP-control (the historical default).
    Heuristic,
    /// Model selection: emit iff the polymorphic model beats the monomorphic one by BIC.
    Bic,
    /// Freebayes-style SFS-prior marginal: emit iff `−10·log10 P(monomorphic)` clears a
    /// QUAL floor. See [`crate::ssr::cohort::freebayes_emit`].
    Freebayes,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EmCfg {
    /// Maximum EM iterations.
    pub(crate) max_iters: usize,
    /// Convergence tolerance on the largest `π` change between iterations.
    pub(crate) tol: f64,
    /// Outlier mixing weight `λ` (uniform junk term).
    pub(crate) lambda: f64,
    /// Seed inbreeding coefficient `F` (the autozygous-branch weight) used **only** by the
    /// [`run_locus_em`] convenience wrapper (tests / standalone). The production path
    /// (`driver::genotype_locus` → [`run_locus_em_with`]) supplies an explicit per-sample
    /// `F` frozen by the burn-in and never reads this field — so setting it and calling
    /// `run_locus_em_with` directly has no effect (review Mi11).
    pub(crate) inbreeding_f: f64,
    /// Maximum per-locus refinement rounds (I1 shape + I2 level). `0` keeps the frozen
    /// `θ_period` seed shape and group level (no per-locus adaptation).
    pub(crate) refit_max_rounds: usize,
    /// Pseudo-count weight shrinking `θ_locus` toward the frozen `θ_period` seed — the
    /// anti-oscillation knob: a thin locus cannot overcome it and stays at the seed.
    pub(crate) theta_shrink: f64,
    /// Convergence tolerance on the largest `θ_locus` coefficient change between rounds.
    pub(crate) theta_tol: f64,
    /// Pseudo-count weight (in slipped-read units) shrinking the per-locus stutter-rate
    /// multiplier toward the group level — the level analogue of `theta_shrink` (I2).
    pub(crate) level_shrink: f64,
    /// Convergence tolerance on the per-locus level multiplier between rounds.
    pub(crate) level_tol: f64,
    /// Genotype-prior model for the π-EM. `false` (default) = the plug-in
    /// HWE(π̂) prior; `true` = the marginalized Dirichlet-multinomial prior with
    /// leave-one-out cohort sharing (the improved SNP-path prior, seeded by the
    /// same mode-centred `G₀`). Off by default so SSR output is byte-identical
    /// until the benchmark (Phase 3.5) decides whether to flip it. See
    /// [`marginalized_genotype_log_priors`].
    pub(crate) marginalized_prior: bool,
    /// The emission model (`PVC_SSR_EMIT_MODEL`). Default `Heuristic`. Read inside
    /// [`run_locus_em_with`] only to decide whether to compute the freebayes SFS marginal
    /// (the count-vector enumeration is not free); the emit/drop *decision* is dispatched
    /// in `driver::genotype_locus`. Orthogonal to `marginalized_prior`.
    pub(crate) emit_model: EmitModel,
    /// Population-scaled diversity `θ` for the freebayes SFS prior (read only when
    /// `emit_model == Freebayes`). See [`crate::ssr::cohort::freebayes_emit::SFS_THETA`].
    pub(crate) sfs_theta: f64,
    /// (`PVC_SSR_NULL_FROM_HOMS`): estimate the per-locus stutter (shape θ *and* rate) from the
    /// plants that are unambiguously homozygous by a read-fraction rule — their off-modal reads
    /// are stutter by definition, with no genotype-call ambiguity — pooled across the whole
    /// cohort, then judge emission (**both** BIC and freebayes) against that genotype-free null,
    /// so an alt must clear the honest stutter. The EM's genotypes/`variable` are untouched
    /// (emission-only). Off by default → byte-identical.
    pub(crate) null_from_homs: bool,
    /// Min reads before a plant is eligible to be a clean-hom stutter donor for
    /// [`Self::null_from_homs`]. Defaults to [`DEFAULT_NULL_HOM_MIN_DEPTH`].
    pub(crate) null_hom_min_depth: u64,
    /// Modal-length read fraction at/above which a plant counts as clean-hom for
    /// [`Self::null_from_homs`]. Defaults to [`DEFAULT_NULL_HOM_FRAC`].
    pub(crate) null_hom_frac: f64,
}

/// Clean-hom stutter-donor thresholds for [`EmCfg::null_from_homs`] (`attribute_clean_homs`).
/// Mirror the silver-standard scorer (`benchmarks/ssr_tomato1/scripts/silver_standard.py`,
/// `MIN_DEPTH` / `HOM_FRAC`) so the caller's clean-hom donor set matches the benchmark that
/// evaluates it. Overridable per run via `PVC_SSR_NULL_HOM_MINDEPTH` / `PVC_SSR_NULL_HOM_FRAC`.
pub(crate) const DEFAULT_NULL_HOM_MIN_DEPTH: u64 = 4;
pub(crate) const DEFAULT_NULL_HOM_FRAC: f64 = 0.75;

impl EmCfg {
    pub(crate) fn dev_default() -> Self {
        Self {
            max_iters: 50,
            tol: 1e-6,
            lambda: 0.01,
            inbreeding_f: 0.0,
            refit_max_rounds: 3,
            theta_shrink: 50.0,
            theta_tol: 1e-3,
            level_shrink: 20.0,
            level_tol: 1e-3,
            marginalized_prior: false,
            emit_model: EmitModel::Heuristic,
            sfs_theta: crate::ssr::cohort::freebayes_emit::SFS_THETA,
            null_from_homs: false,
            null_hom_min_depth: DEFAULT_NULL_HOM_MIN_DEPTH,
            null_hom_frac: DEFAULT_NULL_HOM_FRAC,
        }
    }
}

/// One sample's genotype call.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SampleCall {
    /// The called allele copies as candidate indices (ascending; empty = no-call).
    pub(crate) allele_indices: Vec<usize>,
    /// The called allele lengths in repeat units (parallel to `allele_indices`).
    pub(crate) genotype_units: Vec<u16>,
    /// Posterior probability of the called genotype.
    pub(crate) posterior: f64,
    /// Phred genotype quality, capped at 99.
    pub(crate) gq: u8,
    /// Deconvolved per-called-allele read responsibility, parallel to `allele_indices`:
    /// `n_a = Σ_reads count · Qᵣ(obs | a) / Σ_{a'∈G} Qᵣ(obs | a')` for the called genotype
    /// `G`. This is the sequence-aware input to the allele-balance FP term (arch §3), which
    /// splits reads by *composition* (so a same-length het is deconvolved) rather than by
    /// tract length. Empty on a no-call and until P1.4 fills it (Phase 1's issue-2 fix).
    pub(crate) allele_support: SmallVec<[f64; 2]>,
}

impl SampleCall {
    pub(crate) fn no_call() -> Self {
        Self {
            allele_indices: Vec::new(),
            genotype_units: Vec::new(),
            posterior: 0.0,
            gq: 0,
            allele_support: SmallVec::new(),
        }
    }
}

/// The result of genotyping one locus.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct LocusCall {
    /// Per present sample (parallel to `locus.samples`).
    pub(crate) calls: Vec<SampleCall>,
    /// The converged candidate allele frequencies.
    pub(crate) pi: Vec<f64>,
    /// Per present sample: the posterior probability the sample is homozygous (the
    /// E1 `F` reduce reads this against the HWE expectation `Σ π_i²`).
    pub(crate) posterior_hom: Vec<f64>,
    /// The site admission verdict (drives the VCF FILTER).
    pub(crate) admit: Admission,
}

/// A diploid genotype as an ordered pair of candidate indices `i ≤ j`.
#[derive(Debug, Clone, Copy)]
struct Genotype {
    i: usize,
    j: usize,
}

fn enumerate_diploid_genotypes(k: usize) -> Vec<Genotype> {
    let mut genotypes = Vec::with_capacity(k * (k + 1) / 2);
    for i in 0..k {
        for j in i..k {
            genotypes.push(Genotype { i, j });
        }
    }
    genotypes
}

/// HWE genotype prior with an IBD-`F` autozygous branch.
fn genotype_prior(g: Genotype, pi: &[f64], f: f64) -> f64 {
    if g.i == g.j {
        f * pi[g.i] + (1.0 - f) * pi[g.i] * pi[g.i]
    } else {
        (1.0 - f) * 2.0 * pi[g.i] * pi[g.j]
    }
}

// ------------------------------------------------------------------------
// Marginalized Dirichlet-multinomial genotype prior (Phase 3.2).
//
// The improved SNP-path prior, ported to SSR: instead of plugging a point
// frequency `π` into HWE ([`genotype_prior`]), average the genotype probability
// over the frequency uncertainty using the Dirichlet-multinomial, seeded by the
// **mode-centred `G₀`** concentration (§Q1 — SSR keeps its own hypervariable
// seed, it does NOT inherit the SNP reference-is-common seed). Wright-`F` IBD is
// mixed on top exactly as the SNP engine does. This is additive: the plug-in path
// above is untouched, and the SNP engine is untouched (this mirrors its
// `e_step_cohort_loo`, it does not extract from it — protecting the SNP SIMD
// gains). Gated on behind an `EmCfg` toggle in a later step.
// ------------------------------------------------------------------------

/// The flat Dirichlet-multinomial inputs for a diploid candidate set of size `k`:
/// `genotype_allele_counts` (row-major `n_genotypes × k`) and the per-genotype
/// `log_multinomial_coeffs` (`ln 1 = 0` for a homozygote, `ln 2` for a het), in
/// `genotypes` order. These are exactly what
/// [`crate::genetics::dirichlet_multinomial_log_priors`] consumes.
fn diploid_dm_inputs(genotypes: &[Genotype], k: usize) -> (Vec<u32>, Vec<f64>) {
    let mut counts = vec![0u32; genotypes.len() * k];
    let mut log_coeffs = vec![0.0_f64; genotypes.len()];
    for (g_idx, g) in genotypes.iter().enumerate() {
        // A candidate index `≥ k` would write into a later genotype's row and
        // silently corrupt it (the flat buffer stays in-bounds), so guard loudly
        // rather than let a shape bug escape as a wrong prior.
        debug_assert!(
            g.i < k && g.j < k,
            "genotype ({}, {}) out of range k={k}",
            g.i,
            g.j
        );
        let base = g_idx * k;
        if g.i == g.j {
            counts[base + g.i] = 2;
            // ln(2!/2!) = ln 1 = 0 (already).
        } else {
            counts[base + g.i] = 1;
            counts[base + g.j] = 1;
            log_coeffs[g_idx] = std::f64::consts::LN_2; // ln(2!/(1!1!)) = ln 2
        }
    }
    (counts, log_coeffs)
}

/// Leave-one-out concentration `α'_s = G₀ + max(0, E[cohort copies] − E[own copies])`
/// for one sample, mirroring the SNP engine's `e_step_cohort_loo`. `g0` is the
/// mode-centred seed; `cohort_expected` / `own_expected` are posterior-weighted
/// allele copies over the whole cohort and over this sample. The `max(0, …)`
/// guards float noise on `cohort − own` (own is one non-negative addend of the
/// total, so the true difference is `≥ 0`). All three slices are length `k`.
fn leave_one_out_alpha(g0: &[f64], cohort_expected: &[f64], own_expected: &[f64]) -> Vec<f64> {
    debug_assert_eq!(g0.len(), cohort_expected.len());
    debug_assert_eq!(g0.len(), own_expected.len());
    g0.iter()
        .zip(cohort_expected)
        .zip(own_expected)
        .map(|((&g0_a, &total_a), &own_a)| g0_a + (total_a - own_a).max(0.0))
        .collect()
}

/// Per-genotype log-prior under the marginalized Dirichlet-multinomial with a
/// Wright-`F` IBD branch, for a diploid locus of `k` candidates. `alpha` is the
/// per-candidate concentration (the sample's `α'_s` from [`leave_one_out_alpha`]);
/// `f` is the sample's inbreeding coefficient. Returns one value per genotype in
/// `genotypes` order, log-priors up to a shared additive constant (softmax-ready)
/// — matching the SNP engine's convention and [`dirichlet_multinomial_log_priors`].
///
/// The IBD mixture mirrors the SNP `e_step`: a homozygote `(i, i)` is
/// `logsumexp((1−f)·DM,  f·(α_i/Σα))`; a heterozygote is `(1−f)·DM`.
fn marginalized_genotype_log_priors(
    genotypes: &[Genotype],
    k: usize,
    alpha: &[f64],
    f: f64,
) -> Vec<f64> {
    debug_assert_eq!(alpha.len(), k);
    let (counts, log_coeffs) = diploid_dm_inputs(genotypes, k);
    let log_indep =
        crate::genetics::dirichlet_multinomial_log_priors(&counts, &log_coeffs, k, alpha);

    let sum_alpha: f64 = alpha.iter().sum();
    let log_sum_alpha = sum_alpha.ln();
    let log_one_minus_f = (1.0 - f).ln(); // -∞ at f = 1
    let log_f = f.ln(); // -∞ at f = 0

    genotypes
        .iter()
        .zip(&log_indep)
        .map(|(g, &log_indep_g)| {
            if g.i == g.j {
                // IBD marginal allele log-frequency log(α_i / Σα).
                let log_p_effective_i = alpha[g.i].ln() - log_sum_alpha;
                log_sum_exp(&[log_one_minus_f + log_indep_g, log_f + log_p_effective_i])
            } else {
                log_one_minus_f + log_indep_g
            }
        })
        .collect()
}

/// The `G₀` decay for the locus period (the shared coded fallback if absent — review M3).
fn period_decay(params: &ParamSet, period: usize) -> G0PseudocountDecay {
    const FALLBACK: G0PseudocountDecay = G0PseudocountDecay {
        p: DEFAULT_G0_FALLBACK_P,
    };
    params
        .pseudocount_decay_per_loci_group
        .get(&(period as u8))
        .copied()
        .unwrap_or(FALLBACK)
}

/// `ε` and the per-length stutter level for a present sample (by its sample group).
/// `level_per_group` is supplied separately so the E1 outer loop can refit it without
/// rebuilding the params.
///
/// PANIC-FREE: decision E guarantees every present sample resolves to a frozen sample
/// group (`build_param_set` hard-errors `UnresolvedSamples` otherwise), and the per-group
/// vectors are sized per group — so all three lookups are in range. We index loudly rather
/// than fabricate group-0 / default-ε / default-level chemistry, which would silently
/// genotype a sample on wrong parameters under a broken invariant — exactly the
/// no-silent-default failure `UnresolvedSamples` exists to prevent (review M2).
fn sample_chemistry(
    locus: &CohortLocus,
    present_k: usize,
    params: &ParamSet,
    level_per_group: &[StutterLevel],
) -> (f64, StutterLevel) {
    let global = locus.present[present_k] as usize;
    let group = *params
        .group_of_sample
        .get(global)
        .expect("decision E: every present sample has a frozen sample group");
    let eps = params
        .error_per_sample_group
        .get(group.0 as usize)
        .expect("every sample group has a frozen ε")
        .0;
    let level = *level_per_group
        .get(group.0 as usize)
        .expect("every sample group has a frozen stutter level");
    (eps, level)
}

/// The per-read stutter **level** for one candidate: the group level line
/// `baseline + slope·units`, scaled by the per-locus rate `multiplier` and — the Phase-2
/// addition — the candidate's **purity** factor `purity.level_factor(interruptions)`, clamped
/// to `[0, 1]`. A pure allele (0 interruptions) or a neutral `purity` (`per_interruption_factor
/// = 1`) leaves the pre-Phase-2 arithmetic exactly. Shared by the data likelihood and the
/// allele-balance deconvolution so both score on the identical level (arch §3).
fn candidate_level(
    level: StutterLevel,
    units: u16,
    multiplier: f64,
    purity: &PurityLevel,
    interruptions: u32,
) -> f64 {
    ((level.baseline + level.slope * f64::from(units))
        * multiplier
        * purity.level_factor(interruptions))
    .clamp(0.0, 1.0)
}

/// Distinct observed sequences across the cohort at this locus (`D`, the outlier
/// term's denominator).
fn distinct_sequence_count(locus: &CohortLocus) -> usize {
    let mut seqs: Vec<&[u8]> = locus
        .samples
        .iter()
        .flat_map(|ev| ev.seq_counts.iter().map(|(s, _)| s.as_ref()))
        .collect();
    seqs.sort_unstable();
    seqs.dedup();
    seqs.len().max(1)
}

/// log-sum-exp of a slice (stable).
fn log_sum_exp(xs: &[f64]) -> f64 {
    let max = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    max + xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln()
}

/// Cohort marginal read log-likelihoods for the per-locus **emission** model choice
/// (spec `ssr_cohort_recurrence_aware_gates.md`, model-selection). Computed at the settled
/// read model so the emission decision reuses exactly the likelihoods that produced the calls.
///
/// - [`Self::ln_marginal`] — Σ over present samples of `log Σ_G P(G | π̂, F) P(reads | G)`, the
///   maximised marginal likelihood under the **polymorphic** candidate model (the alt alleles
///   segregate; genotypes integrated out, frequency `π̂` fit by the EM).
/// - [`Self::ln_monomorphic`] — Σ of `P(reads | hom aa)` for the single fixed allele `a` that
///   best explains the cohort (no segregating alt; stutter alone accounts for off-length reads).
///
/// A BIC / likelihood-ratio emission gate compares the two: the extra allele must earn its
/// parameter cost before the locus is emitted, which rejects a systematic stutter shoulder
/// (it barely lifts the likelihood) while keeping a real segregating allele (it lifts it a lot).
/// The freebayes-style SFS marginal ([`Self::freebayes_ln_p_mono`]) is carried alongside the
/// BIC pair so **both** emission models decide from the same one read-likelihood pass
/// (`data_ll`) — only the decision differs, never the reads. It is `Some` only when
/// `emit_model == Freebayes` (the count-vector marginal is computed on demand).
#[derive(Debug, Clone, Copy)]
pub(crate) struct EmissionEvidence {
    pub(crate) ln_marginal: f64,
    pub(crate) ln_monomorphic: f64,
    /// Freebayes-style `ln P(monomorphic | data)` (≤ 0) from the SFS-prior marginal over the
    /// same `data_ll`; `Some` only for `EmitModel::Freebayes`, `None` otherwise / filtered loci.
    pub(crate) freebayes_ln_p_mono: Option<f64>,
}

/// Genotype one locus on the supplied parameters.
///
/// PLOIDY CONTRACT: v1 is diploid-only and **panics** on `ploidy != 2`. This is a
/// deliberate loud failure (the no-silent-fallback invariant, plan §5) — a
/// non-diploid call in a diploid build is a programming error, not something to
/// quietly no-call. Ploidy generalization (genotype enumeration + dosage priors) is
/// a documented follow-up.
pub(crate) fn run_locus_em(
    locus: &CohortLocus,
    rungs: &Rungs,
    candidates: &CandidateSet,
    params: &ParamSet,
    seed: &LocusSeed,
    ploidy: u8,
    cfg: &EmCfg,
) -> LocusCall {
    let f_per_present = vec![cfg.inbreeding_f; locus.samples.len()];
    run_locus_em_with(
        locus,
        rungs,
        candidates,
        params,
        seed,
        ploidy,
        cfg,
        &f_per_present,
        &params.level_seed,
        &HipstrModel,
    )
    .0
}

/// Genotype one locus with an explicit per-sample inbreeding `F` and per-group level
/// (the E1 outer loop drives these); also emits each sample's posterior homozygosity.
// The frozen params + dev-config + the per-sample F / per-group level are threaded
// explicitly (not bundled) so the per-locus call stays a pure function of its inputs — the
// byte-identity contract. Arg count is intentional.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_locus_em_with<M: ReadLikelihoodModel>(
    locus: &CohortLocus,
    rungs: &Rungs,
    candidates: &CandidateSet,
    params: &ParamSet,
    seed: &LocusSeed,
    ploidy: u8,
    cfg: &EmCfg,
    f_per_present: &[f64],
    level_per_group: &[StutterLevel],
    model: &M,
) -> (LocusCall, EmissionEvidence) {
    assert_eq!(ploidy, 2, "C4 EM is diploid-only (v1)");
    let n_present = locus.samples.len();

    // A no-call locus still reports REF + the no-call samples; no emission evidence (a
    // filtered locus is emitted on its FILTER, not on the model choice).
    if candidates.admit != Admission::Pass {
        return (
            LocusCall {
                calls: vec![SampleCall::no_call(); n_present],
                pi: seed.pi0.clone(),
                posterior_hom: vec![0.0; n_present],
                admit: candidates.admit,
            },
            EmissionEvidence {
                ln_marginal: 0.0,
                ln_monomorphic: 0.0,
                freebayes_ln_p_mono: None,
            },
        );
    }

    let period = rungs.period();
    let k = candidates.alleles.len();
    let genotypes = enumerate_diploid_genotypes(k);
    let distinct = distinct_sequence_count(locus);
    let g0 = g0_pseudocounts(candidates, rungs, &period_decay(params, period));

    // Candidate lengths in repeat units (for the per-candidate level and the call).
    let cand_units: Vec<u16> = candidates
        .alleles
        .iter()
        .map(|a| (a.len() / period) as u16)
        .collect();
    // Per-candidate interruption count (Phase 2): the purity that scales its stutter level. A
    // pure allele scores 0 → factor 1 → pre-Phase-2 arithmetic (spec §6.3).
    let cand_interruptions: Vec<u32> = candidates
        .alleles
        .iter()
        .map(|a| interruption_count(a, period))
        .collect();

    // The iteration-invariant locus shaping, built once and reused across the refit rounds
    // (review Mi2).
    let locus_model = LocusModel {
        locus,
        candidates,
        cand_units: &cand_units,
        cand_interruptions: &cand_interruptions,
        genotypes: &genotypes,
        distinct,
    };

    let mut scratch = M::Scratch::default();

    // Per-locus refinement (I1 shape + I2 level): genotype under the frozen `θ_period`
    // seed and group level, then refit a per-locus shape `θ_locus` (shrunk toward the
    // seed) and a per-locus stutter-rate multiplier on the group level (shrunk toward 1)
    // from the called alleles' slips, re-genotyping until both settle. With few / no slips
    // both collapse to the frozen priors (no oscillation); each round recomputes the data
    // log-likelihood for the new shape + rate.
    let mut theta = seed.theta0;
    let mut level_multiplier = 1.0;
    let mut data_ll = compute_data_ll(
        &locus_model,
        params,
        level_per_group,
        &theta,
        level_multiplier,
        cfg.lambda,
        model,
        &mut scratch,
    );
    let (mut pi, mut calls, mut posterior_hom) = genotype_pass(
        &data_ll,
        &genotypes,
        &g0,
        &seed.pi0,
        f_per_present,
        k,
        &cand_units,
        cfg,
    );

    for _ in 0..cfg.refit_max_rounds {
        let fit = attribute_locus(locus, &calls, period, candidates, params, level_per_group);
        let new_theta = refine_theta_locus(&fit.profile, &seed.theta0, cfg.theta_shrink);
        let new_level_multiplier = refit_level_multiplier(&fit, cfg.level_shrink);
        if shapes_close(&new_theta, &theta, cfg.theta_tol)
            && (new_level_multiplier - level_multiplier).abs() < cfg.level_tol
        {
            break;
        }
        theta = new_theta;
        level_multiplier = new_level_multiplier;
        data_ll = compute_data_ll(
            &locus_model,
            params,
            level_per_group,
            &theta,
            level_multiplier,
            cfg.lambda,
            model,
            &mut scratch,
        );
        let (p, c, ph) = genotype_pass(
            &data_ll,
            &genotypes,
            &g0,
            &seed.pi0,
            f_per_present,
            k,
            &cand_units,
            cfg,
        );
        pi = p;
        calls = c;
        posterior_hom = ph;
    }

    // Deconvolve each call's per-allele read support (sequence-aware allele balance, arch §3),
    // using the same settled shape / level that produced these calls.
    fill_allele_support(
        &mut calls,
        &locus_model,
        params,
        level_per_group,
        &theta,
        level_multiplier,
        model,
        &mut scratch,
    );

    // Emission evidence for the model-selection gate, from the settled `data_ll` + `pi`. The
    // BIC pair (cheap) is always computed; the freebayes SFS marginal (a count-vector
    // enumeration) only when that model is selected — both from one `data_ll` so the emission
    // choice, not the reads, is all that differs between models (the fairness invariant).
    //
    // `null_from_homs`: judge emission against a genotype-free per-locus stutter null learned
    // from the clearly-hom plants. Re-estimate θ + rate from their pooled off-modal reads
    // (shrunk toward the frozen group prior when the clean-hom pool is thin), recompute the read
    // likelihoods under that honest null, and derive the emission evidence — the marginal *and*
    // the monomorphic term — from it, so an alt must clear the real stutter. The EM's genotypes
    // are untouched, so `variable`/GTs/recall-of-candidates are unchanged (emission-only).
    let data_ll_clean_hom_null = if cfg.null_from_homs {
        let clean_hom_fit = attribute_clean_homs(
            locus,
            period,
            params,
            level_per_group,
            cfg.null_hom_min_depth,
            cfg.null_hom_frac,
        );
        let clean_hom_theta =
            refine_theta_locus(&clean_hom_fit.profile, &seed.theta0, cfg.theta_shrink);
        let clean_hom_level_multiplier = refit_level_multiplier(&clean_hom_fit, cfg.level_shrink);
        Some(compute_data_ll(
            &locus_model,
            params,
            level_per_group,
            &clean_hom_theta,
            clean_hom_level_multiplier,
            cfg.lambda,
            model,
            &mut scratch,
        ))
    } else {
        None
    };
    // Borrow (no clone) the honest-null matrix when the toggle is on, else the settled
    // genotyping one — with the toggle off this is exactly `&data_ll` (byte-identical).
    let data_ll_emit: &[Vec<f64>] = data_ll_clean_hom_null.as_deref().unwrap_or(&data_ll);
    let mut evidence = emission_evidence(data_ll_emit, &genotypes, &pi, f_per_present, k);
    if cfg.emit_model == EmitModel::Freebayes {
        evidence.freebayes_ln_p_mono = Some(crate::ssr::cohort::freebayes_emit::ln_p_monomorphic(
            data_ll_emit,
            &pi,
            f_per_present,
            k,
            cfg.sfs_theta,
        ));
    }

    (
        LocusCall {
            calls,
            pi,
            posterior_hom,
            admit: candidates.admit,
        },
        evidence,
    )
}

/// Compute the [`EmissionEvidence`] (polymorphic marginal vs best monomorphic) from the settled
/// per-sample per-genotype data log-likelihoods `data_ll`, the fit frequency `pi`, and per-sample
/// `F`. Uses the plug-in HWE genotype prior (`genotype_prior`) — the emission choice is a
/// separate axis from which genotype prior the calls used.
fn emission_evidence(
    data_ll: &[Vec<f64>],
    genotypes: &[Genotype],
    pi: &[f64],
    f_per_present: &[f64],
    k: usize,
) -> EmissionEvidence {
    // Polymorphic model: marginal likelihood integrating genotypes under the fit frequency.
    let mut ln_marginal = 0.0;
    for (sample_ll, &f) in data_ll.iter().zip(f_per_present) {
        let log_joint: Vec<f64> = genotypes
            .iter()
            .zip(sample_ll)
            .map(|(g, ll)| genotype_prior(*g, pi, f).ln() + ll)
            .collect();
        ln_marginal += log_sum_exp(&log_joint);
    }
    // Monomorphic model: every sample homozygous for the single best-fitting allele.
    let mut ln_monomorphic = f64::NEG_INFINITY;
    for a in 0..k {
        if let Some(gi) = genotypes.iter().position(|g| g.i == a && g.j == a) {
            let s: f64 = data_ll.iter().map(|sll| sll[gi]).sum();
            if s > ln_monomorphic {
                ln_monomorphic = s;
            }
        }
    }
    EmissionEvidence {
        ln_marginal,
        ln_monomorphic,
        freebayes_ln_p_mono: None,
    }
}

/// The iteration-invariant per-locus shaping `compute_data_ll` reads: the locus and its
/// candidate set, the candidate lengths in repeat units, the enumerated genotypes, and the
/// distinct-observation count. Built once in [`run_locus_em_with`] and reused across every
/// `θ_locus` / level-multiplier refit round, so the precompute call stays a pure function
/// of `(model, chemistry, shape, rate)` (the byte-identity contract) instead of an 11-arg
/// positional list (review Mi2).
struct LocusModel<'a> {
    locus: &'a CohortLocus,
    candidates: &'a CandidateSet,
    cand_units: &'a [u16],
    cand_interruptions: &'a [u32],
    genotypes: &'a [Genotype],
    distinct: usize,
}

/// Each sample's per-genotype **data** log-likelihood `data_ll[s][g]` under stutter
/// shape `theta` and per-locus rate multiplier `level_multiplier` on the group level (both
/// constant across the π iterations; recomputed when `θ_locus` / the multiplier change).
// The model + its scratch join the (already-bundled) locus shaping, chemistry, and
// per-round shape/level — all distinct inputs threaded explicitly so the precompute
// stays a pure function (the byte-identity contract). One over the lint's bound.
#[allow(clippy::too_many_arguments)]
fn compute_data_ll<M: ReadLikelihoodModel>(
    locus_model: &LocusModel,
    params: &ParamSet,
    level_per_group: &[StutterLevel],
    theta: &StutterShape,
    level_multiplier: f64,
    lambda: f64,
    model: &M,
    scratch: &mut M::Scratch,
) -> Vec<Vec<f64>> {
    let LocusModel {
        locus,
        candidates,
        cand_units,
        cand_interruptions,
        genotypes,
        distinct,
    } = *locus_model;
    let k = candidates.alleles.len();
    let mut data_ll: Vec<Vec<f64>> = Vec::with_capacity(locus.samples.len());
    for (k_present, evidence) in locus.samples.iter().enumerate() {
        let (eps, level) = sample_chemistry(locus, k_present, params, level_per_group);
        // Per-candidate level (group line × per-locus multiplier × purity) — independent of the
        // read, so compute it once per sample rather than per (read × candidate).
        let cand_levels: Vec<f64> = (0..k)
            .map(|c| {
                candidate_level(
                    level,
                    cand_units[c],
                    level_multiplier,
                    &params.purity_level,
                    cand_interruptions[c],
                )
            })
            .collect();
        // Per (distinct obs, candidate) Qᵣ — the read likelihood matrix.
        let obs_qr: Vec<(u32, Vec<f64>)> = evidence
            .seq_counts
            .iter()
            .map(|(obs, count)| {
                let row = (0..k)
                    .map(|c| {
                        let ctx = ReadScoringContext {
                            motif: &locus.motif,
                            shape: theta,
                            level: cand_levels[c],
                            eps,
                        };
                        model.q_r(obs, &candidates.alleles[c], &ctx, scratch)
                    })
                    .collect::<Vec<_>>();
                (*count, row)
            })
            .collect();

        let sample_ll = genotypes
            .iter()
            .map(|g| {
                obs_qr
                    .iter()
                    .map(|(count, row)| {
                        let p = read_given_genotype(row, &[g.i, g.j], lambda, distinct);
                        *count as f64 * p.ln()
                    })
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        data_ll.push(sample_ll);
    }
    data_ll
}

/// The π loop: E-step posteriors over genotypes → `G₀`-regularized M-step on the allele
/// frequencies, to convergence (or `max_iters`).
fn run_pi_em(
    data_ll: &[Vec<f64>],
    genotypes: &[Genotype],
    g0: &[f64],
    pi0: &[f64],
    f_per_present: &[f64],
    cfg: &EmCfg,
) -> Vec<f64> {
    let mut pi = pi0.to_vec();
    for _ in 0..cfg.max_iters {
        let mut expected = g0.to_vec();
        for (k_present, sample_ll) in data_ll.iter().enumerate() {
            let f = f_per_present[k_present];
            let log_joint: Vec<f64> = genotypes
                .iter()
                .zip(sample_ll)
                .map(|(g, ll)| genotype_prior(*g, &pi, f).ln() + ll)
                .collect();
            let norm = log_sum_exp(&log_joint);
            for (g, lj) in genotypes.iter().zip(&log_joint) {
                let post = (lj - norm).exp();
                expected[g.i] += post;
                expected[g.j] += post;
            }
        }
        let total: f64 = expected.iter().sum();
        let new_pi: Vec<f64> = expected.iter().map(|e| e / total).collect();
        let delta = new_pi
            .iter()
            .zip(&pi)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        pi = new_pi;
        if delta < cfg.tol {
            break;
        }
    }
    pi
}

/// Final E-step → per sample: the MAP genotype + GQ + posterior homozygosity.
fn final_calls(
    data_ll: &[Vec<f64>],
    genotypes: &[Genotype],
    pi: &[f64],
    f_per_present: &[f64],
    cand_units: &[u16],
) -> (Vec<SampleCall>, Vec<f64>) {
    let mut calls = Vec::with_capacity(data_ll.len());
    let mut posterior_hom = Vec::with_capacity(data_ll.len());
    for (k_present, sample_ll) in data_ll.iter().enumerate() {
        let f = f_per_present[k_present];
        let log_joint: Vec<f64> = genotypes
            .iter()
            .zip(sample_ll)
            .map(|(g, ll)| genotype_prior(*g, pi, f).ln() + ll)
            .collect();
        let norm = log_sum_exp(&log_joint);
        // PANIC-FREE: a PASS locus enumerates ≥1 genotype (`log_joint` non-empty), and every
        // entry is finite-or-−∞ — `genotype_prior(...).ln()` is finite-or-−∞ and `ll` is a
        // sum of `count·ln(p)` with `p ≥ lambda/D > 0`, so no entry is NaN. `total_cmp` is a
        // total order over f64 (NaN-safe), so the argmax cannot panic (review M4).
        let (best, best_lj) = log_joint
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .expect("a PASS locus enumerates ≥1 genotype");
        let posterior = (best_lj - norm).exp();
        let p_hom: f64 = genotypes
            .iter()
            .zip(&log_joint)
            .filter(|(g, _)| g.i == g.j)
            .map(|(_, lj)| (lj - norm).exp())
            .sum();
        let g = genotypes[best];
        calls.push(SampleCall {
            allele_indices: vec![g.i, g.j],
            genotype_units: vec![cand_units[g.i], cand_units[g.j]],
            posterior,
            gq: phred_gq(posterior),
            // Filled in P1.4 (sequence-aware allele balance); empty here in Phase-1 types.
            allele_support: SmallVec::new(),
        });
        posterior_hom.push(p_hom);
    }
    (calls, posterior_hom)
}

/// The converged state of the marginalized π-EM, carried into
/// [`final_calls_marginalized`] so it can rebuild each sample's leave-one-out
/// concentration `α'_s = G₀ + (E[cohort] − E[own_s])`.
struct MarginalizedFit {
    /// Reported allele frequencies `π = E[cohort copies] / Σ` (raw copies, no `G₀`).
    pi: Vec<f64>,
    /// Converged posterior-weighted allele copies over the whole cohort.
    expected: Vec<f64>,
    /// Per-sample posterior-weighted copies from the final iteration (the
    /// leave-one-out subtrahend). Length `n_samples`, each entry length `k`.
    own_counts: Vec<Vec<f64>>,
}

/// The π loop under the **marginalized Dirichlet-multinomial** prior with
/// leave-one-out cohort sharing — the SSR analogue of the SNP engine's
/// `e_step_cohort_loo`. The first iteration (0-indexed `iteration == 0`) uses a
/// flat prior (`α = G₀`, the seed the whole cohort shares, no leave-one-out term
/// yet); every later iteration gives each sample its own
/// `α'_s = G₀ + (E[cohort] − E[own_s])`. Convergence is tested
/// on the **driver** — the change in the posterior-weighted allele copies
/// `E[cohort]` per chromosome — not on the reported `π` (the Phase-0 discipline:
/// `π` is a `G₀`-rescaled readout, `E[cohort]` is what the leave-one-out prior
/// actually feeds back). Diploid.
fn run_pi_em_marginalized(
    data_ll: &[Vec<f64>],
    genotypes: &[Genotype],
    g0: &[f64],
    f_per_present: &[f64],
    k: usize,
    cfg: &EmCfg,
) -> MarginalizedFit {
    let n_samples = data_ll.len();
    // A PASS locus always has ≥1 present sample. Guard the invariant loudly: with
    // `n_samples == 0` the driver delta and `π` readout below would divide by zero
    // and emit NaN (the plug-in `run_pi_em` seeds `expected = g0` and stays finite;
    // this path seeds raw zeros).
    debug_assert!(n_samples > 0, "marginalized π-EM needs ≥1 present sample");
    let chromosomes = 2.0 * n_samples as f64; // diploid
    let mut expected = vec![0.0_f64; k];
    let mut own_counts = vec![vec![0.0_f64; k]; n_samples];

    for iteration in 0..cfg.max_iters {
        let mut next_expected = vec![0.0_f64; k];
        let mut next_own = vec![vec![0.0_f64; k]; n_samples];
        for (s, sample_ll) in data_ll.iter().enumerate() {
            let f = f_per_present[s];
            // Flat seed on the first iteration (no cohort counts yet); LOO after.
            let alpha = if iteration == 0 {
                g0.to_vec()
            } else {
                leave_one_out_alpha(g0, &expected, &own_counts[s])
            };
            let log_prior = marginalized_genotype_log_priors(genotypes, k, &alpha, f);
            let log_joint: Vec<f64> = log_prior
                .iter()
                .zip(sample_ll)
                .map(|(lp, ll)| lp + ll)
                .collect();
            let norm = log_sum_exp(&log_joint);
            for (g, lj) in genotypes.iter().zip(&log_joint) {
                let post = (lj - norm).exp();
                next_expected[g.i] += post;
                next_expected[g.j] += post;
                next_own[s][g.i] += post;
                next_own[s][g.j] += post;
            }
        }
        // Driver delta on the per-chromosome allele copies. Never converge on the
        // flat first iteration (its counts reflect the seed, not the LOO prior
        // that produces the emitted calls) — the SSR analogue of the SNP guard.
        let delta = next_expected
            .iter()
            .zip(&expected)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max)
            / chromosomes;
        expected = next_expected;
        own_counts = next_own;
        if iteration > 0 && delta < cfg.tol {
            break;
        }
    }

    let total: f64 = expected.iter().sum();
    let pi = expected.iter().map(|e| e / total).collect();
    MarginalizedFit {
        pi,
        expected,
        own_counts,
    }
}

/// Final E-step under the converged marginalized prior → per sample the MAP
/// genotype + GQ + posterior homozygosity, rebuilding each sample's `α'_s` from
/// the converged [`MarginalizedFit`]. Mirrors [`final_calls`] (the plug-in tail),
/// differing only in the per-sample prior.
fn final_calls_marginalized(
    data_ll: &[Vec<f64>],
    genotypes: &[Genotype],
    g0: &[f64],
    fit: &MarginalizedFit,
    f_per_present: &[f64],
    k: usize,
    cand_units: &[u16],
) -> (Vec<SampleCall>, Vec<f64>) {
    let mut calls = Vec::with_capacity(data_ll.len());
    let mut posterior_hom = Vec::with_capacity(data_ll.len());
    for (s, sample_ll) in data_ll.iter().enumerate() {
        let f = f_per_present[s];
        let alpha = leave_one_out_alpha(g0, &fit.expected, &fit.own_counts[s]);
        let log_prior = marginalized_genotype_log_priors(genotypes, k, &alpha, f);
        let log_joint: Vec<f64> = log_prior
            .iter()
            .zip(sample_ll)
            .map(|(lp, ll)| lp + ll)
            .collect();
        let norm = log_sum_exp(&log_joint);
        // PANIC-FREE, as in `final_calls`: a PASS locus enumerates ≥1 genotype so
        // `log_joint` is non-empty, and every entry is finite-or-−∞, so `total_cmp`
        // (a NaN-safe total order) cannot fail the argmax.
        let (best, best_lj) = log_joint
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .expect("a PASS locus enumerates ≥1 genotype");
        let posterior = (best_lj - norm).exp();
        let p_hom: f64 = genotypes
            .iter()
            .zip(&log_joint)
            .filter(|(g, _)| g.i == g.j)
            .map(|(_, lj)| (lj - norm).exp())
            .sum();
        let g = genotypes[best];
        calls.push(SampleCall {
            allele_indices: vec![g.i, g.j],
            genotype_units: vec![cand_units[g.i], cand_units[g.j]],
            posterior,
            gq: phred_gq(posterior),
            allele_support: SmallVec::new(),
        });
        posterior_hom.push(p_hom);
    }
    (calls, posterior_hom)
}

/// One genotyping pass — the π-EM to convergence plus the final per-sample calls
/// — dispatching on [`EmCfg::marginalized_prior`]. The plug-in branch is exactly
/// the previous `run_pi_em` + `final_calls` sequence (byte-identical at the
/// default); the marginalized branch runs the leave-one-out Dirichlet-multinomial
/// prior. `k` is the candidate count.
#[allow(clippy::too_many_arguments)]
fn genotype_pass(
    data_ll: &[Vec<f64>],
    genotypes: &[Genotype],
    g0: &[f64],
    pi0: &[f64],
    f_per_present: &[f64],
    k: usize,
    cand_units: &[u16],
    cfg: &EmCfg,
) -> (Vec<f64>, Vec<SampleCall>, Vec<f64>) {
    if cfg.marginalized_prior {
        let fit = run_pi_em_marginalized(data_ll, genotypes, g0, f_per_present, k, cfg);
        let (calls, posterior_hom) =
            final_calls_marginalized(data_ll, genotypes, g0, &fit, f_per_present, k, cand_units);
        (fit.pi, calls, posterior_hom)
    } else {
        let pi = run_pi_em(data_ll, genotypes, g0, pi0, f_per_present, cfg);
        let (calls, posterior_hom) =
            final_calls(data_ll, genotypes, &pi, f_per_present, cand_units);
        (pi, calls, posterior_hom)
    }
}

/// The largest per-locus stutter-rate multiplier the refit may return (the
/// `level_multiplier` clamp; the resulting per-read level is `[0,1]`-clamped anyway, but
/// this bounds the multiplier itself against a degenerate `expected ≈ 0` denominator).
const LEVEL_MULT_MAX: f64 = 10.0;

/// The per-locus slip attribution: the sufficient statistics for both per-locus refits
/// (I1 shape, I2 level), built by hard-attributing each read to its nearest called allele.
#[derive(Debug, Default)]
struct LocusSlipFit {
    /// Signed slip-magnitude counts → the `θ_locus` shape M-step (I1).
    profile: SlipProfile,
    /// Total slipped reads (`Δ ≠ 0`) → the level numerator (I2).
    slipped: u64,
    /// Reads expected to slip under the *group* level at their parent length → the level
    /// denominator (I2): the multiplier is `observed / expected`.
    expected_slipped: f64,
}

/// IDEA B: pool the stutter sufficient statistics from the plants that are *unambiguously
/// homozygous* by a read-fraction rule (≥ `min_depth` reads, ≥ `hom_frac` of them at one
/// length), keyed to each donor's own modal length. Their off-modal reads are stutter by
/// definition — no genotype call, no allele-ambiguity, and no circularity with the caller we
/// are trying to improve. Pooling across the whole cohort makes the estimate robust even at
/// ~3 reads/plant. Ambiguous plants (thin or heterozygous-looking) contribute nothing.
/// Same [`LocusSlipFit`] shape as [`attribute_locus`], so the same `refine_theta_locus` /
/// `refit_level_multiplier` M-steps consume it. Determinism: integer slip counts + a fixed
/// per-locus sample/observation order → thread-invariant.
fn attribute_clean_homs(
    locus: &CohortLocus,
    period: usize,
    params: &ParamSet,
    level_per_group: &[StutterLevel],
    min_depth: u64,
    hom_frac: f64,
) -> LocusSlipFit {
    let mut fit = LocusSlipFit::default();
    for k_present in 0..locus.samples.len() {
        let seq_counts = &locus.samples[k_present].seq_counts;
        // Length histogram (few distinct lengths → linear merge, no map alloc).
        let mut len_counts: SmallVec<[(usize, u64); 8]> = SmallVec::new();
        let mut total: u64 = 0;
        for (obs, count) in seq_counts {
            let c = *count as u64;
            total += c;
            match len_counts.iter_mut().find(|(l, _)| *l == obs.len()) {
                Some((_, acc)) => *acc += c,
                None => len_counts.push((obs.len(), c)),
            }
        }
        if total < min_depth {
            continue;
        }
        let (modal_len, modal_count) = len_counts
            .iter()
            .copied()
            .max_by_key(|(_, c)| *c)
            .unwrap_or((0, 0));
        if (modal_count as f64) < hom_frac * total as f64 {
            continue; // not clean-hom → ambiguous, excluded
        }
        let modal_units = (modal_len / period) as i32;
        let (_eps, level) = sample_chemistry(locus, k_present, params, level_per_group);
        let parent_level = (level.baseline + level.slope * modal_units as f64).clamp(0.0, 1.0);
        for (obs, count) in seq_counts {
            let delta = (obs.len() / period) as i32 - modal_units;
            let c = *count as u64;
            if delta != 0 {
                fit.profile.add_slip(delta, c);
                fit.slipped += c;
            }
            fit.expected_slipped += parent_level * c as f64;
        }
    }
    fit
}

/// Attribute the locus's reads to their nearest called allele and accumulate the slip
/// sufficient statistics. Integer slip counts ⇒ the shape refit is order-independent;
/// the level's float `expected_slipped` is summed in fixed per-locus order, so the whole
/// per-locus call stays identical on whichever thread runs it. The soft per-read
/// responsibility split is a deferred refinement.
fn attribute_locus(
    locus: &CohortLocus,
    calls: &[SampleCall],
    period: usize,
    candidates: &CandidateSet,
    params: &ParamSet,
    level_per_group: &[StutterLevel],
) -> LocusSlipFit {
    let mut fit = LocusSlipFit::default();
    let mut scratch = HmmScratch::new();
    for (k_present, call) in calls.iter().enumerate() {
        if call.allele_indices.is_empty() {
            continue; // no-call contributes no slips
        }
        let (eps, level) = sample_chemistry(locus, k_present, params, level_per_group);
        // The called alleles as (sequence, length-in-units), parallel to the call — so a read
        // is attributed to its parent allele **by sequence**, same-length ties broken by
        // composition (spec §5.3). Slip counts stay integers → order-free reduce.
        let called: SmallVec<[(&[u8], u16); 2]> = call
            .allele_indices
            .iter()
            .zip(&call.genotype_units)
            .map(|(&idx, &units)| (candidates.alleles[idx].as_ref(), units))
            .collect();
        for (obs, count) in &locus.samples[k_present].seq_counts {
            let read_units = (obs.len() / period) as i32;
            // PANIC-FREE: the loop above skips empty-`allele_indices` calls, and `called` is built
            // parallel to `allele_indices`/`genotype_units`, so it is non-empty here;
            // `nearest_called_by_sequence` returns `None` only for an empty `called` slice.
            let (parent_pos, delta) =
                nearest_called_by_sequence(obs, read_units, &called, eps, &mut scratch)
                    .expect("a non-empty call has ≥1 allele");
            let parent = called[parent_pos].1;
            let c = *count as u64;
            if delta != 0 {
                fit.profile.add_slip(delta, c);
                fit.slipped += c;
            }
            // The group's predicted slip probability at the parent length, for the rate
            // multiplier's "expected slips" denominator.
            let parent_level = (level.baseline + level.slope * parent as f64).clamp(0.0, 1.0);
            fit.expected_slipped += parent_level * c as f64;
        }
    }
    fit
}

/// Fill each sample's [`SampleCall::allele_support`] — the deconvolved per-called-allele read
/// responsibilities `n_a = Σ_reads count · Qᵣ(obs | a) / Σ_{a'∈G} Qᵣ(obs | a')` for the called
/// genotype `G` (Q-I1 / arch §3, Mark-2 §6 amendment). Computed here, at the final E-step,
/// because this is where the read model + candidates + the settled `θ_locus` / level live; the
/// VCF stage then reads `allele_support` directly, with no read model and no length attribution.
///
/// Determinism: `Qᵣ` is a pure function of the frozen `ε` and the settled shape/level, and the
/// responsibility sums run over `seq_counts` in its fixed byte-sorted order (the Stage-1
/// contract), so `(n_A, n_B)` is byte-identical across threads (arch §4). A read that explains
/// none of the called alleles (`Σ Qᵣ = 0` — random junk `λ` owns) contributes nothing.
// The model + its scratch join the locus shaping, chemistry, and settled shape/level — all
// distinct inputs threaded explicitly so this stays a pure function of them (the byte-identity
// contract), like `compute_data_ll`. Two over the lint's bound.
#[allow(clippy::too_many_arguments)]
fn fill_allele_support<M: ReadLikelihoodModel>(
    calls: &mut [SampleCall],
    locus_model: &LocusModel,
    params: &ParamSet,
    level_per_group: &[StutterLevel],
    theta: &StutterShape,
    level_multiplier: f64,
    model: &M,
    scratch: &mut M::Scratch,
) {
    let LocusModel {
        locus,
        candidates,
        cand_units,
        cand_interruptions,
        ..
    } = *locus_model;
    for (k_present, call) in calls.iter_mut().enumerate() {
        if call.allele_indices.is_empty() {
            continue; // a no-call keeps its empty support
        }
        let (eps, level) = sample_chemistry(locus, k_present, params, level_per_group);
        // Per-called-allele level — same formula (purity included) as the data likelihood, so
        // the two never disagree — computed once per sample, not per read.
        let called_levels: SmallVec<[f64; 2]> = call
            .allele_indices
            .iter()
            .map(|&c| {
                candidate_level(
                    level,
                    cand_units[c],
                    level_multiplier,
                    &params.purity_level,
                    cand_interruptions[c],
                )
            })
            .collect();
        let mut support: SmallVec<[f64; 2]> = smallvec![0.0; call.allele_indices.len()];
        for (obs, count) in &locus.samples[k_present].seq_counts {
            // Qᵣ of this observed sequence against each called allele (a handful of distinct
            // sequences × the ploidy called alleles — cheap; arch §3).
            let qrs: SmallVec<[f64; 2]> = call
                .allele_indices
                .iter()
                .zip(&called_levels)
                .map(|(&c, &lvl)| {
                    let ctx = ReadScoringContext {
                        motif: &locus.motif,
                        shape: theta,
                        level: lvl,
                        eps,
                    };
                    model.q_r(obs, &candidates.alleles[c], &ctx, scratch)
                })
                .collect();
            if qrs.iter().sum::<f64>() <= 0.0 {
                continue; // the read explains none of the called alleles
            }
            let resp = allele_responsibilities(&qrs);
            for (s, r) in support.iter_mut().zip(&resp) {
                *s += *count as f64 * r;
            }
        }
        call.allele_support = support;
    }
}

/// The per-locus stutter-rate multiplier on the group level (I2): `observed / expected`
/// slips, shrunk toward `1` with pseudo-count `strength`. With few slips (or a thin
/// locus) it collapses to `1` (the group rate); a genuinely stuttery locus pushes it up.
///
/// HIERARCHY (no double-count): the per-group stutter level is the prior *mean*,
/// estimated upstream by the burn-in's `reduce_level` from the called genotypes (not
/// from this multiplier); `m` is the per-locus *deviation* around it, applied only at
/// genotyping time. At burn-in convergence loci that match the group have `m ≈ 1`, so
/// the group level and `m` model the rate once, hierarchically — not twice.
fn refit_level_multiplier(fit: &LocusSlipFit, strength: f64) -> f64 {
    ((fit.slipped as f64 + strength) / (fit.expected_slipped + strength)).clamp(0.0, LEVEL_MULT_MAX)
}

/// Whether two shapes agree within `tol` on every coefficient (the `θ_locus` stop).
fn shapes_close(a: &StutterShape, b: &StutterShape, tol: f64) -> bool {
    (a.up_rate - b.up_rate).abs() < tol
        && (a.down_rate - b.down_rate).abs() < tol
        && (a.decay - b.decay).abs() < tol
}

/// Phred-scaled genotype quality from the MAP posterior, capped at 99.
fn phred_gq(posterior: f64) -> u8 {
    let err = (1.0 - posterior).max(1e-10);
    let q = (-10.0 * err.log10()).round();
    q.clamp(0.0, 99.0) as u8
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::candidate_set::{CandidateCfg, assemble_candidates};
    use crate::ssr::cohort::em_init::seed_locus;
    use crate::ssr::cohort::param_estimation::{
        G0PseudocountDecay, PerBaseError, PurityLevel, SampleGroupId, StutterShape,
    };
    use crate::ssr::cohort::rung_ladder::{RungCfg, build_rungs};
    use crate::ssr::cohort::sim::{
        SimChemistry, SimCohort, SimCohortSpec, SimGenotype, SimLocus, SimSample,
    };
    use crate::ssr::types::Motif;
    use std::collections::HashMap;

    /// softmax of a log-prior row, for the marginalized-prior tests.
    fn softmax(logs: &[f64]) -> Vec<f64> {
        let m = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let exps: Vec<f64> = logs.iter().map(|&l| (l - m).exp()).collect();
        let z: f64 = exps.iter().sum();
        exps.iter().map(|&e| e / z).collect()
    }

    #[test]
    fn marginalized_prior_at_f_zero_is_the_bare_dm_prior() {
        // With F = 0 the Wright mixture collapses: every genotype's log-prior is
        // just its Dirichlet-multinomial term, so the wrapper must reproduce
        // `dirichlet_multinomial_log_priors` exactly.
        let genotypes = enumerate_diploid_genotypes(3);
        let alpha = [0.6, 0.3, 0.1];
        let (counts, coeffs) = diploid_dm_inputs(&genotypes, 3);
        let bare = crate::genetics::dirichlet_multinomial_log_priors(&counts, &coeffs, 3, &alpha);
        let got = marginalized_genotype_log_priors(&genotypes, 3, &alpha, 0.0);
        assert_eq!(got.len(), bare.len());
        for (g, b) in got.iter().zip(&bare) {
            assert!((g - b).abs() < 1e-12, "got {g}, bare {b}");
        }
    }

    #[test]
    fn marginalized_prior_reduces_to_plugin_hwe_at_high_concentration() {
        // As Σα → ∞ with α ∝ π the Dirichlet-multinomial collapses to the
        // multinomial, i.e. plug-in HWE(π). The softmaxed marginalized prior must
        // approach HWE(π) = [p0², 2·p0·p1, p1²].
        let genotypes = enumerate_diploid_genotypes(2);
        let pi = [0.7_f64, 0.3];
        let big = 1.0e6;
        let alpha = [big * pi[0], big * pi[1]];
        let probs = softmax(&marginalized_genotype_log_priors(
            &genotypes, 2, &alpha, 0.0,
        ));
        // HWE reference at F = 0 via the plug-in prior (already normalised: the
        // three genotype priors sum to 1 for a biallelic locus).
        let hwe: Vec<f64> = genotypes
            .iter()
            .map(|&g| genotype_prior(g, &pi, 0.0))
            .collect();
        for (p, h) in probs.iter().zip(&hwe) {
            assert!((p - h).abs() < 1e-4, "marginalized {p} vs HWE {h}");
        }
    }

    #[test]
    fn marginalized_prior_at_f_one_forbids_heterozygotes() {
        // Full inbreeding: every heterozygote is impossible (log-prior −∞) and
        // the homozygotes carry finite mass at the IBD marginal allele frequency.
        let genotypes = enumerate_diploid_genotypes(2);
        let alpha = [0.6, 0.4];
        let logs = marginalized_genotype_log_priors(&genotypes, 2, &alpha, 1.0);
        for (g, &lp) in genotypes.iter().zip(&logs) {
            if g.i == g.j {
                assert!(
                    lp.is_finite(),
                    "homozygote {g:?} should be finite, got {lp}"
                );
            } else {
                assert!(
                    lp == f64::NEG_INFINITY,
                    "heterozygote {g:?} must be forbidden at F=1, got {lp}",
                );
            }
        }
    }

    #[test]
    fn leave_one_out_alpha_excludes_own_counts_and_floors_at_seed() {
        let g0 = [1.0_f64, 0.5];
        // α'_s = G₀ + (cohort − own); own is excluded.
        let alpha = leave_one_out_alpha(&g0, &[3.0, 2.0], &[1.0, 1.0]);
        assert!((alpha[0] - 3.0).abs() < 1e-12); // 1.0 + (3.0 − 1.0)
        assert!((alpha[1] - 1.5).abs() < 1e-12); // 0.5 + (2.0 − 1.0)
        // Float noise making (cohort − own) slightly negative floors at G₀, never
        // below (keeps the concentration strictly positive).
        let noisy = leave_one_out_alpha(&g0, &[1.0, 1.0], &[1.0 + 1e-15, 1.0 + 1e-15]);
        assert!((noisy[0] - 1.0).abs() < 1e-12);
        assert!((noisy[1] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn marginalized_prior_at_f_one_equals_ibd_marginal_log_frequency() {
        // At F = 1 a homozygote collapses to the IBD term alone, so its value is
        // the hand-computable marginal allele log-frequency ln(α_i / Σα). This
        // pins the −log_sum_alpha normalization the boundary tests can't see (a
        // shared constant washes out under softmax at F=1).
        let genotypes = enumerate_diploid_genotypes(2);
        let alpha = [0.6_f64, 0.4];
        let sum: f64 = alpha.iter().sum();
        let logs = marginalized_genotype_log_priors(&genotypes, 2, &alpha, 1.0);
        for (g, &lp) in genotypes.iter().zip(&logs) {
            if g.i == g.j {
                let expected = (alpha[g.i] / sum).ln();
                assert!(
                    (lp - expected).abs() < 1e-12,
                    "hom {g:?}: got {lp}, want {expected}"
                );
            }
        }
    }

    #[test]
    fn marginalized_prior_mixes_dm_and_ibd_at_intermediate_f() {
        // 0 < F < 1: reconstruct every genotype's log-prior from the documented
        // formula, catching any error that vanishes at the F=0/F=1 boundaries
        // (swapped log_f/log_one_minus_f, a shifted IBD term, or a het branch that
        // wrongly picks up an IBD contribution).
        let genotypes = enumerate_diploid_genotypes(3);
        let alpha = [0.6_f64, 0.3, 0.1];
        let f = 0.25_f64;
        let (counts, coeffs) = diploid_dm_inputs(&genotypes, 3);
        let dm = crate::genetics::dirichlet_multinomial_log_priors(&counts, &coeffs, 3, &alpha);
        let sum: f64 = alpha.iter().sum();
        let got = marginalized_genotype_log_priors(&genotypes, 3, &alpha, f);
        for ((g, &dm_g), &lp) in genotypes.iter().zip(&dm).zip(&got) {
            let expected = if g.i == g.j {
                log_sum_exp(&[(1.0 - f).ln() + dm_g, f.ln() + (alpha[g.i] / sum).ln()])
            } else {
                (1.0 - f).ln() + dm_g
            };
            assert!(
                (lp - expected).abs() < 1e-12,
                "{g:?}: got {lp}, want {expected}"
            );
        }
    }

    #[test]
    fn diploid_dm_inputs_encodes_homozygote_as_double_count_and_het_as_two_singles() {
        // k=2 genotypes: (0,0), (0,1), (1,1). Isolates the flat-input builder from
        // the DM math + mixture so a shape bug localizes here, not downstream.
        let genotypes = enumerate_diploid_genotypes(2);
        let (counts, coeffs) = diploid_dm_inputs(&genotypes, 2);
        assert_eq!(counts, vec![2, 0, 1, 1, 0, 2]);
        assert_eq!(coeffs, vec![0.0, std::f64::consts::LN_2, 0.0]);
    }

    #[test]
    fn marginalized_prior_handles_single_candidate_locus() {
        // k=1 boundary: one genotype (0,0); ln(α/α) = 0 so the homozygote's IBD
        // term is 0 and the prior stays finite.
        let genotypes = enumerate_diploid_genotypes(1);
        let (counts, coeffs) = diploid_dm_inputs(&genotypes, 1);
        assert_eq!(counts, vec![2]);
        assert_eq!(coeffs, vec![0.0]);
        let logs = marginalized_genotype_log_priors(&genotypes, 1, &[0.5], 0.3);
        assert_eq!(logs.len(), 1);
        assert!(
            logs[0].is_finite(),
            "single-candidate prior must be finite, got {}",
            logs[0]
        );
    }

    #[test]
    fn marginalized_prior_returns_empty_on_empty_genotypes() {
        // Degenerate empty candidate set: no genotypes in, no priors out (no
        // index/unwrap on a first genotype, no empty-sum division).
        let logs = marginalized_genotype_log_priors(&[], 2, &[0.6, 0.4], 0.5);
        assert!(logs.is_empty());
    }

    #[test]
    fn emission_evidence_separates_monomorphic_from_polymorphic() {
        // Two candidate alleles; genotype data-likelihoods are built order-independently so the
        // test does not depend on the genotype enumeration order.
        let genotypes = enumerate_diploid_genotypes(2); // (0,0),(0,1),(1,1)
        let pi = [0.5, 0.5];
        let f = [0.0, 0.0];
        let ll = |want_i: usize, want_j: usize| -> Vec<f64> {
            genotypes
                .iter()
                .map(|g| {
                    if g.i == want_i && g.j == want_j {
                        0.0
                    } else {
                        -20.0
                    }
                })
                .collect()
        };
        // Monomorphic: both plants strongly hom-ref → the extra allele earns nothing, so the
        // polymorphic marginal cannot beat the best single-allele model (BIC delta small/negative).
        let mono = vec![ll(0, 0), ll(0, 0)];
        let e = emission_evidence(&mono, &genotypes, &pi, &f, 2);
        assert!(
            2.0 * (e.ln_marginal - e.ln_monomorphic) < 5.0,
            "monomorphic locus should not favour the polymorphic model"
        );
        // Polymorphic: the allele segregates (one plant hom-ref, one hom-alt) → the marginal
        // sits far above the best single-allele model (large positive BIC delta).
        let poly = vec![ll(0, 0), ll(1, 1)];
        let e = emission_evidence(&poly, &genotypes, &pi, &f, 2);
        assert!(
            2.0 * (e.ln_marginal - e.ln_monomorphic) > 10.0,
            "a segregating allele should favour the polymorphic model"
        );
    }

    /// Params for a single low-stutter sample group matching the clean sim chemistry.
    fn clean_params(n_samples: usize) -> ParamSet {
        let mut parent = HashMap::new();
        parent.insert(
            2u8,
            StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
        );
        let mut decay = HashMap::new();
        decay.insert(2u8, G0PseudocountDecay { p: 0.5 });
        ParamSet {
            error_per_sample_group: vec![PerBaseError(0.001)],
            stutter_shape_parent: parent,
            stutter_shape_by_cell: HashMap::new(),
            level_seed: vec![StutterLevel {
                baseline: 0.05,
                slope: 0.0,
            }],
            pseudocount_decay_per_loci_group: decay,
            group_of_sample: vec![SampleGroupId(0); n_samples],
            f0_seed: 0.0,
            purity_level: PurityLevel::none(),
        }
    }

    /// One CA locus, a cohort of homozygotes + separated hets at high depth.
    fn checkpoint_spec() -> SimCohortSpec {
        let chem = SimChemistry {
            error: PerBaseError(0.001),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.05,
                slope: 0.0,
            },
        };
        let genos = [
            SimGenotype::homozygous(8, 2),
            SimGenotype::homozygous(8, 2),
            SimGenotype::diploid(6, 10),
            SimGenotype::homozygous(10, 2),
            SimGenotype::diploid(6, 10),
        ];
        SimCohortSpec {
            seed: 2026,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 40,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 8,
            }],
            groups: vec![chem],
            samples: genos
                .iter()
                .enumerate()
                .map(|(i, g)| SimSample {
                    name: format!("S{i}"),
                    group: SampleGroupId(0),
                    genotypes: vec![Some(g.clone())],
                })
                .collect(),
            depth: 120,
        }
    }

    fn call_locus(cohort: &SimCohort) -> LocusCall {
        let locus = cohort.merger().next().unwrap().expect("one locus");
        let (_, locus) = locus;
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        let params = clean_params(locus.present.len());
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        run_locus_em(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            2,
            &EmCfg::dev_default(),
        )
    }

    /// M3 (review): `null_from_homs` is emission-only. It substitutes a clean-hom stutter null
    /// into the emission evidence but must leave the genotyping `calls`/`pi` untouched (they are
    /// computed from the genotyping `data_ll`, which the toggle never recomputes).
    #[test]
    fn null_from_homs_leaves_genotypes_unchanged() {
        let cohort = crate::ssr::cohort::sim::simulate(&checkpoint_spec());
        let locus = cohort.merger().next().unwrap().expect("one locus").1;
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        let params = clean_params(locus.present.len());
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        let f_per_present = vec![0.0; locus.present.len()];

        let run = |cfg: &EmCfg| {
            run_locus_em_with(
                &locus,
                &rungs,
                &candidates,
                &params,
                &seed,
                2,
                cfg,
                &f_per_present,
                &params.level_seed,
                &HipstrModel,
            )
        };
        let off = EmCfg::dev_default();
        let on = EmCfg {
            null_from_homs: true,
            ..EmCfg::dev_default()
        };
        let (call_off, _) = run(&off);
        let (call_on, _) = run(&on);

        // Genotypes and frequencies are byte-identical; only the emission evidence may move.
        assert_eq!(
            call_off.calls, call_on.calls,
            "null_from_homs must not change GTs"
        );
        assert_eq!(call_off.pi, call_on.pi, "null_from_homs must not change π");
    }

    #[test]
    fn checkpoint_1_called_genotypes_match_truth_at_high_depth() {
        let cohort = crate::ssr::cohort::sim::simulate(&checkpoint_spec());
        let call = call_locus(&cohort);
        assert_eq!(call.admit, Admission::Pass);

        let truth = [[8, 8], [8, 8], [6, 10], [10, 10], [6, 10]];
        for (k, expected) in truth.iter().enumerate() {
            let mut got = call.calls[k].genotype_units.clone();
            got.sort_unstable();
            let mut want = expected.to_vec();
            want.sort_unstable();
            assert_eq!(
                got, want,
                "sample {k}: called {:?} but truth is {:?} (GQ {})",
                call.calls[k].genotype_units, expected, call.calls[k].gq
            );
        }
    }

    #[test]
    fn em_calls_correct_genotypes_under_moderate_stutter() {
        // The same genotypes as checkpoint 1, but a markedly higher stutter level —
        // an EM that only works at near-zero stutter would fail here.
        let baseline = 0.15;
        let chem = SimChemistry {
            error: PerBaseError(0.001),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline,
                slope: 0.0,
            },
        };
        let genos = [
            SimGenotype::homozygous(8, 2),
            SimGenotype::homozygous(8, 2),
            SimGenotype::diploid(6, 10),
            SimGenotype::homozygous(10, 2),
            SimGenotype::diploid(6, 10),
        ];
        let spec = SimCohortSpec {
            seed: 77,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 40,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 8,
            }],
            groups: vec![chem],
            samples: genos
                .iter()
                .enumerate()
                .map(|(i, g)| SimSample {
                    name: format!("S{i}"),
                    group: SampleGroupId(0),
                    genotypes: vec![Some(g.clone())],
                })
                .collect(),
            depth: 200,
        };
        let cohort = crate::ssr::cohort::sim::simulate(&spec);
        let (_, locus) = cohort.merger().next().unwrap().expect("one locus");
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        let mut params = clean_params(locus.present.len());
        params.level_seed[0].baseline = baseline;
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        let call = run_locus_em(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            2,
            &EmCfg::dev_default(),
        );

        let truth = [[8, 8], [8, 8], [6, 10], [10, 10], [6, 10]];
        for (k, expected) in truth.iter().enumerate() {
            let mut got = call.calls[k].genotype_units.clone();
            got.sort_unstable();
            let mut want = expected.to_vec();
            want.sort_unstable();
            assert_eq!(got, want, "sample {k} under moderate stutter");
        }
    }

    #[test]
    fn high_depth_clean_calls_are_confident() {
        let cohort = crate::ssr::cohort::sim::simulate(&checkpoint_spec());
        let call = call_locus(&cohort);
        // Every clean high-depth call should be confident.
        assert!(
            call.calls.iter().all(|c| c.gq >= 20),
            "GQs: {:?}",
            call.calls.iter().map(|c| c.gq).collect::<Vec<_>>()
        );
    }

    #[test]
    fn pi_concentrates_on_the_true_alleles() {
        let cohort = crate::ssr::cohort::sim::simulate(&checkpoint_spec());
        let call = call_locus(&cohort);
        let total: f64 = call.pi.iter().sum();
        assert!((total - 1.0).abs() < 1e-9);
        // The true alleles are 6, 8, 10 units; their combined π should dominate.
        assert!(call.pi.iter().all(|&p| p >= 0.0));
    }

    // ── I1: the per-locus θ_locus shape refit ──

    fn ca(units: u16) -> Box<[u8]> {
        "CA".repeat(units as usize).into_bytes().into_boxed_slice()
    }

    #[test]
    fn attribute_locus_slips_bins_reads_by_signed_distance() {
        use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 0,
                end: 16,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGG".as_slice()),
            ca(8),
        );
        locus.push(
            0,
            SampleEvidence {
                seq_counts: vec![(ca(6), 3), (ca(7), 5), (ca(8), 100), (ca(9), 4)],
                qc: SsrQc::default(),
            },
        );
        // A homozygous-8 call: 8 is faithful (no slip), 7 is Δ=−1, 9 is Δ=+1, 6 is Δ=−2.
        let calls = vec![SampleCall {
            allele_indices: vec![0, 0],
            genotype_units: vec![8, 8],
            posterior: 0.99,
            gq: 40,
            allele_support: SmallVec::new(),
        }];
        let params = clean_params(1);
        let candidates = CandidateSet {
            alleles: vec![ca(8)],
            ref_idx: 0,
            admit: Admission::Pass,
        };
        let fit = attribute_locus(&locus, &calls, 2, &candidates, &params, &params.level_seed);
        assert_eq!(fit.profile.down[0], 5, "Δ=−1 count");
        assert_eq!(fit.profile.up[0], 4, "Δ=+1 count");
        assert_eq!(fit.profile.down[1], 3, "Δ=−2 count");
        assert_eq!(fit.profile.up[1], 0);
        // 5 + 4 + 3 = 12 slipped reads; the faithful 8s are not slips.
        assert_eq!(fit.slipped, 12);
    }

    #[test]
    fn attribute_locus_skips_no_calls() {
        use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 0,
                end: 16,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGG".as_slice()),
            ca(8),
        );
        locus.push(
            0,
            SampleEvidence {
                seq_counts: vec![(ca(7), 5)],
                qc: SsrQc::default(),
            },
        );
        let params = clean_params(1);
        let candidates = CandidateSet {
            alleles: vec![ca(8)],
            ref_idx: 0,
            admit: Admission::Pass,
        };
        let fit = attribute_locus(
            &locus,
            &[SampleCall::no_call()],
            2,
            &candidates,
            &params,
            &params.level_seed,
        );
        assert_eq!(fit.profile.down.iter().sum::<u64>(), 0);
        assert_eq!(fit.profile.up.iter().sum::<u64>(), 0);
        assert_eq!(fit.slipped, 0);
        assert_eq!(fit.expected_slipped, 0.0);
    }

    // ── `attribute_clean_homs` (the `null_from_homs` clean-hom stutter donor pool) ──

    /// Build a one-sample CA locus (ref 8 units) with the given per-length read counts.
    fn clean_hom_locus(seq_counts: Vec<(Box<[u8]>, u32)>) -> CohortLocus {
        use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 0,
                end: 16,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGG".as_slice()),
            ca(8),
        );
        locus.push(
            0,
            SampleEvidence {
                seq_counts,
                qc: SsrQc::default(),
            },
        );
        locus
    }

    #[test]
    fn attribute_clean_homs_bins_off_modal_reads_as_slips() {
        // modal length = 8 units (100/112 ≈ 0.89 ≥ 0.75); 7 → Δ=−1, 9 → Δ=+1, 6 → Δ=−2.
        let locus = clean_hom_locus(vec![(ca(6), 3), (ca(7), 5), (ca(8), 100), (ca(9), 4)]);
        let params = clean_params(1);
        let fit = attribute_clean_homs(&locus, 2, &params, &params.level_seed, 4, 0.75);
        assert_eq!(fit.profile.down[0], 5, "Δ=−1");
        assert_eq!(fit.profile.up[0], 4, "Δ=+1");
        assert_eq!(fit.profile.down[1], 3, "Δ=−2");
        assert_eq!(fit.slipped, 12);
        // `expected_slipped` sums the group level over *every* read (112) in a fixed order —
        // pinning it guards the float sum against a reordering that would break thread-identity.
        // level = baseline 0.05 + slope 0.0 · 8 = 0.05 ⇒ 0.05 · 112 = 5.6.
        assert!((fit.expected_slipped - 5.6).abs() < 1e-9);
    }

    #[test]
    fn attribute_clean_homs_excludes_a_plant_below_hom_frac() {
        // 50/50 het: modal fraction 0.5 < 0.75 → excluded; the true second allele's reads must
        // NOT be counted as stutter (the circularity `attribute_clean_homs` exists to avoid).
        let locus = clean_hom_locus(vec![(ca(8), 6), (ca(10), 6)]);
        let params = clean_params(1);
        let fit = attribute_clean_homs(&locus, 2, &params, &params.level_seed, 4, 0.75);
        assert_eq!(fit.slipped, 0);
        assert_eq!(fit.profile.up.iter().sum::<u64>(), 0);
        assert_eq!(fit.profile.down.iter().sum::<u64>(), 0);
        assert_eq!(fit.expected_slipped, 0.0);
    }

    #[test]
    fn attribute_clean_homs_excludes_a_plant_below_min_depth() {
        // total 3 < min_depth 4 → excluded (guards the `total < min_depth` gate boundary).
        let locus = clean_hom_locus(vec![(ca(8), 2), (ca(7), 1)]);
        let params = clean_params(1);
        let fit = attribute_clean_homs(&locus, 2, &params, &params.level_seed, 4, 0.75);
        assert_eq!(fit.slipped, 0);
        assert_eq!(fit.expected_slipped, 0.0);
    }

    #[test]
    fn attribute_clean_homs_empty_pool_round_trips_to_the_prior() {
        // Every donor filtered (too thin) ⇒ empty fit ⇒ `refine_theta_locus` returns the prior
        // and `refit_level_multiplier` returns 1.0 — i.e. the toggle is a safe no-op with no
        // clean-hom evidence.
        let locus = clean_hom_locus(vec![(ca(8), 1)]);
        let params = clean_params(1);
        let fit = attribute_clean_homs(&locus, 2, &params, &params.level_seed, 4, 0.75);
        assert_eq!(fit.slipped, 0);
        assert_eq!(fit.expected_slipped, 0.0);
        let prior = StutterShape {
            up_rate: 1.0,
            down_rate: 2.0,
            decay: 0.1,
        };
        assert_eq!(refine_theta_locus(&fit.profile, &prior, 50.0), prior);
        assert_eq!(refit_level_multiplier(&fit, 20.0), 1.0);
    }

    #[test]
    fn attribute_clean_homs_tie_break_is_deterministic() {
        // Two lengths tied for modal (5 each); `max_by_key` returns the LAST maximal element
        // (10 units), so the 8-unit reads slip to Δ=−2. `hom_frac = 0.5` admits the tie.
        let locus = clean_hom_locus(vec![(ca(8), 5), (ca(10), 5)]);
        let params = clean_params(1);
        let fit = attribute_clean_homs(&locus, 2, &params, &params.level_seed, 4, 0.5);
        assert_eq!(
            fit.profile.down[1], 5,
            "8-unit reads slip off the last-tied modal 10"
        );
        assert_eq!(fit.slipped, 5);
    }

    #[test]
    fn attribute_clean_homs_counts_far_slip_in_slipped_but_not_profile() {
        // A read |Δ| = 12 units > MAX_SLIP (10): `add_slip` drops it from the shape profile, but
        // `slipped` still counts it — the same divergence `attribute_locus` has.
        let locus = clean_hom_locus(vec![(ca(8), 20), (ca(20), 3)]);
        let params = clean_params(1);
        let fit = attribute_clean_homs(&locus, 2, &params, &params.level_seed, 4, 0.75);
        assert_eq!(fit.slipped, 3, "far read counts toward slipped");
        assert_eq!(
            fit.profile.up.iter().sum::<u64>(),
            0,
            "but is dropped from the profile"
        );
        assert_eq!(fit.profile.down.iter().sum::<u64>(), 0);
    }

    #[test]
    #[should_panic(expected = "frozen sample group")]
    fn sample_chemistry_panics_on_a_sample_missing_its_group() {
        // A ParamSet whose group_of_sample is shorter than the present samples violates
        // decision E; sample_chemistry must fail loud, not fabricate group-0 chemistry (M2).
        use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 0,
                end: 16,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGG".as_slice()),
            ca(8),
        );
        locus.push(
            0,
            SampleEvidence {
                seq_counts: vec![(ca(8), 10)],
                qc: SsrQc::default(),
            },
        );
        let mut params = clean_params(1);
        params.group_of_sample.clear(); // break the decision-E density invariant
        let _ = sample_chemistry(&locus, 0, &params, &params.level_seed);
    }

    #[test]
    fn refit_level_multiplier_collapses_to_one_without_excess_slips() {
        // Observed slips equal the group-expected slips ⇒ multiplier 1 (the group rate).
        let matched = LocusSlipFit {
            profile: SlipProfile::default(),
            slipped: 10,
            expected_slipped: 10.0,
        };
        assert!((refit_level_multiplier(&matched, 20.0) - 1.0).abs() < 1e-9);
        // No data ⇒ also 1 (collapses to the group rate, no oscillation).
        assert!((refit_level_multiplier(&LocusSlipFit::default(), 20.0) - 1.0).abs() < 1e-9);
        // A genuinely stuttery locus (more observed than expected) pushes the rate up,
        // but shrinkage keeps a thin signal modest.
        let stuttery = LocusSlipFit {
            profile: SlipProfile::default(),
            slipped: 60,
            expected_slipped: 20.0,
        };
        let m = refit_level_multiplier(&stuttery, 20.0);
        assert!(m > 1.0 && m < 3.0, "shrunk multiplier {m}");
    }

    #[test]
    fn shapes_close_compares_every_coefficient() {
        let a = StutterShape {
            up_rate: 0.3,
            down_rate: 0.7,
            decay: 0.2,
        };
        assert!(shapes_close(&a, &a, 1e-9));
        let mut b = a;
        b.decay = 0.25; // a 0.05 change exceeds the tolerance
        assert!(!shapes_close(&a, &b, 1e-3));
    }

    #[test]
    fn candidate_level_scales_by_purity_only_for_impure_alleles() {
        let level = StutterLevel {
            baseline: 0.1,
            slope: 0.0,
        };
        let neutral = PurityLevel::none();
        let halving = PurityLevel {
            per_interruption_factor: 0.5,
        };
        // A pure allele (0 interruptions) is unchanged — even under a non-neutral factor —
        // so Phase 2 collapses to Phase 1 wherever alleles are pure.
        assert!((candidate_level(level, 6, 1.0, &neutral, 0) - 0.1).abs() < 1e-12);
        assert!((candidate_level(level, 6, 1.0, &halving, 0) - 0.1).abs() < 1e-12);
        // An impure allele's level is scaled once per interruption (compounding).
        assert!((candidate_level(level, 6, 1.0, &halving, 1) - 0.05).abs() < 1e-12);
        assert!((candidate_level(level, 6, 1.0, &halving, 2) - 0.025).abs() < 1e-12);
        // The result stays clamped to [0, 1].
        let high = StutterLevel {
            baseline: 0.9,
            slope: 0.1,
        };
        assert!(candidate_level(high, 20, 5.0, &neutral, 0) <= 1.0);
    }

    #[test]
    fn level_multiplier_recovers_calls_under_an_underestimated_group_level() {
        // True stutter level 0.22 but the group seed claims a much lower 0.05: the
        // per-locus rate multiplier must scale the level up and still call the truth.
        let chem = SimChemistry {
            error: PerBaseError(0.001),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.22,
                slope: 0.0,
            },
        };
        let genos = [
            SimGenotype::homozygous(8, 2),
            SimGenotype::homozygous(8, 2),
            SimGenotype::diploid(6, 10),
            SimGenotype::homozygous(10, 2),
            SimGenotype::diploid(6, 10),
        ];
        let spec = SimCohortSpec {
            seed: 4040,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 40,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 8,
            }],
            groups: vec![chem],
            samples: genos
                .iter()
                .enumerate()
                .map(|(i, g)| SimSample {
                    name: format!("S{i}"),
                    group: SampleGroupId(0),
                    genotypes: vec![Some(g.clone())],
                })
                .collect(),
            depth: 200,
        };
        let cohort = crate::ssr::cohort::sim::simulate(&spec);
        let (_, locus) = cohort.merger().next().unwrap().expect("one locus");
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        // clean_params keeps the (wrong-for-this-locus) low group level 0.05.
        let params = clean_params(locus.present.len());
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        let call = run_locus_em(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            2,
            &EmCfg::dev_default(),
        );

        let truth = [[8, 8], [8, 8], [6, 10], [10, 10], [6, 10]];
        for (k, expected) in truth.iter().enumerate() {
            let mut got = call.calls[k].genotype_units.clone();
            got.sort_unstable();
            let mut want = expected.to_vec();
            want.sort_unstable();
            assert_eq!(got, want, "sample {k} under an underestimated group level");
        }
    }

    #[test]
    fn refit_max_rounds_zero_reproduces_the_seed_shape_result() {
        // With the θ refit disabled the EM is the pre-I1 (frozen-seed-shape) genotyper;
        // it must still call the clean truth — pinning that the refit is purely additive.
        let cohort = crate::ssr::cohort::sim::simulate(&checkpoint_spec());
        let (_, locus) = cohort.merger().next().unwrap().expect("one locus");
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        let params = clean_params(locus.present.len());
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        let no_refit = EmCfg {
            refit_max_rounds: 0,
            ..EmCfg::dev_default()
        };
        let call = run_locus_em(&locus, &rungs, &candidates, &params, &seed, 2, &no_refit);

        let truth = [[8, 8], [8, 8], [6, 10], [10, 10], [6, 10]];
        for (k, expected) in truth.iter().enumerate() {
            let mut got = call.calls[k].genotype_units.clone();
            got.sort_unstable();
            let mut want = expected.to_vec();
            want.sort_unstable();
            assert_eq!(got, want, "sample {k} with θ refit disabled");
        }
    }

    #[test]
    fn marginalized_prior_recovers_the_clean_truth_end_to_end() {
        // The marginalize+LOO Dirichlet-multinomial prior (EmCfg toggle on) must
        // genotype the clean high-depth checkpoint cohort exactly as the plug-in
        // path does — end-to-end proof that the new π-EM loop + final calls are
        // wired correctly and converge to the right genotypes.
        let cohort = crate::ssr::cohort::sim::simulate(&checkpoint_spec());
        let (_, locus) = cohort.merger().next().unwrap().expect("one locus");
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        let params = clean_params(locus.present.len());
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        let marginalized = EmCfg {
            marginalized_prior: true,
            ..EmCfg::dev_default()
        };
        let call = run_locus_em(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            2,
            &marginalized,
        );

        let truth = [[8, 8], [8, 8], [6, 10], [10, 10], [6, 10]];
        for (k, expected) in truth.iter().enumerate() {
            let mut got = call.calls[k].genotype_units.clone();
            got.sort_unstable();
            let mut want = expected.to_vec();
            want.sort_unstable();
            assert_eq!(got, want, "sample {k} under the marginalized DM prior");
        }
    }

    #[test]
    fn marginalized_and_plugin_agree_on_clean_high_depth() {
        // On clean high-depth data the marginalized DM prior and the plug-in HWE
        // prior should call the *same* genotypes — a drift-detector independent of
        // the hardcoded truth array (a marginalized regression that still happened
        // to match the truth would be caught here).
        let cohort = crate::ssr::cohort::sim::simulate(&checkpoint_spec());
        let (_, locus) = cohort.merger().next().unwrap().expect("one locus");
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        let params = clean_params(locus.present.len());
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        let plugin = run_locus_em(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            2,
            &EmCfg::dev_default(),
        );
        let marginalized_cfg = EmCfg {
            marginalized_prior: true,
            ..EmCfg::dev_default()
        };
        let marginalized = run_locus_em(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            2,
            &marginalized_cfg,
        );
        for (k, (p, m)) in plugin.calls.iter().zip(&marginalized.calls).enumerate() {
            assert_eq!(
                p.genotype_units, m.genotype_units,
                "sample {k}: marginalized and plug-in disagree on clean data",
            );
        }
    }

    #[test]
    fn run_pi_em_marginalized_single_sample_stays_finite_and_collapses_to_g0() {
        // n_samples == 1 boundary: `total = 2` (finite π), and the leave-one-out
        // subtraction removes the sample's whole contribution, so `α'_s == G₀` —
        // `own_counts[0]` must equal the cohort `expected`.
        let genotypes = enumerate_diploid_genotypes(2); // (0,0), (0,1), (1,1)
        let g0 = [1.0, 1.0];
        let data_ll = vec![vec![0.0, -20.0, -40.0]]; // strongly hom-ref (0,0)
        let cfg = EmCfg {
            marginalized_prior: true,
            ..EmCfg::dev_default()
        };
        let fit = run_pi_em_marginalized(&data_ll, &genotypes, &g0, &[0.0], 2, &cfg);
        assert!(
            fit.pi.iter().all(|p| p.is_finite()),
            "π must be finite: {:?}",
            fit.pi
        );
        let sum: f64 = fit.expected.iter().sum();
        assert!(
            (sum - 2.0).abs() < 1e-9,
            "expected copies must sum to 2·n_samples, got {sum}"
        );
        for (own, exp) in fit.own_counts[0].iter().zip(&fit.expected) {
            assert!(
                (own - exp).abs() < 1e-9,
                "single sample: own must equal cohort expected"
            );
        }
    }

    #[test]
    fn run_pi_em_marginalized_leave_one_out_moves_an_ambiguous_sample() {
        // Three samples strongly hom-ref (0,0) and one ambiguous between (0,1) and
        // (1,1). Running only the flat iteration (max_iters = 1, no LOO) vs the full
        // loop must change the cohort `expected`, proving the leave-one-out
        // iterations actually feed the three hom-ref samples' counts into the
        // ambiguous sample's prior. Deterministic via hand-built likelihoods.
        let genotypes = enumerate_diploid_genotypes(2);
        let g0 = [1.0, 1.0];
        let data_ll = vec![
            vec![0.0, -30.0, -60.0],
            vec![0.0, -30.0, -60.0],
            vec![0.0, -30.0, -60.0],
            vec![-5.0, 0.0, -0.5], // (0,1) vs (1,1) close; (0,0) unlikely
        ];
        let f = [0.0; 4];
        let flat = EmCfg {
            marginalized_prior: true,
            max_iters: 1,
            ..EmCfg::dev_default()
        };
        let full = EmCfg {
            marginalized_prior: true,
            ..EmCfg::dev_default()
        };
        let fit_flat = run_pi_em_marginalized(&data_ll, &genotypes, &g0, &f, 2, &flat);
        let fit_full = run_pi_em_marginalized(&data_ll, &genotypes, &g0, &f, 2, &full);
        let moved = fit_flat
            .expected
            .iter()
            .zip(&fit_full.expected)
            .any(|(a, b)| (a - b).abs() > 1e-6);
        assert!(
            moved,
            "leave-one-out iterations did not move the cohort expected counts"
        );
    }

    #[test]
    fn marginalized_inbreeding_is_threaded_and_raises_posterior_homozygosity() {
        // A sample whose reads slightly favour the het (0,1) over the hom (0,0):
        // with F = 0 the het wins, but a strong F engages the IBD branch through
        // the loop AND the final calls, raising the posterior homozygosity. Proves
        // `f_per_present` is threaded end-to-end (not silently dropped to 0).
        let genotypes = enumerate_diploid_genotypes(2);
        let g0 = [1.0, 1.0];
        let cand_units = [5u16, 6u16];
        let data_ll = vec![vec![-1.0, 0.0, -3.0], vec![-1.0, 0.0, -3.0]];
        let posterior_hom_at = |f: f64| {
            let cfg = EmCfg {
                marginalized_prior: true,
                ..EmCfg::dev_default()
            };
            let fpp = [f; 2];
            let fit = run_pi_em_marginalized(&data_ll, &genotypes, &g0, &fpp, 2, &cfg);
            let (_, ph) =
                final_calls_marginalized(&data_ll, &genotypes, &g0, &fit, &fpp, 2, &cand_units);
            ph
        };
        let ph_outbred = posterior_hom_at(0.0);
        let ph_inbred = posterior_hom_at(0.9);
        assert!(
            ph_inbred[0] > ph_outbred[0] + 1e-6,
            "F>0 did not raise posterior homozygosity (F not threaded): {} vs {}",
            ph_inbred[0],
            ph_outbred[0],
        );
    }

    #[test]
    fn capped_refit_returns_calls_consistent_with_the_final_round() {
        // A capped refit (refit_max_rounds: 1) that may not fully converge must still
        // return the calls recomputed at the end of its last round — never stale pre-loop
        // calls. At high depth one round already recovers the truth and agrees with the
        // converged (refit_max_rounds: 3) genotypes (review Mi15).
        let cohort = crate::ssr::cohort::sim::simulate(&checkpoint_spec());
        let (_, locus) = cohort.merger().next().unwrap().expect("one locus");
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        let params = clean_params(locus.present.len());
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        let capped = EmCfg {
            refit_max_rounds: 1,
            ..EmCfg::dev_default()
        };
        let one = run_locus_em(&locus, &rungs, &candidates, &params, &seed, 2, &capped);
        let three = run_locus_em(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            2,
            &EmCfg::dev_default(),
        );

        let truth = [[8, 8], [8, 8], [6, 10], [10, 10], [6, 10]];
        for (k, expected) in truth.iter().enumerate() {
            let mut got = one.calls[k].genotype_units.clone();
            got.sort_unstable();
            let mut want = expected.to_vec();
            want.sort_unstable();
            assert_eq!(got, want, "capped refit sample {k} should still call truth");
            let mut got3 = three.calls[k].genotype_units.clone();
            got3.sort_unstable();
            assert_eq!(got, got3, "capped vs converged disagree at sample {k}");
        }
    }

    #[test]
    fn theta_locus_refit_recovers_calls_under_a_shape_mismatched_seed() {
        // Heavy-tailed stutter (decay 0.5) but a seed shape claiming low decay (0.1):
        // the per-locus θ refit must adapt and still call the truth.
        let chem = SimChemistry {
            error: PerBaseError(0.001),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.5,
            },
            level: StutterLevel {
                baseline: 0.15,
                slope: 0.0,
            },
        };
        let genos = [
            SimGenotype::homozygous(8, 2),
            SimGenotype::homozygous(8, 2),
            SimGenotype::diploid(6, 10),
            SimGenotype::homozygous(10, 2),
            SimGenotype::diploid(6, 10),
        ];
        let spec = SimCohortSpec {
            seed: 909,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 40,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 8,
            }],
            groups: vec![chem],
            samples: genos
                .iter()
                .enumerate()
                .map(|(i, g)| SimSample {
                    name: format!("S{i}"),
                    group: SampleGroupId(0),
                    genotypes: vec![Some(g.clone())],
                })
                .collect(),
            depth: 200,
        };
        let cohort = crate::ssr::cohort::sim::simulate(&spec);
        let (_, locus) = cohort.merger().next().unwrap().expect("one locus");
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        // clean_params carries the (wrong-for-this-locus) low-decay 0.1 seed shape.
        let mut params = clean_params(locus.present.len());
        params.level_seed[0].baseline = 0.15;
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        let call = run_locus_em(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            2,
            &EmCfg::dev_default(),
        );

        let truth = [[8, 8], [8, 8], [6, 10], [10, 10], [6, 10]];
        for (k, expected) in truth.iter().enumerate() {
            let mut got = call.calls[k].genotype_units.clone();
            got.sort_unstable();
            let mut want = expected.to_vec();
            want.sort_unstable();
            assert_eq!(
                got, want,
                "sample {k} under shape-mismatched seed + θ refit"
            );
        }
    }

    #[test]
    fn no_call_locus_reports_no_call_samples() {
        // A thin cohort → LowDepth → every sample no-called.
        let chem = SimChemistry {
            error: PerBaseError(0.001),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.05,
                slope: 0.0,
            },
        };
        let spec = SimCohortSpec {
            seed: 1,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 40,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 8,
            }],
            groups: vec![chem],
            samples: vec![SimSample {
                name: "S0".into(),
                group: SampleGroupId(0),
                genotypes: vec![Some(SimGenotype::homozygous(8, 2))],
            }],
            depth: 2, // below the LowDepth floor
        };
        let cohort = crate::ssr::cohort::sim::simulate(&spec);
        let call = call_locus(&cohort);
        assert_eq!(call.admit, Admission::LowDepth);
        assert!(call.calls.iter().all(|c| c.allele_indices.is_empty()));
    }

    #[test]
    fn run_locus_em_fills_allele_support_from_the_final_e_step() {
        // The Q-I1 data wire: after the EM, every called sample carries a populated
        // `allele_support` (the deconvolved per-allele responsibilities the balance term reads,
        // arch §3), parallel to its genotype and summing to the sample's read depth — the
        // invariant that a genotype's reads are fully apportioned across its called alleles.
        use crate::ssr::cohort::candidate_set::assemble_candidates;
        use crate::ssr::cohort::em_init::seed_locus;
        use crate::ssr::cohort::rung_ladder::{RungCfg, build_rungs};
        use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};

        let pure = ca(6); // CACACACACACA (12 bp)
        let interrupted = {
            let mut s = pure.to_vec();
            s[5] = b'T';
            s.into_boxed_slice()
        };
        let sample = |obs: Vec<(Box<[u8]>, u32)>| {
            let mut seq_counts = obs;
            seq_counts.sort_by(|a, b| a.0.cmp(&b.0));
            SampleEvidence {
                seq_counts,
                qc: SsrQc::default(),
            }
        };
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 40,
                end: 52,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGG".as_slice()),
            pure.clone(),
        );
        let mut evs = Vec::new();
        for _ in 0..3 {
            evs.push(sample(vec![(pure.clone(), 60)]));
        }
        for _ in 0..3 {
            evs.push(sample(vec![(interrupted.clone(), 60)]));
        }
        evs.push(sample(vec![(pure.clone(), 30), (interrupted.clone(), 30)]));
        for (i, ev) in evs.into_iter().enumerate() {
            locus.push(i as u32, ev);
        }

        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        assert_eq!(candidates.alleles.len(), 2, "pure + interrupted candidates");
        let params = clean_params(locus.present.len());
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        let call = run_locus_em(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            2,
            &EmCfg::dev_default(),
        );

        let mut any_call = false;
        for (k, c) in call.calls.iter().enumerate() {
            if c.allele_indices.is_empty() {
                continue;
            }
            any_call = true;
            assert_eq!(
                c.allele_support.len(),
                c.allele_indices.len(),
                "allele_support is parallel to the called genotype (sample {k})"
            );
            let depth: f64 = locus.samples[k]
                .seq_counts
                .iter()
                .map(|(_, n)| f64::from(*n))
                .sum();
            let support_total: f64 = c.allele_support.iter().sum();
            assert!(
                (support_total - depth).abs() < 1e-6,
                "responsibilities apportion all {depth} reads (sample {k}, got {support_total})"
            );
        }
        assert!(
            any_call,
            "the cohort produced at least one non-no-call sample"
        );
    }

    #[test]
    fn fill_allele_support_skips_a_read_matching_no_called_allele() {
        // A read wildly off both called alleles (Σ Qᵣ = 0 — junk the λ term owns) must
        // contribute nothing to allele_support; the matching reads still apportion fully
        // (review §8). The junk read is too sparse to be nominated as a candidate.
        use crate::ssr::cohort::candidate_set::assemble_candidates;
        use crate::ssr::cohort::em_init::seed_locus;
        use crate::ssr::cohort::rung_ladder::{RungCfg, build_rungs};
        use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};

        let pure = ca(6);
        let interrupted = {
            let mut s = pure.to_vec();
            s[5] = b'T';
            s.into_boxed_slice()
        };
        let junk = ca(30); // 60 bp — 24 units off the 6-unit alleles → Qᵣ underflows to 0
        let sample = |obs: Vec<(Box<[u8]>, u32)>| {
            let mut seq_counts = obs;
            seq_counts.sort_by(|a, b| a.0.cmp(&b.0));
            SampleEvidence {
                seq_counts,
                qc: SsrQc::default(),
            }
        };
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 40,
                end: 52,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGG".as_slice()),
            pure.clone(),
        );
        let mut evs = Vec::new();
        for _ in 0..3 {
            evs.push(sample(vec![(pure.clone(), 60)]));
        }
        for _ in 0..3 {
            evs.push(sample(vec![(interrupted.clone(), 60)]));
        }
        // Carrier: hom-pure (60 reads) + 2 junk reads (too few to be a clear peak → not a
        // candidate) at a wildly different length.
        evs.push(sample(vec![(pure.clone(), 60), (junk.clone(), 2)]));
        let carrier = evs.len() - 1;
        for (i, ev) in evs.into_iter().enumerate() {
            locus.push(i as u32, ev);
        }

        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        assert_eq!(
            candidates.alleles.len(),
            2,
            "the junk read is not promoted to a candidate"
        );
        let params = clean_params(locus.present.len());
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        let call = run_locus_em(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            2,
            &EmCfg::dev_default(),
        );

        let c = &call.calls[carrier];
        assert!(!c.allele_indices.is_empty(), "carrier is called");
        let support_total: f64 = c.allele_support.iter().sum();
        // The 2 junk reads must NOT be apportioned: support sums to the 60 matching reads only.
        assert!(
            (support_total - 60.0).abs() < 1e-6,
            "junk read must be skipped (support {support_total}, expected 60)"
        );
    }

    #[test]
    fn low_depth_same_length_het_carrier_is_recovered_by_the_em() {
        // The P1.5 exit gate (spec §5.3/§9): the seed CANNOT express a same-length het (it sees
        // one length peak with one representative → a hom seed), so a genuine same-length het
        // carrier is recovered ONLY if the composition-aware E-step moves posterior weight onto
        // the het. This fixture puts such a carrier at realistic LOW depth (2 pure + 2
        // interrupted reads) alongside common hom carriers of both alleles, and asserts the
        // carrier is called HET — the recovery the seed can't do. (F = 0 here; the low-depth /
        // high-F_IS undercall tail is a measured benchmark concern, not this unit invariant.)
        use crate::ssr::cohort::candidate_set::assemble_candidates;
        use crate::ssr::cohort::em_init::seed_locus;
        use crate::ssr::cohort::rung_ladder::{RungCfg, build_rungs};
        use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};

        let pure = ca(6);
        let interrupted = {
            let mut s = pure.to_vec();
            s[5] = b'T';
            s.into_boxed_slice()
        };
        let sample = |obs: Vec<(Box<[u8]>, u32)>| {
            let mut seq_counts = obs;
            seq_counts.sort_by(|a, b| a.0.cmp(&b.0));
            SampleEvidence {
                seq_counts,
                qc: SsrQc::default(),
            }
        };
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 40,
                end: 52,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGG".as_slice()),
            pure.clone(),
        );
        // 4 hom-pure + 4 hom-interrupted (both alleles common → admitted, π not tiny), then the
        // low-depth same-length HET carrier as the last sample.
        let mut evs = Vec::new();
        for _ in 0..4 {
            evs.push(sample(vec![(pure.clone(), 60)]));
        }
        for _ in 0..4 {
            evs.push(sample(vec![(interrupted.clone(), 60)]));
        }
        evs.push(sample(vec![(pure.clone(), 2), (interrupted.clone(), 2)])); // 4 reads, het
        let carrier = evs.len() - 1;
        for (i, ev) in evs.into_iter().enumerate() {
            locus.push(i as u32, ev);
        }

        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        assert_eq!(
            candidates.alleles.len(),
            2,
            "two same-length candidates admitted"
        );
        let params = clean_params(locus.present.len());
        let seed = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        let call = run_locus_em(
            &locus,
            &rungs,
            &candidates,
            &params,
            &seed,
            2,
            &EmCfg::dev_default(),
        );

        // The locus is variable, and the low-depth carrier is called HET — two DISTINCT
        // candidate indices, i.e. one pure and one interrupted allele.
        let carrier_call = &call.calls[carrier];
        assert_eq!(carrier_call.allele_indices.len(), 2, "carrier is called");
        assert_ne!(
            carrier_call.allele_indices[0], carrier_call.allele_indices[1],
            "the low-depth same-length het is recovered as a het, not collapsed to a hom \
             (indices {:?})",
            carrier_call.allele_indices
        );
        // Its deconvolved support is balanced (the 2/2 split), so the allele-balance term does
        // not wrongly penalise this genuine same-length het.
        let support = &carrier_call.allele_support;
        let total: f64 = support.iter().sum();
        assert!(
            support[0].min(support[1]) / total > 0.3,
            "a balanced same-length het is not flagged imbalanced (support {support:?})"
        );
    }
}
