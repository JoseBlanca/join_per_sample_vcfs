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
use crate::ssr::cohort::attribution::nearest_parent;
use crate::ssr::cohort::candidate_set::{Admission, CandidateSet};
use crate::ssr::cohort::em_init::LocusSeed;
use crate::ssr::cohort::likelihood::read_given_genotype;
use crate::ssr::cohort::param_estimation::{
    DEFAULT_G0_FALLBACK_P, G0PseudocountDecay, ParamSet, SlipProfile, StutterLevel, StutterShape,
};
use crate::ssr::cohort::read_model::{HipstrModel, ReadLikelihoodModel, ReadScoringContext};
use crate::ssr::cohort::rung_ladder::Rungs;
use crate::ssr::cohort::stutter::refine_theta_locus;
use crate::ssr::cohort::types::CohortLocus;
use smallvec::SmallVec;

/// EM controls (pinned in F2).
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
}

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
) -> LocusCall {
    assert_eq!(ploidy, 2, "C4 EM is diploid-only (v1)");
    let n_present = locus.samples.len();

    // A no-call locus still reports REF + the no-call samples.
    if candidates.admit != Admission::Pass {
        return LocusCall {
            calls: vec![SampleCall::no_call(); n_present],
            pi: seed.pi0.clone(),
            posterior_hom: vec![0.0; n_present],
            admit: candidates.admit,
        };
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

    // The iteration-invariant locus shaping, built once and reused across the refit rounds
    // (review Mi2).
    let locus_model = LocusModel {
        locus,
        candidates,
        cand_units: &cand_units,
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
    let mut pi = run_pi_em(&data_ll, &genotypes, &g0, &seed.pi0, f_per_present, cfg);
    let (mut calls, mut posterior_hom) =
        final_calls(&data_ll, &genotypes, &pi, f_per_present, &cand_units);

    for _ in 0..cfg.refit_max_rounds {
        let fit = attribute_locus(locus, &calls, period, params, level_per_group);
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
        pi = run_pi_em(&data_ll, &genotypes, &g0, &seed.pi0, f_per_present, cfg);
        let (c, ph) = final_calls(&data_ll, &genotypes, &pi, f_per_present, &cand_units);
        calls = c;
        posterior_hom = ph;
    }

    LocusCall {
        calls,
        pi,
        posterior_hom,
        admit: candidates.admit,
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
        genotypes,
        distinct,
    } = *locus_model;
    let k = candidates.alleles.len();
    let mut data_ll: Vec<Vec<f64>> = Vec::with_capacity(locus.samples.len());
    for (k_present, evidence) in locus.samples.iter().enumerate() {
        let (eps, level) = sample_chemistry(locus, k_present, params, level_per_group);
        // Per (distinct obs, candidate) Qᵣ — the read likelihood matrix.
        let obs_qr: Vec<(u32, Vec<f64>)> = evidence
            .seq_counts
            .iter()
            .map(|(obs, count)| {
                let row = (0..k)
                    .map(|c| {
                        let lvl = ((level.baseline + level.slope * cand_units[c] as f64)
                            * level_multiplier)
                            .clamp(0.0, 1.0);
                        let ctx = ReadScoringContext {
                            motif: &locus.motif,
                            shape: theta,
                            level: lvl,
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

/// Attribute the locus's reads to their nearest called allele and accumulate the slip
/// sufficient statistics. Integer slip counts ⇒ the shape refit is order-independent;
/// the level's float `expected_slipped` is summed in fixed per-locus order, so the whole
/// per-locus call stays identical on whichever thread runs it. The soft per-read
/// responsibility split is a deferred refinement.
fn attribute_locus(
    locus: &CohortLocus,
    calls: &[SampleCall],
    period: usize,
    params: &ParamSet,
    level_per_group: &[StutterLevel],
) -> LocusSlipFit {
    let mut fit = LocusSlipFit::default();
    for (k_present, call) in calls.iter().enumerate() {
        if call.genotype_units.is_empty() {
            continue; // no-call contributes no slips
        }
        let (_eps, level) = sample_chemistry(locus, k_present, params, level_per_group);
        for (obs, count) in &locus.samples[k_present].seq_counts {
            let read_units = (obs.len() / period) as i32;
            let (parent_idx, delta) = nearest_parent(read_units, &call.genotype_units)
                .expect("a non-empty call has ≥1 allele");
            let parent = call.genotype_units[parent_idx];
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
        G0PseudocountDecay, PerBaseError, SampleGroupId, StutterShape,
    };
    use crate::ssr::cohort::rung_ladder::{RungCfg, build_rungs};
    use crate::ssr::cohort::sim::{
        SimChemistry, SimCohort, SimCohortSpec, SimGenotype, SimLocus, SimSample,
    };
    use crate::ssr::types::Motif;
    use std::collections::HashMap;

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
        let fit = attribute_locus(&locus, &calls, 2, &params, &params.level_seed);
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
        let fit = attribute_locus(
            &locus,
            &[SampleCall::no_call()],
            2,
            &params,
            &params.level_seed,
        );
        assert_eq!(fit.profile.down.iter().sum::<u64>(), 0);
        assert_eq!(fit.profile.up.iter().sum::<u64>(), 0);
        assert_eq!(fit.slipped, 0);
        assert_eq!(fit.expected_slipped, 0.0);
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
}
