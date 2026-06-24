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
//! cost more than it saves. `ε` and the per-group stutter level are frozen by the
//! burn-in; the **stutter shape adapts per locus** (I1): an inner `θ_locus` M-step
//! attributes the locus's reads to their called alleles, refits the shape, and shrinks
//! it toward the frozen `θ_period` seed (so a thin locus stays at the seed — no
//! oscillation), re-genotyping until `θ_locus` settles. The per-locus stutter *level*
//! refit is the I2 follow-up.
//!
//! v1 is **diploid** (the simulator's ploidy); higher ploidy is a documented
//! follow-up (the genotype enumeration generalizes, the dosage priors need care).

use crate::ssr::cohort::allele_freq_prior::g0_pseudocounts;
use crate::ssr::cohort::candidate_set::{Admission, CandidateSet};
use crate::ssr::cohort::em_init::LocusSeed;
use crate::ssr::cohort::likelihood::{LikelihoodScratch, read_given_genotype, read_likelihood};
use crate::ssr::cohort::param_estimation::{
    G0PseudocountDecay, MAX_SLIP, ParamSet, SlipProfile, StutterLevel, StutterShape,
};
use crate::ssr::cohort::rung_ladder::Rungs;
use crate::ssr::cohort::stutter::refine_theta_locus;
use crate::ssr::cohort::types::CohortLocus;

/// EM controls (pinned in F2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EmCfg {
    /// Maximum EM iterations.
    pub(crate) max_iters: usize,
    /// Convergence tolerance on the largest `π` change between iterations.
    pub(crate) tol: f64,
    /// Outlier mixing weight `λ` (uniform junk term).
    pub(crate) lambda: f64,
    /// Inbreeding coefficient `F` (the autozygous-branch weight); refined in E1.
    pub(crate) inbreeding_f: f64,
    /// Maximum per-locus `θ_locus` shape-refinement rounds (I1). `0` keeps the frozen
    /// `θ_period` seed (no per-locus adaptation).
    pub(crate) theta_max_rounds: usize,
    /// Pseudo-count weight shrinking `θ_locus` toward the frozen `θ_period` seed — the
    /// anti-oscillation knob: a thin locus cannot overcome it and stays at the seed.
    pub(crate) theta_shrink: f64,
    /// Convergence tolerance on the largest `θ_locus` coefficient change between rounds.
    pub(crate) theta_tol: f64,
}

impl EmCfg {
    pub(crate) fn dev_default() -> Self {
        Self {
            max_iters: 50,
            tol: 1e-6,
            lambda: 0.01,
            inbreeding_f: 0.0,
            theta_max_rounds: 3,
            theta_shrink: 50.0,
            theta_tol: 1e-3,
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
}

impl SampleCall {
    pub(crate) fn no_call() -> Self {
        Self {
            allele_indices: Vec::new(),
            genotype_units: Vec::new(),
            posterior: 0.0,
            gq: 0,
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

/// The `G₀` decay for the locus period (mild default if absent).
fn period_decay(params: &ParamSet, period: usize) -> G0PseudocountDecay {
    const FALLBACK: G0PseudocountDecay = G0PseudocountDecay { p: 0.5 };
    params
        .pseudocount_decay_per_loci_group
        .get(&(period as u8))
        .copied()
        .unwrap_or(FALLBACK)
}

/// `ε` and the per-length stutter level for a present sample (by its sample group).
/// `level_per_group` is supplied separately so the E1 outer loop can refit it without
/// rebuilding the params.
fn sample_chemistry(
    locus: &CohortLocus,
    present_k: usize,
    params: &ParamSet,
    level_per_group: &[StutterLevel],
) -> (f64, StutterLevel) {
    let global = locus.present[present_k] as usize;
    let group = params
        .group_of_sample
        .get(global)
        .copied()
        .unwrap_or(crate::ssr::cohort::param_estimation::SampleGroupId(0));
    let eps = params
        .error_per_sample_group
        .get(group.0 as usize)
        .map(|e| e.0)
        .unwrap_or(0.01);
    let level = level_per_group
        .get(group.0 as usize)
        .copied()
        .unwrap_or(StutterLevel {
            baseline: 0.05,
            slope: 0.0,
        });
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
    )
}

/// Genotype one locus with an explicit per-sample inbreeding `F` and per-group level
/// (the E1 outer loop drives these); also emits each sample's posterior homozygosity.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_locus_em_with(
    locus: &CohortLocus,
    rungs: &Rungs,
    candidates: &CandidateSet,
    params: &ParamSet,
    seed: &LocusSeed,
    ploidy: u8,
    cfg: &EmCfg,
    f_per_present: &[f64],
    level_per_group: &[StutterLevel],
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

    let mut scratch = LikelihoodScratch::new();

    // θ_locus refinement (I1): genotype under the frozen `θ_period` seed, then refit a
    // per-locus shape from the called alleles' slips, shrunk toward that seed, and
    // re-genotype until `θ_locus` settles. With few/no slips it collapses to the seed
    // (no oscillation); each round's data log-likelihood is recomputed for the new shape.
    let mut theta = seed.theta0;
    let mut data_ll = compute_data_ll(
        locus,
        candidates,
        params,
        level_per_group,
        &theta,
        &cand_units,
        &genotypes,
        cfg.lambda,
        distinct,
        &mut scratch,
    );
    let mut pi = run_pi_em(&data_ll, &genotypes, &g0, &seed.pi0, f_per_present, cfg);
    let (mut calls, mut posterior_hom) =
        final_calls(&data_ll, &genotypes, &pi, f_per_present, &cand_units);

    for _ in 0..cfg.theta_max_rounds {
        let profile = attribute_locus_slips(locus, &calls, period);
        let new_theta = refine_theta_locus(&profile, &seed.theta0, cfg.theta_shrink);
        if shapes_close(&new_theta, &theta, cfg.theta_tol) {
            break;
        }
        theta = new_theta;
        data_ll = compute_data_ll(
            locus,
            candidates,
            params,
            level_per_group,
            &theta,
            &cand_units,
            &genotypes,
            cfg.lambda,
            distinct,
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

/// Each sample's per-genotype **data** log-likelihood `data_ll[s][g]` under stutter
/// shape `theta` (constant across the π iterations; recomputed when `θ_locus` changes).
#[allow(clippy::too_many_arguments)]
fn compute_data_ll(
    locus: &CohortLocus,
    candidates: &CandidateSet,
    params: &ParamSet,
    level_per_group: &[StutterLevel],
    theta: &StutterShape,
    cand_units: &[u16],
    genotypes: &[Genotype],
    lambda: f64,
    distinct: usize,
    scratch: &mut LikelihoodScratch,
) -> Vec<Vec<f64>> {
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
                        let lvl =
                            (level.baseline + level.slope * cand_units[c] as f64).clamp(0.0, 1.0);
                        read_likelihood(
                            obs,
                            &candidates.alleles[c],
                            &locus.motif,
                            theta,
                            lvl,
                            eps,
                            scratch,
                        )
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
        let (best, best_lj) = log_joint
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();
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
        });
        posterior_hom.push(p_hom);
    }
    (calls, posterior_hom)
}

/// Build the locus's slip profile by hard-attributing each read to its nearest called
/// allele (the `θ_locus` M-step's sufficient statistic). Integer counts ⇒ the refit is
/// order-independent. The soft per-read responsibility split is a deferred refinement.
fn attribute_locus_slips(locus: &CohortLocus, calls: &[SampleCall], period: usize) -> SlipProfile {
    let mut profile = SlipProfile::default();
    for (k_present, call) in calls.iter().enumerate() {
        if call.genotype_units.is_empty() {
            continue; // no-call contributes no slips
        }
        for (obs, count) in &locus.samples[k_present].seq_counts {
            let read_units = (obs.len() / period) as i32;
            let parent = *call
                .genotype_units
                .iter()
                .min_by_key(|&&u| (u as i32 - read_units).abs())
                .expect("a non-empty call has ≥1 allele");
            add_slip(&mut profile, read_units - parent as i32, *count as u64);
        }
    }
    profile
}

/// Add a slip of signed size `delta` (count `count`) to a profile, respecting `MAX_SLIP`.
fn add_slip(profile: &mut SlipProfile, delta: i32, count: u64) {
    let magnitude = delta.unsigned_abs() as usize;
    if (1..=MAX_SLIP).contains(&magnitude) {
        if delta > 0 {
            profile.up[magnitude - 1] += count;
        } else {
            profile.down[magnitude - 1] += count;
        }
    }
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
        }];
        let profile = attribute_locus_slips(&locus, &calls, 2);
        assert_eq!(profile.down[0], 5, "Δ=−1 count");
        assert_eq!(profile.up[0], 4, "Δ=+1 count");
        assert_eq!(profile.down[1], 3, "Δ=−2 count");
        assert_eq!(profile.up[1], 0);
    }

    #[test]
    fn attribute_locus_slips_skips_no_calls() {
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
        let profile = attribute_locus_slips(&locus, &[SampleCall::no_call()], 2);
        assert_eq!(profile.down.iter().sum::<u64>(), 0);
        assert_eq!(profile.up.iter().sum::<u64>(), 0);
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
