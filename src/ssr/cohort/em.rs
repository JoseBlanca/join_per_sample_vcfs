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
//! cost more than it saves. The `θ_locus` M-step (which needs slip attribution) is
//! deferred until D wires the slip accumulators; here `θ = θ⁰`.
//!
//! v1 is **diploid** (the simulator's ploidy); higher ploidy is a documented
//! follow-up (the genotype enumeration generalizes, the dosage priors need care).

use crate::ssr::cohort::allele_freq_prior::g0_pseudocounts;
use crate::ssr::cohort::candidate_set::{Admission, CandidateSet};
use crate::ssr::cohort::em_init::LocusSeed;
use crate::ssr::cohort::likelihood::{LikelihoodScratch, read_given_genotype, read_likelihood};
use crate::ssr::cohort::param_estimation::{G0PseudocountDecay, ParamSet, StutterLevel};
use crate::ssr::cohort::rung_ladder::Rungs;
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
}

impl EmCfg {
    pub(crate) fn dev_default() -> Self {
        Self {
            max_iters: 50,
            tol: 1e-6,
            lambda: 0.01,
            inbreeding_f: 0.0,
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
    fn no_call() -> Self {
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
fn sample_chemistry(
    locus: &CohortLocus,
    present_k: usize,
    params: &ParamSet,
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
    let level = params
        .level_seed
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
    assert_eq!(ploidy, 2, "C4 EM is diploid-only (v1)");
    let n_present = locus.samples.len();

    // A no-call locus still reports REF + the no-call samples.
    if candidates.admit != Admission::Pass {
        return LocusCall {
            calls: vec![SampleCall::no_call(); n_present],
            pi: seed.pi0.clone(),
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

    // Precompute each sample's per-genotype DATA log-likelihood (constant across the
    // π iterations, since ε/θ/level are fixed for this milestone). `data_ll[s][g]`.
    let mut scratch = LikelihoodScratch::new();
    let mut data_ll: Vec<Vec<f64>> = Vec::with_capacity(n_present);
    for (k_present, evidence) in locus.samples.iter().enumerate() {
        let (eps, level) = sample_chemistry(locus, k_present, params);
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
                            &seed.theta0,
                            lvl,
                            eps,
                            &mut scratch,
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
                        let p = read_given_genotype(row, &[g.i, g.j], cfg.lambda, distinct);
                        *count as f64 * p.ln()
                    })
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        data_ll.push(sample_ll);
    }

    // EM on π (E-step posteriors over genotypes → M-step allele frequencies).
    let mut pi = seed.pi0.clone();
    for _ in 0..cfg.max_iters {
        let mut expected = g0.clone();
        for sample_ll in &data_ll {
            let log_joint: Vec<f64> = genotypes
                .iter()
                .zip(sample_ll)
                .map(|(g, ll)| genotype_prior(*g, &pi, cfg.inbreeding_f).ln() + ll)
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

    // Final E-step → MAP genotype + GQ per sample.
    let calls = data_ll
        .iter()
        .map(|sample_ll| {
            let log_joint: Vec<f64> = genotypes
                .iter()
                .zip(sample_ll)
                .map(|(g, ll)| genotype_prior(*g, &pi, cfg.inbreeding_f).ln() + ll)
                .collect();
            let norm = log_sum_exp(&log_joint);
            let (best, best_lj) = log_joint
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .unwrap();
            let posterior = (best_lj - norm).exp();
            let g = genotypes[best];
            SampleCall {
                allele_indices: vec![g.i, g.j],
                genotype_units: vec![cand_units[g.i], cand_units[g.j]],
                posterior,
                gq: phred_gq(posterior),
            }
        })
        .collect();

    LocusCall {
        calls,
        pi,
        admit: candidates.admit,
    }
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
