//! The prior-side outer loop: per-individual inbreeding `F_i` + per-group stutter
//! level, wrapping the per-locus EM over all loci (arch `ssr_call_genotyping.md` §5,
//! spec §4.4, implementation plan §4 E1).
//!
//! Each round runs the per-locus EM with the current `F_i` (per sample) and level
//! (per group), then reduces:
//! - **`F_i`** from the posterior homozygosity over *variable* loci vs. the HWE
//!   expectation `Σ π_i²` (the excess-homozygosity estimator), shrunk toward the
//!   cohort mean and clamped to `[0, F_CEILING]`;
//! - **the per-group level** by hard-attributing each read to its called allele and
//!   refitting `baseline + slope·length`.
//!
//! Both reduces go through `FixedPointAccum` / integer counts so they are
//! order-independent (byte-identical across thread counts — the determinism the
//! parallel F1 step proves). The level refit here is the hard-label form; the soft
//! per-allele responsibility reduce is the deferred refinement.

use rayon::prelude::*;

use crate::ssr::cohort::candidate_set::{CandidateCfg, CandidateSet, assemble_candidates};
use crate::ssr::cohort::em::{EmCfg, LocusCall, run_locus_em_with};
use crate::ssr::cohort::em_init::{LocusSeed, seed_locus};
use crate::ssr::cohort::param_estimation::SampleStutterStats;
use crate::ssr::cohort::param_estimation::{FixedPointAccum, ParamSet, StutterLevel};
use crate::ssr::cohort::prepass::{add_bin, fit_level};
use crate::ssr::cohort::read_model::HipstrModel;
use crate::ssr::cohort::rung_ladder::{RungCfg, Rungs, build_rungs};
use crate::ssr::cohort::types::CohortLocus;

/// The hard ceiling on `F` (spec §4.4: an apparent `F_IS` above this is almost
/// always a data artifact, not real inbreeding).
const F_CEILING: f64 = 0.99;

/// Controls for the prior-side outer loop (pinned in F2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct OuterCfg {
    /// Maximum outer rounds.
    pub(crate) max_rounds: usize,
    /// Convergence tolerance on the largest `F_i` change between rounds.
    pub(crate) f_tol: f64,
    /// Convergence tolerance on the largest level coefficient change between rounds.
    pub(crate) level_tol: f64,
    /// Pseudo-count weight shrinking each `F_i` toward the cohort mean.
    pub(crate) f_shrink: f64,
    /// User cap on `F` (≤ [`F_CEILING`]).
    pub(crate) f_cap: f64,
}

impl OuterCfg {
    pub(crate) fn dev_default() -> Self {
        Self {
            max_rounds: 8,
            f_tol: 1e-3,
            level_tol: 1e-3,
            f_shrink: 5.0,
            f_cap: F_CEILING,
        }
    }
}

/// The outer-loop result: the final per-locus calls + the converged `F_i` / level.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CohortCalls {
    /// Final per-locus genotype calls (parallel to the input loci).
    pub(crate) per_locus: Vec<LocusCall>,
    /// Converged per-(global)-sample inbreeding `F`.
    pub(crate) f_per_sample: Vec<f64>,
    /// Converged per-group stutter level.
    pub(crate) level_per_group: Vec<StutterLevel>,
}

/// Per-locus inputs that do not change across outer rounds (they depend on the data
/// and the frozen params, not on `F` / level).
struct Prepared {
    rungs: Rungs,
    candidates: CandidateSet,
    seed: LocusSeed,
}

/// The HWE homozygosity expectation `Σ π_i²` for a locus.
fn expected_homozygosity(pi: &[f64]) -> f64 {
    pi.iter().map(|p| p * p).sum()
}

/// Reduce `F_i` from the per-locus calls (excess homozygosity, shrunk + clamped).
fn reduce_f(
    loci: &[CohortLocus],
    calls: &[LocusCall],
    n_global: usize,
    cfg: &OuterCfg,
    f0: f64,
) -> Vec<f64> {
    // Per sample: accumulate the per-locus `F` estimate over variable loci.
    let mut acc = vec![FixedPointAccum::new(); n_global];
    let mut counts = vec![0u32; n_global];
    for (locus, call) in loci.iter().zip(calls) {
        let exp_hom = expected_homozygosity(&call.pi);
        if exp_hom >= 1.0 - 1e-9 {
            continue; // monomorphic locus carries no inbreeding signal
        }
        for (k, &global) in locus.present.iter().enumerate() {
            let p_hom = call.posterior_hom[k];
            let f_locus = (p_hom - exp_hom) / (1.0 - exp_hom);
            acc[global as usize].add(f_locus);
            counts[global as usize] += 1;
        }
    }

    let raw: Vec<f64> = (0..n_global)
        .map(|s| {
            if counts[s] > 0 {
                acc[s].value() / counts[s] as f64
            } else {
                f0
            }
        })
        .collect();

    // Shrink each toward the cohort mean (over samples with data), then clamp.
    let observed: Vec<f64> = (0..n_global)
        .filter(|&s| counts[s] > 0)
        .map(|s| raw[s])
        .collect();
    let cohort_mean = if observed.is_empty() {
        f0
    } else {
        observed.iter().sum::<f64>() / observed.len() as f64
    };
    let ceiling = cfg.f_cap.min(F_CEILING);
    (0..n_global)
        .map(|s| {
            let n = counts[s] as f64;
            let shrunk = (n * raw[s] + cfg.f_shrink * cohort_mean) / (n + cfg.f_shrink);
            shrunk.clamp(0.0, ceiling)
        })
        .collect()
}

/// Reduce the per-group level by hard-attributing each read to its called allele and
/// refitting the level line per group.
fn reduce_level(
    loci: &[CohortLocus],
    prepared: &[Prepared],
    calls: &[LocusCall],
    params: &ParamSet,
    n_groups: usize,
) -> Vec<StutterLevel> {
    let mut group_stats: Vec<SampleStutterStats> = vec![SampleStutterStats::default(); n_groups];
    for ((locus, prep), call) in loci.iter().zip(prepared).zip(calls) {
        let period = prep.rungs.period();
        for (k, &global) in locus.present.iter().enumerate() {
            let sample_call = &call.calls[k];
            if sample_call.genotype_units.is_empty() {
                continue; // no-call
            }
            let raw_group = params
                .group_of_sample
                .get(global as usize)
                .map(|g| g.0 as usize)
                .unwrap_or(0);
            debug_assert!(
                raw_group < n_groups,
                "sample group {raw_group} out of range for {n_groups} level groups \
                 (group_of_sample / level_seed size mismatch)"
            );
            let group = raw_group.min(n_groups - 1);
            for (obs, count) in &locus.samples[k].seq_counts {
                let read_units = (obs.len() / period) as u16;
                // Attribute to the nearest called allele length.
                let parent = *sample_call
                    .genotype_units
                    .iter()
                    .min_by_key(|&&u| (u as i32 - read_units as i32).abs())
                    .unwrap();
                let (faithful, slipped) = if read_units == parent {
                    (*count as u64, 0)
                } else {
                    (0, *count as u64)
                };
                add_bin(&mut group_stats[group].by_length, parent, faithful, slipped);
            }
        }
    }
    group_stats.iter().map(fit_level).collect()
}

/// Run the prior-side outer loop over a set of loci.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_cohort_em(
    loci: &[CohortLocus],
    params: &ParamSet,
    level_seed: Vec<StutterLevel>,
    ploidy: u8,
    em_cfg: &EmCfg,
    rung_cfg: &RungCfg,
    cand_cfg: &CandidateCfg,
    outer_cfg: &OuterCfg,
) -> CohortCalls {
    let prepared: Vec<Prepared> = loci
        .iter()
        .map(|locus| {
            let rungs = build_rungs(locus, rung_cfg);
            let candidates = assemble_candidates(locus, &rungs, ploidy, cand_cfg);
            let seed = seed_locus(
                locus,
                &rungs,
                &candidates,
                params,
                ploidy,
                rung_cfg.prominence,
            );
            Prepared {
                rungs,
                candidates,
                seed,
            }
        })
        .collect();

    let n_global = params.group_of_sample.len();
    let n_groups = level_seed.len().max(1);
    let mut f_per_sample = vec![params.f0_seed; n_global];
    let mut level_per_group = level_seed;
    let mut per_locus: Vec<LocusCall> = Vec::new();

    for _ in 0..outer_cfg.max_rounds {
        // Per-locus EM is independent and pure; rayon's indexed `collect` preserves
        // locus order, so the result is byte-identical across thread counts (F1).
        per_locus = loci
            .par_iter()
            .zip(prepared.par_iter())
            .map(|(locus, prep)| {
                let f_present: Vec<f64> = locus
                    .present
                    .iter()
                    .map(|&g| f_per_sample[g as usize])
                    .collect();
                run_locus_em_with(
                    locus,
                    &prep.rungs,
                    &prep.candidates,
                    params,
                    &prep.seed,
                    ploidy,
                    em_cfg,
                    &f_present,
                    &level_per_group,
                    &HipstrModel,
                )
            })
            .collect();

        let new_f = reduce_f(loci, &per_locus, n_global, outer_cfg, params.f0_seed);
        let new_level = reduce_level(loci, &prepared, &per_locus, params, n_groups);

        let f_delta = new_f
            .iter()
            .zip(&f_per_sample)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        let level_delta = new_level
            .iter()
            .zip(&level_per_group)
            .map(|(a, b)| {
                (a.baseline - b.baseline)
                    .abs()
                    .max((a.slope - b.slope).abs())
            })
            .fold(0.0, f64::max);
        f_per_sample = new_f;
        level_per_group = new_level;
        if f_delta < outer_cfg.f_tol && level_delta < outer_cfg.level_tol {
            break;
        }
    }

    CohortCalls {
        per_locus,
        f_per_sample,
        level_per_group,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::param_estimation::{
        G0FitCfg, G0PseudocountDecay, PerBaseError, SampleGroupId, StutterShape,
    };
    use crate::ssr::cohort::sim::{
        SimChemistry, SimCohortSpec, SimGenotype, SimLocus, SimSample, simulate,
    };
    use crate::ssr::types::Motif;
    use std::collections::HashMap;

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
            error_per_sample_group: vec![PerBaseError(0.002)],
            stutter_shape_parent: parent,
            stutter_shape_by_cell: HashMap::new(),
            level_seed: vec![StutterLevel {
                baseline: 0.06,
                slope: 0.0,
            }],
            pseudocount_decay_per_loci_group: decay,
            group_of_sample: vec![SampleGroupId(0); n_samples],
            f0_seed: 0.0,
        }
    }

    /// A cohort split into inbred (homozygous everywhere) and outbred (het
    /// everywhere) samples across many biallelic 8/10 loci.
    fn inbreeding_spec() -> (SimCohortSpec, usize) {
        let chem = SimChemistry {
            error: PerBaseError(0.002),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.06,
                slope: 0.0,
            },
        };
        let n_loci = 16;
        let n_inbred = 6;
        let n_outbred = 6;
        let mut samples = Vec::new();
        // Inbred samples: homozygous, alternating 8 and 10 across loci.
        for i in 0..n_inbred {
            let genotypes = (0..n_loci)
                .map(|l| {
                    let allele = if (i + l) % 2 == 0 { 8 } else { 10 };
                    Some(SimGenotype::homozygous(allele, 2))
                })
                .collect();
            samples.push(SimSample {
                name: format!("inbred{i}"),
                group: SampleGroupId(0),
                genotypes,
            });
        }
        // Outbred samples: heterozygous 8/10 at every locus.
        for i in 0..n_outbred {
            samples.push(SimSample {
                name: format!("outbred{i}"),
                group: SampleGroupId(0),
                genotypes: vec![Some(SimGenotype::diploid(8, 10)); n_loci],
            });
        }
        let loci = (0..n_loci)
            .map(|l| SimLocus {
                chrom: "chr1".into(),
                start: 40 + 40 * l as u32,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 8,
            })
            .collect();
        let spec = SimCohortSpec {
            seed: 2718,
            loci,
            groups: vec![chem],
            samples,
            depth: 120,
        };
        (spec, n_inbred)
    }

    fn collect_loci(spec: &SimCohortSpec) -> Vec<CohortLocus> {
        simulate(spec)
            .merger()
            .map(|item| item.expect("locus"))
            .map(|(_, locus)| locus)
            .collect()
    }

    #[test]
    fn f_loop_separates_inbred_from_outbred_samples() {
        let (spec, n_inbred) = inbreeding_spec();
        let loci = collect_loci(&spec);
        let params = clean_params(spec.samples.len());
        let result = run_cohort_em(
            &loci,
            &params,
            params.level_seed.clone(),
            2,
            &EmCfg::dev_default(),
            &RungCfg::dev_default(),
            &CandidateCfg::dev_default(),
            &OuterCfg::dev_default(),
        );

        let inbred_mean: f64 =
            result.f_per_sample[..n_inbred].iter().sum::<f64>() / n_inbred as f64;
        let outbred_mean: f64 = result.f_per_sample[n_inbred..].iter().sum::<f64>()
            / (result.f_per_sample.len() - n_inbred) as f64;
        assert!(
            inbred_mean > outbred_mean + 0.3,
            "F should separate inbred ({inbred_mean:.3}) from outbred ({outbred_mean:.3})"
        );
        assert!(
            result
                .f_per_sample
                .iter()
                .all(|&f| (0.0..=F_CEILING).contains(&f))
        );
    }

    #[test]
    fn full_pipeline_calls_and_emits_a_variant_vcf_line() {
        use crate::ssr::cohort::candidate_set::assemble_candidates;
        use crate::ssr::cohort::param_estimation::{G0PseudocountDecay, PerBaseError};
        use crate::ssr::cohort::prepass::{estimate, run_prepass_stats};
        use crate::ssr::cohort::rung_ladder::build_rungs;
        use crate::ssr::cohort::sample_groups::{ClusterCfg, group_samples};
        use crate::ssr::cohort::vcf_out::{
            FpControlCfg, apply_fp_control, format_vcf_record, is_variable, site_qual,
        };
        use std::collections::HashMap;

        // A clean variant cohort at one locus: 6 homozygotes (8) + 6 separated hets (6/10).
        let chem = SimChemistry {
            error: PerBaseError(0.002),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.06,
                slope: 0.0,
            },
        };
        let mut samples: Vec<SimSample> = (0..6)
            .map(|i| SimSample {
                name: format!("hom{i}"),
                group: SampleGroupId(0),
                genotypes: vec![Some(SimGenotype::homozygous(8, 2))],
            })
            .collect();
        samples.extend((0..6).map(|i| SimSample {
            name: format!("het{i}"),
            group: SampleGroupId(0),
            genotypes: vec![Some(SimGenotype::diploid(6, 10))],
        }));
        let spec = SimCohortSpec {
            seed: 555,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 40,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 8,
            }],
            groups: vec![chem],
            samples,
            depth: 150,
        };
        let loci = collect_loci(&spec);

        // Pre-pass → estimate → cluster.
        let stats = run_prepass_stats(&loci, 2, &RungCfg::dev_default());
        let est = estimate(&stats, &G0FitCfg::dev_default());
        let grouped = group_samples(&stats, &est, &ClusterCfg::dev_default());

        // Assemble a ParamSet from the pre-pass output (the real interface).
        let n = spec.samples.len();
        let mut decay = HashMap::new();
        decay.insert(2u8, G0PseudocountDecay { p: 0.5 });
        let group_of_sample = (0..n as u32)
            .map(|s| {
                grouped
                    .group_of_sample
                    .get(&s)
                    .copied()
                    .unwrap_or(SampleGroupId(0))
            })
            .collect();
        let params = ParamSet {
            error_per_sample_group: grouped
                .eps_per_group
                .iter()
                .map(|&e| PerBaseError(e))
                .collect(),
            stutter_shape_parent: est.shape_by_period.clone(),
            stutter_shape_by_cell: grouped.shape_by_group_period.clone(),
            level_seed: grouped.level_per_group.clone(),
            pseudocount_decay_per_loci_group: decay,
            group_of_sample,
            f0_seed: 0.0,
        };

        // Genotype with the outer loop, then FP-control + emit.
        let result = run_cohort_em(
            &loci,
            &params,
            grouped.level_per_group.clone(),
            2,
            &EmCfg::dev_default(),
            &RungCfg::dev_default(),
            &CandidateCfg::dev_default(),
            &OuterCfg::dev_default(),
        );

        let locus = &loci[0];
        let rungs = build_rungs(locus, &RungCfg::dev_default());
        let cands = assemble_candidates(locus, &rungs, 2, &CandidateCfg::dev_default());
        let fp = FpControlCfg::dev_default();
        let mut call = result.per_locus[0].clone();
        apply_fp_control(&mut call, &fp);

        assert!(is_variable(&call, &cands), "the site is polymorphic");
        let qual = site_qual(&call, &cands, &fp);
        let line = format_vcf_record("chr1", locus, &cands, &call, qual, locus.present.len());
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols[6], "PASS", "line: {line}");
        assert_ne!(cols[4], ".", "should carry an ALT allele");

        // The het samples (indices 6..12) call genotype 6/10.
        for k in 0..call.calls.len() {
            let global = locus.present[k] as usize;
            if global >= 6 {
                let mut units = call.calls[k].genotype_units.clone();
                units.sort_unstable();
                assert_eq!(units, vec![6, 10], "het sample {global} miscalled");
            }
        }
    }

    #[test]
    fn f_loop_is_deterministic() {
        let (spec, _) = inbreeding_spec();
        let loci = collect_loci(&spec);
        let params = clean_params(spec.samples.len());
        let run = || {
            run_cohort_em(
                &loci,
                &params,
                params.level_seed.clone(),
                2,
                &EmCfg::dev_default(),
                &RungCfg::dev_default(),
                &CandidateCfg::dev_default(),
                &OuterCfg::dev_default(),
            )
        };
        assert_eq!(run().f_per_sample, run().f_per_sample);
    }

    #[test]
    fn cohort_em_is_byte_identical_across_thread_counts() {
        let (spec, _) = inbreeding_spec();
        let loci = collect_loci(&spec);
        let params = clean_params(spec.samples.len());
        let run = || {
            run_cohort_em(
                &loci,
                &params,
                params.level_seed.clone(),
                2,
                &EmCfg::dev_default(),
                &RungCfg::dev_default(),
                &CandidateCfg::dev_default(),
                &OuterCfg::dev_default(),
            )
        };
        let pool = |n: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
        };
        let single = pool(1).install(run);
        let multi = pool(4).install(run);
        assert_eq!(
            single, multi,
            "the genotyping outer loop must be byte-identical across thread counts (F1)"
        );
    }

    #[test]
    fn level_refit_recovers_the_injected_level() {
        // A homozygote-only cohort: the level refit should land near the injected 0.06.
        let chem = SimChemistry {
            error: PerBaseError(0.002),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.06,
                slope: 0.0,
            },
        };
        let mut samples = Vec::new();
        for length in 7u16..=11 {
            for r in 0..4 {
                let idx = samples.len();
                samples.push(SimSample {
                    name: format!("S{idx}_{length}_{r}"),
                    group: SampleGroupId(0),
                    genotypes: vec![Some(SimGenotype::homozygous(length, 2))],
                });
            }
        }
        let spec = SimCohortSpec {
            seed: 11,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 60,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 9,
            }],
            groups: vec![chem],
            samples,
            depth: 200,
        };
        let loci = collect_loci(&spec);
        let params = clean_params(spec.samples.len());
        // Start the level seed deliberately wrong; the refit should correct it.
        let bad_seed = vec![StutterLevel {
            baseline: 0.25,
            slope: 0.0,
        }];
        let result = run_cohort_em(
            &loci,
            &params,
            bad_seed,
            2,
            &EmCfg::dev_default(),
            &RungCfg::dev_default(),
            &CandidateCfg::dev_default(),
            &OuterCfg::dev_default(),
        );
        let baseline = result.level_per_group[0].baseline;
        assert!(
            (baseline - 0.06).abs() < 0.04,
            "level refit landed at {baseline}, expected ≈ 0.06"
        );
    }
}
