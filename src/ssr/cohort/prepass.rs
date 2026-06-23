//! The parameter pre-pass: estimate `ε`, the stutter shape, and the stutter level
//! from confident genotypes (arch `ssr_call_parameters.md` §2–§4, spec §4.4,
//! implementation plan §4 D1/D2).
//!
//! D1 — from each **confident genotype** (the B1 heuristic gate: homozygotes ∪
//! well-separated hets), label every read against its nearest parent allele and
//! accumulate two order-independent integer accumulators: a per-period
//! [`SlipProfile`] (the cohort `θ_period` parent) and a per-sample
//! [`SampleStutterStats`] (the level line + `ε` + the clustering input).
//!
//! D2 — **measure**: turn those sufficient statistics into estimates — `ε` from the
//! within-tract base mismatch fraction, the per-period stutter **shape** from the
//! slip profile (direction split + geometric-decay MLE), and the per-sample stutter
//! **level** line from the faithful-vs-slipped counts by length.
//!
//! Scope notes (the refinements this milestone defers): the gate here is B1's
//! **heuristic** resolution (param-free), so the estimator is a single pass — the
//! model-based 1-vs-2-peak BIC gate, the soft full-cohort EM responsibility reduce,
//! the `ℓ_pen`-plateau burn-in iteration with multi-start, and the
//! `FixedPointAccum` float reduce all matter once the gate co-evolves with the
//! params / data is messy, and arrive with the parallel machinery (F1). On clean
//! data the hard-label single pass recovers the injected parameters (checkpoint 2).
//! Het inner-valley reads are hard-assigned to the nearer allele here (a documented
//! approximation; the soft split is the refinement).

use std::collections::HashMap;

use rayon::prelude::*;

use crate::ssr::cohort::param_estimation::{
    MAX_SLIP, SampleStutterStats, SlipProfile, StutterLevel, StutterShape,
};
use crate::ssr::cohort::rung_ladder::{
    Resolution, ResolvedGenotype, RungCfg, Rungs, build_rungs, resolve_confident_genotype,
};
use crate::ssr::cohort::types::CohortLocus;

/// The accumulated confident-genotype sufficient statistics.
#[derive(Debug, Default, Clone, PartialEq)]
pub(crate) struct PrepassStats {
    /// Period → the cohort slip profile (the `θ_period` parent).
    pub(crate) slip_by_period: HashMap<u8, SlipProfile>,
    /// `(global sample, period)` → that sample's slip profile, so D3 can
    /// re-aggregate the slip profile **per sample group** once clustering fixes the
    /// groups (the per-`(group, period)` shape).
    pub(crate) slip_by_sample_period: HashMap<(u32, u8), SlipProfile>,
    /// Global sample index → its stutter stats (level line + `ε` + depth).
    pub(crate) per_sample: HashMap<u32, SampleStutterStats>,
}

impl PrepassStats {
    /// Fold another partial in. All fields are integer counts, so this is
    /// commutative + associative — the merge order does not change the totals (the
    /// determinism the parallel reduce relies on).
    pub(crate) fn merge(&mut self, other: &PrepassStats) {
        for (&period, profile) in &other.slip_by_period {
            merge_slip_profile(self.slip_by_period.entry(period).or_default(), profile);
        }
        for (&key, profile) in &other.slip_by_sample_period {
            merge_slip_profile(self.slip_by_sample_period.entry(key).or_default(), profile);
        }
        for (&sample, stats) in &other.per_sample {
            merge_sample_stats(self.per_sample.entry(sample).or_default(), stats);
        }
    }
}

/// Add `src`'s slip counts into `dst` (integer, order-independent).
fn merge_slip_profile(dst: &mut SlipProfile, src: &SlipProfile) {
    for k in 0..MAX_SLIP {
        dst.up[k] += src.up[k];
        dst.down[k] += src.down[k];
    }
}

/// The measured global parameters (pre-grouping; D3 splits these per sample group).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct EstimatedParams {
    /// Cohort per-base substitution rate `ε`.
    pub(crate) eps: f64,
    /// Per-period stutter shape (the `θ_period` parent).
    pub(crate) shape_by_period: HashMap<u8, StutterShape>,
    /// Per-sample stutter level line.
    pub(crate) level_by_sample: HashMap<u32, StutterLevel>,
}

/// Add `(faithful, slipped)` to the bin for `length`, creating it if new.
pub(crate) fn add_bin(
    by_length: &mut Vec<(u16, u64, u64)>,
    length: u16,
    faithful: u64,
    slipped: u64,
) {
    match by_length.iter_mut().find(|(l, _, _)| *l == length) {
        Some(bin) => {
            bin.1 += faithful;
            bin.2 += slipped;
        }
        None => by_length.push((length, faithful, slipped)),
    }
}

/// Match / mismatch base counts of an observed sequence against a same-length allele.
fn compare_bases(obs: &[u8], allele: &[u8]) -> (u64, u64) {
    let matches = obs.iter().zip(allele).filter(|(a, b)| a == b).count() as u64;
    (matches, obs.len() as u64 - matches)
}

/// A sample's per-base substitution rate `ε` from its faithful-peak base counts.
pub(crate) fn sample_eps(stats: &SampleStutterStats) -> f64 {
    let total = stats.base_match + stats.base_mismatch;
    if total > 0 {
        stats.base_mismatch as f64 / total as f64
    } else {
        0.0
    }
}

/// Fold a sample's per-locus stats into its running accumulator.
pub(crate) fn merge_sample_stats(dst: &mut SampleStutterStats, src: &SampleStutterStats) {
    dst.base_match += src.base_match;
    dst.base_mismatch += src.base_mismatch;
    dst.read_depth += src.read_depth;
    for &(length, faithful, slipped) in &src.by_length {
        add_bin(&mut dst.by_length, length, faithful, slipped);
    }
}

/// D1 — accumulate one locus's confident-genotype statistics into `stats`.
///
/// BIAS NOTE: `ε` is taken off a peak's *faithful* reads by comparing each to the
/// peak's representative allele. For a het this is mildly biased high — the other
/// allele's rare same-length (`≥2`-step) stutter and any impure same-length variant
/// both count as base mismatch. Negligible on homozygote-dominated cohorts; the
/// deferred soft full-cohort EM reduce removes it by fractional attribution. Slips
/// with `|Δ| > MAX_SLIP` are dropped (the spec cap; rare in practice).
pub(crate) fn accumulate_locus(
    stats: &mut PrepassStats,
    locus: &CohortLocus,
    rungs: &Rungs,
    ploidy: u8,
    cfg: &RungCfg,
) {
    let period = rungs.period();
    for (k, evidence) in locus.samples.iter().enumerate() {
        let global = locus.present[k];
        let Resolution::Confident(ResolvedGenotype::Peaks(peaks)) =
            resolve_confident_genotype(evidence, rungs, ploidy, cfg)
        else {
            continue;
        };

        // Accumulate locally first (avoids holding two &mut into `stats` at once).
        let mut local = SampleStutterStats::default();
        let mut slips: Vec<(i32, u32)> = Vec::new();
        for (obs, count) in &evidence.seq_counts {
            let read_units = (obs.len() / period) as u16;
            let peak = peaks
                .iter()
                .min_by_key(|p| (p.repeat_len as i32 - read_units as i32).abs())
                .expect("a confident genotype has ≥1 peak");
            let delta = read_units as i32 - peak.repeat_len as i32;
            if delta == 0 {
                let (m, mm) = compare_bases(obs, &peak.allele);
                local.base_match += m * *count as u64;
                local.base_mismatch += mm * *count as u64;
                add_bin(&mut local.by_length, peak.repeat_len, *count as u64, 0);
            } else {
                add_bin(&mut local.by_length, peak.repeat_len, 0, *count as u64);
                slips.push((delta, *count));
            }
            local.read_depth += *count as u64;
        }

        merge_sample_stats(stats.per_sample.entry(global).or_default(), &local);
        for (delta, count) in &slips {
            add_slip(
                stats.slip_by_period.entry(period as u8).or_default(),
                *delta,
                *count,
            );
            add_slip(
                stats
                    .slip_by_sample_period
                    .entry((global, period as u8))
                    .or_default(),
                *delta,
                *count,
            );
        }
    }
}

/// Add a slip of signed size `delta` (count `count`) to a profile, respecting the
/// `MAX_SLIP` cap.
fn add_slip(profile: &mut SlipProfile, delta: i32, count: u32) {
    let magnitude = delta.unsigned_abs() as usize;
    if (1..=MAX_SLIP).contains(&magnitude) {
        if delta > 0 {
            profile.up[magnitude - 1] += count as u64;
        } else {
            profile.down[magnitude - 1] += count as u64;
        }
    }
}

/// A mild fallback shape when a period has no observed slips.
const FALLBACK_SHAPE: StutterShape = StutterShape {
    up_rate: 0.5,
    down_rate: 0.5,
    decay: 0.1,
};

/// Estimate a stutter shape from a slip profile: the direction split + the
/// geometric-decay MLE `decay = (mean magnitude − 1) / mean magnitude`.
fn estimate_shape(profile: &SlipProfile) -> StutterShape {
    let up: u64 = profile.up.iter().sum();
    let down: u64 = profile.down.iter().sum();
    let total = up + down;
    if total == 0 {
        return FALLBACK_SHAPE;
    }
    let up_fraction = up as f64 / total as f64;
    let mut weighted_magnitude = 0.0;
    for k in 0..MAX_SLIP {
        weighted_magnitude += (k + 1) as f64 * (profile.up[k] + profile.down[k]) as f64;
    }
    let mean_magnitude = weighted_magnitude / total as f64;
    let decay = ((mean_magnitude - 1.0) / mean_magnitude).clamp(0.0, 0.999);
    StutterShape {
        up_rate: up_fraction,
        down_rate: 1.0 - up_fraction,
        decay,
    }
}

/// Fit a sample's stutter level line `baseline + slope·length` from its
/// faithful/slipped counts by length (weighted least squares; constant fit when a
/// single length is observed).
///
/// The returned line is a **seed** and is NOT clamped here: a weighted-LS fit can
/// extrapolate to a negative or `>1` rate at lengths outside the observed range.
/// Every consumer clamps `level` into `[0,1]` at evaluation (`em.rs`, `sim.rs`), so
/// the raw line is intentionally preserved for the outer-loop refit (E1).
pub(crate) fn fit_level(stats: &SampleStutterStats) -> StutterLevel {
    let points: Vec<(f64, f64, f64)> = stats
        .by_length
        .iter()
        .filter(|(_, faithful, slipped)| faithful + slipped > 0)
        .map(|(length, faithful, slipped)| {
            let total = (faithful + slipped) as f64;
            (*length as f64, *slipped as f64 / total, total)
        })
        .collect();

    if points.len() < 2 {
        let total: u64 = stats.by_length.iter().map(|(_, f, s)| f + s).sum();
        let slipped: u64 = stats.by_length.iter().map(|(_, _, s)| *s).sum();
        let baseline = if total > 0 {
            slipped as f64 / total as f64
        } else {
            0.0
        };
        return StutterLevel {
            baseline,
            slope: 0.0,
        };
    }

    // Weighted least squares for rate ~ baseline + slope·length.
    let (mut sw, mut swx, mut swy, mut swxx, mut swxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (x, y, w) in points {
        sw += w;
        swx += w * x;
        swy += w * y;
        swxx += w * x * x;
        swxy += w * x * y;
    }
    let denom = sw * swxx - swx * swx;
    if denom.abs() < f64::EPSILON {
        return StutterLevel {
            baseline: swy / sw,
            slope: 0.0,
        };
    }
    let slope = (sw * swxy - swx * swy) / denom;
    let baseline = (swy - slope * swx) / sw;
    StutterLevel { baseline, slope }
}

/// D2 — measure the parameters from the accumulated statistics.
pub(crate) fn estimate(stats: &PrepassStats) -> EstimatedParams {
    let (mut matches, mut mismatches) = (0u64, 0u64);
    for sample in stats.per_sample.values() {
        matches += sample.base_match;
        mismatches += sample.base_mismatch;
    }
    let eps = if matches + mismatches > 0 {
        mismatches as f64 / (matches + mismatches) as f64
    } else {
        0.0
    };

    let shape_by_period = stats
        .slip_by_period
        .iter()
        .map(|(period, profile)| (*period, estimate_shape(profile)))
        .collect();
    let level_by_sample = stats
        .per_sample
        .iter()
        .map(|(sample, st)| (*sample, fit_level(st)))
        .collect();

    EstimatedParams {
        eps,
        shape_by_period,
        level_by_sample,
    }
}

/// Accumulate the pre-pass sufficient statistics over the cohort loci (parallel).
///
/// Each thread accumulates a per-thread [`PrepassStats`] partial; the partials are
/// reduced with [`PrepassStats::merge`]. The statistics are integer counts, so the
/// reduce is order-independent — the result is byte-identical across thread counts
/// (F1), and `estimate` reads off it deterministically regardless of any internal
/// `Vec`/`HashMap` ordering.
pub(crate) fn run_prepass_stats(loci: &[CohortLocus], ploidy: u8, cfg: &RungCfg) -> PrepassStats {
    loci.par_iter()
        .fold(PrepassStats::default, |mut acc, locus| {
            let rungs = build_rungs(locus, cfg);
            accumulate_locus(&mut acc, locus, &rungs, ploidy, cfg);
            acc
        })
        .reduce(PrepassStats::default, |mut a, b| {
            a.merge(&b);
            a
        })
}

/// Run the pre-pass over the cohort loci and return the measured parameters.
pub(crate) fn run_prepass(loci: &[CohortLocus], ploidy: u8, cfg: &RungCfg) -> EstimatedParams {
    estimate(&run_prepass_stats(loci, ploidy, cfg))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::param_estimation::{PerBaseError, SampleGroupId};
    use crate::ssr::cohort::sim::{
        SimChemistry, SimCohortSpec, SimGenotype, SimLocus, SimSample, simulate,
    };
    use crate::ssr::types::Motif;

    /// A cohort of homozygotes spread over lengths 7..=12 (4 samples each, so every
    /// allele is cohort-recurrent), under a known chemistry.
    fn recovery_spec(eps: f64, up: f64, down: f64, decay: f64, level: f64) -> SimCohortSpec {
        let chem = SimChemistry {
            error: PerBaseError(eps),
            shape: StutterShape {
                up_rate: up,
                down_rate: down,
                decay,
            },
            level: StutterLevel {
                baseline: level,
                slope: 0.0,
            },
        };
        let mut samples = Vec::new();
        for length in 7u16..=12 {
            for replicate in 0..4 {
                let idx = samples.len();
                samples.push(SimSample {
                    name: format!("S{idx}_L{length}_{replicate}"),
                    group: SampleGroupId(0),
                    genotypes: vec![Some(SimGenotype::homozygous(length, 2))],
                });
            }
        }
        SimCohortSpec {
            seed: 4242,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 60,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 10,
            }],
            groups: vec![chem],
            samples,
            depth: 200,
        }
    }

    fn collect_loci(spec: &SimCohortSpec) -> Vec<CohortLocus> {
        simulate(spec)
            .merger()
            .map(|item| item.expect("locus"))
            .map(|(_, locus)| locus)
            .collect()
    }

    fn run(spec: &SimCohortSpec) -> EstimatedParams {
        run_prepass(&collect_loci(spec), 2, &RungCfg::dev_default())
    }

    /// A cohort of cohort-recurrent separated hets (6/10), under a known chemistry.
    fn het_spec() -> SimCohortSpec {
        let chem = SimChemistry {
            error: PerBaseError(0.004),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.10,
                slope: 0.0,
            },
        };
        let samples = (0..8)
            .map(|i| SimSample {
                name: format!("H{i}"),
                group: SampleGroupId(0),
                genotypes: vec![Some(SimGenotype::diploid(6, 10))],
            })
            .collect();
        SimCohortSpec {
            seed: 99,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 60,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 10,
            }],
            groups: vec![chem],
            samples,
            depth: 200,
        }
    }

    #[test]
    fn separated_hets_contribute_two_length_bins() {
        // A 6/10 het must deposit faithful/slipped stats at BOTH allele lengths.
        let loci = collect_loci(&het_spec());
        let mut stats = PrepassStats::default();
        for locus in &loci {
            let rungs = build_rungs(locus, &RungCfg::dev_default());
            accumulate_locus(&mut stats, locus, &rungs, 2, &RungCfg::dev_default());
        }
        // Every resolved het sample carries bins at lengths 6 and 10.
        assert!(!stats.per_sample.is_empty(), "hets should resolve");
        for sample in stats.per_sample.values() {
            let mut lengths: Vec<u16> = sample.by_length.iter().map(|(l, _, _)| *l).collect();
            lengths.sort_unstable();
            assert_eq!(
                lengths,
                vec![6, 10],
                "a 6/10 het should bin at both alleles"
            );
        }
        // ε is still recovered off the (cleaner outer) het skirts.
        let est = estimate(&stats);
        assert!((est.eps - 0.004).abs() < 0.003, "ε recovered {}", est.eps);
    }

    fn with_threads<T: Send>(n: usize, f: impl FnOnce() -> T + Send) -> T {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(f)
    }

    #[test]
    fn prepass_is_byte_identical_across_thread_counts() {
        let loci = collect_loci(&recovery_spec(0.004, 1.0, 2.0, 0.1, 0.10));
        let single = with_threads(1, || run_prepass(&loci, 2, &RungCfg::dev_default()));
        let multi = with_threads(4, || run_prepass(&loci, 2, &RungCfg::dev_default()));
        assert_eq!(
            single, multi,
            "the pre-pass must be byte-identical across thread counts (F1)"
        );
    }

    #[test]
    fn run_prepass_is_deterministic() {
        let loci = collect_loci(&recovery_spec(0.004, 1.0, 2.0, 0.1, 0.10));
        let a = run_prepass(&loci, 2, &RungCfg::dev_default());
        let b = run_prepass(&loci, 2, &RungCfg::dev_default());
        assert_eq!(
            a, b,
            "the pre-pass must be order-independent / reproducible"
        );
    }

    #[test]
    fn checkpoint_2_recovers_epsilon_shape_and_level() {
        let (eps, up, down, decay, level) = (0.004, 1.0, 2.0, 0.1, 0.10);
        let est = run(&recovery_spec(eps, up, down, decay, level));

        // ε within a few parts in 1000.
        assert!(
            (est.eps - eps).abs() < 0.002,
            "ε recovered {} vs injected {eps}",
            est.eps
        );

        // Shape: contraction bias (down > up), direction split ≈ 1:2, decay ≈ 0.1.
        let shape = est.shape_by_period[&2];
        assert!(
            shape.down_rate > shape.up_rate,
            "should recover contraction bias"
        );
        let injected_up_frac = up / (up + down); // 1/3
        assert!(
            (shape.up_rate - injected_up_frac).abs() < 0.08,
            "up fraction recovered {} vs {injected_up_frac}",
            shape.up_rate
        );
        assert!(
            (shape.decay - decay).abs() < 0.06,
            "decay recovered {} vs {decay}",
            shape.decay
        );

        // Level: the mean per-sample baseline ≈ the injected level.
        let mean_baseline: f64 = est
            .level_by_sample
            .values()
            .map(|l| l.baseline)
            .sum::<f64>()
            / est.level_by_sample.len() as f64;
        assert!(
            (mean_baseline - level).abs() < 0.03,
            "mean level baseline recovered {mean_baseline} vs {level}"
        );
    }

    #[test]
    fn recovers_a_higher_epsilon() {
        let est = run(&recovery_spec(0.02, 1.0, 2.0, 0.1, 0.08));
        assert!((est.eps - 0.02).abs() < 0.005, "ε recovered {}", est.eps);
    }

    #[test]
    fn recovers_a_higher_stutter_level() {
        let est = run(&recovery_spec(0.004, 1.0, 2.0, 0.1, 0.25));
        let mean: f64 = est
            .level_by_sample
            .values()
            .map(|l| l.baseline)
            .sum::<f64>()
            / est.level_by_sample.len() as f64;
        assert!((mean - 0.25).abs() < 0.04, "level recovered {mean}");
    }

    #[test]
    fn estimate_is_empty_without_confident_genotypes() {
        // A single thin sample → no confident genotype → no per-sample stats.
        let chem = SimChemistry {
            error: PerBaseError(0.004),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay: 0.1,
            },
            level: StutterLevel {
                baseline: 0.1,
                slope: 0.0,
            },
        };
        let spec = SimCohortSpec {
            seed: 1,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 60,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 10,
            }],
            groups: vec![chem],
            samples: vec![SimSample {
                name: "S0".into(),
                group: SampleGroupId(0),
                genotypes: vec![Some(SimGenotype::homozygous(8, 2))],
            }],
            depth: 3, // below the resolution depth floor
        };
        let est = run(&spec);
        assert_eq!(est.eps, 0.0);
        assert!(est.level_by_sample.is_empty());
    }
}
