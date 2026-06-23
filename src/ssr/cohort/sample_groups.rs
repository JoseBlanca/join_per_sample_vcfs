//! Sample-group clustering + per-`(group, period)` shape + the `ε`-freeze check
//! (arch `ssr_call_parameters.md` §4, spec §4.4, implementation plan §4 D3 — the M3
//! milestone).
//!
//! With the per-sample `(ε, level)` from D2, this:
//! 1. **clusters** samples into data-driven sample groups — a deterministic
//!    distance-based agglomeration (union-find on close neighbours in scaled
//!    `(ε, level)` space; no k-means, no random init; the group *count* falls out of
//!    the threshold, not a preset K). A single-protocol cohort collapses to one
//!    group; a two-protocol cohort separates;
//! 2. with the groups fixed, fits the per-`(group, period)` stutter **shape** by
//!    re-aggregating the per-sample slip profiles within each group and shrinking
//!    toward the cohort `θ_period` parent (`refine_theta_locus`), plus the per-group
//!    `ε` and level;
//! 3. runs the **`ε`-freeze check** — a binomial BIC comparison of frozen-per-group
//!    `ε` vs. per-sample `ε`; if the richer per-sample model is not justified, the
//!    freeze stands.
//!
//! Scope (deferred to F): uncertainty-scaled distances beyond the basic
//! depth weighting, the BIC-selected group *count* (here a calibrated threshold),
//! and the singleton-group shrinkage corner. The clustering is deterministic — ties
//! resolve on the sample catalog index.

use std::collections::HashMap;

use crate::ssr::cohort::param_estimation::{
    SampleGroupId, SampleStutterStats, SlipProfile, StutterLevel, StutterShape,
};
use crate::ssr::cohort::prepass::{
    EstimatedParams, PrepassStats, fit_level, merge_sample_stats, sample_eps,
};
use crate::ssr::cohort::stutter::refine_theta_locus;

/// Clustering + per-group fit controls (pinned in F2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ClusterCfg {
    /// Characteristic `ε` scale for the distance metric.
    pub(crate) eps_scale: f64,
    /// Characteristic level scale for the distance metric.
    pub(crate) level_scale: f64,
    /// Two samples merge when their scaled distance is below this.
    pub(crate) distance_threshold: f64,
    /// Reference depth at which the precision weighting saturates (above it, full
    /// distance is used; below it, distances deflate so a thin sample merges readily).
    pub(crate) precision_reference_depth: f64,
    /// Shrinkage strength of the per-`(group, period)` shape toward `θ_period`.
    pub(crate) shape_shrink_strength: f64,
}

impl ClusterCfg {
    pub(crate) fn dev_default() -> Self {
        Self {
            eps_scale: 0.005,
            level_scale: 0.05,
            distance_threshold: 1.5,
            precision_reference_depth: 100.0,
            shape_shrink_strength: 50.0,
        }
    }
}

/// The grouped chemistry: the sample → group map plus per-group parameters.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GroupedParams {
    /// Global sample index → its sample group.
    pub(crate) group_of_sample: HashMap<u32, SampleGroupId>,
    /// Number of sample groups found.
    pub(crate) n_groups: usize,
    /// Per-group `ε` (indexed by group id).
    pub(crate) eps_per_group: Vec<f64>,
    /// Per-group stutter level line (indexed by group id).
    pub(crate) level_per_group: Vec<StutterLevel>,
    /// Per-`(group, period)` stutter shape, shrunk toward the `θ_period` parent.
    pub(crate) shape_by_group_period: HashMap<(SampleGroupId, u8), StutterShape>,
    /// Whether freezing `ε` per group (vs. per sample) is justified by a binomial BIC.
    pub(crate) eps_freeze_justified: bool,
}

/// Scaled distance between two samples in `(ε, level)` space, deflated by the less
/// certain sample's precision (a thin/wobbly sample is "close to many").
fn scaled_distance(a: (f64, f64, u64), b: (f64, f64, u64), cfg: &ClusterCfg) -> f64 {
    let (eps_a, level_a, depth_a) = a;
    let (eps_b, level_b, depth_b) = b;
    let de = (eps_a - eps_b) / cfg.eps_scale;
    let dl = (level_a - level_b) / cfg.level_scale;
    let raw = (de * de + dl * dl).sqrt();
    // Precision ∝ √depth, normalized to the reference depth. The less-certain (thinner)
    // sample's precision deflates the distance, so a wobbly thin sample is "close to
    // many things" and merges readily; above the reference depth it saturates at 1.0.
    let precision = (depth_a.min(depth_b) as f64).sqrt() / cfg.precision_reference_depth.sqrt();
    raw * precision.min(1.0)
}

/// Union-find root with path halving.
fn find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

/// Cluster the per-sample features into groups (deterministic union-find over
/// sorted sample indices). Returns `(group_of_sample, n_groups)`.
///
/// LIMITATION: this is single-linkage — a smooth `(ε, level)` continuum whose
/// neighbours are each within the threshold chains into one group even when the
/// extremes are far apart. The penalized-likelihood split test that breaks an
/// over-merged group (spec §4.4) is deferred (F); for well-separated protocols (the
/// common case) single-linkage is correct and cheap.
fn cluster(
    features: &[(u32, f64, f64, u64)],
    cfg: &ClusterCfg,
) -> (HashMap<u32, SampleGroupId>, usize) {
    let n = features.len();
    let mut parent: Vec<usize> = (0..n).collect();
    for i in 0..n {
        for j in (i + 1)..n {
            let di = (features[i].1, features[i].2, features[i].3);
            let dj = (features[j].1, features[j].2, features[j].3);
            if scaled_distance(di, dj, cfg) < cfg.distance_threshold {
                let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                if ri != rj {
                    // Union toward the smaller root index (deterministic).
                    let (lo, hi) = (ri.min(rj), ri.max(rj));
                    parent[hi] = lo;
                }
            }
        }
    }

    // Assign group ids by ascending root index (features are already sorted by
    // sample index, so this is a stable, catalog-order labelling).
    let mut root_to_group: HashMap<usize, u16> = HashMap::new();
    let mut group_of_sample = HashMap::new();
    for (idx, feature) in features.iter().enumerate() {
        let root = find(&mut parent, idx);
        let next = root_to_group.len() as u16;
        let group = *root_to_group.entry(root).or_insert(next);
        group_of_sample.insert(feature.0, SampleGroupId(group));
    }
    let n_groups = root_to_group.len();
    (group_of_sample, n_groups)
}

/// Binomial log-likelihood of `(match, mismatch)` base counts under rate `eps`.
fn binomial_loglik(matches: u64, mismatches: u64, eps: f64) -> f64 {
    let p_match = (1.0 - eps).clamp(1e-12, 1.0);
    let p_mismatch = eps.clamp(1e-12, 1.0);
    matches as f64 * p_match.ln() + mismatches as f64 * p_mismatch.ln()
}

/// Group the samples and fit the per-group chemistry + per-`(group, period)` shape.
pub(crate) fn group_samples(
    stats: &PrepassStats,
    est: &EstimatedParams,
    cfg: &ClusterCfg,
) -> GroupedParams {
    // Per-sample features in ascending sample-index order (determinism).
    let mut samples: Vec<u32> = stats.per_sample.keys().copied().collect();
    samples.sort_unstable();
    let features: Vec<(u32, f64, f64, u64)> = samples
        .iter()
        .map(|&s| {
            let st = &stats.per_sample[&s];
            let level = est
                .level_by_sample
                .get(&s)
                .map(|l| l.baseline)
                .unwrap_or_else(|| fit_level(st).baseline);
            (s, sample_eps(st), level, st.read_depth)
        })
        .collect();

    let (group_of_sample, n_groups) = cluster(&features, cfg);

    // Aggregate per group: stutter stats (→ ε, level) and slip profiles by period.
    let mut group_stats: Vec<SampleStutterStats> = vec![SampleStutterStats::default(); n_groups];
    let mut group_period_slips: HashMap<(SampleGroupId, u8), SlipProfile> = HashMap::new();
    for &s in &samples {
        let group = group_of_sample[&s];
        merge_sample_stats(&mut group_stats[group.0 as usize], &stats.per_sample[&s]);
    }
    for (&(sample, period), profile) in &stats.slip_by_sample_period {
        if let Some(&group) = group_of_sample.get(&sample) {
            let entry = group_period_slips.entry((group, period)).or_default();
            for k in 0..profile.up.len() {
                entry.up[k] += profile.up[k];
                entry.down[k] += profile.down[k];
            }
        }
    }

    let eps_per_group: Vec<f64> = group_stats.iter().map(sample_eps).collect();
    let level_per_group: Vec<StutterLevel> = group_stats.iter().map(fit_level).collect();

    // Per-(group, period) shape: the group slip profile shrunk toward θ_period.
    let shape_by_group_period = group_period_slips
        .iter()
        .map(|(&(group, period), profile)| {
            let parent = est
                .shape_by_period
                .get(&period)
                .copied()
                .unwrap_or(StutterShape {
                    up_rate: 0.5,
                    down_rate: 0.5,
                    decay: 0.1,
                });
            let shape = refine_theta_locus(profile, &parent, cfg.shape_shrink_strength);
            ((group, period), shape)
        })
        .collect();

    // ε-freeze check: per-group frozen ε vs. per-sample ε via binomial BIC.
    let eps_freeze_justified = eps_freeze_check(stats, &group_of_sample, &eps_per_group, n_groups);

    GroupedParams {
        group_of_sample,
        n_groups,
        eps_per_group,
        level_per_group,
        shape_by_group_period,
        eps_freeze_justified,
    }
}

/// Whether freezing `ε` per group beats a per-sample `ε` model on a binomial BIC
/// (true = the freeze is justified, i.e. the richer model does not earn its params).
fn eps_freeze_check(
    stats: &PrepassStats,
    group_of_sample: &HashMap<u32, SampleGroupId>,
    eps_per_group: &[f64],
    n_groups: usize,
) -> bool {
    let mut frozen_ll = 0.0;
    let mut per_sample_ll = 0.0;
    let mut total_bases = 0u64;
    for (&sample, st) in &stats.per_sample {
        let group = group_of_sample[&sample];
        frozen_ll += binomial_loglik(
            st.base_match,
            st.base_mismatch,
            eps_per_group[group.0 as usize],
        );
        per_sample_ll += binomial_loglik(st.base_match, st.base_mismatch, sample_eps(st));
        total_bases += st.base_match + st.base_mismatch;
    }
    let n_samples = stats.per_sample.len();
    let ln_n = (total_bases.max(1) as f64).ln();
    let bic_frozen = -2.0 * frozen_ll + n_groups as f64 * ln_n;
    let bic_per_sample = -2.0 * per_sample_ll + n_samples as f64 * ln_n;
    bic_frozen <= bic_per_sample
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::param_estimation::PerBaseError;
    use crate::ssr::cohort::prepass::{estimate, run_prepass_stats};
    use crate::ssr::cohort::rung_ladder::RungCfg;
    use crate::ssr::cohort::sim::{
        SimChemistry, SimCohortSpec, SimGenotype, SimLocus, SimSample, simulate,
    };
    use crate::ssr::cohort::types::CohortLocus;
    use crate::ssr::types::Motif;

    fn chem(eps: f64, level: f64, decay: f64) -> SimChemistry {
        SimChemistry {
            error: PerBaseError(eps),
            shape: StutterShape {
                up_rate: 1.0,
                down_rate: 2.0,
                decay,
            },
            level: StutterLevel {
                baseline: level,
                slope: 0.0,
            },
        }
    }

    /// Homozygous samples spread over recurrent lengths, split across `groups` by a
    /// round-robin over `group_assignment`.
    fn grouped_spec(groups: Vec<SimChemistry>, group_of: &[usize]) -> SimCohortSpec {
        let mut samples = Vec::new();
        // For each group, four samples at each of lengths 7..=11 (cohort-recurrent).
        for (sample_group, _) in group_of.iter().enumerate() {
            for length in 7u16..=11 {
                for _ in 0..3 {
                    let idx = samples.len();
                    samples.push(SimSample {
                        name: format!("S{idx}"),
                        group: SampleGroupId(group_of[sample_group] as u16),
                        genotypes: vec![Some(SimGenotype::homozygous(length, 2))],
                    });
                }
            }
        }
        SimCohortSpec {
            seed: 2024,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 60,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 9,
            }],
            groups,
            samples,
            depth: 250,
        }
    }

    fn prepare(spec: &SimCohortSpec) -> (PrepassStats, EstimatedParams) {
        let loci: Vec<CohortLocus> = simulate(spec)
            .merger()
            .map(|item| item.expect("locus"))
            .map(|(_, locus)| locus)
            .collect();
        let stats = run_prepass_stats(&loci, 2, &RungCfg::dev_default());
        let est = estimate(&stats);
        (stats, est)
    }

    #[test]
    fn single_protocol_cohort_collapses_to_one_group() {
        let spec = grouped_spec(vec![chem(0.004, 0.10, 0.1)], &[0]);
        let (stats, est) = prepare(&spec);
        let grouped = group_samples(&stats, &est, &ClusterCfg::dev_default());
        assert_eq!(grouped.n_groups, 1, "one chemistry → one group");
    }

    #[test]
    fn m3_two_protocol_cohort_recovers_per_group_shapes() {
        // Group 0: low ε, low level, shallow decay. Group 1: high ε, high level,
        // heavy-tailed decay. Distinct in BOTH chemistry and shape.
        let spec = grouped_spec(
            vec![chem(0.002, 0.06, 0.05), chem(0.012, 0.22, 0.45)],
            &[0, 1],
        );
        let (stats, est) = prepare(&spec);
        let grouped = group_samples(&stats, &est, &ClusterCfg::dev_default());

        assert_eq!(grouped.n_groups, 2, "two protocols → two groups");

        // The two groups' ε differ in the injected direction.
        let mut eps = grouped.eps_per_group.clone();
        eps.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!(eps[0] < 0.006 && eps[1] > 0.008, "group ε: {eps:?}");

        // The per-(group, period) shapes recover the divergent decay: one shallow,
        // one heavy-tailed.
        let mut decays: Vec<f64> = grouped
            .shape_by_group_period
            .values()
            .map(|s| s.decay)
            .collect();
        decays.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(decays.len(), 2, "one shape per group at period 2");
        assert!(
            decays[0] < 0.2 && decays[1] > 0.3,
            "per-group decays should stay divergent: {decays:?}"
        );
    }

    #[test]
    fn group_samples_is_deterministic() {
        let spec = grouped_spec(
            vec![chem(0.002, 0.06, 0.05), chem(0.012, 0.22, 0.45)],
            &[0, 1],
        );
        let (stats, est) = prepare(&spec);
        let a = group_samples(&stats, &est, &ClusterCfg::dev_default());
        let b = group_samples(&stats, &est, &ClusterCfg::dev_default());
        assert_eq!(a, b, "clustering + per-group fit must be reproducible");
    }

    #[test]
    fn frozen_epsilon_is_justified_for_a_single_protocol() {
        let spec = grouped_spec(vec![chem(0.004, 0.10, 0.1)], &[0]);
        let (stats, est) = prepare(&spec);
        let grouped = group_samples(&stats, &est, &ClusterCfg::dev_default());
        assert!(
            grouped.eps_freeze_justified,
            "one ε for one homogeneous group should beat a per-sample ε model"
        );
    }
}
