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

use crate::ssr::cohort::attribution::nearest_parent;
use crate::ssr::cohort::param_estimation::{
    AlleleSpreadAccum, G0FitCfg, G0PseudocountDecay, LengthBin, MAX_SLIP, SampleStutterStats,
    SlipProfile, StutterLevel, StutterShape,
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
    /// Period → the count-weighted spread of confident alleles around the per-locus
    /// mode, over **variable loci only** — the sufficient statistic for the `G₀`
    /// decay fit (arch `ssr_call_driver.md` §9).
    pub(crate) allele_spread_by_period: HashMap<u8, AlleleSpreadAccum>,
}

impl PrepassStats {
    /// Fold another partial in. All fields are integer counts, so this is
    /// commutative + associative — the merge order does not change the totals (the
    /// determinism the parallel reduce relies on).
    ///
    /// DETERMINISM LAYER: the *values* (and everything `estimate` derives from them —
    /// `ε`, the per-period shape, the per-sample level) are merge-order-independent and
    /// byte-identical across thread counts. The *struct itself* is **not** a canonical
    /// form: `SampleStutterStats.by_length` is a `Vec` whose order follows insertion /
    /// merge order, so two `PrepassStats` reduced at different thread counts can differ
    /// under derived `PartialEq`. Determinism tests must compare `EstimatedParams` /
    /// `CohortCalls`, not the raw stats.
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
        for (&period, accum) in &other.allele_spread_by_period {
            let dst = self.allele_spread_by_period.entry(period).or_default();
            dst.total_copies += accum.total_copies;
            dst.distance_weighted += accum.distance_weighted;
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
    /// Per loci group (= period) `G₀` decay `p`, fit from the confident-allele spread
    /// (arch `ssr_call_driver.md` §9). One entry per period observed among the
    /// variable loci; thin periods carry the coded fallback.
    pub(crate) g0_by_period: HashMap<u8, G0PseudocountDecay>,
}

/// Add `(faithful, slipped)` to the bin for `length`, creating it if new.
pub(crate) fn add_bin(by_length: &mut Vec<LengthBin>, length: u16, faithful: u64, slipped: u64) {
    match by_length.iter_mut().find(|bin| bin.length == length) {
        Some(bin) => {
            bin.faithful += faithful;
            bin.slipped += slipped;
        }
        None => by_length.push(LengthBin {
            length,
            faithful,
            slipped,
        }),
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
    for bin in &src.by_length {
        add_bin(&mut dst.by_length, bin.length, bin.faithful, bin.slipped);
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
    // Cohort-wide confident-allele tally for this locus (length in repeat units →
    // chromosome copies), used for the G₀ allele-spread fit below.
    let mut locus_alleles: Vec<AlleleCopies> = Vec::new();
    for (k, evidence) in locus.samples.iter().enumerate() {
        let global = locus.present[k];
        let Resolution::Confident(ResolvedGenotype::Peaks(peaks)) =
            resolve_confident_genotype(evidence, rungs, ploidy, cfg)
        else {
            continue;
        };

        // Each peak is one distinct called allele; a 1-peak (homozygous) call carries
        // `ploidy` chromosome copies, a `ploidy`-peak call one copy each. This integer
        // split is exact only when the peak count divides the ploidy — true for the
        // diploid contract (1 or 2 peaks); a future polyploid with a non-dividing peak
        // count would need fractional/dosage attribution, so guard it here.
        debug_assert!(
            (ploidy as usize).is_multiple_of(peaks.len()),
            "peak count {} does not divide ploidy {ploidy} (diploid-only contract)",
            peaks.len()
        );
        let copies_per_allele = ploidy as u64 / peaks.len() as u64;
        for peak in &peaks {
            add_allele_copies(&mut locus_alleles, peak.repeat_len, copies_per_allele);
        }

        // Accumulate locally first (avoids holding two &mut into `stats` at once).
        let mut local = SampleStutterStats::default();
        let mut slips: Vec<(i32, u32)> = Vec::new();
        let peak_units: Vec<u16> = peaks.iter().map(|p| p.repeat_len).collect();
        for (obs, count) in &evidence.seq_counts {
            let read_units = (obs.len() / period) as i32;
            let (peak_idx, delta) =
                nearest_parent(read_units, &peak_units).expect("a confident genotype has ≥1 peak");
            let peak = &peaks[peak_idx];
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
            stats
                .slip_by_period
                .entry(period as u8)
                .or_default()
                .add_slip(*delta, *count as u64);
            stats
                .slip_by_sample_period
                .entry((global, period as u8))
                .or_default()
                .add_slip(*delta, *count as u64);
        }
    }

    // G₀ spread: only loci that VARY (≥2 distinct confident alleles cohort-wide)
    // inform the prior — "given a locus varies, how spread are its alleles." Distance
    // is measured from the per-locus modal length, the same reference `g0_pseudocounts`
    // uses at call time.
    if locus_alleles.len() >= 2
        && let Some(mode) = rungs.modal_length()
    {
        let accum = stats
            .allele_spread_by_period
            .entry(period as u8)
            .or_default();
        for entry in &locus_alleles {
            let distance = (entry.length as i32 - mode as i32).unsigned_abs() as u64;
            accum.total_copies += entry.copies;
            accum.distance_weighted += entry.copies * distance;
        }
    }
}

/// One repeat-length entry of a locus's cohort-wide confident-allele tally: how many
/// chromosome copies carry `length` repeat units. Feeds the `G₀` allele-spread fit; a
/// named struct over a `(u16, u64)` tuple for the same reason as [`LengthBin`] (Mi5).
#[derive(Debug, Clone, Copy)]
struct AlleleCopies {
    length: u16,
    copies: u64,
}

/// Add `copies` chromosome copies of repeat length `length` to a per-locus allele
/// tally, creating the entry if new.
fn add_allele_copies(tally: &mut Vec<AlleleCopies>, length: u16, copies: u64) {
    match tally.iter_mut().find(|entry| entry.length == length) {
        Some(entry) => entry.copies += copies,
        None => tally.push(AlleleCopies { length, copies }),
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
        .filter(|bin| bin.faithful + bin.slipped > 0)
        .map(|bin| {
            let total = (bin.faithful + bin.slipped) as f64;
            (bin.length as f64, bin.slipped as f64 / total, total)
        })
        .collect();

    if points.len() < 2 {
        let total: u64 = stats
            .by_length
            .iter()
            .map(|bin| bin.faithful + bin.slipped)
            .sum();
        let slipped: u64 = stats.by_length.iter().map(|bin| bin.slipped).sum();
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

/// Fit one period's `G₀` decay `p` from its confident-allele spread (arch
/// `ssr_call_driver.md` §9): the closed-form MLE of the symmetric geometric `p^|Δ|`,
/// `p = (√(1 + K̄²) − 1) / K̄` with `K̄` the mean distance-from-mode. A thin period
/// (too few copies) takes the coded fallback; the fit is clamped no steeper than
/// `min_p`. With `K̄ = 0` (every allele on the mode) the formula's limit is `0`, so it
/// clamps straight to `min_p`.
fn fit_g0_decay(accum: &AlleleSpreadAccum, cfg: &G0FitCfg) -> G0PseudocountDecay {
    // `.max(1)` keeps the divisor positive below even if a caller sets `min_copies = 0`
    // — `total_copies` is then guaranteed ≥ 1, so `K̄` can never be a `0/0` NaN.
    if accum.total_copies < cfg.min_copies.max(1) {
        return G0PseudocountDecay { p: cfg.fallback_p };
    }
    let k_bar = accum.distance_weighted as f64 / accum.total_copies as f64;
    let p = if k_bar > 0.0 {
        ((1.0 + k_bar * k_bar).sqrt() - 1.0) / k_bar
    } else {
        0.0
    };
    G0PseudocountDecay {
        p: p.clamp(cfg.min_p, 1.0),
    }
}

/// D2 — measure the parameters from the accumulated statistics.
pub(crate) fn estimate(stats: &PrepassStats, g0_cfg: &G0FitCfg) -> EstimatedParams {
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
    let g0_by_period = stats
        .allele_spread_by_period
        .iter()
        .map(|(period, accum)| (*period, fit_g0_decay(accum, g0_cfg)))
        .collect();

    EstimatedParams {
        eps,
        shape_by_period,
        level_by_sample,
        g0_by_period,
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
pub(crate) fn run_prepass(
    loci: &[CohortLocus],
    ploidy: u8,
    cfg: &RungCfg,
    g0_cfg: &G0FitCfg,
) -> EstimatedParams {
    estimate(&run_prepass_stats(loci, ploidy, cfg), g0_cfg)
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
        run_prepass(
            &collect_loci(spec),
            2,
            &RungCfg::dev_default(),
            &G0FitCfg::dev_default(),
        )
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
            let mut lengths: Vec<u16> = sample.by_length.iter().map(|bin| bin.length).collect();
            lengths.sort_unstable();
            assert_eq!(
                lengths,
                vec![6, 10],
                "a 6/10 het should bin at both alleles"
            );
        }
        // ε is still recovered off the (cleaner outer) het skirts.
        let est = estimate(&stats, &G0FitCfg::dev_default());
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
        let single = with_threads(1, || {
            run_prepass(&loci, 2, &RungCfg::dev_default(), &G0FitCfg::dev_default())
        });
        let multi = with_threads(4, || {
            run_prepass(&loci, 2, &RungCfg::dev_default(), &G0FitCfg::dev_default())
        });
        assert_eq!(
            single, multi,
            "the pre-pass must be byte-identical across thread counts (F1)"
        );
    }

    #[test]
    fn run_prepass_is_deterministic() {
        let loci = collect_loci(&recovery_spec(0.004, 1.0, 2.0, 0.1, 0.10));
        let a = run_prepass(&loci, 2, &RungCfg::dev_default(), &G0FitCfg::dev_default());
        let b = run_prepass(&loci, 2, &RungCfg::dev_default(), &G0FitCfg::dev_default());
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

    // ── G1: the G₀ decay fit ──

    #[test]
    fn g0_fit_recovers_the_closed_form_from_mean_distance() {
        // K̄ = 16/32 = 0.5 → p = (√(1+K̄²) − 1)/K̄. min_p = 0 so the clamp can't mask it.
        let cfg = G0FitCfg {
            min_copies: 1,
            fallback_p: 0.5,
            min_p: 0.0,
        };
        let accum = AlleleSpreadAccum {
            total_copies: 32,
            distance_weighted: 16,
        };
        let k = 0.5_f64;
        let expected = ((1.0 + k * k).sqrt() - 1.0) / k;
        assert!((fit_g0_decay(&accum, &cfg).p - expected).abs() < 1e-12);
    }

    #[test]
    fn g0_fit_falls_back_for_a_thin_period() {
        let cfg = G0FitCfg::dev_default(); // min_copies = 30
        let accum = AlleleSpreadAccum {
            total_copies: 10,
            distance_weighted: 30,
        };
        assert_eq!(fit_g0_decay(&accum, &cfg).p, cfg.fallback_p);
    }

    #[test]
    fn g0_fit_clamps_against_over_tightening() {
        let cfg = G0FitCfg {
            min_copies: 1,
            fallback_p: 0.5,
            min_p: 0.1,
        };
        // A very tight spread (K̄ ≈ 0.01) would fit p ≈ 0.005; clamped up to min_p.
        let tight = AlleleSpreadAccum {
            total_copies: 100,
            distance_weighted: 1,
        };
        assert_eq!(fit_g0_decay(&tight, &cfg).p, 0.1);
        // Every allele on the mode (K̄ = 0) also clamps rather than collapsing to 0.
        let on_mode = AlleleSpreadAccum {
            total_copies: 100,
            distance_weighted: 0,
        };
        assert_eq!(fit_g0_decay(&on_mode, &cfg).p, 0.1);
    }

    #[test]
    fn g0_spread_counts_only_variable_loci() {
        // A variable cohort (homozygotes spread over 7..=12) populates the spread; a
        // monomorphic cohort (all the same length) contributes nothing.
        let variable = run_prepass_stats(
            &collect_loci(&recovery_spec(0.004, 1.0, 2.0, 0.1, 0.10)),
            2,
            &RungCfg::dev_default(),
        );
        let spread = variable.allele_spread_by_period[&2];
        assert!(
            spread.total_copies > 0,
            "variable loci should populate spread"
        );

        let mono = SimCohortSpec {
            seed: 7,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 60,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 10,
            }],
            groups: vec![SimChemistry {
                error: PerBaseError(0.004),
                shape: StutterShape {
                    up_rate: 1.0,
                    down_rate: 2.0,
                    decay: 0.1,
                },
                level: StutterLevel {
                    baseline: 0.06,
                    slope: 0.0,
                },
            }],
            // All eight samples homozygous for the SAME length → one distinct allele.
            samples: (0..8)
                .map(|i| SimSample {
                    name: format!("S{i}"),
                    group: SampleGroupId(0),
                    genotypes: vec![Some(SimGenotype::homozygous(10, 2))],
                })
                .collect(),
            depth: 200,
        };
        let mono_stats = run_prepass_stats(&collect_loci(&mono), 2, &RungCfg::dev_default());
        assert!(
            !mono_stats.allele_spread_by_period.contains_key(&2),
            "a monomorphic locus must not feed the G₀ spread"
        );
    }

    #[test]
    fn g0_by_period_is_fit_for_a_variable_cohort() {
        let est = run(&recovery_spec(0.004, 1.0, 2.0, 0.1, 0.10));
        let p = est.g0_by_period[&2].p;
        assert!(
            (0.0..=1.0).contains(&p) && p >= G0FitCfg::dev_default().min_p,
            "fitted G₀ p = {p} should be a clamped decay in (0,1]"
        );
    }

    #[test]
    fn g0_spread_counts_het_alleles_as_one_copy_each() {
        // 8 separated 6/10 hets → alleles 6 and 10, ONE chromosome copy each per sample
        // (the 2-peak multiplicity branch): 8 × (1 + 1) = 16 copies. The two alleles sit
        // 0 and 4 units from the read modal length (which is one of 6 / 10), so the
        // distance-weighted sum is 8·0 + 8·4 = 32 either way.
        let stats = run_prepass_stats(&collect_loci(&het_spec()), 2, &RungCfg::dev_default());
        let spread = stats.allele_spread_by_period[&2];
        assert_eq!(
            spread.total_copies, 16,
            "het = one copy per allele per sample"
        );
        assert_eq!(spread.distance_weighted, 32);
    }

    #[test]
    fn g0_spread_measures_distance_from_an_asymmetric_mode() {
        // A skewed cohort: 6 homozygotes at the modal length 10, 3 at 13 (both ≥
        // recurrence_k). The mode is 10 (more reads), so the minority allele lands at
        // distance 3: copies = 6·2 + 3·2 = 18; weighted = 12·0 + 6·3 = 18 ⇒ K̄ = 1.0.
        let chem = SimChemistry {
            error: PerBaseError(0.004),
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
                name: format!("modal{i}"),
                group: SampleGroupId(0),
                genotypes: vec![Some(SimGenotype::homozygous(10, 2))],
            })
            .collect();
        samples.extend((0..3).map(|i| SimSample {
            name: format!("minor{i}"),
            group: SampleGroupId(0),
            genotypes: vec![Some(SimGenotype::homozygous(13, 2))],
        }));
        let spec = SimCohortSpec {
            seed: 31,
            loci: vec![SimLocus {
                chrom: "chr1".into(),
                start: 60,
                motif: Motif::new(b"CA").unwrap(),
                ref_units: 10,
            }],
            groups: vec![chem],
            samples,
            depth: 200,
        };
        let stats = run_prepass_stats(&collect_loci(&spec), 2, &RungCfg::dev_default());
        let spread = stats.allele_spread_by_period[&2];
        assert_eq!(spread.total_copies, 18);
        assert_eq!(spread.distance_weighted, 18);
    }
}
