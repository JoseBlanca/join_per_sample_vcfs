//! The per-locus hidden-paralog scoring function (Q3).
//!
//! [`score_locus_for_paralogy`] asks, at one locus, which of two stories
//! better explains all the samples' coverage and allele counts:
//!
//! - **H1 — a real, single-copy variant.** The allele frequency `p` is
//!   marginalised under the folded site-frequency-spectrum prior
//!   `∝ 1/(p(1−p))` on `[1/2N, 1−1/2N]`; per sample the genotype prior is the
//!   Wright inbreeding-adjusted HWE with the sample's `F`, the allele factor
//!   is binomial over `vaf ∈ {ε, ½, 1−ε}`, and coverage is `~ Normal(1, σ₀)`
//!   *independent of genotype*.
//! - **H2 — a hidden (reference-collapsed) paralog.** Per carrier
//!   configuration `(T, m)` (total copies `T ∈ {3,4,6,8}`, keeping only the
//!   single-PSV `m = 1` and balanced `m ≈ T/2`), a carrier has coverage
//!   `~ Normal(T/2, σ₀·√(T/2))` and `vaf = m/T < 1`; a non-carrier is
//!   single-copy REF. Carrier vs non-carrier follows the Wright dosage HWE in
//!   the carrier frequency `q`, marginalised over `configuration × q`.
//!
//! The result is `paralog_log_likelihood_ratio = logL(H2) − logL(H1)`, a
//! **pure** likelihood ratio (nothing added — MQDiff is shared with
//! introgression and stays a VCF INFO field, spec §2–§3), so the prior/FDR
//! step (Q5) can treat it as a Bayes factor. Both hypotheses are *marginal*
//! (log-sum-exp over their grids, not maximised — averaging charges each
//! hypothesis for its flexibility). Because every H2 carrier has `vaf < 1`, a
//! confidently hom-alt sample is near-impossible under H2 — the hom-alt veto
//! falls out for free.
//!
//! Faithful to the prototype `benchmarks/tomato2/src/build_paralog_lr.py`,
//! with the production divergences settled in arch Premise 2 / spec §5:
//! per-sample σ₀ (from each sample's [`SingleCopyCoverageModel`]), the
//! folded-SFS `p` prior, and MQDiff excluded from the score.
//!
//! [`SingleCopyCoverageModel`]: crate::paralog::SingleCopyCoverageModel

use super::ParalogModelParams;

/// One sample's evidence at a locus — built for *every* locus scored, not
/// only known paralogs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SampleObservation {
    /// Window depth expressed in copies (`1.0 = one copy`; from the
    /// [`SingleCopyCoverageModel`]). Winsorised to
    /// `params.max_relative_copy_number` inside the scorer.
    ///
    /// [`SingleCopyCoverageModel`]: crate::paralog::SingleCopyCoverageModel
    pub relative_copy_number: f64,
    /// Reads supporting the ALT allele(s) (`AD` alt).
    pub alt_reads: u32,
    /// Total reads at the site (`AD` sum).
    pub total_reads: u32,
    /// The sample's inbreeding coefficient `F` (selfer ↔ outbred), `[0, 0.99]`.
    pub inbreeding_coefficient: f64,
}

/// All samples' evidence at one locus. `samples[i] == None` means sample `i`
/// has no usable data here; the slice length is the cohort size `N` (used
/// for the SFS floor `1/2N`).
#[derive(Debug, Clone, Copy)]
pub struct LocusObservations<'a> {
    /// Per-cohort-sample evidence, `None` where unusable.
    pub samples: &'a [Option<SampleObservation>],
}

/// The paralog verdict for one locus.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParalogScore {
    /// `logL(H2) − logL(H1)`: a *pure* likelihood ratio. `> 0` favours a
    /// hidden paralog, `< 0` a real single-copy variant.
    pub paralog_log_likelihood_ratio: f64,
    /// Samples with usable data that entered the score.
    pub samples_used: usize,
    /// Confident hom-alt carriers (`vaf > homalt_vaf_threshold` and
    /// `total_reads >= homalt_min_depth`) — the hom-alt veto signal.
    pub confident_homalt_carriers: usize,
    /// `logL(H1)` — log P(data | real variant). Exposed for the parity check.
    pub log_likelihood_real_variant: f64,
    /// `logL(H2)` — log P(data | hidden paralog). Exposed for the parity check.
    pub log_likelihood_hidden_paralog: f64,
}

impl ParalogScore {
    /// The neutral verdict for a locus with no usable samples: a zero
    /// likelihood ratio (neither hypothesis is favoured).
    fn neutral() -> Self {
        Self {
            paralog_log_likelihood_ratio: 0.0,
            samples_used: 0,
            confident_homalt_carriers: 0,
            log_likelihood_real_variant: 0.0,
            log_likelihood_hidden_paralog: 0.0,
        }
    }
}

/// `0.5·ln(2π)` — the constant term of a Normal log-density.
const LN_SQRT_2PI: f64 = 0.918_938_533_204_672_7;

/// A tiny floor clamping an analytically-positive probability away from
/// `ln(0) = −∞`. The prototype instead *adds* `1e-300`; the clamp is
/// equivalent to f64 precision for any realistic probability and likewise
/// never bites on the default grids.
const PROB_FLOOR: f64 = 1e-300;

/// A usable sample paired with its one-copy relative-depth SD (σ₀). Formed
/// once in [`score_locus_for_paralogy`] after validation (coverage
/// winsorised, `alt_reads ≤ total_reads`, σ₀ finite and `> 0`) so both
/// likelihood functions see a consistent, already-checked pair.
#[derive(Debug, Clone, Copy)]
struct UsableSample {
    observation: SampleObservation,
    single_copy_depth_sd: f64,
}

/// One H2 carrier configuration `(T, m)`: the coverage mean `T/2`, its σ
/// scale `√(T/2)`, and the log allele factors for `vaf = m/T`.
#[derive(Debug, Clone, Copy)]
struct CarrierConfig {
    coverage_mean: f64,
    coverage_sigma_factor: f64,
    log_vaf: f64,
    log1m_vaf: f64,
}

/// Score one locus for paralogy. `single_copy_depth_sd` holds each sample's
/// σ₀ (one-copy relative-depth SD), indexed parallel to
/// `observations.samples`. Pure: same inputs → same output, no state.
pub fn score_locus_for_paralogy(
    observations: &LocusObservations,
    single_copy_depth_sd: &[f64],
    params: &ParalogModelParams,
) -> ParalogScore {
    let cohort_size = observations.samples.len();
    let configs = enumerate_carrier_configs(params);
    // A σ₀ slice out of step with the observations would silently truncate the
    // cohort via the `zip` below (the SFS floor `1/2N` still counts the dropped
    // sample), so a mismatch is a hard precondition failure → neutral, not a
    // plausible-but-wrong score.
    if cohort_size == 0 || configs.is_empty() || single_copy_depth_sd.len() != cohort_size {
        return ParalogScore::neutral();
    }

    // Winsorise coverage and pair each usable sample with its σ₀, dropping
    // samples whose σ₀ is degenerate (an unfit coverage model → `≤ 0` or
    // non-finite would make the Normal term `NaN`).
    let cmax = params.max_relative_copy_number;
    let mut usable: Vec<UsableSample> = Vec::new();
    let mut confident_homalt_carriers = 0usize;
    for (obs, &sigma0) in observations.samples.iter().zip(single_copy_depth_sd.iter()) {
        let Some(mut sample) = *obs else { continue };
        if !(sigma0.is_finite() && sigma0 > 0.0) {
            continue;
        }
        // A malformed record with `alt_reads > total_reads` would underflow the
        // `total_reads − alt_reads` ref count; clamp once so both likelihoods
        // see a coherent pair.
        sample.alt_reads = sample.alt_reads.min(sample.total_reads);
        sample.relative_copy_number = sample.relative_copy_number.clamp(0.0, cmax);
        // `total_reads > 0` guards the `homalt_min_depth == 0` config (with the
        // default `5` it is already implied).
        if sample.total_reads >= params.homalt_min_depth
            && sample.total_reads > 0
            && f64::from(sample.alt_reads) / f64::from(sample.total_reads)
                > params.homalt_vaf_threshold
        {
            confident_homalt_carriers += 1;
        }
        usable.push(UsableSample {
            observation: sample,
            single_copy_depth_sd: sigma0,
        });
    }
    if usable.is_empty() {
        return ParalogScore::neutral();
    }

    let log_likelihood_real_variant = h1_log_likelihood(&usable, cohort_size, params);
    let log_likelihood_hidden_paralog = h2_log_likelihood(&usable, &configs, params);

    ParalogScore {
        paralog_log_likelihood_ratio: log_likelihood_hidden_paralog - log_likelihood_real_variant,
        samples_used: usable.len(),
        confident_homalt_carriers,
        log_likelihood_real_variant,
        log_likelihood_hidden_paralog,
    }
}

/// H1 — real, single-copy variant. Coverage is `Normal(1, σ₀)` independent of
/// genotype (a constant added once); the allele/genotype term is marginalised
/// over `p` under the folded-SFS prior.
fn h1_log_likelihood(
    usable: &[UsableSample],
    cohort_size: usize,
    params: &ParalogModelParams,
) -> f64 {
    let eps = params.pseudocount_vaf;
    // H1 genotype VAFs: hom-ref ε, het ½, hom-alt 1−ε.
    let log_vaf = [eps.ln(), 0.5f64.ln(), (1.0 - eps).ln()];
    let log1m_vaf = [(1.0 - eps).ln(), 0.5f64.ln(), eps.ln()];

    // Per-sample allele factor per genotype, and the genotype-independent
    // coverage term summed once.
    let mut coverage_term = 0.0;
    let allele_by_genotype: Vec<[f64; 3]> = usable
        .iter()
        .map(|sample| {
            let obs = &sample.observation;
            coverage_term += ln_normal(obs.relative_copy_number, 1.0, sample.single_copy_depth_sd);
            let k = f64::from(obs.alt_reads);
            let nk = f64::from(obs.total_reads - obs.alt_reads);
            [
                k * log_vaf[0] + nk * log1m_vaf[0],
                k * log_vaf[1] + nk * log1m_vaf[1],
                k * log_vaf[2] + nk * log1m_vaf[2],
            ]
        })
        .collect();

    // Marginalise p under the folded SFS `1/(p(1−p))` on `[1/2N, 1−1/2N]`:
    //   logL1 = LSE_p( summed1[p] + logw_p ) + coverage,
    // with `logw_p` the normalised prior. We accumulate the numerator LSE
    // (`summed1[p] + raw_logw_p`) and the normaliser LSE (`raw_logw_p`)
    // online, so `LSE_num − LSE_norm` is the p-marginal without a p vector.
    let n_points = params.allele_freq_prior.n_points;
    let inv2n = 1.0 / (2.0 * cohort_size as f64);
    let mut lse_num = LogSumExp::new();
    let mut lse_norm = LogSumExp::new();
    for i in 0..n_points {
        let p = sfs_grid_point(i, n_points, inv2n);
        let raw_logw = -(p * (1.0 - p)).ln();
        let mut summed = 0.0;
        for (sample, allele) in usable.iter().zip(allele_by_genotype.iter()) {
            let (log_homref, log_het, log_homalt) =
                wright_genotype_log_priors(p, sample.observation.inbreeding_coefficient);
            summed += log_sum_exp3(
                log_homref + allele[0],
                log_het + allele[1],
                log_homalt + allele[2],
            );
        }
        lse_num.push(summed + raw_logw);
        lse_norm.push(raw_logw);
    }

    lse_num.value() - lse_norm.value() + coverage_term
}

/// H2 — hidden paralog. Per `(configuration, q)`, each sample is a mixture of
/// non-carrier (single-copy REF) and carrier `(T, m)`; marginalised flat over
/// `configuration × q`.
fn h2_log_likelihood(
    usable: &[UsableSample],
    configs: &[CarrierConfig],
    params: &ParalogModelParams,
) -> f64 {
    let eps = params.pseudocount_vaf;
    let (log_eps, log1m_eps) = (eps.ln(), (1.0 - eps).ln());
    let n_configs = configs.len();

    // Per-sample non-carrier base term (coverage 1, VAF ε) and the carrier
    // branch per configuration (coverage T/2, VAF m/T). `carrier_branch` is a
    // row-major `[sample][config]` table with stride `n_configs`.
    let mut noncarrier_base = Vec::with_capacity(usable.len());
    let mut carrier_branch = vec![0.0f64; usable.len() * n_configs];
    for (s, sample) in usable.iter().enumerate() {
        let obs = &sample.observation;
        let sigma0 = sample.single_copy_depth_sd;
        let c = obs.relative_copy_number;
        let k = f64::from(obs.alt_reads);
        let nk = f64::from(obs.total_reads - obs.alt_reads);
        noncarrier_base.push(ln_normal(c, 1.0, sigma0) + k * log_eps + nk * log1m_eps);
        for (j, config) in configs.iter().enumerate() {
            let cov = ln_normal(
                c,
                config.coverage_mean,
                sigma0 * config.coverage_sigma_factor,
            );
            let allele = k * config.log_vaf + nk * config.log1m_vaf;
            carrier_branch[s * n_configs + j] = cov + allele;
        }
    }

    // Marginalise over configuration × q (flat prior): logL2 = LSE(g) − ln|g|,
    // g[config,q] = Σ_s logaddexp(P(noncarrier|q)+base, P(carrier|q)+carrier).
    let grid = params.carrier_freq_grid;
    let mut log_noncarrier = vec![0.0f64; usable.len()];
    let mut log_carrier = vec![0.0f64; usable.len()];
    let mut lse = LogSumExp::new();
    let mut cells = 0usize;
    for iq in 0..grid.n_points {
        let q = linspace_point(iq, grid.n_points, grid.lo, grid.hi);
        // Wright dosage HWE: P(non-carrier) = (1−q)² + F·q(1−q). The floor
        // guards a non-default grid reaching `q → 1` (where P0 → 0); it is
        // inert on the default grid where `min P0 = 0.16` (prototype floors
        // only the carrier branch, but this is harmless there).
        for (s, sample) in usable.iter().enumerate() {
            let f = sample.observation.inbreeding_coefficient;
            let p_noncarrier = ((1.0 - q) * (1.0 - q) + f * q * (1.0 - q)).max(PROB_FLOOR);
            log_noncarrier[s] = p_noncarrier.ln();
            log_carrier[s] = (1.0 - p_noncarrier).max(PROB_FLOOR).ln();
        }
        for j in 0..n_configs {
            let mut acc = 0.0;
            for s in 0..usable.len() {
                let a = log_noncarrier[s] + noncarrier_base[s];
                let b = log_carrier[s] + carrier_branch[s * n_configs + j];
                acc += log_add_exp(a, b);
            }
            lse.push(acc);
            cells += 1;
        }
    }

    lse.value() - (cells as f64).ln()
}

/// Enumerate the kept H2 carrier configurations: for each `T`, the single-PSV
/// (`m = 1`) and balanced (`m ≈ T/2`) copy configurations with `1 ≤ m ≤ T−2`.
fn enumerate_carrier_configs(params: &ParalogModelParams) -> Vec<CarrierConfig> {
    let mut out = Vec::new();
    for &t in &params.carrier_copy_numbers {
        if t < 3 {
            continue; // needs 1 ≤ m ≤ T−2, so T ≥ 3
        }
        let coverage_mean = f64::from(t) / 2.0;
        let coverage_sigma_factor = coverage_mean.sqrt();
        let mut ms = [1u32, t / 2];
        ms.sort_unstable();
        let mut last = None;
        for &m in &ms {
            if Some(m) == last || m < 1 || m > t - 2 {
                continue;
            }
            last = Some(m);
            let vaf = f64::from(m) / f64::from(t);
            out.push(CarrierConfig {
                coverage_mean,
                coverage_sigma_factor,
                log_vaf: vaf.ln(),
                log1m_vaf: (1.0 - vaf).ln(),
            });
        }
    }
    out
}

/// The `i`-th of `n` folded-SFS grid points on `[inv2n, 1−inv2n]`, uniform in
/// `p`. A single-point grid collapses to the midpoint `0.5`.
fn sfs_grid_point(i: usize, n: usize, inv2n: f64) -> f64 {
    linspace_point(i, n, inv2n, 1.0 - inv2n)
}

/// The `i`-th of `n` points on the inclusive linear grid `[lo, hi]`. `n <= 1`
/// yields the midpoint (a degenerate grid still returns a usable point).
fn linspace_point(i: usize, n: usize, lo: f64, hi: f64) -> f64 {
    if n <= 1 {
        return 0.5 * (lo + hi);
    }
    lo + (hi - lo) * (i as f64) / ((n - 1) as f64)
}

/// The Wright inbreeding-adjusted HWE genotype log-priors `(hom-ref, het,
/// hom-alt)` at ALT frequency `p` and inbreeding `F`:
/// `P(het) = 2pq(1−F)`, homozygotes `q²+Fpq` / `p²+Fpq` (`q = 1−p`).
fn wright_genotype_log_priors(p: f64, f: f64) -> (f64, f64, f64) {
    let q = 1.0 - p;
    let het = 2.0 * p * q * (1.0 - f);
    let hom_ref = q * q + f * p * q;
    let hom_alt = p * p + f * p * q;
    (
        hom_ref.max(PROB_FLOOR).ln(),
        het.max(PROB_FLOOR).ln(),
        hom_alt.max(PROB_FLOOR).ln(),
    )
}

/// A Normal log-density `ln N(x; μ, σ)`. `σ` must be `> 0`.
fn ln_normal(x: f64, mu: f64, sigma: f64) -> f64 {
    let z = (x - mu) / sigma;
    -0.5 * z * z - sigma.ln() - LN_SQRT_2PI
}

/// `ln(exp(a) + exp(b))`, stable and `−∞`-safe.
fn log_add_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    hi + (lo - hi).exp().ln_1p()
}

/// `ln(exp(a) + exp(b) + exp(c))`, stable.
fn log_sum_exp3(a: f64, b: f64, c: f64) -> f64 {
    log_add_exp(log_add_exp(a, b), c)
}

/// Online log-sum-exp accumulator: fold values in one at a time without
/// materialising the vector, tracking the running max and rescaled sum.
#[derive(Debug, Clone, Copy)]
struct LogSumExp {
    max: f64,
    sum: f64,
}

impl LogSumExp {
    fn new() -> Self {
        Self {
            max: f64::NEG_INFINITY,
            sum: 0.0,
        }
    }

    fn push(&mut self, x: f64) {
        if x == f64::NEG_INFINITY {
            return;
        }
        if x > self.max {
            self.sum = self.sum * (self.max - x).exp() + 1.0;
            self.max = x;
        } else {
            self.sum += (x - self.max).exp();
        }
    }

    fn value(&self) -> f64 {
        if self.max == f64::NEG_INFINITY {
            f64::NEG_INFINITY
        } else {
            self.max + self.sum.ln()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SIGMA0: f64 = 0.26;

    fn obs(rel: f64, alt: u32, total: u32, f: f64) -> Option<SampleObservation> {
        Some(SampleObservation {
            relative_copy_number: rel,
            alt_reads: alt,
            total_reads: total,
            inbreeding_coefficient: f,
        })
    }

    fn score(samples: &[Option<SampleObservation>]) -> ParalogScore {
        let sds = vec![SIGMA0; samples.len()];
        score_locus_for_paralogy(
            &LocusObservations { samples },
            &sds,
            &ParalogModelParams::default(),
        )
    }

    /// A cohort of confident hom-alt samples at normal coverage: H2 can never
    /// produce `VAF = 1` (every carrier config has `vaf < 1`), so the pure LR
    /// is driven strongly negative — the hom-alt veto — and the diagnostic
    /// counter records every carrier.
    #[test]
    fn confident_homalt_forces_negative_lr() {
        let samples: Vec<_> = (0..30).map(|_| obs(1.0, 20, 20, 0.9)).collect();
        let s = score(&samples);
        assert!(
            s.paralog_log_likelihood_ratio < 0.0,
            "LR = {}",
            s.paralog_log_likelihood_ratio
        );
        assert_eq!(s.confident_homalt_carriers, 30);
        assert_eq!(s.samples_used, 30);
    }

    /// Excess-coverage carriers (2× depth, balanced VAF) mixed with
    /// single-copy REF non-carriers: H1 must pay the Normal(1, σ₀) coverage
    /// penalty for every 2× window, which H2 explains for free — so the LR is
    /// positive. Coverage is the introgression-safe primary signal.
    #[test]
    fn excess_coverage_carriers_force_positive_lr() {
        let mut samples: Vec<_> = (0..15).map(|_| obs(2.0, 10, 20, 0.0)).collect();
        samples.extend((0..15).map(|_| obs(1.0, 0, 20, 0.0)));
        let s = score(&samples);
        assert!(
            s.paralog_log_likelihood_ratio > 0.0,
            "LR = {}",
            s.paralog_log_likelihood_ratio
        );
        assert_eq!(s.confident_homalt_carriers, 0);
    }

    /// Excess heterozygosity at *normal* coverage (an introgression profile):
    /// H2 would need 2× coverage to call the hets carriers, which the depth
    /// contradicts, so the LR stays negative regardless of the allele
    /// pattern — the model does not flag introgressions.
    #[test]
    fn normal_coverage_excess_het_stays_negative() {
        let mut samples: Vec<_> = (0..15).map(|_| obs(1.0, 10, 20, 0.0)).collect();
        samples.extend((0..15).map(|_| obs(1.0, 0, 20, 0.0)));
        let s = score(&samples);
        assert!(
            s.paralog_log_likelihood_ratio < 0.0,
            "LR = {}",
            s.paralog_log_likelihood_ratio
        );
    }

    /// Coverage is the primary signal: the same allele pattern at 2× depth
    /// scores strictly higher than at 1× depth.
    #[test]
    fn higher_coverage_raises_lr() {
        let het_1x: Vec<_> = (0..30).map(|_| obs(1.0, 10, 20, 0.0)).collect();
        let het_2x: Vec<_> = (0..30).map(|_| obs(2.0, 10, 20, 0.0)).collect();
        let lr_1x = score(&het_1x).paralog_log_likelihood_ratio;
        let lr_2x = score(&het_2x).paralog_log_likelihood_ratio;
        assert!(lr_2x > lr_1x, "1x LR {lr_1x} should be < 2x LR {lr_2x}");
    }

    /// The reported LR is exactly `logL(H2) − logL(H1)` and both diagnostic
    /// likelihoods are finite.
    #[test]
    fn lr_equals_difference_of_diagnostics() {
        let samples: Vec<_> = (0..25)
            .map(|i| obs(1.0 + (i % 2) as f64, 5, 20, 0.1))
            .collect();
        let s = score(&samples);
        assert!(s.log_likelihood_real_variant.is_finite());
        assert!(s.log_likelihood_hidden_paralog.is_finite());
        let diff = s.log_likelihood_hidden_paralog - s.log_likelihood_real_variant;
        assert!((s.paralog_log_likelihood_ratio - diff).abs() < 1e-12);
    }

    /// The scoring function is deterministic — identical inputs give bit-equal
    /// outputs (a prerequisite for the reproducible calibration pass).
    #[test]
    fn scoring_is_deterministic() {
        let samples: Vec<_> = (0..20).map(|i| obs(1.0, i, 20, 0.2)).collect();
        assert_eq!(score(&samples), score(&samples));
    }

    /// A locus with no usable samples (empty, or all `None`) returns the
    /// neutral verdict rather than a `NaN`.
    #[test]
    fn no_usable_samples_is_neutral() {
        assert_eq!(score(&[]).paralog_log_likelihood_ratio, 0.0);
        let none: Vec<Option<SampleObservation>> = vec![None; 5];
        let s = score(&none);
        assert_eq!(s.paralog_log_likelihood_ratio, 0.0);
        assert_eq!(s.samples_used, 0);
    }

    /// `None` samples are skipped and `samples_used` counts only usable ones;
    /// the SFS floor `1/2N` still uses the full cohort length.
    #[test]
    fn none_samples_are_skipped() {
        let samples = vec![obs(1.0, 10, 20, 0.0), None, obs(2.0, 5, 20, 0.0), None];
        let s = score(&samples);
        assert_eq!(s.samples_used, 2);
    }

    /// The hom-alt veto counter respects both thresholds: a high VAF below
    /// `homalt_min_depth`, or a moderate VAF at high depth, does not count.
    #[test]
    fn confident_homalt_counter_respects_thresholds() {
        let samples = vec![
            obs(1.0, 4, 4, 0.0),   // VAF 1.0 but depth 4 < 5 → not counted
            obs(1.0, 19, 20, 0.0), // VAF 0.95 ≥ 0.9, depth 20 → counted
            obs(1.0, 15, 20, 0.0), // VAF 0.75 < 0.9 → not counted
            obs(1.0, 20, 20, 0.0), // VAF 1.0, depth 20 → counted
        ];
        let s = score(&samples);
        assert_eq!(s.confident_homalt_carriers, 2);
    }

    /// The default carrier configurations are the single-PSV + balanced set
    /// `(3,1),(4,1),(4,2),(6,1),(6,3),(8,1),(8,4)` = 7 configs, with the
    /// expected `T/2` means and `m/T` VAFs.
    #[test]
    fn default_carrier_configs_are_the_seven_expected() {
        let configs = enumerate_carrier_configs(&ParalogModelParams::default());
        assert_eq!(configs.len(), 7);
        let means: Vec<f64> = configs.iter().map(|c| c.coverage_mean).collect();
        assert_eq!(means, vec![1.5, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0]);
        // VAF = exp(log_vaf); check the balanced m≈T/2 configs are ½.
        let vafs: Vec<f64> = configs.iter().map(|c| c.log_vaf.exp()).collect();
        assert!((vafs[1] - 0.25).abs() < 1e-12); // (4,1)
        assert!((vafs[2] - 0.5).abs() < 1e-12); // (4,2)
        assert!((vafs[6] - 0.5).abs() < 1e-12); // (8,4)
    }

    /// `ln N(x; μ, σ)` matches the closed form at the mean and one σ out.
    #[test]
    fn ln_normal_matches_closed_form() {
        // At the mean: −ln(σ√2π).
        assert!(
            (ln_normal(1.0, 1.0, 0.5) - (-(0.5f64 * (2.0 * std::f64::consts::PI).sqrt()).ln()))
                .abs()
                < 1e-12
        );
        // One σ out drops by exactly ½.
        let d = ln_normal(1.0, 1.0, 0.5) - ln_normal(1.5, 1.0, 0.5);
        assert!((d - 0.5).abs() < 1e-12);
    }

    /// `log_add_exp` is `−∞`-safe and matches the naive form where stable.
    #[test]
    fn log_add_exp_is_stable_and_neg_inf_safe() {
        assert_eq!(log_add_exp(f64::NEG_INFINITY, -3.0), -3.0);
        assert_eq!(log_add_exp(-3.0, f64::NEG_INFINITY), -3.0);
        let naive = (1.0f64.exp() + 2.0f64.exp()).ln();
        assert!((log_add_exp(1.0, 2.0) - naive).abs() < 1e-12);
    }

    /// The online `LogSumExp` matches a batched log-sum-exp and returns `−∞`
    /// for an empty accumulation.
    #[test]
    fn log_sum_exp_accumulator_matches_batch() {
        let xs = [-2.0, 0.5, 3.0, -10.0];
        let mut lse = LogSumExp::new();
        for &x in &xs {
            lse.push(x);
        }
        let max = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let batch = max + xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
        assert!((lse.value() - batch).abs() < 1e-12);
        assert_eq!(LogSumExp::new().value(), f64::NEG_INFINITY);
    }

    /// A malformed record with `alt_reads > total_reads` is clamped at the
    /// source (no `u32` underflow), so the LR stays finite.
    #[test]
    fn score_stays_finite_when_alt_exceeds_total() {
        let samples: Vec<_> = (0..25).map(|_| obs(1.0, 30, 20, 0.0)).collect();
        let s = score(&samples);
        assert!(
            s.paralog_log_likelihood_ratio.is_finite(),
            "LR = {}",
            s.paralog_log_likelihood_ratio
        );
    }

    /// A degenerate σ₀ (`0.0` / non-finite) from an unfit sample is dropped
    /// rather than propagating a `NaN` through the Normal coverage term.
    #[test]
    fn degenerate_sigma0_sample_is_dropped_not_nan() {
        let samples: Vec<_> = (0..25).map(|_| obs(1.0, 10, 20, 0.0)).collect();
        let mut sds = vec![SIGMA0; samples.len()];
        sds[0] = 0.0;
        sds[1] = f64::NAN;
        let s = score_locus_for_paralogy(
            &LocusObservations { samples: &samples },
            &sds,
            &ParalogModelParams::default(),
        );
        assert!(s.paralog_log_likelihood_ratio.is_finite());
        assert_eq!(s.samples_used, 23); // two degenerate σ₀ samples dropped
    }

    /// A σ₀ slice out of step with the observations is a precondition failure:
    /// the score is neutral, not a silently-truncated cohort.
    #[test]
    fn mismatched_sigma0_length_is_neutral() {
        let samples: Vec<_> = (0..25).map(|_| obs(1.0, 10, 20, 0.0)).collect();
        let sds = vec![SIGMA0; samples.len() - 1]; // one short
        let s = score_locus_for_paralogy(
            &LocusObservations { samples: &samples },
            &sds,
            &ParalogModelParams::default(),
        );
        assert_eq!(s.paralog_log_likelihood_ratio, 0.0);
        assert_eq!(s.samples_used, 0);
    }

    /// Absolute-value numerical parity anchor: a single sample on collapsed
    /// grids (one SFS point → p = ½, one q point → q = ½, T = 4) makes both
    /// marginals hand-computable. Pins the *magnitude* of logL1/logL2, which
    /// the sign-only tests cannot — a constant offset in either would slip by.
    #[test]
    fn reduced_grid_matches_hand_computed_likelihoods() {
        use crate::paralog::{GridSpec, SfsPriorSpec};
        let params = ParalogModelParams {
            allele_freq_prior: SfsPriorSpec { n_points: 1 },
            carrier_freq_grid: GridSpec {
                lo: 0.5,
                hi: 0.5,
                n_points: 1,
            },
            carrier_copy_numbers: vec![4],
            ..Default::default()
        };
        // One sample: rel 1.0, 5 alt of 10, F = 0, σ₀ = 0.26.
        let samples = [obs(1.0, 5, 10, 0.0)];
        let s =
            score_locus_for_paralogy(&LocusObservations { samples: &samples }, &[SIGMA0], &params);
        // Hand-derived from ln N + binomial-allele + Wright priors (see the
        // module's H1/H2 formulas): logL1 ≈ −7.1965, logL2 ≈ −11.3157.
        assert!(
            (s.log_likelihood_real_variant - (-7.196484)).abs() < 2e-3,
            "logL1 = {}",
            s.log_likelihood_real_variant
        );
        assert!(
            (s.log_likelihood_hidden_paralog - (-11.315711)).abs() < 2e-3,
            "logL2 = {}",
            s.log_likelihood_hidden_paralog
        );
        assert!(
            (s.paralog_log_likelihood_ratio - (-4.119227)).abs() < 2e-3,
            "LR = {}",
            s.paralog_log_likelihood_ratio
        );
    }

    /// Config enumeration on odd and sub-3 `T`: `T = 2` is skipped (needs
    /// `1 ≤ m ≤ T−2`), `T = 5 → m ∈ {1,2}`, `T = 7 → m ∈ {1,3}`.
    #[test]
    fn enumerate_carrier_configs_handles_odd_and_small_t() {
        let params = ParalogModelParams {
            carrier_copy_numbers: vec![2, 5, 7],
            ..Default::default()
        };
        let configs = enumerate_carrier_configs(&params);
        let vafs: Vec<f64> = configs.iter().map(|c| c.log_vaf.exp()).collect();
        assert_eq!(configs.len(), 4);
        assert!((vafs[0] - 1.0 / 5.0).abs() < 1e-12);
        assert!((vafs[1] - 2.0 / 5.0).abs() < 1e-12);
        assert!((vafs[2] - 1.0 / 7.0).abs() < 1e-12);
        assert!((vafs[3] - 3.0 / 7.0).abs() < 1e-12);
    }

    /// A degenerate (`n <= 1`) grid returns the interval midpoint rather than
    /// dividing by `n − 1 = 0`.
    #[test]
    fn linspace_point_degenerate_grid_returns_midpoint() {
        assert_eq!(linspace_point(0, 1, 0.2, 0.8), 0.5);
        assert_eq!(linspace_point(0, 0, 0.2, 0.8), 0.5);
        assert_eq!(sfs_grid_point(0, 1, 0.01), 0.5);
    }

    /// The Wright genotype priors are a proper distribution: `hom-ref + het +
    /// hom-alt = 1` for every `(p, F)` — a coefficient typo would break this.
    #[test]
    fn wright_genotype_priors_sum_to_one() {
        for &p in &[0.01, 0.2, 0.5, 0.9] {
            for &f in &[0.0, 0.3, 0.99] {
                let (a, b, c) = wright_genotype_log_priors(p, f);
                let sum = a.exp() + b.exp() + c.exp();
                assert!((sum - 1.0).abs() < 1e-12, "p={p} F={f} sum={sum}");
            }
        }
    }
}
