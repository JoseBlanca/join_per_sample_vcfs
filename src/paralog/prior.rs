//! The empirical-Bayes prior + false-discovery cut (Q5).
//!
//! The per-locus score ([`super::score_locus_for_paralogy`]) is a pure
//! likelihood ratio — which explanation fits better, but not how *likely* the
//! locus is to be a paralog. That needs the genome-wide **prior probability**
//! that a locus is a hidden duplication (`π`), estimated from the data itself
//! by EM, after which each locus's LR becomes a posterior and a
//! false-discovery q-value:
//!
//! ```text
//! posterior(LR) = σ(LR + logit π)          (a Bayes-factor update)
//! q(LR)         = E[ 1 − posterior | LR' ≥ LR ]   (tail FDR)
//! ```
//!
//! Faithful to the prototype `benchmarks/tomato2/src/build_paralog_eb.py`
//! (`π` by EM `r_i = σ(LR_i + logit π); π = mean(r_i)`; q by the running mean
//! of the local FDR over loci ranked by posterior). The production divergence
//! (arch Premise 6): both `π` and the `q(LR)` curve are derived from a
//! **fixed-size LR histogram**, never a genome-wide per-locus vector, so RAM
//! stays flat in variant count. Because the posterior is strictly monotone in
//! the LR, ranking by posterior is ranking by LR, so the histogram (bins
//! ordered by LR) reproduces the tail FDR exactly up to bin resolution.

/// Default EM start for `π` (paralogs are rare; the fixed point is unique so
/// the start only affects iteration count — prototype `pi = 0.03`).
pub const DEFAULT_EM_START: f64 = 0.03;
/// Default EM convergence tolerance on `|Δπ|` (prototype `1e-9`).
pub const DEFAULT_EM_TOL: f64 = 1e-9;
/// Default EM iteration cap (prototype `500`).
pub const DEFAULT_EM_MAX_ITER: usize = 500;

/// Default lower edge of the LR histogram. LRs below it saturate into the
/// bottom bin, where the posterior is ≈ 0 regardless — so the exact edge does
/// not matter as long as it is well below the sigmoid transition.
pub const DEFAULT_LR_HISTOGRAM_LO: f64 = -100.0;
/// Default upper edge of the LR histogram (LRs above saturate to ≈ posterior 1).
pub const DEFAULT_LR_HISTOGRAM_HI: f64 = 100.0;
/// Default LR-histogram bin count (width `0.1` over `[-100, 100]`): fine
/// enough around the sigmoid transition (`LR ≈ −logit π`) to pin `π` and the
/// FDR cut, ~16 KB of counts.
pub const DEFAULT_LR_HISTOGRAM_BINS: usize = 2000;

/// EM configuration for [`ParalogPrior::estimate`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EmConfig {
    /// Initial `π`.
    pub start: f64,
    /// Stop when `|π_new − π| < tol`.
    pub tol: f64,
    /// Iteration cap.
    pub max_iter: usize,
}

impl Default for EmConfig {
    fn default() -> Self {
        Self {
            start: DEFAULT_EM_START,
            tol: DEFAULT_EM_TOL,
            max_iter: DEFAULT_EM_MAX_ITER,
        }
    }
}

/// A uniform binning of the LR axis into `n_bins` bins over `[lo, hi]`,
/// saturating out-of-range values into the end bins. It is the single
/// definition of the bin geometry that both [`ParalogLrHistogram`] and
/// [`ParalogFdrCurve`] share, so the "the curve's bins are the histogram's
/// bins" invariant (load-bearing for the FDR being faithful) is held in one
/// place rather than by two hand-synced copies.
#[derive(Debug, Clone, Copy, PartialEq)]
struct UniformBinning {
    lo: f64,
    hi: f64,
    n_bins: usize,
}

impl UniformBinning {
    /// Returns `None` for a degenerate binning (`n_bins == 0`, a non-finite
    /// edge, or `lo >= hi`) — mirroring [`super::GridSpec::new`]'s
    /// construction-boundary convention.
    fn new(lo: f64, hi: f64, n_bins: usize) -> Option<Self> {
        (n_bins >= 1 && lo.is_finite() && hi.is_finite() && lo < hi).then_some(Self {
            lo,
            hi,
            n_bins,
        })
    }

    /// Bin index for `lr`, saturating to `[0, n_bins − 1]`. A non-finite `lr`
    /// must be screened by the caller (`NaN` would not clamp).
    fn bin_index(&self, lr: f64) -> usize {
        let t = (lr - self.lo) / (self.hi - self.lo);
        let scaled = t * self.n_bins as f64;
        (scaled.floor() as isize).clamp(0, self.n_bins as isize - 1) as usize
    }

    /// The representative LR of bin `i` (its centre).
    fn bin_center(&self, i: usize) -> f64 {
        self.lo + (i as f64 + 0.5) * (self.hi - self.lo) / self.n_bins as f64
    }
}

/// A fixed-size histogram of per-locus likelihood ratios: the bounded-RAM
/// sufficient statistic from which `π` and the FDR curve are derived. LRs
/// outside `[lo, hi]` saturate into the end bins (the posterior saturates
/// there anyway).
#[derive(Debug, Clone, PartialEq)]
pub struct ParalogLrHistogram {
    binning: UniformBinning,
    counts: Vec<u64>,
    total: u64,
}

impl ParalogLrHistogram {
    /// A new empty histogram of `n_bins` uniform bins over `[lo, hi]`.
    /// Returns `None` for a degenerate range (`n_bins == 0`, a non-finite
    /// edge, or `lo >= hi`) — the same construction-boundary convention as
    /// the sibling [`super::GridSpec::new`].
    pub fn new(lo: f64, hi: f64, n_bins: usize) -> Option<Self> {
        UniformBinning::new(lo, hi, n_bins).map(|binning| Self {
            binning,
            counts: vec![0; n_bins],
            total: 0,
        })
    }

    /// A new empty histogram with the default range and resolution (the
    /// trusted compile-time defaults, so construction cannot fail).
    pub fn with_defaults() -> Self {
        Self::new(
            DEFAULT_LR_HISTOGRAM_LO,
            DEFAULT_LR_HISTOGRAM_HI,
            DEFAULT_LR_HISTOGRAM_BINS,
        )
        .expect("default LR-histogram range is valid")
    }

    /// Fold one likelihood ratio into its bin (saturating at the ends). A
    /// non-finite LR is dropped (it cannot be a valid Bayes factor).
    pub fn push(&mut self, lr: f64) {
        if !lr.is_finite() {
            return;
        }
        let bin = self.binning.bin_index(lr);
        self.counts[bin] += 1;
        self.total += 1;
    }

    /// Total LRs folded in.
    pub fn total(&self) -> u64 {
        self.total
    }

    /// Whether no LRs have been folded in.
    pub fn is_empty(&self) -> bool {
        self.total == 0
    }
}

/// The genome-wide prior probability that a locus is a hidden duplication,
/// estimated from the data (the `π` of empirical Bayes). Applying it to a
/// locus's LR yields that locus's posterior.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParalogPrior {
    /// Prior `P(locus is a hidden duplication)`.
    pub prior_probability: f64,
    /// Whether the EM reached the tolerance within `max_iter` (and the
    /// histogram was non-empty). `false` marks a seed-valued / unconverged
    /// `prior_probability` the caller should treat with suspicion — it is
    /// not distinguishable from a real estimate by the value alone.
    pub converged: bool,
}

impl ParalogPrior {
    /// Estimate `π` from the LR histogram by EM:
    /// `π ← Σ_bin (count_bin / total) · σ(LR_bin + logit π)`, iterated to
    /// `|Δπ| < tol` (the fixed point is unique, so the start only sets the
    /// iteration count). An empty histogram returns `π = cfg.start` with
    /// `converged = false`; hitting `max_iter` returns the last iterate, also
    /// with `converged = false`.
    pub fn estimate(histogram: &ParalogLrHistogram, cfg: &EmConfig) -> Self {
        if histogram.is_empty() {
            return Self {
                prior_probability: cfg.start,
                converged: false,
            };
        }
        let total = histogram.total as f64;
        let mut pi = cfg.start;
        let mut converged = false;
        for _ in 0..cfg.max_iter {
            let offset = logit(clamp_probability(pi));
            let mut acc = 0.0;
            for (i, &count) in histogram.counts.iter().enumerate() {
                if count == 0 {
                    continue;
                }
                acc += count as f64 * sigmoid(histogram.binning.bin_center(i) + offset);
            }
            let pi_new = acc / total;
            let done = (pi_new - pi).abs() < cfg.tol;
            pi = pi_new;
            if done {
                converged = true;
                break;
            }
        }
        Self {
            prior_probability: pi,
            converged,
        }
    }

    /// Posterior `P(paralog | data) = σ(LR + logit π)`.
    pub fn paralog_posterior(&self, lr: f64) -> f64 {
        sigmoid(lr + logit(clamp_probability(self.prior_probability)))
    }

    /// The LR at which the posterior crosses `0.5`:
    /// `log((1 − π) / π) = −logit π` — positive, since hidden duplications are
    /// rare (the cut is NOT at `LR = 0`).
    pub fn half_posterior_ratio(&self) -> f64 {
        -logit(clamp_probability(self.prior_probability))
    }
}

/// The monotone tail-FDR-as-a-function-of-LR curve, precomputed per histogram
/// bin so a locus's q-value is a lookup `q_of_lr(LR_i)` — no genome-wide sort.
#[derive(Debug, Clone, PartialEq)]
pub struct ParalogFdrCurve {
    /// The bin geometry, copied from the source histogram so `q_of_lr` indexes
    /// against the same binning the tail FDR was computed over.
    binning: UniformBinning,
    /// Per-bin q-value: `q_of_bin[i] = E[1 − posterior | LR' ≥ bin_i]`.
    q_of_bin: Vec<f64>,
}

impl ParalogFdrCurve {
    /// Build the tail-FDR curve from the histogram and estimated prior. The
    /// tail is accumulated from the highest-LR bin down; because the posterior
    /// is monotone in LR, the running mean of the local FDR `1 − posterior` is
    /// monotone, so the curve is monotone decreasing in LR by construction.
    pub fn from_histogram(histogram: &ParalogLrHistogram, prior: &ParalogPrior) -> Self {
        let n = histogram.binning.n_bins;
        let mut q_of_bin = vec![0.0; n];
        let mut tail_count = 0u64;
        let mut tail_lfdr_sum = 0.0;
        for i in (0..n).rev() {
            let local_lfdr = 1.0 - prior.paralog_posterior(histogram.binning.bin_center(i));
            tail_count += histogram.counts[i];
            tail_lfdr_sum += histogram.counts[i] as f64 * local_lfdr;
            // Empty tail (above all data): the limiting FDR is the local one.
            q_of_bin[i] = if tail_count > 0 {
                tail_lfdr_sum / tail_count as f64
            } else {
                local_lfdr
            };
        }
        Self {
            binning: histogram.binning,
            q_of_bin,
        }
    }

    /// The tail-FDR q-value at likelihood ratio `lr` (its bin's value).
    pub fn q_of_lr(&self, lr: f64) -> f64 {
        self.q_of_bin[self.binning.bin_index(lr)]
    }

    /// The least-stringent LR threshold whose q-value is `≤ target_fdr`, i.e.
    /// the smallest LR at which flagging everything above it holds the FDR.
    /// `None` if no bin achieves the target. The returned LR is the crossing
    /// bin's centre.
    pub fn lr_threshold_for_fdr(&self, target_fdr: f64) -> Option<f64> {
        // q is monotone decreasing in LR, so the first bin (lowest LR) at or
        // below the target is the least-stringent cut.
        (0..self.q_of_bin.len())
            .find(|&i| self.q_of_bin[i] <= target_fdr)
            .map(|i| self.binning.bin_center(i))
    }
}

/// The logistic function `σ(x) = 1 / (1 + e^{−x})`, `−∞`/`+∞`-safe.
fn sigmoid(x: f64) -> f64 {
    if x >= 0.0 {
        1.0 / (1.0 + (-x).exp())
    } else {
        let e = x.exp();
        e / (1.0 + e)
    }
}

/// The logit `ln(p / (1 − p))`. Callers clamp `p` into the open unit interval
/// first via [`clamp_probability`].
fn logit(p: f64) -> f64 {
    (p / (1.0 - p)).ln()
}

/// Clamp a probability into the open unit interval so `logit` stays finite.
fn clamp_probability(p: f64) -> f64 {
    p.clamp(1e-12, 1.0 - 1e-12)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() <= eps
    }

    /// EM recovers the planted paralog fraction on a well-separated two-point
    /// mixture: reals at LR −15 (posterior ≈ 0), paralogs at LR +15
    /// (posterior ≈ 1), so `π → fraction of paralogs`.
    #[test]
    fn em_recovers_pi_on_two_component_mixture() {
        let pi_true = 0.1;
        let n = 10_000u64;
        let mut h = ParalogLrHistogram::with_defaults();
        for _ in 0..((1.0 - pi_true) * n as f64) as u64 {
            h.push(-15.0);
        }
        for _ in 0..(pi_true * n as f64) as u64 {
            h.push(15.0);
        }
        let prior = ParalogPrior::estimate(&h, &EmConfig::default());
        assert!(
            approx(prior.prior_probability, pi_true, 1e-3),
            "π = {}",
            prior.prior_probability
        );
    }

    /// The 50%-posterior cut is `−logit π > 0` (never `LR = 0`), and the
    /// posterior at that LR is exactly `0.5`.
    #[test]
    fn half_posterior_cut_is_positive_and_consistent() {
        let prior = ParalogPrior {
            prior_probability: 0.09,
            converged: true,
        };
        let half = prior.half_posterior_ratio();
        assert!(half > 0.0, "half cut = {half}");
        assert!(approx(prior.paralog_posterior(half), 0.5, 1e-12));
        // Posterior is monotone increasing in LR.
        assert!(prior.paralog_posterior(0.0) < prior.paralog_posterior(5.0));
    }

    /// `q_of_lr` is monotone non-increasing in LR (higher LR ⇒ smaller tail
    /// FDR) across a sweep.
    #[test]
    fn q_of_lr_is_monotone_decreasing() {
        let mut h = ParalogLrHistogram::with_defaults();
        // A spread of LRs so many bins are populated.
        for i in 0..2000 {
            h.push(-30.0 + 0.03 * i as f64);
        }
        let prior = ParalogPrior::estimate(&h, &EmConfig::default());
        let curve = ParalogFdrCurve::from_histogram(&h, &prior);
        let mut prev = f64::INFINITY;
        let mut lr = -40.0;
        while lr <= 40.0 {
            let q = curve.q_of_lr(lr);
            assert!(q <= prev + 1e-12, "q not monotone at LR={lr}: {q} > {prev}");
            prev = q;
            lr += 0.5;
        }
    }

    /// A higher target FDR admits a lower (less-stringent) LR threshold.
    #[test]
    fn lr_threshold_relaxes_with_higher_fdr() {
        let mut h = ParalogLrHistogram::with_defaults();
        for i in 0..5000 {
            h.push(-25.0 + 0.01 * i as f64);
        }
        let prior = ParalogPrior::estimate(&h, &EmConfig::default());
        let curve = ParalogFdrCurve::from_histogram(&h, &prior);
        let t01 = curve.lr_threshold_for_fdr(0.01).expect("1% achievable");
        let t20 = curve.lr_threshold_for_fdr(0.20).expect("20% achievable");
        assert!(t20 <= t01, "20% cut {t20} should be <= 1% cut {t01}");
    }

    /// The histogram-derived `π` and per-locus q-values match a brute-force
    /// full-vector computation (the prototype's algorithm) within bin
    /// resolution — the bounded-RAM form is faithful.
    #[test]
    fn histogram_matches_full_vector_reference() {
        // A deterministic mixture of LRs spanning the sigmoid transition.
        let mut lrs: Vec<f64> = Vec::new();
        for i in 0..900 {
            lrs.push(-12.0 + 0.02 * (i % 200) as f64); // reals, ~[-12, -8]
        }
        for i in 0..100 {
            lrs.push(6.0 + 0.05 * (i % 100) as f64); // paralogs, ~[6, 11]
        }

        // Full-vector EM (prototype).
        let mut pi = DEFAULT_EM_START;
        for _ in 0..DEFAULT_EM_MAX_ITER {
            let off = logit(clamp_probability(pi));
            let pi_new = lrs.iter().map(|&x| sigmoid(x + off)).sum::<f64>() / lrs.len() as f64;
            if (pi_new - pi).abs() < DEFAULT_EM_TOL {
                pi = pi_new;
                break;
            }
            pi = pi_new;
        }

        // Histogram EM.
        let mut h = ParalogLrHistogram::with_defaults();
        for &x in &lrs {
            h.push(x);
        }
        let prior = ParalogPrior::estimate(&h, &EmConfig::default());
        assert!(
            approx(prior.prior_probability, pi, 2e-3),
            "hist π {} vs full-vector π {pi}",
            prior.prior_probability
        );

        // Full-vector tail FDR at a few LR thresholds vs the curve.
        let curve = ParalogFdrCurve::from_histogram(&h, &prior);
        let off = logit(clamp_probability(pi));
        for &threshold in &[-5.0, 0.0, 4.0, 8.0] {
            let tail: Vec<f64> = lrs
                .iter()
                .filter(|&&x| x >= threshold)
                .map(|&x| 1.0 - sigmoid(x + off))
                .collect();
            if tail.is_empty() {
                continue;
            }
            let full_q = tail.iter().sum::<f64>() / tail.len() as f64;
            assert!(
                approx(curve.q_of_lr(threshold), full_q, 5e-3),
                "q_of_lr({threshold}) {} vs full {full_q}",
                curve.q_of_lr(threshold)
            );
        }
    }

    /// Non-finite LRs are dropped, and an empty histogram estimates `π =
    /// start` rather than dividing by zero.
    #[test]
    fn empty_and_nonfinite_are_handled() {
        let mut h = ParalogLrHistogram::with_defaults();
        h.push(f64::NAN);
        h.push(f64::INFINITY);
        assert!(h.is_empty());
        let prior = ParalogPrior::estimate(&h, &EmConfig::default());
        assert_eq!(prior.prior_probability, DEFAULT_EM_START);
    }

    /// Out-of-range LRs saturate into the end bins rather than panicking.
    #[test]
    fn out_of_range_lrs_saturate() {
        let mut h = ParalogLrHistogram::new(-10.0, 10.0, 100).expect("valid range");
        h.push(-1e9);
        h.push(1e9);
        assert_eq!(h.total(), 2);
        assert_eq!(h.counts[0], 1);
        assert_eq!(h.counts[99], 1);
    }

    /// `ParalogLrHistogram::new` rejects a degenerate range (matching the
    /// `GridSpec::new` convention) rather than panicking.
    #[test]
    fn histogram_new_rejects_degenerate_range() {
        assert!(ParalogLrHistogram::new(-10.0, 10.0, 100).is_some());
        assert!(ParalogLrHistogram::new(10.0, -10.0, 100).is_none()); // lo > hi
        assert!(ParalogLrHistogram::new(-10.0, 10.0, 0).is_none()); // no bins
        assert!(ParalogLrHistogram::new(f64::NAN, 10.0, 100).is_none()); // non-finite
    }

    /// EM that cannot reach tolerance within `max_iter` returns the last
    /// iterate flagged `converged = false`; a normal run flags `true`.
    #[test]
    fn em_flags_non_convergence() {
        let mut h = ParalogLrHistogram::with_defaults();
        for _ in 0..1000 {
            h.push(-15.0);
        }
        for _ in 0..100 {
            h.push(15.0);
        }
        let one_iter = EmConfig {
            max_iter: 1,
            tol: 1e-12,
            ..Default::default()
        };
        assert!(!ParalogPrior::estimate(&h, &one_iter).converged);
        assert!(ParalogPrior::estimate(&h, &EmConfig::default()).converged);
    }

    /// `lr_threshold_for_fdr` returns `None` for an unachievable target and a
    /// finite LR when the target is trivially met.
    #[test]
    fn lr_threshold_boundaries() {
        let mut h = ParalogLrHistogram::with_defaults();
        for i in 0..2000 {
            h.push(-20.0 + 0.02 * i as f64);
        }
        let prior = ParalogPrior::estimate(&h, &EmConfig::default());
        let curve = ParalogFdrCurve::from_histogram(&h, &prior);
        // A negative FDR target is unachievable (q >= 0 everywhere).
        assert!(curve.lr_threshold_for_fdr(-0.1).is_none());
        // Any q ≤ 1 always holds, so a target of 1.0 admits the lowest bin.
        assert!(curve.lr_threshold_for_fdr(1.0).is_some());
    }

    /// A degenerate π (`0` / `1`, reachable via the `pub` field) is clamped
    /// so the posterior and half-cut stay finite.
    #[test]
    fn degenerate_pi_stays_finite() {
        for pi in [0.0, 1.0] {
            let prior = ParalogPrior {
                prior_probability: pi,
                converged: true,
            };
            assert!(prior.paralog_posterior(3.0).is_finite());
            assert!(prior.half_posterior_ratio().is_finite());
        }
    }

    /// The FDR curve of an empty histogram is well-defined everywhere via the
    /// empty-tail `local_lfdr` fallback (no division by a zero tail count).
    #[test]
    fn empty_histogram_curve_is_finite() {
        let h = ParalogLrHistogram::with_defaults();
        let prior = ParalogPrior::estimate(&h, &EmConfig::default());
        let curve = ParalogFdrCurve::from_histogram(&h, &prior);
        for &lr in &[-50.0, 0.0, 50.0] {
            assert!(curve.q_of_lr(lr).is_finite());
        }
    }

    /// `sigmoid` is stable and correct at the extremes and the origin.
    #[test]
    fn sigmoid_is_stable() {
        assert!(approx(sigmoid(0.0), 0.5, 1e-12));
        assert!(approx(sigmoid(1000.0), 1.0, 1e-12));
        assert!(approx(sigmoid(-1000.0), 0.0, 1e-12));
        assert!(sigmoid(-1000.0).is_finite());
    }
}
