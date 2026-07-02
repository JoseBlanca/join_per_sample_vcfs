//! S1 — the up-front per-sample state and the running `Hexp` accumulator.
//!
//! Two containers the two-pass flow sets up before (and during) the main caller
//! pass, holding the quantities the paralog score needs that are *not* per-locus
//! spill payload:
//!
//! - [`ParalogPrePass`] — fit once at reader-open, from each sample's `.psp`
//!   metadata ([`SampleSummary`]): the single-copy [`SingleCopyCoverageModel`]
//!   (Q2) and the observed-heterozygosity rate `obs_het` (Q4), plus the
//!   callable-position count. The inbreeding coefficient `F` is **not** formed
//!   here — it needs the cohort `Hexp`, which is not known until the whole
//!   genome has streamed (arch §3 step 1, §6).
//! - [`HexpAccumulator`] — folds each locus's cohort allele frequencies into the
//!   running expected-heterozygosity sum during the main caller pass, and
//!   finalises to `Hexp = Σ expected-het / callable_ref` on the
//!   **per-callable-position scale** the observed-het rate lives on (the R1
//!   correction; arch §6). `F` is formed from this finished `Hexp` in S4.
//!
//! Both are memory-flat: the pre-pass is `O(samples)` small structs held for the
//! run, and the accumulator is two running scalars.

use crate::paralog::{CoverageFitConfig, SingleCopyCoverageModel, obs_het};
use crate::sample_summary::SampleSummary;

/// One sample's up-front fitted state, formed from its `.psp` metadata before
/// any records stream.
///
/// Carries everything the paralog score needs per sample *except* the
/// inbreeding coefficient `F`, which is formed later from the finished cohort
/// `Hexp` (S4). A sample whose coverage fit is rejected, or whose summary is
/// absent, is not represented by this struct at all — [`ParalogPrePass`] holds
/// it as `None` (absent, not fatal).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ParalogSampleModel {
    /// The fitted single-copy coverage model (Q2): a window's
    /// `(gc, depth)` → relative copy number, plus σ₀.
    pub coverage_model: SingleCopyCoverageModel,
    /// The observed-heterozygosity **rate** `n_het_sites / callable_positions`
    /// (Q4). One half of `F = 1 − obs_het / Hexp`.
    pub obs_het: f64,
    /// The sample's callable-position count (the het-rate denominator, and a
    /// candidate for the cohort `Hexp` denominator).
    pub callable_positions: u64,
}

/// The up-front, per-sample paralog state: one optional [`ParalogSampleModel`]
/// per cohort sample, indexed parallel to the cohort's sample columns.
///
/// A `None` entry marks a sample carried as **absent** — either its `.psp`
/// summary was missing or its coverage fit was rejected (a degenerate,
/// essentially-uncovered sample). Absence is never fatal: the sample simply
/// contributes no coverage evidence to any locus score.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ParalogPrePass {
    /// Per cohort sample, in column order; `None` = absent / rejected.
    samples: Vec<Option<ParalogSampleModel>>,
    /// The analysis/GC window width the coverage histograms were tiled at (the
    /// `pileup --gc-window-bp`, `500` by default). The window-depth accumulator
    /// (S6c) must tile at the **same** width so a locus's window mean depth is
    /// on the scale the coverage model was fit on. Taken from the first present
    /// summary (all cohort samples share one pileup window).
    window_bp: u32,
}

impl ParalogPrePass {
    /// Fit each sample's coverage model + `obs_het` from its `.psp` metadata.
    ///
    /// `summaries` is indexed parallel to the cohort's sample columns; a `None`
    /// entry (no summary for that column) becomes an absent sample. A summary
    /// whose coverage fit is rejected ([`SingleCopyCoverageModel::fit`] errors —
    /// e.g. a degenerate/too-shallow sample) is likewise carried as absent, not
    /// propagated as an error: one bad sample must not abort the whole cohort.
    pub(crate) fn fit(summaries: &[Option<SampleSummary>], cfg: &CoverageFitConfig) -> Self {
        let samples = summaries
            .iter()
            .map(|summary| {
                let summary = summary.as_ref()?;
                let model = SingleCopyCoverageModel::fit(&summary.coverage_by_gc, cfg).ok()?;
                let callable_positions = summary.coverage_by_gc.callable_positions;
                Some(ParalogSampleModel {
                    coverage_model: model,
                    obs_het: obs_het(&summary.heterozygosity, callable_positions),
                    callable_positions,
                })
            })
            .collect();
        // The window is a pileup property shared by the cohort, so read it from
        // the first present summary (even a fit-rejected one carries it);
        // default to the pileup default if no summary is present at all.
        let window_bp = summaries
            .iter()
            .flatten()
            .map(|s| s.coverage_by_gc.window_bp)
            .next()
            .unwrap_or(crate::sample_summary::DEFAULT_GC_WINDOW_BP)
            .max(1);
        Self { samples, window_bp }
    }

    /// The analysis/GC window width (bp) the coverage histograms were tiled at
    /// — the width the window-depth accumulator (S6c) must reuse.
    pub(crate) fn window_bp(&self) -> u32 {
        self.window_bp
    }

    /// The per-sample fitted state, in cohort column order (`None` = absent).
    pub(crate) fn samples(&self) -> &[Option<ParalogSampleModel>] {
        &self.samples
    }

    /// The number of samples that fit (non-absent entries).
    pub(crate) fn fit_count(&self) -> usize {
        self.samples.iter().filter(|s| s.is_some()).count()
    }

    /// The representative callable-position count for the cohort `Hexp`
    /// denominator: the **median** callable count over the samples that fit
    /// (arch §6 / R1). `None` if no sample fit — the cohort has no usable
    /// coverage, so `Hexp` cannot be formed.
    ///
    /// Uses the lower-middle element on an even count (`len / 2` after sort), a
    /// deterministic convention matching the R1 harness.
    ///
    /// Deliberately drops the R1 harness's `.max(1)` fallback (which coerced an
    /// empty cohort to a callable count of `1`): here a no-fit cohort returns
    /// `None`, [`HexpAccumulator::finish`] maps a zero/absent reference to
    /// `Hexp = 0`, and [`crate::paralog::inbreeding_coefficient`] maps a zero
    /// `Hexp` to `F = 0` (outbred). Both routes reach the same defined-outbred
    /// end state, so the `.max(1)` guard is redundant.
    pub(crate) fn callable_reference(&self) -> Option<u64> {
        let mut callables: Vec<u64> = self
            .samples
            .iter()
            .filter_map(|s| s.as_ref().map(|s| s.callable_positions))
            .collect();
        if callables.is_empty() {
            return None;
        }
        callables.sort_unstable();
        Some(callables[callables.len() / 2])
    }

    /// Per-sample σ₀ (one-copy relative-depth SD), indexed parallel to the
    /// cohort columns, with `1.0` for absent samples. This is the slice
    /// [`crate::paralog::score_locus_for_paralogy`] consumes; the `1.0`
    /// placeholder for an absent sample is inert because an absent sample never
    /// contributes a [`crate::paralog::SampleObservation`] at any locus.
    pub(crate) fn single_copy_depth_sd(&self) -> Vec<f64> {
        self.samples
            .iter()
            .map(|s| {
                s.as_ref()
                    .map(|s| s.coverage_model.single_copy_depth_sd())
                    .unwrap_or(1.0)
            })
            .collect()
    }
}

/// The running accumulator for the cohort's expected heterozygosity `Hexp`,
/// folded during the main caller pass and finalised once the genome has
/// streamed.
///
/// It keeps only the running sum of per-locus expected heterozygosity (`Σ`) and
/// the count of variant loci folded — two scalars, memory-flat in variant
/// count. The site expected heterozygosity is `1 − Σ pᵢ²` over the cohort
/// allele frequencies, which for a biallelic locus is exactly `2p(1−p)` (the
/// `Σ2pq` of the architecture doc); monomorphic positions contribute `0`, so
/// only variant loci need be folded.
///
/// [`finish`] divides by a representative callable count, putting `Hexp` on the
/// **per-callable-position scale** the observed-het rate lives on — *not* the
/// mean over variant sites (the bug R1 caught, which saturates every `F` at
/// `0.99` and inflates `π`).
///
/// [`finish`]: HexpAccumulator::finish
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub(crate) struct HexpAccumulator {
    /// Running `Σ (1 − Σ pᵢ²)` over the variant loci folded so far.
    expected_het_sum: f64,
    /// Number of loci folded with a finite contribution (diagnostics; not the
    /// `Hexp` denominator). In the pipeline only variant loci are fed.
    loci: u64,
}

impl HexpAccumulator {
    /// A fresh accumulator with an empty sum.
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Fold one locus's cohort allele frequencies into the running sum.
    ///
    /// `allele_frequencies` is the per-allele estimated frequency `p̂` (the
    /// [`crate::var_calling::posterior_engine::PosteriorRecord::allele_frequencies`]
    /// the caller pass produces). The site expected heterozygosity `1 − Σ pᵢ²`
    /// is floored at `0` (allele frequencies need not sum to exactly `1`, so the
    /// raw value can dip a hair below zero) and non-finite contributions are
    /// dropped rather than poisoning the sum. An empty slice yields `Σ pᵢ² = 0`
    /// → a contribution of `1.0`; the caller pass never emits an alleleless
    /// record, so this corner is inert.
    pub(crate) fn observe(&mut self, allele_frequencies: &[f64]) {
        let sum_sq: f64 = allele_frequencies.iter().map(|p| p * p).sum();
        let site_het = 1.0 - sum_sq;
        // Screen finiteness *before* the `0` floor: `f64::max` treats `NaN` as
        // the missing operand, so `NaN.max(0.0)` is `0.0` — flooring first would
        // silently fold a poisoned locus in as zero rather than dropping it.
        if site_het.is_finite() {
            self.expected_het_sum += site_het.max(0.0);
            self.loci += 1;
        }
    }

    /// The running expected-heterozygosity sum `Σ (1 − Σ pᵢ²)` folded so far.
    pub(crate) fn expected_het_sum(&self) -> f64 {
        self.expected_het_sum
    }

    /// Variant loci folded so far.
    pub(crate) fn loci(&self) -> u64 {
        self.loci
    }

    /// Finalise to `Hexp = Σ expected-het / callable_ref`, on the
    /// per-callable-position scale. Returns `0.0` for a zero (or absent)
    /// callable reference — a cohort with no callable positions has no defined
    /// expected heterozygosity, and a `0` `Hexp` yields `F = 0` (treat every
    /// sample as outbred) rather than a division by zero.
    pub(crate) fn finish(&self, callable_ref: u64) -> f64 {
        if callable_ref == 0 {
            return 0.0;
        }
        self.expected_het_sum / callable_ref as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sample_summary::{CoverageByGcHistogram, HetCounts};

    /// A single-copy coverage-by-GC histogram that [`SingleCopyCoverageModel`]
    /// accepts: a symmetric single-copy depth spread in both GC bins, dense
    /// enough to seed the GC curve. Mirrors the coverage-model unit tests.
    fn single_copy_histogram(callable_positions: u64) -> CoverageByGcHistogram {
        let gc_bins = 2u32;
        let depth_bins = 40u32;
        let width = 1.0f64;
        let row_stride = depth_bins as usize + 1;
        let mut counts = vec![0u32; gc_bins as usize * row_stride];
        let mut n_tiles = 0u64;
        // (gc_fraction, depth, count) tiles: peak at depth-bin 10 in both GC bins.
        for &(gc, depth, count) in &[
            (0.25, 9.5, 100u32),
            (0.25, 10.5, 300),
            (0.25, 11.5, 100),
            (0.75, 9.5, 100),
            (0.75, 10.5, 300),
            (0.75, 11.5, 100),
        ] {
            let gb = ((gc * gc_bins as f64) as usize).min(gc_bins as usize - 1);
            let db = ((depth / width) as usize).min(depth_bins as usize);
            counts[gb * row_stride + db] += count;
            n_tiles += u64::from(count);
        }
        CoverageByGcHistogram {
            window_bp: 500,
            gc_bins,
            depth_bin_width: width,
            depth_bins,
            n_tiles,
            n_skipped_tiles: 0,
            callable_positions,
            counts,
        }
    }

    fn het_counts(n_het: u64) -> HetCounts {
        HetCounts {
            n_het_sites: n_het,
            n_hom_alt_sites: 0,
            n_ambiguous_sites: 0,
            n_variant_sites: n_het,
            min_depth: 4,
            error_rate: 0.02,
            lr_margin: std::f64::consts::LN_10,
        }
    }

    fn summary(n_het: u64, callable_positions: u64) -> SampleSummary {
        SampleSummary {
            version: crate::sample_summary::SAMPLE_SUMMARY_VERSION,
            coverage_by_gc: single_copy_histogram(callable_positions),
            heterozygosity: het_counts(n_het),
        }
    }

    /// A histogram too shallow to anchor (bottom-bin-only) — the coverage fit
    /// rejects it, so the pre-pass must carry the sample as absent.
    fn degenerate_summary() -> SampleSummary {
        let mut s = summary(0, 500);
        // All mass in depth bin 0 → `DepthModeAtBottomBin` rejection.
        let row_stride = s.coverage_by_gc.depth_bins as usize + 1;
        s.coverage_by_gc.counts = vec![0u32; s.coverage_by_gc.gc_bins as usize * row_stride];
        s.coverage_by_gc.counts[0] = 500;
        s.coverage_by_gc.counts[row_stride] = 500;
        s.coverage_by_gc.n_tiles = 1000;
        s.coverage_by_gc.callable_positions = 1000;
        s
    }

    /// A fitting sample yields a populated model, its het rate, and its callable
    /// count; the σ₀ slice reports the fitted SD.
    #[test]
    fn fit_populates_a_healthy_sample() {
        let summaries = [Some(summary(359, 1_000_000))];
        let prepass = ParalogPrePass::fit(&summaries, &CoverageFitConfig::default());
        assert_eq!(prepass.fit_count(), 1);
        let model = prepass.samples()[0].as_ref().expect("fit");
        assert!((model.obs_het - 359e-6).abs() < 1e-12);
        assert_eq!(model.callable_positions, 1_000_000);
        assert!(model.coverage_model.single_copy_depth_sd() > 0.0);
        assert_eq!(prepass.single_copy_depth_sd().len(), 1);
        assert!(prepass.single_copy_depth_sd()[0] > 0.0);
    }

    /// A degenerate sample (coverage fit rejected) and a `None` summary are both
    /// carried as absent — not fatal — and the σ₀ slice reports `1.0` for them.
    #[test]
    fn degenerate_and_missing_samples_are_carried_absent() {
        let summaries = [
            Some(summary(100, 500_000)),
            Some(degenerate_summary()),
            None,
        ];
        let prepass = ParalogPrePass::fit(&summaries, &CoverageFitConfig::default());
        assert_eq!(prepass.fit_count(), 1);
        assert!(prepass.samples()[0].is_some());
        assert!(prepass.samples()[1].is_none(), "fit-rejected → absent");
        assert!(prepass.samples()[2].is_none(), "missing summary → absent");
        let sigma = prepass.single_copy_depth_sd();
        assert_eq!(sigma.len(), 3);
        assert!(sigma[0] > 0.0);
        assert_eq!(sigma[1], 1.0, "absent sample → inert σ₀ placeholder");
        assert_eq!(sigma[2], 1.0);
    }

    /// The callable reference is the median callable count over the fit samples
    /// (lower-middle on an even count), ignoring absent ones.
    #[test]
    fn callable_reference_is_median_over_fit_samples() {
        let summaries = [
            Some(summary(1, 100)),
            Some(summary(1, 300)),
            None, // absent — excluded
            Some(summary(1, 200)),
        ];
        let prepass = ParalogPrePass::fit(&summaries, &CoverageFitConfig::default());
        // Sorted fit callables: [100, 200, 300] → median 200.
        assert_eq!(prepass.callable_reference(), Some(200));
    }

    /// With no sample fitting, there is no callable reference.
    #[test]
    fn callable_reference_none_when_no_sample_fits() {
        let summaries = [Some(degenerate_summary()), None];
        let prepass = ParalogPrePass::fit(&summaries, &CoverageFitConfig::default());
        assert_eq!(prepass.fit_count(), 0);
        assert_eq!(prepass.callable_reference(), None);
    }

    /// `HexpAccumulator` matches an independent `Σ 2pq` over biallelic loci, and
    /// `finish` divides by the callable reference (per-callable-position scale).
    #[test]
    fn hexp_matches_independent_twopq_sum() {
        // Biallelic loci given as [p_ref, p_alt]; expected het = 2·p_alt·(1−p_alt).
        let alt_freqs = [0.1, 0.25, 0.5, 0.02, 0.4];
        let mut acc = HexpAccumulator::new();
        let mut independent_sum = 0.0;
        for &p in &alt_freqs {
            acc.observe(&[1.0 - p, p]);
            independent_sum += 2.0 * p * (1.0 - p);
        }
        assert_eq!(acc.loci(), alt_freqs.len() as u64);
        assert!((acc.expected_het_sum() - independent_sum).abs() < 1e-12);

        let callable_ref = 1_000_000u64;
        let hexp = acc.finish(callable_ref);
        assert!((hexp - independent_sum / callable_ref as f64).abs() < 1e-18);
    }

    /// The general form `1 − Σ pᵢ²` handles a multi-allelic locus (three
    /// alleles) — expected het is `1 − (p₀² + p₁² + p₂²)`.
    #[test]
    fn hexp_handles_multiallelic_site() {
        let freqs = [0.5, 0.3, 0.2];
        let mut acc = HexpAccumulator::new();
        acc.observe(&freqs);
        let expected = 1.0 - (0.25 + 0.09 + 0.04);
        assert!((acc.expected_het_sum() - expected).abs() < 1e-12);
    }

    /// A monomorphic locus (`p = 1`) contributes zero expected heterozygosity;
    /// a non-finite frequency contribution is dropped.
    #[test]
    fn hexp_monomorphic_and_nonfinite_are_handled() {
        let mut acc = HexpAccumulator::new();
        acc.observe(&[1.0]); // monomorphic → 0
        assert_eq!(acc.expected_het_sum(), 0.0);
        assert_eq!(acc.loci(), 1);
        acc.observe(&[f64::NAN, 0.5]); // dropped
        assert_eq!(acc.loci(), 1, "non-finite contribution not counted");
        assert!(acc.expected_het_sum().is_finite());
    }

    /// A zero (or absent) callable reference yields `Hexp = 0` rather than a
    /// division by zero.
    #[test]
    fn hexp_finish_zero_callable_is_zero() {
        let mut acc = HexpAccumulator::new();
        acc.observe(&[0.5, 0.5]);
        assert_eq!(acc.finish(0), 0.0);
    }
}
