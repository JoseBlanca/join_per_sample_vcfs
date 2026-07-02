//! S1 — the up-front per-sample state.
//!
//! [`ParalogPrePass`] is fit once at reader-open, from each sample's `.psp`
//! metadata ([`SampleSummary`]): the single-copy [`SingleCopyCoverageModel`]
//! (Q2), which is the *only* per-sample quantity the paralog score needs from a
//! sample's own `.psp`. The inbreeding coefficient `F` is **not** per-sample: it
//! is one cohort-level number taken from the caller's `--inbreeding-coefficient`
//! (see `doc/devel/architecture/hidden_paralog_single_sample_scoring.md`), so the
//! pre-pass no longer computes an observed heterozygosity or a cohort expected
//! heterozygosity (`Hexp`). The old `HexpAccumulator` — which folded every
//! locus's allele frequencies during the main pass — is gone; nothing global is
//! accumulated for the score any more (only π, over loci, in the calibrate pass).
//!
//! Memory-flat: `O(samples)` small structs held for the run.

use crate::paralog::{CoverageFitConfig, SingleCopyCoverageModel};
use crate::sample_summary::SampleSummary;

/// One sample's up-front fitted state, formed from its `.psp` metadata before
/// any records stream.
///
/// Just the single-copy coverage model — the only per-sample quantity the score
/// needs from a sample's own `.psp`. The inbreeding coefficient `F` is *not*
/// per-sample: it is one cohort-level number taken from the caller's
/// `--inbreeding-coefficient` (see
/// `doc/devel/architecture/hidden_paralog_single_sample_scoring.md`), so no
/// per-sample heterozygosity is stored here. A sample whose coverage fit is
/// rejected, or whose summary is absent, is not represented by this struct at
/// all — [`ParalogPrePass`] holds it as `None` (absent, not fatal).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ParalogSampleModel {
    /// The fitted single-copy coverage model (Q2): a window's
    /// `(gc, depth)` → relative copy number, plus σ₀.
    pub coverage_model: SingleCopyCoverageModel,
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
    /// Fit each sample's coverage model from its `.psp` metadata.
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
                Some(ParalogSampleModel {
                    coverage_model: model,
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

    /// A minimal het-count block. The score no longer reads it (`F` is
    /// cohort-level), but [`SampleSummary`] still carries the field.
    fn het_counts() -> HetCounts {
        HetCounts {
            n_het_sites: 100,
            n_hom_alt_sites: 0,
            n_ambiguous_sites: 0,
            n_variant_sites: 100,
            min_depth: 4,
            error_rate: 0.02,
            lr_margin: std::f64::consts::LN_10,
        }
    }

    fn summary(callable_positions: u64) -> SampleSummary {
        SampleSummary {
            version: crate::sample_summary::SAMPLE_SUMMARY_VERSION,
            coverage_by_gc: single_copy_histogram(callable_positions),
            heterozygosity: het_counts(),
        }
    }

    /// A histogram too shallow to anchor (bottom-bin-only) — the coverage fit
    /// rejects it, so the pre-pass must carry the sample as absent.
    fn degenerate_summary() -> SampleSummary {
        let mut s = summary(500);
        // All mass in depth bin 0 → `DepthModeAtBottomBin` rejection.
        let row_stride = s.coverage_by_gc.depth_bins as usize + 1;
        s.coverage_by_gc.counts = vec![0u32; s.coverage_by_gc.gc_bins as usize * row_stride];
        s.coverage_by_gc.counts[0] = 500;
        s.coverage_by_gc.counts[row_stride] = 500;
        s.coverage_by_gc.n_tiles = 1000;
        s.coverage_by_gc.callable_positions = 1000;
        s
    }

    /// A fitting sample yields a populated coverage model; the σ₀ slice reports
    /// the fitted SD.
    #[test]
    fn fit_populates_a_healthy_sample() {
        let summaries = [Some(summary(1_000_000))];
        let prepass = ParalogPrePass::fit(&summaries, &CoverageFitConfig::default());
        assert_eq!(prepass.fit_count(), 1);
        let model = prepass.samples()[0].as_ref().expect("fit");
        assert!(model.coverage_model.single_copy_depth_sd() > 0.0);
        assert_eq!(prepass.single_copy_depth_sd().len(), 1);
        assert!(prepass.single_copy_depth_sd()[0] > 0.0);
    }

    /// A degenerate sample (coverage fit rejected) and a `None` summary are both
    /// carried as absent — not fatal — and the σ₀ slice reports `1.0` for them.
    #[test]
    fn degenerate_and_missing_samples_are_carried_absent() {
        let summaries = [Some(summary(500_000)), Some(degenerate_summary()), None];
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

    /// The window width is read from the first present summary (all cohort
    /// samples share one pileup window), defaulting when none is present.
    #[test]
    fn window_bp_from_first_present_summary() {
        let with = ParalogPrePass::fit(&[None, Some(summary(1000))], &CoverageFitConfig::default());
        assert_eq!(with.window_bp(), 500);
        let without = ParalogPrePass::fit(&[None], &CoverageFitConfig::default());
        assert_eq!(
            without.window_bp(),
            crate::sample_summary::DEFAULT_GC_WINDOW_BP
        );
    }
}
