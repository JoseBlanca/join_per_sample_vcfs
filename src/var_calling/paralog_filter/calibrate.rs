//! S4 — the per-locus paralog scorer + calibration (LR histogram → π + FDR cut).
//!
//! Two responsibilities after the inline-scoring change (arch:
//! `hidden_paralog_inline_scoring.md`):
//!
//! 1. **The per-locus scorer** — [`score_spilled_locus`], a pure function of a
//!    locus's record + window coverage + the pre-pass. In production the caller
//!    **worker** runs it once per locus through [`ParalogScoringContext`] (built
//!    once from the pre-pass + cohort `F` + params, shared read-only across
//!    workers) and stores the LR on the spill; nothing re-scores downstream.
//! 2. **The calibration** — [`calibrate_from_histogram`] estimates the
//!    empirical-Bayes prior `π` (Q5, with a fixed-prior fallback if the EM does
//!    not converge), builds the tail-FDR `q_of_lr` curve, and resolves the LR
//!    threshold for the operator's target FDR. It runs on the histogram the sink
//!    folds inline during the main pass — no spill read, no scoring.
//!
//! The spill-streaming [`calibrate`] (score every record into a histogram, then
//! calibrate) and its [`score_spill_record`] wrapper survive only as `#[cfg(test)]`
//! helpers that drive the scorer end-to-end from record fixtures.

#[cfg(test)]
use std::io::Read;

use crate::paralog::{
    EmConfig, LocusObservations, ParalogFdrCurve, ParalogLrHistogram, ParalogModelParams,
    ParalogPrior, ParalogScorePrecompute, SampleObservation, score_locus_for_paralogy,
};

use crate::var_calling::posterior_engine::PosteriorRecord;
use crate::var_calling::types::LocusWindowCoverage;

use super::prepass::ParalogPrePass;
#[cfg(test)]
use super::spill::{ParalogSpillReader, ParalogSpillRecord, SpillError};

/// Default fallback `π` used when the EM does not converge — a rare-paralog
/// prior deliberately low so a pathological dataset yields a conservative
/// (few-drops) calibration rather than trusting a junk estimate. Chosen to
/// match [`crate::paralog::prior::DEFAULT_EM_START`] (`0.03`); they are
/// independent constants, so this comment is the only thing coupling them.
/// The operator is warned when this is used (S6).
pub(crate) const DEFAULT_FALLBACK_PARALOG_PRIOR: f64 = 0.03;

/// Tuning knobs for [`calibrate`]. [`Default`] is the production configuration.
///
/// There is deliberately **no minimum-sample gate**: an under-powered locus
/// self-gates in the likelihood ratio (when the coverage can't tell 1× from 2×,
/// H1 and H2 fit the data equally → LR≈0 → kept), so a count threshold would be
/// redundant. See `doc/devel/architecture/hidden_paralog_single_sample_scoring.md`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CalibrationConfig {
    /// EM configuration for the prior estimate.
    pub em: EmConfig,
    /// Fallback `π` when the EM does not converge.
    pub fallback_prior: f64,
}

impl Default for CalibrationConfig {
    fn default() -> Self {
        Self {
            em: EmConfig::default(),
            fallback_prior: DEFAULT_FALLBACK_PARALOG_PRIOR,
        }
    }
}

/// The global calibration derived from the spill: the prior, the FDR curve, and
/// the resolved LR cut for the target FDR.
#[derive(Debug, Clone)]
pub(crate) struct ParalogCalibration {
    /// The empirical-Bayes prior `π` (from the EM, or the fallback). Its
    /// `converged` flag is `false` when the fallback was used.
    pub prior: ParalogPrior,
    /// The tail-FDR curve `q_of_lr`, built over the same histogram bins the LRs
    /// were folded into.
    pub curve: ParalogFdrCurve,
    /// The least-stringent LR at/above which a locus's tail FDR is `<=
    /// target_fdr`; `None` if the target is unachievable (nothing is flagged).
    /// Recorded for VCF-header provenance (S5).
    pub lr_threshold: Option<f64>,
    /// The operator's target false-discovery rate — the operating knob.
    pub target_fdr: f64,
}

impl ParalogCalibration {
    /// Whether a locus with likelihood ratio `lr` is flagged (dropped): its
    /// tail-FDR q-value is at or below the target. Uses the curve directly, so
    /// an unachievable target flags nothing (q never reaches it). Equivalent to
    /// `lr >= lr_threshold` by the curve's monotonicity, but stated on the FDR
    /// so the guarantee is explicit.
    ///
    /// A non-finite `lr` is never flagged: the histogram itself refuses to fold
    /// a non-finite LR (it is not a valid Bayes factor), so it contributed to
    /// neither `π` nor the curve — classifying on it here would be treating the
    /// same value by a different rule than the calibration did. Screened
    /// explicitly because `q_of_lr` would otherwise saturate a `NaN` into bin 0.
    pub(crate) fn flags(&self, lr: f64) -> bool {
        lr.is_finite() && self.curve.q_of_lr(lr) <= self.target_fdr
    }
}

/// The cohort inbreeding coefficient `F`, one value for every sample.
///
/// `F` is *not* per-individual (it isn't identifiable from one sample, and the
/// per-individual het-rate proxy is divergence-contaminated — see
/// `doc/devel/architecture/hidden_paralog_single_sample_scoring.md`): it is the
/// caller's single cohort `--inbreeding-coefficient`, the same one the SNP
/// caller uses in its genotype prior. The write pass (S5) rebuilds this slice
/// identically so its `score_spilled_locus` inputs are bit-identical to the
/// calibrate pass's.
pub(crate) fn cohort_inbreeding(n_samples: usize, f: f64) -> Vec<f64> {
    vec![f; n_samples]
}

/// Whether the record is a biallelic SNP — the class the paralog score is
/// defined over (a single-base REF and a single-base ALT). Other loci (indels,
/// multiallelics) are not scored and are never flagged.
fn is_biallelic_snp(record: &PosteriorRecord) -> bool {
    record.alleles.len() == 2
        && record.alleles[0].seq.len() == 1
        && record.alleles[1].seq.len() == 1
}

/// One sample's paralog observation at a locus, or `None` if the sample is
/// unusable here (absent from the pre-pass, no covered window, or no reads).
///
/// `window` is the locus's per-sample centred-window coverage, read straight
/// from the record spill ([`ParalogSpillRecord::window_coverage`], gathered by
/// the caller from the psp windowed columns). Each sample carries its own GC +
/// mean coverage; a `NaN` in either means the sample had no covered record at
/// the locus and is skipped. Together with the record's per-sample AD they form
/// the score's per-sample observation.
fn build_observation(
    record: &PosteriorRecord,
    window: &LocusWindowCoverage,
    sample_idx: usize,
    prepass: &ParalogPrePass,
    inbreeding: &[f64],
) -> Option<SampleObservation> {
    let model = prepass.samples().get(sample_idx)?.as_ref()?;
    let window_gc = *window.gc.get(sample_idx)?;
    let mean_depth = *window.coverage.get(sample_idx)?;
    // A NaN sentinel (sample absent at the locus) is not a usable observation.
    if !window_gc.is_finite() || !mean_depth.is_finite() {
        return None;
    }
    let rel = model
        .coverage_model
        .relative_copy_number(f64::from(window_gc), f64::from(mean_depth));
    // Biallelic AD from the forwarded scalars: [REF, ALT] `num_obs`.
    let row = record.scalars_row(sample_idx);
    let ref_obs = row.first()?.num_obs;
    let alt_obs = row.get(1)?.num_obs;
    let total = ref_obs.checked_add(alt_obs)?;
    if total == 0 {
        return None;
    }
    Some(SampleObservation {
        relative_copy_number: rel,
        alt_reads: alt_obs,
        total_reads: total,
        // Fallible like every sibling read above, so this free function has no
        // length precondition of its own (the caller guards `inbreeding.len()
        // == n`, but matching the `?` style removes the lone panic vector).
        inbreeding_coefficient: *inbreeding.get(sample_idx)?,
    })
}

/// Score one locus for paralogy from its caller record and joined window
/// coverage, returning its likelihood ratio, or `None` if the locus cannot be
/// scored at all (not a biallelic SNP, a cohort-size mismatch, or *no* usable
/// samples — all "keep, never flag" cases).
///
/// There is **no minimum-sample count**: an under-powered locus is scored and
/// self-gates in the LR (LR≈0 → kept). Only a locus with zero usable samples is
/// left unscored — there is nothing to fold into the histogram.
///
/// `window` is the locus's per-sample centred-window coverage, read straight
/// from the record spill — each sample's own GC + mean coverage (`NaN` where a
/// sample had no covered record there).
///
/// **Pure and deterministic**, so every locus's LR is reproducible from its own
/// inputs — the invariant that lets the sink fold, and the write pass apply, the
/// exact same value. `obs_buf` is a caller scratch buffer (reused across loci to
/// keep the per-locus cost allocation-free and memory `O(samples)`).
///
/// Private to this module: production reaches it only through
/// [`ParalogScoringContext::score`]; the `#[cfg(test)]` [`calibrate`] /
/// [`score_spill_record`] and the unit tests are its only other callers.
#[allow(clippy::too_many_arguments)]
fn score_spilled_locus(
    record: &PosteriorRecord,
    window: &LocusWindowCoverage,
    prepass: &ParalogPrePass,
    inbreeding: &[f64],
    single_copy_depth_sd: &[f64],
    precompute: &ParalogScorePrecompute,
    obs_buf: &mut Vec<Option<SampleObservation>>,
) -> Option<f64> {
    if !is_biallelic_snp(record) {
        return None;
    }
    let n = record.n_samples;
    // The σ₀ slice must line up with the cohort, or the scorer would silently
    // truncate via its `zip`; a mismatch is a hard skip, not a wrong score.
    if single_copy_depth_sd.len() != n || inbreeding.len() != n {
        return None;
    }

    obs_buf.clear();
    obs_buf.resize(n, None);
    for (sample_idx, slot) in obs_buf.iter_mut().enumerate() {
        *slot = build_observation(record, window, sample_idx, prepass, inbreeding);
    }
    // No usable sample → nothing to score (a locus with all-`None` observations
    // would fold a meaningless LR≈0 into the histogram); one or more → scored,
    // and the LR self-gates if the coverage can't discriminate.
    if obs_buf.iter().all(|o| o.is_none()) {
        return None;
    }

    let score = score_locus_for_paralogy(
        &LocusObservations { samples: obs_buf },
        single_copy_depth_sd,
        precompute,
    );
    Some(score.paralog_log_likelihood_ratio)
}

/// Score one record-spill record from its own carried window coverage.
///
/// **Test-only after M7**: production scores each locus once in the worker
/// ([`ParalogScoringContext::score`]) and stores the LR on the spill, so no read
/// pass re-scores. Kept to drive the calibrate tests directly from spill records.
/// `None` — unscored, so **kept** — when the scorer declines (not a SNP, no
/// usable samples).
#[cfg(test)]
pub(crate) fn score_spill_record(
    spill: &ParalogSpillRecord,
    prepass: &ParalogPrePass,
    inbreeding: &[f64],
    single_copy_depth_sd: &[f64],
    precompute: &ParalogScorePrecompute,
    obs_buf: &mut Vec<Option<SampleObservation>>,
) -> Option<f64> {
    score_spilled_locus(
        &spill.record,
        &spill.window_coverage,
        prepass,
        inbreeding,
        single_copy_depth_sd,
        precompute,
        obs_buf,
    )
}

/// The cohort-constant inputs the hidden-paralog score needs, built **once**
/// when the filter is on and shared read-only across the caller workers.
///
/// After the single-individual reformulation the LR depends on nothing
/// genome-wide (arch `hidden_paralog_inline_scoring.md`), so a worker can score
/// a locus the moment its own data exists — it needs only this context (the
/// pre-pass coverage models, the cohort `F` slice, the σ₀ slice, and the
/// locus-invariant precompute) plus the locus's record + window coverage + a
/// reusable scratch buffer. Owning these here lets [`score`](Self::score) be a
/// method on shared `&self`, keeping [`crate::var_calling::VariantCaller`]
/// `Sync`.
pub(crate) struct ParalogScoringContext {
    prepass: ParalogPrePass,
    inbreeding: Vec<f64>,
    single_copy_depth_sd: Vec<f64>,
    precompute: ParalogScorePrecompute,
}

impl ParalogScoringContext {
    /// Build the context from the fitted pre-pass, the cohort inbreeding
    /// coefficient `F` (the caller's `--inbreeding-coefficient`, applied to every
    /// sample), and the model params — the same quantities the retired
    /// calibrate/write scoring built per pass, now built once up front.
    pub(crate) fn new(
        prepass: ParalogPrePass,
        inbreeding_coefficient: f64,
        params: &ParalogModelParams,
    ) -> Self {
        let single_copy_depth_sd = prepass.single_copy_depth_sd();
        let inbreeding = cohort_inbreeding(single_copy_depth_sd.len(), inbreeding_coefficient);
        let precompute = ParalogScorePrecompute::new(params, &inbreeding);
        Self {
            prepass,
            inbreeding,
            single_copy_depth_sd,
            precompute,
        }
    }

    /// The hidden-paralog LR for one called locus from its record + per-sample
    /// window coverage, or `f64::NAN` when the locus is unscored (not a biallelic
    /// SNP, or no usable samples) — the "keep, never flag" cases. `NaN` (rather
    /// than `Option`) so the worker stores it straight into the per-record
    /// `paralog_lr` vector and every downstream step treats non-finite as
    /// "unscored → kept" uniformly (the histogram drops it, `flags` returns
    /// `false`). `obs_buf` is a per-worker scratch buffer reused across the
    /// chunk's records.
    pub(crate) fn score(
        &self,
        record: &PosteriorRecord,
        window: &LocusWindowCoverage,
        obs_buf: &mut Vec<Option<SampleObservation>>,
    ) -> f64 {
        score_spilled_locus(
            record,
            window,
            &self.prepass,
            &self.inbreeding,
            &self.single_copy_depth_sd,
            &self.precompute,
            obs_buf,
        )
        .unwrap_or(f64::NAN)
    }
}

/// Estimate `π` and the FDR cut from an **already-folded** LR histogram — the
/// histogram the sink builds inline during the main pass (arch
/// `hidden_paralog_inline_scoring.md`), so calibration needs no spill read and
/// no scoring. Pure: same histogram → same calibration.
///
/// On EM non-convergence the prior falls back to `cfg.fallback_prior` (the
/// caller inspects `prior.converged` to warn).
pub(crate) fn calibrate_from_histogram(
    histogram: &ParalogLrHistogram,
    fdr_target: f64,
    cfg: &CalibrationConfig,
) -> ParalogCalibration {
    let estimated = ParalogPrior::estimate(histogram, &cfg.em);
    let prior = if estimated.converged {
        estimated
    } else {
        // Fallback: the EM did not converge (or the histogram was empty), so
        // trust a documented default π over a junk estimate. `converged =
        // false` signals the caller to warn.
        ParalogPrior {
            prior_probability: cfg.fallback_prior,
            converged: false,
        }
    };
    let curve = ParalogFdrCurve::from_histogram(histogram, &prior);
    let lr_threshold = curve.lr_threshold_for_fdr(fdr_target);
    ParalogCalibration {
        prior,
        curve,
        lr_threshold,
        target_fdr: fdr_target,
    }
}

/// Stream the spill once, scoring every locus and folding its LR into a
/// fixed-size histogram, then estimate `π` and the FDR cut — all RAM-flat (the
/// histogram is a few KB; no per-locus vector).
///
/// **Test-only after M7**: production folds the histogram inline in the main
/// pass ([`ParalogScoringContext::score`] in the worker + the sink) and
/// calibrates via [`calibrate_from_histogram`]. This spill-reading path is kept
/// to exercise the record → LR → π/cut chain end-to-end from record fixtures.
///
/// `inbreeding` is the cohort inbreeding coefficient `F` (the caller's
/// `--inbreeding-coefficient`), applied to every sample; `fdr_target` is the
/// operator's FDR.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn calibrate<R: Read>(
    prepass: &ParalogPrePass,
    inbreeding_coefficient: f64,
    records: &mut ParalogSpillReader<R>,
    params: &ParalogModelParams,
    fdr_target: f64,
    cfg: &CalibrationConfig,
) -> Result<ParalogCalibration, SpillError> {
    let single_copy_depth_sd = prepass.single_copy_depth_sd();
    let inbreeding = cohort_inbreeding(single_copy_depth_sd.len(), inbreeding_coefficient);
    let precompute = ParalogScorePrecompute::new(params, &inbreeding);

    let mut histogram = ParalogLrHistogram::with_defaults();
    let mut obs_buf: Vec<Option<SampleObservation>> = Vec::new();
    while let Some(record) = records.next_record() {
        let record = record?;
        if let Some(lr) = score_spill_record(
            &record,
            prepass,
            &inbreeding,
            &single_copy_depth_sd,
            &precompute,
            &mut obs_buf,
        ) {
            histogram.push(lr);
        }
        // `record` is dropped here — nothing per-locus is retained.
    }

    Ok(calibrate_from_histogram(&histogram, fdr_target, cfg))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::paralog_filter::spill::ParalogSpillRecord;
    use crate::var_calling::paralog_filter::test_support::{
        allele, normal_locus, paralog_locus, prepass, window_coverage, write_spill,
    };
    use std::io::Cursor;

    /// Run [`calibrate`] over the record fixtures (each carries its own window
    /// coverage — no sibling window spill).
    fn run_calibrate(
        prepass: &ParalogPrePass,
        f: f64,
        records: &[ParalogSpillRecord],
        fdr: f64,
        cfg: &CalibrationConfig,
    ) -> ParalogCalibration {
        let mut reader = ParalogSpillReader::new(Cursor::new(write_spill(records)));
        calibrate(
            prepass,
            f,
            &mut reader,
            &ParalogModelParams::default(),
            fdr,
            cfg,
        )
        .expect("calibrate")
    }

    /// Score one fixture locus directly from its record-carried window coverage.
    fn score_fixture(
        record: &ParalogSpillRecord,
        prepass: &ParalogPrePass,
        inbreeding: &[f64],
        sigma0: &[f64],
        buf: &mut Vec<Option<SampleObservation>>,
    ) -> Option<f64> {
        let precompute = ParalogScorePrecompute::new(&ParalogModelParams::default(), inbreeding);
        score_spilled_locus(
            &record.record,
            &record.window_coverage,
            prepass,
            inbreeding,
            sigma0,
            &precompute,
            buf,
        )
    }

    /// Calibration separates paralog-like loci (positive LR, flagged) from
    /// normal loci (negative LR, kept): π lands between 0 and 1, and the FDR cut
    /// flags the paralog-like loci but not the normal ones.
    #[test]
    fn calibrate_separates_paralog_from_normal() {
        let n = 30;
        let prepass = prepass(n);
        let mut records = Vec::new();
        for i in 0..90u32 {
            records.push(normal_locus(1000 + i, n));
        }
        for i in 0..10u32 {
            records.push(paralog_locus(5000 + i, n));
        }

        let cfg = CalibrationConfig::default();
        // Hexp → F ≈ 0 (obs_het 0.005, so F ≈ 0.75; still outbred-ish).
        let cal = run_calibrate(&prepass, 0.02, &records, 0.05, &cfg);

        assert!(
            cal.prior.prior_probability > 0.0 && cal.prior.prior_probability < 1.0,
            "π = {}",
            cal.prior.prior_probability
        );
        // The paralog-like locus is flagged; the normal locus is not.
        let inbreeding = cohort_inbreeding(prepass.samples().len(), 0.02);
        let sigma0 = prepass.single_copy_depth_sd();
        let mut buf = Vec::new();
        let lr_paralog = score_fixture(
            &paralog_locus(5000, n),
            &prepass,
            &inbreeding,
            &sigma0,
            &mut buf,
        )
        .expect("paralog scored");
        let lr_normal = score_fixture(
            &normal_locus(1000, n),
            &prepass,
            &inbreeding,
            &sigma0,
            &mut buf,
        )
        .expect("normal scored");
        assert!(
            lr_paralog > lr_normal,
            "paralog LR {lr_paralog} > normal LR {lr_normal}"
        );
        assert!(
            cal.flags(lr_paralog),
            "paralog-like locus should be flagged"
        );
        assert!(!cal.flags(lr_normal), "normal locus should be kept");
    }

    /// Calibrating a spill with no scorable loci (all indels) leaves the EM
    /// unconverged, so the prior takes the documented fallback.
    #[test]
    fn empty_histogram_takes_fallback_prior() {
        let n = 10;
        let prepass = prepass(n);
        // An indel record (REF "AT") — not a biallelic SNP, never scored.
        let mut indel = normal_locus(100, n);
        indel.record.alleles = vec![allele(b"AT"), allele(b"A")];
        let cal = run_calibrate(
            &prepass,
            0.02,
            &[indel],
            0.01,
            &CalibrationConfig::default(),
        );
        assert!(!cal.prior.converged, "empty histogram → non-converged");
        assert_eq!(cal.prior.prior_probability, DEFAULT_FALLBACK_PARALOG_PRIOR);
    }

    /// `score_spilled_locus` returns `None` only when there is nothing to score —
    /// a non-SNP locus or a locus with no usable samples — never on a
    /// sample-count threshold. A low-n locus with real coverage is still scored
    /// (the LR self-gates in place of the old `min_samples` gate).
    #[test]
    fn unscored_loci_are_none() {
        let n = 4;
        let prepass = prepass(n);
        let inbreeding = cohort_inbreeding(prepass.samples().len(), 0.02);
        let sigma0 = prepass.single_copy_depth_sd();
        let mut buf = Vec::new();

        // No usable sample (every sample's window coverage is `NaN`) → nothing
        // to fold.
        let mut no_depth = normal_locus(1, n);
        no_depth.window_coverage = window_coverage(n, 0.5, &|_| None);
        assert!(score_fixture(&no_depth, &prepass, &inbreeding, &sigma0, &mut buf).is_none());
        // Not a biallelic SNP.
        let mut indel = normal_locus(1, n);
        indel.record.alleles = vec![allele(b"AT"), allele(b"A")];
        assert!(score_fixture(&indel, &prepass, &inbreeding, &sigma0, &mut buf).is_none());
        // A tiny cohort (n=4) is still scored — no count gate.
        assert!(
            score_fixture(
                &normal_locus(1, n),
                &prepass,
                &inbreeding,
                &sigma0,
                &mut buf
            )
            .is_some(),
            "a low-n locus is scored; the LR self-gates instead of a count gate"
        );
    }

    /// A record carrying no usable window coverage (every sample `NaN`) is left
    /// unscored (kept) via [`score_spill_record`] — the record-spill analogue of
    /// the retired "no window for the tile" join miss.
    #[test]
    fn locus_without_window_coverage_is_unscored() {
        let n = 10;
        let prepass = prepass(n);
        let inbreeding = cohort_inbreeding(prepass.samples().len(), 0.02);
        let sigma0 = prepass.single_copy_depth_sd();
        let precompute = ParalogScorePrecompute::new(&ParalogModelParams::default(), &inbreeding);
        let mut buf = Vec::new();
        let mut rec = paralog_locus(5000, n);
        rec.window_coverage = window_coverage(n, 0.5, &|_| None);
        let scored =
            score_spill_record(&rec, &prepass, &inbreeding, &sigma0, &precompute, &mut buf);
        assert!(scored.is_none(), "no usable window coverage → unscored");
    }

    /// `build_observation` skips a sample whose window GC is `NaN` even when its
    /// coverage is finite (and vice-versa). The `window_coverage` fixtures only
    /// NaN the coverage (GC stays 0.5), so this pins the `gc.is_finite()` half of
    /// the guard, which no other test exercises.
    #[test]
    fn build_observation_skips_on_half_nan_window() {
        let n = 1;
        let prepass = prepass(n);
        let inbreeding = cohort_inbreeding(1, 0.02);
        let rec = normal_locus(1, n);

        // NaN GC, finite coverage → skipped.
        let nan_gc = LocusWindowCoverage {
            gc: vec![f32::NAN],
            coverage: vec![20.0],
        };
        assert!(
            build_observation(&rec.record, &nan_gc, 0, &prepass, &inbreeding).is_none(),
            "a NaN window GC must skip the sample"
        );

        // Finite GC, NaN coverage → skipped.
        let nan_cov = LocusWindowCoverage {
            gc: vec![0.5],
            coverage: vec![f32::NAN],
        };
        assert!(
            build_observation(&rec.record, &nan_cov, 0, &prepass, &inbreeding).is_none(),
            "a NaN window coverage must skip the sample"
        );

        // Both finite → a usable observation.
        let finite = LocusWindowCoverage {
            gc: vec![0.5],
            coverage: vec![20.0],
        };
        assert!(
            build_observation(&rec.record, &finite, 0, &prepass, &inbreeding).is_some(),
            "a fully finite window is a usable observation"
        );
    }

    /// A non-finite LR is never flagged — it is screened before the curve
    /// lookup, matching the histogram's refusal to fold a non-finite LR (so the
    /// two passes treat garbage identically: excluded, kept).
    #[test]
    fn flags_returns_false_for_non_finite_lr() {
        let n = 10;
        let prepass = prepass(n);
        // A calibration with a low target so finite high-LR loci would flag.
        let records: Vec<ParalogSpillRecord> = (0..20).map(|i| paralog_locus(i + 1, n)).collect();
        let cal = run_calibrate(&prepass, 0.02, &records, 0.5, &CalibrationConfig::default());
        assert!(!cal.flags(f64::NAN));
        assert!(!cal.flags(f64::INFINITY));
        assert!(!cal.flags(f64::NEG_INFINITY));
    }

    /// A higher target FDR flags at least as many loci as a lower one (the cut
    /// relaxes) — `flags` is monotone in the target.
    #[test]
    fn higher_fdr_flags_more() {
        let n = 20;
        let prepass = prepass(n);
        let mut records = Vec::new();
        for i in 0..80u32 {
            records.push(normal_locus(1000 + i, n));
        }
        for i in 0..20u32 {
            records.push(paralog_locus(5000 + i, n));
        }
        let cfg = CalibrationConfig::default();

        let cal_strict = run_calibrate(&prepass, 0.02, &records, 0.01, &cfg);
        let cal_loose = run_calibrate(&prepass, 0.02, &records, 0.5, &cfg);

        // The loose threshold is at or below the strict one (less stringent).
        match (cal_strict.lr_threshold, cal_loose.lr_threshold) {
            (Some(strict), Some(loose)) => {
                assert!(
                    loose <= strict,
                    "loose cut {loose} should be <= strict {strict}"
                )
            }
            (None, Some(_)) => {} // strict unachievable, loose achievable — fine
            other => panic!("unexpected thresholds: {other:?}"),
        }
    }
}
