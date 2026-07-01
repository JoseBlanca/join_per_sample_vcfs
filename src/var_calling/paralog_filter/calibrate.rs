//! S4 — the calibrate pass: spill → LR → histogram → π + FDR cut.
//!
//! The first of the two spill reads. It streams the spill locus by locus,
//! scores each candidate variant for paralogy, folds the likelihood ratio into
//! a **fixed-size histogram**, and discards the locus — nothing per-locus is
//! retained, so RAM stays flat in variant count (arch §3 step 5, §6). From the
//! finished histogram it estimates the empirical-Bayes prior `π` (Q5, with a
//! fixed-prior fallback if the EM does not converge), builds the tail-FDR
//! `q_of_lr` curve, and resolves the LR threshold for the operator's target
//! FDR.
//!
//! The per-locus scoring lives in [`score_spilled_locus`], a pure function of
//! the spilled inputs and the pre-pass. The **write pass (S5) calls the exact
//! same function** to recompute each LR, so the applied cut matches the
//! histogram it came from bit-for-bit (the non-negotiable invariant in arch §3).

use std::io::Read;

use crate::paralog::{
    EmConfig, LocusObservations, ParalogFdrCurve, ParalogLrHistogram, ParalogModelParams,
    ParalogPrior, SampleObservation, inbreeding_coefficient, score_locus_for_paralogy,
};

use super::prepass::ParalogPrePass;
use super::spill::{ParalogSpillReader, ParalogSpillRecord, SpillError};

/// Default minimum usable samples for a locus to be scored. A locus with fewer
/// usable samples has too little evidence for a trustworthy paralog verdict, so
/// it is left unscored (never flagged). Matches the R1 prototype's `m < 20`
/// skip.
pub(crate) const DEFAULT_MIN_SAMPLES_FOR_SCORE: usize = 20;

/// Default fallback `π` used when the EM does not converge — a rare-paralog
/// prior deliberately low so a pathological dataset yields a conservative
/// (few-drops) calibration rather than trusting a junk estimate. Chosen to
/// match [`crate::paralog::prior::DEFAULT_EM_START`] (`0.03`); they are
/// independent constants, so this comment is the only thing coupling them.
/// The operator is warned when this is used (S6).
pub(crate) const DEFAULT_FALLBACK_PARALOG_PRIOR: f64 = 0.03;

/// Tuning knobs for [`calibrate`]. [`Default`] is the production configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CalibrationConfig {
    /// Minimum usable samples for a locus to be scored (else left unscored).
    pub min_samples: usize,
    /// EM configuration for the prior estimate.
    pub em: EmConfig,
    /// Fallback `π` when the EM does not converge.
    pub fallback_prior: f64,
}

impl Default for CalibrationConfig {
    fn default() -> Self {
        Self {
            min_samples: DEFAULT_MIN_SAMPLES_FOR_SCORE,
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

/// Per-sample inbreeding coefficient `F = 1 − obs_het/Hexp`, indexed parallel to
/// the cohort columns; `0.0` (outbred, inert) for an absent sample.
///
/// The write pass (S5) recomputes this the same way so its `score_spilled_locus`
/// inputs are bit-identical to the calibrate pass's.
pub(crate) fn inbreeding_by_sample(prepass: &ParalogPrePass, hexp: f64) -> Vec<f64> {
    prepass
        .samples()
        .iter()
        .map(|s| match s {
            Some(model) => inbreeding_coefficient(model.obs_het, hexp),
            None => 0.0,
        })
        .collect()
}

/// Whether the spilled record is a biallelic SNP — the class the paralog score
/// is defined over (a single-base REF and a single-base ALT). Other loci
/// (indels, multiallelics) are not scored and are never flagged.
fn is_biallelic_snp(record: &crate::var_calling::posterior_engine::PosteriorRecord) -> bool {
    record.alleles.len() == 2
        && record.alleles[0].seq.len() == 1
        && record.alleles[1].seq.len() == 1
}

/// One sample's paralog observation at a spilled locus, or `None` if the sample
/// is unusable here (absent from the pre-pass, no covered window, or no reads).
fn build_observation(
    spill: &ParalogSpillRecord,
    sample_idx: usize,
    prepass: &ParalogPrePass,
    inbreeding: &[f64],
) -> Option<SampleObservation> {
    let model = prepass.samples().get(sample_idx)?.as_ref()?;
    let mean_depth = (*spill.window_mean_depth.get(sample_idx)?)?;
    let rel = model
        .coverage_model
        .relative_copy_number(f64::from(spill.window_gc), f64::from(mean_depth));
    // Biallelic AD from the forwarded scalars: [REF, ALT] `num_obs`.
    let row = spill.record.scalars_row(sample_idx);
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

/// Score one spilled locus for paralogy, returning its likelihood ratio, or
/// `None` if the locus is not scored (not a biallelic SNP, a cohort-size
/// mismatch, or fewer than `min_samples` usable samples — all "keep, never
/// flag" cases).
///
/// **Pure and deterministic**, so the calibrate pass (S4) and the write pass
/// (S5) get bit-identical LRs from the same spilled inputs — the invariant that
/// lets the cut apply to the histogram it was built from. `obs_buf` is a caller
/// scratch buffer (reused across loci to keep the per-locus cost allocation-free
/// and memory `O(samples)`).
pub(crate) fn score_spilled_locus(
    spill: &ParalogSpillRecord,
    prepass: &ParalogPrePass,
    inbreeding: &[f64],
    single_copy_depth_sd: &[f64],
    params: &ParalogModelParams,
    min_samples: usize,
    obs_buf: &mut Vec<Option<SampleObservation>>,
) -> Option<f64> {
    let record = &spill.record;
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
        *slot = build_observation(spill, sample_idx, prepass, inbreeding);
    }
    let usable = obs_buf.iter().filter(|o| o.is_some()).count();
    if usable < min_samples {
        return None;
    }

    let score = score_locus_for_paralogy(
        &LocusObservations { samples: obs_buf },
        single_copy_depth_sd,
        params,
    );
    Some(score.paralog_log_likelihood_ratio)
}

/// Stream the spill once, scoring every locus and folding its LR into a
/// fixed-size histogram, then estimate `π` and the FDR cut — all RAM-flat (the
/// histogram is a few KB; no per-locus vector).
///
/// `hexp` is the finished cohort expected heterozygosity (from the main pass's
/// [`super::prepass::HexpAccumulator`]); `fdr_target` is the operator's FDR. On
/// EM non-convergence the prior falls back to `cfg.fallback_prior` (the caller
/// inspects `prior.converged` to warn).
pub(crate) fn calibrate<R: Read>(
    prepass: &ParalogPrePass,
    hexp: f64,
    spill: &mut ParalogSpillReader<R>,
    params: &ParalogModelParams,
    fdr_target: f64,
    cfg: &CalibrationConfig,
) -> Result<ParalogCalibration, SpillError> {
    let inbreeding = inbreeding_by_sample(prepass, hexp);
    let single_copy_depth_sd = prepass.single_copy_depth_sd();

    let mut histogram = ParalogLrHistogram::with_defaults();
    let mut obs_buf: Vec<Option<SampleObservation>> = Vec::new();
    while let Some(record) = spill.next_record() {
        let record = record?;
        if let Some(lr) = score_spilled_locus(
            &record,
            prepass,
            &inbreeding,
            &single_copy_depth_sd,
            params,
            cfg.min_samples,
            &mut obs_buf,
        ) {
            histogram.push(lr);
        }
        // `record` is dropped here — nothing per-locus is retained.
    }

    let estimated = ParalogPrior::estimate(&histogram, &cfg.em);
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
    let curve = ParalogFdrCurve::from_histogram(&histogram, &prior);
    let lr_threshold = curve.lr_threshold_for_fdr(fdr_target);

    Ok(ParalogCalibration {
        prior,
        curve,
        lr_threshold,
        target_fdr: fdr_target,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::paralog_filter::test_support::{
        allele, normal_locus, paralog_locus, prepass, write_spill,
    };
    use std::io::Cursor;

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
        let bytes = write_spill(&records);

        let cfg = CalibrationConfig {
            min_samples: 5,
            ..Default::default()
        };
        let mut reader = ParalogSpillReader::new(Cursor::new(bytes));
        let cal = calibrate(
            &prepass,
            0.02, // Hexp → F ≈ 0 (obs_het 0.005, so F ≈ 0.75; still outbred-ish)
            &mut reader,
            &ParalogModelParams::default(),
            0.05,
            &cfg,
        )
        .expect("calibrate");

        assert!(
            cal.prior.prior_probability > 0.0 && cal.prior.prior_probability < 1.0,
            "π = {}",
            cal.prior.prior_probability
        );
        // The paralog-like locus is flagged; the normal locus is not.
        let inbreeding = inbreeding_by_sample(&prepass, 0.02);
        let sigma0 = prepass.single_copy_depth_sd();
        let mut buf = Vec::new();
        let lr_paralog = score_spilled_locus(
            &paralog_locus(5000, n),
            &prepass,
            &inbreeding,
            &sigma0,
            &ParalogModelParams::default(),
            5,
            &mut buf,
        )
        .expect("paralog scored");
        let lr_normal = score_spilled_locus(
            &normal_locus(1000, n),
            &prepass,
            &inbreeding,
            &sigma0,
            &ParalogModelParams::default(),
            5,
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
        let bytes = write_spill(&[indel]);

        let mut reader = ParalogSpillReader::new(Cursor::new(bytes));
        let cal = calibrate(
            &prepass,
            0.02,
            &mut reader,
            &ParalogModelParams::default(),
            0.01,
            &CalibrationConfig {
                min_samples: 5,
                ..Default::default()
            },
        )
        .expect("calibrate");
        assert!(!cal.prior.converged, "empty histogram → non-converged");
        assert_eq!(cal.prior.prior_probability, DEFAULT_FALLBACK_PARALOG_PRIOR);
    }

    /// `score_spilled_locus` returns `None` for a locus below `min_samples`
    /// (kept, never flagged) and for a non-SNP locus.
    #[test]
    fn unscored_loci_are_none() {
        let n = 4;
        let prepass = prepass(n);
        let inbreeding = inbreeding_by_sample(&prepass, 0.02);
        let sigma0 = prepass.single_copy_depth_sd();
        let params = ParalogModelParams::default();
        let mut buf = Vec::new();

        // Too few usable samples (require 100, only 4 exist).
        assert!(
            score_spilled_locus(
                &normal_locus(1, n),
                &prepass,
                &inbreeding,
                &sigma0,
                &params,
                100,
                &mut buf
            )
            .is_none()
        );
        // Not a biallelic SNP.
        let mut indel = normal_locus(1, n);
        indel.record.alleles = vec![allele(b"AT"), allele(b"A")];
        assert!(
            score_spilled_locus(&indel, &prepass, &inbreeding, &sigma0, &params, 1, &mut buf)
                .is_none()
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
        let records: Vec<ParalogSpillRecord> = (0..20).map(|i| paralog_locus(i, n)).collect();
        let bytes = write_spill(&records);
        let cal = calibrate(
            &prepass,
            0.02,
            &mut ParalogSpillReader::new(Cursor::new(bytes)),
            &ParalogModelParams::default(),
            0.5,
            &CalibrationConfig {
                min_samples: 5,
                ..Default::default()
            },
        )
        .unwrap();
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
        let bytes = write_spill(&records);
        let params = ParalogModelParams::default();
        let cfg = CalibrationConfig {
            min_samples: 5,
            ..Default::default()
        };

        let cal_strict = calibrate(
            &prepass,
            0.02,
            &mut ParalogSpillReader::new(Cursor::new(bytes.clone())),
            &params,
            0.01,
            &cfg,
        )
        .unwrap();
        let cal_loose = calibrate(
            &prepass,
            0.02,
            &mut ParalogSpillReader::new(Cursor::new(bytes)),
            &params,
            0.5,
            &cfg,
        )
        .unwrap();

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
