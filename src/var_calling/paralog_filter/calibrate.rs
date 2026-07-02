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
    ParalogPrior, ParalogScorePrecompute, SampleObservation, score_locus_for_paralogy,
};

use crate::var_calling::posterior_engine::PosteriorRecord;

use super::prepass::ParalogPrePass;
use super::spill::{ParalogSpillReader, SpillError, WindowJoin};

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
/// `window_gc` and `window_depths` are the locus's joined window coverage (from
/// the sibling window spill): the shared reference GC and each sample's window
/// mean depth. Together with the record's per-sample AD they form the score's
/// per-sample observation.
fn build_observation(
    record: &PosteriorRecord,
    window_gc: f32,
    window_depths: &[Option<f32>],
    sample_idx: usize,
    prepass: &ParalogPrePass,
    inbreeding: &[f64],
) -> Option<SampleObservation> {
    let model = prepass.samples().get(sample_idx)?.as_ref()?;
    let mean_depth = (*window_depths.get(sample_idx)?)?;
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
/// `window_gc` + `window_depths` are the locus's window coverage joined from the
/// window spill by tile key: the window's shared reference GC and each sample's
/// window mean depth (`None` where a sample had no covered window there).
///
/// **Pure and deterministic**, so the calibrate pass (S4) and the write pass
/// (S5) get bit-identical LRs from the same joined inputs — the invariant that
/// lets the cut apply to the histogram it was built from. `obs_buf` is a caller
/// scratch buffer (reused across loci to keep the per-locus cost allocation-free
/// and memory `O(samples)`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn score_spilled_locus(
    record: &PosteriorRecord,
    window_gc: f32,
    window_depths: &[Option<f32>],
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
        *slot = build_observation(record, window_gc, window_depths, sample_idx, prepass, inbreeding);
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

/// Score one record-spill record by joining its window from `windows` (by tile
/// key `(chrom_id, (pos−1)/window_bp)`) and running the pure scorer. Shared by
/// the calibrate pass (S4) and the write pass (S5) so their LRs are
/// bit-identical. `Ok(None)` — unscored, so **kept** — when the window spill
/// holds no window for the locus's tile (the producer never marked it variant)
/// or the scorer declines (not a SNP, no usable samples). Advances `windows` past
/// any window sorting before the locus; the locus stream must be monotone in
/// `(chrom_id, tile)` (it is — the record spill is genomic-ordered).
#[allow(clippy::too_many_arguments)]
pub(crate) fn score_joined_locus<WR: Read>(
    record: &PosteriorRecord,
    windows: &mut WindowJoin<WR>,
    window_bp: u32,
    prepass: &ParalogPrePass,
    inbreeding: &[f64],
    single_copy_depth_sd: &[f64],
    precompute: &ParalogScorePrecompute,
    obs_buf: &mut Vec<Option<SampleObservation>>,
) -> Result<Option<f64>, SpillError> {
    let tile = record.locus.start.saturating_sub(1) / window_bp.max(1);
    let Some(window) = windows.window_for(record.locus.chrom_id, tile)? else {
        return Ok(None);
    };
    Ok(score_spilled_locus(
        record,
        window.gc,
        &window.depths,
        prepass,
        inbreeding,
        single_copy_depth_sd,
        precompute,
        obs_buf,
    ))
}

/// Stream the spill once, scoring every locus and folding its LR into a
/// fixed-size histogram, then estimate `π` and the FDR cut — all RAM-flat (the
/// histogram is a few KB; no per-locus vector).
///
/// `inbreeding` is the cohort inbreeding coefficient `F` (the caller's
/// `--inbreeding-coefficient`), applied to every sample; `fdr_target` is the
/// operator's FDR. On EM non-convergence the prior falls back to
/// `cfg.fallback_prior` (the caller inspects `prior.converged` to warn).
#[allow(clippy::too_many_arguments)]
pub(crate) fn calibrate<R: Read, WR: Read>(
    prepass: &ParalogPrePass,
    inbreeding_coefficient: f64,
    records: &mut ParalogSpillReader<R>,
    windows: &mut WindowJoin<WR>,
    window_bp: u32,
    params: &ParalogModelParams,
    fdr_target: f64,
    cfg: &CalibrationConfig,
) -> Result<ParalogCalibration, SpillError> {
    let single_copy_depth_sd = prepass.single_copy_depth_sd();
    let inbreeding = cohort_inbreeding(single_copy_depth_sd.len(), inbreeding_coefficient);
    // The per-pass locus-invariant scoring tables (Wright genotype priors +
    // carrier-frequency probs + configs), built once from the params + cohort
    // inbreeding and reused for every locus. The write pass builds an identical
    // one, so the recomputed LRs stay bit-identical.
    let precompute = ParalogScorePrecompute::new(params, &inbreeding);

    let mut histogram = ParalogLrHistogram::with_defaults();
    let mut obs_buf: Vec<Option<SampleObservation>> = Vec::new();
    while let Some(record) = records.next_record() {
        let record = record?;
        if let Some(lr) = score_joined_locus(
            &record.record,
            windows,
            window_bp,
            prepass,
            &inbreeding,
            &single_copy_depth_sd,
            &precompute,
            &mut obs_buf,
        )? {
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
    use crate::var_calling::paralog_filter::spill::{
        ParalogSpillRecord, WindowJoin, WindowSpillReader, WindowSpillRecord,
    };
    use crate::var_calling::paralog_filter::test_support::{
        TEST_WINDOW_BP, allele, normal_locus, normal_window, paralog_locus, paralog_window, prepass,
        snp_window, write_spill, write_window_spill,
    };
    use std::io::Cursor;

    /// Run [`calibrate`] over parallel record / window fixtures (the record spill
    /// and its sibling window spill, joined by tile key).
    fn run_calibrate(
        prepass: &ParalogPrePass,
        f: f64,
        records: &[ParalogSpillRecord],
        windows: &[WindowSpillRecord],
        fdr: f64,
        cfg: &CalibrationConfig,
    ) -> ParalogCalibration {
        let mut reader = ParalogSpillReader::new(Cursor::new(write_spill(records)));
        let mut join =
            WindowJoin::new(WindowSpillReader::new(Cursor::new(write_window_spill(windows))))
                .unwrap();
        calibrate(
            prepass,
            f,
            &mut reader,
            &mut join,
            TEST_WINDOW_BP,
            &ParalogModelParams::default(),
            fdr,
            cfg,
        )
        .expect("calibrate")
    }

    /// Score one fixture locus directly (its record + its window).
    fn score_fixture(
        record: &ParalogSpillRecord,
        window: &WindowSpillRecord,
        prepass: &ParalogPrePass,
        inbreeding: &[f64],
        sigma0: &[f64],
        buf: &mut Vec<Option<SampleObservation>>,
    ) -> Option<f64> {
        let precompute = ParalogScorePrecompute::new(&ParalogModelParams::default(), inbreeding);
        score_spilled_locus(
            &record.record,
            window.gc,
            &window.depths,
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
        let mut windows = Vec::new();
        for i in 0..90u32 {
            records.push(normal_locus(1000 + i, n));
            windows.push(normal_window(1000 + i, n));
        }
        for i in 0..10u32 {
            records.push(paralog_locus(5000 + i, n));
            windows.push(paralog_window(5000 + i, n));
        }

        let cfg = CalibrationConfig::default();
        // Hexp → F ≈ 0 (obs_het 0.005, so F ≈ 0.75; still outbred-ish).
        let cal = run_calibrate(&prepass, 0.02, &records, &windows, 0.05, &cfg);

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
            &paralog_window(5000, n),
            &prepass,
            &inbreeding,
            &sigma0,
            &mut buf,
        )
        .expect("paralog scored");
        let lr_normal = score_fixture(
            &normal_locus(1000, n),
            &normal_window(1000, n),
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
        assert!(cal.flags(lr_paralog), "paralog-like locus should be flagged");
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
            &[normal_window(100, n)],
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

        // No usable sample (every window depth is `None`) → nothing to fold.
        let no_depth = snp_window(1, n, 0.5, &|_| None);
        assert!(
            score_fixture(
                &normal_locus(1, n),
                &no_depth,
                &prepass,
                &inbreeding,
                &sigma0,
                &mut buf
            )
            .is_none()
        );
        // Not a biallelic SNP.
        let mut indel = normal_locus(1, n);
        indel.record.alleles = vec![allele(b"AT"), allele(b"A")];
        assert!(
            score_fixture(
                &indel,
                &normal_window(1, n),
                &prepass,
                &inbreeding,
                &sigma0,
                &mut buf
            )
            .is_none()
        );
        // A tiny cohort (n=4) is still scored — no count gate.
        assert!(
            score_fixture(
                &normal_locus(1, n),
                &normal_window(1, n),
                &prepass,
                &inbreeding,
                &sigma0,
                &mut buf
            )
            .is_some(),
            "a low-n locus is scored; the LR self-gates instead of a count gate"
        );
    }

    /// A record whose tile has no window in the window spill is left unscored
    /// (kept), not scored on stale coverage — the join returns `None`.
    #[test]
    fn locus_without_a_window_is_unscored() {
        let n = 10;
        let prepass = prepass(n);
        let inbreeding = cohort_inbreeding(prepass.samples().len(), 0.02);
        let sigma0 = prepass.single_copy_depth_sd();
        let precompute = ParalogScorePrecompute::new(&ParalogModelParams::default(), &inbreeding);
        let mut buf = Vec::new();
        // Window spill holds only tile for pos 100; the record at pos 5000 has no
        // window → the join yields None → unscored.
        let mut join = WindowJoin::new(WindowSpillReader::new(Cursor::new(write_window_spill(&[
            paralog_window(100, n),
        ]))))
        .unwrap();
        let rec = paralog_locus(5000, n);
        let scored = score_joined_locus(
            &rec.record,
            &mut join,
            TEST_WINDOW_BP,
            &prepass,
            &inbreeding,
            &sigma0,
            &precompute,
            &mut buf,
        )
        .unwrap();
        assert!(scored.is_none(), "no window for the tile → unscored");
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
        let windows: Vec<WindowSpillRecord> = (0..20).map(|i| paralog_window(i + 1, n)).collect();
        let cal = run_calibrate(
            &prepass,
            0.02,
            &records,
            &windows,
            0.5,
            &CalibrationConfig::default(),
        );
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
        let mut windows = Vec::new();
        for i in 0..80u32 {
            records.push(normal_locus(1000 + i, n));
            windows.push(normal_window(1000 + i, n));
        }
        for i in 0..20u32 {
            records.push(paralog_locus(5000 + i, n));
            windows.push(paralog_window(5000 + i, n));
        }
        let cfg = CalibrationConfig::default();

        let cal_strict = run_calibrate(&prepass, 0.02, &records, &windows, 0.01, &cfg);
        let cal_loose = run_calibrate(&prepass, 0.02, &records, &windows, 0.5, &cfg);

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
