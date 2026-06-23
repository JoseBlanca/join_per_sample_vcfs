//! The per-locus EM seeds `π⁰` and `θ⁰` (arch `ssr_call_genotyping.md` §5, spec
//! §4.3/§5.5, implementation plan §3, C3).
//!
//! Computed **in Phase 3** at each locus's EM start (Q-P3): the pre-pass runs on a
//! subset and cannot seed loci it never visits, so it supplies only the *global*
//! params (frozen `ε`, shape, level, `G₀` decay) and the per-locus seeds are built
//! here off the rungs + candidates.
//!
//! - `π⁰` = normalized(`putative-genotype counts` + `G₀` pseudocounts). Every sample
//!   contributes exactly `ploidy` allele-copies off its clear peaks (a homozygote
//!   repeats its top peak); with no confident sample anywhere, `π⁰` falls back to
//!   the normalized `G₀` (the conservative small-N behaviour).
//! - `θ⁰` = the locus period's stutter shape from the supplied params.

use crate::ssr::cohort::allele_freq_prior::g0_pseudocounts;
use crate::ssr::cohort::candidate_set::CandidateSet;
use crate::ssr::cohort::param_estimation::{ParamSet, StutterShape};
use crate::ssr::cohort::rung_ladder::{Rungs, sample_clear_peaks};
use crate::ssr::cohort::types::CohortLocus;

/// The EM starting point for one locus.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct LocusSeed {
    /// Initial candidate allele frequencies (normalized, parallel to candidates).
    pub(crate) pi0: Vec<f64>,
    /// Initial stutter shape for the locus (the period's supplied shape).
    pub(crate) theta0: StutterShape,
}

/// A default shape if the params carry none for the period (a degenerate fallback —
/// the pre-pass normally supplies a parent for every period present).
const FALLBACK_SHAPE: StutterShape = StutterShape {
    up_rate: 1.0,
    down_rate: 1.0,
    decay: 0.1,
};

/// The first candidate index whose allele has `length` repeat units, if any.
///
/// SIMPLIFICATION: a rung can legitimately hold several same-length sequences
/// (substitution / interruption variants, first-class per spec §5/§7); this maps a
/// nominated length to the *first* such candidate. Exact for the pure tracts that
/// dominate; multi-variant-per-length attribution is a follow-up tied to S2 impure
/// alleles.
fn candidate_of_length(candidates: &CandidateSet, period: usize, length: u16) -> Option<usize> {
    candidates
        .alleles
        .iter()
        .position(|seq| (seq.len() / period) as u16 == length)
}

/// Build the per-locus EM seed.
pub(crate) fn seed_locus(
    locus: &CohortLocus,
    rungs: &Rungs,
    candidates: &CandidateSet,
    params: &ParamSet,
    ploidy: u8,
    prominence: u32,
) -> LocusSeed {
    let period = rungs.period();

    // Start the tally from the G₀ pseudocounts (so every candidate stays > 0).
    let mut counts = g0_pseudocounts(candidates, rungs, period_decay(params, period));

    for evidence in &locus.samples {
        let mut peaks = sample_clear_peaks(evidence, period, prominence);
        if peaks.is_empty() {
            continue; // no putative genotype for this sample → prior only
        }
        peaks.sort_by(|a, b| {
            b.support
                .cmp(&a.support)
                .then_with(|| a.repeat_len.cmp(&b.repeat_len))
        });

        // Exactly `ploidy` allele-copies: the top peaks get one each; any shortfall
        // (an under-resolved / homozygous sample) is absorbed by the top peak.
        let resolved = peaks.len().min(ploidy as usize);
        let mut copies_left = ploidy as usize;
        for peak in &peaks[..resolved] {
            if let Some(idx) = candidate_of_length(candidates, period, peak.repeat_len) {
                counts[idx] += 1.0;
                copies_left -= 1;
            }
        }
        if copies_left > 0
            && let Some(idx) = candidate_of_length(candidates, period, peaks[0].repeat_len)
        {
            counts[idx] += copies_left as f64;
        }
    }

    let total: f64 = counts.iter().sum();
    let pi0 = counts.iter().map(|c| c / total).collect();

    let theta0 = params
        .stutter_shape_parent
        .get(&(period as u8))
        .copied()
        .unwrap_or(FALLBACK_SHAPE);

    LocusSeed { pi0, theta0 }
}

/// The `G₀` decay for the locus period, or a mild default if the params carry none.
fn period_decay(
    params: &ParamSet,
    period: usize,
) -> &crate::ssr::cohort::param_estimation::G0PseudocountDecay {
    use crate::ssr::cohort::param_estimation::G0PseudocountDecay;
    const FALLBACK_DECAY: G0PseudocountDecay = G0PseudocountDecay { p: 0.5 };
    params
        .pseudocount_decay_per_loci_group
        .get(&(period as u8))
        .unwrap_or(&FALLBACK_DECAY)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::candidate_set::{Admission, CandidateCfg, assemble_candidates};
    use crate::ssr::cohort::param_estimation::{G0PseudocountDecay, SampleGroupId};
    use crate::ssr::cohort::rung_ladder::{RungCfg, build_rungs};
    use crate::ssr::cohort::types::{CohortLocus, LocusId, SampleEvidence, SsrQc};
    use crate::ssr::types::Motif;
    use std::collections::HashMap;

    fn ca_seq(units: u16) -> Box<[u8]> {
        std::iter::repeat_n(*b"CA", units as usize)
            .flatten()
            .collect()
    }

    fn ca_evidence(bins: &[(u16, u32)]) -> SampleEvidence {
        let mut seq_counts: Vec<(Box<[u8]>, u32)> =
            bins.iter().map(|&(u, c)| (ca_seq(u), c)).collect();
        seq_counts.sort_by(|(a, _), (b, _)| a.cmp(b));
        SampleEvidence {
            seq_counts,
            qc: SsrQc::default(),
        }
    }

    fn locus_with(ref_units: u16, samples: Vec<SampleEvidence>) -> CohortLocus {
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 20,
                end: 20 + 2 * ref_units as u32,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGG".as_slice()),
            ca_seq(ref_units),
        );
        for (idx, ev) in samples.into_iter().enumerate() {
            locus.push(idx as u32, ev);
        }
        locus
    }

    fn params_with_decay(p: f64) -> ParamSet {
        let mut decay = HashMap::new();
        decay.insert(2u8, G0PseudocountDecay { p });
        ParamSet {
            error_per_sample_group: vec![],
            stutter_shape_parent: HashMap::new(),
            stutter_shape_by_cell: HashMap::new(),
            level_seed: vec![],
            pseudocount_decay_per_loci_group: decay,
            group_of_sample: vec![SampleGroupId(0)],
            f0_seed: 0.0,
        }
    }

    fn seed(locus: &CohortLocus) -> LocusSeed {
        let rungs = build_rungs(locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(locus, &rungs, 2, &CandidateCfg::dev_default());
        seed_locus(locus, &rungs, &candidates, &params_with_decay(0.5), 2, 3)
    }

    #[test]
    fn pi0_is_a_normalized_distribution() {
        let sample = ca_evidence(&[(5, 5), (6, 50), (7, 5)]);
        let locus = locus_with(6, vec![sample.clone(), sample]);
        let s = seed(&locus);
        let sum: f64 = s.pi0.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
        assert!(s.pi0.iter().all(|&p| p > 0.0));
    }

    #[test]
    fn pi0_concentrates_on_the_observed_homozygous_allele() {
        // Two homozygous-6 samples → most of π⁰ on the length-6 candidate.
        let sample = ca_evidence(&[(5, 5), (6, 50), (7, 5)]);
        let locus = locus_with(6, vec![sample.clone(), sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        let s = seed_locus(&locus, &rungs, &candidates, &params_with_decay(0.5), 2, 3);

        let six = candidate_of_length(&candidates, 2, 6).unwrap();
        let top = (0..s.pi0.len())
            .max_by(|&a, &b| s.pi0[a].partial_cmp(&s.pi0[b]).unwrap())
            .unwrap();
        assert_eq!(top, six, "π⁰ should peak on the homozygous allele");
        assert!(s.pi0[six] > 0.5);
    }

    #[test]
    fn pi0_falls_back_to_normalized_g0_without_confident_samples() {
        // A single thin, structureless sample: no clear peak → π⁰ = normalized G₀.
        let sample = ca_evidence(&[(6, 2), (8, 2), (10, 2)]); // flat, no prominence
        let locus = locus_with(6, vec![sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = CandidateSet {
            alleles: vec![ca_seq(6), ca_seq(8)],
            ref_idx: 0,
            admit: Admission::Pass,
        };
        let s = seed_locus(&locus, &rungs, &candidates, &params_with_decay(0.5), 2, 3);
        // No genotype tally added, so π⁰ ∝ G₀ = [p^|Δ from mode|]; mode is 6/8/10 tie →
        // smallest (6). Candidate 6 (Δ=0) must carry at least as much as candidate 8.
        assert!(s.pi0[0] >= s.pi0[1]);
        let sum: f64 = s.pi0.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn theta0_uses_the_period_parent_shape_when_present() {
        let sample = ca_evidence(&[(6, 40)]);
        let locus = locus_with(6, vec![sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        let mut params = params_with_decay(0.5);
        let parent = StutterShape {
            up_rate: 0.2,
            down_rate: 0.8,
            decay: 0.3,
        };
        params.stutter_shape_parent.insert(2, parent);
        let s = seed_locus(&locus, &rungs, &candidates, &params, 2, 3);
        assert_eq!(s.theta0, parent);
    }
}
