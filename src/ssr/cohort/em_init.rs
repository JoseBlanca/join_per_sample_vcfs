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
use crate::ssr::cohort::pair_hmm::{HmmScratch, align_subst};
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

/// The candidate a sample's peak seeds its allele-copies onto: the candidate whose sequence
/// **equals** the sample's representative `target` at `length` (so a same-length homozygote
/// seeds *its own* composition, not just the first candidate at that length — the P1.3 fix),
/// else — when that representative is not itself an admitted candidate — the closest
/// same-length candidate by `align_subst` (the substitution metric the read likelihood uses,
/// arch §2.4), ties broken on the lower candidate index. `None` iff no candidate sits at
/// `length`.
///
/// `eps` is evaluated **lazily**, only on the fallback path (a sample whose plurality
/// sequence is not an admitted candidate), so the common exact-match case pays nothing and
/// never touches per-sample chemistry.
///
/// SAME-LENGTH HET LIMITATION (spec §5.3): a same-length het shows one length peak with one
/// representative, so **both** allele copies seed the *majority* composition here (a hom
/// seed); the minority same-length allele starts only at its `G₀` floor. Recovering the het
/// is the EM's job (P1.4), not the seed's — documented, not fixed, in Phase 1.
fn candidate_for_sequence(
    candidates: &CandidateSet,
    target: &[u8],
    period: usize,
    length: u16,
    eps: impl FnOnce() -> f64,
    scratch: &mut HmmScratch,
) -> Option<usize> {
    if let Some(idx) = candidates
        .alleles
        .iter()
        .position(|seq| seq.as_ref() == target)
    {
        return Some(idx);
    }
    let eps = eps();
    candidates
        .alleles
        .iter()
        .enumerate()
        .filter(|(_, seq)| (seq.len() / period) as u16 == length)
        .map(|(idx, seq)| (idx, align_subst(target, seq, eps, scratch)))
        // Strict `>` keeps the first (lowest-index) maximum on a score tie (arch §2.4).
        .reduce(|best, cur| if cur.1 > best.1 { cur } else { best })
        .map(|(idx, _)| idx)
}

/// The frozen per-base error `ε` for present sample `k_present` (its sample group's value).
/// Reached only on the seed's rare fallback path, so it indexes **loudly** — matching
/// `em::sample_chemistry`'s decision-E contract — to surface a malformed `ParamSet` rather
/// than fabricate default chemistry.
fn sample_eps(locus: &CohortLocus, k_present: usize, params: &ParamSet) -> f64 {
    let global = locus.present[k_present] as usize;
    let group = *params
        .group_of_sample
        .get(global)
        .expect("decision E: every present sample has a frozen sample group");
    params
        .error_per_sample_group
        .get(group.0 as usize)
        .expect("every sample group has a frozen ε")
        .0
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

    // Reused across samples for the (rare) composition fallback in `candidate_for_sequence`.
    let mut scratch = HmmScratch::new();

    for (k_present, evidence) in locus.samples.iter().enumerate() {
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
        // (an under-resolved / homozygous sample) is absorbed by the top peak. Each peak
        // seeds the candidate matching its **representative sequence** (not just its length),
        // so two same-length candidates are seeded by composition (spec §5.3).
        let resolved = peaks.len().min(ploidy as usize);
        let mut copies_left = ploidy as usize;
        for peak in &peaks[..resolved] {
            if let Some(idx) = candidate_for_sequence(
                candidates,
                peak.allele.as_ref(),
                period,
                peak.repeat_len,
                || sample_eps(locus, k_present, params),
                &mut scratch,
            ) {
                counts[idx] += 1.0;
                copies_left -= 1;
            }
        }
        if copies_left > 0
            && let Some(idx) = candidate_for_sequence(
                candidates,
                peaks[0].allele.as_ref(),
                period,
                peaks[0].repeat_len,
                || sample_eps(locus, k_present, params),
                &mut scratch,
            )
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
    use crate::ssr::cohort::param_estimation::{
        G0PseudocountDecay, PerBaseError, PurityLevel, SampleGroupId,
    };
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
            purity_level: PurityLevel::none(),
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

        let six = candidates
            .alleles
            .iter()
            .position(|seq| seq.len() / 2 == 6)
            .unwrap();
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

    // ── P1.3: sequence-aware seed (same-length composition) ──

    /// `CA×6` (12 bp) — the pure same-length allele (candidate index 0 in these fixtures).
    fn pure6() -> Box<[u8]> {
        ca_seq(6)
    }

    /// A same-length (12 bp) interruption sibling of `pure6`: one interior base flipped.
    fn interrupted6() -> Box<[u8]> {
        let mut s = pure6().to_vec();
        s[5] = b'T';
        s.into_boxed_slice()
    }

    /// A sample observing the given `(sequence, count)` pairs (byte-sorted, zeros dropped).
    fn seq_sample(entries: &[(Box<[u8]>, u32)]) -> SampleEvidence {
        let mut seq_counts: Vec<(Box<[u8]>, u32)> =
            entries.iter().filter(|(_, c)| *c > 0).cloned().collect();
        seq_counts.sort_by(|(a, _), (b, _)| a.cmp(b));
        SampleEvidence {
            seq_counts,
            qc: SsrQc::default(),
        }
    }

    /// Params carrying a real per-base error (for the seed's composition fallback path).
    fn params_with_eps(p_decay: f64, eps: f64) -> ParamSet {
        let mut params = params_with_decay(p_decay);
        params.error_per_sample_group = vec![PerBaseError(eps)];
        params
    }

    /// The two same-length candidates `[pure6, interrupted6]` (ref = pure at index 0).
    fn same_length_candidates() -> CandidateSet {
        CandidateSet {
            alleles: vec![pure6(), interrupted6()],
            ref_idx: 0,
            admit: Admission::Pass,
        }
    }

    #[test]
    fn same_length_homozygote_seeds_its_own_composition() {
        // A sample homozygous for the *interrupted* allele must seed the interrupted candidate
        // (index 1), not the first-at-length (pure, index 0) the old length-only seed loaded.
        // G₀ starts both same-length candidates equal (both at the modal length), so any
        // asymmetry in π⁰ comes purely from the composition-aware seed (spec §5.3).
        let sample = seq_sample(&[(interrupted6(), 50)]);
        let locus = locus_with(6, vec![sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let s = seed_locus(
            &locus,
            &rungs,
            &same_length_candidates(),
            &params_with_decay(0.5),
            2,
            3,
        );
        assert!(
            s.pi0[1] > s.pi0[0],
            "the interrupted allele's own composition must dominate the seed"
        );
        assert!(s.pi0[1] > 0.5, "both copies seed the interrupted candidate");
    }

    #[test]
    fn same_length_het_seeds_a_hom_for_the_majority_documented_limitation() {
        // A same-length het (majority interrupted, minority pure) shows ONE length peak with
        // ONE representative, so the seed cannot express the het: both copies load the
        // majority (interrupted), and the minority (pure) starts only at its G₀ floor. This
        // is the documented Phase-1 limitation (spec §5.3) — a correct het seed would be
        // ~0.5/0.5; recovering it is the EM's job (P1.4), asserted in P1.5.
        let sample = seq_sample(&[(pure6(), 10), (interrupted6(), 40)]);
        let locus = locus_with(6, vec![sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let s = seed_locus(
            &locus,
            &rungs,
            &same_length_candidates(),
            &params_with_decay(0.5),
            2,
            3,
        );
        // Both copies on the majority → π⁰ ≈ [0.25, 0.75], NOT the ~[0.5, 0.5] of a het seed.
        assert!(
            s.pi0[1] > 0.7,
            "both allele copies seed the majority (a hom seed), not a het"
        );
        assert!(
            s.pi0[0] < 0.3,
            "the minority same-length allele gets only its G₀ floor"
        );
    }

    #[test]
    fn seed_falls_back_to_the_composition_closest_candidate() {
        // A sample whose plurality sequence is a private double-substitution variant that is
        // NOT an admitted candidate. The seed's exact match fails, so it falls back to the
        // closest same-length candidate by `align_subst`: the variant differs from
        // interrupted6 by one base and from pure6 by two, so the interrupted candidate (index
        // 1) must win — exercising the fallback + `sample_eps` path.
        let variant6: Box<[u8]> = {
            let mut s = interrupted6().to_vec();
            s[9] = b'G'; // a second substitution → 1 base from interrupted6, 2 from pure6
            s.into_boxed_slice()
        };
        assert_eq!(variant6.len(), 12, "still a same-length (12 bp) variant");
        let sample = seq_sample(&[(variant6, 50)]);
        let locus = locus_with(6, vec![sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let s = seed_locus(
            &locus,
            &rungs,
            &same_length_candidates(),
            &params_with_eps(0.5, 0.01),
            2,
            3,
        );
        assert!(
            s.pi0[1] > s.pi0[0],
            "the composition-closest candidate (interrupted) must take the seed mass"
        );
    }

    #[test]
    #[should_panic(expected = "frozen sample group")]
    fn sample_eps_panics_on_a_sample_missing_its_group() {
        // Decision-E loud-fail on the seed's composition-fallback path: a ParamSet whose
        // `group_of_sample` is shorter than the present samples must panic, not fabricate
        // default-ε chemistry — mirroring `em::sample_chemistry`'s regression test (review Mi5).
        let sample = seq_sample(&[(interrupted6(), 50)]);
        let locus = locus_with(6, vec![sample]);
        let mut params = params_with_eps(0.5, 0.01);
        params.group_of_sample.clear();
        let _ = sample_eps(&locus, 0, &params);
    }

    #[test]
    fn candidate_for_sequence_returns_none_when_no_candidate_sits_at_the_length() {
        // No candidate at the target length → None (the fallback must not return a wrong-length
        // candidate). review §8.
        let candidates = same_length_candidates(); // both at length 6
        let mut scratch = HmmScratch::new();
        let idx =
            candidate_for_sequence(&candidates, ca_seq(9).as_ref(), 2, 9, || 0.01, &mut scratch);
        assert_eq!(idx, None);
    }
}
