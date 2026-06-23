//! The `G₀` geometric pseudocount prior on candidate allele frequencies `π` (arch
//! `ssr_call_genotyping.md` §5, spec §4.3/§5.5, implementation plan §3, C3).
//!
//! Every candidate gets a small prior "count" that decays geometrically with its
//! **unit offset** `Δ` from the per-locus cohort *modal* allele — `p^|Δ|`, floored
//! at a tiny positive constant so a far candidate (a masked het the seed missed)
//! keeps `π > 0` and stays recoverable rather than falling into the `π = 0`
//! absorbing trap (verify-fix #4). It re-enters **every** EM M-step as a small-N
//! regularizer + false-positive control. The decay `p` is fit per loci group in the
//! pre-pass (D); here it is supplied.

use crate::ssr::cohort::candidate_set::CandidateSet;
use crate::ssr::cohort::param_estimation::G0PseudocountDecay;
use crate::ssr::cohort::rung_ladder::Rungs;

/// The pseudocount floor: any representable `> 0` so `p^|Δ|` cannot underflow to
/// exactly `0.0` for a far candidate over a long tract (verify-fix #4).
const G0_FLOOR: f64 = 1e-12;

/// The `G₀` pseudocount vector, parallel to `candidates.alleles`: `max(p^|Δ|, FLOOR)`
/// where `Δ` is each candidate's unit offset from the cohort modal length.
///
/// With no observed data (the only candidate is the seeded reference) the modal
/// length falls back to the reference's own length, so `Δ = 0` there.
pub(crate) fn g0_pseudocounts(
    candidates: &CandidateSet,
    rungs: &Rungs,
    decay: &G0PseudocountDecay,
) -> Vec<f64> {
    let period = rungs.period();
    let modal = rungs
        .modal_length()
        .unwrap_or_else(|| (candidates.alleles[candidates.ref_idx].len() / period) as u16);
    candidates
        .alleles
        .iter()
        .map(|seq| {
            let units = (seq.len() / period) as i32;
            let delta = (units - modal as i32).abs();
            decay.p.powi(delta).max(G0_FLOOR)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::candidate_set::Admission;
    use crate::ssr::cohort::rung_ladder::{RungCfg, build_rungs};
    use crate::ssr::cohort::types::{CohortLocus, LocusId, SampleEvidence, SsrQc};
    use crate::ssr::types::Motif;

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

    fn locus_with(samples: Vec<SampleEvidence>) -> CohortLocus {
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 20,
                end: 32,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGG".as_slice()),
            ca_seq(6),
        );
        for (idx, ev) in samples.into_iter().enumerate() {
            locus.push(idx as u32, ev);
        }
        locus
    }

    #[test]
    fn pseudocounts_decay_with_distance_from_the_mode() {
        // Cohort mode at 6 units; candidates at 4, 6, 9.
        let sample = ca_evidence(&[(5, 5), (6, 50), (7, 5)]);
        let locus = locus_with(vec![sample.clone(), sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let candidates = CandidateSet {
            alleles: vec![ca_seq(6), ca_seq(4), ca_seq(9)],
            ref_idx: 0,
            admit: Admission::Pass,
        };
        let g0 = g0_pseudocounts(&candidates, &rungs, &G0PseudocountDecay { p: 0.5 });
        assert!((g0[0] - 1.0).abs() < 1e-12); // Δ=0 → p^0 = 1
        assert!((g0[1] - 0.5f64.powi(2)).abs() < 1e-12); // Δ=2
        assert!((g0[2] - 0.5f64.powi(3)).abs() < 1e-12); // Δ=3
        assert!(g0[0] > g0[1] && g0[1] > g0[2]);
    }

    #[test]
    fn a_far_candidate_stays_positive_via_the_floor() {
        let sample = ca_evidence(&[(6, 60)]);
        let locus = locus_with(vec![sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        // A very distant candidate with a steep decay would underflow without the floor.
        let candidates = CandidateSet {
            alleles: vec![ca_seq(6), ca_seq(200)],
            ref_idx: 0,
            admit: Admission::Pass,
        };
        let g0 = g0_pseudocounts(&candidates, &rungs, &G0PseudocountDecay { p: 0.1 });
        assert!(g0[1] > 0.0, "the floor must keep a far candidate positive");
        assert_eq!(g0[1], G0_FLOOR);
    }
}
