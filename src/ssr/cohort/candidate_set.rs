//! Candidate assembly (arch `ssr_call_genotyping.md` §2, spec §5,
//! implementation plan §3.3, C1).
//!
//! [`assemble_candidates`] turns the cohort rung ladder into the per-locus
//! [`CandidateSet`]: a **locus-admission** check (the adjacent observed-length
//! spacing must be motif-dominated), per-sample **nomination** (top-ploidy clear
//! peaks, or a `±1` neighbour rescue when a sample under-resolves), the **union**
//! across samples, the **reference allele seeded unconditionally**, and a cap on the
//! allele count. Recall is the goal — a kept-but-spurious candidate is the EM's job
//! to drive to `π ≈ 0`, but a dropped one is gone (spec §5).

use std::collections::BTreeMap;

use crate::ssr::cohort::rung_ladder::{Rungs, sample_clear_peaks};
use crate::ssr::cohort::types::CohortLocus;

/// Why a locus was (or was not) admitted to genotyping — the per-site FILTER
/// reason (arch §6).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Admission {
    /// The locus is admitted.
    Pass,
    /// Adjacent-rung spacing did not match the motif length — not a periodic
    /// ladder (`notPeriodic`).
    NotPeriodic,
    /// The candidate count exceeded the per-locus cap (`tooManyAlleles`).
    TooManyAlleles,
    /// Too little depth to resolve candidates (`lowDepth`).
    LowDepth,
}

/// The candidate alleles for one locus, plus its admission verdict.
///
/// Each entry of [`Self::alleles`] is an independent candidate **tract sequence**
/// (impure peaks kept first-class); [`Self::ref_idx`] indexes the reference allele,
/// which is seeded unconditionally.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CandidateSet {
    /// Candidate tract sequences (each an independent allele).
    pub(crate) alleles: Vec<Box<[u8]>>,
    /// Index of the reference allele within [`Self::alleles`].
    pub(crate) ref_idx: usize,
    /// The site's admission verdict → FILTER.
    pub(crate) admit: Admission,
}

/// Thresholds for candidate assembly (pinned in F2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CandidateCfg {
    /// Prominence floor for a clear local maximum (shares the rung definition).
    pub(crate) prominence: u32,
    /// Cohort total depth below which the locus is `LowDepth` (no-call).
    pub(crate) min_cohort_depth: u32,
    /// Maximum number of candidate alleles before the locus is `TooManyAlleles`.
    pub(crate) max_candidate_alleles: usize,
    /// Minimum number of distinct observed lengths needed to *run* the periodicity
    /// test; below this the locus is admitted by default (too little to judge).
    pub(crate) min_lengths_for_admission: usize,
}

impl CandidateCfg {
    /// Working defaults (recalibrated in F2).
    pub(crate) fn dev_default() -> Self {
        Self {
            prominence: 3,
            min_cohort_depth: 10,
            max_candidate_alleles: 24,
            min_lengths_for_admission: 3,
        }
    }
}

/// The total reads observed across the cohort at a locus.
fn cohort_depth(locus: &CohortLocus) -> u32 {
    locus
        .samples
        .iter()
        .flat_map(|ev| ev.seq_counts.iter().map(|(_, c)| *c))
        .sum()
}

/// Whether the locus behaves like an SSR: the **mode** of the adjacent
/// observed-length differences (in bytes) equals the motif period. Loci with too
/// few distinct lengths to judge are admitted by default (spec §5 level 3).
fn is_periodic(locus: &CohortLocus, cfg: &CandidateCfg) -> bool {
    let period = locus.motif.period();
    let mut lengths: Vec<usize> = locus
        .samples
        .iter()
        .flat_map(|ev| ev.seq_counts.iter().map(|(seq, _)| seq.len()))
        .collect();
    lengths.sort_unstable();
    lengths.dedup();
    if lengths.len() < cfg.min_lengths_for_admission {
        return true;
    }
    let mut diff_counts: BTreeMap<usize, u32> = BTreeMap::new();
    for window in lengths.windows(2) {
        *diff_counts.entry(window[1] - window[0]).or_insert(0) += 1;
    }
    // The dominant adjacent spacing must be the motif period.
    diff_counts
        .iter()
        .max_by(|(da, ca), (db, cb)| ca.cmp(cb).then_with(|| db.cmp(da)))
        .map(|(diff, _)| *diff == period)
        .unwrap_or(true)
}

/// The cohort's most-supported sequence at rung `length`, if that rung is occupied.
fn cohort_representative(rungs: &Rungs, length: u16) -> Option<Box<[u8]>> {
    rungs
        .seqs_at(length)?
        .iter()
        .max_by(|(a_seq, a_c), (b_seq, b_c)| a_c.cmp(b_c).then_with(|| b_seq.cmp(a_seq)))
        .map(|(seq, _)| seq.clone())
}

/// Assemble the candidate set for one locus from its rung ladder (spec §5).
pub(crate) fn assemble_candidates(
    locus: &CohortLocus,
    rungs: &Rungs,
    ploidy: u8,
    cfg: &CandidateCfg,
) -> CandidateSet {
    // The reference allele is always candidate 0 (REF must exist for the VCF).
    let mut alleles: Vec<Box<[u8]>> = vec![locus.ref_tract.clone()];
    let ref_idx = 0;

    // Admission: depth first (cheapest no-call), then periodicity.
    if cohort_depth(locus) < cfg.min_cohort_depth {
        return CandidateSet {
            alleles,
            ref_idx,
            admit: Admission::LowDepth,
        };
    }
    if !is_periodic(locus, cfg) {
        return CandidateSet {
            alleles,
            ref_idx,
            admit: Admission::NotPeriodic,
        };
    }

    let period = rungs.period();
    let occupied = |length: u16| rungs.cohort_support(length) > 0;

    let push_unique = |alleles: &mut Vec<Box<[u8]>>, seq: Box<[u8]>| {
        if !alleles.contains(&seq) {
            alleles.push(seq);
        }
    };

    for evidence in &locus.samples {
        let mut peaks = sample_clear_peaks(evidence, period, cfg.prominence);
        peaks.sort_by(|a, b| {
            b.support
                .cmp(&a.support)
                .then_with(|| a.repeat_len.cmp(&b.repeat_len))
        });

        // Lengths this sample nominates: its top-ploidy peaks, plus the ±1 observed
        // neighbours when it under-resolves (a possible merged/hidden allele).
        let mut nominated: Vec<u16> = peaks
            .iter()
            .take(ploidy as usize)
            .map(|p| p.repeat_len)
            .collect();
        if peaks.len() < ploidy as usize {
            for peak in &peaks {
                for neighbour in [
                    peak.repeat_len.checked_sub(1),
                    peak.repeat_len.checked_add(1),
                ]
                .into_iter()
                .flatten()
                {
                    if occupied(neighbour) {
                        nominated.push(neighbour);
                    }
                }
            }
        }

        for length in nominated {
            if let Some(seq) = cohort_representative(rungs, length) {
                push_unique(&mut alleles, seq);
            }
        }
    }

    let admit = if alleles.len() > cfg.max_candidate_alleles {
        Admission::TooManyAlleles
    } else {
        Admission::Pass
    };
    CandidateSet {
        alleles,
        ref_idx,
        admit,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::rung_ladder::{RungCfg, build_rungs};
    use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};
    use crate::ssr::types::Motif;

    fn ca_evidence(bins: &[(u16, u32)]) -> SampleEvidence {
        let mut seq_counts: Vec<(Box<[u8]>, u32)> = bins
            .iter()
            .map(|&(units, count)| {
                let mut seq = Vec::new();
                for _ in 0..units {
                    seq.extend_from_slice(b"CA");
                }
                (seq.into_boxed_slice(), count)
            })
            .collect();
        seq_counts.sort_by(|(a, _), (b, _)| a.cmp(b));
        SampleEvidence {
            seq_counts,
            qc: SsrQc::default(),
        }
    }

    /// A CA-locus cohort with a reference tract of `ref_units` units.
    fn ca_cohort(ref_units: u16, samples: Vec<SampleEvidence>) -> CohortLocus {
        let ref_tract: Vec<u8> = std::iter::repeat_n(*b"CA", ref_units as usize)
            .flatten()
            .collect();
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 20,
                end: 20 + ref_tract.len() as u32,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGG".as_slice()), // ref_frame unused by C1
            ref_tract.into_boxed_slice(),
        );
        for (idx, evidence) in samples.into_iter().enumerate() {
            locus.push(idx as u32, evidence);
        }
        locus
    }

    fn ca_seq(units: u16) -> Box<[u8]> {
        std::iter::repeat_n(*b"CA", units as usize)
            .flatten()
            .collect()
    }

    fn assemble(locus: &CohortLocus) -> CandidateSet {
        let rungs = build_rungs(locus, &RungCfg::dev_default());
        assemble_candidates(locus, &rungs, 2, &CandidateCfg::dev_default())
    }

    #[test]
    fn reference_allele_is_seeded_at_index_zero() {
        let sample = ca_evidence(&[(5, 5), (6, 50), (7, 5)]);
        let locus = ca_cohort(6, vec![sample.clone(), sample]);
        let cs = assemble(&locus);
        assert_eq!(cs.admit, Admission::Pass);
        assert_eq!(&*cs.alleles[cs.ref_idx], &*ca_seq(6));
    }

    #[test]
    fn separated_het_nominates_both_alleles_without_rescue() {
        let sample = ca_evidence(&[(3, 1), (4, 30), (5, 4), (8, 4), (9, 28), (10, 1)]);
        let locus = ca_cohort(6, vec![sample.clone(), sample]);
        let cs = assemble(&locus);
        assert_eq!(cs.admit, Admission::Pass);
        assert!(cs.alleles.iter().any(|a| **a == *ca_seq(4)));
        assert!(cs.alleles.iter().any(|a| **a == *ca_seq(9)));
    }

    #[test]
    fn under_resolved_sample_rescues_the_adjacent_neighbour() {
        // One clear peak at 6 (ploidy 2 → under-resolved), with an occupied rung at
        // 7: the ±1 rescue must add the length-7 sequence as a candidate.
        let sample = ca_evidence(&[(5, 4), (6, 50), (7, 20)]);
        let locus = ca_cohort(6, vec![sample.clone(), sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        // Sanity: 6 is a clear peak, 7 is not (it sits on 6's shoulder).
        let peaks = sample_clear_peaks(&locus.samples[0], 2, 3);
        assert_eq!(
            peaks.iter().map(|p| p.repeat_len).collect::<Vec<_>>(),
            vec![6]
        );

        let cs = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
        assert!(
            cs.alleles.iter().any(|a| **a == *ca_seq(7)),
            "the ±1 rescue should add the occupied neighbour at length 7"
        );
    }

    #[test]
    fn non_periodic_locus_is_not_admitted() {
        // Observed byte lengths 10, 13, 17 → diffs 3, 4 (mode ≠ period 2).
        let sample = SampleEvidence {
            seq_counts: vec![
                (Box::from(&b"AAAAAAAAAA"[..]), 20),        // len 10
                (Box::from(&b"AAAAAAAAAAAAA"[..]), 20),     // len 13
                (Box::from(&b"AAAAAAAAAAAAAAAAA"[..]), 20), // len 17
            ],
            qc: SsrQc::default(),
        };
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 20,
                end: 30,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGG".as_slice()),
            Box::from(b"CACACA".as_slice()),
        );
        locus.push(0, sample.clone());
        locus.push(1, sample);
        let cs = assemble(&locus);
        assert_eq!(cs.admit, Admission::NotPeriodic);
    }

    #[test]
    fn thin_locus_is_low_depth() {
        let sample = ca_evidence(&[(6, 3)]);
        let locus = ca_cohort(6, vec![sample]);
        assert_eq!(assemble(&locus).admit, Admission::LowDepth);
    }

    #[test]
    fn too_many_alleles_is_flagged() {
        // Many distinct well-separated alleles across samples exceed the cap.
        let cfg = CandidateCfg {
            max_candidate_alleles: 3,
            ..CandidateCfg::dev_default()
        };
        let samples: Vec<SampleEvidence> = (4..12)
            .map(|u| ca_evidence(&[(u - 2, 1), (u, 40), (u + 2, 1)]))
            .collect();
        let locus = ca_cohort(6, samples);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());
        let cs = assemble_candidates(&locus, &rungs, 2, &cfg);
        assert_eq!(cs.admit, Admission::TooManyAlleles);
    }

    #[test]
    fn candidate_set_ref_idx_points_into_alleles() {
        let cs = CandidateSet {
            alleles: vec![
                Box::from(b"ATAT".as_slice()),
                Box::from(b"ATATAT".as_slice()),
            ],
            ref_idx: 0,
            admit: Admission::Pass,
        };
        assert_eq!(&*cs.alleles[cs.ref_idx], b"ATAT");
        assert_eq!(cs.admit, Admission::Pass);
    }

    #[test]
    fn admission_variants_are_distinct() {
        assert_ne!(Admission::Pass, Admission::NotPeriodic);
        assert_ne!(Admission::TooManyAlleles, Admission::LowDepth);
    }
}
