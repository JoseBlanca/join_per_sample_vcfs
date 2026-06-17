//! Per-read analysis — the worker's core step (arch §2 revision, realign-everything).
//!
//! Composes the per-read pipeline: classify the read against the locus
//! ([`super::triage`]), and for a spanning read realign it against the candidate
//! rungs ([`super::candidate_generation`] → [`super::pair_hmm`]) to produce its
//! dense `Qᵣ` (one log-likelihood per candidate). The locus aggregator
//! (`locus_record`, pending the container schema) then sparsifies and tallies
//! these per-read results.
//!
//! **Off-ladder candidate generation is not yet wired here.** A genuinely
//! off-ladder read is currently scored against the rungs only (it votes for the
//! nearest rung). Off-ladder is empty at the vast majority of loci (§11), and
//! wiring it needs the read's *observed tract* isolated from the extracted region
//! — well-defined only for both-flanks-aligned reads (§5.8) — a focused follow-up.

use crate::bam::alignment_input::MappedRead;
use crate::ssr_mark1::types::{Allele, Locus};

use super::candidate_generation::{CandidateAllele, build_rungs};
use super::pair_hmm::{HmmModel, PairHmmScratch, score_candidates};
use super::triage::{TriageResult, triage_read};

/// The per-read outcome the locus aggregator consumes. `Spanning` carries the
/// dense `Qᵣ` — one `(allele, forward log-likelihood)` per candidate rung;
/// `Flanking` / `InRepeat` carry no length evidence and are tallied for QC only
/// (arch §3.3).
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum ReadOutcome {
    Spanning(Vec<(Allele, f64)>),
    Flanking,
    InRepeat,
}

/// Analyze one read against one locus. Reuses the caller's scratch buffers — the
/// candidate set and the pair-HMM rows — across reads, so there is no per-read
/// allocation beyond the returned `Qᵣ` vec. `window` is the candidate half-width
/// (`STUTTER_WINDOW_UNITS`).
pub(crate) fn analyze_read(
    read: &MappedRead,
    locus: &Locus,
    window: u16,
    model: &HmmModel,
    candidates: &mut Vec<CandidateAllele>,
    hmm_scratch: &mut PairHmmScratch,
) -> ReadOutcome {
    match triage_read(read, locus) {
        TriageResult::Spanning(spanning) => {
            candidates.clear();
            build_rungs(locus, spanning.observed_count, window, candidates);
            let region = spanning.region;
            let region_seq = &read.seq[region.clone()];
            let region_qual = &read.qual[region];
            let scores = score_candidates(region_seq, region_qual, candidates, hmm_scratch, model);
            ReadOutcome::Spanning(scores)
        }
        TriageResult::Flanking => ReadOutcome::Flanking,
        TriageResult::InRepeat => ReadOutcome::InRepeat,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::walker::CigarOp;
    use crate::ssr_mark1::pileup::candidate_generation::STUTTER_WINDOW_UNITS;
    use crate::ssr_mark1::types::Motif;

    /// CA locus with 6 bp flanks: GGGGGG | CACACA | TTTTTT, tract ref [16,22).
    fn locus6() -> Locus {
        Locus::new(
            "chr1".into(),
            16,
            22,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGGGGCACACATTTTTT").into(),
            10,
        )
        .unwrap()
    }

    fn mapped_read(pos: u64, cigar: Vec<CigarOp>, seq: &[u8]) -> MappedRead {
        MappedRead {
            qname: b"r".to_vec(),
            flag: 0,
            ref_id: 0,
            pos,
            mapq: 60,
            cigar,
            seq: seq.to_vec(),
            qual: vec![40; seq.len()],
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            source_file_index: 0,
        }
    }

    #[test]
    fn spanning_read_scores_highest_against_its_own_rung() {
        let locus = locus6();
        let model = HmmModel::new();
        let mut cands = Vec::new();
        let mut hmm = PairHmmScratch::new();
        // A clean 3-unit read == rung 3's candidate (GGGGGG + CACACA + TTTTTT).
        let read = mapped_read(11, vec![CigarOp::Match(18)], b"GGGGGGCACACATTTTTT");
        match analyze_read(
            &read,
            &locus,
            STUTTER_WINDOW_UNITS,
            &model,
            &mut cands,
            &mut hmm,
        ) {
            ReadOutcome::Spanning(scores) => {
                // observed_count 3 ± 3 → rungs 0..=6.
                assert_eq!(scores.len(), 7);
                let best = scores
                    .iter()
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                    .unwrap();
                assert_eq!(best.0, Allele::OnLadder { units: 3 });
            }
            other => panic!("expected Spanning, got {other:?}"),
        }
    }

    #[test]
    fn longer_allele_read_scores_highest_against_the_recovered_rung() {
        // A soft-clipped 6-unit allele (mapper aligned 3 + clip): triage recovers
        // observed_count 6, and the read realigns best against rung 6.
        let locus = locus6();
        let model = HmmModel::new();
        let mut cands = Vec::new();
        let mut hmm = PairHmmScratch::new();
        let read = mapped_read(
            11,
            vec![CigarOp::Match(12), CigarOp::SoftClip(12)],
            b"GGGGGGCACACACACACATTTTTT",
        );
        match analyze_read(
            &read,
            &locus,
            STUTTER_WINDOW_UNITS,
            &model,
            &mut cands,
            &mut hmm,
        ) {
            ReadOutcome::Spanning(scores) => {
                let best = scores
                    .iter()
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                    .unwrap();
                assert_eq!(best.0, Allele::OnLadder { units: 6 });
            }
            other => panic!("expected Spanning, got {other:?}"),
        }
    }

    #[test]
    fn non_spanning_reads_carry_no_scores() {
        let locus = locus6();
        let model = HmmModel::new();
        let mut cands = Vec::new();
        let mut hmm = PairHmmScratch::new();

        let flanking = mapped_read(6, vec![CigarOp::Match(12)], b"GGGGGCACACAC");
        assert_eq!(
            analyze_read(
                &flanking,
                &locus,
                STUTTER_WINDOW_UNITS,
                &model,
                &mut cands,
                &mut hmm
            ),
            ReadOutcome::Flanking
        );

        let in_repeat = mapped_read(18, vec![CigarOp::Match(4)], b"CACA");
        assert_eq!(
            analyze_read(
                &in_repeat,
                &locus,
                STUTTER_WINDOW_UNITS,
                &model,
                &mut cands,
                &mut hmm
            ),
            ReadOutcome::InRepeat
        );
    }

    #[test]
    fn scratch_buffers_are_reused_across_reads() {
        // The candidate buffer is cleared per read, so scoring a second read after
        // a first gives the first read's own result, not leftovers.
        let locus = locus6();
        let model = HmmModel::new();
        let mut cands = Vec::new();
        let mut hmm = PairHmmScratch::new();
        let r1 = mapped_read(11, vec![CigarOp::Match(18)], b"GGGGGGCACACATTTTTT");
        let _ = analyze_read(
            &r1,
            &locus,
            STUTTER_WINDOW_UNITS,
            &model,
            &mut cands,
            &mut hmm,
        );
        let r2 = mapped_read(11, vec![CigarOp::Match(18)], b"GGGGGGCACACATTTTTT");
        match analyze_read(
            &r2,
            &locus,
            STUTTER_WINDOW_UNITS,
            &model,
            &mut cands,
            &mut hmm,
        ) {
            ReadOutcome::Spanning(scores) => assert_eq!(scores.len(), 7),
            other => panic!("expected Spanning, got {other:?}"),
        }
    }
}
