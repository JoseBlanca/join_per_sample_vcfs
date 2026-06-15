//! Candidate generation — the translation layer between read length-estimates
//! and typed alleles (arch §6).
//!
//! Produces the candidate haplotype sequences the pair-HMM scores (arch §5):
//! the `observed_count ± W` on-ladder rungs and — when a read earns one — a
//! single read-derived off-ladder candidate. The forward never constructs a
//! sequence; it only scores what this module hands it (arch §5.7).
//!
//! Each candidate is built from the **catalog** locus (its embedded reference
//! bytes + motif), so a given rung is the byte-identical sequence in every
//! sample — which is what lets Stage 2 recognise "the L-unit allele" as the same
//! allele across the cohort by its count alone (arch §6 Job 1).
//!
//! **Built so far: Job 1 (the on-ladder rungs).** Job 2 — the read-derived
//! off-ladder candidate and its canonical normalization (`normalize_offladder`,
//! the adapter onto [`crate::norm_seqs`]) — is deferred until its input contract
//! is pinned: what the triage stage hands in as the observed tract (raw bytes
//! vs. an alignment carrying the indel position) and what `NormalizedSeq`
//! canonicalizes to operationally. Driving the indel-norm kernel needs the
//! variant *bounds*, which raw tract bytes alone do not give; guessing risks the
//! cross-sample-identity invariant the normalization exists to protect.

use crate::ssr::types::{Allele, Locus};

/// One candidate allele the read is scored against: its full local DNA sequence
/// (what the forward aligns the read to) + which [`Allele`] that sequence
/// encodes (the key the `Qᵣ` result is recorded under). `candidate_seq` is
/// `left_flank + tract + right_flank`, where the tract is the motif tiled
/// `units` times for an on-ladder rung.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CandidateAllele {
    pub(crate) candidate_seq: Vec<u8>,
    pub(crate) allele: Allele,
}

/// Half-width of the candidate-`L` window: a slow-path read's forward is
/// evaluated over `observed_count ± STUTTER_WINDOW_UNITS` rungs (arch §5).
///
/// Named for stutter even though Stage 1 is stutter-free (arch §10 naming note):
/// the window must be wide enough that Stage 2 has support points to convolve
/// the stutter kernel against. The value is a **calibration** parameter (arch
/// §14), a placeholder until tuned; callers pass it explicitly to `build_rungs`.
pub(crate) const STUTTER_WINDOW_UNITS: u16 = 3;

/// Build the on-ladder rungs for a read's window — one [`CandidateAllele`] per
/// repeat count `L` in `[observed_count − w, observed_count + w]`, each a clean
/// tiling `left_flank + (motif × L) + right_flank` read off the catalog locus
/// (no FASTA). Appends to `out` (a reused per-worker scratch); does not clear it.
///
/// The low end is clamped at `0` (a repeat count cannot be negative), so a read
/// near the bottom of the ladder yields fewer than `2w + 1` rungs rather than
/// wrapping. The reference rung (`L = ref units`) reproduces the locus's
/// embedded `ref_bytes` exactly for a perfect locus — the property that makes a
/// rung sample-independent.
pub(crate) fn build_rungs(
    locus: &Locus,
    observed_count: u16,
    w: u16,
    out: &mut Vec<CandidateAllele>,
) {
    let lo = observed_count.saturating_sub(w);
    let hi = observed_count.saturating_add(w);
    let left = locus.left_flank();
    let right = locus.right_flank();

    for units in lo..=hi {
        let allele = Allele::OnLadder { units };
        let tract = allele.to_sequence(locus);
        let mut candidate_seq = Vec::with_capacity(left.len() + tract.len() + right.len());
        candidate_seq.extend_from_slice(left);
        candidate_seq.extend_from_slice(&tract);
        candidate_seq.extend_from_slice(right);
        out.push(CandidateAllele {
            candidate_seq,
            allele,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::types::Motif;

    /// A perfect dinucleotide locus: 3bp left flank + (CA)×3 tract + 3bp right
    /// flank, tract at genomic [13, 19); `ref_bytes` = `GGGCACACATTT`.
    fn sample_locus() -> Locus {
        Locus::new(
            "chr1".into(),
            13,
            19,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGCACACATTT").into(),
            10,
        )
        .unwrap()
    }

    #[test]
    fn builds_a_full_window_of_rungs() {
        let locus = sample_locus();
        let mut out = Vec::new();
        build_rungs(&locus, 5, 2, &mut out); // L ∈ 3..=7 → 5 rungs
        assert_eq!(out.len(), 5);
        let counts: Vec<u16> = out
            .iter()
            .map(|c| match c.allele {
                Allele::OnLadder { units } => units,
                _ => panic!("rungs are on-ladder"),
            })
            .collect();
        assert_eq!(counts, vec![3, 4, 5, 6, 7]);
    }

    #[test]
    fn each_rung_is_left_flank_plus_tiling_plus_right_flank() {
        let locus = sample_locus();
        let mut out = Vec::new();
        build_rungs(&locus, 2, 0, &mut out); // single rung L=2
        assert_eq!(out.len(), 1);
        // GGG + CACA + TTT
        assert_eq!(out[0].candidate_seq, b"GGGCACATTT");
        assert_eq!(out[0].allele, Allele::OnLadder { units: 2 });
    }

    #[test]
    fn reference_rung_reproduces_the_locus_ref_bytes() {
        // A perfect locus's reference rung (L = ref units = 3) rebuilds the
        // embedded ref_bytes exactly — the sample-independence property.
        let locus = sample_locus();
        let mut out = Vec::new();
        build_rungs(&locus, 3, 0, &mut out);
        assert_eq!(out[0].candidate_seq, locus.ref_bytes());
    }

    #[test]
    fn low_end_clamps_at_zero_rather_than_wrapping() {
        let locus = sample_locus();
        let mut out = Vec::new();
        build_rungs(&locus, 1, 3, &mut out); // would be -2..=4; clamps to 0..=4
        let counts: Vec<u16> = out
            .iter()
            .map(|c| match c.allele {
                Allele::OnLadder { units } => units,
                _ => unreachable!(),
            })
            .collect();
        assert_eq!(counts, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn zero_unit_rung_is_just_the_flanks() {
        let locus = sample_locus();
        let mut out = Vec::new();
        build_rungs(&locus, 0, 0, &mut out);
        // L=0 → empty tract → left + right flank only.
        assert_eq!(out[0].candidate_seq, b"GGGTTT");
        assert_eq!(out[0].allele, Allele::OnLadder { units: 0 });
    }

    #[test]
    fn appends_without_clearing_the_scratch() {
        let locus = sample_locus();
        let mut out = Vec::new();
        build_rungs(&locus, 2, 0, &mut out);
        build_rungs(&locus, 5, 0, &mut out);
        assert_eq!(out.len(), 2); // reused buffer accumulates
    }
}
