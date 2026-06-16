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
//! **The off-ladder canonical form (Job 2) is the full tract, stored verbatim
//! (contract A + B1).** Triage hands in the raw observed-tract bytes; the stored
//! [`NormalizedSeq`] is the whole tract sequence (so [`Allele::to_sequence`]
//! means the same thing — the tract — for both allele kinds). No left-alignment
//! is applied, and that is **correct, not a shortcut**: a literal alt *sequence*
//! is invariant under left-alignment (left-alignment canonicalizes an indel's
//! *position in a `(ref, alt)` description*, never the alt bytes), and Stage 0
//! guarantees clean flanks (arch §5.9) so the tract is deterministically
//! delimited — hence two reads of the same molecule yield the same tract bytes,
//! already canonical. The shared [`crate::norm_seqs`] kernel would only be needed
//! for a *trimmed-variant* representation (the rejected delta form, shared-types
//! §4); revisit if we ever adopt that or relax the clean-flank guarantee.

use crate::ssr::pileup::count_repeats::pure_tiling_units;
use crate::ssr::types::{Allele, Locus, NormalizedSeq};

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

/// Canonicalize an observed off-ladder tract into its [`NormalizedSeq`] key
/// (Job 2). The canonical form **is** the tract bytes verbatim — see the module
/// doc: with the full-tract representation (B1) and Stage-0 clean flanks, the
/// sequence is already left-alignment-invariant, so there is nothing to shift.
/// Kept as a named function so the canonicalization contract lives at one place
/// (and is the single change point if a trimmed-variant form is ever adopted).
pub(crate) fn normalize_offladder(observed_tract: &[u8]) -> NormalizedSeq {
    NormalizedSeq::new(observed_tract.to_vec())
}

/// Build the single read-derived off-ladder candidate for an observed tract, or
/// `None` if the tract is a clean tiling (then it is an on-ladder rung, not
/// off-ladder — the degenerate case, arch §5.8). The candidate is `left_flank +
/// tract + right_flank` with the [`Allele::OffLadder`] key carrying the
/// canonical tract.
///
/// Whether a read *earns* an off-ladder candidate at all is gated upstream on
/// its slow-path reason (only impure / interior-indel reads, arch §5.8); this
/// builds the candidate once that gate has fired, so it takes just the locus and
/// the observed tract.
pub(crate) fn build_offladder(locus: &Locus, observed_tract: &[u8]) -> Option<CandidateAllele> {
    // Degenerate case: a pure tiling is a rung, already covered by `build_rungs`.
    if pure_tiling_units(observed_tract, &locus.motif()).is_some() {
        return None;
    }

    let normalized = normalize_offladder(observed_tract);
    let left = locus.left_flank();
    let right = locus.right_flank();
    let mut candidate_seq =
        Vec::with_capacity(left.len() + normalized.as_bytes().len() + right.len());
    candidate_seq.extend_from_slice(left);
    candidate_seq.extend_from_slice(normalized.as_bytes());
    candidate_seq.extend_from_slice(right);

    Some(CandidateAllele {
        candidate_seq,
        allele: Allele::OffLadder(normalized),
    })
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

    // --- Job 2: off-ladder candidates ----------------------------------------

    #[test]
    fn normalize_offladder_is_verbatim() {
        // B1: the canonical form is the tract bytes themselves.
        let seq = normalize_offladder(b"CACAACA");
        assert_eq!(seq.as_bytes(), b"CACAACA");
    }

    #[test]
    fn build_offladder_wraps_an_impure_tract_with_the_flanks() {
        let locus = sample_locus(); // motif CA, flanks GGG / TTT
        // "CACAACA" is not a clean (CA) tiling (the AA breaks it) → off-ladder.
        let cand = build_offladder(&locus, b"CACAACA").unwrap();
        assert_eq!(cand.candidate_seq, b"GGGCACAACATTT");
        assert_eq!(
            cand.allele,
            Allele::OffLadder(NormalizedSeq::new(b"CACAACA".to_vec()))
        );
    }

    #[test]
    fn build_offladder_rejects_a_pure_tiling_as_a_rung() {
        // A clean (CA) tiling is an on-ladder rung, not off-ladder (§5.8).
        let locus = sample_locus();
        assert!(build_offladder(&locus, b"CACACA").is_none());
        assert!(build_offladder(&locus, b"CACACACA").is_none());
        // The empty tract is the zero-unit rung, also a tiling.
        assert!(build_offladder(&locus, b"").is_none());
    }

    #[test]
    fn off_ladder_candidate_round_trips_through_to_sequence() {
        // The OffLadder allele's to_sequence is the canonical tract — the same
        // bytes build_offladder put between the flanks.
        let locus = sample_locus();
        let cand = build_offladder(&locus, b"CACAACA").unwrap();
        assert_eq!(cand.allele.to_sequence(&locus), b"CACAACA");
    }

    #[test]
    fn same_tract_two_reads_gives_an_identical_key() {
        // Cross-sample identity: the verbatim full-tract key is byte-equal for
        // two reads of the same molecule (there is only one spelling of a
        // literal sequence — the property that makes left-alignment unnecessary).
        let locus = sample_locus();
        let a = build_offladder(&locus, b"CACAACA").unwrap();
        let b = build_offladder(&locus, b"CACAACA").unwrap();
        assert_eq!(a.allele, b.allele);
    }
}
