//! Synthetic hot-path harness for the `ssr-pileup` perf benches/profiles.
//!
//! The realignment core of Stage 1 — triage → `build_rungs` → `score_candidates`
//! (the pair-HMM forward over `2·window + 1` rungs per spanning read) — is all
//! `pub(crate)`, and so are the domain types it runs on (`Locus`, `MappedRead`).
//! A `benches/` crate (or an `examples/` profile driver) can only touch the
//! public API, so this module is the *single* `#[doc(hidden)] pub` seam they use:
//! it builds a synthetic per-locus workload from primitive parameters and runs
//! [`analyze_read`] over it, exactly as [`super::driver::process_locus`] does
//! after the fetch/reservoir step.
//!
//! It deliberately covers only the **CPU-bound realignment**, not the indexed
//! fetch (that is I/O-bound and needs real alignment files — measured
//! separately). Everything here is gated behind `#[doc(hidden)]` and is not part
//! of the supported surface; it exists only so the perf review can measure the
//! hot path without making the whole SSR domain public.

use crate::bam::alignment_input::MappedRead;
use crate::pileup::walker::CigarOp;
use crate::ssr::types::{Locus, Motif};

use super::driver::LocusScratch;
use super::pair_hmm::HmmModel;
use super::read_analysis::{ReadOutcome, analyze_read};

/// A built synthetic locus plus its spanning reads, ready to feed the
/// realignment core. Opaque on purpose: the bench holds it as a token across the
/// untimed setup / timed measurement boundary.
#[doc(hidden)]
pub struct SyntheticLocusWorkload {
    locus: Locus,
    reads: Vec<MappedRead>,
    window: u16,
}

/// Build a synthetic SSR locus and `depth` spanning reads over it.
///
/// - `period` — motif length (1 = homopolymer, 2 = dinucleotide, …); the motif
///   is the first `period` bytes cycling through `ACGT…` (so it tiles cleanly).
/// - `units` — true repeat count of the embedded tract.
/// - `flank_bp` — clean flank each side of the tract (the catalog default is 50).
/// - `depth` — number of spanning reads built (the per-locus realignment count).
/// - `window` — pair-HMM candidate half-width; each read scores `2·window + 1`
///   rungs (the `ssr-pileup` default is `DEFAULT_WINDOW` = 6).
/// - `clip_units` — if non-zero, every read is given a longer true allele
///   (`units + clip_units` copies) whose extra tract + far flank sit in a
///   trailing soft-clip, exercising the soft-clip recovery + longer-haplotype DP.
///
/// All reads are Q30. Reads are clean tilings of the (possibly extended) tract,
/// so triage classifies them `Spanning` and the realignment runs the full DP for
/// every rung — the realign-everything common case.
#[doc(hidden)]
pub fn build_synthetic_workload(
    period: usize,
    units: u16,
    flank_bp: usize,
    depth: usize,
    window: u16,
    clip_units: u16,
) -> SyntheticLocusWorkload {
    assert!((1..=6).contains(&period), "period must be 1..=6");
    let motif_bytes: Vec<u8> = b"ACGTAC"[..period].to_vec();
    let motif = Motif::new(&motif_bytes).unwrap();

    let left: Vec<u8> = b"GC".iter().copied().cycle().take(flank_bp).collect();
    let right: Vec<u8> = b"AT".iter().copied().cycle().take(flank_bp).collect();
    let tract: Vec<u8> = motif_bytes
        .iter()
        .copied()
        .cycle()
        .take(period * units as usize)
        .collect();

    // Reference window (catalog ref_bytes): left + tract(units) + right, placed
    // at an arbitrary contig offset so coordinates are non-degenerate.
    let ref_bytes_start: u32 = 1_000;
    let mut ref_bytes = Vec::with_capacity(flank_bp * 2 + tract.len());
    ref_bytes.extend_from_slice(&left);
    ref_bytes.extend_from_slice(&tract);
    ref_bytes.extend_from_slice(&right);
    let start = ref_bytes_start + flank_bp as u32;
    let end = start + tract.len() as u32;

    let locus = Locus::new(
        "chr1".into(),
        start,
        end,
        motif,
        1.0,
        ref_bytes.clone().into_boxed_slice(),
        ref_bytes_start,
    )
    .unwrap();

    // The read's own (possibly longer) tract.
    let read_tract: Vec<u8> = motif_bytes
        .iter()
        .copied()
        .cycle()
        .take(period * (units + clip_units) as usize)
        .collect();
    let mut read_seq = Vec::with_capacity(flank_bp * 2 + read_tract.len());
    read_seq.extend_from_slice(&left);
    read_seq.extend_from_slice(&read_tract);
    read_seq.extend_from_slice(&right);

    // CIGAR: a clean read is one Match over its whole length, aligned at the
    // window start. A clipped read aligns left flank + `units` tract copies and
    // soft-clips the remainder (the extra tract + far flank), so triage must
    // recover the long allele from the clip.
    let cigar = if clip_units == 0 {
        vec![CigarOp::Match(read_seq.len() as u32)]
    } else {
        let aligned = flank_bp + period * units as usize;
        vec![
            CigarOp::Match(aligned as u32),
            CigarOp::SoftClip((read_seq.len() - aligned) as u32),
        ]
    };
    // pos is 1-based; the window opens at ref_bytes_start (0-based).
    let pos = u64::from(ref_bytes_start) + 1;

    let reads = (0..depth)
        .map(|i| MappedRead {
            qname: format!("r{i}").into_bytes(),
            flag: 0,
            ref_id: 0,
            pos,
            mapq: 60,
            cigar: cigar.clone(),
            seq: read_seq.clone(),
            qual: vec![30u8; read_seq.len()],
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            source_file_index: 0,
        })
        .collect();

    SyntheticLocusWorkload {
        locus,
        reads,
        window,
    }
}

/// Run the realignment core over every read in the workload, reusing one
/// [`LocusScratch`] and one [`HmmModel`] (exactly as `process_locus` does).
/// Returns a checksum of the produced log-likelihoods so the optimizer cannot
/// elide the work.
#[doc(hidden)]
pub fn analyze_workload(workload: &SyntheticLocusWorkload) -> f64 {
    let model = HmmModel::default();
    let mut scratch = LocusScratch::new();
    let mut sink = 0.0f64;
    for read in &workload.reads {
        let (cands, hmm) = scratch.parts_mut();
        let outcome = analyze_read(read, &workload.locus, workload.window, &model, cands, hmm);
        if let ReadOutcome::Spanning(scores) = outcome {
            for (_, ll) in &scores {
                sink += *ll;
            }
        }
    }
    sink
}
