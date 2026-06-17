//! Read triage — classify reads near a locus and prepare spanning ones for the
//! pair-HMM (arch §3, realign-everything revision §2).
//!
//! Everything downstream trusts this stage's sorting. Per the realign-everything
//! design (§2 revision), triage does **not** run a fast/slow gate or trust the
//! CIGAR's indel placement: it classifies *coverage* from the read's footprint
//! (mapping position + clip lengths), and for a spanning read extracts the
//! locus-region bases and centres the candidate window with the content
//! pre-probe. Every spanning read is then realigned (§5/§6); the direct-count
//! fast path ([`super::count_repeats`]) is the deferred, measured optimization.
//!
//! - [`find_longest_stretch`] — the content pre-probe (longest contiguous motif
//!   run + total copies, arch §3.2); centres the window and recovers soft-clipped
//!   long alleles without trusting the CIGAR.
//! - [`triage_read`] — coverage classification → region extract → window centre,
//!   over the shared `MappedRead` ([`crate::bam::alignment_input`]) the SNP pileup
//!   already yields (the sanctioned reader reuse, arch §3.1).

use std::ops::Range;

use crate::bam::alignment_input::MappedRead;
use crate::pileup::walker::CigarOp;
use crate::ssr_mark1::types::{Locus, Motif};

/// Minimum flank each side (in reference bases) a read's footprint must bracket
/// for the read to count as spanning the tract (arch §9). A **calibration**
/// placeholder (arch §14).
pub(crate) const MIN_FLANK_BP: u32 = 5;

/// Result of the content pre-probe: a cheap `O(read len)` scan of a read's own
/// repeat content (arch §3.2). Used **only** when a flank lives in a soft-clip,
/// to estimate the (possibly large) repeat count without trusting the CIGAR.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ProbeHit {
    /// Longest run of **contiguous** (adjacent, non-overlapping) motif copies.
    /// This centres the pair-HMM window (`count ± W`) — the robust count at short
    /// periods, where incidental copies elsewhere inflate the total but not the
    /// contiguous run (arch §3.2).
    pub(crate) longest_run: u32,
    /// Total motif copies anywhere in the read (non-overlapping). A junk bound /
    /// slow-path trigger (against `MIN_MOTIF_RUN_BY_PERIOD`, §10), **not** the
    /// window centre.
    pub(crate) total_copies: u32,
}

/// Scan `seq` for its repeat content against `motif`, GangSTR-style — the
/// longest contiguous run of motif copies and the total copy count
/// (port of GangSTR `realignment.cpp` `find_longest_stretch`, arch §3.2).
///
/// Copies are counted **non-overlapping** (a match at `i` advances the scan by
/// `period`); a run is a maximal stretch of copies that are immediately
/// adjacent. The longest run centres the window; the total is the junk bound.
///
/// **One deliberate deviation from the C++ source:** GangSTR's loop bound
/// (`i < len − period`) misses a copy ending exactly at the read's end; we scan
/// every valid start (`i + period ≤ len`) so the final copy is counted. The
/// thresholds this feeds (`MIN_MOTIF_RUN_BY_PERIOD`) are uncalibrated placeholders
/// (arch §14), so the off-by-one is not a parity concern.
pub(crate) fn find_longest_stretch(seq: &[u8], motif: &Motif) -> ProbeHit {
    let unit = motif.as_bytes();
    let period = unit.len();
    if period == 0 || seq.len() < period {
        return ProbeHit {
            longest_run: 0,
            total_copies: 0,
        };
    }

    let mut longest_run = 0u32;
    let mut current_run = 0u32;
    let mut total_copies = 0u32;
    let mut i = 0usize;
    while i + period <= seq.len() {
        if &seq[i..i + period] == unit {
            total_copies += 1;
            current_run += 1;
            longest_run = longest_run.max(current_run);
            i += period; // non-overlapping; the next copy must abut to extend the run
        } else {
            current_run = 0;
            i += 1;
        }
    }

    ProbeHit {
        longest_run,
        total_copies,
    }
}

// --- Read classification (realign-everything: coverage → extract → centre) ----

/// Outcome of triaging one read against one locus (arch §3.3) — how many of the
/// tract's two ends the read's footprint brackets (2/1/0). Only `Spanning` feeds
/// the likelihood; the other two are tallied for QC, then dropped.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum TriageResult {
    /// Footprint brackets the tract + `MIN_FLANK_BP` on both ends → realign.
    Spanning(SpanningRead),
    /// Brackets one end only — counted (`n_flanking`), not used in v1.
    Flanking,
    /// Brackets neither end (read buried in a long tract) — `n_frr`, not used in v1.
    InRepeat,
}

/// A spanning read ready for the pair-HMM (arch §2 revision: v1 realigns every
/// spanning read). `region` is the read-coordinate span (into `seq`/`qual`)
/// covering the locus's embedded reference window (flanks + tract, soft-clips
/// included); `observed_count` is the content pre-probe's longest contiguous run,
/// which centres the `count ± W` candidate window.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SpanningRead {
    pub(crate) region: Range<usize>,
    pub(crate) observed_count: u16,
}

/// The read's reference footprint, from the CIGAR's op lengths + soft-clips only
/// — never its indel *placement* (arch §2 revision). Coordinates are 0-based;
/// `ref_end` is exclusive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Footprint {
    ref_start: u32,
    ref_end: u32,
    leading_clip: u32,
    trailing_clip: u32,
}

/// First soft-clip length at this end (skipping a hard-clip), or `0` if the end
/// is an aligned op. A `find_map` predicate: `None` = keep scanning past a
/// hard-clip, `Some` = stop.
fn end_soft_clip(op: &CigarOp) -> Option<u32> {
    match op {
        CigarOp::HardClip(_) => None,
        CigarOp::SoftClip(n) => Some(*n),
        _ => Some(0),
    }
}

fn read_footprint(cigar: &[CigarOp], pos: u64) -> Footprint {
    let ref_start = pos.saturating_sub(1) as u32;
    let ref_span: u32 = cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n)
            | CigarOp::Deletion(n)
            | CigarOp::Skip(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n) => *n,
            _ => 0,
        })
        .sum();
    Footprint {
        ref_start,
        ref_end: ref_start + ref_span,
        leading_clip: cigar.iter().find_map(end_soft_clip).unwrap_or(0),
        trailing_clip: cigar.iter().rev().find_map(end_soft_clip).unwrap_or(0),
    }
}

/// `(left_bracketed, right_bracketed)`: does the footprint — aligned span
/// extended by the soft-clip on each side (the optimistic long-allele reach) —
/// reach ≥ `MIN_FLANK_BP` past the tract on that side? Signed arithmetic so a
/// locus near a contig end never underflows.
fn brackets(fp: Footprint, locus: &Locus) -> (bool, bool) {
    let flank = MIN_FLANK_BP as i64;
    let left = (fp.ref_start as i64 - fp.leading_clip as i64) <= (locus.start() as i64 - flank);
    let right = (fp.ref_end as i64 + fp.trailing_clip as i64) >= (locus.end() as i64 + flank);
    (left, right)
}

/// Read coordinate of a reference position `target` that lies within the aligned
/// span `[ref_start, ref_end]`. Dual-cursor CIGAR walk (the `indel_norm`
/// `count_mismatches` shape); a `target` inside a deletion maps to the read
/// position the deletion sits at.
fn ref_to_read(cigar: &[CigarOp], ref_start: u32, leading_clip: u32, target: u32) -> usize {
    let mut ref_cur = ref_start;
    let mut read_cur = leading_clip as usize;
    for op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                if target < ref_cur + n {
                    return read_cur + (target - ref_cur) as usize;
                }
                ref_cur += n;
                read_cur += *n as usize;
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                if target < ref_cur + n {
                    return read_cur;
                }
                ref_cur += n;
            }
            CigarOp::Insertion(n) => read_cur += *n as usize,
            CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
        }
    }
    read_cur
}

/// The read-coordinate span covering the locus's embedded reference window
/// (`ref_bytes`: flanks + tract). Where the window extends beyond the aligned
/// span, the whole soft-clip on that side is included — that is where a long
/// allele's extra tract + far flank live (arch §3.2), and the pair-HMM realigns
/// within the grabbed bases.
fn extract_region(
    cigar: &[CigarOp],
    fp: Footprint,
    read_len: usize,
    locus: &Locus,
) -> Range<usize> {
    let w_start = locus.ref_bytes_start();
    let w_end = w_start + locus.ref_bytes().len() as u32;

    let r_start = if w_start < fp.ref_start {
        0 // window opens left of the alignment → take the full leading clip
    } else {
        ref_to_read(cigar, fp.ref_start, fp.leading_clip, w_start)
    };
    let r_end = if w_end > fp.ref_end {
        read_len // window closes right of the alignment → take the full trailing clip
    } else {
        ref_to_read(cigar, fp.ref_start, fp.leading_clip, w_end)
    };
    r_start..r_end.min(read_len).max(r_start)
}

/// Classify a read against a locus and, if it spans, prepare it for the pair-HMM
/// (arch §2 revision — realign-everything). Coverage is decided from the
/// footprint (position + clip lengths, never indel placement); a spanning read's
/// locus-region bases are extracted and the window centred by the content
/// pre-probe.
pub(crate) fn triage_read(read: &MappedRead, locus: &Locus) -> TriageResult {
    let fp = read_footprint(&read.cigar, read.pos);
    match brackets(fp, locus) {
        (true, true) => {
            let region = extract_region(&read.cigar, fp, read.seq.len(), locus);
            let probe = find_longest_stretch(&read.seq[region.clone()], &locus.motif());
            TriageResult::Spanning(SpanningRead {
                region,
                observed_count: probe.longest_run.min(u16::MAX as u32) as u16,
            })
        }
        (false, false) => TriageResult::InRepeat,
        _ => TriageResult::Flanking,
    }
}

/// The fetcher's cheap coordinate-reach admission gate (arch §3.1): could this
/// read plausibly *span* the locus, decided from the footprint alone — no
/// sequence scan? The depth cap (arch §8.3) runs in the fetcher's single pass,
/// so the reads it admits should already be plausible spanning evidence;
/// otherwise the cap budget is spent on reads the worker will discard and a
/// messy high-depth locus loses its real spanning reads to eviction.
///
/// - A **soft-clipped** read is **always** admitted — the clip may carry a long
///   allele's far flank that the aligned span does not reach, and only the
///   content scan ([`find_longest_stretch`], in the worker) can tell. Staying
///   conservative here means the cap never evicts a possible long allele.
/// - Otherwise, admit only a read whose aligned footprint brackets **both**
///   tract ends (`MIN_FLANK_BP` each side) — it could be spanning. A non-clipped
///   read reaching one flank only (flanking) or neither (buried) clearly cannot
///   span and would just burn a reservoir slot.
///
/// This is *not* the full classification — spanning confirmation, soft-clip
/// recovery, and the flanking / in-repeat split stay in [`triage_read`] in the
/// worker. Nor is it a QC filter: QC counts are tallied over *all* reads
/// regardless of this gate (arch §3.3).
pub(crate) fn reaches_locus(read: &MappedRead, locus: &Locus) -> bool {
    let fp = read_footprint(&read.cigar, read.pos);
    if fp.leading_clip > 0 || fp.trailing_clip > 0 {
        return true;
    }
    let (left, right) = brackets(fp, locus);
    left && right
}

#[cfg(test)]
mod tests {
    use super::*;

    fn motif(bytes: &[u8]) -> Motif {
        Motif::new(bytes).unwrap()
    }

    #[test]
    fn pure_repeat_run_counts_all_copies() {
        let hit = find_longest_stretch(b"CACACACA", &motif(b"CA"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 4,
                total_copies: 4
            }
        );
    }

    #[test]
    fn flanked_repeat_counts_only_the_tract() {
        // GGG | CACACA | TTT — three contiguous CA copies, nothing in the flanks.
        let hit = find_longest_stretch(b"GGGCACACATTT", &motif(b"CA"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 3,
                total_copies: 3
            }
        );
    }

    #[test]
    fn interruption_splits_the_run_but_total_counts_both_halves() {
        // CACA | GG | CACA — two runs of 2, four copies total.
        let hit = find_longest_stretch(b"CACAGGCACA", &motif(b"CA"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 2,
                total_copies: 4
            }
        );
    }

    #[test]
    fn incidental_copies_inflate_total_not_the_contiguous_run() {
        // The arch §3.2 point at short periods: isolated CA copies bump the total
        // but the longest contiguous run stays 1, so centring on the run is robust.
        let hit = find_longest_stretch(b"CAGGGGCAGGGGCA", &motif(b"CA"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 1,
                total_copies: 3
            }
        );
    }

    #[test]
    fn homopolymer_counts_each_base() {
        let hit = find_longest_stretch(b"AAAAA", &motif(b"A"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 5,
                total_copies: 5
            }
        );
    }

    #[test]
    fn trinucleotide_run() {
        let hit = find_longest_stretch(b"CAGCAGCAGCAG", &motif(b"CAG"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 4,
                total_copies: 4
            }
        );
    }

    #[test]
    fn read_shorter_than_period_finds_nothing() {
        assert_eq!(
            find_longest_stretch(b"C", &motif(b"CA")),
            ProbeHit {
                longest_run: 0,
                total_copies: 0
            }
        );
        assert_eq!(
            find_longest_stretch(b"", &motif(b"CA")),
            ProbeHit {
                longest_run: 0,
                total_copies: 0
            }
        );
    }

    #[test]
    fn single_copy_at_the_read_end_is_counted() {
        // The deliberate fix vs GangSTR's `i < len - period` bound: a lone copy
        // exactly filling the read is found (GangSTR would miss it).
        let hit = find_longest_stretch(b"CA", &motif(b"CA"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 1,
                total_copies: 1
            }
        );
    }

    // --- Classification (footprint / coverage / extract / triage_read) --------

    /// A CA locus with 6 bp flanks (≥ MIN_FLANK_BP): GGGGGG | CACACA | TTTTTT,
    /// tract at ref [16, 22), embedded window [10, 28).
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
    fn footprint_plain_match() {
        let fp = read_footprint(&[CigarOp::Match(18)], 11);
        assert_eq!(
            fp,
            Footprint {
                ref_start: 10,
                ref_end: 28,
                leading_clip: 0,
                trailing_clip: 0
            }
        );
    }

    #[test]
    fn footprint_records_soft_clips_both_ends() {
        let fp = read_footprint(
            &[
                CigarOp::SoftClip(3),
                CigarOp::Match(15),
                CigarOp::SoftClip(4),
            ],
            11,
        );
        assert_eq!(
            fp,
            Footprint {
                ref_start: 10,
                ref_end: 25,
                leading_clip: 3,
                trailing_clip: 4
            }
        );
    }

    #[test]
    fn footprint_deletion_extends_ref_span_only() {
        // M5 D2 M5: ref span 12, read consumes 10.
        let fp = read_footprint(
            &[CigarOp::Match(5), CigarOp::Deletion(2), CigarOp::Match(5)],
            11,
        );
        assert_eq!((fp.ref_start, fp.ref_end), (10, 22));
    }

    #[test]
    fn footprint_leading_clip_skips_a_hard_clip() {
        let fp = read_footprint(
            &[
                CigarOp::HardClip(2),
                CigarOp::SoftClip(3),
                CigarOp::Match(15),
            ],
            11,
        );
        assert_eq!(fp.leading_clip, 3);
    }

    #[test]
    fn brackets_a_fully_aligned_spanning_read() {
        let fp = Footprint {
            ref_start: 10,
            ref_end: 28,
            leading_clip: 0,
            trailing_clip: 0,
        };
        assert_eq!(brackets(fp, &locus6()), (true, true));
    }

    #[test]
    fn brackets_one_side_only_is_flanking() {
        // Reaches the left flank (10 ≤ 11) but the right end (20) misses 27.
        let fp = Footprint {
            ref_start: 10,
            ref_end: 20,
            leading_clip: 0,
            trailing_clip: 0,
        };
        assert_eq!(brackets(fp, &locus6()), (true, false));
    }

    #[test]
    fn brackets_neither_side_when_buried_in_a_long_tract() {
        let fp = Footprint {
            ref_start: 17,
            ref_end: 21,
            leading_clip: 0,
            trailing_clip: 0,
        };
        assert_eq!(brackets(fp, &locus6()), (false, false));
    }

    #[test]
    fn brackets_soft_clips_extend_the_reach() {
        // Aligned [13,25) misses both flank windows, but the clips reach past them.
        let fp = Footprint {
            ref_start: 13,
            ref_end: 25,
            leading_clip: 3,
            trailing_clip: 3,
        };
        assert_eq!(brackets(fp, &locus6()), (true, true));
    }

    #[test]
    fn region_of_a_full_alignment_is_the_window() {
        let cigar = vec![CigarOp::Match(18)];
        let fp = read_footprint(&cigar, 11);
        assert_eq!(extract_region(&cigar, fp, 18, &locus6()), 0..18);
    }

    #[test]
    fn region_grabs_the_leading_clip_when_the_window_opens_left_of_the_alignment() {
        // Aligned at ref [16,28) (M12) with the left flank + units in a 6 bp clip.
        let cigar = vec![CigarOp::SoftClip(6), CigarOp::Match(12)];
        let fp = read_footprint(&cigar, 17);
        assert_eq!(extract_region(&cigar, fp, 18, &locus6()), 0..18);
    }

    #[test]
    fn triage_spanning_clean_read_centres_on_the_tract() {
        let read = mapped_read(11, vec![CigarOp::Match(18)], b"GGGGGGCACACATTTTTT");
        match triage_read(&read, &locus6()) {
            TriageResult::Spanning(s) => {
                assert_eq!(s.region, 0..18);
                assert_eq!(s.observed_count, 3); // CACACA = 3 CA copies
            }
            other => panic!("expected Spanning, got {other:?}"),
        }
    }

    #[test]
    fn triage_in_repeat_read_is_dropped() {
        // Buried mid-tract (ref [17,21)), no clips reaching a flank.
        let read = mapped_read(18, vec![CigarOp::Match(4)], b"CACA");
        assert_eq!(triage_read(&read, &locus6()), TriageResult::InRepeat);
    }

    #[test]
    fn triage_flanking_read_is_dropped() {
        // Reaches the left flank, runs off the right without bracketing it.
        let read = mapped_read(6, vec![CigarOp::Match(12)], b"GGGGGCACACAC");
        assert_eq!(triage_read(&read, &locus6()), TriageResult::Flanking);
    }

    #[test]
    fn reach_gate_admits_a_spanning_read() {
        let read = mapped_read(11, vec![CigarOp::Match(18)], b"GGGGGGCACACATTTTTT");
        assert!(reaches_locus(&read, &locus6()));
    }

    #[test]
    fn reach_gate_skips_a_one_end_unclipped_read() {
        // Brackets the left flank only (Flanking), no clip on the missing
        // right side → clearly cannot span → skipped.
        let read = mapped_read(6, vec![CigarOp::Match(12)], b"GGGGGCACACAC");
        assert!(!reaches_locus(&read, &locus6()));
    }

    #[test]
    fn reach_gate_skips_a_buried_unclipped_read() {
        // Buried mid-tract (ref [17,21)), brackets neither end, no clips → skip.
        let read = mapped_read(18, vec![CigarOp::Match(4)], b"CACA");
        assert!(!reaches_locus(&read, &locus6()));
    }

    #[test]
    fn reach_gate_always_admits_a_soft_clipped_read() {
        // Same buried aligned span as above, but a soft-clip the gate must not
        // second-guess (the clip may carry a long allele) → admitted.
        let read = mapped_read(
            18,
            vec![CigarOp::SoftClip(3), CigarOp::Match(4)],
            b"NNNCACA",
        );
        assert!(reaches_locus(&read, &locus6()));
    }

    #[test]
    fn triage_recovers_a_soft_clipped_long_allele_via_the_clip() {
        // Aligner placed the left flank + 3 tract units (M12 over ref [10,22)) and
        // soft-clipped the rest (3 more units + right flank). The trailing clip
        // extends the right reach → Spanning, and the pre-probe over the
        // clip-included region recovers the full 6-unit count — the realign-
        // everything win over trusting the CIGAR's 3.
        let read = mapped_read(
            11,
            vec![CigarOp::Match(12), CigarOp::SoftClip(12)],
            b"GGGGGGCACACACACACATTTTTT",
        );
        match triage_read(&read, &locus6()) {
            TriageResult::Spanning(s) => {
                assert_eq!(s.region, 0..24); // whole read, clip included
                assert_eq!(s.observed_count, 6);
            }
            other => panic!("expected Spanning, got {other:?}"),
        }
    }
}
