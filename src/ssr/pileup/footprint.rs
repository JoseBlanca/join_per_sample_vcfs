//! Read-vs-locus footprint geometry — where a read sits relative to a locus,
//! derived from its mapping position and CIGAR op lengths (never its indel
//! *placement*, which the mapper gets wrong inside repeats).
//!
//! Two consumers:
//! - the fetcher's cheap **admission gate** ([`reaches_locus`]): could this read
//!   plausibly span the locus? — a footprint/clip test, no sequence scan, used to
//!   protect the reservoir budget (arch §3.1/§8.3);
//! - the delimiter's **region extraction** ([`extract_region`]): the
//!   read-coordinate slice covering the locus's embedded reference window (clips
//!   included) that the pair-HMM realigns within (`alignment.rs`).
//!
//! Lifted from the Mark-1 `triage.rs`; the Mark-1 rung-window machinery
//! (`triage_read`, `SpanningRead`, the `find_longest_stretch` content pre-probe)
//! is dropped — Mark-2 delimits by alignment, not a rung window.

use std::ops::Range;

use crate::bam::alignment_input::MappedRead;
use crate::pileup::walker::CigarOp;
use crate::ssr::types::Locus;

/// Minimum flank each side (in reference bases) a read's footprint must bracket
/// for the read to count as plausibly spanning the tract (arch §9). A
/// **calibration** placeholder (arch §14).
pub(crate) const MIN_FLANK_BP: u32 = 5;

/// The read's reference footprint, from the CIGAR's op lengths + soft-clips only
/// — never its indel *placement*. Coordinates are 0-based; `ref_end` is exclusive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Footprint {
    pub(crate) ref_start: u32,
    pub(crate) ref_end: u32,
    pub(crate) leading_clip: u32,
    pub(crate) trailing_clip: u32,
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

/// The read's reference footprint from its CIGAR + mapping position.
pub(crate) fn read_footprint(cigar: &[CigarOp], pos: u64) -> Footprint {
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
pub(crate) fn extract_region(
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

/// The fetcher's cheap coordinate-reach admission gate (arch §3.1): could this
/// read plausibly *span* the locus, decided from the footprint alone — no
/// sequence scan? The depth cap (arch §8.3) runs in the fetcher's single pass,
/// so the reads it admits should already be plausible spanning evidence;
/// otherwise the cap budget is spent on reads the worker will discard and a
/// messy high-depth locus loses its real spanning reads to eviction.
///
/// - A **soft-clipped** read is **always** admitted — the clip may carry a long
///   allele's far flank that the aligned span does not reach, and only the
///   later realignment can tell. Staying conservative here means the cap never
///   evicts a possible long allele.
/// - Otherwise, admit only a read whose aligned footprint brackets **both**
///   tract ends (`MIN_FLANK_BP` each side) — it could be spanning.
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
    use crate::ssr::types::Motif;

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
}
