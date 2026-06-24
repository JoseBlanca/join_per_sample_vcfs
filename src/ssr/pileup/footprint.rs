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

/// The read length the CIGAR implies — the sum of its read-consuming op lengths
/// (M / I / S / `=` / X). A well-formed record has `cigar_read_len(cigar) ==
/// seq.len()`; [`fetch_locus_reads`](super::fetch_reads::fetch_locus_reads)
/// drops reads where the two disagree (a truncated/malformed record), because
/// the delimiter slices `seq` / `qual` by ranges derived from `seq.len()`.
pub(crate) fn cigar_read_len(cigar: &[CigarOp]) -> usize {
    cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n)
            | CigarOp::Insertion(n)
            | CigarOp::SoftClip(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n) => *n as usize,
            _ => 0,
        })
        .sum()
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
    // Clamp both ends into `[0, read_len]`: `ref_to_read` is bounded by the
    // CIGAR's read-consumption, so a CIGAR consistent with `seq.len()` keeps
    // `r_start <= read_len`. The clamp is defense-in-depth — the fetcher already
    // drops length-inconsistent records — so an out-of-range index can never
    // reach the `seq`/`qual` slice in `process_locus`.
    let r_start = r_start.min(read_len);
    r_start..r_end.min(read_len).max(r_start)
}

/// Whether a delimited tract looks **truncated by the extraction window** rather than by
/// the read genuinely ending — the long-allele recovery trigger (plan
/// `ssr_pileup_long_allele_window_recovery.md`). A read carrying an allele longer than the
/// reference, whose mapper represented it as all-Match (no soft-clip, no insertion), has
/// its far flank pushed out of the ref-sized window; the pair-HMM then can't anchor that
/// flank and collapses the allele toward the reference length.
///
/// A side is **window-bounded** (more read bases exist there, the window just stopped short)
/// when [`extract_region`] did not reach the read edge on that side: `region.start > 0`
/// (left) / `region.end < read_len` (right). The tract is suspicious when a window-bounded
/// side carries fewer flank bytes than the locus declares — the displaced tract ate into the
/// flank. `tract` is **region-relative** (the [`Delimited::Region`](super::alignment::Delimited)
/// span). A side that reached the read edge is not flagged: there are no more read bases to
/// recover, so it is the genuine "allele ≥ read length" case the delimiter already reports.
pub(crate) fn flank_truncated(
    region: &Range<usize>,
    tract: &Range<usize>,
    read_len: usize,
    left_flank_len: usize,
    right_flank_len: usize,
) -> bool {
    let region_len = region.end - region.start;
    let left_flank_bytes = tract.start; // region-relative tract start = bytes before the tract
    let right_flank_bytes = region_len - tract.end;
    let left_window_bounded = region.start > 0;
    let right_window_bounded = region.end < read_len;
    (left_window_bounded && left_flank_bytes < left_flank_len)
        || (right_window_bounded && right_flank_bytes < right_flank_len)
}

/// Widen a read-coordinate region by one full reference flank on each side, clamped to the
/// read (plan `ssr_pileup_long_allele_window_recovery.md`). Extending in **read** rather
/// than reference coordinates sidesteps the unreliable CIGAR of a mis-aligned long-allele
/// read; the margin is the locus's own flank length, which the catalog's bundle guarantee
/// (`bundle_threshold ≥ flank_bp` of clean sequence past each flank) makes safe — the wider
/// window cannot reach a neighbouring tract.
pub(crate) fn widen_region(
    region: Range<usize>,
    read_len: usize,
    left_flank_len: usize,
    right_flank_len: usize,
) -> Range<usize> {
    let start = region.start.saturating_sub(left_flank_len);
    let end = (region.end + right_flank_len).min(read_len);
    start..end
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
    fn cigar_read_len_sums_read_consuming_ops() {
        // M / I / S / = / X consume read bases; D / N / H / P do not.
        assert_eq!(cigar_read_len(&[CigarOp::Match(30)]), 30);
        assert_eq!(
            cigar_read_len(&[
                CigarOp::SoftClip(3),
                CigarOp::Match(10),
                CigarOp::Insertion(2),
            ]),
            15
        );
        // A deletion extends the reference span but consumes no read base.
        assert_eq!(
            cigar_read_len(&[CigarOp::Match(5), CigarOp::Deletion(4), CigarOp::Match(5)]),
            10
        );
        // A hard clip is not present in the read sequence.
        assert_eq!(
            cigar_read_len(&[CigarOp::HardClip(4), CigarOp::Match(8)]),
            8
        );
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
    fn extract_region_maps_window_across_an_internal_deletion() {
        // M8 D2 M10 from pos 11: aligned ref [10,30) brackets the window [10,28);
        // the 2 bp deletion means the read covers those 18 ref bases in 16 read
        // bases, so the window maps to read[0..16].
        let cigar = vec![CigarOp::Match(8), CigarOp::Deletion(2), CigarOp::Match(10)];
        let fp = read_footprint(&cigar, 11);
        assert_eq!(extract_region(&cigar, fp, 18, &locus6()), 0..16);
    }

    #[test]
    fn extract_region_maps_window_across_an_internal_insertion() {
        // M8 I2 M10 from pos 11: aligned ref [10,28) == the window; the 2 inserted
        // read bases fall inside it, so the window maps to the whole read[0..20].
        let cigar = vec![CigarOp::Match(8), CigarOp::Insertion(2), CigarOp::Match(10)];
        let fp = read_footprint(&cigar, 11);
        assert_eq!(extract_region(&cigar, fp, 20, &locus6()), 0..20);
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

    // ── long-allele window recovery: flank_truncated / widen_region ──

    #[test]
    fn flank_truncated_false_for_a_full_flanked_tract() {
        // region = GGGGGG | CACACA | TTTTTT (6 bp flanks both sides), read continues
        // past both ends → window-bounded but flanks complete → not suspicious.
        assert!(!flank_truncated(&(5..33), &(6..12), 42, 6, 6));
    }

    #[test]
    fn flank_truncated_true_when_the_far_flank_is_eaten_by_a_long_allele() {
        // The collapse case: region 28 bp, tract delimited to [10,26) leaving only 2 bp of
        // right flank (region_len - 26 = 2 < 6), and the right side is window-bounded
        // (region.end 33 < read_len 42) → suspicious.
        assert!(flank_truncated(&(5..33), &(10..26), 42, 6, 6));
    }

    #[test]
    fn flank_truncated_false_when_the_short_side_reached_the_read_end() {
        // Same short right flank, but the region ends at the read end (33 == read_len):
        // no more read bases to recover → the genuine allele-≥-read-length case, not a
        // window truncation → not flagged.
        assert!(!flank_truncated(&(5..33), &(10..26), 33, 6, 6));
    }

    #[test]
    fn widen_region_extends_one_flank_each_side_clamped_to_the_read() {
        // Interior region widens by the full flank on both sides.
        assert_eq!(widen_region(20..40, 100, 6, 6), 14..46);
        // Clamps at the read bounds (start floors at 0, end caps at read_len).
        assert_eq!(widen_region(3..38, 40, 6, 6), 0..40);
    }
}
