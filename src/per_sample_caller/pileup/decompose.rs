//! CIGAR decomposition: walk a `PreparedRead`'s CIGAR once and
//! produce a sorted list of `ReadEvent`s — the read's per-position
//! contributions to the pileup. See `ia/specs/pileup_walker.md`
//! §"Read decomposition" for the algorithm, edge cases, and rules.

use super::CigarOp;
use super::PreparedRead;

/// One contribution from a read at a single reference position.
/// Variants carry the 1-based reference position they apply to;
/// the walker collects events at `walker_pos` by matching that field.
#[derive(Debug, Clone, PartialEq)]
pub enum ReadEvent {
    /// `M`/`=`/`X`: at `ref_pos`, the read has `base` aligned with
    /// per-base BAQ-capped quality `bq_baq`.
    Match { ref_pos: u32, base: u8, bq_baq: u8 },
    /// Insertion of `seq` anchored at `anchor_ref_pos` (the
    /// reference position immediately before the inserted run).
    /// `bq_proxy` is the min over an `l + 2` BAQ window centred on
    /// the indel.
    Insertion {
        anchor_ref_pos: u32,
        seq: Vec<u8>,
        bq_proxy: u8,
    },
    /// Deletion of `deleted_len` reference bases starting at
    /// `anchor_ref_pos + 1`, anchored at `anchor_ref_pos`.
    /// `bq_proxy` is the min over an `l + 2` BAQ window centred on
    /// the indel (clamped at the read ends).
    Deletion {
        anchor_ref_pos: u32,
        deleted_len: u32,
        bq_proxy: u8,
    },
}

impl ReadEvent {
    /// The reference position at which this event is collected by
    /// the walker (the *anchor* for indels, the position itself
    /// for matches).
    pub fn anchor_pos(&self) -> u32 {
        match self {
            ReadEvent::Match { ref_pos, .. } => *ref_pos,
            ReadEvent::Insertion { anchor_ref_pos, .. } => *anchor_ref_pos,
            ReadEvent::Deletion { anchor_ref_pos, .. } => *anchor_ref_pos,
        }
    }

    /// REF-span footprint length: how many reference bases the
    /// event covers when forming an open record. SNP/INS produce
    /// span 1 (the anchor base only); DEL produces `len + 1`
    /// (anchor + deleted bases).
    pub fn footprint_span(&self) -> u32 {
        match self {
            ReadEvent::Match { .. } => 1,
            ReadEvent::Insertion { .. } => 1,
            ReadEvent::Deletion { deleted_len, .. } => *deleted_len + 1,
        }
    }
}

/// Bases of flanking sequence on each side of an indel that go
/// into the BQ-proxy `min` window. `INDEL_BQ_PROXY_PAD = 1`
/// yields freebayes' `l + 2` total window length. Source:
/// `freebayes/src/AlleleParser.cpp:1626` and the architecture
/// spec's indel-BQ-proxy commitment.
const INDEL_BQ_PROXY_PAD: usize = 1;

/// Walk `read.cigar` once and emit the read's events in
/// ref-position order.
pub fn decompose(read: &PreparedRead) -> Vec<ReadEvent> {
    let mut events: Vec<ReadEvent> = Vec::with_capacity(read.cigar.len());
    let mut ref_pos: u32 = read.alignment_start;
    let mut read_pos: usize = 0;
    let n_ops = read.cigar.len();

    for (op_idx, op) in read.cigar.iter().enumerate() {
        let is_first = op_idx == 0;
        let is_last = op_idx + 1 == n_ops;

        match *op {
            CigarOp::Match(len) | CigarOp::SeqMatch(len) | CigarOp::SeqMismatch(len) => {
                for k in 0..len {
                    events.push(ReadEvent::Match {
                        ref_pos: ref_pos + k,
                        base: read.seq[read_pos + k as usize],
                        bq_baq: read.bq_baq[read_pos + k as usize],
                    });
                }
                ref_pos += len;
                read_pos += len as usize;
            }
            CigarOp::Insertion(len) => {
                // Drop the indel only when (a) it is the first or
                // last CIGAR op (no flanking — freebayes' rule), or
                // (b) its anchor would fall before position 1 (off
                // the chromosome's start edge).
                if !is_first && !is_last && ref_pos > 1 {
                    let anchor = ref_pos - 1;
                    let seq = read.seq[read_pos..read_pos + len as usize].to_vec();
                    let bq = indel_bq_proxy_insertion(&read.bq_baq, read_pos, len);
                    events.push(ReadEvent::Insertion {
                        anchor_ref_pos: anchor,
                        seq,
                        bq_proxy: bq,
                    });
                }
                read_pos += len as usize;
                // ref_pos unchanged: insertions do not consume reference.
            }
            CigarOp::Deletion(len) => {
                if !is_first && !is_last && ref_pos > 1 {
                    let anchor = ref_pos - 1;
                    let bq = indel_bq_proxy_deletion(&read.bq_baq, read_pos, len);
                    events.push(ReadEvent::Deletion {
                        anchor_ref_pos: anchor,
                        deleted_len: len,
                        bq_proxy: bq,
                    });
                }
                ref_pos += len;
                // read_pos unchanged: deletions do not consume read bases.
            }
            CigarOp::Skip(len) => {
                // RNA splice — produces no event. Both ends produce
                // matches independently before/after.
                ref_pos += len;
            }
            CigarOp::SoftClip(len) => {
                // Clipped read bases — drop them.
                read_pos += len as usize;
            }
            CigarOp::HardClip(_) => {
                // H bases are not in `seq` to begin with — nothing
                // to advance.
            }
            CigarOp::Padding(_) => {
                // P: silent deletion from a padded reference. No-op
                // here — almost no real aligner emits this.
            }
        }
    }

    events
}

/// Indel BQ proxy for an insertion: min over a `l + 2 * PAD`
/// window centred on the inserted run, edge-clamped at read ends.
/// Centring matches freebayes' insertion code path
/// (window covers the inserted bases plus ±PAD flanking read
/// bases — `freebayes/src/AlleleParser.cpp:1714`).
fn indel_bq_proxy_insertion(bq_baq: &[u8], read_pos: usize, op_len: u32) -> u8 {
    let lo = read_pos.saturating_sub(INDEL_BQ_PROXY_PAD);
    let hi = (read_pos + op_len as usize + INDEL_BQ_PROXY_PAD).min(bq_baq.len());
    bq_baq[lo..hi].iter().copied().min().unwrap_or(0)
}

/// Indel BQ proxy for a deletion: deleted bases have no read-side
/// BQ, so we take min over an `l + 2 * PAD` window centred on the
/// read-position cursor (which sits at the base immediately after
/// the deletion's preceding `M` end). Edge-clamped at read ends.
/// Source: `freebayes/src/AlleleParser.cpp:1638-1659`.
fn indel_bq_proxy_deletion(bq_baq: &[u8], read_pos: usize, deleted_len: u32) -> u8 {
    // For a deletion of length `l`, freebayes uses `L = l + 2`
    // bases centred on the cursor `read_pos - L/2`. Edge-clamped.
    let l = deleted_len as usize + 2 * INDEL_BQ_PROXY_PAD;
    let centre = read_pos.saturating_sub(l / 2);
    let lo = centre;
    let hi = (centre + l).min(bq_baq.len());
    let lo = lo.min(hi.saturating_sub(1));
    bq_baq[lo..hi].iter().copied().min().unwrap_or(0)
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    fn make_read(
        cigar: Vec<CigarOp>,
        alignment_start: u32,
        seq: &[u8],
        qual: &[u8],
    ) -> PreparedRead {
        // Compute alignment_end from the cigar (consumed reference span).
        let mut span = 0u32;
        for op in &cigar {
            match *op {
                CigarOp::Match(n)
                | CigarOp::SeqMatch(n)
                | CigarOp::SeqMismatch(n)
                | CigarOp::Deletion(n)
                | CigarOp::Skip(n) => span += n,
                _ => {}
            }
        }
        let alignment_end = if span == 0 {
            alignment_start
        } else {
            alignment_start + span - 1
        };
        PreparedRead {
            chrom_id: 0,
            alignment_start,
            alignment_end,
            cigar,
            seq: seq.to_vec(),
            bq_baq: qual.to_vec(),
            mq_log_err: -3.0,
            is_reverse_strand: false,
            qname: Arc::from("r"),
            is_first_mate: true,
            has_mate: false,
        }
    }

    #[test]
    fn pure_match_yields_one_match_event_per_position() {
        let r = make_read(vec![CigarOp::Match(3)], 100, b"ACG", &[30, 30, 30]);
        let events = decompose(&r);
        assert_eq!(events.len(), 3);
        assert_eq!(
            events[0],
            ReadEvent::Match {
                ref_pos: 100,
                base: b'A',
                bq_baq: 30
            }
        );
        assert_eq!(
            events[2],
            ReadEvent::Match {
                ref_pos: 102,
                base: b'G',
                bq_baq: 30
            }
        );
    }

    #[test]
    fn snp_pattern_via_eq_x() {
        let r = make_read(
            vec![
                CigarOp::SeqMatch(2),
                CigarOp::SeqMismatch(1),
                CigarOp::SeqMatch(2),
            ],
            100,
            b"ACGTA",
            &[30; 5],
        );
        let events = decompose(&r);
        assert_eq!(events.len(), 5);
        assert_eq!(events[2].anchor_pos(), 102);
        match &events[2] {
            ReadEvent::Match { base, .. } => assert_eq!(*base, b'G'),
            _ => panic!("expected Match"),
        }
    }

    #[test]
    fn insertion_in_middle_is_anchored_at_preceding_m() {
        // CIGAR: 5M2I3M starting at 100.
        // M positions 100..104, then 2 inserted bases ("GG"),
        // anchored at 104, then M positions 105..107.
        let r = make_read(
            vec![CigarOp::Match(5), CigarOp::Insertion(2), CigarOp::Match(3)],
            100,
            b"ACGTAGGCGT",
            &[30; 10],
        );
        let events = decompose(&r);
        // 5 M + 1 I + 3 M = 9 events
        assert_eq!(events.len(), 9);
        // Event 5 should be the insertion.
        match &events[5] {
            ReadEvent::Insertion {
                anchor_ref_pos,
                seq,
                ..
            } => {
                assert_eq!(*anchor_ref_pos, 104);
                assert_eq!(seq, b"GG");
            }
            other => panic!("expected Insertion, got {:?}", other),
        }
        // Last M event should be at ref_pos 107.
        assert_eq!(events.last().unwrap().anchor_pos(), 107);
    }

    #[test]
    fn deletion_in_middle_is_anchored_at_preceding_m() {
        // CIGAR: 5M3D5M starting at 100.
        // Deletion of bases 105..107 anchored at 104. ref_pos
        // continues at 108 after the D.
        let r = make_read(
            vec![CigarOp::Match(5), CigarOp::Deletion(3), CigarOp::Match(5)],
            100,
            b"ACGTAACGTA",
            &[30; 10],
        );
        let events = decompose(&r);
        // 5 M + 1 D + 5 M = 11 events
        assert_eq!(events.len(), 11);
        match &events[5] {
            ReadEvent::Deletion {
                anchor_ref_pos,
                deleted_len,
                ..
            } => {
                assert_eq!(*anchor_ref_pos, 104);
                assert_eq!(*deleted_len, 3);
            }
            other => panic!("expected Deletion, got {:?}", other),
        }
        assert_eq!(events.last().unwrap().anchor_pos(), 112);
    }

    #[test]
    fn first_op_indel_is_dropped_but_rest_of_read_keeps_events() {
        // CIGAR: 2I5M starting at 100. The leading I has no flanking
        // M on the left, so it is dropped. The trailing 5M still
        // contributes 5 Match events at positions 100..104.
        let r = make_read(
            vec![CigarOp::Insertion(2), CigarOp::Match(5)],
            100,
            b"GGACGTA",
            &[30; 7],
        );
        let events = decompose(&r);
        assert_eq!(events.len(), 5);
        assert_eq!(events[0].anchor_pos(), 100);
        assert_eq!(events.last().unwrap().anchor_pos(), 104);
        assert!(
            !events
                .iter()
                .any(|e| matches!(e, ReadEvent::Insertion { .. }))
        );
    }

    #[test]
    fn last_op_indel_is_dropped() {
        // CIGAR: 5M3D starting at 100. Trailing D dropped.
        let r = make_read(
            vec![CigarOp::Match(5), CigarOp::Deletion(3)],
            100,
            b"ACGTA",
            &[30; 5],
        );
        let events = decompose(&r);
        assert_eq!(events.len(), 5);
        assert!(
            !events
                .iter()
                .any(|e| matches!(e, ReadEvent::Deletion { .. }))
        );
    }

    #[test]
    fn soft_clip_then_indel_at_alignment_start_one_drops_indel() {
        // CIGAR: 1S2I5M with alignment_start = 1. The leading 1S
        // doesn't advance ref_pos, so the I sits with ref_pos = 1,
        // anchor would be 0 — off the chromosome. Drop the I.
        let r = make_read(
            vec![
                CigarOp::SoftClip(1),
                CigarOp::Insertion(2),
                CigarOp::Match(5),
            ],
            1,
            b"XGGACGTA",
            &[30; 8],
        );
        let events = decompose(&r);
        assert_eq!(
            events.len(),
            5,
            "5 Match events from the trailing M, no insertion"
        );
        assert!(
            !events
                .iter()
                .any(|e| matches!(e, ReadEvent::Insertion { .. }))
        );
        assert_eq!(events[0].anchor_pos(), 1);
        assert_eq!(events.last().unwrap().anchor_pos(), 5);
    }

    #[test]
    fn soft_clip_then_indel_away_from_chrom_start_keeps_indel() {
        // Same CIGAR shape but alignment_start = 100. ref_pos = 100
        // when we hit the I; anchor = 99 — valid. Indel survives.
        let r = make_read(
            vec![
                CigarOp::SoftClip(1),
                CigarOp::Insertion(2),
                CigarOp::Match(5),
            ],
            100,
            b"XGGACGTA",
            &[30; 8],
        );
        let events = decompose(&r);
        // 1 I + 5 M = 6 events.
        assert_eq!(events.len(), 6);
        assert!(matches!(
            events[0],
            ReadEvent::Insertion {
                anchor_ref_pos: 99,
                ..
            }
        ));
    }

    #[test]
    fn n_op_skips_reference_and_produces_no_event() {
        // CIGAR: 3M100N3M starting at 100. The read covers ref
        // positions 100,101,102 and 203,204,205 — 6 events total,
        // 100 reference bases skipped silently.
        let r = make_read(
            vec![CigarOp::Match(3), CigarOp::Skip(100), CigarOp::Match(3)],
            100,
            b"ACGTAC",
            &[30; 6],
        );
        let events = decompose(&r);
        assert_eq!(events.len(), 6);
        assert_eq!(events[0].anchor_pos(), 100);
        assert_eq!(events[2].anchor_pos(), 102);
        assert_eq!(events[3].anchor_pos(), 203);
        assert_eq!(events.last().unwrap().anchor_pos(), 205);
    }

    #[test]
    fn padding_and_hard_clip_are_no_ops() {
        let r = make_read(
            vec![
                CigarOp::HardClip(5),
                CigarOp::Match(3),
                CigarOp::Padding(2),
                CigarOp::Match(2),
            ],
            100,
            b"ACGTA",
            &[30; 5],
        );
        let events = decompose(&r);
        // P consumes neither ref nor query, so it doesn't even
        // shift positions; H is already absent from `seq`.
        assert_eq!(events.len(), 5);
        assert_eq!(events[0].anchor_pos(), 100);
        assert_eq!(events.last().unwrap().anchor_pos(), 104);
    }

    #[test]
    fn indel_bq_proxy_takes_min_of_window() {
        // 5M2I3M with the inserted bases at low quality flanked by
        // higher quality matches. min over the (l+2)=4-byte window
        // centred on the inserted bases is the lower of the two
        // inserted-base BQs.
        let mut qual = vec![30, 30, 30, 30, 30, 5, 7, 30, 30, 30];
        // read_pos cursor when we hit the I is 5 (5 M consumed),
        // window covers read_pos-1..read_pos+2+1 = [4, 8) → quals
        // 30, 5, 7, 30 → min 5.
        let r = make_read(
            vec![CigarOp::Match(5), CigarOp::Insertion(2), CigarOp::Match(3)],
            100,
            b"ACGTAXXACG",
            &qual,
        );
        let events = decompose(&r);
        let ins = events
            .iter()
            .find(|e| matches!(e, ReadEvent::Insertion { .. }))
            .unwrap();
        match ins {
            ReadEvent::Insertion { bq_proxy, .. } => assert_eq!(*bq_proxy, 5),
            _ => unreachable!(),
        }
        // exercise the quality buffer reference so the linter
        // doesn't think it's unused after the assertion.
        qual.clear();
    }

    #[test]
    fn deletion_bq_proxy_uses_flanking_bases() {
        // 5M2D5M; deleted bases have no read-side BQ.
        // Window centre = read_pos - L/2 = 5 - 2 = 3.
        // Window = read_pos[3..7] = quals 30, 30, 4, 30 → min 4.
        let qual = vec![30, 30, 30, 30, 30, 30, 30, 30, 30, 30];
        // Tweak position 5 to a low value so we know which slice
        // got picked.
        let mut qual = qual;
        qual[5] = 4;
        let r = make_read(
            vec![CigarOp::Match(5), CigarOp::Deletion(2), CigarOp::Match(5)],
            100,
            b"ACGTAACGTA",
            &qual,
        );
        let events = decompose(&r);
        let del = events
            .iter()
            .find(|e| matches!(e, ReadEvent::Deletion { .. }))
            .unwrap();
        match del {
            ReadEvent::Deletion { bq_proxy, .. } => assert_eq!(*bq_proxy, 4),
            _ => unreachable!(),
        }
    }

    #[test]
    fn consecutive_indels_produce_two_events() {
        // 3M2I3D2M: I anchored at 102, D anchored at 102 (same ref
        // position because the I doesn't consume reference).
        let r = make_read(
            vec![
                CigarOp::Match(3),
                CigarOp::Insertion(2),
                CigarOp::Deletion(3),
                CigarOp::Match(2),
            ],
            100,
            b"ACGGGCG",
            &[30; 7],
        );
        let events = decompose(&r);
        // 3 M + 1 I + 1 D + 2 M = 7 events.
        assert_eq!(events.len(), 7);
        // Find the I and D and check their anchors.
        let ins = events
            .iter()
            .find(|e| matches!(e, ReadEvent::Insertion { .. }));
        let del = events
            .iter()
            .find(|e| matches!(e, ReadEvent::Deletion { .. }));
        match ins {
            Some(ReadEvent::Insertion { anchor_ref_pos, .. }) => assert_eq!(*anchor_ref_pos, 102),
            _ => panic!("expected an Insertion"),
        }
        match del {
            Some(ReadEvent::Deletion { anchor_ref_pos, .. }) => assert_eq!(*anchor_ref_pos, 102),
            _ => panic!("expected a Deletion"),
        }
    }
}
