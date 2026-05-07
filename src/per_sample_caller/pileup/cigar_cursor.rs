//! Lazy CIGAR cursor — emits per-position events on demand from a
//! `PreparedRead`'s CIGAR. Designed to replace the eager
//! `Vec<ReadEvent>` per active read and the per-walker-step clones
//! that drive the O(L²) hot path measured in
//! [`benches/pileup_walker_scaling.rs`](../../../benches/pileup_walker_scaling.rs).
//!
//! Commit 1 of the migration in
//! `ia/feature_implementation_plans/pileup_lazy_cigar.md`: the
//! cursor lands here behind the existing `decompose`, with parity
//! tests asserting it reproduces `decompose`'s output. The walker
//! still uses `decompose` until commit 2.
//!
//! The cursor itself holds only an offset table — one entry per
//! CIGAR op. Every event-fetching call takes the read back as a
//! parameter, so the cursor stays a small self-contained struct
//! with no borrow on the read's buffers.

use super::CigarOp;
use super::PreparedRead;
use super::decompose::{ReadEvent, indel_bq_proxy_deletion, indel_bq_proxy_insertion};

/// Reference and read offsets at the start of one CIGAR op. The
/// table built from these is the cursor's only persistent state.
#[derive(Debug, Clone, Copy)]
struct OpOffset {
    /// 1-based reference position at the op's start.
    ref_pos: u32,
    /// 0-based read offset at the op's start.
    read_pos: u32,
}

/// Pre-walked CIGAR for a single read. Holds an `OpOffset` per
/// CIGAR op plus a sentinel at index `cigar.len()` carrying the
/// post-last-op (ref_pos, read_pos). For typical short and long
/// reads the table is well under 100 entries.
#[derive(Debug, Clone)]
pub(super) struct CigarCursor {
    offsets: Vec<OpOffset>,
}

impl CigarCursor {
    /// Build the offset table by walking the cigar once. The
    /// caller's `alignment_start` seeds `offsets[0].ref_pos`.
    pub(super) fn new(cigar: &[CigarOp], alignment_start: u32) -> Self {
        let mut offsets = Vec::with_capacity(cigar.len() + 1);
        let mut ref_pos = alignment_start;
        let mut read_pos: u32 = 0;
        for op in cigar {
            offsets.push(OpOffset { ref_pos, read_pos });
            match *op {
                CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                    ref_pos += n;
                    read_pos += n;
                }
                CigarOp::Insertion(n) => {
                    read_pos += n;
                }
                CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                    ref_pos += n;
                }
                CigarOp::SoftClip(n) => {
                    read_pos += n;
                }
                CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
            }
        }
        offsets.push(OpOffset { ref_pos, read_pos });
        Self { offsets }
    }

    /// Events whose anchor falls in the half-open ref-position
    /// range `[lo, hi)`. Walks the cigar left-to-right and emits
    /// in cigar order — same order `decompose` produces. Indels
    /// are dropped under the same first/last-op rule as
    /// `decompose`: an indel anchored at position 0 (off the
    /// chromosome's start) or carried by the first/last cigar op
    /// (no flanking match to anchor against) is silently skipped.
    pub(super) fn events_in_window(&self, lo: u32, hi: u32, read: &PreparedRead) -> Vec<ReadEvent> {
        let mut out = Vec::new();
        if hi <= lo {
            return out;
        }
        let n_ops = read.cigar.len();
        for (i, op) in read.cigar.iter().enumerate() {
            let off = self.offsets[i];
            let is_first = i == 0;
            let is_last = i + 1 == n_ops;
            match *op {
                CigarOp::Match(len) | CigarOp::SeqMatch(len) | CigarOp::SeqMismatch(len) => {
                    let op_lo = off.ref_pos;
                    let op_hi = op_lo + len;
                    let emit_lo = op_lo.max(lo);
                    let emit_hi = op_hi.min(hi);
                    if emit_lo < emit_hi {
                        let span = emit_hi - emit_lo;
                        for k in 0..span {
                            let ref_pos = emit_lo + k;
                            let read_off = off.read_pos + (ref_pos - op_lo);
                            out.push(ReadEvent::Match {
                                ref_pos,
                                base: read.seq[read_off as usize],
                                bq_baq: read.bq_baq[read_off as usize],
                            });
                        }
                    }
                }
                CigarOp::Insertion(len) => {
                    if is_first || is_last || off.ref_pos <= 1 {
                        continue;
                    }
                    let anchor = off.ref_pos - 1;
                    if anchor < lo || anchor >= hi {
                        continue;
                    }
                    let read_off = off.read_pos as usize;
                    let seq = read.seq[read_off..read_off + len as usize].to_vec();
                    let bq_proxy = indel_bq_proxy_insertion(&read.bq_baq, read_off, len);
                    out.push(ReadEvent::Insertion {
                        anchor_ref_pos: anchor,
                        seq,
                        bq_proxy,
                    });
                }
                CigarOp::Deletion(len) => {
                    if is_first || is_last || off.ref_pos <= 1 {
                        continue;
                    }
                    let anchor = off.ref_pos - 1;
                    if anchor < lo || anchor >= hi {
                        continue;
                    }
                    let read_off = off.read_pos as usize;
                    let bq_proxy = indel_bq_proxy_deletion(&read.bq_baq, read_off, len);
                    out.push(ReadEvent::Deletion {
                        anchor_ref_pos: anchor,
                        deleted_len: len,
                        bq_proxy,
                    });
                }
                CigarOp::Skip(_)
                | CigarOp::SoftClip(_)
                | CigarOp::HardClip(_)
                | CigarOp::Padding(_) => {}
            }
        }
        out
    }

    /// Convenience: events anchored exactly at `walker_pos`.
    /// Equivalent to `events_in_window(walker_pos, walker_pos + 1, read)`.
    /// At most one Match plus one indel can sit at any single
    /// anchor (e.g. last base of one M op = anchor of the I/D
    /// that follows), so the returned vec carries at most 2.
    pub(super) fn events_at(&self, walker_pos: u32, read: &PreparedRead) -> Vec<ReadEvent> {
        self.events_in_window(walker_pos, walker_pos + 1, read)
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::per_sample_caller::pileup::CigarOp;
    use crate::per_sample_caller::pileup::decompose::decompose;

    fn make_read(
        cigar: Vec<CigarOp>,
        alignment_start: u32,
        seq: &[u8],
        qual: &[u8],
    ) -> PreparedRead {
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

    /// All the CIGAR shapes covered by `decompose`'s own tests,
    /// plus an explicit consecutive-indels case. Used by the
    /// parity oracle below to assert the cursor reproduces
    /// `decompose`'s output byte-for-byte.
    fn pattern_corpus() -> Vec<(&'static str, PreparedRead)> {
        let dq_with_low_at_5 = {
            let mut q = vec![30u8; 10];
            q[5] = 4;
            q
        };
        vec![
            (
                "pure_match",
                make_read(vec![CigarOp::Match(3)], 100, b"ACG", &[30, 30, 30]),
            ),
            (
                "seq_match_mismatch",
                make_read(
                    vec![
                        CigarOp::SeqMatch(2),
                        CigarOp::SeqMismatch(1),
                        CigarOp::SeqMatch(2),
                    ],
                    100,
                    b"ACGTA",
                    &[30; 5],
                ),
            ),
            (
                "insertion_in_middle",
                make_read(
                    vec![CigarOp::Match(5), CigarOp::Insertion(2), CigarOp::Match(3)],
                    100,
                    b"ACGTAGGCGT",
                    &[30; 10],
                ),
            ),
            (
                "deletion_in_middle",
                make_read(
                    vec![CigarOp::Match(5), CigarOp::Deletion(3), CigarOp::Match(5)],
                    100,
                    b"ACGTAACGTA",
                    &[30; 10],
                ),
            ),
            (
                "first_op_indel_dropped",
                make_read(
                    vec![CigarOp::Insertion(2), CigarOp::Match(5)],
                    100,
                    b"GGACGTA",
                    &[30; 7],
                ),
            ),
            (
                "last_op_indel_dropped",
                make_read(
                    vec![CigarOp::Match(5), CigarOp::Deletion(3)],
                    100,
                    b"ACGTA",
                    &[30; 5],
                ),
            ),
            (
                "soft_clip_indel_at_chrom_start_dropped",
                make_read(
                    vec![
                        CigarOp::SoftClip(1),
                        CigarOp::Insertion(2),
                        CigarOp::Match(5),
                    ],
                    1,
                    b"XGGACGTA",
                    &[30; 8],
                ),
            ),
            (
                "soft_clip_indel_away_from_start_kept",
                make_read(
                    vec![
                        CigarOp::SoftClip(1),
                        CigarOp::Insertion(2),
                        CigarOp::Match(5),
                    ],
                    100,
                    b"XGGACGTA",
                    &[30; 8],
                ),
            ),
            (
                "n_skip",
                make_read(
                    vec![CigarOp::Match(3), CigarOp::Skip(100), CigarOp::Match(3)],
                    100,
                    b"ACGTAC",
                    &[30; 6],
                ),
            ),
            (
                "hard_clip_padding",
                make_read(
                    vec![
                        CigarOp::HardClip(5),
                        CigarOp::Match(3),
                        CigarOp::Padding(2),
                        CigarOp::Match(2),
                    ],
                    100,
                    b"ACGTA",
                    &[30; 5],
                ),
            ),
            (
                "indel_bq_proxy_insertion",
                make_read(
                    vec![CigarOp::Match(5), CigarOp::Insertion(2), CigarOp::Match(3)],
                    100,
                    b"ACGTAXXACG",
                    &[30, 30, 30, 30, 30, 5, 7, 30, 30, 30],
                ),
            ),
            (
                "indel_bq_proxy_deletion",
                make_read(
                    vec![CigarOp::Match(5), CigarOp::Deletion(2), CigarOp::Match(5)],
                    100,
                    b"ACGTAACGTA",
                    &dq_with_low_at_5,
                ),
            ),
            (
                "consecutive_indels",
                make_read(
                    vec![
                        CigarOp::Match(3),
                        CigarOp::Insertion(2),
                        CigarOp::Deletion(3),
                        CigarOp::Match(2),
                    ],
                    100,
                    b"ACGGGCG",
                    &[30; 7],
                ),
            ),
        ]
    }

    /// Parity oracle: events_in_window covering the entire read
    /// produces the same Vec<ReadEvent> as `decompose` for every
    /// pattern in the corpus.
    #[test]
    fn cursor_matches_decompose_on_pattern_corpus() {
        for (name, read) in pattern_corpus() {
            let cursor = CigarCursor::new(&read.cigar, read.alignment_start);
            let cursor_events = cursor.events_in_window(0, u32::MAX, &read);
            let decompose_events = decompose(&read);
            assert_eq!(
                cursor_events, decompose_events,
                "pattern {name}: cursor and decompose disagree (cursor={cursor_events:?}, decompose={decompose_events:?})",
            );
        }
    }

    /// `events_at(pos)` returns exactly the events whose anchor
    /// is `pos`, equivalent to filtering decompose's output by
    /// anchor.
    #[test]
    fn events_at_picks_only_walker_pos_events() {
        for (name, read) in pattern_corpus() {
            let cursor = CigarCursor::new(&read.cigar, read.alignment_start);
            let decompose_events = decompose(&read);
            // Walk a generous range so we exercise both
            // event-bearing and silent positions.
            for walker_pos in read.alignment_start..=read.alignment_end {
                let want: Vec<ReadEvent> = decompose_events
                    .iter()
                    .filter(|e| e.anchor_pos() == walker_pos)
                    .cloned()
                    .collect();
                let got = cursor.events_at(walker_pos, &read);
                assert_eq!(
                    got, want,
                    "pattern {name}, walker_pos {walker_pos}: cursor disagrees with decompose",
                );
            }
        }
    }

    /// `events_in_window(lo, hi)` returns exactly the events whose
    /// anchor falls in `[lo, hi)`.
    #[test]
    fn events_in_window_subsets_correctly() {
        for (name, read) in pattern_corpus() {
            let cursor = CigarCursor::new(&read.cigar, read.alignment_start);
            let decompose_events = decompose(&read);
            let mid = read.alignment_start + (read.alignment_end - read.alignment_start) / 2;
            let windows = [
                (read.alignment_start, read.alignment_start + 1),
                (mid, mid + 1),
                (read.alignment_start, mid),
                (mid, read.alignment_end + 1),
                (0, u32::MAX),
                // Empty window: always yields nothing.
                (mid, mid),
            ];
            for (lo, hi) in windows {
                let want: Vec<ReadEvent> = decompose_events
                    .iter()
                    .filter(|e| {
                        let a = e.anchor_pos();
                        a >= lo && a < hi
                    })
                    .cloned()
                    .collect();
                let got = cursor.events_in_window(lo, hi, &read);
                assert_eq!(
                    got, want,
                    "pattern {name}, window [{lo}, {hi}): cursor disagrees with decompose",
                );
            }
        }
    }

    /// The cursor is stateless under repeated queries — calling
    /// `events_in_window` or `events_at` twice with the same args
    /// returns identical output, and order of mixed calls does
    /// not matter.
    #[test]
    fn cursor_queries_are_stateless() {
        let (_, read) = pattern_corpus().into_iter().next().unwrap();
        let cursor = CigarCursor::new(&read.cigar, read.alignment_start);
        let a = cursor.events_in_window(0, u32::MAX, &read);
        let b = cursor.events_at(read.alignment_start, &read);
        let c = cursor.events_in_window(0, u32::MAX, &read);
        let d = cursor.events_at(read.alignment_start, &read);
        assert_eq!(a, c);
        assert_eq!(b, d);
    }
}
