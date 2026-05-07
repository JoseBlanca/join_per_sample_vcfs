//! Lazy CIGAR cursor — emits per-position events on demand from a
//! `PreparedRead`'s CIGAR. Replaces the eager `Vec<ReadEvent>` per
//! active read and the per-walker-step clones that drove the
//! super-linear scaling measured in
//! [`benches/pileup_walker_scaling.rs`](../../../benches/pileup_walker_scaling.rs).
//!
//! The cursor holds only an offset table — one entry per CIGAR
//! op. Every event-fetching call takes the read back as a
//! parameter, so the cursor stays a small self-contained struct
//! with no borrow on the read's buffers.
//!
//! Two query shapes:
//!
//! - `events_at(walker_pos, read)`: events whose **anchor** is
//!   `walker_pos`. The walker uses this once per step to identify
//!   contributing reads.
//! - `events_overlapping(lo, hi, read)`: events whose **footprint**
//!   intersects `[lo, hi)`. The open-record fold uses this to
//!   compute the haplotype string a read presents under a given
//!   record's footprint.
//!
//! The two are deliberately not equivalent — a deletion anchored
//! before the window whose deleted run reaches into it overlaps
//! the window but is not anchored at any position inside.
//!
//! Two implementation strategies, selected automatically at
//! construction based on the CIGAR's op count:
//!
//! - **Linear** (default for ≤ `BINARY_SEARCH_OP_THRESHOLD` ops):
//!   walk all ops left-to-right with an early break once we pass
//!   the query window. Constant per-call setup; cheap when only
//!   a handful of ops exist (typical Illumina with 1–5 ops).
//! - **Binary search** (for > `BINARY_SEARCH_OP_THRESHOLD` ops):
//!   `partition_point` the offset table to skip the prefix of
//!   ops that can't contribute, then walk forward. Pays a small
//!   fixed setup cost per call, amortised by skipping tens or
//!   hundreds of ops on long reads with many small indels.
//!
//! The crossover is empirical — see the bench discussion in
//! `ia/feature_implementation_plans/pileup_lazy_cigar.md` and
//! the implementation report. Auto-selection keeps the choice
//! internal to the cursor so callers don't need to reason about
//! it.
//!
//! Equivalence with the (now test-only) `decompose` reference is
//! asserted by the parity oracle in this module's test submodule
//! — both modes are exercised against the same corpus.

use super::CigarOp;
use super::PreparedRead;
use super::decompose::{ReadEvent, indel_bq_proxy_deletion, indel_bq_proxy_insertion};

/// CIGAR op-count threshold above which the cursor switches from
/// linear-walk-with-early-break to binary-search-on-offsets. The
/// value is chosen from the multi-op benchmark runs in
/// `benches/pileup_walker_scaling.rs`: at 5 ops linear was ~16 %
/// faster than binary; at 19 ops binary was ~14 % faster than
/// linear; the crossover sits between. Set conservatively so
/// reads in the ambiguous middle ground prefer linear (whose
/// per-call cost is fixed and small) over binary (whose
/// per-call partition_point setup amortises only at higher op
/// counts).
const BINARY_SEARCH_OP_THRESHOLD: usize = 16;

/// Which inner implementation the cursor dispatches its public
/// queries to. Set once at construction and never changed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum CursorMode {
    Linear,
    BinarySearch,
}

impl CursorMode {
    /// Pick a strategy based on the CIGAR's op count. Threshold
    /// at `BINARY_SEARCH_OP_THRESHOLD`.
    fn auto_select(num_ops: usize) -> Self {
        if num_ops > BINARY_SEARCH_OP_THRESHOLD {
            Self::BinarySearch
        } else {
            Self::Linear
        }
    }
}

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
    mode: CursorMode,
}

impl CigarCursor {
    /// Build the offset table by walking the cigar once. The
    /// caller's `alignment_start` seeds `offsets[0].ref_pos`.
    /// Mode is auto-selected from the CIGAR's op count via
    /// `CursorMode::auto_select`.
    pub(super) fn new(cigar: &[CigarOp], alignment_start: u32) -> Self {
        Self::with_mode(cigar, alignment_start, CursorMode::auto_select(cigar.len()))
    }

    /// Construct a cursor with an explicitly-chosen mode. Used by
    /// the parity tests so the same input can be exercised under
    /// both implementations.
    pub(super) fn with_mode(cigar: &[CigarOp], alignment_start: u32, mode: CursorMode) -> Self {
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
        Self { offsets, mode }
    }

    /// Events whose **footprint** intersects the half-open
    /// ref-position range `[lo, hi)`. A Match's footprint is the
    /// single base at its anchor; an Insertion's is also one base
    /// (the anchor base, with the inserted run sitting after it);
    /// a Deletion's footprint covers the anchor plus the
    /// `deleted_len` deleted bases — so a deletion anchored
    /// before `lo` whose deleted run reaches into `[lo, hi)` is
    /// returned.
    ///
    /// This matches the predicate the open-record fold uses
    /// today (`open_record::process_position`'s `s < rec_end &&
    /// s + span > rec_pos`) so the fold can replace the
    /// in-memory `full_window_events` filter with one cursor
    /// call.
    ///
    /// Indels are dropped under the same first/last-op rule as
    /// `decompose`: an indel anchored at position 0 (off the
    /// chromosome's start) or carried by the first/last cigar op
    /// (no flanking match to anchor against) is silently skipped.
    pub(super) fn events_overlapping(
        &self,
        lo: u32,
        hi: u32,
        read: &PreparedRead,
    ) -> Vec<ReadEvent> {
        match self.mode {
            CursorMode::Linear => self.events_overlapping_linear(lo, hi, read),
            CursorMode::BinarySearch => self.events_overlapping_binary(lo, hi, read),
        }
    }

    /// Linear-walk implementation: iterate every op left-to-right
    /// with an early break once we pass `hi`. Cheap setup, fixed
    /// cost per op visited; preferred when the CIGAR has few ops
    /// (typical Illumina with 1–5 ops). Selected when
    /// `mode == CursorMode::Linear`.
    ///
    /// The per-op match arms are **inlined manually** rather than
    /// delegated to `emit_event_for_op_overlapping`. The bench
    /// showed a ~15-20% regression on Linear-mode workloads when
    /// the per-op work was extracted, even with `#[inline]` /
    /// `#[inline(always)]` on the helper — rustc apparently
    /// doesn't fully fuse the helper into this hot path. The
    /// duplication with `emit_event_for_op_overlapping` is
    /// deliberate.
    fn events_overlapping_linear(&self, lo: u32, hi: u32, read: &PreparedRead) -> Vec<ReadEvent> {
        let mut out = Vec::new();
        if hi <= lo {
            return out;
        }
        let n_ops = read.cigar.len();
        for i in 0..n_ops {
            let off = self.offsets[i];
            if off.ref_pos > hi {
                break;
            }
            let is_first = i == 0;
            let is_last = i + 1 == n_ops;
            match read.cigar[i] {
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
                    if anchor >= hi || anchor + 1 <= lo {
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
                    let footprint_end = anchor + len + 1;
                    if anchor >= hi || footprint_end <= lo {
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

    /// Binary-search implementation: `partition_point` the offset
    /// table to skip the prefix of ops that can't contribute,
    /// then walk forward with the same upper-bound break. Used
    /// when the CIGAR has many ops; the fixed setup cost
    /// amortises over the saved op visits.
    fn events_overlapping_binary(&self, lo: u32, hi: u32, read: &PreparedRead) -> Vec<ReadEvent> {
        let mut out = Vec::new();
        if hi <= lo {
            return out;
        }
        let n_ops = read.cigar.len();
        // Start at the last op whose start ref_pos is `<= lo`.
        // Earlier ops are safe to skip: they cover only positions
        // before `lo`, and even a long deletion at an earlier op
        // cannot reach into `[lo, hi)` because the CIGAR walk
        // guarantees its footprint ends at exactly
        // `offsets[j+1].ref_pos`, which is `<= lo` by construction
        // (j+1 ≤ i_start, so its ref_pos is also `<= lo`).
        let i_start = self
            .offsets
            .partition_point(|o| o.ref_pos <= lo)
            .saturating_sub(1);
        for i in i_start..n_ops {
            let off = self.offsets[i];
            if off.ref_pos > hi {
                break;
            }
            self.emit_event_for_op_overlapping(i, lo, hi, &mut out, read);
        }
        out
    }

    /// Per-op emission used by both `events_overlapping_*`
    /// implementations. Pulled into a helper so the linear and
    /// binary-search variants only differ in the iteration
    /// strategy, not in the per-op decision logic. Marked
    /// `#[inline]` because the bench showed the function-call
    /// overhead dominating the linear path on short-CIGAR reads
    /// (~27 % slowdown on the L=150 multi-op fixture without
    /// inlining).
    #[inline(always)]
    fn emit_event_for_op_overlapping(
        &self,
        i: usize,
        lo: u32,
        hi: u32,
        out: &mut Vec<ReadEvent>,
        read: &PreparedRead,
    ) {
        let off = self.offsets[i];
        let n_ops = read.cigar.len();
        let is_first = i == 0;
        let is_last = i + 1 == n_ops;
        match read.cigar[i] {
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
                    return;
                }
                let anchor = off.ref_pos - 1;
                // Insertion footprint is 1 base (the anchor),
                // identical to Match's overlap test.
                if anchor >= hi || anchor + 1 <= lo {
                    return;
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
                    return;
                }
                let anchor = off.ref_pos - 1;
                // Deletion footprint is `len + 1` bases:
                // anchor + the deleted run.
                let footprint_end = anchor + len + 1;
                if anchor >= hi || footprint_end <= lo {
                    return;
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

    /// Events whose **anchor** is exactly `walker_pos`. The
    /// walker uses this at every step to identify the events the
    /// active reads contribute at the current position. At most
    /// one Match plus one indel can share an anchor (the last
    /// base of one M op equals the anchor of the I/D that
    /// follows), so the returned vec carries at most 2.
    ///
    /// Note this is **not** equivalent to
    /// `events_overlapping(pos, pos + 1)`: the latter would
    /// include a deletion anchored before `pos` whose footprint
    /// reaches `pos`, which the walker has already processed at
    /// an earlier step.
    pub(super) fn events_at(&self, walker_pos: u32, read: &PreparedRead) -> Vec<ReadEvent> {
        match self.mode {
            CursorMode::Linear => self.events_at_linear(walker_pos, read),
            CursorMode::BinarySearch => self.events_at_binary(walker_pos, read),
        }
    }

    /// Linear-walk implementation: iterate every op left-to-right
    /// with an early break once we pass `walker_pos + 1`. The
    /// per-op work is just an integer-range check; selected when
    /// the CIGAR has few ops so the constant per-call partition
    /// setup of the binary version isn't worth paying.
    ///
    /// Per-op match arms are inlined manually — see the same
    /// note on `events_overlapping_linear`.
    fn events_at_linear(&self, walker_pos: u32, read: &PreparedRead) -> Vec<ReadEvent> {
        let mut out = Vec::new();
        let n_ops = read.cigar.len();
        // Anything that can anchor at `walker_pos` lives in an op
        // whose start `ref_pos` is at most `walker_pos + 1` (the
        // op containing `walker_pos` for a Match, or the op whose
        // first ref position is one past `walker_pos` for an
        // indel anchored at the previous op's end).
        let break_threshold = walker_pos.saturating_add(1);
        for i in 0..n_ops {
            let off = self.offsets[i];
            if off.ref_pos > break_threshold {
                break;
            }
            let is_first = i == 0;
            let is_last = i + 1 == n_ops;
            match read.cigar[i] {
                CigarOp::Match(len) | CigarOp::SeqMatch(len) | CigarOp::SeqMismatch(len) => {
                    let op_lo = off.ref_pos;
                    let op_hi = op_lo + len;
                    if walker_pos >= op_lo && walker_pos < op_hi {
                        let read_off = off.read_pos + (walker_pos - op_lo);
                        out.push(ReadEvent::Match {
                            ref_pos: walker_pos,
                            base: read.seq[read_off as usize],
                            bq_baq: read.bq_baq[read_off as usize],
                        });
                    }
                }
                CigarOp::Insertion(len) => {
                    if is_first || is_last || off.ref_pos <= 1 {
                        continue;
                    }
                    let anchor = off.ref_pos - 1;
                    if anchor != walker_pos {
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
                    if anchor != walker_pos {
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

    /// Binary-search implementation: `partition_point` the offset
    /// table to land directly on the candidate op(s), inspecting
    /// at most one Match op plus the (typically one or two)
    /// indel ops whose anchor lands at `walker_pos`. Selected
    /// when the CIGAR has many ops so skipping the prefix is
    /// worth the per-call setup.
    fn events_at_binary(&self, walker_pos: u32, read: &PreparedRead) -> Vec<ReadEvent> {
        let mut out = Vec::new();
        let n_ops = read.cigar.len();
        if n_ops == 0 {
            return out;
        }
        // `i_after` is the first op index whose start `ref_pos`
        // is greater than `walker_pos`. The op at `i_after - 1`
        // is the unique candidate for a Match anchored at
        // `walker_pos`. The op(s) at `i_after` onwards (sharing
        // `ref_pos == walker_pos + 1`) are candidates for indels
        // anchored at `walker_pos`. Indels can stack (an
        // Insertion immediately followed by a Deletion both
        // anchor at the same ref position because an Insertion
        // does not advance ref_pos), so the forward-walk
        // continues while ref_pos matches.
        let i_after = self.offsets.partition_point(|o| o.ref_pos <= walker_pos);

        // Match candidate.
        if i_after > 0 {
            let i_lo = i_after - 1;
            // i_lo could equal n_ops if walker_pos is past every
            // real op start (i.e., walker_pos >= the post-last
            // sentinel offset). The sentinel has no matching op
            // in `read.cigar`.
            if i_lo < n_ops
                && let CigarOp::Match(len) | CigarOp::SeqMatch(len) | CigarOp::SeqMismatch(len) =
                    read.cigar[i_lo]
            {
                let off = self.offsets[i_lo];
                let op_lo = off.ref_pos;
                let op_hi = op_lo + len;
                if walker_pos >= op_lo && walker_pos < op_hi {
                    let read_off = off.read_pos + (walker_pos - op_lo);
                    out.push(ReadEvent::Match {
                        ref_pos: walker_pos,
                        base: read.seq[read_off as usize],
                        bq_baq: read.bq_baq[read_off as usize],
                    });
                }
            }
        }

        // Indel candidates: walk forward from `i_after` while
        // `ref_pos == walker_pos + 1`.
        let target_ref_pos = walker_pos.saturating_add(1);
        let mut j = i_after;
        while j < n_ops {
            let off = self.offsets[j];
            if off.ref_pos != target_ref_pos {
                break;
            }
            let is_first = j == 0;
            let is_last = j + 1 == n_ops;
            match read.cigar[j] {
                CigarOp::Insertion(len) => {
                    if !is_first && !is_last && off.ref_pos > 1 {
                        let anchor = off.ref_pos - 1;
                        let read_off = off.read_pos as usize;
                        let seq = read.seq[read_off..read_off + len as usize].to_vec();
                        let bq_proxy = indel_bq_proxy_insertion(&read.bq_baq, read_off, len);
                        out.push(ReadEvent::Insertion {
                            anchor_ref_pos: anchor,
                            seq,
                            bq_proxy,
                        });
                    }
                }
                CigarOp::Deletion(len) => {
                    if !is_first && !is_last && off.ref_pos > 1 {
                        let anchor = off.ref_pos - 1;
                        let read_off = off.read_pos as usize;
                        let bq_proxy = indel_bq_proxy_deletion(&read.bq_baq, read_off, len);
                        out.push(ReadEvent::Deletion {
                            anchor_ref_pos: anchor,
                            deleted_len: len,
                            bq_proxy,
                        });
                    }
                }
                CigarOp::Match(_)
                | CigarOp::SeqMatch(_)
                | CigarOp::SeqMismatch(_)
                | CigarOp::Skip(_)
                | CigarOp::SoftClip(_)
                | CigarOp::HardClip(_)
                | CigarOp::Padding(_) => {}
            }
            j += 1;
        }
        out
    }

    /// Per-op emission used by `events_at_linear`. Mirrors
    /// `emit_event_for_op_overlapping`'s shape but checks for an
    /// exact anchor match against `walker_pos` rather than a
    /// footprint overlap. Marked `#[inline]` for the same reason
    /// — the linear path is hot and the function call overhead
    /// dominates without inlining.
    #[inline(always)]
    fn emit_event_for_op_at(
        &self,
        i: usize,
        walker_pos: u32,
        out: &mut Vec<ReadEvent>,
        read: &PreparedRead,
    ) {
        let off = self.offsets[i];
        let n_ops = read.cigar.len();
        let is_first = i == 0;
        let is_last = i + 1 == n_ops;
        match read.cigar[i] {
            CigarOp::Match(len) | CigarOp::SeqMatch(len) | CigarOp::SeqMismatch(len) => {
                let op_lo = off.ref_pos;
                let op_hi = op_lo + len;
                if walker_pos >= op_lo && walker_pos < op_hi {
                    let read_off = off.read_pos + (walker_pos - op_lo);
                    out.push(ReadEvent::Match {
                        ref_pos: walker_pos,
                        base: read.seq[read_off as usize],
                        bq_baq: read.bq_baq[read_off as usize],
                    });
                }
            }
            CigarOp::Insertion(len) => {
                if is_first || is_last || off.ref_pos <= 1 {
                    return;
                }
                let anchor = off.ref_pos - 1;
                if anchor != walker_pos {
                    return;
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
                    return;
                }
                let anchor = off.ref_pos - 1;
                if anchor != walker_pos {
                    return;
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

    /// Parity oracle: `events_overlapping(0, MAX)` returns every
    /// event `decompose` would emit, in the same order. Footprint
    /// overlap with `[0, MAX)` is the universal case — every
    /// event qualifies — so this is the cleanest single-shot
    /// equivalence test. Run under **both** modes so each
    /// implementation is checked against `decompose`.
    #[test]
    fn cursor_matches_decompose_on_pattern_corpus() {
        for mode in [CursorMode::Linear, CursorMode::BinarySearch] {
            for (name, read) in pattern_corpus() {
                let cursor = CigarCursor::with_mode(&read.cigar, read.alignment_start, mode);
                let cursor_events = cursor.events_overlapping(0, u32::MAX, &read);
                let decompose_events = decompose(&read);
                assert_eq!(
                    cursor_events, decompose_events,
                    "mode {mode:?}, pattern {name}: cursor and decompose disagree",
                );
            }
        }
    }

    /// `events_at(pos)` returns exactly the events whose anchor
    /// is `pos`, equivalent to filtering decompose's output by
    /// anchor. Run under both modes.
    #[test]
    fn events_at_picks_only_walker_pos_events() {
        for mode in [CursorMode::Linear, CursorMode::BinarySearch] {
            for (name, read) in pattern_corpus() {
                let cursor = CigarCursor::with_mode(&read.cigar, read.alignment_start, mode);
                let decompose_events = decompose(&read);
                for walker_pos in read.alignment_start..=read.alignment_end {
                    let want: Vec<ReadEvent> = decompose_events
                        .iter()
                        .filter(|e| e.anchor_pos() == walker_pos)
                        .cloned()
                        .collect();
                    let got = cursor.events_at(walker_pos, &read);
                    assert_eq!(
                        got, want,
                        "mode {mode:?}, pattern {name}, walker_pos {walker_pos}: cursor disagrees with decompose",
                    );
                }
            }
        }
    }

    /// `events_overlapping(lo, hi)` returns exactly the events
    /// whose footprint intersects `[lo, hi)` — the same predicate
    /// `open_record::process_position` applies to
    /// `full_window_events` today. Run under both modes.
    #[test]
    fn events_overlapping_subsets_correctly() {
        for mode in [CursorMode::Linear, CursorMode::BinarySearch] {
            for (name, read) in pattern_corpus() {
                let cursor = CigarCursor::with_mode(&read.cigar, read.alignment_start, mode);
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
                            // Empty windows can never overlap.
                            if lo >= hi {
                                return false;
                            }
                            let s = e.anchor_pos();
                            let span = e.footprint_span();
                            s < hi && s + span > lo
                        })
                        .cloned()
                        .collect();
                    let got = cursor.events_overlapping(lo, hi, &read);
                    assert_eq!(
                        got, want,
                        "mode {mode:?}, pattern {name}, window [{lo}, {hi}): cursor disagrees with decompose",
                    );
                }
            }
        }
    }

    /// Auto-selection picks `Linear` for short CIGARs and
    /// `BinarySearch` for long ones. Threshold is
    /// `BINARY_SEARCH_OP_THRESHOLD`; we test just below and just
    /// above to pin the boundary.
    #[test]
    fn auto_select_picks_mode_by_op_count() {
        // 1 op — well below threshold → Linear.
        let one_op = vec![CigarOp::Match(150)];
        let cursor = CigarCursor::new(&one_op, 1);
        assert_eq!(cursor.mode, CursorMode::Linear);

        // Exactly threshold ops (≤ threshold uses Linear).
        let at_threshold: Vec<CigarOp> = (0..BINARY_SEARCH_OP_THRESHOLD)
            .map(|_| CigarOp::Match(50))
            .collect();
        let cursor = CigarCursor::new(&at_threshold, 1);
        assert_eq!(
            cursor.mode,
            CursorMode::Linear,
            "exactly {BINARY_SEARCH_OP_THRESHOLD} ops should still pick Linear (> threshold, not >=)",
        );

        // One above threshold → BinarySearch.
        let above_threshold: Vec<CigarOp> = (0..BINARY_SEARCH_OP_THRESHOLD + 1)
            .map(|_| CigarOp::Match(50))
            .collect();
        let cursor = CigarCursor::new(&above_threshold, 1);
        assert_eq!(cursor.mode, CursorMode::BinarySearch);
    }

    /// Deletion anchored just before a window: footprint reaches
    /// in, so `events_overlapping` includes it; `events_at`
    /// (anchor-only) does not. Pinning the semantic split since
    /// the open-record fold relies on it.
    #[test]
    fn deletion_anchored_before_window_overlaps_but_not_anchored() {
        // 5M5D5M starting at 100. Deletion anchor at 104, footprint
        // spans [104, 110) (anchor + 5 deleted bases). A window
        // [105, 108) does not contain the anchor 104 but the
        // deletion's footprint reaches in.
        let read = make_read(
            vec![CigarOp::Match(5), CigarOp::Deletion(5), CigarOp::Match(5)],
            100,
            b"ACGTAACGTA",
            &[30; 10],
        );
        let cursor = CigarCursor::new(&read.cigar, read.alignment_start);

        let overlapping = cursor.events_overlapping(105, 108, &read);
        assert!(
            overlapping.iter().any(|e| matches!(
                e,
                ReadEvent::Deletion {
                    anchor_ref_pos: 104,
                    ..
                }
            )),
            "deletion should be returned by events_overlapping when its footprint reaches the window",
        );
        let anchored = cursor.events_at(104, &read);
        assert!(
            anchored.iter().any(|e| matches!(
                e,
                ReadEvent::Deletion {
                    anchor_ref_pos: 104,
                    ..
                }
            )),
            "deletion should be returned by events_at(104)",
        );
        let anchored_inside = cursor.events_at(105, &read);
        assert!(
            !anchored_inside
                .iter()
                .any(|e| matches!(e, ReadEvent::Deletion { .. })),
            "deletion must NOT be returned by events_at(105) — that anchor is the deletion's interior, not its anchor",
        );
    }

    /// The cursor is stateless under repeated queries — calling
    /// the same method twice with the same args returns identical
    /// output, and order of mixed calls does not matter.
    #[test]
    fn cursor_queries_are_stateless() {
        let (_, read) = pattern_corpus().into_iter().next().unwrap();
        let cursor = CigarCursor::new(&read.cigar, read.alignment_start);
        let a = cursor.events_overlapping(0, u32::MAX, &read);
        let b = cursor.events_at(read.alignment_start, &read);
        let c = cursor.events_overlapping(0, u32::MAX, &read);
        let d = cursor.events_at(read.alignment_start, &read);
        assert_eq!(a, c);
        assert_eq!(b, d);
    }
}
