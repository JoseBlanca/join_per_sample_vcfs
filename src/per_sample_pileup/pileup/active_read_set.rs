//! Active-set bookkeeping: holds reads currently overlapping the
//! walker's position, plus a secondary `read_id → vec_index` map
//! for O(1) mate lookup. See `ia/specs/pileup_walker.md` §"Active
//! read table" and §"Active-set bookkeeping".

use ahash::AHashMap;

use super::PreparedRead;
use super::chain_id_allocator::{ChainId, ChainIdAllocator};
use super::cigar_cursor::CigarCursor;
use super::errors::WalkerError;

/// One read currently in the walker's active set. The `read` field
/// is the owned `PreparedRead`; `cursor` carries the per-op offset
/// table the walker queries on demand for events at the current
/// position or in an open record's footprint window.
#[derive(Debug, Clone)]
pub struct ActiveRead {
    pub read_id: u32,
    pub read: PreparedRead,
    pub cursor: CigarCursor,
    pub chain_id: ChainId,
    /// `Some(other.read_id)` when the partner mate has also been
    /// admitted. Filled by the active-set admission code on the
    /// second-mate path.
    pub mate_read_id: Option<u32>,
}

#[derive(Debug)]
pub struct ActiveReads {
    /// Primary container. Iteration order is admission order; we
    /// rely on this for deterministic record-formation tiebreaks
    /// (smaller `read_id` first).
    reads: Vec<ActiveRead>,
    /// Secondary index: `read_id → index in reads`. Maintained in
    /// lockstep with `reads`.
    by_read_id: AHashMap<u32, usize>,
    /// Monotonically-increasing local id, allocated on admission.
    /// Wrap is not a concern at any realistic input size (`u32`
    /// covers 4 billion reads).
    next_read_id: u32,
}

impl ActiveReads {
    pub fn new() -> Self {
        Self {
            reads: Vec::new(),
            by_read_id: AHashMap::new(),
            next_read_id: 0,
        }
    }

    /// Reset all in-flight state. Called at chromosome boundaries.
    pub fn reset(&mut self) {
        self.reads.clear();
        self.by_read_id.clear();
        // next_read_id keeps advancing — read ids stay unique
        // across the whole run, which makes log messages
        // unambiguous.
    }

    pub fn is_empty(&self) -> bool {
        self.reads.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &ActiveRead> {
        self.reads.iter()
    }

    /// Look up a read by `read_id` via the secondary index. Used
    /// by `open_record::process_position` to query each
    /// contributor's `CigarCursor` for window events at the open
    /// record's footprint, and by the active-set tests.
    pub fn get_by_read_id(&self, read_id: u32) -> Option<&ActiveRead> {
        let idx = *self.by_read_id.get(&read_id)?;
        Some(&self.reads[idx])
    }

    /// Test-only helper: number of reads currently in the set.
    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.reads.len()
    }

    /// Admit a read: assign a `read_id`, decompose its CIGAR,
    /// allocate a phase-chain id via `chain_ids`, and (if this is
    /// the second mate of a pair) cross-link `mate_read_id` with
    /// the first mate.
    pub fn admit(
        &mut self,
        read: PreparedRead,
        chain_ids: &mut ChainIdAllocator,
    ) -> Result<u32, WalkerError> {
        let read_id = self.next_read_id;
        // `checked_add` guards against the 4-billion-read case
        // where `next_read_id` would otherwise wrap silently in
        // release mode and collide with an active-set entry via
        // the `by_read_id` index. Implausible at routine
        // workloads but cheap to defend against.
        self.next_read_id =
            self.next_read_id
                .checked_add(1)
                .ok_or_else(|| WalkerError::Internal {
                    detail: format!("active-set next_read_id overflowed u32 at {read_id}"),
                    qname: read.qname.to_string(),
                    chrom_id: read.chrom_id,
                    pos: read.alignment_start,
                })?;

        let qname_for_register = read.qname.clone();
        let (chain_id, partner_read_id) = chain_ids.allocate_for_read(&read)?;

        // Tell the chain-id allocator this fresh first-mate's read_id
        // (no-op when this is a second mate or solo read, since the
        // qname isn't in `pending_mates` in those cases).
        chain_ids.register_first_mate_read_id(&qname_for_register, read_id);

        let cursor = CigarCursor::new(&read.cigar, read.alignment_start);

        let active = ActiveRead {
            read_id,
            read,
            cursor,
            chain_id,
            mate_read_id: partner_read_id,
        };

        let new_index = self.reads.len();
        self.reads.push(active);
        self.by_read_id.insert(read_id, new_index);

        // If this was the second mate, also stitch the back-link on
        // the first mate's `ActiveRead` so either side can find the
        // other.
        if let Some(partner) = partner_read_id
            && let Some(partner_idx) = self.by_read_id.get(&partner).copied()
        {
            self.reads[partner_idx].mate_read_id = Some(read_id);
        }

        Ok(read_id)
    }

    /// Drop reads whose alignment ends before `walker_pos`. For
    /// each dropped read, release its phase-chain bookkeeping and
    /// update the secondary index.
    pub fn expire_passed(
        &mut self,
        walker_pos: u32,
        chain_ids: &mut ChainIdAllocator,
    ) -> Result<(), WalkerError> {
        // Iterate from the end so swap_remove indices stay valid.
        let mut i = self.reads.len();
        while i > 0 {
            i -= 1;
            if self.reads[i].read.alignment_end < walker_pos {
                // Single-fold-per-(record, read) is enforced by
                // `OpenPileupRecord.folded_reads`, not by per-read
                // state here — the lazy CIGAR cursor is stateless.
                //
                // The `pending_mates` entry (if any) stays alive
                // here: pair tracking is governed by
                // `mate_lookup_window`, not by active-set
                // residence. The chain-id allocator's
                // `evict_stale_pending` walks the map on every
                // subsequent `allocate_for_read` and drops entries
                // whose `seen_at` is past the window.
                let chrom_id = self.reads[i].read.chrom_id;
                chain_ids.note_read_exit(chrom_id, walker_pos)?;
                self.swap_remove(i);
            }
        }
        Ok(())
    }

    /// Drop every read unconditionally — used at chromosome
    /// boundaries and end-of-input. Releases each read's
    /// active-count slot. `walker_pos` is the walker's current
    /// locus, used only to populate `WalkerError::Internal`'s
    /// context fields if a release ever hits the panic-free path
    /// (Mi2).
    pub fn flush_all(
        &mut self,
        chain_ids: &mut ChainIdAllocator,
        walker_pos: u32,
    ) -> Result<(), WalkerError> {
        while let Some(active) = self.reads.pop() {
            // The secondary index is rebuilt fresh after a flush.
            // `pending_mates` is cleaned up by the chain-id allocator's
            // `reset()` call (which `walker::flush_chromosome_into`
            // makes right after this method); we don't touch it
            // here.
            let chrom_id = active.read.chrom_id;
            chain_ids.note_read_exit(chrom_id, walker_pos)?;
        }
        self.by_read_id.clear();
        Ok(())
    }

    fn swap_remove(&mut self, idx: usize) {
        let removed_read_id = self.reads[idx].read_id;
        self.by_read_id.remove(&removed_read_id);
        // swap_remove moves the last element (if any) into `idx`.
        let last_idx = self.reads.len() - 1;
        if idx != last_idx {
            let moved_read_id = self.reads[last_idx].read_id;
            self.by_read_id.insert(moved_read_id, idx);
        }
        self.reads.swap_remove(idx);
    }
}

impl Default for ActiveReads {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::per_sample_pileup::pileup::CigarOp;

    fn solo_read(qname: &str, chrom_id: u32, alignment_start: u32, span: u32) -> PreparedRead {
        PreparedRead {
            chrom_id,
            alignment_start,
            alignment_end: alignment_start + span - 1,
            cigar: vec![CigarOp::Match(span)],
            seq: vec![b'A'; span as usize],
            bq_baq: vec![30; span as usize],
            mq_log_err: -3.0,
            is_reverse_strand: false,
            qname: Arc::from(qname),
            mate_role: super::super::MateRole::Solo,
            adaptor_boundary: None,
        }
    }

    fn paired_read(
        qname: &str,
        is_first_mate: bool,
        alignment_start: u32,
        span: u32,
    ) -> PreparedRead {
        let mut r = solo_read(qname, 0, alignment_start, span);
        r.mate_role = if is_first_mate {
            super::super::MateRole::FirstOfPair
        } else {
            super::super::MateRole::SecondOfPair
        };
        r
    }

    #[test]
    fn admit_assigns_increasing_read_ids() {
        let mut s = ActiveReads::new();
        let mut a = ChainIdAllocator::new();
        let r0 = s.admit(solo_read("a", 0, 100, 50), &mut a).unwrap();
        let r1 = s.admit(solo_read("b", 0, 110, 50), &mut a).unwrap();
        let r2 = s.admit(solo_read("c", 0, 120, 50), &mut a).unwrap();
        assert_eq!((r0, r1, r2), (0, 1, 2));
    }

    #[test]
    fn secondary_index_maps_read_id_to_correct_entry() {
        let mut s = ActiveReads::new();
        let mut a = ChainIdAllocator::new();
        s.admit(solo_read("a", 0, 100, 50), &mut a).unwrap();
        let r1 = s.admit(solo_read("b", 0, 110, 50), &mut a).unwrap();
        let entry = s.get_by_read_id(r1).expect("must find by id");
        assert_eq!(entry.read.qname.as_ref(), "b");
    }

    #[test]
    fn expire_passed_drops_only_reads_behind_walker() {
        let mut s = ActiveReads::new();
        let mut a = ChainIdAllocator::new();
        let _r0 = s.admit(solo_read("short", 0, 100, 10), &mut a).unwrap(); // ends 109
        let r1 = s.admit(solo_read("long", 0, 100, 200), &mut a).unwrap(); // ends 299
        // Move walker to 150 — short is past, long still active.
        s.expire_passed(150, &mut a).unwrap();
        assert_eq!(s.len(), 1);
        assert!(s.get_by_read_id(r1).is_some());
    }

    #[test]
    fn expire_uses_swap_remove_and_keeps_index_consistent() {
        // Admit three reads. Expire the middle one (walker between
        // the second's end and the third's end). Verify the index
        // still resolves the remaining two correctly.
        let mut s = ActiveReads::new();
        let mut a = ChainIdAllocator::new();
        let r0 = s.admit(solo_read("a", 0, 100, 1000), &mut a).unwrap(); // ends 1099
        let _r1 = s.admit(solo_read("b", 0, 100, 50), &mut a).unwrap(); // ends 149
        let r2 = s.admit(solo_read("c", 0, 100, 1000), &mut a).unwrap(); // ends 1099

        s.expire_passed(200, &mut a).unwrap();
        assert_eq!(s.len(), 2, "middle read expired");
        assert_eq!(
            s.get_by_read_id(r0).map(|r| r.read.qname.as_ref()),
            Some("a")
        );
        assert_eq!(
            s.get_by_read_id(r2).map(|r| r.read.qname.as_ref()),
            Some("c")
        );
    }

    #[test]
    fn paired_reads_get_mate_read_id_cross_links() {
        let mut s = ActiveReads::new();
        let mut a = ChainIdAllocator::new();
        let m1 = s.admit(paired_read("p", true, 100, 50), &mut a).unwrap();
        let m2 = s.admit(paired_read("p", false, 130, 50), &mut a).unwrap();
        // Both mates should now reference each other.
        assert_eq!(s.get_by_read_id(m1).unwrap().mate_read_id, Some(m2));
        assert_eq!(s.get_by_read_id(m2).unwrap().mate_read_id, Some(m1));
    }

    #[test]
    fn flush_all_drops_every_active_read() {
        let mut s = ActiveReads::new();
        let mut a = ChainIdAllocator::new();
        for i in 0..5 {
            s.admit(solo_read(&format!("r{i}"), 0, 100, 50), &mut a)
                .unwrap();
        }
        assert_eq!(s.len(), 5);
        s.flush_all(&mut a, 200).unwrap();
        assert!(s.is_empty());
    }
}
