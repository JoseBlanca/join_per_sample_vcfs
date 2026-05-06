//! Phase-chain slot allocator: hands out small integer ids to
//! reads (or pairs of mates), recycles them once both mates have
//! exited the active set, and reports lifecycle markers so emitted
//! records carry `new_chains` / `expired_chains` correctly.
//!
//! The allocator is decoupled from the walker — it only knows that
//! reads enter and exit, and that paired reads share a slot. See
//! `ia/specs/pileup_walker.md` §"Phase chain slot allocator".

use std::sync::Arc;

use ahash::AHashMap;

use super::errors::WalkerError;
use super::{MATE_LOOKUP_WINDOW, PreparedRead};

/// Phase-chain slot identifier. `u16` gives ~16× headroom over the
/// default `MAX_ACTIVE_SLOTS` cap (4096); `u8`'s 256-slot ceiling
/// would risk silent overflow at higher coverages. The cap, not the
/// type, is the binding constraint.
pub type SlotId = u16;

/// Hard cap on the number of concurrently-active phase-chain slots.
/// Exceeding it is a hard error rather than silent slot reuse.
pub const MAX_ACTIVE_SLOTS: u32 = 4096;

/// State the allocator holds for a first mate whose partner has not
/// yet been admitted: which slot was issued, the first mate's
/// `read_id` so the active-set admission code can link
/// `mate_read_id`, and where the first mate was seen so we can
/// evict stale entries past the lookup window.
#[derive(Debug, Clone, Copy)]
pub struct PendingMate {
    pub chain_slot_id: SlotId,
    pub first_mate_read_id: u32,
    pub seen_at: u32,
}

#[derive(Debug)]
pub struct SlotAllocator {
    /// Recycled slot ids, kept sorted descending so `pop()` returns
    /// the smallest free id. The "lowest free id" rule keeps slot
    /// numbers compact in the typical case.
    free: Vec<SlotId>,
    /// Next never-used-before slot id; advances when `free` is empty.
    next_fresh: u32,
    /// Per-slot reference count: 1 for solo reads, 2 once both mates
    /// are admitted. The slot is released when this hits 0 again.
    slot_refcount: Vec<u8>,
    /// First mates whose partner has not yet arrived. Sized for
    /// genuinely-pending pairs only — solo reads (no `has_mate`)
    /// never enter this map.
    pending_mates: AHashMap<Arc<str>, PendingMate>,
    /// Slots that started since the previous emitted record.
    /// Drained on each emission.
    new_marks: Vec<SlotId>,
    /// Slots that ended since the previous emitted record. Drained
    /// on each emission.
    expired_marks: Vec<SlotId>,
    /// Bookkeeping for the run summary.
    counters: SlotAllocatorCounters,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct SlotAllocatorCounters {
    pub slot_allocations: u64,
    pub slot_high_water: u32,
    pub mate_lookup_evictions: u64,
}

impl SlotAllocator {
    pub fn new() -> Self {
        Self {
            free: Vec::new(),
            next_fresh: 0,
            slot_refcount: vec![0u8; MAX_ACTIVE_SLOTS as usize],
            pending_mates: AHashMap::new(),
            new_marks: Vec::new(),
            expired_marks: Vec::new(),
            counters: SlotAllocatorCounters::default(),
        }
    }

    /// Reset to the initial state. Called at chromosome boundaries
    /// (chains do not span chromosomes per
    /// `ia/specs/pileup_walker.md` §"Chromosome boundaries").
    pub fn reset(&mut self) {
        self.free.clear();
        self.next_fresh = 0;
        for slot in &mut self.slot_refcount {
            *slot = 0;
        }
        self.pending_mates.clear();
        self.new_marks.clear();
        self.expired_marks.clear();
        // counters are cumulative across the whole run, do not reset
    }

    /// Allocate a slot for an entering read.
    ///
    /// - First mate of a pair (`has_mate` set, qname not in
    ///   `pending_mates`): allocate a fresh or recycled slot and
    ///   register the qname. Refcount is pre-set to 2,
    ///   anticipating the future second mate; this guarantees the
    ///   slot stays held even if the first mate exits the active
    ///   set before the second arrives.
    /// - Second mate (qname found in `pending_mates`): reuse the
    ///   first mate's slot (same chain). Refcount is unchanged
    ///   (already at 2 from the first-mate pre-allocation).
    /// - Solo read (`has_mate` unset): allocate a slot with
    ///   refcount 1.
    ///
    /// Returns `(slot_id, first_mate_read_id_if_pairing)` — the
    /// second value is `Some` only on the second-mate path, so the
    /// caller can fill `mate_read_id` cross-references in the
    /// active set.
    pub fn allocate_for_read(
        &mut self,
        read: &PreparedRead,
    ) -> Result<(SlotId, Option<u32>), WalkerError> {
        // Garbage-collect stale pending-mate entries before any
        // allocation; otherwise a long pileup of unmatched first
        // mates can grow the map until memory pressure forces a
        // hard error.
        self.evict_stale_pending(read.alignment_start);

        if read.has_mate
            && let Some(pending) = self.pending_mates.remove(&read.qname)
        {
            // Second mate of a known pair — reuse the slot. Refcount
            // is already 2 (pre-bumped at first-mate admission), so
            // we don't bump again here.
            self.counters.slot_allocations += 1;
            return Ok((pending.chain_slot_id, Some(pending.first_mate_read_id)));
        }

        let slot = self.acquire_slot(read)?;
        // Pre-bump to 2 if a partner is expected (so the slot
        // stays held if the first mate exits before the second
        // mate arrives), otherwise 1 for solo reads.
        let initial_refcount: u8 = if read.has_mate { 2 } else { 1 };
        self.set_refcount(slot, initial_refcount);
        self.counters.slot_allocations += 1;
        self.new_marks.push(slot);

        if read.has_mate {
            self.pending_mates.insert(
                read.qname.clone(),
                PendingMate {
                    chain_slot_id: slot,
                    first_mate_read_id: 0, // caller fills via `register_first_mate_read_id`
                    seen_at: read.alignment_start,
                },
            );
        }

        Ok((slot, None))
    }

    /// Update a freshly-registered first mate's `PendingMate` entry
    /// with its assigned `read_id`. Called by the active-set
    /// admission code right after it has decided the new read's
    /// `read_id`.
    pub fn register_first_mate_read_id(&mut self, qname: &Arc<str>, first_mate_read_id: u32) {
        if let Some(entry) = self.pending_mates.get_mut(qname) {
            entry.first_mate_read_id = first_mate_read_id;
        }
    }

    /// Decrement a slot's refcount on read exit. When it reaches
    /// zero the slot is recycled and stamped into `expired_marks`.
    pub fn release_slot(&mut self, slot: SlotId) -> Result<(), WalkerError> {
        let idx = slot as usize;
        if idx >= self.slot_refcount.len() || self.slot_refcount[idx] == 0 {
            return Err(WalkerError::Internal {
                detail: format!("release_slot on non-active slot {slot}"),
                qname: String::new(),
                chrom_id: 0,
                pos: 0,
            });
        }
        self.slot_refcount[idx] -= 1;
        if self.slot_refcount[idx] == 0 {
            // Sorted-descending insert so `pop()` yields the
            // smallest free id.
            let pos = self
                .free
                .binary_search_by(|s| slot.cmp(s))
                .unwrap_or_else(|p| p);
            self.free.insert(pos, slot);
            self.expired_marks.push(slot);
        }
        Ok(())
    }

    /// Drain the pending lifecycle markers into a `(new, expired)`
    /// pair, deduplicating slots that appear in both lists (a chain
    /// allocated and released between two emissions never became
    /// observable, so we suppress it from both).
    pub fn drain_lifecycle_marks(&mut self) -> (Vec<SlotId>, Vec<SlotId>) {
        // Find the intersection (slots in both lists) and remove
        // them from both. Sort + unique on each side first so the
        // emitted lists are deterministic.
        sort_dedup(&mut self.new_marks);
        sort_dedup(&mut self.expired_marks);

        let mut new_out = Vec::new();
        let mut exp_out = Vec::new();
        let mut i = 0;
        let mut j = 0;
        while i < self.new_marks.len() && j < self.expired_marks.len() {
            match self.new_marks[i].cmp(&self.expired_marks[j]) {
                std::cmp::Ordering::Less => {
                    new_out.push(self.new_marks[i]);
                    i += 1;
                }
                std::cmp::Ordering::Greater => {
                    exp_out.push(self.expired_marks[j]);
                    j += 1;
                }
                std::cmp::Ordering::Equal => {
                    // present in both — suppress
                    i += 1;
                    j += 1;
                }
            }
        }
        new_out.extend(self.new_marks.drain(i..));
        exp_out.extend(self.expired_marks.drain(j..));
        self.new_marks.clear();
        self.expired_marks.clear();
        (new_out, exp_out)
    }

    pub fn counters(&self) -> SlotAllocatorCounters {
        self.counters
    }

    fn set_refcount(&mut self, slot: SlotId, count: u8) {
        let idx = slot as usize;
        self.slot_refcount[idx] = count;
        let active = self.active_count();
        if active > self.counters.slot_high_water {
            self.counters.slot_high_water = active;
        }
    }

    fn acquire_slot(&mut self, read: &PreparedRead) -> Result<SlotId, WalkerError> {
        if let Some(slot) = self.free.pop() {
            return Ok(slot);
        }
        if self.next_fresh >= MAX_ACTIVE_SLOTS {
            return Err(WalkerError::SlotExhausted {
                cap: MAX_ACTIVE_SLOTS,
                chrom_id: read.chrom_id,
                pos: read.alignment_start,
            });
        }
        let slot = self.next_fresh as SlotId;
        self.next_fresh += 1;
        Ok(slot)
    }

    fn active_count(&self) -> u32 {
        self.slot_refcount.iter().filter(|&&c| c > 0).count() as u32
    }

    fn evict_stale_pending(&mut self, walker_pos: u32) {
        // Collect first; remove after, so we don't iterate while
        // mutating.
        let mut stale: Vec<(Arc<str>, SlotId)> = Vec::new();
        for (qname, entry) in &self.pending_mates {
            if entry.seen_at + MATE_LOOKUP_WINDOW < walker_pos {
                stale.push((qname.clone(), entry.chain_slot_id));
            }
        }
        for (qname, slot) in stale {
            self.pending_mates.remove(&qname);
            self.counters.mate_lookup_evictions += 1;
            // The first mate's slot was pre-bumped to refcount 2
            // anticipating the partner. Now that we've given up on
            // the partner, decrement to refcount 1 so the slot
            // releases when the surviving first mate exits.
            let idx = slot as usize;
            if idx < self.slot_refcount.len() && self.slot_refcount[idx] >= 2 {
                self.slot_refcount[idx] -= 1;
            }
        }
    }
}

fn sort_dedup(v: &mut Vec<SlotId>) {
    v.sort_unstable();
    v.dedup();
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_read(qname: &str, has_mate: bool, alignment_start: u32) -> PreparedRead {
        PreparedRead {
            chrom_id: 0,
            alignment_start,
            alignment_end: alignment_start + 99,
            cigar: Vec::new(),
            seq: Vec::new(),
            bq_baq: Vec::new(),
            mq_log_err: -3.0,
            is_reverse_strand: false,
            qname: Arc::from(qname),
            is_first_mate: true,
            has_mate,
        }
    }

    #[test]
    fn solo_read_allocates_slot_and_does_not_register_qname() {
        let mut a = SlotAllocator::new();
        let read = make_read("solo", false, 100);
        let (slot, mate_id) = a.allocate_for_read(&read).expect("solo allocates");
        assert_eq!(slot, 0);
        assert_eq!(mate_id, None);
        assert!(a.pending_mates.is_empty());
    }

    #[test]
    fn first_mate_registers_then_second_mate_reuses_slot() {
        let mut a = SlotAllocator::new();
        let m1 = make_read("pair1", true, 100);
        let (slot1, mate_id1) = a.allocate_for_read(&m1).expect("first mate");
        a.register_first_mate_read_id(&m1.qname, 42);

        let m2 = make_read("pair1", true, 200);
        let (slot2, mate_id2) = a.allocate_for_read(&m2).expect("second mate");
        assert_eq!(slot1, slot2, "mates must share a slot");
        assert_eq!(mate_id1, None, "first mate has no partner registered yet");
        assert_eq!(mate_id2, Some(42), "second mate sees first mate's read_id");
        assert!(
            a.pending_mates.is_empty(),
            "pending entry consumed when second mate arrives"
        );
    }

    #[test]
    fn slot_releases_only_after_both_mates_exit() {
        let mut a = SlotAllocator::new();
        let m1 = make_read("p", true, 100);
        let (slot, _) = a.allocate_for_read(&m1).unwrap();
        let m2 = make_read("p", true, 200);
        a.allocate_for_read(&m2).unwrap();

        a.release_slot(slot).unwrap();
        let (new, expired) = a.drain_lifecycle_marks();
        assert_eq!(new, vec![slot]);
        assert!(expired.is_empty(), "slot still held by other mate");

        a.release_slot(slot).unwrap();
        let (new, expired) = a.drain_lifecycle_marks();
        assert!(new.is_empty());
        assert_eq!(expired, vec![slot]);
    }

    #[test]
    fn solo_read_emits_new_then_expired_across_two_drains() {
        // The lifecycle markers are drained per emission. Allocate,
        // drain (simulating an emission between admit and exit),
        // release, drain again — only then do we see the expired
        // mark unsuppressed.
        let mut a = SlotAllocator::new();
        let read = make_read("solo", false, 100);
        let (slot, _) = a.allocate_for_read(&read).unwrap();
        let (new, expired) = a.drain_lifecycle_marks();
        assert_eq!(new, vec![slot]);
        assert!(expired.is_empty());
        a.release_slot(slot).unwrap();
        let (new, expired) = a.drain_lifecycle_marks();
        assert!(new.is_empty());
        assert_eq!(expired, vec![slot]);
    }

    #[test]
    fn drain_lifecycle_marks_suppresses_slots_present_in_both() {
        let mut a = SlotAllocator::new();
        let r = make_read("transient", false, 100);
        let (slot, _) = a.allocate_for_read(&r).unwrap();
        a.release_slot(slot).unwrap();
        let (new, expired) = a.drain_lifecycle_marks();
        // The slot existed only between two emissions and was never
        // visible: don't surface it.
        assert!(new.is_empty(), "got {:?}", new);
        assert!(expired.is_empty(), "got {:?}", expired);
    }

    #[test]
    fn lowest_free_slot_is_reused_first() {
        let mut a = SlotAllocator::new();
        let r0 = make_read("a", false, 100);
        let r1 = make_read("b", false, 100);
        let r2 = make_read("c", false, 100);
        let (s0, _) = a.allocate_for_read(&r0).unwrap();
        let (s1, _) = a.allocate_for_read(&r1).unwrap();
        let (s2, _) = a.allocate_for_read(&r2).unwrap();
        assert_eq!((s0, s1, s2), (0, 1, 2));
        a.release_slot(s1).unwrap();
        let _ = a.drain_lifecycle_marks();
        // releasing slot 1 puts it in the free pool; the next
        // allocation should reuse it (smallest free id).
        let r3 = make_read("d", false, 100);
        let (s3, _) = a.allocate_for_read(&r3).unwrap();
        assert_eq!(s3, 1);
    }

    #[test]
    fn slot_exhausted_returns_hard_error() {
        let mut a = SlotAllocator::new();
        // Fill all slots; we don't release anything.
        for i in 0..MAX_ACTIVE_SLOTS {
            let r = make_read(&format!("r{i}"), false, 100);
            a.allocate_for_read(&r).expect("under cap");
        }
        let r = make_read("overflow", false, 100);
        let err = a.allocate_for_read(&r).expect_err("must hit cap");
        assert!(
            matches!(err, WalkerError::SlotExhausted { .. }),
            "got {err:?}"
        );
    }

    #[test]
    fn stale_pending_mates_are_evicted_after_window() {
        let mut a = SlotAllocator::new();
        let m1 = make_read("orphan", true, 100);
        let _ = a.allocate_for_read(&m1).unwrap();
        assert_eq!(a.pending_mates.len(), 1);

        // Bring the walker far past the lookup window.
        let r2 = make_read("later", false, 100 + MATE_LOOKUP_WINDOW + 1);
        a.allocate_for_read(&r2).unwrap();
        assert!(
            a.pending_mates.is_empty(),
            "orphan first mate should have been evicted"
        );
        assert_eq!(a.counters().mate_lookup_evictions, 1);
    }

    #[test]
    fn reset_clears_state_but_preserves_counters() {
        let mut a = SlotAllocator::new();
        let r = make_read("x", false, 100);
        a.allocate_for_read(&r).unwrap();
        let allocs_before = a.counters().slot_allocations;
        a.reset();
        assert_eq!(a.next_fresh, 0);
        assert!(a.pending_mates.is_empty());
        // Counters are cumulative across the whole run.
        assert_eq!(a.counters().slot_allocations, allocs_before);
    }
}
