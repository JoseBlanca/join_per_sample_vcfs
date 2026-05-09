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

/// Active-slot count at which the allocator emits a one-shot soft
/// warning, well below `MAX_ACTIVE_SLOTS`. Lets the user spot a
/// pathological-coverage region while the run is still alive instead
/// of being greeted with `SlotExhausted` and no context. 75% of the
/// hard cap leaves room for a steep climb to still be visible before
/// the run fails.
pub const HIGH_WATER_WARN_THRESHOLD: u32 = MAX_ACTIVE_SLOTS * 3 / 4;

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
    /// Slots released since the most recent `drain_lifecycle_marks`
    /// call. They cannot be reused yet — moving them straight to
    /// `free` would let a new admission hand them out within the
    /// same emission window, putting `expired[S]` and `new[S]` in
    /// one record's marks and silently merging two distinct phase
    /// chains. `drain_lifecycle_marks` migrates this list to
    /// `free` after emitting the marks for the just-closed window.
    pending_free: Vec<SlotId>,
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
    /// Set the first time the active-slot count reaches
    /// `HIGH_WATER_WARN_THRESHOLD`. Idempotent within a run — and
    /// because `reset` (called at chromosome boundaries) deliberately
    /// preserves it, the warning fires at most once per run rather
    /// than once per chromosome.
    high_water_warned: bool,
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
            pending_free: Vec::new(),
            next_fresh: 0,
            slot_refcount: vec![0u8; MAX_ACTIVE_SLOTS as usize],
            pending_mates: AHashMap::new(),
            new_marks: Vec::new(),
            expired_marks: Vec::new(),
            counters: SlotAllocatorCounters::default(),
            high_water_warned: false,
        }
    }

    /// Reset to the initial state. Called at chromosome boundaries
    /// (chains do not span chromosomes per
    /// `ia/specs/pileup_walker.md` §"Chromosome boundaries").
    pub fn reset(&mut self) {
        self.free.clear();
        self.pending_free.clear();
        self.next_fresh = 0;
        for slot in &mut self.slot_refcount {
            *slot = 0;
        }
        self.pending_mates.clear();
        self.new_marks.clear();
        self.expired_marks.clear();
        // counters and high_water_warned are cumulative across the
        // whole run, do not reset
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
        self.maybe_warn_high_water(read.chrom_id, read.alignment_start);
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
    /// zero the slot is parked in `pending_free` (not `free`) and
    /// its id is stamped into `expired_marks`. The slot only
    /// becomes available for reuse at the next
    /// `drain_lifecycle_marks` call, guaranteeing that any
    /// `expired[S]` mark surfaces before any subsequent
    /// `new[S]` for the same id.
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
            self.pending_free.push(slot);
            self.expired_marks.push(slot);
        }
        Ok(())
    }

    /// If `qname` still has a `pending_mates` entry — meaning this
    /// read was admitted as a first mate and its partner never
    /// arrived — release the partner-pre-bumped refcount on the
    /// slot. Called by the active set when a paired read exits, so
    /// the orphan's slot reaches refcount 0 in the same step as the
    /// read's own `release_slot` call (rather than getting stuck at
    /// refcount 1 and never emitting `expired_marks[S]`).
    ///
    /// The call is paired with — and ordered before — the read's own
    /// `release_slot`, so an orphan's two refs collapse to zero
    /// inside a single `expire_passed` step. The expired mark
    /// therefore fires while records are still being drained at this
    /// walker_pos and gets stamped onto the next aged record, just
    /// like a normal slot expiry.
    ///
    /// No-op when `qname` is solo (never registered) or partnered
    /// (the second-mate path consumed the entry on partner arrival).
    /// Counted as a `mate_lookup_eviction` so summary diagnostics
    /// treat the give-up the same as a windowed eviction.
    ///
    /// Regression for finding M1 in `ia/reviews/pileup_2026-05-09.md`.
    pub fn release_pending_partner_ref_if_present(
        &mut self,
        qname: &Arc<str>,
    ) -> Result<(), WalkerError> {
        if let Some(entry) = self.pending_mates.remove(qname) {
            self.counters.mate_lookup_evictions += 1;
            self.release_slot(entry.chain_slot_id)?;
        }
        Ok(())
    }

    /// Drain the pending lifecycle markers into a `(new, expired)`
    /// pair. Both lists are sorted-deduplicated; slots present in
    /// *both* lists at the same drain are emitted in *both* (they
    /// can only be genuinely-transient slots — allocated and
    /// released within one emission window with no allele tagged —
    /// because real reuse is forced into a different drain by the
    /// `pending_free` deferral). The consumer applies `new_chains`
    /// before processing and `expired_chains` after, so a transient
    /// slot has no observable effect on the alive set.
    ///
    /// At end-of-drain, slots parked in `pending_free` are
    /// migrated into `free` and become available for reuse.
    pub fn drain_lifecycle_marks(&mut self) -> (Vec<SlotId>, Vec<SlotId>) {
        sort_dedup(&mut self.new_marks);
        sort_dedup(&mut self.expired_marks);
        let new_out = std::mem::take(&mut self.new_marks);
        let exp_out = std::mem::take(&mut self.expired_marks);

        // Migrate pending_free → free. `free` is kept sorted
        // descending so `pop()` yields the smallest available id.
        while let Some(slot) = self.pending_free.pop() {
            let pos = self
                .free
                .binary_search_by(|s| slot.cmp(s))
                .unwrap_or_else(|p| p);
            self.free.insert(pos, slot);
        }

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

    /// Emit a one-shot soft warning the first time the active-slot
    /// count crosses `HIGH_WATER_WARN_THRESHOLD`. Surfaces a
    /// pathological-coverage region while the run is still alive,
    /// before it potentially trips `MAX_ACTIVE_SLOTS` and dies with
    /// `WalkerError::SlotExhausted`.
    fn maybe_warn_high_water(&mut self, chrom_id: u32, pos: u32) {
        if !self.high_water_warned && self.counters.slot_high_water >= HIGH_WATER_WARN_THRESHOLD {
            self.high_water_warned = true;
            eprintln!(
                "warning: pileup walker reached {}/{} active phase-chain slots at \
                 chrom_id={} pos={}; if usage exceeds {} the run will fail with \
                 SlotExhausted (consider raising --max-active-chain-slots or \
                 pre-filtering this region)",
                self.counters.slot_high_water, MAX_ACTIVE_SLOTS, chrom_id, pos, MAX_ACTIVE_SLOTS,
            );
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
            adaptor_boundary: None,
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
    fn drain_lifecycle_marks_emits_both_for_transient_slot_within_one_drain() {
        // A slot allocated and released between two consecutive
        // drains surfaces in *both* lists of the next drain. The
        // previous "suppress same-id-in-both" rule was wrong: it
        // also masked same-emission reuse (slot 0 released for
        // r1 then reacquired for r2), silently merging two
        // distinct phase chains. With the slot-deferral fix
        // (release goes to `pending_free`, not `free`), reuse
        // can't collide on the same drain — so any same-id-both
        // case is genuinely transient and harmless to surface
        // (the consumer applies new before expired and ends up
        // with no net alive change).
        let mut a = SlotAllocator::new();
        let r = make_read("transient", false, 100);
        let (slot, _) = a.allocate_for_read(&r).unwrap();
        a.release_slot(slot).unwrap();
        let (new, expired) = a.drain_lifecycle_marks();
        assert_eq!(new, vec![slot]);
        assert_eq!(expired, vec![slot]);
    }

    #[test]
    fn same_emission_reuse_does_not_collide_on_slot_id() {
        // r1 acquires slot, releases. Before any drain, r2
        // acquires a slot. The walker MUST NOT hand out r1's old
        // id to r2 within the same emission window — doing so
        // would put the same id in both `new_marks` and
        // `expired_marks` of the next drain, and a downstream
        // consumer could not tell the resulting "chain ends and
        // new chain starts" apart from a transient slot, silently
        // merging the two reads' alleles into a single chain.
        //
        // The slot-deferral fix moves released slots to
        // `pending_free`; only `drain_lifecycle_marks` migrates
        // them to `free`. So r2 here gets a fresh id.
        let mut a = SlotAllocator::new();
        let r1 = make_read("r1", false, 100);
        let (s1, _) = a.allocate_for_read(&r1).unwrap();
        a.release_slot(s1).unwrap();
        let r2 = make_read("r2", false, 100);
        let (s2, _) = a.allocate_for_read(&r2).unwrap();
        assert_ne!(
            s1, s2,
            "same-emission reuse must not collide; got s1={s1}, s2={s2}"
        );

        // After a drain, the freed slot becomes available again.
        let _ = a.drain_lifecycle_marks();
        a.release_slot(s2).unwrap();
        let _ = a.drain_lifecycle_marks();
        let r3 = make_read("r3", false, 100);
        let (s3, _) = a.allocate_for_read(&r3).unwrap();
        assert!(
            s3 == s1 || s3 == s2,
            "post-drain reuse may take any freed slot; got s3={s3}"
        );
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

    #[test]
    fn high_water_warning_fires_once_at_threshold_and_then_stays_set() {
        // Allocate up to one slot below the threshold — the flag
        // must still be clear, since the warning is supposed to
        // fire exactly when the count reaches the threshold.
        let mut a = SlotAllocator::new();
        for i in 0..(HIGH_WATER_WARN_THRESHOLD - 1) {
            let r = make_read(&format!("r{i}"), false, 100);
            a.allocate_for_read(&r).expect("under threshold");
        }
        assert!(
            !a.high_water_warned,
            "must not fire before crossing the threshold"
        );

        // The next allocation lifts active_count to the threshold;
        // the flag flips and the eprintln has fired (side effect
        // not captured here — the flag is the test surface).
        let r = make_read("threshold", false, 100);
        a.allocate_for_read(&r).expect("at threshold");
        assert!(a.high_water_warned, "must fire on reaching threshold");

        // Subsequent allocations must not re-fire. The flag stays
        // set; we'd otherwise spam stderr per allocation in deep
        // regions.
        let r = make_read("post-threshold", false, 100);
        a.allocate_for_read(&r).expect("past threshold");
        assert!(a.high_water_warned, "must remain set (one-shot)");
    }

    #[test]
    fn reset_preserves_high_water_warning_flag() {
        // The warning is one-shot per *run*, not per chromosome —
        // `reset` is called at chromosome boundaries and must not
        // re-arm the flag, otherwise a deep contig would re-warn on
        // every following contig we visit.
        let mut a = SlotAllocator::new();
        for i in 0..HIGH_WATER_WARN_THRESHOLD {
            let r = make_read(&format!("r{i}"), false, 100);
            a.allocate_for_read(&r).unwrap();
        }
        assert!(a.high_water_warned);
        a.reset();
        assert!(
            a.high_water_warned,
            "reset must preserve the one-shot warning flag across chromosomes"
        );
    }

    #[test]
    fn release_pending_partner_paired_with_release_slot_emits_expired_for_orphan() {
        // The exact call sequence the active set's `expire_passed`
        // makes for a paired read exit: partner-ref release first,
        // then the read's own slot release. For an orphan, refcount
        // collapses 2 → 1 → 0 in this single sequence and the
        // expired mark fires on the second call.
        let mut a = SlotAllocator::new();
        let m1 = make_read("orphan", true, 100);
        let (slot, _) = a.allocate_for_read(&m1).unwrap();
        let (new, _) = a.drain_lifecycle_marks();
        assert_eq!(new, vec![slot]);

        a.release_pending_partner_ref_if_present(&m1.qname).unwrap();
        a.release_slot(slot).unwrap();

        let (_, expired) = a.drain_lifecycle_marks();
        assert_eq!(
            expired,
            vec![slot],
            "orphan: 2→1 from partner release, 1→0 from own release, expired emitted",
        );
        assert_eq!(
            a.counters().mate_lookup_evictions,
            1,
            "the abandoned partner counts as one eviction",
        );
    }

    #[test]
    fn release_pending_partner_is_noop_for_solo_read() {
        // Solo reads never enter pending_mates, so the partner-release
        // step is a no-op and the read's own release alone emits
        // expired (refcount 1 → 0).
        let mut a = SlotAllocator::new();
        let solo = make_read("solo", false, 100);
        let (slot, _) = a.allocate_for_read(&solo).unwrap();

        a.release_pending_partner_ref_if_present(&solo.qname).unwrap();
        a.release_slot(slot).unwrap();

        let (_, expired) = a.drain_lifecycle_marks();
        assert_eq!(expired, vec![slot]);
        assert_eq!(
            a.counters().mate_lookup_evictions,
            0,
            "solo path must not be counted as an eviction",
        );
    }

    #[test]
    fn release_pending_partner_is_noop_for_partnered_pair() {
        // Once the second mate arrives, pending_mates is consumed.
        // Either mate's exit calls release_pending_partner_ref_if_present
        // and finds nothing to do — both mates' release_slot calls then
        // collapse the slot 2 → 1 → 0 normally.
        let mut a = SlotAllocator::new();
        let m1 = make_read("p", true, 100);
        let m2 = make_read("p", true, 200);
        let (slot, _) = a.allocate_for_read(&m1).unwrap();
        let _ = a.allocate_for_read(&m2).unwrap();
        assert!(
            a.pending_mates.is_empty(),
            "second-mate path consumed the pending entry",
        );

        // First-to-exit: partner-release is a no-op; release_slot 2→1.
        a.release_pending_partner_ref_if_present(&m1.qname).unwrap();
        a.release_slot(slot).unwrap();
        let (_, expired_mid) = a.drain_lifecycle_marks();
        assert!(expired_mid.is_empty());

        // Second-to-exit: same shape; release_slot 1→0, expired fires.
        a.release_pending_partner_ref_if_present(&m2.qname).unwrap();
        a.release_slot(slot).unwrap();
        let (_, expired_final) = a.drain_lifecycle_marks();
        assert_eq!(expired_final, vec![slot]);
        assert_eq!(
            a.counters().mate_lookup_evictions,
            0,
            "partnered pair must not be counted as an eviction",
        );
    }
}
