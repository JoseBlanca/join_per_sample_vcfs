//! Phase-chain identifier allocator: hands out a fresh, never-reused
//! `u64` to each read (or read pair).
//!
//! See `ia/specs/phase_chain.md`. One read (or one read pair, with
//! mates collapsed via `pending_mates`) gets one identifier. The
//! identifier space is per-`.psp` file, monotonically allocated from
//! zero, and never recycled. Overflow at `u64::MAX + 1` is caught and
//! surfaced as [`WalkerError::ChainIdSpaceExhausted`] — see
//! `ia/feature_implementation_plans/unique_chain_ids.md` §"Overflow
//! guard".
//!
//! With unique ids there is no slot lifecycle to track downstream:
//! the per-record `new_chains` / `expired_chains` markers and the
//! writer-side active-set bookkeeping are gone. The allocator's
//! remaining responsibilities are (a) minting the id, (b) collapsing
//! a read pair onto one id via the `pending_mates` map, and
//! (c) defending against pathological depth via the
//! `max_active_reads` cap that lives on the walker side.
//!
//! The cap on *concurrent active reads* (`WalkerConfig::max_active_slots`)
//! is enforced here as a defensive check: if the number of active reads
//! ever exceeds the cap, [`WalkerError::ActiveReadsExhausted`] surfaces.
//! The cap is independent of, and far smaller than, the `u64`
//! identifier space.

use std::sync::Arc;

use ahash::AHashMap;

use super::PreparedRead;
use super::errors::WalkerError;

/// Phase-chain identifier. `u64` gives ~1.8 × 10¹⁹ values per file —
/// well beyond any realistic read count. Overflow is still caught
/// (see [`SlotAllocator::allocate_for_read`]).
pub type SlotId = u64;

/// Default value for [`WalkerConfig::max_active_slots`]: hard cap on
/// the number of concurrently-active reads. The name is kept for
/// API compatibility; it now bounds *concurrent active reads*, not
/// a slot id namespace (slot ids are unique `u64`s, never recycled).
///
/// [`WalkerConfig::max_active_slots`]: super::WalkerConfig::max_active_slots
pub const DEFAULT_MAX_ACTIVE_SLOTS: u32 = 4096;

/// Hard cap on the number of entries the `pending_mates` map will
/// hold at any time. Exceeding it surfaces as
/// [`WalkerError::PendingMatesExhausted`].
///
/// The map's natural bound under `evict_stale_pending` is
/// `(orphans per bp) × mate_lookup_window`. At realistic 30× coverage
/// with ~5 % orphan rate and a 10 kb window, peak is ~100 entries.
/// Even on pessimistic 300× coverage with 50 % orphans the peak sits
/// around 10 k. This cap (10 k) catches truly malformed inputs (e.g.
/// a BAM where every paired read carries the FirstOfPair flag and no
/// second mate ever arrives) without firing on real data. Bumping
/// the constant if `mate_lookup_window` is configured much higher
/// for long-insert libraries is a one-line change.
pub const MAX_PENDING_MATES: usize = 10_000;

/// Fraction of the active-read cap at which the allocator emits a
/// one-shot soft warning. Lets the user spot a pathological-coverage
/// region before the run trips the cap.
fn high_water_warn_threshold(cap: u32) -> u32 {
    cap.saturating_mul(3) / 4
}

/// State the allocator holds for a first mate whose partner has not
/// yet been admitted: which chain id was issued, the first mate's
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
    /// Next never-used chain id. Incremented on every fresh
    /// allocation via [`checked_add`](u64::checked_add) so silent
    /// wrap to zero is impossible — overflow surfaces as
    /// [`WalkerError::ChainIdSpaceExhausted`].
    next_id: u64,
    /// Count of currently-active reads. Bumped on every successful
    /// fresh allocation, decremented in [`Self::note_read_exit`].
    /// Enforced against `max_active_reads` to surface
    /// `ActiveReadsExhausted` before pathological-depth regions
    /// blow up memory.
    active_count: u32,
    /// First mates whose partner has not yet arrived. Sized for
    /// genuinely-pending pairs only — solo reads
    /// (`MateRole::Solo`) never enter this map.
    pending_mates: AHashMap<Arc<str>, PendingMate>,
    /// Bookkeeping for the run summary.
    counters: SlotAllocatorCounters,
    /// Set the first time the active-read count reaches
    /// `high_water_warn_threshold(max_active_reads)`. Idempotent
    /// within a run — and because `reset` (called at chromosome
    /// boundaries) deliberately preserves it, the warning fires
    /// at most once per run rather than once per chromosome.
    high_water_warned: bool,
    /// Hard cap on concurrent active reads. The name "max_active_slots"
    /// is preserved on the public config struct for backwards-
    /// compatible callers; here we store the same value.
    max_active_reads: u32,
    /// Per-instance pending-mate lookup window, mirrors
    /// `WalkerConfig::mate_lookup_window`.
    mate_lookup_window: u32,
    /// Per-instance pending-mates cap. Initialised from
    /// [`MAX_PENDING_MATES`]; tests can lower it via the
    /// `#[cfg(test)]` constructor so the cap path is reachable
    /// without inserting 10 000 qnames.
    pending_mates_cap: usize,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct SlotAllocatorCounters {
    /// Number of fresh chain ids minted (i.e., `next_id`
    /// increments). Second-mate allocations do not count.
    pub slot_allocations: u64,
    /// Peak observed value of `active_count` during the run.
    pub slot_high_water: u32,
    /// Number of first mates whose partner never arrived within
    /// `mate_lookup_window` and were evicted from `pending_mates`.
    pub mate_lookup_evictions: u64,
}

impl SlotAllocator {
    /// Construct with the default cap and lookup window. Used by
    /// tests; production code calls [`SlotAllocator::with_caps`]
    /// from `walker::run`.
    #[cfg(test)]
    pub fn new() -> Self {
        Self::with_caps(DEFAULT_MAX_ACTIVE_SLOTS, super::DEFAULT_MATE_LOOKUP_WINDOW)
    }

    /// Construct with explicit caps. `max_active_reads` is the
    /// hard cap on concurrent reads; `mate_lookup_window` is the
    /// pending-mate eviction distance.
    pub fn with_caps(max_active_reads: u32, mate_lookup_window: u32) -> Self {
        Self {
            next_id: 0,
            active_count: 0,
            pending_mates: AHashMap::new(),
            counters: SlotAllocatorCounters::default(),
            high_water_warned: false,
            max_active_reads,
            mate_lookup_window,
            pending_mates_cap: MAX_PENDING_MATES,
        }
    }

    /// Test-only constructor that seeds [`Self::next_id`] near the
    /// `u64` ceiling so the overflow-guard tests can exercise the
    /// `ChainIdSpaceExhausted` path without literally allocating
    /// 2⁶⁴ chains.
    #[cfg(test)]
    pub(crate) fn with_next_id_for_testing(start: u64) -> Self {
        Self {
            next_id: start,
            ..Self::new()
        }
    }

    /// Test-only constructor that lowers the pending-mates cap to
    /// `cap` so the `PendingMatesExhausted` path is reachable
    /// without inserting `MAX_PENDING_MATES` qnames.
    #[cfg(test)]
    pub(crate) fn with_pending_mates_cap_for_testing(cap: usize) -> Self {
        Self {
            pending_mates_cap: cap,
            ..Self::new()
        }
    }

    /// Reset to the initial state. Called at chromosome boundaries
    /// (chains do not span chromosomes per
    /// `ia/specs/calling_pipeline_architecture.md` §"Phase chain
    /// identifiers"). `next_id` is *not* reset — chain ids stay
    /// unique across the whole file.
    pub fn reset(&mut self) {
        self.active_count = 0;
        self.pending_mates.clear();
        // counters, high_water_warned, max_active_reads,
        // mate_lookup_window, and crucially next_id are preserved
        // across chromosome resets — they are file-scoped.
    }

    /// Allocate a chain id for an entering read.
    ///
    /// - First mate of a pair (`mate_role.is_paired()` true, qname
    ///   not in `pending_mates`): mint a fresh id and register the
    ///   qname. `active_count` is incremented once.
    /// - Second mate (qname found in `pending_mates`): reuse the
    ///   first mate's id; the pair shares one chain. `active_count`
    ///   is bumped a second time because both mates are now
    ///   simultaneously in the active set.
    /// - Solo read (`MateRole::Solo`): mint a fresh id; `active_count`
    ///   is incremented once.
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

        if read.mate_role.is_paired()
            && let Some(pending) = self.pending_mates.remove(&read.qname)
        {
            // Second mate of a known pair — reuse the id. We do
            // bump `active_count` here because both mates now sit
            // simultaneously in the active set; the cap protects
            // against pathological per-position depth, not against
            // distinct chain ids.
            self.bump_active_count(read.chrom_id, read.alignment_start)?;
            self.maybe_warn_high_water(read.chrom_id, read.alignment_start);
            return Ok((pending.chain_slot_id, Some(pending.first_mate_read_id)));
        }

        // Fresh allocation: bump active_count first (so the cap
        // surfaces before we mint a new id we'd then have to roll
        // back), then mint.
        self.bump_active_count(read.chrom_id, read.alignment_start)?;
        let slot = self.next_id;
        // Bump next_id with overflow check. `2^64` chain ids per
        // .psp file is astronomically beyond any realistic
        // workload, but silent wrap-around would silently merge
        // two distinct molecules into one chain id with no way to
        // detect the collision after the fact — see
        // `ia/feature_implementation_plans/unique_chain_ids.md`
        // §"Overflow guard".
        self.next_id = self
            .next_id
            .checked_add(1)
            .ok_or(WalkerError::ChainIdSpaceExhausted {
                chrom_id: read.chrom_id,
                pos: read.alignment_start,
            })?;
        self.counters.slot_allocations += 1;
        self.maybe_warn_high_water(read.chrom_id, read.alignment_start);

        if read.mate_role.is_paired() {
            // Defensive cap on the pending-mates map size. The
            // `evict_stale_pending` call at the top of every
            // `allocate_for_read` already bounds the map via the
            // mate-lookup window, so under realistic workloads
            // this check never fires. It catches truly malformed
            // inputs (every paired read flagged FirstOfPair with
            // no SecondOfPair ever arriving) before they consume
            // meaningful memory.
            if self.pending_mates.len() >= self.pending_mates_cap {
                return Err(WalkerError::PendingMatesExhausted {
                    cap: self.pending_mates_cap,
                    chrom_id: read.chrom_id,
                    pos: read.alignment_start,
                });
            }
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

    /// Note that a read has exited the active set: decrement the
    /// active-read count. With unique chain ids the act of expiring
    /// a read no longer touches a slot id namespace — the id stays
    /// unique for the file — so this is just refcount-style
    /// bookkeeping for the active-read cap.
    ///
    /// `chrom_id` and `pos` are only consulted on the error path:
    /// they populate `WalkerError::Internal`'s context fields if the
    /// active-read count is asked to go below zero.
    pub fn note_read_exit(&mut self, chrom_id: u32, pos: u32) -> Result<(), WalkerError> {
        if self.active_count == 0 {
            return Err(WalkerError::Internal {
                detail: "note_read_exit called with active_count already zero".to_string(),
                qname: String::new(),
                chrom_id,
                pos,
            });
        }
        self.active_count -= 1;
        Ok(())
    }

    pub fn counters(&self) -> SlotAllocatorCounters {
        self.counters
    }

    /// Bump the active-read count and surface
    /// [`WalkerError::ActiveReadsExhausted`] if the new value would
    /// exceed the cap. Track the high-water mark.
    fn bump_active_count(&mut self, chrom_id: u32, pos: u32) -> Result<(), WalkerError> {
        if self.active_count >= self.max_active_reads {
            return Err(WalkerError::ActiveReadsExhausted {
                cap: self.max_active_reads,
                chrom_id,
                pos,
            });
        }
        self.active_count += 1;
        if self.active_count > self.counters.slot_high_water {
            self.counters.slot_high_water = self.active_count;
        }
        Ok(())
    }

    /// Emit a one-shot soft warning the first time `active_count`
    /// crosses 75 % of the configured cap. Surfaces a pathological-
    /// coverage region while the run is still alive, before it
    /// potentially trips the cap and dies with
    /// `WalkerError::ActiveReadsExhausted`.
    fn maybe_warn_high_water(&mut self, chrom_id: u32, pos: u32) {
        let threshold = high_water_warn_threshold(self.max_active_reads);
        if !self.high_water_warned && self.counters.slot_high_water >= threshold {
            self.high_water_warned = true;
            eprintln!(
                "warning: pileup walker reached {}/{} active reads at \
                 chrom_id={} pos={}; if usage exceeds {} the run will fail with \
                 ActiveReadsExhausted (raise WalkerConfig::max_active_slots or \
                 pre-filter this region)",
                self.counters.slot_high_water,
                self.max_active_reads,
                chrom_id,
                pos,
                self.max_active_reads,
            );
        }
    }

    fn evict_stale_pending(&mut self, walker_pos: u32) {
        // `AHashMap::retain` lets us mutate the map in a single
        // pass without per-entry `Arc::clone`. The eviction only
        // clears the `pending_mates` entry; `active_count` stays
        // attached to the first mate's still-live active-read entry
        // and is released when that read exits.
        let counters = &mut self.counters;
        let mate_lookup_window = self.mate_lookup_window;
        self.pending_mates.retain(|_qname, entry| {
            // `saturating_add` guards `seen_at` near `u32::MAX` on
            // multi-Gbp chromosomes.
            if walker_pos > entry.seen_at.saturating_add(mate_lookup_window) {
                counters.mate_lookup_evictions += 1;
                false // drop
            } else {
                true // keep
            }
        });
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::super::MateRole;
    use super::*;

    fn make_read(qname: &str, mate_role: MateRole, alignment_start: u32) -> PreparedRead {
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
            mate_role,
            adaptor_boundary: None,
        }
    }

    #[test]
    fn solo_read_allocates_a_fresh_id_starting_at_zero() {
        let mut a = SlotAllocator::new();
        let read = make_read("solo", MateRole::Solo, 100);
        let (slot, mate_id) = a.allocate_for_read(&read).expect("solo allocates");
        assert_eq!(slot, 0);
        assert_eq!(mate_id, None);
        assert!(a.pending_mates.is_empty());
    }

    #[test]
    fn ids_are_unique_and_monotonically_increasing_across_solo_reads() {
        let mut a = SlotAllocator::new();
        let r0 = make_read("a", MateRole::Solo, 100);
        let r1 = make_read("b", MateRole::Solo, 100);
        let r2 = make_read("c", MateRole::Solo, 100);
        let (s0, _) = a.allocate_for_read(&r0).unwrap();
        let (s1, _) = a.allocate_for_read(&r1).unwrap();
        let (s2, _) = a.allocate_for_read(&r2).unwrap();
        assert_eq!((s0, s1, s2), (0, 1, 2));
    }

    #[test]
    fn first_mate_registers_then_second_mate_reuses_id() {
        let mut a = SlotAllocator::new();
        let m1 = make_read("pair1", MateRole::FirstOfPair, 100);
        let (slot1, mate_id1) = a.allocate_for_read(&m1).expect("first mate");
        a.register_first_mate_read_id(&m1.qname, 42);

        let m2 = make_read("pair1", MateRole::FirstOfPair, 200);
        let (slot2, mate_id2) = a.allocate_for_read(&m2).expect("second mate");
        assert_eq!(slot1, slot2, "mates must share a chain id");
        assert_eq!(mate_id1, None, "first mate has no partner registered yet");
        assert_eq!(mate_id2, Some(42), "second mate sees first mate's read_id");
        assert!(
            a.pending_mates.is_empty(),
            "pending entry consumed when second mate arrives"
        );
    }

    #[test]
    fn released_ids_are_not_recycled() {
        // Same-id reuse is the entire bug class this design exists to
        // remove: a freed chain id must *not* surface again, even
        // after the slot allocator has been told the read exited.
        let mut a = SlotAllocator::new();
        let r1 = make_read("r1", MateRole::Solo, 100);
        let (s1, _) = a.allocate_for_read(&r1).unwrap();
        a.note_read_exit(0, 100).unwrap();
        let r2 = make_read("r2", MateRole::Solo, 100);
        let (s2, _) = a.allocate_for_read(&r2).unwrap();
        assert_ne!(s1, s2, "released ids must never be recycled");
        assert_eq!((s1, s2), (0, 1));
    }

    #[test]
    fn active_reads_cap_returns_hard_error() {
        let mut a = SlotAllocator::new();
        for i in 0..DEFAULT_MAX_ACTIVE_SLOTS {
            let r = make_read(&format!("r{i}"), MateRole::Solo, 100);
            a.allocate_for_read(&r).expect("under cap");
        }
        let r = make_read("overflow", MateRole::Solo, 100);
        let err = a.allocate_for_read(&r).expect_err("must hit cap");
        assert!(
            matches!(err, WalkerError::ActiveReadsExhausted { .. }),
            "got {err:?}"
        );
    }

    #[test]
    fn stale_pending_mates_are_evicted_after_window() {
        let mut a = SlotAllocator::new();
        let m1 = make_read("orphan", MateRole::FirstOfPair, 100);
        let _ = a.allocate_for_read(&m1).unwrap();
        assert_eq!(a.pending_mates.len(), 1);

        let r2 = make_read(
            "later",
            MateRole::Solo,
            100 + super::super::DEFAULT_MATE_LOOKUP_WINDOW + 1,
        );
        a.allocate_for_read(&r2).unwrap();
        assert!(
            a.pending_mates.is_empty(),
            "orphan first mate should have been evicted"
        );
        assert_eq!(a.counters().mate_lookup_evictions, 1);
    }

    #[test]
    fn reset_clears_active_count_but_preserves_next_id_and_counters() {
        // `reset` is called at chromosome boundaries. With unique ids
        // the *id counter* must stay monotonic across the boundary —
        // otherwise the same id would surface on two chromosomes.
        let mut a = SlotAllocator::new();
        let r = make_read("x", MateRole::Solo, 100);
        a.allocate_for_read(&r).unwrap();
        let allocs_before = a.counters().slot_allocations;
        let next_id_before = a.next_id;
        a.reset();
        assert_eq!(a.active_count, 0);
        assert!(a.pending_mates.is_empty());
        assert_eq!(
            a.next_id, next_id_before,
            "next_id must persist across reset",
        );
        assert_eq!(a.counters().slot_allocations, allocs_before);
    }

    #[test]
    fn high_water_warning_fires_once_at_threshold_and_then_stays_set() {
        let mut a = SlotAllocator::new();
        let threshold = high_water_warn_threshold(DEFAULT_MAX_ACTIVE_SLOTS);
        for i in 0..(threshold - 1) {
            let r = make_read(&format!("r{i}"), MateRole::Solo, 100);
            a.allocate_for_read(&r).expect("under threshold");
        }
        assert!(
            !a.high_water_warned,
            "must not fire before crossing the threshold"
        );

        let r = make_read("threshold", MateRole::Solo, 100);
        a.allocate_for_read(&r).expect("at threshold");
        assert!(a.high_water_warned, "must fire on reaching threshold");

        let r = make_read("post-threshold", MateRole::Solo, 100);
        a.allocate_for_read(&r).expect("past threshold");
        assert!(a.high_water_warned, "must remain set (one-shot)");
    }

    #[test]
    fn reset_preserves_high_water_warning_flag() {
        let mut a = SlotAllocator::new();
        let threshold = high_water_warn_threshold(DEFAULT_MAX_ACTIVE_SLOTS);
        for i in 0..threshold {
            let r = make_read(&format!("r{i}"), MateRole::Solo, 100);
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
    fn note_read_exit_on_empty_active_set_errors_with_locus_context() {
        let mut a = SlotAllocator::new();
        let err = a
            .note_read_exit(/* chrom_id */ 7, /* pos */ 12345)
            .expect_err("note_read_exit with empty active set must error");
        match err {
            WalkerError::Internal {
                detail,
                chrom_id,
                pos,
                ..
            } => {
                assert_eq!(chrom_id, 7);
                assert_eq!(pos, 12345);
                assert!(
                    detail.contains("active_count"),
                    "detail must mention the active_count invariant; got: {detail}",
                );
            }
            other => panic!("expected WalkerError::Internal, got {other:?}"),
        }
    }

    #[test]
    fn pending_mate_entry_survives_first_mate_release_for_window_matching() {
        // Pair-tracking is governed by `mate_lookup_window`, not by
        // active-set residence. After the first mate exits the
        // active set, its `pending_mates` entry must survive so a
        // later second-mate arrival within the window can still
        // match.
        let mut a = SlotAllocator::new();
        let m1 = make_read("pair", MateRole::FirstOfPair, 100);
        let (slot1, _) = a.allocate_for_read(&m1).unwrap();
        // First mate exits the active set.
        a.note_read_exit(0, 200).unwrap();
        assert_eq!(a.active_count, 0);
        assert_eq!(
            a.pending_mates.len(),
            1,
            "pending entry must outlive the first mate's release"
        );
        // Second mate arrives later (within the window); reuses
        // the first mate's chain id.
        let m2 = make_read("pair", MateRole::SecondOfPair, 200);
        let (slot2, partner) = a.allocate_for_read(&m2).unwrap();
        assert_eq!(slot1, slot2, "mates share the same chain id");
        assert_eq!(partner, Some(0), "second-mate path saw the first mate");
        assert!(
            a.pending_mates.is_empty(),
            "second-mate consumption clears the pending entry"
        );
    }

    // --- Overflow guard ------------------------------------------------

    #[test]
    fn allocator_errors_when_chain_id_counter_would_overflow() {
        // Seed `next_id` to `u64::MAX` via the test constructor. The
        // very next fresh allocation must error rather than wrap to
        // zero (which would silently produce colliding chain ids).
        let mut a = SlotAllocator::with_next_id_for_testing(u64::MAX);
        let r = make_read("at-max", MateRole::Solo, 100);
        let err = a
            .allocate_for_read(&r)
            .expect_err("must surface ChainIdSpaceExhausted at u64::MAX");
        match err {
            WalkerError::ChainIdSpaceExhausted { chrom_id, pos } => {
                assert_eq!(chrom_id, 0);
                assert_eq!(pos, 100);
            }
            other => panic!("expected ChainIdSpaceExhausted, got {other:?}"),
        }
    }

    #[test]
    fn allocator_succeeds_at_chain_id_counter_one_below_max() {
        // The off-by-one boundary: with `next_id = u64::MAX - 1` the
        // walker must successfully admit one more read (returning
        // id `u64::MAX - 1`) and only then refuse on the next
        // allocation. Pins that the `checked_add` is at the right
        // place and the test couldn't accidentally relax to
        // `wrapping_add` without failing this test.
        let mut a = SlotAllocator::with_next_id_for_testing(u64::MAX - 1);
        let r = make_read("at-max-minus-one", MateRole::Solo, 100);
        let (slot, _) = a
            .allocate_for_read(&r)
            .expect("must succeed at u64::MAX - 1");
        assert_eq!(slot, u64::MAX - 1);

        // The next fresh allocation would mint `u64::MAX` and try to
        // bump `next_id` to `u64::MAX + 1`, which overflows — error.
        let r2 = make_read("would-be-max", MateRole::Solo, 100);
        // Active_count is now 1; release first so the active-cap check
        // doesn't pre-empt the overflow check.
        a.note_read_exit(0, 100).unwrap();
        let err = a.allocate_for_read(&r2).expect_err("must surface overflow");
        assert!(
            matches!(err, WalkerError::ChainIdSpaceExhausted { .. }),
            "got {err:?}"
        );
    }

    // --- Pending-mates cap --------------------------------------------

    #[test]
    fn allocator_errors_when_pending_mates_cap_exceeded() {
        // Lower the pending-mates cap to a tiny number via the
        // test-only constructor so the cap path is reachable without
        // 10 000 qnames. With cap=2, the third orphan first-mate
        // admission must fail with PendingMatesExhausted.
        let mut a = SlotAllocator::with_pending_mates_cap_for_testing(2);
        let r0 = make_read("a", MateRole::FirstOfPair, 100);
        let r1 = make_read("b", MateRole::FirstOfPair, 101);
        a.allocate_for_read(&r0).unwrap();
        a.allocate_for_read(&r1).unwrap();
        assert_eq!(a.pending_mates.len(), 2);
        let r2 = make_read("c", MateRole::FirstOfPair, 102);
        let err = a
            .allocate_for_read(&r2)
            .expect_err("must surface PendingMatesExhausted at the cap");
        match err {
            WalkerError::PendingMatesExhausted { cap, chrom_id, pos } => {
                assert_eq!(cap, 2);
                assert_eq!(chrom_id, 0);
                assert_eq!(pos, 102);
            }
            other => panic!("expected PendingMatesExhausted, got {other:?}"),
        }
    }

    #[test]
    fn pending_mates_cap_does_not_block_second_mate_arrival() {
        // The cap check is on `pending_mates.insert`. A second-mate
        // arrival takes the existing-entry path (consumes the entry,
        // does not insert), so the cap does not fire there — even
        // with `pending_mates.len() == cap`, a second mate can still
        // resolve cleanly.
        let mut a = SlotAllocator::with_pending_mates_cap_for_testing(1);
        let r0 = make_read("pair", MateRole::FirstOfPair, 100);
        a.allocate_for_read(&r0).unwrap();
        assert_eq!(a.pending_mates.len(), 1);
        let r1 = make_read("pair", MateRole::SecondOfPair, 200);
        a.allocate_for_read(&r1)
            .expect("second-mate path must not hit the cap");
        assert!(a.pending_mates.is_empty());
    }
}
