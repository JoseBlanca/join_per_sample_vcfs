//! The walker loop. Single-threaded; drives all the building
//! blocks (active set, slot allocator, open-record table) through
//! the closure rule and yields finalised `PileupRecord`s through
//! a pull-shaped `Iterator`.
//!
//! The walker is a state machine: each call to `Iterator::next()`
//! advances it until at least one record is ready (or end-of-input
//! is reached), then yields one record at a time. A single walker
//! tick may emit 0, 1, or many records; the iterator buffers them
//! in a small `VecDeque` and drains across successive `next()`
//! calls. Tick-internal order is preserved so the lifecycle-mark
//! attachment rule — marks stamped onto the first record of the
//! tick's batch — survives the change.

use std::collections::VecDeque;
use std::iter::Peekable;

use ahash::AHashMap;

use super::active_read_set::ActiveReads;
use super::decompose::ReadEvent;
use super::errors::WalkerError;
use super::open_record::{
    OpenPileupRecord, OpenPileupRecordTable, ReadContribution, process_position,
};
use super::slot_allocator::{SlotAllocator, SlotAllocatorCounters, SlotId};
use super::{PileupRecord, PreparedRead, ReadLengthError, RefSeqFetcher, WalkerConfig};

/// Construct a [`PileupWalker`] over a coordinate-sorted stream of
/// prepared reads. The walker is an `Iterator<Item = Result<PileupRecord,
/// WalkerError>>`; callers drive it by repeatedly calling `next()`
/// (or by collecting / for-looping). After iteration ends, the
/// run's cumulative counters are available via
/// [`PileupWalker::summary`].
///
/// Coordinate-order invariant: every read pulled from `reads` must
/// have `(chrom_id, alignment_start)` non-decreasing relative to
/// the previous one. A regression is a hard error — stale or
/// malformed input shouldn't pass silently.
pub fn run<R, F>(reads: R, ref_fetcher: F, config: &WalkerConfig) -> PileupWalker<R::IntoIter, F>
where
    R: IntoIterator<Item = PreparedRead>,
    F: RefSeqFetcher,
{
    PileupWalker::new(reads.into_iter(), ref_fetcher, config)
}

/// Pull-shaped walker over a coordinate-sorted stream of prepared
/// reads. See [`run`] for the convenience constructor.
pub struct PileupWalker<I, F>
where
    I: Iterator<Item = PreparedRead>,
    F: RefSeqFetcher,
{
    reads: Peekable<I>,
    ref_fetcher: F,
    state: WalkerState,
    /// Records produced by walker ticks but not yet consumed by
    /// `Iterator::next`. A single tick may emit 0–many records
    /// (e.g. a wide deletion at an earlier anchor unblocks several
    /// narrower records simultaneously); they're appended here in
    /// emission order and drained via `pop_front`.
    pending: VecDeque<PileupRecord>,
    /// `true` once end-of-input has been flushed *or* a `next()`
    /// call has returned an error. Both terminal states stop the
    /// iterator from doing further work.
    done: bool,
}

impl<I, F> PileupWalker<I, F>
where
    I: Iterator<Item = PreparedRead>,
    F: RefSeqFetcher,
{
    pub fn new(reads: I, ref_fetcher: F, config: &WalkerConfig) -> Self {
        let mut reads = reads.peekable();
        let mut state = WalkerState::new(*config);
        // Initial chromosome anchor: the first peeked read sets
        // `chrom_id`, `walker_pos = 1`, and
        // `last_admitted_chrom_id`. Subsequent re-anchors happen
        // inside the chromosome-boundary block of `fill_pending`,
        // right after `flush_chromosome_into`.
        if let Some(first) = reads.peek() {
            state.enter_chrom(first.chrom_id);
        }
        Self {
            reads,
            ref_fetcher,
            state,
            pending: VecDeque::new(),
            done: false,
        }
    }

    /// Cumulative counters for the run so far. Safe to call
    /// mid-stream; the final summary is the value observed after
    /// `next()` has returned `None`.
    pub fn summary(&self) -> RunSummary {
        self.state.summary()
    }

    /// Drive the walker until at least one record is ready in
    /// `pending`, or until end-of-input. End-of-input also flushes
    /// the remaining chromosome and sets `done = true`.
    fn fill_pending(&mut self) -> Result<(), WalkerError> {
        loop {
            // Terminal condition: no more reads to pull and the
            // active set is empty. Stopping at "reads empty" alone
            // would leak open records whose anchors sit ahead of
            // the walker but inside the active set's coverage.
            if self.reads.peek().is_none() && self.state.active_reads.is_empty() {
                self.state.flush_chromosome_into(&mut self.pending)?;
                self.done = true;
                return Ok(());
            }

            // Chromosome boundary: the next peeked read sits on a
            // new chromosome. Finalise everything still in flight
            // from the previous chromosome and re-anchor.
            //
            // The forward-direction check has to run *before* the
            // flush so a backward chromosome change errors out
            // without first emitting the previous chromosome's
            // records into `pending`.
            //
            // We pull the chrom_id and qname out of the peek
            // borrow into locals so the subsequent calls into
            // `self.state` aren't blocked by the peek borrow.
            let chrom_transition: Option<u32> = {
                let peeked = self.reads.peek();
                match (peeked, self.state.last_admitted_chrom_id) {
                    (Some(p), Some(prev)) if prev != p.chrom_id => {
                        if p.chrom_id < prev {
                            let prev_pos = self
                                .state
                                .last_admitted_locus
                                .map(|l| l.pos)
                                .unwrap_or(self.state.walker_pos);
                            return Err(WalkerError::OutOfOrder {
                                qname: p.qname.to_string(),
                                prev_chrom_id: prev,
                                prev_pos,
                                chrom_id: p.chrom_id,
                                pos: p.alignment_start,
                            });
                        }
                        Some(p.chrom_id)
                    }
                    _ => None,
                }
            };
            if let Some(new_chrom) = chrom_transition {
                self.state.flush_chromosome_into(&mut self.pending)?;
                self.state.enter_chrom(new_chrom);
            }

            // Pull every read with alignment_start ≤ walker_pos
            // (only on the current chromosome; reads on later
            // chromosomes wait for the chromosome flush above).
            while let Some(peeked_read) = self.reads.peek() {
                if peeked_read.chrom_id != self.state.chrom_id {
                    break;
                }
                if peeked_read.alignment_start > self.state.walker_pos {
                    break;
                }
                // PANIC-FREE: `peek()` returned Some on the loop
                // condition above, and `self.reads` has not been
                // advanced between then and here.
                let r = self.reads.next().expect("peek matched");
                self.state.admit_read(r)?;
            }

            // Process events at walker_pos, expire passed reads,
            // and close aged records. Order matters: expire BEFORE
            // close so any slot lifecycle marks generated by expiry
            // land on records emitted at this same walker step
            // rather than getting deferred to the next one.
            self.state.process_position(&self.ref_fetcher)?;
            self.state.expire_passed_reads()?;
            self.state.close_aged_records_into(&mut self.pending);

            // Advance walker_pos to the next interesting position:
            // walker_pos+1 if any active read still has events at
            // or beyond it, otherwise jump to the next read's
            // alignment_start (skip uncovered span).
            self.state.advance(self.reads.peek())?;

            if !self.pending.is_empty() {
                return Ok(());
            }
        }
    }
}

impl<I, F> Iterator for PileupWalker<I, F>
where
    I: Iterator<Item = PreparedRead>,
    F: RefSeqFetcher,
{
    type Item = Result<PileupRecord, WalkerError>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(record) = self.pending.pop_front() {
            return Some(Ok(record));
        }
        if self.done {
            return None;
        }
        match self.fill_pending() {
            Ok(()) => self.pending.pop_front().map(Ok),
            Err(e) => {
                // Terminal-on-first-error: stop yielding after this.
                // The previous push-based shape stopped emission as
                // soon as `run` returned `Err`; preserve that.
                self.done = true;
                Some(Err(e))
            }
        }
    }
}

/// Cumulative counters reported back by `run` so callers can log a
/// per-sample summary.
#[derive(Debug, Default, Clone, Copy)]
#[non_exhaustive]
pub struct RunSummary {
    pub reads_admitted: u64,
    pub records_emitted: u64,
    pub record_widen_events: u64,
    pub mate_overlap_positions: u64,
    pub slot_allocations: u64,
    pub slot_high_water: u32,
    pub mate_lookup_evictions: u64,
    /// Number of columns where the contributor list was truncated
    /// because depth exceeded the applicable per-column cap (see
    /// `WalkerConfig::max_snp_column_depth` /
    /// `max_indel_column_depth`). A non-zero value flags
    /// pathologically deep regions; QC pipelines may want to look
    /// at those samples / regions specifically.
    pub column_depth_truncations: u64,
}

impl RunSummary {
    /// Fold the slot allocator's counters into this summary at
    /// run-end. The walker tracks reads/records/widens/overlap
    /// itself; the allocator owns slot bookkeeping.
    fn merge_slot_counters(mut self, slot: SlotAllocatorCounters) -> Self {
        self.slot_allocations = slot.slot_allocations;
        self.slot_high_water = slot.slot_high_water;
        self.mate_lookup_evictions = slot.mate_lookup_evictions;
        self
    }
}

/// A genomic locus: a position on a specific chromosome. The
/// `Ord` derive gives lexicographic ordering — chromosomes
/// compared first, position within a chromosome second — which is
/// exactly the walker's coordinate-sort invariant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Locus {
    chrom_id: u32,
    /// 1-based reference position.
    pos: u32,
}

struct WalkerState {
    chrom_id: u32,
    walker_pos: u32,
    /// `None` until the first read is admitted. Tracks the
    /// chromosome the walker has been processing so the
    /// flush-on-chromosome-change logic knows whether to flush
    /// before re-anchoring.
    last_admitted_chrom_id: Option<u32>,
    /// Last admitted locus for the coordinate-order invariant.
    last_admitted_locus: Option<Locus>,
    active_reads: ActiveReads,
    slots: SlotAllocator,
    open_records: OpenPileupRecordTable,
    summary: RunSummary,
    config: WalkerConfig,
    /// Reusable per-step buffer for the contributors list. Hoisted
    /// here so the per-walker-step `Vec<ReadContribution>` is
    /// allocated once and reused via `clear()` between steps. L6
    /// in `ia/reviews/perf_pileup_2026-05-10.md`.
    contributors_buf: Vec<ReadContribution>,
    /// Reusable per-step buffer for `close_aged_records`'s drained
    /// records. Paired with `OpenPileupRecordTable::closing_keys_buf`;
    /// together they remove the two per-walker-step `Vec` allocations
    /// `drain_aged` paid in round-1. H2 in
    /// `ia/reviews/perf_pileup_2026-05-12.md`.
    drained_buf: Vec<OpenPileupRecord>,
}

impl WalkerState {
    fn new(config: WalkerConfig) -> Self {
        Self {
            chrom_id: 0,
            walker_pos: 1,
            last_admitted_chrom_id: None,
            last_admitted_locus: None,
            active_reads: ActiveReads::new(),
            slots: SlotAllocator::with_caps(config.max_active_slots, config.mate_lookup_window),
            open_records: OpenPileupRecordTable::with_cap(config.max_record_span),
            summary: RunSummary::default(),
            config,
            contributors_buf: Vec::new(),
            drained_buf: Vec::new(),
        }
    }

    /// Anchor the walker to a chromosome. Called twice in the
    /// run lifecycle: once before the loop with the first peeked
    /// read's chrom, and again inside the boundary block right
    /// after `flush_chromosome` when the next peeked read sits on
    /// a new chrom. Walker_pos resets per chromosome — the walker
    /// only emits within a chromosome, so position numbering
    /// restarts from 1.
    ///
    /// `last_admitted_locus` is deliberately preserved across
    /// chromosome boundaries (Mi14 in
    /// `ia/reviews/pileup_2026-05-11.md`): the per-read tuple
    /// comparison in `admit_read` correctly admits a forward
    /// chrom change (`(new_chrom, _) > (old_chrom, _)` holds
    /// whenever `new_chrom > old_chrom`), and keeping the locus
    /// sticky lets the outer chrom regression's error message
    /// report the *actual* last admitted (chrom, pos) instead of
    /// falling back to a misleading `walker_pos`.
    fn enter_chrom(&mut self, chrom_id: u32) {
        self.chrom_id = chrom_id;
        self.walker_pos = 1;
        self.last_admitted_chrom_id = Some(chrom_id);
    }

    fn admit_read(&mut self, read: PreparedRead) -> Result<(), WalkerError> {
        // Order-invariant check.
        let read_locus = Locus {
            chrom_id: read.chrom_id,
            pos: read.alignment_start,
        };
        if let Some(prev) = self.last_admitted_locus
            && read_locus < prev
        {
            return Err(WalkerError::OutOfOrder {
                qname: read.qname.to_string(),
                prev_chrom_id: prev.chrom_id,
                prev_pos: prev.pos,
                chrom_id: read_locus.chrom_id,
                pos: read_locus.pos,
            });
        }
        // Zero-ref-span check.
        if read.alignment_end < read.alignment_start {
            return Err(WalkerError::ZeroRefSpan {
                qname: read.qname.to_string(),
                chrom_id: read.chrom_id,
                pos: read.alignment_start,
            });
        }
        // Length invariants the cursor relies on. See
        // `PreparedRead::length` for the rationale.
        read.length()
            .map_err(|e| malformed_read_from_length_err(e, &read))?;
        self.last_admitted_locus = Some(read_locus);
        self.active_reads.admit(read, &mut self.slots)?;
        self.summary.reads_admitted += 1;
        Ok(())
    }

    fn process_position<F: RefSeqFetcher>(&mut self, ref_fetcher: &F) -> Result<(), WalkerError> {
        // Step 1: query each active read's cursor for events
        // anchored at walker_pos. Reads with no event here are
        // silent (deletion interior or N-skip), so they are not
        // added as contributors at all.
        let walker_pos = self.walker_pos;
        // Hoisted buffer; cleared per step. L6.
        self.contributors_buf.clear();
        let contributors = &mut self.contributors_buf;

        for active_read in self.active_reads.iter() {
            let events_at_pos = active_read.cursor.events_at(walker_pos, &active_read.read);

            if events_at_pos.is_empty() {
                continue;
            }

            // Fallback to BQ=0 when this contributor has only indel events at
            // walker_pos (no Match). BQ=0 → ln(P_err)=0 in `phred_to_ln_perr`,
            // so the contributor's Match-side q_sum contribution is zero; the
            // indel BQ itself flows through events_at_pos and is folded
            // separately. Not a recovered error — there is no Match BQ here.
            let bq_at_walker = events_at_pos
                .iter()
                .find_map(|e| match e {
                    ReadEvent::Match { bq_baq, .. } => Some(*bq_baq),
                    ReadEvent::Insertion { .. } | ReadEvent::Deletion { .. } => None,
                })
                .unwrap_or(0);

            contributors.push(ReadContribution {
                read_id: active_read.read_id,
                chain_slot_id: active_read.chain_slot_id,
                events_at_pos,
                bq_baq_at_walker_pos: bq_at_walker,
                mq_log_err: active_read.read.mq_log_err,
                is_reverse_strand: active_read.read.is_reverse_strand,
                alignment_start: active_read.read.alignment_start,
                mate_role: active_read.read.mate_role,
                bq_zero_in_window: false,
                bq_override_at_walker_pos: None,
            });
        }

        // Step 2: resolve mate overlap on events at this walker
        // position. For each pair of contributors whose reads
        // share a chain_slot_id, compare BQ; the lower-BQ side
        // has its `bq_baq_at_walker_pos` zeroed (still one
        // observation, contributing ln(1)=0 log-likelihood mass)
        // and is flagged so any window event the fold pulls from
        // its cursor also gets BQ-zeroed.
        resolve_mate_overlap_at_pos(contributors, &mut self.summary);

        // Step 2b: per-column depth cap. Adopted from samtools'
        // mpileup (see `WalkerConfig` doc-comment). Apply *after*
        // mate-overlap so the cap counts genuine post-collapse
        // observations, not per-mate. Detect "indel column" from
        // the post-collapse contributor events — any Insertion or
        // Deletion at this anchor flips the column to the tighter
        // indel cap. Reads in the active set are not
        // allele-correlated in iteration order, so a deterministic
        // truncate-to-first-N is approximately unbiased and avoids
        // the random-sample machinery a per-allele clip would
        // require.
        let cap = column_depth_cap(contributors, &self.config);
        if contributors.len() > cap {
            contributors.truncate(cap);
            self.summary.column_depth_truncations += 1;
        }

        // Step 3–6: fold contributors into the records affected
        // at this walker_pos. The fold queries each contributor's
        // cursor through `&self.active`; the returned outcome
        // counts only the records that *widened* during this
        // call — fresh opens and re-finds against an
        // already-large-enough footprint do not count.
        let outcome = process_position(
            &mut self.open_records,
            walker_pos,
            self.chrom_id,
            contributors,
            &self.active_reads,
            ref_fetcher,
        )?;
        self.summary.record_widen_events += outcome.widen_count;

        Ok(())
    }

    /// Finalise any open records whose footprint is fully behind
    /// the walker, append them to `out` in emission order, and
    /// stamp the slot-allocator's pending lifecycle marks onto the
    /// first record of this tick's batch.
    ///
    /// "First of this tick's batch" is indexed via `out`'s pre-call
    /// length: marks land on `out[batch_start]`, which is the same
    /// position the previous push-based code stamped onto
    /// (`records[0]` of the per-tick local `Vec`). Returns without
    /// touching `out` if there are no aged records to drain.
    fn close_aged_records_into(&mut self, out: &mut VecDeque<PileupRecord>) {
        self.open_records
            .drain_aged_into(self.walker_pos, &mut self.drained_buf);
        if self.drained_buf.is_empty() {
            return;
        }
        let batch_start = out.len();
        // `finalise()` consumes each `OpenPileupRecord` by value, so
        // we drain the hoisted buffer rather than `into_iter()`ing
        // it; the backing `Vec` stays allocated and reusable.
        for open in self.drained_buf.drain(..) {
            out.push_back(open.finalise());
            self.summary.records_emitted += 1;
        }
        let (new, expired) = self.slots.drain_lifecycle_marks();
        if let Some(first) = out.get_mut(batch_start) {
            first.new_chains = new;
            first.expired_chains = expired;
        }
    }

    fn expire_passed_reads(&mut self) -> Result<(), WalkerError> {
        // Reads whose alignment_end < walker_pos can no longer
        // contribute. Their last event (if any) was processed at
        // their alignment_end position; once walker advances past
        // that they're done.
        debug_assert!(
            self.walker_pos > 0,
            "walker_pos starts at 1 and never decreases below 1",
        );
        self.active_reads
            .expire_passed(self.walker_pos, &mut self.slots)
    }

    fn advance(&mut self, next_pulled: Option<&PreparedRead>) -> Result<(), WalkerError> {
        // Default: one position forward, so any active read's REF
        // contribution gets folded at every position it sits on.
        // `checked_add` guards against the (extreme) case of a
        // chromosome longer than `u32::MAX` bp — silent wrap to 0
        // would corrupt all subsequent record positions. Realistic
        // mammal/plant genomes don't approach this, but a few
        // salamander/lungfish genomes do (≥ 4 Gbp).
        let mut next_pos = self
            .walker_pos
            .checked_add(1)
            .ok_or_else(|| WalkerError::Internal {
                detail: format!("walker_pos overflowed u32 at {}", self.walker_pos),
                qname: String::new(),
                chrom_id: self.chrom_id,
                pos: self.walker_pos,
            })?;

        // If the active set is empty and the next pulled read
        // starts past the walker, skip the uncovered span.
        if self.active_reads.is_empty()
            && let Some(peeked_read) = next_pulled
            && peeked_read.chrom_id == self.chrom_id
            && peeked_read.alignment_start > self.walker_pos
        {
            next_pos = peeked_read.alignment_start;
        }

        self.walker_pos = next_pos;
        Ok(())
    }

    /// Finalise everything still in flight at end-of-chromosome
    /// (or end-of-input), appending the records to `out` in
    /// emission order and stamping the slot-allocator's pending
    /// lifecycle marks onto the first record of the batch.
    fn flush_chromosome_into(
        &mut self,
        out: &mut VecDeque<PileupRecord>,
    ) -> Result<(), WalkerError> {
        // Drain remaining open records (anything that was still
        // open at end-of-chromosome is by definition ready to
        // close — there are no future reads on this chromosome).
        let batch_start = out.len();
        for open in self.open_records.drain_all() {
            out.push_back(open.finalise());
            self.summary.records_emitted += 1;
        }
        // Release any active-set reads first, then drain marks
        // once. flush_all only emits `expired_marks` (release_slot
        // never touches `new_marks`), so we don't need to merge
        // a pre-flush drain with a post-flush drain.
        self.active_reads
            .flush_all(&mut self.slots, self.walker_pos)?;
        let (new, expired) = self.slots.drain_lifecycle_marks();
        if let Some(first) = out.get_mut(batch_start) {
            first.new_chains = new;
            first.expired_chains = expired;
        }
        // Reset chromosome-scoped state. `self.open_records.reset()`
        // keeps the perf-hoisted `allele_seq_buf` capacity across the
        // chromosome boundary (Mi11).
        self.slots.reset();
        self.active_reads.reset();
        self.open_records.reset();
        Ok(())
    }

    fn summary(&self) -> RunSummary {
        self.summary.merge_slot_counters(self.slots.counters())
    }
}

/// Wrap a [`ReadLengthError`] (which carries only raw lengths) into
/// a [`WalkerError::MalformedRead`] with the offending read's locus
/// context attached.
fn malformed_read_from_length_err(err: ReadLengthError, read: &PreparedRead) -> WalkerError {
    let reason = match err {
        ReadLengthError::SeqBqMismatch {
            seq_len,
            bq_baq_len,
        } => format!("seq.len ({seq_len}) != bq_baq.len ({bq_baq_len})"),
        ReadLengthError::CigarSeqMismatch {
            cigar_consumed,
            seq_len,
        } => format!("CIGAR consumes {cigar_consumed} read bases but seq.len = {seq_len}"),
    };
    WalkerError::MalformedRead {
        reason,
        qname: read.qname.to_string(),
        chrom_id: read.chrom_id,
        pos: read.alignment_start,
    }
}

/// Resolve mate-overlap at the current walker position.
///
/// Two regimes, distinguished by whether either side carries an
/// indel anchored at this position:
///
/// - **Match-only overlap.** Both mates have only `Match` events
///   at this anchor. The lower-BQ side has its event BQs zeroed
///   in the local fold (so its `q_sum` contribution becomes
///   `ln(1) = 0`); both still count as observations and both
///   contribute the shared chain slot. Tie-break: first mate
///   keeps its BQ.
///
/// - **Indel overlap.** Either both mates report an indel at the
///   same anchor, or one reports an indel and the other a clean
///   Match (mates disagree on indel presence). The pair collapses
///   to a single observation: the loser is removed from the
///   contributor list at this walker step, so it contributes
///   nothing to the anchor record. Tie-break on `bq_baq_at_walker_pos`
///   (Match BQ where present, indel `bq_proxy` mapped through 0
///   when the loser carries no Match here); ties go to the first
///   mate.
// `&mut Vec<_>` is intentional: the function `swap_remove`s
// contributors on indel-overlap, which requires the owning `Vec`,
// not a slice.
#[allow(clippy::ptr_arg)]
fn resolve_mate_overlap_at_pos(contributors: &mut Vec<ReadContribution>, summary: &mut RunSummary) {
    // Fast path: mate overlap requires two contributors at this
    // walker_pos sharing a chain_slot_id. Solo-read inputs never
    // hit it; in paired-end inputs the geometry of insert sizes
    // means most positions also don't. Detecting the no-pair case
    // up front lets us skip the AHashMap allocation that would
    // otherwise fire on every walker step (1.5 M allocations on
    // the `pileup_walker_multi_op/L=5000` fixture, dhat 2026-05-10).
    // O(n²) early-break duplicate check; at typical n ≤ ~30
    // contributors per column it beats AHashMap construction +
    // n probes, and short-circuits as soon as any pair is found.
    let n = contributors.len();
    let any_shared = (0..n).any(|i| {
        let s = contributors[i].chain_slot_id;
        ((i + 1)..n).any(|j| contributors[j].chain_slot_id == s)
    });
    if !any_shared {
        return;
    }

    // Build a small index: chain_slot_id → list of contributor
    // indices. Anything with a list length >= 2 is a candidate.
    // ahash::AHashMap matches the rest of the module — std HashMap's
    // RandomState would make iteration non-deterministic between runs
    // and is also slower for this hot path. Mi4 in
    // `ia/reviews/pileup_2026-05-09.md`.
    let mut by_slot: AHashMap<SlotId, Vec<usize>> = AHashMap::new();
    for (i, c) in contributors.iter().enumerate() {
        by_slot.entry(c.chain_slot_id).or_default().push(i);
    }

    // Indices to discard outright (indel-overlap losers).
    let mut to_remove: Vec<usize> = Vec::new();
    // (idx, new_bq_at_walker_pos, zero_in_window) — applied to
    // each contributor of a match-only overlap pair (S7). Agree-
    // case keeper gets the summed BQ (capped at 200, zero_in_window
    // = false); disagree-case winner gets `0.8 * bq` truncated
    // (zero_in_window = false). Other / loser gets new_bq=0 with
    // zero_in_window=true. The fold honours `bq_zero_in_window`
    // (zeros every window event from this contributor's cursor)
    // and `bq_override_at_walker_pos` (rewrites walker_pos events'
    // BQ on top of the cursor pull).
    let mut bq_updates: Vec<(usize, u8, bool /*zero_in_window*/)> = Vec::new();

    for indices in by_slot.values() {
        if indices.len() < 2 {
            continue;
        }
        // Spec invariant: only mate pairs share a slot, so at most
        // two contributors per slot. Assert here so a future
        // change that admits a third reader of the same slot
        // (e.g. supplementary alignments slipping past upstream
        // filters) surfaces in tests instead of in production.
        debug_assert!(
            indices.len() <= 2,
            "more than two contributors share chain_slot_id {:?}",
            indices,
        );
        // All-pairs comparison so a future relaxation of the
        // invariant doesn't silently miss the (i, j>i+1) cases
        // that `indices.windows(2)` skips.
        for i in 0..indices.len() {
            for j in (i + 1)..indices.len() {
                let (a, b) = (indices[i], indices[j]);
                summary.mate_overlap_positions += 1;
                let any_indel_here = pair_has_indel(&contributors[a], &contributors[b]);
                if any_indel_here {
                    // Indel on either side at this walker_pos:
                    // collapse to a single observation by removing
                    // the loser entirely. Tie-break: BQ first,
                    // then first-of-pair, then alignment_start.
                    let loser_idx = pick_overlap_loser(contributors, a, b);
                    to_remove.push(loser_idx);
                } else {
                    // Match-only mate overlap (S7): apply
                    // samtools-style BQ math.
                    //
                    // PANIC-FREE: inside the !any_indel_here branch, every
                    // event at walker_pos on either side is a Match by
                    // definition of `pair_has_indel`, so `match_base_at_pos`
                    // returns Some.
                    let base_a = match_base_at_pos(&contributors[a])
                        .expect("match-only overlap: each side has a Match event at walker_pos");
                    let base_b = match_base_at_pos(&contributors[b])
                        .expect("match-only overlap: each side has a Match event at walker_pos");
                    if base_a == base_b {
                        // Agree case: sum BQs (cap 200), keeper
                        // takes the sum, other is zeroed.
                        let combined_bq = sum_bq_capped_at_200(
                            contributors[a].bq_baq_at_walker_pos,
                            contributors[b].bq_baq_at_walker_pos,
                        );
                        let keeper_idx = pick_agree_keeper(contributors, a, b);
                        let other_idx = if keeper_idx == a { b } else { a };
                        bq_updates.push((keeper_idx, combined_bq, false));
                        bq_updates.push((other_idx, 0, true));
                    } else {
                        // Disagree case: higher-BQ side keeps its
                        // BQ scaled by 0.8 (samtools' "we trust
                        // this less" haircut); loser zeroed.
                        let winner_idx = pick_disagree_winner(contributors, a, b);
                        let loser_idx = if winner_idx == a { b } else { a };
                        let scaled_bq =
                            scale_bq_by_0_8(contributors[winner_idx].bq_baq_at_walker_pos);
                        bq_updates.push((winner_idx, scaled_bq, false));
                        bq_updates.push((loser_idx, 0, true));
                    }
                }
            }
        }
    }

    // Apply bq updates in place. The fold honours both
    // `bq_zero_in_window` (zeroing every window event from this
    // contributor's cursor) and `bq_override_at_walker_pos`
    // (rewriting walker_pos events' BQ on top of the cursor
    // pull). Update the local contribution's `bq_baq_at_walker_pos`
    // and `events_at_pos` for consistency with the override.
    for (idx, new_bq, zero_in_window) in bq_updates {
        contributors[idx].bq_baq_at_walker_pos = new_bq;
        for ev in contributors[idx].events_at_pos.iter_mut() {
            set_match_event_bq(ev, new_bq);
        }
        if zero_in_window {
            contributors[idx].bq_zero_in_window = true;
        } else {
            contributors[idx].bq_override_at_walker_pos = Some(new_bq);
        }
    }

    // Drop indel-overlap losers from the contributor list.
    // Sort descending so swap_remove keeps earlier indices valid.
    to_remove.sort_unstable();
    to_remove.dedup();
    for idx in to_remove.into_iter().rev() {
        contributors.swap_remove(idx);
    }
}

/// Loser-selection for the indel-overlap case. BQ first, then
/// first-of-pair, then `alignment_start`. Matches the pre-S7
/// semantics (which combined match-only and indel paths behind
/// the same loser-selection).
fn pick_overlap_loser(contributors: &[ReadContribution], a: usize, b: usize) -> usize {
    let bq_a = contributors[a].bq_baq_at_walker_pos;
    let bq_b = contributors[b].bq_baq_at_walker_pos;
    match bq_a.cmp(&bq_b) {
        std::cmp::Ordering::Less => a,
        std::cmp::Ordering::Greater => b,
        std::cmp::Ordering::Equal => {
            let a_first = contributors[a].mate_role.is_first_of_pair();
            let b_first = contributors[b].mate_role.is_first_of_pair();
            match (a_first, b_first) {
                (true, false) => b,
                (false, true) => a,
                _ => {
                    if contributors[a].alignment_start <= contributors[b].alignment_start {
                        b
                    } else {
                        a
                    }
                }
            }
        }
    }
}

/// Keeper-selection for the agree case (S7). The choice is
/// statistically irrelevant — the surviving side carries the
/// summed BQ regardless — but must be deterministic. samtools
/// uses a qname hash; we mirror our existing tie-break logic
/// (first-of-pair, then `alignment_start`).
fn pick_agree_keeper(contributors: &[ReadContribution], a: usize, b: usize) -> usize {
    let a_first = contributors[a].mate_role.is_first_of_pair();
    let b_first = contributors[b].mate_role.is_first_of_pair();
    match (a_first, b_first) {
        (true, false) => a,
        (false, true) => b,
        _ => {
            if contributors[a].alignment_start <= contributors[b].alignment_start {
                a
            } else {
                b
            }
        }
    }
}

/// Winner-selection for the disagree case (S7): higher BQ wins;
/// ties fall back to first-of-pair, then `alignment_start`
/// (samtools uses a qname hash on ties).
fn pick_disagree_winner(contributors: &[ReadContribution], a: usize, b: usize) -> usize {
    let bq_a = contributors[a].bq_baq_at_walker_pos;
    let bq_b = contributors[b].bq_baq_at_walker_pos;
    match bq_a.cmp(&bq_b) {
        std::cmp::Ordering::Greater => a,
        std::cmp::Ordering::Less => b,
        std::cmp::Ordering::Equal => pick_agree_keeper(contributors, a, b),
    }
}

/// Extract the base from the `Match` event in `events_at_pos`.
/// In a match-only mate-overlap (no indel anchored at walker_pos
/// on either side), each contributor has exactly one Match event
/// at walker_pos.
fn match_base_at_pos(c: &ReadContribution) -> Option<u8> {
    c.events_at_pos.iter().find_map(|e| match e {
        ReadEvent::Match { base, .. } => Some(*base),
        ReadEvent::Insertion { .. } | ReadEvent::Deletion { .. } => None,
    })
}

/// `min(a + b, 200)` in u8 space without overflow. Cap from
/// samtools (`tweak_overlap_quality` in
/// [`htslib/sam.c:5919-5921`](../../../htslib/sam.c#L5919-L5921)) —
/// quality values above ~Q200 are effectively meaningless.
fn sum_bq_capped_at_200(a: u8, b: u8) -> u8 {
    let sum = (a as u16) + (b as u16);
    sum.min(200) as u8
}

/// `(bq * 0.8)` truncated to u8, matching samtools' C `0.8 *
/// uint8_t` cast at
/// [`htslib/sam.c:5927`](../../../htslib/sam.c#L5927) (truncation,
/// not rounding).
fn scale_bq_by_0_8(bq: u8) -> u8 {
    (bq as f64 * 0.8) as u8
}

/// In-place BQ rewrite on a `Match` event. No-op on indel
/// events — the S7 BQ math only applies to match-only overlaps.
fn set_match_event_bq(ev: &mut ReadEvent, bq: u8) {
    if let ReadEvent::Match { bq_baq, .. } = ev {
        *bq_baq = bq;
    }
}

/// True iff at least one of the two contributors has an
/// Insertion or Deletion anchored at the current walker_pos.
fn pair_has_indel(a: &ReadContribution, b: &ReadContribution) -> bool {
    let has_indel = |c: &ReadContribution| {
        c.events_at_pos
            .iter()
            .any(|e| matches!(e, ReadEvent::Insertion { .. } | ReadEvent::Deletion { .. },))
    };
    has_indel(a) || has_indel(b)
}

/// Per-column depth cap. Returns the lower indel cap if any
/// contributor reports an Insertion or Deletion at this anchor;
/// otherwise the SNP/REF cap.
fn column_depth_cap(contributors: &[ReadContribution], config: &WalkerConfig) -> usize {
    let any_indel = contributors.iter().any(|c| {
        c.events_at_pos
            .iter()
            .any(|e| matches!(e, ReadEvent::Insertion { .. } | ReadEvent::Deletion { .. },))
    });
    if any_indel {
        config.max_indel_column_depth as usize
    } else {
        config.max_snp_column_depth as usize
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::super::cigar_cursor::EventsAt;
    use super::*;

    fn contribution(
        bq: u8,
        is_first_mate: bool,
        alignment_start: u32,
        events: EventsAt,
    ) -> ReadContribution {
        ReadContribution {
            read_id: 0,
            chain_slot_id: 0,
            events_at_pos: events,
            bq_baq_at_walker_pos: bq,
            mq_log_err: -3.0,
            is_reverse_strand: false,
            alignment_start,
            mate_role: if is_first_mate {
                super::super::MateRole::FirstOfPair
            } else {
                super::super::MateRole::SecondOfPair
            },
            bq_zero_in_window: false,
            bq_override_at_walker_pos: None,
        }
    }

    fn match_evs(pos: u32, base: u8, bq: u8) -> EventsAt {
        let mut v = EventsAt::new();
        v.push(ReadEvent::Match {
            ref_pos: pos,
            base,
            bq_baq: bq,
        });
        v
    }

    fn indel_ins_evs(anchor: u32, bq: u8) -> EventsAt {
        let mut v = EventsAt::new();
        v.push(ReadEvent::Match {
            ref_pos: anchor,
            base: b'A',
            bq_baq: bq,
        });
        v.push(ReadEvent::Insertion {
            anchor_ref_pos: anchor,
            seq: b"A".to_vec(),
            bq_proxy: bq,
        });
        v
    }

    // --- M19: pick_* tertiary tie-break tests --------------------
    //
    // All three pick functions fall through to comparing
    // `alignment_start` when both contributors' first-of-pair bits
    // agree. A flipped comparison here would silently change
    // tie-break determinism — exactly the kind of bug that surfaces
    // as "different VCF on the same input".

    #[test]
    fn pick_agree_keeper_breaks_remaining_tie_by_earlier_alignment_start() {
        // Both first-mate, BQ tie (no BQ check in this function).
        // Earlier alignment_start (50) wins over later (100).
        let c = vec![
            contribution(30, true, 100, match_evs(1, b'A', 30)),
            contribution(30, true, 50, match_evs(1, b'A', 30)),
        ];
        assert_eq!(pick_agree_keeper(&c, 0, 1), 1);
        // Swap order: now index 0 has the earlier alignment_start.
        let c = vec![
            contribution(30, true, 50, match_evs(1, b'A', 30)),
            contribution(30, true, 100, match_evs(1, b'A', 30)),
        ];
        assert_eq!(pick_agree_keeper(&c, 0, 1), 0);
    }

    #[test]
    fn pick_overlap_loser_breaks_bq_and_first_mate_tie_by_alignment_start() {
        // BQ tie + first-mate tie → the loser is the one with the
        // larger alignment_start.
        let c = vec![
            contribution(30, true, 100, match_evs(1, b'A', 30)),
            contribution(30, true, 50, match_evs(1, b'A', 30)),
        ];
        // a=0 (start 100), b=1 (start 50) → loser is a (later start).
        assert_eq!(pick_overlap_loser(&c, 0, 1), 0);
    }

    #[test]
    fn pick_disagree_winner_on_bq_tie_delegates_to_pick_agree_keeper() {
        // BQ tie + first-mate tie → falls back to alignment_start.
        // The winner under the agree-keeper rule is the earlier
        // alignment_start.
        let c = vec![
            contribution(30, true, 100, match_evs(1, b'A', 30)),
            contribution(30, true, 50, match_evs(1, b'A', 30)),
        ];
        assert_eq!(pick_disagree_winner(&c, 0, 1), 1);
    }

    // --- Boundary tests for the samtools-C parity helpers --------

    #[test]
    fn sum_bq_capped_at_200_caps_exactly_at_200() {
        assert_eq!(sum_bq_capped_at_200(0, 0), 0);
        assert_eq!(sum_bq_capped_at_200(100, 100), 200);
        assert_eq!(sum_bq_capped_at_200(150, 100), 200);
        assert_eq!(sum_bq_capped_at_200(255, 255), 200);
        assert_eq!(sum_bq_capped_at_200(99, 100), 199);
    }

    #[test]
    fn scale_bq_by_0_8_truncates_not_rounds() {
        assert_eq!(scale_bq_by_0_8(0), 0);
        assert_eq!(scale_bq_by_0_8(5), 4); // 4.0 exact
        assert_eq!(scale_bq_by_0_8(7), 5); // 5.6 → trunc 5 (round would give 6)
        assert_eq!(scale_bq_by_0_8(30), 24); // 24.0
    }

    // --- column_depth_cap: any-indel rule ------------------------

    #[test]
    fn column_depth_cap_returns_indel_cap_when_only_some_contributors_have_indel() {
        // Mixed SNP + one indel contributor at the same anchor.
        // The "any" rule must flip the column to the indel cap.
        let cfg = WalkerConfig {
            max_snp_column_depth: 8000,
            max_indel_column_depth: 2,
            ..WalkerConfig::default()
        };
        let v = vec![
            contribution(30, true, 1, match_evs(1, b'A', 30)),
            contribution(30, true, 1, indel_ins_evs(1, 30)),
            contribution(30, true, 1, match_evs(1, b'A', 30)),
            contribution(30, true, 1, match_evs(1, b'A', 30)),
            contribution(30, true, 1, match_evs(1, b'A', 30)),
        ];
        assert_eq!(column_depth_cap(&v, &cfg), 2);
    }
}
