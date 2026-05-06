//! Open-record formation, merging, and widening — the
//! walker's central operation. See `ia/specs/pileup_walker.md`
//! §"Open-record formation and merging".
//!
//! An open record is the in-flight version of a `PileupRecord`:
//! same shape, but its REF span and allele list grow as the walker
//! folds events into it. When the walker confirms no future event
//! can touch the record (its footprint is fully behind the walker
//! per the closure rule), it converts the open record into a
//! finalised `PileupRecord` and pushes it through the channel.

use std::collections::BTreeMap;

use ahash::AHashMap;

use super::decompose::ReadEvent;
use super::errors::WalkerError;
use super::slot_allocator::SlotId;
use super::{AlleleObservation, FiveScalars, MAX_RECORD_SPAN, PileupRecord, RefBaseFetcher};

/// One in-flight allele bucket inside an `OpenPileupRecord`.
#[derive(Debug, Clone)]
pub struct OpenAllele {
    pub seq: Vec<u8>,
    pub scalars: FiveScalars,
    /// Sorted, deduplicated.
    pub chain_slots: Vec<SlotId>,
}

impl OpenAllele {
    fn new(seq: Vec<u8>) -> Self {
        Self {
            seq,
            scalars: FiveScalars::zero(),
            chain_slots: Vec::new(),
        }
    }
}

/// One in-flight per-position record. Its REF span is `ref_seq.len()`;
/// no separate `ref_span` field is stored. `alleles[0]` is always
/// REF (`seq == ref_seq`).
#[derive(Debug, Clone)]
pub struct OpenPileupRecord {
    pub chrom_id: u32,
    /// 1-based anchor position.
    pub pos: u32,
    pub ref_seq: Vec<u8>,
    pub alleles: Vec<OpenAllele>,
    /// Per-read fold state — the contribution this record currently
    /// holds for each read that has folded into it. Used to enforce
    /// "fold each (record, read) pair exactly once" across walker
    /// steps: at re-fold time the previous contribution is
    /// subtracted from its old bucket before the new contribution is
    /// added to the new bucket. Without this state, a read with
    /// Match events at every position inside an open record's
    /// footprint would be re-folded once per walker step inside
    /// that footprint, multiplying every five-scalar value by
    /// `ref_span` (B1 in `ia/reviews/pileup_2026-05-06.md`).
    folded_reads: AHashMap<u32, FoldedReadState>,
}

/// What a single read currently contributes to one bucket of one
/// open record. Carries enough state to subtract the contribution
/// cleanly when the read re-folds (e.g. on widening that grows the
/// haplotype seq under the record).
#[derive(Debug, Clone, Copy)]
struct FoldedReadState {
    allele_index: usize,
    contribution: FiveScalars,
}

impl OpenPileupRecord {
    /// Open a fresh record at `pos` with the given initial REF
    /// sequence. The REF allele bucket is created up front (with
    /// zero observations) so the `alleles[0] == REF` invariant
    /// holds from the very start.
    fn new(chrom_id: u32, pos: u32, ref_seq: Vec<u8>) -> Self {
        let ref_allele = OpenAllele::new(ref_seq.clone());
        Self {
            chrom_id,
            pos,
            ref_seq,
            alleles: vec![ref_allele],
            folded_reads: AHashMap::new(),
        }
    }

    pub fn ref_span(&self) -> u32 {
        self.ref_seq.len() as u32
    }

    /// Footprint end (exclusive), in 1-based coordinates: the
    /// position one past the last reference base this record
    /// covers.
    pub fn footprint_end_exclusive(&self) -> u32 {
        self.pos + self.ref_span()
    }

    /// Convert into a finalised `PileupRecord`. The slot lifecycle
    /// markers (`new_chains` / `expired_chains`) are filled by the
    /// caller from the `SlotAllocator`'s drain; this method only
    /// converts the per-allele state.
    pub fn finalise(self) -> PileupRecord {
        let alleles = self
            .alleles
            .into_iter()
            .map(|a| AlleleObservation {
                seq: a.seq,
                scalars: a.scalars,
                chain_slots: a.chain_slots,
            })
            .collect();
        PileupRecord {
            chrom_id: self.chrom_id,
            pos: self.pos,
            new_chains: Vec::new(),
            expired_chains: Vec::new(),
            alleles,
        }
    }
}

/// The set of currently-open records, keyed by anchor position.
/// Range queries (find records overlapping a given event span)
/// use the BTreeMap's ordered structure.
#[derive(Debug, Default)]
pub struct OpenPileupRecordTable {
    /// 1-based anchor position → record.
    records: BTreeMap<u32, OpenPileupRecord>,
}

impl OpenPileupRecordTable {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&u32, &OpenPileupRecord)> {
        self.records.iter()
    }

    /// Drain every record whose footprint is fully behind the
    /// walker (`pos + ref_span ≤ walker_pos`), in coordinate
    /// order. Used by the walker's `close_aged_records` step.
    pub fn drain_aged(&mut self, walker_pos: u32) -> Vec<OpenPileupRecord> {
        let mut closing_keys: Vec<u32> = Vec::new();
        for (&pos, rec) in self.records.iter() {
            if rec.footprint_end_exclusive() <= walker_pos {
                closing_keys.push(pos);
            } else {
                // BTreeMap iteration is sorted; once we see a
                // record whose footprint is still ahead of the
                // walker, no later key can be aged either, since
                // their `pos` is even larger.
                break;
            }
        }
        let mut out = Vec::with_capacity(closing_keys.len());
        for pos in closing_keys {
            if let Some(rec) = self.records.remove(&pos) {
                out.push(rec);
            }
        }
        out
    }

    /// Drain everything unconditionally (chromosome boundary or
    /// end-of-input). Records come out in coordinate order.
    pub fn drain_all(&mut self) -> Vec<OpenPileupRecord> {
        let mut out = Vec::with_capacity(self.records.len());
        let keys: Vec<u32> = self.records.keys().copied().collect();
        for k in keys {
            if let Some(r) = self.records.remove(&k) {
                out.push(r);
            }
        }
        out
    }

    /// Find the open record (if any) whose footprint overlaps the
    /// half-open interval `[event_start, event_end)`. "Overlap"
    /// here is non-empty interval intersection — touching
    /// intervals are not overlapping. Returns the anchor position
    /// of the matched record.
    pub fn find_overlapping(&self, event_start: u32, event_end: u32) -> Option<u32> {
        // Candidates are records whose anchor `Q ≤ event_start`
        // (any record opened to the right of the event's start
        // would have its footprint start ≥ event_end > event_start
        // — they can't overlap). Scan from the largest Q ≤
        // event_start downward and return the first whose footprint
        // reaches into the event.
        for (&q, rec) in self.records.range(..=event_start).rev() {
            if rec.footprint_end_exclusive() > event_start && q < event_end {
                return Some(q);
            }
            // If this record's footprint ends at or before
            // event_start, no earlier record's footprint can reach
            // event_start either (footprints are bounded by
            // MAX_RECORD_SPAN; in practice we'd still need to walk
            // a couple back to be safe, but for correctness we
            // just terminate the search).
            if rec.footprint_end_exclusive() <= event_start {
                break;
            }
        }
        None
    }

    /// Widen the record at `key` so its REF span covers up to
    /// `new_end_exclusive`, fetching the additional reference
    /// bases. Existing alleles are rewritten by appending the
    /// new reference bases (alleles whose previous coverage
    /// already ran to the end of the old span); deletion alleles
    /// keep their existing length, expressing "more bases deleted
    /// relative to the wider REF" — see
    /// `ia/specs/pileup_walker.md` §"Step 4.2".
    fn widen(
        &mut self,
        key: u32,
        new_end_exclusive: u32,
        fasta: &dyn RefBaseFetcher,
    ) -> Result<(), WalkerError> {
        let rec = self
            .records
            .get_mut(&key)
            .expect("widen called on absent record");
        let old_end = rec.footprint_end_exclusive();
        if new_end_exclusive <= old_end {
            return Ok(());
        }
        let extra_len = new_end_exclusive - old_end;
        if (new_end_exclusive - rec.pos) > MAX_RECORD_SPAN {
            return Err(WalkerError::RecordTooWide {
                chrom_id: rec.chrom_id,
                pos: rec.pos,
                span: new_end_exclusive - rec.pos,
                cap: MAX_RECORD_SPAN,
            });
        }
        let extra_bases = fasta
            .fetch(rec.chrom_id, old_end, extra_len)
            .map_err(|source| WalkerError::Fasta {
                chrom_id: rec.chrom_id,
                start: old_end,
                start_plus_len: new_end_exclusive,
                source,
            })?;
        rec.ref_seq.extend_from_slice(&extra_bases);

        // Rewrite each existing allele.
        for allele in &mut rec.alleles {
            // If the allele was ref-aligned over the OLD span (i.e.
            // its seq length equals the old ref span), append the
            // new reference bases. Otherwise — the allele was a
            // deletion (shorter than old span) or insertion
            // (longer than old span by the inserted bases) — the
            // rewrite rule preserves the event:
            //
            //  - DEL allele: keep its current length. The widened
            //    REF now has more bases past the deletion's end,
            //    but the deletion's seq stays "anchor only" or
            //    "what remains after the deletion was applied to
            //    the OLD REF span." We still need to append the
            //    new ref bases, because the ref bases past the
            //    deletion are present in the read (and therefore
            //    in the haplotype).
            //  - INS allele: same — append new ref bases past the
            //    insertion's REF coverage.
            //
            // In all cases, append the new REF bases to the end
            // of the allele's seq, because the new reference
            // stretch is conserved (no event in this allele claims
            // those bases as modified).
            allele.seq.extend_from_slice(&extra_bases);
        }
        Ok(())
    }

    /// Open a fresh record at `pos` with REF span `span`, fetching
    /// the reference bases from `fasta`.
    fn open_new(
        &mut self,
        chrom_id: u32,
        pos: u32,
        span: u32,
        fasta: &dyn RefBaseFetcher,
    ) -> Result<&mut OpenPileupRecord, WalkerError> {
        if span > MAX_RECORD_SPAN {
            return Err(WalkerError::RecordTooWide {
                chrom_id,
                pos,
                span,
                cap: MAX_RECORD_SPAN,
            });
        }
        let ref_seq = fasta
            .fetch(chrom_id, pos, span)
            .map_err(|source| WalkerError::Fasta {
                chrom_id,
                start: pos,
                start_plus_len: pos + span,
                source,
            })?;
        let rec = OpenPileupRecord::new(chrom_id, pos, ref_seq);
        self.records.insert(pos, rec);
        Ok(self.records.get_mut(&pos).expect("just inserted"))
    }
}

/// Apply a list of events from one read to the open record's REF
/// sequence to compute the haplotype string this read presents
/// under that record. Events are assumed to be inside the record's
/// footprint and sorted by anchor.
pub fn apply_events_to_ref(record_pos: u32, ref_seq: &[u8], events: &[&ReadEvent]) -> Vec<u8> {
    // Walk the ref_seq in offset order, applying events as we
    // pass their anchor positions. Each event is positioned by an
    // offset = (event.anchor_pos - record_pos).
    let mut out: Vec<u8> = Vec::with_capacity(ref_seq.len() + 8);

    // Sort events by anchor for safety (caller should have, but we
    // re-sort since the cost is tiny and the contract is clearer).
    let mut sorted: Vec<&ReadEvent> = events.to_vec();
    sorted.sort_by_key(|e| e.anchor_pos());

    // Skip indices already consumed by an event so we don't
    // double-emit reference bases that an event has overridden.
    let mut consumed_until: u32 = 0; // ref offset (exclusive) consumed by the last event
    let mut ref_cursor: u32 = 0;

    for ev in sorted {
        let offset = ev.anchor_pos().checked_sub(record_pos).unwrap_or(0);

        // Emit any REF bases between ref_cursor and the event's
        // offset that haven't been consumed by a previous event.
        while ref_cursor < offset {
            if ref_cursor >= consumed_until && (ref_cursor as usize) < ref_seq.len() {
                out.push(ref_seq[ref_cursor as usize]);
            }
            ref_cursor += 1;
        }

        match ev {
            ReadEvent::Match { base, .. } => {
                if (offset as usize) < ref_seq.len() {
                    out.push(*base);
                }
                ref_cursor = offset + 1;
                consumed_until = consumed_until.max(offset + 1);
            }
            ReadEvent::Insertion { seq, .. } => {
                // Insertion sits AFTER the anchor base; the anchor
                // base itself is unchanged (in the ref-matching
                // sense). Emit the anchor base if not already
                // emitted, then append the inserted bases.
                if (offset as usize) < ref_seq.len() && offset >= consumed_until {
                    out.push(ref_seq[offset as usize]);
                }
                out.extend_from_slice(seq);
                ref_cursor = offset + 1;
                consumed_until = consumed_until.max(offset + 1);
            }
            ReadEvent::Deletion { deleted_len, .. } => {
                // DEL: keep the anchor base, drop the next
                // `deleted_len` reference bases.
                if (offset as usize) < ref_seq.len() && offset >= consumed_until {
                    out.push(ref_seq[offset as usize]);
                }
                let skip_until = offset + 1 + *deleted_len;
                ref_cursor = skip_until;
                consumed_until = consumed_until.max(skip_until);
            }
        }
    }

    // Tail: emit remaining REF bases.
    while (ref_cursor as usize) < ref_seq.len() {
        if ref_cursor >= consumed_until {
            out.push(ref_seq[ref_cursor as usize]);
        }
        ref_cursor += 1;
    }

    out
}

/// Find or create the allele bucket inside a record matching `seq`.
/// Returns the bucket index. Linear scan is fine — records
/// typically carry ≤ a few alleles. Returning the index rather
/// than `&mut OpenAllele` keeps the `rec` borrow short, so callers
/// can index into `alleles` and update sibling state (e.g.
/// `folded_reads`) in the same scope.
pub fn find_or_create_allele_index(rec: &mut OpenPileupRecord, seq: Vec<u8>) -> usize {
    if let Some(idx) = rec.alleles.iter().position(|a| a.seq == seq) {
        idx
    } else {
        rec.alleles.push(OpenAllele::new(seq));
        rec.alleles.len() - 1
    }
}

/// Process all events at `walker_pos` from the given list of
/// per-read contributions. Each contributor MUST have at least
/// one event at `walker_pos` (anchor == walker_pos); silent reads
/// (e.g. inside their own deletion or an N-skip at this pos) are
/// already filtered out by the caller.
///
/// 1. Step 3 — identify candidate records per event. Each event
///    either merges into an existing overlapping record or opens
///    a fresh one. Record keys touched here go into `affected`.
///
/// 2. Step 4-6 — for each affected record, fold every contributor
///    whose events overlap the record's footprint exactly once
///    into the matching allele bucket. The fold uses the
///    contributor's full event list (not just events_at_pos) so
///    that compound alleles spanning the whole record's footprint
///    collapse into a single allele bucket.
///
/// Records that already existed at walker_pos but were *not*
/// affected at this step are not re-folded — those reads were
/// folded at the walker step where the record was created or
/// widened, and folding again here would double-count.
pub fn process_position(
    open: &mut OpenPileupRecordTable,
    walker_pos: u32,
    chrom_id: u32,
    contributors: &[ReadContribution],
    fasta: &dyn RefBaseFetcher,
) -> Result<(), WalkerError> {
    let mut affected: Vec<u32> = Vec::new();

    // Step 3: each event either lands in an existing record (and
    // possibly widens it) or opens a fresh one.
    for contrib in contributors {
        for ev in &contrib.events_at_pos {
            let event_start = ev.anchor_pos();
            let event_end = event_start + ev.footprint_span();

            let key = if let Some(k) = open.find_overlapping(event_start, event_end) {
                let cur_end = open
                    .records
                    .get(&k)
                    .expect("just located")
                    .footprint_end_exclusive();
                if event_end > cur_end {
                    open.widen(k, event_end, fasta)?;
                }
                k
            } else {
                let new = open.open_new(chrom_id, event_start, ev.footprint_span(), fasta)?;
                new.pos
            };
            if !affected.contains(&key) {
                affected.push(key);
            }
        }
    }

    // Step 4-6: for each affected record (in coordinate order),
    // fold each contributor that has events overlapping the
    // record's footprint. Each (record, contributor) pair folds
    // exactly once across the record's lifetime: re-folds at
    // later walker steps subtract the prior contribution from
    // the old bucket before adding the new one.
    affected.sort_unstable();
    for key in affected {
        let (rec_pos, rec_ref_seq) = {
            let rec = open.records.get(&key).expect("affected key must exist");
            (rec.pos, rec.ref_seq.clone())
        };
        let rec_end = rec_pos + rec_ref_seq.len() as u32;
        let _ = walker_pos;

        for contrib in contributors {
            let window_events: Vec<&ReadEvent> = contrib
                .full_window_events
                .iter()
                .filter(|e| {
                    let s = e.anchor_pos();
                    let span = e.footprint_span();
                    s < rec_end && (s + span) > rec_pos
                })
                .collect();

            // A contributor only folds into a record if it has
            // events overlapping the record's footprint. (No
            // events overlapping = the read doesn't observe this
            // record's REF stretch at all — it shouldn't fold.)
            if window_events.is_empty() {
                continue;
            }

            let allele_seq = apply_events_to_ref(rec_pos, &rec_ref_seq, &window_events);
            let bq = ln_bq_for_read(&window_events, contrib.bq_baq_at_walker_pos);
            let ln_q = bq.max(contrib.mq_log_err);
            let new_contribution = FiveScalars {
                num_obs: 1,
                q_sum: ln_q,
                fwd: u32::from(!contrib.is_reverse_strand),
                placed_left: u32::from(contrib.alignment_start < rec_pos),
                placed_start: u32::from(contrib.alignment_start == rec_pos),
            };

            let rec = open.records.get_mut(&key).expect("affected key must exist");
            // Subtract any prior contribution this read had to
            // this record. Re-folds happen when a record widens
            // and the read's haplotype under the wider span
            // either changes bucket or stays in the same bucket
            // with a possibly-different ln_q.
            if let Some(prev) = rec.folded_reads.remove(&contrib.read_id) {
                subtract_contribution(
                    &mut rec.alleles[prev.allele_index].scalars,
                    &prev.contribution,
                );
            }
            let new_index = find_or_create_allele_index(rec, allele_seq);
            add_contribution(&mut rec.alleles[new_index].scalars, &new_contribution);
            insert_sorted_unique(
                &mut rec.alleles[new_index].chain_slots,
                contrib.chain_slot_id,
            );
            rec.folded_reads.insert(
                contrib.read_id,
                FoldedReadState {
                    allele_index: new_index,
                    contribution: new_contribution,
                },
            );
        }
    }

    Ok(())
}

fn add_contribution(scalars: &mut FiveScalars, c: &FiveScalars) {
    scalars.num_obs += c.num_obs;
    scalars.q_sum += c.q_sum;
    scalars.fwd += c.fwd;
    scalars.placed_left += c.placed_left;
    scalars.placed_start += c.placed_start;
}

fn subtract_contribution(scalars: &mut FiveScalars, c: &FiveScalars) {
    // Saturating-style: an internal-bookkeeping bug that produced
    // a negative would otherwise wrap silently. Saturate to zero
    // and rely on the upper invariant (num_obs over the record
    // total still adds up) to surface mistakes through tests.
    scalars.num_obs = scalars.num_obs.saturating_sub(c.num_obs);
    scalars.q_sum -= c.q_sum;
    scalars.fwd = scalars.fwd.saturating_sub(c.fwd);
    scalars.placed_left = scalars.placed_left.saturating_sub(c.placed_left);
    scalars.placed_start = scalars.placed_start.saturating_sub(c.placed_start);
}

/// Per-read BQ for an allele's quality contribution: min over
/// the read's events in the record's footprint, converted to
/// `ln(P_err)`. Confirmed against freebayes
/// `AlleleParser.cpp:3151-3155` (haplotype-allele construction
/// uses `min(quality)`). For a contributor with no events in the
/// window (a clean REF read), we fall back to the read's BQ at
/// the walker_pos, which is a `Match` event's `bq_baq` already
/// stamped on the contribution.
fn ln_bq_for_read(window_events: &[&ReadEvent], fallback_bq: u8) -> f64 {
    let min_bq = window_events
        .iter()
        .filter_map(|e| match e {
            ReadEvent::Match { bq_baq, .. } => Some(*bq_baq),
            ReadEvent::Insertion { bq_proxy, .. } => Some(*bq_proxy),
            ReadEvent::Deletion { bq_proxy, .. } => Some(*bq_proxy),
        })
        .min()
        .unwrap_or(fallback_bq);
    phred_to_ln_perr(min_bq)
}

/// `Q -> ln(P_err)` where `P_err = 10^(-Q/10)`.
/// For Q = 0 the error probability is 1, so `ln(1) = 0`.
fn phred_to_ln_perr(q: u8) -> f64 {
    if q == 0 {
        return 0.0;
    }
    -(q as f64) * std::f64::consts::LN_10 / 10.0
}

fn insert_sorted_unique(v: &mut Vec<SlotId>, slot: SlotId) {
    match v.binary_search(&slot) {
        Ok(_) => {} // already present
        Err(pos) => v.insert(pos, slot),
    }
}

/// One read's contribution to the current walker position. The
/// walker assembles these from its active set before calling
/// `process_position`.
#[derive(Debug, Clone)]
pub struct ReadContribution {
    /// Active-set local id of the contributing read. Keys the
    /// per-record `folded_reads` map so re-folding a read into the
    /// same record at a later walker step subtracts its prior
    /// contribution before adding the new one.
    pub read_id: u32,
    pub chain_slot_id: SlotId,
    /// Events whose anchor *is* this walker_pos (used by step 3
    /// to identify candidate records).
    pub events_at_pos: Vec<ReadEvent>,
    /// All this read's events in a window around the walker_pos,
    /// for the haplotype-string computation in step 5. In the
    /// simple non-merged case this is identical to events_at_pos.
    pub full_window_events: Vec<ReadEvent>,
    /// BAQ-capped BQ at this walker position (the Match-event's
    /// quality, used as the fallback when window_events is empty
    /// — a clean REF read).
    pub bq_baq_at_walker_pos: u8,
    pub mq_log_err: f64,
    pub is_reverse_strand: bool,
    pub alignment_start: u32,
}

/// Apply slot lifecycle markers from the slot allocator to the
/// most-recently-emitted record. Used by the walker's
/// close_aged_records step.
pub fn stamp_lifecycle_marks(
    records: &mut [PileupRecord],
    new_chains: Vec<SlotId>,
    expired_chains: Vec<SlotId>,
) {
    if let Some(first) = records.first_mut() {
        first.new_chains = new_chains;
        first.expired_chains = expired_chains;
    }
    // If multiple records emit at once, only the first carries
    // the lifecycle markers — Stage 5's recv loop applies them
    // before reading allele observations, so attaching to the
    // first record is enough.
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::per_sample_caller::pileup::tests::MockFasta;

    fn fa(s: &str) -> MockFasta {
        MockFasta::new(s)
    }

    #[test]
    fn open_new_creates_record_with_ref_allele_zero_obs() {
        let mut t = OpenPileupRecordTable::new();
        let f = fa("ACGTAC");
        let rec = t.open_new(0, 1, 1, &f).unwrap();
        assert_eq!(rec.pos, 1);
        assert_eq!(rec.ref_seq, b"A");
        assert_eq!(rec.alleles.len(), 1);
        assert_eq!(rec.alleles[0].seq, b"A");
        assert_eq!(rec.alleles[0].scalars.num_obs, 0);
    }

    #[test]
    fn widen_extends_ref_seq_and_existing_alleles() {
        let mut t = OpenPileupRecordTable::new();
        let f = fa("ACGTAC");
        // Open at pos 1 with span 1 ("A").
        t.open_new(0, 1, 1, &f).unwrap();
        // Now widen to span 3 ("ACG").
        t.widen(1, 4, &f).unwrap();
        let rec = t.records.get(&1).unwrap();
        assert_eq!(rec.ref_seq, b"ACG");
        assert_eq!(rec.alleles[0].seq, b"ACG");
    }

    #[test]
    fn drain_aged_emits_in_coordinate_order() {
        let mut t = OpenPileupRecordTable::new();
        let f = fa("ACGTACGTAC");
        t.open_new(0, 1, 1, &f).unwrap();
        t.open_new(0, 5, 1, &f).unwrap();
        t.open_new(0, 8, 1, &f).unwrap();
        // Walker at pos 6 — record at 1 (ends 2) and at 5 (ends 6) are aged out.
        let drained = t.drain_aged(6);
        assert_eq!(drained.len(), 2);
        assert_eq!(drained[0].pos, 1);
        assert_eq!(drained[1].pos, 5);
        assert_eq!(t.records.len(), 1);
    }

    #[test]
    fn find_overlapping_returns_record_when_event_falls_inside_footprint() {
        let mut t = OpenPileupRecordTable::new();
        let f = fa("AAAACCCCGGGG");
        // Open a deletion-shaped record at pos 1 with span 5
        // (footprint [1, 6)).
        t.open_new(0, 1, 5, &f).unwrap();
        // Event at pos 3 with span 1 (a SNP) — overlaps [1, 6).
        let key = t.find_overlapping(3, 4);
        assert_eq!(key, Some(1));
        // Event at pos 6 with span 1 — does NOT overlap (touching, not overlapping).
        let key = t.find_overlapping(6, 7);
        assert_eq!(key, None);
    }

    #[test]
    fn apply_events_pure_match_yields_unchanged_ref() {
        let ref_seq = b"ACGTA";
        let out = apply_events_to_ref(100, ref_seq, &[]);
        assert_eq!(out, b"ACGTA");
    }

    #[test]
    fn apply_events_snp_replaces_one_base() {
        let ref_seq = b"ACGTA";
        let snp = ReadEvent::Match {
            ref_pos: 102,
            base: b'X',
            bq_baq: 30,
        };
        let out = apply_events_to_ref(100, ref_seq, &[&snp]);
        assert_eq!(out, b"ACXTA");
    }

    #[test]
    fn apply_events_deletion_drops_bases_after_anchor() {
        let ref_seq = b"ACGTA";
        let del = ReadEvent::Deletion {
            anchor_ref_pos: 100,
            deleted_len: 2,
            bq_proxy: 30,
        };
        let out = apply_events_to_ref(100, ref_seq, &[&del]);
        // Anchor at 100 ("A"), then 2 bases deleted, then "TA".
        assert_eq!(out, b"ATA");
    }

    #[test]
    fn apply_events_insertion_appends_inserted_bases_after_anchor() {
        let ref_seq = b"ACGTA";
        let ins = ReadEvent::Insertion {
            anchor_ref_pos: 100,
            seq: b"XX".to_vec(),
            bq_proxy: 30,
        };
        let out = apply_events_to_ref(100, ref_seq, &[&ins]);
        // Anchor "A", inserted "XX", then ref "CGTA".
        assert_eq!(out, b"AXXCGTA");
    }

    #[test]
    fn find_or_create_allele_returns_same_bucket_on_match() {
        let mut rec = OpenPileupRecord::new(0, 100, b"ACG".to_vec());
        let idx1 = find_or_create_allele_index(&mut rec, b"ACT".to_vec());
        rec.alleles[idx1].scalars.num_obs = 1;
        let idx2 = find_or_create_allele_index(&mut rec, b"ACT".to_vec());
        assert_eq!(idx1, idx2);
        assert_eq!(rec.alleles[idx2].scalars.num_obs, 1);
        // REF + ACT = 2 buckets total.
        assert_eq!(rec.alleles.len(), 2);
    }
}
