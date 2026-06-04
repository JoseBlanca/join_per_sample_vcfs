//! Cohort producer — section 1 (appendix §B).
//!
//! *(today: `var_calling::driver` `BlockIterator`, `two_pass`, `loader`)*
//!
//! `CohortChunkIntegrator` owns the N `SamplePspReader`s, the `DustAheadPool`,
//! and the REF fetcher. It **streams `CohortPileupRecord`s, sliced into chunks
//! at safe gaps** (the producer algorithm, §2.2):
//!
//! 1. lockstep-read one psp segment per sample (light columns only) up to the
//!    watermark = `min(peek_next_span)`;
//! 2. merge across samples by position → variable positions (AC / `min_alt_obs`),
//!    then apply dust;
//! 3. build one `CohortPileupRecord` per variable position (heavy columns via
//!    the readers' `take_*` getters, variable rows only);
//! 4. accumulate records; cut at a safe gap (`find_block_cut`) — between whole
//!    records, so nothing is split;
//! 5. fetch the chunk's REF span; ship the chunk tagged `chunk_order`.
//!
//! **Memory invariant:** the cohort-wide footprint is *only* the consolidated
//! records — variable positions, AC / `min_alt_obs` already applied. Never
//! materialise a full-coverage cohort structure.
//!
//! `CohortPerPositionMerge` is the cohort join — the revived
//! [`per_position_merger`](crate::var_calling_new::per_position_merger), here
//! in the producer.
//!
//! Phase 2 builds this module, in two steps:
//! - **2a** ([`CohortSpanFold`]): the cohort keep/cut math — the revived
//!   `two_pass::WindowSummary`, folding over the reader's *light* columns
//!   ([`SamplePspChunk`](crate::var_calling_new::sample_reader::SamplePspChunk)
//!   `positions`/`ref_spans`/`nonref_obs`) instead of a `SampleColumns`.
//!   Plus [`drop_dust_masked`], the dust step. Byte-identity-critical, so the
//!   arithmetic is copied verbatim from `two_pass.rs`.
//! - **2b** ([`CohortChunkIntegrator`]): the streaming integrator over **one
//!   covered interval** — segment buffering to the watermark, per-position
//!   record building (columnar→record via
//!   [`SamplePspChunk::records_for`](crate::var_calling_new::sample_reader::SamplePspChunk::records_for),
//!   merged by the revived [`PerPositionMerger`](crate::var_calling_new::per_position_merger::PerPositionMerger)),
//!   safe-gap cut, REF fetch (closure-injected), `chunk_order` stamping.
//!   The chromosome / interval iteration, per-chromosome REF fetcher, and
//!   `DustAheadPool` wiring are the next step (2b-wiring); `produce_chunk`
//!   takes the REF fetch as a closure and the dust mask per
//!   [`begin_interval`](CohortChunkIntegrator::begin_interval) so the engine
//!   stays decoupled and unit-testable.

use std::ops::Range;

/// Reach of a position given its (max) ref span — `pos + max(span, 1) - 1`,
/// saturating. Copied from `two_pass::reach`; the grouping arithmetic must
/// match byte-for-byte.
#[inline]
fn reach(pos: u32, span: u32) -> u32 {
    pos.saturating_add(span.max(1)).saturating_sub(1)
}

/// Compact per-position cohort summary over a span — the revived
/// `two_pass::WindowSummary`, folding the readers' **light** columns
/// (`positions` / `ref_spans` / `nonref_obs`) one sample at a time so only
/// the cohort summary (plus the one sample being merged) is ever resident.
///
/// Parallel sorted arrays keyed by [`positions`](Self::positions). The keep /
/// cut logic ([`derive_is_kept`](Self::derive_is_kept) /
/// [`find_cut`](Self::find_cut) / [`chunk_cuts`](Self::chunk_cuts)) is copied
/// verbatim from `two_pass.rs` — it is the byte-identity core. The only change
/// is the fold's input: light-column slices from a
/// [`SamplePspChunk`](crate::var_calling_new::sample_reader::SamplePspChunk)
/// rather than `SampleColumns` accessors (the values are identical:
/// `ref_spans[i]` == `ref_span_at(i)`, `nonref_obs[i]` == `non_ref_obs_sum_at(i)`).
#[derive(Debug, Default)]
pub struct CohortSpanFold {
    /// Sorted, unique 1-based positions present in any folded sample.
    positions: Vec<u32>,
    /// Parallel to `positions`: max `ref_span` across samples with a record
    /// at the position (drives the grouping reach).
    max_ref_span: Vec<u32>,
    /// Parallel to `positions`: max over samples of that sample's summed
    /// non-REF observations at the position. `0` iff no sample carries an
    /// observed non-REF allele there.
    max_nonref_obs: Vec<u32>,
    // Double buffers for the merge (kept to avoid a per-call allocation).
    tmp_positions: Vec<u32>,
    tmp_max_ref_span: Vec<u32>,
    tmp_max_nonref_obs: Vec<u32>,
}

impl CohortSpanFold {
    /// A fresh, empty fold.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reset to empty, preserving allocated capacity for reuse across spans.
    pub fn clear(&mut self) {
        self.positions.clear();
        self.max_ref_span.clear();
        self.max_nonref_obs.clear();
    }

    /// Number of distinct positions folded so far.
    pub fn n_positions(&self) -> usize {
        self.positions.len()
    }

    /// The folded positions (sorted, unique).
    pub fn positions(&self) -> &[u32] {
        &self.positions
    }

    /// The per-position max ref span (parallel to [`positions`](Self::positions)).
    pub fn max_ref_span(&self) -> &[u32] {
        &self.max_ref_span
    }

    /// Merge one sample's light columns into the summary.
    ///
    /// `positions` must be position-sorted (the `SamplePspChunk` invariant)
    /// and the three slices parallel (one entry per record). Aggregation is
    /// `max` on ref-span and non-REF obs, with a position-union on the keys —
    /// order-independent, so folding the samples in any order yields the same
    /// summary. Mirrors `WindowSummary::fold_sample`.
    pub fn fold_sample_light(&mut self, positions: &[u32], ref_spans: &[u32], nonref_obs: &[u32]) {
        debug_assert_eq!(positions.len(), ref_spans.len());
        debug_assert_eq!(positions.len(), nonref_obs.len());
        let n = positions.len();
        self.tmp_positions.clear();
        self.tmp_max_ref_span.clear();
        self.tmp_max_nonref_obs.clear();

        let mut i = 0usize; // cursor into the current summary
        let mut j = 0usize; // cursor into the sample's records
        let m = self.positions.len();
        while i < m && j < n {
            let ps = self.positions[i];
            let pc = positions[j];
            if ps < pc {
                self.push_tmp(ps, self.max_ref_span[i], self.max_nonref_obs[i]);
                i += 1;
            } else if pc < ps {
                self.push_tmp(pc, ref_spans[j], nonref_obs[j]);
                j += 1;
            } else {
                self.push_tmp(
                    ps,
                    self.max_ref_span[i].max(ref_spans[j]),
                    self.max_nonref_obs[i].max(nonref_obs[j]),
                );
                i += 1;
                j += 1;
            }
        }
        while i < m {
            self.push_tmp(
                self.positions[i],
                self.max_ref_span[i],
                self.max_nonref_obs[i],
            );
            i += 1;
        }
        while j < n {
            self.push_tmp(positions[j], ref_spans[j], nonref_obs[j]);
            j += 1;
        }

        std::mem::swap(&mut self.positions, &mut self.tmp_positions);
        std::mem::swap(&mut self.max_ref_span, &mut self.tmp_max_ref_span);
        std::mem::swap(&mut self.max_nonref_obs, &mut self.tmp_max_nonref_obs);
    }

    #[inline]
    fn push_tmp(&mut self, pos: u32, max_ref_span: u32, max_nonref_obs: u32) {
        self.tmp_positions.push(pos);
        self.tmp_max_ref_span.push(max_ref_span);
        self.tmp_max_nonref_obs.push(max_nonref_obs);
    }

    /// Derive `is_kept` (parallel to [`positions`](Self::positions)): `true`
    /// for every position in a variant group whose summed `max_nonref_obs`
    /// reaches `min_alt_obs`. Groups form by the overlapping-reach rule and
    /// are kept or dropped whole. **Dust is *not* applied here** — it is
    /// folded in (masked positions still contribute obs+reach to the group
    /// decision) and dropped afterward by [`drop_dust_masked`], matching the
    /// old `is_kept` (pre-dust) → materialise → partition-skips-masked order.
    /// Copied verbatim from `WindowSummary::derive_is_kept`.
    pub fn derive_is_kept(&self, min_alt_obs: u32, out: &mut Vec<bool>) {
        let n = self.positions.len();
        out.clear();
        out.resize(n, false);
        let threshold = u64::from(min_alt_obs.max(1));

        let mut i = 0usize;
        while i < n {
            let mut group_end = reach(self.positions[i], self.max_ref_span[i]);
            let mut obs_sum = u64::from(self.max_nonref_obs[i]);
            let mut j = i + 1;
            while j < n && self.positions[j] <= group_end {
                group_end = group_end.max(reach(self.positions[j], self.max_ref_span[j]));
                obs_sum += u64::from(self.max_nonref_obs[j]);
                j += 1;
            }
            if obs_sum >= threshold {
                out[i..j].iter_mut().for_each(|slot| *slot = true);
            }
            i = j;
        }
    }

    /// Merge another fold into `self` (`max` aggregation, position union) —
    /// the reduce step for a parallel per-sample fold. Copied from
    /// `WindowSummary::merge`.
    pub fn merge(&mut self, other: &CohortSpanFold) {
        self.tmp_positions.clear();
        self.tmp_max_ref_span.clear();
        self.tmp_max_nonref_obs.clear();

        let (a_pos, a_mrs, a_mno) = (&self.positions, &self.max_ref_span, &self.max_nonref_obs);
        let (b_pos, b_mrs, b_mno) = (&other.positions, &other.max_ref_span, &other.max_nonref_obs);
        let (mut i, mut j) = (0usize, 0usize);
        let (m, n) = (a_pos.len(), b_pos.len());
        while i < m && j < n {
            if a_pos[i] < b_pos[j] {
                self.tmp_positions.push(a_pos[i]);
                self.tmp_max_ref_span.push(a_mrs[i]);
                self.tmp_max_nonref_obs.push(a_mno[i]);
                i += 1;
            } else if b_pos[j] < a_pos[i] {
                self.tmp_positions.push(b_pos[j]);
                self.tmp_max_ref_span.push(b_mrs[j]);
                self.tmp_max_nonref_obs.push(b_mno[j]);
                j += 1;
            } else {
                self.tmp_positions.push(a_pos[i]);
                self.tmp_max_ref_span.push(a_mrs[i].max(b_mrs[j]));
                self.tmp_max_nonref_obs.push(a_mno[i].max(b_mno[j]));
                i += 1;
                j += 1;
            }
        }
        self.tmp_positions.extend_from_slice(&a_pos[i..]);
        self.tmp_max_ref_span.extend_from_slice(&a_mrs[i..]);
        self.tmp_max_nonref_obs.extend_from_slice(&a_mno[i..]);
        self.tmp_positions.extend_from_slice(&b_pos[j..]);
        self.tmp_max_ref_span.extend_from_slice(&b_mrs[j..]);
        self.tmp_max_nonref_obs.extend_from_slice(&b_mno[j..]);

        std::mem::swap(&mut self.positions, &mut self.tmp_positions);
        std::mem::swap(&mut self.max_ref_span, &mut self.tmp_max_ref_span);
        std::mem::swap(&mut self.max_nonref_obs, &mut self.tmp_max_nonref_obs);
    }

    /// The clean group boundary at or before `watermark`: the start of the
    /// still-open group, or `watermark + 1` if every group is closed. Copied
    /// from `WindowSummary::find_cut` / `loader::find_block_cut`.
    pub fn find_cut(&self, watermark: u32) -> u32 {
        let n = self.positions.len();
        if n == 0 {
            return watermark.saturating_add(1);
        }
        let mut group_start = self.positions[0];
        let mut group_end = reach(self.positions[0], self.max_ref_span[0]);
        for k in 1..n {
            let pos = self.positions[k];
            if pos <= group_end {
                group_end = group_end.max(reach(pos, self.max_ref_span[k]));
            } else {
                group_start = pos;
                group_end = reach(pos, self.max_ref_span[k]);
            }
        }
        if group_end > watermark {
            group_start
        } else {
            watermark.saturating_add(1)
        }
    }

    /// Count kept (variable) positions with `position < cut`, by the same
    /// group/threshold rule as [`derive_is_kept`](Self::derive_is_kept).
    /// `cut` should land on a clean group boundary (e.g. from
    /// [`find_cut`](Self::find_cut)) so no group straddles it. Copied from
    /// `WindowSummary::count_kept_below`. Used by the read-ahead to size a
    /// chunk to ~`target_variants` kept positions.
    pub fn count_kept_below(&self, cut: u32, min_alt_obs: u32) -> u32 {
        let n = self.positions.len();
        let threshold = u64::from(min_alt_obs.max(1));
        let mut kept = 0u32;
        let mut i = 0usize;
        while i < n && self.positions[i] < cut {
            let mut group_end = reach(self.positions[i], self.max_ref_span[i]);
            let mut obs_sum = u64::from(self.max_nonref_obs[i]);
            let mut j = i + 1;
            while j < n && self.positions[j] <= group_end {
                group_end = group_end.max(reach(self.positions[j], self.max_ref_span[j]));
                obs_sum += u64::from(self.max_nonref_obs[j]);
                j += 1;
            }
            if obs_sum >= threshold {
                kept += (j - i) as u32;
            }
            i = j;
        }
        kept
    }

    /// Partition `[interval_start, interval_end_exclusive)` into chunk
    /// boundaries: each chunk accumulates ~`target_variants` kept positions
    /// and ends at the next variant group's start (a clean boundary). Copied
    /// from `WindowSummary::chunk_cuts`.
    pub fn chunk_cuts(
        &self,
        min_alt_obs: u32,
        target_variants: u32,
        interval_start: u32,
        interval_end_exclusive: u32,
    ) -> Vec<u32> {
        let n = self.positions.len();
        let threshold = u64::from(min_alt_obs.max(1));
        let target = target_variants.max(1);
        let mut cuts = vec![interval_start];
        let mut kept_in_chunk = 0u32;
        let mut i = 0usize;
        while i < n {
            let mut group_end = reach(self.positions[i], self.max_ref_span[i]);
            let mut obs_sum = u64::from(self.max_nonref_obs[i]);
            let mut j = i + 1;
            while j < n && self.positions[j] <= group_end {
                group_end = group_end.max(reach(self.positions[j], self.max_ref_span[j]));
                obs_sum += u64::from(self.max_nonref_obs[j]);
                j += 1;
            }
            if obs_sum >= threshold {
                kept_in_chunk += (j - i) as u32;
            }
            i = j;
            if kept_in_chunk >= target {
                let cut = if i < n {
                    self.positions[i]
                } else {
                    interval_end_exclusive
                };
                if cut > *cuts.last().expect("seeded with interval_start") {
                    cuts.push(cut);
                }
                kept_in_chunk = 0;
            }
        }
        if *cuts.last().expect("seeded with interval_start") < interval_end_exclusive {
            cuts.push(interval_end_exclusive);
        }
        cuts
    }
}

/// Apply the dust mask to an already-derived `is_kept` (parallel to
/// `positions`): clear every position that falls inside a masked interval.
///
/// `mask` is sorted, non-overlapping, **half-open** `[start, end)` genomic
/// intervals (the chunk-level DUST mask in genomic coordinates — exactly the
/// shape `partition_window` consumes). Run *after*
/// [`CohortSpanFold::derive_is_kept`], so masked positions still counted
/// toward their group's keep decision but are not themselves emitted; the
/// caller's grouper bridges the gaps they leave (their reach was never used —
/// see `partition.rs` §DUST). Byte-identical to the old skip-in-partition.
pub fn drop_dust_masked(positions: &[u32], is_kept: &mut [bool], mask: &[Range<u32>]) {
    debug_assert_eq!(positions.len(), is_kept.len());
    if mask.is_empty() {
        return;
    }
    // Merge-walk: positions ascending, mask intervals ascending.
    let mut m = 0usize;
    for (p_idx, &pos) in positions.iter().enumerate() {
        // Advance past mask intervals that end at or before `pos`.
        while m < mask.len() && mask[m].end <= pos {
            m += 1;
        }
        if m >= mask.len() {
            break;
        }
        if mask[m].start <= pos {
            // `pos` ∈ [start, end) — masked.
            is_kept[p_idx] = false;
        }
    }
}

// ===========================================================================
// 2b — the streaming integrator (single covered interval)
// ===========================================================================

use crate::pileup_record::PileupRecord;
use crate::psp::PspReadError;
use crate::var_calling_new::per_position_merger::{PerPositionMerger, PerPositionMergerError};
use crate::var_calling_new::sample_reader::{SamplePspChunk, SamplePspReader};
use crate::var_calling_new::types::{CohortPileupRecord, PileupCohortChunk, RefSpan};
use std::io::{Read, Seek};

/// Errors the cohort producer can surface.
#[derive(Debug, thiserror::Error)]
pub enum ProducerError {
    /// A per-sample `.psp` decode failed.
    #[error("psp decode failed: {0}")]
    Decode(#[from] PspReadError),
    /// The per-position cohort merge failed (pathological input only — the
    /// in-memory kept-record streams are monotone by construction).
    #[error("per-position merge failed: {0}")]
    Merge(Box<PerPositionMergerError>),
    /// The reference-fetch closure failed.
    #[error("reference fetch failed: {0}")]
    Ref(String),
}

/// Wrap an owned [`PileupRecord`] as the `Result` item the
/// [`PerPositionMerger`] consumes. A named `fn` (not a closure) so every
/// per-sample iterator has the *same* type and they collect into a `Vec<I>`.
fn ok_record(r: PileupRecord) -> Result<PileupRecord, PspReadError> {
    Ok(r)
}

type KeptRecordIter = std::iter::Map<
    std::vec::IntoIter<PileupRecord>,
    fn(PileupRecord) -> Result<PileupRecord, PspReadError>,
>;

/// Streaming cohort producer over **one covered interval** (appendix §B).
///
/// Owns the N [`SamplePspReader`]s. [`begin_interval`](Self::begin_interval)
/// points them at a `[start, end)` interval (+ its dust mask);
/// [`produce_chunk`](Self::produce_chunk) then yields one
/// [`PileupCohortChunk`] of variable [`CohortPileupRecord`]s at a time, cut at
/// safe gaps, until the interval is drained (`Ok(None)`). The chromosome /
/// interval iteration, per-chromosome REF fetcher, and dust computation are
/// the wiring layer (Phase 2b-wiring); `produce_chunk` takes the REF fetch as
/// a closure so it stays decoupled and unit-testable.
///
/// Per chunk: lockstep-read segments to the watermark → fold the **light**
/// columns → `derive_is_kept` → drop dust → cut at a safe gap → build records
/// (columnar→record via [`SamplePspChunk::records_for`], merged per-position
/// by the revived [`PerPositionMerger`]) → fetch the REF span → stamp
/// `chunk_order`. Only the variable positions are ever built (the memory
/// invariant); buffered columnar segments straddling a cut survive to the next
/// chunk (records_for is non-consuming).
pub struct CohortChunkIntegrator<R: Read + Seek> {
    readers: Vec<SamplePspReader<R>>,
    sample_names: Vec<String>,
    n: usize,
    min_alt_obs: u32,
    target_variants: u32,

    // Current-interval state.
    chrom_id: u32,
    interval_end: u32, // exclusive
    next_chunk_start: u32,
    /// Per sample: buffered columnar segments (carryover + read-ahead), in
    /// genomic order; dropped once fully below `next_chunk_start`.
    buffers: Vec<Vec<SamplePspChunk>>,
    /// Per sample: reader exhausted within the current interval.
    exhausted: Vec<bool>,
    dust_mask: Vec<Range<u32>>,

    /// Monotonic genomic-order chunk counter (the writer's reorder ticket).
    /// **Not** reset per interval.
    chunk_order: u64,

    // Reusable scratch.
    fold: CohortSpanFold,
    is_kept: Vec<bool>,
}

impl<R: Read + Seek> CohortChunkIntegrator<R> {
    /// Build a producer over the N per-sample readers (`sample_names`
    /// parallel to `readers`). `min_alt_obs` is the cohort keep threshold;
    /// `target_variants` the soft per-chunk kept-position target (`0` ⇒ 1).
    pub fn new(
        readers: Vec<SamplePspReader<R>>,
        sample_names: Vec<String>,
        min_alt_obs: u32,
        target_variants: u32,
    ) -> Self {
        let n = readers.len();
        Self {
            readers,
            sample_names,
            n,
            min_alt_obs,
            target_variants: target_variants.max(1),
            chrom_id: 0,
            interval_end: 0,
            next_chunk_start: 0,
            buffers: vec![Vec::new(); n],
            exhausted: vec![false; n],
            dust_mask: Vec::new(),
            chunk_order: 0,
            fold: CohortSpanFold::new(),
            is_kept: Vec::new(),
        }
    }

    /// Point every reader at `[interval.start, interval.end)` of `chrom_id`
    /// (reusing decode buffers) and install the interval's dust mask. Clears
    /// the per-interval buffers; preserves `chunk_order`.
    pub fn begin_interval(
        &mut self,
        chrom_id: u32,
        interval: Range<u32>,
        dust_mask: Vec<Range<u32>>,
    ) {
        debug_assert!(interval.start < interval.end, "empty interval");
        let region_end = interval.end - 1; // reader region_end is inclusive
        for r in &mut self.readers {
            r.reset(chrom_id, interval.start, region_end);
        }
        self.chrom_id = chrom_id;
        self.interval_end = interval.end;
        self.next_chunk_start = interval.start;
        for b in &mut self.buffers {
            b.clear();
        }
        for e in &mut self.exhausted {
            *e = false;
        }
        self.dust_mask = dust_mask;
    }

    /// The last position every sample's buffered data reaches (`interval_end`
    /// for exhausted / empty-buffer samples) — `min` over samples is the
    /// cohort watermark: positions `≤ W` have complete cohort data.
    fn coverage(&self, s: usize) -> u32 {
        if self.exhausted[s] {
            self.interval_end
        } else {
            self.buffers[s]
                .last()
                .and_then(|c| c.positions().last().copied())
                .unwrap_or(0)
        }
    }

    fn watermark(&self) -> u32 {
        (0..self.n)
            .map(|s| self.coverage(s))
            .min()
            .unwrap_or(self.interval_end)
    }

    fn all_exhausted(&self) -> bool {
        self.exhausted.iter().all(|&e| e)
    }

    /// Read one more segment for sample `s` (or mark it exhausted).
    fn read_one(&mut self, s: usize) -> Result<(), ProducerError> {
        match self.readers[s].next_chunk()? {
            Some(c) => self.buffers[s].push(c),
            None => self.exhausted[s] = true,
        }
        Ok(())
    }

    /// Rebuild [`self.fold`] over the buffered light columns in
    /// `[next_chunk_start, w]` (all samples' data is complete there).
    fn rebuild_fold(&mut self, w: u32) {
        let mut fold = std::mem::take(&mut self.fold);
        fold.clear();
        let lo_bound = self.next_chunk_start;
        for buf in &self.buffers {
            for chunk in buf {
                let pos = chunk.positions();
                let lo = pos.partition_point(|&p| p < lo_bound);
                let hi = pos.partition_point(|&p| p <= w);
                if lo < hi {
                    fold.fold_sample_light(
                        &pos[lo..hi],
                        &chunk.ref_spans()[lo..hi],
                        &chunk.nonref_obs()[lo..hi],
                    );
                }
            }
        }
        self.fold = fold;
    }

    /// Read-ahead: pull segments (advancing the laggard samples) until the
    /// chunk has ~`target_variants` kept positions *and* a safe gap to cut at,
    /// or the interval is fully read. Leaves [`self.fold`] valid over
    /// `[next_chunk_start, watermark]`.
    fn fill_to_target(&mut self) -> Result<(), ProducerError> {
        loop {
            // Seed any non-exhausted sample that has no buffered segment.
            for s in 0..self.n {
                if !self.exhausted[s] && self.buffers[s].is_empty() {
                    self.read_one(s)?;
                }
            }
            let w = self.watermark();
            self.rebuild_fold(w);
            let kept = self
                .fold
                .count_kept_below(w.saturating_add(1), self.min_alt_obs);
            let has_safe_cut =
                self.all_exhausted() || self.fold.find_cut(w) > self.next_chunk_start;
            if (kept >= self.target_variants && has_safe_cut)
                || self.all_exhausted()
                || w >= self.interval_end
            {
                return Ok(());
            }
            // Advance the laggards (coverage == W) to raise the watermark.
            let mut advanced = false;
            for s in 0..self.n {
                if !self.exhausted[s] && self.coverage(s) == w {
                    self.read_one(s)?;
                    advanced = true;
                }
            }
            if !advanced {
                return Ok(());
            }
        }
    }

    /// Produce the next chunk of the current interval, or `Ok(None)` when the
    /// interval is drained. `fetch_ref` returns the REF bases for a 1-based
    /// `(start, length)` span (errors stringified into [`ProducerError::Ref`]).
    pub fn produce_chunk(
        &mut self,
        fetch_ref: &mut dyn FnMut(u32, u32) -> Result<Vec<u8>, String>,
    ) -> Result<Option<PileupCohortChunk>, ProducerError> {
        loop {
            if self.next_chunk_start >= self.interval_end {
                return Ok(None);
            }

            self.fill_to_target()?;
            let w = self.watermark();
            let cut = if self.all_exhausted() {
                self.interval_end
            } else {
                self.fold.find_cut(w).min(self.interval_end)
            };
            debug_assert!(
                cut > self.next_chunk_start,
                "chunk cut must make progress ({} <= {})",
                cut,
                self.next_chunk_start,
            );

            // Keep decision over the (complete) fold, then drop dust.
            self.fold
                .derive_is_kept(self.min_alt_obs, &mut self.is_kept);
            drop_dust_masked(self.fold.positions(), &mut self.is_kept, &self.dust_mask);
            let variable: Vec<u32> = self
                .fold
                .positions()
                .iter()
                .zip(self.is_kept.iter())
                .filter(|&(&p, &k)| k && p < cut)
                .map(|(&p, _)| p)
                .collect();

            let chunk_start = self.next_chunk_start;

            if variable.is_empty() {
                // Variant-free span — emit nothing, keep chunk_order gapless.
                self.next_chunk_start = cut;
                self.drop_buffers_below(cut);
                continue;
            }

            // Build from the buffers *before* advancing past them.
            let records = self.build_records(&variable, chunk_start, cut)?;
            let ref_span = self.fetch_ref_span(&variable, fetch_ref)?;
            self.next_chunk_start = cut;
            self.drop_buffers_below(cut);

            let chunk_order = self.chunk_order;
            self.chunk_order += 1;
            return Ok(Some(PileupCohortChunk {
                chunk_order,
                records,
                ref_span,
            }));
        }
    }

    /// Drop buffered segments wholly below `cut` (fully consumed).
    fn drop_buffers_below(&mut self, cut: u32) {
        for buf in &mut self.buffers {
            buf.retain(|c| c.positions().last().is_some_and(|&p| p >= cut));
        }
    }

    /// Build the chunk's [`CohortPileupRecord`]s: per sample, reconstruct its
    /// records on the variable positions in `[chunk_start, cut)` (columnar→
    /// record), then merge per-position across samples via [`PerPositionMerger`].
    fn build_records(
        &self,
        variable: &[u32],
        chunk_start: u32,
        cut: u32,
    ) -> Result<Vec<CohortPileupRecord>, ProducerError> {
        let mut iters: Vec<KeptRecordIter> = Vec::with_capacity(self.n);
        for buf in &self.buffers {
            let mut recs: Vec<PileupRecord> = Vec::new();
            for chunk in buf {
                let pos = chunk.positions();
                let keep: Vec<bool> = pos
                    .iter()
                    .map(|&p| p >= chunk_start && p < cut && variable.binary_search(&p).is_ok())
                    .collect();
                if keep.iter().any(|&k| k) {
                    recs.extend(chunk.records_for(&keep));
                }
            }
            iters.push(
                recs.into_iter()
                    .map(ok_record as fn(PileupRecord) -> Result<PileupRecord, PspReadError>),
            );
        }

        let merger = PerPositionMerger::new(iters, self.sample_names.clone(), Vec::new())
            .map_err(|e| ProducerError::Merge(Box::new(e)))?;
        let mut records = Vec::with_capacity(variable.len());
        for item in merger {
            let pp = item.map_err(|e| ProducerError::Merge(Box::new(e)))?;
            records.push(CohortPileupRecord {
                chrom_id: pp.chrom_id,
                pos: pp.pos,
                per_sample: pp.per_sample,
            });
        }
        Ok(records)
    }

    /// Fetch the chunk's one contiguous REF span `[first_kept, max_reach]`
    /// (monotonic-forward across chunks). `max_reach` uses the cohort
    /// `max_ref_span` (the grouping reach), matching the old per-group fetch.
    fn fetch_ref_span(
        &self,
        variable: &[u32],
        fetch_ref: &mut dyn FnMut(u32, u32) -> Result<Vec<u8>, String>,
    ) -> Result<RefSpan, ProducerError> {
        let first = variable[0];
        let positions = self.fold.positions();
        let spans = self.fold.max_ref_span();
        let mut max_reach = first;
        for &p in variable {
            let idx = positions
                .binary_search(&p)
                .expect("variable position must be in the fold");
            max_reach = max_reach.max(reach(p, spans[idx]));
        }
        let len = max_reach - first + 1;
        let bytes = fetch_ref(first, len).map_err(ProducerError::Ref)?;
        Ok(RefSpan {
            genomic_start: first,
            bytes,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Fold a set of per-sample light columns `(positions, ref_spans,
    /// nonref_obs)` into a fresh [`CohortSpanFold`].
    fn fold(samples: &[(Vec<u32>, Vec<u32>, Vec<u32>)]) -> CohortSpanFold {
        let mut s = CohortSpanFold::new();
        for (pos, rs, nro) in samples {
            s.fold_sample_light(pos, rs, nro);
        }
        s
    }

    #[test]
    fn fold_unions_positions_and_maxes_aggregates() {
        // A: pos10 (REF only, nonref 0), pos20 (nonref 1).
        let a = (vec![10, 20], vec![1, 1], vec![0, 1]);
        // B: pos20 (nonref 2), pos30 (REF only).
        let b = (vec![20, 30], vec![1, 1], vec![2, 0]);
        let s = fold(&[a.clone(), b.clone()]);
        assert_eq!(s.positions(), &[10, 20, 30]);
        assert_eq!(s.max_nonref_obs, vec![0, 2, 0]);
        // Order-independent.
        let s2 = fold(&[b, a]);
        assert_eq!(s2.positions(), s.positions());
        assert_eq!(s2.max_nonref_obs, s.max_nonref_obs);
        assert_eq!(s2.max_ref_span, s.max_ref_span);
    }

    #[test]
    fn keep_threshold_one_is_variant_filter() {
        // Three isolated SNPs (ref_span 1): pos10 obs0, pos20 obs1, pos30 obs0.
        let s = fold(&[(vec![10, 20, 30], vec![1, 1, 1], vec![0, 1, 0])]);
        let mut kept = Vec::new();
        s.derive_is_kept(1, &mut kept);
        assert_eq!(kept, vec![false, true, false]);
    }

    #[test]
    fn keep_threshold_two_drops_singletons_keeps_doubletons() {
        let s = fold(&[(vec![20, 30], vec![1, 1], vec![1, 2])]);
        let mut kept = Vec::new();
        s.derive_is_kept(2, &mut kept);
        assert_eq!(kept, vec![false, true]);
    }

    #[test]
    fn over_approximation_is_max_not_sum_across_samples() {
        // Two samples each one ALT obs at the same position: max = 1 < 2.
        let s = fold(&[(vec![20], vec![1], vec![1]), (vec![20], vec![1], vec![1])]);
        let mut kept = Vec::new();
        s.derive_is_kept(2, &mut kept);
        assert_eq!(kept, vec![false]);
    }

    #[test]
    fn multi_position_group_kept_whole() {
        // An MNP at 10 with ref_span 3 reaches 12; the SNP at 12 (obs 0)
        // joins its group and is kept because the group has a variant.
        let s = fold(&[(vec![10, 12], vec![3, 1], vec![2, 0])]);
        let mut kept = Vec::new();
        s.derive_is_kept(1, &mut kept);
        assert_eq!(
            kept,
            vec![true, true],
            "both positions of the variant group kept"
        );
    }

    #[test]
    fn merge_matches_sequential_fold() {
        let a = (vec![10, 20], vec![1, 1], vec![0, 1]);
        let b = (vec![20, 30], vec![1, 1], vec![2, 0]);
        let seq = fold(&[a.clone(), b.clone()]);
        let mut pa = fold(&[a]);
        let pb = fold(&[b]);
        pa.merge(&pb);
        assert_eq!(pa.positions(), seq.positions());
        assert_eq!(pa.max_ref_span, seq.max_ref_span);
        assert_eq!(pa.max_nonref_obs, seq.max_nonref_obs);
    }

    #[test]
    fn find_cut_keeps_open_group_intact() {
        // Variant SNPs at 10 and 100, non-variant at 50. Watermark 60: the
        // group at 100 is still open (reach 100 > 60) ⇒ cut at its start.
        let s = fold(&[(vec![10, 50, 100], vec![1, 1, 1], vec![2, 0, 2])]);
        assert_eq!(s.find_cut(60), 100);
        // Watermark past everything ⇒ all groups closed.
        assert_eq!(s.find_cut(200), 201);
    }

    #[test]
    fn drop_dust_masked_clears_masked_positions() {
        // positions 10, 20, 30; mask [15,25) covers 20 only.
        let positions = [10u32, 20, 30];
        let mut is_kept = vec![true, true, true];
        let r = 15u32..25;
        drop_dust_masked(&positions, &mut is_kept, &[r]);
        assert_eq!(is_kept, vec![true, false, true]);
    }

    #[test]
    fn drop_dust_masked_half_open_boundaries() {
        // [20,30) masks 20..=29; 30 is NOT masked (half-open end).
        let positions = [19u32, 20, 29, 30];
        let mut is_kept = vec![true, true, true, true];
        let r = 20u32..30;
        drop_dust_masked(&positions, &mut is_kept, &[r]);
        assert_eq!(is_kept, vec![true, false, false, true]);
    }

    #[test]
    fn drop_dust_masked_multiple_intervals() {
        let positions = [5u32, 10, 15, 20, 25];
        let mut is_kept = vec![true; 5];
        drop_dust_masked(&positions, &mut is_kept, &[8..12, 22..30]);
        // 10 masked by [8,12); 25 masked by [22,30).
        assert_eq!(is_kept, vec![true, false, true, true, false]);
    }

    // --- 2b: the streaming integrator (cross-checked against a one-shot ref) ---

    use crate::pileup_record::PileupRecord;
    use crate::psp::PspReader;
    use crate::psp::test_fixtures::writer_header;
    use crate::psp::writer::PspWriter;
    use crate::var_calling::test_helpers::{allele, record};
    use crate::var_calling_new::sample_reader::SamplePspReader;
    use crate::var_calling_new::types::CohortPileupRecord;
    use std::io::Cursor;

    const INTERVAL_END: u32 = 4000;

    /// One synthetic sample: a misaligned, multi-allele record stream with a
    /// mix of variant / non-variant positions and the occasional 2-base REF
    /// (MNP reach). `s` shifts coverage + obs so samples disagree per position.
    fn gen_sample(s: u32, n: u32) -> Vec<PileupRecord> {
        (0..n)
            .filter_map(|i| {
                // Per-sample coverage gaps → misaligned position sets.
                if (i + s).is_multiple_of(4) {
                    return None;
                }
                let pos = 10 + i * 5;
                let ref_seq: &[u8] = if i.is_multiple_of(9) { b"AC" } else { b"A" };
                let mut alleles = vec![allele(ref_seq, 5, -1.0, &[])];
                let alt_obs = if (i * 7 + s * 3).is_multiple_of(5) {
                    2
                } else {
                    0
                };
                if alt_obs > 0 {
                    alleles.push(allele(b"T", alt_obs, -1.5, &[u64::from(i) + 1]));
                }
                Some(record(pos, alleles))
            })
            .collect()
    }

    fn psp_bytes(records: &[PileupRecord], block_target: usize) -> Vec<u8> {
        let mut w = PspWriter::new_with_block_target(
            Cursor::new(Vec::new()),
            writer_header(1),
            block_target,
        )
        .expect("writer opens");
        for r in records {
            w.write_record(r).expect("write_record");
        }
        w.finish().expect("finish").into_inner()
    }

    fn names(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("S{i}")).collect()
    }

    /// One-shot reference: fold all samples fully, derive kept positions, drop
    /// dust, assemble per-position cohort records — no streaming / chunking.
    fn reference(
        samples: &[Vec<PileupRecord>],
        min_alt_obs: u32,
        mask: &[Range<u32>],
    ) -> Vec<CohortPileupRecord> {
        let mut fold = CohortSpanFold::new();
        for recs in samples {
            let pos: Vec<u32> = recs.iter().map(|r| r.pos).collect();
            let rs: Vec<u32> = recs.iter().map(|r| r.alleles[0].seq.len() as u32).collect();
            let nro: Vec<u32> = recs
                .iter()
                .map(|r| {
                    r.alleles[1..]
                        .iter()
                        .fold(0u32, |a, al| a.saturating_add(al.support.num_obs))
                })
                .collect();
            fold.fold_sample_light(&pos, &rs, &nro);
        }
        let mut is_kept = Vec::new();
        fold.derive_is_kept(min_alt_obs, &mut is_kept);
        drop_dust_masked(fold.positions(), &mut is_kept, mask);
        let variable: Vec<u32> = fold
            .positions()
            .iter()
            .zip(&is_kept)
            .filter(|&(_, &k)| k)
            .map(|(&p, _)| p)
            .collect();

        let iters: Vec<KeptRecordIter> = samples
            .iter()
            .map(|recs| {
                let kept: Vec<PileupRecord> = recs
                    .iter()
                    .filter(|r| variable.binary_search(&r.pos).is_ok())
                    .cloned()
                    .collect();
                kept.into_iter()
                    .map(ok_record as fn(PileupRecord) -> Result<PileupRecord, PspReadError>)
            })
            .collect();
        let merger = PerPositionMerger::new(iters, names(samples.len()), Vec::new()).unwrap();
        merger
            .map(|x| {
                let pp = x.unwrap();
                CohortPileupRecord {
                    chrom_id: pp.chrom_id,
                    pos: pp.pos,
                    per_sample: pp.per_sample,
                }
            })
            .collect()
    }

    /// Streaming output: drain `produce_chunk` over one interval, asserting
    /// `chunk_order` is gapless and records stay genomically sorted.
    fn streaming(
        samples: &[Vec<PileupRecord>],
        block_targets: &[usize],
        min_alt_obs: u32,
        target_variants: u32,
        mask: Vec<Range<u32>>,
    ) -> Vec<CohortPileupRecord> {
        let readers: Vec<_> = samples
            .iter()
            .zip(block_targets)
            .map(|(recs, &bt)| {
                let bytes = psp_bytes(recs, bt);
                let reader = PspReader::new(Cursor::new(bytes)).unwrap();
                SamplePspReader::new(reader, 0, 1, 1)
            })
            .collect();
        let mut integ =
            CohortChunkIntegrator::new(readers, names(samples.len()), min_alt_obs, target_variants);
        integ.begin_interval(0, 1..INTERVAL_END, mask);

        let mut fetch = |_start: u32, len: u32| Ok::<_, String>(vec![b'N'; len as usize]);
        let mut out: Vec<CohortPileupRecord> = Vec::new();
        let mut expected_order = 0u64;
        while let Some(chunk) = integ.produce_chunk(&mut fetch).unwrap() {
            assert_eq!(chunk.chunk_order, expected_order, "chunk_order gapless");
            expected_order += 1;
            assert!(!chunk.records.is_empty(), "empty chunks are never emitted");
            // RefSpan covers the chunk's first kept position.
            assert_eq!(chunk.ref_span.genomic_start, chunk.records[0].pos);
            out.extend(chunk.records);
        }
        // Globally sorted, strictly ascending positions.
        assert!(
            out.windows(2).all(|w| w[0].pos < w[1].pos),
            "records strictly ascending across chunks",
        );
        out
    }

    fn cohort() -> (Vec<Vec<PileupRecord>>, Vec<usize>) {
        let samples = vec![gen_sample(0, 300), gen_sample(1, 300), gen_sample(2, 300)];
        let block_targets = vec![200usize, 350, 512]; // misaligned blocks
        (samples, block_targets)
    }

    #[test]
    fn streaming_matches_reference_across_chunk_sizes() {
        let (samples, bts) = cohort();
        for &min_alt_obs in &[1u32, 2] {
            let want = reference(&samples, min_alt_obs, &[]);
            assert!(!want.is_empty(), "fixture should yield variants");
            // Same output regardless of how the interval is partitioned.
            for &target in &[1u32, 3, 17, 100_000] {
                let got = streaming(&samples, &bts, min_alt_obs, target, Vec::new());
                assert_eq!(got, want, "min_alt_obs={min_alt_obs} target={target}");
            }
        }
    }

    #[test]
    fn streaming_applies_dust_mask() {
        let (samples, bts) = cohort();
        let min_alt_obs = 1;
        // Mask a chunk of the genome; reference and streaming must agree, and
        // the masked span must actually drop some records.
        let masked_span = 500u32..900;
        let mask = vec![masked_span];
        let want = reference(&samples, min_alt_obs, &mask);
        let unmasked = reference(&samples, min_alt_obs, &[]);
        assert!(want.len() < unmasked.len(), "mask should drop records");
        for &target in &[1u32, 5, 100_000] {
            let got = streaming(&samples, &bts, min_alt_obs, target, mask.clone());
            assert_eq!(got, want, "dust target={target}");
        }
        assert!(
            want.iter().all(|r| !(500..900).contains(&r.pos)),
            "no masked position survives",
        );
    }

    #[test]
    fn interval_with_no_variants_yields_no_chunks() {
        // A cohort with zero ALT obs anywhere → no variant positions.
        let flat: Vec<Vec<PileupRecord>> = (0..2)
            .map(|_| {
                (0..50)
                    .map(|i| record(10 + i * 5, vec![allele(b"A", 4, -1.0, &[])]))
                    .collect()
            })
            .collect();
        let got = streaming(&flat, &[128, 128], 1, 4, Vec::new());
        assert!(got.is_empty(), "no variants ⇒ no chunks");
    }
}
