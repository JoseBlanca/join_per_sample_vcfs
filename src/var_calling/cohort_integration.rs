//! Cohort producer — section 1 (appendix §B).
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
//! [`per_position_merger`](crate::var_calling::per_position_merger), here
//! in the producer.
//!
//! Two parts:
//! - [`CohortSpanFold`]: the cohort keep/cut math, folding over the reader's
//!   *light* columns ([`SamplePspChunk`]
//!   `positions`/`ref_spans`/`nonref_obs`). Plus [`drop_dust_masked`], the dust
//!   step. Byte-identity-critical: the arithmetic is byte-for-byte the
//!   pre-rewrite cohort keep/cut (verified out-of-tree).
//! - [`CohortChunkIntegrator`]: the streaming integrator over **one covered
//!   interval** — segment buffering to the watermark, safe-gap cut, REF fetch
//!   (closure-injected), `chunk_order` stamping. `produce_chunk` takes the REF
//!   fetch as a closure and the dust mask per
//!   [`begin_interval`](CohortChunkIntegrator::begin_interval) so the engine
//!   stays decoupled and unit-testable. The columns→records rebuild + the
//!   per-position merge run on the **caller**
//!   ([`merge_compacted_samples`](crate::var_calling::variant_caller::merge_compacted_samples)).

use std::ops::Range;

/// Reach of a position given its (max) ref span — `pos + max(span, 1) - 1`,
/// saturating. Byte-for-byte the pre-rewrite grouping `reach` (the grouping
/// arithmetic must match exactly; verified out-of-tree).
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
/// [`find_cut`](Self::find_cut)) is copied
/// verbatim from `two_pass.rs` — it is the byte-identity core. The only change
/// is the fold's input: light-column slices from a
/// [`SamplePspChunk`]
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
    /// Byte-for-byte the pre-rewrite cohort keep rule (verified out-of-tree).
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
    /// the reduce step for a parallel per-sample fold (associative +
    /// commutative, so the reduce is order-independent).
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
    /// still-open group, or `watermark + 1` if every group is closed.
    /// Byte-for-byte the pre-rewrite safe-gap cut (verified out-of-tree).
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
    /// [`find_cut`](Self::find_cut)) so no group straddles it. Used by the
    /// read-ahead to size a chunk to ~`target_variants` kept positions.
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

/// Membership mask `out[i] = needles[i] ∈ haystack`, for two ascending,
/// duplicate-free `u32` position slices, via one two-pointer merge-walk —
/// O(|needles| + |haystack|) instead of one binary search per needle (L7).
/// Both inputs are position-sorted by the `SamplePspChunk` / fold invariants;
/// byte-identical to the per-row `haystack.binary_search(&p).is_ok()` it
/// replaces.
fn membership_mask(needles: &[u32], haystack: &[u32]) -> Vec<bool> {
    let mut out = Vec::with_capacity(needles.len());
    let mut j = 0usize;
    for &p in needles {
        while j < haystack.len() && haystack[j] < p {
            j += 1;
        }
        out.push(j < haystack.len() && haystack[j] == p);
    }
    out
}

// ===========================================================================
// 2b — the streaming integrator (single covered interval)
// ===========================================================================

use crate::psp::PspReadError;
use crate::psp::block::new_column_decompressor;
use crate::var_calling::sample_reader::{SamplePspChunk, SamplePspReader, TwoPhaseSegment};
use crate::var_calling::types::{RawCohortChunk, RefSpan};
use std::io::{Read, Seek};
use std::sync::Arc;

/// Boxed error type the producer's REF-fetch closure yields. Boxing keeps
/// the producer decoupled from the concrete fetcher while preserving the
/// typed cause through [`ProducerError::Ref`]'s `source()` chain (the
/// pipeline's closure boxes a [`ChromRefFetchError`](crate::fasta::ChromRefFetchError)
/// directly, no stringification).
pub type RefFetchError = Box<dyn std::error::Error + Send + Sync>;

/// Errors the cohort producer can surface.
#[derive(Debug, thiserror::Error)]
pub enum ProducerError {
    /// A per-sample `.psp` decode failed.
    #[error("psp decode failed: {0}")]
    Decode(#[from] PspReadError),
    /// The reference-fetch closure failed; the fetcher's typed cause is
    /// preserved through `source()`.
    #[error("reference fetch failed")]
    Ref(#[source] RefFetchError),
    /// The chunk cut failed to advance past `next_chunk_start`. A degenerate
    /// fold (every reach collapsing at the boundary) would otherwise spin
    /// `produce_chunk` forever; this is the release-level guard for the
    /// progress invariant (previously a `debug_assert!`).
    #[error("chunk cut {cut} did not advance past next_chunk_start {next_chunk_start}")]
    StalledCut { next_chunk_start: u32, cut: u32 },
}

/// Merge per-block `(first_pos, last_pos)` ranges (across all samples, any
/// order) into the chromosome's **covered intervals** `[start, end)`: two
/// blocks join unless the gap between them exceeds `max_group_span`, so every
/// interval boundary is a safe gap no variant group can span — the producer
/// processes each interval independently.
fn merge_block_ranges(
    ranges: impl IntoIterator<Item = (u32, u32)>,
    max_group_span: u32,
) -> Vec<Range<u32>> {
    let mut v: Vec<(u32, u32)> = ranges.into_iter().collect();
    if v.is_empty() {
        return Vec::new();
    }
    v.sort_unstable();
    let mut out: Vec<Range<u32>> = Vec::new();
    let (mut cur_start, mut cur_last) = v[0];
    for &(s, last) in &v[1..] {
        if s <= cur_last.saturating_add(max_group_span) {
            cur_last = cur_last.max(last);
        } else {
            out.push(cur_start..cur_last.saturating_add(1));
            cur_start = s;
            cur_last = last;
        }
    }
    out.push(cur_start..cur_last.saturating_add(1));
    out
}

/// One chromosome's covered intervals for a slice of readers — every sample's
/// blocks for `chrom_id`, unioned and gap-merged at `max_group_span` (no file
/// I/O, reads only the in-memory block indexes).
pub(crate) fn covered_intervals_for<R: Read + Seek>(
    readers: &[SamplePspReader<R>],
    chrom_id: u32,
    max_group_span: u32,
) -> Vec<Range<u32>> {
    let ranges = readers.iter().flat_map(|r| {
        r.block_index()
            .iter()
            .filter(move |b| b.chrom_id == chrom_id)
            .map(|b| (b.first_pos, b.last_pos))
    });
    merge_block_ranges(ranges, max_group_span)
}

/// A planned-but-not-yet-compacted cohort chunk — the hand-off between the
/// producer's two internal stages (appendix §B, "fold/plan" → "compact").
///
/// The **fold/plan** stage ([`plan_chunk`](CohortChunkIntegrator::plan_chunk))
/// runs the cohort fold over the light columns, derives the variable positions
/// and the safe-gap cut, fetches the chunk's REF span, and captures — per
/// sample — the (still-compressed) [`TwoPhaseSegment`]s overlapping
/// `[chunk_start, cut)` as shared `Arc`s. The **compact** stage
/// ([`compact_plan`]) then inflates only those segments' heavy columns for the
/// variable rows, off the fold stage's critical path. A bounded queue between
/// the two lets the next plan's fold overlap this plan's heavy decode; the
/// `Arc` sharing keeps a segment that straddles two chunks resident once (no
/// copy), so peak memory tracks the queue depth, not duplication.
pub(crate) struct ChunkPlan {
    chunk_order: u64,
    chrom_id: u32,
    chunk_start: u32,
    cut: u32,
    /// Every variable (post-dust) position over the whole current fold, sorted.
    kept_all: Vec<u32>,
    ref_span: RefSpan,
    /// Per sample: the buffered segments overlapping `[chunk_start, cut)`,
    /// shared with the producer's live buffer (a straddler is referenced by two
    /// consecutive plans).
    per_sample_segments: Vec<Vec<Arc<TwoPhaseSegment>>>,
}

/// Compact stage: inflate each plan'd segment's heavy columns for its variable
/// rows in `[chunk_start, cut)` and assemble the cohort chunk. Stateless and
/// `self`-free, so it runs on a separate stage thread (still fanning out across
/// the N samples on the producer rayon pool). Byte-identical regardless of
/// `cache_straddlers`: same per-segment inflate restricted to
/// `kept_all ∩ [chunk_start, cut)`, same genomic row order.
///
/// `cache_straddlers` (the default; `false` under `--low-memory`) decodes a
/// cut-spanning segment's heavy columns once into its cache and slices both
/// overlapping chunks from it. With it `false`, every chunk re-decompresses the
/// straddler instead — a few percent slower for a peak-RSS reduction that scales
/// with the cohort size.
pub(crate) fn compact_plan(
    plan: ChunkPlan,
    cache_straddlers: bool,
) -> Result<RawCohortChunk, ProducerError> {
    use rayon::prelude::*;
    let ChunkPlan {
        chunk_order,
        chrom_id,
        chunk_start,
        cut,
        kept_all,
        ref_span,
        per_sample_segments,
    } = plan;
    let per_sample: Vec<SamplePspChunk> = per_sample_segments
        .par_iter()
        .map(|segs| -> Result<SamplePspChunk, ProducerError> {
            // Per-task zstd context + scratch, created lazily.
            let mut decompressor: Option<zstd::bulk::Decompressor<'static>> = None;
            let (mut cs, mut ds) = (Vec::new(), Vec::new());
            let mut compacted = SamplePspChunk::empty_compacted(chrom_id);
            for seg in segs {
                // Rows that are variable (∈ kept_all) and in this chunk's span.
                let positions = seg.positions();
                let mut keep = membership_mask(positions, &kept_all);
                for (k, &p) in keep.iter_mut().zip(positions.iter()) {
                    *k &= p >= chunk_start && p < cut;
                }
                if !keep.iter().any(|&b| b) {
                    continue;
                }
                // Decode-once (default): a segment that extends past this
                // chunk's cut is a straddler — the next chunk needs it too — so
                // inflate every heavy column once into the segment's cache and
                // slice both chunks from it (no repeat zstd decode). A
                // fully-consumed segment (used only here) takes the memory-frugal
                // one-column-at-a-time `set_variable_rows`. Under `--low-memory`
                // (`cache_straddlers == false`) every segment takes that frugal
                // path — straddlers re-decompress per chunk, no cache held.
                let chunk = if let Some(decoded) = seg.decoded_cache() {
                    seg.compact_rows(decoded, &keep)
                } else {
                    let dz = decompressor.get_or_insert_with(|| {
                        // PANIC-FREE: `zstd::bulk::Decompressor::new` is
                        // infallible on every supported platform.
                        new_column_decompressor().expect(
                            "zstd::bulk::Decompressor::new is infallible on supported platforms",
                        )
                    });
                    if cache_straddlers && seg.last_pos().is_some_and(|lp| lp >= cut) {
                        let decoded = seg.inflate_all(dz, &mut cs, &mut ds)?;
                        let chunk = seg.compact_rows(&decoded, &keep);
                        seg.cache_decoded(decoded);
                        chunk
                    } else {
                        seg.set_variable_rows(&keep, dz, &mut cs, &mut ds)?
                    }
                };
                compacted.append_range(&chunk, chunk_start, cut);
            }
            Ok(compacted)
        })
        .collect::<Result<_, _>>()?;
    Ok(RawCohortChunk {
        chunk_order,
        per_sample,
        ref_span,
    })
}

/// The producer's **read stage** (the second internal queue). Owns the N
/// readers, walks the precomputed `schedule` (same order the fold loop walks),
/// and decodes each sample's light-column segments **ahead** into per-sample
/// bounded channels so the fold stage never blocks on decode — overlapping the
/// light decode (the fold stage's dominant critical-path cost) with folding +
/// compaction.
///
/// Back-pressure without deadlock: each round decodes one segment for every
/// sample whose channel currently has room (`< depth`) in parallel on `pool`,
/// then sends (room was checked and only the fold drains, so the send never
/// blocks). When every non-exhausted channel is full it blocks on `wake` until
/// the fold drains one. A closed channel / `wake` (the fold stage stopped on an
/// error) ends the stage. Resets the readers at each interval boundary, in
/// lockstep with the fold loop's [`begin_interval`](CohortChunkIntegrator::begin_interval),
/// and emits one [`ReadMsg::IntervalEnd`] per sample per interval — so the fold
/// stage sees exactly the segment stream (and exhaustion points) the owned
/// decode would, byte-identically.
pub(crate) fn drive_read_stage<R: Read + Seek + Send>(
    mut readers: Vec<SamplePspReader<R>>,
    schedule: &[(u32, Range<u32>)],
    txs: &[crossbeam_channel::Sender<ReadMsg>],
    wake: &crossbeam_channel::Receiver<()>,
    pool: &rayon::ThreadPool,
    depth: usize,
) -> Result<(), ProducerError> {
    use rayon::prelude::*;
    let n = readers.len();
    let mut which = vec![false; n];
    for (chrom_id, interval) in schedule {
        let region_end = interval.end - 1; // reader region_end is inclusive
        for r in &mut readers {
            r.reset(*chrom_id, interval.start, region_end);
        }
        let mut exhausted = vec![false; n];
        loop {
            // Samples to refill this round: not exhausted, channel has room.
            let mut any = false;
            for s in 0..n {
                which[s] = !exhausted[s] && txs[s].len() < depth;
                any |= which[s];
            }
            if !any {
                if exhausted.iter().all(|&e| e) {
                    break; // interval drained
                }
                // Every non-exhausted channel is full — wait for the fold to
                // drain one (or stop if it's gone).
                if wake.recv().is_err() {
                    return Ok(());
                }
                continue;
            }
            // Decode one segment per selected sample, in parallel.
            let outcomes: Vec<(usize, Result<Option<TwoPhaseSegment>, PspReadError>)> = pool
                .install(|| {
                    readers
                        .par_iter_mut()
                        .enumerate()
                        .filter_map(|(s, r)| which[s].then(|| (s, r.next_two_phase())))
                        .collect()
                });
            for (s, res) in outcomes {
                let msg = match res? {
                    Some(seg) => ReadMsg::Segment(Box::new(seg)),
                    None => {
                        exhausted[s] = true;
                        ReadMsg::IntervalEnd
                    }
                };
                // Room was checked and only the fold drains, so this won't
                // block; an `Err` means the fold stage dropped the receiver.
                if txs[s].send(msg).is_err() {
                    return Ok(());
                }
            }
        }
    }
    Ok(())
}

/// Streaming cohort producer over **one covered interval** (appendix §B).
///
/// Owns the N [`SamplePspReader`]s. [`begin_interval`](Self::begin_interval)
/// points them at a `[start, end)` interval (+ its dust mask);
/// [`produce_chunk`](Self::produce_chunk) then yields one
/// [`RawCohortChunk`] (the chunk's variable rows in columnar per-sample form)
/// at a time, cut at safe gaps, until the interval is drained (`Ok(None)`).
/// `produce_chunk` takes the REF fetch as a closure so it stays decoupled and
/// unit-testable.
///
/// Per chunk: lockstep-read segments to the watermark → fold the **light**
/// columns → `derive_is_kept` → drop dust → cut at a safe gap → compact each
/// sample's variable rows → fetch the REF span → stamp `chunk_order`. Only the
/// variable positions are ever compacted (the memory invariant); buffered
/// columnar segments straddling a cut survive to the next chunk. The
/// columns→records rebuild + per-position merge run on the **caller**
/// ([`merge_compacted_samples`](crate::var_calling::variant_caller::merge_compacted_samples)),
/// not here.
/// One per-sample read result from the producer's **read stage** (the second
/// internal queue). `Segment` carries the next decoded (light columns eager,
/// heavy retained) segment for the sample; `IntervalEnd` marks that sample's
/// reader has no more blocks in the current covered interval. Boxed segment so
/// the enum is small (the channel buffers many).
pub(crate) enum ReadMsg {
    Segment(Box<TwoPhaseSegment>),
    IntervalEnd,
}

/// Where the integrator's [`read_samples`](CohortChunkIntegrator::read_samples)
/// gets each sample's next segment.
///
/// - `Owned`: the integrator owns the readers and decodes inline on the
///   producer pool. The byte-identity oracle path the tests + the test `run`
///   driver use.
/// - `Channel`: a separate **read stage** ([`drive_read_stage`]) owns the
///   readers and decodes ahead into per-sample bounded channels; the fold stage
///   just pulls. This decouples the light-column decode (the fold stage's
///   dominant critical-path cost) from folding so it overlaps fold + compact —
///   the production path.
enum ReadSource<R: Read + Seek> {
    /// Test/`run`-driver path only (the production pipeline always uses the
    /// `Channel` read stage) — kept out of `cfg(test)` so `R` stays used and
    /// the match arms compile unconditionally.
    #[allow(dead_code)]
    Owned(Vec<SamplePspReader<R>>),
    Channel {
        rx: Vec<crossbeam_channel::Receiver<ReadMsg>>,
        /// Pulsed after each receive so the read stage rechecks channel room
        /// (its back-pressure release). Bounded(1) — a pending pulse coalesces.
        wake: crossbeam_channel::Sender<()>,
    },
}

pub struct CohortChunkIntegrator<R: Read + Seek> {
    source: ReadSource<R>,
    n: usize,
    min_alt_obs: u32,
    target_variants: u32,

    // Current-interval state.
    chrom_id: u32,
    interval_end: u32, // exclusive
    next_chunk_start: u32,
    /// Per sample: buffered segments (carryover + read-ahead), in genomic
    /// order; dropped once fully below `next_chunk_start`. Each carries the
    /// light fold columns + the still-compressed heavy columns; held behind an
    /// `Arc` so a [`ChunkPlan`] can share a straddling segment with the live
    /// buffer (compaction inflates the heavy columns on its own stage).
    buffers: Vec<Vec<Arc<TwoPhaseSegment>>>,
    /// Per sample: reader exhausted within the current interval.
    exhausted: Vec<bool>,
    /// `false` until the first [`begin_interval`](Self::begin_interval). Gates
    /// the `Channel`-source leftover-drain: there is no prior interval to drain
    /// on the first call, and draining then would discard the first interval's
    /// real segments.
    started: bool,
    dust_mask: Vec<Range<u32>>,

    /// Monotonic genomic-order chunk counter (the writer's reorder ticket).
    /// **Not** reset per interval.
    chunk_order: u64,

    // Reusable scratch.
    fold: CohortSpanFold,
    is_kept: Vec<bool>,
}

// `R: Send` so the producer can decode the N per-sample readers in parallel
// (`read_samples`); the pipeline's `BufReader<File>` and the tests' `Cursor`
// are both `Send`.
impl<R: Read + Seek + Send> CohortChunkIntegrator<R> {
    /// Build a producer over the N per-sample readers. `min_alt_obs` is the
    /// cohort keep threshold; `target_variants` the soft per-chunk
    /// kept-position target (`0` ⇒ 1). (Sample names are no longer needed
    /// here — the per-position merge that used them now runs on the caller.)
    #[allow(dead_code)] // owned-reader path: tests + the `run` driver only
    pub fn new(readers: Vec<SamplePspReader<R>>, min_alt_obs: u32, target_variants: u32) -> Self {
        let n = readers.len();
        Self::with_source(ReadSource::Owned(readers), n, min_alt_obs, target_variants)
    }

    /// Build a producer that pulls each sample's segments from a per-sample
    /// channel (the production read-stage path; see [`ReadSource`]). `n` is the
    /// sample count (= `receivers.len()`). The caller spawns
    /// [`drive_read_stage`] to feed the channels and drives the fold loop over a
    /// precomputed interval schedule ([`cohort_intervals`]) — `covered_intervals`
    /// is unavailable here (the integrator owns no readers).
    pub(crate) fn new_streaming(
        receivers: Vec<crossbeam_channel::Receiver<ReadMsg>>,
        wake: crossbeam_channel::Sender<()>,
        min_alt_obs: u32,
        target_variants: u32,
    ) -> Self {
        let n = receivers.len();
        Self::with_source(
            ReadSource::Channel {
                rx: receivers,
                wake,
            },
            n,
            min_alt_obs,
            target_variants,
        )
    }

    fn with_source(
        source: ReadSource<R>,
        n: usize,
        min_alt_obs: u32,
        target_variants: u32,
    ) -> Self {
        Self {
            source,
            n,
            min_alt_obs,
            target_variants: target_variants.max(1),
            chrom_id: 0,
            interval_end: 0,
            next_chunk_start: 0,
            buffers: (0..n).map(|_| Vec::new()).collect(),
            exhausted: vec![false; n],
            started: false,
            dust_mask: Vec::new(),
            chunk_order: 0,
            fold: CohortSpanFold::new(),
            is_kept: Vec::new(),
        }
    }

    /// The chromosome's covered intervals — every sample's blocks for
    /// `chrom_id`, unioned and gap-merged at `max_group_span` (no file I/O,
    /// reads only the in-memory block indexes). Only available in `Owned` mode
    /// (the test/`run` path); the production read-stage path precomputes the
    /// whole-cohort schedule from [`covered_intervals_for`] before the readers
    /// move into the stage.
    #[cfg(test)]
    pub fn covered_intervals(&self, chrom_id: u32, max_group_span: u32) -> Vec<Range<u32>> {
        match &self.source {
            ReadSource::Owned(readers) => covered_intervals_for(readers, chrom_id, max_group_span),
            ReadSource::Channel { .. } => unreachable!("covered_intervals is owned-mode only"),
        }
    }

    /// Drive the whole cohort: walk `0..n_chromosomes`, and within each its
    /// covered intervals, emitting every [`RawCohortChunk`] (in genomic
    /// order, `chunk_order` monotonic across the whole run). The REF bytes and
    /// dust mask are **injected** — the chromosome→FASTA fetcher and the sdust
    /// computation are the caller's, keeping this loop pure orchestration.
    ///
    /// Test-only convenience driver: production
    /// ([`pipeline::run_var_calling`](crate::var_calling::pipeline::run_var_calling))
    /// runs the equivalent `covered_intervals` → `begin_interval` →
    /// `produce_chunk` loop itself, inside the producer rayon pool, so the
    /// per-chunk REF/dust fetchers can stay thread-local.
    ///
    /// - `ref_fetch(chrom_id, start_1based, len)` → the REF bases for a chunk's
    ///   span (monotonic-forward within a chromosome);
    /// - `dust_for(chrom_id, &interval)` → that interval's sorted, half-open
    ///   low-complexity mask (empty when complexity filtering is off);
    /// - `emit(chunk)` receives each produced chunk.
    #[cfg(test)]
    pub fn run<RefF, DustF, Emit>(
        &mut self,
        n_chromosomes: u32,
        max_group_span: u32,
        mut ref_fetch: RefF,
        mut dust_for: DustF,
        mut emit: Emit,
    ) -> Result<(), ProducerError>
    where
        RefF: FnMut(u32, u32, u32) -> Result<Vec<u8>, RefFetchError>,
        DustF: FnMut(u32, &Range<u32>) -> Vec<Range<u32>>,
        Emit: FnMut(RawCohortChunk),
    {
        for chrom_id in 0..n_chromosomes {
            let intervals = self.covered_intervals(chrom_id, max_group_span);
            for interval in intervals {
                let mask = dust_for(chrom_id, &interval);
                self.begin_interval(chrom_id, interval, mask);
                loop {
                    let mut rf = |s: u32, l: u32| ref_fetch(chrom_id, s, l);
                    match self.produce_chunk(&mut rf)? {
                        Some(chunk) => emit(chunk),
                        None => break,
                    }
                }
            }
        }
        Ok(())
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
        // `Owned`: reset the readers to the interval. `Channel`: the read stage
        // resets its own readers from the same schedule, in lockstep — the fold
        // stage only resets its per-interval buffer state below.
        match &mut self.source {
            ReadSource::Owned(readers) => {
                let region_end = interval.end - 1; // reader region_end is inclusive
                for r in readers {
                    r.reset(chrom_id, interval.start, region_end);
                }
            }
            // The read stage sends exactly one `IntervalEnd` per sample per
            // interval. The fold loop normally consumes it (a sample is read to
            // exhaustion), but `plan_chunk` can finish an interval *early* when a
            // cut reaches `interval_end` before the reader exhausts — leaving that
            // sample's `IntervalEnd` (and any trailing segments) queued. If they
            // are not drained, the *next* interval's first read pulls the stale
            // marker, marks the sample exhausted immediately, and the whole
            // schedule desyncs (every later interval then folds nothing — the
            // multi-interval staged-path bug). So before starting a new interval,
            // pull each not-yet-exhausted sample up to and including its pending
            // `IntervalEnd`, discarding the leftovers. Skipped on the first
            // interval (`started == false`): there is no prior interval, and the
            // queued messages are this interval's real data.
            ReadSource::Channel { rx, wake } if self.started => {
                for (s, rx_s) in rx.iter().enumerate() {
                    while !self.exhausted[s] {
                        match rx_s.recv() {
                            Ok(ReadMsg::Segment(_)) => {}
                            Ok(ReadMsg::IntervalEnd) | Err(_) => self.exhausted[s] = true,
                        }
                        // Free a slot so a read stage parked on a full channel
                        // can make progress (best-effort, coalesced).
                        let _ = wake.try_send(());
                    }
                }
            }
            ReadSource::Channel { .. } => {}
        }
        self.started = true;
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

    /// Read the next segment of each selected sample (`which[s] == true`) **in
    /// parallel**, pushing it to that sample's buffer (or marking the sample
    /// exhausted).
    ///
    /// The per-sample `.psp` decode (decompress every column of a segment) is
    /// the producer's heavy work; the readers are independent, so this is the
    /// §1-decode parallelism (`main` does the same with rayon-across-samples),
    /// kept inside the producer class. The cohort fold downstream is the serial
    /// barrier; reads are order-independent into its max-aggregation, so this
    /// does not affect byte-identity. (Whether the read-ahead made *progress*
    /// is `which.any()` at the call site — a selected sample always changes
    /// state, pushing a segment **or** exhausting, and an exhausting sample
    /// still raises the watermark, so progress must not hinge on a push.)
    fn read_samples(&mut self, which: &[bool]) -> Result<(), ProducerError> {
        debug_assert_eq!(which.len(), self.n);
        match &mut self.source {
            ReadSource::Owned(readers) => {
                use rayon::prelude::*;
                // Disjoint `&mut readers[s]` via `par_iter_mut().enumerate()`;
                // decode the selected readers concurrently, collect, apply.
                let outcomes: Vec<(usize, Result<Option<TwoPhaseSegment>, PspReadError>)> = readers
                    .par_iter_mut()
                    .enumerate()
                    .filter_map(|(s, reader)| which[s].then(|| (s, reader.next_two_phase())))
                    .collect();
                for (s, outcome) in outcomes {
                    match outcome? {
                        Some(seg) => self.buffers[s].push(Arc::new(seg)),
                        None => self.exhausted[s] = true,
                    }
                }
            }
            ReadSource::Channel { rx, wake } => {
                // Pull one already-decoded segment per selected sample (the
                // read stage decodes ahead concurrently). A closed channel or an
                // `IntervalEnd` marks the sample exhausted for this interval.
                // Each sample's channel is FIFO, so the segments arrive in
                // genomic order and the interval-boundary markers stay in step
                // with the fold loop's schedule walk — byte-identical to the
                // owned decode. Each receive frees a channel slot; pulse `wake`
                // so the read stage rechecks room (best-effort, coalesced).
                for s in 0..self.n {
                    if !which[s] {
                        continue;
                    }
                    match rx[s].recv() {
                        Ok(ReadMsg::Segment(seg)) => self.buffers[s].push(Arc::from(seg)),
                        Ok(ReadMsg::IntervalEnd) | Err(_) => self.exhausted[s] = true,
                    }
                    let _ = wake.try_send(());
                }
            }
        }
        Ok(())
    }

    /// Rebuild [`self.fold`] over the buffered light columns in
    /// `[next_chunk_start, w]` (all samples' data is complete there).
    ///
    /// L1 (perf review): this fold is the producer's serial floor — profiling
    /// pinned it as the dominant non-decode self-time, re-run every read-ahead
    /// iteration over the growing window. The aggregation is per-position
    /// **integer `max`** on ref-span / non-REF obs plus a position union — both
    /// associative and commutative — so each sample folds independently and the
    /// partials reduce via [`CohortSpanFold::merge`] in **any order** for a
    /// **byte-identical** result. Run it across the producer pool (this method
    /// executes inside `producer_pool.install(...)`).
    fn rebuild_fold(&mut self, w: u32) {
        use rayon::prelude::*;
        let lo_bound = self.next_chunk_start;
        let fold = self
            .buffers
            .par_iter()
            .map(|buf| {
                let mut f = CohortSpanFold::new();
                for chunk in buf {
                    let pos = chunk.positions();
                    let lo = pos.partition_point(|&p| p < lo_bound);
                    let hi = pos.partition_point(|&p| p <= w);
                    if lo < hi {
                        f.fold_sample_light(
                            &pos[lo..hi],
                            &chunk.ref_spans()[lo..hi],
                            &chunk.nonref_obs()[lo..hi],
                        );
                    }
                }
                f
            })
            .reduce(CohortSpanFold::new, |mut a, b| {
                a.merge(&b);
                a
            });
        self.fold = fold;
    }

    /// Read-ahead: pull segments (advancing the laggard samples) until the
    /// chunk has ~`target_variants` kept positions *and* a safe gap to cut at,
    /// or the interval is fully read. Leaves [`self.fold`] valid over
    /// `[next_chunk_start, watermark]`.
    fn fill_to_target(&mut self) -> Result<(), ProducerError> {
        loop {
            // Seed any non-exhausted sample that has no buffered segment (parallel).
            let seed: Vec<bool> = (0..self.n)
                .map(|s| !self.exhausted[s] && self.buffers[s].is_empty())
                .collect();
            if seed.iter().any(|&b| b) {
                self.read_samples(&seed)?;
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
            // Advance the laggards (coverage == W) to raise the watermark
            // (parallel). "Progress" is whether there *were* laggards to
            // process — each changes state (push or exhaust), and an exhausting
            // sample still raises W, so the loop must continue when any was
            // selected (matches the old serial `advanced` flag).
            let lag: Vec<bool> = (0..self.n)
                .map(|s| !self.exhausted[s] && self.coverage(s) == w)
                .collect();
            let any_laggard = lag.iter().any(|&b| b);
            self.read_samples(&lag)?;
            if !any_laggard {
                return Ok(());
            }
        }
    }

    /// Produce the next chunk of the current interval, or `Ok(None)` when the
    /// interval is drained — [`plan_chunk`](Self::plan_chunk) followed by the
    /// inline [`compact_plan`]. The production pipeline runs the two as
    /// **separate stages** (a bounded queue between them lets the next plan's
    /// fold overlap this plan's heavy decode); this convenience wrapper keeps
    /// the single-call shape the producer's tests + the test `run` driver use.
    #[cfg(test)]
    pub fn produce_chunk(
        &mut self,
        fetch_ref: &mut dyn FnMut(u32, u32) -> Result<Vec<u8>, RefFetchError>,
    ) -> Result<Option<RawCohortChunk>, ProducerError> {
        // `true`: exercise the straddler-cache (default) path in the
        // byte-identity oracle tests. The output is identical either way.
        match self.plan_chunk(fetch_ref)? {
            Some(plan) => Ok(Some(compact_plan(plan, true)?)),
            None => Ok(None),
        }
    }

    /// Fold/plan stage: advance the cohort fold to the next safe-gap cut, derive
    /// the chunk's variable positions, fetch its REF span, and capture the
    /// (still-compressed) per-sample segments overlapping `[chunk_start, cut)`
    /// as a [`ChunkPlan`]. `Ok(None)` once the interval is drained. The heavy
    /// decode runs later in [`compact_plan`]; this stage stays on the producer's
    /// critical path so REF fetch (monotonic-forward within a contig) and the
    /// serial fold bookkeeping keep their ordering.
    pub fn plan_chunk(
        &mut self,
        fetch_ref: &mut dyn FnMut(u32, u32) -> Result<Vec<u8>, RefFetchError>,
    ) -> Result<Option<ChunkPlan>, ProducerError> {
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
            // Forward progress: the loop sets `next_chunk_start = cut` on a
            // variant-free span and `continue`s, so a non-advancing cut would
            // spin forever. Guard it at release level (was a `debug_assert!`).
            if cut <= self.next_chunk_start {
                return Err(ProducerError::StalledCut {
                    next_chunk_start: self.next_chunk_start,
                    cut,
                });
            }

            // Keep decision over the (complete) fold, then drop dust.
            self.fold
                .derive_is_kept(self.min_alt_obs, &mut self.is_kept);
            drop_dust_masked(self.fold.positions(), &mut self.is_kept, &self.dust_mask);
            // `kept_all` = every variable (post-dust) position over the whole
            // fold — the mask the compact stage uses to slice each segment to
            // variable-only. `variable` is the prefix `< cut`: this chunk's
            // positions (for the REF span and the variant-free check). Both are
            // sorted ascending (the fold is).
            let kept_all: Vec<u32> = self
                .fold
                .positions()
                .iter()
                .zip(self.is_kept.iter())
                .filter(|&(_, &k)| k)
                .map(|(&p, _)| p)
                .collect();
            let variable: Vec<u32> = kept_all.iter().copied().take_while(|&p| p < cut).collect();

            let chunk_start = self.next_chunk_start;

            if variable.is_empty() {
                // Variant-free span — emit nothing, keep chunk_order gapless.
                self.next_chunk_start = cut;
                self.drop_buffers_below(cut);
                continue;
            }

            // Capture the segments overlapping [chunk_start, cut) (shared Arcs)
            // and fetch the REF span — both before dropping consumed buffers.
            let per_sample_segments = self.collect_overlapping(chunk_start, cut);
            let ref_span = self.fetch_ref_span(&variable, fetch_ref)?;
            self.next_chunk_start = cut;
            self.drop_buffers_below(cut);

            let chunk_order = self.chunk_order;
            self.chunk_order += 1;
            return Ok(Some(ChunkPlan {
                chunk_order,
                chrom_id: self.chrom_id,
                chunk_start,
                cut,
                kept_all,
                ref_span,
                per_sample_segments,
            }));
        }
    }

    /// Per sample: the buffered segments overlapping `[chunk_start, cut)`, as
    /// shared `Arc` clones (cheap refcount bump — a straddler is shared with the
    /// live buffer, never copied). The compact stage inflates their heavy
    /// columns off the fold stage's critical path.
    fn collect_overlapping(&self, chunk_start: u32, cut: u32) -> Vec<Vec<Arc<TwoPhaseSegment>>> {
        self.buffers
            .iter()
            .map(|buf| {
                buf.iter()
                    .filter(|seg| {
                        let first = seg.positions().first().copied().unwrap_or(0);
                        let last = seg.last_pos().unwrap_or(0);
                        first < cut && last >= chunk_start
                    })
                    .map(Arc::clone)
                    .collect()
            })
            .collect()
    }

    /// Drop buffered segments wholly below `cut` (fully consumed).
    fn drop_buffers_below(&mut self, cut: u32) {
        for buf in &mut self.buffers {
            buf.retain(|c| c.positions().last().is_some_and(|&p| p >= cut));
        }
    }

    /// Fetch the chunk's one contiguous REF span `[first_kept, max_reach]`
    /// (monotonic-forward across chunks). `max_reach` uses the cohort
    /// `max_ref_span` (the grouping reach), matching the old per-group fetch.
    fn fetch_ref_span(
        &self,
        variable: &[u32],
        fetch_ref: &mut dyn FnMut(u32, u32) -> Result<Vec<u8>, RefFetchError>,
    ) -> Result<RefSpan, ProducerError> {
        let first = variable[0];
        let positions = self.fold.positions();
        let spans = self.fold.max_ref_span();
        let mut max_reach = first;
        for &p in variable {
            // PANIC-FREE: every `variable` position is drawn from
            // `self.fold.positions()` (filtered by `is_kept`), so it is always
            // present in the fold.
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

    #[test]
    fn membership_mask_matches_binary_search() {
        // Byte-identity oracle: the merge-walk must equal the per-row
        // `haystack.binary_search(&p).is_ok()` it replaces (L7), across
        // disjoint, overlapping, prefix/suffix, empty, and singleton cases.
        let cases: &[(&[u32], &[u32])] = &[
            (&[], &[1, 2, 3]),
            (&[1, 2, 3], &[]),
            (&[1, 2, 3, 4, 5], &[2, 4]),
            (&[2, 4], &[1, 2, 3, 4, 5]),
            (&[1, 5, 9], &[1, 5, 9]),
            (&[1, 5, 9], &[2, 6, 10]),
            (&[10, 20, 30], &[5, 15, 30, 40]),
            (&[5], &[5]),
            (&[5], &[6]),
            (&[0, u32::MAX], &[0, 7, u32::MAX]),
        ];
        for (needles, haystack) in cases {
            let got = membership_mask(needles, haystack);
            let want: Vec<bool> = needles
                .iter()
                .map(|&p| haystack.binary_search(&p).is_ok())
                .collect();
            assert_eq!(got, want, "needles={needles:?} haystack={haystack:?}");
        }
    }

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
    // The columns→records rebuild helpers now live on the caller; the producer's
    // tests reuse them to reconstruct the cohort records a caller would build.
    use crate::var_calling::per_position_merger::PerPositionMerger;
    use crate::var_calling::sample_reader::SamplePspReader;
    use crate::var_calling::test_helpers::{allele, record};
    use crate::var_calling::types::CohortPileupRecord;
    use crate::var_calling::variant_caller::{KeptRecordIter, merge_compacted_samples, ok_record};
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
        let mut integ = CohortChunkIntegrator::new(readers, min_alt_obs, target_variants);
        integ.begin_interval(0, 1..INTERVAL_END, mask);

        let mut fetch = |_start: u32, len: u32| Ok::<_, RefFetchError>(vec![b'N'; len as usize]);
        let names = names(samples.len());
        let mut out: Vec<CohortPileupRecord> = Vec::new();
        let mut expected_order = 0u64;
        while let Some(chunk) = integ.produce_chunk(&mut fetch).unwrap() {
            assert_eq!(chunk.chunk_order, expected_order, "chunk_order gapless");
            expected_order += 1;
            // Reconstruct the cohort records the caller would build (the
            // columns→records + per-position merge moved off the producer).
            let records = merge_compacted_samples(&chunk.per_sample, &names).unwrap();
            assert!(!records.is_empty(), "empty chunks are never emitted");
            // RefSpan covers the chunk's first kept position.
            assert_eq!(chunk.ref_span.genomic_start, records[0].pos);
            out.extend(records);
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
    fn low_memory_cache_off_matches_reference() {
        // The `--low-memory` path (`compact_plan(.., cache_straddlers = false)`)
        // re-decompresses straddlers instead of caching them — it must still be
        // byte-identical to the reference (and hence to the default cache path).
        let (samples, bts) = cohort();
        let min_alt_obs = 1u32;
        let want = reference(&samples, min_alt_obs, &[]);
        let readers: Vec<_> = samples
            .iter()
            .zip(&bts)
            .map(|(recs, &bt)| {
                let reader = PspReader::new(Cursor::new(psp_bytes(recs, bt))).unwrap();
                SamplePspReader::new(reader, 0, 1, 1)
            })
            .collect();
        let mut integ = CohortChunkIntegrator::new(readers, min_alt_obs, 3);
        integ.begin_interval(0, 1..INTERVAL_END, Vec::new());
        let mut fetch = |_s: u32, len: u32| Ok::<_, RefFetchError>(vec![b'N'; len as usize]);
        let names = names(samples.len());
        let mut out: Vec<CohortPileupRecord> = Vec::new();
        while let Some(plan) = integ.plan_chunk(&mut fetch).unwrap() {
            let chunk = compact_plan(plan, false).unwrap();
            out.extend(merge_compacted_samples(&chunk.per_sample, &names).unwrap());
        }
        assert_eq!(out, want);
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

    #[test]
    fn merge_block_ranges_gap_merges_and_splits() {
        // Unsorted input; gap ≤ 10 merges, gap > 10 splits.
        let got = merge_block_ranges(vec![(200, 210), (10, 20), (25, 30)], 10);
        // 25 ≤ 20+10 → merge into [10,30]; 200 > 30+10 → new interval.
        assert_eq!(got, vec![10..31, 200..211]);
        assert!(merge_block_ranges(Vec::<(u32, u32)>::new(), 10).is_empty());
    }

    // --- 2b-wiring: the multi-chromosome / multi-interval run loop ---

    const MGS: u32 = 50; // max_group_span (≪ the inter-cluster gap, ≥ ref spans)

    /// Records for one `(chrom_id, sample)` at the given positions, same
    /// variant/coverage pattern as [`gen_sample`].
    fn gen_records(
        chrom_id: u32,
        s: u32,
        positions: impl Iterator<Item = u32>,
    ) -> Vec<PileupRecord> {
        positions
            .enumerate()
            .filter_map(|(k, pos)| {
                let i = k as u32;
                if (i + s).is_multiple_of(4) {
                    return None;
                }
                let ref_seq: &[u8] = if i.is_multiple_of(9) { b"AC" } else { b"A" };
                let mut alleles = vec![allele(ref_seq, 5, -1.0, &[])];
                if (i * 7 + s * 3).is_multiple_of(5) {
                    alleles.push(allele(b"T", 2, -1.5, &[u64::from(i) + 1]));
                }
                Some(PileupRecord::new(chrom_id, pos, alleles))
            })
            .collect()
    }

    /// One sample across 2 chromosomes; chrom 0 has two clusters separated by
    /// a > MGS gap (⇒ two covered intervals), chrom 1 has one.
    fn multichrom_sample(s: u32) -> Vec<PileupRecord> {
        let mut recs = Vec::new();
        recs.extend(gen_records(0, s, (0..80).map(|i| 10 + i * 5))); // chrom0 cluster1
        recs.extend(gen_records(0, s, (0..80).map(|i| 1000 + i * 5))); // chrom0 cluster2 (gap ~600)
        recs.extend(gen_records(1, s, (0..80).map(|i| 10 + i * 5))); // chrom1
        recs
    }

    fn psp_bytes_2chrom(records: &[PileupRecord], block_target: usize) -> Vec<u8> {
        // Small window grid so blocks cut on a fixed genomic grid (the
        // cross-sample segment alignment); the inter-cluster gap then lands on
        // a block boundary and splits the chromosome into two covered intervals.
        let mut w = PspWriter::new_with_block_layout(
            Cursor::new(Vec::new()),
            writer_header(2),
            block_target,
            128,
        )
        .expect("writer opens");
        for r in records {
            w.write_record(r).expect("write_record");
        }
        w.finish().expect("finish").into_inner()
    }

    /// One-shot reference over both chromosomes: fold + assemble per chrom,
    /// concatenated in chrom order (interval boundaries are safe gaps, so a
    /// whole-chrom fold equals the per-interval runs).
    fn reference_2chrom(
        samples: &[Vec<PileupRecord>],
        min_alt_obs: u32,
    ) -> Vec<CohortPileupRecord> {
        let mut out = Vec::new();
        for chrom_id in 0..2u32 {
            let per_chrom: Vec<Vec<PileupRecord>> = samples
                .iter()
                .map(|recs| {
                    recs.iter()
                        .filter(|r| r.chrom_id == chrom_id)
                        .cloned()
                        .collect()
                })
                .collect();
            out.extend(reference(&per_chrom, min_alt_obs, &[]));
        }
        out
    }

    fn run_cohort(
        samples: &[Vec<PileupRecord>],
        block_targets: &[usize],
        min_alt_obs: u32,
        target_variants: u32,
    ) -> Vec<CohortPileupRecord> {
        let readers: Vec<_> = samples
            .iter()
            .zip(block_targets)
            .map(|(recs, &bt)| {
                let bytes = psp_bytes_2chrom(recs, bt);
                SamplePspReader::new(PspReader::new(Cursor::new(bytes)).unwrap(), 0, 1, 1)
            })
            .collect();
        let mut integ = CohortChunkIntegrator::new(readers, min_alt_obs, target_variants);

        let names = names(samples.len());
        let mut out: Vec<CohortPileupRecord> = Vec::new();
        let mut expected_order = 0u64;
        integ
            .run(
                2,
                MGS,
                |_chrom, _start, len| Ok(vec![b'N'; len as usize]),
                |_chrom, _interval| Vec::new(),
                |chunk| {
                    assert_eq!(
                        chunk.chunk_order, expected_order,
                        "chunk_order gapless across run"
                    );
                    expected_order += 1;
                    out.extend(merge_compacted_samples(&chunk.per_sample, &names).unwrap());
                },
            )
            .unwrap();
        out
    }

    #[test]
    fn run_walks_chromosomes_and_intervals() {
        let samples = vec![
            multichrom_sample(0),
            multichrom_sample(1),
            multichrom_sample(2),
        ];
        let bts = vec![200usize, 350, 512];

        // Two covered intervals on chrom 0 (the gap splits it), one on chrom 1.
        {
            let readers: Vec<_> = samples
                .iter()
                .zip(&bts)
                .map(|(recs, &bt)| {
                    let bytes = psp_bytes_2chrom(recs, bt);
                    SamplePspReader::new(PspReader::new(Cursor::new(bytes)).unwrap(), 0, 1, 1)
                })
                .collect();
            let integ = CohortChunkIntegrator::new(readers, 1, 4);
            assert_eq!(
                integ.covered_intervals(0, MGS).len(),
                2,
                "chrom0 splits at the gap"
            );
            assert_eq!(integ.covered_intervals(1, MGS).len(), 1);
        }

        for &min_alt_obs in &[1u32, 2] {
            let want = reference_2chrom(&samples, min_alt_obs);
            assert!(!want.is_empty());
            for &target in &[1u32, 7, 100_000] {
                let got = run_cohort(&samples, &bts, min_alt_obs, target);
                assert_eq!(got, want, "min_alt_obs={min_alt_obs} target={target}");
            }
        }

        // Emit order: all chrom 0 before chrom 1, ascending within each chrom.
        let got = run_cohort(&samples, &bts, 1, 7);
        let boundary = got.iter().position(|r| r.chrom_id == 1).unwrap();
        assert!(got[..boundary].iter().all(|r| r.chrom_id == 0));
        assert!(got[boundary..].iter().all(|r| r.chrom_id == 1));
        assert!(got[..boundary].windows(2).all(|w| w[0].pos < w[1].pos));
        assert!(got[boundary..].windows(2).all(|w| w[0].pos < w[1].pos));
    }

    /// Drive the cohort through the **production staged path**: a
    /// [`drive_read_stage`] thread decodes ahead into per-sample channels and the
    /// fold loop pulls from them via a `Channel`-source integrator — mirroring
    /// the staged branch of `pipeline.rs` (read stage + bounded queues + the
    /// per-interval `begin_interval` / `plan_chunk` walk). Returns the assembled
    /// records in genomic order. `run_cohort` is the owned-reader analogue.
    fn run_cohort_staged(
        samples: &[Vec<PileupRecord>],
        block_targets: &[usize],
        min_alt_obs: u32,
        target_variants: u32,
    ) -> Vec<CohortPileupRecord> {
        const DEPTH: usize = 2; // mirrors pipeline.rs READ_QUEUE_DEPTH
        let readers: Vec<_> = samples
            .iter()
            .zip(block_targets)
            .map(|(recs, &bt)| {
                let bytes = psp_bytes_2chrom(recs, bt);
                SamplePspReader::new(PspReader::new(Cursor::new(bytes)).unwrap(), 0, 1, 1)
            })
            .collect();

        // Whole-cohort interval schedule (both chroms), shared by the read stage
        // and the fold loop — exactly as the pipeline builds it.
        let mut schedule: Vec<(u32, Range<u32>)> = Vec::new();
        for chrom_id in 0..2u32 {
            for iv in covered_intervals_for(&readers, chrom_id, MGS) {
                schedule.push((chrom_id, iv));
            }
        }

        let n = readers.len();
        let (seg_txs, seg_rxs): (
            Vec<crossbeam_channel::Sender<ReadMsg>>,
            Vec<crossbeam_channel::Receiver<ReadMsg>>,
        ) = (0..n)
            .map(|_| crossbeam_channel::bounded::<ReadMsg>(DEPTH))
            .unzip();
        let (wake_tx, wake_rx) = crossbeam_channel::bounded::<()>(1);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();

        let names = names(n);
        let mut out: Vec<CohortPileupRecord> = Vec::new();
        let schedule_ref = &schedule;
        let pool_ref = &pool;
        std::thread::scope(|scope| {
            let read_handle = scope.spawn(move || {
                drive_read_stage(readers, schedule_ref, &seg_txs, &wake_rx, pool_ref, DEPTH)
                    .unwrap();
            });

            let mut integ = CohortChunkIntegrator::<Cursor<Vec<u8>>>::new_streaming(
                seg_rxs,
                wake_tx,
                min_alt_obs,
                target_variants,
            );
            let mut fetch = |_c: u32, len: u32| Ok::<_, RefFetchError>(vec![b'N'; len as usize]);
            let mut expected_order = 0u64;
            for (chrom_id, iv) in schedule_ref {
                integ.begin_interval(*chrom_id, iv.clone(), Vec::new());
                while let Some(plan) = integ.plan_chunk(&mut fetch).unwrap() {
                    assert_eq!(plan.chunk_order, expected_order, "chunk_order gapless");
                    expected_order += 1;
                    let chunk = compact_plan(plan, true).unwrap();
                    out.extend(merge_compacted_samples(&chunk.per_sample, &names).unwrap());
                }
            }
            // Drop the integrator (releases the channel receivers + the wake
            // sender) so a read stage parked on a full queue or on `wake`
            // observes the disconnect and the scope can join.
            drop(integ);
            read_handle.join().unwrap();
        });
        out
    }

    #[test]
    fn staged_channel_path_matches_owned_across_intervals() {
        // Regression guard for the multi-interval staged-path desync bug. The
        // production path (read stage → per-sample channels → `Channel`-source
        // fold) sends one `IntervalEnd` marker per sample per interval. When an
        // interval finishes via a cut reaching `interval_end` *before* the reader
        // exhausts, that marker (and any trailing segments) stays queued; if
        // `begin_interval` doesn't drain it, every *later* interval pulls the
        // stale marker first, folds nothing, and the run emits almost nothing.
        //
        // Single-interval fixtures never trip it (the lone interval always drains
        // to exhaustion), which is why the tomato bench — one ~2 Mb region — and
        // every owned-path test missed it; the human bottle (≈1000 BED intervals)
        // collapsed to 3 records. So this drives the *channel* path across MANY
        // intervals (chrom0 ×2 + chrom1 ×1) and asserts byte-equality with the
        // owned reference. Small `target_variants` forces the early-cut leak;
        // large reads to exhaustion — both must match.
        // A SINGLE-sample cohort is the reliable trigger (and the real-world
        // shape that surfaced this — GIAB HG002): with one sample the cohort
        // watermark is that sample's own coverage, so an interval's last chunk is
        // cut at `interval_end` the moment the watermark reaches the last position
        // — no extra read, so the `IntervalEnd` is never consumed → leak. With
        // several unequal-length samples the `min` watermark forces a laggard read
        // to exhaustion, which consumes the marker and masks the bug; we test that
        // shape too (it must stay correct), but n=1 is what actually fails.
        let single = vec![multichrom_sample(0)];
        let multi = vec![
            multichrom_sample(0),
            multichrom_sample(1),
            multichrom_sample(2),
        ];
        let cases: &[(&[Vec<PileupRecord>], &[usize])] =
            &[(&single, &[200]), (&multi, &[200, 350, 512])];

        for (samples, bts) in cases {
            for &min_alt_obs in &[1u32, 2] {
                let want = reference_2chrom(samples, min_alt_obs);
                assert!(!want.is_empty(), "fixture should yield variants");
                for &target in &[1u32, 3, 7, 100_000] {
                    let got = run_cohort_staged(samples, bts, min_alt_obs, target);
                    // Cross-check the owned path too, so a future divergence in
                    // either direction is caught.
                    let owned = run_cohort(samples, bts, min_alt_obs, target);
                    assert_eq!(owned, want, "owned path sanity: target={target}");
                    assert_eq!(
                        got,
                        want,
                        "staged channel path diverged: n={} min_alt_obs={min_alt_obs} target={target}",
                        samples.len(),
                    );
                }
            }
        }
    }

    #[test]
    fn merge_reduce_tree_is_order_independent_with_ties() {
        // `rebuild_fold` reduces per-sample folds via `merge` over an arbitrary
        // rayon tree, so byte-identity needs `merge` associative + commutative.
        // `merge_matches_sequential_fold` only checks a two-way merge of
        // *distinct* values; this checks two different reduce-tree shapes over
        // folds that **tie** at a shared position (where a non-commutative
        // tie-break bug would hide).
        let a = (vec![10u32, 20], vec![1u32, 2], vec![0u32, 3]);
        let b = (vec![20u32, 30], vec![2u32, 1], vec![3u32, 1]);
        let c = (vec![20u32, 25], vec![2u32, 1], vec![3u32, 0]);
        let d = (vec![15u32, 20], vec![1u32, 2], vec![2u32, 3]);

        use std::slice::from_ref;
        // ((a ∘ b) ∘ c) ∘ d
        let mut left = fold(from_ref(&a));
        left.merge(&fold(from_ref(&b)));
        left.merge(&fold(from_ref(&c)));
        left.merge(&fold(from_ref(&d)));

        // (a ∘ b) ∘ (c ∘ d)
        let mut ab = fold(from_ref(&a));
        ab.merge(&fold(from_ref(&b)));
        let mut cd = fold(from_ref(&c));
        cd.merge(&fold(from_ref(&d)));
        ab.merge(&cd);

        assert_eq!(left.positions(), ab.positions());
        assert_eq!(left.max_ref_span, ab.max_ref_span);
        assert_eq!(left.max_nonref_obs, ab.max_nonref_obs);
        // And both equal the flat sequential fold of all four.
        let seq = fold(&[a, b, c, d]);
        assert_eq!(left.positions(), seq.positions());
        assert_eq!(left.max_ref_span, seq.max_ref_span);
        assert_eq!(left.max_nonref_obs, seq.max_nonref_obs);
    }

    #[test]
    fn produce_chunk_with_zero_samples_yields_none() {
        // Degenerate cohort (no readers): `watermark`/`all_exhausted` must drive
        // the interval straight to `Ok(None)` rather than spin or panic.
        let mut integ =
            CohortChunkIntegrator::new(Vec::<SamplePspReader<Cursor<Vec<u8>>>>::new(), 1, 4);
        integ.begin_interval(0, 1..100, Vec::new());
        let mut fetch = |_s: u32, l: u32| Ok::<_, RefFetchError>(vec![b'N'; l as usize]);
        assert!(integ.produce_chunk(&mut fetch).unwrap().is_none());
    }
}
