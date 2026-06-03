//! Chunk loader — read per-sample [`PileupRecord`] iterators into
//! columnar storage, apply the cohort-wide variant-position filter,
//! and compact the survivors into a [`MaterialisedChunk`] for the
//! pre-pass and worker.

use std::sync::atomic::{AtomicBool, Ordering};

use rayon::prelude::*;
use thiserror::Error;

use crate::pileup_record::PileupRecord;
use crate::psp::PspReadError;
use crate::var_calling::columns::{MaterialisedChunk, SampleColumns};

/// A per-sample source of pileup columns addressed by genomic span.
///
/// The streaming loader pulls one of these per sample. It hides the
/// `.psp` block layout: the loader peeks the next servable span end,
/// takes the cohort-wide `min` (the watermark every sample now covers),
/// and asks each source for the columns up to it. The production
/// implementor is
/// [`ColumnSpanReader`](crate::var_calling::column_span_reader::ColumnSpanReader);
/// tests use an in-memory `Vec<PileupRecord>` source.
pub trait SpanColumnSource {
    /// Inclusive end of the next span this source can serve without a
    /// new decode, or `None` once its region is exhausted. The loader
    /// takes `min` across sources to pick the watermark.
    fn peek_next_span(&self) -> Option<u32>;

    /// Append every not-yet-served record with position `< end` (in
    /// ascending order, within the source's region) to `out`.
    fn read_span(&mut self, end: u32, out: &mut SampleColumns) -> Result<(), PspReadError>;
}

/// Reusable scratch buffers for one iteration of the chunk loader.
///
/// The driver owns one of these and threads `&mut` to [`load_chunk_from_iters`]
/// per chunk; nothing inside is meant to survive across multiple
/// iterations except as memory the next iteration overwrites. All
/// inner buffers are `clear()`-ed on every call so the driver does
/// not have to.
///
/// **Fields.**
/// - [`Self::raw_per_sample`]: one [`SampleColumns`] per sample,
///   pre-sized at construction. Accumulates each sample's
///   pre-filter records (carryover prefix + freshly-loaded records)
///   before the variant filter compacts them into the output chunk.
/// - [`Self::position_union`]: sorted, dedup'd 1-based position
///   timeline across every sample in the current chunk.
/// - [`Self::has_variant_at`]: parallel to [`Self::position_union`];
///   each entry is `true` iff at least one sample at that position
///   carries a non-reference allele with `num_obs > 0`. Used to
///   propagate "variant" up to the group each position belongs to.
/// - [`Self::max_ref_span_at`]: parallel to [`Self::position_union`];
///   the max `ref_span` across samples that have a record at the
///   position. Drives the grouping simulation's reach calculation.
/// - [`Self::is_kept`]: parallel to [`Self::position_union`]; the
///   result of the grouping simulation — `true` iff the position
///   lands in a provisional group that contains at least one variant
///   position. These are the positions whose records survive into
///   the output chunk's `per_sample` columns.
#[derive(Debug)]
pub struct ChunkLoadScratch {
    raw_per_sample: Vec<SampleColumns>,
    position_union: Vec<u32>,
    has_variant_at: Vec<bool>,
    max_ref_span_at: Vec<u32>,
    is_kept: Vec<bool>,
    /// Variant-group index spans over [`Self::position_union`] —
    /// `(start_idx, end_idx)` half-open ranges built by the grouping
    /// simulation. Reused so the ALT-count filter can score whole
    /// groups without re-deriving boundaries.
    group_ranges: Vec<(usize, usize)>,
    /// Parallel to [`Self::position_union`]: the group index each
    /// position belongs to. Only populated when the ALT-count
    /// pushdown is active (`min_alt_obs > 1`); lets each sample's
    /// per-record scan attribute its observations to a group.
    group_of: Vec<u32>,
}

impl ChunkLoadScratch {
    /// Build scratch sized for a cohort of `n_samples` samples.
    pub fn with_n_samples(n_samples: usize) -> Self {
        Self {
            raw_per_sample: (0..n_samples).map(|_| SampleColumns::empty()).collect(),
            position_union: Vec::new(),
            has_variant_at: Vec::new(),
            max_ref_span_at: Vec::new(),
            is_kept: Vec::new(),
            group_ranges: Vec::new(),
            group_of: Vec::new(),
        }
    }

    /// Number of samples this scratch was sized for.
    pub fn n_samples(&self) -> usize {
        self.raw_per_sample.len()
    }

    /// Reset every internal buffer (raw per-sample columns, position
    /// union, variant decisions) while preserving allocated capacity.
    pub fn clear(&mut self) {
        for sample in &mut self.raw_per_sample {
            sample.clear();
        }
        self.position_union.clear();
        self.has_variant_at.clear();
        self.max_ref_span_at.clear();
        self.is_kept.clear();
        self.group_ranges.clear();
        self.group_of.clear();
    }
}

/// M13: per-call load-extent policy bundle for
/// [`load_chunk_from_iters`].
///
/// The four span/variant knobs always travel together —
/// `range_start`, `initial_span`, `target_variants`, `max_span` —
/// and a `u32`-only signature like
/// `load_chunk_from_iters(.., 10, 100, 0, 90, ..)` makes argument-
/// order mistakes silent. Bundling them into a named struct turns
/// the call site into a field-by-field literal and the helper into
/// a single-`extent: ChunkLoadExtent` parameter.
///
/// **Fields.**
/// - `chrom_id` — chromosome the chunk's records belong to.
/// - `range_start` — 1-based inclusive lower bound of the chunk's
///   load span.
/// - `initial_span` — BP span of the loader's first pull attempt.
/// - `target_variants` — soft lower bound on the post-filter variant
///   count (`0` disables the variant-bounded extension loop).
/// - `max_span` — hard cap on the chunk's BP span across every
///   extension iteration. `>= initial_span` is enforced; the loader
///   returns `ChunkLoadError::InvalidRange` otherwise.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChunkLoadExtent {
    pub chrom_id: u32,
    pub range_start: u32,
    pub initial_span: u32,
    pub target_variants: u32,
    pub max_span: u32,
}

impl ChunkLoadExtent {
    /// Build an extent from positional values. The struct-literal
    /// form is preferred at production call sites; this constructor
    /// exists for terse test fixtures.
    pub fn new(
        chrom_id: u32,
        range_start: u32,
        initial_span: u32,
        target_variants: u32,
        max_span: u32,
    ) -> Self {
        Self {
            chrom_id,
            range_start,
            initial_span,
            target_variants,
            max_span,
        }
    }
}

/// Per-chunk diagnostic stats returned by
/// [`load_chunk_from_iters`]. Drivers fold these into their own
/// cumulative counters; callers that don't care can drop the value.
// Mi1: `#[non_exhaustive]` — counter struct; future per-chunk stats
// can be added without breaking out-of-crate consumers.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ChunkLoadStats {
    /// Number of cohort-wide positions kept by the variant filter —
    /// the count of `true` entries in
    /// [`ChunkLoadScratch::is_kept`](ChunkLoadScratch). This is the
    /// per-chunk workload size used by the variant-count-bounded
    /// loading loop (Phase B prerequisite).
    pub variant_count: u32,
}

/// Errors surfaced by the chunk loader. Generic over `E`, the upstream
/// per-sample iterator's error type — typically `PspReadError` for the
/// production glue, or a unit-test error type for fixture-driven tests.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum ChunkLoadError<E> {
    /// A per-sample upstream iterator surfaced an error before the
    /// chunk's range was exhausted.
    #[error("failed to read pileup record for sample {sample_idx}")]
    UpstreamRead {
        sample_idx: usize,
        #[source]
        source: E,
    },

    /// `per_sample_iters.len()` did not match the cohort size the
    /// scratch was built for.
    #[error("per-sample iterator count {got} does not match scratch cohort size {expected}")]
    SampleCountMismatch { expected: usize, got: usize },

    /// `carryover.len()` did not match the cohort size the scratch
    /// was built for.
    #[error("carryover length {got} does not match scratch cohort size {expected}")]
    CarryoverLengthMismatch { expected: usize, got: usize },

    /// The chunk's requested range is empty or reversed
    /// (`start >= end`).
    #[error("chunk range {start}..{end} is empty or reversed")]
    InvalidRange { start: u32, end: u32 },

    /// An upstream iterator yielded a record whose `chrom_id` did
    /// not match the chunk's `chrom_id` — a PSP bug or a wiring
    /// mistake by the caller. The chunk loader does not cross
    /// chromosome boundaries.
    #[error(
        "sample {sample_idx} yielded record on chrom {got_chrom_id} \
         while chunk expects chrom {expected_chrom_id}"
    )]
    UnexpectedChromosome {
        sample_idx: usize,
        expected_chrom_id: u32,
        got_chrom_id: u32,
    },
}

/// Load one chunk from per-sample record iterators, apply the
/// cohort-wide variant-position filter, and compact the survivors
/// into `out`. The chunk's BP span grows adaptively to hit a target
/// post-filter variant count — the Phase B prerequisite that
/// decouples worker workload from PSP block size + variant density.
///
/// **Inputs.**
/// - `scratch`, `out`: caller-owned scratch and output chunk; both
///   are cleared at entry.
/// - `chrom_id`, `range_start`: the chunk's starting locus
///   (inclusive).
/// - `initial_span`: BP span of the loader's first pull attempt.
///   The first attempt pulls records with
///   `range_start <= pos < range_start + initial_span`.
/// - `target_variants`: the post-filter kept-position count the
///   loader tries to reach before stopping. The loader doubles its
///   attempt span until the count is met or `max_span` is reached.
///   `0` disables extension (one attempt only — the legacy
///   behaviour).
/// - `max_span`: hard cap on the chunk's BP span across all
///   extension iterations. The loader never pulls records past
///   `range_start + max_span`. Callers should clamp this to
///   chromosome length.
/// - `per_sample_iters`: one iterator per sample, in cohort order.
///   Iterators are wrapped in [`std::iter::Peekable`] internally
///   so the loader can pause at the current attempt's boundary
///   without consuming the next record, then resume on the next
///   extension iteration with a larger boundary.
/// - `carryover`: per-sample [`SampleColumns`] holding records from
///   the previous chunk's `>= safe_end` tail. These are *prepended*
///   to each sample's raw load. Cleared after consumption.
///
/// **Algorithm.**
/// 1. Drain `carryover[s]` into `scratch.raw_per_sample[s]`.
/// 2. Extension loop. With `attempt_end` starting at
///    `range_start + initial_span` and doubling each iteration up
///    to `max_span`: pull records with `pos < attempt_end` from
///    each iterator into `scratch.raw_per_sample[s]` (records at
///    or above `attempt_end` are left in the iterator for the
///    next pass), then rebuild the cohort-wide position union +
///    per-position `has_variant_at` + `max_ref_span_at`, re-run
///    the grouping simulation into `is_kept`, and count the
///    `true` entries. Break when the count reaches
///    `target_variants` or `attempt_end >= range_start + max_span`.
/// 3. Compact: for each sample, walk its raw rows and the
///    `position_union`/`is_kept` arrays in parallel; push every
///    row whose position is kept into `out.per_sample[s]`.
/// 4. Set `out.chrom_id = chrom_id`,
///    `out.range = range_start..final_attempt_end`,
///    `out.safe_end = final_attempt_end` (the pre-pass may revise
///    it down), leave `out.windows` empty for the pre-pass.
pub fn load_chunk_from_iters<I, E>(
    scratch: &mut ChunkLoadScratch,
    out: &mut MaterialisedChunk,
    extent: ChunkLoadExtent,
    per_sample_iters: Vec<I>,
    carryover: &mut [SampleColumns],
) -> Result<ChunkLoadStats, ChunkLoadError<E>>
where
    I: IntoIterator<Item = Result<PileupRecord, E>>,
    I::IntoIter: Send,
    E: Send,
{
    let ChunkLoadExtent {
        chrom_id,
        range_start,
        initial_span,
        target_variants,
        max_span,
    } = extent;
    let n_samples = scratch.n_samples();
    if per_sample_iters.len() != n_samples {
        return Err(ChunkLoadError::SampleCountMismatch {
            expected: n_samples,
            got: per_sample_iters.len(),
        });
    }
    if carryover.len() != n_samples {
        return Err(ChunkLoadError::CarryoverLengthMismatch {
            expected: n_samples,
            got: carryover.len(),
        });
    }
    if initial_span == 0 {
        return Err(ChunkLoadError::InvalidRange {
            start: range_start,
            end: range_start,
        });
    }
    // M8: validate the load-extent contract instead of silently
    // upgrading a too-small `max_span` to `initial_span` (the previous
    // `initial_span.min(max_span.max(initial_span))` was structurally
    // `initial_span`, masking a caller bug). Today's driver passes
    // `max_load_span >= initial_load_span` (B3), so this branch is
    // unreachable on the production path; the explicit error pins
    // the API contract for future callers.
    if max_span < initial_span {
        return Err(ChunkLoadError::InvalidRange {
            start: range_start,
            end: range_start.saturating_add(max_span),
        });
    }
    let max_attempt_end = range_start.saturating_add(max_span);

    scratch.clear();
    out.clear();
    out.chrom_id = chrom_id;
    // Make sure the output is sized for the cohort: a fresh
    // `MaterialisedChunk` may have come in with no samples slotted.
    if out.per_sample.len() != n_samples {
        out.per_sample.resize_with(n_samples, SampleColumns::empty);
    }

    // ── Step 1: drain carryover (records from the prior chunk that
    //            sat past its safe_end) into the per-sample raw
    //            buffers. ──
    for (sample_idx, carry) in carryover.iter_mut().enumerate() {
        let raw = &mut scratch.raw_per_sample[sample_idx];
        for row_idx in 0..carry.n_records() {
            raw.push_row_from(carry, row_idx);
        }
        carry.clear();
    }

    // ── Step 2: extension loop — pull + filter + variant-count
    //            check, doubling attempt_end until target_variants
    //            is hit or max_span is reached. ──
    let mut peekable_iters: Vec<std::iter::Peekable<I::IntoIter>> = per_sample_iters
        .into_iter()
        .map(|iter| iter.into_iter().peekable())
        .collect();
    let mut attempt_end = range_start.saturating_add(initial_span);
    let mut variant_count: u32;
    loop {
        let clamped_end = attempt_end.min(max_attempt_end);
        // Per-sample decode is independent — each sample has its own
        // iterator and raw buffer — so pull in parallel. At high N this
        // is the dominant cost (N× zstd decompress + record
        // materialise) and was previously serial. Determinism is
        // unaffected: each sample fills its own buffer.
        scratch
            .raw_per_sample
            .par_iter_mut()
            .zip(peekable_iters.par_iter_mut())
            .enumerate()
            .try_for_each(|(sample_idx, (raw, iter))| {
                pull_records_with_pos_under(
                    iter,
                    sample_idx,
                    chrom_id,
                    range_start,
                    clamped_end,
                    raw,
                )
            })?;
        attempt_end = clamped_end;

        // Legacy batch loader keeps the historical "any observed
        // non-REF" criterion (`min_alt_obs = 1`): it is the
        // streaming-vs-batch equivalence oracle and feeds the
        // contamination estimator, neither of which wants the ALT-count
        // pushdown. Only the streaming var-calling path opts in.
        rebuild_position_union_and_is_kept(scratch, 1);
        variant_count = scratch.is_kept.iter().filter(|kept| **kept).count() as u32;

        if variant_count >= target_variants {
            break;
        }
        if attempt_end >= max_attempt_end {
            break;
        }
        // Double the span up to max_attempt_end. Saturating arithmetic
        // keeps us safe at the u32 ceiling (a u32-overflowing span is
        // pathological — chromosomes top out near 2^28).
        let next_span = (attempt_end - range_start).saturating_mul(2);
        attempt_end = range_start.saturating_add(next_span).min(max_attempt_end);
    }

    out.range = range_start..attempt_end;
    out.safe_end = attempt_end;

    // ── Step 3: per-sample compact (run once on final accumulated
    //            data). ──
    // Per-sample compact is independent — each reads its own raw rows
    // plus the shared (immutable) union/is_kept and writes its own
    // column — so run in parallel. `out.per_sample` was cleared +
    // resized to `n_samples` above, so each slot is an empty column
    // ready to receive in row order (per-sample order is preserved, so
    // output is unchanged).
    let position_union = &scratch.position_union;
    let is_kept = &scratch.is_kept;
    out.per_sample
        .par_iter_mut()
        .zip(scratch.raw_per_sample.par_iter())
        .for_each(|(dst, raw)| {
            // Parallel walk over `raw.positions` and `position_union`:
            // both are sorted ascending, so we never need to backtrack.
            let mut union_idx = 0_usize;
            for row_idx in 0..raw.n_records() {
                let pos = raw.position_at(row_idx);
                while union_idx < position_union.len() && position_union[union_idx] < pos {
                    union_idx += 1;
                }
                debug_assert!(
                    union_idx < position_union.len() && position_union[union_idx] == pos,
                    "raw position {pos} missing from union — invariant broken",
                );
                if is_kept[union_idx] {
                    dst.push_row_from(raw, row_idx);
                }
            }
        });

    Ok(ChunkLoadStats { variant_count })
}

/// Rebuild `scratch.position_union`, `scratch.has_variant_at`,
/// `scratch.max_ref_span_at`, and `scratch.is_kept` from the
/// current `scratch.raw_per_sample` accumulator. Idempotent —
/// called once per extension iteration inside the loader's loop.
///
/// `min_alt_obs` is the downstream `min_alt_obs_per_sample` threshold
/// pushed up to load time. `<= 1` reproduces the historical criterion
/// exactly (keep a group iff any position carries any observed non-REF
/// allele — `1` is the `>= 1` special case). `> 1` activates the
/// conservative ALT-count pushdown: a group survives only if some
/// sample's non-REF observations, summed across the group's positions,
/// reach the threshold. This is a safe over-approximation of the
/// worker's per-(sample, unified-allele) filter (a unified allele's
/// obs is a subset-sum of the sample's group non-REF obs), so it never
/// drops a group the worker would have kept — the VCF is byte-identical
/// while the doomed-singleton groups never get materialised.
fn rebuild_position_union_and_is_kept(scratch: &mut ChunkLoadScratch, min_alt_obs: u32) {
    // Destructure for disjoint field borrows so the per-position scan
    // can write `has_variant_at` / `max_ref_span_at` in parallel while
    // reading `raw_per_sample` / `position_union`.
    let ChunkLoadScratch {
        raw_per_sample,
        position_union,
        has_variant_at,
        max_ref_span_at,
        is_kept,
        group_ranges,
        group_of,
    } = scratch;

    // ── Position union ──
    position_union.clear();
    for raw in raw_per_sample.iter() {
        position_union.extend_from_slice(&raw.positions);
    }
    // Parallel sort: at high N this union is N× a single sample's
    // positions. `par_sort_unstable` falls back to sequential for
    // small inputs and produces the same total order, so the dedup'd
    // union — and everything downstream — is unchanged.
    position_union.par_sort_unstable();
    position_union.dedup();

    let n = position_union.len();

    // ── Per-position scan: has_variant_at + max_ref_span_at ──
    // For each cohort position we walk the samples that have a record
    // there exactly once, computing both predicates in one pass:
    //   - `has_variant_at[i]` — at least one sample's record carries a
    //     non-REF allele with `num_obs > 0`. Used by Step 5 to decide
    //     whether a provisional group has any variant evidence.
    //   - `max_ref_span_at[i]` — max `ref_span` across samples with a
    //     record. Drives the grouping simulation's reach calculation.
    // Positions are independent, so the scan (O(positions × samples) —
    // the dominant fold cost at high N) runs in parallel. Every slot is
    // overwritten, so determinism/output is unchanged.
    has_variant_at.resize(n, false);
    max_ref_span_at.resize(n, 0);
    let raw_ref: &[SampleColumns] = raw_per_sample;
    has_variant_at
        .par_iter_mut()
        .zip(max_ref_span_at.par_iter_mut())
        .zip(position_union.par_iter())
        .for_each(|((hv, mrs), &pos)| {
            let mut has_variant = false;
            let mut max_ref_span: u32 = 0;
            for raw in raw_ref {
                if let Ok(row_idx) = raw.binary_search_position(pos) {
                    let ref_span = raw.ref_span_at(row_idx);
                    if ref_span > max_ref_span {
                        max_ref_span = ref_span;
                    }
                    if !has_variant && raw.has_observed_non_ref_allele_at(row_idx) {
                        has_variant = true;
                    }
                }
            }
            *hv = has_variant;
            *mrs = max_ref_span;
        });

    // ── Step 5: grouping simulation. ──
    // Walk the cohort-wide timeline in order, building provisional
    // groups under the same join rule the streaming grouper uses
    // (`pos <= group_end`, where `group_end` is the rolling max of
    // `pos + ref_span - 1` across the group). Mark every position in
    // a group that contains at least one variant position as kept;
    // positions in groups with no variants are dropped (matches the
    // streaming pipeline, where the per-group merger drops pure-REF
    // groups after merging).
    //
    // The point of this simulation — vs. the simpler per-position
    // "has any non-REF" check — is that the per-group merger gathers
    // per-(sample, allele) support from *every* position the group
    // span covers, including positions that are pure-REF in isolation
    // but inside the reach of another sample's MNP/DEL/INS. Dropping
    // those positions here would silently under-count REF evidence
    // for homref samples in multi-position groups.
    // Reset to all-false: the grouping sim only ever *sets* `true`, so
    // a bare `resize` would leak stale `true`s from a previous call when
    // the scratch is reused across rounds (the streaming loader). The
    // batch loader cleared the whole scratch per call and never hit this.
    is_kept.clear();
    is_kept.resize(n, false);

    // Build the group index spans (half-open `[start, end)` over
    // `position_union`) under the overlapping-reach join rule. The
    // group is the unit kept or dropped as a whole — splitting it would
    // under-count REF evidence for homref samples inside a multi-position
    // group (see above).
    group_ranges.clear();
    {
        let mut group_open = false;
        let mut group_start_idx: usize = 0;
        let mut group_end_pos: u32 = 0;
        for i in 0..n {
            let pos = position_union[i];
            let reach = pos.saturating_add(max_ref_span_at[i].max(1)) - 1;
            if group_open && pos <= group_end_pos {
                if reach > group_end_pos {
                    group_end_pos = reach;
                }
            } else {
                if group_open {
                    group_ranges.push((group_start_idx, i));
                }
                group_open = true;
                group_start_idx = i;
                group_end_pos = reach;
            }
        }
        if group_open {
            group_ranges.push((group_start_idx, n));
        }
    }

    let threshold = min_alt_obs.max(1);
    if threshold <= 1 {
        // Historical criterion: keep a group iff any position in it
        // carries an observed non-REF allele. `threshold == 1` is the
        // `group-sum >= 1` special case, so this is exact.
        for &(start, end) in group_ranges.iter() {
            if has_variant_at[start..end].iter().any(|&v| v) {
                is_kept[start..end].iter_mut().for_each(|slot| *slot = true);
            }
        }
    } else {
        // ALT-count pushdown. Score whole groups by each sample's
        // non-REF obs summed across the group; keep the group if any
        // sample reaches the threshold. Parallel over samples (the
        // architecture's existing axis); the per-group flag is an
        // order-independent OR, so the result is deterministic.
        group_of.clear();
        group_of.resize(n, u32::MAX);
        for (g, &(start, end)) in group_ranges.iter().enumerate() {
            group_of[start..end]
                .iter_mut()
                .for_each(|slot| *slot = g as u32);
        }
        let qualifies: Vec<AtomicBool> = (0..group_ranges.len())
            .map(|_| AtomicBool::new(false))
            .collect();
        let group_of_ref: &[u32] = group_of;
        let union_ref: &[u32] = position_union;
        raw_per_sample.par_iter().for_each(|raw| {
            let mut cur_group = u32::MAX;
            let mut sum: u32 = 0;
            for row_idx in 0..raw.n_records() {
                let pos = raw.position_at(row_idx);
                // Every record position contributed to the union, so the
                // search always hits; skip defensively if it ever doesn't.
                let Ok(uidx) = union_ref.binary_search(&pos) else {
                    continue;
                };
                let g = group_of_ref[uidx];
                if g != cur_group {
                    cur_group = g;
                    sum = 0;
                }
                if qualifies[g as usize].load(Ordering::Relaxed) {
                    continue;
                }
                sum = sum.saturating_add(raw.non_ref_obs_sum_at(row_idx));
                if sum >= threshold {
                    qualifies[g as usize].store(true, Ordering::Relaxed);
                }
            }
        });
        for (g, &(start, end)) in group_ranges.iter().enumerate() {
            if qualifies[g].load(Ordering::Relaxed) {
                is_kept[start..end].iter_mut().for_each(|slot| *slot = true);
            }
        }
    }
}

/// Streaming block loader — the memory-bounded replacement for
/// [`load_chunk_from_iters`].
///
/// It reads the cohort forward through [`SpanColumnSource`]s, folding
/// and **compacting incrementally** so that only the open variant group
/// at the working frontier is held in memory (the `pending` columns),
/// not the whole span needed to reach the variant target. Closed groups
/// are compacted into the output block and their raw records dropped on
/// every round — the fix for the "grown span × N samples" footprint of
/// the batch loader.
///
/// One block is produced per [`fill_block`](Self::fill_block) call;
/// `pending` carries the still-open group across calls within a covered
/// interval, and is reset at the interval boundary
/// ([`reset_interval`](Self::reset_interval)).
pub struct StreamingBlockLoader {
    /// `raw_per_sample` holds the *pending* (not-yet-finalised) records
    /// — the open group at the frontier plus the current round's read.
    /// The other fields are the per-round fold scratch.
    scratch: ChunkLoadScratch,
    /// Reused split buffer: receives the records that stay pending
    /// (position `>= cut`) while the closed prefix is compacted out,
    /// then swapped back into `scratch.raw_per_sample`.
    next_pending: Vec<SampleColumns>,
    /// Downstream `min_alt_obs_per_sample` threshold pushed up to the
    /// load-time variant-group filter. See
    /// [`rebuild_position_union_and_is_kept`].
    min_alt_obs: u32,
}

impl StreamingBlockLoader {
    pub fn with_n_samples(n_samples: usize, min_alt_obs: u32) -> Self {
        Self {
            scratch: ChunkLoadScratch::with_n_samples(n_samples),
            next_pending: (0..n_samples).map(|_| SampleColumns::empty()).collect(),
            min_alt_obs,
        }
    }

    /// Drop all pending records — called at a covered-interval /
    /// chromosome boundary, where the gap is wider than any group span
    /// so nothing carries over.
    pub fn reset_interval(&mut self) {
        for sample in &mut self.scratch.raw_per_sample {
            sample.clear();
        }
    }

    /// Produce one block into `out` from `sources`, covering
    /// `[range_start, …]` within `[…, interval_end]` (1-based
    /// inclusive). Accumulates closed variant groups until the kept
    /// (variable) position count reaches `target_variants`, then cuts at
    /// the start of the still-open group (a clean boundary by
    /// construction) and reserves that group as `pending` for the next
    /// call. At the interval's end every group is closed and flushed.
    ///
    /// Returns the number of distinct kept (variable) cohort positions in
    /// the produced block; `0` means the interval is exhausted with nothing
    /// left to emit (`out` is left empty).
    pub fn fill_block<S: SpanColumnSource + Send>(
        &mut self,
        sources: &mut [S],
        out: &mut MaterialisedChunk,
        chrom_id: u32,
        range_start: u32,
        interval_end: u32,
        target_variants: u32,
    ) -> Result<u32, PspReadError> {
        out.clear();
        out.chrom_id = chrom_id;
        let n_samples = sources.len();
        if out.per_sample.len() != n_samples {
            out.per_sample.resize_with(n_samples, SampleColumns::empty);
        }
        debug_assert_eq!(self.scratch.raw_per_sample.len(), n_samples);

        let mut kept_positions: u32 = 0;
        loop {
            let watermark = sources
                .iter()
                .filter_map(SpanColumnSource::peek_next_span)
                .min();
            let Some(watermark) = watermark else {
                // Interval exhausted: every remaining group is closed,
                // so flush all pending (cut past the interval end). A
                // zero return means nothing was left to emit.
                rebuild_position_union_and_is_kept(&mut self.scratch, self.min_alt_obs);
                let cut = interval_end.saturating_add(1);
                kept_positions += self.compact_closed_prefix(out, cut);
                out.range = range_start..cut;
                out.safe_end = cut;
                return Ok(kept_positions);
            };

            // Pull every sample up to the shared watermark (parallel —
            // each fills its own pending column).
            let read_end = watermark.saturating_add(1);
            self.scratch
                .raw_per_sample
                .par_iter_mut()
                .zip(sources.par_iter_mut())
                .try_for_each(|(pending, source)| source.read_span(read_end, pending))?;

            rebuild_position_union_and_is_kept(&mut self.scratch, self.min_alt_obs);
            let cut = find_block_cut(
                &self.scratch.position_union,
                &self.scratch.max_ref_span_at,
                watermark,
            );
            kept_positions += self.compact_closed_prefix(out, cut);

            if kept_positions >= target_variants {
                out.range = range_start..cut;
                out.safe_end = cut;
                return Ok(kept_positions);
            }
        }
    }

    /// Move records with position `< cut` out of `pending`: kept
    /// (variable / in-reach) ones append to `out`, the rest are dropped;
    /// records `>= cut` become the new `pending`. Returns the number of
    /// distinct kept cohort positions in `[…, cut)` (the block's variant
    /// tally). Reuses the fold's `position_union` / `is_kept`, valid for
    /// the current `pending`.
    fn compact_closed_prefix(&mut self, out: &mut MaterialisedChunk, cut: u32) -> u32 {
        let position_union = &self.scratch.position_union;
        let is_kept = &self.scratch.is_kept;

        out.per_sample
            .par_iter_mut()
            .zip(self.next_pending.par_iter_mut())
            .zip(self.scratch.raw_per_sample.par_iter())
            .for_each(|((dst, stays_pending), pending)| {
                stays_pending.clear();
                let mut union_idx = 0_usize;
                for row_idx in 0..pending.n_records() {
                    let pos = pending.position_at(row_idx);
                    if pos >= cut {
                        stays_pending.push_row_from(pending, row_idx);
                        continue;
                    }
                    while union_idx < position_union.len() && position_union[union_idx] < pos {
                        union_idx += 1;
                    }
                    debug_assert!(
                        union_idx < position_union.len() && position_union[union_idx] == pos,
                        "pending position {pos} missing from union — invariant broken",
                    );
                    if is_kept[union_idx] {
                        dst.push_row_from(pending, row_idx);
                    }
                }
            });

        // The split buffer is now the new pending; the old pending
        // columns are recycled as next round's split buffer.
        std::mem::swap(&mut self.scratch.raw_per_sample, &mut self.next_pending);

        position_union
            .iter()
            .zip(is_kept.iter())
            .take_while(|(pos, _)| **pos < cut)
            .filter(|(_, kept)| **kept)
            .count() as u32
    }
}

/// Find the block cut for the streaming loader: the start position of
/// the still-open variant group (the last group whose reach extends
/// past `watermark`), or `watermark + 1` when every group is closed.
///
/// Replays the same overlapping-reach grouping the fold uses
/// ([`rebuild_position_union_and_is_kept`]), so the cut always lands on
/// a clean group boundary — no group straddles it, which is what makes
/// the block split independent of where it falls.
fn find_block_cut(position_union: &[u32], max_ref_span_at: &[u32], watermark: u32) -> u32 {
    let reach = |pos: u32, span: u32| pos.saturating_add(span.max(1)) - 1;
    let mut sites = position_union.iter().zip(max_ref_span_at.iter());
    let Some((&first_pos, &first_span)) = sites.next() else {
        return watermark.saturating_add(1);
    };
    let mut group_start_pos = first_pos;
    let mut group_end_pos = reach(first_pos, first_span);
    for (&pos, &span) in sites {
        if pos <= group_end_pos {
            group_end_pos = group_end_pos.max(reach(pos, span));
        } else {
            group_start_pos = pos;
            group_end_pos = reach(pos, span);
        }
    }
    if group_end_pos > watermark {
        group_start_pos
    } else {
        watermark.saturating_add(1)
    }
}

/// Decision branches surfaced by the peek phase of
/// [`pull_records_with_pos_under`]. Split out so the immutable
/// peek-borrow of `iter` is released before the mutable
/// `iter.next()` consumes the same record. Mirroring the borrow
/// boundary in the type lets the loop body read straightforwardly.
enum PullDecision {
    /// No more records, or the next record is at `pos >=
    /// attempt_end`. Stop pulling without consuming.
    Stop,
    /// Next record is below `range_start`; consume + discard.
    SkipBelowRangeStart,
    /// Next record is in `[range_start, attempt_end)`; consume +
    /// push into `raw`.
    Consume,
    /// `iter.peek()` returned `Some(Err(_))`. Consume to extract
    /// the error and surface it.
    SurfaceUpstreamErr,
    /// Next record's `chrom_id` does not match the expected
    /// chromosome. Carries the observed value for the error.
    ReturnChromMismatch(u32),
}

/// Drain `iter` of every record with `range_start <= pos <
/// attempt_end`, pushing each into `raw`. Records before
/// `range_start` are consumed + skipped (the iterator's start may
/// be slightly before the chunk start because of PSP block
/// alignment). The first record at or above `attempt_end` is left
/// in the iterator — the loader's extension loop (step 3) resumes
/// the pull from there with a larger `attempt_end`.
fn pull_records_with_pos_under<I, E>(
    iter: &mut std::iter::Peekable<I>,
    sample_idx: usize,
    chrom_id: u32,
    range_start: u32,
    attempt_end: u32,
    raw: &mut SampleColumns,
) -> Result<(), ChunkLoadError<E>>
where
    I: Iterator<Item = Result<PileupRecord, E>>,
{
    loop {
        let decision = match iter.peek() {
            None => PullDecision::Stop,
            Some(Err(_)) => PullDecision::SurfaceUpstreamErr,
            Some(Ok(record)) => {
                if record.chrom_id != chrom_id {
                    PullDecision::ReturnChromMismatch(record.chrom_id)
                } else if record.pos >= attempt_end {
                    PullDecision::Stop
                } else if record.pos < range_start {
                    PullDecision::SkipBelowRangeStart
                } else {
                    PullDecision::Consume
                }
            }
        };
        match decision {
            PullDecision::Stop => return Ok(()),
            PullDecision::SkipBelowRangeStart => {
                iter.next();
            }
            PullDecision::Consume => {
                // M25 PANIC-FREE: `decision` was set from
                // `iter.peek()` above, which returned `Some(Ok(_))`
                // for this branch. `Peekable` guarantees the next
                // `next()` returns the same `Some(Ok(_))`; no other
                // mutation occurred (we hold `&mut iter` exclusively).
                let record = iter
                    .next()
                    .and_then(Result::ok)
                    .expect("peek returned Some(Ok(_))");
                raw.push_record(record);
            }
            PullDecision::SurfaceUpstreamErr => {
                // M25 PANIC-FREE: same invariant as Consume — `peek`
                // returned `Some(Err(_))` for this branch, so `next`
                // returns the same `Some(Err(_))`.
                let err = iter
                    .next()
                    .and_then(Result::err)
                    .expect("peek returned Some(Err(_))");
                return Err(ChunkLoadError::UpstreamRead {
                    sample_idx,
                    source: err,
                });
            }
            PullDecision::ReturnChromMismatch(got_chrom_id) => {
                return Err(ChunkLoadError::UnexpectedChromosome {
                    sample_idx,
                    expected_chrom_id: chrom_id,
                    got_chrom_id,
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::test_helpers::{
        allele, record, ref_obs, ref_plus_alt, ref_plus_unobserved_alt, run_loader,
    };
    use std::convert::Infallible;

    // ──────────────────────────────────────────────────────────────
    // Streaming loader equivalence — the streaming loader's kept
    // records (concatenated across all the blocks it emits) must equal
    // the batch loader's compacted output over the same region, for any
    // block granularity / variant target / cross-sample misalignment.
    // ──────────────────────────────────────────────────────────────

    /// In-memory [`SpanColumnSource`] over one sample's pre-sorted
    /// records, chopped into artificial "blocks" of `block_records` to
    /// exercise the watermark/round logic (and, with different sizes
    /// per sample, cross-sample misalignment).
    struct VecColumnSource {
        records: Vec<PileupRecord>,
        cursor: usize,
        region_end: u32,
        block_records: usize,
    }

    impl VecColumnSource {
        fn new(
            records: Vec<PileupRecord>,
            region_start: u32,
            region_end: u32,
            block_records: usize,
        ) -> Self {
            let cursor = records
                .iter()
                .position(|r| r.pos >= region_start)
                .unwrap_or(records.len());
            Self {
                records,
                cursor,
                region_end,
                block_records: block_records.max(1),
            }
        }
    }

    impl SpanColumnSource for VecColumnSource {
        fn peek_next_span(&self) -> Option<u32> {
            let head = self.records.get(self.cursor)?;
            if head.pos > self.region_end {
                return None;
            }
            let last = (self.cursor + self.block_records - 1).min(self.records.len() - 1);
            Some(self.records[last].pos.min(self.region_end))
        }

        fn read_span(&mut self, end: u32, out: &mut SampleColumns) -> Result<(), PspReadError> {
            let upper = end.min(self.region_end.saturating_add(1));
            while let Some(r) = self.records.get(self.cursor) {
                if r.pos >= upper {
                    break;
                }
                out.push_record(r.clone());
                self.cursor += 1;
            }
            Ok(())
        }
    }

    /// Batch reference: fold + compact the whole region in one shot.
    fn batch_kept(
        per_sample: &[Vec<PileupRecord>],
        start: u32,
        end_inclusive: u32,
    ) -> Vec<SampleColumns> {
        let n = per_sample.len();
        let span = end_inclusive - start + 1;
        let mut scratch = ChunkLoadScratch::with_n_samples(n);
        let mut chunk = MaterialisedChunk::with_n_samples(n);
        let iters: Vec<_> = per_sample
            .iter()
            .map(|rs| rs.clone().into_iter().map(Ok::<_, Infallible>))
            .collect();
        let mut carry: Vec<SampleColumns> = (0..n).map(|_| SampleColumns::empty()).collect();
        load_chunk_from_iters(
            &mut scratch,
            &mut chunk,
            ChunkLoadExtent::new(0, start, span, 0, span),
            iters,
            &mut carry,
        )
        .expect("batch load");
        chunk.per_sample
    }

    /// Streaming: emit blocks until the interval is exhausted and
    /// concatenate each sample's kept records across all blocks.
    fn streaming_kept(
        per_sample: &[Vec<PileupRecord>],
        start: u32,
        end_inclusive: u32,
        block_records_per_sample: &[usize],
        target_variants: u32,
    ) -> Vec<SampleColumns> {
        let n = per_sample.len();
        let mut sources: Vec<VecColumnSource> = per_sample
            .iter()
            .zip(block_records_per_sample)
            .map(|(rs, &bs)| VecColumnSource::new(rs.clone(), start, end_inclusive, bs))
            .collect();
        // Threshold 1 = historical criterion, matching `batch_kept`'s
        // oracle (which always passes 1). The ALT-count pushdown
        // (threshold > 1) is validated end-to-end against `main`'s VCF,
        // not by this streaming-vs-batch equivalence test.
        let mut loader = StreamingBlockLoader::with_n_samples(n, 1);
        let mut out = MaterialisedChunk::with_n_samples(n);
        let mut combined: Vec<SampleColumns> = (0..n).map(|_| SampleColumns::empty()).collect();
        let mut range_start = start;
        let mut blocks = 0;
        while loader
            .fill_block(
                &mut sources,
                &mut out,
                0,
                range_start,
                end_inclusive,
                target_variants,
            )
            .expect("fill_block")
            > 0
        {
            for (combined_s, per_sample_s) in combined.iter_mut().zip(&out.per_sample) {
                for row in 0..per_sample_s.n_records() {
                    combined_s.push_row_from(per_sample_s, row);
                }
            }
            range_start = out.safe_end;
            blocks += 1;
            assert!(blocks < 100_000, "fill_block did not terminate");
        }
        combined
    }

    /// Multi-sample fixture: SNPs, a homref sample at a variant
    /// position, a 3-bp MNP whose reach keeps in-span homref positions,
    /// and pure-REF positions that must drop — spread over many
    /// positions so block chopping splits groups across rounds.
    fn streaming_fixture() -> Vec<Vec<PileupRecord>> {
        let mut s0 = Vec::new();
        let mut s1 = Vec::new();
        for i in 0..60u32 {
            let pos = 100 + i * 13;
            match i % 4 {
                0 => {
                    // SNP in both.
                    s0.push(record(pos, ref_plus_alt(8, 4)));
                    s1.push(record(pos, ref_plus_alt(9, 3)));
                }
                1 => {
                    // Variant in s0 only; s1 homref at the same position.
                    s0.push(record(pos, ref_plus_alt(7, 5)));
                    s1.push(record(pos, vec![ref_obs(20)]));
                }
                2 => {
                    // 3-bp MNP in s0 reaching pos..pos+2; s1 plain REF at
                    // all three (in-reach homref must be kept).
                    s0.push(record(
                        pos,
                        vec![allele(b"AAA", 6, -1.0, &[]), allele(b"ACA", 5, -1.0, &[])],
                    ));
                    s1.push(record(pos, vec![ref_obs(7)]));
                    s1.push(record(pos + 1, vec![ref_obs(7)]));
                    s1.push(record(pos + 2, vec![ref_obs(7)]));
                }
                _ => {
                    // Pure REF in both (must drop) + an unobserved ALT.
                    s0.push(record(pos, ref_plus_unobserved_alt(15)));
                    s1.push(record(pos, vec![ref_obs(12)]));
                }
            }
        }
        vec![s0, s1]
    }

    #[test]
    fn streaming_matches_batch_aligned() {
        let f = streaming_fixture();
        let expected = batch_kept(&f, 1, 1_000_000);
        for &target in &[1u32, 4, 64, 100_000] {
            let got = streaming_kept(&f, 1, 1_000_000, &[5, 5], target);
            assert_eq!(got, expected, "target={target}");
        }
    }

    #[test]
    fn streaming_matches_batch_misaligned_blocks() {
        let f = streaming_fixture();
        let expected = batch_kept(&f, 1, 1_000_000);
        // Different artificial block sizes per sample → staggered
        // watermarks, the cross-sample case.
        for sizes in [[1usize, 7], [3, 2], [11, 4]] {
            let got = streaming_kept(&f, 1, 1_000_000, &sizes, 4);
            assert_eq!(got, expected, "sizes={sizes:?}");
        }
    }

    #[test]
    fn streaming_matches_batch_clamped_region() {
        let f = streaming_fixture();
        let (start, end) = (220u32, 600u32);
        let expected = batch_kept(&f, start, end);
        let got = streaming_kept(&f, start, end, &[3, 5], 4);
        assert_eq!(got, expected);
    }

    /// Helper: run the loader and return the post-load
    /// `ChunkLoadStats` along with the chunk. Mirrors `run_loader`
    /// but exposes the variant count and lets callers control the
    /// extension knobs.
    fn run_loader_with_stats(
        chrom_id: u32,
        range_start: u32,
        initial_span: u32,
        target_variants: u32,
        max_span: u32,
        per_sample_records: Vec<Vec<PileupRecord>>,
        mut carryover: Vec<SampleColumns>,
    ) -> (MaterialisedChunk, ChunkLoadStats) {
        let n_samples = per_sample_records.len();
        let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let mut out = MaterialisedChunk::with_n_samples(n_samples);
        let iters: Vec<_> = per_sample_records
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, Infallible>))
            .collect();
        let stats = load_chunk_from_iters(
            &mut scratch,
            &mut out,
            ChunkLoadExtent::new(
                chrom_id,
                range_start,
                initial_span,
                target_variants,
                max_span,
            ),
            iters,
            &mut carryover,
        )
        .expect("load_chunk_from_iters succeeded");
        (out, stats)
    }

    /// M29: `load_chunk_from_iters` returns `SampleCountMismatch`
    /// when the per-sample iter vec length differs from the scratch's
    /// `n_samples()`. Without this regression test a refactor that
    /// drops the check would surface as an `index out of range` panic
    /// inside `scratch.raw_per_sample[sample_idx]` instead of a clean
    /// typed error.
    #[test]
    fn loader_rejects_sample_count_mismatch() {
        let mut scratch = ChunkLoadScratch::with_n_samples(2);
        let mut out = MaterialisedChunk::with_n_samples(2);
        let mut carry = vec![SampleColumns::empty(); 2];
        let iters: Vec<_> = vec![Vec::<PileupRecord>::new()] // 1 iter, scratch is sized for 2
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, Infallible>))
            .collect();
        let err = load_chunk_from_iters(
            &mut scratch,
            &mut out,
            ChunkLoadExtent::new(0, 10, 90, 0, 90),
            iters,
            &mut carry,
        )
        .unwrap_err();
        assert!(matches!(
            err,
            ChunkLoadError::SampleCountMismatch {
                expected: 2,
                got: 1
            }
        ));
    }

    /// M29: same shape, on the carryover slice's length check.
    #[test]
    fn loader_rejects_carryover_length_mismatch() {
        let mut scratch = ChunkLoadScratch::with_n_samples(2);
        let mut out = MaterialisedChunk::with_n_samples(2);
        let mut carry = vec![SampleColumns::empty(); 1]; // mismatched length
        let iters: Vec<_> = (0..2)
            .map(|_| Vec::<PileupRecord>::new())
            .map(|rs| rs.into_iter().map(Ok::<_, Infallible>))
            .collect();
        let err = load_chunk_from_iters(
            &mut scratch,
            &mut out,
            ChunkLoadExtent::new(0, 10, 90, 0, 90),
            iters,
            &mut carry,
        )
        .unwrap_err();
        assert!(matches!(
            err,
            ChunkLoadError::CarryoverLengthMismatch {
                expected: 2,
                got: 1
            }
        ));
    }

    /// M8: `load_chunk_from_iters` rejects `max_span < initial_span`
    /// with `InvalidRange` rather than silently treating it as
    /// `max_span = initial_span` (the previous structurally-no-op
    /// `effective_initial_span = initial_span.min(max_span.max(initial_span))`).
    /// Today's driver passes `max_load_span >= initial_load_span`
    /// (B3), so this branch is unreachable on the production path;
    /// the test pins the API contract for future callers.
    #[test]
    fn loader_rejects_max_span_below_initial_span() {
        let mut scratch = ChunkLoadScratch::with_n_samples(1);
        let mut out = MaterialisedChunk::with_n_samples(1);
        let mut carry = vec![SampleColumns::empty(); 1];
        let iters: Vec<_> = vec![Vec::<PileupRecord>::new()]
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, Infallible>))
            .collect();
        // initial_span = 100, max_span = 50 — the loader cannot honour
        // a 100-BP load when capped at 50 BP.
        let err = load_chunk_from_iters(
            &mut scratch,
            &mut out,
            ChunkLoadExtent::new(0, 10, 100, 0, 50),
            iters,
            &mut carry,
        )
        .unwrap_err();
        assert!(matches!(err, ChunkLoadError::InvalidRange { .. }));
    }

    #[test]
    fn loader_variant_count_matches_kept_position_set() {
        // Sample 0 + 1: variants at pos 10 (sample 0 only) and pos 14
        // (sample 1 only). Pos 12 is REF-only in both samples (drops
        // out of the variant filter). Expected variant_count = 2.
        let s0 = vec![
            record(10, ref_plus_alt(3, 4)),
            record(12, vec![ref_obs(15)]),
        ];
        let s1 = vec![
            record(12, vec![ref_obs(15)]),
            record(14, ref_plus_alt(5, 6)),
        ];
        let (chunk, stats) = run_loader_with_stats(
            0,
            1,
            99,
            0,
            99,
            vec![s0, s1],
            vec![SampleColumns::empty(); 2],
        );
        assert_eq!(stats.variant_count, 2);
        // Sanity-check the chunk shape too.
        assert_eq!(chunk.per_sample[0].n_records(), 1);
        assert_eq!(chunk.per_sample[1].n_records(), 1);
    }

    #[test]
    fn loader_variant_count_is_zero_for_pure_ref_chunk() {
        // No samples carry any non-REF allele; the filter drops
        // every position and variant_count = 0.
        let s0 = vec![record(10, vec![ref_obs(20)])];
        let s1 = vec![record(20, vec![ref_obs(15)])];
        let (chunk, stats) = run_loader_with_stats(
            0,
            1,
            99,
            0,
            99,
            vec![s0, s1],
            vec![SampleColumns::empty(); 2],
        );
        assert_eq!(stats.variant_count, 0);
        assert_eq!(chunk.per_sample[0].n_records(), 0);
        assert_eq!(chunk.per_sample[1].n_records(), 0);
    }

    #[test]
    fn loader_variant_count_includes_homref_positions_inside_mnp_reach() {
        // A 3-bp MNP at position 10 makes positions 10/11/12 all
        // "kept" by the grouping-simulation filter (matching the
        // streaming pipeline's behaviour — see
        // `loader_keeps_homref_position_inside_mnp_reach`). The
        // variant_count therefore reflects the kept set, not the
        // "has any non-REF" set.
        let s0 = vec![record(
            10,
            vec![allele(b"AAA", 5, -1.0, &[]), allele(b"ACA", 4, -1.0, &[])],
        )];
        let s1 = vec![
            record(10, vec![ref_obs(7)]),
            record(11, vec![ref_obs(7)]),
            record(12, vec![ref_obs(7)]),
        ];
        let (_, stats) = run_loader_with_stats(
            0,
            1,
            99,
            0,
            99,
            vec![s0, s1],
            vec![SampleColumns::empty(); 2],
        );
        assert_eq!(stats.variant_count, 3);
    }

    // ──────────────────────────────────────────────────────────────
    // Phase B step 3 — extension loop tests.
    // ──────────────────────────────────────────────────────────────

    #[test]
    fn loader_extension_doubles_span_until_target_variants_reached() {
        // Variants at positions 10, 110, 210, 310, 410, 510. Each
        // variant is one position (ref_span=1). With initial_span=50
        // (covers only the first variant) and target_variants=4,
        // the loader must extend: 50 → 100 (still 1) → 200 (2) →
        // 400 (4 = target). Final attempt_end = range_start + 400 =
        // 410. Final chunk has 4 kept positions (10, 110, 210, 310).
        let positions = [10, 110, 210, 310, 410, 510];
        let s0: Vec<_> = positions
            .iter()
            .map(|&p| record(p, ref_plus_alt(3, 4)))
            .collect();
        let (chunk, stats) = run_loader_with_stats(
            0,
            1,
            50,   // initial_span
            4,    // target_variants
            1000, // max_span (plenty of headroom)
            vec![s0],
            vec![SampleColumns::empty()],
        );
        assert!(
            stats.variant_count >= 4,
            "expected at least 4 variants after extension, got {}",
            stats.variant_count
        );
        assert!(
            chunk.range.end >= 311,
            "expected range.end past variant at 310, got {}",
            chunk.range.end
        );
    }

    #[test]
    fn loader_extension_stops_at_max_span_under_target() {
        // Sparse fixture: only one variant in the whole chromosome.
        // With target=10 the loader cannot satisfy the target; it
        // must stop at max_span and return whatever was loaded
        // (graceful degradation).
        let s0 = vec![record(50, ref_plus_alt(3, 4))];
        let (chunk, stats) = run_loader_with_stats(
            0,
            1,
            20,  // initial_span
            10,  // target_variants (impossible — only 1 variant exists)
            100, // max_span (cap)
            vec![s0],
            vec![SampleColumns::empty()],
        );
        // Variant count is 1 (whatever was loaded), and chunk.range.end
        // == range_start + max_span.
        assert_eq!(stats.variant_count, 1);
        assert_eq!(chunk.range.end, 101);
    }

    #[test]
    fn loader_no_extension_when_target_is_zero() {
        // target_variants = 0 disables extension: the loader runs
        // a single pull attempt covering initial_span and returns
        // whatever was found. Same shape as the legacy behaviour
        // every other test uses.
        let s0 = vec![
            record(10, ref_plus_alt(3, 4)),
            record(60, ref_plus_alt(3, 4)),
            record(110, ref_plus_alt(3, 4)),
        ];
        let (chunk, stats) = run_loader_with_stats(
            0,
            1,
            50,  // initial_span — only first variant lands inside
            0,   // target_variants = 0 disables extension
            500, // max_span (would allow extension if target > 0)
            vec![s0],
            vec![SampleColumns::empty()],
        );
        assert_eq!(stats.variant_count, 1);
        assert_eq!(chunk.range.end, 51);
        assert_eq!(chunk.per_sample[0].n_records(), 1);
    }

    #[test]
    fn loader_no_extension_when_target_already_met_by_initial_span() {
        // Dense initial span already exceeds target — zero
        // extensions needed.
        let s0 = vec![
            record(10, ref_plus_alt(3, 4)),
            record(20, ref_plus_alt(3, 4)),
            record(30, ref_plus_alt(3, 4)),
            record(40, ref_plus_alt(3, 4)),
            record(50, ref_plus_alt(3, 4)),
        ];
        let (chunk, stats) = run_loader_with_stats(
            0,
            1,
            100,  // initial_span covers all 5 variants
            3,    // target_variants — already met after first pass
            1000, // max_span (wouldn't trigger anyway)
            vec![s0],
            vec![SampleColumns::empty()],
        );
        assert_eq!(stats.variant_count, 5);
        // No extension: chunk.range.end == range_start + initial_span.
        assert_eq!(chunk.range.end, 101);
    }

    #[test]
    fn loader_drops_positions_with_no_observed_non_ref_allele() {
        // Sample 0 has REF-only at pos 10; both samples carry REF +
        // unobserved ALT at pos 12; sample 1 has REF + observed ALT
        // at pos 14. Expect pos 14 in the output; pos 10 and pos 12
        // dropped.
        let s0 = vec![
            record(10, vec![ref_obs(20)]),
            record(12, ref_plus_unobserved_alt(15)),
        ];
        let s1 = vec![
            record(12, ref_plus_unobserved_alt(15)),
            record(14, ref_plus_alt(5, 6)),
        ];
        let chunk = run_loader(0, 1..100, vec![s0, s1], vec![SampleColumns::empty(); 2]);
        assert_eq!(chunk.range, 1..100);
        assert_eq!(chunk.safe_end, 100);
        assert_eq!(chunk.per_sample[0].n_records(), 0);
        assert_eq!(chunk.per_sample[1].n_records(), 1);
        assert_eq!(chunk.per_sample[1].position_at(0), 14);
    }

    #[test]
    fn loader_keeps_homref_position_inside_mnp_reach() {
        // Sample 0 has a 3-bp MNP at position 10 (REF=AAA, ALT=ACA).
        // Sample 1 has plain SNP-shape records at 10/11/12 — pure REF.
        // Cohort-wide: position 10 is variant (sample 0 has ALT). But
        // positions 11 and 12 are pure-REF in isolation (sample 0's
        // 3-bp record covers them; sample 1 has REF-only records).
        // The grouping simulation extends the kept set to 11 and 12
        // because they fall inside sample 0's record's reach (10..12).
        let s0 = vec![record(
            10,
            vec![allele(b"AAA", 5, -1.0, &[]), allele(b"ACA", 4, -1.0, &[])],
        )];
        let s1 = vec![
            record(10, vec![ref_obs(7)]),
            record(11, vec![ref_obs(7)]),
            record(12, vec![ref_obs(7)]),
        ];
        let chunk = run_loader(0, 1..100, vec![s0, s1], vec![SampleColumns::empty(); 2]);

        // Sample 0 only has a record at 10. Kept.
        assert_eq!(chunk.per_sample[0].n_records(), 1);
        assert_eq!(chunk.per_sample[0].position_at(0), 10);
        // Sample 1 has records at 10, 11, 12 — all three kept because
        // they're inside sample 0's variant record's reach.
        assert_eq!(chunk.per_sample[1].n_records(), 3);
        assert_eq!(chunk.per_sample[1].position_at(0), 10);
        assert_eq!(chunk.per_sample[1].position_at(1), 11);
        assert_eq!(chunk.per_sample[1].position_at(2), 12);
    }

    #[test]
    fn loader_keeps_homref_sample_at_variant_position() {
        // Position 50: sample 0 homref (REF only), sample 1 has the
        // ALT. The merger needs sample 0's record at 50 for joint
        // likelihoods, so it must survive the filter.
        let s0 = vec![record(50, vec![ref_obs(20)])];
        let s1 = vec![record(50, ref_plus_alt(2, 6))];
        let chunk = run_loader(0, 1..100, vec![s0, s1], vec![SampleColumns::empty(); 2]);
        assert_eq!(chunk.per_sample[0].n_records(), 1);
        assert_eq!(chunk.per_sample[1].n_records(), 1);
        assert_eq!(chunk.per_sample[0].position_at(0), 50);
        assert_eq!(chunk.per_sample[1].position_at(0), 50);
    }

    #[test]
    fn loader_prepends_carryover_in_position_order() {
        // Carryover from the previous chunk: sample 0 has a variant
        // at pos 5; sample 1 has nothing. New chunk range starts at
        // pos 10. Expect output to carry pos 5 + pos 12 in order
        // for sample 0.
        let mut carry_s0 = SampleColumns::empty();
        carry_s0.push_record(record(5, ref_plus_alt(3, 4)));
        let carry = vec![carry_s0, SampleColumns::empty()];
        let s0 = vec![record(12, ref_plus_alt(2, 3))];
        let s1 = vec![record(12, vec![ref_obs(10)])];
        let chunk = run_loader(0, 10..100, vec![s0, s1], carry);

        assert_eq!(chunk.per_sample[0].n_records(), 2);
        assert_eq!(chunk.per_sample[0].position_at(0), 5);
        assert_eq!(chunk.per_sample[0].position_at(1), 12);
        assert_eq!(chunk.per_sample[1].n_records(), 1);
        assert_eq!(chunk.per_sample[1].position_at(0), 12);
    }

    #[test]
    fn loader_clears_carryover_after_draining_it() {
        let mut carry_s0 = SampleColumns::empty();
        carry_s0.push_record(record(5, ref_plus_alt(3, 4)));
        let mut carry = vec![carry_s0, SampleColumns::empty()];
        let s0 = vec![record(12, ref_plus_alt(2, 3))];
        let s1 = vec![record(12, vec![ref_obs(10)])];

        let n_samples = 2;
        let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let mut out = MaterialisedChunk::with_n_samples(n_samples);
        let iters: Vec<_> = vec![s0, s1]
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>))
            .collect();
        load_chunk_from_iters(
            &mut scratch,
            &mut out,
            ChunkLoadExtent::new(0, 10, 90, 0, 90),
            iters,
            &mut carry,
        )
        .unwrap();

        assert_eq!(carry[0].n_records(), 0);
        assert_eq!(carry[1].n_records(), 0);
    }

    #[test]
    fn loader_stops_pulling_at_range_end() {
        let s0 = vec![
            record(15, ref_plus_alt(2, 3)),
            record(60, ref_plus_alt(2, 3)),
            record(150, ref_plus_alt(2, 3)),
        ];
        let chunk = run_loader(0, 10..100, vec![s0], vec![SampleColumns::empty()]);
        assert_eq!(chunk.per_sample[0].n_records(), 2);
        assert_eq!(chunk.per_sample[0].position_at(0), 15);
        assert_eq!(chunk.per_sample[0].position_at(1), 60);
    }

    #[test]
    fn loader_rejects_inverted_range() {
        let n_samples = 1;
        let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let mut out = MaterialisedChunk::with_n_samples(n_samples);
        let mut carry = vec![SampleColumns::empty()];
        let iters: Vec<_> = vec![Vec::<PileupRecord>::new()]
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>))
            .collect();
        let err = load_chunk_from_iters(
            &mut scratch,
            &mut out,
            ChunkLoadExtent::new(0, 100, 0, 0, 0),
            iters,
            &mut carry,
        )
        .expect_err("zero initial_span rejected");
        assert!(matches!(err, ChunkLoadError::InvalidRange { .. }));
    }

    #[test]
    fn loader_rejects_record_from_unexpected_chromosome() {
        let n_samples = 1;
        let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let mut out = MaterialisedChunk::with_n_samples(n_samples);
        let mut carry = vec![SampleColumns::empty()];
        let bogus = PileupRecord::new(99, 50, ref_plus_alt(2, 3));
        let iters: Vec<_> = vec![vec![bogus]]
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>))
            .collect();
        let err = load_chunk_from_iters(
            &mut scratch,
            &mut out,
            ChunkLoadExtent::new(0, 10, 90, 0, 90),
            iters,
            &mut carry,
        )
        .expect_err("wrong-chromosome record rejected");
        assert!(matches!(
            err,
            ChunkLoadError::UnexpectedChromosome {
                expected_chrom_id: 0,
                got_chrom_id: 99,
                ..
            }
        ));
    }

    #[test]
    fn loader_surfaces_upstream_errors() {
        let n_samples = 1;
        let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let mut out = MaterialisedChunk::with_n_samples(n_samples);
        let mut carry = vec![SampleColumns::empty()];
        let iters: Vec<Box<dyn Iterator<Item = Result<PileupRecord, &'static str>> + Send>> =
            vec![Box::new(std::iter::once(Err("psp blew up")))];
        let err = load_chunk_from_iters(
            &mut scratch,
            &mut out,
            ChunkLoadExtent::new(0, 10, 90, 0, 90),
            iters,
            &mut carry,
        )
        .expect_err("upstream error surfaced");
        assert!(matches!(
            err,
            ChunkLoadError::UpstreamRead {
                sample_idx: 0,
                source: "psp blew up"
            }
        ));
    }

    #[test]
    fn loader_clears_output_chunk_between_iterations() {
        // First chunk: one variant at pos 14. Second chunk
        // overwrites the same `out` buffer: range moves to a fresh
        // window with a different chrom_id and a different variant.
        let n_samples = 1;
        let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let mut out = MaterialisedChunk::with_n_samples(n_samples);
        let mut carry = vec![SampleColumns::empty()];

        let iters_a: Vec<_> = vec![vec![record(14, ref_plus_alt(2, 3))]]
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>))
            .collect();
        load_chunk_from_iters(
            &mut scratch,
            &mut out,
            ChunkLoadExtent::new(0, 10, 40, 0, 40),
            iters_a,
            &mut carry,
        )
        .unwrap();
        assert_eq!(out.chrom_id, 0);
        assert_eq!(out.range, 10..50);
        assert_eq!(out.per_sample[0].n_records(), 1);

        let iters_b: Vec<_> = vec![vec![PileupRecord::new(7, 220, ref_plus_alt(3, 7))]]
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>))
            .collect();
        load_chunk_from_iters(
            &mut scratch,
            &mut out,
            ChunkLoadExtent::new(7, 200, 100, 0, 100),
            iters_b,
            &mut carry,
        )
        .unwrap();
        assert_eq!(out.chrom_id, 7);
        assert_eq!(out.range, 200..300);
        assert_eq!(out.per_sample[0].n_records(), 1);
        assert_eq!(out.per_sample[0].position_at(0), 220);
    }
}
