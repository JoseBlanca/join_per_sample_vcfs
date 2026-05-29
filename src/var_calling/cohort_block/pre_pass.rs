//! Pre-pass — pick the chunk's `safe_end`, split records past it
//! into the carryover for the next chunk, and partition
//! `[range.start, safe_end)` into the worker windows the chunk loop
//! will dispatch.

use thiserror::Error;

use crate::var_calling::cohort_block::columns::{MaterialisedChunk, SampleColumns};

/// Reusable scratch for the [`fix_boundaries`] pre-pass — the
/// cohort-wide `(position, max_reach)` timeline plus its prefix-max
/// helper column. The driver owns one of these next to the chunk
/// loader's scratch.
#[derive(Debug, Default)]
pub struct FixBoundariesScratch {
    /// `(position, max_reach)` — one entry per **unique** position in
    /// the post-filter chunk, sorted ascending by position. `max_reach`
    /// is the maximum value of `position + ref_span - 1` across the
    /// samples that have a record at this position.
    timeline: Vec<(u32, u32)>,
    /// Running maximum of `timeline[..=i].1` — `prefix_max_reach[i]`
    /// is the reach of the longest allele at any position ≤
    /// `timeline[i].0`. Lets the safe-gap scanner check the
    /// "no allele below crosses the candidate boundary" rule in O(1)
    /// per gap.
    prefix_max_reach: Vec<u32>,
}

impl FixBoundariesScratch {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn clear(&mut self) {
        self.timeline.clear();
        self.prefix_max_reach.clear();
    }
}

/// Errors surfaced by the [`fix_boundaries`] pre-pass.
#[non_exhaustive]
#[derive(Error, Debug, PartialEq)]
pub enum FixBoundariesError {
    /// The chunk contains records but no gap of width
    /// `> max_group_span` exists between any pair of adjacent
    /// cohort-wide positions, nor between the rightmost record and
    /// `range.end`. The driver's contract is to extend the chunk's
    /// load range (up to the 4× hard cap) and retry; past the cap
    /// this becomes a hard "pathological input" error rather than
    /// silently producing a multi-GB chunk.
    #[error(
        "no safe boundary found inside chunk {chrom_id}:{range_start}..{range_end} \
         under max_group_span={max_group_span}"
    )]
    NoSafeGap {
        chrom_id: u32,
        range_start: u32,
        range_end: u32,
        max_group_span: u32,
    },

    /// `target_window_count` was zero — at least one window is
    /// required if the chunk's logical range is non-empty.
    #[error("target_window_count must be >= 1, got 0")]
    ZeroTargetWindowCount,

    /// `carryover.len()` did not match the chunk's cohort size.
    #[error("carryover length {got} does not match chunk cohort size {expected}")]
    CarryoverLengthMismatch { expected: usize, got: usize },
}

/// Choose the chunk's `safe_end`, split records past it into
/// `carryover`, and partition `[chunk.range.start, safe_end)` into
/// `target_window_count` windows whose boundaries all fall in safe
/// spots.
///
/// **Safe boundary — chunk end (`safe_end`).** A position `B` is a
/// safe chunk end iff (a) the gap to the nearest record below is
/// `> max_group_span` — so the grouper cannot reach across `B` even
/// after a chain of joins all the way to its `max_group_span` cap —
/// and (b) no allele rooted in any record below `B` reaches into
/// `[B, …)`. Rule (a) is needed at chunk boundaries because records
/// past `B` carry over to the next chunk's fresh grouper; both rules
/// together guarantee the carryover side starts as a fresh grouping
/// pass with no in-flight reach from this chunk.
///
/// **Safe boundary — internal window splits.** Internal boundaries
/// do not carry data over; the two sides of an internal split are
/// processed in the same chunk's data and see the same positions
/// across the boundary. Rule (a)'s gap-width condition is therefore
/// not needed — only rule (b) (no earlier reach crosses) is. The
/// simpler check makes the probe-and-slide placement (below) work
/// without enumerating safe gaps.
///
/// **Placement.** After `safe_end` is chosen, `target_window_count`
/// windows are placed in `[range.start, safe_end)` by:
/// 1. picking `target_window_count - 1` evenly-spaced desired
///    positions, then
/// 2. sliding each one *left* until it lands at a position satisfying
///    rule (b).
///
/// Boundaries that slide past `range.start` are dropped (the
/// resulting partition has fewer windows than requested rather than
/// erroring); duplicates are collapsed. Almost all positions on a
/// real genome are non-variant, so the slide terminates within a
/// handful of positions in practice.
///
/// **Carryover.** Records with position `>= safe_end` are *moved*
/// out of `chunk.per_sample[s]` and appended to `carryover[s]`.
/// `carryover` is assumed empty on entry — the loader cleared it.
/// On the chunk that ends at chromosome end, `safe_end == range.end`
/// and no records carry over.
pub fn fix_boundaries(
    chunk: &mut MaterialisedChunk,
    carryover: &mut [SampleColumns],
    scratch: &mut FixBoundariesScratch,
    max_group_span: u32,
    target_window_count: usize,
) -> Result<(), FixBoundariesError> {
    if target_window_count == 0 {
        return Err(FixBoundariesError::ZeroTargetWindowCount);
    }
    if carryover.len() != chunk.n_samples() {
        return Err(FixBoundariesError::CarryoverLengthMismatch {
            expected: chunk.n_samples(),
            got: carryover.len(),
        });
    }

    scratch.clear();
    chunk.windows.clear();

    // ── Step 1: cohort-wide (position, max_reach) timeline. ──
    for sample in &chunk.per_sample {
        for row_idx in 0..sample.n_records() {
            let pos = sample.position_at(row_idx);
            let reach = pos + sample.ref_span_at(row_idx) - 1;
            scratch.timeline.push((pos, reach));
        }
    }
    scratch.timeline.sort_unstable_by_key(|&(p, _)| p);
    // Merge duplicate positions by max-of-reach. Walk in place.
    let mut write = 0_usize;
    let mut read = 0_usize;
    while read < scratch.timeline.len() {
        let (pos, reach) = scratch.timeline[read];
        if write > 0 && scratch.timeline[write - 1].0 == pos {
            let prev_reach = scratch.timeline[write - 1].1;
            scratch.timeline[write - 1].1 = prev_reach.max(reach);
        } else {
            scratch.timeline[write] = (pos, reach);
            write += 1;
        }
        read += 1;
    }
    scratch.timeline.truncate(write);

    // ── Step 2: prefix-max of reaches. ──
    scratch.prefix_max_reach.reserve(scratch.timeline.len());
    let mut running_max = 0_u32;
    for &(_, r) in &scratch.timeline {
        running_max = running_max.max(r);
        scratch.prefix_max_reach.push(running_max);
    }

    // ── Step 3: pick safe_end (chunk boundary). ──
    let safe_end = pick_safe_end(chunk, scratch, max_group_span)?;

    // ── Step 4: split rows >= safe_end out of the chunk's columns
    //            into carryover. Carryover is assumed empty here.
    for (sample_idx, sample) in chunk.per_sample.iter_mut().enumerate() {
        let split_row_idx = match sample.binary_search_position(safe_end) {
            Ok(idx) | Err(idx) => idx,
        };
        sample.drain_rows_from_into(split_row_idx, &mut carryover[sample_idx]);
    }

    chunk.safe_end = safe_end;

    // ── Step 5: emit windows via probe-and-slide. ──
    emit_windows(chunk, scratch, target_window_count);

    Ok(())
}

/// Choose `safe_end`: the latest position in `[range.start, range.end]`
/// where the chunk can be cleanly cut and remaining records carried
/// over to the next chunk's fresh grouper.
fn pick_safe_end(
    chunk: &MaterialisedChunk,
    scratch: &FixBoundariesScratch,
    max_group_span: u32,
) -> Result<u32, FixBoundariesError> {
    // Empty chunk (no variant positions): the whole logical range
    // is safe.
    if scratch.timeline.is_empty() {
        return Ok(chunk.range.end);
    }

    // Case A — the gap from the rightmost record to range.end is
    // already wide enough, and no allele reaches into / past
    // range.end. safe_end = range.end (no carryover).
    let last_pos = scratch.timeline.last().unwrap().0;
    let overall_max_reach = *scratch.prefix_max_reach.last().unwrap();
    if chunk.range.end > last_pos
        && chunk.range.end - last_pos > max_group_span
        && overall_max_reach < chunk.range.end
    {
        return Ok(chunk.range.end);
    }

    // Case B — scan internal gaps right-to-left for the first safe
    // boundary (latest possible safe_end). Same gap-width +
    // reach-aware rule as the chunk-end check above.
    for i in (0..scratch.timeline.len().saturating_sub(1)).rev() {
        let (p_left, _) = scratch.timeline[i];
        let (p_right, _) = scratch.timeline[i + 1];
        if p_right - p_left > max_group_span && scratch.prefix_max_reach[i] < p_right {
            return Ok(p_right);
        }
    }

    Err(FixBoundariesError::NoSafeGap {
        chrom_id: chunk.chrom_id,
        range_start: chunk.range.start,
        range_end: chunk.range.end,
        max_group_span,
    })
}

/// Partition `[range.start, safe_end)` into `target_window_count`
/// windows by probe-and-slide. Pushes the resulting half-open ranges
/// into `chunk.windows`.
///
/// Walk-through:
/// - For `target_window_count == 1` (or `safe_end == range.start`):
///   one window covering the whole logical range.
/// - For `T > 1`: place `T - 1` desired boundaries evenly across
///   `[range.start, safe_end)`, then slide each one *left* until it
///   sits at a [`is_internal_split_safe`]-approved position. Drop
///   boundaries that slide past `range.start`; collapse duplicates.
///   The result is `≤ T` windows (the partition may shrink when
///   safe positions are sparse, but never errors).
fn emit_windows(
    chunk: &mut MaterialisedChunk,
    scratch: &FixBoundariesScratch,
    target_window_count: usize,
) {
    let range_start = chunk.range.start;
    let safe_end = chunk.safe_end;

    if safe_end <= range_start {
        return;
    }
    if target_window_count <= 1 {
        chunk.windows.push(range_start..safe_end);
        return;
    }

    let span = u64::from(safe_end - range_start);
    let t = target_window_count as u64;
    let mut chosen: Vec<u32> = Vec::with_capacity(target_window_count - 1);

    for i in 1..t {
        // Desired position for boundary `i` — even division across
        // the span. `safe_end - range_start` fits in u32, so the u64
        // arithmetic cannot overflow on any realistic chromosome.
        let desired = range_start + ((i * span) / t) as u32;
        // Slide left from `desired` to the nearest safe position.
        // The boundary must be > `range_start` (the first window
        // would otherwise be empty) and > the previous boundary
        // (no duplicates).
        let min_open = chosen.last().copied().unwrap_or(range_start);
        if let Some(b) = slide_left_to_safe(desired, min_open, safe_end, scratch) {
            chosen.push(b);
        }
    }

    let mut start = range_start;
    for &b in &chosen {
        chunk.windows.push(start..b);
        start = b;
    }
    chunk.windows.push(start..safe_end);
}

/// Slide `desired` left until it lands at a safe internal split.
/// Returns `None` if no safe position exists in `(min_open, safe_end]`
/// (the slide hit the floor without finding a non-variant +
/// no-earlier-reach position).
///
/// `min_open` is the latest known unsafe lower bound — either
/// `range_start` (for the first boundary) or the previously chosen
/// boundary (to keep boundaries strictly increasing).
fn slide_left_to_safe(
    desired: u32,
    min_open: u32,
    safe_end: u32,
    scratch: &FixBoundariesScratch,
) -> Option<u32> {
    // Clamp the starting candidate into the open interval
    // `(min_open, safe_end)`. `safe_end` itself is not a valid
    // internal boundary (the last window would be empty); the floor
    // is `min_open + 1` (boundaries must be strictly above the prior
    // window's start).
    let mut b = desired.min(safe_end.saturating_sub(1));
    while b > min_open {
        if is_internal_split_safe(b, scratch) {
            return Some(b);
        }
        b -= 1;
    }
    None
}

/// Is position `b` a safe internal split? Two O(log n) checks:
/// 1. `b` is not itself a variant position (no timeline entry at `b`),
/// 2. no earlier variant's reach extends to or past `b` —
///    `prefix_max_reach[i] < b` where `i` is the largest timeline
///    index with `position < b`.
fn is_internal_split_safe(b: u32, scratch: &FixBoundariesScratch) -> bool {
    // partition_point returns the smallest index `i` with
    // `timeline[i].0 >= b`; positions strictly less than `b` are
    // `timeline[..pp]`.
    let pp = scratch.timeline.partition_point(|&(p, _)| p < b);

    // (1) Reject if `b` is itself a variant position.
    if pp < scratch.timeline.len() && scratch.timeline[pp].0 == b {
        return false;
    }

    // (2) Check earlier reach. No earlier variants ⇒ safe by default.
    if pp == 0 {
        return true;
    }
    scratch.prefix_max_reach[pp - 1] < b
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::cohort_block::test_helpers::{
        loaded_chunk, record, ref_plus_alt, run_pre_pass,
    };

    #[test]
    fn pre_pass_empty_chunk_sets_safe_end_to_range_end() {
        let mut chunk = MaterialisedChunk::with_n_samples(2);
        chunk.chrom_id = 0;
        chunk.range = 1..1000;
        chunk.safe_end = 1000;
        let mut carry = vec![SampleColumns::empty(); 2];
        run_pre_pass(&mut chunk, &mut carry, 100).unwrap();
        assert_eq!(chunk.safe_end, 1000);
        assert_eq!(chunk.windows, vec![1..1000]);
        assert_eq!(carry[0].n_records(), 0);
        assert_eq!(carry[1].n_records(), 0);
    }

    #[test]
    fn pre_pass_no_carryover_when_range_end_is_far_past_last_record() {
        // Single variant at position 50; chunk range goes out to 5000.
        // Gap to range.end is 4950 > max_group_span=100, so safe_end =
        // range.end and there's no carryover.
        let recs = vec![vec![record(50, ref_plus_alt(2, 4))]];
        let (mut chunk, mut carry) = loaded_chunk(0, 1..5000, recs);
        run_pre_pass(&mut chunk, &mut carry, 100).unwrap();
        assert_eq!(chunk.safe_end, 5000);
        assert_eq!(chunk.windows, vec![1..5000]);
        assert_eq!(carry[0].n_records(), 0);
        assert_eq!(chunk.per_sample[0].n_records(), 1);
    }

    #[test]
    fn pre_pass_splits_tail_at_internal_safe_gap() {
        // Variants at 100, 105, then a big gap, then 600, 605, 610.
        // `range.end == 615` keeps the gap-to-range-end small (5 ≤
        // max_group_span=50), so case-A is bypassed and the pre-pass
        // commits to the internal gap (105 → 600). Records at 600+
        // become carryover.
        let s0 = vec![
            record(100, ref_plus_alt(3, 4)),
            record(105, ref_plus_alt(3, 4)),
            record(600, ref_plus_alt(3, 4)),
            record(605, ref_plus_alt(3, 4)),
            record(610, ref_plus_alt(3, 4)),
        ];
        let (mut chunk, mut carry) = loaded_chunk(0, 1..615, vec![s0]);
        run_pre_pass(&mut chunk, &mut carry, 50).unwrap();
        assert_eq!(chunk.safe_end, 600);
        assert_eq!(chunk.windows, vec![1..600]);
        assert_eq!(chunk.per_sample[0].n_records(), 2);
        assert_eq!(chunk.per_sample[0].position_at(0), 100);
        assert_eq!(chunk.per_sample[0].position_at(1), 105);
        assert_eq!(carry[0].n_records(), 3);
        assert_eq!(carry[0].position_at(0), 600);
        assert_eq!(carry[0].position_at(2), 610);
    }

    #[test]
    fn pre_pass_safe_end_equals_range_end_when_gap_to_right_is_wide() {
        // Same record layout as above but with range extended well past
        // the last record — case A (the range-end-is-far branch) must
        // pick safe_end = range.end and produce no carryover.
        let s0 = vec![
            record(100, ref_plus_alt(3, 4)),
            record(610, ref_plus_alt(3, 4)),
        ];
        let (mut chunk, mut carry) = loaded_chunk(0, 1..1000, vec![s0]);
        run_pre_pass(&mut chunk, &mut carry, 50).unwrap();
        assert_eq!(chunk.safe_end, 1000);
        assert_eq!(chunk.windows, vec![1..1000]);
        assert_eq!(carry[0].n_records(), 0);
        assert_eq!(chunk.per_sample[0].n_records(), 2);
    }

    #[test]
    fn pre_pass_returns_no_safe_gap_when_records_are_dense() {
        // Five records, all within max_group_span of each other and
        // less than max_group_span from range.end. No safe boundary
        // exists inside this chunk.
        let s0 = vec![
            record(100, ref_plus_alt(3, 4)),
            record(140, ref_plus_alt(3, 4)),
            record(180, ref_plus_alt(3, 4)),
            record(220, ref_plus_alt(3, 4)),
            record(260, ref_plus_alt(3, 4)),
        ];
        let (mut chunk, mut carry) = loaded_chunk(0, 1..300, vec![s0]);
        let err = run_pre_pass(&mut chunk, &mut carry, 100).unwrap_err();
        assert!(matches!(err, FixBoundariesError::NoSafeGap { .. }));
    }

    #[test]
    fn pre_pass_rejects_zero_target_window_count() {
        let mut chunk = MaterialisedChunk::with_n_samples(1);
        chunk.range = 1..100;
        let mut carry = vec![SampleColumns::empty()];
        let mut scratch = FixBoundariesScratch::new();
        let err = fix_boundaries(&mut chunk, &mut carry, &mut scratch, 50, 0).unwrap_err();
        assert_eq!(err, FixBoundariesError::ZeroTargetWindowCount);
    }

    // ──────────────────────────────────────────────────────────────
    // Phase B — target_window_count > 1 (probe-and-slide).
    // ──────────────────────────────────────────────────────────────

    /// Helper: run with explicit `target_window_count`.
    fn run_pre_pass_with_t(
        chunk: &mut MaterialisedChunk,
        carryover: &mut [SampleColumns],
        max_group_span: u32,
        target_window_count: usize,
    ) -> Result<(), FixBoundariesError> {
        let mut scratch = FixBoundariesScratch::new();
        fix_boundaries(
            chunk,
            carryover,
            &mut scratch,
            max_group_span,
            target_window_count,
        )
    }

    #[test]
    fn pre_pass_t2_with_no_variants_splits_at_midpoint() {
        // Empty chunk + T=2 → one boundary at the midpoint of the
        // logical range; both halves are valid windows.
        let mut chunk = MaterialisedChunk::with_n_samples(1);
        chunk.chrom_id = 0;
        chunk.range = 1..1001;
        chunk.safe_end = 1001;
        let mut carry = vec![SampleColumns::empty(); 1];
        run_pre_pass_with_t(&mut chunk, &mut carry, 100, 2).unwrap();
        assert_eq!(chunk.safe_end, 1001);
        // Midpoint = 1 + 500 = 501; the probe is non-variant +
        // no-earlier-reach (no variants exist), so it's accepted as
        // the boundary.
        assert_eq!(chunk.windows, vec![1..501, 501..1001]);
    }

    #[test]
    fn pre_pass_t4_with_no_variants_splits_evenly() {
        // Empty chunk + T=4 → three boundaries at the 1/4, 2/4, 3/4
        // marks of the logical range.
        let mut chunk = MaterialisedChunk::with_n_samples(1);
        chunk.chrom_id = 0;
        chunk.range = 1..401;
        chunk.safe_end = 401;
        let mut carry = vec![SampleColumns::empty(); 1];
        run_pre_pass_with_t(&mut chunk, &mut carry, 50, 4).unwrap();
        assert_eq!(chunk.windows, vec![1..101, 101..201, 201..301, 301..401]);
    }

    #[test]
    fn pre_pass_t2_slides_left_off_variant_to_nearest_safe_position() {
        // Variant at position 500 with ref_span=1; T=2 picks desired
        // boundary at 501 (midpoint of 1..1001). 501 is non-variant
        // and no earlier reach (variant at 500 reaches 500), so 501
        // is safe and chosen directly. To exercise the slide-left
        // path, place the desired boundary AT a variant position.
        // Range 1..1001 with variant at 501: midpoint 501 is a
        // variant → slide left to 500 (also a variant per the
        // setup) → keep sliding to 499 (non-variant, no reach
        // crosses).
        let s0 = vec![record(501, ref_plus_alt(2, 4))];
        let (mut chunk, mut carry) = loaded_chunk(0, 1..1001, vec![s0]);
        // safe_end picked via case A — gap to 1001 is 499 > 100.
        run_pre_pass_with_t(&mut chunk, &mut carry, 100, 2).unwrap();
        assert_eq!(chunk.safe_end, 1001);
        // Slide-left lands at the first non-variant position below
        // the desired (501): 500.
        assert_eq!(chunk.windows, vec![1..500, 500..1001]);
    }

    #[test]
    fn pre_pass_t2_slides_left_off_deletion_reach() {
        // Crucial test: a deletion at pos 400 with ref_span=110
        // reaches 509. T=2 picks desired boundary at 501 (midpoint
        // of 1..1001). 501 is non-variant but IS inside the
        // deletion's reach — the simple "non-variant only" rule
        // would cut here and split the group. The reach-aware
        // check rejects 501 and slides left to position 400 (the
        // variant position) → reject → 399 → check: pp finds
        // partition_point of `p < 399` → 0 (no earlier variants)
        // → safe.
        //
        // Build the deletion by hand: REF of 110bp, ALT="A" (single
        // base — net length 109 deleted bases, ref_span=110).
        let big_ref = vec![b'A'; 110];
        let s0 = vec![record(
            400,
            vec![
                crate::var_calling::cohort_block::test_helpers::allele(&big_ref, 5, -1.0, &[]),
                crate::var_calling::cohort_block::test_helpers::allele(b"A", 4, -1.0, &[]),
            ],
        )];
        let (mut chunk, mut carry) = loaded_chunk(0, 1..1001, vec![s0]);
        run_pre_pass_with_t(&mut chunk, &mut carry, 200, 2).unwrap();
        assert_eq!(chunk.safe_end, 1001);
        // Slide-left from desired=501. 501 is non-variant but
        // prefix_max_reach[idx for variants < 501] = 509 → unsafe.
        // 500..510 are all inside the deletion's reach → unsafe.
        // 400 is variant → unsafe. 399 is non-variant + no earlier
        // reach → safe.
        assert_eq!(chunk.windows, vec![1..399, 399..1001]);
    }

    #[test]
    fn pre_pass_t2_with_no_safe_internal_position_emits_single_window() {
        // Variant at position 1 with ref_span reaching all the way
        // past the midpoint: only the chunk boundary can act as a
        // split. T=2 desired boundary slides all the way to
        // range_start without finding a safe position; result is
        // a single window.
        let big_ref = vec![b'A'; 600];
        let s0 = vec![record(
            1,
            vec![
                crate::var_calling::cohort_block::test_helpers::allele(&big_ref, 5, -1.0, &[]),
                crate::var_calling::cohort_block::test_helpers::allele(b"A", 4, -1.0, &[]),
            ],
        )];
        let (mut chunk, mut carry) = loaded_chunk(0, 1..1001, vec![s0]);
        run_pre_pass_with_t(&mut chunk, &mut carry, 100, 2).unwrap();
        // Desired boundary midpoint=501 slides left:
        //   - 501..600 inside reach → unsafe.
        //   - 1..1 (range_start) → loop exits without finding a
        //     safe position.
        // Result: single window, no internal split.
        assert_eq!(chunk.windows, vec![1..1001]);
    }

    #[test]
    fn pre_pass_t_above_one_lands_boundaries_around_variant_cluster() {
        // T=4 over a variant cluster near the middle of the range.
        // Variants at 250 and 260 (ref_span=1 each); range 1..501,
        // safe_end=501 (case A — gap to 501 is 241 > max_group_span=100).
        // T=4 desired boundaries at 126, 251, 376.
        //   - 126: non-variant, no earlier reach → safe at 126.
        //   - 251: non-variant; earlier variant at 250 has reach 250
        //     < 251 → safe at 251.
        //   - 376: non-variant; max earlier reach = 260 < 376 → safe
        //     at 376.
        // Result: 4 windows, no slides needed.
        let s0 = vec![
            record(250, ref_plus_alt(3, 4)),
            record(260, ref_plus_alt(3, 4)),
        ];
        let (mut chunk, mut carry) = loaded_chunk(0, 1..501, vec![s0]);
        run_pre_pass_with_t(&mut chunk, &mut carry, 100, 4).unwrap();
        assert_eq!(chunk.windows, vec![1..126, 126..251, 251..376, 376..501]);
    }

    #[test]
    fn pre_pass_t_above_one_drops_boundaries_that_cannot_slide_above_min_open() {
        // A single wide-reach variant covers most of the chunk; T=4
        // desired boundaries at 251 / 501 / 751 all fall inside the
        // reach. The first probe slides left of the variant and
        // commits a boundary; subsequent probes cannot slide below
        // that boundary (the `min_open` floor keeps boundaries
        // strictly increasing) and are dropped — final partition
        // has fewer than T windows rather than zero-width ones.
        //
        // Variant at 200 with ref_span=600 → reach 799. Range
        // 1..1001 with safe_end=1001 (case A: gap to 1001 is 201
        // > 100). prefix_max_reach[0] = 799.
        // T=4 desired boundaries:
        //   - 251 → unsafe (reach 799 > 251); slide left to 199
        //     (199 is non-variant; no earlier variants exist).
        //   - 501 → unsafe; slide left bounded by min_open=199 →
        //     no safe position above 199 → dropped.
        //   - 751 → same → dropped.
        // Final: 2 windows.
        let big_ref = vec![b'A'; 600];
        let s0 = vec![record(
            200,
            vec![
                crate::var_calling::cohort_block::test_helpers::allele(&big_ref, 5, -1.0, &[]),
                crate::var_calling::cohort_block::test_helpers::allele(b"A", 4, -1.0, &[]),
            ],
        )];
        let (mut chunk, mut carry) = loaded_chunk(0, 1..1001, vec![s0]);
        run_pre_pass_with_t(&mut chunk, &mut carry, 100, 4).unwrap();
        assert_eq!(chunk.windows, vec![1..199, 199..1001]);
    }
}
