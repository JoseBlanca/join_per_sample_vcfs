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

    /// Phase A only supports `target_window_count == 1`. Phase B
    /// extends the pre-pass to partition `[range.start, safe_end)`
    /// into T windows whose boundaries also fall in safe gaps.
    #[error("fix_boundaries only supports target_window_count == 1 in Phase A; got {got}")]
    UnsupportedTargetWindowCount { got: usize },

    /// `carryover.len()` did not match the chunk's cohort size.
    #[error("carryover length {got} does not match chunk cohort size {expected}")]
    CarryoverLengthMismatch { expected: usize, got: usize },
}

/// Choose the chunk's `safe_end`, split records past it into
/// `carryover`, and partition `[chunk.range.start, safe_end)` into
/// `target_window_count` windows whose boundaries all fall in safe
/// gaps.
///
/// **Safe boundary.** A position `B` is safe iff (a) the gap to the
/// nearest record below is `> max_group_span` — so the grouper
/// cannot reach across `B` even after a chain of joins all the way
/// to its `max_group_span` cap — and (b) no allele rooted in any
/// record below `B` reaches into `[B, …)`. These are the same two
/// rules the chunk plan §"Pre-pass: safe-boundary placement"
/// describes; they are checked here against the post-variant-filter
/// timeline (non-variant positions are already gone, which makes
/// safe gaps easier to find).
///
/// **Phase A.** Only `target_window_count == 1` is supported; the
/// resulting partition is a single window `[range.start, safe_end)`.
/// Phase B extends this to T-1 internal boundaries placed in safe
/// gaps near evenly-spaced positions.
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
    if target_window_count != 1 {
        return Err(FixBoundariesError::UnsupportedTargetWindowCount {
            got: target_window_count,
        });
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

    // Empty chunk (no variant positions): the whole logical range
    // is safe. No carryover.
    if scratch.timeline.is_empty() {
        chunk.safe_end = chunk.range.end;
        chunk.windows.push(chunk.range.start..chunk.range.end);
        return Ok(());
    }

    // ── Step 2: prefix-max of reaches. ──
    scratch.prefix_max_reach.reserve(scratch.timeline.len());
    let mut running_max = 0_u32;
    for &(_, r) in &scratch.timeline {
        running_max = running_max.max(r);
        scratch.prefix_max_reach.push(running_max);
    }

    // ── Step 3: pick safe_end. ──
    // Case A — the gap from the rightmost record to range.end is
    // already wide enough, and no allele reaches into / past
    // range.end. safe_end = range.end (no carryover).
    let last_pos = scratch.timeline.last().unwrap().0;
    let overall_max_reach = *scratch.prefix_max_reach.last().unwrap();
    if chunk.range.end > last_pos
        && chunk.range.end - last_pos > max_group_span
        && overall_max_reach < chunk.range.end
    {
        chunk.safe_end = chunk.range.end;
        chunk.windows.push(chunk.range.start..chunk.range.end);
        return Ok(());
    }

    // Case B — scan internal gaps right-to-left for the first safe
    // boundary (latest possible safe_end).
    let mut chosen: Option<u32> = None;
    for i in (0..scratch.timeline.len().saturating_sub(1)).rev() {
        let (p_left, _) = scratch.timeline[i];
        let (p_right, _) = scratch.timeline[i + 1];
        if p_right - p_left > max_group_span && scratch.prefix_max_reach[i] < p_right {
            chosen = Some(p_right);
            break;
        }
    }

    let safe_end = match chosen {
        Some(end) => end,
        None => {
            return Err(FixBoundariesError::NoSafeGap {
                chrom_id: chunk.chrom_id,
                range_start: chunk.range.start,
                range_end: chunk.range.end,
                max_group_span,
            });
        }
    };

    // ── Step 4: split rows >= safe_end out of the chunk's columns
    //            into carryover. Carryover is assumed empty here.
    for (sample_idx, sample) in chunk.per_sample.iter_mut().enumerate() {
        let split_row_idx = match sample.binary_search_position(safe_end) {
            Ok(idx) => idx,
            Err(idx) => idx,
        };
        sample.drain_rows_from_into(split_row_idx, &mut carryover[sample_idx]);
    }

    chunk.safe_end = safe_end;
    if safe_end > chunk.range.start {
        chunk.windows.push(chunk.range.start..safe_end);
    }
    Ok(())
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
    fn pre_pass_rejects_target_window_count_above_one_in_phase_a() {
        let mut chunk = MaterialisedChunk::with_n_samples(1);
        chunk.range = 1..100;
        let mut carry = vec![SampleColumns::empty()];
        let mut scratch = FixBoundariesScratch::new();
        let err = fix_boundaries(&mut chunk, &mut carry, &mut scratch, 50, 4).unwrap_err();
        assert!(matches!(
            err,
            FixBoundariesError::UnsupportedTargetWindowCount { got: 4 }
        ));
    }
}
