//! Variant-group partitioning of one worker window.
//!
//! After [`load_chunk_from_iters`](super::loader::load_chunk_from_iters)
//! has filled a chunk and [`fix_boundaries`](super::pre_pass::fix_boundaries)
//! has sliced it into worker windows, this stage takes one window
//! and produces the variant-group work items the worker math will
//! iterate.
//!
//! **Input:** a `&MaterialisedChunk` + one window `Range<u32>` from
//! `chunk.windows`.
//!
//! **Output:** a [`WindowPartition`] — CSR-style columnar list of
//! variant groups. Each group carries:
//! - its 1-based inclusive `[start, end]` genomic span;
//! - the cohort positions inside it (sorted ascending);
//! - per position, the `(sample_idx, row_idx)` pointers the worker
//!   uses to read each sample's record out of `chunk.per_sample`.
//!
//! Two positions are in the same group iff their reference-span
//! footprints overlap transitively, capped at `max_group_span` —
//! the same domain rule the streaming grouper enforced, but
//! evaluated directly on the columnar chunk rather than over a
//! pull-iterator chain.
//!
//! **DUST low-complexity masking** is applied inline: the caller
//! passes a sorted, non-overlapping `&[Range<u32>]` of 1-based
//! half-open genomic intervals (the chunk-level DUST mask,
//! translated to genomic coordinates), and `partition_window`
//! skips every cohort position that falls inside one of those
//! intervals. Masked positions are advanced through (their sample
//! cursors move forward) but their records do not enter any
//! group and do not contribute to any group's reach. The running
//! group state survives masked positions — so an MNP at position
//! P that already extended the group's reach past a later
//! masked position Q still allows the group to continue past Q
//! when an unmasked position arrives. That matches what the
//! streaming `DustFilter → VariantGrouper` chain on `main` does:
//! the masked position is dropped before the grouper sees it, so
//! it never contributes to or breaks the running group.

use std::ops::Range;

use thiserror::Error;

use crate::var_calling::cohort_block::columns::{MaterialisedChunk, u32_from_usize};

/// One worker window's variant-group partitions in CSR-style
/// columnar layout. Symmetric with
/// [`SampleColumns`](super::columns::SampleColumns) — every
/// ragged dimension is encoded with an offsets column + a flat
/// values column.
///
/// **Per-group fixed-width** (length = `n_groups`):
/// - [`Self::group_starts`], [`Self::group_ends`] — 1-based
///   inclusive genomic span the group covers (REF spans accounted
///   for via reach).
///
/// **Per-group CSR into positions** (length = `n_groups + 1`):
/// - [`Self::group_position_offsets`] — group `g`'s cohort
///   positions occupy `positions[group_position_offsets[g]
///   ..group_position_offsets[g + 1]]`.
///
/// **Per-position fixed-width** (length = total positions):
/// - [`Self::positions`] — the cohort positions, sorted ascending,
///   that fall inside the union of all groups.
///
/// **Per-position CSR into sample-handle arrays** (length =
/// `positions.len() + 1`):
/// - [`Self::position_sample_offsets`] — position `p`'s
///   sample-handle pairs occupy
///   `samples_at_pos[position_sample_offsets[p]
///   ..position_sample_offsets[p + 1]]` (likewise for
///   [`Self::rows_at_pos`]).
///
/// **Per-(position, sample-with-record) fixed-width**:
/// - [`Self::samples_at_pos`] — `sample_idx` of every sample that
///   has a record at the corresponding cohort position.
/// - [`Self::rows_at_pos`] — `row_idx` inside
///   `chunk.per_sample[sample_idx]` for that record.
///
/// Samples without a record at a position are simply absent from
/// the sample-handle pairs at that position; the worker treats
/// them as "no evidence at this position".
///
/// **Lifecycle:** the driver holds one persistent `WindowPartition`
/// across every window in every chunk. [`Self::clear`] resets it
/// without freeing the underlying allocations.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowPartition {
    pub chrom_id: u32,
    pub group_starts: Vec<u32>,
    pub group_ends: Vec<u32>,
    pub group_position_offsets: Vec<u32>,
    pub positions: Vec<u32>,
    pub position_sample_offsets: Vec<u32>,
    pub samples_at_pos: Vec<u32>,
    pub rows_at_pos: Vec<u32>,
}

impl WindowPartition {
    /// Empty partition with both CSR offset columns sentineled at
    /// `[0]`, so the invariant `len(offsets) == n + 1` holds at
    /// every step.
    pub fn empty() -> Self {
        Self {
            chrom_id: 0,
            group_starts: Vec::new(),
            group_ends: Vec::new(),
            group_position_offsets: vec![0],
            positions: Vec::new(),
            position_sample_offsets: vec![0],
            samples_at_pos: Vec::new(),
            rows_at_pos: Vec::new(),
        }
    }

    /// Reset every column while preserving allocated capacity. CSR
    /// offset columns are reset to their single-`0` sentinel state.
    pub fn clear(&mut self) {
        self.chrom_id = 0;
        self.group_starts.clear();
        self.group_ends.clear();
        self.group_position_offsets.clear();
        self.group_position_offsets.push(0);
        self.positions.clear();
        self.position_sample_offsets.clear();
        self.position_sample_offsets.push(0);
        self.samples_at_pos.clear();
        self.rows_at_pos.clear();
    }

    /// Number of variant groups currently stored.
    pub fn n_groups(&self) -> usize {
        self.group_starts.len()
    }

    /// Total cohort positions across every group.
    pub fn n_positions(&self) -> usize {
        self.positions.len()
    }

    /// Position-index range inside [`Self::positions`] for group `g`.
    pub fn position_range_for_group(&self, g: usize) -> Range<usize> {
        let lo = self.group_position_offsets[g] as usize;
        let hi = self.group_position_offsets[g + 1] as usize;
        lo..hi
    }

    /// Sample-handle range inside [`Self::samples_at_pos`] /
    /// [`Self::rows_at_pos`] for the position at `position_idx`.
    pub fn sample_range_for_position(&self, position_idx: usize) -> Range<usize> {
        let lo = self.position_sample_offsets[position_idx] as usize;
        let hi = self.position_sample_offsets[position_idx + 1] as usize;
        lo..hi
    }
}

impl Default for WindowPartition {
    fn default() -> Self {
        Self::empty()
    }
}

/// Reusable scratch for one [`partition_window`] call.
///
/// **Fields.**
/// - [`Self::sample_cursors`]: one cursor per sample —
///   `chunk.per_sample[s]`'s current row index during the N-way
///   position merge across the window. Initialised at the first
///   row whose position is `>= window.start`; advanced as each
///   position is emitted.
#[derive(Debug)]
pub struct PartitionScratch {
    sample_cursors: Vec<usize>,
}

impl PartitionScratch {
    /// Build scratch sized for a cohort of `n_samples` samples.
    pub fn with_n_samples(n_samples: usize) -> Self {
        Self {
            sample_cursors: vec![0; n_samples],
        }
    }

    /// Number of samples this scratch was sized for.
    pub fn n_samples(&self) -> usize {
        self.sample_cursors.len()
    }

    /// Reset cursors to 0; capacity preserved.
    pub fn clear(&mut self) {
        for c in &mut self.sample_cursors {
            *c = 0;
        }
    }
}

/// Errors surfaced by [`partition_window`].
#[non_exhaustive]
#[derive(Error, Debug, PartialEq)]
pub enum PartitionError {
    /// `window` is empty, reversed, or escapes the chunk's
    /// `[range.start, safe_end)`.
    #[error(
        "window {window_start}..{window_end} is empty or escapes \
         chunk range {chunk_start}..{chunk_safe_end}"
    )]
    InvalidWindow {
        window_start: u32,
        window_end: u32,
        chunk_start: u32,
        chunk_safe_end: u32,
    },

    /// `scratch.n_samples()` did not match `chunk.n_samples()`.
    #[error("scratch sized for {expected} samples but chunk has {got}")]
    ScratchSampleCountMismatch { expected: usize, got: usize },

    /// A variant group would exceed `max_group_span`. The pre-pass
    /// in [`fix_boundaries`](super::pre_pass::fix_boundaries) makes
    /// this impossible for window-aligned groups — surfacing means
    /// a single-position record's REF span alone exceeds the cap,
    /// or the caller passed a `max_group_span` smaller than the
    /// pre-pass used.
    #[error(
        "variant group at {start}..{attempted_end} exceeds \
         max_group_span={max_group_span}"
    )]
    GroupTooWide {
        start: u32,
        attempted_end: u32,
        max_group_span: u32,
    },
}

/// Partition one worker window of a [`MaterialisedChunk`] into
/// variant groups.
///
/// **Algorithm.** N-way merge across the per-sample position
/// columns clipped to `window`. At each emitted cohort position,
/// every sample with a record there contributes its `ref_span`
/// to the group's running reach. A position joins the open group
/// iff `pos <= current_end` (the rolling max of `pos + ref_span -
/// 1` across the group); otherwise it closes the open group and
/// opens a new one. `attempted_end - start + 1` is checked against
/// `max_group_span`; the pre-pass guarantees window-aligned groups
/// stay under the cap, but the check defends against caller error.
///
/// **DUST.** `masked_intervals` is a sorted, non-overlapping list
/// of half-open 1-based genomic intervals to skip — typically the
/// `sdust_mask` output translated to chromosome coordinates.
/// Cohort positions inside any masked interval are advanced past
/// (sample cursors move forward) but do not enter any group or
/// affect the running reach. Pass an empty slice to disable
/// masking.
///
/// `out` is cleared at entry — callers hold one persistent
/// [`WindowPartition`] across windows. Same scratch-reuse contract
/// as [`load_chunk_from_iters`](super::loader::load_chunk_from_iters)
/// and [`fix_boundaries`](super::pre_pass::fix_boundaries).
pub fn partition_window(
    chunk: &MaterialisedChunk,
    window: &Range<u32>,
    masked_intervals: &[Range<u32>],
    max_group_span: u32,
    scratch: &mut PartitionScratch,
    out: &mut WindowPartition,
) -> Result<(), PartitionError> {
    if window.start >= window.end || window.start < chunk.range.start || window.end > chunk.safe_end
    {
        return Err(PartitionError::InvalidWindow {
            window_start: window.start,
            window_end: window.end,
            chunk_start: chunk.range.start,
            chunk_safe_end: chunk.safe_end,
        });
    }

    let n_samples = chunk.n_samples();
    if scratch.n_samples() != n_samples {
        return Err(PartitionError::ScratchSampleCountMismatch {
            expected: scratch.n_samples(),
            got: n_samples,
        });
    }

    // M19: lock the documented `masked_intervals` invariant
    // ("sorted, non-overlapping list"). The mask-cursor logic below
    // checks only the *current* interval; a malformed mask makes
    // emit/skip decisions silently wrong. `compute_dust_mask_for_chrom`
    // (the only production caller today) returns sorted+non-overlapping
    // intervals by sdust's construction; a future caller or test helper
    // that breaks the invariant trips this assertion under tests.
    debug_assert!(
        masked_intervals.windows(2).all(|w| w[0].end <= w[1].start),
        "partition_window: masked_intervals must be sorted and non-overlapping",
    );

    scratch.clear();
    out.clear();
    out.chrom_id = chunk.chrom_id;

    // Place each sample's cursor at the first row whose position is
    // `>= window.start`. `binary_search_position` returns `Ok(i)` on
    // exact hit and `Err(i)` at the insertion point — either way `i`
    // is the right starting cursor.
    for (s, sample) in chunk.per_sample.iter().enumerate() {
        let cursor = match sample.binary_search_position(window.start) {
            Ok(idx) | Err(idx) => idx,
        };
        scratch.sample_cursors[s] = cursor;
    }

    // Tracking the in-flight group:
    let mut group_open = false;
    let mut group_start: u32 = 0;
    let mut group_end: u32 = 0;

    // Parallel forward cursor over `masked_intervals` — both the
    // emitted positions and the mask are sorted ascending, so one
    // monotonic walk suffices.
    let mut mask_cursor = 0_usize;

    loop {
        // Find the next cohort position to emit — min across all
        // samples whose cursor is still inside the window.
        let mut next_pos: Option<u32> = None;
        for (s, sample) in chunk.per_sample.iter().enumerate() {
            let cursor = scratch.sample_cursors[s];
            if cursor < sample.n_records() {
                let pos = sample.position_at(cursor);
                if pos < window.end {
                    next_pos = Some(match next_pos {
                        Some(m) => m.min(pos),
                        None => pos,
                    });
                }
            }
        }
        let Some(pos) = next_pos else { break };

        // Advance the mask cursor past every interval whose right
        // boundary is `<= pos`. Mask is sorted + non-overlapping,
        // so the only candidate to test is `masked_intervals[mask_cursor]`.
        while mask_cursor < masked_intervals.len() && masked_intervals[mask_cursor].end <= pos {
            mask_cursor += 1;
        }
        let pos_masked = mask_cursor < masked_intervals.len()
            && masked_intervals[mask_cursor].start <= pos
            && pos < masked_intervals[mask_cursor].end;

        if pos_masked {
            // Drop this cohort position: advance every sample's
            // cursor past `pos` without recording handles or
            // touching the running group state. Matches what the
            // streaming pipeline's DustFilter does upstream of the
            // grouper — the masked position is invisible to the
            // grouping decision.
            for (s, sample) in chunk.per_sample.iter().enumerate() {
                let cursor = scratch.sample_cursors[s];
                if cursor < sample.n_records() && sample.position_at(cursor) == pos {
                    scratch.sample_cursors[s] = cursor + 1;
                }
            }
            continue;
        }

        // Record this position's sample handles + compute its reach
        // (max position + ref_span - 1 across contributing samples).
        let samples_at_pos_start_len = out.samples_at_pos.len();
        let mut pos_reach = pos;
        for (s, sample) in chunk.per_sample.iter().enumerate() {
            let cursor = scratch.sample_cursors[s];
            if cursor < sample.n_records() && sample.position_at(cursor) == pos {
                let ref_span = sample.ref_span_at(cursor);
                pos_reach = pos_reach.max(pos + ref_span - 1);
                out.samples_at_pos.push(u32_from_usize(s));
                out.rows_at_pos.push(u32_from_usize(cursor));
                scratch.sample_cursors[s] = cursor + 1;
            }
        }
        debug_assert!(
            out.samples_at_pos.len() > samples_at_pos_start_len,
            "position {pos} emitted with zero sample handles",
        );

        // Group-join decision: extend the current group iff `pos`
        // lies within its current reach. Otherwise close it and
        // start a new one.
        if group_open && pos <= group_end {
            let attempted_end = group_end.max(pos_reach);
            if attempted_end - group_start + 1 > max_group_span {
                return Err(PartitionError::GroupTooWide {
                    start: group_start,
                    attempted_end,
                    max_group_span,
                });
            }
            group_end = attempted_end;
        } else {
            if group_open {
                out.group_starts.push(group_start);
                out.group_ends.push(group_end);
                out.group_position_offsets
                    .push(u32_from_usize(out.positions.len()));
            }
            if pos_reach - pos + 1 > max_group_span {
                return Err(PartitionError::GroupTooWide {
                    start: pos,
                    attempted_end: pos_reach,
                    max_group_span,
                });
            }
            group_open = true;
            group_start = pos;
            group_end = pos_reach;
        }

        out.positions.push(pos);
        out.position_sample_offsets
            .push(u32_from_usize(out.samples_at_pos.len()));
    }

    if group_open {
        out.group_starts.push(group_start);
        out.group_ends.push(group_end);
        out.group_position_offsets
            .push(u32_from_usize(out.positions.len()));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup_record::PileupRecord;
    use crate::var_calling::cohort_block::test_helpers::{
        allele, loaded_chunk, record, ref_plus_alt,
    };

    fn run_partition(
        chunk: &MaterialisedChunk,
        window: &Range<u32>,
        max_group_span: u32,
    ) -> WindowPartition {
        run_partition_masked(chunk, window, &[], max_group_span)
    }

    fn run_partition_masked(
        chunk: &MaterialisedChunk,
        window: &Range<u32>,
        masked_intervals: &[Range<u32>],
        max_group_span: u32,
    ) -> WindowPartition {
        let mut scratch = PartitionScratch::with_n_samples(chunk.n_samples());
        let mut out = WindowPartition::empty();
        partition_window(
            chunk,
            window,
            masked_intervals,
            max_group_span,
            &mut scratch,
            &mut out,
        )
        .expect("partition_window succeeded");
        out
    }

    /// Materialise (sample, row) handles for the position at
    /// `partition.positions[position_idx]` into a Vec for easier
    /// assertion in tests.
    fn handles_at(partition: &WindowPartition, position_idx: usize) -> Vec<(u32, u32)> {
        let range = partition.sample_range_for_position(position_idx);
        partition.samples_at_pos[range.clone()]
            .iter()
            .zip(&partition.rows_at_pos[range])
            .map(|(&s, &r)| (s, r))
            .collect()
    }

    #[test]
    fn empty_window_produces_empty_partition() {
        let s0 = vec![record(100, ref_plus_alt(3, 4))];
        let (chunk, _) = loaded_chunk(0, 1..1000, vec![s0]);
        // Window has no records inside it.
        let partition = run_partition(&chunk, &(200..300), 100);
        assert_eq!(partition.chrom_id, 0);
        assert_eq!(partition.n_groups(), 0);
        assert_eq!(partition.n_positions(), 0);
    }

    #[test]
    fn single_position_yields_one_singleton_group() {
        let s0 = vec![record(100, ref_plus_alt(3, 4))];
        let (chunk, _) = loaded_chunk(0, 1..1000, vec![s0]);
        let partition = run_partition(&chunk, &(1..1000), 100);
        assert_eq!(partition.n_groups(), 1);
        assert_eq!(partition.group_starts, vec![100]);
        assert_eq!(partition.group_ends, vec![100]);
        assert_eq!(partition.positions, vec![100]);
        assert_eq!(partition.group_position_offsets, vec![0, 1]);
        assert_eq!(partition.position_sample_offsets, vec![0, 1]);
        assert_eq!(handles_at(&partition, 0), vec![(0, 0)]);
    }

    #[test]
    fn close_positions_join_into_one_group() {
        // Two SNPs at 100 and 102 — ref_spans of 1 each, reaches
        // 100 and 102. The grouper joins them because 102 <= 100
        // (reach of 100, after rolling) — actually no, 102 > 100.
        // For these to join the SNP at 100 needs ref_span > 2.
        // Use an MNP (ref_span=3) at 100 → reach 102 → SNP at
        // 102 joins.
        let s0 = vec![
            record(
                100,
                vec![allele(b"ACG", 3, -1.0, &[]), allele(b"AAG", 4, -1.0, &[])],
            ),
            record(102, ref_plus_alt(3, 4)),
        ];
        let (chunk, _) = loaded_chunk(0, 1..1000, vec![s0]);
        let partition = run_partition(&chunk, &(1..1000), 100);
        assert_eq!(partition.n_groups(), 1);
        assert_eq!(partition.group_starts, vec![100]);
        assert_eq!(partition.group_ends, vec![102]);
        assert_eq!(partition.positions, vec![100, 102]);
        assert_eq!(partition.group_position_offsets, vec![0, 2]);
    }

    #[test]
    fn far_positions_split_into_separate_groups() {
        let s0 = vec![
            record(100, ref_plus_alt(3, 4)),
            record(200, ref_plus_alt(3, 4)),
        ];
        let (chunk, _) = loaded_chunk(0, 1..1000, vec![s0]);
        let partition = run_partition(&chunk, &(1..1000), 50);
        assert_eq!(partition.n_groups(), 2);
        assert_eq!(partition.group_starts, vec![100, 200]);
        assert_eq!(partition.group_ends, vec![100, 200]);
        assert_eq!(partition.positions, vec![100, 200]);
        assert_eq!(partition.group_position_offsets, vec![0, 1, 2]);
    }

    #[test]
    fn multi_sample_position_carries_every_sample_handle() {
        // Two samples both with a variant at position 50; sample 0
        // also has a variant at position 80, sample 1 also at 90.
        // Window: 1..200.
        let s0 = vec![
            record(50, ref_plus_alt(5, 4)),
            record(80, ref_plus_alt(5, 4)),
        ];
        let s1 = vec![
            record(50, ref_plus_alt(6, 3)),
            record(90, ref_plus_alt(6, 3)),
        ];
        let (chunk, _) = loaded_chunk(0, 1..200, vec![s0, s1]);
        let partition = run_partition(&chunk, &(1..200), 25);
        // 50 / 80 / 90 are >25 apart pairwise → three groups.
        assert_eq!(partition.n_groups(), 3);
        assert_eq!(partition.group_starts, vec![50, 80, 90]);
        assert_eq!(partition.positions, vec![50, 80, 90]);
        // Position 50: both samples; positions 80/90: one sample each.
        assert_eq!(handles_at(&partition, 0), vec![(0, 0), (1, 0)]);
        assert_eq!(handles_at(&partition, 1), vec![(0, 1)]);
        assert_eq!(handles_at(&partition, 2), vec![(1, 1)]);
    }

    #[test]
    fn partition_clears_persisted_buffer_between_windows() {
        // Build two chunks; partition them into the same buffer.
        let s0 = vec![record(100, ref_plus_alt(3, 4))];
        let (chunk_a, _) = loaded_chunk(0, 1..1000, vec![s0]);

        let s0 = vec![
            record(200, ref_plus_alt(3, 4)),
            record(300, ref_plus_alt(3, 4)),
        ];
        let (chunk_b, _) = loaded_chunk(0, 1..1000, vec![s0]);

        let mut scratch = PartitionScratch::with_n_samples(1);
        let mut out = WindowPartition::empty();

        partition_window(&chunk_a, &(1..1000), &[], 50, &mut scratch, &mut out).unwrap();
        assert_eq!(out.n_groups(), 1);
        assert_eq!(out.positions, vec![100]);

        partition_window(&chunk_b, &(1..1000), &[], 50, &mut scratch, &mut out).unwrap();
        assert_eq!(out.n_groups(), 2);
        assert_eq!(out.positions, vec![200, 300]);
        assert_eq!(out.chrom_id, 0);
    }

    #[test]
    fn partition_window_uses_only_records_inside_window() {
        // Window 50..120 picks up only the middle of three records.
        let s0 = vec![
            record(40, ref_plus_alt(2, 3)),
            record(80, ref_plus_alt(2, 3)),
            record(160, ref_plus_alt(2, 3)),
        ];
        let (chunk, _) = loaded_chunk(0, 1..200, vec![s0]);
        let partition = run_partition(&chunk, &(50..120), 50);
        assert_eq!(partition.positions, vec![80]);
        assert_eq!(partition.n_groups(), 1);
    }

    #[test]
    fn partition_rejects_window_escaping_chunk_range() {
        let s0 = vec![record(50, ref_plus_alt(2, 3))];
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        let mut scratch = PartitionScratch::with_n_samples(1);
        let mut out = WindowPartition::empty();
        // Window starts before chunk.range.start.
        let err = partition_window(&chunk, &(0..50), &[], 50, &mut scratch, &mut out).unwrap_err();
        assert!(matches!(err, PartitionError::InvalidWindow { .. }));
        // Window ends after chunk.safe_end.
        let err =
            partition_window(&chunk, &(50..200), &[], 50, &mut scratch, &mut out).unwrap_err();
        assert!(matches!(err, PartitionError::InvalidWindow { .. }));
        // Empty window.
        let err = partition_window(&chunk, &(50..50), &[], 50, &mut scratch, &mut out).unwrap_err();
        assert!(matches!(err, PartitionError::InvalidWindow { .. }));
    }

    #[test]
    fn partition_rejects_scratch_sized_for_wrong_cohort() {
        let s0 = vec![record(50, ref_plus_alt(2, 3))];
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        let mut scratch = PartitionScratch::with_n_samples(7);
        let mut out = WindowPartition::empty();
        let err = partition_window(&chunk, &(1..100), &[], 50, &mut scratch, &mut out).unwrap_err();
        assert!(matches!(
            err,
            PartitionError::ScratchSampleCountMismatch {
                expected: 7,
                got: 1
            }
        ));
    }

    #[test]
    fn partition_rejects_group_exceeding_max_span_via_single_record() {
        // One record with a 10-bp REF span; max_group_span=5 → the
        // first record alone violates the cap.
        let long_ref = b"AAAAAAAAAA"; // 10 bp
        let long_alt = b"AAAAAAAATA"; // 10 bp
        let s0 = vec![record(
            50,
            vec![
                allele(long_ref, 3, -1.0, &[]),
                allele(long_alt, 4, -1.0, &[]),
            ],
        )];
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        let mut scratch = PartitionScratch::with_n_samples(1);
        let mut out = WindowPartition::empty();
        let err = partition_window(&chunk, &(1..100), &[], 5, &mut scratch, &mut out).unwrap_err();
        assert!(matches!(err, PartitionError::GroupTooWide { .. }));
    }

    #[test]
    fn partition_carries_chrom_id_from_chunk() {
        let s0 = vec![record(50, ref_plus_alt(2, 3))];
        // Build a chunk on chromosome 7 by editing the fixture manually
        // (loaded_chunk hard-codes chrom 0; we replace it after).
        let (mut chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        // The records inside the chunk also need to claim chrom 7 for
        // a real load — but partition_window only reads chrom_id from
        // the chunk struct, so editing it directly is sound here.
        chunk.chrom_id = 7;
        let partition = run_partition(&chunk, &(1..100), 50);
        assert_eq!(partition.chrom_id, 7);
    }

    /// Sanity: partitioning a chunk with no records at all gives an
    /// empty partition even when the window covers the whole chunk.
    #[test]
    fn partition_empty_chunk() {
        // Build the empty chunk by hand — `loaded_chunk` always
        // routes through the loader, which expects at least one
        // record per sample path to flow through.
        let mut chunk = MaterialisedChunk::with_n_samples(2);
        chunk.chrom_id = 0;
        chunk.range = 1..100;
        chunk.safe_end = 100;
        let partition = run_partition(&chunk, &(1..100), 50);
        assert_eq!(partition.n_groups(), 0);
        assert_eq!(partition.n_positions(), 0);
        // Sentinels intact.
        assert_eq!(partition.group_position_offsets, vec![0]);
        assert_eq!(partition.position_sample_offsets, vec![0]);
    }

    /// Sanity: dropping a sample with no records in the window
    /// doesn't break the cursor advancement loop.
    #[test]
    fn partition_handles_sample_with_no_records_in_window() {
        let s0 = vec![record(50, ref_plus_alt(2, 3))];
        let s1: Vec<PileupRecord> = Vec::new(); // sample 1 has nothing
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0, s1]);
        assert_eq!(chunk.per_sample[1].n_records(), 0);
        let partition = run_partition(&chunk, &(1..100), 50);
        assert_eq!(partition.n_groups(), 1);
        assert_eq!(partition.positions, vec![50]);
        assert_eq!(handles_at(&partition, 0), vec![(0, 0)]);
    }

    // ── DUST mask tests ──

    #[test]
    fn dust_mask_drops_position_inside_masked_interval() {
        let s0 = vec![
            record(50, ref_plus_alt(2, 3)),
            record(60, ref_plus_alt(2, 3)),
            record(70, ref_plus_alt(2, 3)),
        ];
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        // Mask position 60.
        let partition = run_partition_masked(&chunk, &(1..100), &[60..61], 5);
        assert_eq!(partition.positions, vec![50, 70]);
        // 50 and 70 are >5 apart (max_group_span) → two groups.
        assert_eq!(partition.n_groups(), 2);
    }

    #[test]
    fn dust_mask_drops_every_masked_position_in_interval() {
        let s0 = vec![
            record(50, ref_plus_alt(2, 3)),
            record(55, ref_plus_alt(2, 3)),
            record(58, ref_plus_alt(2, 3)),
            record(70, ref_plus_alt(2, 3)),
        ];
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        // Mask 53..60 (covers 55 and 58 but not 50 or 70).
        let partition = run_partition_masked(&chunk, &(1..100), &[53..60], 50);
        assert_eq!(partition.positions, vec![50, 70]);
    }

    #[test]
    fn dust_mask_does_nothing_when_no_position_falls_inside() {
        // Compare unmasked vs masked output — the mask interval
        // sits entirely between two positions and touches neither,
        // so output must be identical.
        let s0 = vec![
            record(50, ref_plus_alt(2, 3)),
            record(60, ref_plus_alt(2, 3)),
        ];
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        let unmasked = run_partition(&chunk, &(1..100), 50);
        let masked = run_partition_masked(&chunk, &(1..100), &[51..55], 50);
        assert_eq!(unmasked.positions, masked.positions);
        assert_eq!(unmasked.group_starts, masked.group_starts);
        assert_eq!(unmasked.group_ends, masked.group_ends);
    }

    #[test]
    fn dust_mask_with_multiple_intervals_walks_monotonically() {
        // Positions every 10 bp; mask alternates.
        let s0 = vec![
            record(10, ref_plus_alt(2, 3)),
            record(20, ref_plus_alt(2, 3)),
            record(30, ref_plus_alt(2, 3)),
            record(40, ref_plus_alt(2, 3)),
            record(50, ref_plus_alt(2, 3)),
        ];
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        // Drop 20 and 40.
        let partition = run_partition_masked(&chunk, &(1..100), &[20..21, 40..41], 50);
        assert_eq!(partition.positions, vec![10, 30, 50]);
    }

    #[test]
    fn dust_mask_covering_all_positions_yields_empty_partition() {
        let s0 = vec![
            record(50, ref_plus_alt(2, 3)),
            record(60, ref_plus_alt(2, 3)),
        ];
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        let partition = run_partition_masked(&chunk, &(1..100), &[1..100], 50);
        assert_eq!(partition.n_groups(), 0);
        assert_eq!(partition.n_positions(), 0);
        assert_eq!(partition.group_position_offsets, vec![0]);
        assert_eq!(partition.position_sample_offsets, vec![0]);
    }

    /// Boundary check: a half-open interval `[a, b)` masks `a` but
    /// not `b`. So `60..61` masks exactly position 60; `60..60` is
    /// an empty interval and masks nothing.
    #[test]
    fn dust_mask_boundaries_are_half_open() {
        let s0 = vec![
            record(60, ref_plus_alt(2, 3)),
            record(61, ref_plus_alt(2, 3)),
        ];
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        // [60, 61) → drops 60 only.
        let partition = run_partition_masked(&chunk, &(1..100), &[60..61], 50);
        assert_eq!(partition.positions, vec![61]);

        // [60, 60) → empty, drops nothing.
        let partition = run_partition_masked(&chunk, &(1..100), &[60..60], 50);
        assert_eq!(partition.positions, vec![60, 61]);
    }

    /// A masked position that would have started a *new* group must
    /// not affect the group state — the next unmasked position is
    /// what decides whether the running group continues or closes.
    #[test]
    fn dust_mask_preserves_grouping_decisions_around_dropped_position() {
        // An MNP at 100 with ref_span=20 reaches 119. A position at
        // 110 normally joins that group. If we mask 110, the group
        // still ends at 119 (set by 100's ref_span). A later position
        // at 115 still joins, even though 110 was dropped.
        let s0 = vec![
            record(
                100,
                vec![
                    allele(&[b'A'; 20], 3, -1.0, &[]),
                    allele(&[b'C'; 20], 4, -1.0, &[]),
                ],
            ),
            record(110, ref_plus_alt(2, 3)),
            record(115, ref_plus_alt(2, 3)),
        ];
        let (chunk, _) = loaded_chunk(0, 1..200, vec![s0]);
        let partition = run_partition_masked(&chunk, &(1..200), &[110..111], 50);
        assert_eq!(partition.positions, vec![100, 115]);
        assert_eq!(partition.n_groups(), 1);
        assert_eq!(partition.group_starts, vec![100]);
        assert_eq!(partition.group_ends, vec![119]);
    }
}
