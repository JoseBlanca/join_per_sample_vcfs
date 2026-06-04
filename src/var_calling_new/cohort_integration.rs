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
//! - **2b** (`CohortChunkIntegrator`): the streaming integrator — segment
//!   buffering to the watermark, per-position record building via `take_*`,
//!   REF fetch, `chunk_order` stamping, the `run` loop.

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
}
