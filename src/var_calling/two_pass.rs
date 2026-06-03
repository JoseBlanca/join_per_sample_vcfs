//! Low-memory two-pass producer building blocks.
//!
//! The default producer folds **all N samples' window records into
//! memory at once** (`raw_per_sample`) to decide which positions are
//! variable, then compacts. The two-pass producer avoids ever holding
//! every sample's records simultaneously:
//!
//! 1. **Pass 1 — summarise.** Read each sample's window and fold it into
//!    a compact, per-position [`WindowSummary`] (no per-sample records
//!    retained). Its size is O(distinct window positions), independent
//!    of the cohort size N. From it we derive the variant groups and the
//!    kept-position set.
//! 2. **Pass 2 — materialise.** Re-read the same window and append only
//!    the kept positions into the output chunk.
//!
//! This module owns Phase 1: the summary type, its incremental fold, and
//! the keep decision. The producer wiring (the two read passes + the cut)
//! lives in the driver/loader.
//!
//! **Byte-identity.** The keep decision uses a *position-based* upper
//! bound on the downstream `min_alt_obs_per_sample` filter:
//! `max_nonref_obs[pos]` is the max over samples of a sample's summed
//! non-REF observations at that position, so a group's summed
//! `max_nonref_obs` is `>=` any single sample's true group obs-sum. It is
//! therefore a *safe over-approximation* — it never drops a group the
//! worker's exact per-(sample, allele) filter would keep, so the emitted
//! VCF is identical (the worker drops any extra groups that slip
//! through). At `min_alt_obs <= 1` it reproduces the historical "group
//! has any observed non-REF allele" criterion exactly.

use crate::var_calling::columns::SampleColumns;

/// Reach of a position given its (max) ref span, mirroring the loader's
/// grouping arithmetic: `pos + max(span, 1) - 1`, saturating.
#[inline]
fn reach(pos: u32, span: u32) -> u32 {
    pos.saturating_add(span.max(1)).saturating_sub(1)
}

/// Compact per-position cohort summary for the two-pass producer.
///
/// Parallel sorted arrays keyed by [`positions`](Self::positions). Built
/// by [`fold_sample`](Self::fold_sample), one sample at a time, so only
/// the cohort summary (plus the one sample being merged) is ever
/// resident — never all N samples' records.
#[derive(Debug, Default)]
pub struct WindowSummary {
    /// Sorted, unique 1-based positions present in any folded sample.
    positions: Vec<u32>,
    /// Parallel to `positions`: max `ref_span` across samples with a
    /// record at the position (drives the grouping reach calculation).
    max_ref_span: Vec<u32>,
    /// Parallel to `positions`: max over samples of that sample's summed
    /// non-REF observations at the position. `0` iff no sample carries an
    /// observed non-REF allele there. See the module-level byte-identity
    /// note.
    max_nonref_obs: Vec<u32>,
    // Double buffers for the merge in `fold_sample` (kept to avoid a
    // per-call allocation).
    tmp_positions: Vec<u32>,
    tmp_max_ref_span: Vec<u32>,
    tmp_max_nonref_obs: Vec<u32>,
}

impl WindowSummary {
    /// A fresh, empty summary.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reset to empty, preserving allocated capacity for reuse across
    /// windows.
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

    /// Merge one sample's window records into the summary. `cols` must be
    /// position-sorted (the `SampleColumns` invariant). Aggregation is
    /// `max` on `ref_span` and on the per-position non-REF obs sum, and a
    /// position-union on the keys — order-independent, so folding the
    /// samples in any order yields the same summary.
    pub fn fold_sample(&mut self, cols: &SampleColumns) {
        let n = cols.n_records();
        self.tmp_positions.clear();
        self.tmp_max_ref_span.clear();
        self.tmp_max_nonref_obs.clear();

        let mut i = 0usize; // cursor into the current summary
        let mut j = 0usize; // cursor into the sample's records
        let m = self.positions.len();
        while i < m && j < n {
            let ps = self.positions[i];
            let pc = cols.position_at(j);
            if ps < pc {
                self.push_tmp(ps, self.max_ref_span[i], self.max_nonref_obs[i]);
                i += 1;
            } else if pc < ps {
                self.push_tmp(pc, cols.ref_span_at(j), cols.non_ref_obs_sum_at(j));
                j += 1;
            } else {
                self.push_tmp(
                    ps,
                    self.max_ref_span[i].max(cols.ref_span_at(j)),
                    self.max_nonref_obs[i].max(cols.non_ref_obs_sum_at(j)),
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
            self.push_tmp(
                cols.position_at(j),
                cols.ref_span_at(j),
                cols.non_ref_obs_sum_at(j),
            );
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

    /// Derive `is_kept` (parallel to [`positions`](Self::positions)):
    /// `true` for every position in a variant group whose summed
    /// `max_nonref_obs` reaches `min_alt_obs`. Groups are formed by the
    /// same overlapping-reach rule the loader uses, and kept or dropped
    /// whole (so multi-position groups keep their pure-REF-in-isolation
    /// positions, preserving REF evidence for homref samples).
    ///
    /// `min_alt_obs <= 1` is the historical "group has any observed
    /// non-REF allele" criterion exactly; `> 1` is the byte-identity-safe
    /// ALT-count over-approximation (see the module note).
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

    /// Merge another summary into `self` (`max` aggregation, position
    /// union) — the reduce step when pass 1 summarises samples in
    /// parallel into per-thread partials. `other` is consumed
    /// conceptually (left untouched here for borrow simplicity).
    pub fn merge(&mut self, other: &WindowSummary) {
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

    /// The clean group boundary at or before `watermark`: the start of
    /// the still-open group (the last group whose reach extends past
    /// `watermark`), or `watermark + 1` if every group is closed.
    /// Mirrors the loader's `find_block_cut` so chunk cuts never split a
    /// variant group. Used to carve an interval into worker-chunks.
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

    /// Count kept (variable) positions with `position < cut`, using the
    /// same group/threshold rule as [`derive_is_kept`](Self::derive_is_kept).
    /// `cut` must land on a clean group boundary (e.g. from
    /// [`find_cut`](Self::find_cut)) so no group straddles it.
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

    /// Partition the interval `[interval_start, interval_end_exclusive)`
    /// into worker-chunk boundaries: each chunk accumulates ~`target_variants`
    /// kept positions and ends at the start of the next variant group (a
    /// clean boundary — no group straddles a chunk, so the worker's
    /// per-group output is independent of where the chunks fall, and the
    /// VCF stays byte-identical to any other partitioning).
    ///
    /// Returns `[interval_start, c1, c2, …, interval_end_exclusive]`;
    /// chunk `k` spans `[cuts[k], cuts[k+1])`. Always at least
    /// `[interval_start, interval_end_exclusive]` (one chunk).
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::test_helpers::{allele, record};

    /// Build a one-sample `SampleColumns` from `(pos, alleles)` rows.
    fn columns(rows: Vec<(u32, Vec<crate::pileup_record::AlleleObservation>)>) -> SampleColumns {
        let mut c = SampleColumns::empty();
        for (pos, alleles) in rows {
            c.push_record(record(pos, alleles));
        }
        c
    }

    #[test]
    fn fold_unions_positions_and_maxes_aggregates() {
        // Sample A: pos 10 (REF only), pos 20 (REF + ALT obs 1).
        let a = columns(vec![
            (10, vec![allele(b"A", 3, -1.0, &[])]),
            (
                20,
                vec![allele(b"A", 2, -1.0, &[]), allele(b"T", 1, -1.0, &[])],
            ),
        ]);
        // Sample B: pos 20 (REF + ALT obs 2), pos 30 (REF only).
        let b = columns(vec![
            (
                20,
                vec![allele(b"A", 1, -1.0, &[]), allele(b"T", 2, -1.0, &[])],
            ),
            (30, vec![allele(b"A", 4, -1.0, &[])]),
        ]);

        let mut s = WindowSummary::new();
        s.fold_sample(&a);
        s.fold_sample(&b);

        assert_eq!(s.positions(), &[10, 20, 30]);
        // max_nonref_obs: pos10=0, pos20=max(1,2)=2, pos30=0.
        assert_eq!(s.max_nonref_obs, vec![0, 2, 0]);
        // Folding order must not matter.
        let mut s2 = WindowSummary::new();
        s2.fold_sample(&b);
        s2.fold_sample(&a);
        assert_eq!(s2.positions(), s.positions());
        assert_eq!(s2.max_nonref_obs, s.max_nonref_obs);
        assert_eq!(s2.max_ref_span, s.max_ref_span);
    }

    #[test]
    fn keep_threshold_one_is_variant_filter() {
        // Three isolated SNP positions (ref_span 1, so each is its own
        // group): pos10 obs0, pos20 obs1, pos30 obs0.
        let a = columns(vec![
            (10, vec![allele(b"A", 3, -1.0, &[])]),
            (
                20,
                vec![allele(b"A", 2, -1.0, &[]), allele(b"T", 1, -1.0, &[])],
            ),
            (30, vec![allele(b"A", 4, -1.0, &[])]),
        ]);
        let mut s = WindowSummary::new();
        s.fold_sample(&a);
        let mut kept = Vec::new();
        s.derive_is_kept(1, &mut kept);
        // Only the variant position (pos20) survives.
        assert_eq!(kept, vec![false, true, false]);
    }

    #[test]
    fn keep_threshold_two_drops_singletons_keeps_doubletons() {
        // pos20 has a single ALT obs (singleton) → dropped at threshold 2.
        // pos30 has two ALT obs → kept.
        let a = columns(vec![
            (
                20,
                vec![allele(b"A", 2, -1.0, &[]), allele(b"T", 1, -1.0, &[])],
            ),
            (
                30,
                vec![allele(b"A", 2, -1.0, &[]), allele(b"T", 2, -1.0, &[])],
            ),
        ]);
        let mut s = WindowSummary::new();
        s.fold_sample(&a);
        let mut kept = Vec::new();
        s.derive_is_kept(2, &mut kept);
        assert_eq!(kept, vec![false, true]);
    }

    #[test]
    fn merge_matches_sequential_fold() {
        let a = columns(vec![
            (10, vec![allele(b"A", 3, -1.0, &[])]),
            (
                20,
                vec![allele(b"A", 2, -1.0, &[]), allele(b"T", 1, -1.0, &[])],
            ),
        ]);
        let b = columns(vec![
            (
                20,
                vec![allele(b"A", 1, -1.0, &[]), allele(b"T", 2, -1.0, &[])],
            ),
            (30, vec![allele(b"A", 4, -1.0, &[])]),
        ]);
        // Sequential fold of both into one summary.
        let mut seq = WindowSummary::new();
        seq.fold_sample(&a);
        seq.fold_sample(&b);
        // Two per-thread partials reduced via merge.
        let mut pa = WindowSummary::new();
        pa.fold_sample(&a);
        let mut pb = WindowSummary::new();
        pb.fold_sample(&b);
        pa.merge(&pb);
        assert_eq!(pa.positions(), seq.positions());
        assert_eq!(pa.max_ref_span, seq.max_ref_span);
        assert_eq!(pa.max_nonref_obs, seq.max_nonref_obs);
    }

    #[test]
    fn find_cut_and_count_kept_below() {
        // Two isolated variant SNPs (ref_span 1) at 10 and 100, plus a
        // non-variant at 50. With watermark = 60, the group at 100 is
        // still open (reach 100 > 60), so the cut is its start (100).
        let a = columns(vec![
            (
                10,
                vec![allele(b"A", 2, -1.0, &[]), allele(b"T", 2, -1.0, &[])],
            ),
            (50, vec![allele(b"A", 3, -1.0, &[])]),
            (
                100,
                vec![allele(b"A", 2, -1.0, &[]), allele(b"T", 2, -1.0, &[])],
            ),
        ]);
        let mut s = WindowSummary::new();
        s.fold_sample(&a);
        assert_eq!(s.find_cut(60), 100);
        // Below cut=100: only pos 10 is a kept variant (pos 50 is its own
        // non-variant group, dropped).
        assert_eq!(s.count_kept_below(100, 1), 1);
        // Below cut covering everything: pos 10 and pos 100 kept.
        assert_eq!(s.count_kept_below(200, 1), 2);
    }

    #[test]
    fn over_approximation_sums_across_two_samples_in_a_group() {
        // Two samples each contribute one ALT obs at the SAME position;
        // max (not sum across samples) = 1, so at threshold 2 the
        // position-based bound keeps it only if a single sample reaches 2.
        // Here neither does → dropped (matches the worker's per-sample
        // filter, which also requires a single sample to reach 2).
        let a = columns(vec![(
            20,
            vec![allele(b"A", 2, -1.0, &[]), allele(b"T", 1, -1.0, &[])],
        )]);
        let b = columns(vec![(
            20,
            vec![allele(b"A", 2, -1.0, &[]), allele(b"T", 1, -1.0, &[])],
        )]);
        let mut s = WindowSummary::new();
        s.fold_sample(&a);
        s.fold_sample(&b);
        let mut kept = Vec::new();
        s.derive_is_kept(2, &mut kept);
        assert_eq!(kept, vec![false]);
    }
}
