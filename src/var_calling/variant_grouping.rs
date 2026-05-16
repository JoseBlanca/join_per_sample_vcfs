//! Stage 4 — variant grouping over the per-position pileup stream.
//!
//! Given an upstream iterator that yields `Result<PerPositionPileups,
//! MergerError>` items in strictly increasing `(chrom_id, pos)` order
//! (the contract `PerPositionMerger` provides), the [`VariantGrouper`]
//! bundles positions whose REF spans overlap into a single
//! [`OverlappingVarGroup`] and emits a stream of independent groups.
//!
//! Pure-REF positions across every sample with no group currently
//! open are dropped at the iterator boundary without ever being
//! folded into a group — see the "seed-time filter" rationale in
//! `doc/devel/implementation_plans/cohort_variant_grouping.md`. The
//! `max_var_group_span` cap in [`GrouperConfig`] hard-bounds the
//! reach of any emitted group; exceeding it is a
//! [`GrouperError::VarGroupTooWide`] rather than a silent drop.
//!
//! The grouper is generic on the upstream iterator type so synthetic
//! `std::iter`-based fixtures can drive it in unit tests without
//! standing up a real merger.

use thiserror::Error;

use crate::var_calling::per_position_merger::{MergerError, PerPositionPileups};

/// Default value for [`GrouperConfig::max_var_group_span`].
///
/// 2× Stage 1's `MAX_RECORD_SPAN` (5000 bp). High enough that real
/// cohorts have headroom for cross-sample transitive extension
/// chains, low enough to bound memory on adversarial inputs.
pub const DEFAULT_MAX_VAR_GROUP_SPAN: u32 = 10_000;

/// Tunable knobs for the grouper. Follows the per-stage config
/// convention described in
/// `doc/devel/specs/calling_pipeline_architecture.md` §"Configurable
/// parameters".
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct GrouperConfig {
    /// Hard cap on the reference span of any emitted
    /// [`OverlappingVarGroup`]. A group whose `end - start + 1`
    /// would exceed this value triggers
    /// [`GrouperError::VarGroupTooWide`] and the iterator latches
    /// `done`. Defaults to [`DEFAULT_MAX_VAR_GROUP_SPAN`].
    ///
    /// The bound is defensive — Stage 1's `MAX_RECORD_SPAN` already
    /// caps every single record's reach, so in normal operation
    /// the cap is never reached. Users hitting it on real data
    /// should raise it explicitly via the config rather than rely
    /// on the grouper to silently drop the group.
    pub max_var_group_span: u32,
}

impl Default for GrouperConfig {
    fn default() -> Self {
        Self {
            max_var_group_span: DEFAULT_MAX_VAR_GROUP_SPAN,
        }
    }
}

/// One emitted item: a single overlap-bundle of per-position pileups,
/// in strictly increasing `pos` order. `start` and `end` are 1-based
/// inclusive and cover the union of every record's REF span.
#[derive(Debug, Clone, PartialEq)]
pub struct OverlappingVarGroup {
    pub chrom_id: u32,
    pub start: u32,
    pub end: u32,
    pub records: Vec<PerPositionPileups>,
}

impl OverlappingVarGroup {
    pub fn span(&self) -> (u32, u32, u32) {
        (self.chrom_id, self.start, self.end)
    }
}

/// Errors the grouper can emit. Upstream errors come through wrapped
/// in [`GrouperError::Upstream`]; grouper-side failures (currently
/// just [`GrouperError::VarGroupTooWide`]) are their own variants.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum GrouperError {
    #[error("upstream: {0}")]
    Upstream(#[from] MergerError),

    /// A group's reference span would exceed
    /// [`GrouperConfig::max_var_group_span`]. The locus is reported
    /// so the user can either investigate it or raise the cap.
    #[error(
        "variant group at chrom {chrom_id} starting at {start} would span \
         {attempted_span} bp (>{cap} bp cap); raise GrouperConfig::max_var_group_span \
         or investigate the locus"
    )]
    VarGroupTooWide {
        chrom_id: u32,
        start: u32,
        attempted_end: u32,
        attempted_span: u32,
        cap: u32,
    },
}

/// Streaming single-pass overlap bundler over an upstream iterator
/// of per-position pileups.
pub struct VariantGrouper<I>
where
    I: Iterator<Item = Result<PerPositionPileups, MergerError>>,
{
    upstream: I,
    config: GrouperConfig,
    /// First item of the next group, already pulled from upstream
    /// but not yet folded into a group. `None` before the first
    /// emission and after exhaustion.
    pending_seed: Option<PerPositionPileups>,
    /// Latch: once set, all subsequent `next()` calls return `None`.
    done: bool,
}

impl<I> std::fmt::Debug for VariantGrouper<I>
where
    I: Iterator<Item = Result<PerPositionPileups, MergerError>>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Exhaustive destructure so a new field on `VariantGrouper`
        // fails to compile here instead of silently being dropped
        // from the debug rendering.
        let Self {
            upstream: _,
            config,
            pending_seed,
            done,
        } = self;
        f.debug_struct("VariantGrouper")
            .field("config", config)
            .field(
                "pending_seed_at",
                &pending_seed.as_ref().map(|pp| (pp.chrom_id, pp.pos)),
            )
            .field("done", done)
            .finish()
    }
}

impl<I> VariantGrouper<I>
where
    I: Iterator<Item = Result<PerPositionPileups, MergerError>>,
{
    /// Construct a grouper with explicit tuning. Pass
    /// [`GrouperConfig::default()`] for the standard defaults; pass an
    /// explicit value for [`GrouperConfig::max_var_group_span`] to
    /// override.
    pub fn with_config(upstream: I, config: GrouperConfig) -> Self {
        Self {
            upstream,
            config,
            pending_seed: None,
            done: false,
        }
    }

    pub fn config(&self) -> &GrouperConfig {
        &self.config
    }

    /// Phase A — pull from upstream (or pop `pending_seed`) until
    /// we either exhaust, hit an error, or find a variant-bearing
    /// position. Pure-REF positions are dropped here without
    /// allocation.
    fn next_seed(&mut self) -> Result<Option<PerPositionPileups>, GrouperError> {
        loop {
            let candidate = match self.pending_seed.take() {
                Some(pp) => pp,
                None => match self.upstream.next() {
                    None => return Ok(None),
                    Some(Err(e)) => return Err(GrouperError::Upstream(e)),
                    Some(Ok(pp)) => pp,
                },
            };
            if has_variant_observation(&candidate) {
                return Ok(Some(candidate));
            }
            // pure-REF with no open group; drop and keep looking.
        }
    }
}

impl<I> Iterator for VariantGrouper<I>
where
    I: Iterator<Item = Result<PerPositionPileups, MergerError>>,
{
    type Item = Result<OverlappingVarGroup, GrouperError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        // Phase A — find a variant-bearing seed.
        let seed = match self.next_seed() {
            Ok(Some(pp)) => pp,
            Ok(None) => {
                self.done = true;
                return None;
            }
            Err(e) => {
                self.done = true;
                return Some(Err(e));
            }
        };

        // Initialise the in-progress group.
        let chrom_id = seed.chrom_id;
        let start = seed.pos;
        // `max_ref_span` falls back to `1` when every slot is `None`
        // (in production the seed-time variant check screens that
        // case out); `saturating_sub` guards a hypothetical future
        // caller that exposes `ref_span = 0` from blowing up here.
        let seed_end = start + max_ref_span(&seed).saturating_sub(1);

        // Seed-time cap check — guards against a single record's
        // ref_span alone exceeding the cap. Stage 1's MAX_RECORD_SPAN
        // should prevent this in normal operation, but the cap
        // should not rely on that invariant silently.
        if seed_end - start + 1 > self.config.max_var_group_span {
            self.done = true;
            return Some(Err(GrouperError::VarGroupTooWide {
                chrom_id,
                start,
                attempted_end: seed_end,
                attempted_span: seed_end - start + 1,
                cap: self.config.max_var_group_span,
            }));
        }

        let mut end = seed_end;
        let mut records: Vec<PerPositionPileups> = vec![seed];

        // Phase B — extend until a position falls outside the
        // open group's reach (or upstream ends / errors).
        loop {
            match self.upstream.next() {
                None => {
                    self.done = true;
                    return Some(Ok(OverlappingVarGroup {
                        chrom_id,
                        start,
                        end,
                        records,
                    }));
                }
                Some(Err(e)) => {
                    // Drop the in-progress group — partial groups
                    // are never emitted because Stage 5 has no way
                    // to know the group was truncated.
                    self.done = true;
                    return Some(Err(GrouperError::Upstream(e)));
                }
                Some(Ok(pp)) => {
                    if pp.chrom_id != chrom_id || pp.pos > end {
                        // The just-pulled item starts a new group.
                        // Stash it and re-enter Phase A on the next
                        // call (possibly dropping it if it is itself
                        // pure-REF).
                        self.pending_seed = Some(pp);
                        return Some(Ok(OverlappingVarGroup {
                            chrom_id,
                            start,
                            end,
                            records,
                        }));
                    }

                    let pp_ref_span = max_ref_span(&pp);
                    let attempted_end = end.max(pp.pos + pp_ref_span - 1);
                    if attempted_end - start + 1 > self.config.max_var_group_span {
                        self.done = true;
                        return Some(Err(GrouperError::VarGroupTooWide {
                            chrom_id,
                            start,
                            attempted_end,
                            attempted_span: attempted_end - start + 1,
                            cap: self.config.max_var_group_span,
                        }));
                    }
                    end = attempted_end;
                    records.push(pp);
                }
            }
        }
    }
}

/// True if at least one `Some` slot of `pp` has more than one
/// distinct allele (i.e. some sample observed a non-REF allele at
/// this position). Walks slots until the first hit and short-circuits.
fn has_variant_observation(pp: &PerPositionPileups) -> bool {
    pp.per_sample.iter().flatten().any(|r| r.alleles.len() > 1)
}

/// Maximum `ref_span` across all `Some` slots of `pp`. Falls back to
/// `1` if every slot is `None`.
///
/// PANIC-FREE: callers gate this behind `has_variant_observation`,
/// which short-circuits in Phase A on any pp with no `Some` slots, so
/// the fallback branch is unreachable in production. It is kept for
/// future direct callers and to keep the `pos + ref_span - 1`
/// inclusive-end arithmetic well-defined on a hypothetical empty pp.
fn max_ref_span(pp: &PerPositionPileups) -> u32 {
    pp.per_sample
        .iter()
        .flatten()
        .map(|r| r.ref_span())
        .max()
        .unwrap_or(1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::per_sample_pileup::pileup::{AlleleObservation, AlleleSupportStats, PileupRecord};

    type Item = Result<PerPositionPileups, MergerError>;
    type TestIter = std::vec::IntoIter<Item>;

    // ---------- fixture builders ----------

    fn allele(seq: &[u8]) -> AlleleObservation {
        AlleleObservation::new(seq.to_vec(), AlleleSupportStats::default(), Vec::new())
    }

    fn ref_record(chrom_id: u32, pos: u32, ref_seq: &[u8]) -> PileupRecord {
        PileupRecord::new(chrom_id, pos, vec![allele(ref_seq)])
    }

    fn snp_record(chrom_id: u32, pos: u32) -> PileupRecord {
        // REF "A" + ALT "T" — variant with ref_span = 1.
        PileupRecord::new(chrom_id, pos, vec![allele(b"A"), allele(b"T")])
    }

    fn del_record(chrom_id: u32, pos: u32, ref_span: u32) -> PileupRecord {
        // REF anchor + (ref_span-1) deleted bases, ALT = anchor only.
        let ref_seq: Vec<u8> = std::iter::repeat_n(b'A', ref_span as usize).collect();
        PileupRecord::new(chrom_id, pos, vec![allele(&ref_seq), allele(b"A")])
    }

    fn pp_with(
        chrom_id: u32,
        pos: u32,
        n_samples: usize,
        slot: usize,
        record: PileupRecord,
    ) -> PerPositionPileups {
        let mut per_sample: Vec<Option<PileupRecord>> = (0..n_samples).map(|_| None).collect();
        per_sample[slot] = Some(record);
        PerPositionPileups {
            chrom_id,
            pos,
            per_sample,
        }
    }

    fn ref_only_pp(chrom_id: u32, pos: u32, n_samples: usize, slot: usize) -> PerPositionPileups {
        pp_with(
            chrom_id,
            pos,
            n_samples,
            slot,
            ref_record(chrom_id, pos, b"A"),
        )
    }

    fn snp_pp(chrom_id: u32, pos: u32, n_samples: usize, slot: usize) -> PerPositionPileups {
        pp_with(chrom_id, pos, n_samples, slot, snp_record(chrom_id, pos))
    }

    fn del_pp(
        chrom_id: u32,
        pos: u32,
        n_samples: usize,
        slot: usize,
        ref_span: u32,
    ) -> PerPositionPileups {
        pp_with(
            chrom_id,
            pos,
            n_samples,
            slot,
            del_record(chrom_id, pos, ref_span),
        )
    }

    fn empty_pp(chrom_id: u32, pos: u32, n_samples: usize) -> PerPositionPileups {
        PerPositionPileups {
            chrom_id,
            pos,
            per_sample: (0..n_samples).map(|_| None).collect(),
        }
    }

    fn iter_from(items: Vec<Item>) -> TestIter {
        items.into_iter()
    }

    fn ok_iter(pps: Vec<PerPositionPileups>) -> TestIter {
        iter_from(pps.into_iter().map(Ok).collect())
    }

    fn collect_groups<I>(
        grouper: VariantGrouper<I>,
    ) -> Result<Vec<OverlappingVarGroup>, GrouperError>
    where
        I: Iterator<Item = Result<PerPositionPileups, MergerError>>,
    {
        grouper.collect()
    }

    // ---------- empty / pass-through cases ----------

    #[test]
    fn empty_upstream_yields_immediately() {
        let mut grouper =
            VariantGrouper::with_config(iter_from(Vec::new()), GrouperConfig::default());
        assert!(grouper.next().is_none());
        // latches:
        assert!(grouper.next().is_none());
    }

    #[test]
    fn single_trivial_snp() {
        let grouper = VariantGrouper::with_config(
            ok_iter(vec![snp_pp(0, 100, 1, 0)]),
            GrouperConfig::default(),
        );
        let groups = collect_groups(grouper).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].chrom_id, 0);
        assert_eq!(groups[0].start, 100);
        assert_eq!(groups[0].end, 100);
        assert_eq!(groups[0].records.len(), 1);
    }

    // ---------- pure-REF filtering ----------

    #[test]
    fn pure_ref_only_runs_dropped_then_snp_emitted() {
        let upstream = ok_iter(vec![
            ref_only_pp(0, 100, 1, 0),
            ref_only_pp(0, 101, 1, 0),
            ref_only_pp(0, 102, 1, 0),
            snp_pp(0, 103, 1, 0),
        ]);
        let grouper = VariantGrouper::with_config(upstream, GrouperConfig::default());
        let groups = collect_groups(grouper).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].span(), (0, 103, 103));
        assert_eq!(groups[0].records.len(), 1);
    }

    #[test]
    fn long_pure_ref_run_followed_by_one_snp() {
        // 1000 pure-REF positions, then one SNP. Smoke-test that the
        // seed-search loop does not accumulate per-dropped-position
        // state.
        let mut pps: Vec<PerPositionPileups> =
            (100..1100).map(|p| ref_only_pp(0, p, 1, 0)).collect();
        pps.push(snp_pp(0, 1100, 1, 0));
        let grouper = VariantGrouper::with_config(ok_iter(pps), GrouperConfig::default());
        let groups = collect_groups(grouper).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].span(), (0, 1100, 1100));
    }

    #[test]
    fn ref_only_positions_drawn_into_real_group_are_kept() {
        // Deletion at p=100 with ref_span=5 covers p=100..104.
        // Pure-REF records at p=101, p=102 in another sample fall
        // inside that span — they must be folded into the group,
        // *not* dropped.
        let upstream = ok_iter(vec![
            del_pp(0, 100, 2, 0, 5),
            ref_only_pp(0, 101, 2, 1),
            ref_only_pp(0, 102, 2, 1),
            snp_pp(0, 104, 2, 1),
        ]);
        let grouper = VariantGrouper::with_config(upstream, GrouperConfig::default());
        let groups = collect_groups(grouper).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].span(), (0, 100, 104));
        assert_eq!(groups[0].records.len(), 4);
    }

    #[test]
    fn mixed_ref_only_runs_between_variant_groups() {
        // REF-only run, variant group, REF-only run, variant group.
        let upstream = ok_iter(vec![
            ref_only_pp(0, 100, 1, 0),
            ref_only_pp(0, 101, 1, 0),
            snp_pp(0, 102, 1, 0),
            ref_only_pp(0, 103, 1, 0),
            ref_only_pp(0, 104, 1, 0),
            snp_pp(0, 105, 1, 0),
        ]);
        let grouper = VariantGrouper::with_config(upstream, GrouperConfig::default());
        let groups = collect_groups(grouper).unwrap();
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].span(), (0, 102, 102));
        assert_eq!(groups[1].span(), (0, 105, 105));
    }

    #[test]
    fn pending_seed_re_evaluated_drops_pure_ref() {
        // SNP at 100, then SNP at 200 closes the first group; but a
        // pure-REF item at 150 between them would land in
        // pending_seed and must be dropped by Phase A re-entry.
        // Equivalent end-to-end sequence:
        let upstream = ok_iter(vec![
            snp_pp(0, 100, 1, 0),
            ref_only_pp(0, 150, 1, 0),
            snp_pp(0, 200, 1, 0),
        ]);
        let grouper = VariantGrouper::with_config(upstream, GrouperConfig::default());
        let groups = collect_groups(grouper).unwrap();
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].span(), (0, 100, 100));
        assert_eq!(groups[1].span(), (0, 200, 200));
    }

    // ---------- overlap and extension ----------

    #[test]
    fn deletion_draws_in_downstream_snp() {
        // Sample 0 has a deletion at p=100 with ref_span=5 (covers
        // p=100..104). Sample 1 has a SNP at p=103.
        let upstream = ok_iter(vec![del_pp(0, 100, 2, 0, 5), snp_pp(0, 103, 2, 1)]);
        let grouper = VariantGrouper::with_config(upstream, GrouperConfig::default());
        let groups = collect_groups(grouper).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].span(), (0, 100, 104));
        assert_eq!(groups[0].records.len(), 2);
        // Per-sample slot vectors come through verbatim.
        assert!(groups[0].records[0].per_sample[0].is_some());
        assert!(groups[0].records[0].per_sample[1].is_none());
        assert!(groups[0].records[1].per_sample[0].is_none());
        assert!(groups[0].records[1].per_sample[1].is_some());
    }

    #[test]
    fn transitive_chain_extends_through_multiple_records() {
        // The merger guarantees strict (chrom_id, pos) monotonicity,
        // so every upstream item sits at a distinct position; the
        // chain has to be carried by deletions because SNPs alone
        // have ref_span=1 and cannot overlap distinct positions.
        // Del span-4 at 100 → end=103.
        // SNP   at 101 (≤103) → end stays 103.
        // Del span-4 at 102 (≤103) → end = max(103, 102+4-1)=105.
        // SNP   at 105 (≤105) → end stays 105.
        let upstream = ok_iter(vec![
            del_pp(0, 100, 1, 0, 4),
            snp_pp(0, 101, 1, 0),
            del_pp(0, 102, 1, 0, 4),
            snp_pp(0, 105, 1, 0),
        ]);
        let groups = collect_groups(VariantGrouper::with_config(
            upstream,
            GrouperConfig::default(),
        ))
        .unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].span(), (0, 100, 105));
        assert_eq!(groups[0].records.len(), 4);
    }

    #[test]
    fn end_does_not_shrink_on_late_short_extension() {
        // Seed at p=100 with ref_span=10 (end=109). A later
        // record at p=105 with ref_span=1 must not shrink end
        // back to 105.
        let upstream = ok_iter(vec![del_pp(0, 100, 1, 0, 10), snp_pp(0, 105, 1, 0)]);
        let grouper = VariantGrouper::with_config(upstream, GrouperConfig::default());
        let groups = collect_groups(grouper).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].span(), (0, 100, 109));
    }

    #[test]
    fn single_record_end_is_inclusive_one_based() {
        // ref_span = 3 at p=100 → end must be 102, not 103.
        let upstream = ok_iter(vec![del_pp(0, 100, 1, 0, 3)]);
        let groups = collect_groups(VariantGrouper::with_config(
            upstream,
            GrouperConfig::default(),
        ))
        .unwrap();
        assert_eq!(groups[0].span(), (0, 100, 102));
    }

    // ---------- chromosome boundaries ----------

    #[test]
    fn group_never_spans_chromosomes() {
        // Deletion at end of chrom 0 with span 100, then SNP at
        // pos=1 on chrom 1. Numerically pos=1 ≤ end=199 would
        // overlap, but the chrom change forces a close.
        let upstream = ok_iter(vec![del_pp(0, 100, 1, 0, 100), snp_pp(1, 1, 1, 0)]);
        let grouper = VariantGrouper::with_config(upstream, GrouperConfig::default());
        let groups = collect_groups(grouper).unwrap();
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].span(), (0, 100, 199));
        assert_eq!(groups[1].span(), (1, 1, 1));
    }

    #[test]
    fn multi_chromosome_with_no_overlap() {
        let upstream = ok_iter(vec![
            snp_pp(0, 100, 1, 0),
            snp_pp(0, 200, 1, 0),
            snp_pp(1, 50, 1, 0),
            snp_pp(1, 150, 1, 0),
        ]);
        let grouper = VariantGrouper::with_config(upstream, GrouperConfig::default());
        let groups = collect_groups(grouper).unwrap();
        assert_eq!(groups.len(), 4);
        assert_eq!(groups[0].span(), (0, 100, 100));
        assert_eq!(groups[1].span(), (0, 200, 200));
        assert_eq!(groups[2].span(), (1, 50, 50));
        assert_eq!(groups[3].span(), (1, 150, 150));
    }

    // ---------- error propagation ----------

    fn fake_upstream_err(sample_idx: usize) -> MergerError {
        MergerError::OutOfOrder {
            sample_idx,
            sample_name: format!("S{sample_idx}"),
            chrom_id: 0,
            pos: 0,
        }
    }

    #[test]
    fn upstream_error_before_any_group() {
        let items: Vec<Item> = vec![Err(fake_upstream_err(0))];
        let mut grouper = VariantGrouper::with_config(iter_from(items), GrouperConfig::default());
        match grouper.next() {
            Some(Err(GrouperError::Upstream(MergerError::OutOfOrder { .. }))) => {}
            other => panic!("expected upstream out-of-order, got {other:?}"),
        }
        assert!(grouper.next().is_none());
    }

    #[test]
    fn upstream_error_mid_group_drops_in_progress_group() {
        // Build a group that is still being extended (deletion at
        // p=100 with span 5 keeps the group open through p=104),
        // then have upstream error before the group can close.
        // The error must surface; the in-progress group must NOT
        // be emitted.
        let items: Vec<Item> = vec![
            Ok(del_pp(0, 100, 1, 0, 5)),
            Ok(snp_pp(0, 101, 1, 0)),
            Err(fake_upstream_err(0)),
        ];
        let mut grouper = VariantGrouper::with_config(iter_from(items), GrouperConfig::default());
        match grouper.next() {
            Some(Err(GrouperError::Upstream(MergerError::OutOfOrder { .. }))) => {}
            other => panic!("expected upstream error, got {other:?}"),
        }
        assert!(grouper.next().is_none());
    }

    // ---------- max_var_group_span cap ----------

    #[test]
    fn cap_triggered_by_extension() {
        // Configure span cap = 10.
        // Seed: deletion at p=100 with ref_span=8 (end=107, span=8).
        // Within cap, group opens.
        // Extension: deletion at p=107 with ref_span=5 →
        // attempted_end = max(107, 107+5-1)=111, span=12 > 10.
        // The cap must trip; no group is emitted.
        let config = GrouperConfig {
            max_var_group_span: 10,
        };
        let upstream = ok_iter(vec![del_pp(0, 100, 1, 0, 8), del_pp(0, 107, 1, 0, 5)]);
        let mut grouper = VariantGrouper::with_config(upstream, config);
        match grouper.next() {
            Some(Err(GrouperError::VarGroupTooWide {
                chrom_id,
                start,
                attempted_end,
                attempted_span,
                cap,
            })) => {
                assert_eq!(chrom_id, 0);
                assert_eq!(start, 100);
                assert_eq!(attempted_end, 111);
                assert_eq!(attempted_span, 12);
                assert_eq!(cap, 10);
            }
            other => panic!("expected VarGroupTooWide, got {other:?}"),
        }
        assert!(grouper.next().is_none());
    }

    #[test]
    fn cap_triggered_by_seed_alone() {
        // Configure span cap = 4.
        // Seed: deletion at p=100 with ref_span=5 → seed_end=104,
        // span=5 > 4 → error even before any extension.
        let config = GrouperConfig {
            max_var_group_span: 4,
        };
        let upstream = ok_iter(vec![del_pp(0, 100, 1, 0, 5)]);
        let mut grouper = VariantGrouper::with_config(upstream, config);
        match grouper.next() {
            Some(Err(GrouperError::VarGroupTooWide {
                chrom_id,
                start,
                attempted_end,
                attempted_span,
                cap,
            })) => {
                assert_eq!(chrom_id, 0);
                assert_eq!(start, 100);
                assert_eq!(attempted_end, 104);
                assert_eq!(attempted_span, 5);
                assert_eq!(cap, 4);
            }
            other => panic!("expected VarGroupTooWide, got {other:?}"),
        }
        assert!(grouper.next().is_none());
    }

    #[test]
    fn cap_honoured_at_exact_boundary() {
        // Configure span cap = 5.
        // Seed: deletion at p=100 with ref_span=5 → end=104,
        // span=5. Exactly at the cap, not over — must emit.
        let config = GrouperConfig {
            max_var_group_span: 5,
        };
        let upstream = ok_iter(vec![del_pp(0, 100, 1, 0, 5)]);
        let groups = collect_groups(VariantGrouper::with_config(upstream, config)).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].span(), (0, 100, 104));
    }

    // ---------- pathological / defensive ----------

    #[test]
    fn empty_per_sample_dropped_at_seed_time() {
        // A PerPositionPileups with every slot None has no variant
        // observation and must be dropped without panicking.
        // max_ref_span falls back to 1, but we never reach that
        // path because has_variant_observation rejects first.
        let upstream = ok_iter(vec![empty_pp(0, 100, 2), snp_pp(0, 101, 2, 0)]);
        let groups = collect_groups(VariantGrouper::with_config(
            upstream,
            GrouperConfig::default(),
        ))
        .unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].span(), (0, 101, 101));
    }

    // ---------- config / API smoke tests ----------

    #[test]
    fn default_config_has_documented_cap() {
        let cfg = GrouperConfig::default();
        assert_eq!(cfg.max_var_group_span, DEFAULT_MAX_VAR_GROUP_SPAN);
    }

    #[test]
    fn config_accessor_returns_active_config() {
        let cfg = GrouperConfig {
            max_var_group_span: 42,
        };
        let grouper = VariantGrouper::with_config(iter_from(Vec::<Item>::new()), cfg);
        assert_eq!(grouper.config().max_var_group_span, 42);
    }
}
