//! Caller step 1 — grouping (appendix §D.1).
//!
//! Turns a chunk's `CohortPileupRecord`s into `OverlappingPileupRecords`:
//! walk the records, joining ones whose reach overlaps into groups (bridging
//! the gaps left by dropped dust positions), via the
//! [`variant_grouping`](crate::var_calling::variant_grouping) grouper. The
//! grouping is byte-identical by construction (and is not perf-critical).
//!
//! Built and consumed inside the caller — `OverlappingPileupRecords` never
//! crosses a queue.

use crate::var_calling::per_position_merger::{PerPositionMergerError, PerPositionPileups};
use crate::var_calling::types::CohortPileupRecord;
use crate::var_calling::variant_grouping::{
    GrouperConfig, GrouperError, OverlappingVariantGroup, VariantGrouper,
};

/// Group a chunk's variable [`CohortPileupRecord`]s into overlapping variant
/// groups — the record-based grouping (appendix §D.1), promoting the copied
/// [`VariantGrouper`] to production.
///
/// Consumes the chunk's records (already only the variable positions; dust
/// positions were dropped in the producer) and yields one
/// [`OverlappingVariantGroup`] per cluster of reach-overlapping positions. The
/// grouper bridges the gaps left by dropped positions exactly as `main`'s
/// `DustFilter → VariantGrouper` chain does — pure-REF positions inside an open
/// group fold in, pure-REF positions outside any group are dropped — so this is
/// byte-identical by construction.
///
/// Returned lazily (the grouper is a streaming iterator); the
/// [`VariantCaller`](crate::var_calling::variant_caller::VariantCaller)
/// consumes it group-by-group, so the grouped records never cross a queue.
pub fn overlapping_groups(
    records: Vec<CohortPileupRecord>,
    config: GrouperConfig,
) -> impl Iterator<Item = Result<OverlappingVariantGroup, GrouperError>> {
    // The grouper is generic over an upstream yielding `Result<PerPositionPileups, E>`
    // (its production upstreams are the per-position merger / dust filter). Our
    // records are already the per-position cohort pileups — `CohortPileupRecord`
    // is the same shape as `PerPositionPileups` — so we re-wrap them with an
    // infallible `Ok`.
    let upstream = records.into_iter().map(|r| {
        Ok::<PerPositionPileups, PerPositionMergerError>(PerPositionPileups {
            chrom_id: r.chrom_id,
            pos: r.pos,
            per_sample: r.per_sample,
        })
    });
    VariantGrouper::with_config(upstream, config)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup_record::PileupRecord;
    use crate::var_calling::test_helpers::allele;

    /// One single-sample variable position: REF `ref_seq` + an ALT `T` so the
    /// position seeds a group. `ref_seq.len()` is the REF span (the grouping
    /// reach driver).
    fn rec(pos: u32, ref_seq: &[u8]) -> CohortPileupRecord {
        let alleles = vec![allele(ref_seq, 5, -1.0, &[]), allele(b"T", 3, -1.5, &[])];
        CohortPileupRecord {
            chrom_id: 0,
            pos,
            per_sample: vec![Some(PileupRecord::new(0, pos, alleles))],
        }
    }

    #[test]
    fn overlapping_groups_bridges_by_reach_and_splits_on_gap() {
        // pos 10 has a 3-base REF span ⇒ reach [10, 12]; pos 12 falls inside it
        // and joins the same group. pos 1000 is far past the reach ⇒ its own
        // group. Verifies the reach-overlap join + the gap split.
        let records = vec![rec(10, b"AAA"), rec(12, b"A"), rec(1000, b"A")];
        let groups: Vec<_> = overlapping_groups(records, GrouperConfig::default())
            .map(|g| g.expect("grouping"))
            .collect();
        assert_eq!(groups.len(), 2, "two groups: {groups:?}");
        assert_eq!((groups[0].start, groups[0].end), (10, 12));
        assert_eq!(groups[0].records.len(), 2, "pos 10 + 12 bridged");
        assert_eq!((groups[1].start, groups[1].end), (1000, 1000));
        assert_eq!(groups[1].records.len(), 1);
    }

    #[test]
    fn overlapping_groups_empty_input_yields_no_groups() {
        let groups: Vec<_> = overlapping_groups(Vec::new(), GrouperConfig::default()).collect();
        assert!(groups.is_empty());
    }
}
