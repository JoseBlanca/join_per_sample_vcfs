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
