//! Caller step 1 — grouping (appendix §D.1).
//!
//! *(today: the `#[cfg(test)]` `build_overlapping_variant_group` oracle +
//! [`variant_grouping`](crate::var_calling_new::variant_grouping))*
//!
//! Turns a chunk's `CohortPileupRecord`s into `OverlappingPileupRecords`:
//! walk the records, joining ones whose reach overlaps into groups (bridging
//! the gaps left by dropped dust positions). **Record-based** — the row
//! grouper (`OverlappingVariantGroup`) is promoted to production; this is
//! byte-identical by construction (and grouping is not perf-critical).
//!
//! Built and consumed inside the caller — `OverlappingPileupRecords` never
//! crosses a queue.
//!
//! Phase 3 builds this module on top of the copied
//! [`variant_grouping`](crate::var_calling_new::variant_grouping).

// Phase 3.
