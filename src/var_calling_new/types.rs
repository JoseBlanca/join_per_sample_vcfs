//! Shared data types for the re-architected cohort `.psp` → VCF pipeline
//! (appendix §C of the re-architecture design).
//!
//! These are the structures that flow between the three runtime sections
//! (producer → caller → writer):
//!
//! - [`CohortPileupRecord`] — one variable cohort position with every
//!   sample's pileup data at it (the producer's emit unit, the grouper's
//!   input). Revived from `per_position_merger::PerPositionPileups`.
//! - [`PileupCohortChunk`] — the producer→caller work-unit: a safe-gap-bounded
//!   list of [`CohortPileupRecord`]s + its [`RefSpan`], tagged `chunk_order`.
//! - [`RefSpan`] — the chunk's contiguous REF bytes, fetched once on the
//!   producer (monotonic-forward) and sliced per group by the caller.
//! - [`Variant`] — the final emitted record (today's `PosteriorRecord`).
//! - [`CalledChunk`] — the caller→writer payload (today's `BlockResult`).
//!
//! Phase 0: the types land in Phase 1; this stub fixes the module home and
//! the doc map. See
//! `doc/devel/implementation_plans/re_architecture_execution_plan.md` §P1.

// Phase 1 fills this module. The concrete field shapes are pinned in
// appendix §C and reuse the existing per-position record form (cf.
// `PerPositionPileups`), so they are deferred to P1 alongside
// `SamplePspChunk`'s typed getters.
