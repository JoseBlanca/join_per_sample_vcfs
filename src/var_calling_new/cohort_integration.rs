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
//! Phase 2 builds this module.

// Phase 2.
