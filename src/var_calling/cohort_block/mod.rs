//! Cohort-level chunk-at-a-time materialisation.
//!
//! At each step of the chunk loop, the driver loads one genomic
//! chunk × N samples into memory (variant positions only),
//! partitions it into autonomous windows whose boundaries fall in
//! `max_group_span`-wide safe gaps, and runs the worker pipeline
//! on each window to produce final variant records. See
//! `doc/devel/implementation_plans/cohort_within_chromosome_parallel.md`
//! for the architectural rationale and phasing.
//!
//! **Submodule layout.**
//! - [`columns`] — [`SampleColumns`] (per-sample CSR columnar
//!   storage) and [`MaterialisedChunk`] (one chunk × N samples).
//! - [`loader`] — [`load_chunk_from_iters`] +
//!   [`ChunkLoadScratch`] + [`ChunkLoadError`]: read per-sample
//!   record iterators, apply the cohort-wide variant-position
//!   filter, compact survivors into the chunk.
//! - [`chunk_boundaries`] — [`finalise_chunk_boundaries`] +
//!   [`BoundaryFinalisationScratch`] + [`BoundaryFinalisationError`]: pick the
//!   chunk's `safe_end`, partition `[range.start, safe_end)` into
//!   windows, split records past `safe_end` into carryover.

pub mod chunk_boundaries;
pub mod columns;
pub mod driver;
pub mod kernels;
pub mod loader;
pub mod partition;
#[cfg(test)]
pub(crate) mod test_helpers;
pub mod worker;

pub use chunk_boundaries::{
    BoundaryFinalisationError, BoundaryFinalisationScratch, finalise_chunk_boundaries,
};
pub use columns::{MaterialisedChunk, SampleColumns};
pub use driver::{
    ChunkDriverError, ChunkDriverParams, ChunkDriverStats, ChunkSizingParams,
    DEFAULT_CHUNK_GENOMIC_SPAN, DownstreamFilterParams, drive_cohort_chunked,
};
pub use loader::{
    ChunkLoadError, ChunkLoadExtent, ChunkLoadScratch, ChunkLoadStats, load_chunk_from_iters,
};
pub use partition::{PartitionError, PartitionScratch, WindowPartition, partition_window};
pub use worker::{WindowRunStats, WorkerScratch, into_shared_ref_fetcher, run_window};
