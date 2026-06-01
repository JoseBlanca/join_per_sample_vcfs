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
//! - [`loader`] — the batch [`load_chunk_from_iters`] (retained as the
//!   streaming loader's equivalence-test oracle) and the streaming
//!   [`StreamingBlockLoader`], which folds the cohort forward through
//!   [`SpanColumnSource`]s and compacts closed variant groups
//!   incrementally — the production producer's engine.

pub mod column_span_reader;
pub mod columns;
pub mod driver;
pub mod kernels;
pub mod loader;
pub mod partition;
#[cfg(test)]
pub(crate) mod test_helpers;
pub mod worker;

pub use column_span_reader::ColumnSpanReader;
pub use columns::{MaterialisedChunk, SampleColumns};
pub use driver::{
    ChunkDriverError, ChunkDriverParams, ChunkDriverStats, ChunkSizingParams,
    DEFAULT_CHUNK_GENOMIC_SPAN, DownstreamFilterParams, drive_cohort_chunked,
};
pub use loader::{
    ChunkLoadError, ChunkLoadExtent, ChunkLoadScratch, ChunkLoadStats, SpanColumnSource,
    StreamingBlockLoader, load_chunk_from_iters,
};
pub use partition::{PartitionError, PartitionScratch, WindowPartition, partition_window};
pub use worker::{WindowRunStats, WorkerScratch, run_window};
