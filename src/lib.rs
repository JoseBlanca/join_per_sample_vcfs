//! # `pop_var_caller`
//!
//! Multi-sample variant caller — per-sample → cohort merge pipeline.
//! See `ia/specs/calling_pipeline_architecture.md` for the stage
//! breakdown.
//!
//! ## Feature flags
//!
//! - `dhat-heap` — opt-in `dhat::Alloc` global allocator for heap
//!   profiling under benches and examples. Bench/example use only;
//!   not for production builds.
//! - `alloc-mimalloc` — opt-in `mimalloc` global allocator for
//!   benches and examples. Cannot be combined with `dhat-heap` at
//!   the `#[global_allocator]` declaration in a bench/example shim.

#![forbid(unsafe_code)]

pub mod baq;
pub mod fasta;
pub mod iter_ext;
pub mod per_sample_pileup;
pub mod pileup_record;
pub mod pop_var_caller;
pub mod psp;
pub mod var_calling;
