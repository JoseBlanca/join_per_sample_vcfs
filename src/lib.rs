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

pub mod bam;
pub mod baq;
pub mod fasta;
pub mod genetics;
pub mod iter_ext;
pub mod ng;
pub(crate) mod norm_seqs;
pub mod paralog;
pub mod pileup;
pub mod pileup_record;
pub mod pop_var_caller;
pub mod pop_var_caller_exp;
pub mod psp;
pub mod regions;
pub mod sample_summary;
pub mod ssr;
pub mod var_calling;
pub mod vcf;
