//! Cohort `.psp` → multi-sample VCF pipeline (the production record-streaming
//! topology).
//!
//! The runtime is a three-section parallel pipeline — producer
//! ([`cohort_integration`]) → caller ([`variant_caller`]) → writer
//! ([`vcf_writer`]) — wired by [`pipeline`] with two bounded crossbeam
//! channels inside a `std::thread::scope`. The output is byte-identical for
//! any worker count.
//!
//! The architecture and its design rationale are documented in the
//! re-architecture plans under `doc/devel/implementation_plans/`
//! (`re_architecture_streaming_pipeline.md`, `re_architecture_module_outline.md`,
//! `re_architecture_execution_plan.md`).
//!
//! ## Structure modules vs numeric kernels
//!
//! Two groups of modules:
//!
//! - **Structure** — the producer/caller/writer plumbing and the interchange
//!   types: [`types`], [`sample_reader`], [`cohort_integration`],
//!   [`pileup_overlaps`], [`variant_caller`] (the `VariantCaller` worker),
//!   [`vcf_writer`], [`pipeline`].
//! - **Numeric kernels** — the byte-identity-sensitive math (a closed,
//!   columnar-free subset whose only outside dependencies are the shared
//!   `psp` / `fasta` / `pileup_record` crate modules): [`per_group_merger`],
//!   [`posterior_engine`], [`variant_grouping`], [`per_position_merger`],
//!   [`dust_filter`], [`contamination_estimation`]. These carry the
//!   variant-calling arithmetic byte-for-byte from the pre-rewrite pipeline;
//!   byte-identity of the emitted VCF is verified out-of-tree (against the
//!   pre-rewrite pipeline + the GIAB benchmark), so edits to them must hold
//!   that contract.
//!
//! The row [`per_group_merger`] produces the flat SoA `log_likelihoods` buffer
//! the SIMD EM consumes directly, and the closed-form log-likelihood is scalar
//! — so the record-based path runs the SIMD EM with no SIMD loss (no separate
//! columnar likelihood path is needed).

// --- Numeric kernels (byte-identity-sensitive math) ---
pub mod contamination_estimation;
pub mod dust_filter;
pub mod per_group_merger;
pub mod per_position_merger;
pub mod posterior_engine;
pub mod variant_grouping;

// --- Structure (producer / caller / writer plumbing) ---
pub mod cohort_integration;
pub mod pileup_overlaps;
pub mod pipeline;
pub mod sample_reader;
pub mod types;
pub mod variant_caller;
pub mod vcf_writer;

/// Default `--min-qual` (phred) emission gate (CLI default).
pub const DEFAULT_MIN_QUAL_PHRED: f64 = 30.0;
/// Default `--min-alt-obs-per-sample` (CLI default).
pub const DEFAULT_MIN_ALT_OBS_PER_SAMPLE: u32 = 2;
/// Default `--min-mapq-diff-t` Welch's-t threshold (CLI default).
pub const DEFAULT_MIN_MAPQ_DIFF_T: f32 = -3.0;

#[cfg(test)]
mod test_helpers;
