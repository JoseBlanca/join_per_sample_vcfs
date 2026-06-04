//! Re-architected cohort `.psp` → VCF pipeline (parallel build package).
//!
//! This is the **`var_calling_new::`** package from the re-architecture plan —
//! built alongside the shipping [`var_calling`](crate::var_calling) until the
//! one-commit swap (execution-plan §P7). The design docs:
//!
//! - `doc/devel/implementation_plans/re_architecture_streaming_pipeline.md` —
//!   the architecture, constraints, build strategy;
//! - `doc/devel/implementation_plans/re_architecture_module_outline.md` —
//!   the module/type map (the names below are its **target** names);
//! - `doc/devel/implementation_plans/re_architecture_execution_plan.md` —
//!   the ordered, byte-identity-gated phases.
//!
//! ## Structure (rebuilt) vs kernels (copied verbatim)
//!
//! The pipeline structure is rebuilt from scratch; the byte-identity-sensitive
//! numeric kernels are **copied verbatim** from [`var_calling`](crate::var_calling)
//! (plan §7 / principle 2):
//!
//! - **Rebuilt** (Phases 1–4): [`types`], [`sample_reader`],
//!   [`cohort_integration`], [`pileup_overlaps`], [`em_posterior_calc`]
//!   (the `VariantCaller` worker), [`vcf_writer`], [`pipeline`].
//! - **Copied verbatim** (the kernels — a closed, columnar-free subset whose
//!   only outside dependencies are the shared `psp` / `fasta` / `pileup_record`
//!   crate modules): [`per_group_merger`], [`posterior_engine`],
//!   [`variant_grouping`], [`per_position_merger`], [`dust_filter`],
//!   [`contamination_estimation`].
//!
//! The columnar `kernels/` chain (`unify_alleles_columnar` /
//! `project_scalars_columnar` / `compute_log_likelihoods_columnar`) and its
//! `columns.rs` / `partition.rs` / `MaterialisedChunk` dependencies are
//! deliberately **not** transplanted: the row [`per_group_merger`] already
//! produces the flat SoA `log_likelihoods` buffer the SIMD EM consumes, and
//! the closed-form log-likelihood is scalar in both paths — so the record-based
//! path keeps the SIMD EM with no SIMD loss (appendix §D.2 conditional,
//! verified against the code in Phase 0).

// Build-phase only: the package is wired incrementally (P1→P4), so verbatim
// kernels carry pub(crate) items whose only callers were the dropped columnar
// driver/worker. Removed at the P7 swap, when the package becomes
// `var_calling::` and the old consumers are deleted/ported.
#![allow(dead_code)]

// --- Copied verbatim from `var_calling` (the numeric kernels) ---
pub mod contamination_estimation;
pub mod dust_filter;
pub mod per_group_merger;
pub mod per_position_merger;
pub mod posterior_engine;
pub mod variant_grouping;

// --- Rebuilt structure (Phases 1–4) ---
pub mod cohort_integration;
pub mod em_posterior_calc;
pub mod pileup_overlaps;
pub mod pipeline;
pub mod sample_reader;
pub mod types;
pub mod vcf_writer;
