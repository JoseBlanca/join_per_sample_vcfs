//! Caller step 2 — per-group merge + per-record EM (appendix §D.2), plus the
//! `VariantCaller` worker that composes grouping (§D.1) with the math here.
//!
//! *(today: `var_calling::per_group_merger`, `posterior_engine`, `worker`)*
//!
//! `OverlappingPileupRecords` → `Vec<Variant>`, one group at a time:
//!
//! - per-group merge via the copied
//!   [`per_group_merger`](crate::var_calling_new::per_group_merger) — its
//!   row-shape `MergedRecord` already carries the **flat SoA**
//!   `log_likelihoods` buffer (`[sample * n_genotypes + g]`), so it feeds the
//!   SIMD EM directly through `MergedAllelesView` — **no SIMD loss**;
//! - per-record EM via the copied
//!   [`posterior_engine`](crate::var_calling_new::posterior_engine) — the
//!   `InterpUnivariateSimdMath` lane-of-4 `ln`/`exp` backend, kept verbatim.
//!
//! > **Design note (verified Phase 0).** The columnar `kernels/` chain
//! > (`unify_alleles_columnar` / `project_scalars_columnar` /
//! > `compute_log_likelihoods_columnar`) and its `MaterialisedChunk` /
//! > `WindowPartition` dependencies are **not** transplanted: they were a
//! > byte-identical reimplementation of the row `per_group_merger`, whose
//! > closed-form log-likelihood is scalar in *both* paths (the SIMD is only in
//! > the EM iteration, which fills its lanes from the flat SoA buffer either
//! > path produces). This is the design's §D.2 conditional — "row-shape
//! > `MergedRecord` is fine as long as it hands SoA likelihood buffers to the
//! > EM" — with its predicate confirmed against the code.
//!
//! Phase 3 builds the `VariantCaller` worker here.

// Phase 3.
