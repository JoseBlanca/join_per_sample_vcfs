//! Variant-calling stages (Stages 3–6).
//!
//! Sits between Stage 2 (the per-sample `.psp` reader) and the
//! posterior engine: DUST filter, grouping, per-group merger,
//! posterior. Used in cohort mode across many samples, but the same
//! code path applies to a single sample. See
//! `doc/devel/specs/calling_pipeline_architecture.md` for the full
//! stage breakdown.
//!
//! First occupant: [`per_position_merger`], the linear-scan k-way
//! merge over per-sample `.psp` record streams. Stage 4 lands in
//! [`variant_grouping`], the streaming overlap bundler that turns
//! the merger's output into `OverlappingVariantGroup`s for Stage 5.
//! Later stages land in sibling modules as they are implemented.

pub mod cohort_block;
pub mod contamination_estimation;
pub mod dust_filter;
pub mod per_group_merger;
pub mod per_position_merger;
pub mod posterior_engine;
pub mod variant_grouping;
