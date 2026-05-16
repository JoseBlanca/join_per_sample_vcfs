//! Cohort-side stages of the multi-sample SNP caller.
//!
//! Sits between Stage 2 (the per-sample `.psp` reader) and Stages 3-6
//! (DUST filter, grouping, per-group, posterior). See
//! `doc/devel/specs/calling_pipeline_architecture.md` for the full
//! stage breakdown.
//!
//! First occupant: [`per_position_merger`], the linear-scan k-way
//! merge over per-sample `.psp` record streams. Stage 4 lands in
//! [`variant_grouping`], the streaming overlap bundler that turns
//! the merger's output into `OverlappingVarGroup`s for Stage 5.
//! Later stages land in sibling modules as they are implemented.

pub mod per_group_merger;
pub mod per_position_merger;
pub mod variant_grouping;
