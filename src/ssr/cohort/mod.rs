//! Stage 2 — `ssr-call` (the cohort caller). Reads the N per-sample `.ssr.psp`
//! files, k-way-merges them by catalog locus into the [`CohortLocus`](types::CohortLocus)
//! work-item, and (Phases 2–3) genotypes each locus with the per-locus EM.
//!
//! This Phase-1 spine covers reading & merge only; the parameter pre-pass and the
//! genotyping EM land later. Design docs:
//! `doc/devel/architecture/ssr_call_reading.md` (this layer),
//! `ssr_call_parameters.md`, `ssr_call_genotyping.md`; build plan
//! `doc/devel/implementation_plans/ssr_call_reading.md`.

pub mod candidate_set;
pub mod driver;
pub mod merge;
pub mod pair_hmm;
pub mod param_estimation;
pub mod reader;
pub mod rung_ladder;
pub mod stutter;
pub mod types;

#[cfg(test)]
pub(crate) mod sim;
#[cfg(test)]
pub(crate) mod test_support;
