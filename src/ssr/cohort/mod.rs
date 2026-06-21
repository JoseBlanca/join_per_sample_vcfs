//! Stage 2 — `ssr-call` (the cohort caller). Reads the N per-sample `.ssr.psp`
//! files, k-way-merges them by catalog locus into the [`CohortLocus`](types::CohortLocus)
//! work-item, and (Phases 2–3) genotypes each locus with the per-locus EM.
//!
//! This Phase-1 spine covers reading & merge only; the parameter pre-pass and the
//! genotyping EM land later. Design docs:
//! `doc/devel/architecture/ssr_call_reading.md` (this layer),
//! `ssr_call_parameters.md`, `ssr_call_genotyping.md`; build plan
//! `doc/devel/implementation_plans/ssr_call_reading.md`.

pub mod reader;
pub mod types;
