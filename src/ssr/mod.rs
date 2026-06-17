// Net-new Mark-2 SSR module, built incrementally; no production consumers wire
// into it yet (the CLI still drives `ssr_mark1` until the Stage-1 cutover), so
// allow dead code while the tree fills in.
#![allow(dead_code)]
//! SSR/STR genotyping — Mark 2 (the empirical-candidate model).
//!
//! Genotypes short tandem repeats (microsatellites, period ≤ 6) from aligned
//! reads. Candidate alleles are **observed sequences** assembled from the reads,
//! not rungs generated off the reference (the Mark-1 model); the reference is only
//! a coordinate frame. The three stages mirror the SNP path's
//! `pileup → .psp → var-calling`:
//!
//! ```text
//! reference ─► ssr-catalog ─► catalog ─► ssr-pileup ─► .ssr.psp ─► ssr-call ─► VCF
//! ```
//!
//! Design docs: `doc/devel/architecture/ssr_ladder_model.md` (the Mark-2 model),
//! `doc/devel/architecture/ssr_pileup_mark2.md` (Stage 1), and the build plan
//! `doc/devel/implementation_plans/ssr_pileup_mark2.md`.
//!
//! Built so far: [`types`] (`Locus`, `Motif`), [`catalog`] (Stage 0, copied
//! verbatim from the Mark-1 tree — Stage 0 is model-agnostic), and the Stage-1
//! read path in [`pileup`] (footprint geometry + reservoir fetch).

pub mod catalog;
pub mod pileup;
pub mod types;
