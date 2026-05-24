//! Pileup: the per-position read-folding model.
//!
//! Two submodules:
//!
//! - [`walker`] — the sequential walker algorithm. Takes a
//!   coordinate-sorted stream of [`walker::PreparedRead`]s and emits
//!   [`crate::pileup_record::PileupRecord`]s. Algorithmic core; does
//!   not know about CRAM, BAQ, or the on-disk PSP format.
//! - [`per_sample`] — Stage 1 of the multi-sample calling pipeline.
//!   Composes the CRAM input slice, the BAQ glue (`baq_engine`,
//!   `baq_stream`), the walker, and the [`crate::psp`] writer into a
//!   per-sample BAM/CRAM → `.psp` pipeline.
//!
//! The split mirrors the BAQ / FASTA / PSP shape: the walker is a
//! reusable algorithm sitting next to the stage that drives it.

pub mod per_sample;
pub mod walker;
