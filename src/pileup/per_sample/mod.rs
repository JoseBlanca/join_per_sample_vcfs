//! Stage 1 of the multi-sample calling pipeline: per-sample pileup.
//!
//! Reads one or more coordinate-sorted CRAM files belonging to a
//! single sample, validates their headers against each other and
//! against the reference FASTA, merges them into a single
//! coordinate-sorted record stream, BAQ-prepares the reads, drives
//! the walker, and emits a `.psp` artefact. The per-read filter
//! cascade is specified in `doc/devel/specs/per_sample_pileup.md`
//! §"Read filters".
//!
//! Module composition: BAQ glue (`baq_engine`, `baq_stream`) and the
//! pileup→PSP seam (`pileup_to_psp`) live here as siblings; the BAQ
//! algorithm core lives at [`crate::baq`]; the walker algorithm at
//! [`crate::pileup::walker`]; the per-position record data model at
//! [`crate::pileup_record`]; the on-disk `.psp` format at
//! [`crate::psp`]; the reference-FASTA fetchers and `ContigList` at
//! [`crate::fasta`]; and the mapped-read input slice (CRAM reader,
//! typed CRAM-input errors) at [`crate::bam`].

pub mod baq_engine;
pub mod baq_stream;
pub mod pileup_to_psp;
pub mod read_pipeline;
pub mod read_processor;

#[cfg(test)]
mod baq_tests;
#[cfg(test)]
pub(crate) mod cram_files;
#[cfg(test)]
pub(crate) mod record_specs;
