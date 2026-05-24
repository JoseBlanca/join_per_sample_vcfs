//! Stage 1 of the multi-sample calling pipeline: per-sample pileup.
//!
//! This slice of the per-sample pileup reads one or more
//! coordinate-sorted CRAM files belonging to a single sample, validates
//! their headers against each other and against the reference FASTA,
//! merges them into a single coordinate-sorted record stream, and
//! applies the per-read filter cascade specified in
//! `doc/devel/specs/per_sample_pileup.md` §"Read filters".
//!
//! Downstream stages (BAQ, pileup walking, allele extraction, phase
//! chains, `.psp` writing) are separate slices and live in sibling
//! modules.

pub mod baq_engine;
pub mod baq_stream;
pub mod cram_input;
pub mod errors;
pub mod pileup;
pub mod pileup_to_psp;
pub mod psp;
pub mod ref_fetcher;

#[cfg(test)]
mod baq_tests;
#[cfg(test)]
pub(crate) mod cram_files;
#[cfg(test)]
pub(crate) mod record_specs;
