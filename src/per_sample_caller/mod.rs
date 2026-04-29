//! Stage 1 of the multi-sample calling pipeline: per-sample caller.
//!
//! This slice of the per-sample caller reads one or more
//! coordinate-sorted CRAM files belonging to a single sample, validates
//! their headers against each other and against the reference FASTA,
//! merges them into a single coordinate-sorted record stream, and
//! applies the per-read filter cascade specified in
//! `ia/specs/per_sample_caller.md` §"Read filters".
//!
//! Downstream stages (BAQ, pileup walking, allele extraction, phase
//! chains, `.psf` writing) are separate slices and live in sibling
//! modules.

pub mod cram_input;
pub mod errors;

#[cfg(test)]
pub(crate) mod cram_files;
#[cfg(test)]
pub(crate) mod record_specs;
