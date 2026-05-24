//! Mapped-read file input — BAM and CRAM.
//!
//! Named `bam/` (not `cram/`) because BAM is the more common name for
//! mapped-read files and BAM-format support is planned. Today the
//! module only ships the CRAM reader; when BAM lands it will land as
//! a sibling submodule and the shared types (read decoding shape,
//! filter cascade, contig-list cross-checks) will be hoisted into
//! this `mod.rs` or into a `shared` sibling.
//!
//! - [`cram_input`] — CRAM reader. `CramMergedReader` opens a list of
//!   coordinate-sorted CRAMs for one sample, validates their `@SQ`
//!   headers against each other and the reference FASTA, and produces
//!   a single coordinate-sorted [`cram_input::MappedRead`] stream
//!   filtered through [`cram_input::CramMergedReaderConfig`].
//! - [`errors`] — typed I/O errors raised by the CRAM reader. Reused
//!   downstream by the BAQ stream and the cohort CLI's error bridge,
//!   which is why they live here at the module root rather than
//!   inside `cram_input`.
//! - [`index_preflight`] — `.crai` / `.bai` / `.csi` detection +
//!   opt-in auto-build for callers that need per-contig random
//!   access (e.g. the per-chromosome parallel `var-calling-from-bam`
//!   driver).
//!
//! The walker-side input types ([`MappedRead`](cram_input::MappedRead)
//! → [`crate::pileup::walker::PreparedRead`]) are produced here and
//! consumed by the BAQ glue ([`crate::pileup::per_sample::baq_engine`])
//! and then by the walker. CIGAR ops use [`crate::pileup::walker::CigarOp`].

pub mod cram_input;
pub mod errors;
pub mod index_preflight;
