//! Mapped-read file input — BAM and CRAM.
//!
//! Named `bam/` (not `cram/`) because BAM is the more common name for
//! mapped-read files. Today this module hosts the format-agnostic
//! merge / filter / header-validation surface (in
//! [`alignment_input`]) plus a CRAM record-stream decoder; BAM
//! record-stream support lands as a sibling submodule
//! (`bam_input`) that produces the same `OpenAlignmentFile` shape.
//!
//! - [`alignment_input`] — `AlignmentMergedReader` opens a list of
//!   coordinate-sorted alignment files for one sample, validates
//!   their `@SQ` headers against each other and the reference
//!   FASTA, and produces a single coordinate-sorted
//!   [`alignment_input::MappedRead`] stream filtered through
//!   [`alignment_input::AlignmentMergedReaderConfig`]. The merge
//!   loop and filter cascade are format-agnostic; only the
//!   per-input record-stream factories (currently CRAM-only,
//!   inlined here) need to grow when BAM support lands.
//! - [`errors`] — typed I/O errors raised by the alignment-file
//!   readers. Reused downstream by the BAQ stream and the cohort
//!   CLI's error bridge, which is why they live here at the
//!   module root rather than inside `alignment_input`.
//! - [`index_preflight`] — `.crai` / `.bai` / `.csi` detection +
//!   opt-in auto-build for callers that need per-contig random
//!   access (e.g. the per-chromosome parallel `var-calling-from-bam`
//!   driver).
//!
//! The walker-side input types ([`MappedRead`](alignment_input::MappedRead)
//! → [`crate::pileup::walker::PreparedRead`]) are produced here and
//! consumed by the BAQ glue ([`crate::pileup::per_sample::baq_engine`])
//! and then by the walker. CIGAR ops use [`crate::pileup::walker::CigarOp`].

pub mod alignment_input;
pub mod bam_input;
pub mod cram_input;
pub mod errors;
pub mod index_preflight;
