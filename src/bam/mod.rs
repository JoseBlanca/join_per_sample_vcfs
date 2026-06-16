//! Mapped-read file input — BAM and CRAM.
//!
//! Named `bam/` (not `cram/`) because BAM is the more common name for
//! mapped-read files. The module hosts the format-agnostic merge /
//! filter / header-validation surface (in [`alignment_input`]) plus
//! two per-format record-stream decoders ([`cram_input`] and
//! [`bam_input`]), both producing the same `OpenAlignmentFile`
//! shape the merge consumes.
//!
//! - [`alignment_input`] — `AlignmentMergedReader` opens a list of
//!   coordinate-sorted alignment files for one sample, validates
//!   their `@SQ` headers against each other and the reference
//!   FASTA, and produces a single coordinate-sorted
//!   [`alignment_input::MappedRead`] stream filtered through
//!   [`alignment_input::AlignmentMergedReaderConfig`]. The merge
//!   loop and filter cascade are format-agnostic; per-input
//!   record-stream construction dispatches on file extension
//!   into the format-specific sibling modules.
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
// `bam_input` and `cram_input` are the per-format record-stream
// decoders. They expose only `pub(super)` items consumed by
// `alignment_input`; no caller outside `crate::bam` references
// them. `pub(crate)` keeps the path scope honest.
pub(crate) mod bam_input;
pub(crate) mod cram_input;
pub mod errors;
pub mod index_preflight;
// The shared indexed-segment read source: pooled, thread-safe
// per-segment read queries over one BAM/CRAM. Reuses the per-format
// scanners' decode helpers; adds reader ownership + a pool. Consumed
// by the SSR Stage-1 fetcher (and, later, the SNP `--regions` path).
pub(crate) mod segment_reader;
