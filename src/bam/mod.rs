//! Mapped-read file input — BAM and CRAM.
//!
//! Named `bam/` (not `cram/`) because BAM is the more common name for
//! mapped-read files. The module hosts the format-agnostic header-
//! validation + per-read filter surface (in [`alignment_input`]) plus
//! two per-format open helpers ([`cram_input`] and [`bam_input`]); the
//! actual per-segment read fetch lives in [`segment_reader`] /
//! [`segment_merge`].
//!
//! - [`alignment_input`] — validates a list of coordinate-sorted
//!   alignment files for one sample (their `@SQ` headers against each
//!   other, a single sample name, coordinate sort order) via
//!   [`alignment_input::load_pileup_inputs`], and owns the
//!   format-agnostic per-read filter cascade
//!   ([`alignment_input::AlignmentMergedReaderConfig`]) plus the shared
//!   `RecordBuf` → [`alignment_input::MappedRead`] conversion the
//!   segment readers feed through.
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
// The multi-file coordinate merge over `segment_reader` — one sample's
// per-file segment streams k-way-merged into the coordinate-sorted stream
// the pileup walker consumes. Increment #5 (SNP `--regions` retrofit); no
// production consumer until `run_pileup` is flipped onto it.
pub(crate) mod segment_merge;
