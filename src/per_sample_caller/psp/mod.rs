//! `.psp` per-sample pileup byte format — writer and reader.
//!
//! Stage 1 of the pipeline emits per-sample pileup records (one per
//! covered reference position); this module is the bridge between
//! those in-memory records and the on-disk artefact downstream stages
//! consume. The byte layout is specified in
//! `ia/specs/per_sample_pileup_format.md`; this module is its
//! authoritative implementation.
//!
//! Module layout:
//!
//! - [`registry`]: the v1.0 column-tag registry — the single source
//!   of truth shared by writer, reader, and `psp_spec_dump`.
//! - [`varint`]: LEB128 / zig-zag-LEB128 codecs.
//! - [`errors`]: typed error enums.
//! - [`header`], [`block`], [`index`], [`trailer`]: per-section
//!   wire codecs.
//! - [`writer`]: streaming [`writer::PspWriter`] over any `Write` sink.
//!
//! The reader-side public API (`PspReader`) is planned but not yet
//! implemented. The decoder primitives in `block`, `index`,
//! `header`, and `trailer` are already in place and exercised by
//! their respective `#[cfg(test)] mod tests` round-trip suites.

pub mod block;
pub mod errors;
pub mod header;
pub mod index;
pub mod registry;
pub mod trailer;
pub mod varint;
pub mod writer;

pub use errors::{PspReadError, PspWriteError};
