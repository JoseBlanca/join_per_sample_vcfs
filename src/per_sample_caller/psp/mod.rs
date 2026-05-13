//! `.psp` per-sample pileup byte format — writer and reader.
//!
//! Stage 1 of the pipeline emits per-sample pileup records (one per
//! covered reference position); this module is the bridge between
//! those in-memory records and the on-disk artefact downstream stages
//! consume. The byte layout is specified in
//! `ia/specs/per_sample_pileup_format.md`; this module is its
//! authoritative implementation.
//!
//! The first slice of this module covers the shared primitives:
//! - [`registry`]: the v1.0 column-tag registry — the single source
//!   of truth shared by writer, reader, and (later) `psp_spec_dump`.
//! - [`varint`]: LEB128 / zig-zag-LEB128 encoders and decoders used
//!   everywhere in the binary body.
//! - [`errors`]: typed error enums.
//!
//! Block layout, header build/parse, the index, the trailer, and the
//! `PspWriter` / `PspReader` public API land in subsequent slices.

pub mod block;
pub mod errors;
pub mod header;
pub mod index;
pub mod registry;
pub mod trailer;
pub mod varint;

pub use errors::{PspReadError, PspWriteError};
