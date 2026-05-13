//! Typed error enums for the `.psp` writer and reader.
//!
//! Two top-level enums — `PspWriteError` for the writer side,
//! `PspReadError` for the reader — match the producer / consumer
//! split in the spec's §"Header-binary consistency: required reader
//! checks" and §"Block invariants". Sub-enums (e.g. `VarintError`)
//! describe failures intrinsic to a primitive and get wrapped into
//! the top-level enums with context at the call site.
//!
//! Each variant is one concrete failure mode. No `Other(String)`
//! catch-alls — callers that need to react to a specific failure
//! match on the variant; tests assert on the exact variant rather
//! than stringly-comparing messages.
//!
//! New variants land as later modules need them; this file currently
//! covers what the registry / varint / wire-primitive slice requires.

use thiserror::Error;

/// Failures intrinsic to LEB128 / zig-zag-LEB128 decoding. Always
/// wrapped by [`PspReadError::Varint`] (the read side is the only
/// side that decodes — the writer constructs varints by construction
/// and cannot fail mid-encode).
#[derive(Error, Debug, Clone, Copy, PartialEq, Eq)]
pub enum VarintError {
    /// Buffer ran out mid-varint: a continuation bit was set on the
    /// last byte available, but no more bytes were present.
    #[error("varint truncated: continuation bit set on the final available byte")]
    Truncated,

    /// More than 10 LEB128 continuation bytes were consumed. The
    /// spec caps the varint width at 10 bytes (covers `u64`); a
    /// longer encoding cannot represent a valid `u64` and indicates
    /// either corruption or a writer bug.
    #[error("varint overflow: continuation bytes exceeded the 10-byte cap")]
    Overflow,
}

/// Errors the `.psp` reader can emit. Variants land as the
/// corresponding slice ships; the registry / varint slice
/// contributes the [`Self::Varint`] case.
#[derive(Error, Debug)]
pub enum PspReadError {
    #[error("varint decode failed{}: {source}", match context {
        Some(c) => format!(" while reading {c}"),
        None => String::new(),
    })]
    Varint {
        /// What was being decoded when the failure occurred — e.g.
        /// `"block header n_records"`. Optional so callers can pass
        /// the bare primitive error up where context is added higher.
        context: Option<&'static str>,
        #[source]
        source: VarintError,
    },

    #[error("I/O error reading {context}: {source}")]
    Io {
        context: &'static str,
        #[source]
        source: std::io::Error,
    },
}

/// Errors the `.psp` writer can emit. Variants land as the
/// corresponding slice ships.
#[derive(Error, Debug)]
pub enum PspWriteError {
    #[error("I/O error writing {context}: {source}")]
    Io {
        context: &'static str,
        #[source]
        source: std::io::Error,
    },
}
