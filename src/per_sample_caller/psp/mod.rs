//! `.psp` per-sample pileup byte format — writer and reader.
//!
//! Stage 1 of the pipeline emits per-sample pileup records (one per
//! covered reference position); this module is the bridge between
//! those in-memory records and the on-disk artefact downstream stages
//! consume. The byte layout is specified in
//! `ia/specs/per_sample_pileup_format.md`; this module is its
//! authoritative implementation.
//!
//! Module layout (most sub-modules are `pub(crate)` — accessed
//! through the [`writer::PspWriter`] entry point and the public
//! [`header`] type surface):
//!
//! - `registry`: the v1.0 column-tag registry — the single source
//!   of truth shared by writer, reader, and `psp_spec_dump`.
//! - `varint`: LEB128 / zig-zag-LEB128 codecs.
//! - `errors`: typed error enums (re-exported below).
//! - [`header`], `block`, `index`, `trailer`: per-section wire codecs.
//! - [`writer`]: streaming [`writer::PspWriter`] over any `Write` sink.
//!
//! The reader-side public API (`PspReader`) is planned but not yet
//! implemented. The decoder primitives are already in place and
//! exercised by their respective `#[cfg(test)] mod tests` round-trip
//! suites.

// Mi5: most submodules are crate-internal scaffolding for the
// writer/reader implementations. The bench imports `header::*` and
// `writer::PspWriter`; everything else is reachable only through
// the re-exports below or through the writer's API surface.
pub(crate) mod block;
pub(crate) mod errors;
pub mod header;
pub(crate) mod index;
pub(crate) mod registry;
pub(crate) mod trailer;
pub(crate) mod varint;
pub mod writer;

#[cfg(test)]
mod test_fixtures;

// Public re-exports keep the typed error API reachable for
// downstream code without naming the (now-`pub(crate)`) `errors`
// submodule. The `*Kind` sub-enums must be re-exported too because
// the top-level error variants carry them as `#[source]`.
pub use errors::{
    BlockHeaderInvariantKind, InvalidRecordKind, PhaseChainMarkerInconsistencyKind, PspReadError,
    PspWriteError, ScalarDecodeError, TomlSerError, VarintError,
};
