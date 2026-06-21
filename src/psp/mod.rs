//! `.psp` generic columnar container — writer and reader.
//!
//! Originally the SNP per-sample pileup byte format, this module is now
//! a generic columnar container (architecture §10) hosting two schemas
//! — `snp` (the Stage-1 per-sample pileup, one record per covered
//! reference position) and `ssr` (per-locus microsatellite evidence).
//! It is the bridge between those in-memory records and the on-disk
//! artefact downstream stages consume. The byte layout is specified in
//! `ia/specs/per_sample_pileup_format.md`; this module is its
//! authoritative implementation.
//!
//! Module layout (most sub-modules are `pub(crate)` — accessed
//! through the [`writer::PspWriter`] / [`reader::PspReader`] entry
//! points and the public [`header`] type surface):
//!
//! - [`kind`]: the cross-schema [`kind::PspKind`] / `BlockAccumulator`
//!   / `BlockDecoder` trait surface a schema implements.
//! - `registry`: the v1.0 **SNP** column-tag registry, and the
//!   kind→registry dispatch (`columns_for_kind`); `registry_ssr` is the
//!   SSR counterpart. Each is the single source of truth for its
//!   schema's columns, shared by writer, reader, and `psp_spec_dump`.
//! - `varint`: LEB128 / zig-zag-LEB128 codecs.
//! - `errors`: typed error enums (re-exported below).
//! - [`header`], `block`, `index`, `trailer`: per-section wire codecs.
//! - [`writer`]: streaming [`writer::PspWriter`] over any `Write` sink.
//! - [`reader`]: streaming [`reader::PspReader`] over any `Read +
//!   Seek` source — sequential iteration over every record, plus
//!   coordinate-keyed `region_records` queries via the block index.

// Mi5: most submodules are crate-internal scaffolding for the
// writer/reader implementations. The bench imports `header::*` and
// `writer::PspWriter`; everything else is reachable only through
// the re-exports below or through the writer's API surface.
pub(crate) mod block;
pub(crate) mod errors;
pub mod header;
pub(crate) mod index;
pub mod kind;
pub mod reader;
pub(crate) mod registry;
pub(crate) mod registry_ssr;
pub(crate) mod trailer;
pub(crate) mod varint;
pub mod writer;

#[cfg(test)]
pub(crate) mod test_fixtures;

// Public re-exports keep the typed error API reachable for
// downstream code without naming the (now-`pub(crate)`) `errors`
// submodule. The `*Kind` sub-enums must be re-exported too because
// the top-level error variants carry them as `#[source]`.
pub use errors::{
    BlockHeaderInvariantKind, InvalidRecordKind, PspReadError, PspWriteError, ScalarDecodeError,
    TomlSerError, VarintError,
};
pub use index::BlockIndexEntry;
pub use reader::{BlockColumnReader, BlockColumns, OwnedRecordsIter, PspReader, RecordsIter};
// The container-schema abstraction (architecture §10) puts `ColumnDef`
// in the `pub` [`kind::PspKind`] signatures, so it must be reachable at
// at least that visibility. Re-exported here from the (otherwise
// `pub(crate)`) registry rather than widening the whole registry.
pub use registry::ColumnDef;
