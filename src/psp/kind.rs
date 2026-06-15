//! The `.psp` container *schema* abstraction — what the generic
//! container core (blocking, the block index, per-column zstd, the
//! self-describing header) cannot know on its own.
//!
//! A `.psp` "kind" (`snp` today; `ssr` next) supplies the four things
//! the otherwise schema-agnostic machinery needs (architecture §10.2):
//! its **column table**, a **record type**, the **per-block
//! accumulator** that buffers records into column buffers, and the
//! **record → column-bytes** encode. The *decode* half (columns →
//! record) is added in step 1b alongside the reader's typed path; this
//! module is the writer-side surface only.
//!
//! The generic [`PspWriter`](super::writer::PspWriter) is parameterised
//! on a [`PspKind`]; [`SnpKind`](super::writer::SnpKind) is the only
//! implementation in steps 1–3. The flush/encode machinery iterates
//! [`PspKind::columns`] and calls [`PspKind::encode_column`] rather than
//! closing over the hardcoded `V1_0_COLUMNS` + `PileupRecord`, so a
//! second kind is "a table + a record mapping" (§10.4), not a rewrite.

use std::io::Read;

use super::ColumnDef;
use super::block::BlockHeader;
use super::errors::{PspReadError, PspWriteError};

/// One `.psp` kind (`snp` | `ssr`): the schema the generic container
/// core needs beyond its block / index / header machinery.
///
/// The reader-side `record_from_columns` half is deliberately absent —
/// it lands in step 1b with the typed-record reader. Per architecture
/// §10.2 the heavy cross-sample path rides the *columnar* reader (which
/// is already schema-agnostic), so the typed decode can stay thin and
/// per-kind.
pub trait PspKind {
    /// The in-memory record this kind writes and reads (`PileupRecord`
    /// for SNP; the future `SsrLocusRecord` for SSR).
    type Record;

    /// The per-block accumulator that buffers this kind's records into
    /// its column buffers (the **write** side).
    type Block: BlockAccumulator<Record = Self::Record>;

    /// The per-block decoder that turns decompressed columns back into
    /// records (the **read** side — mirror of [`Self::Block`]). Left
    /// **unbounded here** so the `pub` trait surface does not leak the
    /// `pub(crate)` [`BlockDecoder`] (or the wire types in its
    /// signatures); the bound is applied where the decoder is driven
    /// (the reader's `RecordsIter`).
    type Decoder;

    /// Header schema-family tag (architecture §10.3). Written into the
    /// TOML header and used to select the registry on read.
    const KIND: &'static str;

    /// This kind's column registry table — the `&'static [ColumnDef]`
    /// the writer walks to emit columns in tag order (SNP: `V1_0_COLUMNS`).
    fn columns() -> &'static [ColumnDef];

    /// Encode one column from a filled block to its uncompressed wire
    /// bytes, appending into `out` (which the caller has `clear()`ed).
    /// This is the per-schema exhaustive column dispatch — the SNP impl
    /// matches on `def.key` and folds in the writer-side column-size
    /// self-check.
    fn encode_column(
        def: &ColumnDef,
        block: &Self::Block,
        out: &mut Vec<u8>,
    ) -> Result<(), PspWriteError>;

    /// The `(chrom_id, position)` of a record, used by the reader to
    /// apply a region-query window. SNP returns the record's point
    /// position; the interval generalisation (architecture §10.5, a
    /// later step) widens the right edge to `record.end`.
    fn record_coord(record: &Self::Record) -> (u32, u32);
}

/// The per-block accumulator a [`PspKind`] fills as records arrive: its
/// column buffers plus the structural metadata the generic flush and
/// block index read back.
///
/// `chrom_id` / `first_pos` are exposed as accessors (rather than owned
/// by the writer) because the accumulator already holds them — SNP needs
/// `first_pos` for delta-position encoding — so a single source of truth
/// is cheaper than duplicating that state into the writer.
pub trait BlockAccumulator {
    /// The record type this accumulator appends (matches
    /// [`PspKind::Record`]).
    type Record;

    /// Open a new block at `(chrom_id, first_pos)`. The generic writer
    /// calls this when no block is open and the next record arrives.
    fn new_block(chrom_id: u32, first_pos: u32) -> Self;

    /// Append one record's content to the open block's column buffers.
    fn append(&mut self, record: &Self::Record);

    /// The block's chromosome id (constant for the block's lifetime).
    fn chrom_id(&self) -> u32;

    /// The block's first record position — the block-index `first_pos`
    /// and the windowing reference.
    fn first_pos(&self) -> u32;

    /// The block-index `last_pos`: SNP = last record start; SSR =
    /// `max(record.end)` over the block (architecture §10.5).
    fn last_pos(&self) -> u32;

    /// Number of per-record entries buffered (the block header's
    /// `n_records`).
    fn n_records(&self) -> u32;

    /// Number of per-entry entries buffered — the per-allele total for
    /// SNP (the block header's `n_total_alleles`); the CSR entry total
    /// for SSR.
    fn n_entries(&self) -> u32;

    /// Rough projection of the uncompressed byte total, for the
    /// size-driven auto-flush safety cap.
    fn projected_bytes(&self) -> usize;
}

/// The per-block **decoder** a kind drives on the read path — the
/// mirror of [`BlockAccumulator`] (architecture §10.2). The generic
/// reader (`RecordsIter`) owns the block framing and the shared
/// decompression scratch; the decoder owns the schema-specific decode
/// state (its decoded columns + any cross-block reuse buffers) and the
/// per-block record cursor.
///
/// `pub(crate)` on purpose: it names internal wire types
/// ([`BlockHeader`]) and the zstd context in its signatures, so it must
/// not become part of the public API. [`PspKind::Decoder`] is therefore
/// left unbounded in the `pub` trait, and the `Decoder: BlockDecoder`
/// bound is applied only where the reader instantiates it.
pub(crate) trait BlockDecoder {
    /// The record this decoder materialises (matches
    /// [`PspKind::Record`]).
    type Record;

    /// A fresh decoder with empty reuse buffers and no block loaded.
    fn new_decoder() -> Self;

    /// Decode one block's columns into the decoder, replacing any
    /// previously-loaded block and resetting the record cursor. The
    /// generic caller has already positioned `source` at the block
    /// payload (just past its header) and supplies the shared
    /// decompression scratch; `budget` bounds the decompressed size.
    fn decode_block<R: Read>(
        &mut self,
        source: &mut R,
        header: &BlockHeader,
        budget: u64,
        decompressor: &mut zstd::bulk::Decompressor<'static>,
        compressed_scratch: &mut Vec<u8>,
        decompressed_scratch: &mut Vec<u8>,
    ) -> Result<(), PspReadError>;

    /// Materialise the next record of the loaded block, advancing the
    /// cursor. `None` once the block is exhausted (or before any block
    /// is loaded). Mirrors the writer's `append` in reverse.
    fn next_record(&mut self) -> Option<Result<Self::Record, PspReadError>>;
}
