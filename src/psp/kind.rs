//! The `.psp` container *schema* abstraction — what the generic
//! container core (blocking, the block index, per-column zstd, the
//! self-describing header) cannot know on its own.
//!
//! A `.psp` "kind" (`snp` today; `ssr` next) supplies the four things
//! the otherwise schema-agnostic machinery needs (architecture §10.2):
//! its **column table**, a **record type**, the **per-block
//! accumulator** that buffers records into column buffers, and the
//! **record → column-bytes** encode — plus the read-side mirror, the
//! per-block [`BlockDecoder`] (columns → record), which is `pub(crate)`
//! because it names internal wire types.
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

/// Declare a per-schema column-dispatch key enum with a compile-time
/// single source of truth for the `key ↔ tag` pairing (M5).
///
/// Both directions — `tag(self) -> u16` and `from_tag(u16) ->
/// Option<Self>` — are generated from one `Variant = tag` list as
/// exhaustive `match`es, so a new variant *must* appear in the list to
/// compile (it can't be silently omitted from a hand-written lookup
/// array, which is the hole that a separate `from_tag` reintroduced).
/// Tags must be integer literals so they double as `match` patterns.
/// The enum is deliberately **not** `#[non_exhaustive]` — exhaustive
/// dispatch at every call site is the point.
macro_rules! column_key {
    (
        $(#[$enum_meta:meta])*
        $vis:vis enum $name:ident {
            $( $(#[$var_meta:meta])* $variant:ident = $tag:literal ),+ $(,)?
        }
    ) => {
        $(#[$enum_meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        $vis enum $name {
            $( $(#[$var_meta])* $variant ),+
        }

        impl $name {
            /// The column tag for this key. Exhaustive `match Self`, so
            /// a new variant forces an arm here.
            pub const fn tag(self) -> u16 {
                match self { $( Self::$variant => $tag ),+ }
            }

            /// Map a column tag to its dispatch key, or `None` for a
            /// tag this schema does not define. Generated from the same
            /// variant↔tag list as [`Self::tag`], so the pairing is a
            /// compile-time bijection — a forgotten variant cannot
            /// compile.
            pub fn from_tag(tag: u16) -> Option<Self> {
                match tag {
                    $( $tag => Some(Self::$variant), )+
                    _ => None,
                }
            }
        }
    };
}
pub(crate) use column_key;

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

    /// The `(chrom_id, start, end)` interval a record occupies, used by
    /// the reader to apply a region-query window (architecture §10.5).
    /// `end` is **exclusive**. SNP is a degenerate point: `(chrom_id,
    /// pos, pos + 1)`. SSR returns the locus interval. A query `[q_start,
    /// q_end]` (inclusive) keeps a record iff `start <= q_end` and `end >
    /// q_start`.
    fn record_interval(record: &Self::Record) -> (u32, u32, u32);
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

    /// Release the decoded block's per-block column buffers once it is
    /// exhausted, before the generic reader advances to the next block.
    /// Cross-block *reuse* buffers (e.g. the SNP CSR slabs) may be kept
    /// — they are overwritten on the next `decode_block` — but the
    /// large per-block columns should be freed so a held-but-idle
    /// iterator does not pin the last block's memory (the project's
    /// RAM-for-scaling thesis). After `unload`, `next_record` must
    /// return `None` until the next `decode_block`.
    fn unload(&mut self);
}
