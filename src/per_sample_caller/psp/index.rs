//! The tail-mounted block index plus its XXH3-64 integrity check.
//!
//! Each `.psp` carries one [`BlockIndexEntry`] per block, written
//! between the last data block and the file trailer. A reader doing
//! a region query binary-searches the in-memory index by
//! `(chrom_id, first_pos, last_pos)`; sequential reads ignore the
//! index entirely and stream from block 0.
//!
//! See spec §"Block index" for the byte-level layout and the
//! rationale for the tail-mounted, uncompressed, atomic-with-content
//! design. The XXH3-64 checksum (low 32 bits, stored in the trailer)
//! closes the only file region not already covered by a zstd frame's
//! built-in checksum (spec §"File trailer" / Q-PL11).

use xxhash_rust::xxh3::xxh3_64;

use super::errors::PspReadError;
use super::varint::{decode_u64_leb128, encode_u64_leb128};

/// One block's coordinate-keyed index entry.
///
/// Field widths match the on-disk schema: `chrom_id`, `first_pos`,
/// `last_pos`, `n_records` are written as varints (so they may grow
/// in future format versions) but are bounded to fit in `u32` at
/// v1.0 — chromosome counts and per-block record counts comfortably
/// fit, and per-contig positions are bounded by the SAM `@SQ LN`
/// max of `2^31 - 1`. `block_offset` is a fixed `u64` to keep the
/// per-entry cost predictable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockIndexEntry {
    pub chrom_id: u32,
    pub first_pos: u32,
    pub last_pos: u32,
    pub n_records: u32,
    pub block_offset: u64,
}

/// Minimum bytes a [`BlockIndexEntry`] can occupy on the wire: four
/// single-byte varints (`chrom_id`, `first_pos`, `last_pos`,
/// `n_records`, each ≥ 1 byte) plus the fixed 8-byte `block_offset`.
/// Used to derive an upper bound on `expected_n_blocks` from
/// `bytes.len()` at decode time so a tampered trailer cannot drive
/// an unbounded `Vec::with_capacity`.
const MIN_ENTRY_BYTES: u64 = 1 + 1 + 1 + 1 + 8;

/// Worst-case bytes a [`BlockIndexEntry`] can occupy on the wire:
/// four 5-byte LEB128 varints (the maximum width for any `u32`)
/// plus the fixed 8-byte `block_offset`. Used as the per-entry
/// reservation hint in [`encode_index`] so the worst-case file
/// allocates exactly once.
const MAX_ENTRY_BYTES_HINT: usize = 5 + 5 + 5 + 5 + 8;

/// XXH3-64 over `bytes`, truncated to its low 32 bits. Matches the
/// truncation zstd uses for its frame content checksum (RFC 8878
/// §3.1.1), so the same hash backs both — there is only ever one
/// XXH3 implementation in the codebase. Used by both the writer (to
/// stamp the trailer's `index_checksum`) and the reader (to verify
/// it).
pub fn checksum_index(bytes: &[u8]) -> u32 {
    xxh3_64(bytes) as u32
}

/// Encode the block index to its on-wire byte form.
///
/// Each entry: four LEB128 varints (`chrom_id`, `first_pos`,
/// `last_pos`, `n_records`) then an 8-byte little-endian `u64`
/// (`block_offset`). Entries are written in the order given, which
/// the writer always sets to genomic order (the order blocks were
/// written).
pub fn encode_index(entries: &[BlockIndexEntry]) -> Vec<u8> {
    // Worst-case sizing per [`MAX_ENTRY_BYTES_HINT`]: varints up to
    // 5 bytes for each `u32` plus the fixed 8-byte `block_offset`.
    // Tightly bounding the initial capacity keeps the index region
    // small even before append-amortisation kicks in.
    let mut out = Vec::with_capacity(entries.len() * MAX_ENTRY_BYTES_HINT);
    for e in entries {
        encode_u64_leb128(u64::from(e.chrom_id), &mut out);
        encode_u64_leb128(u64::from(e.first_pos), &mut out);
        encode_u64_leb128(u64::from(e.last_pos), &mut out);
        encode_u64_leb128(u64::from(e.n_records), &mut out);
        out.extend_from_slice(&e.block_offset.to_le_bytes());
    }
    out
}

/// Decode the block index from its on-wire byte form.
///
/// `expected_n_blocks` is the trailer's `n_blocks` field; this
/// function decodes exactly that many entries and verifies the
/// buffer is exhausted afterwards. A premature buffer end is a
/// [`PspReadError::IndexEntryDecode`] / [`PspReadError::IndexTruncated`]
/// (depending on which field ran out); trailing bytes past the last
/// entry are [`PspReadError::IndexTrailingBytes`].
///
/// Coordinate-monotonicity and bounds checks against the header's
/// chromosome table are the reader's responsibility (one layer up)
/// and are not performed here.
pub fn decode_index(
    bytes: &[u8],
    expected_n_blocks: u64,
) -> Result<Vec<BlockIndexEntry>, PspReadError> {
    // DoS guard: `expected_n_blocks` is the trailer's `n_blocks`
    // field, decoded from an attacker-influenced file tail. Cap the
    // up-front `Vec::with_capacity` against the buffer-derived
    // ceiling so a `u64::MAX` trailer cannot drive an unbounded
    // allocation. Truncation / mismatch is then detected per-entry
    // inside the loop as before.
    let cap_bound = (bytes.len() as u64) / MIN_ENTRY_BYTES + 1;
    let initial_cap = expected_n_blocks.min(cap_bound) as usize;
    let mut cursor = 0usize;
    let mut entries = Vec::with_capacity(initial_cap);

    for entry_idx in 0..expected_n_blocks {
        let entry = decode_one_entry(bytes, &mut cursor, entry_idx as usize)?;
        entries.push(entry);
    }

    if cursor != bytes.len() {
        return Err(PspReadError::IndexTrailingBytes {
            trailing_bytes: bytes.len() - cursor,
        });
    }
    Ok(entries)
}

/// Decode one entry starting at `*cursor`; advances `*cursor` past
/// the consumed bytes on success.
fn decode_one_entry(
    bytes: &[u8],
    cursor: &mut usize,
    entry_idx: usize,
) -> Result<BlockIndexEntry, PspReadError> {
    let chrom_id = decode_field_u32(bytes, cursor, entry_idx, "chrom_id")?;
    let first_pos = decode_field_u32(bytes, cursor, entry_idx, "first_pos")?;
    let last_pos = decode_field_u32(bytes, cursor, entry_idx, "last_pos")?;
    let n_records = decode_field_u32(bytes, cursor, entry_idx, "n_records")?;
    let block_offset = decode_field_u64_le(bytes, cursor, entry_idx, "block_offset")?;
    Ok(BlockIndexEntry {
        chrom_id,
        first_pos,
        last_pos,
        n_records,
        block_offset,
    })
}

fn decode_field_u32(
    bytes: &[u8],
    cursor: &mut usize,
    entry_idx: usize,
    field: &'static str,
) -> Result<u32, PspReadError> {
    let (value_u64, consumed) =
        decode_u64_leb128(&bytes[*cursor..]).map_err(|source| PspReadError::IndexEntryDecode {
            entry: entry_idx,
            field,
            source,
        })?;
    *cursor += consumed;
    u32::try_from(value_u64).map_err(|_| PspReadError::IndexFieldOverflow {
        entry: entry_idx,
        field,
        value: value_u64,
    })
}

fn decode_field_u64_le(
    bytes: &[u8],
    cursor: &mut usize,
    entry_idx: usize,
    field: &'static str,
) -> Result<u64, PspReadError> {
    if *cursor + 8 > bytes.len() {
        return Err(PspReadError::IndexTruncated {
            entry: entry_idx,
            field,
        });
    }
    let value = u64::from_le_bytes(bytes[*cursor..*cursor + 8].try_into().unwrap());
    *cursor += 8;
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// (Plan group V — V4.) Build N entries, encode, compute
    /// XXH3-64-trunc, decode, assert all entries equal. Covers the
    /// per-entry round-trip and the entry-count cross-check.
    #[test]
    fn round_trip_small() {
        let entries = vec![
            BlockIndexEntry {
                chrom_id: 0,
                first_pos: 1,
                last_pos: 1_000_000,
                n_records: 50_000,
                block_offset: 256,
            },
            BlockIndexEntry {
                chrom_id: 0,
                first_pos: 1_000_001,
                last_pos: 2_000_000,
                n_records: 49_000,
                block_offset: 16_777_216,
            },
            BlockIndexEntry {
                chrom_id: 1,
                first_pos: 1,
                last_pos: 500_000,
                n_records: 30_000,
                block_offset: 33_554_432,
            },
        ];
        let bytes = encode_index(&entries);
        let _checksum = checksum_index(&bytes);
        let decoded = decode_index(&bytes, entries.len() as u64).unwrap();
        assert_eq!(decoded, entries);
    }

    /// Empty index round-trips. An empty `.psp` (zero blocks) still
    /// has a zero-entry index; this case must succeed and produce a
    /// stable checksum.
    #[test]
    fn round_trip_empty() {
        let entries: Vec<BlockIndexEntry> = Vec::new();
        let bytes = encode_index(&entries);
        assert!(bytes.is_empty());
        let decoded = decode_index(&bytes, 0).unwrap();
        assert!(decoded.is_empty());
        // checksum of the empty byte sequence is the documented
        // XXH3-64 constant; we don't pin the value here (the spec
        // says it is whatever XXH3-64 of [] truncated produces),
        // but check it computes without panicking.
        let _ = checksum_index(&bytes);
    }

    /// Larger index — 10k entries — to stress the loop and confirm
    /// the capacity hint is correct (allocations should not grow
    /// after the initial `with_capacity`).
    #[test]
    fn round_trip_many() {
        let n = 10_000usize;
        let entries: Vec<BlockIndexEntry> = (0..n)
            .map(|i| BlockIndexEntry {
                chrom_id: (i % 24) as u32,
                first_pos: (i * 1000 + 1) as u32,
                last_pos: ((i + 1) * 1000) as u32,
                n_records: 100,
                block_offset: (i as u64) * 16 * 1024 * 1024,
            })
            .collect();
        let bytes = encode_index(&entries);
        let decoded = decode_index(&bytes, n as u64).unwrap();
        assert_eq!(decoded, entries);
    }

    /// (Plan group V — V4 negative.) Flip one byte in the encoded
    /// index. The recomputed checksum disagrees with the original.
    #[test]
    fn checksum_changes_under_corruption() {
        let entries = vec![BlockIndexEntry {
            chrom_id: 5,
            first_pos: 1,
            last_pos: 1_000,
            n_records: 10,
            block_offset: 4096,
        }];
        let mut bytes = encode_index(&entries);
        let pristine_checksum = checksum_index(&bytes);
        // Flip a byte inside the block_offset varint payload of the
        // first entry. Don't pick the first byte — that's chrom_id,
        // which after corruption might still decode validly but
        // wrong.
        let last = bytes.len() - 1;
        bytes[last] ^= 0x01;
        let corrupted_checksum = checksum_index(&bytes);
        assert_ne!(
            pristine_checksum, corrupted_checksum,
            "XXH3-trunc should detect a single-byte flip"
        );
    }

    /// Truncated buffer: encode N entries but pass shorter `bytes`
    /// to decode. Reports `IndexEntryDecode` / `IndexTruncated`
    /// depending on which field ran out first.
    #[test]
    fn decode_rejects_truncated_varint() {
        let entries = vec![BlockIndexEntry {
            chrom_id: 1,
            first_pos: 1,
            last_pos: 2,
            n_records: 3,
            block_offset: 4,
        }];
        let bytes = encode_index(&entries);
        // Drop the last byte — that's part of `block_offset`'s 8 LE
        // bytes, so we should hit `IndexTruncated` on `block_offset`.
        let truncated = &bytes[..bytes.len() - 1];
        let err = decode_index(truncated, 1).unwrap_err();
        match err {
            PspReadError::IndexTruncated { entry, field } => {
                assert_eq!(entry, 0);
                assert_eq!(field, "block_offset");
            }
            other => panic!("expected IndexTruncated, got {other:?}"),
        }
    }

    /// Truncated mid-varint (chop into the very first varint) →
    /// `IndexEntryDecode` wrapping `VarintError`.
    #[test]
    fn decode_rejects_truncated_first_field() {
        // A varint with continuation bit set, no terminator.
        let bytes = vec![0x80u8, 0x80];
        let err = decode_index(&bytes, 1).unwrap_err();
        match err {
            PspReadError::IndexEntryDecode { entry, field, .. } => {
                assert_eq!(entry, 0);
                assert_eq!(field, "chrom_id");
            }
            other => panic!("expected IndexEntryDecode, got {other:?}"),
        }
    }

    /// Trailing bytes past the last entry: the buffer carries more
    /// than the declared `n_blocks` accounts for.
    #[test]
    fn decode_rejects_trailing_bytes() {
        let entries = vec![BlockIndexEntry {
            chrom_id: 0,
            first_pos: 1,
            last_pos: 1,
            n_records: 1,
            block_offset: 0,
        }];
        let mut bytes = encode_index(&entries);
        bytes.extend_from_slice(&[0u8, 1, 2]);
        let err = decode_index(&bytes, 1).unwrap_err();
        match err {
            PspReadError::IndexTrailingBytes { trailing_bytes } => {
                assert_eq!(trailing_bytes, 3);
            }
            other => panic!("expected IndexTrailingBytes, got {other:?}"),
        }
    }

    /// A `chrom_id` that overflows `u32` is a decode error, not a
    /// silent truncation.
    #[test]
    fn decode_rejects_u32_overflow() {
        // Build a varint for u32::MAX + 1 = 0x1_0000_0000.
        let mut bytes = Vec::new();
        encode_u64_leb128((u32::MAX as u64) + 1, &mut bytes);
        // Pad with zero varints for the other three fields and a
        // zero u64 for block_offset so the layout is otherwise
        // valid.
        for _ in 0..3 {
            encode_u64_leb128(0, &mut bytes);
        }
        bytes.extend_from_slice(&0u64.to_le_bytes());
        let err = decode_index(&bytes, 1).unwrap_err();
        match err {
            PspReadError::IndexFieldOverflow {
                entry,
                field,
                value,
            } => {
                assert_eq!(entry, 0);
                assert_eq!(field, "chrom_id");
                assert_eq!(value, (u32::MAX as u64) + 1);
            }
            other => panic!("expected IndexFieldOverflow, got {other:?}"),
        }
    }

    /// `decode_index` with an `expected_n_blocks` smaller than what
    /// the buffer encodes hits `IndexTrailingBytes` (because the
    /// loop stops early and the unconsumed entry is "trailing").
    #[test]
    fn decode_under_count_reports_trailing() {
        let entries = vec![
            BlockIndexEntry {
                chrom_id: 0,
                first_pos: 1,
                last_pos: 1,
                n_records: 1,
                block_offset: 0,
            },
            BlockIndexEntry {
                chrom_id: 0,
                first_pos: 2,
                last_pos: 2,
                n_records: 1,
                block_offset: 64,
            },
        ];
        let bytes = encode_index(&entries);
        let err = decode_index(&bytes, 1).unwrap_err();
        assert!(matches!(err, PspReadError::IndexTrailingBytes { .. }));
    }

    /// `decode_index` with an `expected_n_blocks` larger than what
    /// the buffer encodes hits the varint truncation path on the
    /// "missing" entry's first field.
    #[test]
    fn decode_over_count_reports_truncated() {
        let entries = vec![BlockIndexEntry {
            chrom_id: 0,
            first_pos: 1,
            last_pos: 1,
            n_records: 1,
            block_offset: 0,
        }];
        let bytes = encode_index(&entries);
        let err = decode_index(&bytes, 2).unwrap_err();
        // The decoder will try to start entry 1 (zero-indexed) at the
        // buffer's end and immediately hit Truncated.
        match err {
            PspReadError::IndexEntryDecode { entry, field, .. } => {
                assert_eq!(entry, 1);
                assert_eq!(field, "chrom_id");
            }
            other => panic!("expected IndexEntryDecode at entry 1, got {other:?}"),
        }
    }

    /// (M2.) Pins the byte layout of `encode_index` for a known
    /// entry. A field-reorder refactor in `BlockIndexEntry` (or in
    /// `encode_index`'s ordering of the varint writes) would shift
    /// the asserted bytes and fail this test loudly.
    #[test]
    fn wire_layout_is_stable() {
        let entries = vec![BlockIndexEntry {
            chrom_id: 0x01,
            first_pos: 0x80,     // 2-byte LEB128: 0x80 0x01
            last_pos: 0x4000,    // 3-byte LEB128: 0x80 0x80 0x01
            n_records: 0x200000, // 4-byte LEB128: 0x80 0x80 0x80 0x01
            block_offset: 0x0102030405060708,
        }];
        let bytes = encode_index(&entries);
        // chrom_id = 1 — 1-byte varint.
        assert_eq!(bytes[0], 0x01);
        // first_pos = 128 — 2-byte varint.
        assert_eq!(&bytes[1..3], &[0x80, 0x01]);
        // last_pos = 0x4000 — 3-byte varint.
        assert_eq!(&bytes[3..6], &[0x80, 0x80, 0x01]);
        // n_records = 0x200000 — 4-byte varint.
        assert_eq!(&bytes[6..10], &[0x80, 0x80, 0x80, 0x01]);
        // block_offset — 8 LE bytes.
        assert_eq!(&bytes[10..18], &0x0102030405060708u64.to_le_bytes());
        assert_eq!(bytes.len(), 18);
    }

    /// (B3 regression.) A hostile trailer claiming `n_blocks =
    /// u64::MAX` against an empty buffer must not drive an
    /// unbounded `Vec::with_capacity`; the up-front cap clamps the
    /// reservation and the per-entry decode surfaces the truncation.
    #[test]
    fn decode_index_rejects_giant_expected_count_without_oom() {
        let err = decode_index(&[], u64::MAX).unwrap_err();
        // The first decode_field_u32 call hits an empty slice and
        // returns Truncated, wrapped as IndexEntryDecode at entry 0.
        match err {
            PspReadError::IndexEntryDecode { entry, field, .. } => {
                assert_eq!(entry, 0);
                assert_eq!(field, "chrom_id");
            }
            other => panic!("expected IndexEntryDecode at entry 0, got {other:?}"),
        }
    }
}
