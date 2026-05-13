//! Per-block primitives: block header, per-shape column codecs, zstd
//! compression.
//!
//! This module owns the wire-level shape of a single block:
//!
//! - The uncompressed block header (varint stream).
//! - The per-column manifest entries (varint stream).
//! - The per-shape encoders and decoders that turn typed value lists
//!   into bytes and back (`encode_scalar_column` /
//!   `decode_scalar_column`, plus the list and bytes variants).
//! - Thin wrappers over the `zstd` crate for compressing and
//!   decompressing the per-column payloads.
//!
//! What lives one layer up (writer.rs / reader.rs):
//!
//! - The `BlockAccumulator` that buffers records into per-column
//!   buffers and decides when to flush.
//! - The dispatch from `[[column]]` schema entries to the
//!   type-specific codec calls here.
//! - The cross-block phase-chain consistency checks.
//!
//! Each block on the wire:
//!
//! ```text
//! block:
//!   <block_header>            uncompressed varint stream
//!   <column_payloads>         concatenated, each independently
//!                             zstd-compressed; the manifest in the
//!                             block_header gives byte lengths
//! ```
//!
//! See spec §"Block layout" for the canonical definitions.

use std::io;

use super::errors::{PspReadError, ScalarDecodeError, VarintError};
use super::varint::{decode_u64_leb128, encode_u64_leb128};

// ---------------------------------------------------------------------
// WireScalar trait — fixed-width per-entry codecs
// ---------------------------------------------------------------------

/// Encode and decode a single fixed-width scalar value to / from
/// little-endian bytes. Implemented for `u8`, `u16`, `u32`, `u64`,
/// `i32`, `i64`, `f32`, `f64`, `bool`. **Not** implemented for the
/// variable-width `varint` / `svarint` element types — those use
/// the standalone codecs in [`crate::per_sample_caller::psp::varint`]
/// plus the [`encode_varint_column`] / [`decode_varint_column`]
/// wrappers here.
pub trait WireScalar: Sized + Copy {
    /// Encode self appended to `out`. Always pushes exactly
    /// [`Self::FIXED_BYTE_WIDTH`] bytes.
    fn encode_le(self, out: &mut Vec<u8>);

    /// Decode one value from the front of `bytes`. Returns the value
    /// and the number of bytes consumed (always
    /// [`Self::FIXED_BYTE_WIDTH`]).
    fn decode_le(bytes: &[u8]) -> Result<(Self, usize), ScalarDecodeError>;

    /// Byte width of one encoded element. Pinned at the type level
    /// so column-size predictions can be made statically.
    const FIXED_BYTE_WIDTH: usize;
}

macro_rules! impl_wire_scalar_le {
    ($t:ty, $width:literal) => {
        impl WireScalar for $t {
            #[inline]
            fn encode_le(self, out: &mut Vec<u8>) {
                out.extend_from_slice(&self.to_le_bytes());
            }
            #[inline]
            fn decode_le(bytes: &[u8]) -> Result<(Self, usize), ScalarDecodeError> {
                if bytes.len() < $width {
                    return Err(ScalarDecodeError::Truncated);
                }
                let arr: [u8; $width] = bytes[..$width].try_into().unwrap();
                Ok((<$t>::from_le_bytes(arr), $width))
            }
            const FIXED_BYTE_WIDTH: usize = $width;
        }
    };
}

impl_wire_scalar_le!(u16, 2);
impl_wire_scalar_le!(u32, 4);
impl_wire_scalar_le!(u64, 8);
impl_wire_scalar_le!(i32, 4);
impl_wire_scalar_le!(i64, 8);
impl_wire_scalar_le!(f32, 4);
impl_wire_scalar_le!(f64, 8);

impl WireScalar for u8 {
    #[inline]
    fn encode_le(self, out: &mut Vec<u8>) {
        out.push(self);
    }
    #[inline]
    fn decode_le(bytes: &[u8]) -> Result<(Self, usize), ScalarDecodeError> {
        bytes
            .first()
            .copied()
            .map(|b| (b, 1))
            .ok_or(ScalarDecodeError::Truncated)
    }
    const FIXED_BYTE_WIDTH: usize = 1;
}

impl WireScalar for bool {
    #[inline]
    fn encode_le(self, out: &mut Vec<u8>) {
        out.push(u8::from(self));
    }
    #[inline]
    fn decode_le(bytes: &[u8]) -> Result<(Self, usize), ScalarDecodeError> {
        match bytes.first().copied() {
            None => Err(ScalarDecodeError::Truncated),
            Some(0) => Ok((false, 1)),
            Some(1) => Ok((true, 1)),
            Some(other) => Err(ScalarDecodeError::InvalidBool(other)),
        }
    }
    const FIXED_BYTE_WIDTH: usize = 1;
}

// ---------------------------------------------------------------------
// Scalar column codecs — fixed-width
// ---------------------------------------------------------------------

/// Encode a fixed-width scalar column.
///
/// Output: `values.len() * T::FIXED_BYTE_WIDTH` bytes, end to end.
/// The caller is responsible for getting `values.len()` right (it
/// must equal `n_records` for per-record columns and
/// `n_total_alleles` for per-allele columns).
pub fn encode_scalar_column<T: WireScalar>(values: &[T], out: &mut Vec<u8>) {
    out.reserve(values.len() * T::FIXED_BYTE_WIDTH);
    for &v in values {
        v.encode_le(out);
    }
}

/// Decode a fixed-width scalar column from a complete payload
/// buffer. `expected_count` is the cardinality-derived count
/// (`n_records` or `n_total_alleles`).
///
/// Errors:
/// - [`PspReadError::ColumnElementDecode`] if any element fails to
///   decode (buffer truncated mid-element; invalid bool byte).
/// - [`PspReadError::ColumnTrailingBytes`] if bytes remain past the
///   expected-count'th element.
/// - [`PspReadError::ColumnTruncated`] if the buffer ends before
///   `expected_count` elements have been decoded.
pub fn decode_scalar_column<T: WireScalar>(
    bytes: &[u8],
    expected_count: usize,
    column_name: &str,
) -> Result<Vec<T>, PspReadError> {
    let mut values = Vec::with_capacity(expected_count);
    let mut cursor = 0;
    for entry in 0..expected_count {
        if cursor >= bytes.len() && entry < expected_count {
            return Err(PspReadError::ColumnTruncated {
                column: column_name.to_string(),
                decoded: entry,
                expected: expected_count,
            });
        }
        let (v, n) =
            T::decode_le(&bytes[cursor..]).map_err(|source| PspReadError::ColumnElementDecode {
                column: column_name.to_string(),
                entry,
                source,
            })?;
        cursor += n;
        values.push(v);
    }
    if cursor != bytes.len() {
        return Err(PspReadError::ColumnTrailingBytes {
            column: column_name.to_string(),
            trailing: bytes.len() - cursor,
        });
    }
    Ok(values)
}

// ---------------------------------------------------------------------
// Varint scalar column codecs — variable-width
// ---------------------------------------------------------------------

/// Encode a `varint`-element-type scalar column: one LEB128-encoded
/// `u64` per value, packed end to end.
pub fn encode_varint_column(values: &[u64], out: &mut Vec<u8>) {
    // 1 byte for values 0..=127; cap at the 10-byte ceiling per
    // value as a worst-case capacity estimate.
    out.reserve(values.len());
    for &v in values {
        encode_u64_leb128(v, out);
    }
}

/// Decode a `varint`-element-type scalar column. Error shape mirrors
/// [`decode_scalar_column`].
pub fn decode_varint_column(
    bytes: &[u8],
    expected_count: usize,
    column_name: &str,
) -> Result<Vec<u64>, PspReadError> {
    let mut values = Vec::with_capacity(expected_count);
    let mut cursor = 0;
    for entry in 0..expected_count {
        if cursor >= bytes.len() {
            return Err(PspReadError::ColumnTruncated {
                column: column_name.to_string(),
                decoded: entry,
                expected: expected_count,
            });
        }
        let (v, n) =
            decode_u64_leb128(&bytes[cursor..]).map_err(|e| PspReadError::ColumnElementDecode {
                column: column_name.to_string(),
                entry,
                source: ScalarDecodeError::from(e),
            })?;
        cursor += n;
        values.push(v);
    }
    if cursor != bytes.len() {
        return Err(PspReadError::ColumnTrailingBytes {
            column: column_name.to_string(),
            trailing: bytes.len() - cursor,
        });
    }
    Ok(values)
}

// ---------------------------------------------------------------------
// List column codecs
// ---------------------------------------------------------------------

/// Encode a list column whose element type is a fixed-width
/// `WireScalar`.
///
/// Per entry: a `varint` count `k`, then `k` little-endian
/// `element-type` values packed end to end. The list of lists is
/// expected to have one entry per record (per-record list column) or
/// per allele (per-allele list column).
pub fn encode_list_column<T: WireScalar>(lists: &[&[T]], out: &mut Vec<u8>) {
    for list in lists {
        encode_u64_leb128(list.len() as u64, out);
        for &v in list.iter() {
            v.encode_le(out);
        }
    }
}

/// Decode a list column. `expected_count` is the cardinality-derived
/// per-entry count.
pub fn decode_list_column<T: WireScalar>(
    bytes: &[u8],
    expected_count: usize,
    column_name: &str,
) -> Result<Vec<Vec<T>>, PspReadError> {
    let mut out = Vec::with_capacity(expected_count);
    let mut cursor = 0;
    for entry in 0..expected_count {
        if cursor >= bytes.len() {
            return Err(PspReadError::ColumnTruncated {
                column: column_name.to_string(),
                decoded: entry,
                expected: expected_count,
            });
        }
        let (k_u64, n_k) =
            decode_u64_leb128(&bytes[cursor..]).map_err(|e| PspReadError::ColumnElementDecode {
                column: column_name.to_string(),
                entry,
                source: ScalarDecodeError::from(e),
            })?;
        cursor += n_k;
        let k = k_u64 as usize;
        let mut list = Vec::with_capacity(k);
        for _ in 0..k {
            let (v, n_v) = T::decode_le(&bytes[cursor..]).map_err(|source| {
                PspReadError::ColumnElementDecode {
                    column: column_name.to_string(),
                    entry,
                    source,
                }
            })?;
            cursor += n_v;
            list.push(v);
        }
        out.push(list);
    }
    if cursor != bytes.len() {
        return Err(PspReadError::ColumnTrailingBytes {
            column: column_name.to_string(),
            trailing: bytes.len() - cursor,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------
// Bytes column codecs
// ---------------------------------------------------------------------

/// Encode the bytes payload of a `bytes`-shape column: a flat
/// concatenation of every entry's bytes. The lengths-per-entry are
/// emitted separately, by the paired length column (a `(per-X,
/// scalar, varint)` column whose name lives in the schema's
/// `length-column` field).
pub fn encode_bytes_concat(values: &[&[u8]], out: &mut Vec<u8>) {
    let total: usize = values.iter().map(|v| v.len()).sum();
    out.reserve(total);
    for v in values {
        out.extend_from_slice(v);
    }
}

/// Decode the bytes payload of a `bytes`-shape column given the
/// already-decoded `lengths` from the paired length column. Returns
/// one `Vec<u8>` per entry.
///
/// Errors:
/// - [`PspReadError::ColumnTrailingBytes`] if the sum of `lengths`
///   is less than `bytes.len()`.
/// - [`PspReadError::ColumnTruncated`] if the sum of `lengths` is
///   greater than `bytes.len()`.
pub fn decode_bytes_split(
    bytes: &[u8],
    lengths: &[u64],
    column_name: &str,
) -> Result<Vec<Vec<u8>>, PspReadError> {
    let total: u64 = lengths.iter().sum();
    let buf_len = bytes.len() as u64;
    if total > buf_len {
        return Err(PspReadError::ColumnTruncated {
            column: column_name.to_string(),
            decoded: 0,
            expected: lengths.len(),
        });
    }
    if total < buf_len {
        return Err(PspReadError::ColumnTrailingBytes {
            column: column_name.to_string(),
            trailing: (buf_len - total) as usize,
        });
    }
    let mut out = Vec::with_capacity(lengths.len());
    let mut cursor = 0usize;
    for &len in lengths {
        let n = len as usize;
        out.push(bytes[cursor..cursor + n].to_vec());
        cursor += n;
    }
    Ok(out)
}

// ---------------------------------------------------------------------
// zstd helpers
// ---------------------------------------------------------------------

/// Default zstd compression level for `.psp` block column payloads.
/// Spec §"Compression" pins this at level 9. Not on the CLI.
pub const ZSTD_COMPRESSION_LEVEL: i32 = 9;

/// Construct a `zstd::bulk::Compressor` configured for `.psp` column
/// payloads: spec level 9 with frame content checksum enabled. Each
/// writer keeps one of these alive for its lifetime so the CCtx
/// workspace and tables are allocated once instead of per-column —
/// the per-call setup of a fresh encoder is the dominant driver of
/// the `__brk` / `__glibc_morecore` / `systrim` cluster the H3
/// finding identified in the samply profile.
pub fn new_column_compressor() -> io::Result<zstd::bulk::Compressor<'static>> {
    let mut c = zstd::bulk::Compressor::new(ZSTD_COMPRESSION_LEVEL)?;
    c.include_checksum(true)?;
    Ok(c)
}

/// Compress `input` with the persistent `compressor`, writing the
/// resulting zstd frame bytes into `out` (which is cleared first so
/// the same `Vec<u8>` can be reused across many calls without
/// reallocation). Output is one complete zstd frame with content
/// checksum enabled — the XXH64-truncated per-frame CRC the spec
/// relies on for block-level integrity (§"Compression" / Q-PL8).
pub fn zstd_compress_into(
    compressor: &mut zstd::bulk::Compressor<'static>,
    input: &[u8],
    out: &mut Vec<u8>,
) -> io::Result<()> {
    // `Compressor::compress_to_buffer` writes into the destination's
    // `spare_capacity_mut`, so the buffer must have enough free
    // capacity. Reserve the worst-case bound (`compress_bound(N)` is
    // ~N + N/255 + 16 bytes for level 9). For a 16 MiB block payload
    // that is ~16 MiB + ~65 KiB; the reserved capacity is reused on
    // the next call (zstd_compress_into only `clear`s, never
    // shrinks).
    out.clear();
    let bound = zstd::zstd_safe::compress_bound(input.len());
    out.reserve(bound);
    compressor.compress_to_buffer(input, out)?;
    Ok(())
}

/// Compress a column payload with zstd. Single-shot variant for
/// tests and other callers that do not have a persistent compressor.
/// Uses the same frame settings as [`zstd_compress_into`] so the
/// produced bytes decompress under the same reader path.
pub fn zstd_compress(input: &[u8]) -> io::Result<Vec<u8>> {
    let mut compressor = new_column_compressor()?;
    let mut buf = Vec::new();
    zstd_compress_into(&mut compressor, input, &mut buf)?;
    Ok(buf)
}

/// Decompress a zstd-compressed column payload. Failures wrap to
/// [`PspReadError::Zstd`] at the call site.
pub fn zstd_decompress(input: &[u8]) -> io::Result<Vec<u8>> {
    zstd::stream::decode_all(input)
}

// ---------------------------------------------------------------------
// Block header
// ---------------------------------------------------------------------

/// In-memory form of a block header.
///
/// Field widths mirror the wire form: counts are decoded into `u32`
/// to keep the public types ergonomic; the wire form uses varints so
/// future versions could widen.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlockHeader {
    pub chrom_id: u32,
    pub first_pos: u32,
    pub n_records: u32,
    pub n_total_alleles: u32,
    /// `active_chain_slots_at_block_start`: slot ids inherited from
    /// the immediately preceding block on the same chromosome.
    /// **Strictly ascending, no duplicates**, empty whenever this is
    /// the first block on its chromosome. Spec §"Phase-chain state
    /// across blocks".
    pub active_chain_slots: Vec<u16>,
    /// Per-column manifest, **strictly ascending by `column_tag`,
    /// no duplicates** — the writer sorts at emit time, the reader
    /// rejects out-of-order manifests.
    pub manifest: Vec<ColumnManifestEntry>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColumnManifestEntry {
    pub tag: u16,
    pub compressed_len: u32,
    pub uncompressed_len: u32,
}

/// Encode a block header to its varint-stream wire form. Validates
/// the in-memory invariants the decoder will also check, so a writer
/// that hands us a malformed header gets a hard error before any
/// byte is emitted.
pub fn encode_block_header(header: &BlockHeader, out: &mut Vec<u8>) -> Result<(), PspReadError> {
    validate_block_header_invariants(header)?;
    encode_u64_leb128(u64::from(header.chrom_id), out);
    encode_u64_leb128(u64::from(header.first_pos), out);
    encode_u64_leb128(u64::from(header.n_records), out);
    encode_u64_leb128(u64::from(header.n_total_alleles), out);
    encode_u64_leb128(header.active_chain_slots.len() as u64, out);
    for &slot in &header.active_chain_slots {
        out.extend_from_slice(&slot.to_le_bytes());
    }
    encode_u64_leb128(header.manifest.len() as u64, out);
    for entry in &header.manifest {
        encode_u64_leb128(u64::from(entry.tag), out);
        encode_u64_leb128(u64::from(entry.compressed_len), out);
        encode_u64_leb128(u64::from(entry.uncompressed_len), out);
    }
    Ok(())
}

/// Decode a block header from the start of `bytes`. Returns the
/// header and the number of bytes consumed. Validates the same
/// in-memory invariants the encoder enforces.
pub fn decode_block_header(bytes: &[u8]) -> Result<(BlockHeader, usize), PspReadError> {
    let mut cursor = 0usize;
    let chrom_id = decode_u32_field(bytes, &mut cursor, "chrom_id")?;
    let first_pos = decode_u32_field(bytes, &mut cursor, "first_pos")?;
    let n_records = decode_u32_field(bytes, &mut cursor, "n_records")?;
    let n_total_alleles = decode_u32_field(bytes, &mut cursor, "n_total_alleles")?;

    let k = decode_u32_field(bytes, &mut cursor, "active_chain_slots_count")?;
    let mut active_chain_slots = Vec::with_capacity(k as usize);
    for _ in 0..k {
        if cursor + 2 > bytes.len() {
            return Err(PspReadError::BlockHeaderField {
                field: "active_chain_slots",
                source: VarintError::Truncated,
            });
        }
        let arr: [u8; 2] = bytes[cursor..cursor + 2].try_into().unwrap();
        active_chain_slots.push(u16::from_le_bytes(arr));
        cursor += 2;
    }

    let n_columns = decode_u32_field(bytes, &mut cursor, "n_columns")?;
    let mut manifest = Vec::with_capacity(n_columns as usize);
    for _ in 0..n_columns {
        let tag = decode_u16_field(bytes, &mut cursor, "manifest.tag")?;
        let compressed_len = decode_u32_field(bytes, &mut cursor, "manifest.compressed_len")?;
        let uncompressed_len = decode_u32_field(bytes, &mut cursor, "manifest.uncompressed_len")?;
        manifest.push(ColumnManifestEntry {
            tag,
            compressed_len,
            uncompressed_len,
        });
    }

    let header = BlockHeader {
        chrom_id,
        first_pos,
        n_records,
        n_total_alleles,
        active_chain_slots,
        manifest,
    };
    validate_block_header_invariants(&header)?;
    Ok((header, cursor))
}

fn decode_u32_field(
    bytes: &[u8],
    cursor: &mut usize,
    field: &'static str,
) -> Result<u32, PspReadError> {
    let (v, n) = decode_u64_leb128(&bytes[*cursor..])
        .map_err(|source| PspReadError::BlockHeaderField { field, source })?;
    *cursor += n;
    u32::try_from(v).map_err(|_| PspReadError::BlockHeaderInvariant {
        reason: format!("{field}: value {v} exceeds u32 range"),
    })
}

fn decode_u16_field(
    bytes: &[u8],
    cursor: &mut usize,
    field: &'static str,
) -> Result<u16, PspReadError> {
    let (v, n) = decode_u64_leb128(&bytes[*cursor..])
        .map_err(|source| PspReadError::BlockHeaderField { field, source })?;
    *cursor += n;
    u16::try_from(v).map_err(|_| PspReadError::BlockHeaderInvariant {
        reason: format!("{field}: value {v} exceeds u16 range"),
    })
}

fn validate_block_header_invariants(header: &BlockHeader) -> Result<(), PspReadError> {
    if header.n_records == 0 {
        return Err(PspReadError::BlockHeaderInvariant {
            reason: "n_records must be >= 1 (empty blocks are forbidden)".to_string(),
        });
    }
    if header.n_total_alleles < header.n_records {
        return Err(PspReadError::BlockHeaderInvariant {
            reason: format!(
                "n_total_alleles {} < n_records {} (every record has at least one allele)",
                header.n_total_alleles, header.n_records
            ),
        });
    }
    // Active chain slots strictly ascending, no duplicates.
    for w in header.active_chain_slots.windows(2) {
        if w[0] >= w[1] {
            return Err(PspReadError::BlockHeaderInvariant {
                reason: format!(
                    "active_chain_slots not strictly ascending: {} then {}",
                    w[0], w[1]
                ),
            });
        }
    }
    // Manifest tags strictly ascending, no duplicates.
    for w in header.manifest.windows(2) {
        if w[0].tag >= w[1].tag {
            return Err(PspReadError::BlockHeaderInvariant {
                reason: format!(
                    "manifest tags not strictly ascending: {:#x} then {:#x}",
                    w[0].tag, w[1].tag
                ),
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------- WireScalar round-trips --------------------

    /// Round-trip every fixed-width scalar through encode_le /
    /// decode_le, checking the byte width matches the type's
    /// FIXED_BYTE_WIDTH advertisement.
    #[test]
    fn wire_scalar_round_trip_u8() {
        for v in [0u8, 1, 127, 128, 255] {
            let mut buf = Vec::new();
            v.encode_le(&mut buf);
            assert_eq!(buf.len(), u8::FIXED_BYTE_WIDTH);
            let (decoded, n) = u8::decode_le(&buf).unwrap();
            assert_eq!(decoded, v);
            assert_eq!(n, 1);
        }
    }

    #[test]
    fn wire_scalar_round_trip_u16() {
        for v in [0u16, 1, 256, 65535] {
            let mut buf = Vec::new();
            v.encode_le(&mut buf);
            assert_eq!(buf.len(), 2);
            let (decoded, n) = u16::decode_le(&buf).unwrap();
            assert_eq!(decoded, v);
            assert_eq!(n, 2);
        }
    }

    #[test]
    fn wire_scalar_round_trip_u32() {
        for v in [0u32, 1, u32::MAX] {
            let mut buf = Vec::new();
            v.encode_le(&mut buf);
            let (decoded, n) = u32::decode_le(&buf).unwrap();
            assert_eq!(decoded, v);
            assert_eq!(n, 4);
        }
    }

    #[test]
    fn wire_scalar_round_trip_u64() {
        for v in [0u64, 1, u64::MAX] {
            let mut buf = Vec::new();
            v.encode_le(&mut buf);
            let (decoded, n) = u64::decode_le(&buf).unwrap();
            assert_eq!(decoded, v);
            assert_eq!(n, 8);
        }
    }

    #[test]
    fn wire_scalar_round_trip_i32() {
        for v in [0i32, 1, -1, i32::MIN, i32::MAX] {
            let mut buf = Vec::new();
            v.encode_le(&mut buf);
            let (decoded, _) = i32::decode_le(&buf).unwrap();
            assert_eq!(decoded, v);
        }
    }

    #[test]
    fn wire_scalar_round_trip_i64() {
        for v in [0i64, 1, -1, i64::MIN, i64::MAX] {
            let mut buf = Vec::new();
            v.encode_le(&mut buf);
            let (decoded, _) = i64::decode_le(&buf).unwrap();
            assert_eq!(decoded, v);
        }
    }

    #[test]
    fn wire_scalar_round_trip_f32() {
        for v in [0.0f32, 1.5, -3.14, f32::MIN, f32::MAX] {
            let mut buf = Vec::new();
            v.encode_le(&mut buf);
            let (decoded, _) = f32::decode_le(&buf).unwrap();
            assert_eq!(decoded.to_bits(), v.to_bits());
        }
    }

    #[test]
    fn wire_scalar_round_trip_f64() {
        for v in [0.0f64, 1.5, -3.14, f64::MIN, f64::MAX] {
            let mut buf = Vec::new();
            v.encode_le(&mut buf);
            let (decoded, _) = f64::decode_le(&buf).unwrap();
            assert_eq!(decoded.to_bits(), v.to_bits());
        }
    }

    #[test]
    fn wire_scalar_round_trip_bool() {
        for v in [false, true] {
            let mut buf = Vec::new();
            v.encode_le(&mut buf);
            assert_eq!(buf.len(), 1);
            let (decoded, n) = bool::decode_le(&buf).unwrap();
            assert_eq!(decoded, v);
            assert_eq!(n, 1);
        }
        // Other bytes are an error.
        for bad in [2u8, 3, 0xff] {
            let buf = [bad];
            let err = bool::decode_le(&buf).unwrap_err();
            assert_eq!(err, ScalarDecodeError::InvalidBool(bad));
        }
    }

    /// Truncated buffer in every fixed-width type → `Truncated`.
    #[test]
    fn wire_scalar_truncated() {
        assert_eq!(
            u8::decode_le(&[]).unwrap_err(),
            ScalarDecodeError::Truncated
        );
        assert_eq!(
            u16::decode_le(&[0u8]).unwrap_err(),
            ScalarDecodeError::Truncated
        );
        assert_eq!(
            u32::decode_le(&[0u8; 3]).unwrap_err(),
            ScalarDecodeError::Truncated
        );
        assert_eq!(
            u64::decode_le(&[0u8; 7]).unwrap_err(),
            ScalarDecodeError::Truncated
        );
        assert_eq!(
            f64::decode_le(&[0u8; 4]).unwrap_err(),
            ScalarDecodeError::Truncated
        );
    }

    // ------------------- Scalar column round-trips -----------------

    #[test]
    fn scalar_column_round_trip_u32() {
        let values = vec![1u32, 2, 3, 1_000_000, u32::MAX];
        let mut buf = Vec::new();
        encode_scalar_column(&values, &mut buf);
        assert_eq!(buf.len(), values.len() * 4);
        let decoded: Vec<u32> = decode_scalar_column(&buf, values.len(), "test").unwrap();
        assert_eq!(decoded, values);
    }

    #[test]
    fn scalar_column_round_trip_f64() {
        let values = vec![0.0f64, -3.14, 1e308, -1e-308];
        let mut buf = Vec::new();
        encode_scalar_column(&values, &mut buf);
        let decoded: Vec<f64> = decode_scalar_column(&buf, values.len(), "test").unwrap();
        for (a, b) in decoded.iter().zip(values.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    #[test]
    fn scalar_column_round_trip_empty() {
        let values: Vec<u32> = Vec::new();
        let mut buf = Vec::new();
        encode_scalar_column(&values, &mut buf);
        assert!(buf.is_empty());
        let decoded: Vec<u32> = decode_scalar_column(&buf, 0, "test").unwrap();
        assert!(decoded.is_empty());
    }

    /// Truncated scalar-column payload: encoded N elements, told the
    /// decoder there are N+1.
    #[test]
    fn scalar_column_truncation_under_count() {
        let values = vec![1u32, 2, 3];
        let mut buf = Vec::new();
        encode_scalar_column(&values, &mut buf);
        let err: PspReadError = decode_scalar_column::<u32>(&buf, 4, "obs-count").unwrap_err();
        match err {
            PspReadError::ColumnTruncated {
                column,
                decoded,
                expected,
            } => {
                assert_eq!(column, "obs-count");
                assert_eq!(decoded, 3);
                assert_eq!(expected, 4);
            }
            other => panic!("expected ColumnTruncated, got {other:?}"),
        }
    }

    /// Trailing bytes past the last expected element.
    #[test]
    fn scalar_column_trailing_bytes() {
        let values = vec![1u32, 2, 3];
        let mut buf = Vec::new();
        encode_scalar_column(&values, &mut buf);
        buf.extend_from_slice(&[0xaa, 0xbb]);
        let err: PspReadError = decode_scalar_column::<u32>(&buf, 3, "obs-count").unwrap_err();
        match err {
            PspReadError::ColumnTrailingBytes { trailing, .. } => assert_eq!(trailing, 2),
            other => panic!("expected ColumnTrailingBytes, got {other:?}"),
        }
    }

    /// Invalid bool byte mid-column → ColumnElementDecode pointing at
    /// the offending entry.
    #[test]
    fn scalar_column_invalid_bool_byte() {
        // Encoded: [true, false, ?]  with the third byte being 2.
        let buf = vec![1u8, 0, 2];
        let err: PspReadError = decode_scalar_column::<bool>(&buf, 3, "future-bool").unwrap_err();
        match err {
            PspReadError::ColumnElementDecode {
                column,
                entry,
                source,
            } => {
                assert_eq!(column, "future-bool");
                assert_eq!(entry, 2);
                assert_eq!(source, ScalarDecodeError::InvalidBool(2));
            }
            other => panic!("expected ColumnElementDecode, got {other:?}"),
        }
    }

    // ------------------- Varint column round-trips -----------------

    #[test]
    fn varint_column_round_trip() {
        let values = vec![0u64, 1, 127, 128, 16383, 16384, u64::MAX];
        let mut buf = Vec::new();
        encode_varint_column(&values, &mut buf);
        let decoded = decode_varint_column(&buf, values.len(), "delta-pos").unwrap();
        assert_eq!(decoded, values);
    }

    #[test]
    fn varint_column_empty() {
        let mut buf = Vec::new();
        encode_varint_column(&[], &mut buf);
        let decoded = decode_varint_column(&buf, 0, "delta-pos").unwrap();
        assert!(decoded.is_empty());
    }

    /// Truncation at an element boundary: encoded `[1, 2, 3]` (3
    /// single-byte varints), drop the third → the decoder gets to
    /// entry 2 with an empty remaining buffer and returns
    /// `ColumnTruncated`.
    #[test]
    fn varint_column_truncation_at_boundary() {
        let mut buf = Vec::new();
        encode_varint_column(&[1u64, 2, 3], &mut buf);
        let truncated = &buf[..buf.len() - 1];
        let err = decode_varint_column(truncated, 3, "delta-pos").unwrap_err();
        match err {
            PspReadError::ColumnTruncated {
                decoded, expected, ..
            } => {
                assert_eq!(decoded, 2);
                assert_eq!(expected, 3);
            }
            other => panic!("expected ColumnTruncated, got {other:?}"),
        }
    }

    /// Truncation inside a varint: `200` encodes to two bytes; drop
    /// the second one and the decoder hits `VarintError::Truncated`
    /// at entry 1 → `ColumnElementDecode`.
    #[test]
    fn varint_column_truncation_mid_varint() {
        let mut buf = Vec::new();
        encode_varint_column(&[1u64, 200], &mut buf);
        // 1 byte for value 1, 2 bytes for value 200 = 3 bytes total.
        assert_eq!(buf.len(), 3);
        let truncated = &buf[..buf.len() - 1];
        let err = decode_varint_column(truncated, 2, "delta-pos").unwrap_err();
        match err {
            PspReadError::ColumnElementDecode { entry, source, .. } => {
                assert_eq!(entry, 1);
                assert_eq!(source, ScalarDecodeError::Truncated);
            }
            other => panic!("expected ColumnElementDecode mid-varint, got {other:?}"),
        }
    }

    // ------------------- List column round-trips -------------------

    #[test]
    fn list_column_round_trip_u16() {
        let r0: &[u16] = &[];
        let r1: &[u16] = &[7];
        let r2: &[u16] = &[1, 2, 3, 4, 5];
        let lists: Vec<&[u16]> = vec![r0, r1, r2];
        let mut buf = Vec::new();
        encode_list_column(&lists, &mut buf);
        let decoded: Vec<Vec<u16>> = decode_list_column(&buf, lists.len(), "chain-slots").unwrap();
        assert_eq!(decoded.len(), 3);
        assert!(decoded[0].is_empty());
        assert_eq!(decoded[1], vec![7]);
        assert_eq!(decoded[2], vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn list_column_truncation_in_count() {
        let lists: Vec<&[u16]> = vec![&[1, 2, 3]];
        let mut buf = Vec::new();
        encode_list_column(&lists, &mut buf);
        // Tell decoder there are two entries; the second's count
        // varint is missing.
        let err = decode_list_column::<u16>(&buf, 2, "chain-slots").unwrap_err();
        match err {
            PspReadError::ColumnTruncated {
                decoded, expected, ..
            } => {
                assert_eq!(decoded, 1);
                assert_eq!(expected, 2);
            }
            other => panic!("expected ColumnTruncated, got {other:?}"),
        }
    }

    #[test]
    fn list_column_truncation_in_element() {
        // Claim 3 elements but only encode the count byte (3) and
        // two-and-a-half u16s (4 bytes for the first two, 1 byte for
        // a half).
        let mut buf = Vec::new();
        encode_u64_leb128(3, &mut buf); // count = 3
        buf.extend_from_slice(&7u16.to_le_bytes());
        buf.extend_from_slice(&9u16.to_le_bytes());
        buf.push(0x11); // dangling half of u16
        let err = decode_list_column::<u16>(&buf, 1, "chain-slots").unwrap_err();
        match err {
            PspReadError::ColumnElementDecode { entry, source, .. } => {
                assert_eq!(entry, 0);
                assert_eq!(source, ScalarDecodeError::Truncated);
            }
            other => panic!("expected ColumnElementDecode, got {other:?}"),
        }
    }

    // ------------------- Bytes column round-trips ------------------

    #[test]
    fn bytes_column_round_trip() {
        let v0: &[u8] = b"A";
        let v1: &[u8] = b"ACGT";
        let v2: &[u8] = b"N";
        let values: Vec<&[u8]> = vec![v0, v1, v2];
        let lengths: Vec<u64> = values.iter().map(|v| v.len() as u64).collect();
        let mut bytes_buf = Vec::new();
        encode_bytes_concat(&values, &mut bytes_buf);
        let decoded = decode_bytes_split(&bytes_buf, &lengths, "allele-seq").unwrap();
        assert_eq!(decoded[0], b"A");
        assert_eq!(decoded[1], b"ACGT");
        assert_eq!(decoded[2], b"N");
    }

    #[test]
    fn bytes_column_too_short_for_lengths() {
        let bytes_buf = b"AB".to_vec();
        let lengths = vec![1u64, 2, 3]; // sum=6, buf=2
        let err = decode_bytes_split(&bytes_buf, &lengths, "allele-seq").unwrap_err();
        assert!(matches!(err, PspReadError::ColumnTruncated { .. }));
    }

    #[test]
    fn bytes_column_too_long_for_lengths() {
        let bytes_buf = b"ABCD".to_vec();
        let lengths = vec![1u64, 1]; // sum=2, buf=4
        let err = decode_bytes_split(&bytes_buf, &lengths, "allele-seq").unwrap_err();
        match err {
            PspReadError::ColumnTrailingBytes { trailing, .. } => assert_eq!(trailing, 2),
            other => panic!("expected ColumnTrailingBytes, got {other:?}"),
        }
    }

    // ------------------- zstd round-trip --------------------------

    #[test]
    fn zstd_round_trip() {
        let original: Vec<u8> = (0u32..10_000).flat_map(|i| i.to_le_bytes()).collect();
        let compressed = zstd_compress(&original).unwrap();
        // Compression actually reduces the input on this regular
        // pattern.
        assert!(compressed.len() < original.len());
        let decompressed = zstd_decompress(&compressed).unwrap();
        assert_eq!(decompressed, original);
    }

    #[test]
    fn zstd_detects_corruption() {
        let original: Vec<u8> = (0u32..1000).flat_map(|i| i.to_le_bytes()).collect();
        let mut compressed = zstd_compress(&original).unwrap();
        // The content checksum is the last 4 bytes of the frame
        // when `include_checksum(true)` is set. Flip the very last
        // byte — that is guaranteed to be checksum bytes, so
        // decompression must fail.
        let last = compressed.len() - 1;
        compressed[last] ^= 0x01;
        let err = zstd_decompress(&compressed)
            .expect_err("zstd frame checksum should reject this byte-flip");
        // Sanity-check: it's a decompression error, not e.g. EOF.
        assert!(
            format!("{err}").to_lowercase().contains("checksum")
                || err.kind() != io::ErrorKind::UnexpectedEof
        );
    }

    // ------------------- Block header round-trip ------------------

    fn sample_block_header() -> BlockHeader {
        BlockHeader {
            chrom_id: 0,
            first_pos: 1,
            n_records: 100,
            n_total_alleles: 130,
            active_chain_slots: vec![3, 7, 42],
            manifest: vec![
                ColumnManifestEntry {
                    tag: 0x01,
                    compressed_len: 50,
                    uncompressed_len: 100,
                },
                ColumnManifestEntry {
                    tag: 0x02,
                    compressed_len: 60,
                    uncompressed_len: 100,
                },
                ColumnManifestEntry {
                    tag: 0x11,
                    compressed_len: 800,
                    uncompressed_len: 1040,
                },
            ],
        }
    }

    #[test]
    fn block_header_round_trip() {
        let header = sample_block_header();
        let mut buf = Vec::new();
        encode_block_header(&header, &mut buf).unwrap();
        let (decoded, consumed) = decode_block_header(&buf).unwrap();
        assert_eq!(consumed, buf.len());
        assert_eq!(decoded, header);
    }

    #[test]
    fn block_header_round_trip_empty_active_slots() {
        let mut header = sample_block_header();
        header.active_chain_slots.clear();
        let mut buf = Vec::new();
        encode_block_header(&header, &mut buf).unwrap();
        let (decoded, _) = decode_block_header(&buf).unwrap();
        assert_eq!(decoded, header);
    }

    #[test]
    fn block_header_rejects_zero_records() {
        let mut header = sample_block_header();
        header.n_records = 0;
        let mut buf = Vec::new();
        let err = encode_block_header(&header, &mut buf).unwrap_err();
        assert!(matches!(err, PspReadError::BlockHeaderInvariant { .. }));
    }

    #[test]
    fn block_header_rejects_n_total_alleles_less_than_n_records() {
        let mut header = sample_block_header();
        header.n_records = 100;
        header.n_total_alleles = 99;
        let mut buf = Vec::new();
        let err = encode_block_header(&header, &mut buf).unwrap_err();
        assert!(matches!(err, PspReadError::BlockHeaderInvariant { .. }));
    }

    #[test]
    fn block_header_rejects_unsorted_active_slots() {
        let mut header = sample_block_header();
        header.active_chain_slots = vec![7, 3, 42];
        let mut buf = Vec::new();
        let err = encode_block_header(&header, &mut buf).unwrap_err();
        assert!(matches!(err, PspReadError::BlockHeaderInvariant { .. }));
    }

    #[test]
    fn block_header_rejects_duplicate_active_slots() {
        let mut header = sample_block_header();
        header.active_chain_slots = vec![3, 3, 7];
        let mut buf = Vec::new();
        let err = encode_block_header(&header, &mut buf).unwrap_err();
        assert!(matches!(err, PspReadError::BlockHeaderInvariant { .. }));
    }

    #[test]
    fn block_header_rejects_unsorted_manifest() {
        let mut header = sample_block_header();
        header.manifest[1].tag = 0x00;
        let mut buf = Vec::new();
        let err = encode_block_header(&header, &mut buf).unwrap_err();
        assert!(matches!(err, PspReadError::BlockHeaderInvariant { .. }));
    }

    #[test]
    fn block_header_rejects_duplicate_manifest_tags() {
        let mut header = sample_block_header();
        header.manifest[1].tag = 0x01; // same as manifest[0]
        let mut buf = Vec::new();
        let err = encode_block_header(&header, &mut buf).unwrap_err();
        assert!(matches!(err, PspReadError::BlockHeaderInvariant { .. }));
    }

    #[test]
    fn block_header_decode_rejects_truncated_buffer() {
        let header = sample_block_header();
        let mut buf = Vec::new();
        encode_block_header(&header, &mut buf).unwrap();
        // Truncate just past the n_records field.
        let truncated = &buf[..5.min(buf.len())];
        let err = decode_block_header(truncated).unwrap_err();
        assert!(matches!(err, PspReadError::BlockHeaderField { .. }));
    }
}
