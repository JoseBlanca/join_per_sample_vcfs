//! The fixed 32-byte file trailer.
//!
//! The trailer sits at the very end of every `.psp` file. It carries:
//!
//! - `index_offset` (`u64`): absolute file offset of the block
//!   index's first byte.
//! - `index_byte_length` (`u64`): byte length of the block index.
//! - `n_blocks` (`u64`): total number of blocks in the file (must
//!   equal the count of entries the reader decodes from the index).
//! - `index_checksum` (`u32`): low 32 bits of XXH3-64 over the index
//!   bytes — closes the one file region not covered by a zstd frame
//!   checksum (spec §"File trailer" / Q-PL11).
//! - `trailer_magic` (4 B): `PSPE`, placed last so an EOF-aligned
//!   `pread` of the final four bytes can reject a non-`.psp` or
//!   truncated file before the rest of the trailer is touched.
//!
//! Wire layout is fixed-width little-endian per the spec's body
//! encoding conventions; this module only handles the byte
//! marshalling, not the cross-file arithmetic (which lives in the
//! reader once it knows the file size).

use super::errors::PspReadError;

/// Fixed byte length of the trailer.
pub const TRAILER_BYTES: usize = 32;

/// Trailing magic bytes — `PSPE` (ASCII). Different from the head
/// magic so truncation that copies the head pattern would not pass
/// the tail check.
pub const TRAILER_MAGIC: [u8; 4] = *b"PSPE";

/// In-memory form of the file trailer.
///
/// Field order matches the wire layout; encoding writes each field's
/// little-endian bytes in this order back to back.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Trailer {
    pub index_offset: u64,
    pub index_byte_length: u64,
    pub n_blocks: u64,
    pub index_checksum: u32,
}

/// Encode the trailer to its 32-byte wire form.
pub fn encode_trailer(t: &Trailer) -> [u8; TRAILER_BYTES] {
    let mut bytes = [0u8; TRAILER_BYTES];
    bytes[0..8].copy_from_slice(&t.index_offset.to_le_bytes());
    bytes[8..16].copy_from_slice(&t.index_byte_length.to_le_bytes());
    bytes[16..24].copy_from_slice(&t.n_blocks.to_le_bytes());
    bytes[24..28].copy_from_slice(&t.index_checksum.to_le_bytes());
    bytes[28..32].copy_from_slice(&TRAILER_MAGIC);
    bytes
}

/// Decode the trailer from its 32-byte wire form. Verifies the
/// magic; rejects with [`PspReadError::BadTrailerMagic`] on
/// mismatch. Self-consistency arithmetic across the file (such as
/// `index_offset + index_byte_length == trailer_start_offset`) is
/// the reader's responsibility, not this function's.
pub fn decode_trailer(bytes: &[u8; TRAILER_BYTES]) -> Result<Trailer, PspReadError> {
    let mut magic = [0u8; 4];
    magic.copy_from_slice(&bytes[28..32]);
    if magic != TRAILER_MAGIC {
        return Err(PspReadError::BadTrailerMagic {
            got: magic,
            expected: TRAILER_MAGIC,
        });
    }
    Ok(Trailer {
        index_offset: u64::from_le_bytes(bytes[0..8].try_into().unwrap()),
        index_byte_length: u64::from_le_bytes(bytes[8..16].try_into().unwrap()),
        n_blocks: u64::from_le_bytes(bytes[16..24].try_into().unwrap()),
        index_checksum: u32::from_le_bytes(bytes[24..28].try_into().unwrap()),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// (Plan group V — V3.) Build a representative trailer, encode,
    /// decode, assert equality. Covers byte-level marshalling and
    /// magic placement.
    #[test]
    fn round_trip() {
        let t = Trailer {
            index_offset: 1234567,
            index_byte_length: 4096,
            n_blocks: 17,
            index_checksum: 0xDEADBEEF,
        };
        let bytes = encode_trailer(&t);
        let decoded = decode_trailer(&bytes).expect("decode should succeed");
        assert_eq!(decoded, t);
    }

    /// Encoded magic sits at the last four bytes — the EOF-aligned
    /// `pread` that readers do can identify a `.psp` from just those
    /// four bytes.
    #[test]
    fn magic_is_at_the_tail() {
        let t = Trailer {
            index_offset: 0,
            index_byte_length: 0,
            n_blocks: 0,
            index_checksum: 0,
        };
        let bytes = encode_trailer(&t);
        assert_eq!(&bytes[28..32], &TRAILER_MAGIC);
    }

    /// Wrong magic bytes fail decode with a clear error naming the
    /// observed and expected magics. (Plan group V — V3 negative.)
    #[test]
    fn decode_rejects_bad_magic() {
        let mut bytes = encode_trailer(&Trailer {
            index_offset: 0,
            index_byte_length: 0,
            n_blocks: 0,
            index_checksum: 0,
        });
        // Flip one byte of the magic.
        bytes[28] = b'X';
        let err = decode_trailer(&bytes).unwrap_err();
        match err {
            PspReadError::BadTrailerMagic { got, expected } => {
                assert_eq!(got, [b'X', b'S', b'P', b'E']);
                assert_eq!(expected, TRAILER_MAGIC);
            }
            other => panic!("expected BadTrailerMagic, got {other:?}"),
        }
    }

    /// Field order on the wire matches the documented layout: every
    /// field at its expected offset. This is the load-bearing
    /// property — a future refactor that accidentally swaps field
    /// positions would break every existing `.psp`.
    #[test]
    fn wire_layout_is_stable() {
        let t = Trailer {
            index_offset: 0x0102030405060708,
            index_byte_length: 0x1112131415161718,
            n_blocks: 0x2122232425262728,
            index_checksum: 0x31323334,
        };
        let bytes = encode_trailer(&t);
        assert_eq!(&bytes[0..8], &0x0102030405060708u64.to_le_bytes());
        assert_eq!(&bytes[8..16], &0x1112131415161718u64.to_le_bytes());
        assert_eq!(&bytes[16..24], &0x2122232425262728u64.to_le_bytes());
        assert_eq!(&bytes[24..28], &0x31323334u32.to_le_bytes());
        assert_eq!(&bytes[28..32], &TRAILER_MAGIC);
    }

    /// All-zero values are a legal trailer (empty file: header +
    /// zero-entry index + trailer; the index's XXH3-64 of the empty
    /// byte sequence happens to be a specific non-zero constant, but
    /// the trailer struct itself accepts any value here).
    #[test]
    fn all_zero_fields_round_trip() {
        let t = Trailer {
            index_offset: 0,
            index_byte_length: 0,
            n_blocks: 0,
            index_checksum: 0,
        };
        let bytes = encode_trailer(&t);
        let decoded = decode_trailer(&bytes).unwrap();
        assert_eq!(decoded, t);
    }
}
