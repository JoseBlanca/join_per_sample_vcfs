// Mi5: this module is `pub(crate)`; its decoder primitives plus
// the signed-varint codecs are reader-side surface awaiting the
// not-yet-built `PspReader`. They're exercised by tests but are
// not yet reached from any production code path.
#![allow(dead_code)]
//! LEB128 / zig-zag-LEB128 encoders and decoders.
//!
//! The `.psp` body uses unsigned LEB128 for every variable-length
//! integer and zig-zag-LEB128 ("svarint") for the handful of signed
//! deltas. Both encodings are bounded at 10 bytes — the cap for
//! `u64` / `i64`. See the spec's §"Encoding conventions / Body
//! (binary)" rows for the canonical definitions.
//!
//! This module is the only place in the `.psp` codebase that
//! understands those byte sequences. The writer pipes integers
//! through [`encode_u64_leb128`] / [`encode_i64_svarint`] into per-
//! column buffers; the reader pipes the same buffers back through
//! [`decode_u64_leb128`] / [`decode_i64_svarint`].
//!
//! Both decoders return `(value, bytes_consumed)` on success and a
//! [`VarintError`] on a truncated or over-long encoding. The caller
//! adds context (block index, column, record index) when wrapping
//! the error into the relevant [`super::errors::PspReadError`]
//! variant — `IndexEntryDecode`, `BlockHeaderField`, or
//! `ColumnElementDecode`, depending on which file region was being
//! decoded.

use super::errors::VarintError;

/// Maximum number of LEB128 bytes a `u64` can occupy: `ceil(64 / 7) =
/// 10`. The spec pins this cap (§"Encoding conventions / Body
/// (binary)" → `varint` row); a longer encoding cannot represent a
/// valid `u64` and is treated as corruption.
pub const MAX_VARINT_BYTES: usize = 10;

/// Encode `value` as unsigned LEB128 and append the bytes to `out`.
/// Always produces between 1 and 10 bytes. Never fails.
///
/// Split into an inlined fast path (`value < 0x80`, one
/// `Vec::push`) and a `#[cold] #[inline(never)]` multi-byte path,
/// per L4. On the SNP-typical writer workload the vast majority of
/// values (`delta_pos == 1`, `n_alleles == 1`, `allele_seq_len == 1`,
/// per-list-count `== 0`) take the fast path; LLVM can lay the cold
/// body out-of-line so the inlined call site stays tight.
#[inline]
pub fn encode_u64_leb128(value: u64, out: &mut Vec<u8>) {
    if value < 0x80 {
        out.push(value as u8);
    } else {
        encode_u64_leb128_cold(value, out);
    }
}

#[cold]
#[inline(never)]
fn encode_u64_leb128_cold(mut value: u64, out: &mut Vec<u8>) {
    while value >= 0x80 {
        out.push((value as u8) | 0x80);
        value >>= 7;
    }
    out.push(value as u8);
}

/// Decode a single unsigned LEB128 integer from the front of `bytes`.
///
/// Returns the decoded value and the number of bytes consumed.
///
/// Split into an inlined fast path (`bytes[0] < 0x80`, one load + one
/// compare) and a `#[cold] #[inline(never)]` multi-byte path, mirror
/// of [`encode_u64_leb128`]'s L4 split. On the SNP-typical reader
/// workload almost every per-record varint is single-byte
/// (`delta_pos == 1`, `n_alleles == 1`, list-counts == 0), so the
/// fast path covers the dominant case and LLVM keeps the cold body
/// out of the call site's icache footprint. L5 in
/// `ia/reviews/perf_psp_reader_2026-05-13.md`.
///
/// Errors:
/// - [`VarintError::Truncated`] if the buffer ends before a byte
///   without the continuation bit is found.
/// - [`VarintError::Overflow`] if more than [`MAX_VARINT_BYTES`]
///   continuation bytes are consumed.
#[inline]
pub fn decode_u64_leb128(bytes: &[u8]) -> Result<(u64, usize), VarintError> {
    match bytes.first() {
        Some(&b) if b < 0x80 => Ok((b as u64, 1)),
        Some(_) => decode_u64_leb128_cold(bytes),
        None => Err(VarintError::Truncated),
    }
}

#[cold]
#[inline(never)]
fn decode_u64_leb128_cold(bytes: &[u8]) -> Result<(u64, usize), VarintError> {
    let mut value: u64 = 0;
    let mut shift: u32 = 0;
    for (i, &b) in bytes.iter().enumerate().take(MAX_VARINT_BYTES) {
        let data = u64::from(b & 0x7f);
        value |= data << shift;
        if b & 0x80 == 0 {
            return Ok((value, i + 1));
        }
        shift += 7;
    }
    // If we get here, either we consumed exactly MAX_VARINT_BYTES with
    // continuation still set (overflow), or the buffer was shorter
    // than that and ran out mid-varint (truncated).
    if bytes.len() >= MAX_VARINT_BYTES {
        Err(VarintError::Overflow)
    } else {
        Err(VarintError::Truncated)
    }
}

/// Encode an `i64` using zig-zag-LEB128 ("svarint") and append the
/// bytes to `out`. Zig-zag maps small-magnitude signed integers to
/// small-magnitude unsigned integers so the LEB128 byte count is
/// proportional to the *magnitude* of the input, not the
/// sign-extended bit pattern.
pub fn encode_i64_svarint(value: i64, out: &mut Vec<u8>) {
    encode_u64_leb128(zigzag_encode(value), out);
}

/// Decode a single `i64` zig-zag-LEB128 value from the front of
/// `bytes`. Returns the decoded value and the number of bytes
/// consumed; errors as for [`decode_u64_leb128`].
pub fn decode_i64_svarint(bytes: &[u8]) -> Result<(i64, usize), VarintError> {
    let (z, n) = decode_u64_leb128(bytes)?;
    Ok((zigzag_decode(z), n))
}

/// Zig-zag forward: signed → unsigned, small-magnitude-preserving.
/// `n => (n << 1) ^ (n >> 63)` with wrapping shifts so `i64::MIN`
/// and `i64::MAX` do not panic in debug builds.
fn zigzag_encode(n: i64) -> u64 {
    // `n.wrapping_shl(1)` is the signed-magnitude two's-complement
    // shift; `n.wrapping_shr(63)` is the arithmetic right-shift,
    // which fills with the sign bit. XOR turns "negative numbers
    // look big as u64" into "negative numbers map to odd
    // small u64s alongside positive evens".
    (n.wrapping_shl(1) ^ n.wrapping_shr(63)) as u64
}

/// Inverse of [`zigzag_encode`].
fn zigzag_decode(z: u64) -> i64 {
    // `(z >> 1)` recovers magnitude / 2; `-(z & 1)` is 0 for even
    // (originally non-negative) and -1 for odd (originally
    // negative), which when XORed with the magnitude reconstructs
    // the signed value.
    ((z >> 1) as i64) ^ -((z & 1) as i64)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `encode → decode` round-trip for every load-bearing `u64`
    /// value (the LEB128 width boundaries plus the limits).
    /// (Plan group V — V1.)
    #[test]
    fn u64_round_trip() {
        let cases: &[(u64, usize)] = &[
            (0, 1),
            (1, 1),
            (127, 1),
            (128, 2),
            (16383, 2),
            (16384, 3),
            (u32::MAX as u64, 5),
            (u64::MAX, 10),
        ];
        for &(value, expected_len) in cases {
            let mut buf = Vec::new();
            encode_u64_leb128(value, &mut buf);
            assert_eq!(
                buf.len(),
                expected_len,
                "value {value} should encode to {expected_len} bytes, got {}",
                buf.len()
            );
            let (decoded, consumed) = decode_u64_leb128(&buf).expect("decode should succeed");
            assert_eq!(decoded, value, "round-trip value mismatch");
            assert_eq!(consumed, buf.len(), "consumed != encoded length");
        }
    }

    /// Trailing bytes past the varint are not consumed. The decoder
    /// reads exactly the bytes belonging to the value, leaving the
    /// remainder for the next decode call.
    #[test]
    fn u64_decode_leaves_trailing_bytes_alone() {
        let mut buf = Vec::new();
        encode_u64_leb128(300, &mut buf); // multi-byte
        buf.extend_from_slice(b"trailing");
        let (value, consumed) = decode_u64_leb128(&buf).unwrap();
        assert_eq!(value, 300);
        assert!(
            consumed < buf.len(),
            "should not consume the trailing bytes"
        );
        assert_eq!(&buf[consumed..], b"trailing");
    }

    /// A varint with the continuation bit set on the final available
    /// byte must report [`VarintError::Truncated`], not silently
    /// return the partial value.
    #[test]
    fn u64_decode_rejects_truncated() {
        // 0x80 = continuation bit set, no data, no further byte.
        let buf = [0x80u8];
        assert_eq!(decode_u64_leb128(&buf), Err(VarintError::Truncated));

        // Two continuation bytes, no terminator.
        let buf = [0x80u8, 0x80];
        assert_eq!(decode_u64_leb128(&buf), Err(VarintError::Truncated));

        // Empty buffer.
        assert_eq!(decode_u64_leb128(&[]), Err(VarintError::Truncated));
    }

    /// More than 10 continuation bytes is [`VarintError::Overflow`],
    /// not [`VarintError::Truncated`] — the buffer is long enough
    /// that "truncation" doesn't fit the failure shape, the encoding
    /// itself is illegal.
    #[test]
    fn u64_decode_rejects_overflow() {
        // 10 continuation bytes followed by a non-terminator — the
        // decoder bails at the 10-byte cap.
        let buf = [0x80u8; 11];
        assert_eq!(decode_u64_leb128(&buf), Err(VarintError::Overflow));
    }

    /// Reasonable representatives of `i64` round-tripping through
    /// svarint. (Plan group V — V2.)
    #[test]
    fn i64_round_trip() {
        let cases: &[i64] = &[
            0,
            1,
            -1,
            63,
            -63,
            64,
            -64,
            i32::MIN as i64,
            i32::MAX as i64,
            i64::MIN,
            i64::MAX,
        ];
        for &value in cases {
            let mut buf = Vec::new();
            encode_i64_svarint(value, &mut buf);
            let (decoded, consumed) = decode_i64_svarint(&buf).expect("decode should succeed");
            assert_eq!(decoded, value, "round-trip value mismatch for {value}");
            assert_eq!(consumed, buf.len());
        }
    }

    /// Zig-zag encoding maps small-magnitude signed values to small
    /// unsigned values, alternating odd/even for negative/positive.
    /// This is the property the encoding exists for; pin it directly.
    #[test]
    fn zigzag_keeps_small_magnitudes_small() {
        assert_eq!(zigzag_encode(0), 0);
        assert_eq!(zigzag_encode(-1), 1);
        assert_eq!(zigzag_encode(1), 2);
        assert_eq!(zigzag_encode(-2), 3);
        assert_eq!(zigzag_encode(2), 4);
        // Decode side mirrors:
        assert_eq!(zigzag_decode(0), 0);
        assert_eq!(zigzag_decode(1), -1);
        assert_eq!(zigzag_decode(2), 1);
        assert_eq!(zigzag_decode(3), -2);
        assert_eq!(zigzag_decode(4), 2);
    }

    /// `i64::MIN` and `i64::MAX` do not panic in debug builds when
    /// zigzag-encoded — the implementation uses `wrapping_shl` /
    /// `wrapping_shr` for exactly this reason.
    #[test]
    fn zigzag_handles_signed_extremes() {
        let _ = zigzag_encode(i64::MIN);
        let _ = zigzag_encode(i64::MAX);
        assert_eq!(zigzag_decode(zigzag_encode(i64::MIN)), i64::MIN);
        assert_eq!(zigzag_decode(zigzag_encode(i64::MAX)), i64::MAX);
    }

    // ------------------- Property-based round-trips (M13) ----------

    use proptest::prelude::*;

    proptest! {
        /// Round-trip property: `decode(encode(v))` recovers every
        /// `u64`, and the byte count consumed equals the byte count
        /// produced. Catches off-by-one cursor / shift bugs that
        /// example-based tests miss.
        #[test]
        fn proptest_u64_round_trip(v in any::<u64>()) {
            let mut buf = Vec::new();
            encode_u64_leb128(v, &mut buf);
            let (decoded, consumed) = decode_u64_leb128(&buf).expect("decode succeeds");
            prop_assert_eq!(decoded, v);
            prop_assert_eq!(consumed, buf.len());
        }

        #[test]
        fn proptest_i64_round_trip(v in any::<i64>()) {
            let mut buf = Vec::new();
            encode_i64_svarint(v, &mut buf);
            let (decoded, consumed) = decode_i64_svarint(&buf).expect("decode succeeds");
            prop_assert_eq!(decoded, v);
            prop_assert_eq!(consumed, buf.len());
        }
    }
}
