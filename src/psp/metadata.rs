//! The optional, schema-agnostic **metadata section**.
//!
//! A `.psp` may carry one metadata section, written at `finish()`
//! between the block index and the 32-byte trailer:
//!
//! ```text
//! +---------+--------+-------+----------+---------+
//! | header  | blocks | index | metadata | trailer |
//! +---------+--------+-------+----------+---------+
//!                            ^          ^
//!                            index_end  trailer_start
//! ```
//!
//! The container core treats the section as **opaque compressed
//! bytes**: a single zstd frame (content checksum on, like every block
//! payload). The decompressed bytes are a kind-defined document — the
//! `snp` schema puts a TOML document of per-sample summary statistics
//! there (architecture
//! `doc/devel/architecture/hidden_paralog_psp_integration.md`), but
//! this module neither produces nor parses that document.
//!
//! **Location, not a trailer pointer.** The trailer's 32-byte layout is
//! unchanged; the reader recovers the section as the bytes between
//! `index_offset + index_byte_length` and the trailer start. A
//! zero-length gap means *no section*, so a producer that writes none
//! (and every pre-existing `.psp`) stays byte-identical and valid.
//!
//! The section region must be **exactly one frame** — the reader
//! rejects any trailing bytes past the frame
//! ([`PspReadError::MetadataTrailingBytes`]), since the trailer's index
//! checksum stops at `index_end` and the frame checksum covers only the
//! frame payload, leaving the post-frame remainder otherwise unchecked.

use std::io::Read;

use super::block::zstd_compress;
use super::errors::PspReadError;

/// Upper bound on a metadata section's **decompressed** size. The
/// legitimate payload (per-sample coverage histogram + het counts) is
/// kilobytes; this cap exists only so a hostile or corrupt file cannot
/// drive an unbounded allocation through a zip-bomb frame. 64 MiB is
/// far above any real summary and far below an OOM risk.
pub(crate) const MAX_METADATA_DECOMPRESSED_BYTES: usize = 64 * 1024 * 1024;

/// Compress an opaque metadata payload into the on-disk section frame
/// (one zstd frame, content-checksummed). Returns the frame bytes the
/// writer appends after the block index.
pub(crate) fn compress_metadata(raw: &[u8]) -> std::io::Result<Vec<u8>> {
    zstd_compress(raw)
}

/// Decompress a metadata section frame back into its opaque payload
/// bytes. Rejects, before materialising the full output, a frame that
/// decompresses past [`MAX_METADATA_DECOMPRESSED_BYTES`] (zip-bomb
/// guard) or a region that carries bytes past the single zstd frame
/// (`MetadataTrailingBytes`).
pub(crate) fn decompress_metadata(frame: &[u8]) -> Result<Vec<u8>, PspReadError> {
    decompress_metadata_capped(frame, MAX_METADATA_DECOMPRESSED_BYTES)
}

/// [`decompress_metadata`] with an explicit decompressed-size cap. The
/// public wrapper pins the cap to [`MAX_METADATA_DECOMPRESSED_BYTES`];
/// the parameter exists so tests can drive the zip-bomb guard with a
/// small cap instead of a 64 MiB fixture.
fn decompress_metadata_capped(frame: &[u8], cap: usize) -> Result<Vec<u8>, PspReadError> {
    let zstd_err = |source| PspReadError::Zstd {
        context: "metadata section".to_string(),
        source,
    };

    // The section region must be exactly one frame. A streaming decoder
    // stops at the end of the first frame and would silently ignore any
    // trailing bytes; the trailer's index checksum stops at the frame's
    // start and the frame checksum covers only its own payload, so the
    // post-frame remainder is otherwise unchecked. Compare the first
    // frame's compressed length against the whole region. (This also
    // rejects a truncated frame: `find_frame_compressed_size` needs a
    // complete frame and errors otherwise.)
    let frame_len = zstd::zstd_safe::find_frame_compressed_size(frame).map_err(|code| {
        zstd_err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            zstd::zstd_safe::get_error_name(code),
        ))
    })?;
    if frame_len != frame.len() {
        return Err(PspReadError::MetadataTrailingBytes {
            trailing: frame.len() - frame_len,
        });
    }

    // `zstd_decompress` (single-shot `decode_all`) would allocate the
    // whole decompressed output before we could check its size, so a
    // bomb frame would be materialised first. Stream instead and stop
    // one byte past the cap. `saturating_add` keeps the `cap == usize::MAX`
    // edge (never reachable in production) from wrapping the limit to 0.
    let mut decoder = zstd::stream::read::Decoder::new(frame).map_err(zstd_err)?;
    let mut out = Vec::new();
    let read = decoder
        .by_ref()
        .take((cap as u64).saturating_add(1))
        .read_to_end(&mut out)
        .map_err(zstd_err)?;
    if read > cap {
        return Err(PspReadError::MetadataSectionTooLarge { cap });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Round-trip: an arbitrary payload compresses and decompresses
    /// back to the same bytes.
    #[test]
    fn compress_decompress_round_trips() {
        let payload = b"[coverage-gc]\nwindow-bp = 500\n".repeat(40);
        let frame = compress_metadata(&payload).expect("compress");
        let back = decompress_metadata(&frame).expect("decompress");
        assert_eq!(back, payload);
    }

    /// An empty payload still round-trips (the SNP payload is never
    /// empty, but the codec must be total).
    #[test]
    fn empty_payload_round_trips() {
        let frame = compress_metadata(b"").expect("compress");
        let back = decompress_metadata(&frame).expect("decompress");
        assert!(back.is_empty());
    }

    /// `compress_metadata` is deterministic: the same payload always
    /// produces byte-identical frame bytes. Pins the on-disk stability
    /// contract (a future change such as multithreaded zstd that
    /// reordered output would be caught here, not silently shipped).
    #[test]
    fn compress_metadata_is_deterministic() {
        let payload = b"[coverage-gc]\nwindow-bp = 500\n".repeat(7);
        assert_eq!(
            compress_metadata(&payload).expect("compress a"),
            compress_metadata(&payload).expect("compress b"),
        );
    }

    /// A frame that decompresses just under a cap is accepted; one byte
    /// over is rejected with `MetadataSectionTooLarge`. Driven through
    /// `decompress_metadata_capped` with a small cap so the test does
    /// not allocate the 64 MiB production cap.
    #[test]
    fn rejects_payload_over_decompressed_cap() {
        const TEST_CAP: usize = 4096;

        let over = vec![0u8; TEST_CAP + 1];
        let frame = compress_metadata(&over).expect("compress");
        let err = decompress_metadata_capped(&frame, TEST_CAP).expect_err("over-cap must fail");
        assert!(
            matches!(err, PspReadError::MetadataSectionTooLarge { cap } if cap == TEST_CAP),
            "expected MetadataSectionTooLarge, got {err:?}"
        );

        let at_cap = vec![0u8; TEST_CAP];
        let frame = compress_metadata(&at_cap).expect("compress");
        let back = decompress_metadata_capped(&frame, TEST_CAP).expect("at-cap must pass");
        assert_eq!(back.len(), TEST_CAP);
    }

    /// A corrupt frame (mid-frame byte flip) surfaces as `Zstd`, not a
    /// panic.
    #[test]
    fn corrupt_frame_errors_cleanly() {
        let frame = compress_metadata(b"hello world").expect("compress");
        let mut corrupt = frame.clone();
        let mid = corrupt.len() / 2;
        corrupt[mid] ^= 0xFF;
        let err = decompress_metadata(&corrupt).expect_err("corruption must fail");
        assert!(
            matches!(err, PspReadError::Zstd { .. }),
            "expected Zstd, got {err:?}"
        );
    }

    /// A truncated frame (cut short before its end) surfaces as `Zstd`,
    /// not a partial-output success or a hang — the frame-completeness
    /// check rejects it.
    #[test]
    fn truncated_frame_errors_cleanly() {
        let frame = compress_metadata(b"hello world, this is a payload").expect("compress");
        let truncated = &frame[..frame.len() - 4];
        let err = decompress_metadata(truncated).expect_err("truncation must fail");
        assert!(
            matches!(err, PspReadError::Zstd { .. }),
            "expected Zstd, got {err:?}"
        );
    }

    /// Bytes that are not a zstd frame at all surface as `Zstd`.
    #[test]
    fn non_frame_bytes_error_cleanly() {
        let err = decompress_metadata(b"not-a-zstd-frame-at-all").expect_err("non-frame must fail");
        assert!(
            matches!(err, PspReadError::Zstd { .. }),
            "expected Zstd, got {err:?}"
        );
    }

    /// A valid frame followed by trailing bytes is rejected with
    /// `MetadataTrailingBytes` — the region must be exactly one frame,
    /// so the post-frame remainder cannot pass unchecked.
    #[test]
    fn rejects_trailing_bytes_after_frame() {
        let frame = compress_metadata(b"summary payload").expect("compress");
        let mut padded = frame.clone();
        padded.extend_from_slice(b"junk");
        let err = decompress_metadata(&padded).expect_err("trailing bytes must fail");
        assert!(
            matches!(err, PspReadError::MetadataTrailingBytes { trailing } if trailing == 4),
            "expected MetadataTrailingBytes {{ trailing: 4 }}, got {err:?}"
        );
    }
}
