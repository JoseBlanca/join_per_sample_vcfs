use std::fs;
use std::io::Write;
use std::path::PathBuf;

use tempfile::NamedTempFile;

use join_per_sample_vcfs::utils_magic::{are_gzipped_magic_bytes, file_is_gzipped};

#[test]
fn are_gzipped_magic_bytes_detects_gzip_header() {
    let bytes = [0x1f, 0x8b, 0x00, 0x00];
    assert!(are_gzipped_magic_bytes(&bytes));
}

#[test]
fn are_gzipped_magic_bytes_detects_non_gzip_header() {
    let bytes = [0x00, 0x01, 0x02, 0x03];
    assert!(!are_gzipped_magic_bytes(&bytes));
}

#[test]
fn are_gzipped_magic_bytes_returns_false_for_insufficient_bytes() {
    let empty: [u8; 0] = [];
    let one_byte = [0x1f];

    assert!(!are_gzipped_magic_bytes(&empty));
    assert!(!are_gzipped_magic_bytes(&one_byte));
}

#[test]
fn file_is_gzipped_returns_true_for_file_with_gzip_magic_bytes() {
    let mut tmp = NamedTempFile::new().unwrap();

    tmp.write_all(&[0x1f, 0x8b, 0x00, 0x00]).unwrap();
    tmp.flush().unwrap();

    let path = tmp.path();
    let result = file_is_gzipped(path).unwrap();

    assert!(result);
}

#[test]
fn file_is_gzipped_returns_false_for_plain_file() {
    let mut tmp = NamedTempFile::new().unwrap();

    tmp.write_all(&[0x00, 0x01, 0x02, 0x03]).unwrap();
    tmp.flush().unwrap();

    let path = tmp.path();
    let result = file_is_gzipped(path).unwrap();

    assert!(!result);
}

#[test]
fn file_is_gzipped_returns_false_for_short_file() {
    let mut tmp = NamedTempFile::new().unwrap();

    // Only 1 byte; too short to match gzip magic.
    tmp.write_all(&[0x1f]).unwrap();
    tmp.flush().unwrap();

    let path = tmp.path();
    let result = file_is_gzipped(path).unwrap();

    assert!(!result);
}

#[test]
fn file_is_gzipped_propagates_io_error_for_missing_file() {
    // Create a temp file, get its path, then delete it so the path is invalid.
    let tmp = NamedTempFile::new().unwrap();
    let path: PathBuf = tmp.path().to_path_buf();
    drop(tmp);

    let err = file_is_gzipped(&path).unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::NotFound);
}

#[test]
fn file_is_gzipped_works_with_paths_from_fs_module() {
    // Slightly more realistic test using fs::write and PathBuf.
    let tmp_dir = tempfile::tempdir().unwrap();
    let file_path = tmp_dir.path().join("test.gz");

    fs::write(&file_path, &[0x1f, 0x8b, 0x08, 0x00]).unwrap();
    let result = file_is_gzipped(&file_path).unwrap();

    assert!(result);
}
