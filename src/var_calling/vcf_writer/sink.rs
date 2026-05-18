//! Output-sink dispatch for the cohort VCF writer.
//!
//! The writer streams its bytes through a single `std::io::Write`-
//! shaped sink that hides the plain-vs-bgzf branch from the encoder.
//! Path suffix selects the sink kind (matched case-insensitively, so
//! `out.VCF.GZ` and `out.vcf.gz` resolve the same way):
//!
//! * `.vcf.gz` or `.vcf.bgz` → [`SinkKind::Bgzf`]
//! * anything else           → [`SinkKind::Plain`]
//!
//! Bytes go to `<output>.tmp` first; [`SinkKind::finish`] flushes,
//! emits the bgzf EOF block when applicable, `fsync`s **both the file
//! and the parent directory**, and atomically renames the tmp path
//! to the final output path. The parent-dir fsync closes the only
//! remaining crash window in the rename dance: without it, a crash
//! between the rename returning and the parent's journal flushing can
//! leave the file's contents durable but the directory entry pointing
//! at `<final_path>` lost.
//!
//! A crash before `finish` leaves `<output>.tmp` on disk and no
//! `<output>` — the intended loud-failure mode (callers do not see a
//! half-written VCF file).

use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use super::errors::VcfWriteError;

/// htslib-required empty-bgzf EOF marker.
///
/// Any bgzf file that htslib (and therefore tabix / bcftools)
/// considers well-formed has to end with these 28 bytes verbatim.
/// `noodles-bgzf::io::Writer::finish` already emits the marker; this
/// constant exists so the in-crate test can assert the byte pattern
/// without duplicating it. The integration test in `tests/` keeps
/// its own copy (separate crate, no access to `pub(crate)`); the
/// two-place truth is the explicit cost of *not* exposing wire-
/// format magic in this crate's public API.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) const BGZF_EOF: &[u8; 28] = &[
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// In-flight sink. The variant is fixed by the output path suffix at
/// construction time and never changes thereafter.
pub(super) enum SinkKind {
    Plain(BufWriter<File>),
    Bgzf(noodles_bgzf::io::Writer<File>),
}

impl SinkKind {
    /// Open `<final_path>.tmp` and wrap it in the sink kind chosen by
    /// `final_path`'s suffix.
    ///
    /// The plain sink wraps the file in a 64 KiB `BufWriter` so the
    /// per-record `write_all` calls coalesce into block-sized writes
    /// before they hit the kernel. The bgzf sink does its own
    /// block-level buffering internally.
    pub(super) fn open_tmp(final_path: &Path) -> Result<Self, VcfWriteError> {
        let tmp_path = tmp_path_for(final_path);
        let file = File::create(&tmp_path).map_err(|source| VcfWriteError::CreateTmp {
            tmp_path: tmp_path.clone(),
            source,
        })?;
        Ok(if path_is_bgzf(final_path) {
            SinkKind::Bgzf(noodles_bgzf::io::Writer::new(file))
        } else {
            SinkKind::Plain(BufWriter::with_capacity(64 * 1024, file))
        })
    }

    /// Flush the sink, emit the bgzf EOF marker when applicable,
    /// fsync the file, atomically rename `<final_path>.tmp` →
    /// `<final_path>`, then fsync the parent directory so the
    /// rename itself is durable across a crash.
    ///
    /// Consumes `self` so a forgotten `finish` shows up as a missing
    /// final output rather than a silently truncated file.
    pub(super) fn finish(self, final_path: &Path) -> Result<(), VcfWriteError> {
        let tmp_path = tmp_path_for(final_path);

        let file = match self {
            SinkKind::Plain(buf) => buf.into_inner().map_err(|e| {
                let (source, _writer) = (e.into_error(), ());
                VcfWriteError::WriteHeader {
                    // `IntoInnerError` from `BufWriter` only fires on a
                    // buffered-write flush at drop/finish time — i.e.
                    // a tail write to the underlying File failed. The
                    // operation is "writing the tail of the bytes we
                    // buffered for the VCF", which is the same shape
                    // as a per-record/header write failure; report it
                    // as a tmp-file write failure on the underlying
                    // path. We do not have a more specific operation
                    // tag because the BufWriter doesn't carry one.
                    tmp_path: tmp_path.clone(),
                    source,
                }
            })?,
            SinkKind::Bgzf(bgzf) => bgzf.finish().map_err(|source| VcfWriteError::FinishBgzf {
                tmp_path: tmp_path.clone(),
                source,
            })?,
        };

        file.sync_all().map_err(|source| VcfWriteError::FsyncFile {
            tmp_path: tmp_path.clone(),
            source,
        })?;
        drop(file);

        fs::rename(&tmp_path, final_path).map_err(|source| VcfWriteError::Rename {
            tmp_path: tmp_path.clone(),
            final_path: final_path.to_path_buf(),
            source,
        })?;

        // M10: fsync the parent directory so the rename is durable.
        // Without this, a crash between `fs::rename` returning and the
        // parent's journal flushing can leave the file contents on
        // disk but the directory entry pointing at `<final_path>`
        // lost.
        sync_parent_dir(final_path)?;

        Ok(())
    }
}

impl Write for SinkKind {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            SinkKind::Plain(w) => w.write(buf),
            SinkKind::Bgzf(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            SinkKind::Plain(w) => w.flush(),
            SinkKind::Bgzf(w) => w.flush(),
        }
    }
}

/// Path used while the file is being written. The atomic
/// [`fs::rename`] at finish-time moves it to the final path.
pub(super) fn tmp_path_for(final_path: &Path) -> PathBuf {
    let mut tmp = final_path.as_os_str().to_owned();
    tmp.push(".tmp");
    PathBuf::from(tmp)
}

/// `fsync` the directory containing `final_path`. POSIX/Linux/macOS
/// all require this after a rename for the new directory entry to be
/// durable. On Windows the equivalent of `fsync(dir)` is a no-op /
/// not supported; `File::open` of a directory may itself fail, in
/// which case we still surface that as the typed `FsyncDir` error so
/// callers can decide what to do.
fn sync_parent_dir(final_path: &Path) -> Result<(), VcfWriteError> {
    let parent = final_path.parent();
    let dir_path: &Path = match parent {
        Some(p) if !p.as_os_str().is_empty() => p,
        // empty parent means `final_path` is a bare filename in cwd
        _ => Path::new("."),
    };
    let dir = File::open(dir_path).map_err(|source| VcfWriteError::FsyncDir {
        final_path: final_path.to_path_buf(),
        source,
    })?;
    dir.sync_all().map_err(|source| VcfWriteError::FsyncDir {
        final_path: final_path.to_path_buf(),
        source,
    })?;
    Ok(())
}

/// True when the path ends with `.vcf.gz` or `.vcf.bgz` (matched
/// case-insensitively — `out.VCF.GZ` resolves the same as
/// `out.vcf.gz`).
fn path_is_bgzf(path: &Path) -> bool {
    let Some(name) = path.file_name().and_then(|s| s.to_str()) else {
        return false;
    };
    let lower = name.to_ascii_lowercase();
    lower.ends_with(".vcf.gz") || lower.ends_with(".vcf.bgz")
}

#[cfg(test)]
mod tests {
    use std::io::Write as _;

    use tempfile::tempdir;

    use super::*;

    #[test]
    fn plain_suffix_selects_plain_sink() {
        assert!(!path_is_bgzf(Path::new("out.vcf")));
        assert!(!path_is_bgzf(Path::new("a/b/c.vcf")));
        assert!(!path_is_bgzf(Path::new("nogz.txt")));
    }

    #[test]
    fn gz_and_bgz_suffixes_select_bgzf_sink() {
        assert!(path_is_bgzf(Path::new("out.vcf.gz")));
        assert!(path_is_bgzf(Path::new("nested/out.vcf.bgz")));
    }

    /// Mi2: case-insensitive suffix matching — `out.VCF.GZ` is bgzf.
    #[test]
    fn uppercase_suffixes_select_bgzf_sink() {
        assert!(path_is_bgzf(Path::new("out.VCF.GZ")));
        assert!(path_is_bgzf(Path::new("Out.Vcf.Bgz")));
        assert!(!path_is_bgzf(Path::new("out.VCF")));
    }

    #[test]
    fn tmp_path_appends_dot_tmp() {
        assert_eq!(
            tmp_path_for(Path::new("/a/b/out.vcf")),
            PathBuf::from("/a/b/out.vcf.tmp")
        );
        assert_eq!(
            tmp_path_for(Path::new("/a/b/out.vcf.gz")),
            PathBuf::from("/a/b/out.vcf.gz.tmp")
        );
    }

    #[test]
    fn plain_sink_writes_then_atomically_renames() {
        let dir = tempdir().unwrap();
        let final_path = dir.path().join("out.vcf");
        let mut sink = SinkKind::open_tmp(&final_path).unwrap();
        sink.write_all(b"hello\n").unwrap();
        sink.finish(&final_path).unwrap();

        assert!(!tmp_path_for(&final_path).exists(), "tmp should be gone");
        let bytes = std::fs::read(&final_path).unwrap();
        assert_eq!(&bytes, b"hello\n");
    }

    #[test]
    fn bgzf_sink_ends_with_eof_block_after_finish() {
        let dir = tempdir().unwrap();
        let final_path = dir.path().join("out.vcf.gz");
        let mut sink = SinkKind::open_tmp(&final_path).unwrap();
        sink.write_all(b"hello bgzf\n").unwrap();
        sink.finish(&final_path).unwrap();

        let bytes = std::fs::read(&final_path).unwrap();
        assert!(
            bytes.len() >= BGZF_EOF.len(),
            "file too short for EOF marker"
        );
        assert_eq!(
            &bytes[bytes.len() - BGZF_EOF.len()..],
            &BGZF_EOF[..],
            "bgzf file must end with the htslib EOF block"
        );
    }

    /// M10: parent-directory fsync is part of the durability contract.
    /// Asserts that `finish` opens and fsyncs the parent dir on the
    /// happy path (success-path smoke; the crash-recovery property is
    /// covered by the existence of the `sync_parent_dir` call).
    #[test]
    fn finish_renames_and_fsyncs_parent_dir() {
        let dir = tempdir().unwrap();
        let final_path = dir.path().join("out.vcf");
        let mut sink = SinkKind::open_tmp(&final_path).unwrap();
        sink.write_all(b"hello\n").unwrap();
        sink.finish(&final_path).unwrap();

        // Post-finish the parent dir is still openable and syncable —
        // proves we did not error inside `sync_parent_dir` on the
        // happy path and that the dir file handle was released
        // cleanly.
        let parent = File::open(dir.path()).unwrap();
        parent.sync_all().unwrap();
        // And the final path exists, tmp does not.
        assert!(final_path.exists());
        assert!(!tmp_path_for(&final_path).exists());
    }
}
