//! Output-sink dispatch for the cohort VCF writer.
//!
//! The writer streams its bytes through a single `std::io::Write`-
//! shaped sink that hides the plain-vs-bgzf branch from the encoder.
//! Path suffix selects the sink kind:
//!
//! * `.vcf.gz` or `.vcf.bgz` → [`SinkKind::Bgzf`]
//! * anything else           → [`SinkKind::Plain`]
//!
//! Bytes go to `<output>.tmp` first; [`SinkKind::finish`] flushes,
//! emits the bgzf EOF block when applicable, syncs, and atomically
//! renames the tmp path to the final output path. A crash before
//! `finish` leaves `<output>.tmp` on disk and no `<output>`, which
//! is the intended loud-failure mode (callers do not see a
//! half-written VCF file).

use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use super::errors::VcfWriteError;

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
        let file = File::create(&tmp_path)?;
        Ok(if path_is_bgzf(final_path) {
            SinkKind::Bgzf(noodles_bgzf::io::Writer::new(file))
        } else {
            SinkKind::Plain(BufWriter::with_capacity(64 * 1024, file))
        })
    }

    /// Flush the sink, emit the bgzf EOF marker when applicable,
    /// `fsync`, and atomically rename `<final_path>.tmp` →
    /// `<final_path>`.
    ///
    /// Consumes `self` so a forgotten `finish` shows up as a missing
    /// final output rather than a silently truncated file.
    pub(super) fn finish(self, final_path: &Path) -> Result<(), VcfWriteError> {
        let file = match self {
            SinkKind::Plain(buf) => buf
                .into_inner()
                .map_err(|e| VcfWriteError::Io(e.into_error()))?,
            SinkKind::Bgzf(bgzf) => bgzf.finish()?,
        };
        file.sync_all()?;
        drop(file);
        let tmp_path = tmp_path_for(final_path);
        fs::rename(&tmp_path, final_path)?;
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

/// True when the path ends with `.vcf.gz` or `.vcf.bgz` — both are
/// htslib-style bgzipped VCFs.
fn path_is_bgzf(path: &Path) -> bool {
    let name = match path.file_name().and_then(|s| s.to_str()) {
        Some(n) => n,
        None => return false,
    };
    name.ends_with(".vcf.gz") || name.ends_with(".vcf.bgz")
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

    /// The htslib-required empty-bgzf EOF marker. Any bgzf file that
    /// htslib (and therefore tabix / bcftools) considers well-formed
    /// has to end with these 28 bytes verbatim.
    const BGZF_EOF: &[u8] = &[
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

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
            BGZF_EOF,
            "bgzf file must end with the htslib EOF block"
        );
    }
}
