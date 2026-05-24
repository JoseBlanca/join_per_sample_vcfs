//! Pure-Rust bgzf-aware concatenation of complete VCF fragments.
//!
//! The per-chromosome parallel path in
//! [`crate::pop_var_caller::cohort_driver`] drives one
//! [`CohortVcfWriter`](super::CohortVcfWriter) per chromosome, each
//! emitting a self-contained `.vcf` / `.vcf.gz` fragment. After all
//! per-chrom workers join, [`concat_fragments`] assembles those
//! fragments — in contig-table order — into the final cohort VCF:
//!
//! 1. Open `<final_path>.tmp` via [`SinkKind::open_tmp`], inheriting
//!    the writer's plain-vs-bgzf dispatch + atomic-rename discipline.
//! 2. For each fragment, walk lines through the matching reader
//!    ([`noodles_bgzf::io::Reader`] or [`BufReader<File>`]). The first
//!    fragment is copied verbatim; every subsequent fragment has its
//!    VCF header stripped — both the `##`-meta lines and the single
//!    `#CHROM\t…` line.
//! 3. Call [`SinkKind::finish`] — flush + fsync the tmp file, rename
//!    to `<final_path>`, fsync the parent directory.
//!
//! v1 decompresses + re-encodes each bgzf fragment. The cost is
//! bounded by total output size and is I/O-bound in practice; a
//! block-level header-surgery variant is tracked as a v2 follow-up
//! in the implementation plan.
//!
//! ## Why concat, not sort
//!
//! Each per-chrom fragment is internally monotonic by construction
//! (the per-chrom pipeline emits records in `(chrom_id, pos)` order +
//! the writer's per-contig monotonicity check asserts it). Walking
//! the contig table in declared order yields a final VCF that is
//! monotonic without a reorder buffer, in O(total_output_bytes)
//! streaming I/O — `O(n log n)` external-merge sort would pay for
//! something we already have.

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use noodles_bgzf::io::Reader as BgzfReader;

use super::errors::VcfWriteError;
use super::sink::{SinkKind, path_is_bgzf};
use crate::pop_var_caller::common::DEFAULT_BUFFERED_IO_CAPACITY;

/// Concatenate `fragments` (each a complete VCF with header) into a
/// single VCF at `final_path`, in the order they appear in
/// `fragments`. The first fragment's header is preserved; every
/// subsequent fragment's header lines (both `##`-meta and the single
/// `#CHROM\t…` line) are stripped, leaving only the data body.
///
/// The plain-vs-bgzf kind is selected by `final_path`'s suffix
/// (`.vcf.gz` / `.vcf.bgz` → bgzf; anything else → plain). Each
/// fragment **must** be the same kind — the concat helper opens
/// matching readers but does not transcode.
///
/// Atomic-write discipline is inherited from
/// [`SinkKind`](super::SinkKind): bytes go to `<final_path>.tmp`
/// first; on success the sink is flushed + fsynced, atomically
/// renamed into place, and the parent directory is fsynced so the
/// rename is durable across a crash.
///
/// Returns [`VcfWriteError::EmptyFragmentList`] when `fragments` is
/// empty — a downstream pipeline that produces zero per-chrom
/// fragments is a caller bug, not a quiet success.
pub(crate) fn concat_fragments(
    final_path: &Path,
    fragments: &[PathBuf],
) -> Result<(), VcfWriteError> {
    if fragments.is_empty() {
        return Err(VcfWriteError::EmptyFragmentList {
            final_path: final_path.to_path_buf(),
        });
    }

    let mut sink = SinkKind::open_tmp(final_path)?;
    let is_bgzf = path_is_bgzf(final_path);

    for (idx, fragment) in fragments.iter().enumerate() {
        let strip_header = idx > 0;
        let file = File::open(fragment).map_err(|source| VcfWriteError::ReadFragment {
            fragment: fragment.clone(),
            source,
        })?;
        if is_bgzf {
            let mut reader = BgzfReader::new(file);
            append_fragment_body(&mut reader, &mut sink, strip_header, final_path, fragment)?;
        } else {
            let mut reader = BufReader::with_capacity(DEFAULT_BUFFERED_IO_CAPACITY, file);
            append_fragment_body(&mut reader, &mut sink, strip_header, final_path, fragment)?;
        }
    }

    sink.finish(final_path)
}

fn append_fragment_body<R: BufRead>(
    reader: &mut R,
    sink: &mut SinkKind,
    strip_header: bool,
    final_path: &Path,
    fragment: &Path,
) -> Result<(), VcfWriteError> {
    let mut buf = String::new();
    loop {
        buf.clear();
        let n = reader
            .read_line(&mut buf)
            .map_err(|source| VcfWriteError::ReadFragment {
                fragment: fragment.to_path_buf(),
                source,
            })?;
        if n == 0 {
            return Ok(());
        }
        if strip_header && buf.starts_with('#') {
            continue;
        }
        sink.write_all(buf.as_bytes())
            .map_err(|source| VcfWriteError::WriteConcat {
                final_path: final_path.to_path_buf(),
                fragment: fragment.to_path_buf(),
                source,
            })?;
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write as _;

    use noodles_bgzf::io::Writer as BgzfWriter;
    use tempfile::tempdir;

    use super::*;

    /// Two header lines + two body lines, plain text.
    const FRAG_A: &[u8] = b"\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tT\t.\tPASS\t.
1\t200\t.\tG\tC\t.\tPASS\t.
";

    const FRAG_B: &[u8] = b"\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
2\t50\t.\tT\tA\t.\tPASS\t.
2\t300\t.\tC\tG\t.\tPASS\t.
";

    const FRAG_HEADER_ONLY: &[u8] = b"\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

    fn write_plain(path: &Path, bytes: &[u8]) {
        let mut f = File::create(path).unwrap();
        f.write_all(bytes).unwrap();
        f.sync_all().unwrap();
    }

    fn write_bgzf(path: &Path, bytes: &[u8]) {
        let mut w = BgzfWriter::new(File::create(path).unwrap());
        w.write_all(bytes).unwrap();
        w.finish().unwrap();
    }

    fn read_bgzf(path: &Path) -> Vec<u8> {
        use std::io::Read as _;
        let mut buf = Vec::new();
        BgzfReader::new(File::open(path).unwrap())
            .read_to_end(&mut buf)
            .unwrap();
        buf
    }

    #[test]
    fn concat_plain_two_fragments_strips_second_header() {
        let dir = tempdir().unwrap();
        let frag1 = dir.path().join("chr_000_one.vcf");
        let frag2 = dir.path().join("chr_001_two.vcf");
        let final_path = dir.path().join("out.vcf");
        write_plain(&frag1, FRAG_A);
        write_plain(&frag2, FRAG_B);

        concat_fragments(&final_path, &[frag1, frag2]).unwrap();

        let got = std::fs::read(&final_path).unwrap();
        let expected = b"\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tT\t.\tPASS\t.
1\t200\t.\tG\tC\t.\tPASS\t.
2\t50\t.\tT\tA\t.\tPASS\t.
2\t300\t.\tC\tG\t.\tPASS\t.
";
        assert_eq!(&got[..], expected);
    }

    #[test]
    fn concat_bgzf_two_fragments_round_trip() {
        let dir = tempdir().unwrap();
        let frag1 = dir.path().join("chr_000_one.vcf.gz");
        let frag2 = dir.path().join("chr_001_two.vcf.gz");
        let final_path = dir.path().join("out.vcf.gz");
        write_bgzf(&frag1, FRAG_A);
        write_bgzf(&frag2, FRAG_B);

        concat_fragments(&final_path, &[frag1, frag2]).unwrap();

        let got = read_bgzf(&final_path);
        let expected = b"\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tT\t.\tPASS\t.
1\t200\t.\tG\tC\t.\tPASS\t.
2\t50\t.\tT\tA\t.\tPASS\t.
2\t300\t.\tC\tG\t.\tPASS\t.
";
        assert_eq!(&got[..], expected);
    }

    /// N=1 fragment ⇒ identity. The first-fragment branch keeps the
    /// header; the contents of the final output must equal the
    /// original fragment.
    #[test]
    fn concat_single_fragment_is_identity() {
        let dir = tempdir().unwrap();
        let frag = dir.path().join("chr_000_only.vcf");
        let final_path = dir.path().join("out.vcf");
        write_plain(&frag, FRAG_A);

        concat_fragments(&final_path, &[frag]).unwrap();

        let got = std::fs::read(&final_path).unwrap();
        assert_eq!(&got[..], FRAG_A);
    }

    /// A middle fragment with a header but zero body lines must
    /// concatenate as a no-op (the header is stripped, no body is
    /// written). Mimics the empty-chromosome case from the
    /// integration tests.
    #[test]
    fn concat_skips_header_only_middle_fragment() {
        let dir = tempdir().unwrap();
        let frag1 = dir.path().join("chr_000_one.vcf");
        let frag2 = dir.path().join("chr_001_empty.vcf");
        let frag3 = dir.path().join("chr_002_two.vcf");
        let final_path = dir.path().join("out.vcf");
        write_plain(&frag1, FRAG_A);
        write_plain(&frag2, FRAG_HEADER_ONLY);
        write_plain(&frag3, FRAG_B);

        concat_fragments(&final_path, &[frag1, frag2, frag3]).unwrap();

        let got = std::fs::read(&final_path).unwrap();
        let expected = b"\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tT\t.\tPASS\t.
1\t200\t.\tG\tC\t.\tPASS\t.
2\t50\t.\tT\tA\t.\tPASS\t.
2\t300\t.\tC\tG\t.\tPASS\t.
";
        assert_eq!(&got[..], expected);
    }

    /// Empty fragment list is a caller bug — the cohort's contig
    /// table is non-empty by construction, so producing zero
    /// fragments means the orchestrator built the wrong list. Surface
    /// the typed error rather than writing an empty file.
    #[test]
    fn concat_empty_fragments_list_is_typed_error() {
        let dir = tempdir().unwrap();
        let final_path = dir.path().join("out.vcf");

        let err = concat_fragments(&final_path, &[]).unwrap_err();
        match err {
            VcfWriteError::EmptyFragmentList { final_path: fp } => {
                assert_eq!(fp, final_path);
            }
            other => panic!("expected EmptyFragmentList, got {other:?}"),
        }
        // No tmp / final file should have been created.
        assert!(!final_path.exists());
        let tmp = super::super::tmp_path_for(&final_path);
        assert!(!tmp.exists());
    }
}
