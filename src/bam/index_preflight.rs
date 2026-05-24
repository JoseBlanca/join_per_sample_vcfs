//! Alignment-index pre-flight: detect `.crai` next to each input
//! CRAM, and optionally build the index in place if it is missing.
//!
//! Per-chromosome parallelism in
//! [`crate::pop_var_caller::var_calling_from_bam`] requires an
//! alignment index next to every input alignment file so each rayon
//! worker can issue a contig-scoped `noodles` `Reader::query(...)`.
//! This module is the gate that enforces the requirement: callers
//! either pass `build_if_missing = false` (hard-error on the first
//! input without an index, intended for non-interactive pipelines
//! that should never silently create files next to user inputs), or
//! `build_if_missing = true` (build missing indexes in place, one
//! per missing input, with a single-line stderr progress message
//! per build).
//!
//! Design and policy: see
//! `doc/devel/implementation_plans/var_calling_from_bam_per_chromosome.md`
//! §2 "Index pre-flight".
//!
//! Today the only alignment-file format we read is CRAM (matches
//! [`crate::bam::alignment_input::AlignmentMergedReader`]). When BAM input
//! support lands, extend [`AlignmentFileKind`] with a `Bam` variant
//! and add the corresponding `.bai` / `.csi` detection + build
//! branches; the rest of the pre-flight shape stays the same.
//!
//! `noodles_cram::fs::index` stream-walks the source file once to
//! collect container offsets, so an index build costs roughly one
//! extra pass over the CRAM's bytes.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::bam::errors::AlignmentIndexError;

// ---------------------------------------------------------------------
// Pre-loaded alignment index (shared across rayon workers)
// ---------------------------------------------------------------------

/// Pre-loaded alignment index for a single input file, shareable
/// across rayon workers via the `Arc` payload.
///
/// Drivers that need per-chromosome random access call
/// [`load_alignment_index`] once per input at startup, wrap the
/// result in this enum, and pass a `&[AlignmentIndex]` (or
/// `Vec<AlignmentIndex>`) into
/// [`crate::bam::alignment_input::AlignmentMergedReader::query`]. Each rayon
/// worker pays one `Arc::clone` to get its own handle on the
/// already-parsed index; no per-worker disk reads or re-parsing.
///
/// Single-variant today (CRAM-only); when BAM input support lands,
/// add `Bai(Arc<noodles_bam::bai::Index>)` and
/// `Csi(Arc<noodles_csi::Index>)`.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum AlignmentIndex {
    Crai(Arc<noodles_cram::crai::Index>),
}

/// Load an alignment index from disk. The expected location is the
/// canonical sibling path (`<input>.crai` for CRAM); for missing or
/// malformed indexes, prefer running [`preflight_alignment_indexes`]
/// first so the failure surfaces at the typed
/// [`AlignmentIndexError`] layer rather than as a bare `io::Error`
/// here.
pub fn load_alignment_index(input: &Path) -> std::io::Result<AlignmentIndex> {
    match AlignmentFileKind::from_path(input) {
        Some(AlignmentFileKind::Cram) => {
            let index = noodles_cram::crai::fs::read(crai_path_for(input))?;
            Ok(AlignmentIndex::Crai(Arc::new(index)))
        }
        None => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!(
                "unsupported alignment-file extension for '{}'",
                input.display()
            ),
        )),
    }
}

// ---------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------

/// Confirm every input in `inputs` has an alignment index next to it,
/// optionally building any missing index when `build_if_missing` is
/// true.
///
/// Returns `Ok(())` once every input is index-backed.
///
/// # Errors
///
/// - [`AlignmentIndexError::MissingAlignmentIndex`] — `build_if_missing`
///   was false and at least one input has no index. The variant
///   carries the source path and the index path we looked for.
/// - [`AlignmentIndexError::BuildFailed`] — `build_if_missing` was
///   true and an index build (or write) failed.
/// - [`AlignmentIndexError::UnsupportedExtension`] — an input's
///   extension is not `.cram`; we do not guess the file format here.
pub fn preflight_alignment_indexes(
    inputs: &[PathBuf],
    build_if_missing: bool,
) -> Result<(), AlignmentIndexError> {
    // First pass: classify + find missing inputs. Done before any
    // build so a `--build-map-file-index` run with one unbuildable
    // input does not partially build the others before failing.
    let mut missing: Vec<(usize, AlignmentFileKind)> = Vec::new();
    for (idx, input) in inputs.iter().enumerate() {
        let kind = AlignmentFileKind::from_path(input).ok_or_else(|| {
            AlignmentIndexError::UnsupportedExtension {
                path: input.clone(),
            }
        })?;
        if existing_index_for(input, kind).is_none() {
            missing.push((idx, kind));
        }
    }

    if missing.is_empty() {
        return Ok(());
    }

    if !build_if_missing {
        let (idx, kind) = missing[0];
        let input = &inputs[idx];
        return Err(AlignmentIndexError::MissingAlignmentIndex {
            path: input.clone(),
            expected_index_path: target_index_path(input, kind),
        });
    }

    let total = missing.len();
    for (sequence, (idx, kind)) in missing.into_iter().enumerate() {
        let input = &inputs[idx];
        eprintln!(
            "building index for {} ({} of {})...",
            input.display(),
            sequence + 1,
            total,
        );
        build_index(input, kind).map_err(|source| AlignmentIndexError::BuildFailed {
            path: input.clone(),
            source,
        })?;
    }

    Ok(())
}

// ---------------------------------------------------------------------
// Internals
// ---------------------------------------------------------------------

/// Classification of a mapped-read input file, derived from its
/// extension only. We do not sniff the file's magic bytes — an
/// explicit extension is required.
///
/// Single-variant today; the `enum` keeps the dispatch surface
/// stable for when BAM support lands.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AlignmentFileKind {
    Cram,
}

impl AlignmentFileKind {
    fn from_path(path: &Path) -> Option<Self> {
        match path.extension().and_then(|ext| ext.to_str()) {
            Some("cram") => Some(Self::Cram),
            _ => None,
        }
    }
}

/// Resolve the alignment index that already exists for `input`, if
/// any.
fn existing_index_for(input: &Path, kind: AlignmentFileKind) -> Option<PathBuf> {
    match kind {
        AlignmentFileKind::Cram => {
            let p = crai_path_for(input);
            p.exists().then_some(p)
        }
    }
}

/// The index path we would *build* for `input` if asked.
fn target_index_path(input: &Path, kind: AlignmentFileKind) -> PathBuf {
    match kind {
        AlignmentFileKind::Cram => crai_path_for(input),
    }
}

fn crai_path_for(cram: &Path) -> PathBuf {
    append_extension(cram, "crai")
}

/// Append `.<ext>` to a path. Unlike [`Path::with_extension`], this
/// keeps the existing extension intact — `sample.cram` becomes
/// `sample.cram.crai`, not `sample.crai`.
fn append_extension(base: &Path, ext: &str) -> PathBuf {
    let mut s = base.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

fn build_index(input: &Path, kind: AlignmentFileKind) -> std::io::Result<()> {
    match kind {
        AlignmentFileKind::Cram => {
            let index = noodles_cram::fs::index(input)?;
            noodles_cram::crai::fs::write(crai_path_for(input), &index)
        }
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::fs::{self, File};
    use std::io::Write;

    use tempfile::TempDir;

    use super::*;

    /// Write a one-byte placeholder file. Pre-flight only checks
    /// path existence on the no-build path; the file's content does
    /// not need to be a valid CRAM for these tests.
    fn touch(path: &Path) {
        let mut f = File::create(path).expect("create placeholder");
        f.write_all(b"x").expect("write placeholder");
    }

    #[test]
    fn preflight_accepts_existing_crai() {
        let dir = TempDir::new().expect("tempdir");
        let cram = dir.path().join("sample.cram");
        let crai = dir.path().join("sample.cram.crai");
        touch(&cram);
        touch(&crai);

        preflight_alignment_indexes(&[cram], false).expect("preflight ok");
    }

    #[test]
    fn preflight_errors_when_missing_and_flag_unset() {
        let dir = TempDir::new().expect("tempdir");
        let cram = dir.path().join("sample.cram");
        touch(&cram);

        let err = preflight_alignment_indexes(std::slice::from_ref(&cram), false)
            .expect_err("preflight should fail without index");

        match err {
            AlignmentIndexError::MissingAlignmentIndex {
                path,
                expected_index_path,
            } => {
                assert_eq!(path, cram);
                assert_eq!(expected_index_path, dir.path().join("sample.cram.crai"));
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn preflight_errors_on_unsupported_extension() {
        let dir = TempDir::new().expect("tempdir");
        let input = dir.path().join("sample.sam");
        touch(&input);

        let err = preflight_alignment_indexes(std::slice::from_ref(&input), false)
            .expect_err("preflight should reject .sam");
        assert!(matches!(
            err,
            AlignmentIndexError::UnsupportedExtension { path } if path == input
        ));
    }

    #[test]
    fn preflight_idempotent_when_index_exists_and_flag_set() {
        let dir = TempDir::new().expect("tempdir");
        let cram = dir.path().join("sample.cram");
        let crai = dir.path().join("sample.cram.crai");
        touch(&cram);
        touch(&crai);

        let before = fs::metadata(&crai).expect("crai metadata").modified().ok();

        preflight_alignment_indexes(&[cram], true).expect("preflight ok");

        let after = fs::metadata(&crai).expect("crai metadata").modified().ok();
        assert_eq!(before, after, "existing index must not be rebuilt");
    }

    #[test]
    fn preflight_reports_first_missing_input() {
        let dir = TempDir::new().expect("tempdir");
        let cram_a = dir.path().join("a.cram");
        let cram_b = dir.path().join("b.cram");
        let crai_a = dir.path().join("a.cram.crai");
        touch(&cram_a);
        touch(&crai_a);
        touch(&cram_b);
        // b has no .crai.

        let err = preflight_alignment_indexes(&[cram_a, cram_b.clone()], false)
            .expect_err("preflight should fail on missing b");
        match err {
            AlignmentIndexError::MissingAlignmentIndex { path, .. } => {
                assert_eq!(path, cram_b, "first missing input must be reported");
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn target_index_path_appends_to_existing_extension() {
        let cram = Path::new("/tmp/sample.cram");
        let crai = target_index_path(cram, AlignmentFileKind::Cram);
        assert_eq!(crai, PathBuf::from("/tmp/sample.cram.crai"));
    }
}
