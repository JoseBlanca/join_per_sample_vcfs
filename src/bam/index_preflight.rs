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
//! Both CRAM (`.cram` + `.crai`) and BAM (`.bam` + `.csi`/`.bai`)
//! inputs are accepted. The classification + missing-index pass
//! also rejects mixed `.cram` + `.bam` inputs in one invocation —
//! the merge downstream wants one format per sample and the
//! all-or-nothing rule keeps both code and error messages clean.
//!
//! `noodles_cram::fs::index` stream-walks the source file once to
//! collect container offsets, so an index build costs roughly one
//! extra pass over the CRAM's bytes.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::bam::errors::AlignmentIndexError;

// ---------------------------------------------------------------------
// BAM index policy constants
// ---------------------------------------------------------------------

/// Read-side preference order when both `.csi` and `.bai` exist
/// next to a BAM input. `.csi` wins because it has no 512 Mbp
/// per-contig length cap (some plant references exceed this);
/// `.bai` is accepted on read as a fallback for compatibility
/// with existing data sets. See
/// [`existing_index_for`] for the consuming policy and the
/// `BAM_INDEX_BUILD_FORMAT` constant below for the build-side
/// pin.
pub const BAM_INDEX_READ_PREFERENCE: &[&str] = &["csi", "bai"];

/// On-disk format produced by [`preflight_alignment_indexes`]
/// when `--build-map-file-index` is set on a BAM input. Always
/// `.csi`: at `CSI_DEPTH = 6` it addresses contigs up to ~2^32
/// bp (~4.3 Gbp), well past the `.bai` 2^29 (~537 Mbp) cap that
/// `.bai`'s 16 kbp bin grid imposes.
pub const BAM_INDEX_BUILD_FORMAT: &str = "csi";

/// `.csi` min-shift exponent. `14` → 16 kbp leaf bins (matches
/// `.bai` resolution). Bumping widens leaf bins (smaller index,
/// coarser queries); lowering narrows them. Default 14 is the
/// noodles / htslib convention.
pub const CSI_MIN_SHIFT: u8 = 14;

/// `.csi` binning depth. Addressable contig length is
/// approximately `2^(min_shift + 3 * depth)`. At our
/// (`min_shift = 14`, `depth = 6`) the addressable length is
/// ~2^32 ≈ 4.3 Gbp — well past the `.bai` cap, safe for
/// large plant references the project may target in future
/// (wheat ~700 Mbp, lily ~3 Gbp). Increase only if a use case
/// surfaces with contigs > 4 Gbp.
pub const CSI_DEPTH: u8 = 6;

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
/// For BAM inputs we keep the on-disk format we loaded — either
/// `.csi` or `.bai` — because the cohort driver wraps the payload
/// in a per-format BAM-side enum before handing it to the
/// indexed-record-stream opener. Both index formats satisfy
/// noodles' [`noodles_csi::binning_index::BinningIndex`] trait, so
/// the merge cares only about the on-disk file that's available
/// locally.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum AlignmentIndex {
    Crai(Arc<noodles_cram::crai::Index>),
    BamCsi(Arc<noodles_csi::Index>),
    BamBai(Arc<noodles_bam::bai::Index>),
}

impl AlignmentIndex {
    /// User-facing name for an [`AlignmentIndex`] variant, used in
    /// [`crate::bam::errors::AlignmentInputError::AlignmentIndexFormatMismatch`]
    /// messages. Mirror of [`AlignmentFileKind::display_name`].
    pub(crate) fn display_name(&self) -> &'static str {
        match self {
            Self::Crai(_) => "CRAI",
            Self::BamCsi(_) => "CSI",
            Self::BamBai(_) => "BAI",
        }
    }
}

/// Load an alignment index from disk. The expected location is the
/// canonical sibling path (`<input>.crai` for CRAM, `<input>.csi`
/// preferred or `<input>.bai` fallback for BAM); missing or
/// malformed indexes surface as the typed
/// [`AlignmentIndexError`] variants `MissingAlignmentIndex`,
/// `LoadFailed`, or `UnsupportedExtension`.
pub fn load_alignment_index(input: &Path) -> Result<AlignmentIndex, AlignmentIndexError> {
    let kind = AlignmentFileKind::from_path(input).ok_or_else(|| {
        AlignmentIndexError::UnsupportedExtension {
            path: input.to_path_buf(),
        }
    })?;
    // Reuse `existing_index_for` so the `.csi`-preferred /
    // `.bai`-fallback policy lives in exactly one place (M13 +
    // Mi11). For BAM, this returns the actual on-disk path that
    // exists; we then parse it with the matching index reader.
    let index_path = existing_index_for(input, kind).ok_or_else(|| {
        AlignmentIndexError::MissingAlignmentIndex {
            path: input.to_path_buf(),
            expected_index_path: target_index_path(input, kind),
        }
    })?;
    match kind {
        AlignmentFileKind::Cram => {
            let index = noodles_cram::crai::fs::read(&index_path).map_err(|source| {
                AlignmentIndexError::LoadFailed {
                    path: input.to_path_buf(),
                    index_path: index_path.clone(),
                    source,
                }
            })?;
            Ok(AlignmentIndex::Crai(Arc::new(index)))
        }
        AlignmentFileKind::Bam => {
            // `existing_index_for` picks `.csi` when both formats
            // are present (Mi11 / M13); dispatch on the resolved
            // path's extension.
            let ext = index_path.extension().and_then(|e| e.to_str());
            if ext == Some("csi") {
                let index = noodles_csi::fs::read(&index_path).map_err(|source| {
                    AlignmentIndexError::LoadFailed {
                        path: input.to_path_buf(),
                        index_path: index_path.clone(),
                        source,
                    }
                })?;
                Ok(AlignmentIndex::BamCsi(Arc::new(index)))
            } else {
                let index = noodles_bam::bai::fs::read(&index_path).map_err(|source| {
                    AlignmentIndexError::LoadFailed {
                        path: input.to_path_buf(),
                        index_path: index_path.clone(),
                        source,
                    }
                })?;
                Ok(AlignmentIndex::BamBai(Arc::new(index)))
            }
        }
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
///   extension is neither `.cram` nor `.bam`; we do not guess the
///   file format here.
/// - [`AlignmentIndexError::MixedAlignmentFileFormats`] — the
///   input list mixes `.cram` and `.bam` files. One format per
///   invocation only; mixed-format support is an explicit
///   non-goal in `doc/devel/implementation_plans/bam_input_support.md`.
pub fn preflight_alignment_indexes(
    inputs: &[PathBuf],
    build_if_missing: bool,
) -> Result<(), AlignmentIndexError> {
    // First pass: classify + find missing inputs + reject mixed
    // formats. Done before any build so a `--build-map-file-index`
    // run with one unbuildable input does not partially build the
    // others before failing.
    let mut missing: Vec<(usize, AlignmentFileKind)> = Vec::new();
    let mut first_seen: Option<(usize, AlignmentFileKind)> = None;
    for (idx, input) in inputs.iter().enumerate() {
        let kind = AlignmentFileKind::from_path(input).ok_or_else(|| {
            AlignmentIndexError::UnsupportedExtension {
                path: input.clone(),
            }
        })?;
        match first_seen {
            None => first_seen = Some((idx, kind)),
            Some((first_idx, first_kind)) if first_kind != kind => {
                return Err(AlignmentIndexError::MixedAlignmentFileFormats {
                    first_path: inputs[first_idx].clone(),
                    first_format: first_kind.display_name(),
                    other_path: input.clone(),
                    other_format: kind.display_name(),
                });
            }
            Some(_) => {}
        }
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
/// `pub(crate)` so the
/// [`crate::bam::alignment_input::AlignmentMergedReader::new`]
/// dispatch can use the same classifier when opening the per-input
/// record streams (avoids two different "which format is this
/// path" rules drifting apart).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlignmentFileKind {
    Cram,
    Bam,
}

impl AlignmentFileKind {
    pub(crate) fn from_path(path: &Path) -> Option<Self> {
        match path.extension().and_then(|ext| ext.to_str()) {
            Some("cram") => Some(Self::Cram),
            Some("bam") => Some(Self::Bam),
            _ => None,
        }
    }

    /// User-facing format name embedded in error messages. Same
    /// strings appear in [`AlignmentIndexError`] /
    /// [`crate::bam::errors::AlignmentInputError`] mixed-format
    /// messages.
    pub(crate) fn display_name(self) -> &'static str {
        match self {
            Self::Cram => "CRAM",
            Self::Bam => "BAM",
        }
    }
}

/// Resolve the alignment index that already exists for `input`,
/// if any. For BAM, walks [`BAM_INDEX_READ_PREFERENCE`] in order
/// and returns the first format whose file is on disk
/// (`.csi`-preferred, `.bai`-fallback).
fn existing_index_for(input: &Path, kind: AlignmentFileKind) -> Option<PathBuf> {
    match kind {
        AlignmentFileKind::Cram => {
            let p = crai_path_for(input);
            p.exists().then_some(p)
        }
        AlignmentFileKind::Bam => BAM_INDEX_READ_PREFERENCE
            .iter()
            .map(|ext| bam_index_path_for(input, ext))
            .find(|p| p.exists()),
    }
}

/// The index path we would *build* for `input` if asked. For BAM
/// this is always `<input>.<BAM_INDEX_BUILD_FORMAT>` (`.csi`) —
/// see the constant's doc for rationale.
fn target_index_path(input: &Path, kind: AlignmentFileKind) -> PathBuf {
    match kind {
        AlignmentFileKind::Cram => crai_path_for(input),
        AlignmentFileKind::Bam => bam_index_path_for(input, BAM_INDEX_BUILD_FORMAT),
    }
}

fn crai_path_for(cram: &Path) -> PathBuf {
    append_extension(cram, "crai")
}

/// `<bam>.<ext>` for any BAM-side index extension (`csi`, `bai`).
/// Used by both the read-side preference walk and the build-side
/// canonical-path helper so the path-construction rule lives in
/// one place.
fn bam_index_path_for(bam: &Path, ext: &str) -> PathBuf {
    append_extension(bam, ext)
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
        AlignmentFileKind::Bam => build_csi_for_bam(input),
    }
}

/// Build a `.csi` next to `bam_path` by scanning the BAM once and
/// feeding the chunks into a [`noodles_csi::binning_index::Indexer`]
/// parameterised by [`CSI_MIN_SHIFT`] / [`CSI_DEPTH`].
fn build_csi_for_bam(bam_path: &Path) -> std::io::Result<()> {
    use noodles_csi::binning_index::Indexer;
    // BinnedIndex is the CSI on-disk shape; LinearIndex would emit
    // the .bai shape. We only build .csi here — see
    // BAM_INDEX_BUILD_FORMAT / CSI_DEPTH constants.
    use noodles_csi::binning_index::index::reference_sequence::index::BinnedIndex;

    let mut reader = noodles_bam::io::reader::Builder.build_from_path(bam_path)?;
    let header = reader.read_header()?;

    // Bypass `Indexer::default()` (which picks min_shift=14 /
    // depth=5, capping addressable contig length at ~537 Mbp).
    // We use the named constants explicitly so the bump to
    // depth=6 (~4.3 Gbp cap) is visible at the build site
    // rather than hidden behind a no-argument default.
    let mut indexer: Indexer<BinnedIndex> = Indexer::new(CSI_MIN_SHIFT, CSI_DEPTH);
    populate_binning_index(&mut reader, &mut indexer)?;

    let index = indexer.build(header.reference_sequences().len());
    noodles_csi::fs::write(bam_index_path_for(bam_path, BAM_INDEX_BUILD_FORMAT), &index)
}

/// Walk `reader` once and add a chunk per record to `indexer`.
/// Shared body between [`build_csi_for_bam`] (production CSI
/// build) and `crate::bam::bam_input::tests::build_bai_in_memory`
/// (test fixture, LinearIndex parameterisation for `.bai`). Same
/// chunk_start / chunk_end / alignment-context shape across the
/// two callers; the indexer's type parameter is the only thing
/// that differs.
///
/// Integration-test fixtures in `tests/common/mod.rs` carry their
/// own copy because the test crate is outside `crate::bam` and
/// cannot reach `pub(crate)` items. Reducing the production +
/// in-crate-test duplication to one helper is the achievable
/// part of the Mi9 cleanup; the integration-test copy stays.
pub(crate) fn populate_binning_index<I>(
    reader: &mut noodles_bam::io::Reader<noodles_bgzf::io::Reader<std::fs::File>>,
    indexer: &mut noodles_csi::binning_index::Indexer<I>,
) -> std::io::Result<()>
where
    I: noodles_csi::binning_index::index::reference_sequence::Index + Default,
{
    use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles_sam::alignment::Record as _;

    let mut chunk_start = reader.get_ref().virtual_position();
    let mut record = noodles_bam::Record::default();

    while reader.read_record(&mut record)? != 0 {
        let chunk_end = reader.get_ref().virtual_position();
        let alignment_context = match (
            record.reference_sequence_id().transpose()?,
            record.alignment_start().transpose()?,
            record.alignment_end().transpose()?,
        ) {
            (Some(id), Some(start), Some(end)) => {
                let is_mapped = !record.flags().is_unmapped();
                Some((id, start, end, is_mapped))
            }
            _ => None,
        };
        let chunk = Chunk::new(chunk_start, chunk_end);
        indexer.add_record(alignment_context, chunk)?;
        chunk_start = chunk_end;
    }
    Ok(())
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

    // --- BAM-side tests ----------------------------------------------

    #[test]
    fn target_index_path_for_bam_picks_csi_not_bai() {
        // Build-side: always emits .csi (the format with no
        // contig-length cap).
        let bam = Path::new("/tmp/sample.bam");
        let csi = target_index_path(bam, AlignmentFileKind::Bam);
        assert_eq!(csi, PathBuf::from("/tmp/sample.bam.csi"));
    }

    #[test]
    fn preflight_accepts_existing_csi() {
        let dir = TempDir::new().expect("tempdir");
        let bam = dir.path().join("sample.bam");
        let csi = dir.path().join("sample.bam.csi");
        touch(&bam);
        touch(&csi);

        preflight_alignment_indexes(&[bam], false).expect("preflight ok with .csi");
    }

    #[test]
    fn preflight_accepts_existing_bai_as_fallback() {
        let dir = TempDir::new().expect("tempdir");
        let bam = dir.path().join("sample.bam");
        let bai = dir.path().join("sample.bam.bai");
        touch(&bam);
        touch(&bai);

        preflight_alignment_indexes(&[bam], false).expect("preflight ok with .bai");
    }

    #[test]
    fn preflight_prefers_csi_over_bai_when_both_present() {
        let dir = TempDir::new().expect("tempdir");
        let bam = dir.path().join("sample.bam");
        let csi = dir.path().join("sample.bam.csi");
        let bai = dir.path().join("sample.bam.bai");
        touch(&bam);
        touch(&csi);
        touch(&bai);

        // Both exist; existing_index_for must return the .csi one.
        let picked = existing_index_for(&bam, AlignmentFileKind::Bam).expect("some index");
        assert_eq!(picked, csi, ".csi must win over .bai when both exist");
    }

    #[test]
    fn preflight_errors_when_bam_missing_index_and_flag_unset() {
        let dir = TempDir::new().expect("tempdir");
        let bam = dir.path().join("sample.bam");
        touch(&bam);

        let err = preflight_alignment_indexes(std::slice::from_ref(&bam), false)
            .expect_err("preflight should fail without index");
        match err {
            AlignmentIndexError::MissingAlignmentIndex {
                path,
                expected_index_path,
            } => {
                assert_eq!(path, bam);
                // Build path is .csi-canonical (see
                // target_index_path); the missing-index error
                // points the user at the file the build would
                // create.
                assert_eq!(expected_index_path, dir.path().join("sample.bam.csi"));
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn preflight_errors_on_mixed_cram_and_bam() {
        let dir = TempDir::new().expect("tempdir");
        let cram = dir.path().join("first.cram");
        let bam = dir.path().join("second.bam");
        // Touch all four to keep the test focused on the
        // mixed-format check (no missing-index noise).
        touch(&cram);
        touch(&dir.path().join("first.cram.crai"));
        touch(&bam);
        touch(&dir.path().join("second.bam.csi"));

        let err = preflight_alignment_indexes(&[cram.clone(), bam.clone()], false)
            .expect_err("preflight should reject mixed cram + bam");
        match err {
            AlignmentIndexError::MixedAlignmentFileFormats {
                first_path,
                first_format,
                other_path,
                other_format,
            } => {
                assert_eq!(first_path, cram);
                assert_eq!(first_format, "CRAM");
                assert_eq!(other_path, bam);
                assert_eq!(other_format, "BAM");
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn preflight_rejects_bam_first_then_cram_too() {
        // Reverse the input order to confirm the check is
        // order-independent (and that the first-seen rule fires
        // on whichever format appears first).
        let dir = TempDir::new().expect("tempdir");
        let bam = dir.path().join("first.bam");
        let cram = dir.path().join("second.cram");
        touch(&bam);
        touch(&dir.path().join("first.bam.csi"));
        touch(&cram);
        touch(&dir.path().join("second.cram.crai"));

        let err = preflight_alignment_indexes(&[bam.clone(), cram.clone()], false)
            .expect_err("preflight should reject mixed bam + cram");
        match err {
            AlignmentIndexError::MixedAlignmentFileFormats {
                first_path,
                first_format,
                other_path,
                other_format,
            } => {
                assert_eq!(first_path, bam);
                assert_eq!(first_format, "BAM");
                assert_eq!(other_path, cram);
                assert_eq!(other_format, "CRAM");
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    // --- M12: load_alignment_index direct unit tests ----------------
    //
    // The function is `pub`; integration tests only exercise one
    // configuration each. A regression that reversed the csi/bai
    // preference, dropped the bai fallback, or changed the
    // missing-index path would slip past the integration tests
    // and fail diffusely. These tests pin each branch in turn.

    #[test]
    fn load_alignment_index_returns_bam_csi_when_csi_present() {
        let dir = TempDir::new().expect("tempdir");
        let bam = dir.path().join("sample.bam");
        let csi = dir.path().join("sample.bam.csi");
        touch(&bam);
        // The on-disk .csi here is a placeholder one-byte file
        // (the loader will fail on parse). We assert that the
        // loader picked the .csi *path* via the LoadFailed error
        // — the typed surface tells us which file was tried.
        touch(&csi);

        let err = load_alignment_index(&bam).expect_err("parse fails on placeholder bytes");
        match err {
            AlignmentIndexError::LoadFailed {
                path,
                index_path,
                source: _,
            } => {
                assert_eq!(path, bam);
                assert_eq!(
                    index_path, csi,
                    "loader must have picked the .csi path, not the .bai"
                );
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn load_alignment_index_falls_back_to_bai_when_no_csi() {
        let dir = TempDir::new().expect("tempdir");
        let bam = dir.path().join("sample.bam");
        let bai = dir.path().join("sample.bam.bai");
        touch(&bam);
        // Only .bai present; the .csi-preferred policy must
        // still locate .bai as the fallback.
        touch(&bai);

        let err = load_alignment_index(&bam).expect_err("parse fails on placeholder bytes");
        match err {
            AlignmentIndexError::LoadFailed {
                path,
                index_path,
                source: _,
            } => {
                assert_eq!(path, bam);
                assert_eq!(index_path, bai, ".bai must be picked when .csi is absent");
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn load_alignment_index_errors_with_missing_alignment_index_when_no_bam_index() {
        let dir = TempDir::new().expect("tempdir");
        let bam = dir.path().join("sample.bam");
        touch(&bam);
        // No .csi, no .bai. Loader must surface a typed
        // MissingAlignmentIndex naming the .csi-canonical path
        // (the policy build path).
        let err = load_alignment_index(&bam).expect_err("must fail on missing index");
        match err {
            AlignmentIndexError::MissingAlignmentIndex {
                path,
                expected_index_path,
            } => {
                assert_eq!(path, bam);
                assert_eq!(
                    expected_index_path,
                    dir.path().join("sample.bam.csi"),
                    "expected_index_path must be the .csi-canonical build target"
                );
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn load_alignment_index_rejects_unsupported_extension() {
        let dir = TempDir::new().expect("tempdir");
        let input = dir.path().join("sample.sam");
        touch(&input);

        let err = load_alignment_index(&input).expect_err("must reject .sam");
        match err {
            AlignmentIndexError::UnsupportedExtension { path } => {
                assert_eq!(path, input);
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }
}
