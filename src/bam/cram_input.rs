//! CRAM open helpers. The header-validation surface lives in the
//! sibling [`crate::bam::alignment_input`] module; this file owns only
//! the CRAM-specific noodles-cram open seam:
//!
//! - [`open_cram_reader_with_header`] — opens (gating the CRAM version) +
//!   reads the header, used by `load_pileup_inputs` and the pooled
//!   [`crate::bam::segment_reader`].

use std::fs::File;
use std::path::Path;

use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::errors::AlignmentInputError;

// ---------------------------------------------------------------------
// Open helpers
// ---------------------------------------------------------------------

/// Shared CRAM open + file-definition + header-read sequence.
/// One source of truth for the CRAM-version gate (rejects CRAM
/// major-version != 3).
///
/// `repository = None` reads the header without the reference
/// (header reading does not consult it); slice decoding does, so a
/// record-decoding caller passes `Some`.
///
/// Also used by the pooled [`crate::bam::segment_reader`] to open a
/// fresh seekable reader positioned past the file definition + header
/// (the returned header is discarded — the caller already holds the
/// validated one); the version gate still runs on every open.
pub(super) fn open_cram_reader_with_header(
    path: &Path,
    repository: Option<&fasta::Repository>,
) -> Result<(cram::io::Reader<File>, sam::Header), AlignmentInputError> {
    let mut builder = cram::io::reader::Builder::default();
    if let Some(repo) = repository {
        builder = builder.set_reference_sequence_repository(repo.clone());
    }
    let mut noodles_cram_reader =
        builder
            .build_from_path(path)
            .map_err(|source| AlignmentInputError::OpenFailed {
                path: path.to_path_buf(),
                source,
            })?;

    // Read file definition first to gate on CRAM major version
    // before any container is decoded.
    let file_definition = noodles_cram_reader
        .read_file_definition()
        .map_err(|source| AlignmentInputError::OpenFailed {
            path: path.to_path_buf(),
            source,
        })?;
    let version = file_definition.version();
    if version.major() != 3 {
        return Err(AlignmentInputError::UnsupportedCramVersion {
            path: path.to_path_buf(),
            major: version.major(),
            minor: version.minor(),
        });
    }

    let noodles_sam_header = noodles_cram_reader.read_file_header().map_err(|source| {
        AlignmentInputError::OpenFailed {
            path: path.to_path_buf(),
            source,
        }
    })?;
    Ok((noodles_cram_reader, noodles_sam_header))
}

// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::per_sample::cram_files::build_cram_with_major_version;

    /// The shared CRAM open path must fire the CRAM-version gate
    /// (rejecting major-version != 3) before returning a header, so
    /// CRAM 4.x is refused at open time rather than mid-decode.
    #[test]
    fn open_cram_reader_rejects_cram_4x() {
        let (_dir, cram_path) = build_cram_with_major_version(4, 0).expect("build forced cram 4.0");
        let err = open_cram_reader_with_header(&cram_path, None)
            .map(|_| ())
            .expect_err("CRAM 4.0 must be rejected at open time");
        match err {
            AlignmentInputError::UnsupportedCramVersion { major, minor, path } => {
                assert_eq!(major, 4);
                assert_eq!(minor, 0);
                assert_eq!(path, cram_path);
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }
}
