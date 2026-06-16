//! CRAM record-stream decoders. The merge / filter / header-validation
//! surface lives in the sibling [`crate::bam::alignment_input`]
//! module; this file owns only the CRAM-specific noodles-cram seam:
//!
//! - [`OwnedCramRecords`] — linear-order owned iterator over a single
//!   CRAM's records, used by the no-index merge.
//! - [`open_cram_record_stream`] — the open helper `AlignmentMergedReader::new`
//!   calls once per input CRAM.
//! - [`open_cram_reader_with_header`] — opens (gating the CRAM version) +
//!   reads the header, used by `load_pileup_inputs` and the pooled
//!   [`crate::bam::segment_reader`].
//!
//! The shape `Box<dyn Iterator<Item = io::Result<RecordBuf>> + Send>`
//! (aliased as [`super::alignment_input::AlignmentRecordsIter`]) is
//! the format-agnostic seam the merge consumes; BAM support lives
//! in the sibling [`super::bam_input`] module that produces the
//! same shape.

use std::fs::File;
use std::io;
use std::path::Path;

use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::alignment_input::AlignmentRecordsIter;
use super::errors::AlignmentInputError;

// ---------------------------------------------------------------------
// Owned CRAM record iterator (linear order)
// ---------------------------------------------------------------------

/// Owns a CRAM `Reader<File>` plus its `sam::Header` and a clone of the
/// `fasta::Repository`, and yields decoded `RecordBuf`s. Reproduces
/// the work `noodles_cram::io::reader::Records` does, but as an owned
/// iterator (`'static + Send`) so it can be packed into the merge's
/// `Box<dyn Iterator>`.
struct OwnedCramRecords {
    reader: cram::io::Reader<File>,
    header: sam::Header,
    repository: fasta::Repository,
    container: cram::io::reader::Container,
    pending: std::vec::IntoIter<sam::alignment::RecordBuf>,
    /// `true` once `read_container` returns 0 (clean EOF). Subsequent
    /// `next()` calls return `None` without re-entering noodles —
    /// `noodles_cram 0.93`'s `Reader::read_container` is **not**
    /// idempotent at EOF: the second call returns
    /// `Err(InvalidData, TryFromIntError)`, the third returns
    /// `Err(UnexpectedEof, "failed to fill whole buffer")`. Without
    /// the latch, any consumer that polls past the iterator's first
    /// `None` (notably the k-way merge's `BufferedPeekable`, which
    /// peeks every per-input stream on each merge step) re-enters
    /// noodles past EOF and surfaces a spurious `MalformedRecord`
    /// error after every real record has already been emitted.
    /// Single-CRAM consumers never trip this because the outer
    /// iterator stops at the first `None` from the only inner
    /// iterator. Two-or-more-CRAM merges trip it on the trailing
    /// post-last-record poll the merge does to confirm the second
    /// stream is exhausted. Repro: `d1_repro_two_crams_…` in the
    /// `alignment_input` tests.
    eof_latched: bool,
}

impl OwnedCramRecords {
    fn refill(&mut self) -> io::Result<bool> {
        let n = self.reader.read_container(&mut self.container)?;
        if n == 0 {
            return Ok(true);
        }
        let compression_header = self.container.compression_header()?;
        let mut all_records: Vec<sam::alignment::RecordBuf> = Vec::new();
        for slice_result in self.container.slices() {
            let slice = slice_result?;
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
            let cram_records = slice.records(
                self.repository.clone(),
                &self.header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )?;
            for record in &cram_records {
                let rb =
                    sam::alignment::RecordBuf::try_from_alignment_record(&self.header, record)?;
                all_records.push(rb);
            }
        }
        self.pending = all_records.into_iter();
        Ok(false)
    }
}

impl Iterator for OwnedCramRecords {
    type Item = io::Result<sam::alignment::RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(rb) = self.pending.next() {
                return Some(Ok(rb));
            }
            if self.eof_latched {
                return None;
            }
            match self.refill() {
                Ok(true) => {
                    self.eof_latched = true;
                    return None;
                }
                Ok(false) => continue,
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

// ---------------------------------------------------------------------
// Open helpers
// ---------------------------------------------------------------------

/// Open `path` as a CRAM, validate its file definition (rejects
/// CRAM-major-version != 3), read the SAM header, and return the
/// header + an owned record-stream iterator ready to feed into
/// [`super::alignment_input::OpenAlignmentFile`].
///
/// The caller threads the returned `sam::Header` into the
/// format-agnostic cross-file validators in `alignment_input`;
/// nothing here knows about cohort-wide invariants.
pub(super) fn open_cram_record_stream(
    path: &Path,
    repository: fasta::Repository,
) -> Result<(sam::Header, AlignmentRecordsIter), AlignmentInputError> {
    let (noodles_cram_reader, noodles_sam_header) =
        open_cram_reader_with_header(path, Some(&repository))?;

    let records: AlignmentRecordsIter = Box::new(OwnedCramRecords {
        reader: noodles_cram_reader,
        header: noodles_sam_header.clone(),
        repository,
        container: cram::io::reader::Container::default(),
        pending: Vec::new().into_iter(),
        eof_latched: false,
    });
    Ok((noodles_sam_header, records))
}

/// Shared CRAM open + file-definition + header-read sequence.
/// One source of truth for the CRAM-version gate (rejects CRAM
/// major-version != 3), called by `open_cram_record_stream`.
///
/// `repository = None` reads the header without the reference
/// (header reading does not consult it); slice decoding does, so
/// `open_cram_record_stream` always passes `Some`.
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
