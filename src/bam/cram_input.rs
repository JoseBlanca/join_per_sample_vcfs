//! CRAM record-stream decoders. The merge / filter / header-validation
//! surface lives in the sibling [`crate::bam::alignment_input`]
//! module; this file owns only the CRAM-specific noodles-cram seam:
//!
//! - [`OwnedCramRecords`] — linear-order owned iterator over a single
//!   CRAM's records, used by the no-index merge.
//! - [`OwnedIndexedCramRecords`] — per-`.crai`-driven owned iterator
//!   over a single CRAM's records on one target contig, used by the
//!   per-chromosome parallel driver.
//! - [`open_cram_record_stream`] /
//!   [`open_indexed_cram_record_stream`] — open helpers that
//!   `AlignmentMergedReader::{new, query}` call once per input CRAM.
//!
//! The shape `Box<dyn Iterator<Item = io::Result<RecordBuf>> + Send>`
//! (aliased as [`super::alignment_input::AlignmentRecordsIter`]) is
//! the format-agnostic seam the merge consumes; BAM support lives
//! in the sibling [`super::bam_input`] module that produces the
//! same shape.

use std::fs::File;
use std::io::{self, SeekFrom};
use std::path::Path;
use std::sync::Arc;

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
// Indexed (per-contig) owned CRAM record iterator
// ---------------------------------------------------------------------

/// Owned iterator that walks a CRAM via its `.crai` index, yielding
/// only records that align to `target_reference_sequence_id`. Mirror
/// of [`OwnedCramRecords`] but driven by the index instead of by
/// linear file order — the loop seeks to each container whose index
/// entry matches the target ref_id, decodes its records, and filters
/// the decoded set to records whose own `reference_sequence_id`
/// matches the target (a CRAM container may carry records from more
/// than one ref_seq, though slices typically do not).
///
/// The `Arc<crai::Index>` is shared read-only across rayon workers
/// (we hold our own clone here). The reader and the FASTA repository
/// are owned per worker — neither is `Send`-shareable in the
/// noodles surface used here.
///
/// Behaviour matches noodles' own
/// [`noodles_cram::io::reader::Query`](../../../noodles-cram-0.93.0/src/io/reader/query.rs)
/// stripped of the borrow-the-reader lifetime — we own the reader
/// instead so the iterator is `'static` and `Send`, satisfying the
/// `Box<dyn Iterator<Item=…> + Send>` shape the rest of the merger
/// expects.
struct OwnedIndexedCramRecords {
    reader: cram::io::Reader<File>,
    /// Shared parsed header. We use the caller's already-loaded
    /// `Arc<sam::Header>` rather than the per-worker reader's own
    /// `read_header()` result — both are equivalent (every CRAM in
    /// the input set has been cross-validated to carry the same
    /// header by the driver) but the `Arc` skips a re-parse per
    /// worker.
    header: Arc<sam::Header>,
    /// Shared parsed `.crai`. Held as `Arc` so many workers can
    /// share the parse cost; cloning the `Arc` per worker is one
    /// atomic.
    index: Arc<cram::crai::Index>,
    repository: fasta::Repository,
    target_reference_sequence_id: usize,
    /// Cursor into the index. We hold the index by `Arc`, so we walk
    /// it via this integer rather than a `slice::Iter` to keep the
    /// struct free of self-referential borrows.
    next_index_record: usize,
    pending: std::vec::IntoIter<sam::alignment::RecordBuf>,
    /// Same EOF-latch as [`OwnedCramRecords`] — once exhausted,
    /// short-circuit subsequent `next()` calls without re-entering
    /// noodles. The index-cursor exhaustion path is idempotent by
    /// construction (re-reading past `index.len()` simply returns
    /// `Ok(true)` again), but the read_container call on the LAST
    /// matching container is not — and the merger's
    /// `BufferedPeekable` will poll past `None` to confirm
    /// exhaustion. Cheap defensive latch; same shape as the one
    /// `OwnedCramRecords` carries.
    eof_latched: bool,
}

impl OwnedIndexedCramRecords {
    fn refill(&mut self) -> io::Result<bool> {
        loop {
            let record = match self.index.as_slice().get(self.next_index_record) {
                Some(r) => r,
                None => return Ok(true),
            };
            self.next_index_record += 1;

            if record.reference_sequence_id() != Some(self.target_reference_sequence_id) {
                continue;
            }

            self.reader.seek(SeekFrom::Start(record.offset()))?;
            // Allocate a fresh `Container` for each seek rather than
            // reusing one. noodles' own `Query::read_next_container`
            // does the same — the container caches per-instance
            // state (compression header, block buffers) that the
            // sequential reader keeps consistent across adjacent
            // containers but that becomes stale after a non-
            // sequential seek.
            let mut container = cram::io::reader::Container::default();
            let n = self.reader.read_container(&mut container)?;
            if n == 0 {
                return Ok(true);
            }

            let compression_header = container.compression_header()?;
            let mut all_records: Vec<sam::alignment::RecordBuf> = Vec::new();
            for slice_result in container.slices() {
                let slice = slice_result?;
                let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
                let cram_records = slice.records(
                    self.repository.clone(),
                    self.header.as_ref(),
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                )?;
                for record in &cram_records {
                    let rb = sam::alignment::RecordBuf::try_from_alignment_record(
                        self.header.as_ref(),
                        record,
                    )?;
                    if rb.reference_sequence_id() == Some(self.target_reference_sequence_id) {
                        all_records.push(rb);
                    }
                }
            }
            self.pending = all_records.into_iter();
            return Ok(false);
        }
    }
}

impl Iterator for OwnedIndexedCramRecords {
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

/// Open a CRAM, validate the file definition (rejects CRAM
/// major-version != 3), read the SAM header, and discard the
/// reader. Same validation as the first half of
/// [`open_cram_record_stream`]; used by callers that only need
/// the header (e.g. the cohort driver's per-input header harvest)
/// so they pay one open + one header-read per input rather than
/// duplicating the validation sequence inline.
///
/// Does not take a `fasta::Repository` — header-reading does not
/// consult the reference; slice decoding does.
pub(crate) fn read_cram_header_only(path: &Path) -> Result<sam::Header, AlignmentInputError> {
    let (_reader, header) = open_cram_reader_with_header(path, None)?;
    Ok(header)
}

/// Shared CRAM open + file-definition + header-read sequence.
/// One source of truth for the CRAM-version gate; both
/// `open_cram_record_stream` and `read_cram_header_only` call
/// through here.
///
/// `repository = None` is used by `read_cram_header_only`, which
/// does not need the reference (header reading does not consult
/// it). Slice decoding does, so `open_cram_record_stream` always
/// passes `Some`.
fn open_cram_reader_with_header(
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

/// Open `path` as a CRAM, advance past its file definition and SAM
/// header (both results discarded — the caller already holds the
/// validated header via `header`), and return an owned per-contig
/// indexed record-stream iterator.
///
/// `target_reference_sequence_id` must be valid against `header`'s
/// `@SQ` order (the caller is expected to have resolved the contig
/// name into an integer index against the canonical contig list
/// before invoking).
pub(super) fn open_indexed_cram_record_stream(
    path: &Path,
    header: Arc<sam::Header>,
    index: Arc<cram::crai::Index>,
    repository: fasta::Repository,
    target_reference_sequence_id: usize,
) -> Result<AlignmentRecordsIter, AlignmentInputError> {
    let mut noodles_cram_reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .build_from_path(path)
        .map_err(|source| AlignmentInputError::OpenFailed {
            path: path.to_path_buf(),
            source,
        })?;

    // Advance the reader past the CRAM file definition and file
    // header so subsequent seeks land cleanly. Both results are
    // discarded — we already have the validated header via the
    // caller's `Arc<sam::Header>`.
    noodles_cram_reader
        .read_file_definition()
        .map_err(|source| AlignmentInputError::OpenFailed {
            path: path.to_path_buf(),
            source,
        })?;
    noodles_cram_reader
        .read_file_header()
        .map_err(|source| AlignmentInputError::OpenFailed {
            path: path.to_path_buf(),
            source,
        })?;

    let records: AlignmentRecordsIter = Box::new(OwnedIndexedCramRecords {
        reader: noodles_cram_reader,
        header,
        index,
        repository,
        target_reference_sequence_id,
        next_index_record: 0,
        pending: Vec::new().into_iter(),
        eof_latched: false,
    });
    Ok(records)
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::per_sample::cram_files::build_cram_with_major_version;

    /// `read_cram_header_only` must fire the CRAM-version gate
    /// before returning a header. Regression for the bug closed by
    /// the M5 fix (the cohort driver's old inline header loader
    /// skipped this check, letting CRAM 4.x files pass header-load
    /// and only be rejected later by the per-worker `query`).
    #[test]
    fn read_cram_header_only_rejects_cram_4x() {
        let (_dir, cram_path) = build_cram_with_major_version(4, 0).expect("build forced cram 4.0");
        let err = read_cram_header_only(&cram_path)
            .expect_err("CRAM 4.0 must be rejected at header time");
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
