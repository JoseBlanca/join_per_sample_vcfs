//! BAM open helpers. Sibling of [`crate::bam::cram_input`] — the
//! header-validation surface lives in [`crate::bam::alignment_input`];
//! this file owns only the BAM-specific noodles-bam open seam.
//!
//! - [`open_bam_reader_with_header`] — opens + reads the header, used by
//!   `load_pileup_inputs` and the pooled [`crate::bam::segment_reader`].
//! - [`query_interval`] / [`BamIndex`] — the index-query shapes the pooled
//!   segment reader reuses.

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_core::region::Interval;
use noodles_sam as sam;

use super::alignment_input::ContigInterval;
use super::errors::AlignmentInputError;

// ---------------------------------------------------------------------
// BAM-side index payload
// ---------------------------------------------------------------------

/// BAM index payload handed to the indexed open helper. Either index
/// format works at noodles' [`BinningIndex`] layer; the caller (the
/// `index_preflight` layer, downstream) hands us whichever it
/// loaded.
#[derive(Clone)]
#[non_exhaustive]
pub(super) enum BamIndex {
    /// `.bai` — the legacy BAM index. 16 kbp bins; cannot address
    /// contigs above 512 Mbp.
    Bai(Arc<bam::bai::Index>),
    /// `.csi` — coordinate-sorted index. Preferred on read and
    /// always produced on build, because it has no contig-length
    /// limit.
    Csi(Arc<noodles_csi::Index>),
}

// ---------------------------------------------------------------------
// Open helpers
// ---------------------------------------------------------------------

/// Shared BAM open + header-read sequence so the magic-byte check
/// and the error wrapping stay in one place. Also used by the pooled
/// [`crate::bam::segment_reader`] to open a fresh seekable reader
/// (discarding the re-parsed header — the caller already holds the
/// validated one), so the reader is positioned past the header and
/// any chunk seek lands at a record boundary.
pub(super) fn open_bam_reader_with_header(
    path: &Path,
) -> Result<(bam::io::Reader<bgzf::io::Reader<File>>, sam::Header), AlignmentInputError> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .map_err(|source| AlignmentInputError::OpenFailed {
            path: path.to_path_buf(),
            source,
        })?;
    let header = reader
        .read_header()
        .map_err(|source| AlignmentInputError::OpenFailed {
            path: path.to_path_buf(),
            source,
        })?;
    Ok((reader, header))
}

/// The `noodles` query interval for a `region`: the 1-based inclusive
/// range when `Some`, or the whole contig (unbounded) when `None`.
pub(super) fn query_interval(region: Option<ContigInterval>) -> Interval {
    match region {
        Some(iv) => {
            // `iv.start`/`iv.end` are >= 1 by construction, so the
            // `Position::new` calls cannot fail.
            let start = Position::new(iv.start as usize).expect("region start >= 1");
            let end = Position::new(iv.end as usize).expect("region end >= 1");
            Interval::from(start..=end)
        }
        None => Interval::from(..),
    }
}
