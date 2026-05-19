//! Batch-assignment TSV reader for the `estimate-contamination`
//! subcommand.
//!
//! The TSV maps cohort sample names to user-chosen batch identifiers.
//! It is the only way to express multi-batch cohorts to the cohort CLI;
//! the contamination estimator's batch concept (per-batch `q_b`) is
//! invisible to the rest of the pipeline.
//!
//! ## Format
//!
//! Two columns exactly, tab-separated, header row required:
//!
//! ```tsv
//! sample	batch
//! NA12878	lane_3
//! NA12891	lane_3
//! NA12892	lane_7
//! ```
//!
//! - Header row literal: `sample\tbatch` (case-sensitive).
//! - Body rows: `<sample>\t<batch>` — sample and batch are arbitrary
//!   non-empty strings that may not contain `\t` or `\n`.
//! - Sample identifiers must be unique within the file (duplicate
//!   entries surface a hard error rather than a silent last-wins).
//! - Comments / blank lines are not supported. Trailing whitespace is
//!   preserved verbatim — the consumer is expected to feed clean
//!   files.
//! - Sample names that aren't in the cohort being processed are
//!   silently ignored downstream; missing-from-TSV samples fall back
//!   to the default batch `"all_samples"` via [`BatchAssignment::batch_for`].
//!
//! The expected consumer is the contamination subcommand; sample-name
//! reconciliation against the `.psp` inputs lives there, not here.

// The TSV example in the module docs contains literal tab characters —
// faithful to the on-disk format. Clippy's `tabs_in_doc_comments` lint
// would prefer 4-space replacements, but that would misrepresent the
// format. Allowed module-wide.
#![allow(clippy::tabs_in_doc_comments)]

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use thiserror::Error;

/// Default batch identifier assigned to samples not listed in a
/// `--batch-assignment` TSV (or when the TSV is omitted entirely).
/// Documented in the cohort CLI plan §"Batch assignment".
pub const DEFAULT_BATCH_ID: &str = "all_samples";

/// Sample → batch lookup loaded from `--batch-assignment`.
///
/// Construct with [`Self::from_tsv`] for the user-supplied case or
/// [`Self::empty`] when no TSV was passed. Look up a sample's batch
/// with [`Self::batch_for`]; samples not present fall back to
/// [`DEFAULT_BATCH_ID`].
#[derive(Debug, Clone, Default)]
pub struct BatchAssignment {
    by_sample: HashMap<String, String>,
}

impl BatchAssignment {
    /// Empty mapping — every sample falls back to [`DEFAULT_BATCH_ID`].
    /// Used when `--batch-assignment` is omitted.
    pub fn empty() -> Self {
        Self::default()
    }

    /// Read and validate a TSV file at `path`. See the module docs for
    /// the exact format.
    pub fn from_tsv(path: &Path) -> Result<Self, BatchAssignmentError> {
        let contents = fs::read_to_string(path).map_err(|source| BatchAssignmentError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        Self::from_str_with_path(&contents, path)
    }

    /// Same as [`Self::from_tsv`] but takes the file contents as a
    /// string slice. Lets tests cover edge cases without writing to
    /// disk. `display_path` is reported in error messages.
    fn from_str_with_path(
        contents: &str,
        display_path: &Path,
    ) -> Result<Self, BatchAssignmentError> {
        let mut lines = contents.lines().enumerate();
        let (_, header) = lines
            .next()
            .ok_or_else(|| BatchAssignmentError::EmptyFile {
                path: display_path.to_path_buf(),
            })?;
        if header.trim_end() != "sample\tbatch" {
            return Err(BatchAssignmentError::MalformedHeader {
                path: display_path.to_path_buf(),
                got: header.to_string(),
            });
        }
        let mut by_sample: HashMap<String, String> = HashMap::new();
        for (idx, line) in lines {
            let line_no = idx + 1; // human-readable 1-based
            if line.is_empty() {
                return Err(BatchAssignmentError::EmptyLine {
                    path: display_path.to_path_buf(),
                    line_no,
                });
            }
            let mut parts = line.split('\t');
            let sample = parts
                .next()
                .ok_or_else(|| BatchAssignmentError::MalformedRow {
                    path: display_path.to_path_buf(),
                    line_no,
                    got: line.to_string(),
                })?;
            let batch = parts
                .next()
                .ok_or_else(|| BatchAssignmentError::MalformedRow {
                    path: display_path.to_path_buf(),
                    line_no,
                    got: line.to_string(),
                })?;
            if parts.next().is_some() {
                return Err(BatchAssignmentError::MalformedRow {
                    path: display_path.to_path_buf(),
                    line_no,
                    got: line.to_string(),
                });
            }
            if sample.is_empty() {
                return Err(BatchAssignmentError::EmptySample {
                    path: display_path.to_path_buf(),
                    line_no,
                });
            }
            if batch.is_empty() {
                return Err(BatchAssignmentError::EmptyBatch {
                    path: display_path.to_path_buf(),
                    line_no,
                    sample: sample.to_string(),
                });
            }
            if let Some(existing) = by_sample.insert(sample.to_string(), batch.to_string()) {
                return Err(BatchAssignmentError::DuplicateSample {
                    path: display_path.to_path_buf(),
                    line_no,
                    sample: sample.to_string(),
                    first_batch: existing,
                    second_batch: batch.to_string(),
                });
            }
        }
        Ok(Self { by_sample })
    }

    /// Look up the batch a sample belongs to. Falls back to
    /// [`DEFAULT_BATCH_ID`] when the sample is not present in the
    /// loaded TSV.
    pub fn batch_for<'a>(&'a self, sample: &str) -> &'a str {
        self.by_sample
            .get(sample)
            .map(String::as_str)
            .unwrap_or(DEFAULT_BATCH_ID)
    }

    /// Number of explicit (sample, batch) entries loaded. `0` for
    /// [`Self::empty`].
    pub fn len(&self) -> usize {
        self.by_sample.len()
    }

    /// `true` iff no explicit entries were loaded.
    pub fn is_empty(&self) -> bool {
        self.by_sample.is_empty()
    }
}

/// Errors surfaced when reading or validating a batch-assignment TSV.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum BatchAssignmentError {
    /// `read_to_string` failed.
    #[error("batch assignment {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// The file was zero-length.
    #[error("batch assignment {path}: file is empty (expected at least a header row)")]
    EmptyFile { path: PathBuf },

    /// The header row was missing or didn't match `sample\tbatch`.
    #[error("batch assignment {path}: malformed header — expected `sample\\tbatch`, got `{got}`")]
    MalformedHeader { path: PathBuf, got: String },

    /// A body row did not have exactly two tab-separated fields.
    #[error(
        "batch assignment {path}:{line_no}: malformed row \
         (expected exactly `sample\\tbatch`), got `{got}`"
    )]
    MalformedRow {
        path: PathBuf,
        line_no: usize,
        got: String,
    },

    /// A body row was a literal empty line.
    #[error("batch assignment {path}:{line_no}: empty line (blank rows not supported)")]
    EmptyLine { path: PathBuf, line_no: usize },

    /// A body row's `sample` field was empty.
    #[error("batch assignment {path}:{line_no}: sample field is empty")]
    EmptySample { path: PathBuf, line_no: usize },

    /// A body row's `batch` field was empty.
    #[error("batch assignment {path}:{line_no}: batch field for sample `{sample}` is empty")]
    EmptyBatch {
        path: PathBuf,
        line_no: usize,
        sample: String,
    },

    /// Two rows assigned different batches to the same sample.
    #[error(
        "batch assignment {path}:{line_no}: sample `{sample}` appears twice \
         (first as `{first_batch}`, then as `{second_batch}`)"
    )]
    DuplicateSample {
        path: PathBuf,
        line_no: usize,
        sample: String,
        first_batch: String,
        second_batch: String,
    },
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(contents: &str) -> Result<BatchAssignment, BatchAssignmentError> {
        BatchAssignment::from_str_with_path(contents, Path::new("test.tsv"))
    }

    #[test]
    fn empty_falls_back_to_default_batch() {
        let m = BatchAssignment::empty();
        assert_eq!(m.batch_for("any_sample"), DEFAULT_BATCH_ID);
        assert_eq!(m.len(), 0);
        assert!(m.is_empty());
    }

    #[test]
    fn parses_three_sample_two_batch_file() {
        let m =
            parse("sample\tbatch\nNA12878\tlane_3\nNA12891\tlane_3\nNA12892\tlane_7\n").unwrap();
        assert_eq!(m.batch_for("NA12878"), "lane_3");
        assert_eq!(m.batch_for("NA12891"), "lane_3");
        assert_eq!(m.batch_for("NA12892"), "lane_7");
        assert_eq!(m.batch_for("NA00000"), DEFAULT_BATCH_ID);
        assert_eq!(m.len(), 3);
    }

    #[test]
    fn header_with_trailing_newline_is_ok() {
        // File without trailing newline.
        let m = parse("sample\tbatch\nNA12878\tlane_3").unwrap();
        assert_eq!(m.batch_for("NA12878"), "lane_3");
    }

    /// Mi7 verification: `&str::lines()` recognises `\r\n` as a single
    /// line terminator and strips both characters, so a CRLF-authored
    /// (Windows-style) TSV parses identically to LF-only. This test
    /// locks the contract so a future change to a custom splitter
    /// can't silently store `"NA12878\r"` / `"lane_3\r"` as keys.
    #[test]
    fn body_rows_with_trailing_carriage_return_preserve_sample_name() {
        let m = parse("sample\tbatch\r\nNA12878\tlane_3\r\nNA12891\tlane_7\r\n").unwrap();
        assert_eq!(m.batch_for("NA12878"), "lane_3");
        assert_eq!(m.batch_for("NA12891"), "lane_7");
        // And the `\r`-suffixed forms must not exist as keys.
        assert_eq!(m.batch_for("NA12878\r"), DEFAULT_BATCH_ID);
    }

    #[test]
    fn header_with_crlf_is_ok() {
        // CRLF-terminated header — `trim_end()` strips the `\r` so
        // header-equality still holds.
        let m = parse("sample\tbatch\r\nNA12878\tlane_3\r\n").unwrap();
        assert_eq!(m.batch_for("NA12878"), "lane_3");
    }

    #[test]
    fn empty_file_errors() {
        let err = parse("").unwrap_err();
        assert!(matches!(err, BatchAssignmentError::EmptyFile { .. }));
    }

    #[test]
    fn missing_header_errors() {
        let err = parse("NA12878\tlane_3\n").unwrap_err();
        assert!(matches!(err, BatchAssignmentError::MalformedHeader { .. }));
    }

    #[test]
    fn header_with_extra_column_errors() {
        let err = parse("sample\tbatch\textra\nNA12878\tlane_3\n").unwrap_err();
        assert!(matches!(err, BatchAssignmentError::MalformedHeader { .. }));
    }

    #[test]
    fn case_sensitive_header() {
        // `Sample` ≠ `sample`.
        let err = parse("Sample\tBatch\nNA12878\tlane_3\n").unwrap_err();
        assert!(matches!(err, BatchAssignmentError::MalformedHeader { .. }));
    }

    #[test]
    fn empty_batch_field_errors() {
        let err = parse("sample\tbatch\nNA12878\t\n").unwrap_err();
        match err {
            BatchAssignmentError::EmptyBatch { sample, .. } => assert_eq!(sample, "NA12878"),
            other => panic!("expected EmptyBatch, got {other:?}"),
        }
    }

    #[test]
    fn empty_sample_field_errors() {
        let err = parse("sample\tbatch\n\tlane_3\n").unwrap_err();
        assert!(matches!(err, BatchAssignmentError::EmptySample { .. }));
    }

    #[test]
    fn empty_line_in_body_errors() {
        let err = parse("sample\tbatch\nNA12878\tlane_3\n\nNA12891\tlane_7\n").unwrap_err();
        assert!(matches!(err, BatchAssignmentError::EmptyLine { .. }));
    }

    #[test]
    fn extra_columns_in_body_error() {
        let err = parse("sample\tbatch\nNA12878\tlane_3\textra\n").unwrap_err();
        assert!(matches!(err, BatchAssignmentError::MalformedRow { .. }));
    }

    #[test]
    fn missing_tab_in_body_errors() {
        let err = parse("sample\tbatch\nNA12878 lane_3\n").unwrap_err();
        // single field after split('\t') → second `parts.next()` is None.
        assert!(matches!(err, BatchAssignmentError::MalformedRow { .. }));
    }

    #[test]
    fn duplicate_sample_with_different_batch_errors() {
        let err = parse("sample\tbatch\nNA12878\tlane_3\nNA12878\tlane_7\n").unwrap_err();
        match err {
            BatchAssignmentError::DuplicateSample {
                sample,
                first_batch,
                second_batch,
                ..
            } => {
                assert_eq!(sample, "NA12878");
                assert_eq!(first_batch, "lane_3");
                assert_eq!(second_batch, "lane_7");
            }
            other => panic!("expected DuplicateSample, got {other:?}"),
        }
    }

    #[test]
    fn duplicate_sample_with_same_batch_also_errors() {
        // Same-batch duplicate is also rejected — symmetric "first
        // hit wins" rule would silently mask bigger issues if a user
        // appends an already-present row by mistake.
        let err = parse("sample\tbatch\nNA12878\tlane_3\nNA12878\tlane_3\n").unwrap_err();
        assert!(matches!(err, BatchAssignmentError::DuplicateSample { .. }));
    }

    #[test]
    fn from_tsv_io_error_propagates() {
        let err = BatchAssignment::from_tsv(Path::new("/no/such/path.tsv")).unwrap_err();
        assert!(matches!(err, BatchAssignmentError::Io { .. }));
    }
}
