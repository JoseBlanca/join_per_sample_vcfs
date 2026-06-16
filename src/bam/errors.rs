//! Errors produced by the per-sample pileup's alignment-file input
//! slice (CRAM + BAM).
//!
//! See `ia/specs/per_sample_pileup.md` §"Errors" for the catalogue of
//! failure modes this enum covers.

use std::path::PathBuf;

use thiserror::Error;

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum AlignmentInputError {
    #[error("at least one alignment-file input is required")]
    NoInputs,

    #[error("failed to open alignment file '{path}': {source}")]
    OpenFailed {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("FASTA index (.fai) is missing for '{fasta_path}'")]
    MissingFastaIndex { fasta_path: PathBuf },

    #[error(
        "unsupported CRAM major version {major}.{minor} in '{path}' (only CRAM 3.x is supported)"
    )]
    UnsupportedCramVersion { path: PathBuf, major: u8, minor: u8 },

    #[error("alignment file '{path}' is not coordinate-sorted (found SO:'{sort_order}')")]
    NotCoordinateSorted { path: PathBuf, sort_order: String },

    #[error("@SQ list mismatch between '{reference_path}' and '{other_path}': {detail}")]
    ContigListMismatch {
        reference_path: PathBuf,
        other_path: PathBuf,
        detail: String,
    },

    #[error("FASTA '{fasta_path}' disagrees with alignment file '{alignment_file_path}': {detail}")]
    FastaContigMismatch {
        fasta_path: PathBuf,
        alignment_file_path: PathBuf,
        detail: String,
    },

    #[error("@RG '{read_group_id}' in '{path}' has no SM tag")]
    MissingSampleTag {
        path: PathBuf,
        read_group_id: String,
    },

    #[error(
        "multiple sample names across alignment files: '{path_a}' has SM:'{sm_a}', '{path_b}' has SM:'{sm_b}'"
    )]
    MultipleSampleNames {
        path_a: PathBuf,
        sm_a: String,
        path_b: PathBuf,
        sm_b: String,
    },

    #[error(
        "multiple sample names within '{path}': @RG '{rg_a}' has SM:'{sm_a}', @RG '{rg_b}' has SM:'{sm_b}'"
    )]
    MultipleSampleNamesInFile {
        path: PathBuf,
        rg_a: String,
        sm_a: String,
        rg_b: String,
        sm_b: String,
    },

    #[error(
        "out-of-order read in '{path}': QNAME '{qname}' at \
         (ref_id={this_ref_id}, pos={this_pos}) regresses from \
         (ref_id={prev_ref_id}, pos={prev_pos})"
    )]
    OutOfOrderRead {
        path: PathBuf,
        qname: String,
        prev_ref_id: usize,
        prev_pos: u64,
        this_ref_id: usize,
        this_pos: u64,
    },

    #[error(
        "duplicate read across alignment files: QNAME '{qname}' at ref_id={ref_id} pos={pos} appears in both '{path_a}' and '{path_b}'"
    )]
    DuplicateReadAcrossFiles {
        qname: String,
        path_a: PathBuf,
        path_b: PathBuf,
        ref_id: usize,
        pos: u64,
    },

    #[error("malformed @SQ M5 in '{path}' for contig '{contig}': {detail}")]
    MalformedMd5 {
        path: PathBuf,
        contig: String,
        detail: String,
    },

    #[error(
        "malformed record in '{path}'{}: {source}",
        match qname {
            Some(q) => format!(" (qname='{q}')"),
            None => String::new(),
        }
    )]
    MalformedRecord {
        path: PathBuf,
        qname: Option<String>,
        #[source]
        source: std::io::Error,
    },

    #[error("I/O error on '{path}': {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// `AlignmentMergedReader::query` was asked for a contig name that
    /// the canonical contig list does not carry. Programmer error
    /// at the call site (the canonical list comes from a prior
    /// validation pass; the caller is expected to query contigs
    /// from that same list), but surfaced as a typed error rather
    /// than a panic to keep the boundary clean.
    #[error(
        "contig '{contig}' is not in the canonical contig list \
         (cohort has {known_contigs} contigs total)"
    )]
    ContigNotInList {
        contig: String,
        known_contigs: usize,
    },

    /// `AlignmentMergedReader::query` was given a `headers` or `indexes`
    /// slice whose length does not match `inputs.len()`.
    #[error(
        "per-input handle count mismatch: {inputs} inputs, {headers} headers, {indexes} indexes"
    )]
    PerInputHandleCountMismatch {
        inputs: usize,
        headers: usize,
        indexes: usize,
    },

    /// Input file's extension is neither `.cram` nor `.bam`.
    /// `AlignmentMergedReader::{new,query}` need to dispatch on
    /// the source format and we do not sniff magic bytes.
    #[error(
        "input alignment file '{path}' has an unsupported extension \
         (expected .cram or .bam)"
    )]
    UnsupportedExtension { path: PathBuf },

    /// Inputs to a single invocation mixed `.cram` and `.bam`
    /// files. One format per invocation only — same policy the
    /// pre-flight enforces; surfaced here too because not every
    /// caller (e.g. `pileup`) runs through the pre-flight.
    #[error(
        "mixed alignment-file formats are not supported in one \
         invocation: '{first_path}' is {first_format}, '{other_path}' is {other_format}"
    )]
    MixedAlignmentFileFormats {
        first_path: PathBuf,
        first_format: &'static str,
        other_path: PathBuf,
        other_format: &'static str,
    },

    /// `AlignmentMergedReader::query` was handed an `AlignmentIndex`
    /// variant whose on-disk format does not match the input
    /// file's extension (e.g. a `.cram` paired with a
    /// `BamCsi`-typed index, or vice versa). Indicates the
    /// driver's index loader and file-classification logic
    /// disagree — surfaced as a typed error rather than a panic
    /// to keep the boundary clean.
    #[error(
        "alignment-index format does not match input file '{path}': \
         input is {file_format}, index is {index_format}"
    )]
    AlignmentIndexFormatMismatch {
        path: PathBuf,
        file_format: &'static str,
        index_format: &'static str,
    },

    /// Index pre-flight or per-input index load failed (missing index
    /// with no build opt-in, a failed build, etc.). Carries the typed
    /// [`AlignmentIndexError`] so the CLI layer can name the
    /// responsible flag.
    #[error("alignment index: {0}")]
    AlignmentIndex(#[from] AlignmentIndexError),

    /// [`crate::bam::segment_reader::AlignmentFile::from_input`] was
    /// asked to build a CRAM segment reader without a FASTA reference
    /// repository. CRAM slice decoding needs the reference; the caller
    /// must pass one. Programmer error at the call site, surfaced as a
    /// typed error rather than a panic.
    #[error("CRAM input '{path}' requires a FASTA reference repository but none was provided")]
    MissingCramReference { path: PathBuf },

    /// [`crate::bam::segment_reader::AlignmentFile::get_reads_from_segment`]
    /// was given a segment that is not a valid 1-based inclusive range
    /// (`start == 0`, or `start > end`).
    #[error("invalid segment on '{chrom}': start={start}, end={end} (require 1 <= start <= end)")]
    InvalidSegment { chrom: String, start: u32, end: u32 },
}

/// Errors raised by [`crate::bam::index_preflight::preflight_alignment_indexes`].
///
/// The variants here are deliberately CLI-vocabulary-free (no `--flag`
/// names). The CLI layer wraps them into user-facing error variants
/// that name the responsible flag — keeps `crate::bam` independent of
/// any single caller's argument surface.
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum AlignmentIndexError {
    /// An input alignment file has no index next to it.
    #[error(
        "input alignment file '{path}' has no index \
         (looked for '{expected_index_path}')"
    )]
    MissingAlignmentIndex {
        path: PathBuf,
        expected_index_path: PathBuf,
    },

    /// Index construction or write failed (commonly a read-only
    /// directory next to the source file).
    #[error("failed to build alignment index for '{path}': {source}")]
    BuildFailed {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// Loading a previously-built alignment index from disk
    /// failed (corrupted index, permission denied, …). Surfaced
    /// for both CRAM (`.crai`) and BAM (`.csi`/`.bai`) inputs.
    #[error("failed to load alignment index '{index_path}' for '{path}': {source}")]
    LoadFailed {
        path: PathBuf,
        index_path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// Input's extension is neither `.cram` nor `.bam`. Pre-flight
    /// does not sniff file magic — an explicit extension is required.
    #[error(
        "input alignment file '{path}' has an unsupported extension \
         (expected .cram or .bam)"
    )]
    UnsupportedExtension { path: PathBuf },

    /// Inputs to a single invocation mixed `.cram` and `.bam`
    /// files. One format per invocation only; mixed-format support
    /// is an explicit non-goal of the BAM-input plan.
    #[error(
        "mixed alignment-file formats are not supported in one \
         invocation: '{first_path}' is {first_format}, '{other_path}' is {other_format}"
    )]
    MixedAlignmentFileFormats {
        first_path: PathBuf,
        first_format: &'static str,
        other_path: PathBuf,
        other_format: &'static str,
    },
}
