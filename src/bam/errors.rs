//! Errors produced by the per-sample pileup's alignment-file input
//! slice (CRAM today; BAM support tracked in
//! `doc/devel/implementation_plans/bam_input_support.md`).
//!
//! See `ia/specs/per_sample_pileup.md` §"Errors" for the catalogue of
//! failure modes this enum covers.

use std::path::PathBuf;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum AlignmentInputError {
    #[error("at least one CRAM input is required")]
    NoInputs,

    #[error("failed to open CRAM '{path}': {source}")]
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

    #[error("CRAM '{path}' is not coordinate-sorted (found SO:'{sort_order}')")]
    NotCoordinateSorted { path: PathBuf, sort_order: String },

    #[error("@SQ list mismatch between '{reference_path}' and '{other_path}': {detail}")]
    ContigListMismatch {
        reference_path: PathBuf,
        other_path: PathBuf,
        detail: String,
    },

    #[error("FASTA '{fasta_path}' disagrees with CRAM '{cram_path}': {detail}")]
    FastaContigMismatch {
        fasta_path: PathBuf,
        cram_path: PathBuf,
        detail: String,
    },

    #[error("@RG '{read_group_id}' in '{path}' has no SM tag")]
    MissingSampleTag {
        path: PathBuf,
        read_group_id: String,
    },

    #[error(
        "multiple sample names across CRAMs: '{path_a}' has SM:'{sm_a}', '{path_b}' has SM:'{sm_b}'"
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
        "duplicate read across CRAMs: QNAME '{qname}' at ref_id={ref_id} pos={pos} appears in both '{path_a}' and '{path_b}'"
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
    /// slice whose length does not match `crams.len()`.
    #[error("per-input handle count mismatch: {crams} crams, {headers} headers, {indexes} indexes")]
    PerInputHandleCountMismatch {
        crams: usize,
        headers: usize,
        indexes: usize,
    },
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

    /// Input's extension is neither `.cram` nor `.bam`. Pre-flight
    /// does not sniff file magic — an explicit extension is required.
    #[error(
        "input alignment file '{path}' has an unsupported extension \
         (expected .cram or .bam)"
    )]
    UnsupportedExtension { path: PathBuf },
}
