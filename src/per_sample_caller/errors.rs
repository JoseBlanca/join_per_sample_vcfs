//! Errors produced by the per-sample caller's CRAM input slice.
//!
//! See `ia/specs/per_sample_caller.md` §"Errors" for the catalogue of
//! failure modes this enum covers.

use std::path::PathBuf;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum CramInputError {
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

    #[error("malformed record in '{path}' (qname='{qname}'): {source}")]
    MalformedRecord {
        path: PathBuf,
        qname: String,
        #[source]
        source: std::io::Error,
    },

    #[error("I/O error on '{path}': {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
}
