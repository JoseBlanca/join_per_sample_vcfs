//! Error type for the cohort VCF writer.

use std::io;
use std::path::PathBuf;

use thiserror::Error;

/// Errors surfaced by [`CohortVcfWriter`](super::CohortVcfWriter) and
/// its helpers.
///
/// Variants name the *operation that failed*, not the underlying
/// mechanism — the source chain (`std::error::Error::source`) carries
/// the original cause. Each `#[source]` field is the typed root cause;
/// `Display` messages describe what the writer was trying to do and on
/// which input, so chain-walking printers (`anyhow`, `eyre`,
/// `tracing::error!(error = …)`) render a clean "operation: …; caused
/// by: …" trace without doubling the cause text.
///
/// The enum is `#[non_exhaustive]` so future variants (additional
/// per-operation I/O cases, new structural-mismatch reasons) land as
/// non-breaking minor versions once the crate is published.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum VcfWriteError {
    // ---------- sink I/O (per-operation, source = io::Error) ----------
    /// `File::create` of the tmp output failed.
    #[error("failed to create tmp output {tmp_path}")]
    CreateTmp {
        tmp_path: PathBuf,
        #[source]
        source: io::Error,
    },

    /// Writing the VCF header through the line-oriented writer failed.
    #[error("failed to write VCF header to {tmp_path}")]
    WriteHeader {
        tmp_path: PathBuf,
        #[source]
        source: io::Error,
    },

    /// Writing one VCF data record through the line-oriented writer
    /// failed. `chrom_id` / `pos` identify the record under
    /// construction.
    #[error("failed to write VCF record at {chrom_id}:{pos}")]
    WriteRecord {
        chrom_id: u32,
        pos: u32,
        #[source]
        source: io::Error,
    },

    /// The bgzf sink failed to emit its EOF block or flush its tail.
    #[error("failed to flush bgzf sink for {tmp_path}")]
    FinishBgzf {
        tmp_path: PathBuf,
        #[source]
        source: io::Error,
    },

    /// `File::sync_all` on the tmp output file failed.
    #[error("failed to fsync {tmp_path}")]
    FsyncFile {
        tmp_path: PathBuf,
        #[source]
        source: io::Error,
    },

    /// `File::sync_all` on the parent directory failed; without this,
    /// the rename is not durable across a crash.
    #[error("failed to fsync parent directory of {final_path}")]
    FsyncDir {
        final_path: PathBuf,
        #[source]
        source: io::Error,
    },

    /// `fs::rename` from the tmp path to the final path failed.
    #[error("failed to rename {tmp_path} -> {final_path}")]
    Rename {
        tmp_path: PathBuf,
        final_path: PathBuf,
        #[source]
        source: io::Error,
    },

    // ---------- noodles-vcf encode boundary ----------
    /// A noodles encode call failed. `operation` tags the project-side
    /// site (e.g. `"##source key parse"`, `"allele bytes UTF-8"`); the
    /// typed cause is boxed to keep `noodles_*` error types out of
    /// this enum's public surface — the writer is published as a
    /// library and we do not want noodles' fast semver cadence to
    /// force our own major bumps every release. Downstream code that
    /// genuinely needs the underlying error can `source.downcast_ref`
    /// on the chain.
    #[error("VCF encode failed during {operation}")]
    Encode {
        operation: &'static str,
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    // ---------- metadata validation (no underlying source) ----------
    /// `CohortMetadata` failed pre-construction validation.
    #[error("invalid metadata: {0}")]
    InvalidMetadata(String),

    // ---------- per-record contract checks ----------
    /// Records arrived out of order. `(chrom_id, pos)` of the offending
    /// record and the prior accepted record are both recorded.
    #[error(
        "record out of order: record at {chrom_id}:{pos} \
         is not after previous {prev_chrom_id}:{prev_pos}"
    )]
    RecordOutOfOrder {
        chrom_id: u32,
        pos: u32,
        prev_chrom_id: u32,
        prev_pos: u32,
    },

    /// A sample's `best_genotype` value indexes past the canonical
    /// `genotype_order` table for its `(ploidy, n_alleles)`.
    #[error(
        "record at {chrom_id}:{pos}: sample {sample_idx} \
         best_genotype index {got} out of bounds for n_genotypes {n_genotypes}"
    )]
    GenotypeIndexOutOfBounds {
        chrom_id: u32,
        pos: u32,
        sample_idx: usize,
        got: usize,
        n_genotypes: usize,
    },

    /// A decoded genotype refers to an allele index past the record's
    /// allele set. Distinct from
    /// [`GenotypeIndexOutOfBounds`](Self::GenotypeIndexOutOfBounds):
    /// the genotype-order table lookup itself succeeded, but a
    /// referenced allele does not exist on the record.
    #[error(
        "record at {chrom_id}:{pos}: sample {sample_idx} genotype decodes \
         allele index {allele_idx} but record carries only {n_alleles} alleles"
    )]
    AlleleIndexOutOfBounds {
        chrom_id: u32,
        pos: u32,
        sample_idx: usize,
        allele_idx: u8,
        n_alleles: usize,
    },

    /// `RecordLocus.chrom_id` does not index a contig in the cohort
    /// metadata's contig table.
    #[error(
        "record at {chrom_id}:{pos}: chrom_id is out of bounds for the \
         {n_contigs} contig(s) declared in the cohort metadata"
    )]
    UnknownChromId {
        chrom_id: u32,
        pos: u32,
        n_contigs: usize,
    },

    /// The record's `n_samples` does not match the cohort metadata's
    /// sample count.
    #[error(
        "record at {chrom_id}:{pos}: cohort metadata names {expected_samples} \
         samples but the posterior arrays carry {got_samples}"
    )]
    SampleCountMismatch {
        chrom_id: u32,
        pos: u32,
        expected_samples: usize,
        got_samples: usize,
    },

    /// A per-record vector (`best_genotype`, `gq_phred`,
    /// `allele_frequencies`, `scalars`, `posteriors`,
    /// `chain_anchor_flags`, …) does not have the length expected
    /// from the record's declared shape.
    #[error("record at {chrom_id}:{pos}: field `{field}` has length {actual}, expected {expected}")]
    InconsistentRecord {
        chrom_id: u32,
        pos: u32,
        field: &'static str,
        expected: usize,
        actual: usize,
    },

    /// A depth value (per-sample DP or cohort total DP) overflows
    /// `i32`, which is the VCF spec's per-field integer range.
    /// `sample_idx` is `Some(s)` for per-sample overflow and `None`
    /// for the cohort total.
    #[error("record at {chrom_id}:{pos}: depth {depth} overflows i32 (sample {sample_idx:?})")]
    DepthOverflow {
        chrom_id: u32,
        pos: u32,
        sample_idx: Option<usize>,
        depth: u64,
    },

    /// A contig length exceeds `i32::MAX`, which htslib treats as the
    /// `##contig=length=` cap.
    #[error("contig '{name}' length {length} exceeds i32::MAX (VCF contig length limit)")]
    ContigLengthOverflow { name: String, length: u32 },
}
