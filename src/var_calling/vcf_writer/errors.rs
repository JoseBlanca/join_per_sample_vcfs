//! Error type for the cohort VCF writer.

use std::io;

use thiserror::Error;

/// Errors surfaced by [`CohortVcfWriter`](super::CohortVcfWriter) and
/// its helpers.
///
/// `Io` covers sink errors (file create, write, rename). `Encode`
/// wraps anything `noodles-vcf` flags at encode time. `InvalidMetadata`
/// is the constructor-time validation on [`CohortMetadata`](super::CohortMetadata).
/// `RecordOutOfOrder` latches if the upstream surfaces records that
/// regress in `(chrom_id, pos)`. `GenotypeIndexOutOfBounds` and
/// `ContigLengthOverflow` are defensive guards against malformed
/// upstream input.
#[derive(Debug, Error)]
pub enum VcfWriteError {
    #[error("io: {0}")]
    Io(#[from] io::Error),

    #[error("vcf encode: {0}")]
    Encode(String),

    #[error("invalid metadata: {0}")]
    InvalidMetadata(String),

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

    #[error(
        "record at {chrom_id}:{pos}: chrom_id is out of bounds for the \
         {n_contigs} contig(s) declared in the cohort metadata"
    )]
    UnknownChromId {
        chrom_id: u32,
        pos: u32,
        n_contigs: usize,
    },

    #[error(
        "record at {chrom_id}:{pos}: posterior arrays were sized for \
         {expected_samples} samples but the cohort metadata names {got_samples}"
    )]
    SampleCountMismatch {
        chrom_id: u32,
        pos: u32,
        expected_samples: usize,
        got_samples: usize,
    },

    #[error("contig '{name}' length {length} exceeds i32::MAX (VCF contig length limit)")]
    ContigLengthOverflow { name: String, length: u32 },
}
