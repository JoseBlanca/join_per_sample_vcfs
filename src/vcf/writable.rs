//! Trait that captures everything the VCF writer reads from a single
//! per-locus posterior record.
//!
//! [`VcfWritable`] is the contract between [`crate::vcf::CohortVcfWriter`]
//! and whatever upstream produces per-locus posterior data. The pipeline's
//! posterior engine emits
//! [`PosteriorRecord`](crate::var_calling::posterior_engine::PosteriorRecord),
//! which implements this trait directly via accessor methods on its
//! existing fields — no intermediate `VcfRow` struct, no per-record
//! allocation. The writer is generic over `R: VcfWritable`; monomorphisation
//! inlines every dispatch back to a direct field access.
//!
//! Why a trait rather than `&PosteriorRecord` everywhere: the writer is
//! `pub` and `src/vcf/` sits at the same top-level layer as
//! [`crate::psp`] / [`crate::fasta`]. Naming the contract here, instead
//! of reaching into [`crate::var_calling::posterior_engine`], makes the
//! dependency direction honest — the writer no longer back-references a
//! pipeline-stage module just to know what shape it consumes.

use crate::pileup_record::AlleleSupportStats;

/// Read-only view of one cohort-level per-locus posterior record that
/// the VCF writer needs.
///
/// Conventions:
/// - **Locus coordinates** are 0-based `chrom_id` and 1-based
///   `pos_1based`, matching the [`PosteriorRecord`](crate::var_calling::posterior_engine::PosteriorRecord)
///   shape and the VCF format's `POS` column.
/// - **Allele indexing**: `allele_idx == 0` is REF; `1..n_alleles` are
///   ALTs. Bytes are uppercase `{A,C,G,T,N}`.
/// - **Per-sample tables** (`posteriors_row`, `scalars_row`,
///   `chain_anchor_flags_row`) return per-sample slices of fixed width
///   (`n_genotypes`, `n_alleles`, `n_alleles` respectively). The writer
///   calls them once per sample per record.
/// - **Flat-length accessors** (`posteriors_len`, `scalars_len`,
///   `chain_anchor_flags_len`) expose the underlying row-major-table
///   length so the writer can defend against shape disagreements at
///   the record-level.
pub trait VcfWritable {
    // ---- Locus + shape -------------------------------------------------

    /// 0-based contig index into the writer's contig table.
    fn chrom_id(&self) -> u32;
    /// 1-based reference position (the VCF `POS` column).
    fn pos_1based(&self) -> u32;
    /// Cohort-wide ploidy.
    fn ploidy(&self) -> u8;
    /// Number of samples this record carries posteriors for.
    fn n_samples(&self) -> usize;
    /// Number of canonical genotypes per sample —
    /// `genotype_order(ploidy, n_alleles).len()`.
    fn n_genotypes(&self) -> usize;
    /// Number of alleles, including REF.
    fn n_alleles(&self) -> usize;
    /// Uppercase ASCII allele bytes for `allele_idx`. `allele_idx == 0` is REF.
    fn allele_seq(&self, allele_idx: usize) -> &[u8];

    // ---- Cohort-level scalars ------------------------------------------

    /// Site-level QUAL in Phred. `f64::INFINITY` is allowed; the writer
    /// caps to a finite VCF-displayable value.
    fn qual_phred(&self) -> f64;
    /// True iff the EM converged within the iteration cap. Drives the
    /// `EMNoConv` FILTER tag.
    fn converged(&self) -> bool;
    /// Per-allele frequency estimates `p̂`. Length equals `n_alleles`.
    fn allele_frequencies(&self) -> &[f64];
    /// Per-compound cohort frequency estimates. `Some(f)` for compound
    /// alleles, `None` otherwise. Length equals `n_alleles`.
    fn compound_frequencies(&self) -> &[Option<f64>];

    // ---- Per-sample tables ---------------------------------------------

    /// Per-sample argmax genotype, indexed by genotype-enumeration
    /// position. Length equals `n_samples`.
    fn best_genotype(&self) -> &[usize];
    /// Per-sample genotype quality in Phred. Length equals `n_samples`.
    fn gq_phred(&self) -> &[f64];
    /// Per-sample slice of `n_genotypes` posterior probabilities. Each
    /// row sums to 1 within tolerance.
    fn posteriors_row(&self, sample_idx: usize) -> &[f64];
    /// Per-sample slice of `n_alleles` allele-support scalars (read
    /// counts, MAPQ moments, strand / placement counts).
    fn scalars_row(&self, sample_idx: usize) -> &[AlleleSupportStats];
    /// Per-sample slice of `n_alleles` "chain anchor" flags (which
    /// alleles in this record are supported by a chain anchored at this
    /// locus, vs forwarded from an upstream chain).
    fn chain_anchor_flags_row(&self, sample_idx: usize) -> &[bool];

    // ---- Flat-length accessors (validation) ----------------------------

    /// Total length of the row-major `posteriors` table; the writer
    /// defends against `posteriors_len() != n_samples * n_genotypes`.
    fn posteriors_len(&self) -> usize;
    /// Total length of the row-major `scalars` table;
    /// `scalars_len() == n_samples * n_alleles`.
    fn scalars_len(&self) -> usize;
    /// Total length of the row-major `chain_anchor_flags` table;
    /// `chain_anchor_flags_len() == n_samples * n_alleles`.
    fn chain_anchor_flags_len(&self) -> usize;
}
