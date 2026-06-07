//! Per-sample pileup record data model — the unit of interchange
//! between every stage of the pipeline.
//!
//! The pileup walker emits these records; the PSP format serialises
//! and deserialises them; the downstream stages (DUST filter, variant
//! grouper, per-group merger, contamination estimator, posterior
//! engine, VCF writer) consume them. The types live here, in a
//! peer-of-stages module, so no consumer has to reach into any
//! particular stage's namespace just to talk about per-sample
//! pileup data.
//!
//! - [`PileupRecord`] — one record per covered reference position.
//! - [`AlleleObservation`] — one allele bucket inside a record;
//!   `alleles[0]` is unconditionally REF.
//! - [`AlleleSupportStats`] — the fixed-shape per-allele scalar
//!   statistics (read counts, quality sums, strand / placement /
//!   MAPQ moments).
//! - [`ChainId`] — phase-chain identifier minted by the walker's
//!   allocator and threaded through every observation.
//!
//! See `doc/devel/specs/calling_pipeline_architecture.md` §"The five
//! per-allele scalars" for the field-by-field justification of
//! `AlleleSupportStats`, and `doc/devel/specs/phase_chain.md` for the
//! chain-id semantics.

/// Phase-chain identifier. `u64` gives ~1.8 × 10¹⁹ values per file —
/// well beyond any realistic read count. Overflow is caught by the
/// walker's allocator (see the pileup walker's
/// `ChainIdAllocator::allocate_for_read`).
pub type ChainId = u64;

// ---------------------------------------------------------------------
// AlleleSupportStats
// ---------------------------------------------------------------------

/// The fixed set of per-allele support statistics sufficient to
/// reproduce freebayes' likelihood and observation-bias priors
/// exactly. See `doc/devel/specs/calling_pipeline_architecture.md`
/// §"The five per-allele scalars" for the field-by-field
/// justification. Stored compact because every record carries one
/// of these per allele.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[non_exhaustive]
pub struct AlleleSupportStats {
    /// Number of supporting reads for this allele in this record.
    pub num_obs: u32,
    /// `Σ max(ln(P_err_BQ_BAQ), ln(P_err_MQ))` over supporting reads.
    /// Sums in log-error space because freebayes' `prodQout`
    /// reconstruction needs the per-read max combination pre-summed.
    pub q_sum: f64,
    /// Reads on the forward strand among `num_obs`.
    pub fwd: u32,
    /// Reads whose mapped 5′ end is strictly to the left of this
    /// record's anchor position (freebayes' `placedLeft`).
    pub placed_left: u32,
    /// Reads whose mapped 5′ end *is* this record's anchor
    /// position (freebayes' `placedStart`).
    pub placed_start: u32,
    /// Σ mapping quality over supporting reads. Combined with
    /// `num_obs` this yields the mean MAPQ of the allele; with
    /// `mapq_sum_sq` it yields the sample variance and Welch's t.
    pub mapq_sum: u32,
    /// Σ mapq² over supporting reads. With `mapq_sum` and `num_obs`
    /// this gives the sample variance via the identity
    /// `var = (mapq_sum_sq − mapq_sum² / num_obs) / (num_obs − 1)`.
    pub mapq_sum_sq: u64,
}

impl AlleleSupportStats {
    /// Constructor for external code (benches, integration tests).
    /// See [`PileupRecord::new`].
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        num_obs: u32,
        q_sum: f64,
        fwd: u32,
        placed_left: u32,
        placed_start: u32,
        mapq_sum: u32,
        mapq_sum_sq: u64,
    ) -> Self {
        Self {
            num_obs,
            q_sum,
            fwd,
            placed_left,
            placed_start,
            mapq_sum,
            mapq_sum_sq,
        }
    }
}

// ---------------------------------------------------------------------
// AlleleObservation
// ---------------------------------------------------------------------

/// One allele bucket inside a `PileupRecord`. `seq` is the literal
/// allele over `{A,C,G,T,N}`: for SNP/MNP/DEL alleles the length
/// equals the record's REF span (= `alleles[0].seq.len()`);
/// INS-bearing alleles are longer.
///
/// **Observation count.** A `num_obs == 0` entry means "this allele
/// exists in the record but no read in this sample currently has it
/// as its haplotype." Two paths produce this:
/// 1. The REF entry at a position no read covers as REF.
/// 2. A non-REF bucket that *was* observed under an earlier
///    (narrower) footprint and lost its only observer when the
///    record widened and that read re-folded into a different
///    bucket. The empty bucket stays in the record because its
///    seq is still a legal haplotype call across the widened
///    footprint — downstream consumers should treat it as having
///    no evidence in this sample.
///
/// **Chain ids invariant.** `chain_ids` lists exactly the chains of
/// reads currently folded into this bucket — derived from the open-
/// record's fold map at finalise time, not accumulated during
/// folding. Consequences:
/// - `chain_ids.len() <= num_obs as usize` in every bucket (each
///   chain id is backed by at least one observation; paired mates
///   may share a chain id, so the inequality can be strict).
/// - `num_obs == 0 ⇒ chain_ids.is_empty()`.
/// - **The REF bucket (`alleles[0]`) always has `chain_ids.is_empty()`,
///   even when `num_obs > 0`.** The walker deliberately drops REF chain
///   ids at finalise: the only consumer of `chain_ids` is the per-group
///   merger's compound-haplotype check, which links non-reference
///   alleles only and skips allele 0, so REF chain ids are never read.
///   Storing them would be dead weight (~96.6% of all chain ids on real
///   cohorts). See `pileup::walker::open_record::OpenPileupRecord::finalise`
///   and doc/devel/specs/phase_chain.md §8.
///
/// Downstream callers that match constituents by chain-id
/// intersection (e.g. the per-group merger building chain-anchored
/// compounds) can rely on the contrapositive: any bucket reachable
/// via a chain id intersection has `num_obs >= 1`.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct AlleleObservation {
    pub seq: Vec<u8>,
    pub support: AlleleSupportStats,
    /// Distinct phase-chain ids of reads currently folded into this
    /// bucket, sorted ascending and deduplicated.
    pub chain_ids: Vec<ChainId>,
}

impl AlleleObservation {
    /// Constructor for external code (benches, integration tests).
    /// See [`PileupRecord::new`].
    pub fn new(seq: Vec<u8>, support: AlleleSupportStats, chain_ids: Vec<ChainId>) -> Self {
        Self {
            seq,
            support,
            chain_ids,
        }
    }
}

// ---------------------------------------------------------------------
// PileupRecord
// ---------------------------------------------------------------------

/// One emitted per-position record. `alleles[0]` is unconditionally
/// REF — the walker invariant — so `ref_span` is derivable as
/// `alleles[0].seq.len()` and is not stored separately.
///
/// Each `AlleleObservation`'s `chain_ids` field references chain
/// ids minted by the walker's chain-id allocator. Chain ids are
/// unique-per-`.psp`-file `u64` values; there is no recycling, no
/// active-set bookkeeping, and no per-record lifecycle markers.
/// Two allele observations sharing a chain id came from the same
/// read or read-pair in this sample.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct PileupRecord {
    pub chrom_id: u32,
    /// 1-based anchor position.
    pub pos: u32,
    /// At least one entry; `alleles[0]` is always REF.
    pub alleles: Vec<AlleleObservation>,
}

impl PileupRecord {
    /// Number of reference positions this record covers, derived
    /// from REF's `seq` length.
    pub fn ref_span(&self) -> u32 {
        self.alleles[0].seq.len() as u32
    }

    /// Constructor for external code (benches, integration tests).
    /// The walker builds records via struct-update syntax internally
    /// because [`PileupRecord`] is `#[non_exhaustive]`; this `new` is
    /// the supported entry point from outside the crate.
    pub fn new(chrom_id: u32, pos: u32, alleles: Vec<AlleleObservation>) -> Self {
        Self {
            chrom_id,
            pos,
            alleles,
        }
    }
}
