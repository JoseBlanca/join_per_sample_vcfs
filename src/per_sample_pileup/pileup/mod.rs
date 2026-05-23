//! Stage 1's sequential pileup walker — turns a coordinate-sorted
//! stream of prepared (filtered + BAQ-adjusted) reads into the
//! per-position records that Stage 2 encodes to the `.psp` file.
//!
//! See `ia/specs/pileup_walker.md` for the implementation-ready
//! specification. This module ships the public types
//! (`PreparedRead`, `PileupRecord`, `AlleleObservation`, `AlleleSupportStats`,
//! `ChainId`, `MultiChromRefFetcher`, `WalkerError`) and the public entry
//! point `run`. Internal building blocks live in private submodules.

mod active_read_set;
mod chain_id_allocator;
mod cigar_cursor;
mod decompose;
mod errors;
mod open_record;
mod walker;

#[cfg(test)]
pub(crate) mod tests;

use std::sync::Arc;

pub use crate::per_sample_pileup::cram_input::CigarOp;
pub use crate::per_sample_pileup::ref_fetcher::ChromRefFetchError;
pub use chain_id_allocator::{ChainId, DEFAULT_MAX_ACTIVE_READS};
pub use errors::WalkerError;
pub use walker::{PileupWalker, RunSummary, run};

// ---------------------------------------------------------------------
// Defaults / tunables
// ---------------------------------------------------------------------

/// Default for [`WalkerConfig::max_record_span`]. Bounds per-record
/// memory by capping any open record's `ref_span`. Enforced both
/// upstream (the filter stage drops reads whose CIGAR walk plus
/// mate gap consumes more than this) and inside the walker
/// (`widen` / `open_new` error with [`WalkerError::RecordTooWide`]
/// if the cap is reached). Bounds memory rather than serving as a
/// closure trigger; closure is per-record (see
/// `ia/specs/pileup_walker.md` §"Closure rule").
pub const DEFAULT_MAX_RECORD_SPAN: u32 = 5_000;

/// Default for [`WalkerConfig::mate_lookup_window`]. Defensive
/// bound on how far past the first mate's `alignment_start` the
/// walker will keep an entry in `pending_mates` before evicting
/// it and treating the first mate as solo. Sized for typical
/// Illumina paired-end fragments (200–500 bp insert; outer-mate
/// distance under a few kb on jumping libraries) with comfortable
/// headroom. Long-read paired protocols with mate gaps beyond
/// this would need it raised.
pub const DEFAULT_MATE_LOOKUP_WINDOW: u32 = 10_000;

/// Default SNP/REF-column contributor cap, adopted from samtools'
/// `MPLP_MAX_DEPTH`. Pathologically deep regions truncate at this
/// many post-overlap contributors (admission order). See
/// [`WalkerConfig`] for the truncation contract and bias caveats.
pub const DEFAULT_MAX_SNP_COLUMN_DEPTH: u32 = 8_000;

/// Default indel-column contributor cap, adopted from samtools'
/// `MPLP_MAX_INDEL_DEPTH`. Tighter than the SNP cap because
/// homopolymer-context indels saturate beyond what the likelihood
/// math needs.
pub const DEFAULT_MAX_INDEL_COLUMN_DEPTH: u32 = 250;

// ---------------------------------------------------------------------
// WalkerConfig
// ---------------------------------------------------------------------

/// Tunable thresholds the walker reads at run time. Mirrors the
/// per-stage pattern of `CramMergedReaderConfig` (input stage):
/// every value here is a parameter the future CLI binds to. New
/// walker tunables go in this struct rather than as bare `pub
/// const`s; existing consts are migrated opportunistically.
///
/// **Per-column depth caps** (`max_snp_column_depth`,
/// `max_indel_column_depth`) defend against pathologically deep
/// regions. Adopted from samtools' `MPLP_MAX_DEPTH` (8000) and
/// `MPLP_MAX_INDEL_DEPTH` (250) — the indel cap is far tighter
/// because homopolymer-context indels saturate beyond what the
/// likelihood math needs. The cap is per-column, not per-allele:
/// per-allele clipping would silently bias the allele-frequency
/// estimate (a 99/1 column with cap=250 would report ~71/29).
///
/// **Truncation contract.** When the contributor count at a column
/// exceeds the cap, the walker keeps the first `cap` contributors
/// in **admission order** (= the order `cram_input` delivered them,
/// which under normal pipeline operation is coordinate-then-arrival
/// order across the merged CRAMs) and drops the rest. This is *not*
/// a uniform random sample. It is approximately unbiased under the
/// typical assumption that reads at a given coordinate are not
/// allele-correlated — random-fragment shotgun shouldn't cluster
/// alleles by arrival order. Pipelines that batch coverage by
/// source (lane-stratified, region-restricted CRAMs concatenated
/// rather than merged) can show batch-direction bias at over-cap
/// columns; the fix in that case is upstream interleaving, not a
/// random sampler here — the cap exists as a defensive truncation,
/// not a statistical sampler. Pinned by
/// `column_depth_cap_keeps_first_n_of_admission_order`. Mi5 in
/// `ia/reviews/pileup_2026-05-09.md`.
///
/// With the current `MAX_ACTIVE_READS = 4096`, only the indel cap
/// can fire (column depth is already bounded below the SNP cap).
/// The SNP cap is future-proofing for if/when the active-reads cap
/// is raised.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct WalkerConfig {
    /// Maximum contributors to fold at a column whose alleles are
    /// all SNP/REF (no indel events at this anchor). Defaults to
    /// [`DEFAULT_MAX_SNP_COLUMN_DEPTH`].
    pub max_snp_column_depth: u32,
    /// Maximum contributors to fold at a column carrying any
    /// indel observation at this anchor. Defaults to
    /// [`DEFAULT_MAX_INDEL_COLUMN_DEPTH`].
    pub max_indel_column_depth: u32,
    /// Hard cap on per-record reference span, enforced inside the
    /// walker as a defensive bound (`widen` / `open_new` error
    /// with [`WalkerError::RecordTooWide`]). Defaults to
    /// [`DEFAULT_MAX_RECORD_SPAN`]. Upstream filtering drops reads
    /// whose CIGAR walk plus mate gap consumes more than this, so
    /// in normal operation the walker-side cap is never reached.
    pub max_record_span: u32,
    /// Defensive bound on how far past the first mate's
    /// `alignment_start` the walker will keep a `pending_mates`
    /// entry before evicting it and treating the first mate as
    /// solo. Defaults to [`DEFAULT_MATE_LOOKUP_WINDOW`].
    pub mate_lookup_window: u32,
    /// Hard cap on the number of concurrently-active reads.
    /// Exceeding it is a hard error
    /// ([`WalkerError::ActiveReadsExhausted`]) rather than silent
    /// memory blow-up on pathologically deep regions. Defaults to
    /// [`DEFAULT_MAX_ACTIVE_READS`].
    pub max_active_reads: u32,
}

impl Default for WalkerConfig {
    fn default() -> Self {
        Self {
            max_snp_column_depth: DEFAULT_MAX_SNP_COLUMN_DEPTH,
            max_indel_column_depth: DEFAULT_MAX_INDEL_COLUMN_DEPTH,
            max_record_span: DEFAULT_MAX_RECORD_SPAN,
            mate_lookup_window: DEFAULT_MATE_LOOKUP_WINDOW,
            max_active_reads: DEFAULT_MAX_ACTIVE_READS,
        }
    }
}

// ---------------------------------------------------------------------
// MateRole
// ---------------------------------------------------------------------

/// SAM-flag-derived mate role for a [`PreparedRead`]. Collapses the
/// SAM `0x1` (paired) and `0x40` (first segment) bits into a single
/// field so combinations the walker never wants to see — notably
/// "solo + first-of-pair" — cannot be constructed.
///
/// "First" / "Second" here is the SAM-flag sense (which segment of
/// the template the read came from), not the walker arrival-order
/// sense — under coordinate-sorted input the second-of-pair can be
/// the first to reach the walker.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MateRole {
    /// SAM flag `0x1` unset. Treated as solo by the walker and never
    /// registered in `pending_mates`.
    Solo,
    /// SAM flag `0x1` set, flag `0x40` set — first segment in the
    /// template.
    FirstOfPair,
    /// SAM flag `0x1` set, flag `0x40` unset — last segment in the
    /// template.
    SecondOfPair,
}

impl MateRole {
    /// True for both `FirstOfPair` and `SecondOfPair` — i.e. SAM
    /// flag `0x1` is set.
    pub fn is_paired(self) -> bool {
        !matches!(self, MateRole::Solo)
    }

    /// True only for `FirstOfPair`. Used by mate-overlap tie-breaks
    /// on equal-BQ positions; on `Solo` and `SecondOfPair` it is
    /// false.
    pub fn is_first_of_pair(self) -> bool {
        matches!(self, MateRole::FirstOfPair)
    }
}

// ---------------------------------------------------------------------
// PreparedRead
// ---------------------------------------------------------------------

/// A read after the upstream filter and BAQ stages: every field the
/// walker needs is owned and pre-decoded. The walker treats `bq_baq`
/// as opaque per-base BQ — it does not re-apply BAQ — so for tests
/// or until the BAQ stage lands, raw BAM BQ flows through unchanged.
///
/// All bytes are uppercase ASCII over `{A,C,G,T,N}` (`seq`) or Phred
/// 0–93 (`bq_baq`). `qname` is shared as `Arc<str>` so cheap clones
/// can sit in the `pending_mates` map alongside the `ActiveRead`
/// without the bytes being duplicated.
// `#[non_exhaustive]` is deliberately NOT applied to `PreparedRead`.
// PreparedRead is the input contract: callers must populate it with
// concrete bytes per the field-level docs, and adding a new field
// should force every caller (test, bench, production constructor)
// to update its literal explicitly. `#[non_exhaustive]` would push
// callers toward `..Default::default()`, which is exactly the
// silent-absorb hazard the refactor-safety rule guards against.
#[derive(Debug, Clone)]
pub struct PreparedRead {
    /// Index into the merged `ContigList`.
    pub chrom_id: u32,
    /// 1-based reference position of the first reference base the
    /// read covers (CIGAR `M`/`=`/`X` aligned).
    pub alignment_start: u32,
    /// 1-based reference position of the last reference base the
    /// read covers, inclusive. Cached at decode time so the walker
    /// does not re-walk the CIGAR to compute it.
    pub alignment_end: u32,
    pub cigar: Vec<CigarOp>,
    /// Read bases. Length matches the read-consuming CIGAR ops
    /// (M/=/X/I/S).
    pub seq: Vec<u8>,
    /// BAQ-capped per-base quality, Phred. Same length as `seq`.
    /// Outside test contexts this is `min(BQ, BAQ)`; inside tests
    /// (or before BAQ lands) raw BAM BQ flows through unchanged.
    pub bq_baq: Vec<u8>,
    /// `ln(P_err)` derived from MAPQ once at filter time. Cached
    /// because every event-folding step uses it.
    pub mq_log_err: f64,
    /// Raw BAM mapping quality, preserved alongside `mq_log_err` so
    /// per-allele `mapq_sum` / `mapq_sum_sq` can be accumulated in
    /// the pileup fold. Reads that pass `--min-mapq` carry their
    /// original Phred MAPQ here (BWA-MEM caps at 60 in practice).
    pub mapq: u8,
    pub is_reverse_strand: bool,
    pub qname: Arc<str>,
    /// SAM-flag-derived mate role. Collapses the SAM `0x1` (paired)
    /// and `0x40` (first segment) bits into one field so the
    /// previously-unrepresentable "solo + first-of-pair" combination
    /// is gone at the type level. The walker consults
    /// [`MateRole::is_paired`] to decide whether to register the
    /// qname for second-mate lookup and [`MateRole::is_first_of_pair`]
    /// as the deterministic tie-breaker on equal-BQ mate-overlap
    /// positions.
    pub mate_role: MateRole,
    /// 1-based reference position of the first read base that lies
    /// inside the mate-pair adaptor. `None` when the boundary cannot
    /// be reliably computed (single-end, mate unmapped, geometry
    /// inconsistent, or molecule longer than the read).
    ///
    /// On the **forward** strand any read base at `ref_pos >=
    /// adaptor_boundary` was sequenced *through* the molecule's 3′
    /// end into the far adaptor. On the **reverse** strand any read
    /// base at `ref_pos <= adaptor_boundary` was sequenced through
    /// the 5′ end into the near adaptor. The cursor's Match-emit
    /// sites apply this test direction-aware. See finding `G1` in
    /// `ia/reviews/pileup_gatk_comparison_2026-05-08.md`.
    ///
    /// # Default
    /// `None` *disables the G1 adaptor filter* for this read: no
    /// read base is treated as past the boundary. This is the
    /// safe choice when the upstream cannot trust its mate
    /// geometry — a false-positive filter would drop real
    /// evidence. Callers that need strict filtering must compute
    /// and set a boundary themselves. Mi17 in
    /// `ia/reviews/pileup_2026-05-11.md`.
    pub adaptor_boundary: Option<u32>,
}

// ---------------------------------------------------------------------
// ReadLengthError
// ---------------------------------------------------------------------

/// Failure modes for [`PreparedRead::length`]. Carries the raw
/// lengths only — locus context (qname, chrom_id, pos) is added by
/// whichever caller maps this into its own error type (the walker
/// wraps it in [`WalkerError::MalformedRead`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadLengthError {
    /// `seq.len()` and `bq_baq.len()` disagree.
    SeqBqMismatch { seq_len: usize, bq_baq_len: usize },
    /// The CIGAR's read-consuming ops sum to `cigar_consumed` but
    /// the `seq` field has length `seq_len`.
    CigarSeqMismatch { cigar_consumed: u64, seq_len: usize },
}

impl PreparedRead {
    /// The read's length (in `u64`, matching the cigar-sum width),
    /// verified against the three internal representations that
    /// must agree on it:
    ///
    /// 1. `seq.len() == bq_baq.len()`.
    /// 2. The CIGAR's read-consuming ops (M/=/X/I/S) sum to
    ///    `seq.len()`.
    ///
    /// On success the three lengths agree and the common value is
    /// returned. On failure a typed [`ReadLengthError`] reports
    /// which invariant broke.
    ///
    /// M21 in `ia/reviews/pileup_2026-05-09.md`. The cigar cursor
    /// and `decompose` both index `seq[..]` / `bq_baq[..]` using
    /// offsets derived from the CIGAR; a mismatch would otherwise
    /// panic with `slice index out of bounds` and kill the run on
    /// the offending read. The accessor surfaces a typed error
    /// instead so callers can attach locus context and continue.
    pub fn length(&self) -> Result<u64, ReadLengthError> {
        if self.seq.len() != self.bq_baq.len() {
            return Err(ReadLengthError::SeqBqMismatch {
                seq_len: self.seq.len(),
                bq_baq_len: self.bq_baq.len(),
            });
        }
        let cigar_consumed: u64 = self
            .cigar
            .iter()
            .map(|op| match *op {
                CigarOp::Match(n)
                | CigarOp::SeqMatch(n)
                | CigarOp::SeqMismatch(n)
                | CigarOp::Insertion(n)
                | CigarOp::SoftClip(n) => n as u64,
                CigarOp::Deletion(_)
                | CigarOp::Skip(_)
                | CigarOp::HardClip(_)
                | CigarOp::Padding(_) => 0,
            })
            .sum();
        if cigar_consumed != self.seq.len() as u64 {
            return Err(ReadLengthError::CigarSeqMismatch {
                cigar_consumed,
                seq_len: self.seq.len(),
            });
        }
        Ok(cigar_consumed)
    }
}

// ---------------------------------------------------------------------
// AlleleSupportStats
// ---------------------------------------------------------------------

/// The fixed set of per-allele support statistics sufficient to
/// reproduce freebayes' likelihood and observation-bias priors
/// exactly. See `ia/specs/calling_pipeline_architecture.md` §"The
/// five per-allele scalars" for the field-by-field justification.
/// Stored compact because every record carries one of these per
/// allele.
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
// MultiChromRefFetcher
// ---------------------------------------------------------------------

/// Multi-chromosome reference-FASTA fetcher used by the Stage 1
/// pileup walker. Errors are typed (`ChromRefFetchError`) so
/// callers can route I/O failures, range failures, and contract
/// violations distinctly. Single-chromosome consumers (DUST, BAQ,
/// PerGroupMerger) should use [`ChromRefFetcher`] instead — it
/// drops the `chrom_id` parameter and exposes a sliding-buffer
/// contract specifically for monotonic-forward access.
///
/// [`ChromRefFetcher`]: crate::per_sample_pileup::ref_fetcher::ChromRefFetcher
pub trait MultiChromRefFetcher {
    /// Fetch `length` reference bases starting at the 1-based
    /// position `start` on chromosome `chrom_id`. Bytes are
    /// uppercase ASCII over `{A,C,G,T,N}` (canonicalised by the
    /// fetcher implementation).
    ///
    /// # Errors
    ///
    /// - [`ChromRefFetchError::OutOfBounds`] if the requested
    ///   window exceeds the chromosome length.
    /// - [`ChromRefFetchError::InvalidStart`] if `start_1based == 0`.
    /// - [`ChromRefFetchError::Io`] on any underlying FASTA I/O
    ///   failure or unknown `chrom_id`.
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, ChromRefFetchError>;

    /// Forward sequential iterator over every uppercased base of
    /// `chrom_id`'s contig, in 1..=`length` order. Used by the DUST
    /// mask construction (one pass per chrom) to avoid materialising
    /// the whole contig as a single `Vec<u8>`.
    ///
    /// Default impl materialises via `fetch(chrom_id, 1, length)`;
    /// streaming fetchers override to walk a sliding buffer instead.
    /// The boxed iterator costs one heap allocation per chrom plus
    /// one virtual dispatch per byte — the latter is dominated by
    /// the inner sdust scoring work in practice.
    ///
    /// # Errors
    ///
    /// Same failure modes as [`Self::fetch`] (the default impl
    /// just calls `fetch(chrom_id, 1, length)`).
    fn iter_bases<'a>(
        &'a self,
        chrom_id: u32,
        length: u32,
    ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>
    {
        let seq = self.fetch(chrom_id, 1, length)?;
        Ok(Box::new(seq.into_iter().map(Ok)))
    }
}

/// Forwarding impl so callers may pass either an owned fetcher or a
/// shared reference into [`PileupWalker::new`] / [`run`]. The walker
/// only ever calls `&self` methods, so the borrow is sufficient.
impl<T: MultiChromRefFetcher + ?Sized> MultiChromRefFetcher for &T {
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, ChromRefFetchError> {
        (**self).fetch(chrom_id, start_1based, length)
    }

    fn iter_bases<'a>(
        &'a self,
        chrom_id: u32,
        length: u32,
    ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>
    {
        (**self).iter_bases(chrom_id, length)
    }
}
