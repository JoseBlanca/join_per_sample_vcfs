//! Stage 1's sequential pileup walker — turns a coordinate-sorted
//! stream of prepared (filtered + BAQ-adjusted) reads into the
//! per-position records that Stage 2 encodes to the `.psf` file.
//!
//! See `ia/specs/pileup_walker.md` for the implementation-ready
//! specification. This module ships the public types
//! (`PreparedRead`, `PileupRecord`, `AlleleObservation`, `FiveScalars`,
//! `SlotId`, `RefBaseFetcher`, `WalkerError`) and the public entry
//! point `run`. Internal building blocks live in private submodules.

mod active_set;
mod cigar_cursor;
mod decompose;
mod errors;
mod open_record;
mod slot_allocator;
mod walker;

#[cfg(test)]
mod tests;

use std::sync::Arc;

pub use crate::per_sample_caller::cram_input::CigarOp;
pub use errors::WalkerError;
pub use slot_allocator::{DEFAULT_MAX_ACTIVE_SLOTS, SlotId};
pub use walker::{RunSummary, run};

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

/// Default channel buffer between the walker thread and the Stage 2
/// encoder thread. Sized for the worst-case burst of simultaneous
/// closures plus enough headroom for a few subsequent steps to
/// proceed before backpressure kicks in.
pub const DEFAULT_OUTPUT_CHANNEL_CAPACITY: usize = 64;

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
/// With the current `MAX_ACTIVE_SLOTS = 4096`, only the indel cap
/// can fire (column depth is already bounded below the SNP cap).
/// The SNP cap is future-proofing for if/when the slot cap is
/// raised.
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
    /// Hard cap on the number of concurrently-active phase-chain
    /// slots. Exceeding it is a hard error
    /// ([`WalkerError::SlotExhausted`]) rather than silent slot
    /// reuse. Defaults to [`DEFAULT_MAX_ACTIVE_SLOTS`].
    pub max_active_slots: u32,
}

impl Default for WalkerConfig {
    fn default() -> Self {
        Self {
            max_snp_column_depth: DEFAULT_MAX_SNP_COLUMN_DEPTH,
            max_indel_column_depth: DEFAULT_MAX_INDEL_COLUMN_DEPTH,
            max_record_span: DEFAULT_MAX_RECORD_SPAN,
            mate_lookup_window: DEFAULT_MATE_LOOKUP_WINDOW,
            max_active_slots: DEFAULT_MAX_ACTIVE_SLOTS,
        }
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
    pub is_reverse_strand: bool,
    pub qname: Arc<str>,
    /// SAM flag `0x40`. Used as the deterministic tie-breaker on
    /// equal-BQ mate-overlap positions.
    pub is_first_mate: bool,
    /// SAM flag `0x1`. When unset, the read is treated as solo and
    /// is not registered in `pending_mates`.
    pub has_mate: bool,
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
// FiveScalars
// ---------------------------------------------------------------------

/// The fixed set of per-allele scalars sufficient to reproduce
/// freebayes' likelihood and observation-bias priors exactly. See
/// `ia/specs/calling_pipeline_architecture.md` §"The five per-allele
/// scalars" for the field-by-field justification. Stored compact
/// because every record carries one of these per allele.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[non_exhaustive]
pub struct FiveScalars {
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
}

// ---------------------------------------------------------------------
// AlleleObservation
// ---------------------------------------------------------------------

/// One allele bucket inside a `PileupRecord`. `seq` is the literal
/// allele over `{A,C,G,T,N}`: for SNP/MNP/DEL alleles the length
/// equals the record's REF span (= `alleles[0].seq.len()`);
/// INS-bearing alleles are longer.
///
/// A `num_obs == 0` entry is valid and means "this allele exists in
/// the record but no read in this sample observed it" — only the
/// REF entry can be in this state under the current walker (allele
/// buckets for non-REF alleles are created lazily on first event,
/// so they always have `num_obs ≥ 1`).
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct AlleleObservation {
    pub seq: Vec<u8>,
    pub scalars: FiveScalars,
    /// Distinct phase-chain slot ids that contributed to this
    /// allele, sorted ascending and deduplicated.
    pub chain_slots: Vec<SlotId>,
}

// ---------------------------------------------------------------------
// PileupRecord
// ---------------------------------------------------------------------

/// One emitted per-position record. `alleles[0]` is unconditionally
/// REF — the walker invariant — so `ref_span` is derivable as
/// `alleles[0].seq.len()` and is not stored separately.
///
/// `new_chains` and `expired_chains` carry the lifecycle markers
/// for phase-chain slots: ids that started or ended since the
/// previous emitted record. Both lists are deduplicated and may be
/// emitted in any order; Stage 2 may reorder them.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct PileupRecord {
    pub chrom_id: u32,
    /// 1-based anchor position.
    pub pos: u32,
    pub new_chains: Vec<SlotId>,
    pub expired_chains: Vec<SlotId>,
    /// At least one entry; `alleles[0]` is always REF.
    pub alleles: Vec<AlleleObservation>,
}

impl PileupRecord {
    /// Number of reference positions this record covers, derived
    /// from REF's `seq` length.
    pub fn ref_span(&self) -> u32 {
        self.alleles[0].seq.len() as u32
    }
}

// ---------------------------------------------------------------------
// RefBaseFetcher
// ---------------------------------------------------------------------

/// What the walker needs from the reference: a way to fetch the
/// literal bases over a `[start, start + length)` window on a
/// chromosome, 1-based start. Implemented by the production wrapper
/// over `noodles_fasta::Repository`; tests can plug in a mock that
/// reads from an in-memory string.
pub trait RefBaseFetcher {
    /// Fetch `length` reference bases starting at the 1-based
    /// position `start` on chromosome `chrom_id`. Bases are
    /// uppercase ASCII over `{A,C,G,T,N}`. Returns `Err` if the
    /// range exceeds the chromosome's length or the fetch itself
    /// fails (FASTA I/O).
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, std::io::Error>;
}
