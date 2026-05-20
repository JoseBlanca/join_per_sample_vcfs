//! Stage 5 — per-group merger.
//!
//! Consumes `Result<OverlappingVariantGroup, GrouperError>` items from
//! Stage 4 and emits one `MergedRecord` per group: a unified allele
//! set plus the per-sample scalar table and per-sample, per-genotype
//! log-likelihood vector under `c_s = 0` (no contamination
//! correction). See
//! `doc/devel/implementation_plans/cohort_per_group_merger.md` and
//! the algorithmic contract in
//! `doc/devel/specs/calling_pipeline_architecture.md` §"Stage 5".
//!
//! The merger projects every sample-local allele onto the group's
//! reference span, deduplicates by byte equality, admits cross-record
//! compound haplotypes only when chain-id evidence inside at least one
//! sample links the constituents, and reconstructs the freebayes
//! closed-form likelihood from the five per-allele scalars carried in
//! `AlleleSupportStats`. Chain-broken samples at chain-anchored
//! compounds fall back to a constituents-independent per-position
//! likelihood and surface that fact through `chain_anchor_flags` so Stage 6 can
//! use a cohort-derived compound frequency as the prior signal there.
//!
//! The merger pulls a batch of upstream groups, processes them
//! serially, and emits the resulting records one-by-one in input
//! order. Parallelism above the merger lives at the per-chromosome
//! seam in [`crate::pop_var_caller::cohort_driver`]; running an inner
//! `rayon::par_iter` here would only nest under that outer
//! decomposition.

use std::collections::{BTreeMap, VecDeque};
use std::sync::Arc;

use ahash::AHashMap;
use smallvec::SmallVec;
use thiserror::Error;

/// Per-sample list of `(record_idx, local_allele_idx)` pairs used by
/// the unify-alleles bookkeeping. Inline-1 covers the dominant
/// biallelic case (≤ 1 source per sample); chain-anchored compound
/// constituents (2+) spill to the heap, matching the previous `Vec`
/// behaviour for that case. Sized to keep the per-allele
/// `Vec<SourceList>` table cheap to allocate and drop per group —
/// the dominant Stage 5 allocator cost per the perf review.
type SourceList = SmallVec<[(usize, usize); 1]>;

use crate::per_sample_pileup::pileup::{AlleleSupportStats, ChainId, RefSeqFetcher};
use crate::var_calling::variant_grouping::{GrouperError, OverlappingVariantGroup};

/// Maximum number of alleles retained in a single merged record
/// before the cap drops the lowest-cohort-count alleles into the
/// per-sample `other_scalars` pool. Matches GATK's
/// `maxAlternateAlleles` default.
pub const DEFAULT_MAX_ALLELES_PER_RECORD: usize = 6;

/// Default cohort-wide ploidy. Diploid is the common case; the merger
/// supports higher ploidies (triploid, tetraploid) through
/// [`PerGroupMergerConfig::ploidy`]. See
/// `doc/devel/specs/calling_pipeline_architecture.md` §"Stage 5 —
/// per-group processing" for the ploidy contract.
pub const DEFAULT_PLOIDY: u8 = 2;

/// Default number of upstream groups pulled and processed in parallel
/// per `next()` cycle. **Placeholder** — chosen as a reasonable
/// starting point for cache-friendly batches; no comparative benchmark
/// has selected this specific value yet. Raise it on very wide cohorts
/// where per-group work amortises better, and re-evaluate when the
/// `var_calling_per_group_merger/*` criterion bench gets a dedicated
/// batch-size sweep.
pub const DEFAULT_BATCH_SIZE: usize = 32;

/// Tunable knobs for the per-group merger.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct PerGroupMergerConfig {
    /// Cohort-wide ploidy. Mixed per-sample ploidies are out of scope
    /// for v1 (see the implementation plan §"Out-of-scope
    /// follow-ups"). Defaults to [`DEFAULT_PLOIDY`].
    pub ploidy: u8,

    /// Hard cap on the number of alleles in a single merged record.
    /// Excess alleles (lowest cohort-wide count first; REF and any
    /// chain-anchored compound protected) drop into a per-sample
    /// `<OTHER>` scalar pool that only contributes to the error-cost
    /// term. Defaults to [`DEFAULT_MAX_ALLELES_PER_RECORD`].
    pub max_alleles: usize,

    /// Number of upstream groups pulled and processed in parallel per
    /// `next()` refill. Defaults to [`DEFAULT_BATCH_SIZE`].
    pub batch_size: usize,
}

impl Default for PerGroupMergerConfig {
    fn default() -> Self {
        Self {
            ploidy: DEFAULT_PLOIDY,
            max_alleles: DEFAULT_MAX_ALLELES_PER_RECORD,
            batch_size: DEFAULT_BATCH_SIZE,
        }
    }
}

/// Inclusive upper bound on [`PerGroupMergerConfig::ploidy`]. Octoploid
/// wheat is the realistic high-water mark for population cohorts; `u8`'s
/// `255` is mathematically permitted but degenerate in practice. Matches
/// the cohort CLI plan's validation table.
pub const MAX_PLOIDY: u8 = 8;

/// Inclusive upper bound on [`PerGroupMergerConfig::max_alleles`]. Wide
/// enough for the most polyallelic real loci; `16` already exceeds the
/// GATK default (`6`) almost 3×. Matches the cohort CLI plan.
pub const MAX_ALLELES_PER_VAR_CAP: usize = 16;

impl PerGroupMergerConfig {
    /// Validated constructor. Range bounds — matching the cohort CLI
    /// plan's "Config validation" table — are enforced here:
    ///
    /// * `ploidy` ∈ `1..=MAX_PLOIDY` (default `2`).
    /// * `max_alleles` ∈ `2..=MAX_ALLELES_PER_VAR_CAP` (default `6`).
    /// * `batch_size` ≥ `1` (a zero batch would never refill).
    ///
    /// Public fields stay accessible for tests and library use; this
    /// constructor is the validated entry point that the CLI funnels
    /// through.
    pub fn new(
        ploidy: u8,
        max_alleles: usize,
        batch_size: usize,
    ) -> Result<Self, PerGroupMergerConfigError> {
        if !(1..=MAX_PLOIDY).contains(&ploidy) {
            return Err(PerGroupMergerConfigError::InvalidPloidy { got: ploidy });
        }
        if !(2..=MAX_ALLELES_PER_VAR_CAP).contains(&max_alleles) {
            return Err(PerGroupMergerConfigError::InvalidMaxAlleles { got: max_alleles });
        }
        if batch_size == 0 {
            return Err(PerGroupMergerConfigError::InvalidBatchSize { got: batch_size });
        }
        Ok(Self {
            ploidy,
            max_alleles,
            batch_size,
        })
    }
}

/// Config-construction errors for [`PerGroupMergerConfig::new`].
#[non_exhaustive]
#[derive(Error, Debug, PartialEq)]
pub enum PerGroupMergerConfigError {
    /// `ploidy` was out of range `1..=MAX_PLOIDY`.
    #[error("ploidy must be in 1..={MAX_PLOIDY}, got {got}")]
    InvalidPloidy { got: u8 },

    /// `max_alleles` was out of range `2..=MAX_ALLELES_PER_VAR_CAP`.
    #[error("max_alleles must be in 2..={MAX_ALLELES_PER_VAR_CAP}, got {got}")]
    InvalidMaxAlleles { got: usize },

    /// `batch_size` was zero. A zero batch would never refill the
    /// merger.
    #[error("batch_size must be >= 1, got {got}")]
    InvalidBatchSize { got: usize },
}

/// One allele in the merged set.
#[derive(Debug, Clone, PartialEq)]
pub struct MergedAllele {
    /// Allele bytes projected onto the `[start, end]` reference span
    /// of the group. REF (`alleles[0]`) is always the bare reference
    /// sequence over that span.
    pub seq: Vec<u8>,
    /// `true` iff this allele is a cross-record compound — its
    /// non-REF bytes come from more than one of the group's records.
    /// The chain-anchor rule has already cleared this entry: at least
    /// one sample's chain-id intersection across the constituents was
    /// non-empty.
    pub is_compound: bool,
    /// For `is_compound = true`: the constituent record + local
    /// allele indices, taken from the *first* sample that
    /// chain-anchored the compound. Different samples may carry the
    /// same compound under different local allele indices because
    /// each sample's record at a given position has its own allele
    /// ordering; this field names a representative source rather than
    /// a cohort-wide canonical id.
    pub constituents: Vec<CompoundConstituent>,
}

/// One constituent of a compound allele: position in the source
/// `OverlappingVariantGroup.records` plus the local allele index in that
/// record (per the first chain-anchoring sample — see
/// [`MergedAllele::constituents`]).
///
/// **The `local_allele_idx` is a representative, not a cohort-canonical
/// identifier.** Different samples may carry the same compound through
/// different `local_allele_idx` values because of per-sample allele
/// packing within a record. To find a given sample's contribution to a
/// compound, walk that sample's records by `record_idx` and match on
/// allele byte sequence, not on `local_allele_idx` equality.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CompoundConstituent {
    pub record_idx: usize,
    pub local_allele_idx: usize,
}

/// One emitted item: a unified allele set for a single
/// `OverlappingVariantGroup`, the per-sample scalar table, the per-sample
/// chain-anchor flags, and the per-sample, per-genotype log-likelihood
/// vector under `c_s = 0`.
///
/// **Storage layout.** The three per-(sample, allele) and
/// per-(sample, genotype) matrices below are row-major flat `Vec`s
/// (sample-major). Use the accessor methods [`Self::scalars_row`],
/// [`Self::chain_anchor_flag`] / [`Self::chain_anchor_flags_row`],
/// and [`Self::log_likelihoods_row`] to read them rather than indexing
/// the buffers directly. The flattening collapses what was previously
/// `n_samples + 1` small heap allocations per matrix per record into a
/// single allocation — the per-group merger's emit path is on Stage 5's
/// allocator-bound critical path per the 2026-05-16 perf review.
#[derive(Debug, Clone, PartialEq)]
pub struct MergedRecord {
    pub chrom_id: u32,
    /// 1-based inclusive start of the group's reference span.
    pub start: u32,
    /// 1-based inclusive end of the group's reference span.
    pub end: u32,
    /// Merged allele set. `alleles[0]` is always REF.
    pub alleles: Vec<MergedAllele>,
    /// Cohort-wide ploidy used for genotype enumeration. Copied from
    /// [`PerGroupMergerConfig::ploidy`].
    pub ploidy: u8,
    /// Number of samples this record carries scalars / likelihoods for.
    pub n_samples: usize,
    /// Number of genotypes per sample in `log_likelihoods`. Equals
    /// `genotype_order(ploidy, alleles.len()).len()`.
    pub n_genotypes: usize,
    /// Flat row-major scalars table — `scalars[sample_idx * alleles.len() + allele_idx]`
    /// is the projected [`AlleleSupportStats`] for that sample/allele.
    /// A sample with no records in the group has every entry zeroed.
    /// Use [`Self::scalars_row`] for the per-sample slice.
    pub scalars: Vec<AlleleSupportStats>,
    /// `other_scalars[sample_idx]` — pooled scalars for alleles that
    /// were dropped by the `max_alleles` cap. Always counts into the
    /// error-cost term so the likelihoods preserve the pre-cap
    /// evidence balance.
    pub other_scalars: Vec<AlleleSupportStats>,
    /// Flat row-major flag table — `chain_anchor_flags[sample_idx * alleles.len() + allele_idx]`
    /// is `true` iff `allele_idx` is a compound and `sample_idx` is
    /// chain-broken at it (the fallback likelihood path was used).
    /// Stage 6 reads this to route the cohort-derived compound
    /// frequency `f_C` as the prior signal for chain-broken samples.
    /// Use [`Self::chain_anchor_flag`] for individual lookups or
    /// [`Self::chain_anchor_flags_row`] for the per-sample slice.
    pub chain_anchor_flags: Vec<bool>,
    /// Flat row-major log-likelihood table —
    /// `log_likelihoods[sample_idx * n_genotypes + genotype_idx]` is
    /// the natural-log likelihood under `c_s = 0`. Genotype order is
    /// canonical per [`genotype_order`]. Use [`Self::log_likelihoods_row`]
    /// for the per-sample slice.
    pub log_likelihoods: Vec<f64>,
}

impl MergedRecord {
    /// Number of samples this record carries scalars/likelihoods for.
    pub fn n_samples(&self) -> usize {
        self.n_samples
    }

    /// Per-sample slice into the flat row-major `scalars` table.
    pub fn scalars_row(&self, sample_idx: usize) -> &[AlleleSupportStats] {
        let n_alleles = self.alleles.len();
        let start = sample_idx * n_alleles;
        &self.scalars[start..start + n_alleles]
    }

    /// Single chain-anchor flag at `(sample_idx, allele_idx)`.
    pub fn chain_anchor_flag(&self, sample_idx: usize, allele_idx: usize) -> bool {
        self.chain_anchor_flags[sample_idx * self.alleles.len() + allele_idx]
    }

    /// Per-sample slice into the flat row-major `chain_anchor_flags` table.
    pub fn chain_anchor_flags_row(&self, sample_idx: usize) -> &[bool] {
        let n_alleles = self.alleles.len();
        let start = sample_idx * n_alleles;
        &self.chain_anchor_flags[start..start + n_alleles]
    }

    /// Per-sample slice into the flat row-major `log_likelihoods` table.
    pub fn log_likelihoods_row(&self, sample_idx: usize) -> &[f64] {
        let start = sample_idx * self.n_genotypes;
        &self.log_likelihoods[start..start + self.n_genotypes]
    }
}

/// Why a likelihood evaluation hit an unrecoverable value. Surfaced
/// through [`PerGroupMergerError::DegenerateLikelihood`] because the
/// closed-form formula is finite for all valid inputs — a degeneracy
/// signals an internal bug, not a data condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DegeneracyKind {
    NaN,
    PositiveInfinity,
}

/// Identifies which compound-likelihood loop surfaced a
/// [`PerGroupMergerError::ZeroObservationConstituent`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompoundPhase {
    /// The quality-gathering loop in `project_scalars` that builds
    /// `min_mean_q` for a chain-anchored compound.
    QualityGather,
    /// The constituent-subtraction loop in `project_scalars` that
    /// removes the compound's claim from each constituent's scalars.
    ConstituentSubtraction,
}

/// Errors surfaced by the per-group merger.
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum PerGroupMergerError {
    // Single-origin: only one `?` in this module surfaces a
    // `GrouperError`. If a second site is ever added, drop `#[from]`
    // here and convert explicitly at each site so error provenance
    // does not silently collapse into the same variant.
    #[error("upstream: {0}")]
    Upstream(#[from] GrouperError),

    #[error("reference fetch at chrom {chrom_id} {start}-{end}: {source}")]
    RefFetch {
        chrom_id: u32,
        start: u32,
        end: u32,
        #[source]
        source: std::io::Error,
    },

    /// Genotype enumeration produced a numerically degenerate
    /// likelihood (NaN or +∞). The closed-form formula is finite for
    /// every valid input, so this signals an internal bug rather than
    /// a data condition.
    #[error(
        "degenerate likelihood at chrom {chrom_id} {start}-{end} \
         for sample_idx {sample_idx} genotype_idx {genotype_idx}: {kind:?}"
    )]
    DegenerateLikelihood {
        chrom_id: u32,
        start: u32,
        end: u32,
        sample_idx: usize,
        genotype_idx: usize,
        kind: DegeneracyKind,
    },

    /// Compound projection asked for a `(record_idx, local_allele_idx)`
    /// pair that no sample in the group carries. The walker invariant
    /// says any admitted compound constituent has at least one
    /// anchoring sample providing the byte sequence.
    #[error(
        "compound projection at chrom {chrom_id} {start}-{end}: no sample \
         carries record_idx={record_idx} local_allele_idx={local_allele_idx}"
    )]
    MissingCompoundAlleleBytes {
        chrom_id: u32,
        start: u32,
        end: u32,
        record_idx: usize,
        local_allele_idx: usize,
    },

    /// A non-REF allele bucket reported `num_obs = 0` while building
    /// per-sample scalars for a chain-anchored compound. Any bucket
    /// reached here was selected via chain-id intersection, and the
    /// walker's `chain_ids` invariant
    /// (`src/per_sample_pileup/pileup/mod.rs` near `AlleleObservation`)
    /// guarantees `chain_ids` lists only reads currently folded into
    /// the bucket — so a non-empty `chain_ids` implies `num_obs >= 1`.
    /// Hitting this error means that invariant has been violated
    /// upstream (e.g. a buggy walker change, or a hand-crafted
    /// `PileupRecord` in tests / fuzz input); silently continuing
    /// would distort the homogeneous-quality approximation or panic
    /// the subtraction step on `inter / 0`.
    #[error(
        "zero-observation constituent for chain-anchored compound at chrom \
         {chrom_id} {start}-{end} (phase: {phase:?}): sample_idx={sample_idx} \
         allele_idx={allele_idx} record_idx={record_idx} \
         local_allele_idx={local_allele_idx}"
    )]
    ZeroObservationConstituent {
        chrom_id: u32,
        start: u32,
        end: u32,
        phase: CompoundPhase,
        sample_idx: usize,
        allele_idx: usize,
        record_idx: usize,
        local_allele_idx: usize,
    },

    /// A chain-anchored compound (`inter > 0`) produced no usable
    /// quality from any of its constituent sources for this sample —
    /// every source either lacked a record at that position or was
    /// already filtered. Walker invariant says `inter > 0` for a
    /// sample implies at least one usable constituent in that sample.
    #[error(
        "chain-anchored compound at chrom {chrom_id} {start}-{end} \
         (sample_idx={sample_idx} allele_idx={allele_idx}) has \
         chain-anchor count {inter} but no constituent yielded a quality"
    )]
    NoQualityForChainAnchoredCompound {
        chrom_id: u32,
        start: u32,
        end: u32,
        sample_idx: usize,
        allele_idx: usize,
        inter: u32,
    },
}

// ---------------------------------------------------------------------
// Genotype enumeration
// ---------------------------------------------------------------------

/// Canonical, cohort-shared enumeration of genotypes for the given
/// ploidy + allele count. Genotypes are non-decreasing index tuples
/// in lexicographic order on the *reversed* tuple — the VCF / GATK PL
/// convention (e.g. for ploidy 2, 3 alleles: `AA, AB, BB, AC, BC,
/// CC`). All samples in a `MergedRecord` use the same order.
pub fn genotype_order(ploidy: u8, n_alleles: usize) -> Vec<Vec<u8>> {
    let ploidy = ploidy as usize;
    let mut enumerated_genotypes: Vec<Vec<u8>> = Vec::new();
    let mut partial_genotype: Vec<u8> = vec![0; ploidy];
    collect_non_decreasing(
        &mut partial_genotype,
        0,
        0,
        n_alleles,
        ploidy,
        &mut enumerated_genotypes,
    );
    enumerated_genotypes.sort_by(|a, b| {
        for i in (0..ploidy).rev() {
            match a[i].cmp(&b[i]) {
                std::cmp::Ordering::Equal => continue,
                other => return other,
            }
        }
        std::cmp::Ordering::Equal
    });
    enumerated_genotypes
}

fn collect_non_decreasing(
    partial_genotype: &mut [u8],
    pos: usize,
    min_allele: u8,
    n_alleles: usize,
    ploidy: usize,
    enumerated_genotypes: &mut Vec<Vec<u8>>,
) {
    if pos == ploidy {
        enumerated_genotypes.push(partial_genotype.to_vec());
        return;
    }
    for a in min_allele..(n_alleles as u8) {
        partial_genotype[pos] = a;
        collect_non_decreasing(
            partial_genotype,
            pos + 1,
            a,
            n_alleles,
            ploidy,
            enumerated_genotypes,
        );
    }
}

// ---------------------------------------------------------------------
// Public merger type
// ---------------------------------------------------------------------

/// Type-erased reference fetcher handle shared across per-chromosome
/// workers. The `Send + Sync` bounds are mandatory because each
/// per-chrom worker constructs its own `PerGroupMerger` and the
/// shared fetcher `Arc` moves across the worker boundary; the merger
/// itself only ever asks for `&` access, so an immutable shared
/// handle is sufficient.
pub type SharedRefFetcher = Arc<dyn RefSeqFetcher + Send + Sync>;

/// Streaming, internally-parallel per-group merger.
pub struct PerGroupMerger<I>
where
    I: Iterator<Item = Result<OverlappingVariantGroup, GrouperError>>,
{
    upstream: I,
    ref_fetcher: SharedRefFetcher,
    config: PerGroupMergerConfig,
    /// Precomputed `genotype_order(config.ploidy, n)` for `n` in
    /// `0..=config.max_alleles`, indexed by `n_alleles`. Built once at
    /// construction; held behind `Arc` so the cohort driver's per-chrom
    /// workers can each construct a `PerGroupMerger` that shares the
    /// same table without rebuild on each worker boundary (the prior
    /// shape allocated `Vec<Vec<u8>>` per `process_group` call). Safe
    /// because `enforce_max_alleles` caps `unified.alleles.len() <=
    /// config.max_alleles` before `compute_log_likelihoods` runs.
    genotype_tables: Arc<Vec<Vec<Vec<u8>>>>,
    /// FIFO of merged records produced by the most recent batch refill.
    /// The iterator drains this before pulling another batch.
    /// REF-only groups are filtered out inside `refill` before they
    /// land here, so this queue only carries real records and
    /// surfaced errors.
    pending: VecDeque<Result<MergedRecord, PerGroupMergerError>>,
    /// Latched after the first surfaced error or upstream exhaustion.
    done: bool,
}

/// Build the per-(ploidy, n_alleles) genotype table cache used by
/// `PerGroupMerger`. Tests that drive `process_group` directly call
/// this with the test's own config so the cache slice matches.
fn build_genotype_tables(config: &PerGroupMergerConfig) -> Vec<Vec<Vec<u8>>> {
    (0..=config.max_alleles)
        .map(|n_alleles| genotype_order(config.ploidy, n_alleles))
        .collect()
}

impl<I> std::fmt::Debug for PerGroupMerger<I>
where
    I: Iterator<Item = Result<OverlappingVariantGroup, GrouperError>>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Exhaustive destructure: a new field on `PerGroupMerger`
        // fails to compile here so the omission from Debug output is
        // explicit rather than silent.
        let Self {
            upstream: _,
            ref_fetcher: _,
            config,
            genotype_tables: _,
            pending,
            done,
        } = self;
        f.debug_struct("PerGroupMerger")
            .field("config", config)
            .field("pending_len", &pending.len())
            .field("done", done)
            .finish()
    }
}

impl<I> PerGroupMerger<I>
where
    I: Iterator<Item = Result<OverlappingVariantGroup, GrouperError>>,
{
    /// Construct a merger with explicit tuning. Pass
    /// [`PerGroupMergerConfig::default()`] for the standard defaults
    /// (`ploidy = 2`, `max_alleles = 6`, `batch_size = 32`); pass an
    /// explicit value for any field to override.
    pub fn with_config(
        upstream: I,
        ref_fetcher: SharedRefFetcher,
        config: PerGroupMergerConfig,
    ) -> Self {
        let genotype_tables = Arc::new(build_genotype_tables(&config));
        Self {
            upstream,
            ref_fetcher,
            config,
            genotype_tables,
            pending: VecDeque::new(),
            done: false,
        }
    }

    pub fn config(&self) -> &PerGroupMergerConfig {
        &self.config
    }

    /// Pull up to `batch_size` groups from upstream, process them in
    /// parallel, and push the resulting items onto `pending`. Returns
    /// `false` once upstream is fully drained and `pending` is empty.
    fn refill(&mut self) -> bool {
        let batch_size = self.config.batch_size.max(1);
        let mut batch: Vec<OverlappingVariantGroup> = Vec::with_capacity(batch_size);

        while batch.len() < batch_size {
            match self.upstream.next() {
                None => break,
                Some(Err(e)) => {
                    self.pending
                        .push_back(Err(PerGroupMergerError::Upstream(e)));
                    self.done = true;
                    return true;
                }
                Some(Ok(group)) => batch.push(group),
            }
        }

        if batch.is_empty() {
            return false;
        }

        // Collect into a `Vec<Result<_>>` rather than processing
        // inline so the successful prefix in `batch` is emitted in
        // input order before the first error latches `done` below.
        // A short-circuiting fold here would discard already-computed
        // records below the error. Per-chrom outer parallelism lives
        // in `cohort_driver`; this loop runs serially on a single
        // worker thread.
        let results: Vec<Result<Option<MergedRecord>, PerGroupMergerError>> = batch
            .into_iter()
            .map(|group| {
                process_group(
                    group,
                    self.ref_fetcher.as_ref(),
                    &self.config,
                    &self.genotype_tables,
                )
            })
            .collect();

        for result in results {
            match result {
                Ok(Some(record)) => self.pending.push_back(Ok(record)),
                Ok(None) => { /* REF-only group dropped */ }
                Err(e) => {
                    self.pending.push_back(Err(e));
                    self.done = true;
                    return true;
                }
            }
        }

        true
    }
}

impl<I> Iterator for PerGroupMerger<I>
where
    I: Iterator<Item = Result<OverlappingVariantGroup, GrouperError>>,
{
    type Item = Result<MergedRecord, PerGroupMergerError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(item) = self.pending.pop_front() {
                return Some(item);
            }
            if self.done {
                return None;
            }
            if !self.refill() {
                self.done = true;
                return None;
            }
        }
    }
}

// ---------------------------------------------------------------------
// Per-group worker (pure function — borrows shared state by reference)
// ---------------------------------------------------------------------

/// Merge a single `OverlappingVariantGroup` into a `MergedRecord` (or
/// `None` if the group reduces to REF-only after unification).
fn process_group(
    group: OverlappingVariantGroup,
    ref_fetcher: &(dyn RefSeqFetcher + Send + Sync),
    config: &PerGroupMergerConfig,
    genotype_tables: &[Vec<Vec<u8>>],
) -> Result<Option<MergedRecord>, PerGroupMergerError> {
    let chrom_id = group.chrom_id;
    let start = group.start;
    let end = group.end;
    // Defensive: a `ploidy = 0` config would make `genotype_order`
    // return `[[]]` and force a divide-by-zero in the chain-broken
    // fallback. Surface it as a typed error rather than producing a
    // degenerate record.
    if config.ploidy == 0 {
        return Err(PerGroupMergerError::RefFetch {
            chrom_id,
            start,
            end,
            source: std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "PerGroupMergerConfig::ploidy must be >= 1",
            ),
        });
    }
    if end < start {
        // `OverlappingVariantGroup` is a `pub` struct with no
        // constructor invariant. The cohort's own grouper guarantees
        // `start <= end`, but a hand-built fixture or alternative
        // upstream could feed an inverted range; surface it as a
        // typed error rather than wrapping `u32`.
        return Err(PerGroupMergerError::RefFetch {
            chrom_id,
            start,
            end,
            source: std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("group end ({end}) precedes start ({start})"),
            ),
        });
    }
    let span = end - start + 1;

    let ref_seq = ref_fetcher.fetch(chrom_id, start, span).map_err(|source| {
        PerGroupMergerError::RefFetch {
            chrom_id,
            start,
            end,
            source,
        }
    })?;
    if ref_seq.len() != span as usize {
        // The trait contract requires exactly `length` bytes; an
        // alternative fetcher implementation that returns fewer
        // would otherwise blow up downstream in
        // `project_local_allele`'s slice arithmetic.
        return Err(PerGroupMergerError::RefFetch {
            chrom_id,
            start,
            end,
            source: std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                format!("fetcher returned {} bytes for span {}", ref_seq.len(), span,),
            ),
        });
    }

    let n_samples = group
        .records
        .first()
        .map(|pp| pp.per_sample.len())
        .unwrap_or(0);
    debug_assert!(
        group
            .records
            .iter()
            .all(|pp| pp.per_sample.len() == n_samples),
        "OverlappingVariantGroup records disagree on per_sample width",
    );

    let unified = unify_alleles(&group, &ref_seq, n_samples)?;
    let unified = enforce_max_alleles(unified, config.max_alleles);

    // REF-only ⇒ no record. Compound rejection (no entry added) and
    // the cap path (drops into `dropped_other`) both leave `alleles`
    // shorter; the OTHER pool never enters as a real allele.
    if unified.alleles.len() < 2 {
        return Ok(None);
    }

    let projection = project_scalars(&group, &unified, n_samples)?;
    let chain_anchor_flags = build_chain_anchor_flags(&unified, n_samples);
    let (log_likelihoods, n_genotypes) = compute_log_likelihoods(
        &unified,
        &projection,
        &chain_anchor_flags,
        &group,
        &LikelihoodContext {
            chrom_id,
            start,
            end,
            ploidy: config.ploidy,
        },
        genotype_tables,
    )?;

    let alleles: Vec<MergedAllele> = unified
        .alleles
        .into_iter()
        .map(|a| MergedAllele {
            seq: a.seq,
            is_compound: a.is_compound,
            constituents: a.constituents,
        })
        .collect();

    Ok(Some(MergedRecord {
        chrom_id,
        start,
        end,
        alleles,
        ploidy: config.ploidy,
        n_samples,
        n_genotypes,
        scalars: projection.scalars,
        other_scalars: projection.other_scalars,
        chain_anchor_flags,
        log_likelihoods,
    }))
}

// ---------------------------------------------------------------------
// Allele unification (Steps 1 + 2)
// ---------------------------------------------------------------------

/// Mutable scratch entry for the merged allele set during unification.
struct UnifiedAllele {
    seq: Vec<u8>,
    is_compound: bool,
    constituents: Vec<CompoundConstituent>,
    /// Per-sample list of `(record_idx, local_allele_idx)` pairs that
    /// project to this allele in the given sample. Populated for
    /// non-compound alleles (one entry per record/sample that
    /// projects to this allele) and for compound alleles in
    /// chain-anchoring samples (one entry per *constituent* —
    /// allowing the scalar-subtraction pass below to find every
    /// constituent contribution that the compound must claim).
    per_sample_sources: Vec<SourceList>,
    /// Per chain-anchoring sample: the number of chain ids in this
    /// sample whose proposal byte sequence equals this allele's
    /// `seq`. Zero for non-compound alleles and for samples that
    /// don't anchor.
    chain_anchor_counts: Vec<u32>,
    /// Cohort-wide count of supporting reads; drives the `max_alleles`
    /// cap subset.
    cohort_count: u64,
    /// Protected from the `max_alleles` cap. REF (`alleles[0]`) and
    /// chain-anchored compounds are protected.
    cap_protected: bool,
}

/// `alleles[0]` is always REF over the group's reference span.
/// Subsequent entries are ALTs in insertion order; cross-record
/// chain-anchored compounds end the list.
struct UnifiedAlleleSet {
    alleles: Vec<UnifiedAllele>,
    /// Pooled per-sample `(record_idx, local_allele_idx)` sources of
    /// alleles dropped by the `max_alleles` cap. Empty if the cap did
    /// not drop anything. Drives the per-sample `other_scalars` in
    /// the emitted record.
    dropped_other: DroppedOther,
}

/// Per-sample pool of cap-dropped allele sources, plus the cohort-wide
/// count those alleles contributed. Read by `project_scalars` to fold
/// into `other_scalars`.
#[derive(Default)]
struct DroppedOther {
    /// Indexed by `sample_idx`; each inner list holds the
    /// `(record_idx, local_allele_idx)` pairs whose scalars must sum
    /// into the per-sample OTHER pool. Length is either `0` (cap
    /// inactive) or `n_samples`.
    per_sample_sources: Vec<SourceList>,
}

impl DroppedOther {
    fn is_empty(&self) -> bool {
        self.per_sample_sources.is_empty()
    }
}

fn unify_alleles(
    group: &OverlappingVariantGroup,
    ref_seq: &[u8],
    n_samples: usize,
) -> Result<UnifiedAlleleSet, PerGroupMergerError> {
    let (mut alleles, mut byte_index) = project_per_position_alleles(group, ref_seq, n_samples);
    admit_compound_candidates(&mut alleles, &mut byte_index, group, ref_seq, n_samples)?;
    Ok(UnifiedAlleleSet {
        alleles,
        dropped_other: DroppedOther::default(),
    })
}

/// Step 1 of unification: seed `alleles[0]` with REF over the group's
/// reference span, then project every per-sample `(record, local
/// allele)` onto the same span and dedupe by byte equality. Returns
/// the in-progress allele vector and the `byte_key → allele_idx`
/// lookup table that step 2 reuses when admitting compound candidates.
fn project_per_position_alleles(
    group: &OverlappingVariantGroup,
    ref_seq: &[u8],
    n_samples: usize,
) -> (Vec<UnifiedAllele>, AHashMap<Vec<u8>, usize>) {
    // Cheap upper bound on distinct alleles in the group: REF + every
    // per-sample (record, local allele) pair. Most workloads cluster
    // around 2-3 unique alleles per group, so a small with_capacity
    // here avoids the first-insert rehash.
    let max_distinct: usize = 1 + group
        .records
        .iter()
        .map(|pp| {
            pp.per_sample
                .iter()
                .flatten()
                .map(|r| r.alleles.len())
                .sum::<usize>()
        })
        .sum::<usize>();
    let mut alleles: Vec<UnifiedAllele> = Vec::with_capacity(max_distinct);
    let mut byte_index: AHashMap<Vec<u8>, usize> = AHashMap::with_capacity(max_distinct);

    alleles.push(UnifiedAllele {
        seq: ref_seq.to_vec(),
        is_compound: false,
        constituents: Vec::new(),
        per_sample_sources: vec![SourceList::new(); n_samples],
        chain_anchor_counts: vec![0; n_samples],
        cohort_count: 0,
        cap_protected: true,
    });
    byte_index.insert(ref_seq.to_vec(), 0);

    for (record_idx, pp) in group.records.iter().enumerate() {
        let local_offset = (pp.pos - group.start) as usize;
        for (sample_idx, slot) in pp.per_sample.iter().enumerate() {
            let Some(rec) = slot else { continue };
            let local_span = rec.ref_span() as usize;
            for (local_allele_idx, allele) in rec.alleles.iter().enumerate() {
                let projected =
                    project_local_allele(ref_seq, local_offset, local_span, &allele.seq);
                let entry_idx = match byte_index.get(&projected) {
                    Some(&idx) => idx,
                    None => {
                        let idx = alleles.len();
                        byte_index.insert(projected.clone(), idx);
                        alleles.push(UnifiedAllele {
                            seq: projected,
                            is_compound: false,
                            constituents: Vec::new(),
                            per_sample_sources: vec![SourceList::new(); n_samples],
                            chain_anchor_counts: vec![0; n_samples],
                            cohort_count: 0,
                            cap_protected: false,
                        });
                        idx
                    }
                };
                alleles[entry_idx].per_sample_sources[sample_idx]
                    .push((record_idx, local_allele_idx));
                alleles[entry_idx].cohort_count += allele.support.num_obs as u64;
            }
        }
    }

    (alleles, byte_index)
}

/// Step 2 of unification: detect cross-record compound candidates
/// from chain-id evidence, project each onto the group span, and
/// admit them into `alleles` (either extending an existing entry whose
/// byte sequence already matches or appending a new entry). Populates
/// per-sample chain-anchor counts and constituent-source lists used
/// by the scalar-subtraction pass in [`project_scalars`].
fn admit_compound_candidates(
    alleles: &mut Vec<UnifiedAllele>,
    byte_index: &mut AHashMap<Vec<u8>, usize>,
    group: &OverlappingVariantGroup,
    ref_seq: &[u8],
    n_samples: usize,
) -> Result<(), PerGroupMergerError> {
    let compounds = detect_compound_candidates(group, n_samples);
    for candidate in compounds {
        let projected =
            project_compound_onto_group(ref_seq, group, &candidate.constituents_per_first_anchor)?;
        let target_idx = match byte_index.get(&projected) {
            Some(&idx) => {
                // The compound's byte sequence already exists in the
                // set, either as a per-position projection or another
                // candidate. Flag the existing entry as compound and
                // record the constituents.
                if !alleles[idx].is_compound {
                    alleles[idx].is_compound = true;
                    alleles[idx].constituents =
                        sort_constituents(&candidate.constituents_per_first_anchor);
                }
                idx
            }
            None => {
                let idx = alleles.len();
                byte_index.insert(projected.clone(), idx);
                alleles.push(UnifiedAllele {
                    seq: projected,
                    is_compound: true,
                    constituents: sort_constituents(&candidate.constituents_per_first_anchor),
                    per_sample_sources: vec![SourceList::new(); n_samples],
                    chain_anchor_counts: vec![0; n_samples],
                    cohort_count: 0,
                    cap_protected: true,
                });
                idx
            }
        };
        // Record per-sample chain-anchor counts and the per-sample
        // constituent sources (for the scalar-subtraction step).
        for (sample_idx, anchor_evidence) in candidate.per_sample.iter().enumerate() {
            alleles[target_idx].chain_anchor_counts[sample_idx] = anchor_evidence.intersection;
            // The cohort count picks up |chain_id ∩| from each
            // anchoring sample (one observation per spanning read).
            alleles[target_idx].cohort_count += anchor_evidence.intersection as u64;
            if anchor_evidence.intersection > 0 {
                alleles[target_idx].per_sample_sources[sample_idx] =
                    anchor_evidence.constituent_sources.clone();
            }
        }
    }
    Ok(())
}

/// Project a single record's local allele onto the group's reference
/// span. The walker's anchor convention guarantees that the local
/// `allele.seq` substitutes contiguously at `[local_offset,
/// local_offset + local_span)` inside `ref_seq` and produces a
/// well-defined merged byte sequence.
fn project_local_allele(
    ref_seq: &[u8],
    local_offset: usize,
    local_span: usize,
    allele_seq: &[u8],
) -> Vec<u8> {
    let mut out = Vec::with_capacity(ref_seq.len() + allele_seq.len().saturating_sub(local_span));
    out.extend_from_slice(&ref_seq[..local_offset]);
    out.extend_from_slice(allele_seq);
    out.extend_from_slice(&ref_seq[local_offset + local_span..]);
    out
}

/// A compound candidate proposed by chain-id evidence in at least one
/// sample. Constituents are `(record_idx, local_allele_idx)` pairs;
/// the lists are sorted by `record_idx` for determinism.
struct CompoundCandidate {
    /// Constituent indices from the lexicographically first
    /// chain-anchoring sample; used both as the canonical
    /// representative on [`MergedAllele::constituents`] and as the
    /// inputs to compound projection (any anchoring sample's bytes
    /// would give the same projection because the candidate is keyed
    /// by the same byte sequence).
    constituents_per_first_anchor: Vec<CompoundConstituent>,
    /// Per-sample info indexed by `sample_idx`.
    per_sample: Vec<CompoundChainAnchorEvidence>,
}

#[derive(Default, Clone)]
struct CompoundChainAnchorEvidence {
    /// |chain-id intersection across the candidate's constituents in
    /// this sample|. Zero if the sample is chain-broken.
    intersection: u32,
    /// `(record_idx, local_allele_idx)` source pairs for this
    /// sample's constituent alleles in this compound. Only populated
    /// when `intersection > 0`. Used by the scalar-subtraction pass.
    constituent_sources: SourceList,
}

/// Enumerate compound candidates by chain-id evidence. A candidate is
/// keyed by the canonical tuple of its constituent
/// `(record_idx, local_allele_idx)` pairs (sorted by `record_idx`).
/// Under the walker invariant each `(record_idx, local_allele_idx)`
/// names a unique allele byte sequence inside its record, so this
/// tuple is 1-1 with the projected compound byte sequence — keying by
/// it side-steps the byte-collision hazard a `Vec<u8>` key would have
/// for pathological inputs where two distinct local alleles share
/// bytes.
fn detect_compound_candidates(
    group: &OverlappingVariantGroup,
    n_samples: usize,
) -> Vec<CompoundCandidate> {
    let mut candidates: BTreeMap<Vec<(usize, usize)>, CompoundCandidate> = BTreeMap::new();

    for sample_idx in 0..n_samples {
        for proposal in build_chain_proposals(group, sample_idx) {
            // Skip degenerate proposals — chain ids touching only one
            // record do not propose a compound.
            if proposal.constituents.len() < 2 {
                continue;
            }
            let key: Vec<(usize, usize)> = proposal
                .constituents
                .iter()
                .map(|c| (c.record_idx, c.local_allele_idx))
                .collect();

            let entry = candidates.entry(key).or_insert_with(|| CompoundCandidate {
                constituents_per_first_anchor: proposal.constituents.clone(),
                per_sample: vec![CompoundChainAnchorEvidence::default(); n_samples],
            });
            entry.per_sample[sample_idx].intersection += 1;
            if entry.per_sample[sample_idx].constituent_sources.is_empty() {
                entry.per_sample[sample_idx].constituent_sources = proposal
                    .constituents
                    .iter()
                    .map(|c| (c.record_idx, c.local_allele_idx))
                    .collect();
            }
        }
    }

    candidates.into_values().collect()
}

/// One chain id's proposed compound, with the constituent
/// `(record_idx, local_allele_idx)` pairs sorted by `record_idx`.
struct ChainProposal {
    constituents: Vec<CompoundConstituent>,
}

/// Walk the group's records and group chain ids by sample. Each chain
/// id touching ≥ 2 distinct non-REF allele observations across
/// records becomes one proposed compound.
fn build_chain_proposals(group: &OverlappingVariantGroup, sample_idx: usize) -> Vec<ChainProposal> {
    // chain_id → Vec<(record_idx, local_allele_idx)> for this sample.
    let mut by_chain: BTreeMap<ChainId, Vec<CompoundConstituent>> = BTreeMap::new();

    for (record_idx, pp) in group.records.iter().enumerate() {
        let Some(rec) = pp.per_sample.get(sample_idx).and_then(|s| s.as_ref()) else {
            continue;
        };
        for (local_allele_idx, allele) in rec.alleles.iter().enumerate() {
            if local_allele_idx == 0 {
                continue; // REF — never participates in a compound
            }
            for &chain_id in &allele.chain_ids {
                let entry = by_chain.entry(chain_id).or_default();
                // Prevent duplicate (same record, same allele) under
                // pathological mocks. Walker invariants make this
                // unreachable on real `.psp` input.
                if !entry
                    .iter()
                    .any(|c| c.record_idx == record_idx && c.local_allele_idx == local_allele_idx)
                {
                    entry.push(CompoundConstituent {
                        record_idx,
                        local_allele_idx,
                    });
                }
            }
        }
    }

    by_chain
        .into_values()
        .filter(|constituents| constituents.len() >= 2)
        .map(|mut constituents| {
            constituents.sort_by_key(|c| (c.record_idx, c.local_allele_idx));
            ChainProposal { constituents }
        })
        .collect()
}

/// Project a compound's constituent alleles onto the group span. The
/// constituents are processed in `record_idx` order; if they do not
/// overlap in reference span (the typical case), each substitutes
/// independently. Overlapping constituents fall back to a
/// last-write-wins composition over the ref bytes, which the chain-id
/// rule already constrained to be biologically consistent in the
/// anchoring read.
fn project_compound_onto_group(
    ref_seq: &[u8],
    group: &OverlappingVariantGroup,
    constituents: &[CompoundConstituent],
) -> Result<Vec<u8>, PerGroupMergerError> {
    if constituents.is_empty() {
        return Ok(ref_seq.to_vec());
    }

    // We process constituents in record-position order. Each
    // constituent's record may have different allele bytes across
    // samples (different sample at the same record can carry
    // different alleles). For the projection we use the constituent's
    // local_allele_idx as named, against *some* sample that has that
    // record-allele pair (any anchoring sample suffices because by
    // construction the candidate's bytes are the same).
    let mut result: Vec<u8> = Vec::new();
    let mut cursor_offset: usize = 0; // index into ref_seq

    // Every caller currently passes a `Vec<CompoundConstituent>`
    // sorted by `record_idx` (see `build_chain_proposals`,
    // `sort_constituents`); the debug-assert keeps that contract
    // visible.
    debug_assert!(
        constituents
            .windows(2)
            .all(|w| w[0].record_idx <= w[1].record_idx),
        "project_compound_onto_group expects constituents sorted by record_idx",
    );

    for c in constituents {
        let pp = &group.records[c.record_idx];
        let local_offset = (pp.pos - group.start) as usize;
        // Find any sample that has this record's allele at
        // local_allele_idx; that sample's record gives us the bytes.
        // The walker invariant guarantees at least one anchoring
        // sample; missing it is a contract violation, not a data
        // condition.
        let allele_bytes = pp
            .per_sample
            .iter()
            .filter_map(|slot| slot.as_ref())
            .find(|rec| c.local_allele_idx < rec.alleles.len())
            .map(|rec| rec.alleles[c.local_allele_idx].seq.clone())
            .ok_or(PerGroupMergerError::MissingCompoundAlleleBytes {
                chrom_id: group.chrom_id,
                start: group.start,
                end: group.end,
                record_idx: c.record_idx,
                local_allele_idx: c.local_allele_idx,
            })?;
        let local_span = pp
            .per_sample
            .iter()
            .filter_map(|slot| slot.as_ref())
            .map(|rec| rec.ref_span() as usize)
            .next()
            .unwrap_or(1);

        if local_offset < cursor_offset {
            // Overlapping constituent — fall back to last-write-wins:
            // truncate any earlier-written bytes that this constituent
            // overlaps before re-appending.
            let overlap = cursor_offset - local_offset;
            result.truncate(result.len().saturating_sub(overlap));
            cursor_offset = local_offset;
        }

        if local_offset > cursor_offset {
            result.extend_from_slice(&ref_seq[cursor_offset..local_offset]);
        }
        result.extend_from_slice(&allele_bytes);
        cursor_offset = local_offset + local_span;
    }

    if cursor_offset < ref_seq.len() {
        result.extend_from_slice(&ref_seq[cursor_offset..]);
    }

    Ok(result)
}

fn sort_constituents(src: &[CompoundConstituent]) -> Vec<CompoundConstituent> {
    let mut out = src.to_vec();
    out.sort_by_key(|c| (c.record_idx, c.local_allele_idx));
    out
}

// ---------------------------------------------------------------------
// max_alleles cap
// ---------------------------------------------------------------------

fn enforce_max_alleles(mut unified: UnifiedAlleleSet, max_alleles: usize) -> UnifiedAlleleSet {
    if unified.alleles.len() <= max_alleles {
        return unified;
    }

    // Partition into protected (REF, compounds) and prunable.
    let mut prunable_indices: Vec<usize> = unified
        .alleles
        .iter()
        .enumerate()
        .filter_map(|(i, a)| if a.cap_protected { None } else { Some(i) })
        .collect();
    prunable_indices.sort_by_key(|&i| std::cmp::Reverse(unified.alleles[i].cohort_count));

    let protected_count = unified.alleles.iter().filter(|a| a.cap_protected).count();
    let budget_for_prunable = max_alleles.saturating_sub(protected_count);

    if prunable_indices.len() <= budget_for_prunable {
        return unified;
    }

    // Drop prunable indices beyond the budget. Iterate high-to-low so
    // the surviving `remove` calls do not shift earlier indices.
    let mut to_remove: Vec<usize> = prunable_indices.split_off(budget_for_prunable);
    to_remove.sort_unstable();
    to_remove.reverse();
    let removed: Vec<UnifiedAllele> = to_remove
        .into_iter()
        .map(|i| unified.alleles.remove(i))
        .collect();

    // Sum the dropped per-sample sources into the typed OTHER pool.
    // No sentinel allele is pushed onto `alleles`; every downstream
    // helper reads `dropped_other` directly.
    let n_samples = unified.alleles[0].per_sample_sources.len();
    let mut pool: Vec<SourceList> = vec![SourceList::new(); n_samples];
    for dropped in removed {
        for (s, sources) in dropped.per_sample_sources.into_iter().enumerate() {
            pool[s].extend(sources);
        }
    }
    unified.dropped_other = DroppedOther {
        per_sample_sources: pool,
    };
    unified
}

// ---------------------------------------------------------------------
// Scalar projection (Step 3)
// ---------------------------------------------------------------------

struct ScalarProjection {
    /// Flat row-major scalars buffer, length `n_samples * n_kept`.
    /// Indexed as `scalars[sample_idx * n_kept + allele_idx]`. Covers
    /// the *kept* alleles only (compound and per-position) — the OTHER
    /// bucket lives in `other_scalars`.
    scalars: Vec<AlleleSupportStats>,
    /// Width of one sample row in `scalars`.
    n_kept: usize,
    /// Number of samples this projection covers.
    n_samples: usize,
    /// `other_scalars[sample_idx]` — pooled scalars for cap-dropped
    /// alleles. Zero if the cap did not drop anything.
    other_scalars: Vec<AlleleSupportStats>,
}

fn project_scalars(
    group: &OverlappingVariantGroup,
    unified: &UnifiedAlleleSet,
    n_samples: usize,
) -> Result<ScalarProjection, PerGroupMergerError> {
    let n_kept = unified.alleles.len();
    let mut scalars: Vec<AlleleSupportStats> =
        vec![AlleleSupportStats::default(); n_samples * n_kept];

    sum_per_position_scalars(&mut scalars, n_kept, group, unified);
    let other_scalars = pool_dropped_other_scalars(group, &unified.dropped_other, n_samples);
    project_compound_scalars(&mut scalars, n_kept, group, unified)?;
    subtract_compound_from_constituents(&mut scalars, n_kept, group, unified)?;

    Ok(ScalarProjection {
        scalars,
        n_kept,
        n_samples,
        other_scalars,
    })
}

/// Step 1 of scalar projection: sum every non-compound allele's
/// per-sample sources into `scalars[sample_idx][allele_idx]`. The
/// compound entries are filled in by [`project_compound_scalars`].
fn sum_per_position_scalars(
    scalars: &mut [AlleleSupportStats],
    n_kept: usize,
    group: &OverlappingVariantGroup,
    unified: &UnifiedAlleleSet,
) {
    for (allele_idx, allele) in unified.alleles.iter().enumerate() {
        if allele.is_compound {
            continue;
        }
        for (sample_idx, sources) in allele.per_sample_sources.iter().enumerate() {
            let mut bucket = AlleleSupportStats::default();
            for &(record_idx, local_allele_idx) in sources {
                let pp = &group.records[record_idx];
                let Some(rec) = pp.per_sample.get(sample_idx).and_then(|s| s.as_ref()) else {
                    continue;
                };
                let support = rec.alleles[local_allele_idx].support;
                add_support(&mut bucket, &support);
            }
            scalars[sample_idx * n_kept + allele_idx] = bucket;
        }
    }
}

/// Step 2 of scalar projection: pool the per-sample sources of every
/// allele dropped by the `max_alleles` cap into the per-sample OTHER
/// bucket. Those scalars still contribute to the error-cost term in
/// every genotype's likelihood, so they fold into `other_scalars`
/// here rather than disappearing entirely.
fn pool_dropped_other_scalars(
    group: &OverlappingVariantGroup,
    dropped_other: &DroppedOther,
    n_samples: usize,
) -> Vec<AlleleSupportStats> {
    let mut other_scalars: Vec<AlleleSupportStats> = vec![AlleleSupportStats::default(); n_samples];
    if dropped_other.is_empty() {
        return other_scalars;
    }
    for (sample_idx, sources) in dropped_other.per_sample_sources.iter().enumerate() {
        for &(record_idx, local_allele_idx) in sources {
            let pp = &group.records[record_idx];
            let Some(rec) = pp.per_sample.get(sample_idx).and_then(|s| s.as_ref()) else {
                continue;
            };
            let support = rec.alleles[local_allele_idx].support;
            add_support(&mut other_scalars[sample_idx], &support);
        }
    }
    other_scalars
}

/// Step 3 of scalar projection: build the per-sample compound scalar
/// for every chain-anchored compound allele. Chain-evident samples
/// get `(count = |chain_id ∩|, S_s = count × min over constituents of
/// (q_sum / num_obs))`; bias counts (fwd, placed_left, placed_start)
/// are scaled from any constituent's bias by `|chain_id ∩| / num_obs`.
/// See the per-group merger plan §"Step 3" for the derivation.
fn project_compound_scalars(
    scalars: &mut [AlleleSupportStats],
    n_kept: usize,
    group: &OverlappingVariantGroup,
    unified: &UnifiedAlleleSet,
) -> Result<(), PerGroupMergerError> {
    for (allele_idx, allele) in unified.alleles.iter().enumerate() {
        if !allele.is_compound {
            continue;
        }
        for (sample_idx, &inter) in allele.chain_anchor_counts.iter().enumerate() {
            if inter == 0 {
                continue;
            }
            let sources = &allele.per_sample_sources[sample_idx];
            // Gather per-constituent (num_obs, q_sum, fwd, placed_left,
            // placed_start) for the homogeneous-quality approximation.
            let mut min_mean_q: Option<f64> = None;
            let mut bias_fwd = 0_u64;
            let mut bias_left = 0_u64;
            let mut bias_start = 0_u64;
            let mut bias_basis: u64 = 0;
            for &(record_idx, local_allele_idx) in sources {
                let pp = &group.records[record_idx];
                let Some(rec) = pp.per_sample.get(sample_idx).and_then(|s| s.as_ref()) else {
                    continue;
                };
                let stats = rec.alleles[local_allele_idx].support;
                if stats.num_obs == 0 {
                    return Err(PerGroupMergerError::ZeroObservationConstituent {
                        chrom_id: group.chrom_id,
                        start: group.start,
                        end: group.end,
                        phase: CompoundPhase::QualityGather,
                        sample_idx,
                        allele_idx,
                        record_idx,
                        local_allele_idx,
                    });
                }
                let mean_q = stats.q_sum / (stats.num_obs as f64);
                min_mean_q = Some(match min_mean_q {
                    Some(curr) => curr.min(mean_q),
                    None => mean_q,
                });
                bias_fwd += stats.fwd as u64;
                bias_left += stats.placed_left as u64;
                bias_start += stats.placed_start as u64;
                bias_basis += stats.num_obs as u64;
            }
            let count = inter;
            let mean_q =
                min_mean_q.ok_or(PerGroupMergerError::NoQualityForChainAnchoredCompound {
                    chrom_id: group.chrom_id,
                    start: group.start,
                    end: group.end,
                    sample_idx,
                    allele_idx,
                    inter,
                })?;
            let q_sum = mean_q * (count as f64);
            let scale = if bias_basis == 0 {
                0.0
            } else {
                (count as f64) / (bias_basis as f64)
            };
            scalars[sample_idx * n_kept + allele_idx] = AlleleSupportStats {
                num_obs: count,
                q_sum,
                fwd: ((bias_fwd as f64) * scale).round() as u32,
                placed_left: ((bias_left as f64) * scale).round() as u32,
                placed_start: ((bias_start as f64) * scale).round() as u32,
            };
        }
    }
    Ok(())
}

/// Step 4 of scalar projection: in every chain-anchored sample,
/// subtract the compound's claim from each constituent's per-position
/// scalars. A read that supports the compound also appears under each
/// constituent's local allele in the raw record; without this step the
/// constituent would be double-counted against the compound at the
/// group level. The subtraction is scaled by `|chain_id ∩| / num_obs`
/// and clamped per-field to avoid negative bias counts. Uses
/// [`build_source_index`] once for O(1) constituent lookups across
/// the entire pass.
fn subtract_compound_from_constituents(
    scalars: &mut [AlleleSupportStats],
    n_kept: usize,
    group: &OverlappingVariantGroup,
    unified: &UnifiedAlleleSet,
) -> Result<(), PerGroupMergerError> {
    let source_index = build_source_index(unified);
    for (allele_idx, allele) in unified.alleles.iter().enumerate() {
        if !allele.is_compound {
            continue;
        }
        for (sample_idx, &inter) in allele.chain_anchor_counts.iter().enumerate() {
            if inter == 0 {
                continue;
            }
            for &(record_idx, local_allele_idx) in &allele.per_sample_sources[sample_idx] {
                let pp = &group.records[record_idx];
                let Some(rec) = pp.per_sample.get(sample_idx).and_then(|s| s.as_ref()) else {
                    continue;
                };
                let Some(&constituent_idx) =
                    source_index.get(&(sample_idx, record_idx, local_allele_idx))
                else {
                    continue;
                };
                let support = rec.alleles[local_allele_idx].support;
                if support.num_obs == 0 {
                    return Err(PerGroupMergerError::ZeroObservationConstituent {
                        chrom_id: group.chrom_id,
                        start: group.start,
                        end: group.end,
                        phase: CompoundPhase::ConstituentSubtraction,
                        sample_idx,
                        allele_idx,
                        record_idx,
                        local_allele_idx,
                    });
                }
                let dest_idx = sample_idx * n_kept + constituent_idx;
                let scale = (inter as f64) / (support.num_obs as f64);
                let mut to_subtract = AlleleSupportStats {
                    num_obs: inter.min(support.num_obs),
                    q_sum: support.q_sum * scale.min(1.0),
                    fwd: ((support.fwd as f64) * scale.min(1.0)).round() as u32,
                    placed_left: ((support.placed_left as f64) * scale.min(1.0)).round() as u32,
                    placed_start: ((support.placed_start as f64) * scale.min(1.0)).round() as u32,
                };
                to_subtract.fwd = to_subtract.fwd.min(scalars[dest_idx].fwd);
                to_subtract.placed_left =
                    to_subtract.placed_left.min(scalars[dest_idx].placed_left);
                to_subtract.placed_start =
                    to_subtract.placed_start.min(scalars[dest_idx].placed_start);
                subtract_support(&mut scalars[dest_idx], &to_subtract);
            }
        }
    }
    Ok(())
}

/// Build a `(sample_idx, record_idx, local_allele_idx) → allele_idx`
/// map for every non-compound entry of `unified.alleles`. Used by the
/// compound-subtraction pass to find each constituent's kept per-position
/// allele in O(1).
fn build_source_index(unified: &UnifiedAlleleSet) -> AHashMap<(usize, usize, usize), usize> {
    let expected: usize = unified
        .alleles
        .iter()
        .filter(|a| !a.is_compound)
        .map(|a| {
            a.per_sample_sources
                .iter()
                .map(SmallVec::len)
                .sum::<usize>()
        })
        .sum();
    let mut idx: AHashMap<(usize, usize, usize), usize> = AHashMap::with_capacity(expected);
    for (allele_idx, allele) in unified.alleles.iter().enumerate() {
        if allele.is_compound {
            continue;
        }
        for (sample_idx, sources) in allele.per_sample_sources.iter().enumerate() {
            for &(record_idx, local_allele_idx) in sources {
                idx.insert((sample_idx, record_idx, local_allele_idx), allele_idx);
            }
        }
    }
    idx
}

fn add_support(into: &mut AlleleSupportStats, src: &AlleleSupportStats) {
    // Exhaustive destructure so a new field on `AlleleSupportStats`
    // fails to compile here instead of silently being left out of
    // the merged total.
    let AlleleSupportStats {
        num_obs,
        q_sum,
        fwd,
        placed_left,
        placed_start,
    } = *src;
    into.num_obs = into.num_obs.saturating_add(num_obs);
    into.q_sum += q_sum;
    into.fwd = into.fwd.saturating_add(fwd);
    into.placed_left = into.placed_left.saturating_add(placed_left);
    into.placed_start = into.placed_start.saturating_add(placed_start);
}

fn subtract_support(into: &mut AlleleSupportStats, src: &AlleleSupportStats) {
    // `q_sum` is a sum of `ln(P_err) ≤ 0` (always non-positive); the
    // residual after subtracting a less-negative `src.q_sum` must
    // remain `≤ 0`. The clamp guards against over-subtraction taking
    // the residual *positive*; it must not be `.max(0.0)`, which
    // would erase a perfectly valid negative residual.
    let AlleleSupportStats {
        num_obs,
        q_sum,
        fwd,
        placed_left,
        placed_start,
    } = *src;
    into.num_obs = into.num_obs.saturating_sub(num_obs);
    into.q_sum = (into.q_sum - q_sum).min(0.0);
    into.fwd = into.fwd.saturating_sub(fwd);
    into.placed_left = into.placed_left.saturating_sub(placed_left);
    into.placed_start = into.placed_start.saturating_sub(placed_start);
}

// ---------------------------------------------------------------------
// chain-anchor flag table (Step 4)
// ---------------------------------------------------------------------

/// Build the flat `chain_anchor_flags` table — `out[sample_idx * n_alleles + allele_idx]`
/// is `true` iff the allele is a compound and the sample is chain-broken
/// at it. Single heap allocation instead of `n_samples + 1`.
fn build_chain_anchor_flags(unified: &UnifiedAlleleSet, n_samples: usize) -> Vec<bool> {
    let n_alleles = unified.alleles.len();
    let mut out: Vec<bool> = vec![false; n_samples * n_alleles];
    for (allele_idx, allele) in unified.alleles.iter().enumerate() {
        if !allele.is_compound {
            continue;
        }
        for sample_idx in 0..n_samples {
            if allele.chain_anchor_counts[sample_idx] == 0 {
                // Chain-broken at this compound for this sample.
                // (A sample with no record at any constituent position
                // is treated as chain-broken — the fallback path
                // still works because the per-position likelihood
                // degenerates to zero scalars.)
                out[sample_idx * n_alleles + allele_idx] = true;
            }
        }
    }
    out
}

// ---------------------------------------------------------------------
// Per-sample, per-genotype log-likelihood (Step 5)
// ---------------------------------------------------------------------

/// Locus + ploidy context bundle threaded through likelihood
/// computation so callers don't pass six positional `u32`/`u8`s.
struct LikelihoodContext {
    chrom_id: u32,
    start: u32,
    end: u32,
    ploidy: u8,
}

/// Compute the flat row-major log-likelihood table. Returns
/// `(out, n_genotypes)` where `out` has length `n_samples * n_genotypes`
/// and `out[sample_idx * n_genotypes + genotype_idx]` is the natural
/// log-likelihood for that sample/genotype under `c_s = 0`.
fn compute_log_likelihoods(
    unified: &UnifiedAlleleSet,
    projection: &ScalarProjection,
    chain_anchor_flags: &[bool],
    group: &OverlappingVariantGroup,
    ctx: &LikelihoodContext,
    genotype_tables: &[Vec<Vec<u8>>],
) -> Result<(Vec<f64>, usize), PerGroupMergerError> {
    let n_samples = projection.n_samples;
    let n_kept = projection.n_kept;
    let ploidy = ctx.ploidy;
    let n_alleles = unified.alleles.len();
    // `enforce_max_alleles` caps `n_alleles <= config.max_alleles`, and
    // the merger precomputed `genotype_tables` over that range, so the
    // lookup is bounds-safe. Fallback to `genotype_order` keeps direct
    // callers honest if they pass a cache that doesn't cover the input.
    let fallback;
    let genotypes: &[Vec<u8>] = if n_alleles < genotype_tables.len() {
        &genotype_tables[n_alleles]
    } else {
        fallback = genotype_order(ploidy, n_alleles);
        &fallback
    };
    let n_genotypes = genotypes.len();
    let mut out: Vec<f64> = vec![0.0; n_samples * n_genotypes];

    for sample_idx in 0..n_samples {
        let scalars_row = &projection.scalars[sample_idx * n_kept..sample_idx * n_kept + n_kept];
        let flags_row =
            &chain_anchor_flags[sample_idx * n_alleles..sample_idx * n_alleles + n_alleles];
        let other = &projection.other_scalars[sample_idx];
        let out_row_start = sample_idx * n_genotypes;
        for (g_idx, genotype) in genotypes.iter().enumerate() {
            // Chain-broken-compound fallback check.
            let chain_broken_compound = genotype.iter().find_map(|&a| {
                let a = a as usize;
                if unified.alleles[a].is_compound && flags_row[a] {
                    Some(a)
                } else {
                    None
                }
            });

            let value = if let Some(compound_idx) = chain_broken_compound {
                chain_broken_log_likelihood(
                    &unified.alleles[compound_idx],
                    compound_idx as u8,
                    group,
                    sample_idx,
                    genotype,
                    ploidy,
                )
            } else {
                standard_log_likelihood(scalars_row, other, genotype, n_alleles, ploidy)
            };

            if value.is_nan() {
                return Err(PerGroupMergerError::DegenerateLikelihood {
                    chrom_id: ctx.chrom_id,
                    start: ctx.start,
                    end: ctx.end,
                    sample_idx,
                    genotype_idx: g_idx,
                    kind: DegeneracyKind::NaN,
                });
            }
            if value.is_infinite() && value > 0.0 {
                return Err(PerGroupMergerError::DegenerateLikelihood {
                    chrom_id: ctx.chrom_id,
                    start: ctx.start,
                    end: ctx.end,
                    sample_idx,
                    genotype_idx: g_idx,
                    kind: DegeneracyKind::PositiveInfinity,
                });
            }

            out[out_row_start + g_idx] = value;
        }
    }

    Ok((out, n_genotypes))
}

/// Freebayes' closed-form likelihood for `(sample, genotype)` with no
/// chain-broken compound in the genotype.
fn standard_log_likelihood(
    scalars: &[AlleleSupportStats],
    other_scalars: &AlleleSupportStats,
    genotype: &[u8],
    n_alleles: usize,
    ploidy: u8,
) -> f64 {
    // Scratch is sized by `n_alleles ≤ 64`. The `max_alleles` cap
    // defaults to 6 and would have to be raised an order of magnitude
    // before the bitmask shape needs widening — pleasant to keep the
    // allocation off the heap.
    debug_assert!(
        n_alleles <= MAX_BITMASK_ALLELES,
        "standard_log_likelihood assumes n_alleles ≤ {MAX_BITMASK_ALLELES} (u64 bitmask + \
         fixed-size scratch); got {n_alleles}",
    );

    // Bit `a` set iff allele `a` appears in `genotype`. Replaces the
    // previous `BTreeSet<usize>` for the membership check.
    let g_bits: u64 = genotype.iter().fold(0u64, |acc, &a| acc | (1u64 << a));

    // Per-allele copy count in `genotype`. Replaces the previous
    // `BTreeMap<usize, u32>`. Entries past `n_alleles - 1` stay zero.
    let mut p_counts: [u32; MAX_BITMASK_ALLELES] = [0; MAX_BITMASK_ALLELES];
    for &a in genotype {
        p_counts[a as usize] += 1;
    }

    // Error-cost: sum S_s(a) for a ∉ G plus the OTHER pool.
    let mut log_l = 0.0;
    for (a, stats) in scalars.iter().enumerate().take(n_alleles) {
        if (g_bits >> a) & 1 == 0 {
            log_l += stats.q_sum;
        }
    }
    log_l += other_scalars.q_sum;

    // Multinomial: log(N!) − Σ log(n_i!) + Σ n_i log(p_i), with the
    // `xlogy` convention: 0·log 0 = 0.
    let mut n_total: u64 = 0;
    let mut sum_ln_n_fact = 0.0;
    let mut sum_n_log_p = 0.0;
    for a in 0..n_alleles {
        if (g_bits >> a) & 1 == 0 {
            continue;
        }
        let n_i = scalars[a].num_obs as u64;
        n_total += n_i;
        sum_ln_n_fact += ln_factorial(n_i);
        // PANIC-FREE: `a` is filtered by `g_bits`, which is built from
        // `genotype`; `p_counts[a]` is therefore non-zero by construction.
        let p_i = (p_counts[a] as f64) / (ploidy as f64);
        sum_n_log_p += xlogy(n_i as f64, p_i);
    }

    log_l += ln_factorial(n_total) - sum_ln_n_fact + sum_n_log_p;
    log_l
}

/// Width of the per-allele bitmask used by [`standard_log_likelihood`]
/// and its scratch tables. 64 leaves an order-of-magnitude headroom
/// over the default `max_alleles` cap of 6.
const MAX_BITMASK_ALLELES: usize = 64;

/// Chain-broken-compound fallback: decompose to per-position
/// likelihoods at the compound's constituent positions. Non-constituent
/// positions and other compound alleles are not included — the merged
/// set's whole-group scalars don't apply when the haplotype is being
/// reconstructed per-position.
fn chain_broken_log_likelihood(
    compound: &UnifiedAllele,
    compound_slot: u8,
    group: &OverlappingVariantGroup,
    sample_idx: usize,
    genotype: &[u8],
    ploidy: u8,
) -> f64 {
    // Decompose the genotype at each constituent position: the
    // compound's slot in the haplotype takes the constituent's local
    // allele; non-compound slots in the haplotype keep their
    // per-position projection.
    if compound.constituents.is_empty() {
        return 0.0;
    }

    let n_slots = ploidy as usize;
    let mut log_l = 0.0;
    for constituent in &compound.constituents {
        let pp = &group.records[constituent.record_idx];
        let Some(rec) = pp.per_sample.get(sample_idx).and_then(|s| s.as_ref()) else {
            // No record from this sample at this constituent position
            // ⇒ zero scalars; the per-position likelihood reduces to
            // the multinomial of an empty sample, which is 0.0 in log
            // space.
            continue;
        };
        let n_local_alleles = rec.alleles.len();
        // Stack-resident scratch for the per-position multinomial. The
        // previous shape allocated `Vec<u32>` / `Vec<f64>` / `Vec<u32>`
        // per (sample × genotype × constituent) call; for the
        // chain-broken bench shape that was tens-of-millions of small
        // heap allocations on the allocator-bound critical path.
        // Matches the existing `[u32; MAX_BITMASK_ALLELES]` convention
        // in `standard_log_likelihood`.
        debug_assert!(
            n_local_alleles <= MAX_BITMASK_ALLELES,
            "chain_broken_log_likelihood assumes n_local_alleles ≤ {MAX_BITMASK_ALLELES}; got {n_local_alleles}",
        );
        let mut per_pos_counts: [u32; MAX_BITMASK_ALLELES] = [0; MAX_BITMASK_ALLELES];

        // Decode the genotype at this position: each slot that is
        // `compound_slot` takes the constituent's local allele; other
        // slots take REF (slot index 0) as the sentinel for "not
        // observed in the chain-broken interpretation". This matches
        // the plan's "G_decomposed: C/REF ⇒ at p_i, (Y_i, REF_i)"
        // convention.
        for &slot in genotype {
            let local_idx = if slot == compound_slot {
                constituent.local_allele_idx
            } else {
                0 // REF
            };
            per_pos_counts[local_idx] += 1;
        }

        // Read support stats inline from `rec.alleles[a].support` — no
        // gather-into-scratch step is needed because the per-position
        // loops touch each allele exactly once.

        // Error cost at this position: Σ q_sum over local alleles not
        // in per_pos_counts.
        for (a, &count) in per_pos_counts.iter().enumerate().take(n_local_alleles) {
            if count == 0 {
                log_l += rec.alleles[a].support.q_sum;
            }
        }

        // Multinomial at this position.
        let mut n_total: u64 = 0;
        let mut sum_ln_n_fact = 0.0;
        let mut sum_n_log_p = 0.0;
        for (a, &count) in per_pos_counts.iter().enumerate().take(n_local_alleles) {
            if count == 0 {
                continue;
            }
            let n_i = rec.alleles[a].support.num_obs as u64;
            n_total += n_i;
            sum_ln_n_fact += ln_factorial(n_i);
            let p_i = (count as f64) / (n_slots as f64);
            sum_n_log_p += xlogy(n_i as f64, p_i);
        }
        log_l += ln_factorial(n_total) - sum_ln_n_fact + sum_n_log_p;
    }

    log_l
}

/// `0 · log 0 = 0`. Matches the convention used in
/// `genotype_posteriors.rs` and the architecture spec.
fn xlogy(n: f64, p: f64) -> f64 {
    if n == 0.0 { 0.0 } else { n * p.ln() }
}

/// Size of the precomputed `ln(n!)` table. WGS per-allele observation
/// counts at a single position are typically O(10)–O(100); 1024 leaves
/// an order-of-magnitude headroom before the iterative fallback kicks
/// in. Sized for the hot path, not for the absolute worst case.
const LN_FACTORIAL_TABLE_SIZE: usize = 1024;

/// Precomputed `ln(n!)` for `n ∈ [0, LN_FACTORIAL_TABLE_SIZE)`. Built
/// once on first access via [`std::sync::LazyLock`]; subsequent calls
/// to [`ln_factorial`] become an O(1) array load for any `n` in range.
///
/// Trade-off vs the previous iterative implementation: a one-time
/// ~8 KB heap allocation + ~1024 `f64::ln` calls during the first
/// `ln_factorial` invocation, in exchange for O(1) per-call cost on
/// the per-sample × per-genotype × per-allele triple loop in
/// `compute_log_likelihoods`. For depths exceeding the table size the
/// iterative fallback keeps correctness.
static LN_FACTORIAL_TABLE: std::sync::LazyLock<Vec<f64>> = std::sync::LazyLock::new(|| {
    let mut table = Vec::with_capacity(LN_FACTORIAL_TABLE_SIZE);
    table.push(0.0); // ln(0!) = 0
    let mut acc = 0.0;
    for i in 1..LN_FACTORIAL_TABLE_SIZE {
        acc += (i as f64).ln();
        table.push(acc);
    }
    table
});

/// `ln(n!)` — O(1) for `n < LN_FACTORIAL_TABLE_SIZE` via lookup,
/// O(n − table_size) for larger `n` via continued iterative summation
/// from the last tabled value.
fn ln_factorial(n: u64) -> f64 {
    let n_usize = n as usize;
    if n_usize < LN_FACTORIAL_TABLE_SIZE {
        return LN_FACTORIAL_TABLE[n_usize];
    }
    let mut acc = LN_FACTORIAL_TABLE[LN_FACTORIAL_TABLE_SIZE - 1];
    for i in LN_FACTORIAL_TABLE_SIZE..=n_usize {
        acc += (i as f64).ln();
    }
    acc
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::per_sample_pileup::pileup::{AlleleObservation, AlleleSupportStats, PileupRecord};
    use crate::var_calling::per_position_merger::PerPositionPileups;

    // ---------- mock ref fetcher ----------

    #[derive(Clone)]
    struct MockRef {
        seq: Vec<u8>,
        /// 1-based position of `seq[0]`.
        base_offset: u32,
    }

    impl RefSeqFetcher for MockRef {
        fn fetch(
            &self,
            _chrom_id: u32,
            start_1based: u32,
            length: u32,
        ) -> Result<Vec<u8>, std::io::Error> {
            let start_idx = (start_1based - self.base_offset) as usize;
            let end_idx = start_idx + length as usize;
            if end_idx > self.seq.len() {
                return Err(std::io::Error::other("out of range"));
            }
            Ok(self.seq[start_idx..end_idx].to_vec())
        }
    }

    fn fetcher(seq: &[u8], base_offset: u32) -> SharedRefFetcher {
        Arc::new(MockRef {
            seq: seq.to_vec(),
            base_offset,
        })
    }

    // ---------- fixture builders ----------

    fn stats(num_obs: u32, q_sum: f64) -> AlleleSupportStats {
        AlleleSupportStats::new(num_obs, q_sum, num_obs / 2, 0, 0)
    }

    fn allele(seq: &[u8], support: AlleleSupportStats, chain_ids: &[ChainId]) -> AlleleObservation {
        AlleleObservation::new(seq.to_vec(), support, chain_ids.to_vec())
    }

    fn record(chrom_id: u32, pos: u32, alleles: Vec<AlleleObservation>) -> PileupRecord {
        PileupRecord::new(chrom_id, pos, alleles)
    }

    fn pp_one(
        chrom_id: u32,
        pos: u32,
        n_samples: usize,
        slots: Vec<(usize, PileupRecord)>,
    ) -> PerPositionPileups {
        let mut per_sample: Vec<Option<PileupRecord>> = (0..n_samples).map(|_| None).collect();
        for (slot, rec) in slots {
            per_sample[slot] = Some(rec);
        }
        PerPositionPileups {
            chrom_id,
            pos,
            per_sample,
        }
    }

    fn ok_iter(
        groups: Vec<OverlappingVariantGroup>,
    ) -> std::vec::IntoIter<Result<OverlappingVariantGroup, GrouperError>> {
        groups
            .into_iter()
            .map(Ok::<_, GrouperError>)
            .collect::<Vec<_>>()
            .into_iter()
    }

    // ---------- basic merging ----------

    #[test]
    fn stage5_simple_snp_merge_across_two_samples() {
        // Sample 0: A→T at p=100. Sample 1: A→T at p=100.
        // Merged set should be {REF=A, ALT=T} (2 alleles).
        let rec_s0 = record(
            0,
            100,
            vec![
                allele(b"A", stats(8, -1.0), &[1]),
                allele(b"T", stats(6, -3.0), &[2]),
            ],
        );
        let rec_s1 = record(
            0,
            100,
            vec![
                allele(b"A", stats(2, -1.0), &[11]),
                allele(b"T", stats(10, -5.0), &[12]),
            ],
        );
        let pp = pp_one(0, 100, 2, vec![(0, rec_s0), (1, rec_s1)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"A", 100),
            PerGroupMergerConfig::default(),
        );
        let out: Vec<_> = merger.collect();
        assert_eq!(out.len(), 1);
        let record = out.into_iter().next().unwrap().unwrap();
        assert_eq!(record.alleles.len(), 2);
        assert_eq!(record.alleles[0].seq, b"A");
        assert_eq!(record.alleles[1].seq, b"T");
        assert_eq!(record.n_samples, 2);
        assert_eq!(record.scalars_row(0).len(), 2);
        // diploid biallelic ⇒ 3 genotypes
        assert_eq!(record.log_likelihoods_row(0).len(), 3);
    }

    #[test]
    fn stage5_simple_insertion_across_samples() {
        // Sample 0: REF=A, ALT=AT (insertion). Sample 1: REF=A only.
        let rec_s0 = record(
            0,
            100,
            vec![
                allele(b"A", stats(2, -1.0), &[]),
                allele(b"AT", stats(6, -2.0), &[1]),
            ],
        );
        let rec_s1 = record(0, 100, vec![allele(b"A", stats(8, -0.5), &[])]);
        let pp = pp_one(0, 100, 2, vec![(0, rec_s0), (1, rec_s1)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"A", 100),
            PerGroupMergerConfig::default(),
        );
        let out: Vec<_> = merger.collect();
        let record = out.into_iter().next().unwrap().unwrap();
        assert_eq!(record.alleles.len(), 2);
        assert!(record.alleles.iter().any(|a| a.seq == b"AT"));
    }

    #[test]
    fn stage5_simple_deletion_across_samples() {
        // Sample 0: REF=AT, ALT=A (1 bp deletion). Sample 1: only REF=A.
        // ref bases for [100, 101] = "AT".
        let rec_s0 = record(
            0,
            100,
            vec![
                allele(b"AT", stats(2, -1.0), &[]),
                allele(b"A", stats(8, -3.0), &[1]),
            ],
        );
        // Sample 1 covers p=100 only; ref_span=1, REF=A.
        let rec_s1 = record(0, 100, vec![allele(b"A", stats(10, -0.5), &[])]);
        let pp = pp_one(0, 100, 2, vec![(0, rec_s0), (1, rec_s1)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 101,
            records: vec![pp],
        };
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"AT", 100),
            PerGroupMergerConfig::default(),
        );
        let record = merger
            .collect::<Vec<_>>()
            .into_iter()
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(record.alleles[0].seq, b"AT"); // REF over group span
        assert!(record.alleles.iter().any(|a| a.seq == b"A")); // 1 bp deletion
    }

    #[test]
    fn stage5_pure_ref_only_group_drops() {
        // Every sample has only REF allele after unification.
        let rec_s0 = record(0, 100, vec![allele(b"A", stats(8, -1.0), &[])]);
        let pp = pp_one(0, 100, 1, vec![(0, rec_s0)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"A", 100),
            PerGroupMergerConfig::default(),
        );
        assert!(merger.collect::<Vec<_>>().is_empty());
    }

    #[test]
    fn stage5_sample_with_no_record_at_some_position() {
        // Two-position group; sample 1 has no record at the second
        // position. The scalars table must still be the right shape
        // and the likelihood must not panic.
        let rec_s0 = record(
            0,
            100,
            vec![
                allele(b"A", stats(0, 0.0), &[]),
                allele(b"T", stats(8, -2.0), &[1]),
            ],
        );
        let rec_s1 = record(0, 100, vec![allele(b"A", stats(8, -0.5), &[])]);
        let pp0 = pp_one(0, 100, 2, vec![(0, rec_s0), (1, rec_s1)]);
        let rec_s0_p101 = record(
            0,
            101,
            vec![
                allele(b"C", stats(4, -1.0), &[]),
                allele(b"G", stats(5, -2.0), &[2]),
            ],
        );
        let pp1 = pp_one(0, 101, 2, vec![(0, rec_s0_p101)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 101,
            records: vec![pp0, pp1],
        };
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"AC", 100),
            PerGroupMergerConfig::default(),
        );
        let record = merger
            .collect::<Vec<_>>()
            .into_iter()
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(record.n_samples, 2);
        // Sample 1 has no record at p=101 ⇒ its merged scalars for
        // anything that came from p=101 must be zeroed.
        assert!(record.scalars_row(1).iter().all(|s| s.num_obs <= 8));
    }

    // ---------- chain-anchored compounds ----------

    #[test]
    fn stage5_two_overlapping_dels_same_sample_chain_anchored() {
        // Sample 0 has a SNP at p=100 (A→T) and a SNP at p=102 (G→C)
        // sharing chain id 42 — one read supports both.
        let rec_s0_p100 = record(
            0,
            100,
            vec![
                allele(b"A", stats(2, -0.5), &[]),
                allele(b"T", stats(6, -2.0), &[42, 100]),
            ],
        );
        let rec_s0_p102 = record(
            0,
            102,
            vec![
                allele(b"G", stats(2, -0.5), &[]),
                allele(b"C", stats(4, -1.5), &[42, 200]),
            ],
        );
        let pp0 = pp_one(0, 100, 1, vec![(0, rec_s0_p100)]);
        let pp1 = pp_one(0, 102, 1, vec![(0, rec_s0_p102)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 102,
            records: vec![pp0, pp1],
        };
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"AAG", 100),
            PerGroupMergerConfig::default(),
        );
        let record = merger
            .collect::<Vec<_>>()
            .into_iter()
            .next()
            .unwrap()
            .unwrap();
        // Merged set: REF=AAG, T at p=100 (TAG), C at p=102 (AAC),
        // compound (TAC). 4 alleles.
        assert!(
            record
                .alleles
                .iter()
                .any(|a| a.seq == b"TAC" && a.is_compound)
        );
        assert_eq!(record.alleles.iter().filter(|a| a.is_compound).count(), 1);
        // Sample is chain-evident at the compound ⇒ ca_flag false.
        let compound_idx = record.alleles.iter().position(|a| a.is_compound).unwrap();
        assert!(!record.chain_anchor_flag(0, compound_idx));
    }

    #[test]
    fn stage5_no_chain_anchor_anywhere_constituents_only() {
        // Two SNPs at different positions in two samples, but no
        // single read in any sample carries both ⇒ compound rejected.
        let rec_s0_p100 = record(
            0,
            100,
            vec![
                allele(b"A", stats(2, -0.5), &[]),
                allele(b"T", stats(6, -2.0), &[1]),
            ],
        );
        let rec_s0_p102 = record(
            0,
            102,
            vec![
                allele(b"G", stats(2, -0.5), &[]),
                allele(b"C", stats(4, -1.5), &[2]),
            ],
        );
        let pp0 = pp_one(0, 100, 1, vec![(0, rec_s0_p100)]);
        let pp1 = pp_one(0, 102, 1, vec![(0, rec_s0_p102)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 102,
            records: vec![pp0, pp1],
        };
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"AAG", 100),
            PerGroupMergerConfig::default(),
        );
        let record = merger
            .collect::<Vec<_>>()
            .into_iter()
            .next()
            .unwrap()
            .unwrap();
        // No compound entered the merged set.
        assert!(record.alleles.iter().all(|a| !a.is_compound));
        // 3 alleles: REF, T-at-p100, C-at-p102.
        assert_eq!(record.alleles.len(), 3);
    }

    #[test]
    fn stage5_chain_broken_sample_uses_fallback() {
        // Two samples; both have the same two SNPs at p=100 and
        // p=102. Sample 0 has a single read linking them (chain id
        // 42). Sample 1 has them on independent reads. Compound is
        // chain-anchored by sample 0; sample 1 is chain-broken ⇒
        // chain_anchor_flags[1][compound] = true.
        let rec_s0_p100 = record(
            0,
            100,
            vec![
                allele(b"A", stats(2, -0.5), &[]),
                allele(b"T", stats(6, -2.0), &[42]),
            ],
        );
        let rec_s0_p102 = record(
            0,
            102,
            vec![
                allele(b"G", stats(2, -0.5), &[]),
                allele(b"C", stats(4, -1.5), &[42]),
            ],
        );
        let rec_s1_p100 = record(
            0,
            100,
            vec![
                allele(b"A", stats(3, -0.5), &[]),
                allele(b"T", stats(5, -2.0), &[101]),
            ],
        );
        let rec_s1_p102 = record(
            0,
            102,
            vec![
                allele(b"G", stats(3, -0.5), &[]),
                allele(b"C", stats(4, -1.5), &[202]),
            ],
        );
        let pp0 = pp_one(0, 100, 2, vec![(0, rec_s0_p100), (1, rec_s1_p100)]);
        let pp1 = pp_one(0, 102, 2, vec![(0, rec_s0_p102), (1, rec_s1_p102)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 102,
            records: vec![pp0, pp1],
        };
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"AAG", 100),
            PerGroupMergerConfig::default(),
        );
        let record = merger
            .collect::<Vec<_>>()
            .into_iter()
            .next()
            .unwrap()
            .unwrap();
        let compound_idx = record
            .alleles
            .iter()
            .position(|a| a.is_compound)
            .expect("compound is chain-anchored");
        assert!(
            !record.chain_anchor_flag(0, compound_idx),
            "sample 0 anchors"
        );
        assert!(
            record.chain_anchor_flag(1, compound_idx),
            "sample 1 is chain-broken"
        );
    }

    // ---------- genotype enumeration ----------

    #[test]
    fn genotype_order_diploid_three_alleles_matches_vcf_order() {
        let order = genotype_order(2, 3);
        // VCF PL order: 0/0, 0/1, 1/1, 0/2, 1/2, 2/2.
        let expected: Vec<Vec<u8>> = vec![
            vec![0, 0],
            vec![0, 1],
            vec![1, 1],
            vec![0, 2],
            vec![1, 2],
            vec![2, 2],
        ];
        assert_eq!(order, expected);
    }

    #[test]
    fn genotype_order_tetraploid_biallelic_count() {
        let order = genotype_order(4, 2);
        assert_eq!(order.len(), 5); // AAAA, AAAB, AABB, ABBB, BBBB
        assert_eq!(order[0], vec![0, 0, 0, 0]);
        assert_eq!(order[4], vec![1, 1, 1, 1]);
    }

    // ---------- max_alleles cap ----------

    #[test]
    fn max_alleles_cap_drops_lowest_count_alleles() {
        // 5 alts (so 6 alleles total) at p=100; cap = 4 ⇒ 3 alts
        // dropped (only 1 alt kept) — REF protected.
        let alleles = vec![
            allele(b"A", stats(10, -1.0), &[]), // REF
            allele(b"T", stats(1, -1.0), &[1]),
            allele(b"G", stats(2, -1.0), &[2]),
            allele(b"C", stats(3, -1.0), &[3]),
            allele(b"N", stats(4, -1.0), &[4]),
            allele(b"M", stats(20, -1.0), &[5]), // highest count alt, kept
        ];
        let rec = record(0, 100, alleles);
        let pp = pp_one(0, 100, 1, vec![(0, rec)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let config = PerGroupMergerConfig {
            max_alleles: 4,
            ..Default::default()
        };
        let merger = PerGroupMerger::with_config(ok_iter(vec![group]), fetcher(b"A", 100), config);
        let record = merger
            .collect::<Vec<_>>()
            .into_iter()
            .next()
            .unwrap()
            .unwrap();
        assert!(record.alleles.len() <= 4);
        // The highest-count alt "M" must be retained.
        assert!(record.alleles.iter().any(|a| a.seq == b"M"));
        // OTHER pool must have non-zero scalars (the dropped alleles
        // had counts 1+2+3 = 6 at least).
        assert!(record.other_scalars[0].num_obs >= 1);
    }

    // ---------- pure REF group drop after unification ----------

    #[test]
    fn pure_ref_only_after_unification_drops() {
        // Sample 0 has only REF allele at one position.
        let rec_s0 = record(0, 100, vec![allele(b"A", stats(8, 0.0), &[])]);
        let pp = pp_one(0, 100, 1, vec![(0, rec_s0)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"A", 100),
            PerGroupMergerConfig::default(),
        );
        assert!(merger.collect::<Vec<_>>().is_empty());
    }

    // ---------- xlogy edge case ----------

    #[test]
    fn xlogy_zero_count_zero_prob_is_zero() {
        assert_eq!(xlogy(0.0, 0.0), 0.0);
        assert_eq!(xlogy(0.0, 0.5), 0.0);
    }

    #[test]
    fn xlogy_nonzero_count_zero_prob_is_neg_inf() {
        assert!(xlogy(1.0, 0.0).is_infinite() && xlogy(1.0, 0.0) < 0.0);
    }

    // ---------- error surface ----------

    #[test]
    fn ref_fetch_error_surfaces() {
        // MockRef raises out-of-range error when fetched outside seq.
        let rec_s0 = record(
            0,
            100,
            vec![
                allele(b"A", stats(2, 0.0), &[]),
                allele(b"T", stats(8, -2.0), &[1]),
            ],
        );
        let pp = pp_one(0, 100, 1, vec![(0, rec_s0)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        // Empty fetcher seq → out-of-range error.
        let mut merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"", 100),
            PerGroupMergerConfig::default(),
        );
        match merger.next() {
            Some(Err(PerGroupMergerError::RefFetch { .. })) => {}
            other => panic!("expected RefFetch error, got {other:?}"),
        }
        // Latched.
        assert!(matches!(merger.next(), Some(Err(_)) | None));
    }

    // ---------- compound preserved by max_alleles cap ----------

    #[test]
    fn max_alleles_cap_protects_chain_anchored_compound() {
        // 5 SNPs at p=100 plus a chain-anchored compound across p=100
        // and p=102. Cap = 4 ⇒ the compound must survive even though
        // it has lower cohort count than two of the SNPs.
        let many_alts = vec![
            allele(b"A", stats(10, -1.0), &[]),
            allele(b"T", stats(5, -2.0), &[10]),
            allele(b"G", stats(4, -2.0), &[11]),
            allele(b"C", stats(3, -2.0), &[12]),
            allele(b"N", stats(2, -2.0), &[13]),
            allele(b"M", stats(8, -2.0), &[14]),
        ];
        let rec_s0_p100 = record(0, 100, many_alts);
        let rec_s0_p102 = record(
            0,
            102,
            vec![
                allele(b"G", stats(2, -0.5), &[]),
                allele(b"C", stats(1, -1.5), &[10]), // shares chain id with T at p100
            ],
        );
        // p=100 T is chain id 10; p=102 C is chain id 10 ⇒ they form
        // a chain-anchored compound.
        let pp0 = pp_one(0, 100, 1, vec![(0, rec_s0_p100)]);
        let pp1 = pp_one(0, 102, 1, vec![(0, rec_s0_p102)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 102,
            records: vec![pp0, pp1],
        };
        let config = PerGroupMergerConfig {
            max_alleles: 4,
            ..Default::default()
        };
        let merger =
            PerGroupMerger::with_config(ok_iter(vec![group]), fetcher(b"AAG", 100), config);
        let record = merger
            .collect::<Vec<_>>()
            .into_iter()
            .next()
            .unwrap()
            .unwrap();
        assert!(record.alleles.len() <= 4);
        assert!(
            record.alleles.iter().any(|a| a.is_compound),
            "compound must be retained under the cap"
        );
    }

    // ---------- ploidy variation ----------

    #[test]
    fn triploid_genotype_count() {
        // 3 alleles, ploidy 3 ⇒ C(3+3-1, 3-1) = C(5, 2) = 10
        // genotypes.
        let order = genotype_order(3, 3);
        assert_eq!(order.len(), 10);
    }

    #[test]
    fn diploid_biallelic_likelihood_is_finite() {
        // Plain biallelic SNP. Every L should be a finite f64 (modulo
        // a contradicting genotype giving -inf, which is the only
        // allowed infinity).
        let rec_s0 = record(
            0,
            100,
            vec![
                allele(b"A", stats(0, 0.0), &[]),
                allele(b"T", stats(10, -5.0), &[1]),
            ],
        );
        let pp = pp_one(0, 100, 1, vec![(0, rec_s0)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"A", 100),
            PerGroupMergerConfig::default(),
        );
        let record = merger
            .collect::<Vec<_>>()
            .into_iter()
            .next()
            .unwrap()
            .unwrap();
        for ll in record.log_likelihoods_row(0) {
            assert!(
                !(ll.is_nan() || (ll.is_infinite() && *ll > 0.0)),
                "likelihood is {ll}, expected finite or -inf",
            );
        }
    }

    // ---------- multi-group emit shape ----------

    #[test]
    fn two_groups_two_records_emitted() {
        let rec_a = record(
            0,
            100,
            vec![
                allele(b"A", stats(1, 0.0), &[]),
                allele(b"T", stats(5, -1.0), &[1]),
            ],
        );
        let rec_b = record(
            0,
            200,
            vec![
                allele(b"C", stats(1, 0.0), &[]),
                allele(b"G", stats(5, -1.0), &[1]),
            ],
        );
        let pp_a = pp_one(0, 100, 1, vec![(0, rec_a)]);
        let pp_b = pp_one(0, 200, 1, vec![(0, rec_b)]);
        let group_a = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp_a],
        };
        let group_b = OverlappingVariantGroup {
            chrom_id: 0,
            start: 200,
            end: 200,
            records: vec![pp_b],
        };
        // Two separate fetches happen — pass a fetcher covering both
        // positions.
        let ref_seq: Vec<u8> = (100..=200)
            .map(|p| if p == 200 { b'C' } else { b'A' })
            .collect();
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group_a, group_b]),
            Arc::new(MockRef {
                seq: ref_seq,
                base_offset: 100,
            }),
            PerGroupMergerConfig::default(),
        );
        let records: Vec<_> = merger.collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!((records[0].start, records[0].end), (100, 100));
        assert_eq!((records[1].start, records[1].end), (200, 200));
    }

    // ---------- scalar arithmetic helpers ----------

    #[test]
    fn sub_stats_subtraction_preserves_negative_q_sum_residual() {
        // Regression: `subtract_support` previously clamped via `.max(0.0)`,
        // which wiped a legitimate negative residual. `q_sum` is a
        // sum of `ln(P_err) ≤ 0`; subtracting a less-negative `src`
        // from a more-negative `into` must produce a more-negative
        // residual, capped at 0 (over-subtraction).
        let mut into = AlleleSupportStats::new(10, -50.0, 5, 2, 1);
        let src = AlleleSupportStats::new(3, -15.0, 1, 0, 0);
        subtract_support(&mut into, &src);
        assert_eq!(into.num_obs, 7);
        assert!(
            (into.q_sum - (-35.0)).abs() < 1e-9,
            "q_sum = {} (expected -35.0)",
            into.q_sum,
        );
        assert_eq!(into.fwd, 4);
        assert_eq!(into.placed_left, 2);
        assert_eq!(into.placed_start, 1);
    }

    #[test]
    fn sub_stats_over_subtraction_clamps_q_sum_to_zero() {
        // If the scaled subtraction over-runs into.q_sum (a bug or a
        // pathological scaler), the clamp must hold the residual at
        // 0, not let it go positive.
        let mut into = AlleleSupportStats::new(5, -10.0, 2, 1, 0);
        let src = AlleleSupportStats::new(5, -25.0, 2, 1, 0);
        subtract_support(&mut into, &src);
        assert_eq!(into.num_obs, 0);
        assert_eq!(into.q_sum, 0.0);
    }

    #[test]
    fn ln_factorial_small_values_via_table() {
        // ln(0!) = ln(1!) = 0; ln(2!) = ln 2; ln(5!) = ln 120.
        assert_eq!(ln_factorial(0), 0.0);
        assert_eq!(ln_factorial(1), 0.0);
        assert!((ln_factorial(2) - 2.0_f64.ln()).abs() < 1e-12);
        assert!((ln_factorial(5) - 120.0_f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn ln_factorial_iterative_fallback_matches_table_path_at_boundary() {
        // Boundary check: ln_factorial(LN_FACTORIAL_TABLE_SIZE) must
        // equal ln_factorial(LN_FACTORIAL_TABLE_SIZE - 1) + ln(N).
        let table_last = ln_factorial((LN_FACTORIAL_TABLE_SIZE - 1) as u64);
        let just_past = ln_factorial(LN_FACTORIAL_TABLE_SIZE as u64);
        let expected = table_last + (LN_FACTORIAL_TABLE_SIZE as f64).ln();
        assert!(
            (just_past - expected).abs() < 1e-9,
            "table/fallback boundary mismatch: {just_past} vs {expected}",
        );
    }

    #[test]
    fn cap_one_returns_none_no_record_emitted() {
        // Regression: with `max_alleles = 1`, REF stays protected,
        // every ALT folds into the ghost OTHER entry; the merged
        // vector is `[REF, ghost]` (length 2) but the *real* allele
        // count is 1. Must emit no record.
        let rec_s0 = record(
            0,
            100,
            vec![
                allele(b"A", stats(10, -1.0), &[]),
                allele(b"T", stats(5, -2.0), &[1]),
                allele(b"G", stats(4, -2.0), &[2]),
            ],
        );
        let pp = pp_one(0, 100, 1, vec![(0, rec_s0)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let config = PerGroupMergerConfig {
            ploidy: DEFAULT_PLOIDY,
            max_alleles: 1,
            batch_size: DEFAULT_BATCH_SIZE,
        };
        let merger = PerGroupMerger::with_config(ok_iter(vec![group]), fetcher(b"A", 100), config);
        let out: Vec<_> = merger.collect::<Result<Vec<_>, _>>().unwrap();
        assert!(
            out.is_empty(),
            "max_alleles=1 must drop the group, got {} records",
            out.len(),
        );
    }

    #[test]
    fn cap_zero_returns_none_no_record_emitted() {
        let rec_s0 = record(
            0,
            100,
            vec![
                allele(b"A", stats(10, -1.0), &[]),
                allele(b"T", stats(5, -2.0), &[1]),
            ],
        );
        let pp = pp_one(0, 100, 1, vec![(0, rec_s0)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let config = PerGroupMergerConfig {
            ploidy: DEFAULT_PLOIDY,
            max_alleles: 0,
            batch_size: DEFAULT_BATCH_SIZE,
        };
        let merger = PerGroupMerger::with_config(ok_iter(vec![group]), fetcher(b"A", 100), config);
        let out: Vec<_> = merger.collect::<Result<Vec<_>, _>>().unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn add_stats_saturates_on_overflow() {
        let mut into = AlleleSupportStats::new(u32::MAX - 1, -1.0, 0, 0, 0);
        let src = AlleleSupportStats::new(10, -2.0, 0, 0, 0);
        add_support(&mut into, &src);
        assert_eq!(into.num_obs, u32::MAX, "saturating add expected");
        assert!((into.q_sum - (-3.0)).abs() < 1e-9);
    }

    // ---------- defensive boundary tests ----------

    #[test]
    fn process_group_returns_error_on_inverted_start_end() {
        let rec_s0 = record(0, 200, vec![allele(b"A", stats(1, 0.0), &[])]);
        let pp = pp_one(0, 200, 1, vec![(0, rec_s0)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 200,
            end: 100,
            records: vec![pp],
        };
        let mut merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"A", 100),
            PerGroupMergerConfig::default(),
        );
        match merger.next() {
            Some(Err(PerGroupMergerError::RefFetch { .. })) => {}
            other => panic!("expected RefFetch error on inverted span, got {other:?}"),
        }
    }

    #[test]
    fn process_group_returns_error_on_short_fetcher_return() {
        struct ShortRef;
        impl RefSeqFetcher for ShortRef {
            fn fetch(
                &self,
                _chrom_id: u32,
                _start_1based: u32,
                _length: u32,
            ) -> Result<Vec<u8>, std::io::Error> {
                Ok(vec![b'A'])
            }
        }
        let rec_s0 = record(
            0,
            100,
            vec![
                allele(b"AT", stats(1, 0.0), &[]),
                allele(b"A", stats(5, -1.0), &[1]),
            ],
        );
        let pp = pp_one(0, 100, 1, vec![(0, rec_s0)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 101,
            records: vec![pp],
        };
        let mut merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            Arc::new(ShortRef),
            PerGroupMergerConfig::default(),
        );
        match merger.next() {
            Some(Err(PerGroupMergerError::RefFetch { .. })) => {}
            other => panic!("expected RefFetch error on short fetcher, got {other:?}"),
        }
    }

    #[test]
    fn process_group_returns_error_on_zero_ploidy() {
        let rec_s0 = record(
            0,
            100,
            vec![
                allele(b"A", stats(1, 0.0), &[]),
                allele(b"T", stats(5, -1.0), &[1]),
            ],
        );
        let pp = pp_one(0, 100, 1, vec![(0, rec_s0)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let config = PerGroupMergerConfig {
            ploidy: 0,
            max_alleles: DEFAULT_MAX_ALLELES_PER_RECORD,
            batch_size: DEFAULT_BATCH_SIZE,
        };
        let mut merger =
            PerGroupMerger::with_config(ok_iter(vec![group]), fetcher(b"A", 100), config);
        match merger.next() {
            Some(Err(PerGroupMergerError::RefFetch { .. })) => {}
            other => panic!("expected error on ploidy=0, got {other:?}"),
        }
    }

    #[test]
    fn empty_records_group_drops() {
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![],
        };
        let merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"A", 100),
            PerGroupMergerConfig::default(),
        );
        assert!(merger.collect::<Result<Vec<_>, _>>().unwrap().is_empty());
    }

    // ---------- M4 contract-violation errors ----------

    #[test]
    fn project_compound_onto_group_errors_when_no_sample_carries_allele() {
        // Walker invariant says every admitted compound constituent
        // has at least one anchoring sample providing the bytes.
        // Hand-build a constituent referencing a record where the
        // requested local_allele_idx doesn't exist on any sample.
        let rec_s0 = record(0, 100, vec![allele(b"A", stats(2, -0.5), &[])]);
        let pp = pp_one(0, 100, 1, vec![(0, rec_s0)]);
        let group = OverlappingVariantGroup {
            chrom_id: 7,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let constituents = vec![CompoundConstituent {
            record_idx: 0,
            local_allele_idx: 5, // out of range — no sample carries it
        }];
        match project_compound_onto_group(b"A", &group, &constituents) {
            Err(PerGroupMergerError::MissingCompoundAlleleBytes {
                chrom_id,
                record_idx,
                local_allele_idx,
                ..
            }) => {
                assert_eq!(chrom_id, 7);
                assert_eq!(record_idx, 0);
                assert_eq!(local_allele_idx, 5);
            }
            Ok(_) => panic!("expected MissingCompoundAlleleBytes, got Ok"),
            Err(other) => panic!("expected MissingCompoundAlleleBytes, got {other:?}"),
        }
    }

    #[test]
    fn process_group_errors_on_zero_obs_constituent_in_quality_gather() {
        // Two chain-anchored SNPs in one sample — should admit a
        // compound. The p=102 constituent has num_obs = 0 (walker
        // invariant violation). project_scalars must surface the
        // contradiction rather than fabricating an infinite-quality
        // compound.
        let rec_s0_p100 = record(
            0,
            100,
            vec![
                allele(b"A", stats(2, -0.5), &[]),
                allele(b"T", stats(6, -2.0), &[42]),
            ],
        );
        let rec_s0_p102 = record(
            0,
            102,
            vec![
                allele(b"G", stats(2, -0.5), &[]),
                allele(b"C", stats(0, 0.0), &[42]), // zero obs — bogus
            ],
        );
        let pp0 = pp_one(0, 100, 1, vec![(0, rec_s0_p100)]);
        let pp1 = pp_one(0, 102, 1, vec![(0, rec_s0_p102)]);
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 102,
            records: vec![pp0, pp1],
        };
        let mut merger = PerGroupMerger::with_config(
            ok_iter(vec![group]),
            fetcher(b"AAG", 100),
            PerGroupMergerConfig::default(),
        );
        match merger.next() {
            Some(Err(PerGroupMergerError::ZeroObservationConstituent {
                phase: CompoundPhase::QualityGather,
                ..
            })) => {}
            other => panic!("expected ZeroObservationConstituent(QualityGather), got {other:?}"),
        }
    }

    #[test]
    fn project_scalars_errors_when_chain_anchored_compound_has_no_usable_constituent() {
        // inter > 0 for sample 0 but per_sample_sources[0] is empty
        // ⇒ no constituent yields a quality, min_mean_q stays None,
        // the post-loop fallback must surface NoQualityForChainAnchoredCompound.
        let unified = UnifiedAlleleSet {
            alleles: vec![
                UnifiedAllele {
                    seq: b"AAG".to_vec(),
                    is_compound: false,
                    constituents: Vec::new(),
                    per_sample_sources: vec![SourceList::new()],
                    chain_anchor_counts: vec![0],
                    cohort_count: 0,
                    cap_protected: true,
                },
                UnifiedAllele {
                    seq: b"TAC".to_vec(),
                    is_compound: true,
                    constituents: Vec::new(),
                    per_sample_sources: vec![SourceList::new()],
                    chain_anchor_counts: vec![3], // claim 3 chain-anchor obs
                    cohort_count: 3,
                    cap_protected: true,
                },
            ],
            dropped_other: DroppedOther::default(),
        };
        let group = OverlappingVariantGroup {
            chrom_id: 0,
            start: 100,
            end: 102,
            records: vec![pp_one(0, 100, 1, vec![])],
        };
        match project_scalars(&group, &unified, 1) {
            Err(PerGroupMergerError::NoQualityForChainAnchoredCompound {
                sample_idx,
                allele_idx,
                inter,
                ..
            }) => {
                assert_eq!(sample_idx, 0);
                assert_eq!(allele_idx, 1);
                assert_eq!(inter, 3);
            }
            Ok(_) => panic!("expected NoQualityForChainAnchoredCompound, got Ok"),
            Err(other) => panic!("expected NoQualityForChainAnchoredCompound, got {other:?}"),
        }
    }

    // ---------- config / API smoke tests ----------

    #[test]
    fn config_new_accepts_defaults() {
        let cfg = PerGroupMergerConfig::new(
            DEFAULT_PLOIDY,
            DEFAULT_MAX_ALLELES_PER_RECORD,
            DEFAULT_BATCH_SIZE,
        )
        .unwrap();
        assert_eq!(cfg.ploidy, DEFAULT_PLOIDY);
        assert_eq!(cfg.max_alleles, DEFAULT_MAX_ALLELES_PER_RECORD);
        assert_eq!(cfg.batch_size, DEFAULT_BATCH_SIZE);
    }

    #[test]
    fn config_new_accepts_range_boundaries() {
        // Just inside both bounds.
        PerGroupMergerConfig::new(1, 2, 1).unwrap();
        PerGroupMergerConfig::new(MAX_PLOIDY, MAX_ALLELES_PER_VAR_CAP, 1024).unwrap();
    }

    #[test]
    fn config_new_rejects_zero_ploidy() {
        let err = PerGroupMergerConfig::new(0, DEFAULT_MAX_ALLELES_PER_RECORD, DEFAULT_BATCH_SIZE)
            .unwrap_err();
        assert_eq!(err, PerGroupMergerConfigError::InvalidPloidy { got: 0 });
    }

    #[test]
    fn config_new_rejects_too_high_ploidy() {
        let err = PerGroupMergerConfig::new(
            MAX_PLOIDY + 1,
            DEFAULT_MAX_ALLELES_PER_RECORD,
            DEFAULT_BATCH_SIZE,
        )
        .unwrap_err();
        assert_eq!(
            err,
            PerGroupMergerConfigError::InvalidPloidy {
                got: MAX_PLOIDY + 1
            }
        );
    }

    #[test]
    fn config_new_rejects_one_allele() {
        let err = PerGroupMergerConfig::new(DEFAULT_PLOIDY, 1, DEFAULT_BATCH_SIZE).unwrap_err();
        assert_eq!(err, PerGroupMergerConfigError::InvalidMaxAlleles { got: 1 });
    }

    #[test]
    fn config_new_rejects_too_many_alleles() {
        let err = PerGroupMergerConfig::new(
            DEFAULT_PLOIDY,
            MAX_ALLELES_PER_VAR_CAP + 1,
            DEFAULT_BATCH_SIZE,
        )
        .unwrap_err();
        assert_eq!(
            err,
            PerGroupMergerConfigError::InvalidMaxAlleles {
                got: MAX_ALLELES_PER_VAR_CAP + 1
            }
        );
    }

    #[test]
    fn config_new_rejects_zero_batch_size() {
        let err = PerGroupMergerConfig::new(DEFAULT_PLOIDY, DEFAULT_MAX_ALLELES_PER_RECORD, 0)
            .unwrap_err();
        assert_eq!(err, PerGroupMergerConfigError::InvalidBatchSize { got: 0 });
    }
}
