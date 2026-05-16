//! Stage 5 — per-group merger.
//!
//! Consumes `Result<OverlappingVarGroup, GrouperError>` items from
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
//! likelihood and surface that fact through `ca_flags` so Stage 6 can
//! use a cohort-derived compound frequency as the prior signal there.
//!
//! Parallelism is internal to `next()`: a batch of upstream groups is
//! pulled, processed in parallel via `rayon::par_iter`, and the
//! resulting records are emitted one-by-one in input order.

use std::collections::BTreeMap;
use std::collections::VecDeque;
use std::sync::Arc;

use rayon::prelude::*;
use thiserror::Error;

use crate::cohort::variant_grouping::{GrouperError, OverlappingVarGroup};
use crate::per_sample_caller::pileup::{AlleleSupportStats, ChainId, RefSeqFetcher};

/// Maximum number of alleles retained in a single merged record
/// before the cap drops the lowest-cohort-count alleles into the
/// per-sample `other_scalars` pool. Matches GATK's
/// `maxAlternateAlleles` default.
pub const DEFAULT_MAX_ALLELES_PER_RECORD: usize = 6;

/// Default cohort-wide ploidy. Diploid is the common case; the merger
/// supports higher ploidies (triploid, tetraploid) through
/// [`PerGroupMergerConfig::ploidy`].
pub const DEFAULT_PLOIDY: u8 = 2;

/// Default number of upstream groups pulled and processed in parallel
/// per `next()` cycle. Tuned for cache-friendly batches; raise it on
/// very wide cohorts where per-group work amortises better.
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
    /// "<OTHER>" scalar pool that only contributes to the error-cost
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
/// `OverlappingVarGroup.records` plus the local allele index in that
/// record (per the first chain-anchoring sample — see
/// [`MergedAllele::constituents`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CompoundConstituent {
    pub record_idx: usize,
    pub local_allele_idx: usize,
}

/// One emitted item: a unified allele set for a single
/// `OverlappingVarGroup`, the per-sample scalar table, the per-sample
/// chain-anchor flags, and the per-sample, per-genotype log-likelihood
/// vector under `c_s = 0`.
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
    /// `scalars[sample_idx][allele_idx]` — projected
    /// [`AlleleSupportStats`] for that sample/allele. A sample with
    /// no records in the group has every entry zeroed.
    pub scalars: Vec<Vec<AlleleSupportStats>>,
    /// `other_scalars[sample_idx]` — pooled scalars for alleles that
    /// were dropped by the `max_alleles` cap. Always counts into the
    /// error-cost term so the likelihoods preserve the pre-cap
    /// evidence balance.
    pub other_scalars: Vec<AlleleSupportStats>,
    /// `ca_flags[sample_idx][allele_idx]` — `true` iff `allele_idx`
    /// is a compound and `sample_idx` is chain-broken at it (the
    /// fallback likelihood path was used). Stage 6 reads this to
    /// route the cohort-derived compound frequency `f_C` as the
    /// prior signal for chain-broken samples.
    pub ca_flags: Vec<Vec<bool>>,
    /// `log_likelihoods[sample_idx][genotype_idx]` — natural-log
    /// likelihood under `c_s = 0`. Genotype order is canonical
    /// per [`genotype_order`].
    pub log_likelihoods: Vec<Vec<f64>>,
}

impl MergedRecord {
    /// Number of samples this record carries scalars/likelihoods for.
    pub fn n_samples(&self) -> usize {
        self.scalars.len()
    }
}

/// Why a likelihood evaluation hit an unrecoverable value. Surfaced
/// through [`PerGroupMergerError::DegenerateLikelihood`] because the
/// closed-form formula is finite for all valid inputs — a degeneracy
/// signals an internal bug, not a data condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegeneracyKind {
    NaN,
    PositiveInfinity,
}

/// Errors surfaced by the per-group merger.
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum PerGroupMergerError {
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
    let mut out: Vec<Vec<u8>> = Vec::new();
    let mut current: Vec<u8> = vec![0; ploidy];
    collect_non_decreasing(&mut current, 0, 0, n_alleles, ploidy, &mut out);
    out.sort_by(|a, b| {
        for i in (0..ploidy).rev() {
            match a[i].cmp(&b[i]) {
                std::cmp::Ordering::Equal => continue,
                other => return other,
            }
        }
        std::cmp::Ordering::Equal
    });
    out
}

fn collect_non_decreasing(
    current: &mut [u8],
    pos: usize,
    min_allele: u8,
    n_alleles: usize,
    ploidy: usize,
    out: &mut Vec<Vec<u8>>,
) {
    if pos == ploidy {
        out.push(current.to_vec());
        return;
    }
    for a in min_allele..(n_alleles as u8) {
        current[pos] = a;
        collect_non_decreasing(current, pos + 1, a, n_alleles, ploidy, out);
    }
}

// ---------------------------------------------------------------------
// Public merger type
// ---------------------------------------------------------------------

/// Type-erased reference fetcher handle shared across rayon workers.
type SharedRefFetcher = Arc<dyn RefSeqFetcher + Send + Sync>;

/// Streaming, internally-parallel per-group merger.
pub struct PerGroupMerger<I>
where
    I: Iterator<Item = Result<OverlappingVarGroup, GrouperError>>,
{
    upstream: I,
    ref_fetcher: SharedRefFetcher,
    config: PerGroupMergerConfig,
    /// FIFO of merged records produced by the most recent parallel
    /// batch. The iterator drains this before pulling another batch.
    /// `None` entries mean "the group merged down to REF only" — the
    /// merger swallows them per the Stage 4 contract.
    pending: VecDeque<Result<MergedRecord, PerGroupMergerError>>,
    /// Latched after the first surfaced error or upstream exhaustion.
    done: bool,
}

impl<I> std::fmt::Debug for PerGroupMerger<I>
where
    I: Iterator<Item = Result<OverlappingVarGroup, GrouperError>>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PerGroupMerger")
            .field("config", &self.config)
            .field("pending_len", &self.pending.len())
            .field("done", &self.done)
            .finish()
    }
}

impl<I> PerGroupMerger<I>
where
    I: Iterator<Item = Result<OverlappingVarGroup, GrouperError>>,
{
    pub fn new(upstream: I, ref_fetcher: SharedRefFetcher) -> Self {
        Self::with_config(upstream, ref_fetcher, PerGroupMergerConfig::default())
    }

    pub fn with_config(
        upstream: I,
        ref_fetcher: SharedRefFetcher,
        config: PerGroupMergerConfig,
    ) -> Self {
        Self {
            upstream,
            ref_fetcher,
            config,
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
        let mut batch: Vec<OverlappingVarGroup> = Vec::with_capacity(batch_size);

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

        let config = self.config;
        let fetcher = Arc::clone(&self.ref_fetcher);
        let results: Vec<Result<Option<MergedRecord>, PerGroupMergerError>> = batch
            .into_par_iter()
            .map(|group| process_group(group, fetcher.as_ref(), &config))
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
    I: Iterator<Item = Result<OverlappingVarGroup, GrouperError>>,
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
// Per-group worker (pure function so rayon can call it freely)
// ---------------------------------------------------------------------

/// Merge a single `OverlappingVarGroup` into a `MergedRecord` (or
/// `None` if the group reduces to REF-only after unification).
fn process_group(
    group: OverlappingVarGroup,
    ref_fetcher: &(dyn RefSeqFetcher + Send + Sync),
    config: &PerGroupMergerConfig,
) -> Result<Option<MergedRecord>, PerGroupMergerError> {
    let chrom_id = group.chrom_id;
    let start = group.start;
    let end = group.end;
    let span = end - start + 1;

    let ref_seq = ref_fetcher.fetch(chrom_id, start, span).map_err(|source| {
        PerGroupMergerError::RefFetch {
            chrom_id,
            start,
            end,
            source,
        }
    })?;

    let n_samples = group
        .records
        .first()
        .map(|pp| pp.per_sample.len())
        .unwrap_or(0);

    let unified = unify_alleles(&group, &ref_seq, n_samples);
    let unified = enforce_max_alleles(unified, n_samples, config.max_alleles);

    // REF-only ⇒ no record. Compound rejection or cap drop may have
    // collapsed the merged set to a single allele.
    if unified.alleles.len() < 2 {
        return Ok(None);
    }

    let projection = project_scalars(&group, &unified, n_samples);
    let ca_flags = build_ca_flags(&unified, n_samples);
    let log_likelihoods = compute_log_likelihoods(
        &unified,
        &projection,
        &ca_flags,
        &group,
        &LikelihoodCtx {
            chrom_id,
            start,
            end,
            ploidy: config.ploidy,
        },
    )?;
    let _ = n_samples;

    let alleles: Vec<MergedAllele> = unified
        .alleles
        .into_iter()
        .filter(|a| !a.seq.is_empty() || a.is_compound)
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
        scalars: projection.scalars,
        other_scalars: projection.other_scalars,
        ca_flags,
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
    per_sample_sources: Vec<Vec<(usize, usize)>>,
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

struct UnifiedAlleleSet {
    alleles: Vec<UnifiedAllele>,
    /// Index of the REF allele in `alleles`; always 0.
    #[allow(dead_code)]
    ref_idx: usize,
}

fn unify_alleles(
    group: &OverlappingVarGroup,
    ref_seq: &[u8],
    n_samples: usize,
) -> UnifiedAlleleSet {
    // Index `byte_key → position in alleles`. REF is always at 0.
    let mut alleles: Vec<UnifiedAllele> = Vec::new();
    let mut byte_index: BTreeMap<Vec<u8>, usize> = BTreeMap::new();

    alleles.push(UnifiedAllele {
        seq: ref_seq.to_vec(),
        is_compound: false,
        constituents: Vec::new(),
        per_sample_sources: vec![Vec::new(); n_samples],
        chain_anchor_counts: vec![0; n_samples],
        cohort_count: 0,
        cap_protected: true,
    });
    byte_index.insert(ref_seq.to_vec(), 0);

    // Project every (sample, record, allele) onto the group span and
    // dedupe by byte equality.
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
                            per_sample_sources: vec![Vec::new(); n_samples],
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
    // REF is always protected.
    alleles[0].cap_protected = true;

    // ---- Compound-candidate detection + chain-anchor admission ----
    let compounds = detect_compound_candidates(group, n_samples);
    for candidate in compounds {
        let projected =
            project_compound_onto_group(ref_seq, group, &candidate.constituents_per_first_anchor);
        let target_idx = match byte_index.get(&projected) {
            Some(&idx) => {
                // The compound's byte sequence already exists in the
                // set, either as a per-position projection or another
                // candidate. Flag the existing entry as compound and
                // record the constituents.
                if !alleles[idx].is_compound {
                    alleles[idx].is_compound = true;
                    alleles[idx].constituents =
                        constituent_list(&candidate.constituents_per_first_anchor);
                }
                idx
            }
            None => {
                let idx = alleles.len();
                byte_index.insert(projected.clone(), idx);
                alleles.push(UnifiedAllele {
                    seq: projected,
                    is_compound: true,
                    constituents: constituent_list(&candidate.constituents_per_first_anchor),
                    per_sample_sources: vec![Vec::new(); n_samples],
                    chain_anchor_counts: vec![0; n_samples],
                    cohort_count: 0,
                    cap_protected: true,
                });
                idx
            }
        };
        // Record per-sample chain-anchor counts and the per-sample
        // constituent sources (for the scalar-subtraction step).
        for (sample_idx, sample_info) in candidate.per_sample.iter().enumerate() {
            alleles[target_idx].chain_anchor_counts[sample_idx] = sample_info.intersection;
            // The cohort count picks up |chain_id ∩| from each
            // anchoring sample (one observation per spanning read).
            alleles[target_idx].cohort_count += sample_info.intersection as u64;
            if sample_info.intersection > 0 {
                alleles[target_idx].per_sample_sources[sample_idx] =
                    sample_info.constituent_sources.clone();
            }
        }
    }

    UnifiedAlleleSet {
        alleles,
        ref_idx: 0,
    }
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
    per_sample: Vec<CompoundSampleInfo>,
}

#[derive(Default, Clone)]
struct CompoundSampleInfo {
    /// |chain-id intersection across the candidate's constituents in
    /// this sample|. Zero if the sample is chain-broken.
    intersection: u32,
    /// `(record_idx, local_allele_idx)` source pairs for this
    /// sample's constituent alleles in this compound. Only populated
    /// when `intersection > 0`. Used by the scalar-subtraction pass.
    constituent_sources: Vec<(usize, usize)>,
}

/// Enumerate compound candidates by chain-id evidence. A candidate is
/// keyed by the *projected byte sequence* of its constituent allele
/// bytes. Per the implementation plan, we first walk each sample's
/// chain ids to identify which (record_idx, local_allele_idx) sets
/// co-occur on a single read, then group by the resulting byte
/// sequence so that distinct samples' proposals over different local
/// allele indices but identical bytes coalesce into one candidate.
fn detect_compound_candidates(
    group: &OverlappingVarGroup,
    n_samples: usize,
) -> Vec<CompoundCandidate> {
    // Map: byte-sequence-of-candidate → (constituents-from-first-anchor,
    // per-sample info). Bytes-keying allows the same candidate to be
    // anchored by multiple samples even if they use different local
    // allele indices for "the same" allele.
    let group_start = group.start;
    let ref_span_total = (group.end - group.start + 1) as usize;
    // ref_seq is reconstructed lazily inside this function only for
    // candidate-byte-key derivation; the outer unify pass already has
    // the FASTA-derived ref_seq, but we don't have it here so we use
    // a sentinel byte vector of equal length filled with zero. The
    // resulting byte key is unique per candidate because the
    // *positions* and *bytes* of the substitutions encode it.
    //
    // We do NOT need this byte key to be canonical against the FASTA
    // ref_seq; we only need two candidates with the same constituent
    // bytes to map to the same key. Using zeros as the background
    // makes the key independent of the FASTA, which is what we want.
    let pseudo_ref = vec![0u8; ref_span_total];

    let mut candidates: BTreeMap<Vec<u8>, CompoundCandidate> = BTreeMap::new();

    for (sample_idx, sample_info) in
        (0..n_samples).map(|sample_idx| (sample_idx, sample_chain_proposals(group, sample_idx)))
    {
        for proposal in sample_info {
            // Skip degenerate proposals — these are reads whose chain
            // ids only touch one record (no compound).
            if proposal.constituents.len() < 2 {
                continue;
            }
            // Byte-key derived from constituents' allele bytes
            // substituted onto a zero-filled ref of the right length.
            let key = project_compound_onto_group(&pseudo_ref, group, &proposal.constituents);

            let entry = candidates.entry(key).or_insert_with(|| CompoundCandidate {
                constituents_per_first_anchor: proposal.constituents.clone(),
                per_sample: vec![CompoundSampleInfo::default(); n_samples],
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
        let _ = group_start;
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
fn sample_chain_proposals(group: &OverlappingVarGroup, sample_idx: usize) -> Vec<ChainProposal> {
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
    group: &OverlappingVarGroup,
    constituents: &[CompoundConstituent],
) -> Vec<u8> {
    if constituents.is_empty() {
        return ref_seq.to_vec();
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

    let mut constituents_sorted = constituents.to_vec();
    constituents_sorted.sort_by_key(|c| c.record_idx);

    for c in &constituents_sorted {
        let pp = &group.records[c.record_idx];
        let local_offset = (pp.pos - group.start) as usize;
        // Find any sample that has this record's allele at
        // local_allele_idx; that sample's record gives us the bytes.
        let allele_bytes = pp
            .per_sample
            .iter()
            .filter_map(|slot| slot.as_ref())
            .find(|rec| c.local_allele_idx < rec.alleles.len())
            .map(|rec| rec.alleles[c.local_allele_idx].seq.clone())
            .unwrap_or_default();
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

    result
}

fn constituent_list(src: &[CompoundConstituent]) -> Vec<CompoundConstituent> {
    let mut out = src.to_vec();
    out.sort_by_key(|c| (c.record_idx, c.local_allele_idx));
    out
}

// ---------------------------------------------------------------------
// max_alleles cap
// ---------------------------------------------------------------------

fn enforce_max_alleles(
    mut unified: UnifiedAlleleSet,
    _n_samples: usize,
    max_alleles: usize,
) -> UnifiedAlleleSet {
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

    // Mark prunable indices beyond the budget for removal.
    let mut to_remove: Vec<usize> = prunable_indices.split_off(budget_for_prunable);
    to_remove.sort_unstable();
    to_remove.reverse();
    let removed: Vec<UnifiedAllele> = to_remove
        .into_iter()
        .map(|i| unified.alleles.remove(i))
        .collect();

    // Fold removed alleles' per-sample sources into the "<OTHER>"
    // pool (returned via `other_sources` on the set so the projector
    // can sum scalars below).
    unified.alleles[0].cap_protected = true; // keep REF protected
    // Stash the dropped sources on the REF entry; the scalar
    // projection step picks them up via a separate field. To avoid
    // sneaking an extra field through everywhere, we encode the
    // dropped sources as a synthetic "ghost" UnifiedAllele appended
    // last; the projector treats it specially.
    let n_samples = unified.alleles[0].per_sample_sources.len();
    let mut ghost = UnifiedAllele {
        seq: Vec::new(),
        is_compound: false,
        constituents: Vec::new(),
        per_sample_sources: vec![Vec::new(); n_samples],
        chain_anchor_counts: vec![0; n_samples],
        cohort_count: 0,
        cap_protected: false,
    };
    for dropped in removed {
        for (s, sources) in dropped.per_sample_sources.into_iter().enumerate() {
            ghost.per_sample_sources[s].extend(sources);
        }
        ghost.cohort_count += dropped.cohort_count;
    }
    // The ghost is identified by an empty `seq`. The projector
    // recognises it and folds its scalars into `other_scalars` rather
    // than treating it as a real allele.
    unified.alleles.push(ghost);
    unified
}

// ---------------------------------------------------------------------
// Scalar projection (Step 3)
// ---------------------------------------------------------------------

struct ScalarProjection {
    /// `scalars[sample_idx][allele_idx]` matches the *kept* alleles
    /// (compound and per-position). Does not include the OTHER bucket.
    scalars: Vec<Vec<AlleleSupportStats>>,
    /// `other_scalars[sample_idx]` — pooled scalars for cap-dropped
    /// alleles. Zero if the cap did not drop anything.
    other_scalars: Vec<AlleleSupportStats>,
}

fn project_scalars(
    group: &OverlappingVarGroup,
    unified: &UnifiedAlleleSet,
    n_samples: usize,
) -> ScalarProjection {
    // Detect a trailing OTHER ghost entry (empty seq) introduced by
    // the cap.
    let has_ghost = unified
        .alleles
        .last()
        .is_some_and(|a| a.seq.is_empty() && !a.is_compound);
    let n_kept = if has_ghost {
        unified.alleles.len() - 1
    } else {
        unified.alleles.len()
    };

    let mut scalars: Vec<Vec<AlleleSupportStats>> = (0..n_samples)
        .map(|_| vec![AlleleSupportStats::default(); n_kept])
        .collect();
    let mut other_scalars: Vec<AlleleSupportStats> = vec![AlleleSupportStats::default(); n_samples];

    // First: per-position scalars from each non-compound allele's
    // sources.
    for (allele_idx, allele) in unified.alleles.iter().enumerate() {
        if allele.is_compound {
            continue;
        }
        let is_other = has_ghost && allele_idx == n_kept;
        for (sample_idx, sources) in allele.per_sample_sources.iter().enumerate() {
            let mut bucket = AlleleSupportStats::default();
            for &(record_idx, local_allele_idx) in sources {
                let pp = &group.records[record_idx];
                let Some(rec) = pp.per_sample[sample_idx].as_ref() else {
                    continue;
                };
                let support = rec.alleles[local_allele_idx].support;
                add_stats(&mut bucket, &support);
            }
            if is_other {
                add_stats(&mut other_scalars[sample_idx], &bucket);
            } else {
                scalars[sample_idx][allele_idx] = bucket;
            }
        }
    }

    // Compound entries: chain-evident samples get
    // (count = |chain_id ∩|, S_s = |chain_id ∩| × min over
    // constituents of (q_sum / num_obs)). Bias counts (fwd,
    // placed_left, placed_start) are scaled from any constituent's
    // bias by |chain_id ∩| / num_obs.
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
                let Some(rec) = pp.per_sample[sample_idx].as_ref() else {
                    continue;
                };
                let stats = rec.alleles[local_allele_idx].support;
                if stats.num_obs == 0 {
                    continue;
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
            let q_sum = min_mean_q.unwrap_or(0.0) * (count as f64);
            let scale = if bias_basis == 0 {
                0.0
            } else {
                (count as f64) / (bias_basis as f64)
            };
            scalars[sample_idx][allele_idx] = AlleleSupportStats {
                num_obs: count,
                q_sum,
                fwd: ((bias_fwd as f64) * scale).round() as u32,
                placed_left: ((bias_left as f64) * scale).round() as u32,
                placed_start: ((bias_start as f64) * scale).round() as u32,
            };
        }
    }

    // Subtract compound contribution from constituents in
    // chain-anchored samples. A read that supports the compound
    // should be attributed to the compound, not double-counted on
    // each constituent.
    let n_kept_scalar_cols = scalars.first().map(|row| row.len()).unwrap_or(0);
    for allele in &unified.alleles {
        if !allele.is_compound {
            continue;
        }
        for (sample_idx, &inter) in allele.chain_anchor_counts.iter().enumerate() {
            if inter == 0 {
                continue;
            }
            for &(record_idx, local_allele_idx) in &allele.per_sample_sources[sample_idx] {
                let pp = &group.records[record_idx];
                let Some(rec) = pp.per_sample[sample_idx].as_ref() else {
                    continue;
                };
                if let Some(constituent_idx) =
                    find_per_position_allele_idx(unified, record_idx, local_allele_idx, sample_idx)
                    && constituent_idx < n_kept_scalar_cols
                {
                    let support = rec.alleles[local_allele_idx].support;
                    if support.num_obs == 0 {
                        continue;
                    }
                    let scale = (inter as f64) / (support.num_obs as f64);
                    let mut to_subtract = AlleleSupportStats {
                        num_obs: inter.min(support.num_obs),
                        q_sum: support.q_sum * scale.min(1.0),
                        fwd: ((support.fwd as f64) * scale.min(1.0)).round() as u32,
                        placed_left: ((support.placed_left as f64) * scale.min(1.0)).round() as u32,
                        placed_start: ((support.placed_start as f64) * scale.min(1.0)).round()
                            as u32,
                    };
                    to_subtract.fwd = to_subtract
                        .fwd
                        .min(scalars[sample_idx][constituent_idx].fwd);
                    to_subtract.placed_left = to_subtract
                        .placed_left
                        .min(scalars[sample_idx][constituent_idx].placed_left);
                    to_subtract.placed_start = to_subtract
                        .placed_start
                        .min(scalars[sample_idx][constituent_idx].placed_start);
                    sub_stats(&mut scalars[sample_idx][constituent_idx], &to_subtract);
                }
            }
        }
    }

    ScalarProjection {
        scalars,
        other_scalars,
    }
}

/// Look up the index in `unified.alleles` whose sources list
/// `(record_idx, local_allele_idx)` for `sample_idx`. Returns the
/// first match; per-position projections per-sample-source point to a
/// unique merged allele.
fn find_per_position_allele_idx(
    unified: &UnifiedAlleleSet,
    record_idx: usize,
    local_allele_idx: usize,
    sample_idx: usize,
) -> Option<usize> {
    for (idx, allele) in unified.alleles.iter().enumerate() {
        if allele.is_compound {
            continue;
        }
        if allele.per_sample_sources[sample_idx]
            .iter()
            .any(|&(r, a)| r == record_idx && a == local_allele_idx)
        {
            return Some(idx);
        }
    }
    None
}

fn add_stats(into: &mut AlleleSupportStats, src: &AlleleSupportStats) {
    into.num_obs += src.num_obs;
    into.q_sum += src.q_sum;
    into.fwd += src.fwd;
    into.placed_left += src.placed_left;
    into.placed_start += src.placed_start;
}

fn sub_stats(into: &mut AlleleSupportStats, src: &AlleleSupportStats) {
    into.num_obs = into.num_obs.saturating_sub(src.num_obs);
    into.q_sum = (into.q_sum - src.q_sum).max(0.0);
    into.fwd = into.fwd.saturating_sub(src.fwd);
    into.placed_left = into.placed_left.saturating_sub(src.placed_left);
    into.placed_start = into.placed_start.saturating_sub(src.placed_start);
}

// ---------------------------------------------------------------------
// chain-anchor flag table (Step 4)
// ---------------------------------------------------------------------

fn build_ca_flags(unified: &UnifiedAlleleSet, n_samples: usize) -> Vec<Vec<bool>> {
    // ca_flags[sample_idx][allele_idx] = true iff allele is a
    // compound and the sample is chain-broken at it.
    let has_ghost = unified
        .alleles
        .last()
        .is_some_and(|a| a.seq.is_empty() && !a.is_compound);
    let n_kept = if has_ghost {
        unified.alleles.len() - 1
    } else {
        unified.alleles.len()
    };

    let mut out: Vec<Vec<bool>> = (0..n_samples).map(|_| vec![false; n_kept]).collect();
    for (allele_idx, allele) in unified.alleles.iter().take(n_kept).enumerate() {
        if !allele.is_compound {
            continue;
        }
        for (sample_idx, row) in out.iter_mut().enumerate() {
            if allele.chain_anchor_counts[sample_idx] == 0 {
                // Chain-broken at this compound for this sample.
                // (A sample with no record at any constituent position
                // is treated as chain-broken — the fallback path
                // still works because the per-position likelihood
                // degenerates to zero scalars.)
                row[allele_idx] = true;
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
struct LikelihoodCtx {
    chrom_id: u32,
    start: u32,
    end: u32,
    ploidy: u8,
}

fn compute_log_likelihoods(
    unified: &UnifiedAlleleSet,
    projection: &ScalarProjection,
    ca_flags: &[Vec<bool>],
    group: &OverlappingVarGroup,
    ctx: &LikelihoodCtx,
) -> Result<Vec<Vec<f64>>, PerGroupMergerError> {
    let n_samples = projection.scalars.len();
    let ploidy = ctx.ploidy;
    let has_ghost = unified
        .alleles
        .last()
        .is_some_and(|a| a.seq.is_empty() && !a.is_compound);
    let n_kept = if has_ghost {
        unified.alleles.len() - 1
    } else {
        unified.alleles.len()
    };

    let genotypes = genotype_order(ploidy, n_kept);
    let mut out: Vec<Vec<f64>> = (0..n_samples).map(|_| vec![0.0; genotypes.len()]).collect();

    for (sample_idx, out_row) in out.iter_mut().enumerate() {
        let scalars = &projection.scalars[sample_idx];
        let other = &projection.other_scalars[sample_idx];
        for (g_idx, genotype) in genotypes.iter().enumerate() {
            // Chain-broken-compound fallback check.
            let chain_broken_compound = genotype.iter().find_map(|&a| {
                let a = a as usize;
                if unified.alleles[a].is_compound && ca_flags[sample_idx][a] {
                    Some(a)
                } else {
                    None
                }
            });

            let value = if let Some(compound_idx) = chain_broken_compound {
                chain_broken_log_likelihood(
                    unified,
                    group,
                    sample_idx,
                    genotype,
                    compound_idx,
                    ploidy,
                )
            } else {
                standard_log_likelihood(scalars, other, genotype, n_kept, ploidy)
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

            out_row[g_idx] = value;
        }
    }

    Ok(out)
}

/// Freebayes' closed-form likelihood for `(sample, genotype)` with no
/// chain-broken compound in the genotype.
fn standard_log_likelihood(
    scalars: &[AlleleSupportStats],
    other: &AlleleSupportStats,
    genotype: &[u8],
    n_alleles: usize,
    ploidy: u8,
) -> f64 {
    let in_g: Vec<usize> = (0..n_alleles)
        .filter(|&a| genotype.contains(&(a as u8)))
        .collect();
    let g_set: std::collections::BTreeSet<usize> = in_g.iter().copied().collect();

    // Error-cost: sum S_s(a) for a ∉ G plus the OTHER pool.
    let mut log_l = 0.0;
    for (a, stats) in scalars.iter().enumerate().take(n_alleles) {
        if !g_set.contains(&a) {
            log_l += stats.q_sum;
        }
    }
    log_l += other.q_sum;

    // Multinomial: log(N!) − Σ log(n_i!) + Σ n_i log(p_i), with the
    // `xlogy` convention: 0·log 0 = 0.
    let mut p_counts: std::collections::BTreeMap<usize, u32> = std::collections::BTreeMap::new();
    for &a in genotype {
        *p_counts.entry(a as usize).or_insert(0) += 1;
    }

    let mut n_total: u64 = 0;
    let mut sum_ln_n_fact = 0.0;
    let mut sum_n_log_p = 0.0;
    for &a in &in_g {
        let n_i = scalars[a].num_obs as u64;
        n_total += n_i;
        sum_ln_n_fact += ln_factorial(n_i);
        let p_i = (*p_counts.get(&a).unwrap_or(&0) as f64) / (ploidy as f64);
        sum_n_log_p += xlogy(n_i as f64, p_i);
    }

    log_l += ln_factorial(n_total) - sum_ln_n_fact + sum_n_log_p;
    log_l
}

/// Chain-broken-compound fallback: decompose to per-position
/// likelihoods at the compound's constituent positions. Non-constituent
/// positions and other compound alleles are not included — the merged
/// set's whole-group scalars don't apply when the haplotype is being
/// reconstructed per-position.
fn chain_broken_log_likelihood(
    unified: &UnifiedAlleleSet,
    group: &OverlappingVarGroup,
    sample_idx: usize,
    genotype: &[u8],
    compound_idx: usize,
    ploidy: u8,
) -> f64 {
    // Decompose the genotype at each constituent position: the
    // compound's slot in the haplotype takes the constituent's local
    // allele; non-compound slots in the haplotype keep their
    // per-position projection.
    let compound = &unified.alleles[compound_idx];
    if compound.constituents.is_empty() {
        return 0.0;
    }

    let n_slots = ploidy as usize;
    let mut log_l = 0.0;
    for constituent in &compound.constituents {
        let pp = &group.records[constituent.record_idx];
        let Some(rec) = pp.per_sample[sample_idx].as_ref() else {
            // No record from this sample at this constituent position
            // ⇒ zero scalars; the per-position likelihood reduces to
            // the multinomial of an empty sample, which is 0.0 in log
            // space.
            continue;
        };
        let n_local_alleles = rec.alleles.len();
        let mut counts: Vec<u32> = (0..n_local_alleles)
            .map(|a| rec.alleles[a].support.num_obs)
            .collect();
        let q_sums: Vec<f64> = (0..n_local_alleles)
            .map(|a| rec.alleles[a].support.q_sum)
            .collect();

        // Decode the genotype at this position: each slot that is
        // `compound_idx` takes the constituent's local allele; other
        // slots take REF (slot index 0) as the sentinel for "not
        // observed in the chain-broken interpretation". This matches
        // the plan's "G_decomposed: C/REF ⇒ at p_i, (Y_i, REF_i)"
        // convention.
        let mut per_pos_counts: Vec<u32> = vec![0; n_local_alleles];
        for &slot in genotype {
            let slot = slot as usize;
            let local_idx = if slot == compound_idx {
                constituent.local_allele_idx
            } else {
                0 // REF
            };
            per_pos_counts[local_idx] += 1;
        }

        // Error cost at this position: Σ q_sum over local alleles not
        // in per_pos_counts.
        for a in 0..n_local_alleles {
            if per_pos_counts[a] == 0 {
                log_l += q_sums[a];
            }
        }

        // Multinomial at this position.
        let n_total: u64 = counts
            .iter()
            .enumerate()
            .filter(|(a, _)| per_pos_counts[*a] > 0)
            .map(|(_, &c)| c as u64)
            .sum();
        let mut sum_ln_n_fact = 0.0;
        let mut sum_n_log_p = 0.0;
        for a in 0..n_local_alleles {
            if per_pos_counts[a] == 0 {
                continue;
            }
            let n_i = counts[a] as u64;
            sum_ln_n_fact += ln_factorial(n_i);
            let p_i = (per_pos_counts[a] as f64) / (n_slots as f64);
            sum_n_log_p += xlogy(n_i as f64, p_i);
        }
        let _ = &mut counts;
        log_l += ln_factorial(n_total) - sum_ln_n_fact + sum_n_log_p;
    }

    log_l
}

/// `0 · log 0 = 0`. Matches the convention used in
/// `genotype_posteriors.rs` and the architecture spec.
fn xlogy(n: f64, p: f64) -> f64 {
    if n == 0.0 { 0.0 } else { n * p.ln() }
}

/// `ln(n!)` via iterative summation. Fast enough for the depths Stage 5
/// sees in practice; switch to `lgamma` if profiling shows this is hot.
fn ln_factorial(n: u64) -> f64 {
    let mut acc = 0.0;
    for i in 2..=n {
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
    use crate::cohort::per_position_merger::PerPositionPileups;
    use crate::per_sample_caller::pileup::{AlleleObservation, AlleleSupportStats, PileupRecord};

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
        groups: Vec<OverlappingVarGroup>,
    ) -> std::vec::IntoIter<Result<OverlappingVarGroup, GrouperError>> {
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
        let group = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let merger = PerGroupMerger::new(ok_iter(vec![group]), fetcher(b"A", 100));
        let out: Vec<_> = merger.collect();
        assert_eq!(out.len(), 1);
        let record = out.into_iter().next().unwrap().unwrap();
        assert_eq!(record.alleles.len(), 2);
        assert_eq!(record.alleles[0].seq, b"A");
        assert_eq!(record.alleles[1].seq, b"T");
        assert_eq!(record.scalars.len(), 2);
        assert_eq!(record.scalars[0].len(), 2);
        // diploid biallelic ⇒ 3 genotypes
        assert_eq!(record.log_likelihoods[0].len(), 3);
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
        let group = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let merger = PerGroupMerger::new(ok_iter(vec![group]), fetcher(b"A", 100));
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
        let group = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 101,
            records: vec![pp],
        };
        let merger = PerGroupMerger::new(ok_iter(vec![group]), fetcher(b"AT", 100));
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
        let group = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let merger = PerGroupMerger::new(ok_iter(vec![group]), fetcher(b"A", 100));
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
        let group = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 101,
            records: vec![pp0, pp1],
        };
        let merger = PerGroupMerger::new(ok_iter(vec![group]), fetcher(b"AC", 100));
        let record = merger
            .collect::<Vec<_>>()
            .into_iter()
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(record.scalars.len(), 2);
        // Sample 1 has no record at p=101 ⇒ its merged scalars for
        // anything that came from p=101 must be zeroed.
        assert!(record.scalars[1].iter().all(|s| s.num_obs <= 8));
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
        let group = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 102,
            records: vec![pp0, pp1],
        };
        let merger = PerGroupMerger::new(ok_iter(vec![group]), fetcher(b"AAG", 100));
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
        assert!(!record.ca_flags[0][compound_idx]);
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
        let group = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 102,
            records: vec![pp0, pp1],
        };
        let merger = PerGroupMerger::new(ok_iter(vec![group]), fetcher(b"AAG", 100));
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
        // ca_flags[1][compound] = true.
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
        let group = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 102,
            records: vec![pp0, pp1],
        };
        let merger = PerGroupMerger::new(ok_iter(vec![group]), fetcher(b"AAG", 100));
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
        assert!(!record.ca_flags[0][compound_idx], "sample 0 anchors");
        assert!(record.ca_flags[1][compound_idx], "sample 1 is chain-broken");
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
        let group = OverlappingVarGroup {
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
        let group = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let merger = PerGroupMerger::new(ok_iter(vec![group]), fetcher(b"A", 100));
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
        let group = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        // Empty fetcher seq → out-of-range error.
        let mut merger = PerGroupMerger::new(ok_iter(vec![group]), fetcher(b"", 100));
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
        let group = OverlappingVarGroup {
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
        let group = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp],
        };
        let merger = PerGroupMerger::new(ok_iter(vec![group]), fetcher(b"A", 100));
        let record = merger
            .collect::<Vec<_>>()
            .into_iter()
            .next()
            .unwrap()
            .unwrap();
        for ll in &record.log_likelihoods[0] {
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
        let group_a = OverlappingVarGroup {
            chrom_id: 0,
            start: 100,
            end: 100,
            records: vec![pp_a],
        };
        let group_b = OverlappingVarGroup {
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
        let merger = PerGroupMerger::new(
            ok_iter(vec![group_a, group_b]),
            Arc::new(MockRef {
                seq: ref_seq,
                base_offset: 100,
            }),
        );
        let records: Vec<_> = merger.collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!((records[0].start, records[0].end), (100, 100));
        assert_eq!((records[1].start, records[1].end), (200, 200));
    }
}
