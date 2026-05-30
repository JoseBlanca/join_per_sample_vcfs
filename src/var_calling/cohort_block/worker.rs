//! Per-window variant-calling math.
//!
//! Phase A.2 step 2 — column-native end-to-end. The worker runs each
//! partition group through the three columnar kernels
//! ([`unify_alleles_columnar`],
//! [`project_scalars_columnar`],
//! [`compute_log_likelihoods_columnar`]) and feeds the resulting
//! column-native buffers straight into
//! [`run_em_columnar`] through borrowed slices + a [`ColumnarAllelesView`]. No row-shape
//! `MergedRecord` is ever materialised on the production path.
//! `PosteriorRecord` is still row-shape (the VCF writer's input);
//! `Vec<MergedAllele>` is allocated once per group at emit time
//! rather than at EM input time.
//!
//! The Phase A.0 row-shape adapter
//! `build_overlapping_variant_group` is retained `#[cfg(test)]` in
//! `super::test_helpers` (Mi16) as the byte-identity oracle for the
//! kernels' unit tests.

use std::sync::Arc;

use thiserror::Error;

use crate::fasta::fetcher::{ChromRefFetchError, ChromRefFetcher};
use crate::var_calling::cohort_block::columns::MaterialisedChunk;
use crate::var_calling::cohort_block::kernels::compute_log_likelihoods::{
    ComputeLogLikelihoodsError, LogLikelihoodsColumns, compute_log_likelihoods_columnar,
};
use crate::var_calling::cohort_block::kernels::project_scalars::{
    ProjectScalarsError, ProjectScalarsScratch, ProjectedScalarsColumns, project_scalars_columnar,
};
use crate::var_calling::cohort_block::kernels::unify_alleles::{
    UnifiedAllelesColumns, UnifyAllelesError, UnifyAllelesScratch, unify_alleles_columnar,
};
use crate::var_calling::cohort_block::partition::{PartitionScratch, WindowPartition};
use crate::var_calling::per_group_merger::{
    CompoundConstituent, LikelihoodContext, MergedAllele, PerGroupMergerConfig,
    PerGroupMergerError, SharedRefFetcher,
};
use crate::var_calling::posterior_engine::backends::InterpUnivariateSimdMath;
use crate::var_calling::posterior_engine::{
    AllelesView, EmInputs, PosteriorEngineConfig, PosteriorEngineError, PosteriorRecord,
    RecordScratch, run_em_columnar,
};

/// Reusable scratch for the column-native pipeline. The driver owns
/// one of these per chunk-loop instance and threads it into every
/// [`run_window`] call so the per-group columnar buffers retain
/// allocated capacity across windows + chunks. Build empty with
/// [`Self::empty`] and let the first call populate it; subsequent
/// calls only grow buffers when a group exceeds the prior high-water
/// mark.
#[derive(Default)]
pub struct ColumnarPipelineScratch {
    unify: UnifyAllelesScratch,
    unified: UnifiedAllelesColumns,
    project: ProjectScalarsScratch,
    projection: ProjectedScalarsColumns,
    log_likelihoods: LogLikelihoodsColumns,
    /// Per-group chain-anchor flag table —
    /// `chain_anchor_flags[sample × n_alleles + allele]` is `true` iff
    /// the allele is a compound and the sample is chain-broken at it.
    /// Refilled per group; capacity carries across groups.
    ///
    /// Built once and used for two distinct downstream consumers:
    /// 1. EM input validation
    ///    ([`run_em_columnar`](crate::var_calling::posterior_engine::run_em_columnar)
    ///    checks its length).
    /// 2. The emitted [`PosteriorRecord`]'s `chain_anchor_flags`
    ///    field, which the VCF writer reads for the `CA` flag.
    chain_anchor_flags: Vec<bool>,
    /// Per-engine EM scratch — owned at the worker level (rather
    /// than inside a `PosteriorEngine` instance) so the chunk loop
    /// reuses the same EM intermediates across every group.
    record_scratch: RecordScratch,
    /// Math backend instance. `InterpUnivariateSimdMath` is a unit
    /// struct so this holds no data; we keep it as a field for
    /// symmetry with the Phase A.1 design where the backend was a
    /// `PosteriorEngine` parameter.
    math: InterpUnivariateSimdMath,
}

/// M15: per-window run statistics produced by `run_window`. Lives on
/// the `WorkerSlot` as a sibling to `pipeline_scratch` so the driver
/// drains it after each parallel-dispatch round and rolls the
/// per-slot counts into [`ChunkDriverStats`](super::driver::ChunkDriverStats).
///
/// Previously these three counters lived inside
/// `ColumnarPipelineScratch` next to the per-group scratch buffers,
/// mixing two different lifetimes (scratch is per-group; counters are
/// per-`run_window`). Splitting them gives `ColumnarPipelineScratch`
/// a single, scratch-only purpose and lets the driver's stats-rollup
/// code talk to `WindowRunStats` directly instead of through a
/// `take_*` API.
#[derive(Debug, Default, Clone, Copy)]
pub struct WindowRunStats {
    /// Number of groups skipped because
    /// `unified.alleles.len() > cfg.max_alleles_lh_calc`.
    pub lh_cap_groups_skipped: u64,
    /// Sum of `unified.alleles.len()` across every skipped group —
    /// the row-shape merger reports this in `allele_lh_cap_hit:
    /// alleles_total=N`.
    pub lh_cap_alleles_in_skipped: u64,
    /// M23: groups whose post-unification allele set has only REF
    /// (every ALT was absorbed into the OTHER pool by the max-alleles
    /// cap, or none survived the per-position projection dedup). These
    /// are emitted as `Ok(None)` from `build_posterior_record_columnar`
    /// (matching the row-shape merger's REF-only no-op skip).
    pub groups_skipped_post_unify_ref_only: u64,
}

impl WindowRunStats {
    /// Reset every counter to zero while preserving the field shape.
    /// The driver's drain calls `std::mem::take(&mut slot.stats)`
    /// directly; this method exists for tests that want to clear in
    /// place without moving.
    pub fn clear(&mut self) {
        *self = Self::default();
    }
}

/// [`AllelesView`] over [`UnifiedAllelesColumns`]. Phase A.2 step 2:
/// lets the column-native worker feed `run_em_columnar` without
/// allocating a row-shape `Vec<MergedAllele>` for the EM's
/// per-allele queries. The `Vec<MergedAllele>` is still allocated —
/// but only at the end of the group's pipeline, for the emitted
/// [`PosteriorRecord`] that the VCF writer reads.
pub(crate) struct ColumnarAllelesView<'a> {
    columns: &'a UnifiedAllelesColumns,
}

impl<'a> ColumnarAllelesView<'a> {
    pub(crate) fn new(columns: &'a UnifiedAllelesColumns) -> Self {
        Self { columns }
    }
}

impl AllelesView for ColumnarAllelesView<'_> {
    fn len(&self) -> usize {
        self.columns.n_alleles()
    }
    fn ref_len(&self) -> usize {
        // REF is allele 0 by construction; its bytes are the
        // projected-onto-group-span sequence, which equals the
        // group's REF span length.
        self.columns.allele_seq(0).len()
    }
    fn seq_len(&self, allele_idx: usize) -> usize {
        self.columns.allele_seq(allele_idx).len()
    }
    fn is_compound(&self, allele_idx: usize) -> bool {
        self.columns.is_compound[allele_idx]
    }
}

impl ColumnarPipelineScratch {
    pub fn empty() -> Self {
        Self::default()
    }
}

/// Per-block worker scratch: everything needed to partition the block,
/// pre-fetch its REF bytes, run the column-native math, and accumulate
/// the emitted [`PosteriorRecord`]s. The [`BlockIterator`] owns one and
/// reuses it across blocks.
///
/// [`BlockIterator`]: crate::var_calling::cohort_block::driver
pub struct WorkerSlot {
    /// Reusable per-window partition scratch (sample-stride buffers).
    pub partition_scratch: PartitionScratch,
    /// Reusable per-window partition output (CSR group + position
    /// arrays; refilled by
    /// [`partition_window`](super::partition::partition_window)).
    pub partition: WindowPartition,
    /// Pre-fetched REF bytes per group in `partition`. Length is
    /// always `partition.n_groups()` after the driver runs the
    /// [`prefetch_window_ref_bytes`] sequential pre-pass.
    pub pre_fetched_ref_bytes: Vec<Vec<u8>>,
    /// Mi3: column-native pipeline scratch — UnifyAlleles,
    /// ProjectScalars, LogLikelihoods buffers, the EM `RecordScratch`,
    /// etc. Named `pipeline_scratch` to mirror the sibling
    /// `partition_scratch` field; the bare `scratch` name didn't say
    /// *which* scratch this is when the slot already has another.
    pub pipeline_scratch: ColumnarPipelineScratch,
    /// Mi3: per-window emitted records, drained sequentially in
    /// window order after the driver's per-window math finishes. Named
    /// `posterior_records` so the slot drain `for r in slot.posterior_records.drain(..)`
    /// reads as what it is; the previous `output_buf` was a bare
    /// generic name.
    pub posterior_records: Vec<PosteriorRecord>,
    /// M15: per-window run statistics. The driver `std::mem::take`s
    /// this after each `run_window` call and rolls the counts into
    /// `ChunkDriverStats`. Sibling to `pipeline_scratch` (which holds
    /// per-group buffers); the two structs have different lifetimes
    /// (counters are per-window, scratch is per-group), so they live
    /// side by side rather than bundled.
    pub stats: WindowRunStats,
}

impl WorkerSlot {
    pub fn new(n_samples: usize) -> Self {
        Self {
            partition_scratch: PartitionScratch::with_n_samples(n_samples),
            partition: WindowPartition::empty(),
            pre_fetched_ref_bytes: Vec::new(),
            pipeline_scratch: ColumnarPipelineScratch::empty(),
            posterior_records: Vec::new(),
            stats: WindowRunStats::default(),
        }
    }
}

/// Run the worker math on one window's partition.
///
/// **Inputs.**
/// - `chunk`: the loaded `MaterialisedChunk` whose columns the
///   partition's `(sample_idx, row_idx)` handles refer into.
/// - `partition`: the variant-group partition produced by
///   [`partition_window`](super::partition::partition_window) for
///   one of `chunk.windows`. Its `n_groups()` may be 0 — empty
///   partitions are a fast no-op that never reads any REF bytes.
/// - `pre_fetched_ref_bytes`: one entry per variant group in
///   `partition` (`pre_fetched_ref_bytes.len() ==
///   partition.n_groups()`); each entry is the REF bases over
///   `[group_starts[g], group_ends[g]]`. Driver-side pre-fetch
///   lifts the fetcher state out of the per-window math so the
///   worker is `Sync`-friendly for Phase B's parallel dispatch.
/// - `per_group_cfg`: tuning for the column-native kernels
///   (`ploidy`, `max_alleles`, `max_alleles_lh_calc`). The
///   `batch_size` field is ignored — Phase A.1 doesn't batch.
/// - `posterior_cfg`: tuning passed straight to
///   [`PosteriorEngine`](crate::var_calling::posterior_engine::PosteriorEngine).
/// - `scratch`: per-window columnar scratch (
///   [`ColumnarPipelineScratch`]). The driver owns one per worker
///   slot and reuses it across every window in every chunk.
/// - `output`: appended-to record buffer. Cleared on entry so the
///   driver can reuse it across windows.
// Mi11: take `posterior_cfg` by `&` — nothing inside `run_window`
// needs ownership (the function only borrows the config to forward
// into `build_posterior_record_columnar`'s `&PosteriorEngineConfig`
// parameter). The previous by-value shape forced the driver's rayon
// dispatch to `posterior_cfg.clone()` per worker slot, even though
// `PosteriorEngineConfig` is not `Copy` (owns
// `Option<Vec<f64>>` for per-sample fixation overrides).
// `per_group_cfg` stays by value: `PerGroupMergerConfig` is `Copy`.
#[allow(clippy::too_many_arguments)]
pub fn run_window(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    pre_fetched_ref_bytes: &[Vec<u8>],
    per_group_cfg: PerGroupMergerConfig,
    posterior_cfg: &PosteriorEngineConfig,
    scratch: &mut ColumnarPipelineScratch,
    output: &mut Vec<PosteriorRecord>,
    stats: &mut WindowRunStats,
) -> Result<(), PosteriorEngineError> {
    output.clear();
    if partition.n_groups() == 0 {
        return Ok(());
    }
    debug_assert_eq!(
        pre_fetched_ref_bytes.len(),
        partition.n_groups(),
        "one pre-fetched REF buffer required per variant group in the partition",
    );

    for (g, group_ref_bytes) in pre_fetched_ref_bytes.iter().enumerate() {
        if let Some(record) = build_posterior_record_columnar(
            chunk,
            partition,
            g,
            group_ref_bytes,
            &per_group_cfg,
            posterior_cfg,
            scratch,
            stats,
        )? {
            output.push(record);
        }
    }
    Ok(())
}

/// Sequentially populate `out` with one `Vec<u8>` of REF bases per
/// variant group in `partition`. The driver calls this once per
/// window before [`run_window`] so the worker math (parallelised
/// by Phase B's window dispatch) consumes `Sync`-friendly
/// pre-fetched bytes rather than a `!Sync`
/// [`StreamingChromRefFetcher`](crate::fasta::fetcher::StreamingChromRefFetcher).
pub fn prefetch_window_ref_bytes(
    partition: &WindowPartition,
    ref_fetcher: &dyn ChromRefFetcher,
    out: &mut Vec<Vec<u8>>,
) -> Result<(), ChromRefFetchError> {
    // M2: preserve the inner `Vec<u8>` allocations across calls.
    // The previous `out.clear()` shape dropped every inner buffer,
    // making the next pass allocate fresh on every group; this
    // defeated the WorkerSlot scratch-reuse design (one slot per
    // worker × n_groups per window × chunks per chrom).
    //
    // Strategy: shrink the outer Vec to the new `n_groups` without
    // freeing the kept inner buffers (`Vec::truncate` drops the tail
    // — that's expected when the high-water mark drops), then grow
    // with empty buffers if needed, then `clear()` each inner Vec
    // in place and `fetch_into` straight into it. `fetch_into`'s
    // contract is "buffer is cleared first, on success
    // `dst.len() == length`"; subsequent calls reuse the underlying
    // allocator capacity.
    let n_groups = partition.n_groups();
    if out.len() < n_groups {
        out.resize_with(n_groups, Vec::new);
    } else {
        out.truncate(n_groups);
    }

    // M24: pin the invariant that every group's `[start, end]` span
    // is inside the contig the fetcher serves. The walker invariant
    // (record positions ≤ chrom_length, ref_span = 1 at the chrom
    // edge) makes this hold by construction on real input; a
    // regression that emits a group with `group_end > chrom_length`
    // would otherwise surface as a typed `ChromRefFetchError` from
    // `fetch_into` deep on the parallel-dispatch hot path. The
    // debug_assert catches it at the prefetch step under tests, with
    // a focused error message naming the offending group.
    let chrom_length = ref_fetcher.length();
    for (g, bytes) in out.iter_mut().enumerate().take(n_groups) {
        let group_start = partition.group_starts[g];
        let group_end = partition.group_ends[g];
        debug_assert!(
            group_end <= chrom_length,
            "prefetch_window_ref_bytes: group {g} ends at {group_end} \
             past chrom length {chrom_length}",
        );
        let span = group_end - group_start + 1;
        bytes.reserve(span as usize);
        ref_fetcher.fetch_into(group_start, span, bytes)?;
    }
    Ok(())
}

/// Run layers 1+2+3 on group `g`, drive the EM via
/// [`run_em_columnar`], and emit a [`PosteriorRecord`] for the VCF
/// writer. Returns `Ok(None)` for groups the row-shape merger would
/// skip: post-unification `unified.alleles.len() < 2` (REF-only
/// after the max-alleles cap dropped every ALT into the OTHER pool)
/// and `unified.alleles.len() > cfg.max_alleles_lh_calc` (LH-cap
/// bind in pathological low-complexity regions). The skip mirrors
/// `PerGroupMerger::process_group`'s `Ok(None)` returns so the
/// per-category counts stay byte-identical against `main`.
#[allow(clippy::too_many_arguments)]
fn build_posterior_record_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    ref_bytes: &[u8],
    per_group_cfg: &PerGroupMergerConfig,
    posterior_cfg: &PosteriorEngineConfig,
    scratch: &mut ColumnarPipelineScratch,
    stats: &mut WindowRunStats,
) -> Result<Option<PosteriorRecord>, PosteriorEngineError> {
    let chrom_id = chunk.chrom_id;
    let group_start = partition.group_starts[group_idx];
    let group_end = partition.group_ends[group_idx];
    let n_samples = chunk.n_samples();

    // Layer 1 — allele unification reads REF bases pre-fetched by
    // the driver (no fetcher state on the hot path; see
    // [`prefetch_window_ref_bytes`]).
    debug_assert_eq!(
        ref_bytes.len() as u32,
        group_end - group_start + 1,
        "pre-fetched REF bytes length must match the group's reference span",
    );
    unify_alleles_columnar(
        chunk,
        partition,
        group_idx,
        ref_bytes,
        n_samples,
        per_group_cfg.max_alleles,
        &mut scratch.unify,
        &mut scratch.unified,
    )
    .map_err(|e| unify_alleles_error_to_merger(chrom_id, group_start, group_end, e))?;

    let n_alleles = scratch.unified.n_alleles();

    // Row-shape parity — see [`PerGroupMerger::process_group`].
    if n_alleles > per_group_cfg.max_alleles_lh_calc {
        stats.lh_cap_groups_skipped += 1;
        stats.lh_cap_alleles_in_skipped += n_alleles as u64;
        return Ok(None);
    }
    if n_alleles < 2 {
        // M23: surface the silent skip via a counter so operators
        // auditing why a known variant locus produced no record have
        // a breadcrumb. This path fires when post-unification every
        // ALT was absorbed into the OTHER pool by the max-alleles
        // cap, leaving only REF.
        stats.groups_skipped_post_unify_ref_only += 1;
        return Ok(None);
    }

    // Layer 2 — per-(sample, allele) scalar projection.
    project_scalars_columnar(
        chunk,
        partition,
        group_idx,
        &scratch.unified,
        &mut scratch.project,
        &mut scratch.projection,
    )
    .map_err(|e| project_scalars_error_to_merger(chrom_id, group_start, group_end, e))?;

    // Layer 3 — per-(sample, genotype) log-likelihood.
    let lh_ctx = LikelihoodContext {
        chrom_id,
        start: group_start,
        end: group_end,
        ploidy: per_group_cfg.ploidy,
    };
    compute_log_likelihoods_columnar(
        chunk,
        partition,
        group_idx,
        &scratch.unified,
        &scratch.projection,
        &lh_ctx,
        &[],
        &mut scratch.log_likelihoods,
    )
    .map_err(|e| compute_log_likelihoods_error_to_merger(chrom_id, group_start, group_end, e))?;

    // chain_anchor_flags table (row-shape boolean flat): `true` iff
    // the allele is a compound AND the sample's chain_anchor_count is
    // zero. Matches `build_chain_anchor_flags` in `per_group_merger`.
    // Built here for two consumers: `run_em_columnar`'s length-check
    // validation + the emitted `PosteriorRecord.chain_anchor_flags`
    // field that the VCF writer reads for the `CA` flag.
    scratch.chain_anchor_flags.clear();
    scratch
        .chain_anchor_flags
        .resize(n_samples * n_alleles, false);
    for allele_idx in 0..n_alleles {
        if !scratch.unified.is_compound[allele_idx] {
            continue;
        }
        for sample_idx in 0..n_samples {
            if scratch.unified.chain_anchor_count(allele_idx, sample_idx) == 0 {
                scratch.chain_anchor_flags[sample_idx * n_alleles + allele_idx] = true;
            }
        }
    }

    // Layer 4 — EM, mixture pre-pass, posterior summarisation,
    // exact-AF QUAL. Reads the layer 1/2/3 outputs as borrowed
    // slices through the `EmInputs` bundle.
    let n_genotypes = scratch.log_likelihoods.n_genotypes;
    let locus = crate::var_calling::posterior_engine::RecordLocus {
        chrom_id,
        start: group_start,
        end: group_end,
    };
    let em_outputs = {
        let view = ColumnarAllelesView::new(&scratch.unified);
        let inputs = EmInputs {
            locus,
            ploidy: per_group_cfg.ploidy,
            n_samples,
            n_genotypes,
            alleles: &view,
            scalars: &scratch.projection.scalars,
            log_likelihoods: &scratch.log_likelihoods.log_likelihoods,
            chain_anchor_flags_for_validation: &scratch.chain_anchor_flags,
        };
        run_em_columnar(
            inputs,
            posterior_cfg,
            &scratch.math,
            &mut scratch.record_scratch,
        )?
    };

    // Build the row-shape `Vec<MergedAllele>` for the emitted
    // `PosteriorRecord` at emit time, not at EM input time.
    let mut alleles: Vec<MergedAllele> = Vec::with_capacity(n_alleles);
    for a in 0..n_alleles {
        let lo = scratch.unified.constituent_offsets[a] as usize;
        let hi = scratch.unified.constituent_offsets[a + 1] as usize;
        let constituents: Vec<CompoundConstituent> = (lo..hi)
            .map(|k| CompoundConstituent {
                record_idx: scratch.unified.constituent_record_idx[k] as usize,
                local_allele_idx: scratch.unified.constituent_local_allele_idx[k] as usize,
            })
            .collect();
        alleles.push(MergedAllele {
            seq: scratch.unified.allele_seq(a).to_vec(),
            is_compound: scratch.unified.is_compound[a],
            constituents,
        });
    }

    let mut record = PosteriorRecord {
        locus,
        alleles,
        ploidy: per_group_cfg.ploidy,
        n_samples,
        n_genotypes: em_outputs.n_genotypes,
        allele_frequencies: em_outputs.allele_frequencies,
        compound_frequencies: em_outputs.compound_frequencies,
        posteriors: em_outputs.posteriors,
        best_genotype: em_outputs.best_genotype,
        gq_phred: em_outputs.gq_phred,
        qual_phred: em_outputs.qual_phred,
        scalars: scratch.projection.scalars.clone(),
        other_scalars: scratch.projection.other_scalars.clone(),
        chain_anchor_flags: scratch.chain_anchor_flags.clone(),
        diagnostics: em_outputs.diagnostics,
    };
    // Drop candidate ALT alleles no sample's argmax genotype uses, so
    // the VCF never carries AC=0 alleles (and no phantom 0/0 records
    // after `bcftools norm -m`). This is the column-native cohort path's
    // copy of the prune `run_em_for_record` applies on the row-shape /
    // from-bam path — `run_em_columnar` returns the raw `EmOutputs`
    // without it, so the chokepoint is here. No-op for biallelic groups.
    record.prune_unsupported_alleles(posterior_cfg.max_gq_phred);
    Ok(Some(record))
}

fn unify_alleles_error_to_merger(
    _chrom_id: u32,
    _group_start: u32,
    _group_end: u32,
    e: UnifyAllelesError,
) -> PosteriorEngineError {
    let inner = match e {
        UnifyAllelesError::MissingCompoundAlleleBytes {
            chrom_id,
            group_start,
            group_end,
            record_idx,
            local_allele_idx,
        } => PerGroupMergerError::MissingCompoundAlleleBytes {
            chrom_id,
            start: group_start,
            end: group_end,
            record_idx: record_idx as usize,
            local_allele_idx: local_allele_idx as usize,
        },
    };
    inner.into()
}

fn project_scalars_error_to_merger(
    _chrom_id: u32,
    _group_start: u32,
    _group_end: u32,
    e: ProjectScalarsError,
) -> PosteriorEngineError {
    use crate::var_calling::cohort_block::kernels::project_scalars::CompoundPhase as CnPhase;
    use crate::var_calling::per_group_merger::CompoundPhase as RowPhase;
    let inner = match e {
        ProjectScalarsError::ZeroObservationConstituent {
            chrom_id,
            group_start,
            group_end,
            phase,
            sample_idx,
            allele_idx,
            record_idx,
            local_allele_idx,
        } => PerGroupMergerError::ZeroObservationConstituent {
            chrom_id,
            start: group_start,
            end: group_end,
            phase: match phase {
                CnPhase::QualityGather => RowPhase::QualityGather,
                CnPhase::ConstituentSubtraction => RowPhase::ConstituentSubtraction,
            },
            sample_idx,
            allele_idx,
            record_idx,
            local_allele_idx,
        },
        ProjectScalarsError::NoQualityForChainAnchoredCompound {
            chrom_id,
            group_start,
            group_end,
            sample_idx,
            allele_idx,
            inter,
        } => PerGroupMergerError::NoQualityForChainAnchoredCompound {
            chrom_id,
            start: group_start,
            end: group_end,
            sample_idx,
            allele_idx,
            inter,
        },
    };
    inner.into()
}

fn compute_log_likelihoods_error_to_merger(
    chrom_id: u32,
    group_start: u32,
    group_end: u32,
    e: ComputeLogLikelihoodsError,
) -> PosteriorEngineError {
    use crate::var_calling::posterior_engine::RecordLocus;
    match e {
        ComputeLogLikelihoodsError::DegenerateLikelihood {
            chrom_id,
            start,
            end,
            sample_idx,
            genotype_idx,
            kind,
        } => PerGroupMergerError::DegenerateLikelihood {
            chrom_id,
            start,
            end,
            sample_idx,
            genotype_idx,
            kind,
        }
        .into(),
        // B4: surface the upstream-invariant break directly with the
        // real `n_alleles` and group locus, instead of synthesising a
        // `DegenerateLikelihood` with `usize::MAX` placeholders that
        // misled operators into chasing a non-existent sample /
        // genotype.
        ComputeLogLikelihoodsError::NAllelesExceedsBitmask { n_alleles } => {
            PosteriorEngineError::NAllelesExceedsBitmask {
                locus: RecordLocus {
                    chrom_id,
                    start: group_start,
                    end: group_end,
                },
                n_alleles,
            }
        }
    }
}

/// Construct the `SharedRefFetcher` shape `run_window` expects from
/// an owned [`ChromRefFetcher`] implementor. The Arc is still the
/// natural type for the per-chromosome ref fetcher held by the
/// driver — `run_window` derefs into a `&dyn ChromRefFetcher` for
/// the columnar layer 1 fetch.
// M28: verb-named constructor. Bare-noun function names are reserved
// for field-like accessors; converting constructors are
// `from_X` / `into_X` / `to_X`.
pub fn into_shared_ref_fetcher<F>(fetcher: F) -> SharedRefFetcher
where
    F: ChromRefFetcher + Send + 'static,
{
    Arc::new(fetcher)
}

/// Errors from [`unified_alleles_for_group_columnar`]. Wraps the
/// kernel error and the ref-fetcher error so callers see a single
/// surface.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum WorkerUnifyError {
    /// Kernel-side failure (compound projection couldn't find an
    /// anchoring sample, etc.).
    #[error(transparent)]
    Kernel(#[from] UnifyAllelesError),

    /// `fetch_into` failed while pulling the group's REF span — the
    /// kernel needs those bytes as input.
    #[error("reference fetch: {0}")]
    RefFetch(#[from] ChromRefFetchError),
}

/// Phase A.1 sub-step 1.4 — column-native allele unification for one
/// group of a [`WindowPartition`], fetched through the worker's
/// ref-fetcher convention. Retained from layer 1 for direct external
/// access to the unifier without the full pipeline; `run_window`
/// drives the same kernel via [`ColumnarPipelineScratch`].
#[allow(clippy::too_many_arguments)]
pub fn unified_alleles_for_group_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    ref_fetcher: &dyn ChromRefFetcher,
    max_alleles: usize,
    scratch: &mut UnifyAllelesScratch,
    ref_buf: &mut Vec<u8>,
    out: &mut UnifiedAllelesColumns,
) -> Result<(), WorkerUnifyError> {
    let group_start = partition.group_starts[group_idx];
    let group_end = partition.group_ends[group_idx];
    let span = group_end - group_start + 1;
    ref_fetcher.fetch_into(group_start, span, ref_buf)?;
    unify_alleles_columnar(
        chunk,
        partition,
        group_idx,
        ref_buf,
        chunk.n_samples(),
        max_alleles,
        scratch,
        out,
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::cohort_block::partition::{
        PartitionScratch, WindowPartition, partition_window,
    };
    use crate::var_calling::cohort_block::test_helpers::{
        build_overlapping_variant_group, loaded_chunk, record, ref_plus_alt,
    };
    use crate::var_calling::per_group_merger::{PerGroupMerger, SharedRefFetcher};
    use crate::var_calling::posterior_engine::PosteriorEngine;
    use crate::var_calling::variant_grouping::{GrouperError, OverlappingVariantGroup};

    /// Nits (Wave 6): byte-identity oracle shape — `(seq,
    /// is_compound, constituents)`. Naming the triple drops the
    /// `clippy::type_complexity` lint trips at the helper locals.
    type OracleAllele = (Vec<u8>, bool, Vec<(usize, usize)>);

    /// Build a chunk + partition for the worker's adapter tests.
    /// Returns `(chunk, partition)` so callers can probe both
    /// surfaces.
    fn fixture(
        per_sample_records: Vec<Vec<crate::pileup_record::PileupRecord>>,
        window_start: u32,
        window_end: u32,
        max_group_span: u32,
    ) -> (MaterialisedChunk, WindowPartition) {
        let n = per_sample_records.len();
        let (chunk, _) = loaded_chunk(0, 1..1000, per_sample_records);
        let mut scratch = PartitionScratch::with_n_samples(n);
        let mut partition = WindowPartition::empty();
        partition_window(
            &chunk,
            &(window_start..window_end),
            &[],
            max_group_span,
            &mut scratch,
            &mut partition,
        )
        .unwrap();
        (chunk, partition)
    }

    #[test]
    fn adapter_carries_chrom_id_start_end_from_partition() {
        // MNP at 100 (ref_span=10 → reach 109) + SNP at 105.
        // 105 ≤ 109 → joins → one group spanning [100, 109].
        let s0 = vec![
            record(
                100,
                vec![
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"AAAAAAAAAA",
                        3,
                        -1.0,
                        &[],
                    ),
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"AAACAAAAAA",
                        4,
                        -1.0,
                        &[],
                    ),
                ],
            ),
            record(105, ref_plus_alt(3, 4)),
        ];
        let (mut chunk, partition) = fixture(vec![s0], 1, 200, 100);
        chunk.chrom_id = 9;
        let group = build_overlapping_variant_group(&chunk, &partition, 0, chunk.n_samples(), 9);
        assert_eq!(group.chrom_id, 9);
        assert_eq!(group.start, 100);
        assert_eq!(group.end, 109);
        assert_eq!(group.records.len(), 2);
        assert_eq!(group.records[0].pos, 100);
        assert_eq!(group.records[1].pos, 105);
    }

    #[test]
    fn adapter_emits_one_per_sample_slot_per_position_with_some_for_carriers() {
        // Two samples, both at position 50; only sample 1 at
        // position 60.
        let s0 = vec![record(50, ref_plus_alt(2, 3))];
        let s1 = vec![
            record(50, ref_plus_alt(2, 3)),
            record(60, ref_plus_alt(2, 3)),
        ];
        let (chunk, partition) = fixture(vec![s0, s1], 1, 100, 5);
        let n_samples = chunk.n_samples();

        let group_0 = build_overlapping_variant_group(&chunk, &partition, 0, n_samples, 0);
        assert_eq!(group_0.records.len(), 1);
        assert_eq!(group_0.records[0].per_sample.len(), 2);
        assert!(group_0.records[0].per_sample[0].is_some());
        assert!(group_0.records[0].per_sample[1].is_some());

        let group_1 = build_overlapping_variant_group(&chunk, &partition, 1, n_samples, 0);
        assert_eq!(group_1.records.len(), 1);
        assert_eq!(group_1.records[0].per_sample.len(), 2);
        assert!(group_1.records[0].per_sample[0].is_none());
        assert!(group_1.records[0].per_sample[1].is_some());
    }

    #[test]
    fn adapter_round_trips_payload_via_materialised_record() {
        // Sample 0 has a SNP at 50 (REF=A obs=3, ALT=T obs=4). The
        // adapter should materialise the same PileupRecord that
        // `chunk.materialise_record(0, 0)` returns.
        let s0 = vec![record(50, ref_plus_alt(3, 4))];
        let (chunk, partition) = fixture(vec![s0], 1, 100, 5);
        let group = build_overlapping_variant_group(&chunk, &partition, 0, chunk.n_samples(), 0);

        let expected = chunk.materialise_record(0, 0);
        let actual = group.records[0].per_sample[0].as_ref().unwrap();
        assert_eq!(actual, &expected);
    }

    #[test]
    fn run_window_with_empty_partition_clears_output_and_returns_ok() {
        // A real ChromRefFetcher isn't needed for the empty-partition
        // fast path. Build a dangling Arc that would panic on use,
        // and confirm the function never reaches it.
        let s0 = vec![record(50, ref_plus_alt(2, 3))];
        let (chunk, _) = loaded_chunk(0, 1..100, vec![s0]);
        let partition = WindowPartition::empty();
        let mut output: Vec<PosteriorRecord> = Vec::new();
        let mut scratch = ColumnarPipelineScratch::empty();

        // Empty partition ⇒ no groups ⇒ no REF bytes needed.
        let mut stats = WindowRunStats::default();
        let result = run_window(
            &chunk,
            &partition,
            &[],
            PerGroupMergerConfig::default(),
            &PosteriorEngineConfig::default(),
            &mut scratch,
            &mut output,
            &mut stats,
        );
        assert!(result.is_ok());
        assert!(output.is_empty());
    }

    /// B4: `compute_log_likelihoods_error_to_merger` maps the kernel's
    /// `NAllelesExceedsBitmask` directly to
    /// `PosteriorEngineError::NAllelesExceedsBitmask`, preserving the
    /// real `n_alleles` and the group locus. Replaces the prior
    /// behaviour of silently rewriting the upstream-invariant break
    /// as a `DegenerateLikelihood { sample_idx: usize::MAX,
    /// genotype_idx: usize::MAX, kind: NaN }` placeholder.
    #[test]
    fn compute_log_likelihoods_error_to_merger_preserves_n_alleles_and_locus() {
        let err = compute_log_likelihoods_error_to_merger(
            7,
            123,
            456,
            ComputeLogLikelihoodsError::NAllelesExceedsBitmask { n_alleles: 99 },
        );
        match err {
            PosteriorEngineError::NAllelesExceedsBitmask { locus, n_alleles } => {
                assert_eq!(locus.chrom_id, 7);
                assert_eq!(locus.start, 123);
                assert_eq!(locus.end, 456);
                assert_eq!(n_alleles, 99);
            }
            other => panic!("expected NAllelesExceedsBitmask, got {other:?}",),
        }
    }

    /// In-memory `ChromRefFetcher` for the column-native unify
    /// integration test. Mirrors the `MockRef` inside the kernel
    /// tests.
    #[derive(Clone)]
    struct MockRef {
        seq: Vec<u8>,
        /// 1-based position of `seq[0]`.
        base_offset: u32,
    }

    impl crate::fasta::fetcher::sealed::Sealed for MockRef {}

    impl crate::fasta::fetcher::ChromRefFetcher for MockRef {
        fn length(&self) -> u32 {
            self.base_offset.saturating_sub(1) + self.seq.len() as u32
        }
        fn fetch(
            &self,
            start_1based: u32,
            length: u32,
        ) -> Result<Vec<u8>, crate::fasta::fetcher::ChromRefFetchError> {
            let start_idx = (start_1based - self.base_offset) as usize;
            Ok(self.seq[start_idx..start_idx + length as usize].to_vec())
        }
        fn iter_bases<'a>(
            &'a self,
        ) -> Result<
            Box<dyn Iterator<Item = Result<u8, crate::fasta::fetcher::ChromRefFetchError>> + 'a>,
            crate::fasta::fetcher::ChromRefFetchError,
        > {
            Ok(Box::new(self.seq.iter().copied().map(Ok)))
        }
    }

    /// M22: when multiple prunable ALTs have equal `cohort_count`,
    /// the column-native `enforce_max_alleles_columnar` tie-break
    /// must match the row-shape `PerGroupMerger` oracle's tie-break
    /// — otherwise byte-identity against `main` could silently break
    /// on a real-data input that hits the tie. Fixture: 3 ALTs all
    /// at cohort_count=5; `max_alleles = 3` keeps REF + 2 ALTs and
    /// forces the cap to drop exactly one of the three.
    #[test]
    fn unify_max_alleles_ties_match_row_shape_kernel() {
        // Two samples; pos 100 covered by both. Each sample
        // contributes the same three ALTs (T, C, G), each at
        // num_obs=5, so the cohort-pooled count for each is 10. The
        // cap (max_alleles=3 → REF + 2 ALTs) must drop one of the
        // three; the column-native and row-shape kernels must agree
        // on which.
        let alts = || {
            vec![
                crate::var_calling::cohort_block::test_helpers::allele(b"A", 5, -1.0, &[]),
                crate::var_calling::cohort_block::test_helpers::allele(b"T", 5, -1.0, &[]),
                crate::var_calling::cohort_block::test_helpers::allele(b"C", 5, -1.0, &[]),
                crate::var_calling::cohort_block::test_helpers::allele(b"G", 5, -1.0, &[]),
            ]
        };
        let s0 = vec![record(100, alts())];
        let s1 = vec![record(100, alts())];
        let (chunk, partition) = fixture(vec![s0, s1], 1, 200, 5);

        let ref_full = vec![b'A'; 200];
        let fetcher = MockRef {
            seq: ref_full.clone(),
            base_offset: 1,
        };

        // Column-native pipeline.
        let mut cn_scratch = ColumnarPipelineScratch::empty();
        let mut cn_out: Vec<PosteriorRecord> = Vec::new();
        let cfg = PerGroupMergerConfig::new(2, 3, 64, 32).expect("merger config");
        let mut cn_pre_fetched: Vec<Vec<u8>> = Vec::new();
        prefetch_window_ref_bytes(&partition, &fetcher, &mut cn_pre_fetched)
            .expect("prefetch succeeded");
        let mut cn_stats = WindowRunStats::default();
        run_window(
            &chunk,
            &partition,
            &cn_pre_fetched,
            cfg,
            &PosteriorEngineConfig::default(),
            &mut cn_scratch,
            &mut cn_out,
            &mut cn_stats,
        )
        .expect("columnar run_window succeeded");

        // Row-shape oracle.
        let row_fetcher: SharedRefFetcher = Arc::new(fetcher);
        let mut row_out: Vec<PosteriorRecord> = Vec::new();
        let groups: Vec<OverlappingVariantGroup> = (0..partition.n_groups())
            .map(|g| {
                build_overlapping_variant_group(
                    &chunk,
                    &partition,
                    g,
                    chunk.n_samples(),
                    chunk.chrom_id,
                )
            })
            .collect();
        let merger = PerGroupMerger::with_config(
            groups.into_iter().map(Ok::<_, GrouperError>),
            row_fetcher,
            cfg,
        );
        let engine = PosteriorEngine::with_config(merger, PosteriorEngineConfig::default());
        for record in engine {
            row_out.push(record.expect("row-shape pipeline succeeded"));
        }

        // Both pipelines must emit the same set of ALT-allele bytes.
        assert_eq!(cn_out.len(), row_out.len(), "record count mismatch");
        for (cn, row) in cn_out.iter().zip(row_out.iter()) {
            assert!(
                cn.alleles.len() <= 3,
                "max_alleles cap should hold ({} alleles emitted)",
                cn.alleles.len(),
            );
            // Pin the *same* allele set is picked under tie — the
            // strongest possible assertion against tie-break drift.
            assert_eq!(cn, row, "tie-break must match row-shape oracle");
        }
    }

    #[test]
    fn column_native_unify_matches_row_shape_on_worker_fixture() {
        // Same multi-position fixture as the existing adapter test:
        // sample 0 has a 10-bp MNP at pos 100 reaching to 109, plus
        // an SNP at pos 105. The kernel + the row-shape PerGroupMerger
        // should produce identical (seq, is_compound, constituents)
        // tuples.
        let s0 = vec![
            record(
                100,
                vec![
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"AAAAAAAAAA",
                        3,
                        -1.0,
                        &[],
                    ),
                    crate::var_calling::cohort_block::test_helpers::allele(
                        b"AAACAAAAAA",
                        4,
                        -1.0,
                        &[],
                    ),
                ],
            ),
            record(105, ref_plus_alt(3, 4)),
        ];
        let (chunk, partition) = fixture(vec![s0], 1, 200, 100);
        assert_eq!(partition.n_groups(), 1);

        let ref_full = vec![b'A'; 200];
        let fetcher = MockRef {
            seq: ref_full.clone(),
            base_offset: 1,
        };

        // Column-native via the new worker-side entry point.
        let mut scratch = UnifyAllelesScratch::new();
        let mut ref_buf = Vec::new();
        let mut cn_out = UnifiedAllelesColumns::empty();
        unified_alleles_for_group_columnar(
            &chunk,
            &partition,
            0,
            &fetcher,
            16,
            &mut scratch,
            &mut ref_buf,
            &mut cn_out,
        )
        .expect("column-native unify succeeded");

        let cn_alleles: Vec<OracleAllele> = (0..cn_out.n_alleles())
            .map(|i| {
                let lo = cn_out.constituent_offsets[i] as usize;
                let hi = cn_out.constituent_offsets[i + 1] as usize;
                let constituents: Vec<(usize, usize)> = (lo..hi)
                    .map(|k| {
                        (
                            cn_out.constituent_record_idx[k] as usize,
                            cn_out.constituent_local_allele_idx[k] as usize,
                        )
                    })
                    .collect();
                (
                    cn_out.allele_seq(i).to_vec(),
                    cn_out.is_compound[i],
                    constituents,
                )
            })
            .collect();

        // Row-shape oracle via the existing kernel.
        let fetcher_arc: SharedRefFetcher = Arc::new(fetcher);
        let group = build_overlapping_variant_group(&chunk, &partition, 0, chunk.n_samples(), 0);
        let config = PerGroupMergerConfig::new(2, 16, 64, 32).expect("merger config");
        let iter: Vec<Result<OverlappingVariantGroup, GrouperError>> = vec![Ok(group)];
        let mut merger = PerGroupMerger::with_config(iter.into_iter(), fetcher_arc, config);
        let record = merger
            .next()
            .expect("merger yields at least one item")
            .expect("merger succeeded");
        let oracle: Vec<OracleAllele> = record
            .alleles
            .iter()
            .map(|a| {
                (
                    a.seq.clone(),
                    a.is_compound,
                    a.constituents
                        .iter()
                        .map(|c| (c.record_idx, c.local_allele_idx))
                        .collect(),
                )
            })
            .collect();

        assert_eq!(cn_alleles, oracle);
    }

    /// End-to-end byte-identity: same partition + chunk, fed through
    /// the column-native `run_window` and through a row-shape
    /// `PerGroupMerger` → `PosteriorEngine` pipeline. The emitted
    /// `PosteriorRecord` sequences must compare equal.
    #[test]
    fn run_window_columnar_matches_row_shape_pipeline() {
        // Two samples with a SNP at 50 + a multi-position group at
        // 100..109 with a compound candidate, exercising layers
        // 1 (compounds), 2 (compound projection), 3 (chain-broken
        // path on the sample without the compound), 4 (EM).
        let s0_snp = record(50, ref_plus_alt(3, 4));
        let s0_mnp = record(
            100,
            vec![
                crate::var_calling::cohort_block::test_helpers::allele(b"AAAAAAAAAA", 3, -1.0, &[]),
                crate::var_calling::cohort_block::test_helpers::allele(b"AAACAAAAAA", 4, -1.0, &[]),
            ],
        );
        let s0_snp_inside = record(105, ref_plus_alt(3, 4));
        let s1_snp = record(50, ref_plus_alt(5, 6));
        let s1_snp_inside = record(105, ref_plus_alt(4, 5));
        let s0 = vec![s0_snp, s0_mnp, s0_snp_inside];
        let s1 = vec![s1_snp, s1_snp_inside];
        let (chunk, partition) = fixture(vec![s0, s1], 1, 200, 100);
        let ref_full = vec![b'A'; 200];
        let fetcher = MockRef {
            seq: ref_full.clone(),
            base_offset: 1,
        };

        // Column-native path.
        let mut cn_scratch = ColumnarPipelineScratch::empty();
        let mut cn_out: Vec<PosteriorRecord> = Vec::new();
        let cfg = PerGroupMergerConfig::new(2, 16, 64, 32).expect("merger config");
        let mut cn_pre_fetched: Vec<Vec<u8>> = Vec::new();
        prefetch_window_ref_bytes(&partition, &fetcher, &mut cn_pre_fetched)
            .expect("prefetch succeeded");
        let mut cn_stats = WindowRunStats::default();
        run_window(
            &chunk,
            &partition,
            &cn_pre_fetched,
            cfg,
            &PosteriorEngineConfig::default(),
            &mut cn_scratch,
            &mut cn_out,
            &mut cn_stats,
        )
        .expect("columnar run_window succeeded");

        // Row-shape oracle path (Phase A.0).
        let row_fetcher: SharedRefFetcher = Arc::new(fetcher);
        let mut row_out: Vec<PosteriorRecord> = Vec::new();
        let groups: Vec<OverlappingVariantGroup> = (0..partition.n_groups())
            .map(|g| {
                build_overlapping_variant_group(
                    &chunk,
                    &partition,
                    g,
                    chunk.n_samples(),
                    chunk.chrom_id,
                )
            })
            .collect();
        let merger = PerGroupMerger::with_config(
            groups.into_iter().map(Ok::<_, GrouperError>),
            row_fetcher,
            cfg,
        );
        let engine = PosteriorEngine::with_config(merger, PosteriorEngineConfig::default());
        for record in engine {
            row_out.push(record.expect("row-shape pipeline succeeded"));
        }

        assert_eq!(cn_out.len(), row_out.len());
        for (cn, row) in cn_out.iter().zip(row_out.iter()) {
            assert_eq!(cn, row);
        }
    }

    /// Regression test for the cap_protected over-protection bug
    /// that bit on the 3-tomato real-data run before commit (this
    /// file). The pre-fix layer 1.2 admit-compound code marked an
    /// existing per-position match as `cap_protected = true`,
    /// shrinking the prunable pool below what
    /// [`enforce_max_alleles_columnar`] needs — the cap couldn't
    /// trim alleles down to `max_alleles`, yielding emitted records
    /// with more than `max_alleles` ALTs on the columnar path while
    /// the row-shape path emitted only `max_alleles - 1` ALTs.
    /// Byte-identity against the row-shape `PerGroupMerger` oracle
    /// on a fixture that triggers the cap with both compound + bare
    /// per-position alleles present.
    #[test]
    fn run_window_columnar_matches_row_shape_under_cap_with_compounds() {
        // Group span [100, 101]: two cohort positions joined into one
        // group. Many ALTs at pos 100 to force the cap; chain-linked
        // alleles across pos 100 + pos 101 form compound candidates.
        let s0 = vec![
            record(
                100,
                vec![
                    crate::var_calling::cohort_block::test_helpers::allele(b"A", 10, -1.0, &[]),
                    crate::var_calling::cohort_block::test_helpers::allele(b"T", 6, -1.0, &[101]),
                    crate::var_calling::cohort_block::test_helpers::allele(b"G", 2, -1.0, &[]),
                    crate::var_calling::cohort_block::test_helpers::allele(b"C", 1, -1.0, &[]),
                ],
            ),
            record(
                101,
                vec![
                    crate::var_calling::cohort_block::test_helpers::allele(b"A", 10, -1.0, &[]),
                    crate::var_calling::cohort_block::test_helpers::allele(b"G", 5, -1.0, &[101]),
                ],
            ),
        ];
        let s1 = vec![
            record(
                100,
                vec![
                    crate::var_calling::cohort_block::test_helpers::allele(b"A", 8, -1.0, &[]),
                    crate::var_calling::cohort_block::test_helpers::allele(b"G", 3, -1.0, &[]),
                    crate::var_calling::cohort_block::test_helpers::allele(b"T", 1, -1.0, &[]),
                ],
            ),
            record(
                101,
                vec![crate::var_calling::cohort_block::test_helpers::allele(
                    b"A",
                    10,
                    -1.0,
                    &[],
                )],
            ),
        ];
        let (chunk, partition) = fixture(vec![s0, s1], 1, 200, 5);

        let ref_full = vec![b'A'; 200];
        let fetcher = MockRef {
            seq: ref_full.clone(),
            base_offset: 1,
        };

        // `max_alleles = 3` forces the cap; compounds are present.
        let cfg = PerGroupMergerConfig::new(2, 3, 64, 32).expect("merger config");

        let mut cn_scratch = ColumnarPipelineScratch::empty();
        let mut cn_out: Vec<PosteriorRecord> = Vec::new();
        let mut cn_pre_fetched: Vec<Vec<u8>> = Vec::new();
        prefetch_window_ref_bytes(&partition, &fetcher, &mut cn_pre_fetched)
            .expect("prefetch succeeded");
        let mut cn_stats = WindowRunStats::default();
        run_window(
            &chunk,
            &partition,
            &cn_pre_fetched,
            cfg,
            &PosteriorEngineConfig::default(),
            &mut cn_scratch,
            &mut cn_out,
            &mut cn_stats,
        )
        .expect("columnar run_window succeeded");

        let row_fetcher: SharedRefFetcher = Arc::new(fetcher);
        let mut row_out: Vec<PosteriorRecord> = Vec::new();
        let groups: Vec<OverlappingVariantGroup> = (0..partition.n_groups())
            .map(|g| {
                build_overlapping_variant_group(
                    &chunk,
                    &partition,
                    g,
                    chunk.n_samples(),
                    chunk.chrom_id,
                )
            })
            .collect();
        let merger = PerGroupMerger::with_config(
            groups.into_iter().map(Ok::<_, GrouperError>),
            row_fetcher,
            cfg,
        );
        let engine = PosteriorEngine::with_config(merger, PosteriorEngineConfig::default());
        for record in engine {
            row_out.push(record.expect("row-shape pipeline succeeded"));
        }

        assert_eq!(cn_out.len(), row_out.len());
        for (cn, row) in cn_out.iter().zip(row_out.iter()) {
            assert!(
                cn.alleles.len() <= 3,
                "max_alleles cap should hold ({} alleles emitted)",
                cn.alleles.len(),
            );
            assert_eq!(cn, row);
        }
    }
}
