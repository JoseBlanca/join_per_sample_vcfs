//! `#[cfg(test)]`-only fixture builders shared across the
//! [`cohort_block`](super) submodules' test suites.
//!
//! These are crate-visible (`pub(crate)`) so siblings can pull
//! them in without each having its own copy. The shapes match the
//! invariants of [`PileupRecord`] and [`AlleleObservation`] —
//! `alleles[0]` is REF, `chain_ids` is sorted-dedup'd, and so on.

use std::ops::Range;

use crate::pileup_record::{AlleleObservation, AlleleSupportStats, ChainId, PileupRecord};
use crate::var_calling::cohort_block::columns::{MaterialisedChunk, SampleColumns};
use crate::var_calling::cohort_block::loader::{
    ChunkLoadExtent, ChunkLoadScratch, load_chunk_from_iters,
};
use crate::var_calling::cohort_block::partition::WindowPartition;
use crate::var_calling::cohort_block::pre_pass::{
    FixBoundariesError, FixBoundariesScratch, fix_boundaries,
};
use crate::var_calling::per_position_merger::PerPositionPileups;
use crate::var_calling::variant_grouping::OverlappingVariantGroup;

/// Build an [`AlleleObservation`] with a fully-specified
/// `(num_obs, q_sum, chain_ids)` triple. The strand-bias and MAPQ
/// moments are filled in from `num_obs` with stable default shapes
/// — no downstream test on this slice inspects them.
pub(crate) fn allele(
    seq: &[u8],
    num_obs: u32,
    q_sum: f64,
    chain_ids: &[ChainId],
) -> AlleleObservation {
    AlleleObservation::new(
        seq.to_vec(),
        AlleleSupportStats::new(
            num_obs,
            q_sum,
            num_obs,
            0,
            0,
            num_obs * 30,
            u64::from(num_obs) * 900,
        ),
        chain_ids.to_vec(),
    )
}

/// Build a [`PileupRecord`] at `pos` on chromosome 0.
pub(crate) fn record(pos: u32, alleles: Vec<AlleleObservation>) -> PileupRecord {
    PileupRecord::new(0, pos, alleles)
}

/// Build an allele with an arbitrary `num_obs` and a fixed
/// `q_sum = -1.0`. Useful when the variant-filter is the only
/// thing being tested.
pub(crate) fn observed(seq: &[u8], num_obs: u32) -> AlleleObservation {
    allele(seq, num_obs, -1.0, &[])
}

/// REF entry with `num_obs > 0` (i.e. some sample carries REF as
/// its haplotype at this position).
pub(crate) fn ref_obs(num_obs: u32) -> AlleleObservation {
    observed(b"A", num_obs)
}

/// REF + observed ALT — defines a variant position.
pub(crate) fn ref_plus_alt(ref_obs_count: u32, alt_obs_count: u32) -> Vec<AlleleObservation> {
    vec![ref_obs(ref_obs_count), observed(b"T", alt_obs_count)]
}

/// REF + ALT bucket present but with zero observations — the
/// "bucket exists but no read currently carries it" case from
/// the [`AlleleObservation`] doc.
pub(crate) fn ref_plus_unobserved_alt(ref_obs_count: u32) -> Vec<AlleleObservation> {
    vec![ref_obs(ref_obs_count), observed(b"T", 0)]
}

/// Build a chunk by running the loader against synthetic per-sample
/// record vectors. The fixture mirrors what an end-to-end call
/// against PSP readers would produce.
pub(crate) fn run_loader(
    chrom_id: u32,
    range: Range<u32>,
    per_sample_records: Vec<Vec<PileupRecord>>,
    mut carryover: Vec<SampleColumns>,
) -> MaterialisedChunk {
    let n_samples = per_sample_records.len();
    let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
    let mut out = MaterialisedChunk::with_n_samples(n_samples);
    let iters: Vec<_> = per_sample_records
        .into_iter()
        .map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>))
        .collect();
    let span = range.end - range.start;
    // Legacy behaviour: no variant target, no extension (max_span ==
    // initial_span ⇒ single pull attempt).
    load_chunk_from_iters(
        &mut scratch,
        &mut out,
        ChunkLoadExtent {
            chrom_id,
            range_start: range.start,
            initial_span: span,
            target_variants: 0,
            max_span: span,
        },
        iters,
        &mut carryover,
    )
    .expect("load_chunk_from_iters succeeded");
    out
}

/// Build a fresh chunk via [`run_loader`] and pair it with an
/// empty carryover sized to the cohort. Convenient setup for the
/// pre-pass tests.
pub(crate) fn loaded_chunk(
    chrom_id: u32,
    range: Range<u32>,
    per_sample_records: Vec<Vec<PileupRecord>>,
) -> (MaterialisedChunk, Vec<SampleColumns>) {
    let n = per_sample_records.len();
    let chunk = run_loader(
        chrom_id,
        range,
        per_sample_records,
        vec![SampleColumns::empty(); n],
    );
    let carryover = (0..n).map(|_| SampleColumns::empty()).collect();
    (chunk, carryover)
}

/// Run the pre-pass with `target_window_count = 1` and a freshly
/// allocated [`FixBoundariesScratch`].
pub(crate) fn run_pre_pass(
    chunk: &mut MaterialisedChunk,
    carryover: &mut [SampleColumns],
    max_group_span: u32,
) -> Result<(), FixBoundariesError> {
    let mut scratch = FixBoundariesScratch::new();
    fix_boundaries(chunk, carryover, &mut scratch, max_group_span, 1)
}

/// Mi16: byte-identity oracle that materialises one
/// [`OverlappingVariantGroup`] from group `g` of `partition`. Per-
/// position records are materialised on demand from the chunk's
/// columnar slices via
/// [`SampleColumns::materialise_record`](super::columns::SampleColumns::materialise_record).
///
/// Retained from Phase A.0 as the byte-identity oracle in the
/// kernels' unit tests. Not on the production path — `run_window`
/// builds merged records column-natively via
/// `build_posterior_record_columnar`.
///
/// Lifted from `worker.rs` so the cross-module test reach reads as
/// "test helper imported from `test_helpers`" rather than "test
/// helper reached back through the `worker` module surface".
pub(crate) fn build_overlapping_variant_group(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    g: usize,
    n_samples: usize,
    chrom_id: u32,
) -> OverlappingVariantGroup {
    let position_range = partition.position_range_for_group(g);
    let mut records: Vec<PerPositionPileups> = Vec::with_capacity(position_range.len());
    for p in position_range {
        let pos = partition.positions[p];
        let mut per_sample: Vec<Option<PileupRecord>> = vec![None; n_samples];
        let sample_range = partition.sample_range_for_position(p);
        for k in sample_range {
            let sample_idx = partition.samples_at_pos[k] as usize;
            let row_idx = partition.rows_at_pos[k] as usize;
            per_sample[sample_idx] =
                Some(chunk.per_sample[sample_idx].materialise_record(chrom_id, row_idx));
        }
        records.push(PerPositionPileups {
            chrom_id,
            pos,
            per_sample,
        });
    }
    OverlappingVariantGroup {
        chrom_id,
        start: partition.group_starts[g],
        end: partition.group_ends[g],
        records,
    }
}
