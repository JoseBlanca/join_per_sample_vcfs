//! Per-window variant-calling math.
//!
//! Phase A.0 implementation: builds an
//! [`OverlappingVariantGroup`] from each group in a
//! [`WindowPartition`] + the chunk's columnar slices, then pipes
//! the groups through the existing
//! [`PerGroupMerger`](crate::var_calling::per_group_merger::PerGroupMerger)
//! and
//! [`PosteriorEngine`](crate::var_calling::posterior_engine::PosteriorEngine)
//! kernels to produce
//! [`PosteriorRecord`](crate::var_calling::posterior_engine::PosteriorRecord)s.
//!
//! Reusing the kernels lets the cohort-level rewrite reach
//! byte-identical VCFs against `main` without re-deriving the
//! allele-unification, likelihood, EM, and QUAL math from scratch.
//! Phase A.1 (deferred) replaces this adapter with native columnar
//! kernels operating directly on
//! [`SampleColumns`](super::columns::SampleColumns) and
//! [`WindowPartition`].

use std::sync::Arc;

use crate::fasta::fetcher::ChromRefFetcher;
use crate::pileup_record::PileupRecord;
use crate::var_calling::cohort_block::columns::MaterialisedChunk;
use crate::var_calling::cohort_block::partition::WindowPartition;
use crate::var_calling::per_group_merger::{
    PerGroupMerger, PerGroupMergerConfig, SharedRefFetcher,
};
use crate::var_calling::per_position_merger::PerPositionPileups;
use crate::var_calling::posterior_engine::{
    PosteriorEngine, PosteriorEngineConfig, PosteriorEngineError, PosteriorRecord,
};
use crate::var_calling::variant_grouping::{GrouperError, OverlappingVariantGroup};

/// Run the worker math on one window's partition.
///
/// **Inputs.**
/// - `chunk`: the loaded `MaterialisedChunk` whose columns the
///   partition's `(sample_idx, row_idx)` handles refer into.
/// - `partition`: the variant-group partition produced by
///   [`partition_window`](super::partition::partition_window) for
///   one of `chunk.windows`. Its `n_groups()` may be 0 — empty
///   partitions are a fast no-op.
/// - `ref_fetcher`: shared reference fetcher used by the
///   per-group merger to gather REF bases for the group span. The
///   driver typically holds an `Arc<StreamingChromRefFetcher>` per
///   chromosome.
/// - `per_group_config`, `posterior_config`: tuning knobs for the
///   downstream kernels.
/// - `output`: appended-to record buffer. Cleared on entry so the
///   driver can reuse it across windows.
///
/// **Algorithm.** For each variant group `g` in `partition`,
/// materialise its `(start, end, records: Vec<PerPositionPileups>)`
/// from the columnar chunk + partition handles. Feed the resulting
/// `Iterator<Item = Result<OverlappingVariantGroup, GrouperError>>`
/// to `PerGroupMerger`, then to `PosteriorEngine`, then drain into
/// `output`.
pub fn run_window(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    ref_fetcher: SharedRefFetcher,
    per_group_config: PerGroupMergerConfig,
    posterior_config: PosteriorEngineConfig,
    output: &mut Vec<PosteriorRecord>,
) -> Result<(), PosteriorEngineError> {
    output.clear();
    if partition.n_groups() == 0 {
        return Ok(());
    }

    let n_samples = chunk.n_samples();
    let chrom_id = chunk.chrom_id;

    // Materialise every group up front — the kernels expect an
    // iterator over owned `OverlappingVariantGroup` values and the
    // partition is already in memory, so there's no streaming win
    // from a lazy adapter. Phase A.1 will fold this back into a
    // columnar kernel that walks the partition slices directly.
    let groups: Vec<OverlappingVariantGroup> = (0..partition.n_groups())
        .map(|g| build_overlapping_variant_group(chunk, partition, g, n_samples, chrom_id))
        .collect();

    let merger = PerGroupMerger::with_config(
        groups.into_iter().map(Ok::<_, GrouperError>),
        ref_fetcher,
        per_group_config,
    );
    let engine = PosteriorEngine::with_config(merger, posterior_config);
    for record in engine {
        output.push(record?);
    }
    Ok(())
}

/// Build one [`OverlappingVariantGroup`] from group `g` of
/// `partition`. Per-position records are materialised on demand
/// from the chunk's columnar slices via
/// [`SampleColumns::materialise_record`](super::columns::SampleColumns::materialise_record).
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

/// Construct the `SharedRefFetcher` shape `PerGroupMerger` expects
/// from an owned [`ChromRefFetcher`] implementor (typically the
/// streaming fetcher held by the per-chromosome driver). Convenience
/// wrapper so the driver does not need to know the precise `Arc<dyn
/// …>` alias.
pub fn shared_ref_fetcher<F>(fetcher: F) -> SharedRefFetcher
where
    F: ChromRefFetcher + Send + 'static,
{
    Arc::new(fetcher)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::cohort_block::partition::{
        PartitionScratch, WindowPartition, partition_window,
    };
    use crate::var_calling::cohort_block::test_helpers::{loaded_chunk, record, ref_plus_alt};

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

        let fetcher: SharedRefFetcher = Arc::new(NeverCalledFetcher);
        let result = run_window(
            &chunk,
            &partition,
            fetcher,
            PerGroupMergerConfig::default(),
            PosteriorEngineConfig::default(),
            &mut output,
        );
        assert!(result.is_ok());
        assert!(output.is_empty());
    }

    /// A ChromRefFetcher whose every method panics — the
    /// empty-partition fast path must never invoke it.
    struct NeverCalledFetcher;

    impl crate::fasta::fetcher::sealed::Sealed for NeverCalledFetcher {}

    impl crate::fasta::fetcher::ChromRefFetcher for NeverCalledFetcher {
        fn length(&self) -> u32 {
            unreachable!("empty-partition path must not query the fetcher");
        }
        fn fetch(
            &self,
            _start_1based: u32,
            _length: u32,
        ) -> Result<Vec<u8>, crate::fasta::fetcher::ChromRefFetchError> {
            unreachable!("empty-partition path must not query the fetcher");
        }
        fn iter_bases<'a>(
            &'a self,
        ) -> Result<
            Box<dyn Iterator<Item = Result<u8, crate::fasta::fetcher::ChromRefFetchError>> + 'a>,
            crate::fasta::fetcher::ChromRefFetchError,
        > {
            unreachable!("empty-partition path must not query the fetcher");
        }
    }
}
