//! Phase A.1 layer 1 — column-native allele unification.
//!
//! Replaces the row-shape
//! [`unify_alleles`](crate::var_calling::per_group_merger) pipeline
//! (`project_per_position_alleles` → `admit_compound_candidates` →
//! `enforce_max_alleles`) with kernels that walk the columnar chunk
//! and partition directly. The output
//! [`UnifiedAllelesColumns`] is the same semantic content as the
//! existing private `UnifiedAlleleSet`, laid out as CSR/flat-column
//! arrays so layer 2 (per-(sample, allele) stats gather) can iterate
//! without per-emit allocation.
//!
//! **Sub-step progress.** 1.0 (types) and 1.1 (per-position
//! projection — [`project_per_position_alleles_columnar`]) have
//! landed; 1.2 (compound detection) and 1.3 (max-alleles cap) are
//! queued. Until layer 1.4 wires the assembled kernel into the
//! worker, this code is reachable only through its own unit tests;
//! the Phase A.0 row-shape adapter remains the production path.

use ahash::AHashMap;

use crate::pileup_record::ChainId;
use crate::var_calling::cohort_block::columns::MaterialisedChunk;
use crate::var_calling::cohort_block::partition::WindowPartition;

/// Column-native unified allele set for one variant group. Same
/// semantic content as the existing row-shape `UnifiedAlleleSet` in
/// [`per_group_merger`](crate::var_calling::per_group_merger), but
/// laid out for column-native consumption by downstream Phase A.1
/// layers.
///
/// **Per-allele fixed-width** (length = `n_alleles`):
/// - [`Self::is_compound`] — `true` iff the allele's bytes come
///   from a chain-anchored cross-record compound.
/// - [`Self::cap_protected`] — `true` iff the allele is exempt from
///   the `max_alleles` cap (REF + chain-anchored compounds).
/// - [`Self::cohort_count`] — sum of `num_obs` across every
///   (sample, record) source that projects to this allele, plus the
///   per-anchoring-sample chain-id intersection counts for
///   compounds.
///
/// **Per-allele variable-length** (CSR):
/// - [`Self::seq_offsets`] + [`Self::seq_bytes`] — allele byte
///   sequence projected onto the group's REF span. `alleles[0]` is
///   always REF.
/// - [`Self::constituent_offsets`] +
///   [`Self::constituent_record_idx`] +
///   [`Self::constituent_local_allele_idx`] — for compounds, the
///   `(record_idx, local_allele_idx)` pairs from the first
///   chain-anchoring sample. Empty for non-compounds.
///
/// **Per-(allele, sample) fixed-width** (length = `n_alleles ×
/// n_samples`, allele-major):
/// - [`Self::chain_anchor_counts`] — number of chain ids in sample
///   `s` whose proposal byte sequence equals allele `a`'s `seq`.
///   Zero for non-compounds and for samples that don't anchor.
///
/// **Per-(allele, sample) variable-length** (CSR over the same
/// allele-major flat order, length = `n_alleles × n_samples + 1`):
/// - [`Self::source_offsets`] +
///   [`Self::source_record_idx`] +
///   [`Self::source_local_allele_idx`] — `(record_idx,
///   local_allele_idx)` pairs in the group that project to this
///   (allele, sample) cell. Used by layer 2 to gather per-sample
///   `AlleleSupportStats`.
///
/// **OTHER pool** (length = `n_samples + 1` for the offsets):
/// - [`Self::other_offsets`] + [`Self::other_record_idx`] +
///   [`Self::other_local_allele_idx`] — per-sample
///   `(record_idx, local_allele_idx)` pairs that came from alleles
///   dropped by the `max_alleles` cap. Layer 2 folds these into the
///   per-sample `other_scalars` table.
///
/// **Lifecycle.** The driver / worker holds one persistent
/// `UnifiedAllelesColumns` and [`Self::clear`]s it before each
/// group. Sentinels (`[0]` for every CSR offset column) survive the
/// clear so the invariant `len(offsets) == n + 1` holds at every
/// step.
#[derive(Debug, Clone, PartialEq)]
pub struct UnifiedAllelesColumns {
    pub n_samples: usize,
    pub is_compound: Vec<bool>,
    pub cap_protected: Vec<bool>,
    pub cohort_count: Vec<u64>,
    pub seq_offsets: Vec<u32>,
    pub seq_bytes: Vec<u8>,
    pub constituent_offsets: Vec<u32>,
    pub constituent_record_idx: Vec<u32>,
    pub constituent_local_allele_idx: Vec<u32>,
    pub chain_anchor_counts: Vec<u32>,
    pub source_offsets: Vec<u32>,
    pub source_record_idx: Vec<u32>,
    pub source_local_allele_idx: Vec<u32>,
    pub other_offsets: Vec<u32>,
    pub other_record_idx: Vec<u32>,
    pub other_local_allele_idx: Vec<u32>,
}

impl UnifiedAllelesColumns {
    /// Empty columns with every CSR-offset column carrying its
    /// `[0]` sentinel.
    pub fn empty() -> Self {
        Self {
            n_samples: 0,
            is_compound: Vec::new(),
            cap_protected: Vec::new(),
            cohort_count: Vec::new(),
            seq_offsets: vec![0],
            seq_bytes: Vec::new(),
            constituent_offsets: vec![0],
            constituent_record_idx: Vec::new(),
            constituent_local_allele_idx: Vec::new(),
            chain_anchor_counts: Vec::new(),
            source_offsets: vec![0],
            source_record_idx: Vec::new(),
            source_local_allele_idx: Vec::new(),
            other_offsets: vec![0],
            other_record_idx: Vec::new(),
            other_local_allele_idx: Vec::new(),
        }
    }

    /// Number of unified alleles currently held.
    pub fn n_alleles(&self) -> usize {
        self.is_compound.len()
    }

    /// Reset every column to empty while preserving allocated
    /// capacity; CSR offsets revert to their single-`0` sentinel
    /// state.
    pub fn clear(&mut self) {
        self.n_samples = 0;
        self.is_compound.clear();
        self.cap_protected.clear();
        self.cohort_count.clear();
        self.seq_offsets.clear();
        self.seq_offsets.push(0);
        self.seq_bytes.clear();
        self.constituent_offsets.clear();
        self.constituent_offsets.push(0);
        self.constituent_record_idx.clear();
        self.constituent_local_allele_idx.clear();
        self.chain_anchor_counts.clear();
        self.source_offsets.clear();
        self.source_offsets.push(0);
        self.source_record_idx.clear();
        self.source_local_allele_idx.clear();
        self.other_offsets.clear();
        self.other_offsets.push(0);
        self.other_record_idx.clear();
        self.other_local_allele_idx.clear();
    }

    /// Allele bytes for allele `allele_idx` over the group's REF span.
    pub fn allele_seq(&self, allele_idx: usize) -> &[u8] {
        let lo = self.seq_offsets[allele_idx] as usize;
        let hi = self.seq_offsets[allele_idx + 1] as usize;
        &self.seq_bytes[lo..hi]
    }

    /// Number of compound constituents for allele `allele_idx`. Zero
    /// for non-compounds.
    pub fn n_constituents(&self, allele_idx: usize) -> usize {
        let lo = self.constituent_offsets[allele_idx] as usize;
        let hi = self.constituent_offsets[allele_idx + 1] as usize;
        hi - lo
    }

    /// Chain-anchor count at `(allele_idx, sample_idx)`.
    pub fn chain_anchor_count(&self, allele_idx: usize, sample_idx: usize) -> u32 {
        debug_assert!(sample_idx < self.n_samples);
        self.chain_anchor_counts[allele_idx * self.n_samples + sample_idx]
    }

    /// Half-open range of `(record_idx, local_allele_idx)` sources
    /// for `(allele_idx, sample_idx)`. Use with
    /// [`Self::source_record_idx`] / [`Self::source_local_allele_idx`]
    /// to read the pairs.
    pub fn source_range(&self, allele_idx: usize, sample_idx: usize) -> (usize, usize) {
        debug_assert!(sample_idx < self.n_samples);
        let flat = allele_idx * self.n_samples + sample_idx;
        let lo = self.source_offsets[flat] as usize;
        let hi = self.source_offsets[flat + 1] as usize;
        (lo, hi)
    }

    /// Half-open range of `(record_idx, local_allele_idx)` OTHER-pool
    /// sources for `sample_idx`.
    pub fn other_range(&self, sample_idx: usize) -> (usize, usize) {
        debug_assert!(sample_idx < self.n_samples);
        let lo = self.other_offsets[sample_idx] as usize;
        let hi = self.other_offsets[sample_idx + 1] as usize;
        (lo, hi)
    }
}

impl Default for UnifiedAllelesColumns {
    fn default() -> Self {
        Self::empty()
    }
}

/// Reusable scratch buffers for [`project_per_position_alleles_columnar`]
/// and (in later sub-steps) the full
/// `unify_alleles_columnar`. Cleared per group; capacities preserved
/// across groups.
#[derive(Debug, Default)]
pub struct UnifyAllelesScratch {
    /// Per-position projection's byte-sequence → allele-index dedup
    /// map. Cleared at entry to each per-group call.
    pub(crate) byte_index: AHashMap<Vec<u8>, usize>,
    /// Working buffer for projecting a candidate allele's bytes onto
    /// the group's REF span before testing for byte-equality against
    /// existing entries. Re-used inside the per-position projection
    /// loop and (later) by compound projection.
    pub(crate) projection_buf: Vec<u8>,
    /// Working buffer for chain-id intersection across constituent
    /// records during compound detection (sub-step 1.2; unused in
    /// 1.1).
    #[allow(dead_code)]
    pub(crate) chain_id_scratch: Vec<ChainId>,
    /// Per-(allele, sample) working source lists in allele-major flat
    /// order — `intermediate_sources[a * n_samples + s]` holds the
    /// `(record_idx, local_allele_idx)` pairs in the group that
    /// project to allele `a` for sample `s`. Serialised into
    /// [`UnifiedAllelesColumns`]'s source-CSR at the end of each call.
    /// The outer `Vec` may grow past the current group's
    /// `n_alleles × n_samples`; only the first cells are valid, the
    /// rest were cleared by the previous group's serialisation.
    pub(crate) intermediate_sources: Vec<Vec<(u32, u32)>>,
}

impl UnifyAllelesScratch {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn clear(&mut self) {
        self.byte_index.clear();
        self.projection_buf.clear();
        self.chain_id_scratch.clear();
        for cell in self.intermediate_sources.iter_mut() {
            cell.clear();
        }
        // Leave `intermediate_sources` outer length intact — the
        // inner-Vec capacities are preserved for re-use by the next
        // group, and the cells past the current group's logical
        // length are empty so they contribute nothing to the
        // serialised output.
    }
}

/// Phase A.1 sub-step 1.1 — per-position allele projection.
///
/// Native column-native port of the existing
/// `project_per_position_alleles` in
/// [`per_group_merger`](crate::var_calling::per_group_merger). For
/// each cohort position in the partition's group `group_idx`, walks
/// every sample that has a record at that position, projects each
/// allele's local sequence onto the group's REF span via the
/// walker's anchor convention, deduplicates by byte sequence, and
/// records the `(record_idx, local_allele_idx)` source under
/// `(unified_allele_idx, sample_idx)`.
///
/// **`record_idx` semantics.** `record_idx` is the index *within
/// the group's position list* — i.e. the position-iteration order
/// — not a global chunk row index. This matches the existing
/// `OverlappingVariantGroup.records[record_idx]` convention so the
/// outputs of this and the row-shape function can be byte-identity
/// compared.
///
/// **Output ordering.** REF goes at index 0 (its bytes are taken
/// from `ref_seq`). Subsequent ALTs follow in insertion order:
/// outer = position ascending, inner-1 = `sample_idx` ascending,
/// inner-2 = allele declaration order within each record. The
/// row-shape function in `per_group_merger` walks
/// `group.records.iter().enumerate()` → `pp.per_sample.iter().enumerate()`
/// → `rec.alleles.iter().enumerate()`, which produces the same
/// order because the partition stores `samples_at_pos` in
/// sample-index-ascending order (see
/// [`partition_window`](super::super::partition::partition_window)).
///
/// **No compound detection in this sub-step.** Compound alleles
/// land in sub-step 1.2 (`admit_compound_candidates_columnar`).
/// `out.is_compound` is `false` for every entry produced here;
/// `out.constituent_*` and `out.chain_anchor_counts` are emitted as
/// no-op shells (sentinels in the CSR offsets, zeros in the flat
/// grid).
///
/// **Pre-conditions.**
/// - `ref_seq.len() == partition.group_ends[group_idx] -
///   partition.group_starts[group_idx] + 1`.
/// - `n_samples == chunk.n_samples()`.
pub fn project_per_position_alleles_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    ref_seq: &[u8],
    n_samples: usize,
    scratch: &mut UnifyAllelesScratch,
    out: &mut UnifiedAllelesColumns,
) {
    out.clear();
    out.n_samples = n_samples;
    scratch.clear();

    let group_start = partition.group_starts[group_idx];
    let group_end = partition.group_ends[group_idx];
    debug_assert_eq!(
        ref_seq.len(),
        (group_end - group_start + 1) as usize,
        "ref_seq length {} does not match group span {}..{}",
        ref_seq.len(),
        group_start,
        group_end,
    );
    debug_assert_eq!(n_samples, chunk.n_samples());

    // REF goes at index 0.
    start_new_allele(out, ref_seq, false, true, n_samples, scratch);
    scratch.byte_index.insert(ref_seq.to_vec(), 0);

    let position_range = partition.position_range_for_group(group_idx);
    let group_position_lo = position_range.start;
    for p in position_range.clone() {
        let record_idx_in_group = p - group_position_lo;
        let pos = partition.positions[p];
        let local_offset = (pos - group_start) as usize;
        let sample_range = partition.sample_range_for_position(p);
        for k in sample_range {
            let sample_idx = partition.samples_at_pos[k] as usize;
            let row_idx = partition.rows_at_pos[k] as usize;
            let sample = &chunk.per_sample[sample_idx];
            let allele_lo = sample.allele_offsets[row_idx] as usize;
            let allele_hi = sample.allele_offsets[row_idx + 1] as usize;
            let local_span = sample.ref_span_at(row_idx) as usize;
            for k_allele in allele_lo..allele_hi {
                let local_allele_idx = k_allele - allele_lo;
                let seq_lo = sample.allele_seq_offsets[k_allele] as usize;
                let seq_hi = sample.allele_seq_offsets[k_allele + 1] as usize;
                let local_seq = &sample.allele_seq_bytes[seq_lo..seq_hi];

                // Project local seq onto group span: prefix from REF +
                // allele bytes + suffix from REF.
                scratch.projection_buf.clear();
                scratch
                    .projection_buf
                    .extend_from_slice(&ref_seq[..local_offset]);
                scratch.projection_buf.extend_from_slice(local_seq);
                scratch
                    .projection_buf
                    .extend_from_slice(&ref_seq[local_offset + local_span..]);

                let entry_idx = match scratch.byte_index.get(scratch.projection_buf.as_slice()) {
                    Some(&idx) => idx,
                    None => {
                        let idx = out.n_alleles();
                        // Insert into the index by cloning the bytes
                        // out of the projection buffer — the buffer
                        // will be overwritten in the next iteration.
                        scratch
                            .byte_index
                            .insert(scratch.projection_buf.clone(), idx);
                        let seq_copy = scratch.projection_buf.clone();
                        start_new_allele(out, &seq_copy, false, false, n_samples, scratch);
                        idx
                    }
                };

                // Cohort count and per-(allele, sample) source.
                out.cohort_count[entry_idx] += u64::from(sample.allele_num_obs[k_allele]);
                let flat = entry_idx * n_samples + sample_idx;
                scratch.intermediate_sources[flat]
                    .push((record_idx_in_group as u32, local_allele_idx as u32));
            }
        }
    }

    // ── Serialise intermediate_sources into out's CSR. ──
    out.source_offsets.clear();
    out.source_offsets.push(0);
    out.source_record_idx.clear();
    out.source_local_allele_idx.clear();
    let total_cells = out.n_alleles() * n_samples;
    for cell in scratch.intermediate_sources.iter().take(total_cells) {
        for &(rec_idx, local_allele_idx) in cell.iter() {
            out.source_record_idx.push(rec_idx);
            out.source_local_allele_idx.push(local_allele_idx);
        }
        out.source_offsets.push(out.source_record_idx.len() as u32);
    }

    // OTHER pool is empty in sub-step 1.1 — `n_samples` zero ranges
    // so layer-2 consumers can index without an empty-check.
    out.other_offsets.clear();
    out.other_offsets.resize(n_samples + 1, 0);
}

/// Push a fresh allele entry into `out` and the matching working
/// state into `scratch.intermediate_sources`. Used for both REF
/// (at index 0) and every newly-discovered ALT during projection.
fn start_new_allele(
    out: &mut UnifiedAllelesColumns,
    seq: &[u8],
    is_compound: bool,
    cap_protected: bool,
    n_samples: usize,
    scratch: &mut UnifyAllelesScratch,
) {
    out.is_compound.push(is_compound);
    out.cap_protected.push(cap_protected);
    out.cohort_count.push(0);
    out.seq_bytes.extend_from_slice(seq);
    out.seq_offsets.push(out.seq_bytes.len() as u32);
    // No constituents for non-compounds; the offset sentinel just
    // repeats the current end-of-constituents.
    out.constituent_offsets
        .push(out.constituent_record_idx.len() as u32);
    // Append n_samples zeros to the allele-major chain-anchor grid.
    let new_len = out.chain_anchor_counts.len() + n_samples;
    out.chain_anchor_counts.resize(new_len, 0);
    // Grow `intermediate_sources` to cover this allele's per-sample
    // slots. Existing slots past the new length (left over from a
    // bigger previous group) keep their inner-Vec capacity and stay
    // empty — they'll just be skipped during serialisation.
    let needed = out.n_alleles() * n_samples;
    if scratch.intermediate_sources.len() < needed {
        scratch.intermediate_sources.resize_with(needed, Vec::new);
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::fasta::ChromRefFetcher;
    use crate::fasta::fetcher::ChromRefFetchError;
    use crate::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};
    use crate::var_calling::cohort_block::columns::{MaterialisedChunk, SampleColumns};
    use crate::var_calling::cohort_block::partition::{
        PartitionScratch, WindowPartition, partition_window,
    };
    use crate::var_calling::cohort_block::test_helpers::{loaded_chunk, record, ref_plus_alt};
    use crate::var_calling::per_group_merger::{
        PerGroupMerger, PerGroupMergerConfig, SharedRefFetcher,
    };
    use crate::var_calling::variant_grouping::{GrouperError, OverlappingVariantGroup};

    /// In-memory `ChromRefFetcher` for the kernel tests — mirrors the
    /// `MockRef` used inside `per_group_merger`'s test module so this
    /// suite can construct the same `SharedRefFetcher` for byte-
    /// identity comparisons against the existing kernel.
    #[derive(Clone)]
    struct MockRef {
        seq: Vec<u8>,
        /// 1-based position of `seq[0]`.
        base_offset: u32,
    }

    impl crate::fasta::fetcher::sealed::Sealed for MockRef {}
    impl ChromRefFetcher for MockRef {
        fn length(&self) -> u32 {
            self.base_offset.saturating_sub(1) + self.seq.len() as u32
        }
        fn fetch(&self, start_1based: u32, length: u32) -> Result<Vec<u8>, ChromRefFetchError> {
            let start_idx = (start_1based - self.base_offset) as usize;
            Ok(self.seq[start_idx..start_idx + length as usize].to_vec())
        }
        fn iter_bases<'a>(
            &'a self,
        ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>
        {
            Ok(Box::new(self.seq.iter().copied().map(Ok)))
        }
    }

    fn shared_mock(seq: &[u8], base_offset: u32) -> SharedRefFetcher {
        Arc::new(MockRef {
            seq: seq.to_vec(),
            base_offset,
        })
    }

    /// Build a chunk + partition over the supplied per-sample
    /// records, return everything the kernel needs to project one
    /// group's alleles: chunk, partition, group_idx (always 0 in
    /// these fixtures), the group's ref_seq slice over the supplied
    /// reference.
    fn group_fixture(
        per_sample_records: Vec<Vec<PileupRecord>>,
        ref_seq_full: &[u8],
        window: std::ops::Range<u32>,
        max_group_span: u32,
    ) -> (MaterialisedChunk, WindowPartition, Vec<u8>) {
        let n = per_sample_records.len();
        let (chunk, _) = loaded_chunk(0, 1..200, per_sample_records);
        let mut scratch = PartitionScratch::with_n_samples(n);
        let mut partition = WindowPartition::empty();
        partition_window(
            &chunk,
            &window,
            &[],
            max_group_span,
            &mut scratch,
            &mut partition,
        )
        .unwrap();
        assert!(
            partition.n_groups() >= 1,
            "fixture must produce at least one variant group",
        );
        let group_start = partition.group_starts[0];
        let group_end = partition.group_ends[0];
        let ref_seq = ref_seq_full[(group_start - 1) as usize..group_end as usize].to_vec();
        (chunk, partition, ref_seq)
    }

    /// Pull the unified alleles (seq + is_compound) out of the
    /// existing row-shape kernel by feeding one
    /// `OverlappingVariantGroup` through `PerGroupMerger` and
    /// extracting `MergedRecord.alleles`. This is the byte-identity
    /// oracle for the new column-native kernel.
    fn unified_alleles_via_existing_kernel(
        chunk: &MaterialisedChunk,
        partition: &WindowPartition,
        group_idx: usize,
        ref_fetcher: SharedRefFetcher,
    ) -> Vec<(Vec<u8>, bool)> {
        let group = crate::var_calling::cohort_block::worker::build_overlapping_variant_group(
            chunk,
            partition,
            group_idx,
            chunk.n_samples(),
            chunk.chrom_id,
        );
        // Configure the merger with a generous max_alleles so the cap
        // never triggers in these fixtures.
        let config = PerGroupMergerConfig::new(2, 16, 64, 32).expect("merger config");
        let iter: Vec<Result<OverlappingVariantGroup, GrouperError>> = vec![Ok(group)];
        let mut merger = PerGroupMerger::with_config(iter.into_iter(), ref_fetcher, config);
        let item = merger
            .next()
            .expect("merger yields at least one item")
            .expect("merger succeeded");
        item.alleles
            .iter()
            .map(|a| (a.seq.clone(), a.is_compound))
            .collect()
    }

    #[test]
    fn empty_columns_carry_csr_sentinels() {
        let columns = UnifiedAllelesColumns::empty();
        assert_eq!(columns.n_alleles(), 0);
        assert_eq!(columns.n_samples, 0);
        assert_eq!(columns.seq_offsets, vec![0]);
        assert_eq!(columns.constituent_offsets, vec![0]);
        assert_eq!(columns.source_offsets, vec![0]);
        assert_eq!(columns.other_offsets, vec![0]);
    }

    #[test]
    fn clear_restores_sentinels_and_preserves_capacity() {
        let mut columns = UnifiedAllelesColumns::empty();
        // Hand-populate one allele's worth of payload (no real
        // unification yet — this just exercises the storage shape).
        columns.n_samples = 2;
        columns.is_compound.push(false);
        columns.cap_protected.push(true);
        columns.cohort_count.push(7);
        columns.seq_bytes.extend_from_slice(b"AC");
        columns.seq_offsets.push(2);
        columns.chain_anchor_counts.extend_from_slice(&[0, 0]);
        columns.source_offsets.push(1);
        columns.source_offsets.push(1);
        columns.source_record_idx.push(0);
        columns.source_local_allele_idx.push(0);
        columns.other_offsets.push(0);
        columns.other_offsets.push(0);

        let bytes_cap = columns.seq_bytes.capacity();
        columns.clear();
        assert_eq!(columns.n_alleles(), 0);
        assert_eq!(columns.n_samples, 0);
        assert_eq!(columns.seq_offsets, vec![0]);
        assert_eq!(columns.source_offsets, vec![0]);
        assert_eq!(columns.other_offsets, vec![0]);
        assert!(columns.seq_bytes.capacity() >= bytes_cap);
    }

    #[test]
    fn allele_seq_returns_slice_over_csr() {
        let mut columns = UnifiedAllelesColumns::empty();
        columns.n_samples = 0;
        columns.is_compound.push(false);
        columns.cap_protected.push(true);
        columns.cohort_count.push(0);
        columns.seq_bytes.extend_from_slice(b"AC");
        columns.seq_offsets.push(2);
        // Second allele.
        columns.is_compound.push(false);
        columns.cap_protected.push(false);
        columns.cohort_count.push(0);
        columns.seq_bytes.extend_from_slice(b"GTC");
        columns.seq_offsets.push(5);

        assert_eq!(columns.n_alleles(), 2);
        assert_eq!(columns.allele_seq(0), b"AC");
        assert_eq!(columns.allele_seq(1), b"GTC");
    }

    #[test]
    fn source_range_indexes_allele_major_grid() {
        let mut columns = UnifiedAllelesColumns::empty();
        columns.n_samples = 3;
        // Two alleles.
        columns.is_compound.extend_from_slice(&[false, false]);
        columns.cap_protected.extend_from_slice(&[true, false]);
        columns.cohort_count.extend_from_slice(&[0, 0]);
        columns.seq_bytes.extend_from_slice(b"AC");
        columns.seq_offsets.extend_from_slice(&[1, 2]);
        columns.chain_anchor_counts.extend_from_slice(&[0; 6]);
        // Source CSR — flat order (allele_idx * n_samples + sample_idx):
        //   (0,0) → 1 source
        //   (0,1) → 0 sources
        //   (0,2) → 2 sources
        //   (1,0) → 0 sources
        //   (1,1) → 1 source
        //   (1,2) → 0 sources
        columns
            .source_offsets
            .extend_from_slice(&[1, 1, 3, 3, 4, 4]);
        columns.source_record_idx.extend_from_slice(&[5, 9, 11, 17]);
        columns
            .source_local_allele_idx
            .extend_from_slice(&[0, 1, 2, 3]);

        assert_eq!(columns.source_range(0, 0), (0, 1));
        assert_eq!(columns.source_range(0, 1), (1, 1));
        assert_eq!(columns.source_range(0, 2), (1, 3));
        assert_eq!(columns.source_range(1, 0), (3, 3));
        assert_eq!(columns.source_range(1, 1), (3, 4));
        assert_eq!(columns.source_range(1, 2), (4, 4));
    }

    #[test]
    fn chain_anchor_count_indexes_allele_major_grid() {
        let mut columns = UnifiedAllelesColumns::empty();
        columns.n_samples = 3;
        columns.is_compound.extend_from_slice(&[false, true]);
        columns.cap_protected.extend_from_slice(&[true, true]);
        columns.cohort_count.extend_from_slice(&[0, 0]);
        columns.seq_offsets.extend_from_slice(&[0, 0]);
        columns
            .chain_anchor_counts
            .extend_from_slice(&[0, 0, 0, 3, 0, 5]);
        columns
            .source_offsets
            .extend_from_slice(&[0, 0, 0, 0, 0, 0]);

        assert_eq!(columns.chain_anchor_count(0, 0), 0);
        assert_eq!(columns.chain_anchor_count(1, 0), 3);
        assert_eq!(columns.chain_anchor_count(1, 2), 5);
    }

    #[test]
    fn project_single_snp_one_sample_yields_ref_plus_alt() {
        // One sample, single SNP A→T at pos 50.
        // Reference is all 'A' from position 1 onwards.
        let s0 = vec![record(50, ref_plus_alt(3, 4))];
        let ref_full = vec![b'A'; 200];
        let (chunk, partition, ref_seq) = group_fixture(vec![s0], &ref_full, 1..200, 50);

        let mut scratch = UnifyAllelesScratch::new();
        let mut out = UnifiedAllelesColumns::empty();
        project_per_position_alleles_columnar(
            &chunk,
            &partition,
            0,
            &ref_seq,
            chunk.n_samples(),
            &mut scratch,
            &mut out,
        );

        assert_eq!(out.n_alleles(), 2);
        assert_eq!(out.allele_seq(0), b"A");
        assert_eq!(out.allele_seq(1), b"T");
        assert_eq!(out.is_compound, vec![false, false]);
        assert_eq!(out.cap_protected, vec![true, false]);
        // REF cohort_count = 3 (ref_plus_alt(3, 4) puts num_obs=3 on REF
        // and num_obs=4 on ALT).
        assert_eq!(out.cohort_count, vec![3, 4]);
        // Source CSR — REF at sample 0 has one source; ALT at sample 0
        // has one source; both at (record_idx=0, local_allele_idx=0/1).
        assert_eq!(out.source_offsets, vec![0, 1, 2]);
        assert_eq!(out.source_record_idx, vec![0, 0]);
        assert_eq!(out.source_local_allele_idx, vec![0, 1]);
        // OTHER pool empty — one slot per sample (n_samples=1), plus
        // the trailing sentinel that every CSR layout carries.
        assert_eq!(out.other_offsets, vec![0, 0]);
    }

    #[test]
    fn project_snp_shared_across_samples_dedupes_and_sums_cohort_count() {
        // Two samples, both carrying A→T at pos 50 with REF
        // num_obs=2,5 and ALT num_obs=7,9 respectively.
        let s0 = vec![record(50, ref_plus_alt(2, 7))];
        let s1 = vec![record(50, ref_plus_alt(5, 9))];
        let ref_full = vec![b'A'; 200];
        let (chunk, partition, ref_seq) = group_fixture(vec![s0, s1], &ref_full, 1..200, 50);

        let mut scratch = UnifyAllelesScratch::new();
        let mut out = UnifiedAllelesColumns::empty();
        project_per_position_alleles_columnar(
            &chunk,
            &partition,
            0,
            &ref_seq,
            chunk.n_samples(),
            &mut scratch,
            &mut out,
        );

        assert_eq!(out.n_alleles(), 2);
        assert_eq!(out.cohort_count, vec![2 + 5, 7 + 9]);
        // Source CSR — allele-major, so (REF, s0), (REF, s1), (ALT, s0), (ALT, s1).
        assert_eq!(out.source_offsets, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn project_byte_identity_against_existing_kernel_simple_snp() {
        let s0 = vec![record(50, ref_plus_alt(3, 4))];
        let s1 = vec![record(50, ref_plus_alt(5, 6))];
        let ref_full = vec![b'A'; 200];
        let (chunk, partition, ref_seq) = group_fixture(vec![s0, s1], &ref_full, 1..200, 50);

        // Project via the new column-native kernel.
        let mut scratch = UnifyAllelesScratch::new();
        let mut out = UnifiedAllelesColumns::empty();
        project_per_position_alleles_columnar(
            &chunk,
            &partition,
            0,
            &ref_seq,
            chunk.n_samples(),
            &mut scratch,
            &mut out,
        );
        let new_alleles: Vec<(Vec<u8>, bool)> = (0..out.n_alleles())
            .map(|i| (out.allele_seq(i).to_vec(), out.is_compound[i]))
            .collect();

        // Oracle: feed the same data through the row-shape kernel.
        let fetcher = shared_mock(&ref_full, 1);
        let oracle = unified_alleles_via_existing_kernel(&chunk, &partition, 0, fetcher);

        assert_eq!(new_alleles, oracle);
    }

    #[test]
    fn project_byte_identity_against_existing_kernel_multi_position() {
        // Two positions joined by an MNP at pos 50 (ref_span=3),
        // plus a SNP at pos 52. No shared chain_ids → no compounds.
        let mnp_seq = b"ACA";
        let mnp_alt_seq = b"ACG";
        let mnp_record = PileupRecord::new(
            0,
            50,
            vec![
                AlleleObservation::new(
                    mnp_seq.to_vec(),
                    AlleleSupportStats::new(4, -1.0, 4, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    mnp_alt_seq.to_vec(),
                    AlleleSupportStats::new(5, -1.0, 5, 0, 0, 0, 0),
                    Vec::new(),
                ),
            ],
        );
        let snp_record = record(52, ref_plus_alt(2, 3));
        let s0 = vec![mnp_record, snp_record.clone()];
        let s1 = vec![record(50, ref_plus_alt(6, 0))];
        // Reference: A C A from pos 50..52 → "ACA". Fill the rest with 'A'.
        let mut ref_full = vec![b'A'; 200];
        ref_full[49] = b'A';
        ref_full[50] = b'C';
        ref_full[51] = b'A';
        let (chunk, partition, ref_seq) = group_fixture(vec![s0, s1], &ref_full, 1..200, 50);
        assert_eq!(partition.group_starts[0], 50);
        assert_eq!(partition.group_ends[0], 52);
        assert_eq!(ref_seq.as_slice(), b"ACA");

        let mut scratch = UnifyAllelesScratch::new();
        let mut out = UnifiedAllelesColumns::empty();
        project_per_position_alleles_columnar(
            &chunk,
            &partition,
            0,
            &ref_seq,
            chunk.n_samples(),
            &mut scratch,
            &mut out,
        );
        let new_alleles: Vec<(Vec<u8>, bool)> = (0..out.n_alleles())
            .map(|i| (out.allele_seq(i).to_vec(), out.is_compound[i]))
            .collect();

        let fetcher = shared_mock(&ref_full, 1);
        let oracle = unified_alleles_via_existing_kernel(&chunk, &partition, 0, fetcher);
        assert_eq!(new_alleles, oracle);
    }

    #[test]
    fn scratch_clear_drops_contents_preserves_capacity() {
        let mut scratch = UnifyAllelesScratch::new();
        scratch.byte_index.insert(b"AC".to_vec(), 1);
        scratch.projection_buf.extend_from_slice(b"ACGT");
        scratch.chain_id_scratch.push(7);
        let proj_cap = scratch.projection_buf.capacity();
        scratch.clear();
        assert!(scratch.byte_index.is_empty());
        assert!(scratch.projection_buf.is_empty());
        assert!(scratch.chain_id_scratch.is_empty());
        assert!(scratch.projection_buf.capacity() >= proj_cap);
    }
}
