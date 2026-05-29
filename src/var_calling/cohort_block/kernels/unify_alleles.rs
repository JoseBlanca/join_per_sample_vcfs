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
//! **Sub-step progress.** 1.0 (types), 1.1 (per-position
//! projection — [`project_per_position_alleles_columnar`]), 1.2
//! (compound detection + admission), and 1.3 (max-alleles cap —
//! composed in [`unify_alleles_columnar`]) have landed. Layer 1.4
//! wires the assembled kernel into the worker as the production
//! allele-set path; until then, this code is reachable only through
//! its own unit tests, and the Phase A.0 row-shape adapter remains
//! the production path.

use std::collections::BTreeMap;

use ahash::AHashMap;
use thiserror::Error;

use crate::pileup_record::ChainId;
use crate::var_calling::cohort_block::columns::MaterialisedChunk;
use crate::var_calling::cohort_block::partition::WindowPartition;
use crate::var_calling::per_group_merger::CompoundConstituent;

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
// Mi1: `#[non_exhaustive]` — kernel-output columns; future column
// additions land without breaking out-of-crate struct-literal sites.
#[non_exhaustive]
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

/// One working allele entry, mutated during sub-steps 1.1 and 1.2
/// and serialised into [`UnifiedAllelesColumns`] at the end. The
/// intermediate row-shape representation lets sub-step 1.2
/// (`admit_compound_candidates_columnar`) toggle `is_compound` and
/// append `constituents` on an already-written allele without
/// having to re-pack the output CSR offsets.
#[derive(Debug, Default)]
pub(crate) struct WorkingAllele {
    pub(crate) seq: Vec<u8>,
    pub(crate) is_compound: bool,
    pub(crate) cap_protected: bool,
    pub(crate) cohort_count: u64,
    /// Sorted by `(record_idx, local_allele_idx)`. Empty for
    /// non-compounds.
    pub(crate) constituents: Vec<CompoundConstituent>,
    /// Per-sample source pairs from per-position projection. Length
    /// `n_samples`; each inner `Vec` holds the
    /// `(record_idx, local_allele_idx)` pairs that project to this
    /// allele for that sample.
    pub(crate) sources_per_sample: Vec<Vec<(u32, u32)>>,
    /// Per-sample chain-anchor count for compounds. Length
    /// `n_samples`; zero for non-compounds.
    pub(crate) chain_anchor_counts: Vec<u32>,
}

impl WorkingAllele {
    /// Reset this allele in place, preserving every inner `Vec`'s
    /// allocation. Used by the scratch pool when a previously-used
    /// slot is recycled by a new group.
    pub(crate) fn reset(&mut self, n_samples: usize) {
        self.seq.clear();
        self.is_compound = false;
        self.cap_protected = false;
        self.cohort_count = 0;
        self.constituents.clear();
        for cell in self.sources_per_sample.iter_mut() {
            cell.clear();
        }
        if self.sources_per_sample.len() < n_samples {
            self.sources_per_sample.resize_with(n_samples, Vec::new);
        }
        self.chain_anchor_counts.clear();
        self.chain_anchor_counts.resize(n_samples, 0);
    }
}

/// Reusable scratch buffers for the column-native allele unifier.
/// Cleared per group; capacities preserved across groups.
#[derive(Debug, Default)]
pub struct UnifyAllelesScratch {
    /// Per-position projection's byte-sequence → working-allele-index
    /// dedup map. Cleared at entry to each per-group call.
    pub(crate) byte_index: AHashMap<Vec<u8>, usize>,
    /// Working buffer for projecting a candidate allele's bytes onto
    /// the group's REF span before testing for byte-equality against
    /// existing entries.
    pub(crate) projection_buf: Vec<u8>,
    /// Pool of working alleles. The outer `Vec` grows monotonically
    /// across the lifetime of the scratch; [`Self::n_active_alleles`]
    /// is the logical length used by the current group. Slots past
    /// that index were cleared by [`Self::clear`] and stay
    /// allocation-preserved for re-use.
    pub(crate) working_alleles: Vec<WorkingAllele>,
    /// Logical length of [`Self::working_alleles`] for the current
    /// group. Serialisation reads exactly the kept indices
    /// (see [`Self::kept_indices`]).
    pub(crate) n_active_alleles: usize,
    /// Indices into [`Self::working_alleles`] that survive the
    /// max-alleles cap (sub-step 1.3) in original order. Populated by
    /// [`enforce_max_alleles_columnar`]; for groups where the cap
    /// did not fire, this is simply `0..n_active_alleles`. The
    /// serialiser uses this list, not `0..n_active_alleles`, to
    /// write the output columns.
    pub(crate) kept_indices: Vec<usize>,
    /// Per-sample OTHER pool — pooled `(record_idx, local_allele_idx)`
    /// pairs from cap-dropped alleles. Length is `n_samples` after
    /// each call; empty when the cap did not drop anything.
    pub(crate) other_per_sample: Vec<Vec<(u32, u32)>>,
}

impl UnifyAllelesScratch {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn clear(&mut self) {
        self.byte_index.clear();
        self.projection_buf.clear();
        // Reset every previously-used working-allele slot in place so
        // inner `Vec` capacities survive to the next group.
        let prev_active = self.n_active_alleles;
        for slot in self.working_alleles.iter_mut().take(prev_active) {
            slot.reset(0);
        }
        self.n_active_alleles = 0;
        self.kept_indices.clear();
        for cell in self.other_per_sample.iter_mut() {
            cell.clear();
        }
    }

    /// Return a `&mut` to the next working-allele slot, growing the
    /// pool if necessary and resetting the slot to a fresh state.
    /// Increments [`Self::n_active_alleles`] by one.
    pub(crate) fn push_working_allele(&mut self, n_samples: usize) -> &mut WorkingAllele {
        if self.n_active_alleles == self.working_alleles.len() {
            self.working_alleles.push(WorkingAllele::default());
        }
        let idx = self.n_active_alleles;
        self.n_active_alleles += 1;
        let slot = &mut self.working_alleles[idx];
        slot.reset(n_samples);
        slot
    }
}

/// Errors surfaced by [`admit_compound_candidates_columnar`].
#[non_exhaustive]
#[derive(Error, Debug, PartialEq)]
pub enum UnifyAllelesError {
    /// While projecting a chain-anchored compound's bytes onto the
    /// group span, no sample at the constituent's position carried a
    /// record with the named `local_allele_idx`. The walker invariant
    /// guarantees at least one anchoring sample; surfacing this
    /// signals a partition / chunk inconsistency, not a data
    /// condition.
    #[error(
        "compound projection: no sample has constituent (record_idx={record_idx}, \
         local_allele_idx={local_allele_idx}) in group \
         {chrom_id}:{group_start}..{group_end}"
    )]
    MissingCompoundAlleleBytes {
        chrom_id: u32,
        group_start: u32,
        group_end: u32,
        record_idx: u32,
        local_allele_idx: u32,
    },
}

/// Serialise the working alleles in `scratch` into the column-native
/// output. Called at the end of each per-group unification by
/// callers that built up `WorkingAllele` entries during projection
/// + compound admission + the max-alleles cap.
///
/// Iterates [`UnifyAllelesScratch::kept_indices`] (not the full
/// pool) so cap-dropped alleles are excluded from the output;
/// [`UnifyAllelesScratch::other_per_sample`] is written into the
/// OTHER columns at the end. When the cap did not fire,
/// `other_per_sample` is empty and the OTHER offset column gets
/// `n_samples + 1` zero entries so layer-2 consumers can index
/// without an empty-check.
pub(crate) fn serialize_working_to_columns(
    scratch: &UnifyAllelesScratch,
    n_samples: usize,
    out: &mut UnifiedAllelesColumns,
) {
    out.clear();
    out.n_samples = n_samples;
    for &i in &scratch.kept_indices {
        let working = &scratch.working_alleles[i];
        out.is_compound.push(working.is_compound);
        out.cap_protected.push(working.cap_protected);
        out.cohort_count.push(working.cohort_count);
        out.seq_bytes.extend_from_slice(&working.seq);
        out.seq_offsets.push(out.seq_bytes.len() as u32);
        for c in &working.constituents {
            out.constituent_record_idx.push(c.record_idx as u32);
            out.constituent_local_allele_idx
                .push(c.local_allele_idx as u32);
        }
        out.constituent_offsets
            .push(out.constituent_record_idx.len() as u32);
        out.chain_anchor_counts
            .extend_from_slice(&working.chain_anchor_counts);
    }
    // Source CSR — walk (allele, sample) in allele-major order over
    // the kept alleles.
    for &i in &scratch.kept_indices {
        let working = &scratch.working_alleles[i];
        for cell in working.sources_per_sample.iter().take(n_samples) {
            for &(rec, loc) in cell {
                out.source_record_idx.push(rec);
                out.source_local_allele_idx.push(loc);
            }
            out.source_offsets.push(out.source_record_idx.len() as u32);
        }
    }
    // OTHER pool CSR — per-sample sources from cap-dropped alleles.
    if scratch.other_per_sample.is_empty() {
        out.other_offsets.resize(n_samples + 1, 0);
    } else {
        for cell in scratch.other_per_sample.iter().take(n_samples) {
            for &(rec, loc) in cell {
                out.other_record_idx.push(rec);
                out.other_local_allele_idx.push(loc);
            }
            out.other_offsets.push(out.other_record_idx.len() as u32);
        }
    }
}

/// Phase A.1 sub-step 1.3 — max-alleles cap. Populates
/// [`UnifyAllelesScratch::kept_indices`] with the surviving
/// allele indices in original order, and
/// [`UnifyAllelesScratch::other_per_sample`] with the per-sample
/// `(record_idx, local_allele_idx)` pairs pooled from cap-dropped
/// alleles. Mirrors `enforce_max_alleles` in
/// [`per_group_merger`](crate::var_calling::per_group_merger):
/// protected alleles (REF + chain-anchored compounds) always
/// survive; among prunable alleles the top
/// `max_alleles - protected_count` by `cohort_count` survive.
pub(crate) fn enforce_max_alleles_columnar(
    scratch: &mut UnifyAllelesScratch,
    n_samples: usize,
    max_alleles: usize,
) {
    let n_active = scratch.n_active_alleles;

    // Default: keep all in original order.
    scratch.kept_indices.clear();
    scratch.kept_indices.extend(0..n_active);

    if n_active <= max_alleles {
        return;
    }

    // Partition into prunable indices (sorted by descending
    // cohort_count) and protected count.
    let mut prunable: Vec<usize> = (0..n_active)
        .filter(|&i| !scratch.working_alleles[i].cap_protected)
        .collect();
    prunable.sort_by_key(|&i| std::cmp::Reverse(scratch.working_alleles[i].cohort_count));
    let protected_count = n_active - prunable.len();
    let budget_for_prunable = max_alleles.saturating_sub(protected_count);

    if prunable.len() <= budget_for_prunable {
        return;
    }

    // Mark cap-dropped indices.
    let to_remove_sorted = {
        let mut v: Vec<usize> = prunable[budget_for_prunable..].to_vec();
        v.sort_unstable();
        v
    };

    // Rebuild kept_indices in original order, skipping dropped.
    scratch.kept_indices.clear();
    {
        let mut drop_cursor = 0;
        for i in 0..n_active {
            if drop_cursor < to_remove_sorted.len() && to_remove_sorted[drop_cursor] == i {
                drop_cursor += 1;
                continue;
            }
            scratch.kept_indices.push(i);
        }
    }

    // Build OTHER pool from dropped alleles' per-sample sources.
    if scratch.other_per_sample.len() < n_samples {
        scratch.other_per_sample.resize_with(n_samples, Vec::new);
    }
    for cell in scratch.other_per_sample.iter_mut().take(n_samples) {
        cell.clear();
    }
    for &dropped in &to_remove_sorted {
        let working = &scratch.working_alleles[dropped];
        for (s, sources) in working
            .sources_per_sample
            .iter()
            .enumerate()
            .take(n_samples)
        {
            scratch.other_per_sample[s].extend(sources.iter().copied());
        }
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
    scratch.clear();
    project_per_position_into_scratch(chunk, partition, group_idx, ref_seq, n_samples, scratch);
    // No cap in this sub-step — kept_indices is just the full active
    // range so the serialiser writes everything.
    scratch.kept_indices.clear();
    scratch.kept_indices.extend(0..scratch.n_active_alleles);
    serialize_working_to_columns(scratch, n_samples, out);
}

/// Fill `scratch.working_alleles` with the per-position projection
/// (REF at index 0, then ALTs in insertion order). Separated from
/// the serialisation step so sub-step 1.2 can mutate the working
/// alleles before they're written to the output CSR.
pub(crate) fn project_per_position_into_scratch(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    ref_seq: &[u8],
    n_samples: usize,
    scratch: &mut UnifyAllelesScratch,
) {
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
    push_allele_into_scratch(scratch, ref_seq, false, true, n_samples);
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
                        let idx = scratch.n_active_alleles;
                        let seq_copy = scratch.projection_buf.clone();
                        scratch.byte_index.insert(seq_copy.clone(), idx);
                        push_allele_into_scratch(scratch, &seq_copy, false, false, n_samples);
                        idx
                    }
                };

                let working = &mut scratch.working_alleles[entry_idx];
                working.cohort_count += u64::from(sample.allele_num_obs[k_allele]);
                working.sources_per_sample[sample_idx]
                    .push((record_idx_in_group as u32, local_allele_idx as u32));
            }
        }
    }
}

/// Push a fresh working allele entry with the given seq + flags.
/// Inner sample-indexed Vecs are pre-sized to `n_samples`.
fn push_allele_into_scratch(
    scratch: &mut UnifyAllelesScratch,
    seq: &[u8],
    is_compound: bool,
    cap_protected: bool,
    n_samples: usize,
) {
    let working = scratch.push_working_allele(n_samples);
    working.seq.extend_from_slice(seq);
    working.is_compound = is_compound;
    working.cap_protected = cap_protected;
}

// ============================================================
// Sub-step 1.2 — compound detection and admission
// ============================================================

/// Per-sample chain-anchor evidence for one candidate compound.
#[derive(Default, Clone)]
struct ChainAnchorEvidence {
    /// Number of distinct chain ids in this sample whose constituents
    /// match the candidate's constituent tuple.
    intersection: u32,
    /// `(record_idx, local_allele_idx)` pairs for this sample's
    /// constituent alleles in the compound. Populated only when
    /// `intersection > 0`. Used by layer 2's scalar-subtraction
    /// pass.
    constituent_sources: Vec<(u32, u32)>,
}

/// A compound candidate proposed by chain-id evidence in at least
/// one sample. Mirrors the private `CompoundCandidate` in
/// `per_group_merger` so the column-native algorithm is line-by-line
/// comparable.
struct CompoundCandidate {
    /// Constituents from the lex-smallest sample that proposed this
    /// candidate, sorted by `(record_idx, local_allele_idx)`.
    constituents_per_first_anchor: Vec<CompoundConstituent>,
    /// Per-sample evidence indexed by `sample_idx`.
    per_sample: Vec<ChainAnchorEvidence>,
}

/// Sub-step 1.2 — detect compound candidates from chain-id evidence
/// and fold them into the working alleles already produced by
/// sub-step 1.1.
///
/// Walks per-sample chain proposals: a chain id with ≥ 2 constituent
/// non-REF allele observations across different records proposes a
/// compound. Each candidate is keyed by its sorted constituent
/// tuple; matching proposals from multiple samples aggregate into
/// the same `CompoundCandidate`.
///
/// For each candidate, projects its bytes onto the group span. If
/// the bytes match an existing working allele (typically a
/// per-position projection), flag that entry as compound and record
/// the constituents. Otherwise, push a new compound entry.
///
/// **Pre-condition.** Sub-step 1.1
/// (`project_per_position_into_scratch`) must have populated
/// `scratch.working_alleles` before this is called — the candidate
/// admission re-uses the same `byte_index` and works on the
/// scratch's working alleles.
pub(crate) fn admit_compound_candidates_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    ref_seq: &[u8],
    n_samples: usize,
    scratch: &mut UnifyAllelesScratch,
) -> Result<(), UnifyAllelesError> {
    let candidates = detect_compound_candidates_columnar(chunk, partition, group_idx, n_samples);
    let group_start = partition.group_starts[group_idx];
    let group_end = partition.group_ends[group_idx];
    let chrom_id = chunk.chrom_id;

    for candidate in candidates {
        project_compound_onto_group_columnar(
            chunk,
            partition,
            group_idx,
            ref_seq,
            &candidate.constituents_per_first_anchor,
            &mut scratch.projection_buf,
            chrom_id,
            group_start,
            group_end,
        )?;

        let target_idx = match scratch.byte_index.get(scratch.projection_buf.as_slice()) {
            Some(&idx) => {
                let working = &mut scratch.working_alleles[idx];
                if !working.is_compound {
                    // Row-shape parity (per_group_merger.rs:1101–1106).
                    // The compound's bytes match an existing per-position
                    // allele; mark it as a compound and attach the
                    // constituents, but leave `cap_protected` alone — the
                    // entry stays prunable by `enforce_max_alleles` so
                    // its per-position observations can be dropped if
                    // its `cohort_count` is too low.
                    working.is_compound = true;
                    working.constituents = candidate.constituents_per_first_anchor.clone();
                }
                idx
            }
            None => {
                let idx = scratch.n_active_alleles;
                let seq_copy = scratch.projection_buf.clone();
                scratch.byte_index.insert(seq_copy.clone(), idx);
                push_allele_into_scratch(scratch, &seq_copy, true, true, n_samples);
                scratch.working_alleles[idx].constituents =
                    candidate.constituents_per_first_anchor.clone();
                idx
            }
        };

        let working = &mut scratch.working_alleles[target_idx];
        for (sample_idx, anchor) in candidate.per_sample.iter().enumerate() {
            working.chain_anchor_counts[sample_idx] = anchor.intersection;
            working.cohort_count += u64::from(anchor.intersection);
            if anchor.intersection > 0 {
                working.sources_per_sample[sample_idx].clone_from(&anchor.constituent_sources);
            }
        }
    }
    Ok(())
}

/// Build the per-sample chain proposals for one group and aggregate
/// them by constituent tuple into the candidate list. Iterates
/// samples in ascending `sample_idx` so the "first anchoring sample"
/// is deterministically the lex-smallest one.
fn detect_compound_candidates_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    n_samples: usize,
) -> Vec<CompoundCandidate> {
    let mut candidates: BTreeMap<Vec<(usize, usize)>, CompoundCandidate> = BTreeMap::new();
    let mut by_chain: BTreeMap<ChainId, Vec<CompoundConstituent>> = BTreeMap::new();
    for sample_idx in 0..n_samples {
        by_chain.clear();
        build_chain_proposals_columnar(chunk, partition, group_idx, sample_idx, &mut by_chain);
        for constituents in by_chain.values() {
            if constituents.len() < 2 {
                continue;
            }
            // Already sorted by (record_idx, local_allele_idx) below.
            let key: Vec<(usize, usize)> = constituents
                .iter()
                .map(|c| (c.record_idx, c.local_allele_idx))
                .collect();
            let entry = candidates.entry(key).or_insert_with(|| CompoundCandidate {
                constituents_per_first_anchor: constituents.clone(),
                per_sample: vec![ChainAnchorEvidence::default(); n_samples],
            });
            entry.per_sample[sample_idx].intersection += 1;
            if entry.per_sample[sample_idx].constituent_sources.is_empty() {
                entry.per_sample[sample_idx].constituent_sources = constituents
                    .iter()
                    .map(|c| (c.record_idx as u32, c.local_allele_idx as u32))
                    .collect();
            }
        }
    }
    candidates.into_values().collect()
}

/// Walk the group's positions for `sample_idx` and bucket each non-
/// REF allele's `chain_ids` into the `by_chain` map (one entry per
/// distinct chain id). Constituents are sorted by `(record_idx,
/// local_allele_idx)` at the end so the candidate-keying step sees
/// canonical tuples.
fn build_chain_proposals_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    sample_idx: usize,
    by_chain: &mut BTreeMap<ChainId, Vec<CompoundConstituent>>,
) {
    let position_range = partition.position_range_for_group(group_idx);
    let group_position_lo = position_range.start;
    for p in position_range.clone() {
        let record_idx = p - group_position_lo;
        let sample_range = partition.sample_range_for_position(p);
        let samples_slice = &partition.samples_at_pos[sample_range.clone()];
        let rows_slice = &partition.rows_at_pos[sample_range];
        let Some(local_slot) = samples_slice.iter().position(|&s| s as usize == sample_idx) else {
            continue;
        };
        let row_idx = rows_slice[local_slot] as usize;
        let sample = &chunk.per_sample[sample_idx];
        let allele_lo = sample.allele_offsets[row_idx] as usize;
        let allele_hi = sample.allele_offsets[row_idx + 1] as usize;
        // Skip allele 0 (REF) — never participates in a compound.
        for k_allele in (allele_lo + 1)..allele_hi {
            let local_allele_idx = k_allele - allele_lo;
            let chain_lo = sample.allele_chain_ids_offsets[k_allele] as usize;
            let chain_hi = sample.allele_chain_ids_offsets[k_allele + 1] as usize;
            for &chain_id in &sample.allele_chain_ids[chain_lo..chain_hi] {
                let entry = by_chain.entry(chain_id).or_default();
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
    // Sort each chain's constituents — canonical tuple ordering.
    for constituents in by_chain.values_mut() {
        constituents.sort_by_key(|c| (c.record_idx, c.local_allele_idx));
    }
}

/// Project one compound candidate's bytes onto the group's REF
/// span. Constituents are assumed sorted by `record_idx`. The
/// algorithm walks them in order, substituting each one's local
/// allele bytes at its offset within the group span; non-overlapping
/// constituents (the typical case) just concatenate, overlapping
/// ones fall back to a last-write-wins composition (same rule as
/// the existing `project_compound_onto_group`).
#[allow(clippy::too_many_arguments)]
fn project_compound_onto_group_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    ref_seq: &[u8],
    constituents: &[CompoundConstituent],
    out: &mut Vec<u8>,
    chrom_id: u32,
    group_start: u32,
    group_end: u32,
) -> Result<(), UnifyAllelesError> {
    out.clear();
    if constituents.is_empty() {
        out.extend_from_slice(ref_seq);
        return Ok(());
    }

    debug_assert!(
        constituents
            .windows(2)
            .all(|w| w[0].record_idx <= w[1].record_idx),
        "constituents must be sorted by record_idx",
    );

    let position_range = partition.position_range_for_group(group_idx);
    let group_position_lo = position_range.start;
    let mut cursor_offset: usize = 0;

    for c in constituents {
        let p = group_position_lo + c.record_idx;
        let pos = partition.positions[p];
        let local_offset = (pos - group_start) as usize;
        let sample_range = partition.sample_range_for_position(p);

        let mut local_span: Option<usize> = None;
        let mut found_seq: Option<(usize, std::ops::Range<usize>)> = None;
        for k in sample_range {
            let sample_idx = partition.samples_at_pos[k] as usize;
            let row_idx = partition.rows_at_pos[k] as usize;
            let sample = &chunk.per_sample[sample_idx];
            if local_span.is_none() {
                local_span = Some(sample.ref_span_at(row_idx) as usize);
            }
            let allele_lo = sample.allele_offsets[row_idx] as usize;
            let allele_hi = sample.allele_offsets[row_idx + 1] as usize;
            if c.local_allele_idx < (allele_hi - allele_lo) {
                let k_allele = allele_lo + c.local_allele_idx;
                let seq_lo = sample.allele_seq_offsets[k_allele] as usize;
                let seq_hi = sample.allele_seq_offsets[k_allele + 1] as usize;
                found_seq = Some((sample_idx, seq_lo..seq_hi));
                break;
            }
        }
        let (sample_idx, seq_range) =
            found_seq.ok_or(UnifyAllelesError::MissingCompoundAlleleBytes {
                chrom_id,
                group_start,
                group_end,
                record_idx: c.record_idx as u32,
                local_allele_idx: c.local_allele_idx as u32,
            })?;
        let local_span = local_span.unwrap_or(1);

        if local_offset < cursor_offset {
            let overlap = cursor_offset - local_offset;
            out.truncate(out.len().saturating_sub(overlap));
            cursor_offset = local_offset;
        }
        if local_offset > cursor_offset {
            out.extend_from_slice(&ref_seq[cursor_offset..local_offset]);
        }
        out.extend_from_slice(&chunk.per_sample[sample_idx].allele_seq_bytes[seq_range]);
        cursor_offset = local_offset + local_span;
    }

    if cursor_offset < ref_seq.len() {
        out.extend_from_slice(&ref_seq[cursor_offset..]);
    }
    Ok(())
}

/// Full column-native allele unifier — composes sub-steps 1.1
/// (per-position projection), 1.2 (compound detection + admission),
/// and 1.3 (max-alleles cap).
///
/// `max_alleles` is the cap on the kept-allele count. When the cap
/// fires, the lowest-`cohort_count` non-protected alleles are
/// dropped and their per-sample `(record_idx, local_allele_idx)`
/// sources are folded into the OTHER pool (written to the
/// `other_*` columns of `out`). REF and chain-anchored compounds
/// are protected.
#[allow(clippy::too_many_arguments)]
pub fn unify_alleles_columnar(
    chunk: &MaterialisedChunk,
    partition: &WindowPartition,
    group_idx: usize,
    ref_seq: &[u8],
    n_samples: usize,
    max_alleles: usize,
    scratch: &mut UnifyAllelesScratch,
    out: &mut UnifiedAllelesColumns,
) -> Result<(), UnifyAllelesError> {
    scratch.clear();
    project_per_position_into_scratch(chunk, partition, group_idx, ref_seq, n_samples, scratch);
    admit_compound_candidates_columnar(chunk, partition, group_idx, ref_seq, n_samples, scratch)?;
    enforce_max_alleles_columnar(scratch, n_samples, max_alleles);
    serialize_working_to_columns(scratch, n_samples, out);
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::fasta::ChromRefFetcher;
    use crate::fasta::fetcher::ChromRefFetchError;
    use crate::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};
    use crate::var_calling::cohort_block::columns::MaterialisedChunk;
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

    /// Pull the existing kernel's allele set including `constituents`
    /// so we can byte-identity-compare compound detection too.
    fn unified_alleles_with_constituents(
        chunk: &MaterialisedChunk,
        partition: &WindowPartition,
        group_idx: usize,
        ref_fetcher: SharedRefFetcher,
    ) -> Vec<(Vec<u8>, bool, Vec<(usize, usize)>)> {
        let group = crate::var_calling::cohort_block::worker::build_overlapping_variant_group(
            chunk,
            partition,
            group_idx,
            chunk.n_samples(),
            chunk.chrom_id,
        );
        let config = PerGroupMergerConfig::new(2, 16, 64, 32).expect("merger config");
        let iter: Vec<Result<OverlappingVariantGroup, GrouperError>> = vec![Ok(group)];
        let mut merger = PerGroupMerger::with_config(iter.into_iter(), ref_fetcher, config);
        let record = merger
            .next()
            .expect("merger yields at least one item")
            .expect("merger succeeded");
        record
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
            .collect()
    }

    /// Pull the new column-native kernel's allele set including
    /// `constituents` (the same shape the row-shape oracle returns).
    fn columnar_alleles_with_constituents(
        out: &UnifiedAllelesColumns,
    ) -> Vec<(Vec<u8>, bool, Vec<(usize, usize)>)> {
        (0..out.n_alleles())
            .map(|i| {
                let lo = out.constituent_offsets[i] as usize;
                let hi = out.constituent_offsets[i + 1] as usize;
                let constituents: Vec<(usize, usize)> = (lo..hi)
                    .map(|k| {
                        (
                            out.constituent_record_idx[k] as usize,
                            out.constituent_local_allele_idx[k] as usize,
                        )
                    })
                    .collect();
                (out.allele_seq(i).to_vec(), out.is_compound[i], constituents)
            })
            .collect()
    }

    #[test]
    fn unify_byte_identity_against_existing_kernel_with_compound() {
        // Sample 0 has a chain-anchored compound: SNP A→T at pos 100
        // (chain 1) and SNP A→G at pos 102 (chain 1). Same chain
        // means the same read carries both ALTs → a compound.
        // ref_span = 1 (SNPs) but the group spans [100, 102] because
        // we'll attach another MNP somewhere to glue the positions.
        //
        // To get them in one group, sample 1 carries an MNP at pos
        // 100 (REF=ACA, ALT=ACC) with ref_span=3, reaching to 102.
        // That joins 100 and 102 into one group.
        let mnp_record = PileupRecord::new(
            0,
            100,
            vec![
                AlleleObservation::new(
                    b"ACA".to_vec(),
                    AlleleSupportStats::new(4, -1.0, 4, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    b"ACC".to_vec(),
                    AlleleSupportStats::new(3, -1.0, 3, 0, 0, 0, 0),
                    vec![100],
                ),
            ],
        );
        // Sample 0: two SNPs with the same chain_id 1 — proposes a
        // compound.
        let snp_100 = PileupRecord::new(
            0,
            100,
            vec![
                AlleleObservation::new(
                    b"A".to_vec(),
                    AlleleSupportStats::new(5, -1.0, 5, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    b"T".to_vec(),
                    AlleleSupportStats::new(4, -1.0, 4, 0, 0, 0, 0),
                    vec![1],
                ),
            ],
        );
        let snp_102 = PileupRecord::new(
            0,
            102,
            vec![
                AlleleObservation::new(
                    b"A".to_vec(),
                    AlleleSupportStats::new(5, -1.0, 5, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    b"G".to_vec(),
                    AlleleSupportStats::new(4, -1.0, 4, 0, 0, 0, 0),
                    vec![1],
                ),
            ],
        );
        let s0 = vec![snp_100, snp_102];
        let s1 = vec![mnp_record];
        // Reference at positions 100..102 is "ACA" (matches the MNP's
        // REF allele).
        let mut ref_full = vec![b'A'; 200];
        ref_full[99] = b'A';
        ref_full[100] = b'C';
        ref_full[101] = b'A';
        let (chunk, partition, ref_seq) = group_fixture(vec![s0, s1], &ref_full, 1..200, 50);
        assert_eq!(partition.group_starts[0], 100);
        assert_eq!(partition.group_ends[0], 102);
        assert_eq!(ref_seq.as_slice(), b"ACA");

        // Column-native unify (1.1 + 1.2).
        let mut scratch = UnifyAllelesScratch::new();
        let mut out = UnifiedAllelesColumns::empty();
        unify_alleles_columnar(
            &chunk,
            &partition,
            0,
            &ref_seq,
            chunk.n_samples(),
            16, // max_alleles — generous so cap never fires
            &mut scratch,
            &mut out,
        )
        .expect("unify succeeded");
        let new_alleles = columnar_alleles_with_constituents(&out);

        // Oracle (row-shape kernel).
        let fetcher = shared_mock(&ref_full, 1);
        let oracle = unified_alleles_with_constituents(&chunk, &partition, 0, fetcher);

        assert_eq!(new_alleles, oracle);
    }

    #[test]
    fn unify_byte_identity_against_existing_kernel_no_compounds() {
        // Recap: when no chain ids are shared, unify_alleles_columnar
        // should produce the same allele set as the per-position-only
        // sub-step.
        let s0 = vec![record(50, ref_plus_alt(3, 4))];
        let s1 = vec![record(50, ref_plus_alt(5, 6))];
        let ref_full = vec![b'A'; 200];
        let (chunk, partition, ref_seq) = group_fixture(vec![s0, s1], &ref_full, 1..200, 50);

        let mut scratch = UnifyAllelesScratch::new();
        let mut out = UnifiedAllelesColumns::empty();
        unify_alleles_columnar(
            &chunk,
            &partition,
            0,
            &ref_seq,
            chunk.n_samples(),
            16, // max_alleles — generous so cap never fires
            &mut scratch,
            &mut out,
        )
        .expect("unify succeeded");
        let new_alleles = columnar_alleles_with_constituents(&out);

        let fetcher = shared_mock(&ref_full, 1);
        let oracle = unified_alleles_with_constituents(&chunk, &partition, 0, fetcher);

        assert_eq!(new_alleles, oracle);
    }

    #[test]
    fn unify_max_alleles_cap_drops_lowest_cohort_count_prunables() {
        // Four ALTs at the same position with descending cohort_count;
        // max_alleles = 2 (REF + 1 ALT). Only the highest-count ALT
        // survives; the others' (record_idx, local_allele_idx) pairs
        // go into the OTHER pool.
        let alts = vec![
            AlleleObservation::new(
                b"A".to_vec(),
                AlleleSupportStats::new(10, -1.0, 10, 0, 0, 0, 0),
                Vec::new(),
            ),
            AlleleObservation::new(
                b"T".to_vec(),
                AlleleSupportStats::new(8, -1.0, 8, 0, 0, 0, 0),
                Vec::new(),
            ),
            AlleleObservation::new(
                b"C".to_vec(),
                AlleleSupportStats::new(3, -1.0, 3, 0, 0, 0, 0),
                Vec::new(),
            ),
            AlleleObservation::new(
                b"G".to_vec(),
                AlleleSupportStats::new(2, -1.0, 2, 0, 0, 0, 0),
                Vec::new(),
            ),
            AlleleObservation::new(
                b"N".to_vec(),
                AlleleSupportStats::new(1, -1.0, 1, 0, 0, 0, 0),
                Vec::new(),
            ),
        ];
        let s0 = vec![PileupRecord::new(0, 50, alts)];
        let ref_full = vec![b'A'; 200];
        let (chunk, partition, ref_seq) = group_fixture(vec![s0], &ref_full, 1..200, 50);

        let mut scratch = UnifyAllelesScratch::new();
        let mut out = UnifiedAllelesColumns::empty();
        unify_alleles_columnar(
            &chunk,
            &partition,
            0,
            &ref_seq,
            chunk.n_samples(),
            2, // max_alleles — REF + 1 ALT
            &mut scratch,
            &mut out,
        )
        .expect("unify succeeded");

        // Cap survivors: REF (allele 0, A — the first non-REF allele
        // in the record's alleles array is "A" which dedupes into REF)
        // … wait — the first allele in the PileupRecord's alleles
        // array is the REF, and in this fixture I wrote it as "A".
        // The actual REF comes from `ref_seq[0]` which is also "A",
        // so they dedupe.
        //
        // After dedup REF=A counts: REF own (10 in alleles[0]) +
        // alleles[1..] that project to A; here alleles[1]="T" projects
        // to T (different), so no other allele projects to A.
        //
        // Among non-REF projections (T, C, G, N — each with one sample
        // source), the highest is T with cohort_count=8.
        assert_eq!(out.n_alleles(), 2);
        assert_eq!(out.allele_seq(0), b"A");
        assert_eq!(out.cap_protected[0], true);
        assert_eq!(out.allele_seq(1), b"T");
        assert_eq!(out.cap_protected[1], false);

        // OTHER pool: 3 dropped ALTs (C, G, N), each carrying one
        // source pair (record_idx=0, local_allele_idx in 2..5).
        let (lo, hi) = out.other_range(0);
        let other_pairs: Vec<(u32, u32)> = (lo..hi)
            .map(|k| (out.other_record_idx[k], out.other_local_allele_idx[k]))
            .collect();
        // Sort for stable comparison — the pool order isn't externally
        // contracted, only the set is.
        let mut sorted = other_pairs.clone();
        sorted.sort();
        assert_eq!(sorted, vec![(0, 2), (0, 3), (0, 4)]);
    }

    #[test]
    fn unify_max_alleles_cap_protects_compound_even_with_low_count() {
        // Sample 0 has two SNPs sharing chain_id=1 → compound (low
        // cohort_count). Sample 1 has many high-count single SNPs.
        // With max_alleles=3 the compound stays (protected) even
        // though its count is lower than the dropped SNPs.
        //
        // Build via: sample 0 has SNPs at pos 100 (REF=A,ALT=T,
        // chain=1) and pos 102 (REF=A,ALT=G, chain=1). Sample 1 has
        // an MNP at pos 100 ref_span=3 to join positions, plus extra
        // ALTs at pos 100 that compete on cohort_count.
        let mnp_rec = PileupRecord::new(
            0,
            100,
            vec![
                AlleleObservation::new(
                    b"ACA".to_vec(),
                    AlleleSupportStats::new(20, -1.0, 20, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    b"TCA".to_vec(),
                    AlleleSupportStats::new(15, -1.0, 15, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    b"GCA".to_vec(),
                    AlleleSupportStats::new(12, -1.0, 12, 0, 0, 0, 0),
                    Vec::new(),
                ),
            ],
        );
        let snp_100 = PileupRecord::new(
            0,
            100,
            vec![
                AlleleObservation::new(
                    b"A".to_vec(),
                    AlleleSupportStats::new(2, -1.0, 2, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    b"T".to_vec(),
                    AlleleSupportStats::new(1, -1.0, 1, 0, 0, 0, 0),
                    vec![1],
                ),
            ],
        );
        let snp_102 = PileupRecord::new(
            0,
            102,
            vec![
                AlleleObservation::new(
                    b"A".to_vec(),
                    AlleleSupportStats::new(2, -1.0, 2, 0, 0, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    b"G".to_vec(),
                    AlleleSupportStats::new(1, -1.0, 1, 0, 0, 0, 0),
                    vec![1],
                ),
            ],
        );
        let s0 = vec![snp_100, snp_102];
        let s1 = vec![mnp_rec];
        let mut ref_full = vec![b'A'; 200];
        ref_full[99] = b'A';
        ref_full[100] = b'C';
        ref_full[101] = b'A';
        let (chunk, partition, ref_seq) = group_fixture(vec![s0, s1], &ref_full, 1..200, 50);
        assert_eq!(ref_seq.as_slice(), b"ACA");

        // Byte-identity against the row-shape kernel — both should
        // produce the same compound-protected output.
        let mut scratch = UnifyAllelesScratch::new();
        let mut out = UnifiedAllelesColumns::empty();
        unify_alleles_columnar(
            &chunk,
            &partition,
            0,
            &ref_seq,
            chunk.n_samples(),
            3, // max_alleles — tight to force the cap
            &mut scratch,
            &mut out,
        )
        .expect("unify succeeded");
        let new_alleles = columnar_alleles_with_constituents(&out);

        let fetcher = shared_mock(&ref_full, 1);
        let oracle = {
            let group = crate::var_calling::cohort_block::worker::build_overlapping_variant_group(
                &chunk,
                &partition,
                0,
                chunk.n_samples(),
                chunk.chrom_id,
            );
            // Use max_alleles=3 in the oracle config too.
            let config = PerGroupMergerConfig::new(2, 3, 64, 32).expect("merger config");
            let iter: Vec<Result<OverlappingVariantGroup, GrouperError>> = vec![Ok(group)];
            let mut merger = PerGroupMerger::with_config(iter.into_iter(), fetcher, config);
            let record = merger
                .next()
                .expect("merger yields at least one item")
                .expect("merger succeeded");
            record
                .alleles
                .iter()
                .map(|a| {
                    (
                        a.seq.clone(),
                        a.is_compound,
                        a.constituents
                            .iter()
                            .map(|c| (c.record_idx, c.local_allele_idx))
                            .collect::<Vec<(usize, usize)>>(),
                    )
                })
                .collect::<Vec<_>>()
        };

        assert_eq!(new_alleles, oracle);
        // Confirm compound is present.
        assert!(out.is_compound.iter().any(|&x| x));
    }

    #[test]
    fn scratch_clear_drops_contents_preserves_capacity() {
        let mut scratch = UnifyAllelesScratch::new();
        scratch.byte_index.insert(b"AC".to_vec(), 1);
        scratch.projection_buf.extend_from_slice(b"ACGT");
        let proj_cap = scratch.projection_buf.capacity();
        scratch.clear();
        assert!(scratch.byte_index.is_empty());
        assert!(scratch.projection_buf.is_empty());
        assert!(scratch.projection_buf.capacity() >= proj_cap);
    }
}
