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
//! **This sub-step (1.0) is types-only.** The algorithm lands in
//! sub-step 1.1 (per-position projection) and is extended by
//! sub-step 1.2 (compound detection) and sub-step 1.3 (max-alleles
//! cap). Until layer 1.4 wires the kernel into the worker, these
//! types are reachable only through their own unit tests; the
//! Phase A.0 row-shape adapter remains the production path.

use ahash::AHashMap;

use crate::pileup_record::ChainId;

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

/// Reusable scratch buffers for [`unify_alleles_columnar`] (lands
/// in sub-step 1.1). Cleared per group; never reallocated.
#[derive(Debug, Default)]
pub struct UnifyAllelesScratch {
    /// Per-position projection's byte-sequence → allele-index dedup
    /// map. Persists across groups; cleared at entry to each
    /// `unify_alleles_columnar` call.
    pub(crate) byte_index: AHashMap<Vec<u8>, usize>,
    /// Working buffer for projecting a candidate compound's bytes
    /// onto the group's REF span before testing for byte-equality
    /// against existing entries.
    pub(crate) compound_projected: Vec<u8>,
    /// Working buffer for chain-id intersection across constituent
    /// records during compound detection.
    pub(crate) chain_id_scratch: Vec<ChainId>,
}

impl UnifyAllelesScratch {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn clear(&mut self) {
        self.byte_index.clear();
        self.compound_projected.clear();
        self.chain_id_scratch.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn scratch_clear_drops_contents_preserves_capacity() {
        let mut scratch = UnifyAllelesScratch::new();
        scratch.byte_index.insert(b"AC".to_vec(), 1);
        scratch.compound_projected.extend_from_slice(b"ACGT");
        scratch.chain_id_scratch.push(7);
        let proj_cap = scratch.compound_projected.capacity();
        scratch.clear();
        assert!(scratch.byte_index.is_empty());
        assert!(scratch.compound_projected.is_empty());
        assert!(scratch.chain_id_scratch.is_empty());
        assert!(scratch.compound_projected.capacity() >= proj_cap);
    }
}
