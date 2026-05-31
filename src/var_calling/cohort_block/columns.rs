//! Columnar in-memory storage for one chunk × N samples.
//!
//! [`SampleColumns`] is per-sample parallel-array + CSR storage of
//! [`PileupRecord`]-shaped data; [`MaterialisedChunk`] bundles N
//! such columns with the chunk's chromosome, target range, the
//! pre-pass-decided `safe_end`, and the partition into worker
//! windows. Materialising an owned [`PileupRecord`] for one
//! (sample, row) is on demand via [`SampleColumns::materialise_record`].
//!
//! M16: per-allele columns are grouped into three sibling sub-structs
//! ([`PerAlleleFixed`], [`PerAlleleSeq`], [`PerAlleleChainIds`]) so
//! that `clear` / `push_row_from` / `truncate` etc. become "delegate
//! to each group's same-named method" instead of one flat 11-column
//! enumeration per method. Adding a new per-allele column lands on
//! the matching sub-struct, not on every method of `SampleColumns`.

use std::ops::Range;

use crate::pileup_record::{AlleleObservation, AlleleSupportStats, ChainId, PileupRecord};

/// Per-(record, allele) fixed-width scalar columns — the 7 components
/// of [`AlleleSupportStats`]. One entry per allele across the whole
/// sample's chunk; total length equals the parent's
/// `allele_offsets[n_records]`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PerAlleleFixed {
    pub num_obs: Vec<u32>,
    pub q_sum: Vec<f64>,
    pub fwd: Vec<u32>,
    pub placed_left: Vec<u32>,
    pub placed_start: Vec<u32>,
    pub mapq_sum: Vec<u32>,
    pub mapq_sum_sq: Vec<u64>,
}

impl PerAlleleFixed {
    /// Total per-allele cells across all records.
    pub fn len(&self) -> usize {
        self.num_obs.len()
    }

    /// True when no per-allele cells are stored.
    pub fn is_empty(&self) -> bool {
        self.num_obs.is_empty()
    }

    /// Drop every cell while preserving every column's allocated
    /// capacity.
    pub fn clear(&mut self) {
        self.num_obs.clear();
        self.q_sum.clear();
        self.fwd.clear();
        self.placed_left.clear();
        self.placed_start.clear();
        self.mapq_sum.clear();
        self.mapq_sum_sq.clear();
    }

    /// Append the components of one allele's [`AlleleSupportStats`]
    /// — destructured up front so a future field on
    /// [`AlleleSupportStats`] is a compile error here (the columnar
    /// storage must carry every per-allele scalar to preserve byte-
    /// identity downstream).
    pub fn push(&mut self, support: AlleleSupportStats) {
        // M17: no trailing `..` — the 7 named fields cover every field
        // on `AlleleSupportStats` today, and a future field addition
        // should be a compile error here.
        let AlleleSupportStats {
            num_obs,
            q_sum,
            fwd,
            placed_left,
            placed_start,
            mapq_sum,
            mapq_sum_sq,
        } = support;
        self.num_obs.push(num_obs);
        self.q_sum.push(q_sum);
        self.fwd.push(fwd);
        self.placed_left.push(placed_left);
        self.placed_start.push(placed_start);
        self.mapq_sum.push(mapq_sum);
        self.mapq_sum_sq.push(mapq_sum_sq);
    }

    /// Append `src`'s cells in `range` to `self`.
    pub fn extend_from_range(&mut self, src: &Self, range: Range<usize>) {
        self.num_obs.extend_from_slice(&src.num_obs[range.clone()]);
        self.q_sum.extend_from_slice(&src.q_sum[range.clone()]);
        self.fwd.extend_from_slice(&src.fwd[range.clone()]);
        self.placed_left
            .extend_from_slice(&src.placed_left[range.clone()]);
        self.placed_start
            .extend_from_slice(&src.placed_start[range.clone()]);
        self.mapq_sum
            .extend_from_slice(&src.mapq_sum[range.clone()]);
        self.mapq_sum_sq.extend_from_slice(&src.mapq_sum_sq[range]);
    }

    /// Truncate every column to `len` cells.
    pub fn truncate(&mut self, len: usize) {
        self.num_obs.truncate(len);
        self.q_sum.truncate(len);
        self.fwd.truncate(len);
        self.placed_left.truncate(len);
        self.placed_start.truncate(len);
        self.mapq_sum.truncate(len);
        self.mapq_sum_sq.truncate(len);
    }

    /// Build an [`AlleleSupportStats`] from the cell at `allele_idx`.
    pub fn support_at(&self, allele_idx: usize) -> AlleleSupportStats {
        AlleleSupportStats::new(
            self.num_obs[allele_idx],
            self.q_sum[allele_idx],
            self.fwd[allele_idx],
            self.placed_left[allele_idx],
            self.placed_start[allele_idx],
            self.mapq_sum[allele_idx],
            self.mapq_sum_sq[allele_idx],
        )
    }
}

/// Per-allele variable-length sequence-bytes column, nested-CSR-laid-
/// out. `offsets[k]..offsets[k+1]` gives allele `k`'s slice of
/// [`Self::bytes`]; the sentinel `0` at the start keeps the invariant
/// `offsets.len() == n_alleles + 1`.
#[derive(Debug, Clone, PartialEq)]
pub struct PerAlleleSeq {
    pub offsets: Vec<u32>,
    pub bytes: Vec<u8>,
}

impl PerAlleleSeq {
    /// Empty columns ready to receive alleles. `offsets` starts with
    /// the CSR sentinel `0`.
    pub fn empty() -> Self {
        Self {
            offsets: vec![0],
            bytes: Vec::new(),
        }
    }

    /// Drop every cell while preserving allocated capacity. The CSR
    /// offset sentinel is restored so the post-clear state still
    /// honours the invariant.
    pub fn clear(&mut self) {
        self.offsets.clear();
        self.offsets.push(0);
        self.bytes.clear();
    }

    /// Append one allele's sequence bytes and push the resulting CSR
    /// offset.
    pub fn push(&mut self, seq: &[u8]) {
        self.bytes.extend_from_slice(seq);
        self.offsets.push(u32_from_usize(self.bytes.len()));
    }

    /// Append `src`'s alleles in `allele_range` to `self`, copying
    /// their bytes and rewriting offsets to be `self`-relative.
    pub fn extend_from_range(&mut self, src: &Self, allele_range: Range<usize>) {
        let seq_lo = src.offsets[allele_range.start] as usize;
        let seq_hi = src.offsets[allele_range.end] as usize;
        let bytes_base = self.bytes.len();
        self.bytes.extend_from_slice(&src.bytes[seq_lo..seq_hi]);
        for k in (allele_range.start + 1)..=allele_range.end {
            let inner_offset = src.offsets[k] as usize - seq_lo;
            self.offsets.push(u32_from_usize(bytes_base + inner_offset));
        }
    }

    /// Truncate to `n_alleles` alleles (drops their bytes too).
    pub fn truncate(&mut self, n_alleles: usize) {
        let bytes_keep = self.offsets[n_alleles] as usize;
        self.offsets.truncate(n_alleles + 1);
        self.bytes.truncate(bytes_keep);
    }

    /// Allele `allele_idx`'s sequence as a byte slice.
    pub fn slice_at(&self, allele_idx: usize) -> &[u8] {
        let lo = self.offsets[allele_idx] as usize;
        let hi = self.offsets[allele_idx + 1] as usize;
        &self.bytes[lo..hi]
    }
}

/// Per-allele variable-length chain-id column, nested-CSR-laid-out.
/// `offsets[k]..offsets[k+1]` gives allele `k`'s slice of
/// [`Self::ids`]; the sentinel `0` at the start keeps the invariant
/// `offsets.len() == n_alleles + 1`.
#[derive(Debug, Clone, PartialEq)]
pub struct PerAlleleChainIds {
    pub offsets: Vec<u32>,
    pub ids: Vec<ChainId>,
}

impl PerAlleleChainIds {
    /// Empty columns ready to receive alleles. `offsets` starts with
    /// the CSR sentinel `0`.
    pub fn empty() -> Self {
        Self {
            offsets: vec![0],
            ids: Vec::new(),
        }
    }

    /// Drop every cell while preserving allocated capacity. The CSR
    /// offset sentinel is restored so the post-clear state still
    /// honours the invariant.
    pub fn clear(&mut self) {
        self.offsets.clear();
        self.offsets.push(0);
        self.ids.clear();
    }

    /// Append one allele's chain ids and push the resulting CSR offset.
    pub fn push(&mut self, chain_ids: &[ChainId]) {
        self.ids.extend_from_slice(chain_ids);
        self.offsets.push(u32_from_usize(self.ids.len()));
    }

    /// Append `src`'s alleles in `allele_range` to `self`, copying
    /// their ids and rewriting offsets to be `self`-relative.
    pub fn extend_from_range(&mut self, src: &Self, allele_range: Range<usize>) {
        let id_lo = src.offsets[allele_range.start] as usize;
        let id_hi = src.offsets[allele_range.end] as usize;
        let ids_base = self.ids.len();
        self.ids.extend_from_slice(&src.ids[id_lo..id_hi]);
        for k in (allele_range.start + 1)..=allele_range.end {
            let inner_offset = src.offsets[k] as usize - id_lo;
            self.offsets.push(u32_from_usize(ids_base + inner_offset));
        }
    }

    /// Truncate to `n_alleles` alleles (drops their ids too).
    pub fn truncate(&mut self, n_alleles: usize) {
        let ids_keep = self.offsets[n_alleles] as usize;
        self.offsets.truncate(n_alleles + 1);
        self.ids.truncate(ids_keep);
    }

    /// Allele `allele_idx`'s chain ids as a slice.
    pub fn slice_at(&self, allele_idx: usize) -> &[ChainId] {
        let lo = self.offsets[allele_idx] as usize;
        let hi = self.offsets[allele_idx + 1] as usize;
        &self.ids[lo..hi]
    }
}

/// Per-sample, per-chunk columnar storage. Records are sorted by
/// 1-based genomic `position` ascending. Per-allele data lives in
/// three sub-structs that each hold their own CSR-laid-out flat
/// columns indexed via [`Self::allele_offsets`].
///
/// **Why columnar.** PSP's on-disk blocks are already column-oriented;
/// loading is a direct column-by-column copy with no intermediate
/// `PileupRecord` synthesis on the read path. The chunk worker walks
/// these columns directly; consumers that still want an owned
/// [`PileupRecord`] get one via [`Self::materialise_record`].
///
/// **Per-record fixed-width columns:**
/// - [`Self::positions`]: 1-based anchor position.
/// - [`Self::allele_offsets`]: CSR offsets into every per-allele flat
///   column; `[i, i+1)` gives the `i`-th record's allele range, and
///   the final entry equals the total allele count.
///
/// **Per-(record, allele) sub-structs** — each indexed by the same
/// CSR built on top of [`Self::allele_offsets`]:
/// - [`Self::per_allele_fixed`]: 7 scalar columns from
///   [`AlleleSupportStats`] (num_obs, q_sum, fwd, placed_left,
///   placed_start, mapq_sum, mapq_sum_sq).
/// - [`Self::per_allele_seq`]: allele-sequence bytes + nested CSR.
/// - [`Self::per_allele_chain_ids`]: phase-chain ids + nested CSR.
///
/// `chrom_id` is shared across the whole [`MaterialisedChunk`] and is
/// not duplicated here.
// Mi1: `#[non_exhaustive]` so external callers must use the
// `SampleColumns::empty()` / `Default::default()` factories (both
// produce the CSR-sentinel-seeded invariant-honoring state) rather
// than struct-literal construction. M1: hand-written `Default`
// delegating to `empty()` — the derived `Default` would have
// produced empty offset columns missing the `vec![0]` CSR sentinel
// and tripped `n_alleles_total()`'s `expect("CSR sentinel always
// present")`.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct SampleColumns {
    pub positions: Vec<u32>,
    pub allele_offsets: Vec<u32>,
    pub per_allele_fixed: PerAlleleFixed,
    pub per_allele_seq: PerAlleleSeq,
    pub per_allele_chain_ids: PerAlleleChainIds,
}

impl Default for SampleColumns {
    /// M1: delegate to [`Self::empty`] — the derived `Default` would
    /// produce empty offset columns missing the `vec![0]` CSR
    /// sentinel and break every consumer that calls
    /// [`Self::n_alleles_total`] or indexes
    /// `allele_offsets[record_idx + 1]`.
    fn default() -> Self {
        Self::empty()
    }
}

impl SampleColumns {
    /// Empty columns ready to receive records. The top-level CSR
    /// offset column and each per-allele sub-struct start with a
    /// single `0` sentinel so the invariant
    /// `len(offsets) == n + 1` holds at every step.
    pub fn empty() -> Self {
        Self {
            positions: Vec::new(),
            allele_offsets: vec![0],
            per_allele_fixed: PerAlleleFixed::default(),
            per_allele_seq: PerAlleleSeq::empty(),
            per_allele_chain_ids: PerAlleleChainIds::empty(),
        }
    }

    /// Number of records currently stored.
    pub fn n_records(&self) -> usize {
        self.positions.len()
    }

    /// Total per-allele cells across all records.
    pub fn n_alleles_total(&self) -> usize {
        *self
            .allele_offsets
            .last()
            .expect("CSR sentinel always present") as usize
    }

    /// Append one record's columnar data. Records must be pushed in
    /// strictly ascending position order; debug-asserted.
    ///
    /// `record` is consumed; its [`AlleleObservation`] entries are
    /// flattened into the per-(record, allele) sub-structs.
    pub fn push_record(&mut self, record: PileupRecord) {
        if let Some(&prev_pos) = self.positions.last() {
            debug_assert!(
                record.pos > prev_pos,
                "SampleColumns::push_record received non-monotonic position \
                 (previous {prev_pos}, new {})",
                record.pos,
            );
        }
        self.positions.push(record.pos);

        for allele in record.alleles {
            let AlleleObservation {
                seq,
                support,
                chain_ids,
            } = allele;
            self.per_allele_fixed.push(support);
            self.per_allele_seq.push(&seq);
            self.per_allele_chain_ids.push(&chain_ids);
        }
        self.allele_offsets
            .push(u32_from_usize(self.per_allele_fixed.len()));
    }

    /// Append records `[record_lo, record_hi)` of a decoded PSP block
    /// straight from its columns — the same shape [`push_record`] would
    /// build, but with no intermediate row-shape `PileupRecord` (and no
    /// per-allele `Vec`). This is the columnar→columnar fast path the
    /// span reader uses in place of the `columnar → row → columnar`
    /// round-trip.
    ///
    /// `allele_lo` is the cumulative allele offset of `record_lo` in
    /// `cols` (`Σ n_alleles[0..record_lo]`). `abs_positions` carries the
    /// already-decoded absolute 1-based positions of the appended
    /// records — the caller decodes them while it windows the block, so
    /// the delta-varint overflow check lives there rather than here;
    /// `abs_positions.len()` must equal `record_hi - record_lo`.
    ///
    /// The fixed-width per-allele columns are bulk-copied one slice each
    /// (they align 1:1 with the block's columns); the two ragged columns
    /// are pushed per allele, slicing straight from the block's CSR
    /// data; positions and the per-record CSR offsets are bulk / running
    /// appends. No allocation beyond the destination columns' growth.
    pub fn append_block_window(
        &mut self,
        cols: &crate::psp::BlockColumns<'_>,
        record_lo: usize,
        record_hi: usize,
        allele_lo: usize,
        abs_positions: &[u32],
    ) {
        debug_assert_eq!(abs_positions.len(), record_hi - record_lo);
        if record_lo == record_hi {
            return;
        }
        if let (Some(&prev), Some(&first)) = (self.positions.last(), abs_positions.first()) {
            debug_assert!(
                first > prev,
                "append_block_window: non-monotonic position (previous {prev}, new {first})",
            );
        }

        // Cumulative allele count of the appended record range, and the
        // allele base in *this* column set before the append.
        let n_alleles_appended: usize = cols.n_alleles[record_lo..record_hi]
            .iter()
            .map(|&n| n as usize)
            .sum();
        let allele_hi = allele_lo + n_alleles_appended;
        let base = self.per_allele_fixed.len();

        // Per-record positions — bulk copy of the caller's decode.
        self.positions.extend_from_slice(abs_positions);

        // Per-allele fixed-width columns — one bulk slice each.
        let f = &mut self.per_allele_fixed;
        f.num_obs
            .extend_from_slice(&cols.allele_obs_count[allele_lo..allele_hi]);
        f.q_sum
            .extend_from_slice(&cols.allele_q_sum_log[allele_lo..allele_hi]);
        f.fwd
            .extend_from_slice(&cols.allele_fwd_count[allele_lo..allele_hi]);
        f.placed_left
            .extend_from_slice(&cols.allele_placed_left_count[allele_lo..allele_hi]);
        f.placed_start
            .extend_from_slice(&cols.allele_placed_start_count[allele_lo..allele_hi]);
        f.mapq_sum
            .extend_from_slice(&cols.allele_mapq_sum[allele_lo..allele_hi]);
        f.mapq_sum_sq
            .extend_from_slice(&cols.allele_mapq_sum_sq[allele_lo..allele_hi]);

        // Ragged per-allele columns — CSR push straight from the block.
        for j in allele_lo..allele_hi {
            let seq_s = cols.allele_seq_offsets[j] as usize;
            let seq_e = cols.allele_seq_offsets[j + 1] as usize;
            self.per_allele_seq
                .push(&cols.allele_seq_data[seq_s..seq_e]);
            let cid_s = cols.allele_chain_ids_offsets[j] as usize;
            let cid_e = cols.allele_chain_ids_offsets[j + 1] as usize;
            self.per_allele_chain_ids
                .push(&cols.allele_chain_ids_data[cid_s..cid_e]);
        }

        // Per-record CSR offset (cumulative allele end, base-relative).
        let mut cum = 0usize;
        for &n in &cols.n_alleles[record_lo..record_hi] {
            cum += n as usize;
            self.allele_offsets.push(u32_from_usize(base + cum));
        }
    }

    /// 1-based position at row `record_idx`.
    pub fn position_at(&self, record_idx: usize) -> u32 {
        self.positions[record_idx]
    }

    /// Number of alleles at row `record_idx`.
    pub fn n_alleles_at(&self, record_idx: usize) -> usize {
        let lo = self.allele_offsets[record_idx] as usize;
        let hi = self.allele_offsets[record_idx + 1] as usize;
        hi - lo
    }

    /// Binary-search for a 1-based `position`. `Ok(i)` means the row
    /// at index `i` has exactly that position; `Err(i)` is the
    /// insertion point preserving sort order.
    pub fn binary_search_position(&self, position: u32) -> Result<usize, usize> {
        self.positions.binary_search(&position)
    }

    /// Materialise an owned [`PileupRecord`] for the row at
    /// `record_idx`. Per-allele data is copied out of the columnar
    /// storage — `Vec` allocations match what the PSP reader does on
    /// every record emit today, so this is no worse than the
    /// streaming path at the same boundary.
    pub fn materialise_record(&self, chrom_id: u32, record_idx: usize) -> PileupRecord {
        let pos = self.positions[record_idx];
        let allele_lo = self.allele_offsets[record_idx] as usize;
        let allele_hi = self.allele_offsets[record_idx + 1] as usize;
        let mut alleles = Vec::with_capacity(allele_hi - allele_lo);
        for k in allele_lo..allele_hi {
            alleles.push(self.materialise_allele(k));
        }
        PileupRecord::new(chrom_id, pos, alleles)
    }

    /// Drop every record but preserve every column's allocated
    /// capacity. CSR offset columns are reset to their single-`0`
    /// sentinel state. Use this between iterations of the chunk loop
    /// so the same buffer can absorb the next chunk without
    /// reallocation.
    pub fn clear(&mut self) {
        self.positions.clear();
        self.allele_offsets.clear();
        self.allele_offsets.push(0);
        self.per_allele_fixed.clear();
        self.per_allele_seq.clear();
        self.per_allele_chain_ids.clear();
    }

    /// Append the row at `src_row_idx` of `src` to `self` by copying
    /// directly between columns — no [`PileupRecord`] is synthesised
    /// in between. Position monotonicity is debug-asserted as in
    /// [`Self::push_record`].
    ///
    /// Used by the chunk loader's compact step to move filtered
    /// records out of a raw-load buffer into the post-filter chunk
    /// columns without materialising rows.
    pub fn push_row_from(&mut self, src: &SampleColumns, src_row_idx: usize) {
        let pos = src.positions[src_row_idx];
        if let Some(&prev_pos) = self.positions.last() {
            debug_assert!(
                pos > prev_pos,
                "SampleColumns::push_row_from received non-monotonic position \
                 (previous {prev_pos}, new {pos})",
            );
        }
        self.positions.push(pos);

        let allele_lo = src.allele_offsets[src_row_idx] as usize;
        let allele_hi = src.allele_offsets[src_row_idx + 1] as usize;
        let allele_range = allele_lo..allele_hi;

        self.per_allele_fixed
            .extend_from_range(&src.per_allele_fixed, allele_range.clone());
        self.per_allele_seq
            .extend_from_range(&src.per_allele_seq, allele_range.clone());
        self.per_allele_chain_ids
            .extend_from_range(&src.per_allele_chain_ids, allele_range);

        self.allele_offsets
            .push(u32_from_usize(self.per_allele_fixed.len()));
    }

    /// M14: replace this column's contents with a copy of `other`'s,
    /// preserving `self`'s allocated capacity.
    ///
    /// Used by the chunk driver's `NoSafeGap` retry path to snapshot
    /// and restore the carryover the previous chunk handed in (the
    /// loader drains it as part of the raw-load step, so a retry has
    /// to put it back). The single-method form replaces two
    /// hand-written `clear()` + `push_row_from` loops at the call
    /// sites — adding a new column field to `SampleColumns` now
    /// updates one place instead of two parallel loops.
    pub fn clone_from_columns(&mut self, other: &SampleColumns) {
        self.clear();
        for row_idx in 0..other.n_records() {
            self.push_row_from(other, row_idx);
        }
    }

    /// Reference span of the record at `record_idx` — the number of
    /// reference bases the record's REF allele covers. By the
    /// walker invariant `alleles[0]` is REF, so this is the byte
    /// length of allele-index 0's `seq` (recoverable from the
    /// nested-CSR offsets without materialising the bytes).
    pub fn ref_span_at(&self, record_idx: usize) -> u32 {
        let allele_lo = self.allele_offsets[record_idx] as usize;
        let seq_lo = self.per_allele_seq.offsets[allele_lo];
        let seq_hi = self.per_allele_seq.offsets[allele_lo + 1];
        seq_hi - seq_lo
    }

    /// Truncate `self` to its first `keep_n_records` rows, dropping
    /// every CSR-indexed allele payload past the cutoff. Capacities
    /// are preserved — the CSR offset sentinels stay in place and
    /// the underlying allocations can absorb future pushes.
    pub fn truncate(&mut self, keep_n_records: usize) {
        debug_assert!(
            keep_n_records <= self.n_records(),
            "SampleColumns::truncate cannot grow ({keep_n_records} > {})",
            self.n_records(),
        );
        let allele_keep = self.allele_offsets[keep_n_records] as usize;

        self.positions.truncate(keep_n_records);
        self.allele_offsets.truncate(keep_n_records + 1);
        self.per_allele_fixed.truncate(allele_keep);
        self.per_allele_seq.truncate(allele_keep);
        self.per_allele_chain_ids.truncate(allele_keep);
    }

    /// Move rows `split_row_idx..n_records()` into `dst` (appending
    /// in original order) and truncate `self` to `split_row_idx`
    /// rows. Both endpoints preserve allocated capacity.
    pub fn drain_rows_from_into(&mut self, split_row_idx: usize, dst: &mut SampleColumns) {
        for row_idx in split_row_idx..self.n_records() {
            dst.push_row_from(self, row_idx);
        }
        self.truncate(split_row_idx);
    }

    /// True when at least one allele at `record_idx` carries a
    /// non-reference allele with `num_obs > 0`. The walker invariant
    /// `alleles[0] == REF` means alleles at indices `>= 1` are
    /// non-REF by construction; "any with `num_obs > 0`" is the
    /// per-sample-per-position component of the cohort-wide
    /// variant-position filter.
    pub fn has_observed_non_ref_allele_at(&self, record_idx: usize) -> bool {
        let allele_lo = self.allele_offsets[record_idx] as usize;
        let allele_hi = self.allele_offsets[record_idx + 1] as usize;
        // `alleles[0]` is REF; skip it. Walker invariant from
        // `pileup_record.rs`.
        self.per_allele_fixed.num_obs[(allele_lo + 1)..allele_hi]
            .iter()
            .any(|&n| n > 0)
    }

    fn materialise_allele(&self, allele_idx: usize) -> AlleleObservation {
        let support = self.per_allele_fixed.support_at(allele_idx);
        AlleleObservation::new(
            self.per_allele_seq.slice_at(allele_idx).to_vec(),
            support,
            self.per_allele_chain_ids.slice_at(allele_idx).to_vec(),
        )
    }
}

/// One materialised block × N samples worth of [`SampleColumns`] plus
/// the pre-pass-decided right boundary.
///
/// **Lifecycle.**
/// 1. The chunk loader produces a block with `safe_end = range.end`,
///    having applied the cohort-wide variant-position filter.
/// 2. `finalise_chunk_boundaries` walks the loaded block once and picks
///    `safe_end ≤ range.end` at a clean group boundary so no variant
///    group can span the right boundary into the next block; records
///    past `safe_end` are split into the carryover (the reserve).
/// 3. The consumer runs the full per-group pipeline on the whole block
///    `[range.start, safe_end)` and emits final `PosteriorRecord`s.
// Mi1: `#[non_exhaustive]` — external callers construct via
// `MaterialisedChunk::with_n_samples(n)` (or the test-only struct
// literal). Future column additions land without breaking out-of-crate
// struct-literal sites.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialisedChunk {
    pub chrom_id: u32,
    /// `[G_start, G_end)` — the initial target range. May be revised
    /// down to a `safe_end ≤ G_end` by `finalise_chunk_boundaries`.
    pub range: Range<u32>,
    /// Right boundary chosen by the pre-pass. Records at positions
    /// `>= safe_end` are split into the carryover passed to the next
    /// chunk's load.
    pub safe_end: u32,
    /// One [`SampleColumns`] per sample, sample-index aligned with
    /// the cohort's `psp_paths`.
    pub per_sample: Vec<SampleColumns>,
}

impl MaterialisedChunk {
    /// Materialise the row at `record_idx` of `sample_idx`, mounting
    /// the chunk's `chrom_id`.
    pub fn materialise_record(&self, sample_idx: usize, record_idx: usize) -> PileupRecord {
        self.per_sample[sample_idx].materialise_record(self.chrom_id, record_idx)
    }

    /// Number of samples represented by this chunk.
    pub fn n_samples(&self) -> usize {
        self.per_sample.len()
    }

    /// Construct an empty chunk pre-sized for `n_samples`. The
    /// per-sample columns are empty and `range` is `0..0` /
    /// `safe_end` is `0` — the loader and pre-pass will populate them.
    pub fn with_n_samples(n_samples: usize) -> Self {
        Self {
            chrom_id: 0,
            range: 0..0,
            safe_end: 0,
            per_sample: (0..n_samples).map(|_| SampleColumns::empty()).collect(),
        }
    }

    /// Reset every record / range field but preserve all
    /// allocated column capacity. The chunk loader calls this before
    /// loading the next chunk; the driver holds one persistent
    /// `MaterialisedChunk` across iterations.
    ///
    /// Mi2: named `clear()` to match the module-wide convention
    /// (`SampleColumns::clear`, `WindowPartition::clear`,
    /// `ChunkLoadScratch::clear`, etc.). Previously `clear_data` —
    /// the `_data` suffix carried no meaning the type didn't already
    /// imply.
    pub fn clear(&mut self) {
        self.chrom_id = 0;
        self.range = 0..0;
        self.safe_end = 0;
        for sample in &mut self.per_sample {
            sample.clear();
        }
    }
}

/// Convert a `usize` length into a `u32` CSR offset.
///
/// M6: the PSP format caps per-chunk record counts well below
/// `u32::MAX` for any realistic input — but PSP is a disk-derived
/// boundary, so a corrupted or adversarial input could in principle
/// push the offsets past `u32::MAX` and silently corrupt every
/// downstream column index. The previous "release builds wrap
/// silently" contract turned that into wrong likelihoods on the
/// affected chunk. The new contract: panic loudly with a named
/// invariant in both debug and release builds. The cost on the hot
/// path is one branch — well below the cost of a single CSR push.
#[inline]
pub(crate) fn u32_from_usize(value: usize) -> u32 {
    // PANIC-FREE on realistic inputs: PSP per-chunk record counts are
    // bounded by `target_variants_per_chunk × n_samples` (operator
    // input, capped at u32::MAX by the field type itself). A panic
    // here surfaces a corrupted PSP or an upstream invariant break;
    // either way it is a hard error, not a silent miscompute.
    u32::try_from(value).expect("CSR offset exceeds u32::MAX")
}

// Nits (Wave 6): `vec![a..b]` is the natural test-fixture form for
// `MaterialisedChunk.windows` — clippy's `single_range_in_vec_init`
// flags it as ambiguous between "a Vec of one Range" and "a Vec of
// length b filled with value a". Allow at the module level.
#[allow(clippy::single_range_in_vec_init)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::cohort_block::test_helpers::{allele, record, ref_plus_alt};

    #[test]
    fn empty_columns_carry_csr_sentinels() {
        let columns = SampleColumns::empty();
        assert_eq!(columns.n_records(), 0);
        assert_eq!(columns.n_alleles_total(), 0);
        assert_eq!(columns.allele_offsets, vec![0]);
        assert_eq!(columns.per_allele_seq.offsets, vec![0]);
        assert_eq!(columns.per_allele_chain_ids.offsets, vec![0]);
    }

    #[test]
    fn push_record_round_trips_via_materialise() {
        let mut columns = SampleColumns::empty();
        let r1 = record(
            10,
            vec![
                allele(b"A", 12, -1.5, &[1, 7]),
                allele(b"AT", 3, -0.5, &[42]),
            ],
        );
        let r2 = record(15, vec![allele(b"G", 8, -1.1, &[])]);

        columns.push_record(r1.clone());
        columns.push_record(r2.clone());

        assert_eq!(columns.n_records(), 2);
        assert_eq!(columns.n_alleles_total(), 3);
        assert_eq!(columns.n_alleles_at(0), 2);
        assert_eq!(columns.n_alleles_at(1), 1);

        assert_eq!(columns.materialise_record(0, 0), r1);
        assert_eq!(columns.materialise_record(0, 1), r2);
    }

    #[test]
    fn binary_search_position_finds_present_and_absent_rows() {
        let mut columns = SampleColumns::empty();
        for pos in [10_u32, 15, 22, 23, 40] {
            columns.push_record(record(pos, vec![allele(b"A", 1, -1.0, &[])]));
        }
        assert_eq!(columns.binary_search_position(10), Ok(0));
        assert_eq!(columns.binary_search_position(23), Ok(3));
        assert_eq!(columns.binary_search_position(11), Err(1));
        assert_eq!(columns.binary_search_position(50), Err(5));
    }

    #[test]
    #[should_panic(expected = "non-monotonic position")]
    fn push_record_rejects_non_monotonic_positions_in_debug() {
        let mut columns = SampleColumns::empty();
        columns.push_record(record(10, vec![allele(b"A", 1, -1.0, &[])]));
        columns.push_record(record(10, vec![allele(b"A", 1, -1.0, &[])]));
    }

    #[test]
    fn materialised_chunk_dispatches_record_to_per_sample_columns() {
        let mut s0 = SampleColumns::empty();
        s0.push_record(record(100, vec![allele(b"A", 5, -2.0, &[1])]));
        let mut s1 = SampleColumns::empty();
        s1.push_record(record(100, vec![allele(b"T", 7, -1.0, &[2])]));

        let chunk = MaterialisedChunk {
            chrom_id: 3,
            range: 90..200,
            safe_end: 200,
            per_sample: vec![s0, s1],
        };
        assert_eq!(chunk.n_samples(), 2);
        let r0 = chunk.materialise_record(0, 0);
        let r1 = chunk.materialise_record(1, 0);
        assert_eq!(r0.chrom_id, 3);
        assert_eq!(r0.alleles[0].seq, b"A".to_vec());
        assert_eq!(r1.alleles[0].seq, b"T".to_vec());
    }

    #[test]
    fn clear_preserves_csr_sentinels_and_capacity() {
        let mut columns = SampleColumns::empty();
        for pos in [10_u32, 12, 14] {
            columns.push_record(record(pos, vec![allele(b"A", 5, -1.0, &[1])]));
        }
        let positions_cap_before = columns.positions.capacity();
        columns.clear();
        assert_eq!(columns.n_records(), 0);
        assert_eq!(columns.n_alleles_total(), 0);
        assert_eq!(columns.allele_offsets, vec![0]);
        assert_eq!(columns.per_allele_seq.offsets, vec![0]);
        assert_eq!(columns.per_allele_chain_ids.offsets, vec![0]);
        assert!(columns.positions.capacity() >= positions_cap_before);
    }

    #[test]
    fn push_row_from_preserves_record_payload() {
        let mut src = SampleColumns::empty();
        let r1 = record(
            10,
            vec![
                allele(b"A", 12, -1.5, &[1, 7]),
                allele(b"AT", 3, -0.5, &[42]),
            ],
        );
        let r2 = record(15, vec![allele(b"G", 8, -1.1, &[3, 9, 11])]);
        src.push_record(r1.clone());
        src.push_record(r2.clone());

        let mut dst = SampleColumns::empty();
        dst.push_row_from(&src, 0);
        dst.push_row_from(&src, 1);
        assert_eq!(dst.n_records(), 2);
        assert_eq!(dst.materialise_record(0, 0), r1);
        assert_eq!(dst.materialise_record(0, 1), r2);
    }

    #[test]
    fn truncate_drops_rows_and_keeps_csr_consistent() {
        let mut columns = SampleColumns::empty();
        for pos in [10_u32, 12, 14, 16] {
            columns.push_record(record(pos, ref_plus_alt(3, 4)));
        }
        assert_eq!(columns.n_records(), 4);
        columns.truncate(2);
        assert_eq!(columns.n_records(), 2);
        assert_eq!(columns.n_alleles_total(), 4); // 2 alleles × 2 rows
        assert_eq!(columns.materialise_record(0, 0).pos, 10);
        assert_eq!(columns.materialise_record(0, 1).pos, 12);
        // Push more — the truncated buffer must accept ascending
        // positions cleanly.
        columns.push_record(record(20, ref_plus_alt(2, 3)));
        assert_eq!(columns.n_records(), 3);
        assert_eq!(columns.materialise_record(0, 2).pos, 20);
    }

    #[test]
    fn drain_rows_from_into_moves_tail_preserving_round_trip() {
        let mut src = SampleColumns::empty();
        for pos in [10_u32, 12, 14, 16] {
            src.push_record(record(pos, ref_plus_alt(3, 4)));
        }
        let mut dst = SampleColumns::empty();
        src.drain_rows_from_into(2, &mut dst);
        assert_eq!(src.n_records(), 2);
        assert_eq!(dst.n_records(), 2);
        assert_eq!(dst.materialise_record(0, 0).pos, 14);
        assert_eq!(dst.materialise_record(0, 1).pos, 16);
    }

    /// Mi22: `split_row_idx == 0` boundary — drain every row into the
    /// destination, leaving the source empty.
    #[test]
    fn drain_rows_from_into_split_at_zero_moves_everything() {
        let mut src = SampleColumns::empty();
        for pos in [10_u32, 12, 14, 16] {
            src.push_record(record(pos, ref_plus_alt(3, 4)));
        }
        let mut dst = SampleColumns::empty();
        src.drain_rows_from_into(0, &mut dst);
        assert_eq!(src.n_records(), 0);
        assert_eq!(dst.n_records(), 4);
        assert_eq!(dst.materialise_record(0, 0).pos, 10);
        assert_eq!(dst.materialise_record(0, 3).pos, 16);
    }

    /// Mi22: `split_row_idx == n_records()` boundary — drain nothing.
    /// The destination stays empty and the source is unchanged.
    #[test]
    fn drain_rows_from_into_split_at_end_moves_nothing() {
        let mut src = SampleColumns::empty();
        for pos in [10_u32, 12, 14, 16] {
            src.push_record(record(pos, ref_plus_alt(3, 4)));
        }
        let mut dst = SampleColumns::empty();
        src.drain_rows_from_into(src.n_records(), &mut dst);
        assert_eq!(src.n_records(), 4);
        assert_eq!(dst.n_records(), 0);
        assert_eq!(src.materialise_record(0, 0).pos, 10);
        assert_eq!(src.materialise_record(0, 3).pos, 16);
    }
}
