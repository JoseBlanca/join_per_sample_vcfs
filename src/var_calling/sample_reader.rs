//! Per-sample reading (appendix §A) — `SamplePspReader` + `SamplePspChunk`.
//!
//! All per-sample responsibility lives here. **No dust** (the producer is the
//! sole dust consumer, §2.3). One reader is owned per sample; the producer
//! decodes the readers in parallel via disjoint `&mut` (`par_iter_mut`), so
//! `R: Send` is required (appendix §A).
//!
//! ## Two decode paths
//!
//! - **Production — column-selective two-phase.**
//!   [`next_two_phase`](SamplePspReader::next_two_phase) yields one
//!   [`TwoPhaseSegment`] per psp **segment** with the *light* fold columns
//!   ([`positions`](TwoPhaseSegment::positions) /
//!   [`nonref_obs`](TwoPhaseSegment::nonref_obs) /
//!   [`ref_spans`](TwoPhaseSegment::ref_spans)) decoded eagerly and the heavy
//!   columns left deferred-compressed. Once the cohort fold decides which rows
//!   are variable, [`set_variable_rows`](TwoPhaseSegment::set_variable_rows)
//!   inflates only the heavy columns of the *kept* rows into a
//!   [`SamplePspChunk`]. This is the path the producer uses.
//! - **Oracle — eager whole-segment decode.**
//!   [`SamplePspChunk::from_block`] decodes every column of a segment up front
//!   and the typed getters ([`take_seq`](SamplePspChunk::take_seq) /
//!   [`take_chain_ids`](SamplePspChunk::take_chain_ids) /
//!   [`take_scalar`](SamplePspChunk::take_scalar)) move out a `keep` subset.
//!   This is **only** used by the tests, as the byte-identity oracle the
//!   two-phase path is checked against.

use std::io::{Read, Seek};

use crate::pileup_record::{AlleleObservation, AlleleSupportStats, ChainId, PileupRecord};
#[cfg(test)]
use crate::psp::BlockColumns;
use crate::psp::ScalarDecodeError;
use crate::psp::reader::{DecodedColumn, RetainedColumn, TwoPhaseBlock, inflate_retained_column};
use crate::psp::{BlockColumnReader, BlockIndexEntry, PspReadError, PspReader};

// ---------------------------------------------------------------------------
// Heavy per-allele column carriers.
//
// These are the byte-identity-sensitive destination shapes the appendix §A
// take-getters return — the minimal columnar subset the reader needs (the
// full columnar EM machinery is not part of this pipeline).
// ---------------------------------------------------------------------------

/// Per-(record, allele) fixed-width scalar columns — the 7 components of
/// [`AlleleSupportStats`], one cell per allele.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct AlleleScalarColumns {
    pub num_obs: Vec<u32>,
    pub q_sum: Vec<f64>,
    pub fwd: Vec<u32>,
    pub placed_left: Vec<u32>,
    pub placed_start: Vec<u32>,
    pub mapq_sum: Vec<u32>,
    pub mapq_sum_sq: Vec<u64>,
}

impl AlleleScalarColumns {
    /// Total per-allele cells across all records. Test-only (the eager-decode
    /// oracle path); production reads cells via [`support_at`](Self::support_at).
    #[cfg(test)]
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.num_obs.len()
    }

    /// Append `src`'s cells in `range` to `self`.
    pub fn extend_from_range(&mut self, src: &Self, range: std::ops::Range<usize>) {
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

/// Per-allele variable-length sequence-bytes column, nested-CSR.
/// `offsets[k]..offsets[k + 1]` is allele `k`'s slice of [`Self::bytes`];
/// the leading `0` keeps `offsets.len() == n_alleles + 1`.
#[derive(Debug, Clone, PartialEq)]
pub struct AlleleSeqColumns {
    pub offsets: Vec<u32>,
    pub bytes: Vec<u8>,
}

impl AlleleSeqColumns {
    /// Empty columns, CSR sentinel seeded.
    pub fn empty() -> Self {
        Self {
            offsets: vec![0],
            bytes: Vec::new(),
        }
    }

    /// Append one allele's sequence bytes and push the CSR offset. Test-only
    /// (the eager-decode oracle builds columns allele-by-allele); production
    /// fills columns via [`extend_from_range`](Self::extend_from_range).
    #[cfg(test)]
    pub fn push(&mut self, seq: &[u8]) {
        self.bytes.extend_from_slice(seq);
        self.offsets.push(u32_from_usize(self.bytes.len()));
    }

    /// Append `src`'s alleles in `allele_range` to `self`, copying bytes and
    /// rewriting offsets `self`-relative.
    pub fn extend_from_range(&mut self, src: &Self, allele_range: std::ops::Range<usize>) {
        let seq_lo = src.offsets[allele_range.start] as usize;
        let seq_hi = src.offsets[allele_range.end] as usize;
        let bytes_base = self.bytes.len();
        self.bytes.extend_from_slice(&src.bytes[seq_lo..seq_hi]);
        for k in (allele_range.start + 1)..=allele_range.end {
            let inner_offset = src.offsets[k] as usize - seq_lo;
            self.offsets.push(u32_from_usize(bytes_base + inner_offset));
        }
    }

    /// Number of alleles stored. Test-only (the eager-decode oracle).
    #[cfg(test)]
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Allele `allele_idx`'s sequence bytes.
    pub fn slice_at(&self, allele_idx: usize) -> &[u8] {
        let lo = self.offsets[allele_idx] as usize;
        let hi = self.offsets[allele_idx + 1] as usize;
        &self.bytes[lo..hi]
    }
}

/// Per-allele variable-length chain-id column, nested-CSR (mirrors
/// [`AlleleSeqColumns`] for [`ChainId`]s).
#[derive(Debug, Clone, PartialEq)]
pub struct AlleleChainIdColumns {
    pub offsets: Vec<u32>,
    pub ids: Vec<ChainId>,
}

impl AlleleChainIdColumns {
    /// Empty columns, CSR sentinel seeded.
    pub fn empty() -> Self {
        Self {
            offsets: vec![0],
            ids: Vec::new(),
        }
    }

    /// Append one allele's chain ids and push the CSR offset. Test-only (the
    /// eager-decode oracle); production fills via
    /// [`extend_from_range`](Self::extend_from_range).
    #[cfg(test)]
    pub fn push(&mut self, chain_ids: &[ChainId]) {
        self.ids.extend_from_slice(chain_ids);
        self.offsets.push(u32_from_usize(self.ids.len()));
    }

    /// Append `src`'s alleles in `allele_range` to `self`.
    pub fn extend_from_range(&mut self, src: &Self, allele_range: std::ops::Range<usize>) {
        let id_lo = src.offsets[allele_range.start] as usize;
        let id_hi = src.offsets[allele_range.end] as usize;
        let ids_base = self.ids.len();
        self.ids.extend_from_slice(&src.ids[id_lo..id_hi]);
        for k in (allele_range.start + 1)..=allele_range.end {
            let inner_offset = src.offsets[k] as usize - id_lo;
            self.offsets.push(u32_from_usize(ids_base + inner_offset));
        }
    }

    /// Number of alleles stored. Test-only (the eager-decode oracle).
    #[cfg(test)]
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Allele `allele_idx`'s chain ids.
    pub fn slice_at(&self, allele_idx: usize) -> &[ChainId] {
        let lo = self.offsets[allele_idx] as usize;
        let hi = self.offsets[allele_idx + 1] as usize;
        &self.ids[lo..hi]
    }
}

fn u32_from_usize(v: usize) -> u32 {
    debug_assert!(v <= u32::MAX as usize, "column offset {v} exceeds u32");
    v as u32
}

// ---------------------------------------------------------------------------
// SamplePspChunk
// ---------------------------------------------------------------------------

/// One sample's columns for **one psp segment**, clamped to the reader's
/// region (appendix §A `[NEW]`).
///
/// `new` decodes and caches the light fold columns
/// ([`positions`](Self::positions) / [`nonref_obs`](Self::nonref_obs) /
/// [`ref_spans`](Self::ref_spans), one value per record); the heavy
/// per-allele columns are owned and handed out by the move-out take-getters
/// for the `keep` rows. Record indices used by the getters and by `keep` are
/// `0..len()` over this chunk's (region-clamped) records.
#[derive(Debug, Clone)]
pub struct SamplePspChunk {
    chrom_id: u32,
    /// 1-based positions, one per record. Light, cached at `new`.
    positions: Vec<u32>,
    /// Per-record sum of non-REF allele obs (alleles `[1..]`). Light.
    nonref_obs: Vec<u32>,
    /// Per-record REF-allele (`alleles[0]`) reference span in bases. Light.
    ref_spans: Vec<u32>,
    /// CSR record→allele offsets (`len() == n_records + 1`); maps a record
    /// index to its allele range in the heavy columns.
    allele_offsets: Vec<u32>,
    /// Per-record centred-window GC fraction (psp column `WindowedGc`), one
    /// value per record, aligned with `positions`. `NaN` where the pileup had
    /// no window centred here (e.g. an N reference base). Empty for a chunk
    /// built from a legacy `.psp` that predates the windowed columns. Read by
    /// the hidden-paralog score; carried through every per-record maintenance
    /// site alongside `positions`.
    windowed_gc: Vec<f32>,
    /// Per-record centred-window mean coverage (psp column `WindowedCoverage`),
    /// mirror of [`windowed_gc`](Self::windowed_gc).
    windowed_coverage: Vec<f32>,
    /// Heavy columns, owned; moved out (and emptied) by the take-getters.
    scalar: AlleleScalarColumns,
    seq: AlleleSeqColumns,
    chain_ids: AlleleChainIdColumns,
}

/// Byte-identity equality (the `.psp` round-trip oracle compares chunks). The
/// windowed columns hold `NaN` for un-windowed records, so they are compared on
/// their bit pattern — a derived `PartialEq` would make `NaN != NaN` and break
/// reflexivity. Must not derive `Eq`/`Hash` for the same reason (see
/// `PileupRecord`). Exhaustive destructure so a new field cannot silently skip
/// the comparison.
impl PartialEq for SamplePspChunk {
    fn eq(&self, other: &Self) -> bool {
        let Self {
            chrom_id,
            positions,
            nonref_obs,
            ref_spans,
            allele_offsets,
            windowed_gc,
            windowed_coverage,
            scalar,
            seq,
            chain_ids,
        } = self;
        *chrom_id == other.chrom_id
            && *positions == other.positions
            && *nonref_obs == other.nonref_obs
            && *ref_spans == other.ref_spans
            && *allele_offsets == other.allele_offsets
            && windowed_gc.len() == other.windowed_gc.len()
            && windowed_gc
                .iter()
                .zip(&other.windowed_gc)
                .all(|(a, b)| a.to_bits() == b.to_bits())
            && windowed_coverage.len() == other.windowed_coverage.len()
            && windowed_coverage
                .iter()
                .zip(&other.windowed_coverage)
                .all(|(a, b)| a.to_bits() == b.to_bits())
            && *scalar == other.scalar
            && *seq == other.seq
            && *chain_ids == other.chain_ids
    }
}

/// Chain ids to attach to allele index `a` of a record whose allele range
/// is `lo..hi`. The REF allele is always the record's first allele
/// (`a == lo`), and its chain ids are never read downstream — Stage 5's
/// compound-haplotype check (`per_group_merger::build_chain_proposals`)
/// links non-reference alleles only and skips allele 0. So REF gets an
/// empty `Vec` instead of a clone of the column. New `.psp` files no
/// longer store REF chain ids at all (the walker drops them); this guard
/// also covers any legacy `.psp` that still carries them, so they never
/// reach the merger regardless of file age. See phase_chain.md §8.
#[inline]
fn chain_ids_for_allele(chain_ids: &AlleleChainIdColumns, a: usize, lo: usize) -> Vec<ChainId> {
    if a == lo {
        Vec::new()
    } else {
        chain_ids.slice_at(a).to_vec()
    }
}

impl SamplePspChunk {
    /// Decode one block's region-clamped records into an owned chunk — the
    /// eager whole-block decode used as the **test-only byte-identity oracle**
    /// (production decodes column-selectively via [`TwoPhaseSegment`]).
    ///
    /// Records with `pos < floor` or `pos > region_end` are dropped (only
    /// the region-edge blocks carry any). Returns `None` when no record of
    /// the block falls in `[floor, region_end]`.
    #[cfg(test)]
    fn from_block(
        cols: &BlockColumns<'_>,
        chrom_id: u32,
        floor: u32,
        region_end: u32,
    ) -> Result<Option<Self>, PspReadError> {
        let n = cols.n_records as usize;

        // Absolute positions for the whole block (delta-decoded).
        let mut abs = Vec::with_capacity(n);
        let mut p = cols.first_pos;
        for i in 0..n {
            if i > 0 {
                p = advance_pos(p, cols.delta_pos[i], i)?;
            }
            abs.push(p);
        }

        // Region clamp: ascending positions ⇒ a contiguous record range.
        let r_lo = abs.partition_point(|&pos| pos < floor);
        let r_hi = abs.partition_point(|&pos| pos <= region_end);
        if r_lo >= r_hi {
            return Ok(None);
        }

        // Allele range of the clamped record range [r_lo, r_hi): `a_lo` is
        // the cumulative allele count before record `r_lo`.
        let a_lo: usize = cols.n_alleles[..r_lo].iter().map(|&n| n as usize).sum();
        let a_hi = a_lo
            + cols.n_alleles[r_lo..r_hi]
                .iter()
                .map(|&n| n as usize)
                .sum::<usize>();

        // CSR record→allele offsets, chunk-relative.
        let mut allele_offsets = Vec::with_capacity(r_hi - r_lo + 1);
        allele_offsets.push(0u32);
        let mut cum = 0u32;
        for &na in &cols.n_alleles[r_lo..r_hi] {
            cum += na as u32;
            allele_offsets.push(cum);
        }

        // Heavy columns for the clamped allele range.
        let mut scalar = AlleleScalarColumns::default();
        scalar
            .num_obs
            .extend_from_slice(&cols.allele_obs_count[a_lo..a_hi]);
        scalar
            .q_sum
            .extend_from_slice(&cols.allele_q_sum_log[a_lo..a_hi]);
        scalar
            .fwd
            .extend_from_slice(&cols.allele_fwd_count[a_lo..a_hi]);
        scalar
            .placed_left
            .extend_from_slice(&cols.allele_placed_left_count[a_lo..a_hi]);
        scalar
            .placed_start
            .extend_from_slice(&cols.allele_placed_start_count[a_lo..a_hi]);
        scalar
            .mapq_sum
            .extend_from_slice(&cols.allele_mapq_sum[a_lo..a_hi]);
        scalar
            .mapq_sum_sq
            .extend_from_slice(&cols.allele_mapq_sum_sq[a_lo..a_hi]);

        let mut seq = AlleleSeqColumns::empty();
        let mut chain_ids = AlleleChainIdColumns::empty();
        for k in a_lo..a_hi {
            let s = cols.allele_seq_offsets[k] as usize;
            let e = cols.allele_seq_offsets[k + 1] as usize;
            seq.push(&cols.allele_seq_data[s..e]);
            let cs = cols.allele_chain_ids_offsets[k] as usize;
            let ce = cols.allele_chain_ids_offsets[k + 1] as usize;
            chain_ids.push(&cols.allele_chain_ids_data[cs..ce]);
        }

        // Light per-record columns over the clamped range.
        let positions = abs[r_lo..r_hi].to_vec();
        let windowed_gc = cols.windowed_gc[r_lo..r_hi].to_vec();
        let windowed_coverage = cols.windowed_coverage[r_lo..r_hi].to_vec();
        let mut nonref_obs = Vec::with_capacity(r_hi - r_lo);
        let mut ref_spans = Vec::with_capacity(r_hi - r_lo);
        let mut base = a_lo;
        for r in r_lo..r_hi {
            let lo = base;
            let hi = lo + cols.n_alleles[r] as usize;
            base = hi;
            // Non-REF obs: sum over alleles [1..]; REF is allele 0.
            let nro = cols.allele_obs_count[(lo + 1).min(hi)..hi]
                .iter()
                .fold(0u32, |acc, &c| acc.saturating_add(c));
            nonref_obs.push(nro);
            // REF span: byte length of allele 0's sequence.
            let ref_span = cols.allele_seq_offsets[lo + 1] - cols.allele_seq_offsets[lo];
            ref_spans.push(ref_span);
        }

        Ok(Some(Self {
            chrom_id,
            positions,
            nonref_obs,
            ref_spans,
            allele_offsets,
            windowed_gc,
            windowed_coverage,
            scalar,
            seq,
            chain_ids,
        }))
    }

    /// Number of (region-clamped) records.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// 1-based positions, one per record (light, cached). Test-only accessor
    /// (the eager-decode oracle); production reads the field directly.
    #[cfg(test)]
    pub fn positions(&self) -> &[u32] {
        &self.positions
    }

    /// Per-record sum of non-REF allele obs (light, cached). Test-only accessor.
    #[cfg(test)]
    pub fn nonref_obs(&self) -> &[u32] {
        &self.nonref_obs
    }

    /// Per-record REF-allele reference span in bases (light, cached). Test-only
    /// accessor.
    #[cfg(test)]
    pub fn ref_spans(&self) -> &[u32] {
        &self.ref_spans
    }

    /// Number of alleles at record `record_idx`. Test-only accessor.
    #[cfg(test)]
    pub fn n_alleles_at(&self, record_idx: usize) -> usize {
        let lo = self.allele_offsets[record_idx] as usize;
        let hi = self.allele_offsets[record_idx + 1] as usize;
        hi - lo
    }

    /// An empty chunk for `chrom_id`, ready to accumulate compacted records
    /// via [`append_range`](Self::append_range).
    ///
    /// Used by the producer to gather a chunk's variable rows (across the
    /// sample's buffered segments) into one columnar object shipped to the
    /// caller, which rebuilds the records ([`records_all`](Self::records_all)).
    /// The light fold columns (`nonref_obs` / `ref_spans`) are left empty: a
    /// compacted chunk lives *past* the fold, so only the record-building
    /// columns (`positions` / `allele_offsets` / heavy) are populated.
    pub fn empty_compacted(chrom_id: u32) -> Self {
        Self {
            chrom_id,
            positions: Vec::new(),
            nonref_obs: Vec::new(),
            ref_spans: Vec::new(),
            allele_offsets: vec![0],
            windowed_gc: Vec::new(),
            windowed_coverage: Vec::new(),
            scalar: AlleleScalarColumns::default(),
            seq: AlleleSeqColumns::empty(),
            chain_ids: AlleleChainIdColumns::empty(),
        }
    }

    /// Per-record centred-window GC, one value per record (aligned with
    /// [`positions`](Self::positions)); empty for a legacy `.psp` without the
    /// windowed columns. Read by the hidden-paralog score's window gather.
    pub fn windowed_gc(&self) -> &[f32] {
        &self.windowed_gc
    }

    /// Per-record centred-window mean coverage; mirror of
    /// [`windowed_gc`](Self::windowed_gc).
    pub fn windowed_coverage(&self) -> &[f32] {
        &self.windowed_coverage
    }

    /// 1-based positions, one per record — the key the window gather binary-
    /// searches to align a called locus to this sample's windowed values.
    pub fn positions_all(&self) -> &[u32] {
        &self.positions
    }

    /// Append `src`'s records whose position is in `[lo, hi)` onto `self` — the
    /// cohort-chunk slice over a **finalised** (variable-only) two-phase
    /// segment. `src`'s positions are ascending, so the range is contiguous;
    /// copies the heavy columns for that row range (no decode, no filter). The
    /// producer slices each `Ready` segment overlapping `[chunk_start, cut)`
    /// this way to build a cohort chunk's per-sample columns.
    pub fn append_range(&mut self, src: &SamplePspChunk, lo: u32, hi: u32) {
        let r_lo = src.positions.partition_point(|&p| p < lo);
        let r_hi = src.positions.partition_point(|&p| p < hi);
        for r in r_lo..r_hi {
            let a_lo = src.allele_offsets[r] as usize;
            let a_hi = src.allele_offsets[r + 1] as usize;
            self.positions.push(src.positions[r]);
            // Windowed columns stay aligned 1:1 with `positions`; NaN fills a
            // source that predates them (keeps `len() == positions.len()`).
            self.windowed_gc
                .push(src.windowed_gc.get(r).copied().unwrap_or(f32::NAN));
            self.windowed_coverage
                .push(src.windowed_coverage.get(r).copied().unwrap_or(f32::NAN));
            self.scalar.extend_from_range(&src.scalar, a_lo..a_hi);
            self.seq.extend_from_range(&src.seq, a_lo..a_hi);
            self.chain_ids.extend_from_range(&src.chain_ids, a_lo..a_hi);
            let cum = self.allele_offsets.last().copied().unwrap_or(0) + (a_hi - a_lo) as u32;
            self.allele_offsets.push(cum);
        }
    }

    /// Reconstruct every record (a compacted chunk holds only kept rows) —
    /// the columns→records boundary, now run **on the caller**. Equivalent to
    /// `records_for(&vec![true; len()])` without allocating the mask. This is
    /// the per-allele-allocation work the record-building lever moves off the
    /// producer's single thread.
    pub fn records_all(&self) -> Vec<PileupRecord> {
        let mut out = Vec::with_capacity(self.len());
        for r in 0..self.len() {
            let lo = self.allele_offsets[r] as usize;
            let hi = self.allele_offsets[r + 1] as usize;
            let mut alleles = Vec::with_capacity(hi - lo);
            for a in lo..hi {
                alleles.push(AlleleObservation::new(
                    self.seq.slice_at(a).to_vec(),
                    self.scalar.support_at(a),
                    chain_ids_for_allele(&self.chain_ids, a, lo),
                ));
            }
            out.push(PileupRecord::new(self.chrom_id, self.positions[r], alleles));
        }
        out
    }

    /// Reconstruct owned [`PileupRecord`]s for the `keep` records (copies via
    /// `slice_at`/`support_at`). Test-only — part of the eager-decode oracle;
    /// production rebuilds records via [`records_all`](Self::records_all) off the
    /// two-phase path. `keep.len()` must equal [`len`](Self::len).
    #[cfg(test)]
    pub fn records_for(&self, keep: &[bool]) -> Vec<PileupRecord> {
        debug_assert_eq!(keep.len(), self.len(), "keep mask must cover every record");
        let mut out = Vec::new();
        for (r, &k) in keep.iter().enumerate() {
            if !k {
                continue;
            }
            let lo = self.allele_offsets[r] as usize;
            let hi = self.allele_offsets[r + 1] as usize;
            let mut alleles = Vec::with_capacity(hi - lo);
            for a in lo..hi {
                alleles.push(AlleleObservation::new(
                    self.seq.slice_at(a).to_vec(),
                    self.scalar.support_at(a),
                    chain_ids_for_allele(&self.chain_ids, a, lo),
                ));
            }
            out.push(PileupRecord::new(self.chrom_id, self.positions[r], alleles));
        }
        out
    }

    /// Move out the scalar columns for the `keep` records (appendix
    /// §A typed getter). `keep.len()` must equal [`len`](Self::len);
    /// the column is emptied, so a second call yields nothing. Test-only (the
    /// eager-decode oracle's column-selective seam).
    #[cfg(test)]
    pub fn take_scalar(&mut self, keep: &[bool]) -> AlleleScalarColumns {
        debug_assert_eq!(keep.len(), self.len(), "keep mask must cover every record");
        let src = std::mem::take(&mut self.scalar);
        let mut out = AlleleScalarColumns::default();
        for (r, &k) in keep.iter().enumerate() {
            if k {
                let lo = self.allele_offsets[r] as usize;
                let hi = self.allele_offsets[r + 1] as usize;
                out.extend_from_range(&src, lo..hi);
            }
        }
        out
    }

    /// Move out the allele-sequence column for the `keep` records.
    ///
    /// Call-once (per the appendix §A "fetched at most once" contract): the
    /// column is moved out, so a second call on the same chunk is misuse.
    /// Test-only (the eager-decode oracle's column-selective seam).
    #[cfg(test)]
    pub fn take_seq(&mut self, keep: &[bool]) -> AlleleSeqColumns {
        debug_assert_eq!(keep.len(), self.len(), "keep mask must cover every record");
        let src = std::mem::replace(&mut self.seq, AlleleSeqColumns::empty());
        let mut out = AlleleSeqColumns::empty();
        for (r, &k) in keep.iter().enumerate() {
            if k {
                let lo = self.allele_offsets[r] as usize;
                let hi = self.allele_offsets[r + 1] as usize;
                out.extend_from_range(&src, lo..hi);
            }
        }
        out
    }

    /// Move out the allele-chain-id column for the `keep` records. Test-only
    /// (the eager-decode oracle's column-selective seam).
    #[cfg(test)]
    pub fn take_chain_ids(&mut self, keep: &[bool]) -> AlleleChainIdColumns {
        debug_assert_eq!(keep.len(), self.len(), "keep mask must cover every record");
        let src = std::mem::replace(&mut self.chain_ids, AlleleChainIdColumns::empty());
        let mut out = AlleleChainIdColumns::empty();
        for (r, &k) in keep.iter().enumerate() {
            if k {
                let lo = self.allele_offsets[r] as usize;
                let hi = self.allele_offsets[r + 1] as usize;
                out.extend_from_range(&src, lo..hi);
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// TwoPhaseSegment — column-selective decode (light eager, heavy deferred)
// ---------------------------------------------------------------------------

/// Append `full[lo..hi]` for every kept row's allele range to `out` — the
/// per-allele scalar compaction used by [`TwoPhaseSegment::set_variable_rows`].
fn extend_kept<T: Copy>(out: &mut Vec<T>, full: &[T], ranges: &[(usize, usize)]) {
    for &(lo, hi) in ranges {
        out.extend_from_slice(&full[lo..hi]);
    }
}

/// One block's heavy columns inflated to **full** (all rows) — the decode-once
/// cache for a segment that straddles a chunk cut. The cohort producer slices
/// the variable rows for each overlapping chunk out of this
/// ([`compact_rows`](TwoPhaseSegment::compact_rows)) instead of re-inflating the
/// zstd blobs per chunk. `scalar` carries all 7 per-allele stats (including the
/// light `num_obs`); `seq` / `chain_ids` are the full-block ragged columns.
pub(crate) struct DecodedHeavy {
    scalar: AlleleScalarColumns,
    seq: AlleleSeqColumns,
    chain_ids: AlleleChainIdColumns,
    /// Full-block per-record windowed columns (NaN where absent), sliced to the
    /// kept rows by [`compact_rows`](TwoPhaseSegment::compact_rows). Pre-sized to
    /// `n_records` NaN so a `.psp` without the columns still slices cleanly.
    windowed_gc: Vec<f32>,
    windowed_coverage: Vec<f32>,
}

/// A `.psp` segment decoded in two phases (the column-selective producer read):
/// the fold columns are materialised; the per-allele heavy columns stay
/// compressed in `retained` until they are inflated + compacted to the kept
/// (variable) rows, yielding a variable-only [`SamplePspChunk`].
///
/// Two compaction paths (the cohort producer picks per segment):
/// - **fully-consumed** segment (sits wholly below the chunk cut, used once) —
///   [`set_variable_rows`](Self::set_variable_rows) inflates one column at a
///   time and frees it (transient ≈ one inflated column);
/// - **straddler** (extends past the cut, so the *next* chunk needs it too) —
///   [`inflate_all`](Self::inflate_all) decodes every heavy column once into the
///   [`decoded`](Self::decoded) cache; both chunks then slice their rows via
///   [`compact_rows`](Self::compact_rows) with **no** repeat zstd decode.
// Constructed by `SamplePspReader::next_two_phase`, folded over via the light
// accessors, then compacted by the cohort producer's compact stage.
pub(crate) struct TwoPhaseSegment {
    chrom_id: u32,
    positions: Vec<u32>,      // 1-based, full block
    ref_spans: Vec<u32>,      // allele-0 ref span, full block (fold)
    nonref_obs: Vec<u32>,     // Σ non-REF obs, full block (fold)
    allele_offsets: Vec<u32>, // CSR record→allele, full block
    n_total_alleles: usize,
    allele_seq_len: Vec<u64>, // full block — rebuilds the seq/chain CSR offsets
    allele_obs_count: Vec<u32>, // full block (light; supplies the num_obs stat)
    retained: Vec<RetainedColumn>,
    /// Decode-once cache for a straddling segment (see the type doc). Empty for
    /// fully-consumed segments. `OnceLock` so the shared `Arc<TwoPhaseSegment>`
    /// stays `Sync` and the first compaction populates it for the next.
    decoded: std::sync::OnceLock<Box<DecodedHeavy>>,
}

impl TwoPhaseSegment {
    /// Build from a two-phase-decoded block: delta-decode positions and derive
    /// the fold columns (ref span, Σ non-REF obs) from the light columns; keep
    /// the deferred columns compressed. **Whole block, no region clamp** — the
    /// cohort fold's `[next_chunk_start, watermark]` bounds already exclude
    /// out-of-interval rows, and such rows never enter the variable mask.
    pub fn from_two_phase_block(tp: TwoPhaseBlock) -> Result<Self, PspReadError> {
        let n = tp.n_records as usize;
        let mut positions = Vec::with_capacity(n);
        let mut p = tp.first_pos;
        for i in 0..n {
            if i > 0 {
                p = advance_pos(p, tp.delta_pos[i], i)?;
            }
            positions.push(p);
        }
        let mut allele_offsets = Vec::with_capacity(n + 1);
        allele_offsets.push(0u32);
        let mut cum = 0u32;
        for &na in &tp.n_alleles {
            cum += na as u32;
            allele_offsets.push(cum);
        }
        let mut ref_spans = Vec::with_capacity(n);
        let mut nonref_obs = Vec::with_capacity(n);
        for r in 0..n {
            let lo = allele_offsets[r] as usize;
            let hi = allele_offsets[r + 1] as usize;
            ref_spans.push(tp.allele_seq_len[lo] as u32);
            let nro = tp.allele_obs_count[(lo + 1).min(hi)..hi]
                .iter()
                .fold(0u32, |acc, &c| acc.saturating_add(c));
            nonref_obs.push(nro);
        }
        Ok(Self {
            chrom_id: tp.chrom_id,
            positions,
            ref_spans,
            nonref_obs,
            allele_offsets,
            n_total_alleles: tp.n_total_alleles as usize,
            allele_seq_len: tp.allele_seq_len,
            allele_obs_count: tp.allele_obs_count,
            retained: tp.retained,
            decoded: std::sync::OnceLock::new(),
        })
    }

    /// 1-based positions, one per record (the fold's key).
    pub fn positions(&self) -> &[u32] {
        &self.positions
    }
    /// Per-record REF-allele reference span (fold reach).
    pub fn ref_spans(&self) -> &[u32] {
        &self.ref_spans
    }
    /// Per-record Σ non-REF obs (fold AC).
    pub fn nonref_obs(&self) -> &[u32] {
        &self.nonref_obs
    }

    /// Per-record total observed depth = REF-allele obs (allele 0) + Σ non-REF
    /// obs — both **light** columns (`allele_obs_count` and the folded
    /// `nonref_obs`), so this is a couple of additions off already-decoded data,
    /// no heavy-column inflate. Feeds the paralog window-depth accumulator (S6c);
    /// read-only, it touches no fold state. Parallel to [`positions`](Self::positions).
    pub fn depth_at(&self, record_idx: usize) -> u32 {
        let ref_obs = self
            .allele_obs_count
            .get(self.allele_offsets[record_idx] as usize)
            .copied()
            .unwrap_or(0);
        ref_obs.saturating_add(self.nonref_obs[record_idx])
    }
    /// The segment's last 1-based position (its true extent — for the cohort
    /// watermark / finalisation check), or `None` when empty.
    pub fn last_pos(&self) -> Option<u32> {
        self.positions.last().copied()
    }

    /// Phase 2 of the two-phase decode: inflate the deferred columns **one at a
    /// time**, compacting each to the `keep` rows and freeing the inflated
    /// buffer before the next (resident transient ≈ one inflated column).
    /// Returns a variable-only
    /// [`SamplePspChunk`] ready for record building. `keep.len()` must equal
    /// [`len`](Self::len). **Non-consuming** — the producer partial-inflates a
    /// not-yet-finalisable straddling segment with this, then re-inflates it on
    /// the next chunk; a finalisable segment is consumed via [`finalize`](Self::finalize).
    pub fn set_variable_rows(
        &self,
        keep: &[bool],
        decompressor: &mut zstd::bulk::Decompressor<'static>,
        compressed_scratch: &mut Vec<u8>,
        decompressed_scratch: &mut Vec<u8>,
    ) -> Result<SamplePspChunk, PspReadError> {
        debug_assert_eq!(
            keep.len(),
            self.positions.len(),
            "keep mask covers every record"
        );
        let n_records = self.positions.len();

        // Kept positions + re-based CSR + each kept row's full-block allele range.
        // `kept_rows` records the full-block record index of each kept row, so
        // the per-record windowed columns can be compacted the same way.
        let mut positions = Vec::new();
        let mut allele_offsets = vec![0u32];
        let mut ranges: Vec<(usize, usize)> = Vec::new();
        let mut kept_rows: Vec<usize> = Vec::new();
        let mut cum = 0u32;
        for (r, &k) in keep.iter().enumerate() {
            if k {
                let lo = self.allele_offsets[r] as usize;
                let hi = self.allele_offsets[r + 1] as usize;
                positions.push(self.positions[r]);
                ranges.push((lo, hi));
                kept_rows.push(r);
                cum += (hi - lo) as u32;
                allele_offsets.push(cum);
            }
        }

        // `num_obs` is the light obs-count column (already decoded, not deferred).
        let mut scalar = AlleleScalarColumns::default();
        extend_kept(&mut scalar.num_obs, &self.allele_obs_count, &ranges);

        // Inflate + compact each deferred column, freeing between (one resident).
        let mut seq = AlleleSeqColumns::empty();
        let mut chain_ids = AlleleChainIdColumns::empty();
        // Per-record windowed columns, compacted to the kept rows. NaN-filled to
        // `kept_rows.len()` up front so a `.psp` without these columns still
        // yields the `len() == positions.len()` invariant the score relies on.
        let mut windowed_gc = vec![f32::NAN; kept_rows.len()];
        let mut windowed_coverage = vec![f32::NAN; kept_rows.len()];
        // Per-column inflate scratch: seq bytes/offsets, chain-id data/offsets.
        let (mut seq_bytes, mut seq_offsets, mut chain_data, mut chain_offsets) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for rc in &self.retained {
            let col = inflate_retained_column(
                &rc.blob,
                &rc.entry,
                n_records,
                self.n_total_alleles,
                Some(&self.allele_seq_len),
                decompressor,
                compressed_scratch,
                decompressed_scratch,
                &mut seq_bytes,
                &mut seq_offsets,
                &mut chain_data,
                &mut chain_offsets,
            )?;
            match col {
                Some(DecodedColumn::AlleleQSumLog(v)) => {
                    extend_kept(&mut scalar.q_sum, &v, &ranges)
                }
                Some(DecodedColumn::AlleleFwdCount(v)) => extend_kept(&mut scalar.fwd, &v, &ranges),
                Some(DecodedColumn::AllelePlacedLeftCount(v)) => {
                    extend_kept(&mut scalar.placed_left, &v, &ranges)
                }
                Some(DecodedColumn::AllelePlacedStartCount(v)) => {
                    extend_kept(&mut scalar.placed_start, &v, &ranges)
                }
                Some(DecodedColumn::AlleleMapqSum(v)) => {
                    extend_kept(&mut scalar.mapq_sum, &v, &ranges)
                }
                Some(DecodedColumn::AlleleMapqSumSq(v)) => {
                    extend_kept(&mut scalar.mapq_sum_sq, &v, &ranges)
                }
                Some(DecodedColumn::AlleleSeq) => {
                    let full = AlleleSeqColumns {
                        offsets: std::mem::take(&mut seq_offsets),
                        bytes: std::mem::take(&mut seq_bytes),
                    };
                    for &(lo, hi) in &ranges {
                        seq.extend_from_range(&full, lo..hi);
                    }
                }
                Some(DecodedColumn::AlleleChainIds) => {
                    let full = AlleleChainIdColumns {
                        offsets: std::mem::take(&mut chain_offsets),
                        ids: std::mem::take(&mut chain_data),
                    };
                    for &(lo, hi) in &ranges {
                        chain_ids.extend_from_range(&full, lo..hi);
                    }
                }
                // Per-record windowed columns: compact by kept record index
                // (not allele range) into the pre-sized NaN buffers.
                Some(DecodedColumn::WindowedGc(v)) => {
                    for (slot, &r) in windowed_gc.iter_mut().zip(&kept_rows) {
                        *slot = v[r];
                    }
                }
                Some(DecodedColumn::WindowedCoverage(v)) => {
                    for (slot, &r) in windowed_coverage.iter_mut().zip(&kept_rows) {
                        *slot = v[r];
                    }
                }
                // Light columns are never retained; an unknown optional → None.
                _ => {}
            }
        }

        Ok(SamplePspChunk {
            chrom_id: self.chrom_id,
            positions,
            nonref_obs: Vec::new(),
            ref_spans: Vec::new(),
            allele_offsets,
            windowed_gc,
            windowed_coverage,
            scalar,
            seq,
            chain_ids,
        })
    }

    /// The straddler decode-once cache, if already populated.
    pub(crate) fn decoded_cache(&self) -> Option<&DecodedHeavy> {
        self.decoded.get().map(|b| &**b)
    }

    /// Store the full-block decoded heavy columns for reuse by the next chunk's
    /// compaction (no-op if another caller already set it — same data).
    pub(crate) fn cache_decoded(&self, decoded: DecodedHeavy) {
        let _ = self.decoded.set(Box::new(decoded));
    }

    /// Inflate **every** heavy column to full-block form (all rows) — the
    /// decode-once step for a straddling segment. Unlike
    /// [`set_variable_rows`](Self::set_variable_rows) (which inflates one column
    /// at a time and immediately compacts to the kept rows), this keeps the
    /// whole decoded block so [`compact_rows`](Self::compact_rows) can slice it
    /// for each chunk that overlaps the segment with no repeat zstd decode.
    pub(crate) fn inflate_all(
        &self,
        decompressor: &mut zstd::bulk::Decompressor<'static>,
        compressed_scratch: &mut Vec<u8>,
        decompressed_scratch: &mut Vec<u8>,
    ) -> Result<DecodedHeavy, PspReadError> {
        let n_records = self.positions.len();
        let mut scalar = AlleleScalarColumns {
            num_obs: self.allele_obs_count.clone(),
            ..Default::default()
        };
        let mut seq = AlleleSeqColumns::empty();
        let mut chain_ids = AlleleChainIdColumns::empty();
        let mut windowed_gc = vec![f32::NAN; n_records];
        let mut windowed_coverage = vec![f32::NAN; n_records];
        let (mut seq_bytes, mut seq_offsets, mut chain_data, mut chain_offsets) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for rc in &self.retained {
            let col = inflate_retained_column(
                &rc.blob,
                &rc.entry,
                n_records,
                self.n_total_alleles,
                Some(&self.allele_seq_len),
                decompressor,
                compressed_scratch,
                decompressed_scratch,
                &mut seq_bytes,
                &mut seq_offsets,
                &mut chain_data,
                &mut chain_offsets,
            )?;
            match col {
                Some(DecodedColumn::AlleleQSumLog(v)) => scalar.q_sum = v,
                Some(DecodedColumn::AlleleFwdCount(v)) => scalar.fwd = v,
                Some(DecodedColumn::AllelePlacedLeftCount(v)) => scalar.placed_left = v,
                Some(DecodedColumn::AllelePlacedStartCount(v)) => scalar.placed_start = v,
                Some(DecodedColumn::AlleleMapqSum(v)) => scalar.mapq_sum = v,
                Some(DecodedColumn::AlleleMapqSumSq(v)) => scalar.mapq_sum_sq = v,
                Some(DecodedColumn::AlleleSeq) => {
                    seq = AlleleSeqColumns {
                        offsets: std::mem::take(&mut seq_offsets),
                        bytes: std::mem::take(&mut seq_bytes),
                    };
                }
                Some(DecodedColumn::AlleleChainIds) => {
                    chain_ids = AlleleChainIdColumns {
                        offsets: std::mem::take(&mut chain_offsets),
                        ids: std::mem::take(&mut chain_data),
                    };
                }
                Some(DecodedColumn::WindowedGc(v)) => windowed_gc = v,
                Some(DecodedColumn::WindowedCoverage(v)) => windowed_coverage = v,
                // Light columns are never retained; an unknown optional → None.
                _ => {}
            }
        }
        Ok(DecodedHeavy {
            scalar,
            seq,
            chain_ids,
            windowed_gc,
            windowed_coverage,
        })
    }

    /// Slice the `keep` (variable) rows out of an already-decoded full block into
    /// a variable-only [`SamplePspChunk`]. Byte-identical to
    /// [`set_variable_rows`](Self::set_variable_rows) for the same `keep` (same
    /// kept-row CSR + per-allele copies) — it just reads from `decoded` instead
    /// of inflating. `keep.len()` must equal the block's record count.
    pub(crate) fn compact_rows(&self, decoded: &DecodedHeavy, keep: &[bool]) -> SamplePspChunk {
        debug_assert_eq!(keep.len(), self.positions.len(), "keep covers every record");
        let mut positions = Vec::new();
        let mut allele_offsets = vec![0u32];
        let mut ranges: Vec<(usize, usize)> = Vec::new();
        let mut windowed_gc = Vec::new();
        let mut windowed_coverage = Vec::new();
        let mut cum = 0u32;
        for (r, &k) in keep.iter().enumerate() {
            if k {
                let lo = self.allele_offsets[r] as usize;
                let hi = self.allele_offsets[r + 1] as usize;
                positions.push(self.positions[r]);
                ranges.push((lo, hi));
                windowed_gc.push(decoded.windowed_gc[r]);
                windowed_coverage.push(decoded.windowed_coverage[r]);
                cum += (hi - lo) as u32;
                allele_offsets.push(cum);
            }
        }
        let mut scalar = AlleleScalarColumns::default();
        let mut seq = AlleleSeqColumns::empty();
        let mut chain_ids = AlleleChainIdColumns::empty();
        for &(lo, hi) in &ranges {
            scalar.extend_from_range(&decoded.scalar, lo..hi);
            seq.extend_from_range(&decoded.seq, lo..hi);
            chain_ids.extend_from_range(&decoded.chain_ids, lo..hi);
        }
        SamplePspChunk {
            chrom_id: self.chrom_id,
            positions,
            nonref_obs: Vec::new(),
            ref_spans: Vec::new(),
            allele_offsets,
            windowed_gc,
            windowed_coverage,
            scalar,
            seq,
            chain_ids,
        }
    }
}

// ---------------------------------------------------------------------------
// SamplePspReader
// ---------------------------------------------------------------------------

/// Per-sample segment reader bound to one chromosome region (appendix §A
/// `[RENAME ColumnSpanReader]`). Serves the records of `[region_start,
/// region_end]` (1-based inclusive) one **psp segment at a time**, in
/// ascending position order.
pub struct SamplePspReader<R: Read + Seek> {
    blocks: BlockColumnReader<R>,
    chrom_id: u32,
    /// Inclusive upper region bound; records past it are not served.
    region_end: u32,
    /// Inclusive lower region bound; earlier records (region-edge block
    /// only) are skipped.
    floor: u32,
}

impl<R: Read + Seek> SamplePspReader<R> {
    /// Open a reader over `[region_start, region_end]` of `chrom_id`,
    /// taking ownership of `reader`.
    pub fn new(reader: PspReader<R>, chrom_id: u32, region_start: u32, region_end: u32) -> Self {
        let mut blocks = reader.into_column_blocks();
        blocks.seek_to(chrom_id, region_start, region_end);
        Self {
            blocks,
            chrom_id,
            region_end,
            floor: region_start,
        }
    }

    /// Re-point at a new region, reusing the decode buffers (the cohort
    /// producer resets one reader per sample at every covered-interval /
    /// chromosome boundary so the zstd context + column slabs persist).
    pub fn reset(&mut self, chrom_id: u32, region_start: u32, region_end: u32) {
        self.blocks.seek_to(chrom_id, region_start, region_end);
        self.chrom_id = chrom_id;
        self.region_end = region_end;
        self.floor = region_start;
    }

    /// The sample's on-disk block index — the producer unions these across
    /// samples to find the covered intervals per chromosome.
    pub fn block_index(&self) -> &[BlockIndexEntry] {
        self.blocks.index()
    }

    /// Inclusive end of the next segment this reader serves **without a new
    /// decode** (the next block's `last_pos`), clamped to `region_end`.
    /// `None` once the region is exhausted. Test-only — paired with the eager
    /// [`next_chunk`](Self::next_chunk) oracle; production advances the cohort
    /// via the two-phase [`next_two_phase`](Self::next_two_phase) path.
    #[cfg(test)]
    pub fn peek_next_span(&self) -> Option<u32> {
        let entry = self.blocks.peek_block()?;
        if entry.chrom_id != self.chrom_id || entry.first_pos > self.region_end {
            return None;
        }
        Some(entry.last_pos.min(self.region_end))
    }

    /// Decode and return the next in-region segment as a [`SamplePspChunk`],
    /// advancing the cursor past it. `None` once the region is exhausted.
    /// Blocks whose records all fall below `floor` are skipped transparently.
    /// Test-only — the eager whole-segment decode oracle; production uses the
    /// column-selective [`next_two_phase`](Self::next_two_phase).
    #[cfg(test)]
    pub fn next_chunk(&mut self) -> Result<Option<SamplePspChunk>, PspReadError> {
        loop {
            match self.blocks.peek_block() {
                Some(entry)
                    if entry.chrom_id == self.chrom_id && entry.first_pos <= self.region_end => {}
                _ => return Ok(None),
            }
            let loaded = self.blocks.load_current()?;
            debug_assert!(loaded, "peek_block was Some, so the cursor is in range");
            let chunk = {
                let cols = self.blocks.columns().expect("block just loaded");
                SamplePspChunk::from_block(&cols, self.chrom_id, self.floor, self.region_end)?
            };
            self.blocks.advance();
            match chunk {
                Some(c) => return Ok(Some(c)),
                // Whole block was below the region floor — keep scanning.
                None => continue,
            }
        }
    }

    /// Decode and return the next in-region segment as a [`TwoPhaseSegment`]
    /// (column-selective: light columns eager, the rest retained compressed),
    /// advancing the cursor past it. `None` once the region is exhausted.
    ///
    /// Unlike [`next_chunk`](Self::next_chunk) the segment is **not**
    /// region-clamped — the deferred columns can't be sliced while compressed,
    /// and the producer's fold bounds (`[next_chunk_start, watermark]`) already
    /// exclude out-of-interval rows, which therefore never enter a variable mask.
    // Drives the cohort producer's column-selective read path.
    pub(crate) fn next_two_phase(&mut self) -> Result<Option<TwoPhaseSegment>, PspReadError> {
        match self.blocks.peek_block() {
            Some(entry)
                if entry.chrom_id == self.chrom_id && entry.first_pos <= self.region_end => {}
            _ => return Ok(None),
        }
        let tp = self.blocks.decode_current_two_phase()?;
        self.blocks.advance();
        match tp {
            Some(tp) => Ok(Some(TwoPhaseSegment::from_two_phase_block(tp)?)),
            None => Ok(None),
        }
    }
}

/// Advance a delta-encoded position, surfacing a `u32`-overflowing delta as
/// the same typed error the row reader raises (rather than a silent wrap).
fn advance_pos(prev: u32, delta: u64, entry: usize) -> Result<u32, PspReadError> {
    let delta_u32 = u32::try_from(delta).map_err(|_| PspReadError::ColumnElementDecode {
        column: "delta-pos".to_string(),
        entry,
        source: ScalarDecodeError::VarintOverflow,
    })?;
    Ok(prev.saturating_add(delta_u32))
}

#[cfg(test)]
mod tests {
    //! The light accessors and the typed getters must equal a full row
    //! decode (`PspReader::records`) of the same segment — the P1 gate.

    use std::io::Cursor;

    use super::*;
    use crate::pileup_record::PileupRecord;
    use crate::psp::test_fixtures::writer_header;
    use crate::psp::writer::PspWriter;
    use crate::var_calling::test_helpers::{allele, record};

    /// Multi-allele records at ascending positions, dense enough that a
    /// small block target yields several segments.
    fn fixture_records() -> Vec<PileupRecord> {
        (0..200u32)
            .map(|i| {
                let pos = 10 + i * 7;
                match i % 3 {
                    0 => record(pos, vec![allele(b"A", 12, -1.5, &[])]),
                    1 => record(
                        pos,
                        vec![
                            allele(b"A", 9, -1.0, &[]),
                            allele(b"T", 4, -2.0, &[(i as u64) + 1]),
                        ],
                    ),
                    _ => record(
                        pos,
                        vec![
                            allele(b"AC", 7, -0.5, &[]),
                            allele(b"AG", 3, -3.0, &[(i as u64) + 1, (i as u64) + 2]),
                            allele(b"A", 5, -1.2, &[]),
                        ],
                    ),
                }
            })
            .collect()
    }

    fn psp_bytes(records: &[PileupRecord], block_target: usize) -> Vec<u8> {
        let mut writer = PspWriter::new_with_block_target(
            Cursor::new(Vec::new()),
            writer_header(1),
            block_target,
        )
        .expect("writer opens");
        for r in records {
            writer.write_record(r).expect("write_record");
        }
        writer.finish().expect("finish").into_inner()
    }

    /// Step 3b: `TwoPhaseSegment::from_two_phase_block` + `set_variable_rows`
    /// (column-selective: light eager, heavy deferred) must build **identical**
    /// records to the eager `from_block` + `records_for` for the same `keep`.
    #[test]
    fn two_phase_segment_set_variable_rows_matches_eager() {
        let records = vec![
            record(
                10,
                vec![
                    allele(b"A", 5, -1.0, &[]),
                    allele(b"ACGT", 2, -1.5, &[7, 9]),
                ],
            ),
            record(
                25,
                vec![allele(b"T", 6, -1.0, &[3]), allele(b"TT", 1, -2.0, &[])],
            ),
            record(40, vec![allele(b"G", 4, -1.0, &[])]),
            record(
                55,
                vec![allele(b"C", 3, -1.0, &[]), allele(b"CC", 2, -1.0, &[5])],
            ),
        ];
        let bytes = psp_bytes(&records, 1024); // single block
        let keep = vec![true, false, true, true]; // drop row 1

        // Eager reference: from_block + records_for(keep).
        let eager_records = {
            let mut br = PspReader::new(Cursor::new(bytes.clone()))
                .unwrap()
                .into_column_blocks();
            assert!(br.load_current().unwrap());
            let cols = br.columns().unwrap();
            let chunk = SamplePspChunk::from_block(&cols, 0, 1, u32::MAX)
                .unwrap()
                .unwrap();
            chunk.records_for(&keep)
        };

        // Two-phase: from_two_phase_block + set_variable_rows(keep) + records_all.
        let tp_records = {
            let mut br = PspReader::new(Cursor::new(bytes))
                .unwrap()
                .into_column_blocks();
            let tp = br.decode_current_two_phase().unwrap().unwrap();
            let seg = TwoPhaseSegment::from_two_phase_block(tp).unwrap();
            let mut dz = crate::psp::block::new_column_decompressor().unwrap();
            let (mut cs, mut ds) = (Vec::new(), Vec::new());
            let chunk = seg
                .set_variable_rows(&keep, &mut dz, &mut cs, &mut ds)
                .unwrap();
            chunk.records_all()
        };

        assert_eq!(tp_records, eager_records);
        assert_eq!(tp_records.len(), 3); // kept rows 0, 2, 3
        // The light fold accessors match the eager chunk too.
        let seg = {
            let mut br = PspReader::new(Cursor::new(psp_bytes(&records, 1024)))
                .unwrap()
                .into_column_blocks();
            let tp = br.decode_current_two_phase().unwrap().unwrap();
            TwoPhaseSegment::from_two_phase_block(tp).unwrap()
        };
        assert_eq!(seg.positions(), &[10, 25, 40, 55]);
        assert_eq!(seg.ref_spans(), &[1, 1, 1, 1]); // allele-0 lengths: A,T,G,C
        assert_eq!(seg.nonref_obs(), &[2, 1, 0, 2]); // Σ non-REF obs
    }

    /// The per-record windowed columns must survive the round-trip and compact
    /// to the **kept** rows — in both compaction paths (`set_variable_rows` for a
    /// fully-consumed segment, `inflate_all` + `compact_rows` for a straddler),
    /// and byte-identically to the eager `from_block` oracle. Distinct finite
    /// values per row (the shared fixture defaults them to NaN) so a mis-indexed
    /// slice is caught.
    #[test]
    fn windowed_columns_compact_to_kept_rows() {
        let mut records = vec![
            record(10, vec![allele(b"A", 5, -1.0, &[])]),
            record(25, vec![allele(b"T", 6, -1.0, &[3])]),
            record(40, vec![allele(b"G", 4, -1.0, &[])]),
            record(55, vec![allele(b"C", 3, -1.0, &[])]),
        ];
        // Distinct, finite windowed values keyed to the position.
        for (i, r) in records.iter_mut().enumerate() {
            r.windowed_gc = 0.10 + i as f32 * 0.05;
            r.windowed_coverage = 100.0 + i as f32;
        }
        let bytes = psp_bytes(&records, 1024); // single block
        let keep = vec![true, false, true, true]; // drop row 1 (pos 25)
        let want_gc = [0.10f32, 0.20, 0.25]; // rows 0, 2, 3
        let want_cov = [100.0f32, 102.0, 103.0];

        let mut dz = crate::psp::block::new_column_decompressor().unwrap();
        let (mut cs, mut ds) = (Vec::new(), Vec::new());

        // Eager oracle: from_block keeps all rows; compare its whole-block
        // windowed columns to the source (no compaction) as a decode check.
        {
            let mut br = PspReader::new(Cursor::new(bytes.clone()))
                .unwrap()
                .into_column_blocks();
            assert!(br.load_current().unwrap());
            let cols = br.columns().unwrap();
            let chunk = SamplePspChunk::from_block(&cols, 0, 1, u32::MAX)
                .unwrap()
                .unwrap();
            assert_eq!(chunk.windowed_gc(), &[0.10, 0.15, 0.20, 0.25]);
            assert_eq!(chunk.windowed_coverage(), &[100.0, 101.0, 102.0, 103.0]);
        }

        // Path 1 — set_variable_rows (fully-consumed segment).
        {
            let mut br = PspReader::new(Cursor::new(bytes.clone()))
                .unwrap()
                .into_column_blocks();
            let tp = br.decode_current_two_phase().unwrap().unwrap();
            let seg = TwoPhaseSegment::from_two_phase_block(tp).unwrap();
            let chunk = seg
                .set_variable_rows(&keep, &mut dz, &mut cs, &mut ds)
                .unwrap();
            assert_eq!(chunk.windowed_gc(), &want_gc);
            assert_eq!(chunk.windowed_coverage(), &want_cov);
        }

        // Path 2 — inflate_all + compact_rows (straddler decode-once).
        {
            let mut br = PspReader::new(Cursor::new(bytes))
                .unwrap()
                .into_column_blocks();
            let tp = br.decode_current_two_phase().unwrap().unwrap();
            let seg = TwoPhaseSegment::from_two_phase_block(tp).unwrap();
            let decoded = seg.inflate_all(&mut dz, &mut cs, &mut ds).unwrap();
            let chunk = seg.compact_rows(&decoded, &keep);
            assert_eq!(chunk.windowed_gc(), &want_gc);
            assert_eq!(chunk.windowed_coverage(), &want_cov);
        }
    }

    /// Step 4 (reader side): draining `next_two_phase` (+ `set_variable_rows`
    /// keep-all per segment) over a whole-file region rebuilds **the same
    /// records** as draining the eager `next_chunk`, across multiple segments.
    #[test]
    fn sample_reader_next_two_phase_matches_next_chunk() {
        let records = fixture_records();
        let bytes = psp_bytes(&records, 2); // small block target → several segments

        let mut eager = Vec::new();
        {
            let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
            let mut sr = SamplePspReader::new(reader, 0, 1, u32::MAX);
            while let Some(chunk) = sr.next_chunk().unwrap() {
                eager.extend(chunk.records_all());
            }
        }

        let mut tp = Vec::new();
        {
            let reader = PspReader::new(Cursor::new(bytes)).unwrap();
            let mut sr = SamplePspReader::new(reader, 0, 1, u32::MAX);
            let mut dz = crate::psp::block::new_column_decompressor().unwrap();
            let (mut cs, mut ds) = (Vec::new(), Vec::new());
            while let Some(seg) = sr.next_two_phase().unwrap() {
                let keep = vec![true; seg.positions().len()];
                let chunk = seg
                    .set_variable_rows(&keep, &mut dz, &mut cs, &mut ds)
                    .unwrap();
                tp.extend(chunk.records_all());
            }
        }

        assert_eq!(tp, eager);
        assert!(!eager.is_empty());
    }

    /// `append_range(src, lo, hi)` slices a finalised (variable-only) segment by
    /// position — equivalent to `records_for` over the same row range.
    #[test]
    fn append_range_slices_a_compacted_segment_by_position() {
        let records = vec![
            record(
                10,
                vec![allele(b"A", 5, -1.0, &[]), allele(b"AC", 2, -1.5, &[7])],
            ),
            record(25, vec![allele(b"T", 6, -1.0, &[3])]),
            record(40, vec![allele(b"G", 4, -1.0, &[])]),
            record(
                55,
                vec![allele(b"C", 3, -1.0, &[]), allele(b"CC", 1, -2.0, &[])],
            ),
        ];
        let bytes = psp_bytes(&records, 1024);
        let chunk = {
            let mut br = PspReader::new(Cursor::new(bytes))
                .unwrap()
                .into_column_blocks();
            assert!(br.load_current().unwrap());
            let cols = br.columns().unwrap();
            SamplePspChunk::from_block(&cols, 0, 1, u32::MAX)
                .unwrap()
                .unwrap()
        };

        // Slice [25, 55) → positions 25, 40.
        let mut dst = SamplePspChunk::empty_compacted(0);
        dst.append_range(&chunk, 25, 55);
        let got = dst.records_all();

        let want = chunk.records_for(&[false, true, true, false]);
        assert_eq!(got, want);
        assert_eq!(got.iter().map(|r| r.pos).collect::<Vec<_>>(), vec![25, 40]);
    }

    /// Reference: full row decode of `[start, end]`.
    fn rows(bytes: &[u8], chrom: u32, start: u32, end: u32) -> Vec<PileupRecord> {
        let mut reader = PspReader::new(Cursor::new(bytes.to_vec())).expect("reader opens");
        reader
            .region_records(chrom, start, end)
            .map(|r| r.expect("record decodes"))
            .collect()
    }

    /// Drain every chunk of a reader over `[start, end]`, concatenating the
    /// light columns and the `keep`-all heavy columns into flat buffers
    /// matched against the row reference.
    fn drain_all(
        bytes: &[u8],
        chrom: u32,
        start: u32,
        end: u32,
    ) -> (
        Vec<u32>,
        Vec<u32>,
        Vec<u32>,
        AlleleScalarColumns,
        AlleleSeqColumns,
        AlleleChainIdColumns,
    ) {
        let reader = PspReader::new(Cursor::new(bytes.to_vec())).expect("reader opens");
        let mut sr = SamplePspReader::new(reader, chrom, start, end);

        let mut positions = Vec::new();
        let mut nonref_obs = Vec::new();
        let mut ref_spans = Vec::new();
        let mut scalar = AlleleScalarColumns::default();
        let mut seq = AlleleSeqColumns::empty();
        let mut chain = AlleleChainIdColumns::empty();

        while let Some(mut chunk) = sr.next_chunk().expect("next_chunk") {
            positions.extend_from_slice(chunk.positions());
            nonref_obs.extend_from_slice(chunk.nonref_obs());
            ref_spans.extend_from_slice(chunk.ref_spans());
            let keep = vec![true; chunk.len()];
            let cf = chunk.take_scalar(&keep);
            let cs = chunk.take_seq(&keep);
            let cc = chunk.take_chain_ids(&keep);
            let n = cf.len();
            scalar.extend_from_range(&cf, 0..n);
            seq.extend_from_range(&cs, 0..cs.len());
            chain.extend_from_range(&cc, 0..cc.len());
        }
        (positions, nonref_obs, ref_spans, scalar, seq, chain)
    }

    /// Build the expected flat columns straight from the row records.
    fn expected_from_rows(
        recs: &[PileupRecord],
    ) -> (
        Vec<u32>,
        Vec<u32>,
        Vec<u32>,
        AlleleScalarColumns,
        AlleleSeqColumns,
        AlleleChainIdColumns,
    ) {
        let mut positions = Vec::new();
        let mut nonref_obs = Vec::new();
        let mut ref_spans = Vec::new();
        let mut scalar = AlleleScalarColumns::default();
        let mut seq = AlleleSeqColumns::empty();
        let mut chain = AlleleChainIdColumns::empty();
        for r in recs {
            positions.push(r.pos);
            let nro = r.alleles[1..]
                .iter()
                .fold(0u32, |acc, a| acc.saturating_add(a.support.num_obs));
            nonref_obs.push(nro);
            ref_spans.push(r.alleles[0].seq.len() as u32);
            for a in &r.alleles {
                scalar.num_obs.push(a.support.num_obs);
                scalar.q_sum.push(a.support.q_sum);
                scalar.fwd.push(a.support.fwd);
                scalar.placed_left.push(a.support.placed_left);
                scalar.placed_start.push(a.support.placed_start);
                scalar.mapq_sum.push(a.support.mapq_sum);
                scalar.mapq_sum_sq.push(a.support.mapq_sum_sq);
                seq.push(&a.seq);
                chain.push(&a.chain_ids);
            }
        }
        (positions, nonref_obs, ref_spans, scalar, seq, chain)
    }

    #[test]
    fn whole_region_matches_row_decode() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        // Sanity: multi-segment.
        let reader = PspReader::new(Cursor::new(bytes.clone())).unwrap();
        assert!(
            reader.block_index().len() >= 3,
            "fixture should be multi-block"
        );

        let got = drain_all(&bytes, 0, 1, 1_000_000);
        let want = expected_from_rows(&rows(&bytes, 0, 1, 1_000_000));
        assert_eq!(got.0, want.0, "positions");
        assert_eq!(got.1, want.1, "nonref_obs");
        assert_eq!(got.2, want.2, "ref_spans");
        assert_eq!(got.3, want.3, "scalar");
        assert_eq!(got.4, want.4, "seq");
        assert_eq!(got.5, want.5, "chain_ids");
    }

    #[test]
    fn clamped_region_matches_row_decode() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        let (start, end) = (137, 941);
        let got = drain_all(&bytes, 0, start, end);
        let want = expected_from_rows(&rows(&bytes, 0, start, end));
        assert_eq!(got.0, want.0, "positions");
        assert_eq!(got.3, want.3, "scalar");
        assert_eq!(got.4, want.4, "seq");
        assert_eq!(got.5, want.5, "chain_ids");
        // The clamp really dropped records outside the window.
        assert!(got.0.len() < recs.len());
        assert_eq!(*got.0.first().unwrap(), *want.0.first().unwrap());
        assert!(*got.0.first().unwrap() >= start && *got.0.last().unwrap() <= end);
    }

    #[test]
    fn take_with_keep_subset_selects_those_records() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 4096); // single-ish block
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let mut sr = SamplePspReader::new(reader, 0, 1, 1_000_000);

        let mut chunk = sr.next_chunk().expect("next_chunk").expect("a chunk");
        // Keep only even-indexed records of this chunk.
        let keep: Vec<bool> = (0..chunk.len()).map(|i| i % 2 == 0).collect();
        let kept_positions: Vec<u32> = chunk
            .positions()
            .iter()
            .enumerate()
            .filter_map(|(i, &p)| (i % 2 == 0).then_some(p))
            .collect();
        let expected_alleles: usize = (0..chunk.len())
            .filter(|i| i % 2 == 0)
            .map(|i| chunk.n_alleles_at(i))
            .sum();

        let scalar = chunk.take_scalar(&keep);
        let seq = chunk.take_seq(&keep);
        let chain = chunk.take_chain_ids(&keep);
        assert_eq!(scalar.len(), expected_alleles, "scalar allele count");
        assert_eq!(seq.len(), expected_alleles, "seq allele count");
        assert_eq!(chain.len(), expected_alleles, "chain allele count");
        // Cross-check against a fresh decode of just those positions.
        let rows_all = rows(&psp_bytes(&recs, 4096), 0, 1, 1_000_000);
        let want: Vec<u32> = kept_positions.clone();
        let got: Vec<u32> = rows_all
            .iter()
            .filter(|r| want.contains(&r.pos))
            .map(|r| r.pos)
            .collect();
        assert_eq!(got, want);
    }

    #[test]
    fn peek_reports_ascending_segment_ends() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let mut sr = SamplePspReader::new(reader, 0, 1, 1_000_000);
        let mut peeks = Vec::new();
        while let Some(end) = sr.peek_next_span() {
            peeks.push(end);
            sr.next_chunk().unwrap();
        }
        assert!(peeks.len() >= 3, "multi-segment");
        assert!(
            peeks.windows(2).all(|w| w[0] < w[1]),
            "segment ends ascending"
        );
    }

    #[test]
    fn empty_region_yields_no_chunks() {
        let recs = fixture_records();
        let bytes = psp_bytes(&recs, 512);
        // Positions are 10, 17, 24, ...; (11, 16) holds none.
        let reader = PspReader::new(Cursor::new(bytes)).unwrap();
        let mut sr = SamplePspReader::new(reader, 0, 11, 16);
        assert!(sr.next_chunk().unwrap().is_none());
    }
}
