//! Chunk loader — read per-sample [`PileupRecord`] iterators into
//! columnar storage, apply the cohort-wide variant-position filter,
//! and compact the survivors into a [`MaterialisedChunk`] for the
//! pre-pass and worker.

use std::ops::Range;

use thiserror::Error;

use crate::pileup_record::PileupRecord;
use crate::var_calling::cohort_block::columns::{MaterialisedChunk, SampleColumns};

/// Reusable scratch buffers for one iteration of the chunk loader.
///
/// The driver owns one of these and threads `&mut` to [`load_chunk_from_iters`]
/// per chunk; nothing inside is meant to survive across multiple
/// iterations except as memory the next iteration overwrites. All
/// inner buffers are `clear()`-ed on every call so the driver does
/// not have to.
///
/// **Fields.**
/// - [`Self::raw_per_sample`]: one [`SampleColumns`] per sample,
///   pre-sized at construction. Accumulates each sample's
///   pre-filter records (carryover prefix + freshly-loaded records)
///   before the variant filter compacts them into the output chunk.
/// - [`Self::position_union`]: sorted, dedup'd 1-based position
///   timeline across every sample in the current chunk.
/// - [`Self::has_variant_at`]: parallel to [`Self::position_union`];
///   each entry is `true` iff at least one sample at that position
///   carries a non-reference allele with `num_obs > 0`. Used to
///   propagate "variant" up to the group each position belongs to.
/// - [`Self::max_ref_span_at`]: parallel to [`Self::position_union`];
///   the max `ref_span` across samples that have a record at the
///   position. Drives the grouping simulation's reach calculation.
/// - [`Self::is_kept`]: parallel to [`Self::position_union`]; the
///   result of the grouping simulation — `true` iff the position
///   lands in a provisional group that contains at least one variant
///   position. These are the positions whose records survive into
///   the output chunk's `per_sample` columns.
#[derive(Debug)]
pub struct ChunkLoadScratch {
    raw_per_sample: Vec<SampleColumns>,
    position_union: Vec<u32>,
    has_variant_at: Vec<bool>,
    max_ref_span_at: Vec<u32>,
    is_kept: Vec<bool>,
}

impl ChunkLoadScratch {
    /// Build scratch sized for a cohort of `n_samples` samples.
    pub fn with_n_samples(n_samples: usize) -> Self {
        Self {
            raw_per_sample: (0..n_samples).map(|_| SampleColumns::empty()).collect(),
            position_union: Vec::new(),
            has_variant_at: Vec::new(),
            max_ref_span_at: Vec::new(),
            is_kept: Vec::new(),
        }
    }

    /// Number of samples this scratch was sized for.
    pub fn n_samples(&self) -> usize {
        self.raw_per_sample.len()
    }

    /// Reset every internal buffer (raw per-sample columns, position
    /// union, variant decisions) while preserving allocated capacity.
    pub fn clear(&mut self) {
        for sample in &mut self.raw_per_sample {
            sample.clear();
        }
        self.position_union.clear();
        self.has_variant_at.clear();
        self.max_ref_span_at.clear();
        self.is_kept.clear();
    }
}

/// Errors surfaced by the chunk loader. Generic over `E`, the upstream
/// per-sample iterator's error type — typically `PspReadError` for the
/// production glue, or a unit-test error type for fixture-driven tests.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum ChunkLoadError<E> {
    /// A per-sample upstream iterator surfaced an error before the
    /// chunk's range was exhausted.
    #[error("failed to read pileup record for sample {sample_idx}")]
    UpstreamRead {
        sample_idx: usize,
        #[source]
        source: E,
    },

    /// `per_sample_iters.len()` did not match the cohort size the
    /// scratch was built for.
    #[error("per-sample iterator count {got} does not match scratch cohort size {expected}")]
    SampleCountMismatch { expected: usize, got: usize },

    /// `carryover.len()` did not match the cohort size the scratch
    /// was built for.
    #[error("carryover length {got} does not match scratch cohort size {expected}")]
    CarryoverLengthMismatch { expected: usize, got: usize },

    /// The chunk's requested range is empty or reversed
    /// (`start >= end`).
    #[error("chunk range {start}..{end} is empty or reversed")]
    InvalidRange { start: u32, end: u32 },

    /// An upstream iterator yielded a record whose `chrom_id` did
    /// not match the chunk's `chrom_id` — a PSP bug or a wiring
    /// mistake by the caller. The chunk loader does not cross
    /// chromosome boundaries.
    #[error(
        "sample {sample_idx} yielded record on chrom {got_chrom_id} \
         while chunk expects chrom {expected_chrom_id}"
    )]
    UnexpectedChromosome {
        sample_idx: usize,
        expected_chrom_id: u32,
        got_chrom_id: u32,
    },
}

/// Load one chunk from per-sample record iterators, apply the
/// cohort-wide variant-position filter, and compact the survivors
/// into `out`.
///
/// **Inputs.**
/// - `scratch`, `out`: caller-owned scratch and output chunk; both
///   are cleared at entry.
/// - `chrom_id`, `range`: the chunk's genomic window. Only records
///   with `range.start <= pos < range.end` are taken from the
///   upstream iterators.
/// - `per_sample_iters`: one iterator per sample, in cohort order.
///   Iterators are consumed up to the first record at
///   `pos >= range.end` (or until exhausted) — callers that want to
///   keep reading must hand in a fresh iterator next chunk.
/// - `carryover`: per-sample [`SampleColumns`] holding records from
///   the previous chunk's `>= safe_end` tail. These are *prepended*
///   to each sample's raw load. Cleared after consumption so the
///   caller can refill them from the next pre-pass's split.
///
/// **Algorithm.**
/// 1. Drain `carryover[s]` into `scratch.raw_per_sample[s]` (move
///    rows, no row materialisation).
/// 2. Pull records from `per_sample_iters[s]` while
///    `pos < range.end`, pushing into `scratch.raw_per_sample[s]`.
///    Records with `pos < range.start` are dropped (PSP's
///    `region_records` already clamps; this is the defensive
///    fallback against hand-rolled iterators).
/// 3. Build the sorted, dedup'd position union across all samples.
/// 4. For each position in the union, set `is_variant` true iff any
///    sample's record at that position carries a non-REF allele with
///    `num_obs > 0` (the walker's `alleles[0] == REF` invariant
///    makes this the fast path; no reference-base fetch needed).
/// 5. Compact: for each sample, walk its raw rows and the
///    `position_union`/`is_variant` arrays in parallel; push every
///    row whose position is marked variant into `out.per_sample[s]`.
/// 6. Set `out.chrom_id = chrom_id`, `out.range = range`,
///    `out.safe_end = range.end` (the pre-pass may revise it down),
///    leave `out.windows` empty for the pre-pass to populate.
pub fn load_chunk_from_iters<I, E>(
    scratch: &mut ChunkLoadScratch,
    out: &mut MaterialisedChunk,
    chrom_id: u32,
    range: Range<u32>,
    per_sample_iters: Vec<I>,
    carryover: &mut [SampleColumns],
) -> Result<(), ChunkLoadError<E>>
where
    I: IntoIterator<Item = Result<PileupRecord, E>>,
{
    let n_samples = scratch.n_samples();
    if per_sample_iters.len() != n_samples {
        return Err(ChunkLoadError::SampleCountMismatch {
            expected: n_samples,
            got: per_sample_iters.len(),
        });
    }
    if carryover.len() != n_samples {
        return Err(ChunkLoadError::CarryoverLengthMismatch {
            expected: n_samples,
            got: carryover.len(),
        });
    }
    if range.start >= range.end {
        return Err(ChunkLoadError::InvalidRange {
            start: range.start,
            end: range.end,
        });
    }

    scratch.clear();
    out.clear_data();
    out.chrom_id = chrom_id;
    out.range = range.clone();
    out.safe_end = range.end;
    // Make sure the output is sized for the cohort: a fresh
    // `MaterialisedChunk` may have come in with no samples slotted.
    if out.per_sample.len() != n_samples {
        out.per_sample.resize_with(n_samples, SampleColumns::empty);
    }

    // ── Step 1+2: raw load (carryover prefix + iterator pull) ──
    for (sample_idx, (iter, carry)) in per_sample_iters
        .into_iter()
        .zip(carryover.iter_mut())
        .enumerate()
    {
        let raw = &mut scratch.raw_per_sample[sample_idx];
        // Move carryover rows in first (preserve order).
        for row_idx in 0..carry.n_records() {
            raw.push_row_from(carry, row_idx);
        }
        carry.clear();

        for record in iter {
            let record =
                record.map_err(|source| ChunkLoadError::UpstreamRead { sample_idx, source })?;
            if record.chrom_id != chrom_id {
                return Err(ChunkLoadError::UnexpectedChromosome {
                    sample_idx,
                    expected_chrom_id: chrom_id,
                    got_chrom_id: record.chrom_id,
                });
            }
            if record.pos < range.start {
                continue;
            }
            if record.pos >= range.end {
                // First record past the right boundary: stop
                // pulling. Caller is responsible for any rows past
                // `range.end` (re-read next chunk; the iterator is
                // discarded).
                break;
            }
            raw.push_record(record);
        }
    }

    // ── Step 3: cohort-wide position union ──
    for raw in &scratch.raw_per_sample {
        scratch.position_union.extend_from_slice(&raw.positions);
    }
    scratch.position_union.sort_unstable();
    scratch.position_union.dedup();

    // ── Step 4: per-position scan — `has_variant_at` and `max_ref_span_at`. ──
    // For each cohort position we walk the samples that have a record
    // there exactly once, computing both predicates in one pass:
    //   - `has_variant_at[i]` — at least one sample's record carries a
    //     non-REF allele with `num_obs > 0`. Used by Step 5 to decide
    //     whether a provisional group has any variant evidence.
    //   - `max_ref_span_at[i]` — max `ref_span` across samples with a
    //     record. Drives the grouping simulation's reach calculation.
    scratch
        .has_variant_at
        .resize(scratch.position_union.len(), false);
    scratch
        .max_ref_span_at
        .resize(scratch.position_union.len(), 0);
    for (idx, &pos) in scratch.position_union.iter().enumerate() {
        let mut has_variant = false;
        let mut max_ref_span: u32 = 0;
        for raw in &scratch.raw_per_sample {
            if let Ok(row_idx) = raw.binary_search_position(pos) {
                let ref_span = raw.ref_span_at(row_idx);
                if ref_span > max_ref_span {
                    max_ref_span = ref_span;
                }
                if !has_variant && raw.has_observed_non_ref_allele_at(row_idx) {
                    has_variant = true;
                }
            }
        }
        scratch.has_variant_at[idx] = has_variant;
        scratch.max_ref_span_at[idx] = max_ref_span;
    }

    // ── Step 5: grouping simulation. ──
    // Walk the cohort-wide timeline in order, building provisional
    // groups under the same join rule the streaming grouper uses
    // (`pos <= group_end`, where `group_end` is the rolling max of
    // `pos + ref_span - 1` across the group). Mark every position in
    // a group that contains at least one variant position as kept;
    // positions in groups with no variants are dropped (matches the
    // streaming pipeline, where the per-group merger drops pure-REF
    // groups after merging).
    //
    // The point of this simulation — vs. the simpler per-position
    // "has any non-REF" check — is that the per-group merger gathers
    // per-(sample, allele) support from *every* position the group
    // span covers, including positions that are pure-REF in isolation
    // but inside the reach of another sample's MNP/DEL/INS. Dropping
    // those positions here would silently under-count REF evidence
    // for homref samples in multi-position groups.
    scratch.is_kept.resize(scratch.position_union.len(), false);
    let mut group_open = false;
    let mut group_start_idx: usize = 0;
    let mut group_end_pos: u32 = 0;
    let mut group_has_variant = false;
    for i in 0..scratch.position_union.len() {
        let pos = scratch.position_union[i];
        let reach = pos.saturating_add(scratch.max_ref_span_at[i].max(1)) - 1;

        if group_open && pos <= group_end_pos {
            if reach > group_end_pos {
                group_end_pos = reach;
            }
            if scratch.has_variant_at[i] {
                group_has_variant = true;
            }
        } else {
            if group_open && group_has_variant {
                for slot in scratch.is_kept[group_start_idx..i].iter_mut() {
                    *slot = true;
                }
            }
            group_open = true;
            group_start_idx = i;
            group_end_pos = reach;
            group_has_variant = scratch.has_variant_at[i];
        }
    }
    if group_open && group_has_variant {
        for slot in scratch.is_kept[group_start_idx..].iter_mut() {
            *slot = true;
        }
    }

    // ── Step 6: per-sample compact ──
    for (sample_idx, raw) in scratch.raw_per_sample.iter().enumerate() {
        let dst = &mut out.per_sample[sample_idx];
        // Parallel walk over `raw.positions` and `position_union`:
        // both are sorted ascending, so we never need to backtrack.
        let mut union_idx = 0_usize;
        for row_idx in 0..raw.n_records() {
            let pos = raw.position_at(row_idx);
            while union_idx < scratch.position_union.len()
                && scratch.position_union[union_idx] < pos
            {
                union_idx += 1;
            }
            debug_assert!(
                union_idx < scratch.position_union.len()
                    && scratch.position_union[union_idx] == pos,
                "raw position {pos} missing from union — invariant broken",
            );
            if scratch.is_kept[union_idx] {
                dst.push_row_from(raw, row_idx);
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var_calling::cohort_block::test_helpers::{
        allele, record, ref_obs, ref_plus_alt, ref_plus_unobserved_alt, run_loader,
    };

    #[test]
    fn loader_drops_positions_with_no_observed_non_ref_allele() {
        // Sample 0 has REF-only at pos 10; both samples carry REF +
        // unobserved ALT at pos 12; sample 1 has REF + observed ALT
        // at pos 14. Expect pos 14 in the output; pos 10 and pos 12
        // dropped.
        let s0 = vec![
            record(10, vec![ref_obs(20)]),
            record(12, ref_plus_unobserved_alt(15)),
        ];
        let s1 = vec![
            record(12, ref_plus_unobserved_alt(15)),
            record(14, ref_plus_alt(5, 6)),
        ];
        let chunk = run_loader(0, 1..100, vec![s0, s1], vec![SampleColumns::empty(); 2]);
        assert_eq!(chunk.range, 1..100);
        assert_eq!(chunk.safe_end, 100);
        assert_eq!(chunk.per_sample[0].n_records(), 0);
        assert_eq!(chunk.per_sample[1].n_records(), 1);
        assert_eq!(chunk.per_sample[1].position_at(0), 14);
    }

    #[test]
    fn loader_keeps_homref_position_inside_mnp_reach() {
        // Sample 0 has a 3-bp MNP at position 10 (REF=AAA, ALT=ACA).
        // Sample 1 has plain SNP-shape records at 10/11/12 — pure REF.
        // Cohort-wide: position 10 is variant (sample 0 has ALT). But
        // positions 11 and 12 are pure-REF in isolation (sample 0's
        // 3-bp record covers them; sample 1 has REF-only records).
        // The grouping simulation extends the kept set to 11 and 12
        // because they fall inside sample 0's record's reach (10..12).
        let s0 = vec![record(
            10,
            vec![
                allele(b"AAA", 5, -1.0, &[]),
                allele(b"ACA", 4, -1.0, &[]),
            ],
        )];
        let s1 = vec![
            record(10, vec![ref_obs(7)]),
            record(11, vec![ref_obs(7)]),
            record(12, vec![ref_obs(7)]),
        ];
        let chunk = run_loader(0, 1..100, vec![s0, s1], vec![SampleColumns::empty(); 2]);

        // Sample 0 only has a record at 10. Kept.
        assert_eq!(chunk.per_sample[0].n_records(), 1);
        assert_eq!(chunk.per_sample[0].position_at(0), 10);
        // Sample 1 has records at 10, 11, 12 — all three kept because
        // they're inside sample 0's variant record's reach.
        assert_eq!(chunk.per_sample[1].n_records(), 3);
        assert_eq!(chunk.per_sample[1].position_at(0), 10);
        assert_eq!(chunk.per_sample[1].position_at(1), 11);
        assert_eq!(chunk.per_sample[1].position_at(2), 12);
    }

    #[test]
    fn loader_keeps_homref_sample_at_variant_position() {
        // Position 50: sample 0 homref (REF only), sample 1 has the
        // ALT. The merger needs sample 0's record at 50 for joint
        // likelihoods, so it must survive the filter.
        let s0 = vec![record(50, vec![ref_obs(20)])];
        let s1 = vec![record(50, ref_plus_alt(2, 6))];
        let chunk = run_loader(0, 1..100, vec![s0, s1], vec![SampleColumns::empty(); 2]);
        assert_eq!(chunk.per_sample[0].n_records(), 1);
        assert_eq!(chunk.per_sample[1].n_records(), 1);
        assert_eq!(chunk.per_sample[0].position_at(0), 50);
        assert_eq!(chunk.per_sample[1].position_at(0), 50);
    }

    #[test]
    fn loader_prepends_carryover_in_position_order() {
        // Carryover from the previous chunk: sample 0 has a variant
        // at pos 5; sample 1 has nothing. New chunk range starts at
        // pos 10. Expect output to carry pos 5 + pos 12 in order
        // for sample 0.
        let mut carry_s0 = SampleColumns::empty();
        carry_s0.push_record(record(5, ref_plus_alt(3, 4)));
        let carry = vec![carry_s0, SampleColumns::empty()];
        let s0 = vec![record(12, ref_plus_alt(2, 3))];
        let s1 = vec![record(12, vec![ref_obs(10)])];
        let chunk = run_loader(0, 10..100, vec![s0, s1], carry);

        assert_eq!(chunk.per_sample[0].n_records(), 2);
        assert_eq!(chunk.per_sample[0].position_at(0), 5);
        assert_eq!(chunk.per_sample[0].position_at(1), 12);
        assert_eq!(chunk.per_sample[1].n_records(), 1);
        assert_eq!(chunk.per_sample[1].position_at(0), 12);
    }

    #[test]
    fn loader_clears_carryover_after_draining_it() {
        let mut carry_s0 = SampleColumns::empty();
        carry_s0.push_record(record(5, ref_plus_alt(3, 4)));
        let mut carry = vec![carry_s0, SampleColumns::empty()];
        let s0 = vec![record(12, ref_plus_alt(2, 3))];
        let s1 = vec![record(12, vec![ref_obs(10)])];

        let n_samples = 2;
        let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let mut out = MaterialisedChunk::with_n_samples(n_samples);
        let iters: Vec<_> = vec![s0, s1]
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>))
            .collect();
        load_chunk_from_iters(&mut scratch, &mut out, 0, 10..100, iters, &mut carry).unwrap();

        assert_eq!(carry[0].n_records(), 0);
        assert_eq!(carry[1].n_records(), 0);
    }

    #[test]
    fn loader_stops_pulling_at_range_end() {
        let s0 = vec![
            record(15, ref_plus_alt(2, 3)),
            record(60, ref_plus_alt(2, 3)),
            record(150, ref_plus_alt(2, 3)),
        ];
        let chunk = run_loader(0, 10..100, vec![s0], vec![SampleColumns::empty()]);
        assert_eq!(chunk.per_sample[0].n_records(), 2);
        assert_eq!(chunk.per_sample[0].position_at(0), 15);
        assert_eq!(chunk.per_sample[0].position_at(1), 60);
    }

    #[test]
    fn loader_rejects_inverted_range() {
        let n_samples = 1;
        let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let mut out = MaterialisedChunk::with_n_samples(n_samples);
        let mut carry = vec![SampleColumns::empty()];
        let iters: Vec<_> = vec![Vec::<PileupRecord>::new()]
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>))
            .collect();
        let err = load_chunk_from_iters(&mut scratch, &mut out, 0, 100..50, iters, &mut carry)
            .expect_err("inverted range rejected");
        assert!(matches!(err, ChunkLoadError::InvalidRange { .. }));
    }

    #[test]
    fn loader_rejects_record_from_unexpected_chromosome() {
        let n_samples = 1;
        let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let mut out = MaterialisedChunk::with_n_samples(n_samples);
        let mut carry = vec![SampleColumns::empty()];
        let bogus = PileupRecord::new(99, 50, ref_plus_alt(2, 3));
        let iters: Vec<_> = vec![vec![bogus]]
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>))
            .collect();
        let err = load_chunk_from_iters(&mut scratch, &mut out, 0, 10..100, iters, &mut carry)
            .expect_err("wrong-chromosome record rejected");
        assert!(matches!(
            err,
            ChunkLoadError::UnexpectedChromosome {
                expected_chrom_id: 0,
                got_chrom_id: 99,
                ..
            }
        ));
    }

    #[test]
    fn loader_surfaces_upstream_errors() {
        let n_samples = 1;
        let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let mut out = MaterialisedChunk::with_n_samples(n_samples);
        let mut carry = vec![SampleColumns::empty()];
        let iters: Vec<Box<dyn Iterator<Item = Result<PileupRecord, &'static str>>>> =
            vec![Box::new(std::iter::once(Err("psp blew up")))];
        let err = load_chunk_from_iters(&mut scratch, &mut out, 0, 10..100, iters, &mut carry)
            .expect_err("upstream error surfaced");
        assert!(matches!(
            err,
            ChunkLoadError::UpstreamRead {
                sample_idx: 0,
                source: "psp blew up"
            }
        ));
    }

    #[test]
    fn loader_clears_output_chunk_between_iterations() {
        // First chunk: one variant at pos 14. Second chunk
        // overwrites the same `out` buffer: range moves to a fresh
        // window with a different chrom_id and a different variant.
        let n_samples = 1;
        let mut scratch = ChunkLoadScratch::with_n_samples(n_samples);
        let mut out = MaterialisedChunk::with_n_samples(n_samples);
        let mut carry = vec![SampleColumns::empty()];

        let iters_a: Vec<_> = vec![vec![record(14, ref_plus_alt(2, 3))]]
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>))
            .collect();
        load_chunk_from_iters(&mut scratch, &mut out, 0, 10..50, iters_a, &mut carry).unwrap();
        assert_eq!(out.chrom_id, 0);
        assert_eq!(out.range, 10..50);
        assert_eq!(out.per_sample[0].n_records(), 1);

        let iters_b: Vec<_> = vec![vec![PileupRecord::new(7, 220, ref_plus_alt(3, 7))]]
            .into_iter()
            .map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>))
            .collect();
        load_chunk_from_iters(&mut scratch, &mut out, 7, 200..300, iters_b, &mut carry).unwrap();
        assert_eq!(out.chrom_id, 7);
        assert_eq!(out.range, 200..300);
        assert_eq!(out.per_sample[0].n_records(), 1);
        assert_eq!(out.per_sample[0].position_at(0), 220);
    }
}
