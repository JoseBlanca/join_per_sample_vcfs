//! Writer — section 3 (appendix §E).
//!
//! *(today: `var_calling::driver` `drive_blocks_parallel`'s receive loop)*
//!
//! `VcfWriter` consumes `CalledChunk`s and writes VCF, **single consumer**,
//! reordering by `chunk_order` via a `BTreeMap` reorder buffer drained on the
//! *next expected* `chunk_order`.
//!
//! **Gapless invariant:** every chunk yields exactly one `CalledChunk` (empty
//! ones included), or the drain stalls; a debug-assert checks the consumed
//! `chunk_order` sequence is contiguous.
//!
//! Phase 4 builds this module.

use std::collections::BTreeMap;

use crate::vcf::{CohortMetadata, CohortVcfWriter, VcfWriteError, WriterConfig};

use crate::var_calling::types::{CalledChunk, Variant};

/// Post-EM downstream filters applied per record at write time (the
/// `emit_or_drop` decision). Pulled out of `VarCallingArgs` /
/// `ChunkDriverParams.downstream`; the order and predicates are copied verbatim
/// from `driver::emit_or_drop` — they decide *which* records reach the VCF, so
/// they are part of the byte-identity contract.
#[derive(Debug, Clone, Copy)]
pub struct DownstreamFilters {
    pub min_qual_phred: f64,
    pub no_mapq_diff_filter: bool,
    pub min_mapq_diff_t: f32,
}

/// Run-level counters for the run summary (≈ the old `ChunkDriverStats`). The
/// emit-side counters are decided here; the caller-side counters
/// (`lh_cap_*` / `groups_skipped_*` / `records_dropped_low_alt_obs`) are rolled
/// from each [`CalledChunk`]'s [`CallStats`](crate::var_calling::types::CallStats).
/// Not part of the VCF.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct WriterStats {
    pub records_written: u64,
    pub records_dropped_hom_ref: u64,
    pub records_dropped_low_qual: u64,
    pub records_dropped_low_mapq_diff_t: u64,
    pub records_unconverged: u64,
    // Rolled from the callers' per-chunk CallStats.
    pub records_dropped_low_alt_obs: u64,
    pub lh_cap_groups_skipped: u64,
    pub lh_cap_alleles_in_skipped: u64,
    pub groups_skipped_post_unify_ref_only: u64,
}

/// Errors the writer can surface.
#[derive(Debug, thiserror::Error)]
pub enum WriterError {
    #[error("VCF write failed: {0}")]
    Vcf(#[from] VcfWriteError),
    /// `finish` was called with chunks still buffered: a gap in `chunk_order`
    /// stalled the reorder drain (a caller dropped a chunk instead of shipping
    /// exactly one `CalledChunk` per chunk), which would silently truncate the
    /// VCF. Release-level guard for the gapless invariant (was a
    /// `debug_assert!`).
    #[error("{count} chunk(s) never emitted — a gap in chunk_order stalled the writer drain")]
    MissingChunks { count: usize },
}

/// Section 3 — the VCF writer (appendix §E). Single consumer of the
/// caller→writer [`CalledChunk`] stream; reorders by `chunk_order` (a
/// `BTreeMap` drained on the next-expected order, so the emitted VCF is in
/// genomic order regardless of how the parallel callers finish) and applies the
/// downstream `emit_or_drop` filters before handing each surviving [`Variant`]
/// to the reused [`CohortVcfWriter`] (the byte-identity-critical formatting is
/// the shared `crate::vcf` module, untouched).
pub struct VcfWriter {
    inner: CohortVcfWriter,
    filters: DownstreamFilters,
    stats: WriterStats,
    /// Reorder buffer: chunks that arrived before their turn.
    reorder: BTreeMap<u64, CalledChunk>,
    /// The next `chunk_order` to emit (the gapless cursor).
    next_expected: u64,
}

impl VcfWriter {
    /// Open the underlying [`CohortVcfWriter`] (header written here) and start
    /// the reorder buffer at `chunk_order` 0.
    pub fn new(
        metadata: CohortMetadata,
        config: WriterConfig,
        filters: DownstreamFilters,
    ) -> Result<Self, WriterError> {
        Ok(Self {
            inner: CohortVcfWriter::new(metadata, config)?,
            filters,
            stats: WriterStats::default(),
            reorder: BTreeMap::new(),
            next_expected: 0,
        })
    }

    /// Accept one [`CalledChunk`] (in any order) and emit every chunk that is
    /// now contiguously available from `next_expected`.
    pub fn handle(&mut self, chunk: CalledChunk) -> Result<(), WriterError> {
        debug_assert!(
            chunk.chunk_order >= self.next_expected,
            "chunk_order {} already emitted (next_expected {})",
            chunk.chunk_order,
            self.next_expected,
        );
        self.reorder.insert(chunk.chunk_order, chunk);
        while let Some(ready) = self.reorder.remove(&self.next_expected) {
            self.emit_chunk(ready)?;
            self.next_expected += 1;
        }
        Ok(())
    }

    /// Finish writing: every produced chunk must have been consumed (the
    /// gapless invariant — one `CalledChunk` per chunk, empty included), so the
    /// reorder buffer must be empty.
    pub fn finish(self) -> Result<WriterStats, WriterError> {
        if !self.reorder.is_empty() {
            return Err(WriterError::MissingChunks {
                count: self.reorder.len(),
            });
        }
        self.inner.finish()?;
        Ok(self.stats)
    }

    fn emit_chunk(&mut self, chunk: CalledChunk) -> Result<(), WriterError> {
        // Roll the caller-side counters into the run summary.
        self.stats.records_dropped_low_alt_obs += chunk.stats.records_dropped_low_alt_obs;
        self.stats.lh_cap_groups_skipped += chunk.stats.lh_cap_groups_skipped;
        self.stats.lh_cap_alleles_in_skipped += chunk.stats.lh_cap_alleles_in_skipped;
        self.stats.groups_skipped_post_unify_ref_only +=
            chunk.stats.groups_skipped_post_unify_ref_only;
        for record in chunk.records {
            self.emit_or_drop(record)?;
        }
        Ok(())
    }

    /// Per-record filter + write, copied verbatim from `driver::emit_or_drop`
    /// (the `min_alt_obs_per_sample` filter is upstream, in the caller). Order
    /// is load-bearing: hom-ref, then QUAL, then the MAPQ-diff t-test.
    fn emit_or_drop(&mut self, record: Variant) -> Result<(), WriterError> {
        if !record.is_variant_call() {
            self.stats.records_dropped_hom_ref += 1;
            return Ok(());
        }
        if record.qual_phred < self.filters.min_qual_phred {
            self.stats.records_dropped_low_qual += 1;
            return Ok(());
        }
        if !self.filters.no_mapq_diff_filter
            && record_fails_mapq_diff_t(&record, self.filters.min_mapq_diff_t)
        {
            self.stats.records_dropped_low_mapq_diff_t += 1;
            return Ok(());
        }
        if !record.diagnostics.converged {
            self.stats.records_unconverged += 1;
        }
        self.inner.write_record(&record)?;
        self.stats.records_written += 1;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// MAPQ-difference Welch's-t filter (copied verbatim from `driver`).
// ---------------------------------------------------------------------------

const MAPQ_FILTER_MIN_READS_PER_SIDE: u64 = 3;

/// Test alias for the private MAPQ Welch's-t filter
/// ([`record_fails_mapq_diff_t`]) so the integration suite can exercise the
/// decision in isolation (was `driver::record_fails_mapq_diff_t_for_test`).
#[doc(hidden)]
pub fn record_fails_mapq_diff_t_for_test(record: &Variant, threshold: f32) -> bool {
    record_fails_mapq_diff_t(record, threshold)
}

/// Per-allele MAPQ moments pooled across the cohort.
struct PooledMapqMoments {
    n: u64,
    sum: u64,
    sum_of_squares: u128,
}

fn pool_allele_mapq(record: &Variant, allele_idx: usize) -> PooledMapqMoments {
    let mut moments = PooledMapqMoments {
        n: 0,
        sum: 0,
        sum_of_squares: 0,
    };
    for sample_idx in 0..record.n_samples {
        let stats = &record.scalars_row(sample_idx)[allele_idx];
        moments.n += u64::from(stats.num_obs);
        moments.sum += u64::from(stats.mapq_sum);
        moments.sum_of_squares += stats.mapq_sum_sq as u128;
    }
    moments
}

/// Welch's-t MAPQ-difference filter applied post-EM in `emit_or_drop`.
fn record_fails_mapq_diff_t(record: &Variant, threshold: f32) -> bool {
    if !threshold.is_finite() {
        return false;
    }
    let n_alleles = record.alleles.len();
    if n_alleles < 2 {
        return false;
    }
    let ref_moments = pool_allele_mapq(record, 0);
    if ref_moments.n < MAPQ_FILTER_MIN_READS_PER_SIDE {
        return false;
    }
    let mean_ref = ref_moments.sum as f64 / ref_moments.n as f64;
    let var_ref = ((ref_moments.sum_of_squares as f64 - (ref_moments.sum as f64) * mean_ref)
        .max(0.0))
        / ((ref_moments.n - 1) as f64);
    for alt_idx in 1..n_alleles {
        let alt_moments = pool_allele_mapq(record, alt_idx);
        if alt_moments.n < MAPQ_FILTER_MIN_READS_PER_SIDE {
            continue;
        }
        let mean_alt = alt_moments.sum as f64 / alt_moments.n as f64;
        let var_alt = ((alt_moments.sum_of_squares as f64 - (alt_moments.sum as f64) * mean_alt)
            .max(0.0))
            / ((alt_moments.n - 1) as f64);
        let se2 = var_alt / (alt_moments.n as f64) + var_ref / (ref_moments.n as f64);
        if se2 <= 0.0 {
            continue;
        }
        let t = (mean_alt - mean_ref) / se2.sqrt();
        if (t as f32) < threshold {
            return true;
        }
    }
    false
}
