//! Writer — section 3 (appendix §E).
//!
//! `VcfWriter` consumes `CalledChunk`s and writes VCF, **single consumer**,
//! reordering by `chunk_order` via a `BTreeMap` reorder buffer drained on the
//! *next expected* `chunk_order`.
//!
//! **Gapless invariant:** every chunk yields exactly one `CalledChunk` (empty
//! ones included), or the drain stalls; [`finish`](VcfWriter::finish) surfaces a
//! leftover-buffer gap as [`WriterError::MissingChunks`] rather than silently
//! truncating the VCF.

use std::collections::BTreeMap;

use crate::vcf::{CohortMetadata, CohortVcfWriter, VcfWriteError, WriterConfig};

use crate::var_calling::types::{CallStats, CalledChunk, Variant};

/// The MAPQ-difference Welch's-t filter setting: either off, or on with a
/// threshold. A single value so the "off" state cannot carry a stale,
/// silently-ignored threshold (the two old co-dependent fields collapse here).
#[derive(Debug, Clone, Copy)]
pub enum MapqDiffFilter {
    /// `--no-mapq-diff-filter`: the test is never run.
    Off,
    /// Drop a record whose MAPQ-difference Welch's-t falls below `min_t`.
    On { min_t: f32 },
}

/// Post-EM downstream filters applied per record at write time (the
/// `emit_or_drop` decision). They decide *which* records reach the VCF, so
/// they are part of the byte-identity contract.
#[derive(Debug, Clone, Copy)]
pub struct DownstreamFilters {
    pub min_qual_phred: f64,
    pub mapq_diff: MapqDiffFilter,
}

/// Run-level counters for the run summary (≈ the old `ChunkDriverStats`). The
/// emit-side counters are decided here; the caller-side counters
/// (`lh_cap_*` / `groups_skipped_*` / `records_dropped_low_alt_obs`) are rolled
/// from each [`CalledChunk`]'s [`CallStats`].
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
        let CalledChunk {
            chunk_order: _,
            records,
            stats,
        } = chunk;
        // Roll the caller-side counters into the run summary. Destructured
        // exhaustively so a new `CallStats` counter is a compile error here
        // (rather than being silently dropped from the run summary).
        let CallStats {
            lh_cap_groups_skipped,
            lh_cap_alleles_in_skipped,
            groups_skipped_post_unify_ref_only,
            records_dropped_low_alt_obs,
        } = stats;
        self.stats.records_dropped_low_alt_obs += records_dropped_low_alt_obs;
        self.stats.lh_cap_groups_skipped += lh_cap_groups_skipped;
        self.stats.lh_cap_alleles_in_skipped += lh_cap_alleles_in_skipped;
        self.stats.groups_skipped_post_unify_ref_only += groups_skipped_post_unify_ref_only;
        for record in records {
            self.emit_or_drop(record)?;
        }
        Ok(())
    }

    /// Per-record filter + write (the `min_alt_obs_per_sample` filter is
    /// upstream, in the caller). Order is load-bearing: hom-ref, then QUAL, then
    /// the MAPQ-diff t-test — byte-identical to the pre-rewrite emit gate
    /// (verified out-of-tree).
    fn emit_or_drop(&mut self, record: Variant) -> Result<(), WriterError> {
        if !record.is_variant_call() {
            self.stats.records_dropped_hom_ref += 1;
            return Ok(());
        }
        if record.qual_phred < self.filters.min_qual_phred {
            self.stats.records_dropped_low_qual += 1;
            return Ok(());
        }
        if let MapqDiffFilter::On { min_t } = self.filters.mapq_diff
            && record_fails_mapq_diff_t(&record, min_t)
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
// MAPQ-difference Welch's-t filter (byte-identical to the pre-rewrite filter).
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup_record::AlleleSupportStats;
    use crate::psp::header::ParsedChromosome;
    use crate::var_calling::per_group_merger::MergedAllele;
    use crate::var_calling::posterior_engine::{EmDiagnostics, PosteriorRecord, RecordLocus};
    use tempfile::tempdir;

    fn allele(seq: &[u8]) -> MergedAllele {
        MergedAllele {
            seq: seq.to_vec(),
            is_compound: false,
            constituents: Vec::new(),
        }
    }

    fn support(num_obs: u32) -> AlleleSupportStats {
        AlleleSupportStats::new(num_obs, 0.0, 0, 0, 0, 0, 0)
    }

    fn metadata_two_samples() -> CohortMetadata {
        CohortMetadata {
            sample_names: vec!["S0".into(), "S1".into()],
            contigs: vec![ParsedChromosome {
                name: "chr1".into(),
                length: 1_000_000,
                md5: "0123456789abcdef0123456789abcdef".into(),
            }],
            tool_string: "pop_var_caller vcf-writer-test".into(),
            command_line: "pop_var_caller cohort --output out.vcf".into(),
        }
    }

    /// A biallelic A→T SNP at `pos` for a 2-sample cohort. `best` is each
    /// sample's argmax genotype index (0 = hom-ref); a non-zero entry makes it a
    /// variant call. Scalars are inert (MAPQ filtering is disabled in these
    /// tests, so the per-allele moments are unused).
    fn snp_variant(pos: u32, qual: f64, best: [usize; 2], converged: bool) -> Variant {
        PosteriorRecord {
            locus: RecordLocus {
                chrom_id: 0,
                start: pos,
                end: pos,
            },
            alleles: vec![allele(b"A"), allele(b"T")],
            ploidy: 2,
            n_samples: 2,
            n_genotypes: 3,
            allele_frequencies: vec![0.5, 0.5],
            compound_frequencies: vec![None, None],
            posteriors: vec![0.0; 6],
            best_genotype: best.to_vec(),
            gq_phred: vec![50.0, 50.0],
            qual_phred: qual,
            scalars: vec![support(20), support(0), support(20), support(0)],
            other_scalars: vec![],
            chain_anchor_flags: vec![false; 4],
            diagnostics: EmDiagnostics {
                iterations: 1,
                final_max_delta_p: 1e-6,
                converged,
            },
        }
    }

    fn called(chunk_order: u64, records: Vec<Variant>) -> CalledChunk {
        CalledChunk {
            chunk_order,
            records,
            stats: CallStats::default(),
        }
    }

    fn filters_no_mapq() -> DownstreamFilters {
        DownstreamFilters {
            min_qual_phred: 30.0,
            mapq_diff: MapqDiffFilter::Off,
        }
    }

    #[test]
    fn emits_in_chunk_order_when_chunks_arrive_out_of_order() {
        // chunk_order is genomic order (0 ⇒ pos 100, 1 ⇒ 200, 2 ⇒ 300);
        // delivered out of order, the BTreeMap reorder must still emit 100, 200,
        // 300 — the byte-identical-for-any-worker-count mechanism.
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let mut writer = VcfWriter::new(
            metadata_two_samples(),
            WriterConfig::new(out.clone()),
            filters_no_mapq(),
        )
        .unwrap();
        writer
            .handle(called(2, vec![snp_variant(300, 100.0, [1, 0], true)]))
            .unwrap();
        writer
            .handle(called(0, vec![snp_variant(100, 100.0, [1, 0], true)]))
            .unwrap();
        writer
            .handle(called(1, vec![snp_variant(200, 100.0, [1, 0], true)]))
            .unwrap();
        let stats = writer.finish().unwrap();
        assert_eq!(stats.records_written, 3);

        let text = std::fs::read_to_string(&out).unwrap();
        let positions: Vec<&str> = text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').nth(1).unwrap())
            .collect();
        assert_eq!(
            positions,
            ["100", "200", "300"],
            "emitted in chunk_order, not arrival order"
        );
    }

    #[test]
    fn finish_errors_with_missing_chunks_on_gap() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let mut writer = VcfWriter::new(
            metadata_two_samples(),
            WriterConfig::new(out),
            filters_no_mapq(),
        )
        .unwrap();
        writer
            .handle(called(0, vec![snp_variant(100, 100.0, [1, 0], true)]))
            .unwrap();
        // chunk_order 1 is never delivered; chunk 2 stays buffered behind the gap.
        writer
            .handle(called(2, vec![snp_variant(300, 100.0, [1, 0], true)]))
            .unwrap();
        match writer.finish() {
            Err(WriterError::MissingChunks { count }) => assert_eq!(count, 1),
            other => panic!("expected MissingChunks, got {other:?}"),
        }
    }

    #[test]
    fn emit_or_drop_counts_each_bucket() {
        // One chunk carrying one record per outcome (the MAPQ-diff gate is
        // covered in isolation by `record_fails_mapq_diff_t_for_test`):
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let mut writer = VcfWriter::new(
            metadata_two_samples(),
            WriterConfig::new(out),
            filters_no_mapq(),
        )
        .unwrap();
        let records = vec![
            snp_variant(100, 100.0, [0, 0], true), // hom-ref ⇒ dropped_hom_ref
            snp_variant(200, 5.0, [1, 0], true),   // QUAL 5 < 30 ⇒ dropped_low_qual
            snp_variant(300, 100.0, [1, 0], true), // written
            snp_variant(400, 100.0, [0, 1], false), // written + unconverged
        ];
        writer.handle(called(0, records)).unwrap();
        let stats = writer.finish().unwrap();
        assert_eq!(stats.records_written, 2);
        assert_eq!(stats.records_dropped_hom_ref, 1);
        assert_eq!(stats.records_dropped_low_qual, 1);
        assert_eq!(stats.records_unconverged, 1);
    }
}
