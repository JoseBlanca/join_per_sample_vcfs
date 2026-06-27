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

/// The allele-balance filter setting: either off, or on with a threshold and
/// Beta-Binomial concentration. A single value so the "off" state cannot carry
/// a stale, silently-ignored threshold.
#[derive(Debug, Clone, Copy)]
pub enum AlleleBalanceFilter {
    /// `--no-allele-balance-filter`: the test is never run.
    Off,
    /// Drop a biallelic het call whose allele-balance log-LR falls below
    /// `min_log_lr`, evaluated against a Beta-Binomial of the given
    /// `concentration` (see [`crate::var_calling::allele_balance`]).
    On { min_log_lr: f64, concentration: f64 },
}

/// Post-EM downstream filters applied per record at write time (the
/// `emit_or_drop` decision). They decide *which* records reach the VCF, so
/// they are part of the byte-identity contract.
#[derive(Debug, Clone, Copy)]
pub struct DownstreamFilters {
    pub min_qual_phred: f64,
    pub allele_balance: AlleleBalanceFilter,
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
    pub records_dropped_allele_balance: u64,
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
    /// the allele-balance test.
    fn emit_or_drop(&mut self, record: Variant) -> Result<(), WriterError> {
        if !record.is_variant_call() {
            self.stats.records_dropped_hom_ref += 1;
            return Ok(());
        }
        // Gate on the FINAL QUAL we are about to write (engine baseline
        // refined for artifact shape, then clamped) — NOT the pre-refinement
        // baseline. The refinement deflates depth-inflation artifacts (steady
        // low-VAF alt support) to ~0; gating on the baseline let those through
        // as PASS records with a written QUAL of 0. Computing it here, via the
        // same path the writer uses, keeps the filter and the QUAL column in
        // lock-step.
        if f64::from(self.inner.final_qual(&record)) < self.filters.min_qual_phred {
            self.stats.records_dropped_low_qual += 1;
            return Ok(());
        }
        if let AlleleBalanceFilter::On {
            min_log_lr,
            concentration,
        } = self.filters.allele_balance
            && record_fails_allele_balance(&record, min_log_lr, concentration)
        {
            self.stats.records_dropped_allele_balance += 1;
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
// Allele-balance filter (post-EM, in `emit_or_drop`).
// ---------------------------------------------------------------------------

/// Test alias for the private allele-balance filter
/// ([`record_fails_allele_balance`]) so the integration suite can exercise the
/// decision in isolation.
#[doc(hidden)]
pub fn record_fails_allele_balance_for_test(
    record: &Variant,
    min_log_lr: f64,
    concentration: f64,
) -> bool {
    record_fails_allele_balance(record, min_log_lr, concentration)
}

/// Allele-balance drop test. Restricted to **biallelic** records (the SNP
/// artefact class). For each variant-carrying sample:
///
/// - a hom-alt carrier (or any non-het variant genotype) is *untested support*
///   — allele balance can't see hom-alt artefacts — and keeps the site;
/// - a het carrier whose observed (ref, alt) read split fits its expected
///   balance (allele-balance log-LR `≥ min_log_lr`) keeps the site;
/// - only if there is at least one het carrier and **every** variant carrier is
///   a het that fails the test is the record dropped.
///
/// So a clean carrier anywhere rescues the site; the gate fires exactly when
/// all the variant evidence is allele-balance-failing hets.
fn record_fails_allele_balance(record: &Variant, min_log_lr: f64, concentration: f64) -> bool {
    use crate::var_calling::allele_balance::{allele_balance_log_lr, expected_vaf};

    // Biallelic only: REF + one ALT. For this shape `best_genotype[s]` is the
    // sample's alt-copy count (the genotype enumeration orders biallelic
    // genotypes by ascending alt count), so 0 = hom-ref, `ploidy` = hom-alt,
    // anything between = het.
    if record.alleles.len() != 2 {
        return false;
    }
    // SNP only (v1). Indel allele balance from the AD counts is unreliable —
    // reads spanning an indel are assigned ambiguously, so real indel hets show
    // skewed fractions and would be wrongly dropped (measured ~29% indel-recall
    // loss at 300×). A SNP has a single-base REF and ALT over the group span.
    if record.alleles[0].seq.len() != 1 || record.alleles[1].seq.len() != 1 {
        return false;
    }
    let ploidy = record.ploidy;
    let mut saw_failing_het = false;
    for sample_idx in 0..record.n_samples {
        let alt_copies = record.best_genotype[sample_idx];
        if alt_copies == 0 {
            continue; // hom-ref: not a carrier
        }
        if alt_copies as u32 >= u32::from(ploidy) {
            return false; // hom-alt support — allele balance can't challenge it
        }
        let row = record.scalars_row(sample_idx);
        let lr = allele_balance_log_lr(
            row[0].num_obs,
            row[1].num_obs,
            expected_vaf(alt_copies as u8, ploidy),
            concentration,
        );
        if lr >= min_log_lr {
            return false; // a clean het — the variant is real
        }
        saw_failing_het = true;
    }
    saw_failing_het
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

    fn filters_no_ab() -> DownstreamFilters {
        DownstreamFilters {
            min_qual_phred: 30.0,
            allele_balance: AlleleBalanceFilter::Off,
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
            filters_no_ab(),
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
            filters_no_ab(),
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
        // covered in isolation by `record_fails_allele_balance_for_test`):
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let mut writer = VcfWriter::new(
            metadata_two_samples(),
            WriterConfig::new(out),
            filters_no_ab(),
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
