//! [`CohortVcfWriter`] — the runtime that ties sink, header, and
//! record-encoder together.
//!
//! `new` opens `<output>.tmp`, builds the header, writes it, and
//! parks the encoder; `write_record` encodes one `PosteriorRecord`
//! into a noodles `RecordBuf` and feeds it to the underlying VCF
//! writer; `finish` flushes the sink (including the bgzf EOF block
//! for gzipped paths), `fsync`s the file and the parent directory,
//! then atomic-renames into place.
//!
//! Forgetting `finish` leaves `<output>.tmp` on disk and no
//! `<output>` — intended loud failure. The `#[must_use]` attribute
//! on [`CohortVcfWriter`] gives the compiler a chance to catch the
//! forgotten-finish case at the construction site.

use std::collections::HashMap;
use std::path::PathBuf;

use noodles_vcf::Header as VcfHeader;
use noodles_vcf::variant::io::Write as _;
use noodles_vcf::variant::record_buf::samples::Keys;

use super::WriterConfig;
use super::errors::VcfWriteError;
use super::header::{CohortMetadata, build_vcf_header};
use super::record_encode::{build_format_keys, encode};
use super::sink::{SinkKind, tmp_path_for};
use super::writable::VcfWritable;
use crate::psp::header::ParsedChromosome;
use crate::var_calling::per_group_merger::genotype_order;

/// Streams `PosteriorRecord` items to a VCF file.
///
/// Construction opens `<output>.tmp`, builds the noodles header, and
/// writes it immediately so a subsequent `write_record` does not pay
/// header-build cost. Per-record encoding goes through
/// [`super::record_encode::encode`]. The writer's
/// non-decreasing-locus order check is the last line of defence
/// against an upstream that surfaces records out of order.
///
/// **The writer must be finalised with [`CohortVcfWriter::finish`].**
/// Otherwise the output is left at `<output>.tmp` and no `<output>`
/// is produced; the `#[must_use]` attribute warns at the construction
/// site if the value is dropped without `finish`.
#[must_use = "CohortVcfWriter must be finalised by calling `.finish()`; \
              otherwise the output is left at `<output>.tmp` and no \
              `<output>` is produced"]
pub struct CohortVcfWriter {
    /// noodles' line-oriented VCF writer wraps our sink directly. The
    /// sink (plain or bgzf) is owned by the writer and recovered on
    /// `finish`.
    inner: noodles_vcf::io::Writer<SinkKind>,
    /// Built once at construction; passed to noodles per record because
    /// `write_variant_record` takes `&Header`.
    header: VcfHeader,
    /// Cohort metadata frozen at construction. Only the `contigs`
    /// slice is read on the per-record path (to map `chrom_id` →
    /// CHROM string); the sample-name vector is consumed by the
    /// header build at construction and not re-read.
    contigs: Vec<ParsedChromosome>,
    /// Cohort size, re-checked against every record's `n_samples`.
    expected_samples: usize,
    /// Final path; tmp is `<final_path>.tmp`. Held until `finish`.
    final_path: PathBuf,
    /// FORMAT keys list — built once from `config.emit_gp` and reused
    /// per record.
    format_keys: Keys,
    /// Whole-run config — held by value for the per-record encoder and
    /// surfaced via [`CohortVcfWriter::config`] for runtime inspection.
    config: WriterConfig,
    /// Latches the most recent `(chrom_id, start)` so the next
    /// record can be checked for non-decreasing order.
    last_locus: Option<(u32, u32)>,
    /// Mi17: cache of `genotype_order(ploidy, n_alleles)` tables keyed
    /// by the pair. Per-record encoding looks up here instead of
    /// rebuilding the table; in steady state (cohort with constant
    /// ploidy and a small range of `n_alleles`) every record beyond
    /// the first hits the cache.
    genotype_tables: HashMap<(u8, usize), Vec<Vec<u8>>>,
}

impl std::fmt::Debug for CohortVcfWriter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Exhaustive destructure: a new field on the struct fails to
        // compile here so the omission from Debug output is explicit
        // rather than silent.
        let Self {
            inner: _,
            header: _,
            contigs,
            expected_samples,
            final_path,
            format_keys,
            config,
            last_locus,
            genotype_tables,
        } = self;
        f.debug_struct("CohortVcfWriter")
            .field("config", config)
            .field("expected_samples", expected_samples)
            .field("contigs_len", &contigs.len())
            .field("format_keys", format_keys)
            .field("final_path", final_path)
            .field("last_locus", last_locus)
            .field("cached_genotype_tables", &genotype_tables.len())
            .finish_non_exhaustive()
    }
}

impl CohortVcfWriter {
    /// Open `<config.output>.tmp`, build the VCF header, and write
    /// it. If the header write fails the tmp file is removed before
    /// the error bubbles, so a constructor failure is not
    /// indistinguishable from "user forgot to call finish".
    ///
    /// # Errors
    ///
    /// * [`VcfWriteError::InvalidMetadata`] — cohort metadata fails
    ///   validation (empty sample names, duplicate sample/contig
    ///   name).
    /// * [`VcfWriteError::ContigLengthOverflow`] — a contig length
    ///   exceeds `i32::MAX`.
    /// * [`VcfWriteError::Encode`] — noodles refuses a header value
    ///   (`##source` / `##commandline` parse or insert).
    /// * [`VcfWriteError::CreateTmp`] — `File::create` on
    ///   `<output>.tmp` failed.
    /// * [`VcfWriteError::WriteHeader`] — writing the header bytes
    ///   to the sink failed; the tmp file is removed on this path.
    pub fn new(metadata: CohortMetadata, config: WriterConfig) -> Result<Self, VcfWriteError> {
        let header = build_vcf_header(&metadata, &config)?;
        let format_keys = build_format_keys(&config);
        let expected_samples = metadata.sample_names.len();
        // Mi9: move the field out of `metadata` instead of cloning.
        let contigs = metadata.contigs;
        let final_path = config.output.clone();

        let sink = SinkKind::open_tmp(&final_path)?;
        let mut inner = noodles_vcf::io::Writer::new(sink);
        // M13: a header-write failure must remove the tmp file so a
        // constructor crash is not confused with "user forgot to
        // call finish".
        if let Err(write_err) = inner.write_header(&header) {
            let tmp = tmp_path_for(&final_path);
            let _ = std::fs::remove_file(&tmp); // best-effort cleanup
            return Err(VcfWriteError::WriteHeader {
                tmp_path: tmp,
                source: write_err,
            });
        }
        Ok(Self {
            inner,
            header,
            contigs,
            expected_samples,
            final_path,
            format_keys,
            config,
            last_locus: None,
            genotype_tables: HashMap::new(),
        })
    }

    /// The config this writer was constructed with (frozen at `new`).
    /// Useful for runtime inspection ("is GP enabled?", "what path?")
    /// when the caller no longer holds the original `WriterConfig`.
    #[must_use]
    pub fn config(&self) -> &WriterConfig {
        &self.config
    }

    /// Encode and write one record.
    ///
    /// Latches the locus order; an upstream regression surfaces as
    /// [`VcfWriteError::RecordOutOfOrder`] and the writer keeps
    /// running (the next call gets a fresh order check against the
    /// previous accepted record). `last_locus` is *not* updated on
    /// any error path, so a rejected record does not corrupt the
    /// order-check baseline.
    ///
    /// # Errors
    ///
    /// * [`VcfWriteError::RecordOutOfOrder`] — `record.locus` does
    ///   not strictly exceed the previous accepted record's locus.
    /// * [`VcfWriteError::SampleCountMismatch`] — `record.n_samples`
    ///   differs from the cohort metadata's sample count.
    /// * [`VcfWriteError::InconsistentRecord`] — a per-record
    ///   vector length doesn't match the record's declared shape.
    /// * [`VcfWriteError::UnknownChromId`] — `record.locus.chrom_id`
    ///   is out of bounds for the contig table.
    /// * [`VcfWriteError::GenotypeIndexOutOfBounds`] /
    ///   [`VcfWriteError::AlleleIndexOutOfBounds`] — a decoded
    ///   `best_genotype` cell references a non-existent slot.
    /// * [`VcfWriteError::DepthOverflow`] — a depth (per-sample DP,
    ///   per-allele AD, or cohort total) overflows `i32`.
    /// * [`VcfWriteError::Encode`] — noodles refused a value
    ///   (non-UTF-8 allele bytes, invalid `Position`, …).
    /// * [`VcfWriteError::WriteRecord`] — the sink rejected the
    ///   serialised bytes.
    ///
    /// The final (refined + clamped) QUAL this writer would emit for
    /// `record`, computed via the same path as [`write_record`] and reusing
    /// the cached genotype table. The filtering layer calls this so it gates
    /// `--min-qual` on the exact value that lands in the QUAL column, rather
    /// than the engine's pre-refinement baseline.
    ///
    /// [`write_record`]: Self::write_record
    pub fn final_qual<R: VcfWritable>(&mut self, record: &R) -> f32 {
        let table = self
            .genotype_tables
            .entry((record.ploidy(), record.n_alleles()))
            .or_insert_with(|| genotype_order(record.ploidy(), record.n_alleles()));
        super::record_encode::final_qual(record, table)
    }

    pub fn write_record<R: VcfWritable>(&mut self, record: &R) -> Result<(), VcfWriteError> {
        let locus = (record.chrom_id(), record.pos_1based());
        if let Some(prev) = self.last_locus
            && locus <= prev
        {
            return Err(VcfWriteError::RecordOutOfOrder {
                chrom_id: locus.0,
                pos: locus.1,
                prev_chrom_id: prev.0,
                prev_pos: prev.1,
            });
        }

        let table = self
            .genotype_tables
            .entry((record.ploidy(), record.n_alleles()))
            .or_insert_with(|| genotype_order(record.ploidy(), record.n_alleles()))
            .clone();
        let buf = encode(
            record,
            &self.contigs,
            &self.config,
            &self.format_keys,
            self.expected_samples,
            &table,
        )?;
        self.inner
            .write_variant_record(&self.header, &buf)
            .map_err(|source| VcfWriteError::WriteRecord {
                chrom_id: locus.0,
                pos: locus.1,
                source,
            })?;
        self.last_locus = Some(locus);
        Ok(())
    }

    /// Flush, emit the bgzf EOF block when applicable, sync the file
    /// and parent directory, then atomic-rename `<output>.tmp` →
    /// `<output>`.
    ///
    /// Consumes `self`, so a forgotten call leaves `<output>.tmp` on
    /// disk and no `<output>` — see the type-level `#[must_use]`
    /// note.
    ///
    /// # Errors
    ///
    /// * [`VcfWriteError::FinishBgzf`] — the bgzf sink failed to
    ///   emit its EOF block / flush its tail.
    /// * [`VcfWriteError::WriteHeader`] — the `BufWriter` tail flush
    ///   failed on the plain-text path (the operation is "writing
    ///   the buffered tail").
    /// * [`VcfWriteError::FsyncFile`] — `File::sync_all` on the tmp
    ///   file failed.
    /// * [`VcfWriteError::Rename`] — `fs::rename` from the tmp path
    ///   to the final path failed.
    /// * [`VcfWriteError::FsyncDir`] — opening or fsyncing the
    ///   parent directory failed (without this fsync the rename is
    ///   not durable across a crash).
    pub fn finish(self) -> Result<(), VcfWriteError> {
        let sink = self.inner.into_inner();
        sink.finish(&self.final_path)?;
        Ok(())
    }

    /// Abandon the in-progress run: drop the sink (which closes the
    /// tmp file handle) and delete the tmp file at the exact path
    /// this writer used (`tmp_path_for(&self.final_path)`).
    ///
    /// Pair with [`finish`](Self::finish): a driver that returns
    /// early on error should call `abort()` instead of leaving the
    /// `<output>.tmp` file behind. Both consume `self`, so the type
    /// system forces a single terminal call.
    ///
    /// # Errors
    ///
    /// Returns `Err(io::Error)` if `std::fs::remove_file` failed
    /// (typically `NotFound` if the tmp file was never created, or a
    /// permission/IO failure). Callers should generally log and
    /// continue — the original driver error is the operator-facing
    /// one to surface.
    pub fn abort(self) -> std::io::Result<()> {
        let tmp = tmp_path_for(&self.final_path);
        drop(self.inner); // close the file handle
        std::fs::remove_file(tmp)
    }
}

#[cfg(test)]
mod tests {
    use tempfile::tempdir;

    use super::*;
    use crate::pileup_record::AlleleSupportStats;
    use crate::var_calling::per_group_merger::MergedAllele;
    use crate::var_calling::posterior_engine::{EmDiagnostics, PosteriorRecord, RecordLocus};

    fn ref_allele(seq: &[u8]) -> MergedAllele {
        MergedAllele {
            seq: seq.to_vec(),
            is_compound: false,
            constituents: Vec::new(),
        }
    }
    fn alt_allele(seq: &[u8]) -> MergedAllele {
        ref_allele(seq)
    }
    fn support(num_obs: u32) -> AlleleSupportStats {
        AlleleSupportStats {
            num_obs,
            q_sum: 0.0,
            fwd: 0,
            placed_left: 0,
            placed_start: 0,

            mapq_sum: 0,
            mapq_sum_sq: 0,
        }
    }
    fn fixture_metadata() -> CohortMetadata {
        CohortMetadata {
            sample_names: vec!["S0".into(), "S1".into()],
            contigs: vec![ParsedChromosome {
                name: "chr1".into(),
                length: 1_000_000,
                md5: "00000000000000000000000000000001".into(),
            }],
            tool_string: "pop_var_caller test".into(),
            command_line: String::new(),
            paralog_provenance: String::new(),
        }
    }
    fn mk_record(start: u32) -> PosteriorRecord {
        PosteriorRecord {
            locus: RecordLocus {
                chrom_id: 0,
                start,
                end: start,
            },
            alleles: vec![ref_allele(b"A"), alt_allele(b"T")],
            ploidy: 2,
            n_samples: 2,
            n_genotypes: 3,
            allele_frequencies: vec![0.75, 0.25],
            compound_frequencies: vec![None, None],
            posteriors: vec![0.98, 0.01, 0.01, 0.05, 0.90, 0.05],
            best_genotype: vec![0, 1],
            gq_phred: vec![60.0, 40.0],
            qual_phred: 100.0,
            scalars: vec![support(20), support(0), support(10), support(10)],
            other_scalars: vec![],
            chain_anchor_flags: vec![false; 4],
            diagnostics: EmDiagnostics {
                iterations: 5,
                final_max_delta_p: 1e-6,
                converged: true,
            },
            paralog_posterior: None,
        }
    }

    fn cfg_for(out: PathBuf) -> WriterConfig {
        WriterConfig::new(out)
    }

    #[test]
    fn full_round_trip_plain_text() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let metadata = fixture_metadata();

        let mut writer = CohortVcfWriter::new(metadata, cfg_for(out.clone())).unwrap();
        writer.write_record(&mk_record(100)).unwrap();
        writer.write_record(&mk_record(200)).unwrap();
        writer.finish().unwrap();

        let text = std::fs::read_to_string(&out).unwrap();
        assert!(text.contains("##fileformat=VCFv4.4"));
        let data_lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 2);
        assert!(data_lines[0].starts_with("chr1\t100\t"));
        assert!(data_lines[1].starts_with("chr1\t200\t"));
        assert!(!out.with_extension("vcf.tmp").exists());
    }

    /// B1: `abort()` removes the writer's tmp file (the exact path it
    /// used, derived from `tmp_path_for(&final_path)`) and leaves no
    /// final output behind. Compare against `finish()`'s success-path
    /// guarantee (the previous test asserts the tmp does not exist
    /// after `finish` because the rename consumed it).
    #[test]
    fn abort_removes_tmp_file_and_leaves_no_final_output() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let tmp = tmp_path_for(&out);
        let metadata = fixture_metadata();

        let writer = CohortVcfWriter::new(metadata, cfg_for(out.clone())).unwrap();
        // Constructor wrote the header into `tmp` already.
        assert!(tmp.exists(), "tmp file present after constructor");

        writer.abort().expect("abort succeeded");
        assert!(!tmp.exists(), "abort removed the tmp file");
        assert!(!out.exists(), "abort did not rename to final path");
    }

    /// B1: a second `abort()` call (or, equivalently, calling `abort`
    /// after the tmp file was already removed elsewhere) surfaces the
    /// underlying `io::Error` so the driver can log it rather than
    /// silently swallow it.
    #[test]
    fn abort_surfaces_remove_file_error_when_tmp_already_gone() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let tmp = tmp_path_for(&out);
        let metadata = fixture_metadata();

        let writer = CohortVcfWriter::new(metadata, cfg_for(out)).unwrap();
        std::fs::remove_file(&tmp).expect("manual remove succeeded");
        let err = writer.abort().unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::NotFound);
    }

    #[test]
    fn out_of_order_record_errors() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let metadata = fixture_metadata();
        let mut writer = CohortVcfWriter::new(metadata, cfg_for(out)).unwrap();
        writer.write_record(&mk_record(200)).unwrap();
        let err = writer.write_record(&mk_record(100)).unwrap_err();
        assert!(matches!(
            err,
            VcfWriteError::RecordOutOfOrder {
                chrom_id: 0,
                pos: 100,
                prev_chrom_id: 0,
                prev_pos: 200,
            }
        ));
    }

    #[test]
    fn equal_locus_records_also_error() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let metadata = fixture_metadata();
        let mut writer = CohortVcfWriter::new(metadata, cfg_for(out)).unwrap();
        writer.write_record(&mk_record(200)).unwrap();
        let err = writer.write_record(&mk_record(200)).unwrap_err();
        assert!(matches!(err, VcfWriteError::RecordOutOfOrder { .. }));
    }

    /// Mi1: an `RecordOutOfOrder` error does not advance
    /// `last_locus`. A subsequent record at a position past the
    /// *previous accepted* locus must still error against that
    /// locus, not the rejected one.
    #[test]
    fn out_of_order_does_not_advance_last_locus() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let metadata = fixture_metadata();
        let mut writer = CohortVcfWriter::new(metadata, cfg_for(out)).unwrap();
        writer.write_record(&mk_record(200)).unwrap();
        let _ = writer.write_record(&mk_record(100)).unwrap_err();
        // Subsequent record at 150 should still error against last
        // accepted (200), not against the rejected 100.
        let err = writer.write_record(&mk_record(150)).unwrap_err();
        assert!(matches!(
            err,
            VcfWriteError::RecordOutOfOrder {
                prev_pos: 200,
                pos: 150,
                ..
            }
        ));
        // And a record at 201 should succeed.
        writer.write_record(&mk_record(201)).unwrap();
        writer.finish().unwrap();
    }

    /// M12: `config()` accessor returns the frozen config; useful for
    /// runtime inspection without having to keep the original
    /// `WriterConfig` around.
    #[test]
    fn config_accessor_exposes_frozen_settings() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let metadata = fixture_metadata();
        let writer = CohortVcfWriter::new(metadata, cfg_for(out.clone())).unwrap();
        let cfg = writer.config();
        assert_eq!(cfg.output, out);
        assert!(!cfg.emit_gp);
        writer.finish().unwrap();
    }

    /// Mi17: the second record at the same `(ploidy, n_alleles)` is
    /// served from the cache. Asserted indirectly by observing that
    /// the cache grows from 0 → 1 across the first record and then
    /// stays at 1 across the second (same shape).
    #[test]
    fn genotype_table_cache_is_keyed_by_ploidy_and_n_alleles() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let metadata = fixture_metadata();
        let mut writer = CohortVcfWriter::new(metadata, cfg_for(out)).unwrap();
        assert_eq!(writer.genotype_tables.len(), 0);
        writer.write_record(&mk_record(100)).unwrap();
        assert_eq!(writer.genotype_tables.len(), 1, "biallelic table cached");
        writer.write_record(&mk_record(200)).unwrap();
        assert_eq!(
            writer.genotype_tables.len(),
            1,
            "second biallelic record reuses cache"
        );
        writer.finish().unwrap();
    }
}
