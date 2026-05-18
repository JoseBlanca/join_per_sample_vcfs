//! [`CohortVcfWriter`] — the runtime that ties sink, header, and
//! record-encoder together.
//!
//! `new` opens `<output>.tmp`, builds the header, writes it, and
//! parks the encoder; `write_record` encodes one `PosteriorRecord`
//! into a noodles `RecordBuf` and feeds it to the underlying VCF
//! writer; `finish` flushes the sink (including the bgzf EOF block
//! for gzipped paths), `fsync`s, and atomic-renames into place.
//!
//! Forgetting `finish` leaves `<output>.tmp` on disk and no
//! `<output>` — intended loud failure.

use noodles_vcf::Header as VcfHeader;
use noodles_vcf::variant::io::Write as _;
use noodles_vcf::variant::record_buf::samples::Keys;

use super::WriterConfig;
use super::errors::VcfWriteError;
use super::header::{CohortMetadata, build_vcf_header};
use super::record_encode::{build_format_keys, encode};
use super::sink::SinkKind;
use crate::per_sample_pileup::psp::header::ParsedChromosome;
use crate::var_calling::posterior_engine::PosteriorRecord;

/// Streams `PosteriorRecord` items to a VCF file.
///
/// Construction opens `<output>.tmp`, builds the noodles header, and
/// writes it immediately so a subsequent `write_record` does not pay
/// header-build cost. Per-record encoding goes through
/// [`super::record_encode::encode`]. The writer's
/// non-decreasing-locus order check is the last line of defence
/// against an upstream that surfaces records out of order.
pub struct CohortVcfWriter {
    /// noodles' line-oriented VCF writer wraps our sink directly. The
    /// sink (plain or bgzf) is owned by the writer and recovered on
    /// `finish`.
    inner: noodles_vcf::io::Writer<SinkKind>,
    /// Re-encoded by the writer on every record write, but kept here
    /// once because noodles takes `&Header` per call.
    header: VcfHeader,
    /// Cohort metadata frozen at construction. Only the `contigs`
    /// slice is read on the per-record path (to map `chrom_id` →
    /// CHROM string); the sample-name vector is consumed by the
    /// header build at construction and not re-read.
    contigs: Vec<ParsedChromosome>,
    /// Cohort size, re-checked against every record's `n_samples`.
    expected_samples: usize,
    /// Final path; tmp is `<final_path>.tmp`. Held until `finish`.
    final_path: std::path::PathBuf,
    /// FORMAT keys list — built once from `config.emit_gp` and reused
    /// per record.
    format_keys: Keys,
    /// Whole-run config — held by value for the per-record encoder.
    config: WriterConfig,
    /// Latches the most recent `(chrom_id, start)` so the next
    /// record can be checked for non-decreasing order.
    last_locus: Option<(u32, u32)>,
}

impl CohortVcfWriter {
    /// Open `<config.output>.tmp`, build the VCF header, and write it.
    pub fn new(metadata: CohortMetadata, config: WriterConfig) -> Result<Self, VcfWriteError> {
        let header = build_vcf_header(&metadata, &config)?;
        let format_keys = build_format_keys(&config);
        let contigs = metadata.contigs.clone();
        let expected_samples = metadata.sample_names.len();
        let final_path = config.output.clone();

        let sink = SinkKind::open_tmp(&final_path)?;
        let mut inner = noodles_vcf::io::Writer::new(sink);
        inner.write_header(&header)?;
        Ok(Self {
            inner,
            header,
            contigs,
            expected_samples,
            final_path,
            format_keys,
            config,
            last_locus: None,
        })
    }

    /// Encode and write one record. Latches the locus order; an
    /// upstream regression surfaces as
    /// [`VcfWriteError::RecordOutOfOrder`] and the writer keeps
    /// running (the next call gets a fresh order check against the
    /// previous accepted record).
    pub fn write_record(&mut self, record: &PosteriorRecord) -> Result<(), VcfWriteError> {
        let locus = (record.locus.chrom_id, record.locus.start);
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

        let buf = encode(
            record,
            &self.contigs,
            &self.config,
            &self.format_keys,
            self.expected_samples,
        )?;
        self.inner.write_variant_record(&self.header, &buf)?;
        self.last_locus = Some(locus);
        Ok(())
    }

    /// Flush, emit the bgzf EOF block when applicable, sync, and
    /// atomic-rename `<output>.tmp` → `<output>`.
    pub fn finish(self) -> Result<(), VcfWriteError> {
        let sink = self.inner.into_inner();
        sink.finish(&self.final_path)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use tempfile::tempdir;

    use super::*;
    use crate::per_sample_pileup::pileup::AlleleSupportStats;
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
            },
        }
    }

    #[test]
    fn full_round_trip_plain_text() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let metadata = fixture_metadata();
        let config = WriterConfig {
            output: out.clone(),
            default_filter_pass: true,
            emit_gp: false,
        };

        let mut writer = CohortVcfWriter::new(metadata, config).unwrap();
        writer.write_record(&mk_record(100)).unwrap();
        writer.write_record(&mk_record(200)).unwrap();
        writer.finish().unwrap();

        let text = std::fs::read_to_string(&out).unwrap();
        assert!(text.contains("##fileformat=VCFv4.4"));
        let data_lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 2);
        assert!(data_lines[0].starts_with("chr1\t100\t"));
        assert!(data_lines[1].starts_with("chr1\t200\t"));
        // Tmp file is gone.
        assert!(!out.with_extension("vcf.tmp").exists());
    }

    #[test]
    fn out_of_order_record_errors() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("out.vcf");
        let metadata = fixture_metadata();
        let config = WriterConfig {
            output: out,
            default_filter_pass: true,
            emit_gp: false,
        };
        let mut writer = CohortVcfWriter::new(metadata, config).unwrap();
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
        let config = WriterConfig {
            output: out,
            default_filter_pass: true,
            emit_gp: false,
        };
        let mut writer = CohortVcfWriter::new(metadata, config).unwrap();
        writer.write_record(&mk_record(200)).unwrap();
        let err = writer.write_record(&mk_record(200)).unwrap_err();
        assert!(matches!(err, VcfWriteError::RecordOutOfOrder { .. }));
    }
}
