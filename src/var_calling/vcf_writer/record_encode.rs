//! Map a [`PosteriorRecord`] onto a `noodles_vcf::variant::RecordBuf`
//! ready for the writer to serialise.
//!
//! This module owns the field-by-field translation; everything
//! upstream is project types and everything downstream is noodles
//! types. Both directions are pure functions — no I/O, no state.
//!
//! The mapping is documented in
//! `doc/devel/implementation_plans/cohort_vcf_writer.md` §"Record
//! encoding (`record_encode.rs`)". Highlights:
//!
//! * `QUAL` is clamped to `[0, 9999]` (the spec-pinned cap; freebayes
//!   and GATK use the same).
//! * `GT` is unphased — Stage 5 does not currently emit phase — and
//!   the multiset of allele indices is sorted ascending (`0/1` not
//!   `1/0`) per the VCF spec.
//! * `GP` is only assembled when `config.emit_gp` is true.

use std::fmt::Write as _;

use noodles_core::Position;
use noodles_vcf::variant::RecordBuf;
use noodles_vcf::variant::record_buf::{
    AlternateBases, Filters, Ids, Info, Samples, info::field::Value as InfoValue,
    info::field::value::Array as InfoArray, samples::Keys, samples::sample::Value as SampleValue,
};

use super::WriterConfig;
use super::errors::VcfWriteError;
use crate::per_sample_pileup::psp::header::ParsedChromosome;
use crate::var_calling::per_group_merger::genotype_order;
use crate::var_calling::posterior_engine::PosteriorRecord;

/// Hard cap on the QUAL column.
///
/// `PosteriorRecord.qual_phred` is `f64::INFINITY` when every sample
/// is certainly variant. VCF has no sentinel for "infinite quality";
/// we cap at the same value GATK and freebayes use so downstream
/// tooling parses it as a number rather than choking. Search the
/// codebase for `9999` to find this comment.
pub(super) const QUAL_MAX: f32 = 9999.0;

/// Build the four-or-five-key FORMAT-keys list, once per writer
/// instance. Stored on `CohortVcfWriter` and cloned into each
/// `Samples` so we don't rebuild it per record.
pub(super) fn build_format_keys(config: &WriterConfig) -> Keys {
    let mut keys: Vec<String> = vec!["GT".into(), "GQ".into(), "DP".into(), "AD".into()];
    if config.emit_gp {
        keys.push("GP".into());
    }
    keys.into_iter().collect()
}

/// Encode a single posterior record into a noodles `RecordBuf`.
///
/// `contigs` is the same slice the writer holds — used to map
/// `record.locus.chrom_id` back to the contig name. `format_keys` is
/// pre-built via [`build_format_keys`] and shared across records.
pub(super) fn encode(
    record: &PosteriorRecord,
    contigs: &[ParsedChromosome],
    config: &WriterConfig,
    format_keys: &Keys,
    expected_samples: usize,
) -> Result<RecordBuf, VcfWriteError> {
    if record.n_samples != expected_samples {
        return Err(VcfWriteError::SampleCountMismatch {
            chrom_id: record.locus.chrom_id,
            pos: record.locus.start,
            expected_samples,
            got_samples: record.n_samples,
        });
    }

    let chrom_name = contigs
        .get(record.locus.chrom_id as usize)
        .map(|c| c.name.clone())
        .ok_or(VcfWriteError::UnknownChromId {
            chrom_id: record.locus.chrom_id,
            pos: record.locus.start,
            n_contigs: contigs.len(),
        })?;

    let pos_usize = usize::try_from(record.locus.start).map_err(|_| {
        VcfWriteError::Encode(format!(
            "record at {}:{}: 1-based position is not representable as usize",
            record.locus.chrom_id, record.locus.start
        ))
    })?;
    let position = Position::try_from(pos_usize).map_err(|e| {
        VcfWriteError::Encode(format!(
            "record at {}:{}: invalid VCF position: {e}",
            record.locus.chrom_id, record.locus.start
        ))
    })?;

    let ref_bases = allele_to_string(
        &record.alleles[0].seq,
        record.locus.chrom_id,
        record.locus.start,
    )?;
    let alt_strings: Vec<String> = record.alleles[1..]
        .iter()
        .map(|a| allele_to_string(&a.seq, record.locus.chrom_id, record.locus.start))
        .collect::<Result<_, _>>()?;
    let alternate_bases = AlternateBases::from(alt_strings);

    let qual = clamp_qual(record.qual_phred);

    let filters = if config.default_filter_pass {
        Filters::pass()
    } else {
        Filters::default()
    };

    let n_alts = record.alleles.len().saturating_sub(1);
    let table = genotype_order(record.ploidy, record.alleles.len());

    let info = build_info(record, n_alts, &table)?;
    let samples = build_samples(record, format_keys, config, &table)?;

    Ok(RecordBuf::builder()
        .set_reference_sequence_name(chrom_name)
        .set_variant_start(position)
        .set_ids(Ids::default())
        .set_reference_bases(ref_bases)
        .set_alternate_bases(alternate_bases)
        .set_quality_score(qual)
        .set_filters(filters)
        .set_info(info)
        .set_samples(samples)
        .build())
}

/// Convert a `Vec<u8>` allele sequence (guaranteed by Stage 1 to be
/// uppercase `{A,C,G,T,N}`) into an owned `String` for the VCF
/// encoder. An unexpected byte surfaces as an `Encode` error rather
/// than a panic — same defensive shape used at other library /
/// noodles boundaries.
fn allele_to_string(seq: &[u8], chrom_id: u32, pos: u32) -> Result<String, VcfWriteError> {
    std::str::from_utf8(seq)
        .map(|s| s.to_string())
        .map_err(|e| {
            VcfWriteError::Encode(format!(
                "record at {chrom_id}:{pos}: allele bytes are not valid UTF-8: {e}"
            ))
        })
}

/// Clamp QUAL to `[0, 9999]`. Inf maps to 9999 (every sample
/// certainly variant); NaN maps to 0 defensively — the EM math
/// should never produce NaN for an emitted record, so this is a
/// guard not an expected path.
fn clamp_qual(qual_phred: f64) -> f32 {
    if qual_phred.is_nan() {
        return 0.0;
    }
    let q = qual_phred as f32;
    q.clamp(0.0, QUAL_MAX)
}

fn build_info(
    record: &PosteriorRecord,
    n_alts: usize,
    table: &[Vec<u8>],
) -> Result<Info, VcfWriteError> {
    let mut info = Info::default();

    // AF — allele frequencies for ALT alleles only (REF is implied).
    // `Number=A`: empty list when no ALTs.
    let af_values: Vec<Option<f32>> = if n_alts == 0 {
        Vec::new()
    } else {
        record.allele_frequencies[1..]
            .iter()
            .map(|f| Some(*f as f32))
            .collect()
    };
    if !af_values.is_empty() {
        info.as_mut().insert(
            "AF".into(),
            Some(InfoValue::Array(InfoArray::Float(af_values))),
        );
    }

    // AC — per-ALT count of allele copies in argmax genotypes summed
    // across samples; AN — total called allele number.
    let (ac_per_alt, an_total) = tally_called_alleles(record, n_alts, table)?;
    if !ac_per_alt.is_empty() {
        let ac_values: Vec<Option<i32>> = ac_per_alt.iter().map(|&c| Some(c as i32)).collect();
        info.as_mut().insert(
            "AC".into(),
            Some(InfoValue::Array(InfoArray::Integer(ac_values))),
        );
    }
    info.as_mut()
        .insert("AN".into(), Some(InfoValue::Integer(an_total as i32)));

    // DP — total depth across samples = sum of all per-allele
    // observation counts.
    let dp_total: u64 = (0..record.n_samples)
        .map(|s| {
            record
                .scalars_row(s)
                .iter()
                .map(|stats| u64::from(stats.num_obs))
                .sum::<u64>()
        })
        .sum();
    info.as_mut().insert(
        "DP".into(),
        Some(InfoValue::Integer(
            i32::try_from(dp_total).unwrap_or(i32::MAX),
        )),
    );

    // CA — flag set when any sample carries a chain-anchor-broken
    // call. `chain_anchor_flags` is a flat `Vec<bool>` of length
    // `n_samples * n_alleles`.
    if record.chain_anchor_flags.iter().any(|&b| b) {
        info.as_mut().insert("CA".into(), Some(InfoValue::Flag));
    }

    Ok(info)
}

/// Walk each sample's argmax genotype, decode it through the
/// `genotype_order` table, and tally per-ALT counts plus the total
/// called allele number.
fn tally_called_alleles(
    record: &PosteriorRecord,
    n_alts: usize,
    table: &[Vec<u8>],
) -> Result<(Vec<u32>, u32), VcfWriteError> {
    let mut ac_per_alt = vec![0u32; n_alts];
    let mut an_total: u32 = 0;
    for (sample_idx, &gt_idx) in record.best_genotype.iter().enumerate() {
        let gt = table
            .get(gt_idx)
            .ok_or(VcfWriteError::GenotypeIndexOutOfBounds {
                chrom_id: record.locus.chrom_id,
                pos: record.locus.start,
                sample_idx,
                got: gt_idx,
                n_genotypes: table.len(),
            })?;
        for &allele_idx in gt {
            an_total += 1;
            if allele_idx == 0 {
                continue; // REF
            }
            let alt_pos = (allele_idx as usize) - 1;
            // Defensive bounds check; allele_idx > n_alts would mean
            // the upstream produced a genotype indexing an allele
            // beyond the record's allele set.
            if alt_pos < ac_per_alt.len() {
                ac_per_alt[alt_pos] += 1;
            }
        }
    }
    Ok((ac_per_alt, an_total))
}

fn build_samples(
    record: &PosteriorRecord,
    format_keys: &Keys,
    config: &WriterConfig,
    table: &[Vec<u8>],
) -> Result<Samples, VcfWriteError> {
    let n_alleles = record.alleles.len();
    let mut rows: Vec<Vec<Option<SampleValue>>> = Vec::with_capacity(record.n_samples);
    let mut gt_buf = String::new();
    for sample_idx in 0..record.n_samples {
        let gt_idx = record.best_genotype[sample_idx];
        let gt = table
            .get(gt_idx)
            .ok_or(VcfWriteError::GenotypeIndexOutOfBounds {
                chrom_id: record.locus.chrom_id,
                pos: record.locus.start,
                sample_idx,
                got: gt_idx,
                n_genotypes: table.len(),
            })?;
        gt_buf.clear();
        format_gt_unphased(&mut gt_buf, gt);

        let gq = record.gq_phred[sample_idx].round().clamp(0.0, 99.0) as i32;

        let scalars = record.scalars_row(sample_idx);
        let dp_sample: u64 = scalars.iter().map(|stats| u64::from(stats.num_obs)).sum();
        let dp_value = i32::try_from(dp_sample).unwrap_or(i32::MAX);

        let ad_values: Vec<Option<i32>> = (0..n_alleles)
            .map(|a| Some(scalars[a].num_obs as i32))
            .collect();

        let mut row: Vec<Option<SampleValue>> = Vec::with_capacity(format_keys.as_ref().len());
        row.push(Some(SampleValue::String(gt_buf.clone())));
        row.push(Some(SampleValue::Integer(gq)));
        row.push(Some(SampleValue::Integer(dp_value)));
        row.push(Some(SampleValue::from(ad_values)));

        if config.emit_gp {
            let gp_values: Vec<Option<f32>> = record
                .posteriors_row(sample_idx)
                .iter()
                .map(|p| Some(*p as f32))
                .collect();
            row.push(Some(SampleValue::from(gp_values)));
        }

        rows.push(row);
    }

    Ok(Samples::new(format_keys.clone(), rows))
}

/// Format a multiset of allele indices as an unphased `a/b[/c…]`
/// genotype string. The input is assumed already sorted ascending
/// (the `genotype_order` table builds them in canonical order, so a
/// re-sort would be redundant for now — kept as a function so a
/// future phased-emission flag is a localised change).
fn format_gt_unphased(out: &mut String, alleles: &[u8]) {
    for (i, &a) in alleles.iter().enumerate() {
        if i > 0 {
            out.push('/');
        }
        // `write!` into a String never fails; ignore the result.
        let _ = write!(out, "{a}");
    }
}

#[cfg(test)]
mod tests {
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
        MergedAllele {
            seq: seq.to_vec(),
            is_compound: false,
            constituents: Vec::new(),
        }
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

    /// 2-sample biallelic SNP A→T; best_genotype = [0/0, 0/1].
    /// Genotype enumeration at ploidy=2, n_alleles=2 is
    /// [[0,0], [0,1], [1,1]], so indices 0 → 0/0 and 1 → 0/1.
    fn biallelic_two_samples() -> PosteriorRecord {
        let alleles = vec![ref_allele(b"A"), alt_allele(b"T")];
        // Row-major: sample × allele.
        let scalars = vec![
            support(20), // S0, REF
            support(0),  // S0, ALT
            support(10), // S1, REF
            support(10), // S1, ALT
        ];
        let posteriors = vec![
            0.98, 0.01, 0.01, // S0 over 3 genotypes
            0.05, 0.90, 0.05, // S1
        ];
        let best_genotype = vec![0usize, 1usize];
        PosteriorRecord {
            locus: RecordLocus {
                chrom_id: 0,
                start: 1234,
                end: 1234,
            },
            alleles,
            ploidy: 2,
            n_samples: 2,
            n_genotypes: 3,
            allele_frequencies: vec![0.75, 0.25],
            compound_frequencies: vec![None, None],
            posteriors,
            best_genotype,
            gq_phred: vec![60.0, 40.0],
            qual_phred: 200.0,
            scalars,
            other_scalars: vec![],
            chain_anchor_flags: vec![false; 4],
            diagnostics: EmDiagnostics {
                iterations: 5,
                final_max_delta_p: 1e-6,
            },
        }
    }

    fn fixture_contigs() -> Vec<ParsedChromosome> {
        vec![ParsedChromosome {
            name: "chr1".into(),
            length: 1_000_000,
            md5: "00000000000000000000000000000001".into(),
        }]
    }

    /// Drive a `RecordBuf` through the writer to a String; lets tests
    /// inspect the actual rendered bytes.
    fn record_to_line(record: &RecordBuf, header: &noodles_vcf::Header) -> String {
        use noodles_vcf::variant::io::Write as _;
        let mut writer = noodles_vcf::io::Writer::new(Vec::new());
        writer.write_variant_record(header, record).unwrap();
        String::from_utf8(writer.into_inner()).unwrap()
    }

    fn fixture_header(config: &WriterConfig) -> noodles_vcf::Header {
        use super::super::header::{CohortMetadata, build_vcf_header};
        let metadata = CohortMetadata {
            sample_names: vec!["S0".into(), "S1".into()],
            contigs: fixture_contigs(),
            tool_string: "test 0".into(),
            command_line: String::new(),
        };
        build_vcf_header(&metadata, config).unwrap()
    }

    #[test]
    fn biallelic_snp_default_config_renders_expected_line() {
        let record = biallelic_two_samples();
        let contigs = fixture_contigs();
        let config = WriterConfig::default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);

        let buf = encode(&record, &contigs, &config, &keys, 2).unwrap();
        let line = record_to_line(&buf, &header);

        // The single line ends with a newline; one record only.
        let line = line.trim_end();
        let fields: Vec<&str> = line.split('\t').collect();
        // CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S0 S1
        assert_eq!(fields.len(), 11, "unexpected field count: {line}");
        assert_eq!(fields[0], "chr1");
        assert_eq!(fields[1], "1234");
        assert_eq!(fields[2], ".");
        assert_eq!(fields[3], "A");
        assert_eq!(fields[4], "T");
        // QUAL is f32-rounded; just check the integer prefix.
        assert!(fields[5].starts_with("200"), "QUAL = {}", fields[5]);
        assert_eq!(fields[6], "PASS");
        // INFO contains all of AF/AC/AN/DP and no CA flag.
        let info = fields[7];
        assert!(info.contains("AF=0.25"), "INFO missing AF: {info}");
        assert!(info.contains("AC=1"), "INFO missing AC: {info}");
        // AN: S0=0/0 (2 called) + S1=0/1 (2 called) → 4
        assert!(info.contains("AN=4"), "INFO missing AN: {info}");
        // DP: total = 20+0+10+10 = 40
        assert!(info.contains("DP=40"), "INFO missing DP: {info}");
        assert!(!info.contains("CA"), "INFO carried unexpected CA: {info}");
        // FORMAT key string — default config: no GP.
        assert_eq!(fields[8], "GT:GQ:DP:AD");
        // Per-sample cells: 4 fields each.
        assert_eq!(fields[9], "0/0:60:20:20,0", "S0 cell: {}", fields[9]);
        assert_eq!(fields[10], "0/1:40:20:10,10", "S1 cell: {}", fields[10]);
    }

    #[test]
    fn biallelic_snp_with_emit_gp_adds_gp_column() {
        let record = biallelic_two_samples();
        let contigs = fixture_contigs();
        let config = WriterConfig {
            output: "/dev/null".into(),
            default_filter_pass: true,
            emit_gp: true,
        };
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);

        let buf = encode(&record, &contigs, &config, &keys, 2).unwrap();
        let line = record_to_line(&buf, &header).trim_end().to_string();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[8], "GT:GQ:DP:AD:GP");
        // Per-sample cells now carry 5 colon-separated subfields.
        let s0 = fields[9];
        let s0_parts: Vec<&str> = s0.split(':').collect();
        assert_eq!(s0_parts.len(), 5, "S0 with GP: {s0}");
        // GP for S0: 0.98, 0.01, 0.01 over 3 genotypes
        let gp_floats: Vec<f32> = s0_parts[4].split(',').map(|s| s.parse().unwrap()).collect();
        assert_eq!(gp_floats.len(), 3);
        assert!((gp_floats.iter().sum::<f32>() - 1.0).abs() < 1e-3);
    }

    #[test]
    fn triallelic_snp_emits_two_alts_with_number_a_lists() {
        // Same template but with three alleles: A→T, A→C.
        let mut record = biallelic_two_samples();
        record.alleles.push(alt_allele(b"C"));
        // n_genotypes for ploidy=2, n_alleles=3 is 6:
        // [0,0],[0,1],[1,1],[0,2],[1,2],[2,2]
        record.n_genotypes = 6;
        record.allele_frequencies = vec![0.6, 0.25, 0.15];
        record.compound_frequencies = vec![None, None, None];
        record.posteriors = vec![
            0.7, 0.1, 0.05, 0.1, 0.025, 0.025, // S0
            0.05, 0.05, 0.05, 0.7, 0.1, 0.05, // S1 → best = 3 (0/2)
        ];
        record.best_genotype = vec![0, 3]; // S0=0/0, S1=0/2
        record.scalars = vec![
            support(20),
            support(0),
            support(0), // S0 row
            support(10),
            support(0),
            support(10), // S1 row
        ];
        record.chain_anchor_flags = vec![false; 6];

        let contigs = fixture_contigs();
        let config = WriterConfig::default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);

        let buf = encode(&record, &contigs, &config, &keys, 2).unwrap();
        let line = record_to_line(&buf, &header).trim_end().to_string();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[3], "A");
        assert_eq!(fields[4], "T,C");
        let info = fields[7];
        // AF list has 2 entries (Number=A).
        assert!(info.contains("AF=0.25,0.15"), "INFO: {info}");
        // AC list: ALT1 (T) = 0 copies, ALT2 (C) = 1 copy (S1's 0/2).
        assert!(info.contains("AC=0,1"), "INFO: {info}");
    }

    #[test]
    fn chain_anchor_flag_surfaces_as_ca_info_flag() {
        let mut record = biallelic_two_samples();
        // Flip the (S1, ALT) anchor-broken flag.
        record.chain_anchor_flags[3] = true;
        let contigs = fixture_contigs();
        let config = WriterConfig::default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);

        let buf = encode(&record, &contigs, &config, &keys, 2).unwrap();
        let line = record_to_line(&buf, &header).trim_end().to_string();
        let fields: Vec<&str> = line.split('\t').collect();
        let info = fields[7];
        assert!(info.contains("CA"), "INFO should carry CA: {info}");
    }

    #[test]
    fn infinite_qual_renders_as_9999_cap() {
        let mut record = biallelic_two_samples();
        record.qual_phred = f64::INFINITY;
        let contigs = fixture_contigs();
        let config = WriterConfig::default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);

        let buf = encode(&record, &contigs, &config, &keys, 2).unwrap();
        let line = record_to_line(&buf, &header).trim_end().to_string();
        let fields: Vec<&str> = line.split('\t').collect();
        // QUAL field is the 6th column (index 5).
        let qual: f32 = fields[5].parse().unwrap();
        assert_eq!(qual, QUAL_MAX, "QUAL should be capped at 9999, got {qual}");
    }

    #[test]
    fn nan_qual_renders_as_zero() {
        let mut record = biallelic_two_samples();
        record.qual_phred = f64::NAN;
        let contigs = fixture_contigs();
        let config = WriterConfig::default();
        let keys = build_format_keys(&config);
        let header = fixture_header(&config);

        let buf = encode(&record, &contigs, &config, &keys, 2).unwrap();
        let line = record_to_line(&buf, &header).trim_end().to_string();
        let fields: Vec<&str> = line.split('\t').collect();
        let qual: f32 = fields[5].parse().unwrap();
        assert_eq!(qual, 0.0);
    }

    #[test]
    fn unknown_chrom_id_is_an_error() {
        let mut record = biallelic_two_samples();
        record.locus.chrom_id = 99;
        let contigs = fixture_contigs();
        let config = WriterConfig::default();
        let keys = build_format_keys(&config);
        let err = encode(&record, &contigs, &config, &keys, 2).unwrap_err();
        assert!(matches!(
            err,
            VcfWriteError::UnknownChromId { chrom_id: 99, .. }
        ));
    }

    #[test]
    fn sample_count_mismatch_is_an_error() {
        let record = biallelic_two_samples();
        let contigs = fixture_contigs();
        let config = WriterConfig::default();
        let keys = build_format_keys(&config);
        // Cohort metadata names 3 samples but the record has 2.
        let err = encode(&record, &contigs, &config, &keys, 3).unwrap_err();
        assert!(matches!(
            err,
            VcfWriteError::SampleCountMismatch {
                expected_samples: 3,
                got_samples: 2,
                ..
            }
        ));
    }

    #[test]
    fn format_keys_track_emit_gp() {
        let off = build_format_keys(&WriterConfig::default());
        let on_cfg = WriterConfig {
            emit_gp: true,
            ..Default::default()
        };
        let on = build_format_keys(&on_cfg);

        let off_vec: Vec<&str> = off.as_ref().iter().map(String::as_str).collect();
        let on_vec: Vec<&str> = on.as_ref().iter().map(String::as_str).collect();
        assert_eq!(off_vec, vec!["GT", "GQ", "DP", "AD"]);
        assert_eq!(on_vec, vec!["GT", "GQ", "DP", "AD", "GP"]);
    }

    #[test]
    fn format_gt_unphased_joins_with_slashes() {
        let mut s = String::new();
        format_gt_unphased(&mut s, &[0, 1]);
        assert_eq!(s, "0/1");
        s.clear();
        format_gt_unphased(&mut s, &[1, 1, 2]);
        assert_eq!(s, "1/1/2");
    }
}
