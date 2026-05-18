//! Build a `noodles_vcf::Header` from the project's
//! [`CohortMetadata`].
//!
//! The header carries `##fileformat=VCFv4.4`, optional `##source` and
//! `##commandline` unstructured records, one `##contig=` per entry in
//! `metadata.contigs`, the standard INFO / FORMAT definitions for the
//! fields the writer emits, and the cohort sample-name vector. The
//! `GP` FORMAT declaration is only present when `config.emit_gp` is
//! true so a strict parser doesn't flag a declared-but-unused field.

use std::collections::HashSet;

use noodles_vcf::Header;
use noodles_vcf::header::Builder as HeaderBuilder;
use noodles_vcf::header::FileFormat;
use noodles_vcf::header::record::Value as HeaderValue;
use noodles_vcf::header::record::key::other::ParseError as OtherKeyParseError;
use noodles_vcf::header::record::value::Map;
use noodles_vcf::header::record::value::map::{
    Contig, Filter, Format, Info,
    format::{Number as FormatNumber, Type as FormatType},
    info::{Number as InfoNumber, Type as InfoType},
};

use super::WriterConfig;
use super::errors::VcfWriteError;
use crate::per_sample_pileup::psp::header::ParsedChromosome;

/// Inputs the writer needs at construction time. Built by the CLI
/// from `merger.sample_names()`, `merger.chromosomes()`, and the
/// invocation context.
#[derive(Debug, Clone)]
pub struct CohortMetadata {
    /// Sample names in the cohort, ordered to match the layout of
    /// `PosteriorRecord.posteriors` rows. The writer emits the
    /// `#CHROM ... FORMAT SAMPLE_0 SAMPLE_1 ...` data-header line
    /// in this order, so it must come from the upstream merger's
    /// `sample_names()` (the single source of truth) — not from the
    /// CLI argument order, even when the two look identical.
    pub sample_names: Vec<String>,
    /// Contig table — name, length, md5. Sourced from the per-position
    /// merger's chromosome slice, which in turn was sourced from the
    /// `.psp` headers.
    ///
    /// An empty vector is accepted (the header will carry no
    /// `##contig=` lines); operationally meaningless but structurally
    /// valid.
    pub contigs: Vec<ParsedChromosome>,
    /// Goes into `##source=...`. Conventionally `"pop_var_caller
    /// <version>"`. An empty string skips the `##source` line.
    pub tool_string: String,
    /// Goes into `##commandline=...`. CLIs build this from
    /// `std::env::args()`; library users may pass an empty string to
    /// skip the field.
    pub command_line: String,
}

pub(super) fn build_vcf_header(
    metadata: &CohortMetadata,
    config: &WriterConfig,
) -> Result<Header, VcfWriteError> {
    validate_metadata(metadata)?;

    let mut builder = Header::builder().set_file_format(FileFormat::new(4, 4));

    builder = insert_unstructured(builder, "source", &metadata.tool_string)?;
    builder = insert_unstructured(builder, "commandline", &metadata.command_line)?;

    for entry in &metadata.contigs {
        // Mi4: `usize::try_from(u32)` is infallible on 32-bit and
        // 64-bit targets (the project deploys to x86_64 and aarch64).
        // The only real cap is the VCF spec's `##contig=length=` field
        // which htslib treats as signed 32-bit. One check suffices.
        if entry.length > i32::MAX as u32 {
            return Err(VcfWriteError::ContigLengthOverflow {
                name: entry.name.clone(),
                length: entry.length,
            });
        }
        let length_usize = entry.length as usize;
        let mut contig = Map::<Contig>::new();
        *contig.length_mut() = Some(length_usize);
        *contig.md5_mut() = Some(entry.md5.clone());
        builder = builder.add_contig(entry.name.clone(), contig);
    }

    // PASS filter — every record carries PASS by default; declaring
    // it in the header makes bcftools / strict consumers happy.
    builder = builder.add_filter("PASS", Map::<Filter>::pass());

    // INFO entries.
    builder = builder
        .add_info(
            "AF",
            Map::<Info>::new(
                InfoNumber::AlternateBases,
                InfoType::Float,
                "Allele frequency from EM",
            ),
        )
        .add_info(
            "AC",
            Map::<Info>::new(
                InfoNumber::AlternateBases,
                InfoType::Integer,
                "Allele count in called genotypes",
            ),
        )
        .add_info(
            "AN",
            Map::<Info>::new(
                InfoNumber::Count(1),
                InfoType::Integer,
                "Total number of called alleles",
            ),
        )
        .add_info(
            "DP",
            Map::<Info>::new(
                InfoNumber::Count(1),
                InfoType::Integer,
                "Total depth across samples",
            ),
        )
        .add_info(
            "CA",
            Map::<Info>::new(
                InfoNumber::Count(0),
                InfoType::Flag,
                "At least one sample is a chain-anchor-broken compound call",
            ),
        );

    // FORMAT entries — always-on quartet, plus GP only when enabled.
    builder = builder
        .add_format(
            "GT",
            Map::<Format>::new(FormatNumber::Count(1), FormatType::String, "Genotype"),
        )
        .add_format(
            "GQ",
            Map::<Format>::new(
                FormatNumber::Count(1),
                FormatType::Integer,
                "Genotype quality (Phred)",
            ),
        )
        .add_format(
            "DP",
            Map::<Format>::new(FormatNumber::Count(1), FormatType::Integer, "Read depth"),
        )
        .add_format(
            "AD",
            Map::<Format>::new(
                FormatNumber::ReferenceAlternateBases,
                FormatType::Integer,
                "Allelic depths for REF + ALT",
            ),
        );
    if config.emit_gp {
        // noodles' name for VCF `Number=G` is `Number::Samples` (the
        // crate enum is keyed by the spec letter — see
        // header/record/value/map/format/number.rs).
        builder = builder.add_format(
            "GP",
            Map::<Format>::new(
                FormatNumber::Samples,
                FormatType::Float,
                "Genotype posterior probabilities",
            ),
        );
    }

    for name in &metadata.sample_names {
        builder = builder.add_sample_name(name.clone());
    }

    Ok(builder.build())
}

/// Mi6: insert one unstructured `##key=value` header record. Empty
/// `value` is a no-op (callers contract on it to skip the line). Both
/// the key-parse and the insert failure surface as
/// `VcfWriteError::Encode` with the project-side operation tag —
/// the noodles cause is boxed per the project's "no third-party types
/// in our public error API" decision.
fn insert_unstructured(
    builder: HeaderBuilder,
    key: &'static str,
    value: &str,
) -> Result<HeaderBuilder, VcfWriteError> {
    if value.is_empty() {
        return Ok(builder);
    }
    let parsed_key = key
        .parse::<noodles_vcf::header::record::key::Other>()
        .map_err(|e: OtherKeyParseError| VcfWriteError::Encode {
            operation: header_op_for_key(key, KeyOp::Parse),
            source: Box::new(e),
        })?;
    builder
        .insert(parsed_key, HeaderValue::from(value))
        .map_err(|e| VcfWriteError::Encode {
            operation: header_op_for_key(key, KeyOp::Insert),
            source: Box::new(e),
        })
}

enum KeyOp {
    Parse,
    Insert,
}

/// Produces a `&'static str` operation tag for the `Encode` variant.
/// The set of keys we ever insert is known at compile time, so a
/// small `match` keeps `Encode { operation: &'static str }` honest.
fn header_op_for_key(key: &'static str, op: KeyOp) -> &'static str {
    match (key, op) {
        ("source", KeyOp::Parse) => "##source key parse",
        ("source", KeyOp::Insert) => "##source insert",
        ("commandline", KeyOp::Parse) => "##commandline key parse",
        ("commandline", KeyOp::Insert) => "##commandline insert",
        (_, KeyOp::Parse) => "header key parse",
        (_, KeyOp::Insert) => "header insert",
    }
}

fn validate_metadata(metadata: &CohortMetadata) -> Result<(), VcfWriteError> {
    if metadata.sample_names.is_empty() {
        return Err(VcfWriteError::InvalidMetadata(
            "sample_names is empty; a cohort needs at least one sample".into(),
        ));
    }
    let mut seen_samples = HashSet::with_capacity(metadata.sample_names.len());
    for name in &metadata.sample_names {
        if !seen_samples.insert(name.as_str()) {
            return Err(VcfWriteError::InvalidMetadata(format!(
                "duplicate sample name `{name}`"
            )));
        }
    }

    let mut seen_contigs = HashSet::with_capacity(metadata.contigs.len());
    for c in &metadata.contigs {
        if !seen_contigs.insert(c.name.as_str()) {
            return Err(VcfWriteError::InvalidMetadata(format!(
                "duplicate contig name `{}`",
                c.name
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;

    /// Serialise a header through the actual VCF writer so the test
    /// matches what users will see on disk. `noodles_vcf::Header`
    /// does not implement `Display`; the writer is the rendering
    /// path.
    fn header_to_string(header: &Header) -> String {
        let mut writer = noodles_vcf::io::Writer::new(Vec::new());
        writer.write_header(header).unwrap();
        String::from_utf8(writer.into_inner()).unwrap()
    }

    fn fixture_metadata() -> CohortMetadata {
        CohortMetadata {
            sample_names: vec!["S0".into(), "S1".into(), "S2".into()],
            contigs: vec![
                ParsedChromosome {
                    name: "chr1".into(),
                    length: 1_000,
                    md5: "00000000000000000000000000000001".into(),
                },
                ParsedChromosome {
                    name: "chr2".into(),
                    length: 2_000,
                    md5: "00000000000000000000000000000002".into(),
                },
            ],
            tool_string: "pop_var_caller 0.1.0".into(),
            command_line: "pop_var_caller cohort --output out.vcf".into(),
        }
    }

    fn cfg_emit_gp_off() -> WriterConfig {
        WriterConfig::new(PathBuf::from("/dev/null"))
    }

    fn cfg_emit_gp_on() -> WriterConfig {
        let mut c = WriterConfig::new(PathBuf::from("/dev/null"));
        c.emit_gp = true;
        c
    }

    #[test]
    fn header_round_trips_without_gp_by_default() {
        let metadata = fixture_metadata();
        let header = build_vcf_header(&metadata, &cfg_emit_gp_off()).unwrap();
        let text = header_to_string(&header);
        assert!(text.contains("##fileformat=VCFv4.4"));
        assert!(text.contains("##source=pop_var_caller 0.1.0"));
        assert!(text.contains("##commandline="));
        assert!(text.contains("##contig=<ID=chr1"));
        assert!(text.contains("length=1000"));
        assert!(text.contains("md5=00000000000000000000000000000001"));
        assert!(text.contains("##INFO=<ID=AF"));
        assert!(text.contains("##INFO=<ID=CA,Number=0,Type=Flag"));
        assert!(text.contains("##FORMAT=<ID=GT"));
        assert!(text.contains("##FORMAT=<ID=AD,Number=R,Type=Integer"));
        assert!(
            !text.contains("##FORMAT=<ID=GP"),
            "GP should be absent when emit_gp = false"
        );
        let chrom_line = text.lines().find(|l| l.starts_with("#CHROM")).unwrap();
        assert!(chrom_line.ends_with("\tS0\tS1\tS2"));
    }

    #[test]
    fn header_includes_gp_when_enabled() {
        let metadata = fixture_metadata();
        let header = build_vcf_header(&metadata, &cfg_emit_gp_on()).unwrap();
        let text = header_to_string(&header);
        assert!(text.contains("##FORMAT=<ID=GP,Number=G,Type=Float"));
    }

    #[test]
    fn empty_sample_names_rejected() {
        let mut metadata = fixture_metadata();
        metadata.sample_names.clear();
        let err = build_vcf_header(&metadata, &cfg_emit_gp_off()).unwrap_err();
        assert!(matches!(err, VcfWriteError::InvalidMetadata(_)));
    }

    #[test]
    fn duplicate_sample_name_rejected() {
        let mut metadata = fixture_metadata();
        metadata.sample_names = vec!["S0".into(), "S0".into()];
        let err = build_vcf_header(&metadata, &cfg_emit_gp_off()).unwrap_err();
        match err {
            VcfWriteError::InvalidMetadata(msg) => {
                assert!(msg.contains("duplicate sample name"), "got: {msg}");
            }
            other => panic!("expected InvalidMetadata, got {other:?}"),
        }
    }

    #[test]
    fn duplicate_contig_name_rejected() {
        let mut metadata = fixture_metadata();
        metadata.contigs[1].name = "chr1".into();
        let err = build_vcf_header(&metadata, &cfg_emit_gp_off()).unwrap_err();
        match err {
            VcfWriteError::InvalidMetadata(msg) => {
                assert!(msg.contains("duplicate contig name"), "got: {msg}");
            }
            other => panic!("expected InvalidMetadata, got {other:?}"),
        }
    }

    /// Mi4: contig length over `i32::MAX` is rejected with the typed
    /// variant. Previously had zero coverage.
    #[test]
    fn contig_length_above_i32_max_rejected() {
        let mut metadata = fixture_metadata();
        let bad_length = (i32::MAX as u32) + 1;
        metadata.contigs[0].length = bad_length;
        let err = build_vcf_header(&metadata, &cfg_emit_gp_off()).unwrap_err();
        match err {
            VcfWriteError::ContigLengthOverflow { name, length } => {
                assert_eq!(name, "chr1");
                assert_eq!(length, bad_length);
            }
            other => panic!("expected ContigLengthOverflow, got {other:?}"),
        }
    }

    /// Mi5: empty `tool_string` skips the `##source` line; doc on
    /// `CohortMetadata::tool_string` contracts on this skip.
    #[test]
    fn header_omits_source_when_tool_string_empty() {
        let mut metadata = fixture_metadata();
        metadata.tool_string = String::new();
        let header = build_vcf_header(&metadata, &cfg_emit_gp_off()).unwrap();
        let text = header_to_string(&header);
        assert!(
            !text.contains("##source="),
            "##source should be absent with empty tool_string\n{text}"
        );
    }

    /// Mi5: empty `command_line` skips the `##commandline` line.
    #[test]
    fn header_omits_commandline_when_command_line_empty() {
        let mut metadata = fixture_metadata();
        metadata.command_line = String::new();
        let header = build_vcf_header(&metadata, &cfg_emit_gp_off()).unwrap();
        let text = header_to_string(&header);
        assert!(
            !text.contains("##commandline="),
            "##commandline should be absent with empty command_line\n{text}"
        );
    }

    /// Mi3: empty `contigs` is structurally valid; we accept it and
    /// produce a header with no `##contig=` lines. (Documented in
    /// `CohortMetadata::contigs`.)
    #[test]
    fn empty_contigs_accepted_produces_header_with_no_contig_lines() {
        let mut metadata = fixture_metadata();
        metadata.contigs.clear();
        let header = build_vcf_header(&metadata, &cfg_emit_gp_off()).unwrap();
        let text = header_to_string(&header);
        assert!(
            !text.contains("##contig="),
            "no ##contig lines expected\n{text}"
        );
        // Sample names still flow through.
        let chrom_line = text.lines().find(|l| l.starts_with("#CHROM")).unwrap();
        assert!(chrom_line.ends_with("\tS0\tS1\tS2"));
    }
}
