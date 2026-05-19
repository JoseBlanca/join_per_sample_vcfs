//! Integration tests for the cohort CLI subcommands: `var-calling`,
//! `estimate-contamination`, and `var-calling-from-bam`.
//!
//! Each test builds a tiny synthetic FASTA + CRAM(s), drives the
//! relevant `run_*` orchestrator(s), and verifies the resulting
//! artefacts. Fixture helpers mirror those in
//! `tests/pileup_cli_integration.rs` — duplicated here because
//! integration test crates can't import each other's helpers and the
//! per_sample_pileup test-only helpers are `pub(crate)`. Future
//! cleanup may share them via a `tests/common/` module.

use std::fs::{self, File};
use std::io::Write;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use bstr::BString;
use merge_per_sample_vcfs::per_sample_pileup::baq::{
    DEFAULT_BAQ_CHUNK_SIZE, SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH, SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
    SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
};
use merge_per_sample_vcfs::per_sample_pileup::cram_input::{
    DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MISMATCH_BQ_FLOOR,
};
use merge_per_sample_vcfs::per_sample_pileup::pileup::{
    DEFAULT_MATE_LOOKUP_WINDOW, DEFAULT_MAX_ACTIVE_READS, DEFAULT_MAX_INDEL_COLUMN_DEPTH,
    DEFAULT_MAX_RECORD_SPAN, DEFAULT_MAX_SNP_COLUMN_DEPTH,
};
use merge_per_sample_vcfs::pop_var_caller::cli::shared_args::{CohortPipelineArgs, Stage1Args};
use merge_per_sample_vcfs::pop_var_caller::var_calling::{
    VarCallingArgs, VarCallingCliError, run_var_calling,
};
use merge_per_sample_vcfs::pop_var_caller::var_calling_from_bam::{
    VarCallingFromBamArgs, run_var_calling_from_bam,
};
use merge_per_sample_vcfs::pop_var_caller::{PileupArgs, run_pileup};
use merge_per_sample_vcfs::var_calling::dust_filter::{
    DEFAULT_DUST_THRESHOLD, DEFAULT_DUST_WINDOW,
};
use merge_per_sample_vcfs::var_calling::per_group_merger::{
    DEFAULT_MAX_ALLELES_PER_RECORD, DEFAULT_PLOIDY,
};
use merge_per_sample_vcfs::var_calling::posterior_engine::{
    DEFAULT_COMPOUND_ALT_PSEUDOCOUNT, DEFAULT_CONVERGENCE_THRESHOLD,
    DEFAULT_INBREEDING_COEFFICIENT, DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MAX_GQ_PHRED,
    DEFAULT_MAX_ITERATIONS, DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT,
};
use merge_per_sample_vcfs::var_calling::variant_grouping::DEFAULT_MAX_VARIANT_GROUP_SPAN;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;
use sam::alignment::record::Flags;
use sam::alignment::record::MappingQuality;
use sam::alignment::record::cigar::Op;
use sam::alignment::record::cigar::op::Kind;
use sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence};
use sam::header::record::value::Map;
use sam::header::record::value::map::header::Version;
use sam::header::record::value::map::{Header as HeaderMap, ReadGroup, ReferenceSequence};
use tempfile::TempDir;

// ---------------------------------------------------------------------
// Fixture helpers (mirror tests/pileup_cli_integration.rs)
// ---------------------------------------------------------------------

const CONTIG_NAME: &str = "chr1";
const CONTIG_LEN: usize = 200;
const FAKE_MD5: &str = "0123456789abcdef0123456789abcdef";

fn build_fasta(dir: &Path) -> PathBuf {
    let fasta_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let mut fa = File::create(&fasta_path).expect("create fasta");
    let mut fai = File::create(&fai_path).expect("create fai");
    writeln!(fa, ">{}", CONTIG_NAME).unwrap();
    let header_len = (CONTIG_NAME.len() + 2) as u64;
    let seq = vec![b'A'; CONTIG_LEN];
    fa.write_all(&seq).unwrap();
    fa.write_all(b"\n").unwrap();
    writeln!(
        fai,
        "{}\t{}\t{}\t{}\t{}",
        CONTIG_NAME,
        CONTIG_LEN,
        header_len,
        CONTIG_LEN,
        CONTIG_LEN + 1
    )
    .unwrap();
    fasta_path
}

fn build_sam_header(sample: &str) -> sam::Header {
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::SAMPLE;
    use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;

    let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
    hd.other_fields_mut()
        .insert(SORT_ORDER, BString::from(b"coordinate".as_ref()));

    let length_nz = NonZero::new(CONTIG_LEN).unwrap();
    let mut ref_seq = Map::<ReferenceSequence>::new(length_nz);
    ref_seq
        .other_fields_mut()
        .insert(MD5_CHECKSUM, BString::from(FAKE_MD5.as_bytes()));

    let mut rg = Map::<ReadGroup>::default();
    rg.other_fields_mut()
        .insert(SAMPLE, BString::from(sample.as_bytes()));

    sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(CONTIG_NAME.as_bytes().to_vec(), ref_seq)
        .add_read_group(b"rg0".as_ref(), rg)
        .build()
}

fn build_cram(dir: &Path, fasta_path: &Path, sample: &str, records: &[RecordBuf]) -> PathBuf {
    let cram_path = dir.join(format!("{sample}.cram"));
    let header = build_sam_header(sample);

    let indexed_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta_path)
        .expect("indexed fasta reader");
    let adapter = fasta::repository::adapters::IndexedReader::new(indexed_reader);
    let repository = fasta::Repository::new(adapter);

    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(&cram_path)
        .expect("cram writer");
    writer.write_header(&header).expect("write header");
    for record in records {
        sam::alignment::io::Write::write_alignment_record(&mut writer, &header, record)
            .expect("write record");
    }
    sam::alignment::io::Write::finish(&mut writer, &header).expect("finish cram");
    cram_path
}

fn read_record(qname: &str, pos: usize, seq: &[u8]) -> RecordBuf {
    let mut rb = RecordBuf::default();
    *rb.name_mut() = Some(BString::from(qname.as_bytes()));
    *rb.flags_mut() = Flags::from(0u16);
    *rb.reference_sequence_id_mut() = Some(0);
    *rb.alignment_start_mut() = noodles_core::Position::new(pos);
    *rb.mapping_quality_mut() = MappingQuality::new(60);
    let cigar_ops = vec![Op::new(Kind::Match, seq.len())];
    *rb.cigar_mut() = Cigar::from(cigar_ops);
    *rb.sequence_mut() = Sequence::from(seq.to_vec());
    *rb.quality_scores_mut() = QualityScores::from(vec![30u8; seq.len()]);
    rb
}

fn pileup_args(reference: PathBuf, output: PathBuf, crams: Vec<PathBuf>) -> PileupArgs {
    PileupArgs {
        reference,
        output,
        threads: None,
        crams,
        stage1: Stage1Args {
            min_mapq: DEFAULT_MIN_MAPQ,
            no_baq: true, // skip BAQ — tiny synthetic reads don't need it
            min_read_length: 0,
            keep_qc_fail: false,
            keep_duplicates: false,
            max_read_mismatch_fraction: DEFAULT_MAX_READ_MISMATCH_FRACTION,
            mismatch_bq_floor: DEFAULT_MISMATCH_BQ_FLOOR,
            baq_gap_open_prob: SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
            baq_gap_extend_prob: SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
            baq_band_half_width: SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
            baq_chunk_size: DEFAULT_BAQ_CHUNK_SIZE,
            max_snp_column_depth: DEFAULT_MAX_SNP_COLUMN_DEPTH,
            max_indel_column_depth: DEFAULT_MAX_INDEL_COLUMN_DEPTH,
            max_record_span: DEFAULT_MAX_RECORD_SPAN,
            mate_lookup_window: DEFAULT_MATE_LOOKUP_WINDOW,
            max_active_reads: DEFAULT_MAX_ACTIVE_READS,
        },
    }
}

fn var_calling_args(
    reference: PathBuf,
    output: PathBuf,
    psp_files: Vec<PathBuf>,
    contamination_estimates: Option<PathBuf>,
) -> VarCallingArgs {
    VarCallingArgs {
        reference,
        output,
        threads: None,
        contamination_estimates,
        no_complexity_filter: true, // tiny ref; sdust would mask everything
        psp_files,
        cohort: CohortPipelineArgs {
            ploidy: DEFAULT_PLOIDY,
            complexity_window: DEFAULT_DUST_WINDOW,
            complexity_threshold: DEFAULT_DUST_THRESHOLD,
            var_group_max_span: DEFAULT_MAX_VARIANT_GROUP_SPAN,
            max_alleles_per_var: DEFAULT_MAX_ALLELES_PER_RECORD,
            inbreeding_coefficient: DEFAULT_INBREEDING_COEFFICIENT,
            em_convergence_threshold: DEFAULT_CONVERGENCE_THRESHOLD,
            em_max_iterations: DEFAULT_MAX_ITERATIONS,
            ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT,
            snp_alt_pseudocount: DEFAULT_SNP_ALT_PSEUDOCOUNT,
            indel_alt_pseudocount: DEFAULT_INDEL_ALT_PSEUDOCOUNT,
            compound_alt_pseudocount: DEFAULT_COMPOUND_ALT_PSEUDOCOUNT,
            max_gq_phred: DEFAULT_MAX_GQ_PHRED,
            emit_gp: false,
        },
    }
}

/// Build a `.psp` for one sample by running the real pileup
/// orchestrator on a synthetic CRAM. Returns the .psp path.
fn make_psp_for_sample(
    dir: &Path,
    fasta_path: &Path,
    sample: &str,
    reads: &[RecordBuf],
) -> PathBuf {
    let cram = build_cram(dir, fasta_path, sample, reads);
    let psp = dir.join(format!("{sample}.psp"));
    let args = pileup_args(fasta_path.to_path_buf(), psp.clone(), vec![cram]);
    run_pileup(&args).expect("run_pileup OK");
    psp
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

/// **End-to-end happy path.** Three samples → three `.psp` files →
/// run_var_calling → VCF on disk. Verifies that the whole cohort
/// pipeline glue (PspReader → merger → DustFilter [bypassed] →
/// grouper → per-group merger → posterior engine → VCF writer) wires
/// together without errors and produces a parseable VCF.
#[test]
fn var_calling_happy_path_three_samples() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());

    // Three samples, identical-ish reads covering positions 1..14.
    // Reads are mostly REF (all-A reference); some bases are altered
    // to introduce SNP variants that the grouper will pick up.
    let reads_a = vec![
        read_record("r1", 10, b"AAAAA"),
        read_record("r2", 15, b"AAAAA"),
    ];
    let reads_b = vec![
        read_record("r1", 10, b"ACAAA"), // SNP at pos 11
        read_record("r2", 15, b"AAAAA"),
    ];
    let reads_c = vec![
        read_record("r1", 10, b"AAAAA"),
        read_record("r2", 15, b"AAAAA"),
    ];

    let psp_a = make_psp_for_sample(dir.path(), &fasta, "NA00001", &reads_a);
    let psp_b = make_psp_for_sample(dir.path(), &fasta, "NA00002", &reads_b);
    let psp_c = make_psp_for_sample(dir.path(), &fasta, "NA00003", &reads_c);

    let vcf_out = dir.path().join("cohort.vcf");
    let args = var_calling_args(fasta, vcf_out.clone(), vec![psp_a, psp_b, psp_c], None);
    run_var_calling(&args).expect("run_var_calling OK");

    assert!(vcf_out.exists(), "VCF must exist at the final path");
    let body = fs::read_to_string(&vcf_out).expect("read VCF");
    // Self-checks on the header — we expect a multi-sample VCF.
    assert!(body.contains("##fileformat=VCFv4"), "VCF header missing");
    assert!(body.contains("##source="), "##source missing");
    assert!(
        body.contains("\tNA00001\tNA00002\tNA00003"),
        "sample columns missing"
    );
}

/// **var-calling-from-bam happy path.** One CRAM directly → VCF.
/// No `.psp` is written by the test (the helper internally writes
/// the .psp for the cohort path, but `run_var_calling_from_bam`
/// goes BAM → VCF in a single process per the
/// `feedback_no_silent_intermediates` rule).
#[test]
fn var_calling_from_bam_happy_path() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let reads = vec![
        read_record("r1", 10, b"ACAAA"),
        read_record("r2", 15, b"AAAAA"),
    ];
    let cram = build_cram(dir.path(), &fasta, "NA12878", &reads);
    let vcf_out = dir.path().join("solo.vcf");

    let args = VarCallingFromBamArgs {
        reference: fasta.clone(),
        output: vcf_out.clone(),
        threads: None,
        crams: vec![cram],
        no_complexity_filter: true,
        stage1: Stage1Args {
            min_mapq: DEFAULT_MIN_MAPQ,
            no_baq: true,
            min_read_length: 0,
            keep_qc_fail: false,
            keep_duplicates: false,
            max_read_mismatch_fraction: DEFAULT_MAX_READ_MISMATCH_FRACTION,
            mismatch_bq_floor: DEFAULT_MISMATCH_BQ_FLOOR,
            baq_gap_open_prob: SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
            baq_gap_extend_prob: SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
            baq_band_half_width: SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
            baq_chunk_size: DEFAULT_BAQ_CHUNK_SIZE,
            max_snp_column_depth: DEFAULT_MAX_SNP_COLUMN_DEPTH,
            max_indel_column_depth: DEFAULT_MAX_INDEL_COLUMN_DEPTH,
            max_record_span: DEFAULT_MAX_RECORD_SPAN,
            mate_lookup_window: DEFAULT_MATE_LOOKUP_WINDOW,
            max_active_reads: DEFAULT_MAX_ACTIVE_READS,
        },
        cohort: CohortPipelineArgs {
            ploidy: DEFAULT_PLOIDY,
            complexity_window: DEFAULT_DUST_WINDOW,
            complexity_threshold: DEFAULT_DUST_THRESHOLD,
            var_group_max_span: DEFAULT_MAX_VARIANT_GROUP_SPAN,
            max_alleles_per_var: DEFAULT_MAX_ALLELES_PER_RECORD,
            inbreeding_coefficient: DEFAULT_INBREEDING_COEFFICIENT,
            em_convergence_threshold: DEFAULT_CONVERGENCE_THRESHOLD,
            em_max_iterations: DEFAULT_MAX_ITERATIONS,
            ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT,
            snp_alt_pseudocount: DEFAULT_SNP_ALT_PSEUDOCOUNT,
            indel_alt_pseudocount: DEFAULT_INDEL_ALT_PSEUDOCOUNT,
            compound_alt_pseudocount: DEFAULT_COMPOUND_ALT_PSEUDOCOUNT,
            max_gq_phred: DEFAULT_MAX_GQ_PHRED,
            emit_gp: false,
        },
    };
    run_var_calling_from_bam(&args).expect("run_var_calling_from_bam OK");

    assert!(vcf_out.exists(), "VCF must exist at the final path");
    let body = fs::read_to_string(&vcf_out).expect("read VCF");
    assert!(body.contains("##fileformat=VCFv4"), "VCF header missing");
    // Single-sample VCF — only NA12878 in the column header.
    assert!(body.contains("\tNA12878\n") || body.contains("\tNA12878\r\n"));
    // No other sample names should appear in the data-header line.
    assert!(!body.contains("NA00001"));
}

/// **Sample-name reconciliation: missing sample is a hard error.**
/// Build a `.psp` for `NA00001` but hand var-calling a contamination
/// artefact that only references `OTHER_SAMPLE`. `run_var_calling`
/// must surface `ContaminationSampleMissing`.
#[test]
fn var_calling_rejects_contamination_artefact_missing_sample() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let reads = vec![read_record("r1", 10, b"AAAAA")];
    let psp = make_psp_for_sample(dir.path(), &fasta, "NA00001", &reads);

    // Hand-craft a contamination artefact missing NA00001. The
    // on-disk TOML is the realistic input path here: `run_var_calling`
    // always loads via `ContaminationArtefact::read`, and writing a
    // literal mirrors what a user editing the file by hand would
    // produce. Avoids depending on the (`#[non_exhaustive]`) artefact
    // struct surface from outside the crate.
    let artefact_toml = r#"[provenance]
tool = "pop_var_caller"
version = "0.1.0"
subcommand = "estimate-contamination"
created = 2026-05-19T12:00:00Z

[provenance.inputs]
reference = "ref.fa"
input_psps = ["OTHER_SAMPLE.psp"]

[parameters]

[[batches]]
id = "lane_1"
contaminant_ref_prob = 0.999
contaminant_snp_alt_prob = 0.0008
contaminant_indel_alt_prob = 0.0002

[[samples]]
name = "OTHER_SAMPLE"
batch = "lane_1"
contamination_fraction = 0.01
"#;
    let artefact_path = dir.path().join("estimates.toml");
    fs::write(&artefact_path, artefact_toml).expect("write artefact");

    let vcf_out = dir.path().join("out.vcf");
    let args = var_calling_args(fasta, vcf_out, vec![psp], Some(artefact_path));
    let err = run_var_calling(&args).expect_err("must error on missing sample");
    let chain = format!("{err}");
    let nested = match &err {
        VarCallingCliError::ContamArtefact(inner) => format!("{inner}"),
        _ => chain.clone(),
    };
    assert!(
        nested.contains("NA00001") || chain.contains("NA00001"),
        "error should mention the missing sample: {chain} / {nested}"
    );
}

/// **M13 multi-invocation pattern.** Two consecutive `run_pileup`
/// calls in the same process must both succeed — the cohort CLI's
/// rayon-pool gate
/// ([`configure_rayon_pool`](merge_per_sample_vcfs::pop_var_caller))
/// is idempotent, so library consumers and a future
/// multi-subcommand test runner can chain `run_*` helpers freely.
///
/// `args.threads = None` so this exercises the multi-invocation
/// *pattern* without depending on rayon's process-global state
/// being uninitialised by an earlier test (which would be racy
/// under cargo-test's shared-process model). The
/// `configure_rayon_pool(None)` path returns `Ok(())` without
/// touching `build_global()`, so both calls are trivially clean.
///
/// `#[serial]` keeps this test from interleaving with other tests
/// in the same binary that might exercise rayon. The unit test
/// `configure_rayon_pool_none_is_always_ok` in
/// `pop_var_caller::common::tests` covers the helper's
/// `None`-path contract directly.
#[serial_test::serial]
#[test]
fn run_pileup_can_be_called_back_to_back() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let reads = vec![read_record("r1", 10, b"AAAAA")];
    let cram_a = build_cram(dir.path(), &fasta, "A", &reads);
    let cram_b = build_cram(dir.path(), &fasta, "B", &reads);

    let args1 = pileup_args(fasta.clone(), dir.path().join("a.psp"), vec![cram_a]);
    run_pileup(&args1).expect("first run_pileup OK");

    let args2 = pileup_args(fasta, dir.path().join("b.psp"), vec![cram_b]);
    run_pileup(&args2).expect("second run_pileup OK");

    assert!(dir.path().join("a.psp").exists());
    assert!(dir.path().join("b.psp").exists());
}
