//! Integration tests for `pop_var_caller pileup`.
//!
//! Each test builds a tiny synthetic FASTA + CRAM, calls `run_pileup`
//! programmatically, and re-opens the resulting `.psp` with
//! `PspReader` to verify the artefact round-trips.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use bstr::BString;
use merge_per_sample_vcfs::per_sample_caller::baq::{
    DEFAULT_BAQ_CHUNK_SIZE, SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH, SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
    SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
};
use merge_per_sample_vcfs::per_sample_caller::cram_input::{
    DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH,
    DEFAULT_MISMATCH_BQ_FLOOR,
};
use merge_per_sample_vcfs::per_sample_caller::pileup::{
    DEFAULT_MATE_LOOKUP_WINDOW, DEFAULT_MAX_ACTIVE_SLOTS, DEFAULT_MAX_INDEL_COLUMN_DEPTH,
    DEFAULT_MAX_RECORD_SPAN, DEFAULT_MAX_SNP_COLUMN_DEPTH,
};
use merge_per_sample_vcfs::per_sample_caller::psp::PspReader;
use merge_per_sample_vcfs::per_sample_caller::psp::header::ParameterValue;
use merge_per_sample_vcfs::pop_var_caller::{PileupArgs, run_pileup};
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
// Fixture helpers (self-contained — kept inline since `cram_files.rs`
// is `pub(crate)` test-only and not reachable from an integration
// test crate)
// ---------------------------------------------------------------------

const CONTIG_NAME: &str = "chr1";
const CONTIG_LEN: usize = 200;
// MD5 of an all-A 200-base string; doesn't need to match the FASTA
// content because Stage 1 doesn't cross-check (Stage 3 will, with
// a real FASTA). The pipeline only requires *some* 32-hex-char @SQ M5
// to flow through into the .psp header.
const FAKE_MD5: &str = "0123456789abcdef0123456789abcdef";

fn build_fasta(dir: &Path) -> PathBuf {
    let fasta_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let mut fa = File::create(&fasta_path).expect("create fasta");
    let mut fai = File::create(&fai_path).expect("create fai");
    writeln!(fa, ">{}", CONTIG_NAME).unwrap();
    let header_len = (CONTIG_NAME.len() + 2) as u64; // '>' + name + '\n'
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

fn build_sam_header(include_md5: bool) -> sam::Header {
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::SAMPLE;
    use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;

    let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
    hd.other_fields_mut()
        .insert(SORT_ORDER, BString::from(b"coordinate".as_ref()));

    let length_nz = NonZero::new(CONTIG_LEN).unwrap();
    let mut ref_seq = Map::<ReferenceSequence>::new(length_nz);
    if include_md5 {
        ref_seq
            .other_fields_mut()
            .insert(MD5_CHECKSUM, BString::from(FAKE_MD5.as_bytes()));
    }

    let mut rg = Map::<ReadGroup>::default();
    rg.other_fields_mut()
        .insert(SAMPLE, BString::from(b"NA12878".as_ref()));

    sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(CONTIG_NAME.as_bytes().to_vec(), ref_seq)
        .add_read_group(b"rg0".as_ref(), rg)
        .build()
}

fn build_cram(dir: &Path, fasta_path: &Path, include_md5: bool, records: &[RecordBuf]) -> PathBuf {
    let cram_path = dir.join("sample.cram");
    let header = build_sam_header(include_md5);

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

/// Build a synthetic single-end SNP read at `pos` (1-based), 5bp,
/// MAPQ 60, all-match CIGAR. `seq` is the read bases.
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

/// Default `PileupArgs` instance overridable by the test. We
/// disable the min-read-length filter (default 30) so the tiny
/// synthetic 5-bp reads survive into the pipeline.
fn default_args(reference: PathBuf, output: PathBuf, crams: Vec<PathBuf>) -> PileupArgs {
    let _unused = DEFAULT_MIN_READ_LENGTH; // re-export sanity
    PileupArgs {
        reference,
        output,
        threads: None,
        crams,
        min_mapq: DEFAULT_MIN_MAPQ,
        no_baq: false,
        min_read_length: 0, // admit any length
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
        max_active_reads: DEFAULT_MAX_ACTIVE_SLOTS,
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

/// Happy path: BAQ on, single CRAM with three reads covering
/// positions 1..14 — the walker should emit one record per covered
/// position. We verify:
///  - `run_pileup` returns Ok.
///  - The .psp opens cleanly and yields the expected records.
///  - WriterProvenance carries tool/subcommand/basenames and every
///    CLI knob shows up in `parameters`.
#[test]
fn happy_path_default_config() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let records = vec![
        read_record("r1", 1, b"AAAAA"),
        read_record("r2", 5, b"AAAAA"),
        read_record("r3", 10, b"AAAAA"),
    ];
    let cram = build_cram(dir.path(), &fasta, /*md5=*/ true, &records);
    let output = dir.path().join("out.psp");
    let args = default_args(fasta.clone(), output.clone(), vec![cram.clone()]);

    run_pileup(&args).expect("run_pileup OK");

    // No leftover .tmp once we returned successfully.
    assert!(!output.with_extension("psp.tmp").exists());
    assert!(output.exists(), ".psp must exist at the final path");

    // Round-trip through PspReader.
    let f = File::open(&output).expect("open .psp");
    let reader = PspReader::new(BufReader::new(f)).expect("open psp reader");
    let header = reader.header().clone();

    assert_eq!(header.sample, "NA12878");
    assert_eq!(header.writer.tool, "pop_var_caller");
    assert_eq!(header.writer.subcommand, "pileup");
    assert_eq!(header.writer.input_fasta, "ref.fa");
    assert_eq!(header.writer.input_crams, vec!["sample.cram"]);
    assert_eq!(header.chromosomes.len(), 1);
    assert_eq!(header.chromosomes[0].name, CONTIG_NAME);
    assert_eq!(header.chromosomes[0].md5, FAKE_MD5);

    // Every CLI knob got recorded.
    let params: &BTreeMap<String, ParameterValue> = &header.writer.parameters;
    for key in &[
        "min_mapq",
        "min_read_length",
        "drop_qc_fail",
        "drop_duplicate",
        "max_read_mismatch_fraction",
        "mismatch_bq_floor",
        "baq_enabled",
        "baq_gap_open_prob",
        "baq_gap_extend_prob",
        "baq_band_half_width",
        "baq_chunk_size",
        "max_snp_column_depth",
        "max_indel_column_depth",
        "max_record_span",
        "mate_lookup_window",
        "max_active_reads",
        "threads",
    ] {
        assert!(params.contains_key(*key), "missing param: {key}");
    }
    assert_eq!(params["baq_enabled"], ParameterValue::Boolean(true));

    // Count records emitted — the three reads at 1, 5, 10 over a 5bp
    // each produce records at every position they cover. Reads
    // overlap at 5..=9 (r1 covers 1..=5, r2 covers 5..=9, r3 covers
    // 10..=14), so positions 1..=14 should each have a record.
    let mut reader = PspReader::new(BufReader::new(File::open(&output).unwrap())).unwrap();
    let positions: Vec<u32> = reader
        .records()
        .map(|r| r.expect("decode record").pos)
        .collect();
    assert_eq!(
        positions,
        (1u32..=14).collect::<Vec<_>>(),
        "expected one record at each of positions 1..=14"
    );
}

/// `--no-baq` path: the pipeline runs with BAQ bypassed. We verify
/// `WriterProvenance.parameters["baq_enabled"] == Boolean(false)`
/// and that the .psp still round-trips.
#[test]
fn no_baq_path_writes_valid_psp() {
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let records = vec![
        read_record("r1", 1, b"AAAAA"),
        read_record("r2", 3, b"AAAAA"),
    ];
    let cram = build_cram(dir.path(), &fasta, true, &records);
    let output = dir.path().join("out.psp");
    let mut args = default_args(fasta, output.clone(), vec![cram]);
    args.no_baq = true;

    run_pileup(&args).expect("run_pileup OK");

    let f = File::open(&output).expect("open .psp");
    let reader = PspReader::new(BufReader::new(f)).expect("open psp reader");
    let header = reader.header();
    assert_eq!(
        header.writer.parameters["baq_enabled"],
        ParameterValue::Boolean(false)
    );

    let mut reader = PspReader::new(BufReader::new(File::open(&output).unwrap())).unwrap();
    let count = reader.records().count();
    // Records cover positions 1..=7 (r1 1..=5, r2 3..=7).
    assert_eq!(count, 7, "expect 7 records, got {count}");
}

// NOTE: the missing-@SQ M5 hard-error path is covered at the unit
// level by `build_writer_header_errors_on_missing_md5` in
// `src/pop_var_caller/cli.rs::tests`. We cannot easily reproduce
// the missing-M5 case at the integration level because noodles'
// CRAM writer auto-derives an M5 from the reference repository on
// write, so the produced CRAM always has @SQ M5 present even when
// we omit it from the SAM header. Building such a CRAM would
// require hand-rolling the CRAM bytes, which is out of scope.
