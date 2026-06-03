//! Integration tests for `pop_var_caller pileup`.
//!
//! Each test builds a tiny synthetic FASTA + CRAM, calls `run_pileup`
//! programmatically, and re-opens the resulting `.psp` with
//! `PspReader` to verify the artefact round-trips. Shared fixture
//! helpers come from `tests/common/mod.rs` (Mi20).

mod common;

use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use common::{
    CONTIG_NAME, build_cram, build_cram_with_len, build_fasta, build_fasta_with_len, fixture_md5,
    read_record,
};
use pop_var_caller::bam::alignment_input::{
    DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH,
    DEFAULT_MISMATCH_BQ_FLOOR,
};
use pop_var_caller::baq::{
    SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH, SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
    SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
};
use pop_var_caller::pileup::per_sample::baq_stream::DEFAULT_BAQ_CHUNK_SIZE;
use pop_var_caller::pileup::walker::{
    DEFAULT_MATE_LOOKUP_WINDOW, DEFAULT_MAX_ACTIVE_READS, DEFAULT_MAX_INDEL_COLUMN_DEPTH,
    DEFAULT_MAX_RECORD_SPAN, DEFAULT_MAX_SNP_COLUMN_DEPTH,
};
use pop_var_caller::pop_var_caller::cli::shared_args::Stage1Args;
use pop_var_caller::pop_var_caller::{PileupArgs, run_pileup};
use pop_var_caller::psp::PspReader;
use pop_var_caller::psp::header::ParameterValue;
use pop_var_caller::psp::writer::{
    DEFAULT_BLOCK_WINDOW_BP, MAX_BLOCK_TARGET_BYTES, MIN_BLOCK_TARGET_BYTES, TARGET_BLOCK_BYTES,
};
use tempfile::TempDir;

/// Default `PileupArgs` instance overridable by the test. We
/// disable the min-read-length filter (default 30) so the tiny
/// synthetic 5-bp reads survive into the pipeline.
fn default_args(reference: PathBuf, output: PathBuf, alignment_files: Vec<PathBuf>) -> PileupArgs {
    let _unused = DEFAULT_MIN_READ_LENGTH; // re-export sanity
    PileupArgs {
        reference,
        output,
        threads: None,
        alignment_files,
        block_target_bytes: TARGET_BLOCK_BYTES,
        block_window_bp: DEFAULT_BLOCK_WINDOW_BP,
        stage1: Stage1Args {
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
            max_active_reads: DEFAULT_MAX_ACTIVE_READS,
        },
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
    let cram = build_cram(dir.path(), &fasta, "NA12878", Some(fixture_md5()), &records);
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
    assert_eq!(header.writer.input_crams, vec!["NA12878.cram"]);
    assert_eq!(header.chromosomes.len(), 1);
    assert_eq!(header.chromosomes[0].name, CONTIG_NAME);
    assert_eq!(header.chromosomes[0].md5, fixture_md5());

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
        "block_target_bytes",
    ] {
        assert!(params.contains_key(*key), "missing param: {key}");
    }
    assert_eq!(params["baq_enabled"], ParameterValue::Boolean(true));
    assert_eq!(
        params["block_target_bytes"],
        ParameterValue::Integer(TARGET_BLOCK_BYTES as i64)
    );

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
    let cram = build_cram(dir.path(), &fasta, "NA12878", Some(fixture_md5()), &records);
    let output = dir.path().join("out.psp");
    let mut args = default_args(fasta, output.clone(), vec![cram]);
    args.stage1.no_baq = true;

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

// ---------------------------------------------------------------------
// BAM siblings (commit 6 of the BAM-input plan,
// `doc/devel/implementation_plans/bam_input_support.md`).
// ---------------------------------------------------------------------

/// BAM analogue of [`happy_path_default_config`]. Same record set,
/// same expected `.psp` records, with `.bam` as the input format.
/// Verifies the format-dispatch in `AlignmentMergedReader::new`
/// reaches `open_bam_record_stream` and the rest of the Stage-1
/// pipeline runs through unmodified.
#[test]
fn happy_path_default_config_bam() {
    use common::build_bam;

    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let records = vec![
        read_record("r1", 1, b"AAAAA"),
        read_record("r2", 5, b"AAAAA"),
        read_record("r3", 10, b"AAAAA"),
    ];
    let bam = build_bam(dir.path(), "NA12878", Some(fixture_md5()), &records);
    let output = dir.path().join("out.psp");
    let args = default_args(fasta.clone(), output.clone(), vec![bam.clone()]);

    run_pileup(&args).expect("run_pileup OK on bam");

    assert!(output.exists(), ".psp must exist at the final path");

    let f = File::open(&output).expect("open .psp");
    let reader = PspReader::new(BufReader::new(f)).expect("open psp reader");
    let header = reader.header().clone();

    assert_eq!(header.sample, "NA12878");
    // The serialized `input_crams` field name is kept for backward
    // compatibility with existing `.psp` files; for BAM inputs the
    // basename is just the .bam one.
    assert_eq!(header.writer.input_crams, vec!["NA12878.bam"]);
    assert_eq!(header.chromosomes.len(), 1);
    assert_eq!(header.chromosomes[0].name, CONTIG_NAME);

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

/// Mixed CRAM + BAM in one invocation is rejected before any reader
/// opens. The error message names both offending paths.
#[test]
fn pileup_rejects_mixed_cram_and_bam() {
    use common::build_bam;
    use pop_var_caller::bam::errors::AlignmentInputError;
    use pop_var_caller::pop_var_caller::cli::PileupCliError;

    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let reads = vec![read_record("r1", 1, b"AAAAA")];
    let cram = build_cram(dir.path(), &fasta, "NA00001", Some(fixture_md5()), &reads);
    let bam = build_bam(dir.path(), "NA00002", Some(fixture_md5()), &reads);
    let output = dir.path().join("mixed.psp");

    let args = default_args(fasta, output, vec![cram.clone(), bam.clone()]);
    let err = run_pileup(&args).expect_err("mixed cram + bam must error");

    match err {
        PileupCliError::AlignmentInput(AlignmentInputError::MixedAlignmentFileFormats {
            first_path,
            first_format,
            other_path,
            other_format,
        }) => {
            assert_eq!(first_path, cram);
            assert_eq!(first_format, "CRAM");
            assert_eq!(other_path, bam);
            assert_eq!(other_format, "BAM");
        }
        other => panic!("unexpected error variant: {other:?}"),
    }
}

/// M10 regression: `pileup` does not run through preflight, so the
/// `AlignmentMergedReader::new`-side `UnsupportedExtension` arm is
/// the only gate against an input whose extension is neither
/// `.cram` nor `.bam`. A future refactor that dropped the
/// classify pre-pass would let the open-pass's `.unwrap()` panic
/// instead of surfacing the typed error.
#[test]
fn pileup_rejects_input_with_unknown_extension() {
    use pop_var_caller::bam::errors::AlignmentInputError;
    use pop_var_caller::pop_var_caller::cli::PileupCliError;

    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let unknown = dir.path().join("sample.sam");
    std::fs::write(&unknown, b"x").expect("write placeholder");
    let output = dir.path().join("out.psp");
    let args = default_args(fasta, output, vec![unknown.clone()]);

    let err = run_pileup(&args).expect_err("must reject .sam");
    match err {
        PileupCliError::AlignmentInput(AlignmentInputError::UnsupportedExtension { path }) => {
            assert_eq!(path, unknown);
        }
        other => panic!("unexpected error variant: {other:?}"),
    }
}

/// `--block-target-bytes` threads to the writer and changes the
/// emitted .psp's block layout. We run the same large CRAM fixture
/// with the minimum (16 KiB) and the maximum (16 MiB) accepted
/// values; the small-target run must produce strictly more blocks
/// and a strictly larger file. The parameter map round-trips the
/// chosen value verbatim in each case (so the .psp self-describes
/// the writer setting that produced it).
#[test]
fn block_target_bytes_changes_block_layout_end_to_end() {
    // Contig long enough that the writer flushes several times at
    // the small target. ~40 bytes/record × 4000 records ≈ 160 KB
    // uncompressed → ~10 blocks at 16 KiB target, 1 block at 16 MiB.
    const CONTIG_LEN: usize = 4000;
    const READ_LEN: usize = 5;
    let dir = TempDir::new().expect("tempdir");
    let (fasta, md5) = build_fasta_with_len(dir.path(), CONTIG_LEN);

    // One read per starting position, length 5 — guarantees a
    // pileup record at every covered base (1..=CONTIG_LEN).
    let mut records = Vec::with_capacity(CONTIG_LEN - READ_LEN + 1);
    for i in 1..=(CONTIG_LEN - READ_LEN + 1) {
        records.push(read_record(&format!("r{i}"), i, b"AAAAA"));
    }
    let cram = build_cram_with_len(
        dir.path(),
        &fasta,
        "NA00001",
        Some(md5.as_str()),
        &records,
        CONTIG_LEN,
    );

    let run = |target: usize, suffix: &str| -> (usize, u64) {
        let output = dir.path().join(format!("out_{suffix}.psp"));
        let mut args = default_args(fasta.clone(), output.clone(), vec![cram.clone()]);
        args.block_target_bytes = target;
        run_pileup(&args).expect("run_pileup OK");

        let reader =
            PspReader::new(BufReader::new(File::open(&output).unwrap())).expect("open psp reader");
        // Self-describing parameter round-trips verbatim.
        assert_eq!(
            reader.header().writer.parameters["block_target_bytes"],
            ParameterValue::Integer(target as i64),
            "parameter map must record the chosen block target"
        );
        let n_blocks = reader.block_index().len();
        let on_disk_size = std::fs::metadata(&output).unwrap().len();
        (n_blocks, on_disk_size)
    };

    let (small_blocks, small_bytes) = run(MIN_BLOCK_TARGET_BYTES, "min");
    let (large_blocks, large_bytes) = run(MAX_BLOCK_TARGET_BYTES, "max");

    assert!(
        small_blocks > large_blocks,
        "smaller block target must produce more blocks: small={small_blocks} large={large_blocks}"
    );
    assert_eq!(
        large_blocks, 1,
        "16 MiB target should fit ~160 KB of records in a single block, got {large_blocks}"
    );
    assert!(
        small_bytes > large_bytes,
        "smaller blocks → more zstd-frame overhead → larger file: \
         small={small_bytes} bytes, large={large_bytes} bytes"
    );
}
