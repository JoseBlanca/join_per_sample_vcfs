//! Integration tests for `pop_var_caller psp-to-pileup`.
//!
//! Each test writes a tiny synthetic `.psp` file to a tempdir using
//! `PspWriter` directly (no CRAM round-trip — this is a read-side
//! tool), then calls `run_psp_to_pileup` with an output path in the
//! same tempdir and asserts the produced text line-by-line.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::BufWriter;
use std::path::Path;

use pop_var_caller::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};
use pop_var_caller::pop_var_caller::{PspToPileupArgs, run_psp_to_pileup};
use pop_var_caller::psp::header::{ChromosomeEntry, WriterHeader, WriterProvenance};
use pop_var_caller::psp::writer::PspWriter;
use tempfile::TempDir;

// ---------------------------------------------------------------------
// Fixture helpers
// ---------------------------------------------------------------------

fn writer_header_for(chrom_names: &[&str]) -> WriterHeader {
    WriterHeader {
        format_version: (1, 0),
        sample: "S1".to_string(),
        reference: "ref.fa".to_string(),
        created: "2026-05-15T00:00:00Z".parse().unwrap(),
        chromosomes: chrom_names
            .iter()
            .map(|n| ChromosomeEntry {
                name: (*n).to_string(),
                length: 1_000_000,
                md5: "0".repeat(32),
            })
            .collect(),
        writer: WriterProvenance {
            tool: "pop_var_caller".to_string(),
            version: "0.0.0-test".to_string(),
            subcommand: "pileup".to_string(),
            input_crams: vec!["sample.cram".to_string()],
            input_fasta: "ref.fa".to_string(),
            command_line: String::new(),
            parameters: BTreeMap::new(),
        },
    }
}

fn supp(num_obs: u32, fwd: u32) -> AlleleSupportStats {
    AlleleSupportStats::new(num_obs, -1.234, fwd, 0, 0, 0, 0)
}

fn write_psp(path: &Path, chrom_names: &[&str], records: &[PileupRecord]) {
    let file = File::create(path).expect("create psp");
    let sink = BufWriter::with_capacity(64 * 1024, file);
    let mut writer = PspWriter::new(sink, writer_header_for(chrom_names)).expect("psp writer");
    for record in records {
        writer.write_record(record).expect("write record");
    }
    let buf = writer.finish().expect("finish psp");
    // BufWriter::into_inner flushes; we don't need to sync explicitly
    // because the tempdir's lifetime ensures we read the file before it
    // disappears.
    buf.into_inner()
        .expect("into_inner")
        .sync_all()
        .expect("sync");
}

fn run_dump(input: &Path, output: &Path, region: Option<&str>, show_chains: bool) {
    let args = PspToPileupArgs {
        input: input.to_path_buf(),
        output: output.to_path_buf(),
        region: region.map(str::to_string),
        show_chain_ids: show_chains,
    };
    run_psp_to_pileup(&args).expect("run_psp_to_pileup");
}

fn read_lines(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .expect("read output")
        .lines()
        .map(str::to_string)
        .collect()
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[test]
fn streams_every_record_with_seven_columns_per_line() {
    let dir = TempDir::new().expect("tempdir");
    let psp = dir.path().join("sample.psp");
    let txt = dir.path().join("out.pileup");

    // Three records on chr1: a REF-only SNP position, an insertion-
    // bearing position, and a deletion-bearing position.
    let records = vec![
        PileupRecord::new(
            0,
            10,
            vec![AlleleObservation::new(b"A".to_vec(), supp(5, 3), vec![])],
        ),
        PileupRecord::new(
            0,
            20,
            vec![
                AlleleObservation::new(b"C".to_vec(), supp(0, 0), vec![]),
                AlleleObservation::new(b"CTG".to_vec(), supp(2, 2), vec![]),
            ],
        ),
        PileupRecord::new(
            0,
            30,
            vec![
                AlleleObservation::new(b"GAT".to_vec(), supp(1, 0), vec![]),
                AlleleObservation::new(b"G".to_vec(), supp(3, 2), vec![]),
            ],
        ),
    ];

    write_psp(&psp, &["chr1"], &records);
    run_dump(&psp, &txt, None, false);

    let lines = read_lines(&txt);
    assert_eq!(
        lines.len(),
        3,
        "expect three output lines for three records"
    );
    for line in &lines {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 7, "every line has 7 tab-separated columns");
    }

    // Per-record exact-content checks.
    let cols0: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(cols0[0..4], ["chr1", "10", "A", "5"]);
    assert_eq!(cols0[4], "...,,"); // 3 fwd dots, 2 rev commas
    assert_eq!(cols0[5], "!!!!!"); // Phred 0 placeholder × 5
    assert_eq!(cols0[6], "A:5:3:0:0:-1.234");

    let cols1: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(cols1[0..4], ["chr1", "20", "C", "2"]);
    // Insertion of "TG" — 2 forward reads, no reverse.
    assert_eq!(cols1[4], ".+2TG.+2TG");
    let allele_field = cols1[6];
    assert!(allele_field.starts_with("C:0:0:0:0:"));
    assert!(allele_field.contains(",CTG:2:2:0:0:"));

    let cols2: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(cols2[0..4], ["chr1", "30", "G", "4"]);
    // REF allele "GAT" (1 reverse), then DEL allele "G" (2 fwd, 1 rev),
    // dropping "AT" of the ref span.
    assert_eq!(cols2[4], ",.-2AT.-2AT,-2at");
    assert!(cols2[6].starts_with("GAT:1:0:0:0:"));
    assert!(cols2[6].contains(",G:3:2:0:0:"));
}

#[test]
fn region_flag_clamps_to_one_record() {
    let dir = TempDir::new().expect("tempdir");
    let psp = dir.path().join("sample.psp");
    let txt = dir.path().join("out.pileup");

    let records = vec![
        PileupRecord::new(
            0,
            10,
            vec![AlleleObservation::new(b"A".to_vec(), supp(3, 3), vec![])],
        ),
        PileupRecord::new(
            0,
            20,
            vec![AlleleObservation::new(b"C".to_vec(), supp(4, 2), vec![])],
        ),
        PileupRecord::new(
            0,
            30,
            vec![AlleleObservation::new(b"G".to_vec(), supp(5, 5), vec![])],
        ),
    ];

    write_psp(&psp, &["chr1"], &records);
    run_dump(&psp, &txt, Some("chr1:20-20"), false);

    let lines = read_lines(&txt);
    assert_eq!(lines.len(), 1, "region clamp keeps exactly one record");
    let cols: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(cols[1], "20");
    assert_eq!(cols[2], "C");
}

#[test]
fn unknown_chromosome_in_region_is_a_hard_error() {
    let dir = TempDir::new().expect("tempdir");
    let psp = dir.path().join("sample.psp");
    let txt = dir.path().join("out.pileup");

    let records = vec![PileupRecord::new(
        0,
        10,
        vec![AlleleObservation::new(b"A".to_vec(), supp(1, 1), vec![])],
    )];
    write_psp(&psp, &["chr1"], &records);

    let args = PspToPileupArgs {
        input: psp,
        output: txt,
        region: Some("nonexistent_chrom".into()),
        show_chain_ids: false,
    };
    let err = run_psp_to_pileup(&args).expect_err("should reject unknown chrom");
    let msg = format!("{err}");
    assert!(
        msg.contains("nonexistent_chrom"),
        "error mentions missing chrom: {msg}",
    );
}

#[test]
fn show_chain_ids_appends_chain_field() {
    let dir = TempDir::new().expect("tempdir");
    let psp = dir.path().join("sample.psp");
    let txt = dir.path().join("out.pileup");

    let records = vec![PileupRecord::new(
        0,
        100,
        vec![AlleleObservation::new(
            b"T".to_vec(),
            supp(3, 2),
            vec![5, 7, 9],
        )],
    )];
    write_psp(&psp, &["chr1"], &records);
    run_dump(&psp, &txt, None, true);

    let lines = read_lines(&txt);
    let cols: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(cols[6], "T:3:2:0:0:-1.234:5;7;9");
}
