//! Thread-budget contract guard (plan §6 of the thread-budget single-pool
//! plan). A real `var-calling` run must spawn **≤ N + c** live OS threads (a
//! small constant `c`) — not the historical `~3N`. The binary logs
//! `live_threads=<n>` once its full topology is up (`main` + producer pool +
//! caller threads + writer [+ staged coordinators]); this test runs the binary
//! and asserts the count, so a refactor that silently re-grows the budget
//! (e.g. re-introducing a per-stage `N`-sized pool) trips here.
//!
//! Linux-only: the count comes from `/proc/self/status`, which the binary reads
//! on Linux. The dev-container gate is Linux, so this runs there; on other
//! hosts the whole file compiles away.
#![cfg(target_os = "linux")]

mod common;

use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::process::Command;

use common::{CONTIG_LEN, CONTIG_NAME, build_fasta, fixture_md5};
use pop_var_caller::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};
use pop_var_caller::psp::header::{ChromosomeEntry, WriterHeader, WriterProvenance};
use pop_var_caller::psp::writer::PspWriter;
use tempfile::TempDir;

/// `.psp` header matching the all-`A` fixture FASTA from [`build_fasta`]: same
/// contig name/length and the matching MD5, so the CLI's FASTA-vs-`.psp` MD5
/// cross-check passes and the real pipeline runs.
fn writer_header(sample: &str) -> WriterHeader {
    WriterHeader {
        format_version: (1, 0),
        sample: sample.to_string(),
        reference: "ref.fa".to_string(), // basename of `build_fasta`'s output
        created: "2026-06-09T00:00:00Z".parse().unwrap(),
        chromosomes: vec![ChromosomeEntry {
            name: CONTIG_NAME.to_string(),
            length: CONTIG_LEN as u32,
            md5: fixture_md5().to_string(),
        }],
        writer: WriterProvenance {
            tool: "test".to_string(),
            version: "0".to_string(),
            subcommand: "thread-budget".to_string(),
            input_crams: vec!["x.cram".to_string()],
            input_fasta: "ref.fa".to_string(),
            command_line: String::new(),
            parameters: BTreeMap::new(),
        },
    }
}

/// High-confidence support for `num_obs` observations (mirrors the cohort perf
/// bench): per-obs ~Q40, balanced strand, MAPQ 60 with zero variance.
fn support(num_obs: u32) -> AlleleSupportStats {
    let n = num_obs as f64;
    AlleleSupportStats::new(
        num_obs,
        n * -4.0,
        num_obs / 2,
        num_obs / 2,
        0,
        num_obs * 60,
        u64::from(num_obs) * 3_600,
    )
}

/// Write one sample's `.psp`: REF `A` at every position (matches the all-`A`
/// fixture FASTA), with an ALT every 10th so the cohort has variants to call.
/// The data is incidental — the test only needs the run to reach the producer
/// loop where the thread count is logged.
fn write_psp(path: &Path, sample: &str) {
    let sink = BufWriter::new(File::create(path).expect("create psp"));
    let mut writer = PspWriter::new(sink, writer_header(sample)).expect("psp writer");
    for pos in 1..=(CONTIG_LEN as u32 - 10) {
        let mut alleles = vec![AlleleObservation::new(vec![b'A'], support(30), Vec::new())];
        if pos % 10 == 0 {
            alleles.push(AlleleObservation::new(vec![b'C'], support(14), Vec::new()));
        }
        writer
            .write_record(&PileupRecord::new(0, pos, alleles))
            .expect("write record");
    }
    writer
        .finish()
        .expect("finish psp")
        .into_inner()
        .expect("into_inner")
        .sync_all()
        .expect("sync");
}

#[test]
fn var_calling_thread_count_within_budget() {
    let threads = 8usize;
    let dir = TempDir::new().expect("tempdir");
    let fasta = build_fasta(dir.path());
    let psp: Vec<_> = (0..3)
        .map(|i| {
            let sample = format!("S{i}");
            let path = dir.path().join(format!("{sample}.psp"));
            write_psp(&path, &sample);
            path
        })
        .collect();
    let out = dir.path().join("out.vcf");

    // Run the real binary (cargo builds it for integration tests). A fresh
    // process means `/proc/self` counts exactly the pipeline's threads, with no
    // test-harness threads to subtract.
    let output = Command::new(env!("CARGO_BIN_EXE_pop_var_caller"))
        .arg("var-calling")
        .arg("--reference")
        .arg(&fasta)
        .arg("--output")
        .arg(&out)
        .arg("--threads")
        .arg(threads.to_string())
        .arg("--no-complexity-filter") // all-`A` ref would otherwise DUST-mask out
        .args(&psp)
        .output()
        .expect("spawn var-calling binary");
    assert!(
        output.status.success(),
        "var-calling failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    let live: usize = stderr
        .lines()
        .find_map(|l| l.strip_prefix("var-calling: live_threads="))
        .and_then(|n| n.trim().parse().ok())
        .unwrap_or_else(|| panic!("no `live_threads=` line in stderr:\n{stderr}"));

    // Budget: N worker threads (producer pool P + callers C, P+C=N) + main +
    // writer + at most two staged coordinators ⇒ ≤ N + 4. The lower bound (the
    // N worker threads must exist) guards against measuring a partial topology.
    let max = threads + 4;
    assert!(
        (threads..=max).contains(&live),
        "live_threads={live} is outside [{threads}, {max}] for --threads={threads}; \
         the budget contract is ≤ N + c (c ≤ 4)"
    );
    eprintln!("thread-budget: --threads={threads} live_threads={live} (≤ N+4 = {max})");
}
