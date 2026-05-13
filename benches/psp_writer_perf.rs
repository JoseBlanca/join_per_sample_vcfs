//! `.psp` writer end-to-end throughput.
//!
//! One workload: a typical Stage-1 pileup stream — single chromosome,
//! mostly REF, occasional 0/1 SNP, no phase-chain markers. Drives the
//! writer through validation, per-column buffering, the auto-flush
//! cadence (16 MiB target), zstd compression at level 9, and the
//! framed I/O path.
//!
//! Bench shape: a single `Vec<PileupRecord>` is built once outside
//! the timed region; the bench body re-runs `new + write_record* +
//! finish` against an `io::sink()`, borrowing the records by
//! reference. This keeps the timed region (and a samply
//! `--profile-time` run) dominated by writer work rather than the
//! cost of cloning the fixture per iteration.
//!
//! Sized so a full bench iteration emits at least one auto-flushed
//! block (well past the 16 MiB target uncompressed) — the per-block
//! overhead (column encoding, zstd, manifest, header) is then on the
//! measured path. With ~3.3M records and a typical
//! ~50-byte-per-record uncompressed footprint we land ~6 blocks,
//! which keeps both the per-record append cost and the per-flush
//! cost in the measurement at a representative ratio.

use std::hint::black_box;
use std::io;
use std::time::Duration;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use merge_per_sample_vcfs::per_sample_caller::pileup::{
    AlleleObservation, AlleleSupportStats, PileupRecord,
};
use merge_per_sample_vcfs::per_sample_caller::psp::header::{
    ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance,
};
use merge_per_sample_vcfs::per_sample_caller::psp::writer::PspWriter;
use std::collections::BTreeMap;

/// Number of records per bench iteration. Picked so the writer emits
/// several full blocks (16 MiB uncompressed target each).
const NUM_RECORDS: usize = 3_300_000;

/// Probability (out of 1000) that a record carries a second allele.
/// 1-in-1000 ≈ a realistic SNP density at WGS coverage; the bench
/// stays dominated by single-allele REF positions, like a real run.
const SECOND_ALLELE_PER_1000: u32 = 1;

fn writer_header() -> WriterHeader {
    let mut params = BTreeMap::new();
    params.insert("min-mapq".to_string(), ParameterValue::Integer(30));
    WriterHeader {
        format_version: (1, 0),
        sample: "bench".to_string(),
        reference: "ref.fa".to_string(),
        created: "2026-05-13T11:00:00Z".parse().unwrap(),
        chromosomes: vec![ChromosomeEntry {
            name: "chr1".to_string(),
            length: NUM_RECORDS as u32 + 1024,
            md5: "0".repeat(32),
        }],
        writer: WriterProvenance {
            tool: "bench".to_string(),
            version: "0.0.0".to_string(),
            subcommand: "psp-writer-perf".to_string(),
            input_crams: vec!["a.cram".to_string()],
            input_fasta: "ref.fa".to_string(),
            parameters: params,
        },
    }
}

/// Build `NUM_RECORDS` SNP-typical records. Allele bases cycle through
/// {A, C, G, T} so the bytes column is not all the same value (which
/// would flatter zstd unrealistically). q-sum-log is a small negative
/// scalar; the per-allele scalars vary with `i` so the U32 columns
/// also carry compressible-but-not-trivial patterns.
fn build_records() -> Vec<PileupRecord> {
    let mut records = Vec::with_capacity(NUM_RECORDS);
    let bases = [b'A', b'C', b'G', b'T'];
    for i in 0..NUM_RECORDS {
        let pos = (i as u32) + 1;
        let ref_base = bases[i & 3];
        // The 1-in-1000 second-allele case: pick the next base in
        // the rotation so alt != ref.
        let two_alleles = (i as u32 % 1000) < SECOND_ALLELE_PER_1000;
        let mut alleles = Vec::with_capacity(if two_alleles { 2 } else { 1 });
        alleles.push(AlleleObservation::new(
            vec![ref_base],
            AlleleSupportStats::new(
                28 + ((i as u32) % 4),
                -42.0 - (i as f64).rem_euclid(7.0),
                14 + ((i as u32) % 3),
                5 + ((i as u32) % 3),
                1 + ((i as u32) % 3),
            ),
            Vec::new(),
        ));
        if two_alleles {
            let alt_base = bases[(i + 1) & 3];
            alleles.push(AlleleObservation::new(
                vec![alt_base],
                AlleleSupportStats::new(
                    3 + ((i as u32) % 5),
                    -3.5 - (i as f64).rem_euclid(2.0),
                    1 + ((i as u32) % 2),
                    1,
                    0,
                ),
                Vec::new(),
            ));
        }
        records.push(PileupRecord::new(0, pos, Vec::new(), Vec::new(), alleles));
    }
    records
}

fn write_all(records: &[PileupRecord], header: WriterHeader) -> u64 {
    let mut writer = PspWriter::new(io::sink(), header).expect("writer new");
    let mut flushed = 0u64;
    for r in records {
        flushed += writer.write_record(r).expect("write_record");
    }
    writer.finish().expect("finish");
    flushed
}

fn bench_writer(c: &mut Criterion) {
    let mut group = c.benchmark_group("psp_writer");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(30));

    // Build records and the header template once; the bench body
    // borrows `records` and clones only the small header. The
    // `BlockAccumulator` clones the inner Vecs (new_chains,
    // expired_chains, allele.chain_slots) as it appends, so those
    // allocations are still measured exactly as production would
    // produce them.
    let records = build_records();
    let header = writer_header();
    group.throughput(Throughput::Elements(records.len() as u64));

    group.bench_function("snp_typical_3_3M", |b| {
        b.iter(|| black_box(write_all(&records, header.clone())));
    });

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = bench_writer
}

criterion_main!(benches);
