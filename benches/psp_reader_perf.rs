//! `.psp` reader end-to-end throughput.
//!
//! Mirrors the [`psp_writer_perf`] bench shape: build the same three
//! workloads (SNP-typical, phase-chain-heavy, multi-allele), serialise
//! each to a `Vec<u8>` once outside the timed region, then time the
//! [`PspReader`] open + `records()` walk over the bytes.
//!
//! Workloads:
//!
//! - `psp_reader/snp_typical_3_3M` — single chromosome, mostly REF,
//!   ~0.1 % SNPs, no phase-chain markers. Dense common case; the same
//!   shape the writer's `snp_typical_3_3M` covers.
//! - `psp_reader/phase_chain_heavy_1M` — sustained ~6-active phase set,
//!   each record's first allele carries the active set as
//!   `chain_ids`. Drives the list column through `decode_list_column`
//!   and the `mem::take` materialisation in `materialise_next_record`.
//! - `psp_reader/multi_allele_500k` — 2–4 alleles per record, indel-
//!   shaped allele lengths from 1–10 bytes. Drives the bytes column
//!   (`decode_bytes_split`) and the per-allele scalar columns.
//! - `psp_reader/region_window_chr1_mid` — over the SNP file, 100 k
//!   `region_records(chr1, mid-50k, mid+50k)` queries against the
//!   binary-search + per-block snapshot path.

// Opt-in mimalloc global allocator (cargo bench --features alloc-mimalloc ...).
#[cfg(feature = "alloc-mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::hint::black_box;
use std::io::Cursor;
use std::time::Duration;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use pop_var_caller::per_sample_pileup::pileup::{
    AlleleObservation, AlleleSupportStats, PileupRecord,
};
use pop_var_caller::per_sample_pileup::psp::header::{
    ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance,
};
use pop_var_caller::per_sample_pileup::psp::{PspReader, writer::PspWriter};
use std::collections::BTreeMap;

const NUM_RECORDS_SNP: usize = 3_300_000;
const NUM_RECORDS_PHASE: usize = 1_000_000;
const NUM_RECORDS_MULTI: usize = 500_000;

/// Probability (out of 1000) that a SNP-workload record carries a
/// second allele. ~0.1 % matches realistic WGS SNP density.
const SECOND_ALLELE_PER_1000: u32 = 1;

fn writer_header(n_records: usize, sample: &str) -> WriterHeader {
    let mut params = BTreeMap::new();
    params.insert("min-mapq".to_string(), ParameterValue::Integer(30));
    WriterHeader {
        format_version: (1, 0),
        sample: sample.to_string(),
        reference: "ref.fa".to_string(),
        created: "2026-05-13T11:00:00Z".parse().unwrap(),
        chromosomes: vec![ChromosomeEntry {
            name: "chr1".to_string(),
            length: n_records as u32 + 1024,
            md5: "0".repeat(32),
        }],
        writer: WriterProvenance {
            tool: "bench".to_string(),
            version: "0.0.0".to_string(),
            subcommand: "psp-reader-perf".to_string(),
            input_crams: vec!["a.cram".to_string()],
            input_fasta: "ref.fa".to_string(),
            parameters: params,
        },
    }
}

fn build_snp_records(n: usize) -> Vec<PileupRecord> {
    let mut records = Vec::with_capacity(n);
    let bases = [b'A', b'C', b'G', b'T'];
    for i in 0..n {
        let pos = (i as u32) + 1;
        let ref_base = bases[i & 3];
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
        records.push(PileupRecord::new(0, pos, alleles));
    }
    records
}

fn build_phase_chain_heavy_records(n: usize) -> Vec<PileupRecord> {
    let mut records = Vec::with_capacity(n);
    let bases = [b'A', b'C', b'G', b'T'];
    let mut next_id: u64 = 0;
    let mut active: Vec<u64> = Vec::new();

    for i in 0..n {
        let pos = (i as u32) + 1;
        let ref_base = bases[i & 3];

        if i == 0 {
            for _ in 0..4 {
                active.push(next_id);
                next_id += 1;
            }
        } else {
            if i % 13 == 0 && active.len() > 3 {
                active.remove(0);
            }
            if i % 11 == 0 && active.len() < 12 {
                active.push(next_id);
                next_id += 1;
            }
        }

        let chain_ids = active.clone();
        let alleles = vec![AlleleObservation::new(
            vec![ref_base],
            AlleleSupportStats::new(20, -40.0, 10, 5, 2),
            chain_ids,
        )];

        records.push(PileupRecord::new(0, pos, alleles));
    }
    records
}

fn build_multi_allele_records(n: usize) -> Vec<PileupRecord> {
    let mut records = Vec::with_capacity(n);
    let bases = [b'A', b'C', b'G', b'T'];
    for i in 0..n {
        let pos = (i as u32) + 1;
        let r = bases[i & 3];
        let alleles = match i & 3 {
            0 => vec![
                AlleleObservation::new(
                    vec![r],
                    AlleleSupportStats::new(25, -38.0, 12, 4, 2),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    vec![bases[(i + 1) & 3]],
                    AlleleSupportStats::new(7, -10.0, 3, 1, 1),
                    Vec::new(),
                ),
            ],
            1 => vec![
                AlleleObservation::new(
                    vec![r],
                    AlleleSupportStats::new(22, -36.0, 11, 3, 1),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    vec![r, b'A', b'C'],
                    AlleleSupportStats::new(5, -8.0, 2, 1, 0),
                    Vec::new(),
                ),
            ],
            2 => vec![
                AlleleObservation::new(
                    vec![
                        r,
                        bases[(i + 1) & 3],
                        bases[(i + 2) & 3],
                        bases[(i + 3) & 3],
                        r,
                    ],
                    AlleleSupportStats::new(18, -30.0, 9, 2, 1),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    vec![r],
                    AlleleSupportStats::new(4, -7.0, 2, 0, 0),
                    Vec::new(),
                ),
            ],
            _ => vec![
                AlleleObservation::new(
                    vec![r],
                    AlleleSupportStats::new(15, -25.0, 7, 2, 1),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    vec![r, b'G', b'T'],
                    AlleleSupportStats::new(6, -9.0, 3, 1, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    vec![bases[(i + 1) & 3]],
                    AlleleSupportStats::new(3, -5.0, 1, 0, 0),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    vec![r, b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T', b'N'],
                    AlleleSupportStats::new(2, -4.0, 1, 0, 0),
                    Vec::new(),
                ),
            ],
        };
        records.push(PileupRecord::new(0, pos, alleles));
    }
    records
}

/// Serialise `records` to a `Vec<u8>` via the writer. Done once per
/// workload, outside the timed region — the timed body decodes those
/// bytes through `PspReader`.
fn serialise(records: &[PileupRecord], header: WriterHeader) -> Vec<u8> {
    let sink: Vec<u8> = Vec::with_capacity(8 * 1024 * 1024);
    let mut writer = PspWriter::new(sink, header).expect("writer new");
    for r in records {
        writer.write_record(r).expect("write_record");
    }
    writer.finish().expect("finish")
}

/// Drain every record from the source into a sink count, returning
/// the number of records successfully decoded. The sink is `black_box`
/// per-record so LLVM cannot dead-code the column `mem::take`s.
fn read_all_count(bytes: &[u8]) -> u64 {
    let mut reader = PspReader::new(Cursor::new(bytes)).expect("reader new");
    let mut n = 0u64;
    for r in reader.records() {
        let r = r.expect("record decode");
        black_box(&r);
        n += 1;
    }
    n
}

fn bench_reader(c: &mut Criterion) {
    let mut group = c.benchmark_group("psp_reader");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(30));

    // SNP-typical
    let snp_records = build_snp_records(NUM_RECORDS_SNP);
    let snp_header = writer_header(NUM_RECORDS_SNP, "snp");
    let snp_bytes = serialise(&snp_records, snp_header);
    drop(snp_records);
    eprintln!("snp_typical_3_3M: {} bytes serialised", snp_bytes.len());
    group.throughput(Throughput::Elements(NUM_RECORDS_SNP as u64));
    group.bench_function("snp_typical_3_3M", |b| {
        b.iter(|| black_box(read_all_count(black_box(&snp_bytes))));
    });

    // Phase-chain-heavy
    let phase_records = build_phase_chain_heavy_records(NUM_RECORDS_PHASE);
    let phase_header = writer_header(NUM_RECORDS_PHASE, "phase");
    let phase_bytes = serialise(&phase_records, phase_header);
    drop(phase_records);
    eprintln!(
        "phase_chain_heavy_1M: {} bytes serialised",
        phase_bytes.len()
    );
    group.throughput(Throughput::Elements(NUM_RECORDS_PHASE as u64));
    group.bench_function("phase_chain_heavy_1M", |b| {
        b.iter(|| black_box(read_all_count(black_box(&phase_bytes))));
    });

    // Multi-allele
    let multi_records = build_multi_allele_records(NUM_RECORDS_MULTI);
    let multi_header = writer_header(NUM_RECORDS_MULTI, "multi");
    let multi_bytes = serialise(&multi_records, multi_header);
    drop(multi_records);
    eprintln!("multi_allele_500k: {} bytes serialised", multi_bytes.len());
    group.throughput(Throughput::Elements(NUM_RECORDS_MULTI as u64));
    group.bench_function("multi_allele_500k", |b| {
        b.iter(|| black_box(read_all_count(black_box(&multi_bytes))));
    });

    group.finish();
}

fn bench_reader_region(c: &mut Criterion) {
    let mut group = c.benchmark_group("psp_reader_region");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(15));

    // Reuse the SNP workload: the chromosome covers 1..=NUM_RECORDS_SNP+1024.
    // A 100 k-wide window in the middle exercises one binary-search
    // hit + ~one block decode + per-record clamp.
    let records = build_snp_records(NUM_RECORDS_SNP);
    let header = writer_header(NUM_RECORDS_SNP, "snp_region");
    let bytes = serialise(&records, header);
    drop(records);

    let mid = (NUM_RECORDS_SNP / 2) as u32;
    let start = mid - 50_000;
    let end = mid + 50_000;
    group.throughput(Throughput::Elements(100_000));
    group.bench_function("region_window_chr1_mid_100k", |b| {
        b.iter(|| {
            let mut reader = PspReader::new(Cursor::new(black_box(&bytes))).expect("reader new");
            let mut n = 0u64;
            for r in reader.region_records(0, start, end) {
                let r = r.expect("record decode");
                black_box(&r);
                n += 1;
            }
            black_box(n)
        });
    });

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = bench_reader, bench_reader_region
}

criterion_main!(benches);
