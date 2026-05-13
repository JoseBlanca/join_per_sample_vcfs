//! `.psp` writer end-to-end throughput.
//!
//! Workloads:
//!
//! - `psp_writer/snp_typical_3_3M` — single chromosome, mostly REF,
//!   0.1 % SNPs, no phase-chain markers. Dense common case.
//! - `psp_writer/phase_chain_heavy_1M` — single chromosome, sustained
//!   active phase-chain set (~6 active on average), each record's
//!   first allele carries the full active set as `chain_slots`. Drives
//!   the three `Vec<Vec<SlotId>>` list columns and the `active_slots`
//!   set bookkeeping.
//! - `psp_writer/multi_allele_500k` — single chromosome, 2–4 alleles
//!   per record, indel-shaped allele lengths from 1 to 10 bytes.
//!   Exercises the bytes column, per-allele scalars, and the per-byte
//!   ACGTN check inside `validate_record`.
//! - `psp_writer_io/bufwriter_file_1M_64KiB` — same content as
//!   `snp_typical` (1 M records) but the sink is
//!   `BufWriter::with_capacity(64 KiB, tempfile)`. Exposes the
//!   per-block `write_all` syscall cost that `io::sink()` hides.
//!
//! All four workloads build their `Vec<PileupRecord>` once outside
//! the timed region; the bench body borrows `records` and clones
//! only the small header (or recreates the sink for the file
//! variant). The `BlockAccumulator` clones the inner Vecs
//! (`new_chains`, `expired_chains`, `allele.chain_slots`) as it
//! appends, so those allocations are still measured exactly as
//! production would produce them.

use std::hint::black_box;
use std::io::{self, BufWriter, Write};
use std::time::Duration;

use criterion::{BatchSize, Criterion, Throughput, criterion_group, criterion_main};

use merge_per_sample_vcfs::per_sample_caller::pileup::{
    AlleleObservation, AlleleSupportStats, PileupRecord,
};
use merge_per_sample_vcfs::per_sample_caller::psp::header::{
    ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance,
};
use merge_per_sample_vcfs::per_sample_caller::psp::writer::PspWriter;
use std::collections::BTreeMap;

const NUM_RECORDS_SNP: usize = 3_300_000;
const NUM_RECORDS_PHASE: usize = 1_000_000;
const NUM_RECORDS_MULTI: usize = 500_000;
const NUM_RECORDS_FILE: usize = 1_000_000;

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
            subcommand: "psp-writer-perf".to_string(),
            input_crams: vec!["a.cram".to_string()],
            input_fasta: "ref.fa".to_string(),
            parameters: params,
        },
    }
}

/// SNP-typical: one allele per record (~99.9 %), occasional 2-allele.
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
        records.push(PileupRecord::new(0, pos, Vec::new(), Vec::new(), alleles));
    }
    records
}

/// Phase-chain-heavy: ~6 active phase chains on average. Every record's
/// first allele carries the current active set as `chain_slots`.
/// Markers open / close at offset cadences so `new_chains` /
/// `expired_chains` are non-empty on the boundary records.
///
/// Slot IDs cycle inside `[0, 256)` (the slot pool needs to be wider
/// than `active.len()` at all times). Active is kept sorted
/// ascending so it can be cloned straight into the validators'
/// expected order.
fn build_phase_chain_heavy_records(n: usize) -> Vec<PileupRecord> {
    let mut records = Vec::with_capacity(n);
    let bases = [b'A', b'C', b'G', b'T'];
    let mut active: Vec<u16> = Vec::new();

    fn next_unused(active: &[u16], also_avoid: &[u16]) -> u16 {
        for id in 0..=255u16 {
            if active.binary_search(&id).is_err() && !also_avoid.contains(&id) {
                return id;
            }
        }
        panic!("active set exhausted slot pool");
    }

    for i in 0..n {
        let pos = (i as u32) + 1;
        let ref_base = bases[i & 3];

        let mut new_chains: Vec<u16> = Vec::new();
        let mut expired_chains: Vec<u16> = Vec::new();

        // Bootstrap the first record by opening four slots so subsequent
        // records have something to reference / expire.
        if i == 0 {
            for &seed in &[0u16, 1, 2, 3] {
                new_chains.push(seed);
                let idx = active.partition_point(|&s| s < seed);
                active.insert(idx, seed);
            }
        } else {
            // Close one (oldest) every 13 records when there's enough
            // running to keep ≥ 3 active.
            if i % 13 == 0 && active.len() > 3 {
                let slot = active.remove(0);
                expired_chains.push(slot);
            }
            // Open one (smallest unused) every 11 records when we're
            // under the soft cap. Must avoid both the current active
            // set AND any slot we just closed in this record — the
            // writer validates `new_chains` against the pre-apply
            // active set, so reopening the just-closed slot is a
            // collision.
            if i % 11 == 0 && active.len() < 12 {
                let slot = next_unused(&active, &expired_chains);
                new_chains.push(slot);
                let idx = active.partition_point(|&s| s < slot);
                active.insert(idx, slot);
            }
        }

        // First allele carries the full active set as chain_slots.
        // active is kept sorted ascending, which is what the writer
        // requires (strictly ascending, no duplicates).
        let chain_slots = active.clone();
        let alleles = vec![AlleleObservation::new(
            vec![ref_base],
            AlleleSupportStats::new(20, -40.0, 10, 5, 2),
            chain_slots,
        )];

        records.push(PileupRecord::new(
            0,
            pos,
            new_chains,
            expired_chains,
            alleles,
        ));
    }
    records
}

/// Multi-allele: 2–4 alleles per record, allele lengths from 1 to 10
/// bytes. Cycles through {SNP-only, REF + INS, REF + DEL-shaped MNP,
/// REF + INS + alt-SNP} so the per-allele scalars + bytes columns
/// see varied shapes and the per-byte ACGTN check is exercised on
/// longer alleles.
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
                    vec![r, bases[(i + 1) & 3], bases[(i + 2) & 3], bases[(i + 3) & 3], r],
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
        records.push(PileupRecord::new(0, pos, Vec::new(), Vec::new(), alleles));
    }
    records
}

fn write_all_into<W: Write>(sink: W, records: &[PileupRecord], header: WriterHeader) -> u64 {
    let mut writer = PspWriter::new(sink, header).expect("writer new");
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

    let snp_records = build_snp_records(NUM_RECORDS_SNP);
    let snp_header = writer_header(NUM_RECORDS_SNP, "snp");
    group.throughput(Throughput::Elements(snp_records.len() as u64));
    group.bench_function("snp_typical_3_3M", |b| {
        b.iter(|| {
            black_box(write_all_into(
                io::sink(),
                black_box(&snp_records),
                black_box(snp_header.clone()),
            ))
        });
    });

    let phase_records = build_phase_chain_heavy_records(NUM_RECORDS_PHASE);
    let phase_header = writer_header(NUM_RECORDS_PHASE, "phase");
    group.throughput(Throughput::Elements(phase_records.len() as u64));
    group.bench_function("phase_chain_heavy_1M", |b| {
        b.iter(|| {
            black_box(write_all_into(
                io::sink(),
                black_box(&phase_records),
                black_box(phase_header.clone()),
            ))
        });
    });

    let multi_records = build_multi_allele_records(NUM_RECORDS_MULTI);
    let multi_header = writer_header(NUM_RECORDS_MULTI, "multi");
    group.throughput(Throughput::Elements(multi_records.len() as u64));
    group.bench_function("multi_allele_500k", |b| {
        b.iter(|| {
            black_box(write_all_into(
                io::sink(),
                black_box(&multi_records),
                black_box(multi_header.clone()),
            ))
        });
    });

    group.finish();
}

fn bench_writer_phases(c: &mut Criterion) {
    // Two sub-benches that isolate the per-record append cost from
    // the per-flush encode+compress cost. The bench-side
    // `current_block_projected_bytes` peek (on `PspWriter`) is the
    // only writer-state introspection used here — it is `#[doc(hidden)]`
    // and exists only so benches can align with the auto-flush
    // boundary deterministically.
    const TARGET_BYTES: usize = 16 * 1024 * 1024;

    let mut group = c.benchmark_group("psp_writer_phases");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(15));

    // Bench fixtures are the SNP workload — pre-built so the
    // per-iter setup only writes records, not Vec<PileupRecord>
    // construction.
    let phase_records = build_snp_records(NUM_RECORDS_SNP);
    let phase_header = writer_header(NUM_RECORDS_SNP, "phases");

    // ---- write_record steady-state ----
    // Setup: prime a writer with WARMUP records (all in the first
    // block, no auto-flush). Body: write BODY more (still under
    // TARGET, no auto-flush). Body is pure per-record append cost in
    // a warm accumulator.
    const WARMUP_NO_FLUSH: usize = 100_000;
    const BODY_NO_FLUSH: usize = 100_000;
    group.throughput(Throughput::Elements(BODY_NO_FLUSH as u64));
    group.bench_function("write_record_steady_100k", |b| {
        b.iter_batched(
            || {
                let mut w = PspWriter::new(io::sink(), phase_header.clone())
                    .expect("writer new");
                for r in &phase_records[..WARMUP_NO_FLUSH] {
                    w.write_record(r).expect("warmup write_record");
                }
                w
            },
            |mut w| {
                for r in &phase_records[WARMUP_NO_FLUSH..WARMUP_NO_FLUSH + BODY_NO_FLUSH]
                {
                    black_box(w.write_record(black_box(r)).expect("write_record"));
                }
                // do NOT finish — body must time per-record append
                // only, never the trailing flush + index + trailer.
                drop(w);
            },
            BatchSize::PerIteration,
        );
    });

    // ---- flush_block-only ----
    // Setup: prime up to projected_bytes just under TARGET so the
    // very next write triggers exactly one auto-flush. Body: write
    // that one trigger record (flush + tiny append). Reports one
    // measurement per flush — divide by element count to get
    // ms-per-flush.
    group.throughput(Throughput::Elements(1));
    group.bench_function("flush_block_one", |b| {
        b.iter_batched(
            || {
                let mut w = PspWriter::new(io::sink(), phase_header.clone())
                    .expect("writer new");
                let mut idx = 0usize;
                // Prime up to and including the record that pushes
                // projected_bytes >= TARGET. After the loop, the
                // writer's open block has projected_bytes >= TARGET,
                // so the body's first write_record will trigger
                // exactly one auto-flush (the pre-check at
                // writer.rs:113 fires before append).
                while idx < phase_records.len() {
                    w.write_record(&phase_records[idx]).expect("prime");
                    idx += 1;
                    if w.current_block_projected_bytes().unwrap_or(0) >= TARGET_BYTES {
                        break;
                    }
                }
                (w, idx)
            },
            |(mut w, idx)| {
                black_box(
                    w.write_record(black_box(&phase_records[idx]))
                        .expect("trigger write_record"),
                );
                drop(w);
            },
            BatchSize::PerIteration,
        );
    });

    group.finish();
}

fn bench_writer_io(c: &mut Criterion) {
    let mut group = c.benchmark_group("psp_writer_io");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(20));

    let file_records = build_snp_records(NUM_RECORDS_FILE);
    let file_header = writer_header(NUM_RECORDS_FILE, "file");
    group.throughput(Throughput::Elements(file_records.len() as u64));

    // Per-iteration setup creates a fresh anonymous tempfile (it is
    // unlinked on close, so the file system never sees a stale path).
    // The bench body wraps that file in `BufWriter::with_capacity(64
    // KiB, file)` — the buffer size recommended in `PspWriter::new`'s
    // rustdoc — and runs one full writer lifecycle into it.
    group.bench_function("bufwriter_file_1M_64KiB", |b| {
        b.iter_batched(
            || tempfile::tempfile().expect("tempfile"),
            |file| {
                let buf = BufWriter::with_capacity(64 * 1024, file);
                black_box(write_all_into(
                    buf,
                    black_box(&file_records),
                    black_box(file_header.clone()),
                ));
            },
            BatchSize::PerIteration,
        );
    });

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = bench_writer, bench_writer_phases, bench_writer_io
}

criterion_main!(benches);
