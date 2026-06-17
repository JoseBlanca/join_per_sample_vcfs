//! `ssr-pileup` realignment hot-path throughput.
//!
//! Measures the CPU-bound core of Stage 1 — per spanning read: triage →
//! `build_rungs` → `score_candidates` (the pair-HMM forward over `2·window + 1`
//! candidate rungs). This is `process_locus` minus the indexed fetch (the fetch
//! is I/O-bound and benched separately against real alignment files).
//!
//! The workload is one synthetic locus with `depth` clean spanning reads,
//! exercised through the `bench_harness` seam. Each group reports throughput in
//! reads/s so configurations with different `depth` are comparable.
//!
//! Workloads:
//!
//! - `ssr_realign/period` — homopolymer / di / tri / tetra motif, fixed
//!   units=15, flank=50, depth=30, window=6 (the shipped default). Dinucleotide
//!   is the common case.
//! - `ssr_realign/window` — window 3 vs 6 vs 15 at the dinucleotide default;
//!   isolates the per-read candidate count (`2·window + 1` forward passes).
//! - `ssr_realign/units` — short (8) vs long (30) tract; isolates the per-rung
//!   haplotype length (the DP's inner dimension).
//! - `ssr_realign/softclip` — a long allele recovered from a trailing soft-clip
//!   (longer haplotypes + the clip extraction path) vs the clean baseline.

// Opt-in mimalloc global allocator (cargo bench --features alloc-mimalloc ...).
#[cfg(feature = "alloc-mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::hint::black_box;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

use pop_var_caller::ssr::pileup::bench_harness::{analyze_workload, build_synthetic_workload};

const DEPTH: usize = 30;
const FLANK_BP: usize = 50;

fn bench_period(c: &mut Criterion) {
    let mut g = c.benchmark_group("ssr_realign/period");
    g.throughput(Throughput::Elements(DEPTH as u64));
    for (label, period) in [("homopolymer", 1), ("di", 2), ("tri", 3), ("tetra", 4)] {
        let w = build_synthetic_workload(period, 15, FLANK_BP, DEPTH, 6, 0);
        g.bench_with_input(BenchmarkId::from_parameter(label), &w, |b, w| {
            b.iter(|| black_box(analyze_workload(black_box(w))));
        });
    }
    g.finish();
}

fn bench_window(c: &mut Criterion) {
    let mut g = c.benchmark_group("ssr_realign/window");
    g.throughput(Throughput::Elements(DEPTH as u64));
    for window in [3u16, 6, 15] {
        let w = build_synthetic_workload(2, 15, FLANK_BP, DEPTH, window, 0);
        g.bench_with_input(BenchmarkId::from_parameter(window), &w, |b, w| {
            b.iter(|| black_box(analyze_workload(black_box(w))));
        });
    }
    g.finish();
}

fn bench_units(c: &mut Criterion) {
    let mut g = c.benchmark_group("ssr_realign/units");
    g.throughput(Throughput::Elements(DEPTH as u64));
    for units in [8u16, 30] {
        let w = build_synthetic_workload(2, units, FLANK_BP, DEPTH, 6, 0);
        g.bench_with_input(BenchmarkId::from_parameter(units), &w, |b, w| {
            b.iter(|| black_box(analyze_workload(black_box(w))));
        });
    }
    g.finish();
}

fn bench_softclip(c: &mut Criterion) {
    let mut g = c.benchmark_group("ssr_realign/softclip");
    g.throughput(Throughput::Elements(DEPTH as u64));
    for (label, clip_units) in [("clean", 0u16), ("clip10", 10)] {
        let w = build_synthetic_workload(2, 15, FLANK_BP, DEPTH, 6, clip_units);
        g.bench_with_input(BenchmarkId::from_parameter(label), &w, |b, w| {
            b.iter(|| black_box(analyze_workload(black_box(w))));
        });
    }
    g.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(8));
    targets = bench_period, bench_window, bench_units, bench_softclip
}
criterion_main!(benches);
