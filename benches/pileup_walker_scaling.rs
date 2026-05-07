//! Benchmark validating the motivation for the lazy-CIGAR refactor
//! (`ia/feature_implementation_plans/pileup_lazy_cigar.md`, finding
//! `S2` in `ia/reviews/pileup_samtools_comparison_2026-05-07.md`).
//!
//! Hypothesis: with eager CIGAR decomposition, walker work scales
//! linearly in read length `L` at fixed (window, coverage) because
//! every per-step `ReadContribution` clones its read's full event
//! vector. Doubling `L` should roughly double total walker time.
//!
//! Confirming the linear scaling validates the refactor; flat or
//! sub-linear scaling would mean the clone is not the bottleneck
//! and the plan should be revisited.

use std::hint::black_box;
use std::io;
use std::sync::Arc;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

use criterion::{BatchSize, BenchmarkId, Criterion, criterion_group, criterion_main};

use merge_per_sample_vcfs::per_sample_caller::pileup::{
    CigarOp, PileupRecord, PreparedRead, RefBaseFetcher, run,
};

/// Constant-base reference. The walker only needs `RefBaseFetcher`,
/// so we avoid the FASTA toolchain entirely.
struct ConstFasta {
    len: usize,
}

impl RefBaseFetcher for ConstFasta {
    fn fetch(&self, _chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        let lo = (start_1based - 1) as usize;
        let hi = lo + length as usize;
        if hi > self.len {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "off end"));
        }
        Ok(vec![b'A'; length as usize])
    }
}

/// Build a pure-Match input stream of `num_reads` reads of length
/// `read_len`, evenly spaced across `[1, span]` so that average
/// per-position coverage is approximately `coverage`.
///
/// Pure-Match reads are the worst case for the eager-clone hot
/// path: every walker step the read is alive sees a Match event
/// at `walker_pos`, so the per-step `full_window_events.clone()`
/// fires every position the read covers.
fn build_reads(read_len: u32, span: u32, coverage: u32) -> Vec<PreparedRead> {
    assert!(read_len >= 1 && span >= read_len && coverage >= 1);
    // Number of reads = span * coverage / read_len, rounded down.
    let num_reads = ((span as u64) * (coverage as u64) / (read_len as u64)) as u32;
    // Distribute starts evenly over [1, span - read_len + 1].
    let last_start = span - read_len + 1;
    let mut reads = Vec::with_capacity(num_reads as usize);
    for i in 0..num_reads {
        // Even spacing; multiple reads can share a start position
        // when num_reads exceeds the number of distinct starts.
        let start = 1 + ((i as u64) * (last_start as u64 - 1) / num_reads.max(1) as u64) as u32;
        reads.push(PreparedRead {
            chrom_id: 0,
            alignment_start: start,
            alignment_end: start + read_len - 1,
            cigar: vec![CigarOp::Match(read_len)],
            seq: vec![b'A'; read_len as usize],
            bq_baq: vec![30u8; read_len as usize],
            mq_log_err: -3.0,
            is_reverse_strand: false,
            qname: Arc::from(format!("r{i}").as_str()),
            is_first_mate: true,
            has_mate: false,
        });
    }
    // Walker requires non-decreasing alignment_start.
    reads.sort_by_key(|r| r.alignment_start);
    reads
}

fn run_walker(reads: Vec<PreparedRead>, fasta: &ConstFasta) -> u64 {
    let (tx, rx) = mpsc::sync_channel::<PileupRecord>(64);
    let collector = thread::spawn(move || {
        let mut count = 0u64;
        while rx.recv().is_ok() {
            count += 1;
        }
        count
    });
    run(reads, fasta, &tx).expect("walker run failed");
    drop(tx);
    collector.join().expect("collector panicked")
}

fn bench_walker_read_length(c: &mut Criterion) {
    let mut group = c.benchmark_group("pileup_walker_read_length");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(3));

    // Hold (span, coverage) constant so total work is dominated by
    // per-active-read clone cost, not by the size of the input
    // stream. With span = 50 kb and coverage = 30, total reads
    // shrink as `read_len` grows but each read clones a larger
    // event vector at each walker step.
    let span: u32 = 50_000;
    let coverage: u32 = 30;
    let fasta = ConstFasta {
        len: span as usize + 100,
    };

    // Read lengths swept inside the current `MAX_RECORD_SPAN = 5000`
    // ceiling. This is the same ceiling the lazy-CIGAR plan would
    // raise to enable long-read support; the benchmark is the
    // upper-bound test under today's cap.
    for &read_len in &[150u32, 500, 1500, 5000] {
        group.bench_with_input(
            BenchmarkId::from_parameter(read_len),
            &read_len,
            |b, &read_len| {
                b.iter_batched(
                    || build_reads(read_len, span, coverage),
                    |reads| black_box(run_walker(reads, &fasta)),
                    BatchSize::LargeInput,
                );
            },
        );
    }
    group.finish();
}

fn config() -> Criterion {
    Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(3))
}

criterion_group! {
    name = benches;
    config = config();
    targets = bench_walker_read_length
}

criterion_main!(benches);
