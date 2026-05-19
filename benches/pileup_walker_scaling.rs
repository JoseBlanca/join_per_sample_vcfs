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
use std::time::Duration;

use criterion::{BatchSize, BenchmarkId, Criterion, criterion_group, criterion_main};

use pop_var_caller::per_sample_pileup::pileup::{
    CigarOp, MateRole, PreparedRead, RefSeqFetcher, WalkerConfig, run,
};

/// Constant-base reference. The walker only needs `RefSeqFetcher`,
/// so we avoid the FASTA toolchain entirely.
struct ConstFasta {
    len: usize,
}

impl RefSeqFetcher for ConstFasta {
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
            mate_role: MateRole::Solo,
            adaptor_boundary: None,
        });
    }
    // Walker requires non-decreasing alignment_start.
    reads.sort_by_key(|r| r.alignment_start);
    reads
}

/// Build a CIGAR with periodic short insertions. Each cycle is
/// `M(cycle_len)` followed by `I(ins_len)`; the cigar ends with
/// an `M` op so the cursor's first/last-op-indel drop rule
/// doesn't silently strip a trailing insertion. Total reference
/// span equals `ref_span` exactly; total read seq is longer than
/// `ref_span` by the inserted bases.
///
/// **Why insertions, not deletions.** An earlier draft used
/// `M(50) D(2)` cycles, which at high coverage caused overlapping
/// deletion footprints from different reads to chain-widen the
/// same open record past `MAX_RECORD_SPAN = 5000` and trip the
/// walker's safety cap. Insertions have footprint span 1 (just
/// the anchor base), so they never widen records — the fixture
/// scales op count with `ref_span` without putting the walker
/// into its pathological-coverage regime.
///
/// The shape roughly mimics what long-read aligners produce on
/// noisy data: tens to hundreds of small CIGAR ops per read,
/// scaling with `ref_span`. This is the regime where
/// `CigarCursor`'s op-loop cost matters and where the early-break
/// (and a future binary-search-on-offsets) would show up.
fn build_multi_op_cigar(ref_span: u32, cycle_len: u32, ins_len: u32) -> Vec<CigarOp> {
    assert!(cycle_len > 0 && ins_len > 0);
    let mut ops = Vec::new();
    let mut consumed = 0u32;
    // Leave at least one full M op after the last I so the cigar
    // ends in M, not I.
    while consumed + cycle_len < ref_span {
        ops.push(CigarOp::Match(cycle_len));
        ops.push(CigarOp::Insertion(ins_len));
        // Insertions don't consume reference bases.
        consumed += cycle_len;
    }
    let remainder = ref_span - consumed;
    if remainder > 0 {
        ops.push(CigarOp::Match(remainder));
    }
    ops
}

/// Sum the read-consuming op lengths of a CIGAR — the read seq
/// length needed to back this CIGAR.
fn read_seq_len_for_cigar(cigar: &[CigarOp]) -> u32 {
    cigar
        .iter()
        .map(|op| match *op {
            CigarOp::Match(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Insertion(n)
            | CigarOp::SoftClip(n) => n,
            _ => 0,
        })
        .sum()
}

/// Multi-op variant of `build_reads`. Same span/coverage/start
/// distribution; each read carries the periodic-insertion CIGAR
/// from `build_multi_op_cigar(read_len, 50, 2)`. Approximate op
/// counts:
///
/// |  L | total CIGAR ops |
/// |---:|---:|
/// |  150 |   5 |
/// |  500 |  19 |
/// | 1500 |  59 |
/// | 5000 | 199 |
///
/// All reads in one fixture share the canonical CIGAR (only their
/// alignment_start varies), so we build it once and clone.
fn build_multi_op_reads(read_len: u32, span: u32, coverage: u32) -> Vec<PreparedRead> {
    assert!(read_len >= 60 && span >= read_len && coverage >= 1);
    let num_reads = ((span as u64) * (coverage as u64) / (read_len as u64)) as u32;
    let last_start = span - read_len + 1;
    let cigar = build_multi_op_cigar(read_len, 50, 2);
    let read_consumed = read_seq_len_for_cigar(&cigar);
    let seq = vec![b'A'; read_consumed as usize];
    let bq_baq = vec![30u8; read_consumed as usize];

    let mut reads = Vec::with_capacity(num_reads as usize);
    for i in 0..num_reads {
        let start = 1 + ((i as u64) * (last_start as u64 - 1) / num_reads.max(1) as u64) as u32;
        reads.push(PreparedRead {
            chrom_id: 0,
            alignment_start: start,
            alignment_end: start + read_len - 1,
            cigar: cigar.clone(),
            seq: seq.clone(),
            bq_baq: bq_baq.clone(),
            mq_log_err: -3.0,
            is_reverse_strand: false,
            qname: Arc::from(format!("r{i}").as_str()),
            mate_role: MateRole::Solo,
            adaptor_boundary: None,
        });
    }
    reads.sort_by_key(|r| r.alignment_start);
    reads
}

/// Per-iteration timed body: iterate the pull-shaped walker to
/// exhaustion, counting yielded records. After the rewrite from
/// `SyncSender` push to `Iterator` pull, there's no channel to
/// allocate and no collector thread to spawn, so the previous
/// `iter_batched` setup-closure trick (round-2 L2 in
/// `ia/reviews/perf_pileup_2026-05-12.md`) is no longer needed —
/// the timed body now measures only walker work plus the per-record
/// `next()` dispatch.
fn drive_walker(reads: Vec<PreparedRead>, ref_fetcher: &ConstFasta) -> u64 {
    let mut count: u64 = 0;
    for item in run(reads, ref_fetcher, &WalkerConfig::default()) {
        item.expect("walker yielded error");
        count += 1;
    }
    count
}

fn bench_walker_read_length(c: &mut Criterion) {
    let mut group = c.benchmark_group("pileup_walker_read_length");
    group.sample_size(30);
    group.measurement_time(Duration::from_secs(10));

    // Hold (span, coverage) constant so total work is dominated by
    // per-active-read clone cost, not by the size of the input
    // stream. With span = 50 kb and coverage = 30, total reads
    // shrink as `read_len` grows but each read clones a larger
    // event vector at each walker step.
    let span: u32 = 50_000;
    let coverage: u32 = 30;
    let ref_fetcher = ConstFasta {
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
                    |reads| black_box(drive_walker(reads, &ref_fetcher)),
                    BatchSize::PerIteration,
                );
            },
        );
    }
    group.finish();
}

/// Multi-op variant of `bench_walker_read_length`. Same
/// (span, coverage, L sweep), but each read has a many-op CIGAR
/// (periodic small deletions). Designed as the forward-looking
/// guide for changes that affect the cursor's op-loop cost — the
/// early-break already shipped, and the binary-search-on-offsets
/// idea queued in the plan's "out-of-scope follow-ups". The
/// single-op fixture (`bench_walker_read_length`) cannot
/// distinguish these because its loops run exactly once per query.
fn bench_walker_multi_op(c: &mut Criterion) {
    let mut group = c.benchmark_group("pileup_walker_multi_op");
    group.sample_size(30);
    group.measurement_time(Duration::from_secs(10));

    let span: u32 = 50_000;
    let coverage: u32 = 30;
    let ref_fetcher = ConstFasta {
        len: span as usize + 100,
    };

    for &read_len in &[150u32, 500, 1500, 5000] {
        group.bench_with_input(
            BenchmarkId::from_parameter(read_len),
            &read_len,
            |b, &read_len| {
                b.iter_batched(
                    || build_multi_op_reads(read_len, span, coverage),
                    |reads| black_box(drive_walker(reads, &ref_fetcher)),
                    BatchSize::PerIteration,
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
    targets = bench_walker_read_length, bench_walker_multi_op
}

criterion_main!(benches);
