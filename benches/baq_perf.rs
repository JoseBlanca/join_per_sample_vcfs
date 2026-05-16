//! BAQ stage benchmarks. Two scenarios:
//!
//! 1. `baq_engine_read_length` — single-thread `BaqEngine::process`
//!    over a fixed (span, coverage) workload, sweeping read length.
//!    Measures the per-read cost of the HMM + CIGAR walks + cap pass
//!    at steady state — the engine is reused across all reads in one
//!    iteration so scratch buffers stay warm.
//!
//! 2. `baq_stream_chunk_size` — `BaqStream` end-to-end at L = 150
//!    (typical Illumina), sweeping `chunk_size`. Measures the
//!    rayon-parallel pipeline including chunk dispatch and the serial
//!    drain. Picks the chunk size that wins under default rayon.
//!
//! Fixtures use Match-only CIGARs. A multi-op variant is a follow-up
//! if profiling later points at the CIGAR-walk cost as a bottleneck.

use std::hint::black_box;
use std::io;
use std::time::Duration;

use criterion::{BatchSize, BenchmarkId, Criterion, criterion_group, criterion_main};

use merge_per_sample_vcfs::per_sample_pileup::baq::{
    BaqConfig, BaqEngine, BaqOutcome, BaqStream, DEFAULT_BAQ_CHUNK_SIZE,
};
use merge_per_sample_vcfs::per_sample_pileup::cram_input::{CigarOp, MappedRead};
use merge_per_sample_vcfs::per_sample_pileup::errors::CramInputError;
use merge_per_sample_vcfs::per_sample_pileup::pileup::RefSeqFetcher;

/// Constant-base reference: every fetch returns `vec![b'A'; length]`.
/// Same trick the pileup walker bench uses to avoid the FASTA
/// toolchain.
struct ConstRefFetcher {
    len: usize,
}

impl RefSeqFetcher for ConstRefFetcher {
    fn fetch(&self, _chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        let lo = (start_1based - 1) as usize;
        let hi = lo + length as usize;
        if hi > self.len {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "off end"));
        }
        Ok(vec![b'A'; length as usize])
    }
}

/// Build `coverage * span / read_len` Match-only `MappedRead`s at
/// evenly-spaced 1-based positions across `[1, span - read_len + 1]`.
/// All bases A, all Q30, solo / forward strand.
fn build_mapped_reads(read_len: u32, span: u32, coverage: u32) -> Vec<MappedRead> {
    assert!(read_len >= 1 && span >= read_len && coverage >= 1);
    let num_reads = ((span as u64) * (coverage as u64) / (read_len as u64)).max(1) as u32;
    let last_start = span - read_len + 1;
    let mut reads = Vec::with_capacity(num_reads as usize);
    for i in 0..num_reads {
        let start = 1 + ((i as u64) * (last_start as u64 - 1) / num_reads.max(1) as u64) as u32;
        reads.push(MappedRead {
            qname: format!("r{i}").into_bytes(),
            flag: 0,
            ref_id: 0,
            pos: start as u64,
            mapq: 60,
            cigar: vec![CigarOp::Match(read_len)],
            seq: vec![b'A'; read_len as usize],
            qual: vec![30u8; read_len as usize],
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            source_file_index: 0,
        });
    }
    reads.sort_by_key(|r| r.pos);
    reads
}

fn process_all(engine: &mut BaqEngine, reads: Vec<MappedRead>, fetcher: &ConstRefFetcher) -> u64 {
    let mut capped: u64 = 0;
    for r in reads {
        if let BaqOutcome::Capped(_) = engine.process(r, fetcher) {
            capped += 1;
        }
    }
    capped
}

fn bench_engine_read_length(c: &mut Criterion) {
    let mut group = c.benchmark_group("baq_engine_read_length");
    group.sample_size(20);
    group.measurement_time(Duration::from_secs(10));

    let span: u32 = 50_000;
    let coverage: u32 = 30;
    let fetcher = ConstRefFetcher {
        len: span as usize + 200,
    };

    for &read_len in &[150u32, 500, 1500, 5000] {
        group.bench_with_input(BenchmarkId::from_parameter(read_len), &read_len, |b, _| {
            // The setup builds a fresh `Vec<MappedRead>` for each
            // criterion batch and the bench body moves it into
            // `process_all`, matching how `BaqStream::refill_batch`
            // consumes its chunk in production (par_drain). Avoids
            // the per-read `cigar.clone()` / `seq.clone()` cost that
            // a borrow-based bench would re-introduce on top of L3.
            b.iter_batched(
                || {
                    (
                        BaqEngine::new(BaqConfig::default()),
                        build_mapped_reads(read_len, span, coverage),
                    )
                },
                |(mut engine, reads)| black_box(process_all(&mut engine, reads, &fetcher)),
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn drain_stream<F: RefSeqFetcher + Sync>(
    inputs: Vec<Result<MappedRead, CramInputError>>,
    fetcher: &F,
    chunk_size: usize,
) -> u64 {
    let stream = BaqStream::new(
        inputs.into_iter(),
        BaqConfig::default(),
        fetcher,
        chunk_size,
    );
    let mut capped: u64 = 0;
    for outcome in stream {
        if outcome.is_ok() {
            capped += 1;
        }
    }
    capped
}

fn bench_stream_chunk_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("baq_stream_chunk_size");
    group.sample_size(20);
    group.measurement_time(Duration::from_secs(10));

    let span: u32 = 50_000;
    let coverage: u32 = 30;
    let read_len: u32 = 150;
    let fetcher = ConstRefFetcher {
        len: span as usize + 200,
    };
    let reads = build_mapped_reads(read_len, span, coverage);

    for &chunk_size in &[128usize, 512, DEFAULT_BAQ_CHUNK_SIZE, 4096] {
        group.bench_with_input(
            BenchmarkId::from_parameter(chunk_size),
            &chunk_size,
            |b, &chunk_size| {
                b.iter_batched(
                    || {
                        reads
                            .iter()
                            .cloned()
                            .map(Ok::<_, CramInputError>)
                            .collect::<Vec<_>>()
                    },
                    |inputs| black_box(drain_stream(inputs, &fetcher, chunk_size)),
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
    targets = bench_engine_read_length, bench_stream_chunk_size
}

criterion_main!(benches);
