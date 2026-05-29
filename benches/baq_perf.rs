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

use std::fs::{File, create_dir_all};
use std::hint::black_box;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Duration;

use criterion::{BatchSize, BenchmarkId, Criterion, criterion_group, criterion_main};
use tempfile::TempDir;

use pop_var_caller::bam::alignment_input::MappedRead;
use pop_var_caller::bam::errors::AlignmentInputError;
use pop_var_caller::baq::BaqConfig;
use pop_var_caller::fasta::{ContigEntry, ContigList, ManualEvictChromRefFetcher};
use pop_var_caller::pileup::per_sample::baq_engine::{BaqEngine, BaqOutcome};
use pop_var_caller::pileup::per_sample::baq_stream::{BaqStream, DEFAULT_BAQ_CHUNK_SIZE};
use pop_var_caller::pileup::walker::CigarOp;

/// Write a single-contig FASTA + `.fai` to `dir`. Contig is `length`
/// bases of `b'A'` (every read in the bench is Match-only over A's, so
/// the reference matches the queries exactly).
fn write_fasta(dir: &Path, contig_name: &str, length: usize) -> PathBuf {
    let fasta_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let mut fa = File::create(&fasta_path).expect("create fasta");
    writeln!(fa, ">{contig_name}").expect("fa header");
    let header_len = (contig_name.len() + 2) as u64;
    fa.write_all(&vec![b'A'; length]).expect("fa seq");
    fa.write_all(b"\n").expect("fa nl");
    let mut fai = File::create(&fai_path).expect("create fai");
    writeln!(
        fai,
        "{contig_name}\t{length}\t{header_len}\t{length}\t{}",
        length + 1
    )
    .expect("fai write");
    fasta_path
}

fn build_contig_list(name: &str, length: u64) -> ContigList {
    ContigList {
        entries: vec![ContigEntry {
            name: name.to_string(),
            length,
            md5: None,
        }],
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

fn process_all(
    engine: &mut BaqEngine,
    reads: Vec<MappedRead>,
    fetcher: &mut ManualEvictChromRefFetcher,
) -> u64 {
    let mut capped: u64 = 0;
    for r in reads {
        if let BaqOutcome::Capped(_) = engine.process(r, fetcher, true) {
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
    let contig_len = span as usize + 200;

    create_dir_all("tmp").expect("mkdir tmp");
    let dir = TempDir::new_in("tmp").expect("tempdir");
    let fasta_path = write_fasta(dir.path(), "chr0", contig_len);

    for &read_len in &[150u32, 500, 1500, 5000] {
        group.bench_with_input(BenchmarkId::from_parameter(read_len), &read_len, |b, _| {
            // The setup builds a fresh `Vec<MappedRead>` and a fresh
            // `ManualEvictChromRefFetcher` for each criterion batch.
            // The fetcher is `!Sync` and `engine.process` takes
            // `&mut` on it, so each iteration needs its own — same
            // shape `BaqStream::refill_batch` uses in production
            // (per-worker fetcher via `map_init`). Avoids the
            // per-read `cigar.clone()` / `seq.clone()` cost a
            // borrow-based bench would reintroduce on top of L3.
            b.iter_batched(
                || {
                    (
                        BaqEngine::new(BaqConfig::default()),
                        build_mapped_reads(read_len, span, coverage),
                        ManualEvictChromRefFetcher::for_contig(&fasta_path, "chr0")
                            .expect("for_contig(chr0)"),
                    )
                },
                |(mut engine, reads, mut fetcher)| {
                    black_box(process_all(&mut engine, reads, &mut fetcher))
                },
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn drain_stream(
    inputs: Vec<Result<MappedRead, AlignmentInputError>>,
    fasta_path: PathBuf,
    contigs: ContigList,
    chunk_size: usize,
) -> u64 {
    let stream = BaqStream::new(
        inputs.into_iter(),
        BaqConfig::default(),
        fasta_path,
        contigs,
        chunk_size,
        true,
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
    let contig_len = span as usize + 200;

    create_dir_all("tmp").expect("mkdir tmp");
    let dir = TempDir::new_in("tmp").expect("tempdir");
    let fasta_path = write_fasta(dir.path(), "chr0", contig_len);
    let contigs = build_contig_list("chr0", contig_len as u64);
    let reads = build_mapped_reads(read_len, span, coverage);

    for &chunk_size in &[128usize, 512, DEFAULT_BAQ_CHUNK_SIZE, 4096] {
        group.bench_with_input(
            BenchmarkId::from_parameter(chunk_size),
            &chunk_size,
            |b, &chunk_size| {
                b.iter_batched(
                    || {
                        (
                            reads
                                .iter()
                                .cloned()
                                .map(Ok::<_, AlignmentInputError>)
                                .collect::<Vec<_>>(),
                            fasta_path.clone(),
                            contigs.clone(),
                        )
                    },
                    |(inputs, fasta_path, contigs)| {
                        black_box(drain_stream(inputs, fasta_path, contigs, chunk_size))
                    },
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
