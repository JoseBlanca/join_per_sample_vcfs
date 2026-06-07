//! Cohort `var-calling` end-to-end throughput (`.psp` cohort → VCF).
//!
//! The committed perf surface for the re-architected record-streaming
//! pipeline (producer → caller workers → writer). Replaces the deleted
//! `cohort_e2e_perf` bench (which only exercised the removed streaming
//! `drive_cohort_pipeline`). Built to make the **H1 thread-oversubscription**
//! finding rankable: the wall-time gap that grows with thread count at the
//! production `target-variants-per-chunk = 256`.
//!
//! ## What it measures
//!
//! `var_calling/N{n}/T{t}` — drive
//! [`pop_var_caller::var_calling::pipeline::run_var_calling`] over a fully
//! synthetic N-replicate cohort, sweeping `N ∈ {1, 8, 50}` samples ×
//! `t ∈ {1, 2, 8}` worker threads at the production chunk size. The fixture
//! (one synthetic chromosome of `.psp` records + a matching synthetic FASTA)
//! is built **once per N outside the timed region**; only the
//! `run_var_calling` call is timed. Every run asserts on the returned
//! [`WriterStats::records_written`](pop_var_caller::var_calling::vcf_writer::WriterStats)
//! so a fixture or chunking regression trips the bench instead of silently
//! changing the workload.
//!
//! ## Why the pipeline entry point + a per-combo `install`d pool
//!
//! The CLI-level entry calls `configure_rayon_pool`, which is **once per
//! process** — a criterion run that sweeps the thread count in one process
//! could only set rayon's global pool once, and the producer's `par_iter_mut`
//! decode/compaction would then run on a fixed all-cores pool regardless of
//! `t` (measuring only caller scaling, the *opposite* of the production gap).
//!
//! Instead we drive the **pipeline** entry inside a per-combo
//! `rayon::ThreadPool` of `t` threads via [`rayon::ThreadPool::install`]: the
//! producer runs on the calling (main) thread inside that scope, so its
//! `par_iter_mut` uses the `t`-sized pool — exactly what production's
//! `configure_rayon_pool(t)` does to the global pool. The `t` crossbeam
//! **caller** threads are spawned separately, so each timed run has `t` rayon
//! producer workers **plus** `t` callers live concurrently — the `~2·t`
//! oversubscription on the box's cores that **H1** is about. The pool is built
//! once per combo (outside `b.iter`) so thread spawn cost is not timed.
//!
//! This composes with the upcoming Variant-1 fix: once the production code
//! caps the producer's pool via its *own* inner `install`, that inner pool
//! wins over this outer one, so the same bench then measures the fix (callers
//! still `t`, producer capped) — the high-`t` wall should fall back toward the
//! single-pool curve. The authoritative cross-process "vs `main`" sweep on
//! real data (`examples/profile_cohort_e2e` + hyperfine, perf review §3.3)
//! stays the headline experiment; this is the committed in-process guard.
//!
//! Run (container, production arch):
//! ```text
//! ./scripts/dev.sh cargo bench --bench cohort_var_calling_perf
//! # single combo:
//! ./scripts/dev.sh cargo bench --bench cohort_var_calling_perf -- N8/T2
//! ```

// Opt-in mimalloc global allocator (cargo bench --features alloc-mimalloc ...).
#[cfg(feature = "alloc-mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::collections::BTreeMap;
use std::fs::File;
use std::hint::black_box;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use tempfile::TempDir;

use pop_var_caller::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};
use pop_var_caller::pop_var_caller::cli::shared_args::CohortPipelineArgs;
use pop_var_caller::pop_var_caller::var_calling::VarCallingArgs;
use pop_var_caller::psp::header::{
    ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance,
};
use pop_var_caller::psp::writer::PspWriter;
use pop_var_caller::var_calling::dust_filter::{DEFAULT_DUST_THRESHOLD, DEFAULT_DUST_WINDOW};
use pop_var_caller::var_calling::per_group_merger::{
    DEFAULT_MAX_ALLELES_LH_CALC, DEFAULT_MAX_ALLELES_PER_RECORD, DEFAULT_PLOIDY,
};
use pop_var_caller::var_calling::pipeline::run_var_calling;
use pop_var_caller::var_calling::posterior_engine::{
    DEFAULT_COMPOUND_ALT_PSEUDOCOUNT, DEFAULT_CONVERGENCE_THRESHOLD,
    DEFAULT_INBREEDING_COEFFICIENT, DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MAX_GQ_PHRED,
    DEFAULT_MAX_ITERATIONS, DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT,
};
use pop_var_caller::var_calling::variant_grouping::DEFAULT_MAX_VARIANT_GROUP_SPAN;
use pop_var_caller::var_calling::{
    DEFAULT_MIN_ALT_OBS_PER_SAMPLE, DEFAULT_MIN_MAPQ_DIFF_T, DEFAULT_MIN_QUAL_PHRED,
};

/// Single synthetic contig.
const CONTIG_NAME: &str = "chr1";
/// Positions carrying a `.psp` record (and the synthetic FASTA length is this
/// plus a margin). Kept modest so the heaviest combo (N=50) is a few hundred
/// ms — fast enough to sweep all combos, large enough that the
/// producer→caller overlap (not per-call fixed overhead) dominates.
const RECORD_POSITIONS: u32 = 60_000;
/// One in every `VARIANT_PERIOD` positions carries a second (ALT) allele; the
/// rest are REF-only. ~3.3 % variant density drives the grouper + merger + EM.
const VARIANT_PERIOD: u32 = 30;
/// Production per-chunk variable-position target (the `0` sentinel maps to this
/// in the pipeline; set explicitly so the bench is pinned to the production
/// shape regardless of the sentinel default).
const TARGET_VARIANTS_PER_CHUNK: u32 = 256;
/// Sample-count axis (cohort size).
const SAMPLE_COUNTS: &[usize] = &[1, 8, 50];
/// Worker-thread axis (crossbeam caller pool; see module doc on the rayon
/// caveat).
const THREAD_COUNTS: &[usize] = &[1, 2, 8];

/// Deterministic high-entropy base for FASTA position `i`. Uses a full
/// avalanche integer hash (lowbias32) so the sequence is non-periodic — a
/// naive `i*K`-style mix leaves the low bits periodic, which DUST flags as
/// low-complexity and masks wholesale, dropping every variant. With good
/// entropy the DUST pre-pass produces a near-empty mask, matching the
/// production hot path (most positions unmasked).
fn fasta_base(i: u32) -> u8 {
    const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut x = i.wrapping_add(0x9E37_79B9);
    x = (x ^ (x >> 16)).wrapping_mul(0x7feb_352d);
    x = (x ^ (x >> 15)).wrapping_mul(0x846c_a68b);
    x ^= x >> 16;
    B[(x & 3) as usize]
}

/// Write a single-contig FASTA + `.fai` (one-line contig, mirroring the
/// integration-test fixture layout). Length = `RECORD_POSITIONS + 1024`.
fn build_fasta(dir: &Path) -> (PathBuf, u32) {
    let contig_len = RECORD_POSITIONS + 1024;
    let fasta_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let mut fa = BufWriter::new(File::create(&fasta_path).expect("create fasta"));
    writeln!(fa, ">{CONTIG_NAME}").unwrap();
    let header_len = CONTIG_NAME.len() + 2; // '>' + name + '\n'
    let seq: Vec<u8> = (0..contig_len).map(fasta_base).collect();
    fa.write_all(&seq).unwrap();
    fa.write_all(b"\n").unwrap();
    fa.into_inner().unwrap().sync_all().unwrap();

    let mut fai = File::create(&fai_path).expect("create fai");
    writeln!(
        fai,
        "{}\t{}\t{}\t{}\t{}",
        CONTIG_NAME,
        contig_len,
        header_len,
        contig_len,
        contig_len + 1
    )
    .unwrap();
    fai.sync_all().unwrap();
    (fasta_path, contig_len)
}

/// Realistic high-confidence support for `num_obs` observations: per-obs mean
/// quality ≈ Q40 (`q_sum = num_obs * -4.0`, comfortably above the QUAL-30
/// emit threshold), balanced strand (no strand bias), MAPQ 60 with zero
/// variance (so the MAPQ-diff Welch's-t filter stays quiet).
fn support(num_obs: u32) -> AlleleSupportStats {
    let n = num_obs as f64;
    AlleleSupportStats::new(
        num_obs,
        n * -4.0,
        num_obs / 2,                // fwd: balanced strand
        num_obs / 2,                // placed_left: balanced placement
        0,                          // placed_start
        num_obs * 60,               // mapq_sum: MAPQ 60 each
        u64::from(num_obs) * 3_600, // mapq_sum_sq: MAPQ² each → zero variance
    )
}

/// Synthetic per-sample records: every position `1..=RECORD_POSITIONS` is
/// REF-only except every `VARIANT_PERIOD`-th, which carries an ALT. REF bases
/// match the FASTA so `ref_span` reads are consistent. Identical across samples
/// (replicas), so every variant position is cohort-wide variable → kept by the
/// producer's variant filter and called.
fn build_records() -> Vec<PileupRecord> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut records = Vec::with_capacity(RECORD_POSITIONS as usize);
    for pos in 1..=RECORD_POSITIONS {
        let ref_base = fasta_base(pos - 1);
        let is_variant = pos % VARIANT_PERIOD == 0;
        let mut alleles = Vec::with_capacity(if is_variant { 2 } else { 1 });
        // REF depth ~28-31; ALT depth ~12-16 → a clear heterozygous call.
        alleles.push(AlleleObservation::new(
            vec![ref_base],
            support(28 + (pos % 4)),
            Vec::new(),
        ));
        if is_variant {
            // ALT base distinct from REF.
            let alt = bases.iter().copied().find(|&b| b != ref_base).unwrap();
            alleles.push(AlleleObservation::new(
                vec![alt],
                support(12 + (pos % 5)),
                Vec::new(),
            ));
        }
        records.push(PileupRecord::new(0, pos, alleles));
    }
    records
}

fn writer_header(sample: &str, contig_len: u32) -> WriterHeader {
    let mut params = BTreeMap::new();
    params.insert("min-mapq".to_string(), ParameterValue::Integer(30));
    WriterHeader {
        format_version: (1, 0),
        sample: sample.to_string(),
        reference: "ref.fa".to_string(),
        created: "2026-06-07T00:00:00Z".parse().unwrap(),
        chromosomes: vec![ChromosomeEntry {
            name: CONTIG_NAME.to_string(),
            length: contig_len,
            // The pipeline entry point does not cross-check this MD5 (the CLI
            // wrapper does); any placeholder is fine here.
            md5: "0".repeat(32),
        }],
        writer: WriterProvenance {
            tool: "bench".to_string(),
            version: "0.0.0".to_string(),
            subcommand: "cohort-var-calling-perf".to_string(),
            input_crams: vec!["a.cram".to_string()],
            input_fasta: "ref.fa".to_string(),
            command_line: String::new(),
            parameters: params,
        },
    }
}

fn write_psp(path: &Path, header: WriterHeader, records: &[PileupRecord]) {
    let sink = BufWriter::with_capacity(64 * 1024, File::create(path).expect("create psp"));
    let mut writer = PspWriter::new(sink, header).expect("psp writer");
    for record in records {
        writer.write_record(record).expect("write record");
    }
    let buf = writer.finish().expect("finish psp");
    buf.into_inner()
        .expect("into_inner")
        .sync_all()
        .expect("sync");
}

/// A built cohort fixture. `TempDir` is held so the on-disk `.psp`/FASTA
/// outlive every timed iteration for this N.
struct Fixture {
    _dir: TempDir,
    reference: PathBuf,
    psp_files: Vec<PathBuf>,
    output: PathBuf,
}

fn build_fixture(n_samples: usize) -> Fixture {
    let dir = tempfile::tempdir().expect("tempdir");
    let (reference, contig_len) = build_fasta(dir.path());
    let records = build_records();
    let mut psp_files = Vec::with_capacity(n_samples);
    for s in 0..n_samples {
        let sample = format!("S{s:04}");
        let path = dir.path().join(format!("{sample}.psp"));
        write_psp(&path, writer_header(&sample, contig_len), &records);
        psp_files.push(path);
    }
    let output = dir.path().join("out.vcf");
    Fixture {
        _dir: dir,
        reference,
        psp_files,
        output,
    }
}

/// Build the pipeline args for one timed run. `VarCallingArgs` is constructed
/// fresh per call (cheap; the heavy fixture lives on disk) with the swept
/// `threads`.
fn make_args(fx: &Fixture, threads: usize) -> VarCallingArgs {
    VarCallingArgs {
        reference: fx.reference.clone(),
        output: fx.output.clone(),
        regions: None,
        threads: Some(threads),
        contamination_estimates: None,
        no_complexity_filter: false,
        target_variants_per_chunk: TARGET_VARIANTS_PER_CHUNK,
        low_memory: false,
        psp_files: fx.psp_files.clone(),
        cohort: CohortPipelineArgs {
            ploidy: DEFAULT_PLOIDY,
            complexity_window: DEFAULT_DUST_WINDOW,
            complexity_threshold: DEFAULT_DUST_THRESHOLD,
            var_group_max_span: DEFAULT_MAX_VARIANT_GROUP_SPAN,
            max_alleles_per_var: DEFAULT_MAX_ALLELES_PER_RECORD,
            max_alleles_lh_calc: DEFAULT_MAX_ALLELES_LH_CALC,
            inbreeding_coefficient: DEFAULT_INBREEDING_COEFFICIENT,
            em_convergence_threshold: DEFAULT_CONVERGENCE_THRESHOLD,
            em_max_iterations: DEFAULT_MAX_ITERATIONS,
            ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT,
            snp_alt_pseudocount: DEFAULT_SNP_ALT_PSEUDOCOUNT,
            indel_alt_pseudocount: DEFAULT_INDEL_ALT_PSEUDOCOUNT,
            compound_alt_pseudocount: DEFAULT_COMPOUND_ALT_PSEUDOCOUNT,
            max_gq_phred: DEFAULT_MAX_GQ_PHRED,
            min_qual_phred: DEFAULT_MIN_QUAL_PHRED,
            min_alt_obs_per_sample: DEFAULT_MIN_ALT_OBS_PER_SAMPLE,
            no_mapq_diff_filter: false,
            min_mapq_diff_t: DEFAULT_MIN_MAPQ_DIFF_T,
            emit_gp: false,
        },
    }
}

fn bench_cohort(c: &mut Criterion) {
    let mut group = c.benchmark_group("var_calling");
    // Each run is hundreds of ms at N=50; keep the sweep bounded.
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(12));
    group.warm_up_time(Duration::from_secs(2));

    for &n in SAMPLE_COUNTS {
        let fx = build_fixture(n);
        // One untimed run pins the expected emit count so the in-loop assert
        // catches a fixture/chunking regression (output is thread-independent
        // by design, so any thread count yields the same count).
        let expected = run_var_calling(&make_args(&fx, 1), None)
            .expect("warm pipeline run")
            .records_written;
        assert!(
            expected > 0,
            "fixture produced no calls — bad synthetic data"
        );

        for &t in THREAD_COUNTS {
            // Per-combo producer pool of `t` threads (built once, outside the
            // timed region). `install` makes the producer's `par_iter_mut` use
            // it — reproducing production's `configure_rayon_pool(t)` so the
            // `t` axis drives both the producer pool *and* the `t` callers.
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(t)
                .build()
                .expect("build producer pool");
            group.bench_function(format!("N{n}/T{t}"), |b| {
                b.iter(|| {
                    let args = make_args(&fx, t);
                    let stats = pool
                        .install(|| run_var_calling(black_box(&args), None))
                        .expect("pipeline run");
                    assert_eq!(
                        stats.records_written, expected,
                        "records_written drifted — fixture/chunking regression"
                    );
                    black_box(stats)
                });
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_cohort);
criterion_main!(benches);
