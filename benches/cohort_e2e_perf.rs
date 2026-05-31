//! End-to-end `.psp` → cohort-VCF throughput (Stages 3–6).
//!
//! The four per-stage iterator benches in `var_calling_perf.rs`
//! isolate one stage at a time against synthetic in-memory `Vec`s.
//! This file exists for the complementary case the standing
//! PROJECT_STATUS item "Parallel-optimization integration perf
//! benches" calls out: a single timed region that drives the whole
//! cohort pipeline (DUST → grouper → per-group merger → posterior
//! engine → VCF writer), so the parallelisation-tuning pass has a
//! ground truth for end-to-end scaling artefacts the per-stage
//! benches can't catch.
//!
//! Two bench-group families:
//!
//! - `cohort_e2e_core/*` — times
//!   [`drive_cohort_pipeline`][drive_cohort_pipeline]. Inputs are
//!   pre-decoded per-sample record vectors (cloned per iteration into
//!   per-sample iterators feeding the merger), the reference fetcher
//!   is built once per fixture, and per-stage configs are CLI
//!   defaults. The PSP-reader open + decode and the FASTA-MD5
//!   verification are intentionally outside the timed region — this
//!   group is the one the parallelisation-tuning pass should target.
//!   Sub-groups:
//!     - `scaling_samples` — vary N ∈ {10, 64, 256} at fixed
//!       region length, all cores. Answers H7 from the perf review.
//!     - `scaling_threads` — vary T ∈ {1, 2, 4, max} at fixed N, L.
//!       Implemented via a local `rayon::ThreadPool` per measurement
//!       and `pool.install(|| drive_cohort_pipeline(...))`. This is
//!       the headline lever for the deferred rayon-over-records work.
//!     - `scaling_region` — vary L ∈ {1 000, 5 000, 20 000} at
//!       fixed N, all cores. Reports throughput-per-position.
//!
//! - `cohort_e2e_full/*` — times
//!   [`run_var_calling`][run_var_calling], the public CLI entry point.
//!   Includes PSP open + header validate + FASTA MD5 cross-check +
//!   per-stage construction in addition to the core drive. **Cannot
//!   sweep thread count within one `cargo bench` invocation** because
//!   `configure_rayon_pool` / `rayon::ThreadPoolBuilder::build_global`
//!   is once-per-process; the first call wins and every subsequent
//!   call is a silent no-op (the `AtomicBool` guard in
//!   `pop_var_caller::common::configure_rayon_pool`). To bench the
//!   full path at a different thread count, run `cargo bench --bench
//!   cohort_e2e_perf cohort_e2e_full -- ...` with `RAYON_NUM_THREADS=N`
//!   in the environment for each run.
//!
//! Synthetic-fixture shape: one contig, deterministic random
//! reference bases (xorshift seed), every position covered REF in
//! every sample (depth 30), one het SNP every 50 positions in roughly
//! two-thirds of the cohort. The MD5 in the `.psp` header matches
//! the on-disk FASTA bytes so [`run_var_calling`]'s `M5` check
//! passes. Records are built once and reused: the on-disk `.psp` for
//! the full path, the in-memory vector for the core path.
//!
//! Bench-CPU-governor caveat (memory: `feedback_bench_cpu_governor`):
//! the host runs with the `powersave` governor and numbers drift
//! ~2× across sessions. Trust back-to-back runs over absolute
//! numbers; criterion comparisons (`--save-baseline` + `--baseline`)
//! within a single session are reliable.
//!
//! [drive_cohort_pipeline]: pop_var_caller::var_calling::from_bam::pipeline::drive_cohort_pipeline
//! [run_var_calling]: pop_var_caller::pop_var_caller::var_calling::run_var_calling

// Opt-in mimalloc global allocator (cargo bench --features alloc-mimalloc).
#[cfg(feature = "alloc-mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::collections::BTreeMap;
use std::fs::File;
use std::hint::black_box;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread::available_parallelism;
use std::time::{Duration, Instant};

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use md5::{Digest, Md5};
use tempfile::TempDir;

use pop_var_caller::fasta::StreamingChromRefFetcher;
use pop_var_caller::pileup_record::{AlleleObservation, AlleleSupportStats, PileupRecord};
use pop_var_caller::pop_var_caller::cli::shared_args::CohortPipelineArgs;
use pop_var_caller::var_calling::from_bam::pipeline::{
    CohortPipelineParams, DEFAULT_MIN_ALT_OBS_PER_SAMPLE, DEFAULT_MIN_MAPQ_DIFF_T,
    DEFAULT_MIN_QUAL_PHRED, drive_cohort_pipeline,
};
use pop_var_caller::pop_var_caller::var_calling::{
    VarCallingArgs, VarCallingCliError, run_var_calling,
};
use pop_var_caller::psp::PspReadError;
use pop_var_caller::psp::header::{
    ChromosomeEntry, ParsedChromosome, WriterHeader, WriterProvenance,
};
use pop_var_caller::psp::writer::PspWriter;
use pop_var_caller::var_calling::dust_filter::{
    DEFAULT_DUST_THRESHOLD, DEFAULT_DUST_WINDOW, DustFilterConfig,
};
use pop_var_caller::var_calling::per_group_merger::{
    DEFAULT_BATCH_SIZE, DEFAULT_MAX_ALLELES_LH_CALC, DEFAULT_MAX_ALLELES_PER_RECORD,
    DEFAULT_PLOIDY, PerGroupMergerConfig, SharedRefFetcher,
};
use pop_var_caller::var_calling::per_position_merger::PerPositionMerger;
use pop_var_caller::var_calling::posterior_engine::{
    DEFAULT_COMPOUND_ALT_PSEUDOCOUNT, DEFAULT_CONVERGENCE_THRESHOLD,
    DEFAULT_INBREEDING_COEFFICIENT, DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MAX_GQ_PHRED,
    DEFAULT_MAX_ITERATIONS, DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT,
    PosteriorEngineConfig,
};
use pop_var_caller::var_calling::variant_grouping::{
    DEFAULT_MAX_VARIANT_GROUP_SPAN, GrouperConfig,
};
use pop_var_caller::vcf::{CohortMetadata, WriterConfig};

// ---------------------------------------------------------------------
// Bench dimensions
// ---------------------------------------------------------------------

/// Cohort-size sweep — answers H7 from the 2026-05-16 perf review.
const SAMPLE_SWEEP: &[usize] = &[10, 64, 256];

/// Region-length sweep (positions per .psp).
const REGION_SWEEP: &[u32] = &[1_000, 5_000, 20_000];

/// Fixed cohort size used by `scaling_region` and `scaling_threads`.
const FIXED_N_SAMPLES: usize = 64;

/// Fixed region length used by `scaling_samples` and `scaling_threads`.
const FIXED_N_POSITIONS: u32 = 2_000;

/// One het SNP every `SNP_CADENCE` positions in ~2/3 of the cohort.
/// Smaller values mean more grouper / per-group merger work per
/// position; 50 is loose enough that REF positions dominate but
/// dense enough that variant groups exist in every region.
const SNP_CADENCE: u32 = 50;

const REF_DEPTH: u32 = 30;
const HET_REF_DEPTH: u32 = 15;
const HET_ALT_DEPTH: u32 = 15;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

// ---------------------------------------------------------------------
// Fixture
// ---------------------------------------------------------------------

/// Bundle of everything `drive_cohort_pipeline` consumes by value.
/// Returned by [`CohortFixture::core_inputs`] so the per-iter setup
/// for the `core` bench group has a single named return type.
type CoreInputs = (
    PerPositionMerger<std::vec::IntoIter<Result<PileupRecord, PspReadError>>>,
    CohortPipelineParams,
    CohortMetadata,
    WriterConfig,
);

/// One synthetic cohort, built once per (N, L) cell. Holds:
/// - on-disk FASTA + per-sample `.psp` files (consumed by the `full`
///   bench group through `run_var_calling`);
/// - pre-decoded per-sample `Vec<PileupRecord>` (cloned per iteration
///   in the `core` group as the merger's input iterators).
///
///
/// No shared fetcher: post-H1, `SharedRefFetcher` is
/// `Arc<dyn ChromRefFetcher + Send>` (not `+ Sync`), so it can't be held
/// in the fixture and reused across `pool.install` boundaries. Each
/// `core` iteration builds its own `StreamingChromRefFetcher::for_contig`
/// — same shape as `process_one_chromosome` in production. The per-iter
/// `for_contig` cost (one `open(2)` + small `.fai` parse) is negligible
/// against the pipeline it gates.
///
/// `_dir` is held to keep the `.psp` files and FASTA alive until the
/// fixture is dropped; bench code uses `dir()` for per-iter VCF output
/// paths.
struct CohortFixture {
    _dir: TempDir,
    dir_path: PathBuf,
    fasta: PathBuf,
    psp_files: Vec<PathBuf>,
    records: Vec<Vec<PileupRecord>>,
    sample_names: Vec<String>,
    chromosomes: Vec<ParsedChromosome>,
    n_positions: u32,
    n_samples: usize,
}

impl CohortFixture {
    fn build(n_samples: usize, n_positions: u32) -> Self {
        std::fs::create_dir_all("tmp").expect("mkdir tmp");
        let dir = TempDir::new_in("tmp").expect("tempdir");
        let dir_path = dir.path().to_path_buf();

        let contig_name = "chr1".to_string();
        let ref_seq = random_ref_bases(n_positions);
        let contig_md5 = md5_hex(&ref_seq);

        let fasta = write_fasta(&dir_path, &contig_name, &ref_seq);

        // Per-sample synthetic records — built once, reused across
        // iterations by clone for the core path and serialised to
        // `.psp` once for the full path.
        let records: Vec<Vec<PileupRecord>> = (0..n_samples)
            .map(|s| build_sample_records(n_positions, &ref_seq, s))
            .collect();

        let sample_names: Vec<String> = (0..n_samples).map(|s| format!("S{s:04}")).collect();

        let psp_files: Vec<PathBuf> = sample_names
            .iter()
            .zip(records.iter())
            .map(|(name, recs)| {
                let psp = dir_path.join(format!("{name}.psp"));
                write_psp(&psp, name, &contig_name, n_positions, &contig_md5, recs);
                psp
            })
            .collect();

        let chromosomes = vec![ParsedChromosome {
            name: contig_name,
            length: n_positions,
            md5: contig_md5,
        }];

        Self {
            _dir: dir,
            dir_path,
            fasta,
            psp_files,
            records,
            sample_names,
            chromosomes,
            n_positions,
            n_samples,
        }
    }

    /// Build a fresh per-iter input bundle for the core bench. The
    /// merger is constructed over `Vec::into_iter()` so the bench
    /// can hand it by value to `drive_cohort_pipeline`. Cloning the
    /// `Vec<PileupRecord>` per iteration is the small, deliberate
    /// price the core bench pays to keep the timed region pure.
    fn core_inputs(&self, vcf_out: PathBuf) -> CoreInputs {
        let iters: Vec<std::vec::IntoIter<Result<PileupRecord, PspReadError>>> = self
            .records
            .iter()
            .map(|recs| recs.iter().cloned().map(Ok).collect::<Vec<_>>().into_iter())
            .collect();
        let merger =
            PerPositionMerger::new(iters, self.sample_names.clone(), self.chromosomes.clone())
                .expect("PerPositionMerger::new");

        let params = CohortPipelineParams {
            no_complexity_filter: false,
            dust_cfg: DustFilterConfig::new(DEFAULT_DUST_WINDOW, DEFAULT_DUST_THRESHOLD)
                .expect("dust cfg"),
            grouper_cfg: GrouperConfig::new(DEFAULT_MAX_VARIANT_GROUP_SPAN).expect("grouper cfg"),
            per_group_cfg: PerGroupMergerConfig::new(
                DEFAULT_PLOIDY,
                DEFAULT_MAX_ALLELES_PER_RECORD,
                DEFAULT_MAX_ALLELES_LH_CALC,
                DEFAULT_BATCH_SIZE,
            )
            .expect("per-group cfg"),
            posterior_cfg: posterior_default_cfg(),
            min_qual_phred: DEFAULT_MIN_QUAL_PHRED,
            min_alt_obs_per_sample: DEFAULT_MIN_ALT_OBS_PER_SAMPLE,
            no_mapq_diff_filter: false,
            min_mapq_diff_t: DEFAULT_MIN_MAPQ_DIFF_T,
        };

        let metadata = CohortMetadata {
            sample_names: self.sample_names.clone(),
            contigs: self.chromosomes.clone(),
            tool_string: "pop_var_caller-bench".to_string(),
            command_line: "cargo bench --bench cohort_e2e_perf".to_string(),
        };
        let writer_cfg = WriterConfig::new(vcf_out).with_emit_gp(false);

        (merger, params, metadata, writer_cfg)
    }

    /// Build the `VarCallingArgs` the full bench hands to
    /// `run_var_calling`. Mirrors the integration-test defaults; the
    /// FASTA / `.psp` MD5s line up so the M5 check passes.
    fn full_args(&self, vcf_out: PathBuf) -> VarCallingArgs {
        VarCallingArgs {
            reference: self.fasta.clone(),
            output: vcf_out,
            // None → `configure_rayon_pool` accepts whatever the
            // global pool is (or installs the rayon default on first
            // call). The bench can't override per-iter because
            // `build_global` is once-per-process.
            threads: None,
            contamination_estimates: None,
            no_complexity_filter: false,
            // Legacy single-pull behaviour so bench numbers stay
            // comparable to the pre-rewrite shape.
            target_variants_per_chunk: 0,
            psp_files: self.psp_files.clone(),
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
}

fn posterior_default_cfg() -> PosteriorEngineConfig {
    PosteriorEngineConfig::new()
        .with_convergence_threshold(DEFAULT_CONVERGENCE_THRESHOLD)
        .expect("conv")
        .with_max_iterations(DEFAULT_MAX_ITERATIONS)
        .expect("iter")
        .with_ref_pseudocount(DEFAULT_REF_PSEUDOCOUNT)
        .expect("ref pseudo")
        .with_snp_alt_pseudocount(DEFAULT_SNP_ALT_PSEUDOCOUNT)
        .expect("snp pseudo")
        .with_indel_alt_pseudocount(DEFAULT_INDEL_ALT_PSEUDOCOUNT)
        .expect("indel pseudo")
        .with_compound_alt_pseudocount(DEFAULT_COMPOUND_ALT_PSEUDOCOUNT)
        .expect("compound pseudo")
        .with_fixation_index_default(DEFAULT_INBREEDING_COEFFICIENT)
        .expect("F")
        .with_max_gq_phred(DEFAULT_MAX_GQ_PHRED)
        .expect("max gq")
}

// ---------------------------------------------------------------------
// Fixture-building helpers
// ---------------------------------------------------------------------

/// Deterministic random reference bases via a small LCG so the bench
/// doesn't pull in `rand` and runs are byte-identical across sessions.
fn random_ref_bases(n_positions: u32) -> Vec<u8> {
    let mut state: u32 = 0x1234_5678;
    let mut out = Vec::with_capacity(n_positions as usize);
    for _ in 0..n_positions {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        out.push(BASES[((state >> 16) & 3) as usize]);
    }
    out
}

fn md5_hex(bytes: &[u8]) -> String {
    let digest = Md5::new().chain_update(bytes).finalize();
    let mut out = String::with_capacity(32);
    for b in digest {
        use std::fmt::Write as _;
        write!(&mut out, "{b:02x}").unwrap();
    }
    out
}

fn write_fasta(dir: &Path, contig_name: &str, seq: &[u8]) -> PathBuf {
    let fasta_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let mut fa = File::create(&fasta_path).expect("create fasta");
    writeln!(fa, ">{}", contig_name).unwrap();
    let header_len = (contig_name.len() + 2) as u64;
    fa.write_all(seq).unwrap();
    fa.write_all(b"\n").unwrap();
    let mut fai = File::create(&fai_path).expect("create fai");
    writeln!(
        fai,
        "{}\t{}\t{}\t{}\t{}",
        contig_name,
        seq.len(),
        header_len,
        seq.len(),
        seq.len() + 1
    )
    .unwrap();
    fasta_path
}

fn build_sample_records(n_positions: u32, ref_seq: &[u8], sample_idx: usize) -> Vec<PileupRecord> {
    let mut out = Vec::with_capacity(n_positions as usize);
    for pos in 1..=n_positions {
        let ref_base = ref_seq[(pos - 1) as usize];
        let is_snp = pos.is_multiple_of(SNP_CADENCE) && !sample_idx.is_multiple_of(3);
        let alleles = if is_snp {
            let mut alt = BASES[(sample_idx + pos as usize) & 3];
            if alt == ref_base {
                alt = BASES[(sample_idx + pos as usize + 1) & 3];
            }
            vec![
                AlleleObservation::new(
                    vec![ref_base],
                    AlleleSupportStats::new(
                        HET_REF_DEPTH,
                        f64::from(HET_REF_DEPTH) * -3.0,
                        HET_REF_DEPTH / 2,
                        HET_REF_DEPTH / 4,
                        HET_REF_DEPTH / 8,
                        0,
                        0,
                    ),
                    Vec::new(),
                ),
                AlleleObservation::new(
                    vec![alt],
                    AlleleSupportStats::new(
                        HET_ALT_DEPTH,
                        f64::from(HET_ALT_DEPTH) * -3.0,
                        HET_ALT_DEPTH / 2,
                        HET_ALT_DEPTH / 4,
                        HET_ALT_DEPTH / 8,
                        0,
                        0,
                    ),
                    Vec::new(),
                ),
            ]
        } else {
            vec![AlleleObservation::new(
                vec![ref_base],
                AlleleSupportStats::new(
                    REF_DEPTH,
                    0.0,
                    REF_DEPTH / 2,
                    REF_DEPTH / 4,
                    REF_DEPTH / 8,
                    0,
                    0,
                ),
                Vec::new(),
            )]
        };
        out.push(PileupRecord::new(0, pos, alleles));
    }
    out
}

fn write_psp(
    path: &Path,
    sample_name: &str,
    contig_name: &str,
    n_positions: u32,
    contig_md5: &str,
    records: &[PileupRecord],
) {
    let header = WriterHeader {
        format_version: (1, 0),
        sample: sample_name.to_string(),
        reference: "ref.fa".to_string(),
        created: "2026-05-19T00:00:00Z".parse().unwrap(),
        chromosomes: vec![ChromosomeEntry {
            name: contig_name.to_string(),
            length: n_positions,
            md5: contig_md5.to_string(),
        }],
        writer: WriterProvenance {
            tool: "pop_var_caller-bench".to_string(),
            version: "0.0.0".to_string(),
            subcommand: "pileup".to_string(),
            input_crams: vec![format!("{sample_name}.cram")],
            input_fasta: "ref.fa".to_string(),
            parameters: BTreeMap::new(),
        },
    };
    let file = File::create(path).expect("create psp");
    let sink = BufWriter::with_capacity(64 * 1024, file);
    let mut writer = PspWriter::new(sink, header).expect("psp writer");
    for record in records {
        writer.write_record(record).expect("write record");
    }
    let buf = writer.finish().expect("finish psp");
    buf.into_inner()
        .expect("buf into_inner")
        .sync_all()
        .expect("sync");
}

// ---------------------------------------------------------------------
// Core bench group — drive_cohort_pipeline in isolation
// ---------------------------------------------------------------------

fn bench_core_scaling_samples(c: &mut Criterion) {
    let mut group = c.benchmark_group("cohort_e2e_core/scaling_samples");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(8));

    for &n in SAMPLE_SWEEP {
        let fixture = CohortFixture::build(n, FIXED_N_POSITIONS);
        // Throughput: per-position elements × N samples gives a
        // per-cell-of-work rate that's comparable across N.
        group.throughput(Throughput::Elements(
            u64::from(fixture.n_positions) * fixture.n_samples as u64,
        ));
        group.bench_function(BenchmarkId::from_parameter(n), |b| {
            let mut counter = 0u64;
            b.iter_custom(|iters| run_core_iters(&fixture, &mut counter, iters, None));
        });
    }
    group.finish();
}

fn bench_core_scaling_region(c: &mut Criterion) {
    let mut group = c.benchmark_group("cohort_e2e_core/scaling_region");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(8));

    for &l in REGION_SWEEP {
        let fixture = CohortFixture::build(FIXED_N_SAMPLES, l);
        group.throughput(Throughput::Elements(
            u64::from(fixture.n_positions) * fixture.n_samples as u64,
        ));
        group.bench_function(BenchmarkId::from_parameter(l), |b| {
            let mut counter = 0u64;
            b.iter_custom(|iters| run_core_iters(&fixture, &mut counter, iters, None));
        });
    }
    group.finish();
}

fn bench_core_scaling_threads(c: &mut Criterion) {
    let mut group = c.benchmark_group("cohort_e2e_core/scaling_threads");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(8));

    let max_threads = available_parallelism()
        .map(std::num::NonZero::get)
        .unwrap_or(4);
    let thread_counts: Vec<usize> = [1usize, 2, 4, max_threads]
        .into_iter()
        .filter(|&t| t == 1 || t <= max_threads)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();

    let fixture = CohortFixture::build(FIXED_N_SAMPLES, FIXED_N_POSITIONS);
    let total_positions = u64::from(fixture.n_positions) * fixture.n_samples as u64;
    group.throughput(Throughput::Elements(total_positions));

    for t in thread_counts {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build()
            .expect("rayon pool");
        group.bench_function(BenchmarkId::from_parameter(t), |b| {
            let mut counter = 0u64;
            b.iter_custom(|iters| run_core_iters(&fixture, &mut counter, iters, Some(&pool)));
        });
    }
    group.finish();
}

/// Drive `iters` iterations of `drive_cohort_pipeline` against
/// `fixture`, returning the total *timed* duration. Per-iter setup
/// (cloning records, building merger + configs, picking a unique
/// VCF path) and per-iter cleanup (removing the VCF) are outside
/// the timed region. `pool` selects the rayon pool to install
/// around the drive call; `None` lets the per-group merger's
/// `into_par_iter` pick whatever pool is current (typically the
/// global default).
fn run_core_iters(
    fixture: &CohortFixture,
    counter: &mut u64,
    iters: u64,
    pool: Option<&rayon::ThreadPool>,
) -> Duration {
    let mut total = Duration::ZERO;
    for _ in 0..iters {
        *counter += 1;
        let vcf_out = fixture.dir_path.join(format!("core_{counter:08}.vcf"));
        let (merger, params, metadata, writer_cfg) = fixture.core_inputs(vcf_out.clone());

        // Build the fetcher per-iter inside the closure body. Post-H1
        // `SharedRefFetcher` is `Arc<dyn ChromRefFetcher + Send>` (not
        // `+ Sync`), so the Arc itself isn't `Send` — fetchers must be
        // constructed on the thread that uses them, exactly as
        // `process_one_chromosome` does in production
        // (cohort_driver.rs:472).
        let fasta_path = fixture.fasta.clone();
        let chrom_name = fixture.chromosomes[0].name.clone();
        let vcf_out_for_pipeline = vcf_out.clone();
        let start = Instant::now();
        let written = match pool {
            Some(p) => p.install(move || {
                let streaming = StreamingChromRefFetcher::for_contig(&fasta_path, &chrom_name)
                    .expect("StreamingChromRefFetcher::for_contig");
                #[allow(clippy::arc_with_non_send_sync)]
                let fetcher: SharedRefFetcher = Arc::new(streaming);
                drive_cohort_pipeline::<_, VarCallingCliError>(
                    0u32,
                    merger,
                    params,
                    fetcher,
                    &vcf_out_for_pipeline,
                    metadata,
                    writer_cfg,
                )
                .expect("drive_cohort_pipeline")
            }),
            None => {
                let streaming = StreamingChromRefFetcher::for_contig(&fasta_path, &chrom_name)
                    .expect("StreamingChromRefFetcher::for_contig");
                #[allow(clippy::arc_with_non_send_sync)]
                let fetcher: SharedRefFetcher = Arc::new(streaming);
                drive_cohort_pipeline::<_, VarCallingCliError>(
                    0u32,
                    merger,
                    params,
                    fetcher,
                    &vcf_out_for_pipeline,
                    metadata,
                    writer_cfg,
                )
                .expect("drive_cohort_pipeline")
            }
        };
        total += start.elapsed();
        black_box(written);
        let _ = std::fs::remove_file(&vcf_out);
    }
    total
}

// ---------------------------------------------------------------------
// Full bench group — run_var_calling end-to-end
// ---------------------------------------------------------------------

fn bench_full_scaling_samples(c: &mut Criterion) {
    let mut group = c.benchmark_group("cohort_e2e_full/scaling_samples");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(8));

    for &n in SAMPLE_SWEEP {
        let fixture = CohortFixture::build(n, FIXED_N_POSITIONS);
        group.throughput(Throughput::Elements(
            u64::from(fixture.n_positions) * fixture.n_samples as u64,
        ));
        group.bench_function(BenchmarkId::from_parameter(n), |b| {
            let mut counter = 0u64;
            b.iter_custom(|iters| run_full_iters(&fixture, &mut counter, iters));
        });
    }
    group.finish();
}

fn bench_full_scaling_region(c: &mut Criterion) {
    let mut group = c.benchmark_group("cohort_e2e_full/scaling_region");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(8));

    for &l in REGION_SWEEP {
        let fixture = CohortFixture::build(FIXED_N_SAMPLES, l);
        group.throughput(Throughput::Elements(
            u64::from(fixture.n_positions) * fixture.n_samples as u64,
        ));
        group.bench_function(BenchmarkId::from_parameter(l), |b| {
            let mut counter = 0u64;
            b.iter_custom(|iters| run_full_iters(&fixture, &mut counter, iters));
        });
    }
    group.finish();
}

fn run_full_iters(fixture: &CohortFixture, counter: &mut u64, iters: u64) -> Duration {
    let mut total = Duration::ZERO;
    for _ in 0..iters {
        *counter += 1;
        let vcf_out = fixture.dir_path.join(format!("full_{counter:08}.vcf"));
        let args = fixture.full_args(vcf_out.clone());

        let start = Instant::now();
        run_var_calling(&args).expect("run_var_calling");
        total += start.elapsed();
        let _ = std::fs::remove_file(&vcf_out);
    }
    total
}

// ---------------------------------------------------------------------
// Criterion plumbing
// ---------------------------------------------------------------------
//
// Mi24: regression threshold convention. Criterion's CI machinery is
// the threshold — `--baseline pre-fixes` flips each bench to
// `improved` / `no change` / `regressed` based on its own statistical
// model, no hand-picked percentage on top. The comments below pin the
// **intent** so a reviewer can tell whether a `regressed` verdict on a
// given group is expected or actionable.
//
// REGRESSION THRESHOLD: 10% per group is the team-accepted budget for
// any fix that does *not* itself buy back correctness — see Wave 6's
// performance check policy in the code-review-fixes skill. Correctness
// and security fixes are kept regardless of perf delta (the trade-off
// this fix workflow is explicitly designed to accept).

criterion_group!(
    // REGRESSION THRESHOLD: 10% per group.
    // Core driver wall (`drive_cohort_pipeline` in isolation, PSP open
    // + FASTA verify outside the timed region). Scales over samples,
    // region length, and threads.
    name = e2e_core;
    config = Criterion::default();
    targets =
        bench_core_scaling_samples,
        bench_core_scaling_region,
        bench_core_scaling_threads,
);

criterion_group!(
    // REGRESSION THRESHOLD: 10% per group.
    // Full `run_var_calling` wall — header validate + config build +
    // FASTA MD5 verify + driver. Cannot sweep thread count within
    // one `cargo bench` invocation (`configure_rayon_pool` is
    // once-per-process); use `RAYON_NUM_THREADS=N cargo bench` for
    // separate-process thread sweeps.
    name = e2e_full;
    config = Criterion::default();
    targets =
        bench_full_scaling_samples,
        bench_full_scaling_region,
);

criterion_main!(e2e_core, e2e_full);
