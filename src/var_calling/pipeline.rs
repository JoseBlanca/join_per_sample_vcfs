//! Driver / wiring — section orchestration (appendix §G).
//!
//! Wires the staged pipeline read → fold → compact → caller → writer and
//! exposes the cohort `.psp` → VCF entry point ([`run_var_calling`]) the CLI
//! invokes.
//!
//! **Parallel topology.** A bounded-`crossbeam-channel` pipeline inside a
//! [`std::thread::scope`], all hand-offs count-bounded for back-pressure (peak
//! tracks the queue capacities, not scheduling):
//!
//! 1. **read stage** ([`drive_read_stage`]) — owns the per-sample readers,
//!    decodes each sample's light-column segments ahead on the producer pool
//!    into per-sample bounded queues (depth [`READ_QUEUE_DEPTH`]);
//! 2. **fold/plan stage** (main thread, the [`CohortChunkIntegrator`]) — pulls
//!    decoded segments, folds the light columns to find the chunk's variable
//!    positions + safe-gap cut, fetches its REF span, and ships a `ChunkPlan`
//!    (still-compressed per-sample segments, shared by `Arc`) over a depth-2
//!    queue;
//! 3. **compact stage** ([`compact_plan`]) — inflates only those segments' heavy
//!    columns for the variable rows and ships a [`RawCohortChunk`];
//! 4. `W` **caller threads** ([`VariantCaller`]) pop chunks, call them, push
//!    [`CalledChunk`]s;
//! 5. one **writer thread** ([`VcfWriter`]) drains them, reordering by
//!    `chunk_order` via a `BTreeMap` so the emitted VCF is in genomic order
//!    regardless of how the callers finish.
//!
//! Splitting read / fold / compact into concurrent stages keeps the producer
//! pool fed — the light decode (the fold stage's dominant critical-path cost)
//! overlaps folding + compaction instead of blocking the fold thread. The
//! output is byte-identical for any worker count and any stage timing.

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::ops::Range;

use thiserror::Error;

use crate::fasta::{
    ChromRefFetchError, ChromRefFetcher, ManualEvictChromRefFetcher, StreamingChromRefFetcher,
};
use crate::pop_var_caller::var_calling::VarCallingArgs;
use crate::psp::{PspReadError, PspReader};
use crate::regions::{BedError, ContigBounds, Region, RegionSet};
use crate::var_calling::cohort_integration::{
    ChunkPlan, CohortChunkIntegrator, ProducerError, ReadMsg, RefFetchError, compact_plan,
    covered_intervals_for, drive_read_stage,
};
use crate::var_calling::contamination_estimation::ContaminationEstimates;
use crate::var_calling::dust_filter::{MIN_DUST_HALO, sdust_mask_for_span};
use crate::var_calling::per_group_merger::{
    DEFAULT_BATCH_SIZE, PerGroupMergerConfig, PerGroupMergerConfigError,
};
use crate::var_calling::per_position_merger::{PerPositionMergerError, check_chromosome_agreement};
use crate::var_calling::posterior_engine::{PosteriorEngineConfig, PosteriorEngineConfigError};
use crate::var_calling::sample_reader::SamplePspReader;
use crate::var_calling::types::{CalledChunk, RawCohortChunk};
use crate::var_calling::variant_caller::{CallerError, VariantCaller};
use crate::var_calling::variant_grouping::{GrouperConfig, GrouperConfigError};
use crate::var_calling::vcf_writer::{
    AlleleBalanceFilter, DownstreamFilters, VcfWriter, WriterError, WriterStats,
};
use crate::vcf::{CohortMetadata, WriterConfig};

use std::collections::BTreeMap;

use crate::paralog::ParalogModelParams;
use crate::sample_summary::SampleSummary;
use crate::var_calling::paralog_filter::calibrate::{CalibrationConfig, calibrate};
use crate::var_calling::paralog_filter::prepass::{HexpAccumulator, ParalogPrePass};
use crate::var_calling::paralog_filter::spill::{
    ParalogSpill, ParalogSpillRecord, ParalogSpillWriter, SpillError, WindowJoin,
};
use crate::var_calling::paralog_filter::window_gc::ReferenceWindowGc;
use crate::var_calling::paralog_filter::window_producer::WindowSpillBuilder;
use crate::var_calling::paralog_filter::write_pass::{
    WritePassError, paralog_provenance, run_write_pass,
};

/// Open-file buffer size — the shared default, referenced directly so it can't
/// drift from the canonical value.
const BUFFERED_IO_CAPACITY: usize = crate::pop_var_caller::common::DEFAULT_BUFFERED_IO_CAPACITY;
/// Default per-chunk variable-position target when `--target-variants-per-chunk`
/// is left at its `0` sentinel. Finer chunks load-balance the producer→caller
/// pipeline better (the producer is the wall floor, so smaller work-units keep
/// the callers fed sooner), and peak RSS is flat across chunk size — so the
/// tuning is purely a wall trade. A T=6 sweep (2026-06-08) confirmed the trend
/// continues below 256: at N=1000 wall fell monotonically 256→128→64
/// (208→201→196 s), and N=50 was flat within noise; 128 captures most of the
/// gain without pushing toward the per-chunk overhead (REF fetch / safe-gap
/// search) that eventually dominates at very small targets. Calls are
/// byte-identical regardless (chunks always cut at clean group boundaries).
/// Surfaced in the startup log so the effective value is visible.
const DEFAULT_TARGET_VARIANTS: u32 = 128;
/// Bounded-queue depth per worker on both hand-offs: peak resident chunks ≈
/// `QUEUE_DEPTH_PER_WORKER × n_workers`. `2` keeps every worker fed across one
/// hand-off without unbounded look-ahead — the back-pressure cap.
const QUEUE_DEPTH_PER_WORKER: usize = 2;
/// Depth of the bounded queue between the producer's fold/plan stage and its
/// compact stage. `2` lets the next plan's fold overlap the current plan's
/// heavy decode without reading unboundedly ahead — peak memory tracks this
/// depth (each in-flight plan pins its chunk's still-compressed segments).
const PLAN_QUEUE_DEPTH: usize = 2;
/// Depth of each per-sample queue between the producer's read stage and its
/// fold stage. `2` lets the read stage decode one segment ahead per sample
/// (overlapping decode with folding) without reading unboundedly ahead — peak
/// memory tracks this depth × N samples (each slot a still-compressed segment).
const READ_QUEUE_DEPTH: usize = 2;
/// Minimum producer-thread count at which the staged read→fold→compact pipeline
/// is used. Below it the stage threads' coordination overhead outweighs the
/// decode/fold overlap (the overlap needs spare pool capacity), so a single
/// inline producer is faster — see the `staged` decision in [`run_var_calling`].
const STAGED_MIN_THREADS: usize = 4;
/// DUST sub-span for the resident-buffer bound. The mask is byte-identical for
/// any sub-span (runs are coalesced across boundaries); this just caps RAM.
const DUST_SUBSPAN: u32 = 1_000_000;

/// Resolve the thread budget `N` from `--threads`: the explicit value, else
/// `available_parallelism`, floored at 1. The CLI / bench use this to size the
/// producer pool + caller threads (via [`resolve_split`]); both must agree on
/// the budget, so it lives in one place.
pub fn resolve_thread_budget(threads: Option<usize>) -> usize {
    threads
        .unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
        })
        .max(1)
}

/// Cohort size at or below which the producer (the decode-bound wall limiter)
/// should get extra threads: with few samples the per-chunk EM is cheap, so the
/// producer dominates and wants the budget. Measured (split sweep across
/// threads × cohort size, 2026-06-09): at ≤8 samples the optimum is
/// producer-heavy (≈`2N/3`), at 50 samples it is balanced (`N/2`). The penalty
/// for *over*-nudging a large cohort is steep (N=8/50-sample: P=6 → +67 % vs
/// balanced), so this stays at the largest cohort confirmed to still prefer the
/// nudge — it never extrapolates into the unmeasured 9–49 range (those default
/// to the safe balanced split).
const SMALL_COHORT_SAMPLES: usize = 8;

/// Partition the thread budget `n` into `(producer_pool, caller_threads)`,
/// `P + C = n`. Plan B of the thread-budget single-pool plan: a shared pool
/// let heavy non-preemptible EM tasks starve the decode-bound producer (~2×
/// slower than this partition at the same budget — measured), so the two stages
/// get dedicated threads that genuinely overlap.
///
/// Rule (measured): **balanced** `P = ⌊n/2⌋` for normal cohorts — the optimum
/// for large cohorts, where wall actually matters; the caller gets the rounding
/// since caller-starvation is the cliff there. **Nudged** `P = ⌈2n/3⌉` when
/// `n_samples ≤ `[`SMALL_COHORT_SAMPLES`] — EM is cheap, the producer wants the
/// threads (matches the sweep: N=4→3, N=8→6). `PVC_PRODUCER_THREADS` /
/// `PVC_CALLER_THREADS` override each side (the split sweep / power users).
/// Each side floored at 1, so `n = 1` yields `(1, 1)` (the irreducible
/// producer + caller pair; a true 1-thread total is not expressible).
pub fn resolve_split(n: usize, n_samples: usize) -> (usize, usize) {
    let p_default = if n_samples <= SMALL_COHORT_SAMPLES {
        (2 * n).div_ceil(3) // ⌈2n/3⌉ — producer-heavy
    } else {
        n / 2 // ⌊n/2⌋ — balanced (caller takes the rounding)
    };
    // Keep both sides ≥ 1: for n ≥ 2 that pins P into [1, n-1]; for n = 1 the
    // upper bound collapses to 1 and C floors back up to 1 → (1, 1).
    let p_default = p_default.clamp(1, n.saturating_sub(1).max(1));
    let env = |k: &str| {
        std::env::var(k)
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&v| v > 0)
    };
    let p = env("PVC_PRODUCER_THREADS").unwrap_or(p_default).max(1);
    let c = env("PVC_CALLER_THREADS")
        .unwrap_or_else(|| n.saturating_sub(p))
        .max(1);
    (p, c)
}

/// Log the process's live OS-thread count (Linux only; a no-op elsewhere),
/// called once the full pipeline topology is up: `main`, the producer pool, the
/// caller threads, the writer, and (on the staged path) the two coordinators.
/// The thread-budget contract (plan §6) is that this stays at `N + c` for a
/// small constant `c` (not the historical `~3N`); the
/// `thread_budget_integration` test parses this line and asserts it, guarding
/// against the budget silently drifting back. Read from `/proc/self/status` so
/// it reflects *actual* spawned threads, not a recomputed estimate.
fn log_live_thread_count() {
    #[cfg(target_os = "linux")]
    if let Ok(status) = std::fs::read_to_string("/proc/self/status")
        && let Some(n) = status.lines().find_map(|l| l.strip_prefix("Threads:"))
    {
        eprintln!("var-calling: live_threads={}", n.trim());
    }
}

/// Errors surfaced by the re-architected pipeline entry point.
#[derive(Debug, Error)]
pub enum PipelineError {
    #[error("opening / reading a .psp file: {0}")]
    Io(#[from] std::io::Error),
    #[error("decoding a .psp header: {0}")]
    Psp(#[from] PspReadError),
    #[error("cohort chromosome agreement: {0}")]
    ChromAgreement(#[from] PerPositionMergerError),
    #[error("invalid grouper configuration")]
    GrouperConfig(#[from] GrouperConfigError),
    #[error("invalid per-group-merger configuration")]
    MergerConfig(#[from] PerGroupMergerConfigError),
    #[error("invalid posterior-engine configuration")]
    PosteriorConfig(#[from] PosteriorEngineConfigError),
    #[error("invalid --regions BED")]
    Regions(#[from] BedError),
    #[error("reference fetch: {0}")]
    RefFetch(#[from] ChromRefFetchError),
    #[error("dust mask computation: {0}")]
    Dust(std::io::Error),
    #[error(transparent)]
    Producer(#[from] ProducerError),
    #[error(transparent)]
    Caller(#[from] CallerError),
    #[error(transparent)]
    Writer(#[from] WriterError),
    /// Paralog-filter spill I/O (write or read-back).
    #[error(transparent)]
    ParalogSpill(#[from] SpillError),
    /// The spill sink finished with chunks still buffered — a gap in
    /// `chunk_order` stalled the reorder drain (would silently truncate the
    /// spill). Mirrors [`WriterError::MissingChunks`].
    #[error("{count} chunk(s) never spilled — a gap in chunk_order stalled the paralog spill sink")]
    ParalogSpillGap { count: usize },
    /// The paralog write pass (spill → surviving VCF) failed.
    #[error(transparent)]
    ParalogWrite(#[from] WritePassError),
}

/// Cohort `.psp` → VCF entry point — the production pipeline the CLI's
/// [`run_var_calling`](crate::pop_var_caller::var_calling::run_var_calling)
/// invokes. Takes a [`VarCallingArgs`] (the CLI layer also runs the
/// reference-basename / FASTA-vs-`.psp` MD5 validation before calling here).
///
/// `contamination` (loaded by the CLI from `--contamination-estimates`, or
/// `None`) is wired into the EM's [`PosteriorEngineConfig`]. Returns the
/// run-level [`WriterStats`] for the CLI's stderr run summary.
///
/// `pool` is the producer-side rayon pool (built by the caller, sized to `P` —
/// the producer half of the [`resolve_split`] budget) that runs the decode /
/// fold / compact fan-out; `caller_threads` is `C`, the dedicated EM caller
/// count. `P + C = N` honours the `--threads` budget. The caller (CLI / bench)
/// owns pool construction because the criterion bench re-sizes per iteration.
///
/// (Reference-basename and FASTA-vs-`.psp` MD5 *validation* stay in the CLI
/// wrapper — they guard bad input but do not affect the emitted VCF.)
pub fn run_var_calling(
    args: &VarCallingArgs,
    contamination: Option<ContaminationEstimates>,
    pool: &rayon::ThreadPool,
    caller_threads: usize,
) -> Result<WriterStats, PipelineError> {
    let cohort = &args.cohort;

    // --- Open every .psp; validate chromosome agreement; collect metadata. ---
    let mut psp_readers: Vec<PspReader<BufReader<File>>> = Vec::with_capacity(args.psp_files.len());
    for path in &args.psp_files {
        let file = File::open(path)?;
        let buf = BufReader::with_capacity(BUFFERED_IO_CAPACITY, file);
        psp_readers.push(PspReader::new(buf)?);
    }
    let sample_names: Vec<String> = psp_readers
        .iter()
        .map(|r| r.header().sample.clone())
        .collect();
    let chromosomes = check_chromosome_agreement(&psp_readers)?;
    let chrom_names: Vec<String> = chromosomes.iter().map(|c| c.name.clone()).collect();
    let chrom_lengths: Vec<u32> = chromosomes.iter().map(|c| c.length).collect();
    let n_chromosomes = chromosomes.len() as u32;

    // --- Build every per-stage config from the args. Each builder returns a
    //     typed config error, surfaced through its own `PipelineError` variant
    //     (via `?` / `#[from]`) so the cause is preserved in the `source()`
    //     chain and the failing knob is identifiable. ---
    let grouper_cfg = GrouperConfig::new(cohort.var_group_max_span)?;
    let merger_cfg = PerGroupMergerConfig::new(
        cohort.ploidy,
        cohort.max_alleles_per_var,
        cohort.max_alleles_lh_calc,
        DEFAULT_BATCH_SIZE,
    )?;
    let posterior_cfg = PosteriorEngineConfig::new()
        .with_convergence_threshold(cohort.em_convergence_threshold)?
        .with_max_iterations(cohort.em_max_iterations)?
        .with_ref_pseudocount(cohort.ref_pseudocount)?
        .with_snp_alt_pseudocount(cohort.snp_alt_pseudocount)?
        .with_indel_alt_pseudocount(cohort.indel_alt_pseudocount)?
        .with_compound_alt_pseudocount(cohort.compound_alt_pseudocount)?
        .with_fixation_index_default(cohort.inbreeding_coefficient)?
        .with_max_gq_phred(cohort.max_gq_phred)?
        .with_contamination(contamination)?;

    let min_alt_obs = cohort.min_alt_obs_per_sample;
    let max_group_span = cohort.var_group_max_span;
    let target_variants = if args.target_variants_per_chunk == 0 {
        DEFAULT_TARGET_VARIANTS
    } else {
        args.target_variants_per_chunk
    };

    // --- Optional `--regions` BED: covered intervals get clipped to it
    //     (whole genome when absent — the identical byte-for-byte path). ---
    let region_set = match &args.regions {
        Some(bed) => {
            let contig_bounds: Vec<ContigBounds> = chromosomes
                .iter()
                .map(|c| ContigBounds {
                    name: &c.name,
                    length: c.length,
                })
                .collect();
            Some(RegionSet::from_bed_path(bed, &contig_bounds)?)
        }
        None => None,
    };

    // --- Hidden-paralog filter set-up (S1), while the `PspReader`s are still
    //     open: fit each sample's coverage model + obs_het from its `.psp`
    //     metadata section (the `SampleSummary`). Only when the filter is on —
    //     `paralog_requested` gates every filter cost so the off path is
    //     untouched. A sample whose `.psp` carries no summary (or fails to
    //     parse) is carried absent (S1), not fatal. ---
    let paralog_requested =
        PARALOG_WINDOWS_WIRED && !cohort.no_paralog_filter && cohort.paralog_fdr > 0.0;
    let paralog_prepass = if paralog_requested {
        let summaries: Vec<Option<SampleSummary>> = psp_readers
            .iter()
            .map(|r| {
                r.metadata()
                    .and_then(|bytes| SampleSummary::from_toml_bytes(bytes).ok())
            })
            .collect();
        Some(ParalogPrePass::fit(
            &summaries,
            &crate::paralog::CoverageFitConfig::default(),
        ))
    } else {
        None
    };

    // --- Per-sample segment readers. The (chrom_id, region_start, region_end)
    //     = (0, 1, 1) triple is an inert placeholder: the read stage calls
    //     `reset()` on every reader at the first interval before the first
    //     read, so these values are never observed. ---
    let readers: Vec<SamplePspReader<BufReader<File>>> = psp_readers
        .into_iter()
        .map(|r| SamplePspReader::new(r, 0, 1, 1))
        .collect();
    let n_samples = readers.len();

    // --- Whole-cohort interval schedule, computed once from the in-memory
    //     block indexes (no I/O) and shared by the read stage and the fold
    //     loop so both walk the covered intervals in identical order. With
    //     `--regions` each interval is clipped to the analysis regions (the
    //     byte-identical whole-genome path when absent). ---
    let schedule: Vec<(u32, Range<u32>)> = {
        let mut sched = Vec::new();
        for chrom_id in 0..n_chromosomes {
            let mut intervals = covered_intervals_for(&readers, chrom_id, max_group_span);
            if let Some(rs) = &region_set {
                intervals = restrict_intervals_to_regions(&intervals, rs.regions_for(chrom_id));
            }
            for interval in intervals {
                sched.push((chrom_id, interval));
            }
        }
        sched
    };

    // --- The three sections. ---
    let metadata = CohortMetadata {
        sample_names: sample_names.clone(),
        contigs: chromosomes.clone(),
        tool_string: format!("pop_var_caller {}", env!("CARGO_PKG_VERSION")),
        command_line: crate::pop_var_caller::common::current_command_line(),
        // Filled by the paralog write pass (S6) when the filter runs; empty
        // here keeps the single-pass header byte-identical.
        paralog_provenance: String::new(),
    };
    let writer_config = WriterConfig::new(args.output.clone()).with_emit_gp(cohort.emit_gp);
    let filters = DownstreamFilters {
        min_qual_phred: cohort.min_qual_phred,
        allele_balance: if cohort.no_allele_balance_filter {
            AlleleBalanceFilter::Off
        } else {
            AlleleBalanceFilter::On {
                min_log_lr: cohort.min_allele_balance_log_lr,
                concentration: cohort.allele_balance_concentration,
            }
        },
    };
    // The main-pass sink: the direct VCF writer (single-pass) or the paralog
    // spill sink (two-pass). Only the sink differs — the producer→callers
    // topology below is shared. For two-pass, `metadata` / `writer_config` are
    // kept (not moved) to build the final writer *after* calibration, when they
    // carry the paralog provenance; the single-pass branch clones them so they
    // survive for the borrow checker either way.
    // Two ephemeral spills when the filter is on (arch §5, Approach A): the
    // **record spill** (per-locus caller output, written by the sink thread) and
    // the **window spill** (per variant-window coverage, written inline in the
    // fold loop by `window_builder`). Both live in the output's scratch directory
    // and are deleted when these values drop, on success and on any `?` alike.
    let (paralog_spill, paralog_window_spill) = if paralog_requested {
        let scratch_dir = args
            .output
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| std::path::Path::new("."));
        (
            Some(ParalogSpill::create_in(scratch_dir)?),
            Some(ParalogSpill::create_in(scratch_dir)?),
        )
    } else {
        (None, None)
    };
    // The window-depth accumulator's tile width must match the coverage
    // histogram's (arch §4); the pre-pass carries it. Only meaningful when the
    // filter is on (a pre-pass exists).
    let paralog_window_bp = paralog_prepass
        .as_ref()
        .map(|p| p.window_bp())
        .unwrap_or(crate::sample_summary::DEFAULT_GC_WINDOW_BP);
    // The window-spill builder, fed inline in the fold loop (§4). `None` off the
    // filter path. Owns a writer handle to the window spill file.
    let mut window_builder: Option<WindowSpillBuilder<BufWriter<File>>> =
        match &paralog_window_spill {
            Some(ws) => Some(WindowSpillBuilder::new(
                n_samples,
                paralog_window_bp,
                ws.window_writer()?,
            )),
            None => None,
        };
    let sink: Box<dyn ChunkSink> = match &paralog_spill {
        Some(spill) => Box::new(SpillSink::new(spill.writer()?)),
        None => Box::new(VcfWriterSink(VcfWriter::new(
            metadata.clone(),
            writer_config.clone(),
            filters,
        )?)),
    };
    let caller = VariantCaller::new(
        grouper_cfg,
        merger_cfg,
        posterior_cfg,
        min_alt_obs,
        sample_names.clone(),
    );
    // --- Parallel topology: producer → W callers → writer. ---
    //
    // Two bounded crossbeam hand-offs give back-pressure (peak ≈ queue cap).
    // The producer stays on the main thread (its per-sample readers + REF/dust
    // fetchers are thread-local); the callers share `&caller` (it is `Sync` —
    // stateless `call_chunk`), and the writer owns the (`Send`) `CohortVcfWriter`
    // and reorders the out-of-order `CalledChunk`s by `chunk_order`, so the
    // emitted VCF is byte-identical regardless of worker count.
    //
    // Plan B of the thread-budget single-pool plan: the budget `N` is split into
    // a `P`-thread producer pool (this `pool`, sized by the caller) + `C = N-P`
    // dedicated EM caller threads, `P + C = N` (see [`resolve_split`]). A single
    // shared pool was tried first but let heavy non-preemptible EM tasks starve
    // the decode-bound producer (~2× slower at the same budget); dedicated
    // threads for the two stages overlap genuinely. `producer_threads` (= the
    // pool's size `P`) drives the `staged` decision and the startup log.
    let producer_threads = pool.current_num_threads();
    let n_workers = caller_threads;
    let cap = (QUEUE_DEPTH_PER_WORKER * n_workers).max(1);
    // Straddler decode cache: on by default, off under `--low-memory` (trades
    // the cache's cohort-scaling RSS for re-decompressing cut-spanning segments).
    let cache_straddlers = !args.low_memory;
    // The staged read→fold→compact pipeline (the read + compact stage threads)
    // overlaps decode with folding, but only pays off when the producer pool has
    // spare capacity for the overlap. Below `STAGED_MIN_THREADS` its coordination
    // cost — three stages contending for a tiny pool — outweighs the gain
    // (measured: a regression at 2 threads), so fall back to a single inline
    // producer (decode + fold + compact on the main thread, no stage threads),
    // which matches `main`'s topology and keeps the straddler cache. Output is
    // byte-identical either way.
    let staged = producer_threads >= STAGED_MIN_THREADS;
    // Surface the resolved concurrency / chunk-sizing defaults so an operator
    // can recover what a running instance actually chose (the values default
    // implicitly from `--threads` / the `--target-variants-per-chunk` sentinel).
    eprintln!(
        "var-calling: producer_threads={producer_threads} workers={n_workers} queue_cap={cap} target_variants_per_chunk={target_variants} low_memory={} staged={staged}",
        args.low_memory
    );
    // Common hand-offs: producer → callers → writer (back-pressured).
    let (chunk_tx, chunk_rx) = crossbeam_channel::bounded::<RawCohortChunk>(cap);
    let (called_tx, called_rx) = crossbeam_channel::bounded::<CalledChunk>(cap);
    let caller = &caller;
    let producer_pool = pool;
    let schedule = &schedule;

    let sink_output = std::thread::scope(|scope| -> Result<SinkOutput, PipelineError> {
        // Sink thread: drain CalledChunks and hand them to the sink (single-pass
        // = reorder + filter + write VCF; two-pass = reorder + accumulate Hexp +
        // spill). The callers share `&caller` (it is `Sync` — stateless
        // `call_chunk`); the sink reorders out-of-order `CalledChunk`s by
        // `chunk_order`, so the output is byte-identical regardless of worker
        // count or producer topology.
        let writer_handle = scope.spawn(move || -> Result<SinkOutput, PipelineError> {
            let mut sink = sink;
            for called in called_rx {
                sink.accept(called)?;
            }
            sink.close()
        });

        // Caller threads: pop chunks, call, push results.
        let mut caller_handles = Vec::with_capacity(n_workers);
        for _ in 0..n_workers {
            let chunk_rx = chunk_rx.clone();
            let called_tx = called_tx.clone();
            caller_handles.push(scope.spawn(move || -> Result<(), PipelineError> {
                for chunk in chunk_rx {
                    let called = caller.call_chunk(chunk)?;
                    if called_tx.send(called).is_err() {
                        break; // writer gone (it errored) — stop
                    }
                }
                Ok(())
            }));
        }
        // Main holds no caller/writer channel handles, so they close once the
        // producer side / callers finish.
        drop(chunk_rx);
        drop(called_tx);

        // Producer side — yields the producer-side first error (`None` on
        // success). Two topologies (see the `staged` decision above):
        let mut first_err: Option<PipelineError> = if staged {
            // STAGED: read → fold → compact, three internal stages joined by
            // bounded queues. The read stage decodes light columns ahead into
            // per-sample queues; the fold stage (main) pulls them, finds the
            // chunk's variable positions + cut, and ships a `ChunkPlan` over a
            // depth-2 queue; the compact stage inflates the heavy columns. The
            // queues let each stage's work overlap the others — the win on a
            // pool with spare capacity.
            let (seg_txs, seg_rxs): (
                Vec<crossbeam_channel::Sender<ReadMsg>>,
                Vec<crossbeam_channel::Receiver<ReadMsg>>,
            ) = (0..n_samples)
                .map(|_| crossbeam_channel::bounded::<ReadMsg>(READ_QUEUE_DEPTH))
                .unzip();
            let (wake_tx, wake_rx) = crossbeam_channel::bounded::<()>(1);
            let mut producer = CohortChunkIntegrator::<BufReader<File>>::new_streaming(
                seg_rxs,
                wake_tx,
                min_alt_obs,
                target_variants,
            );
            let (plan_tx, plan_rx) = crossbeam_channel::bounded::<ChunkPlan>(PLAN_QUEUE_DEPTH);

            // Compact stage: inflate each plan'd chunk's heavy columns (fanning
            // out across the N samples on `producer_pool`) and forward the
            // cohort chunk. Owns `chunk_tx`, so the caller queue closes when the
            // plan queue drains.
            let compact_handle = scope.spawn(move || -> Result<(), PipelineError> {
                for plan in plan_rx {
                    let chunk = producer_pool.install(|| compact_plan(plan, cache_straddlers))?;
                    if chunk_tx.send(chunk).is_err() {
                        break; // callers gone; real error surfaces on join
                    }
                }
                Ok(())
            });

            // Read stage: owns the readers, walks `schedule` decoding ahead into
            // the per-sample queues. Owns `seg_txs`, so they close when drained.
            let read_handle = scope.spawn(move || -> Result<(), PipelineError> {
                drive_read_stage(
                    readers,
                    schedule,
                    &seg_txs,
                    &wake_rx,
                    producer_pool,
                    READ_QUEUE_DEPTH,
                )?;
                Ok(())
            });

            // Fold/plan stage on the main thread (its `rebuild_fold` runs on
            // `producer_pool` because the whole loop is inside `install`).
            let mut ref_fetcher: Option<(u32, StreamingChromRefFetcher)> = None;
            let mut dust_fetcher: Option<(u32, ManualEvictChromRefFetcher)> = None;
            // Paralog window-spill state (S6c): the builder is fed inline from
            // each plan; `window_gc_chrom` tracks the contig its GC is built for,
            // `variant_scratch` is the reused variant-position buffer.
            let mut window_gc_chrom: Option<u32> = None;
            let mut variant_scratch: Vec<u32> = Vec::new();
            let produce = producer_pool.install(|| -> Result<(), PipelineError> {
                // Full topology (main + pool + callers + writer [+ coordinators])
                // is live by here — log the OS thread count for the §6 budget
                // guard (≤ N + c). See `log_live_thread_count`.
                log_live_thread_count();
                for (chrom_id, interval) in schedule {
                    let chrom_id = *chrom_id;
                    let mask = dust_mask_for(
                        &mut dust_fetcher,
                        chrom_id,
                        interval,
                        &args.reference,
                        &chrom_names,
                        &chrom_lengths,
                        args.no_complexity_filter,
                        cohort.complexity_window,
                        cohort.complexity_threshold,
                    )?;
                    begin_window_interval(
                        &mut window_builder,
                        &mut window_gc_chrom,
                        chrom_id,
                        &args.reference,
                        &chrom_names,
                        paralog_window_bp,
                    )?;
                    producer.begin_interval(chrom_id, interval.clone(), mask);
                    loop {
                        let plan = {
                            let mut fetch =
                                |start: u32, len: u32| -> Result<Vec<u8>, RefFetchError> {
                                    ref_fetch(
                                        &mut ref_fetcher,
                                        chrom_id,
                                        &args.reference,
                                        &chrom_names,
                                        start,
                                        len,
                                    )
                                };
                            producer.plan_chunk(&mut fetch)?
                        };
                        let Some(plan) = plan else { break };
                        if let Some(builder) = window_builder.as_mut() {
                            feed_window_builder(builder, &plan, &mut variant_scratch)?;
                        }
                        if plan_tx.send(plan).is_err() {
                            return Ok(()); // compact stage gone; error on join
                        }
                    }
                }
                Ok(())
            });
            // Close the plan queue; drop `producer` (holds the read-queue
            // receivers + wake sender) so the read stage is released even if it
            // is parked on `wake` or blocked sending — no hang on join.
            drop(plan_tx);
            drop(producer);

            let mut e = produce.err();
            if let Err(x) = read_handle.join().expect("read thread panicked") {
                e.get_or_insert(x);
            }
            if let Err(x) = compact_handle.join().expect("compact thread panicked") {
                e.get_or_insert(x);
            }
            e
        } else {
            // INLINE: one producer on the main thread — decode + fold + compact
            // in a single loop, shipping straight to the callers. No read /
            // compact stage threads (their coordination loses at low thread
            // counts). `produce_chunk` math, just without the stage queues.
            let mut producer = CohortChunkIntegrator::new(readers, min_alt_obs, target_variants);
            let mut ref_fetcher: Option<(u32, StreamingChromRefFetcher)> = None;
            let mut dust_fetcher: Option<(u32, ManualEvictChromRefFetcher)> = None;
            // Paralog window-spill state (S6c), as in the staged branch.
            let mut window_gc_chrom: Option<u32> = None;
            let mut variant_scratch: Vec<u32> = Vec::new();
            let produce = producer_pool.install(|| -> Result<(), PipelineError> {
                // Full topology (main + pool + callers + writer [+ coordinators])
                // is live by here — log the OS thread count for the §6 budget
                // guard (≤ N + c). See `log_live_thread_count`.
                log_live_thread_count();
                for (chrom_id, interval) in schedule {
                    let chrom_id = *chrom_id;
                    let mask = dust_mask_for(
                        &mut dust_fetcher,
                        chrom_id,
                        interval,
                        &args.reference,
                        &chrom_names,
                        &chrom_lengths,
                        args.no_complexity_filter,
                        cohort.complexity_window,
                        cohort.complexity_threshold,
                    )?;
                    begin_window_interval(
                        &mut window_builder,
                        &mut window_gc_chrom,
                        chrom_id,
                        &args.reference,
                        &chrom_names,
                        paralog_window_bp,
                    )?;
                    producer.begin_interval(chrom_id, interval.clone(), mask);
                    loop {
                        let plan = {
                            let mut fetch =
                                |start: u32, len: u32| -> Result<Vec<u8>, RefFetchError> {
                                    ref_fetch(
                                        &mut ref_fetcher,
                                        chrom_id,
                                        &args.reference,
                                        &chrom_names,
                                        start,
                                        len,
                                    )
                                };
                            producer.plan_chunk(&mut fetch)?
                        };
                        let Some(plan) = plan else { break };
                        if let Some(builder) = window_builder.as_mut() {
                            feed_window_builder(builder, &plan, &mut variant_scratch)?;
                        }
                        let chunk = compact_plan(plan, cache_straddlers)?;
                        if chunk_tx.send(chunk).is_err() {
                            return Ok(()); // callers gone; error on join
                        }
                    }
                }
                Ok(())
            });
            drop(chunk_tx); // close the caller queue
            produce.err()
        };

        // First error wins across the producer side, callers, and writer;
        // otherwise return the writer's run stats. `.join().expect(...)`
        // re-raises a worker-thread panic here (fatal-by-design); it cannot
        // deadlock (a panicking thread drops its channel handles, so the rest
        // observe close, drain, and join). `join()` only returns `Err` on panic.
        for handle in caller_handles {
            if let Err(e) = handle.join().expect("caller thread panicked") {
                first_err.get_or_insert(e);
            }
        }
        let writer_result = writer_handle.join().expect("writer thread panicked");
        match (first_err, writer_result) {
            (Some(e), _) | (None, Err(e)) => Err(e),
            (None, Ok(output)) => Ok(output),
        }
    })?;

    // Flush + close the window-spill builder now the fold has finished (drains
    // each contig's last windows). Done before the read passes so both spills
    // are complete on disk. `None` off the filter path.
    if let Some(builder) = window_builder.take() {
        builder.finish()?;
    }

    // Single-pass produced the VCF directly; two-pass now calibrates from the
    // spill and writes the surviving calls (S4 + S5). Each spill is read twice
    // (calibrate, then write) and deleted when its `ParalogSpill` drops at the
    // end of this function — on success and on any `?` error alike.
    match sink_output {
        SinkOutput::Vcf(stats) => Ok(stats),
        SinkOutput::Spill(hexp) => {
            let spill = paralog_spill.expect("a spill sink implies a record spill file");
            let window_spill =
                paralog_window_spill.expect("a spill sink implies a window spill file");
            let prepass = paralog_prepass.expect("a spill sink implies a pre-pass");
            let window_bp = prepass.window_bp();
            eprintln!(
                "var-calling: paralog filter — {} samples fit; Hexp Σ2pq {:.1} over {} variant \
                 loci; record spill at {}, window spill at {}",
                prepass.fit_count(),
                hexp.expected_het_sum(),
                hexp.loci(),
                spill.path().display(),
                window_spill.path().display(),
            );
            // Hexp on the per-callable-position scale (R1 correction): the
            // representative callable is the cohort median; `.max(1)` guards a
            // degenerate zero (→ Hexp 0 → F 0, treat every sample as outbred).
            let callable_ref = prepass.callable_reference().unwrap_or(1).max(1);
            let hexp = hexp.finish(callable_ref);
            let params = ParalogModelParams::default();
            let cfg = CalibrationConfig::default();

            let calibration = {
                let mut reader = spill.reader()?;
                let mut windows = WindowJoin::new(window_spill.window_reader()?)?;
                calibrate(
                    &prepass,
                    hexp,
                    &mut reader,
                    &mut windows,
                    window_bp,
                    &params,
                    cohort.paralog_fdr,
                    &cfg,
                )?
            };
            if !calibration.prior.converged {
                eprintln!(
                    "var-calling: paralog filter — the prior did not converge (or no scorable \
                     loci); using fallback π = {:.4}",
                    calibration.prior.prior_probability,
                );
            }

            // Struct-update spread so a new `CohortMetadata` field is a compile
            // error here rather than silently defaulting.
            let md = CohortMetadata {
                paralog_provenance: paralog_provenance(&calibration),
                ..metadata
            };
            let writer = VcfWriter::new(md, writer_config, filters)?;
            let mut reader = spill.reader()?;
            let mut windows = WindowJoin::new(window_spill.window_reader()?)?;
            let stats = run_write_pass(
                &mut reader,
                &mut windows,
                window_bp,
                &prepass,
                hexp,
                &calibration,
                &params,
                cfg.min_samples,
                writer,
            )?;
            Ok(stats)
        }
    }
}

/// Whether the producer computes per-window depth yet (S6c). Until it does, the
/// two-pass paralog path would spill records with no coverage signal and score
/// nothing — so the filter is held **off** regardless of `--paralog-fdr`. S6c
/// flips this to `true` (and this const is removed) once the windows are wired.
const PARALOG_WINDOWS_WIRED: bool = false;

/// A consumer of the caller→writer [`CalledChunk`] stream. Lets the main pass
/// drive either the direct VCF writer (single-pass) or the paralog spill sink
/// (two-pass) through the same producer→callers→sink topology; only the sink
/// differs. `Send` because it is moved into the sink thread.
trait ChunkSink: Send {
    /// Accept one chunk (reordering by `chunk_order` internally, as the writer
    /// does), so the sink sees records in genomic order.
    fn accept(&mut self, chunk: CalledChunk) -> Result<(), PipelineError>;
    /// Finish the sink, yielding its output. Takes `Box<Self>` so it can be
    /// called on the boxed trait object the sink thread owns.
    fn close(self: Box<Self>) -> Result<SinkOutput, PipelineError>;
}

/// What a [`ChunkSink`] produced: either the finished VCF stats (single-pass) or
/// the finished `Hexp` accumulator (two-pass; the spill is on disk).
#[derive(Debug)]
enum SinkOutput {
    /// Single-pass: the writer already produced the final VCF.
    Vcf(WriterStats),
    /// Two-pass: the spill is written; carries the finished `Hexp` accumulator
    /// (its `Σ2pq`), which the caller divides by the callable reference.
    Spill(HexpAccumulator),
}

/// Single-pass sink: the existing [`VcfWriter`], unchanged. `accept`/`close`
/// monomorphise to the same `handle`/`finish` the direct path always ran, so the
/// emitted VCF is byte-identical to the pre-filter pipeline.
struct VcfWriterSink(VcfWriter);

impl ChunkSink for VcfWriterSink {
    fn accept(&mut self, chunk: CalledChunk) -> Result<(), PipelineError> {
        self.0.handle(chunk)?;
        Ok(())
    }
    fn close(self: Box<Self>) -> Result<SinkOutput, PipelineError> {
        Ok(SinkOutput::Vcf(self.0.finish()?))
    }
}

/// Two-pass sink: reorder the chunks (as the writer does), accumulate `Hexp`
/// from each record's cohort allele frequencies, and spill each record to the
/// **record spill**. The per-sample window coverage lives in the sibling window
/// spill (built inline in the fold loop by `window_builder`), joined to each
/// record by tile key in the read passes (Approach A, S6c).
struct SpillSink {
    writer: ParalogSpillWriter<BufWriter<File>>,
    hexp: HexpAccumulator,
    /// Reorder buffer: chunks that arrived before their turn (mirrors
    /// [`VcfWriter`]'s reorder so the spill is in genomic order).
    reorder: BTreeMap<u64, CalledChunk>,
    /// The next `chunk_order` to spill (the gapless cursor).
    next_expected: u64,
}

impl SpillSink {
    fn new(writer: ParalogSpillWriter<BufWriter<File>>) -> Self {
        Self {
            writer,
            hexp: HexpAccumulator::new(),
            reorder: BTreeMap::new(),
            next_expected: 0,
        }
    }

    /// Accumulate `Hexp` and spill every record of one (in-order) chunk to the
    /// record spill. The record's window coverage is supplied separately by the
    /// window spill (joined by tile key in the read passes).
    fn spill_chunk(&mut self, chunk: CalledChunk) -> Result<(), PipelineError> {
        for record in chunk.records {
            self.hexp.observe(&record.allele_frequencies);
            self.writer.append(&ParalogSpillRecord { record })?;
        }
        Ok(())
    }
}

impl ChunkSink for SpillSink {
    fn accept(&mut self, chunk: CalledChunk) -> Result<(), PipelineError> {
        // Mirror the writer's reorder invariant: a chunk for an already-drained
        // slot would be silently buffered and never emitted (surfacing only as a
        // spurious gap at close). The producer is monotone, so this cannot fire.
        debug_assert!(
            chunk.chunk_order >= self.next_expected,
            "chunk_order {} already spilled (next_expected {})",
            chunk.chunk_order,
            self.next_expected,
        );
        self.reorder.insert(chunk.chunk_order, chunk);
        while let Some(ready) = self.reorder.remove(&self.next_expected) {
            self.spill_chunk(ready)?;
            self.next_expected += 1;
        }
        Ok(())
    }

    fn close(self: Box<Self>) -> Result<SinkOutput, PipelineError> {
        if !self.reorder.is_empty() {
            return Err(PipelineError::ParalogSpillGap {
                count: self.reorder.len(),
            });
        }
        self.writer.finish()?;
        Ok(SinkOutput::Spill(self.hexp))
    }
}

/// Build one contig's per-window reference GC for the paralog window join (S6c),
/// by streaming its FASTA bases through a fresh [`StreamingChromRefFetcher`] —
/// the same memory-flat serial walk the DUST mask uses (a ~1 MB sliding buffer,
/// never the whole contig). A separate fetcher from the fold's ref-span one,
/// because `iter_bases` resets the sliding buffer (the fetchers must not share).
fn build_window_gc(
    fasta: &std::path::Path,
    chrom_name: &str,
    window_bp: u32,
) -> Result<ReferenceWindowGc, PipelineError> {
    let fetcher = StreamingChromRefFetcher::for_contig(fasta, chrom_name)?;
    // `iter_bases` yields `Result<u8, _>`; capture the first fetch error and
    // stop the walk (map_while) rather than folding a bogus base into the GC.
    let mut fetch_err: Option<ChromRefFetchError> = None;
    let gc = ReferenceWindowGc::from_base_iter(
        fetcher.iter_bases()?.map_while(|b| match b {
            Ok(base) => Some(base),
            Err(e) => {
                fetch_err = Some(e);
                None
            }
        }),
        window_bp,
    );
    match fetch_err {
        Some(e) => Err(e.into()),
        None => Ok(gc),
    }
}

/// On a contig change, build the new contig's reference GC and open the window
/// builder's interval for it (draining the previous contig's windows). A no-op
/// off the filter path (`builder` is `None`) or within a contig. `gc_chrom`
/// tracks the contig the builder's GC is currently built for.
fn begin_window_interval<W: std::io::Write>(
    builder: &mut Option<WindowSpillBuilder<W>>,
    gc_chrom: &mut Option<u32>,
    chrom_id: u32,
    fasta: &std::path::Path,
    chrom_names: &[String],
    window_bp: u32,
) -> Result<(), PipelineError> {
    let Some(builder) = builder.as_mut() else {
        return Ok(());
    };
    if *gc_chrom != Some(chrom_id) {
        let gc = build_window_gc(fasta, &chrom_names[chrom_id as usize], window_bp)?;
        builder.begin_interval(chrom_id, gc)?;
        *gc_chrom = Some(chrom_id);
    }
    Ok(())
}

/// Feed one planned chunk's per-sample covered positions + variant positions
/// into the window-spill builder, then flush the windows now closed for the
/// whole cohort. Runs in the fold loop *before* the plan is compacted/sent, so it
/// reads the plan's light data without touching the caller path (byte-identity of
/// the caller output holds). `variant_scratch` is a reused buffer.
fn feed_window_builder<W: std::io::Write>(
    builder: &mut WindowSpillBuilder<W>,
    plan: &ChunkPlan,
    variant_scratch: &mut Vec<u32>,
) -> Result<(), SpillError> {
    plan.for_each_window_observation(|sample_idx, pos, depth| {
        builder.observe(sample_idx, pos, depth);
    });
    variant_scratch.clear();
    variant_scratch.extend(plan.chunk_variant_positions());
    builder.mark_variant_positions(variant_scratch);
    // Flush up to the cohort-fed frontier (this chunk's exclusive `cut`), not
    // each sample's last-covered tile — keeps `pending` bounded when a sample
    // has a coverage gap (S6c review).
    builder.flush_ready(plan.cut())
}

/// Fetch REF bytes for `chrom_id`, building a fresh per-contig
/// [`StreamingChromRefFetcher`] when the chromosome changes (fetches within a
/// chromosome are monotonic-forward across chunks, satisfying the fetcher's
/// sliding-window contract). The typed [`ChromRefFetchError`] is boxed (not
/// stringified) so the producer's [`ProducerError::Ref`] keeps it in the
/// `source()` chain.
fn ref_fetch(
    cache: &mut Option<(u32, StreamingChromRefFetcher)>,
    chrom_id: u32,
    fasta: &std::path::Path,
    chrom_names: &[String],
    start: u32,
    len: u32,
) -> Result<Vec<u8>, RefFetchError> {
    if cache.as_ref().map(|(c, _)| *c) != Some(chrom_id) {
        let f = StreamingChromRefFetcher::for_contig(fasta, &chrom_names[chrom_id as usize])?;
        *cache = Some((chrom_id, f));
    }
    // UNREACHABLE: `cache` is `Some` on both branches above — either the `if`
    // just set it, or it already matched `chrom_id`.
    let (_, fetcher) = cache.as_mut().expect("just set");
    // `ChromRefFetcher::fetch` (trait) returns an owned `Vec<u8>`.
    Ok(fetcher.fetch(start, len)?)
}

/// Compute the interval's DUST mask (empty when complexity filtering is off),
/// reusing a per-contig [`ManualEvictChromRefFetcher`].
// Threads the per-contig fetcher cache + reference + chrom name/length tables
// + the three complexity knobs; these are distinct concerns, not a cohesive
// struct, and the fn has two call sites that pass different cache slots (Mi6).
#[allow(clippy::too_many_arguments)]
fn dust_mask_for(
    cache: &mut Option<(u32, ManualEvictChromRefFetcher)>,
    chrom_id: u32,
    interval: &Range<u32>,
    fasta: &std::path::Path,
    chrom_names: &[String],
    chrom_lengths: &[u32],
    no_complexity_filter: bool,
    window: u32,
    threshold: u32,
) -> Result<Vec<Range<u32>>, PipelineError> {
    if no_complexity_filter {
        return Ok(Vec::new());
    }
    if cache.as_ref().map(|(c, _)| *c) != Some(chrom_id) {
        let f = ManualEvictChromRefFetcher::for_contig(fasta, &chrom_names[chrom_id as usize])?;
        *cache = Some((chrom_id, f));
    }
    // UNREACHABLE: `cache` is `Some` on both branches above (just set, or it
    // already matched `chrom_id`).
    let (_, fetcher) = cache.as_mut().expect("just set");
    dust_mask_for_interval(
        fetcher,
        interval,
        chrom_lengths[chrom_id as usize],
        window,
        threshold,
        DUST_SUBSPAN,
    )
    .map_err(PipelineError::Dust)
}

/// The interval's low-complexity (sdust) mask. Scans the interval in `subspan`
/// windows (evicting the buffer between them to bound resident RAM) and
/// coalesces runs across the window boundaries, so the result equals a single
/// whole-interval scan (byte-identical regardless of `subspan`).
fn dust_mask_for_interval(
    fetcher: &mut ManualEvictChromRefFetcher,
    interval: &Range<u32>,
    chrom_length: u32,
    window: u32,
    threshold: u32,
    subspan: u32,
) -> std::io::Result<Vec<Range<u32>>> {
    let subspan = subspan.max(1);
    let mut mask: Vec<Range<u32>> = Vec::new();
    let mut sub_start = interval.start;
    while sub_start < interval.end {
        let sub_end = sub_start.saturating_add(subspan).min(interval.end);
        let (sub_mask, buf_start) = sdust_mask_for_span(
            sub_start,
            sub_end,
            chrom_length,
            window,
            threshold,
            MIN_DUST_HALO,
            |s, l| {
                fetcher
                    .fetch(s, l)
                    .map(<[u8]>::to_vec)
                    .map_err(std::io::Error::other)
            },
        )?;
        for run in sub_mask {
            match mask.last_mut() {
                Some(last) if last.end == run.start => last.end = run.end,
                _ => mask.push(run),
            }
        }
        fetcher.evict_before(buf_start);
        sub_start = sub_end;
    }
    Ok(mask)
}

/// Clip the chromosome's covered intervals (1-based half-open) to its analysis
/// `regions` (1-based inclusive `[start, end]`), returning the overlap spans in
/// the half-open convention. Two-pointer sweep over the two sorted lists. With
/// no `--regions` this is never called and the path is byte-identical.
fn restrict_intervals_to_regions(covered: &[Range<u32>], regions: &[Region]) -> Vec<Range<u32>> {
    let mut out = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < covered.len() && j < regions.len() {
        let c = &covered[i];
        let region_end_excl = regions[j].end.saturating_add(1);
        let lo = c.start.max(regions[j].start);
        let hi = c.end.min(region_end_excl);
        if lo < hi {
            out.push(lo..hi);
        }
        if c.end <= region_end_excl {
            i += 1;
        } else {
            j += 1;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    // The `&[lo..hi]` covered-interval literals are intentional single-element
    // arrays of `Range`, not a misformed range expression.
    #![allow(clippy::single_range_in_vec_init)]
    use super::*;

    fn region(start: u32, end: u32) -> Region {
        Region {
            chrom_id: 0,
            start,
            end,
        }
    }

    /// The spill sink reorders out-of-order chunks by `chunk_order`, spills the
    /// records in genomic order, accumulates `Hexp` from each record's cohort
    /// allele frequencies, and yields the finished accumulator — the two-pass
    /// main-pass sink's contract.
    #[test]
    fn spill_sink_reorders_accumulates_hexp_and_spills() {
        use crate::var_calling::paralog_filter::spill::ParalogSpill;
        use crate::var_calling::paralog_filter::test_support::normal_locus;
        use crate::var_calling::types::CallStats;

        let n = 4;
        // Three biallelic-SNP records (allele_frequencies [0.5, 0.5] → 2pq 0.5).
        let rec = |pos: u32| normal_locus(pos, n).record;
        let chunk = |order: u64, pos: u32| CalledChunk {
            chunk_order: order,
            records: vec![rec(pos)],
            stats: CallStats::default(),
        };

        let dir = tempfile::tempdir().unwrap();
        let spill = ParalogSpill::create_in(dir.path()).unwrap();
        let mut sink: Box<dyn ChunkSink> = Box::new(SpillSink::new(spill.writer().unwrap()));

        // Deliver out of order: 2, 0, 1 (genomic order is 100, 200, 300).
        sink.accept(chunk(2, 300)).unwrap();
        sink.accept(chunk(0, 100)).unwrap();
        sink.accept(chunk(1, 200)).unwrap();
        let out = sink.close().unwrap();
        let SinkOutput::Spill(hexp) = out else {
            panic!("spill sink must yield SinkOutput::Spill");
        };

        // Hexp accumulator saw all three records' 2pq (0.5 each).
        assert_eq!(hexp.loci(), 3);
        assert!((hexp.expected_het_sum() - 1.5).abs() < 1e-9);

        // The spill holds the records in genomic order (reorder worked).
        let mut reader = spill.reader().unwrap();
        let positions: Vec<u32> = std::iter::from_fn(|| reader.next_record())
            .map(|r| r.unwrap().record.locus.start)
            .collect();
        assert_eq!(positions, vec![100, 200, 300]);
    }

    /// A spill sink closed with a gap in `chunk_order` (a chunk never arrived)
    /// errors rather than silently truncating the spill — mirrors the writer's
    /// gapless invariant.
    #[test]
    fn spill_sink_close_errors_on_gap() {
        use crate::var_calling::paralog_filter::spill::ParalogSpill;
        use crate::var_calling::paralog_filter::test_support::normal_locus;
        use crate::var_calling::types::CallStats;

        let dir = tempfile::tempdir().unwrap();
        let spill = ParalogSpill::create_in(dir.path()).unwrap();
        let mut sink: Box<dyn ChunkSink> = Box::new(SpillSink::new(spill.writer().unwrap()));
        // chunk_order 0 never delivered; 1 stays buffered behind the gap.
        sink.accept(CalledChunk {
            chunk_order: 1,
            records: vec![normal_locus(100, 4).record],
            stats: CallStats::default(),
        })
        .unwrap();
        match sink.close() {
            Err(PipelineError::ParalogSpillGap { count }) => assert_eq!(count, 1),
            other => panic!("expected ParalogSpillGap, got {other:?}"),
        }
    }

    #[test]
    fn restrict_intervals_clips_region_inside_one_interval() {
        // Region [12,15] (1-based inclusive) inside covered [10,20): the
        // half-open overlap is [12, 16).
        let got = restrict_intervals_to_regions(&[10..20], &[region(12, 15)]);
        assert_eq!(got, vec![12..16]);
    }

    #[test]
    fn restrict_intervals_straddling_region_splits_across_intervals() {
        // Region [18,34] straddles the gap between covered [10,20) and
        // [30,40): clipped to [18,20) ∪ [30,35).
        let got = restrict_intervals_to_regions(&[10..20, 30..40], &[region(18, 34)]);
        assert_eq!(got, vec![18..20, 30..35]);
    }

    #[test]
    fn restrict_intervals_touching_exclusive_end_is_empty() {
        // Region [20,25]: covered [10,20) is half-open, so position 20 is not
        // covered — no overlap.
        let got = restrict_intervals_to_regions(&[10..20], &[region(20, 25)]);
        assert!(got.is_empty(), "got {got:?}");
    }

    #[test]
    fn restrict_intervals_region_left_of_covered_is_empty() {
        let got = restrict_intervals_to_regions(&[10..20], &[region(1, 5)]);
        assert!(got.is_empty(), "got {got:?}");
    }

    #[test]
    fn restrict_intervals_empty_regions_yields_empty() {
        assert!(restrict_intervals_to_regions(&[10..20, 30..40], &[]).is_empty());
    }

    #[test]
    fn restrict_intervals_multiple_regions_in_one_interval() {
        // Two disjoint regions inside one covered interval, each clipped.
        let got = restrict_intervals_to_regions(&[10..50], &[region(12, 15), region(30, 40)]);
        assert_eq!(got, vec![12..16, 30..41]);
    }

    #[test]
    fn restrict_intervals_advances_covered_on_end_equals_region_end_excl() {
        // Mi2: when a covered interval ends exactly at `region_end_excl`,
        // the `c.end <= region_end_excl` branch must advance the *covered*
        // pointer (not the region pointer): region [10,20] → end_excl 21,
        // covered [5,21) overlaps as [10,21); the following covered
        // interval [21,30) is past the region and yields nothing.
        let got = restrict_intervals_to_regions(&[5..21, 21..30], &[region(10, 20)]);
        assert_eq!(got, vec![10..21]);
    }

    #[test]
    fn restrict_intervals_one_region_spanning_three_covered_intervals() {
        // Mi2: a single region wide enough to cover three disjoint covered
        // intervals clips to all three (the covered pointer advances three
        // times against one region).
        let got = restrict_intervals_to_regions(&[10..20, 30..40, 50..60], &[region(5, 100)]);
        assert_eq!(got, vec![10..20, 30..40, 50..60]);
    }
}
