//! Driver / wiring — section orchestration (appendix §G).
//!
//! Wires the three sections producer → caller → writer and exposes the cohort
//! `.psp` → VCF entry point ([`run_var_calling`]) the CLI invokes.
//!
//! **Parallel topology.** A bounded-`crossbeam-channel` pipeline inside a
//! [`std::thread::scope`]: the producer ([`CohortChunkIntegrator`]) runs on
//! the main thread (its per-sample readers + REF/dust fetchers are
//! thread-local) and ships [`RawCohortChunk`]s; `W` caller threads
//! ([`VariantCaller`]) pop chunks, call them, and push [`CalledChunk`]s; one
//! writer thread ([`VcfWriter`]) drains them, reordering by `chunk_order` via a
//! `BTreeMap` so the emitted VCF is in genomic order regardless of how the
//! callers finish. Both hand-offs are count-bounded (depth ≈
//! `QUEUE_DEPTH_PER_WORKER × n_workers`) for back-pressure, so peak tracks the
//! queue capacity rather than scheduling. The output is byte-identical for any
//! worker count.

use std::fs::File;
use std::io::BufReader;
use std::ops::Range;

use thiserror::Error;

use crate::fasta::{
    ChromRefFetchError, ChromRefFetcher, ManualEvictChromRefFetcher, StreamingChromRefFetcher,
};
use crate::pop_var_caller::var_calling::VarCallingArgs;
use crate::psp::{PspReadError, PspReader};
use crate::regions::{BedError, ContigBounds, Region, RegionSet};
use crate::var_calling::cohort_integration::{CohortChunkIntegrator, ProducerError, RefFetchError};
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
    DownstreamFilters, MapqDiffFilter, VcfWriter, WriterError, WriterStats,
};
use crate::vcf::{CohortMetadata, WriterConfig};

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
/// DUST sub-span for the resident-buffer bound. The mask is byte-identical for
/// any sub-span (runs are coalesced across boundaries); this just caps RAM.
const DUST_SUBSPAN: u32 = 1_000_000;

/// Resolve `(producer_threads, caller_threads)` from the requested thread
/// budget `t`. The producer's rayon decode/compaction pool and the crossbeam
/// caller pool are *independent*; left both at `t` they oversubscribe a
/// `t`-core box ~2× (perf finding **H1**). This is the single point that sizes
/// them, so the producer pool can be capped below the caller count.
///
/// Default: both `t` — byte-identical and perf-neutral vs the pre-H1 behaviour
/// (the producer's `par_iter_mut` ran on the global pool sized to `t`). The
/// `PVC_PRODUCER_THREADS` / `PVC_CALLER_THREADS` environment variables override
/// each, used for the H1 split sweep (perf review §3.3); once the winning split
/// is chosen it replaces this default. Each value is floored at 1.
fn resolve_thread_split(t: usize) -> (usize, usize) {
    let from_env = |k: &str| {
        std::env::var(k)
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&n| n > 0)
    };
    let producer = from_env("PVC_PRODUCER_THREADS").unwrap_or(t).max(1);
    let callers = from_env("PVC_CALLER_THREADS").unwrap_or(t).max(1);
    (producer, callers)
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
    #[error("building the producer rayon thread pool")]
    ProducerPool(#[from] rayon::ThreadPoolBuildError),
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
/// (Reference-basename and FASTA-vs-`.psp` MD5 *validation* stay in the CLI
/// wrapper — they guard bad input but do not affect the emitted VCF.)
pub fn run_var_calling(
    args: &VarCallingArgs,
    contamination: Option<ContaminationEstimates>,
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

    // --- Per-sample segment readers. The (chrom_id, region_start, region_end)
    //     = (0, 1, 1) triple is an inert placeholder: the producer's
    //     `begin_interval` calls `reset()` on every reader before the first
    //     read, so these values are never observed. ---
    let readers: Vec<SamplePspReader<BufReader<File>>> = psp_readers
        .into_iter()
        .map(|r| SamplePspReader::new(r, 0, 1, 1))
        .collect();

    // --- The three sections. ---
    let metadata = CohortMetadata {
        sample_names: sample_names.clone(),
        contigs: chromosomes.clone(),
        tool_string: format!("pop_var_caller {}", env!("CARGO_PKG_VERSION")),
        command_line: crate::pop_var_caller::common::current_command_line(),
    };
    let writer_config = WriterConfig::new(args.output.clone()).with_emit_gp(cohort.emit_gp);
    let filters = DownstreamFilters {
        min_qual_phred: cohort.min_qual_phred,
        mapq_diff: if cohort.no_mapq_diff_filter {
            MapqDiffFilter::Off
        } else {
            MapqDiffFilter::On {
                min_t: cohort.min_mapq_diff_t,
            }
        },
    };
    let mut writer = VcfWriter::new(metadata, writer_config, filters)?;
    let caller = VariantCaller::new(
        grouper_cfg,
        merger_cfg,
        posterior_cfg,
        min_alt_obs,
        sample_names.clone(),
    );
    let mut producer = CohortChunkIntegrator::new(readers, min_alt_obs, target_variants);

    // --- Parallel topology: producer (this thread) → W callers → writer. ---
    //
    // Two bounded crossbeam hand-offs give back-pressure (peak ≈ queue cap).
    // The producer stays on the main thread (its per-sample readers + REF/dust
    // fetchers are thread-local); the callers share `&caller` (it is `Sync` —
    // stateless `call_chunk`), and the writer owns the (`Send`) `CohortVcfWriter`
    // and reorders the out-of-order `CalledChunk`s by `chunk_order`, so the
    // emitted VCF is byte-identical regardless of worker count.
    let requested_threads = args
        .threads
        .unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
        })
        .max(1);
    // H1: size the producer (rayon) and caller (crossbeam) pools from one
    // budget so they don't each claim `T` cores and oversubscribe ~2×. The
    // producer's `par_iter_mut` runs inside `producer_pool.install(...)` below,
    // so it uses this pool — not the global one — and the bench (which drives
    // this entry directly) measures the same sizing the CLI gets.
    let (producer_threads, n_workers) = resolve_thread_split(requested_threads);
    let producer_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(producer_threads)
        .build()?;
    let cap = (QUEUE_DEPTH_PER_WORKER * n_workers).max(1);
    // Surface the resolved concurrency / chunk-sizing defaults so an operator
    // can recover what a running instance actually chose (the values default
    // implicitly from `--threads` / the `--target-variants-per-chunk` sentinel).
    eprintln!(
        "var-calling: producer_threads={producer_threads} workers={n_workers} queue_cap={cap} target_variants_per_chunk={target_variants}"
    );
    let (chunk_tx, chunk_rx) = crossbeam_channel::bounded::<RawCohortChunk>(cap);
    let (called_tx, called_rx) = crossbeam_channel::bounded::<CalledChunk>(cap);
    let caller = &caller;

    std::thread::scope(|scope| -> Result<WriterStats, PipelineError> {
        // Writer thread: drain CalledChunks (reorder + filter + write); returns
        // the run-level stats for the summary.
        let writer_handle = scope.spawn(move || -> Result<WriterStats, PipelineError> {
            for called in called_rx {
                writer.handle(called)?;
            }
            Ok(writer.finish()?)
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
        // Main holds no extra channel handles, so the channels close once the
        // producer / callers finish.
        drop(chunk_rx);
        drop(called_tx);

        // Producer on the main thread. Its `par_iter_mut` decode/compaction
        // runs on `producer_pool` (H1: sized independently of the callers)
        // because the whole loop executes inside `producer_pool.install(...)`.
        let mut ref_fetcher: Option<(u32, StreamingChromRefFetcher)> = None;
        let mut dust_fetcher: Option<(u32, ManualEvictChromRefFetcher)> = None;
        let produce = producer_pool.install(|| -> Result<(), PipelineError> {
            for chrom_id in 0..n_chromosomes {
                let mut intervals = producer.covered_intervals(chrom_id, max_group_span);
                if let Some(rs) = &region_set {
                    intervals = restrict_intervals_to_regions(&intervals, rs.regions_for(chrom_id));
                }
                for interval in intervals {
                    let mask = dust_mask_for(
                        &mut dust_fetcher,
                        chrom_id,
                        &interval,
                        &args.reference,
                        &chrom_names,
                        &chrom_lengths,
                        args.no_complexity_filter,
                        cohort.complexity_window,
                        cohort.complexity_threshold,
                    )?;
                    producer.begin_interval(chrom_id, interval, mask);
                    loop {
                        let chunk = {
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
                            producer.produce_chunk(&mut fetch)?
                        };
                        let Some(chunk) = chunk else { break };
                        if chunk_tx.send(chunk).is_err() {
                            // All callers gone (errored); stop. The real error
                            // surfaces from a caller's join below.
                            return Ok(());
                        }
                    }
                }
            }
            Ok(())
        });
        // Signal the callers there are no more chunks.
        drop(chunk_tx);

        // First error wins across producer, callers, and writer; otherwise
        // return the writer's run stats.
        //
        // The `.join().expect(...)` re-raises a worker-thread panic on this
        // (the scope-owning) thread. A caller/writer panic is treated as
        // fatal-by-design: it is a bug, not a recoverable run error. It does
        // not deadlock — a panicking thread unwinds and drops its channel-handle
        // clones, so the surviving threads observe channel close, drain, and
        // join (`scope` then completes). `join()` only returns `Err` on panic.
        let mut first_err = produce.err();
        for handle in caller_handles {
            if let Err(e) = handle.join().expect("caller thread panicked") {
                first_err.get_or_insert(e);
            }
        }
        let writer_result = writer_handle.join().expect("writer thread panicked");
        match (first_err, writer_result) {
            (Some(e), _) | (None, Err(e)) => Err(e),
            (None, Ok(stats)) => Ok(stats),
        }
    })
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
