//! Chunk-loop driver for the cohort var-calling pipeline.
//!
//! Replaces the per-chromosome `rayon::par_iter` from
//! [`run_var_calling`](crate::pop_var_caller::var_calling::run_var_calling)
//! with a sequential per-chromosome outer loop + a chunk-loop inner
//! loop that uses the new columnar machinery in this module
//! (loader, pre-pass, partition, worker). Phase A is sequential
//! within and across chromosomes; Phase B extends to parallel
//! windows.
//!
//! **Why "drive" not "iterate"** — the driver owns the persistent
//! scratch reused across every block and chromosome: the
//! `StreamingBlockLoader` (its pending columns + fold scratch), the
//! `PartitionScratch`, the DUST-mask buffer, and the recycled
//! `ReadyBlock` buffers. Streaming compaction keeps the resident set at
//! ~one block per sample + a group-sized straddle, not the whole span ×
//! N — the bulk of the memory-budget win the rewrite promises.
//!
//! **Downstream filters** (`is_variant_call`, `qual_phred`,
//! `min_alt_obs_per_sample`, MAPQ-diff t-test) are applied
//! post-EM in the driver, mirroring the streaming pipeline's drop
//! semantics in
//! [`drive_cohort_pipeline`](crate::var_calling::from_bam::pipeline::drive_cohort_pipeline).
//! The same set of records is dropped from the final VCF, with
//! the same per-category counters incremented — that's the
//! byte-identity contract for Phase A.

use std::collections::{BTreeMap, HashMap, VecDeque};
use std::fs::File;
use std::io::{self, BufReader, Read, Seek};
use std::num::NonZeroU32;
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::sync::{Arc, Condvar, Mutex};
use std::thread::JoinHandle;
use thiserror::Error;

use crate::fasta::fetcher::ChromRefFetchError;
use crate::fasta::{ManualEvictChromRefFetcher, StreamingChromRefFetcher};
use crate::psp::header::ParsedChromosome;
use crate::psp::{PspReadError, PspReader};
use crate::var_calling::dust_filter::{DustFilterConfig, MIN_DUST_HALO, sdust_mask_for_span};
use crate::var_calling::from_psp::column_span_reader::ColumnSpanReader;
use crate::var_calling::from_psp::columns::MaterialisedChunk;
use crate::var_calling::from_psp::loader::{ChunkLoadError, StreamingBlockLoader};
use crate::var_calling::from_psp::partition::{
    PartitionError, PartitionScratch, WindowPartition, partition_window,
};
use crate::var_calling::from_psp::worker::{
    WindowRunStats, WorkerScratch, prefetch_window_ref_bytes, run_window,
};
use crate::var_calling::per_group_merger::{PerGroupMergerConfig, PerGroupMergerError};
use crate::var_calling::posterior_engine::{
    PosteriorEngineConfig, PosteriorEngineError, PosteriorRecord,
};
use crate::var_calling::variant_grouping::GrouperConfig;
use crate::vcf::{CohortMetadata, CohortVcfWriter, VcfWriteError, WriterConfig, tmp_path_for};

/// Default nominal genomic span per chunk, in bases. Picked to fit
/// ~1–2 PSP blocks per sample at the median record density on
/// tomato; the chunk-loop overhead is tiny so the exact value is
/// not critical. Exposed as a knob in [`ChunkDriverParams`].
pub const DEFAULT_CHUNK_GENOMIC_SPAN: u32 = 100_000;

/// I/O buffer capacity for per-sample PSP `BufReader`s. Mirrors the
/// `pop_var_caller::common::DEFAULT_BUFFERED_IO_CAPACITY` constant
/// — kept local so the driver doesn't reach into the CLI module.
const DRIVER_PSP_BUFFER_BYTES: usize = 64 * 1024;

/// Threshold from [`PerPositionMerger::record_fails_mapq_diff_t`]
/// semantics: a side of the Welch's-t test needs at least this many
/// supporting reads. Duplicated here so the chunk driver doesn't
/// reach across modules.
const MAPQ_FILTER_MIN_READS_PER_SIDE: u64 = 3;

/// Mi19: chunk-loop sizing knobs. Bundled separately from
/// `ChunkDriverParams` so the chunk loader / pre-pass / worker-
/// dispatch axis is named and the three fields can be passed around
/// together without dragging along the per-stage configs and
/// downstream filters.
#[non_exhaustive]
#[derive(Debug, Clone, Copy)]
pub struct ChunkSizingParams {
    /// Nominal BP span of the loader's first pull attempt per chunk.
    /// Acts as a starting size; the loader grows up to
    /// `chunk_genomic_span × MAX_CHUNK_SPAN_GROWTH` when
    /// `target_variants_per_chunk` is not satisfied or the pre-pass
    /// can't find a safe boundary.
    pub chunk_genomic_span: u32,
    /// M9: soft lower bound on the post-filter variant count per
    /// chunk. `None` disables the variant-bounded extension (each
    /// chunk is one pull attempt of `chunk_genomic_span` BP).
    /// `Some(n)` lets the loader grow each chunk's span adaptively
    /// until the kept-position count crosses `n` — decoupling worker
    /// workload from PSP block size + variant density. Previously a
    /// bare `u32` with `0` as the disabling sentinel; the
    /// `Option<NonZero…>` form makes the disabled/enabled split
    /// visible in the type. See
    /// [the Phase B prereq plan](https://github.com/JoseBlanca/join_per_sample_vcfs/blob/main/doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_b1_variant_bounded_chunks.md).
    pub target_variants_per_chunk: Option<NonZeroU32>,
}

/// Mi19: downstream per-record filters. `min_alt_obs_per_sample` is
/// applied **pre-EM** at the merger→EM boundary (in
/// `build_posterior_record_columnar`); the rest (`min_qual_phred`,
/// MAPQ-diff t-test) are applied **post-EM** by `emit_or_drop` before
/// each record reaches the VCF writer. Grouped separately from
/// `ChunkDriverParams` so the per-record filter axis is named.
#[non_exhaustive]
#[derive(Debug, Clone, Copy)]
pub struct DownstreamFilterParams {
    /// Minimum per-allele alt-read support — at the per-group merger
    /// → EM boundary, a record is dropped if no ALT allele has
    /// `max(num_obs across samples) >= min_alt_obs_per_sample`.
    pub min_alt_obs_per_sample: u32,
    /// Minimum site-level `QUAL` (phred) required to emit a record.
    /// Records with `qual_phred < min_qual_phred` are dropped.
    pub min_qual_phred: f64,
    /// Skip the Welch's-t MAPQ-difference drop entirely.
    pub no_mapq_diff_filter: bool,
    /// Drop records whose worst per-ALT Welch's t-statistic for the
    /// MAPQ(alt) − MAPQ(ref) distribution falls below this threshold.
    pub min_mapq_diff_t: f32,
}

/// Per-driver tuning bundle. Carries every per-stage config the
/// chunk loop needs.
///
/// Mi19: the per-axis knobs live in dedicated sub-structs
/// ([`ChunkSizingParams`], [`DownstreamFilterParams`]) so the top-
/// level field set names the axes rather than mixing 12 flat
/// `u32` / `bool` / `f64` fields. Per-stage configs (`dust_cfg`,
/// `grouper_cfg`, `per_group_cfg`, `posterior_cfg`) stay at the top
/// level — each is a single subsystem's config bundle and already
/// carries its own internal grouping.
#[non_exhaustive]
pub struct ChunkDriverParams {
    pub no_complexity_filter: bool,
    pub dust_cfg: DustFilterConfig,
    pub grouper_cfg: GrouperConfig,
    pub per_group_cfg: PerGroupMergerConfig,
    pub posterior_cfg: PosteriorEngineConfig,
    /// Chunk-loop sizing knobs — see [`ChunkSizingParams`].
    pub sizing: ChunkSizingParams,
    /// Post-EM downstream filters — see [`DownstreamFilterParams`].
    pub downstream: DownstreamFilterParams,
}

/// Per-driver counters; the shape mirrors
/// [`CohortDriveStats`](crate::var_calling::from_bam::pipeline::CohortDriveStats)
/// so the CLI summary code can consume either driver's output
/// without a branch.
// Mi1: `#[non_exhaustive]` — driver-level counter struct; future
// per-stage counts can be added without breaking out-of-crate
// consumers.
#[non_exhaustive]
#[derive(Debug, Default, Clone)]
pub struct ChunkDriverStats {
    pub records_written: u64,
    pub records_unconverged: u64,
    pub records_dropped_hom_ref: u64,
    pub records_dropped_low_qual: u64,
    pub records_dropped_low_alt_obs: u64,
    pub records_dropped_low_mapq_diff_t: u64,
    pub lh_cap_groups_skipped: u64,
    pub lh_cap_alleles_in_skipped: u64,
    /// M23: groups skipped because post-unification the allele set
    /// was REF-only (every ALT was absorbed into the OTHER pool by
    /// the max-alleles cap, or none survived the per-position
    /// dedup). Mirrors the row-shape merger's REF-only no-op skip;
    /// surfaced separately from `lh_cap_groups_skipped` because the
    /// trigger is different (post-unify vs cap-trip).
    pub groups_skipped_post_unify_ref_only: u64,
    /// Total chunks loaded — includes both terminal loads and
    /// (in the future) retry attempts. Useful for capacity planning
    /// alongside `lh_cap_*` counters.
    pub chunks_loaded: u64,
    /// Sum of `ChunkLoadStats.variant_count` across every chunk
    /// load. Divided by `chunks_loaded` it gives the average
    /// kept-position count per chunk — a quick sanity check that
    /// `target_variants_per_chunk` is doing what the operator expects.
    pub chunk_variants_total: u64,
}

/// Errors surfaced by the chunk-loop driver.
///
/// M4 / M21 — variants name the *operation* being attempted, not the
/// subsystem. Inner errors are carried via `#[source]` so
/// `source()` produces the canonical cause-chain rendering; the
/// `#[error]` message describes what the driver was doing when the
/// failure occurred. Drop categories:
///
/// - **PSP I/O at construction time**: `OpenPspFile` (`io::Error` on
///   `File::open`), `OpenPspReader` (`PspReadError` on header parse).
/// - **Cohort VCF writer**: `OpenVcfWriter` (constructor),
///   `FinishVcfWriter` (terminal rename + flush), `WriteVcf`
///   (per-record; carries the failing record's `(chrom_id, start, end)`
///   per M21).
/// - **Reference / DUST**: `OpenRefFetcher` (per-chrom fetcher
///   construction), `ComputeDustMask` (chromosome-wide mask scan),
///   `PrefetchRefBytes` (per-window REF pre-fetch).
/// - **Per-block pipeline**: `StreamRead`, `Partition` (each carries
///   `chrom_id`), `RunWindow` (per-block math via `run_window`).
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum ChunkDriverError {
    #[error("failed to open PSP file for sample {sample_idx} at {path}")]
    OpenPspFile {
        sample_idx: usize,
        path: PathBuf,
        #[source]
        source: io::Error,
    },
    #[error("failed to read PSP header for sample {sample_idx} at {path}")]
    OpenPspReader {
        sample_idx: usize,
        path: PathBuf,
        #[source]
        source: PspReadError,
    },
    #[error("failed to open cohort VCF writer at {output}")]
    OpenVcfWriter {
        output: PathBuf,
        #[source]
        source: VcfWriteError,
    },
    #[error("failed to finalise cohort VCF writer at {output}")]
    FinishVcfWriter {
        output: PathBuf,
        #[source]
        source: VcfWriteError,
    },
    /// M21: per-record write failure carries the failing record's
    /// `(chrom_id, start, end)` so operators triaging a mid-run I/O
    /// failure (ENOSPC, broken pipe, etc.) can locate the source
    /// position.
    #[error("failed to write VCF record at {chrom_id}:{start}..{end}")]
    WriteVcf {
        chrom_id: u32,
        start: u32,
        end: u32,
        #[source]
        source: VcfWriteError,
    },
    #[error("failed to open reference fetcher for contig {contig}")]
    OpenRefFetcher {
        contig: String,
        #[source]
        source: ChromRefFetchError,
    },
    #[error("failed to compute DUST mask for contig {contig}")]
    ComputeDustMask {
        contig: String,
        #[source]
        source: ChromRefFetchError,
    },
    #[error("DUST mask scan failed for contig {contig}")]
    DustMaskIo {
        contig: String,
        #[source]
        source: io::Error,
    },
    #[error("failed to pre-fetch REF bytes for window on chromosome {chrom_id}")]
    PrefetchRefBytes {
        chrom_id: u32,
        #[source]
        source: ChromRefFetchError,
    },
    #[error("failed to load chunk on chromosome {chrom_id}")]
    LoadChunk {
        chrom_id: u32,
        #[source]
        source: ChunkLoadError<PspReadError>,
    },
    #[error("failed to read columns on chromosome {chrom_id}")]
    StreamRead {
        chrom_id: u32,
        #[source]
        source: PspReadError,
    },
    /// The DUST-ahead thread ended before serving the covered interval
    /// the producer was entering — the producer-vs-DUST interval
    /// sequences fell out of lockstep, or the thread panicked. Surfaced
    /// rather than silently DUSTing nothing (which would change output).
    #[error("DUST-ahead thread ended before serving an interval on chromosome {chrom_id}")]
    DustAheadGone { chrom_id: u32 },
    #[error("failed to partition window on chromosome {chrom_id}")]
    Partition {
        chrom_id: u32,
        #[source]
        source: PartitionError,
    },
    /// Per-window math (`run_window`) failed: layer 1 / 2 / 3 / EM
    /// surfaced an error. Does not currently carry the failing
    /// chromosome / window — `PosteriorEngineError`'s inner variants
    /// carry their own `RecordLocus` for diagnostic.
    #[error("per-window math failed")]
    RunWindow(#[source] PosteriorEngineError),
    /// Carried for `From<PerGroupMergerError>` compatibility with
    /// callers that still surface that type directly (e.g. the
    /// `var-calling-from-bam` shim before its own migration). Wraps
    /// the error with no extra context — provenance lives in the
    /// inner variants.
    #[error("per-group merger error")]
    PerGroupMerger(#[source] PerGroupMergerError),
}

// `From` impls — kept narrow. Each `From` has exactly one production
// origin in this module, so provenance is preserved by the variant
// itself. `?` sites use explicit `.map_err(|source| …)` to attach the
// chrom/window/path context that the variant fields require.

impl From<PosteriorEngineError> for ChunkDriverError {
    fn from(source: PosteriorEngineError) -> Self {
        Self::RunWindow(source)
    }
}

impl From<PerGroupMergerError> for ChunkDriverError {
    fn from(source: PerGroupMergerError) -> Self {
        Self::PerGroupMerger(source)
    }
}

/// Top-level entry: drive the chunk loop across every chromosome
/// in `chromosomes`, emitting the cohort VCF at `output`.
///
/// Owns the persistent scratch + the one `CohortVcfWriter` shared
/// across chromosomes. Per-sample PSP readers are opened once and
/// re-used across chromosomes via `PspReader::region_records`,
/// which seeks via the block index per call.
#[allow(clippy::too_many_arguments)]
pub fn drive_cohort_chunked(
    psp_paths: &[PathBuf],
    sample_names: Vec<String>,
    chromosomes: Vec<ParsedChromosome>,
    fasta_path: &Path,
    output: &Path,
    metadata: CohortMetadata,
    writer_cfg: WriterConfig,
    params: ChunkDriverParams,
) -> Result<ChunkDriverStats, ChunkDriverError> {
    let n_samples = sample_names.len();
    debug_assert_eq!(psp_paths.len(), n_samples);

    // Mi20: surface every effective `ChunkDriverParams` value on
    // stderr at startup so an operator looking at a run log can
    // recover "what config did this run actually use?" without
    // having to read the CLI invocation. Pairs with the variant-
    // filter and chunk-loop counters that already land in
    // `print_run_summary` at the end of the run.
    eprintln!(
        "var-calling: chunk_genomic_span={} target_variants_per_chunk={} \
         min_qual_phred={} min_alt_obs_per_sample={} \
         no_mapq_diff_filter={} min_mapq_diff_t={} no_complexity_filter={}",
        params.sizing.chunk_genomic_span,
        params
            .sizing
            .target_variants_per_chunk
            .map_or(0, NonZeroU32::get),
        params.downstream.min_qual_phred,
        params.downstream.min_alt_obs_per_sample,
        params.downstream.no_mapq_diff_filter,
        params.downstream.min_mapq_diff_t,
        params.no_complexity_filter,
    );

    // Open one PSP reader per sample. Reused across chromosomes;
    // `region_records` re-seeks via the block index per call.
    let mut psp_readers: Vec<PspReader<BufReader<File>>> = Vec::with_capacity(psp_paths.len());
    for (sample_idx, path) in psp_paths.iter().enumerate() {
        let file = File::open(path).map_err(|source| ChunkDriverError::OpenPspFile {
            sample_idx,
            path: path.clone(),
            source,
        })?;
        let buf = BufReader::with_capacity(DRIVER_PSP_BUFFER_BYTES, file);
        let reader = PspReader::new(buf).map_err(|source| ChunkDriverError::OpenPspReader {
            sample_idx,
            path: path.clone(),
            source,
        })?;
        psp_readers.push(reader);
    }

    let mut writer = CohortVcfWriter::new(metadata, writer_cfg).map_err(|source| {
        ChunkDriverError::OpenVcfWriter {
            output: output.to_path_buf(),
            source,
        }
    })?;

    let mut stats = ChunkDriverStats::default();

    // The block generator owns the PSP readers, the reusable scratch,
    // and the reserve; it yields prepared blocks. The consumer below
    // just pulls a block and runs the per-group math — it never sees
    // how the block was produced (data-shaping vs math decoupled).
    let per_group_cfg = params.per_group_cfg;
    let posterior_cfg = &params.posterior_cfg;
    let mut producer = BlockIterator::new(psp_readers, &chromosomes, fasta_path, &params);

    // Consume the blocks: per-group EM/posterior math + downstream
    // filtering + emit. The math is the dominant cost at realistic N
    // (consume grows ~N^1.8, overtaking produce at N≈24), so it is
    // parallelised across blocks when the pool has > 1 thread — a single
    // producer thread feeds a bounded queue, `n_workers` worker threads
    // map `process_block` in parallel, and a collector re-imposes genomic
    // order before emit (byte-identical to the serial path). With a
    // 1-thread pool the serial path runs. Both paths fold the producer's
    // chunk counters into `stats`.
    let n_workers = rayon::current_num_threads();
    let driver_result: Result<(), ChunkDriverError> = if n_workers <= 1 {
        drive_blocks_serial(
            &mut producer,
            per_group_cfg,
            posterior_cfg,
            &params,
            &mut writer,
            &mut stats,
        )
    } else {
        drive_blocks_parallel(
            producer,
            n_workers,
            per_group_cfg,
            posterior_cfg,
            &params,
            &mut writer,
            &mut stats,
        )
    };

    match driver_result {
        Ok(()) => {
            writer
                .finish()
                .map_err(|source| ChunkDriverError::FinishVcfWriter {
                    output: output.to_path_buf(),
                    source,
                })?;
            Ok(stats)
        }
        Err(e) => {
            if let Err(remove_err) = writer.abort() {
                eprintln!(
                    "var-calling: failed to remove tmp VCF {} during error cleanup: {}",
                    tmp_path_for(output).display(),
                    remove_err,
                );
            }
            Err(e)
        }
    }
}

/// Roll one block's [`WindowRunStats`] counters into the driver's
/// running totals. Additive, so the result is independent of the order
/// blocks finish — the parallel collector and the serial loop produce
/// identical totals.
fn roll_window_stats(stats: &mut ChunkDriverStats, block: &WindowRunStats) {
    stats.lh_cap_groups_skipped += block.lh_cap_groups_skipped;
    stats.lh_cap_alleles_in_skipped += block.lh_cap_alleles_in_skipped;
    stats.groups_skipped_post_unify_ref_only += block.groups_skipped_post_unify_ref_only;
    stats.records_dropped_low_alt_obs += block.records_dropped_low_alt_obs;
}

/// Serial consume: pull each block in genomic order, run the math, emit.
/// The reference path for byte-identity, and the path taken when the
/// rayon pool has a single thread.
fn drive_blocks_serial<R: Read + Seek + Send>(
    producer: &mut BlockIterator<R>,
    per_group_cfg: PerGroupMergerConfig,
    posterior_cfg: &PosteriorEngineConfig,
    params: &ChunkDriverParams,
    writer: &mut CohortVcfWriter,
    stats: &mut ChunkDriverStats,
) -> Result<(), ChunkDriverError> {
    let mut worker = WorkerScratch::new();
    let result: Result<(), ChunkDriverError> = (|| {
        while let Some(res) = producer.next_block() {
            let block = res?;
            process_block(
                &block,
                per_group_cfg,
                posterior_cfg,
                params.downstream.min_alt_obs_per_sample,
                &mut worker,
            )
            .map_err(ChunkDriverError::RunWindow)?;
            roll_window_stats(stats, &worker.stats);
            worker.stats.clear();
            for record in worker.posterior_records.drain(..) {
                emit_or_drop(record, params, writer, stats)?;
            }
            producer.recycle(block);
        }
        Ok(())
    })();
    stats.chunks_loaded += producer.chunks_loaded;
    stats.chunk_variants_total += producer.chunk_variants_total;
    result
}

/// One finished block's math output, handed worker → collector. The
/// records are already in genomic order *within* the block; the
/// collector emits blocks in `seq_idx` order, so the global emit order
/// equals the serial path's — that's the byte-identity argument.
struct BlockResult {
    seq_idx: u64,
    records: Vec<PosteriorRecord>,
    stats: WindowRunStats,
}

/// Parallel consume. A single producer thread feeds prepared blocks
/// into a bounded [`BlockQueue`] (back-pressure bounds resident memory
/// to ~`cap` + `n_workers` blocks); `n_workers` worker threads pop
/// blocks, run [`process_block`] with their own persistent
/// [`WorkerScratch`], and forward owned [`BlockResult`]s while recycling
/// the spent block buffers back to the producer. The collector (this
/// thread) re-imposes `seq_idx` order and emits — byte-identical to
/// [`drive_blocks_serial`].
///
/// Workers are dedicated OS threads, **not** rayon tasks: a worker
/// blocked on `queue.pop()` must not occupy a rayon pool thread, because
/// the producer's per-sample decode uses that pool — parking workers on
/// it could starve or deadlock the decode. The producer thread is the
/// only one that touches the rayon pool (via the loader's `par_iter`).
///
/// On the success path no errors occur and output is identical to
/// serial. On error (a worker's math, a produced block, or a VCF write)
/// the first error wins, the queue is closed to wind every thread down,
/// and that error is returned (the run aborts regardless of which block
/// failed first — errors are pathological-input aborts, not part of the
/// byte-identity contract).
fn drive_blocks_parallel<R: Read + Seek + Send>(
    mut producer: BlockIterator<R>,
    n_workers: usize,
    per_group_cfg: PerGroupMergerConfig,
    posterior_cfg: &PosteriorEngineConfig,
    params: &ChunkDriverParams,
    writer: &mut CohortVcfWriter,
    stats: &mut ChunkDriverStats,
) -> Result<(), ChunkDriverError> {
    // Bounded look-ahead: a few blocks per worker keeps every worker fed
    // without letting the producer run far ahead (memory is ~this many
    // blocks × N samples). 2× workers is enough to hide per-block
    // production jitter behind the math.
    let queue = BlockQueue::new(n_workers * 2);
    let (result_tx, result_rx) = mpsc::channel::<Result<BlockResult, ChunkDriverError>>();
    let (recycle_tx, recycle_rx) = mpsc::channel::<ReadyBlock>();

    std::thread::scope(|scope| {
        // ── Producer: generate blocks in genomic order, recycling spent
        //    buffers drained from `recycle_rx`. Returns the chunk
        //    counters (or the first production error) via join. ──
        let queue_ref = &queue;
        let producer_handle = scope.spawn(move || -> Result<(u64, u64), ChunkDriverError> {
            loop {
                while let Ok(spent) = recycle_rx.try_recv() {
                    producer.recycle(spent);
                }
                match producer.next_block() {
                    None => break,
                    Some(Ok(block)) => {
                        if queue_ref.push(block).is_err() {
                            // Queue closed (shutdown after an error).
                            break;
                        }
                    }
                    Some(Err(e)) => {
                        queue_ref.close();
                        return Err(e);
                    }
                }
            }
            queue_ref.close();
            Ok((producer.chunks_loaded, producer.chunk_variants_total))
        });

        // ── Workers: map `process_block` over blocks with private
        //    scratch; forward results, recycle blocks. ──
        // `u32`, captured by Copy into each worker closure.
        let min_alt_obs_per_sample = params.downstream.min_alt_obs_per_sample;
        for _ in 0..n_workers {
            let queue_ref = &queue;
            let result_tx = result_tx.clone();
            let recycle_tx = recycle_tx.clone();
            scope.spawn(move || {
                let mut worker = WorkerScratch::new();
                while let Some(block) = queue_ref.pop() {
                    match process_block(
                        &block,
                        per_group_cfg,
                        posterior_cfg,
                        min_alt_obs_per_sample,
                        &mut worker,
                    ) {
                        Ok(()) => {
                            let out = BlockResult {
                                seq_idx: block.seq_idx,
                                records: std::mem::take(&mut worker.posterior_records),
                                stats: std::mem::take(&mut worker.stats),
                            };
                            // Recycle first so the producer can reuse the
                            // buffer while we ship the result.
                            let _ = recycle_tx.send(block);
                            if result_tx.send(Ok(out)).is_err() {
                                break;
                            }
                        }
                        Err(e) => {
                            let _ = result_tx.send(Err(ChunkDriverError::RunWindow(e)));
                            break;
                        }
                    }
                }
            });
        }
        // Drop the originals so `result_rx` closes once every worker's
        // clone is dropped (and `recycle_rx` once producer + workers end).
        drop(result_tx);
        drop(recycle_tx);

        // ── Collector (this thread): emit in `seq_idx` order. ──
        let mut next_expected: u64 = 0;
        let mut pending: BTreeMap<u64, BlockResult> = BTreeMap::new();
        let mut first_err: Option<ChunkDriverError> = None;
        for msg in result_rx {
            match msg {
                Err(e) => {
                    if first_err.is_none() {
                        first_err = Some(e);
                        queue.close(); // wind down producer + workers
                    }
                }
                Ok(out) => {
                    if first_err.is_some() {
                        continue; // draining; nothing more is emitted
                    }
                    pending.insert(out.seq_idx, out);
                    while let Some(out) = pending.remove(&next_expected) {
                        roll_window_stats(stats, &out.stats);
                        for record in out.records {
                            if let Err(e) = emit_or_drop(record, params, writer, stats) {
                                first_err = Some(e);
                                queue.close();
                                break;
                            }
                        }
                        if first_err.is_some() {
                            break;
                        }
                        next_expected += 1;
                    }
                }
            }
        }

        // All workers have ended (result_rx closed). Join the producer
        // for its counters / production error.
        let producer_result = producer_handle.join().expect("producer thread panicked");

        if let Some(e) = first_err {
            return Err(e);
        }
        match producer_result {
            Ok((chunks_loaded, chunk_variants_total)) => {
                stats.chunks_loaded += chunks_loaded;
                stats.chunk_variants_total += chunk_variants_total;
                Ok(())
            }
            Err(e) => Err(e),
        }
    })
}

/// A bounded, blocking multi-producer/multi-consumer queue of prepared
/// blocks. The single producer `push`es (blocking while `cap` blocks are
/// already in flight — the back-pressure that bounds resident memory);
/// the worker pool `pop`s (blocking while empty). `close` wakes everyone
/// for shutdown: after it, `push` returns `Err` and `pop` drains then
/// returns `None`.
struct BlockQueue {
    inner: Mutex<BlockQueueInner>,
    not_empty: Condvar,
    not_full: Condvar,
}

struct BlockQueueInner {
    q: VecDeque<ReadyBlock>,
    cap: usize,
    closed: bool,
}

impl BlockQueue {
    fn new(cap: usize) -> Self {
        Self {
            inner: Mutex::new(BlockQueueInner {
                q: VecDeque::with_capacity(cap),
                cap: cap.max(1),
                closed: false,
            }),
            not_empty: Condvar::new(),
            not_full: Condvar::new(),
        }
    }

    /// Enqueue, blocking while full. `Err(block)` if closed (shutdown).
    fn push(&self, block: ReadyBlock) -> Result<(), ReadyBlock> {
        let mut g = self.inner.lock().expect("BlockQueue mutex poisoned");
        while g.q.len() >= g.cap && !g.closed {
            g = self.not_full.wait(g).expect("BlockQueue mutex poisoned");
        }
        if g.closed {
            return Err(block);
        }
        g.q.push_back(block);
        drop(g);
        self.not_empty.notify_one();
        Ok(())
    }

    /// Dequeue, blocking while empty. `None` once closed *and* drained.
    fn pop(&self) -> Option<ReadyBlock> {
        let mut g = self.inner.lock().expect("BlockQueue mutex poisoned");
        loop {
            if let Some(block) = g.q.pop_front() {
                drop(g);
                self.not_full.notify_one();
                return Some(block);
            }
            if g.closed {
                return None;
            }
            g = self.not_empty.wait(g).expect("BlockQueue mutex poisoned");
        }
    }

    /// Signal shutdown: unblock all waiters. Idempotent.
    fn close(&self) {
        let mut g = self.inner.lock().expect("BlockQueue mutex poisoned");
        g.closed = true;
        drop(g);
        self.not_empty.notify_all();
        self.not_full.notify_all();
    }
}

/// Apply the streaming pipeline's downstream filters to `record`
/// in the same order; either write it to `writer` or increment the
/// matching `stats` counter and drop it.
fn emit_or_drop(
    record: PosteriorRecord,
    params: &ChunkDriverParams,
    writer: &mut CohortVcfWriter,
    stats: &mut ChunkDriverStats,
) -> Result<(), ChunkDriverError> {
    // NB: the `min_alt_obs_per_sample` filter is NOT here — it runs
    // pre-EM in `build_posterior_record_columnar` (the merger→EM
    // boundary), matching `drive_cohort_pipeline` so the decision is
    // made on the pre-prune allele set. See
    // `WindowRunStats::records_dropped_low_alt_obs`.
    if !record.is_variant_call() {
        stats.records_dropped_hom_ref += 1;
        return Ok(());
    }
    if record.qual_phred < params.downstream.min_qual_phred {
        stats.records_dropped_low_qual += 1;
        return Ok(());
    }
    if !params.downstream.no_mapq_diff_filter
        && record_fails_mapq_diff_t(&record, params.downstream.min_mapq_diff_t)
    {
        stats.records_dropped_low_mapq_diff_t += 1;
        return Ok(());
    }
    if !record.diagnostics.converged {
        stats.records_unconverged += 1;
    }
    // M21: attach the failing record's locus to the typed write
    // error so operators triaging mid-run I/O failures can locate the
    // source position from the error chain.
    writer
        .write_record(&record)
        .map_err(|source| ChunkDriverError::WriteVcf {
            chrom_id: record.locus.chrom_id,
            start: record.locus.start,
            end: record.locus.end,
            source,
        })?;
    stats.records_written += 1;
    Ok(())
}

/// Integration-test alias for [`record_fails_mapq_diff_t`]. The helper
/// itself is private to the driver (single call site in `emit_or_drop`);
/// the alias lets `cohort_vcf_writer_integration` exercise the decision
/// without exposing the production helper as `pub`.
#[doc(hidden)]
pub fn record_fails_mapq_diff_t_for_test(record: &PosteriorRecord, threshold: f32) -> bool {
    record_fails_mapq_diff_t(record, threshold)
}

/// Welch's-t MAPQ-difference filter applied post-EM in `emit_or_drop`.
fn record_fails_mapq_diff_t(record: &PosteriorRecord, threshold: f32) -> bool {
    if !threshold.is_finite() {
        return false;
    }
    let n_alleles = record.alleles.len();
    if n_alleles < 2 {
        return false;
    }
    let ref_moments = pool_allele_mapq(record, 0);
    if ref_moments.n < MAPQ_FILTER_MIN_READS_PER_SIDE {
        return false;
    }
    let mean_ref = ref_moments.sum as f64 / ref_moments.n as f64;
    let var_ref = ((ref_moments.sum_of_squares as f64 - (ref_moments.sum as f64) * mean_ref)
        .max(0.0))
        / ((ref_moments.n - 1) as f64);
    for alt_idx in 1..n_alleles {
        let alt_moments = pool_allele_mapq(record, alt_idx);
        if alt_moments.n < MAPQ_FILTER_MIN_READS_PER_SIDE {
            continue;
        }
        let mean_alt = alt_moments.sum as f64 / alt_moments.n as f64;
        let var_alt = ((alt_moments.sum_of_squares as f64 - (alt_moments.sum as f64) * mean_alt)
            .max(0.0))
            / ((alt_moments.n - 1) as f64);
        let se2 = var_alt / (alt_moments.n as f64) + var_ref / (ref_moments.n as f64);
        if se2 <= 0.0 {
            continue;
        }
        let t = (mean_alt - mean_ref) / se2.sqrt();
        if (t as f32) < threshold {
            return true;
        }
    }
    false
}

/// Mi5: per-allele MAPQ moments pooled across the cohort. Used by
/// [`record_fails_mapq_diff_t`]'s Welch's-t test to compare an ALT's
/// MAPQ distribution against REF's. Wrapping the previously-anonymous
/// `(u64, u64, u128)` tuple in a named struct turns the call site's
/// positional destructure into field access and pins the meaning of
/// each scalar at the type level.
struct PooledMapqMoments {
    /// Total pooled read count across every sample at the allele.
    n: u64,
    /// Sum of per-read MAPQ scores across all those reads.
    sum: u64,
    /// Sum of squares of per-read MAPQ scores (`u128` to absorb the
    /// `mapq² × n` worst case without saturation).
    sum_of_squares: u128,
}

fn pool_allele_mapq(record: &PosteriorRecord, allele_idx: usize) -> PooledMapqMoments {
    let mut moments = PooledMapqMoments {
        n: 0,
        sum: 0,
        sum_of_squares: 0,
    };
    for sample_idx in 0..record.n_samples {
        let stats = &record.scalars_row(sample_idx)[allele_idx];
        moments.n += u64::from(stats.num_obs);
        moments.sum += u64::from(stats.mapq_sum);
        moments.sum_of_squares += stats.mapq_sum_sq as u128;
    }
    moments
}

/// Merge per-block genomic ranges (1-based **inclusive** `(first_pos,
/// last_pos)`) from every sample into sorted, non-overlapping "covered
/// intervals" (1-based **half-open**). Two ranges are merged unless the
/// record gap between them exceeds `max_group_span` — so every gap
/// *between* the returned intervals is `> max_group_span` and is thus a
/// safe chunk boundary (no variant group can reach across it, even after
/// a full chain of joins to the `max_group_span` cap). The chunk loop
/// can therefore process each interval independently and skip the gaps
/// without changing any group — the byte-identity argument for
/// block-index-driven advancement.
///
/// Input order is arbitrary (block indices from N samples interleave);
/// the function sorts. Empty input → empty output.
fn merge_block_ranges(
    ranges: impl IntoIterator<Item = (u32, u32)>,
    max_group_span: u32,
) -> Vec<std::ops::Range<u32>> {
    let mut v: Vec<(u32, u32)> = ranges.into_iter().collect();
    if v.is_empty() {
        return Vec::new();
    }
    v.sort_unstable();
    let mut out: Vec<std::ops::Range<u32>> = Vec::new();
    let (mut cur_start, mut cur_last) = v[0];
    for &(s, last) in &v[1..] {
        // Merge when the next block's first record is within
        // `max_group_span` of the current interval's last record — i.e.
        // the gap is NOT a safe boundary. `cur_last + max_group_span`
        // can't overflow in practice (positions are << u32::MAX) but
        // saturate to be safe.
        if s <= cur_last.saturating_add(max_group_span) {
            cur_last = cur_last.max(last);
        } else {
            out.push(cur_start..cur_last.saturating_add(1));
            cur_start = s;
            cur_last = last;
        }
    }
    out.push(cur_start..cur_last.saturating_add(1));
    out
}

/// Covered intervals for one chromosome: every block's `[first_pos,
/// last_pos]` for `chrom_id`, across all `psp_readers`, merged via
/// [`merge_block_ranges`]. Reads only the in-memory block indices
/// (decoded from the PSP trailer at open) — no file I/O.
fn covered_intervals_for_chrom<W>(
    sources: &[ColumnSpanReader<W>],
    chrom_id: u32,
    max_group_span: u32,
) -> Vec<std::ops::Range<u32>>
where
    W: Read + Seek + Send,
{
    let ranges = sources.iter().flat_map(|s| {
        s.block_index()
            .iter()
            .filter(move |b| b.chrom_id == chrom_id)
            .map(|b| (b.first_pos, b.last_pos))
    });
    merge_block_ranges(ranges, max_group_span)
}

// ─────────────────────────────────────────────────────────────────────
// Stage 3 — DUST-ahead queue.
//
// DUST depends only on the reference sequence + the region, nothing about
// the PSP data or the fold, and the covered intervals are known up front
// from the block indices. Those intervals are independent (no variant
// group reaches across the gap between them), so a **worker pool** DUSTs
// them in parallel and hands the masks to the producer in genomic order;
// the producer slices out the mask for its current block span instead of
// computing it inline. sdust over the whole genome is ~10 s of
// single-threaded work — after Stage 2 the largest serial floor, and the
// wall floor after the Stage 3 single DUST thread — so parallelising it
// across the (independent) intervals is the lever.
// ─────────────────────────────────────────────────────────────────────

/// Look-ahead the [`DustAheadPool`] grants per worker: the pool may keep
/// up to `n_workers × this` covered intervals in flight or
/// completed-but-undelivered ahead of what the producer has consumed. The
/// back-pressure bound — caps in-flight DUST work and resident masks;
/// small because masks are cheap and each in-flight worker already bounds
/// its own reference buffer to ~one sub-span via [`dust_mask_for_interval`].
const DUST_POOL_LOOKAHEAD_PER_WORKER: usize = 2;

/// Sub-span (bases) the DUST-ahead thread DUSTs at a time within a covered
/// interval. `sdust_mask_for_span`'s per-position verdict is byte-
/// identical to a whole-contig scan for any span past the reset barrier,
/// so chunking a long interval into sub-spans and coalescing the masked
/// runs at the boundaries yields the same mask as one whole-interval scan
/// — while keeping the resident reference buffer to ~one sub-span + halo
/// instead of a whole chromosome arm.
const DUST_AHEAD_SUBSPAN: u32 = 1_000_000;

/// Per-chromosome plan shared by the producer and the DUST-ahead thread:
/// the chromosome's id / name / length plus its block-index-driven
/// covered intervals. Both sides iterate the *same* plan list in the same
/// order, so they stay in lockstep on the interval sequence (the producer
/// enters interval N exactly when the thread has served interval N's
/// mask). Built once from the in-memory PSP block indices — no file I/O.
#[derive(Clone)]
struct DustChromPlan {
    chrom_id: u32,
    name: String,
    length: u32,
    /// Block-index-driven covered intervals (1-based half-open),
    /// non-empty. Same value the producer's [`ChromCursor::covered`] holds.
    covered: Vec<std::ops::Range<u32>>,
}

/// One covered interval's precomputed DUST mask, handed DUST-ahead thread
/// → producer over the bounded queue. `interval` is the 1-based half-open
/// covered interval (matched against the interval the producer is
/// entering — a lockstep sanity check); `mask` is the sorted,
/// non-overlapping low-complexity intervals over it, which the producer
/// slices per block.
#[derive(Debug)]
struct IntervalDustMask {
    chrom_id: u32,
    interval: std::ops::Range<u32>,
    mask: Vec<std::ops::Range<u32>>,
}

/// Build the per-chromosome plans (id / name / length + covered
/// intervals) for every chromosome that carries data, in genomic order.
/// Mirrors the skip rules the producer's chromosome advance uses (drop
/// zero-length contigs and contigs no sample covers) so the plan sequence
/// is exactly the interval sequence the producer walks. Reads only the
/// in-memory block indices.
fn build_dust_plans<W>(
    sources: &[ColumnSpanReader<W>],
    chromosomes: &[ParsedChromosome],
    max_group_span: u32,
) -> Vec<DustChromPlan>
where
    W: Read + Seek + Send,
{
    let mut plans = Vec::new();
    for (chrom_idx, chrom) in chromosomes.iter().enumerate() {
        if chrom.length == 0 {
            continue;
        }
        // M25 PANIC-FREE: chrom_id fits in u32 (PSP file-format invariant).
        let chrom_id = u32::try_from(chrom_idx).expect("chrom_id fits in u32");
        let covered =
            covered_intervals_for_chrom(sources, chrom_id, max_group_span.max(MIN_CHUNK_SPAN));
        if covered.is_empty() {
            continue;
        }
        plans.push(DustChromPlan {
            chrom_id,
            name: chrom.name.clone(),
            length: chrom.length,
            covered,
        });
    }
    plans
}

/// Compute the coalesced DUST mask over the covered `interval` (1-based
/// half-open) by DUSTing [`DUST_AHEAD_SUBSPAN`] sub-spans left-to-right
/// and coalescing masked runs across sub-span boundaries. Byte-identical
/// to a single [`sdust_mask_for_span`] over the whole interval (the
/// per-position verdict is identical past the reset barrier), with the
/// reference buffer bounded to ~one sub-span + halo via `evict_before`.
/// Returned intervals are sorted, non-overlapping, 1-based half-open,
/// clipped to `interval`.
fn dust_mask_for_interval(
    fetcher: &mut ManualEvictChromRefFetcher,
    interval: &std::ops::Range<u32>,
    chrom_length: u32,
    window: u32,
    threshold: u32,
    subspan: u32,
) -> io::Result<Vec<std::ops::Range<u32>>> {
    let subspan = subspan.max(1);
    let mut mask: Vec<std::ops::Range<u32>> = Vec::new();
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
                    .map_err(io::Error::other)
            },
        )?;
        // Coalesce a run that the sub-span boundary split: if the first
        // run of this sub-span abuts the last run carried so far, merge
        // them so the result equals a single whole-interval scan.
        for run in sub_mask {
            match mask.last_mut() {
                Some(last) if last.end == run.start => last.end = run.end,
                _ => mask.push(run),
            }
        }
        // Everything left of this sub-span's buffer is done with; the
        // next sub-span only reads from `sub_end - halo` onward, so this
        // keeps the resident reference window to ~one sub-span + halo.
        fetcher.evict_before(buf_start);
        sub_start = sub_end;
    }
    Ok(mask)
}

/// Spawn the DUST-ahead thread: DUST every covered interval in `plans`
/// (genomic order), pushing each `(interval → mask)` into a bounded queue
/// the producer drains as it enters intervals. Off the producer's
/// critical path. Returns the receiver + join handle; the thread exits
/// when it finishes the plans, hits an error (after queueing it), or the
/// receiver is dropped (a `send` to a gone receiver errs).
/// One covered interval queued for the [`DustAheadPool`]: everything a
/// worker needs to DUST it. The job's index in the pool's job list is its
/// delivery sequence — the order the producer enters intervals.
#[derive(Clone)]
struct IntervalJob {
    chrom_id: u32,
    name: String,
    length: u32,
    interval: std::ops::Range<u32>,
}

/// Flatten the per-chromosome plans into one delivery-ordered job list:
/// chromosome-major, interval-order within a chromosome — exactly the
/// order the producer enters intervals (`advance_chrom` walks the plans in
/// order, `interval_idx` 0..n within each). So job index `k` is the k-th
/// interval the producer enters, which is what lets the pool deliver masks
/// in lockstep with the producer.
fn flatten_to_jobs(plans: &[DustChromPlan]) -> Vec<IntervalJob> {
    let mut jobs = Vec::new();
    for plan in plans {
        for interval in &plan.covered {
            jobs.push(IntervalJob {
                chrom_id: plan.chrom_id,
                name: plan.name.clone(),
                length: plan.length,
                interval: interval.clone(),
            });
        }
    }
    jobs
}

/// The per-interval DUST computation a [`DustAheadPool`] worker runs.
/// Behind a trait object so the coordinator can be unit-tested with a stub
/// in place of the real sdust + reference fetch.
type DustComputeFn =
    dyn Fn(&IntervalJob) -> Result<IntervalDustMask, ChunkDriverError> + Send + Sync;

/// Mutable coordination state of a [`DustAheadPool`], behind one mutex.
///
/// Invariants, held under the lock:
/// `next_deliver <= next_dispatch <= next_deliver + lookahead`, and `done`
/// only ever holds keys in `[next_deliver, next_dispatch)`.
struct DustPoolState {
    /// Index of the next job a worker will claim.
    next_dispatch: usize,
    /// Index of the next job the producer will be handed (in order).
    next_deliver: usize,
    /// Completed-but-undelivered results, keyed by job index.
    done: HashMap<usize, Result<IntervalDustMask, ChunkDriverError>>,
    /// Set on shutdown so workers and `recv_next` wind down.
    abort: bool,
}

/// Shared (cross-thread) part of a [`DustAheadPool`].
struct DustPoolShared {
    jobs: Vec<IntervalJob>,
    /// Max `next_dispatch - next_deliver` — how far the pool may run ahead
    /// of what the producer has consumed.
    lookahead: usize,
    state: Mutex<DustPoolState>,
    /// Workers wait here when they are `lookahead` jobs ahead of delivery.
    not_full: Condvar,
    /// `recv_next` waits here for the next-in-order job to complete.
    ready: Condvar,
}

/// A bounded, ordered pool that DUSTs covered intervals in parallel and
/// hands the masks to the producer in genomic order.
///
/// DUST over independent covered intervals is embarrassingly parallel; the
/// only coordination is (a) re-imposing the producer's genomic order on
/// out-of-order completions and (b) bounding how far the pool runs ahead.
/// Workers claim the lowest undispatched job — so they always work on what
/// the producer needs soonest — compute its mask, and stash it; the
/// producer pulls masks in order via [`recv_next`](Self::recv_next). A
/// 1-worker pool reproduces the prior single-thread DUST-ahead.
struct DustAheadPool {
    shared: Arc<DustPoolShared>,
    handles: Vec<JoinHandle<()>>,
}

impl DustAheadPool {
    /// Spawn `n_workers` threads over `jobs`, each running `compute`.
    /// `lookahead` bounds how many jobs may be dispatched ahead of the
    /// delivery frontier. `compute` is shared across workers and must open
    /// its own per-call reference fetcher (no shared mutable state).
    fn spawn(
        jobs: Vec<IntervalJob>,
        lookahead: usize,
        n_workers: usize,
        compute: Arc<DustComputeFn>,
    ) -> Self {
        let shared = Arc::new(DustPoolShared {
            jobs,
            lookahead: lookahead.max(1),
            state: Mutex::new(DustPoolState {
                next_dispatch: 0,
                next_deliver: 0,
                done: HashMap::new(),
                abort: false,
            }),
            not_full: Condvar::new(),
            ready: Condvar::new(),
        });
        let n_workers = n_workers.max(1);
        let mut handles = Vec::with_capacity(n_workers);
        for _ in 0..n_workers {
            let shared = Arc::clone(&shared);
            let compute = Arc::clone(&compute);
            handles.push(std::thread::spawn(move || {
                dust_pool_worker(&shared, &*compute)
            }));
        }
        Self { shared, handles }
    }

    /// Hand the producer the next interval's mask, in genomic order,
    /// blocking until that specific interval completes. `None` once every
    /// interval has been delivered (or on shutdown).
    fn recv_next(&self) -> Option<Result<IntervalDustMask, ChunkDriverError>> {
        let mut state = self
            .shared
            .state
            .lock()
            .expect("DustAheadPool mutex poisoned");
        loop {
            let want = state.next_deliver;
            if let Some(result) = state.done.remove(&want) {
                state.next_deliver = want + 1;
                drop(state);
                // A delivery slot freed: workers may advance one more job.
                self.shared.not_full.notify_all();
                return Some(result);
            }
            if want >= self.shared.jobs.len() || state.abort {
                return None;
            }
            state = self
                .shared
                .ready
                .wait(state)
                .expect("DustAheadPool mutex poisoned");
        }
    }

    /// Signal shutdown and join every worker. Idempotent. A worker mid-DUST
    /// finishes its current interval first, so this blocks at most one
    /// interval's compute.
    fn shutdown(&mut self) {
        {
            let mut state = self
                .shared
                .state
                .lock()
                .expect("DustAheadPool mutex poisoned");
            state.abort = true;
        }
        self.shared.not_full.notify_all();
        self.shared.ready.notify_all();
        for handle in self.handles.drain(..) {
            let _ = handle.join();
        }
    }
}

impl Drop for DustAheadPool {
    fn drop(&mut self) {
        self.shutdown();
    }
}

/// One DUST-pool worker: claim the lowest undispatched job (respecting the
/// look-ahead bound), compute its mask outside the lock, stash it, and
/// repeat until the jobs are exhausted or shutdown is signalled.
fn dust_pool_worker(shared: &DustPoolShared, compute: &DustComputeFn) {
    loop {
        // ── Claim a job (or exit). ──
        let job_idx = {
            let mut state = shared.state.lock().expect("DustAheadPool mutex poisoned");
            loop {
                if state.abort || state.next_dispatch >= shared.jobs.len() {
                    return;
                }
                if state.next_dispatch - state.next_deliver < shared.lookahead {
                    let idx = state.next_dispatch;
                    state.next_dispatch = idx + 1;
                    break idx;
                }
                // Too far ahead of delivery: wait for a slot to free.
                state = shared
                    .not_full
                    .wait(state)
                    .expect("DustAheadPool mutex poisoned");
            }
        };

        let result = compute(&shared.jobs[job_idx]);

        // Stash + wake the producer (it may be waiting for this seq).
        shared
            .state
            .lock()
            .expect("DustAheadPool mutex poisoned")
            .done
            .insert(job_idx, result);
        shared.ready.notify_all();
    }
}

/// Build the production DUST-ahead pool over `plans`: one job per covered
/// interval, a worker per available thread (capped at the job count), each
/// opening its own reference fetcher and running the same
/// [`dust_mask_for_interval`] the single-thread version used. Off the
/// producer's critical path; masks are byte-identical regardless of how
/// many workers race to compute them.
fn spawn_dust_pool(
    fasta_path: &Path,
    plans: &[DustChromPlan],
    window: u32,
    threshold: u32,
) -> DustAheadPool {
    let jobs = flatten_to_jobs(plans);
    let n_workers = rayon::current_num_threads().clamp(1, jobs.len().max(1));
    let lookahead = (n_workers * DUST_POOL_LOOKAHEAD_PER_WORKER).max(1);
    let fasta_path = fasta_path.to_path_buf();
    let compute: Arc<DustComputeFn> =
        Arc::new(move |job: &IntervalJob| {
            let mut fetcher = ManualEvictChromRefFetcher::for_contig(&fasta_path, &job.name)
                .map_err(|source| ChunkDriverError::OpenRefFetcher {
                    contig: job.name.clone(),
                    source,
                })?;
            let mask = dust_mask_for_interval(
                &mut fetcher,
                &job.interval,
                job.length,
                window,
                threshold,
                DUST_AHEAD_SUBSPAN,
            )
            .map_err(|source| ChunkDriverError::DustMaskIo {
                contig: job.name.clone(),
                source,
            })?;
            Ok(IntervalDustMask {
                chrom_id: job.chrom_id,
                interval: job.interval.clone(),
                mask,
            })
        });
    DustAheadPool::spawn(jobs, lookahead, n_workers, compute)
}

/// Slice the precomputed interval mask to the block span
/// `[span_start, span_end)`, clipping each overlapping run to the span and
/// appending it to `out` (sorted, non-overlapping — the shape
/// [`partition_window`](crate::var_calling::from_psp::partition::partition_window)
/// expects). `sweep` is a per-interval forward cursor over
/// `interval_mask`, advanced past runs ending at or before `span_start`; a
/// run straddling the block boundary stays for the next block. Reproduces
/// the inline per-block `sdust_mask_for_span` output exactly because the
/// interval mask carries the same per-position verdict.
fn slice_interval_mask(
    interval_mask: &[std::ops::Range<u32>],
    sweep: &mut usize,
    span_start: u32,
    span_end: u32,
    out: &mut Vec<std::ops::Range<u32>>,
) {
    while *sweep < interval_mask.len() && interval_mask[*sweep].end <= span_start {
        *sweep += 1;
    }
    let mut idx = *sweep;
    while idx < interval_mask.len() && interval_mask[idx].start < span_end {
        let start = interval_mask[idx].start.max(span_start);
        let end = interval_mask[idx].end.min(span_end);
        if start < end {
            out.push(start..end);
        }
        idx += 1;
    }
}

/// A fully-prepared, **owned** block: the materialised columns, their
/// variant-group partition, and the per-group pre-fetched REF bytes —
/// exactly the inputs [`process_block`] needs, and nothing the producer
/// keeps. Owning the payload (rather than borrowing the iterator's
/// reused scratch) is what lets blocks be handed to parallel workers.
///
/// `seq_idx` records the genomic production order so the collector can
/// re-impose it after workers finish out of order. Buffers are recycled
/// through the producer's free-list ([`BlockIterator::recycle`]) so the
/// per-block allocations (columns, CSR arrays, REF byte buffers) are
/// reused, not re-grown, across the run — the scratch-reuse discipline,
/// preserved through the ownership change.
pub struct ReadyBlock {
    /// Monotonic production index (0, 1, 2, …) assigned by the producer
    /// in genomic order. The collector emits records in `seq_idx` order.
    seq_idx: u64,
    /// The loaded columns for `[range.start, safe_end)`.
    chunk: MaterialisedChunk,
    /// The block's single variant-group partition.
    partition: WindowPartition,
    /// One REF byte buffer per group in `partition`
    /// (`len() == partition.n_groups()`).
    pre_fetched_ref_bytes: Vec<Vec<u8>>,
}

impl ReadyBlock {
    /// An empty block pre-sized for `n_samples`; the producer fills it.
    fn with_n_samples(n_samples: usize) -> Self {
        Self {
            seq_idx: 0,
            chunk: MaterialisedChunk::with_n_samples(n_samples),
            partition: WindowPartition::empty(),
            pre_fetched_ref_bytes: Vec::new(),
        }
    }
}

/// Map one prepared block to its final [`PosteriorRecord`]s — the pure
/// data→math seam. Reads only the block + the (per-worker) `scratch`,
/// writes the emitted records into `scratch.posterior_records` and the
/// run counters into `scratch.stats`. No producer state, no shared
/// mutation: this is the function the parallel workers map over blocks.
pub fn process_block(
    block: &ReadyBlock,
    per_group_cfg: PerGroupMergerConfig,
    posterior_cfg: &PosteriorEngineConfig,
    min_alt_obs_per_sample: u32,
    scratch: &mut WorkerScratch,
) -> Result<(), PosteriorEngineError> {
    run_window(
        &block.chunk,
        &block.partition,
        &block.pre_fetched_ref_bytes,
        per_group_cfg,
        posterior_cfg,
        min_alt_obs_per_sample,
        &mut scratch.pipeline_scratch,
        &mut scratch.posterior_records,
        &mut scratch.stats,
    )
}

/// Per-chromosome production cursor owned by [`BlockIterator`]: the
/// reference fetchers, the covered intervals, and the position cursors
/// for the chromosome currently being read. Rebuilt on each chromosome
/// advance.
struct ChromCursor {
    chrom_id: u32,
    chrom_one_past_end: u32,
    /// Forward-only fetcher feeding the per-block REF prefetch.
    ref_fetcher: StreamingChromRefFetcher,
    /// Block-index-driven covered intervals (data-bearing segments).
    covered: Vec<std::ops::Range<u32>>,
    /// Index into `covered` of the interval currently being produced.
    interval_idx: usize,
    /// Genomic start of the next block within the current interval
    /// (advances to each block's `safe_end`; reset to the interval
    /// start when the interval is entered).
    chunk_range_start: u32,
    /// Exclusive end of the current interval (clamped to the chrom).
    interval_one_past_end: u32,
}

/// The chunk/block generator. Owns the per-sample PSP readers, the
/// reusable scratch + the **reserve** (carryover), and yields prepared
/// blocks of ~`desired_num_variants` variable variants. `next_block`
/// reads .psp chunks (decode in parallel across samples), folds to the
/// cohort's variable sites, accumulates until the variant target is
/// met, cuts at a clean group boundary, then slices the precomputed
/// DUST-ahead mask + partitions + prefetches REF bytes for the block.
///
/// This is the *data-shaping* half. `next_block` yields an owned
/// [`ReadyBlock`]; the consumer runs the per-group math via
/// [`process_block`] and emits — never seeing how the block was made.
/// Genomic windows are gone: one block = one unit of work. Spent blocks
/// are returned via [`recycle`](Self::recycle) so their buffers feed the
/// next production without reallocating.
struct BlockIterator<'a, W: Read + Seek + Send> {
    fasta_path: &'a Path,
    params: &'a ChunkDriverParams,
    /// One span-addressable columnar reader per sample; each owns its
    /// `PspReader`. Re-`reset` to the current covered interval.
    sources: Vec<ColumnSpanReader<W>>,
    n_samples: usize,

    /// Streaming fold+compact engine: holds the per-sample pending
    /// (open-group straddle) and the fold scratch, reused across blocks.
    streaming: StreamingBlockLoader,
    /// Reusable per-block partition scratch.
    partition_scratch: PartitionScratch,
    /// Reusable per-block DUST-mask buffer — the current block's slice of
    /// `cur_interval_mask` fed to `partition_window`.
    mask: Vec<std::ops::Range<u32>>,

    // ── DUST-ahead. ──
    /// Per-chromosome plans (covered intervals) built once at construction
    /// and walked by [`advance_chrom`](Self::advance_chrom); the same list
    /// the DUST-ahead pool's jobs are flattened from, keeping the producer
    /// and the pool in lockstep on the interval sequence.
    chrom_plans: Vec<DustChromPlan>,
    next_plan_idx: usize,
    /// Worker pool DUSTing the covered intervals in parallel, delivering
    /// masks in genomic order one per interval entered. `None` when the
    /// complexity filter is disabled (no pool spawned). Joins its workers
    /// when dropped.
    dust_pool: Option<DustAheadPool>,
    /// The current covered interval's DUST mask (whole interval, coalesced)
    /// and a forward sweep cursor into it advanced as blocks progress.
    cur_interval_mask: Vec<std::ops::Range<u32>>,
    mask_sweep: usize,

    // Free-list of recycled block buffers (columns + CSR arrays + REF
    // byte buffers). `next_block` pops one (or makes a fresh one) to
    // fill; the consumer pushes spent blocks back via `recycle`.
    free_blocks: Vec<ReadyBlock>,
    next_seq_idx: u64,

    // Production cursor + counters.
    cur: Option<ChromCursor>,
    chunks_loaded: u64,
    chunk_variants_total: u64,
}

impl<'a, W: Read + Seek + Send> BlockIterator<'a, W> {
    fn new(
        psp_readers: Vec<PspReader<W>>,
        chromosomes: &'a [ParsedChromosome],
        fasta_path: &'a Path,
        params: &'a ChunkDriverParams,
    ) -> Self {
        let n_samples = psp_readers.len();
        let sources: Vec<ColumnSpanReader<W>> = psp_readers
            .into_iter()
            .map(ColumnSpanReader::detached)
            .collect();

        // Plan the covered intervals once (from the in-memory block
        // indices) so the producer and the DUST-ahead pool agree on the
        // interval sequence. DUST off the critical path on a worker pool,
        // unless the complexity filter is disabled.
        let chrom_plans = build_dust_plans(
            &sources,
            chromosomes,
            params.grouper_cfg.max_variant_group_span,
        );
        let dust_pool = if params.no_complexity_filter {
            None
        } else {
            Some(spawn_dust_pool(
                fasta_path,
                &chrom_plans,
                params.dust_cfg.window(),
                params.dust_cfg.threshold(),
            ))
        };

        Self {
            fasta_path,
            params,
            sources,
            n_samples,
            streaming: StreamingBlockLoader::with_n_samples(n_samples),
            partition_scratch: PartitionScratch::with_n_samples(n_samples),
            mask: Vec::new(),
            chrom_plans,
            next_plan_idx: 0,
            dust_pool,
            cur_interval_mask: Vec::new(),
            mask_sweep: 0,
            free_blocks: Vec::new(),
            next_seq_idx: 0,
            cur: None,
            chunks_loaded: 0,
            chunk_variants_total: 0,
        }
    }

    /// Return a spent block's buffers to the free-list so the next
    /// production reuses them instead of allocating.
    fn recycle(&mut self, block: ReadyBlock) {
        self.free_blocks.push(block);
    }

    /// Advance to the next chromosome that carries data, building its
    /// cursor (covered intervals + REF fetcher) positioned at interval 0.
    /// Returns `Ok(false)` when no chromosomes remain. Walks the
    /// precomputed [`chrom_plans`](Self::chrom_plans), so the interval
    /// sequence matches the DUST-ahead thread's. Does not touch the
    /// sources — [`enter_current_interval`](Self::enter_current_interval)
    /// resets them to interval 0.
    fn advance_chrom(&mut self) -> Result<bool, ChunkDriverError> {
        let Some(plan) = self.chrom_plans.get(self.next_plan_idx) else {
            self.cur = None;
            return Ok(false);
        };
        self.next_plan_idx += 1;
        let ref_fetcher = StreamingChromRefFetcher::for_contig(self.fasta_path, &plan.name)
            .map_err(|source| ChunkDriverError::OpenRefFetcher {
                contig: plan.name.clone(),
                source,
            })?;
        self.cur = Some(ChromCursor {
            chrom_id: plan.chrom_id,
            chrom_one_past_end: plan.length.saturating_add(1),
            ref_fetcher,
            covered: plan.covered.clone(),
            interval_idx: 0,
            // Filled by `enter_current_interval`.
            chunk_range_start: 0,
            interval_one_past_end: 0,
        });
        Ok(true)
    }

    /// Point the per-sample sources and the streaming loader at the
    /// cursor's current covered interval, set the block range to its
    /// start, and pull the interval's precomputed DUST mask off the
    /// DUST-ahead queue. Called when a chromosome or interval is
    /// (re-)entered; the interval boundary is wider than any group span,
    /// so the streaming loader starts with empty pending.
    fn enter_current_interval(&mut self) -> Result<(), ChunkDriverError> {
        let (chrom_id, interval, start, interval_one_past_end) = {
            let cur = self.cur.as_ref().expect("cur set before enter");
            let iv = &cur.covered[cur.interval_idx];
            (
                cur.chrom_id,
                iv.clone(),
                iv.start,
                iv.end.min(cur.chrom_one_past_end),
            )
        };
        let end_inclusive = interval_one_past_end.saturating_sub(1);
        for source in &mut self.sources {
            source.reset(chrom_id, start, end_inclusive);
        }
        self.streaming.reset_interval();

        // Take this interval's DUST mask from the DUST-ahead pool. The pool
        // delivers one mask per covered interval in the same genomic order
        // the producer enters them, so this is the mask for exactly this
        // interval (asserted below).
        if let Some(pool) = self.dust_pool.as_ref() {
            match pool.recv_next() {
                Some(Ok(interval_mask)) => {
                    debug_assert_eq!(interval_mask.chrom_id, chrom_id);
                    debug_assert_eq!(interval_mask.interval, interval);
                    self.cur_interval_mask = interval_mask.mask;
                    self.mask_sweep = 0;
                }
                Some(Err(e)) => return Err(e),
                None => return Err(ChunkDriverError::DustAheadGone { chrom_id }),
            }
        }

        let cur = self.cur.as_mut().expect("cur set before enter");
        cur.chunk_range_start = start;
        cur.interval_one_past_end = interval_one_past_end;
        Ok(())
    }

    /// Yield the next prepared block (owned), or `None` when the genome
    /// is exhausted. The block carries its own columns / partition /
    /// pre-fetched REF bytes and a monotonic `seq_idx` for ordered emit;
    /// recycle it with [`recycle`](Self::recycle) once consumed.
    fn next_block(&mut self) -> Option<Result<ReadyBlock, ChunkDriverError>> {
        loop {
            // Make sure a covered interval is entered.
            if self.cur.is_none() {
                match self.advance_chrom() {
                    Ok(true) => {
                        if let Err(e) = self.enter_current_interval() {
                            return Some(Err(e));
                        }
                    }
                    Ok(false) => return None,
                    Err(e) => return Some(Err(e)),
                }
            }
            // Produce one block from the current interval.
            let mut block = self
                .free_blocks
                .pop()
                .unwrap_or_else(|| ReadyBlock::with_n_samples(self.n_samples));
            match self.produce_block(&mut block) {
                Ok(variant_count) if variant_count > 0 => {
                    block.seq_idx = self.next_seq_idx;
                    self.next_seq_idx += 1;
                    self.chunks_loaded += 1;
                    self.chunk_variants_total += u64::from(variant_count);
                    return Some(Ok(block));
                }
                Ok(_) => {
                    // The current interval is exhausted: advance to the
                    // next interval (then chromosome) and keep going.
                    self.free_blocks.push(block);
                    let advanced_interval = match self.cur.as_mut() {
                        Some(cur) if cur.interval_idx + 1 < cur.covered.len() => {
                            cur.interval_idx += 1;
                            true
                        }
                        _ => false,
                    };
                    if advanced_interval {
                        if let Err(e) = self.enter_current_interval() {
                            return Some(Err(e));
                        }
                    } else {
                        match self.advance_chrom() {
                            Ok(true) => {
                                if let Err(e) = self.enter_current_interval() {
                                    return Some(Err(e));
                                }
                            }
                            Ok(false) => return None,
                            Err(e) => return Some(Err(e)),
                        }
                    }
                }
                Err(e) => {
                    self.free_blocks.push(block);
                    return Some(Err(e));
                }
            }
        }
    }

    /// Stream one block of the current interval into `block`: fold +
    /// compact to the variant target (the open group carries forward as
    /// the loader's pending), then slice the precomputed DUST-ahead mask,
    /// partition, and prefetch the REF bytes. Returns the block's variant
    /// count, or `0` when the interval is exhausted (no block produced).
    fn produce_block(&mut self, block: &mut ReadyBlock) -> Result<u32, ChunkDriverError> {
        let ReadyBlock {
            chunk,
            partition,
            pre_fetched_ref_bytes,
            ..
        } = block;
        let Self {
            params,
            sources,
            streaming,
            partition_scratch,
            mask,
            cur,
            cur_interval_mask,
            mask_sweep,
            ..
        } = self;
        let cur = cur
            .as_mut()
            .expect("produce_block needs a current interval");

        let max_group_span = params.grouper_cfg.max_variant_group_span;
        // `.max(1)` so a `0` (disabled) target still cuts at the first
        // closed group rather than yielding an empty block.
        let target_variants = params
            .sizing
            .target_variants_per_chunk
            .map_or(0, NonZeroU32::get)
            .max(1);
        let chrom_id = cur.chrom_id;
        let range_start = cur.chunk_range_start;
        let interval_end_inclusive = cur.interval_one_past_end.saturating_sub(1);

        // Stream one block of ≥ `target_variants` variable variants (or
        // the interval's tail), cutting at the open group's clean
        // boundary; the open group carries forward as `streaming`'s
        // pending. `0` means the interval is exhausted.
        let variant_count = streaming
            .fill_block(
                sources,
                chunk,
                chrom_id,
                range_start,
                interval_end_inclusive,
                target_variants,
            )
            .map_err(|source| ChunkDriverError::StreamRead { chrom_id, source })?;
        if variant_count == 0 {
            return Ok(0);
        }
        // The next block resumes where this one cut.
        cur.chunk_range_start = chunk.safe_end;

        // ── DUST mask (data-shaping). ──
        // Slice the block's span out of the interval mask the DUST-ahead
        // thread precomputed (taken in `enter_current_interval`), instead
        // of DUSTing inline — byte-identical, but off the critical path.
        mask.clear();
        let span_start = chunk.range.start;
        let span_end = chunk.safe_end;
        if !params.no_complexity_filter && span_start < span_end {
            slice_interval_mask(
                cur_interval_mask.as_slice(),
                mask_sweep,
                span_start,
                span_end,
                mask,
            );
        }
        // ── Partition + REF prefetch (data-shaping). ──
        let window = span_start..span_end;
        partition_window(
            chunk,
            &window,
            mask,
            max_group_span,
            partition_scratch,
            partition,
        )
        .map_err(|source| ChunkDriverError::Partition { chrom_id, source })?;
        prefetch_window_ref_bytes(partition, &cur.ref_fetcher, pre_fetched_ref_bytes)
            .map_err(|source| ChunkDriverError::PrefetchRefBytes { chrom_id, source })?;
        Ok(variant_count)
    }
}

/// Coalescing distance (bases) for block-index-driven covered
/// intervals. Adjacent or cross-sample-misaligned blocks within this
/// distance merge into one chunk segment, keeping the chunk count and
/// per-chunk overhead bounded when coverage is fragmented and absorbing
/// block-boundary misalignment across samples. Larger than any
/// realistic `max_group_span`, so coalescing only ever *merges* (loads
/// a little extra empty reference) and never splits a variant group —
/// it cannot change output. Gaps wider than this become safe chunk
/// boundaries the loop skips entirely.
const MIN_CHUNK_SPAN: u32 = 50_000;

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    //! M20: per-category filter tests for `emit_or_drop` and its
    //! helper predicate `record_fails_mapq_diff_t`. Pins each filter's
    //! boundary behaviour and the post-EM order in which `emit_or_drop`
    //! applies them. (The pre-EM `min_alt_obs_per_sample` predicate
    //! lives in `worker::columnar_passes_min_alt_obs`.)

    use super::*;
    use crate::pileup_record::AlleleSupportStats;
    use crate::psp::header::ParsedChromosome;
    use crate::var_calling::per_group_merger::MergedAllele;
    use crate::var_calling::posterior_engine::{EmDiagnostics, PosteriorRecord, RecordLocus};
    use crate::vcf::CohortMetadata;
    use tempfile::tempdir;

    // ── merge_block_ranges ───────────────────────────────────────────

    #[test]
    fn merge_block_ranges_empty_is_empty() {
        assert!(merge_block_ranges(std::iter::empty(), 50).is_empty());
    }

    #[test]
    fn merge_block_ranges_single_block_is_one_half_open_interval() {
        // (first=100, last=200) inclusive -> [100, 201) half-open.
        assert_eq!(merge_block_ranges([(100u32, 200u32)], 50), vec![100..201]);
    }

    #[test]
    fn merge_block_ranges_merges_within_max_group_span() {
        // gap = next.first - cur.last = 230 - 200 = 30 <= 50 -> merge.
        assert_eq!(
            merge_block_ranges([(100, 200), (230, 300)], 50),
            vec![100..301]
        );
    }

    #[test]
    fn merge_block_ranges_splits_when_gap_exceeds_max_group_span() {
        // gap = 260 - 200 = 60 > 50 -> separate (the gap is a safe boundary).
        assert_eq!(
            merge_block_ranges([(100, 200), (260, 300)], 50),
            vec![100..201, 260..301]
        );
    }

    #[test]
    fn merge_block_ranges_sorts_and_merges_overlapping_cross_sample_blocks() {
        // Unsorted, overlapping ranges (e.g. two samples covering the
        // same region with different block splits) collapse to one.
        assert_eq!(
            merge_block_ranges([(260, 300), (100, 200), (150, 270)], 50),
            vec![100..301]
        );
    }

    /// Boundary: gap exactly `max_group_span` merges; one more splits.
    #[test]
    fn merge_block_ranges_boundary_at_max_group_span() {
        assert_eq!(
            merge_block_ranges([(100, 200), (250, 300)], 50),
            vec![100..301]
        );
        assert_eq!(
            merge_block_ranges([(100, 200), (251, 300)], 50),
            vec![100..201, 251..301]
        );
    }

    // ── Stage 3: slice_interval_mask ─────────────────────────────────

    /// Consecutive blocks tiling an interval each get the runs that fall
    /// inside their span, and the sweep cursor advances monotonically so
    /// later blocks never re-scan consumed runs.
    #[test]
    fn slice_interval_mask_clips_and_advances_across_blocks() {
        let interval_mask = vec![10..20, 40..50, 80..100];
        let mut sweep = 0usize;

        let mut out = Vec::new();
        slice_interval_mask(&interval_mask, &mut sweep, 1, 30, &mut out);
        assert_eq!(out, vec![10..20]);

        out.clear();
        slice_interval_mask(&interval_mask, &mut sweep, 30, 60, &mut out);
        assert_eq!(out, vec![40..50]);

        out.clear();
        slice_interval_mask(&interval_mask, &mut sweep, 60, 200, &mut out);
        assert_eq!(out, vec![80..100]);
    }

    /// A masked run straddling a block boundary is clipped into both
    /// blocks (the cursor does not advance past it until a block fully
    /// clears it) — the case that makes the slice byte-identical to a
    /// per-block inline DUST over the same span.
    #[test]
    fn slice_interval_mask_reemits_run_straddling_block_boundary() {
        let interval_mask = vec![10..20, 45..70];
        let mut sweep = 0usize;

        let mut out = Vec::new();
        slice_interval_mask(&interval_mask, &mut sweep, 1, 50, &mut out);
        assert_eq!(out, vec![10..20, 45..50]);

        out.clear();
        slice_interval_mask(&interval_mask, &mut sweep, 50, 100, &mut out);
        assert_eq!(out, vec![50..70]);
    }

    /// A block whose span misses every run yields an empty slice and
    /// leaves the cursor where the next block can pick up.
    #[test]
    fn slice_interval_mask_empty_when_span_misses_every_run() {
        let interval_mask = vec![10..20, 80..90];
        let mut sweep = 0usize;
        let mut out = Vec::new();
        slice_interval_mask(&interval_mask, &mut sweep, 30, 40, &mut out);
        assert!(out.is_empty());
    }

    // ── Stage 3: dust_mask_for_interval ──────────────────────────────

    /// Write a single-contig unwrapped FASTA (+ sibling `.fai`) holding
    /// `seq` and return the tempdir guard + the FASTA path.
    fn write_unwrapped_fasta(seq: &[u8]) -> (tempfile::TempDir, PathBuf) {
        use std::io::Write;
        let dir = tempdir().expect("tempdir");
        let fa_path = dir.path().join("ref.fa");
        let mut f = std::fs::File::create(&fa_path).expect("create fasta");
        f.write_all(b">chr0\n").expect("write header");
        f.write_all(seq).expect("write seq");
        f.write_all(b"\n").expect("write newline");
        let header_len = b">chr0\n".len();
        std::fs::write(
            dir.path().join("ref.fa.fai"),
            // name, length, seq-offset, line_bases, line_width
            format!(
                "chr0\t{}\t{}\t{}\t{}\n",
                seq.len(),
                header_len,
                seq.len(),
                seq.len() + 1
            ),
        )
        .expect("write fai");
        (dir, fa_path)
    }

    /// Chunking a covered interval into sub-spans and coalescing the
    /// masked runs at the boundaries reproduces a single whole-interval
    /// `sdust_mask_for_span` exactly — for any sub-span size. This is the
    /// byte-identity argument for the DUST-ahead thread's sub-span loop.
    #[test]
    fn dust_mask_for_interval_matches_whole_interval_scan() {
        // High-complexity flanks (unmasked) around a long poly-A tract
        // (masked) that crosses the sub-span boundaries.
        let flank: &[u8] = b"ACGTCAGTACGATCAGTAGCATGCAGTAGCATCAGTACGAGCATCAGCAG";
        let mut seq = Vec::new();
        seq.extend_from_slice(flank);
        seq.extend(std::iter::repeat_n(b'A', 120));
        seq.extend_from_slice(flank);
        let len = seq.len() as u32;
        let (dir, fa_path) = write_unwrapped_fasta(&seq);

        let (window, threshold) = (64u32, 20u32);
        let interval = 1..len + 1;

        // One-shot reference scan over the whole interval.
        let mut whole_fetcher =
            ManualEvictChromRefFetcher::for_contig(&fa_path, "chr0").expect("fetcher");
        let (whole, _) = sdust_mask_for_span(
            interval.start,
            interval.end,
            len,
            window,
            threshold,
            MIN_DUST_HALO,
            |s, l| {
                whole_fetcher
                    .fetch(s, l)
                    .map(<[u8]>::to_vec)
                    .map_err(io::Error::other)
            },
        )
        .expect("whole scan");
        assert!(
            !whole.is_empty(),
            "fixture must mask the poly-A tract to exercise coalescing"
        );

        // Chunked + coalesced over a range of sub-span sizes (small ones
        // force many boundary merges; `len` is the whole-interval case).
        for subspan in [7u32, 13, 41, len] {
            let mut fetcher =
                ManualEvictChromRefFetcher::for_contig(&fa_path, "chr0").expect("fetcher");
            let chunked =
                dust_mask_for_interval(&mut fetcher, &interval, len, window, threshold, subspan)
                    .expect("chunked scan");
            assert_eq!(chunked, whole, "subspan={subspan}");
        }
        drop(dir);
    }

    // ── Stage 4: DustAheadPool coordinator ───────────────────────────
    //
    // These test the bounded-ordered coordinator with a *stub* compute
    // closure (no FASTA / no sdust): ordered delivery despite out-of-order
    // completion, look-ahead back-pressure, in-order error propagation,
    // parallelism-invariance, and clean shutdown.

    /// One stub job; the interval start doubles as the job's identity, so a
    /// delivered mask's `interval.start` reveals which job produced it.
    fn dust_job(idx: u32) -> IntervalJob {
        IntervalJob {
            chrom_id: 0,
            name: format!("c{idx}"),
            length: 1000,
            interval: idx..(idx + 1),
        }
    }

    fn dust_jobs(n: u32) -> Vec<IntervalJob> {
        (0..n).map(dust_job).collect()
    }

    fn ok_dust_mask(job: &IntervalJob) -> IntervalDustMask {
        IntervalDustMask {
            chrom_id: job.chrom_id,
            interval: job.interval.clone(),
            mask: Vec::new(),
        }
    }

    /// Masks are delivered in job (genomic) order even when the workers
    /// complete them in the exact reverse order — the reordering contract.
    #[test]
    fn dust_pool_delivers_in_order_despite_reverse_completion() {
        let njobs = 6u32;
        // Gate: job `s` may only finish after job `s+1` has — forces
        // strict reverse completion. `done[s]` set when job `s` finishes.
        let gate: Arc<(Mutex<Vec<bool>>, Condvar)> =
            Arc::new((Mutex::new(vec![false; njobs as usize + 1]), Condvar::new()));
        let gate_c = Arc::clone(&gate);
        let compute: Arc<DustComputeFn> = Arc::new(move |job: &IntervalJob| {
            let s = job.interval.start as usize;
            let (lock, cv) = &*gate_c;
            let mut done = lock.lock().unwrap();
            while s + 1 < njobs as usize && !done[s + 1] {
                done = cv.wait(done).unwrap();
            }
            done[s] = true;
            cv.notify_all();
            Ok(ok_dust_mask(job))
        });
        // n_workers = njobs and lookahead = njobs so all jobs run at once.
        let pool = DustAheadPool::spawn(dust_jobs(njobs), njobs as usize, njobs as usize, compute);
        for expected in 0..njobs {
            let mask = pool.recv_next().expect("delivered").expect("ok");
            assert_eq!(mask.interval.start, expected, "delivered out of order");
        }
        assert!(pool.recv_next().is_none());
    }

    /// Without the producer consuming, the pool dispatches at most
    /// `lookahead` jobs; each consumed mask frees exactly one slot.
    #[test]
    fn dust_pool_respects_lookahead_backpressure() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let started = Arc::new(AtomicUsize::new(0));
        let started_c = Arc::clone(&started);
        let compute: Arc<DustComputeFn> = Arc::new(move |job: &IntervalJob| {
            started_c.fetch_add(1, Ordering::SeqCst);
            Ok(ok_dust_mask(job))
        });
        let lookahead = 3usize;
        let pool = DustAheadPool::spawn(dust_jobs(20), lookahead, 4, compute);

        std::thread::sleep(std::time::Duration::from_millis(80));
        assert_eq!(
            started.load(Ordering::SeqCst),
            lookahead,
            "pool dispatched beyond the look-ahead before any delivery",
        );

        let _ = pool.recv_next().expect("first").expect("ok");
        std::thread::sleep(std::time::Duration::from_millis(80));
        assert_eq!(
            started.load(Ordering::SeqCst),
            lookahead + 1,
            "consuming one mask should free exactly one look-ahead slot",
        );
    }

    /// A compute error surfaces at its own seq, in order, without
    /// disturbing the masks before or after it.
    #[test]
    fn dust_pool_surfaces_compute_error_in_order() {
        let compute: Arc<DustComputeFn> = Arc::new(|job: &IntervalJob| {
            if job.interval.start == 2 {
                Err(ChunkDriverError::DustMaskIo {
                    contig: "boom".into(),
                    source: io::Error::other("boom"),
                })
            } else {
                Ok(ok_dust_mask(job))
            }
        });
        let pool = DustAheadPool::spawn(dust_jobs(5), 8, 4, compute);
        assert_eq!(pool.recv_next().unwrap().unwrap().interval.start, 0);
        assert_eq!(pool.recv_next().unwrap().unwrap().interval.start, 1);
        assert!(matches!(
            pool.recv_next().unwrap().unwrap_err(),
            ChunkDriverError::DustMaskIo { .. }
        ));
        assert_eq!(pool.recv_next().unwrap().unwrap().interval.start, 3);
        assert_eq!(pool.recv_next().unwrap().unwrap().interval.start, 4);
        assert!(pool.recv_next().is_none());
    }

    /// The delivered sequence is identical whether one worker or many
    /// compute the masks — parallelism changes timing, never order.
    #[test]
    fn dust_pool_delivery_is_parallelism_invariant() {
        let run = |n_workers: usize| -> Vec<u32> {
            let compute: Arc<DustComputeFn> = Arc::new(|job: &IntervalJob| Ok(ok_dust_mask(job)));
            let pool = DustAheadPool::spawn(dust_jobs(12), 24, n_workers, compute);
            let mut out = Vec::new();
            while let Some(r) = pool.recv_next() {
                out.push(r.unwrap().interval.start);
            }
            out
        };
        let one = run(1);
        let many = run(8);
        assert_eq!(one, (0..12).collect::<Vec<_>>());
        assert_eq!(one, many, "delivery order differs by worker count");
    }

    /// Dropping the pool mid-flight (workers still computing) shuts down
    /// and joins cleanly — the test simply completing proves no hang.
    #[test]
    fn dust_pool_shuts_down_cleanly_mid_flight() {
        let compute: Arc<DustComputeFn> = Arc::new(|job: &IntervalJob| {
            std::thread::sleep(std::time::Duration::from_millis(5));
            Ok(ok_dust_mask(job))
        });
        let pool = DustAheadPool::spawn(dust_jobs(50), 8, 4, compute);
        let _ = pool.recv_next().expect("one").unwrap();
        drop(pool); // Drop → shutdown → join; reaching the end == clean.
    }

    fn ref_allele() -> MergedAllele {
        MergedAllele {
            seq: b"A".to_vec(),
            is_compound: false,
            constituents: Vec::new(),
        }
    }

    fn alt_allele() -> MergedAllele {
        MergedAllele {
            seq: b"T".to_vec(),
            is_compound: false,
            constituents: Vec::new(),
        }
    }

    /// Build a 2-sample biallelic `PosteriorRecord` with caller-
    /// supplied per-sample (REF_num_obs, ALT_num_obs) and per-sample
    /// (best_genotype_idx, qual_phred). Other fields are filled with
    /// stable defaults so the test surface stays small.
    #[allow(clippy::too_many_arguments)]
    fn mk_record(
        pos: u32,
        per_sample_num_obs: &[(u32, u32)],
        best_genotype: Vec<usize>,
        qual_phred: f64,
        per_sample_alt_mapq_sum: u32,
        per_sample_alt_mapq_sum_sq: u64,
    ) -> PosteriorRecord {
        let n_samples = per_sample_num_obs.len();
        let scalars: Vec<AlleleSupportStats> = per_sample_num_obs
            .iter()
            .flat_map(|&(r, a)| {
                vec![
                    // REF: keep mapq moments low so REF-vs-ALT
                    // difference is the test knob.
                    AlleleSupportStats::new(r, 0.0, r, 0, 0, r * 30, u64::from(r) * 900),
                    // ALT: caller-supplied moments so the
                    // mapq-diff-t test can be triggered or skipped.
                    AlleleSupportStats::new(
                        a,
                        0.0,
                        a,
                        0,
                        0,
                        per_sample_alt_mapq_sum,
                        per_sample_alt_mapq_sum_sq,
                    ),
                ]
            })
            .collect();
        PosteriorRecord {
            locus: RecordLocus {
                chrom_id: 0,
                start: pos,
                end: pos,
            },
            alleles: vec![ref_allele(), alt_allele()],
            ploidy: 2,
            n_samples,
            n_genotypes: 3,
            allele_frequencies: vec![0.5, 0.5],
            compound_frequencies: vec![None, None],
            posteriors: vec![0.0; n_samples * 3],
            best_genotype,
            gq_phred: vec![60.0; n_samples],
            qual_phred,
            scalars,
            other_scalars: vec![],
            chain_anchor_flags: vec![false; n_samples * 2],
            diagnostics: EmDiagnostics {
                iterations: 5,
                final_max_delta_p: 1e-6,
                converged: true,
            },
        }
    }

    // (The `min_alt_obs_per_sample` predicate now lives pre-EM in
    // `worker::columnar_passes_min_alt_obs`, tested there.)

    // ── record_fails_mapq_diff_t ─────────────────────────────────────

    /// Non-finite threshold (`NaN`/`±INF`) short-circuits to `false`
    /// (filter inert).
    #[test]
    fn record_fails_mapq_diff_t_returns_false_when_threshold_is_nan() {
        let r = mk_record(100, &[(10, 10)], vec![0, 0], 100.0, 100, 4000);
        assert!(!record_fails_mapq_diff_t(&r, f32::NAN));
    }

    /// REF count below the minimum (n_ref < 3) ⇒ filter inert
    /// regardless of ALT moments.
    #[test]
    fn record_fails_mapq_diff_t_returns_false_when_n_ref_below_minimum() {
        let r = mk_record(100, &[(1, 10)], vec![0, 0], 100.0, 100, 4000);
        assert!(!record_fails_mapq_diff_t(&r, -1.0));
    }

    /// Both sides ≥ MAPQ_FILTER_MIN_READS_PER_SIDE; zero pooled
    /// variance (every per-read mapq is identical) ⇒ `se2 == 0.0` ⇒
    /// the helper short-circuits with `continue` and never trips.
    /// Documents the `se2 <= 0.0` early-exit branch the Welch's-t
    /// formula needs against degenerate cohorts.
    #[test]
    fn record_fails_mapq_diff_t_skips_alts_with_zero_pooled_variance() {
        // REF mapq is 30 per read; ALT mapq is 1 per read; both
        // are constant across samples ⇒ variance is zero on both
        // sides ⇒ se2 == 0.0 ⇒ no trip even though means differ
        // sharply.
        let r = mk_record(100, &[(10, 10)], vec![0, 0], 100.0, 10, 10);
        assert!(!record_fails_mapq_diff_t(&r, 0.0));
    }

    // ── emit_or_drop — filter order ──────────────────────────────────

    /// Build a minimal `ChunkDriverParams` with explicit drop
    /// thresholds. Other fields use the in-tree config defaults.
    fn params_for_filters(
        min_qual_phred: f64,
        min_alt_obs_per_sample: u32,
        min_mapq_diff_t: f32,
        no_mapq_diff_filter: bool,
    ) -> ChunkDriverParams {
        use crate::var_calling::dust_filter::DustFilterConfig;
        use crate::var_calling::per_group_merger::PerGroupMergerConfig;
        use crate::var_calling::posterior_engine::PosteriorEngineConfig;
        use crate::var_calling::variant_grouping::GrouperConfig;
        ChunkDriverParams {
            no_complexity_filter: true,
            dust_cfg: DustFilterConfig::new(64, 20).expect("dust"),
            grouper_cfg: GrouperConfig::new(50).expect("grouper"),
            per_group_cfg: PerGroupMergerConfig::default(),
            posterior_cfg: PosteriorEngineConfig::default(),
            sizing: ChunkSizingParams {
                chunk_genomic_span: DEFAULT_CHUNK_GENOMIC_SPAN,
                target_variants_per_chunk: None,
            },
            downstream: DownstreamFilterParams {
                min_alt_obs_per_sample,
                min_qual_phred,
                no_mapq_diff_filter,
                min_mapq_diff_t,
            },
        }
    }

    fn open_writer(out: &Path) -> CohortVcfWriter {
        let metadata = CohortMetadata {
            sample_names: vec!["S0".into(), "S1".into()],
            contigs: vec![ParsedChromosome {
                name: "chr1".into(),
                length: 1_000_000,
                md5: "00000000000000000000000000000001".into(),
            }],
            tool_string: "pop_var_caller test".into(),
            command_line: String::new(),
        };
        let cfg = WriterConfig::new(out.to_path_buf());
        CohortVcfWriter::new(metadata, cfg).expect("writer")
    }

    /// `is_variant_call` returns false for all-hom-REF best_genotype
    /// ⇒ records_dropped_hom_ref bumps.
    #[test]
    fn emit_or_drop_increments_hom_ref_when_is_variant_call_returns_false() {
        let dir = tempdir().unwrap();
        let mut writer = open_writer(&dir.path().join("out.vcf"));
        let params = params_for_filters(0.0, 0, f32::NEG_INFINITY, true);
        let mut stats = ChunkDriverStats::default();
        // All-hom-REF: best_genotype index 0 in every sample.
        let r = mk_record(100, &[(10, 5), (10, 5)], vec![0, 0], 100.0, 0, 0);
        emit_or_drop(r, &params, &mut writer, &mut stats).unwrap();
        assert_eq!(stats.records_dropped_hom_ref, 1);
        assert_eq!(stats.records_dropped_low_alt_obs, 0);
        assert_eq!(stats.records_dropped_low_qual, 0);
        assert_eq!(stats.records_dropped_low_mapq_diff_t, 0);
        let _ = writer.abort();
    }

    /// `qual_phred < min_qual_phred` after the variant-call check
    /// passes ⇒ records_dropped_low_qual bumps.
    #[test]
    fn emit_or_drop_increments_low_qual_when_qual_phred_below_threshold() {
        let dir = tempdir().unwrap();
        let mut writer = open_writer(&dir.path().join("out.vcf"));
        let params = params_for_filters(50.0, 0, f32::NEG_INFINITY, true);
        let mut stats = ChunkDriverStats::default();
        // Variant call (best_genotype != [0,0]); qual = 10 < 50.
        let r = mk_record(100, &[(10, 5), (10, 5)], vec![1, 1], 10.0, 0, 0);
        emit_or_drop(r, &params, &mut writer, &mut stats).unwrap();
        assert_eq!(stats.records_dropped_low_qual, 1);
        assert_eq!(stats.records_dropped_low_alt_obs, 0);
        assert_eq!(stats.records_dropped_hom_ref, 0);
        assert_eq!(stats.records_dropped_low_mapq_diff_t, 0);
        let _ = writer.abort();
    }

    /// Filter order: a record that *would* trip both `is_variant_call`
    /// (hom_ref) and `qual_phred` (low qual) increments **only**
    /// `records_dropped_hom_ref` — the post-EM order is
    /// `is_variant_call → qual_phred → mapq_diff_t`, and hom_ref
    /// short-circuits first. (`min_alt_obs` runs earlier, pre-EM.)
    #[test]
    fn emit_or_drop_filter_order_pins_hom_ref_winning_over_low_qual() {
        let dir = tempdir().unwrap();
        let mut writer = open_writer(&dir.path().join("out.vcf"));
        let params = params_for_filters(50.0, 0, f32::NEG_INFINITY, true);
        let mut stats = ChunkDriverStats::default();
        // hom-REF best genotype AND qual below threshold — both
        // filters would trip; hom_ref's branch is reached first.
        let r = mk_record(100, &[(10, 1), (10, 1)], vec![0, 0], 10.0, 0, 0);
        emit_or_drop(r, &params, &mut writer, &mut stats).unwrap();
        assert_eq!(
            stats.records_dropped_hom_ref, 1,
            "hom_ref wins the order race"
        );
        assert_eq!(stats.records_dropped_low_qual, 0);
        let _ = writer.abort();
    }

    /// All filters pass ⇒ records_written bumps.
    #[test]
    fn emit_or_drop_writes_record_and_increments_records_written_on_pass() {
        let dir = tempdir().unwrap();
        let mut writer = open_writer(&dir.path().join("out.vcf"));
        let params = params_for_filters(0.0, 0, f32::NEG_INFINITY, true);
        let mut stats = ChunkDriverStats::default();
        let r = mk_record(100, &[(10, 5), (10, 5)], vec![1, 1], 100.0, 0, 0);
        emit_or_drop(r, &params, &mut writer, &mut stats).unwrap();
        assert_eq!(stats.records_written, 1);
        assert_eq!(stats.records_dropped_low_alt_obs, 0);
        assert_eq!(stats.records_dropped_hom_ref, 0);
        assert_eq!(stats.records_dropped_low_qual, 0);
        assert_eq!(stats.records_dropped_low_mapq_diff_t, 0);
        // finish() consumes the writer cleanly.
        writer.finish().unwrap();
    }
}
