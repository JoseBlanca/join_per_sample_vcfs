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
//! [`drive_cohort_pipeline`](crate::pop_var_caller::cohort_driver::drive_cohort_pipeline).
//! The same set of records is dropped from the final VCF, with
//! the same per-category counters incremented — that's the
//! byte-identity contract for Phase A.

use std::collections::{BTreeMap, VecDeque};
use std::fs::File;
use std::io::{self, BufReader, Read, Seek};
use std::num::NonZeroU32;
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::sync::{Condvar, Mutex};
use thiserror::Error;

use crate::fasta::fetcher::ChromRefFetchError;
use crate::fasta::{ManualEvictChromRefFetcher, StreamingChromRefFetcher};
use crate::psp::header::ParsedChromosome;
use crate::psp::{PspReadError, PspReader};
use crate::var_calling::cohort_block::column_span_reader::ColumnSpanReader;
use crate::var_calling::cohort_block::columns::MaterialisedChunk;
use crate::var_calling::cohort_block::loader::{ChunkLoadError, StreamingBlockLoader};
use crate::var_calling::cohort_block::partition::{
    PartitionError, PartitionScratch, WindowPartition, partition_window,
};
use crate::var_calling::cohort_block::worker::{
    WindowRunStats, WorkerScratch, prefetch_window_ref_bytes, run_window,
};
use crate::var_calling::dust_filter::{DustFilterConfig, MIN_DUST_HALO, sdust_mask_for_span};
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

/// Mi19: post-EM downstream filters applied by `emit_or_drop` before
/// each record reaches the VCF writer. Grouped separately from
/// `ChunkDriverParams` so the per-record filter axis is named and
/// `emit_or_drop`'s signature can take `&DownstreamFilterParams`
/// instead of the whole driver bundle.
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
/// [`CohortDriveStats`](crate::pop_var_caller::cohort_driver::CohortDriveStats)
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
            process_block(&block, per_group_cfg, posterior_cfg, &mut worker)
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
        for _ in 0..n_workers {
            let queue_ref = &queue;
            let result_tx = result_tx.clone();
            let recycle_tx = recycle_tx.clone();
            scope.spawn(move || {
                let mut worker = WorkerScratch::new();
                while let Some(block) = queue_ref.pop() {
                    match process_block(&block, per_group_cfg, posterior_cfg, &mut worker) {
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
    if !passes_min_alt_obs(&record, params.downstream.min_alt_obs_per_sample) {
        stats.records_dropped_low_alt_obs += 1;
        return Ok(());
    }
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

/// `min_alt_obs_per_sample` filter — same predicate as
/// `drive_cohort_pipeline`'s `.filter_map` between merger and EM.
/// Phase A.0 applies it post-EM (the EM has already run on the
/// doomed record); Phase A.1's native-columnar kernel can push it
/// back to the pre-EM boundary if perf review shows the wasted EM
/// work matters.
fn passes_min_alt_obs(record: &PosteriorRecord, min_obs: u32) -> bool {
    if min_obs <= 1 || record.alleles.len() <= 1 {
        return true;
    }
    let n_alleles = record.alleles.len();
    (1..n_alleles).any(|allele_idx| {
        let mut max_obs = 0_u32;
        for sample_idx in 0..record.n_samples {
            let obs = record.scalars[sample_idx * n_alleles + allele_idx].num_obs;
            if obs > max_obs {
                max_obs = obs;
            }
        }
        max_obs >= min_obs
    })
}

/// Welch's-t MAPQ-difference filter — direct port of the streaming
/// driver's `record_fails_mapq_diff_t` so the chunk driver can
/// match its drop semantics without a cross-module call.
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
    scratch: &mut WorkerScratch,
) -> Result<(), PosteriorEngineError> {
    run_window(
        &block.chunk,
        &block.partition,
        &block.pre_fetched_ref_bytes,
        per_group_cfg,
        posterior_cfg,
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
    chrom_name: String,
    chrom_length: u32,
    chrom_one_past_end: u32,
    /// Forward-only fetcher feeding the per-block REF prefetch.
    ref_fetcher: StreamingChromRefFetcher,
    /// Random-access fetcher feeding the per-block DUST mask; `None`
    /// when the complexity filter is disabled.
    dust_fetcher: Option<ManualEvictChromRefFetcher>,
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
/// met, cuts at a clean group boundary (reserving the tail), then
/// DUST-masks + partitions + prefetches REF bytes for the block.
///
/// This is the *data-shaping* half. `next_block` yields an owned
/// [`ReadyBlock`]; the consumer runs the per-group math via
/// [`process_block`] and emits — never seeing how the block was made.
/// Genomic windows are gone: one block = one unit of work. Spent blocks
/// are returned via [`recycle`](Self::recycle) so their buffers feed the
/// next production without reallocating.
struct BlockIterator<'a, W: Read + Seek + Send> {
    chromosomes: &'a [ParsedChromosome],
    fasta_path: &'a Path,
    params: &'a ChunkDriverParams,
    /// One span-addressable columnar reader per sample; each owns its
    /// `PspReader`. Re-`reset` to the current covered interval.
    sources: Vec<ColumnSpanReader<W>>,
    n_samples: usize,

    /// Streaming fold+compact engine: holds the per-sample pending
    /// (open-group straddle) and the fold scratch, reused across blocks.
    streaming: StreamingBlockLoader,
    /// Reusable per-block partition scratch + DUST-mask buffer.
    partition_scratch: PartitionScratch,
    mask: Vec<std::ops::Range<u32>>,

    // Free-list of recycled block buffers (columns + CSR arrays + REF
    // byte buffers). `next_block` pops one (or makes a fresh one) to
    // fill; the consumer pushes spent blocks back via `recycle`.
    free_blocks: Vec<ReadyBlock>,
    next_seq_idx: u64,

    // Production cursor + counters.
    cur: Option<ChromCursor>,
    next_chrom_idx: usize,
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
        let sources = psp_readers.into_iter().map(ColumnSpanReader::detached).collect();
        Self {
            chromosomes,
            fasta_path,
            params,
            sources,
            n_samples,
            streaming: StreamingBlockLoader::with_n_samples(n_samples),
            partition_scratch: PartitionScratch::with_n_samples(n_samples),
            mask: Vec::new(),
            free_blocks: Vec::new(),
            next_seq_idx: 0,
            cur: None,
            next_chrom_idx: 0,
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
    /// cursor (covered intervals + fetchers) positioned at interval 0.
    /// Returns `Ok(false)` when no chromosomes remain. Does not touch
    /// the sources — [`enter_current_interval`](Self::enter_current_interval)
    /// resets them to interval 0.
    fn advance_chrom(&mut self) -> Result<bool, ChunkDriverError> {
        let max_group_span = self.params.grouper_cfg.max_variant_group_span;
        while self.next_chrom_idx < self.chromosomes.len() {
            let chrom_idx = self.next_chrom_idx;
            self.next_chrom_idx += 1;
            let chrom = &self.chromosomes[chrom_idx];
            let chrom_length = chrom.length;
            if chrom_length == 0 {
                continue;
            }
            // M25 PANIC-FREE: chrom_id fits in u32 (PSP file-format invariant).
            let chrom_id = u32::try_from(chrom_idx).expect("chrom_id fits in u32");
            let covered = covered_intervals_for_chrom(
                &self.sources,
                chrom_id,
                max_group_span.max(MIN_CHUNK_SPAN),
            );
            if covered.is_empty() {
                // No sample carries data here: open no fetchers, read nothing.
                continue;
            }
            let ref_fetcher = StreamingChromRefFetcher::for_contig(self.fasta_path, &chrom.name)
                .map_err(|source| ChunkDriverError::OpenRefFetcher {
                    contig: chrom.name.clone(),
                    source,
                })?;
            let dust_fetcher: Option<ManualEvictChromRefFetcher> =
                if self.params.no_complexity_filter {
                    None
                } else {
                    Some(
                        ManualEvictChromRefFetcher::for_contig(self.fasta_path, &chrom.name)
                            .map_err(|source| ChunkDriverError::OpenRefFetcher {
                                contig: chrom.name.clone(),
                                source,
                            })?,
                    )
                };
            self.cur = Some(ChromCursor {
                chrom_id,
                chrom_name: chrom.name.clone(),
                chrom_length,
                chrom_one_past_end: chrom_length.saturating_add(1),
                ref_fetcher,
                dust_fetcher,
                covered,
                interval_idx: 0,
                // Filled by `enter_current_interval`.
                chunk_range_start: 0,
                interval_one_past_end: 0,
            });
            return Ok(true);
        }
        self.cur = None;
        Ok(false)
    }

    /// Point the per-sample sources and the streaming loader at the
    /// cursor's current covered interval, and set the block range to its
    /// start. Called when a chromosome or interval is (re-)entered; the
    /// interval boundary is wider than any group span, so the streaming
    /// loader starts with empty pending.
    fn enter_current_interval(&mut self) {
        let (chrom_id, start, interval_one_past_end) = {
            let cur = self.cur.as_ref().expect("cur set before enter");
            let iv = &cur.covered[cur.interval_idx];
            (cur.chrom_id, iv.start, iv.end.min(cur.chrom_one_past_end))
        };
        let end_inclusive = interval_one_past_end.saturating_sub(1);
        for source in &mut self.sources {
            source.reset(chrom_id, start, end_inclusive);
        }
        self.streaming.reset_interval();
        let cur = self.cur.as_mut().expect("cur set before enter");
        cur.chunk_range_start = start;
        cur.interval_one_past_end = interval_one_past_end;
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
                    Ok(true) => self.enter_current_interval(),
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
                        self.enter_current_interval();
                    } else {
                        match self.advance_chrom() {
                            Ok(true) => self.enter_current_interval(),
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
    /// the loader's pending), then DUST-mask, partition, and prefetch
    /// the REF bytes. Returns the block's variant count, or `0` when the
    /// interval is exhausted (no block produced).
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
        let chrom_length = cur.chrom_length;
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
        mask.clear();
        let span_start = chunk.range.start;
        let span_end = chunk.safe_end;
        if let Some(dust_fetcher) = cur.dust_fetcher.as_mut()
            && span_start < span_end
        {
            let (m, buf_start) = sdust_mask_for_span(
                span_start,
                span_end,
                chrom_length,
                params.dust_cfg.window(),
                params.dust_cfg.threshold(),
                MIN_DUST_HALO,
                |s, l| {
                    dust_fetcher
                        .fetch(s, l)
                        .map(<[u8]>::to_vec)
                        .map_err(io::Error::other)
                },
            )
            .map_err(|source| ChunkDriverError::DustMaskIo {
                contig: cur.chrom_name.clone(),
                source,
            })?;
            *mask = m;
            dust_fetcher.evict_before(buf_start);
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
    //! helper predicates `passes_min_alt_obs` and
    //! `record_fails_mapq_diff_t`. Pins each filter's boundary
    //! behaviour and the order in which `emit_or_drop` applies them.

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

    // ── passes_min_alt_obs ───────────────────────────────────────────

    /// `min_obs <= 1` short-circuits to `true`.
    #[test]
    fn passes_min_alt_obs_returns_true_when_min_obs_is_zero() {
        let r = mk_record(100, &[(10, 0)], vec![0, 0], 100.0, 0, 0);
        assert!(passes_min_alt_obs(&r, 0));
    }

    /// REF-only record (n_alleles == 1) short-circuits to `true`.
    #[test]
    fn passes_min_alt_obs_returns_true_when_record_has_only_ref_allele() {
        let mut r = mk_record(100, &[(10, 0)], vec![0, 0], 100.0, 0, 0);
        r.alleles = vec![ref_allele()];
        r.scalars = vec![AlleleSupportStats::new(10, 0.0, 10, 0, 0, 0, 0)];
        assert!(passes_min_alt_obs(&r, 5));
    }

    /// All samples have ALT obs < min_obs ⇒ false.
    #[test]
    fn passes_min_alt_obs_returns_false_when_all_alts_below_threshold() {
        let r = mk_record(100, &[(10, 1), (10, 2)], vec![0, 0], 100.0, 0, 0);
        assert!(!passes_min_alt_obs(&r, 5));
    }

    /// At least one sample's ALT obs >= min_obs ⇒ true (max-across-
    /// samples wins).
    #[test]
    fn passes_min_alt_obs_returns_true_when_one_alt_meets_threshold_in_one_sample() {
        let r = mk_record(100, &[(10, 1), (10, 7)], vec![0, 0], 100.0, 0, 0);
        assert!(passes_min_alt_obs(&r, 5));
    }

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

    /// `emit_or_drop` increments `records_dropped_low_alt_obs` (not
    /// any other counter) when the record fails `min_alt_obs`.
    #[test]
    fn emit_or_drop_increments_low_alt_obs_when_record_fails_min_alt_obs() {
        let dir = tempdir().unwrap();
        let mut writer = open_writer(&dir.path().join("out.vcf"));
        let params = params_for_filters(0.0, 5, f32::NEG_INFINITY, true);
        let mut stats = ChunkDriverStats::default();
        // Both samples have ALT obs = 1 < min_alt_obs_per_sample = 5.
        let r = mk_record(100, &[(10, 1), (10, 1)], vec![1, 1], 100.0, 0, 0);
        emit_or_drop(r, &params, &mut writer, &mut stats).unwrap();
        assert_eq!(stats.records_dropped_low_alt_obs, 1);
        assert_eq!(stats.records_dropped_hom_ref, 0);
        assert_eq!(stats.records_dropped_low_qual, 0);
        assert_eq!(stats.records_dropped_low_mapq_diff_t, 0);
        assert_eq!(stats.records_written, 0);
        let _ = writer.abort();
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
    /// `records_dropped_hom_ref` — the order is `min_alt_obs →
    /// is_variant_call → qual_phred → mapq_diff_t`, and hom_ref
    /// short-circuits first.
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
