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
//! scratch buffers across every chunk and every chromosome
//! (`ChunkLoadScratch`, `BoundaryFinalisationScratch`, `PartitionScratch`,
//! the `MaterialisedChunk` + `WindowPartition` buffers, the
//! carryover, the per-window `Vec<PosteriorRecord>`). All of them
//! are `clear()`-ed on every iteration and re-filled, never
//! reallocated. This is the bulk of the memory-budget win the
//! rewrite promises: one chunk × N samples, not T_chrom × N samples.
//!
//! **Downstream filters** (`is_variant_call`, `qual_phred`,
//! `min_alt_obs_per_sample`, MAPQ-diff t-test) are applied
//! post-EM in the driver, mirroring the streaming pipeline's drop
//! semantics in
//! [`drive_cohort_pipeline`](crate::pop_var_caller::cohort_driver::drive_cohort_pipeline).
//! The same set of records is dropped from the final VCF, with
//! the same per-category counters incremented — that's the
//! byte-identity contract for Phase A.

use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek};
use std::num::{NonZeroU32, NonZeroUsize};
use std::path::{Path, PathBuf};
use thiserror::Error;

use crate::fasta::StreamingChromRefFetcher;
use crate::fasta::fetcher::ChromRefFetchError;
use crate::psp::header::ParsedChromosome;
use crate::psp::{PspReadError, PspReader};
use crate::var_calling::cohort_block::chunk_boundaries::{
    BoundaryFinalisationError, BoundaryFinalisationScratch, finalise_chunk_boundaries,
};
use crate::var_calling::cohort_block::columns::{MaterialisedChunk, SampleColumns};
use crate::var_calling::cohort_block::loader::{
    ChunkLoadError, ChunkLoadExtent, ChunkLoadScratch, load_chunk_from_iters,
};
use crate::var_calling::cohort_block::partition::{PartitionError, partition_window};
use crate::var_calling::cohort_block::worker::{WorkerPool, prefetch_window_ref_bytes, run_window};
use crate::var_calling::dust_filter::{DustFilterConfig, sdust_mask_streaming};
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
    /// M7: number of parallel worker windows per chunk.
    /// [`finalise_chunk_boundaries`] places `target_window_count - 1` internal
    /// boundaries inside each chunk's `[range.start, safe_end)` so
    /// the [`WorkerPool`] can dispatch the per-window math
    /// concurrently. `NonZeroUsize::MIN` (1) preserves the
    /// sequential single-window-per-chunk behaviour byte-for-byte;
    /// `0` is rejected at the type level. The previous shape was a
    /// bare `usize` where the driver silently mapped `0` to `1` via
    /// `.max(1)`, bypassing the pre-pass's
    /// `BoundaryFinalisationError::ZeroTargetWindowCount` validation.
    /// See [the Phase B plan](https://github.com/JoseBlanca/join_per_sample_vcfs/blob/main/doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_b_parallel_windows.md).
    pub target_window_count: NonZeroUsize,
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
    /// Mi21: chunks where `finalise_chunk_boundaries` produced fewer windows
    /// than `params.sizing.target_window_count` requested. The pre-pass
    /// silently down-grades when safe positions are sparse; an
    /// operator setting `--worker-windows-per-chunk 8` and seeing
    /// `chunks_with_fewer_windows_than_requested > 0` knows the
    /// requested parallelism didn't fully realise on those chunks.
    pub chunks_with_fewer_windows_than_requested: u64,
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
/// - **Per-chunk pipeline**: `LoadChunk`, `FinaliseChunkBoundaries`, `Partition`
///   (each carries `chrom_id`), `RunWindow` (per-window math via
///   `run_window`).
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
    #[error("failed to commit chunk boundaries on chromosome {chrom_id}")]
    FinaliseChunkBoundaries {
        chrom_id: u32,
        #[source]
        source: BoundaryFinalisationError,
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
         target_window_count={} min_qual_phred={} min_alt_obs_per_sample={} \
         no_mapq_diff_filter={} min_mapq_diff_t={} no_complexity_filter={}",
        params.sizing.chunk_genomic_span,
        params
            .sizing
            .target_variants_per_chunk
            .map_or(0, NonZeroU32::get),
        params.sizing.target_window_count.get(),
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

    // Persistent buffers (scratch reuse across every chunk + every chrom).
    let mut chunk_scratch = ChunkLoadScratch::with_n_samples(n_samples);
    let mut fix_scratch = BoundaryFinalisationScratch::new();
    let mut chunk = MaterialisedChunk::with_n_samples(n_samples);
    let mut carryover: Vec<SampleColumns> =
        (0..n_samples).map(|_| SampleColumns::empty()).collect();
    let mut worker_pool = WorkerPool::empty();

    let mut stats = ChunkDriverStats::default();

    // B1 + M5: drive the per-chrom loop inside a closure that borrows
    // `writer` (does not move it), so the outer scope can route to
    // `writer.finish()` on success or `writer.abort()` on error. The
    // tmp file is removed via `abort()` at the exact path the writer
    // used (single source of truth); a remove failure is logged
    // rather than silently swallowed.
    let driver_result: Result<(), ChunkDriverError> = (|| {
        for (chrom_idx, chrom) in chromosomes.iter().enumerate() {
            // M25 PANIC-FREE: `chrom_idx` comes from
            // `chromosomes.iter().enumerate()` where `chromosomes`
            // originates from the PSP header. PSP `chrom_id`s are
            // `u32` by file-format definition, so the contig table
            // never exceeds `u32::MAX` entries.
            let chrom_id = u32::try_from(chrom_idx).expect("chrom_id fits in u32");
            drive_one_chrom_generic(
                chrom_id,
                chrom,
                fasta_path,
                &mut psp_readers,
                &params,
                &mut writer,
                &mut stats,
                &mut chunk_scratch,
                &mut fix_scratch,
                &mut chunk,
                &mut carryover,
                &mut worker_pool,
            )?;
        }
        Ok(())
    })();

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

/// Compute the DUST low-complexity mask for one chromosome by
/// streaming its reference bases through `sdust_mask_streaming`
/// once, then translating the resulting BED-style 0-based
/// half-open slice intervals into 1-based half-open genomic
/// intervals the partition expects.
fn compute_dust_mask_for_chrom(
    ref_fetcher: &StreamingChromRefFetcher,
    contig: &str,
    chrom_length: u32,
    dust_cfg: &DustFilterConfig,
) -> Result<Vec<std::ops::Range<u32>>, ChunkDriverError> {
    use crate::fasta::fetcher::ChromRefFetcher;

    // B2: stream `iter_bases()` directly into `sdust_mask_streaming`
    // instead of materialising the whole chromosome's REF into a
    // `Vec<u8>`. The previous shape allocated O(chrom_length) bytes
    // (~90 MB for tomato chrom 1) just to re-iterate them — defeating
    // the chunk driver's "per-chunk memory bounded by
    // `target_variants_per_chunk × n_samples`" contract.
    let base_iter =
        ref_fetcher
            .iter_bases()
            .map_err(|source| ChunkDriverError::ComputeDustMask {
                contig: contig.to_string(),
                source,
            })?;
    let intervals = sdust_mask_streaming(
        base_iter.map(|r| r.map_err(io::Error::other)),
        chrom_length,
        dust_cfg.window(),
        dust_cfg.threshold(),
    )
    .map_err(|source| ChunkDriverError::DustMaskIo {
        contig: contig.to_string(),
        source,
    })?;

    Ok(intervals
        .into_iter()
        .map(|(s, e)| (s + 1)..(e + 1))
        .collect())
}

/// Drive the chunk loop over one chromosome. Pure helper for
/// [`drive_cohort_chunked`]; broken out so the per-chrom borrow
/// graph is easier to read and so the function stays generic over
/// the `PspReader`'s underlying `Read + Seek` source for testability.
#[allow(clippy::too_many_arguments)]
fn drive_one_chrom_generic<W>(
    chrom_id: u32,
    chrom: &ParsedChromosome,
    fasta_path: &Path,
    psp_readers: &mut [PspReader<W>],
    params: &ChunkDriverParams,
    writer: &mut CohortVcfWriter,
    stats: &mut ChunkDriverStats,
    chunk_scratch: &mut ChunkLoadScratch,
    fix_scratch: &mut BoundaryFinalisationScratch,
    chunk: &mut MaterialisedChunk,
    carryover: &mut [SampleColumns],
    worker_pool: &mut WorkerPool,
) -> Result<(), ChunkDriverError>
where
    W: Read + Seek,
{
    let chrom_length = chrom.length;
    if chrom_length == 0 {
        return Ok(());
    }

    // Mi6: per-chrom reference fetcher. Naming the binding
    // `ref_fetcher` matches the parameter names in
    // `load_and_run_chunk_with_retry` and
    // `prefetch_window_ref_bytes`; the codebase has multiple kinds
    // of fetchers (PSP, FASTA), so the `ref_` qualifier carries
    // meaning the bare `fetcher` did not.
    let ref_fetcher =
        StreamingChromRefFetcher::for_contig(fasta_path, &chrom.name).map_err(|source| {
            ChunkDriverError::OpenRefFetcher {
                contig: chrom.name.clone(),
                source,
            }
        })?;
    let masked_intervals: Vec<std::ops::Range<u32>> = if params.no_complexity_filter {
        Vec::new()
    } else {
        compute_dust_mask_for_chrom(&ref_fetcher, &chrom.name, chrom_length, &params.dust_cfg)?
    };
    // Mi9: borrow the per-chrom fetcher through the retry helper instead
    // of wrapping it in an `Arc<dyn …>`. The fetcher is owned by this
    // function, has exactly one consumer (`prefetch_window_ref_bytes`'s
    // sequential pre-pass), and never crosses the rayon `par_iter_mut`
    // boundary — the parallel block reads only the per-slot
    // `pre_fetched_ref_bytes`. The previous shape carried an Arc with
    // an `#[allow(clippy::arc_with_non_send_sync)]` justifying it.

    for c in carryover.iter_mut() {
        c.clear();
    }

    let max_group_span = params.grouper_cfg.max_variant_group_span;
    let chrom_one_past_end = chrom_length.saturating_add(1);
    let last_chunk_logical_extension = max_group_span.saturating_add(2);

    let mut chunk_range_start: u32 = 1;
    let mut psp_cursor: u32 = 1;
    // Re-used between load attempts when finalise_chunk_boundaries returns
    // `NoSafeGap` and we have to retry the chunk with a wider range —
    // the load consumes the driver-owned carryover, so we keep a
    // shadow copy before each load and restore it on retry.
    let mut carryover_snapshot: Vec<SampleColumns> = (0..carryover.len())
        .map(|_| SampleColumns::empty())
        .collect();

    // M12: bundle the per-iteration mutable state once so the three
    // phase helpers don't restate the borrow graph in their
    // signatures.
    let mut state = ChunkLoopState {
        chrom_id,
        chrom_length,
        chrom_one_past_end,
        last_chunk_logical_extension,
        max_group_span,
        masked_intervals: &masked_intervals,
        ref_fetcher: &ref_fetcher,
        params,
        psp_readers,
        writer,
        stats,
        chunk_scratch,
        fix_scratch,
        chunk,
        carryover,
        carryover_snapshot: &mut carryover_snapshot,
        worker_pool,
    };

    while chunk_range_start < chrom_one_past_end {
        let nominal_span = state.params.sizing.chunk_genomic_span.max(1);
        let max_span = nominal_span.saturating_mul(MAX_CHUNK_SPAN_GROWTH);

        let safe_end = load_and_run_chunk_with_retry(
            &mut state,
            chunk_range_start,
            psp_cursor,
            nominal_span,
            max_span,
        )?;

        if safe_end >= chrom_one_past_end {
            break;
        }
        chunk_range_start = safe_end;
        psp_cursor = state.chunk.range.end.min(chrom_one_past_end);
    }

    Ok(())
}

/// Hard cap on how much the chunk loader may grow a single chunk's
/// nominal span when [`finalise_chunk_boundaries`] cannot find a safe boundary —
/// the plan's "fail loudly on pathological input" backstop. The cap
/// is the largest multiplier of [`ChunkDriverParams::chunk_genomic_span`]
/// we will attempt before surfacing `NoSafeGap` to the caller.
const MAX_CHUNK_SPAN_GROWTH: u32 = 8;

/// M12: per-chrom mutable state bundle. Carries every long-lived
/// borrow that `load_and_run_chunk_with_retry` and its three phase
/// helpers need so the helpers' signatures stay short and the
/// borrow graph is visible at one place. Constructed in
/// `drive_one_chrom_generic` and reused across every chunk on the
/// chromosome.
///
/// Field grouping mirrors the helpers' phase responsibilities:
/// 1. Inputs (constants for the per-chrom loop): `chrom_id`,
///    `chrom_length`, `chrom_one_past_end`, `last_chunk_logical_extension`,
///    `max_group_span`, `masked_intervals`, `ref_fetcher`, `params`.
/// 2. Owned mutable state: `psp_readers`, `writer`, `stats`,
///    `chunk_scratch`, `fix_scratch`, `chunk`, `carryover`,
///    `carryover_snapshot`, `worker_pool`.
struct ChunkLoopState<'a, W: Read + Seek> {
    chrom_id: u32,
    chrom_length: u32,
    chrom_one_past_end: u32,
    last_chunk_logical_extension: u32,
    max_group_span: u32,
    masked_intervals: &'a [std::ops::Range<u32>],
    ref_fetcher: &'a dyn crate::fasta::fetcher::ChromRefFetcher,
    params: &'a ChunkDriverParams,
    psp_readers: &'a mut [PspReader<W>],
    writer: &'a mut CohortVcfWriter,
    stats: &'a mut ChunkDriverStats,
    chunk_scratch: &'a mut ChunkLoadScratch,
    fix_scratch: &'a mut BoundaryFinalisationScratch,
    chunk: &'a mut MaterialisedChunk,
    carryover: &'a mut [SampleColumns],
    carryover_snapshot: &'a mut [SampleColumns],
    worker_pool: &'a mut WorkerPool,
}

/// M12: orchestrate the three phases for one chunk — load with the
/// NoSafeGap retry, partition+prefetch per window, run the per-
/// window math in parallel + drain in order. Returns the chunk's
/// `safe_end` so the per-chrom loop can advance `chunk_range_start`.
///
/// On retry the carryover the previous chunk handed in is restored
/// from a snapshot the helper takes before the first load attempt;
/// the PSP iterators are re-created via `region_records`, which
/// re-seeks the block index per call.
fn load_and_run_chunk_with_retry<W>(
    state: &mut ChunkLoopState<W>,
    chunk_range_start: u32,
    psp_cursor: u32,
    nominal_span: u32,
    max_span: u32,
) -> Result<u32, ChunkDriverError>
where
    W: Read + Seek,
{
    load_chunk_with_safe_boundary_retry(
        state,
        chunk_range_start,
        psp_cursor,
        nominal_span,
        max_span,
    )?;
    run_chunk_windows_parallel(state)?;
    drain_window_outputs(state)?;
    Ok(state.chunk.safe_end)
}

/// M12 phase 1: snapshot the carryover, then loop over load + pre-pass
/// attempts doubling `attempt_span` on every `NoSafeGap`. Updates
/// `state.chunk` and `state.carryover` in place; bumps the
/// `chunks_loaded` / `chunk_variants_total` / `chunks_with_fewer_windows_than_requested`
/// counters as it goes.
fn load_chunk_with_safe_boundary_retry<W>(
    state: &mut ChunkLoopState<W>,
    chunk_range_start: u32,
    psp_cursor: u32,
    nominal_span: u32,
    max_span: u32,
) -> Result<(), ChunkDriverError>
where
    W: Read + Seek,
{
    debug_assert_eq!(state.carryover.len(), state.carryover_snapshot.len());
    // M14: snapshot the carryover so we can restore it on `NoSafeGap` —
    // `load_chunk_from_iters` drains it as part of the raw-load step
    // and we need the original contents back for the retry attempt.
    for (snap, carry) in state
        .carryover_snapshot
        .iter_mut()
        .zip(state.carryover.iter())
    {
        snap.clone_from_columns(carry);
    }

    let chrom_id = state.chrom_id;
    let chrom_length = state.chrom_length;
    let chrom_one_past_end = state.chrom_one_past_end;
    let last_chunk_logical_extension = state.last_chunk_logical_extension;
    let max_group_span = state.max_group_span;

    let mut attempt_span = nominal_span;
    loop {
        let nominal_end = chunk_range_start.saturating_add(attempt_span);
        let chunk_range_end = if nominal_end >= chrom_one_past_end {
            chrom_one_past_end.saturating_add(last_chunk_logical_extension)
        } else {
            nominal_end
        };

        let initial_load_span = chunk_range_end - chunk_range_start;
        // B3: bound the loader's `max_span` by this attempt's
        // `chunk_range_end` (which already reflects `attempt_span` +
        // the last-chunk extension). See the original B3 commit for
        // the decoupling rationale.
        let max_load_span = initial_load_span;
        let psp_inclusive_end = chunk_range_end.saturating_sub(1).min(chrom_length);

        let iters: Vec<_> = state
            .psp_readers
            .iter_mut()
            .map(|r| r.region_records(chrom_id, psp_cursor, psp_inclusive_end))
            .collect();

        let load_stats = load_chunk_from_iters(
            state.chunk_scratch,
            state.chunk,
            ChunkLoadExtent {
                chrom_id,
                range_start: chunk_range_start,
                initial_span: initial_load_span,
                target_variants: state
                    .params
                    .sizing
                    .target_variants_per_chunk
                    .map_or(0, NonZeroU32::get),
                max_span: max_load_span,
            },
            iters,
            state.carryover,
        )
        .map_err(|source| ChunkDriverError::LoadChunk { chrom_id, source })?;
        state.stats.chunks_loaded += 1;
        state.stats.chunk_variants_total += u64::from(load_stats.variant_count);

        // M7: `target_window_count` is `NonZeroUsize` — no more
        // `.max(1)` driver-side fallback that bypassed the pre-pass's
        // `ZeroTargetWindowCount` validation.
        let target_window_count = state.params.sizing.target_window_count.get();
        match finalise_chunk_boundaries(
            state.chunk,
            state.carryover,
            state.fix_scratch,
            max_group_span,
            target_window_count,
        ) {
            Ok(()) => {
                // Mi21: surface when finalise_chunk_boundaries silently produced
                // fewer windows than the operator asked for.
                if state.chunk.windows.len() < target_window_count {
                    state.stats.chunks_with_fewer_windows_than_requested += 1;
                }
                return Ok(());
            }
            Err(BoundaryFinalisationError::NoSafeGap { .. }) if attempt_span < max_span => {
                // M14: retry with double the span. Restore the carryover
                // the previous chunk handed in (the load drained it).
                for (carry, snap) in state
                    .carryover
                    .iter_mut()
                    .zip(state.carryover_snapshot.iter())
                {
                    carry.clone_from_columns(snap);
                }
                attempt_span = attempt_span
                    .checked_mul(2)
                    .unwrap_or(max_span)
                    .min(max_span);
                continue;
            }
            Err(source) => {
                return Err(ChunkDriverError::FinaliseChunkBoundaries { chrom_id, source });
            }
        }
    }
}

/// M12 phase 2: partition + REF prefetch per window (sequential),
/// then run the per-window math in parallel via rayon.
fn run_chunk_windows_parallel<W>(state: &mut ChunkLoopState<W>) -> Result<(), ChunkDriverError>
where
    W: Read + Seek,
{
    // Mi10: iterate `chunk.windows` directly — the loop body's only
    // `chunk` access is the immutable read `partition_window(chunk,
    // window, …)`, which is a shared reborrow under modern Rust's
    // NLL.
    let n_windows = state.chunk.windows.len();
    let n_samples = state.chunk.n_samples();
    state
        .worker_pool
        .ensure_capacity(n_windows.max(1), n_samples);

    let chrom_id = state.chrom_id;
    let max_group_span = state.max_group_span;
    let masked_intervals = state.masked_intervals;
    let ref_fetcher = state.ref_fetcher;
    // Split-borrow: `chunk` is read-only inside the loop while the
    // worker pool's slots are mutated; rust accepts the disjoint
    // field borrows here as long as we name them separately.
    let chunk: &MaterialisedChunk = state.chunk;
    let worker_pool = &mut *state.worker_pool;

    // ── Phase 1 (sequential): partition + prefetch REF bytes per
    //    window. The fetcher is `!Sync` so the prefetch can't run
    //    in parallel; the partition could, but it's cheap relative
    //    to the math and parallelising it would add rayon overhead
    //    for little gain. ──
    for (window_idx, slot) in worker_pool.slots.iter_mut().take(n_windows).enumerate() {
        let window = &chunk.windows[window_idx];
        partition_window(
            chunk,
            window,
            masked_intervals,
            max_group_span,
            &mut slot.partition_scratch,
            &mut slot.partition,
        )
        .map_err(|source| ChunkDriverError::Partition { chrom_id, source })?;
        prefetch_window_ref_bytes(
            &slot.partition,
            ref_fetcher,
            &mut slot.pre_fetched_ref_bytes,
        )
        .map_err(|source| ChunkDriverError::PrefetchRefBytes { chrom_id, source })?;
    }

    // ── Phase 2 (parallel): the per-window math (layers 1+2+3 + EM
    //    + posterior summarisation + QUAL). Each slot's
    //    `pre_fetched_ref_bytes` + `pipeline_scratch` + `posterior_records`
    //    + `stats` are disjoint, so the dispatch is structurally
    //    Send-safe. rayon's `par_iter_mut().try_for_each` collects
    //    the first error and propagates. ──
    let per_group_cfg = state.params.per_group_cfg;
    let posterior_cfg = &state.params.posterior_cfg;
    worker_pool.slots[..n_windows]
        .par_iter_mut()
        .try_for_each(|slot| {
            run_window(
                chunk,
                &slot.partition,
                &slot.pre_fetched_ref_bytes,
                per_group_cfg,
                posterior_cfg,
                &mut slot.pipeline_scratch,
                &mut slot.posterior_records,
                &mut slot.stats,
            )
        })
        .map_err(ChunkDriverError::RunWindow)?;
    Ok(())
}

/// M12 phase 3: sequential drain over the slots in window order —
/// preserves the prior per-record emit order even when phase 2 runs
/// the workers concurrently. M15: `slot.stats` lives next to
/// `slot.pipeline_scratch` as a separate `WindowRunStats` struct; the
/// drive takes its contents (by `std::mem::take`) so the slot's
/// counters reset for the next chunk while the per-driver
/// `ChunkDriverStats` totals accumulate.
fn drain_window_outputs<W>(state: &mut ChunkLoopState<W>) -> Result<(), ChunkDriverError>
where
    W: Read + Seek,
{
    let n_windows = state.chunk.windows.len();
    // Split-borrow: take the &mut on worker_pool but keep the
    // params/writer/stats refs for emit_or_drop.
    let worker_pool = &mut *state.worker_pool;
    let writer = &mut *state.writer;
    let stats = &mut *state.stats;
    let params = state.params;
    for slot in worker_pool.slots.iter_mut().take(n_windows) {
        let slot_stats = std::mem::take(&mut slot.stats);
        stats.lh_cap_groups_skipped += slot_stats.lh_cap_groups_skipped;
        stats.lh_cap_alleles_in_skipped += slot_stats.lh_cap_alleles_in_skipped;
        stats.groups_skipped_post_unify_ref_only += slot_stats.groups_skipped_post_unify_ref_only;
        for record in slot.posterior_records.drain(..) {
            emit_or_drop(record, params, writer, stats)?;
        }
    }
    Ok(())
}

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
                target_window_count: NonZeroUsize::MIN,
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
