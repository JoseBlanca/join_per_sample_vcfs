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
//! (`ChunkLoadScratch`, `FixBoundariesScratch`, `PartitionScratch`,
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
use std::path::{Path, PathBuf};
use thiserror::Error;

use crate::fasta::StreamingChromRefFetcher;
use crate::fasta::fetcher::ChromRefFetchError;
use crate::psp::header::ParsedChromosome;
use crate::psp::{PspReadError, PspReader};
use crate::var_calling::cohort_block::columns::{MaterialisedChunk, SampleColumns};
use crate::var_calling::cohort_block::loader::{
    ChunkLoadError, ChunkLoadScratch, load_chunk_from_iters,
};
use crate::var_calling::cohort_block::partition::{PartitionError, partition_window};
use crate::var_calling::cohort_block::pre_pass::{
    FixBoundariesError, FixBoundariesScratch, fix_boundaries,
};
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

/// Per-driver tuning bundle. Carries every per-stage config the
/// chunk loop needs.
#[non_exhaustive]
pub struct ChunkDriverParams {
    pub no_complexity_filter: bool,
    pub dust_cfg: DustFilterConfig,
    pub grouper_cfg: GrouperConfig,
    pub per_group_cfg: PerGroupMergerConfig,
    pub posterior_cfg: PosteriorEngineConfig,
    pub min_qual_phred: f64,
    pub min_alt_obs_per_sample: u32,
    pub no_mapq_diff_filter: bool,
    pub min_mapq_diff_t: f32,
    /// Nominal BP span of the loader's first pull attempt per chunk.
    /// Acts as a starting size; the loader grows up to
    /// `chunk_genomic_span × MAX_CHUNK_SPAN_GROWTH` when
    /// `target_variants_per_chunk` is not satisfied or the pre-pass
    /// can't find a safe boundary.
    pub chunk_genomic_span: u32,
    /// Soft lower bound on the post-filter variant count per chunk.
    /// `0` disables the variant-bounded extension (each chunk is one
    /// pull attempt of `chunk_genomic_span` BP). Non-zero values let
    /// the loader grow each chunk's span adaptively until the
    /// kept-position count crosses this target — decoupling worker
    /// workload from PSP block size + variant density. See
    /// [the Phase B prereq plan](https://github.com/JoseBlanca/join_per_sample_vcfs/blob/main/doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_b1_variant_bounded_chunks.md).
    pub target_variants_per_chunk: u32,
    /// Number of parallel worker windows per chunk.
    /// [`fix_boundaries`] places `target_window_count - 1` internal
    /// boundaries inside each chunk's `[range.start, safe_end)` so
    /// the [`WorkerPool`] can dispatch the per-window math
    /// concurrently. `1` (default) preserves the
    /// sequential single-window-per-chunk behaviour byte-for-byte.
    /// See [the Phase B plan](https://github.com/JoseBlanca/join_per_sample_vcfs/blob/main/doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_b_parallel_windows.md).
    pub target_window_count: usize,
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
/// - **Per-chunk pipeline**: `LoadChunk`, `FixBoundaries`, `Partition`
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
    FixBoundaries {
        chrom_id: u32,
        #[source]
        source: FixBoundariesError,
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
    let mut fix_scratch = FixBoundariesScratch::new();
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
    if !passes_min_alt_obs(&record, params.min_alt_obs_per_sample) {
        stats.records_dropped_low_alt_obs += 1;
        return Ok(());
    }
    if !record.is_variant_call() {
        stats.records_dropped_hom_ref += 1;
        return Ok(());
    }
    if record.qual_phred < params.min_qual_phred {
        stats.records_dropped_low_qual += 1;
        return Ok(());
    }
    if !params.no_mapq_diff_filter && record_fails_mapq_diff_t(&record, params.min_mapq_diff_t) {
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
    let base_iter = ref_fetcher
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
    fix_scratch: &mut FixBoundariesScratch,
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
    // Re-used between load attempts when fix_boundaries returns
    // `NoSafeGap` and we have to retry the chunk with a wider range —
    // the load consumes the driver-owned carryover, so we keep a
    // shadow copy before each load and restore it on retry.
    let mut carryover_snapshot: Vec<SampleColumns> = (0..carryover.len())
        .map(|_| SampleColumns::empty())
        .collect();

    while chunk_range_start < chrom_one_past_end {
        let nominal_span = params.chunk_genomic_span.max(1);
        let max_span = nominal_span.saturating_mul(MAX_CHUNK_SPAN_GROWTH);

        let safe_end = load_and_run_chunk_with_retry(
            chrom_id,
            chrom_length,
            chrom_one_past_end,
            last_chunk_logical_extension,
            chunk_range_start,
            psp_cursor,
            nominal_span,
            max_span,
            psp_readers,
            &masked_intervals,
            &ref_fetcher,
            params,
            writer,
            stats,
            chunk_scratch,
            fix_scratch,
            chunk,
            carryover,
            &mut carryover_snapshot,
            worker_pool,
            max_group_span,
        )?;

        if safe_end >= chrom_one_past_end {
            break;
        }
        chunk_range_start = safe_end;
        psp_cursor = chunk.range.end.min(chrom_one_past_end);
    }

    Ok(())
}

/// Hard cap on how much the chunk loader may grow a single chunk's
/// nominal span when [`fix_boundaries`] cannot find a safe boundary —
/// the plan's "fail loudly on pathological input" backstop. The cap
/// is the largest multiplier of [`ChunkDriverParams::chunk_genomic_span`]
/// we will attempt before surfacing `NoSafeGap` to the caller.
const MAX_CHUNK_SPAN_GROWTH: u32 = 8;

/// Load + pre-pass one chunk, retrying with a 2× wider load range if
/// the pre-pass can't find a safe boundary inside the current range.
/// Returns the chunk's `safe_end` on success — the caller advances
/// `chunk_range_start` to that value.
///
/// On retry the carryover the previous chunk handed in is restored
/// from a snapshot the helper takes before the first load attempt; the
/// PSP iterators are re-created via `region_records`, which re-seeks
/// the block index per call.
#[allow(clippy::too_many_arguments)]
fn load_and_run_chunk_with_retry<W>(
    chrom_id: u32,
    chrom_length: u32,
    chrom_one_past_end: u32,
    last_chunk_logical_extension: u32,
    chunk_range_start: u32,
    psp_cursor: u32,
    nominal_span: u32,
    max_span: u32,
    psp_readers: &mut [PspReader<W>],
    masked_intervals: &[std::ops::Range<u32>],
    ref_fetcher: &dyn crate::fasta::fetcher::ChromRefFetcher,
    params: &ChunkDriverParams,
    writer: &mut CohortVcfWriter,
    stats: &mut ChunkDriverStats,
    chunk_scratch: &mut ChunkLoadScratch,
    fix_scratch: &mut FixBoundariesScratch,
    chunk: &mut MaterialisedChunk,
    carryover: &mut [SampleColumns],
    carryover_snapshot: &mut [SampleColumns],
    worker_pool: &mut WorkerPool,
    max_group_span: u32,
) -> Result<u32, ChunkDriverError>
where
    W: Read + Seek,
{
    debug_assert_eq!(carryover.len(), carryover_snapshot.len());
    // M14: snapshot the carryover so we can restore it on `NoSafeGap` —
    // `load_chunk_from_iters` drains it as part of the raw-load step
    // and we need the original contents back for the retry attempt.
    for (snap, carry) in carryover_snapshot.iter_mut().zip(carryover.iter()) {
        snap.clone_from_columns(carry);
    }

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
        // the last-chunk extension). The previous shape passed
        // `extension_cap_end - chunk_range_start` (chrom-wide cap)
        // unconditionally, so the loader's internal variant-bounded
        // extension already ran to chrom-end on the *first* attempt,
        // making NoSafeGap retries produce identical chunks — the
        // outer retry's `attempt_span` doubling was a no-op when
        // `target_variants_per_chunk > 0`. Tying `max_load_span` to
        // `chunk_range_end` decouples the two axes: the loader's
        // variant-bounded extension grows within `attempt_span`, the
        // outer NoSafeGap retry grows `attempt_span` itself. PSP
        // iterators are re-opened per attempt via `region_records`
        // and respect the same `chunk_range_end` cap.
        let max_load_span = initial_load_span;
        let psp_inclusive_end = chunk_range_end.saturating_sub(1).min(chrom_length);

        let iters: Vec<_> = psp_readers
            .iter_mut()
            .map(|r| r.region_records(chrom_id, psp_cursor, psp_inclusive_end))
            .collect();

        let load_stats = load_chunk_from_iters(
            chunk_scratch,
            chunk,
            chrom_id,
            chunk_range_start,
            initial_load_span,
            params.target_variants_per_chunk,
            max_load_span,
            iters,
            carryover,
        )
        .map_err(|source| ChunkDriverError::LoadChunk { chrom_id, source })?;
        stats.chunks_loaded += 1;
        stats.chunk_variants_total += u64::from(load_stats.variant_count);

        match fix_boundaries(
            chunk,
            carryover,
            fix_scratch,
            max_group_span,
            params.target_window_count.max(1),
        ) {
            Ok(()) => break,
            Err(FixBoundariesError::NoSafeGap { .. }) if attempt_span < max_span => {
                // M14: retry with double the span. Restore the carryover
                // the previous chunk handed in (the load drained it).
                for (carry, snap) in carryover.iter_mut().zip(carryover_snapshot.iter()) {
                    carry.clone_from_columns(snap);
                }
                attempt_span = attempt_span
                    .checked_mul(2)
                    .unwrap_or(max_span)
                    .min(max_span);
                continue;
            }
            Err(source) => {
                return Err(ChunkDriverError::FixBoundaries { chrom_id, source });
            }
        }
    }

    // Mi10: iterate `chunk.windows` directly — the loop body's only
    // `chunk` access is the immutable read `partition_window(chunk,
    // window, …)`, which is a shared reborrow under modern Rust's
    // NLL. The previous `chunk.windows.clone()` was a defensive copy
    // that allocated a fresh `Vec<Range<u32>>` per chunk for no
    // borrow-checker reason.
    let n_windows = chunk.windows.len();
    let n_samples = chunk.n_samples();
    worker_pool.ensure_capacity(n_windows.max(1), n_samples);

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
    //    `pre_fetched_ref_bytes` + `scratch` + `output_buf` are
    //    disjoint, so the dispatch is structurally
    //    Send-safe. rayon's `par_iter_mut().try_for_each` collects
    //    the first error and propagates. ──
    let per_group_cfg = params.per_group_cfg;
    let posterior_cfg = &params.posterior_cfg;
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
            )
        })
        .map_err(ChunkDriverError::RunWindow)?;

    // ── Sequential drain in window order — preserves the prior
    //    per-record emit order even when step 4 runs the workers
    //    concurrently. ──
    for slot in worker_pool.slots.iter_mut().take(n_windows) {
        let (cap_g, cap_a) = slot.pipeline_scratch.take_lh_cap_stats();
        stats.lh_cap_groups_skipped += cap_g;
        stats.lh_cap_alleles_in_skipped += cap_a;
        stats.groups_skipped_post_unify_ref_only +=
            slot.pipeline_scratch.take_post_unify_ref_only_count();
        for record in slot.posterior_records.drain(..) {
            emit_or_drop(record, params, writer, stats)?;
        }
    }

    Ok(chunk.safe_end)
}
