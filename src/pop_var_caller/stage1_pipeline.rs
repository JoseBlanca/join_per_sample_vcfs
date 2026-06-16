//! Stage 1 pipeline helper — the BAQ → walker chain `run_pileup`
//! drives once per analysis region.
//!
//! Builds the BAQ → walker borrow chain over an already-opened reader
//! on its own stack and hands the walker into a caller-supplied
//! closure. The chain is self-referential by Rust's standards
//! (BaqStream borrows from the reader, the error-shedding adapter
//! borrows from BaqStream, the walker borrows from the adapter), so
//! the only way to expose it externally is via the callback pattern:
//! the helper owns every layer for the duration of the closure, then
//! snapshots counters after the closure returns and the borrow chain
//! unwinds.
//!
//! `run_pileup` opens one `AlignmentMergedReader::query` reader per
//! region and calls [`with_stage1_chain`] once per region, writing the
//! in-region columns to one shared PSP writer inside the closure.

use std::path::Path;

use noodles_fasta as fasta;

use crate::bam::alignment_input::{AlignmentMergedReader, FilterCounts, MappedRead};
use crate::bam::errors::AlignmentInputError;
use crate::baq::BaqConfig;
use crate::fasta::{ContigList, RepositoryRefFetcher};
use crate::pileup::per_sample::baq_stream::{BaqSkipCounts, BaqStream};
use crate::pileup::per_sample::read_pipeline::{
    RawPacket, ReorderReadIter, ResultPacket, StageCounts, produce_packets, run_worker,
};
use crate::pileup::per_sample::read_processor::ReadProcessingConfig;
use crate::pileup::walker::{self, PileupWalker, PreparedRead, WalkerConfig};

/// At or above this `--threads` count, Stage 1 runs the de-barriered
/// source → workers → reorder → walker pipeline; below it, the inline
/// bulk-synchronous `BaqStream` (whose coordination is cheaper when
/// there are too few cores for the pipeline's three roles to overlap).
///
/// `run_pileup` also reads this to decide whether to size the global
/// rayon pool: only the inline path uses rayon, so sizing it for the
/// pipeline path would spawn `n` idle rayon threads on top of the
/// pipeline's own `n` workers.
pub(crate) const STAGED_MIN_THREADS: usize = 4;

/// Bounded-channel depth per worker for the read pipeline — the same
/// back-pressure knob the cohort producer uses. Peak in-flight packets
/// ≈ `QUEUE_DEPTH_PER_WORKER × n_workers`.
const QUEUE_DEPTH_PER_WORKER: usize = 2;

use super::cli::PileupCliError;
use super::cli::error_bridge::ErrorSheddingAdapter;

/// Boxed walker type — uniform across the no-BAQ and BAQ-on branches.
///
/// The walker's input iterator is type-erased to a `Box<dyn ...>` so a
/// A source of `MappedRead`s for the Stage-1 pipeline: a coordinate-sorted
/// iterator plus the cheap-filter drop tally it accumulated. Implemented by
/// both the legacy [`AlignmentMergedReader`] (whole-file / `query`) and the
/// pooled [`SegmentMergedReads`](crate::bam::segment_merge::SegmentMergedReads),
/// so the pipeline is agnostic to which reader feeds it — the
/// SNP `--regions` retrofit swaps the reader without touching the BAQ /
/// walker plumbing.
///
/// `Send` is required because the staged topology moves the source into a
/// scoped producer thread.
pub trait MappedReadSource:
    Iterator<Item = Result<MappedRead, AlignmentInputError>> + Send
{
    /// The reads dropped by the cheap pre-decode filter so far (read after
    /// the source is drained). Named distinctly from any inherent
    /// `filter_counts` accessor so an inherent method returning
    /// `&FilterCounts` does not shadow this by-value trait method.
    fn filter_drop_counts(&self) -> FilterCounts;
}

impl MappedReadSource for AlignmentMergedReader {
    fn filter_drop_counts(&self) -> FilterCounts {
        *self.filter_counts()
    }
}

/// single closure type signature covers both branches. The boxing
/// costs one indirection per `PreparedRead` — invisible against the
/// HMM / walker bookkeeping the per-record budget already pays.
pub type Stage1Walker<'a> =
    PileupWalker<Box<dyn Iterator<Item = PreparedRead> + 'a>, &'a RepositoryRefFetcher>;

/// Context handed to the [`with_stage1_chain`] callback. Contains
/// the walker (by value — the closure may consume or drive it via
/// `by_ref`) plus borrowed metadata the caller typically needs for
/// header building.
pub struct Stage1PipelineContext<'a> {
    pub walker: Stage1Walker<'a>,
    pub sample_name: &'a str,
    pub contigs: &'a ContigList,
}

/// Counter snapshots collected after the closure returns. The fields
/// are owned values — safe to use past `with_stage1_chain`'s
/// return.
pub struct Stage1RunSummary {
    pub filter_counts: FilterCounts,
    /// `None` when the pipeline was built with `no_baq = true`.
    pub baq_skip_counts: Option<BaqSkipCounts>,
}

/// Bundle returned by [`with_stage1_chain`]. Carries the closure's
/// result `R`, the metadata the caller needs (sample name, contigs),
/// the counter snapshot, and any upstream error stashed by the
/// error-shedding adapter while the walker was running.
///
/// A `stashed_upstream_error` of `Some(_)` paired with `result =
/// Ok(_)` is the standard "walker exhausted cleanly but the source
/// stream had an error" case — the caller is expected to surface the
/// error and tear down whatever side effects the closure produced
/// (e.g. delete a half-written output file).
pub struct Stage1Outputs<R> {
    pub result: R,
    pub sample_name: String,
    pub contigs: ContigList,
    pub run_summary: Stage1RunSummary,
    pub stashed_upstream_error: Option<AlignmentInputError>,
}

/// Build the Stage 1 BAQ → walker chain over an already-opened
/// `reader` on this function's stack and hand the walker to `f`. After
/// `f` returns, snapshots `FilterCounts`, `BaqSkipCounts`, sample name
/// and contigs, plus any stashed upstream error, and returns them
/// bundled with `f`'s result.
///
/// The caller owns reader construction: `run_pileup` opens one
/// `query()` reader per analysis region and calls this once per region,
/// reusing the shared PSP writer inside `f`.
///
/// The `walker_fetcher` is built once by the caller and reused across
/// regions (see `run_pileup`); `reference` is still needed here to build
/// `BaqStream`'s per-worker fetchers on the BAQ-on path.
///
/// `E: From<PileupCliError>` lets an upstream read error stashed by the
/// error-shedding adapter lift into the closure's error type without
/// changing the caller's surface; `run_pileup`'s closure picks
/// `E = PileupCliError`.
#[allow(clippy::too_many_arguments)]
pub fn with_stage1_chain<Rd, R, E, F>(
    mut reader: Rd,
    reference: &Path,
    repository: &fasta::Repository,
    walker_fetcher: &RepositoryRefFetcher,
    baq_cfg: BaqConfig,
    proc_cfg: ReadProcessingConfig,
    walker_cfg: WalkerConfig,
    baq_chunk_size: usize,
    n_threads: usize,
    no_baq: bool,
    sample_name: &str,
    contigs: &ContigList,
    f: F,
) -> Result<Stage1Outputs<R>, E>
where
    Rd: MappedReadSource,
    F: FnOnce(Stage1PipelineContext<'_>) -> Result<R, E>,
    E: From<PileupCliError>,
{
    // The reader is already opened and validated; its sample name and
    // contig list are passed in by the caller (they are the same
    // cross-validated metadata the reader was built from), so the
    // pipeline does not depend on the reader carrying accessors for them.

    // Reference fetcher for the walker — [`RepositoryRefFetcher`] reads
    // the walker's reference windows from the same shared noodles
    // `fasta::Repository` the reader keeps resident (for CRAM decode).
    // It is **injected** by the caller (built once over the shared
    // repository). The read-processing stage builds its own per-worker
    // reference access: an uppercased-window `ManualEvictChromRefFetcher`
    // for BAQ and a raw-byte `RawContigRefCache` (over `repository`) for
    // F3/F1.

    // Run the read-processing stage (G2/F3/F1 + BAQ; passthrough under
    // `--no-baq`) feeding the serial walker, by one of two topologies
    // chosen on the thread budget. Both produce the same `(closure
    // result, stage counts, stashed upstream error, reader-level filter
    // counts)`, so the snapshot tail below is shared and the choice
    // cannot change results.
    let apply_baq = !no_baq;
    let (result, stage_counts, stashed_upstream_error, reader_filter_counts): (
        Result<R, E>,
        StageCounts,
        Option<AlignmentInputError>,
        FilterCounts,
    ) = if n_threads >= STAGED_MIN_THREADS {
        run_pipelined(
            reader,
            reference,
            repository,
            walker_fetcher,
            baq_cfg,
            proc_cfg,
            &walker_cfg,
            baq_chunk_size,
            apply_baq,
            n_threads,
            sample_name,
            contigs,
            f,
        )
    } else {
        // Inline bulk-synchronous path: reader → BaqStream → adapter →
        // walker, all driven on this thread.
        let mut read_stream = BaqStream::new(
            reader.by_ref(),
            baq_cfg,
            proc_cfg,
            apply_baq,
            reference.to_path_buf(),
            repository.clone(),
            contigs.clone(),
            baq_chunk_size,
        );
        let mut adapter = ErrorSheddingAdapter::new(read_stream.by_ref());
        let error_handle = adapter.error_handle();
        let input: Box<dyn Iterator<Item = PreparedRead> + '_> = Box::new(adapter.by_ref());
        let walker = walker::run(input, walker_fetcher, &walker_cfg);
        let ctx = Stage1PipelineContext {
            walker,
            sample_name,
            contigs,
        };
        let result = f(ctx);
        drop(adapter);
        let stash = error_handle.take();
        let stage_counts = StageCounts {
            baq_skips: *read_stream.skip_counts(),
            bad_cigar: read_stream.bad_cigar_count(),
            high_mismatch: read_stream.high_mismatch_count(),
        };
        drop(read_stream);
        let reader_filter_counts = reader.filter_drop_counts();
        (result, stage_counts, stash, reader_filter_counts)
    };

    // The reader no longer counts G2/F1 (those filters moved to the
    // read-processing stage); fold the stage's tallies back in so the
    // run summary's totals match the old serial path exactly. BAQ skip
    // counts only exist when BAQ ran.
    let mut filter_counts = reader_filter_counts;
    filter_counts.bad_cigar += stage_counts.bad_cigar;
    filter_counts.high_mismatch_fraction += stage_counts.high_mismatch;
    let baq_skip_counts: Option<BaqSkipCounts> = apply_baq.then_some(stage_counts.baq_skips);

    // Surface order: if the closure errored AND the source stream had
    // already failed (the `ErrorSheddingAdapter` translated it to
    // end-of-stream), prefer the **upstream** error — it is the root
    // cause; the closure's failure is the downstream symptom of the
    // truncated stream. Logic factored into `prefer_upstream_or_closure`
    // so it can be unit-tested in isolation. On Ok the stash is
    // returned alongside the result so the caller can surface it
    // separately (and tear down any half-written output).
    let (r, stashed_upstream_error) = prefer_upstream_or_closure(result, stashed_upstream_error)?;

    Ok(Stage1Outputs {
        result: r,
        sample_name: sample_name.to_string(),
        contigs: contigs.clone(),
        run_summary: Stage1RunSummary {
            filter_counts,
            baq_skip_counts,
        },
        stashed_upstream_error,
    })
}

/// De-barriered Stage 1 topology: a source thread (owns the reader,
/// packetises coordinate-sorted reads), `n_threads` worker threads (the
/// full per-read fold), and the walker on this thread fed by a
/// sequence-ordered reorder buffer — all wired through bounded
/// `crossbeam` channels so the stages overlap. Returns the closure
/// result, the stage's filter tallies, any stashed upstream read error,
/// and the reader-level filter counts (recovered from the producer
/// thread). Byte-identical to the inline path: the reorder buffer
/// restores global coordinate order before the (serial) walker.
#[allow(clippy::too_many_arguments)]
fn run_pipelined<Rd, R, E, F>(
    mut reader: Rd,
    reference: &Path,
    repository: &fasta::Repository,
    walker_fetcher: &RepositoryRefFetcher,
    baq_cfg: BaqConfig,
    proc_cfg: ReadProcessingConfig,
    walker_cfg: &WalkerConfig,
    baq_chunk_size: usize,
    apply_baq: bool,
    n_threads: usize,
    sample_name: &str,
    contigs: &ContigList,
    f: F,
) -> (
    Result<R, E>,
    StageCounts,
    Option<AlignmentInputError>,
    FilterCounts,
)
where
    Rd: MappedReadSource,
    F: FnOnce(Stage1PipelineContext<'_>) -> Result<R, E>,
{
    let err_cell = std::sync::Mutex::new(None);
    let n_workers = n_threads.max(1);
    let cap = (QUEUE_DEPTH_PER_WORKER * n_workers).max(1);
    let (pkt_tx, pkt_rx) = crossbeam_channel::bounded::<RawPacket>(cap);
    let (res_tx, res_rx) = crossbeam_channel::bounded::<ResultPacket>(cap);

    let err_ref = &err_cell;

    let (result, stage_counts, reader_filter_counts) = std::thread::scope(|scope| {
        // Source thread owns the reader; hands back its final filter
        // counts (the reader is moved in, so the caller can't read them).
        let producer = scope.spawn(move || {
            produce_packets(&mut reader, pkt_tx, err_ref, baq_chunk_size);
            reader.filter_drop_counts()
        });

        // Worker pool: each pulls packets and runs G2/F3/F1 + BAQ.
        for _ in 0..n_workers {
            let prx = pkt_rx.clone();
            let rtx = res_tx.clone();
            scope.spawn(move || {
                run_worker(
                    prx, rtx, baq_cfg, proc_cfg, apply_baq, repository, reference, contigs,
                );
            });
        }
        // Drop the originals so each channel disconnects once the
        // producer / all workers finish.
        drop(pkt_rx);
        drop(res_tx);

        // Consumer on this thread: reorder → walker → writer (via `f`).
        let mut reorder = ReorderReadIter::new(res_rx);
        let input: Box<dyn Iterator<Item = PreparedRead> + '_> = Box::new(reorder.by_ref());
        let walker = walker::run(input, walker_fetcher, walker_cfg);
        let ctx = Stage1PipelineContext {
            walker,
            sample_name,
            contigs,
        };
        let result = f(ctx);
        // Walker dropped → `reorder` borrow released. Snapshot the
        // counts, then drop `reorder` (releasing `res_rx`) *before*
        // joining the producer: if `f` errored without draining, this
        // disconnects the channels so the workers and producer unblock
        // from full sends instead of deadlocking the join.
        let stage_counts = reorder.counts();
        drop(reorder);
        let reader_filter_counts = producer
            .join()
            .expect("read-pipeline producer thread panicked");
        (result, stage_counts, reader_filter_counts)
    });

    let stashed_upstream_error = err_cell
        .into_inner()
        .expect("read-pipeline err_cell poisoned");
    (
        result,
        stage_counts,
        stashed_upstream_error,
        reader_filter_counts,
    )
}

/// Surface-order rule used by [`with_stage1_chain`]: if both the
/// closure failed *and* an upstream CRAM-input error was stashed by
/// the `ErrorSheddingAdapter`, prefer the upstream error — it is the
/// root cause. On Ok the stash is handed back to the caller for
/// separate surfacing. Pure function so it is unit-testable in
/// isolation.
fn prefer_upstream_or_closure<R, E>(
    result: Result<R, E>,
    stashed: Option<AlignmentInputError>,
) -> Result<(R, Option<AlignmentInputError>), E>
where
    E: From<PileupCliError>,
{
    match result {
        Ok(v) => Ok((v, stashed)),
        Err(closure_err) => match stashed {
            Some(upstream) => Err(PileupCliError::AlignmentInput(upstream).into()),
            None => Err(closure_err),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::errors::AlignmentInputError;

    /// On success, the closure result is preserved and the stash (if
    /// any) is handed back to the caller.
    #[test]
    fn prefer_upstream_or_closure_returns_ok_with_stash_when_closure_succeeded() {
        let r: Result<(u8, Option<AlignmentInputError>), PileupCliError> =
            prefer_upstream_or_closure(Ok(42), None);
        assert!(matches!(r, Ok((42, None))));

        let stashed = AlignmentInputError::Io {
            path: std::path::PathBuf::from("sample.cram"),
            source: std::io::Error::other("cram truncated"),
        };
        let r: Result<(u8, Option<AlignmentInputError>), PileupCliError> =
            prefer_upstream_or_closure(Ok(7), Some(stashed));
        let (v, stash) = r.expect("Ok branch must propagate");
        assert_eq!(v, 7);
        assert!(matches!(stash, Some(AlignmentInputError::Io { .. })));
    }

    /// When the closure errors and the stash is empty, the closure
    /// error propagates unchanged.
    #[test]
    fn prefer_upstream_or_closure_returns_closure_err_when_stash_is_none() {
        let closure_err = PileupCliError::TimestampFormat;
        let r: Result<(u8, Option<AlignmentInputError>), PileupCliError> =
            prefer_upstream_or_closure::<u8, _>(Err(closure_err), None);
        assert!(matches!(r, Err(PileupCliError::TimestampFormat)));
    }

    /// When the closure errors AND the stash holds a CRAM-input error,
    /// the upstream error wins — this is the M1 case (the stash was
    /// silently dropped before the fix).
    #[test]
    fn prefer_upstream_or_closure_prefers_stash_when_both_present() {
        let closure_err = PileupCliError::TimestampFormat;
        let stashed = AlignmentInputError::Io {
            path: std::path::PathBuf::from("sample.cram"),
            source: std::io::Error::other("cram truncated"),
        };
        let r: Result<(u8, Option<AlignmentInputError>), PileupCliError> =
            prefer_upstream_or_closure::<u8, _>(Err(closure_err), Some(stashed));
        match r {
            Err(PileupCliError::AlignmentInput(AlignmentInputError::Io { .. })) => (),
            other => panic!("expected upstream AlignmentInput::Io, got {other:?}"),
        }
    }
}
