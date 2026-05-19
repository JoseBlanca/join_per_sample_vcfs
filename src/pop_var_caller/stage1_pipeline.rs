//! Stage 1 pipeline helper — shared core of `run_pileup` and (when
//! Task 8 lands) `run_var_calling_from_bam`.
//!
//! Builds the CRAM → BAQ → walker borrow chain on its own stack and
//! hands the walker into a caller-supplied closure. The chain is
//! self-referential by Rust's standards (BaqStream borrows from the
//! reader, the error-shedding adapter borrows from BaqStream, the
//! walker borrows from the adapter), so the only way to expose it
//! externally is via the callback pattern: the helper owns every
//! layer for the duration of the closure, then snapshots counters
//! after the closure returns and the borrow chain unwinds.
//!
//! No temp `.psp` is ever written — the var-calling-from-bam path
//! consumes the walker's record stream directly. See
//! [feedback_no_silent_intermediates] in the user's memory.

use std::path::{Path, PathBuf};

use crate::per_sample_pileup::baq::{BaqConfig, BaqSkipCounts, BaqStream, prepare_passthrough};
use crate::per_sample_pileup::cram_input::{
    ContigList, CramMergedReader, CramMergedReaderConfig, FilterCounts,
};
use crate::per_sample_pileup::errors::CramInputError;
use crate::per_sample_pileup::pileup::{self, PileupWalker, PreparedRead, WalkerConfig};
use crate::per_sample_pileup::ref_fetcher::{ChromBoundaryRefFetcher, SyncRefFetcher};

use super::cli::PileupCliError;
use super::cli::error_bridge::ErrorSheddingAdapter;

/// Boxed walker type — uniform across the no-BAQ and BAQ-on branches.
///
/// The walker's input iterator is type-erased to a `Box<dyn ...>` so a
/// single closure type signature covers both branches. The boxing
/// costs one indirection per `PreparedRead` — invisible against the
/// HMM / walker bookkeeping the per-record budget already pays.
pub type Stage1Walker<'a> =
    PileupWalker<Box<dyn Iterator<Item = PreparedRead> + 'a>, &'a ChromBoundaryRefFetcher>;

/// Context handed to the [`with_stage1_pipeline`] callback. Contains
/// the walker (by value — the closure may consume or drive it via
/// `by_ref`) plus borrowed metadata the caller typically needs for
/// header building.
pub struct Stage1PipelineContext<'a> {
    pub walker: Stage1Walker<'a>,
    pub sample_name: &'a str,
    pub contigs: &'a ContigList,
}

/// Counter snapshots collected after the closure returns. The fields
/// are owned values — safe to use past `with_stage1_pipeline`'s
/// return.
pub struct Stage1RunSummary {
    pub filter_counts: FilterCounts,
    /// `None` when the pipeline was built with `no_baq = true`.
    pub baq_skip_counts: Option<BaqSkipCounts>,
}

/// Bundle returned by [`with_stage1_pipeline`]. Carries the closure's
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
    pub stashed_upstream_error: Option<CramInputError>,
}

/// Build the Stage 1 pipeline on this function's stack and hand the
/// walker to `f`. After `f` returns, snapshots `FilterCounts`,
/// `BaqSkipCounts`, sample name and contigs, plus any stashed
/// upstream error, and returns them bundled with `f`'s result.
///
/// `E: From<PileupCliError>` lets internal setup errors (reader open,
/// fetcher init) lift into the closure's error type without changing
/// the caller's surface. For `run_pileup` the closure picks
/// `E = PileupCliError` (identity); for `run_var_calling_from_bam`,
/// a manual `From<PileupCliError>` on the from-bam error provides
/// the bridge.
#[allow(clippy::too_many_arguments)]
pub fn with_stage1_pipeline<R, E, F>(
    crams: &[PathBuf],
    reference: &Path,
    cram_cfg: CramMergedReaderConfig,
    baq_cfg: BaqConfig,
    walker_cfg: WalkerConfig,
    baq_chunk_size: usize,
    no_baq: bool,
    f: F,
) -> Result<Stage1Outputs<R>, E>
where
    F: FnOnce(Stage1PipelineContext<'_>) -> Result<R, E>,
    E: From<PileupCliError>,
{
    // 1. Open the merged reader and capture identification metadata.
    let mut reader: CramMergedReader =
        CramMergedReader::new(crams, reference, cram_cfg).map_err(PileupCliError::CramInput)?;
    let sample_name = reader.sample_name().to_string();
    let contigs = reader.contigs().clone();

    // 2. Reference fetchers — BAQ needs Send+Sync, walker is
    //    single-threaded and uses the chrom-boundary-evicting one.
    let baq_fetcher =
        SyncRefFetcher::new(reference, contigs.clone()).map_err(PileupCliError::Io)?;
    let walker_fetcher =
        ChromBoundaryRefFetcher::new(reference, contigs.clone()).map_err(PileupCliError::Io)?;

    // 3. Build the boxed input iterator + handle for upstream errors.
    //    Both branches converge on `Box<dyn Iterator<Item =
    //    PreparedRead>>` so the walker has the same concrete type
    //    either way (see `Stage1Walker`).
    let result: Result<R, E>;
    let stashed_upstream_error: Option<CramInputError>;
    let baq_skip_counts: Option<BaqSkipCounts>;

    if no_baq {
        // No-BAQ: passthrough map directly off the reader.
        let mut adapter = ErrorSheddingAdapter::new(reader.by_ref().map(|r| {
            r.map(|read| {
                let chrom_id = u32::try_from(read.ref_id).expect("ref_id fits u32");
                prepare_passthrough(read, chrom_id)
            })
        }));
        let error_handle = adapter.error_handle();
        let input: Box<dyn Iterator<Item = PreparedRead> + '_> = Box::new(adapter.by_ref());
        let walker = pileup::run(input, &walker_fetcher, &walker_cfg);
        let ctx = Stage1PipelineContext {
            walker,
            sample_name: &sample_name,
            contigs: &contigs,
        };
        result = f(ctx);
        // Closure dropped ctx (and walker inside it) — `&mut adapter`
        // borrow released. Now drop adapter to release `&mut reader`.
        drop(adapter);
        stashed_upstream_error = error_handle.take();
        baq_skip_counts = None;
    } else {
        // BAQ-on: reader → BaqStream → adapter → walker.
        let mut baq_stream = BaqStream::new(reader.by_ref(), baq_cfg, &baq_fetcher, baq_chunk_size);
        let mut adapter = ErrorSheddingAdapter::new(baq_stream.by_ref());
        let error_handle = adapter.error_handle();
        let input: Box<dyn Iterator<Item = PreparedRead> + '_> = Box::new(adapter.by_ref());
        let walker = pileup::run(input, &walker_fetcher, &walker_cfg);
        let ctx = Stage1PipelineContext {
            walker,
            sample_name: &sample_name,
            contigs: &contigs,
        };
        result = f(ctx);
        // Closure dropped walker → `&mut adapter` released. Now drop
        // adapter → `&mut baq_stream` released. Then snapshot.
        drop(adapter);
        stashed_upstream_error = error_handle.take();
        baq_skip_counts = Some(*baq_stream.skip_counts());
        drop(baq_stream);
    }

    let filter_counts = *reader.filter_counts();
    let r = result?;

    Ok(Stage1Outputs {
        result: r,
        sample_name,
        contigs,
        run_summary: Stage1RunSummary {
            filter_counts,
            baq_skip_counts,
        },
        stashed_upstream_error,
    })
}
