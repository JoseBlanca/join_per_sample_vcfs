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

use crate::bam::alignment_input::{AlignmentMergedReader, FilterCounts};
use crate::bam::errors::AlignmentInputError;
use crate::baq::BaqConfig;
use crate::fasta::{ContigList, RepositoryRefFetcher};
use crate::pileup::per_sample::baq_stream::{BaqSkipCounts, BaqStream};
use crate::pileup::per_sample::read_processor::ReadProcessingConfig;
use crate::pileup::walker::{self, PileupWalker, PreparedRead, WalkerConfig};

use super::cli::PileupCliError;
use super::cli::error_bridge::ErrorSheddingAdapter;

/// Boxed walker type — uniform across the no-BAQ and BAQ-on branches.
///
/// The walker's input iterator is type-erased to a `Box<dyn ...>` so a
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
pub fn with_stage1_chain<R, E, F>(
    mut reader: AlignmentMergedReader,
    reference: &Path,
    repository: &fasta::Repository,
    walker_fetcher: &RepositoryRefFetcher,
    baq_cfg: BaqConfig,
    proc_cfg: ReadProcessingConfig,
    walker_cfg: WalkerConfig,
    baq_chunk_size: usize,
    no_baq: bool,
    f: F,
) -> Result<Stage1Outputs<R>, E>
where
    F: FnOnce(Stage1PipelineContext<'_>) -> Result<R, E>,
    E: From<PileupCliError>,
{
    // The reader is already opened and validated (via `new()` for a
    // whole-file stream, or `query()` for one region). Capture the
    // identification metadata it carries.
    let sample_name = reader.sample_name().to_string();
    let contigs = reader.contigs().clone();

    // 2. Reference fetcher for the walker — [`RepositoryRefFetcher`]
    //    reads the walker's reference windows from the same shared
    //    noodles `fasta::Repository` the reader already keeps resident
    //    (for CRAM decode, on both CRAM and BAM). It is **injected** by
    //    the caller (built once over the shared repository) — the walker
    //    no longer opens a second, independent streaming reader over the
    //    same FASTA. Random access over a resident contig, so there is
    //    no monotonic-forward contract to satisfy across region/contig
    //    boundaries.
    //
    //    The read-processing stage (`BaqStream`) builds its own
    //    per-worker reference access: an uppercased-window
    //    `ManualEvictChromRefFetcher` for BAQ and a raw-byte
    //    `RawContigRefCache` (over `repository`) for F3/F1.

    // 3. Build the read-processing stream → error-shedding adapter →
    //    walker chain, drive the closure, then snapshot the counters.
    //    One path covers both BAQ and `--no-baq`: the stream applies
    //    BAQ only when `apply_baq` is set, but runs G2/F3/F1 either way
    //    (F3 left-alignment is mandatory even under `--no-baq`).
    let mut read_stream = BaqStream::new(
        reader.by_ref(),
        baq_cfg,
        proc_cfg,
        !no_baq,
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
        sample_name: &sample_name,
        contigs: &contigs,
    };
    let result = f(ctx);
    // Closure dropped walker → `&mut adapter` released. Now drop adapter
    // → `&mut read_stream` released. Then snapshot.
    drop(adapter);
    let stashed_upstream_error = error_handle.take();
    // BAQ skip counts only exist when BAQ ran; G2/F1 rejections are
    // tallied in the stream in both modes and folded into FilterCounts.
    let baq_skip_counts: Option<BaqSkipCounts> = (!no_baq).then(|| *read_stream.skip_counts());
    let stage_bad_cigar = read_stream.bad_cigar_count();
    let stage_high_mismatch = read_stream.high_mismatch_count();
    drop(read_stream);

    // The reader no longer counts G2/F1 (those filters moved to the
    // read-processing stage); fold the stage's tallies back in so the
    // run summary's totals match the old serial path exactly.
    let mut filter_counts = *reader.filter_counts();
    filter_counts.bad_cigar += stage_bad_cigar;
    filter_counts.high_mismatch_fraction += stage_high_mismatch;

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
        sample_name,
        contigs,
        run_summary: Stage1RunSummary {
            filter_counts,
            baq_skip_counts,
        },
        stashed_upstream_error,
    })
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
