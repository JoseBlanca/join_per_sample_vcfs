//! Shared cohort-pipeline driver — DUST → grouper → per-group
//! merger → posterior engine → VCF writer, plus the tmp-rename +
//! best-effort cleanup discipline. Closes **M11** from the
//! 2026-05-19 cohort CLI review: both `run_var_calling` (multi-
//! sample, `.psp` inputs) and `run_var_calling_from_bam` (single-
//! sample, walker input) previously open-coded this wiring twice.
//!
//! The driver takes a pre-built `PerPositionMerger` upstream
//! because the merger's construction differs between the two
//! subcommands (multi-vs-single, .psp-records-vs-walker), but
//! everything *after* the merger is identical.

use std::path::Path;

use crate::per_sample_pileup::psp::header::ParsedChromosome;
use crate::var_calling::dust_filter::{DustFilter, DustFilterConfig, DustFilterError};
use crate::var_calling::per_group_merger::{
    PerGroupMerger, PerGroupMergerConfig, PerGroupMergerError, SharedRefFetcher,
};
use crate::var_calling::per_position_merger::{PerPositionMergerError, PerPositionPileups};
use crate::var_calling::posterior_engine::{
    PosteriorEngine, PosteriorEngineConfig, PosteriorEngineError,
};
use crate::var_calling::variant_grouping::{GrouperConfig, GrouperError, VariantGrouper};
use crate::var_calling::vcf_writer::{
    CohortMetadata, CohortVcfWriter, VcfWriteError, WriterConfig, tmp_path_for,
};

/// Bundle of per-stage configs + shared resources consumed by
/// [`drive_cohort_pipeline`]. Constructed by the caller; consumed
/// (by-value) when the driver runs.
pub(crate) struct CohortPipelineParams {
    /// Whether to skip the DUST filter (set per-CLI by
    /// `--no-complexity-filter`).
    pub no_complexity_filter: bool,
    pub dust_cfg: DustFilterConfig,
    pub grouper_cfg: GrouperConfig,
    pub per_group_cfg: PerGroupMergerConfig,
    pub posterior_cfg: PosteriorEngineConfig,
    /// Reference fetcher shared between DUST and the per-group
    /// merger. The driver borrows from it for DUST and clones it for
    /// the per-group merger — both ends see the same underlying
    /// `Arc<dyn RefSeqFetcher>`.
    pub fetcher: SharedRefFetcher,
    /// Chromosomes table — consumed by DUST when the filter is on
    /// (no use otherwise).
    pub chromosomes: Vec<ParsedChromosome>,
}

/// Drive DUST → grouper → per-group merger → posterior engine →
/// VCF writer over `merger`, finalising the writer on success and
/// best-effort removing the `<output>.tmp` on driver-loop failure.
///
/// Returns the number of records written on success. Caller is
/// responsible for any post-driver bookkeeping (e.g. surfacing a
/// stashed walker error from the from-bam shim's
/// [`ErrorSheddingAdapter`](crate::pop_var_caller::cli::error_bridge::ErrorSheddingAdapter)).
///
/// Generic over the error type `E` so each subcommand can keep its
/// own typed error enum (`VarCallingCliError` /
/// `VarCallingFromBamCliError`); both enums already carry the five
/// required `From` impls so the bounds are satisfied without further
/// work.
pub(crate) fn drive_cohort_pipeline<M, E>(
    merger: M,
    params: CohortPipelineParams,
    output: &Path,
    metadata: CohortMetadata,
    writer_cfg: WriterConfig,
) -> Result<u64, E>
where
    M: Iterator<Item = Result<PerPositionPileups, PerPositionMergerError>>,
    E: From<DustFilterError>
        + From<GrouperError>
        + From<PerGroupMergerError>
        + From<PosteriorEngineError>
        + From<VcfWriteError>,
{
    let CohortPipelineParams {
        no_complexity_filter,
        dust_cfg,
        grouper_cfg,
        per_group_cfg,
        posterior_cfg,
        fetcher,
        chromosomes,
    } = params;

    let tmp_path = tmp_path_for(output);

    // DUST branch + grouper convergence — both arms hand the grouper
    // a single `Result<_, GrouperError>` iterator so downstream
    // construction is monomorphic. `DustFilter` borrows from
    // `&*fetcher` (the blanket `impl RefSeqFetcher for &T` makes
    // this work for `Arc<dyn RefSeqFetcher>` after one deref). The
    // anonymous `'_` lifetime ties the boxed iterator to the
    // local borrow of `fetcher`, which lives for the rest of this
    // function body.
    let upstream_for_grouper: Box<
        dyn Iterator<Item = Result<PerPositionPileups, GrouperError>> + '_,
    > = if no_complexity_filter {
        Box::new(merger.map(|r| r.map_err(GrouperError::from)))
    } else {
        let dust = DustFilter::new(merger, &*fetcher, chromosomes, dust_cfg);
        Box::new(dust.map(|r| r.map_err(GrouperError::from)))
    };

    let grouper = VariantGrouper::with_config(upstream_for_grouper, grouper_cfg);
    // `fetcher.clone()` is cheap (`Arc::clone`); the per-group
    // merger needs its own owned handle.
    let per_group = PerGroupMerger::with_config(grouper, fetcher.clone(), per_group_cfg);
    let posterior = PosteriorEngine::with_config(per_group, posterior_cfg);

    // Open writer, stream records, finalise. On any failure inside
    // the loop, best-effort remove `<output>.tmp` — `finish()`
    // consumes `self`, so a `?` short-circuit otherwise leaves the
    // tmp on disk. Matches the per-subcommand M6 cleanup
    // discipline.
    let mut writer = CohortVcfWriter::new(metadata, writer_cfg)?;
    let mut records_written: u64 = 0;
    let drive_result: Result<(), E> = (|| {
        for item in posterior {
            let record = item?;
            writer.write_record(&record)?;
            records_written += 1;
        }
        writer.finish()?;
        Ok(())
    })();
    if let Err(e) = drive_result {
        let _ = std::fs::remove_file(&tmp_path); // best-effort cleanup
        return Err(e);
    }
    Ok(records_written)
}
