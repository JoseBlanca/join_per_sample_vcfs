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

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use crate::per_sample_pileup::psp::header::ParsedChromosome;
use crate::per_sample_pileup::psp::{PspReadError, PspReader};
use crate::pop_var_caller::common::DEFAULT_BUFFERED_IO_CAPACITY;
use crate::var_calling::dust_filter::{DustFilter, DustFilterConfig, DustFilterError};
use crate::var_calling::per_group_merger::{
    PerGroupMerger, PerGroupMergerConfig, PerGroupMergerError, SharedRefFetcher,
};
use crate::var_calling::per_position_merger::{
    PerPositionMerger, PerPositionMergerError, PerPositionPileups,
};
use crate::var_calling::posterior_engine::{
    PosteriorEngine, PosteriorEngineConfig, PosteriorEngineError,
};
use crate::var_calling::variant_grouping::{GrouperConfig, GrouperError, VariantGrouper};
use crate::var_calling::vcf_writer::{
    CohortMetadata, CohortVcfWriter, VcfWriteError, WriterConfig, tmp_path_for,
};

/// Default for [`CohortPipelineParams::min_qual_phred`]. Matches
/// GATK HaplotypeCaller's `--standard-min-confidence-threshold-for-calling`
/// (30.0) as a hygiene gate.
///
/// **Caveat.** Our site-level QUAL is currently
/// `-10·log10(Π_s P(GT_s = hom-ref))` — a *product* over samples,
/// not a marginalised cohort-AF posterior the way GATK
/// (`P(AC = 0 | data)`, integrated over a shared AF) and FreeBayes
/// (`1 - P(K = 1 | reads)`) compute it. Under our formula QUAL
/// grows roughly linearly with cohort size even when nothing about
/// the per-sample evidence changes, so a fixed threshold has
/// different bite at different N. A fix to the QUAL semantics is
/// tracked separately; until then `30.0` is a coarse low-end
/// hygiene cutoff (drops the obvious tail) rather than a calibrated
/// quality bar.
pub const DEFAULT_MIN_QUAL_PHRED: f64 = 30.0;

/// Default for [`CohortPipelineParams::min_alt_obs_per_sample`].
/// Records where no ALT allele has `max(num_obs across samples) >= 2`
/// are dropped *before* reaching the EM. Empirically (multichrom
/// 10-duplicate cohort, 2026-05-20) `2` cuts FPs vs GATK by 70%
/// while losing < 1 pp recall; 89.9 % of pileup-level candidates
/// are single-alt-read positions that come from base / alignment
/// error rather than real variation.
///
/// `0` and `1` are both effective no-ops (every candidate has at
/// least one supporting alt read by construction); they're allowed
/// as escape hatches for tests / debugging.
///
/// Closest GATK analogue: `--min-pruning = 2` on the De Bruijn
/// assembly, which prunes assembly paths supported by fewer than
/// two reads in the cohort. Our filter is per-allele max across
/// samples rather than per-path, but the effect on candidate-set
/// size is comparable.
pub const DEFAULT_MIN_ALT_OBS_PER_SAMPLE: u32 = 2;

/// Counts the driver tracks as it streams records from the posterior
/// engine to the VCF writer. Surfaced by [`drive_cohort_pipeline`]
/// and [`process_one_chromosome`] so the orchestrator can build a
/// run summary without reading the emitted VCF back.
#[doc(hidden)]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct CohortDriveStats {
    /// Total records written to the writer (or the per-chrom fragment).
    pub records_written: u64,
    /// Subset of `records_written` whose posterior EM hit the
    /// iteration cap without satisfying the convergence threshold —
    /// emitted with `EmDiagnostics::converged = false` and routed by
    /// the writer to `FILTER=EMNoConv`.
    pub records_unconverged: u64,
    /// Posterior records dropped before reaching the writer because
    /// every sample's argmax genotype was hom-ref (so the cohort-level
    /// `AC` would have been 0 across all ALTs). The EM allele
    /// frequency `p̂` can still be non-zero at such sites — it's a
    /// likelihood-weighted estimate from read counts, not from argmax
    /// genotypes — but emitting them produces non-variant records,
    /// which the VCF spec doesn't intend for sites VCFs.
    pub records_dropped_hom_ref: u64,
    /// Posterior records dropped before reaching the writer because
    /// `qual_phred` was strictly below
    /// [`CohortPipelineParams::min_qual_phred`]. See the constant's
    /// doc for the QUAL-semantics caveat.
    pub records_dropped_low_qual: u64,
    /// Per-group records dropped *before* the EM because no ALT
    /// allele had `max(num_obs across samples) >=
    /// min_alt_obs_per_sample`. These records would otherwise have
    /// reached the EM and emitted a 0/1 call from a single supporting
    /// read; pruning them at the per-group → EM boundary saves the EM
    /// CPU on top of the precision gain. See
    /// [`DEFAULT_MIN_ALT_OBS_PER_SAMPLE`].
    pub records_dropped_low_alt_obs: u64,
}

/// Bundle of per-stage configs + shared resources consumed by
/// [`drive_cohort_pipeline`]. Constructed by the caller; consumed
/// (by-value) when the driver runs.
///
/// `Clone` is derived because the per-chromosome parallel path in
/// [`run_var_calling`](super::var_calling::run_var_calling) builds
/// one parameter set per chromosome from a shared template. Every
/// field is cheap-Clone — the configs are `Copy` or wrap a small
/// `Vec`; `fetcher` is an `Arc`; `chromosomes` is a `Vec<ParsedChromosome>`.
#[doc(hidden)]
#[derive(Clone)]
pub struct CohortPipelineParams {
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
    /// Minimum site-level `QUAL` (phred) required to emit a record.
    /// Records with `qual_phred < min_qual_phred` are dropped before
    /// reaching the writer. `0.0` disables the filter. Default:
    /// [`DEFAULT_MIN_QUAL_PHRED`] — see the constant's doc for the
    /// caveat about how our cohort-product QUAL scales with sample
    /// count.
    pub min_qual_phred: f64,
    /// Minimum per-allele alt-read support — at the per-group merger
    /// → EM boundary, a record is dropped if no ALT allele has
    /// `max(num_obs across samples) >= min_alt_obs_per_sample`. `0`
    /// or `1` disables the filter (every candidate has ≥ 1 supporting
    /// read by construction). Default: [`DEFAULT_MIN_ALT_OBS_PER_SAMPLE`].
    pub min_alt_obs_per_sample: u32,
}

/// Drive DUST → grouper → per-group merger → posterior engine →
/// VCF writer over `merger`, finalising the writer on success and
/// best-effort removing the `<output>.tmp` on driver-loop failure.
///
/// Returns [`CohortDriveStats`] on success — total emitted records
/// + the subset whose EM didn't converge. The caller is responsible
/// for post-driver bookkeeping (e.g. surfacing a stashed walker
/// error from the from-bam shim's
/// [`ErrorSheddingAdapter`](crate::pop_var_caller::cli::error_bridge::ErrorSheddingAdapter)).
///
/// Generic over the error type `E` so each subcommand can keep its
/// own typed error enum (`VarCallingCliError` /
/// `VarCallingFromBamCliError`); both enums already carry the five
/// required `From` impls so the bounds are satisfied without further
/// work.
#[doc(hidden)]
pub fn drive_cohort_pipeline<M, E>(
    merger: M,
    params: CohortPipelineParams,
    output: &Path,
    metadata: CohortMetadata,
    writer_cfg: WriterConfig,
) -> Result<CohortDriveStats, E>
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
        min_qual_phred,
        min_alt_obs_per_sample,
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

    // Min-alt-AD filter at the per-group merger → EM boundary. The
    // filter closure cannot mutably borrow `stats` (that borrow is
    // taken by the drive loop below), so we accumulate into a local
    // `AtomicU64` and fold the result into `stats` once the iterator
    // is fully drained. Single-threaded in practice, but `AtomicU64`
    // sidesteps the lifetime gymnastics of `Cell<u64>` in a `move`
    // closure here.
    let alt_obs_dropped = std::sync::atomic::AtomicU64::new(0);
    let alt_obs_dropped_ref = &alt_obs_dropped;
    let per_group_filtered = per_group.filter_map(move |item| match item {
        Err(e) => Some(Err(e)),
        Ok(record) => {
            if min_alt_obs_per_sample <= 1 || record.alleles.len() <= 1 {
                return Some(Ok(record));
            }
            let n_alleles = record.alleles.len();
            let any_alt_passes = (1..n_alleles).any(|allele_idx| {
                let mut max_obs = 0u32;
                for sample_idx in 0..record.n_samples {
                    let obs = record.scalars[sample_idx * n_alleles + allele_idx].num_obs;
                    if obs > max_obs {
                        max_obs = obs;
                    }
                }
                max_obs >= min_alt_obs_per_sample
            });
            if any_alt_passes {
                Some(Ok(record))
            } else {
                alt_obs_dropped_ref.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                None
            }
        }
    });
    let posterior = PosteriorEngine::with_config(per_group_filtered, posterior_cfg);

    // Open writer, stream records, finalise. On any failure inside
    // the loop, best-effort remove `<output>.tmp` — `finish()`
    // consumes `self`, so a `?` short-circuit otherwise leaves the
    // tmp on disk. Matches the per-subcommand M6 cleanup
    // discipline.
    let mut writer = CohortVcfWriter::new(metadata, writer_cfg)?;
    let mut stats = CohortDriveStats::default();
    let drive_result: Result<(), E> = (|| {
        for item in posterior {
            let record = item?;
            if !record.is_variant_call() {
                stats.records_dropped_hom_ref += 1;
                continue;
            }
            if record.qual_phred < min_qual_phred {
                stats.records_dropped_low_qual += 1;
                continue;
            }
            if !record.diagnostics.converged {
                stats.records_unconverged += 1;
            }
            writer.write_record(&record)?;
            stats.records_written += 1;
        }
        writer.finish()?;
        Ok(())
    })();
    if let Err(e) = drive_result {
        let _ = std::fs::remove_file(&tmp_path); // best-effort cleanup
        return Err(e);
    }
    stats.records_dropped_low_alt_obs +=
        alt_obs_dropped.load(std::sync::atomic::Ordering::Relaxed);
    Ok(stats)
}

/// Run the cohort var-calling chain (DUST → … → VCF writer) over a
/// single chromosome's records and write a complete self-contained
/// VCF fragment to `fragment_path`. Designed as the per-worker body
/// of the per-chromosome parallel partition driven from
/// [`run_var_calling`](super::var_calling::run_var_calling).
///
/// Each worker re-opens its own `PspReader` set against `psp_paths`
/// — opens are cheap (one seek to the trailer + a small index read
/// per file) and per-worker readers avoid cross-thread coordination
/// on the read side. The whole-chromosome window is materialised via
/// [`PspReader::region_records`](crate::per_sample_pileup::psp::PspReader::region_records)
/// with `start = 1`, `end = chrom.length`; chromosomes with zero
/// records on disk yield an empty iterator + a header-only fragment,
/// which [`crate::var_calling::vcf_writer::concat::concat_fragments`]
/// handles by stripping the header on non-first fragments.
///
/// Returns `(chrom_id, CohortDriveStats)` on success. Errors are
/// surfaced through the same typed enum the sequential
/// [`drive_cohort_pipeline`] uses; the per-chrom orchestrator
/// collects worker results and returns the first error after every
/// worker has joined.
#[doc(hidden)]
#[allow(clippy::too_many_arguments)] // shared template + per-worker inputs; bundling adds indirection without clarity gain
pub fn process_one_chromosome<E>(
    chrom_id: u32,
    psp_paths: &[PathBuf],
    sample_names: Vec<String>,
    all_chromosomes: Vec<ParsedChromosome>,
    fragment_path: PathBuf,
    metadata: CohortMetadata,
    writer_cfg_template: WriterConfig,
    pipeline_params: CohortPipelineParams,
) -> Result<(u32, CohortDriveStats), E>
where
    E: From<std::io::Error>
        + From<PspReadError>
        + From<PerPositionMergerError>
        + From<DustFilterError>
        + From<GrouperError>
        + From<PerGroupMergerError>
        + From<PosteriorEngineError>
        + From<VcfWriteError>,
{
    let chrom_length = all_chromosomes
        .get(chrom_id as usize)
        .map(|c| c.length)
        .expect("chrom_id past contig table — caller (run_var_calling) must validate the bound");

    // Open one PspReader per sample. Per-worker ownership: no
    // cross-thread coordination on the read side.
    let mut readers: Vec<PspReader<BufReader<File>>> = Vec::with_capacity(psp_paths.len());
    for path in psp_paths {
        let file = File::open(path)?;
        let buf = BufReader::with_capacity(DEFAULT_BUFFERED_IO_CAPACITY, file);
        readers.push(PspReader::new(buf)?);
    }

    // Per-chrom record iterators. `region_records` seeks to the
    // first overlapping block + clamps the iterator to records on
    // this chrom inside [1, chrom_length].
    let record_iters: Vec<_> = readers
        .iter_mut()
        .map(|r| r.region_records(chrom_id, 1, chrom_length))
        .collect();

    let merger = PerPositionMerger::new(record_iters, sample_names, all_chromosomes)?;

    // Inherit `emit_gp` (and any future flags) from the template;
    // override only the output path. Struct-update syntax keeps this
    // robust against future `WriterConfig` field additions.
    let writer_cfg = WriterConfig {
        output: fragment_path.clone(),
        ..writer_cfg_template
    };

    let stats = drive_cohort_pipeline::<_, E>(
        merger,
        pipeline_params,
        &fragment_path,
        metadata,
        writer_cfg,
    )?;
    Ok((chrom_id, stats))
}
