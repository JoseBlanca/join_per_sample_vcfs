//! Driver / wiring — section orchestration (appendix §G).
//!
//! *(today: `var_calling::driver` `drive_blocks_parallel` + the
//! `pop_var_caller::var_calling::run_var_calling` CLI entry)*
//!
//! Wires the three sections producer → caller → writer and exposes the new
//! cohort `.psp` → VCF entry point ([`run_var_calling`]), the one the
//! byte-identity oracle drives against the old `var_calling::` path.
//!
//! **Phase 4: single-threaded.** The pipeline runs the producer
//! ([`CohortChunkIntegrator`]) inline, calling the caller ([`VariantCaller`])
//! and writer ([`VcfWriter`]) for each chunk in genomic order. This is enough
//! to turn the byte-identity oracle green — the producer already emits in
//! `chunk_order` order, so the writer's reorder is a no-op here. The bounded
//! `crossbeam-channel` parallel topology (producer thread / W caller threads /
//! writer thread; the writer's `BTreeMap` reorder then does real work) is a
//! follow-up perf step that does not change the bytes (the writer reorders by
//! `chunk_order` regardless).

use std::fs::File;
use std::io::BufReader;
use std::ops::Range;

use thiserror::Error;

use crate::fasta::{
    ChromRefFetchError, ChromRefFetcher, ManualEvictChromRefFetcher, StreamingChromRefFetcher,
};
use crate::pop_var_caller::var_calling::VarCallingArgs;
use crate::psp::{PspReadError, PspReader};
use crate::var_calling_new::cohort_integration::{CohortChunkIntegrator, ProducerError};
use crate::var_calling_new::dust_filter::{MIN_DUST_HALO, sdust_mask_for_span};
use crate::var_calling_new::em_posterior_calc::{CallerError, VariantCaller};
use crate::var_calling_new::per_group_merger::{DEFAULT_BATCH_SIZE, PerGroupMergerConfig};
use crate::var_calling_new::per_position_merger::{
    PerPositionMergerError, check_chromosome_agreement,
};
use crate::var_calling_new::posterior_engine::PosteriorEngineConfig;
use crate::var_calling_new::sample_reader::SamplePspReader;
use crate::var_calling_new::variant_grouping::GrouperConfig;
use crate::var_calling_new::vcf_writer::{DownstreamFilters, VcfWriter, WriterError};
use crate::vcf::{CohortMetadata, WriterConfig};

/// Open-file buffer size (matches `pop_var_caller::common::DEFAULT_BUFFERED_IO_CAPACITY`).
const BUFFERED_IO_CAPACITY: usize = 64 * 1024;
/// Default per-chunk variable-position target (matches `DEFAULT_DESIRED_VARIANTS_PER_BLOCK`).
const DEFAULT_TARGET_VARIANTS: u32 = 1024;
/// DUST sub-span for the resident-buffer bound. The mask is byte-identical for
/// any sub-span (runs are coalesced across boundaries); this just caps RAM.
const DUST_SUBSPAN: u32 = 1_000_000;

/// Errors surfaced by the re-architected pipeline entry point.
#[derive(Debug, Error)]
pub enum PipelineError {
    #[error("opening / reading a .psp file: {0}")]
    Io(#[from] std::io::Error),
    #[error("decoding a .psp header: {0}")]
    Psp(#[from] PspReadError),
    #[error("cohort chromosome agreement: {0}")]
    ChromAgreement(#[from] PerPositionMergerError),
    #[error("invalid configuration: {0}")]
    Config(String),
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

/// Re-architected cohort `.psp` → VCF entry point.
///
/// Mirrors the old CLI-level
/// [`run_var_calling`](crate::pop_var_caller::var_calling::run_var_calling) so
/// the byte-identity oracle can drive both pipelines from the same
/// [`VarCallingArgs`] and diff the resulting VCFs. The argument-type coupling
/// to the CLI layer is temporary: at the P7 swap this becomes the production
/// entry the CLI invokes directly.
///
/// (Input *validation* the old entry performs — reference-basename check,
/// FASTA-vs-`.psp` per-contig MD5 cross-check — is intentionally omitted here:
/// it guards against bad input but does not affect the emitted VCF bytes, so
/// it is out of scope for the byte-identity gate. Contamination estimates are
/// likewise not yet wired; the oracle path passes none.)
pub fn run_var_calling(args: &VarCallingArgs) -> Result<(), PipelineError> {
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

    // --- Build every per-stage config from the args. ---
    let cfg_err = |e: &dyn std::fmt::Display| PipelineError::Config(e.to_string());
    let grouper_cfg = GrouperConfig::new(cohort.var_group_max_span).map_err(|e| cfg_err(&e))?;
    let merger_cfg = PerGroupMergerConfig::new(
        cohort.ploidy,
        cohort.max_alleles_per_var,
        cohort.max_alleles_lh_calc,
        DEFAULT_BATCH_SIZE,
    )
    .map_err(|e| cfg_err(&e))?;
    let posterior_cfg = PosteriorEngineConfig::new()
        .with_convergence_threshold(cohort.em_convergence_threshold)
        .map_err(|e| cfg_err(&e))?
        .with_max_iterations(cohort.em_max_iterations)
        .map_err(|e| cfg_err(&e))?
        .with_ref_pseudocount(cohort.ref_pseudocount)
        .map_err(|e| cfg_err(&e))?
        .with_snp_alt_pseudocount(cohort.snp_alt_pseudocount)
        .map_err(|e| cfg_err(&e))?
        .with_indel_alt_pseudocount(cohort.indel_alt_pseudocount)
        .map_err(|e| cfg_err(&e))?
        .with_compound_alt_pseudocount(cohort.compound_alt_pseudocount)
        .map_err(|e| cfg_err(&e))?
        .with_fixation_index_default(cohort.inbreeding_coefficient)
        .map_err(|e| cfg_err(&e))?
        .with_max_gq_phred(cohort.max_gq_phred)
        .map_err(|e| cfg_err(&e))?
        .with_contamination(None)
        .map_err(|e| cfg_err(&e))?;

    let min_alt_obs = cohort.min_alt_obs_per_sample;
    let max_group_span = cohort.var_group_max_span;
    let target_variants = if args.target_variants_per_chunk == 0 {
        DEFAULT_TARGET_VARIANTS
    } else {
        args.target_variants_per_chunk
    };

    // --- Per-sample segment readers (region overridden per interval). ---
    let readers: Vec<SamplePspReader<BufReader<File>>> = psp_readers
        .into_iter()
        .map(|r| SamplePspReader::new(r, 0, 1, 1))
        .collect();

    // --- The three sections. ---
    let metadata = CohortMetadata {
        sample_names: sample_names.clone(),
        contigs: chromosomes.clone(),
        tool_string: format!("pop_var_caller {}", env!("CARGO_PKG_VERSION")),
        command_line: "var-calling".to_string(),
    };
    let writer_config = WriterConfig::new(args.output.clone()).with_emit_gp(cohort.emit_gp);
    let filters = DownstreamFilters {
        min_qual_phred: cohort.min_qual_phred,
        no_mapq_diff_filter: cohort.no_mapq_diff_filter,
        min_mapq_diff_t: cohort.min_mapq_diff_t,
    };
    let mut writer = VcfWriter::new(metadata, writer_config, filters)?;
    let caller = VariantCaller::new(grouper_cfg, merger_cfg, posterior_cfg, min_alt_obs);
    let mut producer =
        CohortChunkIntegrator::new(readers, sample_names, min_alt_obs, target_variants);

    // --- Per-chromosome REF + dust fetchers (built on demand, monotonic). ---
    let mut ref_fetcher: Option<(u32, StreamingChromRefFetcher)> = None;
    let mut dust_fetcher: Option<(u32, ManualEvictChromRefFetcher)> = None;

    // --- Drive: chromosome → covered interval → chunk → call → write. ---
    for chrom_id in 0..n_chromosomes {
        let intervals = producer.covered_intervals(chrom_id, max_group_span);
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
                    let mut fetch = |start: u32, len: u32| -> Result<Vec<u8>, String> {
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
                let called = caller.call_chunk(chunk)?;
                writer.handle(called)?;
            }
        }
    }

    writer.finish()?;
    Ok(())
}

/// Fetch REF bytes for `chrom_id`, building a fresh per-contig
/// [`StreamingChromRefFetcher`] when the chromosome changes (fetches within a
/// chromosome are monotonic-forward across chunks, satisfying the fetcher's
/// sliding-window contract). Errors are stringified for the producer's closure.
fn ref_fetch(
    cache: &mut Option<(u32, StreamingChromRefFetcher)>,
    chrom_id: u32,
    fasta: &std::path::Path,
    chrom_names: &[String],
    start: u32,
    len: u32,
) -> Result<Vec<u8>, String> {
    if cache.as_ref().map(|(c, _)| *c) != Some(chrom_id) {
        let f = StreamingChromRefFetcher::for_contig(fasta, &chrom_names[chrom_id as usize])
            .map_err(|e| e.to_string())?;
        *cache = Some((chrom_id, f));
    }
    let (_, fetcher) = cache.as_mut().expect("just set");
    // `ChromRefFetcher::fetch` (trait) returns an owned `Vec<u8>`.
    fetcher.fetch(start, len).map_err(|e| e.to_string())
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

/// The interval's low-complexity (sdust) mask — copied verbatim from
/// `driver::dust_mask_for_interval`. Scans the interval in `subspan` windows
/// (evicting the buffer between them to bound resident RAM) and coalesces runs
/// across the window boundaries, so the result equals a single whole-interval
/// scan (byte-identical regardless of `subspan`).
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
