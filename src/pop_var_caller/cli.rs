//! `pop_var_caller pileup` — the Stage 1 CLI orchestrator.
//!
//! Glues the per-region read source (`SegmentMergedReads` over pooled,
//! re-seekable `AlignmentFile`s) → `BaqStream` (or BAQ-bypass passthrough)
//! → pileup walker → `.psp` writer. Argument parsing is delegated to
//! clap; everything else is in [`run_pileup`].
//!
//! Plan: `ia/feature_implementation_plans/pop_var_caller_pileup_cli.md`.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};

use clap::{Args, Parser, Subcommand};
use thiserror::Error;

use crate::bam::alignment_input::{
    AlignmentMergedReaderConfig, FilterCounts, build_fasta_repository, load_pileup_inputs,
};
use crate::bam::errors::AlignmentInputError;
use crate::bam::segment_merge::SegmentMergedReads;
use crate::bam::segment_reader::{AlignmentFile, SegmentReadFilter};
use crate::baq::BaqConfig;
use crate::fasta::{ContigList, RepositoryRefFetcher};
use crate::pileup::per_sample::baq_stream::BaqSkipCounts;
use crate::pileup::per_sample::pileup_to_psp::{PileupToPspError, drive_region_into_writer};
use crate::pileup::per_sample::read_processor::ReadProcessingConfig;
use crate::pileup::walker::{RunSummary, WalkerConfig};
use crate::pop_var_caller::common::{
    DEFAULT_BUFFERED_IO_CAPACITY, basename, configure_rayon_pool, current_command_line,
    format_md5_hex, rfc3339_now,
};
use crate::psp::header::{ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance};
use crate::psp::writer::{DEFAULT_BLOCK_WINDOW_BP, MAX_BLOCK_TARGET_BYTES, PspWriter};
use crate::regions::{ContigBounds, RegionSet};

pub mod error_bridge;
pub mod parsers;
pub mod shared_args;

// ---------------------------------------------------------------------
// Clap surface
// ---------------------------------------------------------------------

/// Top-level CLI for the `pop_var_caller` binary.
#[derive(Debug, Parser)]
#[command(name = "pop_var_caller", version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub cmd: PopVarCallerCommand,
}

#[derive(Debug, Subcommand)]
pub enum PopVarCallerCommand {
    /// Stage 1: run BAQ + pileup over one sample's CRAMs, emit a .psp.
    Pileup(PileupArgs),

    /// Stream a .psp as samtools-mpileup-style text plus a trailing
    /// column with per-allele aggregates PSP carries.
    PspToPileup(super::psp_to_pileup::PspToPileupArgs),

    /// Estimate per-sample contamination from cohort `.psp` files and
    /// emit a TOML artefact consumed by `var-calling`.
    EstimateContamination(super::estimate_contamination::EstimateContaminationArgs),

    /// Call SNPs across a cohort: `.psp` files → multi-sample VCF.
    VarCalling(super::var_calling::VarCallingArgs),

    /// SSR Stage 0: detect tandem repeats in a reference and write the
    /// per-genome SSR locus catalog (a bgzip TSV).
    SsrCatalog(super::ssr_catalog::SsrCatalogArgs),

    /// SSR Stage 1: genotype one sample's BAM/CRAM against an SSR catalog,
    /// emitting a per-locus `.ssr.psp` evidence file.
    SsrPileup(super::ssr_pileup::SsrPileupArgs),

    /// SSR Stage 2: merge per-sample `.ssr.psp` evidence across a cohort and
    /// genotype each catalog locus into a multi-sample VCF.
    SsrCall(super::ssr_call::SsrCallArgs),
}

/// Arguments accepted by the `pileup` subcommand. The struct is the
/// authoritative knob list; the orchestrator translates it into module
/// configs. Defaults are pulled from each module's `DEFAULT_*` consts
/// so they stay in sync.
///
/// Stage 1 knobs (CRAM-input filters, BAQ HMM, pileup walker) live in
/// the flattened [`Stage1Args`](shared_args::Stage1Args) sub-struct so
/// the `var-calling-from-bam` subcommand reuses the same surface
/// — M10 from the 2026-05-19 cohort CLI review.
#[derive(Debug, Args, Clone)]
pub struct PileupArgs {
    // ===== Common flags (visible in `-h`) =====================
    /// Reference FASTA. A sibling `.fai` is required.
    #[arg(long)]
    pub reference: PathBuf,

    /// Output .psp path. Must be a regular file; stdout is not
    /// supported (the .psp format is random-access on read).
    #[arg(long)]
    pub output: PathBuf,

    /// BED file restricting analysis to the listed regions; the rest of
    /// the genome is ignored. Coordinates are 0-based half-open
    /// (standard BED). Without this flag the whole genome — every
    /// contig in the reference — is analysed. The reads are seeked to
    /// via the alignment index rather than streamed, so a sparse BED is
    /// cheap.
    #[arg(long)]
    pub regions: Option<PathBuf>,

    /// Build a missing alignment index (`.crai` for CRAM, `.csi` for
    /// BAM) in place instead of erroring. Off by default: a missing
    /// index is a hard error naming the path it looked for. The
    /// reference `.fai` is always required and is never built.
    #[arg(long)]
    pub build_map_file_index: bool,

    /// Worker thread count for the read-processing pipeline. If omitted,
    /// defaults to 4 — a good balance for typical hardware. The
    /// pipeline's scaling flattens by ~4 workers, so raising this past
    /// your machine's performance-core count tends not to help (and can
    /// hurt on hybrid performance/efficiency CPUs).
    #[arg(long)]
    pub threads: Option<usize>,

    /// One or more coordinate-sorted CRAM or BAM file(s) for one
    /// sample. All inputs in one invocation must share the same
    /// format (CRAM-only or BAM-only); mixed CRAM + BAM is rejected
    /// with a typed error.
    #[arg(required = true)]
    pub alignment_files: Vec<PathBuf>,

    /// PSP writer **safety cap** (uncompressed bytes per block): force a
    /// block flush mid-window if projected bytes reach this, guarding
    /// against a pathological window blowing up resident memory. The
    /// primary cut is now `--block-window-bp` (a fixed genomic grid);
    /// this cap should sit well above a typical window's size so it only
    /// fires on extreme density. Range 16 KiB..=16 MiB, default 16 MiB.
    ///
    /// Trade-off — peak memory vs on-disk size. The cohort
    /// var-calling driver opens one `PspReader` per sample per
    /// chromosome worker, each holding one decoded block live, so
    /// peak heap scales as `n_threads × N × per_block`. Smaller
    /// blocks → lower cohort-step memory, slightly worse zstd
    /// compression context (so larger .psp on disk).
    ///
    /// Sweep on tomato1 (N=18, T=4): 16 MiB → 2501 MB peak / 183 MB
    /// on disk; 1 MiB → 261 MB / 206 MB (+13% disk, default);
    /// 256 KiB → 108 MB / 246 MB (+34% disk); 64 KiB → 59 MB /
    /// 321 MB (+75% disk, near the knee). Wall time is flat across
    /// the range — no CPU penalty.
    ///
    /// Dial down (256 KiB, 64 KiB) for large-cohort joint
    /// genotyping that hits a memory ceiling. Dial up (4 MiB,
    /// 16 MiB) for single-sample archival where on-disk size
    /// matters more than cohort-step heap.
    #[arg(
        long,
        default_value_t = MAX_BLOCK_TARGET_BYTES,
        value_parser = parsers::parse_block_target_bytes,
        help_heading = "Advanced — PSP writer",
    )]
    pub block_target_bytes: usize,

    /// PSP block-cut genomic window, in bases (the **primary** cut).
    /// Blocks are cut on a fixed reference grid (`pos / window`), so
    /// every independently-written sample cuts at the same positions.
    /// This alignment is what lets the cohort `var-calling` reader
    /// advance all samples in lockstep — without it, misaligned blocks
    /// inflate the per-block synchronisation count ~N-fold and the read
    /// stops scaling with threads. PSP stores a fixed per-position
    /// summary, so a window ≈ that many records ≈ roughly-constant
    /// memory per block regardless of read depth. Default 20 kb
    /// (~2.5 MB blocks). `--block-target-bytes` is the safety cap for
    /// pathological windows.
    #[arg(
        long,
        default_value_t = DEFAULT_BLOCK_WINDOW_BP,
        help_heading = "Advanced — PSP writer",
    )]
    pub block_window_bp: u32,

    // ===== Stage 1 (shared with var-calling-from-bam) =========
    #[command(flatten)]
    pub stage1: shared_args::Stage1Args,
}

pub(crate) fn parse_mismatch_fraction(s: &str) -> Result<f32, String> {
    let v: f32 = s.parse().map_err(|e| format!("not a number: {e}"))?;
    if !v.is_finite() {
        return Err(format!("must be finite, got `{s}`"));
    }
    if !(0.0..=1.0).contains(&v) {
        return Err(format!("must be in [0.0, 1.0], got `{s}`"));
    }
    Ok(v)
}

// ---------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum PileupCliError {
    #[error("alignment input: {0}")]
    AlignmentInput(#[from] AlignmentInputError),
    #[error("pipeline: {0}")]
    Pipeline(#[from] PileupToPspError),
    #[error("regions BED: {0}")]
    Bed(#[from] crate::regions::BedError),
    #[error(
        "contig '{contig}' has no @SQ M5 checksum in the input alignment file(s); \
         re-create the file with a tool that emits @SQ M5 (e.g. samtools view -t)"
    )]
    MissingMd5 { contig: String },
    #[error("rayon thread pool already initialised — refusing to override")]
    RayonAlreadyConfigured,
    #[error("io: {0}")]
    Io(#[from] io::Error),
    #[error("contig '{name}' length {length} exceeds u32::MAX")]
    ContigLengthOverflow { name: String, length: u64 },
    #[error("internal: failed to format current timestamp as RFC3339")]
    TimestampFormat,
}

// ---------------------------------------------------------------------
// Top-level driver
// ---------------------------------------------------------------------

/// Drive the region-by-region pileup and emit one `.psp`.
///
/// There is a single code path: the pileup always operates on a
/// [`RegionSet`] — the `--regions` BED, or (without it) one full-length
/// span per contig. One pooled, re-seekable reader per input file is
/// opened once for the whole run; for each region the per-file index is
/// used to *seek* to the reads overlapping it, and those per-file streams
/// are k-way-merged ([`SegmentMergedReads`](crate::bam::segment_merge::SegmentMergedReads)).
/// The BAQ → walker chain
/// ([`with_stage1_chain`](super::stage1_pipeline::with_stage1_chain)) runs
/// over the merged reads, and only the columns inside the region are
/// written to one shared [`PspWriter`]. "Whole genome" is not a special
/// case — it is the region set whose every span covers an entire
/// contig.
///
/// Records are written to `<output>.tmp` and atomically renamed to
/// `<output>` on success. Stderr receives a one-shot run-summary block
/// totalling every region's counters.
/// Default Stage 1 worker-thread count when `--threads` is omitted. A
/// deliberately modest, hardware-agnostic value: the pipeline's wall
/// time flattens by ~4 workers, and defaulting to all logical cores
/// oversubscribes (badly on hybrid performance/efficiency CPUs).
const DEFAULT_PILEUP_THREADS: usize = 4;

pub fn run_pileup(args: &PileupArgs) -> Result<(), PileupCliError> {
    let stage1 = &args.stage1;

    // 1. Translate CLI → module configs.
    let alignment_cfg = alignment_config_from_args(stage1);
    let baq_cfg = baq_config_from_args(stage1);
    let proc_cfg = read_processing_config_from_args(stage1);
    let walker_cfg = walker_config_from_args(stage1);

    // 2. Resolve the thread budget. Taken from the flag directly (not
    //    `rayon::current_num_threads()`, which would lazily spin up a
    //    rayon pool we may not use). When omitted we default to
    //    `DEFAULT_PILEUP_THREADS` rather than all logical cores: the
    //    pipeline's scaling flattens by ~4 workers, and grabbing every
    //    core oversubscribes — worse on hybrid (performance + efficiency)
    //    CPUs, where the extra workers land on the slow cores.
    let n_threads = args
        .threads
        .filter(|&t| t > 0)
        .unwrap_or(DEFAULT_PILEUP_THREADS);
    // Only the inline read-processing path uses rayon (`BaqStream`'s
    // `par_drain`); the staged pipeline uses dedicated threads. Sizing
    // the global rayon pool for the pipeline path would leave `n` idle
    // rayon threads alongside the pipeline's own `n` workers (~2n threads
    // for an n-thread request), so size it only when the inline path runs.
    if n_threads < super::stage1_pipeline::STAGED_MIN_THREADS {
        configure_rayon_pool(Some(n_threads))
            .map_err(|_| PileupCliError::RayonAlreadyConfigured)?;
    }

    // 3. Load metadata + per-input handles. Index pre-flight runs here:
    //    a missing alignment index is built (with --build-map-file-index)
    //    or hard-errors. The reference `.fai` is required.
    let inputs = load_pileup_inputs(
        &args.alignment_files,
        &args.reference,
        args.build_map_file_index,
    )?;

    // 4. Writer header (validates per-contig u32 lengths + @SQ M5). Built
    //    before the region set so the latter can borrow the header's
    //    already-u32-validated contig bounds.
    let mut header = build_writer_header(
        &inputs.sample_name,
        &args.reference,
        &inputs.contigs,
        &args.alignment_files,
        stage1,
        args.block_target_bytes,
        args.block_window_bp,
        args.regions.as_deref(),
        n_threads,
    )?;

    // 5. Region set: the BED, or one full-length span per contig. The
    //    contig slice index is the `chrom_id` used by the writer and the
    //    `SegmentMergedReads` reader, so this preserves the reference's
    //    contig order.
    let region_set = {
        let contig_bounds: Vec<ContigBounds> = header
            .chromosomes
            .iter()
            .map(|c| ContigBounds {
                name: &c.name,
                length: c.length,
            })
            .collect();
        match &args.regions {
            Some(bed_path) => RegionSet::from_bed_path(bed_path, &contig_bounds)?,
            None => RegionSet::whole_contigs(&contig_bounds),
        }
        // `contig_bounds` borrows `header`; drop it here so `header` can
        // move into the writer below.
    };
    // Record how many analysis spans the BED resolved to (provenance);
    // omitted for the whole-genome default.
    if args.regions.is_some() {
        header.writer.parameters.insert(
            "regions_count".to_string(),
            ParameterValue::Integer(region_set.len() as i64),
        );
    }
    // Announce the resolved analysis mode on stderr so a run is
    // self-describing without inspecting the .psp header (review Mi4): a
    // BED restricts to N spans; its absence means the whole genome.
    match &args.regions {
        Some(bed_path) => eprintln!(
            "regions   : {} spans from {}",
            region_set.len(),
            bed_path.display()
        ),
        None => eprintln!("regions   : whole genome (no --regions BED)"),
    }

    // 6. Open the shared writer on the .tmp path.
    let tmp_path = tmp_path_for(&args.output);
    let file = File::create(&tmp_path).map_err(PileupCliError::Io)?;
    let buf = BufWriter::with_capacity(DEFAULT_BUFFERED_IO_CAPACITY, file);
    let mut writer = PspWriter::new_with_block_layout(
        buf,
        header,
        args.block_target_bytes,
        args.block_window_bp,
    )
    .map_err(PileupToPspError::from)?;

    // 7. Drive each region into the shared writer, totalling counters.
    //
    //    The FASTA repository (CRAM reference resolver + F1/F3
    //    read-reference source) is built **once** here and shared across
    //    every region, rather than rebuilt per region.
    //    The noodles repository is a whole-contig cache, so a per-region
    //    rebuild reloaded the entire contig for every region on it — the
    //    `--regions` perf/memory regression. `RegionSet` is sorted by
    //    `(chrom_id, start)`, so regions arrive grouped by contig; we
    //    `clear()` the cache on each contig transition to keep just one
    //    contig's sequence resident (matching the whole-genome path's
    //    memory profile). See
    //    `doc/devel/implementation_plans/fasta_reference_reading_unification.md`.
    let repository = build_fasta_repository(&args.reference)?;
    let mut current_chrom_id: Option<u32> = None;

    // The walker reads its reference windows from the same shared
    // `repository` the reader keeps resident (CRAM decode + the per-read
    // F1/F3 fetch already load each contig there, for both CRAM and
    // BAM). So the walker no longer opens a second, independent
    // streaming reader over the same FASTA — it shares the resident
    // contig at zero extra memory, observing the per-contig `clear()`
    // through the shared `Arc`.
    let walker_fetcher = RepositoryRefFetcher::new(repository.clone(), inputs.contigs.clone());

    // Open one pooled, re-seekable reader per input file, once for the
    // whole run. Per region we borrow a segment iterator from each and
    // k-way-merge them (`SegmentMergedReads`) — replacing the old
    // per-region `AlignmentMergedReader::query`, which re-opened every
    // input and re-walked its index on each region (the `--regions` tax).
    // The shared `repository` (cleared per contig transition below) is
    // handed to every file; BAM ignores it, CRAM decodes against it.
    let segment_filter = SegmentReadFilter::from(&alignment_cfg);
    let segment_files: Vec<AlignmentFile> = args
        .alignment_files
        .iter()
        .enumerate()
        .map(|(input_index, path)| {
            AlignmentFile::from_input(
                path.clone(),
                inputs.headers[input_index].clone(),
                inputs.indexes[input_index].clone(),
                Some(repository.clone()),
                segment_filter,
                input_index,
            )
        })
        .collect::<Result<_, _>>()?;

    let mut total_filter = FilterCounts::default();
    let mut total_walker = RunSummary::default();
    let mut total_baq_skip = (!stage1.no_baq).then(BaqSkipCounts::default);
    let mut stashed_upstream: Option<AlignmentInputError> = None;

    for region in region_set.iter() {
        // Contig transition: drop the previous contig's cached sequence.
        // Safe because the sorted RegionSet never revisits a contig, so
        // a cleared contig is never fetched again.
        if current_chrom_id != Some(region.chrom_id) {
            repository.clear();
            current_chrom_id = Some(region.chrom_id);
        }

        let contig_name = inputs.contigs.entries[region.chrom_id as usize]
            .name
            .clone();
        let reader =
            SegmentMergedReads::new(&segment_files, &contig_name, region.start, region.end)?;

        let outputs = super::stage1_pipeline::with_stage1_chain::<_, RunSummary, PileupCliError, _>(
            reader,
            &args.reference,
            &repository,
            &walker_fetcher,
            baq_cfg,
            proc_cfg,
            walker_cfg,
            stage1.baq_chunk_size,
            n_threads,
            stage1.no_baq,
            &inputs.contigs,
            |ctx| {
                drive_region_into_writer(ctx.walker, &mut writer, region.start, region.end)
                    .map_err(PileupCliError::from)
            },
        )?;

        total_walker.merge(&outputs.result);
        total_filter.merge(&outputs.run_summary.filter_counts);
        if let (Some(acc), Some(region_skip)) = (
            total_baq_skip.as_mut(),
            outputs.run_summary.baq_skip_counts.as_ref(),
        ) {
            acc.merge(region_skip);
        }

        // An upstream read error on any region halts the run. Stop the
        // loop and surface it after the writer is torn down.
        if let Some(e) = outputs.stashed_upstream_error {
            stashed_upstream = Some(e);
            break;
        }
    }

    // 8. Finalise: fsync + atomically rename the .tmp into place.
    let buf_sink = writer.finish().map_err(PileupToPspError::from)?;
    finalise_output(buf_sink, &tmp_path, &args.output)?;

    // 9. Surface a stashed upstream read error (the walker exited cleanly
    //    from its own perspective, so this is the only report site).
    if let Some(e) = stashed_upstream {
        let _ = fs::remove_file(&args.output); // best-effort cleanup
        return Err(PileupCliError::AlignmentInput(e));
    }

    // 10. Stderr run-summary block, totalled across regions.
    print_run_summary(
        &inputs.sample_name,
        &total_filter,
        total_baq_skip.as_ref(),
        &total_walker,
    );

    Ok(())
}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

pub(crate) fn alignment_config_from_args(
    args: &shared_args::Stage1Args,
) -> AlignmentMergedReaderConfig {
    AlignmentMergedReaderConfig {
        min_mapq: if args.min_mapq == 0 {
            None
        } else {
            Some(args.min_mapq)
        },
        min_read_length: if args.min_read_length == 0 {
            None
        } else {
            Some(args.min_read_length)
        },
        drop_qc_fail: !args.keep_qc_fail,
        drop_duplicate: !args.keep_duplicates,
        max_read_mismatch_fraction: if args.max_read_mismatch_fraction == 0.0 {
            None
        } else {
            Some(args.max_read_mismatch_fraction)
        },
        mismatch_bq_floor: args.mismatch_bq_floor,
    }
}

/// F1 mismatch-fraction config for the parallel read-processing stage.
/// Drawn from the same `Stage1Args` fields the reader config uses;
/// these filters moved out of the reader into `read_processor`.
pub(crate) fn read_processing_config_from_args(
    args: &shared_args::Stage1Args,
) -> ReadProcessingConfig {
    ReadProcessingConfig {
        max_read_mismatch_fraction: if args.max_read_mismatch_fraction == 0.0 {
            None
        } else {
            Some(args.max_read_mismatch_fraction)
        },
        mismatch_bq_floor: args.mismatch_bq_floor,
    }
}

pub(crate) fn baq_config_from_args(args: &shared_args::Stage1Args) -> BaqConfig {
    BaqConfig {
        gap_open_prob: args.baq_gap_open_prob,
        gap_extend_prob: args.baq_gap_extend_prob,
        band_half_width: args.baq_band_half_width,
    }
}

pub(crate) fn walker_config_from_args(args: &shared_args::Stage1Args) -> WalkerConfig {
    WalkerConfig {
        max_snp_column_depth: args.max_snp_column_depth,
        max_indel_column_depth: args.max_indel_column_depth,
        max_record_span: args.max_record_span,
        mate_lookup_window: args.mate_lookup_window,
        max_active_reads: args.max_active_reads,
    }
}

/// Append `.tmp` to the file name. Works for any path: `foo.psp` →
/// `foo.psp.tmp`, `out` → `out.tmp`.
fn tmp_path_for(output: &Path) -> PathBuf {
    let mut s = output.as_os_str().to_os_string();
    s.push(".tmp");
    PathBuf::from(s)
}

/// Recover the underlying File from the BufWriter the walker→writer
/// seam returned, fsync it, then atomically rename the tmp file
/// into place.
fn finalise_output(
    buf: BufWriter<File>,
    tmp_path: &Path,
    final_path: &Path,
) -> Result<(), PileupCliError> {
    let file = buf.into_inner().map_err(|e| {
        // `IntoInnerError` carries both the writer and the IO error;
        // we just want the io::Error for the user-facing message.
        PileupCliError::Io(e.into_error())
    })?;
    file.sync_all().map_err(PileupCliError::Io)?;
    fs::rename(tmp_path, final_path).map_err(PileupCliError::Io)?;
    Ok(())
}

/// Build a [`WriterHeader`] from the merged reader's metadata + CLI
/// inputs. Hard-errors when a contig has no `@SQ M5`.
#[allow(clippy::too_many_arguments)]
fn build_writer_header(
    sample: &str,
    fasta_path: &Path,
    contigs: &ContigList,
    alignment_paths: &[PathBuf],
    args: &shared_args::Stage1Args,
    block_target_bytes: usize,
    block_window_bp: u32,
    regions_bed: Option<&Path>,
    n_threads: usize,
) -> Result<WriterHeader, PileupCliError> {
    let mut chromosomes = Vec::with_capacity(contigs.entries.len());
    for entry in &contigs.entries {
        let length =
            u32::try_from(entry.length).map_err(|_| PileupCliError::ContigLengthOverflow {
                name: entry.name.clone(),
                length: entry.length,
            })?;
        let md5 = entry
            .md5
            .map(format_md5_hex)
            .ok_or_else(|| PileupCliError::MissingMd5 {
                contig: entry.name.clone(),
            })?;
        chromosomes.push(ChromosomeEntry {
            name: entry.name.clone(),
            length,
            md5,
        });
    }

    let created = rfc3339_now()
        .parse()
        .map_err(|_| PileupCliError::TimestampFormat)?;

    let input_crams: Vec<String> = alignment_paths.iter().map(|p| basename(p)).collect();
    let input_fasta = basename(fasta_path);
    let mut parameters = effective_parameters(args, block_target_bytes, block_window_bp, n_threads);
    // Record the analysis-regions BED basename when one was supplied, so
    // the `.psp` self-describes that it covers only those regions (the
    // full path is also in `command_line`). `regions_count` is added by
    // the caller once the BED has been resolved to a region set.
    if let Some(bed) = regions_bed {
        parameters.insert(
            "regions_bed".to_string(),
            ParameterValue::String(basename(bed)),
        );
    }

    Ok(WriterHeader {
        format_version: (1, 0),
        sample: sample.to_string(),
        reference: input_fasta.clone(),
        created,
        chromosomes,
        writer: WriterProvenance {
            tool: "pop_var_caller".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            subcommand: "pileup".to_string(),
            input_crams,
            input_fasta,
            command_line: current_command_line(),
            parameters,
        },
    })
}

/// Effective parameter map, written into `WriterProvenance.parameters`
/// so the .psp self-describes the run. Every CLI knob — including the
/// ones left at their defaults — is recorded.
fn effective_parameters(
    args: &shared_args::Stage1Args,
    block_target_bytes: usize,
    block_window_bp: u32,
    n_threads: usize,
) -> BTreeMap<String, ParameterValue> {
    let mut p = BTreeMap::new();
    // CRAM input filters
    p.insert(
        "min_mapq".into(),
        ParameterValue::Integer(args.min_mapq as i64),
    );
    p.insert(
        "min_read_length".into(),
        ParameterValue::Integer(args.min_read_length as i64),
    );
    p.insert(
        "drop_qc_fail".into(),
        ParameterValue::Boolean(!args.keep_qc_fail),
    );
    p.insert(
        "drop_duplicate".into(),
        ParameterValue::Boolean(!args.keep_duplicates),
    );
    p.insert(
        "max_read_mismatch_fraction".into(),
        ParameterValue::Float(args.max_read_mismatch_fraction as f64),
    );
    p.insert(
        "mismatch_bq_floor".into(),
        ParameterValue::Integer(args.mismatch_bq_floor as i64),
    );
    // BAQ
    p.insert("baq_enabled".into(), ParameterValue::Boolean(!args.no_baq));
    p.insert(
        "baq_gap_open_prob".into(),
        ParameterValue::Float(args.baq_gap_open_prob as f64),
    );
    p.insert(
        "baq_gap_extend_prob".into(),
        ParameterValue::Float(args.baq_gap_extend_prob as f64),
    );
    p.insert(
        "baq_band_half_width".into(),
        ParameterValue::Integer(args.baq_band_half_width as i64),
    );
    p.insert(
        "baq_chunk_size".into(),
        ParameterValue::Integer(args.baq_chunk_size as i64),
    );
    // Walker
    p.insert(
        "max_snp_column_depth".into(),
        ParameterValue::Integer(args.max_snp_column_depth as i64),
    );
    p.insert(
        "max_indel_column_depth".into(),
        ParameterValue::Integer(args.max_indel_column_depth as i64),
    );
    p.insert(
        "max_record_span".into(),
        ParameterValue::Integer(args.max_record_span as i64),
    );
    p.insert(
        "mate_lookup_window".into(),
        ParameterValue::Integer(args.mate_lookup_window as i64),
    );
    p.insert(
        "max_active_reads".into(),
        ParameterValue::Integer(args.max_active_reads as i64),
    );
    // Threads — the resolved budget passed in by the caller (explicit
    // `--threads`, else the stage's default; pileup's is
    // `DEFAULT_PILEUP_THREADS`, not all logical cores). Read from the
    // flag, not `rayon`, since the staged pipeline doesn't size a rayon
    // pool.
    p.insert("threads".into(), ParameterValue::Integer(n_threads as i64));
    // PSP writer.
    p.insert(
        "block_target_bytes".into(),
        ParameterValue::Integer(block_target_bytes as i64),
    );
    p.insert(
        "block_window_bp".into(),
        ParameterValue::Integer(block_window_bp as i64),
    );
    p
}

fn print_run_summary(
    sample: &str,
    filter_counts: &FilterCounts,
    baq_skip: Option<&BaqSkipCounts>,
    summary: &RunSummary,
) {
    eprintln!("sample: {sample}");
    eprintln!(
        "filter_counts: unmapped={} secondary={} supplementary={} qc_fail={} duplicate={} \
         low_mapq={} too_short={} high_mismatch_fraction={} bad_cigar={} baq_rejected={}",
        filter_counts.unmapped,
        filter_counts.secondary,
        filter_counts.supplementary,
        filter_counts.qc_fail,
        filter_counts.duplicate,
        filter_counts.low_mapq,
        filter_counts.too_short,
        filter_counts.high_mismatch_fraction,
        filter_counts.bad_cigar,
        filter_counts.baq_rejected,
    );
    match baq_skip {
        Some(b) => eprintln!("baq: total_skipped={}", b.total),
        None => eprintln!("baq: disabled"),
    }
    eprintln!(
        "walker: reads_admitted={} records_emitted={} record_widen_events={} \
         mate_overlap_positions={} active_reads_high_water={} mate_lookup_evictions={} \
         column_depth_truncations={}",
        summary.reads_admitted,
        summary.records_emitted,
        summary.record_widen_events,
        summary.mate_overlap_positions,
        summary.active_reads_high_water,
        summary.mate_lookup_evictions,
        summary.column_depth_truncations,
    );
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::{
        DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH,
        DEFAULT_MISMATCH_BQ_FLOOR,
    };
    use crate::baq::{
        SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH, SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
        SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
    };
    use crate::pileup::per_sample::baq_stream::DEFAULT_BAQ_CHUNK_SIZE;
    use crate::pileup::walker::{
        DEFAULT_MATE_LOOKUP_WINDOW, DEFAULT_MAX_ACTIVE_READS, DEFAULT_MAX_INDEL_COLUMN_DEPTH,
        DEFAULT_MAX_RECORD_SPAN, DEFAULT_MAX_SNP_COLUMN_DEPTH,
    };
    use crate::psp::writer::TARGET_BLOCK_BYTES;

    fn default_args(
        reference: PathBuf,
        output: PathBuf,
        alignment_files: Vec<PathBuf>,
    ) -> PileupArgs {
        PileupArgs {
            reference,
            output,
            regions: None,
            build_map_file_index: false,
            threads: None,
            alignment_files,
            block_target_bytes: TARGET_BLOCK_BYTES,
            block_window_bp: DEFAULT_BLOCK_WINDOW_BP,
            stage1: shared_args::Stage1Args {
                min_mapq: DEFAULT_MIN_MAPQ,
                no_baq: false,
                min_read_length: DEFAULT_MIN_READ_LENGTH,
                keep_qc_fail: false,
                keep_duplicates: false,
                max_read_mismatch_fraction: DEFAULT_MAX_READ_MISMATCH_FRACTION,
                mismatch_bq_floor: DEFAULT_MISMATCH_BQ_FLOOR,
                baq_gap_open_prob: SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
                baq_gap_extend_prob: SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
                baq_band_half_width: SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
                baq_chunk_size: DEFAULT_BAQ_CHUNK_SIZE,
                max_snp_column_depth: DEFAULT_MAX_SNP_COLUMN_DEPTH,
                max_indel_column_depth: DEFAULT_MAX_INDEL_COLUMN_DEPTH,
                max_record_span: DEFAULT_MAX_RECORD_SPAN,
                mate_lookup_window: DEFAULT_MATE_LOOKUP_WINDOW,
                max_active_reads: DEFAULT_MAX_ACTIVE_READS,
            },
        }
    }

    #[test]
    fn parse_mismatch_fraction_accepts_in_range() {
        assert_eq!(parse_mismatch_fraction("0.0").unwrap(), 0.0);
        assert_eq!(parse_mismatch_fraction("0.1").unwrap(), 0.1);
        assert_eq!(parse_mismatch_fraction("1.0").unwrap(), 1.0);
    }

    #[test]
    fn parse_mismatch_fraction_rejects_out_of_range() {
        assert!(parse_mismatch_fraction("-0.001").is_err());
        assert!(parse_mismatch_fraction("1.0001").is_err());
        assert!(parse_mismatch_fraction("nan").is_err());
        assert!(parse_mismatch_fraction("inf").is_err());
        assert!(parse_mismatch_fraction("hello").is_err());
    }

    #[test]
    fn tmp_path_appends_suffix() {
        assert_eq!(
            tmp_path_for(Path::new("/tmp/x.psp")),
            PathBuf::from("/tmp/x.psp.tmp")
        );
        assert_eq!(tmp_path_for(Path::new("out")), PathBuf::from("out.tmp"));
    }

    #[test]
    fn effective_parameters_records_every_knob() {
        let args = default_args(PathBuf::from("r.fa"), PathBuf::from("o.psp"), vec![]);
        let p = effective_parameters(
            &args.stage1,
            args.block_target_bytes,
            args.block_window_bp,
            4,
        );
        // Sanity: every knob shows up.
        for key in &[
            "min_mapq",
            "min_read_length",
            "drop_qc_fail",
            "drop_duplicate",
            "max_read_mismatch_fraction",
            "mismatch_bq_floor",
            "baq_enabled",
            "baq_gap_open_prob",
            "baq_gap_extend_prob",
            "baq_band_half_width",
            "baq_chunk_size",
            "max_snp_column_depth",
            "max_indel_column_depth",
            "max_record_span",
            "mate_lookup_window",
            "max_active_reads",
            "threads",
            "block_target_bytes",
            "block_window_bp",
        ] {
            assert!(p.contains_key(*key), "missing parameter: {key}");
        }
        // Default booleans land where we expect.
        assert_eq!(p.get("drop_qc_fail"), Some(&ParameterValue::Boolean(true)));
        assert_eq!(p.get("baq_enabled"), Some(&ParameterValue::Boolean(true)));
        // Block-target default tracks `TARGET_BLOCK_BYTES`.
        assert_eq!(
            p.get("block_target_bytes"),
            Some(&ParameterValue::Integer(TARGET_BLOCK_BYTES as i64))
        );
        // Block-window default tracks `DEFAULT_BLOCK_WINDOW_BP`.
        assert_eq!(
            p.get("block_window_bp"),
            Some(&ParameterValue::Integer(DEFAULT_BLOCK_WINDOW_BP as i64))
        );
    }

    #[test]
    fn build_writer_header_errors_on_missing_md5() {
        use crate::fasta::{ContigEntry, ContigList};
        let contigs = ContigList {
            entries: vec![ContigEntry {
                name: "chr1".into(),
                length: 100,
                md5: None,
            }],
        };
        let args = default_args(PathBuf::from("r.fa"), PathBuf::from("o.psp"), vec![]);
        let err = build_writer_header(
            "s",
            Path::new("r.fa"),
            &contigs,
            &[],
            &args.stage1,
            args.block_target_bytes,
            args.block_window_bp,
            None,
            4,
        )
        .expect_err("must error");
        assert!(matches!(err, PileupCliError::MissingMd5 { ref contig } if contig == "chr1"));
    }

    #[test]
    fn build_writer_header_strips_paths_to_basenames() {
        use crate::fasta::{ContigEntry, ContigList};
        let md5 = [0u8; 16];
        let contigs = ContigList {
            entries: vec![ContigEntry {
                name: "chr1".into(),
                length: 100,
                md5: Some(md5),
            }],
        };
        let args = default_args(
            PathBuf::from("/data/ref/grch38.fa"),
            PathBuf::from("o.psp"),
            vec![],
        );
        let alignment_files = [
            PathBuf::from("/scratch/run1/sample.cram"),
            PathBuf::from("./run2/sample2.cram"),
        ];
        let header = build_writer_header(
            "NA12878",
            Path::new("/data/ref/grch38.fa"),
            &contigs,
            &alignment_files,
            &args.stage1,
            args.block_target_bytes,
            args.block_window_bp,
            None,
            4,
        )
        .expect("build_writer_header");
        assert_eq!(header.writer.input_fasta, "grch38.fa");
        assert_eq!(
            header.writer.input_crams,
            vec!["sample.cram", "sample2.cram"]
        );
        // MD5 round-trips as 32 hex chars.
        assert_eq!(header.chromosomes[0].md5.len(), 32);
        // Subcommand identifier matches the binary surface.
        assert_eq!(header.writer.subcommand, "pileup");
        assert_eq!(header.writer.tool, "pop_var_caller");
    }
}
