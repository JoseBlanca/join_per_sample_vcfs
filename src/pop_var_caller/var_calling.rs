//! `pop_var_caller var-calling` — cohort `.psp` → multi-sample VCF.
//!
//! Wires Stages 3–6 of the pipeline:
//!
//!   PspReader(s) → PerPositionMerger → DustFilter → VariantGrouper
//!   → PerGroupMerger → PosteriorEngine → CohortVcfWriter
//!
//! The DUST filter is bypassed when `--no-complexity-filter` is set;
//! both branches converge on the same downstream chain via a
//! `Box<dyn Iterator>` adapter that lifts the upstream error into
//! [`GrouperError`] — the wiring lives in
//! [`crate::pop_var_caller::cohort_driver::drive_cohort_pipeline`].
//!
//! Contamination plumbing:
//!
//! - With `--contamination-estimates <FILE>`: the artefact is loaded
//!   and reconciled against the cohort sample names (extras are
//!   tolerated, absences are not), then handed to the posterior engine
//!   via [`PosteriorEngineConfig::contamination`].
//! - Without: the engine runs in "no contamination" mode and the
//!   `contamination` field stays `None`.
//!
//! Plan: `doc/devel/implementation_plans/pop_var_caller_cohort_cli.md`
//! §"Subcommand: var-calling".

use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use clap::Args;
use thiserror::Error;

use crate::per_sample_pileup::psp::{PspReadError, PspReader};
use crate::per_sample_pileup::ref_fetcher::SyncRefFetcher;
use crate::pop_var_caller::cohort_driver::{CohortPipelineParams, drive_cohort_pipeline};
use crate::pop_var_caller::common::{DEFAULT_BUFFERED_IO_CAPACITY, basename, current_command_line};
use crate::pop_var_caller::contamination_artefact::{
    ContaminationArtefact, ContaminationArtefactError,
};
use crate::var_calling::contamination_estimation::ContaminationEstimates;
use crate::var_calling::dust_filter::{DustFilterConfig, DustFilterError};
use crate::var_calling::per_group_merger::{
    DEFAULT_BATCH_SIZE, PerGroupMergerConfig, PerGroupMergerError, SharedRefFetcher,
};
use crate::var_calling::per_position_merger::{
    PerPositionMerger, PerPositionMergerError, check_chromosome_agreement,
};
use crate::var_calling::posterior_engine::{
    PosteriorEngineConfig, PosteriorEngineConfigError, PosteriorEngineError,
};
use crate::var_calling::variant_grouping::{GrouperConfig, GrouperConfigError, GrouperError};
use crate::var_calling::vcf_writer::{CohortMetadata, VcfWriteError, WriterConfig};

// ---------------------------------------------------------------------
// CLI surface
// ---------------------------------------------------------------------

/// Arguments accepted by `pop_var_caller var-calling`.
///
/// Cohort-pipeline knobs (DUST, grouper, per-group merger,
/// posterior engine, ploidy, VCF writer) live in the flattened
/// [`CohortPipelineArgs`](crate::pop_var_caller::cli::shared_args::CohortPipelineArgs)
/// sub-struct so the `var-calling-from-bam` subcommand reuses the
/// same surface — M10 from the 2026-05-19 cohort CLI review.
#[derive(Debug, Args, Clone)]
pub struct VarCallingArgs {
    // ===== Common (visible in `-h`) ===========================
    /// Reference FASTA used to produce the `.psp` files. Basename
    /// must match the `reference` field every `.psp` header carries.
    #[arg(long)]
    pub reference: PathBuf,

    /// Output VCF path. Extension picks the sink kind:
    /// `.vcf.gz` / `.vcf.bgz` → bgzf, anything else → plain text.
    /// Written via atomic `<output>.tmp` → `fs::rename`.
    #[arg(long)]
    pub output: PathBuf,

    /// Worker threads for the rayon pool. If omitted, rayon picks
    /// the default (all logical cores).
    #[arg(long)]
    pub threads: Option<usize>,

    /// Contamination-estimates artefact produced by
    /// `pop_var_caller estimate-contamination`. When omitted, every
    /// sample's `c_s` is treated as `0` (no contamination correction).
    #[arg(long)]
    pub contamination_estimates: Option<PathBuf>,

    /// Skip the low-complexity (sdust) filter entirely.
    #[arg(long)]
    pub no_complexity_filter: bool,

    /// One or more cohort `.psp` files.
    #[arg(required = true)]
    pub psp_files: Vec<PathBuf>,

    // ===== Cohort pipeline (shared with var-calling-from-bam) ==
    #[command(flatten)]
    pub cohort: crate::pop_var_caller::cli::shared_args::CohortPipelineArgs,
}

// ---------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------

/// Errors surfaced by [`run_var_calling`]. Wraps every upstream module
/// error so the CLI's `format_error_chain` can render the cause.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum VarCallingCliError {
    #[error("io: {0}")]
    Io(#[from] io::Error),

    #[error("psp reader: {0}")]
    PspReader(#[from] PspReadError),

    #[error("merger: {0}")]
    Merger(#[from] PerPositionMergerError),

    #[error("dust filter: {0}")]
    Dust(#[from] DustFilterError),

    #[error("grouper: {0}")]
    Grouper(#[from] GrouperError),

    #[error("per-group merger: {0}")]
    PerGroup(#[from] PerGroupMergerError),

    #[error("posterior engine: {0}")]
    Posterior(#[from] PosteriorEngineError),

    #[error("vcf writer: {0}")]
    Vcf(#[from] VcfWriteError),

    #[error("contamination artefact: {0}")]
    ContamArtefact(#[from] ContaminationArtefactError),

    #[error("grouper config: {0}")]
    GrouperConfig(#[from] GrouperConfigError),

    #[error("per-group merger config: {0}")]
    PerGroupConfig(#[from] crate::var_calling::per_group_merger::PerGroupMergerConfigError),

    #[error("posterior engine config: {0}")]
    PosteriorConfig(#[from] PosteriorEngineConfigError),

    // `DustFilterConfig::new` also returns `DustFilterError`, which
    // funnels through the `Dust` variant above — no separate
    // `DustConfig` slot needed.
    /// A `.psp` file's `reference` field doesn't match the basename
    /// of `--reference`. Surfaces before any record is read.
    #[error(
        "psp {psp}: reference mismatch — header has `{psp_ref}`, CLI passed `{supplied_ref}` \
         (basename comparison)"
    )]
    ReferenceMismatch {
        psp: PathBuf,
        psp_ref: String,
        supplied_ref: String,
    },

    /// `rayon::ThreadPoolBuilder::build_global()` already ran in this
    /// process. The binary calls it at most once.
    #[error("rayon thread pool already initialised — refusing to override")]
    RayonAlreadyConfigured,
}

// ---------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------

/// Run the cohort var-calling pipeline and emit a multi-sample VCF.
///
/// Pipeline (per the cohort CLI plan §"Subcommand: var-calling"):
///
/// 1. Size the global rayon pool when `--threads` is set.
/// 2. Open every `.psp` with [`PspReader`].
/// 3. Cross-check `header.reference` against the basename of
///    `--reference`. **Basename-only — the FASTA on disk is not
///    content-hashed against the `.psp` per-contig MD5s.** A future
///    slice may wire that check; see the v1 contract note below.
/// 4. Cross-check chromosomes across readers
///    ([`check_chromosome_agreement`]). This compares per-contig MD5s
///    *between readers* only, not against the FASTA — so two `.psp`s
///    produced from the same reference agree on contents, but the
///    supplied `--reference` is trusted to match by basename only.
/// 5. Load `--contamination-estimates` if supplied and reconcile
///    sample names; absence → `c_s = 0` for every sample.
/// 6. Build + validate every per-stage config from the CLI args.
/// 7. Wire the pipeline (merger → optional DUST → grouper →
///    per-group merger → posterior engine).
/// 8. Stream records into the
///    [`CohortVcfWriter`](crate::var_calling::vcf_writer::CohortVcfWriter); finalise via
///    atomic tmp+rename.
/// 9. Print a one-shot stderr run-summary block.
pub fn run_var_calling(args: &VarCallingArgs) -> Result<(), VarCallingCliError> {
    let cohort = &args.cohort;

    // 1. Rayon pool.
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|_| VarCallingCliError::RayonAlreadyConfigured)?;
    }

    // 2. Open every .psp.
    let mut readers: Vec<PspReader<BufReader<File>>> = Vec::with_capacity(args.psp_files.len());
    for path in &args.psp_files {
        let file = File::open(path).map_err(VarCallingCliError::Io)?;
        let buf = BufReader::with_capacity(DEFAULT_BUFFERED_IO_CAPACITY, file);
        let reader = PspReader::new(buf)?;
        readers.push(reader);
    }

    // 3. Cross-check references. Basename comparison only; the FASTA
    //    bytes on disk are not hashed against the .psp per-contig MD5s.
    //    See the function docstring's "v1 contract" note.
    let supplied_ref = basename(&args.reference);
    for (path, reader) in args.psp_files.iter().zip(readers.iter()) {
        if reader.header().reference != supplied_ref {
            return Err(VarCallingCliError::ReferenceMismatch {
                psp: path.clone(),
                psp_ref: reader.header().reference.clone(),
                supplied_ref: supplied_ref.clone(),
            });
        }
    }

    // 4. Sample names + chromosome agreement.
    let sample_names: Vec<String> = readers.iter().map(|r| r.header().sample.clone()).collect();
    let chromosomes = check_chromosome_agreement(&readers)?;

    // 5. Build + validate every per-stage config.
    let dust_cfg = DustFilterConfig::new(cohort.complexity_window, cohort.complexity_threshold)?;
    let grouper_cfg = GrouperConfig::new(cohort.var_group_max_span)?;
    let per_group_cfg = PerGroupMergerConfig::new(
        cohort.ploidy,
        cohort.max_alleles_per_var,
        DEFAULT_BATCH_SIZE,
    )?;
    // 5. Build the posterior-engine config via the named-setter
    //    builder chain. Every setter validates and returns
    //    `Result<Self, _>` so a swap or out-of-range value surfaces
    //    a typed error at construction time; the `contamination`
    //    field is private and reachable only via
    //    `with_contamination(...)`.
    // 6. Load contamination if supplied; threaded into the engine
    //    config through the same `with_contamination` setter.
    let posterior_cfg = PosteriorEngineConfig::new()
        .with_convergence_threshold(cohort.em_convergence_threshold)?
        .with_max_iterations(cohort.em_max_iterations)?
        .with_ref_pseudocount(cohort.ref_pseudocount)?
        .with_snp_alt_pseudocount(cohort.snp_alt_pseudocount)?
        .with_indel_alt_pseudocount(cohort.indel_alt_pseudocount)?
        .with_compound_alt_pseudocount(cohort.compound_alt_pseudocount)?
        .with_fixation_index_default(cohort.inbreeding_coefficient)?
        .with_max_gq_phred(cohort.max_gq_phred)?
        .with_contamination(load_contamination(
            args.contamination_estimates.as_deref(),
            &sample_names,
        )?)?;

    // 7. Build the reference fetcher (shared between DUST + per-group
    //    merger). The .psp header's chrom table is the source of truth.
    let fetcher_concrete = SyncRefFetcher::new(&args.reference, contigs_from_parsed(&chromosomes))
        .map_err(VarCallingCliError::Io)?;
    let fetcher: SharedRefFetcher = Arc::new(fetcher_concrete);

    // 8. Cohort metadata for the writer header. Sourced from the
    //    merger's view of the cohort, not from the CLI argument order.
    let metadata = CohortMetadata {
        sample_names: sample_names.clone(),
        contigs: chromosomes.clone(),
        tool_string: format!("pop_var_caller {}", env!("CARGO_PKG_VERSION")),
        command_line: current_command_line(),
    };
    let writer_cfg = WriterConfig::new(args.output.clone()).with_emit_gp(cohort.emit_gp);

    // 9. Wire the pipeline. The k-way merger sits on top of the
    //    record iterators borrowed from each reader.
    let record_iters: Vec<_> = readers.iter_mut().map(|r| r.records()).collect();
    let merger = PerPositionMerger::new(record_iters, sample_names.clone(), chromosomes.clone())?;

    // 10. Drive the cohort pipeline (DUST → grouper → per-group
    //     merger → posterior engine → VCF writer) via the shared
    //     helper. M11 from the 2026-05-19 review — the wiring used
    //     to be duplicated verbatim against `var_calling_from_bam`'s
    //     `run_cohort_pipeline_for_single_sample` helper.
    let pipeline_params = CohortPipelineParams {
        no_complexity_filter: args.no_complexity_filter,
        dust_cfg,
        grouper_cfg,
        per_group_cfg,
        posterior_cfg,
        fetcher,
        chromosomes,
    };
    let records_written = drive_cohort_pipeline::<_, VarCallingCliError>(
        merger,
        pipeline_params,
        &args.output,
        metadata,
        writer_cfg,
    )?;

    // 11. Stderr summary.
    print_run_summary(&sample_names, records_written);
    Ok(())
}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

/// Load and reconcile a contamination-estimates artefact against the
/// cohort sample order. Returns `None` if `path` is `None`.
fn load_contamination(
    path: Option<&Path>,
    sample_names: &[String],
) -> Result<Option<ContaminationEstimates>, VarCallingCliError> {
    let path = match path {
        Some(p) => p,
        None => return Ok(None),
    };
    let artefact = ContaminationArtefact::read(path)?;
    let refs: Vec<&str> = sample_names.iter().map(String::as_str).collect();
    Ok(Some(artefact.to_estimates_for_samples(&refs)?))
}

/// Convert the merger's chromosome table into the `ContigList` the
/// [`SyncRefFetcher`] constructor expects. `ContigList` lives in the
/// CRAM input module — we re-shape the .psp `ParsedChromosome` into
/// it. Md5 round-trips through hex when present.
fn contigs_from_parsed(
    chromosomes: &[crate::per_sample_pileup::psp::header::ParsedChromosome],
) -> crate::per_sample_pileup::cram_input::ContigList {
    use crate::per_sample_pileup::cram_input::{ContigEntry, ContigList};
    let entries = chromosomes
        .iter()
        .map(|c| ContigEntry {
            name: c.name.clone(),
            length: c.length as u64,
            md5: md5_hex_to_bytes(&c.md5),
        })
        .collect();
    ContigList { entries }
}

/// Hex (32-char lowercase) → 16-byte MD5. Returns `None` when the
/// .psp header carried no md5 string (empty); the .psp writer
/// hard-errors on missing md5 so this branch should not trigger in
/// production, but the conversion is defensive.
fn md5_hex_to_bytes(hex: &str) -> Option<[u8; 16]> {
    if hex.len() != 32 {
        return None;
    }
    let mut out = [0u8; 16];
    for (i, byte) in out.iter_mut().enumerate() {
        let high = hex_digit(hex.as_bytes()[i * 2])?;
        let low = hex_digit(hex.as_bytes()[i * 2 + 1])?;
        *byte = (high << 4) | low;
    }
    Some(out)
}

fn hex_digit(c: u8) -> Option<u8> {
    match c {
        b'0'..=b'9' => Some(c - b'0'),
        b'a'..=b'f' => Some(c - b'a' + 10),
        b'A'..=b'F' => Some(c - b'A' + 10),
        _ => None,
    }
}

fn print_run_summary(sample_names: &[String], records_written: u64) {
    eprintln!(
        "var-calling: n_samples={} records_emitted={}",
        sample_names.len(),
        records_written,
    );
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn md5_hex_round_trips() {
        let bytes = [
            0x6a, 0xef, 0x89, 0x7c, 0x3d, 0x6f, 0xf0, 0xc7, 0x8a, 0xff, 0x06, 0xac, 0x18, 0x91,
            0x78, 0xdd,
        ];
        let hex = "6aef897c3d6ff0c78aff06ac189178dd";
        assert_eq!(md5_hex_to_bytes(hex), Some(bytes));
    }

    #[test]
    fn md5_hex_rejects_wrong_length() {
        assert_eq!(md5_hex_to_bytes(""), None);
        assert_eq!(md5_hex_to_bytes("6aef"), None);
        assert_eq!(md5_hex_to_bytes(&"a".repeat(33)), None);
    }

    #[test]
    fn md5_hex_rejects_non_hex() {
        assert_eq!(md5_hex_to_bytes("zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"), None);
        // 31 valid + 1 invalid.
        let mut bad = "6aef897c3d6ff0c78aff06ac189178d".to_string();
        bad.push('Z');
        assert_eq!(md5_hex_to_bytes(&bad), None);
    }

    #[test]
    fn load_contamination_returns_none_without_path() {
        let names = vec!["NA12878".to_string()];
        let result = load_contamination(None, &names).unwrap();
        assert!(result.is_none());
    }
}
