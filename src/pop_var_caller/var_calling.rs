//! `pop_var_caller var-calling` — cohort `.psp` → multi-sample VCF.
//!
//! Wires Stages 3–6 of the pipeline:
//!
//!   PspReader(s) → PerPositionMerger → [DustFilter] → VariantGrouper
//!   → PerGroupMerger → PosteriorEngine → CohortVcfWriter
//!
//! The DUST filter is bypassed when `--no-complexity-filter` is set;
//! both branches converge on the same downstream chain via a shared
//! [`Box<dyn Iterator>`] adapter that lifts the upstream error into
//! [`GrouperError`].
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

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use clap::Args;
use thiserror::Error;

use crate::per_sample_pileup::psp::{PspReadError, PspReader};
use crate::per_sample_pileup::ref_fetcher::SyncRefFetcher;
use crate::pop_var_caller::cli::parsers;
use crate::pop_var_caller::contamination_artefact::{
    ContaminationArtefact, ContaminationArtefactError,
};
use crate::var_calling::contamination_estimation::ContaminationEstimates;
use crate::var_calling::dust_filter::{
    DEFAULT_DUST_THRESHOLD, DEFAULT_DUST_WINDOW, DustFilter, DustFilterConfig, DustFilterError,
};
use crate::var_calling::per_group_merger::{
    DEFAULT_BATCH_SIZE, DEFAULT_MAX_ALLELES_PER_RECORD, DEFAULT_PLOIDY, PerGroupMerger,
    PerGroupMergerConfig, PerGroupMergerError, SharedRefFetcher,
};
use crate::var_calling::per_position_merger::{
    PerPositionMerger, PerPositionMergerError, PerPositionPileups, check_chromosome_agreement,
};
use crate::var_calling::posterior_engine::{
    DEFAULT_COMPOUND_ALT_PSEUDOCOUNT, DEFAULT_CONVERGENCE_THRESHOLD,
    DEFAULT_INBREEDING_COEFFICIENT, DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MAX_GQ_PHRED,
    DEFAULT_MAX_ITERATIONS, DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT, PosteriorEngine,
    PosteriorEngineConfig, PosteriorEngineConfigError, PosteriorEngineError,
};
use crate::var_calling::variant_grouping::{
    DEFAULT_MAX_VARIANT_GROUP_SPAN, GrouperConfig, GrouperConfigError, GrouperError, VariantGrouper,
};
use crate::var_calling::vcf_writer::{
    CohortMetadata, CohortVcfWriter, DEFAULT_EMIT_GP, VcfWriteError, WriterConfig, tmp_path_for,
};

// ---------------------------------------------------------------------
// CLI surface
// ---------------------------------------------------------------------

/// Arguments accepted by `pop_var_caller var-calling`.
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

    /// Cohort-wide ploidy.
    #[arg(long, default_value_t = DEFAULT_PLOIDY, value_parser = parsers::parse_ploidy)]
    pub ploidy: u8,

    /// One or more cohort `.psp` files.
    #[arg(required = true)]
    pub psp_files: Vec<PathBuf>,

    // ===== Advanced — DUST filter ==============================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_DUST_WINDOW,
        value_parser = parsers::parse_dust_window,
        help_heading = "Advanced — DUST filter",
    )]
    pub complexity_window: u32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_DUST_THRESHOLD,
        value_parser = parsers::parse_dust_threshold,
        help_heading = "Advanced — DUST filter",
    )]
    pub complexity_threshold: u32,

    // ===== Advanced — Variant grouping =========================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_VARIANT_GROUP_SPAN,
        value_parser = parsers::parse_var_group_max_span,
        help_heading = "Advanced — Variant grouping",
    )]
    pub var_group_max_span: u32,

    // ===== Advanced — Per-group merger =========================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_ALLELES_PER_RECORD,
        value_parser = parsers::parse_max_alleles,
        help_heading = "Advanced — Per-group merger",
    )]
    pub max_alleles_per_var: usize,

    // ===== Advanced — Posterior engine =========================
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_INBREEDING_COEFFICIENT,
        value_parser = parsers::parse_inbreeding_coefficient,
        help_heading = "Advanced — Posterior engine",
    )]
    pub inbreeding_coefficient: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_CONVERGENCE_THRESHOLD,
        value_parser = parsers::parse_em_convergence_threshold,
        help_heading = "Advanced — Posterior engine",
    )]
    pub em_convergence_threshold: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_ITERATIONS,
        value_parser = parsers::parse_em_max_iterations,
        help_heading = "Advanced — Posterior engine",
    )]
    pub em_max_iterations: u32,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_REF_PSEUDOCOUNT,
        value_parser = parsers::parse_ref_pseudocount,
        help_heading = "Advanced — Posterior engine",
    )]
    pub ref_pseudocount: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_SNP_ALT_PSEUDOCOUNT,
        value_parser = parsers::parse_snp_alt_pseudocount,
        help_heading = "Advanced — Posterior engine",
    )]
    pub snp_alt_pseudocount: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_INDEL_ALT_PSEUDOCOUNT,
        value_parser = parsers::parse_indel_alt_pseudocount,
        help_heading = "Advanced — Posterior engine",
    )]
    pub indel_alt_pseudocount: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_COMPOUND_ALT_PSEUDOCOUNT,
        value_parser = parsers::parse_compound_alt_pseudocount,
        help_heading = "Advanced — Posterior engine",
    )]
    pub compound_alt_pseudocount: f64,

    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MAX_GQ_PHRED,
        value_parser = parsers::parse_max_gq_phred,
        help_heading = "Advanced — Posterior engine",
    )]
    pub max_gq_phred: f64,

    // ===== Advanced — VCF writer ===============================
    /// Emit `GP` (genotype posteriors) `FORMAT` per sample. Off by
    /// default — `GP` is `Number=G`, so the per-sample cell grows as
    /// `(ploidy + n_alleles − 1) choose ploidy` (21 floats at
    /// ploidy=2, n_alleles=6). Opt in when posteriors are wanted on
    /// disk; a 1000-sample × 100 K-variant cohort with `emit_gp` adds
    /// ~8 GiB to the VCF.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_EMIT_GP,
        help_heading = "Advanced — VCF writer",
    )]
    pub emit_gp: bool,
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
/// 8. Stream records into the [`CohortVcfWriter`]; finalise via
///    atomic tmp+rename.
/// 9. Print a one-shot stderr run-summary block.
pub fn run_var_calling(args: &VarCallingArgs) -> Result<(), VarCallingCliError> {
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
        let buf = BufReader::with_capacity(64 * 1024, file);
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
    let dust_cfg = DustFilterConfig::new(args.complexity_window, args.complexity_threshold)?;
    let grouper_cfg = GrouperConfig::new(args.var_group_max_span)?;
    let per_group_cfg =
        PerGroupMergerConfig::new(args.ploidy, args.max_alleles_per_var, DEFAULT_BATCH_SIZE)?;
    let mut posterior_cfg = PosteriorEngineConfig::new(
        args.em_convergence_threshold,
        args.em_max_iterations,
        args.ref_pseudocount,
        args.snp_alt_pseudocount,
        args.indel_alt_pseudocount,
        args.compound_alt_pseudocount,
        args.inbreeding_coefficient,
        args.max_gq_phred,
    )?;

    // 6. Load contamination if supplied.
    posterior_cfg.contamination =
        load_contamination(args.contamination_estimates.as_deref(), &sample_names)?;

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
    let writer_cfg = WriterConfig::new(args.output.clone()).with_emit_gp(args.emit_gp);

    // 9. Wire the pipeline. The k-way merger sits on top of the
    //    record iterators borrowed from each reader.
    let record_iters: Vec<_> = readers.iter_mut().map(|r| r.records()).collect();
    let merger = PerPositionMerger::new(record_iters, sample_names.clone(), chromosomes.clone())?;

    // DUST filter is optional; both branches converge on a
    // `Box<dyn Iterator<Item = Result<_, GrouperError>>>` so the
    // grouper sees a single type either way. The grouper is generic
    // over upstream error via `GrouperError: From<E>` — see the
    // commit that introduced the generalisation.
    // DustFilter wants `F: RefSeqFetcher` by value. `Arc<dyn ...>`
    // does not itself impl the trait (no blanket impl for Arc), but
    // `&dyn RefSeqFetcher` does via the existing `impl<T: RefSeqFetcher
    // + ?Sized> RefSeqFetcher for &T` blanket. Passing `&*fetcher`
    // satisfies the bound and ties the filter's lifetime to `fetcher`,
    // which lives for the whole pipeline scope.
    let upstream_for_grouper: Box<
        dyn Iterator<Item = Result<PerPositionPileups, GrouperError>> + '_,
    > = if args.no_complexity_filter {
        Box::new(merger.map(|r| r.map_err(GrouperError::from)))
    } else {
        let dust = DustFilter::new(merger, &*fetcher, chromosomes.clone(), dust_cfg);
        Box::new(dust.map(|r| r.map_err(GrouperError::from)))
    };

    let grouper = VariantGrouper::with_config(upstream_for_grouper, grouper_cfg);
    let per_group = PerGroupMerger::with_config(grouper, fetcher.clone(), per_group_cfg);
    let posterior = PosteriorEngine::with_config(per_group, posterior_cfg);

    // 10. Stream records into the writer; tmp-then-rename on success.
    //     On any driver-loop failure we best-effort remove the half-
    //     written `<output>.tmp` so failed runs don't leave stale tmp
    //     files behind — parity with `run_pileup` and
    //     `run_var_calling_from_bam`. The writer's `finish()` consumes
    //     `self`, so a `?` short-circuit in the loop never reaches the
    //     rename and the tmp would otherwise linger.
    let tmp_path = tmp_path_for(&args.output);
    let mut writer = CohortVcfWriter::new(metadata, writer_cfg)?;
    let mut records_written: u64 = 0;
    let drive_result: Result<(), VarCallingCliError> = (|| {
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

fn basename(p: &Path) -> String {
    p.file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| p.to_string_lossy().into_owned())
}

fn current_command_line() -> String {
    std::env::args_os()
        .map(|a: OsString| a.to_string_lossy().into_owned())
        .collect::<Vec<_>>()
        .join(" ")
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
    fn basename_strips_directory() {
        assert_eq!(basename(Path::new("/data/grch38.fa")), "grch38.fa");
        assert_eq!(basename(Path::new("ref.fa")), "ref.fa");
    }

    #[test]
    fn load_contamination_returns_none_without_path() {
        let names = vec!["NA12878".to_string()];
        let result = load_contamination(None, &names).unwrap();
        assert!(result.is_none());
    }
}
