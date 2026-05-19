//! `pop_var_caller estimate-contamination` — Stage-6 side-pass CLI.
//!
//! Takes one or more `.psp` files (one per cohort sample) plus an
//! optional `--batch-assignment` TSV, runs the contamination
//! side-pass, and writes a [`ContaminationArtefact`] TOML to disk.
//!
//! The produced artefact is the input `var-calling` consumes via
//! `--contamination-estimates`. The two subcommands share the
//! artefact format defined in
//! [`crate::pop_var_caller::contamination_artifact`].
//!
//! **Convergence-mode only in v1.** The engine's `StoppingMode` enum
//! also supports `FixedSites { num_sites }`, but v1's CLI surface does
//! not expose `--num-sites` (per the cohort CLI plan's flag list).
//! Library callers that want fixed-N can use
//! [`crate::var_calling::contamination_estimation::estimate_contamination`]
//! directly.
//!
//! Plan: `doc/devel/implementation_plans/pop_var_caller_cohort_cli.md` §"Subcommand: estimate-contamination".

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};

use clap::Args;
use thiserror::Error;
use toml::value::Datetime;

use crate::per_sample_pileup::psp::{PspReadError, PspReader};
use crate::pop_var_caller::batch_assignment::{BatchAssignment, BatchAssignmentError};
use crate::pop_var_caller::cli::parsers;
use crate::pop_var_caller::contamination_artifact::{
    BatchEntry, ContaminationArtefact, ContaminationArtefactError, Provenance, ProvenanceInputs,
    SampleEntry,
};
use crate::var_calling::contamination_estimation::{
    ContaminationEstimateSource, ContaminationEstimates, ContaminationEstimationConfig,
    ContaminationEstimationError, DEFAULT_BLOCK_SIZE, DEFAULT_C_S_INIT,
    DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MIN_BATCH_SIZE_FOR_CONTAMINATION,
    DEFAULT_MIN_COHORT_MINOR_COUNT, DEFAULT_MIN_COHORT_MINOR_FRACTION, DEFAULT_MIN_DEPTH,
    DEFAULT_MIN_MAJOR_FRACTION, DEFAULT_Q_B_INIT_PER_CLASS, DEFAULT_REF_PSEUDOCOUNT,
    DEFAULT_SNP_ALT_PSEUDOCOUNT, DEFAULT_STABILITY_BLOCKS, DEFAULT_STABILITY_TOLERANCE,
    StoppingMode, estimate_contamination,
};
use crate::var_calling::per_position_merger::{
    PerPositionMerger, PerPositionMergerError, check_chromosome_agreement,
};

// ---------------------------------------------------------------------
// CLI surface
// ---------------------------------------------------------------------

/// Arguments accepted by `pop_var_caller estimate-contamination`.
#[derive(Debug, Args, Clone)]
pub struct EstimateContaminationArgs {
    // ===== Common (visible in `-h`) ===========================
    /// Reference FASTA used to produce the `.psp` files. Basename
    /// must match the `reference` field every `.psp` header carries.
    #[arg(long)]
    pub reference: PathBuf,

    /// Output TOML path. Atomic tmp-then-rename — see
    /// [`ContaminationArtefact::write`].
    #[arg(long)]
    pub output: PathBuf,

    /// Worker threads for the rayon pool. If omitted, rayon picks
    /// the default (all logical cores).
    #[arg(long)]
    pub threads: Option<usize>,

    /// Optional sample → batch TSV. See
    /// [`crate::pop_var_caller::batch_assignment`]. When omitted,
    /// every cohort sample is assigned the default batch
    /// `all_samples`.
    #[arg(long)]
    pub batch_assignment: Option<PathBuf>,

    /// One or more cohort `.psp` files.
    #[arg(required = true)]
    pub psp_files: Vec<PathBuf>,

    // ===== Advanced — Stopping mode ============================
    /// Convergence tolerance on `max_s |c_s_new − c_s_prev|`
    /// between consecutive block snapshots.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_STABILITY_TOLERANCE,
        value_parser = parsers::parse_stability_tolerance,
        help_heading = "Advanced — Stopping mode",
    )]
    pub stability_tolerance: f64,

    /// Consecutive within-tolerance snapshots required to declare
    /// convergence.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_STABILITY_BLOCKS,
        value_parser = parsers::parse_stability_blocks,
        help_heading = "Advanced — Stopping mode",
    )]
    pub stability_blocks: u32,

    // ===== Advanced — Online EM ================================
    /// Sites per heartbeat — controls both the parameter refresh
    /// and the stability snapshot.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_BLOCK_SIZE,
        value_parser = parsers::parse_block_size,
        help_heading = "Advanced — Online EM",
    )]
    pub block_size: u32,

    // ===== Advanced — Informative-site cuts ====================
    /// Per-sample minimum read depth at a site to consider the
    /// (sample, site) pair informative.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MIN_DEPTH,
        value_parser = parsers::parse_min_depth,
        help_heading = "Advanced — Informative-site cuts",
    )]
    pub min_depth: u32,

    /// Per-sample minimum observed major-allele fraction to call
    /// hom-major at a site.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MIN_MAJOR_FRACTION,
        value_parser = parsers::parse_min_major_fraction,
        help_heading = "Advanced — Informative-site cuts",
    )]
    pub min_major_fraction: f64,

    /// Cohort-summed minimum minor-allele read count for a site to
    /// be informative.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MIN_COHORT_MINOR_COUNT,
        help_heading = "Advanced — Informative-site cuts",
    )]
    pub min_cohort_minor_count: u32,

    /// Cohort-summed minimum minor-allele fraction for a site to
    /// be informative.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MIN_COHORT_MINOR_FRACTION,
        value_parser = parsers::parse_min_cohort_minor_fraction,
        help_heading = "Advanced — Informative-site cuts",
    )]
    pub min_cohort_minor_fraction: f64,

    /// Batches strictly smaller than this floor get `c_s = 0`
    /// regardless of the EM. Singleton batches are always floored.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_MIN_BATCH_SIZE_FOR_CONTAMINATION,
        value_parser = parsers::parse_min_batch_size,
        help_heading = "Advanced — Informative-site cuts",
    )]
    pub min_batch_size: u32,

    // ===== Advanced — Priors ===================================
    /// Dirichlet pseudocount on REF reads in the `q_b` update.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_REF_PSEUDOCOUNT,
        value_parser = parsers::parse_contam_pseudocount,
        help_heading = "Advanced — Priors",
    )]
    pub ref_pseudocount: f64,

    /// Dirichlet pseudocount on SNP-alt reads in the `q_b` update.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_SNP_ALT_PSEUDOCOUNT,
        value_parser = parsers::parse_contam_pseudocount,
        help_heading = "Advanced — Priors",
    )]
    pub snp_alt_pseudocount: f64,

    /// Dirichlet pseudocount on indel-alt reads in the `q_b` update.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_INDEL_ALT_PSEUDOCOUNT,
        value_parser = parsers::parse_contam_pseudocount,
        help_heading = "Advanced — Priors",
    )]
    pub indel_alt_pseudocount: f64,

    /// Initial per-sample `c_s` seed.
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_C_S_INIT,
        value_parser = parsers::parse_c_s_init,
        help_heading = "Advanced — Priors",
    )]
    pub c_s_init: f64,

    /// Initial per-batch `q_b` allele-class prior (one shared value
    /// per class).
    #[arg(
        long,
        hide_short_help = true,
        default_value_t = DEFAULT_Q_B_INIT_PER_CLASS,
        value_parser = parsers::parse_q_b_init_per_class,
        help_heading = "Advanced — Priors",
    )]
    pub q_b_init_per_class: f64,
}

// ---------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------

/// Errors surfaced by [`run_estimate_contamination`].
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum EstimateContaminationCliError {
    #[error("io: {0}")]
    Io(#[from] io::Error),

    #[error("psp reader: {0}")]
    PspReader(#[from] PspReadError),

    #[error("merger: {0}")]
    Merger(#[from] PerPositionMergerError),

    #[error("contamination side-pass: {0}")]
    Engine(#[from] ContaminationEstimationError),

    #[error("contamination artefact: {0}")]
    Artefact(#[from] ContaminationArtefactError),

    #[error("batch assignment: {0}")]
    Batches(#[from] BatchAssignmentError),

    /// A `.psp` file's `reference` field doesn't match the basename
    /// of `--reference`. Surfaces *before* the side-pass starts.
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

    /// Internal: `SystemTime::now()` formatting failed. Realistically
    /// can only fire if the clock is set before the epoch.
    #[error("internal: failed to format current timestamp as RFC3339")]
    TimestampFormat,
}

// ---------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------

/// Run the contamination side-pass and write the artefact to disk.
///
/// Behaviour, step-by-step:
///
/// 1. Size the global rayon pool when `--threads` is set.
/// 2. Open every `.psp` with [`PspReader`].
/// 3. Cross-check that every header's `reference` basename matches
///    `--reference`.
/// 4. Load `--batch-assignment` if supplied, else use the empty
///    mapping (every sample → `all_samples`).
/// 5. Build a [`ContaminationEstimationConfig`] from CLI args and
///    validate it.
/// 6. Combine the readers into a [`PerPositionMerger`] (k-way merge
///    over the cohort).
/// 7. Run [`estimate_contamination`].
/// 8. Convert the engine-side
///    [`ContaminationEstimates`](crate::var_calling::contamination_estimation::ContaminationEstimates)
///    into a [`ContaminationArtefact`] and write it via atomic
///    tmp+rename.
/// 9. Print a one-shot stderr run-summary block.
pub fn run_estimate_contamination(
    args: &EstimateContaminationArgs,
) -> Result<(), EstimateContaminationCliError> {
    // 1. Rayon pool.
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|_| EstimateContaminationCliError::RayonAlreadyConfigured)?;
    }

    // 2. Open readers + collect sample names.
    let mut readers: Vec<PspReader<BufReader<File>>> = Vec::with_capacity(args.psp_files.len());
    for path in &args.psp_files {
        let file = File::open(path).map_err(EstimateContaminationCliError::Io)?;
        let buf = BufReader::with_capacity(64 * 1024, file);
        let reader = PspReader::new(buf)?;
        readers.push(reader);
    }

    // 3. Cross-check references.
    let supplied_ref = basename(&args.reference);
    for (path, reader) in args.psp_files.iter().zip(readers.iter()) {
        if reader.header().reference != supplied_ref {
            return Err(EstimateContaminationCliError::ReferenceMismatch {
                psp: path.clone(),
                psp_ref: reader.header().reference.clone(),
                supplied_ref: supplied_ref.clone(),
            });
        }
    }

    // 4. Batch assignment.
    let batches = match &args.batch_assignment {
        Some(p) => BatchAssignment::from_tsv(p)?,
        None => BatchAssignment::empty(),
    };

    // 5. Build dense batch index.
    let sample_names: Vec<String> = readers.iter().map(|r| r.header().sample.clone()).collect();
    let DenseBatches {
        sample_to_batch,
        batch_id_for_idx,
    } = build_dense_batches(&sample_names, &batches);
    let n_batches = batch_id_for_idx.len();
    let n_samples = sample_names.len();

    // 6. Build + validate config.
    let cfg = ContaminationEstimationConfig {
        stopping_mode: StoppingMode::Convergence {
            tolerance: args.stability_tolerance,
            stability_blocks: args.stability_blocks,
        },
        block_size: args.block_size,
        min_depth: args.min_depth,
        min_major_fraction: args.min_major_fraction,
        min_cohort_minor_count: args.min_cohort_minor_count,
        min_cohort_minor_fraction: args.min_cohort_minor_fraction,
        min_batch_size_for_contamination: args.min_batch_size,
        ref_pseudocount: args.ref_pseudocount,
        snp_alt_pseudocount: args.snp_alt_pseudocount,
        indel_alt_pseudocount: args.indel_alt_pseudocount,
        c_s_init: args.c_s_init,
        q_b_init_per_class: args.q_b_init_per_class,
    };
    cfg.validate()?;

    // 7. Cross-check chromosomes across readers, then build merger.
    let chromosomes = check_chromosome_agreement(&readers)?;
    let record_iters: Vec<_> = readers.iter_mut().map(|r| r.records()).collect();
    let merger = PerPositionMerger::new(record_iters, sample_names.clone(), chromosomes)?;

    // 8. Run side-pass.
    let estimates = estimate_contamination(
        merger,
        n_samples,
        sample_to_batch.clone(),
        n_batches,
        cfg.clone(),
    )?;

    // 9. Build artefact + atomic write.
    let artefact =
        build_artefact_from_estimates(&estimates, &sample_names, &batch_id_for_idx, args, &cfg)?;
    artefact.write(&args.output)?;

    // 10. Stderr summary.
    print_run_summary(&estimates, &sample_names, &batch_id_for_idx);

    Ok(())
}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

/// Dense batch table built from the cohort's sample names and the
/// (possibly empty) [`BatchAssignment`]. Iteration order on the
/// `batch_id_for_idx` vector is the order each batch is first seen in
/// the cohort — matching the order `sample_to_batch[i]` indexes into.
struct DenseBatches {
    sample_to_batch: Vec<usize>,
    batch_id_for_idx: Vec<String>,
}

fn build_dense_batches(sample_names: &[String], batches: &BatchAssignment) -> DenseBatches {
    let mut batch_idx_for_id: HashMap<String, usize> = HashMap::new();
    let mut batch_id_for_idx: Vec<String> = Vec::new();
    let mut sample_to_batch: Vec<usize> = Vec::with_capacity(sample_names.len());
    for sample in sample_names {
        let batch_label = batches.batch_for(sample);
        let dense_idx = match batch_idx_for_id.get(batch_label) {
            Some(idx) => *idx,
            None => {
                let idx = batch_id_for_idx.len();
                batch_idx_for_id.insert(batch_label.to_string(), idx);
                batch_id_for_idx.push(batch_label.to_string());
                idx
            }
        };
        sample_to_batch.push(dense_idx);
    }
    DenseBatches {
        sample_to_batch,
        batch_id_for_idx,
    }
}

/// Convert the engine-side [`ContaminationEstimates`] back into a
/// [`ContaminationArtefact`] for serialisation. Engine convention:
/// floored batches carry an all-zero `q_b` vector — passed through
/// verbatim (the artefact validator accepts all-zero rows). Floored
/// samples (`c_s_per_sample[i] == None`) get `contamination_fraction = 0.0`.
fn build_artefact_from_estimates(
    estimates: &ContaminationEstimates,
    sample_names: &[String],
    batch_id_for_idx: &[String],
    args: &EstimateContaminationArgs,
    cfg: &ContaminationEstimationConfig,
) -> Result<ContaminationArtefact, EstimateContaminationCliError> {
    let provenance = Provenance {
        tool: "pop_var_caller".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        subcommand: "estimate-contamination".to_string(),
        created: rfc3339_now()
            .parse::<Datetime>()
            .map_err(|_| EstimateContaminationCliError::TimestampFormat)?,
        inputs: ProvenanceInputs {
            reference: basename(&args.reference),
            input_psps: args.psp_files.iter().map(|p| basename(p)).collect(),
            batch_assignment: args.batch_assignment.as_ref().map(|p| basename(p)),
        },
    };

    let parameters = build_parameter_map(args, cfg);

    let batches: Vec<BatchEntry> = (0..batch_id_for_idx.len())
        .map(|batch_idx| {
            // Internal access via the engine's `effective_c_s` /
            // `q_b_for_sample` accessors needs a sample-side index;
            // build `BatchEntry` directly from the engine's pub(crate)
            // shape would mean exposing more. Instead pick any sample
            // in this batch and read its q_b through the accessor.
            // Falls back to all-zero (floored signal) if the engine
            // returned an unused batch.
            let representative_sample = estimates
                .sample_to_batch
                .iter()
                .position(|b| *b == batch_idx);
            let q_b: [f64; 3] = match representative_sample {
                Some(s) => *estimates.q_b_for_sample(s),
                None => [0.0; 3],
            };
            BatchEntry {
                id: batch_id_for_idx[batch_idx].clone(),
                contaminant_ref_prob: q_b[0],
                contaminant_snp_alt_prob: q_b[1],
                contaminant_indel_alt_prob: q_b[2],
            }
        })
        .collect();

    let samples: Vec<SampleEntry> = sample_names
        .iter()
        .enumerate()
        .map(|(i, name)| SampleEntry {
            name: name.clone(),
            batch: batch_id_for_idx[estimates.sample_to_batch[i]].clone(),
            contamination_fraction: estimates.effective_c_s(i),
        })
        .collect();

    Ok(ContaminationArtefact {
        provenance,
        parameters,
        batches,
        samples,
    })
}

fn build_parameter_map(
    args: &EstimateContaminationArgs,
    cfg: &ContaminationEstimationConfig,
) -> BTreeMap<String, toml::Value> {
    let mut p = BTreeMap::new();
    // Stopping mode (v1 = convergence-only).
    p.insert(
        "stopping_mode".into(),
        toml::Value::String("convergence".into()),
    );
    p.insert(
        "stability_tolerance".into(),
        toml::Value::Float(args.stability_tolerance),
    );
    p.insert(
        "stability_blocks".into(),
        toml::Value::Integer(args.stability_blocks as i64),
    );
    p.insert(
        "block_size".into(),
        toml::Value::Integer(args.block_size as i64),
    );
    p.insert(
        "min_depth".into(),
        toml::Value::Integer(args.min_depth as i64),
    );
    p.insert(
        "min_major_fraction".into(),
        toml::Value::Float(args.min_major_fraction),
    );
    p.insert(
        "min_cohort_minor_count".into(),
        toml::Value::Integer(args.min_cohort_minor_count as i64),
    );
    p.insert(
        "min_cohort_minor_fraction".into(),
        toml::Value::Float(args.min_cohort_minor_fraction),
    );
    p.insert(
        "min_batch_size".into(),
        toml::Value::Integer(args.min_batch_size as i64),
    );
    p.insert(
        "ref_pseudocount".into(),
        toml::Value::Float(args.ref_pseudocount),
    );
    p.insert(
        "snp_alt_pseudocount".into(),
        toml::Value::Float(args.snp_alt_pseudocount),
    );
    p.insert(
        "indel_alt_pseudocount".into(),
        toml::Value::Float(args.indel_alt_pseudocount),
    );
    p.insert("c_s_init".into(), toml::Value::Float(args.c_s_init));
    p.insert(
        "q_b_init_per_class".into(),
        toml::Value::Float(args.q_b_init_per_class),
    );
    p.insert(
        "threads".into(),
        toml::Value::Integer(rayon::current_num_threads() as i64),
    );
    // `cfg` is currently a superset that mirrors `args` 1:1; recorded
    // for future-proofing in case the orchestrator ever derives a
    // value that doesn't appear directly on `args`.
    let _ = cfg;
    p
}

fn print_run_summary(
    estimates: &ContaminationEstimates,
    sample_names: &[String],
    batch_id_for_idx: &[String],
) {
    let sites_processed = match &estimates.source {
        ContaminationEstimateSource::SidePass {
            sites_processed, ..
        } => *sites_processed,
        // The contamination subcommand always runs the side-pass, so
        // any other variant would mean the engine reshaped its API
        // and the orchestrator wasn't updated.
        _ => 0,
    };
    let mut floored: usize = 0;
    let mut per_batch_size: HashMap<usize, usize> = HashMap::new();
    for (i, b_idx) in estimates.sample_to_batch.iter().enumerate() {
        if estimates.c_s_per_sample[i].is_none() {
            floored += 1;
        }
        *per_batch_size.entry(*b_idx).or_insert(0) += 1;
    }
    eprintln!(
        "estimate-contamination: n_samples={} n_batches={} sites_processed={} floored_samples={}",
        sample_names.len(),
        batch_id_for_idx.len(),
        sites_processed,
        floored,
    );
    let mut per_batch: Vec<(usize, usize)> = per_batch_size.into_iter().collect();
    per_batch.sort_by_key(|(idx, _)| *idx);
    for (idx, size) in per_batch {
        eprintln!("  batch={} size={}", batch_id_for_idx[idx], size);
    }
}

/// Format the current UTC time as a TOML-compatible RFC3339 string.
/// Mirrors the helper in `pop_var_caller::cli` to avoid coupling the
/// two subcommands' implementations.
fn rfc3339_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = (secs / 86_400) as i64;
    let sod = secs % 86_400;
    let (y, m, d) = civil_from_days(days);
    let h = sod / 3600;
    let min = (sod % 3600) / 60;
    let s = sod % 60;
    format!("{y:04}-{m:02}-{d:02}T{h:02}:{min:02}:{s:02}Z")
}

/// Howard Hinnant's `civil_from_days`. Identical to the helper in
/// `pop_var_caller::cli`; duplicated rather than re-exported to keep
/// the subcommands self-contained.
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let y = if m <= 2 { y + 1 } else { y };
    (y, m, d)
}

fn basename(p: &Path) -> String {
    p.file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| p.to_string_lossy().into_owned())
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dense_batches_no_assignment_assigns_all_samples() {
        let names: Vec<String> = ["a", "b", "c"].iter().map(|s| s.to_string()).collect();
        let DenseBatches {
            sample_to_batch,
            batch_id_for_idx,
        } = build_dense_batches(&names, &BatchAssignment::empty());
        assert_eq!(sample_to_batch, vec![0, 0, 0]);
        assert_eq!(batch_id_for_idx, vec!["all_samples".to_string()]);
    }

    #[test]
    fn dense_batches_groups_by_label_preserving_first_seen_order() {
        // Build a BatchAssignment via the parser to avoid leaking its
        // internal state into this test.
        let tsv = "sample\tbatch\nb\tlane_7\na\tlane_3\nc\tlane_3\n";
        let assignment = {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("b.tsv");
            std::fs::write(&path, tsv).unwrap();
            BatchAssignment::from_tsv(&path).unwrap()
        };
        let names: Vec<String> = ["a", "b", "c"].iter().map(|s| s.to_string()).collect();
        let DenseBatches {
            sample_to_batch,
            batch_id_for_idx,
        } = build_dense_batches(&names, &assignment);
        // First sample is `a` → batch lane_3 → dense_idx 0.
        // Second sample is `b` → batch lane_7 → dense_idx 1.
        // Third sample is `c` → batch lane_3 → dense_idx 0.
        assert_eq!(sample_to_batch, vec![0, 1, 0]);
        assert_eq!(
            batch_id_for_idx,
            vec!["lane_3".to_string(), "lane_7".to_string()]
        );
    }

    #[test]
    fn dense_batches_falls_back_for_missing_samples() {
        // `c` is not in the TSV; falls back to `all_samples`.
        let tsv = "sample\tbatch\na\tlane_3\n";
        let assignment = {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("b.tsv");
            std::fs::write(&path, tsv).unwrap();
            BatchAssignment::from_tsv(&path).unwrap()
        };
        let names: Vec<String> = ["a", "c"].iter().map(|s| s.to_string()).collect();
        let DenseBatches {
            sample_to_batch,
            batch_id_for_idx,
        } = build_dense_batches(&names, &assignment);
        assert_eq!(sample_to_batch, vec![0, 1]);
        assert_eq!(
            batch_id_for_idx,
            vec!["lane_3".to_string(), "all_samples".to_string()]
        );
    }

    #[test]
    fn rfc3339_now_parses_as_toml_datetime() {
        let s = rfc3339_now();
        let _: Datetime = s.parse().expect("rfc3339_now must parse as toml Datetime");
    }

    #[test]
    fn civil_from_days_matches_known_dates() {
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        assert_eq!(civil_from_days(364), (1970, 12, 31));
        assert_eq!(civil_from_days(10957), (2000, 1, 1));
        assert_eq!(civil_from_days(19782), (2024, 2, 29));
    }

    #[test]
    fn basename_strips_directory() {
        assert_eq!(basename(Path::new("/data/grch38.fa")), "grch38.fa");
        assert_eq!(basename(Path::new("ref.fa")), "ref.fa");
    }
}
