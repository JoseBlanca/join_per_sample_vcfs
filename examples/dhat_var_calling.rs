//! Heap-profile the cohort `.psp` → VCF pipeline (`run_var_calling`,
//! Stages 3–6) against a real on-disk fixture.
//!
//! Companion to the `perf_ours_joint.py` scaling experiment — that
//! script measures *peak RSS* via psutil polling, which tells us
//! memory grows ~135 MB per added sample (linear in N=1..18 on the
//! tomato1 fixture) but says nothing about *where* the bytes live.
//! This driver swaps in `dhat::Alloc` so every allocation is
//! attributed to a stack trace; the resulting `dhat-heap.json` ranks
//! sites by total bytes / lifetime bytes / block count and answers
//! "if I run with N=18 samples on a tomato 100-kbp-window cohort,
//! where does the RAM actually go?".
//!
//! Build and run inside the dev container:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example dhat_var_calling \
//!     --features dhat-heap -- \
//!     --psp-dir benchmarks/tomato1/results/ours/cohort/psp \
//!     --n-samples 18 \
//!     --reference $HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa \
//!     --output tmp/dhat_var_calling.vcf \
//!     --threads 4
//! ```
//!
//! Produces `dhat-heap.json` in the working directory. Open it at
//! <https://nnethercote.github.io/dh_view/dh_view.html> to see
//! allocation sites ranked by total bytes / lifetime bytes / blocks.
//!
//! Subset selection mirrors `benchmarks/tomato1/scripts/perf_common.py`
//! `pick_subset`: sorted alphabetical, first N. Cross-checking against
//! the `ours_joint.tsv` row at the same N is therefore one-to-one.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;

use pop_var_caller::pop_var_caller::cli::shared_args::CohortPipelineArgs;
use pop_var_caller::pop_var_caller::var_calling::{VarCallingArgs, run_var_calling};
use pop_var_caller::var_calling::dust_filter::{DEFAULT_DUST_THRESHOLD, DEFAULT_DUST_WINDOW};
use pop_var_caller::var_calling::per_group_merger::{
    DEFAULT_MAX_ALLELES_LH_CALC, DEFAULT_MAX_ALLELES_PER_RECORD, DEFAULT_PLOIDY,
};
use pop_var_caller::var_calling::posterior_engine::{
    DEFAULT_COMPOUND_ALT_PSEUDOCOUNT, DEFAULT_CONVERGENCE_THRESHOLD,
    DEFAULT_INBREEDING_COEFFICIENT, DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MAX_GQ_PHRED,
    DEFAULT_MAX_ITERATIONS, DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT,
};
use pop_var_caller::var_calling::variant_grouping::DEFAULT_MAX_VARIANT_GROUP_SPAN;
use pop_var_caller::var_calling::{DEFAULT_MIN_ALT_OBS_PER_SAMPLE, DEFAULT_MIN_QUAL_PHRED};

#[derive(Parser, Debug)]
#[command(about = "Heap-profile run_var_calling on a real PSP cohort fixture")]
struct Cli {
    /// Directory containing the per-sample `.psp` files. The first
    /// `--n-samples` files in sorted order are picked — same selection
    /// rule the perf_*.py scripts use, so dhat numbers line up
    /// one-to-one with the perf TSVs.
    #[arg(long)]
    psp_dir: PathBuf,

    /// Number of samples to include (subset of `--psp-dir`).
    #[arg(long)]
    n_samples: usize,

    /// Reference FASTA matching the PSPs' per-contig MD5s. `.fai`
    /// sibling required.
    #[arg(long)]
    reference: PathBuf,

    /// Output VCF.
    #[arg(long)]
    output: PathBuf,

    /// Rayon worker count. `None` → rayon default (all logical cores).
    /// Set to `4` to match `perf_ours_joint.py` defaults.
    #[arg(long)]
    threads: Option<usize>,

    /// Bypass the DUST (low-complexity) filter.
    #[arg(long, default_value_t = false)]
    no_complexity_filter: bool,

    /// Restrict to a BED of regions (keeps the profile run small).
    #[arg(long)]
    regions: Option<PathBuf>,

    /// `0` = legacy single-pull batch loader (pre-rewrite shape);
    /// `>0` = streaming `StreamingBlockLoader` with this variant target
    /// per block (production default is 1024). Set this to profile the
    /// *current* memory-bounded path.
    #[arg(long, default_value_t = 0)]
    target_variants_per_chunk: u32,

    /// Profile the low-memory two-pass producer.
    #[arg(long, default_value_t = false)]
    low_memory: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();

    let mut psps: Vec<PathBuf> = std::fs::read_dir(&args.psp_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("psp"))
        .collect();
    psps.sort();
    if args.n_samples > psps.len() {
        return Err(format!(
            "requested --n-samples {} but only {} .psp files in {}",
            args.n_samples,
            psps.len(),
            args.psp_dir.display(),
        )
        .into());
    }
    psps.truncate(args.n_samples);

    eprintln!(
        "dhat var-calling: N={} threads={:?} no_complexity_filter={}",
        args.n_samples, args.threads, args.no_complexity_filter,
    );
    eprintln!("inputs:");
    for p in &psps {
        eprintln!("  {}", p.display());
    }

    let var_calling_args = VarCallingArgs {
        reference: args.reference,
        output: args.output.clone(),
        regions: args.regions.clone(),
        threads: args.threads,
        contamination_estimates: None,
        no_complexity_filter: args.no_complexity_filter,
        target_variants_per_chunk: args.target_variants_per_chunk,
        low_memory: args.low_memory,
        psp_files: psps,
        cohort: CohortPipelineArgs {
            ploidy: DEFAULT_PLOIDY,
            complexity_window: DEFAULT_DUST_WINDOW,
            complexity_threshold: DEFAULT_DUST_THRESHOLD,
            var_group_max_span: DEFAULT_MAX_VARIANT_GROUP_SPAN,
            max_alleles_per_var: DEFAULT_MAX_ALLELES_PER_RECORD,
            max_alleles_lh_calc: DEFAULT_MAX_ALLELES_LH_CALC,
            inbreeding_coefficient: DEFAULT_INBREEDING_COEFFICIENT,
            em_convergence_threshold: DEFAULT_CONVERGENCE_THRESHOLD,
            em_max_iterations: DEFAULT_MAX_ITERATIONS,
            ref_pseudocount: DEFAULT_REF_PSEUDOCOUNT,
            snp_alt_pseudocount: DEFAULT_SNP_ALT_PSEUDOCOUNT,
            indel_alt_pseudocount: DEFAULT_INDEL_ALT_PSEUDOCOUNT,
            compound_alt_pseudocount: DEFAULT_COMPOUND_ALT_PSEUDOCOUNT,
            max_gq_phred: DEFAULT_MAX_GQ_PHRED,
            min_qual_phred: DEFAULT_MIN_QUAL_PHRED,
            min_alt_obs_per_sample: DEFAULT_MIN_ALT_OBS_PER_SAMPLE,
            no_allele_balance_filter: false,
            min_allele_balance_log_lr:
                pop_var_caller::var_calling::allele_balance::DEFAULT_AB_MIN_LOG_LR,
            allele_balance_concentration:
                pop_var_caller::var_calling::allele_balance::DEFAULT_AB_CONCENTRATION,
            emit_gp: false,
            paralog_fdr: 0.0,
            no_paralog_filter: true,
            do_not_drop_dup_artifacts: false,
        },
    };

    // ---------- start dhat measurement ----------
    // Scope dhat to just the run_var_calling call so the subset
    // discovery / arg setup above don't show up in the report.
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let run_start = Instant::now();
    run_var_calling(&var_calling_args)?;
    let elapsed = run_start.elapsed();

    // Profiler drop at end of main writes dhat-heap.json.
    eprintln!(
        "dhat var-calling: wall_time={:.2?} → {}",
        elapsed,
        args.output.display(),
    );

    Ok(())
}
