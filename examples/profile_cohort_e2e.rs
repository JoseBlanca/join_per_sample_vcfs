//! Standalone driver for profiling the cohort `.psp` → VCF pipeline
//! end-to-end (Stages 3–6, `var-calling`).
//!
//! Replicates a single real `.psp` file N times under a scratch
//! directory with rewritten `header.sample` so the cohort reader
//! sees N distinct samples, then calls [`run_var_calling`] and
//! prints elapsed wall time. Intended for one-off profiling with
//! `samply` / `perf` and for the thread-scaling baseline that picks
//! between the next code-level perf levers (rayon-over-records vs.
//! `DEFAULT_BATCH_SIZE` tuning).
//!
//! **Not a committed bench.** The `cohort_e2e_perf` criterion bench
//! is the maintained perf surface and uses fully synthetic input.
//! This example is for diagnostics against real per-sample pileup
//! data the user happens to have on disk; the data itself is never
//! committed.
//!
//! Build and run on the host (rootless podman blocks
//! `perf_event_open`, so samply / perf must be invoked outside the
//! container):
//!
//! ```text
//! ./scripts/dev.sh cargo build --release --example profile_cohort_e2e
//!
//! # Baseline wall time, all cores:
//! ./target-container/release/examples/profile_cohort_e2e \
//!     --psp tmp/SRR7279725_small.psp \
//!     --n-samples 10 \
//!     --reference /home/jose/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa \
//!     --output tmp/cohort.vcf
//!
//! # Thread-scaling sweep:
//! for T in 1 2 4 16; do
//!     ./target-container/release/examples/profile_cohort_e2e \
//!         --psp tmp/SRR7279725_small.psp --n-samples 10 \
//!         --reference /home/jose/.../S_lycopersicum_chromosomes.4.00.fa \
//!         --output tmp/cohort_T${T}.vcf --threads $T;
//! done
//!
//! # Sampling profile (record outside the container; rootless podman
//! # blocks perf_event_open inside):
//! samply record -- ./target-container/release/examples/profile_cohort_e2e \
//!     --psp tmp/SRR7279725_small.psp --n-samples 10 \
//!     --reference /home/jose/.../S_lycopersicum_chromosomes.4.00.fa \
//!     --output tmp/cohort.vcf
//! ```

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;

use pop_var_caller::pileup_record::PileupRecord;
use pop_var_caller::pop_var_caller::cli::shared_args::CohortPipelineArgs;
use pop_var_caller::pop_var_caller::var_calling::{VarCallingArgs, run_var_calling};
use pop_var_caller::psp::PspReader;
use pop_var_caller::psp::header::{
    ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance,
};
use pop_var_caller::psp::writer::PspWriter;
use pop_var_caller::var_calling::dust_filter::{DEFAULT_DUST_THRESHOLD, DEFAULT_DUST_WINDOW};
use pop_var_caller::var_calling::from_bam::pipeline::{
    DEFAULT_MIN_ALT_OBS_PER_SAMPLE, DEFAULT_MIN_MAPQ_DIFF_T, DEFAULT_MIN_QUAL_PHRED,
};
use pop_var_caller::var_calling::per_group_merger::{
    DEFAULT_MAX_ALLELES_LH_CALC, DEFAULT_MAX_ALLELES_PER_RECORD, DEFAULT_PLOIDY,
};
use pop_var_caller::var_calling::posterior_engine::{
    DEFAULT_COMPOUND_ALT_PSEUDOCOUNT, DEFAULT_CONVERGENCE_THRESHOLD,
    DEFAULT_INBREEDING_COEFFICIENT, DEFAULT_INDEL_ALT_PSEUDOCOUNT, DEFAULT_MAX_GQ_PHRED,
    DEFAULT_MAX_ITERATIONS, DEFAULT_REF_PSEUDOCOUNT, DEFAULT_SNP_ALT_PSEUDOCOUNT,
};
use pop_var_caller::var_calling::variant_grouping::DEFAULT_MAX_VARIANT_GROUP_SPAN;

#[derive(Parser, Debug)]
#[command(about = "End-to-end profiling driver for cohort var-calling")]
struct Cli {
    /// Single real `.psp` file to replicate into a synthetic cohort.
    #[arg(long)]
    psp: PathBuf,

    /// Number of samples in the synthesised cohort (= number of
    /// `.psp` copies written to `--cohort-dir`).
    #[arg(long)]
    n_samples: usize,

    /// Reference FASTA matching the input `.psp`'s per-contig MD5s.
    /// Required so [`run_var_calling`]'s `M5` check passes; the
    /// example doesn't carry a fake-fetcher fallback.
    #[arg(long)]
    reference: PathBuf,

    /// Output cohort VCF. Atomic tmp+rename — overwrites if present.
    #[arg(long)]
    output: PathBuf,

    /// Worker threads for the rayon pool. `None` → rayon default
    /// (all logical cores). Note: rayon's global pool is
    /// once-per-process — re-running this binary multiple times in
    /// the same shell, each picks its own pool.
    #[arg(long)]
    threads: Option<usize>,

    /// Bypass the DUST (low-complexity) filter.
    #[arg(long, default_value_t = false)]
    no_complexity_filter: bool,

    /// Directory under which replicated `.psp` files are written.
    /// Defaults to `tmp/cohort_synth/`. Created if missing; existing
    /// `S{:04}.psp` files at the requested N are reused (skip the
    /// rewrite) so back-to-back runs at the same N are fast.
    #[arg(long, default_value = "tmp/cohort_synth")]
    cohort_dir: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();

    // 1. Read the input .psp once into RAM. Records carry no
    //    sample-name info — `header.sample` is the only per-sample
    //    field — so we replicate the records as-is and rewrite only
    //    the header on each copy.
    let read_start = Instant::now();
    let (records, base_header) = read_psp(&args.psp)?;
    eprintln!(
        "input  : {} records, sample=\"{}\", reference=\"{}\" ({} chrom), read in {:.2?}",
        records.len(),
        base_header.sample,
        base_header.reference,
        base_header.chromosomes.len(),
        read_start.elapsed(),
    );

    // 2. Write N replicas under `--cohort-dir` (skip files that
    //    already exist at the right path — lets back-to-back runs at
    //    the same N skip the rewrite).
    std::fs::create_dir_all(&args.cohort_dir)?;
    let write_start = Instant::now();
    let psp_files = write_replicas(&args.cohort_dir, args.n_samples, &records, &base_header)?;
    eprintln!(
        "cohort : wrote {} replicas to {} in {:.2?}",
        psp_files.len(),
        args.cohort_dir.display(),
        write_start.elapsed(),
    );

    // 3. Build the cohort-pipeline knobs from defaults (matches the
    //    CLI's shared_args::CohortPipelineArgs defaults).
    let var_calling_args = VarCallingArgs {
        reference: args.reference,
        output: args.output.clone(),
        threads: args.threads,
        contamination_estimates: None,
        no_complexity_filter: args.no_complexity_filter,
        // Legacy single-pull / single-window behaviour so the
        // profile reflects the pre-rewrite per-chunk shape.
        target_variants_per_chunk: 0,
        psp_files,
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
            no_mapq_diff_filter: false,
            min_mapq_diff_t: DEFAULT_MIN_MAPQ_DIFF_T,
            emit_gp: false,
        },
    };

    // 4. Time the whole run_var_calling call. This is what samply /
    //    perf sees as the dominant self-time region.
    let run_start = Instant::now();
    run_var_calling(&var_calling_args)?;
    let elapsed = run_start.elapsed();

    eprintln!(
        "run    : n_samples={} threads={:?} no_complexity_filter={} wall_time={:.2?} → {}",
        args.n_samples,
        args.threads,
        args.no_complexity_filter,
        elapsed,
        args.output.display(),
    );

    Ok(())
}

/// Read every record from `path` into a `Vec` plus a clone of the
/// header re-shaped into `WriterHeader` (so the per-sample rewrite
/// can re-emit it via `PspWriter` directly).
fn read_psp(
    path: &std::path::Path,
) -> Result<(Vec<PileupRecord>, WriterHeader), Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let mut reader = PspReader::new(BufReader::with_capacity(64 * 1024, file))?;
    let parsed = reader.header().clone();

    let writer_header = WriterHeader {
        format_version: parsed.format_version,
        sample: parsed.sample.clone(),
        reference: parsed.reference.clone(),
        created: parsed.created,
        chromosomes: parsed
            .chromosomes
            .iter()
            .map(|c| ChromosomeEntry {
                name: c.name.clone(),
                length: c.length,
                md5: c.md5.clone(),
            })
            .collect(),
        writer: WriterProvenance {
            tool: parsed.writer.tool.clone(),
            version: parsed.writer.version.clone(),
            subcommand: parsed.writer.subcommand.clone(),
            input_crams: parsed.writer.input_crams.clone(),
            input_fasta: parsed.writer.input_fasta.clone(),
            parameters: parameters_clone(&parsed.writer.parameters),
        },
    };

    let mut records = Vec::with_capacity(1024 * 1024);
    for item in reader.records() {
        records.push(item?);
    }
    Ok((records, writer_header))
}

fn parameters_clone(src: &BTreeMap<String, ParameterValue>) -> BTreeMap<String, ParameterValue> {
    src.iter().map(|(k, v)| (k.clone(), v.clone())).collect()
}

/// Write `n_samples` copies of `records` under `dir`, each with
/// `header.sample = "S{:04}"`. Skips files that already exist at the
/// target path — lets the second run at the same N reuse the cohort
/// dir.
fn write_replicas(
    dir: &std::path::Path,
    n_samples: usize,
    records: &[PileupRecord],
    base_header: &WriterHeader,
) -> Result<Vec<PathBuf>, Box<dyn std::error::Error>> {
    let mut paths = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let sample_name = format!("S{i:04}");
        let path = dir.join(format!("{sample_name}.psp"));
        if !path.exists() {
            let mut header = base_header.clone();
            header.sample = sample_name;
            write_psp(&path, header, records)?;
        }
        paths.push(path);
    }
    Ok(paths)
}

fn write_psp(
    path: &std::path::Path,
    header: WriterHeader,
    records: &[PileupRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(path)?;
    let sink = BufWriter::with_capacity(64 * 1024, file);
    let mut writer = PspWriter::new(sink, header)?;
    for record in records {
        writer.write_record(record)?;
    }
    let buf = writer.finish()?;
    buf.into_inner()?.sync_all()?;
    Ok(())
}
