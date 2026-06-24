//! `ssr-call` subcommand — Stage 2 of the SSR caller. Reads the per-sample
//! `.ssr.psp` evidence files produced by `ssr-pileup`, merges them by catalog locus
//! across the cohort, genotypes each locus with the cohort EM, and writes a
//! multi-sample **VCF**.
//!
//! The struct below is the authoritative knob list; [`run_ssr_call`] translates it
//! into the cohort driver config. The driver is the two-pass streaming pipeline (arch
//! `doc/devel/architecture/ssr_call_driver.md`): a bounded burn-in freezes the
//! cross-locus-pooled parameters, then a chunk-parallel sweep genotypes every locus on
//! `--threads` workers and emits the VCF. The VCF is byte-identical for any `--threads`
//! / `--queue-depth` (each locus is a pure function of its reads + the frozen params).

use std::path::PathBuf;

use clap::Args;
use thiserror::Error;

use crate::ssr::cohort::driver::SsrCallConfig;

/// Arguments for the `ssr-call` subcommand.
#[derive(Debug, Args, Clone)]
pub struct SsrCallArgs {
    /// Per-sample `.ssr.psp` evidence files (one per sample), all genotyped against
    /// the same `--catalog`.
    #[arg(required = true)]
    pub psp_files: Vec<PathBuf>,

    /// The `.ssr.catalog` every input was genotyped against (the authoritative,
    /// ordered master locus list). Its md5 must match every input's header.
    #[arg(long)]
    pub catalog: PathBuf,

    /// Output multi-sample VCF path (one record per emitted SSR locus, `GT:GQ:REPCN`).
    #[arg(long)]
    pub output: PathBuf,

    /// Worker threads for the burn-in pool and the chunk-parallel genotyping sweep. The
    /// output VCF is byte-identical for any value; `0` is treated as single-threaded.
    #[arg(long, default_value_t = DEFAULT_THREADS, help_heading = "Advanced")]
    pub threads: usize,

    /// Loci per genotyping-sweep chunk (the resident-loci + parallelism knob). `0` (the
    /// default) lets the driver pick a chunk size; the VCF is identical for any value.
    #[arg(long, default_value_t = DEFAULT_QUEUE_DEPTH, help_heading = "Advanced")]
    pub queue_depth: usize,
}

/// Default genotyping-sweep worker thread count.
pub const DEFAULT_THREADS: usize = 4;

/// Default chunk size: `0` is the "unset" sentinel — the driver substitutes its own
/// `DEFAULT_SWEEP_CHUNK`. (A nonzero CLI default would force tiny chunks and serialize
/// the sweep.)
pub const DEFAULT_QUEUE_DEPTH: usize = 0;

/// Errors from the `ssr-call` subcommand.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum SsrCallCliError {
    /// The Stage-2 run itself failed. The boxed source carries the typed cohort
    /// driver error (boxed so this public error does not name the crate-internal
    /// type).
    #[error("ssr-call run failed")]
    Run {
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },
}

/// Run the `ssr-call` subcommand: open the cohort, merge per-sample evidence by catalog
/// locus, freeze the pooled parameters in a bounded burn-in, then stream-genotype every
/// locus on `--threads` workers and write the multi-sample VCF to `--output`.
pub fn run_ssr_call(args: &SsrCallArgs) -> Result<(), SsrCallCliError> {
    let config = SsrCallConfig {
        catalog: args.catalog.clone(),
        psp_files: args.psp_files.clone(),
        output: args.output.clone(),
        threads: args.threads,
        queue_depth: args.queue_depth,
    };
    crate::ssr::cohort::driver::run(&config).map_err(|source| SsrCallCliError::Run {
        source: Box::new(source),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    /// A throwaway top-level parser so we can exercise `SsrCallArgs` parsing in
    /// isolation, the way the real `Cli` enum drives it.
    #[derive(Debug, Parser)]
    struct TestCli {
        #[command(flatten)]
        args: SsrCallArgs,
    }

    #[test]
    fn parses_required_and_default_args() {
        let cli = TestCli::try_parse_from([
            "ssr-call",
            "a.ssr.psp",
            "b.ssr.psp",
            "--catalog",
            "x.ssr.catalog",
            "--output",
            "out.vcf",
        ])
        .expect("valid args parse");

        assert_eq!(cli.args.psp_files.len(), 2);
        assert_eq!(cli.args.catalog, PathBuf::from("x.ssr.catalog"));
        assert_eq!(cli.args.output, PathBuf::from("out.vcf"));
        assert_eq!(cli.args.threads, 4);
        assert_eq!(cli.args.queue_depth, DEFAULT_QUEUE_DEPTH);
    }

    #[test]
    fn requires_at_least_one_psp_file() {
        let parsed = TestCli::try_parse_from([
            "ssr-call",
            "--catalog",
            "x.ssr.catalog",
            "--output",
            "out.vcf",
        ]);
        assert!(parsed.is_err(), "at least one .ssr.psp is required");
    }

    #[test]
    fn run_errors_on_missing_inputs() {
        let args = SsrCallArgs {
            psp_files: vec![PathBuf::from("/no/such/a.ssr.psp")],
            catalog: PathBuf::from("/no/such/x.ssr.catalog"),
            output: PathBuf::from("/tmp/unused.vcf"),
            threads: 4,
            queue_depth: DEFAULT_QUEUE_DEPTH,
        };
        assert!(matches!(
            run_ssr_call(&args),
            Err(SsrCallCliError::Run { .. })
        ));
    }
}
