//! `ssr-call` subcommand — Stage 2 of the SSR caller. Reads the per-sample
//! `.ssr.psp` evidence files produced by `ssr-pileup`, merges them by catalog locus
//! across the cohort, and writes the result.
//!
//! The struct below is the authoritative knob list; [`run_ssr_call`] translates it
//! into the cohort driver config. **The genotyping EM + VCF are not built yet** — the
//! driver currently writes a catalog-ordered **TSV dump** of the merged evidence
//! (a reading-layer inspection tool + VCF placeholder), and `--threads` /
//! `--queue-depth` are reserved (the run is single-threaded). Build plan:
//! `doc/devel/implementation_plans/ssr_call_reading.md`.

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

    /// Output path. Currently a catalog-ordered TSV dump of the merged evidence;
    /// the multi-sample VCF lands with the genotyping EM.
    #[arg(long)]
    pub output: PathBuf,

    /// RESERVED — EM worker threads. The genotyping EM is not built yet, so the run is
    /// single-threaded and this value currently has no effect.
    #[arg(long, default_value_t = DEFAULT_THREADS, help_heading = "Advanced")]
    pub threads: usize,

    /// RESERVED — bounded producer→worker queue depth. Not yet wired (single-threaded);
    /// currently has no effect.
    #[arg(long, default_value_t = DEFAULT_QUEUE_DEPTH, help_heading = "Advanced")]
    pub queue_depth: usize,
}

/// Default EM worker thread count (reserved until the EM lands).
pub const DEFAULT_THREADS: usize = 4;

/// Default bounded-queue depth (Q-R3; a starting point to be measured against a real
/// run, arch `ssr_call_reading.md` §5).
pub const DEFAULT_QUEUE_DEPTH: usize = 4;

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

/// Run the `ssr-call` subcommand.
///
/// Phase 3: opens the cohort, merges per-sample evidence by catalog locus, and writes
/// a catalog-ordered TSV dump of the reading layer to `--output`. The genotyping EM +
/// VCF land in a later phase; `--threads` / `--queue-depth` are reserved.
pub fn run_ssr_call(args: &SsrCallArgs) -> Result<(), SsrCallCliError> {
    // The reserved flags currently have no effect; say so loudly rather than letting a
    // user believe `--threads 32` did something.
    if args.threads != DEFAULT_THREADS {
        eprintln!(
            "warning: --threads {} has no effect yet (ssr-call is single-threaded \
             until the genotyping EM lands)",
            args.threads
        );
    }
    if args.queue_depth != DEFAULT_QUEUE_DEPTH {
        eprintln!(
            "warning: --queue-depth {} has no effect yet (the producer→worker queue \
             is not wired until the genotyping EM lands)",
            args.queue_depth
        );
    }

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
            output: PathBuf::from("/tmp/unused.tsv"),
            threads: 4,
            queue_depth: DEFAULT_QUEUE_DEPTH,
        };
        assert!(matches!(
            run_ssr_call(&args),
            Err(SsrCallCliError::Run { .. })
        ));
    }
}
