//! `ssr-call` subcommand — Stage 2 of the SSR caller. Reads the per-sample
//! `.ssr.psp` evidence files produced by `ssr-pileup`, merges them by catalog locus
//! across the cohort, and (Phases 2–3) genotypes each locus into a multi-sample VCF.
//!
//! The struct below is the authoritative knob list; [`run_ssr_call`] will translate
//! it into the cohort driver config. **Phase 1 of the build is scaffolding only** —
//! [`run_ssr_call`] is a stub until the reader/merger lands (build plan
//! `doc/devel/implementation_plans/ssr_call_reading.md`).

use std::path::PathBuf;

use clap::Args;
use thiserror::Error;

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

    /// Output multi-sample VCF path.
    #[arg(long)]
    pub output: PathBuf,

    /// Worker threads for the EM pool. The output is identical for any value
    /// (catalog-ordered, frozen-parameter genotyping).
    #[arg(long, default_value_t = 4, help_heading = "Advanced")]
    pub threads: usize,

    /// Bounded-queue depth between the merge producer and the EM workers — the peak
    /// resident-loci knob (trades throughput for RSS).
    #[arg(long, default_value_t = DEFAULT_QUEUE_DEPTH, help_heading = "Advanced")]
    pub queue_depth: usize,
}

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
/// **Phase 1 stub:** the reader/merger and EM are not built yet, so this currently
/// only validates that the CLI is wired and returns. The cohort driver replaces this
/// body in later phases.
pub fn run_ssr_call(_args: &SsrCallArgs) -> Result<(), SsrCallCliError> {
    Ok(())
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
    fn stub_run_is_ok() {
        let args = SsrCallArgs {
            psp_files: vec![PathBuf::from("a.ssr.psp")],
            catalog: PathBuf::from("x.ssr.catalog"),
            output: PathBuf::from("out.vcf"),
            threads: 4,
            queue_depth: DEFAULT_QUEUE_DEPTH,
        };
        assert!(run_ssr_call(&args).is_ok());
    }
}
