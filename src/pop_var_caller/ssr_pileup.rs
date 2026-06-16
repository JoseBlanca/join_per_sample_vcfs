//! `ssr-pileup` subcommand — Stage 1 of the SSR caller. Genotypes one sample's
//! BAM/CRAM against an SSR catalog ([`crate::ssr::catalog`]) and writes a
//! per-locus `.ssr.psp` evidence file. The struct below is the authoritative
//! knob list; [`run_ssr_pileup`] translates it into an
//! [`SsrPileupConfig`] and drives
//! [`crate::ssr::pileup::driver::run`].

use std::path::PathBuf;

use clap::Args;
use thiserror::Error;

use crate::bam::alignment_input::{DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH};
use crate::bam::segment_reader::SegmentReadFilter;
use crate::ssr::pileup::driver::{self, DEFAULT_WINDOW, SsrPileupConfig};
use crate::ssr::pileup::fetch_reads::MAX_READS_PER_LOCUS;

/// Arguments for the `ssr-pileup` subcommand.
#[derive(Debug, Args, Clone)]
pub struct SsrPileupArgs {
    /// Coordinate-sorted, indexed BAM/CRAM file(s) for one sample.
    #[arg(required = true)]
    pub alignment_files: Vec<PathBuf>,

    /// Reference FASTA (with a sibling `.fai`).
    #[arg(long)]
    pub reference: PathBuf,

    /// The sorted `.ssr.catalog` to genotype against.
    #[arg(long)]
    pub catalog: PathBuf,

    /// Output `.ssr.psp` path.
    #[arg(long)]
    pub output: PathBuf,

    /// Sample name override; defaults to the name cross-validated across inputs.
    #[arg(long)]
    pub sample: Option<String>,

    /// Build a missing alignment index (`.csi`/`.crai`) in place instead of
    /// erroring on the first un-indexed input.
    #[arg(long)]
    pub build_index_if_missing: bool,

    /// Drop reads with MAPQ below this (counted into `n_filtered`).
    #[arg(long, default_value_t = DEFAULT_MIN_MAPQ, help_heading = "Read filter")]
    pub min_mapq: u8,

    /// Drop reads shorter than this many decoded bases.
    #[arg(long, default_value_t = DEFAULT_MIN_READ_LENGTH, help_heading = "Read filter")]
    pub min_read_length: u32,

    /// Keep duplicate-flagged reads (default: drop them).
    #[arg(long, help_heading = "Read filter")]
    pub keep_duplicates: bool,

    /// Keep QC-fail-flagged reads (default: drop them).
    #[arg(long, help_heading = "Read filter")]
    pub keep_qc_fail: bool,

    /// Pair-HMM candidate half-width (rungs) centred on the observed count.
    #[arg(long, default_value_t = DEFAULT_WINDOW, help_heading = "Advanced")]
    pub window: u16,

    /// Per-locus read depth cap (reservoir sampling).
    #[arg(long, default_value_t = MAX_READS_PER_LOCUS, help_heading = "Advanced")]
    pub max_reads_per_locus: usize,

    /// Worker threads for the per-locus pool. The output `.ssr.psp` is identical
    /// for any value (per-locus reservoir seeds + catalog-ordered writes).
    #[arg(long, default_value_t = 4, help_heading = "Advanced")]
    pub threads: usize,
}

/// Errors from the `ssr-pileup` subcommand.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum SsrPileupCliError {
    /// The Stage-1 run itself failed. The boxed source carries the typed
    /// `SsrPileupError` cause (boxed so this public error does not name the
    /// crate-internal type).
    #[error("ssr-pileup run failed")]
    Run {
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },
}

/// Translate the CLI args into the driver config (pure; separated from
/// [`run_ssr_pileup`] so the mapping is unit-testable without doing I/O).
fn to_config(args: &SsrPileupArgs) -> SsrPileupConfig {
    SsrPileupConfig {
        alignment_files: args.alignment_files.clone(),
        reference: args.reference.clone(),
        catalog: args.catalog.clone(),
        output: args.output.clone(),
        filter: SegmentReadFilter {
            min_mapq: Some(args.min_mapq),
            min_read_length: Some(args.min_read_length),
            drop_duplicate: !args.keep_duplicates,
            drop_qc_fail: !args.keep_qc_fail,
        },
        window: args.window,
        cap: args.max_reads_per_locus,
        build_index_if_missing: args.build_index_if_missing,
        sample: args.sample.clone(),
        threads: args.threads,
    }
}

/// Genotype one sample against an SSR catalog, writing a `.ssr.psp`.
pub fn run_ssr_pileup(args: &SsrPileupArgs) -> Result<(), SsrPileupCliError> {
    driver::run(&to_config(args)).map_err(|e| SsrPileupCliError::Run {
        source: Box::new(e),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The subcommand parses through the top-level CLI with its defaults.
    #[test]
    fn parses_through_the_top_level_cli() {
        use crate::pop_var_caller::cli::{Cli, PopVarCallerCommand};
        use clap::Parser;
        let cli = Cli::try_parse_from([
            "pop_var_caller",
            "ssr-pileup",
            "--reference",
            "ref.fa",
            "--catalog",
            "cat.ssr.catalog",
            "--output",
            "out.ssr.psp",
            "sample.bam",
        ])
        .unwrap();
        match cli.cmd {
            PopVarCallerCommand::SsrPileup(a) => {
                assert_eq!(a.alignment_files, vec![PathBuf::from("sample.bam")]);
                assert_eq!(a.min_mapq, DEFAULT_MIN_MAPQ);
                assert_eq!(a.min_read_length, DEFAULT_MIN_READ_LENGTH);
                assert_eq!(a.window, DEFAULT_WINDOW);
                assert_eq!(a.max_reads_per_locus, MAX_READS_PER_LOCUS);
                assert_eq!(a.threads, 4);
                assert!(!a.keep_duplicates && !a.keep_qc_fail);
            }
            other => panic!("expected SsrPileup, got {other:?}"),
        }
    }

    #[test]
    fn config_mapping_inverts_keep_flags_and_copies_knobs() {
        let cli = {
            use clap::Parser;
            crate::pop_var_caller::cli::Cli::try_parse_from([
                "pop_var_caller",
                "ssr-pileup",
                "--reference",
                "ref.fa",
                "--catalog",
                "cat.ssr.catalog",
                "--output",
                "out.ssr.psp",
                "--sample",
                "S1",
                "--keep-duplicates",
                "--min-mapq",
                "30",
                "--max-reads-per-locus",
                "42",
                "a.bam",
                "b.cram",
            ])
            .unwrap()
        };
        let args = match cli.cmd {
            crate::pop_var_caller::cli::PopVarCallerCommand::SsrPileup(a) => a,
            other => panic!("expected SsrPileup, got {other:?}"),
        };
        let cfg = to_config(&args);

        assert_eq!(cfg.alignment_files.len(), 2);
        assert_eq!(cfg.sample.as_deref(), Some("S1"));
        assert_eq!(cfg.cap, 42);
        assert_eq!(cfg.filter.min_mapq, Some(30));
        // --keep-duplicates flips drop_duplicate off; qc-fail still dropped.
        assert!(!cfg.filter.drop_duplicate);
        assert!(cfg.filter.drop_qc_fail);
    }
}
