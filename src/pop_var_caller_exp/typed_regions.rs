//! The `type-regions` subcommand: run step 3's walk over a reference and
//! write the typed-region partition to a file. This module owns the `Args`
//! struct, the `run_typed_regions` driver, and its `#[non_exhaustive]`
//! error enum.
//!
//! **Scaffold (Milestone A).** Today this is the binary skeleton only:
//! `TypedRegionsArgs` carries the required flags so `type-regions --help`
//! runs, and `run_typed_regions` is a stub. Milestone B fleshes out the
//! full knob list and the real error enum; Milestone E wires the driver.
//! See `doc/devel/ng/impl_plan/typed_regions_cli.md`.

use std::path::PathBuf;

use clap::Args;
use thiserror::Error;

/// `type-regions` arguments. Scaffold shape — Milestone B adds the walk's
/// knobs (`--min-period`, `--max-str-len`, `--min-copies`, …) under an
/// "Advanced" help heading; `run_typed_regions` will translate the whole
/// struct into a `TypedRegionConfig`.
#[derive(Debug, Args, Clone)]
pub struct TypedRegionsArgs {
    /// Reference FASTA. Its `.fai` is read, or created if absent (spec T1).
    #[arg(long)]
    pub reference: PathBuf,

    /// Output partition file (a plain-text BED-shaped TSV).
    #[arg(long)]
    pub output: PathBuf,

    /// BED of regions to emit; omit for the whole genome. Choosing a BED
    /// costs scan time, not memory or correctness (spec T10).
    #[arg(long)]
    pub regions: Option<PathBuf>,
}

/// Errors from the `type-regions` subcommand. Scaffold shape — Milestone B
/// replaces this with the full set (`Reference`, `ContigTooLong`, `Bed`,
/// `Walk`, `Output`).
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum TypedRegionsCliError {
    /// The driver is not wired yet (Milestone E).
    #[error("type-regions is not implemented yet")]
    NotImplemented,
}

/// Run the walk and write the partition. Stub until Milestone E.
pub fn run_typed_regions(_args: &TypedRegionsArgs) -> Result<(), TypedRegionsCliError> {
    Err(TypedRegionsCliError::NotImplemented)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pop_var_caller_exp::cli::{Cli, PopVarCallerExpCommand};
    use clap::Parser;

    /// The subcommand parses through the top-level CLI with the required
    /// flags. Mirrors `ssr-catalog`'s `parses_through_the_top_level_cli`.
    #[test]
    fn parses_through_the_top_level_cli() {
        let cli = Cli::try_parse_from([
            "pop_var_caller_exp",
            "type-regions",
            "--reference",
            "r.fa",
            "--output",
            "regions.tsv",
        ])
        .unwrap();
        let PopVarCallerExpCommand::TypeRegions(args) = cli.cmd;
        assert_eq!(args.reference, PathBuf::from("r.fa"));
        assert_eq!(args.output, PathBuf::from("regions.tsv"));
        assert_eq!(args.regions, None);
    }
}
