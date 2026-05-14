//! `pop_var_caller` binary namespace — CLI struct, subcommand
//! dispatch, and Stage 1 orchestration. The per-sample algorithmic
//! slices live under [`crate::per_sample_caller`]; this module is the
//! glue that turns CLI arguments into a configured pipeline run.

pub mod cli;

pub use cli::{Cli, PileupArgs, PileupCliError, PopVarCallerCommand, run_pileup};
