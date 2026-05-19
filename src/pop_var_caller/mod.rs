//! `pop_var_caller` binary namespace — CLI struct, subcommand
//! dispatch, and Stage 1 orchestration. The per-sample algorithmic
//! slices live under [`crate::per_sample_pileup`]; this module is the
//! glue that turns CLI arguments into a configured pipeline run.

pub mod batch_assignment;
pub mod cli;
pub mod psp_to_pileup;

pub use batch_assignment::{BatchAssignment, BatchAssignmentError, DEFAULT_BATCH_ID};
pub use cli::{Cli, PileupArgs, PileupCliError, PopVarCallerCommand, run_pileup};
pub use psp_to_pileup::{PspToPileupArgs, PspToPileupError, run_psp_to_pileup};
