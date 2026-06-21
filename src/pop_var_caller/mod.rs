//! `pop_var_caller` binary namespace — CLI struct, subcommand
//! dispatch, and Stage 1 orchestration. The per-sample algorithmic
//! slices live under [`crate::pileup::per_sample`]; this module is the
//! glue that turns CLI arguments into a configured pipeline run.

pub mod batch_assignment;
pub mod cli;
#[doc(hidden)]
pub(crate) mod common;
pub mod contamination_artefact;
pub mod estimate_contamination;
pub mod psp_to_pileup;
pub mod ssr_call;
pub mod ssr_catalog;
pub mod ssr_pileup;
// Stage-1 pipeline seam — consumed only by `cli::run_pileup` (same crate),
// so its surface is crate-internal (review Mi10).
pub(crate) mod stage1_pipeline;
pub mod var_calling;

pub use batch_assignment::{BatchAssignment, BatchAssignmentError, DEFAULT_BATCH_ID};
pub use cli::{Cli, PileupArgs, PileupCliError, PopVarCallerCommand, run_pileup};
pub use contamination_artefact::{ContaminationArtefact, ContaminationArtefactError};
pub use estimate_contamination::{
    EstimateContaminationArgs, EstimateContaminationCliError, run_estimate_contamination,
};
pub use psp_to_pileup::{PspToPileupArgs, PspToPileupError, run_psp_to_pileup};
pub use ssr_call::{SsrCallArgs, SsrCallCliError, run_ssr_call};
pub use ssr_catalog::{SsrCatalogArgs, SsrCatalogCliError, run_ssr_catalog};
pub use ssr_pileup::{SsrPileupArgs, SsrPileupCliError, run_ssr_pileup};
pub use var_calling::{VarCallingArgs, VarCallingCliError, run_var_calling};
