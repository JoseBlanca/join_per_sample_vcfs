//! `pop_var_caller` binary namespace — CLI struct, subcommand
//! dispatch, and Stage 1 orchestration. The per-sample algorithmic
//! slices live under [`crate::per_sample_pileup`]; this module is the
//! glue that turns CLI arguments into a configured pipeline run.

pub mod batch_assignment;
pub mod cli;
pub mod contamination_artefact;
pub mod estimate_contamination;
pub mod psp_to_pileup;
pub mod stage1_pipeline;
pub mod var_calling;
pub mod var_calling_from_bam;

pub use batch_assignment::{BatchAssignment, BatchAssignmentError, DEFAULT_BATCH_ID};
pub use cli::{Cli, PileupArgs, PileupCliError, PopVarCallerCommand, run_pileup};
pub use contamination_artefact::{ContaminationArtefact, ContaminationArtefactError};
pub use estimate_contamination::{
    EstimateContaminationArgs, EstimateContaminationCliError, run_estimate_contamination,
};
pub use psp_to_pileup::{PspToPileupArgs, PspToPileupError, run_psp_to_pileup};
pub use var_calling::{VarCallingArgs, VarCallingCliError, run_var_calling};
pub use var_calling_from_bam::{
    VarCallingFromBamArgs, VarCallingFromBamCliError, run_var_calling_from_bam,
};
