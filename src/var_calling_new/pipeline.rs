//! Driver / wiring — section orchestration (appendix §G).
//!
//! *(today: `var_calling::driver` `drive_blocks_parallel`)*
//!
//! Spawns the three sections (producer → W callers → writer), owns the two
//! bounded `crossbeam-channel` hand-offs, joins, and propagates the first
//! error. This is the new entry point the CLI repoints to at the swap (P7).
//!
//! - producer → caller queue: MPMC, count-capped, carries `PileupCohortChunk`;
//! - caller → writer queue: single-consumer, carries `CalledChunk`.
//!
//! Phase 4 builds the real wiring; Phase 0 ships only the entry-point stub the
//! byte-identity oracle harness drives (red until P4).

use thiserror::Error;

use crate::pop_var_caller::var_calling::VarCallingArgs;

/// Errors surfaced by the re-architected pipeline entry point.
///
/// Phase 0 carries only the `NotImplemented` arm so the oracle harness has a
/// concrete `Err` to observe; the real error taxonomy lands with the wiring in
/// Phase 4 (cf. `var_calling::ChunkDriverError`).
#[derive(Debug, Error)]
pub enum PipelineError {
    /// The new pipeline is not yet wired end-to-end (Phases 1–4 pending).
    #[error(
        "var_calling_new pipeline not yet implemented (re-architecture phases 1-4 pending); \
         this stub exists so the byte-identity oracle harness compiles and runs red by design"
    )]
    NotImplemented,
}

/// Re-architected cohort `.psp` → VCF entry point.
///
/// Mirrors the old CLI-level
/// [`run_var_calling`](crate::pop_var_caller::var_calling::run_var_calling) so
/// the byte-identity oracle can drive both pipelines from the same
/// [`VarCallingArgs`] and diff the resulting VCFs. The argument-type coupling
/// to the CLI layer is temporary: at the P7 swap this becomes the production
/// entry the CLI invokes directly.
///
/// **Phase 0: stub.** Always returns [`PipelineError::NotImplemented`]; the
/// producer/caller/writer wiring lands in Phase 4, where the oracle turns
/// green.
pub fn run_var_calling(_args: &VarCallingArgs) -> Result<(), PipelineError> {
    Err(PipelineError::NotImplemented)
}
