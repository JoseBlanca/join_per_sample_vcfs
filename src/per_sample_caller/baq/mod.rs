//! Base Alignment Quality — per-read HMM that caps effective base
//! quality at `min(BQ, BAQ)` where `BAQ[i]` is the Phred-scaled
//! posterior probability that base `i` is correctly aligned to its
//! reference position.
//!
//! This module is the Rust port of htslib's `probaln_glocal` plus a
//! per-read driver layered on top. The algorithm and parameter choices
//! are pinned in
//! [`ia/feature_implementation_plans/baq.md`](../../../ia/feature_implementation_plans/baq.md);
//! the architecture motivation is in
//! [`ia/specs/calling_pipeline_architecture.md`](../../../ia/specs/calling_pipeline_architecture.md)
//! §"Per-read likelihood quality".
//!
//! Layer 1 (this file) — `BaqConfig` and the pure `probaln_glocal`.
//! Layer 2 (`BaqEngine`, per-read driver) and the pipeline integration
//! land in later commits per the plan.

mod engine;
pub mod errors;
mod probaln;
mod scratch;
mod stream;

#[cfg(test)]
mod tests;

pub use engine::{BaqEngine, BaqOutcome, BaqSkipReason};
pub use errors::BaqOverflow;
pub use stream::{BaqSkipCounts, BaqStream, DEFAULT_BAQ_CHUNK_SIZE};

/// Parameters for `probaln_glocal`. Mirrors htslib's `probaln_par_t`
/// ([htslib/htslib/hts.h:1435](../../../htslib/htslib/hts.h#L1435)).
///
/// Default values come from samtools' production defaults in
/// [`htslib/realn.c:115`](../../../htslib/realn.c#L115) — see the
/// "Parameters: htslib vs GATK" decision in the plan.
#[derive(Debug, Clone, Copy)]
pub struct BaqConfig {
    /// Gap-open probability per base. samtools default: `1e-3`.
    pub d: f32,
    /// Gap-extension probability. samtools default: `0.1`.
    pub e: f32,
    /// Band half-width for the DP. samtools default: `7`, expanded
    /// per-read at the driver layer if a single indel forces it
    /// (`htslib/realn.c:195-197`).
    pub bw: i32,
}

impl Default for BaqConfig {
    fn default() -> Self {
        // htslib/realn.c:115 — Illumina short-read defaults.
        Self {
            d: 1e-3,
            e: 0.1,
            bw: 7,
        }
    }
}
