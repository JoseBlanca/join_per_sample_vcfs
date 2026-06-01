//! Base Alignment Quality — per-read HMM that caps effective base
//! quality at `min(BQ, BAQ)` where `BAQ[i]` is the Phred-scaled
//! posterior probability that base `i` is correctly aligned to its
//! reference position.
//!
//! This module is the algorithmic core: a Rust port of htslib's
//! `probaln_glocal` HMM plus the supporting scratch buffers and
//! configuration types. It has no dependency on the per-sample-pileup
//! walker — the pileup-specific per-read driver and the rayon stream
//! adapter live in
//! [`crate::pileup::per_sample::baq_engine`] and
//! [`crate::pileup::per_sample::baq_stream`].
//!
//! The algorithm and parameter choices are pinned in
//! `doc/devel/implementation_plans/baq.md`; the architecture
//! motivation is in
//! `doc/devel/specs/calling_pipeline_architecture.md` §"Per-read
//! likelihood quality".
//!
//! - [`probaln_glocal`] — the pure HMM kernel.
//! - [`ProbalnScratch`] / [`Q2P`] — preallocated work buffers and the
//!   Q-to-probability lookup table.
//! - [`BaqConfig`] / [`ProbalnError`] — configuration and structured
//!   error variants.

pub mod errors;
mod probaln;
mod scratch;

#[cfg(test)]
mod tests;

pub use errors::ProbalnError;
pub use probaln::probaln_glocal;
pub use scratch::{ProbalnScratch, Q2P};

// htslib/realn.c:115 — samtools' production gap-open probability for
// Illumina short reads.
pub const SAMTOOLS_ILLUMINA_GAP_OPEN_PROB: f32 = 1e-3;
// htslib/realn.c:115 — samtools' production gap-extension probability.
pub const SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB: f32 = 0.1;
// htslib/realn.c:115 — samtools' production band half-width.
// Widened per-read at the driver layer when a single indel forces it
// (`htslib/realn.c:195-197`).
pub const SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH: i32 = 7;

/// Parameters for `probaln_glocal`. Mirrors htslib's `probaln_par_t`
/// ([htslib/htslib/hts.h:1435](../../htslib/htslib/hts.h#L1435)).
///
/// Default values come from samtools' production defaults in
/// [`htslib/realn.c:115`](../../htslib/realn.c#L115) — see the
/// "Parameters: htslib vs GATK" decision in the plan. Prefer
/// [`BaqConfig::samtools_illumina`] at call sites so the parameter
/// source is visible.
#[derive(Debug, Clone, Copy)]
pub struct BaqConfig {
    /// Gap-open probability per base. samtools default: `1e-3`.
    /// htslib name: `d`.
    pub gap_open_prob: f32,
    /// Gap-extension probability. samtools default: `0.1`.
    /// htslib name: `e`.
    pub gap_extend_prob: f32,
    /// Band half-width for the DP. samtools default: `7`, expanded
    /// per-read at the driver layer if a single indel forces it
    /// (`htslib/realn.c:195-197`). htslib name: `bw`.
    ///
    /// **Callers that invoke `probaln_glocal` directly must perform
    /// the same CIGAR-driven expansion** — the HMM only widens `bw`
    /// up to `(l_ref - l_query).abs()`, not the per-indel expansion
    /// the engine applies.
    pub band_half_width: i32,
}

impl BaqConfig {
    /// Build the samtools Illumina-short-read defaults.
    pub fn samtools_illumina() -> Self {
        Self {
            gap_open_prob: SAMTOOLS_ILLUMINA_GAP_OPEN_PROB,
            gap_extend_prob: SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB,
            band_half_width: SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH,
        }
    }
}

impl Default for BaqConfig {
    fn default() -> Self {
        Self::samtools_illumina()
    }
}
