//! Errors produced by the BAQ stage.
//!
//! See [`doc/devel/implementation_plans/baq.md`](../../doc/devel/implementation_plans/baq.md)
//! §"Algorithm port" for how each variant maps onto htslib's failure modes.

use thiserror::Error;

/// Reasons `probaln_glocal` returns early without producing a usable
/// alignment. Mirrors the `INT_MIN` / `errno` exits in
/// [htslib/probaln.c:85-107](../../htslib/probaln.c#L85-L107).
///
/// Marked `#[non_exhaustive]` so future variants can be added without a
/// semver break for out-of-crate consumers.
#[non_exhaustive]
#[derive(Error, Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProbalnError {
    /// One of the caller-provided slices (`iqual`, `state`, `q`) had a
    /// length disagreeing with `query.len()`. This is an engine
    /// programming bug — the per-read driver sizes these
    /// consistently. Surfaced as an error rather than a panic so a
    /// future caller does not silently violate the contract.
    #[error("BAQ slice lengths disagree (iqual/state/q must match query)")]
    SliceLengthMismatch,

    /// The query or reference sequence does not fit in `i32` (htslib's
    /// internal index type). Practically only triggered by
    /// pathological inputs; signalled distinctly from
    /// `AllocationOverflow` so telemetry can tell "input was
    /// preposterously long" apart from "matrix arithmetic overflowed".
    #[error("BAQ query or reference exceeds i32 range")]
    SequenceTooLong,

    /// `BaqConfig::band_half_width` was negative. htslib's `set_u`
    /// macro produces wrong DP-cell offsets in that case; surface as
    /// a structured error rather than corrupting the matrix.
    #[error("BAQ band_half_width must be non-negative")]
    InvalidBandwidth,

    /// The forward / backward matrix would not fit in memory. htslib
    /// signals this by setting `errno = ENOMEM`; we return a structured
    /// error so the per-read driver can count it as a skip reason.
    #[error("BAQ matrix allocation would overflow")]
    AllocationOverflow,
}
