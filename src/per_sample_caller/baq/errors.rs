//! Errors produced by the BAQ stage.
//!
//! See [`ia/feature_implementation_plans/baq.md`](../../../ia/feature_implementation_plans/baq.md)
//! §"Algorithm port" for how each variant maps onto htslib's failure modes.

use thiserror::Error;

/// Reasons `probaln_glocal` returns early without producing a usable
/// alignment. Mirrors the `INT_MIN` / `errno` exits in
/// [htslib/probaln.c:85-107](../../../htslib/probaln.c#L85-L107).
#[derive(Error, Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaqOverflow {
    /// One of the input arrays had a length the algorithm cannot handle
    /// (e.g. `iqual.len() != query.len()` or `query.len() >= i32::MAX - 2`).
    #[error("BAQ inputs are invalid (length mismatch or out of i32 range)")]
    InvalidInput,

    /// The forward / backward matrix would not fit in memory. htslib
    /// signals this by setting `errno = ENOMEM`; we return a structured
    /// error so the per-read driver can count it as a skip reason.
    #[error("BAQ matrix allocation would overflow")]
    AllocationOverflow,
}
