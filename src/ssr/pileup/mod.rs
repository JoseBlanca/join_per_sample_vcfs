//! Stage 1 — `ssr-pileup` (Mark-2): per-sample read fetch + per-read
//! delimitation + per-locus observed-sequence tally → `.ssr.psp`.
//!
//! - [`footprint`] — read-vs-locus geometry, the cheap reach gate, region extraction.
//! - [`fetch_reads`] — per-locus index query + reservoir depth cap.
//! - [`alignment`] — the Viterbi+traceback delimiter and the Q1 quality gate.
//! - [`locus_tally`] — fold per-read outcomes into the per-locus observed ladder.
//! - [`driver`] — config, header build, and the parallel catalog-walk run loop.

pub mod alignment;
pub mod driver;
pub mod fetch_reads;
pub mod footprint;
pub mod locus_tally;
