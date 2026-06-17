//! Stage 1 — `ssr-pileup` (Mark-2): per-sample read fetch + per-read
//! delimitation + per-locus observed-sequence tally → `.ssr.psp`.
//!
//! Built so far: [`footprint`] (read-vs-locus geometry + the cheap admission
//! gate) and [`fetch_reads`] (the reservoir depth cap + per-locus fetch). Still
//! to land: `alignment` (the delimiter + quality gate), `locus_tally`, and the
//! driver.

pub mod alignment;
pub mod fetch_reads;
pub mod footprint;
