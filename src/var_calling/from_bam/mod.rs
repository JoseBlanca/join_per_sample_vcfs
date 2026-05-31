//! Direct CRAM→VCF path (single-sample, streaming per-record).
//!
//! - [`pipeline`] — the streaming per-position pipeline
//!   (`drive_cohort_pipeline`); also the from_psp byte-identity oracle.
//! - [`driver`] — the per-chromosome from-bam engine.

pub mod driver;
pub mod pipeline;
