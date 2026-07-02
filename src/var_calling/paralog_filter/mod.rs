//! The hidden-paralog filter's var-calling wiring (Milestone S).
//!
//! The pure statistics core lives in [`crate::paralog`] (Q1–Q5) and depends on
//! nothing here; this module is the inward-facing integration that plugs it
//! into the cohort caller `crate::var_calling::pipeline::run_var_calling` and
//! turns each candidate variant into a keep/drop VCF verdict. The design is
//! settled in `doc/devel/architecture/hidden_paralog_varcalling_wiring.md` and
//! built in the order of
//! `doc/devel/implementation_plans/paralog_varcalling_wiring.md`.
//!
//! Because a locus's verdict depends on every other locus (`Hexp` needs all
//! allele frequencies; `π` + the FDR cut need all LRs), the filter cannot fit
//! in the current single streaming pass. When it is on, the run becomes:
//! spill every scored locus to an ephemeral file while accumulating the global
//! quantities, calibrate once at the end, then read the spill back to write the
//! surviving calls — all with memory flat in both variant count and sample
//! count.
//!
//! Build order (each piece standalone-tested before the S6 integration):
//!
//! - **S1 (this file's [`prepass`]):** the up-front per-sample state
//!   ([`ParalogPrePass`]) plus the running [`HexpAccumulator`] for the one
//!   global quantity that comes from the caller pass.
//! - S2 `spill` — the ephemeral per-locus store (written once, read twice).
//! - S3 — inline per-window coverage folded into the existing caller pass.
//! - S4 `calibrate` — spill → LR → histogram → `π` + FDR cut.
//! - S5 — the write pass (spill → recompute LR → apply cut → VCF).
//! - S6 — orchestration + CLI.

// The two-pass flow is fully wired into `run_var_calling` (S6): the pre-pass,
// the record + window spills, the calibrate pass, and the write pass. The window
// modules (S3/S6c) feed the window spill that the read passes join by tile key.
pub(crate) mod calibrate;
pub(crate) mod prepass;
pub(crate) mod spill;
#[cfg(test)]
pub(crate) mod test_support;
pub(crate) mod window_coverage;
pub(crate) mod window_gc;
pub(crate) mod window_producer;
pub(crate) mod write_pass;

/// Default `--paralog-fdr`: the target false-discovery rate for the
/// hidden-paralog filter. ≈ 1 % (introgression-safe), so the filter is **on by
/// default**. `0.0` (or `--no-paralog-filter`) disables it. Pinned here pending
/// the T1 flagged-set profile on tomato2.
pub(crate) const DEFAULT_PARALOG_FDR: f64 = 0.01;
