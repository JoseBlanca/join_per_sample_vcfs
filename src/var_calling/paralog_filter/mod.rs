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

// S1–S5 build standalone pieces (each unit-tested against fixtures) before S6
// wires them into `run_var_calling`; until that consumer lands, their public
// surface is dead to the non-test crate. The allow is removed at S6.
#[allow(dead_code)]
pub(crate) mod prepass;
