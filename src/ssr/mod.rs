//! SSR/STR genotyping — the project's second, independent variant caller.
//!
//! Genotypes short tandem repeats (microsatellites, period ≤ 6) from aligned
//! reads, producing per-individual repeat-allele lengths across a cohort. It is
//! a *standalone* pipeline that shares only low-level alignment I/O and a few
//! numerical kernels with the SNP caller — never its records or its math. The
//! three stages mirror the SNP path's `pileup → .psp → var-calling` shape:
//!
//! ```text
//! reference ─► ssr-catalog ─► catalog ─► ssr-pileup ─► .ssr.psp ─► ssr-call ─► VCF
//! ```
//!
//! Design docs: `doc/devel/specs/ssr_genotyping.md` (the model + why) and
//! `doc/devel/architecture/` (`ssr_genotyping_architecture.md`,
//! `ssr_shared_types.md`, `ssr_catalog.md`).
//!
//! **Build status.** Only the shared types (`types`) are present so far; the
//! stage modules are added in data-flow order as they are built.

// The shared types land before their consumers: the catalog/pileup/call stages
// that read them are built in later phases. Until then the `pub(crate)` surface
// has no in-crate caller, so suppress `dead_code` here rather than scatter
// per-item `#[allow]`s. Remove once `ssr-catalog` (Stage 0) wires these up.
#![allow(dead_code)]

pub mod catalog;
pub mod pileup;
pub mod types;
