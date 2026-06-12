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
//! Design docs: [`doc/devel/specs/ssr_genotyping.md`] (the model + why) and
//! [`doc/devel/architecture/`] (`ssr_genotyping_architecture.md`,
//! `ssr_shared_types.md`, `ssr_catalog.md`).
//!
//! **Build status.** Only the shared types (`types`) are present so far; the
//! stage modules are added in data-flow order as they are built.

pub mod types;
