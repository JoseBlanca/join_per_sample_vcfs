//! ng — the next-generation, step-decomposed variant/STR caller: a single-phase,
//! in-memory research lab (see `doc/devel/ng/`). Winning steps are ported back into the
//! production two-phase engine. Landed so far: the shared type vocabulary
//! ([`types`]); the `RefSeq` reference-sequence accessor ([`ref_seq`]); the step-1
//! read-filtering module ([`read`]); the tandem-repeat scanner primitive
//! ([`tandem_repeat`], a shared sequence primitive — types-and-scaffold stage);
//! step 3's typed-region generator ([`region_typing`] — Milestone A, the
//! segment-criteria port); and the read-alignment module ([`alignment`] — a
//! folder of competing aligners, not a pipeline step).
//!
//! **Production is frozen.** ng is a from-scratch caller: it does not edit
//! `src/ssr/` or `src/regions.rs`, and it does not depend on trf-mod. Where ng
//! needs a production behaviour in a different shape, it **copies the code and
//! changes its own version** (owner, 2026-07-16); reuse is for what costs
//! production nothing. Winning steps are ported back only after the experiments
//! ng exists to run have decided something.

#[cfg(test)]
mod scanner_parity;

pub mod alignment;
pub mod locus_generation;
pub mod raw_chrom_reader;
pub mod read;
pub mod ref_seq;
pub mod reference_info;
pub mod region_typing;
pub mod tandem_repeat;
pub mod types;

pub use ref_seq::{
    ContigTable, EvictableRefSeq, InMemoryRefSeq, RawRefSeq, RefSeq, RefSeqError, ResidentRefSeq,
    WindowedRefSeq,
};
pub use types::{BaseQual, Bp, ContigId, DomainError, MapQual, MismatchFraction};
