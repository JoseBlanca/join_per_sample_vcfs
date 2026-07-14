//! ng — the next-generation, step-decomposed variant/STR caller: a single-phase,
//! in-memory research lab (see `doc/devel/ng/`). Winning steps are ported back into the
//! production two-phase engine. Landed so far: the shared type vocabulary
//! ([`types`]), the `RefSeq` reference-sequence accessor ([`ref_seq`]), and the step-1
//! read-filtering module ([`read`]), and the tandem-repeat scanner primitive
//! ([`tandem_repeat`], a shared sequence primitive — types-and-scaffold stage).

pub mod read;
pub mod ref_seq;
pub mod tandem_repeat;
pub mod types;

pub use ref_seq::{InMemoryRefSeq, RawRefSeq, RefSeq, RefSeqError, ResidentRefSeq, WindowedRefSeq};
pub use types::{BaseQual, Bp, ContigId, DomainError, MapQual, MismatchFraction};
