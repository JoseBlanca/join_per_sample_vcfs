//! ng — the next-generation, step-decomposed variant/STR caller: a single-phase,
//! in-memory research lab (see `doc/devel/ng/`). Winning steps are ported back into the
//! production two-phase engine. This is the first landed piece — the shared type
//! vocabulary and the `RefSeq` reference-sequence accessor.

pub mod ref_seq;
pub mod types;

pub use ref_seq::{InMemoryRefSeq, RawRefSeq, RefSeq, RefSeqError, ResidentRefSeq, WindowedRefSeq};
pub use types::ContigId;
