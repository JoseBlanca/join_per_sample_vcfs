//! Writer — section 3 (appendix §E).
//!
//! *(today: `var_calling::driver` `drive_blocks_parallel`'s receive loop)*
//!
//! `VcfWriter` consumes `CalledChunk`s and writes VCF, **single consumer**,
//! reordering by `chunk_order` via a `BTreeMap` reorder buffer drained on the
//! *next expected* `chunk_order`.
//!
//! **Gapless invariant:** every chunk yields exactly one `CalledChunk` (empty
//! ones included), or the drain stalls; a debug-assert checks the consumed
//! `chunk_order` sequence is contiguous.
//!
//! Phase 4 builds this module.

// Phase 4.
