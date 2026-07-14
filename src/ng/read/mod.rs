//! Read-level pipeline stages — turning a decoded alignment record into locus
//! evidence. Two steps share this module because they share the same reference
//! accessor and read filtering's output (`MappedRead`) is read preparation's
//! input:
//!
//! - **step 1 — [`filtering`]** — the fixed whole-read keep/drop prelude (this
//!   milestone). A single file, no bake-off: it is a fixed prelude with no
//!   competing implementations.
//! - **step 2 — read preparation (`ReadPrep`)** — realign/delimit a kept read
//!   against the reference; its swappable implementations land here as siblings
//!   (`trust_mapper.rs`, `pair_hmm.rs`, …) in a later plan.
//!
//! This is the deliberate, documented deviation from the "one folder per step"
//! rule in `doc/devel/ng/arch/module_layout.md` principle 1 — see
//! `doc/devel/ng/spec/read_filtering.md` §2.1.

pub mod filtering;
