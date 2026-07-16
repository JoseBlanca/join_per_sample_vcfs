//! ng step 3 — the typed-region generator: walk the reference and cut it into
//! consecutive typed regions, each a span plus *what the sequence there is*.
//! Design: `doc/devel/ng/spec/typed_regions.md` (spec) and
//! `doc/devel/ng/arch/typed_regions.md` (types & interfaces).
//!
//! **Build status (incremental).** Milestone A only: [`admission`] — ng's own
//! copy of the STR admission policy. The walk's types (`TypedRegion`,
//! `RegionKind`, `TypedRegionConfig`, …), `GenomeRegions`, and
//! `TypedRegionIterator` land in Milestones C–E.
//!
//! **A folder, not a file, and not because of a bake-off** (there is none —
//! spec §6). The admission port is a second concern with its own dense test
//! suite, so it gets its own module beside the walk.
//!
//! ## Production is frozen; ng owns its copies
//!
//! Step 3 needs an STR admission policy that is windowed, 1-based/`u64`,
//! driven by `RepeatInterval`s, all-knobs, and that hands bundle members back
//! instead of dropping them. `ssr::catalog::postprocess::build_loci` is none of
//! those things, and **reshaping it in place is not on the table** (spec
//! Revision 2026-07-16, owner): production stays exactly as it is, so that it
//! remains an *independent yardstick* for the experiments ng exists to run.
//!
//! So [`admission`] is a **port**: the logic transcribed unchanged, the shape
//! ng's. What sharing one function used to guarantee for free, a test now pins
//! — see [`admission`]'s differential against production (spec §8.0).

pub mod admission;
