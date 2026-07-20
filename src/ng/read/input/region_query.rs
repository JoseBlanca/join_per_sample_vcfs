//! Serving one region: the BAM and CRAM region-query record sources, and the
//! order guard over the filtered stream.
//!
//! Both sources are index-driven siblings of the whole-file
//! [`BamRecordSource`](crate::ng::read::filtering::BamRecordSource) /
//! [`CramRecordSource`](crate::ng::read::filtering::CramRecordSource): they
//! query the **already-parsed** index for candidate chunks, seek a pooled
//! reader, and drop records the index over-returned — chunk-edge slop, which is
//! a reader concern and so is dropped **uncounted**, never charged to a filter
//! drop reason. There is no new trait: BAM and CRAM are two containers for one
//! idea, not competing implementations to bake off
//! (`doc/devel/ng/arch/alignment_file.md` §4).
//!
//! Two things here fail *quietly* rather than loudly, which is why each is
//! built and committed on its own with an independent oracle:
//!
//! - **A missed chunk edge is a wrong genotype, not a crash.** The region
//!   query's oracle is the existing whole-file source: an indexed query must
//!   return exactly what a full linear scan filtered to the same region returns
//!   (spec §7, T5).
//! - **A guard that never fires looks exactly like a guard that works.** The
//!   order check is mutation-verified — removing it must let a planted
//!   regression through (T4a).
//!
//! The order guard's state lives in the per-region iterator, never on the
//! handle, so querying region B and then region A is a new forward scan rather
//! than a spurious regression (spec §3.2).

#[cfg(test)]
mod tests {}
