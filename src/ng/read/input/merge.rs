//! The argmin k-way merge: a sample's per-file streams combined into one
//! coordinate-ordered stream, plus the same-file-twice check.
//!
//! A linear argmin over a handful of files rather than a binary heap — for the
//! k values that occur here a linear scan of a small contiguous array beats
//! `O(log k)` on constants and locality. Ties break to the **lowest file
//! index**, which is what makes output order reproducible; because a sample's
//! files are usually separate experiments covering the same coordinate range,
//! ties are routine rather than incidental
//! (`doc/devel/ng/spec/sample_reads.md` §3.2).
//!
//! **Every read of every sample passes through here, so the per-read budget is
//! the thing to protect.** The comparisons are not the worry — k small-integer
//! compares. Two other costs are, and the layout exists to avoid both:
//!
//! - **Compare keys, not reads.** The head keys live in an array *beside* the
//!   head slots, refreshed only when a head is refilled, so the argmin scans a
//!   few contiguous integers instead of chasing a pointer into each read.
//! - **Move each read exactly once, and never clone it.** A `MappedRead` owns
//!   its sequence, qualities and CIGAR — it is the big object in this pipeline.
//!   A `clone()` in this loop is a defect, not a slow path, and no correctness
//!   test would catch it, which is why the budget gets its own benchmark (T14).
//!
//! The same-file-twice check rides on a comparison the argmin already made: two
//! reads at different positions cannot be the same read, so it runs **only on a
//! tie**, and within a tie compares `flag` before `qname` so the string compare
//! happens only for reads that already agree on position and flags.
//!
//! The merge is deliberately **concrete, not generic**. The one plausible
//! second caller — the cohort layer — probably wants *group-by-position* rather
//! than *interleave*, so the shared abstraction cannot yet be identified; and
//! the same-file-twice check does not generalise. If it is generalised later,
//! the shape is a key extractor `Fn(&T) -> K`, never a trait on the item
//! (spec §5).

#[cfg(test)]
mod tests {}
