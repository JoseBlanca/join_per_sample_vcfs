//! `AlignmentFile::open` — the validate-on-open gate, and the handle it
//! produces.
//!
//! Every check lives inside the function that opens the file, so a file is
//! either opened *and* validated or it is an `Err`: there is no window in which
//! an unvalidated handle exists. The checks run fail-fast in this order —
//! `@HD SO`, `@SQ`↔reference, the index, `@RG SM`
//! (`doc/devel/ng/spec/alignment_file.md` §3.1).
//!
//! **The invariant this establishes** is what every later layer leans on: once
//! the gate passes, `ref_id == ContigId` holds by construction for every record
//! the file can yield. That is what makes it sound for the merge one layer up
//! to compare positions across files without re-checking anything.
//!
//! ng reads `@HD SO`, `@SQ` and `@RG SM` off the noodles `sam::Header` itself
//! rather than reusing production's extractors: those are module-private, and
//! making them visible would be a production edit the ng freeze forbids
//! (`doc/devel/ng/arch/alignment_file.md` §5).
//!
//! The name says BAM but the module opens CRAM too — "BAM" in the everyday
//! sense of "the alignment file" (spec §6).

#[cfg(test)]
mod tests {}
