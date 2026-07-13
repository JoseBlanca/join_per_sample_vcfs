//! Shared ng vocabulary — the domain newtypes cross-step code speaks. It starts as this
//! one file and splits into concept modules (`units`, `locus`, …) as clusters grow
//! (`doc/devel/ng/arch/module_layout.md` principle 3). Seeded here with only what the
//! `RefSeq` reference accessor needs.

/// Which reference sequence a coordinate refers to: an index into the reference contig
/// table ([`crate::fasta::ContigList`]), in `@SQ` / `.fai` order. Unconstrained — any
/// `u32` is a legal index at the type level, and an out-of-range id is caught at fetch
/// time — so the field is public and there is no checked constructor.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct ContigId(pub u32);

impl ContigId {
    #[inline]
    pub fn get(self) -> u32 {
        self.0
    }
}
