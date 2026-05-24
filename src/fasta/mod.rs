//! Reference FASTA I/O — streaming and manual-evict fetchers plus the
//! `MultiChromRefFetcher` trait that abstracts "give me reference
//! bases at a given chromosome and position."
//!
//! This module is a peer of the pipeline-stage modules. Every stage
//! that needs the reference (BAQ, the pileup walker, the DUST filter
//! on the cohort side, the per-group merger when extracting REF
//! sequence for compound alleles) reaches in here.
//!
//! - [`fetcher::StreamingChromRefFetcher`] /
//!   [`fetcher::ManualEvictChromRefFetcher`] /
//!   [`fetcher::MultiChromStreamingRefFetcher`] — the production
//!   implementations; see [`fetcher`] for the per-impl tradeoffs.
//! - [`fetcher::ChromRefFetcher`] — single-chromosome trait
//!   (sliding-buffer contract, used by single-chrom consumers).
//! - [`MultiChromRefFetcher`] — multi-chromosome trait used by the
//!   pileup walker.
//! - [`ChromRefFetchError`] — typed errors raised by every fetcher.
//! - [`ContigEntry`] / [`ContigList`] — the reference's contig table
//!   (name, length, optional MD5). Shared across the FASTA, CRAM, and
//!   PSP formats because every input asserts which reference it is
//!   against.

pub mod fetcher;

pub use fetcher::{
    ChromRefBaseIter, ChromRefFetchError, ChromRefFetcher, ManualEvictChromRefFetcher,
    MultiChromStreamingRefFetcher, STREAMING_REF_BUFFER_BYTES, StreamingChromRefFetcher,
};

// ---------------------------------------------------------------------
// ContigList
// ---------------------------------------------------------------------

/// One reference sequence (`@SQ`) entry: name, length, and optional MD5.
#[derive(Debug, Clone)]
pub struct ContigEntry {
    pub name: String,
    pub length: u64,
    pub md5: Option<[u8; 16]>,
}

impl PartialEq for ContigEntry {
    fn eq(&self, other: &Self) -> bool {
        if self.name != other.name || self.length != other.length {
            return false;
        }
        match (self.md5, other.md5) {
            (Some(a), Some(b)) => a == b,
            // Absent MD5 acts as a wildcard: a CRAM that omits M5 is
            // not contradicting a CRAM that carries one.
            (None, _) | (_, None) => true,
        }
    }
}

impl Eq for ContigEntry {}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ContigList {
    pub entries: Vec<ContigEntry>,
}

impl ContigList {
    /// First-difference report between two lists, intended for error
    /// messages. Returns `Ok(())` when the lists agree; otherwise an
    /// `Err` carrying a short string identifying which field
    /// disagreed and where.
    pub(crate) fn first_disagreement(&self, other: &Self) -> Result<(), String> {
        if self.entries.len() != other.entries.len() {
            return Err(format!(
                "@SQ list length differs ({} vs {})",
                self.entries.len(),
                other.entries.len()
            ));
        }
        for (i, (a, b)) in self.entries.iter().zip(other.entries.iter()).enumerate() {
            if a.name != b.name {
                return Err(format!(
                    "name disagreement at index {} ('{}' vs '{}')",
                    i, a.name, b.name
                ));
            }
            if a.length != b.length {
                return Err(format!(
                    "length disagreement at index {} (contig '{}': {} vs {})",
                    i, a.name, a.length, b.length
                ));
            }
            if let (Some(ma), Some(mb)) = (a.md5, b.md5)
                && ma != mb
            {
                return Err(format!(
                    "md5 disagreement at index {} (contig '{}')",
                    i, a.name
                ));
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------
// MultiChromRefFetcher trait
// ---------------------------------------------------------------------

/// Multi-chromosome reference-FASTA fetcher used by the Stage 1
/// pileup walker. Errors are typed ([`ChromRefFetchError`]) so
/// callers can route I/O failures, range failures, and contract
/// violations distinctly. Single-chromosome consumers (DUST, BAQ,
/// PerGroupMerger) should use [`ChromRefFetcher`] instead — it
/// drops the `chrom_id` parameter and exposes a sliding-buffer
/// contract specifically for monotonic-forward access.
pub trait MultiChromRefFetcher {
    /// Fetch `length` reference bases starting at the 1-based
    /// position `start` on chromosome `chrom_id`. Bytes are
    /// uppercase ASCII over `{A,C,G,T,N}` (canonicalised by the
    /// fetcher implementation).
    ///
    /// # Errors
    ///
    /// - [`ChromRefFetchError::OutOfBounds`] if the requested
    ///   window exceeds the chromosome length.
    /// - [`ChromRefFetchError::InvalidStart`] if `start_1based == 0`.
    /// - [`ChromRefFetchError::Io`] on any underlying FASTA I/O
    ///   failure or unknown `chrom_id`.
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, ChromRefFetchError>;

    /// Forward sequential iterator over every uppercased base of
    /// `chrom_id`'s contig, in 1..=`length` order. Used by the DUST
    /// mask construction (one pass per chrom) to avoid materialising
    /// the whole contig as a single `Vec<u8>`.
    ///
    /// Default impl materialises via `fetch(chrom_id, 1, length)`;
    /// streaming fetchers override to walk a sliding buffer instead.
    /// The boxed iterator costs one heap allocation per chrom plus
    /// one virtual dispatch per byte — the latter is dominated by
    /// the inner sdust scoring work in practice.
    ///
    /// # Errors
    ///
    /// Same failure modes as [`Self::fetch`] (the default impl
    /// just calls `fetch(chrom_id, 1, length)`).
    fn iter_bases<'a>(
        &'a self,
        chrom_id: u32,
        length: u32,
    ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>
    {
        let seq = self.fetch(chrom_id, 1, length)?;
        Ok(Box::new(seq.into_iter().map(Ok)))
    }
}

/// Forwarding impl so callers may pass either an owned fetcher or a
/// shared reference into [`PileupWalker::new`] / [`run`]. The walker
/// only ever calls `&self` methods, so the borrow is sufficient.
///
/// [`PileupWalker::new`]: crate::per_sample_pileup::pileup::PileupWalker::new
/// [`run`]: crate::per_sample_pileup::pileup::run
impl<T: MultiChromRefFetcher + ?Sized> MultiChromRefFetcher for &T {
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, ChromRefFetchError> {
        (**self).fetch(chrom_id, start_1based, length)
    }

    fn iter_bases<'a>(
        &'a self,
        chrom_id: u32,
        length: u32,
    ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>
    {
        (**self).iter_bases(chrom_id, length)
    }
}
