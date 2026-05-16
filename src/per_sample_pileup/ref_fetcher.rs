//! Production `RefSeqFetcher` impl that bounds memory by evicting
//! the previous chromosome's bases when the walker advances to a
//! new one.
//!
//! Background: `noodles_fasta::Repository` keeps a `HashMap<Vec<u8>,
//! Arc<Sequence>>` cache that has no eviction policy. A naive impl
//! that builds one Repository at startup and shares it across
//! chromosomes would accumulate the bases of every visited contig
//! for the lifetime of the run — roughly 3 GB across the 24 human
//! chromosomes. The walker is sequential single-pass per chromosome
//! (per `ia/specs/pileup_walker.md` §"Chromosome boundaries"), so
//! steady-state we only need the *current* contig in memory; the
//! previous one's `Arc<Sequence>` can be dropped the moment the
//! walker advances. See finding S6 in
//! `ia/reviews/pileup_samtools_comparison_2026-05-07.md`.

use std::cell::Cell;
use std::io;
use std::path::Path;

use noodles_fasta as fasta;

use super::cram_input::ContigList;
use super::pileup::RefSeqFetcher;

/// `RefSeqFetcher` that holds a single `noodles_fasta::Repository`
/// for the whole run and clears its cache whenever the walker moves
/// to a new chromosome. Steady-state cache size is exactly one
/// contig regardless of how many chromosomes have been visited.
pub struct ChromBoundaryRefFetcher {
    repository: fasta::Repository,
    contigs: ContigList,
    /// Last `chrom_id` seen by `fetch`. `Cell` (not `RefCell`) is
    /// enough because it only holds a `Copy` `u32`. The walker is
    /// single-threaded per call (see S6 §"Risk"), so non-`Sync`
    /// interior mutability is fine.
    current_chrom: Cell<Option<u32>>,
}

impl ChromBoundaryRefFetcher {
    /// Open the `.fa` (and its sibling `.fai`) at `fasta_path` and
    /// build a fetcher that resolves `chrom_id` → contig name via
    /// `contigs`.
    pub fn new(fasta_path: &Path, contigs: ContigList) -> io::Result<Self> {
        let indexed_reader =
            fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
        let adapter = fasta::repository::adapters::IndexedReader::new(indexed_reader);
        let repository = fasta::Repository::new(adapter);
        Ok(Self {
            repository,
            contigs,
            current_chrom: Cell::new(None),
        })
    }

    /// Number of contigs currently held in the underlying
    /// `Repository`'s cache. With the chromosome-boundary eviction
    /// in place, this is 0 before the first fetch and 1 thereafter.
    /// Exposed so tests can pin the eviction invariant.
    pub fn cached_contig_count(&self) -> usize {
        self.repository.len()
    }
}

impl RefSeqFetcher for ChromBoundaryRefFetcher {
    fn fetch(&self, chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        if self.current_chrom.get() != Some(chrom_id) {
            // Walker moved to a new chromosome (or this is the
            // first fetch of the run). Drop the previous contig's
            // bases before loading the next.
            self.repository.clear();
            self.current_chrom.set(Some(chrom_id));
        }

        fetch_from_repository(
            &self.repository,
            &self.contigs,
            chrom_id,
            start_1based,
            length,
        )
    }
}

// ---------------------------------------------------------------------
// SyncRefFetcher — Sync-safe, no eviction. For the BAQ stage.
// ---------------------------------------------------------------------

/// `RefSeqFetcher` variant for the rayon-parallel BAQ stage, which
/// requires a `Sync` fetcher to share across worker threads.
///
/// Trades eviction for thread-safety: under the hood it is the same
/// `noodles_fasta::Repository` (already `Sync` thanks to its
/// internal `Arc<RwLock<...>>`) but without the per-fetch
/// chromosome-boundary check that makes [`ChromBoundaryRefFetcher`]
/// non-`Sync`. The pipeline uses both:
///
/// - BAQ stage → [`SyncRefFetcher`] (parallel; non-evicting cache).
/// - Pileup walker → [`ChromBoundaryRefFetcher`] (sequential;
///   one-chrom-resident cache).
///
/// Memory: this fetcher's cache grows to hold every chromosome the
/// run visits (~3 GB on a 24-chrom human reference). For Stage 1's
/// per-sample whole-CRAM scan that is the dominant memory user, but
/// it is bounded and predictable. A future slice can introduce a
/// thread-safe chrom-boundary fetcher if the doubled cache footprint
/// becomes a problem; today the two-fetcher split is the simpler
/// design.
pub struct SyncRefFetcher {
    repository: fasta::Repository,
    contigs: ContigList,
}

impl SyncRefFetcher {
    pub fn new(fasta_path: &Path, contigs: ContigList) -> io::Result<Self> {
        let indexed_reader =
            fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
        let adapter = fasta::repository::adapters::IndexedReader::new(indexed_reader);
        let repository = fasta::Repository::new(adapter);
        Ok(Self {
            repository,
            contigs,
        })
    }
}

impl RefSeqFetcher for SyncRefFetcher {
    fn fetch(&self, chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        fetch_from_repository(
            &self.repository,
            &self.contigs,
            chrom_id,
            start_1based,
            length,
        )
    }
}

/// Shared body of both fetchers — read a window from a
/// noodles `Repository` after validating `chrom_id` and 1-based
/// coordinates. Promoted to a free function so the two fetchers
/// agree on the error shapes (test invariants checked against
/// `ChromBoundaryRefFetcher` cover `SyncRefFetcher` by reuse).
fn fetch_from_repository(
    repository: &fasta::Repository,
    contigs: &ContigList,
    chrom_id: u32,
    start_1based: u32,
    length: u32,
) -> Result<Vec<u8>, io::Error> {
    let entry = contigs.entries.get(chrom_id as usize).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "chrom_id {chrom_id} out of range (have {} contigs)",
                contigs.entries.len()
            ),
        )
    })?;

    let seq_arc = match repository.get(entry.name.as_bytes()) {
        None => {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("contig {} not in FASTA", entry.name),
            ));
        }
        Some(Err(e)) => return Err(e),
        Some(Ok(seq)) => seq,
    };
    let bytes = seq_arc.as_ref().as_ref();

    if start_1based == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "start_1based must be >= 1",
        ));
    }
    let start_idx = (start_1based - 1) as usize;
    let end_idx = start_idx
        .checked_add(length as usize)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "fetch range overflow"))?;
    if end_idx > bytes.len() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!(
                "fetch [{}, {}) past contig {} length {}",
                start_1based,
                start_1based + length,
                entry.name,
                bytes.len()
            ),
        ));
    }
    // Soft-masked FASTAs (Ensembl/Gencode default) encode repeat
    // regions as lowercase `acgtn`. The mask carries no information
    // used downstream, and the PSP writer rejects anything outside
    // A/C/G/T/N; uppercase here so every consumer sees canonical bases.
    Ok(bytes[start_idx..end_idx]
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect())
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    use crate::per_sample_pileup::cram_files::{ContigSpec, build_fasta};
    use crate::per_sample_pileup::cram_input::{ContigEntry, ContigList};

    fn contig_list(entries: &[(&str, u64)]) -> ContigList {
        ContigList {
            entries: entries
                .iter()
                .map(|(n, l)| ContigEntry {
                    name: (*n).into(),
                    length: *l,
                    md5: None,
                })
                .collect(),
        }
    }

    #[test]
    fn fetch_returns_bases_at_correct_offset() {
        // build_fasta writes 'A' for every base, so the test is
        // about offsets / range handling, not byte values.
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 100,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 100)])).expect("fetcher");

        let bases = fetcher.fetch(0, 5, 4).expect("fetch");
        assert_eq!(bases, b"AAAA");
        assert_eq!(bases.len(), 4);
    }

    #[test]
    fn fetch_uppercases_soft_masked_bases() {
        // Soft-masked FASTAs (Ensembl/Gencode default) encode repeat
        // regions as lowercase `acgtn`. The fetcher must canonicalise
        // to uppercase so the PSP writer (which rejects lowercase
        // bytes) never sees a soft-masked REF byte.
        use std::fs::File;
        use std::io::Write;

        let dir = tempfile::tempdir().expect("tempdir");
        let fasta_path = dir.path().join("soft.fa");
        let fai_path = dir.path().join("soft.fa.fai");

        let header = b">chr0\n";
        let seq = b"ACGTacgtNn";
        {
            let mut fa = File::create(&fasta_path).expect("fa");
            fa.write_all(header).expect("hdr");
            fa.write_all(seq).expect("seq");
            fa.write_all(b"\n").expect("nl");
        }
        {
            let mut fai = File::create(&fai_path).expect("fai");
            writeln!(
                fai,
                "chr0\t{}\t{}\t{}\t{}",
                seq.len(),
                header.len(),
                seq.len(),
                seq.len() + 1
            )
            .expect("fai");
        }

        let fetcher = ChromBoundaryRefFetcher::new(
            &fasta_path,
            contig_list(&[("chr0", seq.len() as u64)]),
        )
        .expect("fetcher");

        let bases = fetcher.fetch(0, 1, seq.len() as u32).expect("fetch");
        assert_eq!(bases, b"ACGTACGTNN");
    }

    #[test]
    fn cache_evicts_previous_chromosome_on_chrom_change() {
        // The eviction invariant: cache size stays at 1 across
        // chrom changes, regardless of how many chromosomes the
        // walker has visited. Without `Repository::clear()` it
        // would grow unboundedly (a regression S6 was added to
        // prevent).
        let specs = vec![
            ContigSpec {
                name: "chr0".into(),
                length: 100,
            },
            ContigSpec {
                name: "chr1".into(),
                length: 100,
            },
        ];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 100), ("chr1", 100)]))
                .expect("fetcher");

        assert_eq!(
            fetcher.cached_contig_count(),
            0,
            "fresh fetcher should have an empty cache"
        );

        fetcher.fetch(0, 1, 4).expect("fetch chr0");
        assert_eq!(
            fetcher.cached_contig_count(),
            1,
            "chr0 must be cached after first fetch"
        );

        fetcher.fetch(1, 1, 4).expect("fetch chr1");
        assert_eq!(
            fetcher.cached_contig_count(),
            1,
            "chr0 must have been evicted when fetcher saw chr1; \
             cache should hold the single current contig only"
        );

        // A second fetch on the same chromosome must reuse the
        // cached entry — no new eviction, no growth.
        fetcher.fetch(1, 5, 4).expect("re-fetch chr1");
        assert_eq!(fetcher.cached_contig_count(), 1);

        // Going back to chr0 also re-loads (chr0 was evicted) and
        // does not double-cache.
        fetcher.fetch(0, 1, 4).expect("re-fetch chr0");
        assert_eq!(fetcher.cached_contig_count(), 1);
    }

    #[test]
    fn fetch_past_contig_end_returns_unexpected_eof() {
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 10)])).expect("fetcher");

        let err = fetcher.fetch(0, 8, 5).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn unknown_chrom_id_returns_invalid_input() {
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 10)])).expect("fetcher");

        let err = fetcher.fetch(99, 1, 4).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    #[test]
    fn start_1based_zero_is_rejected() {
        // 1-based coordinates: 0 is a contract violation, not a
        // valid offset.
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 10)])).expect("fetcher");

        let err = fetcher.fetch(0, 0, 1).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    #[test]
    fn contig_name_in_fasta_must_match_contig_list() {
        // contig_list says the chrom is "chr0", but the FASTA on
        // disk only knows "different_name". The Repository surfaces
        // the missing contig as an `io::Error`; the exact kind is
        // noodles' choice (currently `InvalidInput` from the
        // indexed-reader adapter), so we assert only that fetch
        // errors out — not which kind.
        let specs = vec![ContigSpec {
            name: "different_name".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 10)])).expect("fetcher");

        fetcher.fetch(0, 1, 4).expect_err("must fail");
    }
}
