//! Multi-way per-position merger over per-sample `.psp` iterators.
//!
//! Given N iterators that each yield `PileupRecord`s in genomic order
//! (the contract `PspReader::records` already enforces per file), the
//! [`PerPositionMerger`] yields one [`PerPositionPileups`] item per
//! `(chrom_id, pos)` any sample covers, with an N-slot vector saying
//! which samples had a record there.
//!
//! Algorithm: linear-scan k-way merge. For each output item the
//! merger does an `O(N)` scan over peeked heads to find the min
//! `(chrom_id, pos)`, pulls successors only from the readers tied
//! at that min, and emits the resulting per-sample slot vector.
//! See `doc/devel/implementation_plans/multi_way_per_position_iterator.md`
//! for the rationale (in particular the "linear vs heap" decision —
//! the WGS access pattern this caller targets has `k ≈ N` at nearly
//! every output position, where the linear scan wins).
//!
//! The merger is generic on the iterator type so synthetic
//! `std::iter`-based fixtures can drive it in unit tests without an
//! on-disk `.psp` round-trip. The chromosome-id-space agreement
//! across input `.psp` files is a *caller-side* precondition; the
//! standalone [`check_chromosome_agreement`] helper performs that
//! check on a slice of opened readers and is what the path-based
//! caller uses before instantiating the merger.

use std::io::{Read, Seek};

use thiserror::Error;

use crate::per_sample_pileup::psp::header::ParsedChromosome;
use crate::per_sample_pileup::psp::{PspReadError, PspReader};
use crate::pileup_record::PileupRecord;

/// One emitted item: a single `(chrom_id, pos)` with a per-sample
/// slot for every sample passed to [`PerPositionMerger::new`].
/// Samples without a record at this position carry `None`.
#[derive(Debug, Clone, PartialEq)]
pub struct PerPositionPileups {
    pub chrom_id: u32,
    /// 1-based anchor position. Matches [`PileupRecord::pos`].
    pub pos: u32,
    /// Indexed by sample order (same as
    /// [`PerPositionMerger::sample_names`]).
    /// `per_sample.len() == PerPositionMerger::n_samples()`.
    pub per_sample: Vec<Option<PileupRecord>>,
}

/// Errors the merger and its companion helpers can emit. Carries
/// enough context (`sample_idx`, `sample_name`, coordinates) to
/// point at the offending input without forcing the caller to keep
/// a side table.
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum PerPositionMergerError {
    /// `readers` and `sample_names` had different lengths. The
    /// caller must hand one name per reader, in the same order.
    /// Defensive check — the rest of the merger's API assumes the
    /// two vectors are length-aligned.
    #[error("readers ({n_readers}) and sample_names ({n_sample_names}) length mismatch")]
    SampleCountMismatch {
        n_readers: usize,
        n_sample_names: usize,
    },

    /// One of the per-sample iterators yielded an error, either on
    /// the initial prefetch in [`PerPositionMerger::new`] or on a
    /// later refill during iteration. The merger latches `done`
    /// after this so subsequent `next()` calls return `None`.
    ///
    /// `source` is boxed because `PspReadError` is the largest enum
    /// variant in this crate; inlining it would bloat
    /// `Result<_, PerPositionMergerError>` past 128 bytes for every per-position
    /// emission, which the merger's hot path returns once per output
    /// item.
    #[error("per-sample reader {sample_idx} ({sample_name}) failed: {source}")]
    PerSampleReader {
        sample_idx: usize,
        sample_name: String,
        #[source]
        source: Box<PspReadError>,
    },

    /// A reader emitted a record whose `(chrom_id, pos)` is `≤` the
    /// last emitted output key. `PspReader::records` already
    /// enforces per-reader monotonicity internally, so in production
    /// this fires only on pathological mocks or on drift in a future
    /// reader variant that relaxes the per-reader guarantee.
    ///
    /// Both the regressing key and the boundary it violated are
    /// carried so the error is interpretable on its own without
    /// reproducing the merger state.
    #[error(
        "reader {sample_idx} ({sample_name}) emitted out-of-order record: \
         {regressing_chrom_id}:{regressing_pos} ≤ last emitted \
         {last_emitted_chrom_id}:{last_emitted_pos}"
    )]
    OutOfOrder {
        sample_idx: usize,
        sample_name: String,
        regressing_chrom_id: u32,
        regressing_pos: u32,
        last_emitted_chrom_id: u32,
        last_emitted_pos: u32,
    },

    /// Two `.psp` files disagree on the per-chromosome metadata for a
    /// specific chromosome — name, length, or md5 differs. Emitted by
    /// [`check_chromosome_agreement`]. `detail` is a human-readable
    /// description of the specific divergence; downstream code
    /// matches on `ChromosomeMismatch { .. }` rather than parsing it.
    ///
    /// See also [`PerPositionMergerError::ChromosomeCountMismatch`]
    /// for the case where the two files disagree on the *count* of
    /// chromosomes (no specific chrom_id is meaningful then).
    #[error(
        "sample {sample_idx} ({sample_name}) chromosome {chrom_id} disagrees with sample 0: {detail}"
    )]
    ChromosomeMismatch {
        sample_idx: usize,
        sample_name: String,
        chrom_id: u32,
        detail: String,
    },

    /// Two `.psp` files disagree on the *number* of chromosomes in
    /// their headers. Emitted by [`check_chromosome_agreement`]
    /// before any per-chromosome comparison. Carries both counts so
    /// the error is self-contained.
    #[error(
        "sample {sample_idx} ({sample_name}) declares {n_other} chromosomes; \
         sample 0 declares {n_baseline}"
    )]
    ChromosomeCountMismatch {
        sample_idx: usize,
        sample_name: String,
        n_baseline: usize,
        n_other: usize,
    },
}

/// Linear-scan k-way merge over per-sample iterators yielding
/// `Result<PileupRecord, PspReadError>`. Owns the iterators and
/// emits [`PerPositionPileups`] items in strictly increasing
/// `(chrom_id, pos)` order.
pub struct PerPositionMerger<I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    readers: Vec<I>,
    /// Peeked next record per reader. `None` once that reader is
    /// exhausted. Same length as `readers`.
    heads: Vec<Option<PileupRecord>>,
    sample_names: Vec<String>,
    chromosomes: Vec<ParsedChromosome>,
    /// `(chrom_id, pos)` of the last emitted item, or `None` before
    /// the first emission. Used to enforce strict monotonicity at
    /// the merge layer.
    last_emitted: Option<(u32, u32)>,
    /// Latch: once set, all subsequent `next()` calls return `None`.
    /// Set on natural exhaustion *and* on the first error emitted —
    /// matches the walker's terminate-on-first-error contract.
    done: bool,
}

impl<I> std::fmt::Debug for PerPositionMerger<I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Exhaustive destructure so a new field on `PerPositionMerger`
        // fails to compile here instead of silently being dropped from
        // the debug rendering.
        let Self {
            readers: _,
            heads: _,
            sample_names,
            chromosomes,
            last_emitted,
            done,
        } = self;
        f.debug_struct("PerPositionMerger")
            .field("sample_names", sample_names)
            .field("n_chromosomes", &chromosomes.len())
            .field("last_emitted", last_emitted)
            .field("done", done)
            .finish()
    }
}

impl<I> PerPositionMerger<I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    /// Construct a merger over `readers`. `sample_names` and
    /// `chromosomes` are taken from the caller-validated headers
    /// (see [`check_chromosome_agreement`] for the chromosome
    /// precondition). Their length must match `readers.len()` for
    /// `sample_names`; `chromosomes` is metadata only and is not
    /// length-checked against the readers.
    ///
    /// Returns [`PerPositionMergerError::SampleCountMismatch`] if
    /// `sample_names.len() != readers.len()`, or
    /// [`PerPositionMergerError::PerSampleReader`] if any reader fails on its first
    /// `next()` (the initial prefetch). Empty `readers` is allowed
    /// and yields an immediately-exhausted iterator.
    pub fn new(
        mut readers: Vec<I>,
        sample_names: Vec<String>,
        chromosomes: Vec<ParsedChromosome>,
    ) -> Result<Self, PerPositionMergerError> {
        if readers.len() != sample_names.len() {
            return Err(PerPositionMergerError::SampleCountMismatch {
                n_readers: readers.len(),
                n_sample_names: sample_names.len(),
            });
        }
        let mut heads: Vec<Option<PileupRecord>> = Vec::with_capacity(readers.len());
        for (sample_idx, reader) in readers.iter_mut().enumerate() {
            match reader.next() {
                Some(Ok(record)) => heads.push(Some(record)),
                None => heads.push(None),
                Some(Err(source)) => {
                    return Err(PerPositionMergerError::PerSampleReader {
                        sample_idx,
                        sample_name: sample_names[sample_idx].clone(),
                        source: Box::new(source),
                    });
                }
            }
        }
        Ok(Self {
            readers,
            heads,
            sample_names,
            chromosomes,
            last_emitted: None,
            done: false,
        })
    }

    pub fn sample_names(&self) -> &[String] {
        &self.sample_names
    }

    pub fn chromosomes(&self) -> &[ParsedChromosome] {
        &self.chromosomes
    }

    pub fn n_samples(&self) -> usize {
        self.sample_names.len()
    }
}

impl<I> Iterator for PerPositionMerger<I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    type Item = Result<PerPositionPileups, PerPositionMergerError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        // Find the minimum (chrom_id, pos) over all peeked heads.
        let mut min_key: Option<(u32, u32)> = None;
        for record in self.heads.iter().flatten() {
            let key = (record.chrom_id, record.pos);
            if min_key.is_none_or(|m| key < m) {
                min_key = Some(key);
            }
        }
        let Some(min_key) = min_key else {
            // All readers exhausted.
            self.done = true;
            return None;
        };

        // Strict-monotonicity check across the emit stream.
        if let Some(last) = self.last_emitted
            && min_key <= last
        {
            // The "offending reader" is the first one whose head
            // sits at the regressing key.
            let offender = self
                .heads
                .iter()
                .position(|h| h.as_ref().is_some_and(|r| (r.chrom_id, r.pos) == min_key))
                .expect("min_key was derived from a non-None head");
            self.done = true;
            return Some(Err(PerPositionMergerError::OutOfOrder {
                sample_idx: offender,
                sample_name: self.sample_names[offender].clone(),
                regressing_chrom_id: min_key.0,
                regressing_pos: min_key.1,
                last_emitted_chrom_id: last.0,
                last_emitted_pos: last.1,
            }));
        }

        // Build the per-sample slot vector; advance tied readers.
        // Loop over `per_sample` (a local) so clippy doesn't flag a
        // raw `0..n` range — we still need parallel access to
        // `self.heads`, `self.readers`, and `self.sample_names` by
        // the same index.
        let mut per_sample: Vec<Option<PileupRecord>> = vec![None; self.n_samples()];
        for (sample_idx, slot) in per_sample.iter_mut().enumerate() {
            let head_at_min = self.heads[sample_idx]
                .as_ref()
                .is_some_and(|r| (r.chrom_id, r.pos) == min_key);
            if !head_at_min {
                continue;
            }
            *slot = self.heads[sample_idx].take();
            match self.readers[sample_idx].next() {
                Some(Ok(record)) => {
                    self.heads[sample_idx] = Some(record);
                }
                None => {
                    // Reader exhausted; head stays `None`.
                }
                Some(Err(source)) => {
                    self.done = true;
                    return Some(Err(PerPositionMergerError::PerSampleReader {
                        sample_idx,
                        sample_name: self.sample_names[sample_idx].clone(),
                        source: Box::new(source),
                    }));
                }
            }
        }

        self.last_emitted = Some(min_key);
        Some(Ok(PerPositionPileups {
            chrom_id: min_key.0,
            pos: min_key.1,
            per_sample,
        }))
    }
}

/// Verify every reader agrees with `readers[0]` on the chromosome
/// id space — same count, names, lengths, and MD5s — and return a
/// clone of the agreed-upon list. Empty `readers` yields an empty
/// list.
///
/// Returns [`PerPositionMergerError::ChromosomeCountMismatch`] if
/// the headers declare a different number of chromosomes, or
/// [`PerPositionMergerError::ChromosomeMismatch`] on the first
/// per-chromosome divergence (name/length/md5). The chromosome list
/// is the operative contract for the merger's correctness;
/// mismatched reference strings that happen to agree on every
/// per-chromosome field are not fatal here.
pub fn check_chromosome_agreement<R: Read + Seek>(
    readers: &[PspReader<R>],
) -> Result<Vec<ParsedChromosome>, PerPositionMergerError> {
    let Some(baseline_reader) = readers.first() else {
        return Ok(Vec::new());
    };
    let baseline = &baseline_reader.header().chromosomes;
    for (sample_idx, reader) in readers.iter().enumerate().skip(1) {
        let other = &reader.header().chromosomes;
        let sample_name = reader.header().sample.clone();
        if other.len() != baseline.len() {
            return Err(PerPositionMergerError::ChromosomeCountMismatch {
                sample_idx,
                sample_name,
                n_baseline: baseline.len(),
                n_other: other.len(),
            });
        }
        for (chrom_id, (base, this)) in baseline.iter().zip(other.iter()).enumerate() {
            if base.name != this.name {
                return Err(PerPositionMergerError::ChromosomeMismatch {
                    sample_idx,
                    sample_name,
                    chrom_id: chrom_id as u32,
                    detail: format!("name {:?} vs {:?}", this.name, base.name),
                });
            }
            if base.length != this.length {
                return Err(PerPositionMergerError::ChromosomeMismatch {
                    sample_idx,
                    sample_name,
                    chrom_id: chrom_id as u32,
                    detail: format!("length {} vs {}", this.length, base.length),
                });
            }
            if base.md5 != this.md5 {
                return Err(PerPositionMergerError::ChromosomeMismatch {
                    sample_idx,
                    sample_name,
                    chrom_id: chrom_id as u32,
                    detail: format!("md5 {:?} vs {:?}", this.md5, base.md5),
                });
            }
        }
    }
    Ok(baseline.clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::per_sample_pileup::psp::header::{ChromosomeEntry, WriterHeader};
    use crate::per_sample_pileup::psp::test_fixtures::writer_header;
    use crate::per_sample_pileup::psp::writer::PspWriter;
    use crate::pileup_record::{AlleleObservation, AlleleSupportStats};
    use std::io::Cursor;

    type Item = Result<PileupRecord, PspReadError>;
    type TestIter = std::vec::IntoIter<Item>;

    fn rec(chrom_id: u32, pos: u32) -> PileupRecord {
        PileupRecord::new(
            chrom_id,
            pos,
            vec![AlleleObservation::new(
                b"A".to_vec(),
                AlleleSupportStats::default(),
                Vec::new(),
            )],
        )
    }

    fn iter_from(items: Vec<Item>) -> TestIter {
        items.into_iter()
    }

    fn fake_err() -> PspReadError {
        PspReadError::IndexChecksum {
            stored: 0,
            computed: 1,
        }
    }

    fn names(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("S{i}")).collect()
    }

    fn finish_empty_writer(header: WriterHeader) -> Vec<u8> {
        let writer = PspWriter::new(Cursor::new(Vec::new()), header).unwrap();
        writer.finish().unwrap().into_inner()
    }

    fn psp_reader_with_header(header: WriterHeader) -> PspReader<Cursor<Vec<u8>>> {
        let bytes = finish_empty_writer(header);
        PspReader::new(Cursor::new(bytes)).expect("test psp file opens")
    }

    // ----- merger -----

    #[test]
    fn empty_cohort_yields_immediately() {
        let mut merger =
            PerPositionMerger::<TestIter>::new(Vec::new(), Vec::new(), Vec::new()).unwrap();
        assert_eq!(merger.n_samples(), 0);
        assert!(merger.next().is_none());
    }

    #[test]
    fn single_reader_identity_pass_through() {
        let records = vec![rec(0, 10), rec(0, 20), rec(0, 30)];
        let iter = iter_from(records.iter().cloned().map(Ok).collect());
        let mut merger = PerPositionMerger::new(vec![iter], names(1), Vec::new()).unwrap();
        for expected in &records {
            let pileups = merger.next().unwrap().unwrap();
            assert_eq!(pileups.chrom_id, expected.chrom_id);
            assert_eq!(pileups.pos, expected.pos);
            assert_eq!(pileups.per_sample.len(), 1);
            assert_eq!(pileups.per_sample[0].as_ref(), Some(expected));
        }
        assert!(merger.next().is_none());
    }

    #[test]
    fn two_readers_fully_overlapping() {
        let a = iter_from(vec![Ok(rec(0, 1)), Ok(rec(0, 2)), Ok(rec(0, 3))]);
        let b = iter_from(vec![Ok(rec(0, 1)), Ok(rec(0, 2)), Ok(rec(0, 3))]);
        let mut merger = PerPositionMerger::new(vec![a, b], names(2), Vec::new()).unwrap();
        for pos in 1..=3 {
            let pileups = merger.next().unwrap().unwrap();
            assert_eq!((pileups.chrom_id, pileups.pos), (0, pos));
            assert!(pileups.per_sample[0].is_some(), "A missing at pos {pos}");
            assert!(pileups.per_sample[1].is_some(), "B missing at pos {pos}");
        }
        assert!(merger.next().is_none());
    }

    #[test]
    fn two_readers_disjoint_positions() {
        let a = iter_from(vec![Ok(rec(0, 1)), Ok(rec(0, 3)), Ok(rec(0, 5))]);
        let b = iter_from(vec![Ok(rec(0, 2)), Ok(rec(0, 4)), Ok(rec(0, 6))]);
        let mut merger = PerPositionMerger::new(vec![a, b], names(2), Vec::new()).unwrap();
        let emissions: Vec<_> = (&mut merger).collect::<Result<_, _>>().unwrap();
        assert_eq!(emissions.len(), 6);
        for (idx, pileups) in emissions.iter().enumerate() {
            assert_eq!(pileups.chrom_id, 0);
            assert_eq!(pileups.pos as usize, idx + 1);
            // Even positions come from B, odd from A.
            let from_a = pileups.per_sample[0].is_some();
            let from_b = pileups.per_sample[1].is_some();
            assert!(from_a ^ from_b, "exactly one slot should be Some");
            assert_eq!(from_a, pileups.pos % 2 == 1);
        }
    }

    #[test]
    fn two_readers_partial_overlap() {
        // A: 1, 2, 4. B: 2, 3, 4.
        let a = iter_from(vec![Ok(rec(0, 1)), Ok(rec(0, 2)), Ok(rec(0, 4))]);
        let b = iter_from(vec![Ok(rec(0, 2)), Ok(rec(0, 3)), Ok(rec(0, 4))]);
        let mut merger = PerPositionMerger::new(vec![a, b], names(2), Vec::new()).unwrap();
        let emissions: Vec<_> = (&mut merger).collect::<Result<_, _>>().unwrap();
        let presence: Vec<(u32, bool, bool)> = emissions
            .iter()
            .map(|p| (p.pos, p.per_sample[0].is_some(), p.per_sample[1].is_some()))
            .collect();
        assert_eq!(
            presence,
            vec![
                (1, true, false),
                (2, true, true),
                (3, false, true),
                (4, true, true),
            ]
        );
    }

    #[test]
    fn multi_chromosome_chrom_is_major_sort_key() {
        // A and B both span chrom 0 then chrom 1; the merger must
        // exhaust chrom 0 before any chrom-1 emission, regardless of
        // pos values across chromosomes.
        let a = iter_from(vec![Ok(rec(0, 100)), Ok(rec(1, 5))]);
        let b = iter_from(vec![Ok(rec(0, 200)), Ok(rec(1, 50))]);
        let mut merger = PerPositionMerger::new(vec![a, b], names(2), Vec::new()).unwrap();
        let keys: Vec<(u32, u32)> = (&mut merger)
            .map(|p| {
                let p = p.unwrap();
                (p.chrom_id, p.pos)
            })
            .collect();
        assert_eq!(keys, vec![(0, 100), (0, 200), (1, 5), (1, 50)]);
    }

    #[test]
    fn three_readers_k2_tied_does_not_consume_unrelated_head() {
        // Readers 0 and 1 are tied at (0,1); reader 2 sits at (0,2)
        // and must NOT be consumed in the first emission.
        let a = iter_from(vec![Ok(rec(0, 1))]);
        let b = iter_from(vec![Ok(rec(0, 1))]);
        let c = iter_from(vec![Ok(rec(0, 2))]);
        let mut merger = PerPositionMerger::new(vec![a, b, c], names(3), Vec::new()).unwrap();

        let first = merger.next().unwrap().unwrap();
        assert_eq!((first.chrom_id, first.pos), (0, 1));
        assert!(first.per_sample[0].is_some());
        assert!(first.per_sample[1].is_some());
        assert!(
            first.per_sample[2].is_none(),
            "reader 2 was not at the min and must not contribute"
        );

        let second = merger.next().unwrap().unwrap();
        assert_eq!((second.chrom_id, second.pos), (0, 2));
        assert!(second.per_sample[0].is_none());
        assert!(second.per_sample[1].is_none());
        assert!(
            second.per_sample[2].is_some(),
            "reader 2's head should now be emitted"
        );
        assert!(merger.next().is_none());
    }

    #[test]
    fn reader_error_mid_stream_surfaces_and_latches_done() {
        // A streams cleanly through 1..=3; B sits at (0,5) and
        // errors on its second record. After three successful
        // emissions from A only, the merger reaches (0,5), advances
        // B, and surfaces B's refill error. Errors during refill
        // abort the in-progress emission by design — see
        // `implementation_plans/multi_way_per_position_iterator.md`.
        let a = iter_from(vec![Ok(rec(0, 1)), Ok(rec(0, 2)), Ok(rec(0, 3))]);
        let b = iter_from(vec![Ok(rec(0, 5)), Err(fake_err())]);
        let mut merger = PerPositionMerger::new(vec![a, b], names(2), Vec::new()).unwrap();

        for expected_pos in 1..=3 {
            let pileups = merger.next().unwrap().unwrap();
            assert_eq!((pileups.chrom_id, pileups.pos), (0, expected_pos));
            assert!(pileups.per_sample[0].is_some());
            assert!(pileups.per_sample[1].is_none());
        }
        let err = merger.next().unwrap().unwrap_err();
        let PerPositionMergerError::PerSampleReader {
            sample_idx,
            sample_name,
            ..
        } = err
        else {
            panic!("expected PerSampleReader error, got {err:?}");
        };
        assert_eq!(sample_idx, 1);
        assert_eq!(sample_name, "S1");

        // Once latched, subsequent calls return None.
        assert!(merger.next().is_none());
        assert!(merger.next().is_none());
    }

    #[test]
    fn reader_error_on_prefetch_aborts_construction() {
        let a = iter_from(vec![Ok(rec(0, 1))]);
        let b = iter_from(vec![Err(fake_err())]);
        let err = PerPositionMerger::new(vec![a, b], names(2), Vec::new()).unwrap_err();
        let PerPositionMergerError::PerSampleReader {
            sample_idx,
            sample_name,
            ..
        } = err
        else {
            panic!("expected PerSampleReader error, got {err:?}");
        };
        assert_eq!(sample_idx, 1);
        assert_eq!(sample_name, "S1");
    }

    #[test]
    fn out_of_order_record_is_detected() {
        // A single synthetic reader regresses: (0,5) then (0,3).
        // PspReader would have rejected this, but the merger's
        // own check defends against pathological mocks.
        let a = iter_from(vec![Ok(rec(0, 5)), Ok(rec(0, 3))]);
        let mut merger = PerPositionMerger::new(vec![a], names(1), Vec::new()).unwrap();
        let first = merger.next().unwrap().unwrap();
        assert_eq!((first.chrom_id, first.pos), (0, 5));
        let err = merger.next().unwrap().unwrap_err();
        let PerPositionMergerError::OutOfOrder {
            sample_idx,
            regressing_chrom_id,
            regressing_pos,
            last_emitted_chrom_id,
            last_emitted_pos,
            ..
        } = err
        else {
            panic!("expected OutOfOrder, got {err:?}");
        };
        assert_eq!(sample_idx, 0);
        assert_eq!((regressing_chrom_id, regressing_pos), (0, 3));
        assert_eq!((last_emitted_chrom_id, last_emitted_pos), (0, 5));
        assert!(merger.next().is_none());
    }

    #[test]
    fn emission_order_is_strictly_increasing() {
        let a = iter_from(vec![
            Ok(rec(0, 1)),
            Ok(rec(0, 4)),
            Ok(rec(1, 2)),
            Ok(rec(1, 9)),
        ]);
        let b = iter_from(vec![
            Ok(rec(0, 2)),
            Ok(rec(0, 4)),
            Ok(rec(1, 2)),
            Ok(rec(1, 7)),
        ]);
        let merger = PerPositionMerger::new(vec![a, b], names(2), Vec::new()).unwrap();
        let keys: Vec<(u32, u32)> = merger
            .map(|p| {
                let p = p.unwrap();
                (p.chrom_id, p.pos)
            })
            .collect();
        for window in keys.windows(2) {
            assert!(
                window[0] < window[1],
                "non-strictly-increasing emission: {:?} then {:?}",
                window[0],
                window[1],
            );
        }
    }

    #[test]
    fn sample_count_mismatch_rejected_at_construction() {
        let a = iter_from(vec![Ok(rec(0, 1))]);
        let b = iter_from(vec![Ok(rec(0, 1))]);
        let err = PerPositionMerger::new(vec![a, b], names(1), Vec::new()).unwrap_err();
        assert!(matches!(
            err,
            PerPositionMergerError::SampleCountMismatch {
                n_readers: 2,
                n_sample_names: 1,
            }
        ));
    }

    // ----- check_chromosome_agreement -----

    #[test]
    fn chromosome_agreement_empty_slice_is_ok() {
        let readers: Vec<PspReader<Cursor<Vec<u8>>>> = Vec::new();
        let chroms = check_chromosome_agreement(&readers).unwrap();
        assert!(chroms.is_empty());
    }

    #[test]
    fn chromosome_agreement_identical_lists_passes() {
        let r0 = psp_reader_with_header(writer_header(2));
        let r1 = psp_reader_with_header(writer_header(2));
        let chroms = check_chromosome_agreement(&[r0, r1]).unwrap();
        assert_eq!(chroms.len(), 2);
        assert_eq!(chroms[0].name, "chr1");
        assert_eq!(chroms[1].name, "chr2");
    }

    #[test]
    fn chromosome_agreement_differing_length_fails() {
        let h0 = writer_header(1);
        let mut h1 = writer_header(1);
        h1.chromosomes[0].length = 42;
        let r0 = psp_reader_with_header(h0);
        let r1 = psp_reader_with_header(h1);
        let err = check_chromosome_agreement(&[r0, r1]).unwrap_err();
        let PerPositionMergerError::ChromosomeMismatch {
            sample_idx,
            chrom_id,
            detail,
            ..
        } = err
        else {
            panic!("expected ChromosomeMismatch, got {err:?}");
        };
        assert_eq!(sample_idx, 1);
        assert_eq!(chrom_id, 0);
        assert!(detail.contains("length"), "detail = {detail:?}");
    }

    #[test]
    fn chromosome_agreement_differing_md5_fails() {
        let h0 = writer_header(1);
        let mut h1 = writer_header(1);
        h1.chromosomes[0].md5 = "1".repeat(32);
        let r0 = psp_reader_with_header(h0);
        let r1 = psp_reader_with_header(h1);
        let err = check_chromosome_agreement(&[r0, r1]).unwrap_err();
        assert!(matches!(
            err,
            PerPositionMergerError::ChromosomeMismatch {
                sample_idx: 1,
                chrom_id: 0,
                ..
            }
        ));
    }

    #[test]
    fn chromosome_agreement_differing_name_fails() {
        let h0 = writer_header(1);
        let mut h1 = writer_header(1);
        h1.chromosomes[0].name = "renamed".to_string();
        let r0 = psp_reader_with_header(h0);
        let r1 = psp_reader_with_header(h1);
        let err = check_chromosome_agreement(&[r0, r1]).unwrap_err();
        let PerPositionMergerError::ChromosomeMismatch { detail, .. } = err else {
            panic!("expected ChromosomeMismatch");
        };
        assert!(detail.contains("name"), "detail = {detail:?}");
    }

    #[test]
    fn chromosome_agreement_differing_count_fails() {
        let h0 = writer_header(2);
        let mut h1 = writer_header(2);
        h1.chromosomes.push(ChromosomeEntry {
            name: "chr3".to_string(),
            length: 1_000_000,
            md5: "0".repeat(32),
        });
        let r0 = psp_reader_with_header(h0);
        let r1 = psp_reader_with_header(h1);
        let err = check_chromosome_agreement(&[r0, r1]).unwrap_err();
        let PerPositionMergerError::ChromosomeCountMismatch {
            sample_idx,
            n_baseline,
            n_other,
            ..
        } = err
        else {
            panic!("expected ChromosomeCountMismatch, got {err:?}");
        };
        assert_eq!(sample_idx, 1);
        assert_eq!(n_baseline, 2);
        assert_eq!(n_other, 3);
    }
}
