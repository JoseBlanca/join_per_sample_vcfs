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

use std::iter::FusedIterator;

use crate::bam::alignment_input::MappedRead;
use crate::ng::ref_seq::RawRefSeq;
use crate::ng::types::{ContigId, GenomePosition, Position};

use super::IngestError;
use super::open_bam::RegionReads;

/// Argmin k-way merge over a sample's per-file region streams.
///
/// **Keys are held beside the heads, not read through them.** `keys[i]` mirrors
/// `heads[i]` and is refreshed only when that head is refilled, so each step
/// scans a small contiguous array of `GenomePosition`s in cache rather than
/// dereferencing k `MappedRead`s. That is spec §3.2's per-read budget expressed
/// in the layout rather than in a comment.
///
/// Deliberately **not** built on `Peekable`: `peek()` hands out a
/// `&Result<MappedRead, _>` and the emit path then moves out of the peek slot,
/// which is easy to write and quietly does the work twice.
pub struct MergedRegionReads<'a, R: RawRefSeq> {
    streams: Vec<RegionReads<'a, R>>,
    /// `None` = that stream is exhausted.
    heads: Vec<Option<MappedRead>>,
    /// `keys[i]` is `heads[i]`'s key; `None` in lockstep.
    keys: Vec<Option<GenomePosition>>,
    /// Cleared once every head has been primed.
    filled: bool,
    /// An error met while refilling *after* a good read was already in hand.
    ///
    /// The read is emitted first and this is yielded on the next call. Without
    /// it the read would be dropped on the floor — a silently lost read, which
    /// is the one outcome this whole layer exists to prevent.
    pending_error: Option<IngestError>,
    /// Set on clean exhaustion or after an error is yielded, so the iterator
    /// fuses.
    done: bool,
}

impl<'a, R: RawRefSeq> MergedRegionReads<'a, R> {
    // Built by `SampleReads::reads_in_region` from step E3; until then only
    // this module's own tests construct one. `expect`, so the build fails
    // naming this line once E3 makes it live.
    #[cfg_attr(
        not(test),
        expect(dead_code, reason = "SampleReads::reads_in_region (E3) builds this")
    )]
    pub(crate) fn new(streams: Vec<RegionReads<'a, R>>) -> Self {
        let k = streams.len();
        Self {
            streams,
            heads: (0..k).map(|_| None).collect(),
            keys: vec![None; k],
            filled: false,
            pending_error: None,
            done: false,
        }
    }

    /// Pull one read into slot `i`, recomputing just that key.
    ///
    /// Returns the wrapped error if the stream yielded one; the caller fuses.
    fn refill(&mut self, i: usize) -> Result<(), IngestError> {
        match self.streams[i].next() {
            None => {
                self.heads[i] = None;
                self.keys[i] = None;
                Ok(())
            }
            Some(Ok(read)) => {
                self.keys[i] = Some(key_of(&read));
                self.heads[i] = Some(read);
                Ok(())
            }
            Some(Err(source)) => {
                self.heads[i] = None;
                self.keys[i] = None;
                Err(IngestError::File {
                    source_file_index: i,
                    source,
                })
            }
        }
    }

    /// The lowest key, **ties to the lowest file index**.
    ///
    /// The tie-break is not incidental: because a sample's files usually come
    /// from different experiments and so cover the same coordinate range, ties
    /// are routine rather than rare, and this is what makes the output order
    /// reproducible run to run.
    fn argmin(&self) -> Option<usize> {
        let mut best: Option<usize> = None;
        for (i, key) in self.keys.iter().enumerate() {
            let Some(key) = key else { continue };
            match best {
                // Strictly less, so an equal key leaves the earlier index in
                // place — that *is* the tie-break.
                Some(b) if key >= &self.keys[b].expect("best always has a key") => {}
                _ => best = Some(i),
            }
        }
        best
    }

    /// The same-file-twice check, run **only on a tie**.
    ///
    /// Two reads at different positions cannot be the same read, so this has
    /// nothing to do unless the argmin already found ≥2 heads sharing a
    /// position — which it learned as a by-product of the comparison it had to
    /// make anyway. Within a tie, `flag` (a `u16`) is compared before `qname`
    /// (a short string), so the string compare happens only for reads that
    /// already agree on position *and* flags.
    ///
    /// This is **not deduplication**. Across files there is no such thing as a
    /// legitimate duplicate — reads from different experiments are different
    /// reads — so a match means the caller passed the same file twice, which is
    /// an input-sanity error and never a silent drop.
    fn same_read_across_files(&self, winner: usize) -> Option<IngestError> {
        let winning_key = self.keys[winner]?;
        let winning_read = self.heads[winner].as_ref()?;

        for (other, key) in self.keys.iter().enumerate().skip(winner + 1) {
            if *key != Some(winning_key) {
                continue;
            }
            let Some(other_read) = self.heads[other].as_ref() else {
                continue;
            };
            if other_read.flag == winning_read.flag && other_read.qname == winning_read.qname {
                return Some(IngestError::DuplicateReadAcrossFiles {
                    qname: winning_read.qname.clone(),
                    key: winning_key,
                    files: (winner, other),
                });
            }
        }
        None
    }

    fn fail(&mut self, error: IngestError) -> Option<Result<MappedRead, IngestError>> {
        self.done = true;
        Some(Err(error))
    }
}

/// A read's position. Sound as a *cross-file* key only because every file's
/// open gate proved its `ref_id`s are the reference's `ContigId`s — without
/// that, contig indices from different files would not be comparable at all.
fn key_of(read: &MappedRead) -> GenomePosition {
    GenomePosition {
        // PANIC-FREE: `ref_id` comes from a 32-bit field in both BAM and CRAM.
        contig: ContigId(u32::try_from(read.ref_id).expect("ref_id fits u32")),
        position: Position(read.pos),
    }
}

impl<R: RawRefSeq> Iterator for MergedRegionReads<'_, R> {
    type Item = Result<MappedRead, IngestError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        // A refill error held over from the previous call, after that call's
        // read was safely emitted.
        if let Some(error) = self.pending_error.take() {
            return self.fail(error);
        }

        if !self.filled {
            self.filled = true;
            for i in 0..self.streams.len() {
                if let Err(error) = self.refill(i) {
                    return self.fail(error);
                }
            }
        }

        let Some(winner) = self.argmin() else {
            self.done = true;
            return None;
        };

        if let Some(error) = self.same_read_across_files(winner) {
            return self.fail(error);
        }

        // One move per read: `take` hands the read out and empties the slot,
        // which is then refilled from its own stream. No clone — a `MappedRead`
        // owns its sequence, qualities and CIGAR, and cloning one here would be
        // a defect rather than a slow path.
        let read = self.heads[winner].take().expect("the winner has a head");
        self.keys[winner] = None;

        // The read in hand is good regardless of what the refill finds, so it
        // is emitted now and any error is held over. Returning the error here
        // instead would discard a read the caller had already earned.
        if let Err(error) = self.refill(winner) {
            self.pending_error = Some(error);
        }

        Some(Ok(read))
    }
}

impl<R: RawRefSeq> FusedIterator for MergedRegionReads<'_, R> {}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use noodles_sam::alignment::RecordBuf;
    use tempfile::TempDir;

    use super::*;
    use crate::ng::read::filtering::ReadFilterConfig;
    use crate::ng::read::input::open_bam::AlignmentFile;
    use crate::ng::read::input::test_fixtures::{
        FIXTURE_CONTIGS, bam_header, fixture_reference, indexed_bam, matching_contigs,
        read_named_with_length,
    };
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::ng::types::GenomeRegion;

    fn reference_bases() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(
            FIXTURE_CONTIGS
                .iter()
                .map(|(_, length)| vec![b'A'; *length])
                .collect(),
        )
    }

    fn whole_first_contig() -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(1),
            end: Position(FIXTURE_CONTIGS[0].1 as u64),
        }
    }

    fn bam_with(records: &[RecordBuf]) -> (TempDir, PathBuf) {
        indexed_bam(&bam_header(&matching_contigs()), records)
    }

    /// Open `paths` as one sample's files and merge their streams over the
    /// whole first contig, returning `(qname, source_file_index)` per read.
    fn merged(paths: &[PathBuf]) -> Result<Vec<(String, usize)>, IngestError> {
        let (_reference_dir, reference) = fixture_reference(false);
        let files: Vec<AlignmentFile> = paths
            .iter()
            .enumerate()
            .map(|(i, path)| {
                AlignmentFile::open_as(path, &reference, ReadFilterConfig::default(), false, i)
                    .expect("opens")
            })
            .collect();

        let streams: Vec<_> = files
            .iter()
            .map(|file| {
                file.reads_in_region(whole_first_contig(), reference_bases())
                    .expect("query")
            })
            .collect();

        MergedRegionReads::new(streams)
            .map(|item| {
                item.map(|read| {
                    (
                        String::from_utf8_lossy(&read.qname).into_owned(),
                        read.source_file_index,
                    )
                })
            })
            .collect()
    }

    /// **T6 — two files interleave in coordinate order, each read tagged with
    /// the file it came from.**
    ///
    /// The tag is not decoration: because the usual reason for several files is
    /// several *experiments*, it is the sample's batch label, and a per-batch
    /// error model keys on it.
    #[test]
    fn t6_two_files_interleave_in_coordinate_order_with_their_file_tags() {
        let (_first_dir, first) = bam_with(&[
            read_named_with_length("a1", 0, 1, 30),
            read_named_with_length("a2", 0, 40, 30),
        ]);
        let (_second_dir, second) = bam_with(&[
            read_named_with_length("b1", 0, 20, 30),
            read_named_with_length("b2", 0, 60, 30),
        ]);

        let reads = merged(&[first, second]).expect("both stream cleanly");

        assert_eq!(
            reads,
            vec![
                ("a1".to_string(), 0),
                ("b1".to_string(), 1),
                ("a2".to_string(), 0),
                ("b2".to_string(), 1),
            ],
            "coordinate order across files, each read keeping its own tag"
        );
    }

    /// **T6, the tie-break.** Reads at the same position break to the **lower
    /// file index**. This is load-bearing rather than cosmetic: per-experiment
    /// files cover the same coordinate range, so ties are routine, and without
    /// a fixed rule the output order would vary run to run.
    #[test]
    fn t6_reads_at_one_position_break_to_the_lower_file_index() {
        let (_first_dir, first) = bam_with(&[read_named_with_length("from_file_0", 0, 10, 30)]);
        let (_second_dir, second) = bam_with(&[read_named_with_length("from_file_1", 0, 10, 30)]);

        let reads = merged(&[first.clone(), second.clone()]).expect("streams");
        assert_eq!(
            reads,
            vec![
                ("from_file_0".to_string(), 0),
                ("from_file_1".to_string(), 1)
            ],
            "the lower file index wins the tie"
        );

        // And swapping the inputs swaps the order — proving the rule is the
        // *index*, not something incidental about the reads themselves.
        let swapped = merged(&[second, first]).expect("streams");
        assert_eq!(
            swapped,
            vec![
                ("from_file_1".to_string(), 0),
                ("from_file_0".to_string(), 1)
            ],
        );
    }

    /// **T6, determinism.** The same inputs must give byte-identical output
    /// every run — which is what the tie-break buys.
    #[test]
    fn t6_the_merged_order_is_identical_across_runs() {
        let mut piled = Vec::new();
        for position in [5u64, 5, 5, 20, 20, 40] {
            piled.push(read_named_with_length(
                &format!("a{position}"),
                0,
                position as usize,
                30,
            ));
        }
        let (_first_dir, first) = bam_with(&piled);
        let (_second_dir, second) = bam_with(&[
            read_named_with_length("b5", 0, 5, 30),
            read_named_with_length("b20", 0, 20, 30),
        ]);

        let once = merged(&[first.clone(), second.clone()]).expect("streams");
        for _ in 0..5 {
            assert_eq!(
                merged(&[first.clone(), second.clone()]).expect("streams"),
                once,
                "the merge must be reproducible, not merely correct"
            );
        }
    }

    /// **T7 — the same file passed twice is an error, at the *first*
    /// collision.**
    ///
    /// Not a deduplication: across files there is no such thing as a legitimate
    /// duplicate, because reads from different experiments are different reads.
    /// A match therefore means duplicated *input*, and since every read
    /// collides, the very first overlapping position trips it — which is why
    /// the check can afford to be this cheap.
    #[test]
    fn t7_the_same_file_twice_errors_at_the_first_collision() {
        let (_dir, path) = bam_with(&[
            read_named_with_length("r1", 0, 10, 30),
            read_named_with_length("r2", 0, 50, 30),
        ]);

        let error = merged(&[path.clone(), path]).expect_err("the same file twice is not a sample");

        match &error {
            IngestError::DuplicateReadAcrossFiles { qname, key, files } => {
                assert_eq!(qname, b"r1", "the *first* colliding read, not the last");
                assert_eq!(key.position.get(), 10);
                assert_eq!(*files, (0, 1));
            }
            other => panic!("expected DuplicateReadAcrossFiles, got {other:?}"),
        }
        assert!(
            error.to_string().contains("read 'r1'"),
            "the message names the read: {error}"
        );
    }

    /// **T7's other half — distinct reads at one locus are not duplicates.**
    ///
    /// Two reads sharing `(ref_id, pos)` and differing only in `qname` is the
    /// ordinary case at any covered position, and it is exactly what the
    /// cheap-first comparison order must not get wrong: `flag` matches, so the
    /// check falls through to `qname`, which does not.
    #[test]
    fn t7_reads_sharing_a_position_but_not_a_name_both_survive() {
        let (_first_dir, first) = bam_with(&[read_named_with_length("alpha", 0, 10, 30)]);
        let (_second_dir, second) = bam_with(&[read_named_with_length("beta", 0, 10, 30)]);

        let reads = merged(&[first, second]).expect("different reads at one locus are ordinary");
        assert_eq!(
            reads,
            vec![("alpha".to_string(), 0), ("beta".to_string(), 1)],
            "same position, different names — both kept"
        );
    }

    /// A file that runs out early must not end the merge: the others keep
    /// streaming.
    #[test]
    fn an_exhausted_file_does_not_end_the_merge() {
        let (_first_dir, first) = bam_with(&[read_named_with_length("early", 0, 1, 30)]);
        let (_second_dir, second) = bam_with(&[
            read_named_with_length("b1", 0, 20, 30),
            read_named_with_length("b2", 0, 40, 30),
            read_named_with_length("b3", 0, 60, 30),
        ]);

        let reads = merged(&[first, second]).expect("streams");
        assert_eq!(
            reads,
            vec![
                ("early".to_string(), 0),
                ("b1".to_string(), 1),
                ("b2".to_string(), 1),
                ("b3".to_string(), 1),
            ]
        );
    }

    /// Three files, to prove the argmin is not accidentally a two-way compare.
    #[test]
    fn three_files_merge_in_coordinate_order() {
        let (_a_dir, a) = bam_with(&[read_named_with_length("a", 0, 1, 30)]);
        let (_b_dir, b) = bam_with(&[read_named_with_length("b", 0, 20, 30)]);
        let (_c_dir, c) = bam_with(&[read_named_with_length("c", 0, 40, 30)]);

        let reads = merged(&[c, a, b]).expect("streams");
        assert_eq!(
            reads
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            vec!["a", "b", "c"],
            "coordinate order, regardless of the order the files were given"
        );
    }
}
