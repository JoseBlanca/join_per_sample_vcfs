//! Multi-file coordinate merge over the segment reader.
//!
//! One sample's pileup may have several input alignment files (e.g.
//! per-lane BAMs). [`AlignmentFile`](crate::bam::segment_reader::AlignmentFile)
//! serves the reads for one file + one segment; [`SegmentMergedReads`]
//! k-way-merges those per-file streams into a single coordinate-sorted
//! `MappedRead` stream — the shape the pileup walker consumes — preserving
//! the guarantees the old `AlignmentMergedReader::query` provided:
//!
//! - **coordinate order** `(ref_id, pos)` across all inputs (argmin merge);
//! - **cross-file duplicate detection** at a locus (errors, not silent drop);
//! - **within-file out-of-order detection** + **fuse-on-error**;
//! - **`FilterCounts`** summed from the per-file fetchers (each already
//!   tallies its own cheap-filter drops).
//!
//! What it does *not* do (already handled upstream in the fetcher): the
//! cheap flag/MAPQ/length filter, the `RecordBuf → MappedRead` conversion,
//! and the per-record decode. This is purely the merge.
//!
//! This is increment #5 of the segment-read-fetcher work (the SNP
//! `--regions` retrofit); see
//! `doc/devel/architecture/segment_reader_snp_regions_retrofit.md`. Until
//! `run_pileup` is flipped onto it, the only consumers are the tests.

// No production consumer until the driver flip; the surface is otherwise
// reachable only from `#[cfg(test)]`. Suppress dead-code in the non-test
// build only, so the test build still flags genuinely-dead helpers.
#![cfg_attr(not(test), allow(dead_code))]

use std::path::PathBuf;

use super::alignment_input::{FilterCounts, MappedRead};
use super::errors::AlignmentInputError;
use super::segment_reader::{AlignmentFile, MappedReadsInSegment};

/// A genome location: `ContigList` index + 1-based position. The merge
/// key and the per-file order-regression key. Mirrors the alias the
/// merged reader uses internally.
type Locus = (usize, u64);

/// The four fields whose equality defines "same read" for cross-file
/// duplicate detection — qname, SAM flags, contig-list index, position.
/// Same definition the merged reader uses (`per_sample_pileup.md`
/// §"Duplicate-read detection across CRAMs").
#[derive(Debug, Clone, PartialEq, Eq)]
struct ReadFingerprint {
    qname: Vec<u8>,
    flag: u16,
    ref_id: usize,
    pos: u64,
}

/// A `ReadFingerprint` paired with the input it came from, so a duplicate
/// error can name *both* the previous acceptance and the colliding one.
#[derive(Debug, Clone)]
struct ReadFingerprintWithSourceFile {
    key: ReadFingerprint,
    source_file_index: usize,
}

/// One per-file segment stream with a one-slot peek buffer. Unlike
/// [`std::iter::Peekable`], this exposes the inner iterator so the merge
/// can read its [`FilterCounts`] after draining (std `Peekable` hides it).
struct PeekableSegmentStream<'a> {
    inner: MappedReadsInSegment<'a>,
    /// Buffered head: `None` = not yet peeked / consumed; `Some(item)` =
    /// the next item, already pulled from `inner`.
    head: Option<Result<MappedRead, AlignmentInputError>>,
}

impl<'a> PeekableSegmentStream<'a> {
    fn new(inner: MappedReadsInSegment<'a>) -> Self {
        Self { inner, head: None }
    }

    fn peek(&mut self) -> Option<&Result<MappedRead, AlignmentInputError>> {
        if self.head.is_none() {
            self.head = self.inner.next();
        }
        self.head.as_ref()
    }

    fn next(&mut self) -> Option<Result<MappedRead, AlignmentInputError>> {
        match self.head.take() {
            Some(item) => Some(item),
            None => self.inner.next(),
        }
    }

    /// The inner fetcher's cheap-filter drop tally so far. Final once the
    /// stream is drained (the peek buffer has been refilled past the end).
    fn filter_counts(&self) -> &FilterCounts {
        self.inner.filter_counts()
    }
}

/// The coordinate-sorted, cross-file-merged read stream for one segment.
/// Borrows the pooled [`AlignmentFile`]s; on drop, each inner
/// [`MappedReadsInSegment`] returns its reader to its file's pool.
pub(crate) struct SegmentMergedReads<'a> {
    /// One peekable stream per input file, in input order.
    streams: Vec<PeekableSegmentStream<'a>>,
    /// Per-input file path, parallel to `streams`, for error messages.
    paths: Vec<PathBuf>,
    /// Previous accepted `Locus` per input — within-file order guard.
    per_file_prev_locus: Vec<Option<Locus>>,
    /// The locus the merge is currently at; the fingerprints below all
    /// share it. Cleared when the merge advances past it.
    current_locus: Option<Locus>,
    /// Reads already accepted at `current_locus` — the cross-file dedup
    /// buffer.
    current_locus_fingerprints: Vec<ReadFingerprintWithSourceFile>,
    /// Set on the first error; every later poll returns `None`.
    fused: bool,
}

impl<'a> SegmentMergedReads<'a> {
    /// Open one segment stream per input file and merge them. The reads
    /// overlap `[start, end]` (1-based inclusive) on `chrom`, are
    /// cheap-filtered by each file's `SegmentReadFilter`, and are emitted
    /// in coordinate order across all inputs.
    ///
    /// # Errors
    ///
    /// Any error opening a per-file segment stream (bad segment, unknown
    /// contig, open/seek/index failure) surfaces here before iteration
    /// begins; per-record decode errors and merge violations (out-of-order,
    /// cross-file duplicate) surface as `Some(Err(_))` during iteration.
    pub(crate) fn new(
        files: &'a [AlignmentFile],
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<Self, AlignmentInputError> {
        let mut streams = Vec::with_capacity(files.len());
        let mut paths = Vec::with_capacity(files.len());
        for file in files {
            // On an error here, streams opened so far drop and return
            // their borrowed readers to their pools.
            streams.push(PeekableSegmentStream::new(
                file.get_reads_from_segment(chrom, start, end)?,
            ));
            paths.push(file.path().to_path_buf());
        }
        let n = files.len();
        Ok(Self {
            streams,
            paths,
            per_file_prev_locus: vec![None; n],
            current_locus: None,
            current_locus_fingerprints: Vec::new(),
            fused: false,
        })
    }

    /// The cheap-filter drop tally summed across every input. Read after
    /// draining for the per-segment totals. No merge-level drops exist —
    /// the merge's own rejections (out-of-order, duplicate) are *errors*,
    /// not filtered reads — so this is exactly the sum of the fetchers'.
    pub(crate) fn filter_counts(&self) -> FilterCounts {
        let mut total = FilterCounts::default();
        for stream in &self.streams {
            total.merge(stream.filter_counts());
        }
        total
    }

    fn fail(&mut self, e: AlignmentInputError) -> Option<Result<MappedRead, AlignmentInputError>> {
        self.fused = true;
        Some(Err(e))
    }

    /// Clear the at-locus dedup buffer when the merge advances to a new
    /// locus.
    fn advance_current_locus_if_needed(&mut self, ref_id: usize, pos: u64) {
        match self.current_locus {
            Some(existing) if existing != (ref_id, pos) => {
                self.current_locus_fingerprints.clear();
                self.current_locus = Some((ref_id, pos));
            }
            None => self.current_locus = Some((ref_id, pos)),
            _ => {}
        }
    }

    /// Index of the stream whose head has the smallest `(ref_id, pos)`,
    /// or `None` if every stream is exhausted. Ties break to the lowest
    /// input index (the `<` keeps the first seen), matching the merged
    /// reader's `argmin_head`.
    fn argmin_head(&mut self) -> Option<usize> {
        let mut best: Option<(usize, Locus)> = None;
        for idx in 0..self.streams.len() {
            let Some(Ok(read)) = self.streams[idx].peek() else {
                continue;
            };
            let candidate = (read.ref_id, read.pos);
            best = match best {
                None => Some((idx, candidate)),
                Some((_, current)) if candidate < current => Some((idx, candidate)),
                Some(existing) => Some(existing),
            };
        }
        best.map(|(idx, _)| idx)
    }
}

impl Iterator for SegmentMergedReads<'_> {
    type Item = Result<MappedRead, AlignmentInputError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.fused {
            return None;
        }
        loop {
            // 1. Surface a peeked error from the lowest-index stream that
            //    has one, and fuse. (Deterministic ordering.)
            for idx in 0..self.streams.len() {
                if matches!(self.streams[idx].peek(), Some(Err(_))) {
                    let err = self.streams[idx]
                        .next()
                        .expect("peeked Some")
                        .expect_err("peeked Err");
                    return self.fail(err);
                }
            }

            // 2. Pick the smallest surviving head; None → all exhausted.
            let chosen = self.argmin_head()?;

            // 3. Read the chosen head's keys (it is `Ok` per argmin).
            let (ref_id, pos, flag, qname) = match self.streams[chosen].peek() {
                Some(Ok(read)) => (read.ref_id, read.pos, read.flag, read.qname.clone()),
                _ => continue,
            };

            // 4. Within-file order check against this file's previous accept.
            //    Accepted divergence from the old reader: a read that is both
            //    too-short and a coordinate regression was dropped (and
            //    counted) by the fetcher before reaching here, so it no longer
            //    triggers this error. Corrupt-input edge only; well-formed
            //    sorted data is unaffected (see the #5 architecture doc §2).
            if let Some((prev_ref, prev_pos)) = self.per_file_prev_locus[chosen]
                && (ref_id, pos) < (prev_ref, prev_pos)
            {
                return self.fail(AlignmentInputError::OutOfOrderRead {
                    path: self.paths[chosen].clone(),
                    qname: String::from_utf8_lossy(&qname).into_owned(),
                    prev_ref_id: prev_ref,
                    prev_pos,
                    this_ref_id: ref_id,
                    this_pos: pos,
                });
            }

            // 5. Advance the locus, clearing the dedup buffer on a move.
            self.advance_current_locus_if_needed(ref_id, pos);

            // 6. Cross-file duplicate check at the current locus.
            let new_key = ReadFingerprint {
                qname: qname.clone(),
                flag,
                ref_id,
                pos,
            };
            if let Some(other) = self
                .current_locus_fingerprints
                .iter()
                .find(|entry| entry.key == new_key)
            {
                return self.fail(AlignmentInputError::DuplicateReadAcrossFiles {
                    qname: String::from_utf8_lossy(&new_key.qname).into_owned(),
                    path_a: self.paths[other.source_file_index].clone(),
                    path_b: self.paths[chosen].clone(),
                    ref_id,
                    pos,
                });
            }

            // 7. Consume the chosen head and emit.
            let mapped = match self.streams[chosen].next() {
                Some(Ok(read)) => read,
                Some(Err(e)) => return self.fail(e),
                None => continue,
            };
            self.current_locus_fingerprints
                .push(ReadFingerprintWithSourceFile {
                    key: new_key,
                    source_file_index: chosen,
                });
            self.per_file_prev_locus[chosen] = Some((mapped.ref_id, mapped.pos));
            return Some(Ok(mapped));
        }
    }
}

impl std::iter::FusedIterator for SegmentMergedReads<'_> {}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::num::NonZero;
    use std::path::Path;
    use std::sync::Arc;

    use noodles_bam as bam;
    use noodles_core::Position;
    use noodles_csi::binning_index::Indexer;
    use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles_csi::binning_index::index::reference_sequence::index::BinnedIndex;
    use noodles_sam as sam;
    use noodles_sam::alignment::Record as _;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::MappingQuality;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::ReferenceSequence;
    use tempfile::TempDir;

    use super::*;
    use crate::bam::index_preflight::AlignmentIndex;
    use crate::bam::segment_reader::SegmentReadFilter;

    const CONTIG_LEN: usize = 200;

    fn permissive() -> SegmentReadFilter {
        SegmentReadFilter {
            min_mapq: None,
            min_read_length: None,
            drop_qc_fail: false,
            drop_duplicate: false,
        }
    }

    fn header() -> sam::Header {
        sam::Header::builder()
            .set_header(Default::default())
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZero::new(CONTIG_LEN).unwrap()),
            )
            .build()
    }

    fn record(qname: &str, start: usize, len: usize, mapq: u8) -> RecordBuf {
        RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_reference_sequence_id(0)
            .set_flags(Flags::default())
            .set_mapping_quality(MappingQuality::new(mapq).expect("mapq"))
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar([Op::new(Kind::Match, len)].into_iter().collect())
            .set_sequence(Sequence::from(vec![b'A'; len]))
            .set_quality_scores(QualityScores::from(vec![30u8; len]))
            .build()
    }

    fn build_csi(bam_path: &Path) -> noodles_csi::Index {
        let mut reader = bam::io::reader::Builder
            .build_from_path(bam_path)
            .expect("open bam");
        let parsed = reader.read_header().expect("header");
        let mut indexer: Indexer<BinnedIndex> = Indexer::default();
        let mut chunk_start = reader.get_ref().virtual_position();
        let mut rec = bam::Record::default();
        while reader.read_record(&mut rec).expect("read") != 0 {
            let chunk_end = reader.get_ref().virtual_position();
            let ctx = match (
                rec.reference_sequence_id().transpose().expect("ref"),
                rec.alignment_start().transpose().expect("start"),
                rec.alignment_end().transpose().expect("end"),
            ) {
                (Some(id), Some(s), Some(e)) => Some((id, s, e, !rec.flags().is_unmapped())),
                _ => None,
            };
            indexer
                .add_record(ctx, Chunk::new(chunk_start, chunk_end))
                .expect("add");
            chunk_start = chunk_end;
        }
        indexer.build(parsed.reference_sequences().len())
    }

    /// Build one BAM-backed `AlignmentFile` over `records`, stamped with
    /// `source_file_index`. Keeps the `TempDir` alive via the handle.
    fn bam_file(
        records: &[RecordBuf],
        filter: SegmentReadFilter,
        source_file_index: usize,
    ) -> (TempDir, AlignmentFile) {
        let hdr = header();
        let dir = TempDir::new().expect("tempdir");
        let path = dir.path().join(format!("in{source_file_index}.bam"));
        let mut writer = bam::io::Writer::new(File::create(&path).expect("create"));
        writer.write_header(&hdr).expect("write header");
        for r in records {
            writer.write_alignment_record(&hdr, r).expect("write");
        }
        writer.try_finish().expect("finish");
        let csi = build_csi(&path);
        let file = AlignmentFile::from_input(
            path,
            Arc::new(hdr),
            AlignmentIndex::BamCsi(Arc::new(csi)),
            None,
            filter,
            source_file_index,
        )
        .expect("from_input");
        (dir, file)
    }

    fn drain(merged: SegmentMergedReads<'_>) -> Vec<(u64, String, usize)> {
        merged
            .map(|r| {
                let read = r.expect("merge ok");
                (
                    read.pos,
                    String::from_utf8_lossy(&read.qname).into_owned(),
                    read.source_file_index,
                )
            })
            .collect()
    }

    #[test]
    fn merges_two_files_in_coordinate_order_with_source_tags() {
        let (_a, file_a) = bam_file(
            &[record("a10", 10, 4, 60), record("a30", 30, 4, 60)],
            permissive(),
            0,
        );
        let (_b, file_b) = bam_file(
            &[record("b20", 20, 4, 60), record("b40", 40, 4, 60)],
            permissive(),
            1,
        );
        let files = [file_a, file_b];

        let merged = SegmentMergedReads::new(&files, "chr1", 1, 200).expect("merge");
        let got = drain(merged);
        assert_eq!(
            got,
            vec![
                (10, "a10".into(), 0),
                (20, "b20".into(), 1),
                (30, "a30".into(), 0),
                (40, "b40".into(), 1),
            ]
        );
    }

    #[test]
    fn single_file_passes_through_sorted() {
        let (_a, file_a) = bam_file(
            &[record("r1", 5, 4, 60), record("r2", 25, 4, 60)],
            permissive(),
            0,
        );
        let files = [file_a];
        let merged = SegmentMergedReads::new(&files, "chr1", 1, 200).expect("merge");
        assert_eq!(
            drain(merged),
            vec![(5, "r1".into(), 0), (25, "r2".into(), 0)]
        );
    }

    #[test]
    fn cross_file_duplicate_at_a_locus_errors() {
        // Same read (qname/flag/ref/pos) in both files at chr1:10.
        let (_a, file_a) = bam_file(&[record("dup", 10, 4, 60)], permissive(), 0);
        let (_b, file_b) = bam_file(&[record("dup", 10, 4, 60)], permissive(), 1);
        let files = [file_a, file_b];

        let mut merged = SegmentMergedReads::new(&files, "chr1", 1, 200).expect("merge");
        // First copy emitted.
        assert!(matches!(merged.next(), Some(Ok(_))));
        // Second copy at the same locus → duplicate error.
        assert!(matches!(
            merged.next(),
            Some(Err(AlignmentInputError::DuplicateReadAcrossFiles {
                ref_id: 0,
                pos: 10,
                ..
            }))
        ));
        // Fused.
        assert!(merged.next().is_none());
    }

    #[test]
    fn distinct_reads_at_the_same_locus_are_not_duplicates() {
        // Same locus, different qname → both kept (lower input index first).
        let (_a, file_a) = bam_file(&[record("ra", 10, 4, 60)], permissive(), 0);
        let (_b, file_b) = bam_file(&[record("rb", 10, 4, 60)], permissive(), 1);
        let files = [file_a, file_b];
        let merged = SegmentMergedReads::new(&files, "chr1", 1, 200).expect("merge");
        assert_eq!(
            drain(merged),
            vec![(10, "ra".into(), 0), (10, "rb".into(), 1)]
        );
    }

    #[test]
    fn filter_counts_are_summed_across_inputs() {
        // File A drops a low-MAPQ read; file B drops a too-short read.
        let cfg = SegmentReadFilter::default(); // min_mapq=20, min_read_length=30
        let (_a, file_a) = bam_file(
            &[record("keepA", 10, 40, 60), record("lowmq", 12, 40, 5)],
            cfg,
            0,
        );
        let (_b, file_b) = bam_file(
            &[record("keepB", 20, 40, 60), record("short", 22, 10, 60)],
            cfg,
            1,
        );
        let files = [file_a, file_b];

        let mut merged = SegmentMergedReads::new(&files, "chr1", 1, 200).expect("merge");
        let kept: Vec<String> = merged
            .by_ref()
            .map(|r| String::from_utf8_lossy(&r.expect("ok").qname).into_owned())
            .collect();
        let counts = merged.filter_counts();

        assert_eq!(kept, vec!["keepA".to_string(), "keepB".to_string()]);
        assert_eq!(counts.low_mapq, 1);
        assert_eq!(counts.too_short, 1);
    }

    #[test]
    fn empty_segment_yields_nothing() {
        let (_a, file_a) = bam_file(&[record("a", 10, 4, 60)], permissive(), 0);
        let (_b, file_b) = bam_file(&[record("b", 12, 4, 60)], permissive(), 1);
        let files = [file_a, file_b];
        let merged = SegmentMergedReads::new(&files, "chr1", 100, 150).expect("merge");
        assert!(drain(merged).is_empty());
    }

    /// The gate for the #5 driver flip: the new merge must produce the
    /// **same** surviving reads (same order) and the **same** `FilterCounts`
    /// as the old `AlignmentMergedReader::query` over the same multi-file,
    /// sub-region input. (The accepted divergence — a read that is both
    /// too-short and out-of-order — is a corrupt-input edge not present in
    /// well-formed sorted fixtures, so it does not surface here.)
    #[test]
    fn byte_identical_to_alignment_merged_reader_query() {
        use std::path::PathBuf;

        use crate::bam::alignment_input::{
            AlignmentMergedReader, AlignmentMergedReaderConfig, ContigInterval,
            build_fasta_repository,
        };
        use crate::fasta::{ContigEntry, ContigList};
        use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};

        // Shared FASTA (chr1, 200) for the repository the old query requires
        // (BAM decode never consults it, but the arg is mandatory).
        let (_fa_dir, fasta_path) = build_fasta(&[ContigSpec {
            name: "chr1".into(),
            length: CONTIG_LEN as u64,
        }])
        .expect("fasta");
        let repository = build_fasta_repository(&fasta_path).expect("repo");

        // Two inputs for one sample: interleaved keepers plus a low-MAPQ drop
        // in A and a too-short drop in B, so counts parity is exercised too.
        let good = 34usize; // > DEFAULT_MIN_READ_LENGTH (30)
        let hdr = Arc::new(header());
        let dir = TempDir::new().expect("dir");
        let write = |name: &str, recs: &[RecordBuf]| -> PathBuf {
            let p = dir.path().join(name);
            let mut w = bam::io::Writer::new(File::create(&p).expect("create"));
            w.write_header(&hdr).expect("hdr");
            for r in recs {
                w.write_alignment_record(&hdr, r).expect("rec");
            }
            w.try_finish().expect("finish");
            p
        };
        let path_a = write(
            "a.bam",
            &[
                record("a_keep1", 5, good, 60),
                record("a_lowmq", 8, good, 5), // dropped: MAPQ 5 < 20
                record("a_keep2", 30, good, 60),
            ],
        );
        let path_b = write(
            "b.bam",
            &[
                record("b_keep1", 12, good, 60),
                record("b_short", 15, 10, 60), // dropped: SEQ 10 < 30
                record("b_keep2", 40, good, 60),
            ],
        );

        let idx_a = AlignmentIndex::BamCsi(Arc::new(build_csi(&path_a)));
        let idx_b = AlignmentIndex::BamCsi(Arc::new(build_csi(&path_b)));
        let region = ContigInterval { start: 1, end: 50 };

        // --- old path: AlignmentMergedReader::query ---
        let old_cfg = AlignmentMergedReaderConfig {
            min_mapq: Some(20),
            min_read_length: Some(30),
            drop_qc_fail: true,
            drop_duplicate: true,
            max_read_mismatch_fraction: None, // the merged reader never applies F1 itself
            mismatch_bq_floor: 0,
        };
        let mut old = AlignmentMergedReader::query(
            &[path_a.clone(), path_b.clone()],
            &repository,
            ContigList {
                entries: vec![ContigEntry {
                    name: "chr1".into(),
                    length: CONTIG_LEN as u64,
                    md5: None,
                }],
            },
            "sample0".into(),
            &[hdr.clone(), hdr.clone()],
            &[idx_a.clone(), idx_b.clone()],
            "chr1",
            Some(region),
            old_cfg,
        )
        .expect("query");
        let old_reads: Vec<MappedRead> = old.by_ref().map(|r| r.expect("old ok")).collect();
        let old_counts = *old.filter_counts();

        // --- new path: SegmentMergedReads over pooled AlignmentFiles ---
        let filter = SegmentReadFilter {
            min_mapq: Some(20),
            min_read_length: Some(30),
            drop_qc_fail: true,
            drop_duplicate: true,
        };
        let files = [
            AlignmentFile::from_input(path_a, hdr.clone(), idx_a, None, filter, 0).expect("file a"),
            AlignmentFile::from_input(path_b, hdr.clone(), idx_b, None, filter, 1).expect("file b"),
        ];
        let mut new =
            SegmentMergedReads::new(&files, "chr1", region.start, region.end).expect("merge");
        let new_reads: Vec<MappedRead> = new.by_ref().map(|r| r.expect("new ok")).collect();
        let new_counts = new.filter_counts();

        // The gate: byte-identical reads (incl. source_file_index) + counts.
        assert_eq!(
            new_reads, old_reads,
            "merged reads must match the old query"
        );
        assert_eq!(
            new_counts, old_counts,
            "filter counts must match the old query"
        );

        // Sanity on the expected content.
        let names: Vec<String> = new_reads
            .iter()
            .map(|r| String::from_utf8_lossy(&r.qname).into_owned())
            .collect();
        assert_eq!(names, vec!["a_keep1", "b_keep1", "a_keep2", "b_keep2"]);
        assert_eq!(new_counts.low_mapq, 1);
        assert_eq!(new_counts.too_short, 1);
    }
}
