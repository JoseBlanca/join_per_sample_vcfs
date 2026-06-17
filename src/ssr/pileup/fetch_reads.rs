//! The per-locus read fetcher (arch §8) — index-queries one locus's reads,
//! reach-gates them, drops length-inconsistent records, and reservoir
//! depth-caps the survivors.
//!
//! - [`Reservoir`] — the read-level depth cap (§8.3), the net-new piece: the SNP
//!   side has only a *column* cap (`MPLP_MAX_DEPTH`), not a read-level reservoir.
//! - [`fetch_locus_reads`] — one locus's reads from a sample's per-worker
//!   [`WorkerReader`]s (the decode-caching CRAM / pooled BAM read source):
//!   segment-query the embedded window, apply the cheap coordinate-reach
//!   admission gate ([`reaches_locus`]), drop reads whose seq/qual/CIGAR lengths
//!   disagree, reservoir-cap the survivors, and surface the readers' filter
//!   drops for the `n_filtered` QC scalar.
//!
//! The catalog-walk loop that calls this per locus, and the worker stage that
//! delimits ([`super::alignment`]) + tallies ([`super::locus_tally::tally`]),
//! live in [`super::driver`].

use crate::bam::alignment_input::{FilterCounts, MappedRead};
use crate::bam::errors::AlignmentInputError;
use crate::bam::segment_reader::WorkerReader;
use crate::ssr::types::Locus;

use super::footprint::{cigar_read_len, reaches_locus};

/// Per-locus cap on **admitted** reads (arch §8.3/§10). A hypervariable,
/// high-depth locus is reservoir-sampled down to this many reads so one locus
/// can't make one worker's task (or one bundle) huge. A **calibration**
/// placeholder (arch §14) — a few hundred spanning reads already genotype a
/// locus; the cap only bites at pathologically deep loci.
pub(crate) const MAX_READS_PER_LOCUS: usize = 1000;

/// A tiny deterministic PRNG (SplitMix64) — seeded per locus so the depth-cap
/// subsample is reproducible and `--threads`-invariant (§8.4), with no external
/// RNG dependency whose stream could shift under us.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e3779b97f4a7615);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
}

/// Deterministic per-locus reservoir seed from `(chrom, start)` (arch §8.3):
/// FNV-1a over the contig name folded with the tract start. Stable for a given
/// locus across runs and thread counts — never derived from wall-clock or
/// thread id. Recorded in the `.ssr.psp` header so the subsample is reproducible.
pub(crate) fn locus_seed(chrom: &str, start: u32) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x100000001b3;
    let mut h = FNV_OFFSET;
    for &b in chrom.as_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h ^= start as u64;
    h.wrapping_mul(FNV_PRIME)
}

/// Reservoir sampler (Algorithm R) — an effectively-uniform sample of up to
/// `capacity` items from a stream of unknown length, in one pass with `O(capacity)`
/// memory (arch §8.3). The eviction index is drawn as `next_u64() % seen`, so the
/// sample carries a modulo bias bounded by `seen / 2^64` (negligible at any real
/// depth); the bias is accepted deliberately because the draw is **deterministic
/// and `--threads`-invariant**, which is the property the design needs (an
/// unbiased Lemire reduction would change the kept set and break that). The caller
/// `offer`s each **admitted** read in a fixed total order (`AlignmentMergedReader`'s
/// `(ref_id, pos)` → source-file → record order); with the deterministic per-locus
/// seed, the kept set is identical on every run and at every `--threads`. The
/// caller must not reorder the stream (§8.4).
pub(crate) struct Reservoir<T> {
    capacity: usize,
    held: Vec<T>,
    /// Admitted items offered so far (the `i` of Algorithm R).
    seen: u64,
    rng: SplitMix64,
}

impl<T> Reservoir<T> {
    pub(crate) fn new(capacity: usize, seed: u64) -> Self {
        Self {
            capacity,
            held: Vec::with_capacity(capacity),
            seen: 0,
            rng: SplitMix64::new(seed),
        }
    }

    /// Offer one admitted item. Keeps the first `capacity`; for the `i`-th item
    /// (`i > capacity`) keeps it with probability `capacity / i`, evicting one
    /// held item uniformly at random if kept.
    pub(crate) fn offer(&mut self, item: T) {
        self.seen += 1;
        if self.held.len() < self.capacity {
            self.held.push(item);
        } else {
            // j uniform in [0, seen); replace held[j] when it lands in-window.
            let j = (self.rng.next_u64() % self.seen) as usize;
            if j < self.capacity {
                self.held[j] = item;
            }
        }
    }

    /// The admitted depth (total items offered) — the reservoir sees only
    /// admitted reads, so this is `n_adm`, not the locus's raw depth.
    pub(crate) fn seen(&self) -> u64 {
        self.seen
    }

    /// Consume the reservoir, yielding the sampled items (≤ `capacity`).
    pub(crate) fn into_held(self) -> Vec<T> {
        self.held
    }
}

/// The reads a sample contributes at one locus, plus the fetch-pass tallies the
/// per-read aggregator cannot recover on its own.
pub(crate) struct LocusReads {
    /// The depth-capped reads admitted by the cheap reach gate
    /// ([`reaches_locus`]), concatenated across the sample's input files. The
    /// reservoir is order-independent, so per-file concatenation is sound
    /// (arch §8). These go to the worker for realignment.
    pub(crate) reads: Vec<MappedRead>,
    /// Reads the readers yielded at this locus — overlapping the window and
    /// past the flag/MAPQ/length filter — i.e. the pool the reach gate drew
    /// from (admitted + reach-discarded). Summed across files.
    pub(crate) yielded: u64,
    /// Reads dropped before reservoir admission because their seq/qual/CIGAR
    /// lengths disagree (empty/short `QUAL`, or a CIGAR whose read-consumption
    /// differs from `seq.len()`) — a truncated/malformed record the delimiter
    /// cannot slice safely. Folded into the `n_filtered` QC scalar.
    pub(crate) malformed: u64,
    /// The readers' flag/MAPQ/length drops at this locus, summed across files —
    /// the raw material for the `n_filtered` QC scalar (arch §3.1/§3.3). The
    /// reader owns the filter, so it owns the count
    /// ([`crate::bam::segment_reader::MappedReadsInSegment::filter_counts`]).
    pub(crate) filtered: FilterCounts,
}

/// Fetch one locus's reads from a sample's [`AlignmentFile`]s, gate them to the
/// plausibly-spanning ones, and depth-cap the survivors (arch §8.1/§8.3).
///
/// The query segment is the locus's embedded reference window (`ref_bytes`:
/// tract ± flank), converted from 0-based half-open to the reader's 1-based
/// inclusive convention. Every overlapping read past the readers' cheap filter
/// is offered to one per-locus [`Reservoir`] (seeded by [`locus_seed`] for a
/// `--threads`-invariant subsample) **iff** it passes [`reaches_locus`]; reads
/// that clearly cannot span are dropped from reservoir admission (the worker
/// never sees them — by design, the cap budget is spent on spanning evidence).
///
/// Multiple files for one sample feed the *same* reservoir (order-independent),
/// queried in input order for determinism.
///
/// # Errors
///
/// Propagates the first [`AlignmentInputError`] from opening/seeking a segment
/// or decoding a record (an unknown contig, a bad segment, or I/O).
pub(crate) fn fetch_locus_reads(
    readers: &mut [WorkerReader<'_>],
    locus: &Locus,
    cap: usize,
) -> Result<LocusReads, AlignmentInputError> {
    // Embedded reference window [ref_bytes_start, +len) is 0-based half-open;
    // the reader wants 1-based inclusive (arch §5/§8.1).
    let win_start = locus.ref_bytes_start();
    let seg_start = win_start + 1;
    let seg_end = win_start + locus.ref_bytes().len() as u32;

    let mut reservoir = Reservoir::new(cap, locus_seed(locus.chrom(), locus.start()));
    let mut yielded = 0u64;
    let mut malformed = 0u64;
    let mut filtered = FilterCounts::default();

    for reader in readers.iter_mut() {
        let (reads, counts) = reader.fetch_mapped_reads(locus.chrom(), seg_start, seg_end)?;
        for read in reads {
            yielded += 1;
            // Drop length-inconsistent records (empty/short `QUAL`, or a CIGAR
            // whose read-consumption ≠ `seq.len()`). The delimiter slices `seq`
            // and `qual` by ranges derived from `seq.len()`, so such a record
            // cannot be analyzed without corrupting the slice (arch §8.1; the
            // SNP path guards the same in `baq_engine`). Counted, not analyzed.
            if read.qual.len() != read.seq.len() || cigar_read_len(&read.cigar) != read.seq.len() {
                malformed += 1;
                continue;
            }
            if reaches_locus(&read, locus) {
                reservoir.offer(read);
            }
        }
        filtered.merge(&counts);
    }

    Ok(LocusReads {
        reads: reservoir.into_held(),
        yielded,
        malformed,
        filtered,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use std::fs::File;
    use std::num::NonZero;
    use std::path::{Path, PathBuf};
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
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record::{Flags, MappingQuality};
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::ReferenceSequence;
    use tempfile::TempDir;

    use crate::bam::index_preflight::AlignmentIndex;
    use crate::bam::segment_reader::{AlignmentFile, SegmentReadFilter};
    use crate::ssr::types::Motif;

    #[test]
    fn keeps_everything_when_offered_at_most_capacity() {
        let mut r = Reservoir::new(5, locus_seed("chr1", 100));
        for x in [10u32, 20, 30] {
            r.offer(x);
        }
        assert_eq!(r.seen(), 3);
        assert_eq!(r.into_held(), vec![10, 20, 30]); // first-K kept, in order
    }

    #[test]
    fn caps_at_capacity_and_counts_all_offers() {
        let mut r = Reservoir::new(10, locus_seed("chr1", 100));
        for x in 1..=100u32 {
            r.offer(x);
        }
        assert_eq!(r.seen(), 100);
        assert_eq!(r.into_held().len(), 10);
    }

    #[test]
    fn is_deterministic_for_a_fixed_seed_and_order() {
        let run = || {
            let mut r = Reservoir::new(10, locus_seed("chr7", 4242));
            for x in 1..=100u32 {
                r.offer(x);
            }
            r.into_held()
        };
        assert_eq!(run(), run());
    }

    #[test]
    fn reservoir_keeps_a_deterministic_subset_when_offered_far_past_capacity() {
        // Eviction branch taken ~10k times: the kept set is a subset of the
        // stream, sized at the cap, and identical across runs (it depends only on
        // seed + offer order — the basis of the cross-thread byte-identity claim).
        let run = || {
            let mut r = Reservoir::new(8, locus_seed("chrX", 7));
            for x in 1..=10_000u32 {
                r.offer(x);
            }
            (r.seen(), r.into_held())
        };
        let (seen, held) = run();
        assert_eq!(seen, 10_000);
        assert_eq!(held.len(), 8);
        assert!(
            held.iter().all(|x| (1..=10_000).contains(x)),
            "kept set is a subset of the stream"
        );
        assert_eq!(run().1, held, "kept set is identical across runs");
    }

    #[test]
    fn different_loci_sample_differently() {
        let sample = |chrom, start| {
            let mut r = Reservoir::new(10, locus_seed(chrom, start));
            for x in 1..=100u32 {
                r.offer(x);
            }
            r.into_held()
        };
        // Two distinct loci over the same stream pick (overwhelmingly likely)
        // different subsets — the seed actually drives the sampling.
        assert_ne!(sample("chr1", 100), sample("chr1", 101));
        assert_ne!(sample("chr1", 100), sample("chr2", 100));
    }

    #[test]
    fn every_item_can_be_selected_no_structural_exclusion() {
        // Over many seeds, the union of sampled items covers the whole stream —
        // the reservoir reaches every position, not a fixed prefix/suffix.
        let mut covered = HashSet::new();
        for seed in 0..500u64 {
            let mut r = Reservoir::new(10, seed);
            for x in 1..=100u32 {
                r.offer(x);
            }
            covered.extend(r.into_held());
        }
        assert_eq!(covered.len(), 100);
    }

    #[test]
    fn locus_seed_is_deterministic_and_distinguishes_loci() {
        assert_eq!(locus_seed("chr1", 100), locus_seed("chr1", 100));
        assert_ne!(locus_seed("chr1", 100), locus_seed("chr1", 101));
        assert_ne!(locus_seed("chr1", 100), locus_seed("chr2", 100));
    }

    // --- fetch_locus_reads (indexed-BAM fixture) ----------------------

    const CONTIG_LEN: usize = 200;

    /// A CA locus: GGGGGG | CACACA | TTTTTT — tract ref [16, 22), embedded
    /// window ref [10, 28) (so the 1-based inclusive query is [11, 28]).
    fn locus6() -> Locus {
        Locus::new(
            "chr1".into(),
            16,
            22,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGGGGCACACATTTTTT").into(),
            10,
        )
        .unwrap()
    }

    fn aln_record(qname: &str, start: usize, len: usize, mapq: u8) -> RecordBuf {
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

    /// As [`aln_record`] but with no quality scores (`QUAL = *`) — `qual` decodes
    /// to an empty buffer while `seq` is `len` long, so `qual.len() != seq.len()`.
    fn aln_record_no_qual(qname: &str, start: usize, len: usize, mapq: u8) -> RecordBuf {
        RecordBuf::builder()
            .set_name(qname.as_bytes())
            .set_reference_sequence_id(0)
            .set_flags(Flags::default())
            .set_mapping_quality(MappingQuality::new(mapq).expect("mapq"))
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar([Op::new(Kind::Match, len)].into_iter().collect())
            .set_sequence(Sequence::from(vec![b'A'; len]))
            .build()
    }

    fn build_csi_in_memory(bam_path: &Path) -> noodles_csi::Index {
        let mut reader = bam::io::reader::Builder
            .build_from_path(bam_path)
            .expect("open bam");
        let parsed_header = reader.read_header().expect("read header");
        let mut indexer: Indexer<BinnedIndex> = Indexer::default();
        let mut chunk_start = reader.get_ref().virtual_position();
        let mut record = bam::Record::default();
        while reader.read_record(&mut record).expect("read") != 0 {
            let chunk_end = reader.get_ref().virtual_position();
            let alignment_context = match (
                record.reference_sequence_id().transpose().expect("ref"),
                record.alignment_start().transpose().expect("start"),
                record.alignment_end().transpose().expect("end"),
            ) {
                (Some(id), Some(start), Some(end)) => {
                    Some((id, start, end, !record.flags().is_unmapped()))
                }
                _ => None,
            };
            indexer
                .add_record(alignment_context, Chunk::new(chunk_start, chunk_end))
                .expect("add");
            chunk_start = chunk_end;
        }
        indexer.build(parsed_header.reference_sequences().len())
    }

    /// Build a single-contig BAM-backed [`AlignmentFile`] with the default
    /// filter. Keeps the `TempDir` alive via the returned handle.
    fn bam_file(records: &[RecordBuf]) -> (TempDir, AlignmentFile) {
        let header = sam::Header::builder()
            .set_header(Default::default())
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZero::new(CONTIG_LEN).unwrap()),
            )
            .build();
        let dir = TempDir::new().expect("tempdir");
        let bam_path: PathBuf = dir.path().join("sample.bam");
        let file = File::create(&bam_path).expect("create bam");
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).expect("write header");
        for record in records {
            writer
                .write_alignment_record(&header, record)
                .expect("write record");
        }
        writer.try_finish().expect("finish");

        let csi = build_csi_in_memory(&bam_path);
        let alignment_file = AlignmentFile::from_input(
            bam_path,
            Arc::new(header),
            AlignmentIndex::BamCsi(Arc::new(csi)),
            None,
            SegmentReadFilter::default(),
            0,
        )
        .expect("alignment file");
        (dir, alignment_file)
    }

    #[test]
    fn fetch_admits_spanning_reads_gates_flanking_and_counts_reader_drops() {
        // Reads are ≥ DEFAULT_MIN_READ_LENGTH (30) so the reader keeps them;
        // the locus window is only 18 bp, so a spanning read overruns it.
        // span:  ref [5,35]  → brackets both tract ends → admitted.
        // flank: ref [17,47] → brackets the right end only, no clip → reach-skipped.
        // lowmq: would span, but MAPQ 5 < 20 → dropped by the reader filter,
        //        never yielded, counted in `filtered.low_mapq`.
        let records = [
            aln_record("span", 6, 30, 60),
            aln_record("lowmq", 6, 30, 5),
            aln_record("flank", 18, 30, 60),
        ];
        let (_dir, file) = bam_file(&records);

        let mut readers = vec![file.worker_reader()];
        let got = fetch_locus_reads(&mut readers, &locus6(), MAX_READS_PER_LOCUS).expect("fetch");

        let admitted: Vec<String> = got
            .reads
            .iter()
            .map(|r| String::from_utf8_lossy(&r.qname).into_owned())
            .collect();
        assert_eq!(admitted, vec!["span".to_string()]);
        // `span` + `flank` cleared the reader filter and overlapped the window;
        // `lowmq` did not (filtered out before it could be yielded).
        assert_eq!(got.yielded, 2);
        assert_eq!(got.filtered.low_mapq, 1);
    }

    #[test]
    fn fetch_drops_length_inconsistent_reads_and_counts_them() {
        // A spanning read with no quality scores (`qual.len() != seq.len()`) is
        // dropped before reservoir admission and tallied in `malformed`, so the
        // delimiter never slices it (B2). The clean read is still admitted.
        let records = [
            aln_record("good", 6, 30, 60),
            aln_record_no_qual("noqual", 6, 30, 60),
        ];
        let (_dir, file) = bam_file(&records);

        let mut readers = vec![file.worker_reader()];
        let got = fetch_locus_reads(&mut readers, &locus6(), MAX_READS_PER_LOCUS).expect("fetch");

        let admitted: Vec<String> = got
            .reads
            .iter()
            .map(|r| String::from_utf8_lossy(&r.qname).into_owned())
            .collect();
        assert_eq!(admitted, vec!["good".to_string()]);
        assert_eq!(got.yielded, 2);
        assert_eq!(got.malformed, 1);
    }

    #[test]
    fn fetch_concatenates_reads_across_a_samples_files() {
        // Two files, each with one spanning read → both reach the bundle.
        let (_d1, f1) = bam_file(&[aln_record("a", 6, 30, 60)]);
        let (_d2, f2) = bam_file(&[aln_record("b", 6, 30, 60)]);

        let mut readers = vec![f1.worker_reader(), f2.worker_reader()];
        let got = fetch_locus_reads(&mut readers, &locus6(), MAX_READS_PER_LOCUS).expect("fetch");

        let mut admitted: Vec<String> = got
            .reads
            .iter()
            .map(|r| String::from_utf8_lossy(&r.qname).into_owned())
            .collect();
        admitted.sort();
        assert_eq!(admitted, vec!["a".to_string(), "b".to_string()]);
        assert_eq!(got.yielded, 2);
    }
}
