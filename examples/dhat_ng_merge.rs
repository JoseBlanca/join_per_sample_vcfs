//! **T14 — the k-way merge's per-read budget.**
//!
//! Build and run inside the dev container:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example dhat_ng_merge --features dhat-heap
//! ```
//!
//! Spec `sample_reads.md` §3.2 holds the merge to *O(k) integer comparisons and
//! one `MappedRead` move per read, allocating nothing per read*. A stray
//! `clone()` or a redundant move there costs throughput **without changing a
//! single output byte**, so no correctness test can see it — which is why this
//! check exists separately from them.
//!
//! ## Why dhat, and why an example
//!
//! The project forbids `unsafe` (`Cargo.toml`, `unsafe_code = "forbid"`), so it
//! cannot install a counting allocator of its own — `GlobalAlloc` is an unsafe
//! trait. `dhat::Alloc` provides one whose unsafe lives inside dhat, which is
//! the route the repo already uses for its other heap profiles
//! (`examples/dhat_*.rs`). Feature-gated, so ordinary builds carry no dhat code;
//! `clippy --all-targets` still compiles this file, so it cannot bit-rot
//! silently.
//!
//! ## How the merge's own cost is isolated
//!
//! Everything below the merge allocates per read *by design*: decoding a record
//! builds a `MappedRead`, which owns its sequence, qualities and CIGAR. So
//! "zero allocations while draining a merged stream" is not a property anything
//! could satisfy, and asserting it would prove nothing.
//!
//! What is measurable is the **difference**. The same reads are drained twice:
//! once from one file, which takes the merge-free `Single` arm, and once split
//! across two files, which takes `Merged`. Both decode the same reads through
//! the same filter, so the delta is the merge and nothing else. Run at two
//! sizes, the delta must stay **flat**: if the merge allocated per read,
//! quadrupling the reads would quadruple it.
//!
//! Comparing shapes rather than absolute numbers keeps this from being a
//! brittle golden-number check that needs editing whenever something below the
//! merge changes.

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::fs::File;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles_bam as bam;
use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::{Flags, MappingQuality};
use noodles_sam::alignment::record_buf::{QualityScores, Sequence};
use sam::header::record::value::Map;
use sam::header::record::value::map::{ReadGroup, ReferenceSequence};

use pop_var_caller::bam::index_preflight::preflight_alignment_indexes;
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::ref_seq::InMemoryRefSeq;
use pop_var_caller::ng::reference_info::{ReferenceInfo, ReferenceSource, read_reference_info};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};

const CONTIG: &str = "chr1";
const CONTIG_LENGTH: usize = 100_000;
const READ_LENGTH: usize = 30;

fn main() {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    #[cfg(not(feature = "dhat-heap"))]
    {
        eprintln!(
            "This example measures allocations and needs the dhat allocator.\n\
             Re-run with: cargo run --release --example dhat_ng_merge --features dhat-heap"
        );
        return;
    }

    #[cfg(feature = "dhat-heap")]
    {
        // Warm-up, so first-touch lazy initialisation is not charged to
        // whichever measurement happens to run first.
        let warm = Fixture::build(64);
        warm.drain_single();
        warm.drain_merged();

        let small = Fixture::build(2_000);
        let (small_reads, single_small) = measured(|| small.drain_single());
        let (small_merged_reads, merged_small) = measured(|| small.drain_merged());
        assert_eq!(
            small_reads, small_merged_reads,
            "both arms must yield the same reads, or the comparison is meaningless"
        );

        let large = Fixture::build(8_000);
        let (large_reads, single_large) = measured(|| large.drain_single());
        let (large_merged_reads, merged_large) = measured(|| large.drain_merged());
        assert_eq!(large_reads, large_merged_reads);
        assert_eq!(
            large_reads,
            small_reads * 4,
            "the large fixture must be 4x the small one for the scaling to mean anything"
        );

        let delta_small = merged_small.saturating_sub(single_small);
        let delta_large = merged_large.saturating_sub(single_large);

        println!("reads:              {small_reads} -> {large_reads} (4x)");
        println!("single-arm blocks:  {single_small} -> {single_large}");
        println!("merged-arm blocks:  {merged_small} -> {merged_large}");
        println!("merge's own blocks: {delta_small} -> {delta_large}");
        println!(
            "per read:           {:.4} -> {:.4}",
            delta_small as f64 / small_reads as f64,
            delta_large as f64 / large_reads as f64
        );

        // A per-read allocation would make the delta grow with the reads. The
        // merge's only allocations are the `Vec`s it builds once per query, so
        // the delta must stay flat; the slack absorbs those and any noise.
        let would_be_per_read = delta_small.max(1) * 3;
        assert!(
            delta_large <= would_be_per_read.max(64),
            "the merge's allocation cost grew with the read count \
             ({delta_small} -> {delta_large} over a 4x increase): it is \
             allocating per read, which spec §3.2 forbids — look for a \
             `clone()` or a redundant move in `MergedRegionReads::next`"
        );

        println!("\nOK: the merge's allocation cost is flat in the read count.");
    }
}

/// Total heap blocks allocated so far, and the value `body` produced.
#[cfg(feature = "dhat-heap")]
fn measured<T>(body: impl FnOnce() -> T) -> (T, u64) {
    let before = dhat::HeapStats::get().total_blocks;
    let value = body();
    let after = dhat::HeapStats::get().total_blocks;
    (value, after - before)
}

struct Fixture {
    _dir: tempfile::TempDir,
    _fasta_dir: tempfile::TempDir,
    single: Vec<PathBuf>,
    pair: Vec<PathBuf>,
    reference: ReferenceInfo,
}

impl Fixture {
    /// One file holding `2 * per_file` reads, and two files holding `per_file`
    /// each — the same total, one taking the merge-free arm and one the merge.
    ///
    /// The two files **interleave**: one takes the even slots, the other the
    /// odd ones, so the merge alternates between them at essentially every
    /// step. That is the shape a sample of several experiments actually
    /// produces, and it is what exercises the argmin rather than draining one
    /// file and then the other.
    fn build(per_file: usize) -> Self {
        let fasta_dir = tempfile::tempdir().expect("tempdir");
        let fasta = write_fasta(fasta_dir.path());
        let reference = read_reference_info(ReferenceSource::Fasta {
            fasta: fasta.clone(),
            fai: None,
        })
        .expect("read reference");

        let dir = tempfile::tempdir().expect("tempdir");
        let first = reads(per_file, 0);
        let second = reads(per_file, 1);
        let mut both = first.clone();
        both.extend(second.clone());
        both.sort_by_key(start_of);

        Self {
            single: vec![write_indexed_bam(dir.path(), "solo", &both)],
            pair: vec![
                write_indexed_bam(dir.path(), "a", &first),
                write_indexed_bam(dir.path(), "b", &second),
            ],
            reference,
            _dir: dir,
            _fasta_dir: fasta_dir,
        }
    }

    fn drain_single(&self) -> usize {
        drain(&self.single, &self.reference)
    }

    fn drain_merged(&self) -> usize {
        drain(&self.pair, &self.reference)
    }
}

fn drain(paths: &[PathBuf], reference: &ReferenceInfo) -> usize {
    let sample = SampleReads::open(paths, reference, ReadFilterConfig::default(), false)
        .expect("the fixture files open");
    let region = GenomeRegion {
        contig: ContigId(0),
        start: Position(1),
        end: Position(CONTIG_LENGTH as u64),
    };
    sample
        .reads_in_region(region, || {
            InMemoryRefSeq::from_contigs(vec![vec![b'A'; CONTIG_LENGTH]])
        })
        .expect("query")
        .inspect(|item| {
            // Unwrapped rather than ignored: an error here would silently make
            // the two arms drain different numbers of reads, and the whole
            // comparison rests on them draining the same ones.
            assert!(item.is_ok(), "no fatal error while draining");
        })
        .count()
}

fn start_of(record: &RecordBuf) -> usize {
    record.alignment_start().map(usize::from).unwrap_or(0)
}

/// `count` coordinate-sorted reads, `offset` picking the even or odd slots.
fn reads(count: usize, offset: usize) -> Vec<RecordBuf> {
    let span = CONTIG_LENGTH - READ_LENGTH - 2;
    let mut records: Vec<RecordBuf> = (0..count)
        .map(|i| {
            let start = 1 + ((i * 2 + offset) % span);
            read_at(&format!("r{offset}_{i}"), start)
        })
        .collect();
    // The open gate requires coordinate order.
    records.sort_by_key(start_of);
    records
}

fn read_at(qname: &str, start: usize) -> RecordBuf {
    RecordBuf::builder()
        .set_name(qname.as_bytes())
        .set_reference_sequence_id(0usize)
        // Explicit: `RecordBuf`'s default flags are UNMAPPED, which filter #1
        // drops — leaving nothing to merge and nothing to measure.
        .set_flags(Flags::empty())
        .set_mapping_quality(MappingQuality::new(60).expect("mapq in range"))
        .set_alignment_start(noodles_core::Position::try_from(start).expect("1-based"))
        .set_cigar([Op::new(Kind::Match, READ_LENGTH)].into_iter().collect())
        .set_sequence(Sequence::from(vec![b'A'; READ_LENGTH]))
        .set_quality_scores(QualityScores::from(vec![30u8; READ_LENGTH]))
        .build()
}

fn write_fasta(dir: &Path) -> PathBuf {
    use std::io::Write as _;

    let fasta = dir.join("ref.fa");
    let fai = dir.join("ref.fa.fai");
    let mut sequence_file = File::create(&fasta).expect("create fasta");
    let mut index_file = File::create(&fai).expect("create fai");

    writeln!(sequence_file, ">{CONTIG}").expect("write header");
    sequence_file
        .write_all(&vec![b'A'; CONTIG_LENGTH])
        .expect("write bases");
    sequence_file.write_all(b"\n").expect("write newline");
    writeln!(
        index_file,
        "{CONTIG}\t{CONTIG_LENGTH}\t{}\t{CONTIG_LENGTH}\t{}",
        CONTIG.len() + 2,
        CONTIG_LENGTH + 1
    )
    .expect("write fai");

    fasta
}

fn header_naming(sample: &str) -> sam::Header {
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::SAMPLE;

    let mut hd = Map::<sam::header::record::value::map::Header>::default();
    hd.other_fields_mut()
        .insert(SORT_ORDER, b"coordinate".as_ref().into());

    let mut read_group = Map::<ReadGroup>::default();
    read_group
        .other_fields_mut()
        .insert(SAMPLE, sample.as_bytes().into());

    sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(
            CONTIG,
            Map::<ReferenceSequence>::new(NonZero::new(CONTIG_LENGTH).expect("non-zero")),
        )
        .add_read_group("rg1", read_group)
        .build()
}

fn write_indexed_bam(dir: &Path, name: &str, records: &[RecordBuf]) -> PathBuf {
    // Every file names the *same* sample: they are one sample's experiments,
    // which is the case the merge exists for.
    let header = header_naming("NA12878");
    let path = dir.join(format!("{name}.bam"));

    let mut writer = bam::io::Writer::new(File::create(&path).expect("create bam"));
    writer.write_header(&header).expect("write header");
    for record in records {
        writer
            .write_alignment_record(&header, record)
            .expect("write record");
    }
    writer.try_finish().expect("finish");

    preflight_alignment_indexes(std::slice::from_ref(&path), true).expect("build index");
    path
}
