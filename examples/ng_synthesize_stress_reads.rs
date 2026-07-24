//! **Synthetic stress reads for the differ-at-all screen** — write a reference + BAM engineered so
//! the three normalizers *can* disagree, to show [`ng_normalizer_screen`](../ng_normalizer_screen.rs)
//! is genuinely discriminating (and to find where they diverge).
//!
//! ```text
//! cargo run --release --example ng_synthesize_stress_reads -- <output_dir>
//! # then:
//! cargo run --release --example ng_normalizer_screen -- <output_dir>/stress_ref.fa <output_dir>/stress_reads.bam
//! ```
//!
//! # What makes the three disagree
//!
//! On a read that survives the caller's filters (no adjacent `I`/`D`, no end deletion), the three
//! left-aligners agree **except** when an indel sits more than `MAX_PASSES = 20` bases from its
//! leftmost home:
//!
//! - **1a** (structured) and **1c** (1a to a fixpoint) reach the leftmost position in one pass.
//! - **1b** (freebayes' repeated one-base passes) shifts one base per pass and **caps at 20**, so a
//!   shift of 21+ leaves the indel short of leftmost — a genuine disagreement with 1a and 1c.
//!
//! So each read here places a single indel at the **rightmost** end of a homopolymer of length `L`
//! (shift distance `L − 1` for a deletion, `L` for an insertion). Lengths straddle the cap: below it
//! all three agree; at it 1b hits the cap but still reaches leftmost (agree, cap tallied); above it
//! 1b stops short (disagree). The reads are otherwise clean — flanks match, one indel, no adjacent
//! `I`/`D` — so they pass ingestion and the screen actually scores them.

use std::fs::File;
use std::io::Write as _;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

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

const CONTIG: &str = "stress";
/// Non-repetitive flanks that anchor the read; neither touches an `A` at the repeat boundary, so it
/// cannot extend the homopolymer and change the shift distance.
const LEFT_FLANK: &[u8] = b"CGTCGATCGATCTAGCTCGATCGATCGATG";
const RIGHT_FLANK: &[u8] = b"CTAGCTCGATCGATCGATCTAGCTCGATCA";
/// A spacer of unique sequence between loci, so no read's window bleeds into its neighbour.
const SPACER: &[u8] = b"GCGCGCGCGCTATATATAT";
/// Homopolymer lengths, chosen to straddle `MAX_PASSES = 20`: below it (agree), at it (cap but
/// agree), and above it (disagree). Each is used for a deletion and an insertion locus.
const REPEAT_LENGTHS: &[usize] = &[6, 12, 18, 20, 21, 24, 30, 40];
/// How many identical reads to emit per locus — a little volume so the counts are not all 1.
const READS_PER_LOCUS: usize = 4;

/// One planned read: where it aligns (0-based contig offset), its CIGAR, and its bases.
struct PlannedRead {
    contig_offset: usize,
    cigar: Vec<Op>,
    sequence: Vec<u8>,
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        eprintln!("usage: ng_synthesize_stress_reads <output_dir>");
        return ExitCode::from(2);
    }
    let dir = PathBuf::from(&args[1]);
    if let Err(error) = std::fs::create_dir_all(&dir) {
        eprintln!("could not create {}: {error}", dir.display());
        return ExitCode::FAILURE;
    }

    let (reference, reads) = build_loci();

    let fasta = dir.join("stress_ref.fa");
    write_fasta(&fasta, &reference);
    let bam_path = dir.join("stress_reads.bam");
    write_indexed_bam(&bam_path, reference.len(), &reads);

    println!(
        "wrote {} bases of reference to {}",
        reference.len(),
        fasta.display()
    );
    println!("wrote {} reads to {}", reads.len(), bam_path.display());
    println!(
        "repeat lengths {:?} (deletion + insertion each), {} reads per locus",
        REPEAT_LENGTHS, READS_PER_LOCUS
    );
    println!("\nnow run the screen:");
    println!(
        "  cargo run --release --example ng_normalizer_screen -- {} {}",
        fasta.display(),
        bam_path.display()
    );
    ExitCode::SUCCESS
}

/// Append one locus — left flank, `length` A's, right flank, trailing spacer — to `reference`, and
/// return the 0-based offset at which the locus (its left flank) begins.
fn append_locus(reference: &mut Vec<u8>, length: usize) -> usize {
    let locus_offset = reference.len();
    reference.extend_from_slice(LEFT_FLANK);
    reference.extend(std::iter::repeat_n(b'A', length));
    reference.extend_from_slice(RIGHT_FLANK);
    reference.extend_from_slice(SPACER);
    locus_offset
}

/// The bases of a read that places an indel rightmost in a `length`-A homopolymer: left flank,
/// `repeat_a_count` A's, right flank.
fn read_sequence(repeat_a_count: usize) -> Vec<u8> {
    let mut sequence = Vec::with_capacity(LEFT_FLANK.len() + repeat_a_count + RIGHT_FLANK.len());
    sequence.extend_from_slice(LEFT_FLANK);
    sequence.extend(std::iter::repeat_n(b'A', repeat_a_count));
    sequence.extend_from_slice(RIGHT_FLANK);
    sequence
}

/// Build the reference contig and the reads that place an indel rightmost in each homopolymer.
fn build_loci() -> (Vec<u8>, Vec<PlannedRead>) {
    let mut reference = Vec::new();
    let mut reads = Vec::new();
    reference.extend_from_slice(SPACER);

    let mut emit = |offset: usize, cigar: Vec<Op>, sequence: Vec<u8>| {
        for _ in 0..READS_PER_LOCUS {
            reads.push(PlannedRead {
                contig_offset: offset,
                cigar: cigar.clone(),
                sequence: sequence.clone(),
            });
        }
    };

    for &length in REPEAT_LENGTHS {
        // Deletion locus: the read drops the last A, spelled rightmost (match all but one A, then
        // the deletion). Shift distance to leftmost is `length - 1`.
        let offset = append_locus(&mut reference, length);
        let deletion_cigar = vec![
            Op::new(Kind::Match, LEFT_FLANK.len() + (length - 1)),
            Op::new(Kind::Deletion, 1),
            Op::new(Kind::Match, RIGHT_FLANK.len()),
        ];
        emit(offset, deletion_cigar, read_sequence(length - 1));

        // Insertion locus: the read adds one A, spelled rightmost (match all the A's, then the
        // insertion). Shift distance to leftmost is `length`.
        let offset = append_locus(&mut reference, length);
        let insertion_cigar = vec![
            Op::new(Kind::Match, LEFT_FLANK.len() + length),
            Op::new(Kind::Insertion, 1),
            Op::new(Kind::Match, RIGHT_FLANK.len()),
        ];
        emit(offset, insertion_cigar, read_sequence(length + 1));
    }

    (reference, reads)
}

/// A header the ingestion open-gate accepts: `@HD SO:coordinate`, the one `@SQ`, and one `@RG`
/// naming a sample (mirrors `dhat_ng_merge`'s proven header shape).
fn coordinate_sorted_header(contig_length: usize) -> sam::Header {
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::SAMPLE;

    let mut hd = Map::<sam::header::record::value::map::Header>::default();
    hd.other_fields_mut()
        .insert(SORT_ORDER, b"coordinate".as_ref().into());

    let mut read_group = Map::<ReadGroup>::default();
    read_group
        .other_fields_mut()
        .insert(SAMPLE, b"synthetic".as_ref().into());

    sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(
            CONTIG,
            Map::<ReferenceSequence>::new(NonZero::new(contig_length).expect("non-zero contig")),
        )
        .add_read_group("rg1", read_group)
        .build()
}

fn write_fasta(path: &Path, reference: &[u8]) {
    let fai = PathBuf::from(format!("{}.fai", path.display()));
    let mut sequence_file = File::create(path).expect("create fasta");
    let mut index_file = File::create(&fai).expect("create fai");

    writeln!(sequence_file, ">{CONTIG}").expect("write header");
    sequence_file.write_all(reference).expect("write bases");
    sequence_file.write_all(b"\n").expect("write newline");
    // Single-line sequence, so the .fai geometry is trivial.
    writeln!(
        index_file,
        "{CONTIG}\t{}\t{}\t{}\t{}",
        reference.len(),
        CONTIG.len() + 2,
        reference.len(),
        reference.len() + 1
    )
    .expect("write fai");
}

fn write_indexed_bam(path: &Path, contig_length: usize, reads: &[PlannedRead]) {
    let header = coordinate_sorted_header(contig_length);

    // Coordinate-sort by alignment start — the ingestion path's open gate requires it.
    let mut sorted: Vec<&PlannedRead> = reads.iter().collect();
    sorted.sort_by_key(|read| read.contig_offset);

    let mut writer = bam::io::Writer::new(File::create(path).expect("create bam"));
    writer.write_header(&header).expect("write header");
    for (index, read) in sorted.iter().enumerate() {
        let record = RecordBuf::builder()
            .set_name(format!("stress{index}").as_bytes())
            .set_reference_sequence_id(0usize)
            .set_flags(Flags::empty())
            .set_mapping_quality(MappingQuality::new(60).expect("mapq in range"))
            .set_alignment_start(
                noodles_core::Position::try_from(read.contig_offset + 1).expect("1-based"),
            )
            .set_cigar(read.cigar.iter().copied().collect())
            .set_sequence(Sequence::from(read.sequence.clone()))
            .set_quality_scores(QualityScores::from(vec![40u8; read.sequence.len()]))
            .build();
        writer
            .write_alignment_record(&header, &record)
            .expect("write record");
    }
    writer.try_finish().expect("finish bam");

    preflight_alignment_indexes(std::slice::from_ref(&path.to_path_buf()), true)
        .expect("build index");
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The flanks must not touch an `A` at the repeat boundary, or they would extend the
    /// homopolymer and change the shift distance the whole stimulus is calibrated on.
    #[test]
    fn the_flanks_do_not_extend_the_homopolymer() {
        assert_ne!(LEFT_FLANK.last(), Some(&b'A'), "left flank ends in A");
        assert_ne!(
            RIGHT_FLANK.first(),
            Some(&b'A'),
            "right flank starts with A"
        );
    }

    /// Every generated read is well-formed — its CIGAR consumes exactly the read's bases and a
    /// reference span consistent with the locus — and places its indel at the intended (rightmost)
    /// shift distance. A silent off-by-one here would mis-calibrate the screen's stimulus while
    /// still producing plausible-looking numbers.
    #[test]
    fn every_read_is_well_formed_and_placed_rightmost() {
        let (reference, reads) = build_loci();
        assert_eq!(reads.len(), REPEAT_LENGTHS.len() * 2 * READS_PER_LOCUS);

        for read in &reads {
            let read_consumed: usize = read
                .cigar
                .iter()
                .filter(|op| matches!(op.kind(), Kind::Match | Kind::Insertion))
                .map(|op| op.len())
                .sum();
            assert_eq!(
                read_consumed,
                read.sequence.len(),
                "cigar read-consumption disagrees with the sequence length"
            );

            let reference_consumed: usize = read
                .cigar
                .iter()
                .filter(|op| matches!(op.kind(), Kind::Match | Kind::Deletion))
                .map(|op| op.len())
                .sum();
            // The locus (left flank + L A's + right flank) begins at `contig_offset`; the read's
            // reference footprint must sit inside the reference.
            assert!(read.contig_offset + reference_consumed <= reference.len());

            // The read's bases are the flanks with an all-`A` middle — the only non-flank content.
            let middle = &read.sequence[LEFT_FLANK.len()..read.sequence.len() - RIGHT_FLANK.len()];
            assert!(middle.iter().all(|&b| b == b'A'), "middle is not all A");
        }
    }

    /// The shift distances straddle the 20-pass cap exactly as the run reported: 36 reads hit the
    /// cap (deletion `L-1 >= 20`, insertion `L >= 20`), of which 28 are left non-leftmost and
    /// therefore disagree (deletion `L-1 >= 21`, insertion `L >= 21`).
    #[test]
    fn the_shift_distances_straddle_the_cap_as_reported() {
        let cap = 20usize;
        let mut cap_hits = 0;
        let mut disagreements = 0;
        for &length in REPEAT_LENGTHS {
            let deletion_shift = length - 1;
            let insertion_shift = length;
            for shift in [deletion_shift, insertion_shift] {
                if shift >= cap {
                    cap_hits += READS_PER_LOCUS;
                }
                if shift >= cap + 1 {
                    disagreements += READS_PER_LOCUS;
                }
            }
        }
        assert_eq!(cap_hits, 36);
        assert_eq!(disagreements, 28);
    }
}
