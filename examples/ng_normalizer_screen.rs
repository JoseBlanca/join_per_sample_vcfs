//! **The differ-at-all screen** — algorithm plan 3 (`alignment_normalization.md`), Milestone D.
//!
//! The cheapest discriminating measurement the normalization comparison has: run all three
//! left-aligners over the **same real reads** and count, per algorithm pair, how many outputs
//! disagree. If that count is near zero, differently-spelled equivalents are not scattering a
//! variant's supporting reads, so normalization *placement* cannot explain a difference in calling
//! — and a whole avenue closes for the price of one run (spec §6, §10.3). This binary produces the
//! count; interpreting it and deciding whether the calling path adopts a different normalizer is a
//! different plan's job.
//!
//! ```text
//! cargo run --release --example ng_normalizer_screen -- <reference.fa> <reads.bam|cram> [reads2 ...]
//! ```
//!
//! The reference's sibling `<reference.fa>.fai` is used (the contig table is read from it — cheap,
//! no whole-genome pass). Reads are taken through ng's real ingestion path (`SampleReads`), so they
//! are the **filtered** reads a caller would actually see, and each read's reference window is
//! fetched on demand (`WindowedRefSeq`) rather than loading the genome into memory.
//!
//! Only reads whose CIGAR carries an indel are screened — a pure-match read is a no-op for every
//! normalizer and would agree trivially, inflating the denominator. The three normalizers:
//!
//! - **1a** [`StructuredLeftAligner`] — production's one structured pass.
//! - **1b** [`RepeatedLeftAligner`] — freebayes' repeated capped passes; a run that hits its cap is
//!   tallied (its output is not known to be leftmost).
//! - **1c** [`FixpointLeftAligner`] — 1a to a fixpoint, which **panics** on non-convergence; the
//!   panic is caught here so one pathological read tallies rather than aborting the whole screen.

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::PathBuf;
use std::process::ExitCode;

use pop_var_caller::ng::alignment::left_align_repeated::{ConvergenceReport, RepeatedLeftAligner};
use pop_var_caller::ng::alignment::left_align_structured::{
    FixpointLeftAligner, StructuredLeftAligner,
};
use pop_var_caller::ng::alignment::{Alignment, AlignmentNormalizer};
use pop_var_caller::ng::read::ReadFilterConfig;
use pop_var_caller::ng::read::input::SampleReads;
use pop_var_caller::ng::ref_seq::{RefSeq, WindowedRefSeq};
use pop_var_caller::ng::reference_info::{ReferenceSource, read_reference_info};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position};
use pop_var_caller::pileup::walker::CigarOp;

/// Everything the screen counts across the whole run.
///
/// Not `Copy`: it is a mutable accumulator, only ever taken `&mut`, and a silent by-value copy of a
/// counter aggregate is exactly the bug that shape invites.
#[derive(Debug, Default, Clone)]
struct Tally {
    /// Reads seen (post-filter, all CIGAR shapes).
    reads: u64,
    /// Reads whose CIGAR carries an indel — the ones a normalizer can move.
    reads_with_indel: u64,
    /// Reads that algorithm 1a actually moved from their input placement — how often the
    /// normalizers were genuinely exercised, so a zero-disagreement result cannot be dismissed as
    /// "the input was already leftmost, so all three no-op'd and agreed trivially".
    moved_by_normalization: u64,
    /// Per algorithm-pair disagreements, over `reads_with_indel`.
    disagree_1a_vs_1b: u64,
    disagree_1a_vs_1c: u64,
    disagree_1b_vs_1c: u64,
    /// 1b runs that hit the pass cap (output not known to be leftmost).
    b_hit_cap: u64,
    /// 1c runs that panicked (1a was not a fixpoint — would be a genuine finding).
    c_panicked: u64,
    /// Reads whose reference window could not be fetched (skipped).
    fetch_errors: u64,
    /// Reads the ingestion path yielded as an error (skipped).
    read_errors: u64,
}

/// The full `Display` chain of an error — top message plus each `source()` — so a failure names its
/// cause (which file, the OS reason) rather than just the outermost wrapper.
fn error_chain(error: &dyn std::error::Error) -> String {
    let mut message = error.to_string();
    let mut source = error.source();
    while let Some(cause) = source {
        message.push_str(": ");
        message.push_str(&cause.to_string());
        source = cause.source();
    }
    message
}

/// Whether a CIGAR carries at least one insertion or deletion.
fn has_indel(cigar: &[CigarOp]) -> bool {
    cigar
        .iter()
        .any(|op| matches!(op, CigarOp::Insertion(_) | CigarOp::Deletion(_)))
}

/// How many reference bases a CIGAR consumes — the length of the reference window the read spans.
fn reference_span(cigar: &[CigarOp]) -> u64 {
    cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Deletion(n)
            | CigarOp::Skip(n) => u64::from(*n),
            CigarOp::Insertion(_)
            | CigarOp::SoftClip(_)
            | CigarOp::HardClip(_)
            | CigarOp::Padding(_) => 0,
        })
        .sum()
}

/// Run all three normalizers on one indel-bearing alignment and fold the result into `tally`.
///
/// The caller has already confirmed the CIGAR carries an indel and counted the read. `alignment`
/// starts at offset 0 of `reference` (which is the read's reference window, from its aligned start);
/// `read` is the read's bases.
fn screen_one(alignment: &Alignment, read: &[u8], reference: &[u8], tally: &mut Tally) {
    let mut structured = alignment.clone();
    StructuredLeftAligner.normalize(&mut structured, read, reference);
    if &structured != alignment {
        tally.moved_by_normalization += 1;
    }

    let mut repeated = alignment.clone();
    let report = RepeatedLeftAligner::new().left_align(&mut repeated, read, reference);
    if report == ConvergenceReport::ExhaustedCap {
        tally.b_hit_cap += 1;
    }

    // 1c fails loudly (panics) on non-convergence — a genuine finding, but one it should never
    // reach, since 1a is a one-pass fixpoint. Caught so the screen tallies rather than aborts.
    let mut fixpoint = alignment.clone();
    let fixpoint_ok = catch_unwind(AssertUnwindSafe(|| {
        FixpointLeftAligner.normalize(&mut fixpoint, read, reference);
    }))
    .is_ok();
    if !fixpoint_ok {
        tally.c_panicked += 1;
    }

    if structured != repeated {
        tally.disagree_1a_vs_1b += 1;
    }
    if fixpoint_ok {
        if structured != fixpoint {
            tally.disagree_1a_vs_1c += 1;
        }
        if repeated != fixpoint {
            tally.disagree_1b_vs_1c += 1;
        }
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: ng_normalizer_screen <reference.fa> <reads.bam|cram> [reads2 ...]\n\
             counts, per normalizer pair, how many real reads the three left-aligners spell \
             differently."
        );
        return ExitCode::from(2);
    }
    let fasta = PathBuf::from(&args[1]);
    let read_paths: Vec<PathBuf> = args[2..].iter().map(PathBuf::from).collect();
    let fai = PathBuf::from(format!("{}.fai", fasta.display()));

    // The contig table from the .fai alone (names/order/lengths) — no whole-genome pass. The @SQ
    // assembly check is deferred and not needed for the screen.
    let reference_info = match read_reference_info(ReferenceSource::Fai(fai)) {
        Ok(info) => info,
        Err(error) => {
            eprintln!("could not read reference index: {}", error_chain(&error));
            return ExitCode::FAILURE;
        }
    };
    let contigs = reference_info.contig_list();

    let sample = match SampleReads::open(
        &read_paths,
        &reference_info,
        ReadFilterConfig::default(),
        false,
    ) {
        Ok(sample) => sample,
        Err(error) => {
            eprintln!("could not open reads: {}", error_chain(&error));
            return ExitCode::FAILURE;
        }
    };

    // One windowed reference for this loop's own per-read window fetches. The read stream's own
    // mismatch filter gets a fresh one per query through `make_reference`.
    let window_reference = WindowedRefSeq::new(fasta.clone(), contigs.clone());

    let mut tally = Tally::default();
    for (index, contig) in contigs.entries.iter().enumerate() {
        let region = GenomeRegion {
            contig: ContigId(index as u32),
            start: Position(1),
            end: Position(contig.length),
        };
        let stream = match sample.reads_in_region(region, || {
            WindowedRefSeq::new(fasta.clone(), contigs.clone())
        }) {
            Ok(stream) => stream,
            Err(error) => {
                eprintln!("skipping contig {}: {}", contig.name, error_chain(&error));
                continue;
            }
        };

        for item in stream {
            let read = match item {
                Ok(read) => read,
                Err(_) => {
                    tally.read_errors += 1;
                    continue;
                }
            };
            tally.reads += 1;
            if !has_indel(&read.cigar) {
                continue;
            }
            let span = reference_span(&read.cigar);
            if span == 0 {
                continue;
            }
            let window = match window_reference.fetch(ContigId(read.ref_id as u32), read.pos, span)
            {
                Ok(window) => window,
                Err(_) => {
                    tally.fetch_errors += 1;
                    continue;
                }
            };
            let alignment = Alignment {
                reference_offset: 0,
                cigar: read.cigar.clone(),
            };
            tally.reads_with_indel += 1;
            screen_one(&alignment, &read.seq, &window, &mut tally);
        }

        if tally.reads_with_indel > 0 {
            eprintln!(
                "... {} ({} reads, {} with an indel so far)",
                contig.name, tally.reads, tally.reads_with_indel
            );
        }
    }

    report(&tally);
    ExitCode::SUCCESS
}

/// Print the screen's result — the disagreement counts and the two "not-known-leftmost" tallies.
fn report(tally: &Tally) {
    let pct = |n: u64| {
        if tally.reads_with_indel == 0 {
            0.0
        } else {
            100.0 * n as f64 / tally.reads_with_indel as f64
        }
    };
    println!("--- ng normalizer differ-at-all screen ---");
    println!("reads seen:            {}", tally.reads);
    println!("reads with an indel:   {}", tally.reads_with_indel);
    println!(
        "moved by normalization: {}  ({:.4}%)  <- how often a normalizer was actually exercised",
        tally.moved_by_normalization,
        pct(tally.moved_by_normalization)
    );
    println!(
        "1a vs 1b disagree:     {}  ({:.4}%)",
        tally.disagree_1a_vs_1b,
        pct(tally.disagree_1a_vs_1b)
    );
    println!(
        "1a vs 1c disagree:     {}  ({:.4}%)",
        tally.disagree_1a_vs_1c,
        pct(tally.disagree_1a_vs_1c)
    );
    println!(
        "1b vs 1c disagree:     {}  ({:.4}%)",
        tally.disagree_1b_vs_1c,
        pct(tally.disagree_1b_vs_1c)
    );
    println!(
        "1b hit its pass cap:   {}  ({:.4}%)",
        tally.b_hit_cap,
        pct(tally.b_hit_cap)
    );
    println!("1c panicked (non-fixpoint 1a): {}", tally.c_panicked);
    println!(
        "skipped (fetch/read errors):   {} / {}",
        tally.fetch_errors, tally.read_errors
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use pop_var_caller::pileup::walker::CigarOp::{Deletion, Insertion, Match};

    fn alignment(cigar: Vec<CigarOp>) -> Alignment {
        Alignment {
            reference_offset: 0,
            cigar,
        }
    }

    #[test]
    fn has_indel_and_reference_span_read_the_cigar() {
        assert!(!has_indel(&[Match(5)]));
        assert!(has_indel(&[Match(4), Deletion(1), Match(1)]));
        assert!(has_indel(&[Match(4), Insertion(1), Match(1)]));
        // Deletion consumes reference, insertion does not.
        assert_eq!(reference_span(&[Match(4), Deletion(1), Match(1)]), 6);
        assert_eq!(reference_span(&[Match(4), Insertion(1), Match(1)]), 5);
    }

    #[test]
    fn a_homopolymer_indel_has_the_three_normalizers_agree() {
        // ref GAAAAT, read GAAAT: a pure deletion the bad-CIGAR filter keeps. All three place it at
        // the leftmost A, so no pair disagrees.
        let mut tally = Tally::default();
        screen_one(
            &alignment(vec![Match(4), Deletion(1), Match(1)]),
            b"GAAAT",
            b"GAAAAT",
            &mut tally,
        );
        assert_eq!(tally.disagree_1a_vs_1b, 0);
        assert_eq!(tally.disagree_1a_vs_1c, 0);
        assert_eq!(tally.disagree_1b_vs_1c, 0);
        assert_eq!(tally.b_hit_cap, 0);
        assert_eq!(tally.c_panicked, 0);
        // The rightmost placement is not leftmost, so 1a moved it — the counter that proves a
        // zero-disagreement result is not a trivial all-no-op.
        assert_eq!(tally.moved_by_normalization, 1);
    }

    #[test]
    fn an_indel_placed_past_1b_s_pass_cap_makes_1b_disagree() {
        // The disagreement the synthetic stress reads are built to trigger, and the one real-read
        // shape that survives the filters: a deletion 24 bases from its leftmost home (a 25-base
        // homopolymer, rightmost placement). 1a and 1c reach leftmost in one pass; 1b shifts one
        // base per pass and caps at 20, stopping short — so it disagrees with both.
        let reference: Vec<u8> = std::iter::once(b'G')
            .chain(std::iter::repeat_n(b'A', 25))
            .chain(std::iter::once(b'T'))
            .collect();
        let read: Vec<u8> = std::iter::once(b'G')
            .chain(std::iter::repeat_n(b'A', 24))
            .chain(std::iter::once(b'T'))
            .collect();
        let mut tally = Tally::default();
        screen_one(
            &alignment(vec![Match(25), Deletion(1), Match(1)]),
            &read,
            &reference,
            &mut tally,
        );
        assert_eq!(tally.disagree_1a_vs_1b, 1, "1b caps short of leftmost");
        assert_eq!(tally.disagree_1a_vs_1c, 0, "1c reaches leftmost like 1a");
        assert_eq!(tally.disagree_1b_vs_1c, 1);
        assert_eq!(tally.b_hit_cap, 1);
        assert_eq!(tally.moved_by_normalization, 1);
    }

    #[test]
    fn an_indel_exactly_at_the_cap_hits_it_but_still_agrees() {
        // The subtle boundary the synthetic run also exercises: an insertion 20 bases from home. 1b
        // shifts 20 times, reaching leftmost on the last pass but with no pass left to *confirm* it,
        // so it reports ExhaustedCap — yet the output IS leftmost, so no pair disagrees. This is why
        // the run's cap-hits (36) exceed its disagreements (28).
        let reference: Vec<u8> = std::iter::once(b'G')
            .chain(std::iter::repeat_n(b'A', 20))
            .chain(std::iter::once(b'T'))
            .collect();
        let read: Vec<u8> = std::iter::once(b'G')
            .chain(std::iter::repeat_n(b'A', 21))
            .chain(std::iter::once(b'T'))
            .collect();
        let mut tally = Tally::default();
        screen_one(
            &alignment(vec![Match(21), Insertion(1), Match(1)]),
            &read,
            &reference,
            &mut tally,
        );
        assert_eq!(
            tally.b_hit_cap, 1,
            "1b needs a 21st pass to confirm, but the cap is 20"
        );
        assert_eq!(tally.disagree_1a_vs_1b, 0, "yet 1b did reach leftmost");
        assert_eq!(tally.disagree_1a_vs_1c, 0);
        assert_eq!(tally.disagree_1b_vs_1c, 0);
    }

    #[test]
    fn a_complex_indel_overlap_makes_1a_and_1b_disagree() {
        // The one shape 1a and 1b spell differently: an overlapping D/I 1a trims and 1b does not.
        // (A real ingestion path drops this via the bad-CIGAR filter, so the screen rarely sees it;
        // exercised here directly to prove the screen detects a genuine disagreement when present.)
        let mut tally = Tally::default();
        screen_one(
            &alignment(vec![
                Match(1),
                Deletion(1),
                Deletion(1),
                Insertion(1),
                Match(1),
                Match(1),
            ]),
            b"AAAC",
            b"AAAAC",
            &mut tally,
        );
        assert_eq!(
            tally.disagree_1a_vs_1b, 1,
            "1a trims the overlap, 1b does not"
        );
        // 1c is 1a to a fixpoint, so it matches 1a and disagrees with 1b the same way.
        assert_eq!(tally.disagree_1a_vs_1c, 0);
        assert_eq!(tally.disagree_1b_vs_1c, 1);
        assert_eq!(tally.c_panicked, 0);
    }
}
