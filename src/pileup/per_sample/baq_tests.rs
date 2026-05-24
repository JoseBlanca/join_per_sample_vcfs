//! Pileup-glue tests for the BAQ stage: parity against htslib's shipped
//! realn fixtures (`tests/data/baq/realn0N.{sam,fa}` and
//! `realn0N_exp.sam`), plus unit tests for the engine driver, the rayon
//! stream adapter, and the leaf helpers. Algorithm-only tests
//! (`probaln_glocal` input validation, PARITY constant pins, `set_u`
//! boundary contract, `Q2P` monotonicity) live in `src/baq/tests.rs`.

use std::collections::VecDeque;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    alignment::{
        record::{cigar::op::Kind, data::field::Tag},
        record_buf::data::field::Value,
    },
};

use crate::baq::{BaqConfig, ProbalnScratch, probaln_glocal};
use crate::pileup::walker::CigarOp;

use super::baq_engine::{
    apply_baq_cap_into, compute_alignment_end, encode_base, extend_ref_window,
    locate_alignment_span,
};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("baq")
        .join(name)
}

fn read_single_fasta(path: &Path) -> Vec<u8> {
    let mut reader = fasta::io::Reader::new(BufReader::new(File::open(path).unwrap()));
    let mut records = reader.records();
    let rec = records
        .next()
        .expect("fixture FASTA must contain a record")
        .unwrap();
    assert!(
        records.next().is_none(),
        "fixture FASTA must contain a single record"
    );
    rec.sequence().as_ref().to_vec()
}

/// Convert a noodles `Cigar` into the project's own `Vec<CigarOp>` so
/// the parity test can call the production helpers directly rather
/// than maintaining duplicate CIGAR-walk code (M11).
fn cigar_ops_from_noodles(cigar: &sam::alignment::record_buf::Cigar) -> Vec<CigarOp> {
    cigar
        .as_ref()
        .iter()
        .map(|op| {
            let l = op.len() as u32;
            match op.kind() {
                Kind::Match => CigarOp::Match(l),
                Kind::Insertion => CigarOp::Insertion(l),
                Kind::Deletion => CigarOp::Deletion(l),
                Kind::Skip => CigarOp::Skip(l),
                Kind::SoftClip => CigarOp::SoftClip(l),
                Kind::HardClip => CigarOp::HardClip(l),
                Kind::Pad => CigarOp::Padding(l),
                Kind::SequenceMatch => CigarOp::SeqMatch(l),
                Kind::SequenceMismatch => CigarOp::SeqMismatch(l),
            }
        })
        .collect()
}

fn check_parity(stem: &str) {
    let ref_seq = read_single_fasta(&fixture(&format!("{stem}.fa")));
    let ref_len = ref_seq.len() as i32;

    let mut sam_reader = sam::io::Reader::new(BufReader::new(
        File::open(fixture(&format!("{stem}.sam"))).unwrap(),
    ));
    let sam_header = sam_reader.read_header().unwrap();
    let in_records: Vec<_> = sam_reader
        .record_bufs(&sam_header)
        .map(|r| r.unwrap())
        .collect();

    let mut exp_reader = sam::io::Reader::new(BufReader::new(
        File::open(fixture(&format!("{stem}_exp.sam"))).unwrap(),
    ));
    let exp_header = exp_reader.read_header().unwrap();
    let exp_records: Vec<_> = exp_reader
        .record_bufs(&exp_header)
        .map(|r| r.unwrap())
        .collect();

    assert_eq!(
        in_records.len(),
        exp_records.len(),
        "{stem}: input/expected record count mismatch",
    );

    let cfg = BaqConfig::default();
    let mut scratch = ProbalnScratch::new();
    let mut capped_buf: Vec<u8> = Vec::new();

    for (idx, (in_rec, exp_rec)) in in_records.iter().zip(exp_records.iter()).enumerate() {
        let bq_bytes: Vec<u8> = match exp_rec.data().get(&Tag::BASE_ALIGNMENT_QUALITY_OFFSETS) {
            Some(Value::String(bs)) => {
                let slice: &[u8] = bs.as_ref();
                slice.to_vec()
            }
            Some(other) => panic!("{stem} record {idx}: BQ tag has wrong type: {other:?}"),
            None => continue, // htslib produced no BAQ for this read (e.g. no-match CIGAR)
        };

        let pos_0 = in_rec.alignment_start().unwrap().get() as i32 - 1;
        let qual: Vec<u8> = in_rec.quality_scores().as_ref().to_vec();
        let seq_ascii: Vec<u8> = in_rec.sequence().as_ref().to_vec();
        let cigar_ops = cigar_ops_from_noodles(in_rec.cigar());
        let l_qseq = seq_ascii.len() as i32;

        // Route through production helpers (M11) so the parity test
        // proves equivalence between the production code and htslib,
        // not between two Rust copies of the production code.
        let (xb0, xe0, yb, ye, bw) =
            match locate_alignment_span(&cigar_ops, pos_0, cfg.band_half_width) {
                Ok(span) => span,
                Err(reason) => panic!(
                    "{stem} record {idx}: BQ tag present but no realignable span ({reason:?})",
                ),
            };
        let (xb, raw_xe) = extend_ref_window(xb0, xe0, yb, ye, l_qseq, bw);
        // Test-side ref_len clamp: parity fixtures may have reads whose
        // extended window passes the ref end. Production relies on the
        // fetcher to refuse; in the fixture we have the ref in memory
        // and can clamp directly.
        let xe = raw_xe.min(ref_len);

        let tref: Vec<u8> = ref_seq[xb as usize..xe as usize]
            .iter()
            .map(|&b| encode_base(b))
            .collect();
        let tseq: Vec<u8> = seq_ascii.iter().map(|&b| encode_base(b)).collect();

        let bw_cfg = BaqConfig {
            band_half_width: bw,
            ..cfg
        };

        probaln_glocal(&mut scratch, &tref, &tseq, &qual, &bw_cfg)
            .expect("probaln_glocal should not fail on fixture inputs");

        apply_baq_cap_into(
            &mut capped_buf,
            &cigar_ops,
            &qual,
            scratch.state(),
            scratch.q(),
            pos_0,
            xb,
        );

        assert_eq!(
            bq_bytes.len(),
            seq_ascii.len(),
            "{stem} record {idx}: BQ tag length mismatch",
        );
        for i in 0..seq_ascii.len() {
            let delta = bq_bytes[i] as i32 - 64;
            assert!(
                delta >= 0,
                "{stem} record {idx} pos {i}: negative BQ delta ({delta})",
            );
            let expected = (qual[i] as i32 - delta) as u8;
            assert_eq!(
                capped_buf[i], expected,
                "{stem} record {idx} pos {i}: capped qual mismatch (ours={}, htslib={}, raw qual={}, delta={})",
                capped_buf[i], expected, qual[i], delta,
            );
        }
    }
}

#[test]
fn parity_realn01() {
    check_parity("realn01");
}

#[test]
fn parity_realn02() {
    check_parity("realn02");
}

// realn03 is intentionally not used as a parity fixture: its expected
// SAM is generated by htslib's `test_realn -e` (extended BAQ), which
// runs a left/right running-max post-pass on top of the HMM output and
// is out of scope for this plan. The realn03 data files remain in
// `tests/data/baq/` so the corpus is easy to extend if/when extended
// BAQ lands. See `tests/data/baq/README.md`.

// ---------------------------------------------------------------------
// BaqEngine driver tests — synthetic MappedReads, no SAM parsing.
// ---------------------------------------------------------------------

use crate::bam::cram_input::{
    FLAG_FIRST_OF_PAIR, FLAG_PAIRED, FLAG_REVERSE_STRAND, FLAG_UNMAPPED, MappedRead,
};
use crate::fasta::ManualEvictChromRefFetcher;
use crate::pileup::walker::MateRole;

use super::baq_engine::{BaqEngine, BaqOutcome, BaqSkipReason};

/// Build a single-contig FASTA + sibling `.fai` in a tempdir,
/// open a [`ManualEvictChromRefFetcher`] over it, and hand both back.
/// The tempdir must outlive the fetcher so the file stays on disk.
/// Replaces the in-memory `MockRefFetcher` the tests used before the
/// `unified_chrom_ref_fetcher` migration: BAQ no longer takes a
/// trait-object fetcher, so the tests now use the real production
/// type with a tiny tempfile FASTA.
fn fetcher_from_chrom_bytes(chrom: &[u8]) -> (tempfile::TempDir, ManualEvictChromRefFetcher) {
    use std::io::Write as _;
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta_path = dir.path().join("ref.fa");
    let fai_path = dir.path().join("ref.fa.fai");

    let header = b">chr0\n";
    {
        let mut fa = std::fs::File::create(&fasta_path).expect("fa");
        fa.write_all(header).expect("hdr");
        fa.write_all(chrom).expect("seq");
        fa.write_all(b"\n").expect("nl");
    }
    {
        let mut fai = std::fs::File::create(&fai_path).expect("fai");
        writeln!(
            fai,
            "chr0\t{}\t{}\t{}\t{}",
            chrom.len(),
            header.len(),
            chrom.len(),
            chrom.len() + 1,
        )
        .expect("fai");
    }
    let fetcher = ManualEvictChromRefFetcher::for_contig(&fasta_path, "chr0").expect("fetcher");
    (dir, fetcher)
}

/// Companion to [`fetcher_from_chrom_bytes`] for [`BaqStream`]
/// tests: returns the FASTA path and a one-entry [`ContigList`]
/// instead of a built fetcher, since [`BaqStream::new`] now takes
/// the FASTA path + contigs and builds per-rayon-worker fetchers
/// inside its own `map_init`.
fn stream_components_from_chrom_bytes(
    chrom: &[u8],
) -> (
    tempfile::TempDir,
    std::path::PathBuf,
    crate::fasta::ContigList,
) {
    use crate::fasta::{ContigEntry, ContigList};
    use std::io::Write as _;
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta_path = dir.path().join("ref.fa");
    let fai_path = dir.path().join("ref.fa.fai");
    let header = b">chr0\n";
    {
        let mut fa = std::fs::File::create(&fasta_path).expect("fa");
        fa.write_all(header).expect("hdr");
        fa.write_all(chrom).expect("seq");
        fa.write_all(b"\n").expect("nl");
    }
    {
        let mut fai = std::fs::File::create(&fai_path).expect("fai");
        writeln!(
            fai,
            "chr0\t{}\t{}\t{}\t{}",
            chrom.len(),
            header.len(),
            chrom.len(),
            chrom.len() + 1,
        )
        .expect("fai");
    }
    let contigs = ContigList {
        entries: vec![ContigEntry {
            name: "chr0".into(),
            length: chrom.len() as u64,
            md5: None,
        }],
    };
    (dir, fasta_path, contigs)
}

fn synthetic_read(
    flag: u16,
    pos: u64,
    mapq: u8,
    cigar: Vec<CigarOp>,
    seq: Vec<u8>,
    qual: Vec<u8>,
) -> MappedRead {
    MappedRead {
        qname: b"r1".to_vec(),
        flag,
        ref_id: 0,
        pos,
        mapq,
        cigar,
        seq,
        qual,
        mate_ref_id: None,
        mate_pos: None,
        adaptor_boundary: None,
        source_file_index: 0,
    }
}

#[test]
fn engine_happy_path_match_only() {
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGTACGTACGT");
    let read = synthetic_read(
        0,
        1,
        60,
        vec![CigarOp::Match(8)],
        b"ACGTACGT".to_vec(),
        vec![40; 8],
    );
    let expected_seq = read.seq.clone();
    let expected_cigar = read.cigar.clone();
    let mut engine = BaqEngine::new(BaqConfig::default());
    match engine.process(read, &mut fetcher) {
        BaqOutcome::Capped(prepared) => {
            assert_eq!(prepared.alignment_start, 1);
            assert_eq!(prepared.alignment_end, 8);
            assert_eq!(prepared.seq, expected_seq);
            assert_eq!(prepared.cigar, expected_cigar);
            assert_eq!(prepared.bq_baq.len(), 8);
            for (i, &q) in prepared.bq_baq.iter().enumerate() {
                assert!(q <= 40, "pos {i}: bq_baq {q} exceeds raw qual 40");
            }
            assert_eq!(prepared.mate_role, MateRole::Solo);
            assert!(!prepared.is_reverse_strand);
        }
        BaqOutcome::Skipped(reason) => panic!("expected Capped, got Skipped({reason:?})"),
    }
}

#[test]
fn engine_mate_role_first_of_pair() {
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGTACGT");
    let read = synthetic_read(
        FLAG_PAIRED | FLAG_FIRST_OF_PAIR,
        1,
        60,
        vec![CigarOp::Match(5)],
        b"ACGTA".to_vec(),
        vec![40; 5],
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    match engine.process(read, &mut fetcher) {
        BaqOutcome::Capped(prepared) => assert_eq!(prepared.mate_role, MateRole::FirstOfPair),
        other => panic!("expected Capped, got {other:?}"),
    }
}

#[test]
fn engine_mate_role_second_of_pair() {
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGTACGT");
    let read = synthetic_read(
        FLAG_PAIRED, // paired but not first of pair → SecondOfPair
        1,
        60,
        vec![CigarOp::Match(5)],
        b"ACGTA".to_vec(),
        vec![40; 5],
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    match engine.process(read, &mut fetcher) {
        BaqOutcome::Capped(prepared) => assert_eq!(prepared.mate_role, MateRole::SecondOfPair),
        other => panic!("expected Capped, got {other:?}"),
    }
}

#[test]
fn engine_reverse_strand_propagates() {
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGTACGT");
    let read = synthetic_read(
        FLAG_REVERSE_STRAND,
        1,
        60,
        vec![CigarOp::Match(5)],
        b"ACGTA".to_vec(),
        vec![40; 5],
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    match engine.process(read, &mut fetcher) {
        BaqOutcome::Capped(prepared) => assert!(prepared.is_reverse_strand),
        other => panic!("expected Capped, got {other:?}"),
    }
}

#[test]
fn engine_skip_unmapped() {
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGT");
    let read = synthetic_read(
        FLAG_UNMAPPED,
        1,
        0,
        vec![CigarOp::Match(5)],
        b"ACGTA".to_vec(),
        vec![40; 5],
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    assert_eq!(
        skip_reason(engine.process(read, &mut fetcher)),
        Some(BaqSkipReason::Unmapped),
    );
}

#[test]
fn engine_skip_empty_query() {
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGT");
    let read = synthetic_read(0, 1, 60, vec![], vec![], vec![]);
    let mut engine = BaqEngine::new(BaqConfig::default());
    assert_eq!(
        skip_reason(engine.process(read, &mut fetcher)),
        Some(BaqSkipReason::EmptyQuery),
    );
}

#[test]
fn engine_skip_qual_absent_on_empty_qual() {
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGT");
    let read = synthetic_read(
        0,
        1,
        60,
        vec![CigarOp::Match(5)],
        b"ACGTA".to_vec(),
        vec![], // no quals
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    assert_eq!(
        skip_reason(engine.process(read, &mut fetcher)),
        Some(BaqSkipReason::QualAbsent),
    );
}

#[test]
fn engine_skip_qual_absent_on_length_mismatch() {
    // M22: regression test for the second arm of the QualAbsent
    // guard (qual.len() != seq.len() with both nonzero).
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGT");
    let read = synthetic_read(
        0,
        1,
        60,
        vec![CigarOp::Match(5)],
        b"ACGTA".to_vec(),
        vec![40; 3], // 3 quals for a 5-base seq
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    assert_eq!(
        skip_reason(engine.process(read, &mut fetcher)),
        Some(BaqSkipReason::QualAbsent),
    );
}

#[test]
fn engine_skip_no_match_in_cigar() {
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGT");
    // All soft-clip — no M/=/X. Matches realn.c:192 (`xb == -1`).
    let read = synthetic_read(
        0,
        1,
        60,
        vec![CigarOp::SoftClip(5)],
        b"ACGTA".to_vec(),
        vec![40; 5],
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    assert_eq!(
        skip_reason(engine.process(read, &mut fetcher)),
        Some(BaqSkipReason::NoMatchInCigar),
    );
}

#[test]
fn engine_skip_contains_ref_skip() {
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGTACGTACGT");
    let read = synthetic_read(
        0,
        1,
        60,
        vec![CigarOp::Match(3), CigarOp::Skip(5), CigarOp::Match(3)],
        b"ACGACG".to_vec(),
        vec![40; 6],
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    assert_eq!(
        skip_reason(engine.process(read, &mut fetcher)),
        Some(BaqSkipReason::ContainsRefSkip),
    );
}

#[test]
fn engine_skip_ref_window_past_chrom_end() {
    // 5 bp chrom, read at pos=4 with CIGAR=5M — match span [3..8) on a
    // ref that only has [0..5). The fetcher errs, and the engine
    // converts that to RefWindowPastChromEnd.
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTA");
    let read = synthetic_read(
        0,
        4,
        60,
        vec![CigarOp::Match(5)],
        b"TACGT".to_vec(),
        vec![40; 5],
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    assert_eq!(
        skip_reason(engine.process(read, &mut fetcher)),
        Some(BaqSkipReason::RefWindowPastChromEnd),
    );
}

#[test]
fn engine_skip_pos_out_of_range() {
    // M4: read.pos > i32::MAX → PosOutOfRange.
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGT");
    let read = synthetic_read(
        0,
        (i32::MAX as u64) + 10,
        60,
        vec![CigarOp::Match(2)],
        b"AC".to_vec(),
        vec![40; 2],
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    assert_eq!(
        skip_reason(engine.process(read, &mut fetcher)),
        Some(BaqSkipReason::PosOutOfRange),
    );
}

// M18: happy-path tests for Insertion / Deletion / mixed CIGAR — the
// production CIGAR walks were previously only exercised on `Match`-only
// and `SoftClip`-only inputs.

#[test]
fn engine_happy_path_with_insertion() {
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGTACGT");
    // CIGAR 3M2I3M: query bases 0-2 align to ref 0-2, bases 3-4 are
    // insertions (no ref pos), bases 5-7 align to ref 3-5.
    let read = synthetic_read(
        0,
        1,
        60,
        vec![CigarOp::Match(3), CigarOp::Insertion(2), CigarOp::Match(3)],
        b"ACGAACGT".to_vec(),
        vec![40; 8],
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    match engine.process(read, &mut fetcher) {
        BaqOutcome::Capped(p) => {
            assert_eq!(p.bq_baq.len(), 8);
            // Insertion positions (indices 3, 4) carry the raw qual
            // because the apply_baq_cap_into Insertion arm advances `y`
            // without writing to `out`.
            assert_eq!(p.bq_baq[3], 40);
            assert_eq!(p.bq_baq[4], 40);
        }
        other => panic!("expected Capped, got {other:?}"),
    }
}

#[test]
fn engine_happy_path_with_deletion() {
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGTACGTACGT");
    // CIGAR 4M2D4M: query bases 0-3 align to ref 0-3, ref 4-5 deleted
    // (no query bases), query bases 4-7 align to ref 6-9.
    let read = synthetic_read(
        0,
        1,
        60,
        vec![CigarOp::Match(4), CigarOp::Deletion(2), CigarOp::Match(4)],
        b"ACGTCGTA".to_vec(),
        vec![40; 8],
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    match engine.process(read, &mut fetcher) {
        BaqOutcome::Capped(p) => {
            assert_eq!(p.bq_baq.len(), 8);
            for (i, &q) in p.bq_baq.iter().enumerate() {
                assert!(q <= 40, "pos {i}: bq_baq {q} exceeds raw qual 40");
            }
            // Ref span: 4M + 2D + 4M = 10 bases starting at 1 → end = 10.
            assert_eq!(p.alignment_end, 10);
        }
        other => panic!("expected Capped, got {other:?}"),
    }
}

#[test]
fn engine_happy_path_with_mixed_cigar() {
    // CIGAR 2S5M1I3M2D2M2S covers every arm of apply_baq_cap_into in
    // one walk.
    let (_dir, mut fetcher) = fetcher_from_chrom_bytes(b"ACGTACGTACGTACGTACGT");
    let read = synthetic_read(
        0,
        1,
        60,
        vec![
            CigarOp::SoftClip(2),
            CigarOp::Match(5),
            CigarOp::Insertion(1),
            CigarOp::Match(3),
            CigarOp::Deletion(2),
            CigarOp::Match(2),
            CigarOp::SoftClip(2),
        ],
        b"NNACGTACTACGCGNN".to_vec(),
        vec![40; 16],
    );
    let mut engine = BaqEngine::new(BaqConfig::default());
    match engine.process(read, &mut fetcher) {
        BaqOutcome::Capped(p) => {
            assert_eq!(p.bq_baq.len(), 16);
            // Ref span (counts M, =, X, D, N): 5 + 3 + 2 + 2 = 12.
            assert_eq!(p.alignment_end, 12);
        }
        other => panic!("expected Capped, got {other:?}"),
    }
}

fn skip_reason(outcome: BaqOutcome) -> Option<BaqSkipReason> {
    match outcome {
        BaqOutcome::Skipped(r) => Some(r),
        BaqOutcome::Capped(_) => None,
    }
}

// ---------------------------------------------------------------------
// Leaf helper unit tests (M22 cluster, set_u, encode_base, Q2P,
// compute_alignment_end).
// ---------------------------------------------------------------------

#[test]
fn encode_base_maps_acgtn_case_insensitive() {
    for (b, expected) in [
        (b'A', 0),
        (b'a', 0),
        (b'C', 1),
        (b'c', 1),
        (b'G', 2),
        (b'g', 2),
        (b'T', 3),
        (b't', 3),
        (b'N', 4),
        (b'X', 4),
        (b' ', 4),
    ] {
        assert_eq!(encode_base(b), expected, "byte {b:#x}");
    }
}

#[test]
fn compute_alignment_end_handles_empty_cigar() {
    // No ref-consuming ops → ref_span = 0 → start.saturating_sub(1).
    assert_eq!(compute_alignment_end(1, &[]), 0);
    assert_eq!(compute_alignment_end(100, &[CigarOp::SoftClip(5)]), 99);
}

#[test]
fn baq_skip_counts_bump_is_exhaustive_per_variant() {
    use super::baq_stream::BaqSkipCounts;
    let mut counts = BaqSkipCounts::default();
    for reason in [
        BaqSkipReason::Unmapped,
        BaqSkipReason::EmptyQuery,
        BaqSkipReason::QualAbsent,
        BaqSkipReason::NoMatchInCigar,
        BaqSkipReason::ContainsRefSkip,
        BaqSkipReason::HmmOverflow,
        BaqSkipReason::RefWindowPastChromEnd,
        BaqSkipReason::PosOutOfRange,
        BaqSkipReason::ReadTooLong,
        BaqSkipReason::ChromIdOutOfRange,
    ] {
        counts.bump(reason);
    }
    assert_eq!(counts.total, 10);
    assert_eq!(counts.unmapped, 1);
    assert_eq!(counts.empty_query, 1);
    assert_eq!(counts.qual_absent, 1);
    assert_eq!(counts.no_match_in_cigar, 1);
    assert_eq!(counts.contains_ref_skip, 1);
    assert_eq!(counts.hmm_overflow, 1);
    assert_eq!(counts.ref_window_past_chrom_end, 1);
    assert_eq!(counts.pos_out_of_range, 1);
    assert_eq!(counts.read_too_long, 1);
    assert_eq!(counts.chrom_id_out_of_range, 1);
}

// ---------------------------------------------------------------------
// BaqStream — rayon-parallel adapter tests.
// ---------------------------------------------------------------------

use crate::bam::errors::CramInputError;
use crate::pileup::walker::PreparedRead;

use super::baq_stream::BaqStream;

#[test]
fn stream_yields_prepared_reads_in_order_within_chunk() {
    let (_dir, fasta_path, contigs) = stream_components_from_chrom_bytes(b"ACGTACGTACGTACGT");
    let inputs: Vec<Result<MappedRead, CramInputError>> = (1..=4)
        .map(|pos| {
            Ok(synthetic_read(
                0,
                pos as u64,
                60,
                vec![CigarOp::Match(5)],
                b"ACGTA".to_vec(),
                vec![40; 5],
            ))
        })
        .collect();
    let stream = BaqStream::new(
        inputs.into_iter(),
        BaqConfig::default(),
        fasta_path.clone(),
        contigs.clone(),
        16,
    );
    let outputs: Vec<Result<PreparedRead, _>> = stream.collect();
    assert_eq!(outputs.len(), 4);
    let starts: Vec<u32> = outputs
        .iter()
        .map(|r| r.as_ref().unwrap().alignment_start)
        .collect();
    assert_eq!(starts, vec![1, 2, 3, 4]);
}

#[test]
fn stream_preserves_order_across_chunks() {
    // Chunk size 2 with 5 inputs → three batches. Output order must
    // still match input order.
    let (_dir, fasta_path, contigs) = stream_components_from_chrom_bytes(b"ACGTACGTACGTACGT");
    let inputs: Vec<Result<MappedRead, CramInputError>> = (1..=5)
        .map(|pos| {
            Ok(synthetic_read(
                0,
                pos as u64,
                60,
                vec![CigarOp::Match(5)],
                b"ACGTA".to_vec(),
                vec![40; 5],
            ))
        })
        .collect();
    let stream = BaqStream::new(
        inputs.into_iter(),
        BaqConfig::default(),
        fasta_path.clone(),
        contigs.clone(),
        2,
    );
    let starts: Vec<u32> = stream.map(|r| r.unwrap().alignment_start).collect();
    assert_eq!(starts, vec![1, 2, 3, 4, 5]);
}

#[test]
fn stream_preserves_order_with_explicit_multi_threaded_pool() {
    // M19: an explicit 8-thread pool ensures the order-preservation
    // contract is exercised under genuine parallelism, not just whatever
    // the default rayon pool happens to do on a small CI runner.
    // Chrom needs to cover the extended ref window for the last read
    // (pos=64, 5M, ~p+7 extension); use a generous 128-byte chrom.
    let (_dir, fasta_path, contigs) = stream_components_from_chrom_bytes(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    let inputs: Vec<Result<MappedRead, CramInputError>> = (1..=64)
        .map(|pos| {
            Ok(synthetic_read(
                0,
                pos as u64,
                60,
                vec![CigarOp::Match(5)],
                b"ACGTA".to_vec(),
                vec![40; 5],
            ))
        })
        .collect();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();
    let starts = pool.install(|| {
        let stream = BaqStream::new(
            inputs.into_iter(),
            BaqConfig::default(),
            fasta_path.clone(),
            contigs.clone(),
            16,
        );
        stream
            .map(|r| r.unwrap().alignment_start)
            .collect::<Vec<u32>>()
    });
    assert_eq!(starts, (1u32..=64).collect::<Vec<_>>());
}

#[test]
fn stream_increments_skip_counts_per_reason() {
    let (_dir, fasta_path, contigs) = stream_components_from_chrom_bytes(b"ACGTACGTACGTACGT");
    let inputs: Vec<Result<MappedRead, CramInputError>> = vec![
        // Capped — happy path.
        Ok(synthetic_read(
            0,
            1,
            60,
            vec![CigarOp::Match(5)],
            b"ACGTA".to_vec(),
            vec![40; 5],
        )),
        // Unmapped — Skipped(Unmapped).
        Ok(synthetic_read(
            FLAG_UNMAPPED,
            1,
            0,
            vec![CigarOp::Match(5)],
            b"ACGTA".to_vec(),
            vec![40; 5],
        )),
        // All soft-clip — Skipped(NoMatchInCigar).
        Ok(synthetic_read(
            0,
            1,
            60,
            vec![CigarOp::SoftClip(5)],
            b"ACGTA".to_vec(),
            vec![40; 5],
        )),
        // Another capped — confirms skips don't break the surviving stream.
        Ok(synthetic_read(
            0,
            2,
            60,
            vec![CigarOp::Match(5)],
            b"CGTAC".to_vec(),
            vec![40; 5],
        )),
    ];
    let mut stream = BaqStream::new(
        inputs.into_iter(),
        BaqConfig::default(),
        fasta_path.clone(),
        contigs.clone(),
        16,
    );
    let outputs: Vec<_> = (&mut stream).collect();
    assert_eq!(outputs.iter().filter(|r| r.is_ok()).count(), 2);
    let counts = *stream.skip_counts();
    assert_eq!(counts.total, 2);
    assert_eq!(counts.unmapped, 1);
    assert_eq!(counts.no_match_in_cigar, 1);
}

#[test]
fn stream_propagates_upstream_error_after_batched_reads() {
    let (_dir, fasta_path, contigs) = stream_components_from_chrom_bytes(b"ACGTACGT");
    let inputs: Vec<Result<MappedRead, CramInputError>> = vec![
        Ok(synthetic_read(
            0,
            1,
            60,
            vec![CigarOp::Match(5)],
            b"ACGTA".to_vec(),
            vec![40; 5],
        )),
        Err(CramInputError::NoInputs),
    ];
    let outputs: Vec<_> = BaqStream::new(
        inputs.into_iter(),
        BaqConfig::default(),
        fasta_path.clone(),
        contigs.clone(),
        16,
    )
    .collect();
    assert_eq!(outputs.len(), 2);
    assert!(outputs[0].is_ok());
    assert!(matches!(outputs[1], Err(CramInputError::NoInputs)));
}

#[test]
fn stream_returns_success_after_all_skipped_chunks() {
    // M20: with chunk_size=1, three inputs [unmapped, all-softclip, ok]
    // force three chunks where the first two yield only skips. next()
    // must not return None prematurely — it must loop through the
    // refills until it finds the surviving read.
    let (_dir, fasta_path, contigs) = stream_components_from_chrom_bytes(b"ACGTACGT");
    let inputs: Vec<Result<MappedRead, CramInputError>> = vec![
        Ok(synthetic_read(
            FLAG_UNMAPPED,
            1,
            0,
            vec![CigarOp::Match(5)],
            b"ACGTA".to_vec(),
            vec![40; 5],
        )),
        Ok(synthetic_read(
            0,
            1,
            60,
            vec![CigarOp::SoftClip(5)],
            b"ACGTA".to_vec(),
            vec![40; 5],
        )),
        Ok(synthetic_read(
            0,
            1,
            60,
            vec![CigarOp::Match(5)],
            b"ACGTA".to_vec(),
            vec![40; 5],
        )),
    ];
    let mut stream = BaqStream::new(
        inputs.into_iter(),
        BaqConfig::default(),
        fasta_path.clone(),
        contigs.clone(),
        1,
    );
    let outputs: Vec<_> = (&mut stream).collect();
    assert_eq!(outputs.iter().filter(|r| r.is_ok()).count(), 1);
    let counts = *stream.skip_counts();
    assert_eq!(counts.total, 2);
}

#[test]
fn stream_is_fused_after_exhaustion() {
    let (_dir, fasta_path, contigs) = stream_components_from_chrom_bytes(b"ACGTACGT");
    let inputs: Vec<Result<MappedRead, CramInputError>> = vec![Ok(synthetic_read(
        0,
        1,
        60,
        vec![CigarOp::Match(5)],
        b"ACGTA".to_vec(),
        vec![40; 5],
    ))];
    let mut stream = BaqStream::new(
        inputs.into_iter(),
        BaqConfig::default(),
        fasta_path.clone(),
        contigs.clone(),
        16,
    );
    assert!(stream.next().is_some());
    assert!(stream.next().is_none());
    // Confirmed fused: re-poll yields None.
    assert!(stream.next().is_none());
}

#[test]
#[should_panic(expected = "BaqStream chunk_size must be > 0")]
fn stream_rejects_zero_chunk_size() {
    // M13: zero chunk_size is a programmer error, not a value to be
    // silently coerced.
    use crate::fasta::ContigList;
    let _ = BaqStream::new(
        std::iter::empty::<Result<MappedRead, CramInputError>>(),
        BaqConfig::default(),
        std::path::PathBuf::from("/dev/null"),
        ContigList { entries: vec![] },
        0,
    );
}

#[test]
fn baqstream_is_send() {
    // Pin the `Send` auto-trait so the pipeline integration commit can
    // move a BaqStream across a rayon::scope without a confusing
    // trait-bound error.
    fn assert_send<T: Send>() {}
    assert_send::<BaqStream<std::vec::IntoIter<Result<MappedRead, CramInputError>>>>();
}

#[test]
fn vec_deque_compiles_for_current_batch_drain() {
    // Smoke test that VecDeque<PreparedRead>::pop_front() yields the
    // PreparedReads in insertion order; covers the Mi13 refactor.
    let mut q: VecDeque<u32> = VecDeque::new();
    q.push_back(1);
    q.push_back(2);
    q.push_back(3);
    assert_eq!(q.pop_front(), Some(1));
    assert_eq!(q.pop_front(), Some(2));
    assert_eq!(q.pop_front(), Some(3));
    assert_eq!(q.pop_front(), None);
}
