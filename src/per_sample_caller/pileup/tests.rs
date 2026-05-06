//! End-to-end tests for the pileup walker plus shared helpers
//! (`MockFasta`, fixture builders) used by inner-module tests.
//!
//! Tests exercise scenarios from `ia/specs/pileup_walker.md`: SNP
//! folding, deletion anchoring, REF widening on overlap, mate
//! overlap, phase-chain lifecycle markers, eager closure, etc.

use std::sync::Arc;
use std::sync::mpsc;
use std::thread;

use super::CigarOp;
use super::PreparedRead;
use super::RefBaseFetcher;
use super::run;

// ---------------------------------------------------------------------
// MockFasta
// ---------------------------------------------------------------------

/// In-memory `RefBaseFetcher` implementation backed by a single
/// chromosome string. Tests inject their reference here directly
/// instead of building a real FASTA file.
#[derive(Debug, Clone)]
pub struct MockFasta {
    /// Reference bases, indexed 0-based for storage. The walker
    /// uses 1-based coordinates externally — `fetch` translates.
    chromosomes: Vec<Vec<u8>>,
}

impl MockFasta {
    pub fn new(chr0: &str) -> Self {
        Self {
            chromosomes: vec![chr0.as_bytes().to_vec()],
        }
    }
}

impl RefBaseFetcher for MockFasta {
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, std::io::Error> {
        let chrom = self
            .chromosomes
            .get(chrom_id as usize)
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::NotFound, "chrom"))?;
        let start_idx = (start_1based - 1) as usize;
        let end_idx = start_idx + length as usize;
        if end_idx > chrom.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                format!(
                    "fetch [{}, {}) past chrom len {}",
                    start_1based,
                    start_1based + length,
                    chrom.len()
                ),
            ));
        }
        Ok(chrom[start_idx..end_idx].to_vec())
    }
}

// ---------------------------------------------------------------------
// Fixture helpers
// ---------------------------------------------------------------------

pub fn snp_read(qname: &str, alignment_start: u32, seq: &[u8], qual: &[u8]) -> PreparedRead {
    let len = seq.len() as u32;
    PreparedRead {
        chrom_id: 0,
        alignment_start,
        alignment_end: alignment_start + len - 1,
        cigar: vec![CigarOp::Match(len)],
        seq: seq.to_vec(),
        bq_baq: qual.to_vec(),
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from(qname),
        is_first_mate: true,
        has_mate: false,
    }
}

pub fn paired_snp_reads(
    qname: &str,
    alignment_start_a: u32,
    alignment_start_b: u32,
    seq: &[u8],
    qual: &[u8],
) -> (PreparedRead, PreparedRead) {
    let mut a = snp_read(qname, alignment_start_a, seq, qual);
    a.has_mate = true;
    a.is_first_mate = true;
    let mut b = snp_read(qname, alignment_start_b, seq, qual);
    b.has_mate = true;
    b.is_first_mate = false;
    (a, b)
}

/// Drive `run` on a fixed input list, collecting emitted records.
/// Stage 2 runs on a separate thread per spec; we mirror that
/// shape here to exercise the channel emission path properly.
pub fn drive_walker(reads: Vec<PreparedRead>, fasta: MockFasta) -> Vec<super::PileupRecord> {
    let (tx, rx) = mpsc::sync_channel::<super::PileupRecord>(64);
    let collector = thread::spawn(move || {
        let mut out = Vec::new();
        while let Ok(r) = rx.recv() {
            out.push(r);
        }
        out
    });
    let summary = run(reads, &fasta, &tx).expect("walker run failed");
    drop(tx);
    let _ = summary;
    collector.join().expect("collector thread panicked")
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[test]
fn pure_ref_pileup_emits_one_record_per_position_with_only_ref_allele() {
    // Reference: ACGTA at positions 1..5
    // Two reads, each spanning all 5 positions, all REF.
    let fa = MockFasta::new("ACGTA");
    let r1 = snp_read("r1", 1, b"ACGTA", &[30; 5]);
    let r2 = snp_read("r2", 1, b"ACGTA", &[30; 5]);
    let records = drive_walker(vec![r1, r2], fa);
    assert_eq!(records.len(), 5);
    for (i, rec) in records.iter().enumerate() {
        assert_eq!(rec.pos, (i + 1) as u32);
        assert_eq!(rec.alleles.len(), 1, "REF only at clean position");
        assert_eq!(rec.alleles[0].scalars.num_obs, 2);
    }
}

#[test]
fn snp_at_one_position_emits_record_with_two_alleles() {
    // Reference ACGTA. r1 ref-everywhere, r2 has SNP G→T at pos 3.
    let fa = MockFasta::new("ACGTA");
    let r1 = snp_read("r1", 1, b"ACGTA", &[30; 5]);
    let r2 = snp_read("r2", 1, b"ACTTA", &[30; 5]);
    let records = drive_walker(vec![r1, r2], fa);
    assert_eq!(records.len(), 5);
    let rec_pos3 = &records[2];
    assert_eq!(rec_pos3.pos, 3);
    assert_eq!(rec_pos3.alleles.len(), 2, "REF + SNP");
    // First is REF (G), supported by 1 read.
    assert_eq!(rec_pos3.alleles[0].seq, b"G");
    assert_eq!(rec_pos3.alleles[0].scalars.num_obs, 1);
    // Second is SNP (T), supported by 1 read.
    assert_eq!(rec_pos3.alleles[1].seq, b"T");
    assert_eq!(rec_pos3.alleles[1].scalars.num_obs, 1);
}

#[test]
fn deletion_record_has_extended_ref_span() {
    // Reference AAAATTTTGG. One read with CIGAR 4M3D3M starting at 1:
    //   AAAA at 1..4, deletion at 5..7, TGG at 8..10.
    // Expected: anchor record at 4 (one-before deletion start),
    // REF "ATTT" (4 bases), DEL allele "A" (anchor only).
    let fa = MockFasta::new("AAAATTTGG");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 9,
        cigar: vec![CigarOp::Match(4), CigarOp::Deletion(3), CigarOp::Match(2)],
        seq: b"AAAAGG".to_vec(),
        bq_baq: vec![30; 6],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r1"),
        is_first_mate: true,
        has_mate: false,
    };
    let records = drive_walker(vec![r], fa);
    let anchor = records
        .iter()
        .find(|r| r.pos == 4)
        .expect("must emit anchor at deletion's preceding base");
    assert_eq!(anchor.ref_span(), 4, "anchor + 3 deleted = 4");
    assert_eq!(anchor.alleles[0].seq, b"ATTT", "REF over the deletion span");
    let del = anchor
        .alleles
        .iter()
        .find(|a| a.seq.as_slice() == b"A")
        .expect("DEL allele = anchor only");
    assert_eq!(del.scalars.num_obs, 1);
}

#[test]
fn deletion_record_does_not_double_count_ref_reads() {
    // Reference: ACGTAC (positions 1..6).
    // r1: pure-Match across 1..5 (5M). All REF.
    // r2: 1M3D1M starting at 2 → Match at 2, Deletion of 3 anchored
    //     at 2, Match at 6.
    //
    // The deletion record at pos=2 widens to span 4 (footprint
    // [2, 6)). r1 spans the whole record's footprint with REF
    // bases, so REF.num_obs should be 1 (one ref-spanning read);
    // before the fix it was 4 because r1 was re-folded once per
    // walker step inside the open record's footprint, multiplying
    // every five-scalar value by ref_span.
    let fa = MockFasta::new("ACGTAC");
    let r1 = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 5,
        cigar: vec![CigarOp::Match(5)],
        seq: b"ACGTA".to_vec(),
        bq_baq: vec![30; 5],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r1"),
        is_first_mate: true,
        has_mate: false,
    };
    let r2 = PreparedRead {
        chrom_id: 0,
        alignment_start: 2,
        alignment_end: 6,
        cigar: vec![CigarOp::Match(1), CigarOp::Deletion(3), CigarOp::Match(1)],
        seq: b"CC".to_vec(),
        bq_baq: vec![30; 2],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r2"),
        is_first_mate: true,
        has_mate: false,
    };
    let records = drive_walker(vec![r1, r2], fa);
    let anchor = records
        .iter()
        .find(|r| r.pos == 2)
        .expect("anchor at deletion's preceding base");
    assert_eq!(anchor.ref_span(), 4, "anchor + 3 deleted = 4");
    let ref_allele = &anchor.alleles[0];
    assert_eq!(
        ref_allele.scalars.num_obs, 1,
        "REF: 1 obs from r1 only; got {}",
        ref_allele.scalars.num_obs
    );
    assert_eq!(ref_allele.scalars.fwd, 1, "REF: forward strand count = 1");
    let del = anchor
        .alleles
        .iter()
        .find(|a| a.seq.as_slice() == b"C")
        .expect("DEL allele = anchor base only");
    assert_eq!(del.scalars.num_obs, 1, "DEL: 1 obs from r2");
}

#[test]
fn insertion_record_has_alt_longer_than_ref() {
    // Reference AAAACGT. Read 1M2I5M = 1 M at pos 1 ("A"), 2-base
    // insertion ("XX"), 5 M from pos 2.
    let fa = MockFasta::new("AAAACGT");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 6,
        cigar: vec![CigarOp::Match(1), CigarOp::Insertion(2), CigarOp::Match(5)],
        seq: b"AXXAAACG".to_vec(),
        bq_baq: vec![30; 8],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r1"),
        is_first_mate: true,
        has_mate: false,
    };
    let records = drive_walker(vec![r], fa);
    let anchor = records.iter().find(|r| r.pos == 1).expect("anchor at 1");
    let ins = anchor
        .alleles
        .iter()
        .find(|a| a.seq.len() > anchor.ref_span() as usize);
    assert!(ins.is_some(), "INS allele should be longer than REF");
    let ins = ins.unwrap();
    assert_eq!(ins.seq, b"AXX", "anchor + 2 inserted bases");
    assert_eq!(ins.scalars.num_obs, 1);
}

#[test]
fn forward_strand_count_recorded_correctly() {
    // Reference ACG. Two reads: forward, reverse. Both REF.
    let fa = MockFasta::new("ACG");
    let mut r1 = snp_read("r1", 1, b"ACG", &[30; 3]);
    r1.is_reverse_strand = false;
    let mut r2 = snp_read("r2", 1, b"ACG", &[30; 3]);
    r2.is_reverse_strand = true;
    let records = drive_walker(vec![r1, r2], fa);
    let rec = &records[0];
    assert_eq!(rec.alleles[0].scalars.num_obs, 2);
    assert_eq!(rec.alleles[0].scalars.fwd, 1);
}

#[test]
fn placed_left_and_placed_start_are_per_record() {
    // Reference ACGTA. Two reads:
    //   r1 starts at pos 1, covers 1..5
    //   r2 starts at pos 3, covers 3..5
    // At record pos 3:
    //   r1 was placed_left (start=1 < 3)
    //   r2 was placed_start (start=3 == 3)
    let fa = MockFasta::new("ACGTA");
    let r1 = snp_read("r1", 1, b"ACGTA", &[30; 5]);
    let r2 = snp_read("r2", 3, b"GTA", &[30; 3]);
    let records = drive_walker(vec![r1, r2], fa);
    let rec3 = records.iter().find(|r| r.pos == 3).unwrap();
    assert_eq!(rec3.alleles[0].scalars.num_obs, 2);
    assert_eq!(rec3.alleles[0].scalars.placed_left, 1);
    assert_eq!(rec3.alleles[0].scalars.placed_start, 1);
}

#[test]
fn uncovered_positions_produce_no_records() {
    // Reference ACGTACGTAC (10 bp). Reads at pos 1..3 and pos 7..9.
    // Positions 4..6 should produce no records.
    let fa = MockFasta::new("ACGTACGTAC");
    let r1 = snp_read("r1", 1, b"ACG", &[30; 3]);
    let r2 = snp_read("r2", 7, b"GTA", &[30; 3]);
    let records = drive_walker(vec![r1, r2], fa);
    let positions: Vec<u32> = records.iter().map(|r| r.pos).collect();
    assert_eq!(positions, vec![1, 2, 3, 7, 8, 9]);
}

#[test]
fn paired_mates_share_chain_slot_id() {
    // Two paired reads at non-overlapping positions, same QNAME.
    // Both mates' contributions to their respective records should
    // carry the same chain_slot_id. Use a reference of all 'A's
    // and read sequences of all 'A' so every fold lands in the REF
    // bucket — keeps the test focused on chain_slot wiring without
    // entangling the allele-identification path.
    let fa = MockFasta::new("AAAAAAAAAAAAAAAAAAAA");
    let (m1, m2) = paired_snp_reads("pair", 1, 10, b"AAA", &[30; 3]);
    let records = drive_walker(vec![m1, m2], fa);
    let rec1 = records.iter().find(|r| r.pos == 1).unwrap();
    let rec10 = records.iter().find(|r| r.pos == 10).unwrap();
    // REF bucket on each side carries the same shared slot id.
    assert_eq!(rec1.alleles[0].chain_slots, vec![0u16]);
    assert_eq!(rec10.alleles[0].chain_slots, vec![0u16]);
}

#[test]
fn mate_overlap_zeroes_lower_bq_contribution() {
    // Two paired mates overlapping at the same position with
    // different BAQ-capped BQs. Higher BQ wins; lower contributes
    // to obs count but ln_bq becomes ln(1) = 0 — so q_sum gets
    // only one meaningful contribution.
    //
    // Mate 1 at pos 1, BQ=30, length 3.
    // Mate 2 at pos 1, BQ=10, length 3 (overlapping).
    let fa = MockFasta::new("ACG");
    let (m1, mut m2) = paired_snp_reads("p", 1, 1, b"ACG", &[30; 3]);
    m2.bq_baq = vec![10; 3];
    let records = drive_walker(vec![m1, m2], fa);
    let rec = &records[0];
    assert_eq!(rec.alleles[0].scalars.num_obs, 2, "both mates count");
    // q_sum: max(ln_BQ, ln_MQ).
    // Mate 1 (kept): max(ln(P_err Q=30), ln(P_err MQ ≈ -3)) ≈ -3
    // Mate 2 (zeroed bq): max(ln(1)=0, ln(P_err MQ ≈ -3)) = 0
    // Sum ≈ -3.
    assert!(
        rec.alleles[0].scalars.q_sum > -4.0 && rec.alleles[0].scalars.q_sum < -2.0,
        "expected q_sum ≈ -3, got {}",
        rec.alleles[0].scalars.q_sum
    );
}

#[test]
fn paired_mate_indel_overlap_yields_single_observation() {
    // Both mates of a pair report the same insertion at the same
    // anchor. Per spec §"Mate overlap on indels": treat as one
    // observation, not two — assign the higher-BQ-proxy event to
    // the bucket and drop the other. The previous walker zeroed
    // the loser's BQ but still folded it as a separate
    // observation.
    //
    // Reference AAAACGT (length 7). Both mates' CIGAR is 1M2I5M
    // → anchor Match at 1, Insertion of 2 bp ("XX") at anchor 1,
    // five Matches over positions 2..6.
    let fa = MockFasta::new("AAAACGT");
    let cigar = vec![CigarOp::Match(1), CigarOp::Insertion(2), CigarOp::Match(5)];
    let mate_a = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 6,
        cigar: cigar.clone(),
        seq: b"AXXAAACG".to_vec(),
        bq_baq: vec![30; 8],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("p"),
        is_first_mate: true,
        has_mate: true,
    };
    let mut mate_b = mate_a.clone();
    mate_b.is_first_mate = false;
    mate_b.bq_baq = vec![20; 8]; // lower BQ → loser

    let records = drive_walker(vec![mate_a, mate_b], fa);
    let anchor = records
        .iter()
        .find(|r| r.pos == 1)
        .expect("anchor record at pos 1");
    let ins = anchor
        .alleles
        .iter()
        .find(|a| a.seq.len() > anchor.ref_span() as usize)
        .expect("INS allele present");
    assert_eq!(
        ins.scalars.num_obs, 1,
        "indel-overlap collapses to one observation; got {}",
        ins.scalars.num_obs
    );
    // Forward-strand count should also reflect a single observation.
    assert_eq!(ins.scalars.fwd, 1);
}

#[test]
fn record_emits_in_coordinate_order_across_reads() {
    // 100 reads at increasing starts. We just want to verify the
    // emitted record stream is monotonically ordered by pos.
    let mut chrom = String::with_capacity(200);
    for _ in 0..40 {
        chrom.push_str("ACGTA");
    }
    let fa = MockFasta::new(&chrom);
    let mut reads = Vec::new();
    for i in 0..50u32 {
        let start = i * 2 + 1;
        reads.push(snp_read(&format!("r{i}"), start, b"AC", &[30; 2]));
    }
    let records = drive_walker(reads, fa);
    // Emitted records' positions must be sorted ascending.
    for w in records.windows(2) {
        assert!(w[0].pos <= w[1].pos, "out-of-order emission");
    }
}

#[test]
fn out_of_order_input_is_a_hard_error() {
    let fa = MockFasta::new("ACGTACGTAC");
    let r1 = snp_read("r1", 5, b"ACG", &[30; 3]);
    let r2 = snp_read("r2", 1, b"ACG", &[30; 3]); // before r1 — invalid
    let (tx, _rx) = mpsc::sync_channel::<super::PileupRecord>(64);
    let result = run(vec![r1, r2], &fa, &tx);
    assert!(result.is_err());
    let err = result.unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("out-of-order"), "got: {msg}");
}

#[test]
fn lifecycle_markers_appear_on_emitted_records() {
    // Single read covering 1..3. The slot is allocated when r
    // enters and released when r exits. Across the 3 emitted
    // records we should see exactly one `new_chains` entry (on
    // the first record) and exactly one `expired_chains` entry
    // (on a later record after r has exited).
    let fa = MockFasta::new("ACG");
    let r = snp_read("r", 1, b"ACG", &[30; 3]);
    let records = drive_walker(vec![r], fa);
    let total_new: usize = records.iter().map(|r| r.new_chains.len()).sum();
    let total_expired: usize = records.iter().map(|r| r.expired_chains.len()).sum();
    assert_eq!(total_new, 1, "exactly one chain start across all records");
    assert_eq!(total_expired, 1, "exactly one chain end across all records");
}
