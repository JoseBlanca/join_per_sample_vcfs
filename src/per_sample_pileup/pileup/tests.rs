//! End-to-end tests for the pileup walker plus shared helpers
//! (`MockFasta`, fixture builders) used by inner-module tests.
//!
//! Tests exercise scenarios from `ia/specs/pileup_walker.md`: SNP
//! folding, deletion anchoring, REF widening on overlap, mate
//! overlap, phase-chain lifecycle markers, eager closure, etc.

use std::sync::Arc;

use super::CigarOp;
use super::MateRole;
use super::MultiChromRefFetcher;
use super::PreparedRead;
use super::WalkerConfig;
use super::run;
use crate::per_sample_pileup::ref_fetcher::ChromRefFetchError;

// ---------------------------------------------------------------------
// MockFasta
// ---------------------------------------------------------------------

/// In-memory `MultiChromRefFetcher` implementation backed by a single
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

    /// Multi-chromosome variant: each entry is the literal bases of
    /// chromosome `i` (`chrom_id == i`). Used by tests that exercise
    /// chromosome-boundary behaviour.
    pub fn with_chromosomes(chroms: &[&str]) -> Self {
        Self {
            chromosomes: chroms.iter().map(|s| s.as_bytes().to_vec()).collect(),
        }
    }
}

impl MultiChromRefFetcher for MockFasta {
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, ChromRefFetchError> {
        let chrom_name = format!("chrom_id={chrom_id}");
        let chrom =
            self.chromosomes
                .get(chrom_id as usize)
                .ok_or_else(|| ChromRefFetchError::Io {
                    chrom_name: chrom_name.clone(),
                    source: std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("unknown chrom_id={chrom_id}"),
                    ),
                })?;
        if start_1based == 0 {
            return Err(ChromRefFetchError::InvalidStart);
        }
        let end_exclusive = start_1based + length;
        if (end_exclusive - 1) as usize > chrom.len() {
            return Err(ChromRefFetchError::OutOfBounds {
                chrom_name,
                chrom_length: chrom.len() as u32,
                start: start_1based,
                end: end_exclusive,
            });
        }
        let start_idx = (start_1based - 1) as usize;
        let end_idx = start_idx + length as usize;
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
        mapq: 60,
        is_reverse_strand: false,
        qname: Arc::from(qname),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,
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
    a.mate_role = MateRole::FirstOfPair;
    let mut b = snp_read(qname, alignment_start_b, seq, qual);
    b.mate_role = MateRole::SecondOfPair;
    (a, b)
}

/// Drive `run` on a fixed input list, collecting emitted records.
pub fn drive_walker(
    reads: Vec<PreparedRead>,
    ref_fetcher: MockFasta,
) -> Vec<crate::pileup_record::PileupRecord> {
    drive_walker_with_summary(reads, ref_fetcher).0
}

/// Same as `drive_walker` but also returns the run's
/// `RunSummary`. Useful for tests that assert on counters.
pub fn drive_walker_with_summary(
    reads: Vec<PreparedRead>,
    ref_fetcher: MockFasta,
) -> (
    Vec<crate::pileup_record::PileupRecord>,
    super::walker::RunSummary,
) {
    drive_walker_with_config(reads, ref_fetcher, &WalkerConfig::default())
}

/// Drive `run` with an explicit `WalkerConfig`. Used by tests that
/// need to override defaults (e.g., the per-column depth caps to
/// trip the truncation path with a small synthetic input).
pub fn drive_walker_with_config(
    reads: Vec<PreparedRead>,
    ref_fetcher: MockFasta,
    config: &WalkerConfig,
) -> (
    Vec<crate::pileup_record::PileupRecord>,
    super::walker::RunSummary,
) {
    let mut walker = run(reads, &ref_fetcher, config);
    let records: Vec<crate::pileup_record::PileupRecord> = (&mut walker)
        .map(|r| r.expect("walker yielded error"))
        .collect();
    let summary = walker.summary();
    (records, summary)
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
        assert_eq!(rec.alleles[0].support.num_obs, 2);
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
    assert_eq!(rec_pos3.alleles[0].support.num_obs, 1);
    // Second is SNP (T), supported by 1 read.
    assert_eq!(rec_pos3.alleles[1].seq, b"T");
    assert_eq!(rec_pos3.alleles[1].support.num_obs, 1);
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
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
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
    assert_eq!(del.support.num_obs, 1);
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
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
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
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    let records = drive_walker(vec![r1, r2], fa);
    let anchor = records
        .iter()
        .find(|r| r.pos == 2)
        .expect("anchor at deletion's preceding base");
    assert_eq!(anchor.ref_span(), 4, "anchor + 3 deleted = 4");
    let ref_allele = &anchor.alleles[0];
    assert_eq!(
        ref_allele.support.num_obs, 1,
        "REF: 1 obs from r1 only; got {}",
        ref_allele.support.num_obs
    );
    assert_eq!(ref_allele.support.fwd, 1, "REF: forward strand count = 1");
    let del = anchor
        .alleles
        .iter()
        .find(|a| a.seq.as_slice() == b"C")
        .expect("DEL allele = anchor base only");
    assert_eq!(del.support.num_obs, 1, "DEL: 1 obs from r2");
}

#[test]
fn refold_after_widen_clears_chain_id_from_old_bucket() {
    // Reference ACGTAC (positions 1..6).
    //
    // R0 (1M2D1M @ pos 1): Match(A @ 1), Del(CG @ 2..3),
    //     Match(T @ 4). seq "AT". Opens record at pos 1
    //     directly at span 3 (footprint [1, 4)).
    // R1 (5M @ pos 1) with a T→C SNP at pos 4: seq "ACGCA".
    //     At walker_pos 1, R1 folds into the record at pos 1
    //     with events overlapping [1, 4): three Matches at
    //     pos 1..3, all REF — R1 lands in the REF bucket "ACG".
    // R3 (1M1D1M @ pos 3): Match(G @ 3), Del(T @ 4),
    //     Match(A @ 5). seq "GA". At walker_pos 3, R3's
    //     deletion has footprint [3, 5), which overlaps the
    //     existing record at pos 1 (footprint [1, 4)) and
    //     extends past it — widens the record to span 4
    //     (footprint [1, 5)).
    //
    // At walker_pos 3, R1 is a contributor (Match @ 3) and
    // re-folds under the now-wider REF "ACGT": its event
    // window now includes Match(4) = C (SNP), so the new
    // allele seq is "ACGC" — a *different* bucket from REF
    // "ACGT" (the REF bucket's seq was extended in-place by
    // `widen()` from "ACG" to "ACGT").
    //
    // Walker invariant under test: after R1 re-folds out of
    // the REF bucket into "ACGC", the REF bucket must not
    // retain R1's chain id. With the open-record bug,
    // `subtract_contribution` zeroes the scalars but never
    // removes the read's chain id from the old
    // `OpenAllele.chain_ids` vector — so REF ends up with
    // num_obs = 0 and chain_ids = [1] (R1's stale id).
    // Downstream, the merger sees a chain-anchored
    // constituent with no observations and (rightly) refuses
    // to compute a quality from it.
    let fa = MockFasta::new("ACGTAC");
    let r0 = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 4,
        cigar: vec![CigarOp::Match(1), CigarOp::Deletion(2), CigarOp::Match(1)],
        seq: b"AT".to_vec(),
        bq_baq: vec![30; 2],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r0"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    let r1 = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 5,
        cigar: vec![CigarOp::Match(5)],
        seq: b"ACGCA".to_vec(),
        bq_baq: vec![30; 5],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r1"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    let r3 = PreparedRead {
        chrom_id: 0,
        alignment_start: 3,
        alignment_end: 5,
        cigar: vec![CigarOp::Match(1), CigarOp::Deletion(1), CigarOp::Match(1)],
        seq: b"GA".to_vec(),
        bq_baq: vec![30; 2],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r3"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    let records = drive_walker(vec![r0, r1, r3], fa);
    let anchor = records
        .iter()
        .find(|r| r.pos == 1)
        .expect("anchor record at pos 1");

    // Universal invariant: in every emitted record, every
    // allele bucket's `chain_ids.len()` must be `<= num_obs`.
    // Each chain id represents at least one observation that
    // landed in this bucket; a chain id with no backing
    // observation is a leftover from a re-fold that the
    // walker forgot to clean up.
    for allele in &anchor.alleles {
        assert!(
            allele.chain_ids.len() <= allele.support.num_obs as usize,
            "allele {:?} has chain_ids={:?} but num_obs={} — \
             stale chain ids exceed observations",
            std::str::from_utf8(&allele.seq).unwrap_or("<non-utf8>"),
            allele.chain_ids,
            allele.support.num_obs,
        );
    }
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
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
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
    assert_eq!(ins.support.num_obs, 1);
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
    assert_eq!(rec.alleles[0].support.num_obs, 2);
    assert_eq!(rec.alleles[0].support.fwd, 1);
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
    assert_eq!(rec3.alleles[0].support.num_obs, 2);
    assert_eq!(rec3.alleles[0].support.placed_left, 1);
    assert_eq!(rec3.alleles[0].support.placed_start, 1);
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
fn paired_mates_with_overlapping_positions_share_chain_id() {
    // Both mates start at the same position; the second admits
    // while the first is still in the active set. They collapse
    // onto the same chain id via the `pending_mates` hash lookup.
    let fa = MockFasta::new("AAAAA");
    let (m1, m2) = paired_snp_reads("pair", 1, 1, b"AAA", &[30; 3]);
    let records = drive_walker(vec![m1, m2], fa);
    let rec1 = records.iter().find(|r| r.pos == 1).unwrap();
    assert_eq!(rec1.alleles[0].chain_ids, vec![0u64]);
}

#[test]
fn paired_mates_within_lookup_window_share_chain_id_across_active_set_exit() {
    // The first mate exits the active set well before the second
    // mate admits. Mate-pair tracking is governed by
    // `mate_lookup_window`, not by active-set residence — the
    // `pending_mates` entry stays alive across the first mate's
    // exit, and the second mate's later arrival (still within the
    // window) reuses the first mate's chain id.
    //
    // Fixture: m1 covers pos 1-3, m2 admits at pos 10. m1's exit
    // at walker_pos=4 must *not* drop the pending entry; the
    // 10 → 1 = 9 bp separation is well inside the default
    // `mate_lookup_window` of 10 000 bp.
    let fa = MockFasta::new("AAAAAAAAAAAAAAAAAAAA");
    let (m1, m2) = paired_snp_reads("pair", 1, 10, b"AAA", &[30; 3]);
    let records = drive_walker(vec![m1, m2], fa);
    let rec1 = records.iter().find(|r| r.pos == 1).unwrap();
    let rec10 = records.iter().find(|r| r.pos == 10).unwrap();
    assert_eq!(rec1.alleles[0].chain_ids, vec![0u64]);
    assert_eq!(
        rec10.alleles[0].chain_ids,
        vec![0u64],
        "the two mates of a single pair must share one chain id"
    );
}

#[test]
fn paired_mates_separated_beyond_lookup_window_get_distinct_chain_ids() {
    // When the second mate arrives more than `mate_lookup_window`
    // bp past the first mate's `alignment_start`, the chain-id
    // allocator's `evict_stale_pending` walk has dropped the
    // pending entry by then. The second mate cannot match, mints
    // a fresh chain id, and is therefore treated as a separate
    // molecule. This is correct: when the pair is that far apart,
    // we have no read-level evidence that they're the same
    // physical fragment within our trustworthy-pairing window.
    //
    // Fixture: m1 at pos 1, m2 at pos 12_001. With default
    // `mate_lookup_window = 10_000`, m2 is past the eviction
    // threshold (1 + 10_000 + 1 = 10 002 first hits it).
    let n = 12_010_usize;
    let fa = MockFasta::new(&"A".repeat(n));
    let (m1, m2) = paired_snp_reads("pair", 1, 12_001, b"AAA", &[30; 3]);
    let records = drive_walker(vec![m1, m2], fa);
    let rec_a = records.iter().find(|r| r.pos == 1).unwrap();
    let rec_b = records.iter().find(|r| r.pos == 12_001).unwrap();
    assert_eq!(rec_a.alleles[0].chain_ids, vec![0u64]);
    assert_eq!(
        rec_b.alleles[0].chain_ids,
        vec![1u64],
        "beyond the lookup window the pair-tracking entry has been evicted; \
         the second mate gets a fresh id"
    );
}

#[test]
fn mate_overlap_bq_tie_prefers_first_mate_not_earlier_position() {
    // Two paired mates overlapping at the same anchor with the
    // SAME alignment_start and the SAME BAQ-capped BQ. Distinguish
    // them by mq_log_err so the kept-mate's contribution to q_sum
    // is identifiable. Per spec §"Tie-breaking and disagreement":
    // BQ tie → prefer mate 1 (the read whose `mate_role` is
    // `FirstOfPair`).
    //
    // The previous tie-break used `alignment_start` and dropped
    // contributor `b` arbitrarily on equal starts, so the
    // first-of-pair distinction was ignored.
    let fa = MockFasta::new("ACG");
    let m_first = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 3,
        cigar: vec![CigarOp::Match(3)],
        seq: b"ACG".to_vec(),
        bq_baq: vec![30; 3],
        mq_log_err: -2.0, // distinct mq_log_err for the kept mate
        is_reverse_strand: false,
        qname: Arc::from("p"),
        mate_role: MateRole::FirstOfPair,
        adaptor_boundary: None,

        mapq: 60,
    };
    let m_second = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 3,
        cigar: vec![CigarOp::Match(3)],
        seq: b"ACG".to_vec(),
        bq_baq: vec![30; 3],
        mq_log_err: -10.0, // distinct mq_log_err for the loser
        is_reverse_strand: false,
        qname: Arc::from("p"),
        mate_role: MateRole::SecondOfPair,
        adaptor_boundary: None,

        mapq: 60,
    };
    // First-mate appears AFTER the second mate in the input stream,
    // even though both have alignment_start = 1, to make sure the
    // tie-break is decided by `MateRole::FirstOfPair` rather than
    // by stream order or alignment_start.
    let records = drive_walker(vec![m_second, m_first], fa);
    let rec = &records[0];
    assert_eq!(rec.alleles[0].support.num_obs, 2);
    // Kept mate's contribution = max(ln_BQ(Q=30), -2.0) ≈ -2.0.
    // Zeroed mate contributes max(ln(1)=0, -10.0) = 0.
    // Sum ≈ -2.0. If the tie-break wrongly kept mate 2, sum would
    // be max(0, -10.0) + max(ln_BQ, -2.0) ≈ -2.0 too — but the
    // distinguishable case here is the chain id count: kept
    // mate's bq is non-zero, ln_q ≈ -2.0; loser's ln_q = 0 from
    // its zeroed BQ AND -10 mq, so its contribution to q_sum = 0.
    // Net q_sum ≈ -2.0, NOT ≈ -10.0 (which would be the case if
    // the tie-break wrongly kept mate 2).
    assert!(
        rec.alleles[0].support.q_sum > -3.0 && rec.alleles[0].support.q_sum < -1.0,
        "q_sum ≈ -2.0 (first mate kept); got {}",
        rec.alleles[0].support.q_sum
    );
}

#[test]
fn mate_overlap_zeroes_lower_bq_contribution() {
    // Two paired mates overlapping at the same position with
    // different BAQ-capped BQs. Match-only agree case under S7:
    // keeper carries the *summed* BQ (capped at 200), other is
    // zeroed. Both mates still count as observations.
    //
    // Mate 1 at pos 1, BQ=30, length 3.
    // Mate 2 at pos 1, BQ=10, length 3 (overlapping).
    let fa = MockFasta::new("ACG");
    let (m1, mut m2) = paired_snp_reads("p", 1, 1, b"ACG", &[30; 3]);
    m2.bq_baq = vec![10; 3];
    let records = drive_walker(vec![m1, m2], fa);
    let rec = &records[0];
    assert_eq!(rec.alleles[0].support.num_obs, 2, "both mates count");
    // q_sum at default mq_log_err = -3.0:
    //   keeper: max(ln_perr(40), -3.0) = -3.0  (MQ dominates)
    //   other:  max(ln(1)=0, -3.0)     = 0
    //   sum ≈ -3.0
    // (the BQ-summing change from S7 is invisible here because MQ
    // dominates; tests at low MQ_log_err pin the BQ math directly).
    assert!(
        rec.alleles[0].support.q_sum > -4.0 && rec.alleles[0].support.q_sum < -2.0,
        "expected q_sum ≈ -3 (MQ-dominated), got {}",
        rec.alleles[0].support.q_sum
    );
}

#[test]
fn mate_overlap_agree_keeper_carries_summed_bq() {
    // S7 agree case: when both mates' bases agree at walker_pos,
    // the surviving observation carries the *sum* of BQs (not the
    // higher mate's BQ as the original walker did). To make BQ
    // dominate q_sum so the change is observable, set a strongly-
    // negative mq_log_err so MQ never wins the max.
    let fa = MockFasta::new("A");
    let make = |is_first: bool, bq: u8| PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 1,
        cigar: vec![CigarOp::Match(1)],
        seq: b"A".to_vec(),
        bq_baq: vec![bq],
        mq_log_err: -100.0,
        is_reverse_strand: false,
        qname: Arc::from("p"),
        mate_role: if is_first {
            MateRole::FirstOfPair
        } else {
            MateRole::SecondOfPair
        },
        adaptor_boundary: None,

        mapq: 60,
    };
    let m1 = make(true, 20);
    let m2 = make(false, 20);
    let records = drive_walker(vec![m1, m2], fa);
    assert_eq!(records.len(), 1);
    let rec = &records[0];
    assert_eq!(rec.alleles[0].support.num_obs, 2);
    // Combined BQ = 40. ln_perr(40) = -40 * ln(10) / 10 ≈ -9.21.
    // Keeper contribution: max(-9.21, -100) = -9.21.
    // Other contribution: max(ln_perr(0)=0, -100) = 0.
    // Total q_sum ≈ -9.21. Pre-S7 (Q=20 unsummed): ≈ -4.61.
    let q = rec.alleles[0].support.q_sum;
    assert!(
        q < -8.5 && q > -10.0,
        "q_sum should reflect summed BQ (≈ ln_perr(40) ≈ -9.21), got {q}",
    );
}

#[test]
fn mate_overlap_agree_combined_bq_caps_at_200() {
    // S7 agree case: the combined BQ is clamped to samtools'
    // MPLP cap of 200 (htslib/sam.c:5919-5921). 150 + 100 = 250
    // → capped at 200.
    let fa = MockFasta::new("A");
    let make = |is_first: bool, bq: u8| PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 1,
        cigar: vec![CigarOp::Match(1)],
        seq: b"A".to_vec(),
        bq_baq: vec![bq],
        mq_log_err: -100.0,
        is_reverse_strand: false,
        qname: Arc::from("p"),
        mate_role: if is_first {
            MateRole::FirstOfPair
        } else {
            MateRole::SecondOfPair
        },
        adaptor_boundary: None,

        mapq: 60,
    };
    let m1 = make(true, 150);
    let m2 = make(false, 100);
    let records = drive_walker(vec![m1, m2], fa);
    let q = records[0].alleles[0].support.q_sum;
    // ln_perr(200) ≈ -46.05. Without the cap it would be
    // ln_perr(250) ≈ -57.56.
    assert!(
        q < -45.0 && q > -47.0,
        "q_sum should reflect cap-200 (≈ -46.05), got {q}",
    );
}

#[test]
fn mate_overlap_disagree_winner_bq_scaled_by_0_8() {
    // S7 disagree case: when mate bases disagree, the higher-BQ
    // mate keeps its BQ scaled by 0.8 (samtools' "we trust this
    // less" haircut at htslib/sam.c:5927); the loser is zeroed.
    let fa = MockFasta::new("A");
    let make = |is_first: bool, base: u8, bq: u8| PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 1,
        cigar: vec![CigarOp::Match(1)],
        seq: vec![base],
        bq_baq: vec![bq],
        mq_log_err: -100.0,
        is_reverse_strand: false,
        qname: Arc::from("p"),
        mate_role: if is_first {
            MateRole::FirstOfPair
        } else {
            MateRole::SecondOfPair
        },
        adaptor_boundary: None,

        mapq: 60,
    };
    // Mate 1 has REF "A" with higher BQ → winner. Mate 2 has SNP
    // "G" with lower BQ → loser, zeroed.
    let m1 = make(true, b'A', 30);
    let m2 = make(false, b'G', 20);
    let records = drive_walker(vec![m1, m2], fa);
    let rec = &records[0];
    let ref_allele = rec
        .alleles
        .iter()
        .find(|a| a.seq.as_slice() == b"A")
        .expect("REF allele present");
    let snp_allele = rec
        .alleles
        .iter()
        .find(|a| a.seq.as_slice() == b"G")
        .expect("SNP allele present");
    assert_eq!(ref_allele.support.num_obs, 1);
    assert_eq!(snp_allele.support.num_obs, 1);
    // Winner BQ = (30 * 0.8) as u8 = 24. ln_perr(24) ≈ -5.53.
    // Pre-S7 (Q=30 unscaled): ≈ -6.91.
    let q_ref = ref_allele.support.q_sum;
    assert!(
        q_ref < -5.0 && q_ref > -6.0,
        "REF allele q_sum should reflect scaled BQ=24 (≈ -5.53), got {q_ref}",
    );
    // Loser BQ zeroed → ln_perr(0) = 0 → max(0, -100) = 0.
    assert_eq!(
        snp_allele.support.q_sum, 0.0,
        "SNP allele's BQ was zeroed; q_sum should be 0",
    );
}

#[test]
fn record_widen_events_counter_only_increments_on_real_widens() {
    // Three pure-Match reads at consecutive positions on a 10-bp
    // reference. Every record opens fresh at span 1 and never
    // widens; the run's `record_widen_events` counter should
    // therefore stay at 0. The previous implementation summed
    // ref_span across all open records before/after each
    // process_position step and incremented when the after-sum
    // grew — but a freshly-opened record also grows the sum, so
    // the counter conflated opens with widens.
    let fa = MockFasta::new("ACGTACGTAC");
    let r1 = snp_read("r1", 1, b"ACG", &[30; 3]);
    let r2 = snp_read("r2", 4, b"TAC", &[30; 3]);
    let r3 = snp_read("r3", 7, b"GTA", &[30; 3]);
    let (_records, summary) = drive_walker_with_summary(vec![r1, r2, r3], fa);
    assert_eq!(
        summary.record_widen_events, 0,
        "no widening occurred; counter should be 0, got {}",
        summary.record_widen_events
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
        mate_role: MateRole::FirstOfPair,
        adaptor_boundary: None,

        mapq: 60,
    };
    let mut mate_b = mate_a.clone();
    mate_b.mate_role = MateRole::SecondOfPair;
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
        ins.support.num_obs, 1,
        "indel-overlap collapses to one observation; got {}",
        ins.support.num_obs
    );
    // Forward-strand count should also reflect a single observation.
    assert_eq!(ins.support.fwd, 1);
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
    let err = first_walker_error(vec![r1, r2], fa);
    let msg = err.to_string();
    assert!(msg.contains("out-of-order"), "got: {msg}");
}

#[test]
fn chromosome_id_regression_is_a_hard_error() {
    // Spec invariant: `(chrom_id, alignment_start)` must be
    // monotonically non-decreasing across the input stream. The
    // within-chromosome regression case is already pinned by
    // `out_of_order_input_is_a_hard_error`; this test pins the
    // *cross-chromosome* case, which used to be silently accepted
    // because `flush_chromosome` reset `last_admitted_locus = None`
    // and the regressing read then admitted as a fresh start on
    // the smaller chrom_id. See finding M2 in
    // `ia/reviews/pileup_2026-05-09.md`.
    let fa = MockFasta::with_chromosomes(&["ACGTA", "ACGTA"]);
    let mut r1 = snp_read("a", 1, b"AC", &[30; 2]);
    r1.chrom_id = 1;
    let mut r2 = snp_read("b", 1, b"AC", &[30; 2]);
    r2.chrom_id = 0;
    let err = first_walker_error(vec![r1, r2], fa);
    let msg = err.to_string();
    assert!(
        msg.contains("out-of-order"),
        "error message should reference out-of-order semantics; got: {msg}",
    );
}

#[test]
fn forward_chromosome_change_is_accepted() {
    // Pin the legitimate forward case so the M2 fix doesn't
    // accidentally reject monotonic chrom_id transitions.
    let fa = MockFasta::with_chromosomes(&["ACG", "TTT"]);
    let mut r1 = snp_read("a", 1, b"AC", &[30; 2]);
    r1.chrom_id = 0;
    let mut r2 = snp_read("b", 1, b"TT", &[30; 2]);
    r2.chrom_id = 1;
    let records = drive_walker(vec![r1, r2], fa);
    assert!(
        records.iter().any(|r| r.chrom_id == 0),
        "chrom 0 records must be emitted",
    );
    assert!(
        records.iter().any(|r| r.chrom_id == 1),
        "chrom 1 records must be emitted",
    );
    // And they must come in chrom order.
    let chrom_ids: Vec<u32> = records.iter().map(|r| r.chrom_id).collect();
    let mut sorted = chrom_ids.clone();
    sorted.sort();
    assert_eq!(chrom_ids, sorted, "records must emit in chrom_id order");
}

#[test]
fn chain_ids_are_unique_and_monotonically_allocated() {
    // Three solo reads admit cleanly; every chain id observed across
    // emitted records must be unique (no recycling) and the ids
    // appear in non-decreasing order of first reference position.
    let fa = MockFasta::new("ACGTACGTAC");
    let reads = vec![
        snp_read("a", 1, b"AC", &[30; 2]),
        snp_read("b", 4, b"AC", &[30; 2]),
        snp_read("c", 7, b"AC", &[30; 2]),
    ];
    let records = drive_walker(reads, fa);

    // Collect every chain id from every allele observation.
    let mut all_ids: Vec<u64> = records
        .iter()
        .flat_map(|r| r.alleles.iter().flat_map(|a| a.chain_ids.iter().copied()))
        .collect();
    let n_total_observations = all_ids.len();
    all_ids.sort_unstable();
    all_ids.dedup();
    assert_eq!(
        all_ids.len(),
        3,
        "three reads must produce exactly three distinct chain ids \
         (observed across {n_total_observations} allele observations)",
    );
    // Monotonic: starts at 0 and increases by 1.
    assert_eq!(all_ids, vec![0, 1, 2]);
}

#[test]
fn paired_mates_share_a_single_chain_id() {
    // First and second mates of a pair must collapse to one chain id.
    let fa = MockFasta::new("ACG");
    let (m1, m2) = paired_snp_reads("p", 1, 1, b"ACG", &[30; 3]);
    let records = drive_walker(vec![m1, m2], fa);
    let mut all_ids: Vec<u64> = records
        .iter()
        .flat_map(|r| r.alleles.iter().flat_map(|a| a.chain_ids.iter().copied()))
        .collect();
    all_ids.sort_unstable();
    all_ids.dedup();
    assert_eq!(
        all_ids.len(),
        1,
        "a paired read produces exactly one chain id shared across both mates",
    );
}

#[test]
fn chain_ids_persist_across_chromosome_boundaries() {
    // Chain ids are per-`.psp`-file unique, not per-chromosome. The
    // walker's chain-id allocator must NOT reset `next_id` on
    // chromosome change.
    let fa = MockFasta::with_chromosomes(&["AC", "AC"]);
    let mut r0 = snp_read("a", 1, b"AC", &[30; 2]);
    r0.chrom_id = 0;
    let mut r1 = snp_read("b", 1, b"AC", &[30; 2]);
    r1.chrom_id = 1;
    let records = drive_walker(vec![r0, r1], fa);
    let mut all_ids: Vec<u64> = records
        .iter()
        .flat_map(|r| r.alleles.iter().flat_map(|a| a.chain_ids.iter().copied()))
        .collect();
    all_ids.sort_unstable();
    all_ids.dedup();
    assert_eq!(
        all_ids,
        vec![0, 1],
        "chain ids must remain unique across chromosomes",
    );
}

#[test]
fn column_depth_cap_truncates_snp_only_column_when_over_cap() {
    // Five SNP-only reads anchored at pos 1, each spanning 5
    // bases. Every column has 5 contributors and only Match
    // events, so the SNP cap applies. With max_snp_column_depth=3
    // we expect every column to truncate to 3 contributors.
    let fa = MockFasta::new("ACGTA");
    let reads: Vec<_> = (0..5)
        .map(|i| snp_read(&format!("r{i}"), 1, b"ACGTA", &[30; 5]))
        .collect();
    let cfg = WalkerConfig {
        max_snp_column_depth: 3,
        max_indel_column_depth: 99,
        ..WalkerConfig::default()
    };
    let (records, summary) = drive_walker_with_config(reads, fa, &cfg);

    // 5 columns, all over-cap → 5 truncations.
    assert_eq!(summary.column_depth_truncations, 5);
    // No column should report num_obs > cap.
    for rec in &records {
        for allele in &rec.alleles {
            assert!(
                allele.support.num_obs <= 3,
                "pos {}: num_obs {} should be capped at 3",
                rec.pos,
                allele.support.num_obs,
            );
        }
    }
}

#[test]
fn column_depth_cap_keeps_first_n_of_admission_order() {
    // The cap is a defensive truncation, not a uniform sampler.
    // Pin the contract so an operator can reason about what
    // survives: the first `cap` contributors in admission order
    // (= upstream coordinate-then-arrival order from `cram_input`)
    // fold into the record; the rest are dropped. Mi5 in
    // `ia/reviews/pileup_2026-05-09.md`.
    //
    // Five SNP reads at pos 1, each with a distinct read base so
    // we can identify which contributors survived. Cap = 2 → only
    // r0 (C) and r1 (G) fold; r2 (T), r3 (A=REF), and r4 (C) are
    // dropped.
    let fa = MockFasta::new("A");
    let r0 = snp_read("r0", 1, b"C", &[30]);
    let r1 = snp_read("r1", 1, b"G", &[30]);
    let r2 = snp_read("r2", 1, b"T", &[30]);
    let r3 = snp_read("r3", 1, b"A", &[30]); // matches REF
    let r4 = snp_read("r4", 1, b"C", &[30]);
    let cfg = WalkerConfig {
        max_snp_column_depth: 2,
        max_indel_column_depth: 99,
        ..WalkerConfig::default()
    };
    let (records, summary) = drive_walker_with_config(vec![r0, r1, r2, r3, r4], fa, &cfg);

    assert_eq!(records.len(), 1, "single column emitted");
    assert_eq!(summary.column_depth_truncations, 1);
    let rec = &records[0];

    // REF "A" survives as the alleles[0] entry but with zero
    // observations, because r3 (the only REF-matching read) was
    // past the cap.
    assert_eq!(rec.alleles[0].seq, b"A", "alleles[0] is REF");
    assert_eq!(
        rec.alleles[0].support.num_obs, 0,
        "no surviving read matched REF",
    );

    // r0's "C" and r1's "G" must be present with one observation each.
    let c = rec
        .alleles
        .iter()
        .find(|a| a.seq.as_slice() == b"C")
        .expect("r0's allele 'C' must survive");
    assert_eq!(c.support.num_obs, 1);
    let g = rec
        .alleles
        .iter()
        .find(|a| a.seq.as_slice() == b"G")
        .expect("r1's allele 'G' must survive");
    assert_eq!(g.support.num_obs, 1);

    // r2's "T" must be absent — past the cap.
    assert!(
        rec.alleles.iter().all(|a| a.seq.as_slice() != b"T"),
        "r2 was past the cap and must not have folded",
    );

    // Total buckets: REF + r0 + r1.
    assert_eq!(rec.alleles.len(), 3, "REF + 2 surviving SNP buckets");
}

#[test]
fn column_depth_cap_uses_indel_cap_when_any_indel_event_present() {
    // Four SNP-only reads + one indel-bearing read, all anchored
    // at pos 1. At column pos 1 the indel-bearing read contributes
    // an Insertion event — the column flips to "indel column" and
    // the tighter indel cap applies. At pos 2 onward only Match
    // events remain, so the (much higher) SNP cap applies and
    // does not fire.
    let fa = MockFasta::new("AAAACGT");
    let mut reads: Vec<PreparedRead> = (0..4)
        .map(|i| snp_read(&format!("snp{i}"), 1, b"AAAACGT", &[30; 7]))
        .collect();
    // Indel read: 1M 2I 5M starting at pos 1; anchor of the
    // insertion is pos 1.
    let indel = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 6,
        cigar: vec![CigarOp::Match(1), CigarOp::Insertion(2), CigarOp::Match(5)],
        seq: b"AXXAAACG".to_vec(),
        bq_baq: vec![30; 8],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("indel"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    reads.push(indel);

    let cfg = WalkerConfig {
        max_snp_column_depth: 99,  // far above 5; SNP-only cols don't fire
        max_indel_column_depth: 2, // below 5; indel col at pos 1 fires
        ..WalkerConfig::default()
    };
    let (_records, summary) = drive_walker_with_config(reads, fa, &cfg);

    // Exactly one column carried an indel event (pos 1), and that
    // column had 5 contributors > indel cap of 2 → one truncation.
    assert_eq!(
        summary.column_depth_truncations, 1,
        "indel cap should fire exactly once at the indel-anchor column",
    );
}

#[test]
fn column_depth_cap_does_not_fire_below_threshold() {
    // Two SNP-only reads, default config (caps 8000 and 250).
    // Far below either cap → no truncation, every contributor
    // folds.
    let fa = MockFasta::new("ACGTA");
    let r1 = snp_read("r1", 1, b"ACGTA", &[30; 5]);
    let r2 = snp_read("r2", 1, b"ACGTA", &[30; 5]);
    let (records, summary) = drive_walker_with_summary(vec![r1, r2], fa);

    assert_eq!(summary.column_depth_truncations, 0);
    for rec in &records {
        assert_eq!(
            rec.alleles[0].support.num_obs, 2,
            "pos {}: both reads should fold (no truncation under default cap)",
            rec.pos,
        );
    }
}

// --- G1: adaptor-region per-base filter, walker integration ---------

#[test]
fn g1_walker_drops_match_observations_past_adaptor_boundary() {
    // Reference: ACGTACGT (positions 1..9 1-based, 8 bases).
    // Two reads on the same molecule, ancient-DNA shape:
    //   - r_fwd: forward-strand, M(8) at pos 1, REF on every base.
    //     Adaptor boundary at 5 (insert size = 4 < seq_len = 8) →
    //     positions 1..4 emit, positions 5..8 are dropped as adaptor.
    //   - r_rev: reverse-strand, M(8) at pos 1, REF on every base.
    //     Adaptor boundary at 4 (mate.start - 1) → positions 1..4
    //     are dropped, positions 5..8 emit.
    //
    // Net result at every position 1..8: exactly one read contributes
    // (the one whose strand has that position outside its adaptor).
    // Without G1, both reads would fold at every position and num_obs
    // would be 2 — this is the regression the filter prevents.
    let fa = MockFasta::new("ACGTACGT");
    let mut r_fwd = snp_read("pair", 1, b"ACGTACGT", &[30; 8]);
    r_fwd.is_reverse_strand = false;
    r_fwd.mate_role = MateRole::FirstOfPair;
    r_fwd.adaptor_boundary = Some(5);
    let mut r_rev = snp_read("pair", 1, b"ACGTACGT", &[30; 8]);
    r_rev.is_reverse_strand = true;
    r_rev.mate_role = MateRole::SecondOfPair;
    r_rev.adaptor_boundary = Some(4);

    let records = drive_walker(vec![r_fwd, r_rev], fa);
    assert_eq!(records.len(), 8, "every covered position emits one record");
    for rec in &records {
        // Each position is covered by exactly one of the two mates
        // after the adaptor filter applies. Without G1, mate-overlap
        // resolution would still cap to 1 (sum-and-cap on agreement
        // applies), but that path mutates BQ; here the filter cleanly
        // removes the adaptor base from the contributor list before
        // overlap resolution sees it.
        assert_eq!(
            rec.alleles[0].support.num_obs, 1,
            "pos {}: exactly one mate is outside adaptor at this position",
            rec.pos,
        );
    }
}

// --- Error-variant coverage --------------------------------------
//
// M14–M17 of `ia/reviews/pileup_2026-05-11.md`: every `WalkerError`
// variant must have a regression test pinning the exact variant
// returned. Without these, a swallowed or mis-mapped error could
// silently degrade to `None` / `Ok(())` (data loss) or to the wrong
// variant (operator triage misled).
//
// (The pre-iterator `run_returns_channel_closed_when_receiver_dropped_mid_stream`
// regression has been removed: the pull-shape walker no longer
// owns or sends through a channel, and the `ChannelClosed` variant
// is gone.)

/// Iterate the walker until it yields its first error and return
/// it. Panics if the walker exhausts without erroring.
fn first_walker_error(reads: Vec<PreparedRead>, fa: MockFasta) -> super::WalkerError {
    run(reads, &fa, &WalkerConfig::default())
        .find_map(|r| r.err())
        .expect("walker did not surface any error")
}

#[test]
fn zero_ref_span_input_is_a_hard_error() {
    // M15. A malformed PreparedRead with alignment_end <
    // alignment_start must hard-error on admission.
    let fa = MockFasta::new("ACGT");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 3,
        alignment_end: 2,
        cigar: vec![CigarOp::Match(0)],
        seq: vec![],
        bq_baq: vec![],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("zero"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    let err = first_walker_error(vec![r], fa);
    assert!(
        matches!(err, super::WalkerError::ZeroRefSpan { .. }),
        "got {err:?}",
    );
}

#[test]
fn open_record_widening_past_max_record_span_errors() {
    // M16 widen path. A single read with a deletion that pushes the
    // open record's footprint past MAX_RECORD_SPAN must error.
    let ref_len = (super::DEFAULT_MAX_RECORD_SPAN as usize) + 50;
    let fa = MockFasta::new(&"A".repeat(ref_len));
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: ref_len as u32,
        cigar: vec![
            CigarOp::Match(1),
            CigarOp::Deletion(super::DEFAULT_MAX_RECORD_SPAN + 1),
            CigarOp::Match(1),
        ],
        seq: b"AA".to_vec(),
        bq_baq: vec![30, 30],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("wide"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    let err = first_walker_error(vec![r], fa);
    match err {
        super::WalkerError::RecordTooWide { cap, .. } => {
            assert_eq!(cap, super::DEFAULT_MAX_RECORD_SPAN);
        }
        other => panic!("expected RecordTooWide, got {other:?}"),
    }
}

// M21: malformed-PreparedRead validation. The walker must reject
// the upstream contract violation as a typed error, not panic on
// out-of-bounds indexing.

#[test]
fn admit_rejects_seq_shorter_than_cigar_consumes() {
    // CIGAR consumes 5 read bases (5M), but seq has only 3.
    let fa = MockFasta::new("ACGTA");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 5,
        cigar: vec![CigarOp::Match(5)],
        seq: b"ACG".to_vec(),
        bq_baq: vec![30; 3],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("short"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    let err = first_walker_error(vec![r], fa);
    match err {
        super::WalkerError::MalformedRead { reason, .. } => {
            assert!(
                reason.contains("CIGAR consumes 5 read bases but seq.len = 3"),
                "got reason: {reason}",
            );
        }
        other => panic!("expected MalformedRead, got {other:?}"),
    }
}

#[test]
fn admit_rejects_seq_bq_length_mismatch() {
    let fa = MockFasta::new("ACGTA");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 5,
        cigar: vec![CigarOp::Match(5)],
        seq: b"ACGTA".to_vec(),
        bq_baq: vec![30, 30, 30], // 3 instead of 5
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("bq_short"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    let err = first_walker_error(vec![r], fa);
    match err {
        super::WalkerError::MalformedRead { reason, .. } => {
            assert!(
                reason.contains("seq.len (5) != bq_baq.len (3)"),
                "got reason: {reason}",
            );
        }
        other => panic!("expected MalformedRead, got {other:?}"),
    }
}

#[test]
fn admit_rejects_cigar_consuming_more_read_bases_than_seq_provides() {
    // CIGAR = 2M + 4I + 2M. Consumes 8 read bases; seq has 5.
    // The cursor would index seq[2..6] for the Insertion, which
    // is OOB. Admit-time check must catch it first.
    let fa = MockFasta::new("ACGTA");
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 4,
        cigar: vec![CigarOp::Match(2), CigarOp::Insertion(4), CigarOp::Match(2)],
        seq: b"ACGTA".to_vec(),
        bq_baq: vec![30; 5],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("cigar_long"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    let err = first_walker_error(vec![r], fa);
    match err {
        super::WalkerError::MalformedRead { reason, .. } => {
            assert!(
                reason.contains("CIGAR consumes 8 read bases but seq.len = 5"),
                "got reason: {reason}",
            );
        }
        other => panic!("expected MalformedRead, got {other:?}"),
    }
}

#[test]
fn fasta_fetch_failure_propagates_as_walker_error_fasta() {
    // M17. Reference shorter than the read — MockFasta::fetch
    // returns UnexpectedEof when the open_new fetch runs off the
    // chromosome's end. Must surface as WalkerError::Fasta with
    // the correct locus.
    let fa = MockFasta::new("AC");
    let r = snp_read("r", 1, b"ACGT", &[30; 4]);
    let err = first_walker_error(vec![r], fa);
    match err {
        super::WalkerError::Fasta {
            chrom_id,
            start,
            start_plus_len,
            ..
        } => {
            assert_eq!(chrom_id, 0);
            assert!(start_plus_len > start);
        }
        other => panic!("expected Fasta, got {other:?}"),
    }
}

#[test]
fn walker_iterator_returns_none_after_yielding_error() {
    // The iterator's terminal-on-first-error contract: once `next()`
    // returns `Err`, every subsequent call returns `None`. Mirrors
    // the previous push-shape behaviour, where `run` returned `Err`
    // and stopped emitting at the same moment.
    let fa = MockFasta::new("AC"); // too short to satisfy a 4-base read.
    let r = snp_read("r", 1, b"ACGT", &[30; 4]);
    let mut walker = run(vec![r], &fa, &WalkerConfig::default());
    // The error may be preceded by zero or more successful records
    // (depending on the order of fasta fetch vs. walker_pos). Drive
    // until we see the first Err, then assert exhaustion.
    let mut saw_error = false;
    for item in &mut walker {
        if item.is_err() {
            saw_error = true;
            break;
        }
    }
    assert!(saw_error, "expected the walker to surface a Fasta error");
    assert!(
        walker.next().is_none(),
        "walker must return None after its first error",
    );
}

// ---------------------------------------------------------------------
// PreparedRead::length — internal length-invariant checks
// ---------------------------------------------------------------------

#[test]
fn prepared_read_length_returns_validated_length_on_well_formed_read() {
    // 4-bp Match: seq.len = bq_baq.len = cigar-consumed = 4.
    let r = snp_read("r", 1, b"ACGT", &[30; 4]);
    assert_eq!(r.length(), Ok(4));
}

#[test]
fn prepared_read_length_with_insertion_includes_inserted_bases() {
    // CIGAR 2M1I2M consumes 5 read bases. seq and bq_baq must
    // match. Confirms the helper sums the I op (read-consuming)
    // and not just the M ops.
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 4,
        cigar: vec![CigarOp::Match(2), CigarOp::Insertion(1), CigarOp::Match(2)],
        seq: b"ACXGT".to_vec(),
        bq_baq: vec![30; 5],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    assert_eq!(r.length(), Ok(5));
}

#[test]
fn prepared_read_length_with_deletion_excludes_deleted_bases() {
    // CIGAR 2M2D2M consumes 4 read bases (D op is reference-only).
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 6,
        cigar: vec![CigarOp::Match(2), CigarOp::Deletion(2), CigarOp::Match(2)],
        seq: b"ACGT".to_vec(),
        bq_baq: vec![30; 4],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    assert_eq!(r.length(), Ok(4));
}

#[test]
fn prepared_read_length_seq_bq_mismatch_returns_typed_error() {
    let mut r = snp_read("r", 1, b"ACGT", &[30; 4]);
    r.bq_baq.pop(); // bq_baq now length 3, seq length 4.
    assert_eq!(
        r.length(),
        Err(super::ReadLengthError::SeqBqMismatch {
            seq_len: 4,
            bq_baq_len: 3,
        })
    );
}

#[test]
fn prepared_read_length_cigar_seq_mismatch_returns_typed_error() {
    // seq.len = bq_baq.len = 4 but cigar consumes 3 read bases.
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 3,
        cigar: vec![CigarOp::Match(3)],
        seq: b"ACGT".to_vec(),
        bq_baq: vec![30; 4],
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    assert_eq!(
        r.length(),
        Err(super::ReadLengthError::CigarSeqMismatch {
            cigar_consumed: 3,
            seq_len: 4,
        })
    );
}

#[test]
fn prepared_read_length_checks_seq_bq_before_cigar() {
    // Both invariants violated; the seq/bq check runs first so
    // SeqBqMismatch wins. Locks in the deterministic ordering.
    let r = PreparedRead {
        chrom_id: 0,
        alignment_start: 1,
        alignment_end: 3,
        cigar: vec![CigarOp::Match(3)],
        seq: b"ACGT".to_vec(), // length 4
        bq_baq: vec![30; 3],   // length 3 ≠ seq.len → SeqBqMismatch
        mq_log_err: -3.0,
        is_reverse_strand: false,
        qname: Arc::from("r"),
        mate_role: MateRole::Solo,
        adaptor_boundary: None,

        mapq: 60,
    };
    assert!(matches!(
        r.length(),
        Err(super::ReadLengthError::SeqBqMismatch { .. })
    ));
}
