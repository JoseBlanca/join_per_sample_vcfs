//! ng step 1 — read filtering: the whole-read keep/drop prelude. It takes the
//! reads of one sample's alignment file (BAM/CRAM) and yields the subset worth
//! carrying forward, plus a running tally of what was dropped and why. Every
//! decision is per-read, locus-independent, and content-preserving — filtering
//! *selects* reads, it never rewrites them.
//!
//! Design: `doc/devel/ng/spec/read_filtering.md` (the "why"),
//! `doc/devel/ng/arch/read_filtering.md` (types & interfaces).
//!
//! Read filtering is a **port** of the production filter stack in
//! [`crate::bam::alignment_input`]: it reuses that module's pure predicates
//! (`read_exceeds_mismatch_fraction`, `cigar_is_bad`), its `FLAG_*` /
//! `DEFAULT_*` constants, and its `RecordBuf → MappedRead` decode path as-is,
//! and supplies only its own driver and config.
//!
//! Milestones A (the step-1-local types) and B (the two-phase `verdict_*`
//! cascade) have landed; the record-source seam (C) and the `ReadFilter`
//! iterator (D) follow.

use crate::bam::alignment_input::{
    DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH,
    DEFAULT_MISMATCH_BQ_FLOOR, FLAG_DUPLICATE, FLAG_QC_FAIL, FLAG_SECONDARY, FLAG_SUPPLEMENTARY,
    FLAG_UNMAPPED, MappedRead, cigar_is_bad, cigar_ref_span, read_exceeds_mismatch_fraction,
};
use crate::ng::ref_seq::{RawRefSeq, RefSeqError};
use crate::ng::types::{BaseQual, Bp, ContigId, MapQual, MismatchFraction};

/// The filtering policy: which filters are active and their thresholds. Minimal
/// by design — one field per active filter, no dormant levers (downsampling,
/// read pooling enter only when they enter the pipeline). [`Default`] is the
/// production policy the lab runs with, its thresholds the reused `DEFAULT_*`
/// constants from [`crate::bam::alignment_input`]. Mirrors the filtering subset
/// of the existing `AlignmentMergedReaderConfig`.
///
/// `Option<T>` is the "no threshold" state, never a sentinel: `None` means *no
/// minimum*, `Some(q)` means *drop below `q`* — so `Some(0)` (drop nothing, but
/// a threshold is set) stays structurally distinct from `None` (no threshold).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReadFilterConfig {
    /// `None` = no minimum; `Some(q)` = drop reads with MAPQ `< q` (filter #2).
    pub min_mapq: Option<MapQual>,
    /// `None` = no minimum; `Some(n)` = drop reads shorter than `n` bp (#7).
    pub min_read_length: Option<Bp>,
    /// Drop reads flagged QC-fail (`FLAG_QC_FAIL`) — filter #6.
    pub drop_qc_fail: bool,
    /// Drop reads flagged as PCR/optical duplicates (`FLAG_DUPLICATE`) — #1.
    pub drop_duplicate: bool,
    /// `None` = filter #8 off (no reference access at all); `Some(x)` = drop a
    /// read whose quality-clearing `M`-op mismatch fraction exceeds `x`.
    pub max_read_mismatch_fraction: Option<MismatchFraction>,
    /// BQ floor below which a mismatch does not count toward filter #8. Only
    /// meaningful when `max_read_mismatch_fraction` is `Some`.
    pub mismatch_bq_floor: BaseQual,
}

impl Default for ReadFilterConfig {
    fn default() -> Self {
        Self {
            min_mapq: Some(MapQual(DEFAULT_MIN_MAPQ)),
            min_read_length: Some(Bp(DEFAULT_MIN_READ_LENGTH)),
            drop_qc_fail: true,
            drop_duplicate: true,
            // PANIC-FREE: the default fraction is a known-good in-range constant
            // (0.10, `DEFAULT_MAX_READ_MISMATCH_FRACTION`), so the checked
            // constructor cannot fail here. Guarded by the
            // `default_config_reproduces_the_production_filter_policy` test, which
            // exercises this exact path.
            max_read_mismatch_fraction: Some(
                MismatchFraction::try_new(DEFAULT_MAX_READ_MISMATCH_FRACTION)
                    .expect("DEFAULT_MAX_READ_MISMATCH_FRACTION is in [0, 1]"),
            ),
            mismatch_bq_floor: BaseQual(DEFAULT_MISMATCH_BQ_FLOOR),
        }
    }
}

/// The verdict for one read. `Keep` carries the read on; `Drop` records which
/// filter fired — the first one, per the hit-rate order — for the tally.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterVerdict {
    Keep,
    Drop(DropReason),
}

/// Which filter dropped a read. Variant names line up 1:1 with the
/// [`ReadFilterCounts`] fields, so a drop maps to exactly one counter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropReason {
    Duplicate,
    LowMapq,
    Supplementary,
    Secondary,
    Unmapped,
    QcFail,
    TooShort,
    HighMismatchFraction,
    BadCigar,
}

/// A per-sample tally of the filtering pass — one counter per drop reason, plus
/// the kept count. The ng port of `FilterCounts`. Surfacing every drop is the
/// "no silent caps" discipline: a read that vanished must be accounted for. It
/// is a **running** tally — readable at any point, final once the input is
/// exhausted.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReadFilterCounts {
    pub kept: u64,
    pub duplicate: u64,
    pub low_mapq: u64,
    pub supplementary: u64,
    pub secondary: u64,
    pub unmapped: u64,
    pub qc_fail: u64,
    pub too_short: u64,
    pub high_mismatch_fraction: u64,
    pub bad_cigar: u64,
}

/// Phase one of the cascade — the flag/MAPQ filters (#1–#6), decided on an
/// undecoded record's `flag` and `mapq` alone. Reference-free and decode-free:
/// `Keep` means "decode and continue to phase two", `Drop` is charged to the
/// first filter that fires. Order is identical to production's
/// `classify_pre_decode` (hit-rate-ordered: duplicate, low-MAPQ, supplementary,
/// secondary, unmapped, QC-fail).
///
/// `mapq` is already resolved: SAM's "unavailable" (`0xFF`) is mapped to
/// `MapQual(0)` by the record source (Milestone C), so a non-zero `min_mapq`
/// drops it — matching production. `flag` is the raw SAM bitfield
/// (`MappedRead.flag`), tested against the reused `FLAG_*` constants.
// Wired into the `ReadFilter` iterator in Milestone D; unit-tested standalone here.
#[cfg_attr(not(test), allow(dead_code))]
fn verdict_pre_decode(flag: u16, mapq: MapQual, config: &ReadFilterConfig) -> FilterVerdict {
    // 1. Duplicate — a PCR/optical copy of another molecule (toggle).
    if config.drop_duplicate && (flag & FLAG_DUPLICATE) != 0 {
        return FilterVerdict::Drop(DropReason::Duplicate);
    }
    // 2. Low MAPQ — the aligner is unsure of the placement (toggle via threshold).
    if let Some(min) = config.min_mapq
        && mapq < min
    {
        return FilterVerdict::Drop(DropReason::LowMapq);
    }
    // 3. Supplementary — a chunk of a chimeric read (unconditional).
    if (flag & FLAG_SUPPLEMENTARY) != 0 {
        return FilterVerdict::Drop(DropReason::Supplementary);
    }
    // 4. Secondary — a duplicate projection of a primary alignment (unconditional).
    if (flag & FLAG_SECONDARY) != 0 {
        return FilterVerdict::Drop(DropReason::Secondary);
    }
    // 5. Unmapped — no placement, so no allele evidence (unconditional).
    if (flag & FLAG_UNMAPPED) != 0 {
        return FilterVerdict::Drop(DropReason::Unmapped);
    }
    // 6. QC fail — the sequencer/pipeline flagged the read (toggle).
    if config.drop_qc_fail && (flag & FLAG_QC_FAIL) != 0 {
        return FilterVerdict::Drop(DropReason::QcFail);
    }
    FilterVerdict::Keep
}

/// Phase two of the cascade — the decode-dependent filters, run on a decoded
/// [`MappedRead`] only after it clears phase one. Evaluated **cheapest-first**:
///
/// 1. **#7 too-short** — one length compare, no reference.
/// 2. **#9 bad-CIGAR** — a pure CIGAR scan, no reference.
/// 3. **#8 high-mismatch** — a reference fetch plus a per-base walk, and only
///    when `max_read_mismatch_fraction` is `Some`.
///
/// This is a deliberate reordering relative to the spec's #7/#8/#9 *table*: it
/// puts the one reference-touching filter (#8) last, so a read dropped for
/// being too short or for a malformed CIGAR never pays the reference fetch and
/// base walk. It honours the spec's stated "cheapest, most-often-firing first"
/// principle (spec §3), charges a both-failing read to the root cause
/// (`BadCigar`) rather than the symptom (`HighMismatchFraction`), and leaves the
/// keep/drop *set* unchanged (the filters are independent — a read failing any
/// is dropped regardless of order). The mismatch fraction is measured on the
/// read's own (un-left-aligned) CIGAR: left-alignment only shifts indels across
/// equal bases, so it does not change the match/mismatch tally, and it is
/// deferred to `pileup/` anyway (spec §6).
///
/// `ref_buf` is a caller-owned scratch buffer the reference bytes are read into
/// (reused across reads, so #8 costs no per-read allocation). It is written only
/// when #8 actually runs. A [`RefSeqError`] from the fetch is **fatal to the
/// run**, propagated rather than swallowed into a drop or a keep. This includes
/// an `OutOfBounds` window running past the contig end: a validly-aligned read
/// never covers reference positions the contig does not have, so an
/// out-of-bounds fetch signals a malformed record — and the fatal model treats
/// it, like a truncated file, as corrupt input to fail loudly on rather than
/// filter around (spec §7).
// Wired into the `ReadFilter` iterator in Milestone D; unit-tested standalone here.
#[cfg_attr(not(test), allow(dead_code))]
fn verdict_post_decode(
    read: &MappedRead,
    reference: &impl RawRefSeq,
    config: &ReadFilterConfig,
    ref_buf: &mut Vec<u8>,
) -> Result<FilterVerdict, RefSeqError> {
    // #7 — too short (decoded SEQ length). Cheapest: no CIGAR walk, no reference.
    if let Some(min) = config.min_read_length
        && (read.seq.len() as u32) < min.get()
    {
        return Ok(FilterVerdict::Drop(DropReason::TooShort));
    }

    // #9 — bad CIGAR (adjacent I/D, or a boundary deletion). Cheap: a pure scan
    // of the aligner's CIGAR, no reference. Run before #8 so a malformed read
    // never pays the reference fetch below.
    if cigar_is_bad(&read.cigar) {
        return Ok(FilterVerdict::Drop(DropReason::BadCigar));
    }

    // #8 — high mismatch fraction. The only reference-dependent filter, and the
    // only post-decode one that allocates work; skipped entirely when disabled.
    if let Some(max) = config.max_read_mismatch_fraction {
        let ref_span = cigar_ref_span(&read.cigar);
        // PANIC-FREE: reference coordinates are `u32` at the RefSeq boundary
        // (ref_seq.md, Decision 3). `ref_id` indexes the `u32` contig table and a
        // mapped read's 1-based position fits `u32` for any real contig
        // (< 4.29 Gbp), so neither conversion truncates on real input. A value
        // that did not fit would be a corrupt record — failing loudly is the
        // intended response under the fatal error model, not a silent `as`
        // truncation that would fetch the wrong window and mis-verdict the read.
        let contig = ContigId(u32::try_from(read.ref_id).expect("ref_id fits u32"));
        let pos = u32::try_from(read.pos).expect("read position fits u32");
        // Reads raw (un-canonicalised) bytes, matching production's
        // `RawContigRefCache` path so the ported filter behaves identically.
        // A zero span yields an empty slice → `read_exceeds_mismatch_fraction`
        // has no comparable bases and keeps the read (same as production).
        reference.fetch_raw_into(contig, pos, ref_span, ref_buf)?;
        if read_exceeds_mismatch_fraction(
            &read.cigar,
            &read.seq,
            &read.qual,
            ref_buf,
            config.mismatch_bq_floor.get(),
            max.get(),
        ) {
            return Ok(FilterVerdict::Drop(DropReason::HighMismatchFraction));
        }
    }

    Ok(FilterVerdict::Keep)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::FLAG_PAIRED;
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::pileup::walker::CigarOp;

    #[test]
    fn default_config_reproduces_the_production_filter_policy() {
        let config = ReadFilterConfig::default();
        assert_eq!(config.min_mapq, Some(MapQual(DEFAULT_MIN_MAPQ)));
        assert_eq!(config.min_read_length, Some(Bp(DEFAULT_MIN_READ_LENGTH)));
        assert!(config.drop_qc_fail);
        assert!(config.drop_duplicate);
        assert_eq!(
            config.max_read_mismatch_fraction.map(MismatchFraction::get),
            Some(DEFAULT_MAX_READ_MISMATCH_FRACTION)
        );
        assert_eq!(
            config.mismatch_bq_floor,
            BaseQual(DEFAULT_MISMATCH_BQ_FLOOR)
        );
    }

    #[test]
    fn counts_default_is_all_zero() {
        // Explicit all-zero literal (no `..`): pins every counter to 0 and forces
        // this test to be revisited if a counter field is ever added.
        assert_eq!(
            ReadFilterCounts::default(),
            ReadFilterCounts {
                kept: 0,
                duplicate: 0,
                low_mapq: 0,
                supplementary: 0,
                secondary: 0,
                unmapped: 0,
                qc_fail: 0,
                too_short: 0,
                high_mismatch_fraction: 0,
                bad_cigar: 0,
            }
        );
    }

    // ----- verdict_pre_decode (#1–#6) -------------------------------------

    /// A record's flag/MAPQ pair; the pre-decode cascade needs nothing else.
    fn pre(flag: u16, mapq: u8, config: &ReadFilterConfig) -> FilterVerdict {
        verdict_pre_decode(flag, MapQual(mapq), config)
    }

    #[test]
    fn pre_decode_keeps_a_clean_primary_read() {
        let cfg = ReadFilterConfig::default();
        assert_eq!(pre(0, 60, &cfg), FilterVerdict::Keep);
    }

    #[test]
    fn low_mapq_boundary_keeps_at_threshold_drops_one_below() {
        let cfg = ReadFilterConfig::default(); // min_mapq = Some(20)
        assert_eq!(pre(0, 20, &cfg), FilterVerdict::Keep);
        assert_eq!(pre(0, 19, &cfg), FilterVerdict::Drop(DropReason::LowMapq));
        // Unavailable MAPQ arrives as 0 (resolved by the record source) → dropped.
        assert_eq!(pre(0, 0, &cfg), FilterVerdict::Drop(DropReason::LowMapq));
    }

    #[test]
    fn no_mapq_minimum_keeps_any_quality() {
        let cfg = ReadFilterConfig {
            min_mapq: None,
            ..ReadFilterConfig::default()
        };
        assert_eq!(pre(0, 0, &cfg), FilterVerdict::Keep);
    }

    #[test]
    fn each_flag_bit_drops_to_its_own_bucket() {
        let cfg = ReadFilterConfig::default();
        assert_eq!(
            pre(FLAG_DUPLICATE, 60, &cfg),
            FilterVerdict::Drop(DropReason::Duplicate)
        );
        assert_eq!(
            pre(FLAG_SUPPLEMENTARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
        assert_eq!(
            pre(FLAG_SECONDARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Secondary)
        );
        assert_eq!(
            pre(FLAG_UNMAPPED, 60, &cfg),
            FilterVerdict::Drop(DropReason::Unmapped)
        );
        assert_eq!(
            pre(FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::QcFail)
        );
    }

    #[test]
    fn duplicate_and_qc_fail_toggles_off_keep_those_reads() {
        let cfg = ReadFilterConfig {
            drop_duplicate: false,
            drop_qc_fail: false,
            ..ReadFilterConfig::default()
        };
        assert_eq!(pre(FLAG_DUPLICATE, 60, &cfg), FilterVerdict::Keep);
        assert_eq!(pre(FLAG_QC_FAIL, 60, &cfg), FilterVerdict::Keep);
        // Supplementary/secondary/unmapped have no toggle — still dropped.
        assert_eq!(
            pre(FLAG_SUPPLEMENTARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
    }

    #[test]
    fn pre_decode_attribution_charges_the_first_firing_filter() {
        let cfg = ReadFilterConfig::default();
        // Duplicate (1) + unmapped (5) both set → charged to duplicate (earlier).
        assert_eq!(
            pre(FLAG_DUPLICATE | FLAG_UNMAPPED, 60, &cfg),
            FilterVerdict::Drop(DropReason::Duplicate)
        );
        // Supplementary (3) + secondary (4) → supplementary (earlier).
        assert_eq!(
            pre(FLAG_SUPPLEMENTARY | FLAG_SECONDARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
        // Low MAPQ (2) beats unmapped (5) when both apply.
        assert_eq!(
            pre(FLAG_UNMAPPED, 5, &cfg),
            FilterVerdict::Drop(DropReason::LowMapq)
        );
    }

    #[test]
    fn pre_decode_charges_filters_in_full_cascade_order() {
        // Every filter would fire; peel them off one at a time and confirm the
        // read is charged to each in the exact cascade order
        // (duplicate → low-MAPQ → supplementary → secondary → unmapped → QC-fail).
        let cfg = ReadFilterConfig::default();
        let all_flags =
            FLAG_DUPLICATE | FLAG_SUPPLEMENTARY | FLAG_SECONDARY | FLAG_UNMAPPED | FLAG_QC_FAIL;
        // mapq 0 also fails #2 throughout, so every stage after #1 has low-MAPQ
        // waiting behind it — proving each earlier filter really wins.
        assert_eq!(
            pre(all_flags, 0, &cfg),
            FilterVerdict::Drop(DropReason::Duplicate)
        );
        assert_eq!(
            pre(all_flags & !FLAG_DUPLICATE, 0, &cfg),
            FilterVerdict::Drop(DropReason::LowMapq)
        );
        // From here MAPQ is fine (60), so the next flag in order wins.
        assert_eq!(
            pre(all_flags & !FLAG_DUPLICATE, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
        assert_eq!(
            pre(FLAG_SECONDARY | FLAG_UNMAPPED | FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::Secondary)
        );
        assert_eq!(
            pre(FLAG_UNMAPPED | FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::Unmapped)
        );
        assert_eq!(
            pre(FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::QcFail)
        );
    }

    // ----- verdict_post_decode (#7, #9, #8) -------------------------------

    /// Build a mapped read at contig 0, pos 1, with the given decoded sequence,
    /// per-base qualities, and CIGAR. Flag/MAPQ are already-passed values.
    fn mapped(seq: &[u8], qual: &[u8], cigar: Vec<CigarOp>) -> MappedRead {
        MappedRead {
            qname: b"read".to_vec(),
            flag: FLAG_PAIRED,
            ref_id: 0,
            pos: 1,
            mapq: 60,
            cigar,
            seq: seq.to_vec(),
            qual: qual.to_vec(),
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            source_file_index: 0,
        }
    }

    /// A single contig of `n` adenines — the reference for the post-decode
    /// tests, where a `Match` read of `T`s mismatches every aligned base.
    fn poly_a_ref(n: usize) -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(vec![vec![b'A'; n]])
    }

    /// Post-decode config isolating one filter at a time: #7 off, toggles off,
    /// #8 set to `max` (`None` disables it).
    fn post_config(max: Option<f32>) -> ReadFilterConfig {
        ReadFilterConfig {
            min_mapq: None,
            min_read_length: None,
            drop_qc_fail: false,
            drop_duplicate: false,
            max_read_mismatch_fraction: max.map(|x| MismatchFraction::try_new(x).unwrap()),
            mismatch_bq_floor: BaseQual(0),
        }
    }

    fn post(
        read: &MappedRead,
        reference: &impl RawRefSeq,
        config: &ReadFilterConfig,
    ) -> FilterVerdict {
        let mut buf = Vec::new();
        verdict_post_decode(read, reference, config, &mut buf).unwrap()
    }

    #[test]
    fn too_short_boundary_keeps_at_threshold_drops_one_below() {
        let reference = poly_a_ref(40);
        let cfg = ReadFilterConfig {
            min_read_length: Some(Bp(30)),
            ..post_config(None) // #8 already disabled by post_config(None)
        };
        let at = mapped(&vec![b'A'; 30], &vec![40; 30], vec![CigarOp::Match(30)]);
        let below = mapped(&vec![b'A'; 29], &vec![40; 29], vec![CigarOp::Match(29)]);
        assert_eq!(post(&at, &reference, &cfg), FilterVerdict::Keep);
        assert_eq!(
            post(&below, &reference, &cfg),
            FilterVerdict::Drop(DropReason::TooShort)
        );
    }

    #[test]
    fn bad_cigar_drops_the_two_ill_formed_shapes() {
        let reference = poly_a_ref(40);
        let cfg = post_config(None);
        // Adjacent insertion/deletion pair.
        let adjacent_indel = mapped(
            b"AAAAAAAA",
            &[40; 8],
            vec![
                CigarOp::Match(4),
                CigarOp::Insertion(1),
                CigarOp::Deletion(1),
                CigarOp::Match(4),
            ],
        );
        // Leading deletion (boundary deletion).
        let boundary_deletion = mapped(
            b"AAAAAAAA",
            &[40; 8],
            vec![CigarOp::Deletion(2), CigarOp::Match(8)],
        );
        assert_eq!(
            post(&adjacent_indel, &reference, &cfg),
            FilterVerdict::Drop(DropReason::BadCigar)
        );
        assert_eq!(
            post(&boundary_deletion, &reference, &cfg),
            FilterVerdict::Drop(DropReason::BadCigar)
        );
    }

    #[test]
    fn high_mismatch_boundary_keeps_at_threshold_drops_above() {
        // ref = 10×A; a Match(10) read with k mismatches has fraction k/10.
        let reference = poly_a_ref(10);
        let cfg = post_config(Some(0.10));
        // 1/10 = 0.10, not > 0.10 → kept (boundary is exclusive).
        let at = mapped(b"TAAAAAAAAA", &[40; 10], vec![CigarOp::Match(10)]);
        // 2/10 = 0.20 > 0.10 → dropped.
        let above = mapped(b"TTAAAAAAAA", &[40; 10], vec![CigarOp::Match(10)]);
        assert_eq!(post(&at, &reference, &cfg), FilterVerdict::Keep);
        assert_eq!(
            post(&above, &reference, &cfg),
            FilterVerdict::Drop(DropReason::HighMismatchFraction)
        );
    }

    #[test]
    fn low_quality_mismatches_do_not_count_toward_the_fraction() {
        let reference = poly_a_ref(10);
        let cfg = ReadFilterConfig {
            mismatch_bq_floor: BaseQual(10),
            ..post_config(Some(0.0))
        };
        // Two mismatches, both below the BQ floor → neither counts → kept even
        // at a zero threshold.
        let read = mapped(
            b"TTAAAAAAAA",
            &[5, 5, 40, 40, 40, 40, 40, 40, 40, 40],
            vec![CigarOp::Match(10)],
        );
        assert_eq!(post(&read, &reference, &cfg), FilterVerdict::Keep);
    }

    #[test]
    fn mismatch_filter_disabled_makes_no_reference_access() {
        // Empty reference: any fetch errors with UnknownContig. With #8 disabled
        // the fetch never happens, so a high-mismatch read is kept.
        let empty = InMemoryRefSeq::from_contigs(Vec::new());
        let read = mapped(b"TTTTTTTTTT", &[40; 10], vec![CigarOp::Match(10)]);

        let disabled = post_config(None);
        assert_eq!(post(&read, &empty, &disabled), FilterVerdict::Keep);

        // Enabling #8 on the same empty reference proves a fetch is attempted.
        let enabled = post_config(Some(0.10));
        let mut buf = Vec::new();
        assert!(matches!(
            verdict_post_decode(&read, &empty, &enabled, &mut buf),
            Err(RefSeqError::UnknownContig(_))
        ));
    }

    #[test]
    fn high_mismatch_fetch_past_contig_end_is_fatal() {
        // A read whose reference span runs past the contig end (Match(10) at
        // pos 1 on a 5-base contig). Under the fatal error model the OutOfBounds
        // fetch propagates as `Err`, not a per-read drop or keep — a
        // validly-aligned read cannot cover positions the contig lacks, so this
        // signals a malformed record.
        let reference = poly_a_ref(5);
        let read = mapped(b"TTTTTTTTTT", &[40; 10], vec![CigarOp::Match(10)]);
        let mut buf = Vec::new();
        assert!(matches!(
            verdict_post_decode(&read, &reference, &post_config(Some(0.10)), &mut buf),
            Err(RefSeqError::OutOfBounds { .. })
        ));
    }

    #[test]
    fn bad_cigar_is_charged_before_high_mismatch() {
        // A read that is BOTH a boundary deletion (#9) and, on its M positions,
        // fully mismatched against the reference (#8 at a zero threshold).
        let reference = poly_a_ref(10);
        let both = mapped(
            b"TTTTTTTT",
            &[40; 8],
            vec![CigarOp::Deletion(2), CigarOp::Match(8)],
        );
        // #9 fires first → charged BadCigar, and #8's reference walk is skipped.
        assert_eq!(
            post(&both, &reference, &post_config(Some(0.0))),
            FilterVerdict::Drop(DropReason::BadCigar)
        );
        // The same sequence with a well-formed CIGAR does reach #8 and drops there,
        // proving the read would have failed #8 too — attribution is the only
        // difference the ordering makes.
        let good_cigar = mapped(b"TTTTTTTT", &[40; 8], vec![CigarOp::Match(8)]);
        assert_eq!(
            post(&good_cigar, &reference, &post_config(Some(0.0))),
            FilterVerdict::Drop(DropReason::HighMismatchFraction)
        );
    }

    #[test]
    fn too_short_is_charged_before_bad_cigar() {
        let reference = poly_a_ref(10);
        let cfg = ReadFilterConfig {
            min_read_length: Some(Bp(30)),
            ..post_config(None)
        };
        // Short AND a boundary deletion → charged TooShort (#7 before #9).
        let read = mapped(
            b"AAAAA",
            &[40; 5],
            vec![CigarOp::Deletion(2), CigarOp::Match(5)],
        );
        assert_eq!(
            post(&read, &reference, &cfg),
            FilterVerdict::Drop(DropReason::TooShort)
        );
    }

    #[test]
    fn zero_reference_span_read_is_kept_by_the_mismatch_filter() {
        // An all-soft-clip read has ref_span 0 → empty slice → no comparable
        // bases → kept (matches production's skip-on-no-span behaviour).
        let reference = poly_a_ref(10);
        let read = mapped(b"TTTT", &[40; 4], vec![CigarOp::SoftClip(4)]);
        assert_eq!(
            post(&read, &reference, &post_config(Some(0.0))),
            FilterVerdict::Keep
        );
    }
}
