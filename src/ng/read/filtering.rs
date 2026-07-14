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
//! This milestone (A) lands the step-1-local types only; the two-phase cascade
//! (B), the record-source seam (C), and the `ReadFilter` iterator (D) follow.

use crate::bam::alignment_input::{
    DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH,
    DEFAULT_MISMATCH_BQ_FLOOR,
};
use crate::ng::types::{BaseQual, Bp, MapQual, MismatchFraction};

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

#[cfg(test)]
mod tests {
    use super::*;

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
}
