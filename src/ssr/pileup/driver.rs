//! The Stage-1 `ssr-pileup` driver — turns the sorted catalog + one sample's
//! alignment files into a per-locus `.ssr.psp` evidence file.
//!
//! The work for **one** locus is a single self-contained unit, [`process_locus`]
//! — fetch the locus's reads, analyze each, fold them into the locus record. It
//! is pure over its shared `&` inputs (and [`AlignmentFile`] is `Sync`), so the
//! eventual parallel step is just `for` → `par_iter().map_init(…)`: each worker
//! reuses its own [`LocusScratch`], nothing inside the unit changes (arch §8.4).
//!
//! **Build status (increment 2a):** the per-locus unit — [`process_locus`],
//! [`LocusScratch`], and the [`QcCounts`] assembly ([`qc_counts`]). The catalog
//! walk / `run()` loop, the name→chrom_id container adapter, the writer-header
//! build, and the CLI land in the following increments.

use crate::bam::errors::AlignmentInputError;
use crate::bam::segment_reader::AlignmentFile;
use crate::ssr::types::Locus;

use super::candidate_generation::CandidateAllele;
use super::fetch_reads::{LocusReads, fetch_locus_reads};
use super::locus_record::{QcCounts, SsrLocusRecord, aggregate};
use super::pair_hmm::{HmmModel, PairHmmScratch};
use super::read_analysis::analyze_read;

/// Errors from the Stage-1 driver. More variants (catalog read, container
/// write, header build) land with the `run()` loop in the next increment.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub(crate) enum SsrPileupError {
    #[error(transparent)]
    Read(#[from] AlignmentInputError),
}

/// Per-locus scratch reused across loci to avoid allocation churn — the
/// candidate-allele buffer and the pair-HMM forward-matrix workspace
/// [`analyze_read`] writes through. **One per worker thread:** the serial driver
/// keeps a single instance; the parallel driver hands each thread its own via
/// `rayon`'s `map_init` (so the reuse survives parallelization).
pub(crate) struct LocusScratch {
    cands: Vec<CandidateAllele>,
    hmm: PairHmmScratch,
}

impl LocusScratch {
    pub(crate) fn new() -> Self {
        Self {
            cands: Vec::new(),
            hmm: PairHmmScratch::new(),
        }
    }
}

/// Assemble the locus's QC scalars from the fetch pass. All three describe
/// **independent primary alignments** at the locus:
///
/// - `depth` — usable primary reads considered (passed the reader's cheap
///   filter and overlapped the window; uncapped, pre reach-gate/reservoir).
/// - `n_filtered` — primary reads the gate dropped: quality (QC-fail / low-MAPQ
///   / too-short) **plus** duplicates (auditable as filtered).
/// - `mapped_reads` — the dup-free primary coverage denominator (`depth` + the
///   quality drops, *excluding* duplicates).
///
/// Secondary / supplementary (non-independent) and unmapped (not at the locus)
/// reads are excluded everywhere. So `mapped_reads ≥ depth`, and `n_filtered`
/// adds duplicates on top — it is deliberately *not* `mapped_reads − depth`.
/// (Decided in `ssr_pileup_driver.md` §4/§7.3.)
pub(crate) fn qc_counts(fetched: &LocusReads) -> QcCounts {
    let d = &fetched.filtered;
    let primary_quality_drops = d.qc_fail + d.low_mapq + d.too_short;
    QcCounts {
        depth: fetched.yielded as u32,
        n_filtered: (primary_quality_drops + d.duplicate) as u32,
        mapped_reads: (fetched.yielded + primary_quality_drops) as u32,
    }
}

/// Process one locus end to end: fetch + depth-cap its reads, analyze each
/// (triage → realign spanning reads → `Qᵣ`), and fold the outcomes + QC into the
/// in-memory (chrom-**name**-keyed) [`SsrLocusRecord`]. The name→chrom_id
/// container adapter is applied by the writer stage (next increment).
///
/// `window` is the `analyze_read` candidate half-width; `cap` is the per-locus
/// reservoir depth cap. Concurrency-safe for different loci: the only `&mut` is
/// the caller-owned `scratch`.
pub(crate) fn process_locus(
    files: &[AlignmentFile],
    locus: &Locus,
    window: u16,
    cap: usize,
    model: &HmmModel,
    scratch: &mut LocusScratch,
) -> Result<SsrLocusRecord, SsrPileupError> {
    let fetched = fetch_locus_reads(files, locus, cap)?;
    let outcomes: Vec<_> = fetched
        .reads
        .iter()
        .map(|read| {
            analyze_read(
                read,
                locus,
                window,
                model,
                &mut scratch.cands,
                &mut scratch.hmm,
            )
        })
        .collect();
    Ok(aggregate(locus, &outcomes, qc_counts(&fetched)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::FilterCounts;

    fn locus_reads(yielded: u64, filtered: FilterCounts) -> LocusReads {
        LocusReads {
            reads: Vec::new(),
            yielded,
            filtered,
        }
    }

    #[test]
    fn qc_counts_excludes_dups_from_coverage_but_keeps_them_filtered() {
        // 10 usable; quality drops 2+1+1=4; 3 dups; plus non-independent /
        // not-at-locus reads that must be ignored everywhere.
        let filtered = FilterCounts {
            qc_fail: 2,
            low_mapq: 1,
            too_short: 1,
            duplicate: 3,
            secondary: 5,
            supplementary: 4,
            unmapped: 2,
            ..FilterCounts::default()
        };
        let qc = qc_counts(&locus_reads(10, filtered));

        assert_eq!(qc.depth, 10); // yielded
        assert_eq!(qc.n_filtered, 4 + 3); // quality drops + duplicates
        assert_eq!(qc.mapped_reads, 10 + 4); // depth + quality drops (dup-free)
        // Secondary / supplementary / unmapped influenced nothing.
    }

    #[test]
    fn qc_counts_coverage_is_at_least_depth_and_dup_free() {
        // Only duplicates filtered → coverage == depth (dups don't count), but
        // n_filtered still reports them.
        let filtered = FilterCounts {
            duplicate: 7,
            ..FilterCounts::default()
        };
        let qc = qc_counts(&locus_reads(20, filtered));

        assert_eq!(qc.depth, 20);
        assert_eq!(qc.mapped_reads, 20); // dups excluded from coverage
        assert_eq!(qc.n_filtered, 7); // but visible as filtered
    }

    #[test]
    fn qc_counts_all_zero_for_an_empty_locus() {
        let qc = qc_counts(&locus_reads(0, FilterCounts::default()));
        assert_eq!((qc.depth, qc.n_filtered, qc.mapped_reads), (0, 0, 0));
    }
}
