//! Per-locus tally — fold a locus's per-read delimitation outcomes into its
//! observed evidence (model Q3, arch `ssr_pileup_mark2.md` §5).
//!
//! Replaces the Mark-1 `locus_record` aggregator: no rung profiles and no
//! renormalization — just the **observed repeat-region sequences + their counts**
//! (the per-sample observed ladder), plus the QC scalars. The counts carry the
//! signal; quality was already gated upstream.

use std::collections::HashMap;

use crate::ssr::types::Locus;

/// The QC counts the fetcher tallies over its full pass — the reads triage never
/// classifies (dropped by the admission gate) — which the per-read outcomes
/// cannot recover on their own (arch §3.1/§8.3). The Mark-2 dropout counts
/// (`n_low_quality`, `n_border_off_end`) are derived from the outcomes instead,
/// in [`tally`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct QcCounts {
    /// Usable primary reads considered at the locus (the reach gate's pool).
    pub(crate) depth: u32,
    /// Reads dropped by the admission gate (low-MAPQ / duplicate / qc-fail / short).
    pub(crate) n_filtered: u32,
    /// Dup-free primary coverage denominator.
    pub(crate) mapped_reads: u32,
}

/// One read's delimitation outcome (the Mark-2 analogue of Mark-1 `ReadOutcome`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum ReadObs {
    /// A usable, quality-passing repeat-region sequence.
    Sequence(Box<[u8]>),
    /// A usable, quality-passing sequence that was only recovered by re-delimiting in a
    /// flank-widened window — the mapper had likely mis-aligned a long allele as all-Match
    /// (plan `ssr_pileup_long_allele_window_recovery.md`). Tallied as an observed sequence
    /// **and** counted (`n_widened`) so the mapper trouble is visible in QC.
    WidenedSequence(Box<[u8]>),
    /// Dropped by the first-quartile quality gate.
    LowQuality,
    /// A flank ran off the read end (allele ≥ read length).
    BorderOffEnd,
    /// The tract was still flank-truncated after widening — the allele exceeds what the
    /// read (or the safe `flank_bp` margin) can hold. Dropped from the tally so a collapsed
    /// allele can never reach the VCF, and counted (`n_window_truncated`) for QC.
    WindowTruncated,
}

/// One sample's observed evidence at one locus — in-memory form (chrom-**name**
/// keyed, 0-based). The driver's adapter maps it to the container's id-keyed,
/// 1-based record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SsrLocusObs {
    pub(crate) chrom: Box<str>,
    pub(crate) start: u32,
    pub(crate) end: u32,
    // QC scalars (the lean set; nothing derivable, model P3).
    pub(crate) depth: u32,
    pub(crate) n_filtered: u32,
    pub(crate) mapped_reads: u32,
    pub(crate) n_low_quality: u32,
    pub(crate) n_border_off_end: u32,
    /// Reads whose long allele was recovered by widening the realignment window (kept in
    /// `observed`; a mapper-misalignment QC signal).
    pub(crate) n_widened: u32,
    /// Reads dropped because the allele stayed flank-truncated even after widening (kept
    /// out of `observed`; a stronger mapper-misalignment / too-long-allele QC signal).
    pub(crate) n_window_truncated: u32,
    /// Distinct repeat-region sequences → observation count, **sorted by bytes**
    /// (deterministic storage order — the cross-thread byte-identity invariant).
    pub(crate) observed: Vec<(Box<[u8]>, u32)>,
}

/// Fold a locus's per-read outcomes + the fetch-pass QC into its record: tally
/// the sequences, count the two dropout classes, and emit `observed` sorted by
/// bytes. The tally is order-independent (commutative counts + a final sort), so
/// the record is identical regardless of the order reads were processed.
pub(crate) fn tally(locus: &Locus, outcomes: &[ReadObs], qc: QcCounts) -> SsrLocusObs {
    let mut counts: HashMap<Box<[u8]>, u32> = HashMap::new();
    let mut n_low_quality = 0u32;
    let mut n_border_off_end = 0u32;
    let mut n_widened = 0u32;
    let mut n_window_truncated = 0u32;
    for outcome in outcomes {
        match outcome {
            ReadObs::Sequence(seq) => *counts.entry(seq.clone()).or_insert(0) += 1,
            ReadObs::WidenedSequence(seq) => {
                *counts.entry(seq.clone()).or_insert(0) += 1;
                n_widened += 1;
            }
            ReadObs::LowQuality => n_low_quality += 1,
            ReadObs::BorderOffEnd => n_border_off_end += 1,
            ReadObs::WindowTruncated => n_window_truncated += 1,
        }
    }
    let mut observed: Vec<(Box<[u8]>, u32)> = counts.into_iter().collect();
    observed.sort_unstable_by(|a, b| a.0.cmp(&b.0));

    SsrLocusObs {
        chrom: locus.chrom().into(),
        start: locus.start(),
        end: locus.end(),
        depth: qc.depth,
        n_filtered: qc.n_filtered,
        mapped_reads: qc.mapped_reads,
        n_low_quality,
        n_border_off_end,
        n_widened,
        n_window_truncated,
        observed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::types::Motif;

    fn ca_locus() -> Locus {
        Locus::new(
            "chr1".into(),
            13,
            19,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGCACACATTT").into(),
            10,
        )
        .unwrap()
    }

    fn seq(bytes: &[u8]) -> ReadObs {
        ReadObs::Sequence(bytes.to_vec().into_boxed_slice())
    }

    fn qc() -> QcCounts {
        QcCounts {
            depth: 30,
            n_filtered: 4,
            mapped_reads: 33,
        }
    }

    #[test]
    fn tallies_sequences_and_counts_dropouts() {
        let outcomes = vec![
            seq(b"CACACA"),
            seq(b"CACACA"),
            seq(b"CACACACA"),
            ReadObs::LowQuality,
            ReadObs::BorderOffEnd,
            ReadObs::BorderOffEnd,
        ];
        let rec = tally(&ca_locus(), &outcomes, qc());
        assert_eq!(rec.n_low_quality, 1);
        assert_eq!(rec.n_border_off_end, 2);
        // Two distinct sequences, sorted by bytes ("CACACA" < "CACACACA").
        assert_eq!(
            rec.observed,
            vec![
                (b"CACACA".to_vec().into_boxed_slice(), 2),
                (b"CACACACA".to_vec().into_boxed_slice(), 1),
            ]
        );
    }

    #[test]
    fn copies_coords_and_qc_scalars() {
        let rec = tally(&ca_locus(), &[], qc());
        assert_eq!((&*rec.chrom, rec.start, rec.end), ("chr1", 13, 19));
        assert_eq!((rec.depth, rec.n_filtered, rec.mapped_reads), (30, 4, 33));
        assert!(rec.observed.is_empty());
        assert_eq!((rec.n_low_quality, rec.n_border_off_end), (0, 0));
        assert_eq!((rec.n_widened, rec.n_window_truncated), (0, 0));
    }

    #[test]
    fn widened_sequences_are_observed_and_counted_truncated_are_dropped_and_counted() {
        let outcomes = vec![
            seq(b"CACACA"),
            ReadObs::WidenedSequence(b"CACACACACACA".to_vec().into_boxed_slice()),
            ReadObs::WidenedSequence(b"CACACACACACA".to_vec().into_boxed_slice()),
            ReadObs::WindowTruncated,
        ];
        let rec = tally(&ca_locus(), &outcomes, qc());
        // The widened tracts join `observed` (kept) and bump n_widened; the truncated read
        // is absent from `observed` and only bumps n_window_truncated.
        assert_eq!(rec.n_widened, 2);
        assert_eq!(rec.n_window_truncated, 1);
        assert_eq!(
            rec.observed,
            vec![
                (b"CACACA".to_vec().into_boxed_slice(), 1),
                (b"CACACACACACA".to_vec().into_boxed_slice(), 2),
            ]
        );
    }

    #[test]
    fn observed_is_sorted_by_bytes() {
        let outcomes = vec![seq(b"GG"), seq(b"AA"), seq(b"CC"), seq(b"AA")];
        let rec = tally(&ca_locus(), &outcomes, qc());
        let keys: Vec<&[u8]> = rec.observed.iter().map(|(s, _)| s.as_ref()).collect();
        assert_eq!(keys, vec![b"AA".as_ref(), b"CC".as_ref(), b"GG".as_ref()]);
    }

    #[test]
    fn tally_is_order_independent() {
        // The same multiset of outcomes in two different orders → identical record
        // (counts are commutative and `observed` is sorted).
        let a = vec![
            seq(b"CACACA"),
            seq(b"CACA"),
            seq(b"CACACA"),
            ReadObs::LowQuality,
        ];
        let b = vec![
            ReadObs::LowQuality,
            seq(b"CACACA"),
            seq(b"CACACA"),
            seq(b"CACA"),
        ];
        assert_eq!(tally(&ca_locus(), &a, qc()), tally(&ca_locus(), &b, qc()));
    }
}
