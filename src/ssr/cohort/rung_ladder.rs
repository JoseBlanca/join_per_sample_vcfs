//! Rung-ladder & confident-genotype resolution — the shared peak primitive both
//! Stage-2 halves call (arch `ssr_call_parameters.md` §2, spec §5,
//! implementation plan §3.4, B1).
//!
//! Three things live here: [`build_rungs`] pools the cohort's observed sequences
//! into the **length-keyed rung ladder** (the cross-sample coordinate system §5
//! reads), and [`resolve_confident_genotype`] runs the **heuristic** confident-
//! genotype gate on one sample — 1..ploidy *clear local maxima* guarded by
//! separation, dosage balance, cohort peak-recurrence, and a depth floor. A
//! *confident genotype* (one clear peak = homozygote, up to ploidy = a separated
//! het) is the labelled seed the pre-pass estimates chemistry from.
//!
//! This is the **heuristic** resolution; the model-based 1-vs-2-peak BIC test (D1)
//! layers the C2 likelihood on top of these peaks. Pure-allele lengths are exact
//! multiples of the period; B1 keys rungs by `seq.len() / period` units (substitution
//! variants share a length, the in-tract-no-indel invariant of §6).

use std::collections::BTreeMap;

use crate::ssr::cohort::types::{CohortLocus, SampleEvidence};

/// A distinct observed sequence at a rung with its cohort support — the element a rung
/// holds.
///
/// `reads` is the total cohort read count; `samples` is the number of **distinct samples**
/// that observed it. The sample count is the recurrence signal the same-length admission
/// bar keys on (spec §5.2): a real interruption haplotype recurs across its carriers, a
/// per-base substitution error is sporadic and lands at a different base each time. It is
/// an integer count → order-free, so it does not disturb cross-thread determinism (arch §4).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct RungSeq {
    /// The distinct tract sequence.
    pub(crate) seq: Box<[u8]>,
    /// Total cohort reads supporting this sequence at its length.
    pub(crate) reads: u32,
    /// Number of distinct samples that observed this sequence (spec §5.2 recurrence).
    pub(crate) samples: u32,
}

/// One resolved peak with its labelled parent allele.
///
/// (Shape, not final — B1 may extend it.) Carries what the dosage-consistency and
/// cohort-recurrence checks need: the parent allele's tract sequence, its length in
/// repeat units, and the read support at the peak.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PeakAllele {
    /// The peak's parent allele (tract sequence).
    pub(crate) allele: Box<[u8]>,
    /// The allele's length in **repeat units**.
    pub(crate) repeat_len: u16,
    /// Read support at the peak.
    pub(crate) support: u32,
}

/// A confidently resolved genotype: 1..=ploidy clear peaks, each a labelled parent
/// allele (homozygote = the 1-peak case; a separated het = 2 peaks).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum ResolvedGenotype {
    /// The resolved peaks (length `1..=ploidy`).
    Peaks(Vec<PeakAllele>),
}

/// Why a (sample, locus) failed confident resolution and is left to the soft EM
/// rather than forced into a label (arch §2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum UnresolvedReason {
    /// Peaks were `< 2` units apart (a merged / masked het).
    Merged,
    /// Peak heights are inconsistent with any integer allele dosage (e.g. a
    /// homozygote-plus-heavy-stutter masquerading as a het).
    DosageInconsistent,
    /// A resolved allele does not recur as a clear peak in enough samples.
    NonRecurrent,
    /// Too little depth to resolve the peaks at all.
    Thin,
}

/// The outcome of confident-genotype resolution for one (sample, locus).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Resolution {
    /// Resolved at full confidence — a chemistry seed.
    Confident(ResolvedGenotype),
    /// Not resolved — left to the soft EM, with the reason recorded.
    Unresolved(UnresolvedReason),
}

/// Thresholds for rung admission + confident-genotype resolution.
///
/// All are deliberately exposed (no hidden defaults); the numbers are pinned on the
/// simulator in F2. [`dev_default`](Self::dev_default) gives the working values the
/// tests + the burn-in seed start from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct RungCfg {
    /// Prominence floor: a length is a *clear* local maximum only when its support
    /// stands **more than** this many reads above each adjacent (±1 unit) length
    /// (spec §5: "> 3 reads above each adjacent rung").
    pub(crate) prominence: u32,
    /// A resolved allele must recur as a *clear peak* in at least this many samples
    /// (the cohort-recurrence guard against a stutter artifact masquerading as an
    /// allele). Counts the sample itself.
    pub(crate) recurrence_k: u32,
    /// Resolved peaks must be at least this many repeat units apart (no merged /
    /// masked pair).
    pub(crate) separation_min: u16,
    /// A sample needs at least this total depth to attempt resolution.
    pub(crate) min_depth: u32,
    /// For a multi-peak genotype, the smallest peak's support must be at least this
    /// fraction of the largest's — the dosage-balance guard that rejects a
    /// homozygote-plus-heavy-stutter masquerading as a het.
    pub(crate) balance_ratio: f64,
}

impl RungCfg {
    /// Working defaults (recalibrated in F2).
    pub(crate) fn dev_default() -> Self {
        Self {
            prominence: 3,
            recurrence_k: 2,
            separation_min: 2,
            min_depth: 10,
            balance_ratio: 0.30,
        }
    }
}

/// The pooled cohort rung ladder at one locus: the occupied repeat-unit lengths and,
/// per length, the distinct sequences (with cohort counts) and how many samples
/// peak there.
///
/// "Length" is `seq.len() / period` repeat units. The ladder is the cross-sample
/// coordinate system S2 reads and candidate assembly (C1) unions over; B1 uses the
/// per-length **peak recurrence** for the resolution guard.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Rungs {
    /// The locus motif period (bytes per repeat unit).
    period: usize,
    /// Ascending occupied lengths (repeat units) — the rung coordinate system.
    lengths: Vec<u16>,
    /// Per length: total reads observed across the cohort.
    cohort_support: BTreeMap<u16, u32>,
    /// Per length: number of samples in which that length is a clear local maximum.
    peak_recurrence: BTreeMap<u16, u32>,
    /// Per length: distinct sequences at that length with their cohort support
    /// (a rung holds a *set* of sequences — substitution / interruption variants).
    seqs_by_length: BTreeMap<u16, Vec<RungSeq>>,
}

impl Rungs {
    /// The ascending occupied rung lengths (repeat units).
    pub(crate) fn lengths(&self) -> &[u16] {
        &self.lengths
    }

    /// Total cohort reads observed at `length` units.
    pub(crate) fn cohort_support(&self, length: u16) -> u32 {
        self.cohort_support.get(&length).copied().unwrap_or(0)
    }

    /// How many samples have a clear local maximum at `length` units.
    pub(crate) fn peak_recurrence(&self, length: u16) -> u32 {
        self.peak_recurrence.get(&length).copied().unwrap_or(0)
    }

    /// The distinct sequences (with cohort support) at `length` units, byte-sorted, if any.
    pub(crate) fn seqs_at(&self, length: u16) -> Option<&[RungSeq]> {
        self.seqs_by_length.get(&length).map(Vec::as_slice)
    }

    /// The locus motif period.
    pub(crate) fn period(&self) -> usize {
        self.period
    }

    /// The cohort modal allele length (the most-supported rung; ties → the shorter
    /// length). `None` only for an empty locus. This is the centre of the `G₀` prior
    /// (spec §5.5 — the mode, *not* the reference).
    pub(crate) fn modal_length(&self) -> Option<u16> {
        self.cohort_support
            .iter()
            .max_by(|(la, sa), (lb, sb)| sa.cmp(sb).then_with(|| lb.cmp(la)))
            .map(|(length, _)| *length)
    }
}

/// One sample's length histogram (repeat units → total supporting reads).
fn sample_histogram(evidence: &SampleEvidence, period: usize) -> BTreeMap<u16, u32> {
    let mut histogram = BTreeMap::new();
    for (seq, count) in &evidence.seq_counts {
        let units = (seq.len() / period) as u16;
        *histogram.entry(units).or_insert(0) += count;
    }
    histogram
}

/// Whether `length` is a *clear* local maximum in `histogram`: its support exceeds
/// each adjacent (±1 unit) length's support by more than `prominence` reads (a
/// missing neighbour counts as zero support).
fn is_clear_peak(histogram: &BTreeMap<u16, u32>, length: u16, prominence: u32) -> bool {
    let support = histogram.get(&length).copied().unwrap_or(0);
    if support == 0 {
        return false;
    }
    let lower = length
        .checked_sub(1)
        .and_then(|l| histogram.get(&l).copied())
        .unwrap_or(0);
    let upper = length
        .checked_add(1)
        .and_then(|l| histogram.get(&l).copied())
        .unwrap_or(0);
    support > lower + prominence && support > upper + prominence
}

/// Pool the cohort's observed sequences into the rung ladder (spec §5 levels 1–2).
pub(crate) fn build_rungs(locus: &CohortLocus, cfg: &RungCfg) -> Rungs {
    let period = locus.motif.period().max(1);
    let mut cohort_support: BTreeMap<u16, u32> = BTreeMap::new();
    let mut peak_recurrence: BTreeMap<u16, u32> = BTreeMap::new();
    let mut seqs_by_length: BTreeMap<u16, Vec<RungSeq>> = BTreeMap::new();

    for evidence in &locus.samples {
        // Each `(seq, count)` in a sample's `seq_counts` is one distinct sample observing
        // that sequence (the Stage-1 contract keys `seq_counts` on distinct sequences), so
        // every match/insert here increments the per-sequence distinct-sample tally by one.
        for (seq, count) in &evidence.seq_counts {
            let units = (seq.len() / period) as u16;
            *cohort_support.entry(units).or_insert(0) += count;
            let bucket = seqs_by_length.entry(units).or_default();
            match bucket.iter_mut().find(|rs| rs.seq.as_ref() == seq.as_ref()) {
                Some(rs) => {
                    rs.reads += count;
                    rs.samples += 1;
                }
                None => bucket.push(RungSeq {
                    seq: seq.clone(),
                    reads: *count,
                    samples: 1,
                }),
            }
        }
        let histogram = sample_histogram(evidence, period);
        for &length in histogram.keys() {
            if is_clear_peak(&histogram, length, cfg.prominence) {
                *peak_recurrence.entry(length).or_insert(0) += 1;
            }
        }
    }

    // Deterministic, distinct-sequence order within each rung (bytes ascending).
    for bucket in seqs_by_length.values_mut() {
        bucket.sort_by(|a, b| a.seq.cmp(&b.seq));
    }

    let lengths = cohort_support.keys().copied().collect();
    Rungs {
        period,
        lengths,
        cohort_support,
        peak_recurrence,
        seqs_by_length,
    }
}

/// One sample's clear local maxima as labelled peaks (shared by the confident-
/// genotype gate and candidate assembly, C1). Ascending by length.
pub(crate) fn sample_clear_peaks(
    evidence: &SampleEvidence,
    period: usize,
    prominence: u32,
) -> Vec<PeakAllele> {
    let histogram = sample_histogram(evidence, period);
    histogram
        .keys()
        .copied()
        .filter(|&length| is_clear_peak(&histogram, length, prominence))
        .map(|length| PeakAllele {
            allele: representative_sequence(evidence, length, period),
            repeat_len: length,
            support: histogram[&length],
        })
        .collect()
}

/// The sample's most-supported sequence at `length` units (its representative
/// parent allele for a resolved peak). `length` is known to be occupied.
fn representative_sequence(evidence: &SampleEvidence, length: u16, period: usize) -> Box<[u8]> {
    evidence
        .seq_counts
        .iter()
        .filter(|(seq, _)| (seq.len() / period) as u16 == length)
        // Pick the most-supported sequence; break count ties on the lexicographically
        // smallest sequence (so `max_by` returns it, the tie-break order is reversed).
        .max_by(|(a_seq, a_c), (b_seq, b_c)| a_c.cmp(b_c).then_with(|| b_seq.cmp(a_seq)))
        .map(|(seq, _)| seq.clone())
        .expect("a resolved peak length is occupied by the sample")
}

/// Run the heuristic confident-genotype gate on one sample (spec §5 level 4 + the
/// arch §2 guards). Returns the labelled peaks when confident, else the reason.
pub(crate) fn resolve_confident_genotype(
    sample: &SampleEvidence,
    rungs: &Rungs,
    ploidy: u8,
    cfg: &RungCfg,
) -> Resolution {
    let depth: u32 = sample.seq_counts.iter().map(|(_, c)| c).sum();
    if depth < cfg.min_depth {
        return Resolution::Unresolved(UnresolvedReason::Thin);
    }

    let mut peaks = sample_clear_peaks(sample, rungs.period, cfg.prominence);

    // No clear maximum despite adequate depth: two adjacent alleles cancel each
    // other's prominence (the 1-apart merged het). NOTE: a structureless / pure-noise
    // locus also lands here — both correctly mean "not a seed", so they share the
    // `Merged` reason in v1; split out a `NoStructure` reason only if D1's diagnostics
    // need to tell them apart.
    if peaks.is_empty() {
        return Resolution::Unresolved(UnresolvedReason::Merged);
    }

    // Keep the top-ploidy peaks (ties broken by the shorter allele for determinism),
    // then order them by length for a stable, dosage-readable genotype.
    peaks.sort_by(|a, b| {
        b.support
            .cmp(&a.support)
            .then_with(|| a.repeat_len.cmp(&b.repeat_len))
    });
    peaks.truncate(ploidy as usize);
    peaks.sort_by_key(|p| p.repeat_len);

    // Separation: peaks pairwise ≥ separation_min units apart.
    if peaks
        .windows(2)
        .any(|w| w[1].repeat_len - w[0].repeat_len < cfg.separation_min)
    {
        return Resolution::Unresolved(UnresolvedReason::Merged);
    }

    // Dosage balance: reject a homozygote-plus-heavy-stutter masquerading as a het.
    if peaks.len() >= 2 {
        let max = peaks.iter().map(|p| p.support).max().unwrap();
        let min = peaks.iter().map(|p| p.support).min().unwrap();
        if (min as f64) < cfg.balance_ratio * max as f64 {
            return Resolution::Unresolved(UnresolvedReason::DosageInconsistent);
        }
    }

    // Cohort recurrence: each resolved allele recurs as a clear peak in ≥ k samples.
    if peaks
        .iter()
        .any(|p| rungs.peak_recurrence(p.repeat_len) < cfg.recurrence_k)
    {
        return Resolution::Unresolved(UnresolvedReason::NonRecurrent);
    }

    Resolution::Confident(ResolvedGenotype::Peaks(peaks))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn peak(allele: &[u8], repeat_len: u16, support: u32) -> PeakAllele {
        PeakAllele {
            allele: Box::from(allele),
            repeat_len,
            support,
        }
    }

    #[test]
    fn confident_homozygote_is_one_peak() {
        let res = Resolution::Confident(ResolvedGenotype::Peaks(vec![peak(b"ATATAT", 3, 40)]));
        match res {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(peaks.len(), 1);
                assert_eq!(peaks[0].repeat_len, 3);
            }
            other => panic!("expected a confident 1-peak genotype, got {other:?}"),
        }
    }

    #[test]
    fn confident_separated_het_carries_two_labelled_peaks() {
        let res = Resolution::Confident(ResolvedGenotype::Peaks(vec![
            peak(b"ATAT", 2, 22),
            peak(b"ATATATATAT", 5, 19),
        ]));
        match res {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(peaks.len(), 2);
                assert_eq!(peaks[0].repeat_len, 2);
                assert_eq!(peaks[1].repeat_len, 5);
            }
            other => panic!("expected a confident 2-peak genotype, got {other:?}"),
        }
    }

    #[test]
    fn unresolved_reasons_are_distinct() {
        assert_ne!(
            Resolution::Unresolved(UnresolvedReason::Merged),
            Resolution::Unresolved(UnresolvedReason::Thin)
        );
    }

    // ── B1: build_rungs + resolve_confident_genotype ──

    use crate::ssr::cohort::types::{CohortLocus, LocusId, SampleEvidence, SsrQc};
    use crate::ssr::types::Motif;

    /// A `CA`-tiled tract of `units` repeat units (12 bp for `units = 6`).
    fn ca_seq(units: u16) -> Box<[u8]> {
        std::iter::repeat_n(*b"CA", units as usize)
            .flatten()
            .collect()
    }

    /// A `SampleEvidence` over a CA (period-2) motif from `(units, count)` bins —
    /// each bin becomes the tiled sequence `CA × units`. Sorted ascending by bytes
    /// (the Stage-1 contract).
    fn ca_evidence(bins: &[(u16, u32)]) -> SampleEvidence {
        let mut seq_counts: Vec<(Box<[u8]>, u32)> = bins
            .iter()
            .map(|&(units, count)| {
                let mut seq = Vec::new();
                for _ in 0..units {
                    seq.extend_from_slice(b"CA");
                }
                (seq.into_boxed_slice(), count)
            })
            .collect();
        seq_counts.sort_by(|(a, _), (b, _)| a.cmp(b));
        SampleEvidence {
            seq_counts,
            qc: SsrQc::default(),
        }
    }

    fn ca_cohort(samples: Vec<SampleEvidence>) -> CohortLocus {
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 20,
                end: 32,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGGCACACATTTTTT".as_slice()),
            Box::from(b"CACACA".as_slice()),
        );
        for (idx, evidence) in samples.into_iter().enumerate() {
            locus.push(idx as u32, evidence);
        }
        locus
    }

    fn resolve_sample0(locus: &CohortLocus, cfg: &RungCfg) -> Resolution {
        let rungs = build_rungs(locus, cfg);
        resolve_confident_genotype(&locus.samples[0], &rungs, 2, cfg)
    }

    #[test]
    fn build_rungs_pools_lengths_support_and_peak_recurrence() {
        // Two homozygous-6 samples (faithful peak at 6, ±1 stutter skirt).
        let sample = ca_evidence(&[(5, 5), (6, 50), (7, 5)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        let rungs = build_rungs(&locus, &RungCfg::dev_default());

        assert_eq!(rungs.lengths(), &[5, 6, 7]);
        assert_eq!(rungs.cohort_support(6), 100); // 50 + 50
        assert_eq!(rungs.peak_recurrence(6), 2); // both samples peak at 6
        assert_eq!(rungs.peak_recurrence(5), 0); // stutter band, never a peak
        let modal = rungs.seqs_at(6).unwrap();
        assert_eq!(modal.len(), 1);
        assert_eq!(modal[0].reads, 100); // 50 + 50 cohort reads
        assert_eq!(modal[0].samples, 2); // both samples observed the length-6 sequence
    }

    #[test]
    fn build_rungs_tallies_distinct_samples_per_same_length_sequence() {
        // Two same-length (6-unit / 12 bp) sequences: a pure `CA×6` and an interrupted
        // variant of equal length. The pure allele recurs in 3 samples, the interrupted in
        // 2, and a singleton substitution noise variant in 1 — the per-sequence sample
        // tally must separate them (the §5.2 recurrence signal).
        let pure = ca_seq(6); // CA×6, 12 bp
        let interrupted: Box<[u8]> = {
            let mut s = pure.to_vec();
            s[5] = b'T'; // flip one interior base → same length, distinct sequence
            s.into_boxed_slice()
        };
        let noise: Box<[u8]> = {
            let mut s = pure.to_vec();
            s[7] = b'G'; // a different sporadic substitution
            s.into_boxed_slice()
        };
        let sample_of = |seqs: &[&Box<[u8]>]| {
            let mut seq_counts: Vec<(Box<[u8]>, u32)> =
                seqs.iter().map(|s| ((*s).clone(), 20)).collect();
            seq_counts.sort_by(|(a, _), (b, _)| a.cmp(b));
            SampleEvidence {
                seq_counts,
                qc: SsrQc::default(),
            }
        };
        let samples = vec![
            sample_of(&[&pure, &interrupted]),
            sample_of(&[&pure, &interrupted]),
            sample_of(&[&pure, &noise]),
        ];
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: 0,
                start: 20,
                end: 32,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGGCACACACACACATTTTTT".as_slice()),
            pure.clone(),
        );
        for (idx, ev) in samples.into_iter().enumerate() {
            locus.push(idx as u32, ev);
        }
        let rungs = build_rungs(&locus, &RungCfg::dev_default());

        let at6 = rungs.seqs_at(6).unwrap();
        let tally = |seq: &[u8]| at6.iter().find(|rs| rs.seq.as_ref() == seq).unwrap();
        assert_eq!(tally(&pure).samples, 3);
        assert_eq!(tally(&interrupted).samples, 2);
        assert_eq!(tally(&noise).samples, 1);
        // Bytes-ascending order within the rung (determinism contract).
        assert!(at6.windows(2).all(|w| w[0].seq <= w[1].seq));
    }

    #[test]
    fn resolve_confident_genotype_calls_a_clean_homozygote() {
        let sample = ca_evidence(&[(4, 1), (5, 5), (6, 50), (7, 5), (8, 1)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        match resolve_sample0(&locus, &RungCfg::dev_default()) {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(peaks.len(), 1);
                assert_eq!(peaks[0].repeat_len, 6);
                assert_eq!(&*peaks[0].allele, b"CACACACACACA"); // CA × 6
            }
            other => panic!("expected a confident homozygote, got {other:?}"),
        }
    }

    #[test]
    fn resolve_confident_genotype_calls_a_separated_het_with_two_labelled_peaks() {
        let sample = ca_evidence(&[(3, 1), (4, 30), (5, 4), (8, 4), (9, 28), (10, 1)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        match resolve_sample0(&locus, &RungCfg::dev_default()) {
            Resolution::Confident(ResolvedGenotype::Peaks(peaks)) => {
                assert_eq!(peaks.len(), 2);
                assert_eq!(peaks[0].repeat_len, 4);
                assert_eq!(peaks[1].repeat_len, 9);
            }
            other => panic!("expected a confident separated het, got {other:?}"),
        }
    }

    #[test]
    fn resolve_confident_genotype_rejects_a_one_apart_merged_het() {
        // Alleles 5 and 6 units (1 apart, balanced): neither clears the other's
        // prominence, so there is no clear peak — the masked het.
        let sample = ca_evidence(&[(4, 5), (5, 30), (6, 28), (7, 5)]);
        let locus = ca_cohort(vec![sample.clone(), sample]);
        assert_eq!(
            resolve_sample0(&locus, &RungCfg::dev_default()),
            Resolution::Unresolved(UnresolvedReason::Merged)
        );
    }

    #[test]
    fn resolve_confident_genotype_rejects_homozygote_plus_heavy_stutter_on_dosage() {
        // A tall peak at 8 with a −2 stutter satellite at 6 that is itself a local
        // maximum: two peaks, but wildly unbalanced (12 vs 50) → dosage-inconsistent.
        let sample = ca_evidence(&[(5, 1), (6, 12), (7, 4), (8, 50), (9, 4)]);
        let locus = ca_cohort(vec![sample]);
        assert_eq!(
            resolve_sample0(&locus, &RungCfg::dev_default()),
            Resolution::Unresolved(UnresolvedReason::DosageInconsistent)
        );
    }

    #[test]
    fn resolve_confident_genotype_rejects_a_non_recurrent_allele() {
        // A clean homozygote, but the allele peaks in only this one sample, so it
        // fails the recurrence_k = 2 guard.
        let sample = ca_evidence(&[(5, 5), (6, 50), (7, 5)]);
        let locus = ca_cohort(vec![sample]);
        assert_eq!(
            resolve_sample0(&locus, &RungCfg::dev_default()),
            Resolution::Unresolved(UnresolvedReason::NonRecurrent)
        );
    }

    #[test]
    fn resolve_confident_genotype_rejects_a_thin_sample() {
        let sample = ca_evidence(&[(6, 3)]); // depth 3 < min_depth 10
        let locus = ca_cohort(vec![sample.clone(), sample]);
        assert_eq!(
            resolve_sample0(&locus, &RungCfg::dev_default()),
            Resolution::Unresolved(UnresolvedReason::Thin)
        );
    }
}
