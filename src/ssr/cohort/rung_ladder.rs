//! Rung-ladder & confident-genotype types — the A1 nouns of the shared peak
//! primitive (arch `ssr_call_parameters.md` §2, implementation plan §3.4).
//!
//! The `build_rungs` / `resolve_confident_genotype` logic is written in B1; this
//! file defines only the shapes the resolution produces. A *confident genotype* is
//! a (sample, locus) whose read-length distribution resolves into 1..ploidy clear,
//! well-separated peaks — the labelled seed the pre-pass estimates chemistry from.

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
}
