//! Candidate-allele types — the A1 nouns of the genotyping candidate set (arch
//! `ssr_call_genotyping.md` §2, implementation plan §3.3).
//!
//! The assembly logic (rungs → nomination → union → ref-seed → locus-admission
//! motif filter) is written in C1; this file defines only the shapes it produces.

/// Why a locus was (or was not) admitted to genotyping — the per-site FILTER
/// reason (arch §6).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Admission {
    /// The locus is admitted.
    Pass,
    /// Adjacent-rung spacing did not match the motif length — not a periodic
    /// ladder (`notPeriodic`).
    NotPeriodic,
    /// The candidate count exceeded the per-locus cap (`tooManyAlleles`).
    TooManyAlleles,
    /// Too little depth to resolve candidates (`lowDepth`).
    LowDepth,
}

/// The candidate alleles for one locus, plus its admission verdict.
///
/// Each entry of [`Self::alleles`] is an independent candidate **tract sequence**
/// (impure peaks kept first-class); [`Self::ref_idx`] indexes the reference allele,
/// which is seeded unconditionally.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CandidateSet {
    /// Candidate tract sequences (each an independent allele).
    pub(crate) alleles: Vec<Box<[u8]>>,
    /// Index of the reference allele within [`Self::alleles`].
    pub(crate) ref_idx: usize,
    /// The site's admission verdict → FILTER.
    pub(crate) admit: Admission,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn candidate_set_ref_idx_points_into_alleles() {
        let cs = CandidateSet {
            alleles: vec![
                Box::from(b"ATAT".as_slice()),
                Box::from(b"ATATAT".as_slice()),
            ],
            ref_idx: 0,
            admit: Admission::Pass,
        };
        assert_eq!(&*cs.alleles[cs.ref_idx], b"ATAT");
        assert_eq!(cs.admit, Admission::Pass);
    }

    #[test]
    fn admission_variants_are_distinct() {
        assert_ne!(Admission::Pass, Admission::NotPeriodic);
        assert_ne!(Admission::TooManyAlleles, Admission::LowDepth);
    }
}
