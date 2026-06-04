//! `#[cfg(test)]`-only fixture builders for the `var_calling` test suites
//! (copied from `var_calling::test_helpers` — only the two builders the new
//! package's tests use, so the package is self-contained for the P7 swap when
//! the old package is deleted).

use crate::pileup_record::{AlleleObservation, AlleleSupportStats, ChainId, PileupRecord};

/// Build an [`AlleleObservation`] with a fully-specified
/// `(num_obs, q_sum, chain_ids)` triple. The strand-bias and MAPQ moments are
/// filled in from `num_obs` with stable default shapes — no test on this slice
/// inspects them.
pub(crate) fn allele(
    seq: &[u8],
    num_obs: u32,
    q_sum: f64,
    chain_ids: &[ChainId],
) -> AlleleObservation {
    AlleleObservation::new(
        seq.to_vec(),
        AlleleSupportStats::new(
            num_obs,
            q_sum,
            num_obs,
            0,
            0,
            num_obs * 30,
            u64::from(num_obs) * 900,
        ),
        chain_ids.to_vec(),
    )
}

/// Build a [`PileupRecord`] at `pos` on chromosome 0.
pub(crate) fn record(pos: u32, alleles: Vec<AlleleObservation>) -> PileupRecord {
    PileupRecord::new(0, pos, alleles)
}
