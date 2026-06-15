//! Per-locus aggregation — fold a locus's per-read outcomes into its evidence
//! record (arch §8.2/§11, realign-everything + all-CSR revision).
//!
//! Each spanning read contributes a **pruned, renormalized `Qᵣ` profile** (its
//! plausible repeat lengths + per-read log-probabilities); non-spanning reads
//! contribute only to the QC tallies. There is no confident/ambiguous split and
//! no length histogram — every spanning read is stored as one profile, and zstd
//! compresses the redundancy on disk (decided 2026-06-15; see the §2 design
//! revision and the `.ssr.psp` size being a measure-first concern).
//!
//! This produces an **in-memory** [`SsrLocusRecord`]; flattening the profiles to
//! the container's CSR columns is the (deferred) container schema's job.
//! Off-ladder evidence is deferred (empty) until off-ladder candidate generation
//! lands ([`super::candidate_generation`]).

use crate::ssr::types::{Allele, Locus};

use super::read_analysis::ReadOutcome;

/// A candidate length is dropped from a read's stored `Qᵣ` profile when its
/// forward log-likelihood is more than this many **nats** below the read's best
/// length (then the survivors are renormalized, §11). A length that survives is
/// within `e^(−AMB_LL_DROP)` of the best. A **calibration** placeholder (arch §14):
/// adjacent rungs differ by ~one gap-open penalty (≈10 nats), so ~4 keeps clean
/// reads to a single length and genuine ambiguity to 2–3.
pub(crate) const AMB_LL_DROP: f64 = 4.0;

/// The QC counts the fetcher tallies over its full pre-triage pass — the reads
/// triage never sees (dropped by the admission gate) — which the aggregator
/// cannot derive from the per-read outcomes alone (arch §3.1/§8.3).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct QcCounts {
    /// Total reads considered at the locus.
    pub(crate) depth: u32,
    /// Reads dropped by the admission gate (low-MAPQ / duplicate / clipped).
    pub(crate) n_filtered: u32,
    /// Mapped reads overlapping the locus window — the normalized-depth denominator.
    pub(crate) mapped_reads: u32,
}

/// One sample's evidence at one locus (arch §8.2/§11). In-memory form: the
/// per-read profiles are kept as `Vec`s; the container flattens them to CSR
/// columns (deferred). `n_flank_indel` is intentionally absent — it was a
/// fast-path-gate-era count that realign-everything no longer produces.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SsrLocusRecord {
    pub(crate) chrom: Box<str>,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) depth: u32,
    pub(crate) n_spanning: u32,
    pub(crate) n_flanking: u32,
    pub(crate) n_frr: u32,
    pub(crate) n_filtered: u32,
    pub(crate) mapped_reads: u32,
    /// One renormalized `Qᵣ` profile per spanning read: `(on-ladder length,
    /// log-probability)` pairs, ascending by length, log-probs summing to 1 in
    /// linear space. Off-ladder evidence is deferred (not stored here yet).
    pub(crate) spanning: Vec<Vec<(u16, f32)>>,
}

/// Fold a locus's per-read outcomes + the fetcher's QC counts into its record.
/// Derives `n_spanning` / `n_flanking` / `n_frr` from the outcomes; copies the
/// rest from `counts` and the `locus`.
pub(crate) fn aggregate(
    locus: &Locus,
    outcomes: &[ReadOutcome],
    counts: QcCounts,
) -> SsrLocusRecord {
    let mut n_spanning = 0u32;
    let mut n_flanking = 0u32;
    let mut n_frr = 0u32;
    let mut spanning = Vec::new();

    for outcome in outcomes {
        match outcome {
            ReadOutcome::Spanning(scores) => {
                n_spanning += 1;
                spanning.push(prune_and_renormalize(scores));
            }
            ReadOutcome::Flanking => n_flanking += 1,
            ReadOutcome::InRepeat => n_frr += 1,
        }
    }

    SsrLocusRecord {
        chrom: locus.chrom().into(),
        start: locus.start(),
        end: locus.end(),
        depth: counts.depth,
        n_spanning,
        n_flanking,
        n_frr,
        n_filtered: counts.n_filtered,
        mapped_reads: counts.mapped_reads,
        spanning,
    }
}

/// Turn a read's dense `Qᵣ` (raw forward log-liks over its candidate window) into
/// the stored sparse profile: drop candidates more than `AMB_LL_DROP` below the
/// per-read max, then renormalize the survivors to log-probabilities (subtract
/// their log-sum-exp). Survivors keep candidate order, which `build_rungs`
/// emits ascending by length.
fn prune_and_renormalize(scores: &[(Allele, f64)]) -> Vec<(u16, f32)> {
    let max = scores
        .iter()
        .map(|(_, ll)| *ll)
        .fold(f64::NEG_INFINITY, f64::max);

    let survivors: Vec<(u16, f64)> = scores
        .iter()
        .filter(|(_, ll)| *ll >= max - AMB_LL_DROP)
        .filter_map(|(allele, ll)| match allele {
            Allele::OnLadder { units } => Some((*units, *ll)),
            // Off-ladder candidates are not generated yet (deferred); reaching
            // here means future wiring forgot to teach the aggregator about them.
            Allele::OffLadder(_) => {
                debug_assert!(false, "off-ladder aggregation not yet wired");
                None
            }
        })
        .collect();

    // Renormalize over the survivors: log_prob_i = ll_i − logsumexp(survivors).
    let z = {
        let sum_exp: f64 = survivors.iter().map(|(_, ll)| (*ll - max).exp()).sum();
        max + sum_exp.ln()
    };
    survivors
        .into_iter()
        .map(|(units, ll)| (units, (ll - z) as f32))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::types::Motif;

    fn locus() -> Locus {
        Locus::new(
            "chr1".into(),
            16,
            22,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGGGGCACACATTTTTT").into(),
            10,
        )
        .unwrap()
    }

    fn on(units: u16, ll: f64) -> (Allele, f64) {
        (Allele::OnLadder { units }, ll)
    }

    /// log-probs in a profile exponentiate to sum 1 (a proper distribution).
    fn assert_normalized(profile: &[(u16, f32)]) {
        let sum: f64 = profile.iter().map(|(_, lp)| (*lp as f64).exp()).sum();
        assert!((sum - 1.0).abs() < 1e-5, "profile mass {sum} != 1");
    }

    #[test]
    fn derives_the_three_classification_counts() {
        let outcomes = vec![
            ReadOutcome::Spanning(vec![on(12, -0.01)]),
            ReadOutcome::Spanning(vec![on(12, -0.01)]),
            ReadOutcome::Flanking,
            ReadOutcome::InRepeat,
            ReadOutcome::InRepeat,
        ];
        let rec = aggregate(&locus(), &outcomes, QcCounts::default());
        assert_eq!(rec.n_spanning, 2);
        assert_eq!(rec.n_flanking, 1);
        assert_eq!(rec.n_frr, 2);
        assert_eq!(rec.spanning.len(), 2);
    }

    #[test]
    fn copies_locus_coords_and_fetcher_counts() {
        let counts = QcCounts {
            depth: 30,
            n_filtered: 4,
            mapped_reads: 33,
        };
        let rec = aggregate(&locus(), &[], counts);
        assert_eq!((&*rec.chrom, rec.start, rec.end), ("chr1", 16, 22));
        assert_eq!((rec.depth, rec.n_filtered, rec.mapped_reads), (30, 4, 33));
    }

    #[test]
    fn a_sharp_read_prunes_to_a_single_length_at_log_prob_zero() {
        // A clean read: 12 dominates, neighbours ~10 nats down (one gap each).
        let scores = vec![on(11, -10.4), on(12, -0.01), on(13, -10.4)];
        let profile = prune_and_renormalize(&scores);
        assert_eq!(profile.len(), 1);
        assert_eq!(profile[0].0, 12);
        assert!((profile[0].1).abs() < 1e-6); // ln(1) = 0
    }

    #[test]
    fn a_bimodal_read_keeps_both_lengths_normalized() {
        // 11 and 12 within AMB_LL_DROP of each other → both survive, 13 dropped.
        let scores = vec![on(11, -2.8), on(12, -2.0), on(13, -12.0)];
        let profile = prune_and_renormalize(&scores);
        let lengths: Vec<u16> = profile.iter().map(|(u, _)| *u).collect();
        assert_eq!(lengths, vec![11, 12]); // ascending, 13 pruned
        assert_normalized(&profile);
        // 12 is the more likely of the two.
        assert!(profile[1].1 > profile[0].1);
    }

    #[test]
    fn pruning_drops_the_far_tail() {
        // Only the peak is within 4 nats; everything else is pruned.
        let scores = vec![
            on(9, -31.0),
            on(10, -20.8),
            on(11, -10.4),
            on(12, -0.01),
            on(13, -10.4),
        ];
        let profile = prune_and_renormalize(&scores);
        assert_eq!(profile, vec![(12, 0.0)]);
    }

    #[test]
    fn profiles_are_normalized_per_read() {
        let outcomes = vec![ReadOutcome::Spanning(vec![on(11, -2.8), on(12, -2.0)])];
        let rec = aggregate(&locus(), &outcomes, QcCounts::default());
        assert_normalized(&rec.spanning[0]);
    }
}
