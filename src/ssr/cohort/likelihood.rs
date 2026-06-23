//! The read likelihood `Qᵣ(obs | candidate)` and the genotype likelihood (arch
//! `ssr_call_genotyping.md` §4, spec §6, implementation plan §3.5, C2).
//!
//! `Qᵣ` is HipSTR's sum-over-slips, combining the B2 kernel and the B3 align:
//!
//! ```text
//! Qᵣ(obs | cand) = Σ_Δ S_θ(Δ) · align(obs | cand ⊕ Δ)
//! align(obs | cand ⊕ Δ) = Σ_v (1/|v|) · align_subst(obs | v)   over placement variants v
//! ```
//!
//! The genotype likelihood mixes `Qᵣ` over the genotype's alleles and adds the
//! uniform **outlier** term `λ·(1/D)` (`D` = distinct sequences at the locus) so a
//! read matching no allele goes to junk rather than inventing a false allele. The
//! **allele-balance** defence is deferred to E2 (it acts on deconvolved E-step
//! responsibilities, not here). `align` is recomputed on demand — no cache in v1
//! (Q-G3).

use crate::ssr::cohort::pair_hmm::{HmmScratch, align_subst};
use crate::ssr::cohort::param_estimation::{MAX_SLIP, StutterShape};
use crate::ssr::cohort::stutter::{PlacementVariant, reach_variants, s_theta};
use crate::ssr::types::Motif;

/// Reused buffers for one read-vs-candidate evaluation (the placement-variant set
/// and the align DP), so scoring a locus's reads allocates once.
#[derive(Debug, Default)]
pub(crate) struct LikelihoodScratch {
    variants: Vec<PlacementVariant>,
    hmm: HmmScratch,
}

impl LikelihoodScratch {
    pub(crate) fn new() -> Self {
        Self::default()
    }
}

/// `Qᵣ(obs | cand)` — the probability the observed sequence arose from `cand` via a
/// (marginalized) stutter slip plus within-tract substitutions.
pub(crate) fn read_likelihood(
    obs: &[u8],
    cand: &[u8],
    motif: &Motif,
    shape: &StutterShape,
    level: f64,
    eps: f64,
    scratch: &mut LikelihoodScratch,
) -> f64 {
    let LikelihoodScratch { variants, hmm } = scratch;
    let mut q = 0.0;
    for delta in -(MAX_SLIP as i32)..=(MAX_SLIP as i32) {
        let s = s_theta(delta, shape, level);
        if s == 0.0 {
            continue;
        }
        reach_variants(cand, motif, delta, variants);
        if variants.is_empty() {
            continue; // this Δ cannot be reached from `cand` (e.g. over-contraction)
        }
        let pr_v = 1.0 / variants.len() as f64;
        let mut align = 0.0;
        for variant in variants.iter() {
            align += pr_v * align_subst(obs, &variant.seq, eps, hmm);
        }
        q += s * align;
    }
    q
}

/// `P(read | genotype)` — the per-allele `Qᵣ` mixed uniformly over the genotype's
/// allele copies, with the uniform outlier term `λ·(1/D)`.
///
/// `allele_qr` holds the precomputed `Qᵣ(obs | allele)` for each *candidate* allele;
/// `genotype` is the multiset of candidate indices making up the genotype (a
/// homozygote repeats an index, a diploid het has two). `distinct_count` is `D`.
pub(crate) fn read_given_genotype(
    allele_qr: &[f64],
    genotype: &[usize],
    lambda: f64,
    distinct_count: usize,
) -> f64 {
    debug_assert!(
        !genotype.is_empty(),
        "a genotype has at least one allele copy"
    );
    let mix: f64 = genotype.iter().map(|&a| allele_qr[a]).sum::<f64>() / genotype.len() as f64;
    let junk = 1.0 / distinct_count.max(1) as f64;
    (1.0 - lambda) * mix + lambda * junk
}

#[cfg(test)]
mod tests {
    use super::*;

    fn motif(bytes: &[u8]) -> Motif {
        Motif::new(bytes).unwrap()
    }

    fn shape() -> StutterShape {
        StutterShape {
            up_rate: 1.0,
            down_rate: 2.0,
            decay: 0.2,
        }
    }

    fn ca(units: usize) -> Vec<u8> {
        std::iter::repeat_n(*b"CA", units).flatten().collect()
    }

    const EPS: f64 = 0.01;
    const LEVEL: f64 = 0.1;

    #[test]
    fn faithful_read_is_dominated_by_the_zero_slip_term() {
        let mut scratch = LikelihoodScratch::new();
        let cand = ca(8);
        let q = read_likelihood(
            &cand,
            &cand,
            &motif(b"CA"),
            &shape(),
            LEVEL,
            EPS,
            &mut scratch,
        );
        // The Δ=0 term alone is (1−level)·(1−ε)^len; the total is a touch more (other
        // slips align weakly), and stays below 1.
        let faithful = (1.0 - LEVEL) * (1.0 - EPS).powi(16);
        assert!(
            q >= faithful,
            "q {q} should be at least the faithful term {faithful}"
        );
        assert!(
            q < faithful * 1.05,
            "q {q} should be dominated by the faithful term"
        );
    }

    #[test]
    fn a_read_prefers_its_own_allele_over_a_distant_one() {
        let mut scratch = LikelihoodScratch::new();
        let obs = ca(8);
        let own = read_likelihood(
            &obs,
            &ca(8),
            &motif(b"CA"),
            &shape(),
            LEVEL,
            EPS,
            &mut scratch,
        );
        let distant = read_likelihood(
            &obs,
            &ca(3),
            &motif(b"CA"),
            &shape(),
            LEVEL,
            EPS,
            &mut scratch,
        );
        assert!(own > distant, "own {own} should beat distant {distant}");
    }

    #[test]
    fn a_minus_one_stutter_read_is_explained_by_its_parent() {
        // obs is one unit short of the candidate — explained by the Δ=−1 slip term.
        let mut scratch = LikelihoodScratch::new();
        let cand = ca(8);
        let obs = ca(7);
        let q = read_likelihood(
            &obs,
            &cand,
            &motif(b"CA"),
            &shape(),
            LEVEL,
            EPS,
            &mut scratch,
        );
        // Approximately S_θ(−1)·(1−ε)^14 (the dominant path).
        let down_share = 2.0 / 3.0; // down_rate / (up+down)
        let expected = LEVEL * down_share * (1.0 - 0.2) * (1.0 - EPS).powi(14);
        assert!(q > 0.0);
        assert!(
            (q / expected - 1.0).abs() < 0.1,
            "q {q} should track the −1 stutter path {expected}"
        );
    }

    #[test]
    fn read_likelihood_is_a_probability_at_most_one() {
        let mut scratch = LikelihoodScratch::new();
        for units in [3usize, 6, 8, 12] {
            let cand = ca(units);
            let q = read_likelihood(
                &cand,
                &cand,
                &motif(b"CA"),
                &shape(),
                LEVEL,
                EPS,
                &mut scratch,
            );
            assert!((0.0..=1.0 + 1e-9).contains(&q), "q = {q} for units {units}");
        }
    }

    #[test]
    fn impure_candidate_sums_over_placements() {
        // An interrupted candidate: read_likelihood runs (Σ_v has >1 term) and the
        // faithful read scores positively.
        let mut scratch = LikelihoodScratch::new();
        let cand = b"CACACATTCACA"; // (CA)3 TT (CA)2
        let q = read_likelihood(
            cand,
            cand,
            &motif(b"CA"),
            &shape(),
            LEVEL,
            EPS,
            &mut scratch,
        );
        assert!(q > 0.0);
        assert!(q <= 1.0 + 1e-9);
    }

    #[test]
    fn genotype_likelihood_mixes_alleles_and_adds_outlier_floor() {
        // Two candidates; a het mixes them, and the outlier floor keeps a
        // no-match read non-zero.
        let allele_qr = [0.8, 0.2];
        let het = read_given_genotype(&allele_qr, &[0, 1], 0.01, 5);
        let hom0 = read_given_genotype(&allele_qr, &[0, 0], 0.01, 5);
        assert!(
            hom0 > het,
            "the read matches allele 0 better, so hom0 > het"
        );

        // A read matching neither candidate (both Qr ≈ 0) still gets λ/D.
        let outlier = read_given_genotype(&[0.0, 0.0], &[0, 1], 0.05, 4);
        assert!((outlier - 0.05 * 0.25).abs() < 1e-12);
    }
}
