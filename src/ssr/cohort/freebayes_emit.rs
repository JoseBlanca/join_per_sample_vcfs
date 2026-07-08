//! Freebayes-style joint marginal-likelihood emission test for `ssr-call`
//! (spec [`doc/devel/specs/ssr_freebayes_marginal_emission.md`]).
//!
//! Instead of the heuristic emit gate (`is_variable` + `apply_fp_control`), decide
//! polymorphism at the **cohort level**: is the whole cohort's read evidence better
//! explained by a polymorphic site than by a fixed one? We port freebayes' structure
//! — a **neutral site-frequency-spectrum (Ewens `θ/k`) prior** over the allele-count
//! configuration, combined with the per-sample read likelihoods — onto OUR stutter
//! read model. The read likelihoods are reused verbatim as `data_ll[s][g]`
//! (`em::compute_data_ll`); the only new machinery is the SFS prior and the marginal.
//!
//! Because each sample's genotype prior is conditioned on the population allele
//! frequency `p` ([`em::genotype_prior`], the Wright form with inbreeding `F`), the
//! samples are **conditionally independent given `p`**. So the site marginal is an
//! exact sum over integer allele-count vectors `n` (with `Σ n_a = 2N` cohort
//! chromosomes) — no coupled combo enumeration, no banded search:
//!
//! ```text
//!   term(n)   = lnEwens(n)  +  Σ_s ln L_s(p = n / 2N)
//!   Z         = logsumexp_n term(n)
//!   P(mono)   = Σ_{n a corner} exp(term(n) − Z)          // fixed for one allele
//!   QUAL      = −10·log10 P(mono)                         // = −10·log10 P(monomorphic)
//! ```
//!
//! `L_s(p) = Σ_G genotype_prior(G, p, F_s)·exp(data_ll[s][G])`. The corners `n = 2N·e_a`
//! are the ways the locus is *fixed for one allele* (the generalised monomorphic
//! hypothesis — fixed-ref OR fixed-alt, matching the selfer where a locus can be fixed
//! for a non-reference allele).
//!
//! The count-vector enumeration is exact but grows as `C(2N + m − 1, m − 1)` in the
//! number of alleles `m`; we bound it to the top [`FREEBAYES_MAX_ALLELES`] alleles by
//! cohort frequency (`pi`). This is a computational bound, not a corroboration knob:
//! the excluded alleles are the rare (stutter) tail the SFS prior would suppress
//! anyway, and a genuinely multi-allelic locus still segregates through its strongest
//! alleles. For `m = 2` this is the exact biallelic marginal.

/// Population-scaled diversity `θ` for the Ewens SFS prior (freebayes' default `-T`).
/// Each distinct allele invoked in a configuration pays a factor `θ`, so a rare
/// "allele" (a stutter shoulder) must earn that cost from the read likelihood — the
/// property that suppresses systematic-stutter false positives. Fixed, not a per-run
/// knob (spec §4.4).
pub(crate) const SFS_THETA: f64 = 0.01;

/// Maximum number of alleles entering the exact count-vector marginal. The enumeration
/// is `C(2N + m − 1, m − 1)` in `m`; `4` keeps it tractable while covering ref + up to
/// three segregating alleles, more than a selfer's loci need. Alleles beyond this
/// (ranked by cohort frequency `pi`) are excluded from the hypothesis space.
pub(crate) const FREEBAYES_MAX_ALLELES: usize = 4;

/// Flat index of the diploid genotype `(a, b)` in `enumerate_diploid_genotypes(k)`
/// order (row-major over `i ≤ j`): row `i` starts at `i·k − i(i−1)/2`, offset `j − i`.
#[inline]
fn gt_index(a: usize, b: usize, k: usize) -> usize {
    let (i, j) = if a <= b { (a, b) } else { (b, a) };
    i * k - i * (i.saturating_sub(1)) / 2 + (j - i)
}

/// Online log-sum-exp accumulator (numerically stable, order-stable): tracks
/// `max` and `sum = Σ exp(x − max)` so the running total never overflows.
#[derive(Clone, Copy)]
struct LogSumExp {
    max: f64,
    sum: f64,
}

impl LogSumExp {
    fn new() -> Self {
        Self {
            max: f64::NEG_INFINITY,
            sum: 0.0,
        }
    }

    fn push(&mut self, x: f64) {
        if x == f64::NEG_INFINITY {
            return;
        }
        if x > self.max {
            self.sum = self.sum * (self.max - x).exp() + 1.0;
            self.max = x;
        } else {
            self.sum += (x - self.max).exp();
        }
    }

    /// `ln Σ exp(x_i)`; `−∞` if nothing was pushed.
    fn value(&self) -> f64 {
        if self.sum == 0.0 {
            f64::NEG_INFINITY
        } else {
            self.max + self.sum.ln()
        }
    }
}

/// The reduced Ewens ln-prior of a *labelled* allele-count vector, up to a constant
/// that is shared by every vector at the locus (so it cancels in the normalised
/// `P(monomorphic)` and is omitted). From freebayes' Ewens' Sampling Formula
/// ([`freebayes/src/Ewens.cpp`]) the per-configuration factor is `∏_{a: n_a>0} θ / n_a`;
/// in ln-space, dropping the shared `M! / (θ·∏(θ+h))` normaliser:
///
/// ```text
///   lnEwens(n) = (#nonzero alleles)·ln θ  −  Σ_{a: n_a>0} ln n_a
/// ```
///
/// The `#nonzero·ln θ` term is the rare-allele penalty: each additional allele
/// multiplies the prior by `θ ≈ 0.01`, so a monomorphic corner (`#nonzero = 1`) is
/// favoured over a biallelic configuration by a factor `θ` unless the read likelihood
/// pays for it. (We keep alleles *labelled* — distinct candidate sequences — and
/// enumerate each count vector once, so freebayes' `1/count_j!` exchangeability
/// correction, which is for unordered partitions, is intentionally absent.)
#[inline]
fn ln_ewens_reduced(counts: &[u32], ln_theta: f64) -> f64 {
    let mut acc = 0.0;
    for &c in counts {
        if c > 0 {
            acc += ln_theta - (c as f64).ln();
        }
    }
    acc
}

/// Per-sample scoring inputs for the marginal: the inbreeding `F_s` and the read-weight
/// matrix `w[a][b] = exp(data_ll[s][gt(a,b)] − max_s)` over the selected alleles
/// (`max_s` = per-sample max over selected genotypes, subtracted for stability; the
/// `Σ_s max_s` offset is shared by every count vector and cancels in `P(mono)`).
struct SampleWeights {
    f: f64,
    /// Row-major `m × m`, symmetric; `w[a*m + b]` is the read weight for genotype (a,b).
    w: Vec<f64>,
    m: usize,
}

impl SampleWeights {
    /// `ln L_s(p) = ln[ F·Σ_a p_a·w[a][a] + (1−F)·pᵀ w p ]` — the read likelihood
    /// marginalised over this sample's genotype under the Wright genotype prior
    /// (`em::genotype_prior` expanded: hom `F·p_a + (1−F)·p_a²`, het `(1−F)·2 p_a p_b`).
    #[inline]
    fn ln_likelihood(&self, p: &[f64]) -> f64 {
        let m = self.m;
        let mut hom = 0.0; // Σ_a p_a w[a][a]
        let mut quad = 0.0; // pᵀ w p
        for (a, &pa) in p.iter().enumerate() {
            let row_slice = &self.w[a * m..a * m + m];
            hom += pa * row_slice[a];
            let row: f64 = p.iter().zip(row_slice).map(|(&pb, &wab)| pb * wab).sum();
            quad += pa * row;
        }
        let l = self.f * hom + (1.0 - self.f) * quad;
        if l > 0.0 { l.ln() } else { f64::NEG_INFINITY }
    }
}

/// `ln P(monomorphic | data)` (≤ 0) from the freebayes-style cohort marginal — the raw
/// quantity the emission gate thresholds and the VCF QUAL is derived from. `data_ll[s][g]`
/// is the per-sample per-genotype read log-likelihood (`em::compute_data_ll`), `pi` the EM
/// cohort allele frequencies (to rank alleles), `f_per_present` the frozen per-sample
/// inbreeding `F`, `k` the candidate count, and `theta` the SFS `θ`. Returns `0.0`
/// (`P(mono) = 1` → uninformative/monomorphic) for `< 2` alleles or no present samples.
///
/// Deterministic and order-stable: samples and count vectors are summed in fixed order,
/// so the result is identical regardless of `--threads`.
pub(crate) fn ln_p_monomorphic(
    data_ll: &[Vec<f64>],
    pi: &[f64],
    f_per_present: &[f64],
    k: usize,
    theta: f64,
) -> f64 {
    let n_present = data_ll.len();
    let m = k.min(FREEBAYES_MAX_ALLELES);
    if m < 2 || n_present == 0 {
        return 0.0;
    }
    debug_assert_eq!(pi.len(), k);
    debug_assert_eq!(f_per_present.len(), n_present);

    // Select the top-`m` alleles by cohort frequency `pi` (ties broken by index for
    // determinism). The corners of the enumerated simplex are the fixed-for-`sel[a]`
    // hypotheses; a truly monomorphic locus concentrates `pi` on its fixed allele, so
    // that allele is selected and its corner captured.
    let mut order: Vec<usize> = (0..k).collect();
    order.sort_by(|&x, &y| {
        pi[y]
            .partial_cmp(&pi[x])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(x.cmp(&y))
    });
    let sel = &order[..m];

    // Per-sample read-weight matrices over the selected alleles.
    let samples: Vec<SampleWeights> = (0..n_present)
        .map(|s| {
            let ll = &data_ll[s];
            // max over the selected genotypes for per-sample stability.
            let mut max_s = f64::NEG_INFINITY;
            for a in 0..m {
                for b in a..m {
                    let v = ll[gt_index(sel[a], sel[b], k)];
                    if v > max_s {
                        max_s = v;
                    }
                }
            }
            let mut w = vec![0.0; m * m];
            for a in 0..m {
                for b in a..m {
                    let weight = (ll[gt_index(sel[a], sel[b], k)] - max_s).exp();
                    w[a * m + b] = weight;
                    w[b * m + a] = weight;
                }
            }
            SampleWeights {
                f: f_per_present[s],
                w,
                m,
            }
        })
        .collect();

    let n_chrom = 2 * n_present; // total cohort chromosomes present at the locus
    let ln_theta = theta.ln();
    let inv_n = 1.0 / n_chrom as f64;

    let mut full = LogSumExp::new();
    let mut corners = LogSumExp::new();
    let mut counts = vec![0u32; m];
    let mut p = vec![0.0f64; m];

    // Enumerate every allele-count vector with `Σ = n_chrom` and accumulate its term.
    enumerate_counts(
        &mut counts,
        &mut p,
        0,
        n_chrom as u32,
        n_chrom,
        inv_n,
        ln_theta,
        &samples,
        &mut full,
        &mut corners,
    );

    let ln_p_mono = corners.value() - full.value();
    // `P(mono) ≤ 1` ⇒ `ln P(mono) ≤ 0`; clamp float noise and guard the empty-marginal NaN.
    if ln_p_mono.is_finite() {
        ln_p_mono.min(0.0)
    } else {
        0.0
    }
}

/// Site `QUAL = −10·log10 P(monomorphic | data)` — a convenience wrapper over
/// [`ln_p_monomorphic`]. Not clamped (the caller applies `qual_cap`); `0.0` when
/// `P(mono) = 1`.
pub(crate) fn site_qual(
    data_ll: &[Vec<f64>],
    pi: &[f64],
    f_per_present: &[f64],
    k: usize,
    theta: f64,
) -> f64 {
    let ln_p_mono = ln_p_monomorphic(data_ll, pi, f_per_present, k, theta);
    // −10·log10(P) = −10·ln(P)/ln(10); `ln P(mono) ≤ 0` ⇒ QUAL ≥ 0.
    let qual = -10.0 * ln_p_mono / std::f64::consts::LN_10;
    if qual.is_finite() { qual.max(0.0) } else { 0.0 }
}

/// Recursively enumerate count vectors summing to `n_chrom` over `m` alleles (fixed
/// order), pushing each configuration's `term = lnEwens + Σ_s ln L_s(p)` to the full
/// accumulator and, if it is a corner (one allele holds all chromosomes), to the corner
/// accumulator.
#[allow(clippy::too_many_arguments)]
fn enumerate_counts(
    counts: &mut [u32],
    p: &mut [f64],
    pos: usize,
    remaining: u32,
    n_chrom: usize,
    inv_n: f64,
    ln_theta: f64,
    samples: &[SampleWeights],
    full: &mut LogSumExp,
    corners: &mut LogSumExp,
) {
    let m = counts.len();
    if pos == m - 1 {
        // Last allele takes whatever is left.
        counts[pos] = remaining;
        for a in 0..m {
            p[a] = counts[a] as f64 * inv_n;
        }
        let mut term = ln_ewens_reduced(counts, ln_theta);
        for s in samples {
            term += s.ln_likelihood(p);
            if term == f64::NEG_INFINITY {
                break;
            }
        }
        full.push(term);
        // A corner = fixed for one allele = one entry equal to n_chrom.
        if counts.iter().any(|&c| c as usize == n_chrom) {
            corners.push(term);
        }
        return;
    }
    for c in 0..=remaining {
        counts[pos] = c;
        enumerate_counts(
            counts,
            p,
            pos + 1,
            remaining - c,
            n_chrom,
            inv_n,
            ln_theta,
            samples,
            full,
            corners,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a `data_ll[s][g]` for `k` alleles from per-sample per-genotype log-likelihoods
    /// given as `(i, j, ll)` triples; unspecified genotypes get a very low likelihood.
    fn make_data_ll(k: usize, per_sample: &[Vec<(usize, usize, f64)>]) -> Vec<Vec<f64>> {
        let n_g = k * (k + 1) / 2;
        per_sample
            .iter()
            .map(|entries| {
                let mut row = vec![-50.0; n_g];
                for &(i, j, ll) in entries {
                    row[gt_index(i, j, k)] = ll;
                }
                row
            })
            .collect()
    }

    #[test]
    fn gt_index_matches_enumeration_order() {
        // Row-major i<=j order for k=3: (0,0)(0,1)(0,2)(1,1)(1,2)(2,2).
        let k = 3;
        let expected = [
            (0, 0, 0),
            (0, 1, 1),
            (0, 2, 2),
            (1, 1, 3),
            (1, 2, 4),
            (2, 2, 5),
        ];
        for &(i, j, idx) in &expected {
            assert_eq!(gt_index(i, j, k), idx);
            assert_eq!(gt_index(j, i, k), idx, "symmetric");
        }
    }

    #[test]
    fn segregating_locus_emits_high_qual() {
        // 8 plants: 4 read-dominant homozygous for the ref allele (0/0), 4 read-dominant
        // homozygous for the alt allele (1/1). A clean fixed difference between inbred
        // lines → strongly polymorphic → P(monomorphic) tiny → high QUAL.
        let k = 2;
        let ref_hom = vec![(0usize, 0usize, 0.0f64)];
        let alt_hom = vec![(1usize, 1usize, 0.0f64)];
        let per_sample: Vec<_> = (0..8)
            .map(|s| {
                if s < 4 {
                    ref_hom.clone()
                } else {
                    alt_hom.clone()
                }
            })
            .collect();
        let data_ll = make_data_ll(k, &per_sample);
        let pi = [0.5, 0.5];
        let f = vec![0.82; 8]; // selfer F
        let qual = site_qual(&data_ll, &pi, &f, k, SFS_THETA);
        assert!(
            qual > 60.0,
            "segregating locus should have high QUAL, got {qual}"
        );
    }

    #[test]
    fn systematic_stutter_locus_stays_low_qual() {
        // 8 plants: every plant a strong ref mode with a thin, consistent one-unit
        // shoulder (alt weakly supported everywhere, never dominant). No segregation →
        // the whole cohort is best explained as fixed-ref → P(monomorphic) ≈ 1 → QUAL ≈ 0.
        let k = 2;
        // ref hom strongly favoured; het/alt-hom much less likely but not impossible.
        let stutter = vec![(0usize, 0usize, 0.0f64), (0, 1, -2.5), (1, 1, -12.0)];
        let per_sample: Vec<_> = (0..8).map(|_| stutter.clone()).collect();
        let data_ll = make_data_ll(k, &per_sample);
        let pi = [0.9, 0.1];
        let f = vec![0.82; 8];
        let qual = site_qual(&data_ll, &pi, &f, k, SFS_THETA);
        assert!(
            qual < 10.0,
            "systematic-stutter locus should stay near QUAL 0, got {qual}"
        );
    }

    #[test]
    fn segregating_beats_stutter() {
        // The discriminative property: a real segregating locus outscores a
        // systematic-stutter locus of the same depth by a wide margin.
        let k = 2;
        let seg: Vec<_> = (0..8)
            .map(|s| {
                if s < 4 {
                    vec![(0usize, 0usize, 0.0f64)]
                } else {
                    vec![(1usize, 1usize, 0.0f64)]
                }
            })
            .collect();
        let stut: Vec<_> = (0..8)
            .map(|_| vec![(0usize, 0usize, 0.0f64), (0, 1, -2.5), (1, 1, -12.0)])
            .collect();
        let q_seg = site_qual(
            &make_data_ll(k, &seg),
            &[0.5, 0.5],
            &vec![0.82; 8],
            k,
            SFS_THETA,
        );
        let q_stut = site_qual(
            &make_data_ll(k, &stut),
            &[0.9, 0.1],
            &vec![0.82; 8],
            k,
            SFS_THETA,
        );
        assert!(q_seg > q_stut + 30.0, "seg {q_seg} vs stutter {q_stut}");
    }

    #[test]
    fn triallelic_segregation_emits() {
        // Three inbred alleles each fixed in a subset of plants (multiallelic path).
        let k = 3;
        let per_sample: Vec<_> = (0..9)
            .map(|s| vec![(s / 3, s / 3, 0.0f64)]) // plants 0-2 → 0/0, 3-5 → 1/1, 6-8 → 2/2
            .collect();
        let data_ll = make_data_ll(k, &per_sample);
        let pi = [0.34, 0.33, 0.33];
        let f = vec![0.82; 9];
        let qual = site_qual(&data_ll, &pi, &f, k, SFS_THETA);
        assert!(
            qual > 60.0,
            "triallelic segregation should emit, got {qual}"
        );
    }

    #[test]
    fn fewer_than_two_alleles_is_monomorphic() {
        let data_ll = make_data_ll(1, &[vec![(0, 0, 0.0)], vec![(0, 0, 0.0)]]);
        assert_eq!(
            site_qual(&data_ll, &[1.0], &[0.82, 0.82], 1, SFS_THETA),
            0.0
        );
    }

    #[test]
    fn deterministic_across_repeated_calls() {
        // The reduction is order-stable: identical inputs → bit-identical output.
        let k = 2;
        let per_sample: Vec<_> = (0..6)
            .map(|s| {
                if s % 2 == 0 {
                    vec![(0usize, 0usize, 0.0f64), (0, 1, -3.0)]
                } else {
                    vec![(1usize, 1usize, 0.0f64), (0, 1, -3.0)]
                }
            })
            .collect();
        let data_ll = make_data_ll(k, &per_sample);
        let a = site_qual(&data_ll, &[0.5, 0.5], &vec![0.5; 6], k, SFS_THETA);
        let b = site_qual(&data_ll, &[0.5, 0.5], &vec![0.5; 6], k, SFS_THETA);
        assert_eq!(a.to_bits(), b.to_bits());
    }
}
