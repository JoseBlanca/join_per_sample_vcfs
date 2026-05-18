//! Per-`(ploidy, n_alleles)` precomputed genotype shape, cached for
//! the EM loop's hot path.
//!
//! Every record processed by [`super::PosteriorEngine`] enumerates
//! genotypes, builds an allele-count table, and computes the
//! multinomial log-coefficients. All three quantities are pure
//! functions of `(ploidy, n_alleles)`, and that pair is heavily
//! concentrated in real workloads (overwhelmingly diploid biallelic
//! SNPs). A small thread-local slot-array cache amortises the build
//! cost across records with the same shape.
//!
//! The cached shape also carries two precomputed views that let
//! `e_step` skip per-cell work:
//!
//! - `nonzero_pairs` — for each genotype, the list of
//!   `(allele_idx, count)` pairs with `count > 0`. Replaces the inner
//!   sum's `if k == 0` skip branch.
//! - `homozygous_allele_for` — for each genotype, `Some(a)` if it's
//!   homozygous for allele `a`, else `None`. Replaces the per-cell
//!   linear `homozygous_allele(gt_counts)` scan.
//!
//! Bit-identical to the unoptimised path. Skipping `k == 0` entries
//! doesn't change float arithmetic (`0.0 * log_p + y == y` exactly
//! under IEEE 754), and the precomputed homozygous-allele lookup
//! returns the same value as the scan it replaces.

use std::cell::RefCell;
use std::sync::Arc;

use crate::var_calling::per_group_merger::genotype_order;

/// Cap the slot-array cache at `ploidy ≤ 8` and `n_alleles ≤ 16`.
/// Out-of-range shapes are built fresh on each call (no caching) so the
/// slot array stays small; real workloads almost never exceed these
/// bounds (Stage 5's `max_alleles` defaults to 6).
const MAX_CACHED_PLOIDY: usize = 8;
const MAX_CACHED_N_ALLELES: usize = 16;
const CACHE_SLOTS: usize = (MAX_CACHED_PLOIDY + 1) * (MAX_CACHED_N_ALLELES + 1);

/// Per-`(ploidy, n_alleles)` artefacts the EM loop needs but doesn't
/// change between records of the same shape.
pub(crate) struct GenotypeShape {
    pub(crate) n_genotypes: usize,
    /// Flat row-major `n_genotypes × n_alleles` allele-count table.
    /// `[g_idx * n_alleles + a]` is the number of copies of allele `a`
    /// in genotype `g_idx`.
    pub(crate) genotype_allele_counts: Vec<u32>,
    /// `ln C(ploidy; counts)` per genotype.
    pub(crate) log_multinomial_coeffs: Vec<f64>,
    /// Concatenated lists of `(allele_idx, count)` pairs with
    /// `count > 0`, indexed via `nonzero_pairs_offsets`. For diploid:
    /// each genotype contributes 1 (homozygous) or 2 (heterozygous)
    /// entries.
    pub(crate) nonzero_pairs: Vec<(u8, u32)>,
    /// `(start, len)` into `nonzero_pairs`, per genotype.
    pub(crate) nonzero_pairs_offsets: Vec<(u32, u8)>,
    /// Per genotype: `Some(a)` if homozygous for allele `a`, else
    /// `None`.
    pub(crate) homozygous_allele_for: Vec<Option<u8>>,
}

thread_local! {
    static SHAPE_CACHE: RefCell<[Option<Arc<GenotypeShape>>; CACHE_SLOTS]>
        = const { RefCell::new([const { None }; CACHE_SLOTS]) };
}

/// Return the cached [`GenotypeShape`] for `(ploidy, n_alleles)`,
/// building and caching it on first use. Shapes outside the cache
/// bounds are built fresh on every call — correct but uncached.
pub(crate) fn shape_for(ploidy: u8, n_alleles: usize) -> Arc<GenotypeShape> {
    let p = ploidy as usize;
    if p > MAX_CACHED_PLOIDY || n_alleles > MAX_CACHED_N_ALLELES {
        return Arc::new(build(ploidy, n_alleles));
    }
    let idx = p * (MAX_CACHED_N_ALLELES + 1) + n_alleles;
    SHAPE_CACHE.with(|cell| {
        let mut slots = cell.borrow_mut();
        if let Some(shape) = &slots[idx] {
            return Arc::clone(shape);
        }
        let shape = Arc::new(build(ploidy, n_alleles));
        slots[idx] = Some(Arc::clone(&shape));
        shape
    })
}

fn build(ploidy: u8, n_alleles: usize) -> GenotypeShape {
    let genotypes = genotype_order(ploidy, n_alleles);
    let n_genotypes = genotypes.len();

    let mut genotype_allele_counts = vec![0_u32; n_genotypes * n_alleles];
    for (g_idx, g) in genotypes.iter().enumerate() {
        let row = &mut genotype_allele_counts[g_idx * n_alleles..(g_idx + 1) * n_alleles];
        for &a in g {
            row[a as usize] += 1;
        }
    }

    let log_multinomial_coeffs: Vec<f64> = genotype_allele_counts
        .chunks_exact(n_alleles)
        .map(|counts| log_multinomial_coefficient(ploidy, counts))
        .collect();

    let mut nonzero_pairs: Vec<(u8, u32)> = Vec::with_capacity(n_genotypes * ploidy as usize);
    let mut nonzero_pairs_offsets: Vec<(u32, u8)> = Vec::with_capacity(n_genotypes);
    for counts in genotype_allele_counts.chunks_exact(n_alleles) {
        let start = nonzero_pairs.len() as u32;
        let mut len = 0_u8;
        for (a, &k) in counts.iter().enumerate() {
            if k != 0 {
                nonzero_pairs.push((a as u8, k));
                len += 1;
            }
        }
        nonzero_pairs_offsets.push((start, len));
    }

    let homozygous_allele_for: Vec<Option<u8>> = genotype_allele_counts
        .chunks_exact(n_alleles)
        .map(|counts| homozygous_allele(counts).map(|a| a as u8))
        .collect();

    GenotypeShape {
        n_genotypes,
        genotype_allele_counts,
        log_multinomial_coeffs,
        nonzero_pairs,
        nonzero_pairs_offsets,
        homozygous_allele_for,
    }
}

/// `lgamma(n + 1)` for small non-negative `n`. The cache only invokes
/// this for `n ≤ ploidy` (≤ 8 in practice), so a tiny iterative
/// product is fastest and exact for all relevant inputs.
fn log_factorial(n: u32) -> f64 {
    let mut acc = 0.0_f64;
    for i in 2..=n {
        acc += (i as f64).ln();
    }
    acc
}

fn log_multinomial_coefficient(ploidy: u8, counts: &[u32]) -> f64 {
    let mut log_coef = log_factorial(ploidy as u32);
    for &k in counts {
        log_coef -= log_factorial(k);
    }
    log_coef
}

/// If every non-zero entry of `counts` is at the same allele index,
/// returns that index. Used to flag genotypes the HWE-with-`F` IBD
/// partition contributes mass to.
fn homozygous_allele(counts: &[u32]) -> Option<usize> {
    let mut found: Option<usize> = None;
    for (a, &k) in counts.iter().enumerate() {
        if k == 0 {
            continue;
        }
        match found {
            None => found = Some(a),
            Some(_) => return None,
        }
    }
    found
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // ---- Cache behaviour ----

    #[test]
    fn shape_for_returns_arc_equal_on_repeated_lookups_within_bounds() {
        let a = shape_for(2, 2);
        let b = shape_for(2, 2);
        assert!(Arc::ptr_eq(&a, &b));
    }

    #[test]
    fn shape_for_returns_distinct_shapes_for_different_ploidy_or_n_alleles() {
        let s22 = shape_for(2, 2);
        let s23 = shape_for(2, 3);
        let s32 = shape_for(3, 2);
        assert!(!Arc::ptr_eq(&s22, &s23));
        assert!(!Arc::ptr_eq(&s22, &s32));
    }

    #[test]
    fn shape_for_returns_fresh_shape_for_out_of_bounds_ploidy() {
        let a = shape_for(9, 2);
        let b = shape_for(9, 2);
        // Out-of-bounds path doesn't cache, so each call allocates.
        assert!(!Arc::ptr_eq(&a, &b));
        // Content still matches.
        assert_eq!(a.n_genotypes, b.n_genotypes);
    }

    #[test]
    fn shape_for_returns_fresh_shape_for_out_of_bounds_n_alleles() {
        let a = shape_for(2, 17);
        let b = shape_for(2, 17);
        assert!(!Arc::ptr_eq(&a, &b));
        assert_eq!(a.n_genotypes, b.n_genotypes);
    }

    // ---- Diploid biallelic content ----

    #[test]
    fn diploid_biallelic_shape_has_three_genotypes_with_expected_counts() {
        let s = shape_for(2, 2);
        assert_eq!(s.n_genotypes, 3);
        // Canonical order is AA, AB, BB.
        assert_eq!(&s.genotype_allele_counts[0..2], &[2, 0]);
        assert_eq!(&s.genotype_allele_counts[2..4], &[1, 1]);
        assert_eq!(&s.genotype_allele_counts[4..6], &[0, 2]);
    }

    #[test]
    fn diploid_biallelic_multinomial_coeffs_are_zero_ln2_zero() {
        let s = shape_for(2, 2);
        assert!(approx(s.log_multinomial_coeffs[0], 0.0, 1e-12));
        assert!(approx(s.log_multinomial_coeffs[1], 2.0_f64.ln(), 1e-12));
        assert!(approx(s.log_multinomial_coeffs[2], 0.0, 1e-12));
    }

    #[test]
    fn diploid_biallelic_nonzero_pairs_have_lengths_one_two_one() {
        let s = shape_for(2, 2);
        let offsets: Vec<_> = s
            .nonzero_pairs_offsets
            .iter()
            .map(|&(_, len)| len)
            .collect();
        assert_eq!(offsets, vec![1, 2, 1]);
    }

    #[test]
    fn diploid_biallelic_homozygous_allele_for_has_some_none_some() {
        let s = shape_for(2, 2);
        assert_eq!(s.homozygous_allele_for, vec![Some(0), None, Some(1)]);
    }

    // ---- Helpers (formerly tested in posterior_engine.rs) ----

    #[test]
    fn log_factorial_returns_zero_for_zero_and_one() {
        assert_eq!(log_factorial(0), 0.0);
        assert_eq!(log_factorial(1), 0.0);
        assert!(approx(log_factorial(5), (120.0_f64).ln(), 1e-12));
    }

    #[test]
    fn log_multinomial_coefficient_matches_closed_form() {
        // C(2;1,1) = 2; ln 2.
        assert!(approx(
            log_multinomial_coefficient(2, &[1, 1]),
            2.0_f64.ln(),
            1e-12
        ));
        // C(2;2,0) = 1; ln 1 = 0.
        assert!(approx(log_multinomial_coefficient(2, &[2, 0]), 0.0, 1e-12));
        // C(4;2,2) = 6; ln 6.
        assert!(approx(
            log_multinomial_coefficient(4, &[2, 2]),
            6.0_f64.ln(),
            1e-12
        ));
    }

    #[test]
    fn homozygous_allele_returns_index_for_pure_genotypes() {
        assert_eq!(homozygous_allele(&[2, 0, 0]), Some(0));
        assert_eq!(homozygous_allele(&[0, 2, 0]), Some(1));
        assert_eq!(homozygous_allele(&[0, 0, 2]), Some(2));
        assert_eq!(homozygous_allele(&[1, 1, 0]), None);
        assert_eq!(homozygous_allele(&[0, 0, 0]), None);
    }
}
