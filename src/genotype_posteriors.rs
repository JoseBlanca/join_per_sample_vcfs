//! Estimates genotype posterior probabilities using an EM algorithm.
//!
//! This implements the same approach as GATK's GenotypeGVCFs: it pools
//! evidence across all samples to estimate allele frequencies, then uses
//! those frequencies to refine each sample's genotype call.

/// Configuration for the Dirichlet prior on allele frequencies.
/// These values encode our expectation of how common variants are
/// before looking at any data.
pub struct PriorConfig {
    /// How common are SNPs in the population (default: 0.001)
    pub snp_heterozygosity: f64,
    /// How common are indels in the population (default: 1.25e-4)
    pub indel_heterozygosity: f64,
    /// Uncertainty in our heterozygosity estimate (default: 0.01)
    pub heterozygosity_std_dev: f64,
}

impl Default for PriorConfig {
    fn default() -> Self {
        PriorConfig {
            snp_heterozygosity: 0.001,
            indel_heterozygosity: 1.25e-4,
            heterozygosity_std_dev: 0.01,
        }
    }
}

/// One possible diploid genotype at a site.
/// For example, with alleles [A, T, C]:
///   - genotype 0/0 (AA) has allele_indices = [0, 0]
///   - genotype 0/1 (AT) has allele_indices = [0, 1]
///   - genotype 1/2 (TC) has allele_indices = [1, 2]
pub struct DiploidGenotype {
    /// The two allele indices (e.g., [0, 1] for a het 0/1).
    pub allele_indices: [usize; 2],
    /// Position of this genotype in the PL array (VCF ordering).
    pub pl_index: usize,
}

/// Result of the EM algorithm for one site across all samples.
pub struct SitePosteriors {
    /// Estimated allele frequencies after convergence (one per allele).
    pub allele_frequencies: Vec<f64>,
    /// Per-sample genotype posteriors, stored flat.
    /// For sample `i` with `g` possible genotypes:
    ///   posteriors[i * g .. (i+1) * g]
    pub genotype_posteriors: Vec<f64>,
    /// Number of possible genotypes (same for all samples at this site).
    pub num_genotypes: usize,
    /// QUAL score: confidence that the site is variable.
    pub qual: f64,
}

/// Runs the full EM algorithm for one variant site.
///
/// # Arguments
/// * `num_alleles` - Total number of alleles (1 ref + N alt)
/// * `sample_pls` - PL values for each sample, flat: sample_pls[i * num_genotypes + g]
/// * `prior` - Prior configuration for allele frequencies
///
/// # Returns
/// The estimated allele frequencies, per-sample genotype posteriors, and QUAL.
pub fn estimate_posteriors(
    num_alleles: usize,
    num_samples: usize,
    sample_pls: &[f64],
    prior: &PriorConfig,
) -> SitePosteriors {
    let genotypes = enumerate_diploid_genotypes(num_alleles);
    let num_genotypes = genotypes.len();

    let prior_pseudocounts = compute_prior_pseudocounts(num_alleles, prior);

    // EM iteration
    let mut allele_frequencies = vec![1.0 / num_alleles as f64; num_alleles];
    let mut genotype_posteriors = vec![0.0; num_samples * num_genotypes];

    for _iteration in 0..50 {
        // E-step: compute each sample's genotype posteriors given current frequencies
        compute_genotype_posteriors(
            &allele_frequencies,
            &genotypes,
            sample_pls,
            num_samples,
            &mut genotype_posteriors,
        );

        // M-step: update allele frequencies from the posteriors
        let new_frequencies = update_allele_frequencies(
            &genotype_posteriors,
            &genotypes,
            num_samples,
            num_alleles,
            &prior_pseudocounts,
        );

        // Check convergence: stop if frequencies barely changed
        let max_change = allele_frequencies
            .iter()
            .zip(new_frequencies.iter())
            .map(|(old, new)| (old - new).abs())
            .fold(0.0_f64, f64::max);

        allele_frequencies = new_frequencies;

        if max_change < 1e-6 {
            break;
        }
    }

    // Final E-step with converged frequencies
    compute_genotype_posteriors(
        &allele_frequencies,
        &genotypes,
        sample_pls,
        num_samples,
        &mut genotype_posteriors,
    );

    let qual = compute_qual(&genotype_posteriors, num_samples, num_genotypes);

    SitePosteriors {
        allele_frequencies,
        genotype_posteriors,
        num_genotypes,
        qual,
    }
}

// === Helper functions (filled in below) ===

/// Lists all possible diploid genotypes for the given number of alleles.
///
/// For 2 alleles (ref=0, alt=1) this returns: 0/0, 0/1, 1/1
/// For 3 alleles (ref=0, alt1=1, alt2=2): 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
///
/// The ordering matches the VCF PL field convention.
fn enumerate_diploid_genotypes(num_alleles: usize) -> Vec<DiploidGenotype> {
    let mut genotypes = Vec::new();
    let mut pl_index = 0;

    for second_allele in 0..num_alleles {
        for first_allele in 0..=second_allele {
            genotypes.push(DiploidGenotype {
                allele_indices: [first_allele, second_allele],
                pl_index,
            });
            pl_index += 1;
        }
    }

    genotypes
}

/// Computes the Dirichlet prior pseudocounts for each allele.
///
/// The reference allele gets a large pseudocount (reflecting that most sites
/// are not variable), while each alt allele gets a small pseudocount
/// (reflecting the expected heterozygosity).
///
/// With defaults: ref gets 10.0, each SNP alt gets 0.01.
fn compute_prior_pseudocounts(num_alleles: usize, prior: &PriorConfig) -> Vec<f64> {
    let ref_pseudocount =
        prior.snp_heterozygosity / (prior.heterozygosity_std_dev * prior.heterozygosity_std_dev);
    let alt_pseudocount = prior.snp_heterozygosity * ref_pseudocount;

    let mut pseudocounts = Vec::with_capacity(num_alleles);
    pseudocounts.push(ref_pseudocount); // allele 0 = reference
    for _ in 1..num_alleles {
        pseudocounts.push(alt_pseudocount); // alt alleles
    }

    pseudocounts
}

/// E-step: for each sample, compute the posterior probability of each genotype.
///
/// For each genotype, the posterior combines three things:
/// 1. The sequencing evidence (from the PL field)
/// 2. The Hardy-Weinberg expectation given current allele frequencies
/// 3. A combinatorial factor for heterozygotes (het has 2 arrangements, hom has 1)
///
/// The results are normalized so that each sample's posteriors sum to 1.
fn compute_genotype_posteriors(
    allele_frequencies: &[f64],
    genotypes: &[DiploidGenotype],
    sample_pls: &[f64],
    num_samples: usize,
    out_posteriors: &mut [f64],
) {
    let num_genotypes = genotypes.len();

    for sample in 0..num_samples {
        let pl_offset = sample * num_genotypes;
        let post_offset = sample * num_genotypes;

        // Compute unnormalized log10-posterior for each genotype
        let mut max_log10_posterior = f64::NEG_INFINITY;

        for (g, genotype) in genotypes.iter().enumerate() {
            let log10_likelihood = -sample_pls[pl_offset + g] / 10.0;

            let log10_hw_prior =
                hardy_weinberg_log10_prior(genotype, allele_frequencies);

            let log10_posterior = log10_likelihood + log10_hw_prior;
            out_posteriors[post_offset + g] = log10_posterior;

            if log10_posterior > max_log10_posterior {
                max_log10_posterior = log10_posterior;
            }
        }

        // Normalize: convert from log10 to probabilities that sum to 1.
        // Subtract max first to avoid numerical overflow.
        let mut sum = 0.0;
        for g in 0..num_genotypes {
            let prob = f64::powf(10.0, out_posteriors[post_offset + g] - max_log10_posterior);
            out_posteriors[post_offset + g] = prob;
            sum += prob;
        }
        for g in 0..num_genotypes {
            out_posteriors[post_offset + g] /= sum;
        }
    }
}

/// Computes the Hardy-Weinberg log10 prior probability for a diploid genotype
/// given allele frequencies.
///
/// For a homozygous genotype (e.g., A/A): P = freq(A)^2
/// For a heterozygous genotype (e.g., A/T): P = 2 * freq(A) * freq(T)
///   (the factor of 2 accounts for the two possible arrangements: A|T and T|A)
fn hardy_weinberg_log10_prior(
    genotype: &DiploidGenotype,
    allele_frequencies: &[f64],
) -> f64 {
    let allele_a = genotype.allele_indices[0];
    let allele_b = genotype.allele_indices[1];
    let freq_a = allele_frequencies[allele_a];
    let freq_b = allele_frequencies[allele_b];

    let is_het = allele_a != allele_b;
    let combinatorial_factor = if is_het { 2.0 } else { 1.0 };

    f64::log10(combinatorial_factor * freq_a * freq_b)
}

/// M-step: update allele frequencies from genotype posteriors.
///
/// For each allele, counts how many copies are expected across all samples
/// (weighted by genotype posteriors), then combines with the Dirichlet prior
/// and normalizes to get new frequency estimates.
fn update_allele_frequencies(
    genotype_posteriors: &[f64],
    genotypes: &[DiploidGenotype],
    num_samples: usize,
    _num_alleles: usize,
    prior_pseudocounts: &[f64],
) -> Vec<f64> {
    let num_genotypes = genotypes.len();

    // Start with prior pseudocounts (our baseline expectation)
    let mut allele_counts: Vec<f64> = prior_pseudocounts.to_vec();

    // Add the expected allele counts from each sample
    for sample in 0..num_samples {
        let post_offset = sample * num_genotypes;

        for (g, genotype) in genotypes.iter().enumerate() {
            let posterior = genotype_posteriors[post_offset + g];

            // Each diploid genotype contributes 2 allele copies.
            // For hom 0/0: both copies go to allele 0.
            // For het 0/1: one copy to allele 0, one to allele 1.
            let allele_a = genotype.allele_indices[0];
            let allele_b = genotype.allele_indices[1];

            allele_counts[allele_a] += posterior;
            allele_counts[allele_b] += posterior;
        }
    }

    // Normalize counts to frequencies
    let total: f64 = allele_counts.iter().sum();
    allele_counts.iter().map(|count| count / total).collect()
}

/// Computes the QUAL score: how confident are we that the site is variable?
///
/// For each sample, we look at the posterior probability of being hom-ref (0/0).
/// The product of all those probabilities is the chance that ALL samples are
/// hom-ref (i.e., the site is not variable). QUAL is the phred-scaled version
/// of that probability: QUAL = -10 * log10(P(all samples are hom-ref)).
///
/// High QUAL means we are confident there is a real variant.
fn compute_qual(
    genotype_posteriors: &[f64],
    num_samples: usize,
    num_genotypes: usize,
) -> f64 {
    // The hom-ref genotype (0/0) is always the first one (pl_index = 0).
    let hom_ref_index = 0;

    let mut log10_prob_all_hom_ref = 0.0;
    for sample in 0..num_samples {
        let post_offset = sample * num_genotypes;
        let prob_hom_ref = genotype_posteriors[post_offset + hom_ref_index];
        log10_prob_all_hom_ref += prob_hom_ref.log10();
    }

    // QUAL = -10 * log10(P(no variant))
    -10.0 * log10_prob_all_hom_ref
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_enumerate_genotypes_biallelic() {
        // 2 alleles -> 3 genotypes: 0/0, 0/1, 1/1
        let gts = enumerate_diploid_genotypes(2);
        assert_eq!(gts.len(), 3);
        assert_eq!(gts[0].allele_indices, [0, 0]); // 0/0
        assert_eq!(gts[1].allele_indices, [0, 1]); // 0/1
        assert_eq!(gts[2].allele_indices, [1, 1]); // 1/1
    }

    #[test]
    fn test_enumerate_genotypes_triallelic() {
        // 3 alleles -> 6 genotypes: 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
        let gts = enumerate_diploid_genotypes(3);
        assert_eq!(gts.len(), 6);
        assert_eq!(gts[3].allele_indices, [0, 2]); // 0/2
        assert_eq!(gts[4].allele_indices, [1, 2]); // 1/2
    }

    #[test]
    fn test_prior_pseudocounts() {
        let prior = PriorConfig::default();
        let counts = compute_prior_pseudocounts(2, &prior);
        assert!((counts[0] - 10.0).abs() < 1e-9); // ref
        assert!((counts[1] - 0.01).abs() < 1e-9); // alt
    }

    #[test]
    fn test_clear_het_gets_high_qual() {
        // Two samples, both clearly het (PL: 0/0=99, 0/1=0, 1/1=99)
        //   PL=99 means "very unlikely", PL=0 means "most likely"
        let pls = vec![
            99.0, 0.0, 99.0, // sample 1: clearly 0/1
            99.0, 0.0, 99.0, // sample 2: clearly 0/1
        ];
        let result = estimate_posteriors(2, 2, &pls, &PriorConfig::default());

        // Both samples are clearly het, so QUAL should be high
        assert!(result.qual > 30.0, "QUAL={} should be > 30", result.qual);
        // Alt frequency is pulled toward ref by the strong prior (ref pseudocount=10),
        // so with only 2 samples it won't reach 0.5, but should be well above zero.
        assert!(
            result.allele_frequencies[1] > 0.05,
            "alt freq={} should be > 0.05",
            result.allele_frequencies[1]
        );
    }

    #[test]
    fn test_all_hom_ref_gets_low_qual() {
        // Two samples, both clearly hom-ref (PL: 0/0=0, 0/1=99, 1/1=99)
        let pls = vec![
            0.0, 99.0, 99.0, // sample 1: clearly 0/0
            0.0, 99.0, 99.0, // sample 2: clearly 0/0
        ];
        let result = estimate_posteriors(2, 2, &pls, &PriorConfig::default());

        // Both samples are hom-ref, so QUAL should be very low (no variant)
        assert!(result.qual < 1.0, "QUAL={} should be < 1", result.qual);
    }

    #[test]
    fn test_one_het_one_ref_moderate_qual() {
        // Sample 1: clearly het, sample 2: clearly hom-ref
        let pls = vec![
            99.0, 0.0, 99.0, // sample 1: clearly 0/1
            0.0, 99.0, 99.0, // sample 2: clearly 0/0
        ];
        let result = estimate_posteriors(2, 2, &pls, &PriorConfig::default());

        // One het out of two -> site is variable
        assert!(result.qual > 10.0, "QUAL={} should be > 10", result.qual);
        // Alt frequency should be positive but pulled down by the strong ref prior
        assert!(
            result.allele_frequencies[1] > 0.01 && result.allele_frequencies[1] < 0.5,
            "alt freq={} should be between 0.01 and 0.5",
            result.allele_frequencies[1]
        );
    }

    #[test]
    fn test_posteriors_sum_to_one() {
        let pls = vec![
            10.0, 0.0, 30.0, // sample 1: probably het
            0.0, 15.0, 40.0, // sample 2: probably hom-ref
        ];
        let result = estimate_posteriors(2, 2, &pls, &PriorConfig::default());

        for sample in 0..2 {
            let offset = sample * result.num_genotypes;
            let sum: f64 = result.genotype_posteriors
                [offset..offset + result.num_genotypes]
                .iter()
                .sum();
            assert!(
                (sum - 1.0).abs() < 1e-9,
                "sample {} posteriors sum to {} instead of 1.0",
                sample,
                sum
            );
        }
    }
}
