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
    /// Wright's fixation index (inbreeding coefficient).
    ///
    /// Controls the expected level of homozygosity in the population:
    ///   F = 0.0 — Hardy-Weinberg equilibrium (outcrossing species, animals)
    ///   F = 0.5 — partial selfing (e.g., some legumes, tomato)
    ///   F = 1.0 — complete selfing (all individuals are homozygous)
    ///
    /// Technically, F is the probability that the two alleles in a diploid
    /// individual are identical by descent.  Higher F means more homozygotes
    /// and fewer heterozygotes than HW expectation.
    ///
    /// For polyploids, F is applied as the probability that all allele copies
    /// in an individual are identical by descent.
    pub fixation_index: f64,
}

impl Default for PriorConfig {
    fn default() -> Self {
        PriorConfig {
            snp_heterozygosity: 0.001,
            indel_heterozygosity: 1.25e-4,
            heterozygosity_std_dev: 0.01,
            fixation_index: 0.0,
        }
    }
}

/// One possible genotype at a site, for any ploidy.
///
/// For diploid with alleles [A, T, C]:
///   - genotype 0/0 (AA) has allele_counts = [2, 0, 0]
///   - genotype 0/1 (AT) has allele_counts = [1, 1, 0]
///   - genotype 1/2 (TC) has allele_counts = [0, 1, 1]
///
/// For tetraploid with alleles [A, T]:
///   - genotype 0/0/0/0 has allele_counts = [4, 0]
///   - genotype 0/0/0/1 has allele_counts = [3, 1]
///   - genotype 0/0/1/1 has allele_counts = [2, 2]
///   - genotype 0/1/1/1 has allele_counts = [1, 3]
///   - genotype 1/1/1/1 has allele_counts = [0, 4]
pub struct Genotype {
    /// How many copies of each allele this genotype has.
    /// Length equals the number of alleles at the site.
    /// The values sum to the ploidy.
    pub allele_counts: Vec<usize>,
    /// Position of this genotype in the PL array (VCF ordering).
    pub pl_index: usize,
    /// Precomputed log10 of the multinomial coefficient for this genotype.
    /// This accounts for how many distinguishable arrangements exist.
    /// For diploid: het=log10(2), hom=log10(1)=0.
    /// For tetraploid AAAT: log10(4!/(3!*1!)) = log10(4).
    pub log10_multinomial_coefficient: f64,
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
    ploidy: usize,
    num_samples: usize,
    sample_pls: &[f64],
    prior: &PriorConfig,
) -> SitePosteriors {
    let genotypes = enumerate_genotypes(num_alleles, ploidy);
    let num_genotypes = genotypes.len();

    let prior_pseudocounts = compute_prior_pseudocounts(num_alleles, prior);

    // EM iteration
    let mut allele_frequencies = vec![1.0 / num_alleles as f64; num_alleles];
    let mut genotype_posteriors = vec![0.0; num_samples * num_genotypes];

    let f = prior.fixation_index;

    for _iteration in 0..50 {
        // E-step: compute each sample's genotype posteriors given current frequencies
        compute_genotype_posteriors(
            &allele_frequencies,
            &genotypes,
            sample_pls,
            num_samples,
            f,
            ploidy,
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
        f,
        ploidy,
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

/// Lists all possible genotypes for the given number of alleles and ploidy.
///
/// Diploid, 2 alleles: 0/0, 0/1, 1/1  (3 genotypes)
/// Diploid, 3 alleles: 0/0, 0/1, 1/1, 0/2, 1/2, 2/2  (6 genotypes)
/// Tetraploid, 2 alleles: 0/0/0/0, 0/0/0/1, 0/0/1/1, 0/1/1/1, 1/1/1/1  (5 genotypes)
///
/// The ordering matches the VCF PL field convention (sorted combinations with
/// replacement).
/// Enumerates genotypes in VCF PL order.
///
/// The VCF spec orders genotypes so that the highest allele index changes
/// slowest.  For diploid with 3 alleles:
///   0/0, 0/1, 1/1, 0/2, 1/2, 2/2
///
/// We generate non-decreasing allele sequences and sort them by the reversed
/// sequence (highest position first) to match this convention.
fn enumerate_genotypes(num_alleles: usize, ploidy: usize) -> Vec<Genotype> {
    let mut allele_sequences: Vec<Vec<usize>> = Vec::new();
    let mut current = vec![0usize; ploidy];

    collect_sorted_combinations(&mut current, 0, 0, num_alleles, ploidy, &mut allele_sequences);

    // Sort by VCF convention: compare reversed sequences lexicographically.
    // This makes the highest allele index the slowest-changing.
    allele_sequences.sort_by(|a, b| {
        for i in (0..ploidy).rev() {
            match a[i].cmp(&b[i]) {
                std::cmp::Ordering::Equal => continue,
                other => return other,
            }
        }
        std::cmp::Ordering::Equal
    });

    allele_sequences
        .into_iter()
        .enumerate()
        .map(|(pl_index, seq)| {
            let mut allele_counts = vec![0usize; num_alleles];
            for &a in &seq {
                allele_counts[a] += 1;
            }
            let log10_coeff = log10_multinomial_coefficient(ploidy, &allele_counts);
            Genotype {
                allele_counts,
                pl_index,
                log10_multinomial_coefficient: log10_coeff,
            }
        })
        .collect()
}

/// Collects all non-decreasing sequences of length `ploidy` from alleles 0..num_alleles.
fn collect_sorted_combinations(
    current: &mut [usize],
    pos: usize,
    min_allele: usize,
    num_alleles: usize,
    ploidy: usize,
    out: &mut Vec<Vec<usize>>,
) {
    if pos == ploidy {
        out.push(current.to_vec());
        return;
    }
    for allele in min_allele..num_alleles {
        current[pos] = allele;
        collect_sorted_combinations(current, pos + 1, allele, num_alleles, ploidy, out);
    }
}

/// Computes log10(ploidy! / (count_1! * count_2! * ...)).
///
/// This is the number of distinguishable arrangements of alleles in a genotype.
/// For diploid: hom (e.g., AA) = 1, het (e.g., AT) = 2.
/// For tetraploid AABT: 4!/(2!*1!*1!) = 12.
fn log10_multinomial_coefficient(ploidy: usize, allele_counts: &[usize]) -> f64 {
    let mut log10_result = log10_factorial(ploidy);
    for &count in allele_counts {
        log10_result -= log10_factorial(count);
    }
    log10_result
}

/// Computes log10(n!) using the simple sum log10(1) + log10(2) + ... + log10(n).
/// For the small values of ploidy we deal with, this is fast enough.
fn log10_factorial(n: usize) -> f64 {
    (1..=n).map(|i| (i as f64).log10()).sum()
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
/// 2. The genotype prior given current allele frequencies and fixation index
/// 3. A combinatorial factor (how many distinguishable allele arrangements)
///
/// The results are normalized so that each sample's posteriors sum to 1.
fn compute_genotype_posteriors(
    allele_frequencies: &[f64],
    genotypes: &[Genotype],
    sample_pls: &[f64],
    num_samples: usize,
    fixation_index: f64,
    ploidy: usize,
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
                genotype_log10_prior(genotype, allele_frequencies, fixation_index, ploidy);

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

/// Computes the log10 genotype prior probability, incorporating the fixation
/// index F (Wright's inbreeding coefficient).
///
/// The prior is a mixture of two components:
///
/// 1. **Hardy-Weinberg component** (weight = 1 - F):
///    The standard random-mating expectation.
///    P_hw = multinomial_coeff * product(freq[a] ^ count[a])
///
/// 2. **Inbred component** (weight = F):
///    All allele copies in the individual are identical by descent.
///    Only fully homozygous genotypes (all copies the same allele) have
///    nonzero probability under this component: P_inbred(all_A) = freq(A).
///
/// The combined prior is:
///   P = (1 - F) * P_hw  +  F * P_inbred
///
/// When F = 0 this reduces to pure Hardy-Weinberg (outcrossing/animals).
/// When F = 1 only homozygous genotypes have nonzero prior (complete selfing).
fn genotype_log10_prior(
    genotype: &Genotype,
    allele_frequencies: &[f64],
    fixation_index: f64,
    ploidy: usize,
) -> f64 {
    // Hardy-Weinberg component
    let mut log10_hw = genotype.log10_multinomial_coefficient;
    for (allele, &count) in genotype.allele_counts.iter().enumerate() {
        if count > 0 {
            log10_hw += count as f64 * allele_frequencies[allele].log10();
        }
    }
    let prob_hw = f64::powf(10.0, log10_hw);

    // Inbred component: nonzero only for fully homozygous genotypes
    // (one allele has count == ploidy, all others are 0)
    let prob_inbred = genotype
        .allele_counts
        .iter()
        .enumerate()
        .find(|(_, count)| **count == ploidy)
        .map(|(allele, _)| allele_frequencies[allele])
        .unwrap_or(0.0);

    // Mix the two components
    let combined = (1.0 - fixation_index) * prob_hw + fixation_index * prob_inbred;

    // Guard against log10(0)
    if combined <= 0.0 {
        -1e10
    } else {
        combined.log10()
    }
}

/// M-step: update allele frequencies from genotype posteriors.
///
/// For each allele, counts how many copies are expected across all samples
/// (weighted by genotype posteriors), then combines with the Dirichlet prior
/// and normalizes to get new frequency estimates.
fn update_allele_frequencies(
    genotype_posteriors: &[f64],
    genotypes: &[Genotype],
    num_samples: usize,
    _num_alleles: usize,
    prior_pseudocounts: &[f64],
) -> Vec<f64> {
    let num_genotypes = genotypes.len();

    // Start with prior pseudocounts (our baseline expectation)
    let mut allele_counts: Vec<f64> = prior_pseudocounts.to_vec();

    // Add the expected allele counts from each sample.
    // Each genotype contributes `ploidy` allele copies, distributed
    // according to its allele_counts (e.g., tetraploid AAAT contributes
    // 3 copies of A and 1 copy of T).
    for sample in 0..num_samples {
        let post_offset = sample * num_genotypes;

        for (g, genotype) in genotypes.iter().enumerate() {
            let posterior = genotype_posteriors[post_offset + g];

            for (allele, &count) in genotype.allele_counts.iter().enumerate() {
                if count > 0 {
                    allele_counts[allele] += posterior * count as f64;
                }
            }
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
    fn test_enumerate_genotypes_diploid_biallelic() {
        // Diploid, 2 alleles -> 3 genotypes: 0/0, 0/1, 1/1
        let gts = enumerate_genotypes(2, 2);
        assert_eq!(gts.len(), 3);
        assert_eq!(gts[0].allele_counts, [2, 0]); // 0/0
        assert_eq!(gts[1].allele_counts, [1, 1]); // 0/1
        assert_eq!(gts[2].allele_counts, [0, 2]); // 1/1
    }

    #[test]
    fn test_enumerate_genotypes_diploid_triallelic() {
        // Diploid, 3 alleles -> 6 genotypes
        let gts = enumerate_genotypes(3, 2);
        assert_eq!(gts.len(), 6);
        assert_eq!(gts[3].allele_counts, [1, 0, 1]); // 0/2
        assert_eq!(gts[4].allele_counts, [0, 1, 1]); // 1/2
    }

    #[test]
    fn test_enumerate_genotypes_tetraploid_biallelic() {
        // Tetraploid, 2 alleles -> 5 genotypes:
        //   0/0/0/0, 0/0/0/1, 0/0/1/1, 0/1/1/1, 1/1/1/1
        let gts = enumerate_genotypes(2, 4);
        assert_eq!(gts.len(), 5);
        assert_eq!(gts[0].allele_counts, [4, 0]); // 0/0/0/0
        assert_eq!(gts[1].allele_counts, [3, 1]); // 0/0/0/1
        assert_eq!(gts[2].allele_counts, [2, 2]); // 0/0/1/1
        assert_eq!(gts[3].allele_counts, [1, 3]); // 0/1/1/1
        assert_eq!(gts[4].allele_counts, [0, 4]); // 1/1/1/1
    }

    #[test]
    fn test_multinomial_coefficients() {
        // Diploid het: 2!/(1!*1!) = 2
        let gts = enumerate_genotypes(2, 2);
        assert!((gts[1].log10_multinomial_coefficient - 2.0_f64.log10()).abs() < 1e-9);
        // Diploid hom: 2!/(2!*0!) = 1
        assert!((gts[0].log10_multinomial_coefficient - 0.0).abs() < 1e-9);

        // Tetraploid 0/0/0/1: 4!/(3!*1!) = 4
        let gts4 = enumerate_genotypes(2, 4);
        assert!((gts4[1].log10_multinomial_coefficient - 4.0_f64.log10()).abs() < 1e-9);
        // Tetraploid 0/0/1/1: 4!/(2!*2!) = 6
        assert!((gts4[2].log10_multinomial_coefficient - 6.0_f64.log10()).abs() < 1e-9);
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
        let result = estimate_posteriors(2, 2, 2, &pls, &PriorConfig::default());

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
        let result = estimate_posteriors(2, 2, 2, &pls, &PriorConfig::default());

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
        let result = estimate_posteriors(2, 2, 2, &pls, &PriorConfig::default());

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
        let result = estimate_posteriors(2, 2, 2, &pls, &PriorConfig::default());

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

    #[test]
    fn test_tetraploid_clear_variant() {
        // Tetraploid, biallelic: 5 genotypes (0/0/0/0 .. 1/1/1/1)
        // Sample 1: clearly 0/0/1/1 (PL index 2)
        // Sample 2: clearly 0/0/0/0 (PL index 0)
        let pls = vec![
            99.0, 99.0, 0.0, 99.0, 99.0, // sample 1: 0/0/1/1
            0.0, 99.0, 99.0, 99.0, 99.0,  // sample 2: 0/0/0/0
        ];
        let result = estimate_posteriors(2, 4, 2, &pls, &PriorConfig::default());

        assert_eq!(result.num_genotypes, 5);
        assert!(result.qual > 10.0, "QUAL={} should be > 10", result.qual);
    }

    #[test]
    fn test_tetraploid_posteriors_sum_to_one() {
        let pls = vec![
            10.0, 5.0, 0.0, 20.0, 40.0, // sample 1
            0.0, 10.0, 20.0, 30.0, 40.0, // sample 2
        ];
        let result = estimate_posteriors(2, 4, 2, &pls, &PriorConfig::default());

        for sample in 0..2 {
            let offset = sample * result.num_genotypes;
            let sum: f64 = result.genotype_posteriors
                [offset..offset + result.num_genotypes]
                .iter()
                .sum();
            assert!(
                (sum - 1.0).abs() < 1e-9,
                "tetraploid sample {} posteriors sum to {} instead of 1.0",
                sample,
                sum
            );
        }
    }

    // --- Fixation index (F) tests ---

    fn prior_with_f(f: f64) -> PriorConfig {
        PriorConfig {
            fixation_index: f,
            ..PriorConfig::default()
        }
    }

    #[test]
    fn test_high_f_favors_homozygotes() {
        // Ambiguous data: PL equally supports het and hom-alt
        let pls = vec![
            30.0, 0.0, 0.0, // sample: het and hom-alt equally likely from reads
        ];

        let result_outcrossing = estimate_posteriors(2, 2, 1, &pls, &prior_with_f(0.0));
        let result_selfing = estimate_posteriors(2, 2, 1, &pls, &prior_with_f(0.95));

        // With high F, the hom-alt posterior should be higher than het
        let het_idx = 1;
        let hom_alt_idx = 2;

        let het_post_outcrossing = result_outcrossing.genotype_posteriors[het_idx];
        let het_post_selfing = result_selfing.genotype_posteriors[het_idx];

        let hom_alt_post_outcrossing = result_outcrossing.genotype_posteriors[hom_alt_idx];
        let hom_alt_post_selfing = result_selfing.genotype_posteriors[hom_alt_idx];

        // Selfing should reduce het probability relative to outcrossing
        assert!(
            het_post_selfing < het_post_outcrossing,
            "het posterior with F=0.95 ({}) should be less than F=0.0 ({})",
            het_post_selfing,
            het_post_outcrossing
        );

        // Selfing should increase hom-alt probability relative to outcrossing
        assert!(
            hom_alt_post_selfing > hom_alt_post_outcrossing,
            "hom-alt posterior with F=0.95 ({}) should be greater than F=0.0 ({})",
            hom_alt_post_selfing,
            hom_alt_post_outcrossing
        );
    }

    #[test]
    fn test_f_zero_equals_hardy_weinberg() {
        // F=0 should give the same result as the default (which is F=0)
        let pls = vec![
            10.0, 0.0, 30.0,
            0.0, 15.0, 40.0,
        ];
        let result_default = estimate_posteriors(2, 2, 2, &pls, &PriorConfig::default());
        let result_f0 = estimate_posteriors(2, 2, 2, &pls, &prior_with_f(0.0));

        for (a, b) in result_default
            .genotype_posteriors
            .iter()
            .zip(result_f0.genotype_posteriors.iter())
        {
            assert!((a - b).abs() < 1e-12, "F=0 should equal default HW");
        }
    }

    #[test]
    fn test_f_one_gives_zero_het_prior() {
        // With F=1 (complete selfing) and clear het data, the prior fights
        // against het, so the het posterior should be much lower than with F=0
        let pls = vec![
            99.0, 0.0, 99.0, // sample: clearly het from reads
        ];

        let result_outcrossing = estimate_posteriors(2, 2, 1, &pls, &prior_with_f(0.0));
        let result_full_selfing = estimate_posteriors(2, 2, 1, &pls, &prior_with_f(1.0));

        let het_idx = 1;
        let het_outcrossing = result_outcrossing.genotype_posteriors[het_idx];
        let het_selfing = result_full_selfing.genotype_posteriors[het_idx];

        assert!(
            het_selfing < het_outcrossing,
            "het posterior with F=1.0 ({}) should be < F=0.0 ({})",
            het_selfing,
            het_outcrossing
        );
    }
}
