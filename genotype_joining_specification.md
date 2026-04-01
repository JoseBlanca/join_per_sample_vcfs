# Genotype joining
The overall objective is to replicate GATK's GenotypeGVCFs that joins the variants found in different gVCFs per sample files into a VCF with all samples.

We want to join the variants found in an OverlappingVariantGroup.
The result of processing an OverlappingVariantGroup could be:

- An error if the variants can not be joined.
- A vector of Variants with one or more joined variants.

In general the phase should be maintained.
If there is per sample data with coverages these data should be also maintained in the final joined genotypes.

The output FORMAT fields should be: GT, GQ, PL, DP, AD.

# Nomenclature
We are going to create two new terms in order to clarify the code:

- vallele: an allele found in a variant. Example: In this variant there are three valleles, AC,A,AG
- gallele: an allele found in a genotype. Example: This sample is diploid and homozygous, so it has two galleles (0/0), but only one vallele.

# Examples of Genotype joining

## Example: simple SNP
gVCF sample1
pos REF ALT GT_sample1
1   A   T   1/1

gVCF sample2
pos REF ALT GT_sample2
1   A   C   0/1

joined VCF
pos REF ALT GT_sample1 GT_sample2
1   A   T,C 1/1        0/2

## Example: Mixed:Deletion + SNP
gVCF sample1
pos REF ALT GT_sample1
1   AC  A   1/1
2   C   .   0/0

gVCF sample2
pos REF ALT GT_sample2
1   A   C   0/1
2   C   .   0/0

joined VCF
pos REF ALT  GT_sample1 GT_sample2
1   AC  A,CC 1/1        0/2
2   C   .    0/0        0/0

## Example: Complex:Deletion + SNP in region covered by deletion
gVCF sample1
pos REF ALT GT_sample1
1   AC  A   1/1
2   C   .   0/0

gVCF sample2
pos REF ALT GT_sample2
1   A   .   0/0
2   C   G   1/1

By default we want to remove these kinds of OverlappingVariantGroup. This is not what GATK does. We will have the option of doing what GATK does, but only if the user requests to do it that way explicitly. The parameter should be remove_complex and by default should be true.

joined VCF GATK style
pos REF ALT  GT_sample1 GT_sample2
1   AC  A    1/1        0/0
2   C   *,G  1/1        2/2

# Algorithm

## 1. Group variants by position

Group the variants within the OverlappingVariantGroup by their start position, because each position should be processed almost independently.
Track which samples have active deletions from earlier positions (deletion state carries forward across positions).

## 2. Complex variant detection

If remove_complex: scan all positions. If the OverlappingVariantGroup covers more than one position and any non-first position has a real (non-*, non-ref) ALT allele, return an Error of kind complex_variant_removed.

The rationale: the `*` (spanning deletion) allele is expected in GATK output for deletions and is not inherently complex. The complexity arises when a real variant (e.g., a SNP) co-occurs at a position that is also covered by a deletion from another sample.

## 3. Process each position

### 3.a Build the unified vallele list

Collect all alleles from all samples at this position:
- Strip `<NON_REF>` from the output allele list (but retain the `<NON_REF>` PL values for use in step 3.b).
- For samples inside an active deletion (tracked from step 1), their vallele is `*`.
- Build the unified vallele list: [ref, alt1, alt2, ...].
  The order of the gVCFs given in the input should be respected when adding new alleles.
- At the first position of a multi-position group (where deletions start), if the ref alleles have different lengths, extend shorter alleles by appending trailing ref bases. For example, if REF=AC and a sample has A->C (SNP), the SNP allele becomes CC (extended to match the 2-base ref).
- Multi-allelic variants within a single sample's gVCF are supported: they just contribute more alleles to the unified set.

### 3.b Remap each sample's PLs to the unified vallele set

This is a critical step. Each sample's original PLs cover only the genotypes for its own allele set. The EM needs PLs for all genotypes of the unified allele set.

For each sample:
- Map the sample's original allele indices to unified allele indices.
- For genotypes involving alleles the sample didn't have, use the `<NON_REF>` likelihood values if available. If no `<NON_REF>` data exists, use a high PL penalty (e.g., the max PL from that sample or a fixed high value).
- Samples with no variant at this position (they are in a reference block or absent from the group): use synthetic hom-ref PLs, or remap from the reference block's PLs if the reference block carries PL data.
- Samples with GT-only (no PL field): generate synthetic PLs from the GT call (as already implemented: 99% confidence on the called genotype, remainder spread uniformly).

### 3.c Assign tentative galleles

Once we have the unified valleles, change the galleles for each sample accordingly (mapping original allele indices to unified indices). These galleles are tentative and will be used as fallback; the final galleles come from the EM posteriors.

### 3.d Run the EM algorithm

Run the EM algorithm on the remapped PLs to estimate allele frequencies and genotype posteriors. This uses the existing `estimate_posteriors()` function.

### 3.e Derive output fields from posteriors

- Normalize the posteriors to phred scale: PP = -10 * log10(posterior). Normalize so that the best (most likely genotype) = 0.
- GT = genotype with PP = 0 (highest posterior).
- GQ = second-lowest PP value (confidence in the GT call).
- PL = the remapped input likelihoods from step 3.b (this is what GATK outputs as PL; it preserves the original evidence).
- QUAL = phred-scaled probability that the site is non-variable (all samples hom-ref).

### 3.f Prune unused valleles

After the EM assigns final genotypes, some valleles might no longer be used by any sample (e.g., the EM reassigned all carriers of a rare allele to a different genotype). These unused valleles should be removed from the allele list, and all gallele indices should be renumbered accordingly. PL values should also be recomputed for the reduced allele set.

### 3.g Skip invariant positions

If the site is non-variable after pruning (no sample has a non-ref GT, i.e., all samples are hom-ref), do not emit a variant for this position.

### 3.h Update deletion tracking

Update the deletion state for each sample/haplotype so the next position knows which samples are inside an active deletion.