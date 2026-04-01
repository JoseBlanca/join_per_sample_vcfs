# Genotype joining
The overall objective is to replicate GATK's GenotypeGVCFs that joins the variants found in different gVCFs per sample files into a VCF with all samples.

We want to join the variants found in an OverlappingVariantGroup.

In general the phase should be maintained.
If there is per sample data with coverages these data should be also maintained in the final joined genotypes.

# Nomenclature
We are going to create two new terms in order to clarify the code:

- valelle: an allele found in a variant. Example: In this variant there are three valleles, AC,A,AG
- galelle: an allele found in a genotye. Example: This sample is diploid and homozygous, so it has two galleles (0/0), but only one vallele.

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

Algorithm
1. Group the variants by position because each position should be almost independently processed.

2. If remove_complex and the OverlappingVariantGroup covers more than one position and there's any alternative in any position other than the first one. This join will return an Error of kind complex_variant_removed.

3. for each position:
	-3.1 Rename valleles taking into account all valleles found in the per sample gVCFs. The order of the gVCFs given in the input should be respected. If we are inside of a deletion the vallele for the deleted position should be *.
    -3.2 Once we have the valleles we should change the galleles for each sample accordingly. These galleles will be tentative, the will be used as the input to the EM calculation.
    -3.3 Calculate new PLs using the EM algorithm.
    -3.4 Normalized the PLs, the best one should be 0.
    -3.5 From the new PLs create the final genotypes.
    -3.6 Calculate the genotype quality by getting the second better PL.
    -3.7 Calculate the variant quality.
