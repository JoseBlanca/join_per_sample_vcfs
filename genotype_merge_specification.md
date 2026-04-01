# Per sample gVCFs

Let's imagine that we have two samples: s1 and s2.

pos 12
ref AG
s1  A-
s2  AC

The gVCF files would be something like:
gVCF s1
pos ref  alt           gt  pl
1   AG   A,<NON_REF>   1/1 200,30,0,200,30,200
2   G    <NON_REF>     0/0 0,0,0

gVCF s2
1 A    <NON_REF>     0/0 0,30,200
2 G    C,<NON_REF>   1/1 200,30,0,200,30,200

# allele, vallele, gallele

vallele: one of the possible alleles found in a Variant
gallele: one of the alleles found in a genotype for a sample, so the allele found in one of the chromosomes

# PL meaning

Each possible genotype has one PL, its likelihood relative to the most likely genotype. (The PL for the most likely genotype is always 0).

"Normalized" Phred-scaled likelihoods of the possible genotypes. For the typical case of a monomorphic site (where there is only one ALT allele) in a diploid organism, the PL field will contain three numbers, corresponding to the three possible genotypes (0/0, 0/1, and 1/1). The PL values are "normalized" so that the PL of the most likely genotype (assigned in the GT field) is 0 in the Phred scale. We use "normalized" in quotes because these are not probabilities. We set the most likely genotype PL to 0 for easy reading purpose. The other values are scaled relative to this most likely genotype.

Keep in mind, if you are not familiar with the statistical lingo, that when we say PL is the "Phred-scaled likelihood of the genotype", we mean it is "How much less likely that genotype is compared to the best one". Have a look at this article (https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs) for an example of how PL is calculated.

The basic formula for calculating PL is:

PL = -10 * \log{P(Genotype | Data)}

where P(Genotype | Data) is the conditional probability of the Genotype given the sequence Data that we have observed. 

The PL calculation from the allele counts done by the HaplotypeCaller is described here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890511-Assigning-per-sample-genotypes-HaplotypeCaller

When you run HaplotypeCaller with -ERC GVCF to produce a gVCF, there is an additional calculation to determine the genotype likelihoods associated with the symbolic <NON-REF> allele (which represents the possibilities that remain once you’ve eliminated the REF allele and any ALT alleles that are being evaluated explicitly).

Calculation of PL and GQ by HaplotypeCaller and GenotypeGVCFs
https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs

I've seen an example of a real gVCF file an when there's a deletion all PLs are set to zero.
SL4.0ch00       592829  .       CG      C,<NON_REF>     107.78  .       .     GT:AD:DP:GQ:PGT:PID:PL:PS:SB    1|1:0,3,0:3:9:0|1:592829_CG_C:121,9,0,121,9,121:592829:0,0,3,0
SL4.0ch00       592830  .       G       <NON_REF>       .       .       .       GT:AD:DP:GQ:PL  0/0:0,3:3:0:0,0,0


# Genotype quality (GQ)

The value of GQ is simply the difference between the second lowest PL and the lowest PL (which is always 0). So, in our example GQ = 20 - 0 = 20. Note that the value of GQ is capped at 99 for practical reasons, so even if the calculated GQ is higher, the value emitted to the VCF will be 99.

# GATK joint VCF

GTAK GenotypeGVCFs (https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs) would create the following joint VCF
VCF
pos ref alt s1   s2
1   AG  A   1/1  0/0
2   G   C,* 2/2  1/1

GenotypeGVCFs creates one output position for each input position.

When a position for one sample is in fact a deletion GenotypeGVCFs puts an * as the allele.
The meaning of the * allele is explained here: https://gatk.broadinstitute.org/hc/en-us/articles/360035530752-What-types-of-variants-can-GATK-tools-detect-or-handle
the * allele used to signify the presence of a spanning deletion, or undefined events like a very large allele or one that's fuzzy and not fully modeled; i.e. there's some event going on here but we don't know what exactly)

# Our desired merged VCF

We want do merge the overlapping alleles into one variant.

pos 12
ref AG
s1  A-
s2  AC

VCF ref alt  s1  s2
1   AG  A,AC 1/1 2/2

# PL calculation of the merged genotypes

The final PL values are calculated from the PL values found in the per sample gVCF files.
We need to keep track of the values associated to each allele in the original per sample gVCF files.

gVCF s1
pos ref  alt           gt  pl
1   AG   A,<NON_REF>   1/1 pl1_s1_00,pl1_s1_01,pl1_s1_11,pl1_s1_0nr,pl1_s1_1nr,pl1_s1_nrnr
2   G    <NON_REF>     0/0 pl2_s1_00,pl2_s1_0nr,pl2_s1_nrnr
s1_11,pl1_s1_0nr = 0
pl2_s1_00=0
pl2_s1_0nr=0
pl2_s1_nrnr=0

var_s1_1.pls is a vector with the values:pl1_s1_00,pl1_s1_01,pl1_s1_11,pl1_s1_0nr,pl1_s1_1nr,pl1_s1_nrnr

gVCF s2
pos ref  alt           gt  pl
1   A    <NON_REF>     0/0 pl1_s2_00,pl1_s2_0nr,pl1_s2_nrnr
2   G    C,<NON_REF>   1/1 pl2_s2_00,pl2_s2_01,pl2_s2_11,pl2_s2_0nr,pl2_s2_1nr,pl2_s2_nrnr
pl1_s2_00=0
pl2_s2_11=0

When we are building the merged genotypes out of the original per gVCF genotypes we have to store the PLs that accompanied every original genotype for each sample.

## Building the galleles

Recipe to create the galleles and genotypes
merged_gallele_s1_c1 = gallele_s1_pos1_c1 + gallele_s1_pos2_c2
merged_gallele_s1_c2 = gallele_s1_pos1_c1 + gallele_s1_pos2_c2
Meaning of the nomenclature: ga_s1_pos1_c1: gallele, sample1, pos 1, chrom 1
If there is a deletion the galleles covering the deletion will be "", meaning that those galleles per sample won't be used to create the merged galleles.

We are assuming known phase or not more than one heterozygotic position, otherwise we won't solve those valleles, they will be missing (. or -1).

The valleles will be the ordered set of the reference vallele and then the galleles for all samples.
galleles = [merged_gallele_ref, merged_gallele_s1_c1, merged_gallele_s1_c2, merged_gallele_s2_c1, merged_gallele_s2_c2]
The merged valleles should be obtained from the vector of galleles removing the repeated ones
valleles = non_repeated(galleles)
vallele_idx = position of vallele in valleles
There's a mapping from gallele to vallele in which several galleles could correspond to the same vallele -> vallele_idx_for_gallele

The sample genotypes for sample x to this point are coded as: merged_gallele_sx_c1/merged_gallele_sx_c2
To get the final merged numeric genotype we have to recode as: vallele_idx_for_gallele[merged_gallele_sx_c1]/vallele_idx_for_gallele[merged_gallele_sx_c12]

A PL should be calculated for each possible combination of valleles in each sample.

1 -> 1:0/0 1:0/1 1:1/1
2 -> 2:0/0 2:0/1 2:1/1

PLs for all possible merged genotypes:
PL1:0/0 + PL2:0/0
PL1:0/0 + PL2:0/1
PL1:0/0 + PL2:1/1
PL1:0/1 + PL2:0/0
PL1:0/1 + PL2:0/1
PL1:0/1 + PL2:1/1
PL1:1/1 + PL2:0/0
PL1:1/1 + PL2:0/1
PL1:1/1 + PL2:1/1
nomenclature: PL1:0/0 means PL for position 0 for genotype 0/0
The genotype qual is calculated from the difference between the minimum merged PL and the second one.
After that we remove all combinations with non-ref.
Next we normalize the values.

#

In each position (i) there is only one genotype.
how much we trust genotype in position i = 1 - prob of any of the all other genotypes
This value is related to GQ
