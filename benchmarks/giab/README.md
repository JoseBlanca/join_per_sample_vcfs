Human benchmarks based on the Genome in a Bottle dataset (GIAB)

In these directories there are files for two different datasets:

1. per_sample: SNP calling of three samples/individual for 100 regions. 
2. mendelian: SNP calling of a trio family for 100 regions.

per_sample dataset:

- The regions were selected at random
- The regions are different for each sample.
- Directory: benchmarks/giab/per_sample
- The vcf dir includes the truth downloaded from GIAB
- In the bam dir there are BAMs with different coverages

mendelian dataset:
- The regions were selected at random
- The regions are the same for all samples.
- Directory: benchmarks/giab/mendelian
- The vcf dir includes the truth downloaded from GIAB
- HG002 is the child and HG003 and HG004 are the parents
