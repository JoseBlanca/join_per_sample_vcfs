================================================================================
ssr_hg002 — GIAB HG002 Tandem-Repeat benchmark (SSR gold standard, single sample)
================================================================================

WHAT THIS IS
------------
The GIAB HG002 Tandem-Repeat benchmark v1.0.1 (GRCh38), sliced to a 50,000-region
"Tier" subset, packaged so we can score the SSR caller (ssr-catalog / ssr-pileup /
ssr-call) against an ORTHOGONAL, assembly-based truth set.

It is ONE human sample (HG002), so it is a genuine gold standard for single-sample
SSR genotyping accuracy, but it does NOT exercise the cohort machinery. See
"WHAT IT CAN AND CANNOT EVALUATE" below and
doc/devel/specs/ssr_calling_models_and_optimization.md (§4, §6, §8).

Truth here is assembly-based (the VCF AD field is "coverage from maternal and
paternal assembly"), i.e. it does NOT share the short-read caller's stutter-vs-
allele failure mode — which is exactly what makes it usable as truth.


DIRECTORY LAYOUT
----------------
  README.txt                     this file
  catalog/                       genome-wide SSR locus catalog (OUR ssr-catalog output; INPUT)
  truth/                         GIAB assembly-based truth genotypes (GOLD; VCF)
  regions/                       benchmark confident regions (scoring mask; BED)
  bam/<cov>x/                    HG002 short reads at a ladder of coverages (INPUT reads)
  src/                           scripts (coverage subsampling, ...)

The reference FASTA is NOT duplicated here. It lives at:
  benchmarks/giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
That file's upper-cased-sequence md5 == 8b84395224048cb683894f8f4e237507, which is
exactly the catalog's `reference_md5` header (verified) — so the catalog binds to
it cleanly. Whole-FILE md5 differs (line-wrap width), but content is identical.


FILE-BY-FILE
------------
catalog/CA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.cat
  OUR ssr-catalog v0.1.0 output (BGZF; header is plain-text `## key: value`).
  Genome-wide: 515,352 loci, motif periods 2-6 (period 2: 139,926; 3: 66,126;
  4: 180,328; 5: 104,689; 6: 24,283). Built with trf-mod 4.10.0, flank_bp=50,
  min_purity=0.8, on the GRCh38 reference above.
  NOTE: this is GENOME-WIDE; the BAM/VCF/BED below cover only the 50k Tier subset,
  so evaluation intersects the catalog with the Tier regions.
  ROLE: locus definition, input to ssr-pileup / ssr-call.

truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz (+ .tbi)
  GIAB HG002 TR benchmark truth: 36,497 phased variant records over the Tier subset
  (bcftools view -R of the full TR VCF). GT is HG002's assembly-derived genotype;
  ~72% heterozygous (0|1 or 1|0: 26,373; 1|1: 10,069) — HG002 is outbred, so this
  set genuinely exercises heterozygote calling (unlike the tomato selfer cohort).
  INFO carries TRF motif/period annotations + population AF/HWE (NS=86 is the AF-
  annotation panel; only HG002 has a genotype column).
  ROLE: GOLD STANDARD to score our calls against.

regions/HG002_GRCh38_TandemRepeats_v1.0.1_Tier_50000.bed
regions/HG002_GRCh38_TandemRepeats_v1.0.1_Tier_50000_sort.bed.gz (+ .tbi)
  50,000 Tier1 high-confidence benchmark intervals (the .bed is the raw slice; the
  _sort.bed.gz is bgzipped+tabix'd for `-R` queries). Columns: chrom start end
  Tier tag ... .
  ROLE: scoring mask — restrict evaluation to where the truth is trustworthy.

bam/300x/HG002_TR_v1.0.1_Tier_300x.bam (+ .bai)
  HG002 300x Illumina HiSeq (148 bp reads, novoalign to GRCh38, PL:ILLUMINA),
  subset to the Tier BED. Coordinate-sorted + indexed (the raw download was
  labelled SO:unsorted and has been re-sorted here).
  ROLE: input reads at full depth.

bam/{50,30,20,15,10,5}x/HG002_TR_v1.0.1_Tier_<cov>x.bam (+ .bai)
  Coverage ladder subsampled from the 300x source (see src/subsample_coverages.sh).
  Fractions are computed from the MEASURED mean depth over the regions, seed 42,
  so the ladders are reproducible and nested (5x reads ⊂ 10x reads ⊂ ...).
  ROLE: probe the depth wall — how caller behaviour degrades as depth drops.
  (Measured source mean depth over the regions = 288.32x across 6,089,411 bp; the
  "300x" label is nominal. Subsample fractions are cov/288.32.)


HOW TO RUN THE PIPELINE (sketch)
--------------------------------
  REF=../giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
  CAT=catalog/CA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.cat
  # Stage 1 (per-sample): pileup reads at each coverage against the catalog
  ssr-pileup --catalog $CAT --reference $REF bam/30x/HG002_TR_v1.0.1_Tier_30x.bam -o hg002_30x.ssr.psp
  # Stage 2 (cohort, here n=1): call, sweeping the emission-model toggles
  ssr-call --catalog $CAT hg002_30x.ssr.psp -o hg002_30x.vcf
  # Score vs truth over the Tier regions (assembly truth = gold)
  # ... compare hg002_30x.vcf to truth/...vcf.gz within regions/...bed


WHAT IT CAN AND CANNOT EVALUATE
-------------------------------
CAN:
  - Absolute single-sample SSR genotype accuracy of the read model (Lr, "Model A")
    and the emission decision, against orthogonal assembly truth — real FP/FN and
    genotype-accuracy numbers the silver/HipSTR proxies cannot give.
  - Heterozygote sensitivity (this sample is ~72% het) — a complement to the
    tomato selfer cohort, which has almost no hets.
  - How accuracy degrades with depth (the coverage ladder).
  - End-to-end pipeline robustness on a second genome (human).

CANNOT:
  - The COHORT options that the optimization plan is actually about: the
    marginalized leave-one-out prior (MARG), the chemistry pre-pass / confident-
    genotype gate, cohort-recurrence, and the freebayes SFS-over-carrier-counts
    emission. At n=1 the marginalized prior is byte-identical to plug-in.
  - The low-depth lone-carrier-het tension of §4 at the population level — this is
    ONE sample. For that we still need either the silver tomato standard or a new
    multi-sample (population) gold standard; none was found (see the task note).

In short: this settles absolute single-sample accuracy of the base model; it does
NOT pick the tomato default regime (that stays a cohort question).


RESULTS + SCRIPTS (added 2026-07-09)
------------------------------------
Pipeline scripts in src/ (run from ssr_hg002/):
  add_m5_to_bams.sh          add @SQ M5 to each BAM header (pileup requires it)
  restrict_catalog_to_bed.py catalog -> catalog/HG002_Tier_restricted.cat (13,272 loci;
                             needed so the ssr-call burn-in samples COVERED loci)
  subsample_coverages.sh     the 5x..50x ladder from the 300x source
  run_ours_coverages.sh      ssr-pileup + ssr-call per coverage -> results/ours/{psp,vcf}
  run_hipstr_coverages.sh    HipSTR per coverage -> results/hipstr/ (--use-unpaired
                             --min-reads 5: the BAM slice orphaned mates; see script)
  run_freebayes_coverages.sh freebayes per coverage -> results/freebayes/ (general
                             Bayesian caller, NOT STR-aware; --targets the Tier regions)
  build_eval_table.py        -> results/eval_table.tsv (per-locus truth vs all 3 callers)
scripts/fp_fn_vs_coverage_dashboard.py   marimo dashboard over eval_table.tsv
  ( uv run marimo run benchmarks/ssr_hg002/scripts/fp_fn_vs_coverage_dashboard.py )

Truth-matching note: GIAB anchors an STR indel at the base BEFORE the tract
(REF=1bp anchor, ALT=anchor+units), so build_eval_table.py extends each truth
interval by the longest allele before overlapping it with a catalog locus —
without that, ~1000 correct calls per caller mis-score as false positives.

Headline (DETECTION: did the caller flag a length change; period-aware truth,
2,653 truth-positive loci of 13,272). Recall / precision at 300x -> 30x -> 10x -> 5x:
  ours       R 0.72 -> 0.55 -> 0.10 -> 0.003    P ~0.94 at every depth
  HipSTR     R 0.88 -> 0.87 -> 0.54 -> 0.17     P ~0.91
  freebayes  R 0.99 -> 0.95 -> 0.76 -> 0.50     P ~0.88 (0.94 under "any-indel" truth)

  Findings:
  - Our caller (run single-sample via the recurrence_k fallback) is PRECISE but
    recall-poor and collapses below ~15x — no cohort segregation signal to lean on.
    This is its single-sample FLOOR, not the cohort strength it is built for.
  - HipSTR (STR-specialist) dominates our single-sample mode on recall.
  - freebayes (general Bayesian caller, NOT STR-aware) has the BEST recall + F1 at
    every depth and degrades most gracefully. Its lower period-aware precision is
    largely a definitional artifact: it also calls real NON-period indels (precision
    rises to ~0.94 under the any-indel truth). The surprise: a general haplotype
    caller is an excellent SSR length-change DETECTOR here.
  CAVEAT: this metric is DETECTION (a length change was flagged), NOT genotype-length
  accuracy (the exact repeat count). freebayes emits an indel that need not equal the
  full tract change, so the STR-aware callers likely close the gap on exact-genotype
  concordance — that comparison is the natural next step.


PROVENANCE
----------
  Reference   GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set (content-md5
              8b84395224048cb683894f8f4e237507)
  Truth       GIAB HG002 TandemRepeats benchmark v1.0.1
  Reads       HG002.GRCh38.300x (novoalign) -> samtools view -L Tier.bed
  Catalog     ssr-catalog 0.1.0, trf-mod 4.10.0, flank 50, min_purity 0.8, 2026-07-09
