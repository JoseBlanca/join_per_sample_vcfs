#!/usr/bin/env bash
# run_hipstr_coverages.sh — run HipSTR (single-sample native) on every coverage
# rung over the SAME loci our caller uses, for an apples-to-apples FP/FN comparison.
#
# Two non-default flags matter here:
#   --use-unpaired : the BAMs were sliced with `samtools view -L Tier.bed`, which
#                    physically removed mates mapping outside the small Tier intervals.
#                    HipSTR discards unpaired reads by default -> ~99% of loci skipped.
#                    --use-unpaired makes it genotype from the (now single-end) reads.
#   --min-reads 5  : HipSTR's default (100) is a cohort threshold; for a single-sample
#                    depth study we lower it so HipSTR attempts genotyping at low depth.
#
# Run from ssr_hg002/:  bash src/run_hipstr_coverages.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
HIPSTR="$ROOT/../../HipSTR/HipSTR"
REF="$ROOT/../giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
REG="$ROOT/results/hipstr/HG002_Tier.hipstr_regions.bed"
MIN_READS="${MIN_READS:-5}"
COVERAGES="${COVERAGES:-300 50 30 20 15 10 5}"
OUT="$ROOT/results/hipstr"; LOG="$ROOT/results/logs"
mkdir -p "$OUT" "$LOG"

for cov in $COVERAGES; do
  bam="$ROOT/bam/${cov}x/HG002_TR_v1.0.1_Tier_${cov}x.bam"
  vcf="$OUT/HG002_${cov}x.str.vcf.gz"
  echo "[hipstr] ${cov}x -> $vcf"
  "$HIPSTR" --bams "$bam" --bam-samps HG002 --bam-libs HG002 \
    --fasta "$REF" --regions "$REG" --use-unpaired --min-reads "$MIN_READS" \
    --str-vcf "$vcf" > "$LOG/hipstr_${cov}x.log" 2>&1 || { echo "  FAILED (see log)"; continue; }
  echo "   records: $(zcat < "$vcf" 2>/dev/null | grep -vc '^#')"
done
echo "done."
