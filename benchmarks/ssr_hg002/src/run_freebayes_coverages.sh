#!/usr/bin/env bash
# run_freebayes_coverages.sh — run freebayes (general Bayesian caller, NOT STR-aware)
# on every coverage rung, restricted to the Tier regions, as a third caller in the
# ours-vs-HipSTR-vs-truth comparison. freebayes represents SSR length changes as
# ordinary indels; the point is to see where its haplotype model wins anyway.
#
# freebayes emits records in --targets order, so targets.bed must be sorted (it is;
# built from the sorted Tier bed.gz). It does NOT need @SQ M5. --ploidy 2 (HG002).
# The low-QUAL tail is left in the raw VCF; the eval table applies a QUAL gate.
#
# Run from ssr_hg002/:  bash src/run_freebayes_coverages.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REF="$ROOT/../giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
TARGETS="$ROOT/results/freebayes/targets.bed"
COVERAGES="${COVERAGES:-300 50 30 20 15 10 5}"
OUT="$ROOT/results/freebayes"; LOG="$ROOT/results/logs"
mkdir -p "$OUT" "$LOG"

[[ -s "$TARGETS" ]] || { echo "missing $TARGETS (build from the sorted Tier bed.gz)" >&2; exit 1; }

for cov in $COVERAGES; do
  bam="$ROOT/bam/${cov}x/HG002_TR_v1.0.1_Tier_${cov}x.bam"
  vcf="$OUT/HG002_${cov}x.fb.vcf"
  echo "[freebayes] ${cov}x -> ${vcf}.gz"
  freebayes -f "$REF" -t "$TARGETS" --ploidy 2 "$bam" > "$vcf" 2>"$LOG/freebayes_${cov}x.log"
  bgzip -f "$vcf" && tabix -f -p vcf "${vcf}.gz"
  echo "   records: $(zcat < "${vcf}.gz" | grep -vc '^#')"
done
echo "done."
