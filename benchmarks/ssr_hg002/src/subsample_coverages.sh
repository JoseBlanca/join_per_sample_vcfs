#!/usr/bin/env bash
#
# subsample_coverages.sh — build a coverage ladder from the 300x source BAM.
#
# The GIAB HG002 Tier BAM is ~300x over the benchmark regions. To probe how the
# SSR caller behaves as depth drops (the depth wall of
# doc/devel/specs/ssr_calling_models_and_optimization.md §8), we subsample it to
# a ladder of lower coverages with `samtools view -s`.
#
# Fractions are computed from the MEASURED mean depth of the source over the
# benchmark regions (not the nominal 300x), so the labels are honest. A single
# fixed seed keeps every run reproducible AND keeps the ladders nested (the 5x
# reads are a subset of the 10x reads, etc.) because samtools hashes read names.
#
# Run from the ssr_hg002/ directory root:  bash src/subsample_coverages.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/bam/300x/HG002_TR_v1.0.1_Tier_300x.bam"
BED="$ROOT/regions/HG002_GRCh38_TandemRepeats_v1.0.1_Tier_50000.bed"
SEED=42
THREADS=8
TARGETS=(50 30 20 15 10 5)

[[ -f "$SRC" ]] || { echo "source BAM missing: $SRC" >&2; exit 1; }

echo ">> measuring source mean depth over the benchmark regions ..."
SRC_COV=$(samtools depth -a -b "$BED" "$SRC" \
  | awk '{s+=$3; n++} END{printf "%.4f", (n? s/n : 0)}')
echo ">> source mean depth = ${SRC_COV}x"

for cov in "${TARGETS[@]}"; do
  out="$ROOT/bam/${cov}x/HG002_TR_v1.0.1_Tier_${cov}x.bam"
  # fraction of the source to keep; samtools -s wants seed.frac as one float.
  frac=$(awk -v t="$cov" -v s="$SRC_COV" 'BEGIN{ f=t/s; if(f>=1) f=0.999999; printf "%.6f", f }')
  echo ">> ${cov}x  (fraction ${frac})  -> $out"
  samtools view -@ "$THREADS" -b -s "${SEED}${frac#0}" -o "$out" "$SRC"
  samtools index -@ "$THREADS" "$out"
done

echo ">> done. resulting ladder:"
for cov in 300 "${TARGETS[@]}"; do
  b="$ROOT/bam/${cov}x/HG002_TR_v1.0.1_Tier_${cov}x.bam"
  [[ -f "$b" ]] && printf "   %4sx  %s\n" "$cov" "$(du -h "$b" | cut -f1)"
done
