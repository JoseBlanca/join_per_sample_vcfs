#!/usr/bin/env bash
# run_ours_coverages.sh — run our SSR caller (ssr-pileup -> ssr-call) on every
# coverage rung, producing one VCF per coverage to score against the GIAB truth.
#
# HG002 is a SINGLE sample, so ssr-call runs as an n=1 cohort at each coverage
# (the cohort pre-pass/MARG machinery degenerates — see README "CANNOT evaluate").
#
# Run from ssr_hg002/:  bash src/run_ours_coverages.sh
# Env: BIN, THREADS, COVERAGES, EMIT_MODEL (PVC_SSR_EMIT_MODEL), MIN_QUAL.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REPO="$(cd "$ROOT/../.." && pwd)"
BIN="${BIN:-$REPO/target/release/pop_var_caller}"
REF="$ROOT/../giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
CAT="$ROOT/catalog/CA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.cat"
THREADS="${THREADS:-8}"
COVERAGES="${COVERAGES:-300 50 30 20 15 10 5}"

PSP_DIR="$ROOT/results/ours/psp"
VCF_DIR="$ROOT/results/ours/vcf"
LOG_DIR="$ROOT/results/logs"
mkdir -p "$PSP_DIR" "$VCF_DIR" "$LOG_DIR"

echo "binary : $BIN"
echo "ref    : $REF"
echo "catalog: $CAT"
echo "emit   : ${PVC_SSR_EMIT_MODEL:-heuristic (default)}"
echo

for cov in $COVERAGES; do
  bam="$ROOT/bam/${cov}x/HG002_TR_v1.0.1_Tier_${cov}x.bam"
  psp="$PSP_DIR/HG002_${cov}x.ssr.psp"
  vcf="$VCF_DIR/HG002_${cov}x.ssr.vcf"
  [[ -f "$bam" ]] || { echo "!! missing $bam — skipping ${cov}x" >&2; continue; }

  if [[ -s "$psp" ]]; then
    echo "[skip pileup] ${cov}x ($psp)"
  else
    echo "[ssr-pileup] ${cov}x -> $psp"
    /usr/bin/time -l "$BIN" ssr-pileup \
      --reference "$REF" --catalog "$CAT" --output "$psp" \
      --threads "$THREADS" --build-index-if-missing "$bam" \
      > "$LOG_DIR/pileup_${cov}x.log" 2>&1
  fi

  echo "[ssr-call]   ${cov}x -> $vcf"
  "$BIN" ssr-call \
    --catalog "$CAT" --output "$vcf" --threads "$THREADS" "$psp" \
    > "$LOG_DIR/call_${cov}x.log" 2>&1

  echo "   records: $(grep -vc '^#' "$vcf" 2>/dev/null || echo '?')"
done
echo "done."
