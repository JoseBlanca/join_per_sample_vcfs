#!/usr/bin/env bash
# Phase A (part 2): joint-call the 59 tomato2 .psp into one cohort VCF.
#
# RESEARCH CALLSET — maximum retention. We are hunting the paralog FP class,
# so we must NOT let the production gates drop the very artefacts we study:
#   --min-qual 0                  keep low/zero-QUAL calls (refined paralog
#                                 QUAL can sit near 0)
#   --no-allele-balance-filter    keep the balanced (~0.5 VAF) paralog hets
#                                 the AB filter can't catch anyway, AND the
#                                 unbalanced point-mismap class, so the full
#                                 candidate set survives for analysis
# The DUST/complexity filter is left ON (default): low-complexity is a
# distinct, well-understood artefact class and would only add noise here.
#
# Usage:
#   benchmarks/tomato2/src/run_varcalling.sh
# Env: THREADS (default 16), REFERENCE, OUT (default results/cohort.vcf.gz)

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$TEST_DIR/../.." && pwd)"

REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"
BED="$TEST_DIR/regions_n160_200kb.bed"
PSP_DIR="$TEST_DIR/psp_files"
OUT_DIR="$TEST_DIR/results"
OUT="${OUT:-$OUT_DIR/cohort.vcf.gz}"
BIN="${POP_VAR_CALLER_BIN:-$PROJECT_ROOT/target/release/pop_var_caller}"
THREADS="${THREADS:-16}"

mkdir -p "$OUT_DIR"

shopt -s nullglob
psps=("$PSP_DIR"/*.psp)
shopt -u nullglob
(( ${#psps[@]} > 0 )) || { echo "no .psp in $PSP_DIR — run build_psp.sh first" >&2; exit 1; }

echo "=== var-calling: ${#psps[@]} samples, THREADS=$THREADS ==="
echo "reference: $REFERENCE"
echo "regions:   $BED"
echo "output:    $OUT"
echo "retention: --min-qual 0 --no-allele-balance-filter (DUST on)"
echo

time "$BIN" var-calling \
    --reference "$REFERENCE" --regions "$BED" \
    --output "$OUT" --threads "$THREADS" \
    --min-qual 0 --no-allele-balance-filter \
    "${psps[@]}"

echo
echo "done -> $OUT"
