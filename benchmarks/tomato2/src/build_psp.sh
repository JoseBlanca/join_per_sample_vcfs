#!/usr/bin/env bash
# Phase A (part 1): build one .psp per tomato2 CRAM, restricted to the
# random 200kb regions BED. Runs on the HOST (pileup needs no container).
#
# Mirrors benchmarks/tomato1/scripts/build_cohort_intermediates.sh but:
#   - input CRAMs live in benchmarks/tomato2/bams/
#   - output .psp go to benchmarks/tomato2/psp_files/
#   - analysis is restricted to regions_n160_200kb.bed via --regions
#
# Usage:
#   benchmarks/tomato2/src/build_psp.sh
# Env:
#   JOBS       parallel pileups (default 12)
#   MIN_MAPQ   pileup --min-mapq (default: caller default, unset)
#   REFERENCE  SL4.0 fasta

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$TEST_DIR/../.." && pwd)"

REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"
BED="$TEST_DIR/regions_n160_200kb.bed"
CRAM_DIR="$TEST_DIR/bams"
PSP_DIR="$TEST_DIR/psp_files"
SUFFIX=".bench.cram"
BIN="${POP_VAR_CALLER_BIN:-$PROJECT_ROOT/target/release/pop_var_caller}"
JOBS="${JOBS:-12}"
MIN_MAPQ="${MIN_MAPQ:-}"

mkdir -p "$PSP_DIR"

shopt -s nullglob
crams=("$CRAM_DIR"/*"$SUFFIX")
shopt -u nullglob
(( ${#crams[@]} > 0 )) || { echo "no *$SUFFIX in $CRAM_DIR" >&2; exit 1; }

echo "=== build psp: ${#crams[@]} sample(s), JOBS=$JOBS ==="
echo "reference: $REFERENCE"
echo "regions:   $BED ($(awk '{s+=$3-$2} END{printf "%.1f Mb / %d regions", s/1e6, NR}' "$BED"))"
echo "min-mapq:  ${MIN_MAPQ:-<default>}"
echo "output:    $PSP_DIR"
echo

build_one() {
    local cram="$1" base out
    base="$(basename "$cram" "$SUFFIX")"
    out="$PSP_DIR/$base.psp"
    if [[ -s "$out" ]]; then echo "[skip] $base"; return 0; fi
    local extra=()
    [[ -n "$MIN_MAPQ" ]] && extra=(--min-mapq "$MIN_MAPQ")
    if "$BIN" pileup --reference "$REFERENCE" --regions "$BED" \
            --output "$out" --threads 1 "${extra[@]}" "$cram" \
            > "$PSP_DIR/$base.build.log" 2>&1; then
        echo "[done] $base ($(du -h "$out" | cut -f1))"
    else
        echo "[FAIL] $base (see $PSP_DIR/$base.build.log)"; rm -f "$out"
    fi
}
export -f build_one
export BIN REFERENCE BED SUFFIX PSP_DIR MIN_MAPQ

printf '%s\n' "${crams[@]}" | xargs -P "$JOBS" -I{} bash -c 'build_one "$@"' _ {}
echo
echo "psp files now: $(ls "$PSP_DIR"/*.psp 2>/dev/null | wc -l | tr -d ' ') / ${#crams[@]}"
