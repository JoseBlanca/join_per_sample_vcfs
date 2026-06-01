#!/usr/bin/env bash
# Build the per-sample cohort intermediates that the §4/§5 perf scripts
# read (they assume these already exist):
#
#   psp   -> results/ours/cohort/psp/<sample>.psp        (pop_var_caller pileup; HOST)
#   gvcf  -> results/gatk/cohort/gvcf/<sample>.g.vcf.gz  (GATK HC -ERC GVCF; CONTAINER)
#
# Processes the first COHORT_N CRAMs (alphabetical — same subset the perf
# scripts' pick_subset selects), up to JOBS in parallel, skipping any
# intermediate that already exists. Re-run after adding samples / bumping
# COHORT_N to fill in the rest.
#
# Usage:
#   # PSPs on the host:
#   benchmarks/tomato1/scripts/build_cohort_intermediates.sh psp
#   # GVCFs in the container (GATK + reference live there):
#   DEV_EXTRA_MOUNT=$HOME/genomes ./scripts/dev.sh \
#       benchmarks/tomato1/scripts/build_cohort_intermediates.sh gvcf
#
# Env:
#   COHORT_N            number of samples to build (default 50)
#   JOBS                parallel builders (default 4)
#   REFERENCE           SL4.0 fasta (.fai sibling; .dict too for gvcf)
#   POP_VAR_CALLER_BIN  pop_var_caller binary (psp mode)
#   GATK_BIN            GATK binary (gvcf mode; default /opt/gatk/gatk)
#   JAVA_HEAP           -Xmx per parallel HaplotypeCaller (gvcf; default 2g)

set -uo pipefail

MODE="${1:-}"
case "$MODE" in
    psp|gvcf) ;;
    *) echo "usage: $0 [psp|gvcf]" >&2; exit 2 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$TEST_DIR/../.." && pwd)"

COHORT_N="${COHORT_N:-50}"
JOBS="${JOBS:-4}"
REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"
SUFFIX=".bench.cram"
CRAM_DIR="$TEST_DIR/crams"
BED="$TEST_DIR/regions.bed"
PSP_DIR="$TEST_DIR/results/ours/cohort/psp"
GVCF_DIR="$TEST_DIR/results/gatk/cohort/gvcf"
POP_VAR_CALLER_BIN="${POP_VAR_CALLER_BIN:-$PROJECT_ROOT/target/release/pop_var_caller}"
GATK_BIN="${GATK_BIN:-/opt/gatk/gatk}"
JAVA_HEAP="${JAVA_HEAP:-2g}"

# First COHORT_N CRAMs, sorted (matches perf pick_subset).
shopt -s nullglob
all_crams=("$CRAM_DIR"/*"$SUFFIX")
shopt -u nullglob
(( ${#all_crams[@]} > 0 )) || { echo "no *$SUFFIX in $CRAM_DIR" >&2; exit 1; }
IFS=$'\n' all_crams=($(printf '%s\n' "${all_crams[@]}" | sort)); unset IFS
n=$(( COHORT_N < ${#all_crams[@]} ? COHORT_N : ${#all_crams[@]} ))
crams=("${all_crams[@]:0:n}")

echo "=== build $MODE: ${#crams[@]} sample(s), JOBS=$JOBS ==="
echo "reference: $REFERENCE"
echo

build_one_psp() {
    local cram="$1" base out
    base="$(basename "$cram" "$SUFFIX")"
    out="$PSP_DIR/$base.psp"
    if [[ -s "$out" ]]; then echo "[skip psp] $base"; return 0; fi
    echo "[psp] $base"
    "$POP_VAR_CALLER_BIN" pileup \
        --reference "$REFERENCE" --output "$out" --threads 1 \
        "$cram" > "$PSP_DIR/$base.build.log" 2>&1 \
        && echo "[done psp] $base" || echo "[FAIL psp] $base (see $base.build.log)"
}

build_one_gvcf() {
    local cram="$1" base out
    base="$(basename "$cram" "$SUFFIX")"
    out="$GVCF_DIR/$base.g.vcf.gz"
    if [[ -s "$out" && -s "$out.tbi" ]]; then echo "[skip gvcf] $base"; return 0; fi
    echo "[gvcf] $base"
    "$GATK_BIN" --java-options "-Xmx${JAVA_HEAP}" HaplotypeCaller \
        --reference "$REFERENCE" --input "$cram" --intervals "$BED" \
        --output "$out" --emit-ref-confidence GVCF \
        --native-pair-hmm-threads 1 > "$GVCF_DIR/$base.build.log" 2>&1 \
        && echo "[done gvcf] $base" || echo "[FAIL gvcf] $base (see $base.build.log)"
}

export -f build_one_psp build_one_gvcf
export POP_VAR_CALLER_BIN GATK_BIN REFERENCE BED SUFFIX PSP_DIR GVCF_DIR JAVA_HEAP

if [[ "$MODE" == psp ]]; then
    mkdir -p "$PSP_DIR"
    printf '%s\n' "${crams[@]}" | xargs -P "$JOBS" -I{} bash -c 'build_one_psp "$@"' _ {}
    echo; echo "psp files now: $(ls "$PSP_DIR"/*.psp 2>/dev/null | wc -l | tr -d ' ')"
else
    mkdir -p "$GVCF_DIR"
    printf '%s\n' "${crams[@]}" | xargs -P "$JOBS" -I{} bash -c 'build_one_gvcf "$@"' _ {}
    echo; echo "gvcf files now: $(ls "$GVCF_DIR"/*.g.vcf.gz 2>/dev/null | wc -l | tr -d ' ')"
fi
