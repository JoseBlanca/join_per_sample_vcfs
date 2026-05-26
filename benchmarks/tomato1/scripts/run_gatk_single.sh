#!/usr/bin/env bash
# Single-sample GATK test: one CRAM -> single-sample VCF via
# HaplotypeCaller in default (non-GVCF) mode, restricted to the
# benchmark BED.
#
# Sample: SRR17274057.p1.bench.cram
# Output: tmp/tomato_cohort_test/results/gatk/single_SRR17274057.vcf
#
# Env overrides:
#   GATK_BIN        binary to invoke (default: /opt/gatk/gatk)
#   REFERENCE       SL4.0 fasta (.fai + .dict siblings required)
#   THREADS         --native-pair-hmm-threads (default: 4 — host
#                   P-core budget, memory: feedback_host_p_core_budget)
#   JAVA_HEAP       -Xmx (default: 4g)
#   EXTRA_ARGS      appended verbatim to the HaplotypeCaller line
#
# Reference requirements vs the other callers:
# - GATK uniquely requires a `.dict` next to the FASTA (Picard sequence
#   dictionary). The tomato SL4.0 reference at the default path already
#   has one — see ../regions.bed for what we tested.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

SAMPLE="SRR17274057"
CRAM="$TEST_DIR/crams/${SAMPLE}.p1.bench.cram"
BED="$TEST_DIR/regions.bed"
REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"
THREADS="${THREADS:-4}"
JAVA_HEAP="${JAVA_HEAP:-4g}"
EXTRA_ARGS=${EXTRA_ARGS:-}

GATK_BIN="${GATK_BIN:-/opt/gatk/gatk}"
[[ -x "$GATK_BIN" ]] || { echo "GATK binary not executable: $GATK_BIN" >&2; exit 1; }

# --- preflight ---
REF_DICT="${REFERENCE%.fa}.dict"   # convention: foo.fa -> foo.dict
for f in "$CRAM" "${CRAM}.crai" "$BED" "$REFERENCE" "${REFERENCE}.fai" "$REF_DICT"; do
    [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done

# --- output ---
OUT_DIR="$TEST_DIR/results/gatk"
OUT_VCF="$OUT_DIR/single_${SAMPLE}.vcf"
LOG="$OUT_DIR/single_${SAMPLE}.log"
mkdir -p "$OUT_DIR"

echo "binary    : $GATK_BIN ($($GATK_BIN --version 2>&1 | head -1))"
echo "sample    : $SAMPLE"
echo "input     : $CRAM"
echo "reference : $REFERENCE"
echo "regions   : $BED ($(wc -l < "$BED") intervals)"
echo "threads   : $THREADS (--native-pair-hmm-threads)"
echo "heap      : $JAVA_HEAP"
echo "output    : $OUT_VCF"
echo "log       : $LOG"
echo

t0=$(date +%s)
"$GATK_BIN" --java-options "-Xmx${JAVA_HEAP}" HaplotypeCaller \
    --reference "$REFERENCE" \
    --input "$CRAM" \
    --intervals "$BED" \
    --output "$OUT_VCF" \
    --native-pair-hmm-threads "$THREADS" \
    $EXTRA_ARGS \
    2> >(tee "$LOG" >&2)
t1=$(date +%s)

echo
echo "elapsed: $((t1 - t0)) s"
echo "records: $(grep -vc '^#' "$OUT_VCF")"
