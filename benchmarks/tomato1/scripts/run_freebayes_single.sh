#!/usr/bin/env bash
# Single-sample freebayes test: one CRAM -> single-sample VCF, restricted
# to the benchmark BED.
#
# Sample: SRR17274057.p1.bench.cram
# Output: tmp/tomato_cohort_test/results/freebayes/single_SRR17274057.vcf
#
# Env overrides:
#   FREEBAYES_BIN   binary to invoke (default: auto-detect `freebayes` on PATH)
#   REFERENCE       SL4.0 fasta (.fai sibling required)
#   PLOIDY          --ploidy (default: 2)
#   MIN_QUAL        drop records with QUAL < this before writing
#                   (default: 30 — freebayes emits a long low-QUAL
#                   tail that the other callers don't, and leaving
#                   it in skews any cross-caller comparison).
#   EXTRA_ARGS      appended verbatim to the freebayes command line
#
# Notes:
# - Freebayes is single-threaded. Wall-time is dominated by the CRAM
#   reader + the per-locus model fit; for the 2 Mb / 1-sample test set
#   this is seconds-to-low-minutes, not worth the freebayes-parallel
#   wrapping. The cohort script handles 18-sample scale.
# - We pass the BED via `-t` so calls outside the slice regions
#   (possible from soft-clipped read overhangs) are excluded from the
#   output up front, matching what we'll do for GATK.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

SAMPLE="SRR17274057"
CRAM="$TEST_DIR/crams/${SAMPLE}.p1.bench.cram"
BED="$TEST_DIR/regions.bed"
REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"
PLOIDY="${PLOIDY:-2}"
MIN_QUAL="${MIN_QUAL:-30}"
EXTRA_ARGS=${EXTRA_ARGS:-}

# --- binary discovery ---
if [[ -z "${FREEBAYES_BIN:-}" ]]; then
    if command -v freebayes >/dev/null; then
        FREEBAYES_BIN="$(command -v freebayes)"
    fi
fi
if [[ -z "${FREEBAYES_BIN:-}" || ! -x "${FREEBAYES_BIN}" ]]; then
    echo "no freebayes binary found." >&2
    echo "install via package manager (apt/brew/conda) or set FREEBAYES_BIN=<path>" >&2
    exit 1
fi

# --- preflight ---
for f in "$CRAM" "${CRAM}.crai" "$BED" "$REFERENCE" "${REFERENCE}.fai"; do
    [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done

# --- output ---
OUT_DIR="$TEST_DIR/results/freebayes"
OUT_VCF="$OUT_DIR/single_${SAMPLE}.vcf"
LOG="$OUT_DIR/single_${SAMPLE}.log"
mkdir -p "$OUT_DIR"

echo "binary    : $FREEBAYES_BIN ($($FREEBAYES_BIN --version 2>&1 | head -1))"
echo "sample    : $SAMPLE"
echo "input     : $CRAM"
echo "reference : $REFERENCE"
echo "regions   : $BED ($(wc -l < "$BED") intervals)"
echo "ploidy    : $PLOIDY"
echo "min QUAL  : $MIN_QUAL (records below are dropped pre-write)"
echo "output    : $OUT_VCF"
echo "log       : $LOG"
echo

# Pre-filter QUAL inline. awk avoids a bcftools dependency (not always
# installed on the host); header lines pass through, data rows are
# kept iff column 6 is numeric and >= MIN_QUAL. `$6 + 0` coerces the
# QUAL string to a number; non-numeric values fall to 0 and get
# dropped, which is the right behaviour for a fairness cap.
FILTER_AWK='/^#/ { print; next } $6 != "." && ($6 + 0) >= min { print }'

t0=$(date +%s)
"$FREEBAYES_BIN" \
    -f "$REFERENCE" \
    -t "$BED" \
    -p "$PLOIDY" \
    $EXTRA_ARGS \
    "$CRAM" \
    2> >(tee "$LOG" >&2) \
  | awk -v min="$MIN_QUAL" "$FILTER_AWK" > "$OUT_VCF"
t1=$(date +%s)

echo
echo "elapsed: $((t1 - t0)) s"
echo "records: $(grep -vc '^#' "$OUT_VCF") (QUAL >= $MIN_QUAL, no other filters)"
