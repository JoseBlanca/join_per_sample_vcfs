#!/usr/bin/env bash
# Cohort freebayes test: all 18 CRAMs -> joint VCF (single invocation).
#
# Unlike pop_var_caller, freebayes has no per-sample intermediate
# (no GVCF-style two-stage). All CRAMs are passed in one call and the
# sample names come from each file's @RG SM tag.
#
# Output: tmp/tomato_cohort_test/results/freebayes/cohort.vcf
#
# Env overrides:
#   FREEBAYES_BIN   binary to invoke (default: auto-detect)
#   REFERENCE       SL4.0 fasta (.fai sibling required)
#   PLOIDY          --ploidy (default: 2)
#   MIN_QUAL        drop records with QUAL < this before writing
#                   (default: 30 — freebayes emits a long low-QUAL
#                   tail that the other callers don't, and leaving
#                   it in skews any cross-caller comparison).
#   EXTRA_ARGS      appended verbatim to the freebayes command line
#   PARALLEL        if "1" and `freebayes-parallel` + GNU `parallel`
#                   are on PATH, fan out per-BED-region. Default 0
#                   (single-thread; matches the simplest comparison
#                   surface vs our caller). For wall-time-sensitive
#                   runs set PARALLEL=1 THREADS=4.
#   THREADS         worker count when PARALLEL=1 (default 4 — host
#                   P-core budget, memory: feedback_host_p_core_budget).
#
# Cohort scale notes:
# - 18 samples × 2 Mb. Single-threaded freebayes is the slowest part
#   of the three-way comparison (expect minutes). The variants are
#   per-locus serial; turning PARALLEL on roughly divides by THREADS.
# - freebayes deduces per-sample id from each input's @RG SM tag.
#   The bench CRAMs were header-fixed by fix_cram_headers.sh which
#   preserved one SM per file, so each input contributes one sample.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

CRAM_DIR="$TEST_DIR/crams"
BED="$TEST_DIR/regions.bed"
REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"
PLOIDY="${PLOIDY:-2}"
MIN_QUAL="${MIN_QUAL:-30}"
EXTRA_ARGS=${EXTRA_ARGS:-}
PARALLEL="${PARALLEL:-0}"
THREADS="${THREADS:-4}"

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
for f in "$BED" "$REFERENCE" "${REFERENCE}.fai"; do
    [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done
shopt -s nullglob
crams=("$CRAM_DIR"/*.bench.cram)
if (( ${#crams[@]} == 0 )); then
    echo "no *.bench.cram in $CRAM_DIR" >&2
    exit 1
fi
for c in "${crams[@]}"; do
    [[ -f "${c}.crai" ]] || { echo "missing index: ${c}.crai" >&2; exit 1; }
done

# --- output ---
OUT_DIR="$TEST_DIR/results/freebayes"
OUT_VCF="$OUT_DIR/cohort.vcf"
LOG="$OUT_DIR/cohort.log"
mkdir -p "$OUT_DIR"

echo "binary    : $FREEBAYES_BIN ($($FREEBAYES_BIN --version 2>&1 | head -1))"
echo "samples   : ${#crams[@]}"
echo "reference : $REFERENCE"
echo "regions   : $BED ($(wc -l < "$BED") intervals)"
echo "ploidy    : $PLOIDY"
echo "min QUAL  : $MIN_QUAL (records below are dropped pre-write)"
echo "mode      : $([[ "$PARALLEL" == "1" ]] && echo "parallel (threads=$THREADS)" || echo "single-threaded")"
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
if [[ "$PARALLEL" == "1" ]]; then
    command -v freebayes-parallel >/dev/null || {
        echo "PARALLEL=1 requested but freebayes-parallel not on PATH" >&2
        exit 1
    }
    command -v parallel >/dev/null || {
        echo "PARALLEL=1 requested but GNU parallel not on PATH" >&2
        exit 1
    }
    # freebayes-parallel takes a regions file (one region per line in
    # `chr:start-end` form, NOT BED) plus a thread count plus all the
    # usual freebayes args. Translate the BED on the fly.
    REGIONS_TXT="$OUT_DIR/regions.fbparallel.txt"
    awk '{ printf "%s:%d-%d\n", $1, $2+1, $3 }' "$BED" > "$REGIONS_TXT"
    freebayes-parallel \
        "$REGIONS_TXT" \
        "$THREADS" \
        -f "$REFERENCE" \
        -p "$PLOIDY" \
        $EXTRA_ARGS \
        "${crams[@]}" \
        2> >(tee "$LOG" >&2) \
      | awk -v min="$MIN_QUAL" "$FILTER_AWK" > "$OUT_VCF"
else
    "$FREEBAYES_BIN" \
        -f "$REFERENCE" \
        -t "$BED" \
        -p "$PLOIDY" \
        $EXTRA_ARGS \
        "${crams[@]}" \
        2> >(tee "$LOG" >&2) \
      | awk -v min="$MIN_QUAL" "$FILTER_AWK" > "$OUT_VCF"
fi
t1=$(date +%s)

echo
echo "elapsed: $((t1 - t0)) s"
echo "records: $(grep -vc '^#' "$OUT_VCF") (QUAL >= $MIN_QUAL, no other filters)"
echo "samples in vcf: $(grep -m1 '^#CHROM' "$OUT_VCF" | awk '{print NF-9}')"
