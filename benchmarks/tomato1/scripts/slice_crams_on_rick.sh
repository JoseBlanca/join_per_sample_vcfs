#!/usr/bin/env bash
# Slice every *.cram in IN_DIR down to the intervals in BED, writing
# region-restricted CRAMs (+ .crai indexes) into OUT_DIR.
#
# Run this ON RICK (where the full CRAMs live). Requires samtools and the
# same reference fasta the CRAMs were aligned against.
#
# Usage:
#   ./slice_crams_on_rick.sh REF BED IN_DIR OUT_DIR
#
# Example:
#   ./slice_crams_on_rick.sh \
#       /home/joxi/refs/S_lycopersicum_chromosomes.4.00.fa \
#       /home/joxi/tmp/benchmark_cohort/regions.bed \
#       /media/tomato25_bams/crams/PRJNA790656 \
#       /home/joxi/tmp/benchmark_cohort/crams
#
# Tunables: THREADS env var (default 4).

set -euo pipefail

REF=${1:?"usage: $0 REF BED IN_DIR OUT_DIR"}
BED=${2:?"usage: $0 REF BED IN_DIR OUT_DIR"}
IN_DIR=${3:?"usage: $0 REF BED IN_DIR OUT_DIR"}
OUT_DIR=${4:?"usage: $0 REF BED IN_DIR OUT_DIR"}
THREADS=${THREADS:-4}

for f in "$REF" "$REF.fai" "$BED"; do
    [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done
[[ -d "$IN_DIR" ]] || { echo "not a directory: $IN_DIR" >&2; exit 1; }

mkdir -p "$OUT_DIR"

shopt -s nullglob
crams=("$IN_DIR"/*.cram)
if (( ${#crams[@]} == 0 )); then
    echo "no *.cram in $IN_DIR" >&2
    exit 1
fi

echo "slicing ${#crams[@]} CRAMs -> $OUT_DIR (threads=$THREADS)"
echo "regions: $(wc -l < "$BED") intervals from $BED"
echo

for cram in "${crams[@]}"; do
    base=$(basename "$cram" .cram)
    out="$OUT_DIR/${base}.bench.cram"
    if [[ -s "$out" && -s "${out}.crai" ]]; then
        echo "[skip] $base (already sliced)"
        continue
    fi
    echo "[slice] $base"
    # -C: CRAM output; -T: reference; -L: BED restriction; --write-index: build .crai
    samtools view \
        -C \
        -T "$REF" \
        -L "$BED" \
        -@ "$THREADS" \
        --write-index \
        -o "$out" \
        "$cram"
done

echo
echo "done. output sizes:"
du -h "$OUT_DIR"/*.bench.cram 2>/dev/null | sort -h
