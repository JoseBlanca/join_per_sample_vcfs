#!/usr/bin/env bash
# Sort the per_sample benchmark BEDs into REFERENCE CONTIG ORDER (in place).
#
#   benchmarks/giab/src/sort_beds.sh
#
# Why this exists: the per_sample regions were selected at random and the
# resulting BEDs are in arbitrary chromosome order (e.g. HG003 starts chr7,
# chr6, chr12, ...). Our caller emits variants in reference order regardless,
# but freebayes (run with --targets BED) emits records in BED order, producing
# a VCF that is NOT coordinate-sorted -> `bcftools index -t` then fails with
# "Chromosome blocks not continuous". Sorting the BEDs once, in the same contig
# order as the reference .fai (numeric: chr2 before chr10, not lexicographic),
# makes every downstream callset coordinate-sorted and indexable.
#
# Sort key: (reference contig index from the .fai, then numeric start, then
# numeric end). Idempotent — re-running on an already-sorted BED is a no-op.
#
# The BEDs live under the gitignored data tree; this script is the committed,
# reproducible record of the normalization.

set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SRC_DIR/.." && pwd)"

FAI="$BENCH_DIR/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai"
BED_DIR="$BENCH_DIR/per_sample/bed"

[[ -f "$FAI" ]] || { echo "missing reference .fai: $FAI" >&2; exit 1; }
[[ -d "$BED_DIR" ]] || { echo "missing BED dir: $BED_DIR" >&2; exit 1; }

sort_one() {
    local bed="$1" tmp
    tmp="$(mktemp "${bed}.sorted.XXXXXX")"
    # Decorate each BED line with its contig's .fai index, sort by
    # (index, start, end), then undecorate. Rows whose contig is absent from
    # the .fai sort last (index defaults to a large sentinel) so nothing is
    # silently dropped.
    awk -v OFS='\t' '
        NR == FNR { idx[$1] = FNR; next }
        {
            i = ($1 in idx) ? idx[$1] : 1e9
            print i, $2, $3, $0
        }
    ' "$FAI" "$bed" \
        | sort -k1,1n -k2,2n -k3,3n \
        | cut -f4- > "$tmp"
    # Sanity: line count must be preserved.
    if [[ "$(wc -l < "$bed")" != "$(wc -l < "$tmp")" ]]; then
        echo "line count changed for $bed — aborting, leaving original" >&2
        rm -f "$tmp"; return 1
    fi
    mv "$tmp" "$bed"
}

shopt -s nullglob
beds=("$BED_DIR"/*_bench_azar_merged_100.bed)
if (( ${#beds[@]} == 0 )); then
    echo "no per_sample BEDs found under $BED_DIR" >&2
    exit 1
fi

for bed in "${beds[@]}"; do
    before="$(cut -f1 "$bed" | head -1)"
    sort_one "$bed"
    after="$(cut -f1 "$bed" | head -1)"
    printf 'sorted %-60s (first contig %s -> %s, %s intervals)\n' \
        "$(basename "$bed")" "$before" "$after" "$(wc -l < "$bed" | tr -d ' ')"
done

echo "done. BEDs are now in reference contig order."
