#!/usr/bin/env bash
# Add @SQ M5 checksums to the per_sample BAM headers.
#
#   benchmarks/giab/src/add_m5_to_bams.sh [COVERAGE]
#
# Our pileup requires every @SQ line to carry an M5 checksum (it stores
# the reference md5 in the .psp header). The GIAB per_sample BAMs
# (HG003/HG004) were written without M5, so pileup rejects them with
#   "contig 'chr1' has no @SQ M5 checksum ...".
# (The HG002 300x CRAM already carries M5, so it is left untouched.)
#
# This computes M5 from the benchmark reference with `samtools dict`,
# merges it into each BAM header @SQ line by sequence name (preserving
# the original @SQ order and any existing tags such as AS), and writes a
# reheadered copy `<name>.m5.bam` (+ index) next to the original. The
# original BAMs are left untouched. `run_ours_per_sample.sh` prefers the
# `.m5.bam` copy when it exists.

set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SRC_DIR/.." && pwd)"

COVERAGE="${1:-300x}"
REFERENCE="$BENCH_DIR/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BAM_DIR="$BENCH_DIR/per_sample/bam/$COVERAGE"

command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "missing reference: $REFERENCE" >&2; exit 1; }

# SN -> M5 map from the reference, written once and reused per BAM.
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
M5_MAP="$TMP/m5_map.tsv"
samtools dict "$REFERENCE" 2>/dev/null \
  | awk -F'\t' '/^@SQ/{sn="";m5=""; for(i=1;i<=NF;i++){if($i~/^SN:/)sn=substr($i,4); if($i~/^M5:/)m5=substr($i,4)} if(sn!=""&&m5!="")print sn"\t"m5}' \
  > "$M5_MAP"
echo "reference M5 entries: $(wc -l < "$M5_MAP")"

shopt -s nullglob
bams=("$BAM_DIR"/*.sorted.bam)
shopt -u nullglob
if (( ${#bams[@]} == 0 )); then
    echo "no *.sorted.bam in $BAM_DIR (nothing to reheader)" >&2
    exit 0
fi

for bam in "${bams[@]}"; do
    [[ "$bam" == *.m5.bam ]] && continue
    out="${bam%.bam}.m5.bam"
    if [[ -s "$out" && "$out" -nt "$bam" ]]; then
        echo "[skip] $(basename "$out") already up to date"
        continue
    fi
    echo "[reheader] $(basename "$bam") -> $(basename "$out")"

    new_header="$TMP/header.sam"
    # Append M5:<checksum> to each @SQ line that lacks one, looked up by SN.
    samtools view -H "$bam" \
      | awk -F'\t' -v OFS='\t' -v map="$M5_MAP" '
          BEGIN { while ((getline line < map) > 0) { split(line, a, "\t"); m5[a[1]] = a[2] } }
          /^@SQ/ && $0 !~ /\tM5:/ {
              sn = ""
              for (i = 1; i <= NF; i++) if ($i ~ /^SN:/) sn = substr($i, 4)
              if (sn in m5) { $0 = $0 "\t" "M5:" m5[sn] }
              else { print "no M5 for contig " sn > "/dev/stderr"; exit 3 }
          }
          { print }
      ' > "$new_header"

    samtools reheader "$new_header" "$bam" > "$out"
    samtools index "$out"
done

echo "done."
