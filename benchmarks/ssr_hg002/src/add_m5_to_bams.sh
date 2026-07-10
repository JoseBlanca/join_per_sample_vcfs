#!/usr/bin/env bash
# add_m5_to_bams.sh — add @SQ M5 checksums to every coverage BAM header, IN PLACE.
#
# ssr-pileup requires each @SQ line to carry an M5 checksum (it stores the
# reference md5 in the .ssr.psp header). The GIAB HG002 novoalign BAM was written
# without M5, so pileup rejects it with "contig 'chr1' has no md5 ...".
#
# The M5 map is computed from the benchmark reference with `samtools dict`, merged
# into each BAM's @SQ lines by sequence name, and the BAM is reheadered in place
# (temp + mv; reads untouched) and re-indexed. Idempotent: a BAM that already has
# M5 on every @SQ is skipped.
#
# Run from ssr_hg002/:  bash src/add_m5_to_bams.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REFERENCE="$ROOT/../giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
[[ -f "$REFERENCE" ]] || { echo "missing reference: $REFERENCE" >&2; exit 1; }

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
M5_MAP="$TMP/m5_map.tsv"
samtools dict "$REFERENCE" 2>/dev/null \
  | awk -F'\t' '/^@SQ/{sn="";m5=""; for(i=1;i<=NF;i++){if($i~/^SN:/)sn=substr($i,4); if($i~/^M5:/)m5=substr($i,4)} if(sn!=""&&m5!="")print sn"\t"m5}' \
  > "$M5_MAP"
echo "reference M5 entries: $(wc -l < "$M5_MAP")"

shopt -s nullglob
bams=("$ROOT"/bam/*x/*.bam)
shopt -u nullglob

for bam in "${bams[@]}"; do
  if ! samtools view -H "$bam" | awk '/^@SQ/ && $0 !~ /\tM5:/{f=1} END{exit !f}'; then
    echo "[skip] $(basename "$bam") already has M5 on every @SQ"
    continue
  fi
  echo "[reheader] $(basename "$bam")"
  new_header="$TMP/header.sam"
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
  samtools reheader "$new_header" "$bam" > "$bam.m5.tmp"
  mv "$bam.m5.tmp" "$bam"
  samtools index "$bam"
done
echo "done."
