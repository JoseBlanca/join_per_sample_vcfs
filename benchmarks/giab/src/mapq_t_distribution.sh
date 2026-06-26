#!/usr/bin/env bash
# Dump the MAPQ-difference Welch's-t value the filter computes for every
# candidate record, labelled TP / FP against the GIAB truth, so we can see the
# distribution of t for real variants vs artefacts and judge where the
# --min-mapq-diff-t threshold should sit.
#
#   benchmarks/giab/src/mapq_t_distribution.sh [VARIANT] [COVERAGE]
#
# Method: run var-calling with the filter DISABLED (--min-mapq-diff-t=-inf) so
# every candidate record is emitted AND reaches the (instrumented) filter,
# which prints `MAPQDIFF\t<chrom_id>\t<pos>\t<alt_idx>\t<n_ref>\t<mean_ref>\t
# <n_alt>\t<mean_alt>\t<t>\t<thr>\t<decision>` to stderr (requires the binary
# built with the PVC_DEBUG_MAPQ_DIFF instrumentation). We map chrom_id->name
# via the reference .fai, classify each record SNP/indel from the emitted VCF,
# and label TP/FP by bcftools isec vs the truth. Output:
#   results/per_sample/<cov>/mapq_t_dist.tsv
#     sample  chrom  pos  class  n_ref  mean_ref  n_alt  mean_alt  t  status

set -euo pipefail
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SRC_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$BENCH_DIR/../.." && pwd)"

VARIANT="${1:-ours_nobaq_nodust}"
COVERAGE="${2:-300x}"
BIN="${BIN:-$PROJECT_ROOT/target/release/pop_var_caller}"

REFERENCE="$BENCH_DIR/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BED_DIR="$BENCH_DIR/per_sample/bed"
TRUTH_DIR="$BENCH_DIR/per_sample/vcf"
RES_DIR="$BENCH_DIR/results/per_sample/$COVERAGE/$VARIANT"
WORK="$RES_DIR/mapq_tdist"
TSV="$RES_DIR/mapq_t_dist.tsv"
SAMPLES=(HG002 HG003 HG004)

command -v bcftools >/dev/null || { echo "bcftools not on PATH" >&2; exit 1; }
[[ -x "$BIN" ]] || { echo "binary not found: $BIN" >&2; exit 1; }
mkdir -p "$WORK"

# chrom_id -> name from the .fai (line order == chrom_id, 0-based).
CHROM_MAP="$WORK/chrom_map.tsv"
awk 'BEGIN{OFS="\t"} {print NR-1, $1}' "${REFERENCE}.fai" > "$CHROM_MAP"

printf 'sample\tchrom\tpos\tclass\tn_ref\tmean_ref\tn_alt\tmean_alt\tt\tstatus\n' > "$TSV"

norm() { # src class out passflag bed
    bcftools view $4 -T "$5" "$1" -Ou \
      | bcftools norm -f "$REFERENCE" -m -any -Ou 2>/dev/null \
      | bcftools view -v "$2" -Oz -o "$3"
    bcftools index -tf "$3"
}

for s in "${SAMPLES[@]}"; do
    psp="$RES_DIR/${s}.psp"
    bed="$BED_DIR/${s}_bench_azar_merged_100.bed"
    truth="$TRUTH_DIR/${s}_GRCh38_1_22_v4.2.1_benchmark.selected_100.vcf.gz"
    [[ -s "$psp" ]] || { echo "missing psp: $psp" >&2; continue; }

    vcf="$WORK/${s}.inf.vcf"
    dump="$WORK/${s}.mapqdiff.tsv"
    raw="$WORK/${s}.stderr.log"
    # Use the finite minimum (-50), NOT -inf: the instrumentation lives after
    # the filter's `!threshold.is_finite()` early-return, so -inf prints nothing.
    # The computed t is independent of the threshold, and -50 keeps essentially
    # every record (so the emitted VCF still labels TP/FP for nearly all rows).
    echo "[$s] var-calling (filter at -50, instrumented) + capturing t dump"
    PVC_DEBUG_MAPQ_DIFF=1 "$BIN" var-calling --reference "$REFERENCE" --regions "$bed" \
        --no-complexity-filter --min-mapq-diff-t=-50 \
        --output "$vcf" --threads 1 "$psp" 2> "$raw" || true
    grep '^MAPQDIFF' "$raw" > "$dump" || true

    # Per-class TP/FP position table (chrom:pos<TAB>status<TAB>class) from isec
    # vs truth. Built straight from bcftools output — no bash assoc arrays
    # (this script must run under macOS bash 3.2).
    STATUS_FILE="$WORK/${s}.status.tsv"
    : > "$STATUS_FILE"
    td="$(mktemp -d)"
    for cls in snps indels; do
        norm "$truth" "$cls" "$td/t.vcf.gz" "-f PASS" "$bed"
        norm "$vcf"   "$cls" "$td/q.vcf.gz" ""       "$bed"
        bcftools isec -p "$td/i" "$td/t.vcf.gz" "$td/q.vcf.gz" 2>/dev/null
        # 0003 = shared, query side (TP); 0001 = query-only (FP)
        bcftools query -f '%CHROM:%POS\tTP\t'"$cls"'\n' "$td/i/0003.vcf" 2>/dev/null >> "$STATUS_FILE"
        bcftools query -f '%CHROM:%POS\tFP\t'"$cls"'\n' "$td/i/0001.vcf" 2>/dev/null >> "$STATUS_FILE"
        rm -rf "$td/i"
    done
    rm -rf "$td"

    # Join the t dump to chrom names + status. The dump pos is the record's
    # 1-based locus start (== VCF POS for SNPs; for left-aligned indels norm
    # may shift POS, so indel status matching is best-effort — SNPs are exact).
    awk -F'\t' -v OFS='\t' -v s="$s" -v cmap="$CHROM_MAP" -v sf="$STATUS_FILE" '
        BEGIN {
            while ((getline line < cmap) > 0) { split(line, a, "\t"); name[a[1]] = a[2] }
            while ((getline line < sf) > 0) { split(line, b, "\t"); st[b[1]] = b[2]; cl[b[1]] = b[3] }
        }
        $1=="MAPQDIFF" {
            chrom = name[$2]; pos = $3; key = chrom ":" pos
            status = (key in st) ? st[key] : "other"
            class  = (key in cl) ? cl[key] : "na"
            # cols: chrom_id pos alt_idx n_ref mean_ref n_alt mean_alt t thr decision
            print s, chrom, pos, class, $5, $6, $7, $8, $9, status
        }
    ' "$dump" >> "$TSV"
done

echo "wrote $TSV"
echo "rows: $(($(wc -l < "$TSV") - 1))"
echo "=== status x class counts ==="
awk -F'\t' 'NR>1{c[$10"\t"$4]++} END{for(k in c) print c[k]"\t"k}' "$TSV" | sort
