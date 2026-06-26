#!/usr/bin/env bash
# Build a per-SNP-call signal table for the 2D TP/FP separation analysis:
# combine the MAPQ-difference signal with depth and allele balance (VAF), each
# labelled TP/FP against the GIAB truth. Feeds mapq_depth_dashboard.py.
#
#   benchmarks/giab/src/signal_table.sh [VARIANT] [COVERAGE]
#
# Everything we need is already in the var-calling VCF INFO/FORMAT — no extra
# instrumentation:
#   INFO/DP      total locus depth
#   INFO/MQRef   mean MAPQ of REF-supporting reads
#   INFO/MQAlt   mean MAPQ of ALT-supporting reads (per ALT)
#   INFO/MQDiff  MQAlt - MQRef  (depth-INDEPENDENT effect size; <0 = alt lower)
#   INFO/MQDiffT Welch's t of the MAPQ difference (depth-SENSITIVE statistic)
#   INFO/AF      EM allele frequency
#   FORMAT/AD    REF,ALT read depths -> VAF = alt/(ref+alt)
#
# We use the variant's per-sample VCFs (e.g. the high-recall preset, which keeps
# every candidate so FPs are present), restrict to SNPs within the sample BED,
# split multiallelics (so INFO per-ALT fields line up), isec vs the truth to
# label TP (0003, query side) / FP (0001), and query the fields. Output:
#   results/per_sample/<cov>/<variant>/signal_table.tsv
#     sample chrom pos status dp mqref mqalt mqdiff mqdifft af ref_ad alt_ad vaf

set -euo pipefail
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SRC_DIR/.." && pwd)"

VARIANT="${1:-high-recall}"
COVERAGE="${2:-300x}"
REFERENCE="$BENCH_DIR/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BED_DIR="$BENCH_DIR/per_sample/bed"
TRUTH_DIR="$BENCH_DIR/per_sample/vcf"
RES_DIR="$BENCH_DIR/results/per_sample/$COVERAGE/$VARIANT"
TSV="$RES_DIR/signal_table.tsv"
SAMPLES=(HG002 HG003 HG004)

command -v bcftools >/dev/null || { echo "bcftools not on PATH" >&2; exit 1; }
[[ -d "$RES_DIR" ]] || { echo "no result dir: $RES_DIR (run the preset first)" >&2; exit 1; }

printf 'sample\tchrom\tpos\tstatus\tgt\tdp\tmqref\tmqalt\tmqdiff\tmqdifft\taf\tref_ad\talt_ad\tvaf\n' > "$TSV"

# Query the per-ALT signal fields from a (normalized, biallelic) VCF, tagging
# each row with the given status. [%AD] expands the single sample's REF,ALT
# depths; awk splits it to compute VAF. Missing fields come through as ".".
emit_rows() { # vcf status sample
    local vcf="$1" status="$2" s="$3"
    bcftools query -f '%CHROM\t%POS\t[%GT]\t%INFO/DP\t%INFO/MQRef\t%INFO/MQAlt\t%INFO/MQDiff\t%INFO/MQDiffT\t%INFO/AF\t[%AD]\n' "$vcf" 2>/dev/null \
      | awk -F'\t' -v OFS='\t' -v s="$s" -v st="$status" '{
            split($10, ad, ",")              # AD = ref,alt
            ref_ad = (ad[1]=="."?"":ad[1]) + 0
            alt_ad = (ad[2]=="."?"":ad[2]) + 0
            tot = ref_ad + alt_ad
            vaf = (tot>0) ? alt_ad/tot : ""
            # cols: CHROM POS GT DP MQRef MQAlt MQDiff MQDiffT AF AD
            print s, $1, $2, st, $3, $4, $5, $6, $7, $8, $9, ref_ad, alt_ad, vaf
        }'
}

for s in "${SAMPLES[@]}"; do
    vcf="$RES_DIR/${s}.vcf"
    bed="$BED_DIR/${s}_bench_azar_merged_100.bed"
    truth="$TRUTH_DIR/${s}_GRCh38_1_22_v4.2.1_benchmark.selected_100.vcf.gz"
    [[ -s "$vcf" ]] || { echo "missing VCF: $vcf" >&2; continue; }
    echo "[$s] labelling SNP calls TP/FP and extracting signals"

    td="$(mktemp -d)"
    # Normalize (BED-restrict, biallelic-split, SNPs only). Keep INFO/FORMAT.
    bcftools view -T "$bed" "$vcf" -Ou \
      | bcftools norm -f "$REFERENCE" -m -any -Ou 2>/dev/null \
      | bcftools view -v snps -Oz -o "$td/q.vcf.gz"
    bcftools index -t "$td/q.vcf.gz"
    bcftools view -f PASS -T "$bed" "$truth" -Ou \
      | bcftools norm -f "$REFERENCE" -m -any -Ou 2>/dev/null \
      | bcftools view -v snps -Oz -o "$td/t.vcf.gz"
    bcftools index -t "$td/t.vcf.gz"

    bcftools isec -p "$td/i" "$td/t.vcf.gz" "$td/q.vcf.gz" 2>/dev/null
    # 0003 = shared, query side (TP);  0001 = query-only (FP)
    emit_rows "$td/i/0003.vcf" TP "$s" >> "$TSV"
    emit_rows "$td/i/0001.vcf" FP "$s" >> "$TSV"
    rm -rf "$td"
done

echo "wrote $TSV"
echo "rows: $(($(wc -l < "$TSV") - 1))"
echo "=== status counts ==="
awk -F'\t' 'NR>1{c[$4]++} END{for(k in c) print k, c[k]}' "$TSV"
echo "=== mean VAF / depth / MQDiff by status (sanity) ==="
# cols: 1 sample 2 chrom 3 pos 4 status 5 gt 6 dp 7 mqref 8 mqalt 9 mqdiff
#       10 mqdifft 11 af 12 ref_ad 13 alt_ad 14 vaf
awk -F'\t' 'NR>1 && $14!=""{n[$4]++; v[$4]+=$14; d[$4]+=$6; m[$4]+=$9}
  END{for(k in n) printf "%s: n=%d meanVAF=%.3f meanDP=%.0f meanMQDiff=%.2f\n", k,n[k],v[k]/n[k],d[k]/n[k],m[k]/n[k]}' "$TSV"