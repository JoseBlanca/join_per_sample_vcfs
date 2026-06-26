#!/usr/bin/env bash
# Sweep the var-calling MAPQ-difference Welch's-t filter threshold
# (--min-mapq-diff-t) and measure SNP + indel concordance vs the GIAB truth
# for the per_sample dataset. Emits a tidy TSV for the dashboard.
#
#   benchmarks/giab/src/mapq_diff_sweep.sh [VARIANT] [COVERAGE]
#
# Re-calls each sample's existing .psp (no re-pileup) at every threshold in
# THRESHOLDS, restricted to the sample BED, with the DUST filter off (matching
# the indel-inclusive baseline). For each (threshold, sample, class) it runs
# the same bcftools isec concordance as the dashboard and appends a row to
#   results/per_sample/<cov>/mapq_diff_sweep.tsv
#
# Env:
#   THRESHOLDS   space-separated list (default "-3 -5 -10 -20 -inf")
#   BIN          caller binary (default target/release/pop_var_caller)

set -euo pipefail
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SRC_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$BENCH_DIR/../.." && pwd)"

VARIANT="${1:-ours_nobaq_nodust}"
COVERAGE="${2:-300x}"
THRESHOLDS="${THRESHOLDS:--3 -5 -10 -20 -inf}"
BIN="${BIN:-$PROJECT_ROOT/target/release/pop_var_caller}"

REFERENCE="$BENCH_DIR/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BED_DIR="$BENCH_DIR/per_sample/bed"
TRUTH_DIR="$BENCH_DIR/per_sample/vcf"
RES_DIR="$BENCH_DIR/results/per_sample/$COVERAGE/$VARIANT"
SWEEP_DIR="$RES_DIR/mapq_sweep"
TSV="$RES_DIR/mapq_diff_sweep.tsv"
SAMPLES=(HG002 HG003 HG004)

command -v bcftools >/dev/null || { echo "bcftools not on PATH" >&2; exit 1; }
[[ -x "$BIN" ]] || { echo "binary not found: $BIN" >&2; exit 1; }
mkdir -p "$SWEEP_DIR"
printf 'threshold\tsample\tclass\ttp\tfp\tfn\tprecision\trecall\tf1\n' > "$TSV"

norm() { # src class out passflag bed
    bcftools view $4 -T "$5" "$1" -Ou \
      | bcftools norm -f "$REFERENCE" -m -any -Ou 2>/dev/null \
      | bcftools view -v "$2" -Oz -o "$3"
    bcftools index -tf "$3"
}

for thr in $THRESHOLDS; do
    tag="${thr//-/m}"   # -inf -> minf, -3 -> m3 (filename-safe)
    for s in "${SAMPLES[@]}"; do
        psp="$RES_DIR/${s}.psp"
        bed="$BED_DIR/${s}_bench_azar_merged_100.bed"
        truth="$TRUTH_DIR/${s}_GRCh38_1_22_v4.2.1_benchmark.selected_100.vcf.gz"
        [[ -s "$psp" ]] || { echo "missing psp: $psp" >&2; continue; }
        q="$SWEEP_DIR/${s}.t${tag}.vcf"
        "$BIN" var-calling --reference "$REFERENCE" --regions "$bed" \
            --no-complexity-filter --min-mapq-diff-t="$thr" \
            --output "$q" --threads 4 "$psp" 2>/dev/null

        td="$(mktemp -d)"
        for cls in snps indels; do
            norm "$truth" "$cls" "$td/t.vcf.gz" "-f PASS" "$bed"
            norm "$q"     "$cls" "$td/q.vcf.gz" ""       "$bed"
            bcftools isec -p "$td/i" "$td/t.vcf.gz" "$td/q.vcf.gz" 2>/dev/null
            tp=$(bcftools view -H "$td/i/0002.vcf" | wc -l | tr -d ' ')
            fp=$(bcftools view -H "$td/i/0001.vcf" | wc -l | tr -d ' ')
            fn=$(bcftools view -H "$td/i/0000.vcf" | wc -l | tr -d ' ')
            awk -v thr="$thr" -v s="$s" -v cls="$cls" -v tp="$tp" -v fp="$fp" -v fn="$fn" 'BEGIN{
                p=(tp+fp)?tp/(tp+fp):0; r=(tp+fn)?tp/(tp+fn):0; f=(p+r)?2*p*r/(p+r):0;
                printf "%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n", thr,s,cls,tp,fp,fn,p,r,f;
            }' >> "$TSV"
            rm -rf "$td/i"
        done
        rm -rf "$td"
    done
done

echo "wrote $TSV"
column -t "$TSV"
