#!/usr/bin/env bash
# Accuracy benchmark: compare caller VCF(s) to a truth set (e.g. the GIAB
# HG002 high-confidence calls) within the benchmark BED, and report
# precision / recall / F1 split by variant class (SNP, indel).
#
#   benchmarks/lib/compare_to_truth.sh <bench.config.sh> <query.vcf>...
#
# With no query VCFs it defaults to the three single-sample caller
# outputs under OUT_ROOT (ours / gatk / freebayes), skipping any absent.
#
# Method (bcftools): both query and truth are restricted to the BED,
# left-aligned and biallelic-split (`norm -m -any`), then class-filtered
# and intersected (`isec`, matching on POS+REF+ALT). So this scores
# ALLELE concordance, not genotype concordance — a 0/1-vs-1/1 mismatch at
# a correctly-called site still counts as a TP. (Genotype-level eval
# needs rtg vcfeval or hap.py; this is the dependency-light analogue of
# tomato's agreement dashboard.) The truth set is restricted to FILTER
# PASS; query records keep their own FILTER (callers differ — note this
# when reading the numbers).
#
# Env:
#   TRUTH_VCF   truth callset (from config; overridable)
#   CLASSES     space-separated subset of "snps indels" (default both)

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

bench_load_config "${1:-}"
shift || true

command -v bcftools >/dev/null || { echo "bcftools not on PATH" >&2; exit 1; }
[[ -n "${TRUTH_VCF:-}" && -f "$TRUTH_VCF" ]] || {
    echo "TRUTH_VCF unset or missing: ${TRUTH_VCF:-<none>}" >&2; exit 1; }
bench_preflight "$BED" "$REFERENCE" "${REFERENCE}.fai"

CLASSES="${CLASSES:-snps indels}"

# Default query set: the three single-sample caller outputs.
queries=("$@")
if (( ${#queries[@]} == 0 )); then
    base="$(bench_sample_base "$(bench_single_cram)")"
    for c in ours gatk freebayes; do
        q="$OUT_ROOT/$c/single_${base}.vcf"
        [[ -f "$q" ]] && queries+=("$q")
    done
    if (( ${#queries[@]} == 0 )); then
        echo "no query VCFs given and none found under $OUT_ROOT/{ours,gatk,freebayes}/" >&2
        echo "run the callers first, or pass VCF paths explicitly." >&2
        exit 1
    fi
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# normalize <src-vcf> <class> <out.vcf.gz> — BED-restrict, left-align,
# biallelic-split, keep only <class>, bgzip + index. `extra_view` lets
# the truth call add `-f PASS`.
normalize() {
    local src="$1" cls="$2" out="$3"; shift 3
    bcftools view "$@" -T "$BED" "$src" -Ou \
      | bcftools norm -f "$REFERENCE" -m -any -Ou 2>/dev/null \
      | bcftools view -v "$cls" -Oz -o "$out"
    bcftools index -t "$out"
}

count_records() { bcftools view -H "$1" 2>/dev/null | wc -l | tr -d ' '; }

echo "truth     : $TRUTH_VCF (FILTER=PASS)"
echo "regions   : $BED ($(wc -l < "$BED") intervals)"
echo "reference : $REFERENCE"
echo "classes   : $CLASSES"
echo "queries   : ${#queries[@]}"
echo

# Machine-readable output for the dashboard (mirrors the perf_*.py ->
# perf_dashboard.py convention). One row per (caller, class).
TSV_OUT="${TSV_OUT:-$OUT_ROOT/comparison/accuracy.tsv}"
mkdir -p "$(dirname "$TSV_OUT")"
printf 'caller\tclass\ttp\tfp\tfn\tprecision\trecall\tf1\n' > "$TSV_OUT"

printf '%-26s %-7s %7s %7s %7s %9s %9s %7s\n' \
    caller class TP FP FN precision recall F1
printf '%s\n' "--------------------------------------------------------------------------------------"

for cls in $CLASSES; do
    truth_norm="$TMP/truth.$cls.vcf.gz"
    normalize "$TRUTH_VCF" "$cls" "$truth_norm" -f PASS
    for q in "${queries[@]}"; do
        name="$(basename "$(dirname "$q")")"   # ours / gatk / freebayes
        [[ "$name" == "$BENCH_NAME" || -z "$name" ]] && name="$(basename "$q" .vcf)"
        q_norm="$TMP/query.$name.$cls.vcf.gz"
        normalize "$q" "$cls" "$q_norm"

        isec_dir="$TMP/isec.$name.$cls"
        bcftools isec -p "$isec_dir" "$truth_norm" "$q_norm" >/dev/null 2>&1
        # 0000=truth-only(FN) 0001=query-only(FP) 0002=shared-from-truth(TP)
        local_fn=$(count_records "$isec_dir/0000.vcf")
        local_fp=$(count_records "$isec_dir/0001.vcf")
        local_tp=$(count_records "$isec_dir/0002.vcf")

        read -r prec rec f1 < <(awk -v tp="$local_tp" -v fp="$local_fp" -v fn="$local_fn" 'BEGIN{
            p = (tp+fp>0) ? tp/(tp+fp) : 0;
            r = (tp+fn>0) ? tp/(tp+fn) : 0;
            f = (p+r>0) ? 2*p*r/(p+r) : 0;
            printf "%.4f %.4f %.4f\n", p, r, f;
        }')
        printf '%-26s %-7s %7d %7d %7d %9s %9s %7s\n' \
            "$name" "$cls" "$local_tp" "$local_fp" "$local_fn" "$prec" "$rec" "$f1"
        printf '%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n' \
            "$name" "$cls" "$local_tp" "$local_fp" "$local_fn" "$prec" "$rec" "$f1" >> "$TSV_OUT"
    done
done

echo
echo "tsv: $TSV_OUT"
