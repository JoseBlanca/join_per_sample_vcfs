#!/usr/bin/env bash
# GATK runner — shared across benchmarks.
#
#   benchmarks/lib/run_gatk.sh <bench.config.sh> [single|cohort]
#
# single : one CRAM -> single-sample VCF via HaplotypeCaller in default
#          (non-GVCF) mode, restricted to the benchmark BED.
# cohort : per-sample GVCFs (HaplotypeCaller -ERC GVCF) -> CombineGVCFs
#          -> GenotypeGVCFs, the GATK best-practices three-stage flow.
#          Stage 1 runs sequentially (each at THREADS native-pair-hmm
#          threads); already-built GVCFs are skipped on re-runs.
#
# GATK uniquely requires a `.dict` next to the FASTA (Picard sequence
# dictionary) — built by benchmarks/lib/prepare_reference.sh.
#
# Env overrides (see common.sh for the rest):
#   GATK_BIN            binary (default: /opt/gatk/gatk — in the container)
#   REFERENCE          FASTA (.fai + .dict siblings required)
#   THREADS            --native-pair-hmm-threads (default 4)
#   GATK_HEAP          -Xmx for HC / GenotypeGVCFs (default 4g)
#   GATK_COMBINE_HEAP  -Xmx for CombineGVCFs (default 8g)
#   HC_EXTRA           appended to each HaplotypeCaller invocation
#   COMBINE_EXTRA      appended to CombineGVCFs
#   GENOTYPE_EXTRA     appended to GenotypeGVCFs
#   DRY_RUN=1          print commands instead of running them

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

CONFIG="${1:-}"
MODE="${2:-single}"
bench_load_config "$CONFIG"
bench_discover_gatk_bin

HC_EXTRA=${HC_EXTRA:-}
COMBINE_EXTRA=${COMBINE_EXTRA:-}
GENOTYPE_EXTRA=${GENOTYPE_EXTRA:-}
REF_DICT="$(bench_ref_dict)"
OUT_DIR="$OUT_ROOT/gatk"
GATK_VERSION="$($GATK_BIN --version 2>&1 | head -1 || true)"

run_single() {
    local cram base out_vcf log
    cram="$(bench_single_cram)"
    base="$(bench_sample_base "$cram")"
    bench_preflight "$cram" "${cram}.crai" "$BED" "$REFERENCE" "${REFERENCE}.fai" "$REF_DICT"

    out_vcf="$OUT_DIR/single_${base}.vcf"
    log="$OUT_DIR/single_${base}.log"
    mkdir -p "$OUT_DIR"

    echo "binary    : $GATK_BIN ($GATK_VERSION)"
    echo "mode      : single ($BENCH_NAME)"
    echo "sample    : $base"
    echo "input     : $cram"
    echo "reference : $REFERENCE"
    echo "regions   : $BED ($(wc -l < "$BED") intervals)"
    echo "threads   : $THREADS (--native-pair-hmm-threads)"
    echo "heap      : $GATK_HEAP"
    echo "output    : $out_vcf"
    echo

    local t0 t1
    t0=$(bench_now)
    # shellcheck disable=SC2086
    bench_run "$log" -- "$GATK_BIN" --java-options "-Xmx${GATK_HEAP}" HaplotypeCaller \
        --reference "$REFERENCE" \
        --input "$cram" \
        --intervals "$BED" \
        --output "$out_vcf" \
        --native-pair-hmm-threads "$THREADS" \
        $HC_EXTRA
    t1=$(bench_now)
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0
    echo
    echo "elapsed: $((t1 - t0)) s"
    echo "records: $(bench_record_count "$out_vcf")"
}

run_cohort() {
    bench_list_crams
    bench_preflight "$BED" "$REFERENCE" "${REFERENCE}.fai" "$REF_DICT"
    local c
    for c in "${BENCH_CRAMS[@]}"; do
        bench_preflight "${c}.crai"
    done

    local out_dir="$OUT_DIR/cohort"
    local gvcf_dir="$out_dir/gvcf"
    local combined_gvcf="$out_dir/combined.g.vcf.gz"
    local cohort_vcf="$out_dir/cohort.vcf"
    local combine_log="$out_dir/combine.log"
    local genotype_log="$out_dir/genotype.log"
    mkdir -p "$gvcf_dir"

    echo "binary    : $GATK_BIN ($GATK_VERSION)"
    echo "mode      : cohort ($BENCH_NAME)"
    echo "samples   : ${#BENCH_CRAMS[@]}"
    echo "reference : $REFERENCE"
    echo "regions   : $BED ($(wc -l < "$BED") intervals)"
    echo "threads   : $THREADS (per HaplotypeCaller)"
    echo "heap      : $GATK_HEAP (HC/genotype), $GATK_COMBINE_HEAP (combine)"
    echo "gvcf dir  : $gvcf_dir"
    echo "vcf out   : $cohort_vcf"
    echo

    # ---- Stage 1: per-sample HaplotypeCaller GVCF (sequential) ----
    local t_s1 t_s1_end cram base gvcf log
    t_s1=$(bench_now)
    for cram in "${BENCH_CRAMS[@]}"; do
        base="$(bench_sample_base "$cram")"
        gvcf="$gvcf_dir/${base}.g.vcf.gz"
        log="$gvcf_dir/${base}.log"
        if [[ -s "$gvcf" && -s "${gvcf}.tbi" ]]; then
            echo "[skip HC] $base (GVCF already present)"
            continue
        fi
        echo "[HC] $base"
        # shellcheck disable=SC2086
        bench_run "$log" -- "$GATK_BIN" --java-options "-Xmx${GATK_HEAP}" HaplotypeCaller \
            --reference "$REFERENCE" \
            --input "$cram" \
            --intervals "$BED" \
            --output "$gvcf" \
            --emit-ref-confidence GVCF \
            --native-pair-hmm-threads "$THREADS" \
            $HC_EXTRA
    done
    t_s1_end=$(bench_now)
    echo
    echo "stage 1 elapsed: $((t_s1_end - t_s1)) s"
    echo

    # ---- Stage 2: CombineGVCFs ----
    shopt -s nullglob
    local gvcfs=("$gvcf_dir"/*.g.vcf.gz)
    shopt -u nullglob
    if [[ "${DRY_RUN:-0}" != "1" && ${#gvcfs[@]} -ne ${#BENCH_CRAMS[@]} ]]; then
        echo "WARNING: expected ${#BENCH_CRAMS[@]} GVCFs, found ${#gvcfs[@]}" >&2
    fi
    local combine_args=() g
    for g in "${gvcfs[@]}"; do
        combine_args+=(--variant "$g")
    done

    echo "[CombineGVCFs] ${#gvcfs[@]} GVCFs -> $combined_gvcf (heap=$GATK_COMBINE_HEAP)"
    local t_s2 t_s2_end
    t_s2=$(bench_now)
    # shellcheck disable=SC2086
    bench_run "$combine_log" -- "$GATK_BIN" --java-options "-Xmx${GATK_COMBINE_HEAP}" CombineGVCFs \
        --reference "$REFERENCE" \
        --intervals "$BED" \
        "${combine_args[@]}" \
        --output "$combined_gvcf" \
        $COMBINE_EXTRA
    t_s2_end=$(bench_now)
    echo "stage 2 elapsed: $((t_s2_end - t_s2)) s"
    echo

    # ---- Stage 3: GenotypeGVCFs ----
    echo "[GenotypeGVCFs] -> $cohort_vcf"
    local t_s3 t_s3_end
    t_s3=$(bench_now)
    # shellcheck disable=SC2086
    bench_run "$genotype_log" -- "$GATK_BIN" --java-options "-Xmx${GATK_HEAP}" GenotypeGVCFs \
        --reference "$REFERENCE" \
        --intervals "$BED" \
        --variant "$combined_gvcf" \
        --output "$cohort_vcf" \
        $GENOTYPE_EXTRA
    t_s3_end=$(bench_now)
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0
    echo
    echo "stage 3 elapsed: $((t_s3_end - t_s3)) s"
    echo "total   elapsed: $((t_s3_end - t_s1)) s"
    echo "records: $(bench_record_count "$cohort_vcf")"
    echo "samples in vcf: $(bench_sample_count "$cohort_vcf")"
}

case "$MODE" in
    single) run_single ;;
    cohort) run_cohort ;;
    *) echo "unknown mode: $MODE (expected single|cohort)" >&2; exit 2 ;;
esac
