#!/usr/bin/env bash
# pop_var_caller SSR runner — microsatellite genotyping, parallel to run_ours.sh.
#
#   benchmarks/lib/run_ours_ssr.sh <bench.config.sh> [single|cohort]
#
# Pipeline (three stages, vs the SNP path's two):
#   ssr-catalog : reference FASTA -> per-genome SSR locus catalog (via trf-mod).
#                 Built ONCE and cached at $CATALOG (shared by both modes and
#                 every sample). Delete it to force a rebuild.
#   ssr-pileup  : CRAM/BAM -> per-sample .ssr.psp evidence (one per sample).
#   ssr-call    : the per-sample .ssr.psp files -> multi-sample VCF.
#
# single : catalog + ssr-pileup(SINGLE_CRAM) + ssr-call on that one sample.
# cohort : catalog + ssr-pileup per sample (sequential, cached) + ssr-call over
#          all samples.
#
# Like the SNP runner, our caller has no region flag — the benchmark CRAMs are
# pre-sliced to the region set, and ssr-catalog/ssr-pileup read the whole input.
#
# Env overrides (see common.sh + bench.config.sh):
#   POP_VAR_CALLER_BIN, THREADS, REFERENCE, CATALOG, TRF_MOD_PATH,
#   CATALOG_EXTRA, PILEUP_EXTRA, CALL_EXTRA, DRY_RUN=1

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

CONFIG="${1:-}"
MODE="${2:-single}"
bench_load_config "$CONFIG"
bench_discover_ours_bin

OUT_DIR="$OUT_ROOT/ours"
# CATALOG is resolved here (OUT_ROOT is only defined after the config is sourced).
CATALOG="${CATALOG:-$OUT_DIR/$BENCH_NAME.ssr.catalog}"
CATALOG_EXTRA=${CATALOG_EXTRA:-}
PILEUP_EXTRA=${PILEUP_EXTRA:-}
CALL_EXTRA=${CALL_EXTRA:-}

# trf-mod path flag, only if the config set one (else ssr-catalog uses PATH).
TRF_FLAG=()
[[ -n "${TRF_MOD_PATH:-}" ]] && TRF_FLAG=(--trf-mod-path "$TRF_MOD_PATH")

# Build the SSR catalog once (cached). Shared by every sample and both modes.
ensure_catalog() {
    bench_preflight "$REFERENCE" "${REFERENCE}.fai"
    if [[ -s "$CATALOG" ]]; then
        echo "[skip catalog] $CATALOG already present ($(bench_record_count "$CATALOG") loci)"
        return 0
    fi
    mkdir -p "$(dirname "$CATALOG")"
    local log="${CATALOG%.ssr.catalog}.catalog.log"
    echo "[ssr-catalog] $REFERENCE -> $CATALOG (trf-mod; one-time, whole genome)"
    local t0 t1
    t0=$(bench_now)
    # shellcheck disable=SC2086
    bench_run "$log" -- "$POP_VAR_CALLER_BIN" ssr-catalog \
        --reference "$REFERENCE" \
        --output "$CATALOG" \
        ${TRF_FLAG[@]+"${TRF_FLAG[@]}"} \
        $CATALOG_EXTRA
    t1=$(bench_now)
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0
    echo "catalog: $((t1 - t0)) s, $(bench_record_count "$CATALOG") loci"
}

# pileup_one <cram> <psp> <log>
pileup_one() {
    local cram="$1" psp="$2" log="$3"
    # shellcheck disable=SC2086
    bench_run "$log" -- "$POP_VAR_CALLER_BIN" ssr-pileup \
        --reference "$REFERENCE" \
        --catalog "$CATALOG" \
        --output "$psp" \
        --threads "$THREADS" \
        $PILEUP_EXTRA \
        "$cram"
}

run_single() {
    local cram base psp out_vcf
    cram="$(bench_single_cram)"
    base="$(bench_sample_base "$cram")"
    bench_preflight "$cram" "${cram}.crai"
    mkdir -p "$OUT_DIR"

    echo "binary    : $POP_VAR_CALLER_BIN"
    echo "mode      : single ($BENCH_NAME) — catalog -> ssr-pileup -> ssr-call"
    echo "sample    : $base"
    echo "reference : $REFERENCE"
    echo "catalog   : $CATALOG"
    echo "threads   : $THREADS"
    echo

    ensure_catalog

    psp="$OUT_DIR/single_${base}.ssr.psp"
    out_vcf="$OUT_DIR/single_${base}.ssr.vcf"
    local t0 t_mid t1
    t0=$(bench_now)
    echo "[ssr-pileup] $base -> $psp"
    pileup_one "$cram" "$psp" "$OUT_DIR/single_${base}.pileup.log"
    t_mid=$(bench_now)
    echo "[ssr-call] $psp -> $out_vcf"
    # shellcheck disable=SC2086
    bench_run "$OUT_DIR/single_${base}.log" -- "$POP_VAR_CALLER_BIN" ssr-call \
        --catalog "$CATALOG" \
        --output "$out_vcf" \
        --threads "$THREADS" \
        $CALL_EXTRA \
        "$psp"
    t1=$(bench_now)
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0
    echo
    echo "ssr-pileup : $((t_mid - t0)) s"
    echo "ssr-call   : $((t1 - t_mid)) s"
    echo "records    : $(bench_record_count "$out_vcf") (single sample)"
}

run_cohort() {
    bench_list_crams
    local out_dir="$OUT_DIR/cohort"
    local psp_dir="$out_dir/psp"
    local cohort_vcf="$out_dir/cohort.ssr.vcf"
    mkdir -p "$psp_dir"

    echo "binary    : $POP_VAR_CALLER_BIN"
    echo "mode      : cohort ($BENCH_NAME)"
    echo "samples   : ${#BENCH_CRAMS[@]}"
    echo "reference : $REFERENCE"
    echo "catalog   : $CATALOG"
    echo "threads   : $THREADS"
    echo "vcf out   : $cohort_vcf"
    echo

    ensure_catalog

    # ---- ssr-pileup per sample (sequential, cached) ----
    local t_s1 t_s1_end cram base psp log
    t_s1=$(bench_now)
    for cram in "${BENCH_CRAMS[@]}"; do
        base="$(bench_sample_base "$cram")"
        psp="$psp_dir/${base}.ssr.psp"
        log="$psp_dir/${base}.log"
        if [[ -s "$psp" ]]; then
            echo "[skip pileup] $base (.ssr.psp already present)"
            continue
        fi
        echo "[ssr-pileup] $base"
        pileup_one "$cram" "$psp" "$log"
    done
    t_s1_end=$(bench_now)
    echo
    echo "stage 1 (pileup) elapsed: $((t_s1_end - t_s1)) s"
    echo

    # ---- ssr-call cohort ----
    shopt -s nullglob
    local psps=("$psp_dir"/*.ssr.psp)
    shopt -u nullglob
    if [[ "${DRY_RUN:-0}" != "1" && ${#psps[@]} -ne ${#BENCH_CRAMS[@]} ]]; then
        echo "WARNING: expected ${#BENCH_CRAMS[@]} .ssr.psp files, found ${#psps[@]}" >&2
    fi
    echo "[ssr-call] ${#psps[@]} samples -> $cohort_vcf"
    local t_s2 t_s2_end
    t_s2=$(bench_now)
    # shellcheck disable=SC2086
    bench_run "$out_dir/cohort.log" -- "$POP_VAR_CALLER_BIN" ssr-call \
        --catalog "$CATALOG" \
        --output "$cohort_vcf" \
        --threads "$THREADS" \
        $CALL_EXTRA \
        "${psps[@]}"
    t_s2_end=$(bench_now)
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0
    echo
    echo "stage 2 (call) elapsed: $((t_s2_end - t_s2)) s"
    echo "total elapsed         : $((t_s2_end - t_s1)) s"
    echo "records: $(bench_record_count "$cohort_vcf") sites, $(bench_sample_count "$cohort_vcf") samples"
}

case "$MODE" in
    single) run_single ;;
    cohort) run_cohort ;;
    *) echo "unknown mode: $MODE (expected single|cohort)" >&2; exit 2 ;;
esac
