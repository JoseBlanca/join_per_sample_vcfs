#!/usr/bin/env bash
# pop_var_caller runner — shared across benchmarks.
#
#   benchmarks/lib/run_ours.sh <bench.config.sh> [single|cohort]
#
# single : one CRAM -> single-sample VCF via the two-stage path —
#          `pileup` (CRAM -> .psp) then `var-calling` (.psp -> VCF). The
#          CRAM is config's SINGLE_CRAM. (The old one-shot direct path,
#          `var-calling-from-bam`, was removed; pileup->psp->var-calling
#          is now the only route to a VCF.)
# cohort : all CRAMs -> per-sample .psp (`pileup`) -> joined multi-sample
#          VCF (`var-calling`). Pileups run sequentially (each at THREADS
#          workers) to stay within the host P-core budget; already-built
#          .psp files are skipped on re-runs.
#
# Our caller has no BED/region flag (verified against the CLI) — the
# benchmark CRAMs are already pre-sliced to the region set, so it reads
# the whole CRAM. GATK and freebayes restrict via the BED for fairness.
#
# Env overrides (see common.sh for the rest):
#   POP_VAR_CALLER_BIN  binary to invoke (default: auto-detect)
#   THREADS             rayon worker count (default: 4)
#   REFERENCE           reference FASTA (.fai sibling required)
#   PILEUP_EXTRA        appended to each `pileup` invocation
#   VARCALL_EXTRA       appended to the `var-calling` invocation
#                       (e.g. --no-complexity-filter to disable the DUST
#                       low-complexity filter for a fair cross-caller cmp)
#   DRY_RUN=1           print commands instead of running them

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

CONFIG="${1:-}"
MODE="${2:-single}"
bench_load_config "$CONFIG"
bench_discover_ours_bin

PILEUP_EXTRA=${PILEUP_EXTRA:-}
VARCALL_EXTRA=${VARCALL_EXTRA:-}
OUT_DIR="$OUT_ROOT/ours"

run_single() {
    local cram base psp out_vcf psp_log log
    cram="$(bench_single_cram)"
    base="$(bench_sample_base "$cram")"
    bench_preflight "$cram" "${cram}.crai" "$REFERENCE" "${REFERENCE}.fai"

    psp="$OUT_DIR/single_${base}.psp"
    psp_log="$OUT_DIR/single_${base}.pileup.log"
    out_vcf="$OUT_DIR/single_${base}.vcf"
    log="$OUT_DIR/single_${base}.log"
    mkdir -p "$OUT_DIR"

    echo "binary    : $POP_VAR_CALLER_BIN"
    echo "mode      : single ($BENCH_NAME) — pileup -> .psp -> var-calling"
    echo "sample    : $base"
    echo "input     : $cram"
    echo "reference : $REFERENCE"
    echo "threads   : $THREADS"
    echo "psp       : $psp"
    echo "output    : $out_vcf"
    echo

    local t0 t_mid t1
    t0=$(bench_now)
    # ---- Stage 1: pileup (CRAM -> .psp) ----
    echo "[pileup] $base -> $psp"
    # shellcheck disable=SC2086
    bench_run "$psp_log" -- "$POP_VAR_CALLER_BIN" pileup \
        --reference "$REFERENCE" \
        --output "$psp" \
        --threads "$THREADS" \
        $PILEUP_EXTRA \
        "$cram"
    t_mid=$(bench_now)

    # ---- Stage 2: var-calling (.psp -> VCF) ----
    echo "[var-calling] $psp -> $out_vcf"
    # shellcheck disable=SC2086
    bench_run "$log" -- "$POP_VAR_CALLER_BIN" var-calling \
        --reference "$REFERENCE" \
        --output "$out_vcf" \
        --threads "$THREADS" \
        $VARCALL_EXTRA \
        "$psp"
    t1=$(bench_now)
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0
    echo
    echo "stage 1 (pileup)      : $((t_mid - t0)) s"
    echo "stage 2 (var-calling) : $((t1 - t_mid)) s"
    echo "total   elapsed       : $((t1 - t0)) s"
    echo "records: $(bench_record_count "$out_vcf") (single sample)"
}

run_cohort() {
    bench_list_crams
    bench_preflight "$REFERENCE" "${REFERENCE}.fai"

    local out_dir="$OUT_DIR/cohort"
    local psp_dir="$out_dir/psp"
    local cohort_vcf="$out_dir/cohort.vcf"
    local cohort_log="$out_dir/cohort.log"
    mkdir -p "$psp_dir"

    echo "binary    : $POP_VAR_CALLER_BIN"
    echo "mode      : cohort ($BENCH_NAME)"
    echo "samples   : ${#BENCH_CRAMS[@]}"
    echo "reference : $REFERENCE"
    echo "threads   : $THREADS (per pileup, and for cohort var-calling)"
    echo "psp dir   : $psp_dir"
    echo "vcf out   : $cohort_vcf"
    echo

    # ---- Stage 1: per-sample pileup (sequential) ----
    local t_s1 t_s1_end cram base psp log
    t_s1=$(bench_now)
    for cram in "${BENCH_CRAMS[@]}"; do
        base="$(bench_sample_base "$cram")"
        psp="$psp_dir/${base}.psp"
        log="$psp_dir/${base}.log"
        if [[ -s "$psp" ]]; then
            echo "[skip pileup] $base (.psp already present)"
            continue
        fi
        echo "[pileup] $base"
        # shellcheck disable=SC2086
        bench_run "$log" -- "$POP_VAR_CALLER_BIN" pileup \
            --reference "$REFERENCE" \
            --output "$psp" \
            --threads "$THREADS" \
            $PILEUP_EXTRA \
            "$cram"
    done
    t_s1_end=$(bench_now)
    echo
    echo "stage 1 elapsed: $((t_s1_end - t_s1)) s"
    echo

    # ---- Stage 2: cohort var-calling ----
    shopt -s nullglob
    local psps=("$psp_dir"/*.psp)
    shopt -u nullglob
    if [[ "${DRY_RUN:-0}" != "1" && ${#psps[@]} -ne ${#BENCH_CRAMS[@]} ]]; then
        echo "WARNING: expected ${#BENCH_CRAMS[@]} .psp files, found ${#psps[@]}" >&2
    fi

    echo "[var-calling] ${#psps[@]} samples -> $cohort_vcf"
    local t_s2 t_s2_end
    t_s2=$(bench_now)
    # shellcheck disable=SC2086
    bench_run "$cohort_log" -- "$POP_VAR_CALLER_BIN" var-calling \
        --reference "$REFERENCE" \
        --output "$cohort_vcf" \
        --threads "$THREADS" \
        $VARCALL_EXTRA \
        "${psps[@]}"
    t_s2_end=$(bench_now)
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0
    echo
    echo "stage 2 elapsed: $((t_s2_end - t_s2)) s"
    echo "total   elapsed: $((t_s2_end - t_s1)) s"
    echo "records: $(bench_record_count "$cohort_vcf") (cohort VCF, all sites)"
}

case "$MODE" in
    single) run_single ;;
    cohort) run_cohort ;;
    *) echo "unknown mode: $MODE (expected single|cohort)" >&2; exit 2 ;;
esac
