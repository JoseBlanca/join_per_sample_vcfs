#!/usr/bin/env bash
# freebayes runner — shared across benchmarks.
#
#   benchmarks/lib/run_freebayes.sh <bench.config.sh> [single|cohort]
#
# single : one CRAM -> single-sample VCF, restricted to the benchmark BED.
# cohort : all CRAMs -> joint VCF in one invocation (freebayes has no
#          per-sample intermediate; sample ids come from each @RG SM tag).
#
# A post-call QUAL >= MIN_QUAL filter is applied inline with awk
# (freebayes emits a long low-QUAL tail the other callers don't; leaving
# it in skews any cross-caller comparison). awk avoids a bcftools
# dependency on the host.
#
# Env overrides (see common.sh for the rest):
#   FREEBAYES_BIN   binary (default: auto-detect `freebayes` on PATH)
#   REFERENCE       FASTA (.fai sibling required)
#   PLOIDY          --ploidy (default 2)
#   MIN_QUAL        QUAL floor applied pre-write (default 30)
#   EXTRA_ARGS      appended verbatim to the freebayes command line
#   PARALLEL=1      cohort only: fan out per-BED-region via
#                   freebayes-parallel + GNU parallel (needs both on PATH)
#   THREADS         worker count when PARALLEL=1 (default 4)
#   DRY_RUN=1       print commands instead of running them

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

CONFIG="${1:-}"
MODE="${2:-single}"
bench_load_config "$CONFIG"
bench_discover_freebayes_bin

EXTRA_ARGS=${EXTRA_ARGS:-}
PARALLEL="${PARALLEL:-0}"
OUT_DIR="$OUT_ROOT/freebayes"
FB_VERSION="$($FREEBAYES_BIN --version 2>&1 | head -1 || true)"

# header lines pass through; data rows kept iff QUAL (col 6) is numeric
# and >= min. `$6 + 0` coerces; non-numeric falls to 0 and is dropped.
FILTER_AWK='/^#/ { print; next } $6 != "." && ($6 + 0) >= min { print }'

run_single() {
    local cram base out_vcf log
    cram="$(bench_single_cram)"
    base="$(bench_sample_base "$cram")"
    bench_preflight "$cram" "${cram}.crai" "$BED" "$REFERENCE" "${REFERENCE}.fai"

    out_vcf="$OUT_DIR/single_${base}.vcf"
    log="$OUT_DIR/single_${base}.log"
    mkdir -p "$OUT_DIR"

    echo "binary    : $FREEBAYES_BIN ($FB_VERSION)"
    echo "mode      : single ($BENCH_NAME)"
    echo "sample    : $base"
    echo "input     : $cram"
    echo "reference : $REFERENCE"
    echo "regions   : $BED ($(wc -l < "$BED") intervals)"
    echo "ploidy    : $PLOIDY"
    echo "min QUAL  : $MIN_QUAL (records below are dropped pre-write)"
    echo "output    : $out_vcf"
    echo

    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        printf 'DRY-RUN:'
        printf ' %q' "$FREEBAYES_BIN" -f "$REFERENCE" -t "$BED" -p "$PLOIDY" $EXTRA_ARGS "$cram"
        printf ' | awk -v min=%q %q > %q\n' "$MIN_QUAL" "$FILTER_AWK" "$out_vcf"
        return 0
    fi

    local t0 t1
    t0=$(bench_now)
    # shellcheck disable=SC2086
    "$FREEBAYES_BIN" \
        -f "$REFERENCE" \
        -t "$BED" \
        -p "$PLOIDY" \
        $EXTRA_ARGS \
        "$cram" \
        2> >(tee "$log" >&2) \
      | awk -v min="$MIN_QUAL" "$FILTER_AWK" > "$out_vcf"
    t1=$(bench_now)
    echo
    echo "elapsed: $((t1 - t0)) s"
    echo "records: $(bench_record_count "$out_vcf") (QUAL >= $MIN_QUAL, no other filters)"
}

run_cohort() {
    bench_list_crams
    bench_preflight "$BED" "$REFERENCE" "${REFERENCE}.fai"
    local c
    for c in "${BENCH_CRAMS[@]}"; do
        bench_preflight "${c}.crai"
    done

    local out_vcf="$OUT_DIR/cohort.vcf"
    local log="$OUT_DIR/cohort.log"
    mkdir -p "$OUT_DIR"

    echo "binary    : $FREEBAYES_BIN ($FB_VERSION)"
    echo "mode      : cohort ($BENCH_NAME)"
    echo "samples   : ${#BENCH_CRAMS[@]}"
    echo "reference : $REFERENCE"
    echo "regions   : $BED ($(wc -l < "$BED") intervals)"
    echo "ploidy    : $PLOIDY"
    echo "min QUAL  : $MIN_QUAL (records below are dropped pre-write)"
    echo "mode      : $([[ "$PARALLEL" == "1" ]] && echo "parallel (threads=$THREADS)" || echo "single-threaded")"
    echo "output    : $out_vcf"
    echo

    local t0 t1
    t0=$(bench_now)
    if [[ "$PARALLEL" == "1" ]]; then
        command -v freebayes-parallel >/dev/null || {
            echo "PARALLEL=1 requested but freebayes-parallel not on PATH" >&2; exit 1; }
        command -v parallel >/dev/null || {
            echo "PARALLEL=1 requested but GNU parallel not on PATH" >&2; exit 1; }
        # freebayes-parallel takes a regions file (chr:start-end, 1-based,
        # NOT BED). Translate the BED on the fly.
        local regions_txt="$OUT_DIR/regions.fbparallel.txt"
        awk '{ printf "%s:%d-%d\n", $1, $2+1, $3 }' "$BED" > "$regions_txt"
        if [[ "${DRY_RUN:-0}" == "1" ]]; then
            printf 'DRY-RUN:'
            printf ' %q' freebayes-parallel "$regions_txt" "$THREADS" \
                -f "$REFERENCE" -p "$PLOIDY" $EXTRA_ARGS "${BENCH_CRAMS[@]}"
            printf ' | awk -v min=%q %q > %q\n' "$MIN_QUAL" "$FILTER_AWK" "$out_vcf"
            return 0
        fi
        # shellcheck disable=SC2086
        freebayes-parallel \
            "$regions_txt" \
            "$THREADS" \
            -f "$REFERENCE" \
            -p "$PLOIDY" \
            $EXTRA_ARGS \
            "${BENCH_CRAMS[@]}" \
            2> >(tee "$log" >&2) \
          | awk -v min="$MIN_QUAL" "$FILTER_AWK" > "$out_vcf"
    else
        if [[ "${DRY_RUN:-0}" == "1" ]]; then
            printf 'DRY-RUN:'
            printf ' %q' "$FREEBAYES_BIN" -f "$REFERENCE" -t "$BED" -p "$PLOIDY" \
                $EXTRA_ARGS "${BENCH_CRAMS[@]}"
            printf ' | awk -v min=%q %q > %q\n' "$MIN_QUAL" "$FILTER_AWK" "$out_vcf"
            return 0
        fi
        # shellcheck disable=SC2086
        "$FREEBAYES_BIN" \
            -f "$REFERENCE" \
            -t "$BED" \
            -p "$PLOIDY" \
            $EXTRA_ARGS \
            "${BENCH_CRAMS[@]}" \
            2> >(tee "$log" >&2) \
          | awk -v min="$MIN_QUAL" "$FILTER_AWK" > "$out_vcf"
    fi
    t1=$(bench_now)
    echo
    echo "elapsed: $((t1 - t0)) s"
    echo "records: $(bench_record_count "$out_vcf") (QUAL >= $MIN_QUAL, no other filters)"
    echo "samples in vcf: $(bench_sample_count "$out_vcf")"
}

case "$MODE" in
    single) run_single ;;
    cohort) run_cohort ;;
    *) echo "unknown mode: $MODE (expected single|cohort)" >&2; exit 2 ;;
esac
