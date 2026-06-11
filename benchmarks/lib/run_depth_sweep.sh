#!/usr/bin/env bash
# Depth sweep for the QUAL-distribution analysis (qual-analysis branch).
#
#   benchmarks/lib/run_depth_sweep.sh <bench.config.sh> <subcommand>
#
# Downsamples the single benchmark CRAM to a ladder of global coverage
# depths, re-calls all three callers at each depth, scores each against
# the truth set, and merges the per-depth QUAL/DP tables into two tidy
# files the dashboard reads. The point: study how TP-vs-FP QUAL
# separation moves with sequencing depth (the native CRAM is ~301x, far
# above realistic WGS, so depth has to be created by subsampling).
#
# Subcommands:
#   crams      samtools subsample the native CRAM to each depth (host)
#   ours       run pop_var_caller (pileup -> psp -> var-calling) per depth (host)
#   freebayes  run freebayes per depth (host)
#   gatk       run GATK HaplotypeCaller per depth (DEV CONTAINER ONLY)
#   compare    compare_to_truth.sh per depth -> per-depth qual_dist/accuracy (host)
#   merge      concatenate per-depth tables with a `depth` column (host)
#   host       crams + ours + freebayes (everything runnable on the host)
#
# GATK is not on the host here; run that phase inside the dev container:
#   DEV_EXTRA_MOUNT=$HOME/genomes ./scripts/dev.sh \
#       bash benchmarks/lib/run_depth_sweep.sh <config> gatk
#
# Typical full run from the host:
#   bash benchmarks/lib/run_depth_sweep.sh <config> host
#   DEV_EXTRA_MOUNT=$HOME/genomes ./scripts/dev.sh \
#       bash benchmarks/lib/run_depth_sweep.sh <config> gatk
#   bash benchmarks/lib/run_depth_sweep.sh <config> compare
#   bash benchmarks/lib/run_depth_sweep.sh <config> merge
#
# Env overrides:
#   DEPTHS       space-separated ladder, "full" = native CRAM
#                (default: "5 10 15 20 30 50 100 full")
#   MEAN_DEPTH   measured native mean depth over the BED, used to turn a
#                target depth into a subsample fraction (default 301.35)
#   SUBSAMPLE_SEED  samtools -s seed (default 42; reproducible draws)
#   SWEEP_ROOT   output root (default $OUT_ROOT/depth_sweep)
# Plus everything the per-caller run scripts honor (THREADS, PILEUP_EXTRA,
# VARCALL_EXTRA, HC_EXTRA, MIN_QUAL, ...).

set -euo pipefail
LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$LIB_DIR/common.sh"

CONFIG="${1:-}"
SUBCMD="${2:-}"
bench_load_config "$CONFIG"

DEPTHS="${DEPTHS:-5 10 15 20 30 50 100 full}"
MEAN_DEPTH="${MEAN_DEPTH:-301.35}"
SUBSAMPLE_SEED="${SUBSAMPLE_SEED:-42}"
SWEEP_ROOT="${SWEEP_ROOT:-$OUT_ROOT/depth_sweep}"
CRAM_DEPTH_DIR="$CRAM_DIR/depth"

NATIVE_CRAM="$(bench_single_cram)"
NATIVE_BASE="$(bench_sample_base "$NATIVE_CRAM")"   # HG002_reads_selected_1000_rg

# CRAM path for a depth ("full" -> native CRAM untouched).
depth_cram() {
    local d="$1"
    if [[ "$d" == "full" ]]; then printf '%s\n' "$NATIVE_CRAM"; return; fi
    printf '%s/HG002_d%sx.cram\n' "$CRAM_DEPTH_DIR" "$d"
}
depth_out_root() { printf '%s/d%s\n' "$SWEEP_ROOT" "$1"; }

# samtools -s argument: integer part = seed, fraction = kept proportion.
# frac = target_depth / native_mean_depth.
subsample_arg() {
    awk -v s="$SUBSAMPLE_SEED" -v d="$1" -v m="$MEAN_DEPTH" \
        'BEGIN{ printf "%.6f", s + d/m }'
}

phase_crams() {
    mkdir -p "$CRAM_DEPTH_DIR"
    for d in $DEPTHS; do
        [[ "$d" == "full" ]] && { echo "[crams] $d -> native CRAM (no subsample)"; continue; }
        local out arg
        out="$(depth_cram "$d")"
        arg="$(subsample_arg "$d")"
        if [[ -s "$out" && -s "${out}.crai" ]]; then
            echo "[skip crams] ${d}x ($out present)"; continue
        fi
        echo "[crams] ${d}x  frac=${arg#${SUBSAMPLE_SEED}}  -> $out"
        # version=3.0: the container's GATK (older htsjdk) rejects CRAM 3.1,
        # which is samtools' default output; the native CRAM is 3.0.
        samtools view --output-fmt cram,version=3.0 -T "$REFERENCE" -s "$arg" -o "$out" "$NATIVE_CRAM"
        samtools index "$out"
        echo "  mean depth: $(samtools depth -a -b "$BED" --reference "$REFERENCE" "$out" \
            | awk '{s+=$3;n++} END{printf "%.1fx over %d bp\n", (n?s/n:0), n}')"
    done
}

# Run a per-caller script for every depth with per-depth SINGLE_CRAM/OUT_ROOT.
run_caller_phase() {
    local script="$1"; shift
    for d in $DEPTHS; do
        local cram root
        cram="$(depth_cram "$d")"; root="$(depth_out_root "$d")"
        [[ -s "$cram" ]] || { echo "MISSING cram for ${d}x: $cram (run 'crams' first)" >&2; exit 1; }
        echo "=================== ${d}x : $(basename "$script") ==================="
        SINGLE_CRAM="$cram" OUT_ROOT="$root" bash "$LIB_DIR/$script" "$CONFIG" single "$@"
    done
}

phase_compare() {
    for d in $DEPTHS; do
        local cram root
        cram="$(depth_cram "$d")"; root="$(depth_out_root "$d")"
        echo "=================== ${d}x : compare_to_truth ==================="
        SINGLE_CRAM="$cram" OUT_ROOT="$root" bash "$LIB_DIR/compare_to_truth.sh" "$CONFIG"
    done
}

# Merge per-depth comparison tables, prefixing a `depth` column. "full"
# is emitted as its numeric native mean so it sorts/plots after 100x.
phase_merge() {
    local qual_out="$SWEEP_ROOT/qual_dist_by_depth.tsv"
    local acc_out="$SWEEP_ROOT/accuracy_by_depth.tsv"
    mkdir -p "$SWEEP_ROOT"
    printf 'depth\tcaller\tclass\tstatus\tqual\tdp\n' > "$qual_out"
    printf 'depth\tcaller\tclass\ttp\tfp\tfn\tprecision\trecall\tf1\n' > "$acc_out"
    for d in $DEPTHS; do
        local root dlabel
        root="$(depth_out_root "$d")"
        if [[ "$d" == "full" ]]; then
            dlabel="$(awk -v m="$MEAN_DEPTH" 'BEGIN{printf "%.0f", m}')"
        else dlabel="$d"; fi
        local q="$root/comparison/qual_dist.tsv" a="$root/comparison/accuracy.tsv"
        [[ -s "$q" ]] && tail -n +2 "$q" | awk -v d="$dlabel" '{print d"\t"$0}' >> "$qual_out"
        [[ -s "$a" ]] && tail -n +2 "$a" | awk -v d="$dlabel" '{print d"\t"$0}' >> "$acc_out"
    done
    echo "merged qual : $qual_out ($(($(wc -l < "$qual_out")-1)) rows)"
    echo "merged acc  : $acc_out  ($(($(wc -l < "$acc_out")-1)) rows)"
}

case "$SUBCMD" in
    crams)     phase_crams ;;
    ours)      run_caller_phase run_ours.sh ;;
    freebayes) run_caller_phase run_freebayes.sh ;;
    gatk)      run_caller_phase run_gatk.sh ;;
    compare)   phase_compare ;;
    merge)     phase_merge ;;
    host)      phase_crams; run_caller_phase run_ours.sh; run_caller_phase run_freebayes.sh ;;
    *) echo "usage: $0 <config> {crams|ours|freebayes|gatk|compare|merge|host}" >&2; exit 2 ;;
esac
