#!/usr/bin/env bash
# Shared scaffolding for the benchmark caller runners.
#
# This file is *sourced* (not executed) by the run_<caller>.sh scripts.
# It holds everything those scripts used to duplicate: binary discovery,
# config loading, preflight checks, timing, record counting, and a
# dry-run-aware command wrapper. A per-benchmark `bench.config.sh`
# supplies the paths and knobs; see benchmarks/*/bench.config.sh.
#
# Contract a config file must satisfy (after bench_load_config):
#   REFERENCE          reference FASTA (.fai sibling required; .dict for GATK)
#   BED                regions BED, relative to BENCH_DIR or absolute
#   CRAM_DIR           dir holding the input CRAMs (default $BENCH_DIR/crams)
#   CRAM_GLOB          glob for cohort CRAMs within CRAM_DIR (e.g. "*.cram")
#   CRAM_SUFFIX        suffix stripped from a CRAM filename to form the
#                      sample/output base (e.g. ".bench.cram" or ".cram")
#   SINGLE_CRAM        basename of the CRAM used by `single` mode
#   PLOIDY             default 2
#   THREADS            default 4 (host P-core budget)
#   MIN_QUAL           freebayes post-call QUAL fairness cap (default 30)
#   GATK_HEAP          -Xmx for HaplotypeCaller / GenotypeGVCFs (default 4g)
#   GATK_COMBINE_HEAP  -Xmx for CombineGVCFs (default 8g)
#   OUT_ROOT           results dir (default $BENCH_DIR/results)
#   TRUTH_VCF          optional truth VCF for compare_to_truth.sh
#
# Globals exported for the runners: BENCH_DIR, PROJECT_ROOT, and every
# config var above (with defaults applied).
#
# Honored env knob: DRY_RUN=1 makes bench_run print the command it would
# run (shell-quoted) and skip execution — used to verify command parity
# without needing every tool installed.

set -euo pipefail

# benchmarks/lib/common.sh -> benchmarks/lib -> benchmarks -> repo root
LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$LIB_DIR/../.." && pwd)"

# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

# bench_load_config <path-to-bench.config.sh>
# Sources the config, derives BENCH_DIR from its location, resolves BED /
# CRAM_DIR / OUT_ROOT relative to BENCH_DIR, and applies defaults. Env
# vars already set win over config defaults because the config files use
# ${VAR:-default} themselves; here we only fill in what neither set.
bench_load_config() {
    local config="${1:-}"
    if [[ -z "$config" || ! -f "$config" ]]; then
        echo "usage: $0 <bench.config.sh> [single|cohort]" >&2
        echo "  config file not found: ${config:-<none>}" >&2
        exit 2
    fi
    BENCH_DIR="$(cd "$(dirname "$config")" && pwd)"
    # shellcheck disable=SC1090
    source "$config"

    : "${BENCH_NAME:=$(basename "$BENCH_DIR")}"
    : "${CRAM_DIR:=$BENCH_DIR/crams}"
    : "${CRAM_GLOB:=*.cram}"
    : "${CRAM_SUFFIX:=.cram}"
    : "${PLOIDY:=2}"
    : "${THREADS:=4}"
    : "${MIN_QUAL:=30}"
    : "${GATK_HEAP:=4g}"
    : "${GATK_COMBINE_HEAP:=8g}"
    : "${OUT_ROOT:=$BENCH_DIR/results}"

    # Resolve BED relative to BENCH_DIR if it isn't absolute.
    if [[ -n "${BED:-}" && "$BED" != /* ]]; then
        BED="$BENCH_DIR/$BED"
    fi
    if [[ -z "${REFERENCE:-}" ]]; then
        echo "config error: REFERENCE is unset (in $config)" >&2
        exit 2
    fi
}

# ---------------------------------------------------------------------------
# Binary discovery (identical logic the per-caller scripts used to repeat)
# ---------------------------------------------------------------------------

# Discover the pop_var_caller binary into POP_VAR_CALLER_BIN. Honors a
# pre-set POP_VAR_CALLER_BIN. Tries container build first (canonical on
# Linux per CLAUDE.md), then host build, verifying each actually runs on
# the current platform (a Linux ELF in target-container/ is +x but
# unusable on a macOS host).
bench_discover_ours_bin() {
    if [[ -z "${POP_VAR_CALLER_BIN:-}" ]]; then
        local candidate
        for candidate in \
            "$PROJECT_ROOT/target-container/release/pop_var_caller" \
            "$PROJECT_ROOT/target/release/pop_var_caller"; do
            if [[ -x "$candidate" ]] && "$candidate" --version >/dev/null 2>&1; then
                POP_VAR_CALLER_BIN="$candidate"
                break
            fi
        done
    fi
    if [[ -z "${POP_VAR_CALLER_BIN:-}" || ! -x "${POP_VAR_CALLER_BIN}" ]]; then
        echo "no pop_var_caller binary found." >&2
        echo "build with: ./scripts/dev.sh cargo build --release" >&2
        echo "or set POP_VAR_CALLER_BIN=<path>" >&2
        exit 1
    fi
}

# Discover GATK into GATK_BIN (default /opt/gatk/gatk — present in the
# dev container). Honors a pre-set GATK_BIN.
bench_discover_gatk_bin() {
    GATK_BIN="${GATK_BIN:-/opt/gatk/gatk}"
    [[ -x "$GATK_BIN" ]] || { echo "GATK binary not executable: $GATK_BIN" >&2; exit 1; }
}

# Discover freebayes into FREEBAYES_BIN. Honors a pre-set FREEBAYES_BIN.
bench_discover_freebayes_bin() {
    if [[ -z "${FREEBAYES_BIN:-}" ]] && command -v freebayes >/dev/null; then
        FREEBAYES_BIN="$(command -v freebayes)"
    fi
    if [[ -z "${FREEBAYES_BIN:-}" || ! -x "${FREEBAYES_BIN}" ]]; then
        echo "no freebayes binary found." >&2
        echo "install via package manager (apt/brew/conda) or set FREEBAYES_BIN=<path>" >&2
        exit 1
    fi
}

# The GATK sequence dictionary path for the reference. Convention:
# foo.fa / foo.fna -> foo.dict (strip the final extension).
bench_ref_dict() {
    printf '%s.dict\n' "${REFERENCE%.*}"
}

# ---------------------------------------------------------------------------
# CRAM discovery
# ---------------------------------------------------------------------------

# Populate the BENCH_CRAMS array with every CRAM matching CRAM_GLOB in
# CRAM_DIR (sorted). Exits if none found.
bench_list_crams() {
    shopt -s nullglob
    # shellcheck disable=SC2206
    BENCH_CRAMS=("$CRAM_DIR"/$CRAM_GLOB)
    shopt -u nullglob
    if (( ${#BENCH_CRAMS[@]} == 0 )); then
        echo "no CRAMs matching '$CRAM_GLOB' in $CRAM_DIR" >&2
        exit 1
    fi
}

# bench_sample_base <cram-path> -> base name with CRAM_SUFFIX stripped.
bench_sample_base() {
    basename "$1" "$CRAM_SUFFIX"
}

# Path to the single-mode CRAM (absolute). Resolves SINGLE_CRAM relative
# to CRAM_DIR if not absolute; defaults to the first cohort CRAM.
bench_single_cram() {
    if [[ -n "${SINGLE_CRAM:-}" ]]; then
        [[ "$SINGLE_CRAM" == /* ]] && { printf '%s\n' "$SINGLE_CRAM"; return; }
        printf '%s\n' "$CRAM_DIR/$SINGLE_CRAM"
        return
    fi
    bench_list_crams
    printf '%s\n' "${BENCH_CRAMS[0]}"
}

# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------

# bench_preflight <file>...  — exit 1 if any is missing.
bench_preflight() {
    local f
    for f in "$@"; do
        [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
    done
}

# ---------------------------------------------------------------------------
# Command execution (dry-run aware) + timing
# ---------------------------------------------------------------------------

# bench_run <logfile> -- <cmd> [args...]
# Tees the command's stderr to <logfile> (matching the old scripts'
# `2> >(tee LOG >&2)` idiom). With DRY_RUN=1 it prints the shell-quoted
# command and skips execution, so command-line parity can be verified
# without the tool installed. Pass "" as logfile to skip teeing.
bench_run() {
    local log="$1"; shift
    [[ "${1:-}" == "--" ]] && shift
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        printf 'DRY-RUN:'
        printf ' %q' "$@"
        printf '\n'
        return 0
    fi
    if [[ -n "$log" ]]; then
        "$@" 2> >(tee "$log" >&2)
    else
        "$@"
    fi
}

# Seconds since the epoch — wraps `date` so timing reads uniformly.
bench_now() { date +%s; }

# Count VCF data records (non-header lines). Works on plain or bgzipped.
bench_record_count() {
    local vcf="$1"
    if [[ "$vcf" == *.gz || "$vcf" == *.bgz ]]; then
        gzip -dc "$vcf" | grep -vc '^#' || true
    else
        grep -vc '^#' "$vcf" || true
    fi
}

# Number of sample columns in a VCF (#CHROM columns beyond the fixed 9).
bench_sample_count() {
    local vcf="$1" line
    if [[ "$vcf" == *.gz || "$vcf" == *.bgz ]]; then
        line="$(gzip -dc "$vcf" | grep -m1 '^#CHROM' || true)"
    else
        line="$(grep -m1 '^#CHROM' "$vcf" || true)"
    fi
    awk '{print NF-9}' <<<"$line"
}
