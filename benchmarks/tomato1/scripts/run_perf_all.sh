#!/usr/bin/env bash
# Run the full tomato1 performance benchmark suite that feeds
# perf_dashboard.py (the four-section pop_var_caller vs freebayes vs GATK
# comparison).
#
# Split by where each tool lives:
#   - HOST  (uv run --script): pop_var_caller + freebayes.
#   - CONTAINER (scripts/dev.sh, reference genome bind-mounted): GATK,
#     which isn't on the host. The genome lives outside the project tree,
#     so it's mounted via DEV_EXTRA_MOUNT.
#
# Sections and the scripts that feed them:
#   §1 single-sample direct (CRAM→VCF): perf_ours_from_bam, perf_freebayes,
#                                        perf_gatk_direct   (all N=1 rows)
#   §2 one intermediate, 4 threads:      perf_ours_psp_4t, perf_gatk_gvcf_4t
#   §3 scaling CRAM→VCF:                 perf_ours_whole_pipeline,
#                                        perf_freebayes, perf_gatk_direct
#   §4 scaling intermediate→VCF:         perf_ours_joint, perf_gatk_joint
#
# §4 assumes the per-sample intermediates already exist:
#   - ours: results/ours/cohort/psp/*.psp
#   - GATK: results/gatk/cohort/gvcf/*.g.vcf.gz
# build them first with the cohort drivers if missing (see the warnings
# this script prints).
#
# Usage:
#   benchmarks/tomato1/scripts/run_perf_all.sh [all|host|container]
#
#   all        run everything (default)
#   host       only the host scripts (pop_var_caller + freebayes)
#   container  only the GATK scripts (inside scripts/dev.sh)
#
# Env:
#   GENOMES            dir bind-mounted into the container for GATK's
#                      reference (default: $HOME/genomes)
#   POP_VAR_CALLER_BIN pop_var_caller release binary (default:
#                      <project>/target/release/pop_var_caller)
#   SIZES, THREADS, REFERENCE, GATK_BIN, JAVA_HEAP, ... forwarded to the
#                      per-caller scripts (see each script's header).

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
GENOMES="${GENOMES:-$HOME/genomes}"
MODE="${1:-all}"

case "$MODE" in
    all|host|container) ;;
    *) echo "usage: $0 [all|host|container]" >&2; exit 2 ;;
esac

# Host scripts in dashboard order; container (GATK) scripts separately.
HOST_SCRIPTS=(
    perf_ours_from_bam.py        # §1
    perf_ours_psp_4t.py          # §2
    perf_ours_whole_pipeline.py  # §3
    perf_ours_joint.py           # §4
    perf_freebayes.py            # §1 + §3
)
CONTAINER_SCRIPTS=(
    perf_gatk_direct.py          # §1 + §3
    perf_gatk_gvcf_4t.py         # §2
    perf_gatk_joint.py           # §4
)

declare -a RESULTS=()

record() { RESULTS+=("$1"); }

run_host() {
    local script="$1"
    echo
    echo "=== [host] $script ==="
    if uv run --script "$SCRIPT_DIR/$script"; then
        record "ok   host/$script"
    else
        record "FAIL host/$script (exit $?)"
    fi
}

run_container() {
    local script="$1"
    echo
    echo "=== [container] $script ==="
    if DEV_EXTRA_MOUNT="$GENOMES" "$PROJECT_ROOT/scripts/dev.sh" \
            python3 "benchmarks/tomato1/scripts/$script"; then
        record "ok   container/$script"
    else
        record "FAIL container/$script (exit $?)"
    fi
}

# --- Pre-flight warnings for §4 (intermediates must be pre-built) -------
psp_dir="$TEST_DIR/results/ours/cohort/psp"
gvcf_dir="$TEST_DIR/results/gatk/cohort/gvcf"
if [[ "$MODE" != container ]] && ! compgen -G "$psp_dir/*.psp" >/dev/null; then
    echo "WARNING: no *.psp in $psp_dir — §4 perf_ours_joint will have no inputs." >&2
    echo "         build them first (cohort pileup), or §4-ours will be skipped/empty." >&2
fi
if [[ "$MODE" != host ]] && ! compgen -G "$gvcf_dir/*.g.vcf.gz" >/dev/null; then
    echo "WARNING: no *.g.vcf.gz in $gvcf_dir — §4 perf_gatk_joint will have no inputs." >&2
    echo "         build them first (cohort HaplotypeCaller), or §4-GATK will be skipped/empty." >&2
fi

# --- Run -----------------------------------------------------------------
if [[ "$MODE" == all || "$MODE" == host ]]; then
    for s in "${HOST_SCRIPTS[@]}"; do run_host "$s"; done
fi
if [[ "$MODE" == all || "$MODE" == container ]]; then
    for s in "${CONTAINER_SCRIPTS[@]}"; do run_container "$s"; done
fi

# --- Summary -------------------------------------------------------------
echo
echo "=== summary ==="
printf '%s\n' "${RESULTS[@]}"
echo
echo "view the dashboard with:"
echo "  uvx marimo edit --sandbox benchmarks/tomato1/scripts/perf_dashboard.py"

# Non-zero overall exit if anything failed.
if printf '%s\n' "${RESULTS[@]}" | grep -q '^FAIL'; then
    exit 1
fi
