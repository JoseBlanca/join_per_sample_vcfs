#!/usr/bin/env bash
# Cohort test: all 18 CRAMs -> per-sample .psp files -> joined multi-sample VCF.
#
# Stage 1 (per sample): `pop_var_caller pileup` over each *.bench.cram.
# Stage 2 (cohort):     `pop_var_caller var-calling` over all .psp files.
#
# Pileups run sequentially (each at THREADS workers) to keep within the
# host P-core budget (memory: feedback_host_p_core_budget). Already-built
# .psp outputs are skipped on re-runs.
#
# Output layout:
#   tmp/tomato_cohort_test/results/ours/cohort/psp/<sample>.psp
#   tmp/tomato_cohort_test/results/ours/cohort/cohort.vcf
#   tmp/tomato_cohort_test/results/ours/cohort/*.log
#
# Env overrides:
#   POP_VAR_CALLER_BIN  binary to invoke (default: auto-detect)
#   THREADS             rayon worker count per stage (default: 4)
#   REFERENCE           tomato SL4.0 fasta (.fai sibling required)
#   PILEUP_EXTRA        appended to each `pileup` invocation
#   VARCALL_EXTRA       appended to the `var-calling` invocation
#                       (e.g. --no-complexity-filter, --contamination-estimates ...)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$TEST_DIR/../.." && pwd)"

REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"
THREADS="${THREADS:-4}"
PILEUP_EXTRA=${PILEUP_EXTRA:-}
VARCALL_EXTRA=${VARCALL_EXTRA:-}

if [[ -z "${POP_VAR_CALLER_BIN:-}" ]]; then
    # Try candidates in order; verify each actually runs on the
    # current platform (the executable bit is no guarantee — e.g. a
    # Linux ELF in target-container/ is +x but unusable on macOS host).
    # Order: container build first (canonical workflow per CLAUDE.md
    # on Linux), then host build (only choice on macOS host).
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

CRAM_DIR="$TEST_DIR/crams"
OUT_DIR="$TEST_DIR/results/ours/cohort"
PSP_DIR="$OUT_DIR/psp"
COHORT_VCF="$OUT_DIR/cohort.vcf"
COHORT_LOG="$OUT_DIR/cohort.log"
mkdir -p "$PSP_DIR"

shopt -s nullglob
crams=("$CRAM_DIR"/*.bench.cram)
if (( ${#crams[@]} == 0 )); then
    echo "no *.bench.cram in $CRAM_DIR" >&2
    exit 1
fi

for f in "$REFERENCE" "${REFERENCE}.fai"; do
    [[ -f "$f" ]] || { echo "missing reference file: $f" >&2; exit 1; }
done

echo "binary    : $POP_VAR_CALLER_BIN"
echo "samples   : ${#crams[@]}"
echo "reference : $REFERENCE"
echo "threads   : $THREADS (per pileup, and for cohort var-calling)"
echo "psp dir   : $PSP_DIR"
echo "vcf out   : $COHORT_VCF"
echo

# ---------- Stage 1: per-sample pileup (sequential) ----------
t_stage1=$(date +%s)
for cram in "${crams[@]}"; do
    base=$(basename "$cram" .bench.cram)  # strips .bench.cram -> "SRR....p1"
    psp="$PSP_DIR/${base}.psp"
    log="$PSP_DIR/${base}.log"
    if [[ -s "$psp" ]]; then
        echo "[skip pileup] $base (.psp already present)"
        continue
    fi
    echo "[pileup] $base"
    "$POP_VAR_CALLER_BIN" pileup \
        --reference "$REFERENCE" \
        --output "$psp" \
        --threads "$THREADS" \
        $PILEUP_EXTRA \
        "$cram" \
        2> >(tee "$log" >&2)
done
t_stage1_end=$(date +%s)
echo
echo "stage 1 elapsed: $((t_stage1_end - t_stage1)) s"
echo

# ---------- Stage 2: cohort var-calling ----------
psps=("$PSP_DIR"/*.psp)
if (( ${#psps[@]} != ${#crams[@]} )); then
    echo "WARNING: expected ${#crams[@]} .psp files, found ${#psps[@]}" >&2
fi

echo "[var-calling] ${#psps[@]} samples -> $COHORT_VCF"
t_stage2=$(date +%s)
"$POP_VAR_CALLER_BIN" var-calling \
    --reference "$REFERENCE" \
    --output "$COHORT_VCF" \
    --threads "$THREADS" \
    $VARCALL_EXTRA \
    "${psps[@]}" \
    2> >(tee "$COHORT_LOG" >&2)
t_stage2_end=$(date +%s)

echo
echo "stage 2 elapsed: $((t_stage2_end - t_stage2)) s"
echo "total   elapsed: $((t_stage2_end - t_stage1)) s"
echo "records: $(grep -vc '^#' "$COHORT_VCF") (cohort VCF, all sites)"
