#!/usr/bin/env bash
# Single-sample test: CRAM -> VCF via the direct path
# (`pop_var_caller var-calling-from-bam`, no .psp intermediate).
#
# Sample: SRR17274057.p1.bench.cram
# Output: tmp/tomato_cohort_test/results/ours/single_SRR17274057.vcf
#
# Env overrides:
#   POP_VAR_CALLER_BIN  binary to invoke (default: auto-detect target/release
#                       then target-container/release, relative to project root)
#   THREADS             rayon worker count (default: 4 — host P-core budget)
#   REFERENCE           tomato SL4.0 fasta (.fai sibling required)
#                       default: $HOME/genomes/s_lycopersicum/4.00/
#                                S_lycopersicum_chromosomes.4.00.fa
#   EXTRA_ARGS          appended to the command line (e.g. --no-complexity-filter)

set -euo pipefail

# --- locate project root + this test dir ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$TEST_DIR/../.." && pwd)"

# --- inputs ---
SAMPLE="SRR17274057"
CRAM="$TEST_DIR/crams/${SAMPLE}.p1.bench.cram"
REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"

# --- binary discovery ---
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

THREADS="${THREADS:-4}"
EXTRA_ARGS=${EXTRA_ARGS:-}

# --- preflight ---
for f in "$CRAM" "${CRAM}.crai" "$REFERENCE" "${REFERENCE}.fai"; do
    [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done

# --- output ---
OUT_DIR="$TEST_DIR/results/ours"
OUT_VCF="$OUT_DIR/single_${SAMPLE}.vcf"
LOG="$OUT_DIR/single_${SAMPLE}.log"
mkdir -p "$OUT_DIR"

echo "binary    : $POP_VAR_CALLER_BIN"
echo "sample    : $SAMPLE"
echo "input     : $CRAM"
echo "reference : $REFERENCE"
echo "threads   : $THREADS"
echo "output    : $OUT_VCF"
echo "log       : $LOG"
echo

t0=$(date +%s)
"$POP_VAR_CALLER_BIN" var-calling-from-bam \
    --reference "$REFERENCE" \
    --output "$OUT_VCF" \
    --threads "$THREADS" \
    $EXTRA_ARGS \
    "$CRAM" \
    2> >(tee "$LOG" >&2)
t1=$(date +%s)

echo
echo "elapsed: $((t1 - t0)) s"
echo "records: $(grep -vc '^#' "$OUT_VCF") (PASS+non-PASS, single sample)"
