#!/usr/bin/env bash
# Build a parallel cohort PSP set using --block-target-bytes 262144
# (256 KiB, 1/4 of the 1 MiB pileup default), so perf_ours_joint_smallblock
# has inputs without disturbing the default-block PSPs that
# perf_ours_joint reads.
#
# Idempotent: skips per-sample PSPs that already exist. Delete the
# target file to force a rebuild.
#
# Env overrides:
#   POP_VAR_CALLER_BIN  binary (default: auto-detect target/release)
#   REFERENCE           SL4.0 fasta (.fai sibling required)
#   THREADS             pileup --threads (default: 4)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$TEST_DIR/../.." && pwd)"

REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"
THREADS="${THREADS:-4}"
BLOCK_TARGET=262144  # 256 KiB

if [[ -z "${POP_VAR_CALLER_BIN:-}" ]]; then
    for candidate in \
        "$PROJECT_ROOT/target-container/release/pop_var_caller" \
        "$PROJECT_ROOT/target/release/pop_var_caller"; do
        if [[ -x "$candidate" ]] && "$candidate" --version >/dev/null 2>&1; then
            POP_VAR_CALLER_BIN="$candidate"
            break
        fi
    done
fi
[[ -n "${POP_VAR_CALLER_BIN:-}" && -x "${POP_VAR_CALLER_BIN}" ]] \
    || { echo "no pop_var_caller binary found" >&2; exit 1; }

CRAM_DIR="$TEST_DIR/crams"
OUT_DIR="$TEST_DIR/results/ours/cohort/psp_smallblock"
mkdir -p "$OUT_DIR"

shopt -s nullglob
crams=("$CRAM_DIR"/*.bench.cram)
(( ${#crams[@]} > 0 )) || { echo "no *.bench.cram in $CRAM_DIR" >&2; exit 1; }

for f in "$REFERENCE" "${REFERENCE}.fai"; do
    [[ -f "$f" ]] || { echo "missing reference file: $f" >&2; exit 1; }
done

echo "binary       : $POP_VAR_CALLER_BIN"
echo "samples      : ${#crams[@]}"
echo "block target : $BLOCK_TARGET bytes (256 KiB)"
echo "out dir      : $OUT_DIR"
echo

for cram in "${crams[@]}"; do
    base=$(basename "$cram" .bench.cram)
    psp="$OUT_DIR/${base}.psp"
    if [[ -s "$psp" ]]; then
        echo "[skip] $base"
        continue
    fi
    echo "[pileup] $base"
    "$POP_VAR_CALLER_BIN" pileup \
        --reference "$REFERENCE" \
        --output "$psp" \
        --threads "$THREADS" \
        --block-target-bytes "$BLOCK_TARGET" \
        "$cram"
done

echo
echo "done. $(ls "$OUT_DIR"/*.psp | wc -l | tr -d ' ') PSPs in $OUT_DIR"
