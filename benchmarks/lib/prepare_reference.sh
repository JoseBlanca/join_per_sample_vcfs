#!/usr/bin/env bash
# Ensure a benchmark's reference FASTA has the index siblings the callers
# need: a `.fai` (samtools faidx — every caller) and a `.dict` (Picard
# sequence dictionary — GATK only). Idempotent: existing, non-empty
# siblings are left alone unless FORCE=1.
#
#   benchmarks/lib/prepare_reference.sh <bench.config.sh>
#
# Uses `samtools` for both, so it needs no GATK/Picard on the host.
# The .dict convention matches the runners: foo.fna -> foo.dict.
#
# Env: FORCE=1 rebuild even if siblings exist.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

bench_load_config "${1:-}"

command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "reference FASTA missing: $REFERENCE" >&2; exit 1; }

FAI="${REFERENCE}.fai"
DICT="$(bench_ref_dict)"
FORCE="${FORCE:-0}"

echo "reference : $REFERENCE"
echo "fai       : $FAI"
echo "dict      : $DICT"
echo

if [[ "$FORCE" != "1" && -s "$FAI" ]]; then
    echo "[skip] .fai present"
else
    echo "[faidx] building $FAI"
    samtools faidx "$REFERENCE"
fi

if [[ "$FORCE" != "1" && -s "$DICT" ]]; then
    echo "[skip] .dict present"
else
    echo "[dict] building $DICT"
    # -o writes to a specific path; samtools dict won't overwrite via
    # redirection ambiguity this way.
    samtools dict "$REFERENCE" -o "$DICT"
fi

echo
echo "done."
