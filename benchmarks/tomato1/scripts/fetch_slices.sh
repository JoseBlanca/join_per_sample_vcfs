#!/usr/bin/env bash
# Pull the sliced CRAMs from rick into this project's tmp/ tree.
#
# Run this ON WULF (from the project root). Uses rsync over ssh.
#
# Usage:
#   ./scripts/benchmark_cohort/fetch_slices.sh [REMOTE] [REMOTE_DIR] [LOCAL_DIR]
#
# Defaults:
#   REMOTE      = joxi@rick
#   REMOTE_DIR  = /home/joxi/tmp/benchmark_cohort/crams
#   LOCAL_DIR   = tmp/benchmark_cohort/crams
#
# Also runs `samtools quickcheck` on the fetched files as a sanity test.

set -euo pipefail

REMOTE=${1:-joxi@rick}
REMOTE_DIR=${2:-/home/joxi/tmp/benchmark_cohort/crams}
LOCAL_DIR=${3:-tmp/benchmark_cohort/crams}

mkdir -p "$LOCAL_DIR"

echo "rsync $REMOTE:$REMOTE_DIR/ -> $LOCAL_DIR/"
rsync -avh --progress \
    --include='*.bench.cram' \
    --include='*.bench.cram.crai' \
    --exclude='*' \
    "$REMOTE:$REMOTE_DIR/" "$LOCAL_DIR/"

echo
echo "verifying with samtools quickcheck"
if command -v samtools >/dev/null; then
    if samtools quickcheck -v "$LOCAL_DIR"/*.bench.cram; then
        echo "all slices OK"
    else
        echo "WARNING: quickcheck reported issues above" >&2
    fi
else
    echo "samtools not on PATH; skipping quickcheck"
fi

echo
echo "fetched:"
ls -lh "$LOCAL_DIR"/*.bench.cram 2>/dev/null
