#!/usr/bin/env bash
# Step 1 of option (c): per-sample mean read coverage in fixed W-bp windows,
# taken from the .psp files (NO BAM access). Depth = column 4 of psp-to-pileup
# (post-filter, deduplicated, mate-overlap-resolved fragment coverage).
#
# Output per sample: results/window_cov/<sample>.w<W>.tsv
#   chrom  win_start(1-based)  n_covered_positions  sum_depth
# (mean depth = sum_depth / n_covered_positions; breadth = n_covered / W)
#
# Usage:  W=2000 JOBS=10 benchmarks/tomato2/src/build_window_coverage.sh
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$TEST_DIR/../.." && pwd)"

W="${W:-2000}"
JOBS="${JOBS:-10}"
PSP_DIR="$TEST_DIR/psp_files"
OUT_DIR="$TEST_DIR/results/window_cov"
BIN="${POP_VAR_CALLER_BIN:-$PROJECT_ROOT/target/release/pop_var_caller}"
mkdir -p "$OUT_DIR"

echo "=== window coverage: W=${W}bp, JOBS=$JOBS -> $OUT_DIR ==="

process() {
    local psp="$1" base out
    base="$(basename "$psp" .psp)"
    out="$OUT_DIR/$base.w$W.tsv"
    if [[ -s "$out" ]]; then echo "[skip] $base"; return 0; fi
    if "$BIN" psp-to-pileup --input "$psp" 2>/dev/null \
        | awk -v W="$W" 'BEGIN{OFS="\t"}
            { w=int(($2-1)/W)*W+1; k=$1 SUBSEP w; s[k]+=$4; n[k]++ }
            END{ for(k in s){ split(k,a,SUBSEP); print a[1],a[2],n[k],s[k] } }' \
        | sort -k1,1 -k2,2n > "$out"; then
        echo "[done] $base ($(wc -l < "$out") windows)"
    else
        echo "[FAIL] $base"; rm -f "$out"
    fi
}
export -f process
export W OUT_DIR BIN

ls "$PSP_DIR"/*.psp | xargs -P "$JOBS" -I{} bash -c 'process "$@"' _ {}
echo "done: $(ls "$OUT_DIR"/*.w$W.tsv 2>/dev/null | wc -l | tr -d ' ') sample files"
