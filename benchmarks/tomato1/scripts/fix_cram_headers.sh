#!/usr/bin/env bash
# Fix the duplicate-PL @RG defect in the bench CRAM headers.
#
# Each upstream CRAM has an @RG line like
#   @RG  ID:...  SM:...  PL:ilumina  PL:PRJNA790656_<srr>
# which is invalid SAM (PL must appear at most once) and noodles refuses
# to parse it. htslib silently warns and keeps the first, which is why
# `samtools view` works but pop_var_caller and many strict readers do
# not.
#
# This script rewrites every @RG line in every bench CRAM in place:
#   - keep only the first PL tag
#   - fix the typo `ilumina` -> `ILLUMINA` so the value is in the SAM
#     spec's PL controlled vocabulary (GATK/freebayes accept it cleanly)
#   - the second `PL:` actually carries the library id (confirmed by
#     the original data producer: `PRJNA790656_<srr>` is the library);
#     rewrite it as `LB:` so the info is preserved under the correct
#     SAM tag instead of being thrown away
#   - leave all other tags and other header lines untouched
#
# Uses `samtools reheader --in-place`, which patches the CRAM header
# block without rewriting records and preserves the existing .crai
# index — so re-running is cheap and the slices stay byte-identical
# in their data portion.
#
# Idempotent: a re-run on already-fixed CRAMs is a no-op (the new
# header equals the old).
#
# Usage:
#   ./fix_cram_headers.sh           # fixes every *.bench.cram in ../crams
#   ./fix_cram_headers.sh -n        # dry-run: show what would change
#
# Env overrides:
#   REFERENCE    SL4.0 fasta (needed by samtools to decode the header)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CRAM_DIR="$TEST_DIR/crams"

REFERENCE="${REFERENCE:-/home/jose/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"

DRY_RUN=0
if [[ "${1:-}" == "-n" || "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=1
fi

command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "missing reference: $REFERENCE" >&2; exit 1; }

# Awk: emit a fixed copy of the input header.
#   - For @RG lines: process tags left-to-right, keep first PL only,
#     rewrite the typo `ilumina` -> `ILLUMINA`.
#   - Pass through every other line unchanged.
read -r -d '' FIX_AWK <<'AWK' || true
BEGIN { FS = "\t"; OFS = "\t" }
/^@RG/ {
    # First pass: detect whether an LB tag already exists. If yes
    # we must not invent a second one when rewriting the spare PL.
    has_lb = 0
    for (i = 2; i <= NF; i++) {
        if (substr($i, 1, 3) == "LB:") { has_lb = 1; break }
    }
    # Second pass: emit tags. Keep first PL (typo-fixed), rewrite
    # any further PL as LB (if free), otherwise drop.
    out = $1
    seen_pl = 0
    for (i = 2; i <= NF; i++) {
        if (substr($i, 1, 3) == "PL:") {
            v = substr($i, 4)
            if (!seen_pl) {
                seen_pl = 1
                if (tolower(v) == "ilumina") v = "ILLUMINA"
                out = out OFS "PL:" v
            } else if (!has_lb) {
                has_lb = 1   # claim the LB slot for this RG
                out = out OFS "LB:" v
            }
            # else: drop a second/third PL we can't relocate
        } else {
            out = out OFS $i
        }
    }
    print out
    next
}
{ print }
AWK

shopt -s nullglob
crams=("$CRAM_DIR"/*.bench.cram)
(( ${#crams[@]} > 0 )) || { echo "no *.bench.cram in $CRAM_DIR" >&2; exit 1; }

for cram in "${crams[@]}"; do
    base=$(basename "$cram" .bench.cram)

    old_header=$(mktemp -t "${base}.old.XXXXXX.sam")
    new_header=$(mktemp -t "${base}.new.XXXXXX.sam")
    trap 'rm -f "$old_header" "$new_header"' EXIT

    samtools view -H "$cram" --reference "$REFERENCE" > "$old_header"
    awk "$FIX_AWK" "$old_header" > "$new_header"

    if cmp -s "$old_header" "$new_header"; then
        echo "[skip] $base (header already clean)"
        rm -f "$old_header" "$new_header"
        continue
    fi

    if (( DRY_RUN )); then
        echo "[dry] $base would change:"
        diff -u "$old_header" "$new_header" | grep -E "^[+-]@RG" | head -4 | sed 's/^/    /'
        rm -f "$old_header" "$new_header"
        continue
    fi

    echo "[fix]  $base"
    samtools reheader --in-place "$new_header" "$cram"
    rm -f "$old_header" "$new_header"
done

# Sanity: re-check PL counts after the rewrite.
echo
echo "post-fix audit:"
for cram in "${crams[@]}"; do
    base=$(basename "$cram" .bench.cram)
    rg=$(samtools view -H "$cram" --reference "$REFERENCE" 2>/dev/null | grep -E "^@RG" | head -1)
    pl=$(echo "$rg" | grep -o $'\tPL:' | wc -l)
    printf "  %-22s PL_tags=%d  %s\n" "$base" "$pl" "$rg"
done
