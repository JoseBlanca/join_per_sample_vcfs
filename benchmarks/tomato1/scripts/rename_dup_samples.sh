#!/usr/bin/env bash
# Make every bench CRAM's sample name (@RG SM:) unique.
#
# The tomato pool has biosamples (SRS…) sequenced in more than one run
# (SRR…); each run's CRAM carries the SAME SM, so a cohort VCF — which is
# keyed by sample name — sees duplicate columns (our caller errors;
# freebayes/GATK silently merge same-SM reads). This rewrites the SM on
# the 2nd+ run of each duplicated biosample to "<SM>_<SRR>", keeping the
# earliest run's SM unchanged. Result: one unique sample name per CRAM,
# still traceable to both biosample and run.
#
# Uses `samtools reheader --in-place` (patches the header block, keeps the
# .crai). Idempotent: once renamed the SMs are unique, so a re-run is a
# no-op. Run AFTER fix_cram_headers.sh (needs single, clean @RG SM:).
#
# Usage:   benchmarks/tomato1/scripts/rename_dup_samples.sh [-n]
#   -n / --dry-run : show the renames without applying them.
#
# Env:  REFERENCE  SL4.0 fasta (needed to read the CRAM header)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CRAM_DIR="$TEST_DIR/crams"
REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"
SUFFIX=".bench.cram"

DRY=0
[[ "${1:-}" == "-n" || "${1:-}" == "--dry-run" ]] && DRY=1

command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "missing reference: $REFERENCE" >&2; exit 1; }

sm_of() {  # first @RG SM: value of a CRAM
    samtools view -H "$1" --reference "$REFERENCE" 2>/dev/null \
        | grep -m1 '^@RG' | grep -o 'SM:[^[:space:]]*' | sed 's/SM://'
}

shopt -s nullglob
crams=("$CRAM_DIR"/*"$SUFFIX")
shopt -u nullglob
(( ${#crams[@]} > 0 )) || { echo "no *$SUFFIX in $CRAM_DIR" >&2; exit 1; }
IFS=$'\n' crams=($(printf '%s\n' "${crams[@]}" | sort)); unset IFS

# bash 3.2 (macOS) has no associative arrays; track seen SMs in a string
# (SM values contain no spaces, so a space-delimited membership test is safe).
seen=" "
renamed=0
for cram in "${crams[@]}"; do
    base="$(basename "$cram" "$SUFFIX")"     # e.g. SRR7279484.p1
    run="${base%%.*}"                         # e.g. SRR7279484
    sm="$(sm_of "$cram")"
    [[ -z "$sm" ]] && { echo "WARN: no SM in $base, skipping" >&2; continue; }
    if [[ "$seen" != *" $sm "* ]]; then
        seen="$seen$sm "
        continue                              # first run of this biosample: keep
    fi
    new="${sm}_${run}"
    echo "[rename] $base : SM:$sm -> SM:$new"
    (( DRY )) && continue
    hdr="$(mktemp -t "${base}.hdr.XXXXXX.sam")"
    samtools view -H "$cram" --reference "$REFERENCE" \
        | awk -v old="$sm" -v new="$new" 'BEGIN{FS=OFS="\t"}
            /^@RG/{ for(i=1;i<=NF;i++) if($i=="SM:" old) $i="SM:" new } {print}' > "$hdr"
    samtools reheader --in-place "$hdr" "$cram"
    rm -f "$hdr"
    renamed=$((renamed+1))
done

echo
if (( DRY )); then
    echo "dry-run: would rename the SMs listed above."
else
    echo "renamed $renamed CRAM(s)."
    echo "unique SMs now: $(for c in "${crams[@]}"; do sm_of "$c"; done | sort -u | wc -l | tr -d ' ') / ${#crams[@]}"
fi
