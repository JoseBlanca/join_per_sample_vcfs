#!/usr/bin/env bash
# Run freebayes over the GIAB `per_sample` dataset, for a freebayes-vs-ours
# accuracy comparison.
#
#   benchmarks/giab/src/run_freebayes_per_sample.sh [COVERAGE] [SAMPLE ...]
#
# Mirrors run_ours_per_sample.sh: the per_sample dataset holds three
# independent single-sample callsets (HG002, HG003, HG004), each with its
# OWN random 100-region BED and its OWN GIAB truth VCF, so every sample is
# called on its own, restricted to its confident regions.
#
# NOTE: freebayes emits records in --targets BED order, so the per_sample BEDs
# must be in reference contig order or the output VCF won't be coordinate-sorted
# (and `bcftools index -t` will fail). Run benchmarks/giab/src/sort_beds.sh once
# to normalize them; it's idempotent.
#
# freebayes is run PER SAMPLE, PER COVERAGE, restricted to that sample's BED
# via --targets, with --fasta-reference and --ploidy 2. Unlike our pileup,
# freebayes does NOT require an @SQ M5 tag, so it reads the original 300x
# alignments directly (it happily reads the M5-reheadered copies too — we
# reuse run_ours_per_sample.sh's file naming for simplicity).
#
# Fairness: freebayes emits a long low-QUAL tail the other callers don't.
# A QUAL >= MIN_QUAL (default 30, matching our caller's --min-qual 30)
# post-filter is applied so precision is comparable. We write BOTH:
#   {sample}.raw.vcf  — every freebayes record (no QUAL filter)
#   {sample}.vcf      — QUAL >= MIN_QUAL (the headline, ours-comparable set)
#
# Args:
#   COVERAGE   bam/ coverage subdir to use (default: 300x)
#   SAMPLE...  subset of {HG002,HG003,HG004} to run (default: all three)
#
# Env overrides:
#   FREEBAYES_BIN  binary to invoke (default: auto-detect `freebayes` on PATH)
#   PLOIDY         --ploidy (default 2)
#   MIN_QUAL       QUAL floor for the gated output (default 30)
#   DRY_RUN=1      print the commands instead of running them

set -euo pipefail

# benchmarks/giab/src -> benchmarks/giab -> benchmarks -> repo root
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SRC_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$BENCH_DIR/../.." && pwd)"

COVERAGE="${1:-300x}"
shift || true
SAMPLES=("$@")
if (( ${#SAMPLES[@]} == 0 )); then
    SAMPLES=(HG002 HG003 HG004)
fi

PLOIDY="${PLOIDY:-2}"
MIN_QUAL="${MIN_QUAL:-30}"
REFERENCE="$BENCH_DIR/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BAM_DIR="$BENCH_DIR/per_sample/bam/$COVERAGE"
BED_DIR="$BENCH_DIR/per_sample/bed"
OUT_DIR="$BENCH_DIR/results/per_sample/$COVERAGE/freebayes"

# --- sample -> alignment filename / BED filename ---------------------------
# Same two naming schemes as run_ours_per_sample.sh (kept in sync):
#   300x  — HG002 is a CRAM; HG003/HG004 are the `.m5.bam` reheadered copies.
#           (freebayes doesn't need M5, but those are the files we keep around;
#           the originals would work just as well if present.)
#   else  — downsampled tiers (5x..50x): `{sample}.{cov}.seed42.bam`.
align_file() {
    local dir="$BAM_DIR"
    if [[ "$COVERAGE" != "300x" ]]; then
        local bam="$dir/${1}.${COVERAGE}.seed42.bam"
        if [[ -f "$bam" ]]; then echo "$bam"; return 0; fi
        echo "missing alignment: $bam" >&2; return 1
    fi
    case "$1" in
        HG002) echo "$dir/HG002_reads_selected_100_rg.cram" ;;
        HG003|HG004)
            local m5="$dir/${1}_bench_azar_merged_100.sorted.m5.bam"
            local orig="$dir/${1}_bench_azar_merged_100.sorted.bam"
            if [[ -f "$m5" ]]; then echo "$m5"
            elif [[ -f "$orig" ]]; then echo "$orig"
            else
                echo "missing 300x BAM for $1: $m5 (or $orig)" >&2
                return 1
            fi ;;
        *) echo "unknown sample: $1" >&2; return 1 ;;
    esac
}
bed_file() {
    case "$1" in
        HG002) echo "HG002_bench_azar_merged_100.bed" ;;
        HG003) echo "HG003_bench_azar_merged_100.bed" ;;
        HG004) echo "HG004_bench_azar_merged_100.bed" ;;
        *) echo "unknown sample: $1" >&2; return 1 ;;
    esac
}

# --- binary discovery ------------------------------------------------------
discover_bin() {
    if [[ -z "${FREEBAYES_BIN:-}" ]]; then
        local candidate
        for candidate in /opt/homebrew/bin/freebayes "$(command -v freebayes 2>/dev/null || true)"; do
            if [[ -n "$candidate" && -x "$candidate" ]]; then
                FREEBAYES_BIN="$candidate"
                break
            fi
        done
    fi
    if [[ -z "${FREEBAYES_BIN:-}" || ! -x "${FREEBAYES_BIN}" ]]; then
        echo "no freebayes binary found (set FREEBAYES_BIN=<path>)" >&2
        exit 1
    fi
}

preflight() {
    local f
    for f in "$@"; do
        [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
    done
}

run() {
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        printf 'DRY-RUN:'; printf ' %q' "$@"; printf '\n'
        return 0
    fi
    "$@"
}

record_count() {
    local vcf="$1"
    [[ -f "$vcf" ]] || { echo "?"; return; }
    grep -vc '^#' "$vcf" || true
}

# header lines pass through; data rows kept iff QUAL (col 6) is numeric and
# >= min. `$6 + 0` coerces; non-numeric ('.') falls to 0 and is dropped.
# awk avoids a bcftools dependency in the runner itself.
FILTER_AWK='/^#/ { print; next } $6 != "." && ($6 + 0) >= min { print }'

discover_bin
preflight "$REFERENCE" "${REFERENCE}.fai"
mkdir -p "$OUT_DIR"

FB_VERSION="$("$FREEBAYES_BIN" --version 2>&1 | head -1 || true)"
echo "binary    : $FREEBAYES_BIN ($FB_VERSION)"
echo "dataset   : per_sample / $COVERAGE"
echo "reference : $REFERENCE"
echo "ploidy    : $PLOIDY"
echo "min QUAL  : $MIN_QUAL (gated output; raw kept alongside)"
echo "samples   : ${SAMPLES[*]}"
echo "out dir   : $OUT_DIR"
echo

for sample in "${SAMPLES[@]}"; do
    aln="$(align_file "$sample")"
    bed="$BED_DIR/$(bed_file "$sample")"
    raw_vcf="$OUT_DIR/${sample}.raw.vcf"
    vcf="$OUT_DIR/${sample}.vcf"
    log="$OUT_DIR/${sample}.log"

    preflight "$aln" "$bed"

    echo "=================================================================="
    echo "sample    : $sample"
    echo "alignment : $aln"
    echo "regions   : $bed ($(wc -l < "$bed") intervals)"
    echo "raw vcf   : $raw_vcf"
    echo "vcf       : $vcf (QUAL >= $MIN_QUAL)"
    echo

    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        printf 'DRY-RUN:'
        printf ' %q' "$FREEBAYES_BIN" --fasta-reference "$REFERENCE" \
            --targets "$bed" --ploidy "$PLOIDY" "$aln"
        printf ' > %q\n' "$raw_vcf"
        printf 'DRY-RUN: awk -v min=%q %q %q > %q\n' "$MIN_QUAL" "$FILTER_AWK" "$raw_vcf" "$vcf"
        echo
        continue
    fi

    t0=$(date +%s)
    echo "[freebayes] $sample -> $raw_vcf"
    "$FREEBAYES_BIN" \
        --fasta-reference "$REFERENCE" \
        --targets "$bed" \
        --ploidy "$PLOIDY" \
        "$aln" \
        > "$raw_vcf" 2> >(tee "$log" >&2)
    t1=$(date +%s)

    # Apply the QUAL gate to produce the ours-comparable headline set.
    awk -v min="$MIN_QUAL" "$FILTER_AWK" "$raw_vcf" > "$vcf"

    echo
    echo "elapsed               : $((t1 - t0)) s"
    echo "records (raw)         : $(record_count "$raw_vcf")"
    echo "records (QUAL >= $MIN_QUAL) : $(record_count "$vcf")"
    echo
done

echo "done. VCFs under: $OUT_DIR"
