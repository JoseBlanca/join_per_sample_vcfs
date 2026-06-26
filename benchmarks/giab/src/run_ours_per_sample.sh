#!/usr/bin/env bash
# Run pop_var_caller (our SNP caller) over the GIAB `per_sample` dataset.
#
#   benchmarks/giab/src/run_ours_per_sample.sh [COVERAGE] [SAMPLE ...]
#
# The per_sample dataset holds three independent single-sample callsets
# (HG002, HG003, HG004). Each sample has its OWN random 100-region BED and
# its OWN GIAB truth VCF — the regions are NOT shared across samples — so
# every sample is called on its own, restricted to its confident regions.
#
# For each sample this runs the two-stage path at 6 threads:
#   1. pileup      : alignment (CRAM/BAM) --regions BED -> .psp
#   2. var-calling : .psp           --regions BED -> single-sample VCF
#
# Restricting BOTH stages to the sample's BED keeps the calls inside the
# GIAB high-confidence regions, which is what the FP/TP/FN/TN evaluation
# (a separate step) scores against the truth VCF.
#
# Args:
#   COVERAGE   bam/ coverage subdir to use (default: 300x)
#   SAMPLE...  subset of {HG002,HG003,HG004} to run (default: all three)
#
# Env overrides:
#   POP_VAR_CALLER_BIN  binary to invoke (default: auto-detect host/container build)
#   THREADS             worker threads for both stages (default: 6)
#   DRY_RUN=1           print the commands instead of running them

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

THREADS="${THREADS:-6}"
REFERENCE="$BENCH_DIR/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BAM_DIR="$BENCH_DIR/per_sample/bam/$COVERAGE"
BED_DIR="$BENCH_DIR/per_sample/bed"

# Feature toggles. Our caller enables BAQ (pileup) and the DUST low-complexity
# filter (var-calling) by default; both suppress true indels, so the GIAB
# indel evaluation runs with NO_BAQ=1 NO_DUST=1. The output dir is tagged with
# the active config so on/off runs don't clobber each other.
NO_BAQ="${NO_BAQ:-0}"
NO_DUST="${NO_DUST:-0}"
# MIN_MAPQ: pileup per-read MAPQ floor (caller default 20). Lower it (e.g. 0)
# to admit ambiguously-mapped reads so the var-calling MAPQ-diff filter, not
# the pileup, decides their fate. Tagged into the variant dir when non-default.
MIN_MAPQ="${MIN_MAPQ:-}"
# NO_AB: set to 1 to disable the var-calling allele-balance filter
# (--no-allele-balance-filter). The filter is ON by default (drops biallelic
# SNP het calls whose ALT read fraction is inconsistent with ~0.5; SNP-only,
# depth-aware so inert below ~30x). AB_LR overrides the threshold (default -5).
NO_AB="${NO_AB:-0}"
AB_LR="${AB_LR:-}"

# PRESET: named bundles that set the toggles above together (see GIAB README).
# The allele-balance filter is ON in BOTH presets — it removes ~94% of
# high-coverage SNP false positives for <1% true-call loss, so it is a near-free
# precision win worth keeping even at maximum sensitivity. The presets differ
# only on BAQ/DUST (the indel-recall vs SNP-strictness axis):
#   high-recall     — BAQ off, DUST off, allele-balance on. Keeps true indels and
#                     paralog-region SNPs; AB cleans up the SNP depth-inflation FPs.
#   high-confidence — BAQ on, DUST on, allele-balance on. Strict: sheds more SNP
#                     FPs (low-complexity/alignment) but loses most indels.
# A preset names the output dir after itself. Env toggles still override.
PRESET="${PRESET:-}"
case "$PRESET" in
    high-recall)     NO_BAQ=1; NO_DUST=1; NO_AB=0 ;;
    high-confidence) NO_BAQ=0; NO_DUST=0; NO_AB=0 ;;
    "") ;;
    *) echo "unknown PRESET: $PRESET (expected high-recall|high-confidence)" >&2; exit 2 ;;
esac

PILEUP_EXTRA=(); VARCALL_EXTRA=(); VARIANT="ours"
if [[ "$NO_BAQ" == "1" ]]; then PILEUP_EXTRA+=(--no-baq); VARIANT="${VARIANT}_nobaq"; fi
if [[ "$NO_DUST" == "1" ]]; then VARCALL_EXTRA+=(--no-complexity-filter); VARIANT="${VARIANT}_nodust"; fi
if [[ -n "$MIN_MAPQ" ]]; then PILEUP_EXTRA+=(--min-mapq "$MIN_MAPQ"); VARIANT="${VARIANT}_mapq${MIN_MAPQ}"; fi
if [[ "$NO_AB" == "1" ]]; then VARCALL_EXTRA+=(--no-allele-balance-filter); VARIANT="${VARIANT}_noab"; fi
if [[ -n "$AB_LR" ]]; then VARCALL_EXTRA+=(--min-allele-balance-log-lr="$AB_LR"); VARIANT="${VARIANT}_ablr${AB_LR}"; fi
# A named preset owns the output dir name outright.
[[ -n "$PRESET" ]] && VARIANT="$PRESET"
OUT_DIR="$BENCH_DIR/results/per_sample/$COVERAGE/$VARIANT"

# --- sample -> alignment filename / BED filename ---------------------------
# Two naming schemes:
#   300x  — the original selection: HG002 is a CRAM (carries @SQ M5);
#           HG003/HG004 are `*_bench_azar_merged_100.sorted` BAMs that LACK M5
#           (our pileup requires it), so we use the `.m5.bam` reheadered copies
#           produced by add_m5_to_bams.sh.
#   else  — the downsampled tiers (5x..50x): uniform `{sample}.{cov}.seed42.bam`
#           for all three samples; these already carry @SQ M5.
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
            if [[ -f "$m5" ]]; then echo "$m5"
            else
                echo "missing M5-reheadered BAM: $m5" >&2
                echo "run: benchmarks/giab/src/add_m5_to_bams.sh $COVERAGE" >&2
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

# --- binary discovery (container build first, then host) -------------------
discover_bin() {
    if [[ -z "${POP_VAR_CALLER_BIN:-}" ]]; then
        local candidate
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

discover_bin
preflight "$REFERENCE" "${REFERENCE}.fai"
mkdir -p "$OUT_DIR"

echo "binary    : $POP_VAR_CALLER_BIN"
echo "dataset   : per_sample / $COVERAGE"
echo "reference : $REFERENCE"
echo "threads   : $THREADS"
echo "preset    : ${PRESET:-<none>}"
echo "config    : BAQ=$([[ $NO_BAQ == 1 ]] && echo off || echo on) DUST=$([[ $NO_DUST == 1 ]] && echo off || echo on) allele-balance=$([[ $NO_AB == 1 ]] && echo off || echo "on(${AB_LR:--5})") (variant=$VARIANT)"
echo "samples   : ${SAMPLES[*]}"
echo "out dir   : $OUT_DIR"
echo

for sample in "${SAMPLES[@]}"; do
    aln="$(align_file "$sample")"
    bed="$BED_DIR/$(bed_file "$sample")"
    psp="$OUT_DIR/${sample}.psp"
    vcf="$OUT_DIR/${sample}.vcf"
    pileup_log="$OUT_DIR/${sample}.pileup.log"
    varcall_log="$OUT_DIR/${sample}.varcall.log"

    preflight "$aln" "$bed"

    echo "=================================================================="
    echo "sample    : $sample"
    echo "alignment : $aln"
    echo "regions   : $bed ($(wc -l < "$bed") intervals)"
    echo "psp       : $psp"
    echo "vcf       : $vcf"
    echo

    t0=$(date +%s)
    echo "[pileup] $sample -> $psp"
    run "$POP_VAR_CALLER_BIN" pileup \
        --reference "$REFERENCE" \
        --regions "$bed" \
        --output "$psp" \
        --threads "$THREADS" \
        ${PILEUP_EXTRA[@]+"${PILEUP_EXTRA[@]}"} \
        "$aln" 2> >(tee "$pileup_log" >&2)
    t_mid=$(date +%s)

    echo "[var-calling] $psp -> $vcf"
    run "$POP_VAR_CALLER_BIN" var-calling \
        --reference "$REFERENCE" \
        --regions "$bed" \
        --output "$vcf" \
        --threads "$THREADS" \
        ${VARCALL_EXTRA[@]+"${VARCALL_EXTRA[@]}"} \
        "$psp" 2> >(tee "$varcall_log" >&2)
    t1=$(date +%s)

    [[ "${DRY_RUN:-0}" == "1" ]] && { echo; continue; }
    echo
    echo "stage 1 (pileup)      : $((t_mid - t0)) s"
    echo "stage 2 (var-calling) : $((t1 - t_mid)) s"
    echo "total   elapsed       : $((t1 - t0)) s"
    echo "records               : $(record_count "$vcf")"
    echo
done

echo "done. VCFs under: $OUT_DIR"
