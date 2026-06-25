#!/usr/bin/env bash
# HipSTR runner — microsatellite genotyping, the comparison caller for the
# ssr_tomato1 bench (parallel to run_freebayes.sh / run_gatk.sh for SNPs).
#
#   benchmarks/lib/run_hipstr.sh <bench.config.sh> [single|cohort]
#
# HipSTR genotypes a fixed set of STR regions (its --regions BED). To compare
# apples-to-apples with our caller, that BED is *derived from our own
# ssr-catalog* (catalog_to_hipstr_bed.py), so both tools genotype the same
# loci. The catalog must already exist — build it (and our .ssr.psp / VCF)
# first with run_ours_ssr.sh; this runner errors if it's missing.
#
# single : one CRAM            -> single-sample bgzipped str-vcf.
# cohort : all CRAMs (one cmd) -> joint bgzipped str-vcf. HipSTR has no
#          per-sample intermediate and is single-threaded (no --threads),
#          so this is one process — the fair "native parallelism" for HipSTR.
#          Sample ids come from each CRAM's @RG SM tag.
#
# HipSTR reads CRAM natively; --fasta MUST be the FASTA the CRAMs were
# aligned to (same constraint our caller's catalog MD5 enforces).
#
# Env overrides (see common.sh + bench.config.sh):
#   HIPSTR_BIN     binary (default: $PROJECT_ROOT/HipSTR/HipSTR, then PATH)
#   REFERENCE      FASTA (.fai sibling required)
#   CATALOG        our ssr-catalog (default $OUT_ROOT/ours/$BENCH_NAME.ssr.catalog)
#   HIPSTR_REGIONS regions BED (default $OUT_ROOT/hipstr/$BENCH_NAME.hipstr_regions.bed;
#                  rebuilt from CATALOG + the bench BED if missing)
#   HIPSTR_EXTRA   appended verbatim to the HipSTR command line
#   DRY_RUN=1      print commands instead of running them

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

CONFIG="${1:-}"
MODE="${2:-single}"
bench_load_config "$CONFIG"

OUT_DIR="$OUT_ROOT/hipstr"
CATALOG="${CATALOG:-$OUT_ROOT/ours/$BENCH_NAME.ssr.catalog}"
HIPSTR_REGIONS="${HIPSTR_REGIONS:-$OUT_DIR/$BENCH_NAME.hipstr_regions.bed}"
HIPSTR_EXTRA=${HIPSTR_EXTRA:-}
# The catalog is whole-genome but the cohort CRAMs are sliced to a region
# set; restrict HipSTR's loci to that slice (loci outside it only no-call).
# The SSR config carries no BED of its own, so default to the sibling SNP
# bench's regions.bed (the slice the shared CRAMs were cut to).
SLICE_BED="${BED:-$BENCH_DIR/../tomato1/regions.bed}"

# --- HipSTR binary discovery (no common.sh helper; HipSTR-specific) ---
bench_discover_hipstr_bin() {
    if [[ -z "${HIPSTR_BIN:-}" ]]; then
        if [[ -x "$PROJECT_ROOT/HipSTR/HipSTR" ]]; then
            HIPSTR_BIN="$PROJECT_ROOT/HipSTR/HipSTR"
        elif command -v HipSTR >/dev/null; then
            HIPSTR_BIN="$(command -v HipSTR)"
        fi
    fi
    if [[ -z "${HIPSTR_BIN:-}" || ! -x "${HIPSTR_BIN}" ]]; then
        echo "no HipSTR binary found." >&2
        echo "build it: cd HipSTR && make HipSTR CPPFLAGS=-I\$(brew --prefix)/include LDFLAGS=-L\$(brew --prefix)/lib" >&2
        echo "or set HIPSTR_BIN=<path>" >&2
        exit 1
    fi
}

# Build the HipSTR --regions BED from our ssr-catalog (so both callers see
# the identical locus set), restricted to the bench BED. Cached; delete to
# force a rebuild.
ensure_regions() {
    if [[ -s "$HIPSTR_REGIONS" ]]; then
        echo "[skip regions] $HIPSTR_REGIONS ($(wc -l < "$HIPSTR_REGIONS") loci)"
        return 0
    fi
    [[ -s "$CATALOG" ]] || {
        echo "ssr-catalog not found: $CATALOG" >&2
        echo "build it first: benchmarks/lib/run_ours_ssr.sh $CONFIG $MODE" >&2
        exit 1
    }
    mkdir -p "$OUT_DIR"
    local conv="$PROJECT_ROOT/benchmarks/ssr_tomato1/scripts/catalog_to_hipstr_bed.py"
    echo "[regions] $CATALOG (+ $SLICE_BED) -> $HIPSTR_REGIONS"
    [[ "${DRY_RUN:-0}" == "1" ]] && { echo "DRY-RUN: uv run $conv --catalog $CATALOG --regions $SLICE_BED --out $HIPSTR_REGIONS"; return 0; }
    uv run --quiet "$conv" --catalog "$CATALOG" --regions "$SLICE_BED" --out "$HIPSTR_REGIONS"
}

# run_hipstr <out_vcf_gz> <log> <cram>...
# Sample/library names are forced to the CRAM filename base (e.g.
# SRR7279481.p1) via --bam-samps/--bam-libs so they match our caller's
# sample names. By default HipSTR labels by the @RG SM tag (a biosample
# accession, SRS...), which differs from ours and even collides when two
# run files share a biosample — overriding keeps every run distinct and
# the two VCFs joinable by sample name.
run_hipstr() {
    local out_vcf="$1" log="$2"; shift 2
    local bams names=() c
    bams="$(IFS=,; echo "$*")"
    for c in "$@"; do names+=("$(bench_sample_base "$c")"); done
    local samps; samps="$(IFS=,; echo "${names[*]}")"
    # shellcheck disable=SC2086
    bench_run "$log" -- "$HIPSTR_BIN" \
        --bams "$bams" \
        --bam-samps "$samps" \
        --bam-libs "$samps" \
        --fasta "$REFERENCE" \
        --regions "$HIPSTR_REGIONS" \
        --str-vcf "$out_vcf" \
        $HIPSTR_EXTRA
}

run_single() {
    local cram base out_vcf log
    cram="$(bench_single_cram)"
    base="$(bench_sample_base "$cram")"
    bench_discover_hipstr_bin
    bench_preflight "$cram" "${cram}.crai" "$REFERENCE" "${REFERENCE}.fai"

    out_vcf="$OUT_DIR/single_${base}.str.vcf.gz"
    log="$OUT_DIR/single_${base}.log"
    mkdir -p "$OUT_DIR"

    echo "binary    : $HIPSTR_BIN ($($HIPSTR_BIN --version 2>&1 | head -1))"
    echo "mode      : single ($BENCH_NAME)"
    echo "sample    : $base"
    echo "input     : $cram"
    echo "reference : $REFERENCE"
    echo "output    : $out_vcf"
    echo

    ensure_regions
    echo "regions   : $HIPSTR_REGIONS ($(wc -l < "$HIPSTR_REGIONS" 2>/dev/null || echo '?') loci)"
    echo

    local t0 t1
    t0=$(bench_now)
    run_hipstr "$out_vcf" "$log" "$cram"
    t1=$(bench_now)
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0
    echo
    echo "elapsed: $((t1 - t0)) s"
    echo "records: $(bench_record_count "$out_vcf")"
}

run_cohort() {
    bench_list_crams
    bench_discover_hipstr_bin
    bench_preflight "$REFERENCE" "${REFERENCE}.fai"
    local c
    for c in "${BENCH_CRAMS[@]}"; do bench_preflight "${c}.crai"; done

    local out_vcf="$OUT_DIR/cohort.str.vcf.gz"
    local log="$OUT_DIR/cohort.log"
    mkdir -p "$OUT_DIR"

    echo "binary    : $HIPSTR_BIN ($($HIPSTR_BIN --version 2>&1 | head -1))"
    echo "mode      : cohort ($BENCH_NAME) — single-threaded (HipSTR has no --threads)"
    echo "samples   : ${#BENCH_CRAMS[@]}"
    echo "reference : $REFERENCE"
    echo "output    : $out_vcf"
    echo

    ensure_regions
    echo "regions   : $HIPSTR_REGIONS ($(wc -l < "$HIPSTR_REGIONS" 2>/dev/null || echo '?') loci)"
    echo

    local t0 t1
    t0=$(bench_now)
    run_hipstr "$out_vcf" "$log" "${BENCH_CRAMS[@]}"
    t1=$(bench_now)
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0
    echo
    echo "elapsed: $((t1 - t0)) s"
    echo "records: $(bench_record_count "$out_vcf") sites, $(bench_sample_count "$out_vcf") samples"
}

case "$MODE" in
    single) run_single ;;
    cohort) run_cohort ;;
    *) echo "unknown mode: $MODE (expected single|cohort)" >&2; exit 2 ;;
esac
