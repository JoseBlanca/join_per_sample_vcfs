#!/usr/bin/env bash
# Cohort GATK test: all 18 CRAMs -> per-sample GVCFs -> combined GVCF
# -> joint VCF, following the GATK best-practices three-stage flow.
#
# Stage 1 (per sample): HaplotypeCaller -ERC GVCF over each *.bench.cram.
# Stage 2 (cohort):     CombineGVCFs over the 18 GVCFs.
# Stage 3 (cohort):     GenotypeGVCFs on the combined GVCF.
#
# Stage 1 runs sequentially (each at THREADS native-pair-hmm threads)
# to stay within the host P-core budget (memory:
# feedback_host_p_core_budget). Already-built GVCFs are skipped on
# re-runs so iteration is cheap.
#
# Output layout:
#   tmp/tomato_cohort_test/results/gatk/cohort/gvcf/<sample>.g.vcf.gz
#   tmp/tomato_cohort_test/results/gatk/cohort/combined.g.vcf.gz
#   tmp/tomato_cohort_test/results/gatk/cohort/cohort.vcf
#   tmp/tomato_cohort_test/results/gatk/cohort/*.log
#
# Env overrides:
#   GATK_BIN        binary to invoke (default: /opt/gatk/gatk)
#   REFERENCE       SL4.0 fasta (.fai + .dict siblings required)
#   THREADS         --native-pair-hmm-threads per HaplotypeCaller call
#                   (default 4)
#   JAVA_HEAP       -Xmx (default 4g for stages 1/3, 8g for stage 2)
#   HC_EXTRA        appended to each HaplotypeCaller invocation
#   COMBINE_EXTRA   appended to the CombineGVCFs invocation
#   GENOTYPE_EXTRA  appended to the GenotypeGVCFs invocation
#
# Notes vs other callers:
# - GATK is the only one of the three with a per-sample intermediate
#   (GVCF) that maps to our caller's PSP. Stage 2 (CombineGVCFs) is the
#   GATK analogue of `pop_var_caller var-calling`'s input merging.
# - For cohorts in the thousands, GenomicsDBImport scales better than
#   CombineGVCFs. At 18 samples CombineGVCFs is fine and simpler.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

CRAM_DIR="$TEST_DIR/crams"
BED="$TEST_DIR/regions.bed"
REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"
THREADS="${THREADS:-4}"
JAVA_HEAP="${JAVA_HEAP:-4g}"
HC_EXTRA=${HC_EXTRA:-}
COMBINE_EXTRA=${COMBINE_EXTRA:-}
GENOTYPE_EXTRA=${GENOTYPE_EXTRA:-}

GATK_BIN="${GATK_BIN:-/opt/gatk/gatk}"
[[ -x "$GATK_BIN" ]] || { echo "GATK binary not executable: $GATK_BIN" >&2; exit 1; }

# --- preflight ---
REF_DICT="${REFERENCE%.fa}.dict"
for f in "$BED" "$REFERENCE" "${REFERENCE}.fai" "$REF_DICT"; do
    [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done
shopt -s nullglob
crams=("$CRAM_DIR"/*.bench.cram)
if (( ${#crams[@]} == 0 )); then
    echo "no *.bench.cram in $CRAM_DIR" >&2
    exit 1
fi
for c in "${crams[@]}"; do
    [[ -f "${c}.crai" ]] || { echo "missing index: ${c}.crai" >&2; exit 1; }
done

# --- output ---
OUT_DIR="$TEST_DIR/results/gatk/cohort"
GVCF_DIR="$OUT_DIR/gvcf"
COMBINED_GVCF="$OUT_DIR/combined.g.vcf.gz"
COHORT_VCF="$OUT_DIR/cohort.vcf"
COMBINE_LOG="$OUT_DIR/combine.log"
GENOTYPE_LOG="$OUT_DIR/genotype.log"
mkdir -p "$GVCF_DIR"

echo "binary    : $GATK_BIN ($($GATK_BIN --version 2>&1 | head -1))"
echo "samples   : ${#crams[@]}"
echo "reference : $REFERENCE"
echo "regions   : $BED ($(wc -l < "$BED") intervals)"
echo "threads   : $THREADS (per HaplotypeCaller)"
echo "heap      : $JAVA_HEAP"
echo "gvcf dir  : $GVCF_DIR"
echo "vcf out   : $COHORT_VCF"
echo

# ---------- Stage 1: per-sample HaplotypeCaller GVCF (sequential) ----------
t_stage1=$(date +%s)
for cram in "${crams[@]}"; do
    base=$(basename "$cram" .bench.cram)   # SRR....p1
    gvcf="$GVCF_DIR/${base}.g.vcf.gz"
    log="$GVCF_DIR/${base}.log"
    if [[ -s "$gvcf" && -s "${gvcf}.tbi" ]]; then
        echo "[skip HC] $base (GVCF already present)"
        continue
    fi
    echo "[HC] $base"
    "$GATK_BIN" --java-options "-Xmx${JAVA_HEAP}" HaplotypeCaller \
        --reference "$REFERENCE" \
        --input "$cram" \
        --intervals "$BED" \
        --output "$gvcf" \
        --emit-ref-confidence GVCF \
        --native-pair-hmm-threads "$THREADS" \
        $HC_EXTRA \
        2> >(tee "$log" >&2)
done
t_stage1_end=$(date +%s)
echo
echo "stage 1 elapsed: $((t_stage1_end - t_stage1)) s"
echo

# ---------- Stage 2: CombineGVCFs ----------
gvcfs=("$GVCF_DIR"/*.g.vcf.gz)
if (( ${#gvcfs[@]} != ${#crams[@]} )); then
    echo "WARNING: expected ${#crams[@]} GVCFs, found ${#gvcfs[@]}" >&2
fi

# Build -V flags
combine_args=()
for g in "${gvcfs[@]}"; do
    combine_args+=(--variant "$g")
done

COMBINE_HEAP="${COMBINE_HEAP:-8g}"
echo "[CombineGVCFs] ${#gvcfs[@]} GVCFs -> $COMBINED_GVCF (heap=$COMBINE_HEAP)"
t_stage2=$(date +%s)
"$GATK_BIN" --java-options "-Xmx${COMBINE_HEAP}" CombineGVCFs \
    --reference "$REFERENCE" \
    --intervals "$BED" \
    "${combine_args[@]}" \
    --output "$COMBINED_GVCF" \
    $COMBINE_EXTRA \
    2> >(tee "$COMBINE_LOG" >&2)
t_stage2_end=$(date +%s)
echo "stage 2 elapsed: $((t_stage2_end - t_stage2)) s"
echo

# ---------- Stage 3: GenotypeGVCFs ----------
echo "[GenotypeGVCFs] -> $COHORT_VCF"
t_stage3=$(date +%s)
"$GATK_BIN" --java-options "-Xmx${JAVA_HEAP}" GenotypeGVCFs \
    --reference "$REFERENCE" \
    --intervals "$BED" \
    --variant "$COMBINED_GVCF" \
    --output "$COHORT_VCF" \
    $GENOTYPE_EXTRA \
    2> >(tee "$GENOTYPE_LOG" >&2)
t_stage3_end=$(date +%s)

echo
echo "stage 3 elapsed: $((t_stage3_end - t_stage3)) s"
echo "total   elapsed: $((t_stage3_end - t_stage1)) s"
echo "records: $(grep -vc '^#' "$COHORT_VCF")"
echo "samples in vcf: $(grep -m1 '^#CHROM' "$COHORT_VCF" | awk '{print NF-9}')"
