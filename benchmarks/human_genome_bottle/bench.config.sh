#!/usr/bin/env bash
# Benchmark config: human_genome_bottle — single-sample GIAB HG002,
# CRAM pre-restricted to the 1000-region benchmark BED (~5 Mb of GRCh38).
# This is an accuracy benchmark: each caller's calls are compared to the
# GIAB truth VCF with benchmarks/lib/compare_to_truth.sh.
#
# Sourced by benchmarks/lib/run_<caller>.sh and the helper scripts.
# Every value is env-overridable.

BENCH_NAME="human_genome_bottle"

# GRCh38 analysis set (no-alt + hs38d1). A sibling .fai and .dict are
# required; build them with benchmarks/lib/prepare_reference.sh.
REFERENCE="${REFERENCE:-$HOME/genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna}"

# High-confidence benchmark regions (1000 intervals).
BED="${BED:-HG002_bench_azar_merged_1000.bed}"

# Only one CRAM here; single and cohort modes produce the same VCF.
CRAM_GLOB="${CRAM_GLOB:-*.cram}"
CRAM_SUFFIX="${CRAM_SUFFIX:-.cram}"
SINGLE_CRAM="${SINGLE_CRAM:-HG002_reads_selected_1000_rg.cram}"

# Reference callset (curated baseline for this benchmark) used by
# compare_to_truth.sh as the truth set.
TRUTH_VCF="${TRUTH_VCF:-$BENCH_DIR/results/HG002_bottle_reference_1000rg.vcf.gz}"

PLOIDY="${PLOIDY:-2}"
THREADS="${THREADS:-4}"

# Fair cross-caller QUAL comparison via a ROC-style sweep. Each caller picks a
# different default QUAL gate (ours --min-qual, GATK -stand-call-conf, freebayes
# a post-call cap), so comparing at those defaults is apples-to-oranges. Instead
# every caller emits its FULL QUAL range (gate = 0) and the dashboard sweeps a
# common QUAL cutoff to trace each caller's precision/recall curve and pick its
# F1-optimal threshold. So the fairness lives in the sweep, not in matched
# pre-filters — keep all three emission gates at 0 here.
MIN_QUAL="${MIN_QUAL:-0}"                       # freebayes: no post-call QUAL cap
HC_EXTRA="${HC_EXTRA:--stand-call-conf 0}"      # GATK HaplotypeCaller: emit low-conf too

# ours: emit every site (--min-qual 0) for the full-range QUAL sweep. BAQ and
# the DUST low-complexity filter are kept ON by default (the caller's real
# defaults) for now — note this penalises our indel/SNP recall vs GATK/freebayes
# in low-complexity regions (DUST masks the microsatellites/homopolymers where
# indels live: ~222 of 311 missed truth indels, plus some true SNPs). Revisit
# as a caller-level class-aware-DUST question rather than a benchmark override.
VARCALL_EXTRA="${VARCALL_EXTRA:---min-qual 0}"
GATK_HEAP="${GATK_HEAP:-4g}"
GATK_COMBINE_HEAP="${GATK_COMBINE_HEAP:-8g}"
