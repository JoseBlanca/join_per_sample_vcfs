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
MIN_QUAL="${MIN_QUAL:-30}"
GATK_HEAP="${GATK_HEAP:-4g}"
GATK_COMBINE_HEAP="${GATK_COMBINE_HEAP:-8g}"
