#!/usr/bin/env bash
# Benchmark config: tomato1 — 26-sample S. lycopersicum cohort, CRAMs
# pre-sliced to the regions in regions.bed (~2 Mb of SL4.0).
#
# Sourced by benchmarks/lib/run_<caller>.sh and the helper scripts.
# Every value is env-overridable. BENCH_DIR, CRAM_DIR, and OUT_ROOT are
# derived by common.sh from this file's location unless set here.

BENCH_NAME="tomato1"

REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"
BED="${BED:-regions.bed}"

# Cohort CRAMs: the bench slices share a .bench.cram suffix; the sample
# base (used for .psp / GVCF / output names) strips that suffix, e.g.
# SRR17274057.p1.bench.cram -> SRR17274057.p1
CRAM_GLOB="${CRAM_GLOB:-*.bench.cram}"
CRAM_SUFFIX="${CRAM_SUFFIX:-.bench.cram}"

# Sample used by `single` mode (matches the original run_*_single.sh).
SINGLE_CRAM="${SINGLE_CRAM:-SRR17274057.p1.bench.cram}"

PLOIDY="${PLOIDY:-2}"
THREADS="${THREADS:-4}"
MIN_QUAL="${MIN_QUAL:-30}"
GATK_HEAP="${GATK_HEAP:-4g}"
GATK_COMBINE_HEAP="${GATK_COMBINE_HEAP:-8g}"
