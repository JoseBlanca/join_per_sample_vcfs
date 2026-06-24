#!/usr/bin/env bash
# Benchmark config: ssr_tomato1 — SSR (microsatellite) genotyping on the same
# ~63-sample S. lycopersicum cohort as the SNP `tomato1` bench. Parallel to the
# SNP benches, but the pipeline is three stages:
#   ssr-catalog (trf-mod, once) -> ssr-pileup (per sample) -> ssr-call (cohort).
# Run it with benchmarks/lib/run_ours_ssr.sh. Every value is env-overridable.
#
# NOTE: the SSR caller has no region flag (like the SNP caller); the cohort
# CRAMs are already pre-sliced to regions.bed (~2 Mb of SL4.0), so there is no
# BED here. ssr-catalog scans the WHOLE reference, so most catalog loci sit
# outside the sliced regions and simply no-call (LowDepth) — harmless.

BENCH_NAME="ssr_tomato1"

# Reference: the SAME SL4.0 the cohort CRAMs were aligned to. The catalog's
# whole-reference MD5 is validated against every .ssr.psp, so the catalog and
# the pileups must use this exact FASTA. A sibling .fai is required.
REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"

# Reuse the SNP tomato1 cohort CRAMs (gitignored local data). Defaults to the
# sibling SNP bench's crams; override CRAM_DIR if they live elsewhere (e.g. in a
# separate working tree).
CRAM_DIR="${CRAM_DIR:-$BENCH_DIR/../tomato1/crams}"
CRAM_GLOB="${CRAM_GLOB:-*.bench.cram}"
CRAM_SUFFIX="${CRAM_SUFFIX:-.bench.cram}"
SINGLE_CRAM="${SINGLE_CRAM:-SRR7279481.p1.bench.cram}"

PLOIDY="${PLOIDY:-2}"
THREADS="${THREADS:-4}"

# --- SSR-specific knobs ---
# Path to trf-mod (ssr-catalog finds it on PATH by default; set to override).
TRF_MOD_PATH="${TRF_MOD_PATH:-}"
# Appended to `ssr-catalog` (e.g. --min-purity / --min-score / --flank-bp).
CATALOG_EXTRA="${CATALOG_EXTRA:-}"
# Appended to each `ssr-pileup`. `--build-index-if-missing` rebuilds a missing
# .crai in place rather than erroring on an un-indexed input.
PILEUP_EXTRA="${PILEUP_EXTRA:---build-index-if-missing}"
# Appended to `ssr-call`.
CALL_EXTRA="${CALL_EXTRA:-}"
# CATALOG defaults to $OUT_ROOT/ours/$BENCH_NAME.ssr.catalog inside the runner
# (OUT_ROOT is resolved after this file is sourced); set CATALOG to override.
