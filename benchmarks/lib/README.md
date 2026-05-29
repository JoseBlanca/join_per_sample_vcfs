# Shared benchmark runners

Caller runners + helpers shared across benchmarks. Each benchmark
(`benchmarks/<name>/`) contributes only a `bench.config.sh` describing
its paths and knobs; the scripts here hold all the common scaffolding
(binary discovery, preflight, timing, record counting).

## Layout

```
benchmarks/
  lib/
    common.sh             # sourced by the runners; config + helpers
    run_ours.sh           # pop_var_caller   (single | cohort)
    run_gatk.sh           # GATK             (single | cohort)
    run_freebayes.sh      # freebayes        (single | cohort)
    prepare_reference.sh  # build .fai (+ .dict for GATK) for a reference
    compare_to_truth.sh   # precision/recall/F1 vs a truth VCF (accuracy)
  <name>/
    bench.config.sh       # per-benchmark paths/knobs (see comments inside)
    crams/                # input CRAMs (+ .crai)
    results/              # caller outputs land here (results/<caller>/...)
```

## Usage

```sh
# 1. one-time: make sure the reference has its index siblings
benchmarks/lib/prepare_reference.sh benchmarks/<name>/bench.config.sh

# 2. run a caller — mode is `single` (one sample) or `cohort` (all)
benchmarks/lib/run_ours.sh      benchmarks/<name>/bench.config.sh single
benchmarks/lib/run_gatk.sh      benchmarks/<name>/bench.config.sh cohort
benchmarks/lib/run_freebayes.sh benchmarks/<name>/bench.config.sh single

# 3. (accuracy benchmarks) score caller VCFs against a truth set.
#    With no VCF args it scores the three single-sample outputs.
benchmarks/lib/compare_to_truth.sh benchmarks/<name>/bench.config.sh
```

`DRY_RUN=1` makes the runners print the exact command they would run
(shell-quoted) instead of executing — handy for checking what a config
resolves to without needing every tool installed. Every config value is
also overridable from the environment (`REFERENCE=… THREADS=8 …`).

## Benchmarks

- **tomato1** — 26-sample *S. lycopersicum* cohort, CRAMs pre-sliced to
  ~2 Mb of SL4.0. Multi-sample; the perf experiments (`tomato1/scripts/perf_*.py`)
  build their PSP/GVCF inputs via the cohort runners.
- **human_genome_bottle** — single-sample GIAB HG002, CRAM restricted to
  the 1000-region (~5 Mb) benchmark BED on GRCh38. Accuracy benchmark:
  compare each caller against the GIAB truth VCF with `compare_to_truth.sh`.

## Notes

- GATK lives at `/opt/gatk/gatk` inside the dev container; override with
  `GATK_BIN=…`. `pop_var_caller` is auto-detected from
  `target-container/release` then `target/release` (override with
  `POP_VAR_CALLER_BIN=…`).
- `pop_var_caller` has no BED/region flag, so `run_ours.sh` processes the
  whole CRAM — fine here because the benchmark CRAMs are already
  pre-sliced to the region set. GATK and freebayes restrict via the BED.
- `compare_to_truth.sh` scores **allele** concordance (POS+REF+ALT), not
  genotype concordance; for GT-level eval use rtg vcfeval or hap.py.
