# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Scaling perf experiment for freebayes (CRAM -> VCF) via per-region
parallelisation.

Freebayes' canonical scaling pattern is "split the genome into
regions, run independent freebayes processes per region, merge the
output VCFs" — that's what `freebayes-parallel` does. On macOS
Homebrew that wrapper isn't packaged with freebayes (vcflib + GNU
parallel would also need separate installs), so we drive the same
strategy from Python: one freebayes process per BED region (each
writing directly via `-v`), pool of up to THREADS concurrent
workers, then concatenate the per-region VCFs at the end.

CPU budget matches perf_ours_whole_pipeline (both use THREADS workers),
so the freebayes vs ours_whole_pipeline numbers are an apples-to-apples
"full cohort from BAMs, parallelism budget = THREADS" comparison.

Wall time = freebayes pool wall + merge time. Peak RSS is the max
sum-across-alive-workers observed by the polling loop in
measure_pool (the merge step is trivial Python, negligible RSS).
The post-call MIN_QUAL=30 filter that run_freebayes_*.sh applies is
*not* re-applied here — that's a separate awk pass that doesn't
affect freebayes' own footprint.

Sample sizes: 1, 2, 4, 8, 12, 16, 20, 24, 26 (full tomato1 cohort).
Inputs:       benchmarks/tomato1/crams/*.bench.cram
              benchmarks/tomato1/regions.bed
Output:       benchmarks/tomato1/results/perf/freebayes.tsv
              benchmarks/tomato1/results/perf/freebayes/N<nn>.vcf
              benchmarks/tomato1/results/perf/freebayes/N<nn>/region<NN>.vcf

Env overrides:
  FREEBAYES_BIN   binary (default: auto-detect `freebayes` on PATH)
  REFERENCE       SL4.0 fasta (.fai sibling required)
  PLOIDY          --ploidy (default: 2)
  THREADS         worker concurrency (default: 4 — matches the pileup
                  concurrency in perf_ours_whole_pipeline so CPU budgets align)
  SIZES           comma-separated sample sizes (default: 1,2,4,8,12,16,20,24,26)

Invoke:
  uv run --script benchmarks/tomato1/scripts/perf_freebayes.py
"""

import os
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    CRAM_DIR, DEFAULT_BED, DEFAULT_REFERENCE, DEFAULT_THREADS, PERF_DIR,
    Measurement, banner, check_exists, list_inputs, measure_pool,
    pick_subset, sizes_from_env, write_tsv,
)

CALLER = "freebayes"
BIN = os.environ.get("FREEBAYES_BIN") or shutil.which("freebayes")
PLOIDY = os.environ.get("PLOIDY", "2")
OUT_DIR = PERF_DIR / CALLER


def read_regions(bed: Path) -> list[tuple[str, int, int]]:
    """Parse BED (0-based, half-open) into (chrom, start1, end1)
    triples in 1-based-inclusive form, ready for freebayes' --region
    flag (which expects 1-based inclusive)."""
    regions = []
    with bed.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            regions.append((chrom, start + 1, end))
    if not regions:
        sys.exit(f"no usable BED rows in {bed}")
    return regions


def merge_region_vcfs(region_vcfs: list[Path], out_vcf: Path) -> None:
    """Concatenate per-region VCFs in input order. Header from the
    first file is preserved; subsequent files have their headers
    stripped. Regions are non-overlapping (they came from a BED), so
    no sort or dedup is needed."""
    with out_vcf.open("w") as dst:
        for i, src in enumerate(region_vcfs):
            with src.open() as fh:
                for line in fh:
                    if i > 0 and line.startswith("#"):
                        continue
                    dst.write(line)


def main() -> int:
    if not BIN:
        sys.exit("no freebayes binary on PATH; install via brew/apt or set FREEBAYES_BIN")
    sizes = sizes_from_env()
    check_exists(
        Path(BIN), DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai"), DEFAULT_BED,
    )
    crams = list_inputs(CRAM_DIR, ".bench.cram")
    regions = read_regions(DEFAULT_BED)
    banner(CALLER, sizes)
    print(f"regions:   {len(regions)} (from {DEFAULT_BED.name})")
    print(f"workers:   {DEFAULT_THREADS}")
    print()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[Measurement] = []
    for n in sizes:
        subset = pick_subset(crams, n)
        n_dir = OUT_DIR / f"N{n:02d}"
        n_dir.mkdir(parents=True, exist_ok=True)
        out_vcf = OUT_DIR / f"N{n:02d}.vcf"

        region_vcfs = []
        cmds = []
        for i, (chrom, start, end) in enumerate(regions):
            region_vcf = n_dir / f"region{i:03d}.vcf"
            region_vcfs.append(region_vcf)
            cmds.append([
                BIN,
                "-f", str(DEFAULT_REFERENCE),
                "-r", f"{chrom}:{start}-{end}",
                "-v", str(region_vcf),
                "-p", PLOIDY,
                *[str(c) for c in subset],
            ])

        concurrency = min(DEFAULT_THREADS, len(cmds))
        print(
            f"[{CALLER}] N={n:>2}: {len(cmds)} per-region freebayes @ "
            f"concurrency={concurrency}"
        )
        pool_wall, peak_bytes, exit_code = measure_pool(cmds, max_concurrent=concurrency)
        peak_mb = peak_bytes / 1024 / 1024
        print(f"  pool   -> {pool_wall:.1f}s, peak {peak_mb:.0f} MB")

        if exit_code != 0:
            rows.append(Measurement(CALLER, n, pool_wall, peak_mb, exit_code))
            print(f"  per-region freebayes FAILED (exit {exit_code}); skipping merge")
            continue

        t_merge0 = time.monotonic()
        merge_region_vcfs(region_vcfs, out_vcf)
        merge_wall = time.monotonic() - t_merge0
        total_wall = pool_wall + merge_wall
        rows.append(Measurement(CALLER, n, total_wall, peak_mb, exit_code))
        print(f"  merge  -> {merge_wall:.2f}s")
        print(f"  TOTAL  -> {total_wall:.1f}s, peak {peak_mb:.0f} MB  [ok]")

    tsv = PERF_DIR / f"{CALLER}.tsv"
    write_tsv(rows, tsv)
    print(f"\nwrote {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
