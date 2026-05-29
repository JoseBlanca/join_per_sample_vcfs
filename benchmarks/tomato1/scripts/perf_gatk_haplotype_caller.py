# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Scaling perf experiment for GATK HaplotypeCaller (CRAM -> GVCF,
per-sample intermediate).

Per-sample step of GATK's joint genotyping flow, measured in
isolation. Runs THREADS HaplotypeCaller processes concurrently, each
with --native-pair-hmm-threads=1 so the CPU budget is split across
processes (matches perf_ours_pileup's strategy — same concurrency,
same threads-per-process).

The (perf_gatk_haplotype_caller, perf_ours_pileup) pair shows how
the two per-sample intermediates scale relative to each other when
built under the same parallelism budget.

Per-process Java heap is capped via JAVA_HEAP (default 2g) so that
THREADS=4 concurrent HC processes fit inside the container's memory
budget (4 x 2g = 8g, comfortably under the 16g `scripts/dev.sh`
allocates to Apple container).

Note on availability: GATK isn't in Homebrew. On macOS run this
script in the container (after the Containerfile additions are
built); on Linux either path works.

Sample sizes: 1, 2, 4, 8, 12, 16, 20, 24, 26 (full tomato1 cohort).
Inputs:       benchmarks/tomato1/crams/*.bench.cram
              benchmarks/tomato1/regions.bed
Output:       benchmarks/tomato1/results/perf/gatk_haplotype_caller.tsv
              benchmarks/tomato1/results/perf/gatk_haplotype_caller/N<nn>/<sample>.g.vcf.gz

Env overrides:
  GATK_BIN    binary (default: /opt/gatk/gatk)
  REFERENCE   SL4.0 fasta (.fai + .dict siblings required)
  JAVA_HEAP   -Xmx per HC process (default: 2g)
  THREADS     HC concurrency (default: 4)
  SIZES       comma-separated sample sizes (default: 1,2,4,8,12,16,20,24,26)

Invoke (container path on macOS):
  DEV_EXTRA_MOUNT=$HOME/genomes ./scripts/dev.sh \\
      python3 benchmarks/tomato1/scripts/perf_gatk_haplotype_caller.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    CRAM_DIR, DEFAULT_BED, DEFAULT_REFERENCE, DEFAULT_THREADS, PERF_DIR,
    Measurement, banner, check_exists, list_inputs, measure_pool,
    pick_subset, sizes_from_env, write_tsv,
)

CALLER = "gatk_haplotype_caller"
BIN = Path(os.environ.get("GATK_BIN", "/opt/gatk/gatk"))
JAVA_HEAP = os.environ.get("JAVA_HEAP", "2g")
OUT_DIR = PERF_DIR / CALLER


def main() -> int:
    sizes = sizes_from_env()
    ref_dict = DEFAULT_REFERENCE.with_suffix(".dict")
    check_exists(
        BIN, DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai"),
        ref_dict, DEFAULT_BED,
    )
    crams = list_inputs(CRAM_DIR, ".bench.cram")
    banner(CALLER, sizes)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[Measurement] = []
    for n in sizes:
        subset = pick_subset(crams, n)
        n_dir = OUT_DIR / f"N{n:02d}"
        n_dir.mkdir(parents=True, exist_ok=True)

        cmds = []
        for cram in subset:
            base = cram.name.removesuffix(".bench.cram")
            gvcf = n_dir / f"{base}.g.vcf.gz"
            cmds.append([
                str(BIN), "--java-options", f"-Xmx{JAVA_HEAP}", "HaplotypeCaller",
                "--reference", str(DEFAULT_REFERENCE),
                "--input", str(cram),
                "--intervals", str(DEFAULT_BED),
                "--output", str(gvcf),
                "--emit-ref-confidence", "GVCF",
                "--native-pair-hmm-threads", "1",
            ])
        concurrency = min(DEFAULT_THREADS, len(cmds))
        print(
            f"[{CALLER}] N={n:>2}: {len(cmds)} HaplotypeCallers @ "
            f"concurrency={concurrency} (heap={JAVA_HEAP} each)"
        )
        wall, peak_bytes, exit_code = measure_pool(cmds, max_concurrent=concurrency)
        peak_mb = peak_bytes / 1024 / 1024
        rows.append(Measurement(CALLER, n, wall, peak_mb, exit_code))
        status = "ok" if exit_code == 0 else f"FAILED (exit {exit_code})"
        print(f"  -> {wall:.1f}s, peak {peak_mb:.0f} MB  [{status}]")

    tsv = PERF_DIR / f"{CALLER}.tsv"
    write_tsv(rows, tsv)
    print(f"\nwrote {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
