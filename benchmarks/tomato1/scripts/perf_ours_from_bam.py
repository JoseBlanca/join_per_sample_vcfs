# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Scaling perf experiment for the DIRECT path (`var-calling-from-bam`:
CRAM -> single-sample VCF, no `.psp` intermediate) — the ours counterpart
to freebayes' one-shot CRAM -> VCF.

freebayes has no internal threading; it parallelises by running region
processes concurrently. The direct path is single-sample, so the analogous
parallelisation is one `var-calling-from-bam` process per sample, run up to
THREADS concurrently (each single-threaded). Wall = total time until every
sample's VCF exists; peak RSS = max over polling ticks of the summed RSS
across the concurrently-alive processes.

CAVEATS — this comparison is deliberately NOT apples-to-apples (documented
so the numbers aren't over-read):
  - No region/BED restriction: the direct path calls GENOME-WIDE, whereas
    freebayes is restricted to regions.bed — so ours does substantially
    more work per sample (~36 s single-threaded vs the ~2 Mb of regions).
  - Single-sample: N processes -> N single-sample VCFs, vs freebayes' one
    joint multi-sample VCF.
  - It uses the OLD streaming per-record architecture, not the new columnar
    cohort path (so it does not reflect the columnar/DUST-pool work).

The (perf_ours_from_bam, perf_freebayes) pair is the "direct CRAM -> VCF"
comparison; the cohort columnar path is measured by perf_ours_whole_pipeline
(parallel pileup + one-process joint) and perf_ours_joint.

Sample sizes: 1, 2, 4, 8, 12, 16, 20, 24, 26 (full tomato1 cohort).
Inputs:       benchmarks/tomato1/crams/*.bench.cram
Output:       benchmarks/tomato1/results/perf/ours_from_bam.tsv
              benchmarks/tomato1/results/perf/ours_from_bam/N<nn>/<sample>.vcf

Env overrides:
  POP_VAR_CALLER_BIN  binary (default: $PROJECT_ROOT/target/release/pop_var_caller)
  REFERENCE           SL4.0 fasta (.fai sibling required)
  THREADS             process concurrency (default: 4)
  SIZES               comma-separated sample sizes (default: 1,2,4,8,12,16,20,24,26)

Invoke:
  uv run --script benchmarks/tomato1/scripts/perf_ours_from_bam.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    CRAM_DIR, DEFAULT_REFERENCE, DEFAULT_THREADS, PERF_DIR, PROJECT_ROOT,
    Measurement, banner, check_exists, list_inputs, measure_pool,
    pick_subset, sizes_from_env, write_tsv,
)

CALLER = "ours_from_bam"
BIN = Path(os.environ.get(
    "POP_VAR_CALLER_BIN",
    str(PROJECT_ROOT / "target" / "release" / "pop_var_caller"),
))
OUT_DIR = PERF_DIR / CALLER


def main() -> int:
    sizes = sizes_from_env()
    check_exists(BIN, DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai"))
    crams = list_inputs(CRAM_DIR, ".bench.cram")
    banner(CALLER, sizes)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[Measurement] = []
    for n in sizes:
        subset = pick_subset(crams, n)
        n_dir = OUT_DIR / f"N{n:02d}"
        n_dir.mkdir(parents=True, exist_ok=True)

        # One var-calling-from-bam process per sample (single-threaded);
        # run up to THREADS concurrently — the per-sample-process analogue
        # of freebayes' region-process parallelism. Wall is end-to-end (all
        # VCFs produced); peak RSS sums the concurrently-alive processes.
        cmds = []
        for cram in subset:
            base = cram.name.removesuffix(".bench.cram")
            vcf = n_dir / f"{base}.vcf"
            cmds.append([
                str(BIN), "var-calling-from-bam",
                "--reference", str(DEFAULT_REFERENCE),
                "--output", str(vcf),
                "--threads", str(DEFAULT_THREADS),
                str(cram),
            ])
        # Our caller parallelises with threads inside one process, so each
        # per-sample var-calling-from-bam runs with --threads N, ONE at a
        # time (never several in parallel). Wall = sum over samples; peak
        # RSS = a single process. (freebayes can only parallelise by
        # running region processes concurrently — it has no threads.)
        print(f"[{CALLER}] N={n:>2}: {len(cmds)} var-calling-from-bam, sequential @ --threads {DEFAULT_THREADS}")
        wall, peak_bytes, exit_code = measure_pool(cmds, max_concurrent=1)
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
