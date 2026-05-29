# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Scaling perf experiment for `pop_var_caller pileup` (CRAM -> PSP,
per-sample intermediate).

Per-sample step of the cohort pipeline, measured in isolation. Runs
THREADS pileup processes concurrently (each single-threaded —
parallelism budget exhausted across processes rather than within
one), matching how stage 1 of perf_ours_whole_pipeline runs and how
a user would parallelise pileups across CRAMs in practice.

The (perf_ours_pileup, perf_gatk_haplotype_caller) pair shows how
the two per-sample intermediates scale relative to each other when
produced under the same parallelism budget — analogous in scope
even though the two tools' per-sample artefacts differ in content
(PSP vs GVCF).

Sample sizes: 1, 2, 4, 8, 12, 16, 20, 24, 26 (full tomato1 cohort).
Inputs:       benchmarks/tomato1/crams/*.bench.cram
Output:       benchmarks/tomato1/results/perf/ours_pileup.tsv
              benchmarks/tomato1/results/perf/ours_pileup/N<nn>/<sample>.psp

Env overrides:
  POP_VAR_CALLER_BIN  binary (default: $PROJECT_ROOT/target/release/pop_var_caller)
  REFERENCE           SL4.0 fasta (.fai sibling required)
  THREADS             pileup concurrency (default: 4)
  SIZES               comma-separated sample sizes (default: 1,2,4,8,12,16,20,24,26)

Invoke:
  uv run --script benchmarks/tomato1/scripts/perf_ours_pileup.py
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

CALLER = "ours_pileup"
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

        cmds = []
        for cram in subset:
            base = cram.name.removesuffix(".bench.cram")
            psp = n_dir / f"{base}.psp"
            cmds.append([
                str(BIN), "pileup",
                "--reference", str(DEFAULT_REFERENCE),
                "--output", str(psp),
                "--threads", "1",
                str(cram),
            ])
        concurrency = min(DEFAULT_THREADS, len(cmds))
        print(f"[{CALLER}] N={n:>2}: {len(cmds)} pileups @ concurrency={concurrency}")
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
