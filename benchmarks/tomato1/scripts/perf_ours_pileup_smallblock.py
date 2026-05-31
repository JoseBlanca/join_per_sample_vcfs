# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Scaling perf experiment for `pop_var_caller pileup` with
`--block-target-bytes 262144` (256 KiB, 1/4 of the 1 MiB default).

Companion to perf_ours_pileup.py — same workload, same parallelism
budget, only the PSP writer's block target differs. Block size has
little effect on pileup-stage memory (it's a writer flush threshold,
not a working-set knob), so the difference between this and
perf_ours_pileup should be small. The interesting effect shows up
in the joint stage: see perf_ours_joint vs perf_ours_joint_smallblock.

Sample sizes: 1, 2, 4, 8, 12, 16, 20, 24, 26 (full tomato1 cohort).
Inputs:       benchmarks/tomato1/crams/*.bench.cram
Output:       benchmarks/tomato1/results/perf/ours_pileup_smallblock.tsv
              benchmarks/tomato1/results/perf/ours_pileup_smallblock/N<nn>/<sample>.psp

Env overrides:
  POP_VAR_CALLER_BIN  binary (default: $PROJECT_ROOT/target/release/pop_var_caller)
  REFERENCE           SL4.0 fasta (.fai sibling required)
  THREADS             pileup concurrency (default: 4)
  SIZES               comma-separated sample sizes (default: 1,2,4,8,12,16,20,24,26)

Invoke:
  uv run --script benchmarks/tomato1/scripts/perf_ours_pileup_smallblock.py
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

CALLER = "ours_pileup_smallblock"
BLOCK_TARGET_BYTES = "262144"  # 256 KiB; 1/4 of the 1 MiB pileup default
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
    print(f"block-target-bytes: {BLOCK_TARGET_BYTES}")
    print()
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
                "--threads", str(DEFAULT_THREADS),
                "--block-target-bytes", BLOCK_TARGET_BYTES,
                str(cram),
            ])
        # One caller process at a time, each internally threaded
        # (--threads N) — never several pileup processes in parallel
        # (contention + summed RSS distorts the per-caller measurement).
        print(f"[{CALLER}] N={n:>2}: {len(cmds)} pileups, sequential @ --threads {DEFAULT_THREADS}")
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
