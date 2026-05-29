# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Scaling perf experiment for `pop_var_caller var-calling` over PSPs
that were built with `--block-target-bytes 262144` (256 KiB, 1/4 of
the 1 MiB default).

Companion to perf_ours_joint.py. The cohort reader holds one decoded
block per `n_threads × N` worker; shrinking the block target should
drop peak heap roughly 4x relative to the default-block run, with
slightly worse zstd compression (so larger .psp on disk). This is
the scenario where the block-size change is expected to matter most.

PSP inputs must already exist at results/ours/cohort/psp_smallblock/
— a separate cohort PSP set built with the small block target so it
doesn't disturb the default-block PSPs that perf_ours_joint reads.
Build them with:

  ./benchmarks/tomato1/scripts/build_smallblock_psps.sh

Sample sizes: 1, 2, 4, 8, 12, 16, 20, 24, 26 (full tomato1 cohort).
Inputs:       benchmarks/tomato1/results/ours/cohort/psp_smallblock/*.psp
Output:       benchmarks/tomato1/results/perf/ours_joint_smallblock.tsv
              benchmarks/tomato1/results/perf/ours_joint_smallblock/N<nn>.vcf

Env overrides:
  POP_VAR_CALLER_BIN  binary (default: $PROJECT_ROOT/target/release/pop_var_caller)
  REFERENCE           SL4.0 fasta (.fai sibling required)
  THREADS             rayon worker count (default: 4)
  SIZES               comma-separated sample sizes (default: 1,2,4,8,12,16,20,24,26)

Invoke:
  uv run --script benchmarks/tomato1/scripts/perf_ours_joint_smallblock.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    DEFAULT_REFERENCE, DEFAULT_THREADS, PERF_DIR, PROJECT_ROOT,
    PSP_SMALLBLOCK_DIR,
    Measurement, banner, check_exists, list_inputs, measure, pick_subset,
    sizes_from_env, write_tsv,
)

CALLER = "ours_joint_smallblock"
BIN = Path(os.environ.get(
    "POP_VAR_CALLER_BIN",
    str(PROJECT_ROOT / "target" / "release" / "pop_var_caller"),
))
OUT_DIR = PERF_DIR / CALLER


def main() -> int:
    sizes = sizes_from_env()
    check_exists(BIN, DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai"))
    psps = list_inputs(PSP_SMALLBLOCK_DIR, ".psp")
    banner(CALLER, sizes)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[Measurement] = []
    for n in sizes:
        subset = pick_subset(psps, n)
        out_vcf = OUT_DIR / f"N{n:02d}.vcf"
        cmd = [
            str(BIN), "var-calling",
            "--reference", str(DEFAULT_REFERENCE),
            "--output", str(out_vcf),
            "--threads", str(DEFAULT_THREADS),
            *[str(p) for p in subset],
        ]
        print(f"[{CALLER}] N={n:>2}: var-calling over {len(subset)} smallblock PSPs")
        wall, peak_bytes, exit_code = measure(cmd)
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
