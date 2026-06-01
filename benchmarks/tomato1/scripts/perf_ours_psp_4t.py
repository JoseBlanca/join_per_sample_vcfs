# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Section-2 perf: build ONE sample's per-sample intermediate (the
`.psp`) with the full thread budget — `pop_var_caller pileup --threads
THREADS` on a single CRAM.

This is the our-caller half of the "create one per-sample intermediate
file, one sample, THREADS threads" bar-chart comparison; the GATK half
is `perf_gatk_gvcf_4t.py` (HaplotypeCaller -ERC GVCF). Unlike
`perf_ours_pileup.py` (which runs N single-threaded pileups
concurrently for the cohort-scaling panel), here a single pileup
process gets all THREADS — the realistic "I have one new sample, build
its intermediate as fast as I can" operation.

Single process, so wall and peak RSS characterise the pileup itself.

Inputs:  benchmarks/tomato1/crams/*.bench.cram (first sample by name)
Output:  benchmarks/tomato1/results/perf/ours_psp_4t.tsv
         benchmarks/tomato1/results/perf/ours_psp_4t/<sample>.psp

Env overrides:
  POP_VAR_CALLER_BIN  binary (default: $PROJECT_ROOT/target/release/pop_var_caller)
  REFERENCE           SL4.0 fasta (.fai sibling required)
  THREADS             pileup --threads (default: 4)

Invoke:
  uv run --script benchmarks/tomato1/scripts/perf_ours_psp_4t.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    CRAM_DIR, DEFAULT_REFERENCE, DEFAULT_THREADS, PERF_DIR, PROJECT_ROOT,
    Measurement, check_exists, list_inputs, measure, write_tsv,
)

CALLER = "ours_psp_4t"
BIN = Path(os.environ.get(
    "POP_VAR_CALLER_BIN",
    str(PROJECT_ROOT / "target" / "release" / "pop_var_caller"),
))
OUT_DIR = PERF_DIR / CALLER


def main() -> int:
    check_exists(BIN, DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai"))
    crams = list_inputs(CRAM_DIR, ".bench.cram")
    cram = crams[0]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"=== {CALLER}: pileup one sample @ --threads {DEFAULT_THREADS} ===")
    print(f"sample:    {cram.name}")
    print(f"reference: {DEFAULT_REFERENCE}")
    print()

    base = cram.name.removesuffix(".bench.cram")
    psp = OUT_DIR / f"{base}.psp"
    cmd = [
        str(BIN), "pileup",
        "--reference", str(DEFAULT_REFERENCE),
        "--output", str(psp),
        "--threads", str(DEFAULT_THREADS),
        str(cram),
    ]
    print(f"[{CALLER}] pileup {cram.name} -> {psp.name}")
    wall, peak_bytes, exit_code = measure(cmd)
    peak_mb = peak_bytes / 1024 / 1024
    status = "ok" if exit_code == 0 else f"FAILED (exit {exit_code})"
    print(f"  -> {wall:.1f}s, peak {peak_mb:.0f} MB  [{status}]")

    rows = [Measurement(CALLER, 1, wall, peak_mb, exit_code)]
    tsv = PERF_DIR / f"{CALLER}.tsv"
    write_tsv(rows, tsv)
    print(f"\nwrote {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
