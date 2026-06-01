# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Thread-scaling perf for the cohort PSP -> VCF caller (`var-calling`),
at the MAXIMUM cohort size (all pre-built per-sample PSPs).

This is the cohort-joint model: ONE `var-calling` process, parallelised
internally via `--threads`. We sweep the thread count and record wall +
peak RSS, to show how the columnar produce / DUST worker pool / parallel
block-consume architecture scales with threads at the largest N.

Single process per measurement (no inter-process contention), so wall and
RSS characterise the caller itself.

Sweep:   THREAD_SWEEP env (comma-separated; default 1,2,3,4,6,8).
Inputs:  benchmarks/tomato1/results/ours/cohort/psp/*.psp (the max cohort).
Output:  benchmarks/tomato1/results/perf/ours_joint_threads.tsv
         columns: threads, wall_seconds, peak_rss_mb, exit_code
         benchmarks/tomato1/results/perf/ours_joint_threads/t<NN>.vcf

Env overrides:
  POP_VAR_CALLER_BIN  binary (default: $PROJECT_ROOT/target/release/pop_var_caller)
  REFERENCE           SL4.0 fasta (.fai sibling required)
  THREAD_SWEEP        comma-separated thread counts (default: 1,2,3,4,6,8)

Invoke:
  uv run --script benchmarks/tomato1/scripts/perf_ours_joint_threads.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    DEFAULT_REFERENCE, PERF_DIR, PROJECT_ROOT, PSP_DIR,
    check_exists, list_inputs, measure,
)

BIN = Path(os.environ.get(
    "POP_VAR_CALLER_BIN",
    str(PROJECT_ROOT / "target" / "release" / "pop_var_caller"),
))
THREAD_SWEEP = [
    int(t) for t in os.environ.get("THREAD_SWEEP", "1,2,3,4,6,8").split(",") if t
]
OUT_DIR = PERF_DIR / "ours_joint_threads"


def main() -> int:
    check_exists(BIN, DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai"))
    psps = list_inputs(PSP_DIR, ".psp")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(
        f"=== thread scaling: var-calling over {len(psps)} PSPs "
        f"(max cohort), threads={THREAD_SWEEP} ==="
    )

    rows: list[tuple[int, float, float, int]] = []
    for t in THREAD_SWEEP:
        out_vcf = OUT_DIR / f"t{t:02d}.vcf"
        cmd = [
            str(BIN), "var-calling",
            "--reference", str(DEFAULT_REFERENCE),
            "--output", str(out_vcf),
            "--threads", str(t),
            *[str(p) for p in psps],
        ]
        print(f"[threads={t:>2}] var-calling over {len(psps)} PSPs")
        wall, peak_bytes, exit_code = measure(cmd)
        peak_mb = peak_bytes / 1024 / 1024
        rows.append((t, wall, peak_mb, exit_code))
        status = "ok" if exit_code == 0 else f"FAILED (exit {exit_code})"
        print(f"  -> {wall:.1f}s, peak {peak_mb:.0f} MB  [{status}]")

    tsv = PERF_DIR / "ours_joint_threads.tsv"
    with tsv.open("w") as fh:
        fh.write("threads\twall_seconds\tpeak_rss_mb\texit_code\n")
        for t, wall, peak_mb, exit_code in rows:
            fh.write(f"{t}\t{wall:.3f}\t{peak_mb:.1f}\t{exit_code}\n")
    print(f"\nwrote {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
