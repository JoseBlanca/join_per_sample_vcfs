# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""2-D memory-scaling sweep for the cohort PSP -> VCF caller (`var-calling`):
peak RSS and wall as a function of BOTH `--threads` and cohort size N.

Motivation. The streaming / bounded-queue architecture predicts that peak
RSS should

- grow **~linearly with `--threads`** — the producer hand-off queues are
  sized `≈ 2 × threads`, and each caller thread holds one in-flight chunk,
  so the resident working set is proportional to the number of workers; and
- **plateau with N (cohort size)** — the producer bounds the in-flight
  working set by `target_variants_per_chunk`, i.e. it never materialises
  all-samples data at once, so per-chunk memory grows only with the columns
  it must carry, not unboundedly with N.

This driver produces the data for the two dashboard panels that test those
two claims directly (RSS-vs-N per thread count; RSS-vs-threads per N). One
`var-calling` process per measurement (cohort-joint model) — single process,
so wall + RSS characterise the caller itself with no inter-process noise.

Sweeps (Cartesian product):
  SAMPLE_SWEEP env (comma-separated; default 1,2,4,8,16,32,50)
  THREAD_SWEEP env (comma-separated; default 1,2,3,4,5)

Inputs:  $PSP_DIR/*.psp  (default: perf_common.PSP_DIR). Override PSP_DIR for
         an ad-hoc cohort, e.g. the real aligned per-sample PSPs.
Output:  benchmarks/tomato1/results/perf/ours_joint_mem_scaling.tsv
         columns: threads, n_samples, wall_seconds, peak_rss_mb, exit_code
         benchmarks/tomato1/results/perf/ours_joint_mem_scaling/n<NNN>_t<NN>.vcf

Env overrides:
  POP_VAR_CALLER_BIN  binary (default: $PROJECT_ROOT/target/release/pop_var_caller)
  REFERENCE           SL4.0 fasta (.fai sibling required) — pass an EXPANDED
                      path (the default contains a literal $HOME).
  PSP_DIR             directory of per-sample .psp inputs
  SAMPLE_SWEEP        comma-separated cohort sizes (default 1,2,4,8,16,32,50)
  THREAD_SWEEP        comma-separated thread counts (default 1,2,3,4,5)

Invoke:
  REFERENCE=$HOME/genomes/.../S_lycopersicum_chromosomes.4.00.fa \
  PSP_DIR=/path/to/psp_dir \
  uv run --script benchmarks/tomato1/scripts/perf_ours_joint_mem_scaling.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    DEFAULT_REFERENCE, PERF_DIR, PROJECT_ROOT, PSP_DIR,
    check_exists, list_inputs, measure, pick_subset,
)

BIN = Path(os.environ.get(
    "POP_VAR_CALLER_BIN",
    str(PROJECT_ROOT / "target" / "release" / "pop_var_caller"),
))
PSP_INPUT_DIR = Path(os.environ.get("PSP_DIR", str(PSP_DIR)))
SAMPLE_SWEEP = [
    int(n) for n in os.environ.get("SAMPLE_SWEEP", "1,2,4,8,16,32,50").split(",") if n
]
THREAD_SWEEP = [
    int(t) for t in os.environ.get("THREAD_SWEEP", "1,2,3,4,5").split(",") if t
]
OUT_DIR = PERF_DIR / "ours_joint_mem_scaling"


def main() -> int:
    check_exists(BIN, DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai"))
    psps = list_inputs(PSP_INPUT_DIR, ".psp")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sizes = [n for n in SAMPLE_SWEEP if n <= len(psps)]
    if not sizes:
        print(f"no usable cohort sizes in {SAMPLE_SWEEP}; "
              f"have only {len(psps)} PSPs in {PSP_INPUT_DIR}")
        return 1
    print(
        f"=== mem scaling: var-calling, threads={THREAD_SWEEP} x N={sizes} "
        f"over {len(psps)} PSPs in {PSP_INPUT_DIR} ==="
    )

    rows: list[tuple[int, int, float, float, int]] = []
    # Outer loop on N (subset is reused across the thread sweep), inner on
    # threads, so each TSV block is one N across all thread counts.
    for n in sizes:
        subset = pick_subset(psps, n)
        for t in THREAD_SWEEP:
            out_vcf = OUT_DIR / f"n{n:03d}_t{t:02d}.vcf"
            cmd = [
                str(BIN), "var-calling",
                "--reference", str(DEFAULT_REFERENCE),
                "--output", str(out_vcf),
                "--threads", str(t),
                *[str(p) for p in subset],
            ]
            print(f"[N={n:>3} threads={t:>2}] var-calling over {n} PSPs")
            wall, peak_bytes, exit_code = measure(cmd)
            peak_mb = peak_bytes / 1024 / 1024
            rows.append((t, n, wall, peak_mb, exit_code))
            status = "ok" if exit_code == 0 else f"FAILED (exit {exit_code})"
            print(f"  -> {wall:.1f}s, peak {peak_mb:.0f} MB  [{status}]")

    tsv = PERF_DIR / "ours_joint_mem_scaling.tsv"
    with tsv.open("w") as fh:
        fh.write("threads\tn_samples\twall_seconds\tpeak_rss_mb\texit_code\n")
        for t, n, wall, peak_mb, exit_code in rows:
            fh.write(f"{t}\t{n}\t{wall:.3f}\t{peak_mb:.1f}\t{exit_code}\n")
    print(f"\nwrote {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
