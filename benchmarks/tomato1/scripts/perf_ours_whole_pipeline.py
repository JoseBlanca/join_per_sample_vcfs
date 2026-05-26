# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Scaling perf experiment for the full pop_var_caller cohort
pipeline starting from CRAMs (pileup per sample -> joint var-calling
over PSPs).

`var-calling-from-bam` is single-sample only, so it can't be the
cohort analogue of freebayes' one-shot CRAMs -> VCF. The fair
counterpart for cohort scaling is the full pipeline you'd actually
run: per-sample pileup followed by joint var-calling. That's what
this script measures.

THREADS is treated as a *parallelism budget*: stage 1 spawns up to
THREADS pileup processes concurrently (each single-threaded — the
budget is exhausted by the parallel processes); stage 2 runs the
joint var-calling with --threads THREADS. Pileups write to a
per-N scratch directory so each measurement starts cold.

Total wall = stage1_wall + stage2_wall.
Peak RSS  = max(stage1_peak_sum_across_alive_pileups, stage2_peak).

The (perf_ours_whole_pipeline, perf_freebayes) pair shows how the two
"full cohort from BAM" workflows scale relative to each other; the
(perf_ours_joint, perf_gatk_joint) pair compares only the joint
genotyping stage (intermediates assumed pre-built).

Sample sizes: 1, 2, 4, 8, 12, 16, 18 (full tomato1 cohort).
Inputs:       benchmarks/tomato1/crams/*.bench.cram
Output:       benchmarks/tomato1/results/perf/ours_whole_pipeline.tsv
              benchmarks/tomato1/results/perf/ours_whole_pipeline/N<nn>/<sample>.psp
              benchmarks/tomato1/results/perf/ours_whole_pipeline/N<nn>.vcf

Env overrides:
  POP_VAR_CALLER_BIN  binary (default: $PROJECT_ROOT/target/release/pop_var_caller)
  REFERENCE           SL4.0 fasta (.fai sibling required)
  THREADS             parallelism budget (default: 4 — drives pileup
                      concurrency in stage 1 and --threads in stage 2)
  SIZES               comma-separated sample sizes (default: 1,2,4,8,12,16,18)

Invoke:
  uv run --script benchmarks/tomato1/scripts/perf_ours_whole_pipeline.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    CRAM_DIR, DEFAULT_REFERENCE, DEFAULT_THREADS, PERF_DIR, PROJECT_ROOT,
    Measurement, banner, check_exists, list_inputs, measure, measure_pool,
    pick_subset, sizes_from_env, write_tsv,
)

CALLER = "ours_whole_pipeline"
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
        out_vcf = OUT_DIR / f"N{n:02d}.vcf"

        # --- Stage 1: parallel per-sample pileup ---
        # Each pileup is single-threaded; we run up to THREADS of them
        # at a time. This exhausts the parallelism budget across
        # processes rather than within one.
        pileup_cmds = []
        psp_paths = []
        for cram in subset:
            base = cram.name.removesuffix(".bench.cram")
            psp = n_dir / f"{base}.psp"
            psp_paths.append(psp)
            pileup_cmds.append([
                str(BIN), "pileup",
                "--reference", str(DEFAULT_REFERENCE),
                "--output", str(psp),
                "--threads", "1",
                str(cram),
            ])
        concurrency = min(DEFAULT_THREADS, len(pileup_cmds))
        print(f"[{CALLER}] N={n:>2}: stage 1 — {len(pileup_cmds)} pileups @ concurrency={concurrency}")
        s1_wall, s1_peak, s1_ec = measure_pool(pileup_cmds, max_concurrent=concurrency)
        s1_peak_mb = s1_peak / 1024 / 1024
        print(f"  stage 1 -> {s1_wall:.1f}s, peak {s1_peak_mb:.0f} MB (sum across alive pileups)")

        if s1_ec != 0:
            rows.append(Measurement(CALLER, n, s1_wall, s1_peak_mb, s1_ec))
            print(f"  pileup FAILED (exit {s1_ec}); skipping stage 2 for this N")
            continue

        # --- Stage 2: joint var-calling over the freshly-built PSPs ---
        varcall_cmd = [
            str(BIN), "var-calling",
            "--reference", str(DEFAULT_REFERENCE),
            "--output", str(out_vcf),
            "--threads", str(DEFAULT_THREADS),
            *[str(p) for p in psp_paths],
        ]
        print(f"[{CALLER}] N={n:>2}: stage 2 — joint var-calling over {len(psp_paths)} PSPs")
        s2_wall, s2_peak, s2_ec = measure(varcall_cmd)
        s2_peak_mb = s2_peak / 1024 / 1024
        print(f"  stage 2 -> {s2_wall:.1f}s, peak {s2_peak_mb:.0f} MB")

        total_wall = s1_wall + s2_wall
        total_peak_mb = max(s1_peak_mb, s2_peak_mb)
        rows.append(Measurement(CALLER, n, total_wall, total_peak_mb, s2_ec))
        status = "ok" if s2_ec == 0 else f"FAILED (exit {s2_ec})"
        print(f"  TOTAL   -> {total_wall:.1f}s, peak {total_peak_mb:.0f} MB  [{status}]")

    tsv = PERF_DIR / f"{CALLER}.tsv"
    write_tsv(rows, tsv)
    print(f"\nwrote {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
