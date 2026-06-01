# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""§3 stage-1 measurement: build the cohort's per-sample `.psp` files
ONCE through a FIFO pool of THREADS workers, and record the build-time
*curve* — how long it takes to have the first N samples' PSPs ready —
plus the (≈constant) peak RSS of the parallel build.

Why one build instead of rebuilding per N: a sample's `.psp` is
independent of cohort size (pileup is per-sample), so the PSP for sample
X is identical whether the cohort is 1 or 50. And with a FIFO pool, a
freed worker always starts the lowest-index pending job, so the schedule
of the first N jobs is independent of jobs > N. Therefore

    makespan(first N) = max(completion_time of samples 0..N-1)

measured in a SINGLE build of all max(SIZES) samples is *exactly* the
wall you'd see building only those N. That gives the whole §3 stage-1
curve from one pass (the old perf_ours_whole_pipeline rebuilt PSPs at
every N — ~Σ N redundant pileups).

The PSPs are written to the canonical cohort dir
(results/ours/cohort/psp), so this run also produces the inputs that
perf_ours_joint (§4) and perf_ours_joint_threads (§5) consume — it is
the cohort PSP build. The dir is cleared first for a cold, honest timing.

§3 ours CRAM→VCF total(N) = THIS makespan(N) + perf_ours_joint wall(N);
peak(N) = max(this build peak, joint peak(N)). The dashboard sums them.

Peak RSS here is the steady-state footprint of THREADS concurrent
single-threaded pileups — a fixed parallelism budget, so it does not
grow with N (that bounded footprint is the point). For N < THREADS the
real build would use fewer concurrent processes (lower peak); the row
still reports the full-build peak, which the dashboard treats as the
budget's footprint.

Inputs:  benchmarks/tomato1/crams/*.bench.cram (first max(SIZES), sorted)
Output:  benchmarks/tomato1/results/perf/ours_pileup_build.tsv
         benchmarks/tomato1/results/ours/cohort/psp/<sample>.psp  (canonical)

Env overrides:
  POP_VAR_CALLER_BIN  binary (default: $PROJECT_ROOT/target/release/pop_var_caller)
  REFERENCE           SL4.0 fasta (.fai sibling required)
  THREADS             parallel pileup workers (default: 4)
  SIZES               comma-separated sample sizes (default: see perf_common)

Invoke:
  uv run --script benchmarks/tomato1/scripts/perf_ours_pileup_build.py
"""

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    CRAM_DIR, DEFAULT_REFERENCE, DEFAULT_THREADS, PERF_DIR, PROJECT_ROOT, PSP_DIR,
    Measurement, banner, check_exists, list_inputs, pick_subset,
    sizes_from_env, write_tsv,
)

import psutil  # noqa: E402

CALLER = "ours_pileup_build"
BIN = Path(os.environ.get(
    "POP_VAR_CALLER_BIN",
    str(PROJECT_ROOT / "target" / "release" / "pop_var_caller"),
))


def run_fifo_pool_timed(cmds, logs, max_concurrent, *, poll_interval=0.05):
    """Run `cmds` (in order) through a FIFO pool of `max_concurrent`
    workers. A freed worker always starts the lowest-index pending job.
    Returns (completion_times, exit_codes, peak_bytes):
      completion_times[i] = wall-seconds since start at which cmds[i] ended,
      exit_codes[i]        = its exit code,
      peak_bytes           = max over polling ticks of summed RSS across
                             all live workers (and descendants)."""
    t0 = time.monotonic()
    n = len(cmds)
    completion = [None] * n
    exits = [None] * n
    nxt = 0
    running = []  # (idx, proc, psutil_proc, log_fh)
    peak = 0
    while nxt < n or running:
        while nxt < n and len(running) < max_concurrent:
            fh = open(logs[nxt], "w")
            proc = subprocess.Popen(cmds[nxt], stdout=fh, stderr=subprocess.STDOUT)
            try:
                p = psutil.Process(proc.pid)
            except psutil.NoSuchProcess:
                p = None
            running.append((nxt, proc, p, fh))
            nxt += 1
        total = 0
        still = []
        for idx, proc, p, fh in running:
            ec = proc.poll()
            if ec is not None:
                completion[idx] = time.monotonic() - t0
                exits[idx] = ec
                proc.wait()
                fh.close()
                continue
            still.append((idx, proc, p, fh))
            if p is None:
                continue
            try:
                rss = p.memory_info().rss
                for child in p.children(recursive=True):
                    try:
                        rss += child.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                total += rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
        peak = max(peak, total)
        running = still
        if running or nxt < n:
            time.sleep(poll_interval)
    return completion, exits, peak


def main() -> int:
    sizes = sizes_from_env()
    check_exists(BIN, DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai"))
    crams = list_inputs(CRAM_DIR, ".bench.cram")
    n_max = max(sizes)
    subset = pick_subset(crams, n_max)  # first n_max, sorted — matches pick_subset
    banner(CALLER, sizes)
    print(f"cold build of {len(subset)} PSPs via FIFO pool of {DEFAULT_THREADS} workers")
    print(f"psp dir:   {PSP_DIR} (cleared for cold timing; canonical cohort PSPs)")
    print()

    # Cold build: clear the canonical PSP dir so the timing is from scratch.
    if PSP_DIR.exists():
        shutil.rmtree(PSP_DIR)
    PSP_DIR.mkdir(parents=True, exist_ok=True)

    cmds, logs, bases = [], [], []
    for cram in subset:
        base = cram.name.removesuffix(".bench.cram")
        bases.append(base)
        cmds.append([
            str(BIN), "pileup",
            "--reference", str(DEFAULT_REFERENCE),
            "--output", str(PSP_DIR / f"{base}.psp"),
            "--threads", "1",
            str(cram),
        ])
        logs.append(str(PSP_DIR / f"{base}.build.log"))

    completion, exits, peak_bytes = run_fifo_pool_timed(
        cmds, logs, max_concurrent=DEFAULT_THREADS,
    )
    peak_mb = peak_bytes / 1024 / 1024
    n_fail = sum(1 for e in exits if e not in (0, None))
    print(f"built {len(subset)} PSPs; peak {peak_mb:.0f} MB; {n_fail} failure(s)")

    # makespan(N) = max completion time among the first N samples (sorted).
    rows = []
    for n in sizes:
        first_n = completion[:n]
        first_n_exits = exits[:n]
        makespan = max(first_n) if first_n else 0.0
        ec = next((e for e in first_n_exits if e not in (0, None)), 0)
        rows.append(Measurement(CALLER, n, makespan, peak_mb, ec))
        print(f"  N={n:>2}: build makespan {makespan:.1f}s  (peak {peak_mb:.0f} MB, exit {ec})")

    tsv = PERF_DIR / f"{CALLER}.tsv"
    write_tsv(rows, tsv)
    print(f"\nwrote {tsv}")
    print(f"canonical PSPs left in {PSP_DIR} for perf_ours_joint / _threads")
    return 0


if __name__ == "__main__":
    sys.exit(main())
