"""Shared helpers for the per-caller scaling experiments.

Each `perf_<caller>.py` is a uv PEP-723 script that imports from this
module to:
  - measure(cmd) -> (wall_s, peak_rss_bytes, exit_code) via psutil
    polling, summed across the process tree (rayon workers, JVM
    child processes, etc).
  - pick a deterministic N-sample subset of the tomato1 input pool
    (first N alphabetically — same samples chosen at each N across
    callers, so cross-caller numbers are over the same data).
  - emit a TSV row per (caller, N) into results/perf/<caller>.tsv.

This module is NOT itself a PEP-723 script; psutil is declared in
each caller script's header so `uv run --script` resolves it from
the script's ephemeral env. Sibling imports work because uv
preserves Python's default "script-dir on sys.path[0]" behaviour.
"""

import dataclasses
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Sequence

import psutil


SCRIPT_DIR = Path(__file__).resolve().parent
TEST_DIR = SCRIPT_DIR.parent           # benchmarks/tomato1/
PROJECT_ROOT = TEST_DIR.parent.parent

DEFAULT_REFERENCE = Path(os.environ.get(
    "REFERENCE",
    str(Path.home() / "genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa"),
))
DEFAULT_BED = TEST_DIR / "regions.bed"
DEFAULT_SIZES = [1, 2, 4, 8, 12, 16, 18]
DEFAULT_THREADS = int(os.environ.get("THREADS", "4"))

PERF_DIR = TEST_DIR / "results" / "perf"
CRAM_DIR = TEST_DIR / "crams"
PSP_DIR = TEST_DIR / "results" / "ours" / "cohort" / "psp"
GVCF_DIR = TEST_DIR / "results" / "gatk" / "cohort" / "gvcf"


@dataclasses.dataclass
class Measurement:
    caller: str
    n_samples: int
    wall_seconds: float
    peak_rss_mb: float
    exit_code: int


def measure(
    cmd: Sequence[str],
    *,
    poll_interval: float = 0.05,
    stdout=None,
    env: dict | None = None,
) -> tuple[float, int, int]:
    """Run `cmd`, polling RSS across the process tree at
    `poll_interval` seconds; return (wall_s, peak_rss_bytes,
    exit_code).

    Summing across descendants is essential — rayon-parallelised
    callers (pop_var_caller) fork workers; GATK's launcher process
    spawns the JVM. The parent's own RSS undercounts in both cases.
    """
    t0 = time.monotonic()
    proc = subprocess.Popen(cmd, stdout=stdout, env=env)
    try:
        parent = psutil.Process(proc.pid)
    except psutil.NoSuchProcess:
        # Process exited sub-tick — record what we can and move on.
        proc.wait()
        return time.monotonic() - t0, 0, proc.returncode

    peak_bytes = 0
    while proc.poll() is None:
        try:
            rss = parent.memory_info().rss
            for child in parent.children(recursive=True):
                try:
                    rss += child.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            peak_bytes = max(peak_bytes, rss)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        time.sleep(poll_interval)
    proc.wait()
    return time.monotonic() - t0, peak_bytes, proc.returncode


def measure_pool(
    cmds: Sequence[Sequence[str]],
    *,
    max_concurrent: int,
    poll_interval: float = 0.05,
) -> tuple[float, int, int]:
    """Run a pool of commands with up to `max_concurrent` running at
    once; return (wall_s, peak_rss_bytes, exit_code).

    Wall is end-to-end (first launch to last exit), NOT the sum of
    individual wall times — this is what a parallel-processes
    scaling experiment cares about. Peak RSS is the maximum, over
    all polling ticks, of (sum of RSS across all currently-running
    processes and their descendants). Exit code is 0 iff every
    process succeeded, else the first non-zero exit observed.

    `max_concurrent` of 1 degenerates to a sequential run.
    """
    t0 = time.monotonic()
    pending = list(cmds)
    running: list[tuple[subprocess.Popen, psutil.Process | None]] = []
    first_failure = 0
    peak_bytes = 0

    while pending or running:
        # Launch new processes up to the concurrency cap.
        while pending and len(running) < max_concurrent:
            cmd = pending.pop(0)
            proc = subprocess.Popen(cmd)
            try:
                p = psutil.Process(proc.pid)
            except psutil.NoSuchProcess:
                p = None
            running.append((proc, p))

        # Poll all currently-running processes.
        total_rss = 0
        still_running = []
        for proc, p in running:
            ec = proc.poll()
            if ec is not None:
                if ec != 0 and first_failure == 0:
                    first_failure = ec
                proc.wait()
                continue
            still_running.append((proc, p))
            if p is None:
                continue
            try:
                rss = p.memory_info().rss
                for child in p.children(recursive=True):
                    try:
                        rss += child.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                total_rss += rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
        peak_bytes = max(peak_bytes, total_rss)
        running = still_running

        if running or pending:
            time.sleep(poll_interval)

    return time.monotonic() - t0, peak_bytes, first_failure


def list_inputs(directory: Path, suffix: str) -> list[Path]:
    """Sorted list of *<suffix> in `directory`; bail with a helpful
    message if empty. Sorted order gives deterministic N-sample
    subsets across runs and across callers."""
    files = sorted(directory.glob(f"*{suffix}"))
    if not files:
        sys.exit(f"no *{suffix} files in {directory}")
    return files


def pick_subset(items: list[Path], n: int) -> list[Path]:
    if n > len(items):
        sys.exit(f"requested N={n} but only {len(items)} inputs available")
    return items[:n]


def check_exists(*paths: Path) -> None:
    missing = [p for p in paths if not Path(p).exists()]
    if missing:
        sys.exit("missing required input(s):\n" + "\n".join(f"  {p}" for p in missing))


def write_tsv(rows: list[Measurement], path: Path) -> None:
    """Overwrites — a clean (caller, N) -> measurement table per run.
    Re-running is the usual workflow (you want fresh numbers anyway,
    not a longitudinal accumulation)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("caller\tn_samples\twall_seconds\tpeak_rss_mb\texit_code\n")
        for m in rows:
            fh.write(
                f"{m.caller}\t{m.n_samples}\t{m.wall_seconds:.3f}\t"
                f"{m.peak_rss_mb:.1f}\t{m.exit_code}\n"
            )


def banner(caller: str, sizes: list[int]) -> None:
    print(f"=== {caller} perf sweep ===")
    print(f"sizes:     {sizes}")
    print(f"reference: {DEFAULT_REFERENCE}")
    print(f"threads:   {DEFAULT_THREADS}")
    print()


def sizes_from_env() -> list[int]:
    """Parse SIZES env override (comma-separated). Falls back to
    DEFAULT_SIZES; clamps requested values against the input pool
    happens later in pick_subset."""
    raw = os.environ.get("SIZES")
    if not raw:
        return DEFAULT_SIZES
    return [int(s) for s in raw.split(",") if s.strip()]
