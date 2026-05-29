# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Head-to-head comparison of `pop_var_caller var-calling` between
two pre-built binaries — typically `main` vs the working branch.

Reads the same PSP cohort (`results/ours/cohort/psp/*.psp`) for
every N in `DEFAULT_SIZES`, runs both binaries against an identical
subset, and records wall time, peak process-tree RSS, exit code,
and an md5 of the VCF body (the volatile `##source=` /
`##commandline=` header lines are stripped before hashing so a true
algorithm-level divergence is what shows up).

Output: a single TSV at
`results/perf/main_vs_branch.tsv` with one row per (caller, N).
Both binaries' VCFs land under
`results/perf/main_vs_branch/<caller>/N{NN}.vcf` so they can be
diff'd by hand if the md5s diverge.

Inputs (env-overridable):
  BIN_MAIN=/path/to/pop_var_caller_main
  BIN_BRANCH=/path/to/pop_var_caller_branch
  SIZES=1,2,4,...
  THREADS=4

Default binary paths:
  BIN_MAIN  defaults to  <repo_root>/tmp/bin/pop_var_caller.main
  BIN_BRANCH defaults to <repo_root>/target/release/pop_var_caller
"""

import hashlib
import os
import sys
from pathlib import Path

from perf_common import (
    DEFAULT_REFERENCE,
    DEFAULT_THREADS,
    PERF_DIR,
    PROJECT_ROOT,
    PSP_DIR,
    Measurement,
    banner,
    check_exists,
    list_inputs,
    measure,
    pick_subset,
    sizes_from_env,
)


BIN_MAIN = Path(
    os.environ.get("BIN_MAIN", PROJECT_ROOT / "tmp" / "bin" / "pop_var_caller.main")
)
BIN_BRANCH = Path(
    os.environ.get(
        "BIN_BRANCH", PROJECT_ROOT / "target" / "release" / "pop_var_caller"
    )
)

OUT_DIR = PERF_DIR / "main_vs_branch"
TSV_PATH = PERF_DIR / "main_vs_branch.tsv"


def vcf_body_md5(vcf_path: Path) -> str:
    """md5 of the VCF body, ignoring `##source=` and `##commandline=`
    (both carry the binary's version + invocation path, which differ
    between callers by construction). Returns "" when the file is
    missing — the run failed and the dashboard surfaces the exit
    code instead.
    """
    if not vcf_path.exists():
        return ""
    h = hashlib.md5()
    with vcf_path.open("rb") as fh:
        for line in fh:
            if line.startswith(b"##source=") or line.startswith(b"##commandline="):
                continue
            h.update(line)
    return h.hexdigest()


def run_one(
    binary: Path,
    caller_label: str,
    n: int,
    subset: list[Path],
    out_dir: Path,
) -> tuple[Measurement, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_vcf = out_dir / f"N{n:02d}.vcf"
    cmd = [
        str(binary), "var-calling",
        "--reference", str(DEFAULT_REFERENCE),
        "--output", str(out_vcf),
        "--threads", str(DEFAULT_THREADS),
        *[str(p) for p in subset],
    ]
    print(f"[{caller_label}] N={n:>2}: var-calling over {len(subset)} PSPs")
    wall, peak_bytes, exit_code = measure(cmd)
    peak_mb = peak_bytes / 1024 / 1024
    md5 = vcf_body_md5(out_vcf) if exit_code == 0 else ""
    measurement = Measurement(caller_label, n, wall, peak_mb, exit_code)
    status = "ok" if exit_code == 0 else f"FAILED (exit {exit_code})"
    print(f"  -> {wall:.1f}s, peak {peak_mb:.0f} MB, md5={md5[:8] or '----'}  [{status}]")
    return measurement, md5


def write_compare_tsv(
    rows: list[tuple[Measurement, str]],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(
            "caller\tn_samples\twall_seconds\tpeak_rss_mb\texit_code\tvcf_md5\n"
        )
        for m, md5 in rows:
            fh.write(
                f"{m.caller}\t{m.n_samples}\t{m.wall_seconds:.3f}\t"
                f"{m.peak_rss_mb:.1f}\t{m.exit_code}\t{md5}\n"
            )


def print_pairwise_summary(rows: list[tuple[Measurement, str]]) -> int:
    """Group by N, print main/branch deltas + md5 match. Returns the
    number of size points where the md5s disagree (callers > 0 ⇒ the
    operator should investigate before reading the timing chart)."""
    by_n: dict[int, dict[str, tuple[Measurement, str]]] = {}
    for m, md5 in rows:
        by_n.setdefault(m.n_samples, {})[m.caller] = (m, md5)
    print()
    print(
        f"{'N':>3}  "
        f"{'main wall':>10}  {'branch wall':>11}  {'speedup':>8}  "
        f"{'main RSS':>9}  {'branch RSS':>10}  {'rss ratio':>9}  "
        f"{'md5 match':>9}"
    )
    print("-" * 84)
    mismatch_count = 0
    for n in sorted(by_n):
        slot = by_n[n]
        main_pair = slot.get("main")
        branch_pair = slot.get("branch")
        if main_pair is None or branch_pair is None:
            print(f"{n:>3}  (missing data for one caller)")
            continue
        m_m, m_md5 = main_pair
        b_m, b_md5 = branch_pair
        speedup = m_m.wall_seconds / b_m.wall_seconds if b_m.wall_seconds else 0.0
        rss_ratio = b_m.peak_rss_mb / m_m.peak_rss_mb if m_m.peak_rss_mb else 0.0
        md5_match = "yes" if m_md5 and b_md5 and m_md5 == b_md5 else "no"
        if md5_match != "yes":
            mismatch_count += 1
        print(
            f"{n:>3}  "
            f"{m_m.wall_seconds:>9.1f}s  {b_m.wall_seconds:>10.1f}s  "
            f"{speedup:>7.2f}x  "
            f"{m_m.peak_rss_mb:>7.0f}MB  {b_m.peak_rss_mb:>8.0f}MB  "
            f"{rss_ratio:>8.2f}x  {md5_match:>9}"
        )
    return mismatch_count


def main() -> int:
    sizes = sizes_from_env()
    check_exists(
        BIN_MAIN, BIN_BRANCH, DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai")
    )
    psps = list_inputs(PSP_DIR, ".psp")

    print(f"=== main vs branch perf sweep ===")
    print(f"main bin:   {BIN_MAIN}")
    print(f"branch bin: {BIN_BRANCH}")
    print(f"sizes:      {sizes}")
    print(f"threads:    {DEFAULT_THREADS}")
    print(f"reference:  {DEFAULT_REFERENCE}")
    print()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[tuple[Measurement, str]] = []
    for n in sizes:
        subset = pick_subset(psps, n)
        main_pair = run_one(BIN_MAIN, "main", n, subset, OUT_DIR / "main")
        branch_pair = run_one(BIN_BRANCH, "branch", n, subset, OUT_DIR / "branch")
        rows.append(main_pair)
        rows.append(branch_pair)

    write_compare_tsv(rows, TSV_PATH)
    print(f"\nwrote {TSV_PATH}")

    mismatches = print_pairwise_summary(rows)
    if mismatches:
        print(
            f"\nWARNING: VCF md5 differs at {mismatches} of {len(set(r[0].n_samples for r in rows))} "
            f"size point(s). Investigate before trusting the timing numbers — "
            f"the two callers may not be computing the same answer."
        )
        return 2
    print("\nAll VCFs byte-identical (header-stripped). Timing comparison is valid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
