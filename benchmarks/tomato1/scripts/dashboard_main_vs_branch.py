# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "matplotlib",
# ]
# ///

# Marimo dashboard — head-to-head plot of `pop_var_caller var-calling`
# performance between two binaries (typically `main` vs the working
# branch). Reads `results/perf/main_vs_branch.tsv` written by
# `perf_main_vs_branch.py` and visualises:
#
#   - wall time per N (both callers)
#   - peak process-tree RSS per N (both callers)
#   - branch-vs-main speedup ratio (main_wall / branch_wall)
#   - branch-vs-main memory ratio (branch_rss / main_rss)
#   - VCF byte-identity status per N (md5 match yes/no)
#
# A non-"yes" in the identity column means the two callers are
# producing different VCFs — the timing comparison is meaningless
# until that's investigated.
#
# Recommended invocation:
#
#   uvx marimo edit --sandbox benchmarks/tomato1/scripts/dashboard_main_vs_branch.py
#
# (or `marimo run` for a read-only view).

import marimo

__generated_with = "0.23.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import csv
    from collections import defaultdict
    from pathlib import Path

    import marimo as mo

    return Path, csv, defaultdict, mo


@app.cell
def _(mo):
    mo.md("""
    # var-calling — main vs branch comparison

    Side-by-side wall time + peak RSS for two `pop_var_caller` binaries
    on the same PSP cohort sweep. Source TSV is written by
    `perf_main_vs_branch.py`.

    - **wall time** — `psutil`-polled, real-clock end-to-end.
    - **peak RSS** — max over all polling ticks of the sum across the
      process tree (parent + rayon workers + any helper subprocesses).
    - **md5 match** — md5 of the VCF body with `##source=` /
      `##commandline=` stripped. If this is `no` for any N, the two
      callers are computing different answers and the timing chart
      is not interpretable.
    """)
    return


@app.cell
def _(Path):
    # benchmarks/tomato1/scripts/dashboard_main_vs_branch.py
    #   -> scripts -> tomato1 -> benchmarks -> repo_root
    repo_root = Path(__file__).resolve().parent.parent.parent.parent
    tsv_path = repo_root / "benchmarks/tomato1/results/perf/main_vs_branch.tsv"
    return repo_root, tsv_path


@app.cell
def _(csv, defaultdict, mo, tsv_path):
    if not tsv_path.exists():
        rows = []
        view = mo.md(
            f"""
            _No TSV found at `{tsv_path}`._

            Run the driver first:

            ```sh
            uv run --script benchmarks/tomato1/scripts/perf_main_vs_branch.py
            ```
            """
        )
    else:
        rows = []
        with tsv_path.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows.append({
                    "caller": r["caller"],
                    "n_samples": int(r["n_samples"]),
                    "wall_seconds": float(r["wall_seconds"]),
                    "peak_rss_mb": float(r["peak_rss_mb"]),
                    "exit_code": int(r["exit_code"]),
                    "vcf_md5": r["vcf_md5"],
                })
        view = mo.md(f"Loaded {len(rows)} rows from `{tsv_path.name}`.")

    by_caller_n = defaultdict(dict)
    for r in rows:
        by_caller_n[r["caller"]][r["n_samples"]] = r

    sizes = sorted({r["n_samples"] for r in rows})
    callers = sorted(by_caller_n.keys())
    view
    return by_caller_n, callers, rows, sizes, view


@app.cell
def _(by_caller_n, mo, rows, sizes):
    if not rows:
        identity_table = mo.md("")
    else:
        lines = ["| N | main md5 | branch md5 | match | main exit | branch exit |",
                 "|---|----------|------------|-------|-----------|-------------|"]
        for n in sizes:
            m_row = by_caller_n.get("main", {}).get(n)
            b_row = by_caller_n.get("branch", {}).get(n)
            m_md5 = (m_row or {}).get("vcf_md5", "") or "—"
            b_md5 = (b_row or {}).get("vcf_md5", "") or "—"
            m_exit = (m_row or {}).get("exit_code", "—")
            b_exit = (b_row or {}).get("exit_code", "—")
            if m_md5 != "—" and b_md5 != "—":
                match = "✓" if m_md5 == b_md5 else "✗ DIVERGES"
            else:
                match = "—"
            lines.append(
                f"| {n} | `{m_md5[:12]}` | `{b_md5[:12]}` | {match} | "
                f"{m_exit} | {b_exit} |"
            )
        identity_table = mo.md(
            "## VCF byte-identity (header-stripped)\n\n" + "\n".join(lines)
        )
    identity_table
    return


@app.cell
def _(by_caller_n, callers, mo, rows, sizes):
    if not rows:
        wall_chart = mo.md("")
    else:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 5))
        for caller in callers:
            xs = []
            ys = []
            for n in sizes:
                r = by_caller_n[caller].get(n)
                if r and r["exit_code"] == 0:
                    xs.append(n)
                    ys.append(r["wall_seconds"])
            ax.plot(xs, ys, marker="o", label=caller)
        ax.set_xlabel("N samples")
        ax.set_ylabel("wall seconds")
        ax.set_title("Wall time vs N samples")
        ax.grid(True, alpha=0.3)
        ax.legend()
        wall_chart = fig
    wall_chart
    return


@app.cell
def _(by_caller_n, callers, mo, rows, sizes):
    if not rows:
        rss_chart = mo.md("")
    else:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 5))
        for caller in callers:
            xs = []
            ys = []
            for n in sizes:
                r = by_caller_n[caller].get(n)
                if r and r["exit_code"] == 0:
                    xs.append(n)
                    ys.append(r["peak_rss_mb"])
            ax.plot(xs, ys, marker="o", label=caller)
        ax.set_xlabel("N samples")
        ax.set_ylabel("peak RSS (MB)")
        ax.set_title("Peak RSS vs N samples")
        ax.grid(True, alpha=0.3)
        ax.legend()
        rss_chart = fig
    rss_chart
    return


@app.cell
def _(by_caller_n, mo, rows, sizes):
    if not rows or "main" not in by_caller_n or "branch" not in by_caller_n:
        ratio_chart = mo.md("")
    else:
        import matplotlib.pyplot as plt

        speedups = []
        rss_ratios = []
        xs = []
        for n in sizes:
            m_row = by_caller_n["main"].get(n)
            b_row = by_caller_n["branch"].get(n)
            if not (m_row and b_row):
                continue
            if m_row["exit_code"] != 0 or b_row["exit_code"] != 0:
                continue
            if b_row["wall_seconds"] <= 0 or m_row["peak_rss_mb"] <= 0:
                continue
            xs.append(n)
            speedups.append(m_row["wall_seconds"] / b_row["wall_seconds"])
            rss_ratios.append(b_row["peak_rss_mb"] / m_row["peak_rss_mb"])

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

        ax1.plot(xs, speedups, marker="o", color="tab:green")
        ax1.axhline(1.0, color="grey", linewidth=0.8, linestyle="--")
        ax1.set_xlabel("N samples")
        ax1.set_ylabel("speedup (main_wall / branch_wall)")
        ax1.set_title("Branch wall-time speedup over main\n(>1 = branch faster)")
        ax1.grid(True, alpha=0.3)

        ax2.plot(xs, rss_ratios, marker="o", color="tab:red")
        ax2.axhline(1.0, color="grey", linewidth=0.8, linestyle="--")
        ax2.set_xlabel("N samples")
        ax2.set_ylabel("RSS ratio (branch_rss / main_rss)")
        ax2.set_title("Branch peak-RSS ratio vs main\n(>1 = branch uses more memory)")
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()
        ratio_chart = fig
    ratio_chart
    return


@app.cell
def _(mo, rows, sizes, by_caller_n):
    if not rows:
        raw_table = mo.md("")
    else:
        lines = [
            "## Raw numbers",
            "",
            "| N | main wall | branch wall | main RSS | branch RSS | speedup | rss ratio |",
            "|---|-----------|-------------|----------|------------|---------|-----------|",
        ]
        for n in sizes:
            m_row = by_caller_n.get("main", {}).get(n)
            b_row = by_caller_n.get("branch", {}).get(n)
            if not (m_row and b_row):
                lines.append(f"| {n} | — | — | — | — | — | — |")
                continue
            if m_row["exit_code"] != 0 or b_row["exit_code"] != 0:
                lines.append(
                    f"| {n} | exit={m_row['exit_code']} | exit={b_row['exit_code']} "
                    f"| — | — | — | — |"
                )
                continue
            speedup = m_row["wall_seconds"] / b_row["wall_seconds"] if b_row["wall_seconds"] else 0.0
            rss_ratio = b_row["peak_rss_mb"] / m_row["peak_rss_mb"] if m_row["peak_rss_mb"] else 0.0
            lines.append(
                f"| {n} | {m_row['wall_seconds']:.1f}s | {b_row['wall_seconds']:.1f}s "
                f"| {m_row['peak_rss_mb']:.0f} MB | {b_row['peak_rss_mb']:.0f} MB "
                f"| {speedup:.2f}× | {rss_ratio:.2f}× |"
            )
        raw_table = mo.md("\n".join(lines))
    raw_table
    return


if __name__ == "__main__":
    app.run()
