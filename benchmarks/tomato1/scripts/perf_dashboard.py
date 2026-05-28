# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "matplotlib",
# ]
# ///

# Marimo performance dashboard — runtime + peak RSS vs N samples
# across the four caller modes (ours-direct, ours-psp, freebayes,
# gatk-joint). Reads per-caller TSVs written by the perf_<caller>.py
# scripts; skips any TSV that doesn't exist yet so partial runs still
# render.
#
# Recommended invocation:
#
#   uvx marimo edit --sandbox benchmarks/tomato1/scripts/perf_dashboard.py
#
# Companion to (not replacement of) dashboard.py — that one is about
# correctness (variant agreement, QUAL distributions); this one is
# about performance scaling.

import marimo

__generated_with = "0.23.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import csv
    from pathlib import Path

    import marimo as mo

    return Path, csv, mo


@app.cell
def _(mo):
    mo.md("""
    # Caller scaling — runtime and peak RAM vs N samples

    Two measurement pairs:

    - **A. Full cohort from BAMs** (caller-vs-caller, whole process):
      `ours_whole_pipeline` (parallel pileup + joint var-calling) vs
      `freebayes` (one-shot CRAMs → VCF). End-to-end time and peak
      RAM starting from CRAM inputs.
    - **B. Per-sample intermediate** (CRAM → per-sample file):
      `ours_pileup` (CRAM → PSP, 1 MiB blocks),
      `ours_pileup_smallblock` (CRAM → PSP, 256 KiB blocks),
      and `gatk_haplotype_caller` (CRAM → GVCF). All run THREADS
      instances concurrently, each single-threaded. Block-size effect
      on this stage is expected to be small (block target is a writer
      flush threshold, not a working-set knob).
    - **C. Joint genotyping stage** (per-sample intermediate → cohort VCF):
      `ours_joint` (PSPs at 1 MiB blocks → VCF),
      `ours_joint_smallblock` (PSPs at 256 KiB blocks → VCF),
      and `gatk_joint` (CombineGVCFs + GenotypeGVCFs). Intermediates
      assumed pre-built; per-sample time is excluded. This is where
      the block-size effect should bite — the cohort reader holds one
      decoded block per `n_threads × N` worker, so 4× smaller blocks
      should roughly 4× lower peak heap.
    - **D. Synthetic far-N scaling** (PSP → VCF, post-fcef495 512 KiB
      block default): `ours_joint` (real cohort, N=1..26) +
      `scaling_synthetic` (one tomato PSP replicated into a synthetic
      cohort, N=50/200/1000). Same shape as Pair C, extended into the
      range where memory linearity and CPU distribution are
      load-bearing for future architecture decisions.

    Inputs: per-caller TSVs at `results/perf/<caller>.tsv`,
    written by `perf_<caller>.py`.
    """)
    return


@app.cell
def _(Path):
    # tmp/dashboard.py -> tmp/  /  scripts/dashboard.py -> scripts/
    test_dir = Path(__file__).resolve().parent.parent
    perf_dir = test_dir / "results" / "perf"
    callers = (
        "ours_whole_pipeline",
        "ours_pileup", "ours_pileup_smallblock",
        "ours_joint",  "ours_joint_smallblock",
        "freebayes", "gatk_haplotype_caller", "gatk_joint",
        "scaling_synthetic",
    )
    tsv_paths = {c: perf_dir / f"{c}.tsv" for c in callers}
    return callers, perf_dir, test_dir, tsv_paths


@app.cell
def _(csv, mo, tsv_paths):
    # Load each TSV; skip missing ones with a friendly summary so a
    # partial run (e.g. you've only run perf_ours_whole_pipeline.py so far)
    # still renders something useful.
    rows_by_caller: dict[str, list[dict]] = {}
    present, missing = [], []
    for _caller, _path in tsv_paths.items():
        if not _path.exists():
            missing.append(_caller)
            continue
        with _path.open() as _fh:
            _reader = csv.DictReader(_fh, delimiter="\t")
            _rows = []
            for _r in _reader:
                _rows.append({
                    "n_samples": int(_r["n_samples"]),
                    "wall_seconds": float(_r["wall_seconds"]),
                    "peak_rss_mb": float(_r["peak_rss_mb"]),
                    "exit_code": int(_r["exit_code"]),
                })
        rows_by_caller[_caller] = _rows
        present.append(_caller)

    lines = [f"- **{c}**: {len(rows_by_caller[c])} rows" for c in present]
    if missing:
        lines.append("")
        lines.append("_Missing (not run yet):_ " + ", ".join(f"`{m}`" for m in missing))
    status_view = mo.md("### TSV inputs\n\n" + "\n".join(lines))
    status_view
    return missing, present, rows_by_caller, status_view


@app.cell
def _(mo):
    # Axis-scale toggles. X is sample count (small integer); log only
    # really helps when the sweep spans a wide range (1..18 is borderline,
    # but log lets us see whether scaling is linear vs super-linear).
    # Y log is the bigger lever for both runtime and RSS when the two
    # callers differ by an order of magnitude.
    xscale = mo.ui.radio(
        options=["linear", "log"], value="linear", label="x-axis scale (N samples)",
    )
    yscale = mo.ui.radio(
        options=["linear", "log"], value="linear", label="y-axis scale (runtime / RSS)",
    )
    mo.hstack([xscale, yscale])
    return xscale, yscale


@app.cell
def _():
    import matplotlib.pyplot as plt

    # Stable colour per caller so the same caller looks the same
    # across the two pair figures.
    PALETTE = {
        "ours_whole_pipeline":     "#1f77b4",  # blue
        "ours_pileup":             "#17becf",  # cyan
        "ours_pileup_smallblock":  "#08527a",  # darker cyan/teal
        "ours_joint":              "#2ca02c",  # green
        "ours_joint_smallblock":   "#1a5d1a",  # darker green
        "freebayes":               "#d62728",  # red
        "gatk_haplotype_caller":   "#e377c2",  # pink
        "gatk_joint":              "#9467bd",  # purple
        "scaling_synthetic":       "#ff7f0e",  # orange
    }

    def plot_pair(rows_by_caller, pair, *, title, xscale_val, yscale_val):
        """Render runtime + RSS panels for the callers in `pair`
        (2-or-more — supports default + smallblock variants alongside
        the cross-tool comparison). Plots whichever callers are
        present; reports the absent ones in a placeholder if NONE are
        present.

        The two panels share the x-axis (N samples) but have
        independent y-axes — runtime and RSS are on incomparable
        scales, and yoking them via sharey=True would squash one."""
        present = [c for c in pair if c in rows_by_caller]
        if not present:
            return None, list(pair)
        fig, (ax_t, ax_m) = plt.subplots(1, 2, figsize=(14, 5), sharex=True)
        for caller in present:
            ok_rows = [r for r in rows_by_caller[caller] if r["exit_code"] == 0]
            if not ok_rows:
                continue
            xs = [r["n_samples"] for r in ok_rows]
            ts = [r["wall_seconds"] for r in ok_rows]
            ms = [r["peak_rss_mb"] for r in ok_rows]
            colour = PALETTE.get(caller)
            ax_t.plot(xs, ts, marker="o", label=caller, color=colour)
            ax_m.plot(xs, ms, marker="o", label=caller, color=colour)
        for ax, ylabel, panel_title in (
            (ax_t, "wall-clock seconds", "Runtime"),
            (ax_m, "peak RSS (MB)",      "Peak memory"),
        ):
            ax.set_xlabel("N samples")
            ax.set_ylabel(ylabel)
            ax.set_title(panel_title)
            ax.set_xscale(xscale_val)
            ax.set_yscale(yscale_val)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=9)
        fig.suptitle(title, fontsize=13, y=1.02)
        fig.tight_layout()
        return fig, []

    return PALETTE, plot_pair, plt


@app.cell
def _(mo, plot_pair, rows_by_caller, xscale, yscale):
    # Pair A — full cohort from BAMs (whole-process, caller vs caller).
    fig_a, missing_a = plot_pair(
        rows_by_caller,
        ("ours_whole_pipeline", "freebayes"),
        title="A. Full cohort from BAMs — ours_whole_pipeline vs freebayes",
        xscale_val=xscale.value,
        yscale_val=yscale.value,
    )
    pair_a_view = fig_a if fig_a is not None else mo.md(
        "_(Pair A unavailable; missing TSV(s): "
        + ", ".join(f"`{c}`" for c in missing_a) + ".)_"
    )
    pair_a_view
    return fig_a, missing_a, pair_a_view


@app.cell
def _(mo, plot_pair, rows_by_caller, xscale, yscale):
    # Pair B — per-sample intermediate (CRAM → per-sample file), both
    # measured with THREADS-way concurrency for an apples-to-apples
    # parallelism budget.
    fig_b, missing_b = plot_pair(
        rows_by_caller,
        ("ours_pileup", "ours_pileup_smallblock", "gatk_haplotype_caller"),
        title="B. Per-sample intermediate — ours_pileup (default + smallblock) vs gatk_haplotype_caller",
        xscale_val=xscale.value,
        yscale_val=yscale.value,
    )
    pair_b_view = fig_b if fig_b is not None else mo.md(
        "_(Pair B unavailable; missing TSV(s): "
        + ", ".join(f"`{c}`" for c in missing_b) + ".)_"
    )
    pair_b_view
    return fig_b, missing_b, pair_b_view


@app.cell
def _(mo, plot_pair, rows_by_caller, xscale, yscale):
    # Pair C — joint genotyping stage only (per-sample intermediate
    # → cohort VCF). Per-sample build time excluded.
    fig_c, missing_c = plot_pair(
        rows_by_caller,
        ("ours_joint", "ours_joint_smallblock", "gatk_joint"),
        title="C. Joint genotyping stage only — ours_joint (default + smallblock) vs gatk_joint",
        xscale_val=xscale.value,
        yscale_val=yscale.value,
    )
    pair_c_view = fig_c if fig_c is not None else mo.md(
        "_(Pair C unavailable; missing TSV(s): "
        + ", ".join(f"`{c}`" for c in missing_c) + ".)_"
    )
    pair_c_view
    return fig_c, missing_c, pair_c_view


@app.cell
def _(mo, plot_pair, rows_by_caller, xscale, yscale):
    # Pair D — far-N scaling. ours_joint covers N=1..26 from the real
    # tomato1 cohort; scaling_synthetic extends the same shape to
    # N=50/200/1000 by replicating one tomato PSP into a synthetic
    # cohort. Both run the joint PSP → VCF stage end-to-end with the
    # post-fcef495 512 KiB block default, so the two curves should
    # extrapolate cleanly. Log-x usually wins here because of the
    # 30× span; flip the toggle above if you prefer linear.
    fig_d, missing_d = plot_pair(
        rows_by_caller,
        ("ours_joint", "scaling_synthetic"),
        title="D. Joint stage, far-N scaling — ours_joint (N=1..26) + scaling_synthetic (N=50..1000)",
        xscale_val=xscale.value,
        yscale_val=yscale.value,
    )
    pair_d_view = fig_d if fig_d is not None else mo.md(
        "_(Pair D unavailable; missing TSV(s): "
        + ", ".join(f"`{c}`" for c in missing_d) + ".)_"
    )
    pair_d_view
    return fig_d, missing_d, pair_d_view


@app.cell
def _(mo, rows_by_caller):
    # Raw numbers table — useful when you want exact values instead of
    # eyeballing the plot. One block per caller; failures (exit_code
    # != 0) get a ⚠ marker.
    if not rows_by_caller:
        table_view = mo.md("")
    else:
        blocks = []
        for _caller, _rows in rows_by_caller.items():
            _md = f"#### {_caller}\n\n"
            _md += "| N | wall (s) | peak RSS (MB) | exit |\n|---:|---:|---:|---:|\n"
            for _r in sorted(_rows, key=lambda x: x["n_samples"]):
                _marker = "" if _r["exit_code"] == 0 else " ⚠"
                _md += (
                    f"| {_r['n_samples']} | {_r['wall_seconds']:.1f} | "
                    f"{_r['peak_rss_mb']:.0f} | {_r['exit_code']}{_marker} |\n"
                )
            blocks.append(_md)
        table_view = mo.md("### Raw measurements\n\n" + "\n".join(blocks))
    table_view
    return (table_view,)


if __name__ == "__main__":
    app.run()
