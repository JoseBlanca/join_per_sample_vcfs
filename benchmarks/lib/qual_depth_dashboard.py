# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "matplotlib",
# ]
# ///

# Marimo QUAL-vs-depth dashboard — how each caller's variant-QUAL
# distribution, and the separation between true and false positives,
# moves with sequencing depth.
#
# The native HG002 benchmark CRAM is ~301x, far above realistic WGS, so
# depth is created by subsampling (benchmarks/lib/run_depth_sweep.sh):
# the CRAM is downsampled to a ladder of global coverage depths, all
# three callers are re-run at each, and every call is scored TP/FP
# against the GIAB truth set. This reads the merged tables that sweep
# writes:
#   results/depth_sweep/qual_dist_by_depth.tsv   (depth,caller,class,status,qual,dp)
#   results/depth_sweep/accuracy_by_depth.tsv    (depth,caller,class,tp,fp,fn,prec,rec,f1)
#
# Companion to comparison_dashboard.py (single-depth TP/FP QUAL + ROC).
# That one answers "how separable are TP and FP at native depth"; this
# one answers "how does that separation change as depth drops".
#
# Recommended invocation:
#
#   uvx marimo edit --sandbox benchmarks/lib/qual_depth_dashboard.py
#
# (or `marimo run` for a read-only view).

import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell
def _():
    import bisect
    import csv
    from pathlib import Path

    import marimo as mo

    return Path, bisect, csv, mo


@app.cell
def _(mo):
    mo.md("""
    # QUAL vs sequencing depth — TP/FP separation per caller

    The native HG002 CRAM (~301x) is downsampled to a ladder of
    global depths; every caller is re-run at each depth and its
    calls scored **TP** (in the GIAB truth set) or **FP** (not).

    - **TP/FP QUAL grid** — the raw distributions, one row per
      depth, one column per caller. A caller whose FP mass sits at
      low QUAL while TP mass sits high is *gateable*: a QUAL cutoff
      cleans it up. Overlap that persists as depth drops is the
      interesting failure mode.
    - **Separation vs depth** — median TP and FP QUAL (with IQR
      bands) and a separation score, traced across the depth ladder.
    - **Accuracy vs depth** — precision / recall / F1 from the same
      runs.
    - **Best QUAL cutoff vs depth** — the F1-optimal threshold and
      the F1 it achieves, per caller, as depth changes.

    Source: `results/depth_sweep/{qual_dist,accuracy}_by_depth.tsv`.
    """)
    return


@app.cell
def _(Path):
    # lib -> benchmarks; find every benchmark with a depth sweep.
    benchmarks_dir = Path(__file__).resolve().parent.parent
    sweeps = {
        p.parent.parent.name: p.parent          # bench name -> .../results/depth_sweep
        for p in sorted(benchmarks_dir.glob("*/results/depth_sweep/qual_dist_by_depth.tsv"))
    }
    return (sweeps,)


@app.cell
def _(mo, sweeps):
    if not sweeps:
        bench_selector = None
        view = mo.md(
            "_No `qual_dist_by_depth.tsv` found under "
            "`benchmarks/*/results/depth_sweep/`. Run the sweep first:_\n\n"
            "    bash benchmarks/lib/run_depth_sweep.sh <config> host\n"
            "    DEV_EXTRA_MOUNT=$HOME/genomes ./scripts/dev.sh \\\n"
            "        bash benchmarks/lib/run_depth_sweep.sh <config> gatk\n"
            "    bash benchmarks/lib/run_depth_sweep.sh <config> compare\n"
            "    bash benchmarks/lib/run_depth_sweep.sh <config> merge\n"
        )
    else:
        bench_selector = mo.ui.dropdown(
            options=list(sweeps.keys()),
            value=next(iter(sweeps)),
            label="benchmark",
        )
        view = bench_selector
    view
    return (bench_selector,)


@app.cell
def _(bench_selector, csv, mo, sweeps):
    mo.stop(bench_selector is None, mo.md(""))
    sweep_dir = sweeps[bench_selector.value]

    # qual rows: (depth:int, caller, class, status, qual:float, dp:int|None)
    qual_rows = []
    with (sweep_dir / "qual_dist_by_depth.tsv").open() as _fh:
        for _r in csv.DictReader(_fh, delimiter="\t"):
            try:
                _q = float(_r["qual"])
            except ValueError:
                continue
            _dp = None
            try:
                _dp = int(_r["dp"])
            except (ValueError, KeyError):
                pass
            qual_rows.append((int(_r["depth"]), _r["caller"], _r["class"],
                              _r["status"], _q, _dp))

    acc_rows = []
    _accp = sweep_dir / "accuracy_by_depth.tsv"
    if _accp.exists():
        with _accp.open() as _fh:
            for _r in csv.DictReader(_fh, delimiter="\t"):
                acc_rows.append({
                    "depth": int(_r["depth"]), "caller": _r["caller"],
                    "class": _r["class"], "tp": int(_r["tp"]),
                    "fp": int(_r["fp"]), "fn": int(_r["fn"]),
                    "precision": float(_r["precision"]),
                    "recall": float(_r["recall"]), "f1": float(_r["f1"]),
                })

    depths = sorted({r[0] for r in qual_rows})
    callers = sorted({r[1] for r in qual_rows})
    classes = sorted({r[2] for r in qual_rows}, reverse=True)  # snps before indels
    mo.md(
        f"**{bench_selector.value}** — {len(qual_rows)} call rows · "
        f"depths: {', '.join(str(d) + 'x' for d in depths)} · "
        f"callers: {', '.join('`' + c + '`' for c in callers)}"
    )
    return acc_rows, callers, classes, depths, qual_rows


@app.cell
def _(classes, mo):
    class_sel = mo.ui.dropdown(
        options=classes, value=classes[0], label="variant class"
    )
    qual_max = mo.ui.slider(
        start=100, stop=15000, step=100, value=3000, label="QUAL clip (x-max)"
    )
    n_bins = mo.ui.slider(start=20, stop=120, step=10, value=60, label="bins")
    density = mo.ui.switch(value=True, label="density (normalise TP/FP areas)")
    logy = mo.ui.switch(value=False, label="log y")
    controls = mo.hstack([class_sel, qual_max, n_bins, density, logy],
                         justify="start", gap=1.5)
    controls
    return class_sel, density, logy, n_bins, qual_max


@app.cell
def _():
    import matplotlib.pyplot as plt

    PALETTE = {
        "ours": "#1f77b4",       # blue
        "gatk": "#9467bd",       # purple
        "freebayes": "#d62728",  # red
    }
    TP_COLOR, FP_COLOR = "#2ca02c", "#d62728"
    return FP_COLOR, PALETTE, TP_COLOR, plt


@app.cell
def _(FP_COLOR, TP_COLOR, plt):
    def plot_qual_grid(qual_rows, callers, depths, cls,
                       qual_max, n_bins, density, logy):
        """Grid of TP-vs-FP QUAL histograms: rows = depth (ascending),
        cols = caller. Shared x (0..qual_max, values clipped into the
        last bin) and bin count so every panel lines up. Density mode
        normalises each histogram to area 1 so TP/FP *shapes* compare
        despite wildly different counts (e.g. freebayes' FP flood)."""
        nrow, ncol = len(depths), len(callers)
        fig, axes = plt.subplots(
            nrow, ncol, figsize=(4.4 * ncol, 2.9 * nrow),
            sharex=True, squeeze=False,
        )
        bin_kw = dict(bins=n_bins, range=(0, qual_max), density=density, alpha=0.55)
        for ri, d in enumerate(depths):
            for ci, c in enumerate(callers):
                ax = axes[ri][ci]
                tp = [min(r[4], qual_max) for r in qual_rows
                      if r[0] == d and r[1] == c and r[2] == cls and r[3] == "TP"]
                fp = [min(r[4], qual_max) for r in qual_rows
                      if r[0] == d and r[1] == c and r[2] == cls and r[3] == "FP"]
                if tp:
                    ax.hist(tp, color=TP_COLOR, label=f"TP {len(tp)}", **bin_kw)
                if fp:
                    ax.hist(fp, color=FP_COLOR, label=f"FP {len(fp)}", **bin_kw)
                if logy:
                    ax.set_yscale("log")
                ax.grid(True, alpha=0.3)
                ax.legend(fontsize=7, loc="upper right")
                if ri == 0:
                    ax.set_title(c, fontsize=11, fontweight="bold")
                if ci == 0:
                    ax.set_ylabel(f"{d}x\n" + ("density" if density else "count"),
                                  fontsize=9)
                if ri == nrow - 1:
                    ax.set_xlabel("QUAL")
        fig.suptitle(
            f"{cls} QUAL — TP vs FP, by depth (rows) × caller (cols); "
            f"x clipped at {qual_max:g}",
            fontsize=13, y=1.005,
        )
        fig.tight_layout()
        return fig

    return (plot_qual_grid,)


@app.cell
def _(
    callers,
    class_sel,
    density,
    depths,
    logy,
    mo,
    n_bins,
    plot_qual_grid,
    qual_max,
    qual_rows,
):
    mo.stop(not qual_rows, mo.md(""))
    plot_qual_grid(qual_rows, callers, depths, class_sel.value,
                   qual_max.value, n_bins.value, density.value, logy.value)
    return


@app.cell
def _(mo):
    mo.md("""
    ## QUAL box plots — tools side by side, per coverage

    For each coverage depth, two panels: **FP** (left) and **TP**
    (right). Within each panel one box per caller, so you can read
    off, at that depth, which tool's FP QUAL sits lowest (most
    gateable) and how the TP QUAL distributions line up. Boxes are
    median + IQR; whiskers span the 5th–95th percentile (outliers
    hidden). FP and TP panels auto-scale independently — FP QUAL is
    far lower than TP QUAL, so a shared scale would flatten the FP
    boxes to a line.

    The controls below govern **both** box-plot sections. Tool/depth
    QUAL scales differ wildly, so turn **auto y** off and drag the
    range slider to zoom — e.g. a tight `0–600` range to read the FP
    boxes, or the full span for the high-depth TP boxes. (The `log y`
    switch in the top control bar also applies here.)
    """)
    return


@app.cell
def _(mo):
    box_yauto = mo.ui.switch(value=True, label="auto y-axis (per panel)")
    box_ylim = mo.ui.range_slider(
        start=0, stop=20000, step=100, value=[0, 3000], show_value=True,
        label="box y-axis range (used when auto is off)",
    )
    box_controls = mo.hstack([box_yauto, box_ylim], justify="start", gap=1.5)
    box_controls
    return box_yauto, box_ylim


@app.function
def draw_boxes(ax, groups, logy, ylim=None):
    """groups: list of (label, colour, values). Skips empty groups
    (boxplot rejects empty data) and keeps x positions contiguous.
    ylim=None auto-scales the panel; otherwise (lo, hi) fixes the y-axis
    (log mode clamps a 0 lower bound up to a small positive)."""
    data, labels, colours = [], [], []
    for label, colour, vals in groups:
        if vals:
            data.append(vals); labels.append(label); colours.append(colour)
    if not data:
        ax.text(0.5, 0.5, "no calls", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="#999")
        ax.set_xticks([])
        return
    bp = ax.boxplot(data, patch_artist=True, showfliers=False,
                    whis=(5, 95), widths=0.6,
                    medianprops=dict(color="black", lw=1.4))
    for patch, colour in zip(bp["boxes"], colours):
        patch.set_facecolor(colour); patch.set_alpha(0.55)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=9)
    if logy:
        ax.set_yscale("log")
    if ylim is not None:
        lo, hi = ylim
        if logy and lo <= 0:
            lo = 1.0
        ax.set_ylim(lo, hi)
    ax.grid(True, axis="y", alpha=0.3)


@app.cell
def _(PALETTE, plt):
    def plot_box_by_coverage(qual_rows, callers, depths, cls, logy, ylim):
        """Rows = depth; cols = [FP, TP]; one box per caller in each
        panel. Reads down a column to compare callers at a fixed depth."""
        nrow = len(depths)
        fig, axes = plt.subplots(nrow, 2, figsize=(11, 2.8 * nrow),
                                 squeeze=False)
        for ri, d in enumerate(depths):
            for ci, status in enumerate(("FP", "TP")):
                ax = axes[ri][ci]
                groups = [
                    (c, PALETTE.get(c, "#555"),
                     [r[4] for r in qual_rows
                      if r[0] == d and r[1] == c and r[2] == cls
                      and r[3] == status])
                    for c in callers
                ]
                draw_boxes(ax, groups, logy, ylim)
                if ri == 0:
                    ax.set_title(status, fontsize=12, fontweight="bold")
                if ci == 0:
                    ax.set_ylabel(f"{d}x\nQUAL", fontsize=9)
        fig.suptitle(f"{cls} QUAL by caller, per coverage "
                     f"(left FP · right TP)", fontsize=13, y=1.003)
        fig.tight_layout()
        return fig

    return (plot_box_by_coverage,)


@app.cell
def _(box_ylim, box_yauto, callers, class_sel, depths, logy, mo,
      plot_box_by_coverage, qual_rows):
    mo.stop(not qual_rows, mo.md(""))
    _ylim = None if box_yauto.value else tuple(box_ylim.value)
    plot_box_by_coverage(qual_rows, callers, depths, class_sel.value,
                         logy.value, _ylim)
    return


@app.cell
def _(mo):
    mo.md("""
    ## QUAL box plots — coverage trend, per tool

    Flipped axis: one row per **caller**, two panels (**FP** left,
    **TP** right), one box per coverage depth. Reading left→right
    within a panel shows how that caller's QUAL distribution shifts
    as depth climbs (the QUAL scale grows roughly linearly with
    depth, so a fixed cutoff won't transfer). Panels auto-scale unless
    you fix the range with the controls above (shared with the
    per-coverage view).
    """)
    return


@app.cell
def _(PALETTE, plt):
    def plot_box_by_tool(qual_rows, callers, depths, cls, logy, ylim):
        """Rows = caller; cols = [FP, TP]; one box per depth in each
        panel. Reads across a row to see the coverage trend per tool."""
        nrow = len(callers)
        fig, axes = plt.subplots(nrow, 2, figsize=(12, 3.0 * nrow),
                                 squeeze=False)
        for ri, c in enumerate(callers):
            colour = PALETTE.get(c, "#555")
            for ci, status in enumerate(("FP", "TP")):
                ax = axes[ri][ci]
                groups = [
                    (f"{d}x", colour,
                     [r[4] for r in qual_rows
                      if r[0] == d and r[1] == c and r[2] == cls
                      and r[3] == status])
                    for d in depths
                ]
                draw_boxes(ax, groups, logy, ylim)
                ax.set_xlabel("depth")
                if ri == 0:
                    ax.set_title(status, fontsize=12, fontweight="bold")
                if ci == 0:
                    ax.set_ylabel(f"{c}\nQUAL", fontsize=10, fontweight="bold")
        fig.suptitle(f"{cls} QUAL vs coverage, per caller "
                     f"(left FP · right TP)", fontsize=13, y=1.003)
        fig.tight_layout()
        return fig

    return (plot_box_by_tool,)


@app.cell
def _(box_ylim, box_yauto, callers, class_sel, depths, logy, mo,
      plot_box_by_tool, qual_rows):
    mo.stop(not qual_rows, mo.md(""))
    _ylim = None if box_yauto.value else tuple(box_ylim.value)
    plot_box_by_tool(qual_rows, callers, depths, class_sel.value,
                     logy.value, _ylim)
    return


@app.cell
def _(mo):
    mo.md("""
    ## Separation vs depth

    Per caller: median TP QUAL and median FP QUAL across the depth
    ladder, shaded with their inter-quartile ranges. The wider and
    cleaner the gap between the green (TP) and red (FP) bands, the
    more a single QUAL cutoff can separate good from bad calls at
    that depth. The lower panel is a **separation score** — the
    fraction of FP calls whose QUAL falls below the *median* TP QUAL
    (1.0 = every FP is below the typical TP, i.e. perfectly
    gateable; 0.5 = FP and TP medians coincide).
    """)
    return


@app.cell
def _(PALETTE, plt):
    def _pct(xs, q):
        if not xs:
            return float("nan")
        xs = sorted(xs)
        i = (len(xs) - 1) * q
        lo, hi = int(i), min(int(i) + 1, len(xs) - 1)
        return xs[lo] + (xs[hi] - xs[lo]) * (i - lo)

    def plot_separation(qual_rows, callers, depths, cls):
        """Two stacked panels. Top: median TP & FP QUAL vs depth per
        caller, IQR shaded. Bottom: separation score (share of FP below
        the median TP QUAL) vs depth."""
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(9, 8.5), sharex=True,
            gridspec_kw={"height_ratios": [2, 1]},
        )
        for c in callers:
            colour = PALETTE.get(c, "#555555")
            tp_med, tp_lo, tp_hi = [], [], []
            fp_med, fp_lo, fp_hi = [], [], []
            sep = []
            for d in depths:
                tp = [r[4] for r in qual_rows
                      if r[0] == d and r[1] == c and r[2] == cls and r[3] == "TP"]
                fp = [r[4] for r in qual_rows
                      if r[0] == d and r[1] == c and r[2] == cls and r[3] == "FP"]
                tp_med.append(_pct(tp, 0.5)); tp_lo.append(_pct(tp, 0.25)); tp_hi.append(_pct(tp, 0.75))
                fp_med.append(_pct(fp, 0.5)); fp_lo.append(_pct(fp, 0.25)); fp_hi.append(_pct(fp, 0.75))
                if tp and fp:
                    cut = _pct(tp, 0.5)
                    sep.append(sum(1 for v in fp if v < cut) / len(fp))
                else:
                    sep.append(float("nan"))
            ax1.plot(depths, tp_med, "-o", color=colour, lw=1.8, ms=4,
                     label=f"{c} TP")
            ax1.fill_between(depths, tp_lo, tp_hi, color=colour, alpha=0.12)
            ax1.plot(depths, fp_med, "--s", color=colour, lw=1.4, ms=4,
                     alpha=0.8, label=f"{c} FP")
            ax1.fill_between(depths, fp_lo, fp_hi, color=colour, alpha=0.06)
            ax2.plot(depths, sep, "-o", color=colour, lw=1.8, ms=4, label=c)
        ax1.set_ylabel("QUAL (median, IQR band)")
        ax1.set_title(f"{cls} — TP (solid) vs FP (dashed) QUAL by depth")
        ax1.grid(alpha=0.3)
        ax1.legend(fontsize=8, ncol=len(callers))
        ax2.axhline(0.5, color="#999", ls=":", lw=1)
        ax2.set_ylabel("FP share below median TP QUAL")
        ax2.set_xlabel("sequencing depth (x)")
        ax2.set_ylim(0, 1.02)
        ax2.set_title(f"{cls} — gateability (1.0 = every FP below typical TP)")
        ax2.grid(alpha=0.3)
        ax2.legend(fontsize=8, ncol=len(callers))
        fig.tight_layout()
        return fig

    return (plot_separation,)


@app.cell
def _(callers, class_sel, depths, mo, plot_separation, qual_rows):
    mo.stop(not qual_rows, mo.md(""))
    plot_separation(qual_rows, callers, depths, class_sel.value)
    return


@app.cell
def _(mo):
    mo.md("""
    ## Accuracy vs depth

    Precision / recall / F1 from the same per-depth runs (allele
    concordance vs the truth set, FILTER kept as each caller emits
    it). These are at each caller's emitted gate (=0 here), so they
    reflect the *unfiltered* call set — the QUAL sweep below is what
    a downstream filter would buy.
    """)
    return


@app.cell
def _(PALETTE, plt):
    def plot_accuracy(acc_rows, callers, depths, cls):
        metrics = ("precision", "recall", "f1")
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharex=True)
        for ax, metric in zip(axes, metrics):
            for c in callers:
                ys = []
                for d in depths:
                    rec = next((r for r in acc_rows
                                if r["depth"] == d and r["caller"] == c
                                and r["class"] == cls), None)
                    ys.append(rec[metric] if rec else float("nan"))
                ax.plot(depths, ys, "-o", color=PALETTE.get(c, "#555"),
                        lw=1.8, ms=4, label=c)
            ax.set_title(f"{cls} — {metric}")
            ax.set_xlabel("depth (x)")
            ax.set_ylabel(metric)
            ax.set_ylim(0, 1.04)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        fig.tight_layout()
        return fig

    return (plot_accuracy,)


@app.cell
def _(acc_rows, callers, class_sel, depths, mo, plot_accuracy):
    mo.stop(not acc_rows, mo.md("_No accuracy_by_depth.tsv._"))
    plot_accuracy(acc_rows, callers, depths, class_sel.value)
    return


@app.cell
def _(mo):
    mo.md("""
    ## Best QUAL cutoff vs depth

    For each (caller, depth) we sweep a QUAL threshold over the
    call set and find the F1-optimal cutoff (recall anchored to the
    full truth set, TP+FN). Left: that **optimal threshold** vs
    depth — how high you must gate to get the best trade-off. Right:
    the **F1 achieved** at it. A threshold that climbs steeply as
    depth rises means the QUAL scale is depth-dependent (a fixed
    cutoff won't generalise across depths).
    """)
    return


@app.cell
def _(PALETTE, bisect, plt):
    def _sweep_best(tpq, fpq, n_positives):
        tps, fps = sorted(tpq), sorted(fpq)
        best = None
        for t in sorted(set(tps) | set(fps)):
            tp = len(tps) - bisect.bisect_left(tps, t)
            fp = len(fps) - bisect.bisect_left(fps, t)
            prec = tp / (tp + fp) if (tp + fp) > 0 else 1.0
            rec = tp / n_positives if n_positives > 0 else 0.0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
            if best is None or f1 > best["f1"]:
                best = {"t": t, "f1": f1, "precision": prec, "recall": rec}
        return best

    def plot_best_cutoff(qual_rows, acc_rows, callers, depths, cls):
        positives = {(r["depth"], r["caller"]): r["tp"] + r["fn"]
                     for r in acc_rows if r["class"] == cls}
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))
        for c in callers:
            ts, f1s = [], []
            for d in depths:
                tpq = [r[4] for r in qual_rows
                       if r[0] == d and r[1] == c and r[2] == cls and r[3] == "TP"]
                fpq = [r[4] for r in qual_rows
                       if r[0] == d and r[1] == c and r[2] == cls and r[3] == "FP"]
                npos = positives.get((d, c), len(tpq))
                b = _sweep_best(tpq, fpq, npos) if (tpq or fpq) else None
                ts.append(b["t"] if b else float("nan"))
                f1s.append(b["f1"] if b else float("nan"))
            colour = PALETTE.get(c, "#555555")
            ax1.plot(depths, ts, "-o", color=colour, lw=1.8, ms=4, label=c)
            ax2.plot(depths, f1s, "-o", color=colour, lw=1.8, ms=4, label=c)
        ax1.set_ylabel("F1-optimal QUAL threshold")
        ax1.set_xlabel("depth (x)")
        ax1.set_title(f"{cls} — best QUAL cutoff vs depth")
        ax1.grid(alpha=0.3)
        ax1.legend(fontsize=8)
        ax2.set_ylabel("F1 at optimal cutoff")
        ax2.set_xlabel("depth (x)")
        ax2.set_ylim(0, 1.04)
        ax2.set_title(f"{cls} — best achievable F1 vs depth")
        ax2.grid(alpha=0.3)
        ax2.legend(fontsize=8)
        fig.tight_layout()
        return fig

    return (plot_best_cutoff,)


@app.cell
def _(acc_rows, callers, class_sel, depths, mo, plot_best_cutoff, qual_rows):
    mo.stop(not qual_rows, mo.md(""))
    plot_best_cutoff(qual_rows, acc_rows, callers, depths, class_sel.value)
    return


if __name__ == "__main__":
    app.run()
