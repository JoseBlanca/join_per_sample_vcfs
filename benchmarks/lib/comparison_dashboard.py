# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "matplotlib",
# ]
# ///

# Marimo accuracy dashboard — precision / recall / F1 per caller vs a
# reference callset, split by variant class (SNP, indel). Reads the
# `accuracy.tsv` written by benchmarks/lib/compare_to_truth.sh, scanning
# every benchmark under benchmarks/ that has one and offering a selector.
# The companion to perf_dashboard.py (performance scaling) — that one is
# about speed/RAM; this one is about correctness vs a reference.
#
# Recommended invocation:
#
#   uvx marimo edit --sandbox benchmarks/lib/comparison_dashboard.py
#
# (or `marimo run` for a read-only view).

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
    # Caller accuracy — precision / recall / F1 vs reference callset

    Each caller's VCF is compared to the benchmark's reference callset
    within the benchmark BED (`compare_to_truth.sh`): variants are
    left-aligned + biallelic-split, then intersected on POS+REF+ALT.
    So this scores **allele** concordance, not genotype concordance.

    - **precision** = TP / (TP + FP) — of what the caller emitted, how
      much is in the reference.
    - **recall** = TP / (TP + FN) — of the reference, how much the
      caller recovered.
    - **F1** = harmonic mean of the two.

    Source: `results/comparison/accuracy.tsv` per benchmark.
    """)
    return


@app.cell
def _(Path):
    # benchmarks/lib/comparison_dashboard.py -> lib -> benchmarks
    benchmarks_dir = Path(__file__).resolve().parent.parent
    tsv_by_bench = {
        p.parent.parent.parent.name: p
        for p in sorted(benchmarks_dir.glob("*/results/comparison/accuracy.tsv"))
    }
    return benchmarks_dir, tsv_by_bench


@app.cell
def _(mo, tsv_by_bench):
    if not tsv_by_bench:
        bench_selector = None
        view = mo.md(
            "_No `accuracy.tsv` found under `benchmarks/*/results/comparison/`. "
            "Run `benchmarks/lib/compare_to_truth.sh <bench.config.sh>` first._"
        )
    else:
        bench_selector = mo.ui.dropdown(
            options=list(tsv_by_bench.keys()),
            value=next(iter(tsv_by_bench)),
            label="benchmark",
        )
        view = bench_selector
    view
    return (bench_selector,)


@app.cell
def _(bench_selector, csv, mo, tsv_by_bench):
    mo.stop(bench_selector is None, mo.md(""))
    tsv_path = tsv_by_bench[bench_selector.value]
    rows = []
    with tsv_path.open() as _fh:
        for _r in csv.DictReader(_fh, delimiter="\t"):
            rows.append({
                "caller": _r["caller"],
                "class": _r["class"],
                "tp": int(_r["tp"]),
                "fp": int(_r["fp"]),
                "fn": int(_r["fn"]),
                "precision": float(_r["precision"]),
                "recall": float(_r["recall"]),
                "f1": float(_r["f1"]),
            })
    callers = sorted({r["caller"] for r in rows})
    classes = sorted({r["class"] for r in rows}, reverse=True)  # snps before indels
    mo.md(f"**{bench_selector.value}** — {len(rows)} rows · callers: "
          + ", ".join(f"`{c}`" for c in callers))
    return callers, classes, rows, tsv_path


@app.cell
def _():
    import bisect
    import math

    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    # Stable colour per caller, shared with the perf dashboard's intent.
    PALETTE = {
        "ours":      "#1f77b4",  # blue
        "gatk":      "#9467bd",  # purple
        "freebayes": "#d62728",  # red
    }

    def _bar_labels(ax, bars):
        for b in bars:
            h = b.get_height()
            ax.annotate(f"{h:.2f}", (b.get_x() + b.get_width() / 2, h),
                        ha="center", va="bottom", fontsize=7,
                        xytext=(0, 1), textcoords="offset points")

    def plot_metrics(rows, callers, classes):
        """One chart per (class, metric): a row of three charts
        (precision, recall, F1) for each variant class. Bars are one per
        caller, coloured by caller; Y fixed to [0, 1] so every panel is
        directly comparable."""
        metrics = ("precision", "recall", "f1")
        fig, axes = plt.subplots(
            len(classes), len(metrics),
            figsize=(4.6 * len(metrics), 4.2 * len(classes)),
            squeeze=False,
        )
        xs = list(range(len(callers)))
        colours = [PALETTE.get(c, "#555555") for c in callers]
        for row_i, cls in enumerate(classes):
            by_caller = {r["caller"]: r for r in rows if r["class"] == cls}
            for col_i, metric in enumerate(metrics):
                ax = axes[row_i][col_i]
                heights = [by_caller.get(c, {}).get(metric, 0) for c in callers]
                bars = ax.bar(xs, heights, color=colours, width=0.6)
                _bar_labels(ax, bars)
                ax.set_xticks(xs)
                ax.set_xticklabels(callers)
                ax.set_ylim(0, 1.08)
                ax.set_title(f"{cls} — {metric}")
                ax.set_ylabel(metric)
                ax.grid(True, axis="y", alpha=0.3)
        fig.suptitle("Precision / recall / F1 by caller and variant class",
                     fontsize=13, y=1.0)
        fig.tight_layout()
        return fig

    def plot_pr_scatter(rows, classes):
        """Precision-recall scatter; one marker per (caller, class).
        The top-right corner is ideal. Marker shape distinguishes class,
        colour distinguishes caller."""
        markers = {cls: m for cls, m in zip(classes, ("o", "s", "^", "D"))}
        fig, ax = plt.subplots(figsize=(6, 6))
        for r in rows:
            ax.scatter(r["recall"], r["precision"],
                       color=PALETTE.get(r["caller"], "#555555"),
                       marker=markers.get(r["class"], "o"), s=90,
                       edgecolor="black", linewidth=0.4, zorder=3)
            ax.annotate(f"{r['caller']}/{r['class']}",
                        (r["recall"], r["precision"]),
                        fontsize=7, xytext=(4, 4), textcoords="offset points")
        ax.set_xlabel("recall")
        ax.set_ylabel("precision")
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.02)
        ax.grid(True, alpha=0.3)
        ax.set_title("Precision vs recall (top-right is best)")
        fig.tight_layout()
        return fig

    def plot_qual(qual_rows, cls, n_bins, qual_max, density, logy):
        """One chart per caller, all sharing the same x-axis span and
        bin count so distributions line up. Each chart overlays two
        histograms of the caller's own QUAL values for class `cls`: TP
        (calls in the reference) vs FP (calls not in it). Values above
        `qual_max` are clipped into the last bin so nothing is dropped.
        Density mode normalises each histogram to area 1, which lets you
        compare TP/FP *shapes* despite very different counts."""
        callers_q = sorted({r[0] for r in qual_rows if r[1] == cls})
        fig, axes = plt.subplots(
            1, len(callers_q), figsize=(4.8 * len(callers_q), 4.4),
            sharex=True, squeeze=False,
        )
        bin_kw = dict(bins=n_bins, range=(0, qual_max), density=density, alpha=0.55)
        for ax, c in zip(axes[0], callers_q):
            tp = [min(r[3], qual_max) for r in qual_rows
                  if r[0] == c and r[1] == cls and r[2] == "TP"]
            fp = [min(r[3], qual_max) for r in qual_rows
                  if r[0] == c and r[1] == cls and r[2] == "FP"]
            if tp:
                ax.hist(tp, color="#2ca02c", label=f"TP ({len(tp)})", **bin_kw)
            if fp:
                ax.hist(fp, color="#d62728", label=f"FP ({len(fp)})", **bin_kw)
            if logy:
                ax.set_yscale("log")
            ax.set_title(c)
            ax.set_xlabel("QUAL")
            ax.set_ylabel("density" if density else "count")
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)
        fig.suptitle(
            f"{cls} QUAL distribution — TP vs FP per caller "
            f"(shared x: 0–{qual_max:g}, {n_bins} bins)",
            fontsize=13, y=1.03,
        )
        fig.tight_layout()
        return fig

    def _lens_area(d, ra, rb):
        """Area of the intersection (lens) of two circles whose centres
        are distance d apart. Monotonically decreasing in d over
        [|ra-rb|, ra+rb]."""
        if d >= ra + rb:
            return 0.0
        if d <= abs(ra - rb):
            return math.pi * min(ra, rb) ** 2
        a1 = ra * ra * math.acos((d * d + ra * ra - rb * rb) / (2 * d * ra))
        a2 = rb * rb * math.acos((d * d + rb * rb - ra * ra) / (2 * d * rb))
        a3 = 0.5 * math.sqrt(max(0.0,
            (-d + ra + rb) * (d + ra - rb) * (d - ra + rb) * (d + ra + rb)))
        return a1 + a2 - a3

    def _solve_distance(ra, rb, target):
        """Centre distance d such that the lens area equals `target`.
        Bisection — no scipy needed."""
        if ra <= 0 or rb <= 0:
            return ra + rb
        amax = math.pi * min(ra, rb) ** 2
        if target >= amax:
            return abs(ra - rb)      # smaller circle fully inside the larger
        if target <= 0:
            return ra + rb           # tangent / disjoint
        lo, hi = abs(ra - rb), ra + rb
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if _lens_area(mid, ra, rb) > target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def plot_venn(rows, callers, classes):
        """A 2-set Venn per (class, caller): the caller's calls vs the
        reference/truth. Circle *areas* are proportional to set sizes
        (caller = TP+FP, truth = TP+FN) and the overlap area to TP, on a
        single global scale — so circles are comparable across every
        panel (the truth circle is identical for all callers of a class).
        The three regions are the confusion counts: caller-only = FP,
        truth-only = FN, overlap = TP."""
        totals = [v for r in rows
                  for v in (r["tp"] + r["fp"], r["tp"] + r["fn"]) if v > 0]
        max_total = max(totals) if totals else 1
        rmax = 1.0
        k = math.pi * rmax * rmax / max_total          # area per element

        def radius(n):
            return rmax * math.sqrt(n / max_total) if n > 0 else 0.0

        # Pass 1: geometry + a common extent so every panel shares scale.
        geom = {}
        half = 0.1
        for r in rows:
            ra, rb = radius(r["tp"] + r["fp"]), radius(r["tp"] + r["fn"])
            d = _solve_distance(ra, rb, r["tp"] * k)
            geom[(r["caller"], r["class"])] = (ra, rb, -d / 2, d / 2)
            half = max(half, abs(-d / 2) + ra, abs(d / 2) + rb)
        half *= 1.18

        fig, axes = plt.subplots(
            len(classes), len(callers),
            figsize=(4.2 * len(callers), 3.9 * len(classes)),
            squeeze=False,
        )
        for row_i, cls in enumerate(classes):
            by_caller = {r["caller"]: r for r in rows if r["class"] == cls}
            for col_i, c in enumerate(callers):
                ax = axes[row_i][col_i]
                rec = by_caller.get(c)
                if rec is None:
                    ax.axis("off")
                    continue
                ra, rb, xa, xb = geom[(c, cls)]
                colour = PALETTE.get(c, "#555555")
                ax.add_patch(Circle((xa, 0), ra, color=colour, alpha=0.45, lw=0))
                ax.add_patch(Circle((xb, 0), rb, color="#7f7f7f", alpha=0.45, lw=0))
                ax.text(xa, ra + 0.06 * half, c, ha="center", va="bottom",
                        fontsize=9, color=colour, fontweight="bold")
                ax.text(xb, rb + 0.06 * half, "truth", ha="center", va="bottom",
                        fontsize=9, color="#555555", fontweight="bold")
                # Region labels: FP in the caller-only lobe, FN in the
                # truth-only lobe, TP near the lens centroid (rA-rB)/2.
                ax.text(xa - ra * 0.5, 0, f"FP\n{rec['fp']}",
                        ha="center", va="center", fontsize=9)
                ax.text((ra - rb) / 2, 0, f"TP\n{rec['tp']}",
                        ha="center", va="center", fontsize=9, fontweight="bold")
                ax.text(xb + rb * 0.5, 0, f"FN\n{rec['fn']}",
                        ha="center", va="center", fontsize=9)
                ax.set_xlim(-half, half)
                ax.set_ylim(-half, half)
                ax.set_aspect("equal")
                ax.axis("off")
                ax.set_title(f"{c} — {cls}", fontsize=10)
        fig.suptitle("Caller vs truth — variant set overlap "
                     "(areas ∝ counts, shared scale)", fontsize=13, y=1.0)
        fig.tight_layout()
        return fig

    def qual_sweep(tp_quals, fp_quals, n_positives):
        """Sweep a QUAL cutoff over a caller's calls. At each candidate
        threshold t (every distinct QUAL value), keep calls with
        QUAL >= t and recompute the confusion counts:
          TP(t) = TP-labelled calls with QUAL >= t
          FP(t) = FP-labelled calls with QUAL >= t
          FN(t) = n_positives - TP(t)   (truth not recovered at t)
        n_positives = TP+FN total (all reference variants for the class),
        so recall is anchored to the full truth, not just emitted calls.
        Returns the list of operating points (low->high threshold)."""
        tps = sorted(tp_quals)
        fps = sorted(fp_quals)
        candidates = sorted(set(tps) | set(fps))
        points = []
        for t in candidates:
            tp = len(tps) - bisect.bisect_left(tps, t)
            fp = len(fps) - bisect.bisect_left(fps, t)
            prec = tp / (tp + fp) if (tp + fp) > 0 else 1.0
            rec = tp / n_positives if n_positives > 0 else 0.0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
            points.append({"t": t, "tp": tp, "fp": fp,
                           "precision": prec, "recall": rec, "f1": f1})
        return points

    def plot_roc(qual_rows, rows, cls):
        """PR curve (parametric in QUAL) + F1-vs-threshold per caller.
        The dot on the PR curve and the dotted vertical line mark each
        caller's F1-optimal threshold. Returns (fig, best) where best
        maps caller -> its max-F1 operating point."""
        callers_r = sorted({r[0] for r in qual_rows if r[1] == cls})
        positives = {r["caller"]: r["tp"] + r["fn"]
                     for r in rows if r["class"] == cls}
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
        best = {}
        for c in callers_r:
            tpq = [r[3] for r in qual_rows
                   if r[0] == c and r[1] == cls and r[2] == "TP"]
            fpq = [r[3] for r in qual_rows
                   if r[0] == c and r[1] == cls and r[2] == "FP"]
            pts = qual_sweep(tpq, fpq, positives.get(c, len(tpq)))
            if not pts:
                continue
            colour = PALETTE.get(c, "#555555")
            ax1.plot([p["recall"] for p in pts], [p["precision"] for p in pts],
                     color=colour, lw=1.6, label=c)
            ax2.plot([p["t"] for p in pts], [p["f1"] for p in pts],
                     color=colour, lw=1.6, label=c)
            bp = max(pts, key=lambda p: p["f1"])
            best[c] = bp
            ax1.scatter([bp["recall"]], [bp["precision"]], color=colour,
                        s=80, edgecolor="black", linewidth=0.5, zorder=4)
            ax2.axvline(bp["t"], color=colour, ls=":", alpha=0.7)
        ax1.set_xlabel("recall")
        ax1.set_ylabel("precision")
        ax1.set_xlim(0, 1.02)
        ax1.set_ylim(0, 1.02)
        ax1.grid(alpha=0.3)
        ax1.legend(fontsize=8, loc="lower left")
        ax1.set_title(f"{cls} — precision-recall (dot = max-F1 cutoff)")
        # Zoom the threshold axis to the region around the optima — the
        # QUAL tail runs to ~20k and would otherwise squash everything.
        if best:
            xmax = max(2.5 * bp["t"] for bp in best.values()) or 1.0
            ax2.set_xlim(0, xmax)
        ax2.set_xlabel("QUAL threshold")
        ax2.set_ylabel("F1")
        ax2.grid(alpha=0.3)
        ax2.legend(fontsize=8, loc="lower right")
        ax2.set_title(f"{cls} — F1 vs QUAL threshold (dotted = optimum)")
        fig.suptitle("QUAL sweep — finding the best per-caller cutoff",
                     fontsize=13, y=1.01)
        fig.tight_layout()
        return fig, best

    return (PALETTE, plot_metrics, plot_pr_scatter, plot_qual, plot_roc,
            plot_venn, qual_sweep, plt)


@app.cell
def _(callers, classes, mo, plot_metrics, rows):
    mo.stop(not rows, mo.md(""))
    plot_metrics(rows, callers, classes)
    return


@app.cell
def _(classes, mo, plot_pr_scatter, rows):
    mo.stop(not rows, mo.md(""))
    plot_pr_scatter(rows, classes)
    return


@app.cell
def _(mo):
    mo.md("""
    ---
    ## Caller vs truth — set overlap (Venn)

    Each caller's call set against the Genome-in-a-Bottle reference, one
    Venn per caller (rows = variant class). The regions are the
    confusion counts: the caller-only slice is **FP**, the truth-only
    slice is **FN**, and the overlap is **TP**.
    """)
    return


@app.cell
def _(callers, classes, mo, plot_venn, rows):
    mo.stop(not rows, mo.md(""))
    plot_venn(rows, callers, classes)
    return


@app.cell
def _(mo):
    mo.md("""
    ---
    ## Quality distributions — TP vs FP per caller

    For each caller, the QUAL distribution of its own calls, split into
    true positives (in the reference) and false positives (not in it).
    All charts share the same x span and bin count so they line up. A
    well-calibrated caller separates the two: FPs concentrated at low
    QUAL, TPs at high — so a QUAL cutoff can drop FPs without losing TPs.
    """)
    return


@app.cell
def _(bench_selector, csv, mo, tsv_by_bench):
    # Per-record QUAL labelled TP/FP, written next to accuracy.tsv.
    mo.stop(bench_selector is None, mo.md(""))
    qual_path = tsv_by_bench[bench_selector.value].with_name("qual_dist.tsv")
    qual_rows = []
    if qual_path.exists():
        with qual_path.open() as _fh:
            for _r in csv.DictReader(_fh, delimiter="\t"):
                qual_rows.append(
                    (_r["caller"], _r["class"], _r["status"], float(_r["qual"]))
                )
    return (qual_rows,)


@app.cell
def _(mo, qual_rows):
    mo.stop(
        not qual_rows,
        mo.md("_No `qual_dist.tsv` for this benchmark — re-run "
              "`compare_to_truth.sh` to regenerate it._"),
    )
    q_classes = sorted({r[1] for r in qual_rows}, reverse=True)
    _allq = sorted(r[3] for r in qual_rows)
    _p99 = _allq[min(len(_allq) - 1, int(0.99 * len(_allq)))]
    _gmax = max(_allq)

    class_sel = mo.ui.radio(
        options=q_classes,
        value="snps" if "snps" in q_classes else q_classes[0],
        label="variant class",
    )
    bins_sel = mo.ui.slider(10, 120, value=50, label="bins")
    qmax_sel = mo.ui.slider(
        10, int(round(_gmax)), value=float(round(_p99)),
        step=max(1, int(_gmax / 200)), label="max QUAL (x-axis clip)",
    )
    density_sel = mo.ui.checkbox(value=True, label="density (compare shapes)")
    logy_sel = mo.ui.checkbox(value=False, label="log y")
    controls = mo.hstack(
        [class_sel, bins_sel, qmax_sel, density_sel, logy_sel], wrap=True,
    )
    controls
    return bins_sel, class_sel, density_sel, logy_sel, qmax_sel


@app.cell
def _(bins_sel, class_sel, density_sel, logy_sel, mo, plot_qual, qmax_sel, qual_rows):
    mo.stop(not qual_rows, mo.md(""))
    plot_qual(
        qual_rows, class_sel.value, bins_sel.value, qmax_sel.value,
        density_sel.value, logy_sel.value,
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ---
    ## QUAL sweep — best per-caller cutoff (ROC-style)

    Sweeping a QUAL threshold over each caller's calls traces a
    precision-recall curve: raising the cutoff drops low-QUAL calls,
    trading recall for precision. The dot marks each caller's
    **F1-optimal** threshold; the right panel shows F1 vs threshold
    directly. Recall is anchored to the full reference (TP+FN), so it
    reflects truth missed, not just calls dropped. Uses the variant
    class selected above.
    """)
    return


@app.cell
def _(class_sel, mo, plot_roc, qual_rows, rows):
    mo.stop(not qual_rows, mo.md(""))
    roc_fig, roc_best = plot_roc(qual_rows, rows, class_sel.value)
    roc_fig
    return (roc_best,)


@app.cell
def _(class_sel, mo, roc_best):
    mo.stop(not roc_best, mo.md(""))
    _md = f"### F1-optimal QUAL threshold — {class_sel.value}\n\n"
    _md += "| caller | best QUAL | TP | FP | precision | recall | F1 |\n"
    _md += "|---|---:|---:|---:|---:|---:|---:|\n"
    for _c, _bp in sorted(roc_best.items()):
        _md += (f"| {_c} | {_bp['t']:.1f} | {_bp['tp']} | {_bp['fp']} | "
                f"{_bp['precision']:.4f} | {_bp['recall']:.4f} | "
                f"{_bp['f1']:.4f} |\n")
    mo.md(_md)
    return


@app.cell
def _(mo, rows):
    # Raw counts table.
    mo.stop(not rows, mo.md(""))
    _md = "### Raw counts\n\n"
    _md += "| caller | class | TP | FP | FN | precision | recall | F1 |\n"
    _md += "|---|---|---:|---:|---:|---:|---:|---:|\n"
    for _r in sorted(rows, key=lambda x: (x["class"], x["caller"])):
        _md += (f"| {_r['caller']} | {_r['class']} | {_r['tp']} | {_r['fp']} | "
                f"{_r['fn']} | {_r['precision']:.4f} | {_r['recall']:.4f} | "
                f"{_r['f1']:.4f} |\n")
    mo.md(_md)
    return


if __name__ == "__main__":
    app.run()
