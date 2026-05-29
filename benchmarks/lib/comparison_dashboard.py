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
    import matplotlib.pyplot as plt

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
        """Grouped bars: precision / recall / F1 per caller, one panel
        per variant class. Y fixed to [0, 1] so callers and classes are
        visually comparable across panels."""
        metrics = ("precision", "recall", "f1")
        fig, axes = plt.subplots(
            1, len(classes), figsize=(6.5 * len(classes), 5), squeeze=False,
        )
        width = 0.25
        for ax, cls in zip(axes[0], classes):
            by_caller = {r["caller"]: r for r in rows if r["class"] == cls}
            xs = list(range(len(callers)))
            for i, metric in enumerate(metrics):
                offset = (i - 1) * width
                heights = [by_caller.get(c, {}).get(metric, 0) for c in callers]
                bars = ax.bar([x + offset for x in xs], heights, width,
                              label=metric)
                _bar_labels(ax, bars)
            ax.set_xticks(xs)
            ax.set_xticklabels(callers)
            ax.set_ylim(0, 1.08)
            ax.set_ylabel("score")
            ax.set_title(cls)
            ax.grid(True, axis="y", alpha=0.3)
            ax.legend(fontsize=9, loc="lower right")
        fig.suptitle("Precision / recall / F1 by caller and variant class",
                     fontsize=13, y=1.02)
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

    return PALETTE, plot_metrics, plot_pr_scatter, plt


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
