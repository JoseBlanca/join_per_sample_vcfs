# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "pandas", "matplotlib"]
# ///
"""FP / FN vs coverage — our SSR caller vs HipSTR, scored against the GIAB HG002
Tandem-Repeat assembly truth, on a single sample over a coverage ladder (5x..300x).

Universe: the 13,272 restricted-catalog loci both callers genotype. Truth is
assembly-based (orthogonal to short reads). Data comes from results/eval_table.tsv
(built by src/build_eval_table.py). Single-sample calling of our caller relies on
the recurrence_k=min(2,n_samples) graceful-degradation fix (see PROJECT_STATUS).

Run:  uv run marimo run   benchmarks/ssr_hg002/scripts/fp_fn_vs_coverage_dashboard.py
Edit: uv run marimo edit  benchmarks/ssr_hg002/scripts/fp_fn_vs_coverage_dashboard.py
"""
import marimo

app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    from pathlib import Path
    import pandas as pd
    import matplotlib.pyplot as plt
    return Path, mo, pd, plt


@app.cell
def _(mo):
    mo.md(
        """
        # SSR calling: FP / FN vs coverage — ours vs HipSTR vs freebayes (GIAB HG002)

        Single sample **HG002**, genotyped at a coverage ladder against the
        **GIAB Tandem-Repeat v1.0.1** assembly truth. Universe = the **13,272**
        restricted-catalog loci. A locus is a **truth positive** when HG002 carries
        a repeat-unit length change there.

        - **FP** = caller calls a length variant where truth is non-variant.
        - **FN** = caller misses a truth length variant (no-call *or* ref-call).
        - Precision = TP/(TP+FP), Recall = TP/(TP+FN).

        *Callers:* **ours** (cohort caller run single-sample via the
        `recurrence_k=min(2,n_samples)` fallback — its floor, not its cohort
        strength); **HipSTR** (STR-specialised, `--use-unpaired`); **freebayes**
        (general Bayesian caller, NOT STR-aware, indels QUAL≥20 overlap-matched).
        """
    )
    return


@app.cell
def _(Path, mo, pd):
    table = Path(__file__).resolve().parent.parent / "results" / "eval_table.tsv"
    mo.stop(
        not table.exists(),
        mo.md(f"**Missing** `{table}`. Build it:\n\n`uv run benchmarks/ssr_hg002/src/build_eval_table.py`"),
    )
    df = pd.read_csv(table, sep="\t")
    _names = {2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}
    df["period_name"] = df["period"].map(lambda p: _names.get(p, f"{p}-mer"))
    return (df,)


@app.cell
def _(df, mo):
    _periods = ["all"] + sorted(df["period_name"].unique(), key=lambda x: len(x))
    period_sel = mo.ui.dropdown(_periods, value="all", label="Motif period")
    truth_sel = mo.ui.radio(
        {"repeat-unit change (period-aware)": "truth_var", "any indel": "truth_var_any"},
        value="repeat-unit change (period-aware)",
        label="Truth-positive definition",
    )
    mo.hstack([period_sel, truth_sel], justify="start", gap=2)
    return period_sel, truth_sel


@app.cell
def _(df, pd, period_sel, truth_sel):
    def compute(frame, period, truth_col):
        sub = frame if period == "all" else frame[frame["period_name"] == period]

        def metrics(g, caller):
            t = g[truth_col] == 1
            v = g[caller] == "var"
            nc = g[caller] == "nocall"
            TP = int((t & v).sum()); FP = int((~t & v).sum())
            FN = int((t & ~v).sum()); TN = int((~t & ~v).sum())
            prec = TP / (TP + FP) if TP + FP else float("nan")
            rec = TP / (TP + FN) if TP + FN else float("nan")
            f1 = 2 * prec * rec / (prec + rec) if (prec and rec) else float("nan")
            fpr = FP / (FP + TN) if FP + TN else float("nan")
            return dict(TP=TP, FP=FP, FN=FN, TN=TN, precision=prec, recall=rec,
                        f1=f1, fp_rate=fpr, nocall=int(nc.sum()),
                        nocall_frac=float(nc.mean()))

        recs, states = [], []
        for cov, g in sub.groupby("cov"):
            for caller in ("ours", "hipstr", "freebayes"):
                recs.append({"cov": int(cov), "caller": caller, **metrics(g, caller)})
                vc = g[caller].value_counts()
                states.append({"cov": int(cov), "caller": caller,
                               "var": int(vc.get("var", 0)), "ref": int(vc.get("ref", 0)),
                               "nocall": int(vc.get("nocall", 0))})
        met = pd.DataFrame(recs).sort_values(["caller", "cov"])
        state = pd.DataFrame(states).sort_values(["caller", "cov"])
        n_loci = len(sub) // max(sub["cov"].nunique(), 1)
        deepest = sub[sub["cov"] == sub["cov"].max()]
        truth_pos = int((deepest[truth_col] == 1).sum())
        return met, state, n_loci, truth_pos

    met, state, n_loci, truth_pos = compute(df, period_sel.value, truth_sel.value)
    return met, n_loci, state, truth_pos


@app.cell
def _(met, mo, n_loci, period_sel, truth_pos):
    mo.md(
        f"**{n_loci:,} loci** in scope (period = {period_sel.value}); "
        f"**{truth_pos:,} truth positives**. Metrics per coverage:"
    )
    return


@app.cell
def _(met, mo):
    def fmt(m):
        s = m.copy()
        for c in ("precision", "recall", "f1", "fp_rate", "nocall_frac"):
            s[c] = s[c].map(lambda x: f"{x:.3f}")
        return s

    mo.ui.table(fmt(met), selection=None, pagination=False)
    return


@app.cell
def _(met, plt):
    def make_fig(m):
        colors = {"ours": "#2b6cb0", "hipstr": "#dd6b20", "freebayes": "#38a169"}
        panels = [("recall", "Recall (sensitivity)", False),
                  ("precision", "Precision", False),
                  ("FN", "False negatives (count)", True),
                  ("FP", "False positives (count)", True)]
        fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
        for ax, (col, title, is_count) in zip(axes.flat, panels):
            for caller in ("ours", "hipstr", "freebayes"):
                d = m[m["caller"] == caller].sort_values("cov")
                ax.plot(d["cov"], d[col], "o-", color=colors[caller], label=caller)
            ax.set_title(title); ax.set_xlabel("coverage (x)")
            ax.set_xscale("log"); ax.set_xticks(sorted(m["cov"].unique()))
            ax.get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
            ax.grid(True, alpha=0.3)
            if not is_count:
                ax.set_ylim(0, 1.02)
            ax.legend()
        fig.suptitle("Ours vs HipSTR vs freebayes (GIAB truth) — depth degradation", fontweight="bold")
        fig.tight_layout()
        return fig

    make_fig(met)
    return


@app.cell
def _(mo, state):
    mo.vstack([mo.md("### Call-state breakdown (var / ref / no-call) per coverage"),
               mo.ui.table(state, selection=None, pagination=False)])
    return


if __name__ == "__main__":
    app.run()
