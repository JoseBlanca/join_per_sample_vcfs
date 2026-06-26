# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "polars",
#     "altair",
#     "pyarrow",
# ]
# ///
"""GIAB 2D TP/FP separation dashboard — MAPQ-difference vs depth / allele balance.

Idea: the MAPQ-difference filter alone cannot separate true SNPs from
paralog/artefact false positives (the t and MAPQ-delta distributions overlap).
Can a SECOND axis pull them apart? This explores the joint 2D distribution of
the MAPQ signal against depth and allele balance (VAF), labelled TP vs FP.

Data: benchmarks/giab/src/signal_table.sh writes, per SNP call,
  status (TP/FP), dp, mqref, mqalt, mqdiff (=MQAlt-MQRef, depth-independent),
  mqdifft (Welch's t, depth-sensitive), af, ref_ad, alt_ad, vaf
into results/per_sample/<cov>/<variant>/signal_table.tsv. Run that first:

    benchmarks/giab/src/run_ours_per_sample.sh 300x        # PRESET=high-recall
    benchmarks/giab/src/signal_table.sh high-recall
    uv run marimo run benchmarks/giab/src/mapq_depth_dashboard.py

Findings baked into the defaults: depth does NOT separate the classes (these
benchmark BAMs are coverage-normalised, so paralog pile-up is hidden), but VAF
does — real hets sit near 0.5 (homs near 1.0) while artefacts cluster at
~0.2 VAF. MAPQ-delta adds a negative tail for the genuine paralog FPs.
"""

import marimo

__generated_with = "0.16.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import os
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import polars as pl

    BENCH_DIR = Path(
        os.environ.get("PVC_GIAB_BENCH_DIR", Path(__file__).resolve().parents[1])
    )
    RESULTS_DIR = BENCH_DIR / "results" / "per_sample"

    # Axis menus. label -> (column, human title).
    X_AXES = {
        "MQ-diff (MQAlt−MQRef, depth-independent)": ("mqdiff", "MAPQ difference (alt − ref)"),
        "MQ-diff Welch's t (depth-sensitive)": ("mqdifft", "MAPQ-diff Welch's t"),
        "mean ALT MAPQ": ("mqalt", "mean ALT MAPQ"),
    }
    Y_AXES = {
        "VAF (allele balance)": ("vaf", "VAF = alt / (ref+alt)"),
        "depth deviation (DP / sample median)": ("depth_dev", "depth deviation"),
        "depth (DP)": ("dp", "locus depth"),
        "mean ALT MAPQ": ("mqalt", "mean ALT MAPQ"),
    }
    return RESULTS_DIR, X_AXES, Y_AXES, Path, alt, mo, pl


@app.cell
def _(RESULTS_DIR, mo):
    covs = sorted(p.name for p in RESULTS_DIR.iterdir() if p.is_dir()) if RESULTS_DIR.is_dir() else []
    mo.stop(not covs, mo.md(f"No results under `{RESULTS_DIR}`."))
    # Variants that actually have a signal_table.tsv.
    pairs = []
    for cov in covs:
        for vdir in sorted((RESULTS_DIR / cov).iterdir()):
            if (vdir / "signal_table.tsv").exists():
                pairs.append(f"{cov}/{vdir.name}")
    mo.stop(
        not pairs,
        mo.md("No `signal_table.tsv` found. Run `benchmarks/giab/src/signal_table.sh <variant>` first."),
    )
    src_dd = mo.ui.dropdown(
        options=pairs,
        value=next((p for p in pairs if p.endswith("high-recall")), pairs[0]),
        label="coverage / variant",
    )
    return (src_dd,)


@app.cell
def _(RESULTS_DIR, mo, pl, src_dd):
    tsv = RESULTS_DIR / src_dd.value / "signal_table.tsv"
    df = pl.read_csv(tsv, separator="\t", null_values=["", "."])
    # Drop rows lacking VAF (no AD), and add per-sample depth deviation.
    df = df.filter(pl.col("vaf").is_not_null())
    df = df.with_columns(
        (pl.col("dp") / pl.col("dp").median().over("sample")).alias("depth_dev")
    )
    n_tp = df.filter(pl.col("status") == "TP").height
    n_fp = df.filter(pl.col("status") == "FP").height
    hdr = mo.md(
        f"## 2D TP/FP separation — `{src_dd.value}`\n\n"
        f"{df.height} SNP calls · **{n_tp} TP**, **{n_fp} FP**. "
        "Each point is one called SNP, coloured by truth status."
    )
    return df, hdr, n_fp, n_tp


@app.cell
def _(X_AXES, Y_AXES, mo):
    x_dd = mo.ui.dropdown(options=list(X_AXES), value=list(X_AXES)[0], label="X axis")
    y_dd = mo.ui.dropdown(options=list(Y_AXES), value=list(Y_AXES)[0], label="Y axis")
    return x_dd, y_dd


@app.cell
def _(X_AXES, Y_AXES, x_dd, y_dd):
    # Resolve the selected axes once; both chart cells consume these.
    xcol, xtitle = X_AXES[x_dd.value]
    ycol, ytitle = Y_AXES[y_dd.value]
    return xcol, xtitle, ycol, ytitle


@app.cell
def _(alt, df, hdr, mo, x_dd, xcol, xtitle, y_dd, ycol, ytitle):
    STATUS = dict(sort=["TP", "FP"], scale=alt.Scale(domain=["TP", "FP"],
                                                      range=["#2ca02c", "#d62728"]))
    base = df.select(["sample", "status", xcol, ycol]).drop_nulls()

    pts = (
        alt.Chart(base).mark_circle(size=34, opacity=0.45).encode(
            x=alt.X(f"{xcol}:Q", title=xtitle, scale=alt.Scale(zero=False)),
            y=alt.Y(f"{ycol}:Q", title=ytitle, scale=alt.Scale(zero=False)),
            color=alt.Color("status:N", title="status", **STATUS),
            tooltip=["sample", "status", f"{xcol}:Q", f"{ycol}:Q"],
        ).properties(width=520, height=420)
    )
    top = (
        alt.Chart(base).mark_area(opacity=0.45, interpolate="step").encode(
            x=alt.X(f"{xcol}:Q", bin=alt.Bin(maxbins=40), title=None, axis=None),
            y=alt.Y("count():Q", title=None, stack=None, axis=None),
            color=alt.Color("status:N", **STATUS, legend=None),
        ).properties(width=520, height=70)
    )
    right = (
        alt.Chart(base).mark_area(opacity=0.45, interpolate="step").encode(
            y=alt.Y(f"{ycol}:Q", bin=alt.Bin(maxbins=40), title=None, axis=None),
            x=alt.X("count():Q", title=None, stack=None, axis=None),
            color=alt.Color("status:N", **STATUS, legend=None),
        ).properties(width=70, height=420)
    )
    # Static render (NOT mo.ui.altair_chart): a concatenated spec + marimo's
    # auto legend-selection would collide on a duplicate signal name.
    joint = alt.vconcat(top, alt.hconcat(pts, right, spacing=4), spacing=4).resolve_scale(
        color="shared"
    )
    mo.vstack([
        hdr,
        mo.hstack([x_dd, y_dd], justify="start", gap=2),
        joint,
        mo.md(
            "_Marginal histograms on top/right show each axis's 1D TP/FP overlap; "
            "the scatter shows whether the **pair** separates what neither does alone._"
        ),
    ])
    return (STATUS,)


@app.cell
def _(alt, df, mo, pl, xcol, xtitle, ycol, ytitle):
    # Purity heatmap: 2D bins coloured by fraction TP (green=pure TP, red=pure FP).
    b = (
        df.select(["status", xcol, ycol])
        .drop_nulls()
        .with_columns((pl.col("status") == "TP").cast(pl.Int8).alias("is_tp"))
    )
    heat = (
        alt.Chart(b).mark_rect().encode(
            x=alt.X(f"{xcol}:Q", bin=alt.Bin(maxbins=24), title=xtitle),
            y=alt.Y(f"{ycol}:Q", bin=alt.Bin(maxbins=24), title=ytitle),
            color=alt.Color("mean(is_tp):Q", title="fraction TP",
                            scale=alt.Scale(scheme="redyellowgreen", domain=[0, 1])),
            tooltip=[alt.Tooltip("count():Q", title="n"),
                     alt.Tooltip("mean(is_tp):Q", title="fraction TP", format=".2f")],
        ).properties(width=520, height=420, title="Cell purity (green = TP-pure, red = FP-pure)")
    )
    mo.vstack([mo.md("### Cell purity map"), mo.ui.altair_chart(heat)])
    return


@app.cell
def _(mo):
    vaf_floor = mo.ui.slider(0.0, 0.5, step=0.01, value=0.35, label="keep VAF ≥", show_value=True)
    mqdiff_floor = mo.ui.slider(-30.0, 5.0, step=0.5, value=-30.0,
                                label="keep MQ-diff ≥ (−30 = off)", show_value=True)
    return mqdiff_floor, vaf_floor


@app.cell
def _(df, mo, mqdiff_floor, pl, vaf_floor):
    # Decision-region explorer: a candidate combined filter on VAF + MQ-diff.
    vf, mf = vaf_floor.value, mqdiff_floor.value
    kept = (pl.col("vaf") >= vf) & (pl.col("mqdiff").fill_null(0.0) >= mf)
    g = df.with_columns(kept.alias("keep"))
    tp_keep = g.filter((pl.col("status") == "TP") & pl.col("keep")).height
    fp_keep = g.filter((pl.col("status") == "FP") & pl.col("keep")).height
    tp_drop = g.filter((pl.col("status") == "TP") & ~pl.col("keep")).height
    fp_drop = g.filter((pl.col("status") == "FP") & ~pl.col("keep")).height
    prec = tp_keep / (tp_keep + fp_keep) if (tp_keep + fp_keep) else 0.0
    # Recall here is relative to TP present in this called set (not truth FN).
    rec = tp_keep / (tp_keep + tp_drop) if (tp_keep + tp_drop) else 0.0

    mo.vstack([
        mo.md("### Combined-filter explorer  ·  keep if VAF ≥ v **and** MQ-diff ≥ m"),
        mo.hstack([vaf_floor, mqdiff_floor], justify="start", gap=2),
        mo.md(
            f"**Kept** {tp_keep} TP / {fp_keep} FP  ·  **dropped** {tp_drop} TP / {fp_drop} FP\n\n"
            f"precision of kept set = **{prec:.3f}**  ·  "
            f"TP retained = **{rec:.3f}**  ·  "
            f"FP removed = **{(fp_drop / (fp_keep + fp_drop) if (fp_keep+fp_drop) else 0):.3f}**"
        ),
        mo.md(
            "_A VAF floor alone (MQ-diff off at −30) already removes most FPs at a "
            "small TP cost — the allele-balance signal the MAPQ-diff filter lacks. "
            "Add a MQ-diff floor to also catch the genuine paralog tail at high VAF._"
        ),
    ])
    return


@app.cell
def _(mo):
    mo.md(
        """
        ---
        **Reading guide.** `mqdiff` = MQAlt − MQRef is the *depth-independent*
        effect size; `mqdifft` is the Welch's t that the production filter
        thresholds (depth-sensitive — the reason it misfires at 300×). `vaf` is
        the observed allele balance from FORMAT/AD. TN omitted; recall in the
        explorer is over called TPs only (truth-level FN are separate). These are
        all SNP calls within each sample's confident BED, labelled by isec vs the
        GIAB truth.
        """
    )
    return


if __name__ == "__main__":
    app.run()
