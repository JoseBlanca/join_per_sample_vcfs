# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "polars",
#     "altair",
#     "pyarrow",
#     "numpy",
#     "scipy",
# ]
# ///
"""Allele-balance separation vs coverage — the production-viability test.

The depth-independent finding at 300× was that VAF separates TP from FP. The
production question is whether a DEPTH-AWARE allele-balance statistic holds up
across 5×–300×: aggressive where it should be (high cov, many FP) and quiet
where the data can't support it (low cov, VAF fluctuates).

Statistic: condition on the called genotype (GT). A het should sit at VAF≈0.5,
a hom-alt at ≈1.0. We model the alt count k out of n=ref_ad+alt_ad as
Beta-Binomial(n, mean=p0(GT), concentration s), where s (overdispersion) is the
ONE knob that makes it depth-honest — estimated by MLE from high-coverage TP
hets, not guessed. The per-call score is a log-likelihood ratio:

    AB_LR = logBB(k | n, mean=p0)  −  logBB(k | n, mean=k/n)

≈ 0  → counts consistent with the called genotype's balance (keep)
≪ 0  → counts inconsistent with any real genotype balance (artifact)

By construction AB_LR is depth-scaled: at low n the BB is wide, so even a
skewed k is not far from p0 in likelihood → AB_LR stays near 0 → no false drop.

Inputs: results/per_sample/<cov>/high-recall/signal_table.tsv
(from signal_table.sh, one row per SNP call with gt, ref_ad, alt_ad, status).

    for c in 5x 10x 15x 20x 30x 50x 300x; do
        benchmarks/giab/src/signal_table.sh high-recall $c; done
    uv run marimo run benchmarks/giab/src/allele_balance_dashboard.py
"""

import marimo

__generated_with = "0.16.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import os
    import re
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import numpy as np
    import polars as pl
    from scipy.special import betaln, gammaln

    BENCH_DIR = Path(
        os.environ.get("PVC_GIAB_BENCH_DIR", Path(__file__).resolve().parents[1])
    )
    RESULTS = BENCH_DIR / "results" / "per_sample"
    VARIANT = os.environ.get("PVC_AB_VARIANT", "high-recall")
    HOM_ALT_MEAN = 0.98  # expected VAF for a 1/1 call (leave room for a little error)
    return (
        HOM_ALT_MEAN,
        Path,
        RESULTS,
        VARIANT,
        alt,
        betaln,
        gammaln,
        mo,
        np,
        pl,
        re,
    )


@app.cell
def _(RESULTS, VARIANT, mo, pl, re):
    # Load every coverage tier's signal_table for the chosen variant.
    frames = []
    if RESULTS.is_dir():
        for cov_dir in sorted(RESULTS.iterdir()):
            tsv = cov_dir / VARIANT / "signal_table.tsv"
            if not tsv.exists():
                continue
            m = re.match(r"(\d+)x", cov_dir.name)
            cov_num = int(m.group(1)) if m else 0
            # infer_schema_length=None: scan all rows (mqdifft is int in early
            # rows, float later — partial inference mistypes it as i64).
            d = pl.read_csv(tsv, separator="\t", null_values=["", "."],
                            infer_schema_length=None)
            d = d.with_columns(
                pl.lit(cov_dir.name).alias("coverage"),
                pl.lit(cov_num).alias("cov_x"),
            )
            frames.append(d)
    mo.stop(
        not frames,
        mo.md(f"No `signal_table.tsv` for variant `{VARIANT}` under `{RESULTS}`. "
              "Run `signal_table.sh high-recall <cov>` for each coverage first."),
    )
    raw = pl.concat(frames, how="vertical_relaxed")
    # Keep calls with usable counts.
    raw = raw.with_columns((pl.col("ref_ad") + pl.col("alt_ad")).alias("n")).filter(
        pl.col("n") > 0
    )
    return (raw,)


@app.cell
def _(HOM_ALT_MEAN, betaln, gammaln, np, raw):
    # --- Beta-Binomial helpers (vectorised over numpy arrays) ---------------
    def bb_logpmf(k, n, mean, s):
        a = mean * s
        b = (1.0 - mean) * s
        return (
            gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
            + betaln(k + a, n - k + b) - betaln(a, b)
        )

    k_all = raw["alt_ad"].to_numpy().astype(float)
    n_all = raw["n"].to_numpy().astype(float)
    is_het = (raw["gt"].to_numpy() == "0/1")
    p0 = np.where(is_het, 0.5, HOM_ALT_MEAN)

    # --- Estimate overdispersion s by MLE on the highest-coverage TP hets ---
    cov_x = raw["cov_x"].to_numpy()
    status = raw["status"].to_numpy()
    top_cov = cov_x.max()
    fit_mask = (cov_x == top_cov) & (status == "TP") & is_het
    kf, nf = k_all[fit_mask], n_all[fit_mask]
    s_grid = np.array([5, 10, 20, 50, 100, 150, 200, 300, 500, 800, 1500, 3000], float)
    ll_by_s = np.array([bb_logpmf(kf, nf, 0.5, s).sum() for s in s_grid])
    s_hat = float(s_grid[int(np.argmax(ll_by_s))])

    # --- Per-call AB log-likelihood ratio (genotype balance vs free fit) ----
    p_free = np.clip(k_all / n_all, 1e-6, 1 - 1e-6)
    ab_lr = bb_logpmf(k_all, n_all, p0, s_hat) - bb_logpmf(k_all, n_all, p_free, s_hat)

    df = raw.with_columns(
        ab_lr=ab_lr,
        p0=p0,
        is_het=is_het,
    )
    return df, s_grid, s_hat, top_cov


@app.cell
def _(VARIANT, mo, s_hat, top_cov):
    mo.md(
        f"""
        # Allele-balance separation vs coverage — `{VARIANT}`

        Overdispersion **s = {s_hat:.0f}** (MLE on {top_cov}× TP hets). Larger s →
        tighter around the genotype balance; this is the depth-honesty knob.
        The score is **AB_LR** = logBB(k | genotype balance) − logBB(k | free fit):
        ~0 means the alt/ref split fits the called genotype; strongly negative
        means it fits no real genotype (artifact).
        """
    )
    return


@app.cell
def _(alt, df, mo, pl):
    # Panel A — the FP problem is concentrated at high coverage.
    counts = (
        df.group_by(["coverage", "cov_x", "status"]).len()
        .sort("cov_x")
    )
    chart_a = (
        alt.Chart(counts).mark_line(point=True).encode(
            x=alt.X("cov_x:Q", title="coverage (×)", scale=alt.Scale(type="log")),
            y=alt.Y("len:Q", title="calls"),
            color=alt.Color("status:N", sort=["TP", "FP"],
                            scale=alt.Scale(domain=["TP", "FP"], range=["#2ca02c", "#d62728"])),
            tooltip=["coverage", "status", "len"],
        ).properties(width=560, height=300, title="TP / FP call counts vs coverage")
    )
    mo.vstack([
        mo.md("## A. Where the FPs are — count vs coverage"),
        mo.ui.altair_chart(chart_a),
        mo.md("_FPs explode with depth (low-VAF artifacts only become callable as "
              "alt reads accumulate); TPs plateau. So the filter is needed most at high cov._"),
    ])
    return


@app.cell
def _(df, mo):
    cov_pick = mo.ui.dropdown(
        options=sorted(df["coverage"].unique().to_list(),
                       key=lambda c: int(c.rstrip("x"))),
        value="300x" if "300x" in df["coverage"].to_list() else None,
        label="coverage for distribution",
    )
    return (cov_pick,)


@app.cell
def _(alt, cov_pick, df, mo, pl):
    # Panel B — AB_LR distribution TP vs FP at the selected coverage.
    sel = df.filter(pl.col("coverage") == cov_pick.value)
    # clip for readability (long negative tail)
    selc = sel.with_columns(pl.col("ab_lr").clip(-60, 5).alias("ab_lr_clip"))
    hist = (
        alt.Chart(selc).mark_bar(opacity=0.55).encode(
            x=alt.X("ab_lr_clip:Q", bin=alt.Bin(maxbins=50), title="AB_LR (clipped at −60)"),
            y=alt.Y("count():Q", stack=None, title="calls"),
            color=alt.Color("status:N", sort=["TP", "FP"],
                            scale=alt.Scale(domain=["TP", "FP"], range=["#2ca02c", "#d62728"])),
            tooltip=["status:N", alt.Tooltip("count():Q")],
        ).properties(width=560, height=300,
                     title=f"AB_LR by status — {cov_pick.value} (left = artifact-like)")
    )
    mo.vstack([
        mo.md("## B. AB_LR distribution by status"),
        cov_pick,
        mo.ui.altair_chart(hist),
    ])
    return


@app.cell
def _(df, mo, np, pl):
    # Panel C — separation quality (ROC AUC of −AB_LR) per coverage, computed
    # rank-based (Mann–Whitney): AUC = P(FP scores more artifact-like than TP).
    def auc_more_negative_is_fp(scores, is_fp):
        # higher "artifact score" = more negative AB_LR -> use -ab_lr
        x = -scores
        order = np.argsort(x, kind="mergesort")
        ranks = np.empty(len(x), float)
        ranks[order] = np.arange(1, len(x) + 1)
        # average ties
        # (simple approach; ties rare in continuous AB_LR)
        n_fp = is_fp.sum()
        n_tp = len(x) - n_fp
        if n_fp == 0 or n_tp == 0:
            return float("nan"), int(n_tp), int(n_fp)
        sum_ranks_fp = ranks[is_fp].sum()
        auc = (sum_ranks_fp - n_fp * (n_fp + 1) / 2) / (n_fp * n_tp)
        return float(auc), int(n_tp), int(n_fp)

    rows_c = []
    for cov_c in sorted(df["coverage"].unique().to_list(), key=lambda c: int(c.rstrip("x"))):
        sc = df.filter(pl.col("coverage") == cov_c)
        a, ntp, nfp = auc_more_negative_is_fp(
            sc["ab_lr"].to_numpy(), (sc["status"].to_numpy() == "FP")
        )
        rows_c.append(dict(coverage=cov_c, cov_x=int(cov_c.rstrip("x")), auc=a, n_tp=ntp, n_fp=nfp))
    auc_df = pl.DataFrame(rows_c).sort("cov_x")
    return (auc_df,)


@app.cell
def _(alt, auc_df, mo):
    chart_c = (
        alt.Chart(auc_df).mark_line(point=True).encode(
            x=alt.X("cov_x:Q", title="coverage (×)", scale=alt.Scale(type="log")),
            y=alt.Y("auc:Q", title="AUC (FP more artifact-like than TP)",
                    scale=alt.Scale(domain=[0.5, 1.0])),
            tooltip=["coverage", alt.Tooltip("auc:Q", format=".3f"), "n_tp", "n_fp"],
        ).properties(width=560, height=300, title="AB_LR separation power vs coverage")
    )
    mo.vstack([
        mo.md("## C. How well does AB_LR separate, by coverage?"),
        mo.ui.altair_chart(chart_c),
        mo.ui.table(auc_df, selection=None),
        mo.md("_AUC≈1 = clean separation. Expect it to fall at low coverage (the "
              "signal genuinely isn't there) and the FP counts there are tiny anyway._"),
    ])
    return


@app.cell
def _(mo):
    thr = mo.ui.slider(-60.0, 0.0, step=1.0, value=-10.0,
                       label="drop if AB_LR <", show_value=True)
    return (thr,)


@app.cell
def _(alt, df, mo, pl, thr):
    # Panel D — apply ONE AB_LR threshold across all coverages; does it hold?
    t = thr.value
    g = df.with_columns((pl.col("ab_lr") < t).alias("drop"))
    rows_d = []
    for cov_d in sorted(g["coverage"].unique().to_list(), key=lambda c: int(c.rstrip("x"))):
        sd = g.filter(pl.col("coverage") == cov_d)
        tp_tot = sd.filter(pl.col("status") == "TP").height
        fp_tot = sd.filter(pl.col("status") == "FP").height
        tp_drop = sd.filter((pl.col("status") == "TP") & pl.col("drop")).height
        fp_drop = sd.filter((pl.col("status") == "FP") & pl.col("drop")).height
        rows_d.append(dict(
            coverage=cov_d, cov_x=int(cov_d.rstrip("x")),
            tp_dropped=tp_drop, tp_total=tp_tot,
            fp_dropped=fp_drop, fp_total=fp_tot,
            tp_loss=round(tp_drop / tp_tot, 4) if tp_tot else 0.0,
            fp_removed=round(fp_drop / fp_tot, 4) if fp_tot else 0.0,
        ))
    res = pl.DataFrame(rows_d).sort("cov_x")
    long = res.unpivot(index="cov_x", on=["tp_loss", "fp_removed"],
                       variable_name="metric", value_name="frac")
    chart_d = (
        alt.Chart(long).mark_line(point=True).encode(
            x=alt.X("cov_x:Q", title="coverage (×)", scale=alt.Scale(type="log")),
            y=alt.Y("frac:Q", title="fraction", scale=alt.Scale(domain=[0, 1])),
            color=alt.Color("metric:N",
                            scale=alt.Scale(domain=["tp_loss", "fp_removed"],
                                            range=["#d62728", "#2ca02c"])),
            tooltip=["cov_x", "metric", alt.Tooltip("frac:Q", format=".3f")],
        ).properties(width=560, height=300,
                     title="One threshold across coverages: FP removed (green) vs TP lost (red)")
    )
    mo.vstack([
        mo.md("## D. Does ONE threshold work across coverages?"),
        thr,
        mo.ui.altair_chart(chart_d),
        mo.ui.table(res, selection=None),
        mo.md("_Production-viable if green (FP removed) stays high where FPs matter "
              "(high cov) while red (TP lost) stays near zero everywhere — especially "
              "at low coverage, where a depth-blind VAF floor would bleed real hets._"),
    ])
    return


@app.cell
def _(mo):
    mo.md(
        """
        ---
        **Reading guide.** `AB_LR` conditions on the called genotype, so hom-alts
        (expected VAF≈1) are judged against ~1.0 and pass; only calls whose alt
        fraction fits *no* genotype are pushed negative. The overdispersion `s` is
        fit once on high-coverage TP hets and reused at all depths — that is what
        makes the low-coverage behaviour conservative rather than trigger-happy.
        Recall here is over *called* TPs (truth-level FN are separate). This is the
        Route-2 (annotation/filter) experiment; if it holds, the same Beta-Binomial
        moves into the genotype likelihood as Route 1.
        """
    )
    return


if __name__ == "__main__":
    app.run()
