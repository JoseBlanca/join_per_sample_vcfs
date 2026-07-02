# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "polars",
#     "numpy",
#     "matplotlib",
# ]
# ///

# Marimo EDA dashboard — paralog vs introgression in the tomato2 cohort.
#
# Reads the per-site feature table built by build_feature_table.py:
#   benchmarks/tomato2/results/snp_features.parquet
# (268,794 biallelic SNPs, 59 samples, WGS sliced to 160x200kb random regions,
#  ~6x median coverage; max-retention callset so paralog FPs are present).
#
# Candidate axes (all "coverage" = read depth, per-sample DP normalised by that
# sample's single-copy coverage scale) and the response (obs het = FP proxy):
#   copy-number  mean_rel_cov    cohort-mean relative coverage (HIGH-AF paralogs)
#                n_excess        # samples with coverage too high vs the cohort
#                                at this site (LOW-AF paralogs the mean misses)
#                max_rel_cov     single most-excess sample (1-sample paralogs)
#   divergence   mqdiff          MQAlt - MQRef (negative => alt reads lower MAPQ)
#   het-excess   fis             1 - Hobs/Hexp (negative => heterozygote excess)
#   co-segregat. cov_het_gap     rel-cov(het samples) - rel-cov(hom samples)
#   response     obs_het         fraction of called samples heterozygous
#
# The paralog flag here is COVERAGE-ONLY (mean OR n_excess) so it is
# introgression-safe by construction; obs_het / fis are kept as held-out
# validation, not part of the flag.
#
# Recommended invocation (no global marimo install needed):
#   uvx marimo edit --sandbox benchmarks/tomato2/src/dashboard.py
#   uvx marimo run  --sandbox benchmarks/tomato2/src/dashboard.py   # read-only

import marimo

__generated_with = "0.23.11"
app = marimo.App(width="medium")


@app.cell
def _():
    from pathlib import Path
    import marimo as mo
    import numpy as np
    import polars as pl
    import matplotlib.pyplot as plt

    return Path, mo, np, pl, plt


@app.cell
def _(Path, pl):
    base = Path(__file__).resolve().parents[1] / "results"
    FEATURES = base / "snp_features.parquet"
    df = pl.read_parquet(FEATURES)
    # join the 500bp windowed, GC-corrected coverage features (option c, the
    # recommended operating point: small window + count threshold) if present
    wf = base / "locus_window_features.w500.parquet"
    if wf.exists():
        wcols = ["chrom", "pos", "n_excess_win", "n_excess_carrier_win",
                 "mean_win_gcrel", "max_win_gcrel", "win_cov_het_gap"]
        df = df.join(pl.read_parquet(wf).select(wcols), on=["chrom", "pos"], how="left")
    # join the paralog likelihood ratio (H2 paralog vs H1 variant) if present
    pr = base / "paralog_lr.parquet"
    if pr.exists():
        prc = ["chrom", "pos", "lr", "n_homalt_conf"]
        prcols = pl.read_parquet(pr).columns
        for extra in ("post", "qval"):   # added by build_paralog_eb.py
            if extra in prcols:
                prc.append(extra)
        df = df.join(pl.read_parquet(pr).select(prc), on=["chrom", "pos"], how="left")
    return FEATURES, df


@app.cell
def _(FEATURES, df, mo):
    mo.md(f"""
    # Paralog vs introgression — tomato2 cohort EDA

    **{df.height:,} biallelic SNPs**, 59 samples, ~6× WGS over 160×200 kb random regions.
    Source: `{FEATURES.name}` (max-retention callset — paralog FPs deliberately kept).

    **Working model:** *paralog* = excess coverage + heterozygote-excess + divergent ALT
    (REMOVE); *introgression* = divergent ALT at **normal** coverage (KEEP). Coverage is the
    introgression-safe discriminator; MQ bias is shared (diagnostic only); `obs_het` is the
    FP proxy.

    Two coverage views: `mean_rel_cov` = the cohort-average coverage level at the locus, and
    `n_excess` = **how many individual samples have coverage too high here vs their own coverage
    distribution** (DP above the sample's upper-1% tail). The mean gives the average level; the
    count says how many samples look anomalously deep. The flag below is **coverage-only**
    (`mean_rel_cov` OR `n_excess`) — introgression-safe; `fis`/`obs_het` are held-out validation.
    """)
    return


@app.cell
def _(mo):
    mo.accordion({
        "📖 Glossary — what each measure means (click to expand)": mo.md(
            r"""
            **Coverage scale (the foundation).** Each *sample's* single-copy coverage = the
            median read depth across its SNP sites. A sample's depth at a locus ÷ this scale =
            **relative coverage**, where **1.0 = single copy**.

            | term | plain meaning |
            |---|---|
            | `mean_rel_cov` | cohort-average relative coverage; >1 ⇒ extra copies on average (**high-frequency** paralogs) |
            | `median_rel_cov` | robust cohort center (not pulled up by a few carriers) |
            | `max_rel_cov` | the single most over-covered sample (**1-sample** paralog) |
            | `n_excess` | **how many samples** have coverage too high here vs **their own** coverage distribution (DP above the sample's upper-`α` empirical tail, `α`=0.01). A per-sample count, complementary to the cohort mean |
            | `n_excess_carrier` | of those over-covered samples, **how many also carry the ALT** (het/hom-alt) = the actual paralog carriers. Sharpest single discriminator (obs_het 0.05→0.69 as it goes 0→8). Uses genotype, so it's a *detection* feature, not held-out |
            | `obs_het` | fraction of called samples that are heterozygous (the **FP proxy**) |
            | `fis` | heterozygote excess vs Hardy–Weinberg, `1 − obs_het/exp_het`; **negative = more hets than expected** (paralog tell) |
            | `mqdiff` | mean MAPQ of ALT-reads − REF-reads; **negative = ALT reads map worse** (divergence axis) |
            | `mqdifft` | same difference as a Welch's *t* (depth-scaled significance) |
            | `cov_het_gap` | rel-cov(het samples) − rel-cov(hom samples); **>0 = the heterozygotes are the over-covered samples** (co-segregation) |
            | `corr_het_cov` | across-sample corr(is-het, rel-cov); same idea as `cov_het_gap` on a −1…+1 scale |
            | `het_alt_balance` | among hets, mean fraction of reads supporting ALT (~0.5 = clean balanced het) |
            | `qual` | caller variant QUAL; paralog miscalls are often *confidently* wrong (high QUAL) |
            | `af` | EM-estimated alternate-allele frequency |
            | `paralog_flag` | **coverage-only** flag: `mean_rel_cov > thr` OR `n_excess ≥ thr` (introgression-safe; het/fis are held-out validation) |

            **Windowed coverage (option c)** — per-sample coverage in 500 bp windows from the `.psp`,
            **GC-corrected per sample** (`gc_rel`, 1.0 = single copy). Windowing gives the per-sample
            power that ~6× per-base lacks, catching moderate-AF 2× paralogs the per-base counts miss.
            Operating point: 500 bp + `n_excess_carrier_win ≥ 2–3` (see the operating-point panel below).

            | term | plain meaning |
            |---|---|
            | `mean_win_gcrel` | cohort-mean window copy number (GC-corrected). 1.0 = single copy |
            | `n_excess_carrier_win` | # ALT-carrier samples whose 500 bp window coverage > 1.5× — the **best windowed discriminator**; use with a count threshold (≥2–3) to clear the noise floor |
            | `win_cov_het_gap` | windowed het−hom coverage gap. Weak here: partial-dup dilution shrinks the difference, so the count beats it |

            **Paralog likelihood ratio** — the unified test (see panel below).

            | term | plain meaning |
            |---|---|
            | `lr` | logL(hidden-paralog) − logL(real-variant), per locus. **>0 → paralog (remove), <0 → real variant (keep).** Per-sample coverage + allele balance via the copy-number model (carrier total copies `T`, mutant copies `m` → coverage `T/2`, VAF `m/T`, coverage up to 4×); per-sample inbreeding; hom-alt veto; **coverage-gated mqdiff** term (divergence adds paralog evidence only where there's also coverage excess → recovers diverging paralogs, introgression-safe per simulation) |
            | `n_homalt_conf` | # samples confidently hom-alt (VAF>0.9, ≥5 reads) — evidence *for* a real variant (impossible under the paralog model) |
            | `post` | empirical-Bayes posterior **P(paralog \| data) = σ(LR + log(π/(1−π)))**, π estimated from the data (≈0.11). LR is the **marginal** (Laplace) ratio — integrated over the model parameters, so flexibility is penalised |
            | `qval` | q-value (tail FDR): expected fraction of *real variants* among loci flagged at this cutoff — set the target in the panel below |
            """
        )
    })
    return


@app.cell
def _(mo):
    # ---- interactive controls ----
    x_axis = mo.ui.dropdown(
        ["mean_rel_cov", "median_rel_cov", "max_rel_cov", "n_excess", "n_excess_carrier",
         "mean_win_gcrel", "n_excess_carrier_win", "win_cov_het_gap"],
        value="mean_rel_cov", label="x-axis",
    )
    y_axis = mo.ui.dropdown(
        ["obs_het", "lr", "fis", "n_excess", "n_excess_carrier", "n_excess_carrier_win",
         "mean_win_gcrel", "win_cov_het_gap", "cov_het_gap", "corr_het_cov",
         "het_alt_balance", "qual"],
        value="obs_het", label="y-axis",
    )
    color_by = mo.ui.dropdown(
        ["mqdiff", "lr", "post", "qual", "fis", "obs_het", "n_excess", "n_excess_carrier",
         "n_excess_carrier_win", "mean_win_gcrel", "cov_het_gap", "af"],
        value="lr", label="color by",
    )
    xmax = mo.ui.slider(2.0, 10.0, value=4.0, step=0.5, label="x clip (coverage axes)")
    min_called = mo.ui.slider(20, 59, value=30, step=1, label="min called samples")
    cov_thr = mo.ui.slider(1.0, 3.0, value=1.5, step=0.05, label="flag: mean_rel_cov >")
    nexc_thr = mo.ui.slider(1, 20, value=3, step=1, label="flag: OR n_excess >=")
    point_alpha = mo.ui.slider(0.02, 0.5, value=0.12, step=0.02, label="point alpha")
    controls = mo.vstack([
        mo.hstack([x_axis, y_axis, color_by]),
        mo.hstack([xmax, min_called, point_alpha]),
        mo.hstack([cov_thr, nexc_thr]),
    ])
    controls
    return color_by, cov_thr, min_called, nexc_thr, point_alpha, x_axis, xmax, y_axis


@app.cell
def _(cov_thr, df, min_called, nexc_thr, pl):
    # filtered working set + COVERAGE-ONLY paralog flag (introgression-safe)
    work = df.filter(pl.col("n_called") >= min_called.value).with_columns(
        ((pl.col("mean_rel_cov") > cov_thr.value) | (pl.col("n_excess") >= nexc_thr.value))
        .alias("paralog_flag")
    )
    n_flag = int(work["paralog_flag"].sum())
    return n_flag, work


@app.cell
def _(color_by, n_flag, np, plt, point_alpha, work, x_axis, xmax, y_axis):
    def _plot():
        w = work
        cov_axes = {"mean_rel_cov", "median_rel_cov", "max_rel_cov"}
        x = w[x_axis.value].to_numpy().astype(float)
        if x_axis.value in cov_axes:
            x = x.clip(0, xmax.value)
        y = w[y_axis.value].to_numpy()
        c = w[color_by.value].to_numpy()
        flag = w["paralog_flag"].to_numpy()

        if color_by.value in ("mqdiff", "fis", "cov_het_gap", "lr"):
            lim = np.nanpercentile(np.abs(c), 95) or 1.0
            vmin, vmax, cmap = -lim, lim, "coolwarm"
        else:
            # count-like fields (n_excess*) are mostly 0 -> use the upper tail
            vmin, vmax = np.nanpercentile(c, [2, 99.5])
            cmap = "viridis"
        if not (vmax > vmin):  # degenerate (e.g. a mostly-zero count) -> widen
            vmin, vmax = float(np.nanmin(c)), float(np.nanmax(c))
            if not (vmax > vmin):
                vmax = vmin + 1

        fig, ax = plt.subplots(figsize=(9, 6))
        sc = ax.scatter(x, y, c=c, s=4, alpha=point_alpha.value, cmap=cmap,
                        vmin=vmin, vmax=vmax, rasterized=True, linewidths=0)
        ax.scatter(x[flag], y[flag], s=14, facecolors="none", edgecolors="k",
                   linewidths=0.4, alpha=0.5, label=f"paralog-flag (n={n_flag})")
        if x_axis.value in cov_axes:
            ax.axvline(1.0, color="grey", ls=":", lw=1)
        ax.set_xlabel(x_axis.value + ("  (1.0 = single-copy)" if x_axis.value in cov_axes else ""))
        ax.set_ylabel(y_axis.value)
        ax.set_title(f"{y_axis.value} vs {x_axis.value}, coloured by {color_by.value}")
        fig.colorbar(sc, ax=ax, label=color_by.value)
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        return fig

    _plot()
    return


@app.cell
def _(mo):
    mo.md("""
    ### Coverage distribution vs heterozygosity (x-aligned)
    **A** (top): per-locus coverage distribution. **B** (bottom): obs_het boxplots for the
    loci in each coverage bin. Het stays flat through single-copy, then the upper quartile
    climbs sharply above ~1.8× — the collapsed-paralog signature. The median barely moves
    (paralogs are a minority even at high coverage); switch the trend line to **mean** or
    **frac het>0.3** to see the tail signal far more starkly.
    """)
    return


@app.cell
def _(mo):
    cov_fig_col = mo.ui.dropdown(
        ["mean_rel_cov", "median_rel_cov", "max_rel_cov"],
        value="mean_rel_cov", label="coverage axis",
    )
    cov_fig_step = mo.ui.slider(0.05, 0.3, value=0.15, step=0.05, label="bin width")
    cov_fig_overlay = mo.ui.dropdown(
        ["median", "mean", "frac het>0.3"], value="median", label="trend line",
    )
    cov_fig_controls = mo.hstack([cov_fig_col, cov_fig_step, cov_fig_overlay])
    cov_fig_controls
    return cov_fig_col, cov_fig_overlay, cov_fig_step


@app.cell
def _(cov_fig_col, cov_fig_overlay, cov_fig_step, np, plt, work):
    def _plot():
        XLO, XHI = 0.3, 3.0
        step = cov_fig_step.value
        col = cov_fig_col.value
        x = work[col].to_numpy().astype(float)
        het = work["obs_het"].to_numpy()
        edges = np.arange(XLO, XHI + step / 2, step)
        centers = (edges[:-1] + edges[1:]) / 2
        xc = np.clip(x, XLO + 1e-9, XHI - 1e-9)   # fold the long tail into the end bins
        binidx = np.digitize(xc, edges) - 1
        n_over = int((x > XHI).sum())
        het_by_bin = [het[(binidx == i) & np.isfinite(het)] for i in range(len(centers))]

        fig, (axA, axB) = plt.subplots(
            2, 1, figsize=(11, 7), sharex=True,
            gridspec_kw={"height_ratios": [1, 1.4], "hspace": 0.07},
        )
        axA.hist(xc, bins=edges, color="steelblue", edgecolor="white", linewidth=0.4)
        axA.axvline(1.0, color="grey", ls=":", lw=1.2)
        axA.set_yscale("log")
        axA.set_ylabel("loci per bin (log)")
        axA.set_title(f"{work.height:,} SNPs; top bin folds in {n_over:,} loci with coverage >{XHI:g}×")

        bp = axB.boxplot(
            het_by_bin, positions=centers, widths=step * 0.8, showfliers=False,
            patch_artist=True, manage_ticks=False, medianprops=dict(color="black"),
        )
        for patch in bp["boxes"]:
            patch.set(facecolor="lightsteelblue", alpha=0.9)
        if cov_fig_overlay.value == "mean":
            trend = [np.mean(h) if len(h) else np.nan for h in het_by_bin]
            lbl = "mean obs_het"
        elif cov_fig_overlay.value == "frac het>0.3":
            trend = [float(np.mean(h > 0.3)) if len(h) else np.nan for h in het_by_bin]
            lbl = "fraction obs_het > 0.3"
        else:
            trend = [np.median(h) if len(h) else np.nan for h in het_by_bin]
            lbl = "median obs_het"
        axB.plot(centers, trend, color="crimson", lw=1.5, marker="o", ms=3, label=lbl)
        axB.axvline(1.0, color="grey", ls=":", lw=1.2)
        axB.set_ylim(0, 1)
        axB.set_ylabel("observed heterozygosity")
        axB.set_xlabel(f"{col}  (1.0 = single copy)")
        axB.legend(loc="upper left", fontsize=9)
        axB.set_xlim(XLO, XHI)
        return fig

    _plot()
    return


@app.cell
def _(mo):
    mo.md("""
    ### Measures vs heterozygosity (x-aligned)
    **A** (top): the obs_het distribution (how many loci at each het level). **B, C, …**:
    each selected measure boxplotted per het bin. Read the trends: mean coverage, `n_excess`,
    `n_excess_carrier` (the sharpest), and co-segregation all climb with het; divergence is
    U-shaped (most negative at *mid* het — the normal-coverage introgression-like loci, not
    the high-het paralogs).
    """)
    return


@app.cell
def _(mo):
    het_step = mo.ui.slider(0.025, 0.1, value=0.05, step=0.025, label="het bin width")
    het_measures = mo.ui.multiselect(
        ["post", "lr", "mean_rel_cov", "median_rel_cov", "max_rel_cov", "n_excess", "n_excess_carrier",
         "mean_win_gcrel", "n_excess_carrier_win", "win_cov_het_gap", "qval",
         "mqdiff", "mqdifft", "cov_het_gap", "corr_het_cov", "het_alt_balance", "qual"],
        value=["post", "lr", "n_excess_carrier", "n_excess_carrier_win", "mqdiff"],
        label="measures (panels B, C, …)",
    )
    mo.hstack([het_step, het_measures])
    return het_measures, het_step


@app.cell
def _(het_measures, het_step, np, plt, work):
    def _plot():
        ylims = {
            "mean_rel_cov": (0, 3), "median_rel_cov": (0, 3), "max_rel_cov": (0, 6),
            "n_excess": (0, 12), "n_excess_carrier": (0, 8),
            "mean_win_gcrel": (0, 3), "n_excess_carrier_win": (0, 8),
            "win_cov_het_gap": (-1, 2), "max_win_gcrel": (0, 4), "lr": (-30, 40),
            "post": (0, 1), "qval": (0, 1),
            "mqdiff": (-25, 10), "mqdifft": (-20, 10),
            "cov_het_gap": (-2, 3), "corr_het_cov": (-1, 1),
            "het_alt_balance": (0, 1), "qual": (0, 5000),
        }
        step = het_step.value
        sel = list(het_measures.value)
        het = work["obs_het"].to_numpy()
        edges = np.arange(0.0, 1.0 + step / 2, step)
        centers = (edges[:-1] + edges[1:]) / 2
        binidx = np.digitize(np.clip(het, 0, 1 - 1e-9), edges) - 1

        if not sel:
            fig, ax = plt.subplots(figsize=(11, 1.5))
            ax.text(0.5, 0.5, "select at least one measure", ha="center", va="center")
            ax.axis("off")
            return fig

        nrows = 1 + len(sel)
        fig, axes = plt.subplots(
            nrows, 1, figsize=(11, 2.3 * nrows), sharex=True,
            gridspec_kw={"height_ratios": [1] + [1.2] * len(sel), "hspace": 0.12},
        )
        axes[0].hist(het, bins=edges, color="seagreen", edgecolor="white", linewidth=0.4)
        axes[0].set_yscale("log")
        axes[0].set_ylabel("loci/bin\n(log)")
        axes[0].set_title(f"{work.height:,} SNPs — measures vs observed heterozygosity")

        for ax, col in zip(axes[1:], sel):
            vals = work[col].to_numpy().astype(float)
            by_bin = [vals[(binidx == i) & np.isfinite(vals)] for i in range(len(centers))]
            bp = ax.boxplot(
                by_bin, positions=centers, widths=step * 0.8, showfliers=False,
                patch_artist=True, manage_ticks=False, medianprops=dict(color="black"),
            )
            for patch in bp["boxes"]:
                patch.set(facecolor="lightsteelblue", alpha=0.9)
            med = [np.median(b) if len(b) else np.nan for b in by_bin]
            ax.plot(centers, med, color="crimson", lw=1.3, marker="o", ms=2.5)
            lo, hi = ylims.get(col) or tuple(np.nanpercentile(vals, [1, 99]))
            ax.set_ylim(lo, hi)
            ax.set_ylabel(col, fontsize=8)
            if lo < 0 < hi:
                ax.axhline(0, color="grey", ls=":", lw=0.8)

        axes[-1].set_xlabel("observed heterozygosity")
        axes[-1].set_xlim(0, 1)
        return fig

    _plot()
    return


@app.cell
def _(mo):
    mo.md("""
    ### Marginal distributions
    Where does each axis sit, and how does the flagged class stand apart (red)?
    """)
    return


@app.cell
def _(np, plt, work):
    def _plot():
        cols = ["mean_rel_cov", "n_excess", "obs_het", "fis", "mqdiff"]
        clips = {"mean_rel_cov": (0, 5), "n_excess": (0, 20), "obs_het": (0, 1),
                 "fis": (-1, 1), "mqdiff": (-25, 25)}
        fig, axes = plt.subplots(1, 5, figsize=(18, 3.0))
        for ax, col in zip(axes, cols):
            lo, hi = clips[col]
            allv = work[col].to_numpy().astype(float)
            allv = allv[np.isfinite(allv)].clip(lo, hi)
            flv = work.filter(work["paralog_flag"])[col].to_numpy().astype(float)
            flv = flv[np.isfinite(flv)].clip(lo, hi)
            bins = np.linspace(lo, hi, 50)
            ax.hist(allv, bins=bins, color="steelblue", alpha=0.6, density=True, label="all")
            if flv.size:
                ax.hist(flv, bins=bins, color="crimson", alpha=0.6, density=True, label="flagged")
            ax.set_title(col)
            ax.set_yticks([])
        axes[0].legend(fontsize=8)
        fig.tight_layout()
        return fig

    _plot()
    return


@app.cell
def _(mo, n_flag, pl, work):
    def _summarise(d):
        return d.select([
            pl.len().alias("n"),
            pl.col("obs_het").drop_nans().mean().alias("obs_het"),
            pl.col("mean_rel_cov").drop_nans().mean().alias("mean_cov"),
            pl.col("n_excess").mean().alias("n_excess"),
            pl.col("fis").drop_nans().mean().alias("fis"),
            pl.col("mqdiff").drop_nans().mean().alias("mqdiff"),
            pl.col("cov_het_gap").drop_nans().mean().alias("cov_het_gap"),
            pl.col("qual").drop_nans().median().alias("med_qual"),
        ])

    _flagged = _summarise(work.filter(pl.col("paralog_flag"))).with_columns(pl.lit("paralog-flag").alias("class"))
    _rest = _summarise(work.filter(~pl.col("paralog_flag"))).with_columns(pl.lit("rest").alias("class"))
    _tbl = pl.concat([_flagged, _rest]).select(
        ["class", "n", "obs_het", "mean_cov", "n_excess", "fis", "mqdiff", "cov_het_gap", "med_qual"]
    )
    mo.vstack([
        mo.md(f"**Flagged paralog-like sites: {n_flag:,}** — flag is coverage-only; "
              "obs_het / fis here are *held-out validation*."),
        mo.ui.table(_tbl, selection=None),
    ])
    return


@app.cell
def _(mo):
    mo.md("""
    ### `n_excess` vs mean coverage — the per-sample count adds a dimension
    `n_excess` (y) vs `mean_rel_cov` (x). The two are related, but circled points sit **left**
    of mean=1.5 yet **above** the dashed line: the cohort-mean coverage looks ~normal, yet
    several individual samples are in their own upper tail — candidate low-frequency duplications
    a cohort average would miss. Coloured by obs_het (held-out).
    """)
    return


@app.cell
def _(nexc_thr, np, plt, work):
    def _plot():
        w = work
        x = w["mean_rel_cov"].to_numpy().clip(0, 4)
        y = w["n_excess"].to_numpy().astype(float)
        c = w["obs_het"].to_numpy().clip(0, 1)
        # mean-view misses: normal mean coverage but n_excess over threshold
        miss = (w["mean_rel_cov"].to_numpy() <= 1.5) & (w["n_excess"].to_numpy() >= nexc_thr.value)

        fig, ax = plt.subplots(figsize=(9, 5.5))
        sc = ax.scatter(x, y, c=c, s=5, alpha=0.2, cmap="viridis", vmin=0, vmax=1,
                        rasterized=True, linewidths=0)
        ax.scatter(x[miss], y[miss], s=16, facecolors="none", edgecolors="crimson",
                   linewidths=0.5, alpha=0.6, label=f"low-AF (mean<=1.5, n_excess>={nexc_thr.value}): {int(miss.sum())}")
        ax.axhline(nexc_thr.value, color="k", ls="--", lw=0.8)
        ax.axvline(1.5, color="grey", ls=":", lw=1)
        ax.set_xlabel("mean_rel_cov  (1.0 = single-copy)")
        ax.set_ylabel("n_excess (samples over-covered vs cohort)")
        ax.set_title("low-AF paralog rescue")
        fig.colorbar(sc, ax=ax, label="obs_het")
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        return fig

    _plot()
    return


@app.cell
def _(mo):
    mo.md("""
    ### Genome view — do paralog sites cluster into regions?
    Relative coverage (top) and obs_het (bottom) along a chromosome; flagged sites in red.
    Tight clusters of red = candidate collapsed-paralog regions (objective 1).
    """)
    return


@app.cell
def _(df, mo):
    chrom = mo.ui.dropdown(
        sorted(df["chrom"].unique().to_list()), value="SL4.0ch01", label="chromosome"
    )
    chrom
    return (chrom,)


@app.cell
def _(chrom, plt, work):
    def _plot():
        cw = work.filter(work["chrom"] == chrom.value).sort("pos")
        pos = cw["pos"].to_numpy() / 1e6
        cov = cw["mean_rel_cov"].to_numpy().clip(0, 6)
        het = cw["obs_het"].to_numpy()
        flag = cw["paralog_flag"].to_numpy()

        fig, (a1, a2) = plt.subplots(2, 1, figsize=(14, 5), sharex=True)
        a1.scatter(pos, cov, s=3, alpha=0.25, color="grey", rasterized=True)
        a1.scatter(pos[flag], cov[flag], s=12, color="crimson", label="flagged")
        a1.axhline(1.0, color="k", ls=":", lw=0.8)
        a1.set_ylabel("rel coverage")
        a1.legend(fontsize=8, loc="upper right")
        a1.set_title(f"{chrom.value}: {cw.height:,} SNP sites, {int(flag.sum())} flagged")
        a2.scatter(pos, het, s=3, alpha=0.25, color="grey", rasterized=True)
        a2.scatter(pos[flag], het[flag], s=12, color="crimson")
        a2.set_ylabel("obs het")
        a2.set_xlabel("position (Mb)")
        fig.tight_layout()
        return fig

    _plot()
    return


@app.cell
def _(mo):
    mo.md("""
    ### High-het universe split by coverage
    Among `obs_het > 0.3` sites (where paralog FPs and legitimate high-het loci both
    live), does coverage split a paralog mode from a normal-coverage keep mode? Note mqdiff is
    *more* negative for the normal-coverage keep class — divergence ≠ paralog signal.
    """)
    return


@app.cell
def _(mo, pl, work):
    _hh = work.filter(pl.col("obs_het") > 0.3)
    _classed = _hh.with_columns(
        pl.when(pl.col("mean_rel_cov") > 1.5).then(pl.lit("excess>1.5 (paralog)"))
        .when(pl.col("mean_rel_cov") > 1.2).then(pl.lit("mild 1.2-1.5 (ambiguous)"))
        .otherwise(pl.lit("normal<=1.2 (keep)")).alias("cov_class")
    )
    _split = (_classed.group_by("cov_class").agg(
        pl.len().alias("n"),
        pl.col("mean_rel_cov").median().alias("mean_cov"),
        pl.col("n_excess").median().alias("n_excess"),
        pl.col("fis").drop_nans().median().alias("fis"),
        pl.col("mqdiff").drop_nans().median().alias("mqdiff"),
        pl.col("cov_het_gap").drop_nans().median().alias("cov_het_gap"),
        pl.col("qual").median().alias("med_qual"),
    ).sort("mean_cov"))
    mo.vstack([
        mo.md(f"**high-het sites: {_hh.height:,}**"),
        mo.ui.table(_split, selection=None),
    ])
    return


@app.cell
def _(mo):
    mo.md("""
    ### GC content vs window coverage
    Per-window coverage boxplotted against the window's GC content (window-level, option c —
    not per-site). **`rel`** = raw per-sample relative coverage → shows the GC-bias hump
    (low-GC windows under-covered); **`gc_rel`** = after per-sample GC correction → should be
    flat at 1.0, confirming the correction worked.
    """)
    return


@app.cell
def _(Path, pl):
    _b = Path(__file__).resolve().parents[1] / "results"
    wtabs = {}
    for _w in (500, 1000, 2000):
        _p = _b / f"window_cov.w{_w}.gcnorm.parquet"
        if _p.exists():
            wtabs[_w] = pl.read_parquet(_p).select(["gc", "rel", "gc_rel", "breadth"])
    return (wtabs,)


@app.cell
def _(mo, wtabs):
    gc_w = mo.ui.dropdown([str(k) for k in sorted(wtabs)],
                          value=str(max(wtabs)) if wtabs else "2000", label="window size (bp)")
    gc_y = mo.ui.dropdown(["rel", "gc_rel"], value="rel",
                          label="y (rel=raw, gc_rel=GC-corrected)")
    gc_binw = mo.ui.slider(0.01, 0.05, value=0.02, step=0.01, label="GC bin width")
    mo.hstack([gc_w, gc_y, gc_binw])
    return gc_binw, gc_w, gc_y


@app.cell
def _(gc_binw, gc_w, gc_y, np, pl, plt, wtabs):
    def _plot():
        w = int(gc_w.value)
        d = wtabs[w].filter(pl.col("breadth") > 0.8)
        gc = d["gc"].to_numpy()
        y = d[gc_y.value].to_numpy()
        step = gc_binw.value
        edges = np.arange(0.15, 0.55 + step / 2, step)
        centers = (edges[:-1] + edges[1:]) / 2
        bi = np.digitize(gc, edges) - 1
        by_bin = [np.clip(y[(bi == i) & np.isfinite(y)], 0, 3) for i in range(len(centers))]
        pos = [c for c, b in zip(centers, by_bin) if len(b) >= 20]
        data = [b for b in by_bin if len(b) >= 20]

        fig, ax = plt.subplots(figsize=(11, 5))
        bp = ax.boxplot(data, positions=pos, widths=step * 0.8, showfliers=False,
                        patch_artist=True, manage_ticks=False, medianprops=dict(color="black"))
        for patch in bp["boxes"]:
            patch.set(facecolor="lightsteelblue", alpha=0.9)
        ax.plot(pos, [np.median(b) for b in data], color="crimson", lw=1.5, marker="o", ms=3, label="median")
        ax.axhline(1.0, color="grey", ls=":", lw=1)
        ax.set_xlabel(f"GC content of {w} bp window")
        ax.set_ylabel(gc_y.value + (" — raw relative coverage" if gc_y.value == "rel"
                                    else " — GC-corrected"))
        ax.set_title(f"window coverage vs GC ({w} bp; n={d.height:,} window×sample)")
        ax.set_xlim(0.15, 0.52)
        ax.set_ylim(0, 2.2)
        ax.legend(fontsize=8)
        return fig

    _plot()
    return


@app.cell
def _(mo):
    mo.md("""
    ### Paralog likelihood ratio + empirical-Bayes calibration — the unified test
    `lr` = logL(hidden-paralog H2) − logL(real-variant H1), per locus (per-sample 500 bp coverage
    **+** allele balance; **copy-number model**: carrier total copies `T`, mutant copies `m` →
    coverage `T/2`, VAF `m/T`, so high copy ⇒ low VAF; **per-sample inbreeding**; hom-alt veto). Empirical Bayes
    turns it into a posterior **P(paralog) = σ(LR + log(π/(1−π)))**, with the paralog fraction π
    estimated from the data — so the cut is **FDR-controlled**, landing at a *positive* LR, not the
    naive 0. Set a target FDR; the table's obs_het / fis / balance / homalt are held-out validation.
    """)
    return


@app.cell
def _(df, mo):
    _has = "qval" in df.columns
    fdr_thr = mo.ui.slider(0.01, 0.5, value=0.05, step=0.01, label="target FDR (flag if qval ≤)")
    (fdr_thr if _has else mo.md("*run build_paralog_lr.py then build_paralog_eb.py*"))
    return (fdr_thr,)


@app.cell
def _(df, fdr_thr, np, pl, plt):
    def _plot():
        if "qval" not in df.columns:
            fig, ax = plt.subplots(figsize=(8, 1)); ax.axis("off"); return fig
        d = df.filter(pl.col("lr").is_not_null())
        lr = d["lr"].to_numpy()
        flagged = d.filter(pl.col("qval") <= fdr_thr.value)
        lrcut = float(flagged["lr"].min()) if flagged.height else float("nan")

        fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 4.5))
        a1.hist(np.clip(lr, -20, 60), bins=180, color="steelblue")
        a1.axvline(0, color="grey", ls=":", label="naive LR=0")
        if np.isfinite(lrcut):
            a1.axvline(min(lrcut, 60), color="crimson", lw=1.5,
                       label=f"FDR {fdr_thr.value:.0%} cut ≈ LR {lrcut:.1f}")
        a1.set_yscale("log")
        a1.set_xlabel("LR = logL(paralog) − logL(variant)")
        a1.set_ylabel("loci (log)")
        a1.legend(fontsize=8)
        a1.set_title(f"LR distribution — {flagged.height:,} flagged at FDR ≤ {fdr_thr.value:.0%}")

        x = np.clip(lr, -20, 60)
        y = d["obs_het"].to_numpy()
        c = np.clip(d["mean_rel_cov"].to_numpy(), 0, 3)
        sc = a2.scatter(x, y, c=c, s=3, alpha=0.12, cmap="viridis", rasterized=True, linewidths=0)
        if np.isfinite(lrcut):
            a2.axvline(min(lrcut, 60), color="crimson", lw=1.5)
        a2.set_xlabel("LR")
        a2.set_ylabel("obs_het")
        a2.set_title("LR vs obs_het (colour = mean_rel_cov)")
        fig.colorbar(sc, ax=a2, label="mean_rel_cov")
        fig.tight_layout()
        return fig

    _plot()
    return


@app.cell
def _(df, fdr_thr, mo, pl):
    def _table():
        if "qval" not in df.columns:
            return mo.md("")
        d = df.filter(pl.col("lr").is_not_null())

        def summ(x, lab):
            return x.select(
                pl.len().alias("n"),
                pl.col("obs_het").mean().round(3).alias("obs_het"),
                pl.col("mean_rel_cov").drop_nans().median().round(2).alias("mean_cov"),
                pl.col("fis").drop_nans().mean().round(2).alias("fis"),
                pl.col("het_alt_balance").drop_nans().median().round(2).alias("balance"),
                pl.col("qual").median().round(0).alias("med_qual"),
                pl.col("n_homalt_conf").mean().round(2).alias("homalt"),
            ).with_columns(pl.lit(lab).alias("class"))

        flagged = summ(d.filter(pl.col("qval") <= fdr_thr.value), f"flagged (FDR≤{fdr_thr.value:.0%})")
        kept = summ(d.filter(pl.col("qval") > fdr_thr.value), "kept")
        tbl = pl.concat([flagged, kept]).select(
            ["class", "n", "obs_het", "mean_cov", "fis", "balance", "med_qual", "homalt"])
        nflag = d.filter(pl.col("qval") <= fdr_thr.value).height
        return mo.vstack([
            mo.md(f"**flagged at FDR ≤ {fdr_thr.value:.0%}: {nflag:,} loci** "
                  f"({100*nflag/d.height:.1f}%) — obs_het / fis / balance / homalt are held-out validation"),
            mo.ui.table(tbl, selection=None),
        ])

    _table()
    return


@app.cell
def _(mo):
    mo.md("""
    ### Operating point — window size × carrier-count threshold
    The windowed carrier count trades recall vs precision. **Capture** = fraction of moderate-het
    (obs_het 0.2–0.5) loci flagged (the paralog class per-base missed); **false-flag** = fraction of
    clean low-het (obs_het < 0.05) loci flagged (noise proxy). Smaller windows raise capture but also
    noise; the **count threshold** is the stronger noise filter, so small window + count ≥2–3 wins.
    Pick a window size and read the curves; the count threshold is the x-axis.
    """)
    return


@app.cell
def _(Path, pl):
    _b = Path(__file__).resolve().parents[1] / "results"
    _feat = pl.read_parquet(_b / "snp_features.parquet").select(["chrom", "pos", "obs_het"])
    opdata = {}
    for _w in (500, 1000, 2000):
        _p = _b / f"locus_window_features.w{_w}.parquet"
        if _p.exists():
            opdata[_w] = _feat.join(
                pl.read_parquet(_p).select(["chrom", "pos", "n_excess_carrier_win"]),
                on=["chrom", "pos"], how="inner")
    return (opdata,)


@app.cell
def _(mo, opdata):
    op_w = mo.ui.dropdown([str(k) for k in sorted(opdata)],
                          value="500" if 500 in opdata else str(min(opdata)),
                          label="window size (bp)")
    op_w
    return (op_w,)


@app.cell
def _(np, op_w, opdata, pl, plt):
    def _plot():
        d = opdata[int(op_w.value)]
        mod = d.filter((pl.col("obs_het") >= 0.2) & (pl.col("obs_het") <= 0.5))
        low = d.filter(pl.col("obs_het") < 0.05)
        thrs = list(range(1, 7))
        cap = [mod.filter(pl.col("n_excess_carrier_win") >= t).height / mod.height for t in thrs]
        fp = [low.filter(pl.col("n_excess_carrier_win") >= t).height / low.height for t in thrs]
        enr = [c / f if f else np.nan for c, f in zip(cap, fp)]

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(thrs, [100 * c for c in cap], "o-", color="steelblue", label="moderate-het capture %")
        ax.plot(thrs, [100 * f for f in fp], "o-", color="crimson", label="low-het false-flag %")
        ax.set_xlabel("flag if  n_excess_carrier_win ≥  (threshold)")
        ax.set_ylabel("% of loci")
        ax.legend(loc="center right", fontsize=8)
        ax2 = ax.twinx()
        ax2.plot(thrs, enr, "s--", color="grey", alpha=0.7, label="enrichment ×")
        ax2.set_ylabel("enrichment (capture / false-flag)  ×")
        ax2.legend(loc="upper right", fontsize=8)
        ax.set_title(f"operating point — {op_w.value} bp windows  "
                     f"(capture={100*cap[1]:.0f}% @≥2, {100*cap[2]:.0f}% @≥3)")
        fig.tight_layout()
        return fig

    _plot()
    return


if __name__ == "__main__":
    app.run()
