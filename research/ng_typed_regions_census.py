# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo", "polars", "altair"]
# ///

import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import altair as alt
    from pathlib import Path

    # The walk's output, which is bulky (tomato alone is ~260 MB) and so lives in
    # the gitignored `tmp/` rather than beside this notebook. Resolved relative to
    # the repo so the notebook is not pinned to one machine.
    REPO = Path(__file__).resolve().parents[1]
    DATA = REPO / "tmp" / "ng_typed_regions"
    return DATA, Path, alt, mo, pl


@app.cell
def _(Path, pl):
    def read_run_params(path: Path) -> dict[str, str]:
        """The `## key: value` header block the walk writes, so a file says what made it."""
        params = {}
        with path.open() as fh:
            for line in fh:
                if not line.startswith("##"):
                    break
                key, _, value = line[2:].strip().partition(":")
                params[key.strip()] = value.strip()
        return params

    def read_regions(path: Path, genome: str) -> pl.DataFrame:
        """One typed-regions TSV -> a frame with a `length` column.

        `comment_prefix='##'` drops the run-parameter block but keeps the `#chrom`
        header row, which is where the column names live.
        """
        return (
            pl.read_csv(
                path,
                separator="\t",
                comment_prefix="##",
                null_values=".",
                schema_overrides={"motif": pl.String},
            )
            .rename({"#chrom": "chrom"})
            .with_columns(
                (pl.col("end") - pl.col("start")).alias("length"),
                pl.lit(genome).alias("genome"),
            )
        )

    return read_regions, read_run_params


@app.cell
def _(DATA, read_regions, read_run_params):
    # Human is added here once its walk finishes; tomato alone for now.
    SOURCES = {
        "tomato SL4.00": DATA / "tomato_SL4.00.tsv",
        "human GRCh38": DATA / "human_GRCh38.tsv",
    }
    available = {g: p for g, p in SOURCES.items() if p.exists()}

    frames = {g: read_regions(p, g) for g, p in available.items()}
    params = {g: read_run_params(p) for g, p in available.items()}
    return frames, params


@app.cell
def _(frames, mo, pl):
    regions = pl.concat(list(frames.values()))
    mo.md(f"**{len(frames)} genome(s) loaded — {regions.height:,} typed regions total**")
    return (regions,)


@app.cell
def _(params):
    # The run parameters, so the numbers below are reproducible and comparable.
    params
    return


@app.cell
def _(pl, regions):
    # Q1: how many segments of each type, and how much of the genome do they cover?
    census = (
        regions.group_by("genome", "kind")
        .agg(
            pl.len().alias("count"),
            pl.col("length").sum().alias("total_bp"),
            pl.col("length").mean().round(1).alias("mean_bp"),
            pl.col("length").median().alias("median_bp"),
            pl.col("length").max().alias("max_bp"),
        )
        .with_columns(
            (100 * pl.col("total_bp") / pl.col("total_bp").sum().over("genome"))
            .round(2)
            .alias("pct_of_genome"),
            (100 * pl.col("count") / pl.col("count").sum().over("genome"))
            .round(2)
            .alias("pct_of_segments"),
        )
        .sort("genome", "total_bp", descending=[False, True])
    )
    census
    return


@app.cell
def _(pl, regions):
    # Q2: the length distribution of each type, as quantiles.
    quantiles = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
    length_quantiles = (
        regions.group_by("genome", "kind")
        .agg(
            [pl.col("length").min().alias("min")]
            + [
                pl.col("length").quantile(q).alias(f"p{int(q * 100)}")
                for q in quantiles
            ]
            + [pl.col("length").max().alias("max")]
        )
        .sort("genome", "kind")
    )
    length_quantiles
    return


@app.cell
def _(alt, pl):
    # Log-spaced histogram, computed in Polars: 4.5M rows is far past Altair's
    # row limit, so only the bin table (a few hundred rows) ever reaches the chart.
    import math

    KIND_HUE = {
        "generic": "#2a78d6",     # categorical slot 1
        "ssr_locus": "#008300",   # slot 2
        "ssr_bundle": "#e87ba4",  # slot 3
        "satellite": "#eda100",   # slot 4
    }


    def log_histogram(df, kind, genome, bins_per_decade=12):
        """Count region lengths into bins evenly spaced in log10(bp)."""
        lengths = df.filter(
            (pl.col("kind") == kind) & (pl.col("genome") == genome)
        ).select("length")
        lo = max(1, lengths["length"].min())
        hi = lengths["length"].max()
        n_bins = max(1, int(math.ceil((math.log10(hi) - math.log10(lo)) * bins_per_decade)))
        edges = [10 ** (math.log10(lo) + i / bins_per_decade) for i in range(n_bins + 1)]
        return (
            lengths.with_columns(
                pl.col("length")
                .log10()
                .sub(math.log10(lo))
                .mul(bins_per_decade)
                .floor()
                .clip(0, n_bins - 1)
                .cast(pl.Int32)
                .alias("bin")
            )
            .group_by("bin")
            .agg(pl.len().alias("count"))
            .with_columns(
                pl.col("bin").map_elements(lambda b: edges[b], return_dtype=pl.Float64).alias("bin_lo"),
                pl.col("bin").map_elements(lambda b: edges[b + 1], return_dtype=pl.Float64).alias("bin_hi"),
            )
            .sort("bin")
        )


    def length_hist_chart(df, kind, genome, title):
        """One kind's length distribution. Single series, so no legend: the title names it."""
        hist = log_histogram(df, kind, genome)
        median = df.filter(
            (pl.col("kind") == kind) & (pl.col("genome") == genome)
        )["length"].median()

        bars = (
            alt.Chart(hist)
            .mark_bar(color=KIND_HUE[kind], cornerRadiusTopLeft=2, cornerRadiusTopRight=2)
            .encode(
                x=alt.X("bin_lo:Q", scale=alt.Scale(type="log"), title="region length (bp, log scale)"),
                x2="bin_hi:Q",
                y=alt.Y("count:Q", title="regions"),
                tooltip=[
                    alt.Tooltip("bin_lo:Q", title="from (bp)", format=".0f"),
                    alt.Tooltip("bin_hi:Q", title="to (bp)", format=".0f"),
                    alt.Tooltip("count:Q", title="regions", format=","),
                ],
            )
        )
        rule = (
            alt.Chart(pl.DataFrame({"median": [median]}))
            .mark_rule(color="#52514e", strokeDash=[4, 3], size=1)
            .encode(x="median:Q", tooltip=alt.Tooltip("median:Q", title="median (bp)"))
        )
        return (
            (bars + rule)
            .properties(width=620, height=260, title=title)
            .configure_axis(grid=True, gridOpacity=0.15, domainOpacity=0.3, tickOpacity=0.3)
            .configure_view(strokeWidth=0)
        )


    return (length_hist_chart,)


@app.cell
def _(length_hist_chart, regions):
    # Generic regions: the sequence with nothing more specific to say about it.
    generic_hist = length_hist_chart(
        regions, "generic", "tomato SL4.00",
        "Generic region lengths — tomato SL4.00 (dashed line = median)",
    )
    generic_hist

    return


@app.cell
def _(length_hist_chart, regions):
    # STR loci: the microsatellites the reference alone hands over as genetic objects.
    ssr_hist = length_hist_chart(
        regions, "ssr_locus", "tomato SL4.00",
        "STR locus lengths — tomato SL4.00 (dashed line = median)",
    )
    ssr_hist

    return


@app.cell
def _(alt, pl, regions):
    # STR loci by period, in 10 bp bins.
    #
    # Exact 1 bp bars left most of the axis structurally empty: a period-p tract can
    # only have a length near a multiple of p, so periods 4-6 populate at most every
    # 4th-6th integer and the chart was more gap than data. A 10 bp bin is wider than
    # every stride in the range (period <= 6), so no bin is empty for want of a
    # reachable length. The comb that binning hides is real but it is arithmetic, not
    # biology — `ssr_by_period_data` below keeps the exact lengths for when it matters.
    SSR_BIN_BP = 10

    ssr_by_period_data = (
        regions.filter(
            (pl.col("kind") == "ssr_locus") & (pl.col("genome") == "tomato SL4.00")
        )
        .group_by("period", "length")
        .agg(pl.len().alias("count"))
        .sort("period", "length")
    )

    ssr_by_period_binned = (
        ssr_by_period_data.with_columns(
            (pl.col("length") // SSR_BIN_BP * SSR_BIN_BP).alias("bin_lo")
        )
        .group_by("period", "bin_lo")
        .agg(pl.col("count").sum().alias("count"))
        .with_columns((pl.col("bin_lo") + SSR_BIN_BP).alias("bin_hi"))
        .sort("period", "bin_lo")
    )

    ssr_by_period = (
        alt.Chart(ssr_by_period_binned)
        .mark_bar(color="#008300", stroke="white", strokeWidth=0.5)
        .encode(
            x=alt.X("bin_lo:Q", scale=alt.Scale(domain=[0, 100]), title="locus length (bp)"),
            x2="bin_hi:Q",
            y=alt.Y("count:Q", scale=alt.Scale(type="symlog"), title="loci (log scale)"),
            tooltip=[
                alt.Tooltip("period:Q", title="period"),
                alt.Tooltip("bin_lo:Q", title="from (bp)"),
                alt.Tooltip("bin_hi:Q", title="to (bp)"),
                alt.Tooltip("count:Q", title="loci", format=","),
            ],
        )
        .properties(width=280, height=150)
        .facet(
            facet=alt.Facet("period:O", title=None, header=alt.Header(labelFontWeight="bold")),
            columns=3,
            title=f"STR locus lengths by motif period, {SSR_BIN_BP} bp bins — tomato SL4.00",
        )
        .configure_axis(grid=True, gridOpacity=0.15, domainOpacity=0.3, tickOpacity=0.3)
        .configure_view(strokeWidth=0)
        .configure_header(labelFontSize=12)
    )
    ssr_by_period

    return


if __name__ == "__main__":
    app.run()
