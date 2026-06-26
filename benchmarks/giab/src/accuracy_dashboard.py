# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "polars",
#     "altair",
#     "pyarrow",
# ]
# ///
"""GIAB per_sample accuracy dashboard — TP / FP / FN split by SNP vs indel.

Run it:

    uv run marimo run   benchmarks/giab/src/accuracy_dashboard.py   # app view
    uv run marimo edit  benchmarks/giab/src/accuracy_dashboard.py   # notebook

Method (bcftools, mirrors benchmarks/lib/compare_to_truth.sh): for each
sample the truth callset and our VCF are both restricted to that sample's
confident BED, left-aligned and biallelic-split (`norm -m -any`), then
class-filtered (snps / indels) and intersected (`isec`, matching on
POS+REF+ALT). Truth is restricted to FILTER=PASS; our records keep their
own FILTER. This scores ALLELE concordance, not genotype concordance
(a 0/1-vs-1/1 mismatch at a correct site still counts as a TP) — the
dependency-light analogue of rtg vcfeval / hap.py.

Counts are per (sample, class):
  TP = sites in both truth and ours   (isec 0002, truth side)
  FP = sites only in ours             (isec 0001)
  FN = sites only in truth            (isec 0000)
True negatives are intentionally omitted — over random confident regions
the no-variant denominator is enormous and uninformative.
"""

import marimo

__generated_with = "0.16.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import os
    import shlex
    import shutil
    import subprocess
    import tempfile
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import polars as pl

    # benchmarks/giab/src/this.py -> benchmarks/giab
    BENCH_DIR = Path(
        os.environ.get("PVC_GIAB_BENCH_DIR", Path(__file__).resolve().parents[1])
    )
    REFERENCE = (
        BENCH_DIR
        / "ref_genome_GRCh38"
        / "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
    )
    BED_DIR = BENCH_DIR / "per_sample" / "bed"
    TRUTH_DIR = BENCH_DIR / "per_sample" / "vcf"
    RESULTS_DIR = BENCH_DIR / "results" / "per_sample"

    SAMPLES = ["HG002", "HG003", "HG004"]
    CLASSES = ["snps", "indels"]

    # sample -> (BED basename, truth VCF basename)
    BED_OF = {s: f"{s}_bench_azar_merged_100.bed" for s in SAMPLES}
    TRUTH_OF = {
        s: f"{s}_GRCh38_1_22_v4.2.1_benchmark.selected_100.vcf.gz" for s in SAMPLES
    }
    return (
        BED_DIR,
        BED_OF,
        BENCH_DIR,
        CLASSES,
        REFERENCE,
        RESULTS_DIR,
        SAMPLES,
        TRUTH_DIR,
        TRUTH_OF,
        Path,
        alt,
        mo,
        pl,
        shlex,
        shutil,
        subprocess,
        tempfile,
    )


@app.cell
def _(mo, shutil):
    bcftools = shutil.which("bcftools")
    mo.stop(
        bcftools is None,
        mo.md("**`bcftools` not found on PATH.** Install it (brew/conda) to compute concordance."),
    )
    mo.md(f"`bcftools`: `{bcftools}`")
    return (bcftools,)


@app.cell
def _(RESULTS_DIR, mo):
    # Discover available (coverage, variant) result dirs, e.g.
    # results/per_sample/300x/ours_nobaq_nodust/
    coverages = sorted(p.name for p in RESULTS_DIR.iterdir() if p.is_dir()) if RESULTS_DIR.is_dir() else []
    mo.stop(
        not coverages,
        mo.md(f"No result dirs under `{RESULTS_DIR}`. Run `run_ours_per_sample.sh` first."),
    )
    cov_dd = mo.ui.dropdown(
        options=coverages,
        value="300x" if "300x" in coverages else coverages[0],
        label="Coverage",
    )
    return cov_dd, coverages


@app.cell
def _(RESULTS_DIR, cov_dd, mo):
    cov_path = RESULTS_DIR / cov_dd.value
    # Any of our caller-output variant dirs: the `ours*` toggles plus the named
    # presets (high-recall / high-confidence). A dir qualifies if it holds at
    # least one per-sample VCF (excludes empty/aux dirs and other callers).
    OTHER_CALLERS = {"gatk", "freebayes", "comparison"}
    variants = sorted(
        p.name
        for p in cov_path.iterdir()
        if p.is_dir()
        and p.name not in OTHER_CALLERS
        and any(p.glob("HG*.vcf"))
    )
    mo.stop(not variants, mo.md(f"No caller variant dirs with VCFs under `{cov_path}`."))
    # Prefer the high-recall preset (BAQ/DUST off, allele-balance on) by default.
    default_variant = "high-recall" if "high-recall" in variants else variants[0]
    var_dd = mo.ui.dropdown(options=variants, value=default_variant, label="Caller variant")
    return cov_path, var_dd, variants


@app.cell
def _(cov_dd, mo, var_dd):
    mo.hstack([cov_dd, var_dd], justify="start", gap=2)
    return


@app.cell
def _(REFERENCE, shlex, subprocess, tempfile, Path):
    def _run(cmd: str):
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def _count(vcf: Path) -> int:
        out = subprocess.run(
            f"bcftools view -H {shlex.quote(str(vcf))}",
            shell=True, check=True, capture_output=True, text=True,
        ).stdout
        return sum(1 for line in out.splitlines() if line.strip())

    def _quals(vcf: Path):
        """QUAL values (floats) of every record in a VCF; skips '.'."""
        out = subprocess.run(
            f"bcftools query -f '%QUAL\\n' {shlex.quote(str(vcf))}",
            shell=True, check=True, capture_output=True, text=True,
        ).stdout
        vals = []
        for line in out.splitlines():
            line = line.strip()
            if line and line != ".":
                try:
                    vals.append(float(line))
                except ValueError:
                    pass
        return vals

    def _normalize(src: Path, cls: str, bed: Path, out: Path, pass_only: bool):
        pf = "-f PASS" if pass_only else ""
        _run(
            f"bcftools view {pf} -T {shlex.quote(str(bed))} {shlex.quote(str(src))} -Ou "
            f"| bcftools norm -f {shlex.quote(str(REFERENCE))} -m -any -Ou 2>/dev/null "
            f"| bcftools view -v {cls} -Oz -o {shlex.quote(str(out))}"
        )
        _run(f"bcftools index -t {shlex.quote(str(out))}")

    def counts_for(truth: Path, query: Path, bed: Path, cls: str):
        """Return (tp, fp, fn, tp_quals, fp_quals) for one (sample, class).

        Counts come from the truth side (TP=0002) and the query-only set
        (FP=0001). QUAL distributions come from the caller's OWN records:
        0003 = shared, query side (our record for each TP); 0001 = FP.
        """
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            truth_n = tdp / f"truth.{cls}.vcf.gz"
            query_n = tdp / f"query.{cls}.vcf.gz"
            _normalize(truth, cls, bed, truth_n, pass_only=True)
            _normalize(query, cls, bed, query_n, pass_only=False)
            isec = tdp / "isec"
            _run(
                f"bcftools isec -p {shlex.quote(str(isec))} "
                f"{shlex.quote(str(truth_n))} {shlex.quote(str(query_n))}"
            )
            fn = _count(isec / "0000.vcf")  # truth-only
            fp = _count(isec / "0001.vcf")  # query-only
            tp = _count(isec / "0002.vcf")  # shared (truth side)
            tp_quals = _quals(isec / "0003.vcf")  # our record for each TP
            fp_quals = _quals(isec / "0001.vcf")  # FP records
            return tp, fp, fn, tp_quals, fp_quals

    return (counts_for,)


@app.cell
def _(
    BED_DIR,
    BED_OF,
    CLASSES,
    SAMPLES,
    TRUTH_DIR,
    TRUTH_OF,
    counts_for,
    cov_path,
    mo,
    pl,
    var_dd,
):
    query_dir = cov_path / var_dd.value
    rows = []
    qual_rows = []
    missing = []
    for s in SAMPLES:
        q = query_dir / f"{s}.vcf"
        if not q.exists():
            missing.append(s)
            continue
        bed = BED_DIR / BED_OF[s]
        truth = TRUTH_DIR / TRUTH_OF[s]
        for cls in CLASSES:
            tp, fp, fn, tp_quals, fp_quals = counts_for(truth, q, bed, cls)
            prec = tp / (tp + fp) if (tp + fp) else 0.0
            rec = tp / (tp + fn) if (tp + fn) else 0.0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
            rows.append(
                dict(sample=s, **{"class": cls}, TP=tp, FP=fp, FN=fn,
                     precision=round(prec, 4), recall=round(rec, 4), f1=round(f1, 4))
            )
            for v in tp_quals:
                qual_rows.append(dict(sample=s, **{"class": cls}, status="TP", qual=v))
            for v in fp_quals:
                qual_rows.append(dict(sample=s, **{"class": cls}, status="FP", qual=v))

    df = pl.DataFrame(rows)
    qual_df = pl.DataFrame(
        qual_rows,
        schema={"sample": pl.Utf8, "class": pl.Utf8, "status": pl.Utf8, "qual": pl.Float64},
    )
    warn = mo.md(f"⚠️ missing VCFs for: {', '.join(missing)}") if missing else mo.md("")
    return df, qual_df, query_dir, warn


@app.cell
def _(df, mo, pl, query_dir, var_dd, warn):
    # Aggregate across samples, per class (sum counts, recompute rates).
    agg = (
        df.group_by("class")
        .agg(pl.col("TP").sum(), pl.col("FP").sum(), pl.col("FN").sum())
        .with_columns(
            precision=(pl.col("TP") / (pl.col("TP") + pl.col("FP"))).round(4),
            recall=(pl.col("TP") / (pl.col("TP") + pl.col("FN"))).round(4),
        )
        .with_columns(
            f1=(2 * pl.col("precision") * pl.col("recall")
                / (pl.col("precision") + pl.col("recall"))).round(4)
        )
        .sort("class")
    )
    mo.vstack([
        mo.md(f"## Accuracy — `{var_dd.value}`\n\n`{query_dir}`"),
        warn,
        mo.md("### Cohort totals (all samples summed, per class)"),
        mo.ui.table(agg, selection=None),
        mo.md("### Per-sample"),
        mo.ui.table(df.sort(["class", "sample"]), selection=None),
    ])
    return (agg,)


@app.cell
def _(BED_DIR, BED_OF, CLASSES, RESULTS_DIR, SAMPLES, TRUTH_DIR, TRUTH_OF,
      counts_for, mo, pl, var_dd):
    # TP/FP/FN per coverage for the selected variant (cohort totals across the
    # three samples). Recomputed for every coverage dir that holds this variant,
    # so it's heavier than the single-coverage panels above.
    variant = var_dd.value

    def _cov_x(name):
        stem = name[:-1] if name.endswith("x") else name
        return int(stem) if stem.isdigit() else None

    def _build_rows():
        # Wrapped in a function so the loop locals don't leak as marimo
        # cell-level globals (which would collide with the per-coverage cell).
        rows = []
        for cov_dir in sorted(RESULTS_DIR.iterdir(), key=lambda p: _cov_x(p.name) or 0):
            qdir = cov_dir / variant
            if not qdir.is_dir():
                continue
            for cls in CLASSES:
                tp = fp = fn = 0
                present = False
                for s in SAMPLES:
                    q = qdir / f"{s}.vcf"
                    if not q.exists():
                        continue
                    present = True
                    t, f, n, _, _ = counts_for(
                        TRUTH_DIR / TRUTH_OF[s], q, BED_DIR / BED_OF[s], cls
                    )
                    tp += t
                    fp += f
                    fn += n
                if not present:
                    continue
                prec = tp / (tp + fp) if (tp + fp) else 0.0
                rec = tp / (tp + fn) if (tp + fn) else 0.0
                f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
                rows.append(dict(
                    coverage=cov_dir.name, cov_x=_cov_x(cov_dir.name),
                    **{"class": cls}, TP=tp, FP=fp, FN=fn,
                    precision=round(prec, 4), recall=round(rec, 4), f1=round(f1, 4),
                ))
        return rows

    cov_table = pl.DataFrame(_build_rows())
    mo.stop(
        cov_table.height == 0,
        mo.md(f"### Accuracy across coverage — no other coverages found for `{variant}`."),
    )
    cov_table = cov_table.sort(["class", "cov_x"])
    mo.vstack([
        mo.md(f"### Accuracy across coverage — `{variant}` (cohort totals, per class)"),
        mo.ui.table(cov_table, selection=None),
        mo.md("_TP/FP/FN summed over HG002/3/4. The allele-balance filter is "
              "depth-aware, so FP falls with coverage while FN stays low._"),
    ])
    return


@app.cell
def _(alt, df, mo, pl):
    # Long form for a grouped bar chart of TP/FP/FN.
    long = df.unpivot(
        index=["sample", "class"],
        on=["TP", "FP", "FN"],
        variable_name="status",
        value_name="count",
    )
    chart = (
        alt.Chart(long)
        .mark_bar()
        .encode(
            x=alt.X("status:N", title=None, sort=["TP", "FP", "FN"]),
            y=alt.Y("count:Q"),
            color=alt.Color("status:N", sort=["TP", "FP", "FN"],
                            scale=alt.Scale(domain=["TP", "FP", "FN"],
                                            range=["#2ca02c", "#d62728", "#ff7f0e"])),
            column=alt.Column("class:N", title=None),
            xOffset="sample:N",
            tooltip=["sample", "class", "status", "count"],
        )
        .properties(width=180, height=260, title="Counts by class (grouped by sample)")
    )
    mo.ui.altair_chart(chart)
    return


@app.cell
def _(CLASSES, mo):
    # Interactive QUAL distribution: pick the variant class to inspect.
    class_radio = mo.ui.radio(
        options=CLASSES, value="snps", label="Variant class", inline=True
    )
    clip_switch = mo.ui.switch(value=True, label="Clip x-axis to 99th pct (readability)")
    return class_radio, clip_switch


@app.cell
def _(alt, class_radio, clip_switch, mo, pl, qual_df):
    mo.stop(
        qual_df.height == 0,
        mo.md("### QUAL distribution\n\n_No QUAL values available._"),
    )

    sel = qual_df.filter(pl.col("class") == class_radio.value)
    tp_n = sel.filter(pl.col("status") == "TP").height
    fp_n = sel.filter(pl.col("status") == "FP").height

    x = alt.X("qual:Q", title="QUAL", bin=alt.Bin(maxbins=60))
    if clip_switch.value and sel.height:
        hi = float(sel.select(pl.col("qual").quantile(0.99)).item() or 0.0)
        if hi > 0:
            x = alt.X("qual:Q", title="QUAL (clipped at p99)",
                      bin=alt.Bin(maxbins=60, extent=[0.0, hi]),
                      scale=alt.Scale(domain=[0.0, hi], clamp=True))

    hist = (
        alt.Chart(sel)
        .mark_bar(opacity=0.55)
        .encode(
            x=x,
            y=alt.Y("count():Q", title="records", stack=None),
            color=alt.Color("status:N", sort=["TP", "FP"],
                            scale=alt.Scale(domain=["TP", "FP"],
                                            range=["#2ca02c", "#d62728"])),
            tooltip=["status:N", alt.Tooltip("count():Q", title="records")],
        )
        .properties(width=560, height=300,
                    title=f"QUAL distribution — {class_radio.value} (TP n={tp_n}, FP n={fp_n}, all samples)")
    )
    mo.vstack([
        mo.md("### QUAL distribution by call status"),
        mo.hstack([class_radio, clip_switch], justify="start", gap=2),
        mo.ui.altair_chart(hist),
        mo.md(
            "Overlaid (semi-transparent) histograms: green = true positives, "
            "red = false positives. A QUAL threshold that cleanly separates the "
            "two is the gateable-FP signal; FPs piled up at low QUAL can be "
            "filtered, FPs overlapping the TP mass cannot."
        ),
    ])
    return


@app.cell
def _(mo):
    mo.md(
        """
        ---
        **Notes.** Allele-level concordance (POS+REF+ALT), not genotype-level.
        Truth = FILTER PASS only; our calls keep their own FILTER. TN omitted
        on purpose. SNPs and indels are scored separately because BAQ + the DUST
        low-complexity filter (both on by default in our caller) suppress true
        indels — the `ours_nobaq_nodust` variant turns both off.
        """
    )
    return


if __name__ == "__main__":
    app.run()
