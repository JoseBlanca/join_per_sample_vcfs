# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "polars",
#     "altair",
#     "pyarrow",
# ]
# ///
"""GIAB per_sample accuracy — freebayes vs ours (high-recall preset).

Run it:

    uv run marimo run  benchmarks/giab/src/freebayes_comparison_dashboard.py   # app view
    uv run marimo edit benchmarks/giab/src/freebayes_comparison_dashboard.py   # notebook

It compares two callers, per coverage tier, split by SNP vs indel:
  ours-high-recall : results/per_sample/<cov>/high-recall/{sample}.vcf
                     (BAQ off, DUST off, allele-balance filter on)
  freebayes        : results/per_sample/<cov>/freebayes/{sample}.vcf
                     (QUAL >= 30, matching our caller's --min-qual default)

Generate the inputs first:
  PRESET=high-recall benchmarks/giab/src/run_ours_per_sample.sh <COV>
  benchmarks/giab/src/run_freebayes_per_sample.sh <COV>

Evaluation is IDENTICAL to accuracy_dashboard.py (so the numbers line up):
for each (sample, class) the truth and query VCFs are both restricted to the
sample's confident BED, left-aligned + biallelic-split (`norm -m -any`),
class-filtered (snps/indels) and intersected (`isec`, POS+REF+ALT match).
Truth is FILTER=PASS only. Counts are summed across the three samples into
cohort totals, then precision/recall/F1 are recomputed from those totals.

  TP = sites in both truth and query   (isec 0002, truth side)
  FP = sites only in query             (isec 0001)
  FN = sites only in truth             (isec 0000)

A tidy TSV of all (coverage, caller, class) rows is written to
results/per_sample/freebayes_comparison.tsv when the notebook runs.
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

    # caller label -> result subdir name under results/per_sample/<cov>/.
    # Both "ours" rows are the high-recall config (BAQ off, DUST off,
    # allele-balance on); they differ ONLY in the hidden-paralog filter, so the
    # pair isolates that filter's effect. It is ON by default and expected to be
    # inert on GIAB (single-sample, no seg-dups; high-depth samples rejected as
    # out-of-range), so the two "ours" rows should coincide — the check that the
    # filter does no harm on clean human data. Generate the OFF row with:
    #   PRESET=high-recall NO_PARALOG=1 benchmarks/giab/src/run_ours_per_sample.sh <COV>
    # A subdir is skipped if it holds no per-sample VCFs.
    CALLERS = {
        "ours (paralog filter on)": "high-recall",
        "ours (paralog filter off)": "high-recall-noparalog",
        "freebayes": "freebayes",
    }

    BED_OF = {s: f"{s}_bench_azar_merged_100.bed" for s in SAMPLES}
    TRUTH_OF = {
        s: f"{s}_GRCh38_1_22_v4.2.1_benchmark.selected_100.vcf.gz" for s in SAMPLES
    }
    return (
        BED_DIR,
        BED_OF,
        BENCH_DIR,
        CALLERS,
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

    def _normalize(src: Path, cls: str, bed: Path, out: Path, pass_only: bool):
        # NB freebayes emits records in `--targets` BED order, and the per_sample
        # BEDs are in random (non-reference) order, so its VCFs are NOT
        # coordinate-sorted; `norm` left-alignment can also reorder. A final
        # `bcftools sort` makes the output indexable by tabix and is harmless on
        # already-sorted inputs (e.g. our caller's VCFs).
        pf = "-f PASS" if pass_only else ""
        _run(
            f"bcftools view {pf} -T {shlex.quote(str(bed))} {shlex.quote(str(src))} -Ou "
            f"| bcftools norm -f {shlex.quote(str(REFERENCE))} -m -any -Ou 2>/dev/null "
            f"| bcftools view -v {cls} -Ou "
            f"| bcftools sort -Oz -o {shlex.quote(str(out))} 2>/dev/null"
        )
        _run(f"bcftools index -t {shlex.quote(str(out))}")

    def counts_for(truth: Path, query: Path, bed: Path, cls: str):
        """Return (tp, fp, fn) for one (sample, class). Mirrors accuracy_dashboard."""
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
            return tp, fp, fn

    return (counts_for,)


@app.cell
def _(
    BED_DIR,
    BED_OF,
    CALLERS,
    CLASSES,
    RESULTS_DIR,
    SAMPLES,
    TRUTH_DIR,
    TRUTH_OF,
    counts_for,
    pl,
):
    def _cov_x(name: str):
        stem = name[:-1] if name.endswith("x") else name
        return int(stem) if stem.isdigit() else None

    def _build():
        # Wrapped in a function so loop locals don't leak as marimo globals.
        rows = []
        if not RESULTS_DIR.is_dir():
            return pl.DataFrame()
        cov_dirs = sorted(
            (p for p in RESULTS_DIR.iterdir() if p.is_dir() and _cov_x(p.name) is not None),
            key=lambda p: _cov_x(p.name),
        )
        for cov_dir in cov_dirs:
            for caller_label, subdir in CALLERS.items():
                qdir = cov_dir / subdir
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
                        t, f, n = counts_for(
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
                        caller=caller_label, **{"class": cls},
                        TP=tp, FP=fp, FN=fn,
                        precision=round(prec, 4), recall=round(rec, 4), f1=round(f1, 4),
                    ))
        if not rows:
            return pl.DataFrame()
        return pl.DataFrame(rows).sort(["class", "cov_x", "caller"])

    cmp_df = _build()
    return (cmp_df,)


@app.cell
def _(RESULTS_DIR, cmp_df, mo):
    # Emit the tidy TSV alongside the result tree.
    if cmp_df.height:
        tsv_path = RESULTS_DIR / "freebayes_comparison.tsv"
        cmp_df.write_csv(str(tsv_path), separator="\t")
        tsv_md = mo.md(f"Wrote tidy TSV → `{tsv_path}`")
    else:
        tsv_md = mo.md("")
    mo.stop(
        cmp_df.height == 0,
        mo.md(
            "**No comparison rows found.** Generate inputs first:\n\n"
            "```\nPRESET=high-recall benchmarks/giab/src/run_ours_per_sample.sh <COV>\n"
            "benchmarks/giab/src/run_freebayes_per_sample.sh <COV>\n```"
        ),
    )
    tsv_md
    return


@app.cell
def _(cmp_df, mo):
    mo.vstack([
        mo.md("# freebayes vs ours (high-recall) — GIAB per_sample"),
        mo.md(
            "Cohort totals (HG002+HG003+HG004 summed), per coverage tier, "
            "split by variant class. precision = TP/(TP+FP), recall = "
            "TP/(TP+FN), F1 = harmonic mean. freebayes is QUAL ≥ 30 gated; "
            "both **ours** rows are the high-recall preset (BAQ off, DUST off, "
            "allele-balance filter on) and differ ONLY in the hidden-paralog "
            "filter — **paralog filter on** (default) vs **off** "
            "(`--no-paralog-filter`). On GIAB the filter is expected to be inert "
            "(single-sample, no seg-dups; high-depth samples rejected as "
            "out-of-range), so the two ours rows should coincide — the check "
            "that the filter does no harm on clean human data."
        ),
        mo.md("## SNPs"),
        mo.ui.table(cmp_df.filter(cmp_df["class"] == "snps"), selection=None),
        mo.md("## Indels"),
        mo.ui.table(cmp_df.filter(cmp_df["class"] == "indels"), selection=None),
    ])
    return


@app.cell
def _(alt, cmp_df, mo):
    # Side-by-side precision/recall/FP vs coverage, faceted by class, colored by
    # caller. Rendered as SEPARATE altair_chart embeds (one per metric) so Vega
    # doesn't hit "Duplicate signal name" when concatenating shared-scale charts.
    #
    # The two "ours" rows are near-identical (the paralog filter is inert on
    # GIAB), so paralog-OFF is a DASHED line drawn LAST (on top) — otherwise it
    # sits exactly under the solid paralog-ON line and is invisible. Draw order
    # follows the colour-scale domain, so ON is listed first and OFF last.
    ORDER = [
        "ours (paralog filter on)",
        "freebayes",
        "ours (paralog filter off)",
    ]
    COLORS = {
        "ours (paralog filter on)": "#1f77b4",
        "freebayes": "#ff7f0e",
        "ours (paralog filter off)": "#2ca02c",
    }
    DASHES = {
        "ours (paralog filter on)": [1, 0],
        "freebayes": [1, 0],
        "ours (paralog filter off)": [6, 3],
    }

    def _line(metric: str, title: str, y_title: str, zero_to_one: bool):
        y = alt.Y(f"{metric}:Q", title=y_title)
        if zero_to_one:
            y = alt.Y(f"{metric}:Q", title=y_title, scale=alt.Scale(domain=[0.0, 1.0]))
        return (
            alt.Chart(cmp_df)
            .mark_line(point=True, opacity=0.9)
            .encode(
                x=alt.X("cov_x:Q", title="coverage (x)",
                        scale=alt.Scale(type="log")),
                y=y,
                color=alt.Color("caller:N",
                                scale=alt.Scale(domain=ORDER, range=[COLORS[c] for c in ORDER])),
                strokeDash=alt.StrokeDash(
                    "caller:N",
                    scale=alt.Scale(domain=ORDER, range=[DASHES[c] for c in ORDER]),
                    legend=None,
                ),
                column=alt.Column("class:N", title=None),
                tooltip=["coverage", "caller", "class", "TP", "FP", "FN",
                         "precision", "recall", "f1"],
            )
            .properties(width=260, height=240, title=title)
        )

    prec_chart = _line("precision", "Precision vs coverage", "precision", True)
    rec_chart = _line("recall", "Recall vs coverage", "recall", True)
    fp_chart = _line("FP", "False positives vs coverage", "FP (count)", False)
    f1_chart = _line("f1", "F1 vs coverage", "F1", True)

    mo.vstack([
        mo.md("## Precision / recall / FP / F1 across coverage"),
        mo.ui.altair_chart(prec_chart),
        mo.ui.altair_chart(rec_chart),
        mo.ui.altair_chart(fp_chart),
        mo.ui.altair_chart(f1_chart),
        mo.md(
            "_x-axis is log-scaled coverage (5–300x). The allele-balance "
            "filter makes ours' SNP precision climb with depth; the headline "
            "comparison is precision/FP at higher coverage and how each caller "
            "trades recall vs precision across the depth range._"
        ),
    ])
    return


@app.cell
def _(cmp_df, mo, pl):
    # Wide head-to-head: F1 delta (ours - freebayes) per coverage/class.
    wide = (
        cmp_df.pivot(values="f1", index=["class", "coverage", "cov_x"], on="caller")
        .sort(["class", "cov_x"])
    )
    cols = wide.columns
    ours_on = "ours (paralog filter on)"
    if ours_on in cols and "freebayes" in cols:
        wide = wide.with_columns(
            (pl.col(ours_on) - pl.col("freebayes")).round(4).alias("F1_delta(ours-fb)")
        )
    if ours_on in cols and "ours (paralog filter off)" in cols:
        wide = wide.with_columns(
            (pl.col(ours_on) - pl.col("ours (paralog filter off)"))
            .round(4)
            .alias("F1_delta(on-off)")
        )
    mo.vstack([
        mo.md("## F1 head-to-head (ours − freebayes, and paralog on − off)"),
        mo.ui.table(wide, selection=None),
        mo.md(
            "_`F1_delta(ours-fb)` > 0 = ours (paralog on) beats freebayes. "
            "`F1_delta(on-off)` should be ~0 everywhere — the paralog filter is "
            "inert on GIAB._"
        ),
    ])
    return


@app.cell
def _(mo):
    mo.md(
        """
        ---
        **Notes.** Allele-level concordance (POS+REF+ALT), not genotype-level.
        Truth = FILTER PASS only. TN omitted on purpose. freebayes is gated at
        QUAL ≥ 30 to match our caller's `--min-qual` default and shed its
        low-QUAL tail; `{sample}.raw.vcf` (ungated) is kept alongside each
        gated VCF if you want to inspect the raw recall. Low-coverage SNP FNs
        are sampling-limited (the alt allele often isn't sequenced below ~10x),
        so both callers carry high FN at 5–10x.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Variant QUAL distribution vs coverage — TP vs FP
    Box-and-whisker of per-variant `QUAL` at each coverage tier — **freebayes**
    (orange, left) vs **ours** (high-recall, blue, right). Columns split **true
    positives** (match the GIAB truth) from **false positives** (query-only);
    rows split **SNPs** / **indels**. Pooled over HG002/3/4 within their confident
    BEDs; box = IQR, line = median, whiskers = 1.5×IQR clamped to the data.
    **Log y-axis, capped at 10 000.** FP boxes are sparse — both callers are
    high-precision (very few FPs, essentially none for indels), and groups with
    < 3 calls are omitted. Take-away: FP QUAL sits far below TP QUAL at every
    depth, and the gap widens with coverage, so QUAL is a usable TP/FP separator.
    """)
    return


@app.cell
def _(
    BED_DIR,
    BED_OF,
    CLASSES,
    Path,
    REFERENCE,
    RESULTS_DIR,
    SAMPLES,
    TRUTH_DIR,
    TRUTH_OF,
    pl,
    shlex,
    subprocess,
    tempfile,
):
    # QUAL box statistics per (coverage, caller, class, TP/FP). TP/FP come from
    # bcftools isec vs the GIAB truth (same normalize as the concordance above):
    # 0003 = query-side of shared (TP), 0001 = query-only (FP).
    _BOX_CALLERS = {"freebayes": "freebayes", "ours": "high-recall"}

    def _covx(name):
        stem = name[:-1] if name.endswith("x") else name
        return int(stem) if stem.isdigit() else None

    def _run(cmd):
        subprocess.run(cmd, shell=True, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def _norm(src, cls, bed, out, pass_only):
        pf = "-f PASS" if pass_only else ""
        _run(
            f"bcftools view {pf} -T {shlex.quote(str(bed))} {shlex.quote(str(src))} -Ou "
            f"| bcftools norm -f {shlex.quote(str(REFERENCE))} -m -any -Ou 2>/dev/null "
            f"| bcftools view -v {cls} -Ou "
            f"| bcftools sort -Oz -o {shlex.quote(str(out))} 2>/dev/null"
        )
        _run(f"bcftools index -t {shlex.quote(str(out))}")

    def _qof(vcf):
        out = subprocess.run(
            f"bcftools query -f '%QUAL\\n' {shlex.quote(str(vcf))}",
            shell=True, capture_output=True, text=True,
        ).stdout
        vals = []
        for ln in out.splitlines():
            ln = ln.strip()
            if ln and ln != ".":
                try:
                    vals.append(float(ln))
                except ValueError:
                    pass
        return vals

    def _box(vals):
        ss = pl.Series(vals)
        q1 = float(ss.quantile(0.25))
        med = float(ss.quantile(0.5))
        q3 = float(ss.quantile(0.75))
        iqr = q3 - q1
        wlo = max(max(float(ss.min()), q1 - 1.5 * iqr), 1.0)
        whi = max(min(float(ss.max()), q3 + 1.5 * iqr), 1.0)
        return dict(q1=q1, med=med, q3=q3, wlo=wlo, whi=whi, n=len(vals))

    def _build_qbox():
        if not RESULTS_DIR.is_dir():
            return pl.DataFrame()
        cov_dirs = sorted(
            (p for p in RESULTS_DIR.iterdir() if p.is_dir() and _covx(p.name) is not None),
            key=lambda p: _covx(p.name),
        )
        rows = []
        for cov_dir in cov_dirs:
            for label, subdir in _BOX_CALLERS.items():
                qdir = cov_dir / subdir
                if not qdir.is_dir():
                    continue
                for cls in CLASSES:
                    acc = {"TP": [], "FP": []}
                    for s in SAMPLES:
                        q = qdir / f"{s}.vcf"
                        if not q.exists():
                            continue
                        with tempfile.TemporaryDirectory() as td:
                            tdp = Path(td)
                            qn = tdp / "q.vcf.gz"
                            tn = tdp / "t.vcf.gz"
                            _norm(q, cls, BED_DIR / BED_OF[s], qn, pass_only=False)
                            _norm(TRUTH_DIR / TRUTH_OF[s], cls, BED_DIR / BED_OF[s], tn, pass_only=True)
                            isec = tdp / "isec"
                            _run(f"bcftools isec -p {shlex.quote(str(isec))} "
                                 f"{shlex.quote(str(tn))} {shlex.quote(str(qn))}")
                            acc["TP"].extend(_qof(isec / "0003.vcf"))
                            acc["FP"].extend(_qof(isec / "0001.vcf"))
                    for status, vals in acc.items():
                        if len(vals) < 3:
                            continue
                        rows.append(dict(
                            coverage=cov_dir.name, cov_x=_covx(cov_dir.name),
                            caller=label, **{"class": cls}, status=status, **_box(vals),
                        ))
        return pl.DataFrame(rows) if rows else pl.DataFrame()

    qbox_df = _build_qbox()
    return (qbox_df,)


@app.cell
def _(alt, mo, pl, qbox_df):
    def _plot():
        cov_order = list(qbox_df.sort("cov_x")["coverage"].unique(maintain_order=True))
        dom, rng = ["freebayes", "ours"], ["#ff7f0e", "#1f77b4"]
        ysc = alt.Scale(type="log", domain=[10, 10000], clamp=True)

        def _panel(status, cls, col_title, ytitle, legend):
            d = qbox_df.filter((pl.col("status") == status) & (pl.col("class") == cls))
            base = alt.Chart(d, title=col_title).encode(
                x=alt.X("coverage:O", sort=cov_order, title="coverage"),
                xOffset=alt.XOffset("caller:N", sort=dom),
                color=alt.Color(
                    "caller:N", scale=alt.Scale(domain=dom, range=rng),
                    legend=(alt.Legend(title="caller") if legend else None),
                ),
                tooltip=["class", "status", "coverage", "caller", "n",
                         alt.Tooltip("med:Q", format=".0f", title="median"),
                         alt.Tooltip("q1:Q", format=".0f"),
                         alt.Tooltip("q3:Q", format=".0f")],
            )
            rule = base.mark_rule().encode(y=alt.Y("wlo:Q", title=ytitle, scale=ysc), y2="whi:Q")
            box = base.mark_bar(size=13).encode(y=alt.Y("q1:Q", scale=ysc), y2="q3:Q")
            med = base.mark_tick(thickness=2, size=13).encode(
                y=alt.Y("med:Q", scale=ysc), color=alt.value("black"))
            return (rule + box + med).properties(width=300, height=190)

        return alt.vconcat(
            alt.hconcat(
                _panel("TP", "snps", "True positives", "snps · QUAL (log)", False),
                _panel("FP", "snps", "False positives", "", True),
            ),
            alt.hconcat(
                _panel("TP", "indels", "", "indels · QUAL (log)", False),
                _panel("FP", "indels", "", "", False),
            ),
        )

    mo.stop(qbox_df.is_empty(), mo.md("_No QUAL data — run the callers first._"))
    _plot()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Genotype concordance vs GIAB (TP variants)
    For the true-positive variants (called *and* in the GIAB truth), the percentage
    whose **genotype** matches the GIAB genotype (phase-insensitive: `1|0` ≡ `0/1`),
    computed per sample then boxed over HG002/3/4. **freebayes** (orange, left) vs
    **ours** (blue, right); SNPs (top) / indels (bottom); coverage on x. Each box is
    the 3-sample spread (n = 3; samples with < 10 TP calls dropped). At low depth our
    SNP genotyping trails freebayes (5×: ~84 % vs ~95 %) but both reach ~99–100 % by
    30×; indel genotyping is lower for both and our indel GTs plateau below freebayes'.
    """)
    return


@app.cell
def _(
    BED_DIR,
    BED_OF,
    CLASSES,
    Path,
    REFERENCE,
    RESULTS_DIR,
    SAMPLES,
    TRUTH_DIR,
    TRUTH_OF,
    pl,
    shlex,
    subprocess,
    tempfile,
):
    # Per-sample genotype concordance for TP variants: isec vs truth, then compare
    # the GT on the query side (0003) to the truth side (0002), phase-insensitive.
    _BOX_CALLERS = {"freebayes": "freebayes", "ours": "high-recall"}

    def _covx(name):
        stem = name[:-1] if name.endswith("x") else name
        return int(stem) if stem.isdigit() else None

    def _run(cmd):
        subprocess.run(cmd, shell=True, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def _norm(src, cls, bed, out, pass_only):
        pf = "-f PASS" if pass_only else ""
        _run(
            f"bcftools view {pf} -T {shlex.quote(str(bed))} {shlex.quote(str(src))} -Ou "
            f"| bcftools norm -f {shlex.quote(str(REFERENCE))} -m -any -Ou 2>/dev/null "
            f"| bcftools view -v {cls} -Ou "
            f"| bcftools sort -Oz -o {shlex.quote(str(out))} 2>/dev/null"
        )
        _run(f"bcftools index -t {shlex.quote(str(out))}")

    def _gtmap(vcf):
        out = subprocess.run(
            f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n' {shlex.quote(str(vcf))}",
            shell=True, capture_output=True, text=True,
        ).stdout
        m = {}
        for ln in out.splitlines():
            p = ln.split("\t")
            if len(p) >= 5:
                m[f"{p[0]}:{p[1]}:{p[2]}:{p[3]}"] = p[4]
        return m

    def _norm_gt(gt):
        a = [x for x in gt.replace("|", "/").split("/") if x not in (".", "")]
        if not a:
            return None
        try:
            return "/".join(map(str, sorted(int(x) for x in a)))
        except ValueError:
            return None

    def _concordance(truth_vcf, query_vcf):
        tm, qm = _gtmap(truth_vcf), _gtmap(query_vcf)
        match = total = 0
        for k, qg in qm.items():
            tg = tm.get(k)
            if tg is None:
                continue
            ng_q, ng_t = _norm_gt(qg), _norm_gt(tg)
            if ng_q is None or ng_t is None:
                continue
            total += 1
            if ng_q == ng_t:
                match += 1
        return (100.0 * match / total, total) if total else (None, 0)

    def _box(vals):
        ss = pl.Series(vals)
        return dict(
            q1=float(ss.quantile(0.25)), med=float(ss.quantile(0.5)),
            q3=float(ss.quantile(0.75)), wlo=float(ss.min()), whi=float(ss.max()), n=len(vals),
        )

    def _build_gtbox():
        if not RESULTS_DIR.is_dir():
            return pl.DataFrame()
        cov_dirs = sorted(
            (p for p in RESULTS_DIR.iterdir() if p.is_dir() and _covx(p.name) is not None),
            key=lambda p: _covx(p.name),
        )
        rows = []
        for cov_dir in cov_dirs:
            for label, subdir in _BOX_CALLERS.items():
                qdir = cov_dir / subdir
                if not qdir.is_dir():
                    continue
                for cls in CLASSES:
                    rates = []
                    for s in SAMPLES:
                        q = qdir / f"{s}.vcf"
                        if not q.exists():
                            continue
                        with tempfile.TemporaryDirectory() as td:
                            tdp = Path(td)
                            qn, tn = tdp / "q.vcf.gz", tdp / "t.vcf.gz"
                            _norm(q, cls, BED_DIR / BED_OF[s], qn, pass_only=False)
                            _norm(TRUTH_DIR / TRUTH_OF[s], cls, BED_DIR / BED_OF[s], tn, pass_only=True)
                            isec = tdp / "isec"
                            _run(f"bcftools isec -p {shlex.quote(str(isec))} "
                                 f"{shlex.quote(str(tn))} {shlex.quote(str(qn))}")
                            rate, tot = _concordance(isec / "0002.vcf", isec / "0003.vcf")
                        if rate is not None and tot >= 10:
                            rates.append(rate)
                    if len(rates) >= 2:
                        rows.append(dict(
                            coverage=cov_dir.name, cov_x=_covx(cov_dir.name),
                            caller=label, **{"class": cls}, **_box(rates),
                        ))
        return pl.DataFrame(rows) if rows else pl.DataFrame()

    gtbox_df = _build_gtbox()
    return (gtbox_df,)


@app.cell
def _(alt, mo, gtbox_df):
    def _plot():
        cov_order = list(gtbox_df.sort("cov_x")["coverage"].unique(maintain_order=True))
        dom, rng = ["freebayes", "ours"], ["#ff7f0e", "#1f77b4"]
        ysc = alt.Scale(domain=[50, 100])
        base = alt.Chart(gtbox_df).encode(
            x=alt.X("coverage:O", sort=cov_order, title="coverage"),
            xOffset=alt.XOffset("caller:N", sort=dom),
            color=alt.Color("caller:N", scale=alt.Scale(domain=dom, range=rng), title="caller"),
            tooltip=["class", "coverage", "caller", "n",
                     alt.Tooltip("med:Q", format=".1f", title="median %"),
                     alt.Tooltip("wlo:Q", format=".1f", title="min %"),
                     alt.Tooltip("whi:Q", format=".1f", title="max %")],
        )
        rule = base.mark_rule().encode(
            y=alt.Y("wlo:Q", title="GT concordance vs GIAB (%)", scale=ysc), y2="whi:Q")
        box = base.mark_bar(size=14).encode(y=alt.Y("q1:Q", scale=ysc), y2="q3:Q")
        med = base.mark_tick(thickness=2, size=14).encode(
            y=alt.Y("med:Q", scale=ysc), color=alt.value("black"))
        return (rule + box + med).properties(width=380, height=210).facet(
            row=alt.Row("class:N", sort=["snps", "indels"], title=None))

    mo.stop(gtbox_df.is_empty(), mo.md("_No genotype data — run the callers first._"))
    _plot()
    return


if __name__ == "__main__":
    app.run()
