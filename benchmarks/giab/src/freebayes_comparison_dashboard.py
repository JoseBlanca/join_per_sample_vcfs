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

    # caller label -> result subdir name under results/per_sample/<cov>/
    CALLERS = {
        "ours-high-recall": "high-recall",
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
        mo.md("# freebayes vs ours-high-recall — GIAB per_sample"),
        mo.md(
            "Cohort totals (HG002+HG003+HG004 summed), per coverage tier, "
            "split by variant class. precision = TP/(TP+FP), recall = "
            "TP/(TP+FN), F1 = harmonic mean. freebayes is QUAL ≥ 30 gated; "
            "ours is the high-recall preset (BAQ off, DUST off, "
            "allele-balance filter on)."
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
    COLORS = {"ours-high-recall": "#1f77b4", "freebayes": "#ff7f0e"}

    def _line(metric: str, title: str, y_title: str, zero_to_one: bool):
        y = alt.Y(f"{metric}:Q", title=y_title)
        if zero_to_one:
            y = alt.Y(f"{metric}:Q", title=y_title, scale=alt.Scale(domain=[0.0, 1.0]))
        return (
            alt.Chart(cmp_df)
            .mark_line(point=True)
            .encode(
                x=alt.X("cov_x:Q", title="coverage (x)",
                        scale=alt.Scale(type="log")),
                y=y,
                color=alt.Color("caller:N",
                                scale=alt.Scale(domain=list(COLORS), range=list(COLORS.values()))),
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
    if "ours-high-recall" in cols and "freebayes" in cols:
        wide = wide.with_columns(
            (pl.col("ours-high-recall") - pl.col("freebayes")).round(4).alias("F1_delta(ours-fb)")
        )
    mo.vstack([
        mo.md("## F1 head-to-head (ours − freebayes)"),
        mo.ui.table(wide, selection=None),
        mo.md("_Positive `F1_delta` = ours wins at that coverage/class._"),
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


if __name__ == "__main__":
    app.run()
