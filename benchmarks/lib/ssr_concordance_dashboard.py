# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "matplotlib",
#     "matplotlib-venn",
# ]
# ///

# Marimo dashboard — SSR genotyping concordance, ours vs HipSTR.
#
# Explores agreement between our SSR caller and HipSTR on a shared locus set
# (no STR truth set exists for tomato). The two tools define an STR allele
# differently — HipSTR calls sequence haplotypes (an in-repeat SNP is an
# allele), we call repeat LENGTH — so "variable" here means LENGTH-variable: a
# locus where some sample carries an allele whose length differs from the
# reference tract. Defined that way, HipSTR's seq-only SNP loci collapse to
# length-0 and drop out, so the comparison is apples-to-apples.
#
# Reads the two cohort VCFs a benchmark run produced:
#   <results*>/ours/cohort/cohort.ssr.vcf      (ours)
#   <results*>/hipstr/cohort.str.vcf.gz        (HipSTR)
# and offers a selector over every results* dir that has both.
#
# Recommended invocation:
#
#   uvx marimo edit --sandbox benchmarks/lib/ssr_concordance_dashboard.py
#
# (or `marimo run` for a read-only view).

import marimo

__generated_with = "0.23.10"
app = marimo.App(width="medium")


@app.cell
def _():
    import gzip
    from collections import defaultdict
    from pathlib import Path

    import marimo as mo

    return Path, defaultdict, gzip, mo


@app.cell
def _(mo):
    mo.md("""
    # SSR genotyping concordance — ours vs HipSTR

    Agreement on a shared locus set (no truth set → agreement, not accuracy).
    **Length-variable** = a locus where some sample carries an allele whose
    length differs from the reference tract. This is the apples-to-apples
    notion: HipSTR's in-repeat-SNP alleles are REF-length, so they collapse
    out and only genuine repeat-length polymorphism is compared.

    Genotype unit: the sorted per-allele bp-difference-from-REF (== HipSTR's
    `GB`); two genotypes agree iff those multisets are equal.
    """)
    return


@app.cell
def _(Path):
    # benchmarks/lib/ssr_concordance_dashboard.py -> lib -> benchmarks
    _benchmarks = Path(__file__).resolve().parent.parent
    _ssr = _benchmarks / "ssr_tomato1"
    # Every results* dir that has both cohort VCFs.
    pairs = {}
    for _rd in sorted(_ssr.glob("results*")):
        _ours = _rd / "ours" / "cohort" / "cohort.ssr.vcf"
        _hip = _rd / "hipstr" / "cohort.str.vcf.gz"
        if _ours.exists() and _hip.exists():
            pairs[_rd.name] = (_ours, _hip)
    return (pairs,)


@app.cell
def _(mo, pairs):
    if not pairs:
        run_sel = None
        _view = mo.md(
            "_No `results*/` dir with both `ours/cohort/cohort.ssr.vcf` and "
            "`hipstr/cohort.str.vcf.gz` found. Run the benchmark first._"
        )
    else:
        # Prefer the richest run (most loci ~ largest name suffix) as default.
        _default = "results_ssr15k" if "results_ssr15k" in pairs else next(iter(pairs))
        run_sel = mo.ui.dropdown(
            options=list(pairs.keys()), value=_default, label="benchmark run"
        )
        _view = run_sel
    _view
    return (run_sel,)


@app.cell
def _(gzip):
    # --- VCF parsing + the length-genotype model (shared by every cell) ---
    PERIOD_NAME = {1: "mono", 2: "di", 3: "tri", 4: "tetra", 5: "penta",
                   6: "hexa", 7: "hepta", 8: "octa", 9: "nona"}

    def period_label(p):
        return PERIOD_NAME.get(p, f"{p}-mer") if p else "?"

    def _open(path):
        p = str(path)
        return gzip.open(p, "rt") if p.endswith((".gz", ".bgz")) else open(p)

    def _info_int(info, key):
        for kv in info.split(";"):
            if kv.startswith(key + "="):
                try:
                    return int(kv.split("=")[1])
                except ValueError:
                    return None
        return None

    def _gt_indices(gt):
        parts = gt.replace("|", "/").split("/")
        if any(p in (".", "") for p in parts):
            return None
        return [int(p) for p in parts]

    def parse_vcf(path):
        """-> (samples, {(chrom,pos): rec}). rec carries period + each sample's
        rel-genotype: the sorted tuple of (len(allele)-len(REF)) over its GT,
        or None if missing."""
        samples = []
        loci = {}
        with _open(path) as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]
                    continue
                f = line.rstrip("\n").split("\t")
                chrom, pos, ref, alt, info, fmt = (
                    f[0], int(f[1]), f[3], f[4], f[7], f[8]
                )
                alleles = [ref] + ([] if alt == "." else alt.split(","))
                rlen = len(ref)
                period = _info_int(info, "PERIOD")
                gi = fmt.split(":").index("GT")
                rel = {}
                for name, col in zip(samples, f[9:]):
                    idx = _gt_indices(col.split(":")[gi])
                    rel[name] = (
                        None if idx is None
                        else tuple(sorted(len(alleles[i]) - rlen for i in idx))
                    )
                loci[(chrom, pos)] = {"period": period, "rel": rel}
        return samples, loci

    def locus_is_length_variable(rec):
        """Some sample carries an allele whose length != REF (rel != 0)."""
        for g in rec["rel"].values():
            if g is not None and any(a != 0 for a in g):
                return True
        return False

    def locus_alleles(rec):
        """Distinct rel-allele lengths observed across called samples (+ REF)."""
        s = {0}
        for g in rec["rel"].values():
            if g is not None:
                s.update(g)
        return s

    return locus_alleles, locus_is_length_variable, parse_vcf, period_label


@app.cell
def _(mo, pairs, parse_vcf, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    _ours_path, _hip_path = pairs[run_sel.value]
    ours_samples, ours_loci = parse_vcf(_ours_path)
    hip_samples, hip_loci = parse_vcf(_hip_path)
    shared_samples = [s for s in ours_samples if s in set(hip_samples)]
    mo.md(
        f"**{run_sel.value}** — ours: {len(ours_loci):,} loci / "
        f"{len(ours_samples)} samples · HipSTR: {len(hip_loci):,} loci / "
        f"{len(hip_samples)} samples · **{len(shared_samples)} paired samples**"
    )
    return hip_loci, ours_loci, shared_samples


@app.cell
def _(hip_loci, locus_is_length_variable, mo, ours_loci, run_sel):
    # Length-variable locus sets per tool, joined on (chrom,pos).
    mo.stop(run_sel is None, mo.md(""))
    ours_var = {k for k, r in ours_loci.items() if locus_is_length_variable(r)}
    hip_var = {k for k, r in hip_loci.items() if locus_is_length_variable(r)}
    # HipSTR genotyped these but found no length variation (every sample at the
    # reference length). Disjoint from hip_var; together they are everything
    # HipSTR emitted. The split separates a real disagreement (ours variant,
    # HipSTR monomorphic) from "HipSTR never genotyped the locus" (skipped).
    hip_mono = set(hip_loci) - hip_var
    both_var = ours_var & hip_var
    ours_only = ours_var - hip_var
    hip_only = hip_var - ours_var
    return both_var, hip_mono, hip_only, hip_var, ours_only, ours_var


@app.cell
def _():
    import matplotlib.pyplot as plt

    # ours = blue, HipSTR-variable = red, HipSTR-monomorphic = grey.
    OURS_C, HIP_C, MONO_C = "#1f77b4", "#d62728", "#999999"
    return HIP_C, MONO_C, OURS_C, plt


@app.cell
def _(plt):
    from matplotlib_venn import venn3

    def venn3_plot(sets, labels, colors, title):
        """Area-weighted 3-set Venn; each region is labelled with its true
        locus count. HipSTR-variable and HipSTR-monomorphic are disjoint, so
        their overlap regions read 0 by construction."""
        fig, ax = plt.subplots(figsize=(8.5, 6.4))
        v = venn3(subsets=[set(s) for s in sets], set_labels=labels,
                  set_colors=colors, alpha=0.5, ax=ax)
        for lbl in v.set_labels or []:
            if lbl:
                lbl.set_fontsize(11)
                lbl.set_fontweight("bold")
        for lbl in v.subset_labels or []:
            if lbl:
                lbl.set_fontsize(11)
        ax.set_title(title, fontsize=12)
        fig.tight_layout()
        return fig

    return (venn3_plot,)


@app.cell
def _(mo):
    mo.md("""
    ---
    ## 1. Length-variable loci — ours vs HipSTR (HipSTR split variable / mono)

    Three circles: **ours variable**, **HipSTR variable**, and **HipSTR
    monomorphic** (loci HipSTR genotyped but found at the reference length in
    every sample). HipSTR's two circles are mutually exclusive — a locus is one
    or the other — so their overlap reads 0. The split tells apart *why* a locus
    is ours-only:

    - **ours ∩ HipSTR-variable** — both agree there is length variation.
    - **ours ∩ HipSTR-monomorphic** — genuine disagreement: we call a variant,
      HipSTR genotyped the same locus and saw none.
    - **ours only** (outside both HipSTR circles) — HipSTR never genotyped the
      locus (skipped: too few reads / tract too long).
    """)
    return


@app.cell
def _(HIP_C, MONO_C, OURS_C, hip_mono, hip_var, mo, ours_var, run_sel, venn3_plot):
    mo.stop(run_sel is None, mo.md(""))
    venn3_plot(
        [ours_var, hip_var, hip_mono],
        ("ours\nvariable", "HipSTR\nvariable", "HipSTR\nmonomorphic"),
        (OURS_C, HIP_C, MONO_C),
        "Length-variable SSR loci — ours vs HipSTR (split by HipSTR's call)",
    )
    return


@app.cell
def _(hip_mono, hip_var, mo, ours_var, run_sel):
    # Exact region counts — always legible even when the big monomorphic
    # circle squashes the small overlaps in the Venn above.
    mo.stop(run_sel is None, mo.md(""))
    hip_emit = hip_var | hip_mono
    agree = len(ours_var & hip_var)
    disagree = len(ours_var & hip_mono)
    ours_skipped = len(ours_var - hip_emit)
    hipvar_only = len(hip_var - ours_var)
    monoonly = len(hip_mono - ours_var)
    _md = "### Venn regions\n\n| region | loci | meaning |\n|---|---:|---|\n"
    _md += f"| ours ∩ HipSTR-variable | {agree:,} | both agree there's length variation |\n"
    _md += (f"| ours ∩ HipSTR-monomorphic | {disagree:,} | **disagreement** — we call a "
            f"variant, HipSTR genotyped it and saw none |\n")
    _md += (f"| ours only (HipSTR absent) | {ours_skipped:,} | HipSTR never genotyped it "
            f"(skipped: too few reads / tract too long) |\n")
    _md += f"| HipSTR-variable only | {hipvar_only:,} | HipSTR variant, ours not length-variable |\n"
    _md += f"| HipSTR-monomorphic only | {monoonly:,} | both agree monomorphic (the boring bulk) |\n"
    mo.md(_md)
    return


@app.cell
def _(mo):
    mo.md("""
    ## 2. Length-variable loci per motif length (mono / di / tri / …)
    """)
    return


@app.cell
def _(
    HIP_C,
    OURS_C,
    both_var,
    defaultdict,
    hip_loci,
    hip_only,
    mo,
    ours_loci,
    ours_only,
    period_label,
    plt,
    run_sel,
):
    mo.stop(run_sel is None, mo.md(""))

    def _period_of(k):
        r = ours_loci.get(k) or hip_loci.get(k)
        return r["period"] if r else None

    by_p = defaultdict(lambda: [0, 0, 0])  # period -> [both, ours_only, hip_only]
    for _k in both_var:
        by_p[_period_of(_k)][0] += 1
    for _k in ours_only:
        by_p[_period_of(_k)][1] += 1
    for _k in hip_only:
        by_p[_period_of(_k)][2] += 1
    periods_sorted = sorted(p for p in by_p if p is not None)

    _fig, _ax = plt.subplots(figsize=(9, 4.8))
    _x = range(len(periods_sorted))
    _both = [by_p[p][0] for p in periods_sorted]
    _oo = [by_p[p][1] for p in periods_sorted]
    _ho = [by_p[p][2] for p in periods_sorted]
    _w = 0.27
    _ax.bar([i - _w for i in _x], _both, _w, label="both", color="#6a3d9a")
    _ax.bar(list(_x), _oo, _w, label="ours only", color=OURS_C)
    _ax.bar([i + _w for i in _x], _ho, _w, label="HipSTR only", color=HIP_C)
    _ax.set_xticks(list(_x))
    _ax.set_xticklabels([f"{period_label(p)}\n(p{p})" for p in periods_sorted])
    _ax.set_ylabel("length-variable loci")
    _ax.set_title("Length-variable loci by motif length")
    _ax.legend()
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return by_p, periods_sorted


@app.cell
def _(by_p, mo, period_label, periods_sorted, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    _md = "### Per motif length\n\n"
    _md += "| motif | period | both | ours only | HipSTR only | ours total | HipSTR total | Jaccard |\n"
    _md += "|---|---:|---:|---:|---:|---:|---:|---:|\n"
    _tb = _to = _th = 0
    for _p in periods_sorted:
        _b, _oo, _ho = by_p[_p]
        _ot, _ht = _b + _oo, _b + _ho
        _union = _b + _oo + _ho
        _j = _b / _union if _union else 0.0
        _tb += _b; _to += _oo; _th += _ho
        _md += (f"| {period_label(_p)} | {_p} | {_b} | {_oo} | {_ho} | "
                f"{_ot} | {_ht} | {_j:.2f} |\n")
    _u = _tb + _to + _th
    _md += (f"| **all** | | **{_tb}** | **{_to}** | **{_th}** | "
            f"**{_tb + _to}** | **{_tb + _th}** | **{_tb / _u if _u else 0:.2f}** |\n")
    mo.md(_md)
    return


@app.cell
def _(mo):
    mo.md("""
    ---
    ## 3. Genotype similarity on common loci
    For loci both tools call length-variable, compare each sample's
    length-genotype. A cell counts only when **both** tools call it.
    """)
    return


@app.cell
def _(both_var, defaultdict, hip_loci, mo, ours_loci, run_sel, shared_samples):
    mo.stop(run_sel is None, mo.md(""))
    # Per-cell concordance over the common length-variable loci.
    cc_total = cc_conc = 0
    per_period_cc = defaultdict(lambda: [0, 0])  # period -> [conc, total]
    diff_hist = defaultdict(int)                  # sum(ours)-sum(hip) -> count
    per_locus_conc = []                           # locus concordance fraction
    for _k in both_var:
        _orec, _hrec = ours_loci[_k], hip_loci[_k]
        _p = _orec["period"]
        _lc = _lt = 0
        for _s in shared_samples:
            _og, _hg = _orec["rel"].get(_s), _hrec["rel"].get(_s)
            if _og is None or _hg is None:
                continue
            cc_total += 1
            _lt += 1
            per_period_cc[_p][1] += 1
            if _og == _hg:
                cc_conc += 1
                _lc += 1
                per_period_cc[_p][0] += 1
            else:
                diff_hist[sum(_og) - sum(_hg)] += 1
        if _lt:
            per_locus_conc.append(_lc / _lt)
    overall = cc_conc / cc_total if cc_total else 0.0
    mo.md(
        f"### Overall: **{cc_conc:,} / {cc_total:,} = {overall:.1%}** "
        f"concordant cells across {len(both_var):,} common length-variable loci"
    )
    return diff_hist, per_locus_conc, per_period_cc


@app.cell
def _(mo, per_period_cc, period_label, plt, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    _ps = sorted(p for p in per_period_cc if p is not None)
    _rates = [per_period_cc[p][0] / per_period_cc[p][1] if per_period_cc[p][1] else 0
              for p in _ps]
    _ns = [per_period_cc[p][1] for p in _ps]
    _fig, _ax = plt.subplots(figsize=(9, 4.4))
    _bars = _ax.bar(range(len(_ps)), _rates, color="#6a3d9a", width=0.6)
    for _b, _n in zip(_bars, _ns):
        _ax.annotate(f"n={_n:,}", (_b.get_x() + _b.get_width() / 2, _b.get_height()),
                     ha="center", va="bottom", fontsize=8,
                     xytext=(0, 1), textcoords="offset points")
    _ax.set_xticks(range(len(_ps)))
    _ax.set_xticklabels([f"{period_label(p)}\n(p{p})" for p in _ps])
    _ax.set_ylim(0, 1.08)
    _ax.set_ylabel("concordance")
    _ax.set_title("Genotype concordance by motif length (n = both-called cells)")
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(diff_hist, mo, plt, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    mo.stop(not diff_hist, mo.md("_No discordant cells._"))
    # Discordance: Σ(ours allele bp) − Σ(hipstr allele bp). Clip extremes to ±20.
    _clip = 20
    _binned = {}
    for _d, _c in diff_hist.items():
        _b = max(-_clip, min(_clip, _d))
        _binned[_b] = _binned.get(_b, 0) + _c
    _xs = sorted(_binned)
    _fig, _ax = plt.subplots(figsize=(10, 4.4))
    _ax.bar(_xs, [_binned[x] for x in _xs], color="#d62728", width=0.9)
    _ax.axvline(0, color="black", lw=0.6)
    _ax.set_xlabel("Σ allele-bp difference  (ours − HipSTR;  ±20 clipped)")
    _ax.set_ylabel("discordant cells")
    _ax.set_title("Where ours and HipSTR disagree (peaks at whole-motif-unit offsets)")
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(mo, per_locus_conc, plt, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    mo.stop(not per_locus_conc, mo.md(""))
    _fig, _ax = plt.subplots(figsize=(9, 4.2))
    _ax.hist(per_locus_conc, bins=20, range=(0, 1), color="#6a3d9a", alpha=0.8)
    _ax.set_xlabel("per-locus concordance (fraction of both-called samples agreeing)")
    _ax.set_ylabel("loci")
    _n_perfect = sum(1 for v in per_locus_conc if v == 1.0)
    _ax.set_title(
        f"Per-locus concordance — {_n_perfect:,}/{len(per_locus_conc):,} "
        f"loci fully concordant"
    )
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(mo):
    mo.md("""
    ---
    ## 4. Other relevant stats
    Call rate, heterozygosity (relevant to the high F_IS), and allelic
    richness — computed on the common length-variable loci.
    """)
    return


@app.cell
def _(
    HIP_C,
    OURS_C,
    both_var,
    hip_loci,
    mo,
    ours_loci,
    plt,
    run_sel,
    shared_samples,
):
    mo.stop(run_sel is None, mo.md(""))

    def _tool_stats(loci):
        called = cells = het = 0
        for _k in both_var:
            _rec = loci[_k]
            for _s in shared_samples:
                cells += 1
                _g = _rec["rel"].get(_s)
                if _g is None:
                    continue
                called += 1
                if len(_g) == 2 and _g[0] != _g[1]:
                    het += 1
        return {
            "call_rate": called / cells if cells else 0.0,
            "het_rate": het / called if called else 0.0,
        }

    _o, _h = _tool_stats(ours_loci), _tool_stats(hip_loci)
    _fig, _axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for _ax, _metric, _title in (
        (_axes[0], "call_rate", "call rate (cells called / cells)"),
        (_axes[1], "het_rate", "heterozygosity (het / called)"),
    ):
        _vals = [_o[_metric], _h[_metric]]
        _bars = _ax.bar(["ours", "HipSTR"], _vals, color=[OURS_C, HIP_C], width=0.6)
        for _b in _bars:
            _ax.annotate(f"{_b.get_height():.1%}",
                         (_b.get_x() + _b.get_width() / 2, _b.get_height()),
                         ha="center", va="bottom", fontsize=9,
                         xytext=(0, 1), textcoords="offset points")
        _ax.set_ylim(0, 1.08)
        _ax.set_title(_title)
        _ax.grid(True, axis="y", alpha=0.3)
    _fig.suptitle("Per-tool call rate + heterozygosity on common length-variable loci",
                  y=1.02)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(
    HIP_C,
    OURS_C,
    both_var,
    hip_loci,
    locus_alleles,
    mo,
    ours_loci,
    plt,
    run_sel,
):
    mo.stop(run_sel is None, mo.md(""))
    # Allelic richness: distinct rel-allele lengths per common locus, per tool.
    _o_rich = [len(locus_alleles(ours_loci[k])) for k in both_var]
    _h_rich = [len(locus_alleles(hip_loci[k])) for k in both_var]
    _hi = max(max(_o_rich, default=1), max(_h_rich, default=1))
    _bins = range(1, _hi + 2)
    _fig, _ax = plt.subplots(figsize=(9, 4.4))
    _ax.hist([_o_rich, _h_rich], bins=_bins, align="left",
             color=[OURS_C, HIP_C], label=["ours", "HipSTR"])
    _ax.set_xlabel("distinct allele lengths at a locus (incl. REF)")
    _ax.set_ylabel("loci")
    _ax.set_title("Allelic richness per common length-variable locus")
    _ax.legend()
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(both_var, hip_var, mo, ours_var, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    _u = len(ours_var | hip_var)
    _md = "### Summary\n\n| metric | value |\n|---|---:|\n"
    _md += f"| ours length-variable loci | {len(ours_var):,} |\n"
    _md += f"| HipSTR length-variable loci | {len(hip_var):,} |\n"
    _md += f"| shared (both) | {len(both_var):,} |\n"
    _md += f"| ours only | {len(ours_var - hip_var):,} |\n"
    _md += f"| HipSTR only | {len(hip_var - ours_var):,} |\n"
    _md += f"| Jaccard (both / union) | {len(both_var) / _u if _u else 0:.3f} |\n"
    mo.md(_md)
    return


if __name__ == "__main__":
    app.run()
