# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "matplotlib",
# ]
# ///

# Marimo dashboard — ssr_tomato1: our SSR caller vs HipSTR.
#
# Answers three questions on the ssr_tomato1 benchmark, in order:
#
#   1. LOCUS MATCH   — do the two callers emit the same loci? Both place a
#                      locus at POS = catalog_start + 1 (the HipSTR --regions
#                      BED is projected straight from our ssr-catalog), so a
#                      locus JOINS on exact (CHROM, POS). No fuzzy matching is
#                      needed or wanted: a near-but-not-equal POS means a
#                      genuinely different tract, not the same one shifted.
#
#   2. GENOTYPE MATCH — on loci BOTH callers PASS-call, does each sample's
#                      genotype agree? The two tools define an allele
#                      differently (we call repeat LENGTH / copy number,
#                      HipSTR calls a sequence haplotype), so a genotype is
#                      reduced to the sorted per-allele bp-difference-from-REF
#                      (== HipSTR's GB). An in-repeat SNP — a HipSTR allele at
#                      REF length — collapses to 0 and matches our length view,
#                      which is the apples-to-apples comparison.
#
#   3. WHY loci diverge — a HipSTR locus we don't agree on falls into one of:
#                      ours never emitted the row at all / emitted it as a
#                      no-call (lowDepth) / emitted it but filtered
#                      (notPeriodic) / emitted+PASS but the genotypes differ.
#                      The split matters: on this bench the dominant term is
#                      "ours never emitted the row", i.e. an EMISSION gap, not
#                      a genotyping disagreement.
#
# Reads the two cohort VCFs a benchmark run produced. The ours VCF lives under
# ours/cohort*/ (cohort/ for most runs, cohort53/ for the 63-sample run where
# 10 samples failed pileup); HipSTR's is always hipstr/cohort.str.vcf.gz.
#
#   uvx marimo edit --sandbox benchmarks/ssr_tomato1/scripts/ssr_vs_hipstr_dashboard.py
#   uvx marimo run  --sandbox benchmarks/ssr_tomato1/scripts/ssr_vs_hipstr_dashboard.py

import marimo

__generated_with = "0.23.10"
app = marimo.App(width="medium")


@app.cell
def _():
    import bisect
    import gzip
    from collections import Counter, defaultdict
    from pathlib import Path

    import marimo as mo

    return Counter, Path, bisect, defaultdict, gzip, mo


@app.cell
def _(mo):
    mo.md(
        """
        # ssr_tomato1 — our SSR caller vs HipSTR

        No STR truth set exists for tomato, so this scores **agreement**, not
        accuracy. Three questions, answered in order:

        1. **Locus match** — are we calling the same loci? (join on exact
           `(CHROM, POS)`; both tools place a locus at `catalog_start + 1`)
        2. **Genotype match** — on loci we both PASS-call, do the per-sample
           genotypes agree? (allele = sorted bp-difference-from-REF)
        3. **Why loci diverge** — emission gap vs filter vs real disagreement.
        """
    )
    return


@app.cell
def _(Path):
    # scripts/ssr_vs_hipstr_dashboard.py -> scripts -> ssr_tomato1
    _bench = Path(__file__).resolve().parent.parent
    # Every results* dir that has both an ours and a HipSTR cohort VCF.
    pairs = {}
    for _rd in sorted(_bench.glob("results*")):
        _hip = _rd / "hipstr" / "cohort.str.vcf.gz"
        _ours = next(iter(sorted(_rd.glob("ours/cohort*/cohort.ssr.vcf"))), None)
        if _ours is not None and _ours.exists() and _hip.exists():
            pairs[_rd.name] = (_ours, _hip)
    return (pairs,)


@app.cell
def _(mo, pairs):
    if not pairs:
        run_sel = None
        _view = mo.md(
            "_No `results*/` dir with both `ours/cohort*/cohort.ssr.vcf` and "
            "`hipstr/cohort.str.vcf.gz` found. Run the benchmark first._"
        )
    else:
        _preferred = ("results_rerun_20260708", "results_ssr15k")
        _default = next((p for p in _preferred if p in pairs), next(iter(pairs)))
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
        """-> (samples, {(chrom,pos): rec}).

        rec = {period, filter, rel}. `filter` is the VCF FILTER field (HipSTR
        writes '.', so it reads as PASS here). `rel` maps sample -> the sorted
        tuple of (len(allele) - len(REF)) over its GT, or None if the cell is
        a no-call.
        """
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
                chrom, pos, ref, alt, flt, info, fmt = (
                    f[0], int(f[1]), f[3], f[4], f[6], f[7], f[8]
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
                flag = "PASS" if flt in (".", "PASS", "") else flt
                loci[(chrom, pos)] = {
                    "period": period, "filter": flag, "rel": rel,
                    "ref": ref, "reflen": rlen,
                }
        return samples, loci

    def is_called(rec):
        """Ours PASS-called AND at least one sample has a genotype."""
        if rec["filter"] != "PASS":
            return False
        return any(g is not None for g in rec["rel"].values())

    def locus_is_length_variable(rec):
        for g in rec["rel"].values():
            if g is not None and any(a != 0 for a in g):
                return True
        return False

    def locus_alleles(rec):
        s = {0}
        for g in rec["rel"].values():
            if g is not None:
                s.update(g)
        return s

    return (
        is_called,
        locus_alleles,
        locus_is_length_variable,
        parse_vcf,
        period_label,
    )


@app.cell
def _(mo, pairs, parse_vcf, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    _ours_path, _hip_path = pairs[run_sel.value]
    ours_samples, ours_loci = parse_vcf(_ours_path)
    hip_samples, hip_loci = parse_vcf(_hip_path)
    shared_samples = [s for s in ours_samples if s in set(hip_samples)]
    mo.md(
        f"**{run_sel.value}** — ours: {len(ours_loci):,} rows / "
        f"{len(ours_samples)} samples · HipSTR: {len(hip_loci):,} rows / "
        f"{len(hip_samples)} samples · **{len(shared_samples)} paired samples**\n\n"
        f"_ours file:_ `{_ours_path.relative_to(_ours_path.parents[3])}` · "
        f"_hipstr file:_ `{_hip_path.relative_to(_hip_path.parents[2])}`"
    )
    return hip_loci, ours_loci, shared_samples


# ==========================================================================
# 1. LOCUS MATCH — are we calling the same loci?
# ==========================================================================
@app.cell
def _(mo):
    mo.md(
        """
        ---
        ## 1. Locus match — are we calling the same loci?

        A HipSTR locus (`FILTER='.'`, always emitted) is checked against ours
        on exact `(CHROM, POS)`. For each HipSTR locus we ask *what ours did
        with the same coordinate*:

        - **PASS-called** — ours emitted it and genotyped ≥1 sample. This is
          the set eligible for the genotype comparison in §2.
        - **filtered (notPeriodic)** — ours emitted it but rejected the motif.
        - **no-call (lowDepth)** — ours emitted the row but called no sample.
        - **not emitted** — ours produced no row for this coordinate at all.
          On this bench these are covered, in-slice loci ours *drops* rather
          than emits — an **emission gap**, not a genotyping disagreement.
        """
    )
    return


@app.cell
def _(Counter, hip_loci, is_called, mo, ours_loci, run_sel):
    mo.stop(run_sel is None, mo.md(""))

    def _ours_status(key):
        rec = ours_loci.get(key)
        if rec is None:
            return "not emitted"
        if rec["filter"] != "PASS":
            return f"filtered ({rec['filter']})"
        if is_called(rec):
            return "PASS-called"
        return "no-call (lowDepth)"

    hip_status = {k: _ours_status(k) for k in hip_loci}
    status_counts = Counter(hip_status.values())
    return hip_status, status_counts


@app.cell
def _(hip_loci, mo, plt, run_sel, status_counts):
    mo.stop(run_sel is None, mo.md(""))
    _order = ["PASS-called", "no-call (lowDepth)", "not emitted"]
    _rest = sorted(k for k in status_counts if k not in _order)
    _labels = _order + _rest
    _labels = [l for l in _labels if status_counts.get(l, 0)]
    _vals = [status_counts[l] for l in _labels]
    _colors = {"PASS-called": "#2ca02c", "no-call (lowDepth)": "#ff7f0e",
               "not emitted": "#d62728"}
    _bar_c = [_colors.get(l, "#9467bd") for l in _labels]
    _fig, _ax = plt.subplots(figsize=(9, 4.4))
    _bars = _ax.barh(range(len(_labels)), _vals, color=_bar_c)
    _tot = len(hip_loci)
    for _b, _v in zip(_bars, _vals):
        _ax.annotate(f"{_v:,}  ({_v / _tot:.0%})",
                     (_b.get_width(), _b.get_y() + _b.get_height() / 2),
                     ha="left", va="center", fontsize=9,
                     xytext=(3, 0), textcoords="offset points")
    _ax.set_yticks(range(len(_labels)))
    _ax.set_yticklabels(_labels)
    _ax.invert_yaxis()
    _ax.set_xlabel("HipSTR loci")
    _ax.set_xlim(0, max(_vals) * 1.18)
    _ax.set_title(f"What ours did with each of HipSTR's {_tot:,} loci")
    _ax.grid(True, axis="x", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(hip_loci, mo, ours_loci, run_sel, status_counts):
    mo.stop(run_sel is None, mo.md(""))
    _n_hip = len(hip_loci)
    _n_ours = len(ours_loci)
    _shared = len(set(ours_loci) & set(hip_loci))
    _both_call = status_counts.get("PASS-called", 0)
    _md = "### Locus overlap\n\n| | count |\n|---|---:|\n"
    _md += f"| HipSTR loci emitted | {_n_hip:,} |\n"
    _md += f"| ours rows emitted (incl. no-calls) | {_n_ours:,} |\n"
    _md += f"| shared coordinates (any ours status) | {_shared:,} |\n"
    _md += f"| **both PASS-called (→ §2 compares these)** | **{_both_call:,}** |\n"
    mo.md(_md)
    return


# ==========================================================================
# 1b. LOSS DIAGNOSIS — WHY we lose the loci we care about
# ==========================================================================
@app.cell
def _(mo):
    mo.md(
        """
        ---
        ## 1b. Loss diagnosis — why do we lose loci?

        The counts above join on exact `(CHROM, POS)`, but that *overcounts*
        loss: when a repeat has an impure prefix (`TGCTGC` before a `(TAC)ₙ`
        tract), HipSTR anchors on it while our catalog trims to the pure core,
        so the same locus sits a few bp apart and the exact join misses it.
        This section joins on **tract overlap** instead, and restricts to the
        loci that actually matter — the ones HipSTR calls **length-variable**
        (some sample carries a whole-motif length change). Monomorphic loci
        are set aside: losing them costs nothing.

        Each HipSTR length-variable locus lands in exactly one bucket:

        - **called (exact)** — ours has a row at the same POS.
        - **called (boundary-shifted)** — ours has an overlapping PASS tract at
          a slightly different POS (the impure-prefix case; a join artifact,
          not a real loss, but worth normalizing).
        - **filtered / no-call** — ours emitted a row but filtered it
          (notPeriodic) or called no sample (lowDepth).
        - **genuine drop** — ours emits no overlapping row at all. These are
          the real losses to fix. They are broken down by motif length and by
          whether the tract carries an interior interruption.
        """
    )
    return


@app.cell
def _(bisect):
    def build_pass_index(loci):
        """chrom -> sorted [(start, end)] over PASS loci with a called sample."""
        idx = {}
        for (ch, pos), rec in loci.items():
            if rec["filter"] != "PASS":
                continue
            if not any(g is not None for g in rec["rel"].values()):
                continue
            idx.setdefault(ch, []).append((pos, pos + rec["reflen"]))
        for ch in idx:
            idx[ch].sort()
        return idx

    def overlaps(idx, ch, s, e):
        v = idx.get(ch)
        if not v:
            return False
        starts = [x[0] for x in v]
        j = bisect.bisect_right(starts, e) - 1
        while j >= 0 and v[j][0] >= s - 500:
            if v[j][0] < e and v[j][1] > s:
                return True
            j -= 1
        return False

    def has_interruption(ref, period):
        """A base in the tract disagrees with the first-frame motif tiling."""
        if not period or period < 1 or len(ref) < 2 * period:
            return False
        motif = ref[:period]
        return any(b != motif[i % period] for i, b in enumerate(ref))

    def is_hip_length_variable(rec):
        p = rec["period"]
        if not p:
            return False
        for g in rec["rel"].values():
            if g is not None and any(d != 0 and d % p == 0 for d in g):
                return True
        return False

    return build_pass_index, has_interruption, is_hip_length_variable, overlaps


@app.cell
def _(
    Counter,
    build_pass_index,
    has_interruption,
    hip_loci,
    is_hip_length_variable,
    mo,
    ours_loci,
    overlaps,
    run_sel,
):
    mo.stop(run_sel is None, mo.md(""))
    _pass_idx = build_pass_index(ours_loci)
    buckets = Counter()
    drop_period = Counter()
    drop_interrupted = 0
    for _k, _hrec in hip_loci.items():
        if not is_hip_length_variable(_hrec):
            continue
        _ch, _pos = _k
        _orec = ours_loci.get(_k)
        if _orec is not None and _orec["filter"] == "PASS" \
                and any(g is not None for g in _orec["rel"].values()):
            buckets["called (exact)"] += 1
        elif _orec is not None and _orec["filter"] != "PASS":
            buckets[f"filtered ({_orec['filter']})"] += 1
        elif _orec is not None:
            buckets["no-call (lowDepth row)"] += 1
        elif overlaps(_pass_idx, _ch, _pos, _pos + _hrec["reflen"]):
            buckets["called (boundary-shifted)"] += 1
        else:
            buckets["genuine drop"] += 1
            drop_period[_hrec["period"]] += 1
            if has_interruption(_hrec["ref"], _hrec["period"]):
                drop_interrupted += 1
    n_var = sum(buckets.values())
    return buckets, drop_interrupted, drop_period, n_var


@app.cell
def _(buckets, mo, n_var, plt, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    mo.stop(not n_var, mo.md("_HipSTR called no length-variable loci here._"))
    _order = ["called (exact)", "called (boundary-shifted)",
              "no-call (lowDepth row)", "genuine drop"]
    _rest = sorted(k for k in buckets if k not in _order)
    _labels = [l for l in _order + _rest if buckets.get(l, 0)]
    _vals = [buckets[l] for l in _labels]
    _cmap = {"called (exact)": "#2ca02c", "called (boundary-shifted)": "#98df8a",
             "no-call (lowDepth row)": "#ff7f0e", "genuine drop": "#d62728"}
    _bc = [_cmap.get(l, "#9467bd") for l in _labels]
    _fig, _ax = plt.subplots(figsize=(9, 4.2))
    _bars = _ax.barh(range(len(_labels)), _vals, color=_bc)
    for _b, _v in zip(_bars, _vals):
        _ax.annotate(f"{_v:,}  ({_v / n_var:.0%})",
                     (_b.get_width(), _b.get_y() + _b.get_height() / 2),
                     ha="left", va="center", fontsize=9,
                     xytext=(3, 0), textcoords="offset points")
    _ax.set_yticks(range(len(_labels)))
    _ax.set_yticklabels(_labels)
    _ax.invert_yaxis()
    _ax.set_xlim(0, max(_vals) * 1.2)
    _ax.set_xlabel("HipSTR length-variable loci")
    _ax.set_title(f"Fate of HipSTR's {n_var:,} length-variable loci in ours")
    _ax.grid(True, axis="x", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(drop_interrupted, drop_period, mo, period_label, plt, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    _tot = sum(drop_period.values())
    mo.stop(not _tot, mo.md("_No genuine drops._"))
    _ps = sorted(p for p in drop_period if p is not None)
    _fig, _ax = plt.subplots(figsize=(9, 4.2))
    _bars = _ax.bar(range(len(_ps)), [drop_period[p] for p in _ps],
                    color="#d62728", width=0.6)
    for _b, _p in zip(_bars, _ps):
        _ax.annotate(f"{drop_period[_p]:,}",
                     (_b.get_x() + _b.get_width() / 2, _b.get_height()),
                     ha="center", va="bottom", fontsize=9,
                     xytext=(0, 1), textcoords="offset points")
    _ax.set_xticks(range(len(_ps)))
    _ax.set_xticklabels([f"{period_label(p)}\n(p{p})" for p in _ps])
    _ax.set_ylabel("genuine-drop loci")
    _ax.set_title(
        f"Genuine drops by motif length — {_tot:,} loci, "
        f"{drop_interrupted:,} ({drop_interrupted / _tot:.0%}) carry an "
        f"interior interruption"
    )
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


# ==========================================================================
# 2. GENOTYPE MATCH — on common PASS-called loci
# ==========================================================================
@app.cell
def _(mo):
    mo.md(
        """
        ---
        ## 2. Genotype match — on loci we both PASS-call

        Restricted to loci ours **PASS-called** and HipSTR emitted. A sample
        cell counts only when *both* tools genotyped it; the genotypes agree
        iff their sorted bp-difference-from-REF multisets are equal.
        """
    )
    return


@app.cell
def _(defaultdict, hip_loci, is_called, mo, ours_loci, run_sel, shared_samples):
    mo.stop(run_sel is None, mo.md(""))
    common = [k for k in hip_loci if is_called(ours_loci.get(k, {"filter": "x", "rel": {}}))]
    cc_total = cc_conc = 0
    per_period_cc = defaultdict(lambda: [0, 0])   # period -> [conc, total]
    diff_hist = defaultdict(int)                   # sum(ours)-sum(hip) -> count
    per_locus_conc = []
    for _k in common:
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
    _overall = cc_conc / cc_total if cc_total else 0.0
    mo.md(
        f"### Overall: **{cc_conc:,} / {cc_total:,} = {_overall:.1%}** "
        f"concordant cells across {len(common):,} common PASS-called loci"
    )
    return common, diff_hist, per_locus_conc, per_period_cc


@app.cell
def _(mo, per_period_cc, period_label, plt, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    _ps = sorted(p for p in per_period_cc if p is not None)
    mo.stop(not _ps, mo.md("_No common PASS-called loci to compare._"))
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
    _clip = 20
    _binned = {}
    for _d, _c in diff_hist.items():
        _b = max(-_clip, min(_clip, _d))
        _binned[_b] = _binned.get(_b, 0) + _c
    _xs = sorted(_binned)
    _fig, _ax = plt.subplots(figsize=(10, 4.2))
    _ax.bar(_xs, [_binned[x] for x in _xs], color="#d62728", width=0.9)
    _ax.axvline(0, color="black", lw=0.6)
    _ax.set_xlabel("Σ allele-bp difference  (ours − HipSTR;  ±20 clipped)")
    _ax.set_ylabel("discordant cells")
    _ax.set_title("Where ours and HipSTR disagree (whole-motif offsets = stutter/slip)")
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(mo, per_locus_conc, plt, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    mo.stop(not per_locus_conc, mo.md(""))
    _fig, _ax = plt.subplots(figsize=(9, 4.0))
    _ax.hist(per_locus_conc, bins=20, range=(0, 1), color="#6a3d9a", alpha=0.85)
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


# ==========================================================================
# 3. Call rate / heterozygosity / allelic richness on common loci
# ==========================================================================
@app.cell
def _(mo):
    mo.md(
        """
        ---
        ## 3. Per-tool stats on common PASS-called loci
        Call rate, heterozygosity (relevant to the high reported F_IS), and
        allelic richness — over the loci both tools call.
        """
    )
    return


@app.cell
def _():
    import matplotlib.pyplot as plt
    OURS_C, HIP_C = "#1f77b4", "#d62728"
    return HIP_C, OURS_C, plt


@app.cell
def _(HIP_C, OURS_C, common, hip_loci, mo, ours_loci, plt, run_sel, shared_samples):
    mo.stop(run_sel is None, mo.md(""))
    mo.stop(not common, mo.md("_No common PASS-called loci._"))

    def _tool_stats(loci):
        called = cells = het = 0
        for _k in common:
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
    _fig, _axes = plt.subplots(1, 2, figsize=(10, 4.0))
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
    _fig.suptitle("Per-tool call rate + heterozygosity on common PASS-called loci",
                  y=1.02)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(HIP_C, OURS_C, common, hip_loci, locus_alleles, mo, ours_loci, plt, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    mo.stop(not common, mo.md(""))
    _o_rich = [len(locus_alleles(ours_loci[k])) for k in common]
    _h_rich = [len(locus_alleles(hip_loci[k])) for k in common]
    _hi = max(max(_o_rich, default=1), max(_h_rich, default=1))
    _bins = range(1, _hi + 2)
    _fig, _ax = plt.subplots(figsize=(9, 4.2))
    _ax.hist([_o_rich, _h_rich], bins=_bins, align="left",
             color=[OURS_C, HIP_C], label=["ours", "HipSTR"])
    _ax.set_xlabel("distinct allele lengths at a locus (incl. REF)")
    _ax.set_ylabel("loci")
    _ax.set_title("Allelic richness per common PASS-called locus")
    _ax.legend()
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


if __name__ == "__main__":
    app.run()
