# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "matplotlib",
# ]
# ///

# Marimo dashboard — ssr_tomato1: estimating false calls WITHOUT a truth set.
#
# There is no gold-standard microsatellite truth for tomato, so we cannot score
# accuracy directly. Instead this dashboard leans on two things we DO know:
#
#   IDEA 1 — the selfing ruler.  Tomato self-pollinates, so the plants are
#     highly inbred and almost every TRUE genotype should be homozygous (both
#     chromosomes the same repeat length). Heterozygous calls are therefore
#     mostly ERRORS. So: (a) a caller's heterozygous-call rate is a
#     false-positive meter, and (b) we can split each disagreement by the SHAPE
#     of the reported variation — a length difference seen as a homozygous,
#     recurring allele looks real (a fixed difference between inbred lines),
#     while a lone heterozygote at a slippage-prone repeat looks fake.
#
#   IDEA 3 — the error fingerprint.  A real heterozygote is supported by roughly
#     half its reads for each allele; a false one (PCR/sequencing slippage,
#     "stutter") is lopsided — a big pile at one length and a thin tail at a
#     neighbouring length one motif unit away. HipSTR reports per-read support
#     (PDP / AB / MALLREADS), so we can measure that lopsidedness directly and
#     estimate what fraction of its heterozygous calls are stutter artifacts.
#
# Our caller emits only GT:GQ:REPCN (no per-read support), so idea 3's
# read-balance view is shown for HipSTR — which is where the "are we missing
# real variants or is HipSTR over-calling?" question actually lives.
#
#   uvx marimo edit --sandbox benchmarks/ssr_tomato1/scripts/ssr_error_signals_dashboard.py
#   uvx marimo run  --sandbox benchmarks/ssr_tomato1/scripts/ssr_error_signals_dashboard.py

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
        # ssr_tomato1 — estimating false calls without a truth set

        No microsatellite truth exists for tomato, so we cannot score accuracy
        directly. This dashboard uses two things we *do* know — the **selfing
        biology** (idea 1) and the **read-support fingerprint of slippage**
        (idea 3) — to estimate where the false positives and false negatives
        are, comparing our caller against HipSTR.
        """
    )
    return


# ── run discovery + dropdown ──────────────────────────────────────────────
@app.cell
def _(Path):
    _bench = Path(__file__).resolve().parent.parent
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
            "_No `results*/` dir with both an ours VCF and a HipSTR VCF found._"
        )
    else:
        _pref = ("results_rerun_20260708", "results_ssr15k")
        _default = next((p for p in _pref if p in pairs), next(iter(pairs)))
        run_sel = mo.ui.dropdown(
            options=list(pairs.keys()), value=_default, label="benchmark run"
        )
        _view = run_sel
    _view
    return (run_sel,)


# ── parsing ───────────────────────────────────────────────────────────────
@app.cell
def _(gzip):
    PERIOD_NAME = {2: "di (2 bp)", 3: "tri (3 bp)", 4: "tetra (4 bp)",
                   5: "penta (5 bp)", 6: "hexa (6 bp)"}

    def period_label(p):
        return PERIOD_NAME.get(p, f"{p} bp" if p else "?")

    def open_vcf(path):
        p = str(path)
        return gzip.open(p, "rt") if p.endswith((".gz", ".bgz")) else open(p)

    def info_int(info, key):
        for kv in info.split(";"):
            if kv.startswith(key + "="):
                try:
                    return int(kv.split("=")[1])
                except ValueError:
                    return None
        return None

    def gt_idx(gt):
        parts = gt.replace("|", "/").split("/")
        if any(p in (".", "") for p in parts):
            return None
        return tuple(int(p) for p in parts)

    def to_float(s):
        try:
            return float(s)
        except (ValueError, TypeError):
            return None

    def minor_frac_pdp(pdp):
        """PDP = 'a|b' weighted reads per haploid allele → min/(a+b)."""
        try:
            a, b = (float(x) for x in pdp.split("|"))
        except (ValueError, AttributeError):
            return None
        tot = a + b
        return min(a, b) / tot if tot > 0 else None

    def parse_mallreads(s):
        """MALLREADS = 'bpdiff|count;bpdiff|count;…' → {bpdiff: count}."""
        d = {}
        if s in (".", "", None):
            return d
        for tok in s.split(";"):
            if "|" not in tok:
                continue
            bp, cnt = tok.split("|")
            try:
                d[int(bp)] = d.get(int(bp), 0) + int(cnt)
            except ValueError:
                pass
        return d

    return (PERIOD_NAME, gt_idx, info_int, minor_frac_pdp, open_vcf,
            parse_mallreads, period_label, to_float)


@app.cell
def _(gt_idx, info_int, minor_frac_pdp, open_vcf, parse_mallreads, to_float):
    def parse_ours(path):
        """(chrom,pos) -> {period, filter, reflen, cells:{sample: cell}}.

        cell = {gt, bpdiff, het, gq, dp, minor_frac}. minor_frac (het only) is
        the minority-allele read fraction from the AD field (per-allele read
        depth) — the same balance signal HipSTR exposes via PDP, now that
        ssr-call emits AD. bp-diff = len(allele)-len(REF).
        """
        samples, loci = [], {}
        for line in open_vcf(path):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, flt, info, fmt = (
                f[0], int(f[1]), f[3], f[4], f[6], f[7], f[8])
            alleles = [ref] + ([] if alt == "." else alt.split(","))
            rlen = len(ref)
            keys = fmt.split(":")
            gi = keys.index("GT")
            qi = keys.index("GQ") if "GQ" in keys else None
            ai = keys.index("AD") if "AD" in keys else None
            di = keys.index("DP") if "DP" in keys else None
            cells = {}
            for name, col in zip(samples, f[9:]):
                parts = col.split(":")
                idx = gt_idx(parts[gi])
                if idx is None:
                    cells[name] = None
                    continue
                bp = tuple(sorted(len(alleles[i]) - rlen for i in idx))
                het = len(set(idx)) > 1
                gq = None
                if qi is not None and qi < len(parts):
                    try:
                        gq = int(parts[qi])
                    except ValueError:
                        gq = None
                dp = None
                if di is not None and di < len(parts):
                    try:
                        dp = int(parts[di])
                    except ValueError:
                        dp = None
                # minor-allele read fraction from AD over the two called alleles.
                minor_frac = None
                if het and ai is not None and ai < len(parts) and parts[ai] != ".":
                    try:
                        ad = [int(x) for x in parts[ai].split(",")]
                        a, b = ad[idx[0]], ad[idx[1]]
                        minor_frac = min(a, b) / (a + b) if (a + b) > 0 else None
                    except (ValueError, IndexError):
                        minor_frac = None
                cells[name] = {"gt": idx, "bpdiff": bp, "het": het, "gq": gq,
                               "dp": dp, "minor_frac": minor_frac}
            loci[(chrom, pos)] = {
                "period": info_int(info, "PERIOD"),
                "filter": ("PASS" if flt in (".", "PASS", "") else flt),
                "reflen": rlen, "cells": cells}
        return samples, loci

    def parse_hipstr(path):
        """(chrom,pos) -> {period, reflen, cells:{sample: cell}}.

        cell adds the read-support fingerprint: minor_frac (from PDP), ab (AB
        log10 bias p-value), dp, mallreads dict, and the two alleles' bp-diffs
        from GB.
        """
        samples, loci = [], {}
        for line in open_vcf(path):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, info, fmt = (
                f[0], int(f[1]), f[3], f[4], f[7], f[8])
            rlen = len(ref)
            keys = fmt.split(":")
            gi = keys.index("GT")
            idx_of = {k: keys.index(k) for k in
                      ("GB", "PDP", "AB", "DP", "MALLREADS") if k in keys}
            cells = {}
            for name, col in zip(samples, f[9:]):
                parts = col.split(":")
                idx = gt_idx(parts[gi])
                if idx is None:
                    cells[name] = None
                    continue

                def get(k):
                    j = idx_of.get(k)
                    return parts[j] if j is not None and j < len(parts) else None

                gb = get("GB")
                bp = None
                if gb and gb not in (".", ""):
                    try:
                        bp = tuple(sorted(int(x) for x in gb.replace("|", "/").split("/")))
                    except ValueError:
                        bp = None
                het = len(set(idx)) > 1
                cells[name] = {
                    "gt": idx, "bpdiff": bp, "het": het,
                    "minor_frac": minor_frac_pdp(get("PDP")) if het else None,
                    "ab": to_float(get("AB")) if het else None,
                    "dp": to_float(get("DP")),
                    "mallreads": parse_mallreads(get("MALLREADS")),
                }
            loci[(chrom, pos)] = {
                "period": info_int(info, "PERIOD"), "reflen": rlen,
                "cells": cells}
        return samples, loci

    return parse_hipstr, parse_ours


# ── shared locus classifiers ──────────────────────────────────────────────
@app.cell
def _():
    def variation_shape(loc):
        """How a locus's WHOLE-UNIT LENGTH variation is expressed across the
        cohort — the comparable set for our length (REPCN) model. HipSTR's
        sequence-only and non-whole-unit changes are outside that model, so they
        are not counted here (we are silent on them by design, not wrong).

        Returns (n_length_variant_samples, n_homalt_samples, n_het_samples). A
        sample is counted if any called allele differs from the reference by a
        whole number of motif units; 'hom-alt' if it is homozygous for such an
        allele; 'het' if its two alleles differ.
        """
        period = loc["period"]
        n_var = n_homalt = n_het = 0
        for cell in loc["cells"].values():
            if cell is None or not cell["bpdiff"]:
                continue
            whole_unit = [d for d in cell["bpdiff"]
                          if d != 0 and period and d % period == 0]
            if not whole_unit:
                continue
            n_var += 1
            if cell["het"]:
                n_het += 1
            elif len(set(cell["bpdiff"])) == 1:
                n_homalt += 1
        return n_var, n_homalt, n_het

    def classify_credibility(loc):
        """Split a variant locus into 'credible' vs 'suspect' by its shape,
        using the selfer expectation (real = a homozygous, ideally recurring
        allele; suspect = het-only, especially a singleton)."""
        n_var, n_homalt, n_het = variation_shape(loc)
        if n_var == 0:
            return "monomorphic"
        if n_homalt >= 1:
            return "credible (hom-alt present)"
        if n_var == 1:
            return "suspect (singleton het)"
        return "suspect (het-only, multi-sample)"

    return classify_credibility, variation_shape


@app.cell
def _(bisect):
    def build_pass_index(ours_loci):
        idx = {}
        for (ch, pos), rec in ours_loci.items():
            if rec["filter"] != "PASS":
                continue
            if not any(c is not None for c in rec["cells"].values()):
                continue
            idx.setdefault(ch, []).append((pos, pos + rec["reflen"]))
        for ch in idx:
            idx[ch].sort()
        return idx

    def ours_covers(idx, ch, s, e):
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

    return build_pass_index, ours_covers


@app.cell
def _(mo, pairs, parse_hipstr, parse_ours, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    _op, _hp = pairs[run_sel.value]
    ours_samples, ours_loci = parse_ours(_op)
    hip_samples, hip_loci = parse_hipstr(_hp)
    shared = [s for s in ours_samples if s in set(hip_samples)]
    mo.md(
        f"**{run_sel.value}** — ours {len(ours_loci):,} rows · "
        f"HipSTR {len(hip_loci):,} rows · {len(shared)} shared samples"
    )
    return hip_loci, ours_loci, shared


# ==========================================================================
# IDEA 1a — heterozygosity as a false-positive ruler
# ==========================================================================
@app.cell
def _(mo):
    mo.md(
        """
        ---
        ## 1a. Heterozygosity as a false-positive ruler

        In a near-complete selfer almost every real genotype is homozygous, so
        heterozygous calls are **mostly errors**. A higher heterozygous-call
        rate — especially at the short, slippage-prone repeats — means more
        manufactured calls. This compares the two callers' het rates overall
        and by motif length.
        """
    )
    return


@app.cell
def _(hip_loci, ours_loci, shared):
    def het_stats(loci, ref_idx=0):
        """-> (overall_rate, {period: (het, called)}). Counts one genotype per
        called sample cell; a het is a cell whose two alleles differ."""
        shared_set = set(shared)
        per_period = {}
        het = called = 0
        for rec in loci.values():
            p = rec["period"]
            for name, cell in rec["cells"].items():
                if name not in shared_set or cell is None:
                    continue
                called += 1
                h = 1 if cell["het"] else 0
                het += h
                a, c = per_period.get(p, (0, 0))
                per_period[p] = (a + h, c + 1)
        rate = het / called if called else 0.0
        return rate, het, called, per_period

    ours_het = het_stats(ours_loci)
    hip_het = het_stats(hip_loci)
    return hip_het, ours_het


@app.cell
def _(hip_het, mo, ours_het, period_label, plt, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    _periods = [2, 3, 4, 5, 6]
    _labels = [period_label(p) for p in _periods]

    def _rates(stats):
        pp = stats[3]
        out = []
        for p in _periods:
            a, c = pp.get(p, (0, 0))
            out.append(a / c if c else 0.0)
        return out

    _ours = _rates(ours_het)
    _hip = _rates(hip_het)
    _fig, _ax = plt.subplots(figsize=(9, 4.6))
    _x = range(len(_periods))
    _w = 0.38
    _ax.bar([i - _w / 2 for i in _x], _ours, _w, label="ours", color="#2ca02c")
    _ax.bar([i + _w / 2 for i in _x], _hip, _w, label="HipSTR", color="#1f77b4")
    _ax.set_xticks(list(_x))
    _ax.set_xticklabels(_labels)
    _ax.set_ylabel("heterozygous-call rate")
    _ax.set_xlabel("motif length")
    _ax.set_title("Het rate by motif length — higher = more likely-false calls "
                  "(selfer ⇒ true het ≈ 0)")
    _ax.legend()
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(hip_het, mo, ours_het):
    mo.stop(ours_het is None, mo.md(""))
    _md = "### Overall heterozygous-call rate\n\n| caller | het cells | called cells | het rate |\n|---|--:|--:|--:|\n"
    for _name, _s in (("ours", ours_het), ("HipSTR", hip_het)):
        _md += f"| {_name} | {_s[1]:,} | {_s[2]:,} | **{_s[0]:.1%}** |\n"
    _md += ("\n_In a selfer the true het rate is near zero, so the gap between "
            "these numbers and zero is a rough false-heterozygote estimate. "
            "The caller with the higher rate is manufacturing more._")
    mo.md(_md)
    return


# ==========================================================================
# IDEA 1b — split the disagreements by genotype shape
# ==========================================================================
@app.cell
def _(mo):
    mo.md(
        """
        ---
        ## 1b. Splitting the disagreements by genotype shape

        Two disagreement sets, each split into **credible** (a homozygous
        allele — a fixed difference between inbred lines, what a real
        polymorphism looks like in a selfer) vs **suspect** (het-only,
        especially a lone singleton — what a slippage artifact looks like):

        - **HipSTR calls it variable, we emit nothing** → the loci we "miss".
          Credible ones are our likely **false negatives**; suspect ones are
          likely **HipSTR false positives** (our silence is probably right).
        - **We call it variable, HipSTR calls it monomorphic** → loci only we
          call. Suspect ones here are our likely **false positives**.
        """
    )
    return


@app.cell
def _(build_pass_index, classify_credibility, hip_loci, ours_covers,
      ours_loci, variation_shape):
    _pass_idx = build_pass_index(ours_loci)

    # Set A: HipSTR variable, ours emits no overlapping row (the "misses").
    miss = []
    for (_ch, _pos), _hrec in hip_loci.items():
        _nv, _, _ = variation_shape(_hrec)
        if _nv == 0:
            continue
        _orec = ours_loci.get((_ch, _pos))
        if _orec is not None and _orec["filter"] == "PASS" and any(
                c is not None for c in _orec["cells"].values()):
            continue
        if _orec is not None:  # filtered / no-call row exists
            continue
        if ours_covers(_pass_idx, _ch, _pos, _pos + _hrec["reflen"]):
            continue
        miss.append(((_ch, _pos), _hrec))

    # Set B: ours PASS-variable, HipSTR monomorphic (or absent) — our unique calls.
    unique = []
    for (_ch, _pos), _orec in ours_loci.items():
        if _orec["filter"] != "PASS":
            continue
        _nv, _, _ = variation_shape(_orec)
        if _nv == 0:
            continue
        _hrec = hip_loci.get((_ch, _pos))
        if _hrec is not None:
            _hv, _, _ = variation_shape(_hrec)
            if _hv > 0:
                continue  # both call it variable — not unique
        unique.append(((_ch, _pos), _orec))

    from collections import Counter as _C
    miss_shape = _C(classify_credibility(h) for _, h in miss)
    unique_shape = _C(classify_credibility(o) for _, o in unique)
    return miss, miss_shape, unique, unique_shape


@app.cell
def _(miss, miss_shape, mo, plt, run_sel, unique, unique_shape):
    mo.stop(run_sel is None, mo.md(""))
    _order = ["credible (hom-alt present)", "suspect (het-only, multi-sample)",
              "suspect (singleton het)"]
    _cmap = {"credible (hom-alt present)": "#2ca02c",
             "suspect (het-only, multi-sample)": "#ff7f0e",
             "suspect (singleton het)": "#d62728"}

    def _panel(ax, shape, total, title):
        labs = [l for l in _order if shape.get(l, 0)]
        vals = [shape[l] for l in labs]
        bars = ax.barh(range(len(labs)), vals, color=[_cmap[l] for l in labs])
        for b, v in zip(bars, vals):
            ax.annotate(f"{v}  ({v/total:.0%})",
                        (b.get_width(), b.get_y() + b.get_height() / 2),
                        ha="left", va="center", fontsize=9,
                        xytext=(3, 0), textcoords="offset points")
        ax.set_yticks(range(len(labs)))
        ax.set_yticklabels([l.replace("suspect ", "suspect\n") for l in labs])
        ax.invert_yaxis()
        ax.set_xlim(0, max(vals) * 1.25 if vals else 1)
        ax.set_title(title, fontsize=10)
        ax.grid(True, axis="x", alpha=0.3)

    _fig, _axes = plt.subplots(2, 1, figsize=(9, 7))
    _panel(_axes[0], miss_shape, max(len(miss), 1),
           f"Loci HipSTR calls but we miss (n={len(miss)}) — "
           f"credible = our likely FN, suspect = likely HipSTR FP")
    _panel(_axes[1], unique_shape, max(len(unique), 1),
           f"Loci only we call (n={len(unique)}) — suspect = our likely FP")
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(miss, miss_shape, mo, unique, unique_shape):
    mo.stop(not miss, mo.md(""))

    def _summ(shape, total):
        cred = shape.get("credible (hom-alt present)", 0)
        susp = total - cred - shape.get("monomorphic", 0)
        return cred, susp

    _mc, _ms = _summ(miss_shape, len(miss))
    _uc, _us = _summ(unique_shape, len(unique))
    mo.md(
        f"""
        ### What the split says (no truth needed)

        - **Of the {len(miss)} loci we miss:** ~**{_mc}** ({_mc/len(miss):.0%})
          are credible real polymorphisms (a homozygous allele) → our likely
          **false negatives**, worth recovering. ~**{_ms}**
          ({_ms/len(miss):.0%}) are suspect het-only/singletons → likely
          **HipSTR false positives**, where our silence is probably correct.
        - **Of the {len(unique)} loci only we call:** ~**{_us}**
          ({_us/max(len(unique),1):.0%}) are suspect → our likely **false
          positives**; ~**{_uc}** ({_uc/max(len(unique),1):.0%}) are credible.

        Idea 3 below stress-tests the "suspect HipSTR" pile directly from the
        read support.
        """
    )
    return


# ==========================================================================
# IDEA 3 — the allele-balance fingerprint of HipSTR's het calls
# ==========================================================================
@app.cell
def _(mo):
    mo.md(
        """
        ---
        ## 3. The read-support fingerprint of heterozygous calls

        A real heterozygote splits its reads ~50/50 between its two alleles; a
        slippage artifact is lopsided (a thin minority allele one motif unit
        from the majority). First: HipSTR's het calls, split by whether the
        locus is one **we also call** (agreed) or one **we stay silent on**
        (disputed) — if the disputed hets skew toward 0 they carry the slippage
        fingerprint, so our silence is defensible. Then, now that `ssr-call`
        emits `AD`, the **same test on our own het calls vs HipSTR's** — a
        caller whose hets are more balanced is calling fewer stutter artifacts.
        """
    )
    return


@app.cell
def _(build_pass_index, hip_loci, ours_covers, ours_loci):
    _pass_idx = build_pass_index(ours_loci)

    def _ours_calls_here(ch, pos, reflen):
        rec = ours_loci.get((ch, pos))
        if rec is not None and rec["filter"] == "PASS" and any(
                c is not None for c in rec["cells"].values()):
            return True
        return ours_covers(_pass_idx, ch, pos, pos + reflen)

    agreed_frac, disputed_frac = [], []
    agreed_ab, disputed_ab = [], []
    for (_ch, _pos), _hrec in hip_loci.items():
        _here = _ours_calls_here(_ch, _pos, _hrec["reflen"])
        for _cell in _hrec["cells"].values():
            if _cell is None or not _cell["het"]:
                continue
            _mf, _ab = _cell["minor_frac"], _cell["ab"]
            if _mf is not None:
                (agreed_frac if _here else disputed_frac).append(_mf)
            if _ab is not None:
                (agreed_ab if _here else disputed_ab).append(_ab)
    return agreed_ab, agreed_frac, disputed_ab, disputed_frac


@app.cell
def _(agreed_frac, disputed_frac, mo, plt, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    _fig, _ax = plt.subplots(figsize=(9, 4.6))
    _bins = [i / 20 for i in range(11)]  # 0..0.5
    _ax.hist([disputed_frac, agreed_frac], bins=_bins,
             label=[f"disputed / we're silent (n={len(disputed_frac)})",
                    f"agreed / we also call (n={len(agreed_frac)})"],
             color=["#d62728", "#2ca02c"], density=True)
    _ax.axvline(0.5, color="grey", ls="--", lw=1)
    _ax.annotate("balanced\n(real het)", (0.47, _ax.get_ylim()[1] * 0.85),
                 ha="right", fontsize=8, color="grey")
    _ax.set_xlabel("minority-allele read fraction (0 = fully lopsided, 0.5 = balanced)")
    _ax.set_ylabel("density")
    _ax.set_title("HipSTR het calls: read balance, disputed vs agreed loci")
    _ax.legend()
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(agreed_frac, disputed_frac, mo):
    mo.stop(not disputed_frac, mo.md(""))

    def _skew(vals, thr=0.20):
        n = len(vals)
        lop = sum(1 for v in vals if v < thr)
        med = sorted(vals)[n // 2]
        return n, lop, lop / n, med

    _dn, _dl, _dr, _dm = _skew(disputed_frac)
    _an, _al, _ar, _am = _skew(agreed_frac)
    mo.md(
        f"""
        ### Read-balance summary

        | HipSTR het calls | count | median minor fraction | lopsided (<0.20) |
        |---|--:|--:|--:|
        | on loci **we're silent** (disputed) | {_dn:,} | {_dm:.2f} | **{_dr:.0%}** |
        | on loci **we also call** (agreed) | {_an:,} | {_am:.2f} | {_ar:.0%} |

        A lopsided het (minority allele under ~20% of reads) carries the
        slippage fingerprint. If the disputed column is markedly more lopsided,
        a large share of the loci we "miss" are HipSTR calling stutter as a real
        allele — i.e. those are HipSTR false positives, not our false negatives.
        """
    )
    return


@app.cell
def _(hip_loci, ours_loci):
    def _het_fracs(loci, require_pass):
        out = []
        for rec in loci.values():
            if require_pass and rec["filter"] != "PASS":
                continue
            for cell in rec["cells"].values():
                if cell is None or not cell["het"]:
                    continue
                mf = cell.get("minor_frac")
                if mf is not None:
                    out.append(mf)
        return out

    ours_hetfrac = _het_fracs(ours_loci, require_pass=True)
    hip_hetfrac = _het_fracs(hip_loci, require_pass=False)
    return hip_hetfrac, ours_hetfrac


@app.cell
def _(hip_hetfrac, mo, ours_hetfrac, plt, run_sel):
    mo.stop(run_sel is None or not ours_hetfrac, mo.md(""))
    _fig, _ax = plt.subplots(figsize=(9, 4.6))
    _bins = [i / 20 for i in range(11)]
    _ax.hist([ours_hetfrac, hip_hetfrac], bins=_bins,
             label=[f"ours (n={len(ours_hetfrac):,})",
                    f"HipSTR (n={len(hip_hetfrac):,})"],
             color=["#2ca02c", "#1f77b4"], density=True)
    _ax.axvline(0.2, color="#d62728", ls="--", lw=1)
    _ax.annotate("lopsided\n(stutter)", (0.19, _ax.get_ylim()[1] * 0.85),
                 ha="right", fontsize=8, color="#d62728")
    _ax.set_xlabel("minority-allele read fraction (0 = fully lopsided, 0.5 = balanced)")
    _ax.set_ylabel("density")
    _ax.set_title("Het-call read balance: ours vs HipSTR (all PASS het calls)")
    _ax.legend()
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(hip_hetfrac, mo, ours_hetfrac):
    mo.stop(not ours_hetfrac, mo.md(""))

    def _s(v):
        n = len(v)
        return n, sum(1 for x in v if x < 0.20) / n, sorted(v)[n // 2]

    _on, _olop, _omed = _s(ours_hetfrac)
    _hn, _hlop, _hmed = _s(hip_hetfrac)
    mo.md(
        f"""
        ### Our het balance vs HipSTR's (symmetric, via `AD`)

        | caller | het calls | median minor fraction | lopsided (<0.20) |
        |---|--:|--:|--:|
        | **ours** | {_on:,} | {_omed:.2f} | **{_olop:.0%}** |
        | **HipSTR** | {_hn:,} | {_hmed:.2f} | {_hlop:.0%} |

        Fewer lopsided hets = fewer stutter artifacts called as real
        heterozygotes. This is the read-level counterpart of the precision gap
        the silver standards found in §4 — now measurable on **both** callers,
        because `ssr-call` emits `AD`.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ---
        ### Caveats

        - "Credible" and "suspect" are **priors from the biology**, not truth.
          A credible-shaped call can still be wrong, and a lone het can be a
          genuine recent outcross. These are estimates to be anchored by the
          orthogonal validation subset (idea 5), not final error rates.
        """
    )
    return


# ==========================================================================
# SILVER STANDARDS — turn the confident calls into a pseudo-truth, count FP/FN
# ==========================================================================
@app.cell
def _(mo):
    mo.md(
        """
        ---
        ## 4. Two silver standards → false positives and false negatives

        With no gold truth, we carve out the calls we're *confident* about and
        treat those as truth (a **silver standard**), accepting that the
        uncertain middle is set aside. Two definitions, built from the reads:

        - **true100** — a non-reference length that is **homozygous in ≥2
          independent plants** (a fixed difference between inbred lines;
          artifacts don't recur), or a balanced heterozygote corroborated by
          another carrier.
        - **false100** — no plant shows a clean non-reference allele; every
          apparent variant is minority stutter → the locus is really
          monomorphic.
        - **doubtful / undetermined** — everything else, and loci too shallow to
          judge. Set aside.

        We build this **twice**: once from **our own pileup** reads
        (`our_reads.tsv`, caller-neutral) and once from **HipSTR's** per-read
        fields. Then each caller is scored: a **false negative** is a true100
        locus it calls monomorphic; a **false positive** is a false100 locus it
        calls variable. The **confident core** (both standards agree) is the
        most trustworthy subset — the user's "false even under HipSTR's own
        data" test.
        """
    )
    return


@app.cell
def _(Counter, defaultdict, gzip, mo, pairs, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    _op, _ = pairs[run_sel.value]
    run_dir = _op.parents[2]
    bench_dir = run_dir.parent
    reads_path = run_dir / "our_reads.tsv"
    cat_path = next(iter(sorted(bench_dir.glob("results*/ours/*.ssr.catalog"))), None)

    reads = defaultdict(lambda: defaultdict(Counter))  # (chrom,pos)->sample->{len:count}
    have_reads = reads_path.exists()
    if have_reads:
        with open(reads_path) as _fh:
            next(_fh)
            for _line in _fh:
                _s, _ch, _st, _seq, _c = _line.rstrip("\n").split("\t")
                reads[(_ch, int(_st))][_s][len(_seq)] += int(_c)

    catalog = {}  # (chrom, pos=start+1) -> (reflen, period)
    if cat_path is not None:
        with gzip.open(cat_path, "rt") as _fh:
            for _line in _fh:
                if _line.startswith("#") or not _line.strip():
                    continue
                _p = _line.rstrip("\n").split("\t")
                if len(_p) < 4:
                    continue
                _c0, _s0, _e0, _m0 = _p[0], int(_p[1]), int(_p[2]), _p[3]
                catalog[(_c0, _s0 + 1)] = (_e0 - _s0, len(_m0))

    _msg = (f"read evidence: {len(reads):,} loci · catalog: {len(catalog):,} loci"
            if have_reads else
            "_`our_reads.tsv` not found for this run — generate it with "
            "`cargo run --release --example ssr_slip_dump -- <run>/ours/cohort/psp/*.ssr.psp "
            "> <run>/our_reads.tsv`. The HipSTR-only standard still works below._")
    mo.md(_msg)
    return catalog, have_reads, reads


@app.cell
def _(mo, run_sel):
    mo.stop(run_sel is None, mo.md(""))
    min_depth = mo.ui.slider(2, 10, value=4, label="min reads/sample")
    hom_frac = mo.ui.slider(0.5, 0.95, value=0.75, step=0.05,
                            label="homozygous modal fraction")
    bal_min, recur, min_assessable = 0.40, 2, 6  # fixed
    mo.vstack([mo.md("**Silver-standard thresholds** (recurrence ≥2 plants, "
                     "balanced-het ≥0.40, ≥6 assessable samples — fixed):"),
               min_depth, hom_frac])
    return bal_min, hom_frac, min_assessable, min_depth, recur


@app.cell
def _(Counter):
    def pileup_signals(persample, reflen, period, min_depth, hom_frac, bal_min):
        homalt, balhet, carrier, n_assess = Counter(), set(), Counter(), 0
        for lens in persample.values():
            total = sum(lens.values())
            if total < min_depth:
                continue
            wu = {L: c for L, c in lens.items()
                  if period and (L - reflen) % period == 0}
            if not wu:
                continue
            n_assess += 1
            top = sorted(wu.items(), key=lambda kv: -kv[1])
            lmod, cmod = top[0]
            dmod, fmod = lmod - reflen, cmod / total
            if dmod != 0 and fmod >= hom_frac:
                homalt[dmod] += 1
                carrier[dmod] += 1
            elif len(top) >= 2:
                l2, c2 = top[1]
                if min(cmod, c2) / (cmod + c2) >= bal_min:
                    for d in (dmod, l2 - reflen):
                        if d != 0:
                            balhet.add(d)
                            carrier[d] += 1
        return homalt, balhet, carrier, n_assess

    def hip_signals(rec, bal_min):
        period = rec["period"]
        homalt, balhet, carrier, n_assess = Counter(), set(), Counter(), 0
        for cell in rec["cells"].values():
            if cell is None or not cell["bpdiff"]:
                continue
            n_assess += 1
            bp = cell["bpdiff"]
            wu = [d for d in bp if d != 0 and period and d % period == 0]
            for d in set(wu):
                carrier[d] += 1
            if not cell["het"] and wu and len(set(bp)) == 1:
                homalt[bp[0]] += 1
            elif cell["het"] and cell.get("minor_frac") is not None \
                    and cell["minor_frac"] >= bal_min and wu:
                for d in wu:
                    balhet.add(d)
        return homalt, balhet, carrier, n_assess

    def classify_locus(homalt, balhet, carrier, n_assess, recur, min_assessable):
        if n_assess < min_assessable:
            return "undetermined"
        if any(n >= recur for n in homalt.values()) \
                or any(carrier.get(d, 0) >= 2 for d in balhet):
            return "true100"
        if not homalt and not balhet:
            return "false100"
        return "doubtful"

    return classify_locus, hip_signals, pileup_signals


@app.cell
def _(bal_min, catalog, classify_locus, have_reads, hip_loci, hip_signals,
      hom_frac, min_assessable, min_depth, mo, pileup_signals, reads, recur,
      run_sel):
    mo.stop(run_sel is None, mo.md(""))
    pileup_cls = {}
    if have_reads:
        for (_ch, _pos), _ps in reads.items():
            _meta = catalog.get((_ch, _pos))
            if _meta is None:
                continue
            _rl, _per = _meta
            _sig = pileup_signals(_ps, _rl, _per, min_depth.value,
                                  hom_frac.value, bal_min)
            pileup_cls[(_ch, _pos)] = classify_locus(*_sig, recur, min_assessable)
    hip_cls = {_k: classify_locus(*hip_signals(_rec, bal_min), recur, min_assessable)
               for _k, _rec in hip_loci.items()}
    return hip_cls, pileup_cls


@app.cell
def _(hip_loci, ours_loci, variation_shape):
    def ours_var(key):
        r = ours_loci.get(key)
        return r is not None and r["filter"] == "PASS" and variation_shape(r)[0] > 0

    def hip_var(key):
        r = hip_loci.get(key)
        return r is not None and variation_shape(r)[0] > 0

    return hip_var, ours_var


@app.cell
def _(Counter, have_reads, hip_cls, hip_var, mo, ours_var, pileup_cls, run_sel):
    mo.stop(run_sel is None, mo.md(""))

    def _score(cls, name, callers):
        out = ""
        for cn, cv in callers:
            t = [k for k, c in cls.items() if c == "true100"]
            fset = [k for k, c in cls.items() if c == "false100"]
            fn = sum(1 for k in t if not cv(k))
            fp = sum(1 for k in fset if cv(k))
            out += (f"| {name} | {cn} | {len(t):,} | {fn} ({fn/max(len(t),1):.0%}) "
                    f"| {len(fset):,} | {fp} ({fp/max(len(fset),1):.1%}) |\n")
        return out

    _callers = [("ours", ours_var), ("HipSTR", hip_var)]
    _shared = set(pileup_cls) & set(hip_cls)
    _at = [k for k in _shared if pileup_cls[k] == "true100" == hip_cls[k]]
    _af = [k for k in _shared if pileup_cls[k] == "false100" == hip_cls[k]]
    _core = {**{k: "true100" for k in _at}, **{k: "false100" for k in _af}}

    _sizes = ""
    if have_reads:
        _pc = Counter(pileup_cls.values())
        _sizes += (f"- **pileup**: true100 {_pc.get('true100', 0):,} · "
                   f"false100 {_pc.get('false100', 0):,} · "
                   f"doubtful {_pc.get('doubtful', 0):,} · "
                   f"undetermined {_pc.get('undetermined', 0):,}\n")
    _hc = Counter(hip_cls.values())
    _sizes += (f"- **HipSTR**: true100 {_hc.get('true100', 0):,} · "
               f"false100 {_hc.get('false100', 0):,} · "
               f"doubtful {_hc.get('doubtful', 0):,}\n")

    _tbl = ("| standard | caller | true100 | FN | false100 | FP |\n"
            "|---|---|--:|--:|--:|--:|\n")
    if have_reads:
        _tbl += _score(pileup_cls, "pileup (our reads)", _callers)
    _tbl += _score(hip_cls, "HipSTR (its reads)", _callers)
    if have_reads:
        _tbl += _score(_core, "**confident core**", _callers)

    mo.md(
        f"### Silver-standard sizes\n\n{_sizes}\n"
        f"### False positives / false negatives\n\n{_tbl}\n"
        "> **Read the FP column, not the FN column, for a fair comparison.** "
        "Each standard's *true* set is partly defined by that caller's own "
        "calls, so FN is biased in that caller's favour (HipSTR shows 0% FN "
        "against its own and the core standards by construction; the pileup "
        "standard leans toward ours). **FP is robust** — calling a truly-"
        "monomorphic locus variable is a clear error regardless of read "
        "extraction, and it is comparable across callers."
    )
    return


@app.cell
def _(Counter, have_reads, hip_cls, hip_loci, hip_var, mo, ours_var,
      period_label, pileup_cls, plt, run_sel):
    mo.stop(run_sel is None or not have_reads, mo.md(""))
    _shared = set(pileup_cls) & set(hip_cls)
    _af = [k for k in _shared if pileup_cls[k] == "false100" == hip_cls[k]]
    _fp_h = [k for k in _af if hip_var(k) and not ours_var(k)]
    _per = Counter(hip_loci[k]["period"] for k in _fp_h)
    _periods = sorted(_per)
    _fig, _ax = plt.subplots(figsize=(8, 3.6))
    _ax.bar([period_label(p) for p in _periods], [_per[p] for p in _periods],
            color="#d62728")
    _ax.set_ylabel("HipSTR false-positive loci")
    _ax.set_title(f"HipSTR's {len(_fp_h)} false positives on the confident-false "
                  f"core (ours correctly silent), by motif length")
    _ax.grid(True, axis="y", alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(mo):
    mo.md("### Threshold sensitivity — is the precision gap stable?")
    return


@app.cell
def _(bal_min, catalog, classify_locus, have_reads, hip_cls, hip_signals,
      hip_var, min_assessable, mo, ours_var, pileup_signals, plt, reads, recur,
      run_sel):
    mo.stop(run_sel is None or not have_reads, mo.md(""))
    _depths = [3, 4, 5, 6]
    _ours_fp, _hip_fp = [], []
    for _md in _depths:
        _pc = {}
        for (_ch, _pos), _ps in reads.items():
            _meta = catalog.get((_ch, _pos))
            if _meta is None:
                continue
            _sig = pileup_signals(_ps, _meta[0], _meta[1], _md, 0.75, bal_min)
            _pc[(_ch, _pos)] = classify_locus(*_sig, recur, min_assessable)
        _af = [k for k in (set(_pc) & set(hip_cls))
               if _pc[k] == "false100" == hip_cls[k]]
        _n = max(len(_af), 1)
        _ours_fp.append(100 * sum(1 for k in _af if ours_var(k)) / _n)
        _hip_fp.append(100 * sum(1 for k in _af if hip_var(k)) / _n)
    _fig, _ax = plt.subplots(figsize=(8, 4))
    _ax.plot(_depths, _hip_fp, "o-", color="#1f77b4", label="HipSTR FP %")
    _ax.plot(_depths, _ours_fp, "o-", color="#2ca02c", label="ours FP %")
    _ax.set_xlabel("min reads/sample threshold")
    _ax.set_ylabel("false-positive rate on confident-false core (%)")
    _ax.set_title("Precision gap vs depth threshold (HOM_FRAC fixed 0.75)")
    _ax.set_xticks(_depths)
    _ax.legend()
    _ax.grid(True, alpha=0.3)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(mo):
    mo.md(
        """
        ### Caveats

        - Per-sample depth is shallow here (median 3 reads), so the standard
          leans on **recurrence across plants**, not deep single-sample calls.
        - It judges only the confident tails; the doubtful/undetermined middle
          is where the **orthogonal validation set (idea 5)** will matter.
        - Thresholds are exposed above — the sensitivity plot shows whether the
          precision gap survives moving them.
        """
    )
    return


@app.cell
def _():
    import matplotlib.pyplot as plt

    return (plt,)


if __name__ == "__main__":
    app.run()
