# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "matplotlib",
#     "matplotlib-venn",
# ]
# ///

# Marimo dashboard comparing the three SNP callers (pop_var_caller,
# freebayes, GATK) on the tomato cohort test set. Toggle between the
# single-sample and cohort variants of the test.
#
# Recommended invocation (no global marimo install needed):
#
#   uvx marimo edit --sandbox tmp/tomato_cohort_test/scripts/dashboard.py
#
# To serve as a read-only app instead of opening the editor:
#
#   uvx marimo run --sandbox tmp/tomato_cohort_test/scripts/dashboard.py
#
# Measures (added one at a time):
#   1. Variant agreement — 3-way Venn + counts table at (CHROM, POS,
#      REF, ALT) granularity. No PASS / QUAL filter applied; multi-
#      allelic records are split per ALT; `*` (gVCF spanning deletion)
#      is skipped. NB: indel calls aren't normalised, so the same
#      indel may show as discordant if callers anchored it on
#      different REF bases. Add `bcftools norm` upstream if that
#      matters for a measure.
#   2. Pairwise QUAL hexbins on shared variants — three 2D density
#      plots (ours-vs-fb, ours-vs-gatk, fb-vs-gatk) with a log colour
#      scale and a y=x reference line; bounded by the QUAL-cap slider.
#   3. Per-caller QUAL distributions — three vertically stacked
#      histograms sharing the same x-axis (0..cap) and bin edges;
#      independent y-axes (callers can differ by an order of magnitude
#      in record count); y-scale linear or log via radio toggle.

import marimo

__generated_with = "0.23.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import gzip
    from pathlib import Path

    import marimo as mo

    return Path, gzip, mo


@app.cell
def _(mo):
    mo.md("""
    # SNP caller comparison — tomato cohort test
    """)
    return


@app.cell
def _(mo):
    mode = mo.ui.radio(
        options=["single", "cohort"],
        value="single",
        label="Test set",
    )
    mode
    return (mode,)


@app.cell
def _(Path, mode):
    # tmp/tomato_cohort_test/scripts/dashboard.py -> tmp/tomato_cohort_test/
    test_dir = Path(__file__).resolve().parent.parent
    results = test_dir / "results"
    sample = "SRR17274057"
    if mode.value == "single":
        vcf_paths = {
            "ours": results / "ours" / f"single_{sample}.vcf",
            "freebayes": results / "freebayes" / f"single_{sample}.vcf",
            "gatk": results / "gatk" / f"single_{sample}.vcf",
        }
    else:
        vcf_paths = {
            "ours": results / "ours" / "cohort" / "cohort.vcf",
            "freebayes": results / "freebayes" / "cohort.vcf",
            "gatk": results / "gatk" / "cohort" / "cohort.vcf",
        }
    return (vcf_paths,)


@app.cell
def _(mo, vcf_paths):
    missing = {name: p for name, p in vcf_paths.items() if not p.exists()}
    if missing:
        miss_md = "\n".join(f"- **{name}**: `{p}`" for name, p in missing.items())
        inputs_view = mo.callout(
            mo.md(f"### Missing VCFs\n\n{miss_md}\n\nRun the corresponding `run_*` script first."),
            kind="danger",
        )
    else:
        present_md = "\n".join(
            f"- **{name}**: `{p}` ({p.stat().st_size / 1e6:.2f} MB)"
            for name, p in vcf_paths.items()
        )
        inputs_view = mo.md(f"### Inputs\n\n{present_md}")
    inputs_view
    return inputs_view, missing


@app.cell
def _(gzip):
    def variant_keys(path):
        """Read a VCF and return a set of (CHROM, POS, REF, ALT) tuples.

        Multi-allelic records are split into one entry per ALT.
        `*` (gVCF spanning deletion marker) is skipped — it is not a
        real variant.  No PASS / QUAL filter is applied here; callers
        of this helper layer filtering on top if needed.
        """
        opener = gzip.open if str(path).endswith(".gz") else open
        keys: set[tuple[str, int, str, str]] = set()
        with opener(path, "rt") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alts = fields[4]
                for alt in alts.split(","):
                    if alt == "*" or alt == ".":
                        continue
                    keys.add((chrom, pos, ref, alt))
        return keys

    def variant_quals(path):
        """Return a dict mapping (CHROM, POS, REF, ALT) -> float QUAL.

        QUAL is a *site*-level statistic, so every ALT of a
        multi-allelic record inherits the same value.  Records with
        QUAL = `.` or an unparsable value are skipped — the hexbin
        plots cannot consume missing/non-numeric coordinates.
        """
        opener = gzip.open if str(path).endswith(".gz") else open
        result: dict[tuple[str, int, str, str], float] = {}
        with opener(path, "rt") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                qual_str = fields[5]
                if qual_str in (".", ""):
                    continue
                try:
                    qual = float(qual_str)
                except ValueError:
                    continue
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alts = fields[4]
                for alt in alts.split(","):
                    if alt == "*" or alt == ".":
                        continue
                    result[(chrom, pos, ref, alt)] = qual
        return result

    return variant_keys, variant_quals


@app.cell
def _(missing, variant_keys, vcf_paths):
    # Cell short-circuits when inputs are missing so downstream cells
    # don't crash on a half-built dashboard.
    if missing:
        sets: dict[str, set] = {}
    else:
        sets = {name: variant_keys(p) for name, p in vcf_paths.items()}
    return (sets,)


@app.cell
def _(sets: dict[str, set]):
    # `ours`, `fb`, `gk` are always defined (possibly empty) so the
    # cell's return contract holds regardless of whether the VCFs
    # were found upstream.
    ours = sets.get("ours", set())
    fb = sets.get("freebayes", set())
    gk = sets.get("gatk", set())
    agreement_counts = {
        "in_all_3": len(ours & fb & gk),
        "ours_and_freebayes_only": len(ours & fb - gk),
        "ours_and_gatk_only": len(ours & gk - fb),
        "freebayes_and_gatk_only": len(fb & gk - ours),
        "only_ours": len(ours - fb - gk),
        "only_freebayes": len(fb - ours - gk),
        "only_gatk": len(gk - ours - fb),
        "total_ours": len(ours),
        "total_freebayes": len(fb),
        "total_gatk": len(gk),
    }
    return agreement_counts, fb, gk, ours


@app.cell
def _(agreement_counts, mo):
    # 7 Venn regions + 3 totals, rendered as a markdown table so it
    # works in app/run mode too.
    c = agreement_counts
    denom = max(c["in_all_3"], 1)
    rows = [
        ("In all 3 callers", c["in_all_3"], "100% of triple-agreement baseline"),
        ("ours ∩ freebayes (not gatk)", c["ours_and_freebayes_only"], f"{c['ours_and_freebayes_only'] / denom:.1%}"),
        ("ours ∩ gatk (not freebayes)", c["ours_and_gatk_only"], f"{c['ours_and_gatk_only'] / denom:.1%}"),
        ("freebayes ∩ gatk (not ours)", c["freebayes_and_gatk_only"], f"{c['freebayes_and_gatk_only'] / denom:.1%}"),
        ("only ours", c["only_ours"], f"{c['only_ours'] / denom:.1%}"),
        ("only freebayes", c["only_freebayes"], f"{c['only_freebayes'] / denom:.1%}"),
        ("only gatk", c["only_gatk"], f"{c['only_gatk'] / denom:.1%}"),
        ("— totals —", "", ""),
        ("ours (any membership)", c["total_ours"], ""),
        ("freebayes (any membership)", c["total_freebayes"], ""),
        ("gatk (any membership)", c["total_gatk"], ""),
    ]
    table_md = "| set | count | ratio vs all-3 |\n|---|---:|---:|\n"
    for name, count, ratio in rows:
        table_md += f"| {name} | {count} | {ratio} |\n"
    counts_view = mo.md(
        f"### 1. Variant agreement\n\n"
        f"_(key: CHROM, POS, REF, ALT; multi-allelic split per ALT)_\n\n"
        f"{table_md}"
    )
    counts_view
    return (counts_view,)


@app.cell
def _(fb, gk, ours):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib_venn import venn3

    labels = ("pop_var_caller", "freebayes", "gatk")
    fig, ax = plt.subplots(figsize=(8, 6))
    v = venn3([ours, fb, gk], set_labels=labels, ax=ax)
    ax.set_title("Variant agreement — (CHROM, POS, REF, ALT)")
    # Build a colour-keyed legend out of the venn patches so each
    # caller's circle colour matches a row in the legend (the set
    # labels matplotlib_venn draws next to the circles can be hard
    # to associate at a glance).
    patch_ids = ("100", "010", "001")  # ours, freebayes, gatk (A,B,C order)
    handles = [
        Patch(
            facecolor=v.get_patch_by_id(pid).get_facecolor(),
            edgecolor="black",
            label=lbl,
        )
        for pid, lbl in zip(patch_ids, labels)
        if v.get_patch_by_id(pid) is not None
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.0))
    fig
    return Patch, ax, fig, handles, labels, patch_ids, plt, v, venn3


@app.cell
def _(missing, variant_quals, vcf_paths):
    # Per-caller QUAL maps. Same shape as `sets` upstream but mapping
    # the variant key to its QUAL float rather than just membership.
    if missing:
        quals: dict[str, dict] = {}
    else:
        quals = {name: variant_quals(p) for name, p in vcf_paths.items()}
    return (quals,)


@app.cell
def _(mo, quals):
    # Slider to cap the QUAL axis range on the hexbin plots, so the
    # dense low/mid-QUAL bulk is readable instead of getting squashed
    # by a handful of very-high-QUAL outliers. Default value is the
    # 95th percentile of pooled QUALs — focuses on the bulk while
    # still showing some tail. Re-drawn live as the slider moves.
    if quals:
        all_quals = sorted(q for caller in quals.values() for q in caller.values())
        data_max = float(all_quals[-1])
        p95 = float(all_quals[int(len(all_quals) * 0.95)])
    else:
        # Inert defaults so the widget still renders when VCFs are
        # missing upstream.
        all_quals = []
        data_max = 100.0
        p95 = 100.0
    # ~500 slider positions across the data range — gives smooth
    # sliding for both small (data_max ~ 100) and large
    # (data_max ~ 10 000) value ranges.
    step = max(1.0, data_max / 500)
    qual_cap = mo.ui.slider(
        start=1.0,
        stop=data_max,
        step=step,
        value=p95,
        label="Max QUAL on hexbin axes",
        show_value=True,
        full_width=True,
    )
    qual_cap
    return all_quals, data_max, p95, qual_cap, step


@app.cell
def _(mo, plt, qual_cap, quals):
    # Three pairwise hexbins of QUAL on the variants both callers
    # called. Hexbin (vs scatter) handles the 10⁴-record cohort
    # output cleanly; log colour scale so dense ridges and sparse
    # tails are both readable. A red y=x reference line shows where
    # the two callers would agree on confidence — for callers using
    # different QUAL conventions (we know freebayes / GATK / ours all
    # define -10·log10 P over slightly different events), the bulk
    # rarely sits on the diagonal; the *slope* and *spread* matter.
    #
    # `plt` is dependency-injected from the Venn cell, which already
    # owns the matplotlib import — marimo enforces single-definition
    # for every name across cells.
    if not quals:
        hex_view = mo.md("_(QUAL hexbins unavailable until all three VCFs exist.)_")
    else:
        cap = float(qual_cap.value)
        pairs = [("ours", "freebayes"), ("ours", "gatk"), ("freebayes", "gatk")]
        hex_fig, hex_axes = plt.subplots(1, 3, figsize=(18, 5.5))
        for hex_ax, (a, b) in zip(hex_axes, pairs):
            qa = quals[a]
            qb = quals[b]
            all_shared = set(qa) & set(qb)
            # Apply the cap symmetrically: both axes need to fall in
            # range, otherwise we'd be looking at variants where one
            # caller is extreme and the other isn't — interesting but
            # for a different view.
            in_cap = [(qa[k], qb[k]) for k in all_shared if qa[k] <= cap and qb[k] <= cap]
            if not in_cap:
                hex_ax.set_title(f"{a} vs {b} — no shared variants ≤ cap")
                hex_ax.set_axis_off()
                continue
            xs = [p[0] for p in in_cap]
            ys = [p[1] for p in in_cap]
            hb = hex_ax.hexbin(
                xs, ys,
                gridsize=40,
                mincnt=1,
                cmap="viridis",
                bins="log",
                extent=(0, cap, 0, cap),
            )
            hex_ax.plot([0, cap], [0, cap], "r--", alpha=0.5, linewidth=1, label="y = x")
            hex_ax.set_xlim(0, cap)
            hex_ax.set_ylim(0, cap)
            hex_ax.set_xlabel(f"{a} QUAL")
            hex_ax.set_ylabel(f"{b} QUAL")
            hex_ax.set_title(
                f"{a} vs {b}  ({len(in_cap)} of {len(all_shared)} ≤ cap)"
            )
            hex_fig.colorbar(hb, ax=hex_ax, label="count (log)")
            hex_ax.legend(loc="upper left", fontsize=8)
        hex_fig.suptitle(
            f"2. Pairwise QUAL agreement on shared variants  (cap = {cap:.0f})",
            y=1.02,
            fontsize=13,
        )
        hex_fig.tight_layout()
        hex_view = hex_fig
    hex_view
    return (hex_view,)


@app.cell
def _(mo):
    # Y-axis scale for the per-caller QUAL distributions below.
    # Default linear; flip to log when a heavy low-QUAL tail dominates
    # the bulk and squashes shape detail.
    qual_yscale = mo.ui.radio(
        options=["linear", "log"],
        value="linear",
        label="QUAL distribution y-scale",
    )
    qual_yscale
    return (qual_yscale,)


@app.cell
def _(mo, plt, qual_cap, qual_yscale, quals):
    # Per-caller QUAL distributions. Three panels stacked vertically,
    # all sharing the same x-axis (0..cap) and bin edges so bar
    # positions line up across rows; y-axes are independent because
    # callers can differ by an order of magnitude in record count
    # (freebayes typically emits many more low-QUAL sites than the
    # others) and a shared y would squash the smaller panels.
    if not quals:
        dist_view = mo.md("_(QUAL distributions unavailable until all three VCFs exist.)_")
    else:
        dist_cap = float(qual_cap.value)
        n_bins = 50
        bin_edges = [dist_cap * i / n_bins for i in range(n_bins + 1)]
        callers = ("ours", "freebayes", "gatk")
        dist_fig, dist_axes = plt.subplots(
            len(callers), 1,
            figsize=(10, 8),
            sharex=True,
            sharey=False,
        )
        for dist_ax, _name in zip(dist_axes, callers):
            _values = [q for q in quals[_name].values() if q <= dist_cap]
            _total = len(quals[_name])
            dist_ax.hist(_values, bins=bin_edges, color="steelblue", edgecolor="white")
            dist_ax.set_yscale(qual_yscale.value)
            dist_ax.set_ylabel(f"{_name}\ncount")
            dist_ax.annotate(
                f"n = {len(_values):,} of {_total:,} ≤ cap",
                xy=(0.99, 0.92),
                xycoords="axes fraction",
                ha="right",
                va="top",
                fontsize=9,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7),
            )
        dist_axes[-1].set_xlabel("QUAL")
        dist_axes[-1].set_xlim(0.0, dist_cap)
        dist_fig.suptitle(
            f"3. Per-caller QUAL distributions  "
            f"(cap = {dist_cap:.0f}, {n_bins} bins, y={qual_yscale.value})",
            y=0.995,
            fontsize=13,
        )
        dist_fig.tight_layout()
        dist_view = dist_fig
    dist_view
    return (dist_view,)


if __name__ == "__main__":
    app.run()
