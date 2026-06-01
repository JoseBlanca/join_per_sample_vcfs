# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "matplotlib",
# ]
# ///

# Marimo performance dashboard — pop_var_caller vs freebayes vs GATK on
# the tomato1 benchmark cohort. Four sections, each reading the per-caller
# TSVs written by the perf_<caller>.py scripts; any missing TSV is skipped
# with a friendly note so partial runs still render.
#
#   1. Single-sample direct calling (CRAM -> VCF) — bar charts (time + RAM)
#      for freebayes, GATK (direct HaplotypeCaller), and pop_var_caller
#      (var-calling-from-bam, the direct single-sample route). All three
#      with a 4-thread budget.
#   2. Build one per-sample intermediate (one sample, 4 threads) — bar
#      charts for GATK (HaplotypeCaller -ERC GVCF) and pop_var_caller
#      (pileup -> .psp).
#   3. Scaling CRAM -> VCF with N samples — line charts. freebayes (region
#      parallel, 4 workers), GATK (one direct multi-sample HaplotypeCaller,
#      4 pair-HMM threads), pop_var_caller (4 parallel pileups then one
#      4-thread joint var-calling).
#   4. Scaling intermediate -> VCF with N samples (the re-run-on-new-samples
#      step) — line charts. GATK (CombineGVCFs + GenotypeGVCFs) and
#      pop_var_caller (.psp -> VCF, 4 threads).
#
# Some measurements are shared across sections: freebayes.tsv and
# gatk_direct.tsv each feed both section 1 (N=1 row) and section 3.
#
# Recommended invocation:
#
#   uvx marimo edit --sandbox benchmarks/tomato1/scripts/perf_dashboard.py
#
# Companion to (not replacement of) dashboard.py — that one is about
# correctness (variant agreement, QUAL distributions); this one is
# about performance.

import marimo

__generated_with = "0.23.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import csv
    from pathlib import Path

    import marimo as mo

    return Path, csv, mo


@app.cell
def _(mo):
    mo.md("""
    # pop_var_caller vs freebayes vs GATK — performance on tomato1

    Wall-clock time and peak RSS (summed across the process tree) on the
    tomato1 cohort (CRAMs pre-sliced to `regions.bed`). Four sections:

    - **1. Single-sample direct calling** (CRAM → VCF, bar charts).
      freebayes vs pop_var_caller (`var-calling-from-bam`, the direct
      single-sample route). 4-thread budget each (ours `--threads 4`;
      freebayes 4 concurrent region workers). _GATK is excluded here — its
      only "direct cram→vcf" path is multi-sample HaplotypeCaller, which
      doesn't scale (see §3 note); GATK appears in §2 and §4._
    - **2. Build one per-sample intermediate** (one sample, 4 threads, bar
      charts). GATK `HaplotypeCaller -ERC GVCF` → GVCF vs pop_var_caller
      `pileup` → `.psp`.
    - **3. Scaling CRAM → VCF, N samples** (line charts). freebayes
      (region-parallel) vs pop_var_caller (4 parallel pileups, then one
      4-thread joint `var-calling`). _GATK's direct multi-sample
      `HaplotypeCaller` is excluded: it re-assembles over pooled deep
      coverage, so wall is super-linear in N (≈3.5× per sample-doubling on
      tomato1) — impractical past a few samples. GATK's scalable route is
      the GVCF flow shown in §4._
    - **4. Scaling intermediate → VCF, N samples** (line charts; the step
      you re-run when new samples arrive). GATK `CombineGVCFs` +
      `GenotypeGVCFs` vs pop_var_caller `.psp` → VCF (4 threads).
    - **5. Thread scaling** (line charts) of pop_var_caller's joint stage
      (`.psp` → VCF) at the *maximum* cohort size, sweeping `--threads` —
      how well one `var-calling` process uses more cores.

    Inputs: per-caller TSVs at `results/perf/<caller>.tsv`, written by
    `perf_<caller>.py`. `freebayes.tsv` feeds both §1 (N=1 row) and §3.
    Section 5 reads `ours_joint_threads.tsv` (`perf_ours_joint_threads.py`).
    """)
    return


@app.cell
def _(Path):
    test_dir = Path(__file__).resolve().parent.parent
    perf_dir = test_dir / "results" / "perf"
    # Every TSV any section needs. freebayes is shared across §1 and §3.
    callers = (
        "ours_from_bam", "freebayes",      # §1 (freebayes also §3)
        "ours_psp_4t", "gatk_gvcf_4t",     # §2
        "ours_pileup_build",               # §3 stage-1 (+ ours_joint, freebayes)
        "ours_joint", "gatk_joint",        # §3 stage-2 + §4
    )
    tsv_paths = {c: perf_dir / f"{c}.tsv" for c in callers}
    return callers, perf_dir, test_dir, tsv_paths


@app.cell
def _(csv, mo, tsv_paths):
    # Load each TSV; skip missing ones with a friendly summary so a partial
    # run (e.g. you've only run perf_freebayes.py so far) still renders.
    rows_by_caller: dict[str, list[dict]] = {}
    present, missing = [], []
    for _caller, _path in tsv_paths.items():
        if not _path.exists():
            missing.append(_caller)
            continue
        with _path.open() as _fh:
            _reader = csv.DictReader(_fh, delimiter="\t")
            _rows = []
            for _r in _reader:
                _rows.append({
                    "n_samples": int(_r["n_samples"]),
                    "wall_seconds": float(_r["wall_seconds"]),
                    "peak_rss_mb": float(_r["peak_rss_mb"]),
                    "exit_code": int(_r["exit_code"]),
                })
        rows_by_caller[_caller] = _rows
        present.append(_caller)

    lines = [f"- **{c}**: {len(rows_by_caller[c])} rows" for c in present]
    if missing:
        lines.append("")
        lines.append("_Missing (not run yet):_ " + ", ".join(f"`{m}`" for m in missing))
    status_view = mo.md("### TSV inputs\n\n" + "\n".join(lines))
    status_view
    return missing, present, rows_by_caller, status_view


@app.cell
def _(csv, perf_dir):
    # Thread-scaling TSV for §5 — a DIFFERENT schema (keyed by `threads`,
    # not `n_samples`), written by perf_ours_joint_threads.py. Loaded
    # separately so the n_samples parsing above stays simple. Empty list
    # if the sweep hasn't been run yet.
    thread_tsv = perf_dir / "ours_joint_threads.tsv"
    thread_rows: list[dict] = []
    if thread_tsv.exists():
        with thread_tsv.open() as _fh:
            for _r in csv.DictReader(_fh, delimiter="\t"):
                thread_rows.append({
                    "threads": int(_r["threads"]),
                    "wall_seconds": float(_r["wall_seconds"]),
                    "peak_rss_mb": float(_r["peak_rss_mb"]),
                    "exit_code": int(_r["exit_code"]),
                })
    return (thread_rows,)


@app.cell
def _(mo):
    # Axis-scale toggles for the scaling (line) sections. Log-y is the
    # bigger lever when callers differ by an order of magnitude; log-x
    # helps read linear-vs-superlinear scaling across the 1..26 span.
    xscale = mo.ui.radio(
        options=["linear", "log"], value="linear", label="x-axis scale (N samples)",
    )
    yscale = mo.ui.radio(
        options=["linear", "log"], value="linear", label="y-axis scale (runtime / RSS)",
    )
    mo.hstack([xscale, yscale])
    return xscale, yscale


@app.cell
def _():
    import matplotlib.pyplot as plt

    # Stable colour per caller; ours = blue/green family, freebayes = red,
    # GATK = purple/brown/pink family.
    PALETTE = {
        "ours_from_bam":        "#1f77b4",  # blue
        "ours_psp_4t":          "#17becf",  # cyan
        "ours_pileup_build":    "#2ca02c",  # green
        "ours_cram_to_vcf":     "#2ca02c",  # green (§3 combined ours line)
        "ours_joint":           "#1a5d1a",  # dark green
        "freebayes":            "#d62728",  # red
        "gatk_gvcf_4t":         "#e377c2",  # pink
        "gatk_joint":           "#8c564b",  # brown
    }
    # Human-friendly labels for legends / bar x-ticks.
    LABELS = {
        "ours_from_bam":        "pop_var_caller\n(direct)",
        "ours_psp_4t":          "pop_var_caller\n(pileup→psp)",
        "ours_pileup_build":    "pop_var_caller\n(PSP build)",
        "ours_cram_to_vcf":     "pop_var_caller\n(pileup build + joint)",
        "ours_joint":           "pop_var_caller\n(psp→vcf)",
        "freebayes":            "freebayes",
        "gatk_gvcf_4t":         "GATK\n(HC→gvcf)",
        "gatk_joint":           "GATK\n(Combine+Genotype)",
    }

    def _row_at(rows, n):
        """The exit-0 measurement at N=n for one caller, or None."""
        for r in rows:
            if r["n_samples"] == n and r["exit_code"] == 0:
                return r
        return None

    def bar_pair(rows_by_caller, callers, *, target_n, title):
        """Two bar panels (runtime + peak RSS) over `callers` at a fixed
        sample count `target_n`. One bar per caller. Returns (fig,
        missing_callers); fig is None only if NO caller has a usable row."""
        picked = [(c, _row_at(rows_by_caller.get(c, []), target_n)) for c in callers]
        usable = [(c, r) for c, r in picked if r is not None]
        absent = [c for c, r in picked if r is None]
        if not usable:
            return None, list(callers)

        labels = [LABELS.get(c, c) for c, _ in usable]
        colours = [PALETTE.get(c) for c, _ in usable]
        ts = [r["wall_seconds"] for _, r in usable]
        ms = [r["peak_rss_mb"] for _, r in usable]
        xs = range(len(usable))

        fig, (ax_t, ax_m) = plt.subplots(1, 2, figsize=(13, 5))
        for ax, vals, ylabel, panel_title, fmt in (
            (ax_t, ts, "wall-clock seconds", "Runtime", "{:.1f}"),
            (ax_m, ms, "peak RSS (MB)",      "Peak memory", "{:.0f}"),
        ):
            bars = ax.bar(xs, vals, color=colours)
            ax.set_xticks(list(xs))
            ax.set_xticklabels(labels, fontsize=9)
            ax.set_ylabel(ylabel)
            ax.set_title(panel_title)
            ax.grid(True, axis="y", alpha=0.3)
            ax.bar_label(bars, labels=[fmt.format(v) for v in vals], padding=3, fontsize=9)
            ax.margins(y=0.15)
        fig.suptitle(title, fontsize=13, y=1.02)
        fig.tight_layout()
        return fig, absent

    def line_pair(rows_by_caller, callers, *, title, xscale_val, yscale_val):
        """Two line panels (runtime + peak RSS) vs N samples over
        `callers`. Returns (fig, missing_callers); fig is None only if NO
        caller is present.

        Runs that exited non-zero (e.g. an OOM-kill at high N) are NOT
        dropped — they're drawn as a distinct ✗ marker at the (N, value)
        the run reached before dying. For a memory-scaling comparison the
        point where a tool falls over is a headline result, not noise."""
        present = [c for c in callers if c in rows_by_caller]
        if not present:
            return None, list(callers)
        fig, (ax_t, ax_m) = plt.subplots(1, 2, figsize=(14, 5), sharex=True)
        failed_any = False
        for caller in present:
            rows = sorted(rows_by_caller[caller], key=lambda r: r["n_samples"])
            ok_rows = [r for r in rows if r["exit_code"] == 0]
            bad_rows = [r for r in rows if r["exit_code"] != 0]
            colour = PALETTE.get(caller)
            label = LABELS.get(caller, caller).replace("\n", " ")
            if ok_rows:
                ax_t.plot([r["n_samples"] for r in ok_rows],
                          [r["wall_seconds"] for r in ok_rows],
                          marker="o", label=label, color=colour)
                ax_m.plot([r["n_samples"] for r in ok_rows],
                          [r["peak_rss_mb"] for r in ok_rows],
                          marker="o", label=label, color=colour)
            for r in bad_rows:
                failed_any = True
                ax_t.scatter([r["n_samples"]], [r["wall_seconds"]],
                             marker="x", s=90, color=colour, zorder=5)
                ax_m.scatter([r["n_samples"]], [r["peak_rss_mb"]],
                             marker="x", s=90, color=colour, zorder=5)
            # Annotate the first failing N on the memory panel (where OOM bites).
            if bad_rows:
                first_bad = bad_rows[0]
                ax_m.annotate(f"{label.split()[0]} ✗ N={first_bad['n_samples']}",
                              (first_bad["n_samples"], first_bad["peak_rss_mb"]),
                              textcoords="offset points", xytext=(6, 6),
                              fontsize=8, color=colour)
        for ax, ylabel, panel_title in (
            (ax_t, "wall-clock seconds", "Runtime"),
            (ax_m, "peak RSS (MB)",      "Peak memory"),
        ):
            ax.set_xlabel("N samples")
            ax.set_ylabel(ylabel)
            ax.set_title(panel_title)
            ax.set_xscale(xscale_val)
            ax.set_yscale(yscale_val)
            ax.grid(True, alpha=0.3)
        # Legends last; append the failed-marker entry so it isn't clobbered.
        handles, _lab = ax_t.get_legend_handles_labels()
        if failed_any:
            from matplotlib.lines import Line2D
            handles = handles + [Line2D([], [], marker="x", linestyle="none",
                                        color="grey", label="run failed (exit≠0, e.g. OOM)")]
        ax_t.legend(handles=handles, fontsize=9)
        ax_m.legend(fontsize=9)
        fig.suptitle(title, fontsize=13, y=1.02)
        fig.tight_layout()
        return fig, []

    def thread_pair(thread_rows, *, title, yscale_val):
        """Two panels (runtime + peak RSS) vs thread count for a single
        caller's thread sweep. x-axis is always linear (threads is a
        small integer range); yscale follows the toggle. Returns the fig,
        or None if there are no usable rows."""
        ok_rows = sorted(
            (r for r in thread_rows if r["exit_code"] == 0),
            key=lambda r: r["threads"],
        )
        if not ok_rows:
            return None
        xs = [r["threads"] for r in ok_rows]
        ts = [r["wall_seconds"] for r in ok_rows]
        ms = [r["peak_rss_mb"] for r in ok_rows]
        colour = PALETTE.get("ours_joint")
        fig, (ax_t, ax_m) = plt.subplots(1, 2, figsize=(14, 5), sharex=True)
        for ax, vals, ylabel, panel_title in (
            (ax_t, ts, "wall-clock seconds", "Runtime"),
            (ax_m, ms, "peak RSS (MB)",      "Peak memory"),
        ):
            ax.plot(xs, vals, marker="o", color=colour)
            ax.set_xlabel("threads (--threads)")
            ax.set_ylabel(ylabel)
            ax.set_title(panel_title)
            ax.set_xticks(xs)
            ax.set_yscale(yscale_val)
            ax.grid(True, alpha=0.3)
        # Ideal-linear-speedup reference on the runtime panel (T=1 wall / T).
        t1 = next((r["wall_seconds"] for r in ok_rows if r["threads"] == 1), None)
        if t1 is not None:
            ax_t.plot(xs, [t1 / x for x in xs], linestyle="--", color="grey",
                      alpha=0.7, label="ideal linear speedup")
            ax_t.legend(fontsize=9)
        fig.suptitle(title, fontsize=13, y=1.02)
        fig.tight_layout()
        return fig

    return LABELS, PALETTE, bar_pair, line_pair, plt, thread_pair


@app.cell
def _(bar_pair, mo, rows_by_caller):
    # Section 1 — single-sample direct calling (CRAM -> VCF), N=1 bars.
    fig_1, missing_1 = bar_pair(
        rows_by_caller,
        ("freebayes", "ours_from_bam"),
        target_n=1,
        title="1. Single-sample direct calling (CRAM → VCF, 4-thread budget)",
    )
    sec1_methods = mo.md("""
    **Materials & methods.** One sample (N=1), CRAM → VCF, restricted to
    `regions.bed` (~2 Mb; CRAMs pre-sliced to it). Wall + peak RSS summed
    across the whole process tree.

    - **freebayes** — one process per BED region (20 regions), **up to 4
      concurrent processes, 1 thread each** (freebayes has no internal
      threading); per-region VCFs concatenated at the end.
    - **pop_var_caller (direct)** — **1 `var-calling-from-bam` process, 4
      threads** (`--threads 4`). Runs genome-wide (no BED flag), but the
      pre-sliced CRAM confines the work to the same regions.

    _GATK omitted: its only direct cram→vcf route is multi-sample
    HaplotypeCaller, excluded for the scaling reason noted in §3._
    """)
    _body_1 = mo.as_html(fig_1) if fig_1 is not None else mo.md(
        "_(Section 1 unavailable; missing TSV(s)/N=1 row: "
        + ", ".join(f"`{c}`" for c in missing_1) + ".)_"
    )
    sec1_view = mo.vstack([sec1_methods, _body_1])
    sec1_view
    return fig_1, missing_1, sec1_view


@app.cell
def _(bar_pair, mo, rows_by_caller):
    # Section 2 — build one per-sample intermediate (one sample, 4 threads).
    fig_2, missing_2 = bar_pair(
        rows_by_caller,
        ("gatk_gvcf_4t", "ours_psp_4t"),
        target_n=1,
        title="2. Build one per-sample intermediate (one sample, 4 threads)",
    )
    sec2_methods = mo.md("""
    **Materials & methods.** Build one sample's per-sample intermediate
    file, `--intervals regions.bed`. A single process per tool gets the
    full 4-thread budget; wall + peak RSS of that one process.

    - **GATK (HC → gvcf)** — **1 `HaplotypeCaller -ERC GVCF` process, 4
      pair-HMM threads** (`--native-pair-hmm-threads 4`).
    - **pop_var_caller (pileup → psp)** — **1 `pileup` process, 4 threads**
      (`--threads 4`).
    """)
    _body_2 = mo.as_html(fig_2) if fig_2 is not None else mo.md(
        "_(Section 2 unavailable; missing TSV(s): "
        + ", ".join(f"`{c}`" for c in missing_2) + ".)_"
    )
    sec2_view = mo.vstack([sec2_methods, _body_2])
    sec2_view
    return fig_2, missing_2, sec2_view


@app.cell
def _(line_pair, mo, rows_by_caller, xscale, yscale):
    # Section 3 — scaling CRAM -> VCF with N samples.
    # ours CRAM→VCF = build the PSPs ONCE (makespan per N) + joint psp→vcf
    # per N. A sample's .psp is cohort-size-independent, so the build is not
    # repeated per N; with a FIFO worker pool, makespan(first N) is exact.
    _build = {r["n_samples"]: r for r in rows_by_caller.get("ours_pileup_build", [])}
    _joint = {r["n_samples"]: r for r in rows_by_caller.get("ours_joint", [])}
    _combined = []
    for _n in sorted(set(_build) & set(_joint)):
        _b, _j = _build[_n], _joint[_n]
        _combined.append({
            "n_samples": _n,
            "wall_seconds": _b["wall_seconds"] + _j["wall_seconds"],
            "peak_rss_mb": max(_b["peak_rss_mb"], _j["peak_rss_mb"]),
            "exit_code": 0 if (_b["exit_code"] == 0 and _j["exit_code"] == 0) else 1,
        })
    _rows3 = dict(rows_by_caller)
    _rows3["ours_cram_to_vcf"] = _combined
    fig_3, missing_3 = line_pair(
        _rows3,
        ("freebayes", "ours_cram_to_vcf"),
        title="3. Scaling CRAM → VCF vs N samples",
        xscale_val=xscale.value,
        yscale_val=yscale.value,
    )
    sec3_methods = mo.md("""
    **Materials & methods.** Cohort CRAM → VCF as N grows, `regions.bed`.
    Parallelism budget = 4 for every tool. Wall + peak RSS over the whole
    process tree. A **✗** marks a run that died at that N (e.g. OOM-killed) —
    for a memory-scaling comparison, the point where a tool falls over is
    itself a headline result (pop_var_caller's design trades RAM for the
    ability to keep scaling to more samples).

    - **freebayes** — per-region processes, **up to 4 concurrent, 1 thread
      each**; each process calls all N samples in its region → per-region
      multi-sample VCFs, concatenated.
    - **pop_var_caller (pileup build + joint)** — plotted wall is
      **PSP-build makespan(N) + joint(N)**. Stage 1 builds the per-sample
      `.psp` files **once** through a FIFO pool of **4 `pileup` workers,
      1 thread each**; `makespan(N)` (wall until the first N are ready) is
      read off that single build — a `.psp` is cohort-size-independent, so
      nothing is re-called per N. Stage 2 is **1 `var-calling` process, 4
      threads** over the first N PSPs. Peak RSS = max(build peak ≈ the fixed
      4-pileup footprint, joint peak(N)).

    _GATK's direct multi-sample `HaplotypeCaller` (one process, all N
    `--input`s → joint VCF) is excluded: it re-assembles haplotypes over
    the pooled deep coverage of every sample per active region, so wall is
    super-linear in N (on tomato1: N=1→0.8 min, 2→2.1, 4→7.7, ≈3.5× per
    doubling) — impractical past a few samples. GATK's scalable path is the
    GVCF flow in §4._
    """)
    _body_3 = mo.as_html(fig_3) if fig_3 is not None else mo.md(
        "_(Section 3 unavailable; missing TSV(s): "
        + ", ".join(f"`{c}`" for c in missing_3) + ".)_"
    )
    sec3_view = mo.vstack([sec3_methods, _body_3])
    sec3_view
    return fig_3, missing_3, sec3_view


@app.cell
def _(line_pair, mo, rows_by_caller, xscale, yscale):
    # Section 4 — scaling intermediate -> VCF (re-run-on-new-samples step).
    fig_4, missing_4 = line_pair(
        rows_by_caller,
        ("gatk_joint", "ours_joint"),
        title="4. Scaling intermediate → VCF vs N samples (GVCF/PSP → cohort VCF)",
        xscale_val=xscale.value,
        yscale_val=yscale.value,
    )
    sec4_methods = mo.md("""
    **Materials & methods.** The re-genotyping step: pre-built per-sample
    intermediates → cohort VCF as N grows, `regions.bed`. Intermediates
    assumed already built (their build time is excluded here).

    - **GATK (Combine + Genotype)** — **2 sequential processes**:
      `CombineGVCFs` then `GenotypeGVCFs`, **1 thread each** (default; these
      tools are not pair-HMM-threaded). Wall = sum of the two; peak RSS =
      max of the two.
    - **pop_var_caller (psp → vcf)** — **1 `var-calling` process, 4
      threads** (`--threads 4`).
    """)
    _body_4 = mo.as_html(fig_4) if fig_4 is not None else mo.md(
        "_(Section 4 unavailable; missing TSV(s): "
        + ", ".join(f"`{c}`" for c in missing_4) + ".)_"
    )
    sec4_view = mo.vstack([sec4_methods, _body_4])
    sec4_view
    return fig_4, missing_4, sec4_view


@app.cell
def _(mo, thread_pair, thread_rows, yscale):
    # Section 5 — thread scaling of the cohort joint stage (PSP -> VCF) at
    # the MAXIMUM cohort size. One var-calling process, --threads swept.
    fig_5 = thread_pair(
        thread_rows,
        title="5. Thread scaling — pop_var_caller joint (PSP → VCF) at max N",
        yscale_val=yscale.value,
    )
    sec5_methods = mo.md("""
    **Materials & methods.** Thread scaling of the joint stage at the
    **maximum cohort** (all 26 PSPs), `.psp` → VCF. **1 `var-calling`
    process**, with `--threads` swept over {1, 2, 3, 4, 6, 8} — one process
    per point, no inter-process contention. Dashed line = ideal linear
    speedup extrapolated from the T=1 wall.
    """)
    _body_5 = mo.as_html(fig_5) if fig_5 is not None else mo.md(
        "_(Section 5 unavailable; run `perf_ours_joint_threads.py` to write "
        "`ours_joint_threads.tsv`.)_"
    )
    sec5_view = mo.vstack([sec5_methods, _body_5])
    sec5_view
    return fig_5, sec5_view


@app.cell
def _(LABELS, mo, rows_by_caller):
    # Raw numbers — exact values per caller. Failures (exit_code != 0) get
    # a ⚠ marker.
    if not rows_by_caller:
        table_view = mo.md("")
    else:
        blocks = []
        for _caller, _rows in rows_by_caller.items():
            _label = LABELS.get(_caller, _caller).replace("\n", " ")
            _md = f"#### {_caller} — {_label}\n\n"
            _md += "| N | wall (s) | peak RSS (MB) | exit |\n|---:|---:|---:|---:|\n"
            for _r in sorted(_rows, key=lambda x: x["n_samples"]):
                _marker = "" if _r["exit_code"] == 0 else " ⚠"
                _md += (
                    f"| {_r['n_samples']} | {_r['wall_seconds']:.1f} | "
                    f"{_r['peak_rss_mb']:.0f} | {_r['exit_code']}{_marker} |\n"
                )
            blocks.append(_md)
        table_view = mo.md("### Raw measurements\n\n" + "\n".join(blocks))
    table_view
    return (table_view,)


if __name__ == "__main__":
    app.run()
