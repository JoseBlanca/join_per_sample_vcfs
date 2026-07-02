#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy", "matplotlib"]
# ///
"""Multi-panel, x-aligned figure with observed heterozygosity on the x-axis:
  A         : distribution of obs_het (how many loci at each het level)
  B, C, ... : boxplots of each candidate paralog measure, binned by obs_het

Reads the per-locus heterozygosity from the het side: for loci grouped by their
observed het, how do coverage / n_excess / divergence / co-segregation behave?
Saves a PNG next to the feature table."""
import polars as pl
import numpy as np
import matplotlib.pyplot as plt

FEATURES = "benchmarks/tomato2/results/snp_features.parquet"
OUT = "benchmarks/tomato2/results/het_measures_figure.png"
STEP = 0.05   # obs_het bin width

# the measures we're evaluating (col, label, y-limits for display)
MEASURES = [
    ("mean_rel_cov", "mean rel coverage\n(1 = single copy)", (0, 3)),
    ("n_excess", "n_excess samples", (0, 12)),
    ("mqdiff", "MQAlt - MQRef\n(divergence)", (-25, 10)),
    ("cov_het_gap", "cov(het) - cov(hom)\n(co-segregation)", (-2, 3)),
]

df = pl.read_parquet(FEATURES).filter(pl.col("n_called") >= 30)
het = df["obs_het"].to_numpy()
edges = np.arange(0.0, 1.0 + STEP / 2, STEP)
centers = (edges[:-1] + edges[1:]) / 2
binidx = np.digitize(np.clip(het, 0, 1 - 1e-9), edges) - 1

nrows = 1 + len(MEASURES)
fig, axes = plt.subplots(
    nrows, 1, figsize=(11, 2.4 * nrows), sharex=True,
    gridspec_kw={"height_ratios": [1] + [1.2] * len(MEASURES), "hspace": 0.12},
)

# --- Panel A: het distribution ---
axA = axes[0]
axA.hist(het, bins=edges, color="seagreen", edgecolor="white", linewidth=0.4)
axA.set_yscale("log")
axA.set_ylabel("loci per bin\n(log)")
axA.set_title(f"tomato2 cohort — {df.height:,} biallelic SNPs (n_called≥30)")

# --- Panels B..: each measure boxplotted per het bin ---
for ax, (col, label, ylim) in zip(axes[1:], MEASURES):
    vals = df[col].to_numpy().astype(float)
    by_bin = [vals[(binidx == i) & np.isfinite(vals)] for i in range(len(centers))]
    bp = ax.boxplot(
        by_bin, positions=centers, widths=STEP * 0.8, showfliers=False,
        patch_artist=True, manage_ticks=False, medianprops=dict(color="black"),
    )
    for patch in bp["boxes"]:
        patch.set(facecolor="lightsteelblue", alpha=0.9)
    med = [np.median(b) if len(b) else np.nan for b in by_bin]
    ax.plot(centers, med, color="crimson", lw=1.3, marker="o", ms=2.5)
    ax.set_ylim(*ylim)
    ax.set_ylabel(label, fontsize=9)
    if col in ("mqdiff", "cov_het_gap"):
        ax.axhline(0, color="grey", ls=":", lw=0.8)

axes[-1].set_xlabel("observed heterozygosity")
axes[-1].set_xlim(0, 1)
fig.savefig(OUT, dpi=130, bbox_inches="tight")
print(f"wrote {OUT}")
# per-bin counts so we know which high-het bins are thin
counts = [int(((binidx == i)).sum()) for i in range(len(centers))]
print("loci per het bin:")
for c, n in zip(centers, counts):
    print(f"  het~{c:.2f}  n={n}")
