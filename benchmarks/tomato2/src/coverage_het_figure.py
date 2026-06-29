#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy", "matplotlib"]
# ///
"""Two-panel, x-aligned figure:
  A (top)    : distribution of per-locus mean relative coverage (1.0 = single copy)
  B (bottom) : boxplot of observed heterozygosity for the loci in each coverage bin

Shows directly how obs het rises once coverage exceeds single-copy — the
collapsed-paralog signature. Saves a PNG next to the feature table."""
import polars as pl
import numpy as np
import matplotlib.pyplot as plt

FEATURES = "benchmarks/tomato2/results/snp_features.parquet"
OUT = "benchmarks/tomato2/results/coverage_het_figure.png"

XCOL = "mean_rel_cov"      # per-locus coverage measure (1.0 = single copy)
XLO, XHI, STEP = 0.3, 3.0, 0.15   # coverage axis; loci >XHI folded into top bin

df = pl.read_parquet(FEATURES).filter(pl.col("n_called") >= 30)
x = df[XCOL].to_numpy()
het = df["obs_het"].to_numpy()

edges = np.arange(XLO, XHI + STEP / 2, STEP)
centers = (edges[:-1] + edges[1:]) / 2
xc = np.clip(x, XLO + 1e-9, XHI - 1e-9)   # fold the long tail into the end bins
binidx = np.digitize(xc, edges) - 1        # 0..len(centers)-1
n_over = int((x > XHI).sum())

# per-bin obs_het for the boxplots
het_by_bin = [het[(binidx == i) & np.isfinite(het)] for i in range(len(centers))]
counts = [len(h) for h in het_by_bin]

fig, (axA, axB) = plt.subplots(
    2, 1, figsize=(11, 7), sharex=True,
    gridspec_kw={"height_ratios": [1, 1.4], "hspace": 0.07},
)

# --- Panel A: coverage distribution ---
axA.hist(xc, bins=edges, color="steelblue", edgecolor="white", linewidth=0.4)
axA.axvline(1.0, color="grey", ls=":", lw=1.2)
axA.set_yscale("log")
axA.set_ylabel("loci per bin (log)")
axA.set_title(
    f"tomato2 cohort — {df.height:,} biallelic SNPs (n_called≥30); "
    f"top bin folds in {n_over:,} loci with coverage >{XHI:g}×"
)

# --- Panel B: obs het per coverage bin (x-aligned boxplots) ---
bp = axB.boxplot(
    het_by_bin, positions=centers, widths=STEP * 0.8,
    showfliers=False, patch_artist=True, manage_ticks=False,
    medianprops=dict(color="black"),
)
for patch in bp["boxes"]:
    patch.set(facecolor="lightsteelblue", alpha=0.9)
# overlay per-bin median trend
medians = [np.median(h) if len(h) else np.nan for h in het_by_bin]
axB.plot(centers, medians, color="crimson", lw=1.5, marker="o", ms=3, label="median obs_het")
axB.axvline(1.0, color="grey", ls=":", lw=1.2)
axB.set_ylim(0, 1)
axB.set_ylabel("observed heterozygosity")
axB.set_xlabel(f"{XCOL}  (mean relative coverage; 1.0 = single copy)")
axB.legend(loc="upper left", fontsize=9)
# annotate bin counts under the axis
for c, n in zip(centers, counts):
    if n:
        axB.text(c, -0.06, str(n), ha="center", va="top", fontsize=6, color="grey",
                 rotation=90, clip_on=False)

axB.set_xlim(XLO, XHI)
fig.savefig(OUT, dpi=130, bbox_inches="tight")
print(f"wrote {OUT}")
print("per-bin median obs_het:")
for c, m, n in zip(centers, medians, counts):
    print(f"  cov~{c:.2f}  n={n:6d}  median_het={m:.3f}")
