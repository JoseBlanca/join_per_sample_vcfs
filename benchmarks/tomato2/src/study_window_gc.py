#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy", "matplotlib"]
# ///
"""Study (before deciding the normalisation): the per-sample window depth-vs-GC
relationship. Each sample's window mean depth is normalised by its own median
(relative depth) so shapes are comparable across samples of different depth.

Panels:
  A  relative depth vs GC, cohort density + binned median (the GC effect)
  B  per-sample GC curves (is the shape consistent across samples?)
  C  relative-depth distribution at window scale (single-copy mode + paralog tail)
  D  GC distribution of the windows (where we have data)
"""
import polars as pl
import numpy as np
import matplotlib.pyplot as plt

W = 2000
df = pl.read_parquet(f"benchmarks/tomato2/results/window_cov.w{W}.parquet")

# per-sample median depth -> relative depth (single-copy normalisation, no GC yet)
med = df.group_by("sample").agg(pl.col("mean_depth").median().alias("smed"))
df = df.join(med, on="sample").with_columns((pl.col("mean_depth") / pl.col("smed")).alias("rel"))
# keep well-covered, non-gap windows for the GC curve
ok = df.filter((pl.col("breadth") > 0.8) & (pl.col("gc") > 0.12) & (pl.col("gc") < 0.55))
print(f"{df.height:,} rows; {ok.height:,} used for GC curve")

gc = ok["gc"].to_numpy()
rel = ok["rel"].to_numpy()

fig, axes = plt.subplots(2, 2, figsize=(14, 9))

# A: relative depth vs GC, density + binned median
axA = axes[0, 0]
axA.hexbin(gc, np.clip(rel, 0, 2.5), gridsize=60, bins="log", cmap="Blues", mincnt=1)
edges = np.linspace(0.15, 0.52, 38)
ctr = (edges[:-1] + edges[1:]) / 2
binidx = np.digitize(gc, edges) - 1
med_by = [np.median(rel[binidx == i]) if (binidx == i).any() else np.nan for i in range(len(ctr))]
q1 = [np.percentile(rel[binidx == i], 25) if (binidx == i).any() else np.nan for i in range(len(ctr))]
q3 = [np.percentile(rel[binidx == i], 75) if (binidx == i).any() else np.nan for i in range(len(ctr))]
axA.plot(ctr, med_by, color="crimson", lw=2, label="median")
axA.plot(ctr, q1, color="crimson", lw=0.8, ls="--")
axA.plot(ctr, q3, color="crimson", lw=0.8, ls="--", label="IQR")
axA.axhline(1.0, color="grey", ls=":", lw=1)
axA.set_xlabel("GC fraction of 2kb window")
axA.set_ylabel("relative depth (/ sample median)")
axA.set_title("A. depth vs GC (cohort)")
axA.legend(fontsize=8)
axA.set_ylim(0, 2.5)

# B: per-sample GC curves
axB = axes[0, 1]
samples = sorted(df["sample"].unique().to_list())[:8]
for s in samples:
    d = ok.filter(pl.col("sample") == s)
    g = d["gc"].to_numpy(); r = d["rel"].to_numpy()
    bi = np.digitize(g, edges) - 1
    m = [np.median(r[bi == i]) if (bi == i).any() else np.nan for i in range(len(ctr))]
    axB.plot(ctr, m, lw=1, alpha=0.8)
axB.axhline(1.0, color="grey", ls=":", lw=1)
axB.set_xlabel("GC fraction")
axB.set_ylabel("relative depth")
axB.set_title(f"B. per-sample GC curves (first {len(samples)} samples)")
axB.set_ylim(0.4, 1.6)

# C: relative-depth distribution at window scale
axC = axes[1, 0]
allrel = df.filter(pl.col("breadth") > 0.8)["rel"].to_numpy()
axC.hist(np.clip(allrel, 0, 3), bins=120, color="steelblue")
axC.axvline(1.0, color="grey", ls=":", lw=1)
axC.set_yscale("log")
axC.set_xlabel("relative depth (window)")
axC.set_ylabel("windows (log)")
axC.set_title("C. window relative-depth distribution")

# D: GC distribution
axD = axes[1, 1]
axD.hist(df.filter(pl.col("breadth") > 0.8)["gc"].to_numpy(), bins=80, color="seagreen")
axD.set_xlabel("GC fraction of window")
axD.set_ylabel("window x sample rows")
axD.set_title("D. window GC distribution")

fig.tight_layout()
out = f"benchmarks/tomato2/results/window_gc_study.w{W}.png"
fig.savefig(out, dpi=130, bbox_inches="tight")
print(f"wrote {out}")

# quantify the GC effect
print("\nmedian relative depth by GC bin (the curve to correct):")
for c, m, n in zip(ctr[::4], med_by[::4], range(len(ctr))[::4]):
    cnt = int((binidx == n).sum())
    if not np.isnan(m):
        print(f"  GC~{c:.2f}  median_rel={m:.2f}  n={cnt}")
lo, hi = np.nanmin(med_by), np.nanmax(med_by)
print(f"\nGC effect span: median rel depth ranges {lo:.2f}–{hi:.2f} ({100*(hi-lo):.0f}% swing) across GC")
print(f"windows with rel>1.5 (breadth>0.8): {100*np.mean(allrel>1.5):.1f}%   rel>2: {100*np.mean(allrel>2):.1f}%")
