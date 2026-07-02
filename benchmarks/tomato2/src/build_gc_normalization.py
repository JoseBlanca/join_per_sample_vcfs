#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy", "matplotlib"]
# ///
"""Normalisation step of option (c): GC-correct the per-sample window coverage.

We tested the "pool the GC shape across the whole dataset" idea and REFUTED it
here: after a pooled correction, samples retain large *systematic* GC-correlated
residuals (median 17%, p90 40%) — their GC shapes genuinely differ (library/PCR
heterogeneity across SRA accessions). And pooling isn't needed for robustness:
each sample has ~16k windows, so its own GC curve is estimable to ~1-2% even at
5x. So we fit a PER-SAMPLE curve (binned median, robust to the paralog tail),
falling back to the pooled curve only in the sparse GC extremes.

  expected_depth(sample s, window w) = scale_s * curve_s(GC_w)
  gc_rel = observed_depth / expected_depth        (1.0 = single copy)

We report both the pooled and per-sample residuals so the choice is evidenced.
Output: results/window_cov.w<W>.gcnorm.parquet  (+ a validation figure)
"""
import argparse
import numpy as np
import polars as pl
import matplotlib.pyplot as plt

_ap = argparse.ArgumentParser()
_ap.add_argument("--w", type=int, default=2000)
W = _ap.parse_args().w
df = pl.read_parquet(f"benchmarks/tomato2/results/window_cov.w{W}.parquet")
well = df.filter(pl.col("breadth") > 0.8)

# per-sample scale = median window depth (GC-neutral: all samples share the same
# windows, hence the same GC composition)
scale = well.group_by("sample").agg(pl.col("mean_depth").median().alias("scale"))
df = df.join(scale, on="sample")
df = df.with_columns((pl.col("mean_depth") / pl.col("scale")).alias("rel"))
well = df.filter(pl.col("breadth") > 0.8)

# pooled GC curve: robust median relative depth per fine GC bin over ALL samples
GLO, GHI, STEP = 0.12, 0.55, 0.005
edges = np.arange(GLO, GHI + STEP / 2, STEP)
centers = (edges[:-1] + edges[1:]) / 2
gc = well["gc"].to_numpy()
rel = well["rel"].to_numpy()
bi = np.digitize(gc, edges) - 1
pooled = np.array([np.median(rel[bi == i]) if (bi == i).sum() >= 50 else np.nan
                   for i in range(len(centers))])
# fill thin-data bins by carrying nearest valid value (clamp at the ends)
valid = ~np.isnan(pooled)
pooled = np.interp(centers, centers[valid], pooled[valid])


def curve(g):
    return np.interp(g, centers, pooled)


# ---- PER-SAMPLE GC curves (primary), pooled only as sparse-GC fallback ----
samples = sorted(df["sample"].unique().to_list())
sample_curves = {}
raw_curves = {}  # un-smoothed, for plotting the "shapes differ" fan
for s in samples:
    d = well.filter(pl.col("sample") == s)
    g = d["gc"].to_numpy(); r = d["rel"].to_numpy()
    b = np.digitize(g, edges) - 1
    cur = np.array([np.median(r[b == i]) if (b == i).sum() >= 50 else np.nan
                    for i in range(len(centers))])
    raw_curves[s] = cur.copy()
    cur = np.where(np.isnan(cur), pooled, cur)          # fallback to pooled where sparse
    sm = np.array([np.median(cur[max(0, i - 2): i + 3]) for i in range(len(cur))])  # light smooth
    sample_curves[s] = sm

# apply each sample's own curve
parts = []
for s in samples:
    sub = df.filter(pl.col("sample") == s)
    exp = np.interp(sub["gc"].to_numpy(), centers, sample_curves[s])
    parts.append(sub.with_columns((pl.col("rel") / pl.Series("e", exp)).alias("gc_rel")))
df = pl.concat(parts)
# pooled correction too, for the side-by-side comparison
df = df.with_columns((pl.col("rel") / pl.Series("c", curve(df["gc"].to_numpy()))).alias("gc_rel_pooled"))
df.write_parquet(f"benchmarks/tomato2/results/window_cov.w{W}.gcnorm.parquet")


def madn(x):  # normalised MAD (robust sd)
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def residual_by_sample(col):
    """per sample: max |median(col) - 1| over GC 0.22-0.45 (= leftover GC trend)."""
    wn = df.filter((pl.col("breadth") > 0.8) & (pl.col("gc") > 0.22) & (pl.col("gc") < 0.45))
    out = {}
    for s in samples:
        d = wn.filter(pl.col("sample") == s)
        g = d["gc"].to_numpy(); r = d[col].to_numpy()
        b = np.digitize(g, edges) - 1
        meds = [np.median(r[b == i]) for i in range(len(centers)) if (b == i).sum() >= 20]
        out[s] = float(np.max(np.abs(np.array(meds) - 1.0))) if meds else np.nan
    return out


sc = df.filter((pl.col("breadth") > 0.8) & (pl.col("rel") > 0.4) & (pl.col("rel") < 1.6))
print(f"single-copy mode robust-SD: raw={madn(sc['rel'].to_numpy()):.3f}  "
      f"pooled-corrected={madn(sc['gc_rel_pooled'].to_numpy()):.3f}  "
      f"PER-SAMPLE={madn(sc['gc_rel'].to_numpy()):.3f}")

dev_pool = np.array(list(residual_by_sample("gc_rel_pooled").values()))
dev_ps = np.array(list(residual_by_sample("gc_rel").values()))
print(f"\nleftover GC trend per sample (max |median-1|, GC 0.22-0.45):")
print(f"  POOLED curve     : median={np.nanmedian(dev_pool):.3f}  p90={np.nanpercentile(dev_pool,90):.3f}")
print(f"  PER-SAMPLE curve : median={np.nanmedian(dev_ps):.3f}  p90={np.nanpercentile(dev_ps,90):.3f}")
print("  => per-sample flattens the residual the pooled curve leaves behind")

# ---- figure: shapes differ (A) -> per-sample correction flattens them (B) -> mode (C) ----
fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
for s in samples:
    ax[0].plot(centers, raw_curves[s], lw=0.7, alpha=0.6)
ax[0].plot(centers, pooled, color="k", lw=2, label="pooled")
ax[0].axhline(1, color="grey", ls=":"); ax[0].legend(fontsize=8)
ax[0].set_title("A. raw GC curves, one per sample\n(they fan out => shapes differ)")
ax[0].set_xlabel("GC"); ax[0].set_ylabel("rel depth"); ax[0].set_xlim(0.2, 0.48); ax[0].set_ylim(0.4, 1.6)

wn = df.filter((pl.col("breadth") > 0.8))
for s in samples:
    d = wn.filter(pl.col("sample") == s)
    g = d["gc"].to_numpy(); r = d["gc_rel"].to_numpy()
    b = np.digitize(g, edges) - 1
    m = [np.median(r[b == i]) if (b == i).sum() >= 20 else np.nan for i in range(len(centers))]
    ax[1].plot(centers, m, lw=0.7, alpha=0.6)
ax[1].axhline(1, color="grey", ls=":")
ax[1].set_title("B. residual after PER-SAMPLE correction\n(flat at 1 => fits)")
ax[1].set_xlabel("GC"); ax[1].set_ylabel("gc_rel median"); ax[1].set_xlim(0.22, 0.45); ax[1].set_ylim(0.7, 1.3)

ax[2].hist(np.clip(sc["rel"].to_numpy(), 0, 2), bins=100, alpha=0.45, label="raw", color="grey")
ax[2].hist(np.clip(sc["gc_rel"].to_numpy(), 0, 2), bins=100, alpha=0.55, label="per-sample GC-corrected", color="steelblue")
ax[2].axvline(1, color="k", ls=":"); ax[2].legend(); ax[2].set_title("C. single-copy mode: before vs after")
ax[2].set_xlabel("relative coverage")
fig.tight_layout()
fig.savefig(f"benchmarks/tomato2/results/gc_normalization.w{W}.png", dpi=130, bbox_inches="tight")
print(f"\nwrote results/window_cov.w{W}.gcnorm.parquet + gc_normalization.w{W}.png")
