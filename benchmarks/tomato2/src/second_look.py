#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy"]
# ///
"""Second look: the separability crux. Among HIGH-HET sites (obs_het>0.3 — where
paralog FPs and legitimate high-het loci both live), does depth cleanly split a
paralog mode (excess depth, REMOVE) from a normal-depth mode (introgression /
balancing-selection-like, KEEP)?"""
import polars as pl
import numpy as np

df = pl.read_parquet("benchmarks/tomato2/results/snp_features.parquet")

def prof(d, name):
    """robust (nan-dropping) profile of a candidate class."""
    n = d.height
    def m(c):
        s = d[c].drop_nans().drop_nulls()
        return float(s.mean()) if s.len() else float("nan")
    def md(c):
        s = d[c].drop_nans().drop_nulls()
        return float(s.median()) if s.len() else float("nan")
    print(f"\n=== {name}: n={n} ({100*n/df.height:.2f}%) ===")
    print(f"  mean_rel_cov  mean={m('mean_rel_cov'):.2f}  med={md('mean_rel_cov'):.2f}")
    print(f"  obs_het         mean={m('obs_het'):.2f}  med={md('obs_het'):.2f}")
    print(f"  fis             mean={m('fis'):.2f}  med={md('fis'):.2f}")
    print(f"  mqdiff          mean={m('mqdiff'):.2f}  med={md('mqdiff'):.2f}")
    print(f"  cov_het_gap   mean={m('cov_het_gap'):.2f}  med={md('cov_het_gap'):.2f}")
    print(f"  corr_het_cov  mean={m('corr_het_cov'):.2f}  med={md('corr_het_cov'):.2f}")
    print(f"  het_alt_balance mean={m('het_alt_balance'):.2f}  med={md('het_alt_balance'):.2f}")
    print(f"  qual            med={md('qual'):.1f}")
    print(f"  af(EM)          med={md('af'):.3f}")

# --- the high-het universe, split by depth ---
hh = df.filter(pl.col("obs_het") > 0.3)
print(f"HIGH-HET sites (obs_het>0.3): {hh.height}  ({100*hh.height/df.height:.2f}% of callset)")

para = hh.filter(pl.col("mean_rel_cov") > 1.5)   # excess depth => paralog
keep = hh.filter(pl.col("mean_rel_cov") <= 1.2)  # normal depth => introgression-like
mid  = hh.filter((pl.col("mean_rel_cov") > 1.2) & (pl.col("mean_rel_cov") <= 1.5))
prof(para, "HIGH-HET + EXCESS depth (>1.5x)  -> paralog (REMOVE)")
prof(mid,  "HIGH-HET + MILD excess (1.2-1.5x) -> ambiguous")
prof(keep, "HIGH-HET + NORMAL depth (<=1.2x) -> introgression/balancing (KEEP)")

# --- is the depth distribution of high-het sites bimodal? ---
print("\nDepth distribution of HIGH-HET sites (mean_rel_cov histogram):")
edges = [0,0.8,1.0,1.2,1.5,2.0,2.5,3.0,5.0,100]
h = hh.select("mean_rel_cov").to_series().to_numpy()
counts, _ = np.histogram(h, bins=edges)
for i,c in enumerate(counts):
    print(f"  {edges[i]:5.1f}-{edges[i+1]:<5.1f} : {c:5d}  {'#'*int(60*c/max(counts))}")

# --- introgression-safety check: do NORMAL-depth divergent sites exist & stay clean? ---
print("\nIntrogression-proxy (normal depth <=1.2x, divergent |mqdiff|>5, het>0.3):")
introg = df.filter((pl.col("mean_rel_cov") <= 1.2) &
                   (pl.col("mqdiff").abs() > 5) &
                   (pl.col("obs_het") > 0.3))
prof(introg, "introgression-proxy")

# --- contrast: what fraction of each depth class is divergent? ---
print("\nDivergence (|mqdiff|>5) rate by depth class among high-het:")
for label, d in [("normal<=1.2", keep), ("excess>1.5", para)]:
    r = d.filter(pl.col("mqdiff").abs() > 5).height / max(d.height,1)
    print(f"  {label:12s}: {100*r:.1f}% divergent")
