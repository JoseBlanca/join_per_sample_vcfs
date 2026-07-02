#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy"]
# ///
"""Quick first-look: are the three axes real and separable, and how do they
relate to obs het (the FP proxy)? Prints, no plots."""
import polars as pl
import numpy as np

df = pl.read_parquet("benchmarks/tomato2/results/snp_features.parquet")
print(f"sites: {df.height}\n")

# --- 1. marginal distributions of the axes ---
def q(col):
    s = df[col].drop_nans().drop_nulls()
    qs = np.quantile(s.to_numpy(), [0, .5, .9, .99, .999, 1])
    return "  ".join(f"{x:8.3f}" for x in qs)
print("quantiles  [min   median   p90    p99    p99.9   max]")
for c in ["mean_rel_cov", "frac_excess", "obs_het", "fis",
          "cov_het_gap", "corr_het_cov", "mqdiff", "mqdifft", "het_alt_balance"]:
    print(f"  {c:16s} {q(c)}")

# --- 2. correlation among axes + response ---
print("\nPearson correlation matrix (key axes vs obs_het / fis):")
cols = ["mean_rel_cov", "frac_excess", "cov_het_gap", "corr_het_cov",
        "mqdiff", "mqdifft", "het_alt_balance", "obs_het", "fis"]
sub = df.select(cols).fill_nan(None).drop_nulls()
m = sub.to_numpy()
C = np.corrcoef(m, rowvar=False)
hdr = "                 " + "".join(f"{c[:8]:>9s}" for c in cols)
print(hdr)
for i, c in enumerate(cols):
    print(f"  {c:14s} " + "".join(f"{C[i,j]:9.2f}" for j in range(len(cols))))

# --- 3. obs_het binned by depth-excess ---
print("\nobs_het and fis by mean_rel_cov bin:")
binned = (df
    .with_columns(pl.col("mean_rel_cov").cut(
        [1.1, 1.25, 1.5, 2.0, 3.0],
        labels=["<1.1", "1.1-1.25", "1.25-1.5", "1.5-2", "2-3", ">3"]).alias("relbin"))
    .group_by("relbin")
    .agg(pl.len().alias("n"),
         pl.col("obs_het").mean().alias("obs_het"),
         pl.col("fis").mean().alias("fis"),
         pl.col("mqdiff").mean().alias("mqdiff"),
         pl.col("cov_het_gap").mean().alias("dhgap"),
         pl.col("qual").median().alias("med_qual"))
    .sort("relbin"))
print(binned)

# --- 4. the extreme depth-excess sites (top paralog candidates) ---
print("\nTop 15 depth-excess sites (mean_rel_cov desc):")
top = (df.sort("mean_rel_cov", descending=True)
         .head(15)
         .select(["chrom","pos","mean_rel_cov","frac_excess","obs_het","fis",
                  "cov_het_gap","corr_het_cov","mqdiff","het_alt_balance","qual","af"]))
with pl.Config(tbl_cols=-1, tbl_width_chars=200, float_precision=2):
    print(top)

# --- 5. cross-tab: combined paralog flag vs het ---
print("\nCombined-signal cross-tab (rel>1.5 AND fis<-0.2 = paralog-like):")
df2 = df.with_columns(
    ((pl.col("mean_rel_cov") > 1.5) & (pl.col("fis") < -0.2)).alias("paralog_like"))
print(df2.group_by("paralog_like").agg(
    pl.len().alias("n"),
    pl.col("obs_het").mean().alias("obs_het"),
    pl.col("mean_rel_cov").mean().alias("rel_depth"),
    pl.col("mqdiff").mean().alias("mqdiff"),
    pl.col("qual").median().alias("med_qual")))
