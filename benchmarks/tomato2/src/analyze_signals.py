#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy", "cyvcf2"]
# ///
"""Probe two questions:
  (1) what does n_excess actually flag (how extreme a sample), and how clean a
      discriminator is it?
  (2) does mqdiff add information *conditional on* n_excess, and what is the
      high-het / mqdiff~0 class?
Plus a candidate new signal: n_excess_carrier = samples that are BOTH over-covered
AND carry the ALT allele (paralog carriers), built from the VCF."""
import polars as pl
import numpy as np
from cyvcf2 import VCF

VCF_PATH = "benchmarks/tomato2/results/cohort.vcf.gz"
ALPHA = 0.01
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3

df = pl.read_parquet("benchmarks/tomato2/results/snp_features.parquet")

# (1) what coverage multiple does n_excess flag?  cutoff / scale
print("=== what n_excess flags ===")
print("per-sample cutoff = 99th pct of coverage; scale = median coverage")
print("=> a sample is 'excess' only above ~cutoff/scale x its own normal\n")

# (2) mqdiff conditional on n_excess (cross-tab)
print("=== mqdiff x n_excess cross-tab (mean obs_het / mean cov / med qual / n) ===")
ct = (df
    .with_columns([
        pl.when(pl.col("n_excess") >= 3).then(pl.lit("n_excess>=3"))
          .when(pl.col("n_excess") >= 1).then(pl.lit("n_excess 1-2"))
          .otherwise(pl.lit("n_excess 0")).alias("ne_cls"),
        pl.when(pl.col("mqdiff") <= -5).then(pl.lit("mqdiff<=-5 (divergent)"))
          .when(pl.col("mqdiff") < -1).then(pl.lit("mqdiff -5..-1"))
          .otherwise(pl.lit("mqdiff ~0")).alias("mq_cls"),
    ])
    .group_by(["ne_cls", "mq_cls"]).agg(
        pl.len().alias("n"),
        pl.col("obs_het").mean().round(3).alias("obs_het"),
        pl.col("mean_rel_cov").median().round(2).alias("cov"),
        pl.col("qual").median().round(0).alias("med_qual"),
    ).sort(["ne_cls", "mq_cls"]))
with pl.Config(tbl_rows=20):
    print(ct)

# the high-het / mqdiff~0 class: nearly-identical duplications?
print("\n=== high-het loci (obs_het>0.5) split by mqdiff ===")
hh = df.filter(pl.col("obs_het") > 0.5)
for lab, d in [("mqdiff ~0 (>-2)", hh.filter(pl.col("mqdiff") > -2)),
               ("mqdiff <=-5", hh.filter(pl.col("mqdiff") <= -5))]:
    print(f"  {lab:18s} n={d.height:5d}  cov={d['mean_rel_cov'].median():.2f}  "
          f"n_excess={d['n_excess'].median():.0f}  qual={d['qual'].median():.0f}  "
          f"bal={d['het_alt_balance'].drop_nans().median():.2f}")

# (3) candidate new signal: genotype-aware carrier count from the VCF
def bsnp(v):
    return len(v.REF) == 1 and len(v.ALT) == 1 and len(v.ALT[0]) == 1

vcf = VCF(VCF_PATH)
cols = [v.format("DP")[:, 0].astype(np.int32) for v in vcf if bsnp(v)]
vcf.close()
mat = np.array(cols)
cutoff = np.array([np.quantile(c[c > 0], 1 - ALPHA) if (c > 0).any() else np.nan
                   for c in mat.T])

rows = []
vcf = VCF(VCF_PATH)
for v in vcf:
    if not bsnp(v):
        continue
    gt = v.gt_types
    dp = v.format("DP")[:, 0].astype(np.float64)
    covered = dp > 0
    if covered.sum() < 30:
        continue
    excess = covered & (dp > cutoff)
    het = gt == HET
    altc = (gt == HET) | (gt == HOM_ALT)
    n_excess = int(excess.sum())
    n_excess_carrier = int((excess & altc).sum())   # over-covered AND carries ALT
    n_called = int(((gt == HOM_REF) | het | (gt == HOM_ALT)).sum())
    rows.append((n_excess, n_excess_carrier, int(het.sum()) / max(n_called, 1)))
vcf.close()
g = pl.DataFrame(rows, schema=["n_excess", "n_excess_carrier", "obs_het"], orient="row")

print("\n=== genotype-aware carrier count: of the over-covered samples, are they ALT carriers? ===")
tot_ex = g["n_excess"].sum()
tot_carrier = g["n_excess_carrier"].sum()
print(f"  total excess-sample observations: {tot_ex:,}")
print(f"  of which carry ALT:               {tot_carrier:,}  ({100*tot_carrier/tot_ex:.0f}%)")
print("\n  obs_het by n_excess_carrier bin:")
gb = (g.with_columns(pl.col("n_excess_carrier").clip(0, 8).alias("nc"))
       .group_by("nc").agg(pl.len().alias("n"), pl.col("obs_het").mean().round(3).alias("obs_het"))
       .sort("nc"))
print(gb)
