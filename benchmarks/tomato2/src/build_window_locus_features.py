#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy", "cyvcf2"]
# ///
"""Option (c) payoff: attach each variant locus to its 2kb window's per-sample
GC-corrected coverage (gc_rel) and recompute the carrier-style counts with the
power windowing gives (threshold ~1.5x, which per-base coverage couldn't support).

Window<->genotype alignment: VCF column i == sorted-SRR psp i (verified, 0
mismatches), so the window matrix (columns sorted by SRR) aligns positionally
with cyvcf2's gt arrays.

Per locus (W=2000):
  n_excess_win          # samples whose WINDOW gc_rel > T (default 1.5)
  n_excess_carrier_win  # of those that also carry the ALT (het/hom-alt)
  mean_win_gcrel        cohort-mean window copy number (1.0 = single copy)
  max_win_gcrel         most-amplified sample's window copy number
Output: results/locus_window_features.parquet  (+ comparison print vs per-base)
"""
import argparse
import numpy as np
import polars as pl
from cyvcf2 import VCF

W = 2000
VCF_PATH = "benchmarks/tomato2/results/cohort.vcf.gz"
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3


def bsnp(v):
    return len(v.REF) == 1 and len(v.ALT) == 1 and len(v.ALT[0]) == 1


def main(T):
    gc = pl.read_parquet(f"benchmarks/tomato2/results/window_cov.w{W}.gcnorm.parquet")
    wide = gc.pivot(values="gc_rel", index=["chrom", "win_start"], on="sample")
    srr_cols = sorted(c for c in wide.columns if c not in ("chrom", "win_start"))
    mat = wide.select(srr_cols).to_numpy()  # (n_win, n_samples) in sorted-SRR order
    keys = list(zip(wide["chrom"].to_list(), wide["win_start"].to_list()))
    window_dict = {k: mat[i] for i, k in enumerate(keys)}
    print(f"{len(window_dict)} windows x {len(srr_cols)} samples")

    vcf = VCF(VCF_PATH)
    assert len(vcf.samples) == len(srr_cols), "sample count mismatch"

    rows = []
    for v in vcf:
        if not bsnp(v):
            continue
        win_start = ((v.POS - 1) // W) * W + 1
        arr = window_dict.get((v.CHROM, win_start))
        gt = v.gt_types
        carrier = (gt == HET) | (gt == HOM_ALT)
        if arr is None:
            rows.append((v.CHROM, int(v.POS), win_start, 0, 0, np.nan, np.nan, np.nan))
            continue
        covered = np.isfinite(arr)
        excess = covered & (arr > T)
        n_excess_win = int(excess.sum())
        n_excess_carrier_win = int((excess & carrier).sum())
        mean_gcrel = float(np.nanmean(arr)) if covered.any() else np.nan
        max_gcrel = float(np.nanmax(arr)) if covered.any() else np.nan
        # windowed coverage GAP: are the het samples the over-covered ones?
        # (GC-corrected, window-level analog of per-base cov_het_gap; no threshold)
        het_m = (gt == HET) & covered
        hom_m = ((gt == HOM_REF) | (gt == HOM_ALT)) & covered
        win_gap = (float(np.mean(arr[het_m])) - float(np.mean(arr[hom_m]))) \
            if (het_m.any() and hom_m.any()) else np.nan
        rows.append((v.CHROM, int(v.POS), win_start, n_excess_win,
                     n_excess_carrier_win, mean_gcrel, max_gcrel, win_gap))
    vcf.close()

    df = pl.DataFrame(rows, schema=["chrom", "pos", "win_start", "n_excess_win",
                                    "n_excess_carrier_win", "mean_win_gcrel", "max_win_gcrel",
                                    "win_cov_het_gap"],
                      orient="row")
    df.write_parquet(f"benchmarks/tomato2/results/locus_window_features.w{W}.parquet")
    print(f"wrote locus_window_features.w{W}.parquet ({df.height} loci, T={T}, W={W})")

    # ---- comparison vs per-base, on the het axis ----
    feat = pl.read_parquet("benchmarks/tomato2/results/snp_features.parquet")
    j = feat.join(df, on=["chrom", "pos"], how="inner")
    print(f"\nobs_het by WINDOWED carrier count (vs per-base in parentheses):")
    print("  the test: does the windowed count light up moderate-het loci the per-base missed?")
    b = (j.with_columns([
            pl.col("n_excess_carrier_win").clip(0, 8).alias("ncw"),
        ]).group_by("ncw").agg(
            pl.len().alias("n"),
            pl.col("obs_het").mean().round(3).alias("obs_het"),
            pl.col("n_excess_carrier").mean().round(2).alias("perbase_carrier"),
            pl.col("mean_win_gcrel").median().round(2).alias("win_cov"),
            pl.col("qual").median().round(0).alias("qual"),
        ).sort("ncw"))
    print(b)

    print("\nMODERATE-HET capture (obs_het 0.2-0.5, the class per-base missed):")
    mod = j.filter((pl.col("obs_het") >= 0.2) & (pl.col("obs_het") <= 0.5))
    print(f"  loci: {mod.height}")
    for lab, expr in [("per-base  n_excess_carrier>=1", pl.col("n_excess_carrier") >= 1),
                      ("windowed  n_excess_carrier_win>=1", pl.col("n_excess_carrier_win") >= 1),
                      ("windowed  mean_win_gcrel>1.3", pl.col("mean_win_gcrel") > 1.3),
                      ("windowed  win_cov_het_gap>0.3", pl.col("win_cov_het_gap") > 0.3)]:
        c = mod.filter(expr).height
        print(f"  flagged by {lab:36s}: {c:5d} ({100*c/mod.height:.0f}%)")

    # windowed coverage GAP vs per-base gap: obs_het by gap bin
    print("\nwindowed coverage GAP (win_cov_het_gap) — obs_het by gap bin (per-base gap for ref):")
    g = (j.with_columns(pl.col("win_cov_het_gap").cut(
            [-0.1, 0.1, 0.3, 0.6, 1.0], labels=["<-.1", "-.1-.1", ".1-.3", ".3-.6", ".6-1", ">1"]
         ).alias("gapbin"))
         .group_by("gapbin").agg(
            pl.len().alias("n"),
            pl.col("obs_het").mean().round(3).alias("obs_het"),
            pl.col("cov_het_gap").drop_nans().mean().round(3).alias("perbase_gap"),
            pl.col("mean_win_gcrel").median().round(2).alias("win_cov"),
            pl.col("qual").median().round(0).alias("qual"),
         ).sort("gapbin"))
    print(g)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=1.5)
    ap.add_argument("--w", type=int, default=2000)
    a = ap.parse_args()
    W = a.w  # override module-global window size
    main(a.threshold)
