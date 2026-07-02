#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy", "scipy", "cyvcf2"]
# ///
"""Python companion to the Rust `examples/paralog_score_parity.rs` (plan R1).

Rust has no parquet reader, so the loose Python<->Rust LR correlation (the
porting sanity-check) is done here: read the Rust TSVs + the prototype's
paralog_lr.parquet, join on (chrom, pos), and report the correlation. Also
report the EB prior pi computed on each LR vector, and the Spearman between the
Rust per-sample inbreeding F and the prototype's cohort F recomputed from the
VCF.

Usage (host, uv):
  benchmarks/tomato2/src/paralog_score_parity.py \
      --rust-lr      tmp/paralog_parity_lr.tsv \
      --rust-samples tmp/paralog_parity_samples.tsv \
      --proto-lr     benchmarks/tomato2/results/paralog_lr.parquet \
      --vcf          benchmarks/tomato2/results/cohort.vcf.gz
"""
import argparse
import numpy as np
import polars as pl
from scipy.special import expit, logit
from scipy.stats import spearmanr, pearsonr

ap = argparse.ArgumentParser()
ap.add_argument("--rust-lr", default="tmp/paralog_parity_lr.tsv")
ap.add_argument("--rust-samples", default="tmp/paralog_parity_samples.tsv")
ap.add_argument("--proto-lr", default="benchmarks/tomato2/results/paralog_lr.parquet")
ap.add_argument("--features", default="benchmarks/tomato2/results/snp_features.parquet")
ap.add_argument("--vcf", default="benchmarks/tomato2/results/cohort.vcf.gz")
args = ap.parse_args()

HOM_REF, HET, HOM_ALT = 0, 1, 3


def em_pi(lr, start=0.03):
    pi = start
    for _ in range(500):
        pi_new = float(np.mean(expit(lr + logit(pi))))
        if abs(pi_new - pi) < 1e-9:
            return pi_new
        pi = pi_new
    return pi


# ---- LR correlation: Rust vs prototype, on shared loci ----
rust = pl.read_csv(args.rust_lr, separator="\t")
proto = pl.read_parquet(args.proto_lr).select(["chrom", "pos", "lr"]).rename({"lr": "lr_proto"})
joined = rust.rename({"lr": "lr_rust"}).join(proto, on=["chrom", "pos"], how="inner")
lr_r = joined["lr_rust"].to_numpy()
lr_p = joined["lr_proto"].to_numpy()

print(f"loci: rust={rust.height}  proto={proto.height}  shared={joined.height}")
finite = np.isfinite(lr_r) & np.isfinite(lr_p)
pear = pearsonr(lr_r[finite], lr_p[finite]).statistic
spear = spearmanr(lr_r[finite], lr_p[finite]).statistic
print(f"LR correlation on shared loci: Pearson={pear:.4f}  Spearman={spear:.4f}")
# Also clipped (the extreme veto tails dominate raw Pearson).
clip = lambda x: np.clip(x, -40, 60)
pear_c = pearsonr(clip(lr_r[finite]), clip(lr_p[finite])).statistic
print(f"  Pearson on clip(-40,60): {pear_c:.4f}")

print(f"prior pi:  rust LRs={em_pi(lr_r[finite]):.4f}  proto LRs={em_pi(lr_p[finite]):.4f}")

# ---- Flagged-set paralog profile: elevated het + coverage excess vs kept ----
try:
    feat = pl.read_parquet(args.features).select(
        ["chrom", "pos", "obs_het", "fis", "mean_rel_cov"])
    prof = rust.join(feat, on=["chrom", "pos"], how="inner")
    print(f"\nFDR<5% flagged set profile (Rust qval, {prof.height} loci w/ features):")
    for lab, d in [("flagged (qval<.05)", prof.filter(pl.col("qval") < 0.05)),
                   ("kept    (qval>=.05)", prof.filter(pl.col("qval") >= 0.05))]:
        print(f"  {lab:20s} n={d.height:7d}  "
              f"obs_het={d['obs_het'].mean():.3f}  "
              f"fis={d['fis'].drop_nans().mean():+.2f}  "
              f"mean_rel_cov={d['mean_rel_cov'].drop_nans().median():.2f}")
except Exception as e:  # noqa: BLE001
    print(f"\n(skipped flagged-set profile: {e})")

# ---- F Spearman: Rust per-sample F vs the prototype's cohort F ----
# Prototype F (build_paralog_lr.py): per-sample het deficit on clean loci
# (cohort obs_het <= 0.2), F = clip(1 - obs/exp, 0, 0.99).
from cyvcf2 import VCF  # noqa: E402

vcf = VCF(args.vcf)
samples = list(vcf.samples)
n = len(samples)
obs = np.zeros(n)
exp = np.zeros(n)
for v in vcf:
    if len(v.REF) != 1 or len(v.ALT) != 1 or len(v.ALT[0]) != 1:
        continue
    gt = v.gt_types
    called = (gt == HOM_REF) | (gt == HET) | (gt == HOM_ALT)
    nc = int(called.sum())
    if nc < 30:
        continue
    nh = int((gt == HET).sum())
    if nh / nc > 0.2:
        continue
    p = (2 * int((gt == HOM_ALT).sum()) + nh) / (2 * nc)
    obs[called] += (gt[called] == HET)
    exp[called] += 2 * p * (1 - p)
proto_f = np.clip(1 - obs / np.maximum(exp, 1e-9), 0.0, 0.99)
proto_f_by_sample = dict(zip(samples, proto_f))

rust_s = pl.read_csv(args.rust_samples, separator="\t")
common = [(row["F"], proto_f_by_sample[row["sample"]])
          for row in rust_s.iter_rows(named=True) if row["sample"] in proto_f_by_sample]
rf = np.array([a for a, _ in common])
pf = np.array([b for _, b in common])
fs = spearmanr(rf, pf).statistic
print(f"\nF: {len(common)} shared samples;  Spearman(rust F, proto F) = {fs:.3f}")
print(f"  rust F  median={np.median(rf):.2f} range=[{rf.min():.2f},{rf.max():.2f}]")
print(f"  proto F median={np.median(pf):.2f} range=[{pf.min():.2f},{pf.max():.2f}]")
