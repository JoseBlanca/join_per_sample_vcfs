#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy", "scipy", "matplotlib"]
# ///
"""Empirical-Bayes calibration of the paralog likelihood ratio.

The LR comes from an explicit two-hypothesis generative model, so (treating
exp(LR) as the Bayes factor — the Occam terms ~cancel: H1 and H2 have the same
parameter count and latent structure) the per-locus posterior is a logistic of
the LR shifted by the prior odds:

    P(paralog | data) = σ( LR + log(π/(1-π)) )

We estimate the prior π (paralog fraction) from the genome-wide LRs by EM:
    E: r_i = σ(LR_i + logit(π));   M: π = mean(r_i);  iterate.

Then the 50%-posterior cut is LR = log((1-π)/π) (positive, NOT 0), and we report
the LR threshold that controls the false-discovery rate (q-value = expected
fraction of flagged loci that are real variants). Writes `post` + `qval` back
into paralog_lr.parquet so the dashboard can use them.
"""
import numpy as np
import polars as pl
from scipy.special import expit, logit

LR_PATH = "benchmarks/tomato2/results/paralog_lr.parquet"

df = pl.read_parquet(LR_PATH)
lr = df["lr"].to_numpy().astype(np.float64)

# --- EM for the prior paralog fraction π ---
pi = 0.03
for it in range(500):
    r = expit(lr + logit(pi))
    pi_new = float(np.mean(r))
    if abs(pi_new - pi) < 1e-9:
        break
    pi = pi_new
print(f"EM converged in {it} iters: paralog fraction π = {pi:.4f} ({100*pi:.2f}% of loci)")

off = logit(pi)
post = expit(lr + off)                    # P(paralog | data)
lr_half = -off                            # LR where posterior = 0.5
print(f"posterior = σ(LR {off:+.2f});  50%% cut at LR = {lr_half:.2f}  (vs the naive 0)")

# --- q-values (tail FDR): sort by posterior desc, running mean of local-fdr ---
lfdr = 1 - post
order = np.argsort(-post)
cum = np.cumsum(lfdr[order]) / np.arange(1, len(post) + 1)
qval = np.empty_like(cum)
qval[order] = cum

df = df.with_columns([pl.Series("post", post), pl.Series("qval", qval)])
df.write_parquet(LR_PATH)

print("\nflagged counts + the LR cut at each FDR level:")
print(f"  {'FDR':>5} {'n_flagged':>10} {'min posterior':>14} {'~LR cut':>9}")
for fdr in (0.01, 0.05, 0.10, 0.20):
    sel = qval <= fdr
    n = int(sel.sum())
    if n:
        minpost = float(post[sel].min())
        lrcut = float(lr[sel].min())
        print(f"  {fdr:>5.0%} {n:>10d} {minpost:>14.3f} {lrcut:>9.1f}")
    else:
        print(f"  {fdr:>5.0%} {0:>10d}")

# --- validate: profile of the FDR<5% flagged set vs the rest ---
feat = pl.read_parquet("benchmarks/tomato2/results/snp_features.parquet").select(
    ["chrom", "pos", "obs_het", "fis", "mean_rel_cov", "qual", "het_alt_balance"])
j = df.join(feat, on=["chrom", "pos"], how="inner")
print("\nFDR<5% flagged set vs rest (held-out features):")
for lab, d in [("flagged (qval<.05)", j.filter(pl.col("qval") < 0.05)),
               ("kept", j.filter(pl.col("qval") >= 0.05))]:
    print(f"  {lab:20s} n={d.height:6d}  obs_het={d['obs_het'].mean():.3f}  "
          f"fis={d['fis'].drop_nans().mean():.2f}  mean_cov={d['mean_rel_cov'].drop_nans().median():.2f}  "
          f"balance={d['het_alt_balance'].drop_nans().median():.2f}  med_qual={d['qual'].median():.0f}")

# --- figure: posterior(LR) sigmoid + FDR curve ---
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
xs = np.linspace(-40, 40, 400)
ax[0].plot(xs, expit(xs + off), color="crimson")
ax[0].axhline(0.5, color="grey", ls=":")
ax[0].axvline(lr_half, color="grey", ls=":")
ax[0].axvline(0, color="k", ls="--", lw=0.8, label="naive LR=0")
ax[0].set_xlabel("LR")
ax[0].set_ylabel("P(paralog | data)")
ax[0].set_title(f"posterior = σ(LR {off:+.1f})   π={pi:.3f}, 50% cut LR={lr_half:.1f}")
ax[0].legend(fontsize=8)

sq = np.sort(qval)
ax[1].plot(np.arange(1, len(sq) + 1), sq, color="steelblue")
for fdr in (0.05, 0.10):
    n = int((qval <= fdr).sum())
    ax[1].axhline(fdr, color="grey", ls=":")
    ax[1].annotate(f"FDR {fdr:.0%}: {n} loci", (n, fdr), fontsize=8, va="bottom")
ax[1].set_xscale("log")
ax[1].set_xlabel("loci ranked by posterior (log)")
ax[1].set_ylabel("q-value (tail FDR)")
ax[1].set_ylim(0, 0.5)
ax[1].set_title("FDR curve")
fig.tight_layout()
fig.savefig("benchmarks/tomato2/results/paralog_eb.png", dpi=130, bbox_inches="tight")
print("\nwrote post+qval into paralog_lr.parquet + paralog_eb.png")
