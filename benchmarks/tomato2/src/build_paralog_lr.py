#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy", "cyvcf2", "matplotlib"]
# ///
"""Prototype: per-locus likelihood ratio, hidden-paralog (H2) vs real-variant (H1).

Per sample s we use coverage c_s (500bp per-sample GC-corrected gc_rel) and the
SNP allele reads (k_s alt, n_s total = AD).

H1  real variant, single copy (param: allele freq p)
    genotype G ~ HWE(p);  coverage ~ Normal(1, σ0)  INDEPENDENT of G
    k_s ~ Binom(n_s, vaf), vaf = ε / 0.5 / 1-ε for homref/het/homalt
H2  hidden paralog, dup dosage (param: dup freq q)
    dosage d ∈ {0,1,2} ~ HWE(q);  coverage ~ Normal(1 + d/2, σ_d)
    k_s ~ Binom(n_s, vaf2), vaf2 = ε / 1/3 / 1/2   (tied to the SAME d)
    -> coverage AND allele balance both rise with d, in a fixed ratio.
    -> VAF capped at 0.5: a confident hom-alt sample is ~impossible under H2.

LR = max_q logL_H2 - max_p logL_H1.   LR>0 => paralog, LR<0 => real variant.
Binomial C(n,k) cancels in the ratio (dropped). Params maximized on a grid.
Output: results/paralog_lr.parquet  (+ distribution / spot-check figure)
"""
import numpy as np
import polars as pl
from cyvcf2 import VCF

W = 500
EPS = 0.01
SIGMA0 = 0.26          # single-copy gc_rel spread at 500bp (measured)
MQ_W = 0.25            # mqdiff term: +MQ_W*(-mqdiff) for divergent ALT reads.
                       # coverage-gated by construction (no excess => H2 stays penalised),
                       # so it lifts diverging paralogs (esp. single-carrier) without
                       # flagging introgressions — validated in paralog_simulation.py
CMAX = 4.0             # winsorise coverage for the Normal tail
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3
VCF_PATH = "benchmarks/tomato2/results/cohort.vcf.gz"

# H1 genotype VAFs; H2 dosage VAFs + coverage means/sigmas
VAF_H1 = np.array([EPS, 0.5, 1 - EPS])              # homref, het, homalt
VAF_H2 = np.array([EPS, 1 / 3, 0.5])               # d = 0, 1, 2
MU_H2 = np.array([1.0, 1.5, 2.0])
SIG_H2 = SIGMA0 * np.sqrt(MU_H2)                    # Poisson-ish: SD grows with mean
LOGVAF_H1, LOG1M_H1 = np.log(VAF_H1), np.log(1 - VAF_H1)
LOGVAF_H2, LOG1M_H2 = np.log(VAF_H2), np.log(1 - VAF_H2)

# extended H2: a carrier has total copies T and m mutant copies -> coverage T/2,
# VAF = m/T (high copy => low VAF; m in 1..T-2 so the 2 orthologous copies stay REF)
T_CARRIER = np.array([3, 4, 6, 8])      # coverage 1.5,2,3,4× (cap 4×)
T_MEAN = T_CARRIER / 2.0
T_SIG = SIGMA0 * np.sqrt(T_MEAN)
# per copy level keep only single-PSV (m=1, VAF 1/T) and balanced (m≈T/2, VAF~.5);
# drop the crowded intermediate low-VAF states that aren't separately identifiable
_TM = [(ti, m) for ti, T in enumerate(T_CARRIER)
       for m in sorted({1, int(T // 2)}) if 1 <= m <= T - 2]
TM_TIDX = np.array([ti for ti, _ in _TM])
TM_VC = np.array([m / T_CARRIER[ti] for ti, m in _TM])
LOG_VC, LOG1M_VC = np.log(TM_VC), np.log(1 - TM_VC)
QEXT = np.linspace(0.004, 0.6, 40)

NGRID = 80
PGRID = np.linspace(0.005, 0.995, NGRID)
QGRID = np.linspace(0.005, 0.995, NGRID)


def inbreeding_logprior(freq, Fs):
    """Wright inbreeding genotype/dosage log-priors, per (grid value, sample).
    freq: (NGRID,), Fs: (n_samp,). Returns (NGRID, n_samp, 3) for
    [hom-ref/d0, het/d1, hom-alt/d2] with P(het)=2pq(1-F), homozygotes get +F*pq.
    Tomato is autogamous, so F is large and het is rare under H1; a dup carrier
    in a selfer is usually homozygous for the dup (d=2), not het (d=1)."""
    p = freq[:, None]                      # (NGRID,1) = ALT freq
    q = 1 - p
    F = Fs[None, :]                        # (1,n_samp)
    het = 2 * p * q * (1 - F)
    hom_ref = q * q + F * p * q            # index 0 pairs with VAF≈ε
    hom_alt = p * p + F * p * q            # index 2 pairs with VAF≈1-ε
    return np.log(np.stack([hom_ref, het, hom_alt], axis=2) + 1e-300)


def lse(a, axis):
    m = np.max(a, axis=axis, keepdims=True)
    return (m.squeeze(axis) + np.log(np.sum(np.exp(a - m), axis=axis)))


def bsnp(v):
    return len(v.REF) == 1 and len(v.ALT) == 1 and len(v.ALT[0]) == 1


def normlog(x, mu, sig):
    return -0.5 * ((x - mu) / sig) ** 2 - np.log(sig * np.sqrt(2 * np.pi))


def main():
    # per-sample window coverage matrix, columns sorted by SRR == VCF order
    gc = pl.read_parquet(f"benchmarks/tomato2/results/window_cov.w{W}.gcnorm.parquet")
    wide = gc.pivot(values="gc_rel", index=["chrom", "win_start"], on="sample")
    srr_cols = sorted(c for c in wide.columns if c not in ("chrom", "win_start"))
    mat = wide.select(srr_cols).to_numpy()
    wd = {k: mat[i] for i, k in enumerate(zip(wide["chrom"].to_list(), wide["win_start"].to_list()))}
    n_samp = len(srr_cols)

    # --- per-sample inbreeding F_s from each sample's OWN genome-wide het ---
    # tomato is autogamous and individuals differ (pure inbreds ~0 het, wild
    # relatives higher), so F is per-sample. Estimated on clean loci only
    # (cohort obs_het <= 0.2) so paralog pseudo-hets don't deflate it.
    v0 = VCF(VCF_PATH)
    obs = np.zeros(n_samp)
    exp = np.zeros(n_samp)
    for v in v0:
        if not bsnp(v):
            continue
        gt = v.gt_types
        called = (gt == HOM_REF) | (gt == HET) | (gt == HOM_ALT)
        nc = int(called.sum())
        if nc < 30:
            continue
        nh = int((gt == HET).sum())
        if nh / nc > 0.2:                          # skip paralog-suspect loci
            continue
        p = (2 * int((gt == HOM_ALT).sum()) + nh) / (2 * nc)
        obs[called] += (gt[called] == HET)
        exp[called] += 2 * p * (1 - p)
    v0.close()
    Fs = np.clip(1 - obs / np.maximum(exp, 1e-9), 0.0, 0.99)
    print(f"per-sample inbreeding F_s: median={np.median(Fs):.2f} "
          f"range=[{Fs.min():.2f},{Fs.max():.2f}]  (per-sample het rate, not cohort-wide)")
    LOGPG_full = inbreeding_logprior(PGRID, Fs)        # (NGRID, n_samp, 3)
    # extended-H2 carrier prior P(non-carrier) per (sample, q): inbreeding HWE
    _q = QEXT[None, :]
    _P0 = (1 - _q) ** 2 + Fs[:, None] * _q * (1 - _q)
    LOGP0_FULL = np.log(_P0)                           # (n_samp, nq) non-carrier
    LOGPC_FULL = np.log(np.maximum(1 - _P0, 1e-300))   # (n_samp, nq) carrier

    vcf = VCF(VCF_PATH)
    assert len(vcf.samples) == len(srr_cols)
    rows = []
    for v in vcf:
        if not bsnp(v):
            continue
        arr = wd.get((v.CHROM, ((v.POS - 1) // W) * W + 1))
        if arr is None:
            continue
        ad = v.format("AD")
        ref, alt = ad[:, 0].astype(np.float64), ad[:, 1].astype(np.float64)
        n = ref + alt
        use = np.isfinite(arr) & (ref >= 0) & (alt >= 0) & (n > 0)
        m = int(use.sum())
        if m < 20:
            continue
        c = np.clip(arr[use], 0, CMAX)
        k = alt[use]
        nn = n[use]

        idx = np.where(use)[0]
        # --- H1: real variant; coverage Normal(1,σ0) independent of genotype ---
        lrg = k[:, None] * LOGVAF_H1[None, :] + (nn - k)[:, None] * LOG1M_H1[None, :]
        lc1 = float(np.sum(normlog(c, 1.0, SIGMA0)))
        h1 = LOGPG_full[:, idx, :] + lrg[None, :, :]       # (NGRID,m,3)
        summed1 = np.sum(lse(h1, axis=2), axis=1)          # (NGRID,) per p
        logL1 = lse(summed1, axis=0) - np.log(NGRID) + lc1  # MARGINAL over p (flat prior)

        # --- extended H2: non-carrier (cov 1, hom-ref) vs carrier (copies T,
        #     mutant m -> cov T/2, VAF m/T). max over q (carrier freq) x (T,m) ---
        nbase = normlog(c, 1.0, SIGMA0) + (k * np.log(EPS) + (nn - k) * np.log(1 - EPS))  # (m,)
        cov_car = normlog(c[:, None], T_MEAN[None, :], T_SIG[None, :])                    # (m,5)
        allele_car = k[:, None] * LOG_VC[None, :] + (nn - k)[:, None] * LOG1M_VC[None, :]  # (m,nTM)
        cb = cov_car[:, TM_TIDX] + allele_car                                             # (m,nTM)
        lp0 = LOGP0_FULL[idx]                                                             # (m,nq)
        lpc = LOGPC_FULL[idx]
        A = lp0[:, None, :] + nbase[:, None, None]         # (m,1,nq) non-carrier
        B = lpc[:, None, :] + cb[:, :, None]               # (m,nTM,nq) carrier
        g = np.logaddexp(A, B).sum(axis=0)                 # (nTM,nq) per (config,q)
        logL2 = float(lse(g.ravel(), axis=0) - np.log(g.size))  # MARGINAL over (config,q)

        af = k / nn
        n_homalt_conf = int(np.sum((af > 0.9) & (nn >= 5)))
        mqd = first(v.INFO.get("MQDiff"))                  # per-locus ALT-vs-REF MAPQ diff
        mqd = mqd if np.isfinite(mqd) else 0.0
        lr_val = float(logL2 - logL1) + MQ_W * (-min(mqd, 0.0))
        rows.append((v.CHROM, int(v.POS), lr_val, m,
                     n_homalt_conf, float(np.mean(c)), float(np.mean(af))))
    vcf.close()

    df = pl.DataFrame(rows, schema=["chrom", "pos", "lr", "n_used",
                                    "n_homalt_conf", "mean_c", "mean_af"], orient="row")
    df.write_parquet("benchmarks/tomato2/results/paralog_lr.parquet")
    print(f"{df.height} loci scored")
    print(df["lr"].describe())

    # spot checks vs the per-site features
    feat = pl.read_parquet("benchmarks/tomato2/results/snp_features.parquet").select(
        ["chrom", "pos", "obs_het", "fis", "mean_rel_cov", "qual", "het_alt_balance"])
    j = df.join(feat, on=["chrom", "pos"], how="inner")
    print("\nLR by obs_het bin (paralogs should have HIGH +LR; real variants LOW/-LR):")
    print(j.with_columns(pl.col("obs_het").cut([0.05, 0.1, 0.3, 0.6],
            labels=["<.05", ".05-.1", ".1-.3", ".3-.6", ">.6"]).alias("hb"))
          .group_by("hb").agg(pl.len().alias("n"), pl.col("lr").median().round(2).alias("med_lr"),
                              pl.col("lr").mean().round(2).alias("mean_lr")).sort("hb"))
    print("\nconfident hom-alt present (should push LR negative = real variant):")
    print(j.group_by(pl.col("n_homalt_conf") > 0).agg(
        pl.len().alias("n"), pl.col("lr").median().round(2).alias("med_lr")))
    print("\nknown paralog class (mean_rel_cov>1.5 & fis<-0.2):")
    para = j.filter((pl.col("mean_rel_cov") > 1.5) & (pl.col("fis") < -0.2))
    rest = j.filter(~((pl.col("mean_rel_cov") > 1.5) & (pl.col("fis") < -0.2)))
    print(f"  paralog-class median LR = {para['lr'].median():.2f} (n={para.height})")
    print(f"  rest        median LR = {rest['lr'].median():.2f} (n={rest.height})")

    print("\nmoderate-het (0.2-0.6) by coverage (HWE-H1 was: normal-cov −20.6, excess +74):")
    mod = j.filter((pl.col("obs_het") >= 0.2) & (pl.col("obs_het") <= 0.6))
    for lab, d in [("cov<=1.15", mod.filter(pl.col("mean_rel_cov") <= 1.15)),
                   ("1.15-1.4", mod.filter((pl.col("mean_rel_cov") > 1.15) & (pl.col("mean_rel_cov") <= 1.4))),
                   ("cov>1.4", mod.filter(pl.col("mean_rel_cov") > 1.4))]:
        if d.height:
            print(f"  {lab:10s} n={d.height:5d}  median_LR={d['lr'].median():8.1f}")

    # figure: LR distribution
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    lr = df["lr"].to_numpy()
    ax[0].hist(np.clip(lr, -20, 60), bins=200, color="steelblue")
    ax[0].axvline(0, color="k", ls=":")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("LR = logL(paralog) - logL(variant)")
    ax[0].set_ylabel("loci (log)")
    ax[0].set_title("genome-wide LR distribution")
    hb = j["obs_het"].to_numpy()
    ax[1].scatter(np.clip(j["lr"].to_numpy(), -20, 60), hb, s=3, alpha=0.1,
                  c=np.clip(j["mean_rel_cov"].to_numpy(), 0, 3), cmap="viridis", rasterized=True)
    ax[1].axvline(0, color="k", ls=":")
    ax[1].set_xlabel("LR")
    ax[1].set_ylabel("obs_het")
    ax[1].set_title("LR vs obs_het (colour = mean_rel_cov)")
    fig.tight_layout()
    fig.savefig("benchmarks/tomato2/results/paralog_lr.png", dpi=130, bbox_inches="tight")
    print("\nwrote paralog_lr.parquet + paralog_lr.png")


if __name__ == "__main__":
    main()
