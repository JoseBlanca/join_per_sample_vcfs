#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["cyvcf2", "numpy", "polars"]
# ///
"""
Phase C: per-site feature table for the paralog-vs-introgression EDA.

Reads the max-retention cohort VCF (biallelic SNPs only for now) and emits one
row per site with the three candidate axes plus the response (observed het / FIS,
our FP proxy). All depth normalisation is per-sample.

Axes (all "coverage" = read depth; per-sample DP normalised by that sample's
single-copy coverage scale = median per-sample DP over biallelic SNP sites,
robust to the thin paralog upper tail; GC-conditioning deferred)
  copy-number : two complementary coverage views, because mean coverage only
                moves when the duplication is at HIGH allele frequency and loses
                LOW-frequency paralogs carried by a few samples:
                  mean_rel_cov  cohort-mean relative coverage (high-AF / fixed
                                collapsed paralogs)
                  n_excess      how many samples have coverage "too high" at the
                                locus vs THEIR OWN coverage distribution: DP above
                                the sample's (1-alpha) empirical quantile cutoff.
                                Counts samples that look deep here, full stop.
                  n_excess_carrier  of those, how many also CARRY the ALT allele
                                (het/hom-alt) = the actual paralog carriers. The
                                sharpest single paralog discriminator we have.
                  max_rel_cov   the single most-excess sample (1-sample paralogs)
  divergence  : INFO MQDiff / MQDiffT (primary ALT).
  co-segregat.: do the *het* samples carry the *excess* coverage? Two ways:
                  - cov_het_gap  = cov(het samples) - cov(hom samples)
                  - corr_het_cov = point-biserial corr(het-indicator, rel cov)

Response (held-out; NOT a predictor here, just the y-axis we correlate against)
  obs_het     : fraction of called samples that are heterozygous
  fis         : 1 - Hobs/Hexp, heterozygote excess vs HWE (negative => excess)

Usage:
  benchmarks/tomato2/src/build_feature_table.py \
      --vcf benchmarks/tomato2/results/cohort.vcf.gz \
      --out benchmarks/tomato2/results/snp_features.parquet
"""
import argparse
import numpy as np
from cyvcf2 import VCF
import polars as pl

# cyvcf2 gt_types encoding
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3


def is_biallelic_snp(v):
    return (
        len(v.REF) == 1
        and len(v.ALT) == 1
        and len(v.ALT[0]) == 1
    )


def sample_coverage_model(vcf_path, alpha):
    """Pass 1: per-sample coverage model from each sample's OWN empirical coverage
    distribution over biallelic SNP sites.
      scale     = median DP (single-copy coverage; used for relative-coverage axes)
      hi_thresh = the (1 - alpha) empirical quantile of DP = the 'coverage too high'
                  cutoff. A sample is flagged 'excess' at a locus when its DP exceeds
                  this, i.e. its coverage there is in its own upper-alpha tail. Using
                  the empirical distribution handles coverage over-dispersion with no
                  Poisson/NB assumption."""
    vcf = VCF(vcf_path)
    n = len(vcf.samples)
    cols = []
    for v in vcf:
        if not is_biallelic_snp(v):
            continue
        cols.append(v.format("DP")[:, 0].astype(np.int32))
    vcf.close()
    mat = np.array(cols)  # (n_sites, n_samples)
    scale = np.empty(n, dtype=np.float64)
    hi_thresh = np.empty(n, dtype=np.float64)
    for j in range(n):
        col = mat[:, j]
        col = col[col > 0]
        if col.size:
            scale[j] = np.median(col)
            hi_thresh[j] = np.quantile(col, 1 - alpha)
        else:
            scale[j] = hi_thresh[j] = np.nan
    return vcf_samples(vcf_path), scale, hi_thresh


def vcf_samples(vcf_path):
    vcf = VCF(vcf_path)
    s = list(vcf.samples)
    vcf.close()
    return s


def first(x):
    """INFO Number=A -> scalar (primary ALT) or nan."""
    if x is None:
        return np.nan
    if isinstance(x, (tuple, list, np.ndarray)):
        return float(x[0])
    return float(x)


def build(vcf_path, out_path, min_called=20, alpha=0.01):
    samples, scale, hi_thresh = sample_coverage_model(vcf_path, alpha)
    n = len(samples)
    print(f"{n} samples; coverage scale median={np.nanmedian(scale):.1f} "
          f"range=[{np.nanmin(scale):.1f},{np.nanmax(scale):.1f}]; "
          f"per-sample 'too high' cutoff = empirical (1-{alpha}) quantile "
          f"(median cutoff={np.nanmedian(hi_thresh):.0f} reads)")

    rows = []
    vcf = VCF(vcf_path)
    kept = skipped = 0
    for v in vcf:
        if not is_biallelic_snp(v):
            continue
        gt = v.gt_types  # (n,)
        dp = v.format("DP")[:, 0].astype(np.float64)
        ad = v.format("AD")  # (n, 2) ref,alt

        called = (gt == HOM_REF) | (gt == HET) | (gt == HOM_ALT)
        n_called = int(called.sum())
        if n_called < min_called:
            skipped += 1
            continue
        het = gt == HET
        homalt = gt == HOM_ALT
        homref = gt == HOM_REF
        n_het = int(het.sum())
        n_homalt = int(homalt.sum())

        AC = 2 * n_homalt + n_het
        AN = 2 * n_called
        p = AC / AN  # alt allele freq
        Hobs = n_het / n_called
        Hexp = 2 * p * (1 - p)
        fis = (1 - Hobs / Hexp) if Hexp > 0 else np.nan

        # per-sample relative coverage (covered samples only)
        rel = dp / scale
        covered = dp > 0
        relc = rel[covered]
        mean_rel = float(np.mean(relc)) if relc.size else np.nan
        median_rel = float(np.median(relc)) if relc.size else np.nan
        max_rel = float(np.max(relc)) if relc.size else np.nan

        # per-sample excess: how many samples have coverage TOO HIGH at this locus
        # vs THEIR OWN coverage distribution. A covered sample is flagged when its
        # DP exceeds its own upper-tail cutoff (the (1-alpha) empirical quantile of
        # its coverage). No cohort/site adjustment -> a locus deep for everyone
        # flags many samples (that is the intended "how many samples look deep").
        excess = covered & (dp > hi_thresh)
        n_cov = int(covered.sum())
        n_excess = int(excess.sum())
        frac_excess = (n_excess / n_cov) if n_cov else np.nan
        # genotype-aware carrier count: over-covered AND carrying the ALT allele
        # (het or hom-alt) = the actual paralog carriers. Sharper paralog signal
        # than n_excess alone, which is mostly deep-but-hom-ref samples.
        n_excess_carrier = int((excess & (het | homalt)).sum())

        # co-segregation: het samples vs hom samples relative coverage
        het_cov = het & covered
        hom_cov = (homref | homalt) & covered
        cov_het = float(np.mean(rel[het_cov])) if het_cov.any() else np.nan
        cov_hom = float(np.mean(rel[hom_cov])) if hom_cov.any() else np.nan
        cov_het_gap = (cov_het - cov_hom) if (het_cov.any() and hom_cov.any()) else np.nan
        # point-biserial corr(het indicator, rel coverage) among called+covered
        cc = called & covered
        corr_het_cov = np.nan
        if cc.sum() >= 5:
            hi = het[cc].astype(np.float64)
            rr = rel[cc]
            if hi.std() > 0 and rr.std() > 0:
                corr_het_cov = float(np.corrcoef(hi, rr)[0, 1])

        # het-sample allele balance (paralog ~0.5)
        bal = np.nan
        if het_cov.any():
            adh = ad[het_cov].astype(np.float64)
            tot = adh[:, 0] + adh[:, 1]
            ok = tot > 0
            if ok.any():
                bal = float(np.mean(adh[ok, 1] / tot[ok]))

        info = v.INFO
        rows.append((
            v.CHROM, int(v.POS), v.REF, v.ALT[0],
            float(v.QUAL) if v.QUAL is not None else np.nan,
            n_called, n_het, n_homalt, AC, AN, p,
            Hobs, Hexp, fis,
            mean_rel, median_rel, n_excess, n_excess_carrier, frac_excess, max_rel,
            cov_het, cov_hom, cov_het_gap, corr_het_cov, bal,
            first(info.get("MQDiff")), first(info.get("MQDiffT")),
            info.get("MQRef") if info.get("MQRef") is not None else np.nan,
            first(info.get("MQAlt")),
            first(info.get("AF")),
            int(info.get("DP")) if info.get("DP") is not None else 0,
        ))
        kept += 1

    vcf.close()
    cols = ["chrom", "pos", "ref", "alt", "qual",
            "n_called", "n_het", "n_homalt", "ac", "an", "p",
            "obs_het", "exp_het", "fis",
            "mean_rel_cov", "median_rel_cov", "n_excess", "n_excess_carrier",
            "frac_excess", "max_rel_cov",
            "cov_het", "cov_hom", "cov_het_gap", "corr_het_cov", "het_alt_balance",
            "mqdiff", "mqdifft", "mqref", "mqalt", "af", "total_cov"]
    df = pl.DataFrame(rows, schema=cols, orient="row")
    df.write_parquet(out_path)
    print(f"kept {kept} sites, skipped {skipped} (<{min_called} called) -> {out_path}")
    print(df.describe())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-called", type=int, default=20)
    ap.add_argument("--alpha", type=float, default=0.01,
                    help="per-sample Poisson upper-tail p threshold for 'coverage too high'")
    a = ap.parse_args()
    build(a.vcf, a.out, a.min_called, a.alpha)
