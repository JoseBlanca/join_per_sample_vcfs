#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["cyvcf2", "numpy", "scipy"]
# ///
"""Is n_excess real low-AF paralog signal, or over-dispersion noise?

Decisive per-sample test: at each locus, the Poisson test flags a set of
"excess-coverage" samples. If n_excess is real paralogy, those *same* samples
should carry the ALT allele (het/hom-alt) far more often than the other covered
samples at the locus. If it's noise, the excess samples are random and carry ALT
at the locus's baseline rate."""
import numpy as np
from cyvcf2 import VCF
from scipy.stats import poisson

VCF_PATH = "benchmarks/tomato2/results/cohort.vcf.gz"
ALPHA = 0.01
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3


def is_bsnp(v):
    return len(v.REF) == 1 and len(v.ALT) == 1 and len(v.ALT[0]) == 1


# pass 1: per-sample single-copy coverage scale
vcf = VCF(VCF_PATH)
cols = []
for v in vcf:
    if is_bsnp(v):
        cols.append(v.format("DP")[:, 0].astype(np.int32))
vcf.close()
mat = np.array(cols)
scale = np.array([np.median(c[c > 0]) if (c > 0).any() else np.nan
                  for c in mat.T])
print(f"scale median={np.nanmedian(scale):.1f}")

# pass 2: pool per-sample outcomes split by excess vs non-excess
# buckets: [all sites] and [low-AF sites: mean_rel_cov < 1.2]
acc = {k: dict(ex_n=0, ex_alt=0, ex_het=0, ne_n=0, ne_alt=0, ne_het=0)
       for k in ("all", "lowAF")}

vcf = VCF(VCF_PATH)
for v in vcf:
    if not is_bsnp(v):
        continue
    gt = v.gt_types
    dp = v.format("DP")[:, 0].astype(np.float64)
    covered = dp > 0
    if covered.sum() < 30:
        continue
    rel = dp / scale
    median_rel = np.median(rel[covered])
    expected = scale[covered] * median_rel
    pen = poisson.sf(dp[covered] - 1, expected) < ALPHA  # excess among covered
    cov_idx = np.where(covered)[0]
    excess_idx = cov_idx[pen]
    nonexc_idx = cov_idx[~pen]
    if excess_idx.size == 0:
        continue

    gt_c = gt
    alt = (gt_c == HET) | (gt_c == HOM_ALT)
    het = gt_c == HET
    mean_rel_cov = float(np.mean(rel[covered]))
    keys = ["all"] + (["lowAF"] if mean_rel_cov < 1.2 else [])
    for k in keys:
        a = acc[k]
        a["ex_n"] += excess_idx.size
        a["ex_alt"] += int(alt[excess_idx].sum())
        a["ex_het"] += int(het[excess_idx].sum())
        a["ne_n"] += nonexc_idx.size
        a["ne_alt"] += int(alt[nonexc_idx].sum())
        a["ne_het"] += int(het[nonexc_idx].sum())
vcf.close()

for k, a in acc.items():
    print(f"\n=== {k} sites (with >=1 excess sample) ===")
    ex_alt = a["ex_alt"] / a["ex_n"] if a["ex_n"] else float("nan")
    ne_alt = a["ne_alt"] / a["ne_n"] if a["ne_n"] else float("nan")
    ex_het = a["ex_het"] / a["ex_n"] if a["ex_n"] else float("nan")
    ne_het = a["ne_het"] / a["ne_n"] if a["ne_n"] else float("nan")
    print(f"  excess samples     n={a['ex_n']:>9d}  ALT-carrier={ex_alt:.3f}  het={ex_het:.3f}")
    print(f"  non-excess samples n={a['ne_n']:>9d}  ALT-carrier={ne_alt:.3f}  het={ne_het:.3f}")
    if ne_alt:
        print(f"  ALT-carrier enrichment (excess / non-excess): {ex_alt/ne_alt:.2f}x")
        print(f"  het         enrichment (excess / non-excess): {ex_het/ne_het:.2f}x")
