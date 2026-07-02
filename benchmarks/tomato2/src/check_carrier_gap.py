#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "cyvcf2"]
# ///
"""Why does n_excess_carrier miss moderate-het (0.2-0.5) loci that look paralog-y?
Hypothesis: at ~6x, a 2x paralog carrier has ~12 reads, BELOW the per-sample
'too high' cutoff (~15 reads = 99th pct = 2.5x), so the carriers are het but
never flagged as over-covered. Check it: at each obs_het level, where do the
HET samples' coverages sit, and what fraction clear their own cutoff?"""
import numpy as np
from cyvcf2 import VCF

V = "benchmarks/tomato2/results/cohort.vcf.gz"
ALPHA = 0.01
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3


def bsnp(v):
    return len(v.REF) == 1 and len(v.ALT) == 1 and len(v.ALT[0]) == 1


vcf = VCF(V)
cols = [v.format("DP")[:, 0].astype(np.int32) for v in vcf if bsnp(v)]
vcf.close()
mat = np.array(cols)
scale = np.array([np.median(c[c > 0]) if (c > 0).any() else np.nan for c in mat.T])
cutoff = np.array([np.quantile(c[c > 0], 1 - ALPHA) if (c > 0).any() else np.nan for c in mat.T])
print(f"scale median={np.nanmedian(scale):.1f}  cutoff median={np.nanmedian(cutoff):.1f} "
      f"(= {np.nanmedian(cutoff)/np.nanmedian(scale):.1f}x single-copy)\n")

bins = [(0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 1.01)]
acc = {b: dict(loci=0, het_n=0, het_relsum=0.0, het_flagged=0,
               het_rel_ge15=0, het_rel_ge20=0) for b in bins}

vcf = VCF(V)
for v in vcf:
    if not bsnp(v):
        continue
    gt = v.gt_types
    dp = v.format("DP")[:, 0].astype(np.float64)
    covered = dp > 0
    n_called = int(((gt == HOM_REF) | (gt == HET) | (gt == HOM_ALT)).sum())
    if n_called < 30:
        continue
    het = (gt == HET) & covered
    obs_het = int(het.sum()) / n_called  # note: het here is het&covered; close enough
    b = next((b for b in bins if b[0] <= obs_het < b[1]), None)
    if b is None:
        continue
    a = acc[b]
    a["loci"] += 1
    idx = np.where(het)[0]
    rel = dp[idx] / scale[idx]
    a["het_n"] += idx.size
    a["het_relsum"] += float(rel.sum())
    a["het_flagged"] += int((dp[idx] > cutoff[idx]).sum())
    a["het_rel_ge15"] += int((rel >= 1.5).sum())
    a["het_rel_ge20"] += int((rel >= 2.0).sum())
vcf.close()

print(f"{'obs_het bin':12s} {'loci':>7s} {'mean cov of':>12s} {'% het samp':>11s} "
      f"{'% het >=1.5x':>12s} {'% het >=2x':>11s}")
print(f"{'':12s} {'':>7s} {'het samples':>12s} {'flagged':>11s} {'':>12s} {'':>11s}")
for b in bins:
    a = acc[b]
    if not a["het_n"]:
        continue
    mean_rel = a["het_relsum"] / a["het_n"]
    pct_flag = 100 * a["het_flagged"] / a["het_n"]
    pct15 = 100 * a["het_rel_ge15"] / a["het_n"]
    pct20 = 100 * a["het_rel_ge20"] / a["het_n"]
    print(f"{b[0]:.1f}-{b[1]:<7.2f} {a['loci']:>7d} {mean_rel:>12.2f} "
          f"{pct_flag:>10.1f}% {pct15:>11.1f}% {pct20:>10.1f}%")
