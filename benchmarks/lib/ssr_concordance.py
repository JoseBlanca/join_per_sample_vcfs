# /// script
# requires-python = ">=3.10"
# ///
"""Ours-vs-HipSTR STR genotype concordance (no truth set).

The two callers genotype the *same* locus set (HipSTR's --regions BED is
projected from our ssr-catalog). With no STR truth for tomato, this scores
*agreement*, not accuracy: at loci/samples both call, do the genotypes match?

Comparison unit: per-allele base-pair difference from REF (the locus's
reference tract). This cancels any constant flank/representation offset
between callers and is exactly what HipSTR's GB field encodes. A genotype
is the sorted multiset of its alleles' (len(allele) - len(REF)) values; two
calls are concordant iff those multisets are equal.

Loci are joined on (CHROM, POS). Both callers put POS at the 1-based repeat
start (== catalog_start + 1), so they align exactly; unmatched loci are
reported. Samples are joined by name (VCF column / @RG SM); if exactly one
sample column exists on each side, they're paired regardless of name.

Usage:
  ssr_concordance.py --ours ours.vcf[.gz] --hipstr hipstr.vcf.gz [--out report.txt]
"""

import argparse
import gzip
import sys
from collections import defaultdict


def opener(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def parse_vcf(path):
    """-> (samples, {(chrom,pos): (ref, [alts], {sample: gt_indices|None})})."""
    samples = []
    loci = {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, fmt = f[0], int(f[1]), f[3], f[4], f[8]
            alts = [] if alt == "." else alt.split(",")
            gt_i = fmt.split(":").index("GT")
            calls = {}
            for name, col in zip(samples, f[9:]):
                gt = col.split(":")[gt_i]
                idx = _gt_indices(gt)
                calls[name] = idx
            loci[(chrom, pos)] = (ref, alts, calls)
    return samples, loci


def _gt_indices(gt):
    """'1/2' or '0|1' -> [1,2]; missing ('.', './.') -> None."""
    parts = gt.replace("|", "/").split("/")
    if any(p in (".", "") for p in parts):
        return None
    return [int(p) for p in parts]


def rel_genotype(ref, alts, idx):
    """Sorted per-allele (len(allele) - len(ref)) for a GT index list."""
    alleles = [ref] + alts
    return tuple(sorted(len(alleles[i]) - len(ref) for i in idx))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ours", required=True)
    ap.add_argument("--hipstr", required=True)
    ap.add_argument("--out", help="write the report here too (default: stdout only)")
    args = ap.parse_args()

    o_samples, ours = parse_vcf(args.ours)
    h_samples, hip = parse_vcf(args.hipstr)

    # Sample pairing: by name, else the lone column on each side.
    shared_samples = [s for s in o_samples if s in set(h_samples)]
    name_map = {s: s for s in shared_samples}
    if not shared_samples and len(o_samples) == 1 and len(h_samples) == 1:
        name_map = {o_samples[0]: h_samples[0]}
        shared_samples = [o_samples[0]]

    shared_loci = sorted(set(ours) & set(hip))
    only_ours = len(set(ours) - set(hip))
    only_hip = len(set(hip) - set(ours))

    # Per-sample-locus tallies over shared loci.
    both_called = concordant = 0
    ours_called = hip_called = 0
    sl_total = len(shared_loci) * len(shared_samples)
    diff_hist = defaultdict(int)  # signed sum-of-allele-length-diff between callers

    for key in shared_loci:
        o_ref, o_alts, o_calls = ours[key]
        h_ref, h_alts, h_calls = hip[key]
        for osamp in shared_samples:
            hsamp = name_map[osamp]
            oi = o_calls.get(osamp)
            hi = h_calls.get(hsamp)
            if oi is not None:
                ours_called += 1
            if hi is not None:
                hip_called += 1
            if oi is None or hi is None:
                continue
            both_called += 1
            og = rel_genotype(o_ref, o_alts, oi)
            hg = rel_genotype(h_ref, h_alts, hi)
            if og == hg:
                concordant += 1
            else:
                diff_hist[sum(og) - sum(hg)] += 1

    def pct(n, d):
        return f"{100 * n / d:.1f}%" if d else "n/a"

    lines = []
    lines.append("=== Ours vs HipSTR — STR genotype concordance (agreement, no truth) ===")
    lines.append(f"ours VCF   : {args.ours}  ({len(o_samples)} samples, {len(ours)} loci)")
    lines.append(f"hipstr VCF : {args.hipstr}  ({len(h_samples)} samples, {len(hip)} loci)")
    lines.append("")
    lines.append(f"shared loci (CHROM,POS join) : {len(shared_loci)}")
    lines.append(f"  loci only in ours          : {only_ours}")
    lines.append(f"  loci only in hipstr        : {only_hip}")
    lines.append(f"paired samples               : {len(shared_samples)}  {shared_samples if len(shared_samples) <= 8 else ''}")
    lines.append("")
    lines.append(f"sample×locus cells (shared)  : {sl_total}")
    lines.append(f"  called by ours             : {ours_called}  ({pct(ours_called, sl_total)})")
    lines.append(f"  called by hipstr           : {hip_called}  ({pct(hip_called, sl_total)})")
    lines.append(f"  called by BOTH             : {both_called}  ({pct(both_called, sl_total)})")
    lines.append("")
    lines.append(f"CONCORDANCE (both-called cells with matching allele-length genotype):")
    lines.append(f"  concordant                 : {concordant} / {both_called}  ({pct(concordant, both_called)})")
    lines.append(f"  discordant                 : {both_called - concordant}")
    if diff_hist:
        lines.append("")
        lines.append("  discordant Σallele-bp-diff (ours − hipstr) histogram:")
        for d in sorted(diff_hist):
            lines.append(f"    {d:+d} bp : {diff_hist[d]}")

    report = "\n".join(lines)
    print(report)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(report + "\n")


if __name__ == "__main__":
    main()
