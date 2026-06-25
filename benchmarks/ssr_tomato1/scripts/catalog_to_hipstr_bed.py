# /// script
# requires-python = ">=3.10"
# ///
"""Project our `ssr-catalog` (bgzip TSV) into a HipSTR `--regions` BED.

Both callers must genotype the *same* loci for an apples-to-apples
concordance comparison. Our catalog is the single source of truth (one
trf-mod scan); this script re-expresses each catalog locus in HipSTR's
region format so HipSTR genotypes exactly the loci our caller does.

Coordinate conventions (verified against region.cpp + our VCF output):
  * our catalog : 0-based half-open [start, end); `motif` gives the period.
                  (catalog start 6461 -> our VCF POS 6462 == start+1)
  * HipSTR BED  : tab-delimited `CHROM START STOP PERIOD NCOPIES [NAME]`,
                  START 1-based (>=1; HipSTR stores start-1 internally),
                  STOP 1-based-inclusive (== our 0-based exclusive end),
                  PERIOD an int in [1, 9], NCOPIES a double.
So: START = catalog_start + 1 ; STOP = catalog_end ; PERIOD = len(motif).

With `--regions`, only loci whose tract is fully contained in one of the
BED intervals are emitted (the cohort CRAMs are sliced to that BED; loci
outside it would only ever no-call). Loci with period outside HipSTR's
[1, 9] are dropped and counted (none in the tomato catalog, but the SNP
benches favour loud, explicit handling).
"""

import argparse
import bisect
import gzip
import sys
from collections import defaultdict


def load_regions(path):
    """chrom -> sorted list of (start, end) 0-based half-open intervals."""
    regs = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            regs[f[0]].append((int(f[1]), int(f[2])))
    for chrom in regs:
        regs[chrom].sort()
    return regs


def contained(regs, chrom, start, end):
    """Is [start, end) fully inside some interval for `chrom`?"""
    lst = regs.get(chrom)
    if not lst:
        return False
    i = bisect.bisect_right([r[0] for r in lst], start) - 1
    if i < 0:
        return False
    rs, re = lst[i]
    return rs <= start and end <= re


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalog", required=True, help="ssr-catalog bgzip TSV")
    ap.add_argument("--regions", help="restrict to loci inside this BED (the sliced region set)")
    ap.add_argument("--out", required=True, help="output HipSTR regions BED")
    args = ap.parse_args()

    regs = load_regions(args.regions) if args.regions else None

    kept = []  # (chrom, hstart_1based, stop, period, ncopies, name)
    total = dropped_period = dropped_outside = 0

    with gzip.open(args.catalog, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            total += 1
            f = line.rstrip("\n").split("\t")
            chrom, start, end, motif = f[0], int(f[1]), int(f[2]), f[3]
            period = len(motif)
            if not (1 <= period <= 9):
                dropped_period += 1
                continue
            if regs is not None and not contained(regs, chrom, start, end):
                dropped_outside += 1
                continue
            ncopies = (end - start) / period
            name = f"{chrom}_{start + 1}_{motif}"
            kept.append((chrom, start + 1, end, period, ncopies, name))

    kept.sort(key=lambda r: (r[0], r[1]))
    with open(args.out, "w") as out:
        for chrom, hstart, stop, period, ncopies, name in kept:
            out.write(f"{chrom}\t{hstart}\t{stop}\t{period}\t{ncopies:.2f}\t{name}\n")

    print(
        f"catalog loci: {total}\n"
        f"  kept (emitted)     : {len(kept)}\n"
        f"  dropped period>9/<1: {dropped_period}\n"
        f"  dropped outside BED: {dropped_outside}"
        + ("" if regs is not None else "  (no --regions; all in-period loci kept)"),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
