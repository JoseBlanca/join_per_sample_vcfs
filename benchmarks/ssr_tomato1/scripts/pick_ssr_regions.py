# /// script
# requires-python = ">=3.10"
# ///
"""Pick an SSR-targeted slice BED from our ssr-catalog.

The SNP `tomato1` regions.bed (80 x 100 kb windows) is a sparse, incidental
sampling of SSR space (~774 SSRs in 8 Mb). For the SSR-vs-HipSTR benchmark we
want crams whose reads sit *on* SSRs, so a fixed data budget covers many more
genotypable loci.

This samples N catalog SSRs (seeded, reproducible), expands each tract by
`--flank` bp on each side (so reads can span the repeat + flanks), merges
overlapping windows, and writes a sorted BED. Feed that BED to
slice_crams_on_rick.sh; the resulting crams carry reads for ~N SSRs in roughly
`N * (2*flank + tract)` bp — e.g. 8000 SSRs at +/-500 bp ~= 8 Mb, the same data
budget as the current SNP slice but ~10x the SSRs.

Coordinates: the catalog is 0-based half-open [start, end); the BED is the
same convention, so windows are written verbatim (clipped at 0). Selection is
uniform over catalog loci by default; --per-period stratifies equally across
motif-period classes so rarer long-period SSRs are not swamped.

Usage:
  pick_ssr_regions.py --catalog c.ssr.catalog --n 8000 --flank 500 \
      --seed 42 --out ssr_regions.bed [--per-period]
"""

import argparse
import gzip
import random
import sys
from collections import defaultdict


def read_catalog(path, exclude=()):
    """-> list of (chrom, start, end, period), skipping excluded chroms."""
    exclude = set(exclude)
    loci = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, start, end, motif = f[0], int(f[1]), int(f[2]), f[3]
            if chrom in exclude:
                continue
            loci.append((chrom, start, end, len(motif)))
    return loci


def merge_windows(wins):
    """wins: list of (chrom, start, end) -> sorted, overlaps/adjacent merged."""
    by_chrom = defaultdict(list)
    for chrom, s, e in wins:
        by_chrom[chrom].append((s, e))
    out = []
    for chrom in sorted(by_chrom):
        iv = sorted(by_chrom[chrom])
        cs, ce = iv[0]
        for s, e in iv[1:]:
            if s <= ce:  # overlap or touch
                ce = max(ce, e)
            else:
                out.append((chrom, cs, ce))
                cs, ce = s, e
        out.append((chrom, cs, ce))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--n", type=int, default=8000, help="number of SSRs to sample")
    ap.add_argument("--flank", type=int, default=500, help="bp added each side of every tract")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-period", action="store_true", help="stratify equally across motif periods")
    ap.add_argument("--exclude-chrom", default="", help="comma-separated chroms to skip (e.g. SL4.0ch00)")
    args = ap.parse_args()

    exclude = [c for c in args.exclude_chrom.split(",") if c]
    loci = read_catalog(args.catalog, exclude=exclude)
    rng = random.Random(args.seed)

    if args.per_period:
        by_period = defaultdict(list)
        for rec in loci:
            by_period[rec[3]].append(rec)
        periods = sorted(by_period)
        per = max(1, args.n // len(periods))
        chosen = []
        for p in periods:
            pool = by_period[p]
            chosen.extend(rng.sample(pool, min(per, len(pool))))
    else:
        chosen = rng.sample(loci, min(args.n, len(loci)))

    wins = [(c, max(0, s - args.flank), e + args.flank) for (c, s, e, _) in chosen]
    merged = merge_windows(wins)

    with open(args.out, "w") as out:
        for chrom, s, e in merged:
            out.write(f"{chrom}\t{s}\t{e}\n")

    total_bp = sum(e - s for _, s, e in merged)
    period_mix = defaultdict(int)
    for _, _, _, p in chosen:
        period_mix[p] += 1
    print(
        f"catalog SSRs           : {len(loci)}\n"
        f"sampled                : {len(chosen)} (seed {args.seed}, flank +/-{args.flank}"
        + (", per-period" if args.per_period else "") + ")\n"
        f"windows after merge    : {len(merged)}\n"
        f"total slice bp         : {total_bp:,} (~{total_bp / 1e6:.1f} Mb)\n"
        f"period mix (sampled)   : " + ", ".join(f"p{p}={period_mix[p]}" for p in sorted(period_mix)),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
