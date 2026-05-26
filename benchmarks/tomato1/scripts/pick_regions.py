#!/usr/bin/env python3
# Pick random non-overlapping regions from the tomato SL4.0 reference and
# emit a BED file for slicing CRAMs and constraining variant callers.
#
# Defaults: 20 regions x 100 kb on chr01..chr12 (skips unplaced chr00).
# Seeded for reproducibility — re-running with the same seed yields the
# same BED, so collaborators can regenerate the slices without committing
# the BED itself.
#
# Usage:
#   ./pick_regions.py \
#       --fai /home/jose/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa.fai \
#       --out tmp/benchmark_cohort/regions.bed
#
# Optional: --n-regions, --region-size, --edge-buffer, --min-gap, --seed.

import argparse
import random
import sys
from pathlib import Path


def parse_fai(path: Path, keep_prefix: str = "SL4.0ch", skip_names=("SL4.0ch00",)):
    """Yield (chrom, length) for chromosomes matching keep_prefix and not in skip_names."""
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        chrom, length = fields[0], int(fields[1])
        if not chrom.startswith(keep_prefix) or chrom in skip_names:
            continue
        yield chrom, length


def allocate_counts(chroms, total_regions):
    """Allocate region counts per chromosome proportional to length."""
    total_len = sum(L for _, L in chroms)
    raw = [(c, L, total_regions * L / total_len) for c, L in chroms]
    # floor + distribute remainder by largest fractional part
    floors = [(c, L, int(r), r - int(r)) for c, L, r in raw]
    assigned = sum(f for _, _, f, _ in floors)
    remainder = total_regions - assigned
    floors.sort(key=lambda t: t[3], reverse=True)
    counts = {}
    for i, (c, L, f, _) in enumerate(floors):
        counts[c] = f + (1 if i < remainder else 0)
    # preserve original chrom order in the returned mapping
    return {c: counts[c] for c, _ in chroms}


def pick_chrom_regions(chrom_len, n, region_size, edge_buffer, min_gap, rng,
                       max_attempts_factor=200):
    """Pick n non-overlapping (>= min_gap apart) region starts on one chromosome."""
    lo = edge_buffer
    hi = chrom_len - edge_buffer - region_size
    if hi <= lo or n <= 0:
        return []
    picked = []
    attempts = 0
    cap = max_attempts_factor * max(n, 1)
    while len(picked) < n and attempts < cap:
        attempts += 1
        start = rng.randrange(lo, hi + 1)
        ok = True
        for s in picked:
            if abs(s - start) < (region_size + min_gap):
                ok = False
                break
        if ok:
            picked.append(start)
    if len(picked) < n:
        print(
            f"  warn: only placed {len(picked)}/{n} regions after {attempts} tries",
            file=sys.stderr,
        )
    picked.sort()
    return picked


def main():
    ap = argparse.ArgumentParser(
        description="Pick random non-overlapping regions and emit a BED file."
    )
    ap.add_argument("--fai", type=Path, required=True,
                    help="reference .fai (S_lycopersicum_chromosomes.4.00.fa.fai)")
    ap.add_argument("--out", type=Path, required=True, help="output BED path")
    ap.add_argument("--n-regions", type=int, default=20)
    ap.add_argument("--region-size", type=int, default=100_000)
    ap.add_argument("--edge-buffer", type=int, default=50_000,
                    help="min distance from chromosome ends")
    ap.add_argument("--min-gap", type=int, default=500_000,
                    help="min gap between picked regions on the same chromosome")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    chroms = list(parse_fai(args.fai))
    if not chroms:
        sys.exit(f"no usable chromosomes in {args.fai}")
    counts = allocate_counts(chroms, args.n_regions)

    print(f"# seed={args.seed} n_regions={args.n_regions} "
          f"region_size={args.region_size}", file=sys.stderr)
    bed_lines = []
    tsv_lines = ["chrom\tstart\tend\tlength"]
    for chrom, length in chroms:
        n = counts[chrom]
        if n == 0:
            print(f"  {chrom}: 0 regions", file=sys.stderr)
            continue
        starts = pick_chrom_regions(
            length, n, args.region_size, args.edge_buffer, args.min_gap, rng,
        )
        print(f"  {chrom}: {len(starts)} region(s)", file=sys.stderr)
        for s in starts:
            e = s + args.region_size
            bed_lines.append(f"{chrom}\t{s}\t{e}")
            tsv_lines.append(f"{chrom}\t{s}\t{e}\t{args.region_size}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(bed_lines) + "\n")
    tsv_path = args.out.with_suffix(".tsv")
    tsv_path.write_text("\n".join(tsv_lines) + "\n")

    total_bp = len(bed_lines) * args.region_size
    print(f"wrote {len(bed_lines)} regions ({total_bp/1e6:.2f} Mb) -> {args.out}",
          file=sys.stderr)
    print(f"wrote summary -> {tsv_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
