#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# ///
"""Restrict a bgzipped ssr-catalog to loci overlapping a BED, preserving the
header (so the reference_md5 binding is unchanged). Output is bgzipped via the
`bgzip` CLI so the Rust CatalogReader (which always expects BGZF) can read it.

Why: the genome-wide catalog (515k loci) vs a region-restricted BAM makes the
ssr-call pre-pass sample the first 20k *uncovered* catalog loci -> zero confident
genotypes -> failure. Restricting to the benchmark regions fixes that AND is the
set we actually score against.

Usage: restrict_catalog_to_bed.py <catalog.cat(bgzf)> <regions.bed> <out.cat>
"""
import bisect, gzip, subprocess, sys
from pathlib import Path

cat_path, bed_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

# BED intervals per chrom (0-based half-open), sorted by start.
spans: dict[str, list[tuple[int, int]]] = {}
with open(bed_path) as fh:
    for line in fh:
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.rstrip("\n").split("\t")
        spans.setdefault(f[0], []).append((int(f[1]), int(f[2])))
for c in spans:
    spans[c].sort()
starts = {c: [s for s, _ in v] for c, v in spans.items()}

def overlaps(chrom: str, s: int, e: int) -> bool:
    v = spans.get(chrom)
    if not v:
        return False
    st = starts[chrom]
    j = bisect.bisect_right(st, e) - 1          # last region starting <= e
    while j >= 0 and v[j][1] > s:               # region end > locus start
        if v[j][0] < e and v[j][1] > s:
            return True
        j -= 1
    return False

kept = total = 0
op = gzip.open if cat_path.endswith((".gz", ".bgz", ".cat")) else open
raw = out_path + ".raw"
with op(cat_path, "rt") as fin, open(raw, "w") as fout:
    for line in fin:
        if line.startswith("#"):        # header (## meta) + column line (#chrom)
            fout.write(line)
            continue
        total += 1
        f = line.split("\t", 3)
        if overlaps(f[0], int(f[1]), int(f[2])):
            fout.write(line)
            kept += 1

subprocess.run(["bgzip", "-f", raw], check=True)
Path(raw + ".gz").replace(out_path)
print(f"catalog loci: kept {kept} / {total} overlapping {bed_path}")
