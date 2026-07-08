# /// script
# requires-python = ">=3.10"
# ///
"""Ours-vs-HipSTR SSR loss analysis on the ssr_tomato1 rerun.

Reproduces the dashboard's §1b loss buckets on the two cohort VCFs and adds
finer breakdowns of the "genuine drop" set: motif length, dinucleotide vs.
other, interior interruption, motif composition, tract length. VCF-only — no
catalog needed. Run:

  uv run --no-project benchmarks/ssr_tomato1/scripts/loss_analysis.py \
    --ours  <ours cohort.ssr.vcf> \
    --hipstr <hipstr cohort.str.vcf.gz>
"""
import argparse
import bisect
import gzip
from collections import Counter


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith((".gz", ".bgz")) else open(p)


def _info_int(info, key):
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            try:
                return int(kv.split("=")[1])
            except ValueError:
                return None
    return None


def _gt(gt):
    parts = gt.replace("|", "/").split("/")
    if any(p in (".", "") for p in parts):
        return None
    return [int(p) for p in parts]


def parse_vcf(path):
    samples, loci = [], {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, flt, info, fmt = (
                f[0], int(f[1]), f[3], f[4], f[6], f[7], f[8])
            alleles = [ref] + ([] if alt == "." else alt.split(","))
            rlen = len(ref)
            period = _info_int(info, "PERIOD")
            gi = fmt.split(":").index("GT")
            rel = {}
            for name, col in zip(samples, f[9:]):
                idx = _gt(col.split(":")[gi])
                rel[name] = None if idx is None else tuple(
                    sorted(len(alleles[i]) - rlen for i in idx))
            flag = "PASS" if flt in (".", "PASS", "") else flt
            loci[(chrom, pos)] = {
                "period": period, "filter": flag, "rel": rel,
                "ref": ref, "reflen": rlen}
    return samples, loci


def build_pass_index(loci):
    idx = {}
    for (ch, pos), rec in loci.items():
        if rec["filter"] != "PASS":
            continue
        if not any(g is not None for g in rec["rel"].values()):
            continue
        idx.setdefault(ch, []).append((pos, pos + rec["reflen"]))
    for ch in idx:
        idx[ch].sort()
    return idx


def overlaps(idx, ch, s, e):
    v = idx.get(ch)
    if not v:
        return False
    starts = [x[0] for x in v]
    j = bisect.bisect_right(starts, e) - 1
    while j >= 0 and v[j][0] >= s - 500:
        if v[j][0] < e and v[j][1] > s:
            return True
        j -= 1
    return False


def has_interruption(ref, period):
    if not period or period < 1 or len(ref) < 2 * period:
        return False
    motif = ref[:period]
    return any(b != motif[i % period] for i, b in enumerate(ref))


def is_hip_length_variable(rec):
    p = rec["period"]
    if not p:
        return False
    for g in rec["rel"].values():
        if g is not None and any(d != 0 and d % p == 0 for d in g):
            return True
    return False


def canon_motif(ref, period):
    """Canonical rotation of the first-frame motif (lowercased, min rotation)."""
    if not period or period < 1 or len(ref) < period:
        return None
    m = ref[:period].upper()
    rots = [m[i:] + m[:i] for i in range(len(m))]
    return min(rots)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True)
    ap.add_argument("--hipstr", required=True)
    args = ap.parse_args()

    _, ours = parse_vcf(args.ours)
    _, hip = parse_vcf(args.hipstr)
    pass_idx = build_pass_index(ours)

    buckets = Counter()
    drops = []  # hip recs that are genuine drops
    for k, hrec in hip.items():
        if not is_hip_length_variable(hrec):
            continue
        ch, pos = k
        orec = ours.get(k)
        if orec is not None and orec["filter"] == "PASS" and any(
                g is not None for g in orec["rel"].values()):
            buckets["called (exact)"] += 1
        elif orec is not None and orec["filter"] != "PASS":
            buckets[f"filtered ({orec['filter']})"] += 1
        elif orec is not None:
            buckets["no-call (lowDepth row)"] += 1
        elif overlaps(pass_idx, ch, pos, pos + hrec["reflen"]):
            buckets["called (boundary-shifted)"] += 1
        else:
            buckets["genuine drop"] += 1
            drops.append((k, hrec))

    n_var = sum(buckets.values())
    print(f"HipSTR length-variable loci (comparable): {n_var}")
    print("\n=== Loss buckets ===")
    for lab, v in buckets.most_common():
        print(f"  {lab:32s} {v:5d}  ({v/n_var:5.1%})")

    # --- genuine-drop breakdowns ---
    n = len(drops)
    print(f"\n=== Genuine drops: {n} ===")

    by_period = Counter(h["period"] for _, h in drops)
    print("\nby motif length (period):")
    for p in sorted(by_period):
        print(f"  period {p}: {by_period[p]:4d}  ({by_period[p]/n:5.1%})")

    interr = sum(1 for _, h in drops if has_interruption(h["ref"], h["period"]))
    print(f"\ninterior interruption (impure tract): {interr}  ({interr/n:.1%})")
    print(f"pure tract                           : {n-interr}  ({(n-interr)/n:.1%})")

    print("\ntop canonical motifs among drops:")
    mot = Counter(canon_motif(h["ref"], h["period"]) for _, h in drops)
    for m, c in mot.most_common(15):
        print(f"  {str(m):8s} {c:4d}  ({c/n:5.1%})")

    # tract length distribution
    tl = [h["reflen"] for _, h in drops]
    tl.sort()
    print("\nREF tract length (bp): "
          f"min {tl[0]}  p25 {tl[len(tl)//4]}  median {tl[len(tl)//2]}  "
          f"p75 {tl[3*len(tl)//4]}  max {tl[-1]}")

    # how many drops have HipSTR variation in ONLY a few samples (rare) vs many?
    def n_var_samples(h):
        p = h["period"]
        return sum(1 for g in h["rel"].values()
                   if g is not None and any(d != 0 and d % p == 0 for d in g))
    nv = [n_var_samples(h) for _, h in drops]
    singleton = sum(1 for x in nv if x == 1)
    print(f"\nHipSTR calls the length-variant in exactly 1 sample: "
          f"{singleton}  ({singleton/n:.1%})  "
          f"(median samples-variant among drops: {sorted(nv)[len(nv)//2]})")

    # cross: dinucleotide + interruption
    di = [(k, h) for k, h in drops if h["period"] == 2]
    di_int = sum(1 for _, h in di if has_interruption(h["ref"], 2))
    print(f"\ndinucleotide drops: {len(di)}  of which interrupted: {di_int} "
          f"({di_int/max(len(di),1):.1%})")


if __name__ == "__main__":
    main()
