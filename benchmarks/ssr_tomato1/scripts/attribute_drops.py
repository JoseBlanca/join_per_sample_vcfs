# /// script
# requires-python = ">=3.10"
# ///
"""Attribute the HipSTR-variant loci our SSR caller drops to the mechanism that
killed each one, by joining the genuine-drop set to the PVC_SSR_FATE_TSV sidecar.

  uv run --no-project attribute_drops.py --ours OURS.vcf --hipstr HIP.vcf.gz --fate FATE.tsv
"""
import argparse, bisect, gzip
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
    return None if any(p in (".", "") for p in parts) else [int(p) for p in parts]


def parse_vcf(path):
    s, loci = [], {}
    for line in _open(path):
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            s = line.rstrip("\n").split("\t")[9:]
            continue
        f = line.rstrip("\n").split("\t")
        ref, alt = f[3], f[4]
        alleles = [ref] + ([] if alt == "." else alt.split(","))
        rlen = len(ref)
        per = _info_int(f[7], "PERIOD")
        gi = f[8].split(":").index("GT")
        rel = {}
        for n, c in zip(s, f[9:]):
            idx = _gt(c.split(":")[gi])
            rel[n] = None if idx is None else tuple(sorted(len(alleles[i]) - rlen for i in idx))
        loci[(f[0], int(f[1]))] = {
            "period": per, "filter": ("PASS" if f[6] in (".", "PASS", "") else f[6]),
            "rel": rel, "ref": ref, "reflen": rlen}
    return s, loci


def hip_var(r):
    p = r["period"]
    return bool(p) and any(
        g is not None and any(d != 0 and d % p == 0 for d in g) for g in r["rel"].values())


def has_interruption(ref, period):
    if not period or period < 1 or len(ref) < 2 * period:
        return False
    m = ref[:period]
    return any(b != m[i % period] for i, b in enumerate(ref))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True)
    ap.add_argument("--hipstr", required=True)
    ap.add_argument("--fate", required=True)
    a = ap.parse_args()

    _, ours = parse_vcf(a.ours)
    _, hip = parse_vcf(a.hipstr)

    # ours PASS index for boundary-shift detection
    idx = {}
    for (c, p), r in ours.items():
        if r["filter"] != "PASS" or not any(g is not None for g in r["rel"].values()):
            continue
        idx.setdefault(c, []).append((p, p + r["reflen"]))
    for c in idx:
        idx[c].sort()

    def ours_overlaps(c, s, e):
        v = idx.get(c)
        if not v:
            return False
        starts = [x[0] for x in v]
        j = bisect.bisect_right(starts, e) - 1
        while j >= 0 and v[j][0] >= s - 500:
            if v[j][0] < e and v[j][1] > s:
                return True
            j -= 1
        return False

    # fate TSV: exact map + per-chrom sorted spans for overlap fallback
    fate = {}
    spans = {}
    for line in open(a.fate):
        if line.startswith("#"):
            continue
        c, pos, reflen, n_alt, admit, fk = line.rstrip("\n").split("\t")
        pos, reflen = int(pos), int(reflen)
        fate[(c, pos)] = {"reflen": reflen, "n_alt": int(n_alt), "admit": admit, "fate": fk}
        spans.setdefault(c, []).append((pos, pos + reflen, fk))
    for c in spans:
        spans[c].sort()

    def fate_lookup(c, pos, tlen):
        if (c, pos) in fate:
            return fate[(c, pos)]["fate"]
        # overlap fallback (impure-prefix POS shift)
        v = spans.get(c)
        if not v:
            return "no_catalog_entry"
        starts = [x[0] for x in v]
        j = bisect.bisect_right(starts, pos + tlen) - 1
        while j >= 0 and v[j][0] >= pos - 50:
            if v[j][0] < pos + tlen and v[j][1] > pos:
                return v[j][2]
            j -= 1
        return "no_catalog_entry"

    # genuine drops = HipSTR length-variant, ours emitted no overlapping row
    drops = []
    for k, h in hip.items():
        if not hip_var(h):
            continue
        c, pos = k
        o = ours.get(k)
        if o is not None and o["filter"] == "PASS" and any(g is not None for g in o["rel"].values()):
            continue
        if o is not None:  # filtered or no-call row exists
            continue
        if ours_overlaps(c, pos, pos + h["reflen"]):
            continue
        drops.append((k, h))

    print(f"genuine drops: {len(drops)}\n")
    tally = Counter(fate_lookup(c, p, h["reflen"]) for (c, p), h in drops)
    n = len(drops)
    print("=== fate of the genuine drops ===")
    for fk, v in tally.most_common():
        print(f"  {fk:24s} {v:4d}  ({v/n:5.1%})")

    # break the two dominant drop mechanisms by period + interruption
    print("\n=== drop_no_alt_candidate by period ===")
    sub = [(k, h) for (k, h) in drops if fate_lookup(k[0], k[1], h["reflen"]) == "drop_no_alt_candidate"]
    bp = Counter(h["period"] for _, h in sub)
    for p in sorted(bp):
        print(f"  period {p}: {bp[p]:4d}")
    ii = sum(1 for _, h in sub if has_interruption(h["ref"], h["period"]))
    print(f"  interrupted: {ii}/{len(sub)} ({ii/max(len(sub),1):.1%})")

    print("\n=== drop_em_hom_ref by period ===")
    sub2 = [(k, h) for (k, h) in drops if fate_lookup(k[0], k[1], h["reflen"]) == "drop_em_hom_ref"]
    bp2 = Counter(h["period"] for _, h in sub2)
    for p in sorted(bp2):
        print(f"  period {p}: {bp2[p]:4d}")
    ii2 = sum(1 for _, h in sub2 if has_interruption(h["ref"], h["period"]))
    print(f"  interrupted: {ii2}/{len(sub2)} ({ii2/max(len(sub2),1):.1%})")

    # how many drops are HipSTR singletons per mechanism
    def nvar(h):
        p = h["period"]
        return sum(1 for g in h["rel"].values() if g is not None and any(d != 0 and d % p == 0 for d in g))
    print("\n=== HipSTR-singleton share (variant in exactly 1 sample) per mechanism ===")
    for fk in ["drop_no_alt_candidate", "drop_em_hom_ref", "drop_fp_control"]:
        s = [(k, h) for (k, h) in drops if fate_lookup(k[0], k[1], h["reflen"]) == fk]
        sing = sum(1 for _, h in s if nvar(h) == 1)
        print(f"  {fk:24s}: {sing}/{len(s)} singletons ({sing/max(len(s),1):.1%})")


if __name__ == "__main__":
    main()
