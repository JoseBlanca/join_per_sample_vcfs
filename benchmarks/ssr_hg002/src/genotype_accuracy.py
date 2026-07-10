#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# ///
"""Genotype-LENGTH accuracy (not just detection): among truth-positive loci a caller
detects, does it get HG002's EXACT allele lengths right? This is where an STR-aware
caller should earn its keep — freebayes can flag "a length changed" without nailing
the repeat count.

Common currency = the genotype as a sorted pair of per-allele bp deltas relative to
the reference genome:
  truth : GIAB is phased; per haplotype, sum (len(ALT_on_that_hap) - len(REF)) over
          the overlapping records -> (hap0_bp, hap1_bp).
  ours  : sorted(len(ALT_i) - len(REF_tract)) over GT  (REF tract == reference).
  HipSTR: the GB field IS bp-diff-from-reference per allele.
  freebayes: sorted bp deltas from the best-overlapping indel record's GT (APPROX —
          freebayes is usually unphased and can split a tract across records).

Reported per coverage, over loci that are truth length-variant (period-aware) AND
have a well-defined phased truth genotype:
  detected      = caller called the locus a length variant (the detection recall)
  gt_correct    = detected AND exact genotype match
  acc | detected = gt_correct / detected  (given you detected it, did you nail it?)

Usage: genotype_accuracy.py   (from ssr_hg002/)
"""
import bisect, gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "catalog/HG002_Tier_restricted.cat"
TRUTH = ROOT / "truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz"
OURS = lambda c: ROOT / f"results/ours/vcf/HG002_{c}x.ssr.vcf"
HIP = lambda c: ROOT / f"results/hipstr/HG002_{c}x.str.vcf.gz"
FB = lambda c: ROOT / f"results/freebayes/HG002_{c}x.fb.vcf.gz"
FB_MINQUAL = 20.0
COVERAGES = [300, 50, 30, 20, 15, 10, 5]
MAX_BACK = 300


def _open(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith((".gz", ".bgz", ".cat")) else open(p)


def phased(gt):
    if "|" in gt:
        parts = gt.split("|")
    else:
        parts = gt.split("/")
    if any(x in (".", "") for x in parts):
        return None
    return [int(x) for x in parts], ("|" in gt)


# ---- catalog loci -----------------------------------------------------------
loci = {}
with _open(CATALOG) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        loci[(f[0], int(f[1]) + 1)] = {"start": int(f[1]), "end": int(f[2]), "period": len(f[3])}


# ---- generic overlap index of (s0, e0, record) ------------------------------
def build_index(path, min_qual=None):
    iv = {}
    if not Path(path).exists():
        return iv
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            if min_qual is not None:
                try:
                    if float(f[5]) < min_qual:
                        continue
                except ValueError:
                    pass
            ref, alt = f[3], f[4]
            alleles = [ref] + ([] if alt == "." else alt.split(","))
            s0 = int(f[1]) - 1
            span = max([len(ref)] + [len(a) for a in alleles[1:]])
            iv.setdefault(f[0], []).append((s0, s0 + span, ref, alleles, f[8], f[9]))
    return {c: (sorted(v), [x[0] for x in sorted(v)]) for c, v in iv.items()}


def overlapping(index, chrom, s, e):
    entry = index.get(chrom)
    if not entry:
        return []
    v, st = entry
    j = bisect.bisect_right(st, e) - 1
    out = []
    while j >= 0 and st[j] >= s - MAX_BACK:
        if v[j][0] < e and v[j][1] > s:
            out.append(v[j])
        j -= 1
    return out


def gfield(fmt, sample, key):
    parts = fmt.split(":")
    if key not in parts:
        return None
    return sample.split(":")[parts.index(key)]


# ---- truth phased genotype per locus (bp deltas) ----------------------------
truth_idx = build_index(TRUTH)


def truth_gt(chrom, s, e, period):
    """sorted (hap0_bp, hap1_bp) or None; also whether it's a period-aware variant."""
    hap = [0, 0]
    seen = False
    ok = True
    for s0, e0, ref, alleles, fmt, sample in overlapping(truth_idx, chrom, s, e):
        g = phased(gfield(fmt, sample, "GT"))
        if g is None:
            continue
        idx, is_phased = g
        if not is_phased or len(idx) != 2:
            ok = False
            continue
        for h in (0, 1):
            hap[h] += len(alleles[idx[h]]) - len(ref)
        seen = True
    if not seen or not ok:
        return None
    gt = tuple(sorted(hap))
    is_var = any(d != 0 and (not period or d % period == 0) for d in gt)
    return gt if is_var else None


# ---- caller genotype per locus (bp deltas) ----------------------------------
def ours_hip_gt(path, is_hipstr):
    """(chrom,pos) -> sorted tuple of bp deltas, or None if no-call."""
    out = {}
    if not Path(path).exists():
        return out
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[6] not in (".", "PASS", ""):
                continue
            g = phased(gfield(f[8], f[9], "GT"))
            if g is None:
                continue
            idx = g[0]
            if is_hipstr:
                gb = gfield(f[8], f[9], "GB")
                if gb and all(x.lstrip("-").isdigit() for x in gb.replace("|", "/").split("/")):
                    out[(f[0], int(f[1]))] = tuple(sorted(int(x) for x in gb.replace("|", "/").split("/")))
                    continue
            alleles = [f[3]] + ([] if f[4] == "." else f[4].split(","))
            out[(f[0], int(f[1]))] = tuple(sorted(len(alleles[i]) - len(f[3]) for i in idx))
    return out


def fb_gt(index, chrom, s, e):
    """APPROX freebayes genotype: bp deltas from the overlapping indel record with the
    largest length change that is in the GT."""
    best = None
    for s0, e0, ref, alleles, fmt, sample in overlapping(index, chrom, s, e):
        g = phased(gfield(fmt, sample, "GT"))
        if g is None:
            continue
        idx = g[0]
        deltas = tuple(sorted(len(alleles[i]) - len(ref) for i in idx))
        if all(d == 0 for d in deltas):
            continue
        mag = max(abs(d) for d in deltas)
        if best is None or mag > best[0]:
            best = (mag, deltas)
    return best[1] if best else None


# ---- compare ----------------------------------------------------------------
print(f"{'cov':>4} {'caller':>9} {'truthpos':>8} {'detected':>8} {'gt_correct':>10} "
      f"{'recall':>7} {'acc|det':>8}")
for cov in COVERAGES:
    ours = ours_hip_gt(OURS(cov), False)
    hip = ours_hip_gt(HIP(cov), True)
    fbi = build_index(FB(cov), min_qual=FB_MINQUAL)
    tot = 0
    stat = {c: [0, 0] for c in ("ours", "hipstr", "freebayes")}  # [detected, correct]
    for (chrom, pos), m in loci.items():
        tg = truth_gt(chrom, m["start"], m["end"], m["period"])
        if tg is None:
            continue
        tot += 1
        for name, g in (("ours", ours.get((chrom, pos))),
                        ("hipstr", hip.get((chrom, pos))),
                        ("freebayes", fb_gt(fbi, chrom, m["start"], m["end"]))):
            if g is not None and any(d != 0 for d in g):
                stat[name][0] += 1
                if g == tg:
                    stat[name][1] += 1
    for name in ("ours", "hipstr", "freebayes"):
        det, cor = stat[name]
        rec = cor / tot if tot else 0
        acc = cor / det if det else 0
        print(f"{cov:>4} {name:>9} {tot:>8} {det:>8} {cor:>10} {rec:>7.3f} {acc:>8.3f}")
