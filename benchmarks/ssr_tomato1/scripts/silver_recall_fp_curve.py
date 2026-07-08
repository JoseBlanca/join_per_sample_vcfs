# /// script
# requires-python = ">=3.10"
# ///
"""Silver-standard recall / false-positive curve for an SSR cohort VCF.

Reproduces the confident-core silver standard from §4 of
`ssr_error_signals_dashboard.py` (read-grounded true100 / false100 label sets, the
confident core = loci where the pileup standard and the HipSTR standard agree), then
scores ANY ours-style cohort VCF at a sweep of QUAL thresholds:

  recall = fraction of the true100 core emitted as variable (PASS, MAP-variable, QUAL>=T)
  FP     = fraction of the false100 core emitted as variable

This is the *shared* scoring contract for comparing emission models (freebayes vs BIC).
The classifier functions are copied verbatim from the dashboard so the core is identical.

  uv run --no-project benchmarks/ssr_tomato1/scripts/silver_recall_fp_curve.py \
    --ours   <cohort.ssr.vcf> \
    --hipstr <hipstr cohort.str.vcf.gz> \
    --reads  <our_reads.tsv> \
    --catalog <ssr_tomato1.ssr.catalog> \
    [--thresholds 0,3,10,20,30,50] [--label NAME]
"""
import argparse
import gzip
from collections import Counter, defaultdict


# ── VCF / catalog parsing ─────────────────────────────────────────────────
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
    return tuple(int(p) for p in parts)


def parse_ours(path):
    """(chrom,pos) -> {period, filter, reflen, qual, cells:{sample: {bpdiff, het}}}."""
    samples, loci = [], {}
    for line in _open(path):
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            continue
        f = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt, qual, flt, info, fmt = (
            f[0], int(f[1]), f[3], f[4], f[5], f[6], f[7], f[8])
        alleles = [ref] + ([] if alt == "." else alt.split(","))
        rlen = len(ref)
        try:
            qv = float(qual)
        except ValueError:
            qv = 0.0
        gi = fmt.split(":").index("GT")
        cells = {}
        for name, col in zip(samples, f[9:]):
            idx = _gt(col.split(":")[gi])
            if idx is None:
                cells[name] = None
                continue
            bp = tuple(sorted(len(alleles[i]) - rlen for i in idx))
            cells[name] = {"bpdiff": bp, "het": len(set(idx)) > 1}
        loci[(chrom, pos)] = {
            "period": _info_int(info, "PERIOD"),
            "filter": ("PASS" if flt in (".", "PASS", "") else flt),
            "reflen": rlen, "qual": qv, "cells": cells}
    return samples, loci


def parse_hipstr(path):
    samples, loci = [], {}
    for line in _open(path):
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            continue
        f = line.rstrip("\n").split("\t")
        chrom, pos, ref, info, fmt = f[0], int(f[1]), f[3], f[7], f[8]
        rlen = len(ref)
        keys = fmt.split(":")
        gi = keys.index("GT")
        idx_of = {k: keys.index(k) for k in ("GB", "PDP") if k in keys}
        cells = {}
        for name, col in zip(samples, f[9:]):
            parts = col.split(":")
            idx = _gt(parts[gi])
            if idx is None:
                cells[name] = None
                continue

            def get(k):
                j = idx_of.get(k)
                return parts[j] if j is not None and j < len(parts) else None

            gb = get("GB")
            bp = None
            if gb and gb not in (".", ""):
                try:
                    bp = tuple(sorted(int(x) for x in gb.replace("|", "/").split("/")))
                except ValueError:
                    bp = None
            het = len(set(idx)) > 1
            mf = None
            pdp = get("PDP")
            if het and pdp:
                try:
                    a, b = (float(x) for x in pdp.split("|"))
                    mf = min(a, b) / (a + b) if (a + b) > 0 else None
                except (ValueError, AttributeError):
                    mf = None
            cells[name] = {"bpdiff": bp, "het": het, "minor_frac": mf}
        loci[(chrom, pos)] = {"period": _info_int(info, "PERIOD"),
                              "reflen": rlen, "cells": cells}
    return samples, loci


# ── silver-standard classifiers (verbatim from dashboard §4) ──────────────
def variation_shape(loc):
    period = loc["period"]
    n_var = 0
    for cell in loc["cells"].values():
        if cell is None or not cell["bpdiff"]:
            continue
        whole_unit = [d for d in cell["bpdiff"]
                      if d != 0 and period and d % period == 0]
        if whole_unit:
            n_var += 1
    return n_var


def pileup_signals(persample, reflen, period, min_depth, hom_frac, bal_min):
    homalt, balhet, carrier, n_assess = Counter(), set(), Counter(), 0
    for lens in persample.values():
        total = sum(lens.values())
        if total < min_depth:
            continue
        wu = {L: c for L, c in lens.items()
              if period and (L - reflen) % period == 0}
        if not wu:
            continue
        n_assess += 1
        top = sorted(wu.items(), key=lambda kv: -kv[1])
        lmod, cmod = top[0]
        dmod, fmod = lmod - reflen, cmod / total
        if dmod != 0 and fmod >= hom_frac:
            homalt[dmod] += 1
            carrier[dmod] += 1
        elif len(top) >= 2:
            l2, c2 = top[1]
            if min(cmod, c2) / (cmod + c2) >= bal_min:
                for d in (dmod, l2 - reflen):
                    if d != 0:
                        balhet.add(d)
                        carrier[d] += 1
    return homalt, balhet, carrier, n_assess


def hip_signals(rec, bal_min):
    period = rec["period"]
    homalt, balhet, carrier, n_assess = Counter(), set(), Counter(), 0
    for cell in rec["cells"].values():
        if cell is None or not cell["bpdiff"]:
            continue
        n_assess += 1
        bp = cell["bpdiff"]
        wu = [d for d in bp if d != 0 and period and d % period == 0]
        for d in set(wu):
            carrier[d] += 1
        if not cell["het"] and wu and len(set(bp)) == 1:
            homalt[bp[0]] += 1
        elif cell["het"] and cell.get("minor_frac") is not None \
                and cell["minor_frac"] >= bal_min and wu:
            for d in wu:
                balhet.add(d)
    return homalt, balhet, carrier, n_assess


def classify_locus(homalt, balhet, carrier, n_assess, recur, min_assessable):
    if n_assess < min_assessable:
        return "undetermined"
    if any(n >= recur for n in homalt.values()) \
            or any(carrier.get(d, 0) >= 2 for d in balhet):
        return "true100"
    if not homalt and not balhet:
        return "false100"
    return "doubtful"


# ── main ──────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True)
    ap.add_argument("--hipstr", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--thresholds", default="0,3,10,20,30,50,75,100")
    ap.add_argument("--label", default="ours")
    # fixed silver-standard thresholds (dashboard defaults)
    ap.add_argument("--min-depth", type=int, default=4)
    ap.add_argument("--hom-frac", type=float, default=0.75)
    args = ap.parse_args()
    bal_min, recur, min_assessable = 0.40, 2, 6

    # read evidence
    reads = defaultdict(lambda: defaultdict(Counter))
    with open(args.reads) as fh:
        next(fh)
        for line in fh:
            s, ch, st, seq, c = line.rstrip("\n").split("\t")
            reads[(ch, int(st))][s][len(seq)] += int(c)

    catalog = {}
    with gzip.open(args.catalog, "rt") as fh:  # catalog is gzip regardless of suffix
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            c0, s0, e0, m0 = p[0], int(p[1]), int(p[2]), p[3]
            catalog[(c0, s0 + 1)] = (e0 - s0, len(m0))

    _, ours = parse_ours(args.ours)
    _, hip = parse_hipstr(args.hipstr)

    pileup_cls = {}
    for k, ps in reads.items():
        meta = catalog.get(k)
        if meta is None:
            continue
        sig = pileup_signals(ps, meta[0], meta[1], args.min_depth,
                             args.hom_frac, bal_min)
        pileup_cls[k] = classify_locus(*sig, recur, min_assessable)
    hip_cls = {k: classify_locus(*hip_signals(r, bal_min), recur, min_assessable)
               for k, r in hip.items()}

    # confident core = both standards agree
    shared = set(pileup_cls) & set(hip_cls)
    core_true = [k for k in shared if pileup_cls[k] == "true100" == hip_cls[k]]
    core_false = [k for k in shared if pileup_cls[k] == "false100" == hip_cls[k]]

    def ours_var(k, thr):
        r = ours.get(k)
        return (r is not None and r["filter"] == "PASS"
                and r["qual"] >= thr and variation_shape(r) > 0)

    # total emission = PASS + MAP-variable records at threshold
    def n_emitted(thr):
        return sum(1 for r in ours.values()
                   if r["filter"] == "PASS" and r["qual"] >= thr
                   and variation_shape(r) > 0)

    print(f"# confident core: {len(core_true)} true100 / {len(core_false)} false100")
    print(f"# scoring VCF: {args.ours}  (label={args.label})")
    print(f"{'QUAL>=':>7} | {'recall':>16} | {'FP rate':>16} | {'emitted':>8}")
    print("-" * 60)
    for thr in [float(x) for x in args.thresholds.split(",")]:
        rec = sum(1 for k in core_true if ours_var(k, thr))
        fp = sum(1 for k in core_false if ours_var(k, thr))
        nt, nf = max(len(core_true), 1), max(len(core_false), 1)
        print(f"{thr:7.0f} | {rec:5d}/{nt:<4d} {rec/nt:6.1%} | "
              f"{fp:5d}/{nf:<5d} {fp/nf:6.2%} | {n_emitted(thr):8d}")


if __name__ == "__main__":
    main()
