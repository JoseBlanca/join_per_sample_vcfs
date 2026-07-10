#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# ///
"""Build a tidy per-locus / per-coverage evaluation table: for each catalog locus
and each coverage, what the GIAB truth says (length-variant or not) and what each
caller (ours, HipSTR) called. The marimo dashboard reads this table.

Universe = the 13,272 restricted-catalog loci (both callers genotype exactly these).
Keyed by (chrom, pos) where pos = catalog_start + 1 (== our VCF POS == HipSTR START).

Truth (GIAB TR benchmark VCF) records sit at arbitrary tract positions, so truth is
matched to a catalog locus by OVERLAP. A locus is truth length-variant when some
overlapping truth record has HG002 carrying an allele whose length differs from its
REF by a whole number of repeat units (a repeat-unit indel — what an STR caller
targets). `truth_var_any` (any indel, incl. off-period) is kept as a looser variant.

Call state per caller: `var` (emitted, PASS, GT carries a nonzero length delta),
`ref` (called but all-zero delta, or ours dropped it as monomorphic), `nocall`
(GT ./. or FILTER!=PASS). Ours dropping a locus = judged non-variant => `ref`.

Usage: build_eval_table.py  (paths are relative to ssr_hg002/)
Writes results/eval_table.tsv
"""
import bisect, gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "catalog/HG002_Tier_restricted.cat"
TRUTH = ROOT / "truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz"
OURS = lambda c: ROOT / f"results/ours/vcf/HG002_{c}x.ssr.vcf"
HIP = lambda c: ROOT / f"results/hipstr/HG002_{c}x.str.vcf.gz"
FB = lambda c: ROOT / f"results/freebayes/HG002_{c}x.fb.vcf.gz"
FB_MINQUAL = 20.0   # freebayes low-QUAL tail; gate for a fair precision comparison
COVERAGES = [300, 50, 30, 20, 15, 10, 5]
OUT = ROOT / "results/eval_table.tsv"


def _open(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith((".gz", ".bgz", ".cat")) else open(p)


def gt_indices(gt):
    parts = gt.replace("|", "/").split("/")
    if any(x in (".", "") for x in parts):
        return None
    return [int(x) for x in parts]


# ---- catalog loci: (chrom,pos) -> (start0, end0, period) --------------------
loci = {}
with _open(CATALOG) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        chrom, s, e, motif = f[0], int(f[1]), int(f[2]), f[3]
        loci[(chrom, s + 1)] = {"start": s, "end": e, "period": len(motif)}
print(f"catalog loci: {len(loci)}")


# ---- generalized variant-interval index (truth AND freebayes) ---------------
# Both are position-anchored callers (a variant at an arbitrary tract position,
# not our catalog key), so they match a locus by OVERLAP. GIAB/freebayes anchor
# an indel at the base BEFORE the affected sequence, so we extend each interval
# by the longest allele so the inserted/affected bases overlap the tract.
MAX_BACK = 300  # bounded back-scan; covers any STR indel span in this benchmark


def variant_intervals(path, min_qual=None):
    """path -> {chrom: (sorted [(s0, s0+span, deltas)], [starts])}.
    HG002 is the single (last) sample column. `min_qual` gates on QUAL."""
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
            chrom, pos, ref, alt, qual, fmt = f[0], int(f[1]), f[3], f[4], f[5], f[8]
            if min_qual is not None:
                try:
                    if float(qual) < min_qual:
                        continue
                except ValueError:
                    pass
            alleles = [ref] + ([] if alt == "." else alt.split(","))
            gi = fmt.split(":").index("GT")
            idx = gt_indices(f[9].split(":")[gi])
            if idx is None:
                continue
            deltas = tuple(sorted(len(alleles[i]) - len(ref) for i in idx))
            s0 = pos - 1
            span = max([len(ref)] + [len(a) for a in alleles[1:]])
            iv.setdefault(chrom, []).append((s0, s0 + span, deltas))
    out = {}
    for c, v in iv.items():
        vs = sorted(v)
        out[c] = (vs, [x[0] for x in vs])
    return out


def overlap_state(index, chrom, s, e, period):
    """(var_period, var_any): does an indexed length-changing variant overlap
    locus [s,e)? period-aware for var_period."""
    entry = index.get(chrom)
    if not entry:
        return False, False
    v, st = entry
    j = bisect.bisect_right(st, e) - 1
    var_p = var_any = False
    while j >= 0 and st[j] >= s - MAX_BACK:            # bounded back-scan
        ts, te, deltas = v[j]
        if ts < e and te > s:                          # overlap
            for d in deltas:
                if d != 0:
                    var_any = True
                    if period and d % period == 0:
                        var_p = True
        j -= 1
    return var_p, var_any


truth_idx = variant_intervals(TRUTH)


# ---- caller VCF -> {(chrom,pos): ("var"|"ref"|"nocall")} ---------------------
def load_calls(path):
    out = {}
    if not Path(path).exists():
        return out
    with _open(path) as fh:
        samples = []
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, ref, alt, flt, fmt = f[0], int(f[1]), f[3], f[4], f[6], f[8]
            alleles = [ref] + ([] if alt == "." else alt.split(","))
            gi = fmt.split(":").index("GT")
            idx = gt_indices(f[9].split(":")[gi])
            called_pass = flt in (".", "PASS", "")
            if idx is None or not called_pass:
                out[(chrom, pos)] = "nocall"
                continue
            deltas = [len(alleles[i]) - len(ref) for i in idx]
            out[(chrom, pos)] = "var" if any(d != 0 for d in deltas) else "ref"
    return out


# ---- emit the tidy table ----------------------------------------------------
rows = []
for cov in COVERAGES:
    ours = load_calls(OURS(cov))
    hip = load_calls(HIP(cov))
    # freebayes: position-anchored indels (not STR-aware) -> overlap-match like truth,
    # QUAL-gated (it emits a long low-QUAL tail the STR callers don't).
    fb_idx = variant_intervals(FB(cov), min_qual=FB_MINQUAL)
    for (chrom, pos), meta in loci.items():
        s, e, per = meta["start"], meta["end"], meta["period"]
        tvp, tva = overlap_state(truth_idx, chrom, s, e, per)
        _, fva = overlap_state(fb_idx, chrom, s, e, per)   # caller-var = any length change
        # ours drops monomorphic loci -> absence == called non-variant ("ref").
        os = ours.get((chrom, pos), "ref")
        hs = hip.get((chrom, pos), "nocall")   # HipSTR absence == no genotype
        fs = "var" if fva else "ref"           # freebayes absence == no variant found
        rows.append(
            (cov, chrom, pos, per, int(tvp), int(tva), os, hs, fs)
        )

with open(OUT, "w") as out:
    out.write("cov\tchrom\tpos\tperiod\ttruth_var\ttruth_var_any\tours\thipstr\tfreebayes\n")
    for r in rows:
        out.write("\t".join(map(str, r)) + "\n")
print(f"wrote {len(rows)} rows -> {OUT}")
