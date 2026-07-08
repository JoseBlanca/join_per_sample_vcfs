# /// script
# requires-python = ">=3.10"
# ///
"""Two read-grounded silver standards for the ssr_tomato1 SSR benchmark, and a
scorer for any `ssr-call` VCF against their **confident core**.

There is no gold microsatellite truth for tomato, so we carve out the loci we can
label confidently from the reads and treat those as truth (a *silver* standard),
accepting that the uncertain middle is set aside. We build it twice — once from
OUR pileup reads (caller-neutral) and once from HipSTR's per-read fields — and
take the loci where both agree as the **confident core**:

  * true100  — a non-reference tract length that is read-dominant (>= HOM_FRAC of
               a plant's reads, plant depth >= MIN_DEPTH) in >= RECUR independent
               plants, or a balanced heterozygote corroborated by another carrier.
  * false100 — no plant is read-dominated by any non-reference length (all stutter).

A caller is then scored on the core: recall = fraction of confident-true loci it
emits as a whole-unit length variant; false-positive rate = fraction of
confident-false loci it emits.

Used to compare emission models (heuristic gates / BIC / freebayes-marginal) on an
identical, caller-independent yardstick. Grounded in the reads, so partly
co-defines "segregation" — the non-circular check is an orthogonal truth set.

Usage:
  # core sizes
  uv run --no-project silver_standard.py
  # score a caller VCF on the confident core
  uv run --no-project silver_standard.py --score <ours.ssr.vcf>
  # override inputs (defaults point at results_rerun_20260708)
  uv run --no-project silver_standard.py --reads R.tsv --catalog C --hipstr H.vcf.gz --score V
"""
import argparse
import gzip
from collections import Counter, defaultdict

# --- tunables (kept in one place; the depth-shallow tomato panel needs these low) ---
MIN_DEPTH = 4       # reads before a plant's fractions are trusted
HOM_FRAC = 0.75     # modal-length fraction that counts as read-dominant
BAL_MIN = 0.40      # balanced-het minor fraction
RECUR = 2           # independent plants for recurrence
MIN_ASSESSABLE = 6  # assessable plants before a locus can be judged

_BASE = "benchmarks/ssr_tomato1"
DEF_READS = f"{_BASE}/results_rerun_20260708/our_reads.tsv"
DEF_CATALOG = f"{_BASE}/results/ours/ssr_tomato1.ssr.catalog"
DEF_HIPSTR = f"{_BASE}/results_ssr15k/hipstr/cohort.str.vcf.gz"


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith((".gz", ".bgz")) else open(p)


def _gt(gt):
    a = gt.replace("|", "/").split("/")
    return None if any(x in (".", "") for x in a) else tuple(int(x) for x in a)


def _info_int(info, key):
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            try:
                return int(kv.split("=")[1])
            except ValueError:
                return None
    return None


def load_catalog(path):
    """bgzip TSV (chrom start end motif ...) -> {(chrom, start+1): (reflen, period)}.
    Read with gzip (BGZF is gzip-framed) so no bgzip binary is needed."""
    cat = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                cat[(p[0], int(p[1]) + 1)] = (int(p[2]) - int(p[1]), len(p[3]))
    return cat


def load_reads(path):
    """our_reads.tsv (ssr_slip_dump) -> {(chrom,pos): {sample: {length: count}}}."""
    d = defaultdict(lambda: defaultdict(Counter))
    with open(path) as fh:
        next(fh)
        for line in fh:
            s, ch, st, seq, cnt = line.rstrip("\n").split("\t")
            d[(ch, int(st))][s][len(seq)] += int(cnt)
    return d


def parse_hipstr(path):
    """(chrom,pos) -> homalt/balhet/carrier/n_assess from HipSTR's own GB/PDP fields."""
    loci = {}
    for line in _open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        period = _info_int(f[7], "PERIOD")
        keys = f[8].split(":")
        gi = keys.index("GT")
        gbp = keys.index("GB") if "GB" in keys else None
        pdpp = keys.index("PDP") if "PDP" in keys else None
        homalt, balhet, carrier, n_assess = Counter(), set(), Counter(), 0
        for col in f[9:]:
            parts = col.split(":")
            idx = _gt(parts[gi])
            if idx is None:
                continue
            gb = parts[gbp] if gbp is not None and gbp < len(parts) else None
            if not gb or gb in (".", ""):
                continue
            try:
                diffs = [int(x) for x in gb.replace("|", "/").split("/")]
            except ValueError:
                continue
            n_assess += 1
            wu = [d for d in diffs if d != 0 and period and d % period == 0]
            for d in set(wu):
                carrier[d] += 1
            het = len(set(idx)) > 1
            if not het and wu:
                homalt[diffs[0]] += 1
            elif het and pdpp is not None:
                try:
                    a, b = (float(x) for x in parts[pdpp].split("|"))
                    mf = min(a, b) / (a + b) if a + b > 0 else 0
                except (ValueError, IndexError):
                    mf = 0
                if mf >= BAL_MIN:
                    for d in wu:
                        balhet.add(d)
        loci[(f[0], int(f[1]))] = (homalt, balhet, carrier, n_assess)
    return loci


def pileup_signals(persample, reflen, period):
    homalt, balhet, carrier, n_assess = Counter(), set(), Counter(), 0
    for lens in persample.values():
        total = sum(lens.values())
        if total < MIN_DEPTH:
            continue
        wu = {L: c for L, c in lens.items() if period and (L - reflen) % period == 0}
        if not wu:
            continue
        n_assess += 1
        top = sorted(wu.items(), key=lambda kv: -kv[1])
        lmod, cmod = top[0]
        dmod, fmod = lmod - reflen, cmod / total
        if dmod != 0 and fmod >= HOM_FRAC:
            homalt[dmod] += 1
            carrier[dmod] += 1
        elif len(top) >= 2:
            l2, c2 = top[1]
            if min(cmod, c2) / (cmod + c2) >= BAL_MIN:
                for d in (dmod, l2 - reflen):
                    if d != 0:
                        balhet.add(d)
                        carrier[d] += 1
    return homalt, balhet, carrier, n_assess


def classify(homalt, balhet, carrier, n_assess):
    if n_assess < MIN_ASSESSABLE:
        return "undetermined"
    if any(n >= RECUR for n in homalt.values()) or any(carrier.get(d, 0) >= 2 for d in balhet):
        return "true100"
    if not homalt and not balhet:
        return "false100"
    return "doubtful"


def build_core(reads_path=DEF_READS, catalog_path=DEF_CATALOG, hipstr_path=DEF_HIPSTR):
    """-> (core_true set, core_false set): loci both silver standards agree on."""
    reads = load_reads(reads_path)
    cat = load_catalog(catalog_path)
    hip = parse_hipstr(hipstr_path)
    pileup_cls = {}
    for (ch, st), ps in reads.items():
        m = cat.get((ch, st))
        if m:
            pileup_cls[(ch, st)] = classify(*pileup_signals(ps, m[0], m[1]))
    hip_cls = {k: classify(*v) for k, v in hip.items()}
    shared = set(pileup_cls) & set(hip_cls)
    core_true = {k for k in shared if pileup_cls[k] == "true100" == hip_cls[k]}
    core_false = {k for k in shared if pileup_cls[k] == "false100" == hip_cls[k]}
    return core_true, core_false


def emitted_whole_unit(path):
    """(chrom,pos) -> True if the caller PASS-emitted a whole-unit length variant."""
    out = {}
    for line in _open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[6] not in (".", "PASS", ""):
            out[(f[0], int(f[1]))] = False
            continue
        ref, alt = f[3], f[4]
        rl = len(ref)
        alleles = [ref] + ([] if alt == "." else alt.split(","))
        period = _info_int(f[7], "PERIOD")
        gi = f[8].split(":").index("GT")
        var = False
        for col in f[9:]:
            idx = _gt(col.split(":")[gi])
            if idx is None:
                continue
            if any((len(alleles[i]) - rl) != 0 and period and (len(alleles[i]) - rl) % period == 0
                   for i in idx):
                var = True
                break
        out[(f[0], int(f[1]))] = var
    return out


def score_vcf(vcf_path, core_true, core_false):
    called = emitted_whole_unit(vcf_path)
    tp = sum(1 for k in core_true if called.get(k, False))
    fp = sum(1 for k in core_false if called.get(k, False))
    return {
        "true": len(core_true), "false": len(core_false),
        "tp": tp, "fn": len(core_true) - tp, "fp": fp,
        "recall": tp / len(core_true) if core_true else 0.0,
        "precision": tp / (tp + fp) if (tp + fp) else 0.0,
        "fp_rate": fp / len(core_false) if core_false else 0.0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", default=DEF_READS)
    ap.add_argument("--catalog", default=DEF_CATALOG)
    ap.add_argument("--hipstr", default=DEF_HIPSTR)
    ap.add_argument("--score", metavar="OURS_VCF", help="score this VCF on the confident core")
    a = ap.parse_args()
    core_true, core_false = build_core(a.reads, a.catalog, a.hipstr)
    print(f"confident core: {len(core_true)} true100, {len(core_false)} false100")
    if a.score:
        s = score_vcf(a.score, core_true, core_false)
        print(f"\n{a.score}")
        print(f"  recall    {s['recall']:6.1%}  ({s['tp']}/{s['true']}; FN {s['fn']})")
        print(f"  precision {s['precision']:6.1%}")
        print(f"  FP rate   {s['fp_rate']:6.2%}  ({s['fp']}/{s['false']})")


if __name__ == "__main__":
    main()
