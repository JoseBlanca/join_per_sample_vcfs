# /// script
# requires-python = ">=3.10"
# ///
"""Ours-vs-HipSTR STR genotype concordance (no truth set).

Both callers genotype the *same* locus set (HipSTR's --regions BED is
projected from our ssr-catalog). With no STR truth for tomato, this scores
*agreement*, not accuracy.

The two tools define an STR allele differently, and the report makes that
explicit instead of penalising it:
  * ours   : repeat LENGTH / copy number (REPCN). Same-length sequences are
             the same allele; length-monomorphic loci are not emitted.
  * HipSTR : full sequence HAPLOTYPE. A SNP inside the repeat is a distinct
             allele, so HipSTR emits "variants" that carry no length change.

So every HipSTR variable locus is classified by its ALT-vs-REF lengths
(period from the INFO PERIOD field):
  * length-poly : >=1 ALT whose length differs from REF by a whole motif
                  unit  -> the COMPARABLE set; our caller should also call it.
  * non-unit    : a length change that is not a whole motif unit (e.g. a 1bp
                  indel in a 5bp repeat) -> our period model flags notPeriodic.
  * seq-only    : every ALT is REF-length (a SNP in the repeat) -> outside
                  our REPCN model; our caller is correctly silent (not a miss).

Concordance unit: per-allele base-pair difference from REF, as a sorted
multiset (== HipSTR's GB). SNP-only HipSTR genotypes collapse to (0,..),
matching our length view. Loci join on (CHROM, POS); samples by name.

Usage:
  ssr_concordance.py --ours ours.vcf[.gz] --hipstr hipstr.vcf.gz [--out report.txt]
"""

import argparse
import gzip
import sys
from collections import Counter, defaultdict


def opener(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def parse_vcf(path):
    """-> (samples, {(chrom,pos): Locus})  with Locus = (ref, [alts], period, {sample: gt|None})."""
    samples = []
    loci = {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, info, fmt = f[0], int(f[1]), f[3], f[4], f[7], f[8]
            alts = [] if alt == "." else alt.split(",")
            period = _info_int(info, "PERIOD")
            gt_i = fmt.split(":").index("GT")
            calls = {name: _gt_indices(col.split(":")[gt_i]) for name, col in zip(samples, f[9:])}
            loci[(chrom, pos)] = (ref, alts, period, calls)
    return samples, loci


def _info_int(info, key):
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            try:
                return int(kv.split("=")[1])
            except ValueError:
                return None
    return None


def _gt_indices(gt):
    parts = gt.replace("|", "/").split("/")
    if any(p in (".", "") for p in parts):
        return None
    return [int(p) for p in parts]


def rel_genotype(ref, alts, idx):
    alleles = [ref] + alts
    return tuple(sorted(len(alleles[i]) - len(ref) for i in idx))


def classify_hipstr_locus(ref, alts, period):
    """'mono' | 'length-poly' | 'non-unit' | 'seq-only' (period unknown -> 'unknown')."""
    if not alts:
        return "mono"
    if not period:
        return "unknown"
    deltas = [len(a) - len(ref) for a in alts]
    if any(d != 0 and d % period == 0 for d in deltas):
        return "length-poly"
    if any(d != 0 for d in deltas):  # length change, but not whole-unit
        return "non-unit"
    return "seq-only"  # all ALTs == REF length


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ours", required=True)
    ap.add_argument("--hipstr", required=True)
    ap.add_argument("--out")
    args = ap.parse_args()

    o_samples, ours = parse_vcf(args.ours)
    h_samples, hip = parse_vcf(args.hipstr)

    shared_samples = [s for s in o_samples if s in set(h_samples)]
    name_map = {s: s for s in shared_samples}
    if not shared_samples and len(o_samples) == 1 and len(h_samples) == 1:
        name_map = {o_samples[0]: h_samples[0]}
        shared_samples = [o_samples[0]]

    # Classify every HipSTR locus.
    hip_class = {key: classify_hipstr_locus(ref, alts, period) for key, (ref, alts, period, _) in hip.items()}
    class_counts = Counter(hip_class.values())

    in_ours = set(ours)
    # The comparable set: HipSTR loci with a whole-unit length polymorphism.
    comparable = [k for k, c in hip_class.items() if c == "length-poly"]
    comp_emitted = [k for k in comparable if k in in_ours]
    comp_absent = [k for k in comparable if k not in in_ours]

    # Concordance on the comparable set, per shared sample-cell called by both.
    both_called = concordant = 0
    diff_hist = defaultdict(int)
    for key in comp_emitted:
        o_ref, o_alts, _, o_calls = ours[key]
        h_ref, h_alts, _, h_calls = hip[key]
        for osamp in shared_samples:
            oi, hi = o_calls.get(osamp), h_calls.get(name_map[osamp])
            if oi is None or hi is None:
                continue
            both_called += 1
            og, hg = rel_genotype(o_ref, o_alts, oi), rel_genotype(h_ref, h_alts, hi)
            if og == hg:
                concordant += 1
            else:
                diff_hist[sum(og) - sum(hg)] += 1

    def pct(n, d):
        return f"{100 * n / d:.1f}%" if d else "n/a"

    L = []
    L.append("=== Ours vs HipSTR — STR concordance (agreement, no truth) ===")
    L.append(f"ours VCF   : {args.ours}  ({len(o_samples)} samples, {len(ours)} loci emitted)")
    L.append(f"hipstr VCF : {args.hipstr}  ({len(h_samples)} samples, {len(hip)} loci emitted)")
    L.append(f"paired samples : {len(shared_samples)}")
    L.append("")
    L.append("HipSTR locus classes (by ALT-vs-REF length; the callers define alleles differently):")
    L.append(f"  mono (ref-only)            : {class_counts.get('mono', 0)}")
    L.append(f"  length-poly (whole-unit)   : {class_counts.get('length-poly', 0)}   <- COMPARABLE to our REPCN model")
    L.append(f"  non-unit indel             : {class_counts.get('non-unit', 0)}   (our model flags notPeriodic)")
    L.append(f"  seq-only (SNP in repeat)   : {class_counts.get('seq-only', 0)}   (outside our model; we are silent, not wrong)")
    if class_counts.get("unknown"):
        L.append(f"  unknown (no PERIOD)        : {class_counts['unknown']}")
    L.append("")
    L.append("COMPARABLE set — HipSTR whole-unit length-polymorphic loci:")
    L.append(f"  total                      : {len(comparable)}")
    L.append(f"  also emitted by ours       : {len(comp_emitted)}  ({pct(len(comp_emitted), len(comparable))})")
    L.append(f"  absent from ours (we call length-mono / no-call / filtered) : {len(comp_absent)}")
    L.append("")
    L.append("CONCORDANCE on comparable loci (sample cells called by BOTH, allele-length genotype):")
    L.append(f"  both-called cells          : {both_called}")
    L.append(f"  concordant                 : {concordant} / {both_called}  ({pct(concordant, both_called)})")
    if diff_hist:
        L.append("  discordant Σallele-bp-diff (ours − hipstr):")
        for d in sorted(diff_hist):
            L.append(f"    {d:+d} bp : {diff_hist[d]}")
    if comp_absent:
        L.append("")
        L.append("  comparable loci HipSTR called but ours did not emit (first 20):")
        for k in comp_absent[:20]:
            L.append(f"    {k[0]}:{k[1]}")

    report = "\n".join(L)
    print(report)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(report + "\n")


if __name__ == "__main__":
    main()
