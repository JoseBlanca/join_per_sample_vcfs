# /// script
# requires-python = ">=3.10"
# ///
"""Offline BIC model-selection prototype for the SSR emission decision.

Per locus, is the cohort's reads better explained by MONOMORPHIC (one length +
stutter) or POLYMORPHIC (an alt length segregates + stutter)? Scored on the silver
confident core (see `silver_standard.py`) to check whether the model-selection axis
separates real loci from systematic stutter — the structural question behind
replacing the heuristic emission gates.

This uses a deliberately SIMPLE per-locus stutter kernel (symmetric geometric in
whole units, rate fit per locus per model). It is a structure check, not the real
thing: the in-caller version uses the caller's per-chemistry-anchored read model,
which explains stutter far better and should sharpen the separation. Result on this
core: real loci median LRT ~+95 (strongly polymorphic), stutter loci median ~-2
(the extra allele earns nothing) — the axis is sound; the toy kernel caps precision.

  uv run --no-project bic_prototype.py
"""
import math
from collections import Counter

from silver_standard import (DEF_CATALOG, DEF_READS, build_core, load_catalog,
                             load_reads)

F_INBREED = 0.82
R_GEOM = 0.5
S_GRID = [0.02, 0.05, 0.08, 0.12, 0.16, 0.20, 0.25, 0.30, 0.38, 0.48, 0.60]
F_GRID = [0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]


def kern_log(d, s):
    if d == 0:
        return math.log(max(1 - s, 1e-9))
    return math.log(max((s / 2) * (1 - R_GEOM) * R_GEOM ** (abs(d) - 1), 1e-12))


def logsumexp(a):
    m = max(a)
    return m + math.log(sum(math.exp(x - m) for x in a))


def locus_dbic(k, reads, cat):
    """ΔBIC (BIC1 - BIC0). < 0 → polymorphic favoured. None if no in-frame alt."""
    meta = cat.get(k)
    if not meta:
        return None
    reflen, period = meta
    pooled, persamp = Counter(), []
    for lens in reads.get(k, {}).values():
        h, tot = Counter(), 0
        for L, c in lens.items():
            if (L - reflen) % period == 0:
                h[(L - reflen) // period] += c
                tot += c
        if tot:
            persamp.append(h)
            pooled.update(h)
    if not persamp:
        return None
    modal = max(pooled, key=lambda u: pooled[u])
    alt = [(u, c) for u, c in pooled.items() if u != modal]
    if not alt:
        return None
    v = max(alt, key=lambda uc: uc[1])[0] - modal
    samples = [{u - modal: c for u, c in h.items()} for h in persamp]
    n_reads = sum(sum(h.values()) for h in samples)

    def hom(hh, a, s):
        return sum(c * kern_log(u - a, s) for u, c in hh.items())

    lnL0 = max(sum(hom(hh, 0, s) for hh in samples) for s in S_GRID)
    best1 = -1e18
    for s in S_GRID:
        g00 = [hom(hh, 0, s) for hh in samples]
        gvv = [hom(hh, v, s) for hh in samples]
        g0v = [sum(c * math.log(0.5 * math.exp(kern_log(u, s)) + 0.5 * math.exp(kern_log(u - v, s)))
                   for u, c in hh.items()) for hh in samples]
        for f in F_GRID:
            lp00 = math.log(F_INBREED * (1 - f) + (1 - F_INBREED) * (1 - f) ** 2)
            lpvv = math.log(F_INBREED * f + (1 - F_INBREED) * f * f)
            lp0v = math.log(max((1 - F_INBREED) * 2 * f * (1 - f), 1e-12))
            tot = sum(logsumexp([lp00 + g00[i], lp0v + g0v[i], lpvv + gvv[i]])
                      for i in range(len(samples)))
            best1 = max(best1, tot)
    return -2 * (best1 - lnL0) + math.log(max(n_reads, 2))  # +1 param under M1


def main():
    core_true, core_false = build_core()
    reads = load_reads(DEF_READS)
    cat = load_catalog(DEF_CATALOG)
    dt = [d for k in core_true if (d := locus_dbic(k, reads, cat)) is not None]
    df = [d for k in core_false if (d := locus_dbic(k, reads, cat)) is not None]
    print(f"confident core scored: {len(dt)} true, {len(df)} false "
          f"(false skipped-as-monomorphic: {len(core_false) - len(df)})\n")
    print("call polymorphic when ΔBIC < -margin (FP over the FULL false core):")
    for margin in (0, 5, 10, 20, 40):
        tp = sum(1 for d in dt if d < -margin)
        fp = sum(1 for d in df if d < -margin)
        print(f"  margin {margin:3d}: recall {tp}/{len(core_true)} ({tp/len(core_true):.0%})   "
              f"FP {fp}/{len(core_false)} ({fp/len(core_false):.2%})")


if __name__ == "__main__":
    main()
