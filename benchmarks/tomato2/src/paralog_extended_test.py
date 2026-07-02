#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy"]
# ///
"""Extended paralog statistic (copies T can exceed 2) + synthetic tests.

H1 (real variant): genotype g~HWE(p,F); coverage~Normal(1,σ0) indep of g; VAF g/2.
H2-old (capped):   dosage d∈{0,1,2}; coverage 1+d/2; VAF d/(2+d)  (caps at 2×).
H2-ext (new):      each sample is non-carrier (cov 1, hom-ref) OR a carrier with total
                   copies T and m mutant copies -> coverage T/2, VAF = m/T. T and m are
                   per-locus (maximised), so high copy ⇒ low VAF for a single PSV, and
                   coverage above 2× is representable.

We print, for several synthetic loci, the LR under the old vs the new H2 and the
empirical-Bayes posterior (using π≈0.094 ⇒ offset −2.27), to see what is flagged.
"""
import numpy as np
from scipy.special import logsumexp, expit

SIG0 = 0.26
EPS = 0.01
F = 0.87                       # representative per-sample inbreeding (autogamous)
PRIOR_OFFSET = -2.27           # log(π/(1-π)) with π≈0.094
QGRID = np.linspace(0.004, 0.6, 50)
PGRID = np.linspace(0.004, 0.6, 50)
T_CARRIER = [3, 4, 5, 6, 8]    # carrier total copies (coverage 1.5,2,2.5,3,4)


def nlog(x, mu, sig):
    return -0.5 * ((x - mu) / sig) ** 2 - np.log(sig * np.sqrt(2 * np.pi))


def blog(k, n, v):             # binomial kernel (C(n,k) cancels in the LR)
    return k * np.log(v) + (n - k) * np.log(1 - v)


def hwe(f, F):                 # [hom-ref, het, hom-alt] inbreeding probs
    return np.array([(1 - f) ** 2 + F * f * (1 - f),
                     2 * f * (1 - f) * (1 - F),
                     f ** 2 + F * f * (1 - f)])


def logL_H1(c, k, n):
    best = -np.inf
    cov = np.sum(nlog(c, 1.0, SIG0))
    vaf = np.array([EPS, 0.5, 1 - EPS])
    for p in PGRID:
        lp = np.log(hwe(p, F))
        per = np.array([logsumexp(lp + blog(k[i], n[i], vaf)) for i in range(len(c))])
        best = max(best, per.sum() + cov)
    return best


def logL_H2_old(c, k, n):
    best = -np.inf
    mu = np.array([1.0, 1.5, 2.0]); sg = SIG0 * np.sqrt(mu)
    vaf = np.array([EPS, 1 / 3, 0.5])
    for q in QGRID:
        lp = np.log(hwe(q, F))           # d=0,1,2 with inbreeding
        per = []
        for i in range(len(c)):
            comp = lp + nlog(c[i], mu, sg) + blog(k[i], n[i], vaf)
            per.append(logsumexp(comp))
        best = max(best, np.sum(per))
    return best


def logL_H2_ext(c, k, n):
    best = -np.inf
    for T in T_CARRIER:
        muc, sgc = T / 2, SIG0 * np.sqrt(T / 2)
        for m in range(1, T):            # mutant copies (VAF = m/T)
            vc = m / T
            for q in QGRID:
                p0 = (1 - q) ** 2 + F * q * (1 - q)        # non-carrier
                pc = 1 - p0                                # carrier (het+hom dup)
                l0 = np.log(p0); lc = np.log(pc)
                tot = 0.0
                for i in range(len(c)):
                    non = l0 + nlog(c[i], 1.0, SIG0) + blog(k[i], n[i], EPS)
                    car = lc + nlog(c[i], muc, sgc) + blog(k[i], n[i], vc)
                    tot += logsumexp([non, car])
                best = max(best, tot)
    return best


def scenario(name, c, k, n):
    h1 = logL_H1(c, k, n)
    lr_old = logL_H2_old(c, k, n) - h1
    lr_new = logL_H2_ext(c, k, n) - h1
    print(f"{name:46s}  LR_old={lr_old:+8.1f} (post {expit(lr_old+PRIOR_OFFSET):.3f})"
          f"   LR_new={lr_new:+8.1f} (post {expit(lr_new+PRIOR_OFFSET):.3f})")


N = 59
def locus(carrier_specs, base_cov=1.0, base_n=12):
    """carrier_specs: list of (cov, vaf, n) for the carriers; rest = base hom-ref."""
    c = np.full(N, base_cov); k = np.zeros(N); n = np.full(N, base_n, float)
    for i, (cov, vaf, nn) in enumerate(carrier_specs):
        c[i] = cov; n[i] = nn; k[i] = round(vaf * nn)
    return c, k, n


print("scenario                                          old model (cap 2×)        extended model")
# Q2: ONE sample, 2× coverage, het (VAF 0.5); rest normal hom-ref
scenario("1 carrier @2× het (VAF.5), 58 normal hom-ref", *locus([(2.0, 0.5, 24)]))
# contrast: ONE sample, NORMAL coverage, het = real singleton variant -> should NOT flag
scenario("1 carrier @1× het (real singleton variant)   ", *locus([(1.0, 0.5, 12)]))
# Q1: high-copy, low-VAF PSV (the case the cap misses): several @4× with VAF 1/8
scenario("8 carriers @4× low-VAF (1/8), 51 hom-ref      ", *locus([(4.0, 0.125, 48)] * 8))
# divergent balanced 2× paralog at moderate freq (both models should catch)
scenario("12 carriers @2× balanced het, 47 hom-ref     ", *locus([(2.0, 0.5, 24)] * 12))
# a single 3× het carrier (between 2× and high): cap can't reach 3×
scenario("1 carrier @3× (VAF 1/6), 58 normal hom-ref   ", *locus([(3.0, 1 / 6, 36)]))
