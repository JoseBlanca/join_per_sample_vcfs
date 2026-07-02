#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy"]
# ///
"""Simulation with planted truth, to tune a mqdiff term in the paralog LR.

Classes (each = 59-sample locus; we generate per-sample coverage c, alt/total
reads k/n, and a per-locus mqdiff):
  real_variant     single-copy, normal MAPQ            -> KEEP
  introgression    single-copy, DIVERGENT allele (mqdiff<0), can be hom-alt -> KEEP
  clean_paralog    full 2x excess, balanced het, mild mqdiff   -> REMOVE (easy)
  diverging_paralog MILD excess (reads partly MAPQ-filtered) + DIVERGENT (mqdiff<0) -> REMOVE (the FN gap)
  single_diverging  same but only ONE carrier (private dup)    -> REMOVE (hard)

We run the marginal LR (production model) with a mqdiff term of weight w
(contribution = w * (-mqdiff), only negative mqdiff helps) and report the flag
rate per class vs w, to read the diverging-paralog-recall / introgression-loss
tradeoff off known truth.
"""
import numpy as np
from scipy.special import logsumexp, expit

rng = np.random.default_rng(0)
N = 59
DEPTH = 6
SIG0 = 0.26
EPS = 0.01
F = 0.87
LRCUT = 2.08                      # production 50%-posterior cut
PGRID = np.linspace(0.004, 0.6, 50)
QEXT = np.linspace(0.004, 0.6, 40)
T_CARRIER = np.array([3, 4, 6, 8])
T_MEAN = T_CARRIER / 2.0
T_SIG = SIG0 * np.sqrt(T_MEAN)
_TM = [(ti, m) for ti, T in enumerate(T_CARRIER) for m in sorted({1, int(T // 2)}) if 1 <= m <= T - 2]
TM_TIDX = np.array([ti for ti, _ in _TM])
TM_VC = np.array([m / T_CARRIER[ti] for ti, m in _TM])
LOG_VC, LOG1M_VC = np.log(TM_VC), np.log(1 - TM_VC)


def nlog(x, mu, sig):
    return -0.5 * ((x - mu) / sig) ** 2 - np.log(sig * np.sqrt(2 * np.pi))


def hwe(f):
    return np.array([(1 - f) ** 2 + F * f * (1 - f), 2 * f * (1 - f) * (1 - F), f ** 2 + F * f * (1 - f)])


_VAF1 = np.array([EPS, 0.5, 1 - EPS])
_LOGV1, _LOG1V1 = np.log(_VAF1), np.log(1 - _VAF1)
_LP_GRID = np.log(np.array([hwe(p) for p in PGRID]))        # (nP,3)


def logL_H1(c, k, n):                         # marginal over p (vectorised)
    cov = np.sum(nlog(c, 1.0, SIG0))
    kern = k[:, None] * _LOGV1[None, :] + (n - k)[:, None] * _LOG1V1[None, :]   # (N,3)
    term = _LP_GRID[:, None, :] + kern[None, :, :]          # (nP,N,3)
    summed = logsumexp(term, axis=2).sum(axis=1)            # (nP,)
    return logsumexp(summed) - np.log(len(PGRID)) + cov


def logL_H2(c, k, n):                         # marginal over q x (T,m)
    nbase = nlog(c, 1.0, SIG0) + (k * np.log(EPS) + (n - k) * np.log(1 - EPS))
    cov_car = nlog(c[:, None], T_MEAN[None, :], T_SIG[None, :])
    allele = k[:, None] * LOG_VC[None, :] + (n - k)[:, None] * LOG1M_VC[None, :]
    cb = cov_car[:, TM_TIDX] + allele                       # (N, nTM)
    P0 = (1 - QEXT) ** 2 + F * QEXT * (1 - QEXT)
    lp0, lpc = np.log(P0), np.log(1 - P0)                   # (nq,)
    A = lp0[None, None, :] + nbase[:, None, None]           # (N,1,nq)
    B = lpc[None, None, :] + cb[:, :, None]                 # (N,nTM,nq)
    g = np.logaddexp(A, B).sum(axis=0)                      # (nTM,nq)
    return logsumexp(g.ravel()) - np.log(g.size)


def base_lr(c, k, n):
    return logL_H2(c, k, n) - logL_H1(c, k, n)


# ---- generators ----
def _reads(c):
    n = np.maximum(1, np.round(c * DEPTH)).astype(int)
    return n


def gen_real(divergent):
    p = rng.uniform(0.05, 0.4)
    g = rng.choice(3, size=N, p=hwe(p))
    c = np.clip(rng.normal(1.0, SIG0, N), 0.2, None)
    n = _reads(c)
    vaf = np.array([EPS, 0.5, 1 - EPS])[g]
    k = rng.binomial(n, vaf)
    mq = rng.normal(-15, 4) if divergent else rng.normal(0, 3)
    return c, k.astype(float), n.astype(float), float(mq)


def gen_paralog(diverging, single=False):
    q = (1.0 / N) if single else rng.uniform(0.1, 0.4)
    d = rng.choice(3, size=N, p=hwe(q))             # 0 non-carrier, 1/2 carrier
    carrier = d >= 1
    c = np.clip(rng.normal(1.0, SIG0, N), 0.2, None)
    if diverging:
        cov_car = rng.uniform(1.25, 1.6)            # MILD excess (reads partly filtered)
        vaf_car = rng.uniform(0.22, 0.35)           # partial alt
        mq = rng.normal(-15, 4)
    else:                                           # clean collapsed paralog
        cov_car = 2.0
        vaf_car = 0.5
        mq = rng.normal(-3, 4)
    c[carrier] = np.clip(rng.normal(cov_car, SIG0 * np.sqrt(cov_car), carrier.sum()), 0.2, None)
    n = _reads(c)
    k = rng.binomial(n, np.where(carrier, vaf_car, EPS))
    return c, k.astype(float), n.astype(float), float(mq)


CLASSES = {
    "real_variant   (KEEP)": lambda: gen_real(False),
    "introgression  (KEEP)": lambda: gen_real(True),
    "clean_paralog  (REMOVE)": lambda: gen_paralog(False),
    "diverging_par  (REMOVE)": lambda: gen_paralog(True),
    "single_diverg  (REMOVE)": lambda: gen_paralog(True, single=True),
}

M = 300
print(f"{M} loci/class, flag = LR>{LRCUT}.  flag-rate per class vs mqdiff weight w:\n")
sims = {name: [fn() for _ in range(M)] for name, fn in CLASSES.items()}
# cache base LR + mqdiff per locus once
cache = {name: [(base_lr(c, k, n), mq) for (c, k, n, mq) in loci] for name, loci in sims.items()}
header = "w".ljust(6) + "".join(name.split()[0].ljust(16) for name in CLASSES)
print(header)
for w in (0.0, 0.1, 0.2, 0.3, 0.5):
    row = f"{w:<6.1f}"
    for name, loci in cache.items():
        flags = sum((b + w * (-min(mq, 0.0))) > LRCUT for (b, mq) in loci)
        row += f"{100*flags/M:>5.0f}%".ljust(16)
    print(row)
print("\nwant: KEEP rows LOW, REMOVE rows HIGH. mqdiff (w>0) should lift diverging_par/single_diverg")
print("WITHOUT flagging introgression too much -> pick w at that knee.")
