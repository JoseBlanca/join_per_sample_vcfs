# Algorithm 3 vs 4 — the two-penalty delimiter comparison (Milestone D2)

**Date:** 2026-07-23
**Plan:** [alignment_best_path.md](../../ng/impl_plan/alignment_best_path.md) D2; spec
[alignment.md](../../ng/spec/alignment.md) §4.2, §10.3
**Harness:** [`examples/ssr_delimiter_comparison.rs`](../../../../examples/ssr_delimiter_comparison.rs)
(`cargo run --release --example ssr_delimiter_comparison`) — deterministic, seeded.

## The question

Both delimiters *measure* a read's repeat length. Algorithm 3 (`SsrFlatGapAligner`) charges one
flat per-base gap anywhere inside the tract; algorithm 4 (`SsrUnitSlipAligner`) prices a whole-unit
slip from the stutter model instead. Which measures better, and — the §4.2 worry — does the
two-penalty ruler **round interruptions away**, pulling the measurement toward a tidy in-frame
length? The plan (§10.3) requires the comparison be run on synthetic reads with known truth, split
by period, and at period 1 split by whether an indel is the repeat's own base or a foreign one.

## Method

A seeded simulator builds a read whose true tract length is known, from one of five scenarios, adds
1% per-base sequencing error, and runs both aligners. 20,000 reads per cell, reference tract 6
units, contraction-biased HipSTR-shaped parameters. Three numbers per cell:

- **accuracy** — fraction measured exactly right;
- **bias** — mean *signed* error `measured − true`. For a point-estimate *ruler* this is the analog
  of calibration: a non-zero bias is a systematic pull, which is exactly the §4.2 concern. (The
  marginal sense of calibration does not apply — a best-path delimiter emits a length, not a
  probability. §10.3's "measure calibration" is answered here as bias.)
- **spread** — mean *absolute* error.

Plus the **disagreement rate**: how often the two aligners measure different lengths.

## Results

Reproduced across three independent seeds (values within ~0.3 pp; the period-1 different-base cell
was 92.6% / 99.5% here, 92.4% / 99.5% and 92.5% / 99.4% on the other two). All figures in base
pairs. The flanks are chosen to have no *tiling collision* at either junction (see the harness) — an
earlier pick whose right flank started with period 3's own `motif[0]` let the boundary slide a whole
unit into the flank and collapsed period-3 clean accuracy to 83%, which is a flank artifact, not a
property of either aligner.

```
scenario                      per |  algo 3 (flat gap) | algo 4 (unit slip) |  disagree
                                  |    acc  bias  sprd |    acc  bias  sprd |      rate
------------------------------------------------------------------------------------------
clean                           1 |  99.2% +0.00  0.01 |  99.7% +0.00  0.00 |     0.56%
own-base indel (p1)             1 |  99.2% +0.00  0.01 |  99.6% +0.01  0.01 |     0.45%
different-base indel (p1)       1 |  92.6% -0.13  0.14 |  99.5% -0.00  0.01 |     7.03%
clean                           2 |  99.1% -0.01  0.02 |  99.8% +0.00  0.00 |     0.76%
substitution                    2 |  96.9% -0.07  0.08 |  99.6% +0.00  0.01 |     2.79%
out-of-frame indel              2 |  96.9% -0.06  0.07 |  98.5% -0.01  0.02 |     2.31%
clean                           3 |  99.1% -0.01  0.02 |  99.9% +0.00  0.00 |     0.85%
substitution                    3 |  97.2% -0.05  0.08 |  99.7% -0.00  0.01 |     2.57%
out-of-frame indel              3 |  97.3% -0.05  0.07 |  99.1% -0.00  0.01 |     2.21%
clean                           4 |  99.1% -0.01  0.03 |  99.9% +0.00  0.00 |     0.80%
substitution                    4 |  98.7% -0.06  0.07 |  99.8% -0.00  0.00 |     1.13%
out-of-frame indel              4 |  98.4% -0.05  0.06 |  99.3% -0.00  0.01 |     1.21%
```

Essentially no reads were unmeasured (both aligners anchor every read; at most one read of 20,000 in
a cell, and the same read for both), so every cell is a like-for-like measurement comparison.

## Findings

**1. Algorithm 4 is uniformly the more accurate and less biased ruler** — every period, every
scenario, higher accuracy and bias at or near zero.

**2. The interesting result is what algorithm 3 does, and it is the recorded production failure.**
Algorithm 3 carries a small but consistent **negative** bias: it systematically *under*-measures,
pulling the length toward the reference. This is the ~1.02 tract-transition inconsistency the port
documents (the un-normalised model penalises the length-changing path by ~1.02 per gapped base, a
pull toward the reference length; see `ssr_best_path_flat_gap.rs`) — surfacing here as measurement
bias. It is the same effect as the recorded production failure, "stutter reads pulled to reference."
Algorithm 4 prices slips by the stutter geometric rather than the flat gap, and does not carry it:
its bias is at or below 0.01 bp everywhere.

**3. The §4.2 fear does not materialise in measurement.** The worry was that the two-penalty ruler
would round interruptions toward tidy lengths. It does not: algorithm 4's bias is essentially zero,
including on the substitution and out-of-frame scenarios that carry an interruption. If either ruler
rounds toward the reference, it is **algorithm 3**, not 4.

**4. The spec's predicted period-1 cancellation does not occur — and the reason is instructive.**
The spec expected the own-base indel to favour algorithm 4 (direction asymmetry) and the
different-base indel to favour algorithm 3 (algorithm 4 "mis-routes" a foreign-base insertion as a
slip), so the two would cancel when averaged. Instead **algorithm 4 wins both**, and the
different-base cell is its *largest* win (92.6% → 99.5%). The mis-routing the spec worried about is a
*scoring* concern — algorithm 4 pays a base-error price for the foreign base rather than an
out-of-frame one — but for a *delimiter* the measured length is the same however the foreign base is
routed, so the mis-routing costs no accuracy.

**The mechanism, verified by a hand-read sweep** (D2 review): algorithm 3 is *not* uniformly wrong on
this case. It measures the +1 foreign insertion **correctly at every interior position** of the
homopolymer, and errs **only when the foreign base lands at the left tract boundary** (position 0),
where it collapses the measurement toward the reference length by a jitter-scaled 1–3 bp rather than
by a uniform one. That boundary-collapse model reproduces the aggregate −0.13 bias to two decimals; a
"reinterpret every insertion as a substitution" story does not, and an earlier draft of this report
told the wrong story. There is also a genuine definitional edge at position 0: the foreign base abuts
the flank, and "is it the tract's or the flank's?" has no crisp answer — the harness counts it as
tract, which is one defensible convention.

**This is why the split was mandatory (§10.3):** averaging own-base (a 0.7 pp gap) with different-base
(a 6.9 pp gap) would report a middling ~3.8 pp period-1 difference and hide that the whole effect
lives in the foreign-base case.

## Caveats, stated plainly

- **Algorithm 4 has no parity oracle** (arch §5). Its correctness rests on the differential against
  algorithm 3 (3666 clean cases, byte-identical) and its property tests — not on an external truth.
  This comparison inherits that: it measures *relative* behaviour, not absolute ground truth beyond
  the simulator's own.
- **The simulator builds motif-pure tracts**, which fits algorithm 4's motif-aware slip model well.
  Real repeats are less pure; the substitution and different-base scenarios inject impurity and
  algorithm 4 still wins there, but a genome of ragged, interrupted repeats could narrow the gap.
  This measures the model on clean-ish synthetic data, which is what §10.3 asks for and no more.
- **This is the *measurement* comparison, not the genotype comparison.** Whether the better ruler
  changes called genotypes spans this module and the genotyping that composes it (§5.1, §10.3), and
  that harness lives with the genotyping, not here.
- **Both rulers are already very accurate** (algorithm 3’s worst cell is 92.6%, most are >97%). The
  question this answers is not "does algorithm 3 work" — it does — but "is algorithm 4's stutter
  pricing a measurable improvement," and the answer is a consistent yes, concentrated exactly where
  the biology is (interruptions, and period-1 foreign-base insertions).

## Bottom line

On synthetic reads with known truth, **algorithm 4 is the better delimiter**: uniformly more
accurate, and free of the small reference-pull bias that algorithm 3 inherits from its flat-gap
model. The §4.2 hypothesis that a two-penalty ruler rounds interruptions away is **not** borne out —
the pull, where it exists, is algorithm 3's. The decision to *adopt* algorithm 4 as the STR
delimiter is not this module's to make alone (it turns on the genotype comparison, §10.3), but the
measurement evidence points one way.
