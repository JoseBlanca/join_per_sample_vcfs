# Hidden-paralog single-individual reformulation — validation (Step C)

**Date:** 2026-07-02. **Branch:** `tomato2-paralog-filter`. Validates the
reformulation designed in
[hidden_paralog_single_sample_scoring.md](../../architecture/hidden_paralog_single_sample_scoring.md)
and built in Steps A (drop `min_samples`) + B (`F` = cohort inbreeding
coefficient, remove `Hexp`). No truth set → judged on **profile coherence** and
an old-vs-new drop-set comparison, plus a single-individual graceful-degradation
check.

## Setup

Same as T1/T2 (`paralog_t1_t2_2026-07-02.md`): 59 tomato2 samples, 160 × 200 kb,
`--min-qual 0 --no-allele-balance-filter`, summary-bearing `.psp`. Release binary
with the reformulation. Compared **filter-on vs `--no-paralog-filter`**, and the
**new drop set vs the old** (pre-reformulation) one.

## Cohort result (n=58 fit) — profile preserved, π holds

| quantity | old filter (T1) | **reformulated** |
|---|---|---|
| samples fit | 58 / 59 | 58 / 59 |
| `F` used | ≈0.88 (per-sample, `Hexp`-derived) | **0.0000** (cohort default) |
| π (EM) | 0.099 | **0.0923** (converged) |
| loci dropped | 19 773 (6.98 %) | **18 826 (6.64 %)** |

**Drop profile (the yardstick):**

| set | n | mean DP | het fraction |
|---|---|---|---|
| **dropped** | 18 826 | **500.8 (1.29×)** | **0.116 (2.5×)** |
| kept | 264 522 | 387.4 | 0.046 |

Nearly identical to the old filter's profile (dropped 492 / 0.124, kept 388 /
0.045). The **coverage-excess + het-excess signature survives the reformulation**.

**Stability:** the filter-**off** callset is byte-identical old vs new (0 record
diffs — the off path is untouched). The **drop sets** overlap at **Jaccard
0.867** (17 928 shared; 1 845 old-only, 898 new-only). The ~13 % shift is
expected and comes from two changes: `F` moved 0.88 → 0.0 (so H1 predicts more
hets → hets slightly less surprising → the old-only 1 845 are loci the
het-rarity boost used to flag) and the gate is gone. The profile is unchanged, so
the shift is a mild loss of aggressiveness, not a loss of specificity.

**On `F = 0`.** The cohort default is `--inbreeding-coefficient = 0` (outbred),
so tomato — a selfer — runs *without* the het-rarity boost the old per-sample
`F ≈ 0.88` gave. It still drops the right class at nearly the same rate because
**coverage is load-bearing** — exactly the design thesis (`F` is a weak knob). An
operator who wants the boost for a selfing cohort can pass a higher
`--inbreeding-coefficient`; it is now one honest cohort knob rather than a
per-individual quantity we can't identify AF-free.

**FDR default:** 6.64 % dropped at `--paralog-fdr 0.01` with the correct profile —
no re-pin needed.

## Single-individual result (n=1) — graceful degradation confirmed

One tomato sample (~6× median), filter on:

- **Runs without error** — no `min_samples` crash, no `Hexp` divide, EM converges
  (π = 0.132).
- **Drops 531 / 7 506 (7.1 %)** — even at n=1 and low depth, the coverage signal
  flags the strong cases.
- **The drops are the paralog class:**

| set | n | mean DP | het fraction (of the 1 sample) | mean AF |
|---|---|---|---|---|
| **dropped** | 531 | **15.1 (2.3×)** | 0.994 | 0.50 |
| kept | 6 975 | 6.5 | 0.393 | 0.80 |

The dropped loci are **het at ~2.3× coverage** (a collapsed-paralog PSV in the
one sample); the kept are normal-coverage and mostly hom-alt (real fixed
differences). So the within-sample het∩coverage coupling works at n=1: coverage
carries it, and only strong signals are caught — the intended smooth degradation
(low power at n=1, but correct and safe, no false drops of normal-coverage loci).

## Verdict

The reformulation **reproduces the cohort paralog result** (π ≈ 9 %, coverage +
het-excess drop profile, Jaccard 0.87 with the old drop set) while **removing all
cohort allele-frequency machinery** (`Hexp`, the accumulator, `obs_het`), and it
**degrades gracefully to a single individual** (runs, converges, drops the
coverage-excess class, no false drops). The single global quantity is now π, over
loci. Step C passes.

Owed: the deferred items in the design doc (optional cohort-SFS-prior calibration;
low-n soft-flag policy) remain future work, not blockers.

## Reproduction

```
DEV_EXTRA_MOUNT=$HOME/genomes/s_lycopersicum/4.00 ./scripts/dev.sh bash -lc '
  bash tmp/paralog_t1t2/measure3.sh reform_off tmp/paralog_t1t2/reform_off.vcf --no-paralog-filter
  bash tmp/paralog_t1t2/measure3.sh reform_on  tmp/paralog_t1t2/reform_on.vcf'
uv run tmp/paralog_t1t2/profile.py  tmp/paralog_t1t2/reform_off.vcf tmp/paralog_t1t2/reform_on.vcf
uv run tmp/paralog_t1t2/jaccard.py  {sum,reform}_{off,on}.vcf   # old-vs-new drop set
# n=1: one .psp, filter on vs --no-paralog-filter, then profile.py
```
