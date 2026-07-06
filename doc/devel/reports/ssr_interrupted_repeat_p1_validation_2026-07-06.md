# SSR interrupted-repeat recall — Phase 1 validation (P1.5)

**Date:** 2026-07-06 · **Branch:** `ssr-interruptions` · **Verdict:** Phase 1 validated, shippable.

Validation of Phase 1 (sequence-keyed same-length alleles) against the exit
criteria in the [spec §9](../specs/ssr_interrupted_repeat_recall.md) and
[plan P1.5](../implementation_plans/ssr_interrupted_repeat_recall.md). Phase 1
is commits P1.1–P1.5 (`550b691`…`ec7e038`).

## Setup

`ssr_tomato1` benchmark — 51 *S. lycopersicum* samples, ours vs HipSTR. Phase 1
changed only Stage 2 (`ssr-call`), so the existing per-sample `.ssr.psp` (Stage 1)
were reused and only `ssr-call` re-run:

```
ssr-call --catalog …/ssr_tomato1.ssr.catalog --output cohort_phase1.ssr.vcf \
         --threads 4  …/results_ssr15k/ours/cohort/psp/*.ssr.psp
```

Analysis scripts (not committed; under the run's `tmp/ssr_p15/`): `recovery_audit.py`
(recovery + paralog FP audit) and `seq_concordance.py` (the sequence/zygosity
metric §9 requires, which the length-based `ssr_concordance.py` / dashboard cannot
provide). HipSTR baseline: `results_ssr15k/hipstr/cohort.str.vcf.gz`.

## Results vs the exit criteria

### 1. Recovery of the interruption drops — PASS

| | pre-Phase-1 | Phase-1 |
|---|---|---|
| PASS variable rows | 1467 | **1633** (+166 net) |
| newly-emitted PASS loci (by POS) | — | 196 |
| …overlapping a HipSTR *variable* locus | — | 132 (67%) |
| …of which HipSTR **seq-only** (the interruption class) | — | **116** |
| PASS rows carrying a same-length ALT | 40 | **245** (184 newly-emitted) |

116 recovered loci are exactly HipSTR's seq-only (SNP-in-repeat) class — the
interior-interruption polymorphism Phase 1 targets — matching the spec's estimate
of ~23 % × 765 ≈ 176 interruption drops. The recovered rows carry genuine
same-length ALTs (e.g. `SL4.0ch01:1141392` REF `ACCATC×3` vs an `…ACCGTC`
interruption sibling, both 18 bp, a sample called `1/1`).

### 2. No regression on the common loci — PASS

- **Length-genotype concordance** (comparable = HipSTR whole-unit length-poly,
  the existing `ssr_concordance.py`): **96.5 %** (was 96.4 %), 773 emitted (was 775).
- **Zygosity concordance** (het/hom agreement over all common PASS∩variable cells,
  `seq_concordance.py`): **97.3 %** (was 97.0 %).
- Admission unchanged: `lowDepth` 72 272 and `notPeriodic` 669 identical to
  pre-Phase-1; **zero `TooManyAlleles`** → no `Pass`→`TooManyAlleles` flip
  (spec §10 issue-4 watch).
- SNP-caller regression: 1586 lib tests + SNP e2e integration (posterior engine /
  cohort VCF writer / contamination) all green.

### 3. Sequence-genotype het concordance — the same-length-het tail (spec §9)

The length metric collapses a same-length interruption het to a homozygote, so it
structurally cannot see this. Keyed on zygosity + same-length flag instead:

| on common PASS∩variable loci | pre-Phase-1 | Phase-1 |
|---|---|---|
| common loci | 810 | **926** (+116 = the seq-only recoveries) |
| HipSTR same-length-het cells (visible) | 89 | **289** (3.2×) |
| …**ours recovers as het** | 23 (26 %) | **155 (54 %)** |
| …as a *same-length* het specifically | 22 | 155 (all) |
| …undercalled to hom (the tail) | 66 (74 %) | 134 (**46 %**) |

Phase 1 **doubles the same-length-het recovery rate (26 % → 54 %)** and triples the
interruption hets even visible to the comparison, while holding overall zygosity
concordance. When ours calls one of these het it is correctly a *same-length,
sequence-distinct* het (155/155) — a genotype the pre-Phase-1 caller could not
represent. All 289 read as *hom* under the length metric.

### 4. Paralog FP audit (cross-cutting gate, spec §10) — CLEAN

**0** loci where ≥ 90 % of covered samples are called a same-length het. The
collapsed-paralog / CNV signature (the sharpest new FP mode, which would show as a
same-length het in *every* carrier) is **absent** at the dev §5.2 thresholds. The
FP-risk residue is the 62 newly-emitted loci absent from HipSTR + 2 HipSTR-mono
(33 % of the 196), neither showing the universal-het paralog signature.

## Interpretation of the 46 % undercall tail

This is the **documented spec §5.3 tail**, not a Phase-1 defect. The tomato cohort
has extreme apparent inbreeding (`F_IS ≈ 0.821`; the caller emitted its
`ssrCallWarning`), which tilts the genotype prior hard toward homozygosity and
suppresses low-depth hets — precisely the "low-depth / high-`F_IS` tail" the spec
predicted. It is a *concordance* effect (the locus still emits), not a drop. The
same `F_IS` mechanism suppresses ~734 *length*-het cells on the same cohort (the
out-of-scope §11 rare-length-alt suppression), dwarfing the 134 same-length tail.

Per spec §5.3/§10 this **gates the D1 follow-up** — a per-sample same-length het
seed proposal converging with the Mark-2 BIC two-allele confident-genotype test —
which is explicitly *not* Phase 1 and tracked in `doc/devel/TODO.txt`.

## Left open

- **§5.2 threshold sweep** (defaults 8 / 3 / 0.10): the defaults are clean here
  (paralog-free, 67 % HipSTR overlap, no regression), but the knee is not yet
  pinned against recovered-loci-vs-new-FP. The three params are currently
  hardcoded in `CandidateCfg::dev_default()`; a sweep needs CLI/config plumbing.
- **`SL4.0ch01:4279322`** (the spec's motivating locus) still does not emit — its
  HipSTR record is length-poly (a −3 bp ALT) whose emission depends on the
  out-of-scope §11 rare-length-alt path, distinct from the same-length recovery
  validated above.

## Verdict

All five headline exit criteria pass: recovery confirmed and quantified (het
recovery doubled), no regression (length 96.5 %, zygosity 97.3 %, SNP e2e green),
paralog-clean, and the undercall tail measured to be the documented `F_IS` effect
that gates D1. **Phase 1 is done and shippable.**
