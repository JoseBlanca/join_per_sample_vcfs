# SSR emission-gap drop attribution (ssr_tomato1 vs HipSTR)

**Date:** 2026-07-08 · **Branch:** `ssr-drop-fate` · **Bench:** `benchmarks/ssr_tomato1`
(51-sample *S. lycopersicum* selfer cohort, ours-on-`main` vs HipSTR v0.6.2)

## Question

On `ssr_tomato1`, our SSR caller genotypes concordantly with HipSTR (96.5% of
comparable cells) but **emits only 47% of HipSTR's length-polymorphic loci**.
Prior analysis localised the loss: of 1,578 HipSTR length-variable comparable
loci, **763 are "genuine drops"** — our caller emits no row at all. All 763 are
in our catalog, and only ~3% of the total loss is filtering. So the drops are
loci we admit (`Admission::Pass`) but call **length-monomorphic** across the
cohort → silently dropped (`driver.rs`, emit-iff-variable rule). *Which*
mechanism collapses each real variant to monomorphic was untraced.

## Instrumentation

Added a zero-cost-when-off diagnostic sidecar, `PVC_SSR_FATE_TSV=<path>`
([`src/ssr/cohort/driver.rs`](../../../src/ssr/cohort/driver.rs)). For every
catalog locus in the genotyping sweep it writes one row
`chrom  pos  reflen  n_alt  admit  fate`, where `fate` splits the silent
monomorphic-PASS drop into its three mechanisms:

- `drop_no_alt_candidate` — the candidate set held **no** non-reference allele
  (the slip/variant tract was never nominated; candidate assembly / §5.2 gates).
- `drop_em_hom_ref` — a non-ref candidate existed, but the **EM** called every
  sample hom-ref (the stutter model out-voted the slip).
- `drop_fp_control` — the EM produced a non-ref call, but **FP-control**
  (`min_allele_balance` 0.20 / `no_call_gq` 15) no-called the ALT carrier(s),
  leaving the site monomorphic.

`emitted_variable` and `filtered` cover the non-drop fates. The observation is
read-only: **the VCF is byte-identical with the flag on or off** (verified by
`diff` against the baseline cohort VCF). All 240 `ssr::cohort` tests pass.

## Result — fate of the 763 genuine drops

| Mechanism | count | share |
|---|--:|--:|
| **`drop_fp_control`** (FP-control no-called the ALT) | **374** | **49.0%** |
| **`drop_em_hom_ref`** (EM collapsed to hom-ref) | **316** | **41.4%** |
| `drop_no_alt_candidate` (never nominated) | 65 | 8.5% |
| filtered (boundary-shift join slop) | 8 | 1.0% |

**The headline: the emission gap is downstream genotyping conservatism, not
candidate nomination.** 90% of the loss is FP-control no-calling (49%) + EM
hom-ref collapse (41%). Candidate nomination — the previously-suspected
interruption/§5.2 path — is the *smallest* term at 8.5%.

Breakdowns:

- **`drop_em_hom_ref` is a dinucleotide story:** 226/316 (72%) are period-2.
  The stutter model treats the short-motif slip as noise and votes hom-ref.
- **`drop_fp_control`** spans periods more evenly; these are loci where we *did*
  call a variant and then killed it on allele balance / GQ — the knobs
  `min_allele_balance` and `no_call_gq` directly control this trade.
- Interruptions are a minor slice *within* each mechanism (20–32%), not a
  driver.

**Accuracy caveat (unchanged):** HipSTR is `F`-blind and tomato is a selfer
(F_IS≈0.82); short repeats are where HipSTR over-calls stutter as hets. The
HipSTR-singleton share (variant in exactly 1 of 51 samples) is **~21–23% and
uniform across all three mechanisms** — so an estimated fifth of these "drops"
are likely HipSTR false-polymorphisms, spread evenly, not concentrated in any
one mechanism. Without STR truth we cannot separate the rest.

## Levers (ranked by drops recoverable)

1. **FP-control (≈374 loci).** `min_allele_balance` (0.20) and `no_call_gq` (15)
   are dev-config knobs. Loosening them recovers real imbalanced hets — but also
   readmits HipSTR-style false hets. Sweep against a recovered-vs-new-FP curve.
2. **EM hom-ref at dinucleotides (≈316 loci).** The stutter level / genotype
   prior at period-2 is too aggressive toward hom-ref. Read-model / prior work,
   dinucleotide-specific. Higher risk of manufacturing FPs.
3. **Candidate nomination (≈65 loci).** Smallest lever; the §5.2 promotion gates
   (`min_same_length_*`) are largely *not* the bottleneck here.

## Reproduce

```
# Stage 2 with the sidecar on (reuses existing .ssr.psp; ~9 s):
PVC_SSR_FATE_TSV=/path/fate.tsv \
  pop_var_caller ssr-call <psp>/*.ssr.psp --catalog <cat> --output out.vcf

# Attribute the drops:
uv run --no-project benchmarks/ssr_tomato1/scripts/attribute_drops.py \
  --ours out.vcf --hipstr <hipstr>/cohort.str.vcf.gz --fate /path/fate.tsv
```
