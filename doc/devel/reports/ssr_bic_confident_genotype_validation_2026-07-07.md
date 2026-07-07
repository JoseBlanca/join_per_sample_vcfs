# SSR pre-pass D1 (BIC confident-genotype gate) — validation on ssr_tomato1

*Date: 2026-07-07. Branch `ssr-interruptions` (D1a `82de9ae`, D1b `d0ad7a8`, D1c `35dd6d4`).
Validates the model-based BIC confident-genotype gate + sequence-aware attribution
against the exit criteria of
[`implementation_plans/ssr_bic_confident_genotype.md`](../implementation_plans/ssr_bic_confident_genotype.md)
(step D1d) and the honest success metric of
[`specs/ssr_bic_confident_genotype.md`](../specs/ssr_bic_confident_genotype.md) §6.*

## Setup

Stage 1 is unaffected by D1, so only **`ssr-call`** was re-run on the existing 51-sample
`.ssr.psp` cohort (`benchmarks/ssr_tomato1/results_ssr15k/ours/cohort/psp/`, catalog
`results/ours/ssr_tomato1.ssr.catalog`), in-container. Two measurements:

1. **Frozen chemistry** — the ignored diagnostic `driver::tests::d1_dump_frozen_chemistry`
   runs the pre-pass over the 20 000-locus burn-in subset and prints `ε`, the base
   match/mismatch counts, and the confident-genotype composition. The **baseline** is the
   same diagnostic run at the pre-D1 commit `cc1f401` (Phase-1 tip, the length-histogram
   gate + length-only attribution), in a throwaway worktree.
2. **VCF concordance** — the full `ssr-call` → VCF (`tmp/ssr_d1/cohort_d1.ssr.vcf`),
   compared to the Phase-1 baseline VCF (`tmp/ssr_p15/cohort_postfix.ssr.vcf`) against
   the HipSTR cohort (`results_ssr15k/hipstr/cohort.str.vcf.gz`) with the P1.5 scripts
   (`benchmarks/lib/ssr_concordance.py` length metric; `tmp/ssr_p15/seq_concordance.py`
   sequence/zygosity metric).

## Results vs the exit criteria

### 1. ε de-inflation (the direct, committed win) — PASS

| pre-pass metric (20 000 burn-in loci) | Phase-1 baseline | D1 |
|---|---:|---:|
| frozen **`ε`** | 0.001452 | **0.001115** (**−23.2 %**) |
| within-tract `base_mismatch` | 1202 | **924** (−23 %) |
| `base_match` | 826 585 | 827 844 |
| samples with stats | 51 | 51 |

`ε` falls by **23 %** on the tomato cohort. The drop is exactly the same-length-het
minority reads that the old length-only attribution scored as substitutions against the
pure allele (the `BIAS NOTE`) and now attribute faithfully to their own composition.

### 2. Same-length hets enter the confident set (the enabling win) — PASS

| confident-genotype composition | Phase-1 baseline | D1 |
|---|---:|---:|
| homozygotes | 4621 | 4580 |
| length-separated hets | 1 | 1 |
| **same-length hets** | **0** | **49** |

The length gate resolved **zero** same-length hets (it mislabels them homozygotes); the
BIC gate resolves **49**, and their reads now de-contaminate `ε`/`θ` instead of poisoning
them. (The ~41 fewer homozygotes ≈ the same-length hets no longer collapsed to homs.)

### 3. No regression on the common loci — PASS

| VCF metric (ours vs HipSTR) | Phase-1 baseline | D1 |
|---|---:|---:|
| **length-genotype concordance** (`ssr_concordance.py`) | 96.5 % (27886/28903) | **96.5 %** (27899/28914) |
| comparable loci emitted | 773 | **774** |
| **zygosity concordance** (`seq_concordance.py`) | 97.3 % (35669/36675) | **97.3 %** (35678/36677) |
| FP hets (ours het, HipSTR hom) | 138 | **137** |
| total records emitted | 74 574 | 74 573 |

Both concordance metrics are **flat**, FP hets tick down by one, and D1 emits one
*more* comparable locus. Of the ~74 k records, 3131 (~4.2 %) changed content — the
expected genotype/GQ tweaks from the slightly-shifted frozen params (ε 0.00145→0.00112)
feeding the genotyping EM; none move the concordance.

### 4. SNP caller end-to-end tests — PASS

All SNP integration suites green: `cohort_cli_integration` (12), `posterior_engine_integration`
(11), `pileup_cli_integration` (11), `contamination_estimation_integration` (5). D1
touches only `src/ssr/cohort/`, so the SNP path is untouched by construction.

*(One pre-existing `proptest` unit test, `var_calling::posterior_engine::tests::
larger_ref_pseudocount_cannot_increase_p_alt`, fails on a saved regression seed by a
1e-9 margin. `posterior_engine.rs` is byte-identical to committed HEAD — the failure is
in unmodified SNP code, independent of D1, and is not an end-to-end test. Flagged, not
owned by this change.)*

## Interpretation — what D1 delivers, and what it deliberately does not

**D1's committed target is met and measured:** `ε` de-inflates 23 %, and 49 same-length
hets enter the confident chemistry seed, with no regression on either concordance metric.

**The downstream het-recovery tail is unchanged — as predicted (spec §6).**

| same-length-het recovery (spec §9) | Phase-1 baseline | D1 |
|---|---:|---:|
| HipSTR same-length-het cells (common loci) | 289 | 288 |
| **ours recovers as het** | 155 (54 %) | 155 (54 %) |
| undercalled to hom (the tail) | 134 (46 %) | 133 (46 %) |

The 46 % undercall tail does **not** close under D1, and this is the honest, predicted
outcome, not a shortfall: the P1.5 report attributes that tail to the cohort's extreme
apparent inbreeding (`F_IS ≈ 0.82`, whose `ssrCallWarning` D1 still emits) tilting the
**genotyping-time** prior toward homozygosity — a mechanism in the genotyping EM / outer
`F` loop (roadmap E1/E2), **not** the pre-pass. D1 was scoped to de-contaminate the
pre-pass `ε`/`θ` (done, −23 %) and to let the pre-pass *see* same-length hets (done,
0→49); closing the downstream tail is a separate, out-of-scope lever. Spec §6 stated this
in advance, and the measurement confirms it: the pre-pass win is real and the downstream
tail moves independently of it.

Two forces would make the `ε` win larger in absolute terms than the 23 % seen here: a
cohort that is **less homozygote-dominated** (tomato's burn-in is 4580 homs vs 49
same-length hets, so the contaminating reads are a small slice) and one that is **more
interruption-rich**. On tomato the effect is a clear *relative* 23 % with a small
*absolute* magnitude; the mechanism is validated and will scale with same-length-het
prevalence.

## Verdict

**All three D1 exit criteria pass:** ε de-inflation (−23 %, direct win), same-length-het
confident-set placement (0→49, enabling win), and no regression (length 96.5 %→96.5 %,
zygosity 97.3 %→97.3 %, SNP e2e green). The downstream het-recovery tail is unchanged, as
the spec predicted — it is the F_IS genotyping-prior lever (E1/E2), not D1's. D1 is
validated and shippable.

## Reproduce

```
# frozen ε + confident composition (D1 and, at cc1f401, the baseline):
DEV_EXTRA_MOUNT=…/benchmarks ./scripts/dev.sh env \
  SSR_P20_CATALOG=…/results/ours/ssr_tomato1.ssr.catalog \
  SSR_P20_PSP_DIR=…/results_ssr15k/ours/cohort/psp \
  cargo test --lib d1_dump_frozen_chemistry -- --ignored --nocapture

# full VCF + concordance:
ssr-call <psp…> --catalog <cat> --output cohort_d1.ssr.vcf
uv run benchmarks/lib/ssr_concordance.py --ours cohort_d1.ssr.vcf --hipstr <hipstr.vcf.gz>
uv run tmp/ssr_p15/seq_concordance.py cohort_d1.ssr.vcf <hipstr.vcf.gz>
```
