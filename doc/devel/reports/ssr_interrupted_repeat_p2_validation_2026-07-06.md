# SSR interrupted-repeat recall — Phase 2 validation (P2.3)

**Date:** 2026-07-06 · **Branch:** `ssr-interruptions` · **Verdict:** Phase 2 shipped correct +
guarded; **inert on `ssr_tomato1`** (fits factor 1.0 = Phase 1). Owner-accepted.

Validation of Phase 2 (per-allele purity → stutter-**level**, the minimal form the P2.0 gate
earned) against the [plan P2.3](../implementation_plans/ssr_interrupted_repeat_recall.md) exit
criterion: *calibration on impure loci improves, with no regression on pure loci and no loss of
the Phase 1 recall gain.* Phase 2 is commits P2.1–P2.2b (`84d24f1`…`4563665`).

## What was built

- **Fit** (`prepass.rs::fit_purity_level`, P2.2b): a cohort-global `per_interruption_factor`,
  read-weighted least-squares of `ln(impure_rate/pure_rate)` on interruption count through the
  origin, contrasting pure vs impure **at fixed length** (§6.5a) over confident-genotype slip
  stats. Guarded: capped at 1 (an interruption never *amplifies* stutter), floored at 0.05,
  `none()` (Phase 1) with no clean contrast; a fixed-order reduce (byte-identical across threads).
- **Apply** (`em.rs::candidate_level`, P2.2a): the per-candidate level is scaled by
  `purity.level_factor(interruption_count(candidate))`, at both scoring sites (data likelihood +
  allele-balance deconvolution). A pure allele or neutral factor is the pre-Phase-2 arithmetic.

Both are unit-proven: `fit_purity_level_recovers_a_lower_impure_slip_rate` (a 0.05/0.20 contrast
→ factor 0.25), `candidate_level_scales_by_purity_only_for_impure_alleles`, and the neutrality /
guard tests.

## Result on `ssr_tomato1` — INERT

Re-ran `ssr-call` (Phase-2 binary) on the 51-sample cohort. The VCF is **byte-identical to the
Phase-1 VCF** — the pre-pass fit `factor = 1.0`. Instrumented, the fit sees the data (14 pure
cells, 48 impure cells, 2659 impure reads) but returns 1.0 because the confident-genotype
contrast is not `impure < pure`: at a relaxed read floor the slope is **+0.36** (impure appears
to slip *more*), and the cap-at-1 guard correctly refuses it.

So the exit criterion is **partially met**: **no regression, no loss of Phase-1 recall**
(identical output), but **no calibration improvement** — because the fit honestly finds no
robust purity effect in its data source and declines to invent one.

## The population discrepancy (the substantive finding)

| measurement | population | attribution | result |
|---|---|---|---|
| reads-only proxy | all, hom-dominant | length + inferred period | weak |
| **P2.0** | all called alleles (incl. hets) | **sequence** (`nearest_called_by_sequence`) | impure **0.35×** (strong) |
| **P2.3 pre-pass fit** | **confident homozygotes** (cleanest stutter data) | length (`nearest_parent`) | no clean effect / wrong-signed |

The confident-homozygote population — the cleanest stutter measurement, with no attribution
ambiguity and no het contamination — does **not** confirm the purity→level effect P2.0 measured
on the broad caller-attributed population. Either the confident population is genuinely
different, or P2.0's 0.35× was partly an artifact of sequence-attribution / het inclusion in the
broad population. Unresolved; the three measurements disagree and the cleanest is the least
supportive.

## Decision (owner, 2026-07-06)

**Ship Phase 2 as-is — correct, guarded, and inert on tomato.** Forcing it active would apply a
wrong-signed estimate from the confident population; the guards correctly prevent that. Phase 2
activates only on a cohort whose confident-genotype data carries a clean fixed-length purity
contrast; the synthetic fixtures prove the machinery. The population discrepancy and the option
to **refit the factor from the broad caller genotypes** (the P2.0 population) are tracked in
`doc/devel/TODO.txt` as the follow-up if impure-loci GQ calibration becomes a priority.

## Verdict

Phase 2 is complete, correct, and safely guarded. Its exit criterion is met on the *no-regression
/ no-recall-loss* axes and honestly not on the *calibration-improvement* axis on this cohort,
because the effect is not recoverable from the confident-genotype pre-pass here. The
interrupted-repeat recall feature is **done: Phase 1 delivers the recall (validated, shipped);
Phase 2 is the correct, guarded refinement that stays inert until the data supports it.**
