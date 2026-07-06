# SSR interrupted-repeat recall — P2.0 measurement gate (per-allele stutter)

**Date:** 2026-07-06 · **Branch:** `ssr-interruptions` · **Gate outcome:** build the *minimal*
Phase 2 — a **purity → stutter-level** covariate; **not** the per-allele hierarchical shrinkage,
**not** a decay term.

The [spec §6.5](../specs/ssr_interrupted_repeat_recall.md) / [plan P2.0](../implementation_plans/ssr_interrupted_repeat_recall.md)
gate: with Phase 1 calling impure alleles, measure whether purity actually moves the stutter
rate *before* writing any Phase-2 type. Three questions: (i) does purity lower the rate at
**fixed length** (identified free of the length confound, §6.5a); (ii) does it move the
**level** or the **decay** (§6.5b); (iii) is there per-allele residual beyond the covariates.

## Measurement

Two passes. A reads-only proxy first (inferred period, homozygote-fraction filter) suggested a
**weak** effect — but that was an artefact of its length confounding and period inference. The
authoritative pass is **in-caller**: the true catalog period (`locus.motif`) and the caller's own
sequence-aware attribution (`nearest_called_by_sequence`), bucketed per called allele across the
whole `ssr_tomato1` cohort. Infrastructure: the test-only `driver::tests::measure_allele_slips`
(+ the `#[ignore]`d `p20_dump_per_allele_slip_stats`) and `examples/ssr_slip_dump`; analysis in
`tmp/ssr_p15/p20_analyze_precise.py` (not committed). 5027 alleles measured; 3319 with ≥30 reads
and ≥3 units (1971 pure, 1348 impure).

## Results

**(i) Fixed-length purity contrast** — 48 (period, length) cells holding both a pure and an
impure allele with ≥100 reads each side:

- **41/48 cells (85 %): impure slips less than pure** of the same length.
- median within-cell slip rate: **pure 0.064 vs impure 0.018** → **median ratio 0.35**
  (impure alleles slip at ~⅓ the rate of pure alleles of the same length).

**The Simpson's paradox that vindicates §6.5a.** Pooled across all lengths, impure alleles slip
*more* (read-weighted 0.0265 vs 0.0207, ratio 1.28) — but only because impure alleles are
**longer** (interruptions accumulate with tract length → more slip). Control for length and the
sign reverses to ~3× *less*. This is exactly the confound §6.5a predicted would mask the effect,
and why the fixed-length contrast — the identifying variation Phase 1 newly supplies — is the one
that must be read, not the length-marginal.

**(ii) Level vs decay** — the |Δ| slip-size distributions are essentially the same shape (share
of |Δ| ≥ 3: **0.33 pure vs 0.35 impure**; both peak at |1| ≈ 0.54). No decay steepening. Purity
moves the **level**, not the decay — contra the §6.5b hypothesis that an interruption suppresses
multi-unit slips.

**(iii) Residual** — the reads-only pass found negligible within-(period,length,purity) per-allele
spread (~0.5 pp ≈ the effect size). No residual survives to justify a per-allele deviation term.

## Gate decision

| §6.5 question | answer |
|---|---|
| (i) purity separable from length? | **yes, strongly** — ~3× at fixed length (41/48 cells) → earns a covariate |
| (ii) level or decay? | **level** (no decay steepening) |
| (iii) per-allele residual? | **no** → no per-allele deviation term |

**Build the minimal Phase 2:** a **purity → stutter-level** covariate, fit in the pre-pass and
applied per-candidate at the read-model level site. Explicitly **not** the elaborate covariate-
informed per-allele hierarchical shrinkage of §6.3 (no residual to earn it), and **not** a decay/
shape term (no decay effect). This is the "level, not shape" lever §6.3 provisionally proposed —
now confirmed by measurement rather than assumed.

Rationale for building it (vs stopping): a ~3× slip-rate reduction for impure alleles at fixed
length is material for genotype-quality calibration on impure loci; the pooled rate the current
per-locus `θ_locus` sees is length-confounded and mis-serves impure alleles.
