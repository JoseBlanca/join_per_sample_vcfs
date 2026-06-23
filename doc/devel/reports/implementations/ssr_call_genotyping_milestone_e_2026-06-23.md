# Implementation report — `ssr-call` genotyping+pre-pass, Milestone E (completeness)

**Date:** 2026-06-23 · **Branch:** `ssr-cohort` · **Plan:**
[ssr_call_genotyping_and_parameters.md](../../implementation_plans/ssr_call_genotyping_and_parameters.md)
(Milestone E: E1 the F/level outer loop, E2 FP control + full VCF) · **Skill:** rust-feature-implementation

## 1. Plan

- **E1** — the prior-side outer loop: per-individual `F_i` + per-group stutter level,
  wrapping the per-locus EM over all loci.
- **E2** — FP control + output semantics: the allele-balance defence, emit-iff-variable,
  site QUAL, SSR FILTER reasons, the apparent-`F_IS` warning.

## 2. Assumptions / decisions

- **EM surgery (minimal ripple):** `run_locus_em_with(f_per_present, level_per_group)`
  is the new core (emits per-sample `posterior_hom`); `run_locus_em` is a thin wrapper,
  so C4 call sites are unchanged. `LocusCall` gains `posterior_hom`.
- **`F_i` estimator = excess homozygosity** (posterior homozygosity vs HWE `Σ π_i²`
  over variable loci), accumulated via `FixedPointAccum` (order-independent), shrunk to
  the cohort mean, clamped to `[0, F_CEILING = 0.99]`.
- **Level refit = hard-attribution** (each read to its called allele); the soft
  per-allele responsibility reduce is the deferred refinement.
- **Allele-balance defence** = the minor-allele read fraction of a het; below the floor
  the het is a depth-inflated false het (a hom + stutter/error — the SNP-path blind
  spot), so its GQ is scaled down and it no-calls if it drops below `no_call_gq`. This
  is the headline FP control.
- **Site QUAL is a posterior proxy** — `−10·log₁₀ Π_s P(s hom-modal)`. The exact-AF
  convolution kernel (verify-fix #7b) is an F2 refinement; documented as such.
- **`F_IS` warning** fires when the cohort mean `F` exceeds a threshold (a data
  artifact, not real inbreeding).

## 3. Changes made

- [inbreeding.rs](../../../../src/ssr/cohort/inbreeding.rs) (new) — `run_cohort_em`,
  `reduce_f`, `reduce_level`, `OuterCfg`, `CohortCalls` (E1).
- [em.rs](../../../../src/ssr/cohort/em.rs) — `run_locus_em_with` + per-sample `F` +
  `posterior_hom`; `run_locus_em` wrapper; `LocusCall.posterior_hom`;
  `SampleCall::no_call` made `pub(crate)`.
- [vcf_out.rs](../../../../src/ssr/cohort/vcf_out.rs) — `FpControlCfg`, `allele_balance`,
  `apply_fp_control`, `is_variable`, `site_qual`, `f_is_warning`; `format_vcf_record`
  now emits QUAL (E2).
- [mod.rs](../../../../src/ssr/cohort/mod.rs) — wired `inbreeding`.

## 4. Tests added (8)

- **E1 (3):** `F` separates inbred (homozygous-everywhere) from outbred
  (het-everywhere) samples across 16 loci; the loop is deterministic; the level refit
  corrects a deliberately-wrong seed (0.25 → ≈ 0.06).
- **E2 (5):** allele-balance flags an imbalanced het / keeps a balanced one; FP control
  no-calls a depth-inflated false het; `is_variable` true with an ALT / false all-ref;
  site QUAL high for a confident variant, low for monomorphic; the `F_IS` warning fires
  only when mean `F` is implausible.

## 5. Validation results

- `cargo fmt --check` → pass · `cargo clippy --all-targets --all-features -- -D warnings`
  → pass.
- `cargo test --all-features` → **1250 lib pass** (+11 over Milestone D), integration +
  doctests green.
- Pre-existing unrelated `benches/psp_writer_perf.rs:386` panic under `--all-targets`.

## 6. Tradeoffs and follow-ups

- The soft per-allele level responsibility reduce, the exact-AF site QUAL kernel, and
  overdispersion (beta-binomial) beyond the allele-balance heuristic are documented
  refinements (F2 / future).
- The driver still emits the Phase-1 TSV dump; wiring `run_cohort_em` + `apply_fp_control`
  + `format_vcf_record` into a full VCF (header + per-locus records) is the remaining
  output wiring (F / a small driver step).
- **Milestone F** (parallelism + byte-identity, calibration) is the last milestone.
