# Implementation report — `ssr-call` genotyping+pre-pass, Milestone C (walking skeleton → checkpoint 1)

**Date:** 2026-06-23 · **Branch:** `ssr-cohort` · **Plan:**
[ssr_call_genotyping_and_parameters.md](../../implementation_plans/ssr_call_genotyping_and_parameters.md)
(Milestone C: C1 candidates, C2 likelihood, C3 prior+seeds, C4 EM+VCF) ·
**Skill:** rust-feature-implementation

## 1. Plan

The genotyping walking skeleton on **supplied** parameters, end to end, reaching
**checkpoint 1** (first VCF / called genotypes match truth):

- **C1** candidate assembly (committed earlier, `19b1c61`).
- **C2** the read likelihood `Qᵣ = Σ_Δ S_θ·Σ_v align_subst` + the genotype likelihood.
- **C3** the `G₀` geometric pseudocount prior + the `π⁰/θ⁰` EM seeds.
- **C4** the per-locus EM → first VCF.

## 2. Assumptions / decisions

- **Q-G2 decided: a slim SSR-specific EM**, not a `posterior_engine` graft. The SSR
  genotype model (size-ploidy candidate multisets, HWE+IBD-`F` prior, `G₀`-regularized
  `π` M-step) is small and direct; with `ε`/`θ`/level fixed for this milestone each
  genotype's data log-likelihood is a one-time precompute and only `π` iterates.
  Bending the SNP engine (chain anchors, class pseudocounts, contamination) would cost
  more than it saves.
- **Diploid only (v1)** — the simulator's ploidy; higher ploidy is a documented
  follow-up (`run_locus_em` asserts ploidy == 2).
- **`θ_locus` M-step deferred** to D (it needs slip attribution accumulators); C4
  holds `θ = θ⁰`. Supplied-parameter genotyping (this milestone) doesn't need it.
- **C2 `align` recomputed, no cache** (Q-G3). The genotype data-likelihood precompute
  makes the EM cheap regardless.
- **`G₀` floor** `1e-12` keeps a far candidate `> 0` (verify-fix #4).
- **VCF is minimal** (`GT:GQ:REPCN`, FILTER from admission); site QUAL = `.`, the
  full FILTER vocabulary + no-call refinement are E2. The contig name is supplied by
  the caller (the work-item carries only a cohort-global id).
- **Reading-layer extension** (in `19b1c61`): `CohortLocus.ref_tract` — the REF allele
  the genotyper + VCF need (the Mi4-anticipated "EM consumes the cohort surface").

## 3. Changes made

- [candidate_set.rs](../../../../src/ssr/cohort/candidate_set.rs) — `assemble_candidates`
  + `CandidateCfg` (C1, `19b1c61`).
- [likelihood.rs](../../../../src/ssr/cohort/likelihood.rs) (new) — `read_likelihood`,
  `read_given_genotype`, `LikelihoodScratch`.
- [allele_freq_prior.rs](../../../../src/ssr/cohort/allele_freq_prior.rs) (new) —
  `g0_pseudocounts`; `Rungs::modal_length` added to `rung_ladder.rs`.
- [em_init.rs](../../../../src/ssr/cohort/em_init.rs) (new) — `LocusSeed`, `seed_locus`
  (π⁰ tally + θ⁰).
- [em.rs](../../../../src/ssr/cohort/em.rs) (new) — `run_locus_em`, `EmCfg`,
  `SampleCall`, `LocusCall`, the diploid genotype EM.
- [vcf_out.rs](../../../../src/ssr/cohort/vcf_out.rs) (new) — `format_vcf_record`.
- [mod.rs](../../../../src/ssr/cohort/mod.rs) — wired the new modules.

## 4. Tests added (this report covers C2–C4; C1's 8 are in `19b1c61`)

- **C2 (6):** faithful read dominated by Δ=0; prefers own allele; −1 stutter explained
  by parent; `Qᵣ ≤ 1`; impure Σ_v runs; genotype mix + outlier floor.
- **C3 (6):** `G₀` decays with distance; floor keeps a far candidate `> 0`; π⁰
  normalized; π⁰ concentrates on the homozygous allele; π⁰ → normalized `G₀` without
  confident samples; θ⁰ uses the period parent shape.
- **C4 (6):** **checkpoint 1** — called genotypes match truth (homs + separated hets)
  end-to-end through the real merger; high-depth clean calls are confident (GQ ≥ 20);
  π normalized; no-call locus reports no-call samples; VCF record formats `GT:GQ:REPCN`
  + FILTER; no-call/filtered VCF record.

## 5. Validation results

- `cargo fmt --check` → pass.
- `cargo clippy --all-targets --all-features -- -D warnings` → pass.
- `cargo test --all-features` → **1231 lib pass** (+24), integration + doctests green.
- ▶ **Checkpoint 1 holds:** simulate → merge → rungs → candidates → seed → EM →
  genotypes equal the simulator's truth at high depth.
- Pre-existing unrelated `benches/psp_writer_perf.rs:386` panic under `--all-targets`.

## 6. Tradeoffs and follow-ups

- `θ_locus` M-step, the `F`/level outer loop (E1), allele-balance + site QUAL (E2),
  and parallelism (F1) are the remaining depth.
- The EM is single-threaded and recomputes `align`; the precompute keeps it cheap.
- The driver still writes the Phase-1 TSV dump; wiring `vcf_out` into a full VCF
  (header + the driver) is E2/F.
- **Checkpoint 1 is a human gate** — the loop halts here for sign-off before D.
