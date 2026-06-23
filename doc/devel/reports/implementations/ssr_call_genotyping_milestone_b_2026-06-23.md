# Implementation report — `ssr-call` genotyping+pre-pass, Milestone B (Shared locus primitives)

**Date:** 2026-06-23 · **Branch:** `ssr-cohort` · **Plan:**
[ssr_call_genotyping_and_parameters.md](../../implementation_plans/ssr_call_genotyping_and_parameters.md)
(Milestone B: B1 rung ladder, B2 stutter kernel, B3 align) · **Skill:** rust-feature-implementation

## 1. Plan

The three real algorithms both Stage-2 halves call:

- **B1 `rung_ladder.rs`** — `build_rungs` (pool the cohort into the length-keyed
  rung ladder) + `resolve_confident_genotype` (the heuristic 1..ploidy-clear-peak
  gate: separation, dosage balance, cohort peak-recurrence, depth floor).
- **B2 `stutter.rs`** — `s_theta` (the kernel `S_θ(Δ) = level·shape(Δ)`),
  `reach_variants` (placement-distinct realizations of `cand ⊕ Δ`),
  `refine_theta_locus` (the per-locus shape M-step shrunk to the prior).
- **B3 `pair_hmm.rs`** — `align_subst` (banded forward, flat-`ε`, in-tract
  substitutions only / gaps in flanks only, exact-match fast path).

## 2. Assumptions (spec gaps resolved here)

- **B1 is the heuristic gate.** The model-based 1-vs-2-peak BIC resolution test
  (Q-P7) is D1, which layers the C2 likelihood on these peaks; B1 delivers the
  threshold-parameterized heuristic the plan specifies for this step.
- **Merged-het detection without a 2nd peak.** A balanced 1-apart het produces *zero*
  clear peaks (the alleles cancel each other's prominence), so "adequate depth but no
  clear maximum" is classified `Unresolved(Merged)`. This is the heuristic signal a
  true homozygote (one dominant peak) does not trip.
- **Rung length = `seq.len() / period` units.** Substitution variants share a length
  (the §6 in-tract-no-indel invariant); exact for the pure tracts the simulator emits.
- **B2 kernel matches the simulator forward model** (the plan §B contract): `s_theta`
  is the scoring counterpart of `sim.rs::slip_length` — faithful mass `1−level`,
  direction split `up:down`, magnitude geometric in `decay` truncated at `MAX_SLIP`
  with the tail absorbed — so `Σ_Δ S_θ = 1` exactly.
- **B3 regimes:** exact match → `(1−ε)^len`; equal length → substitution closed form
  `(1−ε)^match·(ε/3)^mismatch`; unequal length → banded forward with gaps confined to
  a `FLANK_SLOP = 2` window each end (interior indels not competed). `gap = ε` and
  `FLANK_SLOP` are provisional, pinned in F2.
- **Period-1 flagging is deferred to where the period is known** (C2 / `vcf_out`), not
  `align_subst` (which has no period argument) — it is a locus property (spec §6).

## 3. Changes made

- [rung_ladder.rs](../../../../src/ssr/cohort/rung_ladder.rs) — `RungCfg`
  (`dev_default`), `Rungs` (+ accessors), `build_rungs`, `resolve_confident_genotype`,
  `sample_histogram`/`is_clear_peak`/`representative_sequence` helpers.
- [stutter.rs](../../../../src/ssr/cohort/stutter.rs) — `s_theta`,
  `geometric_magnitude`, `Segment`/`segment`/`render`, `reach_variants`,
  `refine_theta_locus`.
- [pair_hmm.rs](../../../../src/ssr/cohort/pair_hmm.rs) (new) — `HmmScratch`,
  `align_subst`, `substitution_product`, `banded_forward`, `in_flank`.
- [mod.rs](../../../../src/ssr/cohort/mod.rs) — wired `pair_hmm`.

## 4. Tests added (22)

- **B1 (7):** rungs pool lengths/support/peak-recurrence; clean homozygote → 1 peak;
  separated het → 2 labelled peaks; 1-apart het → Merged; hom+heavy-stutter →
  DosageInconsistent; non-recurrent allele → NonRecurrent; thin → Thin.
- **B2 (9):** kernel sums to 1; faithful mass = 1−level; contraction bias + decay
  ordering + beyond-MAX = 0; pure allele → 1 variant; Δ=0 → candidate; one
  interruption → 2 placements; un-contractable run → no variant; shrinkage collapses
  to prior with no data; refine moves toward observed direction/decay.
- **B3 (6):** exact match = `(1−ε)^len`; one substitution = exact·(ε/3)/(1−ε); flank
  indel scored; interior indel ≪ flank indel; probability ≤ 1; scratch reuse
  bit-identical.

## 5. Validation results

- `cargo fmt --check` → pass.
- `cargo clippy --all-targets --all-features -- -D warnings` → pass (fixed a
  `type_complexity` lint with a `SeqCount` alias).
- `cargo test --all-features` → **1206 lib pass** (+22), integration + doctests green.
- Pre-existing unrelated `benches/psp_writer_perf.rs:386` panic under `--all-targets`
  (out of scope).

## 6. Tradeoffs and follow-ups

- `reach_variants` run decomposition is greedy left-to-right motif matching; correct
  for the pure + single-interruption cases tested. Pathological multi-interruption
  tracts are exercised in C1/D once real candidates flow.
- `align_subst`'s `gap` cost (= ε) and `FLANK_SLOP` are provisional (F2).
- B1/B2/B3 are single-threaded and pure; parallelism is F1.
- These primitives are now ready for C1 (candidate assembly) and C2 (the likelihood
  `Σ_Δ s_theta · Σ_v align_subst`).
