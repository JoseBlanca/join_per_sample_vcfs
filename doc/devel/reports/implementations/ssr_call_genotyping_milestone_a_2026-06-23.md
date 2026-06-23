# Implementation report — `ssr-call` genotyping+pre-pass, Milestone A (Foundations)

**Date:** 2026-06-23 · **Branch:** `ssr-cohort` · **Plan:**
[ssr_call_genotyping_and_parameters.md](../../implementation_plans/ssr_call_genotyping_and_parameters.md)
(Milestone A: A1 core types, A2 simulator) · **Skill:** rust-feature-implementation

## 1. Plan

Deliver the two A-milestone steps that the rest of Stage-2 builds on:

- **A1 — core types & data model** (nouns only, no logic): the parameter set the
  pre-pass freezes and the EM consumes, the sufficient-statistic accumulators, the
  determinism primitive, and the candidate / confident-genotype / placement shapes.
- **A2 — simulator**: generate synthetic cohorts with **known** genotypes, stutter
  (shape × level), `ε`, and sample groups as the `.ssr.psp` bytes the existing
  reader consumes, plus a queryable truth table — so every later step is a
  recover-what-we-put-in test.

Types seeded in their destination files (plan §2 layout) so milestones B/C add
behavior without moving anything.

## 2. Assumptions (spec gaps resolved here)

- **`MAX_SLIP = 10`** — a provisional value; it is a compile-time array bound for
  `SlipProfile`, so it needs a concrete value now. Recalibrated in F2.
- **`PlacementVariant`** shape is tied to B2's `reach_variants`; defined minimally
  (one realized tract sequence). B2 may extend it.
- **`PeakAllele`** carries `{allele bytes, repeat_len, support}` — the minimum the
  B1 dosage-consistency + cohort-recurrence checks need.
- **Floating types** (`StutterShape`/`StutterLevel`/`PerBaseError`/`ParamSet`/
  `SimChemistry`) derive `PartialEq` but **not** `Eq` (f64).
- **Simulator forward model** (documented in `sim.rs` header, to be matched by the
  B2 kernel): per read — slip with prob `level(units)=baseline+slope·units`; on a
  slip, sign `+` with prob `up/(up+down)`, magnitude geometric with continuation
  prob `decay` (capped at `MAX_SLIP`); then within-tract substitutions at rate `ε`.
  Length clamped to `≥ 1` unit.
- **Simulator is `#[cfg(test)]`** — test/bench-only, no release-binary bloat; if a
  `benches/` harness needs it later, relax the gate behind a `#[doc(hidden)]` seam.
- **Dependency-free PRNG** — an inline SplitMix64 rather than adding `rand`, so the
  simulator's determinism rests on no external crate.

## 3. Changes made

New files under [src/ssr/cohort/](../../../../src/ssr/cohort/):

- [param_estimation.rs](../../../../src/ssr/cohort/param_estimation.rs) — `ParamSet`
  (the frozen pre-pass→EM interface), `PerBaseError`/`StutterShape`/`StutterLevel`/
  `G0PseudocountDecay`/`SampleGroupId`, `SlipProfile`/`SampleStutterStats`
  accumulators, `MAX_SLIP`, and `FixedPointAccum` (the `scale→round→i128`
  order-independent float reduce — the only behavior in A1, because its
  associativity is the stage's determinism guarantee).
- [candidate_set.rs](../../../../src/ssr/cohort/candidate_set.rs) — `CandidateSet` +
  `Admission` (site FILTER reasons).
- [rung_ladder.rs](../../../../src/ssr/cohort/rung_ladder.rs) — `Resolution` /
  `ResolvedGenotype` / `PeakAllele` / `UnresolvedReason`.
- [stutter.rs](../../../../src/ssr/cohort/stutter.rs) — `PlacementVariant`.
- [sim.rs](../../../../src/ssr/cohort/sim.rs) — the A2 simulator: `SplitMix64`,
  `SimChemistry`/`SimGenotype`/`SimLocus`/`SimSample`/`SimCohortSpec`, `simulate`,
  `SimCohort` (with a `merger()` that drives the real reader path), and
  `TruthTable`.
- [mod.rs](../../../../src/ssr/cohort/mod.rs) — wired the four production modules +
  the `#[cfg(test)]` `sim` module.

## 4. Tests added

- `param_estimation`: `ParamSet` round-trip; `SlipProfile` default zeroed;
  `FixedPointAccum` value independent of add order (bit-identical); merge matches
  sequential and is grouping-independent.
- `candidate_set`: `ref_idx` indexes into `alleles`; `Admission` variants distinct.
- `rung_ladder`: confident homozygote = 1 peak; separated het = 2 labelled peaks;
  unresolved reasons distinct.
- `stutter`: `PlacementVariant` holds its sequence.
- `sim`: generated cohort round-trips through the merger (present-sample sparsity,
  motif/coords); modal observation = the homozygote tract; truth table queryable;
  same seed ⇒ byte-identical, different seed diverges; mononucleotide (C1) locus
  round-trips; higher stutter level ⇒ lower faithful fraction.

## 5. Validation results

- `cargo fmt --check` → pass (after `cargo fmt`).
- `cargo clippy --all-targets --all-features -- -D warnings` → pass (clean).
- `cargo test --all-features` → **1181 lib tests pass** (+16 over the pre-A1
  baseline of 1165), integration + doctests green.
- Pre-existing, unrelated failure: `cargo test --all-targets` trips
  `benches/psp_writer_perf.rs:386` (index out of bounds in the bench harness,
  untouched by this work) — flagged for a separate fix.

## 6. Tradeoffs and follow-ups

- The simulator's forward model is the **contract B2's kernel must match**; if they
  diverge, D-milestone recovery tests will catch it — but B2 should be written
  against the `sim.rs` header model deliberately.
- `MAX_SLIP`, the geometric forward parameterization, and the substitution model are
  provisional and get pinned on the simulator in F2.
- A1 type shapes are "shape, not final" where noted (`PlacementVariant`,
  `PeakAllele`); B/C may extend fields.
