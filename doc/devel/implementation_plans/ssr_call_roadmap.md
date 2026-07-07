# SSR Stage-2 (`ssr-call`) — implementation roadmap

**Status:** draft, 2026-06-21, branch `ssr-cohort`. The master step list for building
Stage-2 (`ssr-call`). **Each step below gets its own architecture sub-document**
(`doc/devel/architecture/`) **and implementation plan** (`doc/devel/implementation_plans/`)
written *before* it's coded — both TBD, tracked with the checkboxes.

Design is settled in the spec [ssr_cohort_mark2.md](../specs/ssr_cohort_mark2.md) and the
three architecture sketches ([reading](../architecture/ssr_call_reading.md),
[parameters](../architecture/ssr_call_parameters.md),
[genotyping](../architecture/ssr_call_genotyping.md)); all seven Phase-2 questions were
resolved 2026-06-21 (parameters §9). This roadmap turns that design into build order.

---

## How the order was chosen (principles)

- **Types first, then implementation**, within every step (project rule).
- **Genotyping consumes parameters as an interface; the pre-pass *depends on* the
  genotyping likelihood/EM.** The dependency runs **one way**: genotyping takes the frozen
  parameters (ε, stutter shape, level, `G₀`) as *input*, but the pre-pass is **built on the
  genotyping machinery** — the confident-genotype seed gate (CG-seed) reuses the EM scoring (D1 →
  C2), the soft-EM estimator and the `ε`-freeze check run the genotyper itself (D2/D3 → C2/C4),
  all via the shared primitives (Q-P1). So this is **not** a clean two-way decoupling — the
  pre-pass is downstream of genotyping. That is exactly why we **build supplied-parameter
  genotyping first** (Milestone C) and the pre-pass after (Milestone D): it gets a working
  end-to-end VCF early *and* respects the real dependency (the alternative — building D1/D2
  before C2 exists — is unbuildable, since they call the EM).
- **Walking skeleton early.** Reach a runnable end-to-end `ssr-call` (crude,
  single-threaded, parameters supplied) as soon as possible — Milestone C — then deepen.
- **Simulator-driven.** The simulator (Milestone A) generates cohorts with **known
  truth**, so every later step is testable ("recover what we put in") and the calibration
  numbers have a ground truth.
- **Shared primitives once** (Q-P1): `pool → rungs`, the stutter kernel, `align`, and the
  candidate code are written once and used by *both* the pre-pass and genotyping.
- **Determinism throughout; parallelism last.** Build single-threaded and correct first;
  add the parallel pools — and prove byte-identity across thread counts — once the maths
  is settled.

## Preconditions (already in place)

The **reading layer (Phase 1)** exists as a single-threaded skeleton — locus cursor,
catalog-driven k-way merger, and an end-to-end driver (commits up to `2799005`). It feeds
`CohortLocus` work-items. Its parallelism and the in-memory-vs-reread question (reading
Q-R5) fold into **F1**; the pre-pass's seeded subset can come from one-pass reservoir
sampling over the existing stream, so D-milestone work does **not** block on new reader
features.

---

## The steps

### Milestone A — Foundations

**A1. Core types & data model.**  ☐ arch ☐ plan
Delivers the structs/enums for the whole stage — `CohortLocus`, observed sequences,
candidate/allele types, the parameter set (ε, per-`(group,period)` stutter shape shrunk to a per-period parent `θ_period` — M3, level
coefficients, sample group), seeds, and the sufficient-statistic accumulators. Nouns only,
no logic. *Depends:* reading-layer types. *Source:* spec §2/§4.2; parameters §0/§3.

**A2. Simulator + forward stutter/error model.**  ☐ arch ☐ plan
Generates synthetic cohorts with **known** genotypes, stutter (shape × level), ε, and
sample groups → reads/`.psp` the reader can consume. The generative side of the same model
the kernel scores. Must support **per-group shape** divergence (not just per-group level) so
the M3 per-`(group,period)` axis can be tested for recovery (and the mononucleotide/C1 case).
(Crate/test module, not a subcommand.) *Depends:* A1. *Source:* spec §9 (simulator);
parameters §1.

### Milestone B — Shared locus primitives

**B1. `rungs.rs` — pool → rungs → clear-maxima.**  ☐ arch ☐ plan
The shared peak primitive both phases call (Q-P1); pure, threshold-parameterized.
*Depends:* A1. *Source:* parameters §2; spec §5.

**B2. Stutter kernel `S_θ`.**  ☐ arch ☐ plan
`S_θ(Δ) = level × shape(Δ)`; the geometric shape `(u,d,ρ)`; the cheap re-weight interface.
*Depends:* A1. *Source:* spec §5.2; parameters §0/§1.

**B3. Alignment (pair-HMM), recomputed on demand.**  ☐ arch ☐ plan
Flat-emission `align(obs | cand ⊕ Δ)` with ε: in-tract = **substitution closed-form** +
exact-match fast path (C1), banded HMM for flanks/impure; re-weighted by `S_θ`. **No align
cache in v1** — it's a deferred, measure-first optimization (Q-G3): `align` is a pure
function of `(obs, cand⊕Δ, ε)`, recomputed when needed (this also moots the M2 cache-keying
question). A per-locus memo, then a persistent table, are added later *only if* measurement
shows `align` is the wall. *Depends:* A1, B2. *Source:* spec §6; genotyping §4/Q-G3.

### Milestone C — Genotyping walking skeleton (parameters supplied)

**C1. Candidate assembly (S1) + reachability (S2).**  ☐ arch ☐ plan
Rungs → candidate set (±1 rescue, locus-admission motif filter); per-allele stutter
reachability. *Depends:* B1, B2. *Source:* spec §5/§7; genotyping.

**C2. Read likelihood `Qᵣ` (S3).**  ☐ arch ☐ plan
`Qᵣ(obs|cand) = Σ_Δ S_θ(Δ)·align(…)` (sum-over-slips), `align` recomputed on demand (no
cache, B3/Q-G3). *Depends:* B3, C1.
*Source:* spec §6; genotyping.

**C3. `G₀` base measure + seeds (π⁰/θ⁰).**  ☐ arch ☐ plan
Geometric pseudocount vector per candidate; the per-locus seeds computed **in Phase 3**
(Q-P3) off rungs + the supplied globals. *Depends:* B1, C1. *Source:* spec §4.3/§5.5;
parameters §5/§6.

**C4. Per-locus EM → first end-to-end VCF.**  ☐ arch ☐ plan
The per-locus EM (π + `θ_locus` M-step) reusing `posterior_engine` (the extend-vs-fork
call is made here); a minimal VCF. **Milestone: end-to-end `ssr-call` on simulated data
with *supplied* parameters.** *Depends:* C2, C3. *Source:* spec §5.4/§4.2; genotyping §5.

### Milestone D — Parameter pre-pass (estimate the parameters)

**D1. Confident-genotype resolution test (1..ploidy-peak) + per-locus fitting.**  ☐ arch ☐ plan
The peak-count resolution test (Q-P7, generalized by CG-seed) that gates which (sample, locus)
feed estimation — admitting **confident homozygotes ∪ well-separated hets** (peaks ≥ 2 units
apart, dosage-consistent, each allele cohort-recurrent; the diploid core is the one-vs-two-allele
test, extended to *p* peaks); the per-locus parameter fit off those confident genotypes (a het
contributes two labelled outer skirts; polyploids that won't resolve lean on coded priors).
*Depends:* B1, B3, C2. *Source:* parameters §2; spec §4.3/§4.4 (CG-seed).
The heuristic stand-in shipped in Milestone D; the **BIC form** (the model-based part that
fixes the same-length-het ε contamination) is designed in its own trio — spec
[ssr_bic_confident_genotype.md](../specs/ssr_bic_confident_genotype.md), arch
[../architecture/ssr_bic_confident_genotype.md](../architecture/ssr_bic_confident_genotype.md),
plan [ssr_bic_confident_genotype.md](ssr_bic_confident_genotype.md) — stacked on Phase-1
(`ssr-interruptions`), pending sign-off.

**D2. Burn-in loop + measure → freeze parameters.**  ☐ arch ☐ plan
The adaptive burn-in (seeded batches, frozen-params map, barrier reduce, update — batch 32
fixed) + measure (per-locus distributions → averages + shape diagnostics). The estimator is
the **soft full-cohort EM responsibility reduce**, with the confident-**genotype** gate
(homs ∪ separated hets — CG-seed) as a **seed** (C2). **Stop on the penalized marginal
log-likelihood plateau** (`Δℓ_pen/|ℓ_pen| < tol`; `ℓ_pen` = the E-step normalizer, ~free,
summed by fixed-point integer accum — verify-fix #1), **not** "calls/params don't move" (M4); with
**multi-start** (best `ℓ_pen`, divergent basins flagged). Produces **frozen ε**, the **cohort-per-period shape
parent** `θ_period`, and **per-sample** shape/level estimates (the clustering input). The
outer-loop level seed `level⁰` is **per group** and is fit in D3 *after* clustering fixes the
groups (D2's per-sample lines feed that clustering); it is then refined in E1, not frozen (C2). (The
per-`(group,period)` shape is fit in D3, once clustering fixes the groups — M3.) Seeded
subset via reservoir sampling over the reader's stream. **Milestone: recover known
parameters on simulated data.** *Depends:* D1. *Source:* parameters §3/§4 (Q-P4/Q-P5).

**D3. Sample-group clustering + per-(group,period) shape + ε-freeze check.**  ☐ arch ☐ plan
Distance-based grouping of close neighbours (uncertainty-scaled) + shrinkage (Q-P6); then,
with groups fixed, **fit the per-`(group,period)` stutter shape shrunk to the cohort-per-period
parent** (M3 — the invariance-diagnostic spread falls out here); the per-run ε-freeze
sensitivity check (Q-P2). **Milestone (M3): on a simulator with a deliberately
group-divergent shape, recover the per-group shapes** (and confirm they collapse to the
shared shape when truly invariant). *Depends:* D2, C4 (the ε-check runs calls). *Source:*
parameters §4/§9 (Q-P2/Q-P6); spec §4.4 (M3).

### Milestone E — Genotyping completeness

**E1. Outer loop — per-individual `F_i` + per-group stutter level.**  ☐ arch ☐ plan
The prior-side outer loop: the `F_i` reduce (mean-IBD-responsibility over variable loci,
shrink to cohort mean, hard ceiling 0.99) **and** the per-group stutter-level reduce (soft
per-allele responsibilities, re-fit `level_baseline + level_slope·length`; C2). Both are
**fixed-order barrier reduces** (deterministic; no `align` rebuild). *Depends:* C4, D2 (the
level seed). *Source:* spec §4.4; genotyping §5.

**E2. FP control + full VCF semantics.**  ☐ arch ☐ plan
λ outlier + recurrence + allele-balance/overdispersion; emit-variable, site QUAL =
Phred(variable), SSR FILTER reasons, per-sample no-calls. *Depends:* C4. *Source:* spec
§4.5/§6; genotyping.

### Milestone F — Scale & calibrate

**F1. Parallelism & determinism.**  ☐ arch ☐ plan
The Phase-2 batched map-reduce-with-barrier pool and the Phase-3 streaming workers; the
reading layer's parallelism + Q-R5; **byte-identity across thread counts** proven.
*Depends:* C4, D3. *Source:* parameters §4/§7; reading Q-R5.

**F2. Calibration & validation.**  ☐ arch ☐ plan
Set every calibration number on the simulator (the Q-P5 bundle: defaults, the `ℓ_pen`-plateau
settle tolerance + multi-start count (M4), measurement window, shrinkage strengths, `G₀`
decay; plus the Q-P7 purity tuning). Then the **correctness gates** (M4 — external anchors,
not "calls don't move"): simulator ground-truth recovery, the known-protocol positive
control, multi-start `ℓ_pen` agreement, the discrete-vs-continuum (ε, level) cloud, and the
**batch-size invariance check** (m3 — vary the batch size, confirm the `ℓ_pen` plateau and
recovered chemistry agree; *prove it, don't assert it*). *Depends:* all. *Source:*
parameters §9 "Remaining"; spec §9.

---

## Two integration checkpoints to hold the line

1. **After C4** — a real (if crude) VCF comes out of simulated data with supplied
   parameters. Proves the reading → candidates → likelihood → EM → VCF spine end to end.
2. **After D2** — the pre-pass recovers the simulator's known parameters. Proves the
   estimation half before it's wired in front of genotyping.

Everything after these is depth (F, the F-loop, FP control) and scale (parallelism),
not new spine.
