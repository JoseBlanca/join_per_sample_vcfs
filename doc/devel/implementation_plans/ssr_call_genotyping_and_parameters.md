# SSR Stage 2 — `ssr-call` genotyping + parameter pre-pass (implementation plan)

**Status:** draft, 2026-06-23, branch `ssr-cohort`. **One plan for Phases 2 + 3 of
`ssr-call`**, deliberately fused because they **share most of their code** (Q-P1/Q-G1):
the pre-pass that estimates the chemistry parameters is *built on* the genotyping
likelihood/EM, so building them apart would mean building the pre-pass against a
likelihood that does not yet exist. This plan covers everything **after** the reading
layer — candidate assembly, the read likelihood, the per-locus EM, the cohort VCF
(Phase 3) **and** the burn-in/measure/cluster pre-pass that freezes `ε` / fits the
stutter shape + level + `G₀` (Phase 2).

**Specs & architecture this implements:**
- [spec ssr_cohort_mark2.md](../specs/ssr_cohort_mark2.md) — §4.2/§4.3/§4.4/§4.5,
  §5/§6/§7 (model framing, seeding, shared-parameter estimation, output).
- [architecture: ssr_call_genotyping.md](../architecture/ssr_call_genotyping.md) — Phase 3
  module/struct shape (authoritative for the calling half).
- [architecture: ssr_call_parameters.md](../architecture/ssr_call_parameters.md) — Phase 2
  module/struct shape (authoritative for the pre-pass half; **§0 is the canonical
  vocabulary** — read this plan against it).
- [roadmap: ssr_call_roadmap.md](ssr_call_roadmap.md) — the A1→F2 step list this plan
  fleshes out. **Milestone/step IDs below are the roadmap's** (A1, B2, C4, D2, …) so the
  checkboxes there track this plan.

**Companion plan (done):** [ssr_call_reading.md](ssr_call_reading.md) — Phase 1, the
restartable `CohortLocus` merge stream this plan consumes. It is a **precondition**
(below), not part of this plan.

**Build philosophy** (project rules + roadmap principles):
- **Types first, then implementation**, within every step.
- **Walking skeleton early** — a runnable end-to-end `ssr-call` on *supplied* parameters
  (Milestone C) before any estimation; then the pre-pass in front of it (Milestone D).
- **Shared primitives written once** (`rung_ladder`, `stutter`, `pair_hmm`, candidate code,
  `allele_freq_prior`, `em_init`) and called by *both* halves.
- **Single-threaded & correct first; parallelism + byte-identity last** (Milestone F).
- **Simulator-driven** — known truth (Milestone A) makes every later step a
  recover-what-we-put-in test.
- Each milestone compiles, tests, and is reviewable on its own; **pause between
  milestones**.

---

## 0. What exists to build on (verified surfaces)

| surface | path | what we use |
|---|---|---|
| `CohortLocus` / `SampleEvidence` / `LocusId` / `SsrQc` | [src/ssr/cohort/types.rs](../../../src/ssr/cohort/types.rs) | the Phase-1 work-item; `samples[k]`↔`present[k]`, `seq_counts: Vec<(Box<[u8]>, u32)>`, `motif`, `ref_frame` (tract+flanks) |
| restartable merge + driver | [src/ssr/cohort/merge.rs](../../../src/ssr/cohort/merge.rs), [driver.rs](../../../src/ssr/cohort/driver.rs), [reader.rs](../../../src/ssr/cohort/reader.rs) | catalog-driven `CohortLocus` stream, monotonic locus-seq, `restart()` for the two-pass (pre-pass then genotyping); producer→bounded-queue→worker→reorder-by-seq topology (currently a stub worker) |
| `Motif` / `Locus` | [src/ssr/types.rs](../../../src/ssr/types.rs) | `Motif::period()`, `as_bytes()`; `Locus::ref_tract()`/`left_flank()`/`right_flank()`/`purity_fraction()` |
| Stage-1 SSR pair-HMM scratch pattern | [src/ssr/pileup/alignment.rs](../../../src/ssr/pileup/alignment.rs) | `HmmModel`, `ViterbiScratch` (banded rolling rows, `resize(m,n)`). **NB: Stage-1 is Viterbi (max-path) — Phase 3 needs a *forward sum* (probability), so `pair_hmm.rs` is new code reusing the banded/scratch *pattern*, not the function.** |
| `posterior_engine` EM + IBD-`F` prior | [src/var_calling/posterior_engine.rs](../../../src/var_calling/posterior_engine.rs) | reuse the E-step + `F·π + (1−F)·πᵖ` prior + `log_f_per_sample` vector; **replace** class pseudocounts (`ref/snp/indel/compound_pseudocount`, `classify_allele`) with the per-candidate `G₀` vector; **bypass** SNP-only machinery (chain anchors, contamination) |
| exact-AF QUAL kernel | [posterior_engine.rs:3010 `compute_qual_via_exact_af`](../../../src/var_calling/posterior_engine.rs#L3010), [:3132 `convolve_ac_linear`](../../../src/var_calling/posterior_engine.rs#L3132) | the convolution kernel for site QUAL — **generalize off allele-0** (parametrize the collapsed allele, per-`i` Beta from `G₀`), don't call unchanged (verify-fix #7b) |
| allele-balance / depth-inflation defense | [src/vcf/qual_refine.rs](../../../src/vcf/qual_refine.rs) | pattern for the §4 allele-balance penalty on deconvolved per-allele responsibilities |
| CLI / VCF shell | [src/pop_var_caller/cli.rs](../../../src/pop_var_caller/cli.rs), [src/vcf/](../../../src/vcf/) | `SsrCallArgs` already registered (reading layer); extend with param-source flags; reuse REF/ALT/REPCN/BPDIFFS/`prune_unsupported_alleles`/header detection |

**Does not exist yet:** any of `rung_ladder.rs`, `stutter.rs`, `pair_hmm.rs`, `likelihood.rs`,
`candidate_set.rs`, `allele_freq_prior.rs`, `em_init.rs`, `em.rs`, `prior.rs`, `inbreeding.rs`,
`prepass.rs`, `sample_groups.rs`, `vcf_out.rs`; the simulator; the parameter/accumulator
types. The Phase-1 worker is a stub we replace.

---

## 1. The dependency that sets the order (why one plan, this sequence)

The pre-pass does the **same per-(sample, locus) work** as genotyping (pool → rungs →
peaks → align reads vs. candidates with `S_θ` + `ε`) and differs only in the **reduce**:
genotyping → EM → calls; the pre-pass → sufficient statistics → params. So the dependency
is **one-way and downstream**: genotyping takes the frozen params as an *input interface*,
but the pre-pass is *built on* the genotyping scoring (the confident-genotype gate scores
the 1-vs-2-peak models with the EM machinery; the soft estimator **is** the genotyper's
own EM run on a subset). Therefore:

```
A (types + simulator)
  └─ B (shared primitives: rungs, stutter, pair_hmm)
       └─ C (genotyping on SUPPLIED params → first VCF)   ◄── checkpoint 1
            └─ D (pre-pass: estimate the params, on top of C's EM)   ◄── checkpoint 2
                 └─ E (outer F + level loop, FP control, full VCF)
                      └─ F (parallelism + byte-identity, calibration)
```

Building D1/D2 before C2 exists is **unbuildable** (they call the likelihood/EM). The
walking skeleton (C) also gets us a real, if crude, VCF early.

---

## 2. Module layout (combined; under the existing `src/ssr/cohort/`)

```
src/ssr/cohort/
  # ── Phase 1, exists ──
  types.rs        # CohortLocus, SampleEvidence, LocusId, SsrQc   (extend in A1)
  reader.rs merge.rs driver.rs mod.rs test_support.rs

  # ── A1 — parameter & accumulator types ──
  param_estimation.rs # PerBaseError (ε), StutterShape (θ_period/θ_group/θ_locus), StutterLevel,
                      #   G0PseudocountDecay, SampleGroupId; the frozen ParamSet the EM consumes;
                      #   SlipProfile / SampleStutterStats accumulators; FixedPointAccum
  sim.rs              # A2 simulator (test/bench module, NOT a subcommand)

  # ── B — shared locus primitives (both halves call) ──
  rung_ladder.rs      # pool → rungs → clear-maxima; confident-genotype RESOLUTION (1..ploidy
                      #   peaks, separation+dosage+recurrence) → ResolvedGenotype | Unresolved
  stutter.rs      # S_θ(Δ) = level × shape(Δ); reachability ⊕ (placement-variant SET for
                  #   impure A⊕Δ); θ_locus M-step; shape shrinkage θ_locus←(group,period)←θ_period
  pair_hmm.rs     # align_subst(obs | v): banded forward, flat-ε emission, in-tract
                  #   substitutions-only, gaps in flanks only; exact-match fast path

  # ── C — genotyping (supplied params) ──
  candidate_set.rs# S1: nominate + union + ref-seed + locus-admission motif filter → CandidateSet
  likelihood.rs   # S3: Qr(obs|cand)=Σ_Δ S_θ(Δ)·align(obs|cand⊕Δ); align = Σ_v over placements
                  #   + λ outlier + allele-balance term
  allele_freq_prior.rs # G₀ geometric pseudocounts: per-loci-group decay fit + per-candidate vector
  em_init.rs         # π⁰ tally (putative genotypes + G₀) and θ⁰ (globals × locus length); F⁰ default
  prior.rs        # thin adapter onto posterior_engine IBD-F prior; inject per-candidate G₀ vector
  em.rs           # per-locus EM: engine E-step + G₀ M-step + θ_locus M-step

  # ── D — parameter pre-pass ──
  prepass.rs      # burn-in (seeded batch → frozen-param map → barrier reduce → update) →
                  #   measure (stratified) → freeze ε / shape parent / per-sample lines;
                  #   ℓ_pen plateau stop; coded-default fallback + m2(a) zero-flag
  sample_groups.rs# cluster per-sample (ε, level) → sample groups (uncertainty-scaled distance,
                  #   BIC group count); per-(group,period) shape fit; ε-freeze BIC check

  # ── E — completeness ──
  inbreeding.rs   # outer loop: F_i reduce + per-group level reduce (both fixed-point i128)
  vcf_out.rs      # CALL: emit-variable, P(monomorphic) QUAL, SSR FILTER/GQ, F_IS warning
```

**Module names are canonical in *this plan*** (build against these). Three shared files and
one new types file use clearer, self-explaining names than the arch §8 sketches; the arch
docs remain authoritative for *struct/signature shape*, this plan for *file names*. Map when
cross-reading:

| this plan (canonical) | parameters arch §8 | genotyping arch §8 | note |
|---|---|---|---|
| `rung_ladder.rs` | `rungs.rs` | `rungs.rs` | shared (Q-P1/Q-G1) |
| `allele_freq_prior.rs` | `base_measure.rs` | — | shared; the `G₀` Dirichlet base measure |
| `em_init.rs` | `seed.rs` | — | the `π⁰/θ⁰/F⁰` seeds (computed in Phase 3, Q-P3) |
| `param_estimation.rs` | (types scattered) | (types scattered) | **new** — consolidates the A1 types/accumulators in one file |
| `stutter.rs` `pair_hmm.rs` `candidate_set.rs` `likelihood.rs` `prior.rs` `em.rs` `prepass.rs` `sample_groups.rs` `inbreeding.rs` `vcf_out.rs` | same | same | unchanged |

---

## 3. Core types & signatures (the load-bearing sketches)

> All sketches are **shape, not final**. `pub(crate)` throughout (internal module).

### 3.1 Parameters — the interface between the two halves (`param_estimation.rs`, A1)

```rust
/// A data-driven cluster of samples sharing chemistry+provenance (params §0).
pub(crate) struct SampleGroupId(pub u16);

/// Per-base within-tract SUBSTITUTION rate (C1). Frozen per sample group.
pub(crate) struct PerBaseError(pub f64);   // ε

/// Stutter SHAPE — where a slip lands. Period-keyed parent, per-(group,period) cell,
/// per-locus refinement, all the same geometric form (u/d rates + decay ρ).
pub(crate) struct StutterShape {
    pub up_rate: f64,     // mass on +Δ
    pub down_rate: f64,   // mass on −Δ (contraction bias ⇒ usually > up_rate)
    pub decay: f64,       // geometric decay ρ per extra unit
}

/// Stutter LEVEL — how often a slip happens; per sample group, linear in repeat length.
pub(crate) struct StutterLevel { pub baseline: f64, pub slope: f64 }  // level(L)=baseline+slope·L

/// G₀ geometric pseudocount decay, per loci group (= period in v1).
pub(crate) struct G0PseudocountDecay { pub p: f64 }   // pseudocount(Δ) = max(p^|Δ|, FLOOR)

/// The FROZEN output of the pre-pass that genotyping consumes (the §1 interface).
/// `level` is the SEED only — refined in the outer loop (E1), not frozen (C2).
pub(crate) struct ParamSet {
    pub error_per_sample_group: Vec<PerBaseError>,      // ε, frozen per sample group
    pub stutter_shape_parent: HashMap<u8 /*period*/, StutterShape>,        // θ_period
    pub stutter_shape_by_cell: HashMap<(SampleGroupId, u8), StutterShape>, // θ_(group,period)
    pub level_seed:  Vec<StutterLevel>,                 // per sample group (seed → E1)
    pub pseudocount_decay_per_loci_group: HashMap<u8 /*period = loci group*/, G0PseudocountDecay>,
    pub group_of_sample: Vec<SampleGroupId>,            // sample idx → its sample group
    pub f0_seed:     f64,                               // F⁰ seed (estimated in E1)
}
```

### 3.2 Sufficient-statistic accumulators (`param_estimation.rs`, A1; arch parameters §3)

```rust
const MAX_SLIP: usize = 10;   // provisional. It is a compile-time ARRAY BOUND (below), so A1
                              //   must pick a concrete value now; F2 re-calibrates the number.

/// Per (group, period) — the conditional profile of a slip. Integer counts ⇒
/// commutative reduce ⇒ order-independent (the determinism trick).
pub(crate) struct SlipProfile { pub down: [u64; MAX_SLIP], pub up: [u64; MAX_SLIP] }

/// Per sample — feeds (a) the per-group level line, (b) the clustering, and
/// (c) — once groups fixed — the (group,period) shape key.
pub(crate) struct SampleStutterStats {
    pub by_length: Vec<(u16 /*repeat len*/, u64 /*faithful*/, u64 /*slipped*/)>,
    pub base_match: u64, pub base_mismatch: u64,   // → ε (substitutions, C1)
    pub read_depth: u64,                           // clustering PRECISION weight (1/√depth) only,
                                                   //   never a feature; dup-free by construction (verify-fix #5)
}

/// Order-independent float reduce (the M1/verify-fix-#1 trick): scale → round → i128 sum.
/// Used for level/F responsibilities AND the decision floats (ℓ_pen, BIC log-liks).
pub(crate) struct FixedPointAccum { acc: i128 }     // SCALE = 2^40
impl FixedPointAccum {
    pub fn add(&mut self, x: f64) { self.acc += (x * SCALE).round() as i128; }
    pub fn merge(&mut self, o: &Self) { self.acc += o.acc; }   // associative+commutative
    pub fn value(&self) -> f64 { self.acc as f64 / SCALE }
}
```

### 3.3 Candidates (`candidate_set.rs`, C1; arch genotyping §2)

```rust
pub(crate) enum Admission { Pass, NotPeriodic, TooManyAlleles, LowDepth }

pub(crate) struct CandidateSet {
    pub alleles: Vec<Box<[u8]>>,   // candidate tract sequences (each an independent allele)
    pub ref_idx: usize,            // reference allele, seeded unconditionally
    pub admit:   Admission,        // → site FILTER (§6)
}
```

### 3.4 Shared primitive signatures (`rung_ladder.rs` / `stutter.rs` / `pair_hmm.rs`, B)

```rust
// rung_ladder.rs — pool → rungs → clear local maxima (Q-P1: ONE peak definition for both phases)
pub(crate) struct Rungs { /* length-keyed; each rung holds a set of sequences */ }
pub(crate) fn build_rungs(locus: &CohortLocus, cfg: &RungCfg) -> Rungs;

// confident-genotype resolution (CG-seed, params §2) — used by the pre-pass gate
pub(crate) enum ResolvedGenotype {
    Peaks(Vec<PeakAllele>),   // 1..=ploidy peaks, each with its labelled parent allele
}
pub(crate) enum Resolution { Confident(ResolvedGenotype), Unresolved(UnresolvedReason) }
pub(crate) fn resolve_confident_genotype(
    sample: &SampleEvidence, rungs: &Rungs, params: &ParamSet, ploidy: u8, cohort_recurrence: &Recurrence,
) -> Resolution;   // Unresolved = merged(<2 apart) | dosage-inconsistent | non-recurrent | thin

// stutter.rs — the kernel and reachability
pub(crate) fn s_theta(delta: i32, shape: &StutterShape, level: f64) -> f64;  // = level × shape(Δ)
/// Placement-distinct variants of `cand ⊕ Δ` (runs between interruptions; pure ⇒ 1).
pub(crate) fn reach_variants(cand: &[u8], motif: &Motif, delta: i32, out: &mut Vec<PlacementVariant>);
/// θ_locus M-step: regularized shape update shrunk to the (group,period) prior.
pub(crate) fn refine_theta_locus(stats: &LocusSlipStats, prior: &StutterShape, strength: f64) -> StutterShape;

// pair_hmm.rs — one fixed placement variant, FORWARD probability (new; reuses banded scratch)
pub(crate) fn align_subst(obs: &[u8], variant: &[u8], eps: f64, scratch: &mut HmmScratch) -> f64;
//   exact-match fast path: obs == variant ⇒ (1−ε)^len  (no DP); else banded forward,
//   in-tract substitutions only, per-base gaps confined to flanks.
```

### 3.5 The read likelihood (`likelihood.rs`, C2 / S3) — the involved core, pseudocoded

```rust
/// Qr(obs_seq | candidate) = Σ_Δ S_θ(Δ) · align(obs | candidate ⊕ Δ)   (arch genotyping §4)
/// align(obs | cand⊕Δ) = Σ_{v ∈ placements(cand,Δ)} Pr(v) · align_subst(obs | v)
pub(crate) fn read_likelihood(
    obs: &[u8], cand: &[u8], motif: &Motif,
    shape: &StutterShape, level: f64, eps: f64,
    scratch: &mut Scratch,
) -> f64 {
    let mut q = 0.0;
    for delta in -(MAX_SLIP as i32)..=(MAX_SLIP as i32) {
        let s = s_theta(delta, shape, level);
        if s == 0.0 { continue; }
        // Σ_v over run-collapsed placement variants (pure ⇒ 1 term; impure ⇒ ~2–3).
        reach_variants(cand, motif, delta, &mut scratch.variants);
        let pr_v = 1.0 / scratch.variants.len() as f64;        // uniform position prior
        let mut align = 0.0;
        for v in &scratch.variants {
            align += pr_v * align_subst(obs, v, eps, &mut scratch.hmm);
        }
        q += s * align;
    }
    q
}

/// Genotype likelihood with the two FP-defense terms (arch genotyping §4):
///   P(read | G) = (1−λ)·[mix over G's alleles of Qr] + λ·(1/D)
/// allele-balance penalty runs LATER, on deconvolved E-step responsibilities (qual_refine pattern).
```

### 3.6 The per-locus EM (`em.rs`, C4) — reuse engine, replace base measure

```rust
pub(crate) struct LocusEm<'a> { engine: &'a PosteriorEngine, /* … */ }
pub(crate) struct LocusCall { pub genotypes: Vec<Genotype>, pub posteriors: Vec<f64>,
                              pub pi: Vec<f64>, pub theta_locus: StutterShape, pub l_pen: f64 }

/// One locus → calls. Reuses posterior_engine E-step + IBD-F prior; replaces the class
/// pseudocounts with the per-candidate G₀ vector; adds the θ_locus M-step.
pub(crate) fn run_locus_em(
    locus: &CohortLocus, cands: &CandidateSet, params: &ParamSet,
    f_per_sample: &[f64], level_per_group: &[StutterLevel],
    seed: &LocusSeed,           // π⁰/θ⁰ from em_init.rs
    cfg: &EmCfg,
) -> LocusCall;
//  loop {
//    E-step: prior(F_i·π + (1−F_i)·π^p) × Qr  → responsibilities   (align recomputed; cache deferred)
//    M-step π:  expected_counts + G₀ pseudocount  (floored > 0, log-domain)
//    M-step θ_locus:  refine_theta_locus(...)  → re-weight S_θ (no align rebuild)
//  } until Δℓ_pen non-increasing / tol
```

---

## 4. Milestone phases (build order, with per-step deliverable + tests)

Each step = roadmap ID. **Pause + review between milestones.**

### Milestone A — Foundations

**A1. Core types & data model** *(param_estimation.rs; extend types.rs)*
Deliver §3.1–§3.4 structs — `ParamSet`, `PerBaseError`/`StutterShape`/`StutterLevel`/
`G0PseudocountDecay`, `SampleGroupId`, the accumulators, `FixedPointAccum`, `CandidateSet`,
`Admission`, the `rung_ladder`/`stutter` enums (`ResolvedGenotype`, `PlacementVariant`).
Nouns only, no logic. *Depends:* reading types. *Tests:* construction/round-trip;
`FixedPointAccum` add/merge associativity (permute → identical `value()`); `s_theta`/
`level()` aren't here yet (B2).

**A2. Simulator + forward stutter/error model** *(sim.rs; test/bench module)*
Generate synthetic cohorts with **known** genotypes, stutter (shape × level), `ε`, and
sample groups → `.ssr.psp` the reader consumes. Must support **per-group shape
divergence** (not just per-group level) so D3/M3 recovery is testable, plus the
mononucleotide/C1 case and well-separated-het cohorts (CG-seed). *Depends:* A1, Stage-1
`.ssr.psp` fixture helpers. *Tests:* a generated cohort round-trips through the reader to
`CohortLocus`; recorded truth table is queryable for later recovery asserts.

### Milestone B — Shared locus primitives

**B1. `rung_ladder.rs`** — `build_rungs` + `resolve_confident_genotype` (1..ploidy peaks;
separation ≥2 units, dosage-consistent heights, cohort recurrence ≥ *k*, depth floor).
Pure, threshold-parameterized. *Depends:* A1. *Tests:* homozygote → 1 peak; separated het
→ 2 peaks both labelled; 1-apart het → `Unresolved(Merged)`; hom+heavy-stutter → rejected
by dosage/recurrence; thin → `Unresolved(Thin)`.

**B2. Stutter kernel `S_θ`** *(stutter.rs)* — `s_theta`, `reach_variants` (placement-
variant set: pure ⇒ 1, impure ⇒ runs-between-interruptions), `refine_theta_locus`, and
the `θ_locus ← (group,period) ← θ_period` shrinkage. *Depends:* A1. *Tests:* pure allele
⇒ single variant; one interruption ⇒ ~2; kernel normalizes; shrinkage collapses to parent
with no data.

**B3. Alignment forward `align_subst`** *(pair_hmm.rs)* — banded forward, flat-`ε`
emission, **in-tract substitutions only / gaps in flanks only**, exact-match fast path
returning `(1−ε)^len`. **New code** reusing the Stage-1 banded-scratch *pattern* (Stage-1
is Viterbi/max — this is a sum). *Depends:* A1, B2. *Tests:* exact match = `(1−ε)^len`;
one substitution scales by `ε/3`-ish; flank indel scored, in-tract indel **not** competed;
period-1 flagged.

### Milestone C — Genotyping walking skeleton (parameters supplied) → **checkpoint 1**

**C1. Candidate assembly (S1) + reachability (S2)** *(candidate_set.rs, stutter.rs)* —
rungs → per-sample nomination (clear maxima ≥ ploidy → top-ploidy; < ploidy → add ±1
neighbours) → union → ref-seed → locus-admission motif filter (adjacent-rung-Δ mode ==
motif length, else `NotPeriodic`) → cap (`MAX_CANDIDATE_ALLELES`, else `TooManyAlleles`).
Impure peaks kept first-class. *Depends:* B1, B2. *Tests:* merged-het rescue adds the
neighbour; non-periodic locus → `NotPeriodic`; cap → `TooManyAlleles`.

**C2. Read likelihood `Qᵣ` (S3)** *(likelihood.rs)* — §3.5: `Σ_Δ S_θ·align`, `align`=`Σ_v`,
λ outlier term. `align` **recomputed on demand** (no cache — Q-G3). *Depends:* B3, C1.
*Tests:* pure-allele single-term path; impure `Σ_v` ~2–3 terms; clean read dominated by
its own faithful term; `Σ_Δ` sums to a probability ≤ 1.

**C3. `G₀` allele-freq prior + EM seeds** *(allele_freq_prior.rs, em_init.rs)* — geometric pseudocount
vector per candidate (centred on cohort modal allele, floored `> 0` / log-domain,
verify-fix #4); per-locus `π⁰` tally (putative genotypes + `G₀`) and `θ⁰` (globals ×
locus length) — **computed in Phase 3 at EM start** (Q-P3), fused with the rungs pass.
`F⁰` = supplied default. *Depends:* B1, C1. *Tests:* far candidate keeps `π⁰ > 0`; no
confident samples ⇒ `π⁰` = normalized `G₀`; underflow floor holds for long tracts.

**C4. Per-locus EM → first end-to-end VCF** *(em.rs, prior.rs, minimal vcf_out.rs)* —
§3.6: engine E-step + IBD-`F` prior + **per-candidate `G₀` M-step** + `θ_locus` M-step;
**the extend-vs-fork engine decision (Q-G2) is made and recorded here.** Minimal VCF
(REF/ALT/REPCN/GT). *Depends:* C2, C3. **▶ Checkpoint 1: end-to-end `ssr-call` on
simulated data with *supplied* parameters** — proves reading→candidates→likelihood→EM→VCF.
*Tests:* on a clean simulated cohort the called genotypes match truth at high-depth loci.

### Milestone D — Parameter pre-pass (estimate the parameters) → **checkpoint 2**

**D1. Confident-genotype resolution test + per-locus fitting** *(prepass.rs + rung_ladder.rs)* —
the 1-vs-2(..p)-peak BIC resolution test (Q-P7, CG-seed): score reads under best 1-peak vs
2-peak model using the **C2 likelihood**, admit at the peak count that earns its parameters;
from each resolved genotype **label the skirts** (hom: one; separated het: two outer skirts
hard-labelled, inner valley → soft EM) and accumulate `SlipProfile` + `SampleStutterStats`.
*Depends:* B1, B3, C2. *Tests:* on simulated confident homs/separated-hets the labelled
skirts match the injected stutter; a 1-apart het contributes nothing.

**D2. Burn-in loop + measure → freeze parameters** *(prepass.rs)* — the adaptive loop:

```rust
// burn-in (params §4); deterministic via seed + barrier, NOT per-locus overwrite.
let mut params = dev_defaults();
let mut l_pen_prev = f64::NEG_INFINITY;
loop {
    let batch = draw_batch(&catalog, seed, BATCH_SIZE /*=32 fixed, Q-P4*/);
    // parallel map: pure fn of (locus, FROZEN params) → per-thread partials
    let partials: Vec<Partial> = batch.par_map(|loc| score_and_accumulate(loc, &params));
    // barrier → reduce in FIXED catalog order (commutative; FixedPointAccum for floats)
    let reduced = reduce_fixed_order(partials);
    params = update_params(&reduced);                 // soft full-cohort EM reduce (C2)
    let l_pen = reduced.l_pen;                          // E-step normalizer, ~free
    if (l_pen - l_pen_prev).abs() / l_pen.abs() < TOL { break; }   // ℓ_pen plateau (M4)
    if iter >= MAX_ITER { flag_never_settled(); break; }
    l_pen_prev = l_pen;
}
// multi-start guard: run from several seeds, keep best ℓ_pen, flag divergent basins (M4).
```

Then **measure** over a stratified loci sample → frozen `ε`, the cohort-per-period shape
parent `θ_period`, and **per-sample** shape/level lines (the clustering input). `ℓ_pen` and
BIC log-liks summed by `FixedPointAccum` (byte-identical across threads). *Depends:* D1.
**▶ Checkpoint 2: recover known parameters on simulated data.** *Tests:* recovered `ε`/
shape/level within tolerance of injected truth; non-monotone `ℓ_pen` fails loudly; m2(a)
all-merged-het cohort → coded-default fallback + loud warning.

**D3. Sample-group clustering + per-(group,period) shape + ε-freeze check** *(sample_groups.rs)* —
distance-based grouping of close neighbours in (`ε`, level) space, **uncertainty-scaled**
(precision ∝ 1/√depth), deterministic tie-break on sample catalog index; group **count**
from a BIC `ℓ_pen(split) − ℓ_pen(merged)` comparison (M4), not preset K. Then, with groups
fixed, **fit the per-`(group,period)` shape shrunk to `θ_period`** (the invariance spread
is the M3 diagnostic) and the **per-group** level line. The `ε`-freeze check: BIC of
frozen-per-group `ε` vs. `ε + covariate`/per-sample (Q-P2, non-circular) — richer model
wins ⇒ warn + point at the upgrade. *Depends:* D2, C4 (the `ε`-check runs calls).
**▶ M3 milestone: on a deliberately group-divergent-shape simulator, recover the per-group
shapes** (and confirm collapse-to-shared when truly invariant). *Tests:* two-protocol
simulated cohort → two groups tracking protocol not depth; single-protocol → one group.

### Milestone E — Genotyping completeness

**E1. Outer loop — `F_i` + per-group stutter level** *(inbreeding.rs)* — the prior-side
loop wrapping the per-locus EM over all loci:

```rust
let (mut f_i, mut level_g) = (vec![param_set.f0_seed; n_samples], param_set.level_seed.clone());
loop {
    let calls = all_loci.map(|loc| run_locus_em(loc, .., &f_i, &level_g, ..));   // align recomputed
    // F_i reduce over VARIABLE loci (≥2 alleles): mean autozygous-branch responsibility
    //   raw → shrink to cohort mean → clamp ≤ user cap → clamp ≤ F_CEILING=0.99
    f_i = reduce_f(&calls);                       // FixedPointAccum (order-independent)
    // level reduce per SAMPLE GROUP from soft per-allele responsibilities → refit baseline+slope·L
    level_g = reduce_level(&calls);               // FixedPointAccum; re-weight S_θ, NO align rebuild
    if converged(Δf, Δlevel) || rounds >= MAX { break; }
}
// genotype calls = the FINAL per-locus E-step
```

*Depends:* C4, D2 (level seed). *Tests:* `F` recovers injected inbreeding on a structured
sim; both reduces byte-identical across thread counts; level refinement overrules a
deliberately-biased seed.

**E2. FP control + full VCF semantics** *(likelihood.rs allele-balance, vcf_out.rs)* —
λ outlier + recurrence + **allele-balance/overdispersion** on deconvolved responsibilities
(qual_refine pattern; feeds GQ + site QUAL → depth-inflated false het → low GQ → no-call);
**emit-iff-variable**; **site QUAL = `−10·log₁₀ P(monomorphic)`** via the generalized
exact-AF kernel (sum-over-alleles, collapse-allele parametrized off allele-0, per-`i` Beta
from `G₀`, REF-anchored wrapper bypassed — verify-fix #7b); SSR FILTER reasons
(`notPeriodic`/`tooManyAlleles`/`lowDepth`); `./.` no-calls; **apparent-`F_IS` warning** in
header+stderr. *Depends:* C4. *Tests:* depth-inflated false het → low GQ (the SNP-path
blind spot, regression-guarded on sim); monomorphic locus dropped; biallelic QUAL exact vs.
brute force; ≥3-allele QUAL within documented approx bound.

### Milestone F — Scale & calibrate

**F1. Parallelism & determinism** — Phase-2 **batched map-reduce-with-barrier** pool +
Phase-3 streaming workers replacing the driver stub; wire the reading layer's parallelism +
two-pass restart (`restart()`, reading Phase 4). **Prove byte-identity across `--threads`**
for both the pre-pass outputs (params, group assignments, `ℓ_pen` stop iteration) and the
VCF. *Depends:* C4, D3. *Tests:* `ssr-call --threads 1..K` → byte-identical VCF + identical
frozen params + identical group memberships.

**F2. Calibration & validation** — fix every Q-P5/Q-G6 number on the simulator
(`dev_defaults`, `MAX_SLIP`, `TOL`, multi-start count, batch size, stratification window,
shrinkage strengths, `G₀` decay, λ, allele-balance strength, rung recurrence `k`,
prominence, `MAX_CANDIDATE_ALLELES`, separation threshold). Then the **correctness gates**
(external anchors, not "calls don't move"): ground-truth recovery, known-protocol positive
control, multi-start `ℓ_pen` agreement, discrete-vs-continuum (`ε`, level) cloud,
**batch-size invariance** (prove, don't assert). *Depends:* all.

---

## 5. Cross-cutting invariants (re-asserted by every threading step)

- **Determinism / byte-identity across thread counts** is the headline invariant. Integer
  sufficient statistics (`SlipProfile`/`SampleStutterStats`) are commutative; **all decision
  floats** (`ℓ_pen`, BIC log-liks, `F`/level responsibilities) go through `FixedPointAccum`
  (scale→round→`i128`); clustering ties break on **sample catalog index**; barrier reduces
  in **fixed catalog order** (never per-locus running overwrite). The burn-in is reproducible
  via `(seed, batch size)` + pure per-batch map.
- **`align` is a pure function of `(obs, cand⊕Δ, ε)`** — recomputed, **no cache in v1**
  (Q-G3); this moots the M2 cache-keying question and removes per-thread cache state.
- **Recall in assembly, precision in the EM** — never silently drop a sequence; a filtered
  locus is emitted with its reason; far candidates kept alive by the `G₀` floor.
- **No silent fallback** — m2(a) (no confident genotype of either kind) and a `ε`-freeze
  loss both **warn loudly** (VCF header + stderr).
- **Shared primitives are written once** — if a peak/kernel/align behaviour differs between
  the two halves, that is a bug, not a feature.

## 6. Open decisions carried in (resolved upstream; calibration-only here)

| Q | status for this plan |
|---|---|
| **Q-P1 / Q-G1** shared `rung_ladder`/`stutter`/`pair_hmm`/`allele_freq_prior`/`em_init` | settled — Milestone B writes them once (§2 layout) |
| **Q-P2** ε-freeze = penalized-likelihood BIC (non-circular) | settled — D3 |
| **Q-P3 / Q-G5** seeds computed in Phase 3, not the pre-pass | settled — C3 |
| **Q-P4** dedicated batched pool, batch size 32 fixed | settled — D2/F1 (size = F2 calibration on big machines) |
| **Q-P6** distance-based grouping (not k-means), uncertainty-scaled, BIC group count | settled — D3 |
| **Q-P7** confident genotype = 1..ploidy-peak BIC resolution test (CG-seed) | settled — D1 (penalty + separation = F2 calibration) |
| **Q-G2** engine extend-in-place vs slim SSR fork of the π-loop | **decided in C4** (struct-shape call; math boundary fixed either way) |
| **Q-G3 / Q-G3b** defer align cache (recompute); `Σ_v` explicit run-collapsed sum | settled — B3/C2 (per-locus memo only if F1 measurement shows align is the wall) |
| **Q-G4** P(monomorphic) = generalized exact-AF kernel | settled — E2 (`lowDepth` threshold + polymorphism flag = F2 calibration) |

## 7. Out of scope (explicit)

- The **reading layer** (Phase 1) — done; this plan consumes its restartable stream and
  only *uses* its `restart()` for the two-pass (F1).
- A persistent `align` **cache** (Q-G3 deferred; per-locus memo is the first lever *if*
  measured necessary), per-sample `ε`, the Gaussian `G₀` form, Dirichlet-multinomial
  overdispersion, per-locus null alleles — all **documented upgrades**, not v1.
- Polyploid dosage resolution beyond the cleanly-resolvable cases (D1 leans on coded
  priors otherwise — a documented approximation).
- Any Stage-1 change (the in-block-locus invariant is *checked*, not modified).

---

## 8. The execution loop (how a group of agents builds this plan)

The plan is executed **one loop iteration per milestone** (A → B → C → D → E → F): a
milestone bundles its steps (B = B1+B2+B3, etc.); the loop drives all of a milestone's steps
to completion, then reviews and hardens the *whole milestone* before the next one starts.
This is the operational form of "pause + review between milestones" (build philosophy) and
of the roadmap's two integration checkpoints.

**Before the first iteration:** ensure the `ssr-call` Stage-2 feature block in
[PROJECT_STATUS.md](../../../PROJECT_STATUS.md) covers the genotyping/pre-pass work (it
currently tracks the reading layer; the feature_implementation skill extends it on first
touch). All three skills read and update that block — that is how lifecycle state
(`implemented` → `reviewed` → `fixes-applied`) is carried between the loop's stages.

### The per-milestone loop

For each milestone in order:

1. **Implement to the milestone boundary** — feature_implementation skill
   ([ia/skills/feature_implementation_skill.md](../../../ia/skills/feature_implementation_skill.md)).
   Drive *every* step of the milestone, **types first then implementation** within each
   step, tests alongside. Honor the skill's interactive-planning gate and its
   **no-silent-assumptions** rule (record every default the spec left open). Output: code +
   tests, an implementation report under `doc/devel/reports/implementations/`, and a
   PROJECT_STATUS.md update to `implemented`. *Exit criteria below must hold before step 2.*
2. **Commit** the implementation. `feat(ssr): <milestone> — <steps>` (e.g.
   `feat(ssr): milestone B — rung ladder, stutter kernel, pair-HMM forward`).
3. **Code review** — code_review skill
   ([ia/skills/code_review_skill.md](../../../ia/skills/code_review_skill.md)), **scoped to
   the milestone's own commits** (this milestone's diff, not the whole `ssr-cohort` branch;
   prior milestones are already reviewed). Output: a review report under
   `doc/devel/reports/reviews/`, PROJECT_STATUS.md → `reviewed`.
4. **Commit** the review report. `docs(ssr): code review of <milestone>`.
5. **Apply fixes** — code_review_fixes skill
   ([ia/skills/code_review_fixes.md](../../../ia/skills/code_review_fixes.md)). Every finding
   gets a terminal status; **Apply** only High-confidence Blocker/Major by default,
   **Ask/Defer** the rest (the skill forbids silently resolving ambiguity). Output: a
   fix-application report, PROJECT_STATUS.md → `fixes-applied` (or `shipped` if nothing
   material remains).
6. **Commit** the fixes + report. `fix(ssr): apply <milestone> review` *(this step is the
   extension of the bare loop: step 5 produces source changes, which must be committed before
   the next milestone builds on them; a re-review is warranted only if a fix was large or
   risky).*

Then advance to the next milestone.

### Two mandatory human gates — do not auto-advance

The loop **halts for human sign-off** (does not roll straight into the next milestone) after
the two integration checkpoints:

- **after C4** — checkpoint 1: first end-to-end `ssr-call` VCF on *supplied* parameters;
- **after D2** — checkpoint 2: the pre-pass recovers the simulator's known parameters.

These are the load-bearing "is the spine / the estimator actually correct" gates. An agent
must **not** declare them passed on its own — the recovery and end-to-end assertions are
exactly the claims a human PM signs off before the next half is built on top.

### Milestone exit criteria (the gate step 1 must clear before step 2)

A milestone is `implemented` only when **all** hold:

- every step's deliverable in §4 exists, and every test listed for those steps passes;
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, and
  `cargo test --all-targets --all-features` are green (the skills' mandated validation);
- the milestone's own checkpoint assertion holds on simulated data where it has one — **C4**
  (called genotypes match truth at high depth), **D2** (recovered ε/shape/level within
  tolerance), **D3/M3** (per-group shapes recovered; collapse-to-shared when invariant);
- every assumption made against an under-specified part of the spec is written down in the
  implementation report (no silent choices buried in code).

### Execution is sequential — one agent at a time, step by step

**Agents work one after another, not in parallel.** Within a milestone, drive the steps in
their dependency order (the §1 chain), finishing and validating one before starting the next
— A1 → A2; B1 → B2 → B3; C1 → C2 → C3 → C4; D1 → D2 → D3; E1 → E2; F1 → F2. Each step is
**types first, then implementation, then its tests**, and the milestone's exit criteria
(above) gate the move from step 1 to step 2 of the loop. There is no fan-out to coordinate
and no cross-agent merge: the next agent picks up the committed state the previous one left.

### Working environment — the dev container (agents have full permissions there)

All build/test/validation runs go through the project's **dev container**, launched by
[scripts/dev.sh](../../../scripts/dev.sh) (podman on Linux, Apple `container` on macOS — see
the dev-container notes). **Agents have full permissions inside it and inside the project
working tree / worktree** — read, write, and edit files, and run `cargo` / git / file
commands, **without asking for permission for routine file or build operations**. Reserve
questions for genuine decisions (an unspecified contract, an ambiguous default, a public-API
or behavior change) — not for reading or writing files in the project. Heed the standing
caution against full image rebuilds (layer onto the existing image; `container build` is slow
and hard to kill) — but that is a *cost* note, not a permission gate.
