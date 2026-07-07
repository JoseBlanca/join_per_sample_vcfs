# One EM machine for SNP and SSR genotyping

**Status:** design draft, 2026-07-07, branch `em-convergence-criterion`. A design
study — *why* the SNP and SSR genotypers can share one expectation–maximisation
(EM) core, *what* that shared core and its per-domain plug-in points look like,
and *how* to get there in reviewable steps without disturbing either
freshly-validated caller. **No code yet.**

This is the umbrella for the arc whose Phase 0 already landed:
[em_convergence_criterion.md](em_convergence_criterion.md) fixed the SNP cohort
EM to test convergence on the quantity the iteration actually drives. That fix is
the first brick here — it put the SNP engine on the same principle the unified
core is built around.

**Scope.** The genotype-calling EM only: the per-record/per-locus loop that turns
read likelihoods + a genotype prior into posteriors, allele frequencies, and
per-sample calls. It does **not** touch pileup, `.psp` I/O, candidate discovery,
the paralog filter, or (except as a plug-in) the SSR stutter read-model. It is
deliberately staged so the code-sharing steps are **byte-identical** for both
callers; the one step that changes SSR science is isolated and optional (§6, §7).

---

## 0. Vocabulary

- **EM (expectation–maximisation).** Alternate an **E-step** (posterior over
  genotypes given current parameters) with an **M-step** (re-estimate parameters
  from the posteriors) until the parameters stop moving.
- **Genotype prior.** How probable each genotype is *before* reading this
  site's data. Both callers build it from an allele-frequency model + a
  Wright inbreeding-`F` term.
- **Dirichlet–multinomial (DM).** The distribution you get when allele
  frequencies are themselves uncertain (drawn from a Dirichlet with concentration
  vector `α`) and genotypes are then multinomial draws. Its **concentration**
  `Σα` controls how *sharp* the frequency belief is: large `Σα` = "I know the
  frequency, it is `α/Σα`"; small `Σα` = "the frequency is very uncertain."
- **Plug-in vs marginalize.** *Plug-in*: estimate one frequency `π̂` and apply
  Hardy–Weinberg as if `π̂` were fact. *Marginalize*: average the genotype
  probability over the whole frequency belief. These are the **two ends of the DM
  concentration axis** (§2) — the single fact that makes unification real.
- **Sufficient statistic.** The one summary of the posteriors the next iteration
  actually needs. For both callers it is the **posterior-weighted allele counts**
  (`expected_counts` in the SNP engine).

---

## 1. The goal

Today there are two genotype-calling EMs — `run_em_loop` in
`src/var_calling/posterior_engine.rs` (SNP/indel) and `run_pi_em` /
`run_locus_em_with` in `src/ssr/cohort/em.rs` (SSR). They were written
independently, carry different pseudocount structure, and drifted (the Phase-0
bug existed in one and not the other for a *structural* reason). The goal is a
**single EM core** both callers ride, so:

- there is one place where the E/M/convergence logic lives and is tested;
- the "test the driver, not the readout" discipline is enforced by construction
  for any future caller;
- domain differences (allele space, read model, frequency prior) are explicit
  **plug-in points**, not divergent copies of the loop.

---

## 2. The unifying theory — one DM family, one concentration knob

The two machines look different but are the **same statistical object** at
different settings of one knob.

- **SSR today = plug-in HWE(π).** `run_pi_em` estimates a point frequency `π`,
  applies Hardy–Weinberg at `π`, and feeds `π` back into the next E-step. In DM
  terms this is the **large-concentration limit** (`Σα → ∞` with `α/Σα = π̂`): the
  frequency belief is a spike at `π̂`, so DM collapses to plain multinomial =
  HWE(π̂).
- **SNP today = marginalized DM.** `run_em_loop` never commits to a point
  frequency. Its per-sample prior is DM with concentration
  `α'_s = α_species(θ̂) + leave-one-out cohort counts` — a **finite** concentration
  (small in small cohorts), i.e. a genuinely spread frequency belief that is
  *averaged over*, not plugged in. (This is exactly the plug-in→marginalize move
  the [sfs_genotype_prior.md](sfs_genotype_prior.md) work made to fix the
  low-coverage `1/1→0/1` error.)

So **both priors are DM(`α`)**; they differ only in how `α` is built and how
sharp it is. That is the backbone of the whole design: unification is not bolting
two callers together — it is recognising they are one family and exposing the
concentration/`α`-builder as the plug-in point. It also fixes the *direction* of
unification: **toward the marginalized core** (finite `α`), because that is the
statistically-correct end that the SNP work already validated. Copying SSR's
plug-in limit into the SNP path would re-introduce the low-coverage bug — a
non-starter (see em_convergence_criterion.md §Scope).

---

## 3. The two machines today, side by side

The structural facts the shared skeleton has to accommodate:

| Aspect | SNP engine (`posterior_engine.rs`) | SSR engine (`ssr/cohort/em.rs`) |
|---|---|---|
| Loop entry | `run_em_loop` | `run_pi_em`, wrapped by `run_locus_em_with` |
| Genotype space | general `(ploidy, n_alleles)` via `GenotypeShape` | diploid-only (v1), `enumerate_diploid_genotypes` |
| Prior form | DM(`α`), Wright-`F` mixture, cohort leave-one-out | plug-in HWE(`π`), Wright-`F` (`genotype_prior`) |
| Frequency `α`/`π` builder | `α_species(θ̂)` + LOO counts | `G₀` (mode-centred) + counts, normalised to `π` |
| Feedback variable (driver) | `expected_counts` (Phase 0) | `π` |
| Read likelihood → genotype LL | precomputed input; optional contamination mixture pre-pass | `compute_data_ll` via `ReadLikelihoodModel` (stutter) |
| Convergence test | `max|Δ(expected_counts/chromosomes)| < 1e-3` | `max|Δπ| < 1e-6` |
| Outer refinement | none | refit `θ_locus` (stutter shape) + stutter-rate multiplier, `refit_max_rounds` |
| Final outputs | best GT, GQ, site QUAL via exact-AF convolution | best GT, GQ, posterior homozygosity |
| Perf specialisation | biallelic-diploid fast path, SIMD E-step | none (locus-parallel outside the loop) |

The important observation: the **driver is the same sufficient statistic in
both** — posterior-weighted allele counts. SSR's `π` is just those counts
normalised (`π = (G₀ + counts)/total`); the SNP engine keeps them unnormalised
and adds `α_species`. Same information, two conventions. Phase 0 already made both
converge on *that* quantity.

---

## 4. The shared skeleton and its plug-in points

The generic core is small. Everything caller-specific becomes a hook:

```
loop:
    posteriors        = e_step(prior_params, read_lik)     # hook A + hook C
    suff_stats        = accumulate_expected_counts(posteriors)   # SHARED
    prior_params_next = m_step(suff_stats)                 # hook B
    delta             = max_abs_diff(driver(prior_params),
                                     driver(prior_params_next))   # SHARED (Phase 0 rule)
    if delta < tol: break
finalise(posteriors) -> calls, GQ, site-level summary       # hook D
```

The hooks (the interface each caller implements):

- **Hook A — prior builder.** From `prior_params` (+ static `α_species`/`G₀`,
  per-sample `F`) produce the per-genotype log-prior. SNP: DM(`α`) with LOO +
  Wright-`F`. SSR: HWE(`π`) with Wright-`F`. *Both are DM(`α`)* (§2), so this hook
  can be a single DM evaluator parameterised by an `α`-builder.
- **Hook B — parameter update (M-step).** From `suff_stats` produce the next
  `prior_params`. SNP: `α = α_species + LOO(suff_stats)`. SSR:
  `π = (G₀ + suff_stats)/total`.
- **Hook C — read-likelihood model.** Per-`(sample, genotype)` log-likelihood.
  SNP: precomputed (with the contamination mixture as a pre-pass). SSR: the
  `ReadLikelihoodModel` stutter model.
- **Hook D — finalisation.** Best genotype + GQ are shared; the site-level
  summary differs (SNP exact-AF QUAL vs SSR posterior homozygosity).

Two things stay **shared and non-negotiable**: `accumulate_expected_counts`
(the sufficient statistic) and the **convergence rule** (test the driver — the
Phase-0 discipline, now a property of the core rather than of one caller).

The SSR **outer refit loop** (`θ_locus` + stutter-rate) wraps the whole skeleton:
it is an outer driver that re-parameterises Hook C between full EM runs. It stays
SSR-only — modelled as "run the core to convergence, refit read-model params,
repeat" — not pushed into the shared loop.

---

## 5. The hard seams (honest design tensions)

1. **Different feedback topology.** SNP feeds back `expected_counts`; SSR feeds
   back `π`. Resolved by §3: same sufficient statistic, different normalisation.
   The core carries `expected_counts`; each M-step hook maps it to its own
   `prior_params`. No real conflict.
2. **Prior *form* looks different** (DM-LOO vs plug-in HWE). Resolved by §2: both
   are DM(`α`); plug-in is the `Σα→∞` limit. Whether SSR keeps the plug-in limit
   or moves to finite `α` is the §6 (b) decision, not a blocker for sharing code.
3. **Genotype generality.** SNP is general `(ploidy, n_alleles)`; SSR is
   diploid-only v1. The shared genotype-shape machinery (`GenotypeShape`,
   `accumulate_row_counts`) already generalises; SSR adopts it and gets a
   polyploid path "for free" when it needs one.
4. **Read-model asymmetry.** SNP likelihoods are precomputed upstream; SSR
   computes them from the stutter model inside the loop. Hook C absorbs this —
   the core asks for `LL[sample][genotype]`; how it is produced is the hook's
   business.
5. **SNP-only extras** — contamination mixture pre-pass, exact-AF site QUAL,
   SIMD/biallelic fast paths. These stay SNP-side (pre-pass before the core, and
   Hook D), not forced onto SSR.
6. **Byte-identity is the safety rope.** Phases 2–3 must reproduce each caller's
   current output exactly (modulo the Phase-0 change already merged). That is the
   test oracle: extract the skeleton, prove both callers still emit identical
   VCFs, *then* consider the model change.

---

## 6. What is shared vs swappable (decided)

The decision (2026-07-07): **share the SNP path's prior *machinery*; keep the
frequency *seed* swappable.** The SNP genotype prior is two separable layers, and
they go to different places:

- **Shared machinery** — marginalize-over-frequency + leave-one-out cohort
  borrowing + Wright-`F` mixture. This is the nuanced, validated part (it fixed
  the low-coverage `1/1→0/1` error). SSR adopts all of it.
- **Swappable seed** — the concentration recipe (what allele frequencies to
  expect before data). **SNP** seeds "the reference allele is common, alts are
  rare" (neutral SFS from θ̂). **SSR must NOT inherit that** — microsatellites are
  hypervariable and the reference-length allele is often not the common one, so
  SSR keeps its **mode-centred `G₀`** seed. Read model and final-summary are
  swappable too.

The SNP engine already separates these layers (the seed is computed once and
handed to the machinery), so this boundary matches the existing code.

**Consequences for byte-identity:**

- **SNP** rides machinery it already runs → **byte-identical** (pure refactor).
- **SSR** moving from plug-in to the marginalized machinery **changes SSR calls**.
  This is accepted (the change is toward *more correct*), so the old "code-only,
  byte-identical SSR (Phase 3)" and "optional model change (Phase 4)" **merge into
  one benchmark-gated step**. Guardrails: (i) validate on ssr_tomato1
  concordance-vs-HipSTR — the change should help or be neutral at thin-coverage
  loci and vanish where coverage is high; (ii) keep the seed mode-centred (do not
  regress toward the reference allele).

---

## 7. Phasing

- **Phase 0 — SNP convergence fix.** *Done* (`a36f7dd`). Put the SNP engine on the
  driver-based convergence rule; benchmark-validated.
- **Phase 1 — this design doc.** The DM-concentration framing, the shared
  skeleton, the hooks, the seams. No code.
- **Phase 2 — extract the shared skeleton** from the SNP engine (generic over
  Hooks A–D), SNP riding it **byte-identically**. A companion architecture doc
  (`doc/devel/architecture/unified_genotype_em.md`) pins the trait/struct shapes
  before code.
- **Phase 3 — migrate SSR onto the shared marginalized machinery** (§6). SSR keeps
  its **mode-centred seed**, stutter read-model, and outer refit loop as swappable
  parts, but adopts the shared marginalize + LOO + `F` machinery. This **changes
  SSR calls** and is gated on ssr_tomato1 concordance-vs-HipSTR. (Merges the old
  code-only + model-change phases per the §6 decision.)

**Decided:**

- Unify toward the **marginalized DM core**, not the plug-in limit.
- Shared and fixed: the sufficient statistic (`expected_counts`), the
  driver-based convergence rule, and the marginalize + LOO + `F` **machinery**.
- Swappable per caller: the frequency **seed** (SSR stays mode-centred, does not
  inherit the SNP reference-is-common seed), the read model, and finalisation.
- SNP refactor is byte-identical; the SSR change is accepted and benchmark-gated.
- **The SSR stutter-refit outer loop stays SSR-side, wrapping the shared engine as
  a black box** (the engine genotypes once to convergence and returns the calls;
  re-preparing read likelihoods across rounds is the caller's concern). The engine
  is *not* generalised to know about read-model refitting. The outer loop already
  tests convergence on the stutter parameters it feeds back, so it needs no
  Phase-0-style change.
- **SSR adopts the shared general (ploidy-aware) genotype bookkeeping** and retires
  its diploid-only enumerator; genotype enumeration is a *shared* part, not a
  per-caller swappable one. This already covers SSR's diploid-multiallelic case and
  removes ploidy as a blocker — a deliberate step *toward* future polyploid SSR
  (which still needs polyploid work in the SSR stutter read-model and finalisation;
  the bookkeeping just no longer stands in the way).
- **Share via compile-time generics (monomorphisation), not dynamic dispatch** —
  the same technique the `MathBackend` already uses, so the SNP copy compiles to
  today's machine code. Perf specialisations stay **shape-driven** (the biallelic
  `(2,3)` fast path fires on shape, so SNP biallelic sites still hit it; SSR's
  multi-allelic loci take the general path). Guardrail: the tomato1 cohort timing
  (~66 s) must not regress in Phase 2; fallback if it does is to keep the one hot
  step SNP-specific. The EM is only ~3% of cohort wall time, so the risk is small.

### SIMD outlook for SSR (decided)

Two layers, opposite SIMD-friendliness:

- **Genotyping E-step — the free win, taken.** SSR presents read likelihoods as the
  same `f64` [sample × genotype] table as SNP, so it **inherits the 4-samples-at-once
  SIMD E-step for free** by riding the shared engine (general path, not a hardcoded
  fixed-shape kernel). No SSR-specific work.
- **Stutter read-model (`q_r`) — deferred.** Computing the stutter likelihood is a
  sequence-alignment-style dynamic program over variable-length reads — the classic
  *irregular* case that resists simple lane-parallel SIMD. Vectorising it is a
  specialised, SSR-only effort, and `q_r` is largely a *pileup*-stage cost (~75% of
  pileup self-time, re-evaluated per stutter-refit round in Stage 2). **Deferred**:
  pursue only if a profile shows the read model dominates SSR's genotyping runtime.

**Open (for the Phase-2 architecture doc):**

- The exact trait/struct surface for Hooks A–D, and where the SNP contamination
  pre-pass and exact-AF QUAL sit relative to the core.
