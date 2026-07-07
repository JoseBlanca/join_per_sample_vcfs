# Architecture: the shared genotype-EM core

**Status:** draft, 2026-07-07, branch `em-convergence-criterion`. The *how* for
[../specs/unified_genotype_em.md](../specs/unified_genotype_em.md) (the *why*).
Pins the trait surface, the module layout, and the SNP↔SSR mapping so the code
extraction (Phase 2) and the SSR migration (Phase 3) have a fixed target.

**Phase 2 contract:** extract a shared, hook-based EM loop from
`src/var_calling/posterior_engine.rs`, with the SNP caller riding it
**byte-identically**. No SSR work, no behaviour change.

---

## 1. The trait surface

One trait captures the per-record caller-specific behaviour the loop needs. The
SNP implementor (`SnpModel`) is a zero-sized type — all per-record state already
lives in `EmContext` (read-only shape/config) and `RecordScratch` (mutable
buffers), so the model itself carries nothing in Phase 2.

```rust
/// The caller-specific half of the genotype EM. The shared loop
/// (`run_em_loop`) owns iteration, the sufficient statistic, and the
/// convergence rule; the model supplies the E-step, the M-step, and the
/// convergence driver.
pub(crate) trait GenotypeEmModel {
    /// Fill `scratch.posteriors` from the current frequency parameters
    /// and the per-(sample,genotype) read log-likelihoods. `phase`
    /// selects the first-iteration (flat) vs steady-state prior for a
    /// cohort; single-sample ignores it.
    fn e_step<M: MathBackend>(
        &self,
        ctx: EmContext<'_>,
        math: &M,
        log_likelihoods: &[f64],
        scratch: &mut RecordScratch,
        phase: EmStepPhase,
    ) -> Result<(), PosteriorEngineError>;

    /// Advance the frequency parameters from the current posteriors:
    /// accumulate the posterior-weighted allele counts (the shared
    /// sufficient statistic) into `scratch.expected_counts`, then write
    /// the next parameters (`scratch.p_hat_next`, and any auxiliaries).
    fn m_step(&self, ctx: EmContext<'_>, scratch: &mut RecordScratch);

    /// The convergence delta for this iteration: the max change in the
    /// model's *driver* (the quantity the E-step feeds back), on a
    /// dimensionless scale comparable to `convergence_threshold`. SNP:
    /// cohort → `max|Δ(expected_counts/chromosomes)|`; single-sample →
    /// `max|Δp̂|`. This is the Phase-0 discipline, now a per-model hook.
    fn convergence_delta(&self, ctx: EmContext<'_>, scratch: &mut RecordScratch) -> f64;
}
```

Notes:

- **Generic method, not `dyn`.** `e_step<M: MathBackend>` makes the trait not
  object-safe, which is intended: the loop is monomorphised over `(M, Model)`
  exactly like today's `run_em_loop<M>`, so the SNP copy compiles to the current
  machine code (the Phase-1 §Q4 decision).
- **Finalisation stays out of the trait for Phase 2.** `summarise_posteriors`
  (best genotype + GQ) and `compute_qual_via_exact_af` (site QUAL) run *after* the
  loop in `run_em_columnar`; they are not part of the iterate-to-convergence core.
  A `finalize` hook is deferred to Phase 3 (§4), when SSR's different site summary
  (posterior homozygosity) needs it. Keeping them out now minimises the Phase-2
  diff.
- **Contamination stays out too.** The mixture pre-pass rewrites the read
  likelihoods *before* the loop; the loop only ever sees a `log_likelihoods`
  slice. No hook needed.

## 2. What the shared loop owns

`run_em_loop` becomes `run_em_loop<M: MathBackend, Model: GenotypeEmModel>` and
owns exactly:

- the iteration scaffold (iterate to `max_iterations`, the first-iteration
  cohort guard, the `converged` diagnostics);
- calling `model.e_step` → `model.m_step` → `model.convergence_delta` in order;
- the threshold comparison and the `p_hat`/`p_hat_next` buffer swap.

The sufficient statistic (`accumulate_expected_counts`) is called from within
`m_step`, so it is shared by being the one implementation every model's M-step
uses. The convergence *rule* (compare the model's driver delta to the threshold,
never converge on the cohort flat step) is owned by the loop; the *driver* is the
model's to define.

## 3. Module layout

**Phase 2 (SNP-only): the core lives under `var_calling`.**
`src/var_calling/genotype_em/` — a new submodule holding the `GenotypeEmModel`
trait and the generic `run_em_loop`. `SnpModel` and all SNP specifics
(`e_step_*`, `m_step_*`, pseudocount classification, contamination, QUAL) stay in
`posterior_engine.rs`.

**Why not crate-level yet.** The loop reads and writes `RecordScratch` and
`EmContext`, which today carry SNP-and-QUAL-shaped fields (`f_hat_compound`, the
mixture buffers, the exact-AF buffers). Hoisting those to a crate-level "shared"
module in Phase 2 would drag SNP specifics into the shared home *before* SSR
exists to justify the split — a large, byte-identity-risking type migration for
no Phase-2 benefit. So Phase 2 keeps the core under `var_calling` (only SNP
consumes it), and **Phase 3 hoists it to crate-level `src/genotype_em/`** while
splitting `RecordScratch` into a shared core-state struct + a per-model
extension, at the moment SSR forces the boundary to be real. *(This refines the
Phase-1 "crate-level" lean: same destination, but reached in Phase 3 to avoid a
premature migration. Flagged rather than silently changed.)*

## 4. SNP ↔ SSR mapping (Phase 3 preview, not built yet)

| Trait method | `SnpModel` (Phase 2) | `SsrModel` (Phase 3) |
|---|---|---|
| `e_step` | `dispatch_e_step` (flat / LOO / SIMD / single) | HWE(π)→marginalized DM prior over its mode-centred seed + stutter read-lik |
| `m_step` | `m_step_p_hat` + `m_step_f_hat_compound` | `α = G₀ + LOO(expected_counts)` |
| `convergence_delta` | cohort `expected_counts`/chromosomes; single `p̂` | same driver (`expected_counts`) |
| finalize (Phase 3) | best GT + GQ + exact-AF QUAL | best GT + GQ + posterior homozygosity |

The SSR stutter-refit outer loop wraps the whole thing on the SSR side (Phase-1
§Q2); the shared loop is unaware of it.

---

## 5. Phase-2 step plan

Each step is byte-identical for SNP and ends with the full suite green plus the
tomato1 cohort VCF diff at zero against the pre-Phase-2 binary. See
[../implementation_plans/unified_genotype_em.md](../implementation_plans/unified_genotype_em.md).

1. **2.1** — this doc + the plan (design; no code). *Done.*
2. **2.2** — introduce `GenotypeEmModel` + `SnpModel`; route `m_step` +
   `convergence_delta` through the trait. *Done, byte-identical.*
3. **2.3** — move the E-step behind the trait (`SnpModel::e_step` delegates to the
   existing `dispatch_e_step`). *Done, byte-identical. Phase 2 ends here.*
4. **2.4 / 2.5 — deferred to Phase 3.** The module relocation lands at its real
   crate-level home (`src/genotype_em/`, §3) together with the `RecordScratch`
   split, rather than an intermediate move into `var_calling/genotype_em/` that
   would be redone. The trait boundary + generic loop (2.2–2.3) already deliver
   Phase 2's value; Phase 3 does the physical move once.

## 6. Open (resolved in §7 for Phase 3)

- Splitting `RecordScratch` into shared-core vs per-model state, and the
  crate-level hoist to `src/genotype_em/`.
- The `finalize` hook (SSR posterior-homozygosity summary vs SNP exact-AF QUAL).
- Where the SSR mode-centred seed (`G₀`) and stutter read-model attach as model
  state (the `SsrModel` fields).

---

## 7. Phase 3 — generalize the trait to host SSR, and share the DM prior

Phase 3 does two things the user asked for: **(a)** SSR adopts the improved
(marginalized Dirichlet-multinomial) prior, and **(b)** SSR rides the *sensible*
shared parts of the loop. It does **not** merge SSR's outer stutter-refit loop or
its stutter read-model — those stay SSR-side (Phase-1 §Q2). It respects the SSR
`Q-G2` decision ("a slim SSR EM, not a `posterior_engine` graft") by sharing a
thin loop skeleton + a prior *function*, not the SNP engine's baggage (chain
anchors, class pseudocounts, contamination, exact-AF QUAL all stay `SnpModel`-side).

### 7.1 The prior machinery is already shared

`crate::genetics::dirichlet_multinomial_log_priors(genotype_allele_counts,
log_multinomial_coeffs, n_alleles, alpha)` is a crate-level primitive taking flat
arrays — no back-reference into either caller. It **is** "the machinery." What
differs per caller is only the **concentration seed** `α`:

- **SNP:** `alpha_from_diversity(n_alleles, θ̂)` = REF-privileged neutral SFS.
- **SSR:** the mode-centred `g0_pseudocounts(...)` vector (hypervariable, no
  privileged reference — Phase-1 §Q1 guardrail).

Both then add leave-one-out cohort counts (`α'_s = α_seed + (E[cohort] − E[own_s])`)
and the Wright-`F` IBD mixture on top. SSR's current `genotype_prior` (plug-in
HWE(π), no LOO) is replaced by this.

### 7.2 Lightweight sharing — decided 2026-07-07

**Share the meaningful bits, keep two slim loops. The SNP E-step is not
touched.** The heavyweight alternative (one generic `run_em` both callers ride,
via an associated-`Scratch` trait + model-carried inputs) was rejected: it forces
a byte-identity/perf-risky SNP refactor and wraps SSR's deliberately-slim loop in
trait machinery — the very thing the SSR `Q-G2` note pushed back on, and a threat
to the SNP SIMD gains (Phase-1 §Q4). "The parts that make sense to share" are the
prior and the convergence *discipline*, not a five-line `for` loop.

Concretely:

- **Prior machinery — already shared.** `dirichlet_multinomial_log_priors`
  (§7.1). SSR gets a *new, SSR-side* prior wrapper that builds `α'_s = G₀ +
  (E[cohort] − E[own_s])` (marginalize **+ leave-one-out**, decided) and applies
  the Wright-`F` IBD mixture, then calls that shared primitive. This wrapper is
  written in SSR code (a scalar per-locus path — SSR has no SIMD to protect),
  **mirroring** SNP's `e_step_cohort_loo` rather than extracting from it, so the
  SNP E-step (scalar *and* SIMD) is left byte-for-byte alone. (If a genuinely
  perf-neutral, byte-identical extraction of the scalar LOO+`F` fragment turns out
  clean, it can be shared later — but only gated on the tomato1 timing *and* diff,
  never at the cost of the SIMD path.)
- **Convergence discipline — a tiny shared helper.** The Phase-0 rule "test the
  driver (`expected_counts`), not a rescaled readout" becomes a small shared
  function both loops call. It matters for the new SSR prior: with the marginalize
  +LOO DM, SSR's E-step feeds back `expected_counts` (via the LOO `α`), so the new
  path must converge on `expected_counts`, not on the reported `π` — exactly the
  Phase-0 lesson, now applied to SSR. (SSR's *current* plug-in already tests its
  own driver `π` correctly, so this only bites the new prior path.)
- **Loops stay separate.** SNP keeps `run_em_loop`; SSR keeps `run_pi_em`. The
  Phase-2 `GenotypeEmModel` trait stays a tidy SNP-internal abstraction (not
  stretched over SSR).

### 7.3 What is shared vs per-caller (lightweight)

| Shared | Per-caller |
|---|---|
| `genetics::dirichlet_multinomial_log_priors` (already) | the E/M loop (SNP `run_em_loop`; SSR `run_pi_em`) |
| a small convergence-discipline helper (test the driver) | the frequency **seed** (SNP θ̂; SSR `G₀`) |
| | the LOO+`F` prior wrapper (SNP has one; SSR gets its own, mirrored) |
| | read-likelihood model, finalize, SSR outer stutter-refit loop |

The SNP E-step (SIMD + fast paths) is **untouched** — the safest guarantee for the
Phase-1 §Q4 perf gains.

### 7.4 Phase-3 steps (lightweight)

1. **3.1** — this design + plan (no code). *Done.*
2. **3.2** — build the SSR marginalize+LOO+`F` DM prior wrapper (mode-centred `G₀`
   seed, via `dirichlet_multinomial_log_priors`), as additive SSR code with unit
   tests: the high-concentration limit ≈ plug-in HWE, hand-computed small cases,
   and the LOO exclusion. **SNP untouched.**
3. **3.3** — extract the small shared convergence-discipline helper; SSR's *new*
   prior path converges on `expected_counts`. SNP `run_em_loop` adopts the same
   helper only if byte-identical (else left as-is).
4. **3.4** — wire the new prior into `run_pi_em` behind an `EmCfg` toggle
   (default = current plug-in). SSR byte-identical at default (ssr_tomato1
   unchanged); new prior opt-in.
5. **3.5** — benchmark ssr_tomato1 vs HipSTR (plug-in vs marginalized). Watch
   specifically for the `G₀`-as-concentration effect (a too-diffuse or too-tight
   prior), not just overall concordance. Flip the default only if it
   holds/improves; else keep opt-in and record the result.

3.2 is additive (SNP untouched, SSR unchanged until the toggle). 3.4 is
toggle-gated (SSR byte-identical at default). Only 3.5 decides the SSR call
change, on the benchmark.
