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

### 7.2 The generalized trait

The Phase-2 trait took concrete `RecordScratch` / `EmContext` / `math` / `phase`.
To host SSR, generalize it: an **associated `Scratch` type**, and the **model
carries its own inputs** (math backend, read likelihoods, context) so the loop is
fully generic and touches nothing caller-shaped.

```rust
pub(crate) trait GenotypeEmModel {
    /// Per-record/per-locus mutable EM buffers (SNP: the core of today's
    /// RecordScratch; SSR: its own π / expected-count buffers).
    type Scratch;

    /// Fill the posteriors from the current frequency parameters + read
    /// likelihoods. `iteration` lets the model pick its own prior variant
    /// (SNP cohort: flat on iter 1, leave-one-out after).
    fn e_step(&self, iteration: u32, scratch: &mut Self::Scratch) -> Result<(), EmError>;

    /// Advance the frequency parameters from the posteriors (accumulate the
    /// shared sufficient statistic, write the next params, rotate buffers).
    fn m_step(&self, scratch: &mut Self::Scratch);

    /// Max change in the model's driver — the Phase-0 discipline. Called once
    /// per iteration, after `m_step`.
    fn convergence_delta(&self, iteration: u32, scratch: &mut Self::Scratch) -> f64;

    /// Whether this iteration is allowed to declare convergence. SNP cohort
    /// overrides to `false` on iteration 1 (the flat-prior seed); default `true`.
    fn allow_convergence(&self, _iteration: u32) -> bool { true }
}
```

The shared loop is then genuinely generic and owns only the iteration scaffold +
the convergence rule:

```rust
pub(crate) fn run_em<Model: GenotypeEmModel>(
    model: &Model, scratch: &mut Model::Scratch, max_iters: u32, tol: f64,
) -> Result<EmOutcome, EmError> {
    let mut delta = f64::INFINITY;
    for iteration in 1..=max_iters {
        model.e_step(iteration, scratch)?;
        model.m_step(scratch);
        delta = model.convergence_delta(iteration, scratch);
        if delta < tol && model.allow_convergence(iteration) {
            return Ok(EmOutcome::converged(iteration, delta));
        }
    }
    Ok(EmOutcome::capped(max_iters, delta))
}
```

`SnpModel<'a, M>` carries `math`, `log_likelihoods`, `EmContext`, and `config` as
borrows; its `Scratch` is the SNP EM buffers. The `p_hat`/`p_hat_next` swap moves
*into* `SnpModel::m_step` (model-internal), so the loop no longer names SNP
buffers. SNP stays **byte-identical** (same calls, reorganised behind the trait).

### 7.3 What is shared vs per-model ("the parts that make sense")

| Shared (in `src/genotype_em/`) | Per-model (caller-side) |
|---|---|
| `run_em` loop skeleton + convergence rule | E-step / M-step bodies |
| `EmOutcome` (converged? iters, delta) | `Scratch` type + context |
| `genetics::dirichlet_multinomial_log_priors` (already shared) | frequency **seed** (SNP θ̂ vs SSR `G₀`) |
| the LOO α-update + Wright-`F` mixture *pattern* (extract a helper both call) | read-likelihood model (SNP precomputed; SSR stutter) |
| | finalize (SNP exact-AF QUAL; SSR posterior-hom) |
| | SSR outer stutter-refit loop (wraps `run_em`) |

### 7.4 Phase-3 steps

1. **3.1** — this design + plan update (no code).
2. **3.2** — generalize the trait (associated `Scratch`, model-carried inputs,
   `allow_convergence`) + generic `run_em`; `SnpModel` rides it **byte-identical**.
3. **3.3** — hoist the shared trait + `run_em` + `EmOutcome` to crate-level
   `src/genotype_em/`; split the SNP `RecordScratch` so its EM-core fields are
   reachable. SNP byte-identical.
4. **3.4** — extract the shared **LOO α-update + Wright-`F` mixture** helper from
   `SnpModel`'s E-step (SNP byte-identical), so SSR can reuse it.
5. **3.5** — build the SSR marginalized-DM prior (G₀ seed + LOO + Wright-`F` via
   the shared helper + `dirichlet_multinomial_log_priors`), unit-tested against
   the plug-in in the high-concentration limit and against hand-computed values.
6. **3.6** — `SsrModel` implements `GenotypeEmModel`; wire `run_pi_em` onto
   `run_em` with the new prior **behind a config toggle (default = current
   plug-in)**, so SSR output stays byte-identical by default.
7. **3.7** — benchmark ssr_tomato1 vs HipSTR (plug-in vs marginalized); flip the
   default only if concordance holds/improves. Else keep opt-in.

Steps 3.2–3.4 are SNP-byte-identical refactors (same guardrail as Phase 2).
3.5–3.6 are additive + toggle-gated (SSR byte-identical at default). Only 3.7
decides the SSR call change, on the benchmark.
