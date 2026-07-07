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

## 6. Open (for Phase 3)

- Splitting `RecordScratch` into shared-core vs per-model state, and the
  crate-level hoist to `src/genotype_em/`.
- The `finalize` hook (SSR posterior-homozygosity summary vs SNP exact-AF QUAL).
- Where the SSR mode-centred seed (`G₀`) and stutter read-model attach as model
  state (the `SsrModel` fields).
