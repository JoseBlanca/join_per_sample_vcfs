# Code Review: ssr-call θ_locus shape refit (Step I1)

**Date:** 2026-06-24
**Reviewer:** rust-code-review skill (orchestrator; inline synthesis)
**Scope:** commit `46df902` on branch `ssr-cohort` — the per-locus `θ_locus` shape M-step in the genotype EM
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** PR diff, commit `46df902` (Step I1).
- **In-scope files:** [em.rs](../../../../src/ssr/cohort/em.rs) — the `θ_locus` loop in `run_locus_em_with`; the extracted `compute_data_ll` / `run_pi_em` / `final_calls`; new `attribute_locus_slips` / `add_slip` / `shapes_close`; `EmCfg` fields; 4 new tests.
- **Out of scope:** `refine_theta_locus` / `SlipProfile` (reviewed in B2); the burn-in / sweep callers (they call `run_locus_em_with` unchanged).
- **Categories considered:** reliability, errors, naming, idiomatic, defaults, smells, refactor_safety, extras (determinism, hot path).

## 2. Verdict

**Approve-with-changes.** The θ-loop is correct, the helper extraction is behaviour-preserving (the `theta_max_rounds = 0` path reproduces the pre-I1 EM), shrinkage gives the anti-oscillation guarantee, and determinism / byte-identity hold (integer slip counts; the existing cross-thread test still passes). No Blocker/Major. Minors: `add_slip` duplicates `prepass::add_slip`, and a per-locus compute-cost increase worth noting.

## 3. Execution status

- `cargo fmt --check` clean · `cargo clippy --lib -D warnings` clean · `cargo test --lib` = **1275 passed, 0 failed, 2 ignored** (+5 I1 tests; all prior EM/byte-identity tests unchanged).
- Sub-agent fan-out not used (repeated overload); reviewed inline. "Needs verification": 0.

## 4. Open questions and assumptions

1. **Hard-label attribution at the valley.** `attribute_locus_slips` assigns each read to its nearest called allele via `min_by_key`, which on an exact tie picks the first in `genotype_units` (candidate-index order). Deterministic (so byte-identity holds) and consistent with `inbreeding.rs::reduce_level`, but it is the documented hard-label approximation; the soft per-read responsibility split is the deferred refinement. No action — flagged as the known approximation.

## 5. Top 3 priorities
1. **Mi1** — lift one shared `add_slip` (DRY with `prepass::add_slip`).
2. **Mi2** — note (and optionally bound) the per-locus cost: data-ll is recomputed and π restarted from `pi0` each θ round.
3. Pin the `theta_max_rounds = 0` ⇒ pre-I1-behaviour invariant with a test.

## 6. Findings

### Minor

**Mi1: [em.rs](../../../../src/ssr/cohort/em.rs) `add_slip` duplicates [prepass.rs](../../../../src/ssr/cohort/prepass.rs) `add_slip`.**
The same `MAX_SLIP`-bounded signed-delta binning now exists in two private copies. Confidence: High.
*Fix:* lift one `pub(crate) fn add_slip(&mut SlipProfile, delta: i32, count: u64)` next to `SlipProfile` in `param_estimation.rs` and call it from both. Small, clearly-correct DRY.

**Mi2: [em.rs](../../../../src/ssr/cohort/em.rs) — per-locus cost grows with the θ loop.**
Each θ round recomputes the full `data_ll` (the `obs × candidate × genotype` read-likelihood matrix) and `run_pi_em` restarts π from `seed.pi0` rather than warm-starting from the previous round's π. Bounded by `theta_max_rounds = 3` (and early `shapes_close` break — clean loci converge in one round), so worst case ≈ 2–4× the pre-I1 per-locus EM. Confidence: High (by construction).
*Fix:* none required for correctness; note it, and consider warm-starting π across θ rounds as a later optimization (it would not change the converged result materially but would cut iterations). The cost is on the (currently single-threaded) sweep — relevant to Milestone J's throughput.

### Nits
- `compute_data_ll` takes 10 args (`#[allow(clippy::too_many_arguments)]`); acceptable for an internal extraction, but a small `LocusEmContext` borrow-struct would read better if it grows.
- The initial `data_ll`/`pi`/`calls` are always computed before the loop, so `theta_max_rounds = 0` exactly reproduces the pre-I1 EM — worth a one-line test to pin it (section 8).

## 7. Out of scope observations
None.

## 8. Missing tests to add now
Under `em.rs`:
- `theta_max_rounds_zero_reproduces_the_seed_shape_result` — run the moderate-stutter cohort with `EmCfg { theta_max_rounds: 0, .. }` and with the default; assert both call the truth (pins that the refit is purely additive and the disabled path is the old behaviour). A strict "refit ≥ no-refit GQ" assertion is fragile on clean data — skip.

## 9. What's good
- The extraction is genuinely behaviour-preserving: `compute_data_ll` is the verbatim former inline block parametrised on `theta`, so the θ-loop is a clean wrapper and the `theta_max_rounds = 0` path is provably the old EM. [em.rs](../../../../src/ssr/cohort/em.rs)
- Shrinkage toward `seed.theta0` plus the `shapes_close` stop gives the anti-oscillation property *for free* — a thin locus's sparse profile collapses to the seed (`refine_theta_locus` returns the prior on an empty profile), so no special-casing of low-depth loci is needed.
- Integer slip counts keep the M-step order-independent, so the per-locus refit does not threaten the cross-thread byte-identity the stage guarantees — and the existing determinism test confirms it.

## 10. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --lib --all-features -- -D warnings` · `cargo test --lib ssr::cohort`

### Author response convention
Address Mi1, Mi2, Nits by id; answer open question 1 first.
