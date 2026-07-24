# ng alignment — marginal aligners, Milestone A (`LogProb` + the `MarginalAligner` trait)

**Date:** 2026-07-24
**Plan:** [alignment_marginal.md](../../ng/impl_plan/alignment_marginal.md) (plan 2 of 3), steps A1–A2, at **Checkpoint A**.
**Design authority:** spec [alignment.md](../../ng/spec/alignment.md) §5, §7; arch [alignment.md](../../ng/arch/alignment.md) §1, §3.
**Method:** plan-driven implementation — one implement → review → apply-fixes → commit loop per step. Review = per-category sub-agents on the step's working-tree diff; per-category outputs left in the gitignored `tmp/review_2026-07-24_ng-marginal-a{1,2}/` as the audit trail; this report consolidates the milestone.

---

## 1. Plan

Milestone A lands the vocabulary and the interface for the marginal aligners, with **no algorithm** (types-first, project rule):

- **A1** — `LogProb` in `types.rs`.
- **A2** — the `MarginalAligner` trait, no implementations.

Both were built and committed one loop at a time. Algorithms 5 and 6 are Milestones B and C.

## 2. Assumptions / reconciliations with what plan 1 shipped

The arch's signatures are explicitly illustrative; plan 1 (best-path) diverged from them during implementation, and A2 reconciles against **what shipped**, not the arch sketch. Recorded deviations, all licensed latitude (they match the shipped `BestPathAligner` and keep the design intent):

- **`MarginalAligner::Context` is a GAT (`type Context<'a>: Copy`), passed by value** — the arch §3 sketch has a plain `type Context` taken by `&Self::Context`. Algorithm 6 (C1) needs the repeat's geometry, and the shipped `RepeatContext<'a>` *borrows*; a plain associated type could only name it as `RepeatContext<'static>`, forcing every caller to own its locus data forever. This is exactly the problem plan 1 solved on `BestPathAligner`, the identical shape is used here. `()` serves algorithm 5.
- **A `Sized` supertrait** — the arch §3 sketch has none on `MarginalAligner`, but arch §4 decides "static dispatch only, never `dyn`" for the whole module, and plan 1 added `: Sized` to `BestPathAligner` to make that a compile error rather than a convention. Applied identically here.
- **`read`/`reference` are bare `&[u8]`** (not `ReadBases`) — deliberate and *unchanged* from the arch: the marginal is quality-blind, scoring with one flat error rate (spec §5.1, §7). This is the load-bearing difference from `BestPathAligner`, which carries `ReadBases`.
- **No `Output` associated type** — a marginal is always one `LogProb`, so there is nothing to vary; the asymmetry with `BestPathAligner` (which spans `Alignment`/`RepeatSpan`) is intentional.

No stop-and-ask arose; none of these changes the design.

## 3. Changes made

- **A1** — `src/ng/types.rs`: `pub struct LogProb(pub f64)` beside `Bp`, with `.get()`. A probability held as its natural logarithm; `f64::NEG_INFINITY` is a legal value (ln 0 = impossible), which is the whole point — an impossible line-up reaches a finite sentinel a caller can see, where linear space reaches `0` (indistinguishable from underflow). Unconstrained → public field, no checked constructor (matching `Bp`, unlike the checked `MismatchFraction`). Derives `Copy, Clone, PartialEq, PartialOrd, Debug` — not `Eq`/`Ord`/`Hash`, since `f64` has no total order (matching the float newtype `MismatchFraction`).
- **A2** — `src/ng/alignment/mod.rs`: the `MarginalAligner` trait (`Scratch: Default`, GAT `Context<'a>: Copy`, `Sized` supertrait, `marginal_probability(&self, read, reference, context, scratch) -> LogProb`), plus the `LogProb` import and a module-doc note. No implementations; two `#[cfg(test)]` stand-in aligners anchor the trait shape.

## 4. Tests added / updated

- **A1** — `log_prob_carries_any_logarithm_including_negative_infinity`; and (review-applied) `log_prob_partialord_is_not_a_total_order_for_nan` (pins why the type is `PartialOrd`, not `Ord`: two NaNs `partial_cmp` to `None`) and `log_prob_carries_positive_infinity_out_of_domain` (out-of-domain floats carried verbatim).
- **A2** — `marginal_aligner_returns_the_negative_infinity_sentinel_for_an_unreachable_line_up` (the `()` context + the `-∞` sentinel), `marginal_aligner_supports_a_borrowed_repeat_context` (the GAT proof: a borrowed `RepeatContext<'a>` is implementable), and (review-applied) `marginal_aligner_copy_context_and_default_scratch_bind_through_a_generic_caller` with the helper `marginal_reusing_context<A: MarginalAligner>` — which makes the `Context: Copy` and `Scratch: Default` bounds *load-bearing* (a concrete impl compiles regardless; only a generic caller building scratch via `Default` and reusing a moved `Copy` context binds them).

## 5. Review outcomes

- **A1** (5 category sub-agents): 0 Blocker / 0 Major / **1 Minor** (reliability: `+∞`/`NaN` boundary classes untested; `NaN` load-bearing). Four agents independently confirmed the NaN representability is a sanctioned unconstrained-design choice, not a defect. **Applied** — two boundary tests. naming/defaults/idiomatic/smells/errors/refactor_safety: no findings.
- **A2** (4 category sub-agents): 0 Blocker / 0 Major / **2 Minor + 1 Nit**.
  - reliability (2 Minor): the anchor test claimed to guard `Context: Copy` and `Scratch: Default` but did not — it used concrete types, and both bounds only bite through a generic caller (`RepeatContext` derives `Copy` independently; scratch was `()`). The classic "test that cannot fail." **Applied** — the generic helper above pins both.
  - naming (1 Nit): the two stand-ins were named on different axes. **Applied** — `RepeatAwareMarginal` → `RepeatContextMarginal`, parallel with `UnitContextMarginal` (both name their `Context` type).
  - idiomatic / refactor_safety / defaults / smells / errors: no findings (GAT, `Sized`, dropped-`Output` asymmetry all confirmed correct).

Every finding across both steps was `Applied`; none disputed or deferred.

## 6. Validation (container, via `./scripts/dev.sh`)

- `cargo fmt --check` → clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` → clean.
- `cargo test --lib` → **2226 passed, 0 failed**, 4 ignored.

The two project-wide validation commands that are red for pre-existing, unrelated reasons (`cargo test --all-targets` bench panic at `benches/psp_writer_perf.rs:386`; `cargo doc --no-deps` intra-doc links) are the documented Standing items; not touched by this work.

## 7. Tradeoffs / follow-ups

- The `MarginalAligner` contract has properties no test can reach until an implementation exists (the actual summed-probability behaviour, banding, parity against `align_subst`). Those are owed by B1–B3 (algorithm 5) and C1 (algorithm 6) — the same "TODO on the trait until a producer exists" shape plan 1 used for `BestPathAligner`.
- Reporting cadence for this plan: full per-step context lives in each commit message; the per-category review files stay in gitignored `tmp/`; this milestone-consolidated report + the `PROJECT_STATUS` block are the durable record (matching plan 1's milestone-level reporting).

**Checkpoint A** — `LogProb` exists in the shared vocabulary and the trait compiles with no implementations. Hard pause for human review before Milestone B.
