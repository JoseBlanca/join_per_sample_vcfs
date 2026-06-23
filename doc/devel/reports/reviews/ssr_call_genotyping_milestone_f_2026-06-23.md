# Code Review: ssr-call genotyping+pre-pass — Milestone F (parallelism + byte-identity)
**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator, focused inline pass)
**Scope:** Milestone F F1 (commit `0c2cd7b`)
**Status:** Approve-with-changes

---

## 1. Scope
- **In-scope:** the `rayon` parallelisation in
  [prepass.rs](../../../../src/ssr/cohort/prepass.rs) (`run_prepass_stats`,
  `PrepassStats::merge`) and [inbreeding.rs](../../../../src/ssr/cohort/inbreeding.rs)
  (`run_cohort_em` per-locus EM), plus the two byte-identity tests.
- **Categories:** reliability, **unsafe_concurrency** (rayon, shared state),
  refactor_safety, idiomatic, smells. F2 calibration is empirical (out of scope for code
  review).

## 2. Verdict
**Approve-with-changes.** The parallelisation is correct: the closures are pure, all
shared state is read-only `&` (no interior mutability, no locks), and the reduces are
integer/`FixedPointAccum` so order-independence is *structural*, not hoped-for. Both
byte-identity tests (1 vs 4 threads) pass, which is the proof the milestone exists to
produce. One documentation Minor.

## 3. Execution status
- `cargo fmt --check` → pass · `cargo clippy --all-targets --all-features -- -D warnings`
  → pass · `cargo test --all-features` → **1253 lib pass**.
- Needs-verification findings: 0.

## 4. Concurrency review (unsafe_concurrency)
- **No `unsafe`, no locks, no atomics.** `run_prepass_stats` uses rayon `fold`/`reduce`
  over a pure per-locus accumulation; `run_cohort_em` uses `par_iter().zip().collect()`
  where the closure captures `&f_per_sample` / `&params` / `&level_per_group`
  immutably and `run_locus_em_with` is a pure function. `f_per_sample` is mutated only
  **between** rounds (sequential), never during the parallel map — no data race.
- **Determinism is structural:** `PrepassStats` fields are `u64` counts (commutative +
  associative merge); the `F` reduce uses `FixedPointAccum`; rayon's indexed `collect`
  preserves locus order. The byte-identity tests confirm it empirically.

## 5. Findings

### Minor

- `src/ssr/cohort/prepass.rs` (`PrepassStats::merge`) — **[Minor]** struct order is not canonical, only the derived params are
- **Confidence:** High
- **Problem:** `SampleStutterStats.by_length` is a `Vec` whose order depends on
  insertion / merge order, so two `PrepassStats` reduced via different thread counts are
  **not** guaranteed `==` as structs (a derived `PartialEq` compares `Vec` order). The
  byte-identity *property* holds for the **derived** `EstimatedParams` (tested), because
  `estimate` reads order-independent sums — but a future test that asserts
  `PrepassStats == PrepassStats` across thread counts would flake.
- **Why it matters:** Prevents someone from writing a determinism test at the wrong layer
  and concluding the reduce is broken when it is not.
- **Suggested fix:** document on `merge` that determinism is guaranteed for the *values*
  (and the derived estimates), not the internal `Vec` ordering; compare `EstimatedParams`
  / `CohortCalls`, not the raw stats, in determinism tests.

## 6. Out of scope observations
- `benches/psp_writer_perf.rs:386` — pre-existing bench panic, unchanged.
- F2 calibration constants are provisional `dev_default`s (a real-data follow-up, by
  design — see the milestone report).

## 7. Missing tests to add now
- None required — the two byte-identity tests cover the parallelised surfaces; the
  end-to-end test (Milestone E) covers composition.

## 8. What's good
- The determinism was designed in from A1 (`FixedPointAccum`) and B/D (integer reduces,
  catalog-index tie-breaks), so F1 parallelised the iteration *without touching the
  algorithm* — and the byte-identity tests prove it
  ([prepass.rs](../../../../src/ssr/cohort/prepass.rs),
  [inbreeding.rs](../../../../src/ssr/cohort/inbreeding.rs)).
- `PrepassStats::merge` is a clean commutative integer fold; no shared mutable state.
- F2's empirical limits are documented honestly rather than fabricated on the simulator.

## 9. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --all-targets --all-features -- -D warnings` ·
  `cargo test --all-features`

### Author response convention
Address the finding by ID (Mi1).
