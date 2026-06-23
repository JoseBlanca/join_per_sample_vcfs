# Code Review: ssr-call genotyping+pre-pass — Milestone D (D1+D2)
**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator, focused inline pass)
**Scope:** Milestone D D1+D2 (commit `4c3d756`, `prepass.rs`)
**Status:** Approve-with-changes

---

## 1. Scope
- **In-scope:** [prepass.rs](../../../../src/ssr/cohort/prepass.rs) (new) + the `mod.rs`
  wiring.
- **Out of scope:** D3 (post-gate), the rest of the cohort (reviewed).
- **Categories:** reliability, errors, naming, idiomatic, refactor_safety, smells, extras
  (numerical correctness / determinism). Skipped unsafe_concurrency, tooling.

## 2. Verdict
**Approve-with-changes.** Checkpoint 2 holds (parameters recovered), and the sufficient
statistics are integer counts so the single-pass reduce is order-independent by
construction (`ε`, the per-period shape, and the per-sample level are all deterministic
despite `HashMap` iteration, since each is a keyed/integer reduction). Findings are
documentation Minors, one unclamped-output Minor, and two missing tests.

## 3. Execution status
- `cargo fmt --check` → pass · `cargo clippy --all-targets --all-features -- -D warnings`
  → pass · `cargo test --all-features` → **1236 lib pass**.
- Needs-verification findings: 0.

## 4. Top 3 priorities
1. **Mi1** — `ε` taken off a het's faithful peak is mildly contaminated (the other
   allele's rare same-length stutter + impure same-length variants counted as mismatch);
   document the small upward bias.
2. **Mi2** — `fit_level` can return a `baseline`/`slope` line that extrapolates outside
   `[0,1]`; document that the consumer clamps (or clamp `baseline`).
3. **MT-1/MT-2** — add a separated-het contribution test and a determinism test.

## 5. Findings

### Minor

- `src/ssr/cohort/prepass.rs` (`accumulate_locus`) — **[Minor]** `ε` off het peaks is slightly biased high
- **Confidence:** Medium
- **Problem:** For a separated het, faithful reads at allele `A`'s length are compared
  to `A`'s representative sequence. Two small contaminants count as mismatch: the *other*
  allele's rare `≥2`-step stutter landing at `A`'s length, and any impure same-length
  variant (a different real sequence at the same length). Both nudge the base-mismatch
  fraction up.
- **Why it matters:** A small upward `ε` bias on het-heavy cohorts; negligible for the
  homozygote-dominated recovery test, real to document before the soft-EM reduce replaces
  the hard labels.
- **Suggested fix:** document the bias on `accumulate_locus` / `compare_bases`; the soft
  full-cohort EM reduce (deferred) removes it by fractional attribution.

- `src/ssr/cohort/prepass.rs` (`fit_level`) — **[Minor]** the fitted level line is unclamped
- **Confidence:** High
- **Problem:** Weighted-LS `baseline + slope·length` can extrapolate to a negative or
  `>1` rate at lengths outside the observed range; `StutterLevel` stores it raw.
- **Why it matters:** A negative/`>1` stutter probability is unphysical; it happens to be
  clamped at evaluation (`em.rs`, `sim.rs` clamp `level` to `[0,1]`), so it is latent.
- **Suggested fix:** document that the line is a seed clamped at evaluation; or clamp
  `baseline` into `[0,1)` at fit time. (Documented as the chosen resolution.)

### Nits
- `src/ssr/cohort/prepass.rs` (`accumulate_locus`) — the `rungs: &crate::ssr::cohort::rung_ladder::Rungs`
  parameter uses a fully-qualified inline path; import `Rungs` for readability.
- `accumulate_locus` silently drops slips with `|Δ| > MAX_SLIP` (spec-capped); a one-line
  comment noting the cap is deliberate would help (no counter needed in v1).

## 6. Out of scope observations
- `benches/psp_writer_perf.rs:386` — pre-existing bench panic, unchanged.

## 7. Missing tests to add now
- `separated_hets_contribute_two_length_bins` — **input:** a cohort of cohort-recurrent
  separated hets (e.g. `6/10`). **Bug it catches:** an accumulation that mis-attributes
  het reads or only records one allele's bin would skew the level line. **Body:** run the
  pre-pass and assert the per-sample stats carry support at both allele lengths and `ε` is
  still recovered.
- `run_prepass_is_deterministic` — **input:** the recovery cohort, run twice. **Bug it
  catches:** an accidental order-dependence (e.g. a float sum over `HashMap` values).
  **Body:** assert the two `EstimatedParams` are equal.

## 8. What's good
- The estimators are the exact inverses of the simulator's forward draws (decay MLE ↔ the
  geometric magnitude; `ε` ↔ the per-base substitution), and checkpoint 2 pins all three
  ([prepass.rs](../../../../src/ssr/cohort/prepass.rs)).
- Local-accumulate-then-merge sidesteps the two-`&mut`-into-`stats` borrow cleanly.
- The module header is explicit about which heavier-machinery refinements are deferred and
  why — no hidden under-delivery.

## 9. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --all-targets --all-features -- -D warnings` ·
  `cargo test --all-features`

### Author response convention
Address each finding by ID (Mi1, Mi2, Nits, MT-1, MT-2).
