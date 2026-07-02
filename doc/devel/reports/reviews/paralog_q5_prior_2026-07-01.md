# Code review — hidden-paralog filter Q5 (`ParalogPrior` — empirical-Bayes prior + FDR)

## 1. Scope
- **Reviewed:** the new `src/paralog/prior.rs` (LR histogram, EM prior, tail-FDR
  curve) + a re-export line in `src/paralog/mod.rs`.
- **Against:** branch `tomato2-paralog-filter`, uncommitted diff.

## 2. Method
Two parallel review sub-agents per `ai/skills/rust-code-review`: (a) reliability
+ refactor_safety with a **numerical-correctness / prototype-faithfulness**
emphasis; (b) naming + errors + idiomatic + defaults + smells + module_structure.
Verification (dev container): `cargo test --lib paralog` 65 passed (was 60; +5
challenge/boundary tests); `cargo clippy --lib --tests -D warnings` clean;
`cargo fmt --check` clean; 1463 lib tests pass.

## 3. Verdict
**Numerically faithful and sound after fixes.** The numerical reviewer confirmed
the histogram EM equals the prototype's full-vector `π = mean(σ(LR + logit π))`
up to Jensen-order `O(binwidth²)`, the tail-FDR-from-the-top reproduces the
prototype's rank-by-posterior running mean exactly (posterior strictly monotone
in LR ⇒ LR order = rank order), and the saturation / sigmoid-stability /
logit-clamp / `−logit π` half-cut are all correct. No algorithm bug; findings
were an API-convention inconsistency, a duplication, a silent-signal gap, and
missing boundary tests.

## 4. Findings (all fixed)

- **M — constructor-contract inconsistency.** `ParalogLrHistogram::new` panicked
  via `assert!` while its same-module sibling `GridSpec::new` returns `Option`.
  → **Fixed:** `new` now returns `Option<Self>` (matching the module
  convention); `with_defaults` keeps the trusted compile-time-consts path.
- **M — EM non-convergence was silent.** `estimate` returned the last iterate
  with no signal, indistinguishable from a converged one; the `max_iter`
  branch was untested. → **Fixed:** added `ParalogPrior.converged: bool` (also
  `false` for an empty histogram), and a test driving the non-convergent branch
  (`max_iter = 1`).
- **Mi (both agents) — duplicated `bin_index`/`bin_center`** across the
  histogram and the curve, the load-bearing "same binning" invariant held by two
  hand-synced copies. → **Fixed:** extracted a shared private `UniformBinning
  { lo, hi, n_bins }`; both types delegate, and the curve carries the
  histogram's `binning` by value (documented independence).
- **Mi — `paralog_posterior` param naming drift** (`paralog_log_likelihood_ratio`
  vs `lr` everywhere else). → **Fixed:** renamed to `lr`.
- **Mi — empty-histogram silent default** (π = start). → **Fixed:** surfaced via
  `converged = false` (inspectable) rather than a silent seed value.
- **Deferred/skipped:** the `1e-12` clamp literal (local, self-documenting — a
  nit); a `tracing` log on non-convergence (the crate has no logging facade; the
  `converged` flag is the inspectable signal the wiring will act on).

## 5. Tests added (5)
`new` rejects a degenerate range; EM flags non-convergence (`max_iter = 1`) vs a
converged run; `lr_threshold_for_fdr` `None`/all-qualify boundaries; degenerate
π (`0`/`1`) stays finite through the posterior + half-cut; the empty-histogram
FDR curve is finite everywhere (the empty-tail `local_lfdr` fallback). These
join the existing EM-recovery, histogram-vs-full-vector parity, monotonicity,
and threshold-relaxation tests.

## 6. What's good
The histogram/EM/FDR split is a clean bounded-RAM realisation of arch Premise 6;
documented `pub const` defaults; module depends on nothing outside itself (no
`use crate::` lines); `sigmoid`/`logit`/`clamp_probability` are correct,
`−∞`-safe pure helpers.

Audit trail: `tmp/review_2026-07-01_paralog-q5/{reliability,naming_etc}.md`.
