# Code Review: diversity_estimate (SFS genotype prior — Step A)

**Date:** 2026-07-04
**Reviewer:** rust-code-review skill (orchestrator, 8 parallel category agents)
**Scope:** the new `src/var_calling/diversity.rs` pure estimator + its one-line
`mod.rs` registration (Step A of
[../../implementation_plans/sfs_genotype_prior.md](../../implementation_plans/sfs_genotype_prior.md)).
**Status:** Approve-with-changes

---

## 1. Scope

- Reviewed: diff on branch `sfs-genotype-prior` — new module `src/var_calling/diversity.rs`, registration in `src/var_calling/mod.rs`. Read-only context: `src/sample_summary/mod.rs`, `src/paralog/inbreeding.rs`.
- Categories dispatched (all clean-or-minor): reliability, errors, naming, idiomatic, refactor_safety, smells, defaults, module_structure. Skipped: unsafe_concurrency (no unsafe/threads), tooling (no Cargo.toml change), extras (not a parser/hot-path/public crate).
- Audit trail: `tmp/review_2026-07-04_diversity/`.

## 2. Verdict

**Approve-with-changes.** The estimator is correct on the estimated path, well-documented, and its 8 tests cover the inbreeding-free property, human-scale scaling, pooled-vs-mean, the fallbacks, and the F-clamp corners. One Major (an unenforced documented invariant) and a handful of Minor cleanups.

## 3. Execution status

- `cargo test --lib var_calling::diversity` → 8 passed, 0 failed.
- `rustfmt --check src/var_calling/diversity.rs` → clean.
- `cargo clippy --lib` → no warnings on `diversity.rs`. (Crate-wide `-D warnings` is blocked by a pre-existing `needless_lifetimes` in `posterior_engine/shape.rs:384` and rustfmt drift in `write_pass.rs`/`vcf_writer.rs` — untouched code, container clippy/rustfmt 1.95; out of scope.)

## 4. Findings

### M1 — `theta`'s documented "Always >= 0" invariant is unenforced on the override/fallback paths
**Categories:** reliability, errors, smells (convergent). `diversity.rs` — the estimated path is provably non-negative, but `CliOverride` and `PriorFallback` pass `cli_override` / `prior_theta` straight into `theta` with no finiteness/sign check, so a `NaN`/`-inf`/negative value silently seeds the downstream SFS prior. No test exercises malformed override/prior input.
**Fix:** state the caller precondition on the two `f64` inputs in the doc, add `debug_assert!(… .is_finite() && … >= 0.0)` guards, and add a finiteness regression test. (Full CLI-layer validation lands in Step D; this makes the pure estimator's contract explicit and debug-caught now.)

### Mi1 — `inbreeding: Option<Vec<f64>>` is always `Some` (speculative optionality)
**Categories:** idiomatic, smells, refactor_safety (convergent). All return paths set `Some(...)`, so `None` is unreachable, yet every reader must `.as_ref().unwrap()`. The deferred CLI-`F` path the doc cites does not exist yet and would be better modeled when it lands.
**Fix:** change the field to `Vec<f64>`; drop the `Option` wrapping and the `.unwrap()`s in the tests.

### Mi2 — summary reads are not field-addition-safe
**Category:** refactor_safety. `HetCounts`/`CoverageByGcHistogram`/`SampleSummary` are `#[non_exhaustive]` (fields land additively). A future copy-bearing genotype class or a renamed callable denominator would compile silently and bias `θ̂`/`F`.
**Fix:** exhaustively destructure the summary structs at the two read sites (named `_` for ignored fields), so a new field becomes a compile error at the spot that must decide whether it participates.

### Mi3 — missing boundary / invariant / ordering tests
**Category:** reliability. No test pins the `< MIN_COPIES_FOR_ESTIMATE` boundary (an off-by-one `<`→`<=` would pass), the `theta` finiteness invariant (M1), or per-sample `F` ordering across a multi-sample cohort (existing F tests are single-sample or `.len()` only).
**Fix:** add `from_summaries_estimates_at_exactly_min_copies`, `from_summaries_theta_is_finite_and_nonnegative_for_all_sources`, `from_summaries_preserves_inbreeding_order_across_samples`.

### Mi4 — names carry the symbol, not the value
**Category:** naming. Per the project's clear-writing convention (names say what the value *is*; single Greek letters stay in derivations, not the API): `theta` → `nucleotide_diversity`, `inbreeding` → `inbreeding_coefficients`, and the function `per_sample_inbreeding` (a verb-less noun phrase whose `per_sample` prefix misreads its single-sample input) → `estimate_inbreeding`. θ stays in the docs.
**Fix:** rename; keep the symbol in prose.

### Mi5 — shared `F`-ceiling constant is a cross-module coupling
**Categories:** defaults, module_structure (convergent). Reusing `paralog::inbreeding::MAX_INBREEDING_COEFFICIENT` is sound single-source-of-truth (both are the same `[0, 0.99]` inbreeding range and should not drift), but the ownership sits in the paralog algorithm module — a mild layering surprise.
**Fix:** announce the deliberate reuse with a one-line comment at the `use` site. (Hoisting the constant to a neutral home beside `HetCounts` is a reasonable follow-up but reaches into the paralog module — deferred.)

### Mi6 — module filed under the wrong `mod.rs` section
**Category:** module_structure. `pub mod diversity;` sits under the "Numeric kernels (byte-identity-sensitive math)" header, whose contract (byte-carried math, deps limited to `psp`/`fasta`/`pileup_record`) it does not meet — it is a pre-calling estimator depending on `sample_summary`/`paralog`.
**Fix:** move the registration to a fitting section (a new "pre-calling estimation" grouping or the general block).

### Nits
- The diploid factor `2` at the copy-count and chromosome-denominator sites is unnamed (and means two different things — alt copies per hom-alt vs chromosomes per diploid position); clarify with inline comments rather than a single conflating const.
- Test module redundantly re-imports `HetCounts` (already in scope via `use super::*`).
- Test fixture uses `version: 0` (reserved-invalid); use `SAMPLE_SUMMARY_VERSION`.
- Loop binding `s` → `summary`.
- `source` → `theta_source` was suggested; kept as `source` (idiomatic field-named-after-type, unambiguous as `estimate.source`).

## 5. Out of scope observations

- Pre-existing container toolchain drift (rustfmt/clippy 1.95) fails `cargo fmt --check` and `clippy -D warnings` on untouched files (`write_pass.rs`, `vcf_writer.rs`, `posterior_engine/shape.rs`). Separate cleanup; not introduced here.

## 6. What's good

- The copy-counting estimator is genuinely inbreeding-free and the module documents *why* (spec §5 reasoning inline) — `theta_is_inbreeding_free` proves it as a test.
- `MIN_COPIES_FOR_ESTIMATE` is a named `pub const` with a `1/√copies` rationale, referenced by name in the API docs so value and doc cannot drift.
- `prior_theta` as a required parameter (no hidden default) and `DiversitySource` provenance make the fallback behaviour visible at the call site.
- New public enums correctly left *without* `#[non_exhaustive]` (internal crate — keeps future `match` arms compile-forced).
