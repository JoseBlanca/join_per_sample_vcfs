# Code Review: genetics primitives extraction (SFS prior — Step B0)

**Date:** 2026-07-04
**Reviewer:** rust-code-review skill (orchestrator, proportionate 3-category dispatch)
**Scope:** verbatim extraction of shared population-genetics primitives from
`src/paralog/locus_score.rs` into a new `src/genetics.rs` (Step B0 —
prerequisite for `SfsGenotypePrior`). Branch `sfs-genotype-prior`.
**Status:** Approve-with-changes

---

## 1. Scope & category triage

Verbatim move, no new logic, byte-identity-gated. Per the skill's triage, the
risk surface is narrow, so three categories were dispatched (not the full eight):
**module_structure** (new module placement / imports / cycles), **refactor_safety**
(behavior preservation / byte-identity), **naming** (the now-public API). The
"always" categories with no new logic to examine (errors/idiomatic/smells/defaults)
were skipped as inapplicable to a pure move. Audit trail: `tmp/review_2026-07-04_genetics/`.

Moved: `PROBABILITY_FLOOR` (was `PROB_FLOOR`), `linear_grid_point` (was
`linspace_point`), `sfs_grid_point`, `wright_genotype_log_priors`, plus two unit
tests. `paralog/locus_score.rs` imports them back.

## 2. Verdict

**Approve-with-changes** — two categories clean; naming found API-surface Minors.

## 3. Execution status

- `cargo test --lib` → 1533 passed. `genetics` → 2 passed. `paralog` → 98 passed, incl. `reduced_grid_matches_hand_computed_likelihoods` (pins grid+Wright numerics — the byte-identity gate). All integration suites pass. `rustfmt --edition 2024 --check src/genetics.rs` → clean; no clippy warnings on the changed files.

## 4. Findings

- **module_structure — No findings.** Correct top-level peer (consumed by two modules across a boundary), honest non-circular dependency direction (genetics imports nothing), flat `.rs`, imports and `lib.rs` registration correct, no orphaned references.
- **refactor_safety — No findings.** The move is character-for-character verbatim (formula, operation order, the `n <= 1` midpoint branch, the `.max(floor).ln()` clamp, all float literals). Only visibility (`fn`→`pub fn`), doc wording, and test relocation changed. Rust applies no fast-math/FMA contraction, so the private→pub cross-module move cannot perturb IEEE results — byte-identity justified; the existing numerical regression test adequately guards it.
- **naming — 2 Minor + 3 Nits.**
  - **Mi1 (module name `genetics` too broad).** *Disputed / won't-fix* — the module name was chosen by the project owner ("a shared genetics module"). Kept.
  - **Mi2 (`linspace_point` is NumPy jargon).** *Applied* — renamed to `linear_grid_point`.
  - **Nit (`PROB_FLOOR` abbreviates "probability").** *Applied* — `PROBABILITY_FLOOR`.
  - **Nit (`inv2n` cryptic; doc calls it "the grid inset").** *Applied* — parameter renamed `grid_inset`.
  - **Nit (`sfs_grid_point`'s `SFS` acronym).** *Won't-fix* — standard population-genetics vocabulary, expanded in the module doc.
  - Cross-category (noun-form fn names, single-letter `p`/`f`) — *won't-fix*: acceptable under the clear-writing skill's math-symbol exception (these are the standard pop-gen symbols inside a derivation-style function).

## 5–9. Out of scope / tests / what's good

- Pre-existing container toolchain drift (rustfmt/clippy 1.95 on untouched files) — unchanged, out of scope.
- Tests: the two moved unit tests (`wright_genotype_priors_sum_to_one`, `linear_grid_point_degenerate_grid_returns_midpoint`) carry the primitives' coverage; the paralog scorer's `reduced_grid_matches_hand_computed_likelihoods` remains the numerical byte-identity anchor.
- What's good: the module doc expands SFS / Wright / Hardy–Weinberg for the geneticist reader and each fn leads with plain English before the formula.
