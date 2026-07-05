# Fix Application Report: diversity_estimate_2026-07-04.md

**Date:** 2026-07-04
**Source review:** `doc/devel/reports/reviews/diversity_estimate_2026-07-04.md`
**Source state reviewed against:** branch `sfs-genotype-prior`, Step A (`src/var_calling/diversity.rs`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 1 (M1) · Minors: 6 (Mi1–Mi6) · Nits: 5

### Outcome totals
- Applied: 11 · Disputed: 1 (Nit: `source` rename) · all others resolved

### Validation summary
- `rustfmt --edition 2024 --check src/var_calling/diversity.rs` → 0, clean (the mod-tree walk from `mod.rs` reports only pre-existing sibling drift in `write_pass.rs`/`vcf_writer.rs`, unchanged here).
- `cargo test --lib var_calling::diversity` → 0, 10 passed (debug); release build runs the added finiteness guard, 1 passed.
- `cargo test --lib` → 0, **1533 passed, 0 failed**, 2 ignored.
- `cargo clippy --lib` → no warnings on `diversity.rs`. Crate-wide `-D warnings` remains blocked by a pre-existing `needless_lifetimes` in `posterior_engine/shape.rs:384` (untouched; container clippy 1.95).

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Sev | Title | Decision | Final status | Files | Validation |
|---|---|---|---|---|---|---|
| M1 | Major | `theta` "Always ≥0" invariant unenforced on override/fallback | Apply | Applied | diversity.rs | Pass |
| Mi1 | Minor | `inbreeding: Option<Vec>` always `Some` | Apply | Applied | diversity.rs | Pass |
| Mi2 | Minor | Non-field-safe summary reads (`#[non_exhaustive]`) | Apply | Applied | diversity.rs | Pass |
| Mi3 | Minor | Missing boundary/finiteness/ordering tests | Apply | Applied | diversity.rs | Pass |
| Mi4 | Minor | Names carry the symbol, not the value | Apply | Applied | diversity.rs | Pass |
| Mi5 | Minor | Shared F-ceiling cross-module coupling | Apply | Applied | diversity.rs | Pass |
| Mi6 | Minor | Module under wrong `mod.rs` section | Apply | Applied | mod.rs | Pass |
| Nit-a | Nit | Diploid `2` unnamed | Apply | Applied | diversity.rs | Pass |
| Nit-b | Nit | Redundant `HetCounts` test import | Apply | Applied | diversity.rs | Pass |
| Nit-c | Nit | Fixture `version: 0` (reserved) | Apply | Applied | diversity.rs | Pass |
| Nit-d | Nit | Loop binding `s` | Apply | Applied | diversity.rs | Pass |
| Nit-e | Nit | `source` → `theta_source` | Dispute | Disputed | — | — |

## 3. Questions asked and answers
None (all findings were clear `Apply`, one Nit disputed with rationale).

## 4. Per-finding log (condensed)

- **M1 — theta finiteness/sign.** Applied. Documented the caller precondition on `prior_theta`/`cli_override`, added `debug_assert!(is_finite() && >= 0.0)` on both, softened the field doc to state the guarantee holds under those preconditions, and added `theta_is_finite_and_nonnegative_for_all_sources` (release-only — debug asserts fire first). Kept the estimator infallible: the inputs are a compile-time constant and a CLI-validated flag (validation lands in Step D), so a `Result` would be over-engineering for a pure numeric fn.
- **Mi1 — Option→Vec.** Applied. `inbreeding: Option<Vec<f64>>` → `inbreeding_coefficients: Vec<f64>`; removed the five `.as_ref().unwrap()` in tests. The deferred cohort-wide-CLI-`F` path will be modeled when it lands (Step C/D), not pre-represented.
- **Mi2 — field-safety.** Applied. Added `alt_allele_copies(&HetCounts)` and `callable_position_count(&CoverageByGcHistogram)`, each **exhaustively destructuring** (named `_` for ignored fields); `estimate_inbreeding` does the same for `HetCounts`. A new copy-bearing count or renamed callable field is now a compile error at the read site.
- **Mi3 — tests.** Applied. Added `estimates_at_exactly_min_copies` (pins the `<` boundary), `theta_is_finite_and_nonnegative_for_all_sources` (M1), `preserves_inbreeding_order_across_samples` (multi-sample F ordering). 10 debug + 1 release now.
- **Mi4 — naming.** Applied. `theta` → `nucleotide_diversity`; `inbreeding` → `inbreeding_coefficients`; fn `per_sample_inbreeding` → `estimate_inbreeding` (verb). θ symbol kept in the docs/derivations per the clear-writing skill.
- **Mi5 — F-ceiling coupling.** Applied. Added a comment at the reuse site stating the deliberate single-source-of-truth intent (paralog scorer + this estimator share one `[0, 0.99]` clamp). Hoisting the constant to a neutral home is noted as a follow-up (reaches into the paralog module).
- **Mi6 — mod.rs placement.** Applied. Moved `pub mod diversity;` out of the "Numeric kernels (byte-identity-sensitive math)" block into a new "Pre-calling cohort estimation" section with a note that it is not part of the byte-identity math contract.
- **Nit-a.** Applied. Introduced `CHROMOSOMES_PER_DIPLOID` and `ALT_COPIES_PER_HOM_ALT` (deliberately two constants despite the equal value — they mean different things: chromosomes per position vs alt copies per hom-alt genotype).
- **Nit-b/c/d.** Applied. Dropped the redundant `HetCounts` test import; fixture now uses `SAMPLE_SUMMARY_VERSION`; loop binding `s` → `summary`.
- **Nit-e — `source` → `theta_source`.** Disputed. Kept `source`: it is the idiomatic field-named-after-its-type pattern (`source: DiversitySource`), unambiguous as `estimate.source`. Renaming would add noise without clarity.

## 5. Deferred findings to carry forward
- Hoisting `MAX_INBREEDING_COEFFICIENT` to a neutral home (Mi5 follow-up) — deferred; reaches into the committed paralog module, out of scope for Step A.

## 6. Disputed findings to return to reviewer
- Nit-e (`source` rename) — kept as idiomatic; see per-finding log.

## 7–8. Failed-validation / blocked-by-context-mismatch
None.

## 9. Performance check
Skipped — Step A is a pure pre-calling estimator not reachable from any `benches/` harness, and it has no caller yet.

## 10. Commands run
- `rustfmt --edition 2024 --check src/var_calling/diversity.rs`
- `cargo test --lib var_calling::diversity` (debug + `--release`)
- `cargo test --lib`
- `cargo clippy --lib --all-features`

## 12. Notes
- Pre-existing container toolchain drift (rustfmt/clippy 1.95) fails whole-crate `fmt --check` / `clippy -D warnings` on untouched files; distinguished from this change throughout and left alone.
