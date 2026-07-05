# Fix Application Report: sfs_prior_2026-07-04.md

**Date:** 2026-07-04
**Source review:** `doc/devel/reports/reviews/sfs_prior_2026-07-04.md`
**Source state reviewed against:** branch `sfs-genotype-prior`, Step B (`SfsGenotypePrior`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

- Review totals: 0 Blocker · 2 Major · ~8 Minor · 3 Nit (across 5 categories).
- Outcomes: Applied 10 (both Majors + naming ×3 + visibility + dead-branch + defaults ×2 + unused-derive) · Deferred 2 (GenotypeLogPriors struct; proptest) · Declined 1 (`folded_sfs` verb nit).
- Validation: `cargo test --lib var_calling::sfs_prior` → 12 passed; `cargo test --lib` → 1545 passed, 0 failed; `rustfmt --edition 2024 --check` clean; no clippy warnings on the module.
- Unresolved high-priority findings: none (both Majors applied).

## 2. Findings table

| ID | Sev | Title | Decision | Final status | Files |
|---|---|---|---|---|---|
| R-M1 | Major | negative/NaN θ slips guard → silent NaN | Apply | Applied | sfs_prior.rs |
| R-M2 | Major | missing edge tests | Apply | Applied | sfs_prior.rs |
| N-Mi1 | Minor | `theta` param → `nucleotide_diversity` | Apply | Applied | sfs_prior.rs |
| N-Mi2 | Minor | `f` param → `inbreeding_coefficient` | Apply | Applied | sfs_prior.rs |
| N-Mi3 | Minor | `invariant_mass` → `invariant_site_mass` | Apply | Applied | sfs_prior.rs |
| S/I-Mi | Minor | `FrequencyGrid`/`folded_sfs` should be private | Apply | Applied | sfs_prior.rs |
| S/I-Mi | Minor | dead `else { 0.5 }` inset branch | Apply | Applied | sfs_prior.rs |
| D-Mi1 | Minor | `DEFAULT_SFS_GRID_POINTS` unenforced "matches" doc | Apply | Applied | sfs_prior.rs |
| D-Mi2 | Minor | `DEFAULT_INVARIANT_SITE_MASS` provisional | Apply | Applied | sfs_prior.rs |
| I-Mi | Minor | `[f64;3]` → `GenotypeLogPriors` struct | Defer | Deferred | — |
| R-Mi | Minor | proptest for algebraic laws | Defer | Deferred | — |
| Nit | Nit | unused `PartialEq` on float structs | Apply | Applied | sfs_prior.rs |
| Nit | Nit | `folded_sfs` verb-name | Dispute | Disputed | — |

## 3. Per-finding notes

- **R-M1:** `debug_assert` on `nucleotide_diversity`/`invariant_site_mass` (`new`) and `inbreeding_coefficient` (method), documented preconditions, guard hardened `total <= 0.0` → `!(total > 0.0)` with a corrected comment. Kept infallible per the reviewer (machine-generated hyperparameters; θ≥0 guaranteed by `DiversityEstimate`).
- **R-M2:** 5 tests added (θ=0, mass=0, total=0 guard, F=1 flooring, grid symmetry).
- **Naming/visibility/dead-branch/defaults/derive:** applied as described in the review report §4.
- **Deferred:** the `GenotypeLogPriors` struct (the `[f64;3]` is the engine-consumed shape — revisit at Step D) and the `proptest` (no dep; deterministic parametrized tests cover the laws).
- **Disputed:** `folded_sfs` is an idiomatic constructor name.

## 4. Performance check
Skipped — no bench-reachable code (no consumer yet; pure prior).

## 5. Notes
- Pre-existing container toolchain drift on untouched files distinguished throughout and left alone.
