# Fix Application Report: genetics_extract_2026-07-04.md

**Date:** 2026-07-04
**Source review:** `doc/devel/reports/reviews/genetics_extract_2026-07-04.md`
**Source state reviewed against:** branch `sfs-genotype-prior`, Step B0 (genetics extraction)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

- Review totals: 0 Blocker · 0 Major · 2 Minor · 3 Nit (2 of 3 categories: No findings).
- Outcome: Applied 3 · Disputed 2 (kept with rationale) · No findings 2 categories.
- Validation: `cargo test --lib` → 1533 passed; `genetics` → 2; `paralog` → 98 (incl. the numerical byte-identity anchor); integration suites pass; `rustfmt --edition 2024 --check src/genetics.rs` clean; no clippy warnings on changed files.
- Unresolved high-priority findings: none.

## 2. Findings table

| ID | Sev | Title | Decision | Final status | Files |
|---|---|---|---|---|---|
| module_structure | — | (whole category) | — | No findings | — |
| refactor_safety | — | (whole category) | — | No findings | — |
| Mi1 | Minor | module name `genetics` too broad | Dispute | Disputed | — |
| Mi2 | Minor | `linspace_point` NumPy jargon | Apply | Applied | genetics.rs, locus_score.rs |
| Nit-a | Nit | `PROB_FLOOR` abbreviates "probability" | Apply | Applied | genetics.rs, locus_score.rs |
| Nit-b | Nit | `inv2n` cryptic param | Apply | Applied | genetics.rs |
| Nit-c | Nit | `SFS` acronym | Dispute | Disputed | — |

## 3. Per-finding log (condensed)

- **Mi2 — Applied.** `linspace_point` → `linear_grid_point` in `genetics.rs` (def + `sfs_grid_point` call + test) and the paralog call site + import. Byte-identity-neutral (name only).
- **Nit-a — Applied.** `PROB_FLOOR` → `PROBABILITY_FLOOR` in `genetics.rs` (const + 3 internal uses) and the two paralog use sites + import. Value `1e-300` unchanged.
- **Nit-b — Applied.** `sfs_grid_point` parameter `inv2n` → `grid_inset`, with the doc updated (the doc already described it as "the grid inset"). Purely local.
- **Mi1 — Disputed.** Module named `genetics` per the project owner's explicit request; kept.
- **Nit-c — Disputed.** `SFS` is standard population-genetics vocabulary and is expanded in the module doc; kept.

## 4. Deferred / failed / blocked
None.

## 5. Performance check
Skipped — Step B0 is a verbatim move (no algorithmic change) and no new bench-reachable code.

## 6. Notes
- All three renames are byte-identity-neutral (identifiers only, no float ops changed); the paralog byte-identity anchor (`reduced_grid_matches_hand_computed_likelihoods`) still passes.
- Pre-existing container toolchain drift on untouched files distinguished throughout and left alone.
