# Fix Application Report: ssr_call_level_refit_2026-06-24.md

**Date:** 2026-06-24
**Source review:** `doc/devel/reports/reviews/ssr_call_level_refit_2026-06-24.md`
**Source state reviewed against:** commit `fd53330`, branch `ssr-cohort`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 2 (Mi1, Mi2) · Nits: 3

### Outcome totals
- Applied: 2 (Mi1 rename, Mi2 doc) · No action: 3 Nits (all documented / optional)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, **1278 passed, 0 failed, 2 ignored**
- `cargo doc --no-deps` / `cargo audit` → not run (rename + doc-comment changes; no new dep)
- Performance check → not applicable (per-locus EM not benched)

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | `theta_max_rounds` now bounds both refits | Apply | Applied | `em.rs` | Pass |
| Mi2 | Minor | Document the multiplier ↔ group-level hierarchy | Apply | Applied | `em.rs` | Pass |
| Nit-a | Nit | `expected_slipped` float-sum determinism | — | No action (documented) | None | N/A |
| Nit-b | Nit | `LEVEL_MULT_MAX` magic constant | — | No action (documented) | None | N/A |
| Nit-c | Nit | `attribute_locus` discards `_eps` | — | No action (optional) | None | N/A |

## 3. Questions asked and answers
- **OQ1 (scalar per-locus rate deviation):** accepted as the v1 design — `m` is a length-independent overall-rate nudge; the group line keeps the length-dependence. A per-length per-locus rate is out of scope. No change.

## 4. Per-finding log

### Mi1 — rename `theta_max_rounds` → `refit_max_rounds`
- **Final status:** Applied. Renamed the `EmCfg` field (and the `theta_max_rounds_zero_…` test → `refit_max_rounds_zero_…`) via a scoped identifier rename across `em.rs`; updated the field doc to "I1 shape + I2 level." The per-`θ`/`level` `*_shrink`/`*_tol` keep their name-accurate names.
- **Files changed:** `src/ssr/cohort/em.rs`.
- **Validation:** `cargo clippy --lib -D warnings` 0; `cargo test --lib` 1278 passed; the rename touched only `em.rs` (the field is constructed only in `dev_default` + the test).

### Mi2 — hierarchy doc note
- **Final status:** Applied. Added a `HIERARCHY (no double-count)` note to `refit_level_multiplier`: the per-group level is the prior mean (estimated upstream by `reduce_level` from genotypes, not from `m`), and `m` is the per-locus deviation applied only at genotyping; at convergence matched loci have `m ≈ 1`, so the rate is modelled once, hierarchically.
- **Files changed:** `src/ssr/cohort/em.rs` (doc comment).

### Nits — no action
- **Nit-a:** `expected_slipped` is summed in fixed per-locus order and a locus is never split across threads; byte-identity holds and is already noted on `attribute_locus`.
- **Nit-b:** `LEVEL_MULT_MAX = 10.0` is a documented provisional constant (pinned in F2), consistent with the module's dev defaults.
- **Nit-c:** a `sample_level` helper to avoid computing/ignoring `eps` is cosmetic; deferred.

## 5. Deferred findings to carry forward
- None (Nits are documented/optional).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — the per-locus EM is not reachable from a `benches/` harness.

## 10. Commands run
- `cargo fmt --check`
- `cargo clippy --lib --all-features -- -D warnings`
- `cargo test --lib`

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 1278 passed / 0 failed / 2 ignored

## 12. Notes
- With I1 + I2, the per-locus EM now adapts both stutter shape and rate locally, shrunk toward the frozen cohort priors — the HipSTR-style adaptation, with the cohort-pooled anchor HipSTR lacks. The soft per-read responsibility split remains the deferred refinement.
