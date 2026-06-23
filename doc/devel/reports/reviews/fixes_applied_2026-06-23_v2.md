# Fix Application Report: ssr-call genotyping+pre-pass — Milestone B review

**Date:** 2026-06-23
**Source review:** `doc/devel/reports/reviews/ssr_call_genotyping_milestone_b_2026-06-23.md`
**Source state reviewed against:** commit `d324a0c`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 3 · Nits: 2 · Missing tests: 1

### Outcome totals
- Applied: 6 (Mi1, Mi2, Mi3, both Nits, the missing test) · Deferred/Disputed: 0

### Validation summary
- `cargo fmt --check` → pass
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-features` → 0, **1207 lib pass** (+1), integration + doctests green
- Performance check → Skipped — no `benches/` harness covers Stage-2 cohort code;
  changes are a guard, doc comments, and a test.

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | `is_clear_peak` `length + 1` overflow | Apply | Applied | `rung_ladder.rs` | Pass |
| Mi2 | Minor | `align_subst` normalization claim | Apply | Applied | `pair_hmm.rs` | Pass |
| Mi3 | Minor | `Merged` conflates merged-het + no-structure | Apply | Applied | `rung_ladder.rs` | Pass |
| Nit-1 | Nit | `representative_sequence` tie-break clarity | Apply | Applied | `rung_ladder.rs` | Pass |
| Nit-2 | Nit | `banded_forward` memory not banded | Apply | Applied | `pair_hmm.rs` | Pass |
| MT-1 | — | missing pure-contraction test | Apply | Applied | `stutter.rs` | Pass |

## 3. Questions asked and answers
None — Q4-1 (align normalization) resolved by the Mi2 documentation fix (kept as a v1
behavior with an F2-refinement note, per the review's own recommendation).

## 4. Per-finding log

### Mi1 — `is_clear_peak` `length + 1` overflow
- **Final status:** Applied — `length.checked_add(1).and_then(...).unwrap_or(0)`, mirroring
  the existing `checked_sub(1)` lower-neighbour lookup. No interior behavior change
  (lengths are well below `u16::MAX`); removes the latent debug panic.
- **Files:** [rung_ladder.rs](../../../../src/ssr/cohort/rung_ladder.rs)

### Mi2 — `align_subst` normalization claim
- **Final status:** Applied — documented that the equal-length closed form is exactly
  normalized and the flank-gap branch is an unnormalized boundary-slop correction
  (proper transition normalization = F2). No code change; the math is unchanged.
- **Files:** [pair_hmm.rs](../../../../src/ssr/cohort/pair_hmm.rs)

### Mi3 — `Merged` conflation
- **Final status:** Applied — documented on the zero-peak arm that a structureless locus
  shares the `Merged` reason in v1 (both mean "not a seed"); a `NoStructure` split is
  deferred to D1 diagnostics if needed.
- **Files:** [rung_ladder.rs](../../../../src/ssr/cohort/rung_ladder.rs)

### Nit-1 / Nit-2 — clarity + memory note
- **Final status:** Applied — a comment naming the `representative_sequence` tie-break
  rule; a `banded_forward` note that memory-banding is an F1 concern.

### MT-1 — missing pure-contraction test
- **Final status:** Applied — `reach_variants_pure_contraction_shortens_the_tiling`
  (`(CA)×5` with `Δ = −2` → `(CA)×3`). Confirms the successful contraction path the
  prior tests (`+2`, un-contractable `−3`) did not cover.
- **Files:** [stutter.rs](../../../../src/ssr/cohort/stutter.rs)

## 5–8. Deferred / Disputed / Failed / Blocked
None.

## 9. Performance check
Skipped — no `Apply` touched perf-sensitive (`benches/`-reachable) code.

## 10. Commands run
- `cargo fmt`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-features`

## 11. Command results
- `cargo clippy …` → 0, clean
- `cargo test --all-features` → 0, 1207 lib pass

## 12. Notes
- Pre-existing `benches/psp_writer_perf.rs:386` panic remains out of scope.
