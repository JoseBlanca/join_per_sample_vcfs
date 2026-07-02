# Fix Application Report: coverage_accumulator_2026-06-29.md

**Date:** 2026-06-29
**Source review:** `doc/devel/reports/reviews/coverage_accumulator_2026-06-29.md`
**Source state:** branch `tomato2-paralog-filter`, B2 working diff
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers 0; Majors 3 (M1–M3); Minors 3 (Mi1–Mi3).

### Outcome totals
- Applied: 5 (M1, M2, M3, Mi1, Mi2)
- Deferred: 1 (Mi3 — shared-type unification)

### Validation summary
- `cargo fmt --check` → exit 0
- `cargo clippy --lib -- -D warnings` → clean on scope (pre-existing
  `vcf/writer.rs` lints only)
- `cargo test --lib` → 1379 passed, 0 failed (+11 for B2, +4 vs pre-fix)
- `cargo doc --no-deps` → pre-existing `ClassicStutterModel` failure only
- `cargo audit` → not run (no dependency change)

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files |
|---|---|---|---|---|---|
| M1 | Major | re-entered tile double-counts | Apply | Applied | coverage.rs |
| M2 | Major | new() debug/release divergence | Apply | Applied | coverage.rs |
| M3 | Major | expect() missing PANIC-FREE | Apply | Applied | coverage.rs |
| Mi1 | Minor | depth_sum overflow policy | Apply | Applied | coverage.rs |
| Mi2 | Minor | dead `mean_depth <= 0` branch | Apply | Applied | coverage.rs |
| Mi3 | Minor | scheme/histogram field dup | Defer | Deferred | None |

## 3. Questions asked and answers
None.

## 4. Per-finding log (key items)

### M1 — re-entered tile double-counts
Applied. Added `last_observed: Option<(u32, u32)>` and a `debug_assert!`
that `(chrom_id, pos)` is non-decreasing (the walker's coordinate-order
invariant). Out-of-order input panics in debug; the doc-comment now
states the contract instead of "treated as a new tile". Test
`out_of_order_observe_panics_in_debug`.

### M2 — `new()` debug/release divergence
Applied. Replaced the `debug_assert!` + silent release clamp with
`assert!` in both profiles; a non-positive scheme now fails loudly. The
scheme is CLI-validated config (C1), so an invalid one here is a
programmer error. Test `zero_gc_bins_scheme_panics`; `# Panics` doc added.

### M3 — `expect()` panic discipline
Applied. Added the `// PANIC-FREE:` comment explaining `self.open` is
`Some` on every path out of the preceding match.

### Mi1 / Mi2
- Mi1: `tile.depth_sum = tile.depth_sum.saturating_add(u64::from(depth))`.
- Mi2: removed the unreachable `mean_depth <= 0.0` special case; the
  general `(mean/width) as usize` maps 0 to bin 0 and saturates negatives.

## 5. Deferred findings to carry forward
- **Mi3** — unify `CoverageBinScheme` with the four bin-scheme fields of
  `CoverageByGcHistogram` (e.g. `#[serde(flatten)]`). Deferred: distinct
  layers (input config vs output data), and unification would change B1's
  TOML wire shape for marginal benefit. Revisit if a third copy appears.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 9. Performance check
Skipped — no `benches/`-covered hot path (the accumulator is not yet
wired into the pileup; that is C2, which will get a cost check in D3).

## 10. Commands run
- `./scripts/dev.sh cargo test --lib sample_summary::coverage`
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings [-A clippy::doc_lazy_continuation]`

## 11. Command results
- `cargo test --lib` → 1379 passed, 0 failed
- `cargo fmt --check` → exit 0
- `cargo clippy` → 5 pre-existing vcf/writer.rs errors only

## 12. Notes
- The re-entered-tile contract is now "reject (debug-assert) out-of-order
  input" rather than "silently treat as a new tile" — the safer choice
  for a strictly coordinate-ordered single-thread producer.
