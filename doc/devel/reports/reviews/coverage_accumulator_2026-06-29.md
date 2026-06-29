# Code Review: coverage-by-GC accumulator (B2)

**Date:** 2026-06-29
**Reviewer:** rust-code-review skill (orchestrator + 1 consolidated 7-category sub-agent)
**Scope:** feature B2 — `src/sample_summary/coverage.rs` + `pub mod coverage;` in `src/sample_summary/mod.rs`
**Status:** Approve-with-changes

---

## 1. Scope
- Reviewed: B2 working diff on branch `tomato2-paralog-filter`.
- In-scope: `src/sample_summary/coverage.rs` (new), `mod.rs` (one line).
- Categories (consolidated for low-risk deterministic math): reliability,
  errors, defaults, idiomatic, naming, smells, refactor_safety.

## 2. Verdict
**Approve-with-changes.** Deterministic single-pass histogram fold; logic
correct. Findings about the order-invariant guard, debug/release
divergence, panic-discipline comment, and overflow policy.

## 3. Execution status
- `cargo test --lib sample_summary::coverage` → 7 passed (pre-fix).
- `cargo fmt --check` → exit 0. `cargo clippy --lib` → clean on scope.

## 4. Findings

### Major
- **M1 (reliability): re-entered tile double-counts.** An out-of-order
  stream (`tile0, tile1, tile0`) re-opens tile0 and double-bins/double-
  counts `n_tiles`. The walker guarantees coordinate order, but nothing
  asserted it. Fix: debug-assert non-decreasing `(chrom_id, pos)` +
  regression test (out-of-order panics in debug). *(Applied.)*
- **M2 (defaults): `new()` panic-in-debug / silent-clamp-in-release.**
  Release silently substituted a scheme the caller never chose (one B1's
  own `validate` would reject). Fix: `assert!` in both profiles — an
  invalid scheme is a programmer error (CLI validates first). *(Applied.)*
- **M3 (errors): `expect()` in `observe` lacked the `// PANIC-FREE:`
  comment** the crate requires on safe `expect`s. *(Applied — comment
  added explaining `self.open` is always `Some` after the match.)*

### Minor
- **Mi1 (errors): `depth_sum += u64::from(depth)` had no explicit overflow
  policy.** *(Applied — `saturating_add`.)*
- **Mi2 (reliability): dead `mean_depth <= 0.0` branch in `cell_index`**
  (`covered > 0` ⇒ `mean_depth >= 0`; the general formula already maps 0
  to bin 0). *(Applied — removed; doc notes the saturating cast.)*
- **Mi3 (smells): `CoverageBinScheme` duplicates the four bin-scheme
  fields of `CoverageByGcHistogram`.** **Deferred** — the two serve
  different layers (input config vs output data); unifying via
  `#[serde(flatten)]` would change B1's wire shape for marginal benefit.

### Clean categories
idiomatic, naming, refactor_safety — no findings.

## 8. Missing tests added
- `out_of_order_observe_panics_in_debug` (M1)
- `depth_exactly_on_top_edge_is_overflow` (boundary)
- `pos_zero_saturates_to_tile_zero` (boundary)
- `zero_gc_bins_scheme_panics` (M2)

## 9. What's good
- Single running per-tile accumulator + one count matrix: no per-position
  allocation, O(1) per position.
- `N` exclusion and the skipped-tile path match the architecture exactly,
  with tests pinning each.
- The finished histogram is round-tripped through B1's `validate()` in a
  test, tying the producer to the data-model invariants.

## 10. Commands to re-verify
- `./scripts/dev.sh cargo test --lib sample_summary::coverage`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings -A clippy::doc_lazy_continuation`
