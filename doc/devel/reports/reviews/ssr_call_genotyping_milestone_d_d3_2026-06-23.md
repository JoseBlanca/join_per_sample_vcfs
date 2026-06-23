# Code Review: ssr-call genotyping+pre-pass — Milestone D D3 (clustering, M3)
**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator, focused inline pass)
**Scope:** Milestone D D3 (commit `620d091`, `sample_groups.rs` + the `prepass.rs` extension)
**Status:** Approve-with-changes

---

## 1. Scope
- **In-scope:** [sample_groups.rs](../../../../src/ssr/cohort/sample_groups.rs) (new),
  the `prepass.rs` additions (`slip_by_sample_period`, exposed helpers,
  `run_prepass_stats`).
- **Categories:** reliability, errors, naming, idiomatic, refactor_safety, smells, extras
  (numerical correctness / determinism). Skipped unsafe_concurrency, tooling.

## 2. Verdict
**Approve-with-changes.** M3 holds (per-group divergent shapes recovered, single protocol
collapses), the clustering is deterministic (union toward the smaller root + catalog-index
labelling), and the per-group fit reuses the verified primitives. Findings are a magic
number, a documented-but-worth-flagging chaining limitation, and a determinism test.

## 3. Execution status
- `cargo fmt --check` → pass · `cargo clippy --all-targets --all-features -- -D warnings`
  → pass (fixed a `needless_range_loop` during implementation) · `cargo test --all-features`
  → **1241 lib pass**.
- Needs-verification findings: 0.

## 4. Top 3 priorities
1. **Mi1** — `scaled_distance` hardcodes a `100.0` reference depth (and a `1.0` cap);
   move it into `ClusterCfg` so the precision-weighting has no hidden default.
2. **Mi2** — single-linkage union-find chains a smooth `(ε, level)` continuum into one
   blob; the BIC split test that breaks an over-merged group is deferred — re-flag at the
   `cluster` site, not only the module header.
3. **MT-1** — add a `group_samples` determinism test (two runs → equal `GroupedParams`).

## 5. Findings

### Minor

- `src/ssr/cohort/sample_groups.rs` (`scaled_distance`) — **[Minor]** hidden reference-depth constant
- **Confidence:** High
- **Problem:** `let precision = (depth.min() as f64).sqrt() / 100.0_f64.sqrt();` bakes a
  reference depth of `100` and a `.min(1.0)` cap into the function. These are tuning
  knobs and belong in `ClusterCfg` (the defaults checklist: no behaviorally-significant
  hidden defaults).
- **Why it matters:** F2 calibration can't reach them; a reader can't see them at the
  call site.
- **Suggested fix:** add `precision_reference_depth: f64` to `ClusterCfg` (default 100.0)
  and use it; keep the `1.0` cap with a comment that above the reference depth the
  precision weighting saturates.

- `src/ssr/cohort/sample_groups.rs` (`cluster`) — **[Minor]** single-linkage can chain a continuum
- **Confidence:** High
- **Problem:** Union-find with a fixed threshold is single-linkage: a smooth gradient of
  `(ε, level)` values whose neighbours are each within the threshold merges into one
  group, even when the extremes are far apart (the spec's caveat). The BIC split test
  that breaks an over-merged group is deferred.
- **Why it matters:** A deep continuum cohort would under-segment. The module header notes
  the deferral, but the `cluster` function itself should warn.
- **Suggested fix:** document the single-linkage chaining limitation on `cluster`,
  pointing at the deferred BIC split test (F).

### Nits
- `scaled_distance` — the precision direction (lower depth ⇒ smaller distance ⇒ merges
  more readily) is correct but subtle; a one-line "why" comment helps.

## 6. Out of scope observations
- `benches/psp_writer_perf.rs:386` — pre-existing bench panic, unchanged.

## 7. Missing tests to add now
- `group_samples_is_deterministic` — **input:** the M3 cohort, grouped twice. **Bug it
  catches:** an order-dependence in the union-find or group-id assignment (e.g. iterating
  a `HashMap` for the labelling instead of the sorted features). **Body:** assert the two
  `GroupedParams` are equal.

## 8. What's good
- The clustering is deterministic by construction (sorted features, union toward the
  smaller root, group ids by ascending root) — no random init, reproducible across runs
  ([sample_groups.rs](../../../../src/ssr/cohort/sample_groups.rs)).
- D3 reuses `refine_theta_locus` for the shape shrinkage and `fit_level` / `sample_eps`
  for the per-group fit — one definition of each, shared with D2.
- M3 is a real recovery test (divergent decays), not just "it runs".

## 9. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --all-targets --all-features -- -D warnings` ·
  `cargo test --all-features`

### Author response convention
Address each finding by ID (Mi1, Mi2, Nit, MT-1).
