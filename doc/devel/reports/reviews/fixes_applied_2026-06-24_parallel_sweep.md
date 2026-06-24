# Fix Application Report: ssr_call_parallel_sweep_2026-06-24.md

**Date:** 2026-06-24
**Source review:** `doc/devel/reports/reviews/ssr_call_parallel_sweep_2026-06-24.md`
**Source state reviewed against:** commit `8d7e4d2`, branch `ssr-cohort`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 2 (Mi1, Mi2) · Nits: 2

### Outcome totals
- Applied: 2 (Mi1 default chunk + test; Mi2 arch doc note) · No action: 2 Nits (optional)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, **1280 passed, 0 failed, 2 ignored**
- `cargo doc --no-deps` / `cargo audit` → not run (no public-doc-link or dependency change)
- Performance check → not applicable (the sweep is not covered by a `benches/` harness)

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | `queue_depth = 0` serializes the sweep | Apply | Applied | `driver.rs` | Pass |
| Mi2 | Minor | Document chunked realization vs the channel pipeline | Apply | Applied | `ssr_call_driver.md` | Pass |
| Nit-a | Nit | Per-chunk `Vec<Option<String>>` allocation | — | No action | None | N/A |
| Nit-b | Nit | `write_genotyped_chunk` 8 args | — | No action | None | N/A |

## 3. Questions asked and answers
- **OQ1 (realization vs arch topology):** resolved by documenting the chunk-parallel `par_iter` realization as the settled J design (no `seq`-reorder needed; the channel pipeline is a measure-first follow-up).

## 4. Per-finding log

### Mi1 — default chunk when `queue_depth == 0`
- **Final status:** Applied. Added `DEFAULT_SWEEP_CHUNK = 1024` and `chunk_size = if queue_depth == 0 { DEFAULT } else { queue_depth }`, so an unset queue depth no longer collapses the sweep to one-locus (≈ serial) chunks; an explicit value is still honoured. Updated the stale `threads`/`queue_depth` field docs.
- **Files changed:** `src/ssr/cohort/driver.rs`. **Tests added:** `sweep_output_is_independent_of_chunk_size` (queue_depth 0 vs 1 → identical VCF over 3 variant loci).
- **Validation:** `cargo test --lib ssr::cohort::driver` → 12 passed; full lib 1280 passed.

### Mi2 — document the realization
- **Final status:** Applied. Added a "J realization (settled 2026-06-24)" callout to `ssr_call_driver.md` §4: the topology is realized as chunk-parallel `par_iter` (order-preserving ⇒ no `seq`-reorder), with the trade-off (no merger-read/genotype overlap) and the fully-overlapping channel pipeline noted as a measure-first follow-up.
- **Files changed:** `doc/devel/architecture/ssr_call_driver.md`.

### Nits — no action
- **Nit-a:** the per-chunk `Vec<Option<String>>` + `String`s allocation is fine at chunk granularity; a reusable buffer is a later micro-opt if profiling flags it.
- **Nit-b:** 8 args on `write_genotyped_chunk` under `#[allow]`; a `SweepCtx` borrow-struct is premature.

## 5. Deferred findings to carry forward
- None (Nits optional; the fully-overlapping channel pipeline is recorded as a measure-first follow-up in the arch doc).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — the sweep is not reachable from a `benches/` harness. The parallel speedup is the milestone's purpose; an end-to-end throughput measurement belongs to real-data calibration.

## 10. Commands run
- `cargo fmt --check`
- `cargo clippy --lib --all-features -- -D warnings`
- `cargo test --lib`

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --lib --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 1280 passed / 0 failed / 2 ignored

## 12. Notes
- With Step J, the `ssr-call` plan (G → H → I → J) is complete: a correct, accurate, memory-bounded, thread-scaling cohort SSR caller. Remaining work is real-data calibration (a separate effort) and the recorded measure-first follow-ups.
