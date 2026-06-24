# Code Review: ssr-call chunk-parallel sweep (Step J)

**Date:** 2026-06-24
**Reviewer:** rust-code-review skill (orchestrator; inline synthesis)
**Scope:** commit `8d7e4d2` on branch `ssr-cohort` — parallelize the Pass-2 genotyping sweep
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** PR diff, commit `8d7e4d2` (Step J).
- **In-scope files:** [driver.rs](../../../../src/ssr/cohort/driver.rs) — the chunk-parallel Pass-2 sweep in `run`, `write_genotyped_chunk`, the `rayon::prelude` import, the `queue_depth` doc, the module doc, 1 new test.
- **Out of scope:** `genotype_locus` / the EM (reviewed in H4 / I1 / I2); the burn-in (H4).
- **Categories considered:** reliability, errors, **unsafe_concurrency** (rayon `par_iter` + the shared pool), idiomatic, defaults, smells, refactor_safety, extras (determinism, hot path).

## 2. Verdict

**Approve-with-changes.** The chunk-parallel sweep is correct, bounded, ordered, and byte-identical across thread counts (a new multi-chunk test plus the existing 1-locus determinism test confirm it). The concurrency is safe by construction (pure per-locus function over `Sync` shared refs; no shared mutable state). No Blocker/Major. Minors: a `queue_depth = 0` would serialize the sweep, and the chunked realization is a deliberate simplification of the arch's overlapping producer/worker/writer pipeline — both worth pinning down.

## 3. Execution status

- `cargo fmt --check` clean · `cargo clippy --lib -D warnings` clean · `cargo test --lib` = **1279 passed, 0 failed, 2 ignored** (+1 multi-chunk test).
- Sub-agent fan-out not used (overload); reviewed inline. "Needs verification": 0.

## 4. Open questions and assumptions

1. **Realization vs the arch topology.** The plan (Milestone J) named the reading layer's *producer → bounded-queue → worker pool → writer (reorder by seq)* topology. This commit realizes the same guarantees (parallel, bounded, ordered, byte-identical) more simply via **chunking** + an order-preserving `par_iter`, so no `seq`-reorder is needed. That is a legitimate, lower-machinery choice, but it diverges from the literal topology and does **not** overlap the merger read with genotyping (it reads a chunk, then processes it). Worth recording as the chosen design (see Mi2).

## 5. Top 3 priorities
1. **Mi1** — guard against `queue_depth = 0` (it would make 1-locus chunks ≈ serial).
2. **Mi2** — document the chunked-vs-pipeline decision (and that `seq`-reorder is unnecessary) in the arch/plan + a code note.
3. Nits below.

## 6. Findings

### Minor

**Mi1: [driver.rs](../../../../src/ssr/cohort/driver.rs) — `chunk_size = config.queue_depth.max(1)` serializes when `queue_depth == 0`.**
A `0` (or `1`) queue depth yields one-locus chunks, i.e. no genotyping parallelism — easy to hit if the CLI leaves `queue_depth` unset/zero. Confidence: High.
*Fix:* treat `0` as "use a sensible default" — `let chunk_size = if config.queue_depth == 0 { DEFAULT_SWEEP_CHUNK } else { config.queue_depth };` with a documented `DEFAULT_SWEEP_CHUNK` (e.g. 1024). Keeps an explicit small value honoured (for tests) but prevents an accidental serial sweep.

**Mi2: [driver.rs](../../../../src/ssr/cohort/driver.rs) — record the chunked realization vs the arch's overlapping pipeline.**
The chunk approach reads a full chunk before genotyping it (no producer/worker overlap), where the arch §5/§8 describes a producer thread + bounded channel + worker pool + reorder-by-seq writer. The chunk version is simpler and meets every stated guarantee, but the divergence and its trade-off (a small read-then-process serialization vs. no channel machinery) should be written down so a future reader knows it was a choice. Confidence: Medium (documentation).
*Fix:* a one-paragraph note in the arch §J / plan J that the topology is realized as chunk-parallel `par_iter` (order-preserving ⇒ no `seq`-reorder), with the fully-overlapping channel pipeline as a measure-first follow-up.

### Nits
- `write_genotyped_chunk` allocates a `Vec<Option<String>>` (and the `String`s) per chunk; fine at chunk granularity, but a reusable line buffer is a later micro-opt if profiling flags it.
- 8 parameters on `write_genotyped_chunk` (`#[allow]`); acceptable, or fold the `frozen + cfgs` into a small `SweepCtx<'_>` borrow-struct shared with `genotype_locus` if either grows.

## 7. Out of scope observations
None.

## 8. Missing tests to add now
None required — `chunk_parallel_sweep_orders_records_and_is_deterministic` covers multi-chunk ordering + cross-thread byte-identity, and `run_is_byte_identical_across_thread_counts` covers the single-chunk parallel path. (A `queue_depth = 0` test should accompany the Mi1 fix.)

## 9. What's good
- The whole sweep runs inside the existing `pool.install` scope, so `par_iter` transparently uses the `--threads` pool and the I/O stays on the calling thread — no second pool, no channel plumbing.
- Correctness is structural: `genotype_locus` is a pure function of `(locus, frozen params)` over `Sync` shared references, and `par_iter().collect::<Vec<_>>()` is index-ordered, so byte-identity falls out without any per-thread reasoning — exactly what the two determinism tests assert.
- Memory stays bounded by the chunk (peak = `chunk_size` loci + the merger's lockstep blocks), preserving the cohort-scaling property end to end.

## 10. Commands to re-verify
- `cargo fmt --check` · `cargo clippy --lib --all-features -- -D warnings` · `cargo test --lib ssr::cohort::driver`

### Author response convention
Address Mi1, Mi2, Nits by id; answer open question 1 first.
