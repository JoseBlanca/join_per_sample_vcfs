# DUST worker pool — parallel DUST-ahead (2026-05-31)

Implementation report for
[cohort_dust_worker_pool.md](../../../doc/devel/implementation_plans/cohort_dust_worker_pool.md):
replace the single DUST-ahead thread (streaming-columnar produce Stage 3)
with a **bounded, ordered worker pool** that DUSTs the independent covered
intervals in parallel. Profiling had shown the single DUST thread was the
wall floor — sdust over the whole genome is ~10 s of single-threaded work,
and the producer waited on it.

## Plan (as executed)

Sub-staged per the incremental discipline:

- **Stage 1** — the coordinator: `IntervalJob` + `flatten_to_jobs`, and
  `DustAheadPool` (bounded-ordered state, worker loop, `recv_next`,
  shutdown). Unit-tested in isolation with a *stub* compute closure (no
  FASTA / no sdust).
- **Stage 2** — wire into `BlockIterator`: `spawn_dust_pool` builds the
  production compute closure; the iterator holds a `DustAheadPool` instead
  of a `Receiver`/`JoinHandle`; `enter_current_interval` calls
  `pool.recv_next()`; the pool joins its workers on drop. End-to-end
  byte-identity.

## Domain intent

DUST depends only on the reference + the region, and covered intervals are
independent by construction (no variant group reaches across the gap
between them — the property byte-identity already relies on). So the masks
for different intervals compute concurrently with the only coordination
being (a) re-imposing the producer's genomic order on out-of-order
completions and (b) bounding how far the pool runs ahead. The per-interval
math (`dust_mask_for_interval`) and the producer-side slice
(`slice_interval_mask`) are unchanged, so the masks — and the VCF — are
**byte-identical** to the Stage 3 branch.

## Design

- **Jobs.** `flatten_to_jobs` lays the per-chromosome plans into one
  delivery-ordered list (chromosome-major, interval-order within) —
  exactly the order the producer enters intervals. Job index `k` is the
  k-th interval entered.
- **Bounded-ordered coordinator** (`DustAheadPool`, one `Mutex` + two
  `Condvar`s). Invariants under the lock:
  `next_deliver <= next_dispatch <= next_deliver + lookahead`, and `done`
  only holds keys in `[next_deliver, next_dispatch)`.
  - *Workers* claim the lowest undispatched job (so they always work on
    what the producer needs soonest), but only within `lookahead` of the
    delivery frontier — the back-pressure that bounds in-flight work and
    resident masks. Compute happens outside the lock; the result is
    stashed and `ready` is notified.
  - *`recv_next`* hands the producer `done[next_deliver]` in order,
    blocking on `ready` until that specific interval completes; advancing
    `next_deliver` frees a look-ahead slot (`not_full` notified). Returns
    `None` once all delivered.
- **Compute via trait object.** `type DustComputeFn = dyn Fn(&IntervalJob)
  -> Result<IntervalDustMask, ChunkDriverError> + Send + Sync`, shared
  across workers behind an `Arc`. Production opens a per-call
  `ManualEvictChromRefFetcher` and runs `dust_mask_for_interval`; tests
  pass a stub. This is what made the coordinator unit-testable without a
  FASTA.
- **Worker count / look-ahead.** `n_workers =
  current_num_threads().clamp(1, n_jobs)`; `lookahead = n_workers × 2`.
  Dedicated OS threads (not rayon), like the consume workers. A 1-worker
  pool is exactly the Stage 3 behaviour — a strict generalisation.
- **Shutdown.** `DustAheadPool` joins its workers on `Drop` (set `abort`,
  notify both condvars, join); a worker mid-DUST finishes its current
  interval first. `BlockIterator`'s own manual `Drop` (Stage 3) is gone —
  dropping the `Option<DustAheadPool>` field handles it.

## Assumptions / decisions

- **One job per covered interval** (not per sub-span). Simple and already
  collapses the floor toward the largest single interval's DUST time. If a
  giant chromosome arm later dominates, the documented follow-up is
  sub-span jobs + intra-interval coalescing at the orderer. Deferred.
- **Lockstep preserved.** Jobs are flattened from the same `chrom_plans`
  the producer walks, so the k-th `recv_next` is the k-th interval entered
  — the same `(chrom_id, interval)` `debug_assert` from Stage 3 still
  guards it.
- **Byte-identity is by construction.** Parallelism changes only *when*
  masks are computed, never their content or delivery order.

## Changes made

All in [driver.rs](../../../src/var_calling/cohort_block/driver.rs):

- New `IntervalJob`, `flatten_to_jobs`, `DustComputeFn`, `DustPoolState`,
  `DustPoolShared`, `DustAheadPool` (+ `dust_pool_worker`, `spawn_dust_pool`),
  replacing `spawn_dust_ahead`. `dust_mask_for_interval`,
  `slice_interval_mask`, `IntervalDustMask` (now `#[derive(Debug)]`),
  `build_dust_plans`, `DustChromPlan` reused unchanged.
- `BlockIterator`: `dust_rx`/`dust_handle` → `dust_pool:
  Option<DustAheadPool>`; `new` spawns the pool; `enter_current_interval`
  uses `pool.recv_next()`; the manual `Drop` impl removed (pool's `Drop`
  joins).
- Const `DUST_AHEAD_QUEUE_CAP` → `DUST_POOL_LOOKAHEAD_PER_WORKER`.
- Imports: `Arc`, `HashMap`.

## Tests added

Coordinator unit tests (stub compute, no FASTA), in `driver.rs`:

- `dust_pool_delivers_in_order_despite_reverse_completion` — a gate forces
  strict reverse completion; masks are still delivered in job order.
- `dust_pool_respects_lookahead_backpressure` — without consuming, at most
  `lookahead` jobs dispatch; consuming one frees exactly one slot.
- `dust_pool_surfaces_compute_error_in_order` — a compute error surfaces at
  its own seq without disturbing neighbours.
- `dust_pool_delivery_is_parallelism_invariant` — 1-worker and 8-worker
  pools deliver the identical ordered sequence.
- `dust_pool_shuts_down_cleanly_mid_flight` — dropping mid-flight joins
  without hanging.

The Stage 3 `dust_mask_for_interval_*` and `slice_interval_mask_*` tests
are unchanged (per-interval math and producer slice untouched).

## Validation results

- `cargo test --lib`: **1061 passed**, 0 failed, 1 ignored (1056 + 5 pool).
  `cargo clippy --lib`: clean (one pre-existing, unrelated
  `needless_range_loop` in `dust_filter.rs`). `cargo fmt --check`: clean on
  `driver.rs`.
- **Byte-identical** (drop `^##`, then md5) to the Stage 3 baselines,
  **serial and 8 threads**, N=8/26:

  | N  | md5 (serial = 8 threads) | baseline |
  |----|--------------------------|----------|
  | 8  | `8f117ac7…` | `8f117ac7fef9f83d240644b1d0def3c8` ✓ |
  | 26 | `a3164c38…` | `a3164c386c1b21cf94c6adc5cf4c60d3` ✓ |

- **Wall (real tomato cohort, 8 threads):**

  | N  | Stage 3 (1 DUST thread) | **Stage 4 (pool)** | DUST-off floor | DUST's added cost |
  |----|------------------------:|-------------------:|---------------:|------------------:|
  | 8  | 8.71 s | **4.02 s** | 1.83 s | 6.9 → **2.2 s** |
  | 26 | 9.45 s | **7.60 s** | 5.52 s | 4.0 → **2.1 s** |

- **vs `main`** (`perf_main_vs_branch.py`, 8 threads): the branch now
  **beats `main` on wall *and* RSS at both sizes**:

  | N  | main wall | branch wall | speedup | main RSS | branch RSS |
  |----|----------:|------------:|--------:|---------:|-----------:|
  | 8  | 6.0 s | **4.1 s** | **1.47×** | 142 MB | **79 MB** |
  | 26 | 8.9 s | **7.8 s** | **1.14×** | 403 MB | **291 MB** |

  (The branch-vs-`main` md5 difference is the known, pre-existing
  filter-order divergence — orthogonal; the byte-identity gate above is
  against the branch's own output.)

## Tradeoffs and follow-ups

- **Residual DUST cost (~2 s).** With one job per interval, the floor is
  the largest single interval's DUST time plus the cold first-interval
  stall. The escape hatch — if a giant arm dominates — is **sub-span
  jobs** with intra-interval coalescing at the orderer; deferred until
  measured necessary.
- **Oversubscription** at high N (DUST workers + consume workers + decode
  rayon) is bounded by the look-ahead gate (excess workers idle once
  ahead); the static `n_workers = threads` default could be tuned per N if
  a regression appears, but at N=26 it already beats `main`.
