# Cohort `var-calling` — DUST worker pool (parallel DUST-ahead)

Implementation plan for parallelising the **DUST-ahead** stage of the
`cohort_block` producer. This is the direct follow-up to
[cohort_produce_streaming_columnar.md](cohort_produce_streaming_columnar.md)
Stage 3, which moved the sdust low-complexity mask off the producer's
critical path onto a single background thread. Profiling shows that one
thread is now the wall floor: DUST is ~10 s of *single-threaded* work over
the whole genome, and the producer waits on it.

## Motivating measurement

Real tomato cohort, host release binary, `var-calling` PSP→VCF. Ablation
with vs without the complexity filter (`--no-complexity-filter`):

| N  | threads | DUST on | DUST off | DUST's share of wall |
|----|--------:|--------:|---------:|---------------------:|
| 8  | 8 | 8.71 s | **1.83 s** | **~79 % (6.9 s)** |
| 26 | 8 | 9.45 s | 5.44 s | ~42 % (4.0 s) |
| 8  | 1 | 9.95 s | 5.03 s | ~49 % |
| 26 | 1 | 27.1 s | 26.6 s | ~2 % (hidden behind serial EM) |

Thread-scaling confirms a hard serial floor — N=8 wall flatlines at
~8.7 s from 2 threads on (1.14× from 1→8 threads), and **N=26 also
plateaus at 4 threads** (9.98→9.87 s for 4→8). A samply/`sample` profile
of the floor shows the producer's main thread blocked in
`enter_current_interval → recv()` (waiting on the DUST queue) while the
DUST-ahead thread is busy in `dust_filter::find_perfect` / `shift_window`.

**Conclusion:** the remaining wall floor is the single-threaded
DUST-ahead thread. sdust over the whole tomato genome is ~9.7 s, constant
in N (it scales with analysed span, not sample count). At low N there is
almost nothing else to hide it behind (1.83 s of real produce+consume);
at high N the parallel consume overlaps more of it but it still adds ~4 s.
`main` does not have this floor because it runs DUST *parallel across
chromosomes* via its per-chromosome `rayon::par_iter`; our DUST-ahead runs
on one thread.

## Why this parallelises cleanly

The DUST mask depends only on the reference + the region. The producer
already processes the genome as **covered intervals** — maximal
data-bearing segments separated by gaps wider than any variant group's
reach
([`covered_intervals_for_chrom`](../../../src/var_calling/cohort_block/driver.rs)).
Those intervals are *independent by construction* (the same property
byte-identity already relies on). So the masks for different covered
intervals can be computed **concurrently** with no coordination beyond
re-imposing genomic order on the results — exactly the structure `main`
exploits across chromosomes, but at the finer interval granularity.

The per-interval computation
([`dust_mask_for_interval`](../../../src/var_calling/cohort_block/driver.rs))
is unchanged and reused verbatim; only *who runs it and when* changes. The
masks are bit-for-bit the same, so the change is **byte-identical** to the
current branch (Stage 3) — the same gate.

## Goal / non-goals

**Goal.** Replace the single DUST-ahead thread with a bounded, ordered
**worker pool** that DUSTs covered intervals in parallel and feeds them to
the producer in genomic order. Target: at N=8/8-threads, pull wall from
~8.7 s toward the ~1.8 s DUST-off floor (plus the initial cold-interval
stall); at N=26/8-threads, recover the ~4 s DUST adds. Byte-identical at
N=8/26, serial + 8 threads.

**Non-goals.**
- The fold / decode / partition path (measured cheap: 1.83 s at N=8) —
  no multi-producer over the fold.
- The consume half (already parallel across blocks).
- The known, pre-existing filter-order divergence vs `main`.
- Splitting a single covered interval across workers — see *Granularity*
  below; deferred unless the largest-interval floor proves to matter.

## Design

### Jobs — covered intervals in genomic order

Flatten the existing per-chromosome plans
([`DustChromPlan`](../../../src/var_calling/cohort_block/driver.rs)) into a
single seq-numbered job list, chromosome-major then interval-order within
a chromosome — i.e. **exactly the order the producer enters intervals**
(`advance_chrom` walks the plans in order; `interval_idx` 0..n within
each). So `seq == N` is the Nth interval the producer will enter.

```
struct IntervalJob {
    seq: usize,            // genomic order == producer entry order
    chrom_id: u32,
    name: String,          // contig name, for the per-worker fetcher
    length: u32,
    interval: Range<u32>,  // 1-based half-open covered interval
}
```

### Bounded, ordered worker pool

A `DustAheadPool` owns the worker threads and the shared coordination
state. The producer-facing API stays a single pull, so
`BlockIterator::enter_current_interval` changes only from
`rx.recv()` to `pool.recv_next()` (same `Result` semantics).

```
pool.recv_next() -> Option<Result<IntervalDustMask, ChunkDriverError>>
```

Returns the mask for the next interval **in seq order**, blocking until
that specific interval is done; `None` once all intervals delivered.

Shared state (one `Mutex` + two `Condvar`s):

```
next_dispatch: usize                 // next job a worker will claim
next_deliver:  usize                 // next seq the producer wants
done: HashMap<usize, Result<IntervalDustMask, ChunkDriverError>>  // completed, undelivered
n_jobs: usize
abort: bool
```

- **Worker loop.** Claim the lowest undispatched job, but only if it is
  within `lookahead` of the delivery frontier (back-pressure):
  ```
  lock
  while !abort && next_dispatch < n_jobs && next_dispatch - next_deliver >= lookahead:
      wait(not_full)
  if abort || next_dispatch >= n_jobs: break
  job = jobs[next_dispatch]; next_dispatch += 1
  unlock
  result = dust_mask_for_interval(open fetcher for job.name, job.interval, …)   // Err captured, not panicked
  lock; done.insert(job.seq, result); notify_all(ready)
  ```
  In-order dispatch means workers always pick up the intervals the
  producer will need *soonest*, so the producer rarely waits after the
  first interval.

- **`recv_next` (producer side).**
  ```
  lock
  loop:
      if let Some(r) = done.remove(&next_deliver):
          next_deliver += 1; notify_all(not_full); return Some(r)   // freed a lookahead slot
      if next_deliver >= n_jobs: return None
      wait(ready)
  ```

- **Per-worker fetcher.** A worker may bounce between chromosomes, so it
  opens a `ManualEvictChromRefFetcher::for_contig(fasta, &job.name)` per
  job (optionally caching the last `(name, fetcher)` to skip a reopen when
  consecutive jobs share a contig). Covered intervals total only ~tens
  across the genome, so per-job opens are negligible.

- **Shutdown.** `Drop for BlockIterator` (today: drop receiver, join one
  thread) becomes: set `abort = true`, notify both condvars, join all
  workers. A worker mid-DUST finishes its current interval then exits —
  bounded by one interval's DUST time (same property Stage 3 has).

### Byte-identity

The delivered mask for interval `k` is `dust_mask_for_interval(k)`,
identical to Stage 3's. Ordered delivery (`next_deliver`) hands the
producer interval `0, 1, 2, …` exactly as the single-thread version did.
Parallelism changes only *timing*, never content or order. So the VCF is
bit-identical to the current branch; the gate is the existing header-
stripped md5 against the Stage 3 baselines.

### Memory bound

Resident DUST state ≤ `lookahead` interval masks (tiny — `Vec<Range>`) +
`n_workers` transient sub-span reference buffers (`dust_mask_for_interval`
already caps each at ~one `DUST_AHEAD_SUBSPAN` + halo ≈ 1 MB via
`evict_before`). With `lookahead ≈ 2·n_workers` and `n_workers ≤ 8` that
is a few MB on top of the current ~294 MB — flat in N. The `lookahead`
gate is what keeps a fast pool from DUSTing the whole genome ahead of a
slow producer.

### Worker count

DUST is the floor mainly where consume is cheap and cores are idle (low
N), so the pool wants the free cores there; at high N the consume workers
want them. Start with `n_dust_workers = rayon::current_num_threads()`
clamped to `[1, n_jobs]`, and tune by measurement — mild oversubscription
at high N is acceptable (DUST workers block on the lookahead gate once
ahead). Like the consume workers, these are **dedicated OS threads, not
rayon tasks**, so they don't park on the pool the producer's per-sample
decode uses (the starvation caveat documented in `drive_blocks_parallel`).
A pool of 1 worker reproduces Stage 3 — this is a strict generalisation.

### Granularity (the one real limit)

One job = one covered interval, so a single interval is an indivisible
unit: the floor cannot drop below the *largest* interval's DUST time.
With ~12 chromosomes and coverage merged into a handful of intervals each,
the largest interval is well under the per-genome total, so an interval-
granular pool already collapses the ~9.7 s floor toward
`largest_interval + tail`. **If** profiling after this lands shows one
giant arm dominating, the follow-up is to emit **sub-span jobs** (each
`DUST_AHEAD_SUBSPAN` slice is a job) and have the orderer gather + coalesce
an interval's sub-span masks before delivering it — finest granularity,
full load-balancing. Deferred until measured necessary; noted here so it
isn't rediscovered.

## Files

- [driver.rs](../../../src/var_calling/cohort_block/driver.rs):
  `IntervalJob` + `flatten_to_jobs`; the `DustAheadPool` (shared state,
  worker loop, `recv_next`, shutdown) replacing `spawn_dust_ahead`;
  `BlockIterator` holds a `DustAheadPool` instead of
  `dust_rx`/`dust_handle`; `enter_current_interval` calls
  `pool.recv_next()`; `Drop` shuts the pool down. `dust_mask_for_interval`
  and `IntervalDustMask` are reused unchanged.

## Staging

Single feature, sub-staged for gateable steps (types → impl → wire →
verify), per the incremental-step discipline:

1. **Pool types + coordination** — `IntervalJob`, `flatten_to_jobs`,
   `DustAheadPool` with the bounded-ordered state and worker loop. Unit-
   test the coordinator in isolation (ordered delivery, back-pressure,
   error propagation, shutdown) with a stub compute closure — no FASTA.
2. **Wire into `BlockIterator`** — swap `spawn_dust_ahead` for the pool,
   `recv()` for `recv_next()`, update `Drop`. End-to-end byte-identity.
3. **Tune `n_dust_workers` / `lookahead`** by measurement; pick defaults.

## Gate

- **Byte-identical** (header-stripped md5) at N=8/26, serial + 8 threads,
  against the Stage 3 baselines
  (`8f117ac7fef9f83d240644b1d0def3c8`, `a3164c386c1b21cf94c6adc5cf4c60d3`)
  — DUST content/order is unchanged.
- **Wall:** N=8/8-threads from ~8.7 s toward the ~1.8 s DUST-off floor
  (plus the cold first-interval stall); N=26/8-threads recovers most of
  the ~4 s DUST adds. Re-run the thread-scaling sweep and the
  ablation — the N=8 floor should no longer track DUST-on.
- **Memory:** peak RSS flat in N, within a few MB of the current ~294 MB
  (N=26).
- Tests + clippy clean.

## Tests

- **Coordinator unit tests** (stub compute, no FASTA): (a) ordered
  delivery despite out-of-order completion (workers finish high seqs
  first); (b) back-pressure — with `lookahead = L`, `next_dispatch` never
  exceeds `next_deliver + L`; (c) an error from one job surfaces at its
  seq in order and does not corrupt neighbours; (d) shutdown mid-flight
  joins cleanly.
- **Reuse** the Stage 3 `dust_mask_for_interval_matches_whole_interval_scan`
  and `slice_interval_mask_*` tests unchanged (the per-interval math and
  the producer slice are untouched).
- A 1-worker pool equals the multi-worker pool on the same job set
  (parallelism-invariance of the delivered sequence).

## Risks

- **Concurrency correctness** is the substantive risk (Mutex + two
  Condvars, two frontiers). Mitigate: design the invariants on paper
  (`next_deliver ≤ next_dispatch ≤ next_deliver + lookahead`; `done` only
  holds `[next_deliver, next_dispatch)`), unit-test the coordinator with a
  stub before wiring real DUST, and keep `recv_next`'s blocking semantics
  identical to the old `recv` so the producer loop is unchanged.
- **Oversubscription** at high N (DUST workers + consume workers + decode
  rayon). Measured-tuned worker count; the lookahead gate idles excess
  workers once ahead.
- **Largest-interval floor** — accepted for v1; sub-span jobs are the
  documented escape hatch.
