# Concurrency & contention checklist

**Purpose.** Surface lock-contention bottlenecks, atomic-ordering overhead, and parallelism that costs more than it saves.

**Triggers.** Code uses `Arc`, `Mutex`, `RwLock`, atomics, channels (`mpsc`, `crossbeam-channel`, `tokio::sync`), `rayon::par_iter` / `rayon::scope`, `tokio::spawn`, `std::thread::spawn`, `async fn` with `.await` while holding shared state.

**Skip when.** Single-threaded code with no shared state. If dispatched anyway, write `No findings.`

**Note on overlap with the correctness review.** Some patterns in this category (lock held across `.await`, weakened atomic ordering, manual `unsafe impl Send/Sync`) are correctness defects whose performance signature is a side effect. File the finding here with a measurement plan, but explicitly route it through `unsafe_concurrency` review for the correctness call.

## Rules

- **Lock acquisition does not belong inside a hot loop.** Acquire once, batch the work, drop the guard. A common shape of this bug is a comparator or filter closure that locks shared state on every call from inside a sort or scan; the fix moves the acquisition out of the per-call path. Filed for any `lock().unwrap()` (or equivalent) inside a tight inner loop.
- **Lock granularity follows access patterns, not aesthetics.** A single `Mutex<HashMap<K, V>>` accessed concurrently is a serial bottleneck. Replacements, in increasing order of complexity:
  - `RwLock<HashMap<K, V>>` if reads dominate writes ≫ 10:1 (writer starvation is real; read-heavy is the precondition).
  - `dashmap::DashMap` for general per-bucket sharded access (sharded `RwLock`s; reduces but does not eliminate contention; writes still block readers in the same shard).
  - Per-shard `Vec<Mutex<HashMap<K, V>>>` where shard count and hash function are tuned to the workload. Highest manual cost, fewest surprises.
- **`std::sync::Mutex` vs `parking_lot::Mutex` is a contention-pattern decision, not a default-faster decision.** Shape of the trade-off:
  - Short critical sections, low contention: `std::sync::Mutex` typically wins on throughput. Default to it.
  - Long holds, heavy contention: `std::sync::Mutex`'s "barging" lets active threads keep grabbing the lock ahead of queued waiters, which can starve some threads under sustained pressure. `parking_lot::Mutex` enforces eventual fairness via a periodic handoff timer.
  - Bursty workloads, or risk of one thread monopolizing the lock: `parking_lot` is the safer default — its fairness mechanism bounds worst-case wait.
  - Per-thread fairness or P99 latency targets matter: `parking_lot`.
  - Switch only with a benchmark that shows the contention pattern; file the experiment, not the swap.
- **`RwLock` is not a free upgrade from `Mutex`.** It has more overhead in the uncontended case and is vulnerable to writer starvation. Justified only when reads vastly outnumber writes and a measurement shows `Mutex` blocks readers.
- **Memory orderings are justified by the synchronization relationship.** Default to `SeqCst` and only weaken with a comment naming what loads pair with what stores. `Relaxed` is correct for independent counters and free-running statistics; `Relaxed` for synchronization without a separate fence is a Hot-path correctness finding (route through correctness review).
- **Holding a `std::sync::MutexGuard` across `.await` is a bug, not a perf issue.** Either drop the guard before awaiting (scope it in `{}`) or switch to `tokio::sync::Mutex`. Mark the finding and route it through unsafe/concurrency review — performance is downstream of correctness here.
- **`tokio::sync::Mutex` is slower than `std::sync::Mutex`.** Use it only when the lock is held across `.await`. For brief critical sections in async code, `std::sync::Mutex` scoped to a sync block is faster *and* correct.
- **CPU-bound work on a tokio worker stalls every other task on that thread.** Use `tokio::task::spawn_blocking` for sync CPU work; bound the spawn-blocking pool with a semaphore so unbounded parallelism does not exhaust threads.
- **Channels: bounded by default.** Unbounded channels (`mpsc::unbounded_channel`, `crossbeam_channel::unbounded`) silently buffer until the process OOMs. Filed for any unbounded channel reading from an external or large source.
- **`rayon::par_iter` over short work is slower than the serial version.** Fork-join overhead bites when per-element work is below a few hundred nanoseconds. `par_iter()` over collections of small numeric ops without a chunking strategy or `with_min_len` tuning is a finding. Run the serial-vs-parallel benchmark before merging.
- **Rayon's default splitting assumes roughly uniform per-item cost.** Pre-divided binary-tree ranges leave fast threads idle while one thread finishes the heavy tail when per-item cost varies widely (e.g., variable-coverage pileups, simple vs. complex sites). `with_min_len` controls split *depth*, not runtime balance, so it does not fix this. Symptoms: high wall-clock variance across runs, threads finishing at very different times, `perf` showing futex/sched_yield concentrated near end-of-run. Mitigations to consider: cost-aware partitioning at the producer (split heavy items before handing to Rayon), smaller chunks via `with_min_len`, or a custom adaptive stealing loop. Profile with `perf` before tuning — Rayon overhead is often not the real bottleneck.
- **Atomic counters cause cross-CPU traffic.** A single global `AtomicU64` incremented from many threads becomes the bottleneck on its own cache line. Patterns: per-thread counters summed at end-of-phase, or `crossbeam_utils::CachePadded<AtomicU64>` per shard. Cross-reference: `data_layout.md` on false sharing.
- **`Arc::clone` is one atomic increment, not free.** Gratuitous `Arc::clone` in a hot loop (e.g., capturing into a closure that does not outlive the iteration) shows up. Borrow when possible.
