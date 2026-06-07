# Research: H1 thread oversubscription & collision — measurement, mitigation, rayon/crossbeam best practices

**Date:** 2026-06-07
**For:** perf finding **H1** in
[perf_var_calling_cohort_2026-06-06.md](../reviews/perf_var_calling_cohort_2026-06-06.md)
(the T≥2 wall gap = scheduler oversubscription)
**Method:** deep-research harness (5 search angles, 16 sources fetched, 73 claims
extracted, 25 adversarially verified → 21 confirmed / 4 killed). Sources tagged
`primary` (rayon docs/issues, rust-lang PRs, kernel docs, num_cpus source) vs
`blog`/`forum`/`secondary`.

---

## 0. The one correction the research forced (read this first)

The verified findings are sound, but the agents pattern-matched the profile frame
`bridge_producer_consumer` to rayon's **`par_bridge`** issues (#730, #795). **Our
producer does not use `par_bridge`.** It uses indexed `par_iter_mut()` over the
sample slice ([cohort_integration.rs:745](../../../../src/var_calling/cohort_integration.rs#L745),
[:926](../../../../src/var_calling/cohort_integration.rs#L926)). `bridge_producer_consumer`
is rayon's generic *indexed* parallel-iterator driver — it is the function that
recursively splits the producer and runs the consumer for **`par_iter`/`par_iter_mut`
on slices**, not just `par_bridge`.

So the mechanism behind the 64 %-blocked producer is **not** par_bridge
starvation (an underfed serial iterator spinning workers). It is simpler and
exactly what the perf review said: the producer's **main thread calls
`par_iter_mut().install`-style and blocks on `LockLatch::wait_and_reset`
(`__psynch_cvwait`) waiting for the rayon burst to finish** — while `T` crossbeam
caller threads are simultaneously running the EM. That is `T` rayon workers + `T`
callers + writer + producer-main = up to `2·T + 2` runnable CPU-bound threads on
`T` cores.

This correction *narrows* the fix: the par_bridge-specific remedies (batching the
serial iterator, replacing par_bridge with a manual channel) do **not** apply.
The pool-sizing and pool-unification remedies **do**. Everything below is filtered
through that.

---

## 1. Root cause (confirmed, high confidence)

**Two independent pools each self-size to all detected CPUs with no coordination,
so their counts add to ~2·T runnable threads on T cores.**
[rayon FAQ], [ThreadPoolBuilder docs], [PostHog].

- Rayon's default pool = `RAYON_NUM_THREADS` or the logical-CPU count. Our
  `configure_rayon_pool(args.threads)`
  ([var_calling.rs:238](../../../../src/pop_var_caller/var_calling.rs#L238)) sets it
  to `--threads`.
- The crossbeam caller pool is *separately* sized to `--threads` in
  `pipeline.rs`.
- Neither library knows the other exists. PostHog's production root cause was
  literally this: "both crates spawned as many threads as the number of CPU
  cores… both assume they're the only thing running" → "100 % thread
  oversubscription."
- Matches our benchmark shape exactly: on-par at T=1 (one rayon worker ≈ the main
  thread, no contention), then +6/+12/+15 % at T=2/4/8 as the scheduler
  time-slices ~2·T threads onto T cores.

---

## 2. How to MEASURE it (confirmed)

### Linux (the production container — the place that matters)
- **Involuntary context switches** — `/usr/bin/time -v ./prog` →
  "Involuntary context switches". Oversubscription inflates these; a single-pool
  run at the same T is the baseline. (`perf stat -e context-switches,cs` for the
  same signal with timing.)
- **Run-queue depth** — `perf sched record` / `perf sched latency`, or sampling
  `/proc/loadavg` vs owned-core count. Run-queue length > owned cores while CPU is
  pinned at 100 % is the oversubscription fingerprint.
- **Futex / cv wait time** — `perf record` then look for `futex`/`__psynch`-class
  frames; rising futex wait at higher T with no extra useful work = scheduling
  overhead, not progress.
- **cgroup CFS throttling (distinct failure mode — see §5)** — read
  `cpu.stat`: `nr_throttled` (count of throttle events) and `throttled_time`
  (ns; v2: `throttled_usec`). Rising `nr_throttled` relative to `nr_periods`
  means the bottleneck is *quota throttling*, separate from raw oversubscription —
  and oversubscription **compounds** it. [kernel sched-bwc].

### macOS (the profiling host — where the H1 profile was taken)
- `sample` / Instruments (Time Profiler) / dtrace. The cv-wait/lock-latch frames
  we already captured (`__psynch_cvwait`, `LockLatch::wait_and_reset`,
  `swtch_pri 750`, `__ulock_wait2 309`, `__psynch_mutexwait 156`) **are** the
  oversubscription read-out. macOS has no cgroup throttling, so the macOS profile
  cannot show the container's compounding throttle — confirm that on Linux.

### Reading a CPU profile: oversubscription vs the other scaling losses
The discriminator is **what the threads are doing when blocked, and whether more
threads buys more work**:

| Symptom | Evidence that points to it |
|---|---|
| **Oversubscription** (too many runnable threads) | High *involuntary* context switches + elevated `swtch_pri`/futex/cv-wait that grows with T, while per-thread useful self-time falls; CPU near 100 % but wall worse than fewer threads. Fix: fewer threads → wall improves. |
| **Lock contention** | Self-time concentrated in `__psynch_mutexwait`/`futex` on a *specific* lock; threads serialize on one critical section. (Our H1 has *no* hot-path locks — perf review §7 confirms no Mutex/RwLock/atomic in scope — so this is ruled out.) |
| **False sharing** | `perf c2c` (Linux) shows HITM on cache lines; counters, not wall, are the tell. [easyperf c2c] |
| **Memory-bandwidth saturation** | Scaling flattens with *low* context-switch and *low* lock time; `perf stat` memory-bw counters near peak; adding threads does nothing. |
| **Amdahl serial fraction** | Wall plateaus but CPUs are *idle* (not 100 %); the serial section shows as one thread busy while others wait on an empty channel (`__psynch_cvwait` on `chunk_rx`). This is our L1, *downstream* of H1. |

**Our H1 sits squarely in row 1** (involuntary switches + futex/cv-wait grow with
T, no single hot lock, CPUs busy not idle). L1 (producer serial floor) is row 5
and only becomes visible *after* H1 is fixed.

---

## 3. How to MITIGATE it (the two pools are genuinely both CPU-bound)

Both sides are CPU-bound here (producer = zstd decode + compaction; callers = EM),
so the goal is strictly **rayon_threads + caller_threads ≈ owned cores**, not
2·T. Three options, cheapest first — these mirror the perf review's variants but
are now backed by the sources.

### Option A — cap the rayon pool below `--threads` (cheapest, confirmed supported)
A non-zero `ThreadPoolBuilder::num_threads(k)` *guarantees at most k threads*
([ThreadPoolBuilder docs]). The split is the experiment. The perf review's
starting point:
```rust
let rayon_threads = args.threads.map(|t| t.saturating_sub(t / 2).max(1));
configure_rayon_pool(rayon_threads)?;
```
This makes rayon + callers ≈ T. **Open question (measure it):** the compaction
`par_iter_mut` may be bursty enough that rayon can be capped to 1–2 threads with
the callers taking the rest — sweep T=2/4/8 and confirm the sweet spot is
`rayon + workers ≈ owned cores`. Risk: under-feeding decode if the rayon reserve
is too small.

### Option B — unify onto one pool (structural; removes the bridge block entirely)
Instead of a separate global rayon pool *plus* crossbeam threads, run the
`par_iter_mut` compaction **inside** a single capped `rayon::ThreadPool` via
`ThreadPool::install()` (install scopes `join`/`scope`/parallel-iterators onto
that pool — [ThreadPool docs]). Or replace the crossbeam fan-out with
`rayon::scope`/`join` so there is exactly one work-stealing pool sized to owned
cores. This collapses 2·T → T by construction and is the cleanest long-term shape
for a CPU-bound pipeline. Cost: this is the perf review's Variant 2 territory
(higher byte-identity risk because compaction moves relative to the callers) —
defer behind A unless A leaves residual contention.
[pkolaczk multiple-threadpools], [gendignoux rayon-optimized].

### Option C — serialize the compaction (Variant 3)
Drop `par_iter_mut` on `compact_samples`, run it serially on the producer thread,
let the callers carry all the parallelism. Worth testing because per-sample
compaction may not be worth parallelizing once the callers run the EM in
parallel — and it removes the producer's rayon pool entirely (one pool, T
callers). Small diff.

### CPU affinity / pinning — *not verified, treat as lower confidence*
The research did **not** surface authoritative evidence on `core_affinity` /
`hwlocality` tradeoffs for two CPU-bound pools sharing cores. General principle:
pinning helps cache locality for a *single* steady pool but can **hurt rayon's
work-stealing balance** (stealing assumes free migration). Do not pin as part of
H1; if tried later, A/B it against the context-switch counts. Listed as an open
question, not a recommendation.

### Nested parallelism pitfall (relevant, confirmed by our own code)
Running `rayon::par_iter` *inside* a crossbeam worker would multiply pools again.
Our merger already documents avoiding this
([per_group_merger.rs:27](../../../../src/var_calling/per_group_merger.rs#L27):
"running an inner `rayon::par_iter` here would only…"). Keep that invariant — the
EM workers must stay scalar; the parallelism lives at the chunk/caller level only.
[users.rust-lang nested-parallelism].

---

## 4. Best practices: rayon & crossbeam for CPU-bound work (confirmed)

- **Prefer a custom capped `ThreadPoolBuilder` + `install()` over the global pool**
  when another pool coexists. The global pool silently sizes to all logical CPUs.
  [ThreadPoolBuilder docs], [rayon FAQ].
- **`par_bridge` is a bad fit for an underfed serial iterator** — it pulls items
  one at a time under a lock and spins workers when the producer can't keep up
  (#730: 185 s with par_bridge vs 135 s without; #795: throughput *degrades* past
  8–16 threads on 64 cores, CPUs at 100 %). *We don't use it* (§0), but the rule
  stands: if a producer/consumer bridge is ever needed, batch the feed or use an
  explicit bounded channel, don't `par_bridge` a slow channel.
  [rayon #730], [rayon #795], [ParallelBridge docs].
- **crossbeam scoped threads + bounded channels** is the right shape for our
  pipeline (it bounds run-ahead → bounds RSS; perf review §7 already credits this).
  `std::thread::scope` (stable since 1.63) is now equivalent for the scoped-borrow
  case; crossbeam adds the bounded MPMC channel we rely on. Keep crossbeam for the
  channel.
- **Right thread count for CPU-bound work = owned cores, counted once.** Do not
  let two subsystems each independently claim "all cores."

---

## 5. The container / cgroup pitfall (confirmed — highest-risk default)

This is the no-`--threads` path and it is *more* exposed than the explicit sweeps.

- **`available_parallelism()` and the host core count overcount under a cgroup
  quota.** std docs: it "may overcount… when limited by a process-wide affinity
  mask or cgroup quotas." cgroup awareness landed only gradually: **v2 in Rust
  1.61, v1 in 1.64** (rayon PR #937; rust-lang PR #92697). Older toolchains, or
  reading raw core count, see the *machine*, not the *quota*.
  [std available_parallelism], [rust-lang #92697], [rayon #937].
- **Size from the quota:** effective CPUs = `quota / period`
  (cgroup v2 `cpu.max`; v1 `cpu.cfs_quota_us` / `cpu.cfs_period_us`). quota=period
  → 1 CPU; quota=2·period → 2 CPUs. `num_cpus` uses **ceiling** division (a
  1.5-CPU quota rounds up to 2 → itself a mild overcount). [kernel sched-bwc],
  [num_cpus linux.rs].
- **Why it compounds oversubscription:** CFS limits *amortized* usage. A highly
  parallel burst spends the whole quota in a fraction of the period, then **all
  threads are throttled (slept) until the next period** — "disastrous for tail
  latency" (Dan Luu). Layering 2·T threads on top of a quota the scheduler will
  throttle is the worst case. [kernel sched-bwc], [danluu cgroup-throttling].
- **Detect it:** `cpu.stat` → `nr_throttled` / `throttled_time`. Treat
  >~30 % throttling (perf review's own threshold) as "still oversubscribed."

**Action for our default path:** derive the pool size from the **owned cgroup
quota**, not `available_parallelism()`'s view of the limit, and split *that*
budget across the two pools. When validating on the Linux container, watch
`cpu.stat nr_throttled`/`throttled_time` alongside wall.

---

## 6. What this means for H1 (recommendation)

1. **The fix is pool-sizing/unification, confirmed.** Start with **Option A**
   (cap rayon, sweep the split T=2/4/8), gate on byte-identical calls
   (GT/GQ/AD/AF/AC/FILTER; QUAL exempt), and confirm via a fresh T=8 `sample`
   that `swtch_pri`/`__ulock_wait2`/`__psynch_mutexwait` fall toward T=1 levels —
   exactly the perf review's §3 item 3.
2. **Ignore the par_bridge remedies** the raw research surfaced — we use indexed
   `par_iter_mut`, not `par_bridge` (§0). The blocking is plain pool contention,
   not bridge starvation.
3. **Fix the container default separately:** size from the cgroup quota, not host
   cores, and verify with `cpu.stat`. This is the highest-risk default and the
   measurement plan should add a throttling check on Linux.
4. **L1 (producer serial floor) is downstream** — only re-measure it after A
   lands, because A may move the floor (perf review already gates L1 behind H1).
5. **Affinity/pinning and the crossbeam-vs-thread::scope-vs-rayon::scope
   head-to-head are unverified** — open questions, not recommendations.

---

## 7. Open questions (carried from the verification)

- Optimal asymmetric split (rayon vs caller threads) for *this* workload — is
  compaction bursty enough to cap rayon at 1–2?
- Does unifying onto one `install()`-scoped pool remove the bridge block entirely,
  and at what byte-identity/complexity cost (Variant 2)?
- Does disjoint pinning of the two pools help or hurt vs work-stealing balance?
- In the real container: what do `cpu.stat nr_throttled`/`throttled_time` and
  `available_parallelism()` report vs host cores — is CFS throttling compounding
  the oversubscription seen locally?

---

## Sources

Primary: [rayon FAQ](https://github.com/rayon-rs/rayon/blob/main/FAQ.md) ·
[ThreadPoolBuilder](https://docs.rs/rayon/latest/rayon/struct.ThreadPoolBuilder.html) ·
[ThreadPool](https://docs.rs/rayon/latest/rayon/struct.ThreadPool.html) ·
[ParallelBridge](https://docs.rs/rayon/latest/rayon/struct.ParallelBridge.html) ·
[rayon #730](https://github.com/rayon-rs/rayon/issues/730) ·
[rayon #795](https://github.com/rayon-rs/rayon/issues/795) ·
[rayon #937](https://github.com/rayon-rs/rayon/pull/937) ·
[std available_parallelism](https://doc.rust-lang.org/stable/std/thread/fn.available_parallelism.html) ·
[rust-lang #92697](https://github.com/rust-lang/rust/pull/92697) ·
[kernel sched-bwc](https://docs.kernel.org/scheduler/sched-bwc.html) ·
[num_cpus linux.rs](https://github.com/seanmonstar/num_cpus/blob/master/src/linux.rs)

Secondary: [PostHog — Untangling Tokio and Rayon](https://posthog.com/blog/untangling-rayon-and-tokio) ·
[Dan Luu — cgroup throttling](https://danluu.com/cgroup-throttling/) ·
[easyperf — false sharing/c2c](https://easyperf.net/blog/2019/12/17/Detecting-false-sharing-using-perf) ·
[pkolaczk — multiple thread pools](https://pkolaczk.github.io/multiple-threadpools-rust/) ·
[gendignoux — rayon optimized](https://gendignoux.com/blog/2024/11/18/rust-rayon-optimized.html) ·
[users.rust-lang — nested parallelism](https://users.rust-lang.org/t/nested-parallelism-and-rayon/84210)

**Refuted / killed (do not rely on):** par_bridge "batches reads under a mutex"
mechanism (0-3); manual-channel-beats-par_bridge anecdote (1-2); the Tokio
spawn_blocking three-pool mitigation (0-3, not applicable — no async here); the
30×-oversubscription arithmetic from the Medium post (0-3).
