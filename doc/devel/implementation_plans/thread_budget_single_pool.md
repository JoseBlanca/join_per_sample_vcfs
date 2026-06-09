# Thread budget: honour `--threads N` with a single work-stealing pool

**Branch:** `thread-budget-single-pool` (worktree off `main` @ `ca163df`).
**Status:** design proposal; Phase 0 (this audit) done. Phases 0.5–3 are
`[OPEN]`, each byte-identity- and thread-count-gated. The staged-pipeline
correctness bug it depended on — stale `IntervalEnd` not drained at interval
boundaries → dropped calls on multi-interval input — is **fixed and committed**
(`ab2431b`), so the worktree is unblocked. That bug is also a lesson for this
plan's test matrix (§5.1): it was masked by single-region tomato testing.

## Why we're here

`pop_var_caller` is distributed to the scientific community: it runs on
laptops, shared HPC nodes, and big servers alike. **We must not design for
our 32-thread box.** The governing contract is simple and portable:

> `--threads N` should use **about `N` parallel worker threads plus the main
> thread** — not a multiple of `N`. "main + `N` parallel" (i.e. `N + c` for a
> small constant `c`) is acceptable; `2N`–`3N` is not.

**Why this matters more than wall time.** The audience is largely HPC /
shared-cluster users. Under a batch scheduler (SLURM `--cpus-per-task`, PBS,
LSF) a job that requests `N` cores and spawns `~3N` threads is mis-behaving:
the scheduler pins the cgroup to `N` cores, so `3N` threads thrash on them, or
the job is flagged/killed for exceeding its allocation. Honouring the budget
is therefore a **correctness / good-citizen** requirement, not a performance
optimisation — and that means a **wall-neutral** outcome for this work is
fully acceptable. We are buying portability and well-behaved resource use, not
speed.

Today we violate this badly. Empirically (50 real tomato `.psp`, this
machine):

| `--threads` | threads actually spawned | contract target |
|---|---|---|
| 2 | **8** | ~3 |
| 5 | **19** | ~6 |

The count scales as **≈ 3·N + constant** — three *independent* `N`-sized
thread groups. This is the problem to fix.

This plan also records what the deep-research pass (2026-06-09) concluded,
because the fix must respect it: for CPU-bound work the idiomatic shape is
**one work-stealing pool sized to the core budget**, with bounded-queue
backpressure + work-stealing acting as the load balancer (not hand-assigned
per-stage thread counts); long-latency/**blocking** stages stay on dedicated
threads (rayon has no `ForkJoinPool`-style `ManagedBlocker`, so a pool worker
blocked on a channel starves the pool); and a **runtime auto-rebalancer is
overengineering** for a fixed-batch workload (the SEDA "adaptive beats static"
thesis did not survive verification). Sources: rayon FAQ / ThreadPool docs /
discussion #1193 / issue #319; Welsh's SEDA retrospective; the .NET CLR
ThreadPool controller case study; the TBB autotuning study; the JDK
`ForkJoinPool` javadoc.

## 1. The thread inventory today (the audit)

Every thread the `var-calling` path creates, with its count and lifetime:

| Source | count | lifetime / role |
|---|---|---|
| `configure_rayon_pool` global pool (`build_global`) | `N` | runs the parallel FASTA-MD5 verify at startup, then **parks idle for the whole run** |
| `producer_pool` (`pipeline.rs`) | `producer_threads` (= `N` by default) | read-decode + fold + compact `par_iter`s |
| caller threads (crossbeam scope) | `n_workers` (= `N` by default) | `call_chunk` (per-position merge + grouping + EM) |
| read-stage coordinator | 1 | staged path only (`producer_threads ≥ 4`) |
| compact-stage coordinator | 1 | staged path only |
| writer | 1 | ordered VCF sink, blocks on the called-chunk channel |
| main | 1 | the fold/plan driver |

Totals: **`3N + 4`** (staged, `N ≥ 4`) and **`3N + 2`** (inline, `N < 4`),
matching the measured 19 (`N=5`) and 8 (`N=2`).

**Diagnosis.** The offense is the **three `N`-sized groups**:

1. the global verify pool — does real work *once*, then `N` parked threads;
2. the producer rayon pool — `N`;
3. the caller threads — `N`.

Two of the three are rayon pools; the third is dedicated threads. The
coordinators (+2) and writer (+1) and main (+1) are a fixed constant, not the
problem. (`pileup` is a separate subcommand with its own `build_global` +
per-sample `par_iter` model; **out of scope here**, fix separately.)

## 2. The target contract

- **One** rayon work-stealing pool of `N` threads runs **all CPU work** of the
  pipeline: read/decode, fold, compact, *and* the EM caller. Work-stealing
  balances producer-side vs caller-side dynamically — which also deletes the
  hardware-dependent producer/caller split (the thing that made the optimum
  differ between a 6-core and a 32-thread box).
- Plus a **small fixed overhead**: the `main` thread (driver) and at most one
  writer (ordered I/O sink). Target total **≈ `N + 2`** (`c = 2`), independent
  of `N` and of the machine. (See §7 on whether the writer must fold into
  `main` for a strict `N + 1` — not a goal unless wanted.)
- `--threads` still defaults to `available_parallelism`, but the spawned
  count must *match the request* on every machine.
- **No** runtime controller, **no** topology/P-core tuning — rayon sizes to
  `N` and work-stealing absorbs heterogeneous cores; anything finer is left to
  the OS scheduler.

## 3. Design

### 3.1 One pool, passed by reference (not `build_global`)

Build a single explicit `rayon::ThreadPool` of `N` at startup and thread
`&pool` through the FASTA verify **and** the pipeline. We must *not* use
`build_global`: the criterion bench drives `pipeline::run_var_calling`
directly and sizes threads per iteration, which a process-global pool can't
do (this is exactly why `producer_pool` exists today). So:

- CLI `run_var_calling`: `let pool = ThreadPoolBuilder::new().num_threads(N).build()?;`
  run the verify in `pool.install(...)`; pass `&pool` to the pipeline; drop
  `configure_rayon_pool` from the var-calling path (keep it for `pileup`).
- The pipeline takes `pool: &ThreadPool` instead of building its own.

This alone collapses groups (1)+(2): the *same* `N` threads do the verify
then the pipeline — no second idle pool. (`T=5`: 19 → 14.)

### 3.2 The EM caller becomes pool tasks, not `N` dedicated threads

This is the group-(3) fix and the heart of the plan. Today `N` long-lived
caller threads pull `RawCohortChunk`s off a bounded channel and run EM. They
cannot simply move "into" the pool as long-lived tasks — they block on the
channel, and a blocked pool worker starves the pool. Instead, **dispatch each
chunk's call as a transient task on the one pool** as it is produced, so EM
work-steals alongside read/fold/compact. The writer still reorders by
`chunk_order`, so output stays in genomic order regardless of completion
order (unchanged).

**The main risk — producer starvation (design against it, don't discover
it).** Today the producer pool *always* has threads, so it keeps emitting
chunks no matter how busy the callers are. Collapsing producer work and EM
onto one shared pool removes that guarantee: rayon does **not** preempt, so a
burst of **long** EM tasks (EM can run hundreds of ms) can occupy all `N`
workers, stall the fold/produce, and then starve the callers later — i.e. one
pool can *serialize* a pipeline that is currently overlapped.

**Liveness guard (bake in from the start):** cap in-flight EM tasks at
**`< N`** so at least one worker can always advance the producer. The
bounded-in-flight count we need for memory (peak ≈ in-flight chunks) *is* the
liveness guard — one mechanism, two jobs. This is the hypothesis to validate.

**Plan B if the shared pool proves fragile:** a **budget partition** — `P`
producer threads + `N − P` caller threads, summing to `N` (still honours
`N + c`, keeps the always-a-producer-thread guarantee, just less elegant and
reintroduces a small split). Choose shared-pool-with-cap vs partition **from
the measurement (§0.5 / §5.1), not up front.**

Other care points:
- preserve the per-chunk batching / scratch-reuse the caller threads do today
  without relying on per-thread caller state;
- the EM (`call_chunk`) must not block (it doesn't use rayon internally and
  doesn't do I/O — verified), so it is safe as a pool task.

### 3.3 Trim the staged coordinators

With one pool, re-evaluate the read-stage / compact-stage coordinator threads
(queues 1+2). Keep a bounded queue **only** where a genuine
producer/consumer decoupling earns its keep (Welsh: isolate only
long-latency I/O — i.e. the cold `.psp` read), and express the rest as direct
calls / pool tasks. The inline-vs-staged threshold (`STAGED_MIN_THREADS`) may
become moot once there is a single model that scales from `N=1` upward.

### 3.4 Writer + main

The writer is an ordered, blocking I/O sink — keep it as **one** dedicated
thread, or fold it into `main` (main drains the called-chunk channel and
writes while the pool produces) for a strict `N + 1`. Decide in §7.

## 4. Relationship to the just-merged staged pipeline

Be honest: this **reshapes** the queues-1+2 topology we just landed. Those
optimised *throughput on a machine with spare cores*; the new, higher
priority is the **portable thread-budget contract**. What survives unchanged:

- the **straddler decode cache** (orthogonal CPU win, `--low-memory`-gated);
- **bounded-queue backpressure** as the balancer (kept where justified);
- **column-selective (two-phase) decode** and the fold/compact math.

What changes: the producer/caller **pool split** and the **per-stage
dedicated threads** give way to one pool + tasks. Net effect on wall time is
expected to be neutral-to-positive (our H1 result already showed the
multi-pool oversubscription is benign — the pools *park*, they don't
contend), but the win here is **correctness of the thread budget and
portability**, not speed.

## 5. Phasing (each phase: byte-identical + thread-count gated)

- **Phase 0 — audit.** This document + the inventory table.
- **Phase 0.5 — split-sweep measurement (no code change; run on current
  `main`).** Use the existing `PVC_PRODUCER_THREADS` / `PVC_CALLER_THREADS`
  knobs to sweep the *fixed* producer/caller split across the §5.1 matrix and
  find **where starvation happens** (which side is the wall, at which thread
  counts, for which inputs). **Purpose: decide the architecture, not to fit a
  table** — see §5.2. Primary question it must answer: *does one pool + an
  in-flight cap keep the producer fed across the matrix, or is a partition
  needed?* If one pool suffices, there is no table to write (best outcome).
- **Phase 1 — one rayon pool. [DONE — `b8d1fcc`]** §3.1. CLI builds one
  `ThreadPool` of `N`, runs verify via `pool.install`, passes `&pool` to the
  pipeline (which no longer builds its own producer pool); `build_global`
  dropped from var-calling (kept for pileup). `resolve_thread_split` →
  `resolve_thread_budget` + `resolve_caller_threads`; `PVC_PRODUCER_THREADS`
  retired (pool == budget). Removes the idle verify pool (−`N`): measured
  **3N+c → 2N+c** (peak threads at N=1/2/4/8 = 4/6/12/20; T=6: 22 → 16).
  Byte-identical (md5 `5164f285…`) at every N, both `--low-memory` modes, on the
  multi-interval tomato input; wall-neutral (N=8 5.7s). The remaining `N` is the
  caller threads → Phase 2.
- **Phase 2 — EM as pool tasks → fell back to Plan B partition. [DONE —
  `5f6c1b1`]** §3.2. The one-pool prototype (EM dispatched as transient
  `pool.scope` tasks behind an in-flight cap) was built and measured: it is
  **byte-identical** but **~2× slower than a same-budget partition at every
  thread count** (N=6/12/18: 19.4/13.5/11.1s vs Phase-1 2N's 6.9/5.3/5.4s; the
  same-budget partition was ~8.8s at N=8). Root cause exactly as §3.2 warned —
  heavy *non-preemptible* EM tasks steal the decode-bound producer's threads, so
  the producer (the wall limiter) runs *below* its floor; the in-flight cap does
  not fix it (the producer wants bursty access to *all* threads, which a fixed
  cap can't grant). So per the pre-specified fallback we took the **budget
  partition**: a `P`-thread producer pool + `C` dedicated caller threads,
  `P+C=N` (`resolve_split`), default balanced with a small-cohort nudge (§5.2).
  Result: **3N+c → N+c** (c=2 inline / c=4 staged; N=8: 28 → 12 threads),
  byte-identical at N∈{1,2,4,6,8} on multi-interval input, wall at the
  honest-budget partition optimum. The one-pool code was reverted (kept only in
  history); [[project_thread_budget_split_sweep]] records the measurements.
- **Phase 3 — trim coordinators / revisit staging.** §3.3. Collapse to the
  minimum justified threads; re-evaluate `STAGED_MIN_THREADS`. Target
  `≈ N + 2`. **Bonus:** Phase 2 may let the inline and staged paths converge
  into one model (deleting `STAGED_MIN_THREADS`) — adopt that *only* if
  one-pool matches the staged overlap at higher thread counts (verify at
  `T=16` as a rough big-machine proxy; note this box's E-cores).

### 5.1 Measurement matrix (our current benchmarks don't represent the audience)

We have been measuring warm-cache, `N=50` samples, single-region tomato, on a
6-fast-core box — none of which is typical of the community's machines or
inputs. The multi-interval bug (Status) was *masked* by exactly this gap.
Every phase is judged across:

- **threads: 1, 2, 4, 8** — small counts are the common case (and where the
  staged path already regressed once); do **not** optimise for high `T`.
- **inputs: include a multi-interval / multi-chromosome cohort** (e.g. the
  human genome bottle), not just single-region tomato — the interval-boundary
  sync is a known fragility (Status).
- **cache: warm + at least one cold run** — warm runs hide blocking-read
  starvation, which is exactly where the one-pool / I-O-isolation question
  lives.
- **a couple of cohort sizes** — the producer-vs-EM balance shifts with sample
  count (we saw stale-vs-fresh psp move caller starvation 41%→11%).

Gate (all must hold, every phase): **byte-identical VCF** (md5 of non-`##`
lines) at every `N`, **on multi-interval input**, and both `--low-memory`
modes; **thread count ≤ `N + c`** asserted by an integration test and logged
at startup (this is what stops the budget silently drifting back to `3N`);
**no small-machine wall regression**.

### 5.2 On a static thread-allocation table (decision criteria)

The temptation is to bake a per-`N` lookup table mapping thread count → split.
Resist a *machine-fitted* table, for three reasons:

1. **Portability** — a table measured here encodes *our* hardware's bottleneck
   profile; on a different laptop / disk / cache / core-mix the optimum moves.
   Hard-coding it re-introduces the very problem we are solving.
2. **It re-commits to the split we are dissolving** — a per-`N` split only has
   meaning under the partition (Plan B). The shared-pool target makes the
   split *disappear*; if it works, there is no table.
3. **Under-parameterized** — the balance also depends on cohort size, cache
   warmth, and `--low-memory`; a 1-D table keyed on `--threads` is wrong
   off-axis.

What is acceptable, if a split proves unavoidable: the **simplest robust rule
that generalizes** (a *regime threshold* like the existing
`STAGED_MIN_THREADS`, which encodes a qualitative transition rather than
fitted per-`N` optima), env-/flag-overridable for power users. Order:
**measure → decide one-pool-vs-partition → if partition, derive the simplest
robust rule, not a fitted per-`N` table.** The ideal result is *no table*.

**Chosen rule (`resolve_split`, `5f6c1b1`).** The threads × cohort-size sweep
(2026-06-09) showed the optimal *ratio* is ~`0.5` and **stable across `N`** — so
no per-`N` table — but it **shifts with cohort size**: at `N=8` the best `P` is
4 / 5 / 6 for 50 / 8 / 1 samples (heavy EM at many samples pulls toward
balanced; trivial EM at few samples pulls toward producer). The rule is
therefore a *ratio*, not counts: `P = ⌊N/2⌋` (balanced) for cohorts above
[`SMALL_COHORT_SAMPLES`]` = 8`, nudged to `P = ⌈2N/3⌉` for cohorts at/below it.
This keys on `n_samples` (an **input** property, not hardware), so it stays
portable; the threshold sits at the largest cohort still measured to prefer the
nudge, because over-nudging a large cohort is the steep-penalty direction
(N=8/50-sample P=6 → +67 % vs balanced) while under-nudging a small one costs
only tens of ms on sub-2 s runs. Env knobs override each side.

## 6. Hard constraints

- **Byte-identical VCF** (md5 of non-`##` lines) at every phase, every `N`
  (1, 2, 4, …), **on multi-interval input**, both `--low-memory` modes.
- **Thread count ≤ `N + small_constant`**, asserted and logged at startup
  (an integration test spawns a run and counts threads).
- **No deadlock**: pool tasks must never block on the pool; blocking sinks
  (writer) stay dedicated threads.
- **Scales down**: correct and reasonable at `N = 1` and `N = 2`; no design
  assumption of "spare cores." Do not regress small-machine wall.
- One tunable knob (`--threads`) + a sensible default; nothing auto-tuning.

## 7. Open questions

- **Writer placement:** acceptable as a dedicated `+1` (near-idle I/O sink,
  not a CPU competitor), or must it fold into `main` for a strict `N + 1`?
- **EM batching:** can per-chunk EM run as transient pool tasks while keeping
  the current scratch-reuse / batch behaviour without per-thread state?
- **I/O isolation:** keep one dedicated reader thread for cold `.psp` reads
  (Welsh) or fold reads into pool tasks? Needs a **cold-cache** benchmark.
- **Staging after Phase 2:** is any read/compact stage still worth a bounded
  queue, or does one pool + tasks subsume it?
- **`pileup`:** same `3×N`-ish issue? In a follow-up, document and align it.

## 8. Non-goals

- A runtime auto-rebalancing controller (research-refuted for this workload).
- P-core/E-core/NUMA-aware sizing or pinning (leave to the OS scheduler).
- Maximising 32-thread throughput at the expense of the thread-budget
  contract or small-machine behaviour.
