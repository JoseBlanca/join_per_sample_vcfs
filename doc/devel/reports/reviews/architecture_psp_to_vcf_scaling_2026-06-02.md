# Architecture review — the `.psp` → cohort-VCF path and why it stops scaling with threads

**Date:** 2026-06-02
**Scope:** the cohort `var-calling` path only — from the per-sample `.psp`
artefacts to the final cohort VCF. PSP *creation* (`pileup`, one process per
sample) is explicitly out of scope: it parallelises trivially across samples
and is not a scaling concern here.
**Purpose:** reconcile the remembered high-level architecture against the code
as it actually stands, fill in the details, and organise what we know about the
thread-scaling plateau so the next optimisation pass has a shared map.

All file/line references are against the working tree at the time of writing.

---

## 1. TL;DR

The remembered architecture is essentially correct. The path is a **3-stage
pipeline** wired with hand-rolled queues:

```
                 ┌─────────────────────────────────────────────────┐
   N .psp files  │  PRODUCER (1 std::thread)                        │
   ───────────►  │   BlockIterator::next_block                      │
                 │    • decode N samples' columns  (rayon par_iter)  │
                 │    • fold → cohort variable sites (serial)        │
                 │    • cut at a clean group boundary (serial)       │
                 │    • slice DUST mask / partition / prefetch REF   │
                 └───────────────────┬─────────────────────────────┘
                                     │  BlockQueue (bounded MPMC, cap = 2·W)
                                     ▼
                 ┌─────────────────────────────────────────────────┐
   DUST-ahead    │  WORKERS (W = --threads dedicated std::threads)  │
   pool ─────────►   process_block → run_window                     │
   (W threads,   │    • per-group merger + EM posterior (serial     │
    masks via    │      over groups within the block)               │
    ordered Q)   └───────────────────┬─────────────────────────────┘
                                     │  mpsc::channel<BlockResult>
                                     ▼
                 ┌─────────────────────────────────────────────────┐
   final VCF ◄───│  COLLECTOR (the calling thread)                  │
                 │   reorder by seq_idx → emit_or_drop → VCF writer │
                 └─────────────────────────────────────────────────┘
```

The headline finding: **only the middle (consume) stage is replicated across
`--threads`.** Production is a single thread, and emit is a single thread. So
once you have ~2–3 workers the consume stage stops being the bottleneck and the
wall clock pins to the **single producer thread**. That is exactly what the
thread-sweep benchmark shows (§5): ~2× from 1→2 threads, then flat.

---

## 2. Your mental model vs. the code

| # | What you remembered | Verdict | Where it lives |
|---|---|---|---|
| 1 | Read `.psp`, one block at a time; process and yield processed blocks via an iterator | **Correct**, with one refinement: a "block" is *not* one on-disk PSP block. It is a unit of work cut at a variant-count target on a clean group boundary; producing one decodes PSP chunks from **all N samples in parallel** and folds them to the cohort's variable sites. | `BlockIterator` / `next_block` — [driver.rs:1419-1772](../../../src/var_calling/driver.rs#L1419-L1772) |
| 1b | DUST is done by a pool of workers that send results via a queue consumed by the block-processing code | **Correct, and precise.** `DustAheadPool` runs `sdust` over the (independent) covered intervals in parallel and delivers masks in genomic order through an ordered queue; the producer just *slices* the precomputed mask per block. | `DustAheadPool` — [driver.rs:1119-1303](../../../src/var_calling/driver.rs#L1119-L1303); consumed at [driver.rs:1597-1608](../../../src/var_calling/driver.rs#L1597-L1608) |
| 2 | Blocks go through the EM/posterior machinery; "can we run several in parallel?" | **Yes — across blocks.** `W = rayon::current_num_threads()` dedicated worker threads each run `process_block`/`run_window` on a different block. Within a single block, groups are processed **serially** (the old per-group `par_iter` was deliberately removed — "L1"). | workers — [driver.rs:534-568](../../../src/var_calling/driver.rs#L534-L568); `run_window` group loop is serial — [worker.rs:248-269](../../../src/var_calling/worker.rs#L248-L269) |
| 3 | Resulting variants are written to a VCF | **Correct**, via a single **collector** thread that re-imposes genomic order (`seq_idx`) before emitting, so the VCF is byte-deterministic regardless of `W`. | collector — [driver.rs:574-607](../../../src/var_calling/driver.rs#L574-L607); `emit_or_drop` — [driver.rs:702-743](../../../src/var_calling/driver.rs#L702-L743) |

So: nothing in your recollection is wrong. The detail worth internalising is the
**threading topology** — which of these stages actually multiplies with
`--threads`, and which do not.

---

## 3. The threading topology in detail

The whole pipeline is set up in `drive_blocks_parallel`
([driver.rs:474-625](../../../src/var_calling/driver.rs#L474-L625)) inside a
`std::thread::scope`. Three distinct concurrency mechanisms coexist:

1. **One producer thread.** Spawned once; loops `producer.next_block()` and
   pushes into the `BlockQueue`. This is the *only* thread that touches the
   rayon pool — and it does so *inside* a block (see below). There is exactly
   **one** producer regardless of `--threads`.

2. **`W` worker threads** (`W = rayon::current_num_threads().max(1)`,
   [driver.rs:399](../../../src/var_calling/driver.rs#L399)). These are
   **dedicated `std::thread`s, not rayon tasks** — and the doc comment at
   [driver.rs:462-466](../../../src/var_calling/driver.rs#L462-L466) spells out
   why: a worker blocked on `queue.pop()` must not occupy a rayon pool thread,
   because the producer's per-sample decode *uses* that pool; parking workers on
   it could starve or deadlock the decode. (This is the rayon-vs-crossbeam point
   from our earlier conversation, made concrete.)

3. **The DUST-ahead pool** (`DustAheadPool`,
   [driver.rs:1267-1303](../../../src/var_calling/driver.rs#L1267-L1303)): its
   own set of `std::thread`s (also sized from `rayon::current_num_threads()`),
   running entirely off the critical path, feeding masks through an ordered
   bounded queue.

4. **The collector** is the calling thread itself
   ([driver.rs:574](../../../src/var_calling/driver.rs#L574)).

### Where rayon is actually used

Rayon (`par_iter`) appears **only inside the producer**, within `fill_block`:

- per-sample decode `read_span` across all N samples —
  [loader.rs:638-642](../../../src/var_calling/loader.rs#L638-L642)
- the compaction split across all N samples —
  [loader.rs:670-694](../../../src/var_calling/loader.rs#L670-L694)

The serial spine of the producer (the watermark min-scan
[loader.rs:619-622](../../../src/var_calling/loader.rs#L619-L622),
`rebuild_position_union_and_is_kept`, and `find_block_cut`
[loader.rs:644-649](../../../src/var_calling/loader.rs#L644-L649)) is **not**
parallel — it is the cohort fold, and it runs on the single producer thread.

### Consequence: cores are oversubscribed

At `--threads N` you have, concurrently: 1 producer + N workers + (up to) N
DUST threads + an N-thread rayon pool that the producer drives. The design bets
that most of these are blocked at any instant (workers waiting on the queue,
DUST mostly done early, rayon idle except during a decode burst). That bet holds
only while one stage clearly dominates; near the plateau it means the producer's
decode bursts and the workers' EM are fighting for the same physical cores.

---

## 4. Memory model (for context — it's the project's headline)

The pipeline is bounded, not buffered. The `BlockQueue` is bounded to
`cap = 2·W` blocks ([driver.rs:487](../../../src/var_calling/driver.rs#L487));
the producer blocks on a full queue (back-pressure). Spent block buffers are
**recycled** to the producer's free-list rather than freed
([driver.rs:1535-1537](../../../src/var_calling/driver.rs#L1535-L1537),
`ReadyBlock`). So resident memory is ~`(2·W + W)` blocks × N samples plus the
DUST look-ahead — never the whole span × N. This is the RAM-for-scaling trade
the project is built around, and the thread sweep confirms RSS stays modest
(§5). It is **not** the scaling problem; it is the thing the scaling problem
must not break.

---

## 5. The scaling evidence

`benchmarks/tomato1/results/perf/ours_joint_threads.tsv` — one `var-calling`
process over the **full 50-sample** cohort PSPs, sweeping `--threads`:

| threads | wall (s) | speedup vs T=1 | peak RSS (MB) |
|--------:|---------:|---------------:|--------------:|
| 1 | 60.05 | 1.00× | 498 |
| 2 | 32.99 | 1.82× | 617 |
| 3 | 30.31 | 1.98× | 761 |
| 4 | 29.71 | 2.02× | 862 |
| 6 | 29.61 | 2.03× | 916 |
| 8 | 30.56 | 1.97× | 834 |

The curve flattens hard after T=2. Effective ceiling ≈ **2×**, reached by T=3;
T=8 is fractionally *slower* than T=4 (oversubscription cost). By Amdahl this
implies a serial fraction near **0.5** — about half the work is on a
non-replicated stage.

For scale, the cohort-size sweep at the operating thread count
(`ours_joint.tsv`) is healthy and well ahead of the competition (e.g. N=50:
ours 31.7 s / 831 MB vs `gatk_joint` 447 s / 765 MB vs `freebayes` 1424 s /
9617 MB). **The problem is not absolute speed — it's that we cannot buy more
speed with more cores.**

---

## 6. Diagnosis — what the serial half is

The pipeline wall clock is bounded below by the slowest *non-replicated* stage.
Two stages do not replicate with `--threads`:

- **The single producer.** Even with parallel per-sample decode, blocks are
  produced strictly one at a time on one thread, and the cohort fold
  (watermark scan + position-union rebuild + cut finding) is serial on that
  thread. As `W` grows, the consume stage's throughput (`W / C`) rises until it
  exceeds the producer's `1 / P`; past that point the producer is the wall.
  The plateau at ~30 s is, on this reading, the **producer's standalone
  throughput** at N=50.

- **The single collector + VCF writer.** Reorder-and-emit is one thread. Likely
  cheap today (it's bookkeeping + write), but it is a hard serial tail and
  worth confirming isn't a secondary ceiling.

A profile is needed to *attribute* the ~30 s producer floor between (a) the
serial cohort-fold spine, (b) the parallel decode that still has to happen for
all 50 samples on a shared pool, and (c) partition + REF prefetch. The
2026-05-27 scaling measurement already hinted at the shape: at large N the
consume side grows ~N^1.8 and the produce side (`psp_reader` + fold +
`dust_filter`) is a large N-dependent serial component.

---

## 7. Candidate levers (research backlog, not decisions)

Ordered by how directly they attack the §6 diagnosis. None is committed; each
needs a profile to size and must preserve the byte-identity contract.

1. **Parallelise production across covered intervals.** The covered intervals
   are already provably independent — that independence is the whole
   byte-identity argument, and the DUST pool *already* exploits it. A pool of
   producer threads, each owning a disjoint set of intervals and emitting into
   the ordered collector by a global `seq_idx`, would turn the single producer
   into a replicated stage. This is the most direct attack on the ceiling.
   *Risk:* the rayon-based per-sample decode would then need rethinking (N
   producers each driving an N-wide `par_iter` is pure oversubscription) — decode
   probably becomes serial-per-sample within each producer, with parallelism
   coming from running many producers instead.

2. **Shrink the serial cohort-fold spine.** Profile `fill_block`'s serial
   parts (watermark min-scan, `rebuild_position_union_and_is_kept`,
   `find_block_cut`) and cut their cost, so a single producer simply goes
   faster. Cheaper to attempt than #1, smaller ceiling.

3. **Right-size the thread budgets independently.** Today producer-rayon,
   workers, and DUST all derive from one `--threads` knob and oversubscribe.
   Splitting the budget (e.g. a small decode pool + more consume workers, or
   fewer DUST workers once the masks are computed) may reclaim the T=8
   regression without touching the architecture.

4. **Confirm the collector isn't a secondary tail.** Quick to measure;
   determines whether emit needs its own attention after the producer is fixed.

5. **Revisit per-block granularity.** `target_variants_per_chunk` sets the
   block size. Smaller blocks = finer-grained pipelining (workers start sooner,
   producer hands off sooner) at the cost of more per-block overhead. A sweep
   could shift where the producer/consumer crossover lands.

---

## 8. Open questions for the profiling pass

- At N=50, how does the ~30 s producer floor split between parallel decode vs.
  serial fold vs. partition/prefetch? (Determines whether lever #1 or #2 pays.)
- Is the DUST pool fully done before the producer needs masks, or is
  `recv_next` ever a stall point at low thread counts?
- What is the collector's self-time? Is emit ever back-pressuring the workers
  (mpsc is unbounded, but the queue cap transitively bounds it —
  [driver.rs:488-498](../../../src/var_calling/driver.rs#L488-L498))?
- Does the T=8 regression vanish if the three thread budgets are decoupled
  (lever #3), confirming oversubscription as its cause?

---

## 8b. Measured chokepoints — profiling pass (2026-06-02)

Profiled the **real** `var-calling` path (not the synthetic example) on the
50-sample tomato PSP cohort at `--threads 6`, the exact regime of the §5
plateau. Two lenses: lightweight queue-stall counters added to `BlockQueue`
(`POP_VAR_CALLER_QUEUE_STATS=1`; touched only on the blocking path, zero
hot-path cost) for *which stage is the wall*, and macOS `sample` (1 ms, 30 s)
for *within-stage self-time*. Raw artefacts under `tmp/scaling_profile/`.

Run: wall **31.75 s**, 141 s user + 22 s sys (≈5.1 cores busy avg), peak RSS
920 MB, 5501 blocks produced (avg 1153 variants/block), 185 k records emitted.

### Finding 1 — definitively producer-bound (queue stalls)

| counter | value | reading |
|---|---|---|
| `producer_push_waits` | **3** (0.02 s total) | producer **never** blocks on a full queue |
| `worker_pop_waits` | 5345 (**128.4 s** total) | workers block on an *empty* queue constantly |
| per-worker idle | 21.4 s / 31.75 s = **67 %** | each of 6 workers starves two-thirds of the run |

The workers wait on the producer; the producer never waits on the workers. More
workers cannot help. **Confirms the §1/§6 hypothesis with hard numbers.**

### Finding 2 — the serial fold is NOT the problem (lever #2 is dead)

`rebuild_position_union_and_is_kept` self-time: **115 samples** — negligible.
The producer's serial cohort-fold spine is already cheap. Lever #2 ("shrink the
serial fold") would optimise nearly-dead code; **drop it.**

### Finding 3 — the producer's real cost is decode + rayon dispatch overhead

Self-time, summed across all 20 threads (1 ms samples), bucketed:

| bucket | ~samples | % on-CPU | notes |
|---|--:|--:|---|
| Consume / EM (workers) | 35.0 k | ~43 % | `run_window` 21.3 k, `unify_alleles_columnar` 8.7 k (still hot post-H2/H3), `e_step_simd` 3.3 k — *parallel across 6, and idle 67 %* |
| **rayon plumbing / idle** | 13.5 k | ~16 % | `bridge_producer_consumer::helper` 8.6 k, `rayon::slice::sort` 2.4 k, steal/wait/join 2.5 k |
| **DUST (`sdust`)** | 17.7 k | ~22 % | `spawn_dust_pool` closure — the **2nd-hottest single symbol** |
| **PSP decode (zstd)** | 11.7 k | ~14 % | `ZSTD_decompressSequences` 5.7 k, `decode_list_column_csr` 2.1 k, XXH64 1.6 k, … |
| producer fold / compact | 2.8 k | ~3 % | `push_row_from` 2.6 k, `rebuild_union` 0.1 k |
| collector / VCF write | 0.4 k | ~0.5 % | not a tail — confirms the open question |

The producer drives a **per-block, 50-sample-wide `par_iter`** twice per block
(decode + compact) × 5501 blocks ≈ **275 k tiny rayon tasks**. The rayon
dispatch machinery (13.5 k) burns roughly as much CPU as the decode it
dispatches (11.7 k), and the 6 rayon-pool threads sit in `wait_until_cold`
(idle) most of the run because the producer serialises *between* dispatches.
**The fine-grained par_iter is overhead-bound.**

### Finding 4 — the DUST pool is a hidden CPU hog (thread map)

Thread classification from the call graph: 1 collector (main), **6 rayon-pool
threads** (mostly idle / `wait_until_cold`), **6 worker threads** (`run_window`,
idle 67 %), and **6 DUST threads** that burn ~17.7 core-seconds of `sdust` and
**die ~4 s in**. So during the first ~4 s the machine runs ~20 threads
(6 DUST + 6 rayon + 6 workers + producer + collector) on ~8 cores — heavy
oversubscription at startup. "Off the producer's critical path" logically, but
very much *on* the CPU-contention path early. Modest wall cost on a 31.7 s run,
but a real inefficiency and a complexity burden.

### Finding 5 — within production, READING dominates GROUPING ~2:1

Deterministic timers added to the two logical sub-steps of production
(loader `read_nanos`/`fold_nanos`, driver `partition`/`prefetch`/`mask_slice`),
N=50 T=6, three runs, producer-thread wall summed across all 5501 blocks:

| sub-step | time | share of producer busy |
|---|--:|--:|
| **Step 1 — reading pileup (`read_span`, PSP/zstd decode)** | **~17.4 s** | **~68 %** |
| **Step 2 — creating grouped variants** | **~8.3 s** | **~32 %** |
| &nbsp;&nbsp;· cohort variant-position fold + compaction | 7.0 s | |
| &nbsp;&nbsp;· `partition_window` (variant grouping) | 1.18 s | |
| &nbsp;&nbsp;· REF prefetch | 0.12 s | |

Stable to ±1 % across runs. **Reading is the dominant cost inside block
production, ~2× the grouping.** This is the *easily-parallelisable* case:
decode of one PSP region is independent of the next. The grouping step is real
but secondary (and its bulk is the compaction `push_row_from`, not the
variant-position decision, which is near-free).

**Memory caveat for a read/group split:** the variant-vs-non-variant decision
is *cohort-wide* (a position is kept if **any** sample carries a non-ref allele,
plus MNP reach), so a per-sample reader **cannot** pre-drop non-variant records
on its own — it lacks the other samples' view. The drop must happen right after
a span is read across *all* samples (the cohort fold). So the raw,
not-yet-dropped pileup exists only briefly between "read" and "fold/drop"; its
footprint is bounded by the queue's back-pressure depth, not unbounded. Two
clean shapes both respect the memory budget: (a) parallel readers → a small
fold/drop → grouping; or (b) parallelise read+fold+group *per covered interval*
(the producer-pool design), which keeps today's exact memory model.

### Finding 6 — the read bottleneck is PSP block MISALIGNMENT across samples (root cause)

Follow-up after Finding 5 showed reading dominates but does **not** parallelise
(1→6 threads only 1.2×; flat 2→12; unchanged by output block size, by removing
DUST, and not disk-bound since files are cached).

The read path *does* follow the PSP chunk structure correctly
([column_span_reader.rs](../../../src/var_calling/column_span_reader.rs)): the
producer peeks every sample's next on-disk block boundary, advances to the
**minimum** (`W = min peek`), and each sample decodes at most one new block per
step — every block decompressed exactly once. Each step is one `par_iter` over
all samples **followed by a fold barrier**.

Instrumented step count (`read_calls`), N=50 T=6:

| cohort | `read_calls` (sync steps) | STEP1 read | per-step |
|---|--:|--:|--:|
| **real 50 (independently-written files)** | **49,017** | **17.5 s** | 358 µs |
| **50 identical replicas (aligned blocks)** | **920** | **3.7 s** | — |

The 53× step inflation ≈ the sample count. Because each sample's `.psp` is
written independently with content-driven block cuts, the boundaries **don't
line up across samples**; `W = min peek` advances to the *union* of all
samples' boundaries, so at nearly every step only one sample crosses a real
boundary and the other 49 read a sliver and wait at the barrier. The decode of
one small block is a few µs; the per-step `par_iter` dispatch + barrier is
~hundreds of µs — so with 49 k steps the **barriers dominate and the decode
parallelism is drowned**. Output block size can't fix it (steps are set by PSP
boundaries, not output blocks); thread count can't fix it (the wall is barriers,
not decode throughput).

**Three candidate fixes, increasing invasiveness:**

1. **Coarser consensus span (cheapest).** Drive `read_span` with a larger `end`
   than `min peek` — decode *many* blocks per sample per step, so the barrier is
   paid far less often. `read_span` already accepts any `end`; only the
   producer's watermark choice changes. Trades a little more in-flight memory
   per step for ~50× fewer barriers. **Test this first.**
2. **Align PSP block boundaries at write time.** If `pileup` cut blocks on a
   shared genomic grid, the cohort read would sync once per shared block
   (~920 not ~49 k). Requires a PSP-writer change + regenerating `.psp` files.
3. **Decouple readers from the fold barrier.** Let each sample's reader run
   ahead independently (read many blocks without syncing), folding separately —
   the "parallel readers" shape. Biggest change; works on existing files.

NB this also reframes the producer-pool idea: a producer *pool over intervals*
helps the orchestration, but the **per-step barrier from misalignment is the
real read ceiling** — fix #1/#2/#3 is what unblocks reading, independent of how
many producers run.

### Finding 7 — fix #2 implemented and validated: genomic-window block cut (2.6× wall, byte-identical)

Implemented the aligned block cut: the PSP writer now cuts blocks on a fixed
genomic grid (`pos / block_window_bp`, default **20 kb**) instead of by
accumulated bytes; the byte target became a generous safety cap. New CLI knob
`pileup --block-window-bp`. (Writer: [writer.rs](../../../src/psp/writer.rs)
`new_with_block_layout` + `DEFAULT_BLOCK_WINDOW_BP`; CLI:
[cli.rs](../../../src/pop_var_caller/cli.rs).)

Validated by re-chunking the existing 50 real PSPs to a 20 kb grid (same
records, new boundaries — `examples/psp_rechunk.rs`) and re-running
`var-calling`, N=50:

| metric (T=6) | misaligned (byte-cut) | aligned (20 kb grid) | change |
|---|--:|--:|--:|
| sync steps (`read_calls`) | 49,017 | 3,380 | **14.5× fewer** |
| STEP1 read | 17.5 s | 3.8 s | **4.6× faster** |
| worker idle / worker | 19.1 s (67 %) | 0.23 s (~2 %) | starvation gone |
| **wall** | 29.7 s | **11.35 s** | **2.6× faster** |

- Per-sample blocks became uniform: ~480 (was 618–1854), ~16 k records each
  (the equal-memory goal — record count is now depth-independent and uniform).
- Cross-sample alignment ratio 6.0× → **1.9×** (residual is covered-interval
  *edges*, which differ per sample; interior blocks are perfectly grid-aligned).
- **Output byte-identical**: 185,279 records, md5 match vs the misaligned run.
  The hard correctness contract holds.

**The bottleneck has moved.** Reading now parallelises (1→6 threads: 17.6→3.8 s).
At T=6 the workers are ~98 % busy and the producer *occasionally* waits on
them — produce (~10 s single-thread: 3.8 read + 6.1 fold/compact) and the
6-way EM consume are now roughly balanced at ~10–11 s. The next levers are the
ones from the original redesign: **parallelise production** (a producer pool
over the independent covered intervals — now worth it, since reading scales),
and/or the EM consume side (`unify_alleles`, e-step).

Status: committed on branch `psp-aligned-block-window` — writer grid-cut +
`--block-window-bp` CLI knob + `block_window_bp` PSP provenance + a writer
grid-cut unit test; the temporary instrumentation was reverted. fmt / clippy
(`--all-targets -D warnings`) / 1065 lib + integration tests green. The
calibration tools `examples/psp_block_stats.rs` and `examples/psp_rechunk.rs`
ship alongside. Production `.psp` files must be regenerated via `pileup` to
inherit the 20 kb default (the validation used `psp_rechunk` to avoid re-running
pileup).

**Production benchmark (PSPs regenerated via `pileup` at the 20 kb default,
tomato N=50).** Thread sweep — the ~2× plateau is gone:

| threads | before (misaligned) | after (20 kb) |
|--:|--:|--:|
| 1 | 60.0 s | 56.2 s |
| 2 | 33.0 s (1.8×) | 29.0 s (1.9×) |
| 4 | 29.7 s (2.0×) | 16.8 s (3.4×) |
| 6 | 29.6 s (2.0×) | 12.5 s (4.5×) |
| 8 | 30.6 s (2.0×) | 10.8 s (**5.2×**) |

Cohort-size sweep (T=4): wall N=50 31.7 → 16.6 s; N=8 7.9 → 2.2 s. Peak RSS rose
~6× (~16 → ~100 MB/sample; N=50 831 MB → 5.0 GB) — memory stays **linear** in N
(the no-OOM-blowup thesis holds, just a steeper line) because the bigger blocks
now keep the bounded queue full. **PM decision: keep the 20 kb default — the
memory increase is acceptable for this workload.** TSVs at
`benchmarks/tomato1/results/perf/ours_joint{,_threads}.tsv`.

Next levers (now that reading scales): parallelise the single producer thread
(pool over the independent covered intervals) and/or the EM consume side
(`unify_alleles`, e-step).

### What this means for the redesign

- **Lever #1 (producer pool over intervals) is the right move** — and the
  decisive refinement is confirmed by Finding 3: **drop the per-block
  per-sample `par_iter`; decode serial-per-producer; get parallelism from
  running P producers over independent intervals.** That deletes the ~13.5 k
  samples of rayon dispatch overhead *and* fills the idle worker + idle
  rayon-thread capacity. The idle headroom (workers 67 % idle, rayon pool mostly
  idle) is exactly the capacity a producer pool would consume.
- **Fold DUST into the producers** (each DUSTs its own interval) — removes the
  separate 6-thread pool and its startup oversubscription (Finding 4). The
  `DustAheadPool` + ordered `recv_next` coordination disappears.
- **Lever #2 (serial fold) is dead** (Finding 2). **Collector is not a tail**
  (Finding 3) — leave it.
- The per-producer ceiling is then **PSP zstd decode throughput** (Finding 3,
  11.7 k irreducible-ish) — so P producers ≈ P× until cores saturate (~8 here).
- `unify_alleles_columnar` (8.7 k, consume) remains the hottest consume kernel
  after the producer is fixed — the next target once scaling is unblocked.

Instrumentation note: the `BlockQueue` stall counters were added to
[driver.rs](../../../src/var_calling/driver.rs) for this pass and left in place
(env-gated, no hot-path cost). Not yet run through container fmt/clippy/test.

---

## 9. Key code map (for the next reader)

| Concern | File / symbol |
|---|---|
| Top-level entry, opens readers + writer, sizes `W` | `drive_cohort_chunked` — [driver.rs:308-431](../../../src/var_calling/driver.rs#L308-L431) |
| Pipeline wiring (producer/workers/collector, queues) | `drive_blocks_parallel` — [driver.rs:474-625](../../../src/var_calling/driver.rs#L474-L625) |
| Bounded MPMC block queue (back-pressure) | `BlockQueue` — [driver.rs:633-697](../../../src/var_calling/driver.rs#L633-L697) |
| Block production (decode → fold → cut → DUST/partition/prefetch) | `BlockIterator` — [driver.rs:1419-1772](../../../src/var_calling/driver.rs#L1419-L1772) |
| The cohort fold itself (rayon decode + serial union/cut) | `fill_block` — [loader.rs:600-707](../../../src/var_calling/loader.rs#L600-L707) |
| Parallel DUST, ordered delivery | `DustAheadPool` — [driver.rs:1119-1303](../../../src/var_calling/driver.rs#L1119-L1303) |
| Per-block math (merger + EM), serial over groups | `run_window` — [worker.rs:248](../../../src/var_calling/worker.rs#L248) |
| Ordered emit + downstream filters | `emit_or_drop` — [driver.rs:702-743](../../../src/var_calling/driver.rs#L702-L743) |
| The scaling evidence | `benchmarks/tomato1/results/perf/ours_joint_threads.tsv` |
